## main module for parsing, processing and interacting with data/files of 2p-stim protocols (i.e. NAPARM outputs)
from __future__ import absolute_import
from dataclasses import dataclass
from typing import TypedDict, Optional, MutableMapping
import re
import tifffile as tf
import os
import numpy as np
import xml.etree.ElementTree as ET

from packerlabimaging.utils.utils import path_finder, points_in_circle_np


@dataclass
class naparm:
    naparm_path: str

    def __post_init__(self):
        self.target_areas = None
        self._photostimProcessing()

    def _parseNAPARMxml(self):

        print('\n\----- parsing Naparm xml file...')

        print('loading NAPARM_xml_path:')
        NAPARM_xml_path = path_finder(self.naparm_path, '.xml')[0]

        xml_tree = ET.parse(NAPARM_xml_path)
        root = xml_tree.getroot()

        title = root.get('Name')
        n_trials = int(root.get('Iterations'))

        inter_point_delay = 0
        spiral_duration = 20
        for elem in root:
            if int(elem[0].get('InterPointDelay')) > 0:
                inter_point_delay = int(elem[0].get('InterPointDelay'))
                spiral_duration = int(elem[0].get('Duration'))

        n_groups, n_reps, n_shots = [int(s) for s in re.findall(r'\d+', title)]

        print('\tNumbers of trials:', n_trials, '\n\tNumber of groups:', n_groups, '\n\tNumber of shots:', n_shots,
              '\n\tNumber of sequence reps:', n_reps, '\n\tInter-point delay:', inter_point_delay,
              '\n\tSpiral Duration (ms):', spiral_duration)

        # repetitions = int(root[1].get('Repetitions'))
        # print('Repetitions:', repetitions)

        self.n_groups = n_groups
        self.n_reps = n_reps
        self.n_shots = n_shots
        self.n_trials = n_trials
        self.inter_point_delay = inter_point_delay
        self.single_stim_dur = spiral_duration

    def _parseNAPARMgpl(self):
        spiral_size = -1

        print('\n\----- parsing Naparm gpl file...')

        NAPARM_gpl_path = path_finder(self.naparm_path, '.gpl')[0]
        print('loading NAPARM_gpl_path: ', NAPARM_gpl_path)

        xml_tree = ET.parse(NAPARM_gpl_path)
        root = xml_tree.getroot()

        for elem in root:
            if elem.get('Duration'):
                single_stim_dur = float(elem.get('Duration'))
                spiral_size = float(elem.get('SpiralSize'))
                print('Single stim dur (ms):', elem.get('Duration'))
                break

        for elem in root:
            if elem.get('SpiralSize'):
                spiral_size = float(elem.get('SpiralSize'))
                spiral_size = (spiral_size + 0.005155) / 0.005269  # hard-coded size of spiral from MATLAB code
                print('Spiral size .gpl file:', elem.get('SpiralSize'))
                print('spiral size (um):', int(spiral_size))
                break
        if spiral_size != -1:
            self.spiral_size = np.ceil(spiral_size)
        else:
            raise ValueError(f"could not set spiral_size, unable to get from NAPARM gpl")
        # self.single_stim_dur = single_stim_dur  # not sure why this was previously getting this value from here, but I'm now getting it from the xml file above

    def _photostimProcessing(
            self):  # TODO need to figure out how to handle different types of photostimulation experiment setups - think through in detail with Rob at a later date?

        """
        remember that this is currently configured for only the interleaved photostim, not the really fast photostim of multi-groups
        """

        self._parseNAPARMxml()
        self._parseNAPARMgpl()

        #         single_stim = self.single_stim_dur * self.n_shots
        #         total_single_stim = single_stim + self.inter_point_delay

        #         #self.stim_dur = total_single_stim - self.inter_point_delay

        #         total_multi_stim = total_single_stim * self.n_groups

        #         total_stim = total_multi_stim * self.n_reps

        #         self.stim_dur = total_stim - self.inter_point_delay

        ## PRAJAY EDITS HERE FOR PHOTOSTIMPROCESSING:
        # calculate duration (ms) of each individual trial (for a multi group stimulation trial)
        single_stim = self.single_stim_dur + self.inter_point_delay
        total_single_stim = single_stim * self.n_shots * self.n_groups * self.n_reps

        self.stim_dur = total_single_stim
        self.stim_freq = self.n_shots / self.stim_dur  # TODO Rob is this correct?


        print('Single stim. Duration (ms): ', self.single_stim_dur)
        print('Total stim. Duration per trial (ms): ', self.stim_dur)

    def _readTargetsImage(self):
        print('\n\t\----- Loading up target coordinates...')

        # load naparm targets file for this experiment
        naparm_path = os.path.join(self.naparm_path, 'Targets')

        listdir = os.listdir(naparm_path)


        ## All SLM targets
        target_image = np.empty(shape=[1, 1])
        for path in listdir:
            if 'AllFOVTargets' in path:
                target_file = path
                target_image = tf.imread(os.path.join(naparm_path, target_file))
        if target_image == np.empty(shape=[1, 1]):
            raise FileNotFoundError('AllFOVTargets image not found')

        targets = np.where(target_image > 0)
        
        # find targets by stim groups
        target_files = []
        for path in listdir:
            if 'FOVTargets_00' in path:
                target_files.append(path)

        counter = 0
        for slmgroup in target_files:
            target_image = tf.imread(os.path.join(naparm_path, slmgroup))
            targets = np.where(target_image > 0)
            targetCoordinates = list(zip(targets[1] * scale_factor, targets[0] * scale_factor))
            self.target_coords.append(targetCoordinates)
            self.n_targets += [len(targetCoordinates)]
            print('\t\tNumber of targets (in SLM group %s): ' % (counter + 1), len(targetCoordinates))
            counter += 1
        
        return targets


class Targets(naparm):
    def __init__(self, frame_x, frame_y, pix_sz_x, pix_sz_y):
        naparm.__init__()


    def _findTargetsAreas(self, frame_x, frame_y, pix_sz_x, pix_sz_y):  # TODO needs review by Rob - any bits of code missing under here? for e.g. maybe for target coords under multiple SLM groups?

        '''
        Finds cells that have been targeted for optogenetic photostimulation using Naparm in all-optical type experiments.
        output: coordinates of targets, and circular areas of targets
        Note this is not done by target groups however. So all of the targets are just in one big ls.
        '''

        targets = self._readTargetsImage()

        scale_factor = frame_x / 512  ## TODO need to get this from the NAPARM OUTPUT

        # targets_1 = np.where(target_image_scaled_1 > 0)
        # targets_2 = np.where(target_image_scaled_2 > 0)

        targetCoordinates = list(zip(targets[1] * scale_factor, targets[0] * scale_factor))
        print('\t\tNumber of targets:', len(targetCoordinates))

        # targetCoordinates_1 = ls(zip(targets_1[1], targets_1[0]))
        # print('Number of targets, SLM group #1:', len(targetCoordinates_1))
        #
        # targetCoordinates_2 = ls(zip(targets_2[1], targets_2[0]))
        # print('Number of targets, SLM group #2:', len(targetCoordinates_2))

        self.target_coords_all = targetCoordinates
        # self.target_coords_1 = targetCoordinates_1
        # self.target_coords_2 = targetCoordinates_2
        self.n_targets_total = len(targetCoordinates)

        ## specifying target areas in pixels to use for measuring responses of SLM targets
        radius_px = int(np.ceil(((self.spiral_size / 2) + 0) / pix_sz_x))
        print(f"\t\tspiral size: {self.spiral_size}um")
        print(f"\t\tpix sz x: {pix_sz_x}um")
        print("\t\tradius (in pixels): {:.2f}px".format(radius_px * pix_sz_x))

        for coord in targetCoordinates:
            target_area = ([item for item in points_in_circle_np(radius_px, x0=coord[0], y0=coord[1])])
            self.target_areas.append(target_area)

        ## target_areas that need to be excluded when filtering for nontarget cells
        radius_px = int(np.ceil(((self.spiral_size / 2) + 2.5) / pix_sz_x))
        print("\t\tradius of target exclusion zone (in pixels): {:.2f}px".format(radius_px * pix_sz_x))

        for coord in targetCoordinates:
            target_area = ([item for item in points_in_circle_np(radius_px, x0=coord[0], y0=coord[1])])
            self.target_areas_exclude.append(target_area)

        # # get areas for SLM group #1
        # target_areas_1 = []
        # for coord in targetCoordinates_1:
        #     target_area = ([item for item in pj.points_in_circle_np(radius, x0=coord[0], y0=coord[1])])
        #     target_areas_1.append(target_area)
        # self.target_areas_1 = target_areas_1
        #
        # # get areas for SLM group #2
        # target_areas_2 = []
        # for coord in targetCoordinates_2:
        #     target_area = ([item for item in pj.points_in_circle_np(radius, x0=coord[0], y0=coord[1])])
        #     target_areas_2.append(target_area)
        # self.target_areas_2 = target_areas_2



    def getTargetImage(obj):
        targ_img = np.zeros((obj.frame_x, obj.frame_y), dtype='uint8')
        # all FOV targets
        targ_areas = obj.target_areas
        for targ_area in targ_areas:
            for coord in targ_area:
                targ_img[coord[1], coord[0]] = 3000

        # targ_img_1 = np.zeros((obj.frame_x, obj.frame_y), dtype='uint8')
        # # FOV targets group #1
        # targ_areas = obj.target_areas_1
        # for targ_area in targ_areas:
        #     for coord in targ_area:
        #         targ_img_1[coord[1], coord[0]] = 3000
        #
        # targ_img_2 = np.zeros((obj.frame_x, obj.frame_y), dtype='uint8')
        # # FOV targets group #2
        # targ_areas = obj.target_areas_2
        # for targ_area in targ_areas:
        #     for coord in targ_area:
        #         targ_img_2[coord[1], coord[0]] = 3000

        return targ_img  # , targ_img_1, targ_img_2

    def _cropTargets(self, target_image_scaled):
        '''
        Crop the target image based on scan field centre and FOV size in pixels

        Inputs:
            target_image_scaled - 8-bit image with single pixels of 255 indicating target locations
        '''
        # make a bounding box centred on (512,512)
        x1 = 511 - self.imparams.frame_x / 2
        x2 = 511 + self.imparams.frame_x / 2
        y1 = 511 - self.imparams.frame_y / 2
        y2 = 511 + self.imparams.frame_y / 2

        # scan amplitudes in voltage from PrairieView software for 1x FOV
        ScanAmp_X = 2.62 * 2
        ScanAmp_Y = 2.84 * 2

        # scale scan amplitudes according to the imaging zoom
        ScanAmp_V_FOV_X = ScanAmp_X / self.imparams.zoom
        ScanAmp_V_FOV_Y = ScanAmp_Y / self.imparams.zoom

        # find the voltage per pixel
        scan_pix_y = ScanAmp_V_FOV_Y / 1024
        scan_pix_x = ScanAmp_V_FOV_X / 1024

        # find the offset (if any) of the scan field centre, in pixels
        offset_x = self.imparams.scan_x / scan_pix_x
        offset_y = self.imparams.scan_y / scan_pix_y

        # offset the bounding box by the appropriate number of pixels
        x1, x2, y1, y2 = round(x1 + offset_x), round(x2 + offset_x), round(y1 - offset_y), round(y2 - offset_y)

        # crop the target image using the offset bounding box to put the targets in imaging space
        return target_image_scaled[y1:y2, x1:x2]

    def _euclidDist(self, resp_positions):
        '''
        Find the mean Euclidean distance of all cells from a central point
        as a measure of spread

        Inputs:
            resp_positions - the median coordinates of cells
        Outputs:
            euclid_dist    - the mean Euclidean distance from central point
        '''
        # find centroid of all responding cells
        resp_coords = list(zip(*resp_positions))
        centroidx = np.mean(resp_coords[0])
        centroidy = np.mean(resp_coords[1])
        centroid = [centroidx, centroidy]

        # find distances of each cell from centroid
        dists = np.empty(resp_positions.shape[0])

        for i, resp in enumerate(resp_positions):
            dists[i] = np.linalg.norm(resp - centroid)

        euclid_dist = np.mean(dists)  # take average as a measure of spread

        return euclid_dist

    def _targetSpread(self):
        '''
        Find the mean Euclidean distance of responding targeted cells (trial-wise and trial average)
        '''
        # for each trial find targeted cells that responded
        trial_responders = self.trial_sig_dff[0]
        targeted_cells = np.repeat(self.targeted_cells[..., None],
                                   trial_responders.shape[1], 1)  # [..., None] is a quick way to expand_dims
        targeted_responders = targeted_cells & trial_responders

        cell_positions = np.array(self.cell_med[0])

        dists = np.empty(self.n_trials)

        # for each trial, find the spread of responding targeted cells
        for i, trial in enumerate(range(self.n_trials)):
            resp_cell = np.where(targeted_responders[:, trial])
            resp_positions = cell_positions[resp_cell]

            if resp_positions.shape[0] > 1:  # need more than 1 cell to measure spread...
                dists[i] = self._euclidDist(resp_positions)
            else:
                dists[i] = np.nan

        self.trial_euclid_dist = dists

        # find spread of targets that statistically significantly responded over all trials
        responder = self.sta_sig[0]
        targeted_responders = responder & self.targeted_cells

        resp_cell = np.where(targeted_responders)
        resp_positions = cell_positions[resp_cell]

        if resp_positions.shape[0] > 1:  # need more than 1 cell to measure spread...
            dist = self._euclidDist(resp_positions)
        else:
            dist = np.nan

        self.sta_euclid_dist = dist
