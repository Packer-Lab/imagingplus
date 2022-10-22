## main module for parsing, processing and interacting with cellsdata/files of 2p-stim protocols (i.e. NAPARM outputs)
from __future__ import absolute_import
from dataclasses import dataclass
from typing import TypedDict, Optional, MutableMapping
import re
import tifffile as tf
import os
import numpy as np
import xml.etree.ElementTree as ET

from imagingplus.utils.utils import path_finder, points_in_circle_np
from imagingplus.utils.images import ImportTiff


class naparm:

    def __init__(self, path: str):
        """
        Class for NAPARM photostimulation setup protocol.

        :param path: path to output from NAPARM used for photostimulation protocol of current imaging trial
        """
        self.path = path  #: path to output from NAPARM used for photostimulation protocol of current imaging trial
        self.stim_freq: float  #: frequency of photostim protocol (of a single photostim trial? or time between individual photostim trials?)
        self.single_stim_dur: float  #: duration of a single photostim shot (ms)
        self.inter_point_delay: float  #: duration of the delay between each photostim shot
        self.n_shots: int  #: num of photostim shots in a single photostim trial
        self.stim_duration_frames: int  #: num of imaging frames in a single photostim trial

    def __repr__(self):
        print(f"naparm analysis submodule. Loaded from: {self.path}")

    def _parseNAPARMxml(self):
        print('\n\----- parsing Naparm xml file...')

        print('loading NAPARM_xml_path:')
        NAPARM_xml_path = path_finder(self.path, '.xml')[0]

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

        NAPARM_gpl_path = path_finder(self.path, '.gpl')[0]
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

    def _NaparmProcessing(self):

        """
        Run processing of NAPARM files.
        Note: Currently configured for only the interleaved photostim, not other types of photostim of multi-groups
        """

        self._parseNAPARMxml()
        self._parseNAPARMgpl()

        #         single_stim = self.single_stim_dur * self.n_shots
        #         total_single_stim = single_stim + self.inter_point_delay

        #         #self.stim_dur = total_single_stim - self.inter_point_delay

        #         total_multi_stim = total_single_stim * self.n_groups

        #         total_stim = total_multi_stim * self.n_reps

        #         self.stim_dur = total_stim - self.inter_point_delay

        # calculate duration (ms) of each individual trial (for a multi group stimulation trial)
        single_stim = self.single_stim_dur + self.inter_point_delay
        total_single_stim = single_stim * self.n_shots * self.n_groups * self.n_reps

        self.stim_dur = total_single_stim
        self.stim_freq = self.n_shots / self.stim_dur

        print('Single stim. Duration (ms): ', self.single_stim_dur)
        print('Total stim. Duration per trial (ms): ', self.stim_dur)

    @classmethod
    def importNaparm(cls, path):
        """Alternative constructor to be used to run pipeline for loading and parsing NAPARM photostimulation protocol.

        :param path: path to directory containing NAPARM photostimulation protocol outputs.
        """
        npm = cls(path)
        npm._NaparmProcessing()
        return npm

class Targets(naparm):
    """Data and processing for photostimulation targets"""

    def __init__(self, naparm_path, frame_x, frame_y, pix_sz_x, pix_sz_y):

        # SLM target coords attr's
        # self.n_targets = []  # total number of target coordinates per SLM group
        # self.target_coords = []  # x, y locations of target coordinates per SLM group
        super().__init__(path=naparm_path)
        self._NaparmProcessing()
        self.stim_start_frames = None  #: list of frame numbers at onset of photostimulation
        self.photostim_frames = None  #: list of frame numbers contained during photostimulation trials
        self.target_areas = []  # photostim targeted pixels area coordinates per SLM group
        # self.target_coords_all: list = []  # all SLM target coordinates
        # self.n_targets_total: int = 0  # total number of SLM targets
        self.target_areas_exclude = []  # similar to .target_areas, but area diameter expanded (used in excluding cellsdata from this expanded region)

        self.__frame_x = frame_x
        self.__frame_y = frame_y
        self.__pix_sz_x = pix_sz_x
        self.__pix_sz_y = pix_sz_y

        self.target_coords_all, self.target_coords, self.n_targets, self.n_targets_total = self._readTargetsImage(
            self.__frame_x, self.__frame_y)
        self.target_areas, self.target_areas_exclude = self._findTargetsAreas(self.__pix_sz_x)

    def __repr__(self):
        print(f"naparm.Targets analysis submodule. Loaded from: {self.path}")

    def _readTargetsImage(self, frame_x, frame_y):
        scale_factor_x = frame_x / 512
        scale_factor_y = frame_y / 512

        print('\n\t\----- Loading up target coordinates...')
        scale_factor = frame_x / 512

        # load naparm targets file for this experiment
        naparm_path = os.path.join(self.path, 'Targets')

        listdir = os.listdir(naparm_path)

        ## All SLM targets
        target_file = ''
        for path in listdir:
            if 'AllFOVTargets' in path:
                target_file = path
        if target_file == '':
            raise FileNotFoundError('AllFOVTargets image not found')

        # target_image = tf.imread(os.path.join(naparm_path, target_file))
        target_image = ImportTiff(os.path.join(naparm_path, target_file))
        targets = np.where(target_image > 0)

        # find targets by stim groups
        target_files = []
        for path in listdir:
            if 'FOVTargets_00' in path:
                target_files.append(path)

        target_coords = []
        n_targets = []
        counter = 0
        for slmgroup in target_files:
            # target_image = tf.imread(os.path.join(naparm_path, slmgroup))
            target_image = ImportTiff((os.path.join(naparm_path, slmgroup)))
            targets = np.where(target_image > 0)
            targetCoordinates = list(zip(targets[1] * scale_factor_x, targets[0] * scale_factor_x))
            target_coords.append(targetCoordinates)
            n_targets.append(len(targetCoordinates))
            print('\t\tNumber of targets (in SLM group %s): ' % (counter + 1), len(targetCoordinates))
            counter += 1

        target_coords = np.asarray(target_coords)
        target_coords_all = target_coords.reshape(-1, target_coords.shape[-1])
        n_targets_total = len(target_coords_all)

        return target_coords_all, target_coords, n_targets, n_targets_total

    def _findTargetsAreas(self, pix_sz_x):

        """
        Finds cells that have been targeted for optogenetic photostimulation using Naparm in all-optical type experiments.
        output: coordinates of targets, and circular areas of targets
        Note this is not done by target groups however. So all of the targets are just in one big ls.

        :param pix_sz_x: size of pixel in x direction (assumes that length of x-direction == y-direction)
        :return:
        """

        print('\t\tNumber of targets (total):', len(self.target_coords_all))

        ## specifying target areas in pixels to use for measuring responses of SLM targets
        radius_px = int(np.ceil(((self.spiral_size / 2) + 0) / pix_sz_x))
        print(f"\t\tspiral size: {self.spiral_size}um")
        print(f"\t\tpix sz x: {pix_sz_x}um")
        print("\t\tradius (in pixels): {:.2f}px".format(radius_px * pix_sz_x))

        target_areas = []
        for coord in self.target_coords_all:
            target_area = ([item for item in points_in_circle_np(radius_px, x0=coord[0], y0=coord[1])])
            target_areas.append(target_area)

        ## target_areas that need to be excluded when filtering for nontarget cells
        radius_px = int(np.ceil(((self.spiral_size / 2) + 2.5) / pix_sz_x))
        print("\t\tradius of target exclusion zone (in pixels): {:.2f}px".format(radius_px * pix_sz_x))

        target_areas_exclude = []
        for coord in self.target_coords_all:
            target_area = ([item for item in points_in_circle_np(radius_px, x0=coord[0], y0=coord[1])])
            target_areas_exclude.append(target_area)

        return target_areas, target_areas_exclude

    def euclidDist(self, resp_positions):
        """
        Find the mean Euclidean distance of all cells from a central point
        as a measure of spread

        :param resp_positions: the median coordinates of cells
        :return: the mean Euclidean distance from central point
        """
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
        """
        Find the mean Euclidean distance of responding targeted cells (trial-wise and trial average)
        """
        # for each trial find targeted cells that responded
        trial_responders = self.trial_sig_dff[0]
        targeted_cells = np.repeat(self.targeted_cells[..., None],
                                   trial_responders.shape[1], 1)  # [..., None] is a quick way to expand_dims
        targeted_responders = targeted_cells & trial_responders

        cell_positions = np.array(self.Suite2p.cell_med[0])

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
