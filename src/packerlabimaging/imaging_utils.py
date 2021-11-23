### DRAFT SCRIPT FOR FUNCS

"""
This is the main collection of functions for processing and analysis of calcium imaging data in the Packer lab.

The fundamental object of data is one t-series of imaging, along with its associated .paq file, accessory files generated
from the microscope during data collection, and any user generated files associated with the experiment of this t-series.

"""

import re
import glob
import pandas as pd
import itertools

import os
import sys

from .utils_funcs import SaveDownsampledTiff, subselect_tiff, make_tiff_stack, convert_to_8bit, threshold_detect, s2p_loader  # listing functions from .utils_funcs that are used in this script

sys.path.append('/home/pshah/Documents/code/')
# from Vape.utils.paq2py import *
from matplotlib.colors import ColorConverter
from Vape.utils import STAMovieMaker_noGUI as STAMM
import scipy.stats as stats
import statsmodels.api
import statsmodels as sm
from suite2p.run_s2p import run_s2p
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import xml.etree.ElementTree as ET
import tifffile as tf
import csv
import warnings
import bisect

from utils.paq_utils import paq_read  ## TODO choose between adding paq_read to utils_funcs or make a new .py for paq_utils?
import pickle

from numba import njit


## CLASS DEFINITIONS
class TwoPhotonImaging:
    """Just two photon imaging related functions - currently set up for data collected from Bruker microscopes and
    suite2p processed Ca2+ imaging data """

    def __init__(self, tiff_location_dir: str, tiff_path: str, exp_metainfo: dict, analysis_save_path: str,
                 paq_path: str = None,
                 suite2p_path: str = None, make_downsampled_tiff: bool = False):
        """
        :param tiff_location_dir: parent directory where t-series .tiff is located (contains all of the accessory files for this t-series from the microscope)
        :param tiff_path: path to t-series .tiff
        :param exp_metainfo: dictionary containing any metainfo information field needed for this experiment. At minimum it needs to include prep #, t-series # and date of data collection.
        :param analysis_save_path: path of parent location where processed data objects will save by default
        :param paq_path: (optional) path to .paq file associated with current t-series
        :param suite2p_path: (optional) path to suite2p outputs folder associated with this t-series (plane0 file? or ops file? not sure yet)
        :param make_downsampled_tiff: flag to run generation and saving of downsampled tiff of t-series (saves to the analysis save location)
        """

        print('\n***** CREATING NEW TwoPhotonImaging with the following exp_metainfo: ', exp_metainfo)

        self.tiff_location_dir = tiff_location_dir  ## TODO add code this location from the tiff_path (don't require argument in __init__())
        self.tiff_path = tiff_path
        self.paq_path = paq_path
        self.metainfo = exp_metainfo
        self.analysis_save_path = analysis_save_path
        self.suite2p_path = suite2p_path

        ## TODO add checking path exists for all provided paths

        # create analysis save path location
        if os.path.exists(analysis_save_path):
            self.analysis_save_path = analysis_save_path
        else:
            self.analysis_save_path = analysis_save_path
            print('making analysis save folder at: \n  %s' % self.analysis_save_path)
            os.makedirs(self.analysis_save_path)

        # create pkl path and save expobj to pkl object
        pkl_path = "%s/%s_%s.pkl" % (self.analysis_save_path, exp_metainfo['date'],
                                     exp_metainfo['trial'])  # specify path in Analysis folder to save pkl object
        self.save_pkl(pkl_path=pkl_path)

        # if os.path.exists(self.analysis_save_path):
        #     pass
        # elif os.path.exists(self.analysis_save_path[:-17]):
        #     print('making analysis save folder at: \n  %s' % self.analysis_save_path)
        #     os.mkdir(self.analysis_save_path)
        # else:
        #     raise Exception('cannot find save folder path at: ', self.analysis_save_path[:-17])

        # elif os.path.exists(self.analysis_save_path[:-27]):
        #     print('making analysis save folder at: \n  %s \n and %s' % (self.analysis_save_path[:-17], self.analysis_save_path))
        #     os.mkdir(self.analysis_save_path[:-17])
        #     os.mkdir(self.analysis_save_path)

        self._parsePVMetadata()

        if make_downsampled_tiff:
            stack = self.mean_raw_flu_trace(save_pkl=True)
            SaveDownsampledTiff(stack=stack,
                                save_as=f"{analysis_save_path}/{exp_metainfo['date']}_{exp_metainfo['trial']}_downsampled.tif")

        if self.suite2p_path is not None:
            self.suite2p_completed = True
            self.s2pProcessing(s2p_path=self.suite2p_path)

        if self.paq_path is not None:
            self.paqProcessing(paq_path=self.paq_path)

        self.save_pkl(pkl_path=pkl_path)

    def s2pRun(self, ops, db, user_batch_size):

        '''run suite2p on the experiment object, using the attributes deteremined directly from the experiment object'''

        num_pixels = self.frame_x * self.frame_y
        sampling_rate = self.fps / self.n_planes
        diameter_x = 13 / self.pix_sz_x
        diameter_y = 13 / self.pix_sz_y
        diameter = int(diameter_x), int(diameter_y)
        batch_size = user_batch_size * (
                262144 / num_pixels)  # larger frames will be more RAM intensive, scale user batch size based on num pixels in 512x512 images

        if not db:
            db = {
                'data_path': [self.tiff_path],
                'fs': float(sampling_rate),
                'diameter': diameter,
                'batch_size': int(batch_size),
                'nimg_init': int(batch_size),
                'nplanes': self.n_planes
            }

        print(db)

        opsEnd = run_s2p(ops=ops, db=db)

        self.suite2p_completed = True  ## TODO does this make sense?

    def _getPVStateShard(self, path, key):

        '''
        Used in function PV metadata below
        '''

        value = []
        description = []
        index = []

        xml_tree = ET.parse(path)  # parse xml from a path
        root = xml_tree.getroot()  # make xml tree structure

        pv_state_shard = root.find('PVStateShard')  # find pv state shard element in root

        for elem in pv_state_shard:  # for each element in pv state shard, find the value for the specified key

            if elem.get('key') == key:

                if len(elem) == 0:  # if the element has only one subelement
                    value = elem.get('value')
                    break

                else:  # if the element has many subelements (i.e. lots of entries for that key)
                    for subelem in elem:
                        value.append(subelem.get('value'))
                        description.append(subelem.get('description'))
                        index.append(subelem.get('index'))
            else:
                for subelem in elem:  # if key not in element, try subelements
                    if subelem.get('key') == key:
                        value = elem.get('value')
                        break

            if value:  # if found key in subelement, break the loop
                break

        if not value:  # if no value found at all, raise exception
            raise Exception('ERROR: no element or subelement with that key')

        return value, description, index

    def _parsePVMetadata(self):

        print('\n-----parsing PV Metadata')

        tiff_path = self.tiff_location_dir
        path = []

        try:  # look for xml file in path, or two paths up (achieved by decreasing count in while loop)
            count = 2
            while count != 0 and not path:
                count -= 1
                for file in os.listdir(tiff_path):
                    if file.endswith('.xml'):
                        path = os.path.join(tiff_path, file)
                    if file.endswith('.env'):
                        env_path = os.path.join(tiff_path, file)
                tiff_path = os.path.dirname(tiff_path)

        except:
            raise Exception('ERROR: Could not find or load xml for this acquisition from %s' % tiff_path)

        xml_tree = ET.parse(path)  # parse xml from a path
        root = xml_tree.getroot()  # make xml tree structure

        sequence = root.find('Sequence')
        acq_type = sequence.get('type')

        if 'ZSeries' in acq_type:
            n_planes = len(sequence.findall('Frame'))

        else:
            n_planes = 1

        frame_period = float(self._getPVStateShard(path, 'framePeriod')[0])
        fps = 1 / frame_period

        frame_x = int(self._getPVStateShard(path, 'pixelsPerLine')[0])

        frame_y = int(self._getPVStateShard(path, 'linesPerFrame')[0])

        zoom = float(self._getPVStateShard(path, 'opticalZoom')[0])

        scanVolts, _, index = self._getPVStateShard(path, 'currentScanCenter')
        for scanVolts, index in zip(scanVolts, index):
            if index == 'XAxis':
                scan_x = float(scanVolts)
            if index == 'YAxis':
                scan_y = float(scanVolts)

        pixelSize, _, index = self._getPVStateShard(path, 'micronsPerPixel')
        for pixelSize, index in zip(pixelSize, index):
            if index == 'XAxis':
                pix_sz_x = float(pixelSize)
            if index == 'YAxis':
                pix_sz_y = float(pixelSize)

        env_tree = ET.parse(env_path)
        env_root = env_tree.getroot()

        elem_list = env_root.find('TSeries')
        # n_frames = elem_list[0].get('repetitions') # Rob would get the n_frames from env file
        # change this to getting the last actual index from the xml file

        n_frames = root.findall('Sequence/Frame')[-1].get('index')

        self.fps = fps
        self.frame_x = frame_x
        self.frame_y = frame_y
        self.n_planes = n_planes
        self.pix_sz_x = pix_sz_x
        self.pix_sz_y = pix_sz_y
        self.scan_x = scan_x
        self.scan_y = scan_y
        self.zoom = zoom
        self.n_frames = int(n_frames)

        print('n planes:', n_planes,
              '\nn frames:', int(n_frames),
              '\nfps:', fps,
              '\nframe size (px):', frame_x, 'x', frame_y,
              '\nzoom:', zoom,
              '\npixel size (um):', pix_sz_x, pix_sz_y,
              '\nscan centre (V):', scan_x, scan_y
              )

    def s2pProcessing(self, s2p_path: str, subtract_neuropil: bool = True, subset_frames: tuple = None,
                      save: bool = True):
        """processing of suite2p data from the current t-series
        :param s2p_path: path to the directory containing suite2p outputs
        :param subtract_neuropil: choose to subtract neuropil or not when loading s2p traces
        :param subset_frames: (optional) use to specifiy loading a subset of data from the overall s2p trace
        :param save: choose to save data object or not
        """

        if s2p_path is not None and s2p_path != self.suite2p_path:
            assert os.path.exists(s2p_path), print('ERROR: s2p path provided was not found')
            print(f"|- Updating suite2p path to newly provided path: {s2p_path}")
            self.suite2p_path = s2p_path  # update s2p path if provided different path
        self.s2p_subset_frames = subset_frames
        self.cell_id = []
        self.n_units = []
        self.cell_plane = []
        self.cell_med = []
        self.cell_x = []
        self.cell_y = []
        self.raw = []
        self.mean_img = []
        self.radius = []

        if self.n_planes == 1:
            FminusFneu, spks, self.stat, neuropil = s2p_loader(s2p_path,
                                                               subtract_neuropil)  # s2p_loader() is in Vape/utils_func
            ops = np.load(os.path.join(s2p_path, 'ops.npy'), allow_pickle=True).item()

            if self.s2p_subset_frames is None:
                self.raw = FminusFneu
                self.neuropil = neuropil
                self.spks = spks
            elif self.s2p_subset_frames is not None:
                self.raw = FminusFneu[:, self.s2p_subset_frames[0]:self.s2p_subset_frames[1]]
                self.spks = spks[:, self.s2p_subset_frames[0]:self.s2p_subset_frames[1]]
                self.neuropil = neuropil[:, self.s2p_subset_frames[0]:self.s2p_subset_frames[1]]
            self.mean_img = ops['meanImg']
            cell_id = []
            cell_med = []
            cell_x = []
            cell_y = []
            radius = []

            for cell, s in enumerate(self.stat):
                cell_id.append(s['original_index'])  # stat is an np array of dictionaries!
                cell_med.append(s['med'])
                cell_x.append(s['xpix'])
                cell_y.append(s['ypix'])
                radius.append(s['radius'])

            self.cell_id = cell_id
            self.n_units = len(self.cell_id)
            self.cell_med = cell_med
            self.cell_x = cell_x
            self.cell_y = cell_y
            self.radius = radius


        elif self.n_planes > 1:  ## TODO get Rob to review this loop, no experience with multi-plane analysis
            for plane in range(self.n_planes):
                # s2p_path = os.path.join(self.tiff_path, 'suite2p', 'plane' + str(plane))
                FminusFneu, self.spks, self.stat, self.neuro = s2p_loader(s2p_path,
                                                                          subtract_neuropil)  # s2p_loader() is in utils_func
                ops = np.load(os.path.join(s2p_path, 'ops.npy'), allow_pickle=True).item()

                self.raw.append(FminusFneu)
                self.mean_img.append(ops['meanImg'])
                cell_id = []
                cell_plane = []
                cell_med = []
                cell_x = []
                cell_y = []
                radius = []

                for cell, s in enumerate(self.stat):
                    cell_id.append(s['original_index'])  # stat is an np array of dictionaries!
                    cell_med.append(s['med'])
                    cell_x.append(s['xpix'])
                    cell_y.append(s['ypix'])
                    radius.append(s['radius'])

                self.cell_id.append(cell_id)
                self.n_units.append(len(self.cell_id[plane]))
                self.cell_med.append(cell_med)
                self.cell_x.append(cell_x)
                self.cell_y.append(cell_y)
                self.radius = radius

                cell_plane.extend([plane] * self.n_units[plane])
                self.cell_plane.append(cell_plane)

        if self.n_units > 1:
            print(f'|- Loaded {self.n_units} cells, recorded for {round(self.raw.shape[1] / self.fps, 2)} secs')
        else:
            print('******* SUITE2P DATA NOT LOADED.')

        self.save() if save else None

    def cellAreas(self, x=None, y=None):

        '''not sure what this function does'''  ## TODO Rob check

        self.cell_area = []

        if x:
            for i, _ in enumerate(self.cell_id):
                if self.cell_med[i][1] < x:
                    self.cell_area.append(0)
                else:
                    self.cell_area.append(1)

        if y:
            for i, _ in enumerate(self.cell_id):
                if self.cell_med[i][1] < y:
                    self.cell_area.append(0)
                else:
                    self.cell_area.append(1)

    def _good_cells(self, cell_ids: list, raws: np.ndarray, photostim_frames: list, std_thresh: int,
                    radiuses: list = None,
                    min_radius_pix: int = 2.5, max_radius_pix: int = 10, save=True):
        """
        This function is used for filtering for "good" cells based on detection of Flu deflections that are above some std threshold based on std_thresh.
        Note: a moving averaging window of 4 frames is used to find Flu deflections above std threshold.

        :param cell_ids: ls of cell ids to iterate over
        :param raws: raw flu values corresponding to the cell ls in cell_ids
        :param photostim_frames: frames to delete from the raw flu traces across all cells (because they are photostim frames)
        :param std_thresh: std. factor above which to define reliable Flu events
        :param radiuses: radiuses of the s2p ROI mask of all cells in the same order as cell_ids
        :param min_radius_pix:
        :param max_radius_pix:
        :return:
            good_cells: ls of cell_ids that had atleast 1 event above the std_thresh
            events_loc_cells: dictionary containing locations for each cell_id where the moving averaged Flu trace passes the std_thresh
            flu_events_cells: dictionary containing the dff Flu value corresponding to events_loc_cells
            stds = dictionary containing the dFF trace std value for each cell_id
        """

        print('\n----------------------------------------------------------------')
        print('running finding of good cells for suite2p ROIs ')
        print('----------------------------------------------------------------')

        good_cells = []
        events_loc_cells = {}
        flu_events_cells = {}
        stds = {}  # collect the std values for all filtered cells used later for identifying high and low std cells
        for i in range(len(cell_ids)):
            cell_id = cell_ids[i]

            if i % 100 == 0:  # print the progress once every 100 cell iterations
                print(i, " out of ", len(cell_ids), " cells done", end='\r')

            # print(i, " out of ", len(cell_ids), " cells")
            raw = raws[i]
            if np.mean(raw) > 1:  # exclude all negative raw traces and very small mean raw traces
                raw_ = np.delete(raw, photostim_frames)
                raw_dff = normalize_dff_jit(
                    raw_)  # note that this function is defined in this file a little further down
                std_ = raw_dff.std()

                raw_dff_ = moving_average(raw_dff, n=4)

                thr = np.mean(raw_dff) + std_thresh * std_
                events = np.where(raw_dff_ > thr)
                flu_values = raw_dff[events]

                if radiuses is not None:
                    radius = radiuses[i]
                    if len(events[0]) > 0 and radius > min_radius_pix and radius < max_radius_pix:
                        events_loc_cells[cell_id] = events
                        flu_events_cells[cell_id] = flu_values
                        stds[cell_id] = std_
                        good_cells.append(cell_id)
                elif len(events[0]) > 0:
                    events_loc_cells[cell_id] = events
                    flu_events_cells[cell_id] = flu_values
                    good_cells.append(cell_id)
                    stds[cell_id] = std_

                # if i == 465:  # just using this here if ever need to check back with specific cells if function seems to be misbehaving
                #     print(events, len(events[0]), thr)

        print('# of good cells found: ', len(good_cells), ' (out of ', len(cell_ids), ' ROIs)')
        self.good_cells = good_cells
        self.save() if save else None

        return good_cells, events_loc_cells, flu_events_cells, stds

    def paqProcessing(self, paq_path: str = None):
        """
        Loads .paq file and saves data from individual channels.

        :param paq_path: (optional) path to the .paq file for this data object
        """

        print('\n\n-----processing paq file...')

        if not hasattr(self, 'paq_path'):
            if paq_path is not None:
                self.paq_path = paq_path
            else:
                ValueError(
                    'ERROR: no paq_path defined for data object, please provide paq_path to load .paq file from.')
        elif paq_path is not None and paq_path != self.paq_path:
            assert os.path.exists(paq_path), print('ERROR: paq_path provided was not found')
            print(f"|- Updating paq_path to newly provided path: {paq_path}")
            self.paq_path = paq_path  # update paq_path if provided different path

        print('|- loading paq data from:', self.paq_path)

        paq, _ = paq_read(self.paq_path, plot=True)
        self.paq_rate = paq['rate']

        # find frame times
        clock_idx = paq['chan_names'].index('frame_clock')
        clock_voltage = paq['data'][clock_idx, :]

        frame_clock = threshold_detect(clock_voltage, 1)
        self.frame_clock = frame_clock
        plt.figure(figsize=(10, 5))
        plt.plot(clock_voltage)
        plt.plot(frame_clock, np.ones(len(frame_clock)), '.')
        plt.suptitle('frame clock from paq, with detected frame clock instances as scatter')
        sns.despine()
        plt.show()

        # find start and stop frame_clock times -- there might be multiple 2p imaging starts/stops in the paq trial (hence multiple frame start and end times)
        self.frame_start_times = [self.frame_clock[0]]  # initialize list
        self.frame_end_times = []
        i = len(self.frame_start_times)
        for idx in range(1, len(self.frame_clock) - 1):
            if (self.frame_clock[idx + 1] - self.frame_clock[idx]) > 2e3:
                i += 1
                self.frame_end_times.append(self.frame_clock[idx])
                self.frame_start_times.append(self.frame_clock[idx + 1])
        self.frame_end_times.append(self.frame_clock[-1])

        # handling cases where 2p imaging clock has been started/stopped >1 in the paq trial
        if len(self.frame_start_times) > 1:
            diff = [self.frame_end_times[idx] - self.frame_start_times[idx] for idx in
                    range(len(self.frame_start_times))]
            idx = diff.index(max(diff))
            self.frame_start_time_actual = self.frame_start_times[idx]
            self.frame_end_time_actual = self.frame_end_times[idx]
            self.frame_clock_actual = [frame for frame in self.frame_clock if
                                       self.frame_start_time_actual <= frame <= self.frame_end_time_actual]
        else:
            self.frame_start_time_actual = self.frame_start_times[0]
            self.frame_end_time_actual = self.frame_end_times[0]
            self.frame_clock_actual = self.frame_clock

    def mean_raw_flu_trace(self, plot: bool = False, save: bool = True):
        """
        Collects the raw mean of FOV fluorescence trace across the t-series.

        :param plot: (optional) plot mean fluorescence trace
        :param save: (optional) save data object after collecting mean fluorescence trace
        :return: mean fluorescence trace
        """
        print('\n-----collecting mean raw flu trace from tiff file...')
        print(f"|- loading raw TIFF file from: {self.tiff_path}")
        im_stack = tf.imread(self.tiff_path, key=range(self.n_frames))
        print('|- Loaded experiment tiff of shape: ', im_stack.shape)

        self.meanFluImg = np.mean(im_stack, axis=0)
        self.meanRawFluTrace = np.mean(np.mean(im_stack, axis=1), axis=1)

        self.save() if save else None

        aoplot.plotMeanRawFluTrace(expobj=self, stim_span_color=None, x_axis='frames', figsize=[20, 3],
                                   title='Mean raw Flu trace -') if plot else None

        return im_stack

    def plot_single_frame_tiff(self, frame_num: int = 0, title: str = None):
        """
        plots an image of a single specified tiff frame after reading using tifffile.
        :param frame_num: frame # from 2p imaging tiff to show (default is 0 - i.e. the first frame)
        :param title: (optional) give a string to use as title
        :return: matplotlib imshow plot
        """
        stack = tf.imread(self.tiff_path, key=frame_num)
        plt.imshow(stack, cmap='gray')
        plt.suptitle(title) if title is not None else plt.suptitle('frame num: %s' % frame_num)
        plt.show()
        return stack

    def stitch_reg_tiffs(self, first_frame: int, last_frame: int, force_crop: bool = False, force_stack: bool = False, s2p_run_batch: int = 2000):
        """
        Stitches together registered tiffs outputs from suite2p from the provided imaging frame start and end values.

        :param first_frame: first frame from the overall s2p run to start stitching from
        :param last_frame: last frame from the overall s2p run to end stitching at
        :param force_crop:
        :param force_stack:
        :param s2p_run_batch: batch size for suite2p run (defaults to 2000 because that is usually the batch size while running suite2p processing)
        """
        self.curr_trial_frames = [first_frame, last_frame]

        start = first_frame // s2p_run_batch
        end = last_frame // s2p_run_batch + 1

        tif_path_save = self.analysis_save_path + 'reg_tiff_%s.tif' % self.metainfo['trial']
        tif_path_save2 = self.analysis_save_path + 'reg_tiff_%s_r.tif' % self.metainfo['trial']
        reg_tif_folder = self.suite2p_path + '/reg_tif/'
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        sorted_paths = [reg_tif_folder + tif for tif in reg_tif_list][start:end + 1]

        print(tif_path_save)
        print(sorted_paths)

        if os.path.exists(tif_path_save):
            if force_stack:
                make_tiff_stack(sorted_paths, save_as=tif_path_save)
            else:
                pass
        else:
            make_tiff_stack(sorted_paths, save_as=tif_path_save)

        if not os.path.exists(tif_path_save2) or force_crop:
            with tf.TiffWriter(tif_path_save2, bigtiff=True) as tif:
                with tf.TiffFile(tif_path_save, multifile=False) as input_tif:
                    print('cropping registered tiff')
                    data = input_tif.asarray()
                    print('shape of stitched tiff: ', data.shape)
                reg_tif_crop = data[self.curr_trial_frames[0] - start * s2p_run_batch: self.curr_trial_frames[1] - (
                        self.curr_trial_frames[0] - start * s2p_run_batch)]
                print('saving cropped tiff ', reg_tif_crop.shape)
                tif.save(reg_tif_crop)

    def s2pMeanImage(self, s2p_path: str = None, plot: bool = True):
        """
        Return array of the s2p mean image.
        :param s2p_path: (optional) path to location of s2p data
        :param plot: (optional) option to plot the s2p mean image
        :return:
        """

        if s2p_path is None:
            if hasattr(self, 'suite2p_path'):
                s2p_path = self.suite2p_path
            else:
                ValueError(
                    'ERROR: no suite2p path defined for data object, please provide s2p_path to use for locating s2p data.')

        print(f'Plotting s2p mean image from {s2p_path}')

        os.chdir(s2p_path)

        ops = np.load('ops.npy', allow_pickle=True).item()

        mean_img = ops['meanImg']

        mean_img = np.array(mean_img, dtype='uint16')

        if plot:
            plt.imshow(mean_img, cmap='gray')
            plt.suptitle('s2p mean image')
            plt.show()

        return mean_img

    def save_pkl(self, pkl_path: str = None):
        if pkl_path is None:
            if hasattr(self, 'pkl_path'):
                pkl_path = self.pkl_path
            else:
                raise ValueError(
                    'pkl path for saving was not found in data object attributes, please provide pkl_path to save to')
        else:
            self.pkl_path = pkl_path

        with open(self.pkl_path, 'wb') as f:
            pickle.dump(self, f)
        print("\n\t -- data object saved to %s -- " % pkl_path)

    def save(self):
        self.save_pkl()


class WideFieldImaging:

    """
    WideField imaging data object.

    """

    def __init__(self):
        pass
