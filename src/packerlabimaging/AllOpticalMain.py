import os
import time
import re
import glob

import numpy as np
import pandas as pd
import scipy.stats as stats
import signal
import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api
import statsmodels as sm
import xml.etree.ElementTree as ET
import tifffile as tf

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)
from ._utils import convert_to_8bit, threshold_detect, path_finder, points_in_circle_np, normalize_dff, Utils
from ._paq import paq_read

from . TwoPhotonImagingMain import TwoPhotonImagingTrial
from . import _plotting, _anndata_funcs


# %%

PLANE = 0
BADFRAMESLOC = '/home/pshah/Documents/code/packerlabimaging/tests/'

class AllOpticalTrial(TwoPhotonImagingTrial):
    """This class provides methods for All Optical experiments"""

    def __init__(self, metainfo: dict, naparm_path: str, analysis_save_path: str, microscope: str,
                 analysisOptions: dict = {}, prestim_sec: float = 1.0, poststim_sec: float = 3.0, pre_stim_response_window: float = 0.500,
                 post_stim_response_window: float = 0.500, **kwargs):

        """
        :param tiff_path: path to t-series .tiff
        :param paq_path: path to .paq file associated with current t-series
        :param naparm_path: path to folder containing photostimulation setup built by NAPARM
        :param analysis_save_path: path of where to save the experiment analysis object
        :param metainfo: dictionary containing any metainfo information field needed for this experiment. At minimum it needs to include prep #, t-series # and date of data collection.
        :param microscope: name of microscope used to record imaging (options: "Bruker" (default), "Scientifica", "other")
        :param suite2p_path: (optional) path to suite2p outputs folder associated with this t-series (plane0 file? or ops file? not sure yet)
        :param make_downsampled_tiff: flag to run generation and saving of downsampled tiff of t-series (saves to the analysis save location)
        :param prestim_sec:
        :param poststim_sec:
        :param pre_stim_response_window:
        :param post_stim_response_window:
        :kwargs (optional):
        """

        print(f'\----- CREATING AllOpticalTrial data object for {metainfo["t series id"]}')

        # Initialize Attributes:

        # PHOTOSTIM PROTOCOL
        self.stim_channel = kwargs['stim_channel'] if 'stim_channel' in [
            *analysisOptions] else 'markpoints2packio'  # channel on paq file to read for determining stims

        self.photostim_frames = []  # imaging frames that are classified as overlapping with photostimulation
        self.stim_start_frames = []  # frame numbers from the start of each photostim trial
        self.stim_freq = None  # frequency of photostim protocol (of a single photostim trial? or time between individual photostim trials?)
        self.stim_dur: float = None  # total duration of a single photostim trial (units: msecs)
        self.single_stim_dur: float = None  # duration of a single photostim shot
        self.inter_point_delay: float = None  # duration of the delay between each photostim shot
        self.n_shots: int = None  # num of photostim shots in a single photostim trial
        self.stim_duration_frames: int = None  # num of imaging frames in a single photostim trial

        # paq file attr's
        self.paq_rate: int = None  # PackIO acquisition rate for .paq file
        self.frame_start_times: list = None  # paq clock timestamps of the first imaging acquisition frame of t-series
        self.frame_end_times: list = None  # paq clock timestamps of the last imaging acquisition frame of t-series

        # SLM target coords attr's
        self.n_targets = []  # total number of target coordinates per SLM group
        self.target_coords = []  # x, y locations of target coordinates per SLM group
        self.target_areas = []  # photostim targeted pixels area coordinates per SLM group
        self.target_coords_all: list = []  # all SLM target coordinates
        self.n_targets_total: int = 0  # total number of SLM targets
        self.target_areas_exclude = []  # similar to .target_areas, but area diameter expanded to exclude when filtering for nontarget cells

        # set photostim trial Flu trace snippet collection and response measurement analysis time windows
        self.prestim_sec = prestim_sec  # length of pre stim trace collected (in frames)
        self.poststim_sec = poststim_sec  # length of post stim trace collected (in frames)
        self.pre_stim_response_window_msec = pre_stim_response_window  # time window for collecting pre-stim measurement (units: msec)
        self.post_stim_response_window_msec = post_stim_response_window  # time window for collecting post-stim measurement (units: msec)

        # attr's for statistical analysis of photostim trials responses
        self.photostim_responses = None  # anndata object for collecting photostim responses and associated metadata for experiment and cells

        ######## TODO update comment descriptions
        self.all_trials = []  # all trials for each cell, dff detrended
        self.all_amplitudes = []  # all amplitudes of response between dff test periods

        self.stas = []  # avg of all trials for each cell, dff
        self.sta_amplitudes = []  # avg amplitude of response between dff test periods

        self.prob_response = None  # probability of response of cell to photostim trial; obtained from single trial significance (ROB suggestion)
        self.t_tests = []  # result from related samples t-test between dff test periods
        self.wilcoxons = []  # ROB to update
        self.trial_sig_dff = []  # based on dff increase above std of baseline
        self.trial_sig_dfsf = []  # based on df/std(f) increase in test period post-stim
        self.sta_sig = []  # based on t-test between dff test periods
        self.sta_sig_nomulti = []  # as above, no multiple comparisons correction
        ########

        #### initializing data processing, data analysis and/or results associated attr's
        self.n_trials = None  # number of photostimulation trials TODO change to assigning from array: cells x Flu frames x # of photostim trials

        ## PHOTOSTIM SLM TARGETS
        # TODO add attr's related to numpy array's and pandas dataframes for photostim trials - SLM targets
        self.responses_SLMtargets = []  # dF/prestimF responses for all SLM targets for each photostim trial
        self.responses_SLMtargets_tracedFF = []  # poststim dFF - prestim dFF responses for all SLM targets for each photostim trial

        ## ALL CELLS (from suite2p ROIs)
        # TODO add attr's related to numpy array's and pandas dataframes for photostim trials - suite2p ROIs
        self.photostimFluArray: np.ndarray = None  # array of dFF pre + post stim Flu snippets for each stim and cell [num cells x num peri-stim frames collected x num stims]  TODO need to write a quick func to collect over all planes to allow for multiplane imaging
        self.photostimResponseAmplitudes: pd.DataFrame = None  # photostim response amplitudes in a dataframe for each cell and each photostim

        if os.path.exists(naparm_path):
            self.__naparm_path = naparm_path
        else:
            raise FileNotFoundError(f"path not found, naparm_path: {naparm_path}")

        # initialize object as TwoPhotonImagingTrial
        TwoPhotonImagingTrial.__init__(self, metainfo=metainfo, analysis_save_path=analysis_save_path,
                                       microscope=microscope, **kwargs)

        # continue with photostimulation experiment processing
        self._stimProcessing(stim_channel=self.stim_channel)
        self._findTargetsAreas()
        self._find_photostim_add_bad_framesnpy()

        # extend annotated data object with imaging frames in photostim and stim_start_frames as additional keys in vars
        __frames_in_stim = [False] * self.ImagingParams.n_frames
        __stim_start_frame = [False] * self.ImagingParams.n_frames
        for frame in self.photostim_frames: __frames_in_stim[frame] = True
        for frame in self.stim_start_frames[PLANE]: __stim_start_frame[frame] = True
        self.data.add_variables(var_name='photostim_frame', values=__frames_in_stim)
        self.data.add_variables(var_name='stim_start_frame', values=__stim_start_frame)

        ##
        self.save()

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
            information = self.t_series_name
            return f"({information}) TwoPhotonImagingTrial.alloptical experimental trial object, last saved: {lastmod}"
        else:
            return f" -- unsaved TwoPhotonImagingTrial.alloptical experimental trial object -- "

    @property
    def naparm_path(self):
        """setting location of naparm files that specify the photostimulation experiment setup"""
        if self.__naparm_path[-1] == '/':
            return self.__naparm_path
        else:
            return self.__naparm_path + '/'

    @property
    def pre_stim_response_frames_window(self):
        """num frames for measuring Flu trace before each photostimulation trial during photostim response measurement (units: frames)"""
        return int(self.ImagingParams.fps * self.pre_stim_response_window_msec)

    @property
    def pre_stim_frames(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)"""
        return int(self.prestim_sec * self.ImagingParams.fps)

    @property
    def post_stim_frames(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)"""
        return int(self.poststim_sec * self.ImagingParams.fps)

    @property
    def post_stim_response_frames_window(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)"""
        return int(self.ImagingParams.fps * self.post_stim_response_window_msec)  # length of the post stim response test window (in frames)

    @property
    def timeVector(self):
        "vector of frame times in milliseconds rather than frames"
        n_frames = self.all_trials[0].shape[1]  ## <--- TODO double check this line returns as expected
        self.time = np.linspace(-self.prestim_sec, self.poststim_sec, self.ImagingParams.n_frames)

    def paqProcessing(self, paq_path: str = None, plot: bool = True, **kwargs):

        print('\n\----- processing paq file...')

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

        print(f'\tloading paq data from: {self.paq_path}')

        paq, _ = paq_read(self.paq_path, plot=plot)
        self.paq_rate = paq['rate']

        print(f"\t|- loaded {len(paq['chan_names'])} channels from .paq file: {paq['chan_names']}")

        # find frame times

        clock_idx = paq['chan_names'].index('frame_clock')
        clock_voltage = paq['data'][clock_idx, :]

        frame_clock = threshold_detect(clock_voltage, 1)
        self.__frame_clock = frame_clock
        if plot:
            plt.figure(figsize=(10, 5))
            plt.plot(clock_voltage)
            plt.plot(frame_clock, np.ones(len(frame_clock)), '.')
            plt.suptitle('frame clock from paq, with detected frame clock instances as scatter')
            sns.despine()
            plt.show()

        # find start and stop frame_clock times -- there might be multiple 2p imaging starts/stops in the paq trial (hence multiple frame start and end times)
        self.__frame_start_times = [self.__frame_clock[0]]  # initialize ls
        self.__frame_end_times = []
        i = len(self.__frame_start_times)
        for idx in range(1, len(self.__frame_clock) - 1):
            if (self.__frame_clock[idx + 1] - self.__frame_clock[idx]) > 2e3:
                i += 1
                self.__frame_end_times.append(self.__frame_clock[idx])
                self.__frame_start_times.append(self.__frame_clock[idx + 1])
        self.__frame_end_times.append(self.__frame_clock[-1])

        # handling cases where 2p imaging clock has been started/stopped >1 in the paq trial
        if len(self.__frame_start_times) > 1:
            diff = [self.__frame_end_times[idx] - self.__frame_start_times[idx] for idx in
                    range(len(self.__frame_start_times))]
            idx = diff.index(
                max(diff))  # choose the longest imaging sequence in the paq file as the actual frame clocks for the present trial's acquisition
            self.frame_start_times = self.__frame_start_times[idx]
            self.frame_end_times = self.__frame_end_times[idx]
            self.__frame_clock_actual = [frame for frame in self.__frame_clock if
                                self.frame_start_times <= frame <= self.frame_end_times]  #
        else:
            self.frame_start_times = self.__frame_start_times[0]
            self.frame_end_times = self.__frame_end_times[0]
            self.__frame_clock_actual = self.__frame_clock


        # find stim times
        stim_idx = paq['chan_names'].index(self.stim_channel)
        stim_volts = paq['data'][stim_idx, :]
        stim_times = threshold_detect(stim_volts, 1)
        # self.stim_times = stim_times
        self.stim_start_times = stim_times
        print('# of stims found on %s: %s' % (self.stim_channel, len(self.stim_start_times)))

        # correct this based on txt file
        duration_ms = self.stim_dur
        frame_rate = self.ImagingParams.fps / self.ImagingParams.n_planes
        duration_frames = np.ceil((duration_ms / 1000) * frame_rate)
        self.stim_duration_frames = int(duration_frames)

        if plot:
            plt.figure(figsize=(10, 5))
            plt.plot(stim_volts)
            plt.plot(stim_times, np.ones(len(stim_times)), '.')
            plt.suptitle('stim times')
            sns.despine()
            plt.show()

        # find stim frames

        for plane in range(self.ImagingParams.n_planes):  # TODO need to figure out how to handle multi-plane imaging and retrieving stim frame times
            stim_start_frames = []
            for stim in stim_times:
                # the index of the frame immediately preceeding stim
                stim_start_frame = next(
                    i - 1 for i, sample in enumerate(frame_clock[plane::self.ImagingParams.n_planes]) if sample - stim >= 0)
                stim_start_frames.append(stim_start_frame)

            self.stim_start_frames.append(stim_start_frames)

        # read in and save sparse version of all paq channels (only save data from timepoints at frame clock times)
        self.sparse_paq_data = {}
        for idx, chan in enumerate(self.paq_channels):
            self.sparse_paq_data[chan] = paq['data'][idx, self.__frame_clock_actual]


    ### ALLOPTICAL EXPERIMENT PHOTOSTIM PROTOCOL PROCESSING
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

        self.spiral_size = np.ceil(spiral_size)
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

        self.paqProcessing(plot=False)

    def _stimProcessing(self, stim_channel):
        self.stim_channel = stim_channel
        self._photostimProcessing()

    def _findTargetsAreas(
            self):  # TODO needs review by Rob - any bits of code missing under here? for e.g. maybe for target coords under multiple SLM groups?

        '''
        Finds cells that have been targeted for optogenetic photostimulation using Naparm in all-optical type experiments.
        output: coordinates of targets, and circular areas of targets
        Note this is not done by target groups however. So all of the targets are just in one big ls.
        '''

        print('\n\t\-----Loading up target coordinates...')

        # load naparm targets file for this experiment
        naparm_path = os.path.join(self.naparm_path, 'Targets')

        listdir = os.listdir(naparm_path)

        scale_factor = self.ImagingParams.frame_x / 512

        ## All SLM targets
        for path in listdir:
            if 'AllFOVTargets' in path:
                target_file = path
        target_image = tf.imread(os.path.join(naparm_path, target_file))

        targets = np.where(target_image > 0)
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
        radius_px = int(np.ceil(((self.spiral_size / 2) + 0) / self.ImagingParams.pix_sz_x))
        print(f"\t\tspiral size: {self.spiral_size}um")
        print(f"\t\tpix sz x: {self.ImagingParams.pix_sz_x}um")
        print("\t\tradius (in pixels): {:.2f}px".format(radius_px * self.ImagingParams.pix_sz_x))

        for coord in targetCoordinates:
            target_area = ([item for item in points_in_circle_np(radius_px, x0=coord[0], y0=coord[1])])
            self.target_areas.append(target_area)

        ## target_areas that need to be excluded when filtering for nontarget cells
        radius_px = int(np.ceil(((self.spiral_size / 2) + 2.5) / self.ImagingParams.pix_sz_x))
        print("\t\tradius of target exclusion zone (in pixels): {:.2f}px".format(radius_px * self.ImagingParams.pix_sz_x))

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

    def _findTargetedS2pROIs(self, plot: bool = True):
        """finding s2p cell ROIs that were also SLM targets (or more specifically within the target areas as specified by _findTargetAreas - include 15um radius from center coordinate of spiral)
        Make a binary mask of the targets and multiply by an image of the cells
        to find cells that were targeted

        --- LAST UPDATED NOV 6 2021 - copied over from Vape ---
        """

        print('\n\----- Searching for targeted cells in Suite2p ROIs... [Vape version]')

        ## TODO add necessary edits for multi-plane experiments

        ##### IDENTIFYING S2P ROIs THAT ARE WITHIN THE SLM TARGET SPIRAL AREAS
        # make all target area coords in to a binary mask
        targ_img = np.zeros([self.ImagingParams.frame_x, self.ImagingParams.frame_y], dtype='uint16')
        target_areas = np.array(self.target_areas)
        targ_img[target_areas[:, :, 1], target_areas[:, :, 0]] = 1

        # make an image of every cell area, filled with the index of that cell
        cell_img = np.zeros_like(targ_img)

        cell_y = np.array(self.Suite2p.suite2p_overall.cell_x)
        cell_x = np.array(self.Suite2p.suite2p_overall.cell_y)

        for i, coord in enumerate(zip(cell_x, cell_y)):
            cell_img[coord] = i + 1

        # binary mask x cell image to get the cells that overlap with target areas
        targ_cell = cell_img * targ_img

        targ_cell_ids = np.unique(targ_cell)[1:] - 1  # correct the cell id due to zero indexing
        self.targeted_cells = np.zeros([self.Suite2p.suite2p_overall.n_units], dtype='bool')
        self.targeted_cells[targ_cell_ids] = True
        # self.s2p_cell_targets = [self.cell_id[i] for i, x in enumerate(self.targeted_cells) if x is True]  # get ls of s2p cells that were photostim targetted
        self.s2p_cell_targets = [self.Suite2p.suite2p_overall.cell_id[i] for i in np.where(self.targeted_cells)[0]]  # get ls of s2p cells that were photostim targetted

        self.n_targeted_cells = np.sum(self.targeted_cells)

        print('\t|- Search completed.')
        self.save()
        print('\t|- Number of targeted cells: ', self.n_targeted_cells)

        ##### IDENTIFYING S2P ROIs THAT ARE WITHIN THE EXCLUSION ZONES OF THE SLM TARGETS
        # make all target area coords in to a binary mask
        targ_img = np.zeros([self.ImagingParams.frame_x, self.ImagingParams.frame_y], dtype='uint16')
        target_areas_exclude = np.array(self.target_areas_exclude)
        targ_img[target_areas_exclude[:, :, 1], target_areas_exclude[:, :, 0]] = 1

        # make an image of every cell area, filled with the index of that cell
        cell_img = np.zeros_like(targ_img)

        cell_y = np.array(self.Suite2p.suite2p_overall.cell_x)
        cell_x = np.array(self.Suite2p.suite2p_overall.cell_y)

        for i, coord in enumerate(zip(cell_x, cell_y)):
            cell_img[coord] = i + 1

        # binary mask x cell image to get the cells that overlap with target areas
        targ_cell = cell_img * targ_img

        targ_cell_ids = np.unique(targ_cell)[1:] - 1  # correct the cell id due to zero indexing
        self.exclude_cells = np.zeros([self.Suite2p.suite2p_overall.n_units], dtype='bool')
        self.exclude_cells[targ_cell_ids] = True
        self.s2p_cells_exclude = [self.Suite2p.suite2p_overall.cell_id[i] for i in
                                  np.where(self.exclude_cells)[0]]  # get ls of s2p cells that were photostim targetted

        self.n_exclude_cells = np.sum(self.exclude_cells)

        print('\t|-Search completed.')
        self.save()
        print(f"\t|-Number of exclude Suite2p ROIs: {self.n_exclude_cells}")

        # define non targets from suite2p ROIs (exclude cells in the SLM targets exclusion - .s2p_cells_exclude)
        self.s2p_nontargets = [cell for cell in self.Suite2p.suite2p_overall.cell_id if
                               cell not in self.s2p_cells_exclude]  ## exclusion of cells that are classified as s2p_cell_targets

        print(f"\t|-Number of good, s2p non-target ROIs: {len(self.s2p_nontargets)}")

        if plot:
            fig, ax = plt.subplots(figsize=[6, 6])

            targ_img = np.zeros([self.ImagingParams.frame_x, self.ImagingParams.frame_y], dtype='uint16')
            target_areas = np.array(self.target_areas)
            targ_img[target_areas[:, :, 1], target_areas[:, :, 0]] = 1
            ax.imshow(targ_img, cmap='Greys_r', zorder=0)
            ax.set_title('SLM targets areas')
            # for (x, y) in self.target_coords_all:
            #     ax.scatter(x=x, y=y, edgecolors='white', facecolors='none', linewidths=1.0)
            fig.show()


        # add targets classification as observations annotation to .data anndata
        self.data = Utils.add_observation(self.data, 'photostim_target', values=list(self.targeted_cells))



    def _find_photostim_add_bad_framesnpy(self):
        """finds all photostim frames and saves them into the bad_frames attribute for the exp object and saves bad_frames.npy"""
        print('\n\-----Finding photostimulation frames in imaging frames ...')
        print('# of photostim frames calculated per stim. trial: ', self.stim_duration_frames + 1)

        for j in self.stim_start_frames[PLANE]:
            for i in range(
                    self.stim_duration_frames + 1):  # usually need to remove 1 more frame than the stim duration, as the stim isn't perfectly aligned with the start of the imaging frame
                self.photostim_frames.append(j + i)

        # print(photostim_frames)
        print('\t|- Original # of frames:', self.ImagingParams.n_frames, 'frames')
        print('\t|- # of Photostim frames:', len(self.photostim_frames), 'frames')
        print('\t|- Minus photostim. frames total:', self.ImagingParams.n_frames - len(self.photostim_frames), 'frames')

        if len(self.photostim_frames) > 0:
            print(
                # f'***Saving a total of {len(self.photostim_frames)} photostim frames to bad_frames.npy at: {self.tiff_path_dir}/bad_frames.npy')
                f'***Saving a total of {len(self.photostim_frames)} photostim frames to bad_frames.npy at: {BADFRAMESLOC}/bad_frames.npy')  # TODO replace BADFRAMESLOC with self.tiff_path_dir
            np.save(f'{self.tiff_path_dir}/bad_frames.npy',
                    self.photostim_frames)  # save to npy file and remember to move npy file to tiff folder before running with suite2p

    #### TODO review attr's and write docs from the following functions: // start
    # used for creating tiffs that remove artifacts from alloptical experiments with photostim artifacts
    def rm_artifacts_tiffs(expobj, tiffs_loc, new_tiffs):
        ### make a new tiff file (not for suite2p) with the first photostim frame whitened, and save new tiff
        print('\n-----making processed photostim .tiff from:')
        tiff_path = tiffs_loc
        print(tiff_path)
        im_stack = tf.imread(tiff_path, key=range(expobj.n_frames))
        print('Processing experiment tiff of shape: ', im_stack.shape)

        frames_to_whiten = []
        for j in expobj.stim_start_frames:
            frames_to_whiten.append(j)

        # number of photostim frames with artifacts
        frames_to_remove = []
        for j in expobj.stim_start_frames:
            for i in range(0,
                           expobj.stim_duration_frames + 1):  # usually need to remove 1 more frame than the stim duration, as the stim isn't perfectly aligned with the start of the imaging frame
                frames_to_remove.append(j + i)

        print('# of total photostim artifact frames:', len(frames_to_remove))

        im_stack_1 = im_stack
        a = np.full_like(im_stack_1[0], fill_value=0)
        a[0:100, 0:100] = 5000.
        for frame in frames_to_whiten:
            im_stack_1[frame - 3] = im_stack_1[frame - 3] + a
            im_stack_1[frame - 2] = im_stack_1[frame - 2] + a
            im_stack_1[frame - 1] = im_stack_1[frame - 1] + a
        print('Shape', im_stack_1.shape)

        im_stack_1 = np.delete(im_stack_1, frames_to_remove, axis=0)
        print('After delete shape artifactrem', im_stack_1.shape)

        save_path = (new_tiffs + "_artifactrm.tif")
        tf.imwrite(save_path, im_stack_1, photometric='minisblack')

        del im_stack_1

        # draw areas on top of im_stack_1 where targets are:
        im_stack_2 = im_stack
        print('Shape', im_stack_2.shape)

        for stim in range(expobj.n_groups):
            b = np.full_like(im_stack_2[0], fill_value=0)
            targets = expobj.target_areas[stim]
            for i in np.arange(len(targets)):
                for j in targets[i]:
                    b[j] = 5000

            all_stim_start_frames = []
            for stim_frame in expobj.stim_start_frames[stim::expobj.n_groups]:
                all_stim_start_frames.append(stim_frame)
            for frame in all_stim_start_frames:
                #         im_stack_2[frame-4] = im_stack_2[frame-4]+b
                #         im_stack_2[frame-3] = im_stack_2[frame-3]+b
                #        im_stack_2[frame-2] = im_stack_2[frame-2]+b
                im_stack_2[frame - 1] = im_stack_2[frame - 1] + b

        im_stack_2 = np.delete(im_stack_2, expobj.photostim_frames, axis=0)

        print('After delete shape targetcells', im_stack_2.shape)

        save_path = (new_tiffs + '_targetcells.tif')
        tf.imwrite(save_path, im_stack_2, photometric='minisblack')

        print('done saving to: ', save_path)

        del im_stack_2
        del im_stack

    def s2pMasks(expobj, s2p_path, cell_ids):
        os.chdir(s2p_path)
        stat = np.load('stat.npy', allow_pickle=True)
        ops = np.load('ops.npy', allow_pickle=True).item()
        iscell = np.load('iscell.npy', allow_pickle=True)
        mask_img = np.zeros((ops['Ly'], ops['Lx']), dtype='uint8')
        for n in range(0, len(iscell)):
            if n in cell_ids:
                ypix = stat[n]['ypix']
                xpix = stat[n]['xpix']
                mask_img[ypix, xpix] = 3000

        # s2p targets - all SLM targets
        targets_s2p_img = np.zeros((ops['Ly'], ops['Lx']), dtype='uint8')
        for n in range(0, len(iscell)):
            if n in expobj.s2p_cell_targets:
                ypix = stat[n]['ypix']
                xpix = stat[n]['xpix']
                targets_s2p_img[ypix, xpix] = 3000

        # # s2p targets - SLM group #1 targets
        # targets_s2p_img_1 = np.zeros((ops['Ly'], ops['Lx']), dtype='uint8')
        # for n in range(0, len(iscell)):
        #     if n in obj.s2p_cell_targets_1:
        #         ypix = stat[n]['ypix']
        #         xpix = stat[n]['xpix']
        #         targets_s2p_img_1[ypix, xpix] = 3000
        #
        # # s2p targets - SLM group #2 targets
        # targets_s2p_img_2 = np.zeros((ops['Ly'], ops['Lx']), dtype='uint8')
        # for n in range(0, len(iscell)):
        #     if n in obj.s2p_cell_targets_2:
        #         ypix = stat[n]['ypix']
        #         xpix = stat[n]['xpix']
        #         targets_s2p_img_2[ypix, xpix] = 3000

        return mask_img, targets_s2p_img,  # targets_s2p_img_1, targets_s2p_img_2

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

    def s2pMaskStack(obj, pkl_list, s2p_path, parent_folder, force_redo: bool = False):
        """makes a TIFF stack with the s2p mean image, and then suite2p ROI masks for all cells detected, target cells, and also SLM targets as well?"""

        for pkl in pkl_list:
            expobj = obj

            print('Retrieving s2p masks for:', pkl, end='\r')

            # with open(pkl, 'rb') as f:
            #     expobj = pickle.load(f)

            # ls of cell ids to filter s2p masks by
            # cell_id_list = [ls(range(1, 99999)),  # all
            #                 expobj.photostim_r.cell_id[0],  # cells
            #                 [expobj.photostim_r.cell_id[0][i] for i, b in enumerate(expobj.photostim_r.cell_s1[0]) if
            #                  b == False],  # s2 cells
            #                 [expobj.photostim_r.cell_id[0][i] for i, b in enumerate(expobj.photostim_r.is_target) if
            #                  b == 1],  # pr cells
            #                 [expobj.photostim_s.cell_id[0][i] for i, b in enumerate(expobj.photostim_s.is_target) if
            #                  b == 1],  # ps cells
            #                 ]
            #
            cell_ids = expobj.cell_id

            # empty stack to fill with images
            stack = np.empty((0, expobj.frame_y, expobj.frame_x), dtype='uint8')

            s2p_path = s2p_path

            # mean image from s2p
            mean_img = obj.s2pMeanImage(s2p_path)
            mean_img = np.expand_dims(mean_img, axis=0)
            stack = np.append(stack, mean_img, axis=0)

            # mask images from s2p
            mask_img, targets_s2p_img = obj.s2pMasks(s2p_path=s2p_path, cell_ids=cell_ids)
            mask_img = np.expand_dims(mask_img, axis=0)
            targets_s2p_img = np.expand_dims(targets_s2p_img, axis=0)
            # targets_s2p_img_1 = np.expand_dims(targets_s2p_img_1, axis=0)
            # targets_s2p_img_2 = np.expand_dims(targets_s2p_img_2, axis=0)
            stack = np.append(stack, mask_img, axis=0)
            stack = np.append(stack, targets_s2p_img, axis=0)
            # stack = np.append(stack, targets_s2p_img_1, axis=0)
            # stack = np.append(stack, targets_s2p_img_2, axis=0)

            # # sta images
            # for file in os.listdir(stam_save_path):
            #     if all(s in file for s in ['AvgImage', expobj.photostim_r.tiff_path.split('/')[-1]]):
            #         pr_sta_img = tf.imread(os.path.join(stam_save_path, file))
            #         pr_sta_img = np.expand_dims(pr_sta_img, axis=0)
            #     elif all(s in file for s in ['AvgImage', expobj.photostim_s.tiff_path.split('/')[-1]]):
            #         ps_sta_img = tf.imread(os.path.join(stam_save_path, file))
            #         ps_sta_img = np.expand_dims(ps_sta_img, axis=0)

            # stack = np.append(stack, pr_sta_img, axis=0)
            # stack = np.append(stack, ps_sta_img, axis=0)

            # target images
            targ_img = obj.getTargetImage()
            targ_img = np.expand_dims(targ_img, axis=0)
            stack = np.append(stack, targ_img, axis=0)

            # targ_img_1 = np.expand_dims(targ_img_1, axis=0)
            # stack = np.append(stack, targ_img_1, axis=0)
            #
            # targ_img_2 = np.expand_dims(targ_img_2, axis=0)
            # stack = np.append(stack, targ_img_2, axis=0)

            # stack is now: mean_img, all_rois, all_cells, s2_cells, pr_cells, ps_cells,
            # (whisker,) pr_sta_img, ps_sta_img, pr_target_areas, ps_target_areas
            # c, x, y = stack.shape
            # stack.shape = 1, 1, c, x, y, 1  # dimensions in TZCYXS order

            x_pix = expobj.pix_sz_x
            y_pix = expobj.pix_sz_y

            save_path = os.path.join(parent_folder, pkl.split('/')[-1][:-4] + '_s2p_masks.tif')

            tf.imwrite(save_path, stack, photometric='minisblack')
            print('\ns2p ROI + photostim targets masks saved in TIFF to: ', save_path)

    def _cropTargets(self, target_image_scaled):
        '''
        Crop the target image based on scan field centre and FOV size in pixels

        Inputs:
            target_image_scaled - 8-bit image with single pixels of 255 indicating target locations
        '''
        # make a bounding box centred on (512,512)
        x1 = 511 - self.ImagingParams.frame_x / 2
        x2 = 511 + self.ImagingParams.frame_x / 2
        y1 = 511 - self.ImagingParams.frame_y / 2
        y2 = 511 + self.ImagingParams.frame_y / 2

        # scan amplitudes in voltage from PrairieView software for 1x FOV
        ScanAmp_X = 2.62 * 2
        ScanAmp_Y = 2.84 * 2

        # scale scan amplitudes according to the imaging zoom
        ScanAmp_V_FOV_X = ScanAmp_X / self.ImagingParams.zoom
        ScanAmp_V_FOV_Y = ScanAmp_Y / self.ImagingParams.zoom

        # find the voltage per pixel
        scan_pix_y = ScanAmp_V_FOV_Y / 1024
        scan_pix_x = ScanAmp_V_FOV_X / 1024

        # find the offset (if any) of the scan field centre, in pixels
        offset_x = self.ImagingParams.scan_x / scan_pix_x
        offset_y = self.ImagingParams.scan_y / scan_pix_y

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

    def STAProcessing(self, plane):
        '''
        Make stimulus triggered average (STA) traces and calculate the STA amplitude
        of response

        Input:
            plane - imaging plane n
        '''
        # make stas, [plane x cell x frame]
        stas = np.mean(self.all_trials[plane], axis=2)
        self.stas.append(stas)

        # make sta amplitudes, [plane x cell]
        pre_sta = np.mean(stas[:, self.pre_trial_frames], axis=1)
        post_sta = np.mean(stas[:, self.post_trial_frames], axis=1)
        sta_amplitudes = post_sta - pre_sta
        self.sta_amplitudes.append(sta_amplitudes)

    ####                                                                    // end

    ### ALLOPTICAL PROCESSING OF TRACES
    ## ... no methods determined here yet...

    ### ALLOPTICAL ANALYSIS - FOCUS ON SLM TARGETS RELATED METHODS TODO need to review this whole section
    def collect_traces_from_targets(self, curr_trial_frames: list, reg_tif_folder: str, save: bool = True):
        """uses registered tiffs to collect raw traces from SLM target areas"""

        if reg_tif_folder is None:
            if self.suite2p_path:
                reg_tif_folder = self.suite2p_path + '/reg_tif/'
                print(f"\- trying to load registerred tiffs from: {reg_tif_folder}")
        else:
            raise Exception(f"Must provide reg_tif_folder path for loading registered tiffs")
        if not os.path.exists(reg_tif_folder):
            raise Exception(f"no registered tiffs found at path: {reg_tif_folder}")

        print(
            f'\n\ncollecting raw Flu traces from SLM target coord. areas from registered TIFFs from: {reg_tif_folder}')
        # read in registered tiff
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        start = curr_trial_frames[0] // 2000  # 2000 because that is the batch size for suite2p run
        end = curr_trial_frames[1] // 2000 + 1

        mean_img_stack = np.zeros([end - start, self.ImagingParams.frame_x, self.ImagingParams.frame_y])
        # collect mean traces from target areas of each target coordinate by reading in individual registered tiffs that contain frames for current trial
        targets_trace_full = np.zeros([len(self.target_coords_all), (end - start) * 2000], dtype='float32')
        counter = 0
        for i in range(start, end):
            tif_path_save2 = self.suite2p_path + '/reg_tif/' + reg_tif_list[i]
            with tf.TiffFile(tif_path_save2, multifile=False) as input_tif:
                print('|- reading tiff: %s' % tif_path_save2)
                data = input_tif.asarray()

            targets_trace = np.zeros([len(self.target_coords_all), data.shape[0]], dtype='float32')
            for coord in range(len(self.target_coords_all)):
                target_areas = np.array(
                    self.target_areas)  # TODO update this so that it doesn't include the extra exclusion zone
                x = data[:, target_areas[coord, :, 1], target_areas[coord, :, 0]]  # = 1
                targets_trace[coord] = np.mean(x, axis=1)

            targets_trace_full[:, (i - start) * 2000: ((i - start) * 2000) + data.shape[
                0]] = targets_trace  # iteratively write to each successive segment of the targets_trace array based on the length of the reg_tiff that is read in.

            mean_img_stack[counter] = np.mean(data, axis=0)
            counter += 1

        # final part, crop to the *exact* frames for current trial
        self.raw_SLMTargets = targets_trace_full[:,
                              curr_trial_frames[0] - start * 2000: curr_trial_frames[1] - (start * 2000)]

        self.dFF_SLMTargets = normalize_dff(self.raw_SLMTargets, threshold_pct=10)

        self.meanFluImg_registered = np.mean(mean_img_stack, axis=0)

        self.save() if save else None

    def get_alltargets_stim_traces_norm(self, process: str, targets_idx: int = None, subselect_cells: list = None,
                                        pre_stim=15, post_stim=200, filter_sz: bool = False, stims: list = None):
        """
        primary function to measure the dFF and dF/setdF trace SNIPPETS for photostimulated targets.
        :param stims:
        :param targets_idx: integer for the index of target cell to process
        :param subselect_cells: ls of cells to subset from the overall set of traces (use in place of targets_idx if desired)
        :param pre_stim: number of frames to use as pre-stim
        :param post_stim: number of frames to use as post-stim
        :param filter_sz: whether to filter out stims that are occuring seizures
        :return: lists of individual targets dFF traces, and averaged targets dFF over all stims for each target
        """
        if filter_sz:
            print('\n -- filter_sz active --')

        if stims is None:
            stim_timings = self.stim_start_frames
        else:
            stim_timings = stims

        if process == 'trace raw':  ## specify which data to process (i.e. do you want to process whole trace dFF traces?)
            data_to_process = self.raw_SLMTargets
        elif process == 'trace dFF':
            data_to_process = self.dFF_SLMTargets
        else:
            ValueError('need to provide `process` as either `trace raw` or `trace dFF`')

        if subselect_cells:
            num_cells = len(data_to_process[subselect_cells])
            targets_trace = data_to_process[subselect_cells]  ## NOTE USING .raw traces
        else:
            num_cells = len(data_to_process)
            targets_trace = data_to_process

        # collect photostim timed average dff traces of photostim targets
        targets_dff = np.zeros(
            [num_cells, len(self.stim_start_frames), pre_stim + self.stim_duration_frames + post_stim])
        # SLMTargets_stims_dffAvg = np.zeros([num_cells, pre_stim_sec + post_stim_sec])

        targets_dfstdF = np.zeros(
            [num_cells, len(self.stim_start_frames), pre_stim + self.stim_duration_frames + post_stim])
        # targets_dfstdF_avg = np.zeros([num_cells, pre_stim_sec + post_stim_sec])

        targets_raw = np.zeros(
            [num_cells, len(self.stim_start_frames), pre_stim + self.stim_duration_frames + post_stim])
        # targets_raw_avg = np.zeros([num_cells, pre_stim_sec + post_stim_sec])

        if targets_idx is not None:
            print('collecting stim traces for cell ', targets_idx + 1)
            if filter_sz:
                flu = [targets_trace[targets_idx][stim - pre_stim: stim + self.stim_duration_frames + post_stim] for
                       stim in
                       stim_timings if
                       stim not in self.seizure_frames]
            else:
                flu = [targets_trace[targets_idx][stim - pre_stim: stim + self.stim_duration_frames + post_stim] for
                       stim in
                       stim_timings]
            # flu_dfstdF = []
            # flu_dff = []
            for i in range(len(flu)):
                trace = flu[i]
                mean_pre = np.mean(trace[0:pre_stim])
                if process == 'trace raw':
                    trace_dff = ((trace - mean_pre) / mean_pre) * 100
                elif process == 'trace dFF':
                    trace_dff = (trace - mean_pre)
                else:
                    ValueError('not sure how to calculate peri-stim traces...')
                std_pre = np.std(trace[0:pre_stim])
                dFstdF = (trace - mean_pre) / std_pre  # make dF divided by std of pre-stim F trace

                targets_raw[targets_idx, i] = trace
                targets_dff[targets_idx, i] = trace_dff
                targets_dfstdF[targets_idx, i] = dFstdF
            print(f"shape of targets_dff[targets_idx]: {targets_dff[targets_idx].shape}")
            return targets_raw[targets_idx], targets_dff[targets_idx], targets_dfstdF[targets_idx]

        else:
            for cell_idx in range(num_cells):
                print('collecting stim traces for cell %s' % subselect_cells[cell_idx]) if subselect_cells else None

                if filter_sz:
                    if hasattr(self, 'slmtargets_szboundary_stim') and self.slmtargets_szboundary_stim is not None:
                        flu = []
                        for stim in stim_timings:
                            if stim in self.slmtargets_szboundary_stim.keys():  # some stims dont have sz boundaries because of issues with their TIFFs not being made properly (not readable in Fiji), usually it is the first TIFF in a seizure
                                if cell_idx not in self.slmtargets_szboundary_stim[stim]:
                                    flu.append(targets_trace[cell_idx][
                                               stim - pre_stim: stim + self.stim_duration_frames + post_stim])
                    else:
                        flu = []
                        print('classifying of sz boundaries not completed for this expobj',
                              self.metainfo['animal prep.'], self.metainfo['trial'])
                    # flu = [targets_trace[cell_idx][stim - pre_stim_sec: stim + self.stim_duration_frames + post_stim_sec] for
                    #        stim
                    #        in stim_timings if
                    #        stim not in self.seizure_frames]
                else:
                    flu = [targets_trace[cell_idx][stim - pre_stim: stim + self.stim_duration_frames + post_stim] for
                           stim in stim_timings]

                # flu_dfstdF = []
                # flu_dff = []
                # flu = []
                if len(flu) > 0:
                    for i in range(len(flu)):
                        trace = flu[i]
                        mean_pre = np.mean(trace[0:pre_stim])
                        trace_dff = ((trace - mean_pre) / mean_pre) * 100
                        std_pre = np.std(trace[0:pre_stim])
                        dFstdF = (trace - mean_pre) / std_pre  # make dF divided by std of pre-stim F trace

                        targets_raw[cell_idx, i] = trace
                        targets_dff[cell_idx, i] = trace_dff
                        targets_dfstdF[cell_idx, i] = dFstdF
                        # flu_dfstdF.append(dFstdF)
                        # flu_dff.append(trace_dff)

                # targets_dff.append(flu_dff)  # contains all individual dFF traces for all stim times
                # SLMTargets_stims_dffAvg.append(np.nanmean(flu_dff, axis=0))  # contains the dFF trace averaged across all stim times

                # targets_dfstdF.append(flu_dfstdF)
                # targets_dfstdF_avg.append(np.nanmean(flu_dfstdF, axis=0))

                # SLMTargets_stims_raw.append(flu)
                # targets_raw_avg.append(np.nanmean(flu, axis=0))

            targets_dff_avg = np.mean(targets_dff, axis=1)
            targets_dfstdF_avg = np.mean(targets_dfstdF, axis=1)
            targets_raw_avg = np.mean(targets_raw, axis=1)

            print(f"shape of targets_dff_avg: {targets_dff_avg.shape}")
            return targets_dff, targets_dff_avg, targets_dfstdF, targets_dfstdF_avg, targets_raw, targets_raw_avg

    # calculate reliability of photostim responsiveness of all of the targeted cells (found in s2p output)
    def get_SLMTarget_responses_dff(self, process: str, threshold=10, stims_to_use: list = None):
        """
        calculations of dFF responses to photostimulation of SLM Targets. Includes calculating reliability of slm targets,
        saving success stim locations, and saving stim response magnitudes as pandas dataframe.

        :param threshold: dFF threshold above which a response for a photostim trial is considered a success.
        :param stims_to_use: ls of stims to retrieve photostim trial dFF responses
        :return:
        """
        if stims_to_use is None:
            stims_to_use = range(len(self.stim_start_frames))
            stims_idx = [self.stim_start_frames.index(stim) for stim in stims_to_use]
        elif stims_to_use:
            stims_idx = [self.stim_start_frames.index(stim) for stim in stims_to_use]
        else:
            KeyError('no stims set to analyse [1]')

        # choose between .SLMTargets_stims_dff and .SLMTargets_stims_tracedFF for data to process
        if process == 'dF/prestimF':
            if hasattr(self, 'SLMTargets_stims_dff'):
                targets_traces = self.SLMTargets_stims_dff
            else:
                AssertionError('no SLMTargets_stims_dff attr. [2]')
        elif process == 'trace dFF':
            if hasattr(self, 'SLMTargets_stims_dff'):
                targets_traces = self.SLMTargets_tracedFF_stims_dff
            else:
                AssertionError('no SLMTargets_tracedFF_stims_dff attr. [2]')
        else:
            ValueError('need to assign to process: dF/prestimF or trace dFF')

        # initializing pandas df that collects responses of stimulations
        if hasattr(self, 'SLMTargets_stims_dff'):
            d = {}
            for stim in stims_idx:
                d[stim] = [None] * targets_traces.shape[0]
            df = pd.DataFrame(d, index=range(targets_traces.shape[0]))  # population dataframe
        else:
            AssertionError('no SLMTargets_stims_dff attr. [2]')

        # initializing pandas df for binary showing of success and fails (1= success, 0= fails)
        hits_slmtargets = {}  # to be converted in pandas df below - will contain 1 for every success stim, 0 for non success stims
        for stim in stims_idx:
            hits_slmtargets[stim] = [None] * targets_traces.shape[0]  # start with 0 for all stims
        hits_slmtargets_df = pd.DataFrame(hits_slmtargets,
                                          index=range(targets_traces.shape[0]))  # population dataframe

        reliability_slmtargets = {}  # dict will be used to store the reliability results for each targeted cell

        # dFF response traces for successful photostim trials
        traces_dff_successes = {}
        cell_ids = df.index
        for target_idx in range(len(cell_ids)):
            traces_dff_successes_l = []
            success = 0
            counter = 0
            responses = []
            for stim_idx in stims_idx:
                dff_trace = targets_traces[target_idx][stim_idx]
                response_result = np.mean(dff_trace[self.pre_stim + self.stim_duration_frames + 1:
                                                    self.pre_stim + self.stim_duration_frames +
                                                    self.post_stim_response_frames_window])  # calculate the dF over pre-stim mean F response within the response window
                responses.append(round(response_result, 2))
                if response_result >= threshold:
                    success += 1
                    hits_slmtargets_df.loc[target_idx, stim_idx] = 1
                    traces_dff_successes_l.append(dff_trace)
                else:
                    hits_slmtargets_df.loc[target_idx, stim_idx] = 0

                df.loc[target_idx, stim_idx] = response_result
                counter += 1
            reliability_slmtargets[target_idx] = round(success / counter * 100., 2)
            traces_dff_successes[target_idx] = np.array(traces_dff_successes_l)

        return reliability_slmtargets, hits_slmtargets_df, df, traces_dff_successes

    # retrieves photostim avg traces for each SLM target, also calculates the reliability % for each SLM target
    def calculate_SLMTarget_SuccessStims(self, hits_df, process: str, stims_idx_l: list,
                                         exclude_stims_targets: dict = {}):
        """uses outputs of calculate_SLMTarget_responses_dff to calculate overall successrate of the specified stims

        :param hits_df: pandas dataframe of targets x stims where 1 denotes successful stim response (0 is failure)
        :param stims_idx_l: ls of stims to use for this function (useful when needing to filter out certain stims for in/out of sz)
        :param exclude_stims_targets: dictionary of stims (keys) where the values for each stim contains the targets that should be excluded from counting in the analysis of Success/failure of trial

        :return
            reliability_slmtargets: dict; reliability (% of successful stims) for each SLM target
            traces_SLMtargets_successes_avg: np.array; photostims avg traces for each SLM target (successful stims only)

        """

        # choose between .SLMTargets_stims_dff and .SLMTargets_stims_tracedFF for data to process
        if process == 'dF/prestimF':
            if hasattr(self, 'SLMTargets_stims_dff'):
                targets_traces = self.SLMTargets_stims_dff
            else:
                AssertionError('no SLMTargets_stims_dff attr. [2]')
        elif process == 'trace dFF':
            if hasattr(self, 'SLMTargets_stims_dff'):
                targets_traces = self.SLMTargets_tracedFF_stims_dff
            else:
                AssertionError('no SLMTargets_tracedFF_stims_dff attr. [2]')
        else:
            ValueError('need to assign to process: dF/prestimF or trace dFF')

        traces_SLMtargets_successes_avg_dict = {}
        traces_SLMtargets_failures_avg_dict = {}
        reliability_slmtargets = {}
        for target_idx in hits_df.index:
            traces_SLMtargets_successes_l = []
            traces_SLMtargets_failures_l = []
            success = 0
            counter = 0
            for stim_idx in stims_idx_l:
                if stim_idx in exclude_stims_targets.keys():
                    if target_idx not in exclude_stims_targets[stim_idx]:
                        continu_ = True
                    else:
                        continu_ = False
                else:
                    continu_ = True
                if continu_:
                    counter += 1
                    if hits_df.loc[target_idx, stim_idx] == 1:
                        success += 1
                        dff_trace = targets_traces[target_idx][stim_idx]
                        traces_SLMtargets_successes_l.append(dff_trace)
                    else:
                        success += 0
                        dff_trace = targets_traces[target_idx][stim_idx]
                        traces_SLMtargets_failures_l.append(dff_trace)

            if counter > 0:
                reliability_slmtargets[target_idx] = round(success / counter * 100., 2)
            if success > 0:
                traces_SLMtargets_successes_avg_dict[target_idx] = np.mean(traces_SLMtargets_successes_l, axis=0)
            if success < counter:  # this helps protect against cases where a trial is 100% successful (and there's no failures).
                traces_SLMtargets_failures_avg_dict[target_idx] = np.mean(traces_SLMtargets_failures_l, axis=0)

        return reliability_slmtargets, traces_SLMtargets_successes_avg_dict, traces_SLMtargets_failures_avg_dict

    ### ALLOPTICAL ANALYSIS - FOR ALL CELLS FROM SUITE2P  # good progress on this, almost done reviewing
    #### TEMP - need to ask about these two functions from

    ### TODO ROB: how important are these two functions? I also see that detrending is commented out in makeFluTrials - should we include or not?

    def _baselineFluTrial(self, flu_trial, stim_end):
        '''
        Subtract baseline from dff trials to normalise across cells

        Inputs:
            flu_trial           - [cell x frame] dff trial for all cells
        Outputs:
            baselined_flu_trial - detrended dff trial with zeros replacing stim artifact
        '''
        # baseline the flu_trial using pre-stim period mean flu for each cell
        baseline_flu = np.mean(flu_trial[:, :self.pre_stim_frames], axis=1)
        # repeat the baseline_flu value across all frames for each cell
        baseline_flu_stack = np.repeat(baseline_flu, flu_trial.shape[1]).reshape(flu_trial.shape)
        # subtract baseline values for each cell
        baselined_flu_trial = flu_trial - baseline_flu_stack

        # set stim artifact period to 0
        baselined_flu_trial[:, self.pre_stim_frames:stim_end] = 0

        return baselined_flu_trial

    def _detrendFluTrial(self, flu_trial, stim_end):
        '''
        Detrend dff trials to account for drift of signal over a trial

        Inputs:
            flu_trial           - [cell x frame] dff trial for all cells
            stim_end            - frame n of the stim end
        Outputs:
            detrended_flu_trial - detrended dff trial with zeros replacing stim artifact
        '''
        # set stim artifact period to 0
        flu_trial[:, self.pre_frames:stim_end] = 0

        # detrend and baseline-subtract the flu trial for all cells
        detrended_flu_trial = signal.detrend(self.Suite2p.raw, axis=1)
        baselined_flu_trial = self._baselineFluTrial(detrended_flu_trial)

        return baselined_flu_trial

    #### TEMP // end

    def _makePhotostimTrialFluSnippets(self, plane_flu: np.ndarray, plane: int = 0, stim_frames: list = None) -> np.ndarray:  # base code copied from Vape's _makeFluTrials
        '''
        Make Flu snippets timed on photostimulation, for each cell, for each stim instance. [cells x Flu frames x stims]  # TODO triple check order of this array's dimensions

        Inputs:
            plane_flu   - array of dff traces for all cells for this plane only
            plane       - imaging plane corresponding to plane_flu, default = 0 (for one plane datasets)
            stim_frames - optional, if provided then only use these photostim frames to collect photostim_array
        Outputs:
            photostim_array     - dFF peri-photostim Flu array [cell x Flu frames x trial]
        '''

        print('\n\-Collecting peri-stim traces ...')

        trial_array = []
        _stims = self.stim_start_frames[plane] if stim_frames is None else stim_frames

        assert plane_flu.ndim == 2, 'plane_flu needs to be of ndim: 2'
        assert _stims == self.stim_start_frames[plane], "stims not found in the stim frames list of this plane"

        for i, stim in enumerate(_stims):
            # get frame indices of entire trial from pre-stim start to post-stim end
            trial_frames = np.s_[stim - self.pre_stim_frames: stim + self.post_stim_frames]

            # use trial frames to extract this trial for every cell
            flu_trial = plane_flu[:, trial_frames]
            flu_trial_len = self.pre_stim_frames + self.post_stim_frames
            stim_end = self.pre_stim_frames + self.stim_duration_frames

            # catch timeseries which ended in the middle of an ongoing photostim instance
            if flu_trial.shape[1] == flu_trial_len:
                flu_trial = self._baselineFluTrial(flu_trial, stim_end)
                # only append trials of the correct length - will catch corrupt/incomplete data and not include
                if len(trial_array) == 0:
                    trial_array = flu_trial
                else:
                    trial_array = np.dstack((trial_array, flu_trial))
            else:
                print('**incomplete trial detected and not appended to trial_array**', end='\r')

        print(f'\nFinished collecting peri-stim traces, out shape: {trial_array.shape}')

        return trial_array

    def collectPhotostimResponses(self):

        # create parameters, slices, and subsets for making pre-stim and post-stim arrays to use in stats comparison
        # test_period = self.pre_stim_response_window / 1000  # sec
        # self.test_frames = int(self.ImagingParams.fps * test_period)  # test period for stats

        # mean pre and post stimulus (within post-stim response window) flu trace values for all cells, all trials
        self.__analysis_array = self.photostimFluArray
        self.__pre_array = np.mean(self.__analysis_array[:, self.pre_stim_test_slice, :],
                              axis=1)  # [cells x prestim frames] (avg'd taken over all stims)
        self.__post_array = np.mean(self.__analysis_array[:, self.post_stim_test_slice, :],
                               axis=1)  # [cells x poststim frames] (avg'd taken over all stims)

        # Vape's version for collection photostim response amplitudes
        # calculate amplitude of response for all cells, all trials
        all_amplitudes = self.__post_array - self.__pre_array

        df = pd.DataFrame(index=range(self.Suite2p.n_units), columns=self.stim_start_frames, data=all_amplitudes)

        return df

    def _allCellsPhotostimResponsesAnndata(self):   # NOT TESTED!
        """
        Creates annotated data (see anndata library) object based around the Ca2+ matrix of the imaging trial.

        """

        if self.photostimResponseAmplitudes:
            # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
            # build dataframe for obs_meta from suite2p stat information
            obs_meta = pd.DataFrame(
                columns=['original_index', 'photostim_target', 'photostim_exclusion_zone', 'prob_response', 'sig_responder'], index=range(len(self.Suite2p.n_units)))
            for idx in obs_meta.index:
                obs_meta.loc[idx, 'original_index'] = self.Suite2p.stat[idx]['original_index']
            obs_meta.loc[:, 'photostim_target'] = self.targeted_cells
            obs_meta.loc[:, 'photostim_exclusion_zone'] = self.exclude_cells
            obs_meta.loc[:, 'prob_response'] = self.prob_response
            obs_meta.loc[:, 'sig_responder'] = self.sig_units


            # SETUP THE VARIABLES ANNOTATIONS TO USE IN anndata
            # build dataframe for var annot's from paq file
            var_meta = pd.DataFrame(index=self.paq_channels, columns=range(self.stim_start_frames[PLANE]))
            for fr_idx in range(self.stim_start_frames[PLANE]):
                for index in [*self.sparse_paq_data]:
                    var_meta.loc[index, fr_idx] = self.sparse_paq_data[index][fr_idx]

            # BUILD LAYERS TO ADD TO anndata OBJECT
            layers = {'singleTrialSignificance': np.empty_like(self.photostimResponseAmplitudes) # PLACEHOLDER NOT IMPLEMENTED YET
                      }

            print(f"\n\----- CREATING annotated data object for photostim responses using AnnData:")
            adata = Utils.create_anndata2(X=self.photostimResponseAmplitudes, obs_meta=obs_meta, var_meta=var_meta.T, layers=layers)

            print(f"\t{adata}")
            return adata

        else:
            Warning(
                'could not create anndata. anndata creation only available if .photostimResponseAmplitudes have been collected (run .photostimProcessingAllCells())')


    def photostimProcessingAllCells(self, plane: int = 0):  # NOTE: not setup for multi-plane imaging processing yet...
        '''
        Take dfof trace for entire timeseries and break it up in to individual trials, calculate
        the mean amplitudes of response and statistical significance across all trials

        Inputs:
            plane             - imaging plane n
        '''
        print('\n----------------------------------------------------------------')
        print('running trial Processing for all cells ')
        print('----------------------------------------------------------------')

        # make trial arrays from dff data shape: [cells x stims x frames]

        self.photostimFluArray = self._makePhotostimTrialFluSnippets(plane_flu=normalize_dff(self.Suite2p.raw))
        self.photostimResponseAmplitudes = self.collectPhotostimResponses()

        ## create new anndata object for storing measured photostim responses from data, with other relevant data
        self.photostim_responses = self._allCellsPhotostimResponsesAnndata()

    def statisticalProcessingAllCells(self):
        """Runs statistical processing on photostim response arrays"""
        self.wilcoxons = self._runWilcoxonsTest(array1=self.__pre_array, array2=self.__post_array)
        self.sig_units = self._sigTestAvgResponse(p_vals=self.wilcoxons, alpha=0.1)

    @property
    def pre_stim_test_slice(self):
        "num of prestim frames used for quantification of photostim responses"
        return np.s_[self.pre_stim_frames - self.pre_stim_response_frames_window: self.pre_stim_frames]

    @property
    def post_stim_test_slice(self):
        "num of poststim frames used for quantification of photostim responses"
        stim_end = self.pre_stim_frames + self.stim_duration_frames
        return np.s_[stim_end: stim_end + self.post_stim_response_frames_window]


    #### STATISTICS FOR PHOTOSTIM RESPONSES
    def _runWilcoxonsTest(expobj, array1, array2):  # NOTE: not setup for multiplane cells yet

        # check if the two distributions of flu values (pre/post) are different
        assert array1.shape == array2.shape, 'shapes for .__pre_array and .__post_array need to be the same for wilcoxon test'
        wilcoxons = np.empty(len(array1))  # [cell (p-value)]

        for cell in range(len(array1)):
            wilcoxons[cell] = stats.wilcoxon(array2[cell], array1[cell])[1]

        return wilcoxons

    def _sigTestAvgResponse(expobj, p_vals: list, alpha=0.1):  # NOTE: not setup for multiplane cells yet
        """
        Uses the p values and a threshold for the Benjamini-Hochberg correction to return which
        cells are still significant after correcting for multiple significance testing
        """
        print('\n----------------------------------------------------------------')
        print('running statistical significance testing for nontargets response arrays ')
        print('----------------------------------------------------------------')

        sig_units = np.full_like(p_vals, False, dtype=bool)

        try:
            sig_units, _, _, _ = sm.stats.multitest.multipletests(p_vals, alpha=alpha, method='fdr_bh',
                                                                  is_sorted=False, returnsorted=False)
        except ZeroDivisionError:
            print('no cells responding')

        # p values without bonferroni correction
        no_bonf_corr = [i for i, p in enumerate(p_vals) if p < 0.05]
        expobj.nomulti_sig_units = np.zeros(len(expobj.s2p_nontargets), dtype='bool')
        expobj.nomulti_sig_units[no_bonf_corr] = True

        ## TODO - validate by Rob - commented out in his code, is this repeating the sigunits defined by multipletests call just above?
        # p values after bonferroni correction
        bonf_corr = [i for i,p in enumerate(p_vals) if p < 0.05 / expobj.Suite2p.n_units]
        sig_units = np.zeros(expobj.Suite2p.n_units, dtype='bool')
        sig_units[bonf_corr] = True

        return sig_units


    ## NOT REVIEWED FOR USAGE YET
    def _probResponse(self, plane,
                      trial_sig_calc):  ## FROM VAPE'S CODE, TODO NEED TO CHOOSE BETWEEN HERE AND BOTTOM RELIABILITY CODE
        '''
        Calculate the response probability, i.e. proportion of trials that each cell responded on

        Inputs:
            plane          - imaging plane n
            trial_sig_calc - indicating which calculation was used for significance testing ('dff'/'dfsf')
        '''
        n_trials = self.n_trials

        # get the number of responses for each across all trials
        if trial_sig_calc == 'dff':
            num_respond = np.array(self.trial_sig_dff[plane])  # trial_sig_dff is [plane][cell][trial]
        elif trial_sig_calc == 'dfsf':
            num_respond = np.array(self.trial_sig_dfsf[plane])

        # calculate the proportion of all trials that each cell responded on
        self.prob_response.append(np.sum(num_respond, axis=1) / n_trials)

    def cellStaProcessing(self, test='t_test'):

        if self.stim_start_frames:

            # this is the key parameter for the sta, how many frames before and after the stim onset do you want to use
            self.pre_frames = int(np.ceil(self.ImagingParams.fps * 0.5))  # 500 ms pre-stim period
            self.post_frames = int(np.ceil(self.ImagingParams.fps * 3))  # 3000 ms post-stim period

            # ls of cell pixel intensity values during each stim on each trial
            self.all_trials = []  # ls 1 = cells, ls 2 = trials, ls 3 = dff vector

            # the average of every trial
            self.stas = []  # ls 1 = cells, ls 2 = sta vector

            self.all_amplitudes = []
            self.sta_amplitudes = []

            self.t_tests = []
            self.wilcoxons = []

            for plane in range(self.ImagingParams.n_planes):

                all_trials = []  # ls 1 = cells, ls 2 = trials, ls 3 = dff vector

                stas = []  # ls 1 = cells, ls 2 = sta vector

                all_amplitudes = []
                sta_amplitudes = []

                t_tests = []
                wilcoxons = []

                # loop through each cell
                for i, unit in enumerate(self.raw[plane]):

                    trials = []
                    amplitudes = []
                    df = []

                    # a flat ls of all observations before stim occured
                    pre_obs = []
                    # a flat ls of all observations after stim occured
                    post_obs = []

                    for stim in self.stim_start_frames[plane]:
                        # get baseline values from pre_stim_sec
                        pre_stim_f = unit[stim - self.pre_frames: stim]
                        baseline = np.mean(pre_stim_f)

                        # the whole trial and dfof using baseline
                        trial = unit[stim - self.pre_frames: stim + self.post_frames]
                        trial = [((f - baseline) / baseline) * 100 for f in trial]  # dff calc
                        trials.append(trial)

                        # calc amplitude of response
                        pre_f = trial[: self.pre_frames - 1]
                        pre_f = np.mean(pre_f)

                        avg_post_start = self.pre_frames + (self.stim_duration_frames + 1)
                        avg_post_end = avg_post_start + int(np.ceil(self.ImagingParams.fps * 0.5))  # post-stim period of 500 ms

                        post_f = trial[avg_post_start: avg_post_end]
                        post_f = np.mean(post_f)
                        amplitude = post_f - pre_f
                        amplitudes.append(amplitude)

                        # append to flat lists
                        pre_obs.append(pre_f)
                        post_obs.append(post_f)

                    trials = np.array(trials)
                    all_trials.append(trials)

                    # average amplitudes across trials
                    amplitudes = np.array(amplitudes)
                    all_amplitudes.append(amplitudes)
                    sta_amplitude = np.mean(amplitudes, 0)
                    sta_amplitudes.append(sta_amplitude)

                    # average across all trials
                    sta = np.mean(trials, 0)
                    stas.append(sta)

                    # remove nans from flat lists
                    pre_obs = [x for x in pre_obs if ~np.isnan(x)]
                    post_obs = [x for x in post_obs if ~np.isnan(x)]

                    # t_test and man whit test pre and post stim (any other test could also be used here)
                    t_test = stats.ttest_rel(pre_obs, post_obs)
                    t_tests.append(t_test)

                    wilcoxon = stats.wilcoxon(pre_obs, post_obs)
                    wilcoxons.append(wilcoxon)

                self.all_trials.append(np.array(all_trials))
                self.stas.append(np.array(stas))

                self.all_amplitudes.append(np.array(all_amplitudes))
                self.sta_amplitudes.append(np.array(sta_amplitudes))

                self.t_tests.append(np.array(t_tests))
                self.wilcoxons.append(np.array(wilcoxons))

            plt.figure()
            plt.plot([avg_post_start] * 2, [-1000, 1000])
            plt.plot([avg_post_end] * 2, [-1000, 1000])
            plt.plot([self.pre_frames - 1] * 2, [-1000, 1000])
            plt.plot([0] * 2, [-1000, 1000])
            plt.plot(stas[5])
            plt.plot(stas[10])
            plt.plot(stas[15])
            plt.ylim([-100, 200])

            self.staSignificance(test)
            self.singleTrialSignificance()

    def staSignificance(self, test):

        self.sta_sig = []

        for plane in range(self.ImagingParams.n_planes):

            # set this to true if you want to multiple comparisons correct for the number of cells
            multi_comp_correction = True
            if not multi_comp_correction:
                divisor = 1
            else:
                divisor = self.n_units[plane]

            if test == 't_test':
                p_vals = [t[1] for t in self.t_tests[plane]]
            if test == 'wilcoxon':
                p_vals = [t[1] for t in self.wilcoxons[plane]]

            if multi_comp_correction:
                print('performing t-test on cells with mutliple comparisons correction')
            else:
                print('performing t-test on cells without mutliple comparisons correction')

            sig_units = []

            for i, p in enumerate(p_vals):
                if p < (0.05 / divisor):
                    unit_index = self.cell_id[plane][i]
                    # print('stimulation has significantly changed fluoresence of s2p unit {}, its P value is {}'.format(unit_index, p))
                    sig_units.append(unit_index)  # significant units

            self.sta_sig.append(sig_units)

    def singleTrialSignificance(self):

        self.single_sig = []  # single trial significance value for each trial for each cell in each plane

        for plane in range(self.ImagingParams.n_planes):

            single_sigs = []

            for cell, _ in enumerate(self.cell_id[plane]):

                single_sig = []

                for trial in range(self.n_trials):

                    pre_f_trial = self.all_trials[plane][cell][trial][: self.pre_frames]
                    std = np.std(pre_f_trial)

                    if np.absolute(self.all_amplitudes[plane][cell][trial]) >= 2 * std:
                        single_sig.append(True)
                    else:
                        single_sig.append(False)

                single_sigs.append(single_sig)

            self.single_sig.append(single_sigs)

    ## NOT REVIEWED FOR USAGE YET







    # other useful functions for all-optical analysis
    def whiten_photostim_frame(self, tiff_path, save_as=''):
        im_stack = tf.imread(tiff_path, key=range(self.ImagingParams.n_frames))

        frames_to_whiten = []
        for j in self.stim_start_frames:
            frames_to_whiten.append(j)

        im_stack_1 = im_stack
        a = np.full_like(im_stack_1[0], fill_value=0)
        a[0:100, 0:100] = 5000.
        for frame in frames_to_whiten:
            im_stack_1[frame - 3] = im_stack_1[frame - 3] + a
            im_stack_1[frame - 2] = im_stack_1[frame - 2] + a
            im_stack_1[frame - 1] = im_stack_1[frame - 1] + a

        frames_to_remove = []
        for j in self.stim_start_frames:
            for i in range(0,
                           self.stim_duration_frames + 1):  # usually need to remove 1 more frame than the stim duration, as the stim isn't perfectly aligned with the start of the imaging frame
                frames_to_remove.append(j + i)

        im_stack_1 = np.delete(im_stack, frames_to_remove, axis=0)

        tf.imwrite(save_as, im_stack_1, photometric='minisblack')

    def avg_stim_images(self, peri_frames: int = 100, stim_timings: list = [], save_img=False, to_plot=False,
                        verbose=False, force_redo=False):
        """
        Outputs (either by saving or plotting, or both) images from raw t-series TIFF for a trial around each individual
        stim timings.

        :param peri_frames:
        :param stim_timings:
        :param save_img:
        :param to_plot:
        :param force_redo:
        :param verbose:
        :return:
        """

        if force_redo:
            continu = True
        elif hasattr(self, 'avgstimimages_r'):
            if self.avgstimimages_r is True:
                continu = False
            else:
                continu = True
        else:
            continu = True

        if continu:
            print('making stim images...')
            if hasattr(self, 'stim_images'):
                x = [0 for stim in stim_timings if stim not in self.stim_images.keys()]
            else:
                self.stim_images = {}
                x = [0] * len(stim_timings)
            if 0 in x:
                tiffs_loc = '%s/*Ch3.tif' % self.tiff_path_dir
                tiff_path = glob.glob(tiffs_loc)[0]
                print('working on loading up %s tiff from: ' % self.metainfo['trial'], tiff_path)
                im_stack = tf.imread(tiff_path, key=range(self.ImagingParams.n_frames))
                print('Processing seizures from experiment tiff (wait for all seizure comparisons to be processed), \n '
                      'total tiff shape: ', im_stack.shape)

            for stim in stim_timings:
                message = '|- stim # %s out of %s' % (stim_timings.index(stim), len(stim_timings))
                print(message, end='\r')
                if stim in self.stim_images.keys():
                    avg_sub = self.stim_images[stim]
                else:
                    if stim < peri_frames:
                        peri_frames = stim
                    im_sub = im_stack[stim - peri_frames: stim + peri_frames]
                    avg_sub = np.mean(im_sub, axis=0)
                    self.stim_images[stim] = avg_sub

                if save_img:
                    # save in a subdirectory under the ANALYSIS folder path from whence t-series TIFF came from
                    save_path = self.analysis_save_dir + 'avg_stim_images'
                    save_path_stim = save_path + '/%s_%s_stim-%s.tif' % (
                        self.metainfo['date'], self.metainfo['trial'], stim)
                    if os.path.exists(save_path):
                        print("saving stim_img tiff to... %s" % save_path_stim) if verbose else None
                        avg_sub8 = convert_to_8bit(avg_sub, 0, 255)
                        tf.imwrite(save_path_stim,
                                   avg_sub8, photometric='minisblack')
                    else:
                        print('making new directory for saving images at:', save_path)
                        os.mkdir(save_path)
                        print("saving as... %s" % save_path_stim)
                        avg_sub8 = convert_to_8bit(avg_sub, 0, 255)
                        tf.imwrite(save_path_stim,
                                   avg_sub, photometric='minisblack')

                if to_plot:
                    plt.imshow(avg_sub, cmap='gray')
                    plt.suptitle('avg image from %s frames around stim_start_frame %s' % (peri_frames, stim))
                    plt.show()  # just plot for now to make sure that you are doing things correctly so far

            if hasattr(self, 'pkl_path'):
                self.save_pkl()
            else:
                print('note: pkl not saved yet...')

            self.avgstimimages_r = True

        else:
            print('skipping remaking of avg stim images')

    def run_stamm_nogui(self, numDiffStims, startOnStim, everyXStims, preSeconds=0.75, postSeconds=1.25):
        """run STAmoviemaker for the expobj's trial"""
        qnap_path = os.path.expanduser('/home/pshah/mnt/qnap')

        ## data path
        movie_path = self.tiff_path
        sync_path = self.paq_path

        ## stamm save path
        stam_save_path = os.path.join(qnap_path, 'Analysis', self.metainfo['date'], 'STA_Movies',
                                      '%s_%s_%s' % (self.metainfo['date'],
                                                    self.metainfo['animal prep.'],
                                                    self.metainfo['trial']))
        os.makedirs(stam_save_path, exist_ok=True)

        ##
        assert os.path.exists(stam_save_path)

        print('QNAP_path:', qnap_path,
              '\ndata path:', movie_path,
              '\nsync path:', sync_path,
              '\nSTA movie save path:', stam_save_path)

        # define STAmm parameters
        frameRate = int(self.ImagingParams.fps)

        arg_dict = {'moviePath': movie_path,  # hard-code this
                    'savePath': stam_save_path,
                    'syncFrameChannel': 'frame_clock',
                    'syncStimChannel': 'packio2markpoints',
                    'syncStartSec': 0,
                    'syncStopSec': 0,
                    'numDiffStims': numDiffStims,
                    'startOnStim': startOnStim,
                    'everyXStims': everyXStims,
                    'preSeconds': preSeconds,
                    'postSeconds': postSeconds,
                    'frameRate': frameRate,
                    'averageImageStart': 0.5,
                    'averageImageStop': 1.5,
                    'methodDF': False,
                    'methodDFF': True,
                    'methodZscore': False,
                    'syncPath': sync_path,
                    'zPlanes': 1,
                    'useStimOrder': False,
                    'stimOrder': [],
                    'useSingleTrials': False,
                    'doThreshold': False,
                    'threshold': 0,
                    'colourByTime': False,
                    'useCorrelationImage': False,
                    'blurHandS': False,
                    'makeMaxImage': True,
                    'makeColourImage': False
                    }

        # # run STAmm
        # STAMM.STAMovieMaker(arg_dict)

        # show the MaxResponseImage
        img = glob.glob(stam_save_path + '/*MaxResponseImage.tif')[0]
        plotting.plot_single_tiff(img, frame_num=0)



