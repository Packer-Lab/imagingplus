# TODO update all 2p stim related attr's to naparm submodule
from dataclasses import dataclass
import glob
import os
import signal
import time
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import tifffile as tf

from packerlabimaging import TwoPhotonImaging
from packerlabimaging.main.classes import ImagingMetadata, ImagingData, TemporalData, ImagingTrial, CellAnnotations, \
    Experiment
from packerlabimaging.main.paq import PaqData
from packerlabimaging.utils.io import import_obj
from packerlabimaging.utils.utils import convert_to_8bit
from packerlabimaging.processing.naparm import Targets
from packerlabimaging.utils.classes import UnavailableOptionError
# %%
from packerlabimaging.processing.anndata import AnnotatedData

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)

PLANE = 0
BADFRAMESLOC = '/home/pshah/Documents/code/packerlabimaging/tests/'

prestim_sec: float = 1.0  #: length of pre stim trace collected (secs)
poststim_sec: float = 3.0  #: length of post stim trace collected (secs)
pre_stim_response_window: float = 0.500  #: time window for collecting pre-stim measurement (units: msec)
post_stim_response_window: float = 0.500  #: time window for collecting post-stim measurement (units: msec)


class AllOpticalTrial(TwoPhotonImaging):
    """All Optical Experimental Data Analysis Workflow."""

    def __init__(self, naparm_path, dataPath: str, saveDir: str, date: str, trialID: str, expID: str,
                 expGroup: str = '',
                 comment: str = '', imparams: ImagingMetadata = None, cells: CellAnnotations = None,
                 imdata: ImagingData = None,
                 tmdata: PaqData = None):

        """
        :param metainfo: TypedDict containing meta-information field needed for this experiment. Please see TwoPhotonImagingMetainfo for type hints on accepted keys.
        :param paq_options: TypedDict containing meta-information about .paq file associated with current trial
        :param naparm_path: path to folder containing photostimulation setup built by NAPARM
        :param analysis_save_path: path of where to save the experiment analysis object
        :param microscope: name of microscope used to record imaging (options: "Bruker" (default), "other")
        :param prestim_sec: pre-photostimulation timeperiod for collecting photostimulation timed signals
        :param poststim_sec: post-photostimulation timeperiod for collecting photostimulation timed signals
        :param pre_stim_response_window: pre-photostimulation time window for measurement of photostimulation evoked responses
        :param post_stim_response_window: post-photostimulation time window for measurement of photostimulation evoked responses

        """

        initialization_dict = {'date': date,
                               'trialID': trialID,
                               'dataPath': dataPath,
                               'saveDir': saveDir,
                               'expID': expID,
                               'expGroup': expGroup,
                               'comment': comment}

        # 1) initialize object as a TwoPhotonImagingTrial
        super().__init__(imparams=imparams, cells=cells, tmdata=tmdata, **initialization_dict)

        # Initialize Attributes:

        # PHOTOSTIM PROTOCOL
        self.stim_start_times = None
        self.nomulti_sig_units = None

        # attr's for processing/analysis of photostim experiments
        self.__prestim_sec: Union[float, int] = 1.0  #: length of pre stim trace collected (in secs)
        self.__poststim_sec: Union[float, int] = 3.0  #: length of post stim trace collected (in secs)
        self.__pre_stim_response_window: Union[
            float, int] = 0.500  #: time window for collecting pre-stim measurement (units: msec)
        self.__post_stim_response_window: Union[
            float, int] = 0.500  #: time window for collecting post-stim measurement (units: msec)

        # attr's for statistical analysis of photostim trials responses
        self.photostimResponsesData = None  # anndata object for collecting photostim responses and associated metadata for experiment and cells

        # TODO update comment descriptions
        self.all_trials = []  # all trials for each cell, dff detrended
        self.all_amplitudes = []  # all amplitudes of response between dff test periods

        self.stas = []  # avg of all trials for each cell, dff
        self.sta_amplitudes = []  # avg amplitude of response between dff test periods

        self.prob_response = None  # probability of response of cell to photostim trial; obtained from single trial significance (ROB suggestion)
        self.t_tests = []  # result from related samples t-test between dff test periods
        self.wilcoxons = []  # ROB to update
        self.sig_units = None  # ROB to update
        self.trial_sig_dff = []  # based on dff increase above std of baseline
        self.trial_sig_dfsf = []  # based on df/std(f) increase in test period post-stim
        self.sta_sig = []  # based on t-test between dff test periods
        self.sta_sig_nomulti = []  # as above, no multiple comparisons correction
        ########

        # initializing data processing, data analysis and/or results associated attr's

        # PHOTOSTIM SLM TARGETS
        # TODO add attr's related to numpy array's and pandas dataframes for photostim trials - SLM targets
        self.responses_SLMtargets = []  # dF/prestimF responses for all SLM targets for each photostim trial
        self.responses_SLMtargets_tracedFF = []  # poststim dFF - prestim dFF responses for all SLM targets for each photostim trial

        # ALL CELLS (from suite2p ROIs)
        # TODO add attr's related to numpy array's and pandas dataframes for photostim trials - suite2p ROIs

        # FUNCTIONS TO RUN AFTER init's of ALL ATTR'S

        # 3) process 2p stim protocol and collect imaging frames at stim starts and during photostimulation
        # set naparm path
        self.twopstim, self.twopstim.stim_start_frames, self.twopstim.photostim_frames = self.photostimProcessing(
            naparm_path=naparm_path)

        # 5) todo collect Flu traces from SLM targets - probably leave out of the init right??
        if hasattr(self, 'Suite2p'):
            self.raw_SLMTargets, self.dFF_SLMTargets, self.meanFluImg_registered = self.collect_traces_from_targets(
                curr_trial_frames=self.Suite2p.trial_frames, save=True)
            self.targets_dff, self.targets_dff_avg, self.targets_dfstdF, self.targets_dfstdF_avg, self.targets_raw, self.targets_raw_avg = self.get_alltargets_stim_traces_norm(
                process='trace dFF')
        else:
            Warning('NO Flu traces collected from any SLM targets because Suite2p not found for trial.')

        # 6) Collect Flu traces from Suite2p ROIs:
        #   create:
        #       1) array of dFF pre + post stim Flu snippets for each stim and cell [num cells x num peri-stim frames collected x num stims]
        #       2) photostim response amplitudes in a dataframe for each cell and each photostim
        #       3) save photostim response amplitudes to AnnotatedData
        self.photostimFluArray, self.photostimResponseAmplitudes, self.photostimResponsesData = self.photostimProcessingAllCells()

        # extend annotated imaging data object with imaging frames in photostim and stim_start_frames as additional keys in vars
        __frames_in_stim = [False] * self.imparams.n_frames
        __stim_start_frame = [False] * self.imparams.n_frames
        for frame in self.twopstim.photostim_frames: __frames_in_stim[frame] = True
        for frame in self.twopstim.stim_start_frames: __stim_start_frame[frame] = True
        self.data.add_var(var_name='photostim_frame', values=__frames_in_stim)
        self.data.add_var(var_name='stim_start_frame', values=__stim_start_frame)

        # save final object
        self.save()
        print(f'\----- CREATED AllOpticalTrial data object for {self.t_series_name}')

    def __str__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
            information = self.t_series_name
            return f"({information}) TwoPhotonImagingTrial.alloptical experimental trial object, last saved: {lastmod}"
        else:
            return f" -- unsaved TwoPhotonImagingTrial.alloptical experimental trial object -- "

    def __repr__(self):
        return f"TwoPhotonImagingTrial.alloptical experimental trial object"

    @classmethod
    def AllOpticalTrialfromImagingTrial(cls, naparm_path, imaging_trial: ImagingTrial):
        """Alternative constructor for AllOpticalTrial.

        Creates an all optical trial from an existing imaging trial.
        """

        initialization_dict = {'naparm_path': naparm_path, 'dataPath': imaging_trial.dataPath,
                               'saveDir': imaging_trial.saveDir,
                               'expID': imaging_trial.expID, 'group': imaging_trial.expGroup,
                               'comment': imaging_trial.comment}

        aotrial = cls(**initialization_dict)

        return aotrial

    @property
    def twopstim_path(self):
        """path to folder containing photostimulation protocols output by NAPARM"""
        if self.twopstim:
            return self.twopstim.path

    @property
    def prestim_sec(self):
        """length of pre stim trace collected (secs)"""
        return self.__prestim_sec

    @prestim_sec.setter
    def prestim_sec(self, val):
        assert type(val) == int or type(val) == float, 'can only set prestim_sec with int or float'
        self.__prestim_sec = val

    @property
    def poststim_sec(self):
        """length of post stim trace collected (secs)"""
        return self.__poststim_sec

    @poststim_sec.setter
    def poststim_sec(self, val):
        assert type(val) == int or type(val) == float, 'can only set poststim_sec with int or float'
        self.__prestim_sec = val

    @property
    def pre_stim_response_window(self):
        """time window for collecting pre-stim measurement (units: msec)"""
        return self.__pre_stim_response_window

    @pre_stim_response_window.setter
    def pre_stim_response_window(self, val):
        assert type(val) == int or type(val) == float, 'can only set pre_stim_response_window with int or float'
        self.__pre_stim_response_window = val

    @property
    def pre_stim_response_frames_window(self):
        """num frames for measuring Flu trace before each photostimulation trial during photostim response measurement (units: frames)"""
        return int(self.imparams.fps * self.pre_stim_response_window)

    @property
    def post_stim_response_window(self):
        """time window for collecting post-stim measurement (units: msec)"""
        return self.__post_stim_response_window

    @post_stim_response_window.setter
    def post_stim_response_window(self, val):
        assert type(val) == int or type(val) == float, 'can only set post_stim_response_window with int or float'
        self.__post_stim_response_window = val

    @property
    def post_stim_response_frames_window(self):
        """num frames for measuring Flu response after each photostimulation trial during photostim response measurement (units: frames)"""
        return int(self.imparams.fps * self.post_stim_response_window)

    @property
    def pre_stim_frames(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)"""
        return int(self.__prestim_sec * self.imparams.fps)

    @property
    def post_stim_frames(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)"""
        return int(self.__poststim_sec * self.imparams.fps)

    @property
    def n_stims(self):
        #: TODO set property using anndata responses shape property assigning from array: cells x Flu frames x # of photostim trials
        """number of photostimulation trials """
        return

    @property
    def stim_start_frames(self):
        """Imaging frames corresponding to start of photostimulation trials."""
        # todo use anndata for getting this
        return []

    @property
    def photostim_frames(self):
        """Imaging frames during photostimulation trials."""
        # todo just use anndata for getting this
        return []

    @property
    def stim_duration_frames(self):
        """Duration of photostimulation as number of frames.
        Note: Removing 1 more frame than the stim duration, as the stim isn't perfectly aligned with the start of the imaging frame
        """
        if not hasattr(self, 'twopstim'):
            raise ValueError(
                'cannot retrieve stim_duration_frames. photostimulation analysis module cannot be found under .twopstim')
        else:
            duration_ms = self.twopstim.stim_dur
            frame_rate = self.imparams.fps
            duration_frames = np.ceil((duration_ms / 1000) * frame_rate)
            return int(duration_frames) + 1

    @property
    def pre_stim_test_slice(self):
        """slice representing prestim response frames"""
        return np.s_[self.pre_stim_frames - self.pre_stim_response_frames_window: self.pre_stim_frames]

    @property
    def post_stim_test_slice(self):
        """slice representing poststim response frames"""
        stim_end = self.pre_stim_frames + self.stim_duration_frames
        return np.s_[stim_end: stim_end + self.post_stim_response_frames_window]

    def photostimProcessing(self, naparm_path):
        """
        Processing alloptical trial photostimulation protocol.
        - parse NAPARM protocol and SLM Targets information under .twopstim
        - collect stimulation timing synchronized to imaging data frame timings

        """

        self.twopstim = Targets(naparm_path=naparm_path, frame_x=self.imparams.frame_x,
                                frame_y=self.imparams.frame_y,
                                pix_sz_x=self.imparams.pix_sz_x, pix_sz_y=self.imparams.pix_sz_y)

        # find stim frames
        stim_start_frames = []
        for stim in self.tmdata.data['stim_start_times']:
            # the index of the frame immediately preceeding stim
            stim_start_frame = next(i - 1 for i, sample in enumerate(self.tmdata.frame_times) if sample - stim >= 0)
            stim_start_frames.append(stim_start_frame)

        # find all photostimulation frames in imaging series
        print('\n\----- Finding photostimulation frames in imaging frames ...')
        print('# of photostim frames calculated per stim. trial: ', self.stim_duration_frames)

        photostim_frames = []

        for j in stim_start_frames:
            for i in range(
                    self.stim_duration_frames):
                photostim_frames.append(j + i)

        print('\t|- # of Imaging frames:', self.imparams.n_frames,
              'frames')  #: todo set this to n_frames property from trialdata
        print('\t|- # of Photostimulation frames:', len(self.photostim_frames), 'frames')

        return self.twopstim, stim_start_frames, photostim_frames

    def _findTargetedCells(self, plot: bool = True):
        """finding s2p cell ROIs that were also SLM targets (or more specifically within the target areas as specified by _findTargetAreas - include 15um radius from center coordinate of spiral)
        Make a binary mask of the targets and multiply by an image of the cells
        to find cells that were targeted

        --- LAST UPDATED NOV 6 2021 - copied over from Vape ---
        """

        print('\n\----- Searching for targeted cells in annotated cells...')

        ## TODO add necessary edits for multi-plane experiments

        targets = list(['non-target' for cell in
                        range(self.cells.n_cells)])  # initialize all cells as non-target, add as annotation to .cells

        ##### IDENTIFYING S2P ROIs THAT ARE WITHIN THE SLM TARGET SPIRAL AREAS
        # make all target area coords in to a binary mask
        targ_img = np.zeros([self.imparams.frame_x, self.imparams.frame_y], dtype='uint16')
        target_areas = np.array(self.twopstim.target_areas)
        targ_img[target_areas[:, :, 1], target_areas[:, :, 0]] = 1

        # make an image of every cell area, filled with the index of that cell
        cell_img = np.zeros_like(targ_img)

        cell_x = self.cells.cell_coords[:, 0]
        cell_y = self.cells.cell_coords[:, 1]

        for i, coord in enumerate(zip(cell_x, cell_y)):
            cell_img[coord] = i + 1

        # binary mask x cell image to get the cells that overlap with target areas
        targ_cell = cell_img * targ_img

        targ_cell_ids = np.unique(targ_cell)[1:] - 1  # correct the cell id due to zero indexing
        self.targeted_cells = np.zeros([self.cells.n_cells], dtype='bool')
        self.targeted_cells[targ_cell_ids] = True
        self.cell_targets = [self.cells.cell_id[i] for i in np.where(self.targeted_cells)[
            0]]  # get list of cells that were photostim targetted -- todo turn into property accessing targets annotations of .cells

        self.n_targeted_cells = np.sum(
            self.targeted_cells)  # todo turn into property accessing targets annotations of .cells

        # add targeted cells to targets
        for idx, cell in enumerate(self.cells.cell_id):
            if cell in self.cell_targets:
                targets[idx] = 'target'

        print('\t|- Search completed.')
        self.save()
        print('\t|- Number of targeted cells: ', self.n_targeted_cells)

        # IDENTIFYING S2P ROIs THAT ARE WITHIN THE EXCLUSION ZONES OF THE SLM TARGETS
        # make all target area coords in to a binary mask
        targ_img = np.zeros([self.imparams.frame_x, self.imparams.frame_y], dtype='uint16')
        target_areas_exclude = np.array(self.twopstim.target_areas_exclude)
        targ_img[target_areas_exclude[:, :, 1], target_areas_exclude[:, :, 0]] = 1

        # make an image of every cell area, filled with the index of that cell
        cell_img = np.zeros_like(targ_img)

        cell_x = self.cells.cell_coords[:, 0]
        cell_y = self.cells.cell_coords[:, 1]

        for i, coord in enumerate(zip(cell_x, cell_y)):
            cell_img[coord] = i + 1

        # binary mask x cell image to get the cells that overlap with target areas
        targ_cell = cell_img * targ_img

        targ_cell_ids = np.unique(targ_cell)[1:] - 1  # correct the cell id due to zero indexing
        exclude_cells = np.zeros([self.cells.n_cells], dtype='bool')
        exclude_cells[targ_cell_ids] = True
        cells_exclude = [self.cells.n_cells[i] for i in
                                  np.where(exclude_cells)[0]]  # get ls of s2p cells that were photostim targetted

        self.n_exclude_cells = np.sum(exclude_cells)

        print('\t|-Search completed.')
        self.save()
        print(f"\t|-Number of exclude Suite2p ROIs: {self.n_exclude_cells}")

        # define non targets from suite2p ROIs (exclude cells in the SLM targets exclusion - .s2p_cells_exclude)
        for idx, cell in enumerate(self.cells.cell_id):
            if cell not in cells_exclude:
                targets[idx] = 'exclude'

        if plot:
            fig, ax = plt.subplots(figsize=[6, 6])
            targ_img = np.zeros([self.imparams.frame_x, self.imparams.frame_y], dtype='uint16')
            target_areas = np.array(self.twopstim.target_areas)
            targ_img[target_areas[:, :, 1], target_areas[:, :, 0]] = 1
            ax.imshow(targ_img, cmap='Greys_r', zorder=0)
            ax.set_title('SLM targets areas')
            # for (x, y) in self.twopstim.target_coords_all:
            #     ax.scatter(x=x, y=y, edgecolors='white', facecolors='none', linewidths=1.0)
            fig.show()

        # add targets classification as observations annotation to .data anndata
        self.data.add_obs(obs_name='photostim_class', values=targets)

        self.cells.cellsdata['photostim_class'] = targets

        print(f"\t|- Number of non-target ROIs: {len(self.cells.cellsdata['photostim_class'] == 'non-target')}")

    #### TODO review attr's and write docs from the following functions: // start


    ### ALLOPTICAL PROCESSING OF TRACES
    ## ... no methods determined here yet...

    ### ALLOPTICAL ANALYSIS - FOCUS ON SLM TARGETS RELATED METHODS
    def collect_traces_from_targets(self, curr_trial_frames: list, reg_tif_folder: str = None, save: bool = True):
        """uses registered tiffs to collect raw traces from SLM target areas

        :param curr_trial_frames:
        :param reg_tif_folder:
        :param save:
        :return:
        """

        if reg_tif_folder is None:
            if self.Suite2p.s2pResultsPath:
                reg_tif_folder = self.Suite2p.s2pResultsPath + '/reg_tif/'
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

        mean_img_stack = np.zeros([end - start, self.imparams.frame_x, self.imparams.frame_y])
        # collect mean traces from target areas of each target coordinate by reading in individual registered tiffs that contain frames for current trial
        targets_trace_full = np.zeros([len(self.twopstim.target_coords_all), (end - start) * 2000], dtype='float32')
        counter = 0
        for i in range(start, end):
            tif_path_save2 = self.Suite2p.s2pResultsPath + '/reg_tif/' + reg_tif_list[i]
            with tf.TiffFile(tif_path_save2, multifile=False) as input_tif:
                print('|- reading tiff: %s' % tif_path_save2)
                data = input_tif.asarray()

            targets_trace = np.zeros([len(self.twopstim.target_coords_all), data.shape[0]], dtype='float32')
            for coord in range(len(self.twopstim.target_coords_all)):
                target_areas = np.array(
                    self.twopstim.target_areas)  # TODO update this so that it doesn't include the extra exclusion zone
                x = data[:, target_areas[coord, :, 1], target_areas[coord, :, 0]]  # = 1
                targets_trace[coord] = np.mean(x, axis=1)

            targets_trace_full[:, (i - start) * 2000: ((i - start) * 2000) + data.shape[
                0]] = targets_trace  # iteratively write to each successive segment of the targets_trace array based on the length of the reg_tiff that is read in.

            mean_img_stack[counter] = np.mean(data, axis=0)
            counter += 1

        # final part, crop to the *exact* frames for current trial
        raw_SLMTargets = targets_trace_full[:,
                         curr_trial_frames[0] - start * 2000: curr_trial_frames[1] - (start * 2000)]

        dFF_SLMTargets = self.normalize_dff(raw_SLMTargets, threshold_pct=10)

        meanFluImg_registered = np.mean(mean_img_stack, axis=0)

        self.save() if save else None

        return raw_SLMTargets, dFF_SLMTargets, meanFluImg_registered

    def get_alltargets_stim_traces_norm(self, process: str, targets_idx: int = None, subselect_cells: list = None,
                                        pre_stim=0.5, post_stim=4.0, stims: list = None):
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

        self.prestim_sec = pre_stim
        self.poststim_sec = post_stim

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
            flu = [targets_trace[targets_idx][stim - pre_stim: stim + self.stim_duration_frames + post_stim] for stim in
                   stim_timings]
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

    # calculate photostim. response magnitude of photostim responsiveness of all of the targeted cells TODO need to review this whole section
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
                raise AttributeError('no SLMTargets_stims_dff attr. [2]')
        elif process == 'trace dFF':
            if hasattr(self, 'SLMTargets_stims_dff'):
                targets_traces = self.SLMTargets_tracedFF_stims_dff
            else:
                raise AttributeError('no SLMTargets_tracedFF_stims_dff attr. [2]')
        else:
            ValueError('need to assign to process: dF/prestimF or trace dFF')

        # initializing pandas df that collects responses of stimulations
        if hasattr(self, 'SLMTargets_stims_dff'):
            d = {}
            for stim in stims_idx:
                d[stim] = [None] * targets_traces.shape[0]
            df = pd.DataFrame(d, index=range(targets_traces.shape[0]))  # population dataframe
        else:
            raise AttributeError('no SLMTargets_stims_dff attr. [2]')

        # dFF response traces for successful photostim trials
        cell_ids = df.index
        for target_idx in range(len(cell_ids)):
            responses = []
            for stim_idx in stims_idx:
                dff_trace = targets_traces[target_idx][stim_idx]
                response_result = np.mean(dff_trace[self.pre_stim + self.stim_duration_frames + 1:
                                                    self.pre_stim + self.stim_duration_frames +
                                                    self.post_stim_response_frames_window])  # calculate the dF over pre-stim mean F response within the response window
                responses.append(np.round(response_result, 2))

                df.loc[target_idx, stim_idx] = response_result

        return df


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

    ### ALLOPTICAL ANALYSIS - FOR ALL ROIs  # good progress on this, almost done reviewing
    def _makePhotostimTrialFluSnippets(self, plane_flu: np.ndarray, plane: int = 0,
                                       stim_frames: list = None) -> np.ndarray:  # base code copied from Vape's _makeFluTrials
        """
        Make Flu snippets timed on photostimulation, for each cell, for each stim instance. [cells x Flu frames x stims]  # TODO triple check order of this array's dimensions

        Inputs:
            plane_flu   - array of dff traces for all cells for this plane only
            plane       - imaging plane corresponding to plane_flu, default = 0 (for one plane datasets)
            stim_frames - optional, if provided then only use these photostim frames to collect photostim_array
        Outputs:
            photostim_array     - dFF peri-photostim Flu array [cell x Flu frames x trial]
        """

        print('\n\- Collecting peri-stim traces ...')

        trial_array = []
        _stims = self.stim_start_frames if stim_frames is None else stim_frames

        assert plane_flu.ndim == 2, 'plane_flu needs to be of ndim: 2'
        assert _stims == self.stim_start_frames, "stims not found in the stim frames list of this plane"

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

    def collectPhotostimResponses(self, photostimFluArray):
        """
        TODO docstring

        :param photostimFluArray:
        :return:
        """
        # create parameters, slices, and subsets for making pre-stim and post-stim arrays to use in stats comparison
        # test_period = self.pre_stim_response_window / 1000  # sec
        # self.test_frames = int(self.imparams.fps * test_period)  # test period for stats

        # mean pre and post stimulus (within post-stim response window) flu trace values for all cells, all trials
        self.__analysis_array = photostimFluArray
        self.__pre_array = np.mean(self.__analysis_array[:, self.pre_stim_test_slice, :],
                                   axis=1)  # [cells x prestim frames] (avg'd taken over all stims)
        self.__post_array = np.mean(self.__analysis_array[:, self.post_stim_test_slice, :],
                                    axis=1)  # [cells x poststim frames] (avg'd taken over all stims)

        # Vape's version for collection photostim response amplitudes
        # calculate amplitude of response for all cells, all trials
        all_amplitudes = self.__post_array - self.__pre_array

        df = pd.DataFrame(index=range(self.Suite2p.n_units), columns=self.stim_start_frames, data=all_amplitudes)

        return df

    def _allCellsPhotostimResponsesAnndata(self, photostimResponseAmplitudes: pd.DataFrame):  # NOT TESTED!
        """
        Creates annotated data (see anndata library) object based around the Ca2+ matrix of the imaging trial.

        """

        # try:
        # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
        # build dataframe for obs_meta from suite2p stat information
        obs_meta = pd.DataFrame(
            columns=['original_index', 'photostim_target', 'photostim_exclusion_zone', 'prob_response',
                     'sig_responder'], index=range(self.Suite2p.n_units))
        for idx in obs_meta.index:
            obs_meta.loc[idx, 'original_index'] = self.Suite2p.stat[idx]['original_index']
        obs_meta.loc[:, 'photostim_target'] = self.targeted_cells if hasattr(self, 'targeted_cells') else None
        obs_meta.loc[:, 'photostim_exclusion_zone'] = self.exclude_cells if hasattr(self, 'exclude_cells') else None
        obs_meta.loc[:, 'prob_response'] = self.prob_response if hasattr(self, 'prob_response') else None
        obs_meta.loc[:, 'sig_responder'] = self.sig_units if hasattr(self, 'sig_units') else None

        # SETUP THE VARIABLES ANNOTATIONS TO USE IN anndata
        # build dataframe for var annot's from Paq file
        var_meta = pd.DataFrame(index=self.Paq.paq_channels, columns=self.stim_start_frames)
        for fr_idx, fr in enumerate(self.stim_start_frames):
            for index in [*self.Paq.sparse_paq_data]:
                var_meta.loc[index, fr] = self.Paq.sparse_paq_data[index][fr_idx]
        photostimResponseAmplitudes.columns = var_meta.columns
        # var_meta.columns = photostimResponseAmplitudes.columns

        # BUILD LAYERS TO ADD TO anndata OBJECT
        layers = {'singleTrialSignificance': np.empty_like(photostimResponseAmplitudes)
                  # PLACEHOLDER NOT IMPLEMENTED YET
                  }

        print(f"\n\----- CREATING annotated data object for photostim responses using AnnData:")
        adata = AnnotatedData(X=np.asarray(photostimResponseAmplitudes), obs=obs_meta, var=var_meta.T, layers=layers)

        print(f"\t{adata}")
        return adata

        # except Exception:
        #     raise Warning("could not create anndata. anndata creation only available if .photostimResponseAmplitudes have been collected (run .photostimProcessingAllCells())')")

    def photostimProcessingAllCells(self, plane: int = 0):  # NOTE: not setup for multi-plane imaging processing yet...
        """
        Take dfof trace for entire timeseries and break it up in to individual trials, calculate
        the mean amplitudes of response and statistical significance across all trials

        Inputs:
            plane             - imaging plane n
        """
        print('\n----------------------------------------------------------------')
        print('running trial Processing for all cells ')
        print('----------------------------------------------------------------')

        # make trial arrays from dff data shape: [cells x stims x frames]
        if hasattr(self, 'Suite2p'):
            photostimFluArray = self._makePhotostimTrialFluSnippets(plane_flu=self.normalize_dff(self.Suite2p.raw))
            photostimResponseAmplitudes = self.collectPhotostimResponses(photostimFluArray)

            ## create new anndata object for storing measured photostim responses from data, with other relevant data
            photostim_responses_adata = self._allCellsPhotostimResponsesAnndata(
                photostimResponseAmplitudes=photostimResponseAmplitudes)

            return photostimFluArray, photostimResponseAmplitudes, photostim_responses_adata
        else:
            NotImplementedError('Photostim processing cannot be performed without Suite2p data.')

    def statisticalProcessingAllCells(self):
        """Runs statistical processing on photostim response arrays"""
        from packerlabimaging.processing.stats import AllOpticalStats

        self.wilcoxons = AllOpticalStats.runWilcoxonsTest(array1=self.__pre_array, array2=self.__post_array)
        self.sig_units = AllOpticalStats.sigTestAvgResponse(self=self, p_vals=self.wilcoxons, alpha=0.1)

    ## NOT REVIEWED FOR USAGE YET
    def _probResponse(self, plane,
                      trial_sig_calc):  ## FROM VAPE'S CODE, TODO NEED TO CHOOSE BETWEEN HERE AND BOTTOM RELIABILITY CODE
        """
        Calculate the response probability, i.e. proportion of trials that each cell responded on

        Inputs:
            plane          - imaging plane n
            trial_sig_calc - indicating which calculation was used for significance testing ('dff'/'dfsf')
        """
        n_trials = self.n_trials

        # get the number of responses for each across all trials
        if trial_sig_calc == 'dff':
            num_respond = np.array(self.trial_sig_dff[plane])  # trial_sig_dff is [plane][cell][trial]
        elif trial_sig_calc == 'dfsf':
            num_respond = np.array(self.trial_sig_dfsf[plane])

        # calculate the proportion of all trials that each cell responded on
        self.prob_response.append(np.sum(num_respond, axis=1) / n_trials)

    def cellStaProcessing(self, test='t_test'):
        """
        TODO docstring

        :param test:
        """
        if self.stim_start_frames:

            # this is the key parameter for the sta, how many frames before and after the stim onset do you want to use
            self.pre_frames = int(np.ceil(self.imparams.fps * 0.5))  # 500 ms pre-stim period
            self.post_frames = int(np.ceil(self.imparams.fps * 3))  # 3000 ms post-stim period

            # ls of cell pixel intensity values during each stim on each trial
            self.all_trials = []  # ls 1 = cells, ls 2 = trials, ls 3 = dff vector

            # the average of every trial
            self.stas = []  # ls 1 = cells, ls 2 = sta vector

            self.all_amplitudes = []
            self.sta_amplitudes = []

            self.t_tests = []
            self.wilcoxons = []

            for plane in range(self.imparams.n_planes):

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

                    for stim in self.stim_start_frames:
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
                        avg_post_end = avg_post_start + int(
                            np.ceil(self.imparams.fps * 0.5))  # post-stim period of 500 ms

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
        """
        TODO docstring

        :param test:
        """
        self.sta_sig = []

        for plane in range(self.imparams.n_planes):

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
        """
        TODO docstring
        """
        self.single_sig = []  # single trial significance value for each trial for each cell in each plane

        for plane in range(self.imparams.n_planes):

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
    def run_stamm_nogui(self, numDiffStims, startOnStim, everyXStims, preSeconds=0.75, postSeconds=1.25):
        """
        run STAmoviemaker for current trial

        :param numDiffStims:
        :param startOnStim:
        :param everyXStims:
        :param preSeconds:
        :param postSeconds:
        """
        qnap_path = os.path.expanduser('/home/pshah/mnt/qnap')

        # data path
        movie_path = self.tiff_path
        sync_path = self._paq_path

        # stamm save path
        stam_save_path = os.path.join(qnap_path, 'Analysis', self.metainfo['date'], 'STA_Movies',
                                      '%s_%s_%s' % (self.metainfo['date'],
                                                    self.metainfo['expID'],
                                                    self.metainfo['trialID']))
        os.makedirs(stam_save_path, exist_ok=True)

        ##
        assert os.path.exists(stam_save_path)

        print('QNAP_path:', qnap_path,
              '\ndata path:', movie_path,
              '\nsync path:', sync_path,
              '\nSTA movie save s2pResultsPath:', stam_save_path)

        # define STAmm parameters
        frameRate = int(self.imparams.fps)

        arg_dict = {'moviePath': movie_path,  # hard-code this
                    'savePath': stam_save_path,
                    'syncFrameChannel': "frame_clock",
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
        # plot_single_tiff(img, frame_num=0)

    def _targetSpread(self):
        '''
        Find the mean Euclidean distance of responding targeted cells (trial-wise and trial average)
        '''
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

    ####                                                                    // end

if __name__ == '__main__':
    LOCAL_DATA_PATH = '/Users/prajayshah/data/oxford-data-to-process/'
    REMOTE_DATA_PATH = '/home/pshah/mnt/qnap/Data/'
    BASE_PATH = LOCAL_DATA_PATH

    ExperimentMetainfo = {
        'dataPath': f'{BASE_PATH}/2020-12-19/2020-12-19_t-013/2020-12-19_t-013_Cycle00001_Ch3.tif',
        'saveDir': f'{BASE_PATH}/2020-12-19/',
        'expID': 'RL109',
        'comment': 'two photon imaging + alloptical trials',
    }

    expobj = Experiment(**ExperimentMetainfo)


    def alloptical_trial_fixture():
        initialization_dict = {'naparm_path': f'{BASE_PATH}/2020-12-19/photostim/2020-12-19_RL109_ps_014/',
                               'dataPath': f'{BASE_PATH}/2020-12-19/2020-12-19_t-013/2020-12-19_t-013_Cycle00001_Ch3.tif',
                               'saveDir': f'{BASE_PATH}/2020-12-19/',
                               'date': '2020-12-19',
                               'trialID': 't-013',
                               'expID': 'RL109',
                               'expGroup': 'all optical trial with LFP',
                               'comment': ''}

        return initialization_dict


    def test_AllOpticalClass(alloptical_trial_fixture):
        from packerlabimaging.processing.imagingMetadata import PrairieViewMetadata
        from packerlabimaging.main.paq import PaqData

        paqs_loc = f'{BASE_PATH}/2020-12-19/2020-12-19_RL109_013.paq'  # path to the .paq files for the selected trials
        dataPath = alloptical_trial_fixture['dataPath']

        # parses imaging system data
        imparams = PrairieViewMetadata(pv_xml_dir=os.path.dirname(dataPath), microscope='Bruker 2pPlus')

        # sets the stim start frames
        tmdata = PaqData.paqProcessingAllOptical(paq_path=paqs_loc, frame_channel='frame_clock',
                                                 stim_channel='markpoints2packio')

        # create the trial
        aotrial = AllOpticalTrial(imparams=imparams, tmdata=tmdata, **alloptical_trial_fixture)
        return aotrial


    idict = alloptical_trial_fixture()
    aotrial = test_AllOpticalClass(idict)
