# TODO update all 2p stim related attr's to naparm submodule
import glob
import os
import time
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from imagingplus import TwoPhotonImaging
from imagingplus.main.subcore import ImagingMetadata, CellAnnotations
from imagingplus.main.core import Experiment, ImagingTrial
from imagingplus.processing.paq import PaqData
from imagingplus.processing.naparm import Targets

# %%
from imagingplus.processing.anndata import AnnotatedData

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)

PLANE = 0
BADFRAMESLOC = '/home/pshah/Documents/code/imagingplus/tests/'

prestim_sec: float = 1.0  #: length of pre stim trace collected (secs)
poststim_sec: float = 3.0  #: length of post stim trace collected (secs)
pre_stim_response_window: float = 0.500  #: time window for collecting pre-stim measurement (units: msec)
post_stim_response_window: float = 0.500  #: time window for collecting post-stim measurement (units: msec)


# todo think about creating a custom class to hold directly taken coords targets imaging cellsdata - would work well for SLM targets, might even be able to extend in the naparm--> Targets class

class AllOpticalTrial(TwoPhotonImaging):
    """All Optical Experimental Data Analysis Workflow."""

    def __init__(self, naparm_path, dataPath: str, saveDir: str, date: str, trialID: str, expID: str,
                 expGroup: str = '', comment: str = '', imparams: ImagingMetadata = None, cells: CellAnnotations = None,
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
        self.__prestim_response_window: Union[
            float, int] = 0.500  #: time window for collecting pre-stim measurement (units: msec)
        self.__poststim_response_window: Union[
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
        self.wilcoxons = []  #:
        self.sig_units = None  #:
        self.sta_sig = []  # based on t-test between dff test periods
        self.sta_sig_nomulti = []  # as above, no multiple comparisons correction
        ########

        # initializing cellsdata processing, cellsdata analysis and/or results associated attr's

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

        # 6) Collect Flu traces from Suite2p ROIs:
        #   create:
        #       1) array of dFF pre + post stim Flu snippets for each stim and cell [num cells x num peri-stim frames collected x num stims]
        #       2) photostim response amplitudes in a dataframe for each cell and each photostim
        #       3) save photostim response amplitudes to AnnotatedData
        self.photostimProcessingAllCells()

        # save final object
        self.save()
        print(f'\----- CREATED AllOpticalTrial cellsdata object for {self.t_series_name}')

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
        TODO add parameters
        :param naparm_path:
        :param imaging_trial:
        :return:
        """

        initialization_dict = {'naparm_path': naparm_path, 'dataPath': imaging_trial.dataPath,
                               'saveDir': imaging_trial.saveDir,
                               'expID': imaging_trial.expID, 'group': imaging_trial.expGroup,
                               'comment': imaging_trial.comment}

        aotrial = cls(**initialization_dict)

        return aotrial

    @property
    def twopstim_path(self):
        """path to folder containing photostimulation protocols output by NAPARM
        TODO add parameters
        :return:
        """
        if self.twopstim:
            return self.twopstim.path

    @property
    def prestim_sec(self):
        """length of pre stim trace collected (secs)
        TODO add parameters
        :return:
        """
        return self.__prestim_sec

    @prestim_sec.setter
    def prestim_sec(self, val):
        """
    TODO fill explanation and add parameters
        :param val:
        """
        assert type(val) == int or type(val) == float, 'can only set prestim_sec with int or float'
        self.__prestim_sec = val

    @property
    def poststim_sec(self):
        """length of post stim trace collected (secs)
        TODO add parameters
        :return:
        """
        return self.__poststim_sec

    @poststim_sec.setter
    def poststim_sec(self, val):
        """
        TODO fill explanation and add parameters
        :param val:
        """
        assert type(val) == int or type(val) == float, 'can only set poststim_sec with int or float'
        self.__prestim_sec = val

    @property
    def prestim_response_window(self):
        """time window for collecting pre-stim measurement (units: msec)
        TODO add parameters
        :return:
        """
        return self.__prestim_response_window

    @prestim_response_window.setter
    def prestim_response_window(self, val):
        """
        TODO fill explanation and add parameters
        :param val:
        """
        assert type(val) == int or type(val) == float, 'can only set prestim_response_window with int or float'
        self.__prestim_response_window = val

    @property
    def prestim_response_frames_window(self):
        """num frames for measuring Flu trace before each photostimulation trial during photostim response measurement (units: frames)
        TODO add parameters
        :return:
        """
        return int(self.imparams.fps * self.prestim_response_window)

    @property
    def poststim_response_window(self):
        """time window for collecting post-stim measurement (units: msec)
        TODO add parameters
        :return:
        """
        return self.__poststim_response_window

    @poststim_response_window.setter
    def poststim_response_window(self, val):
        """
        TODO fill explanation and add parameters
        :param val:
        """
        assert type(val) == int or type(val) == float, 'can only set poststim_response_window with int or float'
        self.__poststim_response_window = val

    @property
    def poststim_response_frames_window(self):
        """num frames for measuring Flu response after each photostimulation trial during photostim response measurement (units: frames)
        TODO add parameters
        :return:
        """
        return int(self.imparams.fps * self.poststim_response_window)

    @property
    def prestim_frames(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)
        TODO add parameters
        :return:
        """
        return int(self.__prestim_sec * self.imparams.fps)

    @property
    def poststim_frames(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)
        TODO add parameters
        :return:
        """
        return int(self.__poststim_sec * self.imparams.fps)

    @property
    def n_stims(self):
        #: TODO set property using anndata responses shape property assigning from array: cells x Flu frames x # of photostim trials
        """number of photostimulation trials
        TODO add parameters
        :return:
        """
        return

    @property
    def stim_start_frames(self):
        """Imaging frames corresponding to start of photostimulation trials.
        TODO add parameters
        :return:
        """
        # todo use anndata for getting this
        return []

    @property
    def photostim_frames(self):
        """Imaging frames during photostimulation trials.
        TODO add parameters
        :return:
        """
        # todo just use anndata for getting this
        return []

    @property
    def stim_duration_frames(self):
        """Duration of photostimulation as number of frames.
        Note: Removing 1 more frame than the stim duration, as the stim isn't perfectly aligned with the start of the imaging frame
        TODO add parameters
        :return:
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
    def prestim_test_slice(self):
        """slice representing prestim response frames
        TODO add parameters
        :return:
        """
        return np.s_[self.prestim_frames - self.prestim_response_frames_window: self.prestim_frames]

    @property
    def poststim_test_slice(self):
        """slice representing poststim response frames
        TODO add parameters
        :return:
        """
        stim_end = self.prestim_frames + self.stim_duration_frames
        return np.s_[stim_end: stim_end + self.poststim_response_frames_window]

    def photostimProcessing(self, naparm_path):
        """
        Processing alloptical trial photostimulation protocol.
        - parse NAPARM protocol and SLM Targets information under .twopstim
        - collect stimulation timing synchronized to imaging cellsdata frame timings
        TODO add parameters
        :param naparm_path:
        :return:

        """

        self.twopstim = Targets(naparm_path=naparm_path, frame_x=self.imparams.frame_x, frame_y=self.imparams.frame_y, pix_sz_x=self.imparams.pix_sz_x, pix_sz_y=self.imparams.pix_sz_y)

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

    #### TODO review attr's and write docs from the following functions: // start

    ### ALLOPTICAL PROCESSING OF TRACES
    ## ... no methods determined here yet...

    ### ALLOPTICAL ANALYSIS - FOCUS ON SLM TARGETS RELATED METHODS - TEST all methods below
    # todo test
    def getTargetsStimTraceSnippets(self, targets_idx: Union[list, str] = 'all', pre_stim: Union[float, int] = 0.5,
                                    post_stim: Union[float, int] = 4.0, stims: list = None):
        """
        Collect photostimulation timed snippets of signal from selected targets.

        TODO add parameters
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

        pre_stim = self.prestim_frames
        post_stim = self.poststim_frames

        if stims is None:
            stim_timings = self.stim_start_frames
        else:
            stim_timings = stims

        process = 'trace dFF'
        if targets_idx == 'all':
            data_to_process = self.dFF_SLMTargets
        else:
            assert type(targets_idx) == list, 'provide targets_idx as list on cell indexes to collect.'
            data_to_process = np.asarray([self.dFF_SLMTargets[idx] for idx in targets_idx])

        num_cells = data_to_process.shape[0]
        print(f'collecting stim traces for {num_cells} cells ', num_cells)

        # collect photostim timed trace snippets traces of photostim targets
        flu = np.asarray([data_to_process[:, stim - pre_stim: stim + self.stim_duration_frames + post_stim] for stim in
                          stim_timings])

        print(f"shape photostim. trials trace snippets array: {flu.shape}")
        return flu

    # todo test
    def findTargetedCells(self, plot: bool = True):
        """finding s2p cell ROIs that were also SLM targets (or more specifically within the target areas as specified by _findTargetAreas - include 15um radius from center coordinate of spiral)
        Make a binary mask of the targets and multiply by an image of the cells
        to find cells that were targeted

        --- LAST UPDATED NOV 6 2021 - copied over from Vape ---

        TODO add parameters
        :param plot:
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

        # add targets classification as observations annotation to .cellsdata anndata
        self.data.add_obs(obs_name='photostim_class', values=targets)

        self.cells.cellsdata['photostim_class'] = targets

        print(f"\t|- Number of non-target ROIs: {len(self.cells.cellsdata['photostim_class'] == 'non-target')}")

    # todo code and test
    def _makePhotostimTrialFluSnippets(self, plane_flu: np.ndarray,
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
            trial_frames = np.s_[stim - self.prestim_frames: stim + self.poststim_frames]

            # use trial frames to extract this trial for every cell
            flu_trial = plane_flu[:, trial_frames]
            flu_trial_len = self.prestim_frames + self.stim_duration_frames + self.poststim_frames

            # todo test if needed or better way to implement: catch timeseries which ended in the middle of an ongoing photostim instance
            if flu_trial.shape[1] == flu_trial_len:
                # only append trials of the correct length - will catch corrupt/incomplete cellsdata and not include
                if len(trial_array) == 0:
                    trial_array = flu_trial
                else:
                    trial_array = np.dstack((trial_array, flu_trial))
            else:
                print('**incomplete trial detected and not appended to trial_array**', end='\r')

        print(f'\nFinished collecting peri-stim traces, out shape: {trial_array.shape}')

        return trial_array

    # todo test
    def _normalize_snippets_prestim(self, snippets: np.ndarray = None):
        """
        Normalize each trace snippet to pre-stim period.

        TODO add parameters
        :param snippets:
        :return:
        :return:
        """

        snippets = snippets if snippets else self.targets_snippets
        num_cells = snippets.shape[0]

        targets_dff = np.zeros([num_cells, len(self.stim_start_frames),
                                self.prestim_frames + self.stim_duration_frames + self.poststim_frames])

        # create pre-stim period mean substracted photostim trials trace snippets
        for j, traces in enumerate(snippets):
            for i, trace in enumerate(traces):
                mean_pre = np.mean(trace[0: self.prestim_frames])
                trace_dff = (trace - mean_pre)
                targets_dff[j, i] = trace_dff
        print(f"shape photostim. trials trace snippets array: {targets_dff.shape}")
        return targets_dff

    # todo test
    def calculatePhotoresponses(self, snippets: np.ndarray, stims_to_use: Union[list, str] = 'all') -> pd.DataFrame:
        """
        Calculations of responses (post-stim - pre-stim) to photostimulation of SLM Targets of the provided snippets array.

        TODO add parameters
        :param snippets:
        :param stims_to_use: ls of stims to retrieve photostim trial dFF responses
        :return:
        """

        if stims_to_use is 'all':
            stims_to_use = range(len(self.stim_start_frames))
            stims_idx = [self.stim_start_frames.index(stim) for stim in stims_to_use]
        else:
            stims_idx = [self.stim_start_frames.index(stim) for stim in stims_to_use]

        # method 1) initializing pandas df that collects responses of stimulations
        d = {}
        for stim in stims_idx:
            d[stim] = [None] * snippets.shape[0]
        df = pd.DataFrame(d, index=range(snippets.shape[0]))  # population dataframe

        # calculate photostimulation responses
        cell_ids = df.index
        for target_idx in range(len(cell_ids)):
            responses = []
            for stim_idx in stims_idx:
                dff_trace = snippets[target_idx][stim_idx]
                response_result = np.mean(dff_trace[self.prestim_frames + self.stim_duration_frames + 1:
                                                    self.prestim_frames + self.stim_duration_frames +
                                                    self.poststim_response_frames_window])  # calculate the dF over pre-stim mean F response within the response window
                responses.append(np.round(response_result, 2))

                df.loc[target_idx, stim_idx] = response_result

        # method 2) alternate method for calculating photostim responses
        self.__analysis_array = snippets
        self.__pre_array = np.mean(self.__analysis_array[:, self.prestim_test_slice, :],
                                   axis=1)  # [cells x prestim frames] (avg'd taken over all stims)
        self.__post_array = np.mean(self.__analysis_array[:, self.poststim_test_slice, :],
                                    axis=1)  # [cells x poststim frames] (avg'd taken over all stims)

        # Vape's version for collection photostim response amplitudes
        # calculate amplitude of response for all cells, all trials
        all_amplitudes = self.__post_array - self.__pre_array

        df = pd.DataFrame(index=range(self.Suite2p.n_units), columns=self.stim_start_frames, data=all_amplitudes)

        # todo compare two methods

        return df

    # todo code up
    def _TargetsPhotostimResponsesAnndata(self):
        """
        Create an anndata table for photostimulation responses of Targets.
        TODO add parameters
        :return:
        """

    # todo test
    def _CellsPhotostimResponsesAnndata(self, photostimResponseAmplitudes: pd.DataFrame):
        """
        Creates annotated cellsdata (see anndata library) object based around the Ca2+ matrix of the imaging trial.

        TODO add parameters
        :param photostimResponseAmplitudes:
        :return:

        """

        # try:
        # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
        # build dataframe for obs_meta from suite2p stat information
        obs_meta = self.cells.cellsdata

        # SETUP THE VARIABLES ANNOTATIONS TO USE IN anndata
        # build dataframe for var annot's from Paq file
        var_meta = self.tmdata.data
        # var_meta.columns = photostimResponseAmplitudes.columns

        print(f"\n\----- CREATING annotated cellsdata object for photostim responses using AnnData:")
        photostimResponseAmplitudes.columns = var_meta.columns
        adata = AnnotatedData(X=np.asarray(photostimResponseAmplitudes), obs=obs_meta, var=var_meta.T)

        print(f"\t{adata}")
        return adata

    def photostimProcessingAllCells(self, plane: int = 0):  # NOTE: not setup for multi-plane imaging processing yet...
        """
        Take dfof trace for entire timeseries and break it up in to individual trials, calculate
        the mean amplitudes of response and statistical significance across all trials

        TODO add parameters
        Inputs:
        :param plane:             - imaging plane n
        """
        print('\n----------------------------------------------------------------')
        print('running trial Processing for all cells ')
        print('----------------------------------------------------------------')

        # make trial arrays from dff cellsdata shape: [cells x stims x frames]
        if self.Suite2p != None:
            photostimFluArray = self._makePhotostimTrialFluSnippets(plane_flu=self.normalize_dff(self.imdata.imdata))
            photostimResponseAmplitudes = self.calculatePhotoresponses(photostimFluArray)

            ## create new anndata object for storing measured photostim responses from cellsdata, with other relevant cellsdata
            self.data = self._CellsPhotostimResponsesAnndata(
                photostimResponseAmplitudes=photostimResponseAmplitudes)

            # extend annotated imaging cellsdata object with imaging frames in photostim and stim_start_frames as additional keys in vars
            __frames_in_stim = [False] * self.imparams.n_frames
            __stim_start_frame = [False] * self.imparams.n_frames
            for frame in self.twopstim.photostim_frames: __frames_in_stim[frame] = True
            for frame in self.twopstim.stim_start_frames: __stim_start_frame[frame] = True
            self.data.add_var(var_name='photostim_frame', values=__frames_in_stim)
            self.data.add_var(var_name='stim_start_frame', values=__stim_start_frame)

            self.photostimFluArray, self.photostimResponseAmplitudes = photostimFluArray, photostimResponseAmplitudes
        else:
            print('Photostim. processing of fluorescence signal responses cannot be performed without Suite2p cellsdata.')


    def statisticalProcessingAllCells(self):
        """Runs statistical processing on photostim response arrays"""

        from imagingplus._archive.stats import runWilcoxonsTest
        self.wilcoxons = runWilcoxonsTest(array1=self.__pre_array, array2=self.__post_array)
        from imagingplus._archive.stats import sigTestAvgResponse
        self.sig_units = sigTestAvgResponse(trial=self, p_vals=self.wilcoxons, alpha=0.1)

    def staSignificance(self, test):
        """
        TODO docstring
        TODO add parameters
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
        TODO fill explanation
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
        TODO add parameters
        :param numDiffStims:
        :param startOnStim:
        :param everyXStims:
        :param preSeconds:
        :param postSeconds:
        """
        qnap_path = os.path.expanduser('/home/pshah/mnt/qnap')

        # cellsdata path
        movie_path = self.data_path
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
              '\ncellsdata path:', movie_path,
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
        """
        :return:
        """
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
        """
        :param alloptical_trial_fixture:
        :return:
        """
        from imagingplus.processing.imagingMetadata import PrairieViewMetadata
        from imagingplus.processing.paq import PaqData

        paqs_loc = f'{BASE_PATH}/2020-12-19/2020-12-19_RL109_013.paq'  # path to the .paq files for the selected trials
        dataPath = alloptical_trial_fixture['dataPath']

        # parses imaging system cellsdata
        imparams = PrairieViewMetadata(pv_xml_dir=os.path.dirname(dataPath), microscope='Bruker 2pPlus')

        # sets the stim start frames
        tmdata = PaqData.paqProcessingAllOptical(paq_path=paqs_loc, frame_channel='frame_clock',
                                                 stim_channel='markpoints2packio')

        # create the trial
        aotrial = AllOpticalTrial(imparams=imparams, tmdata=tmdata, **alloptical_trial_fixture)
        return aotrial


    idict = alloptical_trial_fixture()
    aotrial = test_AllOpticalClass(idict)
