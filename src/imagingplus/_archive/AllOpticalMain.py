# TODO update all 2p stim related attr's to naparm submodule
import glob
import os
import signal
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import tifffile as tf

from packerlabimaging import TwoPhotonImaging
from packerlabimaging._archive.classes import ImagingMetadata, TemporalData, ImagingTrial, CellAnnotations, \
    Experiment
from packerlabimaging.utils.utils import convert_to_8bit
from packerlabimaging.processing.naparm import Targets
from packerlabimaging.utils.classes import UnavailableOptionError
# %%
from packerlabimaging.processing.anndata import AnnotatedData

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)

PLANE = 0
BADFRAMESLOC = '/home/pshah/Documents/code/imagingplus/tests/'


def getTargetImage():  # TODO write this function and probably add to the Targets class
    pass


class AllOpticalTrial(TwoPhotonImaging):
    """All Optical Experimental Data Analysis Workflow."""
    prestim_sec: float = 1.0  #: length of pre stim trace collected (in frames)
    poststim_sec: float = 3.0  #: length of post stim trace collected (in frames)
    pre_stim_response_window: float = 0.500  #: time window for collecting pre-stim measurement (units: msec)
    post_stim_response_window: float = 0.500  #: time window for collecting post-stim measurement (units: msec)

    def __init__(self, naparm_path, dataPath: str, saveDir: str, date: str, trialID: str, expID: str,
                 expGroup: str = '',
                 comment: str = '', imparams: ImagingMetadata = None, cells: CellAnnotations = None,
                 tmdata: TemporalData = None):

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


        # Paq file attr's
        self.paq_rate: int = -1  # PackIO acquisition rate for .Paq file
        self.frame_start_times: list = [None]  # Paq clock timestamps of the first imaging acquisition frame of t-series
        self.frame_end_times: list = [None]  # Paq clock timestamps of the last imaging acquisition frame of t-series

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

        # initializing cellsdata processing, cellsdata analysis and/or results associated attr's
        self.n_trials = None  # number of photostimulation trials TODO change to assigning from array: cells x Flu frames x # of photostim trials

        # PHOTOSTIM SLM TARGETS
        # TODO add attr's related to numpy array's and pandas dataframes for photostim trials - SLM targets
        self.responses_SLMtargets = []  # dF/prestimF responses for all SLM targets for each photostim trial
        self.responses_SLMtargets_tracedFF = []  # poststim dFF - prestim dFF responses for all SLM targets for each photostim trial

        # ALL CELLS (from suite2p ROIs)
        # TODO add attr's related to numpy array's and pandas dataframes for photostim trials - suite2p ROIs

        # FUNCTIONS TO RUN AFTER init's of ALL ATTR'S

        # 3) process 2p stim protocol
        # set naparm path
        self.__naparm_path = naparm_path if os.path.exists(naparm_path) else FileNotFoundError(
            f"naparm path not found, naparm_path: {naparm_path}")
        self.Targets, self.stim_duration_frames = self._stimProcessing(protocol='naparm')

        # 4) determine bad frames in imaging cellsdata that correspond to photostim frames
        self.photostim_frames = self._find_photostim_add_bad_framesnpy()

        # 5) collect Flu traces from SLM targets
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

        # extend annotated imaging cellsdata object with imaging frames in photostim and stim_start_frames as additional keys in vars
        __frames_in_stim = [False] * self.imparams.n_frames
        __stim_start_frame = [False] * self.imparams.n_frames
        for frame in self.photostim_frames: __frames_in_stim[frame] = True
        for frame in self.stim_start_frames: __stim_start_frame[frame] = True
        self.data.add_var(var_name='photostim_frame', values=__frames_in_stim)
        self.data.add_var(var_name='stim_start_frame', values=__stim_start_frame)

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
        """

        initialization_dict = {'naparm_path': naparm_path, 'dataPath': imaging_trial.dataPath,
                               'saveDir': imaging_trial.saveDir,
                               'expID': imaging_trial.expID, 'group': imaging_trial.expGroup,
                               'comment': imaging_trial.comment}

        aotrial = cls(**initialization_dict)

        return aotrial

    @property
    def naparm_path(self):
        """path to folder containing photostimulation protocols output by NAPARM"""
        if self.__naparm_path[-1] == '/':
            return self.__naparm_path
        else:
            return self.__naparm_path + '/'

    @property
    def pre_stim_response_frames_window(self):
        """num frames for measuring Flu trace before each photostimulation trial during photostim response measurement (units: frames)"""
        return int(self.imparams.fps * self.pre_stim_response_window_msec)

    @property
    def pre_stim_frames(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)"""
        return int(self.prestim_sec * self.imparams.fps)

    @property
    def post_stim_frames(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)"""
        return int(self.poststim_sec * self.imparams.fps)

    @property
    def post_stim_response_frames_window(self):
        """num frames for collecting Flu trace after each photostimulation trial (units: frames)"""
        return int(
            self.imparams.fps * self.post_stim_response_window_msec)  # length of the post stim response test window (in frames)

    @property
    def timeVector(self):
        """vector of frame times in milliseconds rather than frames"""
        time_vector = np.linspace(-self.prestim_sec, self.poststim_sec, self.imparams.n_frames)
        return time_vector

    # ALLOPTICAL EXPERIMENT PHOTOSTIM PROTOCOL PROCESSING ##############################################################

    def _paqProcessingAllOptical(self, stim_channel: str):
        """
        TODO need to update!

        Process .paq file with relation to all optical processing. In particular, this returns frame numbers timed to
        individual photostimulation trial synced to the specified stim_channel.

        """
        paqdata, _, _ = self.Paq.paq_read()
        stim_start_frames, stim_start_times = self.Paq.paq_alloptical_stims(paq_data=paqdata,
                                                                            frame_clock=self.Paq.frame_times,
                                                                            stim_channel=stim_channel)
        return stim_start_frames

    def _stimProcessing(self, protocol: str = 'naparm'):
        _available_protocols = ['naparm']
        if protocol == 'naparm':
            targets = Targets(naparm_path=self.naparm_path, frame_x=self.imparams.frame_x,
                              frame_y=self.imparams.frame_y,
                              pix_sz_x=self.imparams.pix_sz_x, pix_sz_y=self.imparams.pix_sz_y)

            # correct the following based on txt file
            duration_ms = targets.stim_dur
            frame_rate = self.imparams.fps / self.imparams.n_planes
            duration_frames = np.ceil((duration_ms / 1000) * frame_rate)
            stim_duration_frames = int(duration_frames)

            return targets, stim_duration_frames
        else:
            raise UnavailableOptionError(
                f"{protocol} is not available for 2p stim processing. Available options are: {_available_protocols}")

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
        targ_img = np.zeros([self.imparams.frame_x, self.imparams.frame_y], dtype='uint16')
        target_areas = np.array(self.Targets.target_areas)
        targ_img[target_areas[:, :, 1], target_areas[:, :, 0]] = 1

        # make an image of every cell area, filled with the index of that cell
        cell_img = np.zeros_like(targ_img)

        cell_y = np.array(self.Suite2p.cell_x)
        cell_x = np.array(self.Suite2p.cell_y)

        for i, coord in enumerate(zip(cell_x, cell_y)):
            cell_img[coord] = i + 1

        # binary mask x cell image to get the cells that overlap with target areas
        targ_cell = cell_img * targ_img

        targ_cell_ids = np.unique(targ_cell)[1:] - 1  # correct the cell id due to zero indexing
        self.targeted_cells = np.zeros([self.Suite2p.n_units], dtype='bool')
        self.targeted_cells[targ_cell_ids] = True
        # self.s2p_cell_targets = [self.cell_id[i] for i, x in enumerate(self.targeted_cells) if x is True]  # get ls of s2p cells that were photostim targetted
        self.s2p_cell_targets = [self.Suite2p.cell_id[i] for i in
                                 np.where(self.targeted_cells)[0]]  # get ls of s2p cells that were photostim targetted

        self.n_targeted_cells = np.sum(self.targeted_cells)

        print('\t|- Search completed.')
        self.save()
        print('\t|- Number of targeted cells: ', self.n_targeted_cells)

        # IDENTIFYING S2P ROIs THAT ARE WITHIN THE EXCLUSION ZONES OF THE SLM TARGETS
        # make all target area coords in to a binary mask
        targ_img = np.zeros([self.imparams.frame_x, self.imparams.frame_y], dtype='uint16')
        target_areas_exclude = np.array(self.Targets.target_areas_exclude)
        targ_img[target_areas_exclude[:, :, 1], target_areas_exclude[:, :, 0]] = 1

        # make an image of every cell area, filled with the index of that cell
        cell_img = np.zeros_like(targ_img)

        cell_y = np.array(self.Suite2p.cell_x)
        cell_x = np.array(self.Suite2p.cell_y)

        for i, coord in enumerate(zip(cell_x, cell_y)):
            cell_img[coord] = i + 1

        # binary mask x cell image to get the cells that overlap with target areas
        targ_cell = cell_img * targ_img

        targ_cell_ids = np.unique(targ_cell)[1:] - 1  # correct the cell id due to zero indexing
        self.exclude_cells = np.zeros([self.Suite2p.n_units], dtype='bool')
        self.exclude_cells[targ_cell_ids] = True
        self.s2p_cells_exclude = [self.Suite2p.cell_id[i] for i in
                                  np.where(self.exclude_cells)[0]]  # get ls of s2p cells that were photostim targetted

        self.n_exclude_cells = np.sum(self.exclude_cells)

        print('\t|-Search completed.')
        self.save()
        print(f"\t|-Number of exclude Suite2p ROIs: {self.n_exclude_cells}")

        # define non targets from suite2p ROIs (exclude cells in the SLM targets exclusion - .s2p_cells_exclude)
        self.Suite2p.s2p_nontargets = [cell for cell in self.Suite2p.cell_id if
                                       cell not in self.s2p_cells_exclude]  ## exclusion of cells that are classified as s2p_cell_targets

        print(f"\t|-Number of good, s2p non-target ROIs: {len(self.Suite2p.s2p_nontargets)}")

        if plot:
            fig, ax = plt.subplots(figsize=[6, 6])

            targ_img = np.zeros([self.imparams.frame_x, self.imparams.frame_y], dtype='uint16')
            target_areas = np.array(self.Targets.target_areas)
            targ_img[target_areas[:, :, 1], target_areas[:, :, 0]] = 1
            ax.imshow(targ_img, cmap='Greys_r', zorder=0)
            ax.set_title('SLM targets areas')
            # for (x, y) in self.Targets.target_coords_all:
            #     ax.scatter(x=x, y=y, edgecolors='white', facecolors='none', linewidths=1.0)
            fig.show()

        # add targets classification as observations annotation to .cellsdata anndata
        self.data.add_observation(self.data, 'photostim_target', values=list(self.targeted_cells))

    def _find_photostim_add_bad_framesnpy(self):
        """finds all photostim frames and saves them into the bad_frames attribute for the exp object and saves bad_frames.npy"""
        print('\n\----- Finding photostimulation frames in imaging frames ...')
        print('# of photostim frames calculated per stim. trial: ', self.stim_duration_frames + 1)

        photostim_frames = []

        for j in self.stim_start_frames:
            for i in range(
                    self.stim_duration_frames + 1):  # usually need to remove 1 more frame than the stim duration, as the stim isn't perfectly aligned with the start of the imaging frame
                photostim_frames.append(j + i)

        # print(photostim_frames)
        print('\t|- Original # of frames:', self.imparams.n_frames, 'frames')
        print('\t|- # of Photostim frames:', len(self.photostim_frames), 'frames')
        print('\t|- Minus photostim. frames total:', self.imparams.n_frames - len(photostim_frames), 'frames')

        # if using Suite2p then add photostim frames to bad_frames.npy for current Experiment
        if self.Suite2p:
            if len(self.photostim_frames) > 0:
                bad_frames = self.Suite2p.bad_frames
                bad_frames.extend(self.photostim_frames)
                bad_frames = list(np.unique(bad_frames))
                print(
                    f'***Added a total of {len(self.photostim_frames)} photostim frames to bad_frames.npy at: {self.data_path_dir}/bad_frames.npy \n\t total bad_frames: {len(bad_frames)}')
                # f'***Saving a total of {len(photostim_frames)} photostim frames to bad_frames.npy at: {BADFRAMESLOC}/bad_frames.npy')  # TODO replace BADFRAMESLOC with self.pv_xml_dir
                np.save(f'{self.data_path_dir}/bad_frames.npy',
                        bad_frames)  # save to npy file and remember to move npy file to tiff folder before running with suite2p

        return photostim_frames

    #### TODO review attr's and write docs from the following functions: // start
    # used for creating tiffs that remove artifacts from alloptical experiments with photostim artifacts
    def rm_artifacts_tiffs(self, tiffs_loc, new_tiffs):
        """
        TODO docstring
        :param tiffs_loc:
        :param new_tiffs:
        """
        # make a new tiff file (not for suite2p) with the first photostim frame whitened, and save new tiff
        print('\n-----making processed photostim .tiff from:')
        tiff_path = tiffs_loc
        print(tiff_path)
        im_stack = tf.imread(tiff_path, key=range(self.n_frames))
        print('Processing experiment tiff of shape: ', im_stack.shape)

        frames_to_whiten = []
        for j in self.stim_start_frames:
            frames_to_whiten.append(j)

        # number of photostim frames with artifacts
        frames_to_remove = []
        for j in self.stim_start_frames:
            for i in range(0,
                           self.stim_duration_frames + 1):  # usually need to remove 1 more frame than the stim duration, as the stim isn't perfectly aligned with the start of the imaging frame
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

        for stim in range(self.Targets.n_groups):
            b = np.full_like(im_stack_2[0], fill_value=0)
            targets = self.Targets.target_areas[stim]
            for i in np.arange(len(targets)):
                for j in targets[i]:
                    b[j] = 5000

            all_stim_start_frames = []
            for stim_frame in self.stim_start_frames[stim::self.Targets.n_groups]:
                all_stim_start_frames.append(stim_frame)
            for frame in all_stim_start_frames:
                im_stack_2[frame - 1] = im_stack_2[frame - 1] + b

        im_stack_2 = np.delete(im_stack_2, self.photostim_frames, axis=0)

        print('After delete shape targetcells', im_stack_2.shape)

        save_path = (new_tiffs + '_targetcells.tif')
        tf.imwrite(save_path, im_stack_2, photometric='minisblack')

        print('done saving to: ', save_path)

        del im_stack_2
        del im_stack

    def s2pMasks(self, s2p_path, cell_ids):
        """
        Returns arrays that adds targets images to suite2p images.

        :param s2p_path:
        :param cell_ids:
        :return:
        """
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
            if n in self.s2p_cell_targets:
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

        return mask_img, targets_s2p_img

    def s2pMaskStack(self, pkl_list, s2p_path, parent_folder, force_redo: bool = False):
        """makes a TIFF stack with the s2p mean image, and then suite2p ROI masks for all cells detected,
        target cells, and also SLM targets as well?

        :param pkl_list:
        :param s2p_path:
        :param parent_folder:
        :param force_redo:
        """

        for pkl in pkl_list:
            print('Retrieving s2p masks for:', pkl, end='\r')

            # with open(pkl, 'rb') as f:
            #     self = pickle.load(f)

            # ls of cell ids to filter s2p masks by
            # cell_id_list = [ls(range(1, 99999)),  # all
            #                 self.photostim_r.cell_id[0],  # cells
            #                 [self.photostim_r.cell_id[0][i] for i, b in enumerate(self.photostim_r.cell_s1[0]) if
            #                  b == False],  # s2 cells
            #                 [self.photostim_r.cell_id[0][i] for i, b in enumerate(self.photostim_r.is_target) if
            #                  b == 1],  # pr cells
            #                 [self.photostim_s.cell_id[0][i] for i, b in enumerate(self.photostim_s.is_target) if
            #                  b == 1],  # ps cells
            #                 ]
            #
            cell_ids = self.Suite2p.cell_id

            # empty stack to fill with images
            stack = np.empty((0, self.imparams.frame_y, self.imparams.frame_x), dtype='uint8')

            s2p_path = s2p_path

            # mean image from s2p
            mean_img = self.Suite2p.s2pMeanImage(s2p_path)
            mean_img = np.expand_dims(mean_img, axis=0)
            stack = np.append(stack, mean_img, axis=0)

            # mask images from s2p
            mask_img, targets_s2p_img = self.s2pMasks(s2p_path=s2p_path, cell_ids=cell_ids)
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
            #     if all(s in file for s in ['AvgImage', self.photostim_r.tiff_path.split('/')[-1]]):
            #         pr_sta_img = tf.imread(os.path.join(stam_save_path, file))
            #         pr_sta_img = np.expand_dims(pr_sta_img, axis=0)
            #     elif all(s in file for s in ['AvgImage', self.photostim_s.tiff_path.split('/')[-1]]):
            #         ps_sta_img = tf.imread(os.path.join(stam_save_path, file))
            #         ps_sta_img = np.expand_dims(ps_sta_img, axis=0)

            # stack = np.append(stack, pr_sta_img, axis=0)
            # stack = np.append(stack, ps_sta_img, axis=0)

            # target images
            targ_img = getTargetImage()
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

            x_pix = self.imparams.pix_sz_x
            y_pix = self.imparams.pix_sz_y

            save_path = os.path.join(parent_folder, pkl.split('/')[-1][:-4] + '_s2p_masks.tif')

            tf.imwrite(save_path, stack, photometric='minisblack')
            print('\ns2p ROI + photostim targets masks saved in TIFF to: ', save_path)

    def STAProcessing(self, plane):
        """
        Make stimulus triggered average (STA) traces and calculate the STA amplitude
        of response

        Input:
            plane - imaging plane n
        """
        # make stas, [plane x cell x frame]
        stas = np.mean(self.all_trials[plane], axis=2)
        self.stas.append(stas)

        # make sta amplitudes, [plane x cell]
        pre_sta = np.mean(stas[:, self.pre_stim_frames], axis=1)
        post_sta = np.mean(stas[:, self.post_stim_frames], axis=1)
        sta_amplitudes = post_sta - pre_sta
        self.sta_amplitudes.append(sta_amplitudes)

    ####                                                                    // end

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
        targets_trace_full = np.zeros([len(self.Targets.target_coords_all), (end - start) * 2000], dtype='float32')
        counter = 0
        for i in range(start, end):
            tif_path_save2 = self.Suite2p.s2pResultsPath + '/reg_tif/' + reg_tif_list[i]
            with tf.TiffFile(tif_path_save2, multifile=False) as input_tif:
                print('|- reading tiff: %s' % tif_path_save2)
                data = input_tif.asarray()

            targets_trace = np.zeros([len(self.Targets.target_coords_all), data.shape[0]], dtype='float32')
            for coord in range(len(self.Targets.target_coords_all)):
                target_areas = np.array(
                    self.Targets.target_areas)  # TODO update this so that it doesn't include the extra exclusion zone
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
                                        pre_stim=15, post_stim=200, stims: list = None):
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

        if stims is None:
            stim_timings = self.stim_start_frames
        else:
            stim_timings = stims

        if process == 'trace raw':  ## specify which cellsdata to process (i.e. do you want to process whole trace dFF traces?)
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

    # calculate reliability of photostim responsiveness of all of the targeted cells (found in s2p output) TODO need to review this whole section
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

        # choose between .SLMTargets_stims_dff and .SLMTargets_stims_tracedFF for cellsdata to process
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

        # choose between .SLMTargets_stims_dff and .SLMTargets_stims_tracedFF for cellsdata to process
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
        """
        Subtract baseline from dff trials to normalise across cells

        Inputs:
            flu_trial           - [cell x frame] dff trial for all cells
        Outputs:
            baselined_flu_trial - detrended dff trial with zeros replacing stim artifact
        """
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
        """
        Detrend dff trials to account for drift of signal over a trial

        Inputs:
            flu_trial           - [cell x frame] dff trial for all cells
            stim_end            - frame n of the stim end
        Outputs:
            detrended_flu_trial - detrended dff trial with zeros replacing stim artifact
        """
        # set stim artifact period to 0
        flu_trial[:, self.pre_frames:stim_end] = 0

        # detrend and baseline-subtract the flu trial for all cells
        detrended_flu_trial = signal.detrend(self.Suite2p.raw, axis=1)
        baselined_flu_trial = self._baselineFluTrial(detrended_flu_trial)

        return baselined_flu_trial

    #### TEMP // end

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
                # only append trials of the correct length - will catch corrupt/incomplete cellsdata and not include
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
        Creates annotated cellsdata (see anndata library) object based around the Ca2+ matrix of the imaging trial.

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

        print(f"\n\----- CREATING annotated cellsdata object for photostim responses using AnnData:")
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

        # make trial arrays from dff cellsdata shape: [cells x stims x frames]
        if hasattr(self, 'Suite2p'):
            photostimFluArray = self._makePhotostimTrialFluSnippets(plane_flu=self.normalize_dff(self.Suite2p.raw))
            photostimResponseAmplitudes = self.collectPhotostimResponses(photostimFluArray)

            ## create new anndata object for storing measured photostim responses from cellsdata, with other relevant cellsdata
            photostim_responses_adata = self._allCellsPhotostimResponsesAnndata(
                photostimResponseAmplitudes=photostimResponseAmplitudes)

            return photostimFluArray, photostimResponseAmplitudes, photostim_responses_adata
        else:
            NotImplementedError('Photostim processing cannot be performed without Suite2p cellsdata.')

    def statisticalProcessingAllCells(self):
        """Runs statistical processing on photostim response arrays"""
        from packerlabimaging._archive.stats import AllOpticalStats

        self.wilcoxons = AllOpticalStats.runWilcoxonsTest(array1=self.__pre_array, array2=self.__post_array)
        self.sig_units = AllOpticalStats.sigTestAvgResponse(self=self, p_vals=self.wilcoxons, alpha=0.1)

    @property
    def pre_stim_test_slice(self):
        """num of prestim frames used for quantification of photostim responses"""
        return np.s_[self.pre_stim_frames - self.pre_stim_response_frames_window: self.pre_stim_frames]

    @property
    def post_stim_test_slice(self):
        """num of poststim frames used for quantification of photostim responses"""
        stim_end = self.pre_stim_frames + self.stim_duration_frames
        return np.s_[stim_end: stim_end + self.post_stim_response_frames_window]

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
    def whiten_photostim_frame(self, tiff_path, save_as=''):
        """
        TODO docstring

        :param tiff_path:
        :param save_as:
        """
        im_stack = tf.imread(tiff_path, key=range(self.imparams.n_frames))

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
                tiffs_loc = '%s/*Ch3.tif' % self.data_path_dir
                tiff_path = glob.glob(tiffs_loc)[0]
                print('working on loading up %s tiff from: ' % self.metainfo['trialID'], tiff_path)
                im_stack = tf.imread(tiff_path, key=range(self.imparams.n_frames))
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
                    save_path = self.saveDir + 'avg_stim_images'
                    save_path_stim = save_path + '/%s_%s_stim-%s.tif' % (
                        self.metainfo['date'], self.metainfo['trialID'], stim)
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
        """
        run STAmoviemaker for current trial

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
        from packerlabimaging.processing.paq import PaqData

        paqs_loc = f'{BASE_PATH}/2020-12-19/2020-12-19_RL109_013.paq'  # path to the .paq files for the selected trials
        dataPath = alloptical_trial_fixture['dataPath']

        # parses imaging system cellsdata
        imparams = PrairieViewMetadata(pv_xml_dir=os.path.dirname(dataPath), microscope='Bruker 2pPlus')

        # sets the stim start frames
        tmdata = PaqData.paqProcessingAllOptical(paq_path=paqs_loc, frame_channel='frame_clock',
                                                 stim_channel='markpoints2packio???')

        # create the trial
        aotrial = AllOpticalTrial(imparams=imparams, tmdata=tmdata, **alloptical_trial_fixture)
        return aotrial


    idict = alloptical_trial_fixture()
    aotrial = test_AllOpticalClass(idict)
