### DRAFT SCRIPT FOR FUNCS

### TODO consider creating an independent submodule for suite2p_data, and image_processing

"""
This is the main collection of functions for processing and analysis of calcium imaging data in the Packer lab.

The fundamental object of data is one t-series of imaging, along with its associated .paq file, accessory files generated
from the microscope during data collection, and any user generated files associated with the experiment of this t-series.

"""

import re
import glob
import pandas as pd
import time
import datetime
import os
import sys

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)
from packerlabimaging.utils_funcs import SaveDownsampledTiff, subselect_tiff, make_tiff_stack, convert_to_8bit, threshold_detect, \
    s2p_loader, path_finder, points_in_circle_np, moving_average, normalize_dff, paq_read, _check_path_exists

from matplotlib.colors import ColorConverter
import scipy.stats as stats
import statsmodels.api
import statsmodels as sm
from suite2p.run_s2p import run_s2p
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import xml.etree.ElementTree as ET
import tifffile as tf

from packerlabimaging import plotting_utils as plotting
import pickle

###### UTILITIES

# dictionary of terms, phrases, etc. that are used in the processing and analysis of imaging data
terms_dictionary = {
    'dFF': "normalization of datatrace for a given imaging ROI by subtraction and division of a given baseline value",
    'ROI': "a single ROI from the imaging data"
}

def define_term(x):
    try:
        print(f"{x}:    {terms_dictionary[x]}") if type(x) is str else print(
            'ERROR: please provide a string object as the key')
    except KeyError:
        print('input not found in dictionary')

# import .pkl'd objects
def import_obj(pkl_path):
    if not os.path.exists(pkl_path):
        raise FileNotFoundError(f'pkl path NOT found: {pkl_path}' )
    with open(pkl_path, 'rb') as f:
        print(f'\- loading {pkl_path} ... ', end='\r')
        try:
            obj = pickle.load(f)
        except pickle.UnpicklingError:
            raise pickle.UnpicklingError(f"\n** FAILED IMPORT from {pkl_path}\n")
        print(f'|- Loaded {obj.__repr__()}')

    return obj


# option for pre-loading s2p results, and loading Experiment with some trials included in s2p results and some not. -- use s2p_use arg in trialsInformation dict
## TODO need to figure out where to set option for providing input for s2p_batch_size if pre-loading s2p results

## CLASS DEFINITIONS
from dataclasses import dataclass

@dataclass
class Experiment:
    date: str
    comments: str
    dataPath: str
    analysisSavePath: str  # main dir where the experiment object and the trial objects will be saved to
    microscope: str
    expID: str
    trialsInformation: dict # {'trial ID': {'trialType': None, 'tiff_path': None, 'expGroup': None}}
    useSuite2p: bool = False
    # s2pResultsPath: bool = False
    s2pResultsPath: str = None  ## path to the parent directory containing the ops.npy file
    def __post_init__(self):
        print(f'CREATING new Experiment: \n\t{self.__repr__()}')


        ## need to check that the required keys are provided in trialsInformation
        for i in ['trialType', 'tiff_path', 'expGroup']:
            for trial in [*self.trialsInformation]:
                assert i in [*self.trialsInformation[trial]], f"must provide {i} as a key in the trialsInformation dictionary for trial: {trial}"

        # start processing for experiment imaging:

        # start suite2p action:
        if self.useSuite2p or self.s2pResultsPath:
            self.__trialsSuite2p = []
            for trial in [*self.trialsInformation]:
                assert 's2p_use' in [*self.trialsInformation[trial]], 'when trying to utilize suite2p , must provide value for `s2p_use` ' \
                             'in trialsInformation[trial] for each trial to specify if to use trial for this suite2p associated with this experiment'
                self.__trialsSuite2p.append(trial) if self.trialsInformation[trial]['s2p_use'] else None

            if self.s2pResultsPath:  # if s2pResultsPath provided then try to find and pre-load results from provided path, raise error if cannot find results

                # search for suite2p results items in self.suite2pPath, and auto-assign s2pRunComplete --> True if found successfully
                __suite2p_path_files = os.listdir(self.s2pResultsPath)
                self.__s2pResultExists = False
                for filepath in __suite2p_path_files:
                    if 'ops.npy' in filepath:
                        self.__s2pResultExists = True
                        break
                if self.__s2pResultExists:
                    self.Suite2p = Suite2pResultsExperiment(s2pResultsPath=self.s2pResultsPath, trialsSuite2p = self.__trialsSuite2p)
                else:
                    raise ValueError(f"suite2p results could not be found. `suite2pPath` provided was: {self.suite2pPath}")
            elif self.useSuite2p:  # no s2pResultsPath provided, so initialize without pre-loading any results
                self.__s2pResultExists = False
                self.__suite2p_save_path = self.analysisSavePath + '/suite2p/'
                self.Suite2p = Suite2pResultsExperiment(trialsSuite2p = self.__trialsSuite2p)

        # create individual trial objects
        self._runExpTrialsProcessing()


        # save Experiment object:
        # set and/or create analysis save path directory
        if self.analysisSavePath[-4:] == '.pkl':
            self.__pkl_path = analysis_save_path
            self.analysisSavePath = self.analysisSavePath[:[(s.start(), s.end()) for s in re.finditer('/', self.analysisSavePath)][-1][0]]
        else:
            self.analysisSavePath = self.analysisSavePath + '/' if self.analysisSavePath[-1] != '/' else self.analysisSavePath
            self.__pkl_path = f"{self.analysisSavePath}{self.expID}_analysis.pkl"
        os.makedirs(self.analysisSavePath, exist_ok=True)
        self.save_pkl(pkl_path=self.pkl_path)

        # Attributes:
        self.n_frames = None  # total number of imaging frames in the


        print(f"\n\n\nNEW Experiment object created: ")
        print(self)



    def __repr__(self):
        return f"Experiment object (date: {self.date}, expID: {self.expID})"

    def __str__(self):
        lastsaved = time.ctime(os.path.getmtime(self.pkl_path))
        __return_information = f"Experiment object (last saved: {lastsaved}), date: {self.date}, expID: {self.expID}, microscope: {self.microscope}"
        __return_information = __return_information + f"\npkl path: {self.pkl_path}"
        __return_information = __return_information + f"\ntrials in Experiment object:"
        for trial in self.trialsInformation:
            # print(f"\t{trial}: {self.trialsInformation[trial]['trialType']} {self.trialsInformation[trial]['expGroup']}")
            __return_information = __return_information + f"\n\t{trial}: {self.trialsInformation[trial]['trialType']}, {self.trialsInformation[trial]['expGroup']}"
        return __return_information

    def _runExpTrialsProcessing(self):

        total_frames_stitched = 0  # used in calculating # of frames from a single trial in the overall suite2p run
        for trial in self.trialsInformation:
            print(f"\n\n\- PROCESSING trial: {trial}, expID: ({self.expID})")
            _metainfo = {
                'animal prep.': self.expID,
                'trial': trial,
                'date': self.date,
                't series id': f"{self.expID} {trial}",
                'trialsInformation': self.trialsInformation[trial]
            }

            # initialize TwoPhotonImagingTrial
            if _metainfo['trialsInformation']['trialType'] == 'TwoPhotonImagingTrial':
                if 'tiff_path' not in [*_metainfo['trialsInformation']]:
                    raise ValueError('TwoPhotonImagingTrial experiment trial requires `tiff_path` field defined in .trialsInformation dictionary for each trial')

                trial_obj = TwoPhotonImagingTrial(metainfo=_metainfo, analysis_save_path=self.analysisSavePath,
                                                  microscope=self.microscope)

                # update self.trialsInformation using the information from new trial_obj
                self.trialsInformation[trial]['analysis Object Information'] = {'series ID': trial_obj.t_series_name,
                                                                                'repr': trial_obj.__repr__(),
                                                                                'pkl path': trial_obj.pkl_path}

            # initialize AllOpticalTrial
            elif _metainfo['trialsInformation']['trialType'] == 'AllOpticalTrial':
                if 'tiff_path' not in [*_metainfo['trialsInformation']] or 'paq_path' not in [*_metainfo['trialsInformation']] \
                        or 'naparm_path' not in [*_metainfo['trialsInformation']]:
                    raise ValueError(f'AllOpticalTrial experiment trial requires `tiff_path`, `paq_path` and `naparm_path` fields defined in .trialsInformation dictionary for each alloptical trial. '
                                     f'\n{self.trialsInformation[trial]}')
                trial_obj = AllOpticalTrial(microscope=self.microscope, naparm_path=_metainfo['trialsInformation']['naparm_path'],
                                            analysis_save_path=self.analysisSavePath, metainfo=_metainfo,
                                            pre_stim= 1.0, post_stim= 3.0, pre_stim_response_window = 0.500, post_stim_response_window = 0.500)
                # update self.trialsInformation using the information from new trial_obj
                self.trialsInformation[trial]['analysis Object Information'] = {'series ID': trial_obj.t_series_name,
                                                                                'repr': trial_obj.__repr__(),
                                                                                'pkl path': trial_obj.pkl_path}
            else:
                raise ValueError(f"unsupported trial type for trial: {trial}. All trials must be of trialType: TwoPhotonImagingTrial or AllOpticalTrial")

            # initialize suite2p for trial objects
            if trial in self.__trialsSuite2p:
                print(f"\n\----- ADDING Suite2p class to Trial object ... ")
                trial_obj.Suite2p = Suite2pResultsTrial(suite2p_experiment_obj=self.Suite2p,
                                                        trial_frames=(total_frames_stitched, total_frames_stitched + trial_obj.n_frames))  # use trial obj's current trial frames
                total_frames_stitched += trial_obj.n_frames

    @property
    def tiff_path_dir(self):
        # __first_trial_in_experiment = list(self.trialsInformation.keys())[0]
        __first_trial_in_experiment = [*self.trialsInformation][0]
        __tiff_path_first_trial = self.trialsInformation[__first_trial_in_experiment]['tiff_path']
        return __tiff_path_first_trial[:[(s.start(), s.end()) for s in re.finditer('/', __tiff_path_first_trial)][-1][0]]  # this is the directory where the Bruker xml files associated with the 2p imaging TIFF are located


    @property
    def pkl_path(self):
        "path in Analysis folder to save pkl object"
        # return f"{self.analysis_save_dir}{self.metainfo['date']}_{self.metainfo['trial']}.pkl"
        return self.__pkl_path

    @pkl_path.setter
    def pkl_path(self, path: str):
        self.__pkl_path = path

    def save_pkl(self, pkl_path: str = None):
        if pkl_path:
            print(f'\nsaving new pkl object at: {pkl_path}')
            self.pkl_path = pkl_path

        with open(self.pkl_path, 'wb') as f:
            pickle.dump(self, f)
        print("\t -- Experiment analysis object saved to %s -- " % self.pkl_path)

    def save(self):
        self.save_pkl()

    ## suite2p methods
    def s2pRun(self, user_batch_size=2000, trialsSuite2P: list = None):
        """run suite2p on the experiment object, using trials specified in current experiment object, using the attributes
        determined directly from the experiment object."""

        ops = self.Suite2p.ops
        db = self.Suite2p.db

        self.Suite2p.trials = trialsSuite2P if trialsSuite2P else self.Suite2p.trials
        self.__trialsSuite2p = trialsSuite2P if trialsSuite2P else self.__trialsSuite2p

        tiffs_to_use_s2p = []
        for trial in self.Suite2p.trials:
            tiffs_to_use_s2p.append(self.trialsInformation[trial]['tiff_path'])

        sampling_rate = self.fps / self.n_planes
        diameter_x = 13 / self.pix_sz_x
        diameter_y = 13 / self.pix_sz_y
        diameter = int(diameter_x), int(diameter_y)
        self.Suite2p.user_batch_size = user_batch_size
        batch_size = self.Suite2p.user_batch_size * (262144 / (self.frame_x * self.frame_y))  # larger frames will be more RAM intensive, scale user batch size based on num pixels in 512x512 images

        if not db:
            db = {
                'fs': float(sampling_rate),
                'diameter': diameter,
                'batch_size': int(batch_size),
                'nimg_init': int(batch_size),
                'nplanes': self.n_planes
            }

        # specify tiff list to suite2p and data path
        db['tiff_list'] = tiffs_to_s2p
        db['data_path'] = self.dataPath  ## this is where the bad_frames.npy file will be stored for suite2p to use.
        db['save_folder'] = self.suite2pPath

        print(db)

        opsEnd = run_s2p(ops=ops, db=db)

        ## TODO update Experiment attr's and Trial attr's to reflect completion of the suite2p RUN
        self.__s2pResultExists = True
        self.s2pResultsPath = self.__suite2p_save_path + '/plane0/'

class Suite2pResultsExperiment:
    """used to run and further process suite2p processed data, and analysis associated with suite2p processed data."""

    def __init__(self, trialsSuite2p: list, s2pResultsPath: str = None, subtract_neuropil: bool = True):
        print(f"\- ADDING Suite2p class to Experiment object ... ")

        # set trials to run together in suite2p for Experiment
        self.trials = trialsSuite2p
        self.subtract_neuropil = subtract_neuropil
        assert len(self.trials) > 0, "no trials found to run suite2p, option available to provide list of trial IDs in `trialsSuite2P`"

        if s2pResultsPath is None:
            ## initialize needed variables and attr's for future calling of s2pRun
            self.__s2pResultsPath = None
            self.ops = {}
            self.db = {}
        else:
            self.__s2pResultsPath = s2pResultsPath
            try:
                self.s2pProcessing(self.__s2pResultsPath, subtract_neuropil, subset_frames, save)
            except:
                raise ValueError(f'suite2p processed data could not be loaded from the provided `s2pResultsPath`: {s2pResultsPath}')

        # Attributes
        self.n_frames = None  # total number of imaging frames in the Suite2p run

    def s2pProcessing(self, s2p_path: str, subtract_neuropil: bool = True, save: bool = True):
        """processing of suite2p data from the current t-series
        :param s2p_path: path to the directory containing suite2p outputs
        :param subtract_neuropil: choose to subtract neuropil or not when loading s2p traces
        :param save: choose to save data object or not
        """

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
                # s2pResultsPath = os.path.join(self.tiff_path, 'suite2p', 'plane' + str(plane))
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


    def stitch_reg_tiffs(self, first_frame: int, last_frame: int, reg_tif_folder: str = None, force_crop: bool = False,
                         s2p_run_batch: int = 2000):  # TODO refactor as method for trial object
        """
        Stitches together registered tiffs outputs from suite2p from the provided imaging frame start and end values.

        :param first_frame: first frame from the overall s2p run to start stitching from
        :param last_frame: last frame from the overall s2p run to end stitching at
        :param force_crop:
        :param force_stack:
        :param s2p_run_batch: batch size for suite2p run (defaults to 2000 because that is usually the batch size while running suite2p processing)
        """

        if reg_tif_folder is None:
            if self.suite2p_path:
                reg_tif_folder = self.suite2p_path + '/reg_tif/'
                print(f"\- trying to load registerred tiffs from: {reg_tif_folder}")
        else:
            raise Exception(f"Must provide reg_tif_folder path for loading registered tiffs")
        if not os.path.exists(reg_tif_folder):
            raise Exception(f"no registered tiffs found at path: {reg_tif_folder}")

        frame_range = [first_frame, last_frame]

        start = first_frame // s2p_run_batch
        end = last_frame // s2p_run_batch + 1

        # set tiff paths to save registered tiffs:
        tif_path_save = self.analysis_save_dir + 'reg_tiff_%s.tif' % self.metainfo['trial']
        tif_path_save2 = self.analysis_save_dir + 'reg_tiff_%s_r.tif' % self.metainfo['trial']
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        sorted_paths = [reg_tif_folder + tif for tif in reg_tif_list][start:end + 1]

        print(tif_path_save)
        print(sorted_paths)

        if os.path.exists(tif_path_save):
            make_tiff_stack(sorted_paths, save_as=tif_path_save)

        if not os.path.exists(tif_path_save2) or force_crop:
            with tf.TiffWriter(tif_path_save2, bigtiff=True) as tif:
                with tf.TiffFile(tif_path_save, multifile=False) as input_tif:
                    print('cropping registered tiff')
                    data = input_tif.asarray()
                    print('shape of stitched tiff: ', data.shape)
                reg_tif_crop = data[frame_range[0] - start * s2p_run_batch: frame_range[1] - (
                            frame_range[0] - start * s2p_run_batch)]
                print('saving cropped tiff ', reg_tif_crop.shape)
                tif.save(reg_tif_crop)

    def s2pMeanImage(self, plot: bool = True):
        """
        Return array of the s2p mean image.
        :param s2p_path: (optional) path to location of s2p data
        :param plot: (optional) option to plot the s2p mean image
        :return:
        """

        print(f'Plotting s2p mean image .. {self.__s2pResultsPath}')

        os.chdir(self.__s2pResultsPath)

        ops = np.load('ops.npy', allow_pickle=True).item()

        mean_img = ops['meanImg']

        mean_img = np.array(mean_img, dtype='uint16')

        if plot:
            plt.imshow(mean_img, cmap='gray')
            plt.suptitle('s2p mean image')
            plt.show()

        return mean_img

    ### TODO add methods for processing suite2p ROIs

class Suite2pResultsTrial:
    """used to collect suite2p processed data for one trial - out of overall experiment."""

    def __init__(self, suite2p_experiment_obj: Suite2pResultsExperiment, trial_frames: tuple = None):
        self.trial_frames = trial_frames  # tuple of first and last frame (out of the overall suite2p run) corresponding to the present trial

        ## TODO adding attributes to collect Suite2p results for this specific trial

        self.__overall_suite2p = suite2p_experiment_obj

        # __s2pResultsPath = suite2p_experiment_obj.__s2pResultsPath if hasattr(suite2p_experiment_obj, '__s2pResultsPath') else None
        # Suite2pResultsExperiment.__init__(self, trialsSuite2p = suite2p_experiment_obj.trials, s2pResultsPath=__s2pResultsPath,
        #                                   subtract_neuropil=suite2p_experiment_obj.subtract_neuropil)

    def stitch_reg_tiffs(self, first_frame: int, last_frame: int, reg_tif_folder: str = None, force_crop: bool = False,
                         s2p_run_batch: int = 2000):

        """
        Stitches together registered tiffs outputs from suite2p from the provided imaging frame start and end values.

        :param first_frame: first frame from the overall s2p run to start stitching from
        :param last_frame: last frame from the overall s2p run to end stitching at
        :param force_crop:
        :param force_stack:
        :param s2p_run_batch: batch size for suite2p run (defaults to 2000 because that is usually the batch size while running suite2p processing)
        """

        assert hasattr(self, 'user_batch_size'), 'No user_batch_size set for Suite2pResultsTrial, please create new attr for .Suite2p.user_batch_size, before continuing'

        if reg_tif_folder is None:
            if self.s2pResultsPath:
                reg_tif_folder = self.s2pResultsPath + '/reg_tif/'
                print(f"\- trying to load registerred tiffs from: {reg_tif_folder}")
        else:
            raise Exception(f"Must provide reg_tif_folder path for loading registered tiffs")
        if not os.path.exists(reg_tif_folder):
            raise Exception(f"no registered tiffs found at path: {reg_tif_folder}")

        first_frame = self.trial_frames[0]
        last_frame = self.trial_frames[-1]
        frame_range = [first_frame, last_frame]

        start = first_frame // s2p_run_batch
        end = last_frame // s2p_run_batch + 1

        # set tiff paths to save registered tiffs:
        tif_path_save = self.analysis_save_dir + 'reg_tiff_%s.tif' % self.metainfo['trial']
        tif_path_save2 = self.analysis_save_dir + 'reg_tiff_%s_r.tif' % self.metainfo['trial']
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        sorted_paths = [reg_tif_folder + tif for tif in reg_tif_list][start:end + 1]

        print(tif_path_save)
        print(sorted_paths)

        if os.path.exists(tif_path_save):
            make_tiff_stack(sorted_paths, save_as=tif_path_save)

        if not os.path.exists(tif_path_save2) or force_crop:
            with tf.TiffWriter(tif_path_save2, bigtiff=True) as tif:
                with tf.TiffFile(tif_path_save, multifile=False) as input_tif:
                    print('cropping registered tiff')
                    data = input_tif.asarray()
                    print('shape of stitched tiff: ', data.shape)
                reg_tif_crop = data[frame_range[0] - start * s2p_run_batch: frame_range[1] - (
                            frame_range[0] - start * s2p_run_batch)]
                print('saving cropped tiff ', reg_tif_crop.shape)
                tif.save(reg_tif_crop)

class TwoPhotonImagingTrial:
    """Just two photon imaging related functions - currently set up for data collected from Bruker microscopes and
    suite2p processed Ca2+ imaging data """

    def __init__(self, metainfo: dict, analysis_save_path: str, microscope: str = 'Bruker'):
        """
        TODO update function docstring for approp args
        :param tiff_path: path to t-series .tiff
        :param analysis_save_path: path of where to save the experiment analysis object
        :param microscope: name of microscope used to record imaging (options: "Bruker" (default), "Scientifica", "other")
        :param metainfo: dictionary containing any metainfo information field needed for this experiment. At minimum it needs to include prep #, t-series # and date of data collection.
        :param paq_path: (optional) path to .paq file associated with current t-series
        :param suite2p_path: (optional) path to suite2p outputs folder associated with this t-series (plane0 file? or ops file? not sure yet)
        :param make_downsampled_tiff: flag to run generation and saving of downsampled tiff of t-series (saves to the analysis save location)
        """

        print(f'\----- CREATING TwoPhotonImagingTrial for trial: {metainfo["trial"]},  {metainfo["t series id"]}')

        # Initialize Attributes:
        self.n_frames = None  # number of imaging frames in the current trial


        if 'date' in metainfo.keys() and 'trial' in metainfo.keys() and 'animal prep.' in metainfo.keys() and 't series id' in metainfo.keys(): self.__metainfo = metainfo
        else: raise ValueError("dev error: __metainfo argument must contain the minimum fields: 'date', 'trial', 'animal prep.' and 't series id'")
        if os.path.exists(metainfo['trialsInformation']['tiff_path']): self.tiff_path = metainfo['trialsInformation']['tiff_path']
        else: raise FileNotFoundError(f"tiff_path does not exist: {metainfo['trialsInformation']['tiff_path']}")
        if 'paq_path' in [*metainfo['trialsInformation']]:
            if os.path.exists(metainfo['trialsInformation']['paq_path']): self.paq_path = metainfo['trialsInformation']['paq_path']
            else: raise FileNotFoundError(f"paq_path does not exist: {metainfo['trialsInformation']['paq_path']}")

        # set and create analysis save path directory
        self.analysis_save_dir = analysis_save_path
        self.__pkl_path = f"{self.analysis_save_dir}{metainfo['date']}_{metainfo['trial']}.pkl"

        self.save_pkl(pkl_path=self.pkl_path)  # save experiment object to pkl_path

        self._parsePVMetadata() if 'Bruker' in microscope else Warning(f'retrieving data-collection metainformation from '
                                                                       f'{microscope} microscope has not been implemented yet')

        # run paq processing if paq_path is provided for trial
        TwoPhotonImagingTrial.paqProcessing(self, paq_path=self.paq_path, plot=False) if hasattr(self, 'paq_path') else None

        self.save()





    def __repr__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
        else:
            lastmod = "(unsaved pkl object)"
        return repr(f"({self.t_series_name}) TwoPhotonImagingTrial experimental data object, last saved: {lastmod}")

    @property
    def fig_save_path(self):
        today_date = datetime.today().strftime('%Y-%m-%d')
        return self.analysis_save_dir + f'Results_fig/{today_date}/'

    @fig_save_path.setter
    def fig_save_path(self, value: str):
        "set new default fig save path for data object"
        self.fig_save_path = value

    @property
    def date(self):
        "date of the experiment datacollection"
        return self.__metainfo['date']

    @property
    def prep(self):
        return self.__metainfo['animal prep.']

    @property
    def trial(self):
        return self.__metainfo['trial']

    @property
    def t_series_name(self):
        if 't series id' in self.__metainfo.keys():
            return f"{self.__metainfo['t series id']}"
        elif "animal prep." in self.__metainfo.keys() and "trial" in self.__metainfo.keys():
            return f'{self.__metainfo["animal prep."]} {self.__metainfo["trial"]}'
        else:
            raise ValueError('no information found to retrieve t series id')

    @property
    def tiff_path_dir(self):
        return self.tiff_path[:[(s.start(), s.end()) for s in re.finditer('/', self.tiff_path)][-1][0]]  # this is the directory where the Bruker xml files associated with the 2p imaging TIFF are located

    @property
    def pkl_path(self):
        "path in Analysis folder to save pkl object"
        # return f"{self.analysis_save_dir}{self.__metainfo['date']}_{self.__metainfo['trial']}.pkl"
        return self.__pkl_path

    @pkl_path.setter
    def pkl_path(self, path: str):
        self.__pkl_path = path

    def _getPVStateShard(self, root, key):
        '''
        Find the value, description and indices of a particular parameter from an xml file

        Inputs:
            path        - path to xml file
            key         - string corresponding to key in xml tree
        Outputs:
            value       - value of the key
            description - unused
            index       - index that the key was found at
        '''
        value = []
        description = []
        index = []

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
        '''
        Parse all of the relevant acquisition metadata from the PrairieView xml file for this recording
        '''

        print('\n\----- Parsing PV Metadata for Bruker microscope...')


        tiff_path = self.tiff_path_dir  # starting path
        xml_path = []  # searching for xml path

        try:  # look for xml file in path, or two paths up (achieved by decreasing count in while loop)
            count = 2
            while count != 0 and not xml_path:
                count -= 1
                for file in os.listdir(tiff_path):
                    if file.endswith('.xml'):
                        xml_path = os.path.join(tiff_path, file)
                tiff_path = os.path.dirname(tiff_path)  # re-assign tiff_path as next folder up

        except:
            raise Exception('ERROR: Could not find xml for this acquisition, check it exists')

        xml_tree = ET.parse(xml_path)  # parse xml from a path
        root = xml_tree.getroot()  # make xml tree structure

        sequence = root.find('Sequence')
        acq_type = sequence.get('type')

        if 'ZSeries' in acq_type:
            n_planes = len(sequence.findall('Frame'))
        else:
            n_planes = 1

        frame_branch = root.findall('Sequence/Frame')[-1]
        #         frame_period = float(self._getPVStateShard(root,'framePeriod')[0])
        frame_period = float(self._getPVStateShard(frame_branch, 'framePeriod')[0])
        fps = 1 / frame_period

        frame_x = int(self._getPVStateShard(root, 'pixelsPerLine')[0])
        frame_y = int(self._getPVStateShard(root, 'linesPerFrame')[0])
        zoom = float(self._getPVStateShard(root, 'opticalZoom')[0])

        scan_volts, _, index = self._getPVStateShard(root, 'currentScanCenter')
        for scan_volts, index in zip(scan_volts, index):
            if index == 'XAxis':
                scan_x = float(scan_volts)
            if index == 'YAxis':
                scan_y = float(scan_volts)

        pixel_size, _, index = self._getPVStateShard(root, 'micronsPerPixel')
        for pixel_size, index in zip(pixel_size, index):
            if index == 'XAxis':
                pix_sz_x = float(pixel_size)
            if index == 'YAxis':
                pix_sz_y = float(pixel_size)

        if n_planes == 1:
            n_frames = root.findall('Sequence/Frame')[-1].get('index')  # use suite2p output instead later
        else:
            n_frames = root.findall('Sequence')[-1].get('cycle')

        extra_params = root.find('Sequence/Frame/ExtraParameters')
        last_good_frame = extra_params.get('lastGoodFrame')

        self.fps = fps / n_planes
        self.frame_x = frame_x
        self.frame_y = frame_y
        self.n_planes = n_planes
        self.pix_sz_x = pix_sz_x
        self.pix_sz_y = pix_sz_y
        self.scan_x = scan_x
        self.scan_y = scan_y
        self.zoom = zoom
        self.n_frames = int(n_frames)
        self.last_good_frame = last_good_frame

        print('\tn planes:', n_planes,
              '\n\tn frames:', int(n_frames),
              '\n\tfps:', fps,
              '\n\tframe size (px):', frame_x, 'x', frame_y,
              '\n\tzoom:', zoom,
              '\n\tpixel size (um):', pix_sz_x, pix_sz_y,
              '\n\tscan centre (V):', scan_x, scan_y
              )

    def paqProcessing(self, paq_path: str = None, plot: bool = False):
        """
        Loads .paq file and saves data from individual channels.

        :param paq_path: (optional) path to the .paq file for this data object
        """

        print('\n\----- Processing paq file ...')

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

        ## TODO print the paq channels that were loaded. and some useful metadata about the paq channels.

        # find frame times
        clock_idx = paq['chan_names'].index('frame_clock')
        clock_voltage = paq['data'][clock_idx, :]

        frame_clock = threshold_detect(clock_voltage, 1)
        self.frame_clock = frame_clock

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

    def meanRawFluTrace(self, plot: bool = False, save: bool = True):
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

        plotting.plotMeanRawFluTrace(expobj=self, stim_span_color=None, x_axis='frames', figsize=[20, 3],
                                     title='Mean raw Flu trace -') if plot else None

        return im_stack

    def makeDownsampledTiff(self):
        stack = self.meanRawFluTrace(save_pkl=True)
        SaveDownsampledTiff(stack=stack,
                            save_as=f"{self.analysis_save_dir}/{metainfo['date']}_{metainfo['trial']}_downsampled.tif")

    def plotSingleImageFrame(self, frame_num: int = 0, title: str = None):
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

    def save_pkl(self, pkl_path: str = None):
        ## commented out after setting pkl_path as a @property
        # if pkl_path is None:
        #     if hasattr(self, 'pkl_path'):
        #         pkl_path = self.pkl_path
        #     else:
        #         raise ValueError(
        #             'pkl path for saving was not found in data object attributes, please provide pkl_path to save to')
        # else:
        #     self.pkl_path = pkl_path
        if pkl_path:
            print(f'saving new trial object to: {pkl_path}')
            self.pkl_path = pkl_path

        with open(self.pkl_path, 'wb') as f:
            pickle.dump(self, f)
        print("\n\t -- data object saved to %s -- " % self.pkl_path)

    def save(self):
        self.save_pkl()

class AllOpticalTrial(TwoPhotonImagingTrial):
    """This class provides methods for All Optical experiments"""

    def __init__(self, metainfo: dict, naparm_path: str, analysis_save_path: str, microscope: str, analysisOptions: dict = {},
                 pre_stim: float = 1.0, post_stim: float = 3.0, pre_stim_response_window: float = 0.500,
                 post_stim_response_window: float = 0.500):

        # TODO update providing arg values for all optical experiment analysis hyperparameters

        """
        :param tiff_path: path to t-series .tiff
        :param paq_path: path to .paq file associated with current t-series
        :param naparm_path: path to folder containing photostimulation setup built by NAPARM
        :param analysis_save_path: path of where to save the experiment analysis object
        :param metainfo: dictionary containing any metainfo information field needed for this experiment. At minimum it needs to include prep #, t-series # and date of data collection.
        :param microscope: name of microscope used to record imaging (options: "Bruker" (default), "Scientifica", "other")
        :param suite2p_path: (optional) path to suite2p outputs folder associated with this t-series (plane0 file? or ops file? not sure yet)
        :param make_downsampled_tiff: flag to run generation and saving of downsampled tiff of t-series (saves to the analysis save location)
        :param pre_stim:
        :param post_stim:
        :param pre_stim_response_window:
        :param post_stim_response_window:
        :kwargs (optional):
        """

        print(f'\----- CREATING AllOpticalTrial data object for {metainfo["t series id"]}')

        # setting initial attr's
        self.stim_channel = kwargs['stim_channel'] if 'stim_channel' in [*analysisOptions] else 'markpoints2packio'  # stim_channel = channel on paq file to read for determining stims


        #### initializing data processing, data analysis and/or results associated attr's
        self.n_trials = None  # number of photostimulation trials TODO change to assigning from array: cells x Flu frames x # of photostim trials


        ## PHOTOSTIM SLM TARGETS
        self.responses_SLMtargets = []  # dF/prestimF responses for all SLM targets for each photostim trial
        self.responses_SLMtargets_tracedFF = []  # poststim dFF - prestim dFF responses for all SLM targets for each photostim trial

        # TODO add attr's related to numpy array's and pandas dataframes for photostim trials - SLM targets

        ## NON PHOTOSTIM SLM TARGETS

        # TODO add attr's related to numpy array's and pandas dataframes for photostim trials - non SLM suite2p ROIs


        # set photostim analysis time windows
        self.__pre_stim = pre_stim
        self.__post_stim = post_stim
        self.__pre_stim_response_window = pre_stim_response_window
        self.__post_stim_response_window = post_stim_response_window
        # self._set_photostim_windows(pre_stim, post_stim, pre_stim_response_window, post_stim_response_window)  -- delete line

        if os.path.exists(naparm_path): self.__naparm_path = naparm_path
        else: raise FileNotFoundError(f"path not found, naparm_path: {naparm_path}")


        # initialize object as TwoPhotonImagingTrial
        TwoPhotonImagingTrial.__init__(self, metainfo=metainfo, analysis_save_path=analysis_save_path,
                                       microscope=microscope)


        # continue with photostimulation experiment processing
        self._stimProcessing(stim_channel=self.stim_channel)
        self._findTargetsAreas()
        self.photostim_frames = ['not-yet-processed']
        self._find_photostim_add_bad_framesnpy()



        ##
        self.save()

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
        else:
            lastmod = "(unsaved pkl object)"
        if not hasattr(self, 'metainfo'):
            information = f"uninitialized"
        else:
            information = self.t_series_name

        return repr(f"({information}) TwoPhotonImagingTrial.alloptical experimental data object, last saved: {lastmod}")

    @property
    def naparm_path(self):
        """setting location of naparm files that specify the photostimulation experiment setup"""
        if self.__naparm_path[-1] == '/':
            return self.__naparm_path
        else:
            return self.__naparm_path + '/'

    @property
    def pre_stim(self):
        """setting time window for collecting Flu trace before each photostimulation trial"""
        return self.__pre_stim

    @pre_stim.setter
    def pre_stim(self, value):
        self.__pre_stim = value

    @property
    def post_stim(self):
        """setting time window for collecting Flu trace after each photostimulation trial"""
        return self.__post_stim

    @post_stim.setter
    def post_stim(self, value):
        self.__post_stim = value

    @property
    def pre_stim_response_window(self):
        """setting time window for measuring Flu trace before each photostimulation trial"""
        return self.__pre_stim_response_window

    @pre_stim_response_window.setter
    def pre_stim_response_window(self, value):
        self.__pre_stim_response_window = value

    @property
    def post_stim_response_window(self):
        """setting time window for measuring Flu trace after each photostimulation trial"""
        return self.__post_stim_response_window

    @post_stim_response_window.setter
    def post_stim_response_window(self, value):
        self.__post_stim_response_window = value

    @property
    def pre_stim_response_frames_window(self):
        return int(self.fps * self.pre_stim_response_window)  # length of the pre stim response test window (in frames)

    @property
    def post_stim_response_frames_window(self):
        return int(self.fps * self.post_stim_response_frames_window)  # length of the post stim response test window (in frames)

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
        self.frame_clock = frame_clock
        if plot:
            plt.figure(figsize=(10, 5))
            plt.plot(clock_voltage)
            plt.plot(frame_clock, np.ones(len(frame_clock)), '.')
            plt.suptitle('frame clock from paq, with detected frame clock instances as scatter')
            sns.despine()
            plt.show()

        # find start and stop frame_clock times -- there might be multiple 2p imaging starts/stops in the paq trial (hence multiple frame start and end times)
        self.frame_start_times = [self.frame_clock[0]]  # initialize ls
        self.frame_end_times = []
        i = len(self.frame_start_times)
        for idx in range(1, len(self.frame_clock) - 1):
            if (self.frame_clock[idx + 1] - self.frame_clock[idx]) > 2e3:
                i += 1
                self.frame_end_times.append(self.frame_clock[idx])
                self.frame_start_times.append(self.frame_clock[idx + 1])
        self.frame_end_times.append(self.frame_clock[-1])

        # for frame in self.frame_clock[1:]:
        #     if (frame - self.frame_start_times[i - 1]) > 2e3:
        #         i += 1
        #         self.frame_start_times.append(frame)
        #         self.frame_end_times.append(self.frame_clock[np.where(self.frame_clock == frame)[0] - 1][0])
        # self.frame_end_times.append(self.frame_clock[-1])

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

        # find stim times
        stim_idx = paq['chan_names'].index(self.stim_channel)
        stim_volts = paq['data'][stim_idx, :]
        stim_times = threshold_detect(stim_volts, 1)
        # self.stim_times = stim_times
        self.stim_start_times = stim_times
        print('# of stims found on %s: %s' % (self.stim_channel, len(self.stim_start_times)))

        # correct this based on txt file
        duration_ms = self.stim_dur
        frame_rate = self.fps / self.n_planes
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

        self.stim_start_frames = []

        for plane in range(self.n_planes):

            stim_start_frames = []

            for stim in stim_times:
                # the index of the frame immediately preceeding stim
                stim_start_frame = next(
                    i - 1 for i, sample in enumerate(frame_clock[plane::self.n_planes]) if sample - stim >= 0)
                stim_start_frames.append(stim_start_frame)

            self.stim_start_frames = stim_start_frames

    ### KEY FUNCTIONS TO REVIEW WITH ROB FOR ALLOPTICAL WORKFLOW

    def _photostimProcessing(self):  ## TODO need to figure out how to handle different types of photostimulation experiment setups

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
        print('Single stim. Duration (ms): ', self.single_stim_dur)
        print('Total stim. Duration per trial (ms): ', self.stim_dur)

        self.paqProcessing(plot=False)

    def _stimProcessing(self, stim_channel):
        self.stim_channel = stim_channel
        self._photostimProcessing()

    def _find_photostim_add_bad_framesnpy(self):
        """finds all photostim frames and saves them into the bad_frames attribute for the exp object and saves bad_frames.npy"""
        print('\n\-----Finding photostimulation frames in imaging frames ...')
        print('# of photostim frames calculated per stim. trial: ', self.stim_duration_frames + 1)

        photostim_frames = []
        for j in self.stim_start_frames:
            for i in range(
                    self.stim_duration_frames + 1):  # usually need to remove 1 more frame than the stim duration, as the stim isn't perfectly aligned with the start of the imaging frame
                photostim_frames.append(j + i)

        self.photostim_frames = photostim_frames
        # print(photostim_frames)
        print('\t|- Original # of frames:', self.n_frames, 'frames')
        print('\t|- # of Photostim frames:', len(photostim_frames), 'frames')
        print('\t|- Minus photostim. frames total:', self.n_frames - len(photostim_frames), 'frames')

        if len(self.photostim_frames) > 0:
            print(f'***Saving a total of {len(self.photostim_frames)} photostim frames to bad_frames.npy at: {self.tiff_path_dir}/bad_frames.npy')
            np.save(f'{self.tiff_path_dir}/bad_frames.npy', self.photostim_frames)  # save to npy file and remember to move npy file to tiff folder before running with suite2p

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

    def _findTargetsAreas(self):

        '''
        Finds cells that have been targeted for optogenetic photostimulation using Naparm in all-optical type experiments.
        output: coordinates of targets, and circular areas of targets
        Note this is not done by target groups however. So all of the targets are just in one big ls.
        '''

        print('\n\t\-----Loading up target coordinates...')

        self.n_targets = []
        self.target_coords = []
        self.target_areas = []

        # load naparm targets file for this experiment
        naparm_path = os.path.join(self.naparm_path, 'Targets')

        listdir = os.listdir(naparm_path)

        scale_factor = self.frame_x / 512

        ## All SLM targets
        for path in listdir:
            if 'AllFOVTargets' in path:
                target_file = path
        target_image = tf.imread(os.path.join(naparm_path, target_file))

        # # SLM group#1: FOVTargets_001
        # for path in listdir:
        #     if 'FOVTargets_001' in path:
        #         target_file_1 = path
        # target_image_slm_1 = tf.imread(os.path.join(naparm_path, target_file_1))
        #
        # # SLM group#2: FOVTargets_002
        # for path in listdir:
        #     if 'FOVTargets_002' in path:
        #         target_file_2 = path
        # target_image_slm_2 = tf.imread(os.path.join(naparm_path, target_file_2))

        ## idk if there a reason for this...but just keeping it here since Rob had used it for this analysis
        # n = np.array([[0, 0], [0, 1]])
        # target_image_scaled = np.kron(target_image, n)

        # target_image_scaled_1 = target_image_slm_1;
        # del target_image_slm_1
        # target_image_scaled_2 = target_image_slm_2;
        # del target_image_slm_2

        # if frame_x < 1024 or frame_y < 1024:
        #     pass
        # #             # bounding box coords
        # #             x1 = 511 - frame_x / 2
        # #             x2 = 511 + frame_x / 2
        # #             y1 = 511 - frame_y / 2
        # #             y2 = 511 + frame_y / 2
        #
        # #             # calc imaging galvo offset between BAFOV and t-series
        # #             zoom = self.zoom
        # #             scan_x = self.scan_x  # scan centre in V
        # #             scan_y = self.scan_y
        #
        # #             ScanAmp_X = 2.62 * 2
        # #             ScanAmp_Y = 2.84 * 2
        #
        # #             ScanAmp_V_FOV_X = ScanAmp_X / zoom
        # #             ScanAmp_V_FOV_Y = ScanAmp_Y / zoom
        #
        # #             scan_pix_y = ScanAmp_V_FOV_Y / 1024
        # #             scan_pix_x = ScanAmp_V_FOV_X / 1024
        #
        # #             offset_x = scan_x / scan_pix_x  # offset between image centres in pixels
        # #             offset_y = scan_y / scan_pix_y
        #
        # #             # offset the bounding box
        # #             x1, x2, y1, y2 = round(x1 + offset_x), round(x2 + offset_x), round(y1 - offset_y), round(y2 - offset_y)
        #
        # #             # crop the target image using the offset bounding box to get the targets in t-series imaging space
        # #             target_image_scaled = target_image_scaled[y1:y2, x1:x2]
        # #             tf.imwrite(os.path.join(naparm_path, 'target_image_scaled.tif'), target_image_scaled)
        # else:
        #     #             image_8bit = convert_to_8bit(target_image_scaled, np.unit8)
        #     #             tf.imwrite(os.path.join(naparm_path, 'target_image_scaled.tif'), image_8bit)
        #     tf.imwrite(os.path.join(naparm_path, 'target_image_scaled.tif'), target_image_scaled)

        targets = np.where(target_image > 0)
        # targets_1 = np.where(target_image_scaled_1 > 0)
        # targets_2 = np.where(target_image_scaled_2 > 0)

        targetCoordinates = list(zip(targets[1] * scale_factor, targets[0] * scale_factor))
        print('\tNumber of targets:', len(targetCoordinates))

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
        radius_px = int(np.ceil(((self.spiral_size / 2) + 0) / self.pix_sz_x))
        print(f"\tspiral size: {self.spiral_size}um")
        print(f"\tpix sz x: {self.pix_sz_x}um")
        print("\tradius (in pixels): {:.2f}px".format(radius_px * self.pix_sz_x))

        target_areas = []
        for coord in targetCoordinates:
            target_area = ([item for item in points_in_circle_np(radius_px, x0=coord[0], y0=coord[1])])
            target_areas.append(target_area)
        self.target_areas = target_areas

        ## target_areas that need to be excluded when filtering for nontarget cells
        radius_px = int(np.ceil(((self.spiral_size / 2) + 2.5) / self.pix_sz_x))
        print("\tradius of target exclusion zone (in pixels): {:.2f}px".format(radius_px * self.pix_sz_x))

        target_areas = []
        for coord in targetCoordinates:
            target_area = ([item for item in points_in_circle_np(radius_px, x0=coord[0], y0=coord[1])])
            target_areas.append(target_area)
        self.target_areas_exclude = target_areas

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

        self.target_coords = []
        counter = 0
        for slmgroup in target_files:
            target_image = tf.imread(os.path.join(naparm_path, slmgroup))
            targets = np.where(target_image > 0)
            targetCoordinates = list(zip(targets[1] * scale_factor, targets[0] * scale_factor))
            self.target_coords.append(targetCoordinates)
            print('\tNumber of targets (in SLM group %s): ' % (counter + 1), len(targetCoordinates))
            counter += 1


    def _findTargetedS2pROIs(self, plot: bool = True):
        """finding s2p cell ROIs that were also SLM targets (or more specifically within the target areas as specified by _findTargetAreas - include 15um radius from center coordinate of spiral)
        --- LAST UPDATED NOV 6 2021 - copied over from Rob's ---
        """

        '''
        Make a binary mask of the targets and multiply by an image of the cells
        to find cells that were targeted

        --- COPIED FROM ROB'S VAPE INTERAREAL_ANALYSIS.PY ON NOV 5 2021 ---
        '''

        print('searching for targeted cells... [Rob version]')
        ##### IDENTIFYING S2P ROIs THAT ARE WITHIN THE SLM TARGET SPIRAL AREAS
        # make all target area coords in to a binary mask
        targ_img = np.zeros([self.frame_x, self.frame_y], dtype='uint16')
        target_areas = np.array(self.target_areas)
        targ_img[target_areas[:, :, 1], target_areas[:, :, 0]] = 1

        # make an image of every cell area, filled with the index of that cell
        cell_img = np.zeros_like(targ_img)

        cell_y = np.array(self.cell_x)
        cell_x = np.array(self.cell_y)

        for i, coord in enumerate(zip(cell_x, cell_y)):
            cell_img[coord] = i + 1

        # binary mask x cell image to get the cells that overlap with target areas
        targ_cell = cell_img * targ_img

        targ_cell_ids = np.unique(targ_cell)[1:] - 1  # correct the cell id due to zero indexing
        self.targeted_cells = np.zeros([self.n_units], dtype='bool')
        self.targeted_cells[targ_cell_ids] = True
        # self.s2p_cell_targets = [self.cell_id[i] for i, x in enumerate(self.targeted_cells) if x is True]  # get ls of s2p cells that were photostim targetted
        self.s2p_cell_targets = [self.cell_id[i] for i in
                                 np.where(self.targeted_cells)[0]]  # get ls of s2p cells that were photostim targetted

        self.n_targeted_cells = np.sum(self.targeted_cells)

        print('------- Search completed.')
        self.save()
        print('Number of targeted cells: ', self.n_targeted_cells)

        ##### IDENTIFYING S2P ROIs THAT ARE WITHIN THE EXCLUSION ZONES OF THE SLM TARGETS
        # make all target area coords in to a binary mask
        targ_img = np.zeros([self.frame_x, self.frame_y], dtype='uint16')
        target_areas_exclude = np.array(self.target_areas_exclude)
        targ_img[target_areas_exclude[:, :, 1], target_areas_exclude[:, :, 0]] = 1

        # make an image of every cell area, filled with the index of that cell
        cell_img = np.zeros_like(targ_img)

        cell_y = np.array(self.cell_x)
        cell_x = np.array(self.cell_y)

        for i, coord in enumerate(zip(cell_x, cell_y)):
            cell_img[coord] = i + 1

        # binary mask x cell image to get the cells that overlap with target areas
        targ_cell = cell_img * targ_img

        targ_cell_ids = np.unique(targ_cell)[1:] - 1  # correct the cell id due to zero indexing
        self.exclude_cells = np.zeros([self.n_units], dtype='bool')
        self.exclude_cells[targ_cell_ids] = True
        self.s2p_cells_exclude = [self.cell_id[i] for i in
                                  np.where(self.exclude_cells)[0]]  # get ls of s2p cells that were photostim targetted

        self.n_exclude_cells = np.sum(self.exclude_cells)

        print('------- Search completed.')
        self.save()
        print(f"Number of exclude cells: {self.n_exclude_cells}")

        # define non targets from suite2p ROIs (exclude cells in the SLM targets exclusion - .s2p_cells_exclude)
        self.s2p_nontargets = [cell for cell in self.good_cells if
                               cell not in self.s2p_cells_exclude]  ## exclusion of cells that are classified as s2p_cell_targets

        print(f"Number of good, s2p non_targets: {len(self.s2p_nontargets)}")

        if plot:
            fig, ax = plt.subplots(figsize=[6, 6])

            targ_img = np.zeros([self.frame_x, self.frame_y], dtype='uint16')
            target_areas = np.array(self.target_areas)
            targ_img[target_areas[:, :, 1], target_areas[:, :, 0]] = 1
            ax.imshow(targ_img, cmap='Greys_r', zorder=0)
            ax.set_title('Targets areas')
            # for (x, y) in self.target_coords_all:
            #     ax.scatter(x=x, y=y, edgecolors='white', facecolors='none', linewidths=1.0)
            fig.show()

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

        mean_img_stack = np.zeros([end - start, self.frame_x, self.frame_y])
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
        self.raw_SLMTargets = targets_trace_full[:, curr_trial_frames[0] - start * 2000: curr_trial_frames[1] - (start * 2000)]

        self.dFF_SLMTargets = normalize_dff(self.raw_SLMTargets, threshold_pct=10)

        self.meanFluImg_registered = np.mean(mean_img_stack, axis=0)

        self.save() if save else None


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
            AssertionError('no stims set to analyse [1]')

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

    def _makeNontargetsStimTracesArray(self, stim_timings, normalize_to='pre-stim', save=True):
        """
        primary function to retrieve photostimulation trial timed Fluorescence traces for non-targets (ROIs taken from suite2p).
        :param self: alloptical experiment object
        :param normalize_to: str; either "baseline" or "pre-stim" or "whole-trace"
        :return: plot of avg_dFF of 100 randomly selected nontargets
        """
        print('\nCollecting peri-stim traces ')

        # collect photostim timed average dff traces of photostim targets
        dff_traces = []
        dff_traces_avg = []

        dfstdF_traces = []
        dfstdF_traces_avg = []

        raw_traces = []
        raw_traces_avg = []

        for cell in self.s2p_nontargets:
            # print('considering cell # %s' % cell)
            cell_idx = self.cell_id.index(cell)
            flu_trials = [self.raw[cell_idx][stim - self.pre_stim: stim + self.stim_duration_frames + self.post_stim]
                          for stim in stim_timings]

            dff_trace = normalize_dff(self.raw[cell_idx],
                                      threshold_pct=50)  # normalize trace (dFF) to mean of whole trace

            if normalize_to == 'baseline':  # probably gonna ax this anyways
                flu_dff = []
                mean_spont_baseline = np.mean(self.baseline_raw[cell_idx])
                for i in range(len(flu_trials)):
                    trace_dff = ((flu_trials[i] - mean_spont_baseline) / mean_spont_baseline) * 100

                    # add nan if cell is inside sz boundary for this stim
                    if hasattr(self, 'slmtargets_szboundary_stim'):
                        if self.is_cell_insz(cell=cell, stim=stim_timings[i]):
                            trace_dff = [np.nan] * len(flu_trials[i])

                    flu_dff.append(trace_dff)

            elif normalize_to == 'whole-trace':
                print('s2p neu. corrected trace statistics: mean: %s (min: %s, max: %s, std: %s)' %
                      (np.mean(self.raw[cell_idx]), np.min(self.raw[cell_idx]), np.max(self.raw[cell_idx]),
                       np.std(self.raw[cell_idx], ddof=1)))
                # dfstdf_trace = (self.raw[cell_idx] - np.mean(self.raw[cell_idx])) / np.std(self.raw[cell_idx], ddof=1)  # normalize trace (dFstdF) to std of whole trace
                flu_dfstdF = []
                flu_dff = []
                flu_dff_ = [dff_trace[stim - self.pre_stim: stim + self.stim_duration_frames + self.post_stim] for
                            stim in stim_timings if
                            stim not in self.seizure_frames]

                for i in range(len(flu_dff_)):
                    trace = flu_dff_[i]
                    mean_pre = np.mean(trace[0:self.pre_stim])
                    trace_dff = trace - mean_pre  # correct dFF of this trial to mean of pre-stim dFF
                    std_pre = np.std(trace[0:self.pre_stim], ddof=1)
                    dFstdF = trace_dff / std_pre  # normalize dFF of this trial by std of pre-stim dFF

                    flu_dff.append(trace_dff)
                    flu_dfstdF.append(dFstdF)

            elif normalize_to == 'pre-stim':
                flu_dff = []
                flu_dfstdF = []
                # print('|- splitting trace by photostim. trials and correcting by pre-stim period')
                for i in range(len(flu_trials)):
                    trace = flu_trials[i]
                    mean_pre = np.mean(trace[0:self.pre_stim])

                    std_pre = np.std(trace[0:self.pre_stim], ddof=1)
                    # dFstdF = (((trace - mean_pre) / mean_pre) * 100) / std_pre  # make dF divided by std of pre-stim F trace
                    dFstdF = (trace - mean_pre) / std_pre  # make dF divided by std of pre-stim F trace

                    if mean_pre < 1:
                        # print('risky cell here at cell # %s, trial # %s, mean pre: %s [1.1]' % (cell, i+1, mean_pre))
                        trace_dff = [np.nan] * len(trace)
                        dFstdF = [np.nan] * len(
                            trace)  # - commented out to test if we need to exclude cells for this correction with low mean_pre since you're not dividing by a bad mean_pre value
                    else:
                        # trace_dff = ((trace - mean_pre) / mean_pre) * 100
                        trace_dff = normalize_dff(trace, threshold_val=mean_pre)
                        # std_pre = np.std(trace[0:self.pre_stim], ddof=1)
                        # # dFstdF = (((trace - mean_pre) / mean_pre) * 100) / std_pre  # make dF divided by std of pre-stim F trace
                        # dFstdF = (trace - mean_pre) / std_pre  # make dF divided by std of pre-stim F trace

                    # # add nan if cell is inside sz boundary for this stim -- temporarily commented out for a while
                    # if 'post' in self.metainfo['trialType']:
                    #     if hasattr(self, 'slmtargets_szboundary_stim'):
                    #         if self.is_cell_insz(cell=cell, stim=stim_timings[i]):
                    #             trace_dff = [np.nan] * len(trace)
                    #             dFstdF = [np.nan] * len(trace)
                    #     else:
                    #         AttributeError(
                    #             'no slmtargets_szboundary_stim attr, so classify cells in sz boundary hasnot been saved for this expobj')

                    flu_dff.append(trace_dff)
                    flu_dfstdF.append(dFstdF)

            else:
                TypeError('need to specify what to normalize to in get_targets_dFF (choose "baseline" or "pre-stim")')

            dff_traces.append(flu_dff)  # contains all individual dFF traces for all stim times
            dff_traces_avg.append(np.nanmean(flu_dff, axis=0))  # contains the dFF trace averaged across all stim times

            dfstdF_traces.append(flu_dfstdF)
            dfstdF_traces_avg.append(np.nanmean(flu_dfstdF, axis=0))

            raw_traces.append(flu_trials)
            raw_traces_avg.append(np.nanmean(flu_trials, axis=0))

        if normalize_to == 'baseline':
            print(
                '\nCompleted collecting pre to post stim traces -- normalized to spont imaging as baseline -- for %s cells' % len(
                    dff_traces_avg))
            self.dff_traces = dff_traces
            self.dff_traces_avg = dff_traces_avg
            # return dff_traces, dff_traces_avg
        elif normalize_to == 'pre-stim' or normalize_to == 'whole-trace':
            print(
                f'\nCompleted collecting pre to post stim traces -- normalized to pre-stim period or maybe whole-trace -- for {len(dff_traces_avg)} cells, {len(flu_trials)} stims')
            self.dff_traces = np.asarray(dff_traces)
            self.dff_traces_avg = np.asarray([i for i in dff_traces_avg])
            self.dfstdF_traces = np.asarray(dfstdF_traces)
            self.dfstdF_traces_avg = np.asarray([i for i in dfstdF_traces_avg])
            self.raw_traces = np.asarray(raw_traces)
            self.raw_traces_avg = np.asarray([i for i in raw_traces_avg])

        print('\nFinished collecting peri-stim traces ')

        self.save() if save else None

        # return dff_traces, dff_traces_avg, dfstdF_traces, dfstdF_traces_avg, raw_traces, raw_traces_avg

    def _trialProcessing_nontargets(expobj, normalize_to='pre-stim', save=True):
        '''
        Uses dfstdf traces for individual cells and photostim trials, calculate the mean amplitudes of response and
        statistical significance across all trials for all cells

        Inputs:
            plane             - imaging plane n
        '''

        print('\n----------------------------------------------------------------')
        print('running trial Processing for nontargets ')
        print('----------------------------------------------------------------')

        # make trial arrays from dff data shape: [cells x stims x frames]
        expobj._makeNontargetsStimTracesArray(stim_timings=expobj.stim_start_frames, normalize_to=normalize_to,
                                              save=False)

        # create parameters, slices, and subsets for making pre-stim and post-stim arrays to use in stats comparison
        # test_period = expobj.pre_stim_response_window / 1000  # sec
        # expobj.test_frames = int(expobj.fps * test_period)  # test period for stats
        expobj.pre_stim_frames_test = np.s_[expobj.pre_stim - expobj.pre_stim_response_frames_window: expobj.pre_stim]
        stim_end = expobj.pre_stim + expobj.stim_duration_frames
        expobj.post_stim_frames_slice = np.s_[stim_end: stim_end + expobj.post_stim_response_frames_window]

        # mean pre and post stimulus (within post-stim response window) flu trace values for all cells, all trials
        expobj.analysis_array = expobj.dfstdF_traces  # NOTE: USING dF/stdF TRACES
        expobj.pre_array = np.mean(expobj.analysis_array[:, :, expobj.pre_stim_frames_test],
                                   axis=1)  # [cells x prestim frames] (avg'd taken over all stims)
        expobj.post_array = np.mean(expobj.analysis_array[:, :, expobj.post_stim_frames_slice],
                                    axis=1)  # [cells x poststim frames] (avg'd taken over all stims)

        # ar2 = expobj.analysis_array[18, :, expobj.post_stim_frames_slice]
        # ar3 = ar2[~np.isnan(ar2).any(axis=1)]
        # assert np.nanmean(ar2) == np.nanmean(ar3)
        # expobj.analysis_array = expobj.analysis_array[~np.isnan(expobj.analysis_array).any(axis=1)]

        # measure avg response value for each trial, all cells --> return array with 3 axes [cells x response_magnitude_per_stim (avg'd taken over response window)]
        expobj.post_array_responses = []  ### this and the for loop below was implemented to try to root out stims with nan's but it's likley not necessary...
        for i in np.arange(expobj.analysis_array.shape[0]):
            a = expobj.analysis_array[i][~np.isnan(expobj.analysis_array[i]).any(axis=1)]
            responses = a.mean(axis=1)
            expobj.post_array_responses.append(responses)

        expobj.post_array_responses = np.mean(expobj.analysis_array[:, :, expobj.post_stim_frames_slice], axis=2)
        expobj.wilcoxons = expobj._runWilcoxonsTest()

        expobj.save() if save else None


    #### STATISTICS FOR PHOTOSTIM RESPONSES
    def cellStaProcessing(self, test='t_test'):

        if self.stim_start_frames:

            # this is the key parameter for the sta, how many frames before and after the stim onset do you want to use
            self.pre_frames = int(np.ceil(self.fps * 0.5))  # 500 ms pre-stim period
            self.post_frames = int(np.ceil(self.fps * 3))  # 3000 ms post-stim period

            # ls of cell pixel intensity values during each stim on each trial
            self.all_trials = []  # ls 1 = cells, ls 2 = trials, ls 3 = dff vector

            # the average of every trial
            self.stas = []  # ls 1 = cells, ls 2 = sta vector

            self.all_amplitudes = []
            self.sta_amplitudes = []

            self.t_tests = []
            self.wilcoxons = []

            for plane in range(self.n_planes):

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
                        avg_post_end = avg_post_start + int(np.ceil(self.fps * 0.5))  # post-stim period of 500 ms

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

        for plane in range(self.n_planes):

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

        for plane in range(self.n_planes):

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

    def _runWilcoxonsTest(expobj, array1=None, array2=None):

        if array1 is None and array2 is None:
            array1 = expobj.pre_array;
            array2 = expobj.post_array

        # check if the two distributions of flu values (pre/post) are different
        assert array1.shape == array2.shape, 'shapes for expobj.pre_array and expobj.post_array need to be the same for wilcoxon test'
        wilcoxons = np.empty(len(array1))  # [cell (p-value)]

        for cell in range(len(array1)):
            wilcoxons[cell] = stats.wilcoxon(array2[cell], array1[cell])[1]

        return wilcoxons

        # expobj.save() if save else None

    def _sigTestAvgResponse_nontargets(expobj, p_vals=None, alpha=0.1, save=True):
        """
        Uses the p values and a threshold for the Benjamini-Hochberg correction to return which
        cells are still significant after correcting for multiple significance testing
        """
        print('\n----------------------------------------------------------------')
        print('running statistical significance testing for nontargets response arrays ')
        print('----------------------------------------------------------------')

        # p_vals = expobj.wilcoxons
        sig_units = np.full_like(p_vals, False, dtype=bool)

        try:
            sig_units, _, _, _ = sm.stats.multitest.multipletests(p_vals, alpha=alpha, method='fdr_bh',
                                                                  is_sorted=False, returnsorted=False)
        except ZeroDivisionError:
            print('no cells responding')

        # # p values without bonferroni correction
        # no_bonf_corr = [i for i, p in enumerate(p_vals) if p < 0.05]
        # expobj.nomulti_sig_units = np.zeros(len(expobj.s2p_nontargets), dtype='bool')
        # expobj.nomulti_sig_units[no_bonf_corr] = True

        # expobj.save() if save else None

        # p values after bonferroni correction
        #         bonf_corr = [i for i,p in enumerate(p_vals) if p < 0.05 / expobj.n_units[plane]]
        #         sig_units = np.zeros(expobj.n_units[plane], dtype='bool')
        #         sig_units[bonf_corr] = True

        return sig_units

    # other useful functions for all-optical analysis
    def whiten_photostim_frame(self, tiff_path, save_as=''):
        im_stack = tf.imread(tiff_path, key=range(self.n_frames))

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
                im_stack = tf.imread(tiff_path, key=range(self.n_frames))
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
        frameRate = int(self.fps)

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




class WideFieldImaging:
    """
    WideField imaging data object.

    """

    def __init__(self, tiff_path, paq_path, exp_metainfo, pkl_path):
        self.tiff_path = tiff_path
        self.paq_path = paq_path
        self.metainfo = exp_metainfo
        # set and create analysis save path directory
        self.analysis_save_dir = self.pkl_path[:[(s.start(), s.end()) for s in re.finditer('/', self.pkl_path)][-1][0]]
        if not os.path.exists(self.analysis_save_dir):
            print('making analysis save folder at: \n  %s' % self.analysis_save_dir)
            os.makedirs(self.analysis_save_dir)

        self.save_pkl(pkl_path=pkl_path)  # save experiment object to pkl_path

    def __repr__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
        else:
            lastmod = "(unsaved pkl object)"
        if not hasattr(self, 'metainfo'):
            information = f"uninitialized"
        else:
            information = self.t_series_name

        return repr(f"({information}) WidefieldImaging experimental data object, last saved: {lastmod}")

    @property
    def t_series_name(self):
        if 't series id' in self.metainfo.keys():
            return f"{self.metainfo['t series id']}"
        elif "animal prep." in self.metainfo.keys() and "trial" in self.metainfo.keys():
            return f'{self.metainfo["animal prep."]} {self.metainfo["trial"]}'
        else:
            raise ValueError('no information found to retrieve t series id')

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
