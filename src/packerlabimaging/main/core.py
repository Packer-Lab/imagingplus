# this file contains the two fundamental class types (Trial and Experiment) needed to construct an experiment in packerlabimaging
from __future__ import absolute_import
from dataclasses import dataclass, field
from typing import Optional, MutableMapping, Union, TypedDict, List, Dict, Any

import numpy as np
import pandas as pd

from packerlabimaging.main.subcore import ImagingMetadata, ImagingData, CellAnnotations, TemporalData
from packerlabimaging.processing.anndata import AnnotatedData

# from packerlabimaging.processing.imagingMetadata import ImagingMetadata
from packerlabimaging.processing.suite2p import Suite2pResultsTrial
from packerlabimaging.utils.classes import UnavailableOptionError, TrialMetainfo

import os
import time
import re
import tifffile as tf

import pickle


# add new class for temporal synchronization of additional 1d dataarrays with imaging cellsdata - parent of Paq
# todo - test new paqdata child class.

# add new class for cell annotations

# TODO test out new Trial class and new workflows
# todo test anndata creation from Trial

# TODO add function in Experiment object for merging of Trials across time axis (assert they have the same cell annotations).

# TODO [ ]  thinking about restructuring TwoPhotonImaging trial methods to more general trial type
#    [ ]  then making TwoPhotonImaging as an independent workflow, allowing the general trial type to retain the methods that are *currently* in TwoPhotonImaging


@dataclass
class SingleImage:
    dataPath: str
    date: str = None
    imgID: str = None
    expGroup: str = None
    comment: str = None
    imparams: ImagingMetadata = None

    def __post_init__(self):
        self.data = tf.imread(self.dataPath)

    def plotImg(self, **kwargs):
        from packerlabimaging.plotting.plotting import SingleTiffFrame
        SingleTiffFrame(imstack=self.data, **kwargs)


# noinspection DuplicatedCode
@dataclass
class Experiment:
    """Overall experiment. It can collate all imaging trials that are part of a single field-of-view (FOV) of imaging"""

    expID: str  #: given identification name for experiment
    dataPath: str  #: main dir where the imaging cellsdata is contained
    saveDir: str  #: main dir where the experiment object and the trial objects will be saved to
    comment: str = ''  #: notes related to experiment
    singleImages: Dict[str, SingleImage] = None  #: contains single image frames from each experiment

    def __post_init__(self):
        print(f'***********************')
        print(f'CREATING new Experiment: (expID: {self.expID})')
        print(f'***********************\n\n')

        self.singleImages = {}

        # self.TrialsInformation: MutableMapping[str, Union[str, TrialsInformation, PaqInfo]] = {}  #: dictionary of metadata information about each trial. Gets filled while adding each trial.
        self.TrialsInformation: MutableMapping[
            str, TrialMetainfo] = {}  #: dictionary of metadata information about each trial. Gets filled while adding each trial.

        # suite2p related attrs initialization
        from packerlabimaging.processing.suite2p import Suite2pExperiment
        self.Suite2p: Suite2pExperiment = None
        self._suite2p_save_path = None
        # self._trialsTiffsSuite2p = {}  #: dictionary of trial IDs and their respective .tiff paths for each trial that will be used in Suite2p processing for current experiment
        # self._s2pResultExists = False  #: flag for whether suite2p results exist for current experiment

        # save Experiment object
        self._pkl_path = self._get_save_location()
        os.makedirs(self.saveDir, exist_ok=True)
        self.save_pkl(pkl_path=self.pkl_path)

        print(f'\n\n\n******************************')
        print(f"NEW Experiment created: \n")
        print(self.__str__())

    def __repr__(self):
        return f"packerlabimaging.Experiment object (expID: {self.expID})"

    def __str__(self):
        lastsaved = time.ctime(os.path.getmtime(self.pkl_path))
        __return_information = f"packerlabimaging Experiment object (last saved: {lastsaved}), expID: {self.expID}"
        __return_information = __return_information + f"\nfile path: {self.pkl_path}"

        if len(self.TrialsInformation) > 0:
            __return_information = __return_information + f"\n\ntrials in Experiment object:"
            for trial in self.TrialsInformation:
                __return_information = __return_information + self.get_trial_infor(trialID=trial)
            return f"{__return_information}\n"
        else:
            return f"{__return_information}\n"

    def _get_save_location(self):
        """specify location to create pickled experiment object."""
        if self.saveDir[-4:] == '.pkl':
            _pkl_path = self.saveDir
            self.saveDir = self.saveDir[
                           :[(s.start(), s.end()) for s in re.finditer('/', self.saveDir)][-1][0]]
        else:
            self.saveDir = self.saveDir + '/' if self.saveDir[
                                                     -1] != '/' else self.saveDir
            _pkl_path = f"{self.saveDir}{self.expID}_analysis.pkl"
        os.makedirs(self.saveDir, exist_ok=True)
        return _pkl_path

    def get_trial_infor(self, trialID: str):
        if trialID in [*self.TrialsInformation]:
            return f"\n\t{trialID}: {self.TrialsInformation[trialID]['expGroup']}"
        else:
            ValueError(f"{trialID} not found in Experiment.")

    # def add_trial(self, trialobj: ImagingTrial):
    #     """
    #     Add trial object to the experiment. This will add metainformation about the trial to the experiment.
    #
    #     :param trialobj: ImagingTrial instance.
    #     """
    #
    #     print(f"\n\n\- ADDING trial: {trialobj.trialID} to {self.expID} experiment", end='\r')
    #
    #     self.TrialsInformation[trialobj.trialID] = trialobj.metainfo
    #     self.TrialsInformation[trialobj.trialID]['expGroup'] = trialobj.metainfo['expGroup']
    #     for attr in dir(trialobj):
    #         if 'path' in attr:  # add any attr's containing 'path' from trialobj to the TrialsInformation dict
    #             self.TrialsInformation[trialobj.trialID]['paths'][attr] = getattr(trialobj, attr)
    #
    #     print(f"|- ADDED trial: {trialobj.trialID} to {self.expID} experiment")

    def add_trial(self, trialID, trialobj=None, **kwargs):
        """
        Add trial object to the experiment. This will add metainformation about the trial to the experiment.

        :param trialID:
        :param trialobj: ImagingTrial instance.
        """

        print(f"\n\n\- ADDING trial: {trialID} to {self.expID} experiment", end='\r')

        if trialobj is None:
            trialobj = ImagingTrial.newImagingTrialfromExperiment(experiment=self, trialID=trialID, **kwargs)

        self.TrialsInformation[trialobj.trialID] = trialobj.metainfo
        self.TrialsInformation[trialobj.trialID]['expGroup'] = trialobj.metainfo['expGroup']
        for attr in dir(trialobj):
            if 'path' in attr:  # add any attr's containing 'path' from trialobj to the TrialsInformation dict
                self.TrialsInformation[trialobj.trialID]['paths'][attr] = getattr(trialobj, attr)

        self.save()
        print(f"|- ADDED trial: {trialobj.trialID} to {self.expID} experiment")

    def add_single(self, singleimg: SingleImage, **kwargs):
        """

        :param singleimg:
        :param kwargs:
        """
        setattr(self, singleimg.imgID, singleimg)

    def combine_trials(self):
        """todo: Combine anndata table of trials with same cells."""

    @property
    def trialIDs(self):
        return [*self.TrialsInformation]

    @staticmethod
    def tiff_path_dir(tiff_path):
        """retrieve the parent directory of the provided tiff_path"""
        return os.path.dirname(tiff_path)
        # return tiff_path[:[(s.start(), s.end()) for s in re.finditer('/', tiff_path)][-1][
        #     0]]  # this is the directory where the Bruker xml files associated with the 2p imaging TIFF are located

    @property
    def pkl_path(self):
        """path in Analysis folder to save pkl object"""
        return self._pkl_path

    # noinspection PyAttributeOutsideInit
    @pkl_path.setter
    def pkl_path(self, path: str):
        self._pkl_path = path

    def save_pkl(self, pkl_path: str = None):
        """
        saves packerlabimaging Experiment object using Pickle library (.pkl file extension).

        :param pkl_path: full .pkl s2pResultsPath for saving the Experiment object.
        """
        if pkl_path:
            print(f'\nsaving new pkl object at: {pkl_path}')
            self.pkl_path = pkl_path

        with open(self.pkl_path, 'wb') as f:
            pickle.dump(self, f)
        print(f"\n\t|- Experiment analysis object saved to {self.pkl_path} -- ")

    def save(self):
        self.save_pkl()

    def load_trial(self, trialID: str):
        """
        method for importing individual trial objects from Experiment instance using the trial id for a given trial

        :param trialID:
        :return:


        """
        try:
            trial_pkl_path = self.TrialsInformation[trialID]['paths']['pkl_path']
            from packerlabimaging import import_obj
            trialobj = import_obj(trial_pkl_path)
            return trialobj
        except KeyError:
            raise KeyError(f"trialID: {trialID} not found in Experiment instance.")

    def add_suite2p(self, s2p_trials: Union[list, str] = 'all', s2pResultsPath: Optional[str] = None):
        """
        todo: try moving function to suite2p file.
        Wrapper for adding suite2p results to Experiment. Can only be run after adding trials to Experiment.

        :param s2p_trials: list of trials to use for Suite2p processing/analysis pipeline, default = 'all' to use all trials of Experiment
        :param s2pResultsPath: optional, if suite2p already run then here provide path to plane0 folder of existing suite2p results
        """

        print(f'\- Adding suite2p module to experiment. Located under .Suite2p')
        from packerlabimaging.processing.suite2p import Suite2pExperiment, Suite2pResultsTrial

        if not len([*self.TrialsInformation]) > 0:
            raise UnavailableOptionError(
                'need to add at least 1 trial to Experiment before adding Suite2p functionality.')

        _trialsTiffsSuite2p = {}

        if s2p_trials == 'all': s2p_trials = self.trialIDs
        assert type(s2p_trials) == list and len(s2p_trials) > 0, 'no s2p trials given in list form.'
        for trial in s2p_trials: _trialsTiffsSuite2p[trial] = self.TrialsInformation[trial]['paths']['dataPath']

        if s2pResultsPath:  # if s2pResultsPath provided then try to find and pre-load results from provided s2pResultsPath, raise error if cannot find results
            # search for suite2p results items in self.suite2pPath, and auto-assign s2pRunComplete --> True if found successfully
            __suite2p_path_files = os.listdir(s2pResultsPath)
            self._s2pResultExists = False
            for filepath in __suite2p_path_files:
                if 'ops.npy' in filepath:
                    self._s2pResultExists = True
                    break
            if self._s2pResultExists:
                self._suite2p_save_path = s2pResultsPath
                self.Suite2p = Suite2pExperiment(trialsTiffsSuite2p=_trialsTiffsSuite2p,
                                                 s2pResultsPath=s2pResultsPath)
            else:
                raise ValueError(
                    f"suite2p results could not be found. `suite2pPath` provided was: {s2pResultsPath}")
        else:  # no s2pResultsPath provided, so initialize without pre-loading any results
            self.Suite2p = Suite2pExperiment(trialsTiffsSuite2p=_trialsTiffsSuite2p)

        # print(self.Suite2p)
        # adding of suite2p trial level as well in this function as well
        total_frames = 0
        for trial in s2p_trials:
            trialobj = self.load_trial(trialID=trial)
            # print(f'n_frames', self.Suite2p.n_frames)
            trialobj.Suite2p = Suite2pResultsTrial(s2pExp=self.Suite2p, trial_frames=(
                total_frames, total_frames + trialobj.n_frames))  # use trial obj's current trial key_frames
            trialobj.save()
            total_frames += trialobj.n_frames

        print(f'|- Finished adding suite2p module to experiment and trials. Located under .Suite2p')
        self.save()

    @property
    def suite2p_save_path(self):
        return self._suite2p_save_path


@dataclass
class ImagingTrial:
    dataPath: str
    saveDir: str  #: main dir where the experiment object and the trial objects will be saved to
    date: str
    trialID: str
    expID: str
    expGroup: str = ''
    comment: str = ''
    imparams: ImagingMetadata = None
    data: AnnotatedData = None
    imdata: ImagingData = None
    cells: CellAnnotations = None
    tmdata: TemporalData = None
    Suite2p: Suite2pResultsTrial = None

    def __post_init__(self):
        self.metainfo = TrialMetainfo(date=self.date, trialID=self.trialID, expID=self.expID, expGroup=self.expGroup,
                                      comment=self.comment, paths={})  # , microscope=self.microscope)

        if os.path.exists(self.dataPath):
            self.metainfo['paths']['dataPath'] = self.dataPath
        else:
            raise FileNotFoundError(f"dataPath does not exist: {self.dataPath}")

        # processing collect mean FOV Trace -- after collecting imaging params and Paq timing info
        im_stack = self.importTrialTiff()
        self._n_frames = im_stack.shape[0]
        self.meanFluImg, self.meanFovFluTrace = self.meanRawFluTrace(
            im_stack)  #: mean image and mean FOV fluorescence trace

        # set, create analysis save path directory and create pkl object
        os.makedirs(self.saveDir, exist_ok=True)
        self.metainfo['paths']['pkl_path'] = f"{self.saveDir}{self.date}_{self.trialID}.pkl"
        self.save_pkl(pkl_path=self.pkl_path)  # save experiment object to pkl_path
        self.metainfo['paths']['pkl_path'] = self.pkl_path

        # make anndata
        if self.imdata and self.cells and self.tmdata:
            self.create_anndata()

    def __repr__(self):
        # todo test repr
        return repr(f"{self.t_series_name} (ImagingTrial)")

    # todo maybe ? - add alternative constructor to handle construction if temporal cellsdata or cell annotations or imaging cellsdata is provided
    # might be useful for providing like submodules (especially things like Suite2p)
    @classmethod
    def newImagingTrialfromExperiment(cls, experiment: Experiment, trialID, dataPath, date, comment=''):
        """
        Creates a new ImagingTrial from an Experiment.

        :param trialID:
        :param dataPath:
        :param date:
        :param group:
        :param comment:
        :param experiment:
        :param kwargs: Same args as ImagingTrial
        :return: ImagingTrial
        """

        trialobj: ImagingTrial = cls(dataPath=dataPath,
                                     saveDir=experiment.saveDir,
                                     date=date,
                                     trialID=trialID,
                                     expID=experiment.expID,
                                     comment=comment)
        experiment.add_trial(trialID=trialID, trialobj=trialobj)
        return trialobj

    # @property
    # def date(self):
    #     """date of the experiment datacollection"""
    #     return self.metainfo['date']

    # @property
    # def microscope(self):
    #     """name of imaging cellsdata acquisition microscope system"""
    #     return self.metainfo['microscope']

    # @property
    # def expID(self):
    #     """experiment ID of current trial object"""
    #     return self.metainfo['expID']

    # @property
    # def trialID(self):
    #     """trial ID of current trial object"""
    #     return self.metainfo['trialID']

    @property
    def tiff_path(self):
        """tiff path of current trial object"""
        return self.metainfo['paths']['dataPath']

    @property
    def t_series_name(self):
        if "expID" in self.metainfo.keys() and "trialID" in self.metainfo.keys():
            return f'{self.metainfo["expID"]} {self.metainfo["trialID"]}'
        else:
            raise ValueError('no information found to retrieve t series id')

    @property
    def tiff_path_dir(self):
        """Parent directory of .tiff cellsdata path"""
        return os.path.dirname(self.tiff_path)

    @property
    def pkl_path(self):
        """path in Analysis folder to save pkl object"""
        return self.metainfo['paths']['pkl_path']

    @pkl_path.setter
    def pkl_path(self, path: str):
        self.metainfo['paths']['pkl_path'] = path

    def save_pkl(self, pkl_path: str = None):
        """
        Alias method for saving current object to pickle file.

        :param pkl_path: (optional) provide path to save object to pickle file.
        """
        if pkl_path:
            parent = os.path.relpath(os.path.join(pkl_path, os.pardir))
            if os.path.exists(parent):
                print(f'saving new trial object to: {pkl_path}')
                self.pkl_path = pkl_path
            else:
                raise FileNotFoundError(f"Parent directory path: `{parent}` was not found for saving .pkl")

        with open(self.pkl_path, 'wb') as f:
            pickle.dump(self, f)
        print("\n\t -- data object saved to %s -- " % self.pkl_path)

    def save(self):
        """
        Alias method for saving current object as pickle file.
        """
        self.save_pkl()

    def plotSingleImageFrame(self, frame_num: int = 0, title: str = None):
        """
        plots an image of a single specified tiff frame after reading using tifffile.
        :param frame_num: frame # from 2p imaging tiff to show (default is 0 - i.e. the first frame)
        :param title: (optional) give a string to use as title
        :return: matplotlib imshow plot
        """
        from packerlabimaging.plotting.plotting import SingleTiffFrame
        stack = SingleTiffFrame(tiff_path=self.tiff_path, frame_num=frame_num, title=title)

        # stack = tf.imread(self.tiff_path, key=frame_num)
        # plt.imshow(stack, cmap='gray')
        # plt.suptitle(title) if title is not None else plt.suptitle(f'frame num: {frame_num}')
        # plt.show()
        return stack

    def addImagingMetadata(self, microscope: str = '', imaging_metadata: ImagingMetadata = None):
        """Add the imaging metadata submodule to the current ImagingTrial."""
        # imaging metadata
        if 'Bruker' in microscope:
            from packerlabimaging.processing.imagingMetadata import PrairieViewMetadata
            self.imparams = PrairieViewMetadata(pv_xml_dir=self.tiff_path_dir)
        elif imaging_metadata:
            self.imparams = imaging_metadata
        else:
            Warning(
                f"NO imaging microscope parameters set. follow imagingMetadata to create a custom imagingMicroscopeMetadata class.")

    ## below properties/methods have pre-requisite processing steps
    @property
    def n_frames(self):
        """number of imaging frames from tiff cellsdata.
        """
        return self._n_frames

    def importTrialTiff(self) -> np.ndarray:
        """
        Import current trial's imaging tiff in full.

        :return: imaging tiff as numpy array
        """
        print(f"\n\- loading raw TIFF file from: {self.tiff_path}", end='\r')
        im_stack = tf.imread(self.tiff_path)
        print('|- Loaded experiment tiff of shape: ', im_stack.shape)

        return im_stack

    def meanRawFluTrace(self, im_stack: np.ndarray = None):
        """
        Collects the raw mean of FOV fluorescence trace across the t-series.

        :return: mean fluorescence trace
        """
        try:
            im_stack = self.importTrialTiff() if im_stack is None else im_stack

            print('\n-----collecting mean raw flu trace from tiff file...')
            mean_flu_img = np.mean(im_stack, axis=0)
            mean_fov_flutrace = np.mean(np.mean(im_stack, axis=1), axis=1)

            return mean_flu_img, mean_fov_flutrace
        except AssertionError:
            raise Warning('no imaging metadata module found for trial. set imaging metadata to `.imparams`')

    def makeDownsampledTiff(self):
        """Import current trial tiff, create downsampled tiff and save in default analysis directory."""

        stack = self.importTrialTiff()
        from packerlabimaging.utils.utils import SaveDownsampledTiff
        SaveDownsampledTiff(stack=stack, save_as=f"{self.saveDir}/{self.date}_{self.trialID}_downsampled.tif")

    def create_anndata(self, imdata_type: str = None, layers=False):
        """
        Creates annotated cellsdata (see anndata library for more information on AnnotatedData) object based around the Ca2+ matrix of the imaging trial.

        """
        if not self.imdata and self.cells and self.tmdata:
            raise ValueError(
                'cannot create anndata table. anndata creation only available if experiments have ImagingData (.imdata), CellAnnotations (.cells) and TemporalData (.tmdata)')

        # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
        assert hasattr(self.cells, 'cellsdata'), 'missing cellsdata attr from .cells'
        obs_meta = self.cells.cellsdata

        # SETUP THE VARIABLES ANNOTATIONS TO USE IN anndata
        assert hasattr(self.tmdata, 'cellsdata'), 'missing cellsdata attr from .tmdata'
        var_meta = self.tmdata.data

        print(f"\n\----- CREATING annotated cellsdata object using AnnData:")
        assert hasattr(self.imdata, 'cellsdata'), 'missing cellsdata attr from .imdata'
        _data_type = [*self.imdata.imdata][0] if not imdata_type else imdata_type
        primary_data = self.imdata.imdata[_data_type]

        if layers and self.imdata.n_data > 0:
            layers = {}
            for layer in self.imdata.data_labels:
                if layer == _data_type:
                    layers[layer] = self.imdata.imdata[layer]
        else:
            layers = None

        anndata_setup = {'X': primary_data, 'data_label': _data_type, 'obs': obs_meta, 'var': var_meta,
                         'obs_m': self.cells.multidim_data if self.cells.multidim_data else None, 'layers': layers}

        adata = AnnotatedData(**anndata_setup)

        print(f"\n{adata}")
        return adata

    # todo test
    def collectSignalFromCoords(self, curr_trial_frames: Union[tuple, list], target_coords_masks: np.ndarray,
                                reg_tif_folder: str = None, save: bool = True):
        """uses registered tiffs to collect raw traces from provided target masks target areas

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

        num_coords = len(target_coords_masks)
        target_areas = target_coords_masks

        print(
            f'\n\ncollecting raw Flu traces from SLM target coord. areas from registered TIFFs from: {reg_tif_folder}')
        # read in registered tiff
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        start = curr_trial_frames[0] // 2000  # 2000 because that is the batch size for suite2p run
        end = curr_trial_frames[1] // 2000 + 1

        mean_img_stack = np.zeros([end - start, self.imparams.frame_x, self.imparams.frame_y])
        # collect mean traces from target areas of each target coordinate by reading in individual registered tiffs that contain frames for current trial
        targets_trace_full = np.zeros([num_coords, (end - start) * 2000], dtype='float32')
        counter = 0
        for i in range(start, end):
            tif_path_save2 = self.Suite2p.s2pResultsPath + '/reg_tif/' + reg_tif_list[i]
            with tf.TiffFile(tif_path_save2, multifile=False) as input_tif:
                print('|- reading tiff: %s' % tif_path_save2)
                data = input_tif.asarray()

            targets_trace = np.zeros([num_coords, data.shape[0]], dtype='float32')
            for coord in range(num_coords):
                x = data[:, target_areas[coord, :, 1], target_areas[coord, :, 0]]
                targets_trace[coord] = np.mean(x, axis=1)

            targets_trace_full[:, (i - start) * 2000: ((i - start) * 2000) + data.shape[
                0]] = targets_trace  # iteratively write to each successive segment of the targets_trace array based on the length of the reg_tiff that is read in.

            mean_img_stack[counter] = np.mean(data, axis=0)
            counter += 1

        # final part, crop to the *exact* frames for current trial
        raw_coords_signal = targets_trace_full[:,
                         curr_trial_frames[0] - start * 2000: curr_trial_frames[1] - (start * 2000)]

        dFF_coords_signal = self.normalize_dff(raw_coords_signal, threshold_pct=10)

        meanFluImg_registered = np.mean(mean_img_stack, axis=0)

        self.save() if save else None

        return raw_coords_signal, dFF_coords_signal, meanFluImg_registered


    @staticmethod
    def normalize_dff(arr, threshold_pct=20, threshold_val=None):
        """normalize given array (cells x time) to the mean of the fluorescence values below given threshold. Threshold
        will refer to the that lower percentile of the given trace."""

        if arr.ndim == 1:
            if threshold_val is None:
                a = np.percentile(arr, threshold_pct)
                mean_ = arr[arr < a].mean()
            else:
                mean_ = threshold_val
            new_array = ((arr - mean_) / mean_) * 100
            if np.isnan(new_array).any() == True:
                Warning('Cell (unknown) contains nan, normalization factor: %s ' % mean_)

        else:
            new_array = np.empty_like(arr)
            for i in range(len(arr)):
                if threshold_val is None:
                    a = np.percentile(arr[i], threshold_pct)
                else:
                    a = threshold_val
                mean_ = np.mean(arr[i][arr[i] < a])
                new_array[i] = ((arr[i] - mean_) / abs(mean_)) * 100

                if np.isnan(new_array[i]).any() == True:
                    print('Warning:')
                    print('Cell %d: contains nan' % (i + 1))
                    print('      Mean of the sub-threshold for this cell: %s' % mean_)

        return new_array
