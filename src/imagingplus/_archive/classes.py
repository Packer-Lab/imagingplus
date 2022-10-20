### TO BE ARCHIVED -- REFACTORED INTO CORE.PY AND SUBCORE.PY ON MAY 10 2022
# this file contains the two fundamental class types (Trial and Experiment) needed to construct an experiment in imagingplus
from __future__ import absolute_import
from dataclasses import dataclass, field
from typing import Optional, MutableMapping, Union, TypedDict, List, Dict, Any

import numpy as np
import pandas as pd
from packerlabimaging.processing.anndata import AnnotatedData

# from imagingplus.processing.imagingMetadata import ImagingMetadata
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

# noinspection DuplicatedCode
@dataclass
class Experiment:
    """Overall experiment. It can collate all imaging trials that are part of a single field-of-view (FOV) of imaging"""

    expID: str  #: given identification name for experiment
    dataPath: str  #: main dir where the imaging cellsdata is contained
    saveDir: str  #: main dir where the experiment object and the trial objects will be saved to
    comment: str = ''  #: notes related to experiment

    def __post_init__(self):
        print(f'***********************')
        print(f'CREATING new Experiment: (expID: {self.expID})')
        print(f'***********************\n\n')

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
        return f"imagingplus.Experiment object (expID: {self.expID})"

    def __str__(self):
        lastsaved = time.ctime(os.path.getmtime(self.pkl_path))
        __return_information = f"imagingplus Experiment object (last saved: {lastsaved}), expID: {self.expID}"
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
        saves imagingplus Experiment object using Pickle library (.pkl file extension).

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
        """method for importing individual trial objects from Experiment instance using the trial id for a given trial"""
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
class TemporalData:
    """1-D time series datsets corresponding with an imaging trial."""
    file_path: str  #: path to cellsdata file
    sampling_rate: float  #: rate of cellsdata collection (Hz)
    channels: List[str]  #: list of cellsdata channel names.
    data: pd.DataFrame  #: N columns x Time array of N x 1D cellsdata channels collected over Time at the same sampling rate
    frame_times: Union[
        list, np.ndarray] = None  #: timestamps representing imaging frame times. must be of same time duration as imaging dataset.
    sparse_data: Dict[
        str, np.ndarray] = None  #: dictionary of timeseries channels containing cellsdata keyed at imaging frames for each timeseries channel

    def __post_init__(self):
        print(f"Created new TemporalData of {self.n_channels} x {self.n_timepoints} (sampled at {self.sampling_rate}")
        pass

    def __repr__(self):
        return f"TemporalData module, containing {self.n_timepoints} timepoints and {self.n_channels} channels."

    def __str__(self):
        return f"TemporalData module, containing {self.n_timepoints} timepoints and {self.n_channels} channels: \n{self.channels}"

    @property
    def n_frames(self):
        """number of timed imaging frames"""
        return len(self.frame_times)

    @property
    def n_timepoints(self):
        """number of cellsdata collection timepoints"""
        return self.data.shape[0]

    @property
    def n_channels(self):
        """number of cellsdata collection channels"""
        return self.data.shape[1]

    def get_sparse_data(self, frame_times: Union[list, np.ndarray] = None):
        """
        todo: need to probably average the original signal between key_frames if collected at a higher rate than frame_times.
            - why not just use a downsampling algorithm to downsample to the same frame rate and num datapoints as the imaging key_frames????
        Returns dictionary of numpy array keyed on channels from paq_data timed to 2photon imaging frame_times.
        In effect this works as a downsampling algorithm.

        :param frame_times:
        :return:
        """

        # todo insert test to check that original signal has been collected at a rate higher than imaging. if not then need to handle differently.

        assert hasattr(self,
                       'frame_times') or frame_times, 'no frame_times given to retrieve cellsdata from those timestamps.'

        frame_times = self.frame_times if frame_times is None else frame_times

        print(f"\n\t\- Getting imaging key frames timed cellsdata from {len(frame_times)} frames ... ")

        # read in and save sparse version of all cellsdata channels (only save cellsdata from timepoints at frame clock times)
        sparse_data = {}
        for idx, chan in enumerate(self.channels):
            print(f'\t\t\- Adding sparse cellsdata for channel: {chan} ')
            data = self.data.iloc[frame_times, idx]
            sparse_data[chan] = data

        sparse_data = pd.DataFrame(sparse_data)

        return sparse_data

    def cropData(self, begin: int, end: int, channels: List[str] = 'all', replace: bool = False):
        """
        Crops saved temporal cellsdata channels to the timestamps of begin and end.

        :param channels: channels to crop cellsdata, default is 'all' (all channels in dataset).
        :param begin: timestamp to begin cropping at
        :param end: timestamp to end cropping at
        :param replace: if True, then replace .cellsdata attribute with newly cropped cellsdata. if False, return the cropped dataframe.
        """

        channels = self.channels if channels == 'all' else channels
        print(f"\- cropping {channels} to {begin} and {end} paq clock times.")

        assert self.data.shape[0] >= (
                end - begin), f'Not enough time series samples in cellsdata to crop between the provided times.'
        data_ = self.data.loc[begin: end, channels]

        if replace:
            self.data = data_
        else:
            return data_

        # older approach...
        # for channel in channels:
        #     print(f"\- cropping {channel} to {begin} and {end} paq clock times.")
        #     cellsdata = self.cellsdata[channel].to_numpy()
        #     assert len(cellsdata) >= (end - begin), f'{channel} paq cellsdata is not long enough to crop between the provided clock times.'
        #     cropdata = cellsdata[begin: end]
        #     setattr(self, channel, cropdata)


@dataclass
class CellAnnotations:
    """Annotations of cells in an imaging trial."""

    def __init__(self, cells_array: Union[List[int], pd.Index, pd.RangeIndex, np.ndarray], annotations: Union[List[str], pd.Index],
                 cellsdata: Union[pd.DataFrame, pd.Series], multidim_data: Dict[str, List[Any]] = None):
        """

        :param cells_array: ID of all cells in imaging dataset. must be of same cell length as imaging dataset.
        :param annotations: list of names of annotations
        :param cellsdata: M x Cells array of an arbritrary number (M) 1D annotations channels collected for all Cells. must contain same number of cells as cells_array.
        :param multidim_data: annotations with cellsdata of unconstrained dimensions for all cells. Structured as dictionary with keys corresponding to annotation name and a list of the length of cells containing cellsdata in any format.
        """

        self.cells_array = cells_array
        self.annotations = annotations
        self.cellsdata = cellsdata
        self.multidim_data = multidim_data

        print(f'\- added CellAnnotations module. consisting of {self.n_annotations} annotations.')
        if self.multidim_data:
            for label, data in self.multidim_data.items():
                if not len(data) == self.n_cells:
                    raise ValueError(f"length of {label} of multidimdata does not match number of cells.")

    def __repr__(self):
        return f"CellAnnotations module, containing {self.n_cells} cells and {self.n_annotations} annotations."

    def __str__(self):
        return f"CellAnnotations module, containing {self.n_cells} cells and {self.n_annotations} annotations: \n{self.annotations}"

    # constructors:
    # @classmethod
    # def s2p_to_CellAnnotations(cls, s2pTrial):
    #     """alternative constructor for CellAnnotations from suite2p results."""
    #     cellsdata = s2pTrial.setCellsAnnotations()
    #     obj = cls(cells_array=cellsdata['original_index'], annotations=cellsdata.columns, cellsdata=cellsdata)
    #     return obj

    # properties:
    @property
    def n_cells(self):
        """number of cells"""
        return len(self.cells_array)

    @property
    def n_annotations(self):
        """number of annotations"""
        return len(self.annotations)

    # todo ?create properties for accessing individual annotations

    @property
    def cell_id(self):
        """ID of cells"""
        assert 'cell_id' in self.cellsdata, 'cell_id cannot be found in cells annotations under cellsdata'
        return list(self.cellsdata['cell_id'])

    # functions:


@dataclass
class ImagingData:
    """Imaging dataset."""

    def __init__(self, imdata: Union[np.ndarray, pd.DataFrame], **kwargs):
        """
        Imaging dataset.

        :param imdata: N rois x num_frames, table of imaging cellsdata of cells (N) collected over time (Frames)
        """
        self.imdata = imdata
        for i, values in kwargs.items():
            setattr(self, i, values)
        print(f'\- added ImagingData module. consisting of {self.imdata.shape} ROIs x Frames.')


    @property
    def n_frames(self):
        """number of imaging frames in imaging cellsdata"""
        return self.imdata.shape[1]

    @property
    def n_rois(self):
        """number of ROIs in imaging cellsdata"""
        return self.imdata.shape[0]


class ImagingMetadata:
    """Metadata about imaging system parameters."""

    def __init__(self, microscope, n_frames, fps, frame_x, frame_y, n_planes, pix_sz_x, pix_sz_y, **kwargs):
        self.microscope = microscope  #: given name of microscope
        self.n_frames = n_frames  #: number of imaging frames in the current trial
        self.fps = fps  #: rate of imaging acquisition (frames per second)
        self.frame_x = frame_x  #: num of pixels in the x direction of a single frame
        self.frame_y = frame_y  #: num of pixels in the y direction of a single frame
        self.n_planes = n_planes  #: num of FOV planes in imaging acquisition
        self.pix_sz_x = pix_sz_x  #: size of a single imaging pixel in x direction (microns)
        self.pix_sz_y = pix_sz_y  #: size of a single imaging pixel in y direction (microns)
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return f'ImagingMetadata for imaging cellsdata collected with {self.microscope}.'


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
    Suite2p = None

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
        from packerlabimaging.plotting.plotting import SingleFrame
        stack = SingleFrame(tiff_path=self.tiff_path, frame_num=frame_num, title=title)

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
        from packerlabimaging.utils.images import SaveDownsampledTiff
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
