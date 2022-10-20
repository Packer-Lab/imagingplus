# this file contains the two fundamental class types (Trial and Experiment) needed to construct an experiment in imagingplus
from __future__ import absolute_import
from dataclasses import dataclass, field
from typing import Optional, MutableMapping, Union, TypedDict, List, Dict, Any
import numpy as np
import pandas as pd

from imagingplus.main.subcore import ImagingMetadata, ImagingData, CellAnnotations, TemporalData
from imagingplus.processing.anndata import AnnotatedData

# from imagingplus.processing.imagingMetadata import ImagingMetadata
from imagingplus.processing.denoising import Deepinterpolation
from imagingplus.processing.suite2p import Suite2pResultsTrial
from imagingplus.utils.classes import UnavailableOptionError, TrialMetainfo

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
from imagingplus.utils.utils import findClosest
from imagingplus.utils.images import ImportTiff, SaveDownsampledTiff

"""

NOTES:
- need to figure out how to handle temporal data that is collected at a different sampling rate!!!!
    - could have cases where you have different hardware collecitng temporal data!

- merging of multiple ImagingTrials + TemporalData across the same CellAnnotations into a single ImagingTrial

"""


@dataclass
class SingleImage:
    dataPath: str                       #: path to source image
    date: str = None                    #: date associated with processed image
    imgID: str = None                   #: id of image
    expGroup: str = None                #: experimental group of image
    comment: str = None                 #: notes regarding image
    imparams: ImagingMetadata = None    #: image collection parameters of image

    def __post_init__(self):
        # self.data = tf.imread(self.dataPath)
        self.data = ImportTiff(self.dataPath)

    def plotImg(self, **kwargs):
        """
        Plot single scan image data.

        :param kwargs: refer to imagingplus.plotting.plotting.SingleFrame
        """
        from imagingplus.plotting.plotting import SingleFrame
        SingleFrame(imstack=self.data, **kwargs)


@dataclass
class ImagingTrial:
    dataPath: str       #: path to tiff file containing data to create Imaging Trial
    saveDir: str        #: main dir where the experiment object and the trial objects will be saved to
    date: str           #: date to assign to Imaging Trial
    trialID: str        #: trial ID of Imaging Trial
    expID: str          #: experiment ID associated with current Imaging Trial
    expGroup: str = ''  #: group assignment within experiment
    comment: str = ''   #: notes related to Imaging Trial
    imparams: ImagingMetadata = None        #: ImagingMetadata related to current Imaging Trial
    data: AnnotatedData = None              #: anndata structure for current Imaging Trial
    imdata: ImagingData = None              #: Imaging data for current Imaging Trial
    cells: CellAnnotations = None           #: Cells Annotations data for current Imaging Trial
    tmdata: TemporalData = None             #: Temporal Data for current Imaging Trial
    Suite2p: Suite2pResultsTrial = None     #: Suite2p Results for current Imaging Trial
    Deepinterpolation: Deepinterpolation = None     #: Deepinterpolation (denoising) for current Imaging Trial
    optical_channel_name: str = 'Imaging Channel'   #: a name of the optical channel collected for current trial (default = Imaging Channel)
    imaging_plane_name: str = 'Plane 0'   #: a name of the imaging plane collected for current trial (default = Plane 0)

    def __post_init__(self):
        self.metainfo = TrialMetainfo(date=self.date, trialID=self.trialID, expID=self.expID, expGroup=self.expGroup,
                                      comment=self.comment, paths={})

        if os.path.exists(self.dataPath):
            self.metainfo['paths']['dataPath'] = self.dataPath
        else:
            raise FileNotFoundError(f"dataPath does not exist: {self.dataPath}")

        # processing collect mean FOV Trace -- after collecting imaging params and Paq timing info
        im_stack = ImportTiff(tiff_path=self.data_path)
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
        return repr(f"{self.t_series_name} (ImagingTrial)")

    # todo maybe ? - add alternative constructor to handle construction if temporal cellsdata or cell annotations or imaging cellsdata is provided
    # might be useful for providing like submodules (especially things like Suite2p)
    @classmethod
    def newImagingTrialfromExperiment(cls, experiment, trialID, dataPath, date, comment=''):
        """
        Creates a new ImagingTrial under an Experiment.

        :param experiment: Experiment object
        :param trialID: ID of new Imaging Trial to be added
        :param dataPath: path to imaging data tiff file of new Imaging Trial
        :param date: associated date of new Imaging Trial
        :param comment: associated notes for new Imaging Trial

        :return: ImagingTrial instance
        """

        trialobj: ImagingTrial = cls(dataPath=dataPath,
                                     saveDir=experiment.saveDir,
                                     date=date,
                                     trialID=trialID,
                                     expID=experiment.expID,
                                     comment=comment)
        experiment.add_imaging_trial(trialID=trialID, trialobj=trialobj)
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

    def frameNum(self, time: float) -> int:
        """
        Calculate and return frame number of time point, uses the temporally captured frame timing signals to

        :param time: Time point to retrieve frame number from imaging series.
        :return: Frame number at input time
        """
        # subtract the timestamp of the desired frame from the first frame clock timestamp, and convert to secs.
        time_stamp = time * self.tmdata.sampling_rate
        frame_num = findClosest(self.tmdata.frame_times, self.tmdata.frame_times[0] + time_stamp)[1]
        return int(frame_num)

    def timePoint(self, frame):
        """
        Calculate and return timepoint of input frame, uses the temporally captured frame timing signals.

        :param frame: Input frame to retrieve timepoint from imaging series time signals
        :return: Timepoint
        """
        # subtract the timestamp of the desired frame from the first frame clock timestamp, and convert to secs.
        assert hasattr(self.tmdata,
                       'frame_times'), 'cannot retrieve frame timepoint without `frame_times` defined in .tmdata'
        assert frame in range(0, len(self.tmdata.frame_times)), 'frame not in indexes of frame times'
        frame_time = (self.tmdata.frame_times[frame] - self.tmdata.frame_times[0]) / self.tmdata.sampling_rate
        return np.round(frame_time, 5)

    @property
    def data_path(self):
        """data path of current trial object"""
        return self.metainfo['paths']['dataPath']

    @property
    def t_series_name(self):
        """
        Return convience identification label of trial ID and exp ID.
        :return:
        """
        if "expID" in self.metainfo.keys() and "trialID" in self.metainfo.keys():
            return f'{self.metainfo["expID"]} {self.metainfo["trialID"]}'
        else:
            raise ValueError('missing `expID` or `trialID` required to retrieve t series id')

    @property
    def data_path_dir(self):
        """Parent directory of data path"""
        return os.path.dirname(self.data_path)

    @property
    def pkl_path(self):
        """path in Analysis folder to save pkl object of current Imaging Trial"""
        return self.metainfo['paths']['pkl_path']

    @pkl_path.setter
    def pkl_path(self, path: str):
        """
        Setter for pkl path for current Imaging Trial.
        :param path:
        """
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
        from imagingplus.plotting.plotting import SingleFrame
        stack = SingleFrame(tiff_path=self.data_path, frame_num=frame_num, title=title)

        # stack = tf.imread(self.tiff_path, key=frame_num)
        # plt.imshow(stack, cmap='gray')
        # plt.suptitle(title) if title is not None else plt.suptitle(f'frame num: {frame_num}')
        # plt.show()
        return stack

    def addImagingMetadata(self, imaging_metadata: ImagingMetadata = None, Bruker: bool = False):
        """Add the imaging metadata submodule to the current ImagingTrial.

        :param Bruker: if True, then set ImagingMetadata as PrairieViewMetadata
        :param imaging_metadata: input ImagingMetadata to add to current ImagingTrial
        """
        # imaging metadata
        if Bruker:
            from imagingplus.processing.imagingMetadata import PrairieViewMetadata
            self.imparams = PrairieViewMetadata(pv_xml_dir=self.data_path_dir)
        elif imaging_metadata:
            self.imparams = imaging_metadata
        else:
            Warning(
                f"NO imaging microscope parameters set. follow imagingMetadata to create a custom imagingMicroscopeMetadata class.")

    ## below properties/methods have pre-requisite processing steps
    @property
    def n_frames(self):
        """number of imaging frames from tiff cellsdata.
        :return:
        """
        return self._n_frames

    def importTrialTiff(self) -> np.ndarray:
        """
        Import current trial's imaging tiff in full.

        :return: imaging tiff as numpy array
        """
        print(f"\n\- loading raw TIFF file from: {self.data_path}", end='\r')
        im_stack = ImportTiff(tiff_path=self.data_path)
        # im_stack = tf.imread(self.tiff_path)
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
        SaveDownsampledTiff(stack=stack, save_as=f"{self.saveDir}/{self.date}_{self.trialID}_downsampled.tif")

    def create_anndata(self, imdata_type: str, imdata: ImagingData=None, cells: CellAnnotations=None, tmdata: TemporalData=None):
        """
        Creates annotated cellsdata (see anndata library for more information on AnnotatedData) object based around the Ca2+ matrix of the imaging trial.

        Note: the tmdata is used as the source of truth for imaging frames time signals. If there is a mismatch, then tmdata timestamps for extraneous imaging frames will not be included in adata table.

        :param imdata_type: label to describe primary dataset in anndata structure.
        :return:

        """
        # if not self.imdata and self.cells and self.tmdata:
        #     raise ValueError(
        #         'cannot create anndata table. anndata creation only available if experiments have ImagingData (.imdata), CellAnnotations (.cells) and TemporalData (.tmdata)')

        imdata = self.imdata if imdata is None else imdata
        assert hasattr(imdata, 'imdata') or tmdata != None, 'no `imdata` under .imdata'

        # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
        obs = self.cells if cells is None else cells
        assert hasattr(obs, 'cellsdata'), 'no `data` under .cells'
        assert obs.cellsdata.shape[0] == imdata.imdata.shape[0], 'incorrect # of cells passed to anndata creation.'

        # SETUP THE VARIABLES ANNOTATIONS TO USE IN anndata
        vars = self.tmdata.data if tmdata is None else tmdata
        if vars.shape[0] != imdata.imdata.shape[1]:
            print(f'WARNING: mismatch of timestamps in tmdata ({vars.shape[0]}) and frames in imdata ({imdata.n_frames}) passed to anndata creation.')
            vars = vars[:imdata.imdata.shape[1]]

        print(f"\n\----- CREATING annotated cellsdata object using AnnData:")

        anndata_setup = {'X': imdata.imdata, 'data_label': imdata_type, 'obs': obs.cellsdata, 'var': vars,
                         'obsm': obs.multidim_data if obs.multidim_data else None}

        adata = AnnotatedData(**anndata_setup)

        print(f"\n{adata}")
        return adata

    @staticmethod
    def collectSignalFromCoords(tiff_paths: Union[tuple, list], target_coords_masks: np.ndarray):
        """
        Collect and return mean signal from target mask areas for each cell from provided tiffs.

        :param tiff_paths: tiff paths to load movies to use for collecting signals
        :param target_coords_masks: array of masks per cell to use for collecting signals from loaded tiffs

        """

        num_coords = len(target_coords_masks)
        target_areas = target_coords_masks

        targets_trace_full = np.zeros([num_coords, 0], dtype='float32')

        for tiff_path in tiff_paths:
            print(f'|- reading tiff: {tiff_path}')

            from imagingplus.utils.images import ImportTiff
            data = ImportTiff(tiff_path)
            print(f'\t\- shape: {data.shape}')

            targets_trace = np.zeros([num_coords, data.shape[0]], dtype='float32')
            for coord in range(num_coords):
                x = data[:, target_areas[coord, 1], target_areas[coord, 0]]
                targets_trace[coord] = np.mean(x, axis=1)

            targets_trace_full = np.concatenate((targets_trace_full, targets_trace), axis=1)

        return targets_trace_full

    def run_deepinterpolation(self, input_file: str, output_file: str, model_path: str, generator_param: dict = {},
                              inferrence_param: dict = {}):
        self.Deepinterpolation = Deepinterpolation.deepinterpolate(input_file, output_file, model_path, generator_param,
                                                                   inferrence_param)


# noinspection DuplicatedCode
@dataclass
class Experiment:
    """Overall experiment. It can collate all imaging trials that are part of a single field-of-view (FOV) of imaging"""

    expID: str  #: given identification name for experiment
    dataPath: str  #: main dir where the imaging data is contained
    saveDir: str  #: main dir where the experiment object and the trial objects will be saved to
    comment: str = ''  #: notes related to experiment
    singleImages: Dict[str, SingleImage] = None  #: contains single image frames from each experiment
    experimenter = ''  #: experimenter performing pertinent data collection or data analysis
    lab = ''  #: research group with with experiments were performed
    institution = ''  #: insititution where experiments were performed

    def __post_init__(self):
        print(f'***********************')
        print(f'CREATING new Experiment: (expID: {self.expID})')
        print(f'***********************\n\n')

        self.singleImages = {}

        # self.TrialsInformation: MutableMapping[str, Union[str, TrialsInformation, PaqInfo]] = {}  #: dictionary of metadata information about each trial. Gets filled while adding each trial.
        self.TrialsInformation: MutableMapping[
            str, TrialMetainfo] = {}  #: dictionary of metadata information about each trial. Gets filled while adding each trial.

        # suite2p related attrs initialization
        from imagingplus.processing.suite2p import Suite2pExperiment
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
                __return_information = __return_information + self.get_trial_info(trialID=trial)
            return f"{__return_information}\n"
        else:
            return f"{__return_information}\n"

    def _get_save_location(self):
        """
        Specify location to create pickled experiment object.

        :return:
        """
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

    def get_trial_info(self, trialID: str):
        """
        Retrieve meta-information stored about trial (trialID) in the current experiment.

        :param trialID: trial to retrieve meta-information
        :return:
        """
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

    def add_imaging_trial(self, trialID: str = None, trialobj: ImagingTrial = None, **kwargs):
        """
        Add trial to current experiment by adding trial meta-information.

        :param trialID: name ID of trialobj
        :param trialobj: ImagingTrial object.
        """

        print(f"\n\n\- ADDING trial: {trialobj.trialID} to {self.expID} experiment", end='\r')

        if trialobj is None:
            trialobj = ImagingTrial.newImagingTrialfromExperiment(experiment=self, trialID=trialID, **kwargs)

        self.TrialsInformation[trialobj.trialID] = trialobj.metainfo
        self.TrialsInformation[trialobj.trialID]['expGroup'] = trialobj.metainfo['expGroup']
        for attr in dir(trialobj):
            if 'path' in attr:  # add any attr's containing 'path' from trialobj to the TrialsInformation dict
                self.TrialsInformation[trialobj.trialID]['paths'][attr] = getattr(trialobj, attr)

        self.save()
        print(f"|- ADDED trial: {trialobj.trialID} to {self.expID} experiment")

    def add_single(self, singleimg: SingleImage):
        """
        Add SingleImage to current Experiment (e.g. for adding a high quality single scan of current Experiment's FOV.

        :param singleimg: SingleImage object
        """
        setattr(self, singleimg.imgID, singleimg)

    def combine_trials(self):
        """todo: Combine anndata table of trials with same cells.
        """

    @property
    def trialIDs(self):
        """
        Return list of trial IDs contained in Experiment.
        """
        return [*self.TrialsInformation]

    @property
    def pkl_path(self):
        """path in Analysis folder to save pkl object
        :return:
        """
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
        """
        Shortcut for saving current Experiment to its designated pkl path.
        """
        self.save_pkl()

    def load_trial(self, trialID: str):
        """
        Import individual trial objects from Experiment instance using the trial id for a given trial

        :param trialID: ID of trial to import.
        """
        try:
            trial_pkl_path = self.TrialsInformation[trialID]['paths']['pkl_path']
            from imagingplus import import_obj
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
        from imagingplus.processing.suite2p import Suite2pExperiment, Suite2pResultsTrial

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
            __s2pResultExists = False
            for filepath in __suite2p_path_files:
                if 'ops.npy' in filepath:
                    __s2pResultExists = True
                    break
            if __s2pResultExists:
                # self._suite2p_save_path = s2pResultsPath
                self.Suite2p = Suite2pExperiment(trialsTiffsSuite2p=_trialsTiffsSuite2p,
                                                 s2pResultsPath=s2pResultsPath)
            else:
                raise ValueError(
                    f"suite2p results could not be found. `suite2pPath` provided was: {s2pResultsPath}")
        else:  # no s2pResultsPath provided, so initialize without pre-loading any results
            self.Suite2p = Suite2pExperiment(trialsTiffsSuite2p=_trialsTiffsSuite2p)

        # print(f"SHAPE OF LOADED SUITE2P DATA: {self.Suite2p.raw.shape}")

        # print(self.Suite2p)
        # adding of suite2p trial level as well in this function as well
        total_frames = 0  # TODO need to use approach of cross-checking trial ID with trials run in 2p imaging outputs to select correct frame numbers relative to the s2p output
        for trial in self.Suite2p.trials:
            trialobj = self.load_trial(trialID=trial)
            trialobj.Suite2p = Suite2pResultsTrial(s2pExp=self.Suite2p, trial_frames=(
                total_frames, total_frames + trialobj.n_frames))  # use trial obj's current trial key_frames
            assert trialobj.n_frames == (trialobj.Suite2p.imdata.shape[
                1]), f'mismatch of trial frames lengths in loading Suite2p results for trial: {trialobj.trialID}'
            total_frames += trialobj.Suite2p.imdata.shape[1]
            # print(f'frame count: ', total_frames)
            trialobj.save()

        print(f'|- Finished adding suite2p module to experiment and trials. Located under .Suite2p')
        self.save()

    @property
    def suite2p_save_path(self):
        """
        Save path for suite2p outputs/results.
        """
        return self._suite2p_save_path
