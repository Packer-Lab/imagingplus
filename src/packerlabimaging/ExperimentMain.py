### MAIN SCRIPT
"""
This is the main collection of functions for processing and analysis of calcium imaging data in the Packer lab.

The fundamental object of data is one t-series of imaging, along with its associated .Paq file, accessory files generated
from the microscope during data collection, and any user generated files associated with the experiment of this t-series.

"""

# option for pre-loading s2p results, and loading Experiment with some trials included in s2p results and some not. -- use s2p_use arg in TrialsInformation dict
## TODO need to figure out where to set option for providing input for s2p_batch_size if pre-loading s2p results
## TODO consider providing arg values for all optical experiment analysis hyperparameters

from __future__ import absolute_import
from dataclasses import dataclass
from typing import Optional, MutableMapping
from .utils.classes import TrialsInformation

import os
import time
import re

import pickle

from suite2p.run_s2p import run_s2p

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)

from .utils import io

# UTILITIES

# dictionary of terms, phrases, etc. that are used in the processing and analysis of imaging data
terms_dictionary = {
    'dFF': "normalization of datatrace for a given imaging ROI by subtraction and division of a given baseline value",
    'ROI': "a single ROI from the imaging data"
}


def define_term(x: str):
    try:
        print(f"{x}:\t{terms_dictionary[x]}") if type(x) is str else print(
            'ERROR: please provide a string object as the key')
    except KeyError:
        print(f'input - {x} - not found in dictionary')


# TEMP VARIABLES FOR DEVELOPMENT USAGES
N_PLANES = 1
NEUROPIL_COEFF = 0.7


# CLASS DEFINITIONS



@dataclass
class Experiment:
    """A class to initialize and store data of an imaging experiment. This class acts as a bucket to contain
    information about individual trial objects. """
    date: str
    comments: str
    dataPath: str
    analysisSavePath: str  # main dir where the experiment object and the trial objects will be saved to
    microscope: str
    expID: str
    TrialsInformation: MutableMapping[str, TrialsInformation]
    useSuite2p: bool = False
    s2pResultsPath: Optional[str] = None  ## path to the parent directory containing the ops.npy file

    def __post_init__(self):

        trial_descr = f""
        for trial in self.trialIDs:
            trial_descr += f"\n\t{trial}, {self.TrialsInformation[trial]['trialType']}, {self.TrialsInformation[trial]['expGroup']}"
        print(f'***********************')
        print(f'CREATING new Experiment: (date: {self.date}, expID: {self.expID}), with trials: {trial_descr}')
        print(f'***********************\n\n')

        # need to check that the required keys are provided in TrialsInformation
        for i in ['trialType', 'tiff_path', 'expGroup']:
            for trial in [*self.TrialsInformation]:
                assert i in [*self.TrialsInformation[
                    trial]], f"must provide {i} as a key in the TrialsInformation dictionary for trial: {trial}"

        # START PROCESSING FOR EXPERIMENT:

        # 1) start suite2p:
        self._trialsSuite2p = []
        if self.useSuite2p or self.s2pResultsPath: self._add_suite2p_experiment()

        # 2) create individual trial objects
        self._runExpTrialsProcessing()

        # 3) save Experiment object:
        self._get_save_location()
        os.makedirs(self.analysisSavePath, exist_ok=True)
        self.save_pkl(pkl_path=self.pkl_path)

        # FINISH
        print(f'\n\n\n******************************')
        print(f"NEW Experiment created: \n")
        print(self.__str__())

    def _get_save_location(self):
        if self.analysisSavePath[-4:] == '.pkl':
            self.__pkl_path = self.analysisSavePath
            self.analysisSavePath = self.analysisSavePath[
                                    :[(s.start(), s.end()) for s in re.finditer('/', self.analysisSavePath)][-1][0]]
        else:
            self.analysisSavePath = self.analysisSavePath + '/' if self.analysisSavePath[
                                                                       -1] != '/' else self.analysisSavePath
            self.__pkl_path = f"{self.analysisSavePath}{self.expID}_analysis.pkl"
        os.makedirs(self.analysisSavePath, exist_ok=True)

    def _get_trial_infor(self, trialID: str):
        return f"\n\t{trialID}: {self.TrialsInformation[trialID]['trialType']}, {self.TrialsInformation[trialID]['expGroup']}"

    def __repr__(self):
        return f"packerlabimaging.Experiment object (date: {self.date}, expID: {self.expID})"

    def __str__(self):
        lastsaved = time.ctime(os.path.getmtime(self.pkl_path))
        __return_information = f"packerlabimaging Experiment object (last saved: {lastsaved}), date: {self.date}, expID: {self.expID}, microscope: {self.microscope}"
        __return_information = __return_information + f"\nfile path: {self.pkl_path}"
        __return_information = __return_information + f"\n\ntrials in Experiment object:"
        for trial in self.TrialsInformation:
            # print(f"\t{trial}: {self.TrialsInformation[trial]['trialType']} {self.TrialsInformation[trial]['expGroup']}")
            __return_information = __return_information + self._get_trial_infor(trialID=trial)
        return f"{__return_information}\n"

    def _add_suite2p_experiment(self):
        from .processing import suite2p
        from .processing.suite2p import Suite2pResultsExperiment

        for trial in self.trialIDs:
            if trial not in self._trialsSuite2p:
                assert 's2p_use' in [*self.TrialsInformation[
                    trial]], 'when trying to utilize suite2p , must provide value for `s2p_use` ' \
                             'in TrialsInformation[trial] for each trial to specify if to use trial for this suite2p associated with this experiment'
                self._trialsSuite2p.append(trial) if self.TrialsInformation[trial]['s2p_use'] else None

        if self.s2pResultsPath:  # if s2pResultsPath provided then try to find and pre-load results from provided s2pResultsPath, raise error if cannot find results
            # search for suite2p results items in self.suite2pPath, and auto-assign s2pRunComplete --> True if found successfully
            __suite2p_path_files = os.listdir(self.s2pResultsPath)
            self._s2pResultExists = False
            for filepath in __suite2p_path_files:
                if 'ops.npy' in filepath:
                    self._s2pResultExists = True
                    break
            if self._s2pResultExists:
                self.Suite2p = suite2p.Suite2pResultsExperiment(s2pResultsPath=self.s2pResultsPath,
                                                                trialsSuite2p=self._trialsSuite2p)
            else:
                raise ValueError(
                    f"suite2p results could not be found. `suite2pPath` provided was: {self.s2pResultsPath}")
        elif self.useSuite2p:  # no s2pResultsPath provided, so initialize without pre-loading any results
            self._s2pResultExists = False
            self._suite2p_save_path = self.analysisSavePath + '/suite2p/'
            self.Suite2p = Suite2pResultsExperiment(trialsSuite2p=self._trialsSuite2p)

    def update_suite2p(self, trialID: str = None, s2pResultsPath: str = None):  # TODO need to review structure and usage of this func
        """
        Add trial to suite2p trial list, and/or update an existing instance of the Suite2p submodule in Experiment.

        :param trialID: trial name to add to suite2p list
        :param s2pResultsPath: path to output of suite2p to use to update .Suite2p submodule
        """
        from .processing.suite2p import Suite2pResultsExperiment

        assert self.Suite2p is not None, 'No existing Suite2p sub-module found in Experiment.'

        if trialID is not None and trialID not in self._trialsSuite2p:
            self._trialsSuite2p.append(trialID)
            self.TrialsInformation[trialID]['s2p_use'] = True
        try:
            if s2pResultsPath:  # if s2pResultsPath provided then try to find and pre-load results from provided s2pResultsPath, raise error if cannot find results
                # search for suite2p results items in self.suite2pPath, and auto-assign s2pRunComplete --> True if found successfully
                __suite2p_path_files = os.listdir(self.s2pResultsPath)
                _s2pResultExists = False
                for filepath in __suite2p_path_files:
                    if 'ops.npy' in filepath:
                        _s2pResultExists = True
                        break
                if _s2pResultExists:
                    self.Suite2p = Suite2pResultsExperiment(s2pResultsPath=s2pResultsPath,
                                                            trialsSuite2p=self._trialsSuite2p)
                    self.s2pResultsPath = s2pResultsPath
                else:
                    raise ValueError(
                        f"suite2p results could not be found. `suite2pPath` provided was: {self.s2pResultsPath}")
            elif self.useSuite2p:  # no s2pResultsPath provided, so initialize without pre-loading any results
                self._s2pResultExists = False
                self._suite2p_save_path = self.analysisSavePath + '/suite2p/'
                self.Suite2p = Suite2pResultsExperiment(trialsSuite2p=self._trialsSuite2p)
        except Exception:
            raise ValueError(f"something went wrong. could not update suite2p from: {s2pResultsPath}.")

    def add_trial(self, trial_id: str, total_frames_stitched: int, trials_information: TrialsInformation = None):  # TODO need to figure out if there's a way around providing total_frames_stitched as arg! also add example of .add_trial() to tutorial!
        """
        Create and add trial object to the experiment. Return trial object.

        :param trial_id:
        :param total_frames_stitched:
        :param trials_information:
        :return: trial object

        """
        from .TwoPhotonImagingMain import TwoPhotonImagingTrial
        from .AllOpticalMain import AllOpticalTrial
        from .onePstimMain import OnePhotonStim

        print(f"\n\n\- ADDING trial: {trial_id}, expID: ({self.expID})")
        __trialsInformation = trials_information if trials_information else self.TrialsInformation[trial_id]
        _metainfo = {
            'exp_id': self.expID,
            'trial_id': trial_id,
            'date': self.date,
            't_series_id': f"{self.expID} {trial_id}",  # TODO consider removing occurrence especially since its redundant as a property in twophoton imaging
            'TrialsInformation': __trialsInformation
        }

        if _metainfo['TrialsInformation']['trialType'] == 'TwoPhotonImagingTrial':
            if self.TrialsInformation[trial_id]['s2p_use']:  # TODO could use switch statements in the 2p imaging trial class...
                trial_obj = TwoPhotonImagingTrial(metainfo=_metainfo, analysis_save_path=self.analysisSavePath,
                                                  paq_options=_metainfo['TrialsInformation']['PaqInfoTrial'],
                                                  microscope=self.microscope,
                                                  total_frames_stitched=total_frames_stitched,
                                                  suite2p_experiment_obj=self.Suite2p)
            else:
                trial_obj = TwoPhotonImagingTrial(metainfo=_metainfo, analysis_save_path=self.analysisSavePath,
                                                  microscope=self.microscope,
                                                  paq_options=_metainfo['TrialsInformation']['PaqInfoTrial'])

        elif _metainfo['TrialsInformation']['trialType'] == 'AllOpticalTrial':
            if self.TrialsInformation[trial_id]['s2p_use']:
                trial_obj = AllOpticalTrial(metainfo=_metainfo,
                                            naparm_path=_metainfo['TrialsInformation']['naparm_path'],
                                            analysis_save_path=self.analysisSavePath, microscope=self.microscope,
                                            paq_options=_metainfo['TrialsInformation']['PaqInfoTrial'], prestim_sec=1.0,
                                            poststim_sec=3.0, pre_stim_response_window=0.500,
                                            post_stim_response_window=0.500,
                                            total_frames_stitched=total_frames_stitched,
                                            suite2p_experiment_obj=self.Suite2p)
            else:
                trial_obj = AllOpticalTrial(metainfo=_metainfo,
                                            naparm_path=_metainfo['TrialsInformation']['naparm_path'],
                                            analysis_save_path=self.analysisSavePath, microscope=self.microscope,
                                            paq_options=_metainfo['TrialsInformation']['PaqInfoTrial'], prestim_sec=1.0,
                                            poststim_sec=3.0, pre_stim_response_window=0.500,
                                            post_stim_response_window=0.500)

        elif _metainfo['TrialsInformation']['trialType'] == 'OnePhotonAllOpticalTrial':
            trial_obj = OnePhotonStim(metainfo=_metainfo, analysis_save_path=self.analysisSavePath,
                                      microscope=self.microscope)
        else:
            raise NotImplementedError(
                f"the trial type ({_metainfo['TrialsInformation']['trialType']}) for this trial is not implemented yet."
                f"Compatible trialType options are: 'TwoPhotonImagingTrial', 'AllOpticalTrial', 'OnePhotonStim' ")

        # update self.TrialsInformation using the information from new trial_obj
        self.TrialsInformation[trial_id]['analysis_object_information'] = {'series ID': trial_obj.t_series_name,
                                                                          'repr': trial_obj.__repr__(),
                                                                          'pkl path': trial_obj.pkl_path}

        # initialize suite2p for trial objects
        if 's2p_use' in [*__trialsInformation] and __trialsInformation['s2p_use'] is True:
            self.update_suite2p(trialID=trial_id)

        return trial_obj

    def _runExpTrialsProcessing(self):
        """Runs processing of individual Trials that will contribute to the overall Experiment.

        Processing of individual Trials is carried out based on the contents of self.TrialsInformation.
        """
        total_frames_stitched = 0  # used in calculating # of frames from a single trial in the overall suite2p run
        for trial in self.trialIDs:
            trial_obj = self.add_trial(trial_id=trial, total_frames_stitched=total_frames_stitched)
            total_frames_stitched += trial_obj.imparams.n_frames

    @property
    def trialIDs(self):
        return [*self.TrialsInformation]

    @property
    def tiff_path_dir(self):
        "retrieve the parent directory of the provided tiff_path"
        __first_trial_in_experiment = [*self.TrialsInformation][0]
        __tiff_path_first_trial = self.TrialsInformation[__first_trial_in_experiment]['tiff_path']
        return __tiff_path_first_trial[:[(s.start(), s.end()) for s in re.finditer('/', __tiff_path_first_trial)][-1][
            0]]  # this is the directory where the Bruker xml files associated with the 2p imaging TIFF are located

    @property
    def pkl_path(self):
        "path in Analysis folder to save pkl object"
        return self.__pkl_path

    @pkl_path.setter
    def pkl_path(self, path: str):
        self.__pkl_path = path

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
        print("\t -- Experiment analysis object saved to %s -- " % self.pkl_path)

    def save(self):
        self.save_pkl()

    def load_trial(self, trialID: str):
        "method for importing individual trial objects from Experiment instance using the trial id for a given trial"
        try:
            trial_pkl_path = self.TrialsInformation[trialID]['analysis_object_information']['pkl path']
            trialobj = _io.import_obj(trial_pkl_path)
            return trialobj
        except KeyError:
            raise KeyError("trial_id not found in Experiment instance.")


class WideFieldImaging:
    """
    WideField imaging data object.

    """

    def __init__(self, tiff_path, paq_path, exp_metainfo, pkl_path):
        self.tiff_path = tiff_path
        self.paq_path = paq_path
        self.metainfo = exp_metainfo
        # set and create analysis save s2pResultsPath directory
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
        if 't_series_id' in [*self.metainfo]:
            return f"{self.metainfo['t_series_id']}"
        elif "exp_id" in [*self.metainfo] and "trial_id" in [*self.metainfo]:
            return f'{self.metainfo["exp_id"]} {self.metainfo["trial_id"]}'
        else:
            raise ValueError('no information found to retrieve t series id')

    def save_pkl(self, pkl_path: str = None):
        if pkl_path is None:
            if hasattr(self, 'pkl_path'):
                pkl_path = self.pkl_path
            else:
                raise ValueError(
                    'pkl s2pResultsPath for saving was not found in data object attributes, please provide pkl_path to save to')
        else:
            self.pkl_path = pkl_path

        with open(self.pkl_path, 'wb') as f:
            pickle.dump(self, f)
        print("\n\t -- data object saved to %s -- " % pkl_path)

    def save(self):
        self.save_pkl()
