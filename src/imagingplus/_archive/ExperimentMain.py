### MAIN SCRIPT
"""
This is the main collection of functions for processing and analysis of calcium imaging cellsdata in the Packer lab.

The fundamental object of cellsdata is one t-series of imaging, along with its associated .Paq file, accessory files generated
from the microscope during cellsdata collection, and any user generated files associated with the experiment of this t-series.

"""

# option for pre-loading s2p results, and loading Experiment with some trials included in s2p results and some not. -- use s2p_use arg in TrialsInformation dict

from __future__ import absolute_import
from dataclasses import dataclass
from typing import Optional, MutableMapping, Union

from packerlabimaging.utils.classes import TrialsInformation, PaqInfo

import os
import time
import re

import pickle

from suite2p.run_s2p import run_s2p
from packerlabimaging.utils import io



# TEMP VARIABLES FOR DEVELOPMENT
N_PLANES = 1
NEUROPIL_COEFF = 0.7


@dataclass
class Experiment:
    """A class to initialize and store cellsdata of an imaging experiment. This class acts as a bucket to contain
    information about individual trial objects. """
    date: str  #: date of experiment cellsdata collection
    dataPath: str   #: main dir where the imaging cellsdata is contained
    saveDir: str  #: main dir where the experiment object and the trial objects will be saved to
    expID: str  #: given identification name for experiment
    comment: str = ''  #: notes related to experiment

    def __post_init__(self):
        print(f'***********************')
        print(f'CREATING new Experiment: (date: {self.date}, expID: {self.expID})')
        print(f'***********************\n\n')

        self.TrialsInformation: MutableMapping[str, Union[TrialsInformation, PaqInfo]] = {}  #: dictionary of metadata information about each trial

        # suite2p related attrs initialization
        self._trialsTiffsSuite2p = {}  #: dictionary of trial IDs and their respective .tiff paths for each trial that will be used in Suite2p processing for current experiment
        self._s2pResultExists = False  #: flag for whether suite2p results exist for current experiment
        self._suite2p_save_path = self.saveDir + '/suite2p/'  #: default location to save Suite2p output results of current experiment
        self.Suite2p = None  #: suite2p submodule

        # save Experiment object
        self._get_save_location()
        os.makedirs(self.saveDir, exist_ok=True)
        self.save_pkl(pkl_path=self.pkl_path)


        print(f'\n\n\n******************************')
        print(f"NEW Experiment created: \n")
        print(self.__str__())

    def __repr__(self):
        return f"imagingplus.Experiment object (date: {self.date}, expID: {self.expID})"

    def __str__(self):
        lastsaved = time.ctime(os.path.getmtime(self.pkl_path))
        __return_information = f"imagingplus Experiment object (last saved: {lastsaved}), date: {self.date}, expID: {self.expID}"
        __return_information = __return_information + f"\nfile path: {self.pkl_path}"

        if len(self.TrialsInformation) > 0:
            __return_information = __return_information + f"\n\ntrials in Experiment object:"
            for trial in self.TrialsInformation:
                __return_information = __return_information + self.get_trial_infor(trialID=trial)
            return f"{__return_information}\n"
        else:
            return f"{__return_information}\n"

    def _get_save_location(self):
        if self.saveDir[-4:] == '.pkl':
            self._pkl_path = self.saveDir
            self.saveDir = self.saveDir[
                                    :[(s.start(), s.end()) for s in re.finditer('/', self.saveDir)][-1][0]]
        else:
            self.saveDir = self.saveDir + '/' if self.saveDir[
                                                                       -1] != '/' else self.saveDir
            self._pkl_path = f"{self.saveDir}{self.expID}_analysis.pkl"
        os.makedirs(self.saveDir, exist_ok=True)

    def get_trial_infor(self, trialID: str):
        if trialID in [*self.TrialsInformation]:
            return f"\n\t{trialID}: {self.TrialsInformation[trialID]['repr']}, {self.TrialsInformation[trialID]['expGroup']}"
        else:
            ValueError(f"{trialID} not found in Experiment.")

    def add_trial(self, trialobj):
        """
        Add trial object to the experiment.

        :param trialobj:
        :param trialID:
        :param total_frames_stitched:
        :param trials_information:
        :return: trial object

        """

        print(f"\n\n\- ADDING trial: {trialobj.t_series_name}", end='\r')

        self.TrialsInformation[trialobj.trialID] = {}
        for attr in dir(trialobj):
            if 'path' in attr:
                self.TrialsInformation[trialobj.trialID][attr] = getattr(trialobj, attr)
        for key, val in trialobj._metainfo.items():
            self.TrialsInformation[trialobj.trialID][key] = val

        self.TrialsInformation[trialobj.trialID]['tiff_path'] = trialobj.data_path  # this should be redundant but just keeping up until this is confirmed.
        # update self.TrialsInformation using the information from new trial_obj
        self.TrialsInformation[trialobj.trialID]['series ID'] = trialobj.t_series_name
        self.TrialsInformation[trialobj.trialID]['repr'] = trialobj.__repr__()

        print(f"|- ADDED trial: {trialobj.t_series_name}")

    @property
    def trialIDs(self):
        return [*self.TrialsInformation]

    @staticmethod
    def tiff_path_dir(tiff_path):
        """retrieve the parent directory of the provided tiff_path"""
        return tiff_path[:[(s.start(), s.end()) for s in re.finditer('/', tiff_path)][-1][
            0]]  # this is the directory where the Bruker xml files associated with the 2p imaging TIFF are located

    @property
    def pkl_path(self):
        """path in Analysis folder to save pkl object"""
        return self._pkl_path

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
            trial_pkl_path = self.TrialsInformation[trialID]['pkl_path']
            from packerlabimaging import import_obj
            trialobj = import_obj(trial_pkl_path)
            return trialobj
        except KeyError:
            raise KeyError(f"trialID: {trialID} not found in Experiment instance.")


    def add_suite2p(self, s2p_trials: Union[list, str] = 'all', s2pResultsPath: Optional[str] = None):
        """Wrapper for adding suite2p results to Experiment. Can only be run after adding trials to Experiment.

        :param s2p_trials: list of trials to use for Suite2p processing/analysis pipeline, default = 'all' to use all trials of Experiment
        :param s2pResultsPath: optional, if suite2p already run then here provide path to plane0 folder of existing suite2p results
        """

        print(f'\- Adding suite2p module to experiment. Located under .Suite2p')

        from .processing import suite2p
        from .processing.suite2p import Suite2pResultsExperiment

        assert len([*self.TrialsInformation]) > 0, 'need to add at least 1 trial to Experiment before adding Suite2p functionality.'

        if s2p_trials is 'all': s2p_trials = self.trialIDs
        assert type(s2p_trials) == list and len(s2p_trials) > 0, 'no s2p trials given in list form.'
        for trial in s2p_trials: self._trialsTiffsSuite2p[trial] = self.TrialsInformation[trial]['tiff_path']

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
                self.Suite2p = suite2p.Suite2pResultsExperiment(trialsTiffsSuite2p=self._trialsTiffsSuite2p,
                                                                s2pResultsPath=s2pResultsPath)
            else:
                raise ValueError(
                    f"suite2p results could not be found. `suite2pPath` provided was: {s2pResultsPath}")
        else:  # no s2pResultsPath provided, so initialize without pre-loading any results
            self.Suite2p = Suite2pResultsExperiment(trialsTiffsSuite2p=self._trialsTiffsSuite2p)

        # print(self.Suite2p)
        # adding of suite2p trial level as well in this function as well
        total_frames = 0
        for trial in s2p_trials:
            trialobj = self.load_trial(trialID=trial)
            # print(f'n_frames', self.Suite2p.n_frames)
            trialobj.Suite2p = suite2p.Suite2pResultsTrial(s2pExp=self.Suite2p, trial_frames=(total_frames, total_frames + trialobj.n_frames))  # use trial obj's current trial frames
            trialobj.save()
            total_frames += trialobj.n_frames


        print(f'|- Finished adding suite2p module to experiment and trials. Located under .Suite2p')
        self.save()

    @property
    def suite2p_save_path(self):
        return self._suite2p_save_path











