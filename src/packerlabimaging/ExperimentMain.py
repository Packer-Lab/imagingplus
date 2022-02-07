### MAIN SCRIPT
"""
This is the main collection of functions for processing and analysis of calcium imaging data in the Packer lab.

The fundamental object of data is one t-series of imaging, along with its associated .paq file, accessory files generated
from the microscope during data collection, and any user generated files associated with the experiment of this t-series.

"""

# option for pre-loading s2p results, and loading Experiment with some trials included in s2p results and some not. -- use s2p_use arg in trialsInformation dict
## TODO need to figure out where to set option for providing input for s2p_batch_size if pre-loading s2p results
## TODO consider providing arg values for all optical experiment analysis hyperparameters

from __future__ import absolute_import
from dataclasses import dataclass
from typing import TypedDict, Optional, MutableMapping

import os
import time
import re

import pickle

from suite2p.run_s2p import run_s2p

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)

from .TwoPhotonImagingMain import TwoPhotonImagingTrial
from .onePstimMain import OnePhotonStim
from .AllOpticalMain import AllOpticalTrial
from . import _io
from .processing import suite2p

## UTILITIES

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


## TEMP VARIABLES FOR DEVELOPMENT USAGES
N_PLANES = 1
NEUROPIL_COEFF = 0.7

## CLASS DEFINITIONS


class trialsInformation(TypedDict, total=False):
    trialType: str
    tiff_path: str
    expGroup: str
    paq_path: Optional[str]
    s2p_use: Optional[str]
    naparm_path: Optional[str]
    analysis_object_information: Optional[TypedDict('analysis_object_information',
                                                    {'series ID': str, 'repr': str, 'pkl path': str})] = None


@dataclass
class Experiment:
    date: str
    comments: str
    dataPath: str
    analysisSavePath: str  # main dir where the experiment object and the trial objects will be saved to
    microscope: str
    expID: str
    trialsInformation: MutableMapping[str, trialsInformation]
    useSuite2p: bool = False
    s2pResultsPath: Optional[str] = None  ## path to the parent directory containing the ops.npy file
    def __post_init__(self):

        trialDescr = f""
        for trial in self.trialIDs:
            trialDescr += f"\n\t{trial}, {self.trialsInformation[trial]['trialType']}, {self.trialsInformation[trial]['expGroup']}"
        print(f'***********************')
        print(f'CREATING new Experiment: (date: {self.date}, expID: {self.expID}), with trials: {trialDescr}')
        print(f'***********************\n\n')




        ## need to check that the required keys are provided in trialsInformation
        for i in ['trialType', 'tiff_path', 'expGroup']:
            for trial in [*self.trialsInformation]:
                assert i in [*self.trialsInformation[trial]], f"must provide {i} as a key in the trialsInformation dictionary for trial: {trial}"

        # start processing for experiment imaging:

        # start suite2p action:
        if self.useSuite2p or self.s2pResultsPath: self.add_suite2p()

        # create individual trial objects
        self._runExpTrialsProcessing()


        # save Experiment object:
        # set and/or create analysis save path directory
        self._get_save_location()
        os.makedirs(self.analysisSavePath, exist_ok=True)
        self.save_pkl(pkl_path=self.pkl_path)

        # Attributes:
        print(f'\n\n\n******************************')
        print(f"NEW Experiment created: \n")
        print(self.__str__())

    def _get_save_location(self):
        if self.analysisSavePath[-4:] == '.pkl':
            self.__pkl_path = analysis_save_path
            self.analysisSavePath = self.analysisSavePath[:[(s.start(), s.end()) for s in re.finditer('/', self.analysisSavePath)][-1][0]]
        else:
            self.analysisSavePath = self.analysisSavePath + '/' if self.analysisSavePath[-1] != '/' else self.analysisSavePath
            self.__pkl_path = f"{self.analysisSavePath}{self.expID}_analysis.pkl"
        os.makedirs(self.analysisSavePath, exist_ok=True)


    def _get_trial_infor(self, trialID: str):
        return f"\n\t{trialID}: {self.trialsInformation[trialID]['trialType']}, {self.trialsInformation[trialID]['expGroup']}"

    def __repr__(self):
        return f"packerlabimaging.Experiment object (date: {self.date}, expID: {self.expID})"

    def __str__(self):
        lastsaved = time.ctime(os.path.getmtime(self.pkl_path))
        __return_information = f"packerlabimaging Experiment object (last saved: {lastsaved}), date: {self.date}, expID: {self.expID}, microscope: {self.microscope}"
        __return_information = __return_information + f"\nfile path: {self.pkl_path}"
        __return_information = __return_information + f"\n\ntrials in Experiment object:"
        for trial in self.trialsInformation:
            # print(f"\t{trial}: {self.trialsInformation[trial]['trialType']} {self.trialsInformation[trial]['expGroup']}")
            __return_information = __return_information + self._get_trial_infor(trialID=trial)
        return f"{__return_information}\n"



    def add_suite2p(self):  # TODO use add_suite2p method in __post_init__ (or refactor that bit out and create new method for run Suite2p)
        self._trialsSuite2p = []
        for trial in self.trialIDs:
            assert 's2p_use' in [*self.trialsInformation[trial]], 'when trying to utilize suite2p , must provide value for `s2p_use` ' \
                         'in trialsInformation[trial] for each trial to specify if to use trial for this suite2p associated with this experiment'
            self._trialsSuite2p.append(trial) if self.trialsInformation[trial]['s2p_use'] else None

        if self.s2pResultsPath:  # if s2pResultsPath provided then try to find and pre-load results from provided path, raise error if cannot find results
            # search for suite2p results items in self.suite2pPath, and auto-assign s2pRunComplete --> True if found successfully
            __suite2p_path_files = os.listdir(self.s2pResultsPath)
            self._s2pResultExists = False
            for filepath in __suite2p_path_files:
                if 'ops.npy' in filepath:
                    self._s2pResultExists = True
                    break
            if self._s2pResultExists:
                self.Suite2p = suite2p.Suite2pResultsExperiment(s2pResultsPath=self.s2pResultsPath, trialsSuite2p = self._trialsSuite2p)
            else:
                raise ValueError(f"suite2p results could not be found. `suite2pPath` provided was: {self.suite2pPath}")
        elif self.useSuite2p:  # no s2pResultsPath provided, so initialize without pre-loading any results
            self._s2pResultExists = False
            self._suite2p_save_path = self.analysisSavePath + '/suite2p/'
            self.Suite2p = _suite2p.Suite2pResultsExperiment(trialsSuite2p = self._trialsSuite2p)



    def add_trial(self, trialID: str = None, trialsInformation: trialsInformation = None):  # TODO use the add_trial method in _runExpTrialsProcessing! also add example of .add_trial() to tutorial!

        print(f"\n\n\- ADDING trial: {trial}, expID: ({self.expID})")
        __trialsInformation = trialsInformation if trialsInformation else self.trialsInformation[trialID]
        _metainfo = {
            'animal prep.': self.expID,
            'trial': trialID,
            'date': self.date,
            't series id': f"{self.expID} {trialID}",
            'trialsInformation': self.trialsInformation[trialID]
        }

        if _metainfo['trialsInformation']['trialType'] == 'TwoPhotonImagingTrial':
            if self.trialsInformation[trialID]['s2p_use']:  # TODO could use switch statements in the 2p imaging trial class...
                trial_obj = TwoPhotonImagingTrial(metainfo=_metainfo, analysis_save_path=self.analysisSavePath,
                                                  microscope=self.microscope, total_frames_stitched=total_frames_stitched, suite2p_experiment_obj=self.Suite2p)
            else:
                trial_obj = TwoPhotonImagingTrial(metainfo=_metainfo, analysis_save_path=self.analysisSavePath, microscope=self.microscope)

        elif _metainfo['trialsInformation']['trialType'] == 'AllOpticalTrial':
            if self.trialsInformation[trialID]['s2p_use']:
                trial_obj = AllOpticalTrial(metainfo=_metainfo,
                                            naparm_path=_metainfo['trialsInformation']['naparm_path'],
                                            analysis_save_path=self.analysisSavePath, microscope=self.microscope,
                                            prestim_sec=1.0, poststim_sec=3.0, pre_stim_response_window=0.500,
                                            post_stim_response_window=0.500, total_frames_stitched=total_frames_stitched,
                                            suite2p_experiment_obj=self.Suite2p)
            else:
                trial_obj = AllOpticalTrial(metainfo=_metainfo,
                                            naparm_path=_metainfo['trialsInformation']['naparm_path'],
                                            analysis_save_path=self.analysisSavePath, microscope=self.microscope,
                                            prestim_sec=1.0, poststim_sec=3.0, pre_stim_response_window=0.500,
                                            post_stim_response_window=0.500)

        elif _metainfo['trialsInformation']['trialType'] == 'AllOpticalTrial':
            trial_obj = OnePhotonStim(metainfo=_metainfo, analysis_save_path=self.analysisSavePath, microscope=self.microscope)
        else:
            raise NotImplementedError(f"the trial type ({_metainfo['trialsInformation']['trialType']}) for this trial is not implemented yet."
                                      f"Compatible trialType options are: 'TwoPhotonImagingTrial', 'AllOpticalTrial', 'OnePhotonStim' ")


    def _runExpTrialsProcessing(self):
        """Runs processing of individual Trials that will contribute to the overall Experiment.

        Processing of individual Trials is carried out based on the contents of self.trialsInformation.
        """
        total_frames_stitched = 0  # used in calculating # of frames from a single trial in the overall suite2p run
        for trial in self.trialIDs:
            print(f"\n\n\- PROCESSING trial: {trial}, expID: ({self.expID}) {'*'*20}")
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

                if trial in self._trialsSuite2p:  # TODO could use switch statements in the 2p imaging trial class...
                    trial_obj = TwoPhotonImagingTrial(metainfo=_metainfo, analysis_save_path=self.analysisSavePath,
                                                      microscope=self.microscope, total_frames_stitched=total_frames_stitched, suite2p_experiment_obj=self.Suite2p)
                else:
                    trial_obj = TwoPhotonImagingTrial(metainfo=_metainfo, analysis_save_path=self.analysisSavePath, microscope=self.microscope)


                # update self.trialsInformation using the information from new trial_obj
                self.trialsInformation[trial]['analysis_object_information'] = {'series ID': trial_obj.t_series_name,
                                                                                'repr': trial_obj.__repr__(),
                                                                                'pkl path': trial_obj.pkl_path}

            # initialize AllOpticalTrial
            elif _metainfo['trialsInformation']['trialType'] == 'AllOpticalTrial':
                if 'tiff_path' not in [*_metainfo['trialsInformation']] or 'paq_path' not in [*_metainfo['trialsInformation']] \
                        or 'naparm_path' not in [*_metainfo['trialsInformation']]:
                    raise ValueError(f'AllOpticalTrial experiment trial requires `tiff_path`, `paq_path` and `naparm_path` fields defined in .trialsInformation dictionary for each alloptical trial. '
                                     f'\n{self.trialsInformation[trial]}')
                if trial in self._trialsSuite2p:  # TODO could use switch statements in the 2p imaging trial class...
                    trial_obj = AllOpticalTrial(metainfo=_metainfo, naparm_path=_metainfo['trialsInformation']['naparm_path'],
                                                analysis_save_path=self.analysisSavePath, microscope=self.microscope, prestim_sec=1.0,
                                                poststim_sec=3.0, pre_stim_response_window=0.500, post_stim_response_window=0.500,
                                                total_frames_stitched=total_frames_stitched, suite2p_experiment_obj=self.Suite2p)
                else:
                    trial_obj = AllOpticalTrial(metainfo=_metainfo, naparm_path=_metainfo['trialsInformation']['naparm_path'],
                                                analysis_save_path=self.analysisSavePath, microscope=self.microscope, prestim_sec=1.0,
                                                poststim_sec=3.0, pre_stim_response_window=0.500, post_stim_response_window=0.500)

                # update self.trialsInformation using the information from new trial_obj
                self.trialsInformation[trial]['analysis_object_information'] = {'series ID': trial_obj.t_series_name,
                                                                                'repr': trial_obj.__repr__(),
                                                                                'pkl path': trial_obj.pkl_path  }
            else:
                raise ValueError(f"unsupported trial type for trial: {trial}. All trials must be of trialType: TwoPhotonImagingTrial or AllOpticalTrial")

            # initialize suite2p for trial objects
            if trial in self._trialsSuite2p:
                # trial_obj.Suite2p = _suite2p.Suite2pResultsTrial(suite2p_experiment_obj=self.Suite2p,
                #                                         trial_frames=(total_frames_stitched, total_frames_stitched + trial_obj.n_frames))  # use trial obj's current trial frames
                total_frames_stitched += trial_obj.imparams.n_frames
                trial_obj.save()


    @property
    def trialIDs(self):
        return [*self.trialsInformation]

    @property
    def tiff_path_dir(self):
        "retrieve the parent directory of the provided tiff_path"
        __first_trial_in_experiment = [*self.trialsInformation][0]
        __tiff_path_first_trial = self.trialsInformation[__first_trial_in_experiment]['tiff_path']
        return __tiff_path_first_trial[:[(s.start(), s.end()) for s in re.finditer('/', __tiff_path_first_trial)][-1][0]]  # this is the directory where the Bruker xml files associated with the 2p imaging TIFF are located

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

        :param pkl_path: full .pkl path for saving the Experiment object.
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
            trial_pkl_path = self.trialsInformation[trialID]['analysis_object_information']['pkl path']
            trialobj = _io.import_obj(trial_pkl_path)
            return trialobj
        except KeyError:
            raise KeyError("trialID not found in Experiment instance.")

    ## suite2p methods - refactored currently to _utils.Utils !!!!!
    def s2pRun(expobj, user_batch_size=2000, trialsSuite2P: list = None):  ## TODO gotta specify # of planes somewhere here
        """run suite2p on the experiment object, using trials specified in current experiment object, using the attributes
        determined directly from the experiment object."""

        ops = expobj.Suite2p.ops
        db = expobj.Suite2p.db

        expobj.Suite2p.trials = trialsSuite2P if trialsSuite2P else expobj.Suite2p.trials
        expobj._trialsSuite2p = trialsSuite2P if trialsSuite2P else expobj._trialsSuite2p

        tiffs_to_use_s2p = []
        for trial in expobj.Suite2p.trials:
            tiffs_to_use_s2p.append(expobj.trialsInformation[trial]['tiff_path'])

        sampling_rate = expobj.fps / expobj.n_planes
        diameter_x = 13 / expobj.pix_sz_x
        diameter_y = 13 / expobj.pix_sz_y
        diameter = int(diameter_x), int(diameter_y)
        expobj.Suite2p.user_batch_size = user_batch_size
        batch_size = expobj.Suite2p.user_batch_size * (262144 / (expobj.frame_x * expobj.frame_y))  # larger frames will be more RAM intensive, scale user batch size based on num pixels in 512x512 images

        if not db:
            db = {
                'fs': float(sampling_rate),
                'diameter': diameter,
                'batch_size': int(batch_size),
                'nimg_init': int(batch_size),
                'nplanes': expobj.n_planes
            }

        # specify tiff list to suite2p and data path
        db['tiff_list'] = tiffs_to_use_s2p
        db['data_path'] = expobj.dataPath  ## this is where the bad_frames.npy file will be stored for suite2p to use.
        db['save_folder'] = expobj._suite2p_save_path

        print(db)

        opsEnd = run_s2p(ops=ops, db=db)

        ## TODO update Experiment attr's and Trial attr's to reflect completion of the suite2p RUN
        expobj._s2pResultExists = True
        expobj.s2pResultsPath = expobj._suite2p_save_path + '/plane0/'  ## need to further debug that the flow of the suite2p path makes sense








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
