import os.path

import numpy as np
import pandas as pd
import pytest

from packerlabimaging.processing.imagingMetadata import PrairieViewMetadata
from packerlabimaging.main.paq import PaqData
from packerlabimaging.utils.io import import_obj

SUITE2P_FRAMES_SPONT_t005t006 = [0, 14880]
SUITE2P_FRAMES_t013 = 103525

@pytest.fixture(scope="session")
def twophoton_imaging_trial_new_noPreDoneSuite2p_fixture():
    "mar 28 2022 - new no pipeline structure"

    ExperimentMetainfo = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'saveDir': '/home/pshah/Documents/code/packerlabimaging/tests/',
        "expID": 'RL109',
        'date': '2020-12-19',
        'comments': 'testing out analysis workflow',
    }

    trials_info = {}


    trials_list_spont = ['t-005', 't-006']
    for idx, trial in enumerate(trials_list_spont):
        data_path_base = '/home/pshah/mnt/qnap/Data/2020-12-19'
        animal_prep = ExperimentMetainfo['expID']
        date = data_path_base[-10:]

        ## everything below should autopopulate
        paqs_loc = '%s/%s_%s_%s.paq' % (
        data_path_base, date, animal_prep, trial[2:])  # path to the .paq files for the selected trials
        tiffs_loc = f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'

        TwoPhotonImagingMetainfo = {
            'date': ExperimentMetainfo['date'],
            'trialID': trial,
            'expID': ExperimentMetainfo['expID'],
            'microscope': 'Bruker 2pPlus',
            'tiff_path': tiffs_loc,
            'saveDir': ExperimentMetainfo['saveDir'],
            'expGroup': "pre 4ap 2p spont imaging",
            'PaqInfo': {'paq_path': paqs_loc,
                             'frame_channel': 'frame_clock'}
        }

        trials_info[trial] = TwoPhotonImagingMetainfo

    return trials_info


@pytest.fixture(scope="session")
def twophoton_imaging_multitrial_noPreDoneSuite2p_fixture():
    prep = 'HF113'
    date = '2021-01-31'

    trials_paq = {'t-001': '001.paq',
                  't-002': '002.paq',
                  't-003': '003.paq',
                  't-004': '004.paq'}

    TwoPhotonImagingMultiTrial = {}
    for trial, paq in trials_paq.items():
        data_path_base = f'/home/pshah/mnt/qnap/Data/{date}'  # same as above, unclear what this is supposed to be.. (need more clear argument name)

        paqs_loc = f'{data_path_base}/{date}_{prep}_{paq}'  # path to the .paq files for the selected trials
        dataPath = f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'

        imparams = PrairieViewMetadata(pv_xml_dir=os.path.dirname(dataPath), microscope='Bruker 2pPlus')
        tmdata = PaqData.paqProcessingTwoPhotonImaging(paq_path=paqs_loc, frame_channel='frame_clock')


        TwoPhotonImagingMultiTrial[trial] = {'date': date,
                                    'trialID': trial,
                                    'expID': prep,
                                    'imparams': imparams,
                                    'tmdata': tmdata,
                                    'saveDir': f'/home/pshah/mnt/qnap/Analysis/{date}/{prep}/',
                                    'dataPath': dataPath,
                                    'expGroup': "awake spont. 2p imaging + LFP",
                                    }



    return TwoPhotonImagingMultiTrial

@pytest.fixture(scope="session")
def twophoton_imaging_trial_fixture():
    expobj = import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')

    initialization_dict = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'saveDir': '/home/pshah/Documents/code/packerlabimaging/tests/',
        'microscope': "Bruker",
        "expID": 'RL109',
        'date': '2020-12-19',
        'comments': 'testing out analysis workflow',
        'TrialsInformation': {},  # NOTE: this dictionary is populated in the code cells below.
        # 'useSuite2p': True,
        # 'useSuite2p': False,
        's2pResultsPath': "/home/pshah/mnt/qnap/Analysis/2020-12-19/suite2p/alloptical-2p-1x-alltrials/plane0"
    }

    # add information about each trial in experiment to TrialsInformation field of the initialization_dict
    trials_list_spont = ['t-005', 't-006']
    for idx, trial in enumerate(trials_list_spont):
        data_path_base = '/home/pshah/mnt/qnap/Data/2020-12-19'
        animal_prep = initialization_dict['expID']
        date = data_path_base[-10:]

        # create dictionary containing necessary information for initialization
        initialization_dict["TrialsInformation"][trial] = {'trialType': 'TwoPhotonImagingTrial',
                                                           'tiff_path': f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif',
                                                           's2p_use': True,
                                                           'expGroup': "pre 4ap 2p spont imaging",
                                                           'PaqInfo': {
                                                               'frame_channel': "frame_clock",
                                                               'paq_path': f'{data_path_base}/{date}_{animal_prep}_{trial[2:]}.paq'
                                                               # path to the .paq files for the selected trials
                                                           }
                                                           }
        _metainfo = {
            'expID': initialization_dict['expID'],
            'trialID': trial,
            'date': initialization_dict['date'],
            'TrialsInformation': initialization_dict["TrialsInformation"][trial]
        }
        initialization_dict['metainfo'] = _metainfo
        initialization_dict['analysis_save_path'] = initialization_dict['saveDir']
        initialization_dict['suite2p_experiment_obj'] = expobj.Suite2p
        initialization_dict['paq_options'] = _metainfo['TrialsInformation']['PaqInfo']
        initialization_dict['total_frames_stitched'] = SUITE2P_FRAMES_SPONT_t005t006[idx]

    return initialization_dict


@pytest.fixture(scope="session")
def alloptical_trial_fixture():
    expobj = import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')

    initialization_dict = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'saveDir': '/home/pshah/Documents/code/packerlabimaging/tests/',
        'microscope': "Bruker",
        "expID": 'RL109',
        'date': '2020-12-19',
        'comments': 'testing out analysis workflow',
        'TrialsInformation': {},  # NOTE: this dictionary is populated in the code cells below.
        # 'useSuite2p': True,
        # 'useSuite2p': False,
        's2pResultsPath': "/home/pshah/mnt/qnap/Analysis/2020-12-19/suite2p/alloptical-2p-1x-alltrials/plane0"
    }

    # add information about each trial in experiment to TrialsInformation field of the initialization_dict
    trials_list_alloptical = ['t-013']
    naparm_paths = {'t-013': '/home/pshah/mnt/qnap/Data/2020-12-19/photostim/2020-12-19_RL109_ps_014/'}
    for idx, trial in enumerate(trials_list_alloptical):
        data_path_base = '/home/pshah/mnt/qnap/Data/2020-12-19'
        animal_prep = initialization_dict['expID']
        date = data_path_base[-10:]

        ## create dictionary containing necessary information for initialization
        initialization_dict["TrialsInformation"][trial] = {'trialType': 'AllOpticalTrial',
                                                           'tiff_path': f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif',
                                                           's2p_use': True,
                                                           'expGroup': "pre 4ap 2p all optical",
                                                           'PaqInfo': {'frame_channel': "frame_clock",
                                                                            'paq_path': f'{data_path_base}/{date}_{animal_prep}_{trial[2:]}.paq',
                                                                            # path to the .paq files for the selected trials
                                                                            'stim_channel': 'markpoints2packio'
                                                                            },
                                                           'naparm_path': naparm_paths[trial]
                                                           }

        _metainfo = {
            'expID': initialization_dict['expID'],
            'trialID': trial,
            'date': initialization_dict['date'],
            'TrialsInformation': initialization_dict["TrialsInformation"][trial]
        }

        initialization_dict['metainfo'] = _metainfo
        initialization_dict['naparm_path'] = initialization_dict["TrialsInformation"][trial]['naparm_path']
        initialization_dict['analysis_save_path'] = initialization_dict['saveDir']
        initialization_dict['suite2p_experiment_obj'] = expobj.Suite2p
        initialization_dict['paq_options'] = _metainfo['TrialsInformation']['PaqInfo']
        initialization_dict['total_frames_stitched'] = SUITE2P_FRAMES_t013

    return initialization_dict


@pytest.fixture(scope="session")
def experiment_fixture(alloptical_trial_fixture, twophoton_imaging_trial_fixture):
    initialization_dict = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'saveDir': '/home/pshah/Documents/code/packerlabimaging/tests/',
        'microscope': "Bruker",
        "expID": 'RL109',
        'date': '2020-12-19',
        'comments': 'testing out analysis workflow',
        'TrialsInformation': {},  # NOTE: this dictionary is populated in the code cells below.
        # 'useSuite2p': True,
        # 'useSuite2p': False,
        's2pResultsPath': "/home/pshah/mnt/qnap/Analysis/2020-12-19/suite2p/alloptical-2p-1x-alltrials/plane0"
    }

    for trial in [*twophoton_imaging_trial_fixture['TrialsInformation']]:
        initialization_dict['TrialsInformation'][trial] = twophoton_imaging_trial_fixture['TrialsInformation'][trial]

    for trial in [*alloptical_trial_fixture['TrialsInformation']]:
        initialization_dict['TrialsInformation'][trial] = alloptical_trial_fixture['TrialsInformation'][trial]

    return initialization_dict


@pytest.fixture(scope="session")
def experimentnew_fixture():
    ExperimentMetainfo = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'saveDir': '/home/pshah/Documents/code/packerlabimaging/tests/',
        "expID": 'RL109',
        'date': '2020-12-19',
        'comments': 'testing out analysis workflow',
    }

    return ExperimentMetainfo




@pytest.fixture(scope="session")
def existing_trialobj_twophotonimaging_fixture():
    expobj = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    trialobj1 = expobj.load_trial(trialID='t-005')
    trialobj2 = expobj.load_trial(trialID='t-006')
    return trialobj1, trialobj2


@pytest.fixture(scope="session")
def existing_expobj_fixture():
    expobj = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    trialsSuite2p = ['t-005', 't-006', 't-013']
    s2pResultsPath = '/home/pshah/mnt/qnap/Analysis/2020-12-19/suite2p/alloptical-2p-1x-alltrials/plane0'
    return expobj, trialsSuite2p, s2pResultsPath


@pytest.fixture(scope="session")
def existing_trialobj_alloptical_fixture():
    expobj = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    trialobj = expobj.load_trial(trialID='t-013')
    return trialobj


# @pytest.fixture(scope="session")
def existing_expobj_nopredones2p_fixture():
    expobj = import_obj(pkl_path='/home/pshah/mnt/qnap/Analysis/2021-01-25/PS12/PS12_analysis.pkl')
    return expobj


def anndata_trial_data():
    import pandas as pd
    import numpy as np

    # number of observations
    n_obs = 1000
    # say we measure the time of observing the data points
    # add them to a dataframe for storing some annotation
    obs = pd.DataFrame()
    obs['group'] = np.random.choice(['day 1', 'day 2', 'day 4', 'day 8'], n_obs)
    obs['group2'] = np.random.choice(['day 3', 'day 5', 'day 7'], n_obs)
    # set the names of variables/features to the following
    # ['A', 'B', 'C', ..., 'AA', 'BB', 'CC', ..., 'AAA', ...]
    from string import ascii_uppercase
    var_names = [i * letter for i in range(1, 10) for letter in ascii_uppercase]
    # number of variables
    n_vars = len(var_names)
    var_group = {'var_group_1': np.random.choice(['group A', 'group B', 'group C', 'group D'], n_vars),
                 'var_group_2': np.random.choice(['group A', 'group B', 'group C', 'group D'], n_vars)}
    # dataframe for annotating the variables
    var = pd.DataFrame(var_group, index=var_names)
    # the data matrix of shape n_obs x n_vars
    # X = np.arange(n_obs * n_vars).reshape(n_obs, n_vars)
    X = np.random.random(n_obs * n_vars).reshape(n_obs, n_vars)

    return X, var, obs

@pytest.fixture(scope='session')
def existing_anndata():
    expobj = import_obj(
        pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    trialobj = expobj.load_trial(trialID=expobj.trialIDs[0])

    print(trialobj.data)  # this is the anndata object for this trial

    var_meta = pd.DataFrame({
        'exp_group': np.random.choice(['A', 'B', 'C'], trialobj.n_frames),
    },
        index=np.arange(trialobj.n_frames, dtype=int).astype(str),  # these are the same IDs of observations as above!
    )
    # var_meta

    trialobj.data.add_var(var_name='exp_group', values=list(var_meta['exp_group']))

    print(trialobj.data)




































































