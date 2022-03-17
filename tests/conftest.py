import pytest
from packerlabimaging.utils.io import import_obj

SUITE2P_FRAMES_SPONT_t005t006 = [0, 14880]
SUITE2P_FRAMES_t013 = 103525

@pytest.fixture(scope="session")
def twophoton_imaging_trial_noPreDoneSuite2p_fixture():
    prep = 'PS12'
    date = '2021-01-25'

    initialization_dict = {
        'dataPath': f'/home/pshah/mnt/qnap/Data/{date}',
        # todo this seems very vauge, maybe add very specific documentation about what this is supposed to be, or just say tiff path?
        'analysisSavePath': f'/home/pshah/mnt/qnap/Analysis/{date}/{prep}/',
        'microscope': "Bruker",
        "expID": prep,
        'date': date,
        'comments': f'{prep} - interneuron gcamp imaging + LFP pre- and post-4ap',
        'TrialsInformation': {},  # NOTE: this dictionary is populated in the code cells below.
        'useSuite2p': True,
        's2pResultsPath': None  # todo unclear what to put here when suite2p hasn't been done yet...
    }

    # add information about each trial in experiment to trialsInformation field of the initialization_dict
    trials_list_pre4ap = ['t-001', 't-002', 't-003']
    # todo - add functionality to add longer detailed comments for each trial (e.g. t-001: 30 mins spont, t-002: 30 mins spont + LFP, etc.) (other than expGroup)

    for idx, trial in enumerate(trials_list_pre4ap):
        data_path_base = f'/home/pshah/mnt/qnap/Data/{date}'  # same as above, unclear what this is supposed to be.. (need more clear argument name)
        animal_prep = initialization_dict['expID']
        date = data_path_base[-10:]

        ## everything below should autopopulate
        paqs_loc = f'{data_path_base}/{date}_{animal_prep}_{trial[2:]}.paq'  # path to the .paq files for the selected trials
        tiffs_loc = f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'

        initialization_dict["TrialsInformation"][trial] = {'trialType': 'TwoPhotonImagingTrial',
                                                           # need some way of explaining the possible arguments accepted in trialType
                                                           'tiff_path': f"{tiffs_loc}",
                                                           's2p_use': True,
                                                           'expGroup': "pre 4ap 2p imaging",
                                                           'PaqInfoTrial': {'paq_path': paqs_loc,
                                                                            'frame_channel': 'frame_clock'}
                                                           }

    trials_list_post4ap = ['t-006', 't-007', 't-008', 't-009']
    for idx, trial in enumerate(trials_list_post4ap):
        data_path_base = f'/home/pshah/mnt/qnap/Data/{date}'  # same as above, unclear what this is supposed to be.. (need more clear argument name)
        animal_prep = initialization_dict['expID']
        date = data_path_base[-10:]

        ## everything below should autopopulate
        paqs_loc = f'{data_path_base}/{date}_{animal_prep}_{trial[2:]}.paq'  # path to the .paq files for the selected trials
        tiffs_loc = f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'

        initialization_dict["TrialsInformation"][trial] = {'trialType': 'TwoPhotonImagingTrial',
                                                           # need some way of explaining the possible arguments accepted in trialType
                                                           'tiff_path': f"{tiffs_loc}",
                                                           's2p_use': True,
                                                           'expGroup': "post 4ap 2p imaging",
                                                           'PaqInfoTrial': {'paq_path': paqs_loc,
                                                                            'frame_channel': 'frame_clock'}
                                                           }

    return initialization_dict


@pytest.fixture(scope="session")
def twophoton_imaging_trial_fixture():
    expobj = import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')

    initialization_dict = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'analysisSavePath': '/home/pshah/Documents/code/packerlabimaging/tests/',
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
                                                           'PaqInfoTrial': {
                                                               'frame_channel': 'frame_clock',
                                                               'paq_path': f'{data_path_base}/{date}_{animal_prep}_{trial[2:]}.paq'
                                                               # path to the .paq files for the selected trials
                                                           }
                                                           }
        _metainfo = {
            'exp_id': initialization_dict['expID'],
            'trial_id': trial,
            'date': initialization_dict['date'],
            't_series_id': f"{initialization_dict['expID']} {trial}",
            'TrialsInformation': initialization_dict["TrialsInformation"][trial]
        }
        initialization_dict['metainfo'] = _metainfo
        initialization_dict['analysis_save_path'] = initialization_dict['analysisSavePath']
        initialization_dict['suite2p_experiment_obj'] = expobj.Suite2p
        initialization_dict['paq_options'] = _metainfo['TrialsInformation']['PaqInfoTrial']
        initialization_dict['total_frames_stitched'] = SUITE2P_FRAMES_SPONT_t005t006[idx]

    return initialization_dict


@pytest.fixture(scope="session")
def alloptical_trial_fixture():
    expobj = import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')

    initialization_dict = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'analysisSavePath': '/home/pshah/Documents/code/packerlabimaging/tests/',
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
                                                           'PaqInfoTrial': {'frame_channel': 'frame_clock',
                                                                            'paq_path': f'{data_path_base}/{date}_{animal_prep}_{trial[2:]}.paq',
                                                                            # path to the .paq files for the selected trials
                                                                            'stim_channel': 'markpoints2packio'
                                                                            },
                                                           'naparm_path': naparm_paths[trial]
                                                           }

        _metainfo = {
            'exp_id': initialization_dict['expID'],
            'trial_id': trial,
            'date': initialization_dict['date'],
            't_series_id': f"{initialization_dict['expID']} {trial}",
            'TrialsInformation': initialization_dict["TrialsInformation"][trial]
        }

        initialization_dict['metainfo'] = _metainfo
        initialization_dict['naparm_path'] = initialization_dict["TrialsInformation"][trial]['naparm_path']
        initialization_dict['analysis_save_path'] = initialization_dict['analysisSavePath']
        initialization_dict['suite2p_experiment_obj'] = expobj.Suite2p
        initialization_dict['paq_options'] = _metainfo['TrialsInformation']['PaqInfoTrial']
        initialization_dict['total_frames_stitched'] = SUITE2P_FRAMES_t013

    return initialization_dict


@pytest.fixture(scope="session")
def experiment_fixture(alloptical_trial_fixture, twophoton_imaging_trial_fixture):
    initialization_dict = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'analysisSavePath': '/home/pshah/Documents/code/packerlabimaging/tests/',
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

