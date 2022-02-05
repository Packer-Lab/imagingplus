import pytest
from packerlabimaging._io import import_obj


@pytest.fixture(scope="session")
def twophoton_imaging_trial_fixture():
    expobj = import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    SUITE2P_FRAMES = 0

    initialization_dict = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'analysisSavePath': '/home/pshah/Documents/code/packerlabimaging/tests/',
        'microscope': "Bruker",
        "expID": 'RL109',
        'date': '2020-12-19',
        'comments': 'testing out analysis workflow',
        'trialsInformation': {},  # NOTE: this dictionary is populated in the code cells below.
        # 'useSuite2p': True,
        # 'useSuite2p': False,
        's2pResultsPath': "/home/pshah/mnt/qnap/Analysis/2020-12-19/suite2p/alloptical-2p-1x-alltrials/plane0"
    }

    # add information about each trial in experiment to trialsInformation field of the initialization_dict
    trials_list_spont = ['t-005', 't-006']
    for idx, trial in enumerate(trials_list_spont):
        data_path_base = '/home/pshah/mnt/qnap/Data/2020-12-19'
        animal_prep = initialization_dict['expID']
        date = data_path_base[-10:]

        ## everything below should autopopulate and run automatically
        paqs_loc = '%s/%s_%s_%s.paq' % (
        data_path_base, date, animal_prep, trial[2:])  # path to the .paq files for the selected trials
        tiffs_loc = f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'

        initialization_dict["trialsInformation"][trial] = {'trialType': 'TwoPhotonImagingTrial',
                                                           'tiff_path': f"{tiffs_loc}",
                                                           's2p_use': True,
                                                           'expGroup': "pre 4ap 2p spont imaging",
                                                           'paq_path': paqs_loc
                                                           }
        _metainfo = {
            'animal prep.': initialization_dict['expID'],
            'trial': trial,
            'date': initialization_dict['date'],
            't series id': f"{initialization_dict['expID']} {trial}",
            'trialsInformation': initialization_dict["trialsInformation"][trial]
        }
        initialization_dict['metainfo'] = _metainfo
        initialization_dict['analysis_save_path'] = initialization_dict['analysisSavePath']
        initialization_dict['suite2p_experiment_obj'] = expobj.Suite2p
        initialization_dict['total_frames_stitched'] = SUITE2P_FRAMES


    return initialization_dict

@pytest.fixture(scope="session")
def alloptical_trial_fixture():
    expobj = import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    SUITE2P_FRAMES = 0

    initialization_dict = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'analysisSavePath': '/home/pshah/Documents/code/packerlabimaging/tests/',
        'microscope': "Bruker",
        "expID": 'RL109',
        'date': '2020-12-19',
        'comments': 'testing out analysis workflow',
        'trialsInformation': {},  # NOTE: this dictionary is populated in the code cells below.
        # 'useSuite2p': True,
        # 'useSuite2p': False,
        's2pResultsPath': "/home/pshah/mnt/qnap/Analysis/2020-12-19/suite2p/alloptical-2p-1x-alltrials/plane0"
    }

    # add information about each trial in experiment to trialsInformation field of the initialization_dict
    trials_list_alloptical = ['t-013']
    naparm_paths = {'t-013': '/home/pshah/mnt/qnap/Data/2020-12-19/photostim/2020-12-19_RL109_ps_014/'}
    for idx, trial in enumerate(trials_list_alloptical):
        data_path_base = '/home/pshah/mnt/qnap/Data/2020-12-19'
        animal_prep = initialization_dict['expID']
        date = data_path_base[-10:]

        ## everything below should autopopulate and run automatically
        paqs_loc = '%s/%s_%s_%s.paq' % (
        data_path_base, date, animal_prep, trial[2:])  # path to the .paq files for the selected trials
        tiffs_loc = f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'

        initialization_dict["trialsInformation"][trial] = {'trialType': 'AllOpticalTrial',
                                                           'tiff_path': f"{tiffs_loc}",
                                                           's2p_use': True,
                                                           'expGroup': "pre 4ap 2p all optical",
                                                           'paq_path': paqs_loc,
                                                           'naparm_path': naparm_paths[trial]
                                                           }

        _metainfo = {
            'animal prep.': initialization_dict['expID'],
            'trial': trial,
            'date': initialization_dict['date'],
            't series id': f"{initialization_dict['expID']} {trial}",
            'trialsInformation': initialization_dict["trialsInformation"][trial]
        }

        initialization_dict['metainfo'] = _metainfo
        initialization_dict['naparm_path'] = initialization_dict["trialsInformation"][trial]['naparm_path']
        initialization_dict['analysis_save_path'] = initialization_dict['analysisSavePath']
        initialization_dict['suite2p_experiment_obj'] = expobj.Suite2p
        initialization_dict['total_frames_stitched'] = SUITE2P_FRAMES

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
        'trialsInformation': {},  # NOTE: this dictionary is populated in the code cells below.
        # 'useSuite2p': True,
        # 'useSuite2p': False,
        's2pResultsPath': "/home/pshah/mnt/qnap/Analysis/2020-12-19/suite2p/alloptical-2p-1x-alltrials/plane0"
    }
    
    for trial in [*twophoton_imaging_trial_fixture['trialsInformation']]:
        initialization_dict['trialsInformation'][trial] = twophoton_imaging_trial_fixture['trialsInformation'][trial]
    
    for trial in [*alloptical_trial_fixture['trialsInformation']]:
        initialization_dict['trialsInformation'][trial] = alloptical_trial_fixture['trialsInformation'][trial]
    
    return initialization_dict

@pytest.fixture(scope="session")
def existing_trialobj_twophotonimaging_fixture():
    expobj = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    trialobj = expobj.load_trial(trialID='t-005')
    return trialobj