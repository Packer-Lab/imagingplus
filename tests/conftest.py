import os.path

import numpy as np
import pandas as pd
import pytest

from packerlabimaging.utils.io import import_obj

LOCAL_DATA_PATH = '/Users/prajayshah/data/oxford-data-to-process/'
REMOTE_DATA_PATH = '/home/pshah/mnt/qnap/Data/'
BASE_PATH = LOCAL_DATA_PATH


SUITE2P_FRAMES_SPONT_t005t006 = [0, 14880]
SUITE2P_FRAMES_t013 = 103525


@pytest.fixture(scope="session")
def twophoton_imaging_trial_new_noPreDoneSuite2p_fixture():
    """apr 30 2022 - newest v0.2.0 structure"""

    date = '2021-01-31'
    prep = 'HF113'

    # add information about each trial in experiment to trialsInformation field of the initialization_dict
    trials_paq = {'t-001': '001.paq',
                  't-002': '002.paq',
                  't-003': '003.paq',
                  't-004': '004.paq'}

    return trials_paq


@pytest.fixture(scope="session")
def twophoton_imaging_multitrial_noPreDoneSuite2p_fixture():
    prep = 'HF113'
    date = '2021-01-31'

    trials_paq = {'t-001': '001.paq',
                  't-002': '002.paq',
                  't-003': '003.paq',
                  't-004': '004.paq'
                  }

    info = {
        'prep': prep,
        'date': date,
        'trials_paq': trials_paq
    }

    return info


@pytest.fixture(scope="session")
def twophoton_imaging_trial_fixture():
    expobj = import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')

    initialization_dict = {
        'dataPath': '/home/pshah/mnt/qnap/Data/2020-12-19',
        'saveDir': '/home/pshah/Documents/code/packerlabimaging/tests/',
        'microscope': "Bruker",
        "expID": 'RL109',
        'date': '2020-12-19',
        'comment': 'testing out analysis workflow',
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


# @pytest.fixture(scope="session")
def alloptical_trial_fixture():
    expobj = import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')

    initialization_dict = {'naparm_path': f'{BASE_PATH}/2020-12-19/photostim/2020-12-19_RL109_ps_014/',
                           'dataPath': f'{BASE_PATH}/2020-12-19/2020-12-19_t-013/2020-12-19_t-013_Cycle00001_Ch3.tif',
                           'saveDir': f'{BASE_PATH}/2020-12-19/',
                           'date': '2020-12-19',
                           'trialID': 't-013',
                           'expID': 'RL109',
                           'expGroup': 'all optical trial with LFP',
                           'comment': ''}

    return initialization_dict


@pytest.fixture(scope="session")
def experiment_fixture():
    ExperimentMetainfo = {
        'dataPath': '/mnt/qnap_share/Data/packerlabimaging-example/packerlabimaging-test-cellsdata',
        'saveDir': '/mnt/qnap_share/Data/packerlabimaging-example/packerlabimaging-test-analysis',
        "expID": 'HF113',
        'comment': 'two photon imaging + LFP dataset',
    }
    return ExperimentMetainfo


@pytest.fixture(scope="session")
def existing_trialobj_twophotonimaging_fixture():
    expobj = import_obj(
        pkl_path='/mnt/qnap_share/Data/packerlabimaging-example/packerlabimaging-test-analysis/HF113_analysis.pkl')
    trialobj1 = expobj.load_trial(expobj.trialIDs[0])
    return trialobj1


# @pytest.fixture(scope="session")
def existing_trialobj_alloptical_fixture():
    expobj = import_obj(pkl_path='/mnt/qnap_share/Data/packerlabimaging-example/RL109_analysis.pkl')
    trialobj = expobj.load_trial('t-013')
    return trialobj


@pytest.fixture(scope="session")
def existing_expobj_fixture():
    expobj = import_obj(pkl_path='/mnt/qnap_share/Data/packerlabimaging-example/RL109_analysis.pkl')
    return expobj


@pytest.fixture(scope="session")
def suite2p_results_fixture():
    expobj = import_obj(pkl_path='/home/pshah/mnt/qnap/Analysis/2021-01-31/HF113/HF113_analysis.pkl')
    s2p_path = f'/home/pshah/mnt/qnap/Analysis/2021-01-31/HF113//suite2p//plane0/'
    assert os.path.exists(s2p_path), 's2p path not found...'
    from packerlabimaging.processing.suite2p import s2p_loader
    FminusFneu, spks, stat, neuropil = s2p_loader(s2p_path=s2p_path)
    return FminusFneu, spks, stat, neuropil@pytest.fixture(scope="session")


# @pytest.fixture(scope="session")
# def existing_trialobj_alloptical_fixture():
#     expobj = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
#     trialobj = expobj.load_trial(trialID='t-013')
#     return trialobj


@pytest.fixture(scope="session")
def existing_expobj_nopredones2p_fixture():
    expobj = import_obj(pkl_path='/mnt/qnap_share/Data/packerlabimaging-example/RL109_analysis.pkl')
    return expobj

@pytest.fixture(scope="session")
def anndata_trial_data():
    aotrial = expobj.load_trial(trialID = 't-013')
    aotrial.data = aotrial.create_anndata(imdata=aotrial.Suite2p,
                                          cells=aotrial.Suite2p,
                                          tmdata=aotrial.tmdata,
                                          imdata_type='suite2p raw - neuropil corrected')

@pytest.fixture(scope='session')
def existing_anndata():
    expobj = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    trialobj = expobj.load_trial(trialID=expobj.trialIDs[0])

    print(trialobj.cellsdata)  # this is the anndata object for this trial

    var_meta = pd.DataFrame({
        'exp_group': np.random.choice(['A', 'B', 'C'], trialobj.n_frames),
    },
        index=np.arange(trialobj.n_frames, dtype=int).astype(str),  # these are the same IDs of observations as above!
    )

    trialobj.cellsdata.add_var(var_name='exp_group', values=list(var_meta['exp_group']))

    print(trialobj.cellsdata)


@pytest.fixture(scope='session')
def s_tiff_path_fixture():
    return '/home/pshah/mnt/qnap/Data/2021-01-31/2021-01-31_s-007/2021-01-31_s-007_Cycle00001_Ch2_000001.ome.tif'


@pytest.fixture(scope='session')
def tiff_path_fixture():
    # return '/home/pshah/mnt/qnap/Data/2021-01-31/2021-01-31_t-006/2021-01-31_t-006_Cycle00001_Ch3.tif'
    # return '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_t-004/2021-01-28_PS14_t-004_Cycle00001_Ch3.tif'
    return '/home/pshah/mnt/qnap/Data/2020-03a/2020-03-03/2020-03-03_t-001/2020-03-03_t-001_Cycle00001_Ch3.tif'

@pytest.fixture(scope='session')
def image_registration_fixture():
    ref_image_path = '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_s-013/2021-01-28_PS14_s-013_Cycle00001_Ch2_000001.ome.tif'
    multi_tiff_path = '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_t-004/2021-01-28_PS14_t-004_Cycle00001_Ch3.tif'
    return (ref_image_path, multi_tiff_path)



