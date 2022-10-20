import numpy as np
import pandas as pd

from conftest import existing_expobj_nopredones2p_fixture
from imagingplus import Experiment, import_obj
from imagingplus._archive.TwoPhotonImagingMain import TwoPhotonImagingTrial
from imagingplus.processing import suite2p
from imagingplus.processing.suite2p import Suite2pExperiment, Suite2pResultsTrial

# todo this test is likely breaking!!
from imagingplus.workflows import TwoPhotonImaging


def test_Suite2pResultsExperiment(existing_expobj_fixture):
    expobj: Experiment = existing_expobj_fixture
    expobj.save()

def test_Suite2pResultsTrial(existing_trialobj_twophotonimaging_fixture, existing_expobj_fixture,):
                             # existing_trialobj_alloptical_fixture):
    trialobj = existing_trialobj_twophotonimaging_fixture
    # alloptical_trialobj = existing_trialobj_alloptical_fixture
    expobj: Experiment = existing_expobj_fixture

    frames = 0
    for n_obj in [trialobj]:
        from imagingplus.processing.suite2p import Suite2pExperiment
        s2p_expobj: Suite2pExperiment = expobj.Suite2p
        n_obj.Suite2p = suite2p.Suite2pResultsTrial(s2pExp=s2p_expobj, trial_frames=(frames, n_obj.n_frames))  # use trial obj's current trial key_frames
        frames += n_obj.n_frames
        n_obj.save()


def test_Suite2pExp():
    # TODO run test for suite2p without predone results path
    _trialsSuite2p = ['t-005', 't-006']

    from imagingplus import import_obj
    expobj = import_obj(pkl_path='/home/pshah/mnt/qnap/Analysis/2021-01-25/PS12/PS12_analysis.pkl')

    expobj.Suite2p = Suite2pExperiment(trialsTiffsSuite2p=expobj.Suite2p.tiff_paths_to_use_s2p,
                                       s2pResultsPath=expobj.Suite2p.s2pResultsPath,
                                       subtract_neuropil=expobj.Suite2p.subtract_neuropil)
    expobj.save()
# test_Suite2pExp()


def test_Suite2pTrial():
    # TODO run tests for suite2p trial creation without predone results path
    trialsSuite2p: list = ['t-005', 't-006']

    trial_frames: tuple
    s2ptrial = Suite2pResultsTrial(trialsTiffsSuite2p=trialsSuite2p, trial_frames=(0, 1000))
    assert s2ptrial._s2pResultExists == False  # - doesnt seem to be setting to False as expected during the super() call!

# test_Suite2pTrial()

def test_add_bad_frames():
    expobj: Experiment = existing_expobj_nopredones2p_fixture()

    expobj.Suite2p.bad_frames = []
    for trial in expobj.trialIDs:
        if 'post' in expobj.TrialsInformation[trial]['expGroup']:
            trialobj: TwoPhotonImagingTrial = expobj.load_trial(trialID=trial)
            print(len(expobj.Suite2p.bad_frames))

            expobj.Suite2p.add_bad_frames(
                frames=np.arange(trialobj.Suite2p.trial_frames[0], trialobj.Suite2p.trial_frames[1]),
                bad_frames_npy_loc=expobj.dataPath)
    assert len(expobj.Suite2p.bad_frames) > 0
    expobj.save()


def test_update_ops():
    # TODO run test
    expobj: Experiment = existing_expobj_nopredones2p_fixture()
    trialobj: TwoPhotonImagingTrial = expobj.load_trial(trialID='t-001')


    pix_sz_x = trialobj.imparams.pix_sz_x
    pix_sz_y = trialobj.imparams.pix_sz_y
    diameter_x = 13 / pix_sz_x
    diameter_y = 13 / pix_sz_y
    diameter = (int(diameter_x), int(diameter_y)) if int(diameter_y) != int(diameter_x) else int(diameter_x)

    # set imaging parameters using defaults or kwargs if provided
    frame_x = trialobj.imparams.frame_x
    frame_y = trialobj.imparams.frame_y

    # new values for ops dictionary
    new_ops = {
        'fs': trialobj.imparams.fps,
        'tau': 1.15,
        'num_workers': 50,
        'diameter': diameter,
        'delete_bin': False,  # temporary not deleting binaries in case needed for further testing!
        'batch_size': 2000 * (262144 / (
            frame_x * frame_y)),  # larger key_frames will be more RAM intensive, scale user batch size based on num pixels in 512x512 images
        'save_folder': expobj.suite2p_save_path
    }

    expobj.Suite2p.update_ops(new_ops)
    expobj.save()

    from imagingplus import import_obj
    expobj = import_obj(pkl_path='/home/pshah/mnt/qnap/Analysis/2021-01-25/PS12/PS12_analysis.pkl')
    assert expobj.Suite2p.ops['tau'] == 1.15
    print('ops:\n', expobj.Suite2p.ops)
    print('db:\n', expobj.Suite2p.db)
# test_update_ops()


def test_s2pRun():
    expobj: Experiment = existing_expobj_nopredones2p_fixture()
    expobj.Suite2p.s2pRun(expobj=expobj)
    if len(expobj.Suite2p.bad_frames) == 0: test_add_bad_frames()

    # TODO run suite2p - but probably want to select only a subset of trials to run for testing!
    expobj.Suite2p.s2pRun(expobj=expobj)

# test_s2pRun()


def test_FrameAverage():
    from imagingplus.utils.io import import_obj
    expobj = import_obj(pkl_path='/home/pshah/mnt/qnap/Analysis/2021-01-25/PS12/PS12_analysis.pkl')
    trialobj = expobj.load_trial(trialID=expobj.trialIDs[0])

    trialobj.Suite2p.FrameAverage(key_frames=[110, 510], peri_frames=100,
                                  save_path=trialobj.saveDir + f'/export/avg_frames/', to_plot=True)


def test_makeDownSampledTiff():
    from imagingplus.utils.io import import_obj
    expobj = import_obj(pkl_path='/home/pshah/mnt/qnap/Analysis/2021-01-25/PS12/PS12_analysis.pkl')
    reg_tiff_folder = '/home/pshah/mnt/qnap/Analysis/2021-01-25/PS12/suite2p/reg_tif'
    trialobj = expobj.load_trial(trialID=expobj.trialIDs[0])

    trialobj.Suite2p.makeDownSampledTiff(key_frames=[110, 510], peri_frames=100,
                                  save_path=trialobj.saveDir + f'/export/avg_frames/', to_plot=True)

# test_makeFrameAverageTiff()

def test_stitch_s2p_reg_tiff(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImaging = existing_trialobj_twophotonimaging_fixture
    assert hasattr(trialobj, 'Suite2p')
    trialobj.Suite2p


def getCellsAnnotations(self: TwoPhotonImaging):
    # if self.s2pResultExists:
    # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
    # build dataframe for obs_meta from suite2p stat information
    obs_meta = pd.DataFrame(
        columns=['original_index', 'footprint', 'mrs', 'mrs0', 'compact', 'med', 'npix', 'radius',
                 'aspect_ratio', 'npix_norm', 'skew', 'std'], index=range(len(self.stat)))
    for idx, __stat in enumerate(self.stat):
        for __column in obs_meta:
            obs_meta.loc[idx, __column] = __stat[__column]

    obs_m = {'ypix': [],
             'xpix': []}
    for col in [*obs_m]:
        for idx, __stat in enumerate(self.stat):
            obs_m[col].append(__stat[col])
        obs_m[col] = np.asarray(obs_m[col])

    return obs_meta, obs_m

    # else:
    #     raise ValueError('cannot set cell annotations. no s2p results found in trial.')

def test_getCellsAnnotations(existing_trialobj_twophotonimaging_fixture):
    # FminusFneu, spks, stat, neuropil = suite2p_results_fixture
    expobj = import_obj(pkl_path='/home/pshah/mnt/qnap/Analysis/2021-01-31/HF113/HF113_analysis.pkl')
    s2p_path = f'/home/pshah/mnt/qnap/Analysis/2021-01-31/HF113//suite2p//plane0/'

    trialobj: TwoPhotonImaging = existing_trialobj_twophotonimaging_fixture

    trialobj.Suite2p.getCellsAnnotations()


