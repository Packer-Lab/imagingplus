import numpy as np
from conftest import existing_expobj_nopredones2p_fixture
from packerlabimaging import Experiment, TwoPhotonImagingTrial
from packerlabimaging.processing import suite2p
from packerlabimaging.processing.suite2p import Suite2pResultsExperiment, Suite2pResultsTrial


def test_Suite2pResultsExperiment(existing_expobj_fixture):
    expobj: Experiment = existing_expobj_fixture[0]
    trialSuite2p = existing_expobj_fixture[1]
    s2pResultsPath = existing_expobj_fixture[2]

    expobj.Suite2p = suite2p.Suite2pResultsExperiment(trialsTiffsSuite2p=trialSuite2p,
                                                      s2pResultsPath=s2pResultsPath)


def test_Suite2pResultsTrial(existing_trialobj_twophotonimaging_fixture, existing_trialobj_alloptical_fixture,
                             existing_expobj_fixture):
    trialobj, trialobj_ = existing_trialobj_twophotonimaging_fixture
    alloptical_trialobj = existing_trialobj_alloptical_fixture
    expobj: Experiment = existing_expobj_fixture[0]

    for n_obj in [trialobj, trialobj_, alloptical_trialobj]:
        from packerlabimaging.processing.suite2p import Suite2pResultsExperiment
        s2p_expobj: Suite2pResultsExperiment = expobj.Suite2p
        n_obj.Suite2p = suite2p.Suite2pResultsTrial(trialsTiffsSuite2p=s2p_expobj.tiff_paths_to_use_s2p,
                                                    s2pResultsPath=s2p_expobj.s2pResultsPath,
                                                    subtract_neuropil=s2p_expobj.subtract_neuropil,
                                                    trial_frames=n_obj.Suite2p.trial_frames)  # use trial obj's current trial frames

        n_obj.save()

def test_Suite2pExp():
    # TODO run test for suite2p without predone results path
    _trialsSuite2p = ['t-005', 't-006']

    from packerlabimaging import import_obj
    expobj = import_obj(pkl_path='/home/pshah/mnt/qnap/Analysis/2021-01-25/PS12/PS12_analysis.pkl')

    expobj.Suite2p = Suite2pResultsExperiment(trialsTiffsSuite2p=expobj.Suite2p.tiff_paths_to_use_s2p,
                                              s2pResultsPath=expobj.Suite2p.s2pResultsPath,
                                              subtract_neuropil=expobj.Suite2p.subtract_neuropil)
    expobj.save()
# test_Suite2pExp()


def test_Suite2pTrial():
    # TODO run tests for suite2p trial creation without predone results path
    trialsSuite2p: list = ['t-005', 't-006']

    trial_frames: tuple
    s2ptrial = Suite2pResultsTrial(trialsTiffsSuite2p=trialsSuite2p, trial_frames=(0, 1000))
    assert s2ptrial.__s2pResultExists == False  # - doesnt seem to be setting to False as expected during the super() call!

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
            frame_x * frame_y)),  # larger frames will be more RAM intensive, scale user batch size based on num pixels in 512x512 images
        'save_folder': expobj.suite2p_save_path
    }

    expobj.Suite2p.update_ops(new_ops)
    expobj.save()

    from packerlabimaging import import_obj
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


# ## creating test for Suite2p.s2pRun():
# # TODO need to add custom bad_frames creation function
# expobj: Experiment = existing_expobj_nopredones2p_fixture()
#
# expobj.Suite2p.bad_frames = []
# for trial in expobj.trialIDs:
#     if 'post' in expobj.TrialsInformation[trial]['expGroup']:
#         trialobj: TwoPhotonImagingTrial = expobj.load_trial(trialID=trial)
#         print(len(expobj.Suite2p.bad_frames))
#
#         expobj.Suite2p.add_bad_frames(frames=np.arange(trialobj.Suite2p.trial_frames[0], trialobj.Suite2p.trial_frames[1]), bad_frames_npy_loc = expobj.dataPath)
#
# expobj.save()
#
# expobj.Suite2p.s2pRun(expobj=expobj)


def test_makeFrameAverageTiff():
    from packerlabimaging.utils.io import import_obj
    expobj = import_obj(pkl_path='/home/pshah/mnt/qnap/Analysis/2021-01-25/PS12/PS12_analysis.pkl')
    trialobj = expobj.load_trial(trialID=expobj.trialIDs[0])

    trialobj.Suite2p.makeFrameAverageTiff(reg_tif_dir='/home/pshah/mnt/qnap/Analysis/2021-01-25/PS12/suite2p/reg_tif/',
                                          frames=[110, 510],
                                          peri_frames=100, save_dir=trialobj.saveDir + f'/export/avg_frames/',
                                          to_plot=True)

test_makeFrameAverageTiff()

