from packerlabimaging.TwoPhotonImagingMain import TwoPhotonImagingTrial

from packerlabimaging import Experiment

from packerlabimaging.processing import suite2p


def test_Suite2pResultsExperiment(existing_expobj_twophotonimaging_fixture):
    expobj: Experiment = existing_expobj_twophotonimaging_fixture[0]
    trialSuite2p = existing_expobj_twophotonimaging_fixture[1]
    s2pResultsPath = existing_expobj_twophotonimaging_fixture[2]

    expobj.Suite2p = suite2p.Suite2pResultsExperiment(trialsSuite2p = trialSuite2p,
                                                      s2pResultsPath = s2pResultsPath)
    # expobj.save()

def test_Suite2pResultsTrial(existing_trialobj_twophotonimaging_fixture, existing_expobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture
    expobj: Experiment = existing_expobj_twophotonimaging_fixture[0]

    n_obj = trialobj
    n_obj.Suite2p = suite2p.Suite2pResultsTrial(suite2p_experiment_obj = expobj.Suite2p, trial_frames = trialobj.Suite2p.trial_frames)
    n_obj.save()
