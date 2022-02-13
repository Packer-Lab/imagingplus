from packerlabimaging.TwoPhotonImagingMain import TwoPhotonImagingTrial

from packerlabimaging import Experiment

from packerlabimaging.processing import suite2p


def test_Suite2pResultsExperiment(existing_expobj_fixture):
    expobj: Experiment = existing_expobj_fixture[0]
    trialSuite2p = existing_expobj_fixture[1]
    s2pResultsPath = existing_expobj_fixture[2]

    expobj.Suite2p = suite2p.Suite2pResultsExperiment(trialsSuite2p = trialSuite2p,
                                                      s2pResultsPath = s2pResultsPath)

def test_Suite2pResultsTrial(existing_trialobj_twophotonimaging_fixture, existing_trialobj_alloptical_fixture, existing_expobj_fixture):
    trialobj, trialobj_ = existing_trialobj_twophotonimaging_fixture
    alloptical_trialobj = existing_trialobj_alloptical_fixture
    expobj: Experiment = existing_expobj_fixture[0]

    for n_obj in [trialobj, trialobj_, alloptical_trialobj]:
        n_obj.Suite2p = suite2p.Suite2pResultsTrial(suite2p_experiment_obj = expobj.Suite2p, trial_frames = trialobj.Suite2p.trial_frames)
        n_obj.save()
