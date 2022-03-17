from packerlabimaging.TwoPhotonImagingMain import TwoPhotonImagingTrial

from packerlabimaging import Experiment

from packerlabimaging.processing import suite2p
from packerlabimaging.processing.suite2p import Suite2pResultsExperiment, Suite2pResultsTrial


def test_Suite2pResultsExperiment(existing_expobj_fixture):
    expobj: Experiment = existing_expobj_fixture[0]
    trialSuite2p = existing_expobj_fixture[1]
    s2pResultsPath = existing_expobj_fixture[2]

    expobj.Suite2p = suite2p.Suite2pResultsExperiment(trialsSuite2p=trialSuite2p,
                                                      s2pResultsPath=s2pResultsPath)


def test_Suite2pResultsTrial(existing_trialobj_twophotonimaging_fixture, existing_trialobj_alloptical_fixture,
                             existing_expobj_fixture):
    trialobj, trialobj_ = existing_trialobj_twophotonimaging_fixture
    alloptical_trialobj = existing_trialobj_alloptical_fixture
    expobj: Experiment = existing_expobj_fixture[0]

    for n_obj in [trialobj, trialobj_, alloptical_trialobj]:
        from packerlabimaging.processing.suite2p import Suite2pResultsExperiment
        s2p_expobj: Suite2pResultsExperiment = expobj.Suite2p
        n_obj.Suite2p = suite2p.Suite2pResultsTrial(trialsSuite2p=s2p_expobj.trials, s2pResultsPath=s2p_expobj.s2pResultsPath,
                                                    subtract_neuropil=s2p_expobj.subtract_neuropil,
                                                    trial_frames=n_obj.Suite2p.trial_frames)  # use trial obj's current trial frames

        n_obj.save()

def test_Suite2pExp():
    # TODO run test for suite2p without predone results path
    _trialsSuite2p = ['t-005', 't-006']
    # trialIDs = [*twophoton_imaging_trial_fixture_noPreDoneSuite2p['TrialsInformation']]
    # for trial in trialIDs:
    #         assert 's2p_use' in [*twophoton_imaging_trial_fixture_noPreDoneSuite2p.TrialsInformation[trial]], 'when trying to utilize suite2p , must provide value for `s2p_use` ' \
    #                      'in TrialsInformation[trial] for each trial to specify if to use trial for this suite2p associated with this experiment'
    #         _trialsSuite2p.append(trial) if twophoton_imaging_trial_fixture_noPreDoneSuite2p['TrialsInformation'][trial]['s2p_use'] else None

    Suite2p = Suite2pResultsExperiment(trialsSuite2p=_trialsSuite2p)


# test_Suite2pExp()


def test_Suite2pTrial():
    # TODO run tests for suite2p trial creation without predone results path
    trialsSuite2p: list = ['t-005', 't-006']
    trial_frames: tuple
    s2ptrial = Suite2pResultsTrial(trialsSuite2p=trialsSuite2p, trial_frames=(0, 1000))
    assert s2ptrial._s2pResultExists == False  # - doesnt seem to be setting to False as expected during the super() call!

# test_Suite2pTrial()

