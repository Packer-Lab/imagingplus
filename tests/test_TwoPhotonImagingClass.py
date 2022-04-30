import pytest

from conftest import twophoton_imaging_trial_new_noPreDoneSuite2p_fixture
from packerlabimaging._archive import TwoPhotonImagingMain
from packerlabimaging.workflows.TwoPhotonImaging import TwoPhotonImaging

@pytest.mark.skip
def test_TwoPhotonImagingTrial_new(twophoton_imaging_trial_new_noPreDoneSuite2p_fixture):
    dict_ = twophoton_imaging_trial_new_noPreDoneSuite2p_fixture()
    return TwoPhotonImagingMain.TwoPhotonImagingTrial(**dict_['t-005'])

# trialobj = test_TwoPhotonImagingTrial_new(twophoton_imaging_trial_new_noPreDoneSuite2p_fixture)


def test_TwoPhotonImagingTrial(twophoton_imaging_multitrial_noPreDoneSuite2p_fixture):
    """test for TwoPhotonImaging workflow"""
    trials_ = twophoton_imaging_multitrial_noPreDoneSuite2p_fixture()
    for trial, info in trials_.items()[0]:
        TwoPhotonImaging(**info)

@pytest.mark.skip
def test_meanRawFluTrace(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingMain.TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    trialobj.meanFluImg, trialobj.meanFovFluTrace = trialobj.meanRawFluTrace()
    trialobj.save()

# TODO add tests for main public functions - especially those in tutorial 1




