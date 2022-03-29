from packerlabimaging import TwoPhotonImagingMain


def test_TwoPhotonImagingTrial(twophoton_imaging_trial_fixture):
    TwoPhotonImagingMain.TwoPhotonImagingTrial(**twophoton_imaging_trial_fixture)


def test_meanRawFluTrace(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingMain.TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    trialobj.meanFluImg, trialobj.meanFovFluTrace = trialobj.meanRawFluTrace()
    trialobj.save()

# TODO add tests for main public functions - especially those in tutorial 1