from packerlabimaging.TwoPhotonImagingMain import TwoPhotonImagingTrial

def test_paq_func(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture
    trialobj.Paq.storePaqChannel(chan_name='voltage')

# test_paq_func(existing_trialobj_twophotonimaging_fixture)