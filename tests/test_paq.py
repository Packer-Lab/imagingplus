from packerlabimaging.processing.paq import PaqData

from packerlabimaging.TwoPhotonImagingMain import TwoPhotonImagingTrial

def test_import_paqdata(paqpath = '/home/pshah/mnt/qnap/Data/2021-01-25/2021-01-25_PS12_009.paq'):
    paq_data_obj, paqdata = PaqData.import_paqdata(paq_path=paqpath)

# test_import_paqdata()

def test_paq_func(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture
    trialobj.Paq.storePaqChannel(chan_name='voltage')

# test_paq_func(existing_trialobj_twophotonimaging_fixture)