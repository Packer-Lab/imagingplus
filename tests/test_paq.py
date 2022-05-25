from packerlabimaging.plotting.plotting import plot__paq_channel
from packerlabimaging._archive.paq import PaqData

from packerlabimaging._archive.TwoPhotonImagingMain import TwoPhotonImagingTrial

def test_import_paqdata(paqpath = '/home/pshah/mnt/qnap/Data/2021-01-25/2021-01-25_PS12_009.paq'):
    paq_data_obj, paqdata = PaqData.import_paqdata(paq_path=paqpath)

# test_import_paqdata()

def test_paq_func(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture
    trialobj.Paq.storePaqChannel(chan_name='voltage')

# test_paq_func(existing_trialobj_twophotonimaging_fixture)


def test_plot__paq_channel(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture
    plot__paq_channel(trialobj.Paq, channel='voltage', x_axis='Time (secs)', x_tick_secs=120)


