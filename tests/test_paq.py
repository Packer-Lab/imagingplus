from imagingplus import TwoPhotonImaging
from imagingplus.processing.paq import PaqData


def test_import_paqdata(paqpath = '/home/pshah/mnt/qnap/Data/2021-01-25/2021-01-25_PS12_009.paq'):
    paq_data_obj, paqdata = PaqData.import_paqdata(file_path=paqpath)


def test_cropData(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImaging = existing_trialobj_twophotonimaging_fixture
    trialobj.tmdata.cropData(begin=trialobj.tmdata.frame_times[0], end=trialobj.tmdata.frame_times[-1], channels='all')

