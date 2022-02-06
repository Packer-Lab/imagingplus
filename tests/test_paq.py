##
from packerlabimaging import Experiment

import packerlabimaging as pkg
from packerlabimaging.TwoPhotonImagingMain import TwoPhotonImagingTrial
from conftest import existing_trialobj_twophotonimaging_fixture

# trying to use pytest -- seems to work (!)
def test_paq_func(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture
    trialobj.paq.storePaqChannel(chan_name='voltage')

# test_paq_func(existing_trialobj_twophotonimaging_fixture)

# pytest fixture not working for whatever reason
# def test_paq_func():
#     expobj: Experiment = pkg.import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    trialobj: TwoPhotonImagingTrial = pkg.import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/2020-12-19_t-013.pkl')
#     print(expobj.trialsInformation['t-005']['analysis_object_information']['pkl path'])
#     trialobj = expobj.load_trial(trialID='t-005')
#
#     trialobj.paq.storePaqChannel(chan_name='voltage')
