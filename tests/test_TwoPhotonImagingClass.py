import os

import pytest

from imagingplus.processing.paq import PaqData
from imagingplus.processing.microscopes import PrairieViewMetadata
from imagingplus.workflows.TwoPhotonImaging import TwoPhotonImaging



@pytest.mark.skip
def test_TwoPhotonImagingTrial(twophoton_imaging_multitrial_noPreDoneSuite2p_fixture, existing_expobj_fixture):
    # passing
    """test for TwoPhotonImaging workflow"""
    info = twophoton_imaging_multitrial_noPreDoneSuite2p_fixture
    expobj = existing_expobj_fixture
    date = info['date']
    prep = info['prep']
    trials_paq = info['trials_paq']

    for trial, paq in trials_paq.items():
        if trial == 't-001':
            data_path_base = f'/mnt/qnap_share/Data/imagingplus-example/test-data'

            paqs_loc = f'{data_path_base}/{date}_{prep}_{paq}'  # path to the .paq files for the selected trials
            dataPath = f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'

            imparams = PrairieViewMetadata(pv_xml_dir=os.path.dirname(dataPath), microscope='Bruker 2pPlus')
            tmdata = PaqData.paqProcessingTwoPhotonImaging(paq_path=paqs_loc, frame_channel='frame_clock')

            trialobj = TwoPhotonImaging(date=date, trialID= trial, expID= prep, imparams =  imparams, tmdata= tmdata,
                                        saveDir=f'/mnt/qnap_share/Data/imagingplus-example/imagingplus-test-analysis/',
                                        dataPath= dataPath, expGroup= "awake spont. 2p imaging + LFP")

            expobj.add_imaging_trial(trialID=trial, trialobj=trialobj)
        else:
            print('skipping remaining trials becaues those raw data have not been copied over to the qnap_share location.')

@pytest.mark.skip
def test_meanRawFluTrace(existing_trialobj_twophotonimaging_fixture):
    # passing
    trialobj: TwoPhotonImaging = existing_trialobj_twophotonimaging_fixture[0]
    trialobj.meanFluImg, trialobj.meanFovFluTrace = trialobj.meanRawFluTrace()
    trialobj.save()

# TODO add tests for main public functions - especially those in tutorial 1





# archive

# @pytest.mark.skip
# def test_TwoPhotonImagingTrial_new(twophoton_imaging_trial_new_noPreDoneSuite2p_fixture):
#     """archived test im pretty sure......"""
#     trials_paq = twophoton_imaging_trial_new_noPreDoneSuite2p_fixture
#
#     date = '2021-01-31'
#     prep = 'HF113'
#
#
#     for trial, paq in trials_paq.items():
#         data_path_base = f'/home/pshah/mnt/qnap/Data/{date}'  # same as above, unclear what this is supposed to be.. (need more clear argument name)
#
#         paqs_loc = f'{data_path_base}/{date}_{prep}_{paq}'  # path to the .paq files for the selected trials
#         dataPath = f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'
#
#         from imagingplus.processing.imagingMetadata import PrairieViewMetadata
#         imparams = PrairieViewMetadata(pv_xml_dir=os.path.dirname(dataPath), microscope='Bruker 2pPlus')
#         from imagingplus.main.paq import PaqData
#         tmdata = PaqData.paqProcessingTwoPhotonImaging(paq_path=paqs_loc, frame_channel='frame_clock')
#
#         imaging_info = {'date': date,
#                         'trialID': trial,
#                         'expID': prep,
#                         'saveDir': f'/home/pshah/mnt/qnap/Analysis/{date}/{prep}/',
#                         'dataPath': dataPath,
#                         'expGroup': "awake spont. 2p imaging + LFP",
#                         }
#
#         TwoPhotonImaging(imparams=imparams, tmdata=tmdata, **imaging_info)
