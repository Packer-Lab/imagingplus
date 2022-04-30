import os

import pytest

from conftest import twophoton_imaging_trial_new_noPreDoneSuite2p_fixture
from packerlabimaging._archive import TwoPhotonImagingMain
from packerlabimaging.main.paq import PaqData
from packerlabimaging.processing.imagingMetadata import PrairieViewMetadata
from packerlabimaging.workflows.TwoPhotonImaging import TwoPhotonImaging

@pytest.mark.skip
def test_TwoPhotonImagingTrial_new(twophoton_imaging_trial_new_noPreDoneSuite2p_fixture):
    dict_ = twophoton_imaging_trial_new_noPreDoneSuite2p_fixture()
    return TwoPhotonImagingMain.TwoPhotonImagingTrial(**dict_['t-005'])

# trialobj = test_TwoPhotonImagingTrial_new(twophoton_imaging_trial_new_noPreDoneSuite2p_fixture)


def test_TwoPhotonImagingTrial(twophoton_imaging_multitrial_noPreDoneSuite2p_fixture):
    """test for TwoPhotonImaging workflow"""
    info = twophoton_imaging_multitrial_noPreDoneSuite2p_fixture

    date = info['date']
    prep = info['prep']
    trials_paq = info['trials_paq']

    for trial, paq in trials_paq.items():
        data_path_base = f'/home/pshah/mnt/qnap/Data/{date}'  # same as above, unclear what this is supposed to be.. (need more clear argument name)

        paqs_loc = f'{data_path_base}/{date}_{prep}_{paq}'  # path to the .paq files for the selected trials
        dataPath = f'{data_path_base}/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'

        imparams = PrairieViewMetadata(pv_xml_dir=os.path.dirname(dataPath), microscope='Bruker 2pPlus')
        tmdata = PaqData.paqProcessingTwoPhotonImaging(paq_path=paqs_loc, frame_channel='frame_clock')

        TwoPhotonImaging(date=date, trialID= trial, expID= prep, imparams =  imparams, tmdata= tmdata,
                         saveDir=f'/home/pshah/mnt/qnap/Analysis/{date}/{prep}/', dataPath= dataPath, expGroup= "awake spont. 2p imaging + LFP"
        )



@pytest.mark.skip
def test_meanRawFluTrace(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingMain.TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    trialobj.meanFluImg, trialobj.meanFovFluTrace = trialobj.meanRawFluTrace()
    trialobj.save()

# TODO add tests for main public functions - especially those in tutorial 1




