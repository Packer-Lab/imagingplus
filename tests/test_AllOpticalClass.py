import os

import imagingplus as pli
from conftest import alloptical_trial_fixture
from imagingplus import AllOpticalTrial

LOCAL_DATA_PATH = '/Users/prajayshah/data/oxford-data-to-process/'
REMOTE_DATA_PATH = '/home/pshah/mnt/qnap/Data/'
BASE_PATH = LOCAL_DATA_PATH


# %%


# %%

def test_AllOpticalClass(alloptical_trial_fixture, existing_expobj_nopredones2p_fixture):
    """
    TODO fill explaination add parameters
    :param alloptical_trial_fixture:
    :return:
    """
    from imagingplus.processing.microscopes import PrairieViewMetadata
    from imagingplus.processing.paq import PaqData

    paqs_loc = f'{BASE_PATH}/2020-12-19/2020-12-19_RL109_013.paq'  # path to the .paq files for the selected trials
    dataPath = alloptical_trial_fixture['dataPath']

    # parses imaging system cellsdata
    imparams = PrairieViewMetadata(pv_xml_dir=os.path.dirname(dataPath), microscope='Bruker 2pPlus')

    # sets the stim start frames
    tmdata = PaqData.paqProcessingAllOptical(paq_path=paqs_loc, frame_channel='frame_clock', stim_channel='markpoints2packio')

    # create the trial
    aotrial = AllOpticalTrial(imparams=imparams, tmdata=tmdata, **alloptical_trial_fixture)
    return aotrial

# idict = alloptical_trial_fixture()
# aotrial = test_AllOpticalClass(idict)


def test_processing_targets_stims(existing_trialobj_alloptical_fixture):
    self = existing_trialobj_alloptical_fixture
    self.twopstim, self.stim_duration_frames = self.photostimProcessing()

    self.raw_SLMTargets, self.dFF_SLMTargets, self.meanFluImg_registered = self.collectSignalFromCoords(
        curr_trial_frames=self.Suite2p.trial_frames, save=True)
    self.targets_dff, self.targets_dff_avg, self.targets_dfstdF, self.targets_dfstdF_avg, \
    self.targets_raw, self.targets_raw_avg = self.getTargetsStimTraceSnippets()
    self.save()

# %%
# import imagingplus as pkg
#
# expobj = pkg.import_obj('/home/pshah/Documents/code/imagingplus/tests/RL109_analysis.pkl')
# trial = expobj.load_trial('t-013')
#
# trial.cellsdata
#
