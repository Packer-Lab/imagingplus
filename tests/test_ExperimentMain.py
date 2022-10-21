import pytest

# pytest framework
import imagingplus as pli

from imagingplus.main.core import Experiment


# test passing
# @pytest.mark.skip
from imagingplus.utils.io import import_obj


def test_Experiment(experiment_fixture):
    expobj = Experiment(**experiment_fixture)
    print(expobj)


@pytest.mark.skip
def test_add_suite2p(existing_expobj_nopredones2p_fixture):
    # passing
    """test for adding suite2p trials without any results loaded."""
    expobj: Experiment = existing_expobj_nopredones2p_fixture
    expobj.add_suite2p(s2p_trials='all')
    return expobj


# @pytest.mark.skip
def test_add_suite2p_results(existing_expobj_fixture):
    # test passing
    """test for adding suite2p trials with precompleted results. """
    expobj = existing_expobj_fixture
    s2pResultsPath="/home/pshah/mnt/qnap/Analysis/2021-01-31/HF113/suite2p/plane0"

    # expobj.add_suite2p(s2p_trials=['t-005', 't-006', 't-013'], s2pResultsPath=s2pResultsPath)

    # # temp fix ///
    # s2p_trials = expobj.trialIDs
    # for trial in s2p_trials: expobj.TrialsInformation[trial]['paths']['dataPath'] = expobj.TrialsInformation[trial]['paths']['data_path']
    # expobj.save()
    expobj.add_suite2p(s2p_trials='all', s2pResultsPath=s2pResultsPath)

