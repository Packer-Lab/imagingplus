import pytest

# pytest framework
from packerlabimaging.main.core import Experiment


# test passing
@pytest.mark.skip
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


def test_add_suite2p_results(existing_expobj_fixture):
    # todo run test
    """test for adding suite2p trials with precompleted results. """

    expobj, date = existing_expobj_fixture
    s2p_path = f'/home/pshah/mnt/qnap/Analysis/{date}/{expobj.expID}//suite2p//plane0/'

    expobj.add_suite2p(s2p_trials=[expobj.trialIDs[0]], s2pResultsPath=s2p_path)
    # t005 = expobj.load_trial(trialID='t-005')

# TODO write tests for public Experiment methods - especially now for add_trial
