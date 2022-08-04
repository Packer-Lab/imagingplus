import pytest

# pytest framework
import packerlabimaging as pli

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

    expobj = existing_expobj_fixture
    s2pResultsPath="/home/pshah/mnt/qnap/Analysis/2020-12-19/suite2p/alloptical-2p-1x-alltrials/plane0/"

    # expobj.add_suite2p(s2p_trials=['t-005', 't-006', 't-013'], s2pResultsPath=s2pResultsPath)
    expobj.add_suite2p(s2p_trials=['t-013'], s2pResultsPath=s2pResultsPath)


expobj = pli.import_obj(pkl_path='/mnt/qnap_share/Data/packerlabimaging-example/RL109_analysis.pkl')

test_add_suite2p_results(expobj)

# TODO write tests for public Experiment methods - especially now for add_trial
