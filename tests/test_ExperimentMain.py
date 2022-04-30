import pytest

import packerlabimaging as pkg


# pytest framework
from packerlabimaging._archive.ExperimentMain import Experiment

@pytest.mark.skip
def test_ExperimentClass(experiment_fixture):
    # print(experiment_fixture)
    expobj = pkg.Experiment(**experiment_fixture)
    print(expobj)

# TODO need to run tests on new experiment/trial structure
def test_Experiment_new(experimentnew_fixture):
    expobj = Experiment(**experimentnew_fixture)
    print(expobj)

def test_add_suite2p():
    """test for adding suite2p trials without any results loaded."""
    expobj = pkg.import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    expobj.add_suite2p(s2p_trials='all')
    return expobj

def test_add_suite2p_results():
    """test for adding suite2p trials with precompleted results. """
    s2p_results_path = '/home/pshah/mnt/qnap/Analysis/2020-12-19/suite2p/alloptical-2p-1x-alltrials/plane0'

    expobj = pkg.import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
    # expobj.add_suite2p(s2p_trials='all', s2pResultsPath=s2p_results_path)
    expobj.add_suite2p(s2p_trials=['t-005'], s2pResultsPath=s2p_results_path)
    t005 = expobj.load_trial(trialID='t-005')

test_add_suite2p_results()


# TODO write tests for public Experiment methods - especially now for add_trial


