import pytest

import packerlabimaging as pkg


# pytest framework
from conftest import experimentnew_fixture
from packerlabimaging.ExperimentMain import Experiment

@pytest.mark.skip
def test_ExperimentClass(experiment_fixture):
    # print(experiment_fixture)
    expobj = pkg.Experiment(**experiment_fixture)
    print(expobj)

# TODO need to run tests on new experiment/trial structure
def test_Experiment_new(experimentnew_fixture):
    expobj = Experiment(**experimentnew_fixture)
    print(expobj)


# TODO write tests for public Experiment methods - especially now for add_trial and add_suite2p_experiment


