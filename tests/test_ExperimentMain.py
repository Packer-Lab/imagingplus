import packerlabimaging as pkg


# pytest framework
def test_ExperimentClass(experiment_fixture):
    # print(experiment_fixture)
    expobj = pkg.Experiment(**experiment_fixture)
    print(expobj)

# TODO write tests for public Experiment methods

