import numpy as np
from packerlabimaging.AllOpticalMain import AllOpticalTrial

from packerlabimaging._io import import_obj

from packerlabimaging.TwoPhotonImagingMain import TwoPhotonImagingTrial

from packerlabimaging.ExperimentMain import Experiment

import packerlabimaging.plotting as plotting
from packerlabimaging.plotting._utils import heatmap_options, image_frame_options
from packerlabimaging.plotting.plotting import makeSuite2pPlots, plot_flu_trace, plotMeanFovFluTrace, \
    plot_photostim_traces_overlap

expobj: Experiment = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
trialobj: TwoPhotonImagingTrial = expobj.load_trial(trialID='t-005')


def test_makeSuite2pPlots(existing_expobj_fixture):
    expobj: Experiment = existing_expobj_fixture[0]
    makeSuite2pPlots(expobj)


def test_plot_flu_trace(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    plot_flu_trace(trialobj=trialobj, cell=10, to_plot = 'raw', linewidth = 0.10,
                            x_lims=None, y_lims= None)

def test_plotMeanFovFluTrace(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    plotMeanFovFluTrace(trialobj=trialobj)


def test_plot_photostim_traces_overlap(existing_trialobj_alloptical_fixture):
    trialobj: AllOpticalTrial = existing_trialobj_alloptical_fixture
    plot_photostim_traces_overlap()
