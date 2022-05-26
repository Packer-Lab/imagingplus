import pytest

from packerlabimaging.workflows.AllOptical import AllOpticalTrial

from packerlabimaging._archive.TwoPhotonImagingMain import TwoPhotonImagingTrial

from packerlabimaging._archive.ExperimentMain import Experiment

from packerlabimaging.plotting.plotting import makeSuite2pPlots, plot_flu_trace, plotMeanFovFluTrace, \
    plot_photostim_traces_overlap, plot_periphotostim_avg, plotRoiLocations, plot_SLMtargets_Locs, MeanProject, \
    SingleFrame, FrameAverage


# expobj: Experiment = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
# trialobj: TwoPhotonImagingTrial = expobj.load_trial(trialID='t-005')




def test_plotRoiLocations(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    plotRoiLocations(trialobj=trialobj, suite2p_rois='all', background=None)
    plotRoiLocations(trialobj=trialobj, suite2p_rois='all', background=None, scalebar=True)


def test_makeSuite2pPlots(existing_expobj_fixture, existing_trialobj_twophotonimaging_fixture):
    expobj: Experiment = existing_expobj_fixture[0]
    makeSuite2pPlots(obj=expobj)

    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    makeSuite2pPlots(obj=trialobj, scalebar=True)


def test_plot_flu_trace(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    plot_flu_trace(trialobj=trialobj, cell=10, to_plot='raw', linewidth=0.10,
                   x_lims=None, y_lims=None)


def test_plotMeanFovFluTrace(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    plotMeanFovFluTrace(trialobj=trialobj)


def test_plot_photostim_traces_overlap(existing_trialobj_alloptical_fixture):
    trialobj: AllOpticalTrial = existing_trialobj_alloptical_fixture
    plot_photostim_traces_overlap(array=trialobj.targets_dff_avg, trialobj=trialobj)


def test_plot_periphotostim_avg(existing_trialobj_alloptical_fixture):
    trialobj: AllOpticalTrial = existing_trialobj_alloptical_fixture
    plot_periphotostim_avg(arr=trialobj.targets_dff_avg, trialobj=trialobj)


def test_plot_SLMtargets_Locs(existing_trialobj_alloptical_fixture):
    trialobj: AllOpticalTrial = existing_trialobj_alloptical_fixture
    plot_SLMtargets_Locs(trialobj=trialobj)

# @pytest.skip('pass')
def test_FrameAverageTiff(tiff_path_fixture, existing_trialobj_twophotonimaging_fixture):
    # FrameAverage(key_frames=[500, 1000, 1500, 2000, 3000], tiff_path=tiff_path_fixture, peri_frames=100, plot=True,
    #              scalebar_um=100, trialobj=existing_trialobj_twophotonimaging_fixture)
    FrameAverage(key_frames=27520, tiff_path=tiff_path_fixture, peri_frames=100, plot=True, title='custom title')


def test_SingleTiffFrame(s_tiff_path_fixture, existing_trialobj_twophotonimaging_fixture):
    SingleFrame(tiff_path=s_tiff_path_fixture, trialobj=existing_trialobj_twophotonimaging_fixture, scalebar_um=100)


