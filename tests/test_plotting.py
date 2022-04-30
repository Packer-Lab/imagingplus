from packerlabimaging.workflows.AllOptical import AllOpticalTrial

from packerlabimaging._archive.TwoPhotonImagingMain import TwoPhotonImagingTrial

from packerlabimaging._archive.ExperimentMain import Experiment

from packerlabimaging.plotting.plotting import makeSuite2pPlots, plot_flu_trace, plotMeanFovFluTrace, \
    plot_photostim_traces_overlap, plot_periphotostim_avg, plotRoiLocations, plot_SLMtargets_Locs, makeAverageTiff

# expobj: Experiment = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
# trialobj: TwoPhotonImagingTrial = expobj.load_trial(trialID='t-005')


def test_plotRoiLocations(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    plotRoiLocations(trialobj=trialobj, suite2p_rois='all', background = None)
    plotRoiLocations(trialobj=trialobj, suite2p_rois='all', background = None, scalebar=True)


def test_makeSuite2pPlots(existing_expobj_fixture, existing_trialobj_twophotonimaging_fixture):
    expobj: Experiment = existing_expobj_fixture[0]
    makeSuite2pPlots(obj=expobj)

    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    makeSuite2pPlots(obj=trialobj, scalebar=True)


def test_plot_flu_trace(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture[0]
    plot_flu_trace(trialobj=trialobj, cell=10, to_plot = 'raw', linewidth = 0.10,
                            x_lims=None, y_lims= None)

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

def test_makeAverageTiff(tiff_path, save_path):
    makeAverageTiff(tiff_path=tiff_path, save_path=save_path)

test_makeAverageTiff(tiff_path='/home/pshah/mnt/qnap/Data/2020-03a/2020-03-03/2020-03-03_t-001/2020-03-03_t-001_Cycle00001_Ch3_downsampled.tif', save_path='/home/pshah/mnt/qnap/Analysis/analysis_export/2020-03-03_t-001_Cycle00001_Ch3_downsampled_MEAN.tif')