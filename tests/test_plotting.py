import numpy as np
from packerlabimaging._io import import_obj

from packerlabimaging.TwoPhotonImagingMain import TwoPhotonImagingTrial

from packerlabimaging.ExperimentMain import Experiment

import packerlabimaging.plotting as plotting
from packerlabimaging.plotting._utils import heatmap_options, image_frame_options
from packerlabimaging.plotting.plotting import makeSuite2pPlots


def test_makeSuite2pPlots(existing_expobj_twophotonimaging_fixture):
    expobj: Experiment = existing_expobj_twophotonimaging_fixture[0]
    plotting.makeSuite2pPlots(expobj)

# expobj: Experiment = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
# makeSuite2pPlots(expobj)

def test_plot_flu_trace(existing_trialobj_twophotonimaging_fixture):
    trialobj: TwoPhotonImagingTrial = existing_trialobj_twophotonimaging_fixture
    plotting.plot_flu_trace(trialobj, cell=np.random.choice(len(trialobj.Suite2p.cell_id)), to_plot = 'raw', linewidth = 0.10,
                            x_lims=None, y_lims= None)




# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import numpy as np
#
# mpl.rcParams.update({
#     'axes.spines.left': False,
#     'axes.spines.bottom': False,
#     'axes.spines.top': False,
#     'axes.spines.right': False,
#     'legend.frameon': False,
#     'figure.subplot.wspace': .01,
#     'figure.subplot.hspace': .01,
#     'ytick.major.left': False,
#     'xtick.major.bottom': False
# })
#
# expobj: Experiment = import_obj(pkl_path='/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
# obj = expobj
#
# image_frame_options()
# heatmap_options()
#
# plt.figure(figsize=[10, 3])
# plt.subplot(1, 4, 1)
# plt.imshow(obj.Suite2p.output_op['max_proj'], cmap='gray')
# plt.title("Registered Image, Max Projection", wrap=True)
#
# plt.subplot(1, 4, 2)
# plt.imshow(np.nanmax(obj.Suite2p.im, axis=0), cmap='jet')
# plt.title("All ROIs Found")
#
#
#
# plt.subplot(1, 4, 3)
# plt.imshow(np.nanmax(obj.Suite2p.im[~obj.Suite2p.iscell], axis=0, ), cmap='jet')
# plt.title("All Non-Cell ROIs")
#
# plt.subplot(1, 4, 4)
# plt.imshow(np.nanmax(obj.Suite2p.im[obj.Suite2p.iscell], axis=0), cmap='jet')
# plt.title("All Cell ROIs")
#
# plt.show()

