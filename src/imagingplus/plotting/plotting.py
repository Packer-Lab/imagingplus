# library of convenience plotting funcs that are used for making various plots for all optical photostimulation/imaging experiments

# imports
import os
from typing import Union

import mpl_point_clicker
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import tifffile as tf
from imagingplus.workflows.TwoPhotonImaging import TwoPhotonImaging

from imagingplus import Experiment
from imagingplus.main.subcore import TemporalData
from imagingplus.utils.utils import save_to_csv
from imagingplus.utils.images import ImportTiff

from imagingplus.workflows.AllOptical import AllOpticalTrial

from imagingplus.plotting._utils import plotting_decorator, make_random_color_array, _add_scalebar, \
    image_frame_options, dataplot_frame_options, dataplot_ax_options, plot_coordinates, heatmap_options, image_frame_ops
from imagingplus.utils.classes import ObjectClassError


# DATA ANALYSIS PLOTTING FUNCS

@plotting_decorator(figsize=(6, 6))
def plotImg(img: np.ndarray, **kwargs):
    """Plot image in grayscale.
    
    :param img: input image to show
    :param kwargs:
        :trialobj: ImagingTrial or SingleImage; object associated with input image.
        :scalebar_um: int; size of scalebar to plot on image (in um); must provide trialobj parameter.
    """
    assert img.ndim == 2, 'img to plot must only have 2 dimensions.'
    fig, ax = kwargs['fig'], kwargs['ax']
    ax.imshow(img, cmap='gray')
    if 'scalebar_um' in kwargs:
        fig.tight_layout(pad=0.2)
        if 'trialobj' in kwargs:
            _add_scalebar(ax=ax, **kwargs)
        else:
            raise ValueError('must provide trialobj parameter to make scalebar.')


# suite2p cellsdata
# simple plot of the location of the given cell(s) against a black FOV
@mpl.rc_context(image_frame_ops)
@plotting_decorator(figsize=(5, 5))
def plotRoiLocations(trialobj: TwoPhotonImaging, suite2p_rois: Union[list, str] = 'all',
                     background: np.ndarray = None, **kwargs):
    """
    plots an image of the FOV to show the locations of cells given in cells ls.

    :param trialobj: alloptical or 2p imaging object
    :param background: either 2dim numpy array to use as the backsplash or None (default; which is a black background)
    :param suite2p_rois: list of ROIs (suite2p cell IDs) to show coord location, default is 'all' (which is all ROIs)
    :param edgecolor: specify edgecolor of the scatter plot for cells
    :param facecolor: specify facecolor of the scatter plot for cells
    :param cells: list of cells to plot
    :param title: title for plot
    :param color_float_list: if given, it will be used to color the cells according a colormap
    :param cmap: cmap to be used in conjuction with the color_float_array argument
    :param show_s2p_targets: if True, then will prioritize coloring of cell points based on whether they were photostim targets
    :param invert_y: if True, invert the reverse the direction of the y axis
    :param fig: a matplotlib.Figure instance, if provided use this fig for plotting
    :param ax: a matplotlib.Axes.axes instance, if provided use this ax for plotting
    :param show: if False, do not display plot (used when the necessity is to return the fig and ax objects to futher manipulation)
    :returns fig, ax: returns fig and ax object, if show is False
    """
    # image_frame_options()

    # fig = kwargs['fig']
    # suptitle = kwargs['suptitle'] if 'suptitle' in kwargs else None
    image_frame_options() if 'apply_image_frame_options' not in kwargs or kwargs[
        'apply_image_frame_options'] else None

    ax = kwargs['ax']
    kwargs.pop('ax')

    facecolors = kwargs['facecolors'] if 'facecolors' in kwargs else 'none'
    edgecolors = kwargs['edgecolors'] if 'edgecolors' in kwargs else 'orange'

    if suite2p_rois == 'all':
        suite2p_rois = trialobj.Suite2p.cell_id

    for cell in suite2p_rois:
        y, x = trialobj.Suite2p.stat[trialobj.Suite2p.cell_id.index(cell)]['med']
        ax.scatter(x=x, y=y, edgecolors=edgecolors, facecolors=facecolors, linewidths=0.8)

    if background is None:
        black = np.zeros((trialobj.imparams.frame_x, trialobj.imparams.frame_y), dtype='uint16')
        ax.imshow(black, cmap='Greys_r', zorder=0)
        ax.set_xlim(0, trialobj.imparams.frame_x)
        ax.set_ylim(0, trialobj.imparams.frame_y)
    else:
        ax.imshow(background, cmap='Greys_r', zorder=0)

    ax.invert_yaxis()

    _add_scalebar(trialobj=trialobj, ax=ax) if 'scalebar' in kwargs and kwargs['scalebar'] is True else None
    mpl.pyplot.rcdefaults()


@mpl.rc_context(image_frame_ops)
@plotting_decorator(figsize=(15, 5), nrows=1, ncols=4)
def makeSuite2pPlots(obj: Union[Experiment, TwoPhotonImaging], **kwargs):
    """
    Makes four plots that are created by Suite2p output.

    credit: run_s2p tutorial on Mouseland/Suite2p on github

    :param obj: Experiment or Trial object to use for creating Suite2p plots.
    :param axs: a matplotlib.Axes.axes instance (4 cols x 1 rows), if provided use these axes for plotting
    :param show: if False, do not display plot (used when the necessity is to return the fig and ax objects to futher manipulation)
    :returns fig, axs: returns fig and ax object, if show is False

    """

    heatmap_options()

    axs = kwargs['axs']
    kwargs.pop('axs')

    # f, axs = plt.subplots(figsize=[15, 5], nrows=1, ncols=4)

    # plt.subplot(1, 4, 1)
    axs[0].imshow(obj.Suite2p.output_ops['meanImgE'], cmap='gray')
    axs[0].set_title("Registered Image, Mean Enhanced", wrap=True)

    # plt.subplot(1, 4, 2)
    axs[1].imshow(np.nanmax(obj.Suite2p.im, axis=0), cmap='jet')
    axs[1].set_title("All ROIs Found", wrap=True)

    # plt.subplot(1, 4, 3)
    axs[2].imshow(np.nanmax(obj.Suite2p.im[~obj.Suite2p.iscell], axis=0, ), cmap='jet')
    axs[2].set_title("All Non-Cell ROIs", wrap=True)

    # plt.subplot(1, 4, 4)
    axs[3].imshow(np.nanmax(obj.Suite2p.im[obj.Suite2p.iscell], axis=0), cmap='jet')
    axs[3].set_title("All Cell ROIs", wrap=True)

    _add_scalebar(trialobj=obj, ax=axs[3]) if 'scalebar' in kwargs and kwargs['scalebar'] is True else None
    mpl.pyplot.rcdefaults()


@plotting_decorator(figsize=(20, 3))
def plot_flu_trace(trialobj: TwoPhotonImaging, cell: int, to_plot='raw', **kwargs):
    """
    Plot individual cell's flu or dFF trace, with photostim. timings for that cell. Accesses the object's anndata data storage.

    :param trialobj: TwoPhotonImaging (or derived) object containing anndata data storage from which to select cell to plot flu trace
    :param cell: cell index number
    :param to_plot: select data layer to plot (choose 'raw' to plot main raw data layer, otherwise choose name of the desired layer)
    :param fig: a matplotlib.Figure instance, if provided use this fig for plotting
    :param ax: a matplotlib.Axes.axes instance, if provided use this ax for plotting
    :param show: if False, do not display plot (used when the necessity is to return the fig and ax objects to futher manipulation)
    :returns fig, ax: returns fig and ax object, if show is False
    """
    dataplot_frame_options()
    ax = kwargs['ax']
    kwargs.pop('ax')
    for i in ['lw', 'linewidth']:
        try:
            lw = kwargs[i]
        except KeyError:
            lw = 1

    idx = trialobj.Suite2p.cell_id.index(cell)

    if to_plot == 'raw':
        data_to_plot = trialobj.data.X[idx, :]
    else:
        if to_plot in [*trialobj.data.layers]:
            data_to_plot = trialobj.data.layers[to_plot][idx, :]
        else:
            raise KeyError(f"to_plot processed cellsdata is not found in trialobj.cellsdata.layers. ")

    # make the plot either as just the raw trace or as a dFF trace with the std threshold line drawn as well.
    ax.plot(data_to_plot, linewidth=lw)

    dataplot_ax_options(ax=ax, **kwargs)
    mpl.pyplot.rcdefaults()


@plotting_decorator()
def plot__channel(tmdata: TemporalData, channel: str, **kwargs):
    """
    Plot the saved signal from the specified channel from a PaqData submodule.

    :param tmdata: TemporalData object
    :param channel: channel to plot
    :param kwargs:
        :x_axis: str; x axis label, if 'time' or "Time" found in x_axis label, will convert x axis to time domain.
        :x_tick_secs: int; interval to plot x axis ticks
        :ax: matplotlib axis object to use for plotting.
    """
    assert channel in tmdata.channels, f'{channel} not found in .Paq module cellsdata.'

    # set any kwargs provided
    ax = kwargs['ax']
    kwargs.pop('ax')

    kwargs['x_tick_secs'] = 120 if 'x_tick_secs' not in kwargs else kwargs['x_tick_secs']

    lw = 0.5 if 'lw' not in kwargs else kwargs['lw']
    color = 'black' if 'color' not in kwargs else kwargs['color']

    # collect cellsdata to plot
    data = tmdata.data[channel].to_numpy()

    # make plot
    ax.plot(data, lw=lw, color=color)
    ax.set_title(f"{channel}")

    ax.set_xlabel('tmdata clock')

    # set axis options
    dataplot_ax_options(ax=ax, collection_hz=tmdata.sampling_rate, **kwargs)
    kwargs['ax'] = ax

def MeanProject(tiff_path: str = None, frames: tuple = None, save_path: str = None, plot=True,
                imstack: np.ndarray = None, **kwargs) -> np.ndarray:
    """
    Mean Projection function.
    Creates, plots and (optionally) saves a Mean Projection image of the loaded tiff. If frames are provided, tiff is cropped to those frames.

    Note: saves as 8-bit grayscale image.

    :param tiff_path: path to tiff file to load (if imstack not provided)
    :param save_path: path to save output image
    :param plot: if True, show the output image
    :param imstack: data stack to use as input (if tiff_path not provided)
    :param frames: crop input stack to specified frames before calculating mean projection.
    :param kwargs:
        :trialobj: ImagingTrial or SingleImage; object associated with input image.
        :scalebar_um: int; size of scalebar to plot on image (in um); must provide trialobj parameter.
    :return: mean projection image

    """

    if imstack is None and tiff_path is not None:
        # read tiff
        print(f'\t\- Creating mean projection, from tiff: {tiff_path}')
        im_stack = ImportTiff(tiff_path=tiff_path, frames=frames)

        # if frames:
        #     im_stack = tf.imread(tiff_path, key=range(frames[0], frames[1]))
        # else:
        #     im_stack = tf.imread(tiff_path)
    elif imstack is not None and tiff_path is None:
        im_stack = imstack
        assert im_stack.ndim == 3, 'can only project 3D image stacks (implicit dimensions are: Frames x Xpixels x Ypixels)'
    else:
        raise ValueError(
            'values provided for both tiff_path and imstack. Unclear where to source image cellsdata from. Provide only one please.')

    if frames:
        print(f'\t collecting average image for frames: {frames}')
        im_stack = im_stack[frames[0]: frames[1]]

    print(f'\t\- Averaging {len(im_stack)} total frames.')
    img = np.mean(im_stack, axis=0)

    # convert to 8-bit
    from imagingplus.utils.images import convert_to_8bit
    img = convert_to_8bit(img, 0, 255)

    if save_path:
        if '.tif' in save_path:
            save_dir = os.path.dirname(save_path) + '/'
            os.makedirs(save_dir, exist_ok=True)
        else:
            raise ValueError(f'please provide full .tif path to save to. provided value was: {tiff_path}')
        print(f"\t\- Writing tiff to: {save_path}")
        tf.imwrite(save_path, img, photometric='minisblack')

    if plot:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(img, cmap='gray')
        fig.suptitle(f'Mean Project')
        fig.show()  # just plot for now to make sure that you are doing things correctly so far
        if 'scalebar_um' in kwargs and kwargs['scalebar_um']:
            if 'trialobj' in kwargs:
                _add_scalebar(trialobj=kwargs['trialobj'], ax=ax, **kwargs)

    return img


def MaxProject(tiff_path: str = None, frames: tuple = None, save_path: str = None, plot=True,
               imstack: np.ndarray = None, **kwargs) -> np.ndarray:
    """
    Maximum projection function.
    Creates, plots and (optionally) saves an average image of the loaded tiff. If frames are provided, tiff is cropped to those frames.

    Note: saves as 8-bit grayscale image.

    :param tiff_path: path to tiff file to load (if imstack not provided)
    :param save_path: path to save output image
    :param plot: if True, show the output image
    :param imstack: data stack to use as input (if tiff_path not provided)
    :param frames: crop input stack to specified frames before calculating maximum projection.
    :param kwargs:
        :trialobj: ImagingTrial or SingleImage; object associated with input image.
        :scalebar_um: int; size of scalebar to plot on image (in um); must provide trialobj parameter.
    :return: maximum projection image

    """

    if imstack is None and tiff_path is not None:
        # read tiff
        print(f'\t\- Creating max projection, from tiff: {tiff_path}')
        im_stack = ImportTiff(tiff_path=tiff_path, frames=frames)

        # if frames:
        #     im_stack = tf.imread(tiff_path, key=range(frames[0], frames[1]))
        # else:
        #     im_stack = tf.imread(tiff_path)

    elif imstack is not None and tiff_path is None:
        im_stack = imstack
        assert im_stack.ndim == 3, 'can only project 3D image stacks (implicit dimensions are: Frames x Xpixels x Ypixels)'
    else:
        raise ValueError(
            'values provided for both tiff_path and imstack. Unclear where to source image cellsdata from. Provide only one please.')

    if frames:
        print(f'\t collecting average image for frames: {frames}')
        im_stack = im_stack[frames[0]: frames[1]]

    print(f'\t\- Max Project {len(im_stack)} total frames.')
    img = np.sum(im_stack, axis=0)

    # convert to 8-bit
    from imagingplus.utils.images import convert_to_8bit
    img = convert_to_8bit(img, 0, 255)

    if save_path:
        if '.tif' in save_path:
            save_dir = os.path.dirname(save_path) + '/'
            os.makedirs(save_dir, exist_ok=True)
        else:
            raise ValueError(f'please provide full .tif path to save to. provided value was: {tiff_path}')
        print(f"\t\- Writing tiff to: {save_path}")
        tf.imwrite(save_path, img, photometric='minisblack')

    if plot:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(img, cmap='gray')
        fig.suptitle(f'Max Project')
        fig.show()  # just plot for now to make sure that you are doing things correctly so far
        if 'scalebar' in kwargs and kwargs['scalebar']:
            if 'trialobj' in kwargs:
                _add_scalebar(trialobj=kwargs['trialobj'], ax=ax, **kwargs)

    return img


def StdevProject(tiff_path: str = None, frames: tuple = None, save_path: str = None, plot=True,
                 imstack: np.ndarray = None, **kwargs):
    """
    Creates, plots and (optionally) saves an average image of the loaded tiff. If frames are provided, tiff is cropped to those frames.

    Note: saves as 8-bit grayscale image.

    :param tiff_path: path to tiff file to load (if imstack not provided)
    :param imstack: data stack to use as input (if tiff_path not provided)
    :param frames: crop input stack to specified frames before calculating stdev projection.
    :param save_path: path to save output image
    :param plot: if True, show the output image
    :param kwargs:
        :trialobj: ImagingTrial or SingleImage; object associated with input image.
        :scalebar_um: int; size of scalebar to plot on image (in um); must provide trialobj parameter.
    :return: stdev projection image

    """

    if imstack is None and tiff_path is not None:
        # read tiff
        print(f'\t\- Creating std projection, from tiff: {tiff_path}')
        im_stack = ImportTiff(tiff_path=tiff_path, frames=frames)
        # if frames:
        #     im_stack = tf.imread(tiff_path, key=range(frames[0], frames[1]))
        # else:
        #     im_stack = tf.imread(tiff_path)
    elif imstack is not None and tiff_path is None:
        im_stack = imstack
        assert im_stack.ndim == 3, 'can only project 3D image stacks (implicit dimensions are: Frames x Xpixels x Ypixels)'
    else:
        raise ValueError(
            'values provided for both tiff_path and imstack. Unclear where to source image cellsdata from. Provide only one please.')

    if frames:
        print(f'\t collecting image for frames: {frames}')
        im_stack = im_stack[frames[0]: frames[1]]

    print(f'\t\- Std Project {len(im_stack)} total frames.')
    img = np.std(im_stack, axis=0)

    # convert to 8-bit
    from imagingplus.utils.images import convert_to_8bit
    img = convert_to_8bit(img, 0, 255)

    if save_path:
        if '.tif' in save_path:
            save_dir = os.path.dirname(save_path) + '/'
            os.makedirs(save_dir, exist_ok=True)
        else:
            raise ValueError(f'please provide full .tif path to save to. provided value was: {tiff_path}')
        print(f"\t\- Writing tiff to: {save_path}")
        tf.imwrite(save_path, img, photometric='minisblack')

    if plot:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(img, cmap='gray')
        fig.suptitle(f'Std Project')
        fig.show()  # just plot for now to make sure that you are doing things correctly so far
        if 'scalebar' in kwargs and kwargs['scalebar']:
            if 'trialobj' in kwargs:
                _add_scalebar(trialobj=kwargs['trialobj'], ax=ax, **kwargs)

    return img


def InspectTiff(tiff_path: str = None, frames: tuple = None, imstack: np.ndarray = None, **kwargs):
    """
    Plot stdev, max and mean project of input image.

    :param tiff_path: path to tiff file to load (if imstack not provided)
    :param imstack: data stack to use as input (if tiff_path not provided)
    :param frames: crop input stack to specified frames before calculating stdev projection.
    :param kwargs:
        :trialobj: ImagingTrial or SingleImage; object associated with input image.
        :scalebar_um: int; size of scalebar to plot on image (in um); must provide trialobj parameter.

    """
    if imstack is None and tiff_path is not None:
        # read tiff
        print(f'\t\- Creating projections from tiff: {tiff_path}')
        print(f'\t\- Collecting average image for frames: {frames}')
        loaded = ImportTiff(tiff_path=tiff_path, frames=frames)

        # if frames:
        #     loaded = tf.imread(tiff_path, key=range(frames[0], frames[1]))
        # else:
        #     loaded = tf.imread(tiff_path)
    elif imstack and tiff_path is None:
        loaded = imstack
        assert loaded.ndim == 3, 'can only project 3D image stacks (implicit dimensions are: Frames x Xpixels x Ypixels)'
    else:
        raise ValueError(
            'values provided for both tiff_path and imstack. Unclear where to source image cellsdata from. Provide only one please.')

    if frames:
        image = loaded[frames[0]: frames[1]]
    else:
        image = loaded

    _ = StdevProject(imstack=image, **kwargs)
    _ = MaxProject(imstack=image, **kwargs)
    _ = MeanProject(imstack=image, **kwargs)


def FrameAverage(key_frames: Union[int, list], tiff_path: str = None, imstack: np.ndarray = None,
                 peri_frames: int = 100, save_path: str = None, plot=True, **kwargs) -> np.ndarray:
    """
    Creates, plots and/or saves an average image of the specified number of peri-key_frames around the given frame from a multipage imaging TIFF file.

    :param tiff_path: path to tiff file to load (if imstack not provided)
    :param imstack: data stack to use as input (if tiff_path not provided)
    :param key_frames: key frames around which to create peri-frame averages
    :param peri_frames: number of frames around the key frame to calculate average from
    :param plot: if True, show the output image
    :param save_path: path to save output image
    :return: array of images
    """

    print('\nMaking peri-frame avg image...')

    if type(key_frames) == int:
        key_frames = [key_frames]

    if imstack is None and tiff_path:
        # read tiff
        print(f'\t\- Creating avg img for frame: {key_frames}, from tiff: {tiff_path}')
        frames_collect = ((key_frames[0] - peri_frames // 2), (key_frames[-1] + peri_frames // 2))
        im_stack = ImportTiff(tiff_path=tiff_path, frames=frames_collect)
    elif imstack and tiff_path is None:
        im_stack = imstack
    else:
        raise ValueError(
            'values provided for both tiff_path and imstack. Unclear where to source image cellsdata from. Provide only one please.')

    key_frames_adjusted = [frame - (key_frames[0] - peri_frames // 2) for frame in key_frames]

    images = []
    for frame in key_frames_adjusted:
        if len(key_frames_adjusted) > 1:
            if frame not in range(im_stack.shape[0]):
                raise ValueError(f'{frame} not in range of loaded tiff.')

        if frame < peri_frames // 2:
            peri_frames_low = frame
        else:
            peri_frames_low = peri_frames // 2
        peri_frames_high = peri_frames // 2
        im_sub = im_stack[frame - peri_frames_low: frame + peri_frames_high]

        if save_path:
            if '.tif' in save_path:
                save_path = os.path.dirname(save_path) + '/'
            save_path = save_path + f'/{frame}_frame_avg.tif'
            os.makedirs(save_path, exist_ok=True)

            print(f"\t\- Saving averaged tiff for frame: {frame}")
            avg_sub = MeanProject(save_path=save_path, imstack=im_sub)

        else:
            avg_sub = MeanProject(imstack=im_sub, plot=False)

        if plot:
            kwargs['title'] = f'{peri_frames} peri-key_frames avg from frame {frame}' if 'title' not in kwargs else kwargs['title']
            @plotting_decorator(figsize=(6, 6))
            def make_plot(**kwargs):
                fig, ax = kwargs['fig'], kwargs['ax']
                ax.imshow(avg_sub, cmap='gray')
                if 'scalebar_um' in kwargs:
                    fig.tight_layout(pad=0.2)
                    if 'trialobj' in kwargs:
                        _add_scalebar(ax=ax, **kwargs)
                    else:
                        raise ValueError('must provide trialobj parameter to access for making scalebar.')

            # make_plot(**kwargs)
            plotImg(img=avg_sub, **kwargs)
        images.append(avg_sub)

    return np.ndarray(images)


@plotting_decorator(figsize=(6, 6))
def SingleFrame(tiff_path: str = None, frame_num: int = 0, title: str = None, imstack: np.array = None, **kwargs):
    """
    plots an image of a single specified tiff frame after reading using tifffile.

    :param tiff_path: path to .tiff file to loads
    :param frame_num: frame # from 2p imaging tiff to show (default is 0 - i.e. the first frame)
    :param title: (optional) give a string to use as title
    :return: matplotlib imshow plot
    """
    # set any kwargs provided
    ax = kwargs['ax']
    fig = kwargs['fig']

    if imstack is None:
        assert tiff_path is not None, 'please provide a tiff path or input to imstack to use for plotting image.'
        stack = ImportTiff(tiff_path=tiff_path, frames=frame_num)
    else:
        stack = imstack
    ax.imshow(stack, cmap='gray')
    ax.set_title(title) if title is not None else ax.set_title(f'frame num: {frame_num}')
    fig.tight_layout(pad=0.2)

    if 'scalebar_um' in kwargs:
        if 'trialobj' in kwargs:
            _add_scalebar(scalebar_um=kwargs['scalebar_um'], **kwargs)
        else:
            raise ValueError('must provide trialobj parameter to access for making scalebar.')

    fig.show()

    return stack, fig, ax


# plots the raw trace for the Flu mean of the FOV (similar to the ZProject in Fiji)
@plotting_decorator(figsize=(10, 3))
def plotMeanFovFluTrace(trialobj: TwoPhotonImaging, **kwargs):
    """
    Make plot of mean fluorescence trace averaged over the whole FOV

    NOTE: this function will use the underlying timing that is found under trialobj.tmdata (if that is defined). If there is a mismatch in the number of frame timestamps collected and
    the number of imaging frames data, the shorter length will be used for plotting (with no warning printed, as it is often the case that there are dropped frames from the microscope).

    :param trialobj: TwoPhotonImaging (or derived) object
    :param kwargs:
        :ax: matplotlib axis object to use for plotting.
    """

    if not hasattr(trialobj, 'meanFovFluTrace'):
        raise AttributeError('Cannot make plot. Missing meanFovFluTrace attribute.')

    else:
        dataplot_frame_options()
        ax = kwargs['ax']
        kwargs.pop('ax')
        for i in ['lw', 'linewidth']:
            try:
                lw = kwargs[i]
            except KeyError:
                lw = 1

        data_to_plot = trialobj.meanFovFluTrace
        if trialobj.tmdata:
            if len(data_to_plot) < len(trialobj.tmdata.frame_times):
                x_range = trialobj.tmdata.frame_times[:len(data_to_plot)] - trialobj.tmdata.frame_times[0]
            elif len(data_to_plot) > len(trialobj.tmdata.frame_times):
                x_range = trialobj.tmdata.frame_times - trialobj.tmdata.frame_times[0]
                data_to_plot = data_to_plot[:len(x_range)]
            else:
                x_range = trialobj.tmdata.frame_times - trialobj.tmdata.frame_times[0]
            time_res = trialobj.tmdata.sampling_rate
        else:
            x_range = np.arange(0, len(data_to_plot))
            time_res = trialobj.imparams.fps

        print(f"\t \- PLOTTING mean raw flu trace ... ")
        ax.plot(x_range, data_to_plot, c='forestgreen', linewidth=lw)

        ax.set_xlabel('frame #s')
        ax.set_ylabel('Flu (a.u.)')

        # set axis options
        dataplot_ax_options(ax=ax, collection_hz=time_res, **kwargs)


@plotting_decorator(figsize=(10, 6))
def plot_photostim_traces_overlap(array, trialobj: AllOpticalTrial, exclude_id: list = None, y_spacing_factor=1,
                                  title='', x_axis='Time (seconds)', **kwargs):
    """
    Plot fluorescence traces of alloptical stimulated cells on one plot.

    :param array: fluorescence traces of alloptical stimulated cells to plot.
    :param trialobj: AllOptical object to select cells from
    :param exclude_id: list of cell id's to exclude from plot
    :param y_spacing_factor: increase to add more space between plotted traces
    :param title: title of plot
    :param x_axis: x-axis of plot in units of Frames or Time (seconds)
    :param kwargs:
        :ax: matplotlib axis object to use for plotting.


    """

    len_ = len(array)
    dataplot_frame_options()
    ax = kwargs['ax']
    kwargs.pop('ax')

    for i in range(len_):
        if i not in exclude_id:
            if 'linewidth' in kwargs.keys():
                linewidth = kwargs['linewidth']
            else:
                linewidth = 1
            ax.plot(array[i] + i * 40 * y_spacing_factor, linewidth=linewidth)
    for j in trialobj.stim_start_frames:
        if j <= array.shape[1]:
            ax.axvline(x=j, c='gray', alpha=0.3)

    ax.set_xlim([0, trialobj.n_frames - 3000])

    ax.margins(0)
    # change x axis ticks to seconds
    if 'Time' in x_axis or 'time' in x_axis:
        # change x axis ticks to every 30 seconds
        labels = list(range(0, int(trialobj.n_frames // trialobj.imparams.fps), 30))
        ax.set_xticks(ticks=[(label * trialobj.imparams.fps) for label in labels])
        ax.set_xticklabels(labels)
        ax.set_xlabel('Time (secs)')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel(x_axis)

    if 'y_lim' in kwargs.keys():
        ax.set_ylim(kwargs['y_lim'])
    else:
        y_max = np.mean(array[-1] + len_ * 40 * y_spacing_factor) + 3 * np.mean(array[-1])
        ax.set_ylim(0, y_max)

    ax.set_title((title + ' - %s' % len_ + ' cells'), horizontalalignment='center', verticalalignment='top', pad=20,
                 fontsize=10, wrap=True)

    dataplot_ax_options(ax=ax, **kwargs)



### photostim analysis - PLOT avg over photostim. trials traces for the provided traces
# todo decide which photostim avg plot func to keep
@plotting_decorator(figsize=(5, 5.5))
def plot_periphotostim_avg2(dataset, fps=None, legend_labels=None, colors=None, avg_with_std=False,
                            title='high quality plot', pre_stim_sec=None, ylim=None, fig=None, ax=None, **kwargs):
    """
TODO fill documentation and add parameters
    :param dataset:
    :param fps:
    :param legend_labels:
    :param colors:
    :param avg_with_std:
    :param title:
    :param pre_stim_sec:
    :param ylim:
    :param fig:
    :param ax:
    :param kwargs:
    """
    # if 'fig' in kwargs.keys():
    #     fig = kwargs['fig']
    #     ax = kwargs['ax']
    # else:
    #     if 'figsize' in kwargs.keys():
    #         fig, ax = plt.subplots(figsize=kwargs['figsize'])
    #     else:
    #         fig, ax = plt.subplots(figsize=[5, 4])

    meantraces = []
    stdtraces = []
    if type(dataset) == list and len(dataset) > 1:
        assert len(legend_labels) == len(dataset), print('please provide same number of legend labels as dataset')
        if colors is None:
            colors = ['black', make_random_color_array(len(dataset) - 1)]
        assert len(colors) == len(legend_labels)
        avg_only = True
        print('-- plotting average +/- std fill for each dataset')
        for i in range(len(dataset)):
            meanst = np.mean(dataset[i], axis=0)
            std = np.std(dataset[i], axis=0, ddof=1)
            meantraces.append(meanst)
            stdtraces.append(std)
            if not meanst.shape == meantraces[i - 1].shape:
                print(
                    f"|--- length mismatch in mean traces of datasets... {title}, shape0 {meanst.shape} and shape1 {meantraces[i - 1].shape}")
            if not std.shape == stdtraces[i - 1].shape:
                print(
                    f"|--- length mismatch in std traces of datasets...{title}, shape0 {std.shape} and shape1 {stdtraces[i - 1].shape}")

    elif type(dataset) is not list or len(dataset) == 1:
        dataset = list(dataset)
        meanst = np.mean(dataset[0], axis=0)
        std = np.std(dataset[0], axis=0, ddof=1)
        meantraces.append(meanst)
        stdtraces.append(std)
        colors = ['black']
    else:
        AttributeError(
            'please provide the cellsdata to plot in a ls format, each different cellsdata group as a ls item...')

    if 'xlabel' not in kwargs or kwargs['xlabel'] is None or 'key_frames' not in kwargs['xlabel'] or 'Frames' not in \
            kwargs[
                'xlabel']:
        ## change xaxis to time (secs)
        if fps is not None:
            if pre_stim_sec is not None:
                x_range = np.linspace(0, len(meantraces[0]) / fps, len(
                    meantraces[
                        0])) - pre_stim_sec  # x scale, but in time domain (transformed from key_frames based on the provided fps)
                if 'xlabel' in kwargs.keys():
                    ax.set_xlabel(kwargs['xlabel'])
                else:
                    ax.set_xlabel('Time post stim (secs)')
            else:
                AttributeError('need to provide a pre_stim_sec value to the function call!')
        else:
            AttributeError('need to provide fps value to convert xaxis in units of time (secs)')
    elif 'key_frames' in kwargs['xlabel'] or 'Frames' in kwargs['xlabel']:
        x_range = range(len(meanst[0]))
        ax.set_xlabel('Frames')

    for i in range(len(meantraces)):
        if avg_with_std:
            if len(meantraces[i]) < len(x_range):  ## TEMP FIX
                mismatch = (len(x_range) - len(meantraces[i]))
                meantraces[i] = np.append(meantraces[i], [0] * mismatch)
                stdtraces[i] = np.append(stdtraces[i], [0] * mismatch)
                print(
                    f'|------ adding {mismatch} zeros to mean and std-fill traces to make the arrays the same length, new length of plot array: {meantraces[i].shape} ')

            ax.plot(x_range, meantraces[i], color=colors[i], lw=2)
            ax.fill_between(x_range, meantraces[i] - stdtraces[i], meantraces[i] + stdtraces[i], alpha=0.15,
                            color=colors[i])
        else:
            ax.plot(x_range, meantraces[i], color=colors[i], lw=2)
            for trace in dataset[i]:
                ax.plot(x_range, trace, color=colors[i], alpha=0.3, lw=2)

    if legend_labels:
        if 'fontsize' in kwargs.keys():
            fontsize = kwargs['fontsize']
        else:
            fontsize = 'medium'
        ax.legend(legend_labels, fontsize=fontsize)
    if ylim:
        ax.set_ylim([ylim[0], ylim[1]])
    if 'ylabel' in kwargs.keys():
        ax.set_ylabel(kwargs['ylabel'])
    else:
        pass

    if 'savepath' in kwargs.keys():
        plt.savefig(kwargs['savepath'])

    # # finalize plot, set title, and show or return axes
    # ax.set_title(title)
    # if 'show' in kwargs.keys():
    #     fig.show() if kwargs['show'] else None
    # else:
    #     fig.show()
    #
    # if 'fig' in kwargs.keys():
    #     return fig, ax


### photostim analysis - PLOT avg over all photstim. trials traces from PHOTOSTIM TARGETTED cells
@plotting_decorator(figsize=(5, 5.5))
def plot_periphotostim_avg(arr: np.ndarray, trialobj: AllOpticalTrial, pre_stim_sec=1.0, post_stim_sec=3.0, title='',
                           avg_only: bool = False, x_label=None, y_label=None, pad=20, **kwargs):
    """
    plot trace across all stims
    :param arr: Flu traces to plot (will be plotted as individual traces unless avg_only is True) dimensions should be cells x stims x key_frames
    :param trialobj: instance object of AllOpticalTrial
    :param pre_stim_sec: seconds of array to plot for pre-stim period
    :param post_stim_sec: seconds of array to plot for post-stim period
    :param title: title to use for plot
    :param avg_only: if True, only plot the mean trace calculated from the traces provided in arr
    :param x_label: x axis label
    :param y_label: y axis label
    :param kwargs:
        options include:
            'stim_duration': photostimulation duration in secs
            'y_lims': tuple, y min and max of the plot
            'edgecolor': str, edgecolor of the individual traces behind the mean trace
            'savepath': str, path to save plot to

    :param fig: a matplotlib.Figure instance, if provided use this fig for plotting
    :param ax: a matplotlib.Axes.axes instance, if provided use this ax for plotting
    :param show: if False, do not display plot (used when the necessity is to return the fig and ax objects to futher manipulation)
    :returns fig, ax: returns fig and ax object, if show is False

    """

    fps = trialobj.imparams.fps  # key_frames per second rate of the imaging cellsdata collection for the cellsdata to be plotted
    exp_prestim = trialobj.pre_stim_frames  # key_frames of pre-stim cellsdata collected for each trace for this trialobj (should be same as what's under trialobj.pre_stim_sec)
    if 'stim_duration' in kwargs.keys():
        stim_duration = kwargs['stim_duration']
    else:
        stim_duration = trialobj.twopstim.stim_dur / 1000  # seconds of stimulation duration

    dataplot_frame_options()
    ax = kwargs['ax']
    kwargs.pop('ax')

    x = list(range(arr.shape[1]))
    # x range in time (secs)
    x_time = np.linspace(0, arr.shape[1] / fps, arr.shape[
        1]) - pre_stim_sec  # x scale, but in time domain (transformed from key_frames based on the provided fps)

    len_ = len(arr)
    flu_avg = np.mean(arr, axis=0)

    # ax.margins(x=0.07)

    if 'alpha' in kwargs.keys():
        alpha = kwargs['alpha']
    else:
        alpha = 0.2

    if x_label is None or not 'Frames' in x_label or 'key_frames' in x_label:
        x = x_time  # set the x plotting range
        if x_label is not None:
            x_label = x_label + 'post-stimulation relative'
        else:
            x_label = 'Time (secs post-stimulation)'

        if avg_only is True:
            # ax.axvspan(exp_prestim/fps, (exp_prestim + stim_duration + 1) / fps, alpha=alpha, color='plum', zorder = 3)
            ax.axvspan(0 - 1 / fps, 0 + stim_duration + 1 / fps, alpha=alpha, color='plum',
                       zorder=3)  # note that we are setting 0 as the stimulation time
        else:
            ax.axvspan(0 - 1 / fps, 0 - 1 / fps + stim_duration, alpha=alpha, color='plum',
                       zorder=3)  # note that we are setting 0 as the stimulation time
    else:
        ax.axvspan(exp_prestim, exp_prestim + int(stim_duration * fps), alpha=alpha, color='tomato')

    if not avg_only:
        for cell_trace in arr:
            if 'color' in kwargs.keys():
                ax.plot(x, cell_trace, linewidth=1, alpha=0.6, c=kwargs['color'], zorder=1)
            else:
                if arr.shape[0] > 50:
                    alpha = 0.1
                else:
                    alpha = 0.5
                ax.plot(x, cell_trace, linewidth=1, alpha=alpha, zorder=1)

    ax.plot(x, flu_avg, color='black', linewidth=2.3, zorder=2)  # plot average trace

    if 'y_lims' in kwargs.keys():
        ax.set_ylim(kwargs['y_lims'])
    if pre_stim_sec and post_stim_sec:
        if x_label is None or not 'Frames' in x_label or 'key_frames' in x_label:
            ax.set_xlim(-pre_stim_sec,
                        stim_duration + post_stim_sec)  # remember that x axis is set to be relative to the stim time (i.e. stim is at t = 0)
        else:
            ax.set_xlim(exp_prestim - int(pre_stim_sec * fps),
                        exp_prestim + int(stim_duration * fps) + int(post_stim_sec * fps) + 1)

    # set axis labels
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if 'savepath' in kwargs.keys():
        plt.savefig(kwargs['savepath'])

    if title is not None:
        ax.set_title((title + ' - %s' % len_ + ' traces'), horizontalalignment='center', verticalalignment='top',
                     pad=pad, fontsize=10, wrap=True)

    dataplot_ax_options(ax=ax, **kwargs)


# alloptical trial
### plot the location of all SLM targets, along with option for plotting the mean img of the current trial

@mpl.rc_context(image_frame_ops)
@plotting_decorator(figsize=(5, 5))
def plot_SLMtargets_Locs(trialobj: AllOpticalTrial, targets_coords: Union[list, str] = 'all',
                         background: np.ndarray = None, **kwargs):
    """
    plot SLM target coordinate locations

    :param trialobj:
    :param targets_coords: ls containing (x,y) coordinates of targets to plot
    :param background:
    :param kwargs:
    :param fig: a matplotlib.Figure instance, if provided use this fig for plotting
    :param ax: a matplotlib.Axes.axes instance, if provided use this ax for plotting
    :param show: if False, do not display plot (used when the necessity is to return the fig and ax objects to futher manipulation)
    :returns fig, ax: returns fig and ax object, if show is False

    """

    if not type(trialobj) == AllOpticalTrial:
        raise ObjectClassError(function='plot_SLMtargets_Locs', valid_class=[AllOpticalTrial],
                               invalid_class=type(trialobj))
    ax = kwargs['ax']
    fig = kwargs['fig']
    kwargs.pop('ax')

    if background is None:
        background = np.zeros((trialobj.imparams.frame_x, trialobj.imparams.frame_y), dtype='uint16')
        ax.imshow(background, cmap='gray')
    else:
        ax.imshow(background, cmap='gray')

    colors = make_random_color_array(len(trialobj.twopstim.target_coords))
    if targets_coords is 'all':
        if len(trialobj.twopstim.target_coords) > 1:
            for i in range(len(trialobj.twopstim.target_coords)):
                ax.scatter(x=trialobj.twopstim.target_coords[i][:, 0], y=trialobj.twopstim.target_coords[i][:, 1],
                           edgecolors=colors[i], facecolors='none', linewidths=2.0, label=f'SLM Group {i}')
                # for (x, y) in trialobj.twopstim.target_coords[i]:
                #     ax.scatter(x=x, y=y, edgecolors=colors[i], facecolors='none', linewidths=2.0)
        else:
            if 'edgecolors' in kwargs.keys():
                edgecolors = kwargs['edgecolors']
            else:
                edgecolors = 'yellowgreen'
            for (x, y) in trialobj.twopstim.target_coords_all:
                ax.scatter(x=x, y=y, edgecolors=edgecolors, facecolors='none', linewidths=2.0)
    elif targets_coords:
        if 'edgecolors' in kwargs.keys():
            edgecolors = kwargs['edgecolors']
        else:
            edgecolors = 'yellowgreen'
        plot_coordinates(coords=targets_coords, frame_x=trialobj.imparams.frame_x, frame_y=trialobj.imparams.frame_y,
                         edgecolors=edgecolors,
                         background=background, fig=fig, ax=ax)

    ax.legend(loc='upper right', labelcolor='white', frameon=False)
    ax.margins(0)

    ax = _add_scalebar(trialobj=trialobj, ax=ax)

    # fig.tight_layout()

    if 'title' in kwargs.keys():
        if kwargs['title'] is not None:
            ax.set_title(kwargs['title'])
        else:
            pass
    else:
        ax.set_title(f'SLM targets location - {trialobj.t_series_name}')


# alloptical trial plotting


def plot_s2pMasks():
    """ Creates some images of SLM targets to suite2p cells.

    """
    pass


# todo need to find new place for this:
def export_klicker_to_csv(klicker: mpl_point_clicker.clicker, csv_path):
    """Save point cellsdata from mplpointclicker klicker to .csv
    :param klicker: mplpointclicker klicker object
    :param csv_path: path to save csv file of klicker points
    """
    _columns = [*klicker.get_positions()]
    data = {}
    for i in [*klicker.get_positions()]:
        data[i] = klicker.get_positions()[i][:,
                  0]  # take just the x coord of the clicked point from klicker (i.e. time value)
    df = pd.DataFrame(data)
    save_to_csv(df, savepath=csv_path)
