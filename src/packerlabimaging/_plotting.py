# library of convenience plotting funcs that are used for making various plots for all optical photostimulation/imaging experiments

# imports
import os

import matplotlib
import numpy as np
import functools
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import tifffile as tf
from .AllOpticalMain import AllOpticalTrial

from .TwoPhotonImagingMain import TwoPhotonImagingTrial
from ._utils import normalize_dff

# wrapper for piping plots in and out of figures
def plot_piping_decorator(figsize=(5,5)):
    def plot_piping_decorator_(plotting_func):
        """
        Wrapper to help simplify creating plots from matplotlib.pyplot

        :param plotting_func: plotting convenience function to be wrapped
        :return: fig and ax objects if requested, and/or shows the figure as specified

        Examples:
        ---------
        1) in this example the fig and ax will be taken directly from the kwargs inside the inner wrapper
        >>> @plot_piping_decorator
        >>> def example_decorated_plot(title='', **kwargs):
        >>>     fig, ax = kwargs['fig'], kwargs['ax']  # need to add atleast this line to each plotting function's definition
        >>>     print(f'|-kwargs inside example_decorated_plot definition: {kwargs}')
        >>>     ax.plot(np.random.rand(10))
        >>>     ax.set_title(title)

        >>> example_decorated_plot()


        """

        @functools.wraps(plotting_func)
        def inner(**kwargs):
            # print(f'perform fig, ax creation')
            # print(f'|-original kwargs {kwargs}')
            return_fig_obj = False

            # set number of rows, cols and figsize
            if 'nrows' in kwargs.keys():
                nrows = kwargs['nrows']
            else:
                nrows = 1

            if 'ncols' in kwargs.keys():
                ncols = kwargs['ncols']
            else:
                ncols = 1

            if 'figsize' in kwargs.keys():
                figsize_ = kwargs['figsize']
            else:
                figsize_ = figsize

            # create or retrieve the fig, ax objects --> end up in kwargs to use into the plotting func call below
            if 'fig' in kwargs.keys() and 'ax' in kwargs.keys():
                if kwargs['fig'] is None or kwargs['ax'] is None:
                    # print('\-creating fig, ax [1]')
                    kwargs['fig'], kwargs['ax'] = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize_)
            else:
                # print('\-creating fig, ax [2]')
                kwargs['fig'], kwargs['ax'] = plt.subplots(figsize=figsize_)

            # print(f"\nnew kwargs {kwargs}")

            print(f'\- executing plotting_func')
            res = plotting_func(**kwargs)  # these kwargs are the original kwargs defined at the respective plotting_func call + any additional kwargs defined in inner()

            # print(f'\nreturn fig, ax or show figure as called for')
            kwargs['fig'].suptitle(kwargs['suptitle'], wrap=True) if 'suptitle' in kwargs.keys() else None

            kwargs['fig'].tight_layout(pad=1.8)

            if 'show' in kwargs.keys():
                if kwargs['show'] is True:
                    # print(f'\- showing fig of size {figsize_}...[3]')
                    kwargs['fig'].show()
                    # print(f"*res right now: {res}")
                    return res
                else:
                    # print(f"\- not showing, but returning fig_obj of size {figsize_}[4]")
                    if res is not None:
                        return (kwargs['fig'], kwargs['ax'], res)
                    else:
                        return (kwargs['fig'], kwargs['ax'])
            else:
                # print(f'\- showing fig of size {figsize_}...[5]')
                kwargs['fig'].show()
                return res

            # # print(f"|-value of return_fig_obj is {return_fig_obj} [5]")
            # print(f"\- returning fig_obj [4]") if return_fig_obj else None
            # return (kwargs['fig'], kwargs['ax']) if return_fig_obj else None

        return inner
    return plot_piping_decorator_

# %% UTILITY FUNCS
def save_figure(fig, save_path_full: str = None):
    print(f'\nsaving figure to: {save_path_full}')
    os.makedirs(save_path_full)
    fig.savefig(save_path_full)

# custom colorbar for heatmaps
from matplotlib.colors import LinearSegmentedColormap
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return LinearSegmentedColormap('CustomMap', cdict)


# generate an array of random colors
def _get_random_color(pastel_factor=0.5):
    return [(x + pastel_factor) / (1.0 + pastel_factor) for x in [random.uniform(0, 1.0) for i in [1, 2, 3]]]


def _color_distance(c1, c2):
    return sum([abs(x[0] - x[1]) for x in zip(c1, c2)])


def _generate_new_color(existing_colors, pastel_factor=0.5):
    max_distance = None
    best_color = None
    for i in range(0, 100):
        color = _get_random_color(pastel_factor=pastel_factor)
        if not existing_colors:
            return color
        best_distance = min([_color_distance(color, c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    return best_color


def make_random_color_array(n_colors):
    """
    Generates a list of random colors for an input number of colors required.

    :param n_colors: # of colors required
    :return: list of colors in RGB
    """
    colors = []
    for i in range(0, n_colors):
        colors.append(_generate_new_color(colors, pastel_factor=0.2))
    return colors


def add_scalebar(trialobj: TwoPhotonImagingTrial, ax: plt.Axes, scale_bar_um: float = 100, **kwargs):
    """add scalebar to the image being plotted on the a single matplotlib.axes.Axes object using the TwoPhotonImaging object information.
    Option to specify scale bar um length to add to plot.

    """
    numpx = scale_bar_um/trialobj.imparams.pix_sz_x

    lw = 5 if 'lw' not in [*kwargs] else kwargs['lw']
    color = 'white' if 'color' not in [*kwargs] else kwargs['color']
    right_offset = 50 if 'right_offset' not in [*kwargs] else kwargs['right_offset']
    bottom_offset = 50 if 'bottom_offset' not in [*kwargs] else kwargs['bottom_offset']

    ax.plot(np.linspace(trialobj.imparams.frame_x - right_offset - numpx, trialobj.imparams.frame_x - right_offset, 40),
            [trialobj.imparams.frame_y - bottom_offset]*40,
            color=color, lw=lw)
    return ax

# Figure Style settings for notebook.
def image_frame_options():
    mpl.rcParams.update({
        'axes.spines.left': False,
        'axes.spines.bottom': False,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.frameon': False,
        'figure.subplot.wspace': .01,
        'figure.subplot.hspace': .01,
        'ytick.major.left': False,
        'xtick.major.bottom': False
    })

def heatmap_options():
    jet = mpl.cm.get_cmap('jet')
    jet.set_bad(color='k')

# %% GENERAL PLOTTING FUNCS
### plot the location of provided coordinates
@plot_piping_decorator(figsize=(5,5))
def plot_coordinates(coords: list,  frame_x: int, frame_y: int, background: np.ndarray = None, fig=None, ax=None, **kwargs):
    """
    plot coordinate locations

    :param targets_coords: ls containing (x,y) coordinates of targets to plot
    :param background: np.array on which to plot coordinates, default is black background (optional)
    :param kwargs:
    """

    if background is None:
        background = np.zeros((frame_x, frame_y), dtype='uint16')
        ax.imshow(background, cmap='gray')
    else:
        ax.imshow(background, cmap='gray')

    if 'edgecolors' in kwargs.keys():
        edgecolors = kwargs['edgecolors']
    else:
        edgecolors = 'yellowgreen'
    for (x, y) in coords:
        ax.scatter(x=x, y=y, edgecolors=edgecolors, facecolors='none', linewidths=2.0)

    ax.margins(0)
    fig.tight_layout()

    if 'title' in kwargs.keys():
        if kwargs['title'] is not None:
            ax.set_title(kwargs['title'])
        else:
            pass


# %% PROJECT SPECIFIC PLOTTING FUNCS


# suite2p data
# simple plot of the location of the given cell(s) against a black FOV
@plot_piping_decorator(figsize=(5,5))
def plot_cells_loc(expobj: TwoPhotonImagingTrial, cells: list, title=None, background: np.array = None, scatter_only: bool = False,
                   show_s2p_targets: bool = True, color_float_list: list = None, cmap: str = 'Reds', invert_y = True, fig=None, ax=None, **kwargs):
    """
    plots an image of the FOV to show the locations of cells given in cells ls.
    :param background: either 2dim numpy array to use as the backsplash or None (where black backsplash will be created)
    :param expobj: alloptical or 2p imaging object
    :param edgecolor: str to specify edgecolor of the scatter plot for cells
    :param cells: ls of cells to plot
    :param title: str title for plot
    :param color_float_list: if given, it will be used to color the cells according a colormap
    :param cmap: cmap to be used in conjuction with the color_float_array argument
    :param show_s2p_targets: if True, then will prioritize coloring of cell points based on whether they were photostim targets
    :param kwargs: optional arguments
            invert_y: if True, invert the reverse the direction of the y axis
            show: if True, show the plot
            fig: a fig plt.subplots() instance, if provided use this fig for making figure
            ax: a ax plt.subplots() instance, if provided use this ax for plotting
    """

    # # if there is a fig and ax provided in the function call then use those, otherwise start anew
    # if 'fig' in kwargs.keys():
    #     fig = kwargs['fig']
    #     ax = kwargs['ax']
    # else:
    #     fig, ax = plt.subplots()

    x_list = []
    y_list = []
    for cell in cells:
        y, x = expobj.stat[expobj.cell_id.index(cell)]['med']
        x_list.append(x)
        y_list.append(y)

        if show_s2p_targets:
            if hasattr(expobj, 's2p_cell_targets'):
                if cell in expobj.s2p_cell_targets:
                    color_ = '#F02A71'
                else:
                    color_ = 'none'
            else:
                color_ = 'none'
            ax.scatter(x=x, y=y, edgecolors=None, facecolors=color_, linewidths=0.8)
        elif color_float_list:
            # ax.scatter(x=x, y=y, edgecolors='none', c=color_float_list[cells.index(cell)], linewidths=0.8,
            #            cmap=cmap)
            pass
        else:
            if 'edgecolors' in kwargs.keys():
                edgecolors = kwargs['edgecolors']
            else:
                edgecolors = 'yellowgreen'
            ax.scatter(x=x, y=y, edgecolors=edgecolors, facecolors='none', linewidths=0.8)

    if color_float_list:
        ac = ax.scatter(x=x_list, y=y_list, edgecolors='none', c=color_float_list, linewidths=0.8,
                   cmap=cmap, zorder=1)

        plt.colorbar(ac, ax=ax)

    if not scatter_only:
        if background is None:
            black = np.zeros((expobj.frame_x, expobj.frame_y), dtype='uint16')
            ax.imshow(black, cmap='Greys_r', zorder=0)
            ax.set_xlim(0, expobj.frame_x)
            ax.set_ylim(0, expobj.frame_y)
        else:
            ax.imshow(background, cmap='Greys_r', zorder=0)

    if title is not None:
        plt.suptitle(title, wrap=True)

    if 'text' in kwargs.keys():
        if kwargs['text'] is not None:
            ax.text(0.99, 0.95, kwargs['text'],
                    verticalalignment='top', horizontalalignment='right',
                    transform=ax.transAxes, fontweight='bold',
                    color='white', fontsize=10)

    if 'hide_axis_labels' in kwargs.keys():
        ax.set_xticks(ticks=[])
        ax.set_xticklabels([])
        ax.set_yticks(ticks=[])
        ax.set_yticklabels([])


    if 'invert_y' in kwargs.keys():
        if kwargs['invert_y']:
            ax.invert_yaxis()

    ax = add_scalebar(trialobj=expobj, ax=ax) if 'scalebar' in [*kwargs] and kwargs['scalebar'] is True else None


# plot to show s2p ROIs location, colored as specified
def s2pRoiImage(expobj: TwoPhotonImagingTrial, save_fig: str = None):
    """
    plot to show the classification of each cell as the actual's filling in the cell's ROI pixels.

    :param expobj: expobj associated with data
    :param df: pandas dataframe (cell_id x stim frames)
    :param clim: color limits
    :param plot_target_coords: bool, if True plot the actual X and Y coords of all photostim cell targets
    :param save_fig: where to save the save figure (optional)
    :return:
    """
    fig, ax = plt.subplots(figsize=(5, 5))
    if expobj.frame_x == 512:
        s = 0.003 * (1024/expobj.frame_x * 4)
    else:
        s = 0.003
    ##### targets areas image
    targ_img = np.zeros([expobj.frame_x, expobj.frame_y], dtype='float')
    target_areas_exclude = np.array(expobj.target_areas_exclude)
    targ_img[target_areas_exclude[:, :, 1], target_areas_exclude[:, :, 0]] = 1
    x = np.asarray(list(range(expobj.frame_x)) * expobj.frame_y)
    y = np.asarray([i_y for i_y in range(expobj.frame_y) for i_x in range(expobj.frame_x)])
    img = targ_img.flatten()
    im_array = np.array([x, y], dtype=np.float)
    ax.scatter(im_array[0], im_array[1], c=img, cmap='gray', s=s, zorder=0, alpha=1)

    ##### suite2p ROIs areas image - nontargets
    for n in expobj.s2p_nontargets:
        idx = expobj.cell_id.index(n)
        ypix = expobj.stat[idx]['ypix']
        xpix = expobj.stat[idx]['xpix']
        ax.scatter(xpix, ypix, c='lightsteelblue', s=s, zorder=1, alpha=1)

    ##### suite2p ROIs areas image - exclude cells
    for n in expobj.s2p_cells_exclude:
        idx = expobj.cell_id.index(n)
        ypix = expobj.stat[idx]['ypix']
        xpix = expobj.stat[idx]['xpix']
        ax.scatter(xpix, ypix, c='yellow', s=s, zorder=2, alpha=1)

    ##### suite2p ROIs areas image - targeted cells
    for n in expobj.s2p_cell_targets:
        idx = expobj.cell_id.index(n)
        ypix = expobj.stat[idx]['ypix']
        xpix = expobj.stat[idx]['xpix']
        ax.scatter(xpix, ypix, c='red', s=s, zorder=3, alpha=1)

    ax.set_xlim([0, expobj.frame_x])
    ax.set_ylim([0, expobj.frame_y])

    ax = add_scalebar(expobj=expobj, ax=ax)

    plt.margins(x=0, y=0)
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    # fig.show()

    plt.suptitle(f"{expobj.metainfo['animal prep.']} {expobj.metainfo['trial']} - s2p nontargets (blue), exclude (yellow), targets (red); target_areas (white)",
                 y=0.97, fontsize=7)
    plt.show()
    save_figure(fig, save_path_suffix=f"{save_fig}") if save_fig else None

# (full) plot individual cell's flu or dFF trace, with photostim. timings for that cell
def plot_flu_trace(expobj: TwoPhotonImagingTrial, cell, x_lims=None, slm_group=None, to_plot='raw', figsize=(20, 3), linewidth=0.10, show=True):
    idx = expobj.cell_id.index(cell)
    raw = expobj.raw[idx]
    raw_ = np.delete(raw, expobj.photostim_frames)  # this is very problematic for the dFF plotting with stim frames if you're deleting ALL of the photostim frames!?!!!
    raw_dff = normalize_dff(raw_)
    std_dff = np.std(raw_dff, axis=0)
    std = np.std(raw_, axis=0)

    # find locations along time when the trace rises above 2.5std of the mean
    x = []
    # y = []
    for j in np.arange(len(raw_dff), step=4):
        avg = np.mean(raw_dff[j:j + 4])
        if avg > np.mean(raw_dff) + 2 * std_dff:
            x.append(j)
            # y.append(0)

    if to_plot == 'raw':
        to_plot_ = raw
        to_thresh = std
    elif to_plot == 'dff':
        to_plot_ = raw_dff
        to_thresh = std_dff
    else:
        AttributeError('specify to_plot as either "raw" or "dff"')

    # make the plot either as just the raw trace or as a dFF trace with the std threshold line drawn as well.
    plt.figure(figsize=figsize)
    plt.plot(to_plot_, linewidth=linewidth)
    if to_plot == 'raw':
        plt.suptitle(('raw flu for cell #%s' % expobj.cell_id[idx]), horizontalalignment='center',
                     verticalalignment='top',
                     fontsize=15, y=1.00)
    elif to_plot == 'dff':
        plt.scatter(x, y=[0] * len(x), c='r', linewidth=0.1)
        plt.axhline(y=np.mean(to_plot_) + 2.5 * to_thresh, c='green')  # plot threshold line
        plt.suptitle(('dff flu for cell #%s' % expobj.cell_id[idx]), horizontalalignment='center',
                     verticalalignment='top', fontsize=15, y=0.95)

    if slm_group is not None:
        for i in expobj.stim_start_frames[slm_group::expobj.n_groups]:  # select SLM group specific stim trigger frames (may not exist by each individual SLM group though...)
            plt.axvline(x=i - 1, c='gray', alpha=0.1)
    else:
        for i in expobj.stim_start_frames:  # select all stim trigger frames from the trial
            plt.axvline(x=i - 1, c='gray', alpha=0.1)

    if expobj.seizure_frames:
        plt.scatter(expobj.seizure_frames, y=[-20] * len(expobj.seizure_frames), c='g', linewidth=0.10)

    if x_lims:
        plt.xlim(x_lims)

    # plt.ylim(0, 300)
    plt.show() if show else None

# plots the raw trace for the Flu mean of the FOV (similar to the ZProject in Fiji)
@plot_piping_decorator(figsize=(10,3))
def plotMeanRawFluTrace(expobj: TwoPhotonImagingTrial, stim_span_color='white', stim_lines: bool = True, title='raw Flu trace', x_axis='time', shrink_text=1,
                        fig=None, ax=None, **kwargs):
    """make plot of mean Ca trace averaged over the whole FOV"""

    print(f"\t \- PLOTTING mean raw flu trace ... ")

    # # if there is a fig and ax provided in the function call then use those, otherwise start anew
    # if 'fig' in kwargs.keys():
    #     fig = kwargs['fig']
    #     ax = kwargs['ax']
    # else:
    #     if 'figsize' in kwargs.keys():
    #         fig, ax = plt.subplots(figsize=kwargs['figsize'])
    #     else:
    #         if 'xlims' in kwargs.keys():
    #             fig, ax = plt.subplots(figsize=[10 * (kwargs['xlims'][1] - kwargs['xlims'][0]) / 2000, 3])
    #             fig.tight_layout(pad=0)
    #         else:
    #             fig, ax = plt.subplots(figsize=[10 * len(expobj.meanRawFluTrace) / 2000, 3])
    #             fig.tight_layout(pad=0)

    # change linewidth
    if 'linewidth' in kwargs:
        lw = kwargs['linewidth']
    else:
        lw = 2

    ax.plot(expobj.meanRawFluTrace, c='forestgreen', zorder=1, linewidth=lw)
    ax.margins(0)
    if stim_span_color is not None:
        if hasattr(expobj, 'shutter_frames'):
            for start, end in zip(expobj.shutter_start_frames[0], expobj.shutter_end_frames[0]):
                ax.axvspan(start-4.5, end, color=stim_span_color, zorder=2)
        else:
            for stim in expobj.stim_start_frames:
                ax.axvspan(stim-2, 1 + stim + expobj.stim_duration_frames, color=stim_span_color, zorder=2)
    if stim_lines:
        if stim_span_color is not None:
            for line in expobj.stim_start_frames:
                ax.axvline(x=line, color='black', linestyle='--', linewidth=0.6, zorder=2)
        else:
            for line in expobj.stim_start_frames:
                ax.axvline(x=line, color='black', linestyle='--', linewidth=0.6, zorder=0)
    if 'time' in x_axis or 'Time' in x_axis:
        # change x axis ticks to every 30 seconds
        labels = list(range(0, int(len(expobj.meanRawFluTrace) // expobj.fps), 30))
        ax.set_xticks(ticks=[(label * expobj.fps) for label in labels])
        # for item in labels:
        #     labels[labels.index(item)] = int(round(item / expobj.fps))
        # ticks_loc = ax.get_xticks().tolist()
        # ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax.set_xticklabels(labels)
        ax.set_xlabel('Time (secs)')
    elif x_axis == 'frames':
        ax.set_xlabel('frame #s')
    ax.set_ylabel('Flu (a.u.)')
    if 'xlims' in kwargs.keys() and kwargs['xlims'] is not None:
        ax.set_xlim(kwargs['xlims'])

    if 'ylims' in kwargs.keys() and kwargs['ylims'] is not None:
        ax.set_ylim(kwargs['ylims'])

    ax.set_title('%s %s %s %s' % (title, expobj.metainfo['exptype'], expobj.metainfo['animal prep.'], expobj.metainfo['trial']))

    return None
    # # add title
    # if not 'fig' in kwargs.keys():
    #     ax.set_title(
    #         '%s %s %s %s' % (title, expobj.metainfo['exptype'], expobj.metainfo['animal prep.'], expobj.metainfo['trial']))
    #
    # if 'show' in kwargs.keys():
    #     plt.show() if kwargs['show'] else None
    # else:
    #     plt.show()
    #
    # if 'fig' in kwargs.keys():
    #     # adding text because adding title doesn't seem to want to work when piping subplots
    #     ax.text(0.98, 0.97, f"Avg. FOV Flu Trace - {expobj.metainfo['exptype']} {expobj.metainfo['animal prep.']} {expobj.metainfo['trial']}",
    #             verticalalignment='top', horizontalalignment='right',
    #             transform=ax.transAxes, fontweight='bold',
    #             color='black', fontsize=10 * shrink_text)
    #     return fig, ax

def plot_s2p_raw(expobj, cell_id):
    plt.figure(figsize=(50, 3))
    plt.plot(expobj.baseline_raw[expobj.cell_id.index(cell_id)], linewidth=0.5, c='black')
    plt.xlim(0, len(expobj.baseline_raw[0]))
    plt.show()

# suite2p data


# LFP
@plot_piping_decorator(figsize=(10,3))
def plotLfpSignal(expobj, stim_span_color='powderblue', downsample: bool = True, stim_lines: bool = True, sz_markings: bool = False,
                  title='LFP trace', x_axis='time', hide_xlabel=False, fig=None, ax=None, **kwargs):
    """make plot of LFP with also showing stim locations
    NOTE: ONLY PLOTTING LFP SIGNAL CROPPED TO 2P IMAGING FRAME START AND END TIMES - SO CUTTING OUT THE LFP SIGNAL BEFORE AND AFTER"""

    print(f"\t \- PLOTTING LFP Signal trace ... ")

    # # if there is a fig and ax provided in the function call then use those, otherwise start anew
    # if 'fig' in kwargs.keys():
    #     fig = kwargs['fig']
    #     ax = kwargs['ax']
    # else:
    #     if 'figsize' in kwargs.keys():
    #         fig, ax = plt.subplots(figsize=kwargs['figsize'])
    #     else:
    #         fig, ax = plt.subplots(figsize=[60 * (expobj.stim_start_times[-1] + 1e5 - (expobj.stim_start_times[0] - 1e5)) / 1e7, 3])

    if 'alpha' in kwargs:
        alpha = kwargs['alpha']
    else:
        alpha = 1

    # plot LFP signal
    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = 'steelblue'

    # option for downsampling of data plot trace
    x = range(len(expobj.lfp_signal[expobj.frame_start_time_actual: expobj.frame_end_time_actual]))
    signal = expobj.lfp_signal[expobj.frame_start_time_actual: expobj.frame_end_time_actual]
    if downsample:
        labels = list(range(0, int(len(signal) / expobj.paq_rate * 1), 15))[::2]  # set x ticks at every 30 secs
        down = 1000
        signal = signal[::down]
        x = x[::down]
        assert len(x) == len(signal), print('something went wrong with the downsampling')

    # change linewidth
    if 'linewidth' in kwargs:
        lw = kwargs['linewidth']
    else:
        lw = 0.4

    ax.plot(x, signal, c=color, zorder=1, linewidth=lw)  ## NOTE: ONLY PLOTTING LFP SIGNAL RELATED TO
    ax.margins(0.02)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.spines['left'].set_visible(False)


    # plot stims
    if stim_span_color != '':
        for stim in expobj.stim_start_times:
            stim = (stim - expobj.frame_start_time_actual)
            ax.axvspan(stim - 8, 1 + stim + expobj.stim_duration_frames / expobj.fps * expobj.paq_rate, color=stim_span_color, zorder=0, alpha=alpha)
    else:
        if stim_lines:
            for line in expobj.stim_start_times:
                line = (line - expobj.frame_start_time_actual)
                ax.axvline(x=line+2, color='black', linestyle='--', linewidth=0.6, zorder=0)

    # plot seizure onset and offset markings
    if sz_markings:
        if hasattr(expobj, 'seizure_lfp_onsets'):
            for sz_onset in expobj.seizure_lfp_onsets:
                ax.axvline(x=expobj.frame_clock_actual[sz_onset] - expobj.frame_start_time_actual, color='black', linestyle='--', linewidth=1.0, zorder=0)
            for sz_offset in expobj.seizure_lfp_offsets:
                ax.axvline(x=expobj.frame_clock_actual[sz_offset] - expobj.frame_start_time_actual, color='gray', linestyle='--', linewidth=1.0, zorder=0)

    # change x axis ticks to seconds
    labels_ = kwargs['labels'] if 'labels' in [*kwargs] else labels
    labels_ = [int(i) for i in labels_]
    if 'time' in x_axis or 'Time' in x_axis:
        ax.set_xticks(ticks=[(label * expobj.paq_rate) for label in labels_])
        ax.set_xticklabels(labels_)
        ax.tick_params(axis='both', which='both', length=3)
        if not hide_xlabel:
            ax.set_xlabel('Time (secs)')
    # elif 'frame' or "frames" or "Frames" or "Frame" in x_axis:
    #     x_ticks = range(0, expobj.n_frames, 2000)
    #     x_clocks = [x_fr*expobj.paq_rate for x_fr in x_ticks]  ## convert to paq clock dimension
    #     ax.set_xticks(x_clocks)
    #     ax.set_xticklabels(x_ticks)
    #     if not hide_xlabel:
    #         ax.set_xlabel('Frames')
    else:
        ax.set_xlabel('paq clock')
    ax.set_ylabel('Voltage')
    # ax.set_xlim([expobj.frame_start_time_actual, expobj.frame_end_time_actual])  ## this should be limited to the 2p acquisition duration only

    # set ylimits:
    if 'ylims' in kwargs:
        ax.set_ylim(kwargs['ylims'])
    else:
        ax.set_ylim([np.mean(expobj.lfp_signal) - 10, np.mean(expobj.lfp_signal) + 10])

    # set xlimits:
    if 'xlims' in kwargs and kwargs['xlims'] is not None:
        ax.set_xlim(kwargs['xlims'])


    # add title
    ax.set_title(
        '%s - %s %s %s' % (title, expobj.metainfo['exptype'], expobj.metainfo['animal prep.'], expobj.metainfo['trial']))

    # return None
    # if not 'fig' in kwargs.keys():
    #     ax.set_title(
    #         '%s - %s %s %s' % (title, expobj.metainfo['exptype'], expobj.metainfo['animal prep.'], expobj.metainfo['trial']))
    #
    # # options for showing plot or returning plot
    # if 'show' in kwargs.keys():
    #     plt.show() if kwargs['show'] else None
    # else:
    #     plt.show()

    # return fig, ax if 'fig' in kwargs.keys() else None
# LFP


# alloptical trial
### plot the location of all SLM targets, along with option for plotting the mean img of the current trial
@plot_piping_decorator(figsize=(5,5))
def plot_SLMtargets_Locs(expobj: AllOpticalTrial, targets_coords: list = None, background: np.ndarray = None, fig=None, ax=None, **kwargs):
    """
    plot SLM target coordinate locations

    :param expobj:
    :param targets_coords: ls containing (x,y) coordinates of targets to plot
    :param background:
    :param kwargs:
    :return:
    """

    # if 'fig' in kwargs.keys():
    #     fig = kwargs['fig']
    #     ax = kwargs['ax']
    # else:
    #     if 'figsize' in kwargs.keys():
    #         fig, ax = plt.subplots(figsize=kwargs['figsize'])
    #     else:
    #         fig, ax = plt.subplots()

    if background is None:
        background = np.zeros((expobj.frame_x, expobj.frame_y), dtype='uint16')
        ax.imshow(background, cmap='gray')
    else:
        ax.imshow(background, cmap='gray')

    colors = make_random_color_array(len(expobj.target_coords))
    if targets_coords is None:
        if len(expobj.target_coords) > 1:
            for i in range(len(expobj.target_coords)):
                for (x, y) in expobj.target_coords[i]:
                    ax.scatter(x=x, y=y, edgecolors=colors[i], facecolors='none', linewidths=2.0)
        else:
            if 'edgecolors' in kwargs.keys():
                edgecolors = kwargs['edgecolors']
            else:
                edgecolors = 'yellowgreen'
            for (x, y) in expobj.target_coords_all:
                ax.scatter(x=x, y=y, edgecolors=edgecolors, facecolors='none', linewidths=2.0)
    elif targets_coords:
        if 'edgecolors' in kwargs.keys():
            edgecolors = kwargs['edgecolors']
        else:
            edgecolors = 'yellowgreen'
        plot_coordinates(coords=targets_coords, frame_x=expobj.frame_x, frame_y=expobj.frame_y, edgecolors=edgecolors,
                            background=background, fig=fig, ax=ax)

    ax.margins(0)

    ax = add_scalebar(trialobj=expobj, ax=ax)

    fig.tight_layout()

    if 'title' in kwargs.keys():
        if kwargs['title'] is not None:
            ax.set_title(kwargs['title'])
        else:
            pass
    else:
        ax.set_title('SLM targets location - %s %s' % (expobj.metainfo['animal prep.'], expobj.metainfo['trial']))

    # if 'show' in kwargs.keys():
    #     if kwargs['show'] is True:
    #         fig.show()
    #     else:
    #         return fig, ax
    # else:
    #     fig.show()
    #
    # return fig, ax if 'fig' in kwargs.keys() else None


### plot entire trace of individual targeted cells as super clean subplots, with the same y-axis lims
def plot_photostim_traces(array, expobj: AllOpticalTrial, title='', y_min=None, y_max=None, x_label=None,
                          y_label=None, save_fig=None, **kwargs):
    """

    :param array:
    :param expobj:
    :param title:
    :param y_min:
    :param y_max:
    :param x_label:
    :param y_label:
    :param save_fig:
    :param kwargs:
        options include:
            hits: ls; a ls of 1s and 0s that is used to add a scatter point to the plot at stim_start_frames indexes at 1s
    :return:
    """
    # make rolling average for these plots
    w = 30
    array = [(np.convolve(trace, np.ones(w), 'valid') / w) for trace in array]

    len_ = len(array)
    fig, axs = plt.subplots(nrows=len_, sharex=True, figsize=(20, 3 * len_))
    for i in range(len(axs)):
        axs[i].plot(array[i], linewidth=1, color='black', zorder=2)
        if y_min != None:
            axs[i].set_ylim([y_min, y_max])
        for j in expobj.stim_start_frames:
            axs[i].axvline(x=j, c='gray', alpha=0.7, zorder=1)
        if 'scatter' in kwargs.keys():
            x = expobj.stim_start_frames[kwargs['scatter'][i]]
            y = [0] * len(x)
            axs[i].scatter(x, y, c='chocolate', zorder=3)
        if len_ == len(expobj.s2p_cell_targets):
            axs[i].set_title('Cell # %s' % expobj.s2p_cell_targets[i])
        if 'line_ids' in kwargs:
            axs[i].legend(['Target %s' % kwargs['line_ids'][i]], loc='upper left')


    axs[0].set_title((title + ' - %s' % len_ + ' cells'), loc='left', verticalalignment='top', pad=20,
                     fontsize=15)
    axs[0].set_xlabel(x_label)
    axs[0].set_ylabel(y_label)

    if save_fig is not None:
        plt.savefig(save_fig)

    fig.show()


@plot_piping_decorator(figsize=(10,6))
def plot_photostim_traces_overlap(array, expobj: AllOpticalTrial, exclude_id=[], y_spacing_factor=1, title='',
                                  x_axis='Time (seconds)', save_fig=None, fig=None, ax=None, **kwargs):
    '''
    :param array:
    :param expobj:
    :param spacing: a multiplication factor that will be used when setting the spacing between each trace in the final plot
    :param title:
    :param y_min:
    :param y_max:
    :param x_label:
    :param save_fig:
    :output: matplotlib plot
    '''
    # make rolling average for these plots
    # w = 30
    # array = np.asarray([(np.convolve(trace, np.ones(w), 'valid') / w) for trace in array])

    len_ = len(array)

    if 'fig' in kwargs.keys():
        fig = kwargs['fig']
        ax = kwargs['ax']
    else:
        if 'figsize' in kwargs.keys():
            fig, ax = plt.subplots(figsize=kwargs['figsize'])
        else:
            fig, ax = plt.subplots(figsize=[20, 10])

    for i in range(len_):
        if i not in exclude_id:
            if 'linewidth' in kwargs.keys():
                linewidth=kwargs['linewidth']
            else:
                linewidth=1
            ax.plot(array[i] + i * 40 * y_spacing_factor, linewidth=linewidth)
    for j in expobj.stim_start_frames:
        if j <= array.shape[1]:
            ax.axvline(x=j, c='gray', alpha=0.3)

    ax.set_xlim([0, expobj.n_frames-3000])

    ax.margins(0)
    # change x axis ticks to seconds
    if 'Time' in x_axis or 'time' in x_axis:
        # change x axis ticks to every 30 seconds
        labels = list(range(0, int(expobj.n_frames // expobj.fps), 30))
        ax.set_xticks(ticks=[(label * expobj.fps) for label in labels])
        ax.set_xticklabels(labels)
        ax.set_xlabel('Time (secs)')

        # labels = [item for item in ax.get_xticks()]
        # for item in labels:
        #     labels[labels.index(item)] = int(round(item / expobj.fps))
        # ax.set_xticklabels(labels)
        # ax.set_xlabel('Time (secs.)')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel(x_axis)

    if 'y_lim' in kwargs.keys():
        ax.set_ylim(kwargs['y_lim'])
    else:
        y_max = np.mean(array[-1] + len_ * 40 * y_spacing_factor) + 3 * np.mean(array[-1])
        ax.set_ylim(0, y_max)

    if save_fig is not None:
        plt.savefig(save_fig)

    ax.set_title((title + ' - %s' % len_ + ' cells'), horizontalalignment='center', verticalalignment='top', pad=20,
                 fontsize=10, wrap=True)

    # # finalize plot, set title, and show or return axes
    # if 'fig' in kwargs.keys():
    #     ax.set_title((title + ' - %s' % len_ + ' cells'), horizontalalignment='center', verticalalignment='top', pad=20,
    #                  fontsize=10, wrap=True)
    #     # ax.title.set_text((title + ' - %s' % len_ + ' cells'), wrap=True)
    #     return fig, ax
    # else:
    #     ax.set_title((title + ' - %s' % len_ + ' cells'), horizontalalignment='center', verticalalignment='top', pad=20,
    #                  fontsize=10, wrap=True)
    # if 'show' in kwargs.keys():
    #     if kwargs['show'] is True:
    #         fig.show()
    #     else:
    #         pass
    # else:
    #     fig.show()


### photostim analysis - PLOT avg over photostim. trials traces for the provided traces
@plot_piping_decorator(figsize=(5,5.5))
def plot_periphotostim_avg2(dataset, fps=None, legend_labels=None, colors=None, avg_with_std=False,
                            title='high quality plot', pre_stim_sec=None, ylim=None, fig=None, ax=None, **kwargs):

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
                print(f"|--- length mismatch in mean traces of datasets... {title}, shape0 {meanst.shape} and shape1 {meantraces[i - 1].shape}")
            if not std.shape == stdtraces[i - 1].shape:
                print(f"|--- length mismatch in std traces of datasets...{title}, shape0 {std.shape} and shape1 {stdtraces[i - 1].shape}")

    elif type(dataset) is not list or len(dataset) == 1:
        dataset = list(dataset)
        meanst = np.mean(dataset[0], axis=0)
        std = np.std(dataset[0], axis=0, ddof=1)
        meantraces.append(meanst)
        stdtraces.append(std)
        colors = ['black']
    else:
        AttributeError('please provide the data to plot in a ls format, each different data group as a ls item...')

    if 'xlabel' not in kwargs or kwargs['xlabel'] is None or 'frames' not in kwargs['xlabel'] or 'Frames' not in kwargs['xlabel']:
        ## change xaxis to time (secs)
        if fps is not None:
            if pre_stim_sec is not None:
                x_range = np.linspace(0, len(meantraces[0]) / fps, len(
                    meantraces[0])) - pre_stim_sec  # x scale, but in time domain (transformed from frames based on the provided fps)
                if 'xlabel' in kwargs.keys():
                    ax.set_xlabel(kwargs['xlabel'])
                else:
                    ax.set_xlabel('Time post stim (secs)')
            else:
                AttributeError('need to provide a pre_stim_sec value to the function call!')
        else:
            AttributeError('need to provide fps value to convert xaxis in units of time (secs)')
    elif 'frames' in kwargs['xlabel'] or 'Frames' in kwargs['xlabel']:
        x_range = range(len(meanst[0]))
        ax.set_xlabel('Frames')

    for i in range(len(meantraces)):
        if avg_with_std:
            if len(meantraces[i]) < len(x_range):  ## TEMP FIX
                mismatch = (len(x_range) - len(meantraces[i]))
                meantraces[i] = np.append(meantraces[i], [0] * mismatch)
                stdtraces[i] = np.append(stdtraces[i], [0] * mismatch)
                print(f'|------ adding {mismatch} zeros to mean and std-fill traces to make the arrays the same length, new length of plot array: {meantraces[i].shape} ')

            ax.plot(x_range, meantraces[i], color=colors[i], lw=2)
            ax.fill_between(x_range, meantraces[i] - stdtraces[i], meantraces[i] + stdtraces[i], alpha=0.15, color=colors[i])
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
        ax.set_ylim([ylim[0],ylim[1]])
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
@plot_piping_decorator(figsize=(5,5.5))
def plot_periphotostim_avg(arr: np.ndarray, expobj: AllOpticalTrial, pre_stim_sec=1.0, post_stim_sec=3.0, title='',
                           avg_only: bool = False, x_label=None, y_label=None, ax=None, pad=20, **kwargs):
    """
    plot trace across all stims
    :param arr: Flu traces to plot (will be plotted as individual traces unless avg_only is True) dimensions should be cells x stims x frames
    :param expobj: instance object of AllOpticalTrial
    :param pre_stim_sec: seconds of array to plot for pre-stim period
    :param post_stim_sec: seconds of array to plot for post-stim period
    :param title: title to use for plot
    :param avg_only: if True, only plot the mean trace calculated from the traces provided in arr
    :param x_label: x axis label
    :param y_label: y axis label
    :param ax: matplotlib.pyplot.Axes object
    :param kwargs:
        options include:
            'stim_duration': photostimulation duration in secs
            'y_lims': tuple, y min and max of the plot
            'edgecolor': str, edgecolor of the individual traces behind the mean trace
            'savepath': str, path to save plot to
            'show': bool = to show the plot or not
    :return: ls containing some items about the traces
    """

    fps = expobj.fps  # frames per second rate of the imaging data collection for the data to be plotted
    exp_prestim = expobj.pre_stim  # frames of pre-stim data collected for each trace for this expobj (should be same as what's under expobj.pre_stim_sec)
    if 'stim_duration' in kwargs.keys():
        stim_duration = kwargs['stim_duration']
    else:
        stim_duration = expobj.stim_dur / 1000  # seconds of stimulation duration

    x = list(range(arr.shape[1]))
    # x range in time (secs)
    x_time = np.linspace(0, arr.shape[1]/fps, arr.shape[1]) - pre_stim_sec  # x scale, but in time domain (transformed from frames based on the provided fps)

    len_ = len(arr)
    flu_avg = np.mean(arr, axis=0)

    # ax.margins(x=0.07)

    if 'alpha' in kwargs.keys():
        alpha = kwargs['alpha']
    else:
        alpha = 0.2

    if x_label is None or not 'Frames' in x_label or 'frames' in x_label:
        x = x_time  # set the x plotting range
        if x_label is not None:
            x_label = x_label + 'post-stimulation relative'

        if avg_only is True:
            # ax.axvspan(exp_prestim/fps, (exp_prestim + stim_duration + 1) / fps, alpha=alpha, color='plum', zorder = 3)
            ax.axvspan(0 - 1/fps, 0 + stim_duration + 1 / fps, alpha=alpha, color='plum', zorder = 3)  # note that we are setting 0 as the stimulation time
        else:
            ax.axvspan(0 - 1/fps, 0 - 1/fps + stim_duration, alpha=alpha, color='plum', zorder = 3)  # note that we are setting 0 as the stimulation time
    else:
        ax.axvspan(exp_prestim, exp_prestim + int(stim_duration*fps), alpha=alpha, color='tomato')

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
        if x_label is None or not 'Frames' in x_label or 'frames' in x_label:
            ax.set_xlim(-pre_stim_sec, stim_duration + post_stim_sec)  # remember that x axis is set to be relative to the stim time (i.e. stim is at t = 0)
        else:
            ax.set_xlim(exp_prestim - int(pre_stim_sec * fps), exp_prestim + int(stim_duration * fps) + int(post_stim_sec * fps) + 1)

    # spine options
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # set axis labels
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if 'savepath' in kwargs.keys():
        plt.savefig(kwargs['savepath'])

    if title is not None:
        ax.set_title((title + ' - %s' % len_ + ' traces'), horizontalalignment='center', verticalalignment='top',
                     pad=pad, fontsize=10, wrap=True)

    if avg_only:
        return flu_avg
# alloptical trial


