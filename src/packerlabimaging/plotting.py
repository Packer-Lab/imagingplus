# library of functions that are used for making various plots for all optical photostimulation/imaging experiments

# imports
import os
import numpy as np
import functools
import random
import matplotlib.pyplot as plt
import tifffile as tf

from .utils import normalize_dff

# wrapper for piping plots in and out of figures
def plot_piping_decorator(figsize=(5,5)):
    def plot_piping_decorator_(plotting_func):
        """
        Wrapper to help simplify creating plots from matplotlib.pyplot

        :param plotting_func: function to be wrapper
        :return: fig and ax objects if requested, and/or shows the figure as specified

        Examples:
        ---------
        1) in this example the fig and ax will be taken directly from the kwargs inside the inner wrapper
        @plot_piping_decorator
        def example_decorated_plot(title='', **kwargs):
            fig, ax = kwargs['fig'], kwargs['ax']  # need to add atleast this line to each plotting function's definition
            print(f'|-kwargs inside example_decorated_plot definition: {kwargs}')
            ax.plot(np.random.rand(10))
            ax.set_title(title)
        def example_decorated_plot(fig=None, ax=None, title='', **kwargs):
            print(f'|-kwargs inside example_decorated_plot definition: {kwargs}')
            ax.plot(np.random.rand(10))
            ax.set_title(title)

        example_decorated_plot()


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

## UTILITY FUNCS
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
### plot the location of all SLM targets, along with option for plotting the mean img of the current trial
@plot_piping_decorator(figsize=(5,5))
def plot_SLMtargets_Locs(expobj, targets_coords: list = None, background: np.ndarray = None, fig=None, ax=None, **kwargs):
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


# simple plot of the location of the given cell(s) against a black FOV
@plot_piping_decorator(figsize=(5,5))
def plot_cells_loc(expobj, cells: list, title=None, background: np.array = None, scatter_only: bool = False,
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

    # if 'show' in kwargs.keys():
    #     if kwargs['show'] is True:
    #         fig.show()
    #     else:
    #         pass
    # else:
    #     fig.show()
    #
    # return fig, ax if 'fig' in kwargs.keys() else None


# plot to show s2p ROIs location, colored as specified
def s2pRoiImage(expobj, save_fig: str = None):
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
    plt.margins(x=0, y=0)
    plt.gca().invert_yaxis()
    # plt.gca().invert_xaxis()
    # fig.show()

    plt.suptitle(f"{expobj.metainfo['animal prep.']} {expobj.metainfo['trial']} - s2p nontargets (blue), exclude (yellow), targets (red); target_areas (white)",
                 y=0.97, fontsize=7)
    plt.show()
    save_figure(fig, save_path_suffix=f"{save_fig}") if save_fig else None

# (full) plot individual cell's flu or dFF trace, with photostim. timings for that cell
def plot_flu_trace(expobj, cell, x_lims=None, slm_group=None, to_plot='raw', figsize=(20, 3), linewidth=0.10, show=True):
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
def plotMeanRawFluTrace(expobj, stim_span_color='white', stim_lines: bool = True, title='raw Flu trace', x_axis='time', shrink_text=1,
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


# imshow gray plot for a single frame tiff
def plot_single_tiff(tiff_path: str, title: str = None, frame_num: int = 0):
    """
    plots an image of a single tiff frame after reading using tifffile.
    :param tiff_path: path to the tiff file
    :param title: give a string to use as title (optional)
    :return: imshow plot
    """
    stack = tf.imread(tiff_path, key=frame_num)
    plt.imshow(stack, cmap='gray')
    if title is not None:
        plt.suptitle(title)
    else:
        plt.suptitle('frame num: %s' % frame_num)
    plt.show()
    return stack


