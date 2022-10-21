# library of convenience plotting funcs that are used for making various plots for all optical photostimulation/imaging experiments

# imports
from typing import Union

import numpy as np
import functools
import random
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from imagingplus.main.core import ImagingTrial, SingleImage

# global plotting params
params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
plt.rcParams.update(params)
sns.set()
sns.set_style('white')


# %% UTILITY FUNCS

def plotting_decorator(figsize=(3, 3), nrows=1, ncols=1, dpi=300):
    """
    Wrapper for piping plots in and out of figures
    """
    def plotting_decorator(plotting_func):
        """
        Wrapper to help simplify creating plots from matplotlib.pyplot

        :param plotting_func: plotting convenience function to be wrapped
        :return: fig and ax objects if requested, and/or shows the figure as specified

        Examples:
        ---------
        1) in this example the fig and ax will be taken directly from the kwargs inside the inner wrapper
        >>> @plotting_decorator
        >>> def example_decorated_plot(title='', **kwargs):
        >>>     fig, ax = kwargs['fig'], kwargs['ax']  # need to add atleast this line to each plotting function's definition
        >>>     print(f'|-kwargs inside example_decorated_plot definition: {kwargs}')
        >>>     ax.plot(np.random.rand(10))
        >>>     ax.set_title(title)

        >>> example_decorated_plot()

        """

        @functools.wraps(plotting_func)
        def inner(*args, **kwargs):
            return_fig_obj = False

            # set number of rows, cols and figsize
            if 'nrows' in [*kwargs]:
                nrows_ = kwargs['nrows']
            else:
                nrows_ = nrows

            if 'ncols' in [*kwargs]:
                ncols_ = kwargs['ncols']
            else:
                ncols_ = ncols

            if 'figsize' in [*kwargs]:
                figsize_ = kwargs['figsize']
            else:
                figsize_ = figsize

            # create or retrieve the fig, ax objects --> end up in kwargs to use into the plotting func call below
            if 'fig' in [*kwargs] and 'ax' in [*kwargs]:
                if kwargs['fig'] is None or kwargs['ax'] is None:
                    # print('\-creating fig, ax [1]')
                    kwargs['fig'], kwargs['axs'] = plt.subplots(nrows=nrows_, ncols=ncols_, figsize=figsize_)
                    kwargs['fig'].tight_layout(pad=1.8)
                else:
                    pass
            elif ncols_ > 1 or nrows_ > 1:
                kwargs['fig'], kwargs['axs'] = plt.subplots(nrows=nrows_, ncols=ncols_, figsize=figsize_, dpi=dpi)
                kwargs['fig'].tight_layout(pad=1.8)
            else:
                kwargs['fig'], kwargs['ax'] = plt.subplots(figsize=figsize_, dpi=dpi)
                kwargs['fig'].tight_layout(pad=1.8)

            print(f'\- executing plotting function: {plotting_func.__name__}')
            res = plotting_func(
                **kwargs)  # these kwargs are the original kwargs defined at the respective plotting_func call + any additional kwargs defined in inner()

            kwargs['ax'].set_title(kwargs['title'], wrap=True) if 'title' in [*kwargs] else None


            if 'show' in [*kwargs]:
                if kwargs['show'] is True:
                    kwargs['fig'].show()
                    return res
                else:
                    if res is not None:
                        return (kwargs['fig'], kwargs['ax'], res)
                    else:
                        return (kwargs['fig'], kwargs['ax'])
            else:
                kwargs['fig'].show()
                return res

        return inner

    return plotting_decorator


# custom colorbar for heatmaps
from matplotlib.colors import LinearSegmentedColormap


def _make_colormap(seq):
    """Return a LinearSegmentedColormap
    :param seq: a sequence of floats and RGB-tuples. The floats should be increasing and in the interval (0,1).
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


def _get_random_color(pastel_factor: float = 0.5):
    """
    Generate an array of random colors
    """
    return [(x + pastel_factor) / (1.0 + pastel_factor) for x in [random.uniform(0, 1.0) for i in [1, 2, 3]]]


def _color_distance(c1, c2):
    """
    Calculate distance between two colors (input in RGB values).
    """
    return sum([abs(x[0] - x[1]) for x in zip(c1, c2)])


def _generate_new_color(existing_colors: Union[list, set], pastel_factor=0.5):
    """
    Generate a new color based on existing colors as input.
    """
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
    :return: list of colors in RGB values
    """
    colors = []
    for i in range(0, n_colors):
        colors.append(_generate_new_color(colors, pastel_factor=0.2))
    return colors


def _add_scalebar(trialobj: Union[ImagingTrial, SingleImage], ax: plt.Axes, **kwargs):
    """add scalebar to the image being plotted on the a single matplotlib.axes.Axes object using the TwoPhotonImaging object information.
    Option to specify scale bar um length to add to plot.

    :param trialobj: ImagingTrial or SingleImage; object associated with input image.
    :param ax: plot object to add scale bar on
    :param kwargs:
        :scalebar_um: int; size of scalebar to plot on image (in um); must provide trialobj parameter.
        :right_offset: int; move scalebar left/right
        :bottom_offset: int; move scalebar up/down
        :color: str; color of scalebar to show on plot
        :lw: int, float; linewidth of scalebar to show on plot
    :return:

    """

    # REMOVED FOR THE TIME BEING WHILE FIGURING OUT HOW TO ACCEPT CHILD CLASSES OF VALID CLASSES
    # if type(trialobj) not in [TwoPhotonImagingTrial, AllOpticalTrial]:
    #     raise ObjectClassError(function='_add_scalebar', valid_class=[TwoPhotonImagingTrial, AllOpticalTrial], invalid_class=type(trialobj))
    # else:
    assert trialobj.imparams, 'could not find imaging metadata (.imparams) for trialobj.'
    scalebar_um = 100 if 'scalebar_um' not in kwargs else kwargs['scalebar_um']
    assert type(scalebar_um) == int, 'scale bar length must be provided as integer (units: microns)'
    numpx = scalebar_um / trialobj.imparams.pix_sz_x

    lw = 5 if 'lw' not in [*kwargs] else kwargs['lw']
    color = 'white' if 'color' not in [*kwargs] else kwargs['color']
    right_offset = 50 if 'right_offset' not in [*kwargs] else kwargs['right_offset']
    bottom_offset = 50 if 'bottom_offset' not in [*kwargs] else kwargs['bottom_offset']

    print(f'\- adding scalebar of length {scalebar_um} microns.')

    ax.plot(np.linspace(trialobj.imparams.frame_x - right_offset - numpx, trialobj.imparams.frame_x - right_offset, 40),
            [trialobj.imparams.frame_y - bottom_offset] * 40, color=color, lw=lw, solid_capstyle='butt')

    return ax


image_frame_ops = {
    'axes.spines.left': False,
    'axes.spines.bottom': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.frameon': False,
    'figure.subplot.wspace': .01,
    'figure.subplot.hspace': .01,
    'figure.figsize': (18, 13),
    'ytick.major.left': False,
    'xtick.major.bottom': False}


# Figure Style settings for notebook.
def image_frame_options(fig, ax):
    """
    Matplotlib settings for plotting an image.
    """
    # mpl.pyplot.rcdefaults()

    # mpl.rcParams.update({
    #     'axes.spines.left': False,
    #     'axes.spines.bottom': False,
    #     'axes.spines.top': False,
    #     'axes.spines.right': False,
    #     'legend.frameon': False,
    #     'figure.subplot.wspace': .01,
    #     'figure.subplot.hspace': .01,
    #     'figure.figsize': (18, 13),
    #     'ytick.major.left': False,
    #     'xtick.major.bottom': False
    #
    # })

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

    ax.spines.left = False
    ax.spines.bottom = False
    ax.spines.top = False
    ax.spines.right = False
    fig.subplots_adjust(hspace=0.01)
    fig.subplots_adjust(wspace=0.01)
    ax.set_xticks([])
    ax.set_yticks([])
    # ax.xaxis.set_major_locator([])
    # ax.yaxis.set_major_locator([])
    # ax.xaxis.set_minor_locator([])
    # ax.yaxis.set_minor_locator([])


def dataplot_frame_options():
    """
    Matplotlib frame-level settings for plotting a dataplot.
    """
    import matplotlib as mpl

    mpl.rcParams.update({
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.fontsize': 'x-large',
        'axes.labelsize': 'x-large',
        'axes.titlesize': 'x-large',
        'xtick.labelsize': 'x-large',
        'ytick.labelsize': 'x-large',
        'legend.frameon': False,
        'figure.subplot.wspace': .01,
        'figure.subplot.hspace': .01,
    })
    sns.set()
    sns.set_style('white')


def dataplot_ax_options(ax, **kwargs):
    """
    Matplotlib axis-level settings for plotting a dataplot.
    :param
        **kwargs:
            :x_axis: x axis label, if specify Time or time in x_axis then convert x_axis to time domain
            :collection_hz: cellsdata collection rate (in Hz)
            :x_tick_secs: interval (in secs) for plotting x ticks when converting x axis to time domain
            :xlims: set xlimits for plot
            :ylims: set ylimits for plot

    """
    if ax:
        ax.margins(0)
        ax.grid(True)

        # set x and y axis limits
        if 'ylims' in [*kwargs]: ax.set_ylim([ylim for ylim in kwargs['ylims']])

        # set x_axis label
        if 'x_axis' in [*kwargs]:
            x_axis = kwargs['x_axis']
            # change x-axis to time (secs) if time is requested
            if ('time' in x_axis or 'Time' in x_axis) and 'collection_hz' in [*kwargs]:
                if 'xlims' in [*kwargs]:
                    ax.set_xlim([xlim * kwargs['collection_hz'] for xlim in kwargs['xlims']])
                    x_tick_secs = int((kwargs['xlims'][1] - kwargs['xlims'][0]) // 2) if 'x_tick_secs' not in [*kwargs] else kwargs['x_tick_secs']
                    labels = list(np.arange(kwargs['xlims'][0], kwargs['xlims'][1], x_tick_secs)) if (kwargs['xlims'][1] - kwargs['xlims'][0]) > x_tick_secs else list(np.arange(kwargs['xlims'][0], kwargs['xlims'][1]))
                else:
                    x_tick_secs = 30 if 'x_tick_secs' not in [*kwargs] else kwargs['x_tick_secs']
                    start, end = 0, ax.get_xlim()[1]
                    labels = list(
                        range(0, int(end // kwargs['collection_hz']), x_tick_secs))

                x_tick_locations = [(label * kwargs['collection_hz']) for label in labels]
                ax.set_xticks(ticks=x_tick_locations)

                ax.set_xticklabels(labels)
                ax.set_xlabel('Time (secs)')
            else:
                if 'xlims' in [*kwargs]: ax.set_xlim([xlim for xlim in kwargs['xlims']])
                ax.set_xlabel(x_axis)

        # set y_axis label
        ax.set_ylabel(kwargs['y_axis']) if 'y_axis' in [*kwargs] else None


    else:
        pass


def heatmap_options():
    """
    Matplotlib settings for plotting a heatmap.
    """
    jet = mpl.cm.get_cmap('jet')
    jet.set_bad(color='k')


# %% PLOTTING UTILITIES

@plotting_decorator(figsize=(5, 5))
def plot_coordinates(coords: list, frame_x: int, frame_y: int, background: np.ndarray = None, fig=None, ax=None,
                     **kwargs):
    """
    Plot coordinate locations using matplotlib's imshow function.

    :param coords: ls containing (x,y) coordinates of targets to plot
    :param frame_x: x-axis size of background plot
    :param frame_y: y-axis size of background plot
    :param background: np.array on which to plot coordinates, default is black background (optional)
    :param fig: matplotlib figure object to plot
    :param ax: matplotlib axis object to plot
    :param kwargs:
        :edgecolors: edgecolor of coordinates that are plotted
    """

    if background is None:
        background = np.zeros((frame_x, frame_y), dtype='uint16')
        ax.imshow(background, cmap='gray')
    else:
        ax.imshow(background, cmap='gray')

    if 'edgecolors' in [*kwargs]:
        edgecolors = kwargs['edgecolors']
    else:
        edgecolors = 'yellowgreen'
    for (x, y) in coords:
        ax.scatter(x=x, y=y, edgecolors=edgecolors, facecolors='none', linewidths=2.0)

    ax.margins(0)
    fig.tight_layout()

