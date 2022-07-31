# TODO need to update this file to remove duplicates that have been refactored...

import sys
from pathlib import Path
from typing import Union

import matplotlib.pyplot
import numpy as np
import pandas as pd
import tifffile as tf
import os
import matplotlib.pyplot as plt
import csv
import math
import copy

from scipy import signal
from statsmodels import stats

from packerlabimaging.utils import io


def define_term(term: str):
    """
    Search terms dictionary and retrieve definition of input term if found.
    :param term: term to search for definition in dictionary.
    """
    try:
        from packerlabimaging.utils.terms_dictionary import terms_dictionary
        print(f"{term}:\t{terms_dictionary[term]}") if type(term) is str else print(
            'ERROR: please provide a string object as the key')
    except KeyError:
        print(f'input - {term} - not found in dictionary')


# report sizes of variables
def _sizeof_fmt(num, suffix='B'):
    """ by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified"""
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


def print_size_of(var):
    """
    Print size of input variable in current memory.
    :param var: input variable
    """
    print(_sizeof_fmt(sys.getsizeof(var)))


def print_size_vars():
    """
    Print size of all current variables in memory.
    """
    for name, size in sorted(((name, sys.getsizeof(value)) for name, value in locals().items()),
                             key=lambda x: -x[1])[:10]:
        print("{:>30}: {:>8}".format(name, _sizeof_fmt(size)))


def save_figure(fig: matplotlib.pyplot.Figure, save_path_full: str):
    """
    Save a matplotlib generated figure to the path provided. Also, makes necessary directories in the path.

    :param fig: matplotlib generated figure object
    :param save_path_full: path to save figure to.
    """
    print(f'\n\- saving figure to: {save_path_full}', end="\r")
    os.makedirs(os.path.dirname(save_path_full), exist_ok=True)
    fig.savefig(save_path_full)
    print(f'\n|- saved figure to: {save_path_full}')


def save_to_csv(df: pd.DataFrame, savepath: Path = None):
    """
    Save pandas dataframe to csv at savepath.

    :param df: dataframe to save.
    :param savepath: path to save to.
    """
    savepath.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(savepath)
    print(f"|- saved dataframe to {savepath}")


def filterDfBoolCol(df, true_cols=[], false_cols=[]):
    '''Filter indices in a pandas dataframe using logical operations
    on columns with Boolean values

        :param df:         -- dataframe
        :param true_cols:  -- columns where True should be filtered
        :param false_cols: -- columns where False should be filtered
    
    Outputs:
        indices of the dataframe where the logical operation is true
    '''
    if true_cols:
        true_rows = df[true_cols].all(axis='columns')

    if false_cols:
        false_rows = (~df[false_cols]).all(axis='columns')

    if true_cols and false_cols:
        filtered_df = df[true_rows & false_rows]
    elif true_cols:
        filtered_df = df[true_rows]
    elif false_cols:
        filtered_df = df[false_rows]

    return filtered_df.index


def findClosest(arr: Union[np.ndarray, list], input):
    """find the closest value in a list or 1-D array, to the given input

    :param arr: list or array of values
    :param input: key input to find closest value with in `arr`
    :return: value of the closest item from `arr`, index of closest item
    """
    if type(arr) == list: arr = np.array(arr)
    if not arr.ndim == 1: raise ValueError('input `arr` can only be of 1-dimension.')
    subtract = arr - input
    positive_values = abs(subtract)
    # closest_value = min(positive_values) + input
    index = np.where(positive_values == min(positive_values))[0][0]
    closest_value = arr[index]

    return closest_value, index


# calculates average over sliding window for an array
def moving_average(arr: Union[list, np.ndarray], n=4):
    """
    Return moving average over sliding window for an array
    
    :param arr: input array
    :param n: length of sliding window
    :return: array with moving average sliding window applied
    """
    if type(arr) == list: arr = np.array(arr)
    ret = np.cumsum(arr)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


# finding paths to files with a certain extension
def path_finder(umbrella, *args, is_dir=False):
    '''
    Returns the path to the single item in the umbrella folder containing the string names in each arg
    
    :param umbrella: parent directory to search inside
    :param args: key string names to search in umbrella directory
    :param is_dir: True if args is a collection of directories
    :return:
    '''
    # ls of bools, has the function found each argument?
    # ensures two folders / files are not found
    found = [False] * len(args)
    # the paths to the args
    paths = [None] * len(args)

    if is_dir:
        for root, dirs, files in os.walk(umbrella):
            for folder in dirs:
                for i, arg in enumerate(args):
                    if arg in folder:
                        assert not found[i], 'found at least two paths for {},' \
                                             'search {} to find conflicts' \
                            .format(arg, umbrella)
                        paths[i] = os.path.join(root, folder)
                        found[i] = True

    elif not is_dir:
        for root, dirs, files in os.walk(umbrella):
            for file in files:
                for i, arg in enumerate(args):
                    if arg in file:
                        assert not found[i], 'found at least two paths for {},' \
                                             'search {} to find conflicts' \
                            .format(arg, umbrella)
                        paths[i] = os.path.join(root, file)
                        found[i] = True

    print(paths)
    for i, arg in enumerate(args):
        if not found[i]:
            raise ValueError('could not find path to {}'.format(arg))

    return paths


def points_in_circle_np(radius, x0=0, y0=0):
    """
    Return an array of 2-D Euclidean coordinates for a circle centered on x0, y0 with the input radius.

    :param radius: radius of circle to be calculated
    :param x0: x-center of circle
    :param y0: y-center of circle
    """
    x_ = np.arange(x0 - radius - 1, x0 + radius + 1, dtype=int)
    y_ = np.arange(y0 - radius - 1, y0 + radius + 1, dtype=int)
    x, y = np.where((x_[:, np.newaxis] - x0) ** 2 + (y_ - y0) ** 2 <= radius ** 2)
    for x, y in zip(x_[x], y_[y]):
        yield x, y


def threshold_detect(signal: Union[np.ndarray, list], threshold):
    """
    Returns indexes where the input signal reaches above threshold.

    lloyd russell

    :param signal: input signal
    :param threshold: threhsold value
    :return:
    """
    thresh_signal = signal > threshold
    thresh_signal[1:][thresh_signal[:-1] & thresh_signal[1:]] = False
    frames = np.where(thresh_signal)
    return frames[0]


def listdirFullpath(directory: str, string: str = ''):
    """Return full path of all files in directory containing specified string

    :param directory:  path to directory
    :param string:  sequence to be found in file name
    :return: string
    """
    return [os.path.join(directory, file) for file in os.listdir(directory) if string in file]


def read_fiji(csv_path) -> np.ndarray:
    """Reads the csv file saved through plot z axis profile in Fiji

    :param csv_path: path to input csv file
    """

    data = []

    with open(csv_path, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for i, row in enumerate(spamreader):
            if i == 0:
                continue
            data.append(float(row[0].split(',')[1]))

    return np.array(data)


def intersect(lst1, lst2) -> list:
    """
    Find intersection of two input lists.

    :param lst1: list 1
    :param lst2: list 2
    """
    return list(set(lst1) & set(lst2))


# calculate distance between 2 points on a cartesian plane
def calc_distance_2points(p1: tuple, p2: tuple):
    """
    Uses the hypothenus method to calculate the straight line distance between two given points on a 2d cartesian plane.
    :param p1: point 1
    :param p2: point 2
    :return:
    """
    return math.hypot(p2[0] - p1[0], p2[1] - p1[1])