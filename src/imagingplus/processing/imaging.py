"""
Functions for various imaging related processing tasks. 

"""

import os
from typing import Union

import numpy as np
from imagingplus.main.core import ImagingTrial
from imagingplus.utils.utils import points_in_circle_np


def normalize_dff(arr: np.ndarray, normalize_pct: int = 20, normalize_val=None):
    """
    Normalize given array (cells x time) to the mean of the fluorescence values below given threshold. Threshold
    will refer to the that lower percentile of the given trace. Calculate dFF of traces.

    :param arr: numpy array of fluorescence values (cells x time)
    :param normalize_pct: percentile to normalize each cell trace to, default = 20 (20th percentile)
    :param normalize_val: value to normalize each cell trace to
    :return: normalized array of fluorescence values (cells x time)
    """

    if arr.ndim == 1:
        if normalize_val is None:
            a = np.percentile(arr, normalize_pct)
            mean_ = arr[arr < a].mean()
        else:
            mean_ = normalize_val
        new_array = ((arr - mean_) / mean_) * 100
        if np.isnan(new_array).any() == True:
            Warning('Cell (unknown) contains nan, normalization factor: %s ' % mean_)

    else:
        new_array = np.empty_like(arr)
        for i in range(len(arr)):
            if normalize_val is None:
                a = np.percentile(arr[i], normalize_pct)
            else:
                a = normalize_val
            mean_ = np.mean(arr[i][arr[i] < a])
            new_array[i] = ((arr[i] - mean_) / abs(mean_)) * 100

            if np.isnan(new_array[i]).any() == True:
                print('Warning:')
                print('Cell %d: contains nan' % (i + 1))
                print('      Mean of the sub-threshold for this cell: %s' % mean_)

    return new_array



def calcAnnulus(trial: ImagingTrial, coord: tuple, inner_distance: float = 10, outer_distance: float = 15):
    """
    Creates an annulus around a given coordinate (coord) of the specified diameter that is considered the exclusion zone around the coordinate.
    Also return a slice object for the annulus around the

    # use the SLM targets exclusion zone areas as the annulus around each SLM target
    # -- for each SLM target: create a numpy array slice that acts as an annulus around the target

    :param coord:
    :param distance: distance (in microns) from the edge of the spiral to extend the target exclusion zone
    :return:
    """

    radius_px_inner = int(inner_distance) / trial.imparams.pix_sz_x
    inner_area = tuple([item for item in points_in_circle_np(radius_px_inner, x0=coord[0], y0=coord[1])])

    # create annulus by subtracting SLM spiral target pixels
    radius_px_outer = int(outer_distance) / trial.imparams.pix_sz_x

    outer_area = tuple([item for item in points_in_circle_np(radius_px_outer, x0=coord[0], y0=coord[1])])
    area_annulus = np.array([list(coord_) for i, coord_ in enumerate(outer_area) if coord_ not in inner_area])

    coords = []
    for coord in area_annulus:
        if (0 < coord[0] < trial.imparams.frame_x) and (0 < coord[1] < trial.imparams.frame_y):
            coords.append(coord)

    area_annulus_final = np.array(coords)


    # testing the areas created
    # import matplotlib.pyplot as plt
    # frame = np.zeros(shape=(trial.imparams.frame_x, trial.imparams.frame_y), dtype=int)

    # add to frame_array towards creating a plot
    # for x, y in outer_area:
    #     frame[x, y] = 10

    # for x, y in inner_area:
    #     frame[x, y] = 0

    # for x, y in area_annulus:
    #     frame[x, y] = 10

    # plt.figure(figsize=(4, 4), dpi=600)
    # plt.imshow(frame, cmap='bwr')
    # plt.colorbar()
    # plt.show()

    annulus_slice_obj_ = np.s_[area_annulus_final[:, 0], area_annulus_final[:, 1]]

    return area_annulus_final, annulus_slice_obj_


def _collect_annulus_flu(trial, annulus_slice_obj):
    """
    Read in raw registered tiffs, then use slice object to collect individual targets' annulus raw traces directly from the tiffs

    :param annulus_slice_obj: list of len(n_targets) containing the numpy slice object for SLM targets
    """

    print('\n\ncollecting raw Flu traces from SLM target coord. areas from registered TIFFs')

    # read in registered tiff
    reg_tif_folder = trial.s2p_path + '/reg_tif/'
    reg_tif_list = os.listdir(reg_tif_folder)
    reg_tif_list.sort()
    start = trial.curr_trial_frames[0] // 2000  # 2000 because that is the batch size for suite2p run
    end = trial.curr_trial_frames[1] // 2000 + 1

    mean_img_stack = np.zeros([end - start, trial.frame_x, trial.frame_y])
    # collect mean traces from target areas of each target coordinate by reading in individual registered tiffs that contain frames for current trial
    targets_annulus_traces = np.zeros([len(trial.slmtargets_ids), (end - start) * 2000], dtype='float32')
    for i in range(start, end):
        tif_path_save2 = trial.s2p_path + '/reg_tif/' + reg_tif_list[i]
        with tf.TiffFile(tif_path_save2, multifile=False) as input_tif:
            print('\t reading tiff: %s' % tif_path_save2)
            data = input_tif.asarray()

        target_annulus_trace = np.zeros([len(trial.target_coords_all), data.shape[0]], dtype='float32')
        for idx, coord in enumerate(trial.target_coords_all):
            # target_areas = np.array(trial.target_areas)
            # x = data[:, target_areas[coord, :, 1], target_areas[coord, :, 0]]
            x = data[:, annulus_slice_obj[idx][0], annulus_slice_obj[idx][1]]
            target_annulus_trace[idx] = np.mean(x, axis=1)

        targets_annulus_traces[:, (i - start) * 2000: ((i - start) * 2000) + data.shape[
            0]] = target_annulus_trace  # iteratively write to each successive segment of the targets_trace array based on the length of the reg_tiff that is read in.

    # final part, crop to the exact frames for current trial
    trial.raw_SLMTargets_annulus = targets_annulus_traces[:,
                                   trial.curr_trial_frames[0] - start * 2000: trial.curr_trial_frames[1] - (
                                           start * 2000)]

    return trial.raw_SLMTargets_annulus
