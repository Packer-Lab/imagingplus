import os
import time
from pathlib import Path
from typing import Dict, Any, Optional, List, Union, Tuple

import numpy as np
import pandas as pd
import suite2p
# todo add type hinting controls for expobj
import matplotlib.pyplot as plt
from imagingplus.main.subcore import CellAnnotations, ImagingData
from suite2p import run_s2p
import tifffile as tf

from imagingplus.utils.classes import UnavailableOptionError
from imagingplus.utils.images import ImportTiff, SaveDownsampledTiff, WriteTiff, make_tiff_stack


# suite2p methods
def s2pRun(expobj, trialsSuite2P: Union[list, str] = 'all'):
    """run suite2p for an Experiment object, using trials specified in current experiment object, using the attributes
    determined directly from the experiment object.

    :param expobj: a imagingplus.main.Experiment object
    :param trialsSuite2P: list of trialIDs from experiment to use in running suite2p
    """

    _suite2p_save_path = expobj.saveDir + '/suite2p/'  #: default location to save Suite2p output results of current experiment

    expobj.Suite2p.trials = trialsSuite2P if trialsSuite2P != 'all' else expobj.Suite2p.trials

    tiffs_paths_to_use_s2p = []
    for trial in expobj.Suite2p.trials:
        tiffs_paths_to_use_s2p.append(expobj.TrialsInformation[trial]['tiff_path'])

    expobj.Suite2p.tiffs_paths_to_use_s2p = tiffs_paths_to_use_s2p

    # setup db dict
    expobj.Suite2p.ops['batch_size'] = int(expobj.Suite2p.ops['batch_size'])
    db = {'fs': float(expobj.Suite2p.ops['fs']), 'diameter': expobj.Suite2p.ops['diameter'],
          'batch_size': expobj.Suite2p.ops['batch_size'],
          'nimg_init': expobj.Suite2p.ops['batch_size'], 'nplanes': expobj.Suite2p.ops['nplanes'],
          'nchannels': expobj.Suite2p.ops['nchannels'],
          'tiff_list': list(tiffs_paths_to_use_s2p), 'data_path': expobj.dataPath, 'save_folder': _suite2p_save_path
          }

    expobj.Suite2p.db = db

    print('RUNNING Suite2p:')
    print(f'-- .Suite2p.db: ')  # \n\t{expobj.Suite2p.db}\n\n')

    t1 = time.time()
    opsEnd = run_s2p(ops=expobj.Suite2p.ops, db=expobj.Suite2p.db)
    t2 = time.time()
    print('Total time to run suite2p was {}'.format(t2 - t1))

    # update expobj.Suite2p.ops and db
    expobj.Suite2p.ops = opsEnd

    expobj.Suite2p.s2pResultExists = True
    expobj.Suite2p.s2pResultsPath = _suite2p_save_path

    print(f'\n\nCompleted Suite2p run. \n\t Results are saved to: {expobj.Suite2p.s2pResultsPath}')

    expobj.save()


def s2p_loader(s2p_path, subtract_neuropil=True, neuropil_coeff=0.7):
    """
    Load suite2p results from file path at `s2p_path`.

    :param s2p_path: file path from which to load suite2p results
    :param subtract_neuropil: if True, use neuropil cooefficient to subtract neuropil signal from raw fluorescence signal.
    :param neuropil_coeff: neuropil coefficient to subtract neuropil signal.
    :return: Raw signal array (optionally neuropil corrected), deconvolved spikes, stat output, neuropil signal array
    """
    found_stat = False

    for root, dirs, files in os.walk(s2p_path):

        for file in files:

            if file == 'F.npy':
                all_cells = np.load(os.path.join(root, file), allow_pickle=True)
            elif file == 'Fneu.npy':
                neuropil = np.load(os.path.join(root, file), allow_pickle=True)
            elif file == 'iscell.npy':
                is_cells = np.load(os.path.join(root, file),
                                   allow_pickle=True)[:, 0]
                is_cells = np.ndarray.astype(is_cells, 'bool')
                print('Loading {} traces labelled as cells'
                      .format(sum(is_cells)))
            elif file == 'spks.npy':
                spks = np.load(os.path.join(root, file), allow_pickle=True)
            elif file == 'stat.npy':
                stat = np.load(os.path.join(root, file), allow_pickle=True)
                found_stat = True

    if not found_stat:
        raise FileNotFoundError('Could not find stat, '
                                'this is likely not a suit2p folder, path: ', s2p_path)

    # for i, s in enumerate(stat):
    #     s['original_index'] = i

    all_cells = all_cells[is_cells, :]
    neuropil = neuropil[is_cells, :]
    spks = spks[is_cells, :]
    stat = stat[is_cells]

    if not subtract_neuropil:
        return all_cells, spks, stat, neuropil

    else:
        print('Subtracting neuropil with a coefficient of {}'
              .format(neuropil_coeff))
        neuropil_corrected = all_cells - neuropil * neuropil_coeff
        return neuropil_corrected, spks, stat, neuropil


# utility funcs
# quick workaround (patch of suite2p code) because of suite2p internal error for these methods in ROI class
def from_stat_dict(stat: Dict[str, Any]) -> suite2p.ROI:
    """
    Overloading function from suite2p source code.

    :param stat: suite2p results output stat dict
    """
    return suite2p.ROI(ypix=stat['ypix'], xpix=stat['xpix'], lam=stat['lam'], med=stat['med'], do_crop=False)


def to_array_(self, Ly: int, Lx: int) -> np.ndarray:
    """Returns a 2D boolean array of shape (Ly x Lx) indicating where the roi is located.
    Overloading function from suite2p source code.
    """
    arr = np.zeros((Ly, Lx), dtype=float)
    arr[self.ypix, self.xpix] = 1
    return arr


def stats_dicts_to_3d_array_(stat_dict, output_ops):
    """
    Overloading function from suite2p source code.
    """
    arrays = []
    for i, stat in enumerate(stat_dict):
        array = to_array_(from_stat_dict(stat=stat), Ly=output_ops['Ly'], Lx=output_ops['Lx'])
        array *= i + 1
        arrays.append(array)
    return np.stack(arrays)


# todo have suite2p run return a cells annotations (and maybe even imaging cellsdata set??) that can get added to trials' cells and imdata attr's immediately

class Suite2pExperiment:
    """used to run and further process suite2p processed cellsdata, and analysis associated with suite2p processed cellsdata."""

    # default ops dict for suite2p experiment run
    ops = {
        'batch_size': 2000,  # reduce if running out of RAM
        'delete_bin': True,  # whether to delete binary file after processing
        # main settings
        'nplanes': 1,  # each tiff has these many planes in sequence
        'nchannels': 1,  # each tiff has these many channels per plane
        'functional_chan': 1,  # this channel is used to extract functional ROIs (1-based)
        'diameter': 12,
        # this is the main parameter for cell detection, 2-dimensional if Y and X are different (e.g. [6 12])
        # 'tau': 1.26,  # this is the main parameter for deconvolution (1.25-1.5 for gcamp6s)
        # 'fs': 30,  # sampling rate (total across planes)
        # output settings
        'save_mat': True,  # whether to save output as matlab files
        'combined': True,  # combine multiple planes into a single result /single canvas for GUI
        # parallel settings
        'num_workers': 0,  # 0 to select num_cores, -1 to disable parallelism, N to enforce value
        'num_workers_roi': 0,  # 0 to select number of planes, -1 to disable parallelism, N to enforce value
        # registration settings
        'do_registration': True,  # whether to register cellsdata
        'nimg_init': 200,  # subsampled key_frames for finding reference image
        'maxregshift': 0.1,  # max allowed registration shift, as a fraction of frame max(width and height)
        'align_by_chan': 1,  # when multi-channel, you can align by non-functional channel (1-based)
        'reg_tif': False,  # whether to save registered tiffs
        'subpixel': 10,  # precision of subpixel registration (1/subpixel steps)
        # cell detection settings
        'connected': True,  # whether or not to keep ROIs fully connected (set to 0 for dendrites)
        'navg_frames_svd': 5000,  # max number of binned key_frames for the SVD
        'nsvd_for_roi': 1000,  # max number of SVD components to keep for ROI detection
        'max_iterations': 20,  # maximum number of iterations to do cell detection
        'ratio_neuropil': 6.,  # ratio between neuropil basis size and cell radius
        'ratio_neuropil_to_cell': 3,  # minimum ratio between neuropil radius and cell radius
        'tile_factor': 1.,  # use finer (>1) or coarser (<1) tiles for neuropil estimation during cell detection
        'threshold_scaling': 1.,  # adjust the automatically determined threshold by this scalar multiplier
        'max_overlap': 0.75,  # cells with more overlap than this get removed during triage, before refinement
        'inner_neuropil_radius': 2,  # number of pixels to keep between ROI and neuropil donut
        'outer_neuropil_radius': np.inf,  # maximum neuropil radius
        'min_neuropil_pixels': 350,  # minimum number of pixels in the neuropil
        # deconvolution settings
        'baseline': 'maximin',  # baselining mode
        'win_baseline': 60.,  # window for maximin
        'sig_baseline': 10.,  # smoothing constant for gaussian filter
        'prctile_baseline': 8.,  # optional (whether to use a percentile baseline)
        'neucoeff': .7,  # neuropil coefficient
    }

    @classmethod
    def update_ops(cls, new_ops_entries: dict):
        """
        Add or update existing ops fields before running suite2p pipeline.

        :param new_ops_entries: dict specificying new or updated ops fields.
        """
        print(f'Updating {[*new_ops_entries]} keys in ops prior to Suite2p run.')
        for key in new_ops_entries:
            cls.ops[key] = new_ops_entries[key]

    def __init__(self, trialsTiffsSuite2p: dict, s2pResultsPath: Optional[str] = None):
        """
        Submodule for connecting an Experiment to Suite2p for the specified trials/data_tiffs.

        :param trialsTiffsSuite2p:
        :param s2pResultsPath: optional, if suite2p already run then here provide path to plane0 folder of existing suite2p results
        :param subtract_neuropil:
        :param dataPath:
        """
        self.__s2pResultExists = False
        self.bad_frames: list = []  #: list of key_frames to discard during Suite2p ROI detection procedure
        print(f"\- ADDING .Suite2p module to Experiment object ... ", end='\r')

        ## initialize attr's
        self.n_frames: int = 0  # total number of imaging key_frames in the Suite2p run

        self.cell_id = []  # cell ROI # id per plane
        self.n_units = []  # num of ROIs in suite2p result per plane
        self.cell_med = []  # y, x pixel location of each cell in the imaging frame
        self.cell_x = []  # x pixel location of cell
        self.cell_y = []  # y pixel location of cell
        self.raw = []  # [cells x num key_frames], raw traces for each cell from suite2p
        self.stat = []  # stat output from suite2p
        self.spks = []  # deconvolved spikes output from suite2p
        self.neuropil = []  # neuropil output from suite2p
        self.mean_img = []  # mean img output from suite2p
        self.mean_imgE = []  # mean img Enhanced output from suite2p
        self.radius = []  # cell radius from suite2p
        self.xoff = []  # motion correction info
        self.yoff = []  # motion correction info

        self.neuropil_coeff = 0

        if s2pResultsPath is None:
            # initialize needed variables and attr's for future calling of s2pRun
            self.__s2pResultsPath = None
            self.s2pResultExists = False
        else:
            self.__s2pResultsPath = s2pResultsPath
            try:
                self._retrieveSuite2pData(self.s2pResultsPath, neuropil_coeff=self.neuropil_coeff)
            except Exception:
                raise Exception(
                    f'Something went wrong while trying to load suite2p processed cellsdata from: {s2pResultsPath}')
            self.s2pResultExists = True

        # set trials to run together in suite2p for Experiment
        self.trials = []
        tiff_paths_to_use_s2p: dict = trialsTiffsSuite2p
        for trial, path in tiff_paths_to_use_s2p.items():
            if path is not None:
                if self.s2pResultExists and hasattr(self, 'output_ops'):
                    if path not in self.output_ops['filelist']:
                        print(f'\tWARNING: Path for {trial}, ({path}) was not found in suite2p results output_ops filelist.')
                    self.trials.append(trial)
            else:
                tiff_paths_to_use_s2p.pop(trial)

        self.tiff_paths_to_use_s2p: dict = tiff_paths_to_use_s2p
        assert len(self.trials) > 0, "Failed adding trials to suite2p module. \n" \
                                     "\t - Ensure that correct trials are provided." \
                                     "\t - Ensure that file paths provided for each trial correspond to filepaths in output_ops." \

        self.db: dict = {}

        # finish
        print(f"|- ADDED .Suite2p module to Experiment object. ")

    @property
    def s2pResultsPath(self):
        """
        Path to the saved suite2p results output.
        """
        return self.__s2pResultsPath

    @s2pResultsPath.setter
    def s2pResultsPath(self, val):
        """
        Setter for path to the saved suite2p results output.
        :param val: file path
        """
        self.__s2pResultsPath = val


    @property
    def s2pResultExists(self):
        """
        Return true if s2p results exist for current trial.
        """
        return self.__s2pResultExists

    @s2pResultExists.setter
    def s2pResultExists(self, val):
        """
        Set if s2p results exist for current trial.
        :param val:
        """
        self.__s2pResultExists = val

    def __repr__(self):
        if self.s2pResultExists:
            return f'Suite2p Results (Experiment level) Object, containing trials: \n\t{self.trials}'
        else:
            return f'Suite2p Results (Experiment level) Object, containing trials: \n\t{self.trials}. No Suite2p results loaded.'

    # noinspection PyTypeChecker
    def _retrieveSuite2pData(self, s2p_path: str = None, neuropil_coeff: float = 0.7):
        """processing of suite2p cellsdata from the current t-series
        :param s2p_path: s2pResultsPath to the directory containing suite2p outputs
        :param neuropil_coeff: choose to subtract neuropil or not when loading s2p traces
        :param save: choose to save cellsdata object or not
        """

        print(f'\----- Adding Suite2p results to .Suite2p module ...')

        s2p_path = self.s2pResultsPath if s2p_path is None else s2p_path

        # extract suite2p stat.npy and neuropil-subtracted F
        subtract_neuropil = True if neuropil_coeff > 0 else False
        FminusFneu, spks, stat, neuropil = s2p_loader(self.s2pResultsPath, subtract_neuropil=subtract_neuropil,
                                                      neuropil_coeff=neuropil_coeff)
        self.raw.append(FminusFneu)  # raw F of each suite2p ROI (neuropil corrected if neuropil_coeff > 0)
        self.spks.append(spks)  # deconvolved spikes each suite2p ROI
        self.neuropil.append(neuropil)  # neuropil value of each suite2p ROI
        self.stat.append(stat)  # stat dictionary for each suite2p ROI
        self.n_frames = spks.shape[1]

        self.ops: dict = np.load(os.path.join(s2p_path, 'ops.npy'), allow_pickle=True).item()
        self.mean_img.append(self.ops['meanImg'])
        self.mean_imgE.append(self.ops['meanImgE'])  # enhanced mean image from suite2p
        self.xoff = self.ops['xoff']  # motion correction info
        self.yoff = self.ops['yoff']

        cell_id = []
        cell_med = []
        cell_x = []
        cell_y = []

        # stat is an np array of dictionaries
        for cell, s in enumerate(stat):
            cell_id.append(cell)
            cell_med.append(s['med'])

            cell_x.append(s['xpix'])
            cell_y.append(s['ypix'])

        self.cell_id.append(cell_id)
        self.n_units.append(len(self.cell_id))
        self.cell_med.append(cell_med)
        self.cell_x.append(cell_x)
        self.cell_y.append(cell_y)

        num_units = FminusFneu.shape[0]

        print(
            f'|- Loaded {self.n_units} suite2p classified cells from plane 0, recorded for {round(self.raw[0].shape[1] / self.ops["fs"], 2)} secs total, {self.n_frames} frames total')

        # TODO consider replacing this and use returning properties
        # print(f'*** plane 0 cellsdata ***')
        self.raw = self.raw[0]
        self.spks = self.spks[0]
        self.neuropil = self.neuropil[0]
        self.mean_img = self.mean_img[0]
        self.mean_imgE = self.mean_imgE[0]
        self.xoff = self.xoff[0]
        self.yoff = self.yoff[0]
        self.cell_id = self.cell_id[0]
        self.n_units = self.n_units[0]
        self.cell_med = self.cell_med[0]
        self.cell_x = self.cell_x[0]
        self.cell_y = self.cell_y[0]
        self.n_frames = self.spks.shape[1]
        # print(f'n_frames', self.n_frames)

        # read in other files
        self.output_ops: dict = np.load(Path(self.s2pResultsPath).joinpath('ops.npy'), allow_pickle=True).item()
        stats_file = Path(self.s2pResultsPath).joinpath('stat.npy')
        self.iscell = np.load(Path(self.s2pResultsPath).joinpath('iscell.npy'), allow_pickle=True)[:, 0].astype(bool)
        self.stat = np.load(stats_file, allow_pickle=True)

        self.im = stats_dicts_to_3d_array_(stat_dict=self.stat, output_ops=self.output_ops)
        self.im[self.im == 0] = np.nan

    def stitch_reg_tiffs(self, first_frame: int, last_frame: int, reg_tif_folder: str = None,
                         s2p_run_batch: int = 2000):
        """
        Stitches together registered tiffs outputs from suite2p from the provided imaging frame start and end values.

        :param first_frame: first frame from the overall s2p run to start stitching from
        :param last_frame: last frame from the overall s2p run to end stitching at
        :param s2p_run_batch: batch size for suite2p run (defaults to 2000 because that is usually the batch size while running suite2p processing)
        """

        if reg_tif_folder is None:
            if self.s2pResultsPath:
                reg_tif_folder = self.s2pResultsPath + '/reg_tif/'
                print(f"|____ trying to load registered tiffs from: {reg_tif_folder}")
        else:
            raise Exception(f"Must provide reg_tif_folder s2pResultsPath for loading registered tiffs")
        if not os.path.exists(reg_tif_folder):
            raise Exception(f"no registered tiffs found at s2pResultsPath: {reg_tif_folder}")

        frame_range = [first_frame, last_frame]

        start = first_frame // s2p_run_batch
        end = last_frame // s2p_run_batch + 1

        # set tiff paths to save registered tiffs:
        tif_path_save = self.s2pResultsPath + 'reg_tiff.tif'
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        sorted_paths = [reg_tif_folder + tif for tif in reg_tif_list][start:end + 1]

        print(tif_path_save)
        print(sorted_paths)

        if os.path.exists(tif_path_save):
            make_tiff_stack(sorted_paths, save_as=tif_path_save)

        self.reg_tiff_path = tif_path_save

    def s2pMeanImage(self, plot: bool = True):
        """
        Return array of the s2p mean image.
        :param s2p_path: (optional) s2pResultsPath to location of s2p cellsdata
        :param plot: (optional) option to plot the s2p mean image
        :return:
        """

        print(f'Plotting s2p mean image .. {self.s2pResultsPath}')

        os.chdir(self.s2pResultsPath)

        ops = np.load('ops.npy', allow_pickle=True).item()

        mean_img = ops['meanImg']

        mean_img = np.array(mean_img, dtype='uint16')

        if plot:
            plt.imshow(mean_img, cmap='gray')
            plt.suptitle('s2p mean image')
            plt.show()

        return mean_img

    @classmethod
    def subSuite2p(cls, trialsTiffs: dict, s2pResultsPath):
        """
        Alternative constructor for Suite2pExperiment class. Inputs a specific # of trials to use for loading up Suite2p.

        :param trialsTiffs: dictionary of trial ID and associated tiff path to use for creating new Suite2p instance.
        :param s2pResultsPath: path to save suite2p results outputs.
        :return: Suite2p instance
        """
        return cls(trialsTiffsSuite2p=trialsTiffs, s2pResultsPath=s2pResultsPath)

    def add_bad_frames(self, frames, bad_frames_npy_loc) -> None:
        """
        Add frames to bad_frames.npy file that will be excluded by Suite2p during ROI segmentation.

        :param frames: frames to add as bad_frames.
        :param bad_frames_npy_loc: directory to save bad_frame.npy file
        """

        self.bad_frames.extend(frames)
        self.bad_frames = list(np.unique(self.bad_frames))

        if len(self.bad_frames) > 0:
            print('\- Appending a total of ', len(frames), 'to bad_frames.npy',
                  f"\n\t total bad_frames: {len(self.bad_frames)}")
            np.save('%s/bad_frames.npy' % bad_frames_npy_loc,
                    self.bad_frames)  # save to npy file and remember to move npy file to tiff folder before running with suite2p

    def setupForSuite2p(self, trialsSuite2P: list, TrialsInformation: dict, **kwargs):
        """
        Run housekeeping tasks to organize fields prior to running Suite2p pipeline.
        
        :param trialsSuite2P: list of trials to use for Suite2p setup
        :param TrialsInformation: dictionary mapping trial ID to tiff_path to use for Suite2p processing of that trial
        :param kwargs:
            :fs:
            :nplanes:
            :pix_sz_x:
            :pix_sz_y:
            :frame_x:
            :frame_y:
            :n_channels:
        """
        if trialsSuite2P != self.trials:

            trialsTiffs = {}
            for trial in trialsSuite2P:
                trialsTiffs[trial] = TrialsInformation[trial]['tiff_path']

            Suite2p_obj = self.subSuite2p(trialsTiffs=trialsTiffs, s2pResultsPath=self.s2pResultsPath)
        else:
            Suite2p_obj = self
        # self.trials = trialsSuite2P if trialsSuite2P else self.trials

        for trial in Suite2p_obj.trials:
            Suite2p_obj.tiff_paths_to_use_s2p[trial] = TrialsInformation[trial]['tiff_path']

        # load the first tiff in Suite2p_obj.trials to collect default metainformation about imaging setup parameters
        trial = Suite2p_obj.trials[0]
        from imagingplus import import_obj
        from imagingplus.main.core import ImagingTrial
        trialobj: ImagingTrial = import_obj(TrialsInformation[trial]['analysis_object_information']['pkl path'])

        # set imaging parameters using defaults or kwargs if provided
        fps = trialobj.imparams.fps if 'fs' not in kwargs else kwargs['fs']
        n_planes = trialobj.imparams.n_planes if 'n_planes' not in kwargs else kwargs['n_planes']
        pix_sz_x = trialobj.imparams.pix_sz_x if 'pix_sz_x' not in kwargs else kwargs['pix_sz_x']
        pix_sz_y = trialobj.imparams.pix_sz_y if 'pix_sz_y' not in kwargs else kwargs['pix_sz_y']
        frame_x = trialobj.imparams.frame_x if 'frame_x' not in kwargs else kwargs['frame_x']
        frame_y = trialobj.imparams.frame_y if 'frame_y' not in kwargs else kwargs['frame_y']
        n_channels = kwargs['n_channels'] if 'n_channels' in [*kwargs] else 1  # default is 1 channel imaging in .tiffs for suite2p

        # setup ops dictionary
        Suite2p_obj.ops['fs'] = fps / n_planes
        diameter_x = 13 / pix_sz_x
        diameter_y = 13 / pix_sz_y
        Suite2p_obj.ops['diameter'] = int(diameter_x), int(diameter_y) if diameter_y != diameter_x else diameter_x

        # set other ops parameters if provided in kwargs:
        for key in kwargs:
            if key in [*Suite2p_obj.ops]:
                Suite2p_obj.ops[key] = kwargs[key]

        # calculate batch size to use in Suite2p run
        batch_size = Suite2p_obj.ops['batch_size'] * (262144 / (
                frame_x * frame_y))  # larger key_frames will be more RAM intensive, scale user batch size based on num pixels in 512x512 images

        # setup db dict
        Suite2p_obj.db = {'fs': float(Suite2p_obj.ops['fs']), 'diameter': Suite2p_obj.ops['diameter'],
                          'batch_size': int(batch_size),
                          'nimg_init': int(Suite2p_obj.ops['batch_size']), 'nplanes': n_planes, 'nchannels': n_channels,
                          'tiff_list': list(Suite2p_obj.tiff_paths_to_use_s2p.values()),
                          'data_path': os.path.dirname(trialobj.dataPath),
                          'save_folder': Suite2p_obj.s2pResultsPath}

    def s2pRun(self, expobj, trialsSuite2P: Union[list, str] = 'all'):
        """
        Run suite2p pipeline.
        
        :param expobj: experiment object to run Suite2p on.
        :param trialsSuite2P: list of trials to run for Suite2p
        """
        s2pRun(expobj, trialsSuite2P=trialsSuite2P)

    def s2pROIsTiff(self, save_path, cell_ids: Union[str, list] = 'all'):
        """
        Save a TIFF image of the suite2p ROIs masks.
        # todo test function

        :param save_path: path to save to
        :param cell_ids: IDs of Suite2p ROI masks to show in image (use 'all' (default) to plot all ROI masks)
        """

        mask_img = np.zeros((self.output_ops['Ly'], self.output_ops['Lx']), dtype='uint8')

        cell_ids = self.cell_id if cell_ids == 'all' else cell_ids
        for n in range(0, len(self.iscell)):
            if n in cell_ids:
                ypix = self.stat[n]['ypix']
                xpix = self.stat[n]['xpix']
                mask_img[ypix, xpix] = np.random.randint(10, 255)

        WriteTiff(save_path=save_path, stack=mask_img)


# noinspection DuplicatedCode
class Suite2pResultsTrial(CellAnnotations, ImagingData):
    """Class to collect and store suite2p processed cellsdata for one trial - out of overall experiment."""

    def __init__(self, s2pExp: Suite2pExperiment, trial_frames: tuple):
        """
        Connecting Suite2p results from a specific trial (which spans trial_frames out of the overall Suite2p run) to that trial.

        :param s2pExp:
        :param trial_frames:


        """
        self.__s2pResultExists = False
        self.s2pResultsPath = None  #: path where s2p results were loaded from
        # self.raw: np.ndarray  #: array of raw Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
        # self.spks: np.ndarray  #: array of deconvolved spks values from suite2p output [num cells x length of imaging acquisition], one per plane
        # self.neuropil: np.ndarray  #: array of neuropil Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
        self.im: np.ndarray
        self.iscell = None
        self.n_units = None
        self.output_ops = None
        self.stat = None
        # self.cell_id = None

        print(f"\n\----- ADDING .Suite2p module to trial ... ", end='\r')

        self.trial_frames = trial_frames  # tuple of first and last frame (out of the overall suite2p run) corresponding to the present trial

        print(
            f"\t|- current trial key_frames: {trial_frames} out of {s2pExp.n_frames} total key_frames processed through suite2p")
        if s2pExp.s2pResultsPath:
            raw, spks, neuropil = self._get_suite2pResults(s2pExp)
        else:
            raise ValueError('cannot create s2p results trial without existing results in the input s2pExperiment.')

        cells_data, cells_multidim = self.getCellsAnnotations()
        # super(Suite2pResultsTrial, self).__init__(cells_array=cells_data.index, annotations=cells_data.columns, cellsdata=cells_data, multidim_data=cells_multidim)

        CellAnnotations.__init__(self, cells_array=cells_data.index[self.iscell], annotations=cells_data.columns,
                                 cellsdata=cells_data[self.iscell], multidim_data=cells_multidim)

        ImagingData.__init__(self, imdata=raw, spks=spks, neuropil=neuropil, data_label='')

        print(f"\n\----- ADDED .Suite2p module to trial. ", end='\r')

    def __repr__(self):
        if self.s2pResultExists:
            return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} key_frames x {self.n_rois} s2p ROIs'
        else:
            return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} key_frames. No Suite2p Results loaded.'

    def __str__(self):
        if self.s2pResultExists:
            return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} key_frames x {self.n_rois} s2p ROIs'
        else:
            return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} key_frames. No Suite2p Results loaded.'

    @property
    def s2pResultExists(self):
        """
        Return true if suite2p results exist for current trial.
        """
        return self.__s2pResultExists

    @s2pResultExists.setter
    def s2pResultExists(self, val):
        """
        Set true if suite2p results exist for current trial.
        :param val:
        """
        self.__s2pResultExists = val

    def _get_suite2pResults(self, s2pExp: Suite2pExperiment):
        """
        Get suite2p results for current trial's key_frames.
        
        :param s2pExp:
        :return:
        """

        # self.cell_id = s2pExp.cell_id
        self.stat = s2pExp.stat
        self.output_ops = s2pExp.output_ops
        self.n_units = s2pExp.n_units
        self.iscell = s2pExp.iscell

        raw = s2pExp.raw[:, self.trial_frames[0]:self.trial_frames[
            1]]
        spks = s2pExp.spks[:, self.trial_frames[0]:self.trial_frames[
            1]]
        neuropil = s2pExp.neuropil[:, self.trial_frames[0]:self.trial_frames[
            1]]

        self.s2pResultExists = True

        self.im = stats_dicts_to_3d_array_(stat_dict=s2pExp.stat, output_ops=s2pExp.output_ops)
        self.im[self.im == 0] = np.nan
        self.s2pResultsPath = s2pExp.s2pResultsPath

        return np.asarray(raw), np.asarray(spks), np.asarray(neuropil)

    def makeFrameAverageTiff(self, reg_tif_dir: str, frames: Union[int, list, tuple], peri_frames: int = 100,
                             save_dir: str = None, to_plot=False, **kwargs):
        """Creates, plots and/or saves an average image of the specified number of peri-key_frames around the given frame from the suite2p registered TIFF file.

        :param reg_tif_dir: registered tiff directory (from Suite2p results directory)
        :param frames: key frames to plot average image
        :param peri_frames: number of frames to collect pre- and post- from key frame
        :param save_dir: directory to save tiff images to
        :param to_plot: if true, show peri-frame average image
        :param kwargs: see kwargs under imagingplus.plotting.plotting.plotImg
        :return: peri-frame averaged array images
        """

        # read in registered tiff
        reg_tif_folder = reg_tif_dir
        assert os.path.exists(reg_tif_dir), '`reg_tif_dir` could not be found.'
        reg_tif_list = os.listdir(reg_tif_dir)
        reg_tif_list.sort()

        if type(frames) == int:
            frames = [frames]

        frames_s2prun_num = [(frame + self.trial_frames[0]) for frame in frames]

        imgs = []
        for idx, frame in enumerate(frames_s2prun_num):
            batch_num = frame // self.output_ops['batch_size']

            tif_path = reg_tif_folder + reg_tif_list[batch_num]

            # im_batch_reg = tf.imread(tif_path, key=range(0, self.output_ops['batch_size']))
            im_batch_reg = ImportTiff(tif_path, frames=(0, self.output_ops['batch_size']))

            frame_num_batch = frame - (batch_num * self.output_ops['batch_size'])

            if frame_num_batch < peri_frames // 2:
                peri_frames_low = frame_num_batch
            else:
                peri_frames_low = peri_frames // 2
            peri_frames_high = peri_frames // 2
            im_sub_reg = im_batch_reg[frame_num_batch - peri_frames_low: frame_num_batch + peri_frames_high]

            avg_sub = np.mean(im_sub_reg, axis=0)

            # convert to 8-bit
            from imagingplus.utils.images import convert_to_8bit
            avg_sub = convert_to_8bit(avg_sub, 0, 255)

            if save_dir:
                if '.tif' in save_dir:
                    save_dir = os.path.dirname(save_dir) + '/'
                save_path = save_dir + f'/{frames[idx]}_s2preg_frame_avg.tif'
                os.makedirs(save_dir, exist_ok=True)

                print(f"\t\- Saving averaged s2p registered tiff for frame: {frames[idx]}, to: {save_path}")
                tf.imwrite(save_path, avg_sub, photometric='minisblack')

            if to_plot:  # todo replace with proper plotting average tiff frame code
                kwargs['title'] = f'{peri_frames} key_frames avg from s2p reg tif, frame: {frames[idx]}' if not kwargs['title'] else kwargs['title']
                from imagingplus.plotting.plotting import plotImg
                plotImg(img=avg_sub, **kwargs)

            imgs.append(avg_sub)

        return np.asarray(imgs)

    def makeDownSampledTiff(self, save_path: str, group_by: int = 4, reg_tif_folder=None, frameNum: Union[str, tuple, list, int]='all'):
        """
        Makes downsampled TIFF of current suite2p trial imaging series
        
        :param save_path: path to save downsampled TIFF.
        :param group_by: number of frames to created grouped average output.
        :param reg_tif_folder: registered tiff directory
        :param frameNum: frame numbers (start and end) to create downsampled tiff for (use 'all' (default) to create downsampled tiff for all frames).
        """
        reg_tif_folder = reg_tif_folder if reg_tif_folder else self.s2pResultsPath + '/reg_tif/'

        sorted_paths, first_tiff_offset, last_tiff_offset = self.getRegTiffPaths(reg_tif_folder=reg_tif_folder, frameNum=frameNum)
        data = make_tiff_stack(sorted_paths, save_path=None)
        trial_frames_cropped = data[first_tiff_offset: -last_tiff_offset]
        SaveDownsampledTiff(stack=trial_frames_cropped, group_by=group_by, save_as=save_path)

    def getRegTiffPaths(self, reg_tif_folder=None, frameNum: Union[str, int, tuple, list] = 'all'):
        """
        Return current imaging trial's suite2p output registered tiff paths for a certain frame number (default is for all trial frames). note that suite2p uses batched registration and creates batched motion registered tiffs.

        :param reg_tif_folder: registered tiff directory
        :param frameNum: frame numbers (start and end) to create downsampled tiff for (use 'all' (default) to create downsampled tiff for all frames).
        """
        reg_tif_folder = reg_tif_folder if reg_tif_folder else self.s2pResultsPath + '/reg_tif/'
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        tif_list = [file for file in reg_tif_list if '.tif' in file]

        curr_trial_frames = self.trial_frames
        batch_size = self.output_ops['batch_size']

        if type(frameNum) is int:
            s2p_frame = curr_trial_frames[0] + frameNum

            idx = s2p_frame // batch_size
            tiff_path = reg_tif_folder + tif_list[idx]

            offsetFrameNum = int(s2p_frame % batch_size)

            return tiff_path, offsetFrameNum
        elif type(frameNum) is list or type(frameNum) is tuple:
            start = int(np.floor((curr_trial_frames[0] + frameNum[0]) / batch_size))
            end = int(np.ceil((curr_trial_frames[0] + frameNum[1]) / batch_size))

            if start == end:  # i.e. if the frames selected for crop are in the same batch .tif
                stacks = [start]
                last_tiff_offset = batch_size - frameNum[1]
            else:
                stacks = range(int(start), int(end)) if end - start > 1 else [start, end]
                last_tiff_offset = batch_size - int((curr_trial_frames[0] + frameNum[1]) % batch_size)

            reg_tif_list = os.listdir(reg_tif_folder)
            reg_tif_list.sort()
            tif_list = [file for file in reg_tif_list if '.tif' in file]

            tiff_paths_list = [reg_tif_folder + tif_list[i] for i in stacks]

            first_tiff_offset = int(curr_trial_frames[0] + frameNum[0] - (start * batch_size))

            print(tiff_paths_list, first_tiff_offset, last_tiff_offset)
            return tiff_paths_list, first_tiff_offset, last_tiff_offset
        elif frameNum == 'all':
            reg_tif_list = os.listdir(reg_tif_folder)
            reg_tif_list.sort()
            tif_list = [file for file in reg_tif_list if '.tif' in file]
            start = curr_trial_frames[0] // batch_size
            end = curr_trial_frames[1] // batch_size

            tiff_paths_list = [reg_tif_folder + tif_list[i] for i in range(start, end)]

            first_tiff_offset = int(curr_trial_frames[0] - (start * batch_size))
            last_tiff_offset = int(curr_trial_frames[1] % batch_size)
            return tiff_paths_list, first_tiff_offset, last_tiff_offset

        else:
            raise ValueError(f'frameNum must be `all` or an int of a frame number.')

    def getCellsAnnotations(self):
        """
        Get cell annotations from Suite2p output to use in creating anndata structure.
        """
        if self.s2pResultExists:
            # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
            # build dataframe for obs_meta from suite2p stat information
            obs_meta = pd.DataFrame(columns=[*self.stat[0]], index=range(len(self.stat)))[self.iscell]

            for idx, __stat in enumerate(self.stat):
                for __column in [*__stat]:
                    if __column in obs_meta.columns:
                        obs_meta.loc[idx, __column] = __stat[__column]

            obs_m = {'ypix': [],
                     'xpix': []}
            for col in [*obs_m]:
                for idx, __stat in enumerate(self.stat):
                    obs_m[col].append(__stat[col])
                obs_m[col] = np.asarray(obs_m[col])[self.iscell]

            return obs_meta, obs_m

        else:
            raise ValueError('cannot get cell annotations. no s2p results found in trial.')

    def stitch_s2p_reg_tiff(self, save_path: str = None, reg_tif_folder=None):
        """
        Concatenate suite2p registered tiffs into a single stack.

        :param save_path: path to save stack.
        :param reg_tif_folder: folder to retrieve suite2p registered tiffs from.
        """
        # TODO test out

        assert self.s2pResultExists, UnavailableOptionError('stitch_s2p_reg_tiff')

        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        tif_path_save2 = save_path if save_path else None

        start = self.trial_frames[0] // self.output_ops[
            'batch_size']  # 2000 because that is the batch size for suite2p run
        end = self.trial_frames[1] // self.output_ops['batch_size'] + 1

        if reg_tif_folder is None:
            if self.s2pResultsPath:
                reg_tif_folder = self.s2pResultsPath + '/reg_tif/'
                print(f"|____ trying to load registered tiffs from: {reg_tif_folder}")
        else:
            raise Exception(f"Must provide reg_tif_folder s2pResultsPath for loading registered tiffs")
        if not os.path.exists(reg_tif_folder):
            raise Exception(f"no registered tiffs found at s2pResultsPath: {reg_tif_folder}")

        # set tiff paths to save registered tiffs:
        tif_path_save = tif_path_save2
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        sorted_paths = [reg_tif_folder + tif for tif in reg_tif_list][start:end + 1]

        print(tif_path_save)
        print(sorted_paths)

        if os.path.exists(tif_path_save):
            data = make_tiff_stack(sorted_paths, save_path=tif_path_save)

        if not os.path.exists(tif_path_save2):
            with tf.TiffWriter(tif_path_save2, bigtiff=True) as tif:
                with tf.TiffFile(tif_path_save2, multifile=False) as input_tif:
                    print('cropping registered tiff')
                    data = input_tif.asarray()
                    print('shape of stitched tiff: ', data.shape)

                # todo test and debug and fix
                reg_tif_crop = data[self.trial_frames[0] - start * self.output_ops['batch_size']:
                                    self.trial_frames[1] - (
                                            self.trial_frames - start * self.output_ops['batch_size'])]
                print('saving cropped tiff ', reg_tif_crop.shape)
                tif.write(reg_tif_crop)

    @property
    def cell_coords(self):
        """Return X and Y coordinates of cells
        """
        assert 'med' in self.cellsdata, 'med cannot be found in cells annotations under cellsdata'
        coordinates = np.empty(shape=[self.n_cells, 2])
        for cell, value in enumerate(self.cellsdata['med']):
            coordinates[cell, 0] = value[1]
            coordinates[cell, 1] = value[0]
        return coordinates

    def collectTiffList(self, curr_trial_frames: Union[tuple, list], reg_tif_folder: str = None):
        """
        Collects list of tiffs from registered tiff output folder for current trial of suite2p processed dataset.

        :param curr_trial_frames: frame start and end for current trial (relative to total suite2p experiment results).
        :param reg_tif_folder: folder to retrieve suite2p registered tiffs from.
        :return: list of tiffs
        """

        if not os.path.exists(reg_tif_folder):
            raise Exception(f"no registered tiffs found at path: {reg_tif_folder}")

        # collect list of registered tiffs associated with current trial
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        start = curr_trial_frames[0] // self.output_ops['batch_size']
        end = curr_trial_frames[1] // self.output_ops['batch_size'] + 1

        tiffs_list = []
        for i in range(start, end):
            tiff_path = reg_tif_folder + reg_tif_list[i]
            tiffs_list.append(tiff_path)

        return tiffs_list

    def collectSignalFromCoords(self, target_coords_masks: np.ndarray, reg_tif_folder: str = None) -> np.ndarray:
        """
        Collect average fluorescence signal from input coordinates (and their masks).

        uses registered tiffs to collect raw traces from provided target masks target areas
        :param target_coords_masks: array of masks of coordinates to collect signal from (for each cell).
        :param reg_tif_folder: folder to retrieve suite2p registered tiffs from.
        """

        reg_tif_folder = reg_tif_folder if reg_tif_folder else self.s2pResultsPath
        batch_size = self.output_ops['batch_size']

        if not os.path.exists(reg_tif_folder):
            raise Exception(f"no registered tiffs found at path: {reg_tif_folder}")

        num_coords = len(target_coords_masks)
        target_areas = target_coords_masks

        print(
            f'\n\ncollecting raw Flu traces from SLM target coord. areas from registered TIFFs from: {reg_tif_folder}')
        # read in registered tiff
        reg_tif_list = os.listdir(reg_tif_folder)
        reg_tif_list.sort()
        tif_list = [file for file in reg_tif_list if '.tif' in file]
        start = self.trial_frames[0] // batch_size
        end = self.trial_frames[1] // batch_size

        # collect mean traces from target areas of each target coordinate by reading in individual registered tiffs that contain frames for current trial
        # targets_trace_full = np.zeros([num_coords, (end - start) * batch_size], dtype='float32')
        targets_trace_full = np.zeros([num_coords, 0], dtype='float32')

        for i in range(start, end):
            try:
                tiff_path = reg_tif_folder + tif_list[i]
                print(f'|- reading tiff: {tiff_path}')
            except IndexError:
                break

            import cv2
            ret, images = cv2.imreadmulti(tiff_path, [], cv2.IMREAD_ANYCOLOR)
            data = np.asarray(images)

            print(f'\t\- shape: {data.shape}')
            targets_trace = np.zeros([num_coords, data.shape[0]], dtype='float32')
            for coord in range(num_coords):
                x = data[:, target_areas[coord, 1], target_areas[coord, 0]]
                targets_trace[coord] = np.mean(x, axis=1)

            targets_trace_full = np.concatenate((targets_trace_full, targets_trace), axis=1)

            # targets_trace_full[:, (i - start) * batch_size: ((i - start) * batch_size) + data.shape[
            #     0]] = targets_trace  # iteratively write to each successive segment of the targets_trace array based on the length of the reg_tiff that is read in.

        # final part, crop to the *exact* frames for current trial
        start_crop = self.trial_frames[0] - (start * batch_size)
        end_crop = self.trial_frames[0] - (start * batch_size) + (self.trial_frames[1] - self.trial_frames[0])
        raw_coords_signal = targets_trace_full[:, start_crop: end_crop]

        return raw_coords_signal

#### archiving away for now - trying to switch to an approach that doesn't inherit from parent suite2p obj.
# class Suite2PTrial_(Suite2pExperiment):
#     """used to collect and store suite2p processed cellsdata for one trial - out of overall experiment."""
#
#     def __init__(self, trialsTiffsSuite2p: dict, trial_frames: tuple, s2pResultsPath: Optional[str] = None):
#         """
#         Connecting Suite2p results from a specific trial (which spans trial_frames out of the overall Suite2p run) to that trial.
#
#         :param trialsTiffsSuite2p:
#         :param dataPath:
#         :param trial_frames:
#         :param s2pResultsPath:
#         :param subtract_neuropil:
#         """
#         # todo need to reconsider this heirarchy; it saves in the entire suite2p result for every trial....
#         super().__init__(trialsTiffsSuite2p,
#                          s2pResultsPath)  # - TODO it really is confusing to be passing in all trials for a s2p results obj that should be restricted to just one trial
#
#         print(f"\n\----- ADDING .Suite2p module to trial ... ", end='\r')
#
#         ## initializing attributes to collect Suite2p results for this specific trial
#         # self.dfof: list(np.ndarray) = []  # array of dFF normalized Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
#         # self.__raw: List[np.ndarray] = self.raw
#         # self.__spks: List[np.ndarray] = self.spks
#         # self.__neuropil: List[np.ndarray] = self.neuropil
#
#         self.trial_frames = trial_frames  # tuple of first and last frame (out of the overall suite2p run) corresponding to the present trial
#
#         # self.suite2p_overall = suite2p_experiment_obj
#         print(
#             f"\t|- current trial key_frames: {trial_frames} out of {self.n_frames} total key_frames processed through suite2p")
#         self._get_suite2pResults() if self.s2pResultsPath else None
#
#         # s2pResultsPath = suite2p_experiment_obj.s2pResultsPath if hasattr(suite2p_experiment_obj, 's2pResultsPath') else None
#         # Suite2pExperiment.__init__(self, trialsSuite2p = suite2p_experiment_obj.trials, s2pResultsPath=s2pResultsPath,
#         #                                   subtract_neuropil=suite2p_experiment_obj.subtract_neuropil)
#
#         print(f"\n\----- ADDED .Suite2p module to trial. ", end='\r')
#
#     def __repr__(self):
#         if self.s2pResultExists:
#             return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} key_frames x {self.n_units} s2p ROIs'
#         else:
#             return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} key_frames. No Suite2p Results loaded.'
#
#     def _get_suite2pResults(self):  # TODO complete code for getting suite2p results for trial
#         """crop suite2p cellsdata for key_frames only for the present trial"""
#
#         # for attr in ['n_units', 'cell_id', 'cell_plane', 'cell_x', 'cell_y', 'xoff', 'yoff', 'raw', 'spks', 'neuropil', 'stat']:
#         #     try:
#         #         setattr(self, attr, getattr(self.suite2p_overall, attr))
#         #     except AttributeError:
#         #         pass
#
#         if self.n_planes == 1:
#             self.cell_id = self.cell_id
#             self.stat = self.stat[0]
#
#             self.raw = self.raw[:, self.trial_frames[0]:self.trial_frames[
#                 1]]  # array of raw Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
#             self.spks = self.spks[:, self.trial_frames[0]:self.trial_frames[
#                 1]]  # array of deconvolved spks values from suite2p output [num cells x length of imaging acquisition], one per plane
#             self.neuropil = self.neuropil[:, self.trial_frames[0]:self.trial_frames[
#                 1]]  # array of neuropil Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
#
#             self.s2pResultExists = True
#
#         else:
#             for plane in range(self.n_planes):
#                 self.raw.append(self.raw[plane][:, self.trial_frames[0]:self.trial_frames[1]])
#                 self.spks.append(self.spks[plane][:, self.trial_frames[0]:self.trial_frames[1]])
#                 self.neuropil.append(self.neuropil[plane][:, self.trial_frames[0]:self.trial_frames[1]])
#
#                 # self.dfof.append(normalize_dff(self.raw[plane]))  # calculate df/f based on relevant key_frames
#                 self.s2pResultExists = True
#
#     def makeFrameAverageTiff(self, reg_tif_dir: str, frames: Union[int, list], peri_frames: int = 100,
#                              save_dir: str = None,
#                              to_plot=False):
#         """Creates, plots and/or saves an average image of the specified number of peri-key_frames around the given frame from the suite2p registered TIFF file.
#         """
#
#         # read in registered tiff
#         reg_tif_folder = reg_tif_dir
#         assert os.path.exists(reg_tif_dir), '`reg_tif_dir` could not be found.'
#         reg_tif_list = os.listdir(reg_tif_dir)
#         reg_tif_list.sort()
#
#         if type(frames) == int:
#             frames = [frames]
#
#         frames_s2prun_num = [(frame + self.trial_frames[0]) for frame in frames]
#
#         for idx, frame in enumerate(frames_s2prun_num):
#             batch_num = frame // self.ops['batch_size']
#
#             tif_path = reg_tif_folder + reg_tif_list[batch_num]
#
#             im_batch_reg = tf.imread(tif_path, key=range(0, self.ops['batch_size']))
#
#             frame_num_batch = frame - (batch_num * self.ops['batch_size'])
#
#             if frame_num_batch < peri_frames // 2:
#                 peri_frames_low = frame_num_batch
#             else:
#                 peri_frames_low = peri_frames // 2
#             peri_frames_high = peri_frames // 2
#             im_sub_reg = im_batch_reg[frame_num_batch - peri_frames_low: frame_num_batch + peri_frames_high]
#
#             avg_sub = np.mean(im_sub_reg, axis=0)
#
#             # convert to 8-bit
#             from imagingplus.utils.utils import convert_to_8bit
#             avg_sub = convert_to_8bit(avg_sub, 0, 255)
#
#             if save_dir:
#                 if '.tif' in save_dir:
#                     from imagingplus.utils.utils import return_parent_dir
#                     save_dir = return_parent_dir(save_dir) + '/'
#                 save_path = save_dir + f'/{frames[idx]}_s2preg_frame_avg.tif'
#                 os.makedirs(save_dir, exist_ok=True)
#
#                 print(f"\t\- Saving averaged s2p registered tiff for frame: {frames[idx]}, to: {save_path}")
#                 tf.imwrite(save_path, avg_sub, photometric='minisblack')
#
#             if to_plot:
#                 plt.imshow(avg_sub, cmap='gray')
#                 plt.suptitle(f'{peri_frames} key_frames avg from s2p reg tif, frame: {frames[idx]}')
#                 plt.show()  # just plot for now to make sure that you are doing things correctly so far


# noinspection DuplicatedCode
def add_suite2p_results(expobj, s2p_trials: Union[list, str] = 'all', s2pResultsPath: Optional[str] = None):
    """
    todo: probably want to split up the bottom
    todo: test moving function to suite2p file.
    Wrapper for adding suite2p results to Experiment. Can only be run after adding trials to Experiment.

    :param expobj: Experiment object to add suite2p results to
    :param s2p_trials: list of trials to use for Suite2p processing/analysis pipeline, default = 'all' to use all trials of Experiment
    :param s2pResultsPath: optional, if suite2p already run then here provide path existing suite2p results
    """

    print(f'\- Adding suite2p module to experiment. Located under .Suite2p')

    if not len([*expobj.TrialsInformation]) > 0:
        from imagingplus.utils.classes import UnavailableOptionError
        raise UnavailableOptionError('need to add at least 1 trial to Experiment before adding Suite2p functionality.')

    if s2p_trials == 'all': s2p_trials = expobj.trialIDs
    assert type(s2p_trials) == list and len(s2p_trials) > 0, 'no s2p trials given in list form.'
    _trialsTiffsSuite2p = {}
    for trial in s2p_trials: _trialsTiffsSuite2p[trial] = expobj.TrialsInformation[trial]['paths']['dataPath']

    if s2pResultsPath:  # if s2pResultsPath provided then try to find and pre-load results from provided s2pResultsPath, raise error if cannot find results
        # search for suite2p results items in expobj.suite2pPath, and auto-assign s2pRunComplete --> True if found successfully
        __suite2p_path_files = os.listdir(s2pResultsPath)
        expobj.Suite2p.s2pResultExists = False
        for filepath in __suite2p_path_files:
            if 'ops.npy' in filepath:
                expobj.Suite2p.s2pResultExists = True
                break
        if expobj.Suite2p.s2pResultExists:
            expobj._suite2p_save_path = s2pResultsPath
            # todo - why is this re'initing the whole s2p module??? this action should just be to add results to pre-existing module
            expobj.Suite2p = Suite2pResultsExperiment(trialsTiffsSuite2p=expobj.Suite2p._trialsTiffsSuite2p,
                                                      s2pResultsPath=s2pResultsPath)
        else:
            raise ValueError(
                f"suite2p results could not be found. `suite2pPath` provided was: {s2pResultsPath}")
    else:  # no s2pResultsPath provided, so initialize without pre-loading any results
        expobj.Suite2p = Suite2pResultsExperiment(trialsTiffsSuite2p=expobj.Suite2p._trialsTiffsSuite2p)

    # print(expobj.Suite2p)
    # adding of suite2p trial level as well in this function as well
    total_frames = 0
    for trial in s2p_trials:
        trialobj = expobj.load_trial(trialID=trial)
        # print(f'n_frames', expobj.Suite2p.n_frames)
        trialobj.Suite2p = Suite2pResultsTrial(s2pExp=expobj.Suite2p, trial_frames=(
            total_frames, total_frames + trialobj.n_frames))  # use trial obj's current trial key_frames
        trialobj.save()
        total_frames += trialobj.n_frames

    print(f'|- Finished adding suite2p module to experiment and trials. Located under .Suite2p')
    expobj.save()
