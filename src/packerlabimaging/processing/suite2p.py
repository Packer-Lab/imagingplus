import os
from pathlib import Path
from typing import Dict, Any, Optional, List

import numpy as np
import suite2p
import tifffile as tf
import matplotlib.pyplot as plt

from packerlabimaging.ExperimentMain import TrialsInformation
from packerlabimaging.utils.utils import make_tiff_stack, s2p_loader, normalize_dff, Utils

# TEMP VARIABLES FOR DEVELOPMENT USAGES
N_PLANES = 1
NEUROPIL_COEFF = 0.7


class Suite2pResultsExperiment:
    """used to run and further process suite2p processed data, and analysis associated with suite2p processed data."""

    def __init__(self, trialsSuite2p: list, s2pResultsPath: Optional[str] = None, subtract_neuropil: bool = True):
        self.reg_tiff_path = None
        self.ops_end = None
        self.s2pResultsPath = None
        self._s2pResultExists = None
        self.db = None
        self.ops = None
        print(f"\- ADDING Suite2p Results to Experiment object ... ")

        ## initialize attr's
        # TODO add attr comments
        self.n_planes: int = N_PLANES

        self.cell_id = []  # cell ROI # id per plane
        self.n_units = []  # num of ROIs in suite2p result per plane
        self.cell_plane = []  # the corresponding imaging plane for each cell
        self.cell_med = []  # y, x location of each cell in the imaging frame
        self.cell_x = []  # TODO ROB
        self.cell_y = []  # TODO ROB
        self.raw = []  # [cells x num frames], raw traces for each cell from suite2p
        self.stat = []
        self.spks = []
        self.neuropil = []
        self.mean_img = []  # mean img output from suite2p
        self.mean_imgE = []  # mean img Enhanced output from suite2p
        self.radius = []  # cell radius from suite2p
        self.xoff = []  # motion correction info
        self.yoff = []  # motion correction info

        # set trials to run together in suite2p for Experiment
        self.trials = trialsSuite2p
        self.subtract_neuropil = subtract_neuropil
        assert len(
            self.trials) > 0, "no trials found to run suite2p, option available to provide list of trial IDs in " \
                              "`trialsSuite2P` "

        if s2pResultsPath is None:
            # initialize needed variables and attr's for future calling of s2pRun
            self.s2pResultsPath = None
        else:
            self.s2pResultsPath = s2pResultsPath
            try:
                neuropil_coeff = NEUROPIL_COEFF if subtract_neuropil else 0
                self._retrieveSuite2pData(self.s2pResultsPath, neuropil_coeff=neuropil_coeff)
            except:
                raise ValueError(
                    f'Something went wrong while trying to load suite2p processed data from: {s2pResultsPath}')

        # Attributes
        self.n_frames: int = 0  # total number of imaging frames in the Suite2p run

    def __repr__(self):
        return 'Suite2p Results (Experiment level) Object'

    def _retrieveSuite2pData(self, s2p_path: str = None, neuropil_coeff: float = 0.7):
        """processing of suite2p data from the current t-series
        :param s2p_path: s2pResultsPath to the directory containing suite2p outputs
        :param neuropil_coeff: choose to subtract neuropil or not when loading s2p traces
        :param save: choose to save data object or not
        """

        s2p_path = self.s2pResultsPath if s2p_path is None else s2p_path

        for plane in range(self.n_planes):  # TODO really don't know how planes are collected and fed into suite2p

            # extract suite2p stat.npy and neuropil-subtracted F
            substract_neuropil = True if neuropil_coeff > 0 else False
            FminusFneu, spks, stat, neuropil = s2p_loader(self.s2pResultsPath, subtract_neuropil=substract_neuropil,
                                                          neuropil_coeff=neuropil_coeff)
            self.raw.append(FminusFneu)  # raw F of each suite2p ROI
            self.spks.append(spks)  # deconvolved spikes each suite2p ROI
            self.neuropil.append(neuropil)  # neuropil value of each suite2p ROI
            self.stat.append(stat)  # stat dictionary for each suite2p ROI

            self.ops = np.load(os.path.join(s2p_path, 'ops.npy'), allow_pickle=True).item()
            self.mean_img.append(self.ops['meanImg'])
            self.mean_imgE.append(self.ops['meanImgE'])  # enhanced mean image from suite2p
            self.xoff = self.ops['xoff']  # motion correction info
            self.yoff = self.ops['yoff']

            cell_id = []
            cell_plane = []
            cell_med = []
            cell_x = []
            cell_y = []

            # stat is an np array of dictionaries
            for cell, s in enumerate(stat):
                cell_id.append(s['original_index'])
                cell_med.append(s['med'])

                cell_x.append(s['xpix'])
                cell_y.append(s['ypix'])

            self.cell_id.append(cell_id)
            self.n_units.append(len(self.cell_id[plane]))
            self.cell_med.append(cell_med)
            self.cell_x.append(cell_x)
            self.cell_y.append(cell_y)

            num_units = FminusFneu.shape[0]
            cell_plane.extend([plane] * num_units)
            self.cell_plane.append(cell_plane)

            print(
                f'|- Loaded {self.n_units} suite2p classified cells from plane {plane}, recorded for {round(self.raw[plane].shape[1] / self.ops["fs"], 2)} secs total')

        # consider replacing this and use returning properties
        if self.n_planes == 1:
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
            self.cell_plane = self.cell_plane[0]
            self.cell_x = self.cell_x[0]
            self.cell_y = self.cell_y[0]

        # read in other files
        self.output_op = np.load(Path(self.s2pResultsPath).joinpath('ops.npy'), allow_pickle=True).item()
        stats_file = Path(self.s2pResultsPath).joinpath('stat.npy')
        self.iscell = np.load(Path(self.s2pResultsPath).joinpath('iscell.npy'), allow_pickle=True)[:, 0].astype(bool)
        self.stats = np.load(stats_file, allow_pickle=True)

        # quick workaround (patch of suite2p code) because of suite2p internal error for these methods in ROI class
        def from_stat_dict(stat: Dict[str, Any]) -> suite2p.ROI:
            return suite2p.ROI(ypix=stat['ypix'], xpix=stat['xpix'], lam=stat['lam'], med=stat['med'], do_crop=False)

        def to_array_(self, Ly: int, Lx: int) -> np.ndarray:
            """Returns a 2D boolean array of shape (Ly x Lx) indicating where the roi is located."""
            arr = np.zeros((Ly, Lx), dtype=float)
            arr[self.ypix, self.xpix] = 1
            return arr

        def stats_dicts_to_3d_array_():
            arrays = []
            for i, stat in enumerate(self.stats):
                array = to_array_(from_stat_dict(stat=stat), Ly=self.output_op['Ly'], Lx=self.output_op['Lx'])
                array *= i + 1
                arrays.append(array)
            return np.stack(arrays)

        self.im = stats_dicts_to_3d_array_()
        self.im[self.im == 0] = np.nan

    # consider use of properties
    # @property
    # def raw(self):
    #     if self.n_planes == 1:
    #         return self.raw[0]
    #
    # @property
    # def spks(self):
    #     if self.n_planes == 1:
    #         return self.spks[0]
    #
    # @property
    # def neuropil(self):
    #     if self.n_planes == 1:
    #         return self.neuropil[0]

    def stitch_reg_tiffs(self, first_frame: int, last_frame: int, reg_tif_folder: str = None, force_crop: bool = False,
                         s2p_run_batch: int = 2000):  # TODO review and refactor as method for trial object
        """
        Stitches together registered tiffs outputs from suite2p from the provided imaging frame start and end values.

        :param first_frame: first frame from the overall s2p run to start stitching from
        :param last_frame: last frame from the overall s2p run to end stitching at
        :param force_crop:
        :param force_stack:
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
        :param s2p_path: (optional) s2pResultsPath to location of s2p data
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

    # suite2p methods - refactored currently to _utils.Utils !!!!!
    def s2pRun(self, ops, db):  # TODO gotta specify # of planes somewhere here
        """
        run suite2p on the experiment object, using trials specified in current experiment object, using the attributes
        determined directly from the experiment object.

        :param ops:
        :param db:
        """

        self.ops = ops
        self.db = db

        print(f'Starting Suite2p run for : \n\t{db}')

        from suite2p import run_s2p
        self.ops_end = run_s2p(ops=ops, db=db)

        self._s2pResultExists = True
        self.s2pResultsPath = self.ops_end['save_path']


class Suite2pResultsTrial(Suite2pResultsExperiment):
    """used to collect suite2p processed data for one trial - out of overall experiment."""

    def __init__(self, trialsSuite2p: list, trial_frames: tuple, s2pResultsPath: Optional[str] = None,
                 subtract_neuropil: bool = True):

        super().__init__(trialsSuite2p, s2pResultsPath, subtract_neuropil)

        print(f"\n\----- ADDING Suite2pResultsTrial ... ")

        ## initializing attributes to collect Suite2p results for this specific trial
        # self.dfof: list(np.ndarray) = []  # array of dFF normalized Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
        self.__raw: List[np.ndarray] = self.raw
        self.__spks: List[np.ndarray] = self.spks
        self.__neuropil: List[np.ndarray] = self.neuropil

        self.trial_frames = trial_frames  # tuple of first and last frame (out of the overall suite2p run) corresponding to the present trial

        # self.suite2p_overall = suite2p_experiment_obj
        print(
            f"\t|- current trial frames: {trial_frames} out of {self.n_frames} total frames processed through suite2p")
        self._get_suite2pResults() if self.s2pResultsPath else None

        # s2pResultsPath = suite2p_experiment_obj.s2pResultsPath if hasattr(suite2p_experiment_obj, 's2pResultsPath') else None
        # Suite2pResultsExperiment.__init__(self, trialsSuite2p = suite2p_experiment_obj.trials, s2pResultsPath=s2pResultsPath,
        #                                   subtract_neuropil=suite2p_experiment_obj.subtract_neuropil)

    def __repr__(self):
        return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} frames x {self.n_units} s2p ROIs'

    def _get_suite2pResults(self):  # TODO complete code for getting suite2p results for trial
        """crop suite2p data for frames only for the present trial"""

        # for attr in ['n_units', 'cell_id', 'cell_plane', 'cell_x', 'cell_y', 'xoff', 'yoff', 'raw', 'spks', 'neuropil', 'stat']:
        #     try:
        #         setattr(self, attr, getattr(self.suite2p_overall, attr))
        #     except AttributeError:
        #         pass

        if self.n_planes == 1:
            self.cell_id = self.cell_id
            self.stat = self.stat[0]

            self.raw = self.__raw[:, self.trial_frames[0]:self.trial_frames[
                1]]  # array of raw Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
            self.spks = self.__spks[:, self.trial_frames[0]:self.trial_frames[
                1]]  # array of deconvolved spks values from suite2p output [num cells x length of imaging acquisition], one per plane
            self.neuropil = self.__neuropil[:, self.trial_frames[0]:self.trial_frames[
                1]]  # array of neuropil Flu values from suite2p output [num cells x length of imaging acquisition], one per plane

            self.__s2pResultExists = True

        else:
            for plane in range(self.n_planes):
                self.raw.append(self.raw[plane][:, self.trial_frames[0]:self.trial_frames[1]])
                self.spks.append(self.spks[plane][:, self.trial_frames[0]:self.trial_frames[1]])
                self.neuropil.append(self.neuropil[plane][:, self.trial_frames[0]:self.trial_frames[1]])

                # self.dfof.append(normalize_dff(self.raw[plane]))  # calculate df/f based on relevant frames
                self.__s2pResultExists = True

    @property
    def _s2pResultExists(self):
        return self.__s2pResultExists

    # @property
    # def stat(self):
    #     if self.n_planes == 1:
    #         return self.stat[0]
    #     else:
    #         return self.stat
