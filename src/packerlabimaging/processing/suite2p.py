import os
import time
from pathlib import Path
from typing import Dict, Any, Optional, List, Union

import numpy as np
import suite2p
import matplotlib.pyplot as plt
from suite2p import run_s2p
import tifffile as tf

from packerlabimaging.utils.utils import make_tiff_stack, save_array_to_tiff

# TEMP VARIABLES FOR DEVELOPMENT USAGES
N_PLANES = 1

# suite2p methods
def s2pRun(expobj, trialsSuite2P: Union[list, str] = 'all'):  ## TODO gotta specify # of planes somewhere here
    """run suite2p for an Experiment object, using trials specified in current experiment object, using the attributes
    determined directly from the experiment object.

    :param expobj: a packerlabimaging.main.Experiment object
    :param trialsSuite2P: list of trialIDs from experiment to use in running suite2p
    """

    expobj.Suite2p.trials = trialsSuite2P if trialsSuite2P != 'all' else expobj.Suite2p.trials
    # expobj._trialsTiffsSuite2p = trialsSuite2P if trialsSuite2P else expobj._trialsTiffsSuite2p

    tiffs_paths_to_use_s2p = []
    for trial in expobj.Suite2p.trials:
        tiffs_paths_to_use_s2p.append(expobj.TrialsInformation[trial]['tiff_path'])

    expobj.Suite2p.tiffs_paths_to_use_s2p = tiffs_paths_to_use_s2p

    # # load the first tiff in expobj.Suite2p.trials to collect default metainformation about imaging setup parameters
    # trial = expobj.Suite2p.trials[0]
    # from packerlabimaging import TwoPhotonImagingTrial
    # trialobj: TwoPhotonImagingTrial = import_obj(expobj.TrialsInformation[trial]['analysis_object_information']['pkl path'])

    # # set imaging parameters using defaults or kwargs if provided
    # frame_x = trialobj.imparams.frame_x if 'frame_x' not in [*kwargs] else kwargs['frame_x']
    # frame_y = trialobj.imparams.frame_y if 'frame_y' not in [*kwargs] else kwargs['frame_y']

    # # setup ops dictionary
    # # expobj.Suite2p.ops['fs'] = fps / n_planes
    # batch_size = expobj.Suite2p.ops['batch_size'] * (262144 / (
    #         frame_x * frame_y))  # larger frames will be more RAM intensive, scale user batch size based on num pixels in 512x512 images

    # # set other ops parameters if provided in kwargs:
    # for key in [*kwargs]:
    #     if key in [*expobj.Suite2p.ops]:
    #         expobj.Suite2p.ops[key] = kwargs[key]


    # setup db dict
    expobj.Suite2p.ops['batch_size'] = int(expobj.Suite2p.ops['batch_size'])
    db = {'fs': float(expobj.Suite2p.ops['fs']), 'diameter': expobj.Suite2p.ops['diameter'], 'batch_size': expobj.Suite2p.ops['batch_size'],
          'nimg_init': expobj.Suite2p.ops['batch_size'], 'nplanes': expobj.Suite2p.ops['nplanes'], 'nchannels': expobj.Suite2p.ops['nchannels'],
          'tiff_list': list(tiffs_paths_to_use_s2p), 'data_path': expobj.dataPath, 'save_folder': expobj._suite2p_save_path
          }

    # db = {'fs': float(expobj.Suite2p.ops['fs']), 'diameter': expobj.Suite2p.ops['diameter'], 'batch_size': int(batch_size),
    #       'nimg_init': int(batch_size), 'nplanes': n_planes, 'nchannels': n_channels,
    #       'tiff_list': tiffs_paths_to_use_s2p, 'data_path': expobj.dataPath,
    #       'save_folder': expobj._suite2p_save_path}

    expobj.Suite2p.db = db
    expobj.save()

    print('RUNNING Suite2p:')
    print(f'\- db: ') #\n\t{expobj.Suite2p.db}\n\n')

    t1 = time.time()
    opsEnd = run_s2p(ops=expobj.Suite2p.ops, db=expobj.Suite2p.db)
    t2 = time.time()
    print('Total time to run suite2p was {}'.format(t2 - t1))

    # update expobj.Suite2p.ops and db
    expobj.Suite2p.ops = opsEnd

    expobj.__s2pResultExists = True
    expobj.s2pResultsPath = expobj.suite2p_save_path + '/plane0/'  ## need to further debug that the flow of the suite2p s2pResultsPath makes sense

    expobj.save()

def s2p_loader(s2p_path, subtract_neuropil=True, neuropil_coeff=0.7):
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
    for i, s in enumerate(stat):
        s['original_index'] = i

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
    return suite2p.ROI(ypix=stat['ypix'], xpix=stat['xpix'], lam=stat['lam'], med=stat['med'], do_crop=False)

def to_array_(self, Ly: int, Lx: int) -> np.ndarray:
    """Returns a 2D boolean array of shape (Ly x Lx) indicating where the roi is located."""
    arr = np.zeros((Ly, Lx), dtype=float)
    arr[self.ypix, self.xpix] = 1
    return arr

def stats_dicts_to_3d_array_(stat_dict, output_ops):
    arrays = []
    for i, stat in enumerate(stat_dict):
        array = to_array_(from_stat_dict(stat=stat), Ly=output_ops['Ly'], Lx=output_ops['Lx'])
        array *= i + 1
        arrays.append(array)
    return np.stack(arrays)

class Suite2pResultsExperiment:
    """used to run and further process suite2p processed data, and analysis associated with suite2p processed data."""

    # default ops dict for suite2p
    ops = {
        'batch_size': 2000,  # reduce if running out of RAM
        'fast_disk': os.path.expanduser('/mnt/sandbox/pshah/suite2p_tmp'), # used to store temporary binary file, defaults to save_path0 (set as a string NOT a list)
        # 'save_path0': '', # stores results, defaults to first item in data_path
        'delete_bin': True,  # whether to delete binary file after processing
        # main settings
        'nplanes': 1,  # each tiff has these many planes in sequence
        'nchannels': 1,  # each tiff has these many channels per plane
        'functional_chan': 1,  # this channel is used to extract functional ROIs (1-based)
        'diameter': 12,
        # this is the main parameter for cell detection, 2-dimensional if Y and X are different (e.g. [6 12])
        'tau': 1.26,  # this is the main parameter for deconvolution (1.25-1.5 for gcamp6s)
        'fs': 30,  # sampling rate (total across planes)
        # output settings
        'save_mat': True,  # whether to save output as matlab files
        'combined': True,  # combine multiple planes into a single result /single canvas for GUI
        # parallel settings
        'num_workers': 50,  # 0 to select num_cores, -1 to disable parallelism, N to enforce value
        'num_workers_roi': 0,  # 0 to select number of planes, -1 to disable parallelism, N to enforce value
        # registration settings
        'do_registration': True,  # whether to register data
        'nimg_init': 200,  # subsampled frames for finding reference image
        'maxregshift': 0.1,  # max allowed registration shift, as a fraction of frame max(width and height)
        'align_by_chan': 1,  # when multi-channel, you can align by non-functional channel (1-based)
        'reg_tif': True,  # whether to save registered tiffs
        'subpixel': 10,  # precision of subpixel registration (1/subpixel steps)
        # cell detection settings
        'connected': True,  # whether or not to keep ROIs fully connected (set to 0 for dendrites)
        'navg_frames_svd': 5000,  # max number of binned frames for the SVD
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
        # self.reg_tiff_path = None
        # self.ops_end = None
        # self.s2pResultsPath = None
        # self.__s2pResultExists = False
        # self.ops = None
        self.bad_frames: list = []  #: list of frames to discard during Suite2p ROI detection procedure
        print(f"\- ADDING .Suite2p module to Experiment object ... ", end='\r')

        ## initialize attr's
        # TODO add attr comments
        self.n_frames: int = 0  # total number of imaging frames in the Suite2p run

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

        self.neuropil_coeff = 0

        # set trials to run together in suite2p for Experiment
        self.trials = []
        tiff_paths_to_use_s2p: dict = trialsTiffsSuite2p
        for trial, path in tiff_paths_to_use_s2p.items():
            if path is not None:
                if os.path.exists(path):
                    self.trials.append(trial)
            else:
                tiff_paths_to_use_s2p.pop(trial)

        # self.trials = [*trialsTiffsSuite2p]
        self.tiff_paths_to_use_s2p: dict = tiff_paths_to_use_s2p
        assert len(
            self.trials) > 0, "no trials found to run suite2p, option available to provide list of trial IDs in " \
                              "`trialsSuite2P` "

        if s2pResultsPath is None:
            # initialize needed variables and attr's for future calling of s2pRun
            self.s2pResultsPath = None
            self.__s2pResultExists = False
        else:
            self.s2pResultsPath = s2pResultsPath
            try:
                self._retrieveSuite2pData(self.s2pResultsPath, neuropil_coeff=self.neuropil_coeff)
            except Exception:
                raise Exception(
                    f'Something went wrong while trying to load suite2p processed data from: {s2pResultsPath}')
            self.__s2pResultExists = True

        self.db: dict = {'fs': float(self.ops['fs']), 'diameter': self.ops['diameter'], 'batch_size': self.ops['batch_size'],
          'nimg_init': self.ops['batch_size'], 'nplanes': self.ops['nplanes'], 'nchannels': self.ops['nchannels'],
          'tiff_list': list(self.tiff_paths_to_use_s2p.values()),
          'save_folder': self.s2pResultsPath}


        # finish
        print(f"|- ADDED .Suite2p module to Experiment object. ")

    @property
    def _s2pResultExists(self):
        return self.__s2pResultExists

    # @_s2pResultExists.setter
    # def _s2pResultExists(self, val):
    #     self.__s2pResultExists = val


    def __repr__(self):
        if self._s2pResultExists:
            return f'Suite2p Results (Experiment level) Object, containing trials: \n\t{self.trials}'
        else:
            return f'Suite2p Results (Experiment level) Object, containing trials: \n\t{self.trials}. No Suite2p results loaded.'

    # noinspection PyTypeChecker
    def _retrieveSuite2pData(self, s2p_path: str = None, neuropil_coeff: float = 0.7):
        """processing of suite2p data from the current t-series
        :param s2p_path: s2pResultsPath to the directory containing suite2p outputs
        :param neuropil_coeff: choose to subtract neuropil or not when loading s2p traces
        :param save: choose to save data object or not
        """

        print(f'\----- Adding Suite2p results to .Suite2p module ...')

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
            self.n_frames = len(spks)

            self.ops: dict = np.load(os.path.join(s2p_path, 'ops.npy'), allow_pickle=True).item()
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
                f'|- Loaded {self.n_units} suite2p classified cells from plane {plane}, recorded for {round(self.raw[plane].shape[1] / self.ops["fs"], 2)} secs total, {self.n_frames} frames total')

        # consider replacing this and use returning properties
        if self.n_planes == 1:
            # print(f'*** plane 0 data ***')
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
            self.n_frames = self.spks.shape[1]
            # print(f'n_frames', self.n_frames)

        # read in other files
        self.output_ops: dict = np.load(Path(self.s2pResultsPath).joinpath('ops.npy'), allow_pickle=True).item()
        stats_file = Path(self.s2pResultsPath).joinpath('stat.npy')
        self.iscell = np.load(Path(self.s2pResultsPath).joinpath('iscell.npy'), allow_pickle=True)[:, 0].astype(bool)
        self.stat = np.load(stats_file, allow_pickle=True)

        self.im = stats_dicts_to_3d_array_(stat_dict=self.stat, output_ops=self.output_ops)
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

    @classmethod
    def subSuite2p(cls, trialsTiffs, s2pResultsPath):
        """
        Alternative constructor for Suite2pResultsExperiment class.
        """
        return cls(trialsTiffsSuite2p=trialsTiffs, s2pResultsPath=s2pResultsPath)

    def add_bad_frames(self, frames, bad_frames_npy_loc) -> None:
        """
        Add frames to bad_frames.npy file that will be used by Suite2p to ignore these frames during ROI detection.
        """

        self.bad_frames.extend(frames)
        self.bad_frames = list(np.unique(self.bad_frames))

        if len(self.bad_frames) > 0:
            print('\- Appending a total of ', len(frames), 'to bad_frames.npy', f"\n\t total bad_frames: {len(self.bad_frames)}")
            np.save('%s/bad_frames.npy' % bad_frames_npy_loc,
                    self.bad_frames)  # save to npy file and remember to move npy file to tiff folder before running with suite2p


    def setupForSuite2p(self, trialsSuite2P: list, TrialsInformation: dict, **kwargs):
        """

        :param trialsSuite2P: list of trials to use for Suite2p setup
        :param TrialsInformation: dictionary mapping trial ID to tiff_path to use for Suite2p processing of that trial
        :param kwargs:
        """
        if trialsSuite2P != self.trials:

            trialsTiffs = {}
            for trial in trialsSuite2P:
                trialsTiffs[trial] = TrialsInformation[trial]['tiff_path']

            Suite2p_obj = self.subSuite2p(trialsTiffs=trialsSuite2P, s2pResultsPath=self.s2pResultsPath, subtract_neuropil=self.subtract_neuropil)
        else:
            Suite2p_obj = self
        # self.trials = trialsSuite2P if trialsSuite2P else self.trials



        for trial in Suite2p_obj.trials:
            Suite2p_obj.tiff_paths_to_use_s2p[trial] = TrialsInformation[trial]['tiff_path']

        # load the first tiff in Suite2p_obj.trials to collect default metainformation about imaging setup parameters
        trial = Suite2p_obj.trials[0]
        from packerlabimaging import import_obj
        trialobj = import_obj(TrialsInformation[trial]['analysis_object_information']['pkl path'])

        # set imaging parameters using defaults or kwargs if provided
        fps = trialobj.fps if 'fs' not in [*kwargs] else kwargs['fs']
        n_planes = trialobj.n_planes if 'n_planes' not in [*kwargs] else kwargs['n_planes']
        pix_sz_x = trialobj.pix_sz_x if 'pix_sz_x' not in [*kwargs] else kwargs['pix_sz_x']
        pix_sz_y = trialobj.pix_sz_y if 'pix_sz_y' not in [*kwargs] else kwargs['pix_sz_y']
        frame_x = trialobj.frame_x if 'frame_x' not in [*kwargs] else kwargs['frame_x']
        frame_y = trialobj.frame_y if 'frame_y' not in [*kwargs] else kwargs['frame_y']
        n_channels = kwargs['n_channels'] if 'n_channels' in [
            *kwargs] else 1  # default is 1 channel imaging in .tiffs for suite2p

        # setup ops dictionary
        Suite2p_obj.ops['fs'] = fps / n_planes
        diameter_x = 13 / pix_sz_x
        diameter_y = 13 / pix_sz_y
        Suite2p_obj.ops['diameter'] = int(diameter_x), int(diameter_y) if diameter_y != diameter_x else diameter_x

        # set other ops parameters if provided in kwargs:
        for key in [*kwargs]:
            if key in [*Suite2p_obj.ops]:
                Suite2p_obj.ops[key] = kwargs[key]

        # calculate batch size to use in Suite2p run
        batch_size = Suite2p_obj.ops['batch_size'] * (262144 / (
                frame_x * frame_y))  # larger frames will be more RAM intensive, scale user batch size based on num pixels in 512x512 images

        # setup db dict
        Suite2p_obj.db = {'fs': float(Suite2p_obj.ops['fs']), 'diameter': Suite2p_obj.ops['diameter'], 'batch_size': int(batch_size),
              'nimg_init': int(Suite2p_obj.ops['batch_size']), 'nplanes': n_planes, 'nchannels': n_channels,
              'tiff_list': list(Suite2p_obj.tiff_paths_to_use_s2p.values()), 'data_path': Suite2p_obj.dataPath,  # TODO need to figure out how to more appropriately bring in the dataPath here
              'save_folder': Suite2p_obj.s2pResultsPath}


    # # suite2p methods
    # def s2pRun(self, user_batch_size=2000, trialsSuite2P: list = None,
    #            **kwargs):  ## TODO gotta specify # of planes somewhere here
    #     """run suite2p for an Experiment object, using trials specified in current experiment object, using the attributes
    #     determined directly from the experiment object.
    #
    #     :param self: a packerlabimaging.main.Experiment object
    #     :param user_batch_size: batch size for suite2p registration stage - adjust if running into MemmoryError while running suite2p
    #     :param trialsSuite2P: list of trialIDs from experiment to use in running suite2p
    #     """
    #
    #     print(f'db: \n\t{self.db}')
    #
    #
    #     self.trials = trialsSuite2P if trialsSuite2P else self.trials
    #     #
    #     tiffs_paths_to_use_s2p = []
    #     for trial in self.trials:
    #         tiffs_paths_to_use_s2p.append(self.TrialsInformation[trial]['tiff_path'])
    #     #
    #     # # load the first tiff in self.trials to collect default metainformation about imaging setup parameters
    #     # trial = self.trials[0]
    #     # from packerlabimaging import import_obj
    #     # trialobj = import_obj(self.TrialsInformation[trial]['analysis_object_information']['pkl path'])
    #     #
    #     # # set imaging parameters using defaults or kwargs if provided
    #     # fps = trialobj.fps if 'fs' not in [*kwargs] else kwargs['fs']
    #     # n_planes = trialobj.n_planes if 'n_planes' not in [*kwargs] else kwargs['n_planes']
    #     # pix_sz_x = trialobj.pix_sz_x if 'pix_sz_x' not in [*kwargs] else kwargs['pix_sz_x']
    #     # pix_sz_y = trialobj.pix_sz_y if 'pix_sz_y' not in [*kwargs] else kwargs['pix_sz_y']
    #     # frame_x = trialobj.frame_x if 'frame_x' not in [*kwargs] else kwargs['frame_x']
    #     # frame_y = trialobj.frame_y if 'frame_y' not in [*kwargs] else kwargs['frame_y']
    #     # n_channels = kwargs['n_channels'] if 'n_channels' in [
    #     #     *kwargs] else 1  # default is 1 channel imaging in .tiffs for suite2p
    #     #
    #     # # setup ops dictionary
    #     # ops = self.ops
    #     #
    #     # ops['fs'] = fps / n_planes
    #     # diameter_x = 13 / pix_sz_x
    #     # diameter_y = 13 / pix_sz_y
    #     # ops['diameter'] = int(diameter_x), int(diameter_y) if diameter_y != diameter_x else diameter_x
    #     # self.user_batch_size = user_batch_size
    #     # ops['batch_size'] = self.user_batch_size
    #     # batch_size = self.user_batch_size * (262144 / (
    #     #         frame_x * frame_y))  # larger frames will be more RAM intensive, scale user batch size based on num pixels in 512x512 images
    #     #
    #     # # set other ops parameters if provided in kwargs:
    #     # for key in [*kwargs]:
    #     #     if key in [*ops]:
    #     #         ops[key] = kwargs[key]
    #     #
    #     # # setup db dict
    #     # db = {'fs': float(ops['fs']), 'diameter': ops['diameter'], 'batch_size': int(batch_size),
    #     #       'nimg_init': int(batch_size), 'nplanes': n_planes, 'nchannels': n_channels,
    #     #       'tiff_list': tiffs_paths_to_use_s2p, 'data_path': self.dataPath,
    #     #       'save_folder': self.s2pResultsPath}
    #     #
    #     # print(f'db: \n\t{db}')
    #
    #
    #     t1 = time.time()
    #     opsEnd = run_s2p(ops=self.ops, db=self.db)
    #     t2 = time.time()
    #     print('Total time to run suite2p was {}'.format(t2 - t1))
    #
    #     # update self.ops
    #     self.ops = opsEnd
    #
    #     self._s2pResultExists = True
    #     self.s2pResultsPath = self.s2pResultsPath + '/plane0/'  ## need to further debug that the flow of the suite2p s2pResultsPath makes sense

    @staticmethod
    def s2pRun(expobj, trialsSuite2P: Union[list, str] = 'all'):
        s2pRun(expobj, trialsSuite2P=trialsSuite2P)

    def s2pROIsTiff(self, save_path, cell_ids: Union[str, list] = 'all'):
        """save a TIFF image of the suite2p ROIs masks."""
        # todo test function
        # os.chdir(self.s2pResultsPath)
        # stat = np.load('stat.npy', allow_pickle=True)
        # ops = np.load('ops.npy', allow_pickle=True).item()
        # iscell = np.load('iscell.npy', allow_pickle=True)

        mask_img = np.zeros((self.output_ops['Ly'], self.output_ops['Lx']), dtype='uint8')

        cell_ids = self.cell_id if cell_ids == 'all' else cell_ids
        for n in range(0, len(self.iscell)):
            if n in cell_ids:
                ypix = self.stat[n]['ypix']
                xpix = self.stat[n]['xpix']
                mask_img[ypix, xpix] = np.random.randint(10, 255)

        save_array_to_tiff(save_path=save_path, data=mask_img)


    # PROCESSING OF SUITE2P OUTPUT
    def filterROIs(self, overlapthreshold: int = None, skewthreshold: float = None, footprintthreshold: float = None, npixthreshold: float = None,
                   aspectratiothreshold: float = None, classifierthreshold: float = None, boundarybox: np.array = None):

        """todo
        Filter Suite2p ROIs based on a variety of metrics of each ROI. Filters can be combined by providing a value to each filter.

        :param skewthreshold:
        :param footprintthreshold:
        :param npixthreshold:
        :param aspectratiothreshold:
        :param classifierthreshold:
        :param boundarybox:
        :param overlapthreshold: percentage threshold for filtering out overlapping ROIs. if overlap between two or more ROIs exceeds a specified percentage threshold, then remove the ROI with smaller npix.

        """



class Suite2pResultsTrial:
    """used to collect and store suite2p processed data for one trial - out of overall experiment."""

    def __init__(self,  s2pExp: Suite2pResultsExperiment, trial_frames: tuple):
        """
        Connecting Suite2p results from a specific trial (which spans trial_frames out of the overall Suite2p run) to that trial.

        :param s2pExp:
        :param trial_frames:


        """

        print(f"\n\----- ADDING .Suite2p module to trial ... ", end ='\r')

        ## initializing attributes to collect Suite2p results for this specific trial
        # self.dfof: list(np.ndarray) = []  # array of dFF normalized Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
        # self.__raw: List[np.ndarray] = self.raw
        # self.__spks: List[np.ndarray] = self.spks
        # self.__neuropil: List[np.ndarray] = self.neuropil

        self.trial_frames = trial_frames  # tuple of first and last frame (out of the overall suite2p run) corresponding to the present trial

        # self.suite2p_overall = suite2p_experiment_obj
        print(
            f"\t|- current trial frames: {trial_frames} out of {s2pExp.n_frames} total frames processed through suite2p")
        if s2pExp.s2pResultsPath: self._get_suite2pResults(s2pExp)
        else: self.__s2pResultExists = False

        # s2pResultsPath = suite2p_experiment_obj.s2pResultsPath if hasattr(suite2p_experiment_obj, 's2pResultsPath') else None
        # Suite2pResultsExperiment.__init__(self, trialsSuite2p = suite2p_experiment_obj.trials, s2pResultsPath=s2pResultsPath,
        #                                   subtract_neuropil=suite2p_experiment_obj.subtract_neuropil)

        print(f"\n\----- ADDED .Suite2p module to trial. ", end ='\r')

    def __repr__(self):
        if self.s2pResultExists:
            return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} frames x {self.n_units} s2p ROIs'
        else:
            return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} frames. No Suite2p Results loaded.'


    @property
    def s2pResultExists(self):
        return self.__s2pResultExists

    def _get_suite2pResults(self, s2pExp: Suite2pResultsExperiment):  # TODO complete code for getting suite2p results for trial
        """crop suite2p data for frames only for the present trial"""

        # for attr in ['n_units', 'cell_id', 'cell_plane', 'cell_x', 'cell_y', 'xoff', 'yoff', 'raw', 'spks', 'neuropil', 'stat']:
        #     try:
        #         setattr(self, attr, getattr(self.suite2p_overall, attr))
        #     except AttributeError:
        #         pass

        if s2pExp.n_planes == 1:
            self.cell_id = s2pExp.cell_id
            self.stat = s2pExp.stat[0]
            self.output_ops = s2pExp.output_ops
            self.n_units = s2pExp.n_units
            self.iscell = s2pExp.iscell

            self.raw = s2pExp.raw[:, self.trial_frames[0]:self.trial_frames[
                1]]  # array of raw Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
            self.spks = s2pExp.spks[:, self.trial_frames[0]:self.trial_frames[
                1]]  # array of deconvolved spks values from suite2p output [num cells x length of imaging acquisition], one per plane
            self.neuropil = s2pExp.neuropil[:, self.trial_frames[0]:self.trial_frames[
                1]]  # array of neuropil Flu values from suite2p output [num cells x length of imaging acquisition], one per plane

            self.__s2pResultExists = True

        else:
            for plane in range(s2pExp.n_planes):
                self.raw.append(s2pExp.raw[plane][:, self.trial_frames[0]:self.trial_frames[1]])
                self.spks.append(s2pExp.spks[plane][:, self.trial_frames[0]:self.trial_frames[1]])
                self.neuropil.append(s2pExp.neuropil[plane][:, self.trial_frames[0]:self.trial_frames[1]])

                # self.dfof.append(normalize_dff(self.raw[plane]))  # calculate df/f based on relevant frames
                self.__s2pResultExists = True

        self.im = stats_dicts_to_3d_array_(stat_dict=s2pExp.stat, output_ops=s2pExp.output_ops)
        self.im[self.im == 0] = np.nan


    def makeFrameAverageTiff(self, reg_tif_dir: str, frames: Union[int, list], peri_frames: int = 100,
                             save_dir: str = None, to_plot=False):
        """Creates, plots and/or saves an average image of the specified number of peri-frames around the given frame from the suite2p registered TIFF file.
        """

        # read in registered tiff
        reg_tif_folder = reg_tif_dir
        assert os.path.exists(reg_tif_dir), '`reg_tif_dir` could not be found.'
        reg_tif_list = os.listdir(reg_tif_dir)
        reg_tif_list.sort()

        if type(frames) == int:
            frames = [frames]

        frames_s2prun_num = [(frame + self.trial_frames[0]) for frame in frames]

        for idx, frame in enumerate(frames_s2prun_num):
            batch_num = frame // self.output_ops['batch_size']

            tif_path = reg_tif_folder + reg_tif_list[batch_num]

            im_batch_reg = tf.imread(tif_path, key=range(0, self.output_ops['batch_size']))

            frame_num_batch = frame - (batch_num * self.output_ops['batch_size'])

            if frame_num_batch < peri_frames // 2:
                peri_frames_low = frame_num_batch
            else:
                peri_frames_low = peri_frames // 2
            peri_frames_high = peri_frames // 2
            im_sub_reg = im_batch_reg[frame_num_batch - peri_frames_low: frame_num_batch + peri_frames_high]

            avg_sub = np.mean(im_sub_reg, axis=0)

            # convert to 8-bit
            from packerlabimaging.utils.utils import convert_to_8bit
            avg_sub = convert_to_8bit(avg_sub, 0, 255)

            if save_dir:
                if '.tif' in save_dir:
                    from packerlabimaging.utils.utils import return_parent_dir
                    save_dir = return_parent_dir(save_dir) + '/'
                save_path = save_dir + f'/{frames[idx]}_s2preg_frame_avg.tif'
                os.makedirs(save_dir, exist_ok=True)

                print(f"\t\- Saving averaged s2p registered tiff for frame: {frames[idx]}, to: {save_path}")
                tf.imwrite(save_path, avg_sub, photometric='minisblack')

            if to_plot:
                plt.imshow(avg_sub, cmap='gray')
                plt.suptitle(f'{peri_frames} frames avg from s2p reg tif, frame: {frames[idx]}')
                plt.show()  # just plot for now to make sure that you are doing things correctly so far






#### archiving away for now - trying to switch to an approach that doesn't inherit from parent suite2p obj.
class Suite2pResultsTrial_(Suite2pResultsExperiment):
    """used to collect and store suite2p processed data for one trial - out of overall experiment."""

    def __init__(self, trialsTiffsSuite2p: dict, trial_frames: tuple, s2pResultsPath: Optional[str] = None):
        """
        Connecting Suite2p results from a specific trial (which spans trial_frames out of the overall Suite2p run) to that trial.

        :param trialsTiffsSuite2p:
        :param dataPath:
        :param trial_frames:
        :param s2pResultsPath:
        :param subtract_neuropil:
        """
        # todo need to reconsider this heirarchy; it saves in the entire suite2p result for every trial....
        super().__init__(trialsTiffsSuite2p, s2pResultsPath)  # - TODO it really is confusing to be passing in all trials for a s2p results obj that should be restricted to just one trial

        print(f"\n\----- ADDING .Suite2p module to trial ... ", end ='\r')

        ## initializing attributes to collect Suite2p results for this specific trial
        # self.dfof: list(np.ndarray) = []  # array of dFF normalized Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
        # self.__raw: List[np.ndarray] = self.raw
        # self.__spks: List[np.ndarray] = self.spks
        # self.__neuropil: List[np.ndarray] = self.neuropil

        self.trial_frames = trial_frames  # tuple of first and last frame (out of the overall suite2p run) corresponding to the present trial

        # self.suite2p_overall = suite2p_experiment_obj
        print(
            f"\t|- current trial frames: {trial_frames} out of {self.n_frames} total frames processed through suite2p")
        self._get_suite2pResults() if self.s2pResultsPath else None

        # s2pResultsPath = suite2p_experiment_obj.s2pResultsPath if hasattr(suite2p_experiment_obj, 's2pResultsPath') else None
        # Suite2pResultsExperiment.__init__(self, trialsSuite2p = suite2p_experiment_obj.trials, s2pResultsPath=s2pResultsPath,
        #                                   subtract_neuropil=suite2p_experiment_obj.subtract_neuropil)

        print(f"\n\----- ADDED .Suite2p module to trial. ", end ='\r')

    def __repr__(self):
        if self._s2pResultExists:
            return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} frames x {self.n_units} s2p ROIs'
        else:
            return f'Suite2p Results (trial level) Object, {self.trial_frames[1] - self.trial_frames[0]} frames. No Suite2p Results loaded.'


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

            self.raw = self.raw[:, self.trial_frames[0]:self.trial_frames[
                1]]  # array of raw Flu values from suite2p output [num cells x length of imaging acquisition], one per plane
            self.spks = self.spks[:, self.trial_frames[0]:self.trial_frames[
                1]]  # array of deconvolved spks values from suite2p output [num cells x length of imaging acquisition], one per plane
            self.neuropil = self.neuropil[:, self.trial_frames[0]:self.trial_frames[
                1]]  # array of neuropil Flu values from suite2p output [num cells x length of imaging acquisition], one per plane

            self.__s2pResultExists = True

        else:
            for plane in range(self.n_planes):
                self.raw.append(self.raw[plane][:, self.trial_frames[0]:self.trial_frames[1]])
                self.spks.append(self.spks[plane][:, self.trial_frames[0]:self.trial_frames[1]])
                self.neuropil.append(self.neuropil[plane][:, self.trial_frames[0]:self.trial_frames[1]])

                # self.dfof.append(normalize_dff(self.raw[plane]))  # calculate df/f based on relevant frames
                self.__s2pResultExists = True

    def makeFrameAverageTiff(self, reg_tif_dir: str, frames: Union[int, list], peri_frames: int = 100,
                             save_dir: str = None,
                             to_plot=False):
        """Creates, plots and/or saves an average image of the specified number of peri-frames around the given frame from the suite2p registered TIFF file.
        """

        # read in registered tiff
        reg_tif_folder = reg_tif_dir
        assert os.path.exists(reg_tif_dir), '`reg_tif_dir` could not be found.'
        reg_tif_list = os.listdir(reg_tif_dir)
        reg_tif_list.sort()

        if type(frames) == int:
            frames = [frames]

        frames_s2prun_num = [(frame + self.trial_frames[0]) for frame in frames]

        for idx, frame in enumerate(frames_s2prun_num):
            batch_num = frame // self.ops['batch_size']

            tif_path = reg_tif_folder + reg_tif_list[batch_num]

            im_batch_reg = tf.imread(tif_path, key=range(0, self.ops['batch_size']))

            frame_num_batch = frame - (batch_num * self.ops['batch_size'])

            if frame_num_batch < peri_frames // 2:
                peri_frames_low = frame_num_batch
            else:
                peri_frames_low = peri_frames // 2
            peri_frames_high = peri_frames // 2
            im_sub_reg = im_batch_reg[frame_num_batch - peri_frames_low: frame_num_batch + peri_frames_high]

            avg_sub = np.mean(im_sub_reg, axis=0)

            # convert to 8-bit
            from packerlabimaging.utils.utils import convert_to_8bit
            avg_sub = convert_to_8bit(avg_sub, 0, 255)

            if save_dir:
                if '.tif' in save_dir:
                    from packerlabimaging.utils.utils import return_parent_dir
                    save_dir = return_parent_dir(save_dir) + '/'
                save_path = save_dir + f'/{frames[idx]}_s2preg_frame_avg.tif'
                os.makedirs(save_dir, exist_ok=True)

                print(f"\t\- Saving averaged s2p registered tiff for frame: {frames[idx]}, to: {save_path}")
                tf.imwrite(save_path, avg_sub, photometric='minisblack')

            if to_plot:
                plt.imshow(avg_sub, cmap='gray')
                plt.suptitle(f'{peri_frames} frames avg from s2p reg tif, frame: {frames[idx]}')
                plt.show()  # just plot for now to make sure that you are doing things correctly so far

