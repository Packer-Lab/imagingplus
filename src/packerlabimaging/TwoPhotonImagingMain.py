import os
import time
import datetime
import re
from typing import Optional, TypedDict

import pickle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import xml.etree.ElementTree as ET
import tifffile as tf

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)
from packerlabimaging.utils.utils import SaveDownsampledTiff, normalize_dff
from packerlabimaging.utils.classes import UnavailableOptionError, TrialsInformation
from packerlabimaging.processing.paq import PaqData
from .processing import suite2p, anndata as ad
from .utils.imagingMetadata import PrairieViewMetadata, ImagingMetadata


class TwoPhotonImagingMetainfo(TypedDict, total=False):
    data: str
    trial_id: str
    exp_id: str
    t_series_id: str
    TrialsInformation: TrialsInformation


class TwoPhotonImagingTrial:
    """Two Photon Imaging Experiment Data Analysis Workflow."""

    def __init__(self, metainfo: dict, analysis_save_path: str, microscope: str = 'Bruker', paq_options: dict = None,
                 **kwargs):
        """
        TODO update function docstring for approp args
        :param metainfo: TypedDict containing meta-information field needed for this experiment. Please see TwoPhotonImagingMetainfo for type hints on accepted keys.
        :param paq_options: TypedDict containing meta-information about .paq file associated with current trial
        :param analysis_save_path: path of where to save the experiment analysis object
        :param microscope: name of microscope used to record imaging (options: "Bruker" (default), "other")
        :param imagingMicroscopeMetadata: provide ImagingMetadata object (see ImagingMetadata class).
        :param suite2p_experiment_obj: provide Suite2p Experiment Object as variable in order to process Suite2p data for current trial
        :param total_frames_stitched: provide frame number on which current trial starts in Suite2p Experiment Object
        """

        # Initialize Attributes:

        # # Imaging acquisition attr's: - refactored to _prairieview
        # self.n_frames: int = 0  # number of imaging frames in the current trial
        # self.fps = None  # rate of imaging acquisition (frames per second)
        # self.frame_x = None  # num of pixels in the x direction of a single frame
        # self.frame_y = None  # num of pixels in the y direction of a single frame
        # self.n_planes = None  # num of FOV planes in imaging acquisition
        # self.pix_sz_x = None  # size of a single imaging pixel in x direction (microns)
        # self.pix_sz_y = None  # size of a single imaging pixel in y direction (microns)
        # self.scan_x = None  # TODO ROB - not sure what the comment for this is
        # self.scan_y = None  # TODO ROB - not sure what the comment for this is
        # self.zoom: float = 0.0 # zoom level on Bruker microscope
        # self.last_good_frame = None  # indicates when the last good frame was during the t-series recording, if nothing was wrong the value is 0, otherwise it is >0 and that indicates that PV is not sure what happened after the frame listed, but it could be corrupt data

        if 'date' in [*metainfo] and 'trial_id' in [*metainfo] and 'exp_id' in [*metainfo] and 't_series_id' in [
            *metainfo] and 'TrialsInformation' in [*metainfo]:
            self.__metainfo = metainfo
        else:
            raise ValueError(
                "dev error: metainfo argument must contain the minimum fields: 'date', 'trial_id', 'exp_id', 't_series_id', 'trialInformation dict")

        if os.path.exists(metainfo['TrialsInformation']['tiff_path']):
            self.tiff_path = metainfo['TrialsInformation']['tiff_path']  #: path to the tiff for current trial
        else:
            raise FileNotFoundError(f"tiff_path does not exist: {metainfo['TrialsInformation']['tiff_path']}")


        if os.path.exists(metainfo['TrialsInformation']['tiff_path']):
            self.tiff_path = metainfo['TrialsInformation']['tiff_path']  #: path to the tiff for current trial
        else:
            raise FileNotFoundError(f"tiff_path does not exist: {metainfo['TrialsInformation']['tiff_path']}")

        if 'paq_path' in [*paq_options]:
            paq_path = paq_options['paq_path']
            if os.path.exists(paq_path):
                self._use_paq = True
            else:
                raise FileNotFoundError(f"paq_path does not exist: {paq_path}")
        else:
            self._use_paq = False


        # set and create analysis save path directory
        self.save_dir = analysis_save_path  #: path to the directory to save outputs from analysis
        os.makedirs(self.save_dir, exist_ok=True)
        self._pkl_path = f"{self.save_dir}{metainfo['date']}_{metainfo['trial_id']}.pkl"

        self.save_pkl(pkl_path=self.pkl_path)  # save experiment object to pkl_path

        # get imaging setup parameters
        if 'Bruker' in microscope:
            self.imparams = self._getImagingParameters(
                microscope='Bruker')  #: Imaging Microscope Metadata submodule for trial
        elif 'imagingMicroscopeMetadata' in [*kwargs]:
            self.imparams = self._getImagingParameters(metadata=kwargs['imagingMicroscopeMetadata'])
        else:
            Warning(f"NO imaging microscope parameters set. ")

        # run Paq processing if paq_path is provided for trial
        if self._use_paq:
            frame_channel = paq_options['frame_channel'] if 'frame_channel' in [
                *paq_options] else KeyError(
                'No frame_channel specified for .paq processing')  # channel on Paq file to read for determining stims
            self.Paq = self._paqProcessingTwoPhotonImaging(paq_path=paq_options['paq_path'],
                                                           frame_channel=frame_channel)  #: Paq data submodule for trial

        # collect mean FOV Trace
        self.meanFluImg, self.meanFovFluTrace = self.meanRawFluTrace()  #: mean image and mean FOV fluorescence trace

        # if provided Suite2p information provided, add Suite2p results for trial; + run further processing.
        if 'suite2p_experiment_obj' in [*kwargs] and 'total_frames_stitched' in [*kwargs]:
            from packerlabimaging.processing.suite2p import Suite2pResultsExperiment
            s2p_expobj: Suite2pResultsExperiment = kwargs['suite2p_experiment_obj']
            #: Suite2p results submodule for trial
            self.Suite2p = suite2p.Suite2pResultsTrial(trialsTiffsSuite2p=s2p_expobj.tiff_paths_to_use_s2p,
                                                       s2pResultsPath=s2p_expobj.s2pResultsPath,
                                                       subtract_neuropil=s2p_expobj.subtract_neuropil,
                                                       trial_frames=(kwargs['total_frames_stitched'], kwargs[
                                                           'total_frames_stitched'] + self.imparams.n_frames))  # use trial obj's current trial frames

            # normalize dFF for raw Flu
            self.dFF = self.dfof()  #: dFF normalized Suite2p data for trial

            # create annotated data object
            self.data = self.create_anndata()  #: anndata storage submodule

        # SAVE Trial OBJECT
        self.save()
        print(f'\----- CREATING TwoPhotonImagingTrial for trial: {metainfo["trial_id"]},  {metainfo["t_series_id"]}')

    def __str__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
        else:
            lastmod = "(unsaved pkl object)"
        return repr(f"({self.t_series_name}) TwoPhotonImagingTrial experimental data object, last saved: {lastmod}")

    def __repr__(self):
        return repr(f"TwoPhotonImagingTrial experimental data object")

    @property
    def fig_save_path(self):
        """create path for saving figures"""
        today_date = datetime.today().strftime('%Y-%m-%d')
        return self.save_dir + f'Results_fig/{today_date}/'

    @fig_save_path.setter
    def fig_save_path(self, value: str):
        """set new default fig save path for data object"""
        self.fig_save_path = value

    @property
    def date(self):
        """date of the experiment datacollection"""
        return self.__metainfo['date']

    @property
    def expID(self):
        """experiment ID of current trial object"""
        return self.__metainfo['exp_id']

    @property
    def trialID(self):
        """trial ID of current trial object"""
        return self.__metainfo['trial_id']

    @property
    def t_series_name(self):
        if 't_series_id' in self.__metainfo.keys():
            return f"{self.__metainfo['t_series_id']}"
        elif "exp_id" in self.__metainfo.keys() and "trial_id" in self.__metainfo.keys():
            return f'{self.__metainfo["exp_id"]} {self.__metainfo["trial_id"]}'
        else:
            raise ValueError('no information found to retrieve t series id')

    @property
    def tiff_path_dir(self):
        return self.tiff_path[:[(s.start(), s.end()) for s in re.finditer('/', self.tiff_path)][-1][
            0]]  # this is the directory where the Bruker xml files associated with the 2p imaging TIFF are located

    @property
    def pkl_path(self):
        """path in Analysis folder to save pkl object"""
        return self._pkl_path

    @pkl_path.setter
    def pkl_path(self, path: str):
        self._pkl_path = path

    def _getImagingParameters(self, metadata: Optional[dict] = None, microscope: Optional[str] = 'Bruker'):
        """retrieves imaging metadata parameters. If using Bruker microscope and PrairieView, then _prairieview module is used to collect this data.

        :param microscope: name of the microscope, currently the only supported microscope for parsing metadata directly is Bruker/PrairieView imaging setup.
        """
        if 'Bruker' in microscope:
            return PrairieViewMetadata(tiff_path_dir=self.tiff_path_dir)
        else:
            try:
                return ImagingMetadata(**metadata)
            except TypeError:
                Exception('required key not present in metadata')

    @property
    def frame_clock(self):
        if hasattr(self.Paq, 'frame_times'):
            return self.Paq.frame_times
        else:
            raise ValueError('Frame clock timings couldnt be retrieved from .Paq submodule.')

    @property
    def n_frames(self):
        try:
            return self.imparams.n_frames
        except AttributeError:
            return -1

    def _paqProcessingTwoPhotonImaging(self, paq_path, frame_channel):
        """
        Add and further process paq data for current trial.
        :param paq_path: path to .paq file
        :param frame_channel: channel to use for measuring frame times from .paq data

        :return: PAQ data object
        """

        print(f"\n\- ADDING PAQ MODULE ... ")
        paq_data_obj, paqdata = PaqData.import_paqdata(paq_path=paq_path)
        print(f"\n\- PROCESSING PAQDATA ... ")
        assert frame_channel in paq_data_obj.paq_channels, print(f"{frame_channel} not found in channels in .paq data.")
        paq_data_obj.frame_times_channame = frame_channel
        paq_data_obj.frame_times = paq_data_obj.paq_frame_times(frame_channel=frame_channel)
        paq_data_obj.sparse_paq_data = paq_data_obj.get_sparse_paq(frame_clock=paq_data_obj.frame_times)

        return paq_data_obj

    def stitch_s2p_reg_tiff(self): ## TODO refactoring in new code from the Suite2p class script?
        assert self.Suite2p._s2pResultExists, UnavailableOptionError('stitch_s2p_reg_tiff')

        tif_path_save2 = self.save_dir + f'reg_tiff_{self.t_series_name}_r.tif'

        start = self.Suite2p.trial_frames[0] // 2000  # 2000 because that is the batch size for suite2p run
        end = self.Suite2p.trial_frames[1] // 2000 + 1

        if not os.path.exists(tif_path_save2):
            with tf.TiffWriter(tif_path_save2, bigtiff=True) as tif:
                with tf.TiffFile(self.Suite2p.reg_tiff_path, multifile=False) as input_tif:
                    print('cropping registered tiff')
                    data = input_tif.asarray()
                    print('shape of stitched tiff: ', data.shape)
                reg_tif_crop = data[self.Suite2p.trial_frames[0] - start * self.Suite2p.s2p_run_batch: self.Suite2p.trial_frames[1] - (
                        self.Suite2p.trial_frames - start * self.Suite2p.s2p_run_batch)]
                print('saving cropped tiff ', reg_tif_crop.shape)
                tif.write(reg_tif_crop)

    def create_anndata(self):
        """
        Creates annotated data (see anndata library for more information on AnnotatedData) object based around the Ca2+ matrix of the imaging trial.

        """

        if self.Suite2p._s2pResultExists and self.Paq:
            # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
            # build dataframe for obs_meta from suite2p stat information
            obs_meta = pd.DataFrame(
                columns=['original_index', 'footprint', 'mrs', 'mrs0', 'compact', 'med', 'npix', 'radius',
                         'aspect_ratio', 'npix_norm', 'skew', 'std'], index=range(len(self.Suite2p.stat)))
            for idx, __stat in enumerate(self.Suite2p.stat):
                for __column in obs_meta:
                    obs_meta.loc[idx, __column] = __stat[__column]

            # build numpy array for multidimensional obs metadata
            obs_m = {'ypix': [],
                     'xpix': []}
            for col in [*obs_m]:
                for idx, __stat in enumerate(self.Suite2p.stat):
                    obs_m[col].append(__stat[col])
                obs_m[col] = np.asarray(obs_m[col])

            # SETUP THE VARIABLES ANNOTATIONS TO USE IN anndata
            # build dataframe for var annot's from Paq file
            var_meta = pd.DataFrame(index=[self.Paq.frame_times_channame], columns=range(self.imparams.n_frames))
            for fr_idx in range(self.imparams.n_frames):
                var_meta.loc[self.Paq.frame_times_channame, fr_idx] = self.Paq.frame_times[fr_idx]

            # BUILD LAYERS TO ADD TO anndata OBJECT
            layers = {'dFF': self.dFF
                      }

            print(f"\n\----- CREATING annotated data object using AnnData:")
            _data_type = 'Suite2p Raw (neuropil substracted)' if self.Suite2p.subtract_neuropil else 'Suite2p Raw'
            adata = ad.AnnotatedData(X=self.Suite2p.raw, obs=obs_meta, var=var_meta.T, obsm=obs_m, layers=layers,
                                     data_label=_data_type)

            print(f"\n{adata}")
            return adata

        else:
            Warning(
                'did not create anndata. anndata creation only available if experiments were processed with suite2p and .Paq file(s) provided for temporal synchronization')

    def dfof(self):
        """(delta F)/F normalization of raw Suite2p data of trial."""
        if self.Suite2p._s2pResultExists:
            dFF = normalize_dff(self.Suite2p.raw)
            return dFF

    def importTrialTiff(self) -> np.ndarray:
        """
        Import current trial's microscope imaging tiff in full.

        :return: imaging tiff as numpy array
        """
        print(f"\n\- loading raw TIFF file from: {self.tiff_path}", end='\r')
        im_stack = tf.imread(self.tiff_path, key=range(self.imparams.n_frames))
        print('|- Loaded experiment tiff of shape: ', im_stack.shape)

        return im_stack

    def meanRawFluTrace(self):
        """
        Collects the raw mean of FOV fluorescence trace across the t-series.

        :return: mean fluorescence trace
        """
        im_stack = self.importTrialTiff()

        print('\n-----collecting mean raw flu trace from tiff file...')
        mean_flu_img = np.mean(im_stack, axis=0)
        mean_fov_flutrace = np.mean(np.mean(im_stack, axis=1), axis=1)

        return mean_flu_img, mean_fov_flutrace

    def makeDownsampledTiff(self):
        """Import current trial tiff, create downsampled tiff and save in default analysis directory."""

        stack = self.importTrialTiff()
        SaveDownsampledTiff(stack=stack, save_as=f"{self.save_dir}/{self.date}_{self.trialID}_downsampled.tif")

    def plotSingleImageFrame(self, frame_num: int = 0, title: str = None):
        """
        plots an image of a single specified tiff frame after reading using tifffile.
        :param frame_num: frame # from 2p imaging tiff to show (default is 0 - i.e. the first frame)
        :param title: (optional) give a string to use as title
        :return: matplotlib imshow plot
        """
        stack = tf.imread(self.tiff_path, key=frame_num)
        plt.imshow(stack, cmap='gray')
        plt.suptitle(title) if title is not None else plt.suptitle(f'frame num: {frame_num}')
        plt.show()
        return stack

    @staticmethod
    def normalize_dff(arr, threshold_pct=20, threshold_val=None):
        """normalize given array (cells x time) to the mean of the fluorescence values below given threshold. Threshold
        will refer to the that lower percentile of the given trace."""

        if arr.ndim == 1:
            if threshold_val is None:
                a = np.percentile(arr, threshold_pct)
                mean_ = arr[arr < a].mean()
            else:
                mean_ = threshold_val
            # mean_ = abs(arr[arr < a].mean())
            new_array = ((arr - mean_) / mean_) * 100
            if np.isnan(new_array).any() == True:
                Warning('Cell (unknown) contains nan, normalization factor: %s ' % mean_)

        else:
            new_array = np.empty_like(arr)
            for i in range(len(arr)):
                if threshold_val is None:
                    a = np.percentile(arr[i], threshold_pct)
                else:
                    a = threshold_val
                mean_ = np.mean(arr[i][arr[i] < a])
                new_array[i] = ((arr[i] - mean_) / abs(mean_)) * 100

                if np.isnan(new_array[i]).any() == True:
                    print('Warning:')
                    print('Cell %d: contains nan' % (i + 1))
                    print('      Mean of the sub-threshold for this cell: %s' % mean_)

        return new_array

    def save_pkl(self, pkl_path: str = None):
        """
        Alias method for saving current object to pickle file.

        :param pkl_path: (optional) provide path to save object to pickle file.
        """
        if pkl_path:
            parent = os.path.relpath(os.path.join(pkl_path, os.pardir))
            if os.path.exists(parent):
                print(f'saving new trial object to: {pkl_path}')
                self.pkl_path = pkl_path
            else:
                raise FileNotFoundError(f"Parent directory path: `{parent}` was not found for saving .pkl")


        with open(self.pkl_path, 'wb') as f:
            pickle.dump(self, f)
        print("\n\t -- data object saved to %s -- " % self.pkl_path)

    def save(self):
        """
        Alias method for saving current object as pickle file.
        """
        self.save_pkl()




