# this file tries an approach of building twophoton imaging trial on top of a general Trial class

import os
import time
import re
from typing import Optional, TypedDict

import pickle

import numpy as np
import pandas as pd

import tifffile as tf

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)
from packerlabimaging.main.classes import ImagingTrial
from packerlabimaging.processing.imagingMetadata import PrairieViewMetadata, ImagingMetadata
from packerlabimaging.utils.utils import SaveDownsampledTiff
from packerlabimaging.utils.classes import UnavailableOptionError, PaqInfo
from packerlabimaging.main.paq import PaqData
from packerlabimaging.processing import anndata as ad


class TwoPhotonImaging(ImagingTrial):
    """Two Photon Imaging Experiment Data Analysis Workflow."""

    def __init__(self, date: str = None, trialID: str = None, expID: str = None, tiff_path: str = None,
                 microscope: str = '',
                 expGroup: str = None, saveDir: str = None, PaqInfo: PaqInfo = None,
                 ImagingMetadata: ImagingMetadata = None,
                 comments: str = ''):

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

        ImagingTrial(date=date, trialID=trialID, expID=expID, dataPath=tiff_path, group=expGroup, comment=comments,
                     saveDir=saveDir)

        # ADD MODULES -
        self.Suite2p = None  #: Suite2p analysis sub-module  # todo consider adding wrapper method for attaching Suite2p to trial object (like might just need to refactor over from the experiment main file)

        print(f'\----- CREATING TwoPhotonImagingTrial for trial: \n\t{self.trialID}')

        # imaging metadata
        if 'Bruker' in microscope:
            self.imparams = PrairieViewMetadata(pv_xml_dir=self.tiff_path_dir)
        elif ImagingMetadata:
            self.imparams = ImagingMetadata
        else:
            Warning(
                f"NO imaging microscope parameters set. follow imagingMetadata to create a custom imagingMicroscopeMetadata class.")

        # temporal synchronization data from .Paq
        if PaqInfo:
            paq_path = PaqInfo['paq_path']
            if os.path.exists(paq_path):
                self._use_paq = True
            else:
                raise FileNotFoundError(f"paq_path does not exist: {paq_path}")

            frame_channel = PaqInfo['frame_channel'] if 'frame_channel' in [*PaqInfo] else KeyError(
                'No frame_channel specified for .paq processing')  # channel on Paq file to read for determining stims
            self.Paq = self._paqProcessingTwoPhotonImaging(paq_path=PaqInfo['paq_path'],
                                                           frame_channel=frame_channel)  #: Paq data submodule for trial

        # processing collect mean FOV Trace -- after collecting imaging params and Paq timing info
        self.meanFluImg, self.meanFovFluTrace = self.meanRawFluTrace()  #: mean image and mean FOV fluorescence trace

        self.data = None  #: anndata storage submodule
        self.dFF = None  #: dFF normalized traces of cells

        # SAVE Trial OBJECT
        self.save()

    def __str__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
        else:
            lastmod = "(unsaved pkl object)"
        return repr(f"ID: {self.t_series_name} (TwoPhotonImagingTrial experimental data object, last saved: {lastmod})")

    def __repr__(self):
        return repr(f"{self.t_series_name} (TwoPhotonImagingTrial experimental data object)")

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
    def paq_frame_clock(self):
        if hasattr(self.Paq, "frame_times"):
            return self.Paq.frame_times
        else:
            raise ValueError('Frame clock timings couldnt be retrieved from .Paq submodule.')

    def _paqProcessingTwoPhotonImaging(self, paq_path, frame_channel):
        """
        Add and further process paq data for current trial.
        :param paq_path: path to .paq file
        :param frame_channel: channel to use for measuring frame times from .paq data

        :return: PAQ data object
        """

        paq_data_obj = PaqData.import_paqdata(file_path=paq_path, plot=False)
        assert frame_channel in paq_data_obj.paq_channels, f"frame_channel argument: '{frame_channel}', not found in channels in .paq data."
        paq_data_obj.frame_times_channame = frame_channel
        paq_data_obj.frame_times = paq_data_obj.getPaqFrameTimes(frame_channel=frame_channel)
        paq_data_obj.sparse_paq_data = paq_data_obj.get_sparse_data(frame_times=paq_data_obj.frame_times)

        return paq_data_obj

    def stitch_s2p_reg_tiff(self):  ## TODO refactoring in new code from the Suite2p class script?
        assert self.Suite2p._s2pResultExists, UnavailableOptionError('stitch_s2p_reg_tiff')

        tif_path_save2 = self.saveDir + f'reg_tiff_{self.t_series_name}_r.tif'

        start = self.Suite2p.trial_frames[0] // 2000  # 2000 because that is the batch size for suite2p run
        end = self.Suite2p.trial_frames[1] // 2000 + 1

        if not os.path.exists(tif_path_save2):
            with tf.TiffWriter(tif_path_save2, bigtiff=True) as tif:
                with tf.TiffFile(self.Suite2p.reg_tiff_path, multifile=False) as input_tif:
                    print('cropping registered tiff')
                    data = input_tif.asarray()
                    print('shape of stitched tiff: ', data.shape)
                reg_tif_crop = data[self.Suite2p.trial_frames[0] - start * self.Suite2p.s2p_run_batch:
                                    self.Suite2p.trial_frames[1] - (
                                            self.Suite2p.trial_frames - start * self.Suite2p.s2p_run_batch)]
                print('saving cropped tiff ', reg_tif_crop.shape)
                tif.write(reg_tif_crop)

    def dfof(self):
        """(delta F)/F normalization of raw Suite2p data of trial."""
        assert hasattr(self,
                       'Suite2p'), 'no Suite2p module found. dfof function implemented to just normalize raw traces from Suite2p ROIs.'
        if self.Suite2p._s2pResultExists:
            dFF = self.normalize_dff(self.Suite2p.raw)
            return dFF

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