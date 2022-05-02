# this file tries an approach of building twophoton imaging trial on top of a general Trial class

import os
import time

import numpy as np

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)
from packerlabimaging.main.classes import ImagingTrial, CellAnnotations
from packerlabimaging.main.paq import PaqData
from packerlabimaging.processing.imagingMetadata import ImagingMetadata


class TwoPhotonImaging(ImagingTrial):
    """Two Photon Imaging Experiment Data Analysis Workflow."""

    def __init__(self, date: str = None, trialID: str = None, expID: str = None, dataPath: str = None,
                 expGroup: str = None, saveDir: str = None, tmdata: PaqData = None,
                 imparams: ImagingMetadata = None, cells: CellAnnotations = None, comment: str = ''):

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

        print(f'\----- CREATING TwoPhotonImagingTrial for trial: {trialID}')

        super().__init__(date=date, trialID=trialID, expID=expID, dataPath=dataPath, expGroup=expGroup, comment=comment,
                         saveDir=saveDir, imparams=imparams, cells=cells, tmdata=tmdata)

        # processing collect mean FOV Trace -- after collecting imaging params and Paq timing info
        self.meanFluImg, self.meanFovFluTrace = self.meanRawFluTrace()  #: mean image and mean FOV fluorescence trace

        # SAVE two photon trial object
        self.save()

    def __str__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
        else:
            lastmod = "(unsaved pkl object)"
        return repr(f"ID: {self.t_series_name} (TwoPhotonImagingTrial experimental data object, last saved: {lastmod})")

    def __repr__(self):
        return repr(f"{self.t_series_name} (TwoPhotonImagingTrial experimental data object)")

    # def _getImagingParameters(self, metadata: Optional[dict] = None, microscope: Optional[str] = 'Bruker'):
    #     """retrieves imaging metadata parameters. If using Bruker microscope and PrairieView, then _prairieview module is used to collect this data.
    #
    #     :param microscope: name of the microscope, currently the only supported microscope for parsing metadata directly is Bruker/PrairieView imaging setup.
    #     """
    #     if 'Bruker' in microscope:
    #         return PrairieViewMetadata(tiff_path_dir=self.tiff_path_dir)
    #     else:
    #         try:
    #             return ImagingMetadata(**metadata)
    #         except TypeError:
    #             Exception('required key not present in metadata')
    #
    # @property
    # def frame_clock(self):
    #     if hasattr(self.tmdata, "frame_times"):
    #         return self.tmdata.frame_times
    #     else:
    #         raise ValueError('Frame clock timings couldnt be retrieved from .tmdata submodule.')

    # def stitch_s2p_reg_tiff(self):  ## TODO refactoring in new code from the Suite2p class script?
    #     assert self.Suite2p._s2pResultExists, UnavailableOptionError('stitch_s2p_reg_tiff')
    #
    #     tif_path_save2 = self.saveDir + f'reg_tiff_{self.t_series_name}_r.tif'
    #
    #     start = self.Suite2p.trial_frames[0] // 2000  # 2000 because that is the batch size for suite2p run
    #     end = self.Suite2p.trial_frames[1] // 2000 + 1
    #
    #     if not os.path.exists(tif_path_save2):
    #         with tf.TiffWriter(tif_path_save2, bigtiff=True) as tif:
    #             with tf.TiffFile(self.Suite2p.reg_tiff_path, multifile=False) as input_tif:
    #                 print('cropping registered tiff')
    #                 data = input_tif.asarray()
    #                 print('shape of stitched tiff: ', data.shape)
    #             reg_tif_crop = data[self.Suite2p.trial_frames[0] - start * self.Suite2p.s2p_run_batch:
    #                                 self.Suite2p.trial_frames[1] - (
    #                                         self.Suite2p.trial_frames - start * self.Suite2p.s2p_run_batch)]
    #             print('saving cropped tiff ', reg_tif_crop.shape)
    #             tif.write(reg_tif_crop)

    # def dfof(self):
    #     """(delta F)/F normalization of raw Suite2p data of trial."""
    #     assert hasattr(self,
    #                    'Suite2p'), 'no Suite2p module found. dfof function implemented to just normalize raw traces from Suite2p ROIs.'
    #     if self.Suite2p._s2pResultExists:
    #         dFF = self.normalize_dff(self.Suite2p.raw)
    #         return dFF

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
