# this file tries an approach of building twophoton imaging trial on top of a general Trial class

import os
import time

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)
from imagingplus.main.core import ImagingTrial
from imagingplus.main.subcore import CellAnnotations, ImagingData
from imagingplus.processing.paq import PaqData
from imagingplus.processing.imagingMetadata import ImagingMetadata


class TwoPhotonImaging(ImagingTrial):
    """Two Photon Imaging Experiment Data Analysis Workflow."""

    description = 'Imaging Channel'

    def __init__(self, date: str = None, trialID: str = None, expID: str = None, dataPath: str = None,
                 expGroup: str = None, saveDir: str = None, tmdata: PaqData = None, imdata: ImagingData = None,
                 imparams: ImagingMetadata = None, cells: CellAnnotations = None, comment: str = ''):

        """
        TODO update function docstring for approp args
        :param metainfo: TypedDict containing meta-information field needed for this experiment. Please see TwoPhotonImagingMetainfo for type hints on accepted keys.
        :param paq_options: TypedDict containing meta-information about .paq file associated with current trial
        :param analysis_save_path: path of where to save the experiment analysis object
        :param microscope: name of microscope used to record imaging (options: "Bruker" (default), "other")
        :param imagingMicroscopeMetadata: provide ImagingMetadata object (see ImagingMetadata class).
        :param suite2p_experiment_obj: provide Suite2p Experiment Object as variable in order to process Suite2p cellsdata for current trial
        :param total_frames_stitched: provide frame number on which current trial starts in Suite2p Experiment Object
        """

        print(f'\----- CREATING TwoPhotonImagingTrial for trial: {trialID}')

        super().__init__(date=date, trialID=trialID, expID=expID, dataPath=dataPath, expGroup=expGroup, comment=comment,
                         saveDir=saveDir, imparams=imparams, cells=cells, tmdata=tmdata)

        # SAVE two photon trial object
        self.save()

    def __str__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
        else:
            lastmod = "(unsaved pkl object)"
        return repr(f"ID: {self.t_series_name} (TwoPhotonImagingTrial experimental object, last saved: {lastmod})")

    def __repr__(self):
        return repr(f"{self.t_series_name} (TwoPhotonImagingTrial experimental object)")

    # def _getImagingParameters(self, metadata: Optional[dict] = None, microscope: Optional[str] = 'Bruker'):
    #     """retrieves imaging metadata parameters. If using Bruker microscope and PrairieView, then _prairieview module is used to collect this cellsdata.
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
    #                 cellsdata = input_tif.asarray()
    #                 print('shape of stitched tiff: ', cellsdata.shape)
    #             reg_tif_crop = cellsdata[self.Suite2p.trial_frames[0] - start * self.Suite2p.s2p_run_batch:
    #                                 self.Suite2p.trial_frames[1] - (
    #                                         self.Suite2p.trial_frames - start * self.Suite2p.s2p_run_batch)]
    #             print('saving cropped tiff ', reg_tif_crop.shape)
    #             tif.write(reg_tif_crop)

    # def dfof(self):
    #     """(delta F)/F normalization of raw Suite2p cellsdata of trial."""
    #     assert hasattr(self,
    #                    'Suite2p'), 'no Suite2p module found. dfof function implemented to just normalize raw traces from Suite2p ROIs.'
    #     if self.Suite2p._s2pResultExists:
    #         dFF = self.normalize_dff(self.Suite2p.raw)
    #         return dFF

