# todo:
#  make tutorial for NWB conversion
#  add compatibility for SingleImage frames


import os.path
from datetime import datetime
from dateutil import tz

import numpy as np
from imagingplus.utils.images import ImportTiff

from imagingplus.main.subcore import TemporalData

from imagingplus.main.core import ImagingTrial
from pynwb.file import Subject
from pynwb import NWBFile, TimeSeries, NWBHDF5IO
from pynwb.image import ImageSeries
from pynwb.ophys import TwoPhotonSeries, OpticalChannel, ImageSegmentation, \
    Fluorescence, CorrectedImageStack, MotionCorrection, RoiResponseSeries

from imagingplus import TwoPhotonImaging, Experiment


def add_temporal_series(nwbfile: NWBFile, tmdata: TemporalData):
    """add temporal series data to the current nwb file using an imaging+ TemporalData input"""
    print(f"[NWB processing]: Adding temporal series to nwb file ...")
    for i, channel in enumerate(tmdata.channels):
        unit = tmdata.units[i] if tmdata.units else ''

        test_ts = TimeSeries(
            name=channel,
            data=tuple(tmdata.data[channel]),
            unit=unit,
            timestamps=np.linspace(tmdata.crop_offset_time,
                                   (tmdata.n_timepoints / tmdata.sampling_rate) + tmdata.crop_offset_time,
                                   tmdata.n_timepoints)
            # account for any cropped time (from the start of the recording)
        )

        nwbfile.add_acquisition(test_ts)

def save_nwb(nwbfile, path: str = None):
    """save current nwb file. specify path if given
    :param nwbfile: nwb file to save
    :param path: .nwb file path to save nwb filepath
    """
    print(f'\n\t ** Saving nwb file to: {path}')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with NWBHDF5IO(path, 'w') as io:
        io.write(nwbfile)

def newImagingNWB(expobj: Experiment, nwb_subject: Subject, trialobj: ImagingTrial, add_raw_tiff=False, save=True,
                 **kwargs):
    """
    Create a new NWB file from imaging+ data analysis pipeline.

    :param expobj:
    :param nwb_subject:
    :param trialobj:
    :param add_raw_tiff: set True to add the raw imaging tiffstack to the new NWB file
    :param save: set True to save the new NWB file
    :param kwargs:
        location: str; location of imaging in brain
        indicator: str; indicator used for imaging in brain

    """
    location = kwargs['location'] if 'location' in kwargs else ''
    indicator = kwargs['indicator'] if 'indicator' in kwargs else ''

    # initialize the NWB file
    nwbfile = NWBFile(
        session_description=expobj.comment,
        identifier=trialobj.expID,
        session_start_time=datetime.now(),
        experimenter=expobj.experimenter,
        lab=expobj.lab,
        institution=expobj.institution,
        subject=nwb_subject
    )

    device = nwbfile.create_device(
        name=trialobj.imparams.microscope
    )

    emission_lambda = trialobj.imparams.emission_lambda if hasattr(trialobj.imparams, 'emission_lambda') else 0.0
    excitation_lambda = trialobj.imparams.excitation_lambda if hasattr(trialobj.imparams,
                                                                       'excitation_lambda') else 0.0

    optical_channel = OpticalChannel(
        name=trialobj.optical_channel_name,
        description=trialobj.optical_channel_name,
        emission_lambda=emission_lambda

    )

    imaging_plane = nwbfile.create_imaging_plane(
        name=trialobj.imaging_plane_name,
        optical_channel=optical_channel,
        imaging_rate=trialobj.imparams.fps,
        description=trialobj.comment,
        excitation_lambda=excitation_lambda,
        location=location,
        indicator=indicator,
        device=device,
        grid_spacing=[trialobj.imparams.pix_sz_x, trialobj.imparams.pix_sz_y],  # spacing between pixels
        grid_spacing_unit=trialobj.imparams.PIXEL_SIZE_UNITS,
    )

    # add imaging data
    if add_raw_tiff:
        raw_tiff = ImportTiff(trialobj.data_path)
        image_series = TwoPhotonSeries(
            name='ImagingTrial',
            dimension=[trialobj.imparams.frame_x, trialobj.imparams.frame_y],
            imaging_plane=imaging_plane,
            rate=trialobj.imparams.fps,
            data=raw_tiff,
            unit=trialobj.imparams.IMAGING_SIGNAL_UNITS
        )

    else:
        image_series = TwoPhotonSeries(
            name='ImagingTrial',
            dimension=[trialobj.imparams.frame_x, trialobj.imparams.frame_y],
            external_file=[trialobj.data_path],
            imaging_plane=imaging_plane,
            rate=trialobj.imparams.fps,
        )

    print(f"[NWB processing]: Adding 2photon imaging series to nwb file ...")
    nwbfile.add_acquisition(image_series)

    # add Suite2p pre-processed data to nwb file
    if trialobj.Suite2p:
        ophys_module = nwbfile.create_processing_module(name='Pre-processed imaging data',
                                                     description='Pre-processed imaging data (ROIs and Fluorescence)')

        # create image segmentation object
        img_seg = ImageSegmentation()

        ps = img_seg.create_plane_segmentation(
            name='PlaneSegmentation',
            description=f'Image segmentation output from {trialobj.imaging_plane_name}',
            imaging_plane=imaging_plane,
        )
        ophys_module.add(img_seg)

        # ADD ROIs from Suite2p ROI segmentation
        for i, roi in enumerate(list(trialobj.Suite2p.cells_array)):
            pixel_mask = []
            for iy in trialobj.Suite2p.multidim_data['ypix'][i]:
                for ix in trialobj.Suite2p.multidim_data['xpix'][i]:
                    pixel_mask.append((iy, ix, 1))
            ps.add_roi(pixel_mask=pixel_mask)

        # add pre-processed data to nwb file
        if trialobj.Suite2p.imdata is not None:
            rt_region = ps.create_roi_table_region(
                region=list(np.arange(0, trialobj.Suite2p.n_cells)),
                description=f'Segmented ROIs'
            )

            roi_resp_series = RoiResponseSeries(
                name='ROIs x extracted fluorescence signal',
                data=trialobj.Suite2p.imdata,
                rois=rt_region,
                rate=trialobj.imparams.fps,
                unit=trialobj.imparams.IMAGING_SIGNAL_UNITS
            )

            fl = Fluorescence(roi_response_series=roi_resp_series)
            ophys_module.add(fl)

    # add temporal data
    if trialobj.tmdata: add_temporal_series(nwbfile=nwbfile, tmdata=trialobj.tmdata)

    # save newly generated nwb file
    _nwb_save_path = trialobj.pkl_path[:-4] + '.nwb'

    if save:
        save_nwb(nwbfile=nwbfile, path=_nwb_save_path)
        # print(f'\n\t ** Saving nwb file to: {_nwb_save_path}')
        # os.makedirs(os.path.dirname(_nwb_save_path), exist_ok=True)
        # with NWBHDF5IO(_nwb_save_path, 'w') as io:
        #     io.write(nwbfile)

    return nwbfile





def readImagingNWB(filename):
    """
    Reading in imaging+ generated NWB files.

    """
    print(f"[NWB processing]: Reading NWB file from:\n\t {filename}")

    with NWBHDF5IO(filename, 'r') as io:
        nwbfile: NWBFile = io.read()

    # if filename != nwbfile.nwb_path: nwbfile.nwb_path = filename

    return nwbfile
