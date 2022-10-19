from datetime import datetime
from dateutil.tz import tzlocal

import numpy as np
from packerlabimaging.main.core import ImagingTrial
from pynwb.file import Subject
from pynwb import NWBFile, TimeSeries, NWBHDF5IO
from pynwb.image import ImageSeries
from pynwb.ophys import TwoPhotonSeries, OpticalChannel, ImageSegmentation, \
    Fluorescence, CorrectedImageStack, MotionCorrection, RoiResponseSeries

import matplotlib.pyplot as plt

from packerlabimaging import TwoPhotonImaging, Experiment

#### FROM TUTORIAL
# when we create the imaging trial, we also add modules for parsing simultaneously temporal data (.paq files) and imaging microscope's parameters for this imaging trial
from packerlabimaging.processing.paq import PaqData
from packerlabimaging.processing.imagingMetadata import PrairieViewMetadata

# create trial obj for each trial of experiment

trials_list_spont = ['t-005', 't-006']
for idx, trial in enumerate(trials_list_spont):
    date = '2020-12-19'
    prep = 'RL113'

    # define paths to retrieve data from disc
    paqs_loc = f'/home/pshah/mnt/qnap/Data/2020-12-19/{date}_{prep}_{trial[-3:]}.paq'  # path to the .paq files for the selected trials
    dataPath = f'/home/pshah/mnt/qnap/Data/2020-12-19/{date}_{trial}/{date}_{trial}_Cycle00001_Ch3.tif'

    # initialize microscope meta-data for current imaging trial using pre-built module for Bruker-PrairieView xml parser
    imparams = PrairieViewMetadata(pv_xml_dir=os.path.dirname(dataPath), microscope='Bruker 2pPlus')

    # initialize temporal data using pre-built module for .paq parsing
    tmdata = PaqData.paqProcessingTwoPhotonImaging(paq_path=paqs_loc, frame_channel='frame_clock', plot=False)

    # feed in information into `TwoPhotonImaging` to create a TwoPhotonImaging object
    trialobj = TwoPhotonImaging(date=date, trialID=trial, expID=prep, imparams=imparams, tmdata=tmdata,
                                    saveDir=f'/mnt/qnap_share/Data/packerlabimaging-example/packerlabimaging-test-analysis/',
                                    dataPath=dataPath, expGroup="awake spont. 2p imaging + LFP")

    # add each Trial to the overall Experiment using the trialobj

#### END

# example:
subject = Subject(age='P60D',
                  subject_id='RL113',
                  genotype='CamkIIa-GCaMP6s/Niell'
                  )



class WriteImagingNWB(NWBFile):
    """
    Container for creating a new NWB file.
    Works by extending the base NWBFile container for Imaging+ related requirements and usages.

    """

    def __init__(self, expobj: Experiment, nwb_subject: Subject, trialobj: ImagingTrial, **kwargs):
        # initialize the NWB file
        super().__init__(
            session_description=expobj.comment,
            identifier=trialobj.expID,
            session_start_time=datetime.now(tzlocal()),
            experimenter=expobj.experimenter,
            lab=expobj.lab,
            institution=expobj.institution,
            subject=nwb_subject
        )

        device = self.create_device(
            name=imparams.microscope
        )

        optical_channel = OpticalChannel(
            name=trialobj.optical_channel_name,
        )

        imaging_plane = self.create_imaging_plane(
            name=trialobj.imaging_plane_name,
            optical_channel=optical_channel,
            imaging_rate=trialobj.imparams.fps,
            description=trialobj.comment,
            device=device,
            grid_spacing=[trialobj.imparams.pix_sz_x, trialobj.imparams.pix_sz_y],  # spacing between pixels
            grid_spacing_unit=trialobj.imparams.PIXEL_SIZE_UNITS,
            **kwargs
        )

        image_series = TwoPhotonSeries(
            name='ImagingTrial',
            dimension=[trialobj.imparams.frame_x, trialobj.imparams.frame_y],
            external_file=[trialobj.data_path],
            imaging_plane=imaging_plane,
            rate=trialobj.imparams.fps
        )

        # todo add temporal data from paq processing

        self.add_acquisition(image_series)

        # save newly generated nwb file
        self.__nwb_save_path = trialobj.pkl_path[:-4] + '.nwb'
        # self.save()

    @property
    def nwb_path(self):
        """save path for current nwb file"""
        return self.__nwb_save_path

    @nwb_path.setter
    def nwb_path(self, path):
        """set new save path for current nwb file"""
        self.__nwb_save_path = path

    def save(self):
        print(f'\n\t ** Saving nwb file to: {self.nwb_path}')
        with NWBHDF5IO(self.nwb_path, 'w') as io:
            io.write(self)


def readImagingNWB(filename):
    """
    Reading in imaging+ generated NWB files.

    """
    with NWBHDF5IO(filename, 'r') as io:
        nwbfile = io.read()

    return nwbfile



