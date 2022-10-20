"""
Tests for nwb compatibility and interfacing functions.
"""
from packerlabimaging.utils.nwb import WriteImagingNWB, readImagingNWB


def test_WriteImagingNWB(new_imaging_nwb_fixture):
    subject = new_imaging_nwb_fixture[0]
    expobj = new_imaging_nwb_fixture[1]
    trialobj = new_imaging_nwb_fixture[2]
    expobj.experimenter = 'P. Shah'
    expobj.lab = 'packer'
    expobj.institution = 'oxford'
    trialobj.optical_channel_name = 'gcamp imaging channel'

    WriteImagingNWB(expobj=expobj, nwb_subject=subject, trialobj=trialobj, save=True, add_raw_tiff=True)


def test_readImagingNWB(existing_imaging_nwb_path_fixture):
    f = readImagingNWB(existing_imaging_nwb_path_fixture)
    print(f)



