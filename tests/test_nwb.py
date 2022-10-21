"""
Tests for nwb compatibility and interfacing functions.
"""
from imagingplus.utils.nwb import newImagingNWB, readImagingNWB


def test_WriteImagingNWB(new_imaging_nwb_fixture):
    subject = new_imaging_nwb_fixture[0]
    expobj = new_imaging_nwb_fixture[1]
    trialobj = new_imaging_nwb_fixture[2]
    expobj.experimenter = 'P. Shah'
    expobj.lab = 'packer'
    expobj.institution = 'oxford'
    trialobj.optical_channel_name = 'gcamp imaging channel'

    inwb = newImagingNWB(expobj=expobj, nwb_subject=subject, trialobj=trialobj, save=True, add_raw_tiff=True)
    print(inwb)


def test_readImagingNWB(existing_imaging_nwb_path_fixture):
    f = readImagingNWB(existing_imaging_nwb_path_fixture)
    print(f)



