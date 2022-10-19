"""
Tests for nwb compatibility and interfacing functions.
"""
from packerlabimaging.utils.nwb import WriteImagingNWB


def test_WriteImagingNWB(new_imaging_nwb_fixture):
    subject = new_imaging_nwb_fixture[0]
    expobj = new_imaging_nwb_fixture[1]
    trialobj = new_imaging_nwb_fixture[2]

    WriteImagingNWB(expobj=expobj, nwb_subject=subject, trialobj=trialobj, save=True)



