import pytest


def test_SaveDownsampledTiff(tiff_path_fixture):
    from packerlabimaging.utils.utils import ImportTiff
    start, end = 0, 4000
    stack = ImportTiff(tiff_path_fixture, frames=(start,end))
    from packerlabimaging.utils.utils import SaveDownsampledTiff
    SaveDownsampledTiff(stack=stack, save_as=f'/home/pshah/mnt/qnap/Analysis/2020-03-03/2020-03-03_t-001_{start}-{end}fr_pre4ap_downsampled.tif',
                        group_by=8)