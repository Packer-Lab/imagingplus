from packerlabimaging.utils.imagingMetadata import PrairieViewMetadata


def test_PrairieViewMetadata():
    # TODO run tests for suite2p trial creation without predone results path
    testdir = '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_t-002'
    pvdata = PrairieViewMetadata(tiff_path_dir=testdir)

test_PrairieViewMetadata()

