from packerlabimaging.processing.imagingMetadata import PrairieViewMetadata


def test_PrairieViewMetadata():
    testdir = '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_t-002'
    pvdata = PrairieViewMetadata(pv_xml_dir=testdir)

test_PrairieViewMetadata()

