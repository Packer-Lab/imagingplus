import numpy as np
import pandas as pd
import pytest

from packerlabimaging import TwoPhotonImaging

import packerlabimaging as pli
from packerlabimaging.processing import anndata
from conftest import anndata_trial_data, existing_trialobj_alloptical_fixture


def test_AnnotatedData(existing_trialobj_alloptical_fixture):
    aotrial = existing_trialobj_alloptical_fixture
    aotrial.data = aotrial.create_anndata(imdata=aotrial.Suite2p,
                                          cells=aotrial.Suite2p,
                                          tmdata=aotrial.tmdata,
                                          imdata_type='suite2p raw - neuropil corrected')

test_AnnotatedData(existing_trialobj_alloptical_fixture)

@pytest.mark.skip
def test_convert_to_df(existing_anndata):
    'todo continue testing --- not working for everything.'
    adata = existing_anndata()
    df = adata.convert_to_df()
    print(df)
    return df



