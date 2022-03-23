import numpy as np
import pandas as pd
from packerlabimaging import TwoPhotonImagingTrial

import packerlabimaging as pli
from packerlabimaging.processing import anndata
from conftest import anndata_trial_data


def test_AnnotatedData():
    X, var, obs = anndata_trial_data()
    adata = anndata.AnnotatedData(X=X, obs=obs, var=var, data_label='test_data')
    print(adata)
    return adata

adata = test_AnnotatedData()


def test_convert_to_df(existing_anndata):
    'todo continue testing --- not working for everything.'
    adata = existing_anndata()
    df = adata.convert_to_df()
    print(df)
    return df

df = adata.convert_to_df()

# %%
import seaborn as sns
import matplotlib.pyplot as plt

# %%
sns.scatterplot(data=df, x="group", y="test_data")
plt.show()

def make_seaborn_plots(df):
    """todo - make seaborn plot from long-form dataframe to use as example in Anndata tutorial."""



    pass





