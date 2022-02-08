import numpy as np
import pandas as pd
from packerlabimaging.processing.anndata import AnnotatedData


# %%
# TEST
n_obs, n_vars = 100, 10000
X = np.random.random((n_obs, n_vars))
X2 = np.random.random((n_obs, n_vars))



df2 = pd.DataFrame(X2, columns=list('ABCDEFGHIJ')*1000, index=np.arange(n_obs, dtype=int).astype(str))

obs_meta = pd.DataFrame({
        'time_yr': np.random.choice([0, 2, 4, 8], n_obs),
        'subject_id': np.random.choice(['subject 1', 'subject 2', 'subject 4', 'subject 8'], n_obs),
        'instrument_type': np.random.choice(['type a', 'type b'], n_obs),
        'site': np.random.choice(['site x', 'site y'], n_obs),
    },
    index=np.arange(n_obs, dtype=int).astype(str),    # these are the same IDs of observations as above!
)

var_meta = pd.DataFrame({
    'var1': np.random.choice(['a', 'b', 'c'], n_vars)
}, index=df2.columns)

adata_ = AnnotatedData(X=df2, obs=obs_meta, var=var_meta, data_label='df2')

adata_.add_observation('new_obs', values=['a']*n_obs)

print(adata_)