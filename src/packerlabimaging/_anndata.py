import numpy as np
import pandas as pd
import anndata

# trial level class
def create_anndata(trialobj):
    """
    Creates annotated data (see anndata library) object based around the Ca2+ matrix of the imaging trial.

    :param trialobj: TwoPhotonImagingTrial (or derivative) object containing data for creating anndata object.

    """

    if trialobj.Suite2p._s2pResultExists and trialobj.paq_channels:
        # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
        # build dataframe for obs_meta from suite2p stat information
        obs_meta = pd.DataFrame(columns=['original_index', 'footprint', 'mrs', 'mrs0', 'compact', 'med', 'npix', 'radius',
                                         'aspect_ratio', 'npix_norm', 'skew', 'std'], index=range(len(trialobj.Suite2p.stat)))
        for idx, __stat in enumerate(trialobj.Suite2p.stat):
            for __column in obs_meta:
                obs_meta.loc[idx, __column] = __stat[__column]

        # build numpy array for multidimensional obs metadata
        obs_m = {'ypix': [],
                 'xpix': []}
        for col in [*obs_m]:
            for idx, __stat in enumerate(trialobj.Suite2p.stat):
                obs_m[col].append(__stat[col])
            obs_m[col] = np.asarray(obs_m[col])

        # SETUP THE VARIABLES ANNOTATIONS TO USE IN anndata
        # build dataframe for var annot's from paq file
        var_meta = pd.DataFrame(index=trialobj.paq_channels, columns=range(trialobj.n_frames))
        for fr_idx in range(trialobj.n_frames):
            for index in [*trialobj.sparse_paq_data]:
                var_meta.loc[index, fr_idx] = trialobj.sparse_paq_data[index][fr_idx]

        # BUILD LAYERS TO ADD TO anndata OBJECT
        layers = {'dFF': trialobj.dFF
                  }

        print(f"\----- CREATING annotated data object using AnnData:")
        adata = anndata.AnnData(X=trialobj.Suite2p.raw, obs=obs_meta, var=var_meta.T, obsm=obs_m,
                                layers=layers)

        print(f"\t{adata}")
        return adata
    else:
        Warning('did not create anndata. anndata creation only available if experiments were processed with suite2p and .paq file(s) provided for temporal synchronization')




def extend_anndata(adata_obj, additional_adata, axis: int = 0):
    adata = anndata.concat([adata_obj, additional_adata], axis=axis)
    return adata