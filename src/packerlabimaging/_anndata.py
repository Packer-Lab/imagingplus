import pandas as pd
import anndata

# trial level method
def convert_to_anndata(self):

    if self.Suite2p._s2pResultExists and self.paq_channels:
        # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
        # build dataframe for obs_meta from suite2p stat information
        obs_meta = pd.DataFrame(columns=['original_index', 'footprint', 'mrs', 'mrs0', 'compact', 'med', 'npix', 'radius',
                                         'aspect_ratio', 'npix_norm', 'skew', 'std'], index=range(len(self.Suite2p.stat)))
        for idx, __stat in enumerate(self.Suite2p.stat):
            for __column in obs_meta:
                obs_meta.loc[idx, __column] = __stat[__column]

        # SETUP THE VARIABLES ANNOTATIONS TO USE IN anndata
        # build dataframe for var annot's from paq file
        var_meta = pd.DataFrame(index=self.paq_channels, columns=range(self.n_frames))
        for fr_idx in range(self.n_frames):
            for index in [*self.sparse_paq_data]:
                var_meta.loc[index, fr_idx] = self.sparse_paq_data[index][fr_idx]

        adata = anndata.AnnData(X=self.Suite2p.raw, obs=obs_meta, var=var_meta.T)

        print(adata)
        return adata
    else:
        Warning('did not create anndata, only available if experiments were processed with suite2p and .paq file(s) provided for temporal synchronization')
