"NOTE: THESE FUNCTIONS REFACTORED TO _utils.Utils !!! - jan 18 2022"

import anndata as ad
from typing import Optional

class AnnotatedData(ad.AnnData):
    """Creates annotated data (see anndata library for more information on AnnotatedData) object based around the Ca2+ matrix of the imaging trial."""

    def __init__(self, X, obs, var: Optional=None, data_label=None, **kwargs):
        __adata_dict = {'X': X, 'obs': obs, 'var': var}
        for key in [*kwargs]:
            __adata_dict[key] = kwargs[key]

        super().__init__(**__adata_dict)
        self.data_label = data_label if data_label else None


    def _gen_repr(self, n_obs, n_vars) -> str:
        "modification of the default anndata _repr_"
        if self.filename:
            backed_at = f" backed at {str(self.filename)!r}"
        else:
            backed_at = ""

        descr = f"Annotated Data of n_obs (# ROIs) × n_vars (# Frames) = {n_obs} × {n_vars} {backed_at}"
        descr += f"\navailable attributes: "

        descr += f"\n\t.X (primary datamatrix, with data_label): {str(self.data_label)}" if self.data_label else f"\n\t.X (primary datamatrix)"
        descr += f"\n\t.obs (ROIs metadata) keys: {str(list(self.obs.keys()))[1:-1]}"
        descr += f"\n\t.var (frames metadata) keys: {str(list(self.var.keys()))[1:-1]}"
        for attr in [
            # "obs",
            # "var",
            ".uns",
            ".obsm",
            ".varm",
            ".layers",
            ".obsp",
            ".varp",
        ]:
            keys = getattr(self, attr[1:]).keys()
            if len(keys) > 0:
                descr += f"\n\t{attr} keys: {str(list(keys))[1:-1]}"
        return descr


    def add_observation(self, obs_name: str, values: list):
        """adds values to the observations of an anndata object, under the key obs_name"""
        assert len(values) == self.obs.shape[0], f"# of values to add doesn't match # of observations in anndata"
        self.obs[obs_name] = values

    def add_variables(self, var_name: str, values: list):
        """adds values to the variables of an anndata object, under the key var_name"""
        assert len(values) == self.var.shape[0], f"# of values to add doesn't match # of observations in anndata"
        self.var[var_name] = values

    def extend_anndata(self, additional_adata: ad.AnnData, axis: int = 0):
        """
        :param adata_obj: an anndata object of dimensions n obs x m var
        :param additional_adata: an anndata object of dimensions n obs x # var or, # obs x m var (depending on which axis to extend)
        """
        adata = ad.concat([self, additional_adata], axis=axis)
        return adata


def extend_anndata(adata_obj: ad.AnnData, additional_adata: ad.AnnData, axis: int = 0):
    """
    :param adata_obj: an anndata object of dimensions n obs x m var
    :param additional_adata: an anndata object of dimensions n obs x # var or, # obs x m var (depending on which axis to extend)
    """
    adata = ad.concat([adata_obj, additional_adata], axis=axis)
    return adata

def add_observation(adata_obj: ad.AnnData, obs_name: str, values: list):
    """adds values to the observations of an anndata object, under the key obs_name"""
    assert len(values) == adata_obj.obs.shape[0], f"# of values to add doesn't match # of observations in anndata"
    adata_obj.obs[obs_name] = values
    return adata_obj

def add_variables(adata_obj: ad.AnnData, var_name: str, values: list):
    """adds values to the variables of an anndata object, under the key var_name"""
    assert len(values) == adata_obj.var.shape[0], f"# of values to add doesn't match # of observations in anndata"
    adata_obj.var[var_name] = values
    return adata_obj
