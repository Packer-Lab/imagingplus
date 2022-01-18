import numpy as np
import pandas as pd
import anndata


def extend_anndata(adata_obj: anndata.AnnData, additional_adata: anndata.AnnData, axis: int = 0):
    """
    :param adata_obj: an anndata object of dimensions n obs x m var
    :param additional_adata: an anndata object of dimensions n obs x # var or, # obs x m var (depending on which axis to extend)
    """
    adata = anndata.concat([adata_obj, additional_adata], axis=axis)
    return adata

def add_observation(adata_obj: anndata.AnnData, obs_name: str, values: list):
    """adds values to the observations of an anndata object, under the key obs_name"""
    assert len(values) == adata_obj.obs.shape[0], f"# of values to add doesn't match # of observations in anndata"
    adata_obj.obs[obs_name] = values
    return adata_obj

def add_variables(adata_obj: anndata.AnnData, var_name: str, values: list):
    """adds values to the variables of an anndata object, under the key var_name"""
    assert len(values) == adata_obj.var.shape[0], f"# of values to add doesn't match # of observations in anndata"
    adata_obj.var[var_name] = values
    return adata_obj
