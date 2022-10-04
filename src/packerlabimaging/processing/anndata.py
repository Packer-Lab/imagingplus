import anndata as ad
from typing import Optional, Literal

import numpy as np
import pandas as pd


class AnnotatedData(ad.AnnData):
    """Creates annotated cellsdata (see anndata library for more information on AnnotatedData) object based around the Ca2+ matrix of the imaging trial."""

    def __init__(self, X, obs, var=None, data_label=None, **kwargs):
        adata_dict = {'X': X, 'obs': obs, 'var': var}
        for key in [*kwargs]:
            adata_dict[key] = kwargs[key]

        ad.AnnData.__init__(self, **adata_dict)
        self.data_label = data_label if data_label else None

        print(f"Created AnnData object: \n\t{self.__repr__()}")

    def __str__(self):
        "extensive information about the AnnotatedData cellsdata structure"
        if self.filename:
            backed_at = f" backed at {str(self.filename)!r}"
        else:
            backed_at = ""

        descr = f"Annotated Data of n_obs × n_vars = {self.n_obs} × {self.n_vars} {backed_at}"
        descr += f"\navailable attributes: "

        descr += f"\n\t.X (primary datamatrix) of .data_label: \n\t\t|- {str(self.data_label)}" if self.data_label else f"\n\t.X (primary datamatrix)"
        descr += f"\n\t.obs (obs metadata): \n\t\t|- {str(list(self.obs.keys()))[1:-1]}"
        descr += f"\n\t.var (vars metadata): \n\t\t|- {str(list(self.var.keys()))[1:-1]}"
        for attr in [
            ".uns",
            ".obsm",
            ".varm",
            ".layers",
            ".obsp",
            ".varp",
        ]:
            keys = getattr(self, attr[1:]).keys()
            if len(keys) > 0:
                descr += f"\n\t{attr}: \n\t\t|- {str(list(keys))[1:-1]}"
        return descr

    def _gen_repr(self, n_obs, n_vars) -> str:  # overriding base method from AnnData
        """overrides the default anndata _gen_repr_() method for imaging cellsdata usage.

        :param n_obs: number of observations in anndata table
        :param n_vars: number of variables in anndata table
        :return:
        """

        return f"Annotated Data of n_obs (# ROIs) × n_vars (# Frames) = {n_obs} × {n_vars}"

    def add_obs(self, obs_name: str, values: list):
        """adds values to the observations of an anndata object, under the key obs_name

        :param obs_name: name of new observation field to add to anndata table
        :param values: list of data values to add under new observation field
        """
        assert len(values) == self.obs.shape[0], f"# of values to add doesn't match # of observations in anndata array"
        self.obs[obs_name] = values

    def del_obs(self, obs_name: str):
        """removes a key from observations from an anndata object, of the key obs_name

        :param obs_name: name of observation to remove
        """
        _ = self.obs.pop(obs_name)

    def add_var(self, var_name: str, values: list):
        """adds values to the variables of an anndata object, under the key var_name

        :param var_name: name of new variable field to add to anndata table
        :param values: list of data values to add under new variable field
        """
        assert len(values) == self.var.shape[0], f"# of values to add doesn't match # of observations in anndata array"
        self.var[var_name] = values

    def del_var(self, var_name: str):
        """removes a key from variables from an anndata object, of the key var_name

        :param var_name: name of variable to remove
        """
        _ = self.var.pop(var_name)

    def extend_anndata(self, additional_adata: ad.AnnData, axis: Literal[0, 1] = 0):
        """
        Extend with an additional anndata object. Specify axis to extend. Ensure that new anndata object matches number of variables or observations depending on
        which axis is being extended.

        :param additional_adata: an anndata object with matching observations or variables (depending on which axis is being extended).
        :param axis: axis on which to extend
        """
        adata = ad.concat([self, additional_adata], axis=axis)
        return adata

    def convert_to_df(self) -> pd.DataFrame:
        """
        Convert anndata object into a long-form pandas dataframe. primary purpose is to allow access to pandas and seaborn functionality more directly.

        - overall seems to be working well. just need to test with a dataset with >1 obs and var keys(), and to test with the whole larger dataset.
        :return: long-form pandas dataframe

        """

        print(f"\n\- converting anndata cellsdata matrix to long-form pandas dataframe ... [in progress]")

        cols = [self.obs_keys()[0], self.var_keys()[0]]
        cols.extend(self.obs_keys()[1:])
        cols.extend(self.var_keys()[1:])
        cols.extend([self.data_label]) if self.data_label is not None or not '' else cols.extend('data_values')

        df = pd.DataFrame(columns=cols)
        index = 0
        for idi in range(self.n_obs):
            for idj in range(self.n_vars):
                dict_ = {}
                for col in self.obs_keys():
                    dict_[str(col)] = self.obs[col][idi]

                for col in self.var_keys():
                    dict_[str(col)] = self.var[col][idj]

                dict_[cols[-1]] = self.X[idi, idj]
                df = pd.concat([df, pd.DataFrame(dict_, index=[index])])
                index += 1

        print(f"\n|- converting anndata cellsdata matrix to long-form pandas dataframe ... [finished]")

        return df

    # @classmethod
    # def create_anndata(cls, trial: ImagingTrial):
    #     """
    #     Alternative constructor to create anndata object using ImagingTrial as input.
    #     Creates annotated cellsdata (see anndata library for more information on AnnotatedData) object based around the Ca2+ matrix of the imaging trial.
    #
    #     """
    #     if trial.cells and trial.tmdata and trial.imdata:
    #         # SETUP THE OBSERVATIONS (CELLS) ANNOTATIONS TO USE IN anndata
    #         # build dataframe for obs_meta from suite2p stat information
    #         obs_meta = trial.cells.cellsdata
    #
    #         var_meta = trial.tmdata.cellsdata
    #
    #         assert obs_meta.shape[0] == trial.imdata.cellsdata.shape[1], '.cells.cellsdata shape does not match .imdata.cellsdata shape that are being set together.'
    #         if var_meta.shape[0] == trial.imdata.cellsdata.shape[0]:
    #             var_meta = var_meta.T
    #         elif var_meta.shape[1] != trial.imdata.cellsdata.shape[0]:
    #             raise ValueError('.tmdata.cellsdata shape does not match .imdata.cellsdata shape that are being set together.')
    #
    #
    #         print(f"\n\----- CREATING annotated cellsdata object using AnnData:")
    #         adata = cls(X=trial.imdata.cellsdata, obs=obs_meta, var=var_meta.T)
    #
    #         print(f"\n{adata}")
    #         return adata
    #
    #     else:
    #         Warning(
    #             'did not create anndata. anndata creation only available if experiments were processed with suite2p and .Paq file(s) provided for temporal synchronization')
