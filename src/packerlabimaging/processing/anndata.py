import anndata as ad
from typing import Optional, Literal

import pandas as pd


class AnnotatedData(ad.AnnData):
    """Creates annotated data (see anndata library for more information on AnnotatedData) object based around the Ca2+ matrix of the imaging trial."""

    def __init__(self, X, obs, var=None, data_label=None, **kwargs):
        adata_dict = {'X': X, 'obs': obs, 'var': var}
        for key in [*kwargs]:
            adata_dict[key] = kwargs[key]

        ad.AnnData.__init__(self, **adata_dict)
        self.data_label = data_label if data_label else None

        print(f"Created AnnData object: \n\t{self.__repr__()}")

    def __str__(self):
        "extensive information about the AnnotatedData data structure"
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
                descr += f"\n\t{attr}: \n\t\t|- {str(list(keys))[1:-1]}"
        return descr


    def _gen_repr(self, n_obs, n_vars) -> str:  # overriding base method from AnnData
        """overrides the default anndata _gen_repr_() method for imaging data usage."""

        return f"Annotated Data of n_obs (# ROIs) × n_vars (# Frames) = {n_obs} × {n_vars}"

        # if self.filename:
        #     backed_at = f" backed at {str(self.filename)!r}"
        # else:
        #     backed_at = ""
        #
        # descr = f"Annotated Data of n_obs (# ROIs) × n_vars (# Frames) = {n_obs} × {n_vars} {backed_at}"
        # descr += f"\navailable attributes: "
        #
        # descr += f"\n\t.X (primary datamatrix, with .data_label): \n\t\t|-- {str(self.data_label)}" if self.data_label else f"\n\t.X (primary datamatrix)"
        # descr += f"\n\t.obs (ROIs metadata) keys: \n\t\t|-- {str(list(self.obs.keys()))[1:-1]}"
        # descr += f"\n\t.var (frames metadata) keys: \n\t\t|-- {str(list(self.var.keys()))[1:-1]}"
        # for attr in [
        #     # "obs",
        #     # "var",
        #     ".uns",
        #     ".obsm",
        #     ".varm",
        #     ".layers",
        #     ".obsp",
        #     ".varp",
        # ]:
        #     keys = getattr(self, attr[1:]).keys()
        #     if len(keys) > 0:
        #         descr += f"\n\t{attr} keys: \n\t\t|-- {str(list(keys))[1:-1]}"
        # return descr


    def add_obs(self, obs_name: str, values: list):
        """adds values to the observations of an anndata object, under the key obs_name"""
        assert len(values) == self.obs.shape[0], f"# of values to add doesn't match # of observations in anndata array"
        self.obs[obs_name] = values

    def del_obs(self, obs_name: str): # TODO
        """removes a key from observations from an anndata object, of the key obs_name"""
        _ = self.obs.pop(obs_name)


    def add_var(self, var_name: str, values: list):
        """adds values to the variables of an anndata object, under the key var_name"""
        assert len(values) == self.var.shape[0], f"# of values to add doesn't match # of observations in anndata array"
        self.var[var_name] = values

    def del_var(self, obs_name: str): # TODO
        """removes a key from variables from an anndata object, of the key var_name"""
        _ = self.var.pop(obs_name)


    def extend_anndata(self, additional_adata: ad.AnnData, axis: Literal[0,1] = 0):
        """
        :param additional_adata: an anndata object of dimensions n obs x # var or, # obs x m var (depending on which axis to extend)
        :param axis:
        """
        adata = ad.concat([self, additional_adata], axis=axis)
        return adata

    def convert_to_df(self) -> pd.DataFrame:
        """
        convert anndata object into a long-form pandas dataframe. primary purpose is to allow access to pandas and seaborn functionality more directly.

        - overall seems to be working well. just need to test with a dataset with >1 obs and var keys(), and to test with the whole larger dataset.
        :return: long-form pandas dataframe

        """

        print(f"\n\- converting anndata data matrix to long-form pandas dataframe ... [in progress]")

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

        print(f"\n|- converting anndata data matrix to long-form pandas dataframe ... [finished]")

        return df














































































