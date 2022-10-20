import os
import pickle

import pandas as pd

pd.options.display.max_rows = 100
pd.options.display.max_columns = 10


# import .pkl'd objects
def import_obj(pkl_path):
    """
    Imports objects that have been saved as .pkl files.

    :param pkl_path: path to .pkl file to load object from.
    :return: object loaded from .pkl file

    Example:
    >>> import imagingplus as pkg
    >>> trialobj = pkg.import_obj(pkl_path='/home/pshah/Documents/code/imagingplus/tests/2020-12-19_t-013.pkl')
    """

    if not os.path.exists(pkl_path):
        raise FileNotFoundError(f'pkl path NOT found: {pkl_path}')
    with open(pkl_path, 'rb') as f:
        print(f'\n\- loading {pkl_path} ... ', end='\r')
        try:
            obj = pickle.load(f)
        except pickle.UnpicklingError:
            raise pickle.UnpicklingError(f"\n** FAILED IMPORT from {pkl_path}\n")
        print(f'|- Loaded {obj.__repr__()}\n')

    return obj
