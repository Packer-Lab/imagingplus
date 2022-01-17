import os
import pickle

# import .pkl'd objects
def import_obj(pkl_path):
    if not os.path.exists(pkl_path):
        raise FileNotFoundError(f'pkl path NOT found: {pkl_path}' )
    with open(pkl_path, 'rb') as f:
        print(f'\- loading {pkl_path} ... ', end='\r')
        try:
            obj = pickle.load(f)
        except pickle.UnpicklingError:
            raise pickle.UnpicklingError(f"\n** FAILED IMPORT from {pkl_path}\n")
        print(f'|- Loaded {obj.__repr__()}')

    return obj
