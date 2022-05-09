import os
import pickle
import re
import time


class WideFieldImaging:
    """
    WideField imaging cellsdata object.

    """

    def __init__(self, tiff_path, paq_path, exp_metainfo, pkl_path):
        self.tiff_path = tiff_path
        self.paq_path = paq_path
        self.metainfo = exp_metainfo
        # set and create analysis save s2pResultsPath directory
        self.analysis_saveDir = self.pkl_path[:[(s.start(), s.end()) for s in re.finditer('/', self.pkl_path)][-1][0]]
        if not os.path.exists(self.analysis_saveDir):
            print('making analysis save folder at: \n  %s' % self.analysis_saveDir)
            os.makedirs(self.analysis_saveDir)

        self.save_pkl(pkl_path=pkl_path)  # save experiment object to pkl_path

    def __repr__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
        else:
            lastmod = "(unsaved pkl object)"
        if not hasattr(self, 'metainfo'):
            information = f"uninitialized"
        else:
            information = self.t_series_name

        return repr(f"({information}) WidefieldImaging experimental cellsdata object, last saved: {lastmod}")

    @property
    def t_series_name(self):
        if "exp_id" in [*self.metainfo] and "trial_id" in [*self.metainfo]:
            return f'{self.metainfo["exp_id"]} {self.metainfo["trial_id"]}'
        else:
            raise ValueError('no information found to retrieve t series id')

    def save_pkl(self, pkl_path: str = None):
        if pkl_path is None:
            if hasattr(self, 'pkl_path'):
                pkl_path = self.pkl_path
            else:
                raise ValueError(
                    'pkl s2pResultsPath for saving was not found in cellsdata object attributes, please provide pkl_path to save to')
        else:
            self.pkl_path = pkl_path

        with open(self.pkl_path, 'wb') as f:
            pickle.dump(self, f)
        print("\n\t -- cellsdata object saved to %s -- " % pkl_path)

    def save(self):
        self.save_pkl()