# add electrophysiology data to the trial

# retrieving and processing on LFP recordings from the .Paq file
import os.path

from packerlabimaging.main.classes import TemporalData
from packerlabimaging.processing.paq import paq2py


class LFP:
    def __init__(self, data, **kwargs):
        # self.lfp_from_tmdata(chan_name, kwargs['_paq_path']) if '_paq_path' in [*kwargs] else KeyError('no `paq_path` provided to load LFP from.')
        self.data = data  #: 1-D LFP data


    @classmethod
    def lfp_from_tmdata(cls, tmdata: TemporalData, channel: str = 'voltage'):
        """Alternative constructor for LFP. """
        print('\n----- Creating new LFP analysis submodule from TemporalData ...')

        if channel not in tmdata.channels:
            raise ValueError('channel not found in provided tmdata')




