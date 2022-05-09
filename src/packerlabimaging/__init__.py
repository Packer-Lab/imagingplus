"""

packerlabimaging: an essential package for getting started with processing and analysis of calcium imaging cellsdata in the Packer lab.

The fundamental cellsdata object is one t-series of imaging cellsdata, along with its associated .Paq file, accessory files generated
from the microscope during cellsdata collection, and any user generated files associated with the experiment of this t-series.

"""

import warnings

warnings.filterwarnings("ignore")

from ._version import __version__
from .utils.io import import_obj
from packerlabimaging.main.classes import Experiment
from packerlabimaging.workflows.TwoPhotonImaging import TwoPhotonImaging
from packerlabimaging.workflows.AllOptical import AllOpticalTrial
# from packerlabimaging.workflows.OnePhotonStimImaging import OnePhotonStim
from .plotting import plotting
from .utils.utils import define_term

print(f"\nimported packerlabimaging successfully\n\tversion: {__version__}\n")

LOCAL_DATA_PATH = '/Users/prajayshah/data/oxford-data-to-process/'
REMOTE_DATA_PATH = '/home/pshah/mnt/qnap/Data/'
BASE_PATH = LOCAL_DATA_PATH



