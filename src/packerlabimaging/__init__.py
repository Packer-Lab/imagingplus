"""

packerlabimaging: an essential data processing and analysis package for getting started with processing and analysis of calcium imaging data.

"""

import warnings

warnings.filterwarnings("ignore")

from ._version import __version__
from .utils.io import import_obj
from packerlabimaging.main.core import SingleImage
from packerlabimaging.utils import utils
from packerlabimaging.utils import images
from packerlabimaging.main.core import Experiment
from packerlabimaging.processing import imaging
from packerlabimaging.workflows.TwoPhotonImaging import TwoPhotonImaging
from packerlabimaging.workflows.AllOptical import AllOpticalTrial
from .plotting import plotting
from .utils.utils import define_term

print(f"\nimported packerlabimaging successfully\n\tversion: {__version__}\n")

LOCAL_DATA_PATH = '/Users/prajayshah/data/oxford-data-to-process/'
REMOTE_DATA_PATH = '/home/pshah/mnt/qnap/Data/'
BASE_PATH = LOCAL_DATA_PATH



