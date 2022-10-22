"""

imagingplus: integrated tool-suite of essential processing and analysis of imaging+ experiments

"""

import warnings

warnings.filterwarnings("ignore")

from ._version import __version__
from .utils.io import import_obj
from imagingplus.main.core import SingleImage
from imagingplus.utils import utils
from imagingplus.utils import PrairieLink
from imagingplus.utils import images
from imagingplus.main.core import Experiment
from imagingplus.processing import imaging
from imagingplus.workflows.TwoPhotonImaging import TwoPhotonImaging
from imagingplus.workflows.AllOptical import AllOpticalTrial
from .plotting import plotting
from .utils.utils import define_term

print(f"\nimported imagingplus successfully\n\tversion: {__version__}\n")

LOCAL_DATA_PATH = '/Users/prajayshah/data/oxford-data-to-process/'
REMOTE_DATA_PATH = '/home/pshah/mnt/qnap/Data/'
BASE_PATH = LOCAL_DATA_PATH



