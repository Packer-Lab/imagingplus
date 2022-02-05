"""
This is the main package for getting started with processing and analysis of calcium imaging data in the Packer lab.

The fundamental data object is one t-series of imaging data, along with its associated .paq file, accessory files generated
from the microscope during data collection, and any user generated files associated with the experiment of this t-series.

"""

import warnings
warnings.filterwarnings("ignore")

from ._version import __version__
from packerlabimaging._io import import_obj
from .ExperimentMain import Experiment, define_term
from packerlabimaging.utils._utils import Utils

print(f"import packerlabimaging\n\tversion: {__version__}\n")