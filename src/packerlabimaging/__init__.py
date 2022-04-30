"""

packerlabimaging: an essential package for getting started with processing and analysis of calcium imaging data in the Packer lab.

The fundamental data object is one t-series of imaging data, along with its associated .Paq file, accessory files generated
from the microscope during data collection, and any user generated files associated with the experiment of this t-series.

"""

import warnings

warnings.filterwarnings("ignore")

from ._version import __version__
from .utils.io import import_obj
from packerlabimaging._archive.ExperimentMain import Experiment
from packerlabimaging.workflows.AllOptical import AllOpticalTrial
from packerlabimaging.workflows.OnePhotonStimImaging import OnePhotonStim
from .plotting import plotting
from .utils.utils import define_term

print(f"\nimported packerlabimaging successfully\n\tversion: {__version__}\n")
