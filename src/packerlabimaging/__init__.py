"""
This is the main package for getting started with processing and analysis of calcium imaging data in the Packer lab.

The fundamental data object is one t-series of imaging data, along with its associated .paq file, accessory files generated
from the microscope during data collection, and any user generated files associated with the experiment of this t-series.

"""

from .version import __version__
from packerlabimaging.io import import_obj
from .packerlabimaging import Experiment, define_term
# from .packerlabimaging import Experiment, import_obj, define_term