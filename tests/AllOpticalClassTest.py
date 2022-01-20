import packerlabimaging as pkg
from packerlabimaging._utils import normalize_dff

self = pkg.import_obj('/home/pshah/Documents/code/packerlabimaging/tests/2020-12-19_t-013.pkl')

# self.photostimFluArray = self._makePhotostimTrialFluSnippets(plane_flu=normalize_dff(self.Suite2p.raw))

self.photostimProcessingAllCells()
