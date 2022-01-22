import packerlabimaging as pkg

expobj = pkg.import_obj('/home/pshah/Documents/code/packerlabimaging/tests/RL109_analysis.pkl')
self = expobj.load_trial('t-013')

self.data

# self.photostimFluArray = self._makePhotostimTrialFluSnippets(plane_flu=normalize_dff(self.Suite2p.raw))

# self.photostimProcessingAllCells()
