# Record for meeting agendas, notes and todo's 

## Apr 19 2022 - TODO list:

**Major code base related tasks:**

- [ ]  create the three fundamental data analysis *axes as fundamental classes . then plug other stuff into those axes.*
    - like Temporal axis
        - plug paq into temporal axes
        - allows for other temporal signals from outside paq to plug in, but Paq becomes a child of this class
    - data axis
        - suite2p analysis
        - this might be fundamentally just the parent trial class
    - cell annotations axis
        - suite2p analysis
        - can also include manual labelling of cells as well (e.g. all optical SLM targets or
        - allows for ROIs from other algorithms
- [ ]  testing and fleshing out of AllOpticalTrial code
- [ ]  making trial subtypes child objects of general parent `Trial` type

**Major packaging related tasks:**

- [ ]  Tutorial addition: Interactive plotting and data exporting using mpl_point_clicker
- [x]  Add instructions to the README.md for installation of the package
- [x]  implementing tests

**Less important tasks/considerations:**


## Feb 28 2022 - TODO list transferred over:
**Major code base related tasks:**

- [x]  need to debug `_findTargetsAreas` and `targetCoordinates`(already ipr)
- [x]  add 2pstim module - for NAPARM related funcs and methods
    - [x]  update 2p stim related attr’s to use the naparm submodule
- [x]  refactor out stats testing
- [x]  add new class for `STAMovieMaker_nogui` - or add as method for `AllOpticalTrial`?
- [x]  adding plotting module
- [x]  add onePstim module
    - [ ]  need to edit `__init__` to fit into the package pipeline
- [x]  consider making anndata extension funcs as staticmethods, so they are more easily (?potentially) accessible widely - no need...
- [x]  finalize plotting sub-module
    - [x]  plotROIlocations test/debug
- [x]  update Two-photon processing and Alloptical processing tutorials with processing steps (e.g. dFF normalization, calculating photostim responses, sig diff testing), and also plots
    - all optical first pass through should be done to review with Adam?
    - twophoton first pass through also pretty solid

**Major packaging related tasks:**

- [ ]  Tutorial addition: Interactive plotting and data exporting using mpl_point_clicker
- [x]  Add instructions to the README.md for installation of the package
- [ ]  figure out how to release an alpha version of the package
    - [x]  add documentation for installation of package
- [ ]  writing out documentation for user-facing functions/methods
- [x]  implementing tests

**Less important tasks/considerations:**

- ~~[ ]  try out making trial objs children of Experiments~~
    ~~- allows loading of trial objects from: `expobj.t-001`~~
    ~~- [ ]  need to remove requirement for providing `suite2p_experiment_obj`, and `total_frames_stitched`  args for `TwoPhotonImagingTrial` (probably part of making trial objs children of Experiment)~~




## Jan 21 2022 - Meeting with Adam
- tested suite2p integration and results loading
- added integration of anndata data org package
- anndata is a nice, useful data structure for handling the multi-dimensional nature of our experiments
- very basic jupyter notebook tutorial for the experiment class

TODO:

Meeting with Adam:
- show Adam the basic workflow
- talk about the AnnData framework
- concerns over code being too ecosystem specific - i.e Bruker imaging, Bruker all optical, NAPARM for photostim protocols, PackIO for synchronization?
  - I think it's quite difficult to avoid the code becoming tied to the ecosystem...
  - but how much should we try to stay open to other ecosystems?
- which plots do we want to include convenience functions for?


## Jan 14 2022 - Summary for Rob
- lots of work on smoothing out and cleaning up workflow of Experiment class and trial objects, and their Suite2p classes and individual Twophoton imaging and alloptical trial classes
- working on getting attr's in place in the __init__ for the TwoPhotonImaging and AllOptical Trial objects (and also implementing code for some of them)


Items done of note:
- copy pasted + integrated Rob's code for:
  - _parsePVMetadata and _getPVStateShard
  - _retrieveSuite2pData, _makeFluTrials, (_detrendFluTrial, baselined_flu_trial) - for last two see Question below
  - added the .time attr as a property under Alloptical

Questions for Rob:
- (_detrendFluTrial, baselined_flu_trial) - ASK: how important are these functions? also see that detrending is commented out in makeFluTrials - should we include or not?

- Todos:
- update comments descriptions for statistical analysis of photostim attr's (under AllOpticalTrial class) - tagged Rob


## Jan 7 2022
- discuss implementation of parent experiment class for handling multiple sub-experiment objects, and used to run suite2p from
- main todo with Rob: go through and ensure that all the alloptical functions are in place, correct and uptodate with his versions
  (can even create a pipeline to quickly jot down the flow of functions that will be used for an alloptical experiment)
-- then I can handle putting the structure of everything together.

Next big tasks:
[] test linking of Suite2p class to experiment object 
[] link Suite2p class to trial objects
[] review Rob's notes on attr's and functions for alloptical and 2p imaging classes


## Dec 22 2021
- installation works well for Rob - happy with the 
- talking about experiment metadata to extract and include in the expobj
- csv export for experimet metadata - idea 

TODO:
- bad_frames.npy processing - how to best setup bad_frames.npy when you have multiple tiffs
- setting up parent object to run suite2p for multiple tiffs
- how to best setup code to run suite2p for multiple tiffs 
  - two approaches ideas:
    1) independent class for running suite2p and setting up experiment association - Rob likes the sound of this
    2) associating all related expobj inside each expobj - and then using these associations to run suite2p from the expobj
- continue working on correcting attr's for the alloptical object

Rob:
- look through attr's for alloptical object and add to __init__ to initialize the attr's



## Dec 6 2021
- get Rob setup with development installation of the package
  - installation in new environment was successful
  - had problems with dependencies (need to figure out how to set dependencies)


- discuss first example Jupyter notebook of simple standard TwoPhoton imaging only workflow - thoughts, suggestions, ideas for improvements?
- discuss next workflows to work on.


Next TODOs:
  - add necessary dependencies to the setup.py 
    - also CLEAN UP IMPORT STATEMENTS and do testing with new dependencies
  - Prajay will work on creating a barebones alloptical initial object
    - next time: work on going through the code together and see alloptical  