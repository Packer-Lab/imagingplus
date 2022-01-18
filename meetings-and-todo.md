# Record for meeting agendas, notes and todo's 

## Jan 21 2022
- tested suite2p integration and results loading
- added integration of anndata data org package
- anndata is a nice, useful data structure for handling the multi-dimensional nature of our experiments

TODO:
- jupyter notebook tutorial for the experiment class

## Jan 14 2022
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