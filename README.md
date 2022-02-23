# packerlabimaging package

packerlabimaging is a simple Python package for essential processing and analysis of 2photon imaging data collected in the Packer Lab. 
Especially, there is a fully implemented pipeline for data structuring, processing, analysis and plotting of 2photon Ca2+ imaging experiments, and experiments based around 2photon imaging such as all optical experiments (i.e. 2photon optogenetic stim or 1photon optogenetic stim, with combined 2photon Ca2+ imaging).

![Overall packerlabimaging package Flow Diagram](https://github.com/Packer-Lab/packerlabimaging/blob/7e16cf76588fa3fa34f634b9b455d9f386c54226/files/Overall%20Package%20Flow%20Diagram.drawio.png "Overall Flow Diagram")

This package is designed to work best with imaging experiments performed using PackIO, a Bruker 2pPlus microscope and using Suite2p for Ca2+ imaging 
data processing for ROI detection. Ultimately, the goal of this package is to jump-start your own analysis of your awesome experiment. 
It should be completely usable and understandable for anyone with the correct data in hand and basic knowledge of Python. There are tutorials to help along the way.  
We hope that it provides a useful structure to organize your experimental data, and some functionality to interact and process your data in an efficient manner. 

## Installation instructions

1. Clone this github repository using `git clone https://github.com/Packer-Lab/packerlabimaging.git` in the terminal. 
2. Install the conda environment provided in this repository (`myenv.yml`) using `conda env create -f myenv.yml` from the terminal. Note: The package is installable as a stand-alone python package. You can install the package into an existing conda environment, or you may choose to skip using conda environment all together. 
3. Activate the conda environment `conda activate plitest`.
4. `cd` to the parent directory of where this repo was downloaded to.
5. From this parent directly, run `pip install -e packerlabimaging` from terminal to install this package `packerlabimaging`.

#### Test `packerlabimaging` is successfully installed (import package in python):
1. Ensure that the conda environment from which `packerlabimaging` was installed is activated.
2. start python from command line using: `python`, or start python from the same conda environment in your preferred method (e.g. jupyter notebook or IDE).
3. Import the package: `import packerlabimaging` or `import packerlabimaging as pli`.


## Overview
### Data Analysis Structure and Organization

The data organization of the package follows object-oriented programming in Python.
There are two general types of objects within packerlabimaging: the *Experiment* object and the *Trial* object.
In general, it follows that a typical 2photon imaging and/or all optical experiment contains multiple imaging trials of varying durations and/or combined 2photon imaging + optogenetic stim. 
The primary entry point of the package is the Experiment object. 
Multiple individual Trial objects are collected in a single Experiment object. 

The definition of an "experiment" is purposefully loose to allow for accommodation across different experimental designs.  
On the other hand, a "trial" is strictly a single continuous time series of 2photon imaging data collection. 
Each trial may correspond to any of the implemented trial types. Currently, the implemented trial types are: 2-photon imaging only, all optical, one photon stim + 2p imaging.  

### Experiment class

Each imaging experiment is associated with a super-class Experiment object. Each Experiment object has optional access to Suite2p processing associated with that Experiment. 
Given this structure, each Experiment is meant to be built from all individual imaging trials that are closely related to each other based on their imaging characteristics. 
Most presciently, all imaging trials that are run together in Suite2p can be part of a single Experiment. 


### TwoPhotonImagingTrial class

`trialobj: TwoPhotonImagingTrial` - packerlabimaging —> `TwoPhotonImagingTrial` trial class

Attributes:
* this class is the parent for all other trial classes

    `trialobj.imparams` - stores metadata retrieved from microscope regarding imaging collection 

    `trialobj.paq` - stores data from PackIO .paq files

    `trialobj.suite2p` - stores data and methods related to Suite2p (library for processing calcium imaging data)

    `trialobj.lfp` - stores local field potential data (read in from .paq files)

    `trialobj.data` - annotated data object (based on AnnData library) for centralized storage of raw and processed data related to the experiment. contains methods for further modifying the stored annotated data object.

### AllOpticalTrial class

`trialobj: AllOpticalTrial` - packerlabimaging —> `TwoPhotonImagingTrial` —> `AllOpticalTrial` trial class

Attributes:

* this class is a child of the `TwoPhotonImaging` Trial (i.e. inherits all attributes and methods of `TwoPhotonImagingTrial`)

    `trialobj._2pstim` - stores data and methods related to processing of 2p stim protocols

    `trialobj.data` - additional 2p stim related variables (e.g. `photostim_frame` and `stim_start_frame` ) to `.data`


## Data Structure Inside Trial objects

Numpy arrays are heavily used to collect/store raw and processed data inside each Trial object. 
There are multiple Trial object attributes which store data, with each attribute dedicated to a specific method of data processing.

### ***Annotated data object for storing multi-functional data***

The AnnData library is used to store data in an efficient, multi-functional format. This is stored under: `trialobject.data` . The AnnData object is built around the raw Flu matrix of each `trialobject` . In keeping with AnnData conventions, the data structure is organized in *n* observations (obs) x *m* variables (var), where observations are suite2p ROIs and variables are imaging frame timepoints. The rest of the AnnData data object is built according to these dimensions. For instance, the metadata for each suite2p ROI stored in Suite2p’s stat.npy output is added to `trialobject.data` under `obs` and `obsm` (1D and >1-D observations annotations, respectively). And the temporal synchronization data of the experiment collected in .paq output is added to the variables annotations under `var`.

![AnnData integration for packerlabimaging package](https://github.com/Packer-Lab/packerlabimaging/blob/0967f1fb5f04a8407f5e574c2ebc8afbdf3f14de/files/packerlabimaging-anndata-integration-01.jpg "anndata for data storage")

Further processing on the raw data is added to the AnnData object as `layers`. For instance, dFF normalization of the raw data is added as the `dFF` layer to the existing AnnData object. 

The primary benefit of anndata is that it enforces an intuitive data structure and allows a workflow that maintains relationships between different data sources and types, including processed data. Lastly, AnnData objects are highly scalable, allowing individual users to further modify and add observations and variables of their own as their experiment dictates. Note: since AnnData is an independent library, they have an implementation for unstructured data annotations, however usage of this feature is not recommended within the packerlabimaging ecosystem as unstructured data elements should be built as attributes of the `trialobject`.


### ***Flow of the TwoPhotonImaging Experiment class***

- TwoPhotonImagingTrial
    1. paqProcessing
        - retrieves frame clock times, including adjustments for multiple starts/stops of imaging acquisition (chooses the longest imaging sequence in the paq file as the actual frame clocks for the present trial's acquisition, will not work if this assumption fails)
    2. adding Suite2p results
    3. dFF normalization
    4. creating  annotated data object using AnnData

![TwoPhoton Imaging Workflow #2.jpg](https://github.com/Packer-Lab/packerlabimaging/blob/4dd9ee035df2fd2e9ac7b1f3b82a7e7606d38492/files/TwoPhoton%20Imaging%20Workflow%20%232.jpg)


### ***Flow of the AllOptical Experiment class***

![alloptical-workflow-1.drawio.png](https://github.com/Packer-Lab/packerlabimaging/blob/4dd9ee035df2fd2e9ac7b1f3b82a7e7606d38492/files/alloptical-workflow-1.drawio.png)

- `AllOpticalTrial.__init__()`: setting up the alloptical trial object
    1. 2p imaging `__init__()`
    2. `_paqProcessingAllOptical()`
        1. unpacks paq File
        2. locates stim frames based on the specified `stim_channel`
        3. returns stim timings as frame numbers
    3. `_stimProcessing()`
        1. processes information about the 2p stim protocol
            1. setup to use `naparm` and `Targets(naparm)` submodules to unpack naparm outputs
            2. 
        
    4. `_find_photostim_add_bad_framesnpy()`
        - finds all imaging frame that are overlapping with photostimulation trials, and creates bad_frames.npy file using these photostim frames for suite2p to exclude
    - 5. SLM targets processing+analysis: running processing and analysis specific to data collected from coordinates from registered movie
        1. `collect_traces_from_targets`
            1. *uses registered tiffs to collect raw traces from SLM target areas*
        2. `get_alltargets_stim_traces_norm`
            1. *primary function to measure the dFF and dF/setdF trace SNIPPETS for photostimulated targets*
        - [ ]  create new anndata object for storing measured photostim responses from data, with other relevant data
            
            `_allCellsPhotostimResponsesAnndata`
            
            - [ ]  create new anndata object for storing measured photostim responses from data, with other relevant data for SLM targets
                
                anndata module for alloptical trial object: `.photostimResponsesData`
                
                metadata to add to anndata object:
                
                - [ ]  SLM groups
            
            anndata metadata: - function in place, not tested really yet though. 
            
            - [ ]  SLM targets or not
            - [ ]  SLM photostim exclusion region
            - [ ]  `prob_response`
                - [ ]  - need to review code for calculating prob response
        
    - 6. All cells processing+analysis: running processing and analysis specific to data collected from all Suite2p ROIs: `photostimProcessingAllCells()`
        1. finding Suite2p ROIs that are also targets: `_findTargetedS2pROIs`
            - [ ]  test out finding suite2p ROI targets - save in `.targeted_cells`
            - [ ]  adding obs annotations of SLM group IDs to Suite2p ROIs targets
        
        `photostimProcessingAllCells()`:
        
        1. making trace snippets from all cells from suite2p:  `_makePhotostimTrialFluSnippets()`
        - [x]  measure photostim dFF responses: create dataframe of suite2p cells x stims containing responses: no independent method
        - [x]  statistical analysis of responses: `_runWilcoxonsTest()`
            - [ ]  singleTrialSignificance stats measurements?
            
        - [ ]  create new anndata object for storing measured photostim responses from data, with other relevant data
            
            `_allCellsPhotostimResponsesAnndata`
            
            - [ ]  create new anndata object for storing measured photostim responses from data, with other relevant data for SLM targets
                
                anndata module for alloptical trial object: `.photostimResponsesData`
                
                metadata to add to anndata object:
                
                - [ ]  SLM groups
            
            anndata metadata: - function in place, not tested really yet though. 
            
            - [ ]  SLM targets or not
            - [ ]  SLM photostim exclusion region
            - [ ]  `prob_response`
                - [ ]  - need to review code for calculating prob response




## **Plotting of data/analysis**

### functions

1) `plot_SLMtargets_Locs()`

2) `plot_cells_loc()`

3) `s2pRoiImage()`

4) `plot_flu_trace()`

5) `plotMeanRawFluTrace()`

6) `plot_s2p_raw()`

7) `plotLfpSignal()`

8) `plot_SLMtargets_Locs()`

9) `plot_photostim_traces()`

10) `plot_photostim_traces_overlap()`

11) `plot_periphotostim_avg2()`

12) `plot_periphotostim_avg()`



## TODO items:

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

- [ ]  Interactive plotting using mpl_point_clicker
- [ ]  Add instructions to the README.md for installation of the package
- [ ]  figure out how to release an alpha version of the package
    - [ ]  add documentation for installation of package
- [ ]  writing out documentation for user-facing functions/methods
- [ ]  implementing tests
- [ ]  figuring out how to cache during running tests??

**Less important tasks/considerations:**

- [ ]  try out making trial objs children of Experiments
    - allows loading of trial objects from: `expobj.t-001`
    - [ ]  need to remove requirement for providing `suite2p_experiment_obj`, and `total_frames_stitched`  args for `TwoPhotonImagingTrial` (probably part of making trial objs children of Experiment)
