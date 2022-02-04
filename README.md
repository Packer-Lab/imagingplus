# packerlabimaging package

packerlabimaging is a simple Python package for essential processing and analysis of 2photon imaging data collected in the Packer Lab. 
Especially, there is a fully complete code pipeline for analysis and plotting of standard Ca2+ imaging data, and of standard all optical experiments (both 2photon optogenetic stim and 1photon optogenetic stim, with combined 2photon Ca2+ imaging).

![Overall packerlabimaging package Flow Diagram](https://github.com/Packer-Lab/packerlabimaging/blob/7e16cf76588fa3fa34f634b9b455d9f386c54226/files/Overall%20Package%20Flow%20Diagram.drawio.png "Overall Flow Diagram")

This package is designed to work best with imaging experiments performed using PackIO, a Bruker 2pPlus microscope and using Suite2p for Ca2+ imaging 
data processing for ROI detection. Ultimately, the goal of this package is to jump-start your own analysis of your awesome experiment. 
It should be completely usable and understandable for anyone with the correct data in hand and basic knowledge of Python. There are tutorials to help along the way.  
We hope that it provides a useful structure to organize your experimental data, and some functionality to interact and process your data in an efficient manner. 

### Data Analysis Structure and Organization

The data organization of the package follows object-oriented programming in Python.
There are two general types of objects within packerlabimaging: the *Experiment* object and the *Trial* object.
In general, it follows that a typical 2photon imaging and/or all optical experiment contains multiple imaging trials of varying durations and/or combined 2photon imaging + optogenetic stim. 
The primary entry point of the package is the Experiment object. 
Multiple individual Trial objects are collected in a single Experiment object. 

The definition of an "experiment" is purposefully loose to allow for accommodation across different experimental designs.  
On the other hand, a "trial" is strictly a single time series of 2photon imaging data collection. 
Each trial may correspond to either a 2photon imaging only or an all optical trial. 

### Experiment class

Each imaging experiment is associated with a super-class Experiment object. Each Experiment object has optional access to Suite2p processing associated with that Experiment. 
Given this structure, each Experiment is meant to be built from all individual imaging trials that are closely related to each other based on their imaging characteristics. 
Most presciently, all imaging trials that are run together in Suite2p can be part of a single Experiment. 


***TwoPhotonImagingTrial class***

`trialobj: TwoPhotonImagingTrial` - packerlabimaging —> `TwoPhotonImagingTrial` trial class

Attributes:
* this class is the parent for all other trial classes

    `trialobj.imparams` - stores metadata retrieved from microscope regarding imaging collection 

    `trialobj.paq` - stores data from PackIO .paq files

    `trialobj.suite2p` - stores data and methods related to Suite2p (library for processing calcium imaging data)

    `trialobj.lfp` - stores local field potential data (read in from .paq files)

    `trialobj.data` - annotated data object (based on AnnData library) for centralized storage of raw and processed data related to the experiment. contains methods for further modifying the stored annotated data object.

***AllOpticalTrial class***

`trialobj: AllOpticalTrial` - packerlabimaging —> `TwoPhotonImagingTrial` —> `AllOpticalTrial` trial class

Attributes:

* this class is a child of the `TwoPhotonImaging` Trial (i.e. inherits all attributes and methods of `TwoPhotonImagingTrial`)

    `trialobj._2pstim` - stores data and methods related to processing of 2p stim protocols

    `trialobj.data` - additional 2p stim related variables (e.g. `photostim_frame` and `stim_start_frame` ) to `.data`


## Data Structure Inside Trial objects

Numpy arrays are heavily used to collect/store raw and processed data inside each Trial object. 
There are multiple Trial object attributes which store data, with each attribute dedicated to a specific method of data processing.

***Annotated data object for storing multi-functional data***

The AnnData library is used to store data in an efficient, multi-functional format. This is stored under: `trialobject.data` . The AnnData object is built around the raw Flu matrix of each `trialobject` . In keeping with AnnData conventions, the data structure is organized in *n* observations (obs) x *m* variables (var), where observations are suite2p ROIs and variables are imaging frame timepoints. The rest of the AnnData data object is built according to these dimensions. For instance, the metadata for each suite2p ROI stored in Suite2p’s stat.npy output is added to `trialobject.data` under `obs` and `obsm` (1D and >1-D observations annotations, respectively). And the temporal synchronization data of the experiment collected in .paq output is added to the variables annotations under `var`.

![AnnData integration for packerlabimaging package](https://github.com/Packer-Lab/packerlabimaging/blob/0967f1fb5f04a8407f5e574c2ebc8afbdf3f14de/files/packerlabimaging-anndata-integration-01.jpg "anndata for data storage")

Further processing on the raw data is added to the AnnData object as `layers`. For instance, dFF normalization of the raw data is added as the `dFF` layer to the existing AnnData object. 

The primary benefit of anndata is that it enforces an intuitive data structure and allows a workflow that maintains relationships between different data sources and types, including processed data. Lastly, AnnData objects are highly scalable, allowing individual users to further modify and add observations and variables of their own as their experiment dictates. Note: since AnnData is an independent library, they have an implementation for unstructured data annotations, however usage of this feature is not recommended within the packerlabimaging ecosystem as unstructured data elements should be built as attributes of the `trialobject`.


***Flow of the TwoPhotonImaging Experiment class***

- TwoPhotonImagingTrial
    1. paqProcessing
        - retrieves frame clock times, including adjustments for multiple starts/stops of imaging acquisition (chooses the longest imaging sequence in the paq file as the actual frame clocks for the present trial's acquisition, will not work if this assumption fails)
    2. adding Suite2p results
    3. dFF normalization
    4. creating  annotated data object using AnnData

![TwoPhoton Imaging Workflow #2.jpg](https://github.com/Packer-Lab/packerlabimaging/blob/4dd9ee035df2fd2e9ac7b1f3b82a7e7606d38492/files/TwoPhoton%20Imaging%20Workflow%20%232.jpg)


***Flow of the AllOptical Experiment class***

- `AllOpticalTrial.__init__()`: setting up the alloptical trial object
    1. 2p imaging `__init__()`
    2. `stimProcessing()`
        1. unpacks Naparm outputs xml and gpl
            1. determines stim duration and stim properties
        2. unpacks paq File
            1. locates stim frames based on the specified `stim_channel`
            
    3. `_findTargetsAreas()`
        - Loads target coordinates, and organizes target coordinates per SLM groups
        - creates target_areas - circle of pixels outward from the center of the target
        - creates target_areas_exclude - expanded circle of pixels outward from the center of the target that includes the region for excluding nontarget cells
        - creates images of targets // scaled somehow?? not totally sure - dont remember, has been commented out for a while now for myself - maybe check what Rob suggests adding here?
    4. `_find_photostim_add_bad_framesnpy()`
        - finds all imaging frame that are overlapping with photostimulation trials, and creates bad_frames.npy file using these photostim frames for suite2p to exclude

![alloptical-workflow-1.drawio.png](https://github.com/Packer-Lab/packerlabimaging/blob/4dd9ee035df2fd2e9ac7b1f3b82a7e7606d38492/files/alloptical-workflow-1.drawio.png)

- SLM targets processing+analysis: running processing and analysis specific to data collected from coordinates from registerred movie
    - [ ]  get input coordinates - alloptical workflow: SLM targets areas
    - [ ]  making trace snippets from all coords targets areas
    - [ ]  measure dFF responses across stims
    - [ ]  create new anndata object for storing measured photostim responses from data, with other relevant data for SLM targets
        
        anndata metadata:
        
        - [ ]  SLM groups
    
- All cells processing+analysis: running processing and analysis specific to data collected from all Suite2p ROIs: `photostimProcessingAllCells()`
    1. finding Suite2p ROIs that are also targets: `_findTargetedS2pROIs`
        - [ ]  test out finding suite2p ROI targets - save in `.targeted_cells`
        - [ ]  adding obs annotations of SLM group IDs to Suite2p ROIs targets
    
    `photostimProcessingAllCells()`:
    
    1. making trace snippets from all cells from suite2p:  `_makePhotostimTrialFluSnippets()`
    - [x]  measure photostim dFF responses: create dataframe of suite2p cells x stims containing responses: no independent method
    - [x]  statistical analysis of responses: `_runWilcoxonsTest()`
        - [ ]  singleTrialSignificance stats measurements?
    - [ ]  create new anndata object for storing measured photostim responses from data, with other relevant data
        
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