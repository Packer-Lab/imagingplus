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

### TwoPhotonImagingTrial class


### AllOpticalTrial class


## Data Structure Inside Trial objects

Numpy arrays are heavily used to collect/store raw and processed data inside each Trial object. 
There are multiple Trial object attributes which store data, with each attribute dedicated to a specific method of data processing.

***Annotated data object for storing multi-functional data***

The AnnData library is used to store data in an efficient, multi-functional format. This is stored under: `[trialobject.data](http://trialobject.data)` . The AnnData object is built around the raw Flu matrix of each `trialobject` . In keeping with AnnData conventions, the data structure is organized in *n* observations (obs) x *m* variables (var), where observations are suite2p ROIs and variables are imaging frame timepoints. The rest of the AnnData data object is built according to these dimensions. For instance, the metadata for each suite2p ROI stored in Suite2p’s stat.npy output is added to `trialobject.data` under `obs` and `obsm` (1D and >1-D observations annotations, respectively). And the temporal synchronization data of the experiment collected in .paq output is added to the variables annotations under `var`.

![packerlabimaging-anndata-integration.png](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/65f203b3-61f7-49bc-9955-9f44224d238a/packerlabimaging-anndata-integration.png)

Further processing on the raw data is added to the AnnData object as `layers`. For instance, dFF normalization of the raw data is added as the `dFF` layer to the existing AnnData object. 

The primary benefit of anndata is that it enforces an intuitive data structure and allows a workflow that maintains relationships between different data sources and types, including processed data. Lastly, AnnData objects are highly scalable, allowing individual users to further modify and add observations and variables of their own as their experiment dictates. Note: since AnnData is an independent library, they have an implementation for unstructured data annotations, however usage of this feature is not recommended within the packerlabimaging ecosystem as unstructured data elements should be built as attributes of the `trialobject`.
***Flow of the TwoPhotonImaging Experiment class***

1. paqProcessing
    - retrieves frame clock times, including adjustments for multiple starts/stops of imaging acquisition (chooses the longest imaging sequence in the paq file as the actual frame clocks for the present trial's acquisition, will not work if this assumption fails)

***Flow of the AllOptical Experiment class***

- `AllOpticalTrial.__init__()`: setting up the alloptical trial object
    1. 2p imaging `__init__()`
    2. stimProcessing
        1. unpacks Naparm outputs xml and gpl
            1. determines stim duration and stim properties
        2. unpacks paq File
            1. locates stim frames based on the specified `stim_channel`
            
    3. _findTargetsAreas
        - Loads target coordinates, and organizes target coordinates per SLM groups
        - creates target_areas - circle of pixels outward from the center of the target
        - creates target_areas_exclude - expanded circle of pixels outward from the center of the target that includes the region for excluding nontarget cells
        - creates images of targets // scaled somehow?? not totally sure - dont remember, has been commented out for a while now for myself - maybe check what Rob suggests adding here?
    4. _find_photostim_add_bad_framesnpy
        - finds all imaging frame that are overlapping with photostimulation trials, and creates bad_frames.npy file using these photostim frames for suite2p to exclude

- SLM targets processing+analysis: running processing and analysis specific to data collected from SLM targets
  1. making trace snippets from all SLM targets areas

- All cells processing+analysis: running processing and analysis specific to data collected from SLM targets: `photostimProcessingAllCells()`
    1. making trace snippets from all cells from suite2p:  `_makePhotostimTrialFluSnippets()`
    2. measure photostim dFF responses: create dataframe of cells x stims containing responses: `_collectPhotostimResponses()`
    3. statistical analysis of responses
    

## **Plotting of data/analysis**