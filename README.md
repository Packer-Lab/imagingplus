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
2. Install the conda environment provided in this repository (`plitest.yml`) using `conda env create -f myenv.yml` from the terminal. Note: The package is installable as a stand-alone python package. You can install the package into an existing conda environment, or you may choose to skip using conda environment all together. 
3. Activate the conda environment `conda activate plitest`.
4. `cd` to the parent directory of where this repo was downloaded to.
5. From this parent directly, run `pip install -e packerlabimaging` from terminal to install this package `packerlabimaging`.

#### Test `packerlabimaging` is successfully installed (import package in python):
1. Ensure that the conda environment from which `packerlabimaging` was installed is activated.
2. start python from command line using: `python`, or start python from the same conda environment in your preferred method (e.g. jupyter notebook or IDE).
3. Import the package: `import packerlabimaging` or `import packerlabimaging as pli`.

## Documentation

The documentation for `packerlabimaging` is not currently hosted online. 
Instead, you can access the documentation by opening the `index.html` file found from the following file path of the downloaded git repo on your local computer:

{packerlabimaging-git-repo local copy} > docs > build > index.html

Opening the `index.html` will open the documentation as an HTML file in your web-browser.


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
