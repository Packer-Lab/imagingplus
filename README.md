# packerlabimaging package

packerlabimaging is a simple Python package for essential processing and analysis of 2photon imaging data collected in the Packer Lab. 
Especially, there is a fully implemented pipeline for data structuring, processing, analysis and plotting of 2photon Ca2+ imaging experiments, and experiments based around 2photon imaging such as all optical experiments (i.e. 2photon optogenetic stim or 1photon optogenetic stim, with combined 2photon Ca2+ imaging).

![Overall packerlabimaging package Flow Diagram](https://github.com/Packer-Lab/packerlabimaging/blob/7e16cf76588fa3fa34f634b9b455d9f386c54226/files/Overall%20Package%20Flow%20Diagram.drawio.png "Overall Flow Diagram")

This package is designed to work best with imaging experiments performed using PackIO, a Bruker 2pPlus microscope and using Suite2p for Ca2+ imaging 
data processing for ROI detection. Ultimately, the goal of this package is to jump-start your own analysis of your awesome experiment. 
It should be completely usable and understandable for anyone with the correct data in hand and basic knowledge of Python. There are tutorials to help along the way.  
We hope that it provides a useful structure to organize your experimental data, and some functionality to interact and process your data in an efficient manner. 

## Getting started

Recommended tools:
- python >=3.9
- conda - environments for python package management
  - Run `conda create -n <insert-name> python=3.9` to create a new conda environment with python 3.9
  - Note: `packerlabimaging` should install and function normally within existing conda environments as well. 
- Jupyter - for creating python based notebooks
- VS Code or PyCharm - IDE for python code development

## Installation instructions

Note: The package is installable as a stand-alone python package. You can install the package into an existing conda environment, or you may choose to skip using conda environment all together.

1. Clone this github repository using `git clone https://github.com/Packer-Lab/packerlabimaging.git` in the terminal. 
2. Install the conda environment provided in this repository (`plitest.yml`) using `conda env create -f plitest.yml` from the terminal. 
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

Opening the `index.html` file will open the documentation as an HTML file in your web-browser.

