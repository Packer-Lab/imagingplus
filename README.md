# Imaging+ <img src="docs/source/files/imagingplus-logo.png" width="200" title="imagingplus logo" alt="a python -x- calcium spike" align="right" vspace = "50">

imaging+ is a Python package and tool-suite for essential processing and analysis of 2photon imaging data. 
Especially, there is a fully implemented pipeline for data structuring, processing, analysis and plotting of 2photon Ca2+ imaging experiments, and experiments based around 2photon imaging such as all optical experiments (i.e. 2photon optogenetic stim or 1photon optogenetic stim, with combined 2photon Ca2+ imaging).

## Overview

This package is designed for experiments that follow the general structure of an "imaging+" experiment shown below. We have additionally provided specific sub-modules to suit imaging experiments performed using PackIO, a Bruker 2pPlus microscope and using [Suite2p](https://suite2p.readthedocs.io/en/latest/) for Ca2+ imaging 
data processing for ROI segmentation.

This package integrates the [AnnData](https://anndata.readthedocs.io) framework for native and highly efficient, HDF-5 compatible storage and organization of all levels of the data analysis workflow from preprocessing to highly analyzed downstream forms as well. 

It should be completely usable and understandable for anyone with the correct data in hand and basic knowledge of Python. 
There are tutorials to help along the way. Ultimately, the goal of this package is to jump-start your own analysis of your awesome experiment. 

<img src="docs/source/files/Typical-experiment-apr-22-2022.png" title="Imaging+ experiment diagram" alt="a typical imaging+ neuroscience experiment">


We hope that it provides a useful structure to organize your experimental data, and some functionality to interact and process your data in an efficient manner. 

## Getting started

Recommended tools:
- python >=3.9
- conda - environments for python package management
  - Run `conda create -n <insert-name> python=3.9` to create a new conda environment with python 3.9
  - Note: `imaging+` should install and function normally within existing conda environments as well. 
- Jupyter - for creating python based notebooks
- VS Code or PyCharm - IDE for python code development

## Installation instructions


#### Install pypi release:
`pip install imagingplus` -- not stable yet -- do not use while this message is here. 

#### Install Git version:
Note: The package is installable as a stand-alone python package. You can install the package into an existing conda environment, or you may choose to skip using a conda environment all together (steps 2 and 3).

In your terminal, run the following commands:
1. Clone this Github repository using `git clone https://github.com/Packer-Lab/imagingplus.git`. 
2. `cd` to the parent directory where this repo was downloaded to.
3. From this parent directly, run `pip install -e imagingplus` to install this package `imagingplus`.

#### Test `imaging+` is successfully installed (import package in python):
1. Ensure that the conda environment from which `imaging+` was installed is activated.
2. start python from command line using: `python`, or start python from the same conda environment in your preferred method (e.g. jupyter notebook or IDE).
3. Import the package: `import imagingplus as ip`.

## Documentation

The documentation for `imaging+` is available at https://imagingplus.readthedocs.io/en/latest/ 


