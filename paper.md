---
title: 'An integrated analysis tool-kit in Python for multi-modal, imaging+ neuroscience data'
tags:
  - Python
  - neuroscience
  - calcium imaging
  - all optical
  - optogenetics
authors:
  - name: Prajay Shah
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Rob M. Lees
    affiliation: 2
  - name: James Rowland
    affiliation: 2
  - name: Adam M. Packer
    affiliation: 2
affiliations:
  - name: Institute of Biomedical Engineering, University of Toronto, ON, Canada
    index: 1
  - name: Department of Physiology, Anatomy and Genetics, University of Oxford, UK
    index: 2
date: xx XXX 2022
bibliography: paper.bib

# Summary
Calcium imaging is a mainstay in modern neuroscience. Additionally, optical techniques for stimulation of neurons have also become standard. Packerlabimaging is a Python package which provides a structured analytic framework and an integrated data-to-results pipeline for complex, multimodal neuroscience experiments which are based primarily around calcium imaging. It is designed to be scalable, extendable and modular to dovetail users’ personal data analysis needs with other existing Python libraries and custom home-built libraries. 

# A Statement of Need
Microscope-based neuroscience experiments are often multi-modal, combining electrophysiology, real-time behavioural readouts (e.g. locomotion, eye-tracking for brain state classification) and post-hoc molecular biology analysis (e.g. immunohistological classification of imaged neurons) to collect data from an individual experiment or sample. We refer to the conglomeration of data outputs around a primary source of live-imaging data as “imaging+” data. However, such experiments call for modular software, and there is a need for analysis workflows that match the growing complexity of experimental data. Traditionally, experimenters develop analysis workflows in an ad hoc manner to match the data analysis needs of their experiments. Reusable code libraries (i.e. packages) are generally only designed to implement a specific step along the data processing, analysis and/or results generation track. Ultimately, this leads to a lack of standardisation across the community, contributes to potential difficulties in reproducibility of analyses/results, and entails a long road of code development for neuroscientists as they begin new experiments. This can be a significant challenge for the growing number of beginner/novice programmers who engage in experiments that might ultimately require extensive and bespoke software development process. To facilitate the rapidly growing complexity of imaging+ experiments in neuroscience, we developed a full-service analysis toolkit in Python for analysing data from multi-modal neuroscience experiments that are based centrally around 2-photon Ca2+ imaging data and analysis. 

# Overview
Packerlabimaging was first designed to create a pipeline for the data analytics needs of a complex all-optical experiment that combines 2-photon imaging, 2-photon stimulation and electrophysiology data collection. Packerlabimaging creates the necessary structured framework to enable efficient analysis workflow and results generation from this type of imaging+ data.
Especially, there is a fully complete code pipeline for analysis and plotting of standard Ca2+ imaging data, and of standard all optical experiments (both 2photon optogenetic stim or 1photon optogenetic stim, with simultaneous 2photon Ca2+ imaging). We provide a structure and organisation of data that is designed to be intuitive and easily extended to serve the end-user’s unique needs. We utilise excellent pre-existing Python libraries (e.g. numpy, matplotlib, suite2p, anndata) along with custom code developed over years for our systems to deliver an analysis framework that aims to be widely replicable, modular, and easy to use out of the box. Our framework gives the user a single toolkit that integrates multiple modes of data processing and analyses required for multi-modal experiments. Extensive utilisation of object-oriented programming principles allow for modularity and extensibility with bespoke end-user code.
Out-of-the-box, the package is designed to work completely with imaging experiments performed using a Bruker 2pPlus microscope system, PACKIO for experimental hardware synchronisation and integrates Suite2p for processing Ca2+ imaging data. Ultimately, the goal of this package is to provide a useful structure to organise the user’s experimental data in an intuitive manner, and functionality for data processing, data exploration and data visualisation. With its end-to-end functionality, this package enables a coding experience that can move the user from imaging+ data collection to results with high efficiency.
Organisation and Workflow
By strictly following object-oriented programming principles, packerlabimaging enables modularity of data analysis and, most importantly, easy extensibility by the end-user.
There are two overall structural concepts implemented for the analysis flow of the package. The first relates to the representation of imaging+ data within Pythonic objects, and the second relates to the diverse set of sub-modules which can act upon these objects to run specific data processing and/or analysis procedures.
There are two primary types of objects that are used to represent the totality of experimental data: `Experiment` and `-Trial` objects. `-Trial` objects represent a single series of imaging+ data corresponding to the type of experimental trial performed (e.g. 2-photon imaging only vs. all-optical). The `Experiment` collects any number/types of individual -Trial objects. Together, the `Experiment` and `-Trial` are the primary entry points to analysis with the package. In general, this follows a typical imaging experiment design which might contain an arbitrary number of imaging trials and an arbitrary mixture of trial types.
The definition of an "experiment" is purposefully loose to allow accommodation across different experimental designs. However, abstractly, one Experiment object can be thought of as a collection of imaging trials that are processed together. On the other hand, a "trial" is strictly a single time series of 2photon imaging data collection. packerlabimaging comes with two built-in trial types: TwoPhotonImagingTrial (which represents 2photon imaging only trials) and AllOpticalTrial (which represents 2photon imaging + 2photon stimulation trials).

# Modular Data Analysis

packerlabimaging is intentionally designed to be modular. This means a light-weight data processing workflow when users first use the package, and a subsequent ‘as-per-need’ choice of data analysis workflows from either pre-built analysis submodules or user-built modules. These submodules include excellent pre-existing Python libraries. Two examples are: 1) the suite2p pipeline for calcium imaging movie registration, cell segmentation and spike deconvolution; and 2) anndata, which is a package for handling annotated, multi-dimensional data matrices that was originally developed to handle complex big datasets in single cell genomics. Anndata is used as the primary method of data organisation in packerlabimaging and provides an intuitive and efficient on-disk and on-memory access solution for multi-modal data storage of imaging+ data.

We provide support for .paq data generated from PackIO which is used for temporal synchronization of multiple hardware components during an experiment, including for instance electrophysiological data from an electrode collecting data from a sample. We also provide pre-built support for parsing imaging acquisition metadata from Bruker microscopes, one of the popular microscope systems in neuroscience. Non-Bruker imaging acquisition can be readily handled using a general imaging acquisition metadata class or creation of a custom class for the desired microscope system.

Additionally, there is support for parsing 2-photon stimulation protocols from all-optical experiments generated using NAPARM, and for running the STAMovieMaker algorithm for visual analysis all-optical experiments.

# Results and Plotting



# Citations

# Figures

# Acknowledgements
