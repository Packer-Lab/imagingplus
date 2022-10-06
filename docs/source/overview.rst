.. _overview:

**Core structures**
===================

There are two core types of objects that are used to store experimental data: `Experiment` and `ImagingTrial` objects.

`ImagingTrial` objects represent a single, continuous imaging+ trial.


The `Experiment` acts as a container to collect any number/types of individual `ImagingTrial` objects.

Together, the `Experiment` and `ImagingTrial` are the primary entry points to analysis with the package.
In general, this follows a typical imaging experiment design which might contain an arbitrary number and mixture of imaging+ trials.
The definition of an "experiment" is loose but, abstractly, can be thought of as a collection of imaging+ trials across a common imaging field-of-view.
On the other hand, a "trial" is strictly a single time series of imaging+ data collection.


Each `ImagingTrial` is built out-of three sub-types - `TemporalData`, `ImagingData` and `CellAnnotations` - which represent the three general modes of data outputs of an imaging+ data collection trial (Figure 2). Ultimately, the goal of this data processing flow is to organise data of an imaging+ experiment trial into an anndata data structure table that closely and intuitively represents and integrates multi-modal experimental data in a single Pythonic object (Figure 3). The natively functionality of the anndata library also allows for efficient on-disk and on-memory data access solution for imaging+ data.


**Modular Data Processing, Analysis and Visualisation**
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

Packerlabimaging is designed to be highly modular. This allows an as-per-need data analysis workflow utilising a mixture of pre-built analysis submodules and custom user-built modules. We provide a sub-module to integrate the use of Suite2p, which is currently the most popular calcium imaging region-of-interest segmentation library to generate `ImagingData` and `CellsAnnotations` structures from microscope .tiff fluorescence calcium imaging data. We provide a sub-module for .paq data generated from PackIO which is used for temporal synchronisation of multiple hardware components during an experiment, including for instance electrophysiological data from an electrode collecting data from a sample. We also provide pre-built support for parsing microscope acquisition metadata from Bruker microscopes. Non-Bruker imaging acquisition can be readily handled using a general imaging acquisition metadata class or creation of a custom class for the desired microscope system. An additional sub-module provides support for parsing holographic 2-photon SLM stimulation protocols generated using NAPARM (ref) for all-optical imaging+ trials, and for running the STAMovieMaker algorithm for quick analysis of all-optical experiments.

Finally, we also provide a plotting module which contains a number of convenient methods for typical plotting of imaging+ data.

**Extensibility**
+++++++++++++++++

The design structure of packerlabimaging is purposefully meant to enable direct extensibility on top of the core built-in structures of the package that allow end-users to create custom data processing/analysis workflows while benefiting from pre-engineered structures and workflows. For instance, packerlabimaging comes pre-built with a fully built-out workflow for AllOptical imaging+ trials which extends the core `ImagingTrial`. This workflow integrates individual data processing sub-modules to create a custom type of imaging+ trial designed for an all-optical experiment. Additionally methods for typical analysis of AllOptical data are provided herein.





Next
----

:ref:`Data structure details`

:ref:`tutorials`