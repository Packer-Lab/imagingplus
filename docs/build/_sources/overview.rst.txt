.. _overview:

Overview
========

Data Analysis Structure and Organization
----------------------------------------

The data organization of the package follows object-oriented programming in Python.
There are two general types of objects within packerlabimaging: the *Experiment* object and the *Trial* object.
In general, it follows that a typical 2photon imaging and/or all optical experiment contains multiple imaging trials of varying durations and/or combined 2photon imaging + optogenetic stim.
The primary entry point of the package is the Experiment object.
Multiple individual Trial objects are collected in a single Experiment object.

The definition of an "experiment" is purposefully loose to allow for accommodation across different experimental designs.
On the other hand, a "trial" is strictly a single continuous time series of 2photon imaging data collection.
Each trial may correspond to any of the implemented trial types. Currently, the implemented trial types are: 2-photon imaging only, all optical, one photon stim + 2p imaging.

Experiment class
++++++++++++++++

Each imaging experiment is associated with a super-class Experiment object. Each Experiment object has optional access to Suite2p processing associated with that Experiment.
Given this structure, each Experiment is meant to be built from all individual imaging trials that are closely related to each other based on their imaging characteristics.
Most presciently, all imaging trials that are run together in Suite2p can be part of a single Experiment.


TwoPhotonImagingTrial class
+++++++++++++++++++++++++++

``trialobj: TwoPhotonImagingTrial`` - packerlabimaging —> ``TwoPhotonImagingTrial`` trial class

Attributes:
* this class is the parent for all other trial classes

    ``trialobj.imparams`` - stores metadata retrieved from microscope regarding imaging collection

    ``trialobj.paq`` - stores data from PackIO .paq files

    ``trialobj.suite2p`` - stores data and methods related to Suite2p (library for processing calcium imaging data)

    ``trialobj.lfp`` - stores local field potential data (read in from .paq files)

    ``trialobj.data`` - annotated data object (based on AnnData library) for centralized storage of raw and processed data related to the experiment. contains methods for further modifying the stored annotated data object.

AllOpticalTrial class
+++++++++++++++++++++

``trialobj: AllOpticalTrial`` - packerlabimaging —> ``TwoPhotonImagingTrial`` —> ``AllOpticalTrial`` trial class

Attributes:

* this class is a child of the ``TwoPhotonImaging`` Trial (i.e. inherits all attributes and methods of ``TwoPhotonImagingTrial``)

    ``trialobj._2pstim`` - stores data and methods related to processing of 2p stim protocols

    ``trialobj.data`` - additional 2p stim related variables (e.g. ``photostim_frame`` and ``stim_start_frame`` ) to ``.data``
