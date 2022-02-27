.. _overview:

**Data Analysis Flow and Organization**
=======================================

The structure of the package follows object-oriented programming principles. This design philosophy is used to enable modularity of data analysis and most importantly easy extensibility by the end-user towards their own data analysis needs for custom experiments.

There are two primary types of objects within ``packerlabimaging``: the ``Experiment`` class and the ``Trial`` object. The ``Experiment`` class is the primary entry-point into the ``packerlabimaging`` pipeline and is the super-class under which any number/types of individual ``Trial`` objects can be collected. In general, this follows a typical imaging experiment design which might contain an arbitrary number of imaging trials and an arbitrary mixture of trial types.

The definition of an *"experiment"* is purposefully loose to allow accommodation across different experimental designs. However, abstractly, one ``Experiment`` object can be thought of as a collection of trials that are processed together in Suite2p. On the other hand, a *"trial"* is strictly a single time series of 2photon imaging data collection.

The trials types that ``packerlabimaging`` is currently built for are: ``TwoPhotonImagingTrial`` (which represents 2photon imaging only trials) and ``AllOpticalTrial`` (which represents 2photon imaging + 2photon stimulation trials).



Submodules: Experiment
----------------------

**Experiment class**

Each imaging experiment is associated with a super-class Experiment object. Each Experiment object has optional access to Suite2p processing associated with that Experiment. Given this structure, each Experiment is meant to be built from all individual imaging trials that are closely related to each other based on their imaging characteristics. Most presciently, all imaging trials that are run together in Suite2p can be part of a single Experiment.



Submodules: Trial types
-----------------------

**TwoPhotonImagingTrial class**

``trialobj: TwoPhotonImagingTrial`` - packerlabimaging —> ``TwoPhotonImagingTrial`` trial class


*Attributes:*

    ``trialobj.imparams`` - stores metadata retrieved from microscope regarding imaging collection

    ``trialobj.paq`` - stores data from PackIO .paq files

    ``trialobj.suite2p`` - stores data and methods related to Suite2p (library for processing calcium imaging data)

    ``trialobj.lfp`` - stores local field potential data (read in from .paq files)

    ``trialobj.data`` - annotated data object (based on AnnData library) for centralized storage of raw and processed data related to the experiment. contains methods for further modifying the stored annotated data object.




**AllOpticalTrial class**

``trialobj: AllOpticalTrial`` - packerlabimaging —> ``TwoPhotonImagingTrial`` —> ``AllOpticalTrial`` trial class


*Attributes:*

* this class is a child of the ``TwoPhotonImagingTrial`` (i.e. inherits all attributes and methods of ``TwoPhotonImagingTrial``)

    ``trialobj._2pstim`` - stores data and methods related to processing of 2p stim protocols

    ``trialobj.data`` - additional 2p stim related variables (e.g. ``photostim_frame`` and ``stim_start_frame`` ) to ``.data``

    ``OnePstimTrial class``



Submodules: Data Processing
---------------------------

``anndata.py``

``OnePstim.py``

``paq.py``

``stats.py``

``suite2p.py``

``TwoPstim.py``

- ``_readTargetsImage``
- ``_findTargetsAreas``
    - Loads target coordinates, and organizes target coordinates per SLM groups
    - creates target_areas - circle of pixels outward from the center of the target
    - creates target_areas_exclude - expanded circle of pixels outward from the center of the target that includes the region for excluding nontarget cells
    - creates images of targets // scaled somehow?? not totally sure - dont remember, has been commented out for a while now for myself - maybe check what Rob suggests adding here?



Submodules: Experiment Utils
----------------------------

``imagingMetadata.py``

``PrairieLink.py``

``STAMovieMaker_noGUI.py``

``ThorLink.py``

``utils.py``


Next
----

:ref:`Data structure details`

:ref:`tutorials`