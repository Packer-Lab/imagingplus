.. _Quick start guide:

Quick Start Guide
=================

The easiest way to get started with using ``packerlabimaging`` is to follow `Tutorial 1 - Initializing an Experiment`_.
This tutorial demonstrates how to setup a new experiment inside ``packerlabimaging`` for data processing and analysis.

To get started with plotting and data exploration, check out `Tutorial 3 - Making useful plots`_.


Key Structures
==============

There are two core data objects for storing experimental data:


1) `Experiment` and

2) `ImagingTrial` objects

- The `Experiment` class is built from information provided for each constituent trial of the overall Experiment.



And, there 3 sub-core data types to handle multi-modal data within each imaging trial:

1) `ImagingData` - live imaging data

2) `TemporalData` - 1D time series data temporally synchronized to live imaging data

3) `CellAnnotations` - Annotations and cellular/ROI segmentations of the live imaging tissue sample


An more detailed overview of these modules and other key structures of ``packerlabimaging`` can be found at :ref:`overview`.




.. _Tutorial 1 - Initializing an Experiment: Tutorials/Tutorial-1-Initializing-an-Experiment.ipynb
.. _Tutorial 3 - Making useful plots: Tutorials/Tutorial-3-Making-useful-plots.ipynb



Next
----
:ref:`overview`

:ref:`Data structure details`

:ref:`tutorials`

:ref:`Making plots`



