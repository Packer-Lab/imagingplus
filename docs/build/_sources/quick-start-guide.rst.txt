.. _Quick start guide:

Quick Start Guide
=================

The easiest way to get started with using ``imaging+`` is to follow `Tutorial 1 - Getting started with data`_.
This tutorial demonstrates how to setup a new experiment inside ``imaging+`` for data processing and analysis.

To get started with plotting and data exploration, check out `Tutorial 3 - Making useful plots`_.


Key Structures
++++++++++++++

There are two core data objects for storing experimental data:


1) `Experiment` - built from constituent `ImagingTrial` objects, and represents the overall experiment.

2) `ImagingTrial` - a calcium imaging trial at a single field of view.



And, there 3 sub-core data types to handle multi-modal data within each imaging trial:

1) `ImagingData` - live imaging data

2) `TemporalData` - 1D time series data temporally synchronized to live imaging data

3) `CellAnnotations` - Annotations and cellular/ROI segmentations of the live imaging tissue sample


A more detailed overview of these modules and other key structures of ``imaging+`` can be found at :ref:`overview`.


An overview of all primary modules of ``imaging+`` can be found at :ref:`main modules`.


.. _Tutorial 1 - Getting started with data: Tutorials/Tutorial-1-Overview-of-Workflow.ipynb
.. _Tutorial 3 - Making useful plots: Tutorials/Tutorial-3-plotting-module.ipynb



Next
----
:ref:`overview`

:ref:`Data structure details`

:ref:`tutorials`

:ref:`Making plots`



