.. packerlabimaging documentation master file, created by Prajay Shah
   sphinx-quickstart on Thu Feb 24 22:32:00 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

packerlabimaging
================

Welcome to packerlabimaging's documentation!
++++++++++++++++++++++++++++++++++++++++++++


This documentation contains guides and information about the usage and structure of ``packerlabimaging``.
Please follow the instructions for installation and provided guides to get you started with utilizing the package.


Please open an issue or pull request on the project's `Github <https://github.com/Packer-Lab/packerlabimaging>`_ page if you are encountering any problems with the package.

Introduction
++++++++++++

``packerlabimaging`` is an a integrated tool-suite for essential processing, analysis and data visualization steps for multimodal calcium imaging+ neuroscience experiments in Python.


In particular, ``packerlabimaging`` was first designed to serve a structured data analysis framework for all-optical experiments, which combine 2-photon imaging, 2-photon stimulation and multiple other time-synced data signals (Packer et al., 2015, Nat Methods).
However, the ``packerlabimaging`` framework is more broadly useful for multi-modal imaging+ neuroscience experiments, in which calcium imaging data is collected with other time-synchronized signals.
Especially, there are methods for analysis and visualisation of standard calcium imaging data.

.. sidebar:: Multimodal imaging+ experiment:

     .. Figure:: files/Typical-experiment-apr-22-2022.jpeg
         :width: 600



Currently, ``packerlabimaging`` is most functional for experiments performed using PackIO (REF) for temporal synchronization, a Bruker microscope for 2p imaging and utilizing Suite2p (REF) for cell segmentation.

Follow :ref:`installation` for installation.
Please refer to :ref:`overview` and :ref:`tutorials` to get started and learn the toolbox.


Table of Contents
+++++++++++++++++

.. toctree::
   :maxdepth: 1

   Installation <installation.rst>
   Quick start guide <quick-start-guide.rst>
   Overview <overview.rst>
   Tutorials <Tutorials-reference.rst>
   Main Modules <main-modules.rst>
   Data structure details <data-structure-details.rst>
   Making plots <making-plots.rst>
   API Reference <reference.rst>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
