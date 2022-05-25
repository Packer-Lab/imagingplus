.. packerlabimaging documentation master file, created by
   sphinx-quickstart on Thu Feb 24 22:32:00 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

packerlabimaging
================

Welcome to packerlabimaging's documentation!
++++++++++++++++++++++++++++++++++++++++++++


This documentation contains detailed information about the code base of ``packerlabimaging``.
There are instructions for installation and guides to get you started with utilizing the package.


Please open an issue or pull request on the project's `Github <https://github.com/Packer-Lab/packerlabimaging>`_ page if you are encountering any problems with the package.

Introduction
++++++++++++

`packerlabimaging` is a Python package for essential processing and analysis of 2photon imaging data collected in the Packer Lab.

`packerlabimaging` was first designed to provide a toolbox of data processing methods and a structured data analysis framework for all-optical experiments, which combine 2-photon imaging, 2-photon stimulation and multiple other time-synced data streams (Packer et al., 2015, Nat Meth). However, the packerlabimaging toolbox can be more broadly applied to other configurations of imaging+ neuroscience experiments, such as any combination of imaging + behaviour tracking + sensory stimulation + electrophysiology, etc.
Especially, there are workflows for analysis and visualisation of standard Ca2+ imaging data, and of standard all-optical experiments. Out-of-the-box, the package is designed to work best with imaging experiments performed using a Bruker 2pPlus microscope system, PACKIO (REF) for experimental hardware temporal synchronisation and utilises Suite2p (REF) for processing of Ca2+ imaging data.
Packerlabimaging implements a structure of data storage and processing that is designed to be modular, highly intuitive and readily extendable to serve the end-user’s unique analytic needs. This framework gives the user a single framework to integrate the data processing and analyses needs of a variety of multi-modal experiments.



.. sidebar:: Schematic of overall package flow

     .. Figure:: files/OverallPackageFlowDiagram.png
         :width: 600



This package is designed to work best with imaging experiments performed using PackIO, a Bruker 2pPlus microscope and Suite2p for Ca2+ imaging data processing for ROI detection.
Ultimately, the goal of this package is to jump-start your own analysis of your awesome experiment.
It should be completely usable and understandable for anyone with the correct data in hand and basic knowledge of Python. There are :ref:`tutorials` to help along the way.
We hope that it provides a useful structure to organize your experimental data, and some functionality to interact and process your data in an efficient manner.

Table of Contents
+++++++++++++++++

.. toctree::
   :maxdepth: 1

   Installation <installation.rst>
   Quick start guide <quick-start-guide.rst>
   Overview <overview.rst>
   Tutorials <Tutorials-reference.rst>
   Main Modules <main-modules.rst>
   Data Structure details <data-structure-details.rst>
   Making plots <making-plots.rst>
   API Reference <reference.rst>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
