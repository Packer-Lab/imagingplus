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

packerlabimaging is a Python package for essential processing and analysis of 2photon imaging data collected in the Packer Lab.
Especially, there is a fully implemented pipeline for data structuring, processing, analysis and plotting of 2photon Ca2+ imaging experiments, and experiments based around 2photon imaging such as all optical experiments (i.e. 2photon optogenetic stim or 1photon optogenetic stim, with combined 2photon Ca2+ imaging).


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
   Data Structure details <data-structure-details.rst>
   Making plots <making-plots.rst>
   API Reference <reference.rst>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
