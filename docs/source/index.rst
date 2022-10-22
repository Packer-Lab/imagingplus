.. imaging+ documentation master file, created by Prajay Shah
   sphinx-quickstart on Thu Feb 24 22:32:00 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Imaging+
========

Welcome to the imaging+ documentation!
++++++++++++++++++++++++++++++++++++++


This documentation contains guides and information about the usage and structure of ``imaging+``.
Please follow the instructions for installation and provided guides to get you started with utilizing the package.


Please open an issue or pull request on the project's `Github <https://github.com/Packer-Lab/imagingplus>`_ page if you are encountering any problems with the package.

Introduction
++++++++++++

``imaging+`` is an integrated tool-suite for essential processing, analysis and data visualization steps for multimodal calcium imaging+ neuroscience experiments in Python.


This tool-suite was first designed to serve as a scalable data analysis framework for all-optical experiments that combine 2-photon imaging, 2-photon optogenetic stimulation and multiple time-synced data signals (Packer et al., 2015, Nat Methods).
The ``imaging+`` framework is more broadly useful for multi-modal *imaging+* neuroscience experiments, in which calcium imaging data is collected with other time-synchronized signals.
Especially, there are methods for analysis and visualisation of standard calcium imaging data.

.. sidebar:: Multimodal imaging+ experiment:

     .. Figure:: files/Typical-experiment-apr-22-2022.png
         :width: 500


Currently, ``imaging+`` is most functional for experiments performed using `PackIO <http://apacker83.github.io/PackIO/>`_ for temporal synchronization, a Bruker microscope for 2p imaging and utilizing `Suite2p <https://suite2p.readthedocs.io/en/latest/>`_ for automated post-processing cell segmentation.

Follow :ref:`installation` for installation.
Please refer to :ref:`Quick start guide`,  :ref:`tutorials` and :ref:`overview` to get started and learn the toolbox.


Table of Contents
+++++++++++++++++

.. toctree::
   :maxdepth: 1

   Installation <installation.rst>
   Quick start guide <quick-start-guide.rst>
   Overview <overview.rst>
   Tutorials <Tutorials-reference.rst>
   Main Modules <organization-of-modules.rst>
   Data structure details <data-structure-details.rst>
   Making plots <making-plots.rst>
   API Reference <reference.rst>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
