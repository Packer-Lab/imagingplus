.. _Data structure:


Data Structure Inside ImagingTrial
==================================

*AnnData* sub-module for storing multi-modal data
---------------------------------------------

The AnnData library is used in the `anndata` submodule to store *imaging+*  data in a highly intuitive, efficient and multi-functional format.

The `anndata` table is used to store:

#. `.X`: ROI extracted imaging data in the main data matrix (dimensions: # of ROIs x # of frames)
#. `.var`: Multiple channels of 1-D temporal series data as synchronized variables (dimensions: # of channels x # of frames)
#. `.obs`: annotations of ROIs (dimensions # of ROIs x # of ROI annotations collected)
#. `.layers`: Further processing on the main data is added as additional layers (dimensions: # of ROIs x # of frames). For instance, dFF normalization of the raw data is added as the ``dFF`` layer to the existing data table.

.. image:: files/packerlabimaging-anndata-integration-01.jpg
    :width: 600

In keeping with AnnData conventions, the data structure is organized in *n* observations (obs) x *m* variables (var), where observations are suite2p ROIs and variables are imaging frame timepoints.
The rest of the AnnData data object is built according to these dimensions.

For instance, the metadata for each suite2p ROI stored in Suite2p’s stat.npy output is added to ``trialobject.data`` under ``obs`` and ``obsm`` (1D and >1-D observations annotations, respectively).
And the temporal synchronization data of the experiment collected in .paq output is added to the variables annotations under ``var``.


The primary benefit of anndata is that it enforces an intuitive data structure and allows a workflow that maintains relationships between different data sources and types, including processed data.
Lastly, AnnData objects are highly scalable, allowing individual users to further modify and add observations and variables of their own as their experiment dictates. Note: since AnnData is an independent library, they have an implementation for unstructured data annotations, however usage of this feature is not recommended within the packerlabimaging ecosystem as unstructured data elements should be built as attributes of the ``trialobject``.