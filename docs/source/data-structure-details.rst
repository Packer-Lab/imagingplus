.. _Data structure details:


Data Structure Inside Trial objects
===================================

Annotated data object for storing multi-functional data
-------------------------------------------------------

The AnnData library is used to store data in an efficient, multi-functional format. This is stored under: ``trialobject.data``.
The AnnData object is built around the raw Flu matrix of each ``trialobject``.
In keeping with AnnData conventions, the data structure is organized in *n* observations (obs) x *m* variables (var), where observations are suite2p ROIs and variables are imaging frame timepoints.
The rest of the AnnData data object is built according to these dimensions.
For instance, the metadata for each suite2p ROI stored in Suite2p’s stat.npy output is added to ``trialobject.data`` under ``obs`` and ``obsm`` (1D and >1-D observations annotations, respectively).
And the temporal synchronization data of the experiment collected in .paq output is added to the variables annotations under ``var``.

.. image:: files/packerlabimaging-anndata-integration-01.jpg
    :width: 600


Further processing on the raw data is added to the AnnData object as ``layers``. For instance, dFF normalization of the raw data is added as the ``dFF`` layer to the existing AnnData object.

The primary benefit of anndata is that it enforces an intuitive data structure and allows a workflow that maintains relationships between different data sources and types, including processed data.
Lastly, AnnData objects are highly scalable, allowing individual users to further modify and add observations and variables of their own as their experiment dictates. Note: since AnnData is an independent library, they have an implementation for unstructured data annotations, however usage of this feature is not recommended within the packerlabimaging ecosystem as unstructured data elements should be built as attributes of the ``trialobject``.