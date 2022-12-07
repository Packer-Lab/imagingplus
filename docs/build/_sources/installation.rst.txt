Installation instructions
=========================

Installation
------------
``imaging+`` can be installed by downloading or cloning the GitHub repository and running:

.. code-block:: python

    pip install -e imagingplus

on the command line from the directory where the GitHub repository was downloaded to.


*Note: The package is installable as a stand-alone python package. You can install the package into an existing conda environment, or you may choose to skip using conda environment all together.*

**To setup in a conda environment, we recommend following the steps below (on the command line):**

1. Clone the GitHub repository using Git::

    git clone https://github.com/Packer-Lab/imagingplus.git

2. Create a new conda environment using the ``plitest.yml`` provided in this repository by running::

    conda env create -f plitest.yml


3. Activate the conda environment::

    conda activate plitest

4. ``cd`` to the parent directory of where this repo was downloaded to.
5. From this parent directly, run::

    pip install -e imagingplus

    # this installs the package under developer settings (preferred for current release).


Test/Import Installation
------------------------
Test ``imaging+`` is successfully installed by importing the package in python:

1. Ensure that the conda environment from which `imaging+` was installed is activated.
2. start *python* on the command line using::

    python

or start a python session from the same conda environment in your preferred method (e.g. jupyter notebook or IDE).

3. Import the package::

    import imagingplus
    # or
    import imagingplus as ip

-- Ignore any warnings or printed messages :) --

Dependencies
------------
``imaging+`` requires the use of multiple scientific python packages (e.g ``numpy``, ``pandas``, ``anndata``, ``matplotlib``).
These dependencies should be automatically installed as required when installing ``imaging+``. It is common to have conflicts with dependencies that might require manually uninstalling and re-installing specified versions.
However, if there are non-trivial issues with dependencies during installation, please report this by raising an issue on Github.


Next
----

:ref:`overview`

:ref:`Quick start guide`