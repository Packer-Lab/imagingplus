from pathlib import Path
import os
from setuptools import setup, find_packages

test_deps = [
    "pytest~=7.0.1",
"pytest-cov==2.12.1",
"numpy~=1.19.5",
'pandas~=1.4.0',
'seaborn~=0.11.2',
'matplotlib~=3.5.1',
'packerlabimaging~=0.1',
'tifffile~=2022.2.9',
'scipy~=1.8.0',
'statsmodels~=0.13.2',
'h5py~=3.1.0',
'scikit-image~=0.19.1',
'opencv-python~=4.6.0',
'suite2p~=0.10.3',
'anndata~=0.7.8',
'setuptools~=58.0.4',
'tensorflow~=2.5.0',
'pystackreg~=0.2.5',
'deepinterpolation~=0.1.4']

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='imagingplus',
    version='0.2-beta',
    description='integrated tool-suite of essential processing and analysis of imaging+ experiments',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/Packer-Lab/imagingplus.git',
    author='Packer Lab',
    author_email='adampacker@gmail.com',
    license='MIT',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    install_requires=[
        'suite2p',
        'anndata',
        'tifffile',
        'numpy>=1.18.5',
        'seaborn',
        'matplotlib',
        'pandas',
        'scipy~=1.8.0',
        'statsmodels',
        'scikit-image',
        'opencv-python',
        'pystackreg',
        'mpl_point_clicker',
        'pynwb'
    ],
    tests_require=test_deps,
    zip_safe=False
)
