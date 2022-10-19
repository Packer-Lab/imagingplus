from pathlib import Path
import os
from setuptools import setup, find_packages

TESTS_REQUIRE = Path(f"{os.getcwd()}/requirements_dev.txt").read_text().splitlines()

setup(
    name='packerlabimaging',
    version='0.1-beta',
    description='essential processing and analysis code for imaging experiments in Packer lab',
    url='https://github.com/Packer-Lab/packerlabimaging.git',
    author='Packer Lab',
    author_email='adampacker@gmail.com',
    license='MIT',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    install_requires=[
        'anndata',
        'tifffile',
        'numpy>=1.18.5',
        'seaborn',
        'matplotlib',
        'pandas',
        'scipy',
        'statsmodels',
        'scikit-image',
        'suite2p',
        'opencv-python',
        'pystackreg',
        'mpl_point_clicker',
        'pynwb'
    ],
    tests_require=TESTS_REQUIRE,
    zip_safe=False
)
