from pathlib import Path
import os
from setuptools import setup, find_packages

TESTS_REQUIRE = Path(f"{os.getcwd()}/requirements_dev.txt").read_text().splitlines()


# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='imagingplus',
    version='0.1-beta',
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
