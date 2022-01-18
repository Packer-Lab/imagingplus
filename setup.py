from setuptools import setup, find_packages

setup(
    name='packerlabimaging',
    version='0.1',
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
        'numpy',
        'seaborn',
        'suite2p',
        'matplotlib',
        'pandas',
        'scipy',
        'statsmodels',
        'scikit-image'
    ],
    zip_safe=False
)
