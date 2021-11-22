from setuptools import setup, find_packages

setup(
    name='packerlabimaging_tempold',
    version='0.1',
    description='essential processing and analysis code for imaging experiments in Packer lab',
    url='http://github.com/{.....}',
    author='Packer Lab',
    author_email='adampacker@gmail.com',
    license='MIT',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    zip_safe=False
)
