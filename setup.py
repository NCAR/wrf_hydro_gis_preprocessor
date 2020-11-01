from setuptools import find_packages, setup

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='wrfhydro_gis',
    version='0.1.0',
    packages=find_packages(include=['wrfhydro_gis']),
    url='https://github.com/NCAR/wrf_hydro_gis_preprocessor',
    license='MIT',
    install_requires=[
        'gdal',
        'netcdf4==1.5.3',
        'numpy==1.18.1',
        'pyproj==2.6.0',
        'whitebox==1.2.0'
    ],
    author='Kevin Sampson & Matt Casali',
    author_email='ksampson@ucar.edu',
    description='Geospatial Pre-processing functions for the WRF-Hydro model',
    python_requires='>=3.6.10',
)




