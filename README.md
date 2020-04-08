# WRF-Hydro GIS Preprocessor

## Description:
The WRF-Hydro GIS Pre-processor provides various scripts and tools for building the geospatial input files for running a WRF-Hydro simulation.

## Necessary Python Packages and Installation Tips

The scripts for the WRF-Hydro GIS Pre-processor rely on several python modules a user will need to install such as numpy, gdal, and Whitebox-Tools.  We highly recommend using a python distribution such as Miniconda(https://docs.conda.io/en/latest/miniconda.html). We used the conda package manager to install necessary python modules. Our essential packages and versions used are listed below (Windows 64-bit and Python 3.6.10):

| Package       | Version       | 
| ------------- |--------------:| 
| gdal          | 3.0.4         | 
| netcdf4       | 1.5.3         |
| numpy         | 1.18.1        |
| pyproj        | 2.6.0         |
| python        | 3.6.10        |
| whitebox      | 1.2.0         |

If you are using Anaconda, creating a new, clean 'wrfh_gis' environment with these needed packages can be done easily and simply one of several ways:

* In your conda shell, add one necessary channel and then download the wrfh_gis environment from the Anaconda cloud:
  + conda config --add channels conda-forge
  + conda create -n wrfh_gis2 -c conda-forge python=3.6 gdal netCDF4 numpy pyproj whitebox=1.2.0
  
* Alternatively, using the wrfh_gis.yml file in this repository, in your conda shell cd to the directory where your wrfh_gis.yml file is located and run:
  + conda env create -f environment.yml

* To activate this new environment, type the following at the conda prompt
  + activate wrfh_gis
  
## How to Run Scripts 

### The scripts make use of a function script 'wrfhydro_functions.py' to pass all functions and selected global parameters parameters to the primary script: 

+ [Build_Routing_Stack.py](https://github.com/NCAR/wrf_hydro_gis_preprocessor/blob/master/Build_Routing_Stack.py).

In turn, these scripts rely on a set of functions in [wrfhydro_functions.py](https://github.com/NCAR/wrf_hydro_gis_preprocessor/blob/master/wrfhydro_functions.py). 

### Running Build_Routing_Stack.py to generate the routing grids for a new WRF-Hydro simulation domain

After editing the globals in both [wrfhydro_functions.py] and [Build_Routing_Stack.py] for the desired result, such as the threshold to define a channel cell, the following steps will excecute processes to generate new routing grids for the desired domain. This example assumes use of a Bash shell.
1.

## NCAR Disclaimer
The National Center for Atmospheric Research (NCAR) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use.  NCAR has relinquished control of the information and no longer has responsibility to protect the integrity , confidentiality, or availability of the information.  Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by NCAR.  The NCAR seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by NCAR or the National Science Foundation (NSF).