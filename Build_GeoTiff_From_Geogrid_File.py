# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Copyright UCAR (c) 2019
# University Corporation for Atmospheric Research(UCAR)
# National Center for Atmospheric Research(NCAR)
# Research Applications Laboratory(RAL)
# P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# 24/09/2019
#
# Name:        Build_GeoTiff_From_Geogrid_File.py
# Purpose:
# Author:      Kevin Sampson, NCAR
# Created:     24/09/2019
# Licence:     Reserved
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

'''
This script will export a 2D gridded variable from a WRF-Hydro input file (Geogrid
or Fulldom_hires) to an output raster format, with all spatial and coordinate
system metadata.

If a 3-dimensional variable is selected, individual raster bands will be created
in the output raster for each item in the first dimension.

If a 4-dimensional variable is selected, the first item in the first dimension
will be selected and the variable will be treated as a 3-dimensional variable
described above.
'''

# Import Python Core Modules
import os
import time

# Import Additional Modules
import gdal
import numpy
import netCDF4
from distutils.version import LooseVersion

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import (WRF_Hydro_Grid, flip_grid, RasterDriver)
# --- Global Variables --- #

# Input file
in_nc = r"C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\geo_em.d03.20200327_snow.trim.nc"

# Variable information
Variable = 'HGT_M'

# Output directory
out_dir = r'C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\FOSS_Domain'

# Script options
out_Grid_fmt = RasterDriver                                                     # ['GTiff']

# --- End Global Variables --- #

# --- Main Codeblock --- #
if __name__ == '__main__':
    print('Script initiated at {0}'.format(time.ctime()))
    tic = time.time()

    rootgrp = netCDF4.Dataset(in_nc, 'r')                                       # Establish an object for reading the input NetCDF file
    grid_obj = WRF_Hydro_Grid(rootgrp)                                          # Instantiate a grid object
    print('    Created projection definition from input NetCDF GEOGRID file.')

    if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
        rootgrp.set_auto_mask(False)                                            # Change masked arrays to old default (numpy arrays always returned)

    # Convert 2D or 4D input variables to 3D
    variable = rootgrp.variables[Variable]
    dims = variable.dimensions
    array = variable[:]
    if len(dims) == 2:
        array = array[numpy.newaxis]
    elif len(dims) == 4:
        array = array[0]                                                        # Choose the first index, usually Time=0

    # Flip the grid to become north-to-south, which is the ordering used by GDAL to write rasters.
    if grid_obj.isGeogrid:
        array = flip_grid(array)

    # Build output rasters from numpy array of the GEOGRID variables requested
    OutRaster = grid_obj.numpy_to_Raster(array, nband=array.shape[0])
    rootgrp.close()
    del grid_obj, rootgrp, array

    # Build a geotiff using an input GEOGRID file and variable name
    OutGTiff = os.path.join(out_dir, '{0}.tif'.format(Variable))                # Output raster
    if OutRaster is not None:
        target_ds = gdal.GetDriverByName(out_Grid_fmt).CreateCopy(OutGTiff, OutRaster)
        print('    Created {0} raster: {1}'.format(Variable, OutGTiff))
        target_ds = None

    OutRaster = None
    del Variable, in_nc, out_dir, out_Grid_fmt, OutRaster, OutGTiff
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))