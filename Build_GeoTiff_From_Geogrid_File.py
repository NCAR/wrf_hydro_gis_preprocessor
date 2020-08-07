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

descText = 'This is a program to export >=2D variables from a WRF-Hydro input file ' + \
            '(geogrid or Fulldom_hires) file to an output raster format, with all ' + \
            'spatial and coordinate system metadata. If a 3-dimensional variable is ' + \
            'selected, individual raster bands will be created in the output raster ' + \
            'for each index in the 3rd dimension. If a 4-dimensional variable ' + \
            'is selected, the first index in the 4th dimension will be selected ' + \
            'and the variable will be treated as a 3-dimensional variable described above.'

'''
NOTES:
    Check for the position of certain dimensions?
'''

# Import Python Core Modules
import sys
import os
import time

# Import Additional Modules
import gdal
import numpy
import netCDF4
from distutils.version import LooseVersion
from argparse import ArgumentParser
from pathlib import Path

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import (WRF_Hydro_Grid, flip_grid, RasterDriver)

# --- Global Variables --- #

# Script options
out_Grid_fmt = RasterDriver                                                     # ['GTiff']
defaltGeogrid = 'geo_em.d01.nc'
version_number = '1.0'

# --- End Global Variables --- #

# --- Functions --- #
def build_geogrid_raster(in_nc, Variable, OutGTiff):

    rootgrp = netCDF4.Dataset(in_nc, 'r')                                       # Establish an object for reading the input NetCDF file
    grid_obj = WRF_Hydro_Grid(rootgrp)                                          # Instantiate a grid object
    print('    Created projection definition from input NetCDF GEOGRID file.')

    if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
        rootgrp.set_auto_mask(False)                                            # Change masked arrays to old default (numpy arrays always returned)

    # Convert 2D or 4D input variables to 3D
    if Variable in rootgrp.variables:
        variable = rootgrp.variables[Variable]
    else:
        print('Could not find variable {0} in input netCDF file.'.format(Variable))
        raise SystemExit

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
    if OutRaster is not None:
        target_ds = gdal.GetDriverByName(out_Grid_fmt).CreateCopy(OutGTiff, OutRaster)
        print('    Created {0} raster: {1}'.format(Variable, OutGTiff))
        target_ds = None
    OutRaster = None
    del Variable, in_nc, OutRaster, OutGTiff
# --- End Functions --- #

# --- Main Codeblock --- #
if __name__ == '__main__':

    print('Script initiated at {0}'.format(time.ctime()))
    tic = time.time()

    # Setup the input arguments
    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_nc",
                        default='./{0}'.format(defaltGeogrid),
                        help="Path to WPS geogrid (geo_em.d0*.nc) file or WRF-Hydro Fulldom_hires.nc file. default=./geo_em.d01.nc")
    parser.add_argument("-v",
                        dest="Variable",
                        default='HGT_M',
                        help="Name of the variable in the input netCDF file. default=HGT_M")
    parser.add_argument("-o",
                        dest="out_file",
                        default='./Output_GEOGRID_Raster.tif',
                        help="Output GeoTiff raster file.")

    # If no arguments are supplied, print help message
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    # Handle path of input
    if args.in_nc == all_defaults["in_nc"]:
        print('Using default input geogrid location of: {0}'.format(all_defaults["in_nc"]))
        in_nc = Path.cwd().joinpath(defaltGeogrid)
    else:
        in_nc = args.in_nc

    # Handle printing to user the default variable name
    if args.Variable == all_defaults["Variable"]:
        print('Using default variable name: {0}'.format(all_defaults["Variable"]))
    if args.out_file == all_defaults["out_file"]:
        print('Using default output location: {0}'.format(all_defaults["out_file"]))

    # Input WPS Namelist
    print('Input WPS Geogrid or Fulldom file: {0}'.format(in_nc))
    print('Input netCDF variable name: {0}'.format(args.Variable))
    print('Output raster file: {0}'.format(args.out_file))

    build_geogrid_raster(in_nc, args.Variable, args.out_file)
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))
# --- End Main Codeblock --- #