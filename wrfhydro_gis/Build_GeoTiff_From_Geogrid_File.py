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
import numpy
import netCDF4
from distutils.version import LooseVersion
from argparse import ArgumentParser
from pathlib import Path
import osgeo

try:
    if LooseVersion(osgeo.__version__) > LooseVersion('3.0.0'):
        from osgeo import gdal
    else:
        import gdal
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import (WRF_Hydro_Grid, RasterDriver, subset_ncVar)

# --- Global Variables --- #
out_fmt = RasterDriver                                                          # Could overwrite the output format. Default is 'GTiff'
defaltGeogrid = 'geo_em.d01.nc'                                                 # Default input geogrid file name if not provided by user
overwrite_output = True                                                         # Option to overwrite the output file if it exists already
# --- End Global Variables --- #

# --- Functions --- #
def build_geogrid_raster(in_nc, Variable, OutGTiff, out_Grid_fmt=out_fmt):
    '''
    Function to build a properly georeferenced raster object from a WRF-Hydro
    input file and variable name.
    '''

    # Check inputs for validity
    if os.path.exists(in_nc):
        rootgrp = netCDF4.Dataset(in_nc, 'r')                                   # Establish an object for reading the input NetCDF file
        grid_obj = WRF_Hydro_Grid(rootgrp)                                      # Instantiate a grid object

        # Change masked arrays to old default (numpy arrays always returned)
        if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
            rootgrp.set_auto_mask(False)
    else:
        print('    The input netCDF file does not exist: {0}'.format(in_nc))
        sys.exit(1)

    # Convert 2D or 4D input variables to 3D
    if Variable not in rootgrp.variables:
        print('    Could not find variable {0} in input netCDF file. Exiting...'.format(Variable))
        sys.exit(1)

    if os.path.exists(OutGTiff):
        if overwrite_output:
            print('    The output file already exists and will be overwritten: {0}'.format(OutGTiff))
        else:
            print('    The output file already exists. Exiting...')
            sys.exit(1)

    # Enhancement to allow n-dimensional arrays, with max dimension size of 3.
    if grid_obj.isGeogrid:
        array = subset_ncVar(rootgrp.variables[Variable], DimToFlip='south_north')
    else:
        array = subset_ncVar(rootgrp.variables[Variable], DimToFlip='')
    print('    Size of array being sent to raster: {0}'.format(array.shape))

    # Export numpy array to raster (up to 3D).
    OutRaster = grid_obj.numpy_to_Raster(array)
    print('    Bands in output raster: {0}'.format(OutRaster.RasterCount))
    rootgrp.close()
    del grid_obj, rootgrp, array

    # Save in-memory raster file to disk
    if OutRaster is not None:
        target_ds = gdal.GetDriverByName(out_Grid_fmt).CreateCopy(OutGTiff, OutRaster)
        print('    Created {0} format raster from {1} variable: {2}'.format(out_Grid_fmt, Variable, OutGTiff))
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
                        required=True,
                        help="Path to WPS geogrid (geo_em.d0*.nc) file or WRF-Hydro Fulldom_hires.nc file.")
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

    # Handle printing to user the default variable name
    if args.Variable == all_defaults["Variable"]:
        print('Using default variable name: {0}'.format(all_defaults["Variable"]))
    if args.out_file == all_defaults["out_file"]:
        print('Using default output location: {0}'.format(all_defaults["out_file"]))

    # Print information to screen
    print('Input WPS Geogrid or Fulldom file: {0}'.format(args.in_nc))
    print('Input netCDF variable name: {0}'.format(args.Variable))
    print('Output raster file: {0}'.format(args.out_file))

    build_geogrid_raster(args.in_nc, args.Variable, args.out_file)
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))
# --- End Main Codeblock --- #