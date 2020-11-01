# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Copyright UCAR (c) 2020
# University Corporation for Atmospheric Research(UCAR)
# National Center for Atmospheric Research(NCAR)
# Research Applications Laboratory(RAL)
# P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
#
# Name:        Testing_DEM_interpolation.py
# Purpose:
# Author:      $ Kevin Sampson(ksampson)
# Created:     2020
# Licence:     <your licence>
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

descText = "This tool takes a WRF/WPS Geogrid file and a higher-resolution DEM raster " +\
            "and interpolates the file to a grid that is either identical or nested into " +\
            "the geogrid file Mass grid (Stagger='M'). This tool may be used to test " +\
            "the effect of different regridding factors and interpolation methods on the " +\
            "output DEM. It may also be used to interpolate other raster files to the model "  +\
            "grid, but caution should be used with regard to interpolating discrete data." +\
            "Default interpolation method is bilinear."

# --- Import Modules --- #

# Import Python core modules
import sys
sys.dont_write_bytecode = True
import time
import os
import copy
from distutils.version import LooseVersion
from argparse import ArgumentParser

# Import Additional Modules
import gdal
import netCDF4

# Import function library into namespace. Must exist in same directory as this script.
#wrfhGIS_lib = r''
sys.path.append(wrfhGIS_lib)
import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox

# Module options
gdal.UseExceptions()                                                            # this allows GDAL to throw Python Exceptions
gdal.PushErrorHandler('CPLQuietErrorHandler')

# --- End Import Modules --- #

# --- Global Variables --- #

# Default values
resampling_method = gdal.GRA_Bilinear                                           # Default regridding method
default_regridFactor = 1                                                        # Default regridding factor

# --- End Global Variables --- #

# --- Functions --- #
def is_valid_file(parser, arg):
    # https://stackoverflow.com/questions/11540854/file-as-command-line-argument-for-argparse-error-message-if-argument-is-not-va
    if not os.path.exists(arg):
        parser.error("The file {0} does not exist!".format(arg))
    else:
        return str(arg)

def interpolate_raster(in_Geogrid, inDEM, cellsize, out_file):
    tic1 = time.time()

    # Georeference geogrid file
    rootgrp = netCDF4.Dataset(in_Geogrid, 'r')                             # Establish an object for reading the input NetCDF files
    globalAtts = rootgrp.__dict__                                               # Read all global attributes into a dictionary
    if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
        rootgrp.set_auto_mask(False)                                            # Change masked arrays to old default (numpy arrays always returned)

    # Build grid object
    coarse_grid = wrfh.WRF_Hydro_Grid(rootgrp)                                  # Instantiate a grid object
    fine_grid = copy.copy(coarse_grid)                                          # Copy the grid object for modification
    fine_grid.regrid(cellsize)                                             # Regrid to the desired nest ratio
    print('    Created projection definition from input NetCDF GEOGRID file.')
    print('    Proj4: {0}'.format(coarse_grid.proj4))                           # Print Proj.4 string to screen
    print('    Coarse grid GeoTransform: {0}'.format(coarse_grid.GeoTransformStr()))                    # Print affine transformation to screen.
    print('    Coarse grid extent [Xmin, Ymin, Xmax, Ymax]: {0}'.format(coarse_grid.grid_extent()))     # Print extent to screen.
    print('    Fine grid extent [Xmin, Ymin, Xmax, Ymax]:   {0}'.format(fine_grid.grid_extent()))       # Print extent to screen.

    # Create high resolution topography layers
    in_DEM = gdal.Open(inDEM, 0)                                           # Open with read-only mode
    mosprj = fine_grid.project_to_model_grid(in_DEM, saveRaster=True, OutGTiff=out_file, resampling=resampling_method)
    in_DEM = mosprj = None
    print('  Output file created: {0}'.format(out_file))
    print('  Interpolation completed in {0:3.2f} seconds.'.format(time.time()-tic1))

# --- End Functions --- #

# --- Main Codeblock --- #
if __name__ == '__main__':
    print('Script initiated at {0}'.format(time.ctime()))
    tic = time.time()

    # Setup the input arguments
    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_Geogrid",
                        required=True,
                        help="Path to WPS geogrid (geo_em.d0*.nc) file or WRF-Hydro Fulldom_hires.nc file.")
    parser.add_argument("-d",
                        dest="inDEM",
                        type=lambda x: is_valid_file(parser, x),
                        default='',
                        required=True,
                        help="Path to input high-resolution elevation raster [REQUIRED]")
    parser.add_argument("-R",
                        dest="cellsize",
                        type=int,
                        default=default_regridFactor,
                        help="Regridding (nest) Factor. default=10")
    parser.add_argument("-o",
                        dest="out_file",
                        default='./',
                        required=True,
                        help="Output raster file (.tif).")

    # If no arguments are supplied, print help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    # Handle printing to user the default variable name
    print('  Parameter values that have not been altered from script default values:')
    if args.cellsize == all_defaults["cellsize"]:
        print('Using default regridding factor of: {0}'.format(all_defaults["cellsize"]))
    if args.out_file == all_defaults["out_file"]:
        print('Using default output location of: {0}'.format(all_defaults["out_file"]))

    # Print information to screen
    print('  Values that will be used in building this routing stack:')
    print('    Input WPS Geogrid file: {0}'.format(args.in_Geogrid))
    print('    Input high-resolution DEM: {0}'.format(args.inDEM))
    print('    Regridding factor: {0}'.format(args.cellsize))
    print('    Output raster file: {0}'.format(args.out_file))
    print('    Resampling method: {0}\n'.format(str(resampling_method)))

    # Run the function to interpolate the raster, given inputs
    interpolate_raster(args.in_Geogrid, args.inDEM, args.cellsize, args.out_file)
    print('Process completed in {0:3.2f} seconds.'.format(time.time()-tic))