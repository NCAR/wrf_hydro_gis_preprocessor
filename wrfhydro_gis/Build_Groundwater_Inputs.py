# # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Copyright UCAR (c) 2019
# University Corporation for Atmospheric Research(UCAR)
# National Center for Atmospheric Research(NCAR)
# Research Applications Laboratory(RAL)
# P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# 24/09/2019
#
# Name:        module1
# Purpose:
# Author:      $ Kevin Sampson
# Created:     24/09/2019
# Licence:     <your licence>
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

descText = "This script will build WRF-Hydro groundwater basin inputs using an input WPS Geogrid " \
           "file and WRF-Hydro routing grid (Fulldom_hires.nc) file as inputs. Three methods" \
           "are currently available for generating groundwater basins. One method, 'Polygon" \
           "Shapefile or Feature Class' currently requires an input polygon shapefile defining" \
           "the groundwater basins."

# Import Python Core Modules
import os
import time
import shutil
import sys
import copy
from distutils.version import LooseVersion

# Import additional modules
import netCDF4
import gdal
from argparse import ArgumentParser
from pathlib import Path

# Import function library into namespace. Must exist in same directory as this script.
#import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox
from wrfhydro_functions import (WRF_Hydro_Grid, GW_nc, GWGRID_nc, dir_d8, streams,
    basinRaster, RasterDriver, build_GW_Basin_Raster, build_GW_buckets, remove_file,
    zipUpFolder)
#import Examine_Outputs_of_GIS_Preprocessor as EO

# --- Global Variables --- #

# --- EDIT BELOW THIS LINE --- #

# Provide the default groundwater basin generation method.
defaultGWmethod = 'FullDom basn_msk variable'
defaultFulldom = 'Fulldom_hires.nc'
defaultGeogrid = 'geo_em.d01.nc'

# --- EDIT ABOVE THIS LINE --- #

# --- DO NOT EDIT BELOW THIS LINE --- #

# List of routing-stack files to send to output .zip files
nclist = [GW_nc, GWGRID_nc]

# Save a GeoTiff of the location of the derived basins on the fine and coarse grids
saveBasins_Fine = False
saveBasins_Coarse = False

# --- End Global Variables --- #

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()
    print('Script initiated at {0}'.format(time.ctime()))

    # Setup the input arguments
    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_nc",
                        help="Path to WPS geogrid (geo_em.d0*.nc) file or WRF-Hydro Fulldom_hires.nc file.")
    parser.add_argument("-f",
                        dest="in_fulldom",
                        default='./{0}'.format(defaultFulldom),
                        help="Path to WRF-Hydro Fulldom_hires.nc file.")
    parser.add_argument("-m",
                        dest="GWmethod",
                        default='FullDom LINKID local basins',
                        help="Method to create groundwater basins. Choose from 'FullDom basn_msk variable', "
                             "'FullDom LINKID local basins', 'Polygon Shapefile or Feature Class'"
                             " default='FullDom basn_msk variable'")
    parser.add_argument("-g",
                        dest="in_GWPolys",
                        default='',
                        help="Path to groundwater basin polygon or feature class file.")
    parser.add_argument("-o",
                        dest="out_dir",
                        default='',
                        required=True,
                        help="Output directory.")

    # If no arguments are supplied, print help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    # Handle path of input
    if args.in_nc == all_defaults["in_nc"]:
        print('Using default input geogrid location of: {0}'.format(all_defaults["in_nc"]))

    if args.in_fulldom == all_defaults["in_fulldom"]:
        print('Using default input fulldom location of: {0}'.format(all_defaults["in_fulldom"]))

    if args.GWmethod == all_defaults["GWmethod"]:
        print('Using groundwater method of: {0}'.format(all_defaults["GWmethod"]))

    if args.out_dir == all_defaults["out_dir"]:
        print('Using default output location: {0}'.format(all_defaults["out_file"]))

    # Outputs - permanent
    out_zip = os.path.join(args.out_dir, 'GroundwaterBasins_local.zip')

    # Create scratch directory for temporary outputs
    projdir = os.path.join(args.out_dir, 'gw_scratchdir')
    projdir = os.path.abspath(projdir)
    if os.path.exists(projdir):
        shutil.rmtree(projdir)
    os.makedirs(projdir)

    # Setup temporary output files
    fdir = os.path.join(projdir, dir_d8)
    channelgrid = os.path.join(projdir, streams)
    if saveBasins_Fine:
        basinRaster_File = os.path.join(projdir, 'GWBasins_fine.tif')
        nclist.append('GWBasins_fine.tif')
    if saveBasins_Coarse:
        nclist.append(basinRaster)

    rootgrp1 = netCDF4.Dataset(args.in_nc, 'r')
    rootgrp2 = netCDF4.Dataset(args.in_fulldom, 'r')
    if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
        rootgrp2.set_auto_mask(False)                                            # Change masked arrays to old default (numpy arrays always returned)
    coarse_grid = WRF_Hydro_Grid(rootgrp1)                                 # Instantiate a grid object for the coarse grid
    fine_grid = WRF_Hydro_Grid(rootgrp2)                                   # Instantiate a grid object for the fine grid
    #fine_grid = copy.copy(coarse_grid)                                          # Copy the grid object for modification
    #fine_grid.regrid(4)                                                         # Regrid to the coarse grid

    # Build inputs required for creating groundwater buckets
    flowdir = fine_grid.numpy_to_Raster(rootgrp2.variables['FLOWDIRECTION'][:])
    strm_arr = rootgrp2.variables['CHANNELGRID'][:]
    strm_arr[strm_arr==0] = 1                                                   # Set active channels to 1
    strm = fine_grid.numpy_to_Raster(strm_arr)
    del strm_arr

    # Save to disk for the Groundwater tools to use
    out_ds1 = gdal.GetDriverByName(RasterDriver).CreateCopy(fdir, flowdir)
    out_ds2 = gdal.GetDriverByName(RasterDriver).CreateCopy(channelgrid, strm)
    out_ds1 = out_ds2 = flowdir = strm = None

    # Build groundwater files
    print('  Building Groundwater Basin inputs.')
    GWBasns = build_GW_Basin_Raster(args.in_fulldom, projdir, args.GWmethod, channelgrid, fdir, fine_grid, in_Polys=args.in_GWPolys)
    build_GW_buckets(projdir, GWBasns, coarse_grid, Grid=True, saveRaster=saveBasins_Coarse)
    if saveBasins_Fine:
        out_ds3 = gdal.GetDriverByName(RasterDriver).CreateCopy(basinRaster_File, GWBasns)
        out_ds3 = None
    GWBasns = None
    remove_file(fdir)
    remove_file(channelgrid)
    del GWBasns, coarse_grid, fine_grid

    # zip the folder
    tic1 = time.time()
    zipper = zipUpFolder(projdir, out_zip, nclist)
    print('Built output .zip file: {0}'.format(out_zip))

    # Delete all temporary files
    shutil.rmtree(projdir)
    rootgrp1.close()
    rootgrp2.close()
    del rootgrp1, rootgrp2
    print('Process completed in {0:3.2f} seconds.'.format(time.time()-tic))