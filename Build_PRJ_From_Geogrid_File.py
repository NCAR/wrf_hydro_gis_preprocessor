# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
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

descText = "This tool takes an input WRF Geogrid file in NetCDF format and uses the" \
           " specified variable's projection parameters to produce a projection file."

# Import Modules

# Import Python Core Modules
import os
import time
import sys

# Import additional modules
import netCDF4
from argparse import ArgumentParser
from pathlib import Path

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import WRF_Hydro_Grid
print('Script initiated at {0}'.format(time.ctime()))

# Global Variables
defaultGeogrid = 'geo_em.d01.nc'

# Script options
buildPRJ = True                                                                 # Switch for building output Esri Projection File

# Main Codeblock
if __name__ == '__main__':
    print('Script initiated at {0}'.format(time.ctime()))
    tic = time.time()

    # Setup the input arguments
    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_nc",
                        required=True,
                        help="Path to WPS geogrid (geo_em.d0*.nc) file or WRF-Hydro Fulldom_hires.nc file.")
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

    if args.out_dir == all_defaults["out_dir"]:
        print('Using default output location: {0}'.format(all_defaults["out_file"]))

    # Input and output files and directories
    outPRJ = os.path.join(args.out_dir, os.path.basename(args.in_nc).replace('.nc', '.prj'))

    print('Input WPS Geogrid or Fulldom file: {0}'.format(args.in_nc))
    print('Output prj file: {0}'.format(outPRJ))

    rootgrp = netCDF4.Dataset(args.in_nc, 'r')                                   # Establish an object for reading the input NetCDF file
    coarse_grid = WRF_Hydro_Grid(rootgrp)                                  # Instantiate a grid object
    rootgrp.close()
    del rootgrp
    print('    Created projection definition from input NetCDF GEOGRID file.')

    # Build an ESRI style projection file
    if buildPRJ:
        projEsri = coarse_grid.proj.Clone()                                     # Copy the SRS
        projEsri.MorphToESRI()                                                  # Alter the projection to Esri's representation of a coordinate system
        file = open(outPRJ, 'w')
        file.write(projEsri.ExportToWkt())
        file.close()
        print('    Created ESRI Projection file: {0}'.format(outPRJ))
    del coarse_grid
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))
