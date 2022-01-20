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

descText = "This tool takes an WRF Geogrid file and creates a single polygon shapefile" \
           " that makes up the boundary of the domain of the M-grid (HGT_M, for example)."

# Import Modules

# Import Python Core Modules
import os
import sys
import time

# Import additional modules
import netCDF4
from argparse import ArgumentParser
from pathlib import Path
from distutils.version import LooseVersion
import osgeo

try:
    if LooseVersion(osgeo.__version__) > LooseVersion('3.0.0'):
        from osgeo import ogr
    else:
        import ogr
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import WRF_Hydro_Grid                                   # Function script packaged with this toolbox

# --- Global Variables --- #
outDriverName = 'ESRI Shapefile'                                                # Output vector file format (OGR driver name)
outSHPDefault = 'domain_boundary.shp'                                           # Default output filename if none is provided

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()
    print('Script initiated at {0}'.format(time.ctime()))

    # Setup the input arguments
    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_nc",
                        required=True,
                        help="Path to WPS geogrid (geo_em.d0*.nc) file or WRF-Hydro Fulldom_hires.nc file.")
    parser.add_argument("-o",
                        dest="out_dir",
                        default='./{0}'.format(outSHPDefault),
                        required=True,
                        help="Output directory.")

    # If no arguments are supplied, print help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    # Handle output path
    if args.out_dir == all_defaults["out_dir"]:
        print('Using default output location of: {0}'.format(all_defaults["out_dir"]))

    # Input and output files and directories
    outSHP = os.path.join(args.out_dir, os.path.basename(args.in_nc).replace('.nc', '_boundary.shp'))

    rootgrp = netCDF4.Dataset(args.in_nc, 'r')                                   # Establish an object for reading the input NetCDF file
    coarse_grid = WRF_Hydro_Grid(rootgrp)                                       # Instantiate a grid object
    print('    Created projection definition from input NetCDF GEOGRID file.')

    # Write out domain boundary shapefile
    '''For some unkonwn reason, the OGR.createlayer method will not pass CRS scale factors
    into the output for polar sterographic projections.'''
    geom = coarse_grid.boundarySHP(outSHP, outDriverName)
    geom = None
    rootgrp.close()
    del rootgrp, geom
    print('  Output shapefile: {0}'.format(outSHP))
    print('Process completed in {0:3.2f} seconds.'.format(time.time()-tic))