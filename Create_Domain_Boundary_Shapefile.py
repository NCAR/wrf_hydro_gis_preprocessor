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

# Import Modules

# Import Python Core Modules
import os
import time

# Import additional modules
import netCDF4
import ogr

# Import function library into namespace. Must exist in same directory as this script.
#import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox
from wrfhydro_functions import WRF_Hydro_Grid
print('Script initiated at {0}'.format(time.ctime()))

# Global Variables

# Input and output files and directories
inGeogrid = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\geo_em.d01.nc'
out_dir = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs'
outSHP = os.path.join(out_dir, os.path.basename(inGeogrid).replace('.nc', '_boundary.shp'))

# Script options
outDriverName = 'ESRI Shapefile'                                                # Output vector file format (OGR driver name)

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()

    rootgrp = netCDF4.Dataset(inGeogrid, 'r')                                   # Establish an object for reading the input NetCDF file
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