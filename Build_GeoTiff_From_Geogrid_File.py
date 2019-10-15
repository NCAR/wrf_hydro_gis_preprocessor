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

# Import Additional Modules
#import osr
import gdal
import netCDF4

# Import function library into namespace. Must exist in same directory as this script.
import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox

print('Script initiated at {0}'.format(time.ctime()))

# Global Variables

# Input and output files and directories
inGeogrid = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\geo_em.d01.nc'
out_dir = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs'
Variable = 'HGT_M'
OutGTiff = os.path.join(out_dir, '{0}.tif'.format(Variable))       # Output raster

# Variables derived from function script
out_Grid_fmt = wrfh.RasterDriver

# Script options
buildGeogridRaster = True                                                       # Switch for building output GeoTiff from GEOGRID file

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()

    rootgrp = netCDF4.Dataset(inGeogrid, 'r')                                       # Establish an object for reading the input NetCDF file
    coarse_grid = wrfh.WRF_Hydro_Grid(rootgrp)                              # Instantiate a grid object
    print('    Created projection definition from input NetCDF GEOGRID file.')

    # Build output rasters from numpy array of the GEOGRID variables requested
    OutRaster = wrfh.numpy_to_Raster(wrfh.flip_grid(rootgrp.variables[Variable][0]),
                                        coarse_grid.proj,
                                        coarse_grid.DX,
                                        coarse_grid.DY,
                                        coarse_grid.x00,
                                        coarse_grid.y00)
    rootgrp.close()
    del coarse_grid, rootgrp

    # Build a geotiff using an input GEOGRID file and variable name
    if buildGeogridRaster:

        if OutRaster is not None:
            target_ds = gdal.GetDriverByName(out_Grid_fmt).CreateCopy(OutGTiff, OutRaster)
            print('    Created {0} raster: {1}'.format(Variable, OutGTiff))
            target_ds = None

    OutRaster = None
    del Variable, inGeogrid, out_dir, out_Grid_fmt, buildGeogridRaster, OutRaster, OutGTiff
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))
