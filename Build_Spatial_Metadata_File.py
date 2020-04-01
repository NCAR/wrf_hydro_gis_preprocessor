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
import copy
from distutils.version import LooseVersion

# Import additional modules
import netCDF4
import gdal
import osr
from gdalnumeric import *
from osgeo import gdal_array

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import (WRF_Hydro_Grid, projdict, flip_grid,
    numpy_to_Raster, wgs84_proj4, ReprojectCoords, outNCType, create_CF_NetCDF)

print('Script initiated at {0}'.format(time.ctime()))

# Global Variables

# Input and output files and directories
inGeogrid = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\geo_em.d01.nc'
out_nc = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\Fulldom_RTOUT_Spatial_Metadata.nc'

# Globals
#format_out = "LDASOUT"                                                          # ["LDASOUT", "RTOUT"]
format_out = "RTOUT"                                                          # ["LDASOUT", "RTOUT"]
inFactor = 4.0                                                                    # Regridding factor
latlon_vars = True                                                             # Include LATITUDE and LONGITUDE 2D variables?

# Processing Notes to insert into netCDF global attributes
processing_notes_SM = '''Created: %s''' %time.ctime()                           # Processing notes for Spatial Metdata files

# Script options

# Methods test switches
coordMethod1 = True                                                             # Interpolate GEOGRID latitude and longitude coordinate arrays
coordMethod2 = False                                                            # Transform coordinate pairs at each grid cell from projected to geocentric

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()
    print('Begining processing on {0}'.format(time.ctime()))
    projdir = os.path.dirname(out_nc)

    if format_out == "LDASOUT":
        regridFactor = 1.0
    elif format_out == "RTOUT":
        regridFactor = inFactor

    # Print informational messages
    print('Input Dataset: {0}'.format(inGeogrid))
    print('Output Grid Resolution: {0}'.format(format_out))
    print('Output Regridding Factor: {0}'.format(regridFactor))
    print('Directory to be used for outputs: {0}'.format(projdir))
    print('Output netCDF File: {0}'.format(out_nc))

    # Georeference geogrid file
    rootgrp = netCDF4.Dataset(inGeogrid, 'r')                                   # Establish an object for reading the input NetCDF file
    if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
        rootgrp.set_auto_mask(False)                                            # Change masked arrays to old default (numpy arrays always returned)
    coarse_grid = WRF_Hydro_Grid(rootgrp)                                  # Instantiate a grid object
    print('    Map Projection of GEOGRID: {0}'.format(projdict[coarse_grid.map_pro]))
    print('    PROJ4: {0}'.format(coarse_grid.proj4))
    print('    Input GeoTransform: {0}'.format(coarse_grid.GeoTransform()))                   # Print affine transformation to screen.

    # Build the regridded domain
    fine_grid = copy.copy(coarse_grid)                                          # Copy the grid object for modification
    fine_grid.regrid(regridFactor)                                              # Regrid to the coarse grid
    print('    Output GeoTransform: {0}'.format(fine_grid.GeoTransform()))      # Print affine transformation to screen.
    print('    New Resolution: {0} {1}'.format(fine_grid.DX, -fine_grid.DY))

    # Build latitude and longitude arrays for GEOGRID_LDASOUT spatial metadata file
    latArr = flip_grid(rootgrp.variables['XLAT_M'][0])                 # Extract array of GEOGRID latitude values
    lonArr = flip_grid(rootgrp.variables['XLONG_M'][0])                # Extract array of GEOGRID longitude values
    rootgrp.close()
    del rootgrp

    if latlon_vars:
        # Build latitude and longitude arrays for Fulldom_hires netCDF file

        if coordMethod1:
            print('  Deriving geocentric coordinates on routing grid from bilinear interpolation of geogrid coordinates.')

            # Method 1: Use GEOGRID latitude and longitude fields and resample to routing grid
            latRaster = numpy_to_Raster(latArr, coarse_grid.proj, coarse_grid.DX, coarse_grid.DY, coarse_grid.x00, coarse_grid.y00)      # Build raster out of GEOGRID latitude array
            lonRaster = numpy_to_Raster(lonArr, coarse_grid.proj, coarse_grid.DX, coarse_grid.DY, coarse_grid.x00, coarse_grid.y00)      # Build raster out of GEOGRID latitude array

            if format_out == "RTOUT":
                latRaster = fine_grid.project_to_model_grid(latRaster)          # Regrid from GEOGRID resolution to routing grid resolution
                lonRaster = fine_grid.project_to_model_grid(lonRaster)          # Regrid from GEOGRID resolution to routing grid resolution

            latArr = BandReadAsArray(latRaster.GetRasterBand(1))                  # Read into numpy array
            lonArr = BandReadAsArray(lonRaster.GetRasterBand(1))                  # Read into numpy array
            latRaster = lonRaster = None                                          # Destroy raster objects
            del latRaster, lonRaster

        elif coordMethod2:
            print('  Deriving geocentric coordinates on routing grid from direct transformation geogrid coordinates.')

            # Method 2: Transform each point from projected coordinates to geocentric coordinates
            wgs84_proj = osr.SpatialReference()                                 # Build empty spatial reference object
            wgs84_proj.ImportFromProj4(wgs84_proj4)                        # Imprort from proj4 to avoid EPSG errors (4326)

            if format_out == "RTOUT":
                pass

            xmap, ymap = coarse_grid.getxy()                                    # Get x and y coordinates as numpy array
            latArr, lonArr = ReprojectCoords(xmap, ymap, fine_grid.proj, wgs84_proj)   # Transform coordinate arrays
            del xmap, ymap, wgs84_proj
    else:
        latArr = lonArr = None

    # Create the netCDF file with spatial metadata
    rootgrp2 = netCDF4.Dataset(out_nc, 'w', format=outNCType)
    rootgrp2, grid_mapping = create_CF_NetCDF(fine_grid, rootgrp2, projdir,
            notes=processing_notes_SM, addLatLon=latlon_vars, latArr=latArr, lonArr=lonArr)
    rootgrp2.close()

    del rootgrp2, latArr, lonArr
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))