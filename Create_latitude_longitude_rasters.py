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
import gdal
import osr

# Import function library into namespace. Must exist in same directory as this script.
import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox

print('Script initiated at {0}'.format(time.ctime()))

# Global Variables

# Input and output files and directories
inRaster = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\TOPOGRAPHY_HUGEBACK.tif'
out_dir = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs'

# Script options
RasterDriver = 'GTiff'
wgs84_proj4 = '+proj=longlat +datum=WGS84 +no_defs'

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()

    outLat = os.path.join(out_dir, 'LATITUDE.tif')
    outLon = os.path.join(out_dir, 'LONGITUDE.tif')

    # Open (read-only) input raster
    in_raster = gdal.Open(inRaster, 0)                                          # Open with read-only mode

    # Gather information from input raster projection
    proj = wrfh.get_projection_from_raster(in_raster)
    x00, DX, xskew, y00, yskew, DY  = in_raster.GetGeoTransform()
    del xskew, yskew

    print('  Deriving geocentric coordinates on routing grid from direct transformation geogrid coordinates.')
    # Note that this might not be the same method used in the main pre-processing script
    # because a geogrid file is not required here as input.
    wgs84_proj = osr.SpatialReference()                                 # Build empty spatial reference object
    wgs84_proj.ImportFromProj4(wrfh.wgs84_proj4)                        # Imprort from proj4 to avoid EPSG errors (4326)
    xmap, ymap = wrfh.getxy(in_raster)                                          # Get x and y coordinates as numpy array
    in_raster = None
    lonArr2, latArr2 = wrfh.ReprojectCoords(xmap, ymap, proj, wgs84_proj)  # Transform coordinate arrays
    del wgs84_proj, in_raster, xmap, ymap

    # Convert back to rasters
    xmap = wrfh.numpy_to_Raster(lonArr2, proj, DX, DY, x00, y00)
    ymap = wrfh.numpy_to_Raster(latArr2, proj, DX, DY, x00, y00)
    del DX, DY, x00, y00, proj, latArr2, lonArr2

    # Section below causing a RuntimeError
    out_drv = gdal.GetDriverByName(RasterDriver)
    if ymap is not None:
        try:
            target_ds = out_drv.CreateCopy(outLat, ymap)
            target_ds = None
        except:
            pass
    if xmap is not None:
        try:
            target_ds = out_drv.CreateCopy(outLon, xmap)
            target_ds = None
        except:
            pass
    del out_drv, xmap, ymap

    ##    # Save output rasters
    ##    for OutGTiff, InRaster in zip([outLat, outLon], [ymap, xmap]):
    ##        if InRaster is not None:
    ##            target_ds = gdal.GetDriverByName(RasterDriver).CreateCopy(OutGTiff, InRaster)
    ##            target_ds = None
    ##    del OutGTiff, InRaster, xmap, ymap
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))
