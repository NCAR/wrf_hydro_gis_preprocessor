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

descText = "This tool takes an input raster (most likely produced using the ExportGrid tool)" \
           " and uses that grid to produce latitude and longitude ESRI GRID rasters."

# Import Modules

# Import Python Core Modules
import os
import sys
import time

# Import additional modules
import gdal
from argparse import ArgumentParser

try:
    if sys.version_info >= (3, 0):
        from osgeo import osr
    else:
        import osr
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

# Import function library into namespace. Must exist in same directory as this script.
# Import wrfhydro_functions as wrfh, Function script packaged with this toolbox
from wrfhydro_functions import (get_projection_from_raster, wgs84_proj4, getxy,
    ReprojectCoords, numpy_to_Raster)

# Global Variables

# Script options
RasterDriver = 'GTiff'
wgs84_proj4 = '+proj=longlat +datum=WGS84 +no_defs'

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()
    print('Script initiated at {0}'.format(time.ctime()))

    # Setup the input arguments
    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_raster",
                        default='',
                        required=True,
                        help="Path to input raster.")
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

    if args.in_raster == all_defaults["in_raster"]:
        print('Using input raster location of: {0}'.format(all_defaults["in_raster"]))

    if args.out_dir == all_defaults["out_dir"]:
        print('Using output location of: {0}'.format(all_defaults["out_dir"]))

    outLat = os.path.join(args.out_dir, 'LATITUDE.tif')
    outLon = os.path.join(args.out_dir, 'LONGITUDE.tif')

    # Open (read-only) input raster
    in_raster = gdal.Open(args.in_raster, 0)                                          # Open with read-only mode

    # Gather information from input raster projection
    proj = get_projection_from_raster(in_raster)
    x00, DX, xskew, y00, yskew, DY = in_raster.GetGeoTransform()
    del xskew, yskew

    print('  Deriving geocentric coordinates on routing grid from direct transformation geogrid coordinates.')
    # Note that this might not be the same method used in the main pre-processing script
    # because a geogrid file is not required here as input.
    wgs84_proj = osr.SpatialReference()                                 # Build empty spatial reference object
    wgs84_proj.ImportFromProj4(wgs84_proj4)                        # Import from proj4 to avoid EPSG errors (4326)
    xmap, ymap = getxy(in_raster)                                          # Get x and y coordinates as numpy array
    in_raster = None
    lonArr2, latArr2 = ReprojectCoords(xmap, ymap, proj, wgs84_proj)  # Transform coordinate arrays
    del wgs84_proj, in_raster, xmap, ymap

    # Convert back to rasters
    xmap = numpy_to_Raster(lonArr2, proj, DX, DY, x00, y00)
    ymap = numpy_to_Raster(latArr2, proj, DX, DY, x00, y00)
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
    print('Process completed in {0:3.2f} seconds.'.format(time.time()-tic))
