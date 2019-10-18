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
import shutil

# Import Additional Modules
import osr
import gdal
from gdalnumeric import *
import netCDF4

# Import function library into namespace. Must exist in same directory as this script.
import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox

### Add Proj directory to path
import sys
conda_env_path = os.path.join(os.path.dirname(sys.executable))
internal_datadir = os.path.join(conda_env_path, "Library", "share", "proj")
os.environ["PROJ_LIB"] = internal_datadir

##import pyproj
##print('Script initiated at {0}'.format(time.ctime()))

# Global Variables

# Input and output files and directories
FullDom = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\WRF_Hydro_routing_grids_using_example_case_GEOGRID2\Fulldom_hires.nc'

#out_dir = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs'
out_dir = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\TestGW'

out_zip = os.path.join(out_dir, 'wrf_hydro_routing_grids.zip')
inDEM = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\iowadem.tif'
in_csv = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\Iowa_Gauges_v8.csv'
in_csv = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\WRF_Hydro_routing_grids_using_example_case_GEOGRID2\croton_frxst_pts_csv.csv'
in_lakes = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\WRF_Hydro_routing_grids_using_example_case_GEOGRID2\lakes.shp'

# Variables derived from function script
out_Grid_fmt = wrfh.RasterDriver = 'GTiff'

# Script options
runGEOGRID_STANDALONE = True                                                    # Switch for testing the GEOGRID STANDALONE Pre-processing workflow

# Script parameters
routing = False                                                                 # Build reach-based routing inputs
Lake_routing = False                                                            # Allow gridded lake routing
regridFactor = 4                                                                # Regridding factor
ovroughrtfac_val = 1.0
retdeprtfac_val = 1.0
basin_mask = False
threshold = 200
maskRL = False                                                                  # Allow masking of channels in RouteLink file. May cause WRF-Hydro to crash if True

#outNCType = 'NETCDF3_64BIT'                                                     # Set output netCDF format for spatial metdata files. This was the default before 7/31/2018
outNCType = 'NETCDF4_CLASSIC'                                                   # Define the output netCDF version for RouteLink.nc and LAKEPARM.nc

# Processing Notes to insert into netCDF global attributes
processing_notes_SM = '''Created: %s''' %time.ctime()                           # Processing notes for Spatial Metdata files
processing_notesFD = '''Created: %s''' %time.ctime()                            # Processing notes for the FULLDOM (Routing Grid) file

# List of all possible routing-stack files to keep between the working directory and output .zip files
nclist = [wrfh.LDASFile,
            wrfh.FullDom,
            'gw_basns.nc',
            wrfh.GW_ASCII,
            'gw_basns_geogrid.prj',
            wrfh.RT_nc,
            'Route_Link.csv',
            wrfh.LK_nc,
            'streams.shp', 'streams.shx', 'streams.shp.xml', 'streams.sbx', 'streams.sbn', 'streams.prj', 'streams.dbf',
            'lakes.shp', 'lakes.shx', 'lakes.shp.xml', 'lakes.sbx', 'lakes.sbn', 'lakes.prj', 'lakes.dbf',
            wrfh.GW_nc,
            wrfh.GWGRID_nc]

# Groundwater input options
GW_with_Stack = True                                                            # Switch for building default groundwater inputs with any routing stack
defaultGWmethod = 'FullDom LINKID local basins'                                 # Provide the default groundwater basin generation method. Options ['FullDom basn_msk variable', 'FullDom LINKID local basins', 'Polygon Shapefile or Feature Class']

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()

    if runGEOGRID_STANDALONE:

        # Create scratch directory for temporary outputs
        projdir = os.path.join(os.path.dirname(out_zip), 'scratchdir')
        if os.path.exists(projdir):
            shutil.rmtree(projdir)
        os.makedirs(projdir)

        # Interpret the input for reservoir routing
        if not Lake_routing:
            in_lakes = None

        if GW_with_Stack and defaultGWmethod != 'Polygon Shapefile or Feature Class':
            in_GWPolys = None                                                   # The polygon shapefile to use if defaultGWmethod == 'Polygon Shapefile or Feature Class'

        # Create FULLDOM file CRS info
        rootgrp2 = netCDF4.Dataset(FullDom, 'r+')
        proj1, map_pro, DX2, DY2, x00, y00 = wrfh.georeference_geogrid_file(rootgrp2)
        GeoTransform1 = (x00, DX2*float(regridFactor), 0, y00, 0, DY2*float(regridFactor))  # Build an affine transformation (useful info to add to metadata)

        # Build rasters out of Fulldom_hires variables
        fill = wrfh.numpy_to_Raster(rootgrp2.variables['TOPOGRAPHY'][:], proj1, DX2, DY2, x00, y00)
        flac = wrfh.numpy_to_Raster(rootgrp2.variables['FLOWACC'][:], proj1, DX2, DY2, x00, y00)
        fdir = wrfh.numpy_to_Raster(rootgrp2.variables['FLOWDIRECTION'][:], proj1, DX2, DY2, x00, y00)
        channelgrid = wrfh.numpy_to_Raster(rootgrp2.variables['CHANNELGRID'][:], proj1, DX2, DY2, x00, y00)

        ''' Works great up until here '''


        if 'in_csv' in locals():
            rootgrp2 = wrfh.forecast_points(in_csv, rootgrp2, basin_mask, flac, fdir, channelgrid, out_dir)    # Forecast point processing

        # Moved 10/9/2017 by KMS to allow masking routing files (LINKID, Route_Link, etc.) to forecast points if requested
        if routing:
            print('    Reach-based routing files will be created.')
            order = wrfh.numpy_to_Raster(rootgrp2.variables['STREAMORDER'][:], proj1, DX2, DY2, x00, y00)

            rasterExp = "Value = -9999"                                             # Default: all channels will be represented in reach-based routing file
            if in_csv is None:                                                      # Added 10/10/2017 by KMS to include forecast points in reach-based routing file
                frxst_raster = None                                                 # Default is no forecast points for reach-based routing file
            elif maskRL:
                rasterExp = "Value < 0"                                         # Only channels within the masked basin will be in reach-based routing file
            strm2 = SetNull(channelgrid, channelgrid, rasterExp)                    # Alter channelgrid such that -9999 and -1 to be NoData
            linkid_arr = Routing_Table(projdir, proj1, strm2, fdir, fill, order, gages=frxst_raster)
            rootgrp2.variables['LINKID'][:] = linkid_arr
            print('    Process: LINKID written to output netCDF.')
            strm2 = frxst_raster = order = None
            del linkid_arr, strm2, rasterExp, frxst_raster, order
        else:
            print('    Reach-based routing files will not be created.')

        fdir = None
        del fdir

        if in_lakes is not None:
            # Alter Channelgrid for reservoirs
            channelgrid_arr, outRaster_arr = wrfh.add_reservoirs(channelgrid, in_lakes, flac, projdir, fill, cellsize2, proj1)

            # Process: Output LAKEGRID
            rootgrp2.variables['LAKEGRID'][:] = outRaster_arr
            print('    Process: LAKEGRID written to output netCDF.')

            # Process: Output Channelgrid
            rootgrp2.variables['CHANNELGRID'][:] = channelgrid_arr
            print('    Process: CHANNELGRID written to output netCDF.')

        rootgrp2.close()                                                        # Close Fulldom_hires.nc file
        del rootgrp2

        channelgrid = flac = fill = None
        del channelgrid, flac, fill

        # Step 5 - Build groundwater files -- NOT WORKING YET
        if GW_with_Stack:
            print('  Building Groundwater Basin inputs using default method.')
            GWBasns = wrfh.build_GW_Basin_Raster(out_nc2, projdir, defaultGWmethod, channelgrid, fdir, DX2, DY2, x00, y00, proj1, in_Polys=in_GWPolys)
            wrfh.build_GW_buckets(projdir, GWBasns, DX, DY, proj1, map_pro, GeoTransformStr, Grid=True)
            #GWBasns = None
            #del GWBasns
        wrfh.remove_file(channelgrid)                                           # Delete channelgrid from disk
        wrfh.remove_file(fdir)                                                  # Delete fac from disk
        del proj1, map_pro, GeoTransform1, DX2, DY2, x00, y00

        # zip the folder
        zipper = wrfh.zipUpFolder(projdir, out_zip, nclist)
        print('Built output .zip file after {0: 3.2f} seconds:\n    {1}'.format(time.time()-tic, out_zip))

        # Delete all temporary files
        shutil.rmtree(projdir)
    print('Process completed in {0:3.2f} seconds.'.format(time.time()-tic))