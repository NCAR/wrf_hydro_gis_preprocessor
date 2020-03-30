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
import sys
import time
import shutil
import copy
from distutils.version import LooseVersion

# Import Additional Modules
import osr
import gdal
from gdalnumeric import *
import netCDF4

# Import function library into namespace. Must exist in same directory as this script.
import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox
##from wrfhydro_functions import (RasterDriver, LDASFile, FullDom, GW_ASCII, RT_nc,
##    LK_nc, GW_nc, GWGRID_nc, GWGRID_nc, minDepthCSV, TeeNoFile, WRF_Hydro_Grid,
##    flip_grid, create_CF_NetCDF, wgs84_proj4, ReprojectCoords, FullDom, WB_functions,
##    remove_file, forecast_points, Routing_Table, add_reservoirs, build_GW_Basin_Raster,
##    build_GW_buckets, zipUpFolder)

# Add Proj directory to path
conda_env_path = os.path.join(os.path.dirname(sys.executable))
internal_datadir = os.path.join(conda_env_path, "Library", "share", "proj")
os.environ["PROJ_LIB"] = internal_datadir

##import pyproj
##print('Script initiated at {0}'.format(time.ctime()))

# --- Global Variables --- #

# Input and output files and directories

# CONUS
#inGeogrid = r'C:\Users\ksampson\Desktop\GLERL_Visit_201910\Soils_BNU\geo_em.d01_varsubset.nc'
#inDEM = r"C:\Data\gisdata\Elevation\HydroSheds_mosaic_NWM_30m.vrt"
##inDEM = r'C:\Data\Projects\Gochis\NWM_v1_2\Terrain_Update\Data\Combined_NED_HS_Elevation_2.tif'

# Iowa
inGeogrid = r"C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\geo_em.d02.20200327_snow.trim.nc"
#inDEM = r"C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\FOSS_Domain\HGT_M.tif"
inDEM = r"C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\geo_em.d03.20200327_snow.trim.HGTM_SetNull_Elev0_Mask0.tif"
in_csr = r''
in_lakes = r''
#in_csv = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\Iowa_Gauges_v8.csv'
#in_lakes = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\NWM_v_2_1_Reservoirs_Preliminary_20190510.shp'

in_GWPolys = None                                                               # The polygon shapefile to use if defaultGWmethod == 'Polygon Shapefile or Feature Class'

# Outputs - permanent
out_dir = r'C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\FOSS_Domain\Alaska_WPS_DEM_r4_t25_gridded_WBT1_2_0'
out_zip = os.path.join(out_dir, 'Alaska_WPS_DEM_r4_t25_gridded_noflatfix.zip')

# Script parameters
routing = False                                                                 # Build reach-based routing inputs
Lake_routing = False                                                           # Allow gridded lake routing
regridFactor = 4                                                                # Regridding factor
ovroughrtfac_val = 1.0
retdeprtfac_val = 1.0
basin_mask = True
threshold = 25
lksatfac_val = 1000.0

# Provide the default groundwater basin generation method.
# Options ['FullDom basn_msk variable', 'FullDom LINKID local basins', 'Polygon Shapefile or Feature Class']
defaultGWmethod = 'FullDom LINKID local basins'
GW_with_Stack = True                                                            # Switch for building default groundwater inputs with any routing stack

# --- DO NOT EDIT BELOW THIS LINE --- #

# Methods test switches
coordMethod1 = True                                                             # Interpolate GEOGRID latitude and longitude coordinate arrays
coordMethod2 = False                                                            # Transform coordinate pairs at each grid cell from projected to geocentric

# Default temporary output file names
mosprj_name = 'mosaicprj.tif'                                                   # Default regridded input DEM if saved to disk

# Variables derived from function script
out_Grid_fmt = wrfh.RasterDriver

# Script options
runGEOGRID_STANDALONE = True                                                    # Switch for testing the GEOGRID STANDALONE Pre-processing workflow

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
            wrfh.GWGRID_nc,
            wrfh.minDepthCSV]

'''Pre-defining the variables and populating variable attributes is
a much faster strategry than creating and populating each variable
sequentially, especially for netCDF3 versions. Also, unsigned integer
types are only allowed in NETCDF4.'''
# List of variables to create [<varname>, <vardtype>, <long_name>]
varList2D = [['CHANNELGRID', 'i4', ''],
            ['FLOWDIRECTION', 'i2', ''],
            ['FLOWACC', 'i4', ''],
            ['TOPOGRAPHY', 'f4', ''],
            ['RETDEPRTFAC', 'f4', ''],
            ['OVROUGHRTFAC', 'f4', ''],
            ['STREAMORDER', 'i1', ''],
            ['frxst_pts', 'i4', ''],
            ['basn_msk', 'i4', ''],
            ['LAKEGRID', 'i4', ''],
            ['landuse', 'f4', ''],
            ['LKSATFAC', 'f4', '']]

# --- End Global Variables --- #

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()

    # Configure logging
    logfile = out_zip.replace('.zip', '.log')
    tee = wrfh.TeeNoFile(logfile, 'w')

    if runGEOGRID_STANDALONE:

        # Create scratch directory for temporary outputs
        projdir = os.path.join(os.path.dirname(out_zip), 'scratchdir')
        projdir = os.path.abspath(projdir)
        if os.path.exists(projdir):
            shutil.rmtree(projdir)
        os.makedirs(projdir)

        # Interpret the input for reservoir routing
        if not Lake_routing:
            in_lakes = None
        gridded = not routing                                                   # Flag for gridded routing

        if 'in_csv' in locals():
            AddGages = True
        else:
            AddGages = False

        # Add variables depending on the input options
        if routing:
            varList2D.append(['LINKID', 'i4', ''])

        # Step 1 - Georeference geogrid file
        rootgrp = netCDF4.Dataset(inGeogrid, 'r')                               # Establish an object for reading the input NetCDF files
        if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
            rootgrp.set_auto_mask(False)                                        # Change masked arrays to old default (numpy arrays always returned)
        coarse_grid = wrfh.WRF_Hydro_Grid(rootgrp)                              # Instantiate a grid object
        fine_grid = copy.copy(coarse_grid)                                      # Copy the grid object for modification
        fine_grid.regrid(regridFactor)                                          # Regrid to the coarse grid
        print('    Proj4: {0}'.format(coarse_grid.proj4))                       # Print Proj.4 string to screen
        print('    Created projection definition from input NetCDF GEOGRID file.')

        # Build output raster from numpy array of the GEOGRID variable requested. This will be used as a template later on
        LU_INDEX = coarse_grid.numpy_to_Raster(wrfh.flip_grid(rootgrp.variables['LU_INDEX'][:]))
        print('    GeoTransform: {0}'.format(coarse_grid.GeoTransformStr()))                  # Print affine transformation to screen.

        # Create spatial metadata file for GEOGRID/LDASOUT grids
        out_nc1 = os.path.join(projdir, wrfh.LDASFile)
        rootgrp1 = netCDF4.Dataset(out_nc1, 'w', format=outNCType)              # wrf_hydro_functions.outNCType)
        rootgrp1, grid_mapping = wrfh.create_CF_NetCDF(coarse_grid, rootgrp1, projdir,
                notes=processing_notes_SM) # addLatLon=True, latArr=latArr, lonArr=lonArr)
        rootgrp1.close()
        del rootgrp1

        # Step 3 - Create high resolution topography layers
        in_DEM = gdal.Open(inDEM, 0)                                     # Open with read-only mode
        outDEM = os.path.join(projdir, mosprj_name)
        mosprj = fine_grid.project_to_model_grid(in_DEM, saveRaster=True, OutGTiff=outDEM, resampling=gdal.GRA_Bilinear)
        in_DEM = mosprj = None

        # Build latitude and longitude arrays for Fulldom_hires netCDF file
        if coordMethod1:
            print('  Deriving geocentric coordinates on routing grid from bilinear interpolation of geogrid coordinates.')
            # Build latitude and longitude arrays for GEOGRID_LDASOUT spatial metadata file
            latArr = wrfh.flip_grid(rootgrp.variables['XLAT_M'][:])             # Extract array of GEOGRID latitude values
            lonArr = wrfh.flip_grid(rootgrp.variables['XLONG_M'][:])            # Extract array of GEOGRID longitude values

            # Method 1: Use GEOGRID latitude and longitude fields and resample to routing grid
            latRaster1 = coarse_grid.numpy_to_Raster(latArr)                    # Build raster out of GEOGRID latitude array
            lonRaster1 = coarse_grid.numpy_to_Raster(lonArr)                    # Build raster out of GEOGRID longitude array

            latRaster2 = fine_grid.project_to_model_grid(latRaster1)            # Regrid from GEOGRID resolution to routing grid resolution
            lonRaster2 = fine_grid.project_to_model_grid(lonRaster1)            # Regrid from GEOGRID resolution to routing grid resolution
            latRaster1 = lonRaster1 = None                                      # Destroy rater objects
            latArr2 = BandReadAsArray(latRaster2.GetRasterBand(1))              # Read into numpy array
            lonArr2 = BandReadAsArray(lonRaster2.GetRasterBand(1))              # Read into numpy array
            latRaster2 = lonRaster2 = None                                      # Destroy raster objects
            del latArr, lonArr, latRaster1, lonRaster1, latRaster2, lonRaster2

        elif coordMethod2:
            print('  Deriving geocentric coordinates on routing grid from direct transformation geogrid coordinates.')
            # Method 2: Transform each point from projected coordinates to geocentric coordinates
            wgs84_proj = osr.SpatialReference()                                 # Build empty spatial reference object
            wgs84_proj.ImportFromProj4(wrfh.wgs84_proj4)                        # Imprort from proj4 to avoid EPSG errors (4326)
            xmap, ymap = fine_grid.getxy()                                      # Get x and y coordinates as numpy array
            latArr2, lonArr2 = wrfh.ReprojectCoords(xmap, ymap, coarse_grid.proj, wgs84_proj)  # Transform coordinate arrays
            del xmap, ymap, wgs84_proj

        # Create FULLDOM file
        out_nc2 = os.path.join(projdir, wrfh.FullDom)
        rootgrp2 = netCDF4.Dataset(out_nc2, 'w', format=outNCType)              # wrf_hydro_functions.outNCType)
        rootgrp2, grid_mapping = wrfh.create_CF_NetCDF(fine_grid, rootgrp2, projdir,
                notes=processing_notesFD, addVars=varList2D, addLatLon=True,
                latArr=latArr2, lonArr=lonArr2)
        del latArr2, lonArr2

        # Add some global attribute metadata to the Fulldom file, including relevant WPS attributes for defining the model coordinate system
        rootgrp2.geogrid_used = inGeogrid                                       # Paste path of geogrid file to the Fulldom global attributes
        rootgrp2.DX = fine_grid.DX                                              # Add X resolution as a global attribute
        rootgrp2.DY = -fine_grid.DY                                             # Add Y resolution as a global attribute
        globalAtts = rootgrp.__dict__                                           # Read all global attributes into a dictionary
        for item in ['MAP_PROJ', 'corner_lats', 'corner_lons', 'TRUELAT1', 'TRUELAT2', 'STAND_LON', 'POLE_LAT', 'POLE_LON', 'MOAD_CEN_LAT', 'CEN_LAT']:
            if item in globalAtts:
                rootgrp2.setncattr(item, globalAtts[item])
        rootgrp.close()                                                         # Close input GEOGRID file
        del item, globalAtts, rootgrp

        # Process: Resample LU_INDEX grid to a higher resolution
        LU_INDEX2 = fine_grid.project_to_model_grid(LU_INDEX, fine_grid.DX, fine_grid.DY, resampling=gdal.GRA_NearestNeighbour)
        rootgrp2.variables['landuse'][:] = BandReadAsArray(LU_INDEX2.GetRasterBand(1))          # Read into numpy array
        LU_INDEX = None                                                         # Destroy raster object
        print('    Process: landuse written to output netCDF.')
        del LU_INDEX, LU_INDEX2

        ##        # Step X(a) - Test to match LANDMASK - Only used for areas surrounded by water (LANDMASK=0)
        ##        mosprj2, loglines = wrfh.adjust_to_landmask(mosprj, LANDMASK, coarse_grid.proj, projdir, 'm')
        ##        outtable.writelines("\n".join(loglines) + "\n")
        ##        del LANDMASK

        # Step 4 - Hyrdo processing functions -- Whitebox
        rootgrp2, fdir, fac, channelgrid, fill, order = wrfh.WB_functions(rootgrp2, outDEM,
                projdir, threshold, ovroughrtfac_val, retdeprtfac_val, lksatfac_val)
        wrfh.remove_file(outDEM)                                                # Delete output DEM from disk

        # If the user provides forecast points as a CSV file, alter outputs accordingly
        if AddGages:
            if os.path.exists(in_csv):
                rootgrp2 = wrfh.forecast_points(in_csv, rootgrp2, basin_mask, projdir,
                            fine_grid.DX, fine_grid.WKT, fdir, fac, channelgrid)    # Forecast point processing

        # Moved 10/9/2017 by KMS to allow masking routing files (LINKID, Route_Link, etc.) to forecast points if requested
        if routing:
            print('  Reach-based routing files will be created.')
            rootgrp2 = wrfh.Routing_Table(projdir, rootgrp2, fine_grid, fdir, channelgrid, fill, order, gages=AddGages)
        else:
            print('  Reach-based routing files will not be created.')
        wrfh.remove_file(fill)                                                  # Delete fill from disk
        wrfh.remove_file(order)                                                 # Delete order from disk

        if Lake_routing:
            # Alter Channelgrid for reservoirs and build reservoir inputs
            pass
            rootgrp2 = wrfh.add_reservoirs(rootgrp2, projdir, fac, in_lakes, fine_grid, Gridded=gridded)
        rootgrp2.close()                                                        # Close Fulldom_hires.nc file
        del rootgrp2

        # Build groundwater files
        if GW_with_Stack:
            GWBasns = wrfh.build_GW_Basin_Raster(out_nc2, projdir, defaultGWmethod, channelgrid, fdir, fine_grid, in_Polys=in_GWPolys)
            wrfh.build_GW_buckets(projdir, GWBasns, coarse_grid, Grid=True)
            GWBasns = None
        wrfh.remove_file(fdir)                                                  # Delete fdir from disk
        wrfh.remove_file(fac)                                                   # Delete fac from disk
        wrfh.remove_file(channelgrid)                                           # Delete channelgrid from disk

        # Copmress (zip) the output directory
        tic1 = time.time()
        zipper = wrfh.zipUpFolder(projdir, out_zip, nclist)
        print('Built output .zip file in {0: 3.2f} seconds.'.format(time.time()-tic1))  # Diagnotsitc print statement

        # Delete all temporary files
        shutil.rmtree(projdir)

    tee.close()
    del tee                                                                     # Should do nothing
    print('Process completed in {0:3.2f} seconds.'.format(time.time()-tic))