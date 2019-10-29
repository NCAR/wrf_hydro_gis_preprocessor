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

'''
This script will build WRF-Hydro groundwater basin inputs using an input WPS Geogrid
file and WRF-Hydro routing grid (Fulldom_hires.nc) file as inputs. Three methods
are currently available for generating groundwater basins. One method, 'Polygon
Shapefile or Feature Class' currently requires an input polygon shapefile defining
the groundwater basins.
'''

# Import Python Core Modules
import os
import time
import shutil

# Import additional modules
import netCDF4
import gdal

# Import function library into namespace. Must exist in same directory as this script.
import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox
import Examine_Outputs_of_GIS_Preprocessor as EO

# --- Global Variables --- #

# --- EDIT BELOW THIS LINE --- #

# Input and output files and directories
inGeogrid = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\geo_em.d01.nc'
inFulldom = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\TestGW\Fulldom_hires.nc'
out_dir = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\TestGW'

# Outputs - permanent
out_zip = os.path.join(out_dir, 'GroundwaterBasins_local.zip')

# Provide the default groundwater basin generation method.
# Options ['FullDom basn_msk variable', 'FullDom LINKID local basins', 'Polygon Shapefile or Feature Class']
#defaultGWmethod = 'FullDom basn_msk variable'
#defaultGWmethod = 'FullDom LINKID local basins'
defaultGWmethod = 'Polygon Shapefile or Feature Class'

# If the user selects 'Polygon Shapefile or Feature Class' as input, specify path here. Otherwise, in_GWPolys = None
#in_GWPolys = None
#in_GWPolys = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\GW_Basisn_Boundary_Unit.shp'
in_GWPolys = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\GW_Basisn_NHDPlus.shp'

# --- EDIT ABOVE THIS LINE --- #

# --- DO NOT EDIT BELOW THIS LINE --- #

# List of routing-stack files to send to output .zip files
nclist = [wrfh.GW_nc, wrfh.GWGRID_nc]

# Save a GeoTiff of the location of the derived basins on the fine and coarse grids
saveBasins_Fine = False
saveBasins_Coarse = False

# --- End Global Variables --- #

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()
    print('Script initiated at {0}'.format(time.ctime()))

    # Create scratch directory for temporary outputs
    projdir = os.path.join(out_dir, 'gw_scratchdir')
    projdir = os.path.abspath(projdir)
    if os.path.exists(projdir):
        shutil.rmtree(projdir)
    os.makedirs(projdir)

    # Setup temporary output files
    fdir = os.path.join(projdir, wrfh.dir_d8)
    channelgrid = os.path.join(projdir, wrfh.streams)
    if saveBasins_Fine:
        basinRaster_Fine = os.path.join(projdir, 'GWBasins_fine.tif')
        nclist.append('GWBasins_fine.tif')
    if saveBasins_Coarse:
        nclist.append(wrfh.basinRaster)

    rootgrp1 = netCDF4.Dataset(inGeogrid, 'r')
    rootgrp2 = netCDF4.Dataset(inFulldom, 'r')
    coarse_grid = wrfh.WRF_Hydro_Grid(rootgrp1)                                 # Instantiate a grid object for the coarse grid
    fine_grid = wrfh.WRF_Hydro_Grid(rootgrp2)                                   # Instantiate a grid object for the fine grid

    # Build inputs required for creating groundwater buckets
    flowdir = fine_grid.numpy_to_Raster(rootgrp2.variables['FLOWDIRECTION'][:])
    strm_arr = rootgrp2.variables['CHANNELGRID'][:]
    strm_arr[strm_arr==0] = 1                                                   # Set active channels to 1
    strm = fine_grid.numpy_to_Raster(strm_arr)
    del strm_arr

    # Save to disk for the Groundwater tools to use
    out_ds1 = gdal.GetDriverByName(wrfh.RasterDriver).CreateCopy(fdir, flowdir)
    out_ds2 = gdal.GetDriverByName(wrfh.RasterDriver).CreateCopy(channelgrid, strm)
    out_ds1 = out_ds2 = flowdir = strm = None

    # Build groundwater files
    print('  Building Groundwater Basin inputs.')
    GWBasns = wrfh.build_GW_Basin_Raster(inFulldom, projdir, defaultGWmethod, channelgrid, fdir, fine_grid, in_Polys=in_GWPolys)
    wrfh.build_GW_buckets(projdir, GWBasns, coarse_grid, Grid=True, saveRaster=saveBasins_Coarse)
    if saveBasins_Fine:
        out_ds3 = gdal.GetDriverByName(wrfh.RasterDriver).CreateCopy(basinRaster_Fine, GWBasns)
        out_ds3 = None
    GWBasns = None
    wrfh.remove_file(fdir)
    wrfh.remove_file(channelgrid)
    del GWBasns, coarse_grid, fine_grid

    # zip the folder
    tic1 = time.time()
    zipper = wrfh.zipUpFolder(projdir, out_zip, nclist)
    print('Built output .zip file: {0}'.format(out_zip))

    # Delete all temporary files
    shutil.rmtree(projdir)
    rootgrp1.close()
    rootgrp2.close()
    del rootgrp1, rootgrp2
    print('Process completed in {0:3.2f} seconds.'.format(time.time()-tic))