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

"""
This tool takes the output zip file from the ProcessGeogrid script and creates a raster
from each output NetCDF file. The Input should be a .zip file that was created using the
WRF Hydro pre-processing tools. The tool will create the folder which will contain the
results (out_folder), if that folder does not already exist.


 Ensure that the out_folder exists, but is empty for maximum compatibility.
 """

# Import Python Core Modules
import os
import time
import shutil
from distutils.version import LooseVersion

# Import additional modules
import netCDF4
import osr
import gdal
from osgeo import gdal_array
from gdalnumeric import *                                                       # Assists in using BandWriteArray, BandReadAsArray, and CopyDatasetInfo

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import (LK_nc, RT_nc, GW_nc, LDASFile, crsVar,
    numpy_to_Raster, ZipCompat)

# Global Variables

# Input and output files and directories
#in_zip = r"C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\wrf_hydro_routing_grids_Iowa2.zip"
#out_folder = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\wrf_hydro_routing_grids_Examine'

#in_zip = r"C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\FOSS_Domain\Alaska_WPS_DEM_r4_t25_gridded_breachLC.zip"
in_zip = r"C:\Users\ksampson\Desktop\Youcan_Feng\Routing_stack_lakes\WRF_Hydro_routing_grids_Neuse_V2_nomask.zip"

#out_folder = r"C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\FOSS_Domain\Alaska_WPS_DEM_r4_t25_gridded_breachLC_Examine"
out_folder = r"C:\Users\ksampson\Desktop\Youcan_Feng\Routing_stack_lakes_Examine"

# Script Options
RasterDriver = 'GTiff'                                                          # Driver for output raster format
suffix = '.tif'                                                                 # File extension to use for output rasters
skipfiles = []                                                                  # Files that should not be converted or written to output directory

# --- Functions --- #
def examine_outputs(out_folder, dellist=[], skipfiles=[]):
    '''
    Provide a directory, ideally the unzipped directory of a routing stack produced
    by the WRF-Hydro GIS Pre-processor. Files will be examined and derivatives
    made from 2D netCDF files. Some files will be delted from the input directory.
    '''
    tic1 = time.time()
    dellist = []                            # Initialize list to store files to delete from output directory

    # Iterate through unzipped files and copy to output directory as necessary
    for dirpath, dirnames, filenames in os.walk(out_folder):
        for file in filenames:
            infile = os.path.join(dirpath, file)

            # Copy skipped files over to new directory
            if file in skipfiles:
                dellist.append(infile)
                print('    File NOT Copied: {0}'.format(file))
                del file, infile
                continue

            # Trap to eliminate Parameter tables in NC format from this extraction
            if file in [LK_nc, RT_nc, GW_nc, LDASFile]:
                print('    File Copied: {0}'.format(file))
                del file, infile
                continue

            if file.endswith('.nc'):

                # Establish an object for reading the input NetCDF file
                rootgrp = netCDF4.Dataset(infile, 'r')
                if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
                    rootgrp.set_auto_mask(False)                                # Change masked arrays to old default (numpy arrays always returned)

                # Using netCDF4 library - still need point2, DX, DY
                if 'esri_pe_string' in rootgrp.variables[crsVar].__dict__:
                    PE_string = rootgrp.variables[crsVar].esri_pe_string
                elif 'spatial_ref' in rootgrp.variables[crsVar].__dict__:
                    PE_string = rootgrp.variables[crsVar].spatial_ref
                GT = rootgrp.variables[crsVar].GeoTransform.split(" ")[0:6]
                GT = tuple(float(item) for item in GT)
                print('  GeoTransform: {0}'.format(GT))
                print('  DX: {0}'.format(GT[1]))
                print('  DY: {0}'.format(-GT[5]))

                proj = osr.SpatialReference()                                   # Initiate OSR spatial reference object
                proj.ImportFromWkt(PE_string)
                print('  PROJ.4 string: {0}'.format(proj.ExportToProj4()))
                for variablename, ncvar in rootgrp.variables.items():
                    if ncvar.dimensions==('y', 'x'):
                        OutRaster = numpy_to_Raster(ncvar[:].copy(), proj, GT[1], GT[5], GT[0], GT[3])

                        # Save to disk
                        OutGTiff = os.path.join(out_folder, variablename+suffix)# Output raster

                        ##                        gdaltype = gdal_array.NumericTypeCodeToGDALTypeCode(ncvar.dtype)
                        ##                        rows = ncvar.shape[ncvar.dimensions.index('y')]
                        ##                        cols = ncvar.shape[ncvar.dimensions.index('x')]
                        ##                        target_ds = gdal.GetDriverByName(RasterDriver).Create(OutGTiff, rows, cols, 1, gdaltype)
                        ##                        CopyDatasetInfo(OutRaster, target_ds)
                        ##                        target_ds.FlushCache()                                  #saves to disk!!
                        ##                        target_ds = OutRaster = None

                        try:
                            target_ds = gdal.GetDriverByName(RasterDriver).CreateCopy(OutGTiff, OutRaster)
                            target_ds = OutRaster = None
                            del target_ds
                        except:
                            pass
                        del OutRaster

                        print('    File Created: {0}'.format(OutGTiff))
                        del OutGTiff, variablename
                rootgrp.close()
                dellist.append(infile)
                del file, infile, rootgrp, ncvar
                continue

            if file.split('.')[-1] in ['shp', 'shx', 'xml', 'sbx', 'sbn', 'prj', 'dbf']:
                print('    File Copied: {0}'.format(str(file)))
                del file, infile
                continue

            # These file formats are legacy output files, but allow the tool to work with older routing stacks
            if file.endswith('.csv') or file.endswith('.TBL') or file.endswith('.txt') or file.endswith('.prj'):
                print('    File Copied: {0}'.format(file))
                del file, infile
                continue

            dellist.append(infile)
            del file, infile
            continue
    del dirpath, dirnames, filenames

    # Remove each file from the temporary extraction directory
    for infile in dellist:
        os.remove(infile)
    return

# --- End Functions --- #

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()
    print('Script initiated at {0}'.format(time.ctime()))

    # Create output directory for temporary outputs
    if os.path.exists(out_folder):
        print('Requested output directory already exists. \nPlease specify a non-existant directory as output.')
        raise SystemExit
    else:
        os.makedirs(out_folder)

    # Unzip to a known location (make sure no other nc files live here)
    ZipCompat(in_zip).extractall(out_folder)
    examine_outputs(out_folder, skipfiles=skipfiles)
    print('Extraction of WRF routing grids completed.')
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))