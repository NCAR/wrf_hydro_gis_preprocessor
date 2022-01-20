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

descText = "This tool takes the output zip file from the ProcessGeogrid script and creates a raster " \
           "from each output NetCDF file. The Input should be a .zip file that was created using the" \
           "WRF Hydro pre-processing tools. The tool will create the folder which will contain the" \
           "results (out_folder), if that folder does not already exist. "

# Import Python Core Modules
import os
import sys
import time
import shutil
from distutils.version import LooseVersion

# Import additional modules
import netCDF4
from argparse import ArgumentParser
from pathlib import Path
import osgeo

try:
    if LooseVersion(osgeo.__version__) > LooseVersion('3.0.0'):
        from osgeo import gdal
        from osgeo import osr
        from osgeo.gdal_array import *
    else:
        import gdal
        import osr
        from gdal_array import *
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')
#from gdalnumeric import *                                                      # Assists in using BandWriteArray, BandReadAsArray, and CopyDatasetInfo

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import (LK_nc, RT_nc, GW_nc, LDASFile, crsVar,
    numpy_to_Raster, ZipCompat)

# Global Variables

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

                # Old method which will crash if the CRS variable is not present
                ##                if 'esri_pe_string' in rootgrp.variables[crsVar].__dict__:
                ##                    PE_string = rootgrp.variables[crsVar].esri_pe_string
                ##                elif 'spatial_ref' in rootgrp.variables[crsVar].__dict__:
                ##                    PE_string = rootgrp.variables[crsVar].spatial_ref
                ##                GT = rootgrp.variables[crsVar].GeoTransform.split(" ")[0:6]

                # Added 4/14/2021 to allow for the absence of a coordinate system variable.
                if crsVar in rootgrp.variables:
                    crsNCVar = rootgrp.variables[crsVar]
                    if 'esri_pe_string' in crsNCVar.__dict__:
                        #PE_string = crsNCVar.esri_pe_string
                        PE_string = crsNCVar.esri_pe_string.replace("'", '"')
                    elif 'spatial_ref' in crsNCVar.__dict__:
                        PE_string = crsNCVar.spatial_ref
                    if 'GeoTransform' in crsNCVar.__dict__:
                        GT = crsNCVar.GeoTransform.split(" ")[0:6]
                    else:
                        print('  No GeoTransform attribute found. Setting to default.')
                        GT = [0, 1, 0, 0, 0, -1]
                else:
                    # Create dummy variables to allow the script to continue
                    PE_string = ''
                    GT = [0, 1, 0, 0, 0, -1]

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

    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_zip",
                        default='',
                        required=True,
                        help="Path to WRF Hydro routing grids zip file.")
    parser.add_argument("-o",
                        dest="out_folder",
                        default='',
                        required=True,
                        help="Path to output folder.")

    # If no arguments are supplied, print help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    if args.in_zip == all_defaults["in_zip"]:
        print('Using input zip location of: {0}'.format(all_defaults["in_zip"]))

    if args.out_folder == all_defaults["out_folder"]:
        print('Using output location of: {0}'.format(all_defaults["out_folder"]))

    # Create output directory for temporary outputs
    if os.path.exists(args.out_folder):
        print('Requested output directory already exists. \nPlease specify a non-existant directory as output.')
        raise SystemExit
    else:
        os.makedirs(args.out_folder)

    # Unzip to a known location (make sure no other nc files live here)
    ZipCompat(args.in_zip).extractall(args.out_folder)
    examine_outputs(args.out_folder, skipfiles=skipfiles)
    print('Extraction of WRF routing grids completed.')
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))
