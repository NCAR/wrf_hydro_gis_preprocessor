# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Copyright UCAR (c) 2021
# University Corporation for Atmospheric Research(UCAR)
# National Center for Atmospheric Research(NCAR)
# Research Applications Laboratory(RAL)
# P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# 06/05/2021
#
# Name:        create_wrfinput.py
# Purpose:
# Author:      Kevin Sampson
# Created:     06/05/2021
# Licence:
#
# 10/20/2018:
#    The purpose of this script is to ensure that the dominant soil types more
#    closely match the landmask, such that no water soil types are present where
#    landcover indicates a land type.
#
# Based on work by A. Dugger (NCAR)
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

# Import Python Core Modules
import os
import sys
import getopt
import time
import shutil
from argparse import ArgumentParser

# Import additional Python Modules
import xarray as xr
import numpy

try:
    if sys.version_info >= (3, 0):
        from osgeo import gdal
        from osgeo import gdalconst
        from osgeo import gdalconst
    else:
        import gdal
        import gdalconst
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')
from gdalconst import *

# Import whitebox
from whitebox.WBT.whitebox_tools import WhiteboxTools

# Import functions from the WRF-Hydro GIS Pre-processing scripts
from Build_GeoTiff_From_Geogrid_File import build_geogrid_raster
import wrfhydro_functions as wrfh

# Module configurations
sys.dont_write_bytecode = True

# Screen print in case invalid parameters are given
descText = '''Utility for harmonizing the SCT_DOM and SCB_DOM variables in a WRF/
                WRF-Hydro geogrid file with the LANDMASK variable. This is done
                to ensure that no water soil types are present where the LU_INDEX
                or LANDMASK indicate a land type. The reason a mismatch can happen
                is that Geogrid.exe does not ensure these variables will be consistent.
                Another reason these might not be consistent is if the user substitutes
                their own landcover dataset into LU_INDEX. For the purposes of
                hydrologic simulation, these layers should be consistent when it
                comes to water.'''

'''
To import and run these functions using python, from a custom script or console:

    from Harmonize_Soils_to_LANDMASK import update_geogrid_soils
    update_geogrid_soils(r'./geo_em.d01.nc', r'./geo_em.d01.modifiedSoils.nc')
'''

#######################################################
# Global Variables - update relevant arguments below.
#######################################################

# Soil type to use as a fill value in case conflicts between soil water and land cover water cells:
# If the script encounters a cell that is classified as land in the land use field (LU_INDEX)
# but is classified as a water soil type, it will replace the soil type with the value you
# specify below. Ideally there are not very many of these, so you can simply choose the most
# common soil type in your domain. Alternatively, you can set to a "bad" value (e.g., -8888)
# to see how many of these conflicts there are. If you do this DO NOT RUN THE MODEL WITH THESE
# BAD VALUES. Instead, fix them manually with a neighbor fill or similar fill algorithm.
fillsoiltyp = 3

# Soil water type value from input GEOGRID SCT_DOM and SCB_DOM variables
soil_water_val = 14

# Search distance for filling water cell areas over land with nearest non-water
# neighbor value, in cells. (Cell width is used to determine euclidean distance
search_dist = 3

#######################################################
# Do not update below here.
#######################################################

# Data model for output netCDF data ['NETCDF4', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC']
outNCType = 'NETCDF4'                                                           # Data model for output netCDF data

# Use the GeoTiff raster driver for GDAL, which is the only format accepted by Whitebox Tools
RasterDriver = 'GTiff'

# Controls related to the search distance for gap-filling inland soil-wate classes
distance_calc = True        # Use a search distance threshold for nearest neighbor?
distance = 3                # Distance in cell-widths. Uses GEOGRID DX global attribute

# --- Functions --- #
def is_valid_file(parser, arg):
    # https://stackoverflow.com/questions/11540854/file-as-command-line-argument-for-argparse-error-message-if-argument-is-not-va
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return str(arg)

def update_geogrid_soils(inNC, outFile, distance_calc=True, outNCType=outNCType):
    '''
    This function will populate the arrays in the input xarray DataSet based on
    what is in the input GEOGRID file.

    This may only be needed if you supply and insert your own land-cover layer
    into the LU_INDEX variable in the geogrid file, which will require changes to
    other layers for maximum compatibilty between landcover and soil types.

    '''

    # Create scratch directory for temporary outputs
    projdir = os.path.join(os.path.dirname(inNC), 'scratchdir')
    projdir = os.path.abspath(projdir)
    if os.path.exists(projdir):
        print('Removing existing scratch directory: {0}'.format(projdir))
        shutil.rmtree(projdir)
    os.makedirs(projdir)

    # Whitebox options for running Whitebox in a full workflow
    wbt = WhiteboxTools()
    print('    Using {0}'.format(wbt.version().split('\n')[0]))
    wbt.work_dir = projdir                              # Set working directory
    wbt.verbose = False                                 # Verbose output. [True, False]

    # Open input files for reading. Lazy loading
    print('Input netCDF file: {0}'.format(inNC))
    ncDS = xr.open_dataset(inNC)                    # , decode_cf=False

    # Handle whether or not to modify files in place or create new output files.
    if inNC == outFile:
        print('OVERWRITING INPUT FILE')
        # !!! OVERWRITE exiting input file !!!
        # Load the entire dataset into memory. This allows you to overwrite the existing file
        outFile = inNC
        with xr.open_dataset(inNC) as dataset:
            ncDS = dataset.load()

    # Setup directory for temporary outputs
    OutGTiff_LM = os.path.join(projdir, 'LANDMASK.tif')
    OutGTiff_ST = os.path.join(projdir, 'SCT_DOM.tif')
    OutGTiff_SB = os.path.join(projdir, 'SCB_DOM.tif')

    # Build and then ensure the output file exists:
    build_geogrid_raster(inNC, 'LANDMASK', OutGTiff_LM, out_Grid_fmt=RasterDriver)
    assert(os.path.exists(OutGTiff_LM))

    # Open the Landmask to use as a mask array
    LM_arr, LM_ndv= wrfh.return_raster_array(OutGTiff_LM)                # Landmask array
    wrfh.remove_file(OutGTiff_LM)

    # Iterate over soil layers
    for soil_class_layer, inraster in zip(['SCT_DOM', 'SCB_DOM'], [OutGTiff_ST, OutGTiff_SB]):

        build_geogrid_raster(inNC, soil_class_layer, inraster, out_Grid_fmt=RasterDriver)
        assert(os.path.exists(inraster))
        arr, ndv = wrfh.return_raster_array(inraster)
        wrfh.remove_file(inraster)

        # Create a copy
        modRaster = os.path.join(projdir, '{0}_mod.tif'.format(soil_class_layer))
        target_ds = gdal.GetDriverByName(RasterDriver).CreateCopy(modRaster, gdal.Open(inraster, gdalconst.GA_ReadOnly))
        target_ds = None

        # 1) Create a grid with holes of value 0 where SCT_DOM is water (14):
        ds = gdal.Open(modRaster, gdalconst.GA_Update)                          # Open for writing
        band = ds.GetRasterBand(1)
        arr_mod = band.ReadAsArray()
        arr_mod[arr_mod==soil_water_val] = 0                    # Set all water soil type cells to 0
        band.WriteArray(arr_mod)                            # Write the array to the disk
        stats = band.GetStatistics(0,1)                 # Calculate statistics
        ds = band = stats = arr_mod = None

        # 2) Run the Euclidean Allocation tool to fill in the 0 gaps with nearest neighbor values.
        EA_output = os.path.join(projdir, '{0}_EA.tif'.format(soil_class_layer))
        wbt.euclidean_allocation(modRaster, EA_output)
        EA_arr, EA_ndv= wrfh.return_raster_array(EA_output)                # Euclidean allocation array
        wrfh.remove_file(EA_output)

        # 3) If requested, create a euclidean distance raster to that we can apply a
        # limit to the gap-filling process
        if distance_calc:
            cell_dist = distance * float(ncDS.DX)
            print('  Using cell distance of {0} map units'.format(cell_dist))

            ED_output = os.path.join(projdir, '{0}_ED.tif'.format(soil_class_layer))
            wbt.euclidean_distance(modRaster, ED_output)
            ED_arr, ED_ndv = wrfh.return_raster_array(ED_output)                # Euclidean distance array
            wrfh.remove_file(ED_output)
        wrfh.remove_file(modRaster)

        # Modify the input to create a grid which samples nearest non-water land
        # cell for all inland water classes on the soil category grid

        # Fill all gaps with the result from Euclidean Allocation
        arr[LM_arr==1] = EA_arr[LM_arr==1]

        # If asking for a distance threshold, fill all other gaps with default value
        if distance_calc:
            arr[numpy.logical_and(LM_arr==1, ED_arr>cell_dist)] = fillsoiltyp
            del ED_arr, ED_ndv

        # Set all water areas (accoridng to landmask) to the soil water value
        arr[LM_arr==0] = soil_water_val

        # Reverse the order of any arrays if necessary
        ind = wrfh.flip_dim(['y', 'x'], DimToFlip='y')     # Assumed gdal array order returned by 'return_raster_array' is y,x
        arr = arr[tuple(ind)]                               # Index the input array

        # Replace what is in GOEGRID with these values
        ncDS[soil_class_layer][:] = arr
        del arr, ndv, EA_arr, EA_ndv
    del LM_arr, LM_ndv

    # Delete all temporary files
    shutil.rmtree(projdir)

    # Output file to disk
    encoding = {varname:ncDS[varname].encoding for varname in list(ncDS.variables.keys())}
    for key, val in encoding.items():
        val['_FillValue'] = None

    ncDS.to_netcdf(outFile, mode='w', format=outNCType, encoding=encoding)
    ncDS.close()
    del encoding, ncDS
    print('Output netCDF file: {0}'.format(outFile))
    print('Process completed in {0:3.2f} seconds.'.format(time.time()-tic))

if __name__ == '__main__':
    print('Script initiated at {0}'.format(time.ctime()))
    tic = time.time()

    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_Geogrid",
                        type=lambda x: is_valid_file(parser, x),
                        required=True,
                        help="Path to WPS geogrid (geo_em.d0*.nc) file [REQUIRED]")
    parser.add_argument("-o",
                        dest="out_geogrid",
                        default='./geo_em_altered_soiltypes.nc',
                        required=True,
                        help='Output "geogrid" file.')

    # If no arguments are supplied, print help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    # Resolve relative paths to absolute paths
    args.in_Geogrid = os.path.abspath(args.in_Geogrid)
    args.out_geogrid = os.path.abspath(args.out_geogrid)

    ##    in_Geogrid = r"C:\Users\ksampson\Desktop\NWM\NWM_Alaska\WPS\WPS_Output\geo_em.d02.20210419_nlcd2016_snow.nc"
    ##    out_geogrid = r"C:\Users\ksampson\Desktop\NWM\NWM_Alaska\WPS\WPS_Output\geo_em.d02.20210419_nlcd2016_snow_modifiedSoils.nc"
    ##    update_geogrid_soils(in_Geogrid, out_geogrid, outNCType=outNCType, distance_calc=True)

    update_geogrid_soils(args.in_Geogrid, args.out_geogrid, outNCType=outNCType, distance_calc=True)
    print('  Process completed in {0:3.2f} seconds'.format(time.time()-tic))