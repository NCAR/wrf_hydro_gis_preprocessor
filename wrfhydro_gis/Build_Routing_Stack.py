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

descText = 'This is a program to perform the full routing-stack GIS pre-processing' + \
            'for WRF-Hydro. The inputs will be related to the domain, the desired ' + \
            'routing nest factor, and other options and parameter values. The output ' + \
            'will be a routing stack zip file with WRF-Hydro domain and parameter files. '

# --- Import Modules --- #

# Import Python Core Modules
import os
import sys
import time
import shutil
import copy
from distutils.version import LooseVersion
import argparse
from argparse import ArgumentParser
import platform                                                                 # Added 8/20/2020 to detect OS

# Import Additional Modules
import netCDF4
import osgeo

try:
    if LooseVersion(osgeo.__version__) > LooseVersion('3.0.0'):
        from osgeo import osr
        from osgeo import gdal
        from osgeo.gdal_array import *
    else:
        import osr
        import gdal
        from gdal_array import *
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

# Import function library into namespace. Must exist in same directory as this script.
import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox

# --- End Import Modules --- #

# Add Proj directory to path
conda_env_path = os.path.join(os.path.dirname(sys.executable))
if platform.system() == 'Windows':
    internal_datadir = os.path.join(conda_env_path, "Library", "share", "proj")
elif platform.system() == 'Linux':
    internal_datadir = os.path.join(os.path.dirname(conda_env_path), "share", "proj")
os.environ["PROJ_LIB"] = internal_datadir

# --- Global Variables --- #

# Provide the default groundwater basin generation method.
# Options: ['FullDom basn_msk variable', 'FullDom LINKID local basins', 'Polygon Shapefile or Feature Class']
defaultGWmethod = 'FullDom LINKID local basins'
#defaultGWmethod = 'FullDom basn_msk variable'
GW_with_Stack = True                                                            # Switch for building default groundwater inputs with any routing stack

# Processing Notes to insert into output netCDF global attributes. Provide any documentation here.
processing_notes_SM = '''Created: {0}'''.format(time.ctime())                   # Processing notes for Spatial Metdata files
processing_notesFD = '''Created: {0}'''.format(time.ctime())                    # Processing notes for the FULLDOM (Routing Grid) file

# --- DO NOT EDIT BELOW THIS LINE --- #

# Parameter default values
default_regridFactor = 10                                                       # Regridding factor
default_ovroughrtfac_val = 1.0
default_retdeprtfac_val = 1.0
default_threshold = 200
default_lksatfac_val = wrfh.lksatfac_val

# Script options
runGEOGRID_STANDALONE = True                                                    # Switch for testing the GEOGRID STANDALONE Pre-processing workflow
cleanUp = True                                                                 # Switch to keep all temporary files (for troubleshooting)

# Methods test switches
coordMethod1 = True                                                             # Interpolate GEOGRID latitude and longitude coordinate arrays
coordMethod2 = False                                                            # Transform coordinate pairs at each grid cell from projected to geocentric

# Variables derived from function script
out_Grid_fmt = wrfh.RasterDriver

#outNCType = 'NETCDF3_64BIT'                                                     # Set output netCDF format for spatial metdata files. This was the default before 7/31/2018
outNCType = 'NETCDF4_CLASSIC'                                                   # Define the output netCDF version for RouteLink.nc and LAKEPARM.nc

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

# Default temporary output file names
mosprj_name = 'mosaicprj.tif'                                                   # Default regridded input DEM if saved to disk

# Default name for the output routing stack zip file
outZipDefault = 'WRF_Hydro_routing_grids.zip'                                   # Default output routing stack zip file name if not provided by user
defaltGeogrid = 'geo_em.d01.nc'                                                 # Default input geogrid file name if not provided by user

# --- End Global Variables --- #

# --- Functions --- #
def is_valid_file(parser, arg):
    # https://stackoverflow.com/questions/11540854/file-as-command-line-argument-for-argparse-error-message-if-argument-is-not-va
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return str(arg)

def GEOGRID_STANDALONE(inGeogrid,
                        regridFactor,
                        inDEM,
                        projdir,
                        threshold,
                        out_zip,
                        in_csv = '',
                        basin_mask = False,
                        routing = False,
                        varList2D = [],
                        in_lakes = '',
                        GW_with_Stack = True,
                        in_GWPolys = None,
                        ovroughrtfac_val = 1.0,
                        retdeprtfac_val = 1.0,
                        lksatfac_val = 1000.0,
                        startPts = None):
    '''
    This function will validate input parameters and attempt to run the full routing-
    stack GIS pre-processing for WRF-Hydro. The inputs will be related to the domain,
    the desired routing nest factor, and other options and parameter values. The
    output will be a routing stack zip file with WRF-Hydro domain and parameter
    files.
    '''

    global defaultGWmethod
    tic1 = time.time()

    # Print information provided to this function
    for key, value in locals().items():
        if callable(value) and value.__module__ == __name__:
            print('      {0}: {1}'.format(key, value))

    # Set some switches
    if os.path.exists(in_csv):
        AddGages = True
        print('    Forecast points provided.')
    else:
        AddGages = False

    if routing:
        print('  Reach-based routing files will be created.')
        varList2D.append(['LINKID', 'i4', ''])
    else:
        print('  Reach-based routing files will not be created.')

    # Step 1 - Georeference geogrid file
    rootgrp = netCDF4.Dataset(inGeogrid, 'r')                                   # Establish an object for reading the input NetCDF files
    globalAtts = rootgrp.__dict__                                               # Read all global attributes into a dictionary
    if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
        rootgrp.set_auto_mask(False)                                            # Change masked arrays to old default (numpy arrays always returned)
    coarse_grid = wrfh.WRF_Hydro_Grid(rootgrp)                                  # Instantiate a grid object
    fine_grid = copy.copy(coarse_grid)                                          # Copy the grid object for modification
    fine_grid.regrid(regridFactor)                                              # Regrid to the coarse grid
    print('    Created projection definition from input NetCDF GEOGRID file.')
    print('    Proj4: {0}'.format(coarse_grid.proj4))                           # Print Proj.4 string to screen
    print('    Coarse grid GeoTransform: {0}'.format(coarse_grid.GeoTransformStr()))        # Print affine transformation to screen.
    print('    Coarse grid extent [Xmin, Ymin, Xmax, Ymax]: {0}'.format(coarse_grid.grid_extent()))      # Print extent to screen.
    print('    Fine grid extent [Xmin, Ymin, Xmax, Ymax]:   {0}'.format(fine_grid.grid_extent()))      # Print extent to screen.

    # Build output raster from numpy array of the GEOGRID variable requested. This will be used as a template later on
    LU_INDEX = coarse_grid.numpy_to_Raster(wrfh.flip_grid(rootgrp.variables['LU_INDEX'][:]))

    # Create spatial metadata file for GEOGRID/LDASOUT grids
    out_nc1 = os.path.join(projdir, wrfh.LDASFile)
    rootgrp1 = netCDF4.Dataset(out_nc1, 'w', format=outNCType)                  # wrf_hydro_functions.outNCType)
    rootgrp1, grid_mapping = wrfh.create_CF_NetCDF(coarse_grid, rootgrp1, projdir,
            notes=processing_notes_SM) # addLatLon=True, latArr=latArr, lonArr=lonArr)
    for item in wrfh.Geogrid_MapVars + ['DX', 'DY']:
        if item in globalAtts:
            rootgrp1.setncattr(item, globalAtts[item])
    rootgrp1.close()
    del rootgrp1

    # Step 3 - Create high resolution topography layers
    in_DEM = gdal.Open(inDEM, 0)                                                # Open with read-only mode
    outDEM = os.path.join(projdir, mosprj_name)
    mosprj = fine_grid.project_to_model_grid(in_DEM, saveRaster=True, OutGTiff=outDEM, resampling=gdal.GRA_Bilinear)
    in_DEM = mosprj = None

    # Build latitude and longitude arrays for Fulldom_hires netCDF file
    if coordMethod1:
        print('  Deriving geocentric coordinates on routing grid from bilinear interpolation of geogrid coordinates.')
        # Build latitude and longitude arrays for GEOGRID_LDASOUT spatial metadata file
        latArr = wrfh.flip_grid(rootgrp.variables['XLAT_M'][:])                 # Extract array of GEOGRID latitude values
        lonArr = wrfh.flip_grid(rootgrp.variables['XLONG_M'][:])                # Extract array of GEOGRID longitude values

        # Resolve any remaining issues with masked arrays. Happens in the ArcGIS pre-processing tools for python 2.7.
        if numpy.ma.isMA(lonArr):
            lonArr = lonArr.data
        if numpy.ma.isMA(latArr):
            latArr = latArr.data

        # Method 1: Use GEOGRID latitude and longitude fields and resample to routing grid
        latRaster1 = coarse_grid.numpy_to_Raster(latArr)                        # Build raster out of GEOGRID latitude array
        lonRaster1 = coarse_grid.numpy_to_Raster(lonArr)                        # Build raster out of GEOGRID longitude array

        latRaster2 = fine_grid.project_to_model_grid(latRaster1)                # Regrid from GEOGRID resolution to routing grid resolution
        lonRaster2 = fine_grid.project_to_model_grid(lonRaster1)                # Regrid from GEOGRID resolution to routing grid resolution
        latRaster1 = lonRaster1 = None                                          # Destroy rater objects
        latArr2 = BandReadAsArray(latRaster2.GetRasterBand(1))                  # Read into numpy array
        lonArr2 = BandReadAsArray(lonRaster2.GetRasterBand(1))                  # Read into numpy array
        latRaster2 = lonRaster2 = None                                          # Destroy raster objects
        del latArr, lonArr, latRaster1, lonRaster1, latRaster2, lonRaster2

    elif coordMethod2:
        print('  Deriving geocentric coordinates on routing grid from direct transformation geogrid coordinates.')
        # Method 2: Transform each point from projected coordinates to geocentric coordinates
        wgs84_proj = osr.SpatialReference()                                     # Build empty spatial reference object
        wgs84_proj.ImportFromProj4(wrfh.wgs84_proj4)                            # Imprort from proj4 to avoid EPSG errors (4326)
        xmap, ymap = fine_grid.getxy()                                          # Get x and y coordinates as numpy array
        latArr2, lonArr2 = wrfh.ReprojectCoords(xmap, ymap, coarse_grid.proj, wgs84_proj)  # Transform coordinate arrays
        del xmap, ymap, wgs84_proj

    # Create FULLDOM file
    out_nc2 = os.path.join(projdir, wrfh.FullDom)
    rootgrp2 = netCDF4.Dataset(out_nc2, 'w', format=outNCType)                  # wrf_hydro_functions.outNCType)
    rootgrp2, grid_mapping = wrfh.create_CF_NetCDF(fine_grid, rootgrp2, projdir,
            notes=processing_notesFD, addVars=varList2D, addLatLon=True,
            latArr=latArr2, lonArr=lonArr2)
    del latArr2, lonArr2

    # Add some global attribute metadata to the Fulldom file, including relevant WPS attributes for defining the model coordinate system
    rootgrp2.geogrid_used = inGeogrid                                           # Paste path of geogrid file to the Fulldom global attributes
    rootgrp2.DX = fine_grid.DX                                                  # Add X resolution as a global attribute
    rootgrp2.DY = -fine_grid.DY                                                 # Add Y resolution as a global attribute
    for item in wrfh.Geogrid_MapVars:
        if item in globalAtts:
            rootgrp2.setncattr(item, globalAtts[item])
    rootgrp.close()                                                             # Close input GEOGRID file
    del item, globalAtts, rootgrp

    # Process: Resample LU_INDEX grid to a higher resolution
    LU_INDEX2 = fine_grid.project_to_model_grid(LU_INDEX, fine_grid.DX, fine_grid.DY, resampling=gdal.GRA_NearestNeighbour)
    rootgrp2.variables['landuse'][:] = BandReadAsArray(LU_INDEX2.GetRasterBand(1))          # Read into numpy array
    LU_INDEX = None                                                             # Destroy raster object
    print('    Process: landuse written to output netCDF.')
    del LU_INDEX, LU_INDEX2

    ##        # Step X(a) - Test to match LANDMASK - Only used for areas surrounded by water (LANDMASK=0)
    ##        mosprj2, loglines = wrfh.adjust_to_landmask(mosprj, LANDMASK, coarse_grid.proj, projdir, 'm')
    ##        outtable.writelines("\n".join(loglines) + "\n")
    ##        del LANDMASK

    # Step 4 - Hyrdo processing functions -- Whitebox
    rootgrp2, fdir, fac, channelgrid, fill, order = wrfh.WB_functions(rootgrp2, outDEM,
            projdir, threshold, ovroughrtfac_val, retdeprtfac_val, lksatfac_val, startPts=startPts)
    if cleanUp:
        wrfh.remove_file(outDEM)                                                # Delete output DEM from disk

    # If the user provides forecast points as a CSV file, alter outputs accordingly
    if AddGages:
        if os.path.exists(in_csv):
            rootgrp2 = wrfh.forecast_points(in_csv, rootgrp2, basin_mask, projdir,
                        fine_grid.DX, fine_grid.WKT, fdir, fac, channelgrid)    # Forecast point processing

    # Moved 10/9/2017 by KMS to allow masking routing files (LINKID, Route_Link, etc.) to forecast points if requested
    if routing:
        rootgrp2 = wrfh.Routing_Table(projdir, rootgrp2, fine_grid, fdir, channelgrid, fill, order, gages=AddGages)
    if cleanUp:
        wrfh.remove_file(fill)                                                  # Delete fill from disk
        wrfh.remove_file(order)                                                 # Delete order from disk

    gridded = not routing                                                       # Flag for gridded routing
    if os.path.exists(in_lakes):
        #pass
        # Alter Channelgrid for reservoirs and build reservoir inputs
        print('    Reservoir polygons provided. Lake routing will be activated.')
        rootgrp2 = wrfh.add_reservoirs(rootgrp2, projdir, fac, in_lakes, fine_grid, Gridded=gridded)
    rootgrp2.close()                                                            # Close Fulldom_hires.nc file
    del rootgrp2

    # Build groundwater files
    if GW_with_Stack:
        if in_GWPolys is not None:
            if os.path.exists(in_GWPolys):
                print('    Groundwater basin boundary polygons provided. Delineating groundwater basins from these polygons.')
                defaultGWmethod = 'Polygon Shapefile or Feature Class'
        GWBasns = wrfh.build_GW_Basin_Raster(out_nc2, projdir, defaultGWmethod, channelgrid, fdir, fine_grid, in_Polys=in_GWPolys)
        wrfh.build_GW_buckets(projdir, GWBasns, coarse_grid, Grid=True)
        GWBasns = None

    if cleanUp:
        wrfh.remove_file(fdir)                                                  # Delete fdir from disk
        wrfh.remove_file(fac)                                                   # Delete fac from disk
        wrfh.remove_file(channelgrid)                                           # Delete channelgrid from disk
    if routing:
        wrfh.remove_file(os.path.join(projdir, wrfh.stream_id))

    # Copmress (zip) the output directory
    zipper = wrfh.zipUpFolder(projdir, out_zip, nclist)
    print('Built output .zip file in {0: 3.2f} seconds.'.format(time.time()-tic1))  # Diagnotsitc print statement

    # Delete all temporary files
    if cleanUp:
        shutil.rmtree(projdir)

# --- End Functions --- #

# --- Main Codeblock --- #
if __name__ == '__main__':
    print('Script initiated at {0}'.format(time.ctime()))
    tic = time.time()

    # Setup the input arguments
    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_Geogrid",
                        type=lambda x: is_valid_file(parser, x),
                        required=True,
                        help="Path to WPS geogrid (geo_em.d0*.nc) file [REQUIRED]")
    parser.add_argument("--CSV",
                        dest="in_CSV",
                        type=lambda x: is_valid_file(parser, x),
                        default=None,
                        help="Path to input forecast point CSV file [OPTIONAL]")
    parser.add_argument("-b",
                        dest="basin_mask",
                        type=bool,
                        default=False,
                        help="Mask CHANNELGRID variable to forecast basins? [True/False]. default=False")
    parser.add_argument("-r",
                        dest="RB_routing",
                        type=bool,
                        default=False,
                        help="Create reach-based routing (RouteLink) files? [True/False]. default=False")
    parser.add_argument("-l",
                        dest="in_reservoirs",
                        type=lambda x: is_valid_file(parser, x),
                        default=None,
                        help="Path to reservoirs shapefile or feature class [OPTIONAL]. If -l is TRUE, this is required.")
    parser.add_argument("-d",
                        dest="inDEM",
                        type=lambda x: is_valid_file(parser, x),
                        default='',
                        required=True,
                        help="Path to input high-resolution elevation raster [REQUIRED]")
    parser.add_argument("-R",
                        dest="cellsize",
                        type=int,
                        default=default_regridFactor,
                        help="Regridding (nest) Factor. default=10")
    parser.add_argument("-t",
                        dest="threshold",
                        type=int,
                        default=default_threshold,
                        help="Number of routing grid cells to define stream. default=200")
    parser.add_argument("-o",
                        dest="out_zip_file",
                        default='./{0}'.format(outZipDefault),
                        help="Output routing stack ZIP file")
    parser.add_argument("-O",
                        dest="ovroughrtfac_val",
                        type=float,
                        default=default_ovroughrtfac_val,
                        help="OVROUGHRTFAC value. default=1.0")
    parser.add_argument("-T",
                        dest="retdeprtfac_val",
                        type=float,
                        default=default_retdeprtfac_val,
                        help="RETDEPRTFAC value. default=1.0")
    parser.add_argument("--starts",
                        dest="channel_starts",
                        type=lambda x: is_valid_file(parser, x),
                        default=None,
                        help="Path to channel initiation points feature class. Must be 2D point type. [OPTIONAL]")
    parser.add_argument("--gw",
                        dest="gw_polys",
                        type=lambda x: is_valid_file(parser, x),
                        default=None,
                        help="Path to groundwater polygons feature class [OPTIONAL]")

    # If no arguments are supplied, print help message
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    # Handle printing to user the default variable name
    print('  Parameter values that have not been altered from script default values:')
    if args.basin_mask == all_defaults["basin_mask"]:
        print('    Using default basin mask setting: {0}'.format(all_defaults["basin_mask"]))
    if args.RB_routing == all_defaults["RB_routing"]:
        print('    Using default reach-based routing setting: {0}'.format(all_defaults["RB_routing"]))
    if args.cellsize == all_defaults["cellsize"]:
        print('    Using default regridding factor: {0}'.format(all_defaults["cellsize"]))
    if args.threshold == all_defaults["threshold"]:
        print('    Using default stream initiation threshold: {0}'.format(all_defaults["threshold"]))
    if args.out_zip_file == all_defaults["out_zip_file"]:
        print('    Using default output location: {0}'.format(all_defaults["out_zip_file"]))
    if args.ovroughrtfac_val == all_defaults["ovroughrtfac_val"]:
        print('    Using default OVROUGHRTFAC parameter value: {0}'.format(all_defaults["ovroughrtfac_val"]))
    if args.retdeprtfac_val == all_defaults["retdeprtfac_val"]:
        print('    Using default RETDEPRTFAC parameter value: {0}'.format(all_defaults["retdeprtfac_val"]))

    # Handle unsupported configurations
    if args.RB_routing and args.in_reservoirs is not None:
        print('  Reach-based routing with reservoirs configuration not currently supported in this version of the GIS Pre-processing tools. Try the ArcGIS version.')
        print('Exiting.')
        raise SystemExit

    # This block allows us to continue to check for a valid file path while allowing the script later to avoid a NoneType error.
    args.in_Geogrid = os.path.abspath(args.in_Geogrid)                          # Obtain absolute path for required input file.
    args.inDEM = os.path.abspath(args.inDEM)                                    # Obtain absolute path for required input file.
    args.out_zip_file = os.path.abspath(args.out_zip_file)                      # Obtain absolute path for required output file.
    if not args.in_reservoirs:
        args.in_reservoirs = ''
    else:
        args.in_reservoirs = os.path.abspath(args.in_reservoirs)                # Obtain absolute path for optional input file.
    if args.in_CSV == None:
        args.in_CSV = ''
    else:
        args.in_CSV = os.path.abspath(args.in_CSV)                              # Obtain absolute path for optional input file.
    if args.channel_starts != None:
        args.channel_starts = os.path.abspath(args.channel_starts)              # Obtain absolute path for optional input file.
    if args.gw_polys is not None:
        args.gw_polys = os.path.abspath(args.gw_polys)                          # Obtain absolute path for optional input file.

    if runGEOGRID_STANDALONE:

        # Configure logging
        logfile = args.out_zip_file.replace('.zip', '.log')
        tee = wrfh.TeeNoFile(logfile, 'w')

        # Print information to screen
        print('  Values that will be used in building this routing stack:')
        print('    Input WPS Geogrid file: {0}'.format(args.in_Geogrid))
        print('    Forecast Point CSV file: {0}'.format(args.in_CSV))
        print('    Mask CHANNELGRID variable to forecast basins?: {0}'.format(args.basin_mask))
        print('    Create reach-based routing (RouteLink) files?: {0}'.format(args.RB_routing))
        print('    Lake polygon feature class: {0}'.format(args.in_reservoirs))
        print('    Input high-resolution DEM: {0}'.format(args.inDEM))
        print('    Regridding factor: {0}'.format(args.cellsize))
        print('    Stream initiation threshold: {0}'.format(args.threshold))
        print('    OVROUGHRTFAC parameter value: {0}'.format(args.ovroughrtfac_val))
        print('    RETDEPRTFAC parameter value: {0}'.format(args.retdeprtfac_val))
        print('    Input channel initiation start point feature class: {0}'.format(args.channel_starts))
        print('    Input groundwater basin polygons: {0}'.format(args.gw_polys))
        print('    Output ZIP file: {0}'.format(args.out_zip_file))

        # Create scratch directory for temporary outputs
        projdir = os.path.join(os.path.dirname(args.out_zip_file), 'scratchdir')
        projdir = os.path.abspath(projdir)
        if os.path.exists(projdir):
            shutil.rmtree(projdir)
        os.makedirs(projdir)

        # Run pre-process
        print('  Running Process GEOGRID function')
        GEOGRID_STANDALONE(args.in_Geogrid,
                            args.cellsize,
                            args.inDEM,
                            projdir,
                            args.threshold,
                            args.out_zip_file,
                            in_csv = args.in_CSV,
                            basin_mask = args.basin_mask,
                            routing = args.RB_routing,
                            varList2D = varList2D,
                            in_lakes = args.in_reservoirs,
                            GW_with_Stack = GW_with_Stack,
                            in_GWPolys = args.gw_polys,
                            ovroughrtfac_val = args.ovroughrtfac_val,
                            retdeprtfac_val = args.retdeprtfac_val,
                            lksatfac_val = default_lksatfac_val,
                            startPts = args.channel_starts)
        tee.close()
        del tee
    else:
        print('  Will not run Process GEOGRID function. Set global "runGEOGRID_STANDALONE" to True to run.')                                                                  # Should do nothing
    print('Process completed in {0:3.2f} seconds.'.format(time.time()-tic))
# --- End Main Codeblock --- #