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

# --- Import Modules --- #

# Import Python system modules
import sys
sys.dont_write_bytecode = True

# Handle running Whitebox Tools from Python 2.x
if sys.version_info < (3, 0):
    from StringIO import StringIO
else:
    from io import StringIO

# Import Python core modules
import time
import os
import glob
import csv
import zipfile
from zipfile import ZipFile, ZipInfo
from operator import itemgetter                                                 # Added 03/28/2023 Used in the group_min function
import collections                                                              # Added 03/28/2023 Used in the group_min function
from collections import defaultdict                                             # Added 09/03/2015 Needed for topological sorting algorthm
from itertools import takewhile, count                                          # Added 09/03/2015 Needed for topological sorting algorthm
import platform                                                                 # Added 8/20/2020 to detect OS
from packaging.version import parse as LooseVersion                             # To avoid deprecation warnings

# Change any environment variables here
#os.environ["OGR_WKT_PRECISION"] = "5"                                           # Change the precision of coordinates

# Import Additional Modules
import netCDF4
import numpy
import osgeo

from shapely.geometry import LineString                                         # Added 03/28/2023 to support line midpoint geometry generation for RouteLink
from shapely import wkt                                                         # Added 03/28/2023 to support line midpoint geometry generation for RouteLink

try:
    if LooseVersion(osgeo.__version__) > LooseVersion('3.0'):
        from osgeo import gdal
        from osgeo import ogr
        from osgeo import osr
        from osgeo import gdalconst
        from osgeo.gdal_array import *                                          # Assists in using BandWriteArray, BandReadAsArray, and CopyDatasetInfo
    else:
        import gdal
        import ogr
        import osr
        import gdalconst
        from gdal_array import *                                                # Assists in using BandWriteArray, BandReadAsArray, and CopyDatasetInfo
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

# Import whitebox
#from whitebox.WBT.whitebox_tools import WhiteboxTools
from whitebox.whitebox_tools import WhiteboxTools

# Module options
gdal.UseExceptions()                                                            # this allows GDAL to throw Python Exceptions
gdal.PushErrorHandler('CPLQuietErrorHandler')

# Add Proj directory to path
conda_env_path = os.path.join(os.path.dirname(sys.executable))
if platform.system() == 'Windows':
    internal_datadir = os.path.join(conda_env_path, "Library", "share", "proj")
elif platform.system() in ['Linux', 'Darwin']:
    internal_datadir = os.path.join(os.path.dirname(conda_env_path), "share", "proj")
os.environ["PROJ_LIB"] = internal_datadir

# --- End Import Modules --- #

# --- Global Variables --- #

# GDAL/OGR Raster and Vector outupt driver options
RasterDriver = 'GTiff'
VectorDriver = 'ESRI Shapefile'                                                # Output vector file format (OGR driver name)

# Version numbers toa ppend to metadata
WRFH_Version = 5.2                                                              # WRF-Hydro version to be executed using the outputs of this tool
PpVersion = 'Open-Source WRF-Hydro GIS Pre-Processing Tools v0.2.0 (04/2021)'   # WRF-Hydro ArcGIS Pre-processor version to add to FullDom metadata
CFConv = 'CF-1.5'                                                               # CF-Conventions version to place in the 'Conventions' attribute of RouteLink files

# Output netCDF format
outNCType = 'NETCDF4_CLASSIC'                                                   # Define the output netCDF version for RouteLink.nc and LAKEPARM.nc

###################################################
# Default output file names
FullDom = 'Fulldom_hires.nc'                                                    # Default Full Domain routing grid nc file
LDASFile = 'GEOGRID_LDASOUT_Spatial_Metadata.nc'                                # Defualt LDASOUT domain grid nc file
LK_nc = 'LAKEPARM.nc'                                                           # Default Lake parameter table name [.nc]
LK_tbl = 'LAKEPARM.TBL'                                                         # Default Lake parameter table name [.TBL]
RT_nc = 'Route_Link.nc'                                                         # Default Route Link parameter table name
GW_nc = 'GWBUCKPARM.nc'                                                         # Default groundwater bucket parameter table name
GWGRID_nc = 'GWBASINS.nc'
GW_ASCII = 'gw_basns_geogrid.txt'                                               # Default Groundwater Basins ASCII grid output
GW_TBL = 'GWBUCKPARM.TBL'
LakesSHP = 'lakes.shp'                                                          # Default lakes shapefile name
minDepthCSV = 'Lakes_with_minimum_depth.csv'                                    # Output file containing lakes with minimum depth enforced.
basinRaster = 'GWBasins.tif'                                                    # Output file name for raster grid of groundwater bucket locations
###################################################

###################################################
# Global Variables
NoDataVal = -9999                                                               # Default NoData value for gridded variables
walker = 3                                                                      # Number of cells to walk downstream before gaged catchment delineation
LK_walker = 3                                                                   # Number of cells to walk downstream to get minimum lake elevation
z_limit = 1000.0                                                                # Maximum fill depth (z-limit) between a sink and it's pour point. None or float.
x_limit = None                                                                  # Maximum breach length for breaching depressions, in pixels. None or Int/Float
lksatfac_val = 1000.0                                                           # Default LKSATFAC value (unitless coefficient)
###################################################

###################################################
# Channel Routing default parameters for the RouteLink file.
Qi = 0                                                                          # Initial Flow in link (cms)
MusK = 3600                                                                     # Muskingum routing time (s)
MusX = 0.2                                                                      # Muskingum weighting coefficient
n = 0.035                                                                       # Manning's roughness
ChSlp = 0.05                                                                    # Channel Side Slope (%; drop/length)
BtmWdth = 5                                                                     # Bottom Width of Channel (m)
Kc = 0                                                                          # Channel loss parameter (mm/hour), New for v1.2
minSo = 0.001                                                                   # Minimum slope

# Order-based Mannings N values for Strahler orders 1-10
ManningsOrd = True                                                              # Switch to activate order-based Mannings N values
Mannings_Order = {1:0.096,
                    2:0.076,
                    3:0.060,
                    4:0.047,
                    5:0.037,
                    6:0.030,
                    7:0.025,
                    8:0.021,
                    9:0.018,
                    10:0.022}                                                    # Values based on CONUS JTTI research, originally from LR 7/01/2020, confirmed by JMC 6/18/21

# Order-based Channel Side-Slope values for Strahler orders 1-10
ChSSlpOrd = True                                                                # Switch to activate order-based Channel Side-Slope values
Mannings_ChSSlp = {1:0.03,
                    2:0.03,
                    3:0.03,
                    4:0.04,
                    5:0.04,
                    6:0.04,
                    7:0.04,
                    8:0.04,
                    9:0.05,
                    10:0.10}                                                    # Values from LR 7/01/2020

# Order-based Bottom-width values for Strahler orders 1-10
BwOrd = True                                                                    # Switch to activate order-based Bottom-width values
Mannings_Bw = {1:1.6,
               2:2.4,
               3:3.5,
               4:5.3,
               5:7.4,
               6:11.,
               7:14.,
               8:16.,
               9:26.,
               10:110.}                                                         # Values from LR 7/01/2020
###################################################

###################################################
#Default Lake Routing parameters
OrificeC = 0.1                                                                  # Default orifice coefficient (0=closed, 1=open)
OrificA = 1.0                                                                   # Default orifice area (square meters)
WeirC = 0.4                                                                     # Default weir coefficient (0=closed, 1=open)
WeirL = 10.0                                                                    # New default prescribed by D. Yates 5/11/2017 (10m default weir length). Old default weir length (0.0m).
ifd_Val = 0.90                                                                  # Default initial fraction water depth (90%)
minDepth = 1.0                                                                  # Minimum active lake depth for lakes with no elevation variation
ChannelLakeCheck = True                                                         # Ensure (True) lakes intersect the channel network.
dam_length = 10.0                                                               # Default length of the dam, multiplier on weir length
###################################################

###################################################
# Globals that support lake pre-processing
datestr = time.strftime("%Y_%m_%d")                                             # Date string to append to output files
FLID = "ARCID"                                                                  # Field name for the flowline IDs
LakeAssoc = 'WBAREACOMI'                                                        # Field name containing link-to-lake association
NoDownstream = [None, -1, 0]                                                    # A list of possible values indicating that the downstream segment is invalid (ocean, network endpoint, etc.)
LkNodata = 0                                                                    # Set the nodata value for the link-to-lake association
hydroSeq = 'HydroSeq'                                                           # Fieldname that stores a hydrologic sequence (ascending from downstream to upstream)
save_Lake_Link_Type_arr = True                                                  # Switch for saving the Lake_Link_Type_arr array to CSV
###################################################

###################################################
# Default groundwater bucket (GWBUCKPARM) parameters
coeff = 1.0000                                                                  # Bucket model coefficient
expon = 3.000                                                                   # Bucket model exponent
zmax = 50.00                                                                    # Conceptual maximum depth of the bucket
zinit = 10.0000                                                                 # Initial depth of water in the bucket model
Loss = 0                                                                        # Not intended for Community WRF-Hydro
maskGW_Basins = False                                                           # Option to mask the GWBASINS.nc grid to only active channels
addLoss = False                                                                 # Option to add loss function parameter to groundwater buckets. Not intended for Community WRF-Hydro use.s
###################################################

# Dictionaries of GEOGRID projections and projection names
#   See http://www.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap3.htm#_Description_of_the_1
projdict = {1: 'Lambert Conformal Conic',
            2: 'Polar Stereographic',
            3: 'Mercator',
            6: 'Cylindrical Equidistant'}
CF_projdict = {1: "lambert_conformal_conic",
                2: "polar_stereographic",
                3: "mercator",
                6: "latitude_longitude",
                0: "crs"}

# Unify all coordinate system variables to have the same name ("crs"). Ths makes it easier for WRF-Hydro output routines to identify the variable and transpose it to output files
crsVarname = True                                                               # Switch to make all coordinate system variables = "crs" instead of related to the coordinate system name
crsVar = CF_projdict[0]                                                         # Expose this as a global for other functions in other scripts to use
wgs84_proj4 = '+proj=longlat +datum=WGS84 +no_defs'                             # Proj.4 string used to define WGS84 coordinate systems. Could also use EPSG code 4326 if using ImportFromEPSG.

# Point time-series CF-netCDF file coordinate system
'''Note that the point netCDF files are handled using a separate coordinate system than the grids.
This is because input data are usually in WGS84 or some other global spheroidal datum. We treat
these coordinates as though there is no difference between a sphere and a spheroid with respect
to latitude. Thus, we can choose an output coordinate system for the points, although no
transformation is performed. Properly transforming the points back and forth betwen sphere and
spheroid dramatically increases the runtime of the tools, with clear obvious benefit.'''
pointCF = True                                                                  # Switch to turn on CF-netCDF point time-series metadata attributes
pointSR = 4326                                                                  # The spatial reference system of the point time-series netCDF files (RouteLink, LAKEPARM). NAD83=4269, WGS84=4326

# Global attributes for altering the sphere radius used in computations. Do not alter sphere_radius for standard WRF-Hydro simulations
sphere_radius = 6370000.0                                                       # Radius of sphere to use (WRF Default = 6370000.0m)
#wkt_text = "GEOGCS['GCS_Sphere_CUSTOM',DATUM['D_Sphere',SPHEROID['Sphere',%s,0.0]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.99462786704589E-09;0.001;0.001;IsHighPrecision" %sphere_radius

###################################################
# Temporary output file names for Whitebox outputs
fill_depressions = "fill_depressions.tif"
dir_d8 = "dir_d8.tif"
streams = "streams.tif"
strahler = "strahler.tif"
sub_basins = "sub_basins.tif"
snapPour1 = 'snapped_pour_points_1.shp'                                         # Pour points snapped to nearest grid cell center
snapPour2 = 'snapped_pour_points_2.shp'                                         # Pour points snapped with higher tolerance
watersheds = "watersheds.tif"                                                   # Watersheds delineated above pour points
start_pts_temp = 'Projected_start_points.shp'                                   # Channel initiation points projected to model CRS
stream_id = "stream_id.tif"                                                     # Stream link ID raster
streams_vector = "streams.shp"                                                  # Stream vector shapefile
###################################################

###################################################
# Dimension names to be used to identify certain known dimensions
yDims = ['south_north', 'y']
xDims = ['west_east', 'x']
timeDim = ['Time', 'time']
Geogrid_MapVars = ['MAP_PROJ', 'corner_lats', 'corner_lons', 'TRUELAT1', 'TRUELAT2', 'STAND_LON', 'POLE_LAT', 'POLE_LON', 'MOAD_CEN_LAT', 'CEN_LAT']
###################################################

version_number = 'v5.1.2 (8/2020)'                                              # Pre-processing tool version

# --- End Global Variables --- #

# --- Classes --- #
class ZipCompat(ZipFile):
    def __init__(self, *args, **kwargs):
        ZipFile.__init__(self, *args, **kwargs)

    def extract(self, member, path=None):
        if not isinstance(member, ZipInfo):
            member = self.getinfo(member)
        if path is None:
            path = os.getcwd()
        return self._extract_member(member, path)

    def extractall(self, path=None, members=None, pwd=None):
        if members is None:
            members = self.namelist()
        for zipinfo in members:
            self.extract(zipinfo, path)

    def _extract_member(self, member, targetpath):
        if (targetpath[-1:] in (os.path.sep, os.path.altsep)
            and len(os.path.splitdrive(targetpath)[1]) > 1):
            targetpath = targetpath[:-1]
        if member.filename[0] == '/':
            targetpath = os.path.join(targetpath, member.filename[1:])
        else:
            targetpath = os.path.join(targetpath, member.filename)
        targetpath = os.path.normpath(targetpath)
        upperdirs = os.path.dirname(targetpath)
        if upperdirs and not os.path.exists(upperdirs):
            os.makedirs(upperdirs)
        if member.filename[-1] == '/':
            if not os.path.isdir(targetpath):
                os.mkdir(targetpath)
            return targetpath
        target = open(targetpath, "wb")                                         # 5/31/2019: Supporting Python3
        try:
            target.write(self.read(member.filename))
        finally:
            target.close()
        return targetpath

class TeeNoFile(object):
    '''
    Send print statements to a log file:
    http://web.archive.org/web/20141016185743/https://mail.python.org/pipermail/python-list/2007-May/460639.html
    https://stackoverflow.com/questions/11124093/redirect-python-print-output-to-logger/11124247
    '''
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
    def close(self):
        if self.stdout is not None:
            sys.stdout = self.stdout
            self.stdout = None
        if self.file is not None:
            self.file.close()
            self.file = None
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
    def flush(self):
        self.file.flush()
        self.stdout.flush()
    def __del__(self):
        self.close()

class WRF_Hydro_Grid:
    '''
    Class with which to create the WRF-Hydro grid representation. Provide grid
    information to initiate the class, and use getgrid() to generate a grid mesh
    and index information about the intersecting cells.

    Note:  The i,j index begins with (1,1) in the upper-left corner.
    '''
    def __init__(self, rootgrp):

        '''                                                                                  .
        9/24/2019:
            This function will create a georeferenced raster object as well as projection
            definition from an input WPS GEOGRID (geo_em.d0*.nc) or WRF-Hydro Fulldom_hires.nc
            netCDF file.

            See the WPS Documentation for more information:

            http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap3.htm

        '''
        tic1 = time.time()
        # First step: Import and georeference NetCDF file
        print('    WPS netCDF projection identification initiated...')

        corner_index = 13                                                           # 13 = Upper left of the Unstaggered grid

        # Loop through global variables in NetCDF file to gather projection information
        dimensions = rootgrp.dimensions
        globalAtts = rootgrp.__dict__                                               # Read all global attributes into a dictionary
        self.map_pro = globalAtts['MAP_PROJ']                                            # Find out which projection this GEOGRID file is in
        print('    Map Projection: {0}'.format(projdict[self.map_pro]))

        # Collect grid corner XY and DX DY for creating ascii raster later
        if 'corner_lats' in globalAtts:
            corner_lat = globalAtts['corner_lats'][corner_index].astype(numpy.float64)
        if 'corner_lons' in globalAtts:
            corner_lon = globalAtts['corner_lons'][corner_index].astype(numpy.float64)
        if 'DX' in globalAtts:
            self.DX = globalAtts['DX'].astype(numpy.float32)
        if 'DY' in globalAtts:
            self.DY = -globalAtts['DY'].astype(numpy.float32)

        # Collect necessary information to put together the projection file
        if 'TRUELAT1' in globalAtts:
            standard_parallel_1 = globalAtts['TRUELAT1'].astype(numpy.float64)
        if 'TRUELAT2' in globalAtts:
            standard_parallel_2 = globalAtts['TRUELAT2'].astype(numpy.float64)
        if 'STAND_LON' in globalAtts:
            central_meridian = globalAtts['STAND_LON'].astype(numpy.float64)
        if 'POLE_LAT' in globalAtts:
            pole_latitude = globalAtts['POLE_LAT'].astype(numpy.float64)
        if 'POLE_LON' in globalAtts:
            pole_longitude = globalAtts['POLE_LON'].astype(numpy.float64)
        if 'MOAD_CEN_LAT' in globalAtts:
            print('    Using MOAD_CEN_LAT for latitude of origin.')
            latitude_of_origin = globalAtts['MOAD_CEN_LAT'].astype(numpy.float64)
        elif 'CEN_LAT' in globalAtts:
            print('    Using CEN_LAT for latitude of origin.')
            latitude_of_origin = globalAtts['CEN_LAT'].astype(numpy.float64)

        # Check to see if the input netCDF is a WPS-Generated Geogrid file.
        if 'TITLE' in globalAtts and (('GEOGRID' in globalAtts['TITLE']) or ('WRF' in globalAtts['TITLE'])):
            self.isGeogrid = True
        else:
            self.isGeogrid = False
        del globalAtts

        # Handle expected dimension names from either Geogrid or Fulldom
        if 'south_north' in dimensions:
            self.nrows = len(dimensions['south_north'])
        elif 'y' in dimensions:
            self.nrows = len(dimensions['y'])
        if 'west_east' in dimensions:
            self.ncols = len(dimensions['west_east'])
        elif 'x' in dimensions:
            self.ncols = len(dimensions['x'])
        del dimensions

        # Initiate OSR spatial reference object - See http://gdal.org/java/org/gdal/osr/SpatialReference.html
        proj = osr.SpatialReference()

        if self.map_pro == 1:
            # Lambert Conformal Conic
            if 'standard_parallel_2' in locals():
                print('    Using Standard Parallel 2 in Lambert Conformal Conic map projection.')
                proj.SetLCC(standard_parallel_1, standard_parallel_2, latitude_of_origin, central_meridian, 0, 0)
                #proj.SetLCC(double stdp1, double stdp2, double clat, double clong, double fe, double fn)        # fe = False Easting, fn = False Northing
            else:
                proj.SetLCC1SP(latitude_of_origin, central_meridian, 1, 0, 0)       # Scale = 1???
                #proj.SetLCC1SP(double clat, double clong, double scale, double fe, double fn)       # 1 standard parallell

        elif self.map_pro == 2:
            # Polar Stereographic
            ##            phi1 = standard_parallel_1
            ##
            ##            ### Back out the central_scale_factor (minimum scale factor?) using formula below using Snyder 1987 p.157 (USGS Paper 1395)
            ##            ##phi = math.copysign(float(pole_latitude), float(latitude_of_origin))    # Get the sign right for the pole using sign of CEN_LAT (latitude_of_origin)
            ##            ##central_scale_factor = (1 + (math.sin(math.radians(phi1))*math.sin(math.radians(phi))) + (math.cos(math.radians(float(phi1)))*math.cos(math.radians(phi))))/2
            ##
            ##            # Method where central scale factor is k0, Derivation from C. Rollins 2011, equation 1: http://earth-info.nga.mil/GandG/coordsys/polar_stereographic/Polar_Stereo_phi1_from_k0_memo.pdf
            ##            # Using Rollins 2011 to perform central scale factor calculations. For a sphere, the equation collapses to be much  more compact (e=0, k90=1)
            ##            central_scale_factor = (1 + math.sin(math.radians(abs(phi1))))/2        # Equation for k0, assumes k90 = 1, e=0. This is a sphere, so no flattening
            ##            print('        Central Scale Factor: {0}'.format(central_scale_factor))
            ##
            ##            # Adjusted 8/7/2017 based on changes made 4/4/2017. Example: proj1.SetPS(90, -1.5, 1, 0, 0)
            ##            proj.SetPS(pole_latitude, central_meridian, central_scale_factor, 0, 0)
            ##            #proj.SetPS(double clat, double clong, double scale, double fe, double fn)
            ##
            ##            # Alternate method? untested.
            ##            ##proj.SetStereographic(pole_latitude, central_meridian, central_scale_factor, 0, 0)

            ##            # Added 3/30/2020. Sad compromise here. For some reason, Esri will not
            ##            # read the central_scale_factor parameter when Latitude_Of_Origin is
            ##            # 90 or 90.0. Adjust slightly, accepting some error (meters) to allow
            ##            # the resulting files to be properly geolocated in esri products.
            ##            if pole_latitude > 0:
            ##                pole_latitude-=0.0000001
            ##            else:
            ##                pole_latitude+=0.0000001

            # Added 3/30/2020. Using a WKT string instead of a projection object
            # because ArcGIS cannot interpret the scale_factor parameter when constructing
            # the projection definition using proj.SetPS or proj.SetSeterographic.
            ##            Projection_String = ('PROJCS["Sphere_Stereographic",'
            ##                                    'GEOGCS["GCS_Sphere",'
            ##                                    'DATUM["D_Sphere",'
            ##                                    'SPHEROID["Sphere",' + str(sphere_radius) + ',0.0]],'
            ##                                    'PRIMEM["Greenwich",0.0],'
            ##                                    'UNIT["Degree",0.0174532925199433]],'
            ##                                    'PROJECTION["Stereographic"],'
            ##                                    'PARAMETER["False_Easting",0.0],'
            ##                                    'PARAMETER["False_Northing",0.0],'
            ##                                    'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
            ##                                    'PARAMETER["Scale_Factor",' + str(central_scale_factor) + '],'
            ##                                    'PARAMETER["Latitude_Of_Origin",' + str(pole_latitude) + '],'
            ##                                    'UNIT["Meter",1.0]]')

            # 3/30/2020: Testing out utility of providing "Stereographic_North_Pole"
            # or "Stereographic_South_Pole" definitions. This is one variant of  the
            # supported projections for polar sterographic that does not require scale factor.
            # However, this may not work for all WRF-supported aspects of stereographic CRS.
            # https://desktop.arcgis.com/en/arcmap/latest/map/projections/stereographic.htm
            # https://svn.osgeo.org/gdal/trunk/autotest/osr/osr_esri.py
            if pole_latitude > 0:
                pole_orientation = 'North'
            else:
                pole_orientation = 'South'

            # 10/5/2020: Changed parameter "Standard_Parallel_1" from 0.0 to the value of TRUELAT1
            # Appears to work for WPS polar stereographic domains such as NWM AK
            Projection_String = ('PROJCS["Sphere_Stereographic",'
                                    'GEOGCS["GCS_Sphere",'
                                    'DATUM["D_Sphere",'
                                    'SPHEROID["Sphere",' + str(sphere_radius) + ',0.0]],'
                                    'PRIMEM["Greenwich",0.0],'
                                    'UNIT["Degree",0.0174532925199433]],'
                                    'PROJECTION["Stereographic_' + pole_orientation + '_Pole"],'
                                    'PARAMETER["False_Easting",0.0],'
                                    'PARAMETER["False_Northing",0.0],'
                                    'PARAMETER["Central_Meridian",' + str(central_meridian) + '],'
                                    'PARAMETER["Standard_Parallel_1",' + str(standard_parallel_1) + '],'
                                    'UNIT["Meter",1.0]]')

            proj.ImportFromWkt(Projection_String)
            #proj.ImportFromESRI([Projection_String])

        elif self.map_pro == 3:
            # Mercator Projection
            proj.SetMercator(standard_parallel_1, central_meridian, 1, 0, 0)     # Scale = 1???
            #proj.SetMercator(latitude_of_origin, central_meridian, 1, 0, 0)     # Scale = 1???
            #proj.SetMercator(double clat, double clong, double scale, double fe, double fn)

        elif self.map_pro == 6:
            # Cylindrical Equidistant (or Rotated Pole)
            if pole_latitude != float(90) or pole_longitude != float(0):
                # if pole_latitude, pole_longitude, or stand_lon are changed from thier default values, the pole is 'rotated'.
                print('[PROBLEM!] Cylindrical Equidistant projection with a rotated pole is not currently supported.')
                raise SystemExit
            else:
                proj.SetEquirectangular(latitude_of_origin, central_meridian, 0, 0)
                #proj.SetEquirectangular(double clat, double clong, double fe, double fn)
                #proj.SetEquirectangular2(double clat, double clong, double pseudostdparallellat, double fe, double fn)

        # Set Geographic Coordinate system (datum) for projection
        proj.SetGeogCS('Sphere', 'Sphere', '', sphere_radius, 0.0)              # Could try 104128 (EMEP Sphere) well-known?
        #proj.SetGeogCS(String pszGeogName, String pszDatumName, String pszSpheroidName, double dfSemiMajor, double dfInvFlattening)

        # Set the origin for the output raster (in GDAL, usuall upper left corner) using projected corner coordinates
        wgs84_proj = osr.SpatialReference()
        wgs84_proj.ImportFromProj4(wgs84_proj4)

        # Added 11/19/2020 to allow for GDAL 3.0 changes to the order of coordinates in transform
        if int(osgeo.__version__[0]) >= 3:
            # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
            wgs84_proj.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
            proj.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

        transform = osr.CoordinateTransformation(wgs84_proj, proj)
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint_2D(corner_lon, corner_lat)
        point.Transform(transform)
        self.x00 = point.GetX(0)
        self.y00 = point.GetY(0)
        self.proj = proj
        self.WKT = proj.ExportToWkt()
        self.proj4 = proj.ExportToProj4()
        del point, transform, wgs84_proj
        print('    Geo-referencing step completed without error in {0: 3.2f} seconds.'.format(time.time()-tic1))

    def regrid(self, regrid_factor):
        '''
        Change the grid cell spacing while keeping all other grid parameters
        the same.
        '''
        print('      Building sub-grid of model grid.')
        print('        Original grid spacing dx={0}, dy={1}'.format(self.DX, self.DY))
        print('        Original grid size: rows={0}, cols={1}'.format(self.nrows, self.ncols))
        self.DX = float(self.DX)/float(regrid_factor)
        self.DY = float(self.DY)/float(regrid_factor)
        self.nrows = int(self.nrows*regrid_factor)
        self.ncols = int(self.ncols*regrid_factor)
        print('        New grid spacing: dx={0}, dy={1}'.format(self.DX, self.DY))
        print('        New dimensions: rows={0}, cols={1}'.format(self.nrows, self.ncols))
        return self

    def GeoTransform(self):
        '''
        Return the affine transformation for this grid. Assumes a 0 rotation grid.
        (top left x, w-e resolution, 0=North up, top left y, 0 = North up, n-s pixel resolution (negative value))
        '''
        return (self.x00, self.DX, 0, self.y00, 0, self.DY)

    def GeoTransformStr(self):
        return ' '.join([str(item) for item in self.GeoTransform()])

    def getxy(self):
        """
        This function will use the affine transformation (GeoTransform) to produce an
        array of X and Y 1D arrays. Note that the GDAL affine transformation provides
        the grid cell coordinates from the upper left corner. This is typical in GIS
        applications. However, WRF uses a south_north ordering, where the arrays are
        written from the bottom to the top.

        The input raster object will be used as a template for the output rasters.
        """
        print('    Starting Process: Building to XMap/YMap')

        # Build i,j arrays
        j = numpy.arange(self.nrows) + float(0.5)                              # Add 0.5 to estimate coordinate of grid cell centers
        i = numpy.arange(self.ncols) + float(0.5)                               # Add 0.5 to estimate coordinate of grid cell centers

        # col, row to x, y   From https://www.perrygeo.com/python-affine-transforms.html
        x = (i * self.DX) + self.x00
        y = (j * self.DY) + self.y00
        del i, j

        # Create 2D arrays from 1D
        xmap = numpy.repeat(x[numpy.newaxis, :], y.shape, 0)
        ymap = numpy.repeat(y[:, numpy.newaxis], x.shape, 1)
        del x, y
        print('    Conversion of input raster to XMap/YMap completed without error.')
        return xmap, ymap

    def grid_extent(self):
        '''
        Return the grid bounding extent [xMin, yMin, xMax, yMax]
        '''
        xMax = self.x00 + (float(self.ncols)*self.DX)
        yMin = self.y00 + (float(self.nrows)*self.DY)
        return [self.x00, yMin, xMax, self.y00]

    def numpy_to_Raster(self, in_arr, quiet=True, nband=1):
        '''This funciton takes in an input netCDF file, a variable name, the ouput
        raster name, and the projection definition and writes the grid to the output
        raster. This is useful, for example, if you have a FullDom netCDF file and
        the GEOGRID that defines the domain. You can output any of the FullDom variables
        to raster.
        Adapted to use as input 3D arrays.'''
        try:
            # Set up driver for GeoTiff output
            driver = gdal.GetDriverByName('Mem')                                # Write to Memory
            if driver is None:
                print('    {0} driver not available.'.format('Memory'))

            gdaltype = NumericTypeCodeToGDALTypeCode(in_arr.dtype)
            print('    GDAL Data type derived from input array: {0} ({1})'.format(gdaltype, in_arr.dtype))
            if not gdaltype:
                print('    The input numpy array type does not have a compatible GDAL data type. Assuming int32.')
                gdaltype = 5                                                    # Int32
            DataSet = driver.Create('', in_arr.shape[-1], in_arr.shape[-2], nband, gdaltype)
            DataSet.SetProjection(self.WKT)
            DataSet.SetGeoTransform(self.GeoTransform())

            if in_arr.ndim == 2:
                in_arr = in_arr[numpy.newaxis]

            for band in range(nband):
                DataSet.GetRasterBand(band+1).WriteArray(in_arr[band])        # Write the array
                #BandWriteArray(DataSet.GetRasterBand(band+1), band_arr[band])
                stats = DataSet.GetRasterBand(band+1).GetStatistics(0,1)        # Calculate statistics
                #stats = DataSet.GetRasterBand(band+1).ComputeStatistics(0)     # Force recomputation of statistics
            driver = None
        except RuntimeError:
            print('ERROR: Unable to build output raster from numpy array.')
            raise SystemExit
        return DataSet

    def boundarySHP(self, outputFile, DriverName='ESRI Shapefile'):
        '''Build a single-feature rectangular polygon that represents the boundary
        of the WRF/WRF-Hydro domain. '''

        # Now convert it to a vector file with OGR
        tic1 = time.time()
        drv = ogr.GetDriverByName(DriverName)
        if drv is None:
            print('      %s driver not available.' % DriverName)
        else:
            print('      %s driver is available.' % DriverName)
            datasource = drv.CreateDataSource(outputFile)
        if datasource is None:
            print('      Creation of output file failed.\n')
            raise SystemExit

        # Create output polygon vector file
        layer = datasource.CreateLayer('boundary', self.proj, geom_type=ogr.wkbPolygon)
        if layer is None:
            print('        Layer creation failed.\n')
            raise SystemExit
        LayerDef = layer.GetLayerDefn()                                             # Fetch the schema information for this layer

        # Create polygon object that is fully inside the outer edge of the domain
        [xMin, yMin, xMax, yMax] = self.grid_extent()
        ring = ogr.Geometry(type=ogr.wkbLinearRing)
        ring.AddPoint(xMin, yMax)
        ring.AddPoint(xMax, yMax)
        ring.AddPoint(xMax, yMin)
        ring.AddPoint(xMin, yMin)
        ring.AddPoint(xMin, yMax)                                     #close ring
        geometry = ogr.Geometry(type=ogr.wkbPolygon)
        geometry.AssignSpatialReference(self.proj)
        geometry.AddGeometry(ring)

        # Create the feature
        feature = ogr.Feature(LayerDef)                                     # Create a new feature (attribute and geometry)
        feature.SetGeometry(geometry)                                      # Make a feature from geometry object
        layer.CreateFeature(feature)
        print('      Done producing output vector polygon shapefile in {0: 3.2f} seconds'.format(time.time()-tic1))
        datasource = ring = feature = layer = None        # geometry
        return geometry

    def xy_to_grid_ij(self, x, y):
        '''
        This function converts a coordinate in (x,y) to the correct row and column
        on a grid. Code from: https://www.perrygeo.com/python-affine-transforms.html

        Grid indices are 0-based.
        '''
        # x,y to col,row.
        col = int((x - self.x00) / self.DX)
        row = int((y - self.y00) / self.DY)
        return row, col

    def grid_ij_to_xy(self, col, row):
        '''
        This function converts a 2D grid index (i,j) the grid cell center coordinate
        (x,y) in the grid coordinate system.
        Code from: https://www.perrygeo.com/python-affine-transforms.html

        Grid indices are 0-based.
        '''
        # col, row to x, y
        x = (col * self.DX) + self.x00 + self.DX/2.0
        y = (row * self.DY) + self.y00 + self.DY/2.0
        return x, y

    def project_to_model_grid(self, in_raster, saveRaster=False, OutGTiff=None, resampling=gdal.GRA_Bilinear):
        '''
        The second step creates a high resolution topography raster using a hydrologically-
        corrected elevation dataset.

        grid object extent and coordinate system will be respected.
        '''
        tic1 = time.time()
        print('    Raster resampling initiated...')

        te = self.grid_extent()                                         # Target Extent
        print('    The High-resolution dataset will be {0}m'.format(str(self.DX)))

        # Use Warp command
        outDtype = in_raster.GetRasterBand(1).DataType
        OutRaster = gdal.Warp('', in_raster, format='MEM', xRes=self.DX, yRes=self.DY,
                            outputBounds=te, outputBoundsSRS=self.WKT,
                            resampleAlg=resampling, dstSRS=self.WKT,
                            errorThreshold=0.0,
                            outputType=outDtype)
        # Other options to gdal.Warp: dstSRS='EPSG:32610', dstNodata=1, srcNodata=1, outputType=gdal.GDT_Int16
        #   transformerOptions=[ 'SRC_METHOD=NO_GEOTRANSFORM', 'DST_METHOD=NO_GEOTRANSFORM']
        #   width=Xsize_out, height=Ysize_out, targetAlignedPixels=True
        del te

        # Save to disk
        if saveRaster:
            if OutRaster is not None:
                try:
                    target_ds = gdal.GetDriverByName(RasterDriver).CreateCopy(OutGTiff, OutRaster)
                    target_ds = None
                except:
                    pass
        # Finish
        print('    Projected input raster to model grid in {0: 3.2f} seconds.'.format(time.time()-tic1))
        return OutRaster

    def getgrid(self, envelope, layer):
        '''Function with which to create the grid intersecting grid cells based
        on a feature geometry envelope. Initiate the class, and use getgrid() to
        generate a grid mesh and index information about the intersecting cells.

        Note:  The i,j index begins with (1,1) in the Lower Left corner.
        '''

        """Gridder.getgrid() takes as input an OGR geometry envelope, and will
        compute the grid polygons that intersect the evelope, returning a list
        of grid cell polygons along with other attribute information.

        Cell IDs are numbered 1...n
        I-index values are numbered 1...n from the lower-left corner (left to right).
        J-index values are numbered 1...n from the lower-left corner (bottom to top).
        """
        # Calculate the number of grid cells necessary
        xmin, xmax, ymin, ymax = envelope

        # Find the i and j indices
        i0 = int((xmin-self.x00)/self.DX // 1)                              # Floor the value
        j0 = int((ymax-self.y00)/self.DY // 1)                              # Floor the absolute value
        i1 = int((xmax-self.x00)/self.DX // 1)                              # Floor the value
        j1 = int((ymin-self.y00)/self.DY // 1)                              # Floor the absolute value

        # Create a new field on a layer. Add one attribute
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('i_index', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('j_index', ogr.OFTInteger))
        LayerDef = layer.GetLayerDefn()                                         # Fetch the schema information for this layer

        # Build OGR polygon objects for each grid cell in the intersecting envelope
        for x in range(i0, i1+1):
            if x < 0 or x > self.ncols:
                continue
            for y in reversed(range(j0, j1+1)):
                if y < 0 or y > self.nrows:
                    continue
                id1 = (self.nrows*(x+1))-y                                      # This should give the ID of the cell from the lower left corner (1,1)

                # Calculating each grid cell polygon's coordinates
                x0 = self.x00 + (self.DX*x)
                x1 = x0 + self.DX
                y1 = self.y00 - (abs(self.DY)*y)
                y0 = y1 - abs(self.DY)

                # Create ORG geometry polygon object using a ring
                myRing = ogr.Geometry(type=ogr.wkbLinearRing)
                myRing.AddPoint(x0, y1)
                myRing.AddPoint(x1, y1)
                myRing.AddPoint(x1, y0)
                myRing.AddPoint(x0, y0)
                myRing.AddPoint(x0, y1)
                geometry = ogr.Geometry(type=ogr.wkbPolygon)
                geometry.AddGeometry(myRing)

                # Create the feature
                feature = ogr.Feature(LayerDef)                                     # Create a new feature (attribute and geometry)
                feature.SetField('id', id1)
                feature.SetField('cellsize', geometry.Area())
                feature.SetField('i_index', x+1)
                feature.SetField('j_index', (self.nrows)-y)
                feature.SetGeometry(geometry)                                      # Make a feature from geometry object
                layer.CreateFeature(feature)
                geometry = feature = None
                del x0, x1, y1, y0, id1
        return layer

#gridder_obj = Gridder_Layer(WKT, DX, DY, x00, y00, nrows, ncols)

# --- End Classes --- #

# --- Functions --- #

def zipws(zipfile, path, zip, keep, nclist):
    path = os.path.normpath(path)
    for dirpath, dirnames, filenames in os.walk(path):
        for file in filenames:
            if file in nclist:
                if keep:
                    try:
                        zip.write(os.path.join(dirpath, file), os.path.join(os.sep + os.path.join(dirpath, file)[len(path) + len(os.sep):]))
                    except:
                        print('Exception encountered while trying to write to output zip file.')

def zipUpFolder(folder, outZipFile, nclist):
    try:
        zip = zipfile.ZipFile(outZipFile, 'w', zipfile.ZIP_DEFLATED, allowZip64=True)
        zipws(zipfile, str(folder), zip, 'CONTENTS_ONLY', nclist)
        zip.close()
    except RuntimeError:
        print('Exception encountered while trying to write to output zip file.')
        pass

def remove_file(in_file):
    '''
    Remove any individual file using os.remove().
    '''
    if os.path.exists(in_file):
        os.remove(in_file)
    return

def flip_grid(array):
    '''This function takes a two dimensional array (x,y) and flips it up-down to
    correct for the netCDF storage of these grids.'''
    array = array[:, ::-1]                                                     # Flip 2+D grid up-down
    return array

def subset_ncVar(ncVar, times=slice(None), DimToFlip='south_north'):
    '''
    7/7/2020:
    This function will accept a netCDF4 Dataset variable object, and will attempt
    to identify the time, x, and y dimensions. If requested, a dimension can be
    specified that will be reversed. This is typically "south_north" dimesnion.
    Also, a time index or slice may be provided. The time dimension size may
    only be greater than 1 if the variable has only 3 dimensions (time, x, y),
    for example.
    '''

    # Ensure x and y dimensions are in the dimensions of the input dataset
    dimensions = ncVar.dimensions
    assert all([any([dim in dimensions for dim in dims]) for dims in [yDims, xDims]])

    # Construct slice to index entire array in original order
    ind = [slice(None)] * len(dimensions)

    # Find the index for the y dimension
    xDimIdx = [dimensions.index(dim) for dim in dimensions if dim in xDims][0]
    print("    X-dimension: '{0}'.".format(dimensions[xDimIdx]))

    # Find the index for the y dimension
    yDimIdx = [dimensions.index(dim) for dim in dimensions if dim in yDims][0]
    print("    Y-dimension: '{0}'.".format(dimensions[yDimIdx]))

    # Flip y-dimension if necessary
    if DimToFlip in dimensions:
        flipIdx = dimensions.index(DimToFlip)
        ind[flipIdx] = slice(None,None,-1)
        print("    Reversing order of dimension '{0}'".format(dimensions[flipIdx]))
        del flipIdx
    else:
        print("    Requested dimension for reversal not found '{0}'.".format(DimToFlip))

    # Find the index for the time dimension
    timeDimIdx = [dimensions.index(dim) for dim in dimensions if dim in timeDim]
    if len(timeDimIdx) > 0:
        timeDimIdx = timeDimIdx[0]
        print("    Time dimension found: '{0}'.".format(dimensions[timeDimIdx]))

        # Choose either a specific time, all times, or the first time to avoid
        # having >1 dimensions. We don't want a 4th dimension that is size 1.
        if ncVar.shape[timeDimIdx] == 1:
            print('      Time dimension size = 1.')
            ind[timeDimIdx] = 0
        else:
            print('      Found time dimension != 1 [{0}].'.format(ncVar.shape[timeDimIdx]))
            ind[timeDimIdx] = times

            # Set any additional dimensions to 0
            additionals = [dim for dim in dimensions if dim not in set(yDims + xDims + timeDim)]
            for extraDim in additionals:
                print("    Selecting '{}' = 0.".format(extraDim))
                ind[dimensions.index(extraDim)] = 0
    else:
        print('    No time dimension found.')

    # Read the array as requested, reversing y if necessary, and subsetting in time
    print('    Dimensions and indices or slices on those dimensions:')
    for dim,indslice in zip(list(dimensions),ind):
        print('        {0}: {1}'.format(dim,indslice))
    ncArr = ncVar[ind]
    assert len(ncArr.shape) <= 3
    del dimensions, timeDimIdx, ind
    return ncArr

def boundarySHP(in_file, outputFile='', DriverName='MEMORY'):
    '''Build a single-feature rectangular polygon in memory that represents the
    boundary of an input raster file. '''

    # Get pertinent information from input raster file
    ds = gdal.Open(in_file, gdalconst.GA_ReadOnly)
    dst_proj = get_projection_from_raster(ds)
    ncols = ds.RasterXSize
    nrows = ds.RasterYSize
    xMin, DX, xskew, yMax, yskew, DY = ds.GetGeoTransform()
    xMax = xMin + (float(ncols)*DX)
    yMin = yMax + (float(nrows)*DY)
    ds = None

    # Now convert it to a vector file with OGR
    tic1 = time.time()
    drv = ogr.GetDriverByName(DriverName)
    if drv is None:
        print('            %s driver not available.' % DriverName)
    else:
        print('            %s driver is available.' % DriverName)
        datasource = drv.CreateDataSource(outputFile)
    if datasource is None:
        print('            Creation of output file failed.\n')
        raise SystemExit

    # Create output polygon vector file
    layer = datasource.CreateLayer('boundary', dst_proj, geom_type=ogr.wkbPolygon)
    if layer is None:
        print('            Layer creation failed.\n')
        raise SystemExit
    LayerDef = layer.GetLayerDefn()                                             # Fetch the schema information for this layer

    # Create polygon object that is fully inside the outer edge of the domain
    ring = ogr.Geometry(type=ogr.wkbLinearRing)
    ring.AddPoint(xMin, yMax)
    ring.AddPoint(xMax, yMax)
    ring.AddPoint(xMax, yMin)
    ring.AddPoint(xMin, yMin)
    ring.AddPoint(xMin, yMax)                                     #close ring
    geometry = ogr.Geometry(type=ogr.wkbPolygon)
    geometry.AssignSpatialReference(dst_proj)
    geometry.AddGeometry(ring)

    # Create the feature
    feature = ogr.Feature(LayerDef)                                     # Create a new feature (attribute and geometry)
    feature.SetGeometry(geometry)                                      # Make a feature from geometry object
    layer.CreateFeature(feature)
    print('      Done producing output vector polygon shapefile in {0: 3.2f} seconds'.format(time.time()-tic1))
    datasource = myRing = feature = layer = None        # geometry
    return geometry

def numpy_to_Raster(in_arr, proj_in=None, DX=1, DY=-1, x00=0, y00=0, quiet=True):
    '''This funciton takes in an input netCDF file, a variable name, the ouput
    raster name, and the projection definition and writes the grid to the output
    raster. This is useful, for example, if you have a FullDom netCDF file and
    the GEOGRID that defines the domain. You can output any of the FullDom variables
    to raster.'''

    tic1 = time.time()
    try:
        # Set up driver for GeoTiff output
        driver = gdal.GetDriverByName('Mem')                                # Write to Memory
        if driver is None:
            print('    {0} driver not available.'.format('Memory'))

        # Set up the dataset and define projection/raster info
        gdaltype = NumericTypeCodeToGDALTypeCode(in_arr.dtype)
        DataSet = driver.Create('', in_arr.shape[1], in_arr.shape[0], 1, gdaltype) # the '1' is for band 1.
        if proj_in:
            DataSet.SetProjection(proj_in.ExportToWkt())
        DataSet.SetGeoTransform((x00, DX, 0, y00, 0, DY))                      # (top left x, w-e resolution, 0=North up, top left y, 0 = North up, n-s pixel resolution (negative value))

        # 7/16/2020: Check to make sure input is not a masked array. Happens in ArcGIS python 2.7.
        if numpy.ma.isMA(in_arr):
            in_arr = in_arr.data

        #DataSet.GetRasterBand(1).WriteArray(in_arr)                             # Write the array
        BandWriteArray(DataSet.GetRasterBand(1), in_arr)
        stats = DataSet.GetRasterBand(1).GetStatistics(0,1)                     # Calculate statistics
        #stats = DataSet.GetRasterBand(1).ComputeStatistics(0)                  # Force recomputation of statistics
        driver = None

    #except RuntimeError:
    except Exception as ex:
        print('ERROR: Unable to build output raster from numpy array.')
        print(ex)
        raise SystemExit

    # Clear objects and return
    if not quiet:
        print('      Created raster in-memory from numpy array in {0:3.2f} seconds.'.format(time.time()-tic1))
    return DataSet

def get_projection_from_raster(in_raster):
    ''' Get projection from input raster and return.'''
    proj = osr.SpatialReference()
    proj.ImportFromWkt(in_raster.GetProjectionRef())
    return proj

def save_raster(OutGTiff, in_raster, rows, cols, gdaltype, NoData=None, Driver='GTiff'):

    target_ds = gdal.GetDriverByName(Driver).Create(OutGTiff, cols, rows, 1, gdaltype)

    band = in_raster.GetRasterBand(1)
    arr_out = band.ReadAsArray()                                                #Read the data into numpy array

    target_ds.SetGeoTransform(in_raster.GetGeoTransform())
    target_ds.SetProjection(in_raster.GetProjection())
    target_ds.GetRasterBand(1).WriteArray(arr_out)

    if NoData is not None:
        target_ds.GetRasterBand(1).SetNoDataValue(NoDataVal)                    # Set noData

    stats = target_ds.GetRasterBand(1).GetStatistics(0,1)                       # Calculate statistics
    #stats = target_ds.GetRasterBand(1).ComputeStatistics(0)                    # Force recomputation of statistics

    target_ds.FlushCache()                                                      #saves to disk!!
    target_ds = None
    return

def define_projection(input_file, dest_srs):
    '''
    This function will define a coordinate system for an input vector layer by
    overwriting an existing vector layer. It is a workaround to create a layer
    in-place but with specifying the destingation coordinate system and copyin g
    all features to the new layer.

    From:
        https://gis.stackexchange.com/questions/126705/how-to-set-the-spatial-reference-to-a-ogr-layer-using-the-python-api
    '''
    tic1 = time.time()

    driver = ogr.Open(input_file).GetDriver()
    ds_in = driver.Open(input_file, 0)
    input_layer = ds_in.GetLayer()

    # Create in-memory output layer to store output
    drv = ogr.GetDriverByName('MEMORY')                                         # Other options: 'ESRI Shapefile'
    ds_out = drv.CreateDataSource('')                                      # Create the data source. If in-memory, use '' or some other string as the data source name
    dest_layer = ds_out.CreateLayer('', dest_srs,
                input_layer.GetLayerDefn().GetGeomType(),
                ['OVERWRITE=YES', 'GEOMETRY_NAME=geom', 'DIM=2', 'FID=id'])

    # adding fields to new layer
    layer_definition = ogr.Feature(input_layer.GetLayerDefn())
    for i in range(layer_definition.GetFieldCount()):
        dest_layer.CreateField(layer_definition.GetFieldDefnRef(i))

    # adding the features from input to dest
    for i in range(0, input_layer.GetFeatureCount()):
        feature = input_layer.GetFeature(i)
        dest_layer.CreateFeature(feature)

    ds_in = input_layer = None                                                  # Dereference all objects from input, except driver
    driver.DeleteDataSource(input_file)                                         # Delete input file
    out_ds  = driver.CopyDataSource(ds_out, input_file)                         # Copy datasource to the same filename as input file
    out_ds = drv = ds_out = dest_layer = layer_definition = feature = None
    print('    Projection defined for input vector layer in {0:3.2f} seconds.'.format(time.time()-tic1))
    return

def return_raster_array(in_file):
    '''
    Read a GDAL-compatible raster file from disk and return the array of raster
    values as well as the nodata value.
    '''
    ds = gdal.Open(in_file, gdalconst.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    ndv = band.GetNoDataValue()                                                 # Obtain nodata value
    ds = band = None
    return arr, ndv

def array_to_points(in_arr, dtype, GT, proj, NoDataVal=-9999):
    '''
    Build a point feature class for every grid cell that contains a unique value
    in the input array.

    This function is intended to be used to derive pour points from channel pixels
    such as lake outlets.

    Assumes input is a 2D numpy array. Seems to only work with integer field type
    currently.
    '''

    tic1 = time.time()
    valField = 'VALUE'
    xMin, DX, xskew, yMax, yskew, DY = GT

    # Create in-memory output layer to store projected and/or clipped polygons
    drv = ogr.GetDriverByName('MEMORY')                                         # Other options: 'ESRI Shapefile'
    data_source = drv.CreateDataSource('')                                      # Create the data source. If in-memory, use '' or some other string as the data source name
    outLayer = data_source.CreateLayer('', proj, ogr.wkbPoint)               # Create the layer name. Use '' or some other string as the layer name
    outLayer.CreateField(ogr.FieldDefn(valField, dtype))                         # Add a single field to the new layer
    outlayerDef = outLayer.GetLayerDefn()

    # col, row to x, y   From https://www.perrygeo.com/python-affine-transforms.html
    uniques = numpy.unique(in_arr[in_arr!=NoDataVal])                           # Determine unique values
    for idval in uniques:
        locs = numpy.where(in_arr==idval)
        for j, i in zip(locs[0], locs[1]):
            x = (i * DX) + xMin + float(DX/2)
            y = (j * DY) + yMax + float(DY/2)

            # Build output feature
            outFeature = ogr.Feature(outlayerDef)
            outFeature.SetField(valField, int(idval))                                    # Set pixel value attribute
            wkt = "POINT({0} {1})".format(x, y)                                       # create the WKT for the feature using Python string formatting
            point = ogr.CreateGeometryFromWkt(wkt)                                  # Create the point from the Well Known Txt
            outFeature.SetGeometry(point)                                           # Set the feature geometry using the point
            outLayer.CreateFeature(outFeature)                                            # Create the feature in the layer (shapefile)
            outFeature.Destroy()
            outFeature = point = None                                               # Dereference the feature
    outLayer = None
    return data_source

def ReprojectCoords(xcoords, ycoords, src_srs, tgt_srs):
    '''
    Adapted from:
        https://gis.stackexchange.com/questions/57834/how-to-get-raster-corner-coordinates-using-python-gdal-bindings
     Reproject a list of x,y coordinates.
    '''
    tic1 = time.time()

    # Added 11/19/2020 to allow for GDAL 3.0 changes to the order of coordinates in transform
    if int(osgeo.__version__[0]) >= 3:
        # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
        src_srs.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        tgt_srs.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    # Setup coordinate transform
    transform = osr.CoordinateTransformation(src_srs, tgt_srs)

    ravel_x = numpy.ravel(xcoords)
    ravel_y = numpy.ravel(ycoords)
    trans_x = numpy.zeros(ravel_x.shape, ravel_x.dtype)
    trans_y = numpy.zeros(ravel_y.shape, ravel_y.dtype)

    for num,(x,y) in enumerate(zip(ravel_x, ravel_y)):
        x1,y1,z = transform.TransformPoint(x,y)
        trans_x[num] = x1
        trans_y[num] = y1

    # reshape transformed coordinate arrays of the same shape as input coordinate arrays
    trans_x = trans_x.reshape(*xcoords.shape)
    trans_y = trans_y.reshape(*ycoords.shape)
    print('Completed transforming coordinate pairs [{0}] in {1: 3.2f} seconds.'.format(num, time.time()-tic1))
    return trans_x, trans_y

# Function for using forecast points
def FeatToRaster(InputVector, inRaster, fieldname, dtype, NoData=None):
    '''
    This function will take a point shapefile and rasterize it. The point feature
    class must have a field in it with values of 1, which is an optional input
    to the RasterizeLayer function. Currently, this field is named PURPCODE. The
    result is a raster of NoData and 1 values, which is used as the "Input
    Depression Mask Grid" in TauDEM's PitRemove tool.
    '''
    # Python GDAL_RASTERIZE syntax, adatped from:
    #    https://gis.stackexchange.com/questions/212795/rasterizing-shapefiles-with-gdal-and-python

    # Open Raster input
    ds = gdal.Open(inRaster, gdalconst.GA_ReadOnly)

    # Get shapefile information
    in_vector = ogr.Open(InputVector)
    in_layer = in_vector.GetLayer()
    driver = gdal.GetDriverByName('Mem')                                        # Write to Memory
    target_ds = driver.Create('', ds.RasterXSize, ds.RasterYSize, 1, dtype)

    # Copy input raster info to output (SpatialReference, Geotransform, etc)
    CopyDatasetInfo(ds, target_ds)
    band = target_ds.GetRasterBand(1)
    if NoData is not None:
        band.SetNoDataValue(NoData)
    band.FlushCache()
    gdal.RasterizeLayer(target_ds, [1], in_layer, options=["ATTRIBUTE=%s" %fieldname])
    stats = target_ds.GetRasterBand(1).GetStatistics(0,1)                       # Calculate statistics on new raster
    return target_ds

def project_Features(InputVector, outProj, clipGeom=None, geomType=None):
    '''
    This function is intended to project a polygon geometry to a new coordinate
    system. Optionally, the geometries can be clipped to an extent rectangle. If
    this option is chosen, the geometry will be clipped for each intersecting
    polygon.
    '''
    tic1 = time.time()
    trans = False

    # Get input vector information
    in_vect = ogr.Open(InputVector)                                             # Read the input vector file
    in_layer = in_vect.GetLayer()                                               # Get the 'layer' object from the data source
    in_proj = in_layer.GetSpatialRef()                                          # Obtain the coordinate reference object.
    in_LayerDef = in_layer.GetLayerDefn()

    # Added 11/19/2020 to allow for GDAL 3.0 changes to the order of coordinates in transform
    if int(osgeo.__version__[0]) >= 3:
        # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
        in_proj.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        outProj.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    # Check if a coordinate transformation (projection) must be performed
    if not outProj.IsSame(in_proj):
        print('        Input shapefile projection does not match requested output. Transforming.')
        coordTrans = osr.CoordinateTransformation(in_proj, outProj)
        trans = True

    # Create in-memory output layer to store projected and/or clipped features
    drv = ogr.GetDriverByName('memory')                                         # Other options: 'ESRI Shapefile'
    data_source = drv.CreateDataSource('')                                      # Create the data source. If in-memory, use '' or some other string as the data source name
    if geomType is None:
        geomType = in_layer.GetGeomType()                                       # Use the geometry of the input of not specified
    outLayer = data_source.CreateLayer('', outProj, geomType)                   # Create the layer name. Use '' or some other string as the layer name

    # Get field names
    fieldNames = []
    for i in range(in_LayerDef.GetFieldCount()):
        fieldNames.append(in_LayerDef.GetFieldDefn(i).GetName())

    # adding fields to new layer
    layer_definition = ogr.Feature(in_LayerDef)
    for i in range(layer_definition.GetFieldCount()):
        outLayer.CreateField(layer_definition.GetFieldDefnRef(i))
    layer_defininition = None
    outlayerDef = outLayer.GetLayerDefn()

    # Read all features in layer
    for feature in in_layer:
        geometry = feature.GetGeometryRef()                                     # Get the geometry object from this feature
        if trans:
            geometry.Transform(coordTrans)                                      # Transform the geometry
        if clipGeom:
            if clipGeom.Intersects(geometry):
                geometry = geometry.Intersection(clipGeom)                      # Clip the geometry if requested
            else:
                continue                                                        # Go to the next feature (do not copy)
        if not geometry:
            continue                                                            # Trap because some geometries end up as None
        feature.SetGeometry(geometry)                                           # Set output Shapefile's feature geometry
        outLayer.CreateFeature(feature)
        feature = geometry = None                                  # Clear memory
    in_layer.ResetReading()
    outLayer.ResetReading()
    outFeatCount = outLayer.GetFeatureCount()                                   # Get number of features in output layer
    #outLayer = None                                                            # Clear memory

    print('        Number of output features: {0} of {1}'.format(outFeatCount, in_layer.GetFeatureCount()))
    print('      Completed reprojection and-or clipping in {0:3.2f} seconds.'.format(time.time()-tic1))
    in_vect = inlayerDef = in_layer = in_LayerDef = None
    return data_source, outLayer, fieldNames

def raster_to_polygon(in_raster, in_proj, geom_typ=ogr.wkbPolygon):
    '''
    Convert a raster object to a polygon layer.
    '''
    tic1 = time.time()

    # Create temporary polygon vector layer
    ds = ogr.GetDriverByName('MEMORY').CreateDataSource('')
    Layer = ds.CreateLayer('', geom_type=geom_typ, srs=in_proj)   # Use projection from input raster
    Layer.CreateField(ogr.FieldDefn('RASTERVALU', ogr.OFTReal))

    # Get raster band information
    band = in_raster.GetRasterBand(1)                                           # Get raster band 1
    stats = band.ComputeStatistics(0)                                       # Force recomputation of statistics

    # Polygonize the raster and write features to in-memory vector layer
    result = gdal.Polygonize(band, band, Layer, 0, ["8CONNECTED=8"], callback=None)     # With 8-connectedness
    if result != 0:
        print('Polygonize raster failed')
    else:
        print('  Created polygon from input raster in {0: 3.2f} seconds'.format(time.time()-tic1))

    #feature = Layer.GetNextFeature()
    return ds, Layer

def dissolve_polygon_to_multipolygon(inDS, inLayer, fieldname, quiet=True):
    '''
    This function will dissolve the polygons in an input polygon feature layer
    and provide an output in-memory multipolygon layer with the dissolved geometries.
    '''
    tic1 = time.time()
    in_proj = inLayer.GetSpatialRef()

    # Create temporary polygon vector layer
    ds = ogr.GetDriverByName('MEMORY').CreateDataSource('')
    outLayer = ds.CreateLayer('', geom_type=ogr.wkbMultiPolygon, srs=in_proj)   # Use projection from input raster
    inlayerDef = inLayer.GetLayerDefn()                                        # Obtain the layer definition for this layer

    # Copy fields from input vector layer to output
    fieldNames = []                                                             # Build empty list of field names
    for i in range(inlayerDef.GetFieldCount()):
        fieldDef = inlayerDef.GetFieldDefn(i)                                   # Get the field definition for this field
        fieldName =  fieldDef.GetName()                                         # Get the field name for this field
        outLayer.CreateField(ogr.FieldDefn(fieldName, fieldDef.GetType()))      # Create a field in the output that matches the field in the input layer
        fieldNames.append(fieldName)                                            # Add field name to list of field names
    outlayerDef = outLayer.GetLayerDefn()
    inlayerDef = None

    # Set up list of unique lake IDs over which to dissolve singlepart to multiapart polygons
    valuelist = set([feature.GetField(fieldname) for feature in inLayer])       # Get list of unique IDs
    inLayer.ResetReading()                                                      # Reset layer
    for idval in valuelist:
        inLayer.SetAttributeFilter('"%s" = %s' %(fieldname, idval))         # Select the ID from the layer
        polygeom = ogr.Geometry(ogr.wkbMultiPolygon)
        for feature in inLayer:
            polygeom.AddGeometry(feature.GetGeometryRef())
        #polygeom = polygeom.UnionCascaded()
        if not quiet:
            print('  [{0}] Number of features in the original polygon: {1},  multipolygon: {2}'.format(int(idval), inLayer.GetFeatureCount(), polygeom.GetGeometryCount()))

        # Create output Feature
        outFeature = ogr.Feature(outlayerDef)                                   # Create new feature
        outFeature.SetGeometry(polygeom)                                        # Set output Shapefile's feature geometry

        # Fill in fields. All fields in input will be transferred to output
        for fieldname in fieldNames:
            if fieldname == 'AREASQKM':
                outFeature.SetField(fieldname, float(polygeom.Area()/1000000.0))       # Add an area field to re-calculate area
            else:
                outFeature.SetField(fieldname, feature.GetField(fieldname))
        outLayer.CreateFeature(outFeature)                                      # Add new feature to output Layer
        feature = outFeature = polygeom = None                                  # Clear memory
        inLayer.SetAttributeFilter(None)
    inlayerDef = outlayerDef = outLayer = None
    del idval, fieldNames, valuelist, in_proj
    print('    Done dissolving input layer in {0:3.2f} seconds.'.format(time.time()-tic1))
    return ds

def add_CRS_var(rootgrp, sr, map_pro, CoordSysVarName, grid_mapping, PE_string, GeoTransformStr=None):
    '''
    10/13/2017 (KMS):
        This function was added to generalize the creating of a CF-compliant
        coordinate reference system variable. This was modularized in order to
        create CRS variables for both gridded and point time-series CF-netCDF
        files.
    '''
    tic1 = time.time()

    # Scalar projection variable - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
    proj_var = rootgrp.createVariable(CoordSysVarName, 'S1')                    # (Scalar Char variable)
    proj_var.transform_name = grid_mapping                                      # grid_mapping. grid_mapping_name is an alias for this
    proj_var.grid_mapping_name = grid_mapping                                   # for CF compatibility
    proj_var.esri_pe_string = PE_string                                         # For ArcGIS. Not required if esri_pe_string exists in the 2D variable attributes
    #proj_var.spatial_ref = PE_string                                            # For GDAl
    proj_var.long_name = "CRS definition"                                       # Added 10/13/2017 by KMS to match GDAL format
    proj_var.longitude_of_prime_meridian = 0.0                                  # Added 10/13/2017 by KMS to match GDAL format
    if GeoTransformStr is not None:
        proj_var.GeoTransform = GeoTransformStr                                 # For GDAl - GeoTransform array

    # Projection specific parameters - http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/reference/StandardCoordinateTransforms.html
    if map_pro == 1:
        # Lambert Conformal Conic

        # Required transform variables
        proj_var._CoordinateAxes = 'y x'                                            # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_var._CoordinateTransformType = "Projection"
        proj_var.standard_parallel = sr.GetProjParm("standard_parallel_1"), sr.GetProjParm("standard_parallel_2")     # Double
        proj_var.longitude_of_central_meridian = sr.GetProjParm("central_meridian")     # Double. Necessary in combination with longitude_of_prime_meridian?
        proj_var.latitude_of_projection_origin = sr.GetProjParm("latitude_of_origin")   # Double

        # Optional tansform variable attributes
        proj_var.false_easting = sr.GetProjParm("false_easting")                # Double  Always in the units of the x and y projection coordinates
        proj_var.false_northing = sr.GetProjParm("false_northing")              # Double  Always in the units of the x and y projection coordinates
        proj_var.earth_radius = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        proj_var.semi_major_axis = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        proj_var.inverse_flattening = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 2:
        # Polar Stereographic

        # Required transform variables
        proj_var._CoordinateAxes = 'y x'                                            # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_var._CoordinateTransformType = "Projection"
        proj_var.longitude_of_projection_origin = sr.GetProjParm("longitude_of_origin")   # Double - proj_var.straight_vertical_longitude_from_pole = ''
        proj_var.latitude_of_projection_origin = sr.GetProjParm("latitude_of_origin")     # Double
        proj_var.scale_factor_at_projection_origin = sr.GetProjParm("scale_factor")      # Double

        # Optional tansform variable attributes
        proj_var.false_easting = sr.GetProjParm("false_easting")                         # Double  Always in the units of the x and y projection coordinates
        proj_var.false_northing = sr.GetProjParm("false_northing")                       # Double  Always in the units of the x and y projection coordinates
        proj_var.earth_radius = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        proj_var.semi_major_axis = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        proj_var.inverse_flattening = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 3:
        # Mercator

        # Required transform variables
        proj_var._CoordinateAxes = 'y x'                                            # Coordinate systems variables always have a _CoordinateAxes attribute, optional for dealing with implicit coordinate systems
        proj_var._CoordinateTransformType = "Projection"
        proj_var.longitude_of_projection_origin = sr.GetProjParm("central_meridian")   # Double
        proj_var.latitude_of_projection_origin = sr.GetProjParm("latitude_of_origin")     # Double
        proj_var.standard_parallel = sr.GetProjParm("standard_parallel_1")                # Double
        proj_var.earth_radius = sphere_radius                                   # OPTIONAL. Parameter not read by Esri. Default CF sphere: 6371.229 km.
        proj_var.semi_major_axis = sphere_radius                                # Added 10/13/2017 by KMS to match GDAL format
        proj_var.inverse_flattening = float(0)                                  # Added 10/13/2017 by KMS to match GDAL format: Double - optional Lambert Conformal Conic parameter

    elif map_pro == 6:
        # Cylindrical Equidistant or rotated pole

        #http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#appendix-grid-mappings
        # Required transform variables
        #proj_var.grid_mapping_name = "latitude_longitude"                      # or "rotated_latitude_longitude"

        #print('        Cylindrical Equidistant projection not supported.')
        #raise SystemExit
        pass                                                                    # No extra parameters needed for latitude_longitude

    # Added 10/13/2017 by KMS to accomodate alternate datums
    elif map_pro == 0:
        proj_var._CoordinateAxes = 'lat lon'
        proj_var.semi_major_axis = sr.GetSemiMajor()
        proj_var.semi_minor_axis =  sr.GetSemiMinor()
        proj_var.inverse_flattening = sr.GetInvFlattening()
        pass

    # Global attributes related to CF-netCDF
    rootgrp.Conventions = CFConv                                                # Maybe 1.0 is enough?
    return rootgrp

def create_CF_NetCDF(grid_obj, rootgrp, projdir, addLatLon=False, notes='', addVars=[], latArr=None, lonArr=None):
    """This function will create the netCDF file with CF conventions for the grid
    description. Valid output formats are 'GEOGRID', 'ROUTING_GRID', and 'POINT'.
    The output NetCDF will have the XMAP/YMAP created for the x and y variables
    and the LATITUDE and LONGITUDE variables populated from the XLAT_M and XLONG_M
    variables in the GEOGRID file or in the case of the routing grid, populated
    using the getxy function."""

    tic1 = time.time()
    print('      Creating CF-netCDF File.')

    # Build Esri WKT Projection string to store in CF netCDF file
    projEsri = grid_obj.proj.Clone()                                            # Copy the SRS
    projEsri.MorphToESRI()                                                      # Alter the projection to Esri's representation of a coordinate system
    PE_string = projEsri.ExportToWkt().replace("'", '"')                        # INVESTIGATE - this somehow may provide better compatability with Esri products?
    print('        Esri PE String: {0}'.format(PE_string))

    # Find name for the grid mapping
    if CF_projdict.get(grid_obj.map_pro) is not None:
        grid_mapping = CF_projdict[grid_obj.map_pro]
        print('        Map Projection of input raster : {0}'.format(grid_mapping))
    else:
        grid_mapping = 'crs'                                                    # Added 10/13/2017 by KMS to generalize the coordinate system variable names
        print('        Map Projection of input raster (not a WRF projection): {0}'.format(grid_mapping))

    # Create Dimensions
    dim_y = rootgrp.createDimension('y', grid_obj.nrows)
    dim_x = rootgrp.createDimension('x', grid_obj.ncols)

    # Create coordinate variables
    var_y = rootgrp.createVariable('y', 'f8', 'y')                              # (64-bit floating point)
    var_x = rootgrp.createVariable('x', 'f8', 'x')                              # (64-bit floating point)

    # Must handle difference between ProjectionCoordinateSystem and LatLonCoordinateSystem
    if grid_obj.proj.IsGeographic():
        if crsVarname:
            CoordSysVarName = crsVar
        else:
            CoordSysVarName = "LatLonCoordinateSystem"

        # Set variable attributes
        #var_y.standard_name = ''
        #var_x.standard_name = ''
        var_y.long_name = "latitude coordinate"
        var_x.long_name = "longitude coordinate"
        var_y.units = "degrees_north"
        var_x.units = "degrees_east"
        var_y._CoordinateAxisType = "Lat"
        var_x._CoordinateAxisType = "Lon"

    elif grid_obj.proj.IsProjected():
        if crsVarname:
            CoordSysVarName = crsVar
        else:
            CoordSysVarName = "ProjectionCoordinateSystem"
        #proj_units = sr.linearUnitName.lower()                                  # sr.projectionName wouldn't work for a GEOGCS
        proj_units = 'm'                                                        # Change made 11/3/2016 by request of NWC

        # Set variable attributes
        var_y.standard_name = 'projection_y_coordinate'
        var_x.standard_name = 'projection_x_coordinate'
        var_y.long_name = 'y coordinate of projection'
        var_x.long_name = 'x coordinate of projection'
        var_y.units = proj_units                                                # was 'meter', now 'm'
        var_x.units = proj_units                                                # was 'meter', now 'm'
        var_y._CoordinateAxisType = "GeoY"                                      # Use GeoX and GeoY for projected coordinate systems only
        var_x._CoordinateAxisType = "GeoX"                                      # Use GeoX and GeoY for projected coordinate systems only
        var_y.resolution = float(abs(grid_obj.DY))                              # Added 11/3/2016 by request of NWC
        var_x.resolution = float(grid_obj.DX)                                   # Added 11/3/2016 by request of NWC

        # Build coordinate reference system variable
        rootgrp = add_CRS_var(rootgrp, grid_obj.proj, grid_obj.map_pro, CoordSysVarName, grid_mapping, PE_string, grid_obj.GeoTransformStr())

    # For prefilling additional variables and attributes on the same 2D grid, given as a list [[<varname>, <vardtype>, <long_name>],]
    for varinfo in addVars:
        ncvar = rootgrp.createVariable(varinfo[0], varinfo[1], ('y', 'x'))
        ncvar.esri_pe_string = PE_string
        ncvar.grid_mapping = CoordSysVarName
        #ncvar.long_name = varinfo[2]
        #ncvar.units = varinfo[3]

    # Get x and y variables for the netCDF file
    xmap, ymap = grid_obj.getxy()                                               # Get coordinates as numpy array
    var_y[:] = ymap[:,0]                                                        # Assumes even spacing in y across domain
    var_x[:] = xmap[0,:]                                                        # Assumes even spacing in x across domain
    ymap = xmap = None
    del ymap, xmap

    if addLatLon == True:
        print('        Proceeding to add LATITUDE and LONGITUDE variables after {0: 8.2f} seconds.'.format(time.time()-tic1))

        # Populate this file with 2D latitude and longitude variables
        # Latitude and Longitude variables (WRF)
        lat_WRF = rootgrp.createVariable('LATITUDE', 'f4', ('y', 'x'))          # (32-bit floating point)
        lon_WRF = rootgrp.createVariable('LONGITUDE', 'f4', ('y', 'x'))         # (32-bit floating point)
        lat_WRF.long_name = 'latitude coordinate'                               # 'LATITUDE on the WRF Sphere'
        lon_WRF.long_name = 'longitude coordinate'                              # 'LONGITUDE on the WRF Sphere'
        lat_WRF.units = "degrees_north"
        lon_WRF.units = "degrees_east"
        lat_WRF._CoordinateAxisType = "Lat"
        lon_WRF._CoordinateAxisType = "Lon"
        lat_WRF.grid_mapping = CoordSysVarName                                  # This attribute appears to be important to Esri
        lon_WRF.grid_mapping = CoordSysVarName                                  # This attribute appears to be important to Esri
        lat_WRF.esri_pe_string = PE_string
        lon_WRF.esri_pe_string = PE_string

        # Missing value attribute not needed yet
        #missing_val = numpy.finfo(numpy.float32).min                            # Define missing data variable based on numpy
        #lat_WRF.missing_value = missing_val                                     # Float sys.float_info.min?
        #lon_WRF.missing_value = missing_val                                     # Float sys.float_info.min?

        '''Adding the Esri PE String in addition to the CF grid mapping attributes
        is very useful. Esri will prefer the PE string over other CF attributes,
        allowing a spherical datum to be defined. Esri can interpret the coordinate
        system variable alone, but will assume the datum is WGS84. This cannot be
        changed except when using an Esri PE String.'''

        ##    # Create a new coordinate system variable
        ##    LatLonCoordSysVarName = "LatLonCoordinateSystem"
        ##    latlon_var = rootgrp.createVariable(LatLonCoordSysVarName, 'S1')            # (Scalar Char variable)
        ##    latlon_var._CoordinateAxes = 'LATITUDE LONGITUDE'                           # Coordinate systems variables always have a _CoordinateAxes attribute

        # Data variables need _CoodinateSystems attribute
        lat_WRF._CoordinateAxisType = "Lat"
        lon_WRF._CoordinateAxisType = "Lon"
        lat_WRF._CoordinateSystems = CoordSysVarName
        lon_WRF._CoordinateSystems = CoordSysVarName
        ##    lat_WRF._CoordinateSystems = "%s %s" %(CoordSysVarName, LatLonCoordSysVarName)        # For specifying more than one coordinate system
        ##    lon_WRF._CoordinateSystems = "%s %s" %(CoordSysVarName, LatLonCoordSysVarName)        # For specifying more than one coordinate system

        # Populate netCDF variables using input numpy arrays
        lat_WRF[:] = latArr
        lon_WRF[:] = lonArr

    # Global attributes
    rootgrp.GDAL_DataType = 'Generic'
    rootgrp.Source_Software = 'WRF-Hydro GIS Pre-processor {0}'.format(PpVersion)
    rootgrp.proj4 = grid_obj.proj4                                              # Added 3/16/2018 (KMS) to avoid a warning in WRF-Hydro output
    rootgrp.history = 'Created {0}'.format(time.ctime())
    rootgrp.processing_notes = notes
    rootgrp.spatial_ref = PE_string                                             # For GDAl
    print('      netCDF global attributes set after {0: 3.2f} seconds.'.format(time.time()-tic1))

    # Return the netCDF file to the calling script
    return rootgrp, grid_mapping

def build_GW_Basin_Raster(in_nc, projdir, in_method, strm, fdir, grid_obj, in_Polys=None):
    '''
    10/10/2017:
    This function was added to build the groundwater basins raster using a variety
    of methods. The result is a raster on the fine-grid which can be used to create
    groundwater bucket parameter tables in 1D and 2D for input to WRF-Hydro.
    '''

    #(in_nc, in_method, strm, grid_obj, in_Polys) = (inFulldom, defaultGWmethod, channelgrid, fine_grid, in_GWPolys)
    tic1 = time.time()
    print('    Beginning to build 2D groundwater basin inputs')
    print('      Building groundwater inputs using {0}'.format(in_method))

    # Determine which method will be used to generate groundwater bucket grid
    if in_method == 'FullDom basn_msk variable':
        print('        Reading Fulldom_hires for basn_msk variable.')

        # Create a raster layer from the netCDF
        rootgrp = netCDF4.Dataset(in_nc, 'r')                                      # Read-only on FullDom file
        basn_arr = rootgrp.variables['basn_msk'][:]
        if numpy.unique(basn_arr).shape[0] == 1:
            print('    WARNING: No basins are present in Fulldom basn_msk variable.')
        GWBasns = grid_obj.numpy_to_Raster(rootgrp.variables['basn_msk'][:])
        rootgrp.close()
        del rootgrp, basn_arr

    elif in_method == 'FullDom LINKID local basins':
        print('        Generating LINKID grid for building local sub-basins.')

        # Whitebox options for running Whitebox in a full workflow
        wbt = WhiteboxTools()
        esri_pntr = True
        wbt.verbose = False
        wbt.work_dir = projdir

        # Temporary outputs
        streams = os.path.basename(strm)
        dir_d8 = os.path.basename(fdir)
        sub_basins_file = os.path.join(projdir, sub_basins)
        stream_id_file = os.path.join(projdir, stream_id)

        # See if the link ID grid has already been created.
        if not os.path.exists(stream_id_file):
            wbt.stream_link_identifier(fdir, strm, stream_id, esri_pntr=True, zero_background=False)

        # Build sub-basins, one for each reach
        wbt.subbasins(dir_d8, streams, sub_basins, esri_pntr=esri_pntr)

        # Create raster object from output
        Sub_Basins = gdal.Open(sub_basins_file, gdalconst.GA_ReadOnly)
        GWBasns = gdal.GetDriverByName('Mem').CreateCopy('', Sub_Basins)
        Sub_Basins = None
        remove_file(sub_basins_file)

    elif in_method == 'Polygon Shapefile or Feature Class':
        print('        Groundwater  polygon shapefile input: {0}'.format(in_Polys))

        # Project the input polygons to the output coordinate system
        Poly_FC = os.path.join(projdir, 'Projected_GW_Basins.shp')
        poly_ds, poly_layer, fieldNames = project_Features(in_Polys, grid_obj.proj)
        out_ds = ogr.GetDriverByName(VectorDriver).CopyDataSource(poly_ds, Poly_FC) # Copy to file on disk
        poly_layer = poly_ds = None

        # Assign a new ID field for basins, numbered 1...n. Add field to store this information if necessary
        print('        Adding auto-incremented basin ID field (1...n)')
        poly_layer = out_ds.GetLayer()
        basinID = "newID"
        if basinID not in fieldNames:
            poly_layer.CreateField(ogr.FieldDefn(basinID, ogr.OFTInteger))
        for num,feature in enumerate(poly_layer):
            feature.SetField(basinID, num+1)      # Add an area field to re-calculate area
            poly_layer.SetFeature(feature)
        poly_layer.ResetReading()
        feature = out_ds = poly_layer = None
        del fieldNames

        # Convert from polygon features to raster object.
        GWBasns = FeatToRaster(Poly_FC, strm, basinID, gdal.GDT_Int32, NoData=NoDataVal)
        ogr.GetDriverByName(VectorDriver).DeleteDataSource(Poly_FC)             # Delete temporary shapefile

    print('    Finished building fine-grid groundwater basin grids in {0: 3.2f} seconds'.format(time.time()-tic1))
    return GWBasns

def build_GWBASINS_nc(GW_BUCKS, out_dir, grid_obj):
    '''
    5/17/2017: This function will build the 2-dimensional groundwater bucket
    grid in netCDF format.
    '''

    tic1 = time.time()
    out_file = os.path.join(out_dir, GWGRID_nc)
    varList2D = [['BASIN', 'i4', 'Basin ID corresponding to GWBUCKPARM table values']]

    # Build output 2D GWBASINS netCDF file
    rootgrp = netCDF4.Dataset(out_file, 'w', format=outNCType)
    rootgrp, grid_mapping = create_CF_NetCDF(grid_obj, rootgrp, out_dir, addLatLon=False, notes='', addVars=varList2D)
    del grid_mapping

    # Array size check
    GWBasns_arr = BandReadAsArray(GW_BUCKS.GetRasterBand(1))
    print('    NC dimensions: {0}, {1}'.format(len(rootgrp.dimensions['y']), len(rootgrp.dimensions['x'])))
    print('    GWBUCKS array dimensions: {0}, {1}'.format(GWBasns_arr.shape[0], GWBasns_arr.shape[1]))

    # Add groundwater buckets to the file (flip UP-DOWN?)
    ncvar = rootgrp.variables[varList2D[0][0]]                                  # 'BASIN'
    ncvar[:] = GWBasns_arr[:]                                                   # Populate variable with groundwater bucket array
    rootgrp.close()
    print('    Process: {0} completed without error'.format(out_file))
    print('    Finished building groundwater grid file in {0: 3.2f} seconds'.format(time.time()-tic1))
    del GW_BUCKS, GWBasns_arr, varList2D, ncvar, rootgrp
    return

def build_GWBUCKPARM(out_dir, cat_areas, cat_comids):
    '''
    5/17/2017: This function will build the groundwater bucket parameter table.
               Currently, only netCDF output format is available.
    '''
    tic1 = time.time()
    Community = True                                                            # Switch to provide Community WRF-Hydro GWBUCKPARM outputs

    # Produce output in NetCDF format (binary and much faster to produce)
    out_file = os.path.join(out_dir, GW_nc)                                     # Groundwater bucket parameter table path and filename
    rootgrp = netCDF4.Dataset(out_file, 'w', format=outNCType)

    # Create dimensions and set other attribute information
    dim1 = 'feature_id'
    dim = rootgrp.createDimension(dim1, len(cat_comids))

    # Create fixed-length variables
    Basins = rootgrp.createVariable('Basin', 'i4', (dim1))                  # Variable (32-bit signed integer)
    coeffs = rootgrp.createVariable('Coeff', 'f4', (dim1))                  # Variable (32-bit floating point)
    Expons = rootgrp.createVariable('Expon', 'f4', (dim1))                  # Variable (32-bit floating point)
    Zmaxs = rootgrp.createVariable('Zmax', 'f4', (dim1))                    # Variable (32-bit floating point)
    Zinits = rootgrp.createVariable('Zinit', 'f4', (dim1))                  # Variable (32-bit floating point)
    Area_sqkms = rootgrp.createVariable('Area_sqkm', 'f4', (dim1))          # Variable (32-bit floating point)
    ComIDs = rootgrp.createVariable('ComID', 'i4', (dim1))                  # Variable (32-bit signed integer)
    if addLoss:
        LossF = rootgrp.createVariable('Loss', 'f4', (dim1))                  # Variable (32-bit signed integer)

    # Set variable descriptions
    Basins.long_name = 'Basin monotonic ID (1...n)'
    coeffs.long_name = 'Coefficient'
    Expons.long_name = 'Exponent'
    Zmaxs.long_name = 'Bucket height'
    Zinits.long_name = 'Initial height of water in bucket'
    Area_sqkms.long_name = 'Basin area in square kilometers'
    if Community:
        ComIDs.long_name = 'Catchment Gridcode'
    else:
        ComIDs.long_name = 'NHDCatchment FEATUREID (NHDFlowline ComID)'     # For NWM
    Zmaxs.units = 'mm'
    Zinits.units = 'mm'
    Area_sqkms.units = 'km2'
    if addLoss:
        LossF.units = '-'
        LossF.long_name = "Fraction of bucket output lost"

    # Fill in global attributes
    rootgrp.featureType = 'point'                                           # For compliance
    rootgrp.history = 'Created {0}'.format(time.ctime())

    # Fill in variables
    Basins[:] = cat_comids                                                  #Basins[:] = numpy.arange(1,cat_comids.shape[0]+1)
    coeffs[:] = coeff
    Expons[:] = expon
    Zmaxs[:] = zmax
    Zinits[:] = zinit
    Area_sqkms[:] = numpy.array(cat_areas)
    ComIDs[:] = numpy.array(cat_comids)
    if addLoss:
        LossF[:] = Loss

    # Close file
    rootgrp.close()
    print('    Created output bucket parameter table (.nc): {0}.'.format(out_file))
    del rootgrp, dim1, dim, Basins, coeffs, Expons, Zmaxs, Zinits, Area_sqkms, ComIDs, out_file

    # Print statements and return
    print('  Finished building groundwater bucket parameter table in {0: 3.2f} seconds.'.format(time.time()-tic1))
    del tic1, Community, cat_areas, cat_comids
    return

def build_GW_buckets(out_dir, GWBasns, grid_obj, Grid=True, saveRaster=False):
    '''
    5/17/2017: This function will build the groundwater bucket grids and parameter
               tables.

    1) A set of basins must be provided. This is a grid of watershed pixels, in
       which each value on the grid corresponds to a basin.

    Build Options:
        1) Build option 1 will biuld the groundwater buckets from ...

    NOTES:
       * Groundwater buckets are currently resolved on the LSM (coarse/GEOGRID)
         grid. In the future this may change.
       * The ID values for groundwater buckets must be numbered 1...n, and will
         not directly reflect the COMID values of individual basins or pour points.
         Much like the LAKEPARM and LAKEGRID values, a mapping must be made between
         input basins/lakes and the outputs using the shapefiles/parameter tables
         output by these tools.
    '''

    # (out_dir, GWBasns, grid_obj, Grid, saveRaster) = (projdir, GWBasns, coarse_grid, True, False)

    tic1 = time.time()
    print('    Beginning to build coarse-grid groundwater basins and parameters')

    # Read basin information from the array
    GWBasns_arr = BandReadAsArray(GWBasns.GetRasterBand(1))                     # Read input raster into array
    ndv = GWBasns.GetRasterBand(1).GetNoDataValue()                             # Obtain nodata value
    UniqueVals = numpy.unique(GWBasns_arr[GWBasns_arr!=ndv])                    # Array to store the basin ID values in the fine-grid groundwater basins
    UniqueVals = UniqueVals[UniqueVals>=0]                                      # Remove NoData, removes potential noData values (-2147483647, -9999)
    print('        Found {0} basins in the watershed grid'.format(UniqueVals.shape[0]))
    del UniqueVals, GWBasns_arr

    # Resample fine-grid groundwater basins to coarse grid
    GW_BUCKS = grid_obj.project_to_model_grid(GWBasns, resampling=gdal.GRA_NearestNeighbour)
    del GWBasns

    # Re-assign basin IDs to 1...n because sometimes the basins get lost when converting to coarse grid
    band = GW_BUCKS.GetRasterBand(1)
    GWBasns_arr2 = BandReadAsArray(band)                                        # Create array from raster
    ndv = band.GetNoDataValue()                                                 # Obtain nodata value
    GWBasns_arr2[GWBasns_arr2==ndv] = NoDataVal                                 # Ensure all non-basin areas are NoData
    UniqueVals2 = numpy.unique(GWBasns_arr2[:])                                 # Get the unique values, including nodata
    GW_BUCKS = band = ndv = None                                                # Destroy the resampled-to-coarse-grid groundwater basin raster
    print('        Found {0} basins (potentially including nodata values) in the file after resampling to the coarse grid.'.format(UniqueVals2.shape[0]))

    '''Because we resampled to the coarse grid, we lost some basins. Thus, we need to
    re-assign basin ID values to conform to the required 1...n groundwater basin
    ID assignment scheme.'''
    # Fast replace loop from https://stackoverflow.com/questions/3403973/fast-replacement-of-values-in-a-numpy-array
    # This method ensures that any nodata values are issued a 0 index in sort_idx
    sort_idx = numpy.argsort(UniqueVals2)                                       # Index of each unique value, 0-based
    # Added .astype on each input because dtype was changing at the argsort step from int32 to int64 on linux!
    idx = numpy.searchsorted(UniqueVals2, GWBasns_arr2, sorter=sort_idx).astype(UniqueVals2.dtype) # 2D array of index values against GWBasns_arr2
    del GWBasns_arr2, sort_idx                                                  # Free up memory

    # This method requires the nodata value to be issued a 0 index in sort_idx and idx
    to_values = numpy.arange(UniqueVals2.size).astype(UniqueVals2.dtype)        # 0..n values to be substituted, 0 in place of NoDataVal
    GWBasns_arr3 = to_values[idx]                                               # Same as to_values[sort_idx][idx]
    if numpy.where(UniqueVals2==NoDataVal)[0].shape[0] > 0:
        new_ndv = int(to_values[numpy.where(UniqueVals2==NoDataVal)[0]][0])     # Obtain the newly-assigned nodatavalue
    else:
        new_ndv = NoDataVal
        GWBasns_arr3+=1                                                         # Add one so that the basin IDs will be 1...n rather than 0...n when there are no nodata values in the grid
    GWBasns_arr3[GWBasns_arr3==new_ndv] = NoDataVal
    del UniqueVals2

    # Build rasters and arrays to create the NC or ASCII outputs
    GW_BUCKS = grid_obj.numpy_to_Raster(GWBasns_arr3)
    del GWBasns_arr3, idx, to_values, new_ndv

    # If requested, create 2D gridded bucket parameter table
    if Grid:
        build_GWBASINS_nc(GW_BUCKS, out_dir, grid_obj)

    # If requested, save a raster object of the 2D groundwater grid
    if saveRaster:
        out_ds = gdal.GetDriverByName(RasterDriver).CreateCopy(os.path.join(out_dir, basinRaster), GW_BUCKS)
        out_ds = None

    # Alternate method to obtain IDs - read directly from raster attribute table
    print('        Calculating size and ID parameters for basin polygons.')
    GW_BUCKS_arr = BandReadAsArray(GW_BUCKS.GetRasterBand(1))
    GW_BUCKS = None
    uniques = numpy.unique(GW_BUCKS_arr, return_counts=True)
    cat_comids = uniques[0].tolist()
    pixel_counts = uniques[1].tolist()
    del uniques, GW_BUCKS, GW_BUCKS_arr
    cat_areas = [float((item*(grid_obj.DX**2))/1000000) for item in pixel_counts]  # Assumes DX is in units of meters
    del pixel_counts

    # Build the groundwater bucket parameter table in netCDF format
    build_GWBUCKPARM(out_dir, cat_areas, cat_comids)

    # Clean up and return
    del cat_comids, cat_areas
    print('    Finished building groundwater parameter files in {0: 3.2f} seconds'.format(time.time()-tic1))
    return

def force_edges_off_grid(fd_arr, ignore_vals=[]):
    '''
    10/23/2020:
        This function is intended to resolve an incompatibility between Whitebox'
        fill_depressions_planchon_and_darboux method and WRF-Hydro. This method
        is known to produce flow direction values of 0 on grid edges where a flow
        direction cannot be determined. Unlike ArcGIS Fill tool, which will optionally
        allow users to set edge cells to flow off the edge of the grid, Whitebox
        tools do not allow this. Thus, those 0-value flow directions must be identified
        and resolved.

        The resolution is to find which edge the cell is on and force the flow direction
        to flow off of the edge of the grid.

        This function uses the Esri flow-direction scheme. All cells on the left edge
        will flow to the right off the right edge, then any remaining edge cells on the
        top edge will flow upward off the top edge, then any cells on the right edge
        will flow right off the right edge, and finally any cells along the bottom
        will flow downward off the bottom edge.
    '''
    tic1 = time.time()

    silent = True       # Supress messages about cells flowing off of the edge.

    # Get array shape
    arr_shape = fd_arr.shape
    rows = arr_shape[0]
    cols = arr_shape[1]

    # Find indices of cells that have 0 flow direction
    zero_index_arr = numpy.argwhere(fd_arr==0)                                  # Get indices of all zero values in input grid
    fd_arr_out = fd_arr.copy()                                                  # Make a copy
    counter = 0
    for entry in zero_index_arr:
        # The indices from numpy.argwhere are y,x
        i = entry[1]                    # x-direction index
        j = entry[0]                    # y-direction index
        if fd_arr_out[j,i] == 0:
            if i == 0:
                # Left edge
                fd_arr_out[j,i] = 16
            elif j == 0:
                # Top edge
                fd_arr_out[j,i] = 64
            elif i == cols-1:
                # Right edge
                fd_arr_out[j,i] = 1
            elif j == rows-1:
                # Bottom edge
                fd_arr_out[j,i] = 4
            else:
                if not silent:
                    print('        Not able to determine which edge 0-value in flow direction grid at index [{0},{1}] (i,j) flows to.'.format(i,j))
            counter+=1
    if fd_arr_out[fd_arr_out==0].shape[0] == 0:
        print('    Coerced {0} 0-value flow direction cells to flow off of the grid.'.format(counter))
    else:
        print('    Could not corece all 0-value flow direction cells to flow off of the grid.')
    return fd_arr_out

def WB_functions(rootgrp, indem, projdir, threshold, ovroughrtfac_val, retdeprtfac_val, lksatfac_val, sink=False, startPts=None, chmask=None):
    """
    This function is intended to produce the hydroglocial DEM corrections and derivitive
    products using the Whitebox tools suite.

    sink: Flag to identify sinks in the input DEM
    """

    tic1 = time.time()
    print('    Terrain processing step initiated...')

    # Whitebox options for running Whitebox in a full workflow
    wbt = WhiteboxTools()
    print('        Using {0}'.format(wbt.version().split('\n')[0]))
    wbt.work_dir = projdir                              # Set working directory
    wbt.verbose = False                                 # Verbose output. [True, False]
    esri_pntr = True                                    # Use the Esri flow direction classification scheme
    fac_type = 'cells'                                  # Output type; one of 'cells', 'sca' (default), and 'ca'

    # Workflow options
    Full_Workflow = False                               # Use the Flow Accumulation Full Workflow tool (fewer options)
    default_Method = False                               # Fill Depressions (Planchon and Darboux), no z-Limit functionality
    fill_deps = True                                    # Option to Fill Depressions with z_limit
    breach_deps = False                                 # Option to Breach Depressions
    breach_deps_LC = False                              # Option to use Breach Depressions (Least Cost)
    zero_background_stream_order = True                 # 2021/09/24 Adding option for specifying zero-background as output of stream order tools
    fill_depth_raster = False                           # 2022/10/05 - For diagnostics, we can opt to create a grid of fill depths.

    # Temporary output files
    flow_acc = "flow_acc.tif"
    fill_pits = "fill_pits.tif"

    # Workflow for diagnosing sink extent and sink depth on input DEM.
    if sink:
        sink_output = "sinks.tif"
        sink_depth = "sink_depth.tif"
        print('      Outputting layer of sink locations: {0}'.format(sink_output))
        wbt.sink(
            indem,
            sink_output,
            zero_background=False)
        print('      Outputting layer of sink depths: {0}'.format(sink_depth))
        wbt.depth_in_sink(
            indem,
            sink_depth,
            zero_background=False)

    # Determine which terrain processing workflow to follow from Whitebox Tools tools.
    if Full_Workflow:
        print('        Algorithm: Whitebox Flow Accumulation Full Workflow.')
        # Perform Fill, Flow Direction, and Flow Accumulation in one step
        wbt.flow_accumulation_full_workflow(
            indem,
            fill_pits,
            dir_d8,
            flow_acc,
            out_type=fac_type,
            esri_pntr=esri_pntr)
    else:
        # Runs each whitebox tool separately
        if fill_deps:
            print('        Depression Filling algorithm: Whitebox Fill Depressions.')

            # Fill Depressions options
            fix_flats = True                                   # Optional flag indicating whether flat areas should have a small gradient applied. [True, False]
            max_depth = z_limit                                # Optional maximum breach depth (default is Inf) [None]
            fill_pits_bool = False                              # Option to fill single cell pits

            # Optional elevation increment applied to flat areas in Breach Depressions and Fill Depressions tools
            # Appropriate values from 0.01 - 0.00001, None
            flat_increment = None                      # 0.0001
            #flat_increment = 0.0001                      # None

            if fill_pits_bool:
                wbt.fill_single_cell_pits(indem, fill_pits)
                indem = fill_pits

            wbt.fill_depressions(
                indem,
                fill_depressions,
                fix_flats=fix_flats,
                flat_increment=flat_increment,
                max_depth=max_depth)
            #remove_file(os.path.join(projdir, fill_pits))                  # Delete temporary file

        elif breach_deps:
            '''
            Note that even if a maximum breach length or maximum breach depth is
            specified, the Whitebox-Tools Breach Depressions tool will currently
            fill all remaining depressions after the breach depth or length are
            exceeded such that the output DEM will be completely depressionless.
            This is not always desirable.
            '''
            print('        Depression Breaching algorithm: Whitebox Breach Depressions (Lindsay, 2016).')

            # Breach Depression options
            max_length = x_limit                               # Optional maximum breach channel length (in grid cells; default is Inf) [None]
            max_depth = z_limit                                # Optional maximum breach depth (default is Inf) [None]
            fill_pits_bool = True                              # In Breach Depressions tool, option to fill single cell pits

            # Optional elevation increment applied to flat areas in Breach Depressions and Fill Depressions tools
            # Appropriate values from 0.01 - 0.00001, None
            flat_increment = None                      # 0.0001

            #wbt.fill_single_cell_pits(indem, fill_pits)
            #wbt.breach_single_cell_pits(indem, fill_pits)

            wbt.breach_depressions(
                indem,
                fill_depressions,
                max_depth=max_depth,
                max_length=max_length,
                fill_pits=fill_pits_bool,
                flat_increment=flat_increment)

            # The below code is slower, either chained or just using the fill_pits option in the breach depressions algorithm.
            #remove_file(os.path.join(projdir, fill_pits))                       # Delete temporary file

        elif breach_deps_LC:
            print('        Depression Breaching algorithm: Whitebox Breach Depressions Least Cost (Lindsay and Dhun, 2015).')

            # Breach Depressions Least Cost Options
            fill_remaining = True                           # Optional flag indicating whether to fill any remaining unbreached depressions
            dist = 100                                       # Undocumented parameter. Serach radius? Breach distance? 1000. Large numbers have a huge effect on processing time
            min_dist_bool = True                            # Optional flag indicating whether to minimize breach distances
            max_cost_val = None                             # Optional maximum breach cost (default is Inf) (None)

            # Optional elevation increment applied to flat areas in Breach Depressions and Fill Depressions tools
            # Appropriate values from 0.01 - 0.00001, None
            flat_increment = None                      # 0.0001

            wbt.breach_depressions_least_cost(
                indem,
                fill_depressions,
                dist,
                max_cost=max_cost_val,
                min_dist=min_dist_bool,
                flat_increment=flat_increment,
                fill=fill_remaining)

        elif default_Method:
            # Fill Depressions (Planchon and Darboux). This method mimicks the
            # Esri Spatial Analyst "Fill" tool with no fill limit.
            # Planchon, O. and Darboux, F., 2002. A fast, simple and versatile algorithm to fill the depressions of digital elevation models. Catena, 46(2-3), pp.159-176.

            print('        Depression Filling algorithm: Planchon and Darboux (2002).')
            # Fill Depressions (Planchon and Darboux) options
            fix_flats = False                               # Optional flag indicating whether flat areas should have a small gradient applied. [True, False]
            flat_increment = None                           # 0.0001

            wbt.fill_depressions_planchon_and_darboux(indem,
                fill_depressions,
                fix_flats=fix_flats,
                flat_increment=flat_increment)

        if fill_depth_raster:
            # Create a fill depth raster for diagnostic purposes
            fill_depth_raster = 'Fill_Depth.tif'
            wbt.subtract(fill_depressions, indem, fill_depth_raster)

        # This is the variable name for the output filled DEM
        fill_pits = fill_depressions

        # Build the flow direction grid (used to populate FLOWDIRECTION in Fulldom_hires.nc)
        wbt.d8_pointer(fill_pits, dir_d8, esri_pntr=esri_pntr)

        # Build the flow accumulation grid (used to populate FLOWACC in Fulldom_hires.nc)
        wbt.d8_flow_accumulation(fill_pits, flow_acc)

    # Process: Write hydrologically pre-processed DEM to Fulldom_hires.nc
    fill_pits_file = os.path.join(projdir, fill_pits)
    fill_arr, ndv = return_raster_array(fill_pits_file)
    fill_arr[fill_arr==ndv] = NoDataVal                                         # Replace raster NoData with WRF-Hydro NoData value
    rootgrp.variables['TOPOGRAPHY'][:] = fill_arr
    print('        Process: TOPOGRAPHY written to output netCDF.')
    del fill_arr, ndv

    # Process: Flow Direction
    dir_d8_file = os.path.join(projdir, dir_d8)
    fdir_arr, ndv = return_raster_array(dir_d8_file)
    #fdir_arr[fdir_arr==ndv] = 255                                               # Replace raster NoData with specific value
    fdir_arr_out = force_edges_off_grid(fdir_arr)                               # Force 0-value cells to flow off edge of grid.
    fdir_arr_out[fdir_arr_out==ndv] = 255                                       # Replace raster NoData with specific value
    rootgrp.variables['FLOWDIRECTION'][:] = fdir_arr_out
    print('        Process: FLOWDIRECTION written to output netCDF.')
    del fdir_arr, fdir_arr_out, ndv

    # Process: Flow Accumulation (intermediate
    flow_acc_file = os.path.join(projdir, flow_acc)
    flac_arr, ndv = return_raster_array(flow_acc_file)
    flac_arr[flac_arr==ndv] = 0                                                 # Set NoData values to 0 on Flow Accumulation grid
    rootgrp.variables['FLOWACC'][:] = flac_arr
    print('        Process: FLOWACC written to output netCDF.')
    del flac_arr, ndv

    # Create stream channel raster
    if not startPts:
        # Create stream channel raster according to threshold
        print('        Flow accumulation will be thresholded to build channel pixels.')
        wbt.extract_streams(flow_acc, streams, threshold, zero_background=zero_background_stream_order)
    else:
        # Added 8/14/2020 to use a vector of points to seed the channelgrid

        # Project input points and clip to domain if necessary
        geom = boundarySHP(dir_d8_file)                                         # Get domain extent for cliping geometry
        pt_ds, pt_layer, fieldNames = project_Features(startPts,
                                                            geom.GetSpatialReference(),
                                                            clipGeom=geom,
                                                            geomType=ogr.wkbPoint)

        # Save to disk
        temp_pts = os.path.join(projdir, start_pts_temp)
        out_ds = ogr.GetDriverByName(VectorDriver).CopyDataSource(pt_ds, temp_pts)
        pt_layer = pt_ds = fieldNames = out_ds = geom = None

        print('        Flow accumulation will be weighted using input channel initiation points.')
        wbt.trace_downslope_flowpaths(temp_pts, dir_d8, streams, esri_pntr=esri_pntr, zero_background=zero_background_stream_order)

        driver = ogr.Open(temp_pts).GetDriver()
        driver.DeleteDataSource(temp_pts)                                       # Delete input file
        del temp_pts, geom

    # Define the location for the stream raster grid
    streams_file = os.path.join(projdir, streams)

    # Added 9/2/2022 - Option to mask the Channelgrid layer to a mask raster (1 or NoData on the routing grid)
    if not chmask:
        print('        No masking of CHANNELGRID will be performed.')
    if chmask is not None:
        print('        Masking CHANNELGRID to user-provided mask raster.')

        # Open streams file in update mode
        ds = gdal.Open(streams_file, gdalconst.GA_Update)
        band = ds.GetRasterBand(1)
        strm_arr = band.ReadAsArray()
        ndv = band.GetNoDataValue()                                             # Obtain nodata value

        # Open channelgrid mask file in read-only mode
        ds_chgrid = gdal.Open(chmask, gdalconst.GA_ReadOnly)
        band2 = ds_chgrid.GetRasterBand(1)
        chmask_arr = band2.ReadAsArray()
        chmask_ndv = band2.GetNoDataValue()                                             # Obtain nodata value
        assert chmask_arr.shape == strm_arr.shape

        # Set all values outside the mask to NoData
        strm_arr[chmask_arr==chmask_ndv] = ndv

        # Write changes back to file
        BandWriteArray(ds.GetRasterBand(1), strm_arr)
        stats = ds.GetRasterBand(1).GetStatistics(0,1)                          # Calculate statistics
        ds = ds_chgrid = band = band1 = stats = None
        del chmask_arr, chmask_ndv, strm_arr, ndv, ds, ds_chgrid, band, band1, stats

    # Below this point are modifications to set CHANNELGRID into Fulldom_hires.nc format
    strm_arr, ndv = return_raster_array(streams_file)
    strm_arr[strm_arr==ndv] = NoDataVal
    if zero_background_stream_order:
        strm_arr[strm_arr == 0] = NoDataVal

        # Must set NoData in this file because it is used later by reach-based routing routine.
        ds = gdal.Open(streams_file, 1)
        ds.GetRasterBand(1).SetNoDataValue(0)                    # Set noData
        ds = None

    # Modify the value of channels from 1 to 0 or from anything other than NoData to 0.
    if not startPts:
        strm_arr[strm_arr==1] = 0
    if startPts is not None:
        if zero_background_stream_order:
            # Added 03/25/2023 to avoid entire grid getting set to 0
            strm_arr[strm_arr>0] = 0
        else:
            # Added 10/6/2020 to reclassify results of trace_downslope_flowpaths
            strm_arr[strm_arr!=ndv] = 0

    # Write Channelgrid layer to Fulldom_hires.nc
    rootgrp.variables['CHANNELGRID'][:] = strm_arr
    print('        Process: CHANNELGRID written to output netCDF.')
    del strm_arr, ndv, flow_acc

    # Process: Stream Order
    wbt.strahler_stream_order(dir_d8, streams, strahler, esri_pntr=esri_pntr, zero_background=zero_background_stream_order)
    strahler_file = os.path.join(projdir, strahler)
    strahler_arr, ndv = return_raster_array(strahler_file)
    if zero_background_stream_order:
        strahler_arr[strahler_arr==0] = NoDataVal

        # Must set NoData in this file because it is used later by reach-based routing routine.
        ds = gdal.Open(strahler_file, 1)
        ds.GetRasterBand(1).SetNoDataValue(0)                    # Set noData
        ds = None

    # -9999 does not fit in the 8-bit types, so it gets put in as -15 by netCDF4 for some reason
    strahler_arr[strahler_arr==ndv] = NoDataVal
    rootgrp.variables['STREAMORDER'][:] = strahler_arr
    print('        Process: STREAMORDER written to output netCDF.')
    del strahler_arr, ndv

    # Create initial constant raster of value retdeprtfac_val
    rootgrp.variables['RETDEPRTFAC'][:] = float(retdeprtfac_val)
    print('        Process: RETDEPRTFAC written to output netCDF.')

    # Create initial constant raster of ovroughrtfac_val
    rootgrp.variables['OVROUGHRTFAC'][:] = float(ovroughrtfac_val)
    print('        Process: OVROUGHRTFAC written to output netCDF.')

    # Create initial constant raster of LKSATFAC
    rootgrp.variables['LKSATFAC'][:] = float(lksatfac_val)
    print('        Process: LKSATFAC written to output netCDF.')

    # We will assume that no forecast points, basin masks, or lakes are provided
    rootgrp.variables['frxst_pts'][:] = NoDataVal
    rootgrp.variables['basn_msk'][:] = NoDataVal
    rootgrp.variables['LAKEGRID'][:] = NoDataVal

    print('    Terrain processing step completed without error in {0: 3.2f} seconds.'.format(time.time()-tic1))
    return rootgrp, dir_d8_file, flow_acc_file, streams_file, fill_pits_file, strahler_file

def CSV_to_SHP(in_csv, DriverName='MEMORY', xVar='LON', yVar='LAT', idVar='FID', toProj=None):
    tic1 = time.time()

    drv = ogr.GetDriverByName(DriverName)
    if drv is None:
        print('      {0} driver not available.'.format(DriverName))
        raise SystemExit
    else:
        data_source = drv.CreateDataSource('')

    # Read the input CSV file
    csv_arr = numpy.genfromtxt(in_csv, delimiter=',', names=True)

    # create the spatial reference for the input point CSV file, WGS84
    srs = osr.SpatialReference()
    srs.ImportFromProj4(wgs84_proj4)

    # Handle coordinate transformation
    if toProj is not None:
        # Create the spatial reference for the output
        out_srs = osr.SpatialReference()
        out_srs.ImportFromWkt(toProj)
        transform = osr.CoordinateTransformation(srs, out_srs)
    else:
        out_srs = srs.Clone()

    # Added 11/19/2020 to allow for GDAL 3.0 changes to the order of coordinates in transform
    if int(osgeo.__version__[0]) >= 3:
        # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
        srs.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        out_srs.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    # Add the fields we're interested in
    layer = data_source.CreateLayer('frxst_FC', out_srs, ogr.wkbPoint)
    layer.CreateField(ogr.FieldDefn(idVar, ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn(yVar, ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn(xVar, ogr.OFTReal))
    featureDefn = layer.GetLayerDefn()

    # Process the text file and add the attributes and features to the shapefile
    for row in csv_arr:
        x = row[xVar]
        y = row[yVar]

        # Set the attributes using the values from the delimited text file
        feature = ogr.Feature(featureDefn)                               # create the feature
        feature.SetField(idVar, row[idVar])
        feature.SetField(yVar, y)
        feature.SetField(xVar, x)

        #create point geometry
        if toProj is not None:
            x,y,z = transform.TransformPoint(x,y)
            pass
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x,y)

        # Create the feature and set values
        feature.SetGeometry(point)
        layer.CreateFeature(feature)
        feature = point = None
        del x, y
    layer = srs = drv = None                                      # Saveand close the data source
    return data_source

def forecast_points(in_csv, rootgrp, bsn_msk, projdir, DX, WKT, fdir, fac, strm):
    # (in_csv, rootgrp, bsn_msk, projdir, template_raster) = (in_csv, rootgrp2, basin_mask, projdir, mosprj)

    tic1 = time.time()

    # Setup whitebox tool object and options
    wbt = WhiteboxTools()
    wbt.verbose = False
    wbt.work_dir = projdir
    esri_pntr = True

    # Setup snap tolerances for snapping forecast points to channel pixels
    snap_dist1 = int(DX)                                                        # This is to put the point on the grid only. One pixel tolerance
    snap_dist2 = int(DX * walker)                                               # This is to search for the within a distance

    # Make feature layer from CSV
    print('    Forecast points provided and basins being delineated.')
    frxst_FC = os.path.join(projdir, 'Temp_Frxst_Pts.shp')
    ds = CSV_to_SHP(in_csv, DriverName='MEMORY', xVar='LON', yVar='LAT', idVar='FID', toProj=WKT)  # In-memory features
    out_ds = ogr.GetDriverByName(VectorDriver).CopyDataSource(ds, frxst_FC)    # Copy to file on disk
    ds = out_ds = None

    # Snap pour points to channel grid within a tolerance
    wbt.jenson_snap_pour_points(frxst_FC, strm, snapPour1, snap_dist2)

    # Convert point shapefile to raster
    frxst_raster = FeatToRaster(frxst_FC, fac, 'FID', gdal.GDT_Int32, NoData=NoDataVal)
    frxst_raster_arr = BandReadAsArray(frxst_raster.GetRasterBand(1))
    frxst_raster = None
    frxst_raster_arr[frxst_raster_arr==0] = NoDataVal                           # Replace 0 with WRF-Hydro NoData value
    rootgrp.variables['frxst_pts'][:] = frxst_raster_arr
    print('    Process: frxst_pts written to output netCDF.')
    del frxst_raster_arr

    # Snap pour points to flow accumulation grid within a tolerance
    wbt.snap_pour_points(frxst_FC, fac, snapPour2, snap_dist2)

    # Delineate above points
    watershed_file = os.path.join(projdir, watersheds)
    wbt.watershed(fdir, snapPour2, watershed_file, esri_pntr=esri_pntr)
    watershed_arr, ndv = return_raster_array(watershed_file)
    watershed_arr[watershed_arr==ndv] = NoDataVal                               # Replace raster NoData with WRF-Hydro NoData value
    rootgrp.variables['basn_msk'][:] = watershed_arr
    print('    Process: basn_msk written to output netCDF.')
    remove_file(watershed_file)                                                 # Delete fac from disk
    del ndv, watershed_file

    # Delete temporary point shapefiles
    ogr.GetDriverByName(VectorDriver).DeleteDataSource(frxst_FC)
    ogr.GetDriverByName(VectorDriver).DeleteDataSource(os.path.join(projdir, snapPour1))
    ogr.GetDriverByName(VectorDriver).DeleteDataSource(os.path.join(projdir, snapPour2))

    # Set mask for future raster output
    if bsn_msk:
        print('    Channelgrid will be masked to basins.')
        channelgrid_arr = rootgrp.variables['CHANNELGRID'][:]

        # Converts channelgrid values inside basins to 0, outside to -1
        channelgrid_arr[numpy.logical_and(watershed_arr>=0, channelgrid_arr!=NoDataVal)] = 0
        channelgrid_arr[numpy.logical_and(watershed_arr<0, channelgrid_arr!=NoDataVal)] = -1
        rootgrp.variables['CHANNELGRID'][:] = channelgrid_arr
        print('    Process: CHANNELGRID written to output netCDF.')

        # Eliminate inactive channel cells by setting all areas outside of basin_msk to NoData
        ds = gdal.Open(strm, 1)                                                 # Read for writing
        band = ds.GetRasterBand(1)
        strm_arr = band.ReadAsArray()
        ndv = band.GetNoDataValue()                                             # Obtain nodata value
        strm_arr[channelgrid_arr!=0] = ndv                                      # Any inactivate channels will become Nodata in stream array
        band.WriteArray(strm_arr)                                               # Write the array to the disk
        stats = band.GetStatistics(0,1)                                         # Calculate statistics
        ds = band = ndv = stats = None
        del strm_arr, channelgrid_arr
    else:
        print('    Channelgrid will not be masked to basins.')
    del watershed_arr

    print('    Built forecast point outputs in {0: 3.2f} seconds.'.format(time.time()-tic1))
    return rootgrp

def build_RouteLink(RoutingNC, order, From_To, NodeElev, NodesLL, NodesXY, Lengths, StrOrder, Slopes, gageDict=None):
    '''
    8/10/2017: This function is designed to build the routiing parameter netCDF file.
                Ideally, this will be the only place the produces the file, and
                all functions wishing to write the file will reference this function.
    '''
    tic1 = time.time()

    # To create a netCDF parameter file
    rootgrp = netCDF4.Dataset(RoutingNC, 'w', format=outNCType)

    # Create dimensions and set other attribute information
    #dim1 = 'linkDim'
    dim1 = 'feature_id'
    dim2 = 'IDLength'
    dim = rootgrp.createDimension(dim1, len(order))
    gage_id = rootgrp.createDimension(dim2, 15)

    # Create fixed-length variables
    ids = rootgrp.createVariable('link', 'i4', (dim1))                          # Variable (32-bit signed integer)
    froms = rootgrp.createVariable('from','i4',(dim1))                          # Variable (32-bit signed integer)
    tos = rootgrp.createVariable('to','i4',(dim1))                              # Variable (32-bit signed integer)
    slons = rootgrp.createVariable('lon', 'f4', (dim1))                         # Variable (32-bit floating point)
    slats = rootgrp.createVariable('lat', 'f4', (dim1))                         # Variable (32-bit floating point)
    selevs = rootgrp.createVariable('alt', 'f4', (dim1))                        # Variable (32-bit floating point)
    orders = rootgrp.createVariable('order','i4',(dim1))                        # Variable (32-bit signed integer)
    Qis = rootgrp.createVariable('Qi', 'f4', (dim1))                            # Variable (32-bit floating point)
    MusKs = rootgrp.createVariable('MusK','f4',(dim1))                          # Variable (32-bit floating point)
    MusXs = rootgrp.createVariable('MusX', 'f4', (dim1))                        # Variable (32-bit floating point)
    Lengthsnc = rootgrp.createVariable('Length', 'f4', (dim1))                  # Variable (32-bit floating point)
    ns = rootgrp.createVariable('n', 'f4', (dim1))                              # Variable (32-bit floating point)
    Sos = rootgrp.createVariable('So', 'f4', (dim1))                            # Variable (32-bit floating point)
    ChSlps = rootgrp.createVariable('ChSlp', 'f4', (dim1))                      # Variable (32-bit floating point)
    BtmWdths = rootgrp.createVariable('BtmWdth','f4',(dim1))                    # Variable (32-bit floating point)
    Times = rootgrp.createVariable('time', 'f4')                                # Scalar Variable (32-bit floating point)
    geo_x = rootgrp.createVariable('x', 'f4', (dim1))                           # Variable (32-bit floating point)
    geo_y = rootgrp.createVariable('y', 'f4', (dim1))                           # Variable (32-bit floating point)
    Kcs = rootgrp.createVariable('Kchan', 'i2', (dim1))                         # Variable (16-bit signed integer)
    Gages = rootgrp.createVariable('gages', 'S1', (dim1, dim2))                 # Variable (string type character) Added 07/27/2015 - 15 character strings
    LakeDis = rootgrp.createVariable('NHDWaterbodyComID', 'i4', (dim1))         # Variable (32-bit signed integer)

    # Add CF-compliant coordinate system variable
    if pointCF:
        sr = osr.SpatialReference()                                             # Build a spatial reference object
        sr.ImportFromEPSG(pointSR)                                              # Define using EPSG code specified in header
        projEsri = sr.Clone()                                                   # Copy the SRS
        projEsri.MorphToESRI()                                                  # Alter the projection to Esri's representation of a coordinate system
        PE_string = projEsri.ExportToWkt().replace("'", '"')                    # INVESTIGATE - this somehow may provide better compatability with Esri products?
        grid_mapping = crsVar
        rootgrp = add_CRS_var(rootgrp, sr, 0, grid_mapping, 'latitude_longitude', PE_string)

    # Set variable descriptions
    ids.long_name = 'Link ID'
    froms.long_name = 'From Link ID'
    tos.long_name = 'To Link ID'
    slons.long_name = 'longitude of the start node'
    slats.long_name = 'latitude of the start node'
    selevs.long_name = 'Elevation in meters at start node'
    orders.long_name = 'Stream order (Strahler)'
    Qis.long_name = 'Initial flow in link (CMS)'
    MusKs.long_name = 'Muskingum routing time (s)'
    MusXs.long_name = 'Muskingum weighting coefficient'
    Lengthsnc.long_name = 'Stream length (m)'
    ns.long_name = "Manning's roughness"
    Sos.long_name = 'Slope (%; drop/length)'
    ChSlps.long_name = 'Channel side slope (%; drop/length)'
    BtmWdths.long_name = 'Bottom width of channel'
    geo_x.long_name = "x coordinate of projection"
    geo_y.long_name = "y coordinate of projection"
    Kcs.long_name = "channel conductivity"
    LakeDis.long_name = 'ID of the lake element that intersects this flowline'
    Gages.long_name = 'Gage ID'

    # Variable attributes for CF compliance
    slons.units = 'degrees_east'                                                # For compliance
    slats.units = 'degrees_north'                                               # For compliance
    slons.standard_name = 'longitude'                                           # For compliance
    slats.standard_name = 'latitude'                                            # For compliance
    Times.standard_name = 'time'                                                # For compliance
    Times.long_name = 'time of measurement'                                     # For compliance
    Times.units = 'days since 2000-01-01 00:00:00'                              # For compliance
    selevs.standard_name = "height"                                             # For compliance
    selevs.units = "m"                                                          # For compliance
    selevs.positive = "up"                                                      # For compliance
    selevs.axis = "Z"                                                           # For compliance
    ids.cf_role = "timeseries_id"                                               # For compliance
    geo_x.standard_name = "projection_x_coordinate"
    geo_y.standard_name = "projection_y_coordinate"
    geo_x.units = "m"
    geo_y.units = "m"
    Kcs.units = "mm h-2"
    slons.standard_name = 'longitude'                                           # For compliance with NCO
    slats.standard_name = 'latitude'                                            # For compliance with NCO

    # Apply grid_mapping and coordinates attributes to all variables
    for varname, ncVar in rootgrp.variables.items():
        if dim1 in ncVar.dimensions and varname not in ['alt', 'lat', 'lon', 'x', 'y']:
            ncVar.setncattr('coordinates', 'lat lon')                           # For CF-compliance
            if pointCF:
                ncVar.setncattr('grid_mapping', grid_mapping)                       # For CF-compliance
        del ncVar, varname

    # Fill in global attributes
    rootgrp.featureType = 'timeSeries'                                          # For compliance
    rootgrp.history = 'Created %s' %time.ctime()

    print('        Starting to fill in routing table NC file.')
    ids[:] = numpy.array(order)                                                 # Fill in id field information

    # Change None values to 0.  Could alternatively use numpy.nan
    froms[:] = numpy.zeros(len(order))
    tos[:] = numpy.array([From_To.get(featID, 0) for featID in order])

    # Fill in other variables
    slons[:] = numpy.array([NodesLL[featID][0] for featID in order])
    slats[:] = numpy.array([NodesLL[featID][1] for featID in order])
    geo_x[:] = numpy.array([NodesXY[featID][0] for featID in order])
    geo_y[:] = numpy.array([NodesXY[featID][1] for featID in order])
    selevs[:] = numpy.array([round(NodeElev[featID], 3) for featID in order])   # Round to 3 decimal places
    Lengthsnc[:] = numpy.array([round(Lengths[featID], 1) for featID in order]) # Round to 1 decimal place

    # Modify order and slope arrays
    orders[:] = numpy.array([StrOrder[featID] for featID in order])
    Sos_ = numpy.round(numpy.array([Slopes[featID] for featID in order]), 3)
    Sos[:] = Sos_[:]

    # Set default arrays
    Qis[:] = Qi
    MusKs[:] = MusK
    MusXs[:] = MusX
    Times[:] = 0
    Kcs[:] = Kc
    #LakeDis[:] = NoDataVal                                                      # Fill with default nodata value. Disabled to keep all values uninitialized

    # Apply order-based Mannings N values according to global dictionary "Mannings_Order"
    if ManningsOrd:
        ns[:] = numpy.array([Mannings_Order[item] for item in orders[:]])
    else:
        ns[:] = n

    # Apply order-based Channel Side-slope values according to global dictionary "Mannings_Order"
    if ChSSlpOrd:
        ChSlps[:] = numpy.array([Mannings_ChSSlp[item] for item in orders[:]])
    else:
        ChSlps[:] = ChSlp

    # Apply order-based bottom-width values according to global dictionary "Mannings_Order"
    if BwOrd:
        BtmWdths[:] = numpy.array([Mannings_Bw[item] for item in orders[:]])
    else:
        BtmWdths[:] = BtmWdth

    # Added 10/10/2017 by KMS to include user-supplied gages in reach-based routing files
    if gageDict is not None:
        Gages[:,:] = numpy.asarray([tuple(str(gageDict[arcid]).rjust(15)) if arcid in gageDict else tuple('               ') for arcid in order])
    else:
        Gages[:,:] = numpy.asarray([tuple('               ') for arcid in order])    # asarray converts from tuple to array
    del gageDict

    # Close file
    rootgrp.close()
    print('        Done writing NC file to disk.')
    print('    Routing table created without error.')
    return

def find_line_midpoint(geom):
    '''
    This function will use shapely to find the midpoint along a line, given an
    OGR geometry object.
    '''
    shapelyLine = LineString(wkt.loads(geom.ExportToWkt()))
    midPoint = shapelyLine.interpolate(shapelyLine.length/2)
    outGeom = ogr.CreateGeometryFromWkt(midPoint.wkt)
    del midPoint, shapelyLine
    return outGeom

def Routing_Table(projdir, rootgrp, grid_obj, fdir, strm, Elev, Strahler, gages=False, Lakes=None):
    """If "Create reach-based routing files?" is selected, this function will create
    the Route_Link.nc table and Streams.shp shapefiles in the output directory."""

    # Stackless topological sort algorithm, adapted from: http://stackoverflow.com/questions/15038876/topological-sort-python
    def sort_topologically_stackless(graph):

        '''This function will navigate through the list of segments until all are accounted
        for. The result is a sorted list of which stream segments should be listed
        first. Simply provide a topology dictionary {Fromnode:[ToNode,...]} and a sorted list
        is produced that will provide the order for navigating downstream. This version
        is "stackless", meaning it will not hit the recursion limit of 1000.'''

        levels_by_name = {}
        names_by_level = defaultdict(set)

        def add_level_to_name(name, level):
            levels_by_name[name] = level
            names_by_level[level].add(name)

        def walk_depth_first(name):
            stack = [name]
            while(stack):
                name = stack.pop()
                if name in levels_by_name:
                    continue

                if name not in graph or not graph[name]:
                    level = 0
                    add_level_to_name(name, level)
                    continue

                children = graph[name]

                children_not_calculated = [child for child in children if child not in levels_by_name]
                if children_not_calculated:
                    stack.append(name)
                    stack.extend(children_not_calculated)
                    continue

                level = 1 + max(levels_by_name[lname] for lname in children)
                add_level_to_name(name, level)

        for name in graph:
            walk_depth_first(name)

        list1 = list(takewhile(lambda x: x is not None, (names_by_level.get(i, None) for i in count())))
        list2 = [item for sublist in list1 for item in sublist][::-1]           # Added by KMS 9/2/2015 to reverse sort the list
        list3 = [x for x in list2 if x is not None]                             # Remove None values from list
        return list3

    print('    Routing table will be created...')
    tic1 = time.time()

    # Setup whitebox tool object and options
    wbt = WhiteboxTools()
    wbt.verbose = False
    wbt.work_dir = projdir
    esri_pntr = True
    zero_background = False
    id_field = 'STRM_VAL'                                                       # Whitebox-assigned stream ID field

    # Setup temporary and other outputs
    stream_id_file = os.path.join(projdir, stream_id)
    streams_vector_file = os.path.join(projdir, streams_vector)
    RoutingNC = os.path.join(projdir, RT_nc)

    # Run Whitebox functions for creating link IDs and vectors

    '''
    The stream_link_identifier appears to output an int16 raster, limiting the number
    of individual stream link IDs possible. Further, it will populate negative values
    in the output, providing both positive and negative IDs. Unfortunately, any
    IDs that are assigned negative values in the output will not be resolved as
    stream vectors in the raster_streams_to_vector routine.
    '''
    wbt.stream_link_identifier(fdir, strm, stream_id, esri_pntr=esri_pntr, zero_background=zero_background)
    wbt.raster_streams_to_vector(stream_id, fdir, streams_vector, esri_pntr=esri_pntr)
    print('        Stream to features step complete.')

    # Read the link IDs as an array from the output file
    strm_link_arr, ndv = return_raster_array(stream_id_file)
    if numpy.unique(strm_link_arr).shape[0] > 32768 or strm_link_arr[strm_link_arr<0].shape[0] > 0:
        print('        Warning: Number of unique IDs exceeds limit of 16-bit unsigned integer type. ' + \
                'Not all reaches may be converted to stream vectors. Check output carefully.')
    strm_link_arr[strm_link_arr==ndv] = NoDataVal                               # Set nodata values to WRF-Hydro nodata value
    strm_link_arr[strm_link_arr<1] = NoDataVal                                  # Remove zeros from background of grid

    # Find any LINKID reach ID values that did not get transferred to the stream vector file.
    # These are typically single-cell channel cells on the edge of the grid.
    ds = ogr.Open(streams_vector_file)
    lyr = ds.GetLayer(0)                                               # Get the 'layer' object from the data source
    vector_reach_IDs = numpy.unique([feature.GetField('STRM_VAL') for feature in lyr]).astype(int)
    print('        Found {0} unique IDs in stream vector layer.'.format(len(vector_reach_IDs)))
    ds = lyr = None

    # Resolve issue where LINKID values are present that do not correspond to a vector ID (10/25/2020)
    grid_reach_IDs = numpy.unique(strm_link_arr[strm_link_arr!=NoDataVal])
    missing_reach_IDs = grid_reach_IDs[~numpy.in1d(grid_reach_IDs, vector_reach_IDs)]
    print('        Eliminating {0} IDs in LINKID grid that could not be resolved in stream vector layer.'.format(missing_reach_IDs.shape[0]))
    print('          {0}'.format(missing_reach_IDs.tolist()))
    channel_arr = rootgrp.variables['CHANNELGRID'][:]
    strorder_arr = rootgrp.variables['STREAMORDER'][:]
    for idVal in missing_reach_IDs:
        arr_mask = strm_link_arr==idVal         # Build a boolean mask for masking all array elements to be changed
        strm_link_arr[arr_mask] = NoDataVal     # Set all linkid values that didn't get resolved in the routelink file to nodata.
        channel_arr[arr_mask] = NoDataVal       # Set all channel values that didn't get resolved in the routelink file to nodata.
        strorder_arr[arr_mask] = NoDataVal      # Set all channel values that didn't get resolved in the routelink file to nodata.
        del arr_mask
    rootgrp.variables['LINKID'][:] = strm_link_arr
    rootgrp.variables['CHANNELGRID'][:] = channel_arr
    rootgrp.variables['STREAMORDER'][:] = strorder_arr
    del channel_arr, strorder_arr, grid_reach_IDs, missing_reach_IDs

    gage_linkID = {}
    if gages:
        print('        Adding forecast points:LINKID association.')
        gage_arr = rootgrp.variables['frxst_pts'][:]
        unique_gages = numpy.unique(gage_arr[gage_arr!=NoDataVal])
        gage_linkID = {gage:strm_link_arr[gage_arr==gage][0] for gage in unique_gages}    # Create blank dictionary so that it exists and can be deleted later
        print('        Found {0} forecast point:LINKID associations.'.format(len(gage_linkID)))
        del unique_gages, gage_arr
    linkID_gage = {val:key for key, val in gage_linkID.items()}                 # Reverse the dictionary
    del strm_link_arr, ndv, gage_linkID

    # Setup coordinate transform for calculating lat/lon from x/y
    wgs84_proj = osr.SpatialReference()
    wgs84_proj.ImportFromProj4(wgs84_proj4)

    # Added 11/19/2020 to allow for GDAL 3.0 changes to the order of coordinates in transform
    if int(osgeo.__version__[0]) >= 3:
        # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
        wgs84_proj.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
    coordTrans = osr.CoordinateTransformation(grid_obj.proj, wgs84_proj)        # Transformation from grid projection to WGS84

    # Initiate dictionaries for storing topology and attribute information
    Lengths = {}                        # Gather the stream feature length
    StrOrder = {}                       # Store stream order for each node
    NodeElev = {}                       # Elevation of the start node
    slope_dic = {}                      # Slope (unitless drop/length)
    NodeLL = {}                         # Dictionary to store geocentric (longitude, latitude) coordinates for every start node
    NodeXY = {}                         # Dictionary to store projectedc (x, y) coordinates for every start node

    # Open shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(streams_vector_file, 1)
    lyr = data_source.GetLayer()
    topology_dic = {}
    coords_dic = {}
    for feature in lyr:
        # Open each feature and get geometry
        flowline_id = int(feature.GetField(id_field))
        geom = feature.GetGeometryRef()
        flowline_length = geom.Length()

        # Get coordinates of first and last point, flow line ID, and flow line length
        first_point = geom.GetPoint(0)
        #mid_point = geom.Centroid()
        mid_point = find_line_midpoint(geom)
        last_point = geom.GetPoint(geom.GetPointCount() - 1)
        first_point_coords = (first_point[0], first_point[1])
        mid_point_coords = (mid_point.GetX(), mid_point.GetY())
        last_point_coords = (last_point[0], last_point[1])

        # Create topology dictionary of 'bottom_point geometry: stream flowline ID'
        try:
            topology_dic[last_point_coords] += [flowline_id]
        except KeyError:
            topology_dic[last_point_coords] = [flowline_id]

        # Create coordinate dictionary of flowline ID: first point, last point, length
        coords_dic[flowline_id] = first_point_coords, last_point_coords, flowline_length, mid_point_coords
        feature = geom = first_point = last_point = None
    lyr.ResetReading()

    # Create to/from dictionary matching bottom point to top point, creating dic of 'from ID: to ID'
    to_from_dic = {}
    for flowline_id, (first_point_coords, last_point_coords, flowline_length, mid_point_coords) in coords_dic.items():
        if first_point_coords in topology_dic:
            #for feature_id in topology_dic[first_point_coords]:
            for feature_id in topology_dic.pop(first_point_coords):
                to_from_dic[feature_id] = flowline_id

    # Add in flowlines with nothing downstream
    for feature_id in coords_dic:
        if feature_id not in to_from_dic:
            to_from_dic[feature_id] = 0
    del topology_dic

    # Get the order of segments according to a simple topological sort
    order = sort_topologically_stackless({key:[val] for key,val in to_from_dic.items()})

    # Open elevation raster
    dem_array = gdal.Open(Elev, 0)
    dem_rb = dem_array.GetRasterBand(1)

    # Open strahler stream order raster
    strahler_array = gdal.Open(Strahler, 0)
    strahler_rb = strahler_array.GetRasterBand(1)

    # Iterate over coordinate dictionary​
    tic2 = time.time()
    for idval, (top_xy, bot_xy, length, mid_xy) in coords_dic.items():

        # Get top/first coordinates values from DEM
        row, col = grid_obj.xy_to_grid_ij(top_xy[0], top_xy[1])
        top_elevation = float(dem_rb.ReadAsArray(col, row, 1, 1))
        strahler_value = int(strahler_rb.ReadAsArray(col, row, 1, 1))

        # Get bottom/last coordinates values from DEM
        row, col = grid_obj.xy_to_grid_ij(bot_xy[0], bot_xy[1])
        bottom_elevation = dem_rb.ReadAsArray(col, row, 1, 1)

        # Fix negative slopes
        drop = top_elevation - bottom_elevation
        slope = drop/length
        if slope < minSo:
            slope = minSo

        # Populate all dictionaries
        slope_dic[idval] = float(slope)
        StrOrder[idval] = strahler_value
        NodeElev[idval] = top_elevation
        Lengths[idval] = length
        NodeXY[idval] = (top_xy[0], top_xy[1])

        point = ogr.Geometry(ogr.wkbPoint)
        #point.AddPoint(top_xy[0], top_xy[1])        # Top of line
        point.AddPoint(mid_xy[0], mid_xy[1])        # Midpoint of line
        point.Transform(coordTrans)                                      # Transform the geometry
        NodeLL[idval] = (point.GetX(), point.GetY())
        point = None
    del coords_dic
    print('  All dictionaries have been created in {0: 3.2f} seconds.'.format(time.time()-tic2))

    # Create new field in shapefile
    field_defn = ogr.FieldDefn("link", ogr.OFTInteger64)
    lyr.CreateField(field_defn)

    field_defn = ogr.FieldDefn("to", ogr.OFTInteger64)
    lyr.CreateField(field_defn)

    field_defn = ogr.FieldDefn("Order_", ogr.OFTInteger)
    lyr.CreateField(field_defn)

    field_defn = ogr.FieldDefn("GageID", ogr.OFTString)
    field_defn.SetWidth(15)
    lyr.CreateField(field_defn)

    field_defn = ogr.FieldDefn("LakeID", ogr.OFTInteger64)
    lyr.CreateField(field_defn)

    field_defn = ogr.FieldDefn("length", ogr.OFTReal)
    lyr.CreateField(field_defn)

    field_defn = ogr.FieldDefn("Slope", ogr.OFTReal)
    lyr.CreateField(field_defn)

    field_defn = ogr.FieldDefn("TopElev", ogr.OFTReal)
    lyr.CreateField(field_defn)

    # Iterate over shapefile to add new values to the newly created field
    for feature in lyr:
        link_id = int(feature.GetField(id_field))
        feature.SetField("link", link_id)
        feature.SetField("to", to_from_dic.get(link_id, 0))
        feature.SetField("Order_", StrOrder[link_id])
        feature.SetField("GageID", str(linkID_gage.get(link_id, None)))
        feature.SetField("LakeID", NoDataVal)
        feature.SetField("length", Lengths[link_id])
        feature.SetField("Slope", slope_dic[link_id])
        feature.SetField("TopElev", NodeElev[link_id])
        lyr.SetFeature(feature)
    data_source = feature = lyr = None
    print('  Fields have been added to the shapefile.')

    # We need to define the projection for the streams file: this is not done automatically by Whitebox.
    define_projection(streams_vector_file, grid_obj.proj)

    # Added 8/16/2020 because a value of 0 exists in the order list
    order.remove(0)

    # Call function to build the netCDF parameter table
    build_RouteLink(RoutingNC, order, to_from_dic, NodeElev, NodeLL, NodeXY, Lengths, StrOrder, slope_dic, gageDict=linkID_gage)
    del linkID_gage, order, to_from_dic, NodeElev, Lengths, StrOrder, NodeLL, NodeXY, slope_dic
    print('Reach-based routing inputs generated in {0:3.2f} seconds.'.format(time.time()-tic1))
    return rootgrp

def build_LAKEPARM(LakeNC, min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals):
    '''
    8/10/2017: This function is designed to build the lake parameter netCDF file.
                Ideally, this will be the only place the produces the file, and
                all functions wishing to write the file will reference this function.
    '''
    tic1 = time.time()
    min_elev_keys = list(min_elevs.keys())                                      # 5/31/2019: Supporting Python3

    # Create Lake parameter file
    print('    Starting to create lake parameter table.')
    print('        Lakes Table: {0} Lakes'.format(len(list(areas.keys()))))

    # Create NetCDF output table
    rootgrp = netCDF4.Dataset(LakeNC, 'w', format=outNCType)

    # Create dimensions and set other attribute information
    dim1 = 'feature_id'
    dim = rootgrp.createDimension(dim1, len(min_elevs))

    # Create coordinate variables
    ids = rootgrp.createVariable('lake_id','i4',(dim1))                         # Variable (32-bit signed integer)
    ids[:] = numpy.array(min_elev_keys)                                         # Variable (32-bit signed integer)

    # Create fixed-length variables
    LkAreas = rootgrp.createVariable('LkArea','f8',(dim1))                      # Variable (64-bit floating point)
    LkMxEs = rootgrp.createVariable('LkMxE', 'f8', (dim1))                      # Variable (64-bit floating point)
    WeirCs = rootgrp.createVariable('WeirC', 'f8', (dim1))                      # Variable (64-bit floating point)
    WeirLs = rootgrp.createVariable('WeirL', 'f8', (dim1))                      # Variable (64-bit floating point)
    OrificeCs = rootgrp.createVariable('OrificeC', 'f8', (dim1))                # Variable (64-bit floating point)
    OrificeAs = rootgrp.createVariable('OrificeA', 'f8', (dim1))                # Variable (64-bit floating point)
    OrificeEs = rootgrp.createVariable('OrificeE', 'f8', (dim1))                # Variable (64-bit floating point)
    lats = rootgrp.createVariable('lat', 'f4', (dim1))                          # Variable (32-bit floating point)
    longs = rootgrp.createVariable('lon', 'f4', (dim1))                         # Variable (32-bit floating point)
    Times = rootgrp.createVariable('time', 'f8', (dim1))                        # Variable (64-bit floating point)
    WeirEs = rootgrp.createVariable('WeirE', 'f8', (dim1))                      # Variable (64-bit floating point)
    AscendOrder = rootgrp.createVariable('ascendingIndex', 'i4', (dim1))        # Variable (32-bit signed integer)
    ifd = rootgrp.createVariable('ifd', 'f4', (dim1))                           # Variable (32-bit floating point)
    if WRFH_Version >= 5.2:
        Dam = rootgrp.createVariable('Dam_Length', 'i4', (dim1))                # Variable (32-bit signed integer)

    # Add CF-compliant coordinate system variable
    if pointCF:
        sr = osr.SpatialReference()                                             # Build empty spatial reference object
        sr.ImportFromProj4(wgs84_proj4)                                    # Imprort from proj4 to avoid EPSG errors (4326)
        projEsri = sr.Clone()
        projEsri.MorphToESRI()                                                      # Alter the projection to Esri's representation of a coordinate system
        PE_string = projEsri.ExportToWkt().replace("'", '"')                        # INVESTIGATE - this somehow may provide better compatability with Esri products?
        grid_mapping = crsVar
        rootgrp = add_CRS_var(rootgrp, sr, 0, grid_mapping, 'latitude_longitude', PE_string)

    # Set variable descriptions
    ids.long_name = 'Lake ID'
    LkAreas.long_name = 'Lake area (sq. km)'
    LkMxEs.long_name = 'Maximum lake elevation (m ASL)'
    WeirCs.long_name = 'Weir coefficient'
    WeirLs.long_name = 'Weir length (m)'
    OrificeCs.long_name = 'Orifice coefficient'
    OrificeAs.long_name = 'Orifice cross-sectional area (sq. m)'
    OrificeEs.long_name = 'Orifice elevation (m ASL)'
    WeirEs.long_name = 'Weir elevation (m ASL)'
    lats.long_name = 'latitude of the lake centroid'
    longs.long_name = 'longitude of the lake centroid'
    AscendOrder.long_name = 'Index to use for sorting IDs (ascending)'
    ifd.long_name = 'Initial fraction water depth'
    longs.units = 'degrees_east'                                                # For compliance
    lats.units = 'degrees_north'                                                # For compliance
    longs.standard_name = 'longitude'                                           # For compliance
    lats.standard_name = 'latitude'                                             # For compliance
    Times.standard_name = 'time'                                                # For compliance
    Times.long_name = 'time of measurement'                                     # For compliance
    Times.units = 'days since 2000-01-01 00:00:00'                              # For compliance. Reference time arbitrary
    WeirEs.units = 'm'
    ids.cf_role = "timeseries_id"                                               # For compliance
    if WRFH_Version >= 5.2:
        Dam.long_name = 'Dam length (multiplier on weir length)'

    # Apply grid_mapping and coordinates attributes to all variables
    for varname, ncVar in rootgrp.variables.items():
        if dim1 in ncVar.dimensions and varname not in ['alt', 'lat', 'lon', 'x', 'y']:
            ncVar.setncattr('coordinates', 'lat lon')                           # For CF-compliance
            if pointCF:
                ncVar.setncattr('grid_mapping', grid_mapping)                       # For CF-compliance
        del ncVar, varname

    # Fill in global attributes
    rootgrp.featureType = 'timeSeries'                                          # For compliance
    rootgrp.history = 'Created %s' %time.ctime()

    print('        Starting to fill in lake parameter table NC file.')
    AscendOrder[:] = numpy.argsort(ids[:])                                  # Use argsort to give the ascending sort order for IDs. Added by KMS 4/4/2017
    LkAreas[:] = numpy.array([float(areas[lkid]) for lkid in min_elev_keys])
    LkMxEs[:] = numpy.array([max_elevs[lkid] for lkid in min_elev_keys])
    WeirCs[:] = WeirC
    WeirLs[:] = WeirL
    OrificeCs[:] = OrificeC
    OrificeAs[:] = OrificA
    Times[:] = 0
    OrificeEs[:] = numpy.array([OrificEs.get(lkid, 0) for lkid in min_elev_keys])   # Orifice elevation is 1/3 between 'min' and max lake elevation.
    lats[:] = numpy.array([cen_lats[lkid] for lkid in min_elev_keys])
    longs[:] = numpy.array([cen_lons[lkid] for lkid in min_elev_keys])
    WeirEs[:] = numpy.array([WeirE_vals.get(lkid,0) for lkid in min_elev_keys])    # WierH is 0.9 of the distance between the low elevation and max lake elevation
    ifd[:] = ifd_Val
    if WRFH_Version >= 5.2:
        Dam[:] = dam_length

    # Close file
    rootgrp.close()
    print('        Done writing {0} table to disk.'.format(LK_nc))
    return

def add_reservoirs(rootgrp, projdir, fac, in_lakes, grid_obj, lakeIDfield=None, Gridded=True):
    """
    This function is intended to add reservoirs into the model grid stack, such
    that the channelgrid and lake grids are modified to accomodate reservoirs and
    lakes.

    This version does not attempt to subset the lakes by a size threshold, nor
    does it filter based on FTYPE.

    2/23/2018:
        Change made to how AREA paramter is calculated in LAKEPARM. Previously, it
        was based on the gridded lake area. Now it is based on the AREASQKM field
        in the input shapefile. This change was made because in NWM, the lakes
        are represented as objects, and are not resolved on a grid.

    10/15/2019:
        Adding a GRIDDED flag. If True, a gridded WRF-Hydro run is assumed, and
        no lakes will be added to LAKEPARM.nc which are not resolved on the
        routing grid. If False, all lakes coincident with the domain boundary
        will be included in LAKEPARM.nc.
    """
    #(rootgrp, projdir, fac, in_lakes, grid_obj, lakeIDfield, Gridded) = (rootgrp2, projdir, fac, in_lakes, fine_grid, None, gridded)
    tic1 = time.time()                                                          # Set timer
    subsetLakes = True                                                          # Option to eliminate lakes that do not intersect channel network
    print('      Adding reservoirs to routing stack.')
    print('      Gridded: {0}'.format(Gridded))

    # Setup Whitebox tools
    wbt = WhiteboxTools()
    wbt.verbose = False
    wbt.work_dir = projdir

    # Outputs
    LakeNC = os.path.join(projdir, LK_nc)
    outshp = os.path.join(projdir, LakesSHP)                                    # Clipped and projected input lakes shapefile
    frxst_FC = os.path.join(projdir, 'Lake_outlets.shp')
    snapPour = 'Lake_snapped_pour_points.shp'                                   # Pour points snapped downstream of lake outlets

    # Setup coordinate transform for calculating lat/lon from x/y
    wgs84_proj = osr.SpatialReference()
    wgs84_proj.ImportFromProj4(wgs84_proj4)

    # Added 11/19/2020 to allow for GDAL 3.0 changes to the order of coordinates in transform
    if int(osgeo.__version__[0]) >= 3:
        # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
        wgs84_proj.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
    coordTrans = osr.CoordinateTransformation(grid_obj.proj, wgs84_proj)        # Transformation from grid projection to WGS84
    coordTrans_inv = osr.CoordinateTransformation(wgs84_proj, grid_obj.proj)    # Transformation from WGS84 to grid projection

    # Use extent of the template raster to add a feature layer of lake polygons
    geom = grid_obj.boundarySHP('', 'MEMORY')                                   # Get domain extent for cliping geometry
    lake_ds, lake_layer, fieldNames = project_Features(in_lakes, grid_obj.proj, clipGeom=geom)
    geom = lake_layer = None

    # Save to disk in order to use in the FeatToRaster function.
    lake_ds  = ogr.GetDriverByName(VectorDriver).CopyDataSource(lake_ds, outshp)
    lake_layer = lake_ds.GetLayer()

    # Add and re-calculate area information for new, clipped and reprojected lake polygons
    lake_layer.CreateField(ogr.FieldDefn('AREASQKM', ogr.OFTReal))              # Add a single field to the new layer
    for feature in lake_layer:
        geometry = feature.GetGeometryRef()                                     # Get the geometry object from this feature
        feature.SetField('AREASQKM', float(geometry.Area()/1000000.0))          # Add an area field to re-calculate area
        lake_layer.SetFeature(feature)
        feature = geometry = None
    lake_layer.ResetReading()

    # Assign a new ID field for lakes, numbered 1...n. Add field to store this information if necessary
    if lakeIDfield is None:
        print('    Adding auto-incremented lake ID field (1...n)')
        lakeID = "newID"
        if lakeID not in fieldNames:
            lake_layer.CreateField(ogr.FieldDefn(lakeID, ogr.OFTInteger))
        for num,feature in enumerate(lake_layer):
            feature.SetField(lakeID, num+1)      # Add an area field to re-calculate area
            lake_layer.SetFeature(feature)
        lake_layer.ResetReading()
        feature = None
    else:
        print('    Using provided lake ID field: {0}'.format(lakeIDfield))
        lakeID = lakeIDfield                                                    # Use existing field specified by 'lakeIDfield' parameter

    # Generate dictionary of areas, and centroid lat/lon for populatingLAKEPARM.nc
    print('    Starting to gather lake centroid and area information.')
    cen_lats = {}
    cen_lons = {}
    areas = {}
    for feature in lake_layer:
        idval = feature.GetField(lakeID)
        areas[idval] = feature.GetField('AREASQKM')
        centroid = feature.GetGeometryRef().Centroid()
        centroid.Transform(coordTrans)                                      # Transform the geometry
        cen_lats[idval] = centroid.GetY()
        cen_lons[idval] = centroid.GetX()
        feature = centroid = None
    lake_layer.ResetReading()
    print('    Done gathering lake centroid information.')
    lakeIDList = list(areas.keys())

    # Convert lake geometries to raster geometries on the model grid
    LakeRaster = FeatToRaster(outshp, fac, lakeID, gdal.GDT_Int32, NoData=NoDataVal)
    Lake_arr = BandReadAsArray(LakeRaster.GetRasterBand(1))                     # Read raster object into numpy array
    LakeRaster = None
    Lake_arr[Lake_arr==0] = NoDataVal                                           # Convert 0 to WRF-Hydro NoData

    # Code-block to eliminate lakes that do not coincide with active channel cells
    strm_arr = rootgrp.variables['CHANNELGRID'][:]                              # Read channel grid array from Fulldom
    lake_uniques = numpy.unique(Lake_arr[Lake_arr!=NoDataVal])
    if subsetLakes:
        Lk_chan = {lake:strm_arr[numpy.logical_and(Lake_arr==lake, strm_arr==0)].shape[0]>0 for lake in lake_uniques}   # So slow...
        old_Lk_count = lake_uniques.shape[0]
        lake_uniques = numpy.array([lake for lake,val in Lk_chan.items() if val])   # New set of lakes to use
        new_Lk_count = lake_uniques.shape[0]
        Lake_arr[~numpy.isin(Lake_arr, lake_uniques)] = NoDataVal               # Remove lakes from Lake Array that are not on channels
        print('    Found {0} lakes on active channels. Lost {1} lakes that were not on active channels.'.format(new_Lk_count, old_Lk_count-new_Lk_count))
        del Lk_chan, old_Lk_count, new_Lk_count

        # Reset the 1...n index and eliminate lakes from shapefile that were eliminated here
        #num = 1                                                                 # Initialize the lake ID counter
        print('    Removing lakes not on gridded channel network')
        for feature in lake_layer:
            idval = feature.GetField(lakeID)
            if idval not in lake_uniques.tolist():
                #print('      Removing lake: {0}'.format(idval))
                lake_layer.DeleteFeature(feature.GetFID())
            else:
                #feature.SetField(lakeID, num)      # Add an area field to re-calculate area
                #lake_layer.SetFeature(feature)
                #num += 1
                pass
        lake_layer.ResetReading()
    lake_ds = lake_layer = None

    # Save the gridded lake array to the Fulldom file
    if Gridded:
        rootgrp.variables['LAKEGRID'][:] = Lake_arr                             # Write array to output netCDF file
    else:
        rootgrp.variables['LAKEGRID'][:] = NoDataVal                            # No need to populate LAKEGRID if not using gridded lakes
    print('    Process: LAKEGRID written to output netCDF.')

    # Find the maximum flow accumulation value for each lake
    flac_arr = rootgrp.variables['FLOWACC'][:]                                  # Read flow accumulation array from Fulldom
    flac_max = {lake:flac_arr[Lake_arr==lake].max() for lake in lake_uniques}

    # Iterate over lakes, assigning the outlet pixel to the lake ID in channelgrid
    strm_arr[Lake_arr>0] = NoDataVal                                            # Set all lake areas to WRF-Hydro NoData value under these lakes
    for lake,maxfac in flac_max.items():
        strm_arr[numpy.logical_and(Lake_arr==lake, flac_arr==maxfac)] = lake    # Set the lake outlet to the lake ID in Channelgrid
    del flac_arr
    if Gridded:
        rootgrp.variables['CHANNELGRID'][:] = strm_arr
    print('    Process: CHANNELGRID written to output netCDF.')

    # Now march down a set number of pixels to get minimum lake elevation
    strm_arr[strm_arr<1] = NoDataVal                                            # Remove channels (active and inactive)
    ds = array_to_points(strm_arr, ogr.OFTInteger, grid_obj.GeoTransform(), grid_obj.proj)
    out_ds = ogr.GetDriverByName(VectorDriver).CopyDataSource(ds, frxst_FC)    # Copy to file on disk
    ds = out_ds = None
    del strm_arr, ds, out_ds

    tolerance = grid_obj.DX * LK_walker                                         # Snap distance is the horizontal cellsize multiplied by number of pixels to 'walk' downstream
    snapPourFile = os.path.join(projdir, snapPour)
    wbt.snap_pour_points(frxst_FC, fac, snapPour, tolerance)                    # Snap pour points to flow accumulation grid within a tolerance

    # Gathering maximum elevation from input DEM
    fill_arr = rootgrp.variables['TOPOGRAPHY'][:]                               # Read elevation array from Fulldom
    max_elevs = {lake:fill_arr[Lake_arr==lake].max() for lake in lake_uniques}
    del Lake_arr, lake_uniques

    # Read the shapefile from previous Snap Pour Points and extract the values directly from the grid
    snap_ds = ogr.Open(snapPourFile, 0)
    pointlyr = snap_ds.GetLayer()                                               # Get the 'layer' object from the data source
    min_elevs = {}
    for feature in pointlyr:
        idval = feature.GetField('VALUE')
        point = feature.GetGeometryRef()
        row, col = grid_obj.xy_to_grid_ij(point.GetX(), point.GetY())
        min_elevs[idval] = fill_arr[row, col]
        feature = point = None
        del idval, row, col
    snap_ds = pointlyr = None
    ogr.GetDriverByName(VectorDriver).DeleteDataSource(snapPourFile)
    ogr.GetDriverByName(VectorDriver).DeleteDataSource(frxst_FC)

    # Only add in 'missing' lakes if this is a reach-routing simulation and lakes
    # don't need to be resolved on the grid.
    if not Gridded:
        # 2/23/2018: Find the missing lakes and sample elevation at their centroid.
        min_elev_keys = list(min_elevs.keys())
        print('    Lakes in minimum elevation dict: {0}'.format(len(min_elev_keys))) # Delete later
        MissingLks = [item for item in lakeIDList if item not in min_elev_keys]     # 2/23/2018: Find lakes that were not resolved on the grid
        if len(MissingLks) > 0:
            print('    Found {0} lakes that could not be resolved on the grid: {1}\n      Sampling elevation from the centroid of these features.'.format(len(MissingLks), str(MissingLks)))
            centroidElev = {}
            for idval in MissingLks:
                centroid = ogr.Geometry(ogr.wkbPoint)
                centroid.AddPoint(cen_lons[idval], cen_lats[idval])
                centroid.Transform(coordTrans_inv)                              # Transform the geometry
                row, col = grid_obj.xy_to_grid_ij(centroid.GetX(), centroid.GetY())
                centroidElev[idval] = fill_arr[row, col]
                centroid = None
                del row, col, centroid, idval

            # Update dictionaries with information on the lakes that were not resolved on the grid
            max_elevs.update(centroidElev)                                      # Add single elevation value as max elevation
            min_elevs.update(centroidElev)                                      # Make these lakes the minimum depth
            del centroidElev
        del fill_arr, lakeIDList, MissingLks

    # Give a minimum active lake depth to all lakes with no elevation variation
    elevRange = {key:max_elevs[key]-val for key,val in min_elevs.items()}   # Get lake depths
    noDepthLks = {key:val for key,val in elevRange.items() if val<minDepth}     # Make a dictionary of these lakes
    if len(noDepthLks) > 0:
        print('    Found {0} lakes with no elevation range. Providing minimum depth of {1}m for these lakes.'.format(len(noDepthLks), minDepth))
        min_elevs.update({key:max_elevs[key]-minDepth for key,val in noDepthLks.items() if val==0 }) # Give these lakes a minimum depth
        noDepthFile = os.path.join(projdir, minDepthCSV)
        with open(noDepthFile,'w') as f:
            w = csv.writer(f)
            for item in noDepthLks.items():
                w.writerows(noDepthLks.items())
            del noDepthFile
    min_elev_keys = list(min_elevs.keys())
    del elevRange, noDepthLks

    # Calculate the Orifice and Wier heights
    # Orifice elevation is 1/3 between the low elevation and max lake elevation
    OrificEs = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x])/3)) for x in min_elev_keys}

    # WierH is 0.9 of the distance between the low elevation and max lake elevation
    WeirE_vals = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x])*0.9)) for x in min_elev_keys}
    del min_elev_keys

    ##    # Dissolve multiple features to multipolygon feature layer by field value.
    ##    ds1, Layer = raster_to_polygon(LakeRaster, grid_obj.proj, geom_typ=ogr.wkbMultiPolygon)
    ##    ds2 = dissolve_polygon_to_multipolygon(ds1, Layer, 'RASTERVALU')
    ##    out_shp = ogr.GetDriverByName(VectorDriver).CopyDataSource(ds2, out_lake_raster_shp)
    ##    ds1 = Layer = out_shp = None

    # Call function to build lake parameter netCDF file
    build_LAKEPARM(LakeNC, min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals)
    print('    Lake parameter table created without error in {0: 3.2f} seconds.'.format(time.time()-tic1))
    return rootgrp, lakeID

def getxy(ds):
    """
    This function will use the affine transformation (GeoTransform) to produce an
    array of X and Y 1D arrays. Note that the GDAL affine transformation provides
    the grid cell coordinates from the upper left corner. This is typical in GIS
    applications. However, WRF uses a south_north ordering, where the arrays are
    written from the bottom to the top.
    The input raster object will be used as a template for the output rasters.
    """
    print('    Starting Process: Building to XMap/YMap')
    nrows = ds.RasterYSize
    ncols = ds.RasterXSize
    xMin, DX, xskew, yMax, yskew, DY = ds.GetGeoTransform()
    del ds, xskew, yskew
    # Build i,j arrays
    j = numpy.arange(nrows) + float(0.5)                                        # Add 0.5 to estimate coordinate of grid cell centers
    i = numpy.arange(ncols) + float(0.5)                                        # Add 0.5 to estimate coordinate of grid cell centers
    # col, row to x, y   From https://www.perrygeo.com/python-affine-transforms.html
    x = (i * DX) + xMin
    y = (j * DY) + yMax
    del i, j, DX, DY, xMin, yMax
    # Create 2D arrays from 1D
    xmap = numpy.repeat(x[numpy.newaxis, :], y.shape, 0)
    ymap = numpy.repeat(y[:, numpy.newaxis], x.shape, 1)
    del x, y
    print('    Conversion of input raster to XMap/YMap completed without error.')
    return xmap, ymap

def flip_dim(array_dimensions, DimToFlip='south_north'):
    '''
    Function to flip a dimension based on provided dimension names.

        array_dimensions - A list of dimension names for the input dataset
        DimToFlip - The dimension to reverse.
    '''

    # Determine how to slice the array in order to fit into the netCDF
    ind = [slice(None)] * len(array_dimensions)                                 # Build array slice as default (:)

    # Flip a dimension if necessary
    if DimToFlip in array_dimensions:
        flipIdx = array_dimensions.index(DimToFlip)
        ind[flipIdx] = slice(None,None,-1)
        print("    Reversing order of dimension '{0}'".format(array_dimensions[flipIdx]))
        del flipIdx
    else:
        print("    Requested dimension for reversal not found '{0}'.".format(DimToFlip))
    return ind

def move_downstream(DIRECTION, trim=True, mask=slice(None)):
    '''
    This function accepts a flow direction grid, and will return a grid of indices
    for the downstream cells for each cell on the grid. Index values that fall off
    of the grid will be optionally left the same as the input.
    '''

    # Get the indices of the grid (in Fulldom, the order is (y, x))
    j_size, i_size = DIRECTION.shape
    j, i = numpy.indices(DIRECTION.shape)

    # Mask arrays as necessary
    j = j[mask]
    i = i[mask]
    DIRECTION = DIRECTION[mask]

    # Make copies that can be modified
    downstream_index_grid_j = j.copy()
    downstream_index_grid_i = i.copy()

    # For each flow direction, adjust index up or down, left or right
    #   Remember north is up so -1 on the j-grid index is north (up)
    #   Remember north is up so +1 on the j-grid index is south (down)
    downstream_index_grid_j[numpy.logical_and(DIRECTION>=32, DIRECTION<=128)] -= 1
    downstream_index_grid_j[numpy.logical_and(DIRECTION>=2, DIRECTION<=8)] += 1
    downstream_index_grid_i[numpy.logical_and(DIRECTION>=8, DIRECTION<=32)] -= 1
    downstream_index_grid_i[(DIRECTION==2) | (DIRECTION==1) | (DIRECTION==128)] += 1

    # Calculate any invalid indices (cells that flow off the grid)
    valid_mask = numpy.full(DIRECTION.shape, True)
    valid_mask = numpy.logical_and(valid_mask, downstream_index_grid_j>=0)
    valid_mask = numpy.logical_and(valid_mask, downstream_index_grid_j<=j_size-1)
    valid_mask = numpy.logical_and(valid_mask, downstream_index_grid_i>=0)
    valid_mask = numpy.logical_and(valid_mask, downstream_index_grid_i<=i_size-1)

    # Trim the edges so that nothing is beyond the indices of the initial grid
    if trim:
        downstream_index_grid_i[~valid_mask] = i[~valid_mask]
        downstream_index_grid_j[~valid_mask] = j[~valid_mask]

    del j_size, i_size, j, i
    return downstream_index_grid_j, downstream_index_grid_i, valid_mask

def get_tot_chan_and_lakes(CH_NETRT, DIRECTION, CH_nodata=-9999):
    '''
    Numpy arrays from top to bottom, so we need to reverse all j indices
    in the function below
    '''
    print('          Obtaining valid channel cells')
    tic1 = time.time()
    error_cells = []

    # Empty grid
    CH_NETLNK = numpy.full(CH_NETRT.shape, CH_nodata)

    # Vectorized approach
    channelgrid_mask = CH_NETRT>=0

    # For each channel grid, look up what is downstream
    downstream_index_grid_j, downstream_index_grid_i, valid_mask = move_downstream(DIRECTION, trim=False, mask=channelgrid_mask)

    # Number of cells that are not flowing off the grid
    counter2 = valid_mask.sum()

    # Find where the valid downstream channels are also channel cells
    cnt = channelgrid_mask[downstream_index_grid_j[valid_mask], downstream_index_grid_i[valid_mask]].sum()

    # Find flow direction errors
    channel_error_mask = numpy.logical_and(channelgrid_mask, DIRECTION==0)
    error_cells += numpy.asarray(numpy.where(channel_error_mask)).T.tolist()

    # Adjust counters to eliminate these error cells
    cnt -= channel_error_mask.sum()
    counter2 -= channel_error_mask.sum()

    print('        found type 0 nodes      {0}'.format(cnt))
    print('            Found {0} cells that are within parameters'.format(counter2))
    print('            total number of channel elements: {0}'.format(cnt))
    print('            Completed cnt calculation in {0:3.2f} seconds'.format(time.time()-tic1))

    # Reduce all arrays to eliminate problem pixels
    valid_mask = valid_mask[~channel_error_mask[channelgrid_mask]]
    downstream_index_grid_j = downstream_index_grid_j[~channel_error_mask[channelgrid_mask]]
    downstream_index_grid_i = downstream_index_grid_i[~channel_error_mask[channelgrid_mask]]
    channelgrid_mask[channel_error_mask] = False

    CH_NETLNK.flat[numpy.flatnonzero(channelgrid_mask)[valid_mask]] = numpy.arange(1, counter2+1)
    CH_NETLNK[channel_error_mask] = CH_nodata                                   # Redundant with previous step?

    # Reverse the valid mask, so that we are only looking at the cells that flow off the grid
    valid_mask2 = ~valid_mask
    counter3 = (valid_mask2).sum()

    # Using multiple levels of boolean indexing
    CH_NETLNK.flat[numpy.flatnonzero(channelgrid_mask)[valid_mask2]] = numpy.arange(cnt+1, (cnt+1)+(counter3+1))
    cnt += counter3
    counter2 += counter3

    # Still need to add in the channels that drain to a non channel
    valid_mask3 = CH_NETRT[downstream_index_grid_j[valid_mask], downstream_index_grid_i[valid_mask]]<0
    valid_mask4 = valid_mask
    valid_mask4[valid_mask] = valid_mask3
    counter4 = (valid_mask4).sum()
    #CH_NETRT.flat[numpy.flatnonzero(channelgrid_mask)[valid_mask4]] = numpy.arange(cnt+1, (cnt+1)+(counter3+1))
    CH_NETLNK.flat[numpy.flatnonzero(channelgrid_mask)[valid_mask4]] = numpy.arange(cnt+1, (cnt+1)+(counter3+1))
    counter2 += counter4
    cnt += counter4
    del downstream_index_grid_j, downstream_index_grid_i, channelgrid_mask

    print('            Found {0} cells that are within parameters'.format(counter3+counter4))
    print('            total number of channel elements: {0}'.format(cnt))
    print('            Completed lake and boundary calculation in {0:3.2f} seconds'.format(time.time()-tic1))
    return cnt, CH_NETLNK, error_cells

def nlinks_checker(rootgrp_FD,
                    FD_linkidVar='LINKID',
                    FD_chgridVar='CHANNELGRID',
                    FD_FD8Var='FLOWDIRECTION',
                    FD_orderVar='STREAMORDER',
                    silent=False):

    '''
    3/27/2023
        This function will perform a check of the CHANNELGRID layer in Fulldom_hires.nc
        for channel connectivity errors. These typically occur as a result of processing
        with WhiteBox Tools in areas around coastlines. This script will perform
        a similar set of checks to WRF-Hydro and resolve the issues by setting
        CHANNELGRID and FLOWDIRECTION (and STREAMORDER) values to appropriate and
        valid values.

        Inputs:
            rootgrp_FD - Fulldom_hires.nc Dataset object from netcdf4 library.
                         Must be open in write mode to make changes.

    '''
    tic1 = time.time()

    # Create a list of the valid Flow Direction valus
    valid_FDs = [64, 128, 1, 2, 4, 8, 16, 32]

    # Define the NoData value for CHANNELGRID cells
    CH_nodata = -9999

    # Flag to fix errors that are found by this script
    fix_CH = True

    # Fix the masking of the input netcdf.Dataset object.
    if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
        rootgrp_FD.set_auto_mask(False)                                         # Change masked arrays to old default (numpy arrays always returned)
    variables_FD = rootgrp_FD.variables

    if FD_linkidVar in variables_FD:
        FD_linkID = variables_FD[FD_linkidVar][:]
        print('        Found {0} {1} cells in input {1} variable'.format(numpy.unique(FD_linkID[:]).shape[0], FD_linkidVar))
        del FD_linkID

    # Gather the channelgrid
    FD_chgrid = variables_FD[FD_chgridVar][:]
    FD_order = variables_FD[FD_orderVar][:]
    strm_order_min = FD_order.min()
    print('        Found {0} channelgrid cells in input {1} variable'.format((FD_chgrid>=0).sum(), FD_chgridVar))

    # Do NLINKS check - Vectorized approach
    NLINKS = 0
    print('          Obtaining total number of channel cells in domain')
    NLINKS = (FD_chgrid>=0).sum()
    print("            NLINKS IS {0}".format(NLINKS))

    # Perform check on number of links
    FD8_grid = variables_FD[FD_FD8Var][:]

    # Determine if any of the Flowdirection values are invalid and set to 0 if so
    invalid_FDs = ~numpy.in1d(FD8_grid, valid_FDs).reshape(FD8_grid.shape)
    FD8_grid[invalid_FDs] = 0
    del invalid_FDs

    # Vectorized method for obtaining errors
    cnt, CH_NETLNK, error_cells = get_tot_chan_and_lakes(FD_chgrid, FD8_grid, CH_nodata=CH_nodata)
    print("        total number of channel elements {0}".format(cnt))
    print("        total number of NLINKS           {0}".format(NLINKS))
    if cnt != NLINKS:
        print('        Apparent error in network topology {0} {1}'.format(cnt, NLINKS))

    # Perfroam check for errors after the lake and edge step
    secondary_check = True
    if secondary_check:
        CH_OUT = CH_NETLNK.copy()
        CH_OUT[CH_OUT!=CH_nodata] = 0
        error_cells2 = numpy.asarray(numpy.where(FD_chgrid!=CH_OUT)).T.tolist()
        for error_cell2 in error_cells:
            if error_cell2 not in error_cells:
                error_cells.append(error_cell2)
        del CH_OUT, error_cells2
    print('        Found {0} channel grid cells that do not match after assigning IDs.'.format(len(error_cells)))

    if len(error_cells) > 0:
        print('        Enumerating location of error cells:')
        for error_j, error_i in error_cells:
            error_lon = variables_FD['LONGITUDE'][error_j, error_i]
            error_lat = variables_FD['LATITUDE'][error_j, error_i]
            if not silent:
                print('          Direction i,j {0},{1} is invalid'.format(error_i,error_j))
                print('            Longitude/Latitude: {0},{1}'.format(error_lon, error_lat))

            if fix_CH:
                # Set the problematic channelgrid cell to NoData
                FD_chgrid[error_j, error_i] = CH_nodata

                # Set other variables related to the channel to nodata
                FD_order[error_j, error_i] = strm_order_min

                # Reset all channel cell IDs to reflect the loss of this cell
                chID = CH_NETLNK[error_j, error_i]
                CH_NETLNK[error_j, error_i] = CH_nodata
                CH_NETLNK[CH_NETLNK>chID] -= 1
            del error_lon, error_lat

    # Alter grids in the input if requested
    if not numpy.array_equal(variables_FD[FD_FD8Var][:], FD8_grid):
        print('        WARNING: The {0} variable in Fulldom_hires.nc file will be altered.'.format(FD_FD8Var))
        variables_FD[FD_FD8Var][:] = FD8_grid
    else:
        print('        The {0} variable in the input Fulldom_hires.nc file will not be altered.'.format(FD_FD8Var))
    if not numpy.array_equal(variables_FD[FD_chgridVar][:],FD_chgrid):
        print('        WARNING: The {0} variable in Fulldom_hires.nc file will be altered.'.format(FD_chgridVar))
        variables_FD[FD_chgridVar][:] = FD_chgrid
        print('        WARNING: The {0} variable in Fulldom_hires.nc file will be altered.'.format(FD_orderVar))
        variables_FD[FD_orderVar][:] = FD_order
    else:
        print('        The {0} variable in Fulldom_hires.nc file will not be altered.'.format(FD_chgridVar))
        print('        The {0} variable in Fulldom_hires.nc file will not be altered.'.format(FD_orderVar))

    # Clean up
    del variables_FD, FD_chgrid, FD_order
    print('        NLINKS checking process completed in {0:3.2f} seconds'.format(time.time()-tic1))
    return rootgrp_FD

def group_min(l, g):
    '''Function for gathering minimum value from a set of groups'''
    groups = defaultdict(int)                                                   # default value of int is 0
    for li, gi in zip(l, g):
        if li <= groups[gi] or groups[gi]==0:
            groups[gi] = li
    return groups

def set_problem(problem_lakes, LakeID, problemstr):
    '''This function is used to add elements to a dictionary where the values are
    a list and the keys may or may not exist. This speeds up adding elements to
    an existing dictionary.'''
    try:
        problem_lakes[LakeID] += [problemstr]
    except KeyError:
        problem_lakes[LakeID] = [problemstr]
    return problem_lakes

def get_inflow_segs(FLWBarr, localLk, FromComIDs, FromSegs, LakeAssociation=LakeAssoc):
    '''This function will list all inflow segments to a particular lake.
    Inputs:
        1) FLWBarr: Array containing flowlines and the associated lake
        2) localLk: The ID of the lake to examine
        3) FromComIDs: Dictionary of the From:To relationship for flowlines
        4) FromSegs: Dictionary of the upstream links for each flowline
    Globals:
        1) LakeAssociation: The field name used to associate flowlines with lakes.
    Output:
        inflows: Array of links that are inflows to the lake
    '''

    Lake_Links = FLWBarr[FLWBarr[LakeAssociation]==localLk]                           # Find all of the links inside the lake

    # Find all headwater segments as start points
    ToSeg_keys = Lake_Links[FLID]                                               # Array of all lake links for this local lake
    ToSeg_vals = numpy.array([FromComIDs.get(key) for key in ToSeg_keys])       # Array of all downstream links for all of the lake links
    ups = ToSeg_keys[~numpy.in1d(ToSeg_keys, ToSeg_vals)].tolist()              # List of all 'headwater' segments for this lake
    del ToSeg_keys, ToSeg_vals, Lake_Links

    # Gather all inflow segments from lake 'headwater' segments
    inflows = [FromSegs[up] for up in ups if up in FromSegs]                    # list of all contributing segments to the headwater segments for this lake
    inflows = numpy.array([item for sublist in inflows for item in sublist])    # Convert from list of lists to flat list
    del ups
    return inflows

def check_downstream(uplink, localLk, FromComIDs, Lake_LinkDict):
    '''This function will keep looking downstream to see if the link will flow
    back into the same lake and will stop if it encounteres a different lake.
    Inputs:
        1) uplink: ID of the flowline to examine
        2) localLk: ID of the lake under examination
        3) FromComIDs: Dictionary of the From:To relationship for flowlines.
        4) Lake_LinkDict: Dictionary of the Flowline:Lake relationships.
    Outputs:
        result: Boolean (True/False) of if there is a flow loop
        downlinks: List of links to re-associate with this lake
        counter: Number of flowlines downstream iterated over to find the loop
        '''

    # Setup initial values
    result = False                                                              # Default condition
    counter = 0                                                                 # Initiate counter
    downlinks = [uplink]                                                        # Initiate list of links including supplied link
    while not result:
        counter += 1                                                            # Advance the counter
        down = FromComIDs.get(uplink)                                           # Move down one segment
        if down in NoDownstream:
            break                                                               # Reached end of flow network. Break without returning True
        else:
            downlinks.append(down)                                              # Add this link to the list
            if down in Lake_LinkDict:
                if Lake_LinkDict.get(down) == localLk:                          # This means that the links flows eventually back into the same lake
                    result = True
                    break
                else:
                    break                                                       # Link is associated with another lake. Break without returning true.
        uplink=down                                                             # Make this segment the upstream segment
    #print '      Looked downstream %s links from %s in lake %s.' %(counter, uplink, localLk)
    return result, downlinks, counter

def get_lake_routing_info(FLWBarr, localLk, Lake_LinkDict, accum_val, FromComIDs, FromSegs, LakeAssociation=LakeAssoc):
    '''This function will list all inflow segments to a "lake" as well as additional
    datasets to assist in routing through the lake.

    Inputs:
        1) FLWBarr: Array containing flowlines and the associated lake
        2) localLk: The ID of the lake to be examined
        3) Lake_LinkDict: Dictionary of all link:lake associations in the domain (pass through to check_downstream function)
        4) accum_val: The value given to each segment before performing lake routed accumulations
        5) FromComIDs: Dictionary of the From:To relationship for flowlines
        6) FromSegs: Dictionary of the upstream links for each flowline
    Globals:
        1) LakeAssociation: A field describing which lake ComID a flowline is associated with
        2) FLID: The flowline ComID
    Outputs:
        1) Lake_LinksList: List of lake link ComIDs
        2) inflows: list of all contributing segments to the headwater segments for this lake
        3) SegVals: A dictionary initializing the accumulation of flow for each link in this lake
        4) SegVals2: A dictionary initializing the accumulation of flow for each link in this lake. This dict is used to find the outlet.
        5) ups: List of all 'headwater' segments for this lake
        6) newLakeLinks: Dictionary to store newly found flowline:lake associations

    If UseAll == True, the script will look outside of the link:lake associations
    to find segments that may exit the lake and then re-enter. If the flow re-enters
    the lake, it will include those segments.

    Input array FLWBarr must have fields [LakeAssociation, FLID], as defined in global variables.
    '''

    # Local variables
    UseAll = True                                                               # Switch to look outside of just the links that spatially intersect the lake

    # Use all links that are associated with lakes (either through spatial join or attribute join)
    Lake_Links = FLWBarr[FLWBarr[LakeAssociation]==localLk]                           # Find all of the links associated with the current lake
    Lake_LinksList = Lake_Links[FLID].tolist()                                  # List of lake link ComIDs

    # Find all headwater segments as start points
    ToSeg_keys = numpy.unique(Lake_Links[FLID])                                 # Array of all unique lake links for this local lake
    ToSeg_vals = numpy.unique(numpy.array([FromComIDs.get(key) for key in ToSeg_keys])) # Array of all unique downstream links for all of the lake links
    ups = ToSeg_keys[~numpy.in1d(ToSeg_keys, ToSeg_vals)].tolist()              # List of all 'headwater' segments for this lake

    # Find any non-assigned lake segments in the lake flow network by iterating down the network
    newLakeLinks = {}                                                           # Dictionary to store new flowline:lake associations
    counter = 0                                                                 # Initiate counter
    seen = []                                                                   # Initiate list of seen links
    if UseAll:

        # Added 12/30/2019 to avoid situations with no upstream or downstream links
        if len(ups) == 0:                                                       # I am not sure this will ever get triggered
            changes = 0
        else:
            changes = 1                                                         # Initiate the change detector in order to move down the network in the first iteration

        # Search down the network of links in this lake to find the outlet
        iteration = 0                                                           # Initiate the iteration counter
        while changes > 0:
            downs = list(set([FromComIDs.get(key) for key in ups]))             # Get unique list of downstream flowline ComIDs

            for key in downs:
                if key in NoDownstream:                                         # Downstream segment is not a valid ComID
                    # This might be where to catch an endorheic basin terminal lake
                    downs.remove(key)

            # Add any link:lake associations that were not already present. This is an attempt to eliminate flow looping out of lakes and back in.
            checklinks = list(set([item for item in downs if item not in Lake_LinksList]))  # Unique list of links that are downstream of lake links but not associated with that lake

            # Check these suspect links for downstream connectivity back to the lake
            if len(checklinks) > 0:
                for uplink in checklinks:
                    if uplink in seen:
                        continue
                    result, downlinks, counter1 = check_downstream(uplink, localLk, FromComIDs, Lake_LinkDict)
                    if result:
                        # If result == True, then this lake flows out and then back into itself
                        #print '      Looked downstream %s links from %s in lake %s.' %(counter1, uplink, localLk)
                        Lake_LinksList += downlinks                             # Add these links to the lake association
                        newLakeLinks.update({item:localLk for item in downlinks})   # Store these new flowline:lake associations
                    else:
                        # This downstream segment never returns to the lake. Disassociate.
                        Lake_LinksList = [item for item in Lake_LinksList if item not in downlinks[1:]] # Remove these segments from association with this lake
                    downs.remove(uplink)                                        # Remove this segment so that the loop can complete.
                    seen.append(uplink)

            downs = [item for item in downs if item in Lake_LinksList]          # Remove downstream segments that are not associated with this lake
            #for key in downs:
            #    if key in NoDownstream:
            #        continue                                                    # Downstream segment is not a valid ComID
            ups = downs                                                         # Change the upstream segments to examine in the next iteration to the downstream ones from the previous iteration
            changes = sum([1 for item in downs if item is not None])            # Number of changes in this iteration
            iteration += 1                                                      # Add one for each level
        counter += 1                                                            # Advance the counter

    Lake_LinksList = list(set(Lake_LinksList))                                  # Remove duplicate link IDs

    # Set local accumulation values to default accum_val (1)
    SegVals = {key:accum_val for key in Lake_LinksList}                         # Set all initial values to accum_val
    SegVals2 = {key:accum_val for key in Lake_LinksList}                        # Set all initial values to accum_val

    # Test to regenerate upstream inflows to account for newly added links (2/26/2018)
    ToSeg_keys = numpy.array(Lake_LinksList)                                    # Array of all unique lake links for this local lake
    ToSeg_vals = numpy.unique(numpy.array([FromComIDs.get(key) for key in ToSeg_keys])) # Array of all unique downstream links for all of the lake links

    # Re-gather the upstream segments for this lake
    ups = ToSeg_keys[~numpy.in1d(ToSeg_keys, ToSeg_vals)].tolist()              # List of all 'headwater' segments for this lake
    del ToSeg_keys, ToSeg_vals                                                  # Free up memory

    # Gather all inflow segments from lake 'headwater' segments. This may not account for oCONUS contributing links
    #inflows = [FromSegs[up] for up in ups if up in FromSegs]                    # list of all contributing segments to the headwater segments for this lake
    #inflows = numpy.unique(numpy.array([item for sublist in inflows for item in sublist]))  # Convert from list of lists to unique flat list

    # 2/27/2018: We can more closely replicate WRF-Hydro's TYPE=3 segments if we say inflows are all segments that flow into any lake segment that are not already accounted for
    inflows = []
    for link in Lake_LinksList:
        upsegs = FromSegs.get(link)                                             # All upstream segments for each lake link
        if upsegs is None:
            continue
        upsegs2 = [item for item in upsegs if item not in Lake_LinksList]       # All contributing links that are not already associated with the lake
        inflows += upsegs2
        del upsegs, upsegs2
    inflows = numpy.array(inflows)
    return Lake_LinksList, inflows, SegVals, SegVals2, ups, newLakeLinks

def Lake_Link_Type(FLWBarr, FromComIDs, FLarr, subset=None, LakeAssociation=LakeAssoc):
    '''
    This function will assign a link type to each lake.
            3 = Lake Inflow Link
            2 = Internal Lake Link
            1 = Lake Outflow Link (based on accumulated flow in the lake)

    This function sorts the lakes by the minimum HydrSeq of all links associated
    with that lake. Thus, theoretically, lakes can be visited by increasing HydroSeq
    and no downstream searching for lakes should be necessary. ...

    This function will examine each lake and determine the outlets based on an
    accumulation of values through the topology of links within the lake. It will
    identify lakes with (potentially) multiple outlets, as well as internally
    draining lakes. Lakes with multiple outlets according to flow accumulation
    within the lake will have the minor outlets removed.

    2/16/2018: This function finally diagnoses all potential conflicts with WRF-Hydro.
        In the "for LakeID in LakeSeq[FLID]:" loop, any continue statements will
        eliminate that lake from being placed on the flow network. Conditions that
        prevent a lake from functioning in WRF-Hydro are:
            1) The lake has no Lake Link Type of 1 (outlet link). This seems to
               occur when a lake has only one interior link, or when all interior
               links flow directly to a nonexistent outlet link (endorheic).
            2) A lake flows directly into another lake.
            3) Multiple outlet links are defined
    '''

    print('        Starting to gather lake link type information.')

    # Options and defaults
    tic1 = time.time()                                                          # Initiate timer for this function
    #HydroSeq = False                                                            # Switch to use HydroSeq to find the lake outlet
    IterateList = True                                                          # Switch to iterate over lake links to find outlet(s)
    accum_val = 1                                                               # Values given to each segment before accumulation
    debug = False                                                               # Switch to trigger many, many print statements

    # Build default mapping files to store new lake mappings
    problem_lakes = {}                                                          # Problem dictionary
    Old_New_LakeComID = {}                                                      # Dictionary to store the 'old Lake ComID':'new Lake ComID' mapping
    ChainedLakes = {}                                                           # Dictionary to store number of lakes chained together
    Remove_Association = []                                                     # List used to remove a lake association from a flowline (for divergences in lakes)

    # Set output array dtype and field names
    dtype1 = dict(names=(FLID, 'LINK_TYPE', LakeAssociation, 'Accum1', 'Accum2'), formats=('<i4', '<i4', '<i4', '<i4', '<i4'))
    dtype2 = dict(names=(FLID, 'LINK_TYPE', LakeAssociation, 'Accum1', 'Accum2', 'Reason'), formats=('<i4', '<i4', '<i4', '<i4', '<i4', '<i4'))

    # Attempt to use lists instead of arrays
    Lake_Link_Type_COMID = []                                                   # Flowline ComID for flowlines associated with lakes
    Lake_Link_Type = []                                                         # The lake link type for this flowline
    Lake_Link_Type_WBAREACOMI = []                                              # The lake ComID associated with this flowline
    Lake_Link_Type_Accum1 = []                                                  # Added 2/17/2018 to keep track of lake local accumulation
    Lake_Link_Type_Accum2 = []                                                  # Added 2/17/2018 to keep track of lake local accumulation

    # Attempt to keep diagnostic information for lakes that are eliminated by this pre-processor
    Tossed_Lake_Link_Type_arr = numpy.empty(0, dtype=dtype2)                    # Generate new empty array to store tossed lake values

    # Build a dictionary here to pass into the get_lake_routing_info function
    Lake_LinkDict = {item[FLID]:item[LakeAssociation] for item in FLWBarr}            # Generate dictionary of link:lake associations

    # 1) Gather a list of all contributing segments for each segment in the system
    tic2 = time.time()
    FromSegs = {}
    for key,val in FromComIDs.items():
        try:
            FromSegs[val] += [key]
        except KeyError:
            FromSegs[val] = [key]
    print('        Completed FromSegs dictionary in {0:3.2f} seconds.'.format(time.time()-tic2))

    # 2) Create a sorting of lakes that will start with the lake which has the lowest HydroSeq value in it's flowlines
    # Use indexing to grab elements common to both arrays while preserving order of one of the arrays (FLWBarr)
    tic2 = time.time()
    commons = FLarr[numpy.in1d(FLarr[FLID], FLWBarr[FLID])]                     # An array of all of the flowlines in the flowline array that are common to the flowline-waterbody association array
    xsorted = numpy.argsort(commons[FLID])                                      # Find the sorted order for the array that is to be sorted
    ypos = numpy.searchsorted(commons[xsorted][FLID], FLWBarr[numpy.in1d(FLWBarr[FLID], FLarr[FLID])][FLID]) # Search only the common values between arrays
    indices = xsorted[ypos]
    commons2 = commons[indices]                                                 # Re-order the input array to match the comparison array order
    del commons, xsorted, ypos, indices

    # Use the group_min function to find the minimum HydroSeq value for each group of WBAREACOMID values
    group_minDict = group_min(commons2[hydroSeq], FLWBarr[LakeAssociation])     # Find the minimum HydroSeq value for this lake
    del commons2                                                                # Free up memory

    # Construct an array to store the Lake COMID and Minimum Hydrosequence
    dtype = dict(names=(FLID, 'minHydroSeq'), formats=('<i4', '<f8'))
    LakeSeq = numpy.array(list(group_minDict.items()), dtype=dtype)             # Create array of minimum HydroSeq values for each lake
    del group_minDict, dtype                                                    # Free up memory
    LakeSeq = LakeSeq[LakeSeq[FLID]>-9998]                                      # Clip off -9999, -9998
    LakeSeq.sort(order='minHydroSeq')                                           # Sort by minimum HydroSeq
    LakeSeqsize = LakeSeq.shape[0]                                              # Get the number of elements in the array
    print('        Completed sorting lakes by minimum HydroSeq in {0:3.2f} seconds.'.format(time.time()-tic2))

    # 3)  Find outflow link and assign type 1 (but not if it flows into another lake)
    counter = 0                                                                 #
    counter3 = 0                                                                # Counter to keep track of multiple outlet lakes
    counter4 = 0                                                                # Counter to keep track of headwater lakes
    counter5 = 0                                                                # Counter to keep track of terminal lakes
    counter6 = 0                                                                # Counter to keep track of lakes immediately downstream
    counter7 = 0                                                                # Counter to keep track of lakes with outlet that drians to nowhere
    counter8 = 0                                                                # Counter to keep track of lakes with no Type=2 links
    tic2 = time.time()                                                          # Reset the counter to provide progressive print statements
    seen = []                                                                   # Lakes in this list have been 'seen', or reviewed
    if subset is not None:
        LakeSeq = LakeSeq[numpy.in1d(LakeSeq[FLID], subset[FLID])]              # Subset the list of lakes to a subset array if requested
        print('        Subsetted the lakes from {0} to {1} based on a provided subset array.'.format(LakeSeqsize, LakeSeq.shape[0]))

    for LakeID in LakeSeq[FLID]:
        if debug:
            print('          Lake {0}'.format(LakeID))

        if LakeID in seen:
            #print('          Lake {0} has already been examined.'.format(LakeID))
            continue                                                            # This lake has already been examined

        # Get initial lake information, including adding looping outflows
        Lake_LinksList, inflows, SegVals, SegVals2, ups, newLakeLinks = get_lake_routing_info(FLWBarr, LakeID, Lake_LinkDict, accum_val, FromComIDs, FromSegs, LakeAssociation=LakeAssociation)
        if debug:
            print('          Lake {0} has {1} links, {2} inflows, and {3} upstream flowlines.'.format(LakeID, len(Lake_LinksList), len(inflows), len(ups)))
            print('          Lake {0} has {1} inflows: {2}'.format(LakeID, len(inflows), inflows))
            print('          Lake {0} has {1} links: {2}'.format(LakeID, len(Lake_LinksList), Lake_LinksList))

        # First, determine if this is a headwater lake. Must set 'continue statement to remove the lake
        if len(inflows) == 0:
            problem_lakes = set_problem(problem_lakes, LakeID, 'Headwater lake, no inflows')
            counter4 += 1                                                       # Advance the counter
            if debug:
                print('        Headwater lake, no inflows')

        # See if any of the inflow links are on other lakes
        counter2 = 0                                                            # Counter to keep track of lakes immediately upstream for each lake
        num_uplakes = FLWBarr[numpy.in1d(FLWBarr[FLID], inflows)][LakeAssociation].tolist()   # This will be a list of lakes that are associated with this current lake's inflows
        if LakeID in num_uplakes:                                               # Check to make sure lake doesn't flow into itself
            # If the upstream lake is the same ID as the current lake, log the issue
            num_uplakes.remove(LakeID)                                          # If it does, remove that lake to elminate a flow loop
            problem_lakes = set_problem(problem_lakes, LakeID, 'Potential Flow Loop')
            if debug:
                print('        Potential Flow Loop')
        to_alter = list(num_uplakes)                                            # Copy the list to a new list

        # If there are lakes immediately above this one, examine each one, then examine all lakes upstream of those, etc.
        # This will essentially merge any lakes that flow directly into another lake, giving the ID of the most downstream lake to all chained lakes above it.
        while len(num_uplakes) > 0:
            uplake = num_uplakes[0]                                             # Look at the first lake in the list
            if uplake in seen:                                                  # If this lake has already been examined, then skip it
                num_uplakes.remove(uplake)                                      # If this lake has been examined, remove it from this list
                continue
            Old_New_LakeComID[uplake] = LakeID                                  # This will give the upstream lake the ID of the downstream lake it flows directly into
            inflows2 = get_inflow_segs(FLWBarr, uplake, FromComIDs, FromSegs, LakeAssociation=LakeAssociation)   # Get all of the upstream lake's inflows and associate these flowlines with the new LakeID
            uplakes2 = FLWBarr[numpy.in1d(FLWBarr[FLID], inflows2)][LakeAssociation].tolist() # Get all the upstream lakes for this lake
            if uplake in uplakes2:
                uplakes2.remove(uplake)                                         # Remove any flow loops
            num_uplakes += uplakes2                                             # Add these new upstream lakes to the list of upstream lakes
            to_alter += uplakes2
            num_uplakes.remove(uplake)                                          # Remove this lake from the list so that iteration will stop when the list is empty
            counter2 += 1                                                       # Advance the counter
            seen.append(uplake)                                                 # Add this upstream lake to the list of examined lakes
            del inflows2, uplakes2                                              # Free up memory

        # Alter all Waterbody ComIDs at once to match the most downstream lake and regenerate lake information for the new, merged lake
        if len(to_alter) > 0:
            ChainedLakes[LakeID] = to_alter                                     # Add list of upstream lakes to this dictionary for this lake
            print('        Altering {0} lake COMID value(s) to {1} from: {2}'.format(len(to_alter), LakeID, to_alter))
            for item in to_alter:
                problem_lakes = set_problem(problem_lakes, LakeID, 'Upstream segment for {0} belongs to another lake [{1}]'.format(LakeID, item))
            FLWBarr[LakeAssociation][numpy.in1d(FLWBarr[LakeAssociation], numpy.array(to_alter))] = LakeID  # Change all the IDs of the upstream lakes at once to the current lake ID
            Lake_LinksList, inflows, SegVals, SegVals2, ups, newLakeLinks = get_lake_routing_info(FLWBarr, LakeID, Lake_LinkDict, accum_val, FromComIDs, FromSegs, LakeAssociation=LakeAssociation)
        Lake_LinksList2 = Lake_LinksList + inflows.tolist()                     # Add all lake links to all inflows
        if debug:
            print('          Lake {0} now has {1} links: {2}'.format(LakeID, len(Lake_LinksList2), Lake_LinksList2))

        # Assign the initial LINK_TYPE values based on local lake links and inflow links
        #Lake_Link_Type_local = [3 if item in inflows.tolist() and item not in newLakeLinks else 2 for item in Lake_LinksList2]   # Append Default LINK_TYPE (2) for all items in the Lake_LinksList
        Lake_Link_Type_local = [3 if item in inflows.tolist() else 2 for item in Lake_LinksList2]   # Append Default LINK_TYPE (2) for all items in the Lake_LinksList

        # Added 12/30/2019 to avoid situations with no upstream or downstream links
        if len(ups) == 0:
            changes = 0
        else:
            changes = 1                                                         # Initiate the change detector in order to move down the network in the first iteration

        # For each lake, search down the network of links in this lake to find the outlet flowline(s)
        iteration = 0                                                           # Initiate the iteration counter
        while changes > 0:
            downs = list(set([FromComIDs.get(key, 0) for key in ups]))             # Get unique list of downstream flowline ComIDs
            downs = [item for item in downs if item in Lake_LinksList2]         # Remove downstream segments that are not associated with this lake
            for key in downs:
                if key in NoDownstream:
                    # This might be where to catch an endorheic basin terminal lake
                    # This appears to never get triggered, probably because there are no 0 values in Lake_LinksList2
                    problem_lakes = set_problem(problem_lakes, LakeID, 'Drains to ocean or internally')
                    if debug:
                        print('        Link {0} drains to ocean or internally'.format(key))
                    counter5 += 1                                               # Advance the counter
                else:
                    # This downstream segment is a real flowline in the network
                    ups_local = FromSegs[key]                                   # Get all the upstream links for this link
                    vals = [SegVals.get(val) for val in ups_local if val in SegVals]    # Get the accumulated value for each upstream link (usually 1 for each upstream link)
                    newval = 1 + sum(vals)                                      # Add one to the sum of the accumulated values of all upstream links
                    SegVals[key] = newval                                       # Put the accumulated sum in the SegVals dictionary
                    SegVals2[key] = newval                                      # Put the accumulated sum in the SegVals2 dictionary
                    for up in ups_local:
                        SegVals2[up] = 0                                        # Set all upstream values to 0 in the SegVals2 dictionary
            ups = downs                                                         # Move down the network by one link
            #changes = sum([1 for item in downs if item not in NoDownstream])    # Quantify the number of valid downstream links
            changes = sum([0 if item in NoDownstream else 1 for item in downs])    # Quantify the number of valid downstream links
            iteration += 1                                                      # Add one for each level
        seen.append(LakeID)                                                     # Add this lake to the list of lakes that have been seen already
        if debug:
            print('        Lake_LinksList2: {0}'.format(Lake_LinksList2))

        # Record the routed and unrouted accumulation values
        Accum1_local = [0 if item in inflows.tolist() else SegVals[item] for item in Lake_LinksList2]   # Added 2/17/2018 to keep track of lake local accumulation
        Accum2_local = [0 if item in inflows.tolist() else SegVals2[item] for item in Lake_LinksList2]  # Added 2/17/2018 to keep track of lake local accumulation
        if debug:
            print('        Accum1_local: {0}'.format(Accum1_local))
            print('        Accum2_local: {0}'.format(Accum2_local))
            print('SegVals2.keys(): {0}'.format(SegVals2.keys()))

        # Are there multiple outflows? If so, remove all networks contributing to minor outlets
        # Note that this method only chooses the first link ID if multiple links have the same maximum accumulation value
        #maxflow = SegVals2.keys()[list(SegVals2.values()).index(max(SegVals2.values()))]  # Find the flowline COMID with the highest accumulation value
        maxflow = list(SegVals2)[list(SegVals2.values()).index(max(list(SegVals2.values())))]  # Find the flowline COMID with the highest accumulation value

        if debug:
            print('maxflow: {0}'.format(maxflow))

        outflows = [key for key,val in SegVals2.items() if val > 0]             # All links with accumulation still in them
        if debug:
            print('outflows: {0}'.format(outflows))

        if len(outflows) > 1:
            problem_lakes = set_problem(problem_lakes, LakeID, 'Potentially multiple outlets. Secondary Outlets: {0}'.format(outflows))
            if debug:
                print('        Potentially multiple outlets. Secondary Outlets: {0}'.format(outflows))

            # Iterate upstream from the secondary outlet and remove associations.
            non_outlets = [item for item in outflows if item != maxflow]        # Create list of secondary outlets (outlets with fewer contributing flowlines than the maxflow)
            Remove_local = []                                                   # Initiate list of associations to remove
            ups = non_outlets
            while len(ups) > 0:
                Remove_local += ups                                             # Add the segments from the secondary outlet to the association removal list
                #ups = [FromSegs.get(key) for key in ups if key in Lake_LinksList2] # Move upstream in the flow network, only considering the local lake links
                ups = [FromSegs.get(key) for key in ups if key in Lake_LinksList]   # Move upstream in the flow network, only considering the local lake links (not inflows links)
                ups = [item for sublist in ups if sublist is not None for item in sublist]  # Flatten list including None values
            counter3 += len(Remove_local)                                       # Iterate counter to keep track of lake association removals
            Remove_Association += Remove_local                                  # Add these networks that drain to a secondary lake outlet to the association removal list

            # Before eliminating this lake, save the information in a separate table
            arrsize = Tossed_Lake_Link_Type_arr.shape[0]                        # Get current length of the array
            maskList = numpy.array([True if item in Remove_local else False for item in Lake_LinksList2])   # Create a mask list to mask other lists with
            Tossed_Lake_Link_Type_arr.resize(arrsize + maskList.sum(), refcheck=False)              # Add rows to the recarray to store new data
            Tossed_Lake_Link_Type_arr[FLID][arrsize:] = numpy.array(Lake_LinksList2)[maskList]      # Add Lake_LinksList to the list
            Tossed_Lake_Link_Type_arr['LINK_TYPE'][arrsize:] = numpy.array(Lake_Link_Type_local)[maskList] # Add Lake Link Types
            Tossed_Lake_Link_Type_arr[LakeAssociation][arrsize:] = [LakeID for item in Remove_local]      # Add WBAREACOMI lake associations
            Tossed_Lake_Link_Type_arr['Accum1'][arrsize:] = numpy.array(Accum1_local)[maskList]     # Added 2/17/2018 to keep track of lake local accumulation
            Tossed_Lake_Link_Type_arr['Accum2'][arrsize:] = numpy.array(Accum2_local)[maskList]     # Added 2/17/2018 to keep track of lake local accumulation
            Tossed_Lake_Link_Type_arr['Reason'][arrsize:] = 1                   # Reason: 'Potentially multiple outlets. Secondary Outlets: %s' %(outflows)

            # Use a reverse the masklist to remove these minor network elements draining to false outlets from the rest of the data
            Lake_LinksList2 = numpy.array(Lake_LinksList2)[~maskList].tolist()              # Subset the list to exclude the networks draining to a minor outlet
            Lake_Link_Type_local = numpy.array(Lake_Link_Type_local)[~maskList].tolist()    # Subset the list to exclude the networks draining to a minor outlet
            Accum1_local = numpy.array(Accum1_local)[~maskList].tolist()                    # Subset the list to exclude the networks draining to a minor outlet
            Accum2_local = numpy.array(Accum2_local)[~maskList].tolist()                    # Subset the list to exclude the networks draining to a minor outlet

        # Find if this has an outlet or not
        if FromComIDs.get(maxflow) in NoDownstream:
            # This is either a coastline waterbody or an internally draining (endorheic) lake
            # If these lakes are not given any type=1 outlet links, the lake will not be placed into RouteLink
            problem_lakes = set_problem(problem_lakes, LakeID, 'Link with maximum flow drains to ocean or internally')
            if debug:
                print('        Link with maximum flow drains to ocean or internally')
            #Lake_Link_Type_local = [1 if com==maxflow else LT for com,LT in zip(Lake_LinksList2, Lake_Link_Type_local)]    # Test to see if lakes that flow to a nodata point can be kept in WRF-Hydro
            counter7 += 1                                                       # Advance the counter
        elif FLWBarr[FLWBarr[FLID]==FromComIDs.get(maxflow)].shape[0] > 0:      # This is a downstream lake
            # This is a real downstream segment
            problem_lakes = set_problem(problem_lakes, LakeID, 'Downstream of outlet segment is a lake segment')
            if debug:
                print('        Downstream of outlet segment is a lake segment')
            counter6 += 1                                                       # Advance the counter
        else:
            # This is a legitimate downstream flowline. Alter local lake link types to reflect the outlet link type
            Lake_Link_Type_local = [1 if com==maxflow else LT for com,LT in zip(Lake_LinksList2, Lake_Link_Type_local)]

        # Additional check: Does the lake have only Type=3 and Type=1 links in it? If so, elminate. Added 2/15/2018
        if 1 not in list(set(Lake_Link_Type_local)):
            problem_lakes = set_problem(problem_lakes, LakeID, 'This lake has no type=1 (outlet) link in it. Eliminating...')
            if debug:
                print('        This lake has no type=1 (outlet) link in it. Eliminating...')
            counter8 += 1                                                       # Advance the counter

            # Before eliminating this lake, save the information in a separate table
            arrsize = Tossed_Lake_Link_Type_arr.shape[0]                        # Get current length of the array
            Tossed_Lake_Link_Type_arr.resize(arrsize + len(Lake_LinksList2), refcheck=False)    # Add rows to the recarray to store new data
            Tossed_Lake_Link_Type_arr[FLID][arrsize:] = Lake_LinksList2         # Add Lake_LinksList to the list
            Tossed_Lake_Link_Type_arr['LINK_TYPE'][arrsize:] = Lake_Link_Type_local         # Add Lake Link Types
            Tossed_Lake_Link_Type_arr[LakeAssociation][arrsize:] = [LakeID for item in Lake_LinksList2]         # Add WBAREACOMI lake associations
            Tossed_Lake_Link_Type_arr['Accum1'][arrsize:] = Accum1_local        # Added 2/17/2018 to keep track of lake local accumulation
            Tossed_Lake_Link_Type_arr['Accum2'][arrsize:] = Accum2_local        # Added 2/17/2018 to keep track of lake local accumulation
            Tossed_Lake_Link_Type_arr['Reason'][arrsize:] = 2                   # Reason: 'This lake has no type=1 (outlet) link in it. Eliminating...'
            continue                                                            # Added 2/15/2018 as a test to eliminate these lakes from RouteLink

        # Iterate the lake counter and print a statement every so often
        counter += 1                                                            # Advance the main counter
        if counter % 1000 == 0:
            print('        {0} lakes processed in {1:3.2f} seconds.'.format(counter, time.time()-tic2))
            tic2 = time.time()

        # Store link type information in lists
        Lake_Link_Type_COMID += Lake_LinksList2                                 # Add Lake_LinksList to the list
        Lake_Link_Type_WBAREACOMI += [LakeID for item in Lake_LinksList2]       # Append LakeID for all items in the Lake_LinksList
        Lake_Link_Type += Lake_Link_Type_local
        Lake_Link_Type_Accum1 += Accum1_local                                   # Added 2/17/2018 to keep track of lake local accumulation
        Lake_Link_Type_Accum2 += Accum2_local                                   # Added 2/17/2018 to keep track of lake local accumulation
        del Lake_Link_Type_local

    # Subset Lake_Link_Type_arr to just the links in XXX
    Lake_Link_Type_arr = numpy.zeros(len(Lake_Link_Type_COMID), dtype=dtype1)
    Lake_Link_Type_arr[FLID] = Lake_Link_Type_COMID                             # Populate with lake COMIDs
    Lake_Link_Type_arr['LINK_TYPE'] = Lake_Link_Type                            # Add Lake Link Types
    Lake_Link_Type_arr[LakeAssociation] = Lake_Link_Type_WBAREACOMI                   # Add WBAREACOMI lake associations
    Lake_Link_Type_arr['Accum1'] = Lake_Link_Type_Accum1                        # Add unrouted lake accumulation
    Lake_Link_Type_arr['Accum2'] = Lake_Link_Type_Accum2                        # Add routed lake accumulation
    del Lake_Link_Type_COMID, Lake_Link_Type, Lake_Link_Type_WBAREACOMI, dtype1, dtype2, Lake_LinkDict

    # Remove WBAREACOMI association for the multiple outlets (other than that with max flow)
    FLWBarr = FLWBarr[~numpy.in1d(FLWBarr[FLID], numpy.array(Remove_Association))]   # Remove from the array any items that need the association removed
    Lake_Link_Type_arr[LakeAssociation][numpy.in1d(Lake_Link_Type_arr[FLID], numpy.array(Remove_Association))] = 0  # Old Way (left alot of zeros) Remove WBAREACOMI lake associations
    Lake_Link_Type_arr = Lake_Link_Type_arr[Lake_Link_Type_arr[LakeAssociation]!= 0]  # Now remove all lake associations with lake ID = 0 (added 2/14/2018)

    # Clean up and return
    print('        {0} lake associations eliminated due to multiple outlets.'.format(counter3))
    print('        {0} lakes are a headwater lake (no inflows).'.format(counter4))
    print('        {0} lakes have an endorheic lake lake condition.'.format(counter5))
    print('        {0} lakes eliminated due to having a lake immediately downstream.'.format(counter6))
    print('        {0} lakes have an outlet that drains to nowhere.'.format(counter7))
    print('        {0} lakes eliminated due to having no type 1 (outlet) segments associated with it.'.format(counter8))
    print('        Examined {0} lakes in {1:3.2f} seconds.'.format(len(seen), time.time()-tic1))
    return Lake_Link_Type_arr, problem_lakes, seen, ChainedLakes, Old_New_LakeComID, FLWBarr, Remove_Association, Tossed_Lake_Link_Type_arr

def Waterbody_SpatialJoin(in_RL, in_Lakes, link_ID_field, lake_ID_field, quiet=True):
    '''
    3/28/2023
    This function performs a simple intersection between the flowline midpoint
    geometry defined in RouteLink.nc and the waterbodies. The output is a dictionary
    of waterbodies that intersect each flowline midpoint.
    '''

    tic1 = time.time()
    print('    Starting to Intersect flowline network with waterbodies.')

    # Open the lakes file
    lakes_ds = ogr.Open(in_Lakes, 0)
    lakes_lyr = lakes_ds.GetLayerByIndex(0)
    lakes_srs = lakes_lyr.GetSpatialRef()
    lakes_LayerDef = lakes_lyr.GetLayerDefn()
    print('    Input lakes layer has {0} features.'.format(lakes_lyr.GetFeatureCount()))
    lake_fieldNames = [lakes_LayerDef.GetFieldDefn(i).GetName() for i in range(lakes_LayerDef.GetFieldCount())]
    assert lake_ID_field in lake_fieldNames

    # Open the RouteLink as points
    link_ds = ogr.Open(in_RL, 0)
    link_lyr = link_ds.GetLayerByIndex(0)
    link_srs = link_lyr.GetSpatialRef()
    link_LayerDef = link_lyr.GetLayerDefn()
    print('    Input RouteLink layer has {0} features.'.format(link_lyr.GetFeatureCount()))
    link_fieldNames = [link_LayerDef.GetFieldDefn(i).GetName() for i in range(link_LayerDef.GetFieldCount())]
    assert link_ID_field in link_fieldNames

    # Added 11/19/2020 to allow for GDAL 3.0 changes to the order of coordinates in transform
    if int(osgeo.__version__[0]) >= 3:
        # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
        lakes_srs.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        link_srs.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    # Check if a coordinate transformation (projection) must be performed
    if not lakes_srs.IsSame(link_srs):
        print('    Input shapefile projection does not match requested output. Transformatin will be applied to lakes.')
        coordTrans = osr.CoordinateTransformation(lakes_srs, link_srs)
        trans = True

    # Attempt using a grid layer returned by the gridder object
    counter = 0
    WaterbodyDict = {}
    for lake_feature in lakes_lyr:

        lake_id = lake_feature.GetField(lake_ID_field)
        counter += 1
        lake_geometry = lake_feature.GetGeometryRef()
        polygon_area = lake_geometry.GetArea()
        if not quiet:
            print('        Input Lake [{0}] area: {1}.'.format(lake_id, polygon_area))

        if trans:
            lake_geometry.Transform(coordTrans)

        # Set spatial filter to limit how many points are considered
        link_lyr.SetSpatialFilter(lake_geometry)
        if not quiet:
            print('        Input RouteLink layer has {0} features after spatial filter.'.format(link_lyr.GetFeatureCount()))

        # Find intersecting geometries
        point_counter = 0
        for link_feature in link_lyr:
            link_id = link_feature.GetField(link_ID_field)
            link_geom = link_feature.GetGeometryRef()
            if link_geom.Within(lake_geometry):
                point_counter += 1
                WaterbodyDict[link_id] = [lake_id]
            link_feature = link_geom = None
            del link_feature, link_id, link_geom
        if not quiet:
            print('        Found {0} points inside lake {1}'.format(point_counter, lake_id))

        # Remove spatial filter for next iteration
        link_lyr.SetSpatialFilter(None)
        link_lyr.ResetReading()
        lake_feature = lake_geometry = None
        del lake_id, point_counter, polygon_area

    # Clean up
    if trans:
        coordTrans = None
    lakes_ds = lakes_lyr = lakes_srs = link_ds = link_lyr = link_srs = lakes_LayerDef = link_LayerDef = None
    lakes_in_dict = sorted({x for v in WaterbodyDict.values() for x in v})
    print('    Found {0} flowlines with connectivity to {1} lakes.'.format(len(WaterbodyDict), len(set(lakes_in_dict))))
    print('    Finished intersecting flowline network with waterbodies in {0:3.2f}s'.format(time.time()-tic1))
    return WaterbodyDict

def LK_main(outDir, Flowline, Waterbody, link_ID_field, lake_ID_field, Subset_arr=None, datestr=datestr, LakeAssociation=LakeAssoc, update_RL=True, update_LK=True):
    '''
    This is the main lake pre-processing function, but written for open-source GIS pre-processing

    Flowline            The RouteLink.nc file, or None if providing a dictionary the define the mappings
    Waterbody           A Watebody feature class or a dictionary if providing pre-defined flowline:lake mappings
    outDir              A directory to write the output CSV files and other ancillary data.
    link_ID_field       The fieldname in the RouteLink.nc file to identify flowlines.
    lake_ID_field       The fieldname in the input lake shapefil to identify lakes
    Subset_arr          A numpy array with which to subset the lake list
    datestr             A string giving the current date, for file naming
    LakeAssociation     The fielname to use for ???
    '''

    # Setup Logging
    tic1 = time.time()
    print('    Lake module initiated on {0}'.format(time.ctime()))

    if Flowline is not None:
        # This is the normal case. The user wishes to evaluate flowline:waterbody associations for either a flowline feature class
        # which has the associations specified, or perform a spatial join between those flowlines and a Waterbody feature class.
        WaterbodyDict = Waterbody_SpatialJoin(Flowline, Waterbody, link_ID_field, lake_ID_field, quiet=True)
    elif isinstance(Waterbody, dict):
        # This is the case where the user is only wishing to submit a dictionary of flowline:waterbody associations for evaluation
        # Thus, choose Flowline=None and Waterbody=WaterbodyDict in the main() function arguments.
        # Populate dictionary of all flowline/lake intersections
        WaterbodyDict = {item[0]:[item[1]] for item in Waterbody.items()}       # Convert to lists

    # Prepare inputs for the Lake_Link_Type function
    rootgrp = netCDF4.Dataset(Flowline, 'r')
    FromComIDs = {link:to for link,to in zip(rootgrp.variables['link'][:],rootgrp.variables['to'][:])}
    sorted_Flowlinearr = rootgrp.variables['link'][:]
    dtypes = numpy.dtype([(FLID, 'i4'), (hydroSeq, 'i4')])           # Create a numpy dtype object
    order = numpy.empty(len(sorted_Flowlinearr), dtype=dtypes)
    order[FLID] = sorted_Flowlinearr

    # Reverse this range since the sorted array is sorted from uptream to downstream.
    order[hydroSeq] = numpy.arange(len(sorted_Flowlinearr))[::-1]
    rootgrp.close()
    del rootgrp

    # Create an array of all flowlines associated with all lakes from WaterbodyDict
    dtype = dict(names=(FLID, LakeAssociation), formats=('<i4', '<i4'))
    FLWBarr = numpy.array([(item[0], item[1][0]) for item in WaterbodyDict.items()], dtype=dtype)   # Grab the first lake association for any flowline
    print('        Found {0} unique lake ComIDs from flowline association'.format(numpy.unique(FLWBarr[LakeAssociation]).shape[0]))

    # Gather all link lake types
    Lake_Link_Type_arr, problem_lakes, seen, ChainedLakes, Old_New_LakeComID, FLWBarr, Remove_Association, Tossed_Lake_Link_Type_arr = Lake_Link_Type(FLWBarr, FromComIDs, order, subset=Subset_arr, LakeAssociation=LakeAssociation)
    unique_lakes = numpy.unique(Lake_Link_Type_arr[LakeAssociation]).shape[0]
    print('      Found {0} unique lake comID values.'.format(unique_lakes))
    print('      Found {0} outlet flowlines.'.format(Lake_Link_Type_arr[Lake_Link_Type_arr['LINK_TYPE']==1].shape[0]))
    print('      Found {0} internal flowlines.'.format(Lake_Link_Type_arr[Lake_Link_Type_arr['LINK_TYPE']==2].shape[0]))
    print('      Found {0} contributing flowlines.'.format(Lake_Link_Type_arr[Lake_Link_Type_arr['LINK_TYPE']==3].shape[0]))
    print('      Problems: {0}'.format(len(problem_lakes)))
    print('      Chained Lakes: {0}'.format(ChainedLakes.keys()))

    # Write the lake problem file to a CSV format file
    LakeProblemFile = os.path.join(outDir, 'Lake_Problems.csv')
    print('      Writing Lake_Problems dictionary to file: {0}'.format(LakeProblemFile))
    with open(LakeProblemFile,'w') as f:
        w = csv.writer(f)
        w.writerows(problem_lakes.items())

    # Write the dictionary to disk as CSV file that shows which lakes have been merged (added 1/16/2017)
    Old_New_LakeComIDFile = os.path.join(outDir, 'Old_New_LakeComIDs.csv')
    print('      Writing Lake merging dictionary to file: {0}'.format(Old_New_LakeComIDFile))
    with open(Old_New_LakeComIDFile,'w') as f:
        w = csv.writer(f)
        w.writerows(Old_New_LakeComID.items())

    # Save the Lake_Link_Type array
    if save_Lake_Link_Type_arr:
        Lake_Link_Type_File = os.path.join(outDir, 'Lake_Link_Types.csv')
        print('      Writing Lake_Link_Type array to file: {0}'.format(Lake_Link_Type_File))
        numpy.savetxt(Lake_Link_Type_File, Lake_Link_Type_arr, fmt='%i', delimiter=",")

    # Save the tossed Lake_Link_Type array
    Lake_Link_Type_File2 = os.path.join(outDir, 'Tossed_Lake_Link_Types.csv')
    print('      Writing Tossed Lake_Link_Type array to file: {0}'.format(Lake_Link_Type_File2))
    numpy.savetxt(Lake_Link_Type_File2, Tossed_Lake_Link_Type_arr, fmt='%i', delimiter=",")

    # Output to dictionary from new Lake Link Type Array (excluding inflow links)
    WaterbodyDict = {item[FLID]:item[LakeAssociation] for item in Lake_Link_Type_arr[Lake_Link_Type_arr['LINK_TYPE']<3]}

    # Now we can update the input files (LAKEPARM.nc, RouteLink.nc, streams.shp, lakes.shp)
    update_RL = True
    if update_RL and Flowline:
        print('      Updating the RouteLink file to include waterbody associations.')
        if os.path.exists(Flowline):
            in_RL = Flowline
            print('      Using provided RouteLink file: {0}'.format(in_RL))
        else:
            in_RL = os.path.join(outDir, in_RL)
            print('      Using RouteLink file from temporary directory: {0}'.format(in_RL))

        # Add lake association information into RouteLink
        rootgrp_RL = netCDF4.Dataset(in_RL, 'r+')
        variables_RL = rootgrp_RL.variables
        RL_lakes = numpy.array([WaterbodyDict.get(link, NoDataVal) for link in variables_RL['link'][:]])
        rootgrp_RL.variables['NHDWaterbodyComID'][:] = RL_lakes
        rootgrp_RL.close()
        del rootgrp_RL, variables_RL, RL_lakes

        # Add lake association to output streams layer
        streams_vector_file = os.path.join(outDir, streams_vector)
        print('      Updating the streams shapefile to include waterbody associations.')
        if os.path.exists(streams_vector_file):
            data_source = ogr.Open(streams_vector_file, 1)
            lyr = data_source.GetLayer()
            for feature in lyr:
                flowline_id = int(feature.GetField(link_ID_field))
                if flowline_id in WaterbodyDict:
                    feature.SetField("LakeID", int(WaterbodyDict.get(flowline_id, NoDataVal)))
                lyr.SetFeature(feature)
                feature = None
            data_source = lyr = None

    update_LK = True
    if update_LK:
        print('      Subsetting the LAKEPARM (if necessary) to include only valid waterbodies.')
        LakeNC = os.path.join(outDir, LK_nc)
        print('      Using LAKEPARM file from temporary directory: {0}'.format(LakeNC))

        # Subset the LAKEPARM file by building a new one with a subset of the old values
        lk_subsetList = list(set(WaterbodyDict.values()))               # List of lakes to keep in LAKEPARM
        min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals = obtain_LakeParameters(LakeNC,
                                                                                                        subsetList=lk_subsetList)
        build_LAKEPARM(LakeNC, min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals)
        del min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals

        out_lakes = out_lakes = os.path.join(outDir, LakesSHP)
        if os.path.exists(out_lakes):
            # Remove any lakes from the output feature class
            print('      Removing lakes from lakes shapefile that are not on the vector channel network')
            lake_ds = ogr.Open(out_lakes, 1)
            lake_layer = lake_ds.GetLayer()
            for feature in lake_layer:
                idval = feature.GetField(lake_ID_field)
                if idval not in lk_subsetList:
                    lake_layer.DeleteFeature(feature.GetFID())
            lake_layer.ResetReading()
            lake_ds = lake_layer = None
            del lk_subsetList, lake_ds, lake_layer

    # Clean up and return
    del Subset_arr, dtype, problem_lakes, seen, ChainedLakes, Remove_Association, FLWBarr
    print('    Finished building Lake Association and flowline connectivity tables.  Time elapsed: {0:3.2f} seconds.'.format(time.time()-tic1))
    return WaterbodyDict, Lake_Link_Type_arr, Old_New_LakeComID

def obtain_LakeParameters(in_NC, subsetList=None):
    '''
    3/22/2023

    Read an existing LAKEPARM parameter file into a set of dictionaries that can
    be sent into a function such as build_LAKEPARM in the WRF-Hydro GIS Pre-processing
    tools. This function could be handy if you wish to add lakes to an existing
    lake PARAMETER set.

    3/28/2023 - Added ability to subset an input file with a list
    '''

    # Establish an object for reading the input NetCDF file
    rootgrp = netCDF4.Dataset(in_NC, 'r')
    if LooseVersion(netCDF4.__version__) > LooseVersion('1.4.0'):
        # Change masked arrays to old default (numpy arrays always returned)
        rootgrp.set_auto_mask(False)
    ncVars = rootgrp.variables

    # Assemble dictionaries
    IDs = ncVars['lake_id'][:]

    # If no subset list is given, use all IDs
    if not subsetList:
        subsetList = IDs.tolist()

    # Must convert square kilometers in input file back to square meters before going into build_LAKEPARM function
    areas = {key:val*float(1000000) for key,val in zip(IDs, ncVars['LkArea'][:]) if key in subsetList}
    max_elevs = {key:val for key,val in zip(IDs, ncVars['LkMxE'][:]) if key in subsetList}
    OrificEs = {key:val for key,val in zip(IDs, ncVars['OrificeE'][:]) if key in subsetList}
    cen_lats = {key:val for key,val in zip(IDs, ncVars['lat'][:]) if key in subsetList}
    cen_lons = {key:val for key,val in zip(IDs, ncVars['lon'][:]) if key in subsetList}
    WeirE_vals = {key:val for key,val in zip(IDs, ncVars['WeirE'][:]) if key in subsetList}

    # min_elevs is only used for IDs, but we can reconstruct the real value using other parameters
    # Only works for lakes where WeirE or OrificeE is based on Max Elevation and Min Elevation.
    min_elevs = {key:maxE-((maxE-OrE)*(3.0/2.0)) for key,maxE,OrE in zip(IDs, ncVars['LkMxE'][:], ncVars['OrificeE'][:]) if key in subsetList}
    #min_elevs2 = {key:maxE-((maxE-WeirE)*10.) for key,maxE,WeirE in zip(IDs, ncVars['LkMxE'][:], ncVars['WeirE'][:]) if key in subsetList}

    # Clean up and return
    rootgrp.close()
    del IDs, ncVars, rootgrp
    return min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals

# --- End Functions --- #

# --- Main Codeblock --- #
if __name__ == '__main__':
    pass