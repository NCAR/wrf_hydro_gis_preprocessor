#!/usr/bin/env python3

# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# ** Copyright UCAR (c) 2015-2025
# ** University Corporation for Atmospheric Research(UCAR)
# ** National Center for Atmospheric Research(NCAR)
# ** Research Applications Laboratory(RAL)
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

'''
 Name:         WRF_Hydro_Regridding_Spatial_Weights.py
 Author:       Kevin Sampson
               Associate Scientist
               National Center for Atmospheric Research

This script is considered "low memory" because of a few key distinctions. First,
it saves the dictionary created by each node to a pickle file. Then it constructs
the output correspondence netCDF file by reading each pickle file individually,
only holding one core's worth of output at a time.


This script code is intended to generate regridding weights between a
netCDF file (ideally a GEOGRID file) or a high resolution grid in netCDF but
with a GEOGRID file to define the coordinate reference system. The process takes
the following steps:

    1) Use GEOGRID file to define the input coordinate reference system
    2) Save the NetCDF grid to in-memory raster format (optionally save to disk)
    3) Build Gridder object to speed up the spatial intersection process (must be a Cartesian grid!!)
    4) Calculate the spatial intersection between the input grid and the input polygons
    5) Export weights to netCDF.
    6) Optionally output a vector polygon mesh of the analysis grid, in GeoPackage format

2/15/2021:
    Updates to allow functionality with Python3.
    Another update involves using actual FID values rather than a range based on the
    number of features. This is because in GPKG format, it appears that the FID
    values are not always sequential, or in ascending order.
'''

# Import Python Core Modules
import sys
import os
import time
import multiprocessing
import math                                                                     # Only needed for producing coordinate system information
import pickle

# Import Additionale Modules
from osgeo import ogr
from osgeo import osr
from osgeo import gdal
import netCDF4
import numpy
from osgeo import gdalconst

# Module settings
tic1 = time.time()
sys.dont_write_bytecode = True
gdal.UseExceptions()                                                            # this allows GDAL to throw Python Exceptions
gdal.PushErrorHandler('CPLQuietErrorHandler')
multiprocessing.freeze_support()

# ---------- Global variables ---------- #

# Run configurations
ticker = 1000                                                                   # Ticker for how often to print(a message
SaveRaster = True                                                              # Save the NetCDF grid as a GeoTiff file
OutputGridPolys = False                                                         # Save output vector grid polygons
check_geometry = True                                                           # Check geometry to trap any .Area issues
threshold = 50000000000                                                         # Number of square meters above which a polygon gets split up
splits = 10                                                                     # Number of chunks to divide the large polygon into
NC_format = 'NETCDF4'                                                           # NetCDF output format. Others: 'NETCDF4_CLASSIC'

# Multiprocessing parameters
CPU_Count = multiprocessing.cpu_count()                                         # Find the number of CPUs to chunk the data into
Processors = 6                                                                  # To spread the processing over a set number of cores, otherwise use Processors = CPU_Count-1
#Processors = CPU_Count-1
pool_size = Processors*100                                                      # Having a large number of worker processes can help with overly slow cores

# Output Directory
OutDir = r'C:\OutputDirectory'  # Output Directory

# Input netCDF GEOGRID file
isGEOGRID = True                                                               # Is the desired domain a GEOGRID file?
in_geogrid = r"C:\geo_em.d01.nc"
geogridVariable = 'HGT_M'

# Input vector layer (basins, etc)
inDriverName = 'ESRI Shapefile'                                                 # 'GPKG' / 'ESRI Shapefile'
in_basins = r"C:\Basins.shp"
layerName = 'Basins'
fieldname = 'GRIDCODE'   # Unique identifier field

# Output polygon vector file
if OutputGridPolys == True:
    outDriverName = 'ESRI Shapefile'                                                      # 'GPKG' / 'ESRI Shapefile'
    OutGridFile = os.path.join(OutDir, 'grid_polygons.shp')                               # Output GeoPackage

# Output Weights file
regridweightnc = os.path.join(OutDir, 'spatialweights_out.nc')

# If desired netCDF domain is not a GEOGRID file, and is a routing grid
use_inRaster = True
if use_inRaster:
    in_nc = r"C:\DEM.tif"
    Variable = 'HGT_M'
    nest = 1
elif isGEOGRID:
    # This is for input GEOGRID files
    in_nc = in_geogrid
    Variable = geogridVariable
    nest = 1                                                                    # Ratio of high-resolution grid cells to GEOGRID resolution (nest:1 for original)

# Output GeoTiff of the geogrid file
if SaveRaster == True:
    OutGTiff = os.path.join(OutDir, 'Model_grid.tif')

projdict = {1: 'Lambert Conformal Conic',
            2: 'Polar Stereographic',
            3: 'Mercator',
            6: 'Cylindrical Equidistant'}

# Global attributes for altering the sphere radius used in computations. Do not alter sphere_radius for standard WRF simulations
sphere_radius = 6370000.0                                                       # Radius of sphere to use (WRF Default = 6370000.0m)

# Switch for activating old style (used in NWM v1.0, v1.1, and v1.2) of reading CRS from Geogrid file
oldCRS = False
# ---------- End Global variables ---------- #

# ---------- Classes ---------- #
class Gridder_Layer(object):
    '''Class with which to create the grid intersecting grid cells based on a feature
    geometry envelope. Provide grid information to initiate the class, and use getgrid()
    to generate a grid mesh and index information about the intersecting cells.

    Note:  The i,j index begins with (1,1) in the Lower Left corner.'''
    def __init__(self, DX, DY, x00, y00, nrows, ncols):
        self.DX = DX
        self.DY = DY
        self.x00 = x00
        self.y00 = y00
        self.nrows = nrows
        self.ncols = ncols

    def getgrid(self, envelope, layer):
        """Gridder.getgrid() takes as input an OGR geometry envelope, and will
        compute the grid polygons that intersect the evelope, returning a list
        of grid cell polygons along with other attribute information."""
        # Calculate the number of grid cells necessary
        xmin, xmax, ymin, ymax = envelope

        # Find the i and j indices
        i0 = int((xmin-self.x00)/self.DX // 1)                                  # Floor the value
        j0 = int(abs((ymax-self.y00)/self.DY) // 1)                             # Floor the absolute value
        i1 = int((xmax-self.x00)/self.DX // 1)                                  # Floor the value
        j1 = int(abs((ymin-self.y00)/self.DY) // 1)                             # Floor the absolute value

        # Create a new field on a layer. Add one attribute
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('i_index', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('j_index', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('cellsize', ogr.OFTReal))
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
                feature.SetField('j_index', self.nrows-y)
                feature.SetGeometry(geometry)                                      # Make a feature from geometry object
                layer.CreateFeature(feature)
                geometry = feature = None
                del x0, x1, y1, y0, id1
        return layer
# ---------- End Classes ---------- #


# ---------- Functions ---------- #
def checkfield(layer, fieldname, string1):
    '''Check for existence of provided fieldnames'''
    layerDefinition = layer.GetLayerDefn()
    fieldslist = []
    for i in range(layerDefinition.GetFieldCount()):
        fieldslist.append(layerDefinition.GetFieldDefn(i).GetName())
    if fieldname in fieldslist:
        i = fieldslist.index(fieldname)
        field_defn = layerDefinition.GetFieldDefn(i)
    else:
        print('    Field %s not found in input %s. Terminating...' %(fieldname, string1))
        raise SystemExit
    return field_defn, fieldslist


def getfieldinfo(field_defn, fieldname):
    '''Get information about field type for buildng the output NetCDF file later'''
    if field_defn.GetType() == ogr.OFTInteger:
        fieldtype = 'integer'
        print("found ID type of Integer")
    elif field_defn.GetType() == ogr.OFTInteger64:
        fieldtype = 'integer64'
        print("found ID type of Integer64")
    elif field_defn.GetType() == ogr.OFTReal:
        fieldtype = 'float'
        print("found ID type of Float")
        #print("field type: OFTReal not currently supported in output NetCDF file.")
        #raise SystemExit
    elif field_defn.GetType() == ogr.OFTString:
        fieldtype = 'string'
        print("found ID type of String")
    else:
        print("ID Type not found ... Exiting")
        raise SystemExit
    print("    Field Type for field '%s': %s (%s)" %(fieldname, field_defn.GetType(), fieldtype))
    return fieldtype


def loadpickle(fp):
    with open(fp, 'rb') as fh:
        listOfObj = pickle.load(fh)
    return listOfObj


def Read_GEOGRID_for_SRS_old(in_nc, Variable):
    '''Read NetCDF GEOGRID file as a GDAL raster object. Much of the code below was borrowed
    from https://github.com/rveciana/geoexamples/blob/master/python/wrf-NetCDF/read_netcdf.py'''

    tic = time.time()
    try:
        print('Input netCDF GEOGRID file: %s    Variable: %s' %(in_nc, Variable))
        ds_in = gdal.Open(in_nc, gdalconst.GA_ReadOnly)                         # Open input netCDF file using GDAL
        subdatasets = ds_in.GetSubDatasets()                                    # Gather subdatasets from input netCDF file
        variables = [subdataset[1].split(" ")[1] for subdataset in subdatasets] # Gather variables in the input netCDF file
        print('Variables found in input NC file: %s' %(variables))
        if Variable in variables:
            src_ds = gdal.Open('NETCDF:"'+in_nc+'":%s' %(Variable))             # Open using NETCDF driver, file name, and variable
            metadata = ds_in.GetMetadata()                                      # Read metadata
            srcband = src_ds.GetRasterBand(1)                                   # Get raster band
            ncvar = srcband.ReadAsArray()                                       # Read variable as a numpy array
            ds_in = subdatasets = None

            # Initiate dictionaries of GEOGRID projections and parameters
            projdict = {1: 'Lambert Conformal Conic', 2: 'Polar Stereographic', 3: 'Mercator', 6: 'Cylindrical Equidistant'}

            # Read metadata for grid information
            map_pro = int(metadata['NC_GLOBAL#MAP_PROJ'])
            DX = float(metadata['NC_GLOBAL#DX'])
            DY = float(metadata['NC_GLOBAL#DY'])
            corner_lats = metadata['NC_GLOBAL#corner_lats']
            corner_lons = metadata['NC_GLOBAL#corner_lons']

            # Gather corner information [order = (LL center, ULcenter, URcenter, LRcenter, LLLedge, ULLedge, URRedge, LRRedge, LLBedge, ULUedge, URUedge, LRBedge, LLcorner, ULcorner, URcorner, LRcorner)
            corner_latslist = corner_lats.strip('{}').split(',')                # Create list of strings from the corner_lats attribute
            corner_lonslist = corner_lons.strip('{}').split(',')                # Create list of strings from the corner_lons attribute
            ##LLcenter = [corner_lonslist[0], corner_latslist[0]]               # Lower left of the Mass grid staggering
            #ULcenter = [corner_lonslist[1], corner_latslist[1]]               # Upper left of the Mass grid staggering
            ##URcenter = [corner_lonslist[2], corner_latslist[2]]               # Upper right of the Mass grid staggering
            ##LRcenter = [corner_lonslist[3], corner_latslist[3]]               # Lower right of the Mass grid staggering
            ##LLcorner = [corner_lonslist[12], corner_latslist[12]]             # Lower left of the Unstaggered grid
            ULcorner = [corner_lonslist[13], corner_latslist[13]]               # Upper left of the Unstaggered grid
            ##URcorner = [corner_lonslist[14], corner_latslist[14]]             # Upper right of the Unstaggered grid
            ##LRcorner = [corner_lonslist[15], corner_latslist[15]]             # Lower right of the Unstaggered grid

            # Pick a corner or center point to 'hang' the raster from
            lon = float(ULcorner[0])                                            #lon = float(ULcenter[0])
            lat = float(ULcorner[1])                                            #lat = float(ULcenter[1])

            # Read metadata for projection parameters
            if 'NC_GLOBAL#TRUELAT1' in metadata.keys():
                standard_parallel_1 = float(metadata['NC_GLOBAL#TRUELAT1'])
            if 'NC_GLOBAL#TRUELAT2' in metadata.keys():
                standard_parallel_2 = float(metadata['NC_GLOBAL#TRUELAT2'])
            if 'NC_GLOBAL#STAND_LON' in metadata.keys():
                central_meridian = float(metadata['NC_GLOBAL#STAND_LON'])
            if 'NC_GLOBAL#POLE_LAT' in metadata.keys():
                pole_latitude = float(metadata['NC_GLOBAL#POLE_LAT'])
            if 'NC_GLOBAL#POLE_LON' in metadata.keys():
                pole_longitude = float(metadata['NC_GLOBAL#POLE_LON'])
            if 'NC_GLOBAL#CEN_LAT' in metadata.keys():
                latitude_of_origin = float(metadata['NC_GLOBAL#CEN_LAT'])

            # Initiate OSR spatial reference object - See http://gdal.org/java/org/gdal/osr/SpatialReference.html
            proj1 = osr.SpatialReference()

            # ---- NOTE: Below is experimental & untested ---- #

            # Use projection information from global attributes to populate OSR spatial reference object
            # See this website for more information on defining coordinate systems: http://gdal.org/java/org/gdal/osr/SpatialReference.html
            print('    Map Projection: %s' %projdict[int(metadata['NC_GLOBAL#MAP_PROJ'])])
            if map_pro == 1:
                # Lambert Conformal Conic
                if 'standard_parallel_2' in locals():
                    proj1.SetLCC(standard_parallel_1, standard_parallel_2, latitude_of_origin, central_meridian, 0, 0)
                    #proj1.SetLCC(double stdp1, double stdp2, double clat, double clong, double fe, double fn)        # fe = False Easting, fn = False Northing
                else:
                    proj1.SetLCC1SP(latitude_of_origin, central_meridian, 1, 0, 0)       # Scale = 1???
                    #proj1.SetLCC1SP(double clat, double clong, double scale, double fe, double fn)       # 1 standard parallell
            elif map_pro == 2:
                # Polar Stereographic
                phi1 = float(standard_parallel_1)                               # Set up pole latitude
                ### Back out the central_scale_factor (minimum scale factor?) using formula below using Snyder 1987 p.157 (USGS Paper 1395)
                ##phi = math.copysign(float(pole_latitude), float(latitude_of_origin))    # Get the sign right for the pole using sign of CEN_LAT (latitude_of_origin)
                ##central_scale_factor = (1 + (math.sin(math.radians(phi1))*math.sin(math.radians(phi))) + (math.cos(math.radians(float(phi1)))*math.cos(math.radians(phi))))/2
                # Method where central scale factor is k0, Derivation from C. Rollins 2011, equation 1: http://earth-info.nga.mil/GandG/coordsys/polar_stereographic/Polar_Stereo_phi1_from_k0_memo.pdf
                # Using Rollins 2011 to perform central scale factor calculations. For a sphere, the equation collapses to be much  more compact (e=0, k90=1)
                central_scale_factor = (1 + math.sin(math.radians(abs(phi1))))/2        # Equation for k0, assumes k90 = 1, e=0. This is a sphere, so no flattening
                print('        Central Scale Factor: %s' %central_scale_factor)
                #proj1.SetPS(latitude_of_origin, central_meridian, central_scale_factor, 0, 0)    # example: proj1.SetPS(90, -1.5, 1, 0, 0)
                # Adjusted 8/7/2017 based on changes made 4/4/2017 as a result of Monaghan's polar sterographic domain. Example: proj1.SetPS(90, -1.5, 1, 0, 0)
                proj1.SetPS(pole_latitude, central_meridian, central_scale_factor, 0, 0)
                #proj1.SetPS(double clat, double clong, double scale, double fe, double fn)
            elif map_pro == 3:
                # Mercator Projection
                proj1.SetMercator(latitude_of_origin, central_meridian, 1, 0, 0)     # Scale = 1???
                #proj1.SetMercator(double clat, double clong, double scale, double fe, double fn)
            elif map_pro == 6:
                # Cylindrical Equidistant (or Rotated Pole)
                if pole_latitude != float(90) or pole_longitude != float(0):
                    # if pole_latitude, pole_longitude, or stand_lon are changed from thier default values, the pole is 'rotated'.
                    print('[PROBLEM!] Cylindrical Equidistant projection with a rotated pole is not currently supported.')
                    raise SystemExit
                else:
                    proj1.SetEquirectangular(latitude_of_origin, central_meridian, 0, 0)
                    #proj1.SetEquirectangular(double clat, double clong, double fe, double fn)
                    #proj1.SetEquirectangular2(double clat, double clong, double pseudostdparallellat, double fe, double fn)

            # Set Geographic Coordinate system (datum) for projection
            proj1.SetGeogCS('WRF-Sphere', 'Sphere', '', 6370000.0, 0.0)      # Could try 104128 (EMEP Sphere) well-known?
            #proj1.SetGeogCS(String pszGeogName, String pszDatumName, String pszSpheroidName, double dfSemiMajor, double dfInvFlattening)

            variables = src_ds = ncvar = metadata = srcband = ncvar = None
        else:
            print('Could not find variable: %s in file: %s' %(Variable, in_nc))

        # Set the origin for the output raster (in GDAL, usuall upper left corner) using projected corner coordinates
        wgs84_proj = osr.SpatialReference()
        wgs84_proj.ImportFromEPSG(4326)
        transform = osr.CoordinateTransformation(wgs84_proj, proj1)
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint_2D(lon, lat)
        point.Transform(transform)
        x00 = point.GetX(0)
        y00 = point.GetY(0)
        #x00 = point.GetX(0)-(DX/2)                                              # This mimicks the Arcpy method if you use the ULcenter[0]-(DX/2)
        #y00 = point.GetY(0)+abs(DY/2)                                           # This mimicks the Arcpy method if you use the ULcenter[1]+abs(DY/2)

    except RuntimeError as e:
        print('Unable to open %s' %in_nc)
        print(e)
        raise SystemExit
    print('Created projection definition from input NetCDF GEOGRID file %s in %.2fs.' %(in_nc, time.time()-tic))

    # Clear objects and return
    return proj1, DX, DY, x00, y00


def Read_GEOGRID_for_SRS(in_nc):
    """
    The input NetCDF file (Ideally WRF GEOGRID) gets georeferenced and projection
    infromation is created.

    10/6/2017: Proj4 string generation was added, with definitions adapted from
    https://github.com/NCAR/wrf-python/blob/develop/src/wrf/projection.py
    """

    tic = time.time()

    # First step: Import and georeference NetCDF file
    print('  Step 1: NetCDF Conversion initiated...')
    print('    Input netCDF GEOGRID file: %s' %in_nc)

    # Read input WPS GEOGRID file
    # Loop through global variables in NetCDF file to gather projection information
    rootgrp = netCDF4.Dataset(in_nc, 'r')                                       # Establish an object for reading the input NetCDF file
    globalAtts = rootgrp.__dict__                                               # Read all global attributes into a dictionary
    map_pro = globalAtts['MAP_PROJ']                                            # Find out which projection this GEOGRID file is in
    print('    Map Projection: %s' %projdict[map_pro])

    # Collect grid corner XY and DX DY for creating ascii raster later
    if 'corner_lats' in globalAtts:
        corner_lat = globalAtts['corner_lats'][13].astype(numpy.float64)        # Note: The values returned are corner points of the mass grid. 13 = Upper left of the Unstaggered grid
    if 'corner_lons' in globalAtts:
        corner_lon = globalAtts['corner_lons'][13].astype(numpy.float64)        # Note: The values returned are corner points of the mass grid. 13 = Upper left of the Unstaggered grid
    if 'DX' in globalAtts:
        DX = globalAtts['DX'].astype(numpy.float32)
    if 'DY' in globalAtts:
        DY = globalAtts['DY'].astype(numpy.float32)

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
        latitude_of_origin = globalAtts['MOAD_CEN_LAT'].astype(numpy.float64)         # Added 2/26/2017 by KMS
    elif 'CEN_LAT' in globalAtts:
        print('    Using CEN_LAT for latitude of origin.')
        latitude_of_origin = globalAtts['CEN_LAT'].astype(numpy.float64)
    del globalAtts
    rootgrp.close()

    # Initiate OSR spatial reference object - See http://gdal.org/java/org/gdal/osr/SpatialReference.html
    proj1 = osr.SpatialReference()

    # Use projection information from global attributes to populate OSR spatial reference object
    # See this website for more information on defining coordinate systems: http://gdal.org/java/org/gdal/osr/SpatialReference.html
    print('    Map Projection: %s' %projdict[int(map_pro)])
    if map_pro == 1:
        # Lambert Conformal Conic
        if 'standard_parallel_2' in locals():
            proj1.SetLCC(standard_parallel_1, standard_parallel_2, latitude_of_origin, central_meridian, 0, 0)
            #proj1.SetLCC(double stdp1, double stdp2, double clat, double clong, double fe, double fn)        # fe = False Easting, fn = False Northing
        else:
            proj1.SetLCC1SP(latitude_of_origin, central_meridian, 1, 0, 0)       # Scale = 1???
            #proj1.SetLCC1SP(double clat, double clong, double scale, double fe, double fn)       # 1 standard parallell
    elif map_pro == 2:
        # Polar Stereographic
        phi1 = standard_parallel_1
        ### Back out the central_scale_factor (minimum scale factor?) using formula below using Snyder 1987 p.157 (USGS Paper 1395)
        ##phi = math.copysign(float(pole_latitude), float(latitude_of_origin))    # Get the sign right for the pole using sign of CEN_LAT (latitude_of_origin)
        ##central_scale_factor = (1 + (math.sin(math.radians(phi1))*math.sin(math.radians(phi))) + (math.cos(math.radians(float(phi1)))*math.cos(math.radians(phi))))/2
        # Method where central scale factor is k0, Derivation from C. Rollins 2011, equation 1: http://earth-info.nga.mil/GandG/coordsys/polar_stereographic/Polar_Stereo_phi1_from_k0_memo.pdf
        # Using Rollins 2011 to perform central scale factor calculations. For a sphere, the equation collapses to be much  more compact (e=0, k90=1)
        central_scale_factor = (1 + math.sin(math.radians(abs(phi1))))/2        # Equation for k0, assumes k90 = 1, e=0. This is a sphere, so no flattening
        print('        Central Scale Factor: %s' %central_scale_factor)
        #proj1.SetPS(latitude_of_origin, central_meridian, central_scale_factor, 0, 0)    # example: proj1.SetPS(90, -1.5, 1, 0, 0)
        proj1.SetPS(pole_latitude, central_meridian, central_scale_factor, 0, 0)    # Adjusted 8/7/2017 based on changes made 4/4/2017 as a result of Monaghan's polar sterographic domain. Example: proj1.SetPS(90, -1.5, 1, 0, 0)
        #proj1.SetPS(double clat, double clong, double scale, double fe, double fn)
    elif map_pro == 3:
        # Mercator Projection
        proj1.SetMercator(latitude_of_origin, central_meridian, 1, 0, 0)     # Scale = 1???
        #proj1.SetMercator(double clat, double clong, double scale, double fe, double fn)
    elif map_pro == 6:
        # Cylindrical Equidistant (or Rotated Pole)
        if pole_latitude != float(90) or pole_longitude != float(0):
            # if pole_latitude, pole_longitude, or stand_lon are changed from thier default values, the pole is 'rotated'.
            print('[PROBLEM!] Cylindrical Equidistant projection with a rotated pole is not currently supported.')
            raise SystemExit
        else:
            proj1.SetEquirectangular(latitude_of_origin, central_meridian, 0, 0)
            #proj1.SetEquirectangular(double clat, double clong, double fe, double fn)
            #proj1.SetEquirectangular2(double clat, double clong, double pseudostdparallellat, double fe, double fn)

    # Set Geographic Coordinate system (datum) for projection
    proj1.SetGeogCS('WRF_Sphere', 'Sphere', '', sphere_radius, 0.0)      # Could try 104128 (EMEP Sphere) well-known?
    #proj1.SetGeogCS(String pszGeogName, String pszDatumName, String pszSpheroidName, double dfSemiMajor, double dfInvFlattening)

    # Set the origin for the output raster (in GDAL, usuall upper left corner) using projected corner coordinates
    wgs84_proj = osr.SpatialReference()
    wgs84_proj.ImportFromEPSG(4326)
    transform = osr.CoordinateTransformation(wgs84_proj, proj1)
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint_2D(corner_lon, corner_lat)
    point.Transform(transform)
    x00 = point.GetX(0)
    y00 = point.GetY(0)
    print('  Created projection definition from input NetCDF GEOGRID file %s in %.2fs.' %(in_nc, time.time()-tic))
    return proj1, DX, DY, x00, y00


def NetCDF_to_Raster(in_nc, Variable, proj_in=None, DX=1, DY=-1, x00=0, y00=0):
    '''This funciton takes in an input netCDF file, a variable name, the ouput
    raster name, and the projection definition and writes the grid to the output
    raster. This is useful, for example, if you have a FullDom netCDF file and
    the GEOGRID that defines the domain. You can output any of the FullDom variables
    to raster.'''

    tic = time.time()
    try:
        #print('        in_nc: %s    Variable: %s' %(in_nc, Variable))
        ds_in = gdal.Open(in_nc, gdalconst.GA_ReadOnly)                         # Open input netCDF file using GDAL
        subdatasets = ds_in.GetSubDatasets()                                    # Gather subdatasets from input netCDF file
        variables = [subdataset[1].split(" ")[1] for subdataset in subdatasets] # Gather variables in the input netCDF file
        ds_in = subdatasets = None
        if Variable in variables:
            src_ds = gdal.Open('NETCDF:"'+in_nc+'":%s' %(Variable))             # Open using NETCDF driver, file name, and variable
            srcband = src_ds.GetRasterBand(1)                                   # Get raster band
            ncvar = srcband.ReadAsArray()                                       # Read variable as a numpy array

            # Set up driver for GeoTiff output
            driver = gdal.GetDriverByName('Mem')                                # Write to Memory
            if driver is None:
                print('    %s driver not available.' % 'Memory')

            # Set up the dataset and define projection/raster info
            DataSet = driver.Create('', srcband.XSize, srcband.YSize, 1, gdal.GDT_Float32)         # the '1' is for band 1.
            if proj_in is not None:
                DataSet.SetProjection(proj_in.ExportToWkt())
            if DX==1 and DY==-1:
                DataSet.SetGeoTransform((0, DX, 0, 0, 0, DY))                   # Default (top left x, w-e resolution, 0=North up, top left y, 0 = North up, n-s pixel resolution (negative value))
            else:
                DataSet.SetGeoTransform((x00, DX, 0, y00, 0, -DY))              # (top left x, w-e resolution, 0=North up, top left y, 0 = North up, n-s pixel resolution (negative value))
            DataSet.GetRasterBand(1).WriteArray(ncvar)                          # Write the array
            stats = DataSet.GetRasterBand(1).GetStatistics(0,1)                 # Calculate statistics
            #src_ds = srcband = None
            ncvar = driver = None

    except RuntimeError as e:
        print('Unable to open %s' %in_nc)
        print(e)
        raise SystemExit

    # Clear objects and return
    print('Created raster in-memory from input NetCDF file %s in %.2fs.' %(in_nc, time.time()-tic))
    return DataSet


def create_polygons_from_info(gridder_obj, proj1, outputFile, outDriverName, ticker=10000):
    '''This will take the grid info index created in Read_GEOGRID_for_SRS() and
    NetCDF_to_Raster() and produce a GeoPackage of the grid cells.'''

    print('Starting to produce polygon vector file from index')
    tic = time.time()

    # Check if files exist and delete
    if os.path.isfile(outputFile)==True:
        os.remove(outputFile)
        print('      Removed existing file: %s' %outputFile)

    # Now convert it to a vector file with OGR
    drv = ogr.GetDriverByName(outDriverName)
    if drv is None:
        print('      %s driver not available.' % outDriverName)
    else:
        print( '      %s driver is available.' % outDriverName)
        driver = ogr.GetDriverByName(outDriverName)
        datasource = driver.CreateDataSource(outputFile)
    if datasource is None:
        print('      Creation of output file failed.\n')
        raise SystemExit

    # Create output polygon vector file
    layer = datasource.CreateLayer('gridpolys', srs=proj1, geom_type=ogr.wkbPolygon)
    if layer is None:
        print('        Layer creation failed.\n')
        raise SystemExit
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('i_index', ogr.OFTInteger))
    layer.CreateField(ogr.FieldDefn('j_index', ogr.OFTInteger))
    LayerDef = layer.GetLayerDefn()

    point_ref=ogr.osr.SpatialReference()
    point_ref.ImportFromEPSG(4326)                                              # WGS84
    #coordTrans2 = ogr.osr.CoordinateTransformation(proj1, point_ref)            # Create transformation for converting to WGS84

    # Pull info out of th gridder object in order to create a polygon
    ncols = gridder_obj.ncols
    nrows = gridder_obj.nrows
    x00 = gridder_obj.x00
    y00 = gridder_obj.y00
    DX = gridder_obj.DX
    DY = gridder_obj.DY

    # Create polygon object that is fully inside the outer edge of the domain
    myRing = ogr.Geometry(type=ogr.wkbLinearRing)
    myRing.AddPoint(x00+(DX/2), y00-(DY/2))
    myRing.AddPoint(x00+(ncols*DX)-(DX/2), y00-(DY/2))
    myRing.AddPoint(x00+(ncols*DX)-(DX/2), y00-(nrows*DY)+(DY/2))
    myRing.AddPoint(x00+(DX/2), y00-(nrows*DY)+(DY/2))
    myRing.AddPoint(x00+(DX/2), y00-(DY/2))                                     #close ring
    geometry = ogr.Geometry(type=ogr.wkbPolygon)
    geometry.AssignSpatialReference(proj1)
    geometry.AddGeometry(myRing)

    tic2 = time.time()
    layer = gridder_obj.getgrid(geometry.GetEnvelope(), layer)
    print('      %s polygons returned in %s seconds.' %(layer.GetFeatureCount(), time.time()-tic2))
    myRing = geometry = None

    print('Done producing output vector polygon shapefile in %ss' %(time.time()-tic))
    datasource = layer = None


def split_vertical(polygon, peices=2):
    '''Creates a specified number of clipping geometries which are boxes used to
    clip an OGR feature. Returns a list of geometry objects which are verticaly
    split chunks of the original polygon.'''

    tic = time.time()

    # Get polygon geometry information
    polygeom = polygon.GetGeometryRef()
    polygeom.CloseRings()                                                       # Ensure all rings are closed

    # Get min/max
    xmin, xmax, ymin, ymax = polygeom.GetEnvelope()                             # Get individual bounds from bounding envelope
    horizontal_dist = xmax - xmin                                               # Distance across the horizontal plane

    # Create clipping geometries
    clippolys = []           # List of new polygons
    interval = horizontal_dist/peices                                           # Split the horizontal distance using numsplits
    for split in range(peices):

        # Create clip-box bounds
        x0 = xmin+(split*interval)                                              # X-min - changes with section
        x1 = xmin+((split+1)*interval)                                          # X-max - changes with section
        y0 = ymin                                                               # Y-min - always the same
        y1 = ymax                                                               # Y-max - always the same

        # Create geometry for clip box
        myRing = ogr.Geometry(type=ogr.wkbLinearRing)
        myRing.AddPoint(x0, y1)
        myRing.AddPoint(x1, y1)
        myRing.AddPoint(x1, y0)
        myRing.AddPoint(x0, y0)
        myRing.AddPoint(x0, y1)                                                 #close ring
        geometry = ogr.Geometry(type=ogr.wkbPolygon)
        geometry.AddGeometry(myRing)

        # Add to the list of clipping geometries to be returned
        clippolys.append(geometry)
    #print('         Polygon envelope split into %s sections in %.2fs seconds' %(peices, (time.time()-tic)))
    return clippolys


def perform_intersection(gridder_obj, proj1, layer, fieldname, ticker=10000, corenum=0):
    '''This function performs the intersection between two geometries.'''

    # Counter initiate
    counter = 0
    counter2 = 0

    # Test intersection with layer
    tic2 = time.time()
    spatialweights = {}                                                         # This yields the fraction of the key polygon that each overlapping polygon contributes
    regridweights = {}                                                          # This yields the fraction of each overlapping polygon that intersects the key polygon - for regridding
    other_attributes = {}                                                       # This dicitonary stores the i,j indices of the grid cells
    allweights = {}                                                             # This dictionary will store the final returned data

    # Set up coordinate transform from Layer1 to WGS84 for lat/lon
    proj2 = layer.GetSpatialRef()
    coordTrans = osr.CoordinateTransformation(proj2, proj1)                     # Coordinate tansformation from layer to layer1

    # Attempt using a grid layer returned by the gridder object
    #print('[{0: >3}]    Layer feature count: {1: ^8}'.format(corenum, layer.GetFeatureCount()))
    for feature2 in layer:
        id2 = feature2.GetField(fieldname)
        counter2 += 1
        geometry2 = feature2.GetGeometryRef()
        geometry2.Transform(coordTrans)                                         # Transform the geometry from layer CRS to layer1 CRS
        polygon_area = geometry2.GetArea()

        # Check to find incompatible geometry types
        if check_geometry:
            if geometry2==None:
                print('  polygon {0} is NoneType'.format(id2))
                continue
            if not geometry2.IsValid():
                #print('polygon {0} geometry invalid.'.format(id2))
                #geometry2 = geometry2.Union(geometry2)
                geometry2 = geometry2.Buffer(0)
                if not geometry2:
                    print('  polygon {0} is NoneType after Union operation.'.format(id2))
                    continue
                if not geometry2.IsValid():
                    print('  polygon {0} geometry not fixed by performing self-union or Buffer.'.format(id2))
                    continue

        # Split into parts
        if polygon_area > threshold:
            print('[{0: >3}]    Polygon: {1: ^8} area = {2: ^12}. Splitting into {3: ^2} sections.'.format(corenum, id2, polygon_area, splits))
            Areas = []
            inters = 0
            clip_polys = split_vertical(feature2, splits)

            # Create temporary output polygon vector file to store the input feature
            drv1 = ogr.GetDriverByName('Memory')
            in_ds = drv1.CreateDataSource('in_ds')
            inlayer = in_ds.CreateLayer('in_ds', srs=proj1, geom_type=ogr.wkbPolygon)
            LayerDef = inlayer.GetLayerDefn()                                    # Fetch the schema information for this layer
            infeature = ogr.Feature(LayerDef)                                     # Create a new feature (attribute and geometry)
            infeature.SetGeometry(geometry2)                                       # Make a feature from geometry object
            inlayer.CreateFeature(infeature)

            for num,clipgeom in enumerate(clip_polys):
                tic3 = time.time()

                # Create temporary output polygon vector file
                out_ds = drv1.CreateDataSource('out_ds')
                outlayer = out_ds.CreateLayer('out_ds', srs=proj1, geom_type=ogr.wkbPolygon)
                outlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))   # Create a new field on a layer. Add one attribute

                # Create temporary in-memory feature layer to store clipping geometry
                clip_ds = drv1.CreateDataSource('clip_ds')
                cliplayer = clip_ds.CreateLayer('clip_ds', srs=proj1, geom_type=ogr.wkbPolygon)
                LayerDef2 = cliplayer.GetLayerDefn()                                    # Fetch the schema information for this layer
                clipfeature = ogr.Feature(LayerDef2)                                     # Create a new feature (attribute and geometry)
                clipfeature.SetGeometry(clipgeom)                                       # Make a feature from geometry object
                cliplayer.CreateFeature(clipfeature)

                # Perform clip
                inlayer.Clip(cliplayer, outlayer)

                # Read clipped polygon feature
                assert outlayer.GetFeatureCount() == 1                      # Make sure the output has only 1 feature
                feat = outlayer.GetNextFeature()                            # The output should only have one feature
                geometry3 = feat.GetGeometryRef()

                # Create a Layer so that the SetSpatialFilter method can be used (faster for very large geometry2 polygons)
                drv = ogr.GetDriverByName('Memory')
                dst_ds = drv.CreateDataSource('out')
                gridlayer = dst_ds.CreateLayer('out_%s' %corenum, srs=proj1, geom_type=ogr.wkbPolygon)
                gridlayer = gridder_obj.getgrid(geometry3.GetEnvelope(), gridlayer)     # Generate the grid layer
                gridlayer.SetSpatialFilter(geometry3)                                   # Use the SetSpatialFilter method to thin the layer's geometry

                # Check to find incompatible geometry types
                #if check_geometry:
                #    keep_Going = True
                #    copyLayer = False
                #    doSQL = False
                #    for item in gridlayer:
                #        itemgeom = item.geometry()
                #        if itemgeom.Intersection(geometry3) is None:
                #            print('grid polygon %s intersection geometry invalid against polygon %s' %(item.GetField(0), id2))
                #            gridlayer.DeleteFeature(item.GetFID())
                #            #keep_Going = False
                #            copyLayer = True
                #            #doSQL = True
                #            continue
                #    if doSQL:
                #        #gridlayer = dst_ds.ExecuteSQL("select * from out_%s WHERE ST_IsValid(geometry)" %corenum, dialect = "SQLITE")  # Only select valid geometry
                #        gridlayer = dst_ds.ExecuteSQL("select ST_Buffer(geometry, 0), * from out_%s" %corenum, dialect = "SQLITE") # Buffer with 0 distance to create valid geometry
                #    if copyLayer:
                #        driver = ogr.GetDriverByName('ESRI Shapefile')
                #        gridds = driver.CreateDataSource(os.path.join(OutDir, '%s.shp' %id2))
                #        griddslayer = gridds.CopyLayer(gridlayer, 'gridlayer')
                #        gridds = gridslayer = driver = None
                #    if not keep_Going:
                #        continue
                #    gridlayer.ResetReading()

                # First find all intersection areas
                Areas += [[item.GetField(0), item.geometry().Intersection(geometry3).Area(), item.geometry().Area(), item.GetField(1), item.GetField(2)] for item in gridlayer]  # Only iterate over union once
                inters += len(Areas)
                counter += inters                                                       # Advance the counter

                #flush memory
                clipfeature = clipgeom = geometry3 = feat = outlayer = cliplayer = clipfeature = None  # destroy these
                print('[{0: >3}]          Chunk {1: ^2}. Time elapsed: {2: ^4.2f} seconds.'.format(corenum, num+1, (time.time()-tic3)))

            # Collapse all duplicates back down to 1 list
            AreaDict = {}
            for item in Areas:
                try:
                    AreaDict[item[0]][1] += item[1]
                except KeyError:
                    AreaDict[item[0]] = item
            Areas = AreaDict.values()

        else:
            '''This is the normal case where polygons are smaller or more uniform in size.'''

            # Create a Layer so that the SetSpatialFilter method can be used (faster for very large geometry2 polygons)
            drv = ogr.GetDriverByName('Memory')
            dst_ds = drv.CreateDataSource('out')
            gridlayer = dst_ds.CreateLayer('out_%s' %corenum, srs=proj1, geom_type=ogr.wkbPolygon)
            gridlayer = gridder_obj.getgrid(geometry2.GetEnvelope(), gridlayer)     # Generate the grid layer
            gridlayer.SetSpatialFilter(geometry2)                                   # Use the SetSpatialFilter method to thin the layer's geometry

            ## Check to find incompatible geometry types
            #if check_geometry:
            #    keep_Going = True
            #    copyLayer = False
            #    doSQL = False
            #    for item in gridlayer:
            #        itemgeom = item.geometry()
            #        if itemgeom.Intersection(geometry2) is None:
            #            print('grid polygon %s intersection geometry invalid against polygon %s' %(item.GetField(0), id2))
            #            gridlayer.DeleteFeature(item.GetFID())
            #            #keep_Going = False
            #            copyLayer = True
            #            #doSQL = True
            #            continue
            #    if doSQL:
            #        #gridlayer = dst_ds.ExecuteSQL("select * from out_%s WHERE ST_IsValid(geometry)" %corenum, dialect = "SQLITE")  # Only select valid geometry
            #        gridlayer = dst_ds.ExecuteSQL("select ST_Buffer(geometry, 0), * from out_%s" %corenum, dialect = "SQLITE") # Buffer with 0 distance to create valid geometry
            #    if copyLayer:
            #        driver = ogr.GetDriverByName('ESRI Shapefile')
            #        gridds = driver.CreateDataSource(os.path.join(OutDir, '%s.shp' %id2))
            #        griddslayer = gridds.CopyLayer(gridlayer, 'gridlayer')
            #        gridds = griddslayer = driver = None
            #    if not keep_Going:
            #        continue
            #    gridlayer.ResetReading()

            # First find all intersection areas
            Areas = [[item.GetField(0), item.geometry().Intersection(geometry2).Area(), item.geometry().Area(), item.GetField(1), item.GetField(2)] for item in gridlayer]  # Only iterate over union once

            # Use the intersection area to thin the other lists
            inters = len(Areas)
            counter += inters                                                       # Advance the counter

        # Calculate area weights - for averaging
        spatialweights[id2] = [(item[0], (item[1]/polygon_area)) for item in Areas]

        # Calculate regrid weights - for conservative regridding
        regridweights[id2] = [(item[0], (item[1]/item[2])) for item in Areas]

        # Store i,j variables
        other_attributes[id2] = [[item[0], item[3], item[4]] for item in Areas]
        del gridlayer, Areas, inters, dst_ds, drv

    # Counter and printed information below
    print('[{0: >3}]      [{1: ^7} intersections processed in {2: ^4.2f} s] [{3: ^8.2f} features per second] [Processed {4: ^8} features in dest grid]'.format(corenum, counter, time.time()-tic2, (counter/(time.time()-tic2)), counter2))

    # print run information
    #print('[{0: >3}]    Done gathering intersection information between layer 1 and layer 2 in {1: 8.2f} seconds'.format(corenum, time.time()-tic2))
    #print('[{0: >3}]    {1: ^10} polygons processed for intersection with grid. {2: ^10} total polygon intersections processed.'.format(corenum, counter2, counter))
    allweights[0] = spatialweights
    allweights[1] = regridweights
    allweights[2] = other_attributes

    # pickle the dictionary into a file (saves memory)
    allweightsfile = os.path.join(OutDir, "%s.p" %corenum)
    with open(allweightsfile, "wb") as fp:
        pickle.dump(allweights, fp, protocol=pickle.HIGHEST_PROTOCOL)

    # Clean up and return
    del allweights
    return allweightsfile

##def split(a, n):
##    '''Turn a list into a number of chunks. This function is used to split up a job
##    into discrete parts and handle lits of uneven length.'''
##    k, m = len(a) / n, len(a) % n
##    if k == 0:
##        # Fewer items (a) than slots (n)
##        x = ([i,i] for i in a)
##    else:
##        x = (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))   #return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))
##    return [[y[0],y[-1],num+1] for num,y in enumerate(x)]                       # Just use the start and end of each chunk


def split(a, n):
   '''
   Turn a list into a number of chunks. This function is used to split up a job
   into discrete parts and handle lists of uneven length.
   '''
   x = numpy.array_split(a, n)
   return [[y[0], y[-1], num+1] for num,y in enumerate(x) if y.size > 0]


def main(gridder_obj, proj1, ticker, tic1, inDriverName, in_basins, fieldname, OutputGridPolys, chunk):
    '''This script code is intended to take advantage of the above functions in order
    to generate regridding weights between a netCDF file (ideally a GEOGRID file) or
    a high resolution grid in netCDF but with a GEOGRID file to define the coordinate
    reference system. The process takes the following steps:

        1) Use GEOGRID file to define the input coordinate reference system
        2) Save the NetCDF grid to in-memory raster format (optionally save to disk)
        3) Build rtree index to speed up the intersection process
        4) Calculate the spatial intersection between the input grid and the input polygons
        5) Export weights to netCDF.
        '''

    core = chunk[-1]                                                            # Identify the core number

    # Open the basin shapefile file with OGR, read-only access
    driver = ogr.GetDriverByName(inDriverName)                                  # 'GPKG' / 'ESRI Shapefile'
    shp = driver.Open(in_basins, 0)                                            # 0 means read-only. 1 means writeable.
    if shp is None:
        print("Open failed.\n")
        raise SystemExit
    if not layerName:
        layer = shp.GetLayer()
    else:
        layer = shp.GetLayer(layerName)
    #FID_column = layer.GetFIDColumn()
    FID_column = 'FID'
    SQL = '(%s >= %s) AND (%s <= %s)' %(FID_column, chunk[0], FID_column, chunk[1])
    #print('[{0: >3}]    SQL Statement: {1: ^}'.format(core, SQL))
    layer.SetAttributeFilter(SQL)

    # Perform the intersection between the layers and gather the spatial weights
    allweights = perform_intersection(gridder_obj, proj1, layer, fieldname, ticker, corenum=core)
    layer = None
    del driver, shp
    #print('[{0: >3}]    Calculation of spatial correspondences completed in {1: 8.2f} s.'.format(core, (time.time()-tic1)))
    return allweights


def work(chunk):
    '''This is the worker function which wraps all globals into one call to main()
    function. The only argument needed is the chunk, which tells OGR which features
    to grab from layer.'''
    allweights = main(gridder_obj, proj1, ticker, tic1, inDriverName, in_basins, fieldname, OutputGridPolys, chunk)
    return allweights
# ---------- End Functions ---------- #

# ---------- Begin Script Execution ---------- #
print('Script initiated at %s' %time.ctime())
tic1 = time.time()

# Ability to gather geospatial and affine transform information from a raster file
if use_inRaster:
    ds_in = gdal.Open(in_nc, gdalconst.GA_ReadOnly)
    x00, DX, xskew, y00, yskew, DY = ds_in.GetGeoTransform()
    DY = -DY
    ncols = ds_in.RasterXSize
    nrows = ds_in.RasterYSize
    proj1 = ds_in.GetSpatialRef()
    OutRaster = gdal.GetDriverByName('Mem').CreateCopy('', ds_in)
    ds_in = None
    del xskew, yskew

# Get projection information as an OSR SpatialReference object from Geogrid File
else:
    if oldCRS:
        proj1, DX, DY, x00, y00 = Read_GEOGRID_for_SRS_old(in_geogrid, geogridVariable)
    else:
        proj1, DX, DY, x00, y00 = Read_GEOGRID_for_SRS(in_geogrid)

    # Read input netcdf file into GDAL Raster object
    if not isGEOGRID and nest > 1:
        DX = DX/nest
        DY = DY/nest
    OutRaster = NetCDF_to_Raster(in_nc, Variable, proj1, DX, DY, x00, y00)          # To export original Geogrid file HGT_M variable to raster (nest = 1)

# Initiate the gridder class object
ncols = OutRaster.RasterXSize
nrows = OutRaster.RasterYSize
print('Raster will have nrows: %s and ncols: %s' %(nrows, ncols))
gridder_obj = Gridder_Layer(DX, DY, x00, y00, nrows, ncols)                           # Initiate grider class object which will perform the gridding for each basin
del ncols, nrows, DX, DY, x00, y00


if __name__ == '__main__':

    # Option to save geogrid variable to raster format
    if SaveRaster == True:
        target_ds = gdal.GetDriverByName('GTiff').CreateCopy(OutGTiff, OutRaster)
        target_ds = None
    del OutRaster

    if OutputGridPolys == True:
        create_polygons_from_info(gridder_obj, proj1, OutGridFile, outDriverName, 10000)
        print('Created output grid vector file: %s' %(OutGridFile))

    # Open the basins file in order to pull out the necessary basin geometries
    driver = ogr.GetDriverByName(inDriverName)                                  # 'GPKG' / 'ESRI Shapefile'     # Open the basin shapefile file with OGR, read-only access
    shp = driver.Open(in_basins, 0)                                             # 0 means read-only. 1 means writeable.
    if shp is None:
        print("Open failed.\n")
        raise SystemExit
    if 'layerName' in locals():
        layer = shp.GetLayer(layerName)
    else:
        layer = shp.GetLayer()
    numbasins = layer.GetFeatureCount()
    print(f"Number of basins: {numbasins}")

    # The list of IDs from the input feature class
    FID_column = layer.GetFIDColumn()
    feature_IDs = sorted([feature.GetFID() for feature in layer])
    layer.ResetReading()

    # Check for existence of provided fieldnames
    field_defn, fieldslist = checkfield(layer, fieldname, 'shapefile2')
    fieldtype2 = getfieldinfo(field_defn, fieldname)                            # Get information about field types for buildng the output NetCDF file later
    fieldtype1 = fieldtype2
    layer = shp = driver = None
    del driver, shp, layer, field_defn, fieldslist
    print('Found %s polygons in layer. Time elapsed: %s' %(numbasins, time.time()-tic1))

    print('Found %s processors. Distributing processing across %s core(s) using %s worker processes.' %(CPU_Count, Processors, pool_size))
    pool = multiprocessing.Pool(Processors)
    #chunks = list(split(range(1, numbasins+1), pool_size))                      # Split the list of basin ids into x chunks
    chunks = list(split(feature_IDs, pool_size))                                # Split the list of basin ids into x chunks
    del feature_IDs

    # The following command works
    results = pool.map(work, [chunk for chunk in chunks], chunksize=1)          # This farms the function 'work' out to the number of processors defined by 'Processors'
    pool.close()
    pool.join()
    print('Length of results: %s . Time elapsed: %s' %(len(results), time.time()-tic1))

    # Get the size of the dimensions for constructing the netCDF file
    print('Beginning to get the size of the dictionaries.')
    dim1size = 0
    dim2size = 0
    counter = 1
    for dictionary in results:
        #print('%s: Dictionary file: %s' %(counter, dictionary))
        allweights = loadpickle(dictionary)
        dim1size += len(allweights[0])
        dim2size += sum([len(item) for item in allweights[0].values()])
        allweights = None
        print('  [%s] Finished gathering dictionary length information from dictionary: %s' %(counter, dictionary))
        counter += 1
    print('Number of pickled dictionaries found: %s' %(counter-1))

    '''Create a long-vector netCDF file. '''
    # variables for compatability with the code below, which was formerly from a function
    gridflag = 1

    print('Beginning to build weights netCDF file: %s . Time elapsed: %s' %(regridweightnc, time.time()-tic1))
    tic = time.time()

    # Create netcdf file for this simulation
    rootgrp = netCDF4.Dataset(regridweightnc, 'w', format=NC_format)

    # Create dimensions and set other attribute information
    dim1name = 'polyid'
    dim2name = 'data'
    dim1 = rootgrp.createDimension(dim1name, dim1size)
    dim2 = rootgrp.createDimension(dim2name, dim2size)
    print('    Dimensions created after {0: 8.2f} seconds.'.format(time.time()-tic))

    # Handle the data type of the polygon identifier
    if fieldtype1 == 'integer':
        ids = rootgrp.createVariable(dim1name, 'i4', (dim1name))                # Coordinate Variable (32-bit signed integer)
        masks = rootgrp.createVariable('IDmask', 'i4', (dim2name))              # (32-bit signed integer)
    elif fieldtype1 == 'integer64':
        ids = rootgrp.createVariable(dim1name, 'i8', (dim1name))                # Coordinate Variable (64-bit signed integer)
        masks = rootgrp.createVariable('IDmask', 'i8', (dim2name))              # (64-bit signed integer)
    elif fieldtype1 == 'string':
        ids = rootgrp.createVariable(dim1name, str, (dim1name))                 # Coordinate Variable (string type character)
        masks = rootgrp.createVariable('IDmask', str, (dim2name))               # (string type character)
    elif fieldtype1 == 'float':
        ids = rootgrp.createVariable(dim1name, 'f4', (dim1name))                # Coordinate Variable (32-bit float)
        masks = rootgrp.createVariable('IDmask', 'f4', (dim2name))              # (32-bit float)
    print('    Coordinate variable created after {0: 8.2f} seconds.'.format(time.time()-tic))

    # Create fixed-length variables
    overlaps = rootgrp.createVariable('overlaps', 'i4', (dim1name))             # 32-bit signed integer
    weights = rootgrp.createVariable('weight', 'f8', (dim2name))                # (64-bit floating point)
    rweights = rootgrp.createVariable('regridweight', 'f8', (dim2name))         # (64-bit floating point)

    if gridflag == 1:
        iindex = rootgrp.createVariable('i_index', 'i4', (dim2name))            # (32-bit signed integer)
        jindex = rootgrp.createVariable('j_index', 'i4', (dim2name))            # (32-bit signed integer)
        iindex.long_name = 'Index in the x dimension of the raster grid (starting with 1,1 in LL corner)'
        jindex.long_name = 'Index in the y dimension of the raster grid (starting with 1,1 in LL corner)'
    print('    Variables created after {0: 8.2f} seconds.'.format(time.time()-tic))

    # Set variable descriptions
    masks.long_name = 'Polygon ID (polyid) associated with each record'
    weights.long_name = 'fraction of polygon(polyid) intersected by polygon identified by poly2'
    rweights.long_name = 'fraction of intersecting polyid(overlapper) intersected by polygon(polyid)'
    ids.long_name = 'ID of polygon'
    overlaps.long_name = 'Number of intersecting polygons'
    print('    Variable attributes set after {0: 8.2f} seconds.'.format(time.time()-tic))

    # Fill in global attributes
    rootgrp.history = 'Created %s' %time.ctime()
    rootgrp.notes = 'Input polygon file: {0}\nInput grid file: {1}\nInput nest ratio: {2}'.format(in_basins, in_nc, nest)

    # Iterate over dictionaries and begin filling in NC variable arrays
    dim1len = 0
    dim2len = 0
    counter = 1
    for dictionary in results:
        tic2 = time.time()

        # Create dictionaries
        allweights = loadpickle(dictionary)
        spatialweights = allweights[0].copy()                                       # Updates new dictionary with another one
        regridweights = allweights[1].copy()                                        # Updates new dictionary with another one
        other_attributes = allweights[2].copy()                                     # Updates new dictionary with another one
        allweights = None

        # Set dimensions for this slice
        dim1start = dim1len
        dim2start = dim2len
        dim1len += len(spatialweights)
        dim2len += sum([len(item) for item in spatialweights.values()])

        # Start filling in elements
        if fieldtype1 == 'integer':
            ids[dim1start:dim1len] = numpy.array([x[0] for x in spatialweights.items()])    # Test to fix ordering of ID values
        if fieldtype1 == 'integer64':
            #ids[dim1start:dim1len] = numpy.array([x[0] for x in spatialweights.items()], dtype=numpy.long)    # Test to fix ordering of ID values
            ids[dim1start:dim1len] = numpy.array([x[0] for x in spatialweights.items()], dtype=numpy.int64)    # Test to fix ordering of ID values
        elif fieldtype1 == 'string':
            ids[dim1start:dim1len] = numpy.array([x[0] for x in spatialweights.items()], dtype=object)    # Test to fix ordering of ID values
        elif fieldtype1 == 'float':
            ids[dim1start:dim1len] = numpy.array([x[0] for x in spatialweights.items()], dtype=numpy.float32)    # Test to fix ordering of ID values

        overlaps[dim1start:dim1len] = numpy.array([len(x) for x in spatialweights.values()])

        masklist = [[x[0] for y in x[1]] for x in spatialweights.items()]       # Get all the keys for each list of weights
        masks[dim2start:dim2len] = numpy.array([item for sublist in masklist for item in sublist], dtype=object)  # Flatten to 1 list (get rid of lists of lists)
        del masklist

        weightslist = [[item[1] for item in weight] for weight in spatialweights.values()]
        weights[dim2start:dim2len] = numpy.array([item for sublist in weightslist for item in sublist], dtype=object)
        del weightslist

        rweightlist = [[item[1] for item in rweight] for rweight in regridweights.values()]
        rweights[dim2start:dim2len] = numpy.array([item for sublist in rweightlist for item in sublist], dtype=object)
        del rweightlist

        if gridflag == 1:
            iindexlist= [[item[1] for item in attribute] for attribute in other_attributes.values()]
            iindex[dim2start:dim2len] = numpy.array([item for sublist in iindexlist for item in sublist], dtype=object)
            del iindexlist
            jindexlist = [[item[2] for item in attribute] for attribute in other_attributes.values()]
            jindex[dim2start:dim2len] = numpy.array([item for sublist in jindexlist for item in sublist], dtype=object)
            del jindexlist

        spatialweights = regridweights = other_attributes = None
        print('  [%s] Done setting dictionary: %s in %s seconds.' %(counter, os.path.basename(dictionary), time.time()-tic2))
        counter += 1
        os.remove(dictionary)                                                       # Delete the picled dictionary file

    # Close file
    rootgrp.close()
    del fieldtype1
    print('NetCDF correspondence file created in {0: 8.2f} seconds.'.format(time.time()-tic))

    # Build regrid and spatial weights NetCDF file.
    print('Finished building weights netCDF file: %s . Time elapsed: %s' %(regridweightnc, time.time()-tic1))
    # ---------- End Script Execution ---------- #
