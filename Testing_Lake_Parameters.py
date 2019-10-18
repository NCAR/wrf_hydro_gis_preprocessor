#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      ksampson
#
# Created:     10/10/2019
# Copyright:   (c) ksampson 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import time
import os
import ogr
import osr
import gdal
from gdalnumeric import *                                                       # Assists in using BandWriteArray, BandReadAsArray, and CopyDatasetInfo
from osgeo import gdal_array
import netCDF4
import csv

# Import function library into namespace. Must exist in same directory as this script.
import wrfhydro_functions as wrfh                                               # Function script packaged with this toolbox

# Import whitebox.
#from whitebox import whitebox_tools                                             # Required if first-time import
from whitebox.WBT.whitebox_tools import WhiteboxTools

in_lakes = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\NWM_v_2_1_Reservoirs_Preliminary_20190510.shp'
projdir = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\scratchdir'
fulldom = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\scratchdir\Fulldom_hires.nc'
template_raster = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\scratchdir\fill_pits.tif'
fac = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs\scratchdir\flow_acc.tif'

lakeIDfield = None
VectorDriver = 'ESRI Shapefile'                                                # Output vector file format (OGR driver name)
NoDataVal = -9999                                                               # Default NoData value for gridded variables
RasterDriver = 'GTiff'
LK_walker = 3                                                                   # Number of cells to walk downstream to get minimum lake elevation
wgs84_proj4 = '+proj=longlat +datum=WGS84 +no_defs'
minDepth = 1.0                                                                  # Minimum active lake depth for lakes with no elevation variation

# Outputs
snapPour = 'Lake_snapped_pour_points.shp'                                       # Pour points snapped downstream of lake outlets
LK_nc = 'LAKEPARM.nc'                                                           # Default Lake parameter table name [.nc]

rootgrp = netCDF4.Dataset(fulldom, 'r+')
in_raster = gdal.Open(fac, 0)
proj = osr.SpatialReference()
proj.ImportFromWkt(in_raster.GetProjectionRef())
GT = in_raster.GetGeoTransform()
nrows = in_raster.RasterYSize
ncols = in_raster.RasterXSize
in_raster = None

def project_Polygons(InputVector, outProj, clipGeom=None):
    '''
    This function is intended to project a polygon geometry to a new coordinate
    system. Optionally, the geometries can be clipped to an extent rectangle. If
    this option is chosen, the area will be re-calculated for each polygon. The
    assumption is that the linear units are meters.
    '''
    # (InputVector, outProj, clipGeom) = (in_lakes, proj, geom)
    # import ogr
    tic1 = time.time()

    # Get input vector information
    in_vect = ogr.Open(InputVector)                                             # Read the input vector file
    in_layer = in_vect.GetLayer()                                               # Get the 'layer' object from the data source
    in_proj = in_layer.GetSpatialRef()                                          # Obtain the coordinate reference object.
    inlayerDef = in_layer.GetLayerDefn()                                        # Obtain the layer definition for this layer

    # Check if a coordinate transformation (projection) must be performed
    if not outProj.IsSame(in_proj):
        print('    Input shapefile projection does not match requested output. Transforming.')
        coordTrans = osr.CoordinateTransformation(in_proj, outProj)
        trans = True

    # Create in-memory output layer to store projected and/or clipped polygons
    drv = ogr.GetDriverByName('MEMORY')                                         # Other options: 'ESRI Shapefile'
    data_source = drv.CreateDataSource('')                                      # Create the data source. If in-memory, use '' or some other string as the data source name

    # Add the fields we're interested in
    outLayer = data_source.CreateLayer('', outProj, ogr.wkbPolygon)             # Create the layer name. Use '' or some other string as the layer name

    # Build output vector file identical to input with regard to fields
    outLayer.CreateField(ogr.FieldDefn('AREASQKM', ogr.OFTReal))                # Add a single field to the new layer

    # Copy fields from input vector layer to output
    fieldNames = []                                                             # Build empty list of field names
    for i in range(inlayerDef.GetFieldCount()):
        fieldDef = inlayerDef.GetFieldDefn(i)                                   # Get the field definition for this field
        fieldName =  fieldDef.GetName()                                         # Get the field name for this field
        outLayer.CreateField(ogr.FieldDefn(fieldName, fieldDef.GetType()))      # Create a field in the output that matches the field in the input layer
        fieldNames.append(fieldName)                                            # Add field name to list of field names
        #print('    Added field {0}, dtype = {1}'.format(fieldName, fieldDef.GetType()))
    outlayerDef = outLayer.GetLayerDefn()

    # Read all features in layer
    inFeatCount = in_layer.GetFeatureCount()                                    # Get number of input features
    for feature in in_layer:
        geometry = feature.GetGeometryRef()                                     # Get the geometry object from this feature
        if trans:
            geometry.Transform(coordTrans)                                      # Transform the geometry
        if clipGeom:
            if clipGeom.Intersects(geometry):
                geometry = geometry.Intersection(clipGeom)                      # Clip the geometry if requested
            else:
                continue                                                        # Go to the next feature (do not copy)

        # Create output Feature
        outFeature = ogr.Feature(outlayerDef)                                   # Create new feature
        outFeature.SetGeometry(geometry)                                        # Set output Shapefile's feature geometry

        # Fill in fields. All fields in input will be transferred to output
        outFeature.SetField('AREASQKM', float(geometry.Area()/1000000.0))       # Add an area field to re-calculate area
        for fieldname in fieldNames:
            outFeature.SetField(fieldname, feature.GetField(fieldname))
        outLayer.CreateFeature(outFeature)                                      # Add new feature to output Layer
        feature.Destroy()                                                       # Destroy this feature
        feature = outFeature = geometry = None                                  # Clear memory
    outLayer.ResetReading()
    outFeatCount = outLayer.GetFeatureCount()                                   # Get number of features in output layer
    #outLayer = None                                                             # Clear memory

    print('    Number of output polygons: {0} of {1}'.format(outFeatCount, inFeatCount))
    print('  Completed reprojection and-or clipping in {0:3.2f} seconds.'.format(time.time()-tic1))
    in_vect = inlayerDef = in_layer = None
    return data_source, outLayer, fieldNames

def FeatToRaster(InputVector, inRaster, fieldname, dtype, NoData=None):
    '''
    This function will take a point shapefile and rasterize it. The point feature
    class must have a field in it with values of 1, which is an optional input
    to the RasterizeLayer function. Currently, this field is named PURPCODE. The
    result is a raster of NoData and 1 values, which is used as the "Input
    Depression Mask Grid" in TauDEM's PitRemove tool.
    '''
    # Python GDAL_RASTERIZE syntax, adatped from: https://gis.stackexchange.com/questions/212795/rasterizing-shapefiles-with-gdal-and-python

    # Open Raster input
    ds = gdal.Open(inRaster)

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

def raster_to_polygon(in_raster, in_proj):
    '''
    Convert a raster object to a polygon layer.
    '''
    tic1 = time.time()

    # Create temporary polygon vector layer
    ds = ogr.GetDriverByName('MEMORY').CreateDataSource('')
    Layer = ds.CreateLayer('', geom_type=ogr.wkbPolygon, srs=in_proj)   # Use projection from input raster
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


if __name__ == '__main__':

    tic1 = time.time()                                                          # Set timer

    frxst_FC = os.path.join(projdir, 'Lake_outlets.shp')

    # Setup Whitebox tools
    wbt = WhiteboxTools()
    wbt.verbose = False
    wbt.work_dir = projdir

    # Outputs
    outshp = os.path.join(projdir, 'in_lakes_clip.shp')                         # Clipped input lakes shapefile
    LakeRas = os.path.join(projdir, 'LakeGrid.tif')

    # Get domain extent for cliping geometry
    geom = grid_obj.boundarySHP('', 'MEMORY')

    # Use extent of the template raster to add a feature layer of lake polygons
    lake_ds, lake_layer, fieldNames = project_Polygons(in_lakes, proj, clipGeom=geom)
    geom = None

    # Assign a new ID field for lakes, numbered 1...n. Add field to store this information if necessary
    if lakeIDfield is None:
        print('    Adding auto-incremented lake ID field (1...n)')
        lakeID = "newID"
        if lakeID not in fieldNames:
            lake_layer.CreateField(ogr.FieldDefn(lakeID, ogr.OFTInteger))
        for num,feature in enumerate(lake_layer):
            feature.SetField(lakeID, num)      # Add an area field to re-calculate area
            lake_layer.SetFeature(feature)
        lake_layer.ResetReading()
        feature = None
    else:
        print('    Using provided lake ID field: {0}'.format(lakeIDfield))
        lakeID = lakeIDfield                                                    # Use existing field specified by 'lakeIDfield' parameter

    # Gather areas from AREASQKM field
    areas = {}
    for feature in lake_layer:
        areas[feature.GetField(lakeID)] = feature.GetField('AREASQKM')*1000000.0
        feature = None
    lake_layer.ResetReading()
    lakeIDList = list(areas.keys())

    # Save to disk in order to use in the FeatToRaster function.
    out_ds  = ogr.GetDriverByName(VectorDriver).CopyDataSource(lake_ds, outshp)
    out_ds = None

    # Convert lake geometries to raster geometries on the model grid
    LakeRaster = FeatToRaster(outshp, template_raster, lakeID, gdal.GDT_Int32, NoData=NoDataVal)
    Lake_arr = BandReadAsArray(LakeRaster.GetRasterBand(1))                     # Read raster object into numpy array
    Lake_arr[Lake_arr==0] = NoDataVal                                           # Convert 0 to WRF-Hydro NoData
    rootgrp.variables['LAKEGRID'][:] = Lake_arr                                # Write array to output netCDF file
    print('    Process: LAKEGRID written to output netCDF.')

    # Find the maximum flow accumulation value for each lake
    lake_uniques = numpy.unique(Lake_arr[Lake_arr!=NoDataVal])
    flac_arr = rootgrp.variables['FLOWACC'][:]                                  # Read flow accumulation array from Fulldom
    flac_max = {lake:flac_arr[Lake_arr==lake].max() for lake in lake_uniques}

    # Iterate over lakes, assigning the outlet pixel to the lake ID in channelgrid
    strm_arr = rootgrp.variables['CHANNELGRID'][:]                              # Read channel grid array from Fulldom
    for lake,maxfac in flac_max.items():
        strm_arr[numpy.logical_and(Lake_arr==lake, flac_arr==maxfac)] = lake
    del flac_arr
    rootgrp.variables['CHANNELGRID'][:] = strm_arr
    print('    Process: CHANNELGRID written to output netCDF.')

    # Now march down a set number of pixels to get minimum lake elevation
    strm_arr[strm_arr<1] = NoDataVal                                            # Remove channels (active and inactive)
    ds = array_to_points(strm_arr, ogr.OFTInteger, GT, proj)
    out_ds = ogr.GetDriverByName(VectorDriver).CopyDataSource(ds, frxst_FC)    # Copy to file on disk
    ds = out_ds = None
    del strm_arr, ds, out_ds

    tolerance = GT[1] * LK_walker
    snapPourFile = os.path.join(projdir, snapPour)
    wbt.snap_pour_points(frxst_FC, fac, snapPour, tolerance)                  # Snap pour points to flow accumulation grid within a tolerance

    # Rasterize the point shapefile and extract elevations at those locations
    Min_Elev_Raster = FeatToRaster(snapPourFile, fac, 'VALUE', gdal.GDT_Int32, NoData=NoDataVal)
    Min_Elev_arr = BandReadAsArray(Min_Elev_Raster.GetRasterBand(1))
    Min_Elev_Raster = None

    # Gathering minimum elevation reqiures sampling at the location below reservoir outlets
    fill_arr = rootgrp.variables['TOPOGRAPHY'][:]                               # Read elevation array from Fulldom
    min_elevs = {lake:fill_arr[Min_Elev_arr==lake].min() for lake in lake_uniques}  # Slow?
    del Min_Elev_arr
    max_elevs = {lake:fill_arr[Lake_arr==lake].max() for lake in lake_uniques}
    del Lake_arr, lake_uniques

    # Delete temporary point shapefiles
    ogr.GetDriverByName(VectorDriver).DeleteDataSource(frxst_FC)
    ogr.GetDriverByName(VectorDriver).DeleteDataSource(snapPour)

    # 2/23/2018: Find the missing lakes and sample elevation at their centroid.
    min_elev_keys = list(min_elevs.keys())
    print('    Lakes in minimum elevation dict: {0}'.format(len(min_elev_keys))) # Delete later
    MissingLks = [item for item in lakeIDList if item not in min_elev_keys]  # 2/23/2018: Find lakes that were not resolved on the grid
    shapes = {}
    if len(MissingLks) > 0:
        print('    Found {0} lakes that could not be resolved on the grid: {1}\n      Sampling elevation from the centroid of these features.'.format(len(MissingLks), str(MissingLks)))
        ds = ogr.Open(outshp, 0)
        Lakeslyr = ds.GetLayer()
        Lakeslyr.SetAttributeFilter('"%s" IN (%s)' %(lakeID, str(MissingLks)[1:-1]))    # Select the missing lakes from the input shapefile
        centroidElev = {}
        for feature in Lakeslyr:
            idval = feature.GetField(lakeID)
            centroid = feature.GetGeometryRef().Centroid()
            x, y = centroid.GetX(), centroid.GetY()
            row, col = grid_xy_to_ij(x, y, GT)
            centroidElev[idval] = fill_arr[row, col]
            shapes[idval] = (x,y)
        feature = centroid = Lakeslyr = ds = None
    del fill_arr, lakeIDList
    ogr.GetDriverByName(VectorDriver).DeleteDataSource(outshp)

    # Update dictionaries with information on the lakes that were not resolved on the grid
    max_elevs.update(centroidElev)                                              # Add single elevation value as max elevation
    min_elevs.update(centroidElev)                                              # Make these lakes the minimum depth
    del centroidElev, MissingLks

    # Give a minimum active lake depth to all lakes with no elevation variation
    elevRange = {key:max_elevs[key]-val for key,val in min_elevs.items()}   # Get lake depths
    noDepthLks = {key:val for key,val in elevRange.items() if val<minDepth}     # Make a dictionary of these lakes
    if len(noDepthLks) > 0:
        print('    Found {0} lakes with no elevation range. Providing minimum depth of %sm for these lakes.'.format(len(noDepthLks), minDepth))
        min_elevs.update({key:max_elevs[key]-minDepth for key,val in noDepthLks.items() if val==0 }) # Give these lakes a minimum depth
        noDepthFile = os.path.join(projdir, 'Lakes_with_minimum_depth.csv')
        with open(noDepthFile,'w') as f:
            w = csv.writer(f)
            for item in noDepthLks.items():
                w.writerows(noDepthLks.items())
            del noDepthFile
    del elevRange, noDepthLks

    # Calculate the Orifice and Wier heights
    OrificEs = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x])/3)) for x in min_elev_keys}             # Orific elevation is 1/3 between the low elevation and max lake elevation
    WeirE_vals = {x:(min_elevs[x] + ((max_elevs[x] - min_elevs[x]) * 0.9)) for x in min_elev_keys}       # WierH is 0.9 of the distance between the low elevation and max lake elevation

    #  Gather centroid lat/lons
    out_lake_raster_shp = os.path.join(projdir, "out_lake_raster.shp")
    ds, Layer = raster_to_polygon(LakeRaster, proj)
    out_ds  = ogr.GetDriverByName(VectorDriver).CopyDataSource(ds, out_lake_raster_shp)
    out_ds = LakeRaster = None
    #arcpy.Dissolve_management(out_lake_raster_shp, out_lake_raster_dis, "GRIDCODE", "", "MULTI_PART")               # Dissolve to eliminate multipart features
    #arcpy.Delete_management(out_lake_raster_shp)

    # Create a point geometry object from gathered lake centroid points
    print('    Starting to gather lake centroid information.')
    ds = ogr.Open(out_lake_raster_shp, 0)
    Lakeslyr = ds.GetLayer()
    wgs84_proj = osr.SpatialReference()
    wgs84_proj.ImportFromProj4(wgs84_proj4)
    coordTrans = osr.CoordinateTransformation(proj, wgs84_proj)
    cen_lats = {}
    cen_lons = {}
    for feature in Lakeslyr:
        idval = feature.GetField('RASTERVALU')
        centroid = feature.GetGeometryRef().Centroid()
        x, y = centroid.GetX(), centroid.GetY()
        centroid.Transform(coordTrans)                                      # Transform the geometry
        cen_lats[idval] = centroid.GetY()
        cen_lons[idval] = centroid.GetX()
    feature = centroid = Lakeslyr = ds = None
    del x, y, idval
    print('    Done gathering lake centroid information.')

    # Call function to build lake parameter netCDF file
    print('    Starting to create lake parameter table.')
    print('        Lakes Table: {0} Lakes'.format(len(list(areas.keys()))))
    LakeNC = os.path.join(projdir, LK_nc)
    #build_LAKEPARM(LakeNC, min_elevs, areas, max_elevs, OrificEs, cen_lats, cen_lons, WeirE_vals)

    # Delete temporary vector files
    print('    Lake parameter table created without error in {0: 3.2f} seconds.'.format(time.time()-tic1))
    #return rootgrp