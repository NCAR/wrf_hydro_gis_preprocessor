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
import numpy

# Import additional modules
import ogr
import osr

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import CSV_to_SHP                                       # Function script packaged with this toolbox

print('Script initiated at {0}'.format(time.ctime()))

# Global Variables

# Inputs and output directory
inCSV = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\Iowa_Gauges_v8.csv'
inSHP = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Inputs\Temp_Frxst_Pts.shp'
out_dir = r'C:\Users\ksampson\Desktop\WRF_Hydro_GIS_Preprocessor_FOSS\Outputs'

# Otput files
outCSV = os.path.join(out_dir, os.path.basename(inCSV).replace('.shp', '.csv'))
outSHP = os.path.join(out_dir, os.path.basename(inSHP).replace('.csv', '.shp'))

# Coordinate system of all latitude/longitude coordinates: WGS84, EPSG:4326
wgs84_proj4 = '+proj=longlat +datum=WGS84 +no_defs'

Driver = 'ESRI Shapefile'

# Script options
CSV_to_shape = False                                   # Switch for creating a point shapefile from a forecast point CSV file
SHP_to_CSV = True                                  # Switch for creating a CSV file from a point shapefile

# Dictionary to map OGR data types to numpy dtypes - many of these are just a guess,
# with numpy.object for List, string, and binary objects
OGRTypes = {ogr.OFTBinary: numpy.object,
            ogr.OFTDate: numpy.datetime64,
            ogr.OFTDateTime: numpy.datetime64,
            ogr.OFTInteger: numpy.int,
            ogr.OFTInteger64: numpy.int64,
            ogr.OFTInteger64List: numpy.object,
            ogr.OFTIntegerList: numpy.object,
            ogr.OFTReal: numpy.float64,
            ogr.OFTRealList: numpy.object,
            ogr.OFTString: numpy.str,
            ogr.OFTStringList: numpy.object,
            ogr.OFTTime: numpy.datetime64,
            ogr.OFTWideString: numpy.object,
            ogr.OFTWideStringList: numpy.object}

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()

    drv = ogr.GetDriverByName(Driver)
    if drv is None:
        print('      {0} driver not available.'.format(DriverName))
        raise SystemExit

    if CSV_to_shape:
        '''
        This block will create a point shapefile from an input CSV file.
        '''
        ds = CSV_to_SHP(inCSV, DriverName='MEMORY')
        out_ds = drv.CopyDataSource(ds, outSHP)
        out_ds = ds = None

    if SHP_to_CSV:

        '''
        This block will first use any existing fields to fill out the output
        CSV file from the list: ['FID', 'LAT', 'LON']. If any of these are not
        present, it will fill them in with values calculated from the geometry
        (['LAT', 'LON']) or from a 1...n numbering (['FID']).
        '''

        data_source = drv.Open(inSHP, 0)                                 # 0 means read-only. 1 means writeable.
        if data_source is None:
            print('      data source could not be created.')
            raise SystemExit

        layer = data_source.GetLayer()
        featureCount = layer.GetFeatureCount()
        inSR = layer.GetSpatialRef()                                            # Get spatial reference of input
        print('Number of features in {0}: {1}'.format(os.path.basename(inSHP),featureCount))

        # Read shapefile fields
        layerDefinition = layer.GetLayerDefn()
        npDtypes = []
        fields = []
        for i in range(layerDefinition.GetFieldCount()):
            fieldDef = layerDefinition.GetFieldDefn(i)
            fieldName = fieldDef.GetName()
            fieldType = fieldDef.GetType()
            npDtypes.append((fieldName, OGRTypes[fieldType]))
            fields.append(fieldName)

        # Read the input CSV file
        append_dtypes = [('LAT', numpy.float64), ('LON', numpy.float64), ('FID', numpy.int)]
        addDtypes = [item for item in append_dtypes if item[0] not in [item2[0] for item2 in npDtypes]]    # Eliminate redundancy
        dtypes = numpy.dtype(npDtypes + addDtypes)                              # Create numpy dtype object, adding in any required fields
        csv_arr = numpy.empty(featureCount, dtype=dtypes)

        # create the spatial reference for output, WGS84
        outSR = osr.SpatialReference()
        outSR.ImportFromProj4(wgs84_proj4)

        # Create a transformation?
        if outSR.IsSame(inSR):
            print('    No coordinate transformation necessary.')
            mustTransform = False
        else:
            print('    Transforming coordinates from:\n\t\t{0}\n\t  to\n\t\t{1}.'.format(inSR.ExportToProj4(), outSR.ExportToProj4()))
            mustTransform = True
            transform = osr.CoordinateTransformation(inSR, outSR)

        # Build the numpy array using info in the fields or geometries of the shapefile
        for num,feature in enumerate(layer):
            geom = feature.GetGeometryRef().Clone()
            if mustTransform:
                geom.Transform(transform)
            csv_arr[fields][num] = tuple([feature.GetField(field) for field in fields])
            if 'LAT' not in npDtypes:
                csv_arr['LAT'][num] = geom.GetY(0)
            if 'LON' not in npDtypes:
                csv_arr['LON'][num] = geom.GetX(0)
            if 'FID' not in npDtypes:
                csv_arr['FID'][num] = num
        layer.ResetReading()

        numpy.savetxt(outCSV, csv_arr, fmt='%s', delimiter=',', header=','.join(csv_arr.dtype.names), comments='')
    drv = None
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))
