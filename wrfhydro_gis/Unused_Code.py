



# Coastline harmonization with a landmask
def coastlineHarmonize(maskFile, ds, outmaskFile, outDEM, minimum, waterVal=0):
    '''
    This function is designed to take a coastline mask and harmonize elevation
    values to it, such that no elevation values that are masked as water cells
    will have elevation >0, and no land cells will have an elevation < minimum.
    '''
    tic1 = time.time()

    # Read mask file for information
    refDS = gdal.Open(maskFile, gdalconst.GA_ReadOnly)
    target_ds = gdal.GetDriverByName(RasterDriver).Create(outmaskFile, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte)
    DEM_ds = gdal.GetDriverByName(RasterDriver).Create(outDEM, ds.RasterXSize, ds.RasterYSize, 1, ds.GetRasterBand(1).DataType)
    CopyDatasetInfo(ds, target_ds)                                              # Copy information from input to output
    CopyDatasetInfo(ds, DEM_ds)                                                 # Copy information from input to output

    # Resample input to output
    gdal.ReprojectImage(refDS, target_ds, refDS.GetProjection(), target_ds.GetProjection(), gdalconst.GRA_NearestNeighbour)

    # Build numpy array of the mask grid and elevation grid
    maskArr = BandReadAsArray(target_ds.GetRasterBand(1))
    elevArr = BandReadAsArray(ds.GetRasterBand(1))

    # Reassign values
    ndv = ds.GetRasterBand(1).GetNoDataValue()                                  # Obtain nodata value
    mask = maskArr==1                                                           # A boolean mask of True wherever LANDMASK=1
    elevArr[elevArr==ndv] = 0                                                   # Set Nodata cells to 0
    elevArr[mask] += minimum                                                    # For all land cells, add minimum elevation
    elevArr[~mask] = waterVal                                                   # ds.GetRasterBand(1).GetNoDataValue()

    # Write to output
    band = DEM_ds.GetRasterBand(1)
    BandWriteArray(band, elevArr)
    band.SetNoDataValue(ndv)

    # Clean up
    target_ds = refDS = DEM_ds = band = None
    del maskArr, elevArr, ndv, mask
    print('    DEM harmonized with landmask in %3.2f seconds.' %(time.time()-tic1))

def raster_extent(in_raster):
    '''
    Given a raster object, return the bounding extent [xMin, yMin, xMax, yMax]
    '''
    xMin, DX, xskew, yMax, yskew, DY = in_raster.GetGeoTransform()
    Xsize = in_raster.RasterXSize
    Ysize = in_raster.RasterYSize
    xMax = xMin + (float(Xsize)*DX)
    yMin = yMax + (float(Ysize)*DY)
    del Xsize, Ysize, xskew, yskew, DX, DY
    return [xMin, yMin, xMax, yMax]

def alter_GT(GT, regridFactor):
    '''
    This function will alter the resolution of a raster's affine transformation,
    assuming that the extent and CRS remain unchanged.
    '''
    # Georeference geogrid file
    GeoTransform = list(GT)
    DX = GT[1]/float(regridFactor)
    DY = GT[5]/float(regridFactor)
    GeoTransform[1] = DX
    GeoTransform[5] = DY
    GeoTransformStr = ' '.join([str(item) for item in GeoTransform])
    return GeoTransform, GeoTransformStr, DX, DY

# Function to reclassify values in a raster
def reclassifyRaster(array, thresholdDict):
    '''
    Apply a dictionary of thresholds to an array for reclassification.
    This function may be made more complicated as necessary
    '''
    # Reclassify array using bounds and new classified values
    new_arr = array.copy()
    for newval, oldval in thresholdDict.iteritems():
        mask = numpy.where(array==oldval)
        new_arr[mask] = newval
    del array
    return new_arr

# Function to calculate statistics on a raster using gdalinfo command-line
def calcStats(inRaster):
    print('    Calculating statistics on %s' %inRaster)
    subprocess.call('gdalinfo -stats %s' %inRaster, shell=True)

def apply_threshold(array, thresholdDict):
    '''
    Apply a dictionary of thresholds to an array for reclassification.
    This function may be made more complicated as necessary
    '''

    # Reclassify array using bounds and new classified values
    for newval, bounds in thresholdDict.iteritems():
        mask = numpy.where((array > bounds[0]) & (array <= bounds[1]))          # All values between bounds[0] and bounds[1]
        array[mask] = newval
    return array