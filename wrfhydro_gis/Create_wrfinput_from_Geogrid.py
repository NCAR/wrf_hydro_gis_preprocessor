# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Copyright UCAR (c) 2018
# University Corporation for Atmospheric Research(UCAR)
# National Center for Atmospheric Research(NCAR)
# Research Applications Laboratory(RAL)
# P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# 20/10/2018
#
# Name:        create_wrfinput.py
# Purpose:
# Author:      Kevin Sampson
# Created:     20/10/2018
# Licence:
#
# 10/20/2018:
#    The purpose of this script is to port the functinality of the create_wrfinput.R
#    script to Python.
#
# Based on:
#    ############################################################
#    R script to create wrfinput file from geogrid.
#    Usage: Rscript create_Wrfinput.R
#    Developed: 07/09/2017, A. Dugger
#             Mirrors the HRLDAS routines here:
#             https://github.com/NCAR/hrldas-release/blob/release/HRLDAS/HRLDAS_forcing/lib/module_geo_em.F
#             from M. Barlage.
#    ##############################################################
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

# Import Python Core Modules
import os
import sys
import getopt
import time
from argparse import ArgumentParser

# Import additional Python Modules
import numpy

# Module configurations
sys.dont_write_bytecode = True

# Screen print in case invalid parameters are given
descText = '''Utility for converting a WPS Geogrid (geo_em.d0*.nc) into a wrfinput file for initiating
                a WRF-Hydro simulation. An input Geogrid file and month with which to pull LAI is required.'''


'''
To import and run these functions using python, from a custom script or console:

    from Create_wrfinput_from_Geogrid import main_ncdfpy
    main_ncdfpy(r'./geo_em.d01.nc', r'./wrfinput_aug_pytest.nc', lai=8)
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

# The script will attempt to find the dominant SOILCTOP value. If LU_INDEX is water,
# this will be water, otherwise it will be the dominant fractional non-water SOILCTOP
# category. If no dominant land class can be determined, set the soil category to
# use as a fill below.
dom_lc_fill = 8

# User-specified soil layer thickness (dzs). Also used to calculate depth of the center of each layer (zs)
dzs = [0.1, 0.3, 0.6, 1.0]          # Soil layer thickness top layer to bottom (m)

# Missing values to use when defining netcdf variables:
missFloat = -1.e+36
missInt = -9999

#### Number of soil layers (e.g., 4)
nsoil = 4

#######################################################
# Do not update below here.
#######################################################

outNCType = 'NETCDF4'                                                           # Data model for output netCDF data

# Name the dimensions
wedim = 'west_east'                                                             # X dimension name
sndim = 'south_north'                                                           # Y-dimension name
timedim = 'Time'                                                                # Time dimension name
soildim = 'soil_layers_stag'                                                    # Soil layer dimension name
keepDims = [timedim, 'month', sndim, wedim, soildim]                            # Dimensions to transfer from GEOGRID to WRFINPUT

# Alter the names of variable names from GEOGRID (key) to WRFINPUT (value) variable names
mapVars = {'HGT_M': 'HGT',
            'XLAT_M': 'XLAT',
            'XLONG_M': 'XLONG',
            'LU_INDEX': 'IVGTYP'}

# Variables to keep from GEOGRID
keepVars = ['XLAT_M','XLONG_M','HGT_M','LU_INDEX','MAPFAC_MX', 'MAPFAC_MY']

# (Name, units, dimensions, missing value, dtype) for all variables to add to the WRFINPUT file
addVars = [('TMN', 'K', [timedim, sndim, wedim], missFloat, 'f4'),
            ('XLAND', '', [timedim, sndim, wedim], missInt, 'i4'),
            ('SEAICE', '', [timedim, sndim, wedim], missFloat, 'f4'),
            ('ISLTYP', '', [timedim, sndim, wedim], missInt, 'i4'),
            ('SHDMAX', '%', [timedim, sndim, wedim], missFloat, 'f4'),
            ('SHDMIN', '%', [timedim, sndim, wedim], missFloat, 'f4'),
            ('LAI', 'm^2/m^2', [timedim, sndim, wedim], missFloat, 'f4'),
            ('CANWAT', 'kg/m^2', [timedim, sndim, wedim], missFloat, 'f4'),
            ('SNOW', 'kg/m^2', [timedim, sndim, wedim], missFloat, 'f4'),
            ('TSK', 'K', [timedim, sndim, wedim], missFloat, 'f4'),
            ('SMOIS', 'm^3/m^3', [timedim, soildim, sndim, wedim], missFloat, 'f4'),
            ('TSLB', 'K', [timedim, soildim, sndim, wedim], missFloat, 'f4'),
            ('ZS', 'm', [timedim, soildim], missFloat, 'f4'),
            ('DZS', 'm', [timedim, soildim], missFloat, 'f4')]

# Choose the processing method for producing the output. Options = ['netcdf4-python', 'xarray']
method = 'xarray'
if method == 'netcdf4-python':
    import netCDF4
elif method == 'xarray':
    import xarray as xr

# Switch for fixing SoilT values of 0 over water areas
fix_zero_over_water = True

# --- Functions --- #
def is_valid_file(parser, arg):
    # https://stackoverflow.com/questions/11540854/file-as-command-line-argument-for-argparse-error-message-if-argument-is-not-va
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return str(arg)

def fill_wrfinput_ncdfpy(rootgrp_in, rootgrp_out, laimo=8):
    '''
    This function will populate the arrays in the WRFINPUT file based on array and
    attribute values in the input GEOGRID file.
    '''

    # Pull data variables from input GEOGRID file
    soilT = rootgrp_in.variables['SOILTEMP'][:]                                 # Soil temperature in degrees Kelvin
    hgt = rootgrp_in.variables['HGT_M'][:]                                      # Elevation in meters on the mass grid
    use = rootgrp_in.variables['LU_INDEX'][:]                                   # Landcover class
    soil_top_cat = rootgrp_in.variables['SOILCTOP'][:]                          # Soil texture class of the top soil layer
    veg = rootgrp_in.variables['GREENFRAC'][:] * 100.0                          # Green fraction as a percentage (0-100)

    # Define variables using global attributes in GEOGRID file
    ncatts = {key:val for key,val in rootgrp_out.__dict__.items()}
    iswater = ncatts.get('ISWATER')                                             # GEOGRID global attribute describing the LU_INDEX value assigned to water
    isoilwater = ncatts.get('ISOILWATER')                                       # GEOGRID global describing the SOIL class assigned to water

    # SOILTEMP will show 0 value over water. This can cause issues when varying land cover fields
    # from default. Setting to mean non-zero values for now to have something reasonable.
    if fix_zero_over_water:
        soilT_mask = soilT<100                                                      # Create a mask of 'invalid' soil temperature values. All values <100K on SOILTEMP grid
        soilT_mask_Mean = soilT[~soilT_mask].mean()                                 # Calculate the non-NAN soil temperature mean
        soilT[soilT_mask] = soilT_mask_Mean                                         # Replace masked values with the mean of the unmasked values
        if soilT_mask.sum()>0:
            print('    Replaced {0} values in TMN with mean SOILTEMPT value ({1}).'.format(soilT_mask.sum(), soilT_mask_Mean))

    # TMN topographic adjustment
    print('    Performing topographic soil temperature adjustment.')
    tmn = soilT - 0.0065 * hgt                                                  # Topographic soil temperature adjustment
    del soilT, hgt, soilT_mask, soilT_mask_Mean

    # Create land/water mask from LU_INDEX in input GEOGRID file
    msk = use.copy()                                                            # Land use (LU_INDEX from GEOGRID)
    msk[numpy.logical_and(msk>=0, msk!=iswater)] = 1                            # Set non-water cells to 1
    msk[msk==iswater] = 2                                                       # Set all other cells to 2

    # Find the dominant SOILCTOP value. If LU_INDEX is water, this will be water,
    # otherwise it will be the dominant fractional non-water SOILCTOP category
    a = numpy.ma.array(soil_top_cat[0], mask=False)                             # Build a masked array in order to mask the water cateogry
    a.mask[isoilwater-1] = True                                                 # Use a mask to mask the index of the water category
    dominant_value = numpy.amax(a, axis=0).data                                 # Dominant soil category fraction (unmasked)
    dominant_index = numpy.argmax(a, axis=0)+1                                  # Dominant soil category index (1-based index)
    dominant_index[numpy.logical_and(dominant_value < 0.01, msk[0]==1)] = dom_lc_fill     # Set any land pixel values with tiny dominant fractions to a fill value
    if numpy.logical_and(dominant_value < 0.01, msk[0]==1).sum() > 0:
        print('    Replaced {0} values in ISLTYP with {1} because no dominant land class could be determined'.format(numpy.logical_and(dominant_value < 0.01, msk[0]==1).sum(), dom_lc_fill))
    dominant_index[msk[0]==2] = isoilwater                                      # Set any values with landmask == 2 (water) to water
    del a, dominant_value, soil_top_cat

    # Set soils to the dominant type (derived above) and handle water/land soiltype mismatches
    soi = dominant_index[numpy.newaxis]                                         # Add an axis to make this 3D (time, south_north, west_east)
    soi[use==iswater] = isoilwater                                              # Make all water LU_INDEX cells into water soiltypes
    soi[numpy.logical_and(use!=iswater, soi==isoilwater)] = fillsoiltyp         # If the pixel LU_INDEX is land type and SOILCTOP is water type, then make it water soil type
    if numpy.logical_and(use!=iswater, soi==isoilwater).sum() > 0:
        print('    Replaced {0} values in ISLTYP with {1} because of a land landcover type and water soil class'.format(fillsoiltyp, numpy.logical_and(use!=iswater, soi==isoilwater).sum()))
    del dominant_index, use

    # Soil moisture SMOIS 3D array
    smoisArr = numpy.array([0.20, 0.21, 0.25, 0.27])                            # Constant soil moisture with increasing depth by vertical level
    smois = smoisArr[:, None, None] * numpy.ones(msk.shape)                     # Set the soil moisture (SMOIS) array across entire domain by vertical level

    # TSLB 3D array
    tslbArr = numpy.array([285.0, 283.0, 279.0, 277.0])                         # Constant tslb with increasing depth by vertical level
    tslb = tslbArr[:, None, None] * numpy.ones(msk.shape)                       # Set the TSLB array across entire domain by vertical level. tslb = numpy.vstack([msk]*4)

    # Calculate the depths of the center depth of each soil layer based on the
    # layer thicknesses provided in the header. Default is zs = [0.05, 0.25, 0.7, 1.5]
    zs = [item/2 + sum(dzs[:num]) for num,item in enumerate(dzs)]               # Each center depth is half the layer thickness + sum of thicknesses of all levels above

    # Populate output WRFINPUT file variable arrays
    rootgrp_out.variables['TMN'][:] = tmn                                       # Elevation adjusted deep soil temperature
    rootgrp_out.variables['XLAND'][:] = msk                                     # Landmask (1=land, 2=water) from LU_INDEX
    rootgrp_out.variables['SEAICE'][:] = numpy.zeros(msk.shape)                 # zeros
    rootgrp_out.variables['ISLTYP'][:] = soi                                    # Dominant soil type
    rootgrp_out.variables['SHDMAX'][:] = veg.max(axis=1)                        # Maximum GREENFRAC over time dimesion
    rootgrp_out.variables['SHDMIN'][:] = veg.min(axis=1)                        # Minimum GREENFRAC over time dimesion
    rootgrp_out.variables['LAI'][:] = rootgrp_in.variables['LAI12M'][:,laimo-1] # Leaf area index for the user-specified month
    rootgrp_out.variables['CANWAT'][:] = numpy.zeros(msk.shape)                 # zeros
    rootgrp_out.variables['SNOW'][:] = numpy.zeros(msk.shape)                   # zeros
    rootgrp_out.variables['TSK'][:] = numpy.zeros(msk.shape) + 290.0            # Constant value
    rootgrp_out.variables['SMOIS'][:] = smois[numpy.newaxis]                    # Add an axis to make this 4D (time, soil_layer_stag, south_north, west_east)
    rootgrp_out.variables['TSLB'][:] = tslb[numpy.newaxis]                      # Add an axis to make this 4D (time, soil_layer_stag, south_north, west_east)
    rootgrp_out.variables['ZS'][:] = numpy.array(zs)[numpy.newaxis]             # Depths of the center of each soil layer
    rootgrp_out.variables['DZS'][:] = numpy.array(dzs)[numpy.newaxis]           # Thickness of each soil layer
    del msk, veg, ncatts, iswater, isoilwater, soi, smois, smoisArr, tslb, tslbArr, tmn, zs
    return rootgrp_in, rootgrp_out

def main_wrfinput_ncdfpy(geoFile, wrfinFile, lai=8, outNCType='NETCDF4'):
    '''
    This function is designed to build the wrfinput file using only the
    netCDF4-python library.
    '''
    tic1 = time.time()

    # --- Create initial file --- #
    print('  Creating wrfinput file from geogrid file.')
    print('    Input geogrid file: {0}'.format(geoFile))
    print('    Output wrfinput file: {0}'.format(wrfinFile))
    print('    Month selected (1=Januaray, 12=December): {0}'.format(lai))
    rootgrp_in = netCDF4.Dataset(geoFile, 'r')                                  # Read object on input GEOGRID file
    rootgrp_out = netCDF4.Dataset(wrfinFile, 'w', format=outNCType)             # Write object to create output WRFINPUT file

    # Copy dimensions from GEOGRID file
    for dimname, dim in rootgrp_in.dimensions.items():
        if dimname in keepDims:
            rootgrp_out.createDimension(dimname, len(dim))                      # Copy dimensions from the GEOGRID file
    soildimension = rootgrp_out.createDimension(soildim, nsoil)           # Add soil_layers_stag dimension

    # Populate initial file with variables to keep from the input GEOGRID file
    for varname, ncvar in rootgrp_in.variables.items():
        if varname in keepVars:
            varAtts = {key:val for key,val in ncvar.__dict__.items()}
            varDims = tuple(varDim for varDim in ncvar.dimensions)
            if varname in mapVars:
                varname = mapVars[varname]                                      # Alter the variable name
            var = rootgrp_out.createVariable(varname, ncvar.dtype, varDims)
            var.setncatts(varAtts)                                              # Copy variable attributes from GEOGRID file

    # Define new variables based on the addVars list
    for (varname, units, varDims, missing_value, dtype) in addVars:
        var = rootgrp_out.createVariable(varname, dtype, varDims)
        var.setncatts({'units':units, 'missing_value':missing_value})

    # Global Attributes - copy all from GEOGRID file to WRFINPUT
    ncatts = {key:val for key,val in rootgrp_in.__dict__.items()}
    ncatts['Source_Software'] = 'WRF-Hydro {0} script (Python).'.format(sys.argv[0])
    ncatts['creation_time'] = 'Created {0}'.format(time.ctime())
    rootgrp_out.setncatts(ncatts)

    # Add variable values last (makes the script run faster). These are exact copies of existing variables
    for varname, ncvar in rootgrp_in.variables.items():
        if varname in keepVars:
            if varname in mapVars:
                varname = mapVars[varname]
            var = rootgrp_out.variables[varname]
            var[:] = ncvar[:]                                                   # Copy the variable data into the newly created variable

    # Process and populate variables
    rootgrp_in, rootgrp_out = fill_wrfinput_ncdfpy(rootgrp_in, rootgrp_out, laimo=lai)

    # Close and return
    rootgrp_in.close()
    rootgrp_out.close()
    return

def fill_wrfinput_xarray(ds_in, laimo=8):
    '''
    This function will populate the arrays in the WRFINPUT file based on array and
    attribute values in the input GEOGRID file.
    '''

    # Define variables using global attributes in GEOGRID file
    iswater = ds_in.attrs.get('ISWATER')                                        # GEOGRID global attribute describing the LU_INDEX value assigned to water
    isoilwater = ds_in.attrs.get('ISOILWATER')                                  # GEOGRID global describing the SOIL class assigned to water

    # SOILTEMP will show 0 value over water. This can cause issues when varying land cover fields
    # from default. Setting to mean non-zero values for now to have something reasonable.
    hgt = ds_in['HGT'].data                                                     # Elevation in meters on the mass grid
    if fix_zero_over_water:
        soilT = ds_in['SOILTEMP'].data.copy()                                       # Soil temperature in degrees Kelvin
        soilT_mask = soilT<100                                                      # Create a mask of 'invalid' soil temperature values. All values <100K on SOILTEMP grid
        soilT_mask_Mean = soilT[~soilT_mask].mean()                                 # Calculate the non-NAN soil temperature mean
        soilT[soilT_mask] = soilT_mask_Mean                                         # Replace masked values with the mean of the unmasked values
        if soilT_mask.sum()>0:
            print('    Replaced {0} values in TMN with mean SOILTEMPT value ({1}).'.format(soilT_mask.sum(), soilT_mask_Mean))

    # TMN topographic adjustment
    print('    Performing topographic soil temperature adjustment.')
    tmn = soilT - 0.0065 * hgt                                                  # Topographic soil temperature adjustment
    del soilT, hgt, soilT_mask, soilT_mask_Mean

    # Create land/water mask from LU_INDEX in input GEOGRID file
    use = ds_in['IVGTYP'].data                                                  # Landcover class
    msk = use.copy()                                                            # Land use (LU_INDEX from GEOGRID)
    msk[numpy.logical_and(msk>=0, msk!=iswater)] = 1                            # Set non-water cells to 1
    msk[msk==iswater] = 2                                                       # Set all other cells to 2

    # Find the dominant SOILCTOP value. If LU_INDEX is water, this will be water,
    # otherwise it will be the dominant fractional non-water SOILCTOP category
    soil_top_cat = ds_in['SOILCTOP'].data                                       # Soil texture class of the top soil layer
    a = numpy.ma.array(soil_top_cat[0], mask=False)                             # Build a masked array in order to mask the water cateogry
    a.mask[isoilwater-1] = True                                                 # Use a mask to mask the index of the water category
    dominant_value = numpy.amax(a, axis=0).data                                 # Dominant soil category fraction (unmasked)
    dominant_index = numpy.argmax(a, axis=0)+1                                  # Dominant soil category index (1-based index)
    dominant_index[numpy.logical_and(dominant_value < 0.01, msk[0]==1)] = dom_lc_fill     # Set any land pixel values with tiny dominant fractions to a fill value
    if numpy.logical_and(dominant_value < 0.01, msk[0]==1).sum() > 0:
        print('    Replaced {0} values in ISLTYP with {1} because no dominant land class could be determined'.format(numpy.logical_and(dominant_value < 0.01, msk[0]==1).sum(), dom_lc_fill))
    dominant_index[msk[0]==2] = isoilwater                                      # Set any values with landmask == 2 (water) to water
    del a, dominant_value, soil_top_cat

    # Set soils to the dominant type (derived above) and handle water/land soiltype mismatches
    soi = dominant_index[numpy.newaxis]                                         # Add an axis to make this 3D (time, south_north, west_east)
    soi[use==iswater] = isoilwater                                              # Make all water LU_INDEX cells into water soiltypes
    soi[numpy.logical_and(use!=iswater, soi==isoilwater)] = fillsoiltyp         # If the pixel LU_INDEX is land type and SOILCTOP is water type, then make it water soil type
    if numpy.logical_and(use!=iswater, soi==isoilwater).sum() > 0:
        print('    Replaced {0} values in ISLTYP with {1} because of a land landcover type and water soil class'.format(fillsoiltyp, numpy.logical_and(use!=iswater, soi==isoilwater).sum()))
    del dominant_index, use

    # Soil moisture SMOIS 3D array
    smoisArr = numpy.array([0.20, 0.21, 0.25, 0.27])                            # Constant soil moisture with increasing depth by vertical level
    smois = smoisArr[:, None, None] * numpy.ones(msk.shape)                     # Set the soil moisture (SMOIS) array across entire domain by vertical level

    # TSLB 3D array
    tslbArr = numpy.array([285.0, 283.0, 279.0, 277.0])                         # Constant tslb with increasing depth by vertical level
    tslb = tslbArr[:, None, None] * numpy.ones(msk.shape)                       # Set the TSLB array across entire domain by vertical level. tslb = numpy.vstack([msk]*4)

    # Calculate the depths of the center depth of each soil layer based on the
    # layer thicknesses provided in the header. Default is zs = [0.05, 0.25, 0.7, 1.5]
    zs = [item/2 + sum(dzs[:num]) for num,item in enumerate(dzs)]               # Each center depth is half the layer thickness + sum of thicknesses of all levels above

    veg = ds_in['GREENFRAC'].data * 100.0                               # Green fraction as a percentage (0-100)

    # Populate output WRFINPUT file variable arrays
    ds_in.variables['TMN'][:] = tmn                                       # Elevation adjusted deep soil temperature
    ds_in.variables['XLAND'][:] = msk                                     # Landmask (1=land, 2=water) from LU_INDEX
    ds_in.variables['SEAICE'][:] = numpy.zeros(msk.shape)                 # zeros
    ds_in.variables['ISLTYP'][:] = soi                                    # Dominant soil type
    ds_in.variables['SHDMAX'][:] = veg.max(axis=1)                        # Maximum GREENFRAC over time dimesion
    ds_in.variables['SHDMIN'][:] = veg.min(axis=1)                        # Minimum GREENFRAC over time dimesion
    ds_in.variables['LAI'][:] = ds_in.variables['LAI12M'][:,laimo-1] # Leaf area index for the user-specified month
    ds_in.variables['CANWAT'][:] = numpy.zeros(msk.shape)                 # zeros
    ds_in.variables['SNOW'][:] = numpy.zeros(msk.shape)                   # zeros
    ds_in.variables['TSK'][:] = numpy.zeros(msk.shape) + 290.0            # Constant value
    ds_in.variables['SMOIS'][:] = smois[numpy.newaxis]                    # Add an axis to make this 4D (time, soil_layer_stag, south_north, west_east)
    ds_in.variables['TSLB'][:] = tslb[numpy.newaxis]                      # Add an axis to make this 4D (time, soil_layer_stag, south_north, west_east)
    ds_in.variables['ZS'][:] = numpy.array(zs)[numpy.newaxis]             # Depths of the center of each soil layer
    ds_in.variables['DZS'][:] = numpy.array(dzs)[numpy.newaxis]           # Thickness of each soil layer
    del msk, veg, iswater, isoilwater, soi, smois, smoisArr, tslb, tslbArr, tmn, zs
    return ds_in

def main_wrfinput_xarray(geoFile, wrfinFile, lai=8, outNCType='NETCDF4'):
    '''
    This function is designed to build the wrfinput file using the xarray library.
    '''
    tic1 = time.time()

    # --- Create initial file --- #
    print('  Creating wrfinput file from geogrid file.')
    print('    Input geogrid file: {0}'.format(geoFile))
    print('    Output wrfinput file: {0}'.format(wrfinFile))
    print('    Month selected (1=Januaray, 12=December): {0}'.format(lai))

    # Open input files for reading. Lazy loading
    ncDS = xr.open_dataset(geoFile)                         # , decode_cf=False

    # Rename variables
    ncDS = ncDS.rename(mapVars)

    # Add new variables based on the addVars list
    dims = dict(ncDS.dims)
    dims.update({soildim:nsoil})
    newVars = []
    for (varname, units, varDims, missing_value, dtype) in addVars:
        da = xr.DataArray(data=numpy.empty(tuple([dims[dim] for dim in varDims]), dtype=dtype),
            dims=varDims,
            attrs={'units':units, 'missing_value':missing_value})
        ncDS[varname] = da
        newVars.append(varname)

    # Process and populate variables
    ncDS = fill_wrfinput_xarray(ncDS, laimo=lai)

    # Drop dimensions
    dropDims = [item for item in ncDS.dims if item not in keepDims]
    ncDS = ncDS.drop_dims(dropDims)

    # Drop variables first
    keepVars2 = [mapVars.get(item,item) for item in keepVars] + newVars
    dropVars = [item for item in ncDS.variables if item not in keepVars2]
    ncDS = ncDS.drop(dropVars)

    # Add global attributes
    ncDS.attrs['Source_Software'] = 'WRF-Hydro {0} script (Python).'.format(sys.argv[0])
    ncDS.attrs['creation_time'] = 'Created {0}'.format(time.ctime())

    # Output file to disk
    #encoding = {varname:ncDS[varname].encoding for varname in list(ncDS.variables.keys())}
    encoding = {varname:{'_FillValue':None} for varname in list(ncDS.variables.keys())}
    #for key, val in encoding.items():
    #    val['_FillValue'] = None
    ncDS.to_netcdf(wrfinFile, mode='w', format=outNCType, encoding=encoding)
    ncDS.close()
    del encoding, ncDS
    return

# --- End Functions --- #

if __name__ == '__main__':
    print('Script initiated at {0}'.format(time.ctime()))
    tic = time.time()

    parser = ArgumentParser(description=descText, add_help=True)
    parser.add_argument("-i",
                        dest="in_Geogrid",
                        type=lambda x: is_valid_file(parser, x),
                        required=True,
                        help="Path to WPS geogrid (geo_em.d0*.nc) file [REQUIRED]")
    parser.add_argument("-m",
                        dest="LAI_month",
                        type=int,
                        default=8,
                        required=True,
                        help="LAI month for initialization [REQUIRED]. 1=January, 12=December.")
    parser.add_argument("-o",
                        dest="out_wrfinput",
                        default='./wrfinput.nc',
                        required=True,
                        help='Output "wrfinput" file.')

    # If no arguments are supplied, print help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    if args.LAI_month not in range(1,13):
        print('LAI month provided [{0}] is not between 1 and 12. Exiting...'.format(args.LAI_month))
        raise SystemExit

    # Resolve relative paths to absolute paths
    args.in_Geogrid = os.path.abspath(args.in_Geogrid)
    args.out_wrfinput = os.path.abspath(args.out_wrfinput)

    if method == 'netcdf4-python':
        main_wrfinput_ncdfpy(args.in_Geogrid, args.out_wrfinput, lai=args.LAI_month, outNCType=outNCType)
    elif method == 'xarray':
        main_wrfinput_xarray(args.in_Geogrid, args.out_wrfinput, lai=args.LAI_month, outNCType=outNCType)
    print('  Process completed in {0:3.2f} seconds'.format(time.time()-tic))