# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# Copyright UCAR (c) 2018
# University Corporation for Atmospheric Research(UCAR)
# National Center for Atmospheric Research(NCAR)
# Research Applications Laboratory(RAL)
# P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# 14/12/2018
#
# Name:        create_SoilProperties.py
# Purpose:
# Author:      Kevin Sampson, NCAR
# Created:     14/12/2018
# Licence:
#
# Based on:
#    #############################################################
#    R script to create spatial parameter files from TBLs.
#    Usage: Rscript create_SoilProperties.R
#    Developed: 11/11/2016, A. Dugger
#    Updated: 07/23/2017, A.Dugger
#             New capability to handle soil composition
#             fractions and convert to soil parameters. Also
#             new functionality to create HYDRO2DTBL.nc.
#             Pedotransfer and soil texture class functions
#             from M. Barlage (7/2017).
#    #############################################################
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

'''
09/09/2024
'''

# Screen print in case invalid parameters are given
descText = '''Usage: python %s <GEOGRID file> <OUTPUT Directory>
              <GEOGRID file>                -> Input GEOGRID netCDF file
              <Parameter File Directory>    -> Directory containing SOILPARM.TBL, MPTABLE.TBL, GENPARM.TBL, and HYDRO.TBL
              <Output Directory>            -> Directory to write output files to (soil_properties.nc, hydro2dtbl.nc)
              '''

# Import python core modules
import os
import sys
import time
import re
import shutil
from io import StringIO

# Import additional modules
import f90nml
import netCDF4
import numpy as np
import pandas as pd
from argparse import ArgumentParser

# --- Global Variables --- #

#######################################################
# Update relevant arguments below.
#######################################################

#### Land cover classification system?
# Options: ["USGS", "MODIS"]
landClass = "USGS"

# Select which soil parameter set to use from SOILPARM.TBL ['STAS', 'STAS-RUC']
soilparm_set = 'STAS'

# Use SOILCOMP?
useSoilComp = False

# Update texture class in geogrid??
# Note that if TRUE, the script will overwrite the geogrid file specified above.
updateTexture = False

#### Category to fill in for soil class if a cell is water in the soil layer but NOT water in the land cover layer:
# If the script encounters a cell that is classified as land in the land use field (LU_INDEX) but is
# classified as a water soil type, it will replace the soil type with the value you specify below.
# If updateTexture is TRUE, these chages will be propagated to the geogrid. If not, they are just
# used in parameter assignment.
# Ideally there are not very many of these, so you can simply choose the most common soil type in
# your domain. Alternatively, you can set to a "bad" value (e.g., -8888) to see how many of these
# conflicts there are. If you do this DO NOT RUN THE MODEL WITH THESE BAD VALUES. Instead, fix them
# manually with a neighbor fill or similar fill algorithm.
soilFillVal = 3

# Output files to create:
# IMPORTANT: The netcdf files below will be overwritten if they exist!
slpropFile = 'soil_properties.nc'
hyd2dFile = 'hydro2dtbl.nc'

#### Hard-wire urban soil properties in hydro 2d table?
# Some soil parameters are hard-coded to preset values in NoahMP and WRF-Hydro for urban land cover cells.
# If you want to show these in your hyd2dFile parameter file, set this to TRUE. If you want to show
# default parameters, set to FALSE. There should be no answer differences either way.
setUrban = False

#######################################################
# Do not update below this line.
#######################################################

### Number of soil layers (e.g., 4)
# This number should be consistent with the nsoil in the geogrid IF you choose the updateTexture option.
nsoil = 4

# User-specified soil layer thickness (dzs). Also used to calculate depth of the center of each layer (zs)
lyrHt = [0.1, 0.3, 0.6, 1.0]  # Soil layer thickness top layer to bottom (m)

# Data model for output netCDF data
outNCType = 'NETCDF4'

# Default WRF-Hydro nodata value
fillValue = -9999.0

# Input parameter table names and MPTABLE parsing
soilParamFile = 'SOILPARM.TBL'
mpParamFile = 'MPTABLE.TBL'
genParamFile = 'GENPARM.TBL'
if landClass == "USGS":
    hydParamFile = 'HYDRO.TBL'
    mp_params = 'noahmp_usgs_parameters'
elif landClass == "MODIS":
    hydParamFile = 'HYDRO_MODIS.TBL'
    mp_params = 'noahmp_modis_parameters'

# Map variable names from netCDF (lower case, keys) to parameter column headings
# from .TBL files (upper case, values)
nameLookupSoil = dict(smcref="REFSMC",
                      dwsat="SATDW",
                      smcdry="DRYSMC",
                      smcwlt="WLTSMC",
                      bexp="BB",
                      dksat="SATDK",
                      psisat="SATPSI",
                      quartz="QTZ",
                      refdk="REFDK",
                      refkdt="REFKDT",
                      slope="SLOPE",
                      smcmax="MAXSMC",
                      cwpvt="CWPVT",
                      vcmx25="VCMX25",
                      mp="MP",
                      hvt="HVT",
                      mfsno="MFSNO",
                      AXAJ="AXAJ",
                      BXAJ="BXAJ",
                      XXAJ="XXAJ")

# 3D variables in soil_properties file
var3d = ["smcref",
         "dwsat",
         "smcdry",
         "smcwlt",
         "bexp",
         "dksat",
         "psisat",
         "quartz",
         "smcmax"]

# All possible soilParam columns
soilparm_columns = ['solID',
                    'BB',
                    'DRYSMC',
                    'F11',
                    'MAXSMC',
                    'REFSMC',
                    'SATPSI',
                    'SATDK',
                    'SATDW',
                    'WLTSMC',
                    'QTZ',
                    'AXAJ',
                    'BXAJ',
                    'XXAJ',
                    'solName']

## Assumes that each table starts with 'Soil Parameters' and the next line is the name of the parameter set
SOILPARM_Start = 'Soil Parameters'

# Additional parameters to be added from MPTABLE 'noahmp_global_parameters' table
add_vars = ['ssi', 'snowretfac', 'tau0', 'rsurfsnow', 'scamax', 'rsurfexp']

# Default global values for snow variables in case they are not found in MPTABLE
add_var_defaults = {'ssi': 0.03,
                    'snowretfac': 0.00005,
                    'tau0': 1000000.0,
                    'rsurfsnow': 50.0,
                    'scamax': 1.0,
                    'rsurfexp': 5.0}

# Mapping between output variable (key) and input parameter name from MPTABLE (value)
global_mp_dict = {'ssi': 'ssi',
                  'snowretfac': 'snow_ret_fac',
                  'tau0': 'tau0',
                  'rsurfsnow': 'rsurf_snow',
                  'scamax': 'scamax',
                  'rsurfexp': 'rsurf_exp'}

# Hydro 2D Table variable names from netCDF (lower case, keys) to parameter column
# headings from .TBL files (upper case, values)
nameLookupHyd = dict(SMCMAX1="smcmax",
                     SMCREF1="smcref",
                     SMCWLT1="smcwlt",
                     OV_ROUGH2D="OV_ROUGH2D",
                     LKSAT="dksat",
                     NEXP="NEXP")

# Soil layers that will be averaged over when calculating pedo-transfer function parameters.
hydTagList = [0, 1, 2, 3]

# Output netCDF file structure
soilDimName = 'soil_layers_stag'
dims3D = tuple(['Time', 'soil_layers_stag', 'south_north', 'west_east'])
dims2D = tuple(['Time', 'south_north', 'west_east'])

# Pedo-transfer function parameter range limits ([min, max]; set to NULL if you don't want to apply limits):
theta_1500_rng = [0.03, 0.45]
theta_33_rng = [0.07, 0.56]
theta_s33_rng = [0.01, 0.50]
psi_e_rng = [0.1, 30.0]


# --- End Global Variables --- #

# --- Functions --- #

def is_valid_file(parser, arg):
    # https://stackoverflow.com/questions/11540854/file-as-command-line-argument-for-argparse-error-message-if-argument-is-not-va
    if not os.path.exists(arg):
        parser.error(f"The file {arg} does not exist!")
    else:
        return str(arg)


def array_replace(from_values, to_values, array):
    '''
    This function will perform a fast array replacement based on arrays of
    from and to values.
    '''
    sort_idx = np.argsort(from_values)
    idx = np.searchsorted(from_values, array, sorter=sort_idx)
    out = to_values[sort_idx][idx]
    return out


def lyrAvg(arr, hts, inds):
    '''
    This function will calculate a mean from a multidimensional array given
    the list of indexes to average.
    This seems to more heavily weight the inputs by the depth... ?!
    '''
    wtval = 0
    wtsum = 0
    for i in inds:
        wtval += arr[i] * hts[i]
        wtsum += hts[i]
    return wtval / wtsum


def ApplyPedo(sand, clay, orgm=0.0):
    # Calcs
    theta_1500t = (-0.024 * sand) + (0.487 * clay) + (0.006 * orgm) + (0.005 * (sand * orgm)) - (
                0.013 * (clay * orgm)) + (0.068 * sand * clay) + 0.031
    theta_1500 = theta_1500t + ((0.14 * theta_1500t) - 0.02)
    theta_33t = (-0.251 * sand) + (0.195 * clay) + (0.011 * orgm) + (0.006 * (sand * orgm)) - (
                0.027 * (clay * orgm)) + (0.452 * sand * clay) + 0.299
    theta_33 = theta_33t + ((1.283 * theta_33t * theta_33t) - (0.374 * theta_33t) - 0.015)
    theta_s33t = (0.278 * sand) + (0.034 * clay) + (0.022 * orgm) - (0.018 * (sand * orgm)) - (
                0.027 * (clay * orgm)) - (0.584 * sand * clay) + 0.078
    theta_s33 = theta_s33t + ((0.636 * theta_s33t) - 0.107)
    psi_et = (-21.67 * sand) - (27.93 * clay) - (81.97 * theta_s33) + (71.12 * (sand * theta_s33)) + (
                8.29 * (clay * theta_s33)) + (14.05 * sand * clay) + 27.16
    psi_e = psi_et + ((0.02 * psi_et * psi_et) - (0.113 * psi_et) - 0.7)

    # Clip array values to realistic bounds
    if theta_1500_rng:
        theta_1500 = np.clip(theta_1500, theta_1500_rng[0], theta_1500_rng[1])
    if theta_33_rng:
        theta_33 = np.clip(theta_33, theta_33_rng[0], theta_33_rng[0])
    if theta_s33_rng:
        theta_s33 = np.clip(theta_s33, theta_s33_rng[0], theta_s33_rng[0])
    if psi_e_rng:
        psi_e = np.clip(psi_e, psi_e_rng[0], psi_e_rng[1])

    # Almost final param arrays
    smcwlt_3d = theta_1500
    smcref_3d = theta_33
    smcmax_3d = theta_33 + theta_s33 - (0.097 * sand) + 0.043
    bexp_3d = 3.816712826 / (np.log(theta_33) - np.log(theta_1500))
    psisat_3d = psi_e
    dksat_3d = 1930.0 * (smcmax_3d - theta_33) ** (3.0 - (1.0 / bexp_3d))
    quartz_3d = sand

    # Units conversion
    psisat_3d = 0.101997 * psisat_3d  # convert kpa to m
    dksat_3d = dksat_3d / 3600000.0  # convert mm/h to m/s
    dwsat_3d = (dksat_3d * psisat_3d * bexp_3d) / smcmax_3d  # units should be m*m/s
    smcdry_3d = smcwlt_3d
    outDict = dict(smcwlt=smcwlt_3d, smcref=smcref_3d, smcmax=smcmax_3d, bexp=bexp_3d,
                   psisat=psisat_3d, dksat=dksat_3d, quartz=quartz_3d, dwsat=dwsat_3d, smcdry=smcdry_3d)
    return outDict


def demote_dtype(in_df, in_dtype='float64', out_dtype='float32'):
    '''
    Alter dtype of a pandas DataFrame
    '''
    for column in in_df.columns:
        if in_df[column].dtype == in_dtype:
            in_df[column] = in_df[column].astype(out_dtype)
    return in_df


def obtain_soilparams(in_file, table_name='STAS'):
    '''
    Create pandas DataFrame tables from each parameter set in the soil parameter
    files. Assumes that each table starts with 'Soil Parameters' and the next
    line is the name of the parameter set.
    '''

    # Read entire table into a list of lines
    with open(in_file, 'r') as td:
        lines = td.readlines()

    # Find where each table starts and ends
    table_start_lines = [n for n, item in enumerate(lines) if item.startswith(SOILPARM_Start)]
    num_tables = len(table_start_lines)
    print(f'  Found {num_tables} Soil Parameter tables in {os.path.basename(in_file)}')
    table_start_lines += [len(lines)]

    # Build a dictionary of tables
    for n in range(num_tables):
        in_table_name = lines[table_start_lines[n] + 1].strip()
        if in_table_name == table_name:
            print(f"  Found table: {table_name}")

            # Find header inside of quotes and modify header line
            header = re.findall(r"'(.*?)'", lines[table_start_lines[n] + 2].strip(), re.DOTALL)
            header = re.split(r'\s{2,}', header[0])
            header.insert(0, 'solID')
            header = ['solName' if column == '' else column for column in header]
            header += [item for item in soilparm_columns if
                       item not in header]  # Add columns to the end that may be missing

            # Identify the body lines for the sub-table
            body = [line for line in lines[table_start_lines[n] + 3:table_start_lines[n + 1]]]

            # Generate the dataframe and modify data types
            df = pd.read_table(StringIO('\n'.join(body)), sep=',', header=None, names=header)
            df = demote_dtype(df, in_dtype='float64', out_dtype='float32')
            df = demote_dtype(df, in_dtype='i8', out_dtype='i4')
            return df


def obtain_MPparams(in_file, table_name='noahmp_modis_parameters'):
    '''
    Create pandas DataFrame tables from each parameter set in the MP parameter files.
    '''
    # Read the parameter table as a FORTRAN namelist
    namelist_patch = f90nml.read(in_file)

    # Convert top-level namelist to dictionary
    namelist_dict = namelist_patch.todict()
    num_tables = len(namelist_dict)
    if table_name in namelist_dict.keys():
        print(f"  Found table: {table_name}")

        # Read selected table
        table = namelist_dict[table_name]

        # Read the information in the tablef
        val_dict = {}
        list_dict = {}
        for key, val in list(table.items()):
            if type(val) == int:
                val_dict[key] = val
            elif type(val) == list:
                list_dict[key] = val
            elif type(val) == float:
                list_dict[key] = val

        # Handle different returns
        if table_name == 'noahmp_global_parameters':
            return list_dict
        else:

            # set column names to 1...n
            columns = range(1, max([len(x) for x in list_dict.values()]) + 1)

            # Create the dataframe
            df = pd.DataFrame.from_dict(list_dict)
            df.index = columns
            return df


def obtain_GENparams(in_file):
    '''
    Create dictionary of GENPARM parameters. Assumes that all parameters end with '_DATA'.
    '''

    # Read entire table into a list of lines
    GPdict = {}
    with open(in_file, 'r') as td:
        lines = td.readlines()
        lines = [line.strip() for line in lines]

        # Find where the parameters start and end in the file
        param_start_lines = [n for n, item in enumerate(lines) if item.endswith('_DATA')]
        num_params = len(param_start_lines)
        param_start_lines += [len(lines)]
        print(f'  Found {num_params} parameters in {os.path.basename(in_file)}')

        for n in range(num_params):
            param_name = lines[param_start_lines[n]]
            GPdict[param_name] = [float(line) for line in lines[param_start_lines[n] + 1:param_start_lines[n + 1]]]
    return GPdict


def obtain_HYDROparams(in_file):
    # Read entire table into a list of lines
    with open(in_file, 'r') as td:
        lines = td.readlines()
        lines = [line.strip() for line in lines]

    # Find the roughness parameters for landcover types
    num_SFC_ROUGH = int(lines[0].split(' ')[0])
    SFC_ROUGH_startline = [n for n, line in enumerate(lines) if line.startswith('SFC_ROUGH')][0]
    SFC_ROUGH_endline = SFC_ROUGH_startline + 1 + num_SFC_ROUGH
    SFC_ROUGH_vals = [line for line in lines[SFC_ROUGH_startline + 1:SFC_ROUGH_endline]]

    # Create the surface roughness parameter dataframe
    sfcRough_df = pd.read_table(StringIO('\n'.join(SFC_ROUGH_vals)), sep=',', header=None,
                                names=('OV_ROUGH2D', 'descrip'))
    print(f'  Found {len(sfcRough_df)} SFC_ROUGH parameters in {os.path.basename(in_file)}')

    # Find the soil hydro parameters
    soil_param_start = SFC_ROUGH_endline
    num_soil_params = int(lines[soil_param_start].split(',')[0])

    # Find header for soil parameters
    header = re.split(r'\s{2,}', lines[soil_param_start + 1])
    header = ['solName' if column == "'" else column for column in header]

    # Identify the body lines for the sub-table
    body = [line for line in lines[soil_param_start + 2:soil_param_start + 2 + num_soil_params]]

    # Generate the dataframe and modify data types
    hydro_df = pd.read_table(StringIO('\n'.join(body)), sep=',', header=None, names=header)
    hydro_df = demote_dtype(hydro_df, in_dtype='float64', out_dtype='float32')
    hydro_df = demote_dtype(hydro_df, in_dtype='i8', out_dtype='i4')
    return hydro_df, sfcRough_df


def main_soilProp(geoFile,
                  paramDir,
                  outDir,
                  landClass=landClass,
                  soilFillVal=soilFillVal,
                  setUrban=setUrban,
                  soilparm_set=soilparm_set):
    '''This function will process the input geogrid file into the appropriate
    2D soil property and hydro prarameter output variables.

    Inputs:
        geoFile  - WRF WPS GEOGRID file (netCDF).
        paramDir - The directory containing WRF-Hydro parameter files.
        outDir   - Directory to store outputs
    Outputs:
        {slpropFile} - 2-Dimensional soil properties file (netCDF)
        {hyd2dFile}  - 2-Dimensional hydro parameter file (netCDF)
    '''
    tic1 = time.time()

    # Setup input and output files
    slpropF = os.path.join(outDir, slpropFile)
    hydFile = os.path.join(outDir, hyd2dFile)
    soilParamF = os.path.join(paramDir, soilParamFile)
    mpParamF = os.path.join(paramDir, mpParamFile)
    genParamF = os.path.join(paramDir, genParamFile)
    hydParamF = os.path.join(paramDir, hydParamFile)
    for inFile in [geoFile, soilParamF, mpParamF, genParamF]:
        if not os.path.isfile(inFile):
            print(f'Expected input file "{inFile}" not found.\nExiting...')
            raise SystemExit
    print('  Creating initial soil_properties.nc file from the input geogrid file')
    print(f'    Input:  {geoFile}')
    print(f'    Output: {slpropF}')
    print(f'    Output: {hydFile}')
    print(f'  Keyword arguments:')
    print(f'    landClass:    {landClass}')
    print(f'    soilFillVal:  {soilFillVal}')
    print(f'    setUrban:     {setUrban}')
    print(f'    soilparm_set: {soilparm_set}')

    # Instantiate read and write objects on netCDF files
    if updateTexture:
        new_geoFile = os.path.join(outDir, os.path.basename(geoFile).replace('.nc', '_texture.nc'))
        shutil.copyfile(geoFile, new_geoFile)
        print(f'  Copied input geogrild file:\n  \t{geoFile}\n  to:\n  \t{new_geoFile}')
        rootgrp_geo = netCDF4.Dataset(new_geoFile, 'r+')  # Read/Write object on input GEOGRID file
    else:
        rootgrp_geo = netCDF4.Dataset(geoFile, 'r')  # Read object on input GEOGRID file
    rootgrp_SP = netCDF4.Dataset(slpropF, 'w', format=outNCType)  # Write object to create output file

    # Copy dimensions from GEOGRID file
    for dimname, dim in rootgrp_geo.dimensions.items():
        if dimname in dims2D:
            if dimname == 'Time':
                dimlen = None  # Make Time an unlimited dimension
            else:
                dimlen = len(dim)
            rootgrp_SP.createDimension(dimname, dimlen)  # Copy dimensions from the GEOGRID file
    soildim = rootgrp_SP.createDimension(soilDimName, nsoil)  # Add soil_layers_stag dimension
    del dimname, dim, soildim

    # Populate initial file with variables to keep from the input GEOGRID file
    for varname in nameLookupSoil:
        if varname in var3d:
            rootgrp_SP.createVariable(varname, 'f4', dims3D, fill_value=fillValue)
        else:
            rootgrp_SP.createVariable(varname, 'f4', dims2D, fill_value=fillValue)
    del varname

    # Global Attributes - copy all from GEOGRID file to soil_properties
    ncatts = {key: val for key, val in rootgrp_geo.__dict__.items()}
    ncatts['Source_Software'] = f'WRF-Hydro {sys.argv[0]} script (Python) v1.0'
    ncatts['creation_time'] = f'Created {time.ctime()}'
    rootgrp_SP.setncatts(ncatts)

    # Create new hydro2d file with fill values
    dims2D_hyd = [item for item in dims2D if item != 'Time']
    rootgrp_hyd = netCDF4.Dataset(hydFile, 'w', format=outNCType)  # Write object to create output file
    for dimname, dim in rootgrp_geo.dimensions.items():
        if dimname in dims2D_hyd:
            rootgrp_hyd.createDimension(dimname, len(dim))  # Copy dimensions from the GEOGRID file
    for varname in nameLookupHyd.keys():
        rootgrp_hyd.createVariable(varname, 'f4', dims2D_hyd, fill_value=fillValue)
    rootgrp_hyd.setncatts(ncatts)
    del dimname, dim, varname

    # --- Read parameter tables --- #

    # --- Read Soil parameter table --- #

    print(f'Reading {soilParamF}')
    soilTblDf = obtain_soilparams(soilParamF, table_name=soilparm_set)

    # --- End Read Soil parameter table --- #

    # --- Read MP parameter table --- #

    print(f'Reading {mpParamF}')
    mpTblDf = obtain_MPparams(mpParamF, table_name=mp_params)
    mpTblDict = obtain_MPparams(mpParamF, table_name='noahmp_global_parameters')

    # --- End Read MP parameter table --- #

    # --- Read GENPARM parameter table --- #

    print(f'Reading {genParamF}')
    GPdict = obtain_GENparams(genParamF)
    genTab = {'SLOPE': GPdict['SLOPE_DATA'][1],
              'REFKDT': GPdict['REFKDT_DATA'][0],
              'REFDK': GPdict['REFDK_DATA'][0]}

    # --- End Read GENPARM parameter table --- #

    # --- Read HYDRO parameter table --- #

    print(f'Reading {hydParamF}')
    hydroTblDf, sfcRoughDf = obtain_HYDROparams(hydParamF)

    # --- End Read HYDRO parameter table --- #

    # Get 2D fields and global attribute values from input GEOGRID file
    vegmap = rootgrp_geo.variables['LU_INDEX'][0]
    solmap = rootgrp_geo.variables['SCT_DOM'][0]
    maxSoilClass = len(rootgrp_geo.dimensions['soil_cat'])
    vegWater = rootgrp_geo.ISWATER
    vegLake = rootgrp_geo.ISLAKE
    soilWater = rootgrp_geo.ISOILWATER
    vegUrban = rootgrp_geo.ISURBAN
    print(f'Geogrid attributes: vegWater={vegWater} soilWater={soilWater}, maxSoilClass={maxSoilClass}')

    if useSoilComp and 'SOILCOMP' in rootgrp_geo.variables.keys():
        print('Pulling in soil composition data')
        soilc = rootgrp_geo.variables['SOILCOMP'][0]
        soilcFrac = (soilc / 100.0).astype('f8')  # Convert percent to fraction and make Double
        sandFrac3D = soilcFrac[0:4]  # Sand fraction for layers 1-4
        clayFrac3D = soilcFrac[4:8]  # Clay fraction for layers 1-4
        orgm3D = 0.0

        # Create water mask from layers with all 0s
        soilWaterMsk = np.sum(soilc, axis=0)
        soilWaterMsk[soilWaterMsk == 0] = fillValue
        soilWaterMsk[soilWaterMsk > 0] = 0

        # Create layer means
        layerMeans = {}  # Empty dictionary to store data
        layerMeans[f'sand'] = lyrAvg(sandFrac3D, lyrHt, hydTagList)
        layerMeans[f'clay'] = lyrAvg(clayFrac3D, lyrHt, hydTagList)

        # Apply pedotransfer functions
        print('Applying pedotransfer functions')
        soillist3D = ApplyPedo(sandFrac3D, clayFrac3D, orgm3D)
        pedoXfer = {}  # Empty dictionary to store transformed parameters
        pedoXfer[f'soillist'] = ApplyPedo(layerMeans[f'sand'], layerMeans[f'clay'])
        del soilcFrac, sandFrac3D, clayFrac3D
    else:
        # Fill all areas in SCT_DOM where the vegetation is water and soil is water with a fill value
        vegmap[vegmap == vegLake] = vegWater
        solmap[np.logical_and(vegmap != vegWater, solmap == soilWater)] = soilFillVal
        solmap[vegmap == vegWater] = soilWater

    paramList = rootgrp_SP.variables.keys()
    print(f'Updating {slpropFile}')
    for param in paramList:
        paramName = nameLookupSoil[param]
        ncvar = rootgrp_SP.variables[param]
        if paramName in soilTblDf.columns.values:
            # Parameter is in the soil table, map to the categories
            print(f'    Updating SOIL parameter: {param} {paramName}')
            if useSoilComp and 'SOILCOMP' in rootgrp_geo.variables.keys():
                pnew = soillist3D[param][0]

                # Create a water mask
                pnew[vegmap == vegWater] = 0
                # pnew[vegmap!=vegWater] = fillValue
                # pnew[pnew>fillValue] = soilWater
            else:
                pnew = solmap.copy()
                paramVar = ncvar[:]

                # Create 2D replacement matrix
                replaceArr = np.array([soilTblDf.solID.values, soilTblDf[paramName]])
                pnew = array_replace(replaceArr[0, :], replaceArr[1, :], pnew).astype('f4')

                # Build mask array to reclassify out of range values
                maskArr = np.full(pnew.shape, True, dtype=bool)
                maskArr = array_replace(replaceArr[0, :], np.full(soilTblDf.solID.values.shape, False, dtype=bool),
                                        maskArr)
                pnew[maskArr] = fillValue

            # Fill in output netCDF variable with this value
            if 'soil_layers_stag' in ncvar.dimensions:
                dimsize = len(rootgrp_SP.dimensions[ncvar.dimensions[-3]])  # Either 'Time' or 'soil_layers_stag'

                # Write output to NetCDF
                ncvar[0] = np.repeat(pnew[np.newaxis, :, :], dimsize, axis=0)
            else:
                # Write output to NetCDF
                ncvar[0] = pnew

        elif paramName.lower() in mpTblDf.columns:

            # Parameter is in the table, map to the categories.
            print(f'    Updating MP parameter: {param} {paramName}')
            paramVals = mpTblDf[paramName.lower()]
            paramCols = mpTblDf.index

            # Create 2D replacement matrix
            replaceArr = np.array([paramCols.tolist(), paramVals.tolist()])
            pnew = array_replace(replaceArr[0, :], replaceArr[1, :], vegmap).astype('f4')

            # Build mask array to reclassify out of range values
            maskArr = np.full(pnew.shape, True, dtype=bool)
            maskArr = array_replace(replaceArr[0, :], np.full(paramVals.shape, False, dtype=bool), maskArr)
            pnew[maskArr] = fillValue

            # Write output to NetCDF
            ncvar[0] = pnew

        elif paramName in genTab:
            # Parameter is in the general parameter table
            print(f'    Updating GEN parameters: {param} {paramName} with {genTab[paramName]}')

            # Write output to NetCDF
            ncvar[:, :, :] = genTab[paramName]

    # Add additional noahmp global parameters
    print(f'    Adding additional noahmp global parameters to {slpropFile}')
    for add_var in add_vars:
        if global_mp_dict[add_var] in mpTblDict:
            add_val = mpTblDict.get(global_mp_dict[add_var], fillValue)
            print(f'      Adding global parameter {add_var}={add_val} to grid.')
            rootgrp_SP.createVariable(add_var, 'f4', dims2D, fill_value=fillValue)
            rootgrp_SP[add_var][:] = add_val
        else:
            print(f'      Could not find parameter {global_mp_dict[add_var]} in {mpParamFile}.')
            print(f'        Using default value of {add_var}={add_var_defaults[add_var]}.')
            rootgrp_SP.createVariable(add_var, 'f4', dims2D, fill_value=fillValue)
            rootgrp_SP[add_var][:] = add_var_defaults[add_var]

    rootgrp_SP.close()
    del pnew

    # Options to update the input Geogrid file texture classes
    if updateTexture:
        print('Updating texture classes')
        lyrDict = dict(top=dict(indx=0, splitlyr="SOILCTOP", mglyr="SCT_DOM"),
                       bot=dict(indx=-1, splitlyr="SOILCBOT", mglyr="SCB_DOM"))

        for layer, layerOpts in lyrDict.items():

            # Fill the 2D variable (SC*_DOM)
            inLayerName = layerOpts['mglyr']
            splitlyr = layerOpts['splitlyr']
            i = layerOpts['indx']
            ncvar = rootgrp_geo.variables[inLayerName]
            ncvar2 = rootgrp_geo.variables[splitlyr]
            soil_texture = np.full(ncvar[0].shape, fillValue, dtype='f4')

            if useSoilComp and 'SOILCOMP' in rootgrp_geo.variables.keys():
                print(f'  Using SOILCOMP variable to update soils for layer {layerOpts}.')
                # SOILCOMP variable has vertical dimension with sand fraction layers
                # then clay fraction layers. Silt can be calculated as the remainder.
                sand = soilc[0:nsoil]  # Sand fraction for layers 1-4
                clay = soilc[nsoil:]  # Clay fraction for layers 1-4
                silt = 100.0 - sand - clay  # Silt is everything except sand and clay

                # Bin into texture classes
                soil_texture[(silt[i] + 1.5 * clay[i]) < 15] = 1
                soil_texture[((silt[i] + 1.5 * clay[i]) >= 15) & (silt[i] + 2 * clay[i] < 30)] = 2
                soil_texture[((clay[i] >= 7) & (clay[i] < 20)) & (sand[i] > 52) & (((silt[i] + 2 * clay[i]) >= 30) | (
                            (clay[i] < 7) & (silt[i] < 50) & ((silt[i] + 2 * clay[i]) >= 30)))] = 3
                soil_texture[
                    ((clay[i] >= 7) & (clay[i] < 27)) & ((silt[i] >= 28) & (silt[i] < 50)) & (sand[i] <= 52)] = 6
                soil_texture[((silt[i] >= 50) & ((clay[i] >= 12) & (clay[i] < 27))) | (
                            ((silt[i] >= 50) & (silt[i] < 80)) & (clay[i] < 12))] = 4
                soil_texture[((silt[i] >= 80) & (clay[i] < 12))] = 5
                soil_texture[((clay[i] >= 20) & (clay[i] < 35)) & ((silt[i] < 28) & (sand[i] > 45))] = 7
                soil_texture[((clay[i] >= 27) & (clay[i] < 40)) & ((sand[i] > 20) & (sand[i] <= 45))] = 9
                soil_texture[((clay[i] >= 27) & (clay[i] < 40)) & (sand[i] <= 20)] = 8
                soil_texture[(clay[i] >= 35) & (sand[i] > 45)] = 10
                soil_texture[(clay[i] >= 40) & (silt[i] >= 40)] = 11
                soil_texture[(clay[i] >= 40) & (sand[i] <= 45) & (silt[i] < 40)] = 12

                # Apply water masks
                soil_texture[soilWaterMsk == fillValue] = soilFillVal
                soil_texture[vegmap == vegWater] = soilWater
                del soilc, sand, clay, silt
            else:
                print('  SOILCOMP variable not in the input GEOGRID file.')
                print('    Updating Soil Texture in GEOGRID file.')
                # soil_texture = solmap.copy()                       # Will only be SCT_DOM!
                soil_texture = ncvar[0]
                soil_texture[np.logical_and(vegmap != vegWater, soil_texture == soilWater)] = soilFillVal
                soil_texture[vegmap == vegWater] = soilWater

            # Write output to NetCDF
            ncvar[0] = soil_texture

            # Calculate and place split layer (assumes 100% for specified class)
            for stype in range(1, maxSoilClass + 1):
                tmp = soil_texture.copy()
                tmp[soil_texture == stype] = 1
                tmp[soil_texture != stype] = 0
                ncvar2[0, stype - 1] = tmp
            del ncvar, ncvar2, soil_texture, tmp

    # Populate hydro2d file
    print(f'Updating: {hyd2dFile}')
    vegmap = rootgrp_geo.variables['LU_INDEX'][0]
    solmap = rootgrp_geo.variables['SCT_DOM'][0]

    # Loop through params and update
    paramList = rootgrp_hyd.variables
    for param in paramList:
        paramNameHyd = nameLookupHyd.get(param)
        paramNameSoil = nameLookupSoil.get(paramNameHyd)
        print(f'    Processing {param}')
        ncvar = rootgrp_hyd.variables[param]
        if paramNameHyd and (paramNameSoil in soilTblDf.columns.tolist()):
            print(f'      Updating HYDRO soil parameters: {param} {paramNameHyd} {paramNameSoil}')

            pnew = solmap.copy()
            if useSoilComp and 'SOILCOMP' in rootgrp_geo.variables.keys():
                paramVar = pedoXfer[f'soillist'][paramNameHyd]

                # Create a water mask
                pnew[vegmap != vegWater] = fillValue
                pnew[pnew > fillValue] = soilWater
            else:
                paramVar = ncvar[:]

                # Create 2D replacement matrix
                pnew[vegmap == vegWater] = soilWater
                replaceArr = np.array([soilTblDf.solID.values, soilTblDf[paramNameSoil]])
                pnew = array_replace(replaceArr[0, :], replaceArr[1, :], pnew).astype('f4')

            # Added to incorporate pedo-transfer function data
            pnew[pnew < -9998] = paramVar[pnew < -9998]

            # Manually make some changes to urban cells to match hydro code.
            if setUrban:
                for param_to_alter, replaceval in zip(['SMCMAX1', 'SMCREF1', 'SMCWLT1'], [0.45, 0.42, 0.40]):
                    if param == param_to_alter:
                        pnew[np.where((vegmap == vegUrban) & (solmap != soilWater))] = replaceval
                        print(
                            f'        Replaced HYDRO2D parameter {param} with {replaceval} where vegetation is urban but soil is not soilwater')

            # Write output to NetCDF
            ncvar[:, :] = pnew

        elif paramNameHyd in sfcRoughDf:
            print(f'      Updating {paramNameHyd}')
            pnew = vegmap.copy()
            pnew[solmap == soilWater] = vegWater

            # Loop through each vegetation category
            for catTmp in range(len(sfcRoughDf.descrip)):
                pnew[np.where(pnew == int(catTmp + 1))] = float(sfcRoughDf[paramNameHyd][catTmp])

            # Write output to NetCDF
            ncvar[:, :] = pnew
        elif paramNameHyd == "NEXP":
            # Setting this to a global initial value of 1.0
            print(f"Updating HYDRO global parameters: {param}  {paramNameHyd}")
            ncvar[:, :] = 1.0

    rootgrp_hyd.close()
    rootgrp_geo.close()
    del rootgrp_geo, rootgrp_hyd
    print('    Function main_soilProp completed in {0:3.2f} seconds'.format(time.time() - tic1))
    return


# --- End Functions --- #

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
    parser.add_argument("-p",
                        dest="param_dir",
                        type=lambda x: is_valid_file(parser, x),
                        default='./',
                        help="[OPTIONAL] Directory containing SOILPARM.TBL, MPTABLE.TBL, GENPARM.TBL, and HYDRO.TBL. Default is the current directory.")
    parser.add_argument("-o",
                        dest="out_dir",
                        type=lambda x: is_valid_file(parser, x),
                        default='./',
                        help="[OPTIONAL] Directory to write output files to (soil_properties.nc, hydro2dtbl.nc). Default is the current directory.")

    # If no arguments are supplied, print help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    all_defaults = {key: parser.get_default(key) for key in vars(args)}

    # Handle printing to user the default variable name
    print('  Parameter values that have not been altered from script default values:')
    if args.in_Geogrid == all_defaults["in_Geogrid"]:
        print('    Using default basin mask setting: {0}'.format(all_defaults["in_Geogrid"]))
    if args.param_dir == all_defaults["param_dir"]:
        print('    Using default reach-based routing setting: {0}'.format(all_defaults["param_dir"]))
    if args.out_dir == all_defaults["out_dir"]:
        print('    Using default regridding factor: {0}'.format(all_defaults["out_dir"]))

    # This block allows us to continue to check for a valid file path while allowing the script later to avoid a NoneType error.
    geoFile = args.in_Geogrid = os.path.abspath(
        args.in_Geogrid)  # Obtain absolute path for required input geogrid file..
    paramDir = args.param_dir = os.path.abspath(args.param_dir)  # Obtain absolute path for required input directory.
    outDir = args.out_dir = os.path.abspath(args.out_dir)  # Obtain absolute path for required output directory.

    # Print information to screen
    print('  Values that will be used in building this routing stack:')
    print(f'    Input GEOGRID file:  {geoFile}')
    print(f'    Parameter directory: {paramDir}')
    print(f'    Output directory:    {outDir}')

    # Run pre-process
    print('  Running Process: main_soilProp function')
    main_soilProp(geoFile,
                  paramDir,
                  outDir)
    print('Process completed in {0:3.2f} seconds.'.format(time.time() - tic))
