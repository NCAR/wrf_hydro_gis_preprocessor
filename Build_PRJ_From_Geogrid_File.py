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

# Import additional modules
import netCDF4

# Import function library into namespace. Must exist in same directory as this script.
from wrfhydro_functions import WRF_Hydro_Grid
print('Script initiated at {0}'.format(time.ctime()))

# Global Variables

# Input and output files and directories
inGeogrid = r"C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\geo_em.d03.20200327_snow.trim.nc"
out_dir = r'C:\Users\ksampson\Desktop\NWM\NWM_Alaska\HRRR_AK\NWM\FOSS_Domain'
outPRJ = os.path.join(out_dir, os.path.basename(inGeogrid).replace('.nc', '.prj'))

# Script options
buildPRJ = True                                                                 # Switch for building output Esri Projection File

# Main Codeblock
if __name__ == '__main__':
    tic = time.time()

    rootgrp = netCDF4.Dataset(inGeogrid, 'r')                                   # Establish an object for reading the input NetCDF file
    coarse_grid = WRF_Hydro_Grid(rootgrp)                                  # Instantiate a grid object
    rootgrp.close()
    del rootgrp
    print('    Created projection definition from input NetCDF GEOGRID file.')

    # Build an ESRI style projection file
    if buildPRJ:
        projEsri = coarse_grid.proj.Clone()                                     # Copy the SRS
        projEsri.MorphToESRI()                                                  # Alter the projection to Esri's representation of a coordinate system
        file = open(outPRJ, 'w')
        file.write(projEsri.ExportToWkt())
        file.close()
        print('    Created ESRI Projection file: {0}'.format(outPRJ))
    del coarse_grid
    print('Process complted in {0:3.2f} seconds.'.format(time.time()-tic))
