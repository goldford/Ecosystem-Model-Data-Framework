import netCDF4 as nc
import numpy as np
import pandas as pd
import os

from GO_helpers import buildSortableString, saveASCFile, getDataFrame # 2023

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO last edited by GO 2023-05-10
#
# purpose: convert wind files in monthly NC to ASC format for ECOSPACE (with header)
#
# log: 
#    GO - removed functions, put in GO_helpers.py, updated for RDRS
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# /////////////////// MISC ///////////////////////////
startYear = 1980
endyear = startYear + 39
sigdigfmt = '%0.1f' # sig figs

#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39
ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"


# /////////////////// PATHS ///////////////////////////
rootDir = "/project/6006412/goldford/ECOSPACE/DATA/"
greigDir = "/project/6006412/goldford/ECOSPACE/"
saveTemplate = "RDRS21_NEMOgrid_{}_{}_{}.asc"
saveDir = "/project/6006412/goldford/ECOSPACE/DATA/wind_monthly_ASC/"
inputsSubDir = "wind_monthly/"
depthSubDir = "DATA/"
if (os.path.isdir(saveDir) != True):
  os.mkdir(saveDir)


# /////////////////// VARS ///////////////////////////
lstVars = ['u10m','v10m'];
#key - match the variable names in lstVars
#dctFileTemplates = {'msdwswrf':'ERA5_NEMOgrid_light_monthly_{}.nc'}
dctFileTemplates = {'u10m':'RDRS21_NEMOgrid_wind_monthly_{}.nc','v10m':'RDRS21_NEMOgrid_wind_monthly_{}.nc'}

#Unit Conversions - mult for Ecospace parameterization
dctUnitConvert = {'u10m':1, 'v10m':1}

# /////////////////// LAND MASK ///////////////////////////
ecospacegrid_f = os.path.join(greigDir, depthSubDir, "ecospacedepthgrid.asc")
df_dep = getDataFrame(ecospacegrid_f,"-9999.00000000")
dfLandMask = df_dep == 0 #mask is where elev =0


for var in lstVars:
    for iyear in range(startYear, endyear):
    
    # input format: one file per year with daily data
        inputTemplate = dctFileTemplates[var]
        filename = inputTemplate.format(iyear) #fixed GO 20210917
        fullpath = os.path.join(rootDir, inputsSubDir,filename)
        
        
        if(os.path.isfile(fullpath)):
            ds = nc.Dataset(fullpath)
            vars = ds.variables[var]
        else:
            print("Could not find NC file. " + fullpath)
            break

        for imon in range(1,13):
            m = buildSortableString(imon,2)    

            monthlydata = vars[imon-1]

            #Save the monthly means to an .asc file
            savepath = os.path.join(saveDir, var) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)

            savefn = saveTemplate.format(var,iyear, m)
            saveASCFile(os.path.join(savepath, savefn) ,monthlydata,bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe, sigdigfmt, ASCheader, None, dfLandMask)
            print("saved-->", os.path.join(savepath, savefn))