#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by G Oldford, 
# last edited 2024-05-23 -GO
# total revision by GO 2022-08-24
#
# Purpose: convert pre-processed NC files (mean monthly NEMO results in monthly NetCDFs) to ASC for ECOSPACE
# Input: monthly NC files
# Output: monthly ASC files
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

import netCDF4 as nc
import numpy as np
import pandas as pd
import os
from calendar import monthrange

def buildSortableString(number, nZeros): 
    newstring = str(number)
    while (len(newstring) < nZeros) :
        tmpstr = "0" + newstring
        newstring = tmpstr
    return newstring    

# GO2021-12-22 fmt= has big effect on ASC file size, added argument for this
def saveASCFile(filename, data, sigdigfmt):
    #Greig's bounds clipping from "Basemap Converter (Py3).ipynb"
    #I don't know why this works but it does... -JB
    # magic. -GO
    
    trimmedData = []
    i =  bottomleft_row_ewe 
    while i < upperleft_row_ewe:
        trimmedData.append(data[i,upperleft_col_ewe:])
        i += 1
    trimmedData.reverse()

    np.savetxt(filename,trimmedData,fmt=sigdigfmt, delimiter=" ", comments='',header=ASCheader)


# /////////////////////// HEADER //////////////////////////
#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39

ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"

# //////////////////////// PATHS //////////////////////////
greigDir = "/project/6006412/goldford/ECOSPACE/"

NEMOMonthlyNCSubDir = "DATA/SS1500-RUN216/NEMO_monthly_NC/"
modelrun = "216"

startYear = 1979
endyear = 2019
sigdigfmt = '%0.1f' #GO2021-12-23 sets decimal places in out ASC e.g 0.2f is two decimals
saveTemplate = "SalishSea1500-RUN{}_{}_{}_{}.asc"
outputsSubDir = "DATA/SS1500-RUN{}/NEMO_prepped_as_ASC"
outdir = os.path.join(greigDir,outputsSubDir.format(modelrun))
if (os.path.isdir(outdir) != True):
  os.mkdir(outdir)


# //////////////////////// VAR DEF ////////////////////////
#lstVars_T_3D_vmean = ['votemper','vosaline','vozocrtx','vomecrty','vovecrtz'];
lstVars_T_3D_vmean = ['votemper','vosaline','vozocrtx','vomecrty'];
lstVars_T_2D = ['mldkz5', 'mldr10_1'];
lstVars_U = ['vozocrtx'];
lstVars_V = ['vomecrty'];
lstVars_W = ['vovecrtz'];

# variable name conversion
dctVarnameToEwEName_T_3D_vertmean = {
'votemper':'TempVertMean',
'vosaline':'SalinVertMean',
'vozocrtx':'xvelocVertMean',
'vomecrty':'yvelocVertMean',
'vovecrtz':'zvelocVertMean'}
dctVarnameToEwEName_T_3D_vmean10m = {
'votemper':'TempVertMean10m',
'vosaline':'SalinVertMean10m',
'vozocrtx':'xvelocMean10m',
'vomecrty':'yvelocMean10m',
'vovecrtz':'zvelocMean10m'}
dctVarnameToEwEName_T_3D_vmean4m = {
'votemper':'TempVertMean4m',
'vosaline':'SalinVertMean4m',
'vozocrtx':'xvelocMean4m',
'vomecrty':'yvelocMean4m',
'vovecrtz':'zvelocMean4m'}

dctVarnameToEwEName_T_2D = {'mldkz5':'MixingTurboZ','mldr10_1':'MixedLyrZ'}
dctVarnameToEwEName_U = {'vozocrtx':'xvelocMean10m'}
dctVarnameToEwEName_V = {'vomecrty':'yvelocMean10m'}
dctVarnameToEwEName_W = {'vozocrtz':'zvelocMean10m'}

# //////////////// FILE NAME TEMPLATES //////////////////////////
dctFileTemplates_T_3D_vertmean = {
'votemper':'SalishSea1500-RUN{}_MonthlyMean_grid_T_vertmean_y{}m{}.nc', 
'vosaline':'SalishSea1500-RUN{}_MonthlyMean_grid_T_vertmean_y{}m{}.nc',
'vozocrtx':'SalishSea1500-RUN{}_MonthlyMean_grid_U_vertmean_y{}m{}.nc',
'vomecrty':'SalishSea1500-RUN{}_MonthlyMean_grid_V_vertmean_y{}m{}.nc',
'vovecrtz':'SalishSea1500-RUN{}_MonthlyMean_grid_W_vertmean_y{}m{}.nc'
}

dctFileTemplates_T_3D_vmean10m = {
'votemper':'SalishSea1500-RUN{}_MonthlyMean_grid_T_avgtop10m_y{}m{}.nc', 
'vosaline':'SalishSea1500-RUN{}_MonthlyMean_grid_T_avgtop10m_y{}m{}.nc',
'vozocrtx':'SalishSea1500-RUN{}_MonthlyMean_grid_U_avgtop10m_y{}m{}.nc',
'vomecrty':'SalishSea1500-RUN{}_MonthlyMean_grid_V_avgtop10m_y{}m{}.nc',
'vovecrtz':'SalishSea1500-RUN{}_MonthlyMean_grid_W_avgtop10m_y{}m{}.nc'
}

dctFileTemplates_T_3D_vmean4m = {
'votemper':'SalishSea1500-RUN{}_MonthlyMean_grid_T_avgtop4m_y{}m{}.nc', 
'vosaline':'SalishSea1500-RUN{}_MonthlyMean_grid_T_avgtop4m_y{}m{}.nc',
'vozocrtx':'SalishSea1500-RUN{}_MonthlyMean_grid_U_avgtop4m_y{}m{}.nc',
'vomecrty':'SalishSea1500-RUN{}_MonthlyMean_grid_V_avgtop4m_y{}m{}.nc',
'vovecrtz':'SalishSea1500-RUN{}_MonthlyMean_grid_W_avgtop4m_y{}m{}.nc'
}

dctFileTemplates_T_2D = {
'mldkz5':'SalishSea1500-RUN{}_MonthlyMean_grid_T_2D_y{}m{}.nc',
'mldr10_1':'SalishSea1500-RUN{}_MonthlyMean_grid_T_2D_y{}m{}.nc'
}

dctFileTemplates_U = {'vozocrtx':'SalishSea1500-RUN{}_MonthlyMean_grid_U_y{}m{}.nc'}
dctFileTemplates_V = {'vozocrty':'SalishSea1500-RUN{}_MonthlyMean_grid_V_y{}m{}.nc'}
dctFileTemplates_W = {'vozocrtz':'SalishSea1500-RUN{}_MonthlyMean_grid_W_y{}m{}.nc'}
dctFileTemplates_T = {'vosaline':'SalishSea1500-RUN{}_MonthlyMean_grid_T_y{}m{}.nc'}


# /////////////// UNIT CONVERSIONS ///////////////////////
#'vovecrtz':1000 upwelling/downwelling (vertical velocity) 
# is m / s but his makes it really hard to parameterize the EwE UI
#so arbitrarily scale them bigger  - JB
#I changed it so multiplier is 100 to convert from m/s (NEMO) to cm/s (EwE) - GO 2022-02
dctUnitConvert_T_3D_vmean = {'votemper':1,'vosaline':1,'vozocrtx':100,'vomecrty':100,'vovecrtz':10000}
dctUnitConvert_T_2D = {'mldkz5':1,'mldr10_1':1}
dctUnitConvert_U = {'vozocrtx':100}
dctUnitConvert_V = {'vozocrty':100}
dctUnitConvert_W = {'vozocrtz':100}
#dctUnitConvert = {'votemper':1,'vosaline':1}
#dctUnitConvert = {'votemper':1,'vosaline':1,'mldr10_1':1,'mldkz5':1,'vovecrtz':100,'vozocrtx':100,'vomecrty':100}

mos = ['01','02','03','04','05','06','07','08','09','10','11','12']
mo_i = 0 # unused? -GO 2022-08

# ////////////////// CONVERT 3D VMEAN /////////////////
# 3D files that are already vert averaged over full water col
for var in lstVars_T_3D_vmean:
  print(var)
  for iyear in range(startYear, endyear):
    for mo in mos:

      inputTemplate = dctFileTemplates_T_3D_vertmean[var]
      datafullpath = os.path.join(greigDir,NEMOMonthlyNCSubDir,inputTemplate.format(modelrun,iyear,mo))
      dsData = nc.Dataset(datafullpath)
      varsData = dsData.variables[var]
      DataMon = varsData[0,:,:]
      
      #unit conversion
      UnitConvert = dctUnitConvert_T_3D_vmean[var]
      DataMon = DataMon * UnitConvert 
  
      #Save the monthly Data to an .asc file     
      EwEName = dctVarnameToEwEName_T_3D_vertmean[var]
      savepath = os.path.join(greigDir,outputsSubDir.format(modelrun),EwEName)

      if (os.path.isdir(savepath) != True):
          os.mkdir(savepath)
  
      savefn = saveTemplate.format(modelrun,EwEName,iyear,mo)
      savefullpath = os.path.join(savepath, savefn)

      print("Saving file " + savefullpath)
      saveASCFile(savefullpath, DataMon, sigdigfmt)


# ////////////////// CONVERT 3D VMEAN TOP 10m /////////////////
# 3D files that are already vert averaged over top 10 m
for var in lstVars_T_3D_vmean:
  print(var)
  for iyear in range(startYear, endyear):
    for mo in mos:

      inputTemplate = dctFileTemplates_T_3D_vmean10m[var]
      datafullpath = os.path.join(greigDir,NEMOMonthlyNCSubDir,inputTemplate.format(modelrun,iyear,mo))
      dsData = nc.Dataset(datafullpath)
      varsData = dsData.variables[var]
      DataMon = varsData[0,:,:]
      
      #unit conversion
      UnitConvert = dctUnitConvert_T_3D_vmean[var]
      DataMon = DataMon * UnitConvert 
  
      #Save the monthly Data to an .asc file     
      EwEName = dctVarnameToEwEName_T_3D_vmean10m[var]
      savepath = os.path.join(greigDir,outputsSubDir.format(modelrun),EwEName)

      if (os.path.isdir(savepath) != True):
          os.mkdir(savepath)
  
      savefn = saveTemplate.format(modelrun,EwEName,iyear,mo)
      savefullpath = os.path.join(savepath, savefn)

      print("Saving file " + savefullpath)
      saveASCFile(savefullpath, DataMon, sigdigfmt)
      
# ////////////////// CONVERT 3D VMEAN TOP 4m /////////////////
# 3D files that are already vert averaged over top 4 layers
for var in lstVars_T_3D_vmean:
  print(var)
  for iyear in range(startYear, endyear):
    for mo in mos:

      inputTemplate = dctFileTemplates_T_3D_vmean4m[var]
      datafullpath = os.path.join(greigDir,NEMOMonthlyNCSubDir,inputTemplate.format(modelrun,iyear,mo))
      dsData = nc.Dataset(datafullpath)
      varsData = dsData.variables[var]
      DataMon = varsData[0,:,:]
      
      #unit conversion
      UnitConvert = dctUnitConvert_T_3D_vmean[var]
      DataMon = DataMon * UnitConvert 
  
      #Save the monthly Data to an .asc file     
      EwEName = dctVarnameToEwEName_T_3D_vmean4m[var]
      savepath = os.path.join(greigDir,outputsSubDir.format(modelrun),EwEName)

      if (os.path.isdir(savepath) != True):
          os.mkdir(savepath)
  
      savefn = saveTemplate.format(modelrun,EwEName,iyear,mo)
      savefullpath = os.path.join(savepath, savefn)

      print("Saving file " + savefullpath)
      saveASCFile(savefullpath, DataMon, sigdigfmt)
      

# 2D vars that don't need vert avg
mo_i = 0
for var in lstVars_T_2D:
  print(var)
  for iyear in range(startYear, endyear):
    for mo in mos:
      mo_i = mo_i + 1
      
      inputTemplate = dctFileTemplates_T_2D[var]
      datafullpath = os.path.join(greigDir,NEMOMonthlyNCSubDir,inputTemplate.format(modelrun,iyear,mo))
      dsData = nc.Dataset(datafullpath)
      varsData = dsData.variables[var]
      DataMon = varsData[0,:,:]
      
      #unit conversion
      UnitConvert = dctUnitConvert_T_2D[var]
      DataMon = DataMon * UnitConvert #GO2021-12-23 - will this work since it's an array?
  
      #Save the monthly Data to an .asc file     
      EwEName = dctVarnameToEwEName_T_2D[var]
      savepath = os.path.join(greigDir,outputsSubDir.format(modelrun),EwEName)

      if (os.path.isdir(savepath) != True):
          os.mkdir(savepath)
  
      savefn = saveTemplate.format(modelrun,EwEName,iyear,mo)
      savefullpath = os.path.join(savepath, savefn)

      print("Saving file " + savefullpath)
      saveASCFile(savefullpath, DataMon, sigdigfmt)

# ////////////////// CONVERT 3D SINGLE LEV ///////////////////
# 3D fields, individual levels (unfinished code)
#mo_i = 0
#for var in lstVars_T_2D:
#  print(var)
#  for iyear in range(startYear, endyear):
#    for mo in mos:
#      mo_i = mo_i + 1
#      
#      inputTemplate = dctFileTemplates_T_2D[var]
#      datafullpath = os.path.join(greigDir,NEMOMonthlyNCSubDir,inputTemplate.format(modelrun,iyear,mo))
#      dsData = nc.Dataset(datafullpath)
#      varsData = dsData.variables[var]
#      DataMon = varsData[0,:,:]
#      
#      #unit conversion
#      UnitConvert = dctUnitConvert_T_2D[var]
#      DataMon = DataMon * UnitConvert #GO2021-12-23 - will this work since it's an array?
#  
#      #Save the monthly Data to an .asc file     
#      EwEName = dctVarnameToEwEName_T_2D[var]
#      savepath = os.path.join(greigDir,outputsSubDir.format(modelrun),EwEName)
#
#      if (os.path.isdir(savepath) != True):
#          os.mkdir(savepath)
#  
#      savefn = saveTemplate.format(modelrun,EwEName,iyear,mo)
#      savefullpath = os.path.join(savepath, savefn)
#
#      print("Saving file " + savefullpath)
#      saveASCFile(savefullpath, DataMon, sigdigfmt)


#for var in lstVars_T_3D_vmean:
#  for iyear in range(startYear, endyear):
#    
#    #data for year
#    inputTemplate = dctFileTemplates_T_3D_vmean10m[var]
#    
#    datafullpath = os.path.join(greigDir,NEMOMonthlyNCSubDir,inputTemplate.format(iyear))
#    dsData = nc.Dataset(datafullpath)
#    varsData = dsData.variables[var]
#    
#    for imon in range(1,13):
#        
#      #Data for month
#      DataMon = varsData[imon-1,:,:]
#  
#      #Get the unit conversion multiplier
#      UnitConvert = dctUnitConvert[var]
#      DataMon = DataMon * UnitConvert #GO2021-12-23 - will this work since it's an array?
#  
#      #Save the monthly Data to an .asc file
#      
#      EwEName = dctUnitConvert_T_3D_vmean[var]
#      #EwEName = EwEName.format(salinitydepthlayer) #GO ??
#      savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
#      if (os.path.isdir(savepath) != True):
#          os.mkdir(savepath)
#  
#      savefn = saveTemplate.format(modelrun,EwEName,iyear,buildSortableString(imon,2))
#      savefullpath = os.path.join(savepath, savefn)
#  
#      print("Saving file " + savefullpath)
#      saveASCFile(savefullpath, DataMon, sigdigfmt)
