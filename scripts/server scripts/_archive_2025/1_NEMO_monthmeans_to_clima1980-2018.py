import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import os
from calendar import monthrange
import math
from GO_helpers import is_leap_year, buildSortableString, saveASCFile, getDataFrame 

#import matplotlib.pyplot as plt
# module load StdEnv/2020
# module load scipy-stack/2020a

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO, last edited by GO 2024-05-23

# purpose: prep files for ECOSPACE - use clima as initial conditions
#          input: monthly light in annual files from RDRS climate model, atten coef (K) from mixing depth Z and salin at chosen depth
#          output: depth integrated and single levels for key vars output as ASC (daily) for daily NEMO model

#Data in:
# 1/ NEMO outputs in 2D or 3D (var mldkz5 / turbocline depth) as monthly means in annual NC files,
# 2/ Salinity - can be either vertical averages over chosen levels or from a specific depth level, 
# 3/ Light as monthly means interpolated to NEMO grid as annual NC files

# Log:
# - GO-2021-12-22 edits to draw from monthly mean files rather than looping through daily NC data 
# - GO2021-12-23 differentiated 1a from 1b as vertmean of salinity (top x layers) and salinity at x m, respectively
#                Added write of K values as ASC for eval
# - GO 2023-03-23 edits to code using PAR2b script (need to re-run PAR3)
# - GO 2023-04-05 moved functions to GO_helpers.py
# - GO 2023-04-06 found that PAR predicted seemed not to match 2015, 2016 in Pawlowicz
# - GO 2023-05-10 added option to switch to over-depth average rather than just selecting a single dpeth
# - GO 2023-05-12 adding calc for u(I) (growth response curve) based on Merico et al (2004), Steele
# - GO 2024-05-17 modifications to read and write to daily instead of monthly means
# - GO 2024-05-21 edits to write out daily files for non-depth integrated salinity (e.g., just 4 m)
# - GO 2024-05-23 mod for 3-day means
# - GO 2024-05-23 fork to this for annual clima - the climatol of advection is for 3day mean rates!
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# //////////////////////////////////// FUNCTIONS //////////////////////////////////////
# now in GO_helpers.py -2023-04-05


# ################# Basic Params #################
startYear = 1980
endyear = 2018
startMonth = 1
endMonth = 12
NEMO_run = "216" #NEMO run code
out_code = "var" #PAR prep code

# not yet used
process_salt = False
process_temp = False
process_mixing = False
process_light = False
process_advec = True


# set depth and whether it is average over depths (outstanding issue with lev width calc)
salinity_depthlev = 10
temperature_depthlev = 10
advec_depthlev = 10

avg_salt = True
if avg_salt == True: 
  out_code_salt = out_code + "{}_{}mAvg".format("salt", salinity_depthlev)
else:
  out_code_salt = out_code + "{}_{}m".format("salt", salinity_depthlev)

avg_temp = True
if avg_temp == True: 
  out_code_temp = out_code + "{}_{}mAvg".format("temp", temperature_depthlev)
else:
  out_code_temp = out_code + "{}_{}m".format("temp", temperature_depthlev)
  
out_code_mixing = out_code + "mixingdep"

avg_advec = True
if avg_advec == True: 
  out_code_advecX = out_code + "{}_{}mAvg".format("advecX", advec_depthlev)
  out_code_advecY = out_code + "{}_{}mAvg".format("advecY", advec_depthlev)
else:
  out_code_advecX = out_code + "{}_{}m".format("advecX", advec_depthlev)
  out_code_advecY = out_code + "{}_{}m".format("advecY", advec_depthlev)


# ################# ASC Params #################
ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"


# ################# PATHS #################
greigDir = "/project/6006412/goldford/ECOSPACE/"
michaelDir = "/project/6006412/mdunphy/nemo_results/SalishSea1500/"
lightSubDir = "DATA/light_monthly/"
plumeSubDir = "DATA/"
depthSubDir = "DATA/"

# pre processed by G
mixingSubDir = "DATA/SS1500-RUN{}/NEMO_annual_NC".format(NEMO_run)
salinitySubDir = "DATA/SS1500-RUN{}/NEMO_annual_NC".format(NEMO_run)
temperatureSubDir = "DATA/SS1500-RUN{}/NEMO_annual_NC".format(NEMO_run)
advecSubDir = "DATA/SS1500-RUN{}/NEMO_annual_NC".format(NEMO_run)

outputsSubDir ="DATA/SS1500-RUN{}/ECOSPACE_in_climayr_{}-{}_{}".format(NEMO_run, startYear, endyear, out_code)

# (for daily)
#mixingSubDir = "SalishSea1500-RUN{}/CDF".format(NEMO_run)
#salinitySubDir = "SalishSea1500-RUN{}/CDF".format(NEMO_run)
#outputsSubDir ="DATA/SS1500-RUN{}/ECOSPACE_in_annualclima_1980-2018{}".format(NEMO_run, PAR_code)

if (os.path.isdir(os.path.join(greigDir,outputsSubDir)) != True):
  os.mkdir(os.path.join(greigDir,outputsSubDir))

# preprocessed monthlies in yearly files
mixingFileName = 'SalishSea1500-RUN{}_MonthlyMean_grid_T_2D_y{}.nc'
lightFileName = "RDRS21_NEMOgrid_light_monthly_{}.nc"  
salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_y{}.nc'
temperatureFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_y{}.nc'
advectXFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_U_y{}.nc' 
advectYFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_V_y{}.nc' 
#salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_singledepth10m_{}.nc' #depth=9.5m
#salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_avgtop10m_y{}.nc'

# (for daily)
#mixingFileName = 'SalishSea1500-RUN{}_1d_grid_T_2D_y{}m{}.nc'
#lightFileName = "RDRS21_NEMOgrid_light_daily_{}.nc"
#salinityFileName ='SalishSea1500-RUN{}_1d_grid_T_y{}m{}.nc' 


# ################# VARS #################
varNameLight = "solar_rad"
varNameSalinity = 'vosaline'
varNameTemp = 'votemper'
varNamesMixing = 'mldkz5' #turbocline
varNameXadvect = 'vozocrtx' #(U)
varNameYadvect = 'vomecrty' #(V) grid
#multiplier is 100 to convert from m/s (NEMO) to cm/s (EwE) - GO 2022-02
#then for tricking Ecospace to run daily divide by 365 * 3 for 3-day
# 'Advection vectors Xvel(,) are in cm/sec convert to km/year, same units as the mvel()/(Dispersal Rate)
# '[km/year] / [cell length]
# AdScale = 315.36 / Me.EcoSpaceData.CellLength
advec_mult = 100 / 365 * 3 # see AdScale in cEcospace


saveTemplate = "{}_{}-{}.asc" # new ewename, startYear, endYear


#Number of depth strata to use
#Depth indexes are zero based
#So 10 will be the 9 index
#nDepthIndexes = 10 #= 9.001736   10.003407
#nDepthIndexes = 23 #= 31.101034   39.11886
#nDepthIndexes = 26 #= 86.96747   109.73707 


# ################# MAP CLIP #################
#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39


# ################# PLUME MASK #################
#Create the plume mask
PlumeRegionFilename = os.path.join(greigDir, plumeSubDir, "Fraser_Plume_Region.asc")
dfPlumeRegion = getDataFrame(PlumeRegionFilename,"-9999.00000000")
dfPlumeMask = dfPlumeRegion == 1 #Plume is region one


# ################# LAND MASK #################
ecospacegrid_f = os.path.join(greigDir, depthSubDir, "ecospacedepthgrid.asc")
df_dep = getDataFrame(ecospacegrid_f,"-9999.00000000")
dfLandMask = df_dep == 0 #mask is where elev =0


# ################# MAIN LOOP #################

light_allyrs = []
mixing_allyrs = []
salt_allyrs = []
temp_allyrs = []
advecx_allyrs = []
advecy_allyrs = []

for iyear in range(startYear, endyear+1):
    print(iyear)
    #YearSubDirTemplate = str(iyear) # unused -GO

    #LIGHT data (hrly data in annual, or monthly in annual)
    lightfullpath = os.path.join(greigDir,lightSubDir,lightFileName.format(iyear))
    dsLight = nc.Dataset(lightfullpath)
    varsLight = dsLight.variables[varNameLight]
    LightYear1 = varsLight[:,:,:]
    LightYear = np.ma.average(LightYear1, axis = 0)
    LightYear = np.expand_dims(LightYear, axis=0)
    
    #SALINITY data (3D)
    salinityfullpath = os.path.join(greigDir,salinitySubDir,salinityFileName.format(NEMO_run,iyear)) # old (preprocessed)
    dsSalinity = nc.Dataset(salinityfullpath)
    varsSalinity = dsSalinity.variables[varNameSalinity]
    if avg_salt == True:
        SaltYear1 = varsSalinity[:,0:salinity_depthlev,:,:]
        SaltYear2 = np.ma.average(SaltYear1, axis = 1) # avg over depths - will reduce to (2,299,132)
        SaltYear = np.ma.average(SaltYear2, axis = 0) # avg over time - reduce to (299,132)
        #SaltYear = np.expand_dims(SaltYear, axis=0)
    else:
        SaltYear1 = varsSalinity[:,salinity_depthlev,:,:]
        SaltYear = np.ma.average(SaltYear1, axis = 0)
    SaltYear = np.expand_dims(SaltYear, axis=0)
    ar_inf = np.where(np.isinf(SaltYear)) # unused? 

    #TEMPERATURE data (3D)
    temperaturefullpath = os.path.join(greigDir,temperatureSubDir,temperatureFileName.format(NEMO_run,iyear)) # old (preprocessed)
    dsTemperature = nc.Dataset(temperaturefullpath)
    varsTemperature = dsTemperature.variables[varNameTemp]
    if avg_temp == True:
        TempYear1 = varsTemperature[:,0:temperature_depthlev,:,:]
        TempYear2 = np.ma.average(TempYear1, axis = 1) # avg over depths - will reduce to (2,299,132)
        TempYear = np.ma.average(TempYear2, axis = 0) # avg over time - reduce to (299,132)     
    else:
        TempYear1 = varsTemperature[:,temperature_depthlev,:,:]
        TempYear = np.ma.average(TempYear1, axis = 0)
    TempYear = np.expand_dims(TempYear, axis=0)
    ar_inf = np.where(np.isinf(TempYear)) # unused? 
                
    #MIXING data 
    mixingFullPath = os.path.join(greigDir, mixingSubDir, mixingFileName.format(NEMO_run, iyear))
    dsMixing = nc.Dataset(mixingFullPath)
    varsMixing = dsMixing.variables['mldkz5']   
    MixingYear1 = varsMixing[:,:,:]
    MixingYear = np.ma.average(MixingYear1, axis = 0) 
    MixingYear = np.expand_dims(MixingYear, axis=0)    

    # ADVEC data
    advecXFullPath = os.path.join(greigDir, advecSubDir, advectXFileName.format(NEMO_run, iyear))
    advecYFullPath = os.path.join(greigDir, advecSubDir, advectYFileName.format(NEMO_run, iyear))
    dsAdvecX = nc.Dataset(advecXFullPath)
    dsAdvecY = nc.Dataset(advecYFullPath)
    varsAdvecX = dsAdvecX.variables[varNameXadvect]
    varsAdvecY = dsAdvecY.variables[varNameYadvect]
    if avg_advec == True:
        AdvecXYear1 = varsAdvecX[:,0:advec_depthlev,:,:]
        AdvecYYear1 = varsAdvecY[:,0:advec_depthlev,:,:]
        AdvecXYear2 = np.ma.average(AdvecXYear1, axis = 1) # avg over depths - will reduce to (2,299,132)
        AdvecYYear2 = np.ma.average(AdvecYYear1, axis = 1)
        AdvecXYear = np.ma.average(AdvecXYear2, axis = 0) # avg over time - reduce to (299,132)
        AdvecYYear = np.ma.average(AdvecYYear2, axis = 0)
        AdvecXYear = np.expand_dims(AdvecXYear, axis=0)
        AdvecYYear = np.expand_dims(AdvecYYear, axis=0)        
    else:
        AdvecXYear1 = varsAdvecX[:,advec_depthlev,:,:]
        AdvecXYear = np.ma.average(AdvecXYear1, axis = 0)
        AdvecYYear1 = varsAdvecY[:,advec_depthlev,:,:]
        AdvecYYear = np.ma.average(AdvecYYear1, axis = 0)
   

    if iyear == startYear:
        light_allyrs = LightYear
        mixing_allyrs = MixingYear
        salt_allyrs = SaltYear
        temp_allyrs = TempYear
        advecx_allyrs = AdvecXYear
        advecy_allyrs = AdvecYYear
    else:
        light_allyrs = np.concatenate((light_allyrs, LightYear), axis=0)
        mixing_allyrs = np.concatenate((mixing_allyrs, MixingYear), axis=0)
        salt_allyrs = np.concatenate((salt_allyrs, SaltYear), axis=0)
        temp_allyrs = np.concatenate((temp_allyrs, TempYear), axis=0)
        advecx_allyrs = np.concatenate((advecx_allyrs, AdvecXYear), axis=0)
        advecy_allyrs = np.concatenate((advecy_allyrs, AdvecYYear), axis=0)           


print(light_allyrs.shape, salt_allyrs.shape, mixing_allyrs.shape, temp_allyrs.shape, advecx_allyrs.shape, advecy_allyrs.shape) 

light_allyrs = np.ma.average(light_allyrs, axis = 0)
salt_allyrs = np.ma.average(salt_allyrs, axis = 0)
mixing_allyrs = np.ma.average(mixing_allyrs, axis = 0)
temp_allyrs = np.ma.average(temp_allyrs, axis = 0)
advecx_allyrs = np.ma.average(advecx_allyrs, axis = 0)
advecy_allyrs = np.ma.average(advecy_allyrs, axis = 0) 

print(light_allyrs.shape, salt_allyrs.shape, mixing_allyrs.shape, temp_allyrs.shape, advecx_allyrs.shape, advecy_allyrs.shape) 

# salinity
sigdigfmt = '%0.1f'
savepath = os.path.join(greigDir,outputsSubDir,out_code_salt) 
if (os.path.isdir(savepath) != True):
    os.mkdir(savepath)
savefn = saveTemplate.format(out_code_salt, startYear, endyear)
savefullpath = os.path.join(savepath, savefn)
print("Saving file " + savefullpath)
saveASCFile(savefullpath, salt_allyrs, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)


# temperature
sigdigfmt = '%0.1f'
savepath = os.path.join(greigDir, outputsSubDir, out_code_temp) 
if (os.path.isdir(savepath) != True):
    os.mkdir(savepath)
savefn = saveTemplate.format(out_code_temp, startYear, endyear)
savefullpath = os.path.join(savepath, savefn)
print("Saving file " + savefullpath)
saveASCFile(savefullpath, temp_allyrs, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)

# mixing
sigdigfmt = '%0.1f'
savepath = os.path.join(greigDir, outputsSubDir, out_code_mixing) 
if (os.path.isdir(savepath) != True):
    os.mkdir(savepath)
savefn = saveTemplate.format(out_code_mixing, startYear, endyear)
savefullpath = os.path.join(savepath, savefn)
print("Saving file " + savefullpath)
saveASCFile(savefullpath, mixing_allyrs, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)

# X advec
sigdigfmt = '%0.4f'
savepath = os.path.join(greigDir, outputsSubDir, out_code_advecX) 
if (os.path.isdir(savepath) != True):
    os.mkdir(savepath)
savefn = saveTemplate.format(out_code_advecX, startYear, endyear)
savefullpath = os.path.join(savepath, savefn)
print("Saving advec file " + savefullpath)
saveASCFile(savefullpath, advecx_allyrs, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)

# Y advec
sigdigfmt = '%0.4f'
savepath = os.path.join(greigDir, outputsSubDir, out_code_advecY) 
if (os.path.isdir(savepath) != True):
    os.mkdir(savepath)
savefn = saveTemplate.format(out_code_advecY, startYear, endyear)
savefullpath = os.path.join(savepath, savefn)
print("Saving advec file " + savefullpath)
saveASCFile(savefullpath, advecy_allyrs, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)

   
