import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import os
from calendar import monthrange
import math
from GO_helpers import is_leap_year, buildSortableString, saveASCFile, getDataFrame 

#import matplotlib.pyplot as plt

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO, last edited by GO 2024-05-23

# purpose: prep files for ECOSPACE
#          input: daily mean NEMO outputs
#          output: depth integrated and single levels for key vars output as ASC (daily) for daily NEMO model

#Data in:
# 1/ NEMO outputs in 2D Layer data (var mldkz5 / turbocline depth) as daily means in annual NC files, DIRECT from NEMO outputs "xxx_1d_grid_T_2D_1980.nc"
# 2/ left it hea
# 2/ Salinity - can be either vertical averages over chosen levels or from a specific depth level, DIRECT from NEMO outputs "xxx_1d_grid_T_selectedlevs_1980.nc"
# 3/ Light as daily means from climate re-analysis (e.g. ERA5) with near-surface downward radiative flux (var msdwswrf Wm-2) interpolated to NEMO grid as annual NC files# 
# file name e.g. "ERA5_NEMOgrid_light_daily_1980.nc".

# Log:
# - GO 2024-05-21 edits to write out daily files for non-depth integrated salinity (e.g., just 4 m)
# - GO 2024-05-22 note that the depth averaging is not right - needs to be weighted by level width (code exists somewhere to fix), not major issue for shallows
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# //////////////////////////////////// FUNCTIONS //////////////////////////////////////
# now in GO_helpers.py -2023-04-05


# originally CDO was used, it is much faster but I don't have time to fiddle
# ################# Basic Params #################
startYear = 2015 
endyear = 2015
startMonth = 1
endMonth = 2
NEMO_run = "216" #NEMO run code
out_code = "var" #PAR prep code

process_salt = True
process_temp = False
process_mixing = False
process_light = False
process_advec = False


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
#lightSubDir = "DATA/light_monthly/"
lightSubDir = "DATA/light_daily/"
plumeSubDir = "DATA/"
depthSubDir = "DATA/"
# old
#mixingSubDir = "DATA/SS1500-RUN{}/NEMO_annual_NC".format(NEMO_run)
#salinitySubDir = "DATA/SS1500-RUN{}/NEMO_annual_NC".format(NEMO_run)
#outputsSubDir ="DATA/SS1500-RUN{}/ECOSPACE_in_{}".format(NEMO_run, PAR_code)
# new (for daily)
mixingSubDir = "SalishSea1500-RUN{}/CDF".format(NEMO_run)
salinitySubDir = "SalishSea1500-RUN{}/CDF".format(NEMO_run)
temperatureSubDir = "SalishSea1500-RUN{}/CDF".format(NEMO_run)
advecSubDir = "SalishSea1500-RUN{}/CDF".format(NEMO_run)
outputsSubDir ="DATA/SS1500-RUN{}/ECOSPACE_in_daily_vars".format(NEMO_run)

if (os.path.isdir(os.path.join(greigDir,outputsSubDir)) != True):
  os.mkdir(os.path.join(greigDir,outputsSubDir))

# OLD
#mixingFileName = 'SalishSea1500-RUN{}_MonthlyMean_grid_T_2D_y{}.nc'
#lightFileName = "RDRS21_NEMOgrid_light_monthly_{}.nc"  
#salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_y{}.nc' 
#salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_singledepth10m_{}.nc' #depth=9.5m
#salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_avgtop10m_y{}.nc'

# NEW - for daily
mixingFileName = 'SalishSea1500-RUN{}_1d_grid_T_2D_y{}m{}.nc'
lightFileName = "RDRS21_NEMOgrid_light_daily_{}.nc"
salinityFileName ='SalishSea1500-RUN{}_1d_grid_T_y{}m{}.nc'
temperatureFileName ='SalishSea1500-RUN{}_1d_grid_T_y{}m{}.nc' 
advectXFileName ='SalishSea1500-RUN{}_1d_grid_U_y{}m{}.nc' 
advectYFileName ='SalishSea1500-RUN{}_1d_grid_V_y{}m{}.nc' 


# ///////////////////// VARS ////////////////////////
varNameLight = "solar_rad"
varNameSalinity = 'vosaline'
varNameTemp = 'votemper'
varNamesMixing = 'mldkz5' #turbocline
varNameXadvect = 'vozocrtx' #(U)
varNameYadvect = 'vomecrty' #(V) grid
#multiplier is 100 to convert from m/s (NEMO) to cm/s (EwE) - GO 2022-02
#then for tricking Ecospace to run daily divide by 365
# 'Advection vectors Xvel(,) are in cm/sec convert to km/year, same units as the mvel()/(Dispersal Rate)
# '[km/year] / [cell length]
# AdScale = 315.36 / Me.EcoSpaceData.CellLength
advec_mult = 100 / 365 # see AdScale in cEcospace
#lstVars = ['mldkz5']; 

#saveTemplate = "{}_{}_{}.asc" # old
saveTemplate = "{}_{}_{}_{}.asc" # new: ewename, year, mo, day

#Number of depth strata to use
#Depth indexes are zero based
#So 10 will be the 9 index
#nDepthIndexes = 10 #= 9.001736   10.003407
#nDepthIndexes = 23 #= 31.101034   39.11886
#nDepthIndexes = 26 #= 86.96747   109.73707 

# ///////////////// MAP CLIP ///////////////////////
#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39

# ///////////////// PLUME MASK /////////////////////
#Create the plume mask
PlumeRegionFilename = os.path.join(greigDir, plumeSubDir, "Fraser_Plume_Region.asc")
dfPlumeRegion = getDataFrame(PlumeRegionFilename,"-9999.00000000")
dfPlumeMask = dfPlumeRegion == 1 #Plume is region one

# ///////////////// LAND MASK /////////////////////
ecospacegrid_f = os.path.join(greigDir, depthSubDir, "ecospacedepthgrid.asc")
df_dep = getDataFrame(ecospacegrid_f,"-9999.00000000")
dfLandMask = df_dep == 0 #mask is where elev =0

# ///////////////// MAIN LOOP /////////////////////

for iyear in range(startYear, endyear+1):

    #YearSubDirTemplate = str(iyear) # unused -GO

    #LIGHT data (hrly data in annual, or monthly in annual)
    lightfullpath = os.path.join(greigDir,lightSubDir,lightFileName.format(iyear))
    dsLight = nc.Dataset(lightfullpath)
    varsLight = dsLight.variables[varNameLight]


    for imon in range(startMonth,13):
        
        #SALINITY data
        if process_salt:
            salinityfullpath = os.path.join(michaelDir, salinitySubDir, salinityFileName.format(NEMO_run, iyear, buildSortableString(imon,2))) 
            dsSalinity = nc.Dataset(salinityfullpath)
            varsSalinity = dsSalinity.variables[varNameSalinity]
        
        #TEMPERATURE data
        if process_temp:
            temperaturefullpath = os.path.join(michaelDir, temperatureSubDir, temperatureFileName.format(NEMO_run, iyear, buildSortableString(imon,2))) 
            dsTemperature = nc.Dataset(temperaturefullpath)
            varsTemperature = dsTemperature.variables[varNameTemp]
        
        #MIXING data 
        if process_mixing:
            mixingFullPath = os.path.join(michaelDir, mixingSubDir, mixingFileName.format(NEMO_run, iyear, buildSortableString(imon,2)))
            dsMixing = nc.Dataset(mixingFullPath)
            varsMixing = dsMixing.variables['mldkz5']
        
        #ADVEC data
        if process_advec:
            advecXFullPath = os.path.join(michaelDir, advecSubDir, advectXFileName.format(NEMO_run, iyear, buildSortableString(imon,2)))
            advecYFullPath = os.path.join(michaelDir, advecSubDir, advectYFileName.format(NEMO_run, iyear, buildSortableString(imon,2)))
            dsAdvecX = nc.Dataset(advecXFullPath)
            dsAdvecY = nc.Dataset(advecYFullPath)
            varsAdvecX = dsAdvecX.variables[varNameXadvect]
            varsAdvecY = dsAdvecY.variables[varNameYadvect]        
        
        # is leap year?
        if imon == 2 and is_leap_year(iyear):
            num_days = 29
        else:
            num_days = 28 if imon == 2 else 30 if imon in [4, 6, 9, 11] else 31
        
        for iday in range(1,num_days+1):
        
        
            # #############  LIGHT  #############
            # not done
            LightDay = varsLight[iday-1,:,:]          

 
            # #############  SALINITY  #############
            if process_salt:
                sigdigfmt = '%0.1f'
                if avg_salt == True:
                  SalinityDay1 = varsSalinity[iday-1,0:salinity_depthlev,:,:]
                  #print(SalinityDay1.shape) # will be (10,299,132)
                  SalinityDay = np.ma.average(SalinityDay1, axis = 0)
                else:
                  SalinityDay = varsSalinity[iday-1,salinity_depthlev,:,:]
                ar_inf = np.where(np.isinf(SalinityDay)) # unused?
                
                savepath = os.path.join(greigDir,outputsSubDir,out_code_salt) 
                if (os.path.isdir(savepath) != True):
                    os.mkdir(savepath)
                savefn = saveTemplate.format(out_code_salt,iyear,buildSortableString(imon,2),buildSortableString(iday,2))
                savefullpath = os.path.join(savepath, savefn)
                print("Saving file " + savefullpath)
                saveASCFile(savefullpath, SalinityDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
        
        
            # #############  TEMPERATURE  #############
            if process_temp:
                sigdigfmt = '%0.1f'
                if avg_temp == True:
                  TemperDay1 = varsTemperature[iday-1,0:temperature_depthlev,:,:]
                  TemperDay = np.ma.average(TemperDay1, axis = 0)
                else:
                  TemperDay = varsTemperature[iday-1,temperature_depthlev,:,:]
                ar_inf = np.where(np.isinf(TemperDay)) # unused?
                
                savepath = os.path.join(greigDir,outputsSubDir,out_code_temp) 
                if (os.path.isdir(savepath) != True):
                    os.mkdir(savepath)
                savefn = saveTemplate.format(out_code_temp,iyear,buildSortableString(imon,2),buildSortableString(iday,2))
                savefullpath = os.path.join(savepath, savefn)
                print("Saving file " + savefullpath)
                saveASCFile(savefullpath, TemperDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
        
        
            # #############  MIXING  #############
            if process_mixing:
                sigdigfmt = '%0.1f'
                MixingDay = varsMixing[iday-1,:,:] 
                savepath = os.path.join(greigDir,outputsSubDir,out_code_mixing) 
                if (os.path.isdir(savepath) != True):
                    os.mkdir(savepath)
                savefn = saveTemplate.format(out_code_mixing,iyear,buildSortableString(imon,2),buildSortableString(iday,2))
                savefullpath = os.path.join(savepath, savefn)
                print("Saving file " + savefullpath)
                saveASCFile(savefullpath, MixingDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
        
        
            # #############  ADVECTION  #############
            if process_advec:
                sigdigfmt = '%0.4f'
                if avg_advec == True:
                  AdvecXDay1 = varsAdvecX[iday-1,0:advec_depthlev,:,:]
                  AdvecYDay1 = varsAdvecY[iday-1,0:advec_depthlev,:,:]
                  AdvecXDay = np.ma.average(AdvecXDay1, axis = 0)
                  AdvecYDay = np.ma.average(AdvecYDay1, axis = 0)
                else:
                  AdvecXDay = varsAdvecX[iday-1,advec_depthlev,:,:]
                  AdvecYDay = varsAdvecY[iday-1,advec_depthlev,:,:]
                AdvecXDay = AdvecXDay * advec_mult # convert from m/s to cm/s, and s/365 for daily hack
                AdvecYDay = AdvecYDay * advec_mult
                ar_inf = np.where(np.isinf(AdvecXDay)) # unused?

                # X advec
                savepath = os.path.join(greigDir,outputsSubDir,out_code_advecX) 
                if (os.path.isdir(savepath) != True):
                    os.mkdir(savepath)
                savefn = saveTemplate.format(out_code_advecX,iyear,buildSortableString(imon,2),buildSortableString(iday,2))
                savefullpath = os.path.join(savepath, savefn)
                print("Saving advec file " + savefullpath)
                saveASCFile(savefullpath, AdvecXDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
                
                # Y advec
                savepath = os.path.join(greigDir,outputsSubDir,out_code_advecY) 
                if (os.path.isdir(savepath) != True):
                    os.mkdir(savepath)
                savefn = saveTemplate.format(out_code_advecY,iyear,buildSortableString(imon,2),buildSortableString(iday,2))
                savefullpath = os.path.join(savepath, savefn)
                print("Saving advec file " + savefullpath)
                saveASCFile(savefullpath, AdvecYDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
        
    #               #export Sal for debug 2023-03-23
    #               sigdigfmt = '%0.0f'
    #               EwEName = "just_salin"
    #               #EwEName = "RUN102b-Kfromsal_singlelev10m"
    #               #EwEName = "RUN102a-Kfromsal_vertmeantop10m" #GO alternative
    #               #EwEName = EwEName.format(salinitydepthlayer) #GO ??
    #               savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
    #               if (os.path.isdir(savepath) != True):
    #                   os.mkdir(savepath)
    #               savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
    #               savefullpath = os.path.join(savepath, savefn)
    #               print("Saving file " + savefullpath)
    #               saveASCFile(savefullpath, SalinityMon,bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
                
            
        if (imon == endMonth) and (iyear == endyear):
            break

        # why?
        if imon == 12:
            break
