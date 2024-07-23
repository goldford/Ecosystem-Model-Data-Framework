import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import os
from calendar import monthrange
import math
from GO_helpers import is_leap_year, buildSortableString, saveASCFile, getDataFrame 

# module load StdEnv/2020
# module load scipy-stack/2020a

#import matplotlib.pyplot as plt

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO, last edited by GO 2024-05-23

# purpose: prep files for ECOSPACE - 3 day averages 
#          input: daily mean NEMO outputs
#          output: depth integrated and single levels for key vars output as ASC (daily) for daily NEMO model

#Data in:
# 1/ NEMO outputs in 2D or 3D (var mldkz5 / turbocline depth) as daily means in annual NC files, DIRECT from NEMO outputs "xxx_1d_grid_T_2D_1980.nc"
# 2/ Salinity - can be either vertical averages over chosen levels or from a specific depth level, DIRECT from NEMO outputs "xxx_1d_grid_T_selectedlevs_1980.nc"
# 3/ Light as daily means from climate re-analysis (e.g. ERA5) with near-surface downward radiative flux (var msdwswrf Wm-2) interpolated to NEMO grid as annual NC files# 
# file name e.g. "ERA5_NEMOgrid_light_daily_1980.nc".

# Log:
# - GO 2024-05-21 edits to write out daily files for non-depth integrated salinity (e.g., just 4 m)
# - GO 2024-05-22 note that the depth averaging is not right - needs to be weighted by level width (code exists somewhere to fix), not major issue for shallows
# - GO 2024-05-23 Ecospace and Ecosim have max of 500 yr simulations - must therefore change for 2015-2018 validation using PP to 3-day model
#                 Leap years divide by 3 well - 122 'years' per year, others are 121 full time steps and the 0.66 x 3 days get dropped
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# //////////////////////////////////// FUNCTIONS //////////////////////////////////////
# now in GO_helpers.py -2023-04-05


# originally CDO was used, it is much faster but I don't have time to fiddle
# ################# Basic Params #################
startYear = 2015 
endyear = 2018
startMonth = 1
endMonth = 12
NEMO_run = "216" #NEMO run code
out_code = "var" #PAR prep code

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
outputsSubDir ="DATA/SS1500-RUN{}/ECOSPACE_in_3daily_vars".format(NEMO_run)

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
#then for tricking Ecospace to run daily divide by 365 * 3 for 3-day
# 'Advection vectors Xvel(,) are in cm/sec convert to km/year, same units as the mvel()/(Dispersal Rate)
# '[km/year] / [cell length]
# AdScale = 315.36 / Me.EcoSpaceData.CellLength
advec_mult = 100 / 365 * 3 # see AdScale in cEcospace
#lstVars = ['mldkz5']; 

saveTemplate = "{}_{}_{}.asc" # new: ewename, year, dayofyear

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
i = 0 # the hacked 'year' (3 day)
for iyear in range(startYear, endyear+1):

    print(iyear)
    #YearSubDirTemplate = str(iyear) # unused -GO

    #LIGHT data (hrly data in annual, or monthly in annual)
    lightfullpath = os.path.join(greigDir,lightSubDir,lightFileName.format(iyear))
    dsLight = nc.Dataset(lightfullpath)
    varsLight = dsLight.variables[varNameLight]

    varsSalinity_12mo = []
    varsTemperature_12mo = []
    varsMixing_12mo = []
    varsAdvecX_12mo = []
    varsAdvecY_12mo = []

    # read and create stack of all months of data
    for imon in range(startMonth,13):
        print(imon)
        
        #SALINITY data
        if process_salt:
            salinityfullpath = os.path.join(michaelDir, salinitySubDir, salinityFileName.format(NEMO_run, iyear, buildSortableString(imon,2))) 
            dsSalinity = nc.Dataset(salinityfullpath)
            varsSalinity = dsSalinity.variables[varNameSalinity]
            if imon == startMonth:
                varsSalinity_12mo = varsSalinity
            else:
                varsSalinity_12mo = np.concatenate((varsSalinity_12mo, varsSalinity), axis=0)
        
        #TEMPERATURE data
        if process_temp:
            temperaturefullpath = os.path.join(michaelDir, temperatureSubDir, temperatureFileName.format(NEMO_run, iyear, buildSortableString(imon,2))) 
            dsTemperature = nc.Dataset(temperaturefullpath)
            varsTemperature = dsTemperature.variables[varNameTemp]
            if imon == startMonth:
                varsTemperature_12mo = varsTemperature
            else:
                varsTemperature_12mo = np.concatenate((varsTemperature_12mo, varsTemperature), axis=0)
        
        
        #MIXING data 
        if process_mixing:
            mixingFullPath = os.path.join(michaelDir, mixingSubDir, mixingFileName.format(NEMO_run, iyear, buildSortableString(imon,2)))
            dsMixing = nc.Dataset(mixingFullPath)
            varsMixing = dsMixing.variables['mldkz5']
            if imon == startMonth:
                varsMixing_12mo = varsMixing
            else:
                varsMixing_12mo = np.concatenate((varsMixing_12mo, varsMixing), axis=0)
        
        #ADVEC data
        if process_advec:
            advecXFullPath = os.path.join(michaelDir, advecSubDir, advectXFileName.format(NEMO_run, iyear, buildSortableString(imon,2)))
            advecYFullPath = os.path.join(michaelDir, advecSubDir, advectYFileName.format(NEMO_run, iyear, buildSortableString(imon,2)))
            dsAdvecX = nc.Dataset(advecXFullPath)
            dsAdvecY = nc.Dataset(advecYFullPath)
            varsAdvecX = dsAdvecX.variables[varNameXadvect]
            varsAdvecY = dsAdvecY.variables[varNameYadvect]
            if imon == startMonth:
                varsAdvecX_12mo = varsAdvecX
                varsAdvecY_12mo = varsAdvecY
            else:
                varsAdvecX_12mo = np.concatenate((varsAdvecX_12mo, varsAdvecX), axis=0)
                varsAdvecY_12mo = np.concatenate((varsAdvecY_12mo, varsAdvecY), axis=0)

        
        if (imon == endMonth) and (iyear == endyear):
            break

        # why?
        if imon == 12:
            break

    leapTF = is_leap_year(iyear)
    if leapTF:
        num_days = 366
    else:
        num_days = 365
    
    num_days = num_days // 3 # 3 day 'years'
    print("number of 'years' per year")
    print(num_days)
    
    if not leapTF:
        num_days += 1 # this helps catch day 364 and 365
    
    j = 0 # for file naming
    for iday in range(1,num_days+1):

        day_strt = (iday-1)+(j*2)
        day_end = day_strt+2
        middle_day = day_strt+2
        
        # catch if last block of days (364,365) is only 2 long
        if not leapTF:
            if iday == num_days+1:
                day_end = day_strt+1
                middle_day = day_strt+1
            
        # #############  LIGHT  #############
        # not done
        #LightDay = varsLight[day_strt:day_end,:,:]  
        

        # #############  SALINITY  #############
        if process_salt:
            
            if avg_salt == True:
              SalinityDay1 = varsSalinity_12mo[day_strt:day_end,0:salinity_depthlev,:,:]
              #print(SalinityDay1.shape) # #print(SalinityDay1.shape) # will be (2,10,299,132)
              SalinityDay2 = np.ma.average(SalinityDay1, axis = 1) # avg over depths - will reduce to (2,299,132)
              SalinityDay = np.ma.average(SalinityDay2, axis = 0) # avg over time - reduce to (299,132)
            else:          
              SalinityDay1 = varsSalinity_12mo[day_strt:day_end,salinity_depthlev,:,:]
              SalinityDay = np.ma.average(SalinityDay1, axis = 0)
            ar_inf = np.where(np.isinf(SalinityDay)) # unused?
            
            savepath = os.path.join(greigDir,outputsSubDir,out_code_salt) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            savefn = saveTemplate.format(out_code_salt,iyear,buildSortableString(middle_day,2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving file " + savefullpath)
            sigdigfmt = '%0.1f'
            saveASCFile(savefullpath, SalinityDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
    
    
        # #############  TEMPERATURE  #############
        if process_temp:
            sigdigfmt = '%0.1f'
            if avg_temp == True:
              TemperDay1 = varsTemperature_12mo[day_strt:day_end,0:temperature_depthlev,:,:]
              TemperDay2 = np.ma.average(TemperDay1, axis = 1)
              TemperDay = np.ma.average(TemperDay2, axis = 0)
            else:
              TemperDay1 = varsTemperature_12mo[day_strt:day_end,temperature_depthlev,:,:]
              TemperDay = np.ma.average(TemperDay1, axis = 0)
            ar_inf = np.where(np.isinf(TemperDay)) # unused?
            
            savepath = os.path.join(greigDir,outputsSubDir,out_code_temp) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            savefn = saveTemplate.format(out_code_temp,iyear,buildSortableString(middle_day,2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, TemperDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
    
    
        # #############  MIXING  #############
        if process_mixing:
            MixingDay1 = varsMixing_12mo[day_strt:day_end,:,:]
            MixingDay = np.ma.average(MixingDay1, axis = 0)
            savepath = os.path.join(greigDir,outputsSubDir,out_code_mixing) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            savefn = saveTemplate.format(out_code_mixing,iyear,buildSortableString(middle_day,2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving file " + savefullpath)
            sigdigfmt = '%0.1f'
            saveASCFile(savefullpath, MixingDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
    
    
        # #############  ADVECTION  #############
        if process_advec:
            sigdigfmt = '%0.4f'
            if avg_advec == True:
              AdvecXDay1 = varsAdvecX_12mo[day_strt:day_end,0:advec_depthlev,:,:]         
              AdvecYDay1 = varsAdvecY_12mo[day_strt:day_end,0:advec_depthlev,:,:]
              AdvecXDay2 = np.ma.average(AdvecXDay1, axis = 1)
              AdvecYDay2 = np.ma.average(AdvecYDay1, axis = 1)
              AdvecXDay = np.ma.average(AdvecXDay2, axis = 0)
              AdvecYDay = np.ma.average(AdvecYDay2, axis = 0)
            else:
              AdvecXDay1 = varsAdvecX_12mo[day_strt:day_end,advec_depthlev,:,:]
              AdvecYDay1 = varsAdvecY_12mo[day_strt:day_end,advec_depthlev,:,:]
              AdvecXDay = np.ma.average(AdvecXDay1, axis = 0)
              AdvecYDay = np.ma.average(AdvecYDay1, axis = 0)
            AdvecXDay = AdvecXDay * advec_mult # convert from m/s to cm/s, and s/365 for daily hack
            AdvecYDay = AdvecYDay * advec_mult
            ar_inf = np.where(np.isinf(AdvecXDay)) # unused?

            # X advec
            savepath = os.path.join(greigDir,outputsSubDir,out_code_advecX) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            savefn = saveTemplate.format(out_code_advecX,iyear,buildSortableString(middle_day,2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving advec file " + savefullpath)
            saveASCFile(savefullpath, AdvecXDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
            
            # Y advec
            savepath = os.path.join(greigDir,outputsSubDir,out_code_advecY) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            savefn = saveTemplate.format(out_code_advecY,iyear,buildSortableString(middle_day,2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving advec file " + savefullpath)
            saveASCFile(savefullpath, AdvecYDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
    
        j+=1    



    
       
        # # is leap year?
        # if imon == 2 and is_leap_year(iyear):
            # num_days = 29
        # else:
            # num_days = 28 if imon == 2 else 30 if imon in [4, 6, 9, 11] else 31
        
        # num_days = num_days // 3
        # print(num_days)
        
        # j = 0 # 
        # for iday in range(1,num_days+1):
        
            # day_strt = (iday-1)+(j*2)
            # day_end = day_strt+2
        
            # middle_day = day_strt+2 # for file name (middle mnth-day of 3 day block)
        
            # # #############  LIGHT  #############
            # # not done
            # #LightDay = varsLight[day_strt:day_end,:,:]  
            
 
            # # #############  SALINITY  #############
            # if process_salt:
                # sigdigfmt = '%0.1f'
                # if avg_salt == True:
                  # SalinityDay1 = varsSalinity[day_strt:day_end,0:salinity_depthlev,:,:]
                  # #print(SalinityDay1.shape) # #print(SalinityDay1.shape) # will be (2,10,299,132)
                  # SalinityDay2 = np.ma.average(SalinityDay1, axis = 1) # avg over depths - will reduce to (2,299,132)
                  # SalinityDay = np.ma.average(SalinityDay2, axis = 0) # avg over time - reduce to (299,132)
                # else:          
                  # SalinityDay1 = varsSalinity[day_strt:day_end,salinity_depthlev,:,:]
                  # SalinityDay = np.ma.average(SalinityDay1, axis = 0)
                # ar_inf = np.where(np.isinf(SalinityDay)) # unused?
                
                # savepath = os.path.join(greigDir,outputsSubDir,out_code_salt) 
                # if (os.path.isdir(savepath) != True):
                    # os.mkdir(savepath)
                # savefn = saveTemplate.format(out_code_salt,iyear,buildSortableString(imon,2),buildSortableString(middle_day,2))
                # savefullpath = os.path.join(savepath, savefn)
                # print("Saving file " + savefullpath)
                # saveASCFile(savefullpath, SalinityDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
        
        
            # # #############  TEMPERATURE  #############
            # if process_temp:
                # sigdigfmt = '%0.1f'
                # if avg_temp == True:
                  # TemperDay1 = varsTemperature[day_strt:day_end,0:temperature_depthlev,:,:]
                  # TemperDay2 = np.ma.average(TemperDay1, axis = 1)
                  # TemperDay = np.ma.average(TemperDay2, axis = 0)
                # else:
                  # TemperDay1 = varsTemperature[day_strt:day_end,temperature_depthlev,:,:]
                  # TemperDay = np.ma.average(TemperDay1, axis = 0)
                # ar_inf = np.where(np.isinf(TemperDay)) # unused?
                
                # savepath = os.path.join(greigDir,outputsSubDir,out_code_temp) 
                # if (os.path.isdir(savepath) != True):
                    # os.mkdir(savepath)
                # savefn = saveTemplate.format(out_code_temp,iyear,buildSortableString(imon,2),buildSortableString(middle_day,2))
                # savefullpath = os.path.join(savepath, savefn)
                # print("Saving file " + savefullpath)
                # saveASCFile(savefullpath, TemperDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
        
        
            # # #############  MIXING  #############
            # if process_mixing:
                # sigdigfmt = '%0.1f'
                # MixingDay1 = varsMixing[day_strt:day_end,:,:]
                # MixingDay = np.ma.average(MixingDay1, axis = 0)
                # savepath = os.path.join(greigDir,outputsSubDir,out_code_mixing) 
                # if (os.path.isdir(savepath) != True):
                    # os.mkdir(savepath)
                # savefn = saveTemplate.format(out_code_mixing,iyear,buildSortableString(imon,2),buildSortableString(middle_day,2))
                # savefullpath = os.path.join(savepath, savefn)
                # print("Saving file " + savefullpath)
                # saveASCFile(savefullpath, MixingDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
        
        
            # # #############  ADVECTION  #############
            # if process_advec:
                # sigdigfmt = '%0.4f'
                # if avg_advec == True:
                  # AdvecXDay1 = varsAdvecX[day_strt:day_end,0:advec_depthlev,:,:]         
                  # AdvecYDay1 = varsAdvecY[day_strt:day_end,0:advec_depthlev,:,:]
                  # AdvecXDay2 = np.ma.average(AdvecXDay1, axis = 1)
                  # AdvecYDay2 = np.ma.average(AdvecYDay1, axis = 1)
                  # AdvecXDay = np.ma.average(AdvecXDay2, axis = 0)
                  # AdvecYDay = np.ma.average(AdvecYDay2, axis = 0)
                # else:
                  # AdvecXDay1 = varsAdvecX[day_strt:day_end,advec_depthlev,:,:]
                  # AdvecYDay1 = varsAdvecY[day_strt:day_end,advec_depthlev,:,:]
                  # AdvecXDay = np.ma.average(AdvecXDay1, axis = 0)
                  # AdvecYDay = np.ma.average(AdvecYDay1, axis = 0)
                # AdvecXDay = AdvecXDay * advec_mult # convert from m/s to cm/s, and s/365 for daily hack
                # AdvecYDay = AdvecYDay * advec_mult
                # ar_inf = np.where(np.isinf(AdvecXDay)) # unused?

                # # X advec
                # savepath = os.path.join(greigDir,outputsSubDir,out_code_advecX) 
                # if (os.path.isdir(savepath) != True):
                    # os.mkdir(savepath)
                # savefn = saveTemplate.format(out_code_advecX,iyear,buildSortableString(imon,2),buildSortableString(middle_day,2))
                # savefullpath = os.path.join(savepath, savefn)
                # print("Saving advec file " + savefullpath)
                # saveASCFile(savefullpath, AdvecXDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
                
                # # Y advec
                # savepath = os.path.join(greigDir,outputsSubDir,out_code_advecY) 
                # if (os.path.isdir(savepath) != True):
                    # os.mkdir(savepath)
                # savefn = saveTemplate.format(out_code_advecY,iyear,buildSortableString(imon,2),buildSortableString(middle_day,2))
                # savefullpath = os.path.join(savepath, savefn)
                # print("Saving advec file " + savefullpath)
                # saveASCFile(savefullpath, AdvecYDay, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
        
            # j+=1    
            

