import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import os
from calendar import monthrange
import math
from GO_helpers import buildSortableString, saveASCFile, getDataFrame # 2023-04

#import matplotlib.pyplot as plt

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO + JB 2021, last edited by GO 2023-10-11

# purpose: prep files for ECOSPACE
#          input: light from RDRS climate model, atten coef (K) from mixing depth Z and salin
#          output: computes integrated PAR across water depths, u(I) (growth response) following Merico et al (2004), Steele
#          PAR = photosynthetically active radiation
#     Method codes:
#         - PAR0 - fixed mixed depth, fixed K value
#         - PAR1a - variable mixing depth, fixed K value
#         - PAR1b - variable mixed depth, fixed K value
#         - PAR2a - fixed mix depth, linear salinity->turbidity to modify K by cell, function using salinity vertmean 0-10 m
#         - PAR2b - fixed mix depth, linear salinity->turbidity f(n) using salinity at one depth
#     --> - PAR3 -  variable mix depth (turbocline threshold), turbidity from 2b, salin at one depth or avg depth
#         - PAR3b - no longer using 

#Data in:

# 1/ Mixing Layer data (var mldkz5 / turbocline depth) as monthly means in annual NC files, post-processed in other script from NEMO "xxx_MonthlyMeans_grid_T_2D_1980.nc"
# 2/ Salinity - can be either vertical averages over chosen levels or from a specific depth level, post-processed from NEMO in annual NC files "xxx_MonthlyVertMean_grid_T_selectedlevs_1980.nc"
# 3/ Light as monthly means from climate re-analysis (e.g. ERA5) with near-surface downward radiative flux (var msdwswrf Wm-2) interpolated to NEMO grid as annual NC files# 
# file name e.g. "ERA5_NEMOgrid_light_daily_1980.nc".

# Log:
# - GO-2021-12-22 edits to draw from monthly mean files rather than looping through daily NC data 
# - GO2021-12-23 differentiated 1a from 1b as vertmean of salinity (top x layers) and salinity at x m, respectively
#                Added write of K values as ASC for eval
# - GO 2023-03-23 edits to code using PAR2b script (need to re-run PAR3)
# - GO 2023-04-05 moved functions to GO_helpers.py
# - GO 2023-04-06 found that PAR predicted seemed not to match 2015, 2016 in Pawlowicz
# - GO 2023-05-10 added option to switch to over-depth average rather than just selecting a single dpeth
# - GO 2023-05-12 adding calc for u(I) (growth response curve) based on Merico et al (2004), Steele
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# //////////////////////////////////// FUNCTIONS //////////////////////////////////////
# now in GO_helpers.py -2023-04-05


# ////////////////////////////////// Basic Params ///////////////////////////////////////
startYear = 1979 
endyear = startYear + 40
  
NEMO_run = "216" #NEMO run code
PAR_code = "PAR3" #PAR prep code

# declare if single level or pre-calculated mean across depths
avg_salin = True

if avg_salin == True: 
  PAR_code = PAR_code + "_Sal10mAvg"
else:
  PAR_code = PAR_code + "_Sal10m"

# ////////////////////////////////// PAR Ir PARAMS ///////////////////////////////////////
# irradiance PAR across depths
salinitydepthlayer = 4
FixedMixingDepth = 10
b = -0.063 
#a = 0.08529 
a = 1.94


#Variables use for Light attenuation
pPAR = 0.44  #proportion of irrad useful for phyto
alb = 0.067  #albedo reflectance
K = 0.05     #baseline K (clear water)


# not completed
# ////////////////////////////////// Growth Response params ///////////////////////////////////////
# eqn from Merico et al 2004 Appendix A pp. 1822
# Light saturation constants W/m-2
I_sd = 15  #diatoms
I_sf = 15  #flag
I_sdf = 15 #dinoflag
I_sc = 45  #coccolith

# relative max growth rates
#u_max_d = 
#u_max_ 

# ///////////////////////////////////// ASC Params ////////////////////////////////////////
ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"


#////////////////// /////////////////// PATHS /////////////////////////////////////////////
greigDir = "/project/6006412/goldford/ECOSPACE/"
lightSubDir = "DATA/light_monthly/"
plumeSubDir = "DATA/"
depthSubDir = "DATA/"
mixingSubDir = "DATA/SS1500-RUN{}/NEMO_annual_NC".format(NEMO_run)
salinitySubDir = "DATA/SS1500-RUN{}/NEMO_annual_NC".format(NEMO_run)
outputsSubDir ="DATA/SS1500-RUN{}/ECOSPACE_in_{}".format(NEMO_run, PAR_code)

if (os.path.isdir(os.path.join(greigDir,outputsSubDir)) != True):
  os.mkdir(os.path.join(greigDir,outputsSubDir))

mixingFileName = 'SalishSea1500-RUN{}_MonthlyMean_grid_T_2D_y{}.nc'
lightFileName = "RDRS21_NEMOgrid_light_monthly_{}.nc"  


salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_y{}.nc' 
#salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_singledepth10m_{}.nc' #depth=9.5m
#salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_avgtop10m_y{}.nc'


# ///////////////////// VARS ////////////////////////
varNameLight = "solar_rad"
varNameSalinity = 'vosaline'
varNamesMixing = 'mldkz5'
lstVars = ['mldkz5']; 

#Dictionary to convert NEMO names to EwE names
#Used for both Directory name and in the file name
dctVarnameToEwEName = {varNamesMixing:'PAR-VarZ-VarK'}
saveTemplate = "{}_{}_{}.asc"


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
for var in lstVars:
    for iyear in range(startYear, endyear):

        #YearSubDirTemplate = str(iyear) # unused -GO

        #LIGHT data
        lightfullpath = os.path.join(greigDir,lightSubDir,lightFileName.format(iyear))
        dsLight = nc.Dataset(lightfullpath)
        varsLight = dsLight.variables[varNameLight]
        
        #SALINITY data
        salinityfullpath = os.path.join(greigDir,salinitySubDir,salinityFileName.format(NEMO_run,iyear))
        dsSalinity = nc.Dataset(salinityfullpath)
        varsSalinity = dsSalinity.variables[varNameSalinity]
        
        #MIXING data 
        mixingFullPath = os.path.join(greigDir, mixingSubDir, mixingFileName.format(NEMO_run, iyear))
        dsMixing = nc.Dataset(mixingFullPath)
        varsMixing = dsMixing.variables['mldkz5']       
       
        #GO2021-12-23 to-do: fix how shape is declared
        #monthlyMeanPAR = np.empty([varsSalinity.shape[2],varsSalinity.shape[3]])        
        #monthlyMeanK = np.empty([varsSalinity.shape[3],varsSalinity.shape[2]]) 
        #print(monthlyMeanK.shape)
        #print(monthlyMeanPAR.shape)


        for imon in range(1,13):
            
            #Data for month
            
            LightMon = varsLight[imon-1,:,:]
            MixingMon = varsMixing[imon-1,:,:] 
            if avg_salin == True:
              SalinityMon1 = varsSalinity[imon-1,0:salinitydepthlayer,:,:]
              SalinityMon = np.ma.average(SalinityMon1, axis = 0)
            else:
              SalinityMon = varsSalinity[imon-1,salinitydepthlayer,:,:]
            
            ar_inf = np.where(np.isinf(SalinityMon))

            # unused code? 
            #x = lambda a : a + 10
            #SalinityMon2 = x(SalinityMon)
            #print(SalinityMon2[1])
            
            # apply lin regress k w/ salin
            Ksal = a + b * SalinityMon
            Ksal[Ksal < 0.05] = 0.05
           
            #mixingZ = FixedMixingDepth
            
            #Ksal = KLinear(salinity) #Klinear RUN 102a
            #uPAR = L0*(1-math.exp(-Ksal*MixingDepth))/(Ksal*MixingDepth) #vertical avg PAR to mixing lyr depth
            step1 = np.multiply(-Ksal,MixingMon)
            step2 = 1-np.exp(step1)
            step3 = np.multiply(Ksal,MixingMon)
            step4 = np.divide(step2,step3)

            L0 = LightMon * pPAR * (1-alb)
            uPAR = np.multiply(L0,step4) #vertical avg PAR to mixing lyr depth
            uPAR[SalinityMon == 0.0] = 0.0
            Ksal[SalinityMon == 0.0] = 0.0
            uPAR = np.round(uPAR, 2)
            Ksal = np.round(Ksal, 3)

            #export PAR
            sigdigfmt = '%0.1f'
            EwEName = dctVarnameToEwEName[var]
            EwEName = EwEName.format(salinitydepthlayer) #GO ??
            savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, uPAR, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
            
            #export K
            sigdigfmt = '%0.2f'
            EwEName = "RUN" + NEMO_run + "_Kfromsal"
            #EwEName = "RUN102b-Kfromsal_singlelev10m"
            #EwEName = "RUN102a-Kfromsal_vertmeantop10m" #GO alternative
            #EwEName = EwEName.format(salinitydepthlayer) #GO ??
            savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, Ksal,bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)
            
#            #export Sal for debug 2023-03-23
#            sigdigfmt = '%0.0f'
#            EwEName = "just_salin"
#            #EwEName = "RUN102b-Kfromsal_singlelev10m"
#            #EwEName = "RUN102a-Kfromsal_vertmeantop10m" #GO alternative
#            #EwEName = EwEName.format(salinitydepthlayer) #GO ??
#            savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
#            if (os.path.isdir(savepath) != True):
#                os.mkdir(savepath)
#            savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
#            savefullpath = os.path.join(savepath, savefn)
#            print("Saving file " + savefullpath)
#            saveASCFile(savefullpath, SalinityMon,bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader, dfPlumeMask, dfLandMask)

            if imon == 12:
              break
