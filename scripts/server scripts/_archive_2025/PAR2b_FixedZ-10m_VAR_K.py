import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import os
from calendar import monthrange
import math
#import matplotlib.pyplot as plt

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO + JB 2021, last edited by GO 2023-03-23

# purpose: prep files for ECOSPACE diagnostic run 102 (fixed Z, variable K)
#          takes light from climate model, computes K from mixing depth Z, and computes integrated PAR across water depths
#          PAR = photosynthetically active radiation
#     Method codes:
#         - PAR0 - fixed mixed depth, fixed K value
#         - PAR1a - variable mixing depth, fixed K value
#         - PAR1b - variable mixed depth, fixed K value
#         - PAR2a - fixed mix depth, linear salinity->turbidity to modify K by cell, function using salinity vertmean 0-10 m
#     --> - PAR2b - fixed mix depth, linear salinity->turbidity f(n) using salinity at one depth (tested with 4.5 and 10 m)
#         - PAR3 - variable mix depth (turbocline threshold), turbidity from 2b (linear salin -> turbid f(n)
#         - PAR3b - no longer using 

#Data in:

# 1/ Mixing Layer data (var mldkz5 / turbocline depth) as monthly means in annual NC files, post-processed in other script from NEMO "xxx_MonthlyMeans_grid_T_2D_1980.nc"
# 2/ Salinity - can be either vertical averages over chosen levels or from a specific depth level, post-processed from NEMO in annual NC files "xxx_MonthlyVertMean_grid_T_selectedlevs_1980.nc"
# 3/ Light as monthly means from climate re-analysis (e.g. ERA5) with near-surface downward radiative flux (var msdwswrf Wm-2) interpolated to NEMO grid as annual NC files# 
# file name e.g. "ERA5_NEMOgrid_light_daily_1980.nc".

# Log:
# - GO-2021-12-22 edits to draw from monthly mean files rather than looping through daily NC data 
# - GO2021-12-23 differentiated RUN102a from RUN102b as vertmean of salinity (top x layers) and salinity at x m, respectively
#                Added write of K values as ASC for eval
# - GO 2022 removed loops, replaced with np. calls to dramatically increase speed
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# ////////////////// FUNCTIONS ////////////////////
#Open a Dataframe
def getDataFrame(fullPath,NaN):
    if os.path.exists(fullPath):
        nas = [NaN]
        df = pd.read_table(fullPath, skiprows=6, header=None, delim_whitespace=True,na_values=nas)
        return df
        
def buildSortableString(number, nZeros): 
    newstring = str(number)
    while (len(newstring) < nZeros) :
        tmpstr = "0" + newstring
        newstring = tmpstr
    return newstring    


# GO2021-12-22 fmt= has big effect on ASC file size, added argument for this
def saveASCFile(filename, data, dfPlumeMask, sigdigfmt):
    #Greig's bounds clipping from "Basemap Converter (Py3).ipynb"
    #I don't know why this works but it does... -JB
    #(due to how indexing is done in ASC vs NC -GO)
    trimmedData = []
    i =  bottomleft_row_ewe 
    while i < upperleft_row_ewe:
        trimmedData.append(data[i,upperleft_col_ewe:])
        i += 1
    trimmedData.reverse()
    dfEwEGrid = pd.DataFrame(data=trimmedData)
    dfEwEGrid = dfEwEGrid.mask(dfPlumeMask)
    dfEwEGrid.fillna(-9999,inplace=True)

    ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"
    np.savetxt(filename,dfEwEGrid.to_numpy(),fmt=sigdigfmt, delimiter=" ", comments='',header=ASCheader) # GO 20211222 - change from %0.5f



# //////////////// Basic Params /////////////////////
startYear = 1979 
endyear = startYear + 40  
NEMO_run = "203" #NEMO run code
PAR_code = "PAR2b" #PAR prep code

# //////////////// PAR PARAMS /////////////////////
salinitydepthlayer = 10
FixedMixingDepth = 10
a = 0.063 #GO2021-12-23
#a = -0.08529 
b = -1.94

#Variables use for Light attenuation
pPAR = 0.44 #proportion of irrad useful for phyto
alb = 0.067 #albedo reflectance
K = 0.05 #baseline K (clear water)

# /////////////////// PATHS ///////////////////////////
greigDir = "/project/6006412/goldford/ECOSPACE/"
lightSubDir = "DATA/light_monthly/"
plumeSubDir = "DATA/"
mixingSubDir = "DATA/SS1500-RUN{}/CDF/".format(NEMO_run)
salinitySubDir = "DATA/SS1500-RUN{}/NEMO_annual_NC".format(NEMO_run)
outputsSubDir ="DATA/SS1500-RUN{}/ECOSPACE_in_{}".format(NEMO_run, PAR_code)

mixingFileName = 'SalishSea1500-RUN{}_MonthlyMean_grid_T_2D_y{}.nc'
lightFileName = "ERA5_NEMOgrid_light_monthly_{}.nc"                        # = "ERA5_NEMOgrid_light_daily_2019.nc"
salinityFileName ='SalishSea1500-RUN{}_MonthlyMean_grid_T_y{}.nc' #depth=chosenabove
#salinityFileName ='SalishSea1500-RUN202_MonthlyMean_grid_T_singledepth10m_{}.nc' #depth=9.5m
#salinityFileName ='SalishSea1500-RUN202_MonthlyVertMean_grid_T_selectedlevs_{}.nc' #alternative

# /////////////////// VARS ///////////////////////////
varNameLight = "msdwswrf"
varNameSalinity = 'vosaline'
varNamesMixing = 'mldkz5'
lstVars = ['mldkz5']; #GO2021-12-22 ??

#Dictionary to convert NEMO names to EwE names
#Used for both Directory name and in the file name
dctVarnameToEwEName = {varNamesMixing:'LightAttenuation_PAR-FixedZ-VarK'}
saveTemplate = "{}_{}_{}.asc"


#Number of depth strata to use
#Depth indexes are zero based
#So 10 will be the 9 index
#nDepthIndexes = 10 #= 9.001736   10.003407
#nDepthIndexes = 23 #= 31.101034   39.11886
#nDepthIndexes = 26 #= 86.96747   109.73707 

# //////////////// MAP CLIP /////////////////////
#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39

# //////////////// PLUME MASK //////////////////////
#Create the plume mask
PlumeRegionFilename = os.path.join(greigDir, plumeSubDir, "Fraser_Plume_Region.asc")
dfPlumeRegion = getDataFrame(PlumeRegionFilename,"-9999.00000000")
dfPlumeMask = dfPlumeRegion == 1 #Plume is region one



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
        #dsMixing = nc.Dataset(mixingFullPath)
        #varsMixing = dsMixing.variables['mldkz5']       
       
        #GO2021-12-23 to-do: fix how shape is declared
        monthlyMeanPAR = np.empty([varsSalinity.shape[1],varsSalinity.shape[2]])        
        monthlyMeanK = np.empty([varsSalinity.shape[1],varsSalinity.shape[2]]) 


        for imon in range(1,13):
            
            #Data for month
            LightMon = varsLight[imon-1,:,:]
            #MixingMon = varsMixing[imon-1,:,:] 
            SalinityMon = varsSalinity[imon-1,salinitydepthlayer,:,:]
            

# removed loops used in JB's version - should not need them!

            x = lambda a : a + 10
            SalinityMon2 = x(SalinityMon)
            #print(SalinityMon2[1])
            
            # apply lin regress k w/ salin
            Ksal = a * SalinityMon + b
            Ksal[Ksal < 0.05] = 0.05
            mixingZ = FixedMixingDepth
            step1 = -Ksal*mixingZ  
            step2 = 1-np.exp(step1)                
            step3 = Ksal*mixingZ
            step4 = np.divide(step2,step3)
            
            L0 = LightMon * pPAR * (1-alb)
            uPAR = np.multiply(L0,step4) #vertical avg PAR to mixing lyr depth
            uPAR[SalinityMon == 0.0] = 0.0
            Ksal[SalinityMon == 0.0] = 0.0
            uPAR = np.round(uPAR, 2)
            Ksal = np.round(Ksal, 3)
            print(uPAR[1])

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
            saveASCFile(savefullpath, uPAR, dfPlumeMask,sigdigfmt)
            
            #export K
            sigdigfmt = '%0.2f'
            EwEName = "RUN203_Kfromsal_singlelev4m"
            #EwEName = "RUN102b-Kfromsal_singlelev10m"
            #EwEName = "RUN102a-Kfromsal_vertmeantop10m" #GO alternative
            #EwEName = EwEName.format(salinitydepthlayer) #GO ??
            savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, Ksal, dfPlumeMask,sigdigfmt)



