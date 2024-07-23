import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import os
from calendar import monthrange
import math
#import matplotlib.pyplot as plt

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO + JB 2021, last edited by GO 2021-12-23

# purpose: prep files for ECOSPACE diagnostic run 102 (fixed Z, variable K)
#         - RUN102a - linear salinity->turbidity (i.e. K) function using salinity vertmean 0-10 m, RUN102b - linear salinity->turbidity f(n) using 10m depth lyr

#Data in:

# 1/ Mixing Layer data (var mldkz5 / turbocline depth) as monthly means in annual NC files, post-processed from NEMO "xxx_MonthlyMeans_grid_T_2D_1980.nc"

# 2/ Salinity - can be either vertical averages over chosen levels or from a specific depth level, post-processed from NEMO in annual NC files "xxx_MonthlyVertMean_grid_T_selectedlevs_1980.nc"

# 3/ Light as monthly means from climate re-analysis (e.g. ERA5) with near-surface downward radiative flux (var msdwswrf Wm-2) interpolated to NEMO grid as annual NC files# 
# file name e.g. "ERA5_NEMOgrid_light_daily_1980.nc".

# Log:
# - GO-2021-12-22 edits to draw from monthly mean files rather than looping through daily NC data 
# - GO2021-12-23 differentiated RUN102a from RUN102b as vertmean of salinity (top ten layers, 10 m max) and salinity at 10 m, respectively
#                Added write of K values as ASC for eval

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

michaelDir = "/project/6006412/mdunphy/nemo_results/SalishSea1500"
greigDir = "/project/6006412/goldford/ECOSPACE/"

inputsSubDir = "SalishSea1500-RUN202/CDF" # GO-2021-12-22 unused - raw NEMO files are all pre-processed now
YearSubDirTemplate = "" #used if NEMO results are organized into yearly folders

lightSubDir = "DATA/light_monthly/"
plumeSubDir = "DATA/"
mixingSubDir = "DATA/SS1500-RUN202/CDF/" 
salinitySubDir = "DATA/SS1500-RUN202/NEMO_annualNC"
outputsSubDir ="DATA/SS1500-RUN202/ECOSPACE_in_RUN102a"
#outputsSubDir ="DATA/SS1500-RUN202/ECOSPACE_in_RUN102b"

mixingFileName = 'SalishSea1500-RUN202_MonthlyMean_grid_T_2D_{}.nc'
lightFileName = "ERA5_NEMOgrid_light_monthly_{}.nc"                        # = "ERA5_NEMOgrid_light_daily_2019.nc"

#salinityFileName ='SalishSea1500-RUN202_MonthlyMean_grid_T_singledepth4m_{}.nc' #depth=4.5m
#salinityFileName ='SalishSea1500-RUN202_MonthlyMean_grid_T_singledepth10m_{}.nc' #depth=9.5m
salinityFileName ='SalishSea1500-RUN202_MonthlyVertMean_grid_T_selectedlevs_{}.nc' #alternative
#salinityFileName ='SalishSea1500-RUN202_MonthlyVertMean_grid_T_{}.nc' #OVERWRITE TEMPORARY UNTIL KSH RERUN TO FIX NAME 

varNameLight = "msdwswrf"
varNameSalinity = 'vosaline'
varNamesMixing = 'mldkz5'
lstVars = ['mldkz5']; #GO2021-12-22 ??

#Dictionary to convert NEMO names to EwE names
#Used for both Directory name and in the file name
dctVarnameToEwEName = {varNamesMixing:'LightAttenuation_PAR-FixedZ-VarK'}

saveTemplate = "{}_{}_{}.asc"

startYear = 1979# Light data is just 2019 for development
endyear = startYear + 30  

#Number of depth strata to use
#Depth indexes are zero based
#So 10 will be the 9 index
#nDepthIndexes = 10 #= 9.001736   10.003407
#nDepthIndexes = 23 #= 31.101034   39.11886
#nDepthIndexes = 26 #= 86.96747   109.73707 

#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#Variables use for Light attenuation
pPAR = 0.44 #proportion of irrad useful for phyto
alb = 0.067 #albedo
K = 0.05 #baseline K (clear water)
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

salinitydepthlayer = 10 #GO 2021-12-22 unused now
FixedMixingDepth = 10 

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


# GO2021-12-22 fmt= has a lot of effect on ASC file size, added argument for this
def saveASCFile(filename, data, dfPlumeMask, sigdigfmt):
    #Greig's bounds clipping from "Basemap Converter (Py3).ipynb"
    #I don't know why this works but it does... -JB
    #(it's weird but due to how indexing is done in ASC vs NC -GO)
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

# GO2021-12-23 modified to return K
def lightAttenuation(LightWM2, MixingDepth, salinity):
    L0 = LightWM2 * pPAR * (1-alb)
    Ksal = KLinear(salinity) #Klinear RUN 102a
    uPAR = L0*(1-math.exp(-Ksal*MixingDepth))/(Ksal*MixingDepth) #vertical avg PAR to mixing lyr depth
    return uPAR, Ksal 


# GO2021-12-22 function not used anymore
def getSalintyMatrix(Year, iMon, iDay, iSalinityDepthLayer):
    m = buildSortableString(iMon,2)
    d = buildSortableString(iDay,2)
    filename = SalinityTemplate.format(Year,m,d,Year,m,d)
    SalinitySubDir = YearSubDirTemplate.format(Year)
    fullpath = os.path.join(greigDir,salinitySubDir,filename)

    ds = nc.Dataset(fullpath)
    vars = ds.variables[SalinityVarname]
    #for now get the average of the top iSalinityDepthLayer layers
    Salinity = vars[:,0:iSalinityDepthLayer,:].mean(axis=1)   
    return Salinity


def KLinear(Salinity):
    #This is based on Loos et al 2017 -GO
    a = -0.063 #GO2021-12-23
    #a = -0.08529 
    b = 1.94
    k = a * Salinity + b
    if(k<0.05):
        k=0.05
    return k


#Create the plume mask
PlumeRegionFilename = os.path.join(greigDir, plumeSubDir, "Fraser_Plume_Region.asc")
dfPlumeRegion = getDataFrame(PlumeRegionFilename,"-9999.00000000")
#Plume is region one
dfPlumeMask = dfPlumeRegion == 1


# main loops
for var in lstVars:
    for iyear in range(startYear, endyear):

        YearSubDirTemplate = str(iyear) # GO2021-12-22 unused

        #LIGHT data for this year
        lightfullpath = os.path.join(greigDir,lightSubDir,lightFileName.format(iyear))
        dsLight = nc.Dataset(lightfullpath)
        varsLight = dsLight.variables[varNameLight]
        
        #SALINITY data
        salinityfullpath = os.path.join(greigDir,salinitySubDir,salinityFileName.format(iyear))
        dsSalinity = nc.Dataset(salinityfullpath)
        varsSalinity = dsSalinity.variables[varNameSalinity]
        
        #GO2021-12-23 (unused in RUN 102)
        #MIXING data 
        mixingFullPath = os.path.join(greigDir, mixingSubDir, mixingFileName.format(iyear)) # GO20211221
        #dsMixing = nc.Dataset(mixingFullPath)
        #varsMixing = dsMixing.variables['mldkz5']       
       
        #GO2021-12-23 to-do: fix how shape is declared
        monthlyMeanPAR = np.empty([varsSalinity.shape[1],varsSalinity.shape[2]])        
        monthlyMeanK = np.empty([varsSalinity.shape[1],varsSalinity.shape[2]]) 

        for imon in range(1,13):
            
            #Data for month
            LightMon = varsLight[imon-1,:,:]
            #MixingMon = varsMixing[imon-1,:,:] #RUN102 
            SalinityMon = varsSalinity[imon-1,:,:]

            # GO2021-12-22 - get rid of these two loops
            for irow in range(varsSalinity.shape[1]):
              for icol in range(varsSalinity.shape[2]):
                
                mixingZ = FixedMixingDepth #RUN102 - fixed mixing depth (Z)
                #mixingZ = MixingMon[irow,icol]
                sal =  SalinityMon[irow,icol]  # K
                light =  LightMon[irow,icol]

                if (sal != 0): #GO-2021-12-22 - check here should be for land / water mask
                  fMonthlyPAR, k_ = lightAttenuation(light,mixingZ, sal)
                  monthlyMeanPAR[irow,icol] = round(fMonthlyPAR,2)
                  monthlyMeanK[irow,icol] = round(k_,2)
                else:
                  monthlyMeanPAR[irow,icol] = 0.0 #GO2021-12-22 can't use np.empty values on server py env
                  monthlyMeanK[irow,icol] = 0.0
                
                # GO20211222 check
                if (monthlyMeanPAR[irow,icol] > 200):
                  print("Warning: strangely high PAR value")
                  print(monthlyMeanPAR[irow,icol])

            #Save the monthly avgPAR to an .asc file
            sigdigfmt = '%0.1f'
            EwEName = dctVarnameToEwEName[var]
            EwEName = EwEName.format(salinitydepthlayer) #GO ??
            savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            
            savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
            savefullpath = os.path.join(savepath, savefn)

            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, monthlyMeanPAR, dfPlumeMask,sigdigfmt)
            
            #GO2021-12-23 added to export K as ASC, too, for eval
            sigdigfmt = '%0.2f'
            
            #EwEName = "RUN102a-Kfromsal_singlelev4m"
            EwEName = "RUN102a-Kfromsal_vertmean10m"
            #EwEName = EwEName.format(salinitydepthlayer) #GO ??
            savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            
            savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
            savefullpath = os.path.join(savepath, savefn)

            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, monthlyMeanK, dfPlumeMask,sigdigfmt)



