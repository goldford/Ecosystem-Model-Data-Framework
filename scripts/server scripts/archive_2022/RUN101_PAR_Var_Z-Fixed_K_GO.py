import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import os
from calendar import monthrange
import math
#import matplotlib.pyplot as plt #not needed - GO

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO + JB 2021

# purpose: prep files for ECOSPACE diagnostic run 101 (variable Z, fixed K)

#Data in:

# 1/ Mixing Layer data (var mldkz5 / turbocline depth) as monthly means in annual NC files, post-processed from NEMO "xxx_MonthlyMeans_grid_T_2D_1980.nc"

# 2/ Light as monthly means from climate re-analysis (e.g. ERA5) with near-surface downward radiative flux (var msdwswrf Wm-2) interpolated to NEMO grid as annual NC files# 

# file name e.g. "ERA5_NEMOgrid_light_daily_1980.nc".
# notes  
#  - tested with Python 3.8.10 on graham Dec 2021 (first: 'module load python/3.8.10' then 'source py3-8-10_Jul9/bin/activate')
#  - run in virtual env on goldford account
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#
#/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN202/CDF/[year]/
michaelDir = "/project/6006412/mdunphy/nemo_results/SalishSea1500"
#rootDir = "W:\\Sync\\PSF\\EwE\\Georgia Strait 2021\\LTL_model\\DATA"
inputsSubDir = "SalishSea1500-RUN202/NEMO_annualNC"
#InputsSubDir = 'SalishSea1500-RUN202\\CDF'

YearSubDirTemplate = ""
greigDir = "/project/6006412/goldford/ECOSPACE/"
lightSubDir = "DATA/light_monthly/"
plumeSubDir = "DATA/"
#LightSubDir = "ERA5_light_NC"
mixingSubDir = "DATA/SS1500-RUN202/NEMO_annualNC/" # added GO (files created from .ksh)
outputsSubDir ="DATA/SS1500-RUN202/ECOSPACE_in_RUN101"
#outputsSubDir ="SalishSea1500-RUN202\\ECOSPACE_in_RUN101"

lstVars = ['mldkz5'];

MixingFileName = 'SalishSea1500-RUN202_MonthlyMean_grid_T_2D_{}.nc'
MixingVarName = 'mldkz5'

#Dictionary to convert NEMO names to EwE names
#Used for both Directory name and in the file name
dctVarnameToEwEName = {MixingVarName:'LightAttenuation_PAR-VarZ-FixedK'}

fnTempLight = "ERA5_NEMOgrid_light_monthly_{}.nc"# = "ERA5_NEMOgrid_light_daily_2019.nc"
varNameLight = "msdwswrf"

saveTemplate = "{}_{}_{}.asc"

# unused
SalinitySubDir = ''
SalinityTemplate ='SalishSea1500-RUN202_1d_grid_T_{}{}{}-{}{}{}.nc'
SalinityVarname = 'vosaline'

iLightDay = 0

startYear = 1979
endyear = startYear + 40 #  

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
a = 0.067 #albedo reflectance
K = 0.5 #baseline K (0.05 = clearest ocean water) - revised 2022-01-14 to be more realistic (0.5 based on Loos et al 2017 Tab. 3)
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#
salinitydepthlayer = 0 # not used in Ecospace diagnostic run 101
FixedK = 0.05 # this is probably far too low


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


# GO-2021-12-22 fmt= has a lot of effect on ASC file size
def saveASCFile(filename, data, dfMask):
    #Greig's bounds clipping from "Basemap Converter (Py3).ipynb"
    #I don't know why this works but it does...
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
    np.savetxt(filename,dfEwEGrid.to_numpy(),fmt='%0.1f', delimiter=" ", comments='',header=ASCheader)


def lightAttenuation(LightWM2, MixingDepth, salinity):
    L0 = LightWM2 * pPAR * (1-a)
    Ksal = FixedK#KLinear(salinity)
    return L0*(1-math.exp(-Ksal*MixingDepth))/(Ksal*MixingDepth) 


def getSalintyMatrix(Year, iMon, iDay, iSalinityDepthLayer):
    m = buildSortableString(iMon,2)
    d = buildSortableString(iDay,2)
    filename = SalinityTemplate.format(Year,m,d,Year,m,d)
    SalinitySubDir = YearSubDirTemplate.format(Year)
    fullpath = os.path.join(michaelDir, SalinitySubDir,filename)

    ds = nc.Dataset(fullpath)
    vars = ds.variables[SalinityVarname]
    #for now get the average of the top iSalinityDepthLayer layers
    Salinity = vars[:,0:iSalinityDepthLayer,:].mean(axis=1)   
    return Salinity

# not used in this - 
def KLinear(Salinity):
    a = -0.0598 
    b = 1.84
    k = a * Salinity + b
    if(k<0.05):
        k=0.05
    return k

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#For Debugging
#Plot the Salinity Modifier
#plotSalinityModifiers()
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#Create the plume mask
PlumeRegionFilename = os.path.join(greigDir, plumeSubDir,"Fraser_Plume_Region.asc")
#PlumeRegionFilename = os.path.join(rootDir,"Fraser_Plume_Region.asc")
dfPlumeRegion = getDataFrame(PlumeRegionFilename,"-9999.00000000")

#Plume is region one
dfPlumeMask = dfPlumeRegion == 1


#assert False,"Warning Development version! Did you remove the code to match Light year 1980 to the Salinity Year 1979?"

for var in lstVars:
    for iyear in range(startYear, endyear):
        
        #get the light data for this year
        
        # GO-2021-12-22 changed to use monthly light instead of daily (faster)
        lightfullpath = os.path.join(greigDir,lightSubDir, fnTempLight.format(iyear))
        #lightfullpath = os.path.join(lightDir,LightSubDir, fnTempLight.format(iyear)) #path change GO
        dsLight = nc.Dataset(lightfullpath)
        varsLight = dsLight.variables[varNameLight]

        YearSubDirTemplate = str(iyear)
        
        #Mixing Layer 
        mixingFullPath = os.path.join(greigDir, mixingSubDir, MixingFileName.format(iyear)) # GO 20211221
        #mixingFullPath = os.path.join(mixingDir,InputsSubDir, YearSubDirTemplate, MixingFileName.format(iyear))
        dsMixing = nc.Dataset(mixingFullPath)
        varsMixing = dsMixing.variables['mldkz5']       
       
        # GO-2021-12-22
        monthlyMean = np.empty([varsMixing.shape[1],varsMixing.shape[2]])
       
        for imon in range(1,13):
            #Mixing data is the monthly mean so we just need to get one matrix per month
            MixingMon = varsMixing[imon-1,:,:]
            
            #GO-2021-12-22 changed use monthly light means instead of looping daily
            LightMon = varsLight[imon-1,:,:]

            ndays = monthrange(iyear,imon)
            bnewMonth = True 

            # GO-2021-12-22 change this loop below to monthly instead of daily 
            for irow in range(varsMixing.shape[1]):
              for icol in range(varsMixing.shape[2]):

                #get the mixing layer depth for this row col
                mixingDepth =  MixingMon[irow,icol]
                
                #Salinity for this day
                sal = np.nan#Don't need salinity for Fixed K salinity[0,irow,icol]

                if (mixingDepth != 0): #GO-2021-12-22
                  monthlyMean[irow,icol] = lightAttenuation(LightMon[irow,icol],mixingDepth, sal)
                else:
                  monthlyMean[irow,icol] = 0.0 #GO-2021-12-22 can't use np.empty values on server py env
                
                # chec
                if (monthlyMean[irow,icol] > 200):
                  print("Warning: strangely high PAR value")
                  print(monthlyMean[irow,icol])

            #Save the monthly means to an .asc file
            EwEName = dctVarnameToEwEName[var]
            EwEName = EwEName.format(salinitydepthlayer)
            savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            
            savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
            savefullpath = os.path.join(savepath, savefn)

            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, monthlyMean, dfPlumeMask)



