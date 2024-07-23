import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import os
from calendar import monthrange
import math
import matplotlib.pyplot as plt

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#Data requirements

#The Mixing Layer data are the monthly means, one file contain all the data for one year. "xxx_MonthlyMeans_grid_T_2D_1980.nc"

#Salinity are the daily means. One file contains data for one day, with the year month day timestamp in the file name "xxx_1d_grid_T_19791228-19791228.nc"
#We are averaging the salinity over the top 10 layers/meters

#Light are the daily values, one file contains data for one year. The year time stamp is in the file name "ERA5_NEMOgrid_light_daily_1980.nc".
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#
#/project/6006412/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN202/CDF/[year]/
#
#rootDir = "/project/6006412/mdunphy/nemo_results/SalishSea1500"
rootDir = "W:\\Sync\\PSF\\EwE\\Georgia Strait 2021\\LTL_model\\DATA"
InputsSubDir = 'SalishSea1500-RUN202\\CDF'
YearSubDirTemplate = ""
LightSubDir = "ERA5_light_NC"
outputsSubDir ="SalishSea1500-RUN202\\ECOSPACE_in_RUN102"

lstVars = ['mldkz5'];

MixingFileName = 'SalishSea1500-RUN201_MonthlyMeans_grid_T_2D_{}.nc'
MixingVarName = 'mldkz5'

#Dictionary to convert NEMO names to EwE names
#Used for both Directory name and in the file name
dctVarnameToEwEName = {MixingVarName:'LightAttenuation_PAR-FixedZ-VarK'}

fnTempLight = "ERA5_NEMOgrid_light_daily_{}.nc"# = "ERA5_NEMOgrid_light_daily_2019.nc"
varNameLight = "msdwswrf"

saveTemplate = "{}_{}_{}.asc"

SalinityTemplate ='SalishSea1500-RUN200_1d_grid_T_{}{}{}-{}{}{}.nc'
SalinityVarname = 'vosaline'

iLightDay = 0

startYear = 1979# Light data is just 2019 for development
endyear = startYear + 1 # 1980 

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
PAR = 0.44
a = 0.067
K = 0.05
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#
salinitydepthlayer = 10
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


def saveASCFile(filename, data, dfPlumeMask):
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
    np.savetxt(filename,dfEwEGrid.to_numpy(),fmt='%0.5f', delimiter=" ", comments='',header=ASCheader)


def lightAttenuation(LightWM2, MixingDepth, salinity):
    L0 = LightWM2 * PAR * (1-a)
    Ksal = KLinear(salinity)
    return L0*(1-math.exp(-Ksal*MixingDepth))/(Ksal*MixingDepth) 


def getSalintyMatrix(Year, iMon, iDay, iSalinityDepthLayer):
    m = buildSortableString(iMon,2)
    d = buildSortableString(iDay,2)
    filename = SalinityTemplate.format(Year,m,d,Year,m,d)
    SalinitySubDir = YearSubDirTemplate.format(Year)
    fullpath = os.path.join(rootDir,InputsSubDir, SalinitySubDir,filename)

    ds = nc.Dataset(fullpath)
    vars = ds.variables[SalinityVarname]
    #for now get the average of the top iSalinityDepthLayer layers
    Salinity = vars[:,0:iSalinityDepthLayer,:].mean(axis=1)   
    return Salinity


def KLinear(Salinity):
    #This is based on Loos et al 2017
    a = -0.08529 
    b = 2.60
    k = a * Salinity + b
    if(k<0.05):
        k=0.05
    return k

#Create the plume mask
PlumeRegionFilename = os.path.join(rootDir,"Fraser_Plume_Region.asc")
dfPlumeRegion = getDataFrame(PlumeRegionFilename,"-9999.00000000")
#Plume is region one
dfPlumeMask = dfPlumeRegion == 1

#assert False,"Warning Development version! Did you remove the code to match Light year 1980 to the Salinity Year 1979?"

for var in lstVars:
    for iyear in range(startYear, endyear):
        iLightDay=-1
        #get the light data for this year
        lightfullpath = os.path.join(rootDir,LightSubDir, fnTempLight.format(iyear))
        dsLight = nc.Dataset(lightfullpath)
        varsLight = dsLight.variables[varNameLight]

        YearSubDirTemplate = str(iyear)
        #Mixing Layer 
        mixingFullPath = os.path.join(rootDir,InputsSubDir, YearSubDirTemplate, MixingFileName.format(iyear))
        #dsMixing = nc.Dataset(mixingFullPath)
        #varsMixing = dsMixing.variables['mldkz5']       
       
        for imon in range(1,13):
            #Mixing data is the monthly mean so we just need to get one matrix per month
            #MixingMon = varsMixing[imon-1,:,:]

            ndays = monthrange(iyear,imon)
            bnewMonth = True 

            for iday in  range(1,ndays[1]+1):
                iLightDay+=1
                print("Year = "+ str(iyear) + ", Month = " + str(imon)+ ", Day = " + str(iLightDay))

                #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                #HACK for testing
                #Salinity Data is only for one year 1979
                #Mixing Data is for 1980 so tweak the salinity year to match
                #Also 1980 is a leap year! so skip 1980 Feb 29
                #for production data remove the year-1 and Feb 29 check
                #if ((imon != 2) and  (iday != 29)):
                    #get average salinity over the top 10 meters
                   # print("HACK WARNING for debugging Salinity data is only 1979 Light data is 1980 FIX THIS FOR PRODUCTION CODE" )
                salinity = getSalintyMatrix(iyear,imon,iday,salinitydepthlayer)
                #XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

                if(bnewMonth):
                    dailymeans = np.empty([ndays[1],salinity.shape[1],salinity.shape[2]])
                    bnewMonth = False

                for irow in range(salinity.shape[1]):
                    for icol in range(salinity.shape[2]):

                        #get the mixing layer depth for this row col
                        mixingDepth = FixedMixingDepth# MixingMon[irow,icol]
                        #Salinity for this day
                        sal = salinity[0,irow,icol]

                        if (mixingDepth != 0) & (iLightDay < varsLight.shape[0]):
                            dailymeans[iday-1,irow,icol] = lightAttenuation(varsLight[iLightDay,irow,icol],mixingDepth, sal)

            monthlyMean = dailymeans.mean(axis=0)

            #Save the monthly means to an .asc file
            EwEName = dctVarnameToEwEName[var]
            EwEName = EwEName.format(salinitydepthlayer)
            savepath = os.path.join(rootDir,outputsSubDir,EwEName ) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            
            savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
            savefullpath = os.path.join(savepath, savefn)

            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, monthlyMean, dfPlumeMask)



