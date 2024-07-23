import netCDF4 as nc
import numpy as np
import pandas as pd
import os
from calendar import monthrange

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO + JB 2021, last edited by GO 2022-01-20
#
# Purpose: convert pre-processed NC files (mean monthly NEMO results in annual NetCDFs) to ASC for ECOSPACE
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def buildSortableString(number, nZeros): 
    newstring = str(number)
    while (len(newstring) < nZeros) :
        tmpstr = "0" + newstring
        newstring = tmpstr
    return newstring    

# GO2021-12-22 fmt= has a lot of effect on ASC file size, added argument for this
def saveASCFile(filename, data, sigdigfmt):
    #Greig's bounds clipping from "Basemap Converter (Py3).ipynb"
    #I don't know why this works but it does...
    trimmedData = []
    i =  bottomleft_row_ewe 
    while i < upperleft_row_ewe:
        trimmedData.append(data[i,upperleft_col_ewe:])
        i += 1
    trimmedData.reverse()

    np.savetxt(filename,trimmedData,fmt=sigdigfmt, delimiter=" ", comments='',header=ASCheader)

greigDir = "/project/6006412/goldford/ECOSPACE/"
NEMOAnnualNCSubDir = "DATA/SS1500-RUN202/NEMO_annualNC/"
#inputsSubDir = 'CDF/1980' # inserted in loop - GO 20210917
saveTemplate = "SalishSea1500-RUN202_{}_{}_{}.asc"
outputsSubDir = "DATA/SS1500-RUN202/NEMOout_as_ASC"

startYear = 1979
endyear = startYear + 40
sigdigfmt = '%0.1f' #GO2021-12-23 sets decimal places in out ASC e.g 0.2f is two decimals

#Variable names
#For simplicity use the variables names from the NEMO Model
lstVars = ['vozocrtx','vomecrty']; # updating veloc fields - GO 2022-01-26
#lstVars = ['votemper', 'vosaline']; # temp + salin only run for 1m data - GO 2022-01-20
#lstVars = ['votemper', 'vosaline', 'mldr10_1', 'mldkz5'];
#lstVars = ['votemper', 'vosaline','vovecrtz','vozocrtx','vomecrty'];

dctVarnameToEwEName = {'vozocrtx':'xvelocMean10m','vomecrty':'yvelocMean10m'}
#dctVarnameToEwEName = {'votemper':'TempVertMean10m','vosaline':'SalinVertMean10m'} #alt
#dctVarnameToEwEName = {'votemper':'TempAt1m','vosaline':'SalinAt1m'}
#dctVarnameToEwEName = {'votemper':'TempAt4m','vosaline':'SalinAt4m','mldr10_1':'MixedLyrZ','mldkz5':'MixingTurboZ'}
#dctVarnameToEwEName = {'votemper':'TempAt10m','vosaline':'SalinAt10m'}

dctFileTemplates = {'vozocrtx':'SalishSea1500-RUN202_MonthlyVertMean_grid_U_selectedlevs_{}.nc','vozocrty':'SalishSea1500-RUN202_MonthlyVertMean_grid_V_selectedlevs_{}.nc'}
#dctFileTemplates = {'votemper':'SalishSea1500-RUN202_MonthlyMean_grid_T_singledepth1m_{}.nc','vosaline':'SalishSea1500-RUN202_MonthlyMean_grid_T_singledepth1m_{}.nc'}
#dctFileTemplates = {'votemper':'SalishSea1500-RUN202_MonthlyMean_grid_T_singledepth4m_{}.nc','vosaline':'SalishSea1500-RUN202_MonthlyMean_grid_T_singledepth4m_{}.nc','mldr10_1': 'SalishSea1500-RUN202_MonthlyMean_grid_T_2D_{}.nc', 'mldkz5': 'SalishSea1500-RUN202_MonthlyMean_grid_T_2D_{}.nc'}
#dctFileTemplates = {'votemper':'SalishSea1500-RUN202_MonthlyMean_grid_T_singledepth10m_{}.nc','vosaline':'SalishSea1500-RUN202_MonthlyMean_grid_T_singledepth10m_{}.nc'}


#Unit Conversions
#'vovecrtz':1000 upwelling/downwelling (vertical velocity) is in meters per second but his makes it really hard parameterize the EwE UI
#so arbitrarily scale them bigger  - JB
#I changed it so multiplier is 100 to convert from m/s (NEMO) to cm/s (EwE) - GO
dctUnitConvert = {'vozocrtx':100,'vomecrty':100}
#dctUnitConvert = {'votemper':1,'vosaline':1}
#dctUnitConvert = {'votemper':1,'vosaline':1,'mldr10_1':1,'mldkz5':1,'vovecrtz':100,'vozocrtx':100,'vomecrty':100}

#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39

ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"

for var in lstVars:
    for iyear in range(startYear, endyear):
        #inputsSubDir = 'CDF/' + str(iyear)
        
        #data for year
        inputTemplate = dctFileTemplates[var]
        datafullpath = os.path.join(greigDir,NEMOAnnualNCSubDir,inputTemplate.format(iyear))
        dsData = nc.Dataset(datafullpath)
        varsData = dsData.variables[var]
        
        for imon in range(1,13):
            
            #Data for month
            DataMon = varsData[imon-1,:,:]
    
            #Get the unit conversion multiplier
            UnitConvert = dctUnitConvert[var]
            DataMon = DataMon * UnitConvert #GO2021-12-23 - will this work since it's an array?
    
            #Save the monthly Data to an .asc file
            
            EwEName = dctVarnameToEwEName[var]
            #EwEName = EwEName.format(salinitydepthlayer) #GO ??
            savepath = os.path.join(greigDir,outputsSubDir,EwEName ) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)

            savefn = saveTemplate.format(EwEName,iyear,buildSortableString(imon,2))
            savefullpath = os.path.join(savepath, savefn)

            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, DataMon, sigdigfmt)



