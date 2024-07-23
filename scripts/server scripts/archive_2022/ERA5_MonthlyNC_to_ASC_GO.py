import netCDF4 as nc
import numpy as np
import pandas as pd
import os

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# created by GO + JB 2021, last edited by GO 2021-12-23

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

def buildSortableString(number, nZeros): 
    newstring = str(number)
    while (len(newstring) < nZeros) :
        tmpstr = "0" + newstring
        newstring = tmpstr
    return newstring   


def saveASCFile(filename, data):
    #Greig's bounds clipping from "Basemap Converter (Py3).ipynb"
    #I don't know why this works but it does...
    trimmedData = []
    i =  bottomleft_row_ewe 
    while i < upperleft_row_ewe:
        trimmedData.append(data[i,upperleft_col_ewe:])
        i += 1
    trimmedData.reverse()

    np.savetxt(filename,trimmedData,fmt='%0.5f', delimiter=" ", comments='',header=ASCheader)


#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39
ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"

rootDir = "/project/6006412/goldford/ECOSPACE/DATA/"
#inputsSubDir = 'CDF/1980' # inserted in loop - GO 20210917
saveTemplate = "ERA5_NEMOgrid_{}_{}_{}.asc"
saveDir = "/project/6006412/goldford/ECOSPACE/DATA/ASC/"
inputsSubDir = "light_monthly/"

#Variable names
#For simplicity use the variables names from the NEMO Model
#If we want easier to read name i.e. Salintity, UpWelling... we will need to create another dictionary
#with the EwE dumbass names!
lstVars = ['msdwswrf'];
#The key in the dictionary need to match the variable names in lstVars
dctFileTemplates = {'msdwswrf':'ERA5_NEMOgrid_light_monthly_{}.nc'}
#Unit Conversions - a multiplier to help with Ecospace parameterization
#This doesn't matter because we create a response function on the values so that can look like anything
dctUnitConvert = {'msdwswrf':1}

startYear = 1979
endyear = startYear + 41

for var in lstVars:
    for iyear in range(startYear, endyear):
    # input format: one file per year with daily data
        inputTemplate = dctFileTemplates[var]
        filename = inputTemplate.format(iyear) #fixed GO 20210917
        fullpath = os.path.join(rootDir, inputsSubDir,filename)
        print(fullpath)
        
        if(os.path.isfile(fullpath)):
            ds = nc.Dataset(fullpath)
            vars = ds.variables[var]
        else:
            print("Could not find NC file. " + fullpath)

        for imon in range(1,13):
            m = buildSortableString(imon,2)    
            print("month: " + str(imon))

            monthlydata = vars[imon-1]
            print(monthlydata)
            print(monthlydata.shape)

            #Save the monthly means to an .asc file
            savepath = os.path.join(saveDir, var) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)

            savefn = saveTemplate.format(var,iyear, m)
            saveASCFile(os.path.join(savepath, savefn) ,monthlydata)


        




