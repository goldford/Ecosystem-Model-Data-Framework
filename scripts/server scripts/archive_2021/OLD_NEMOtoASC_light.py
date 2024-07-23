import netCDF4 as nc
import numpy as np
import pandas as pd
import os
from calendar import monthrange


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

#Variable names
#For simplicity use the variables names from the NEMO Model
#If we want easier to read name i.e. Salintity, UpWelling... we will need to create another dictionary
#with the EwE dumbass names!
lstVars = ['msdwswrf'];

#The key in the dictionary need to match the variable names in lstVars
dctFileTemplates = {'msdwswrf':'ERA5_NEMOgrid_light_daily_{}.nc'}

#Unit Conversions - a multiplier to help with Ecospace parameterization 
#This doesn't matter because we create a response function on the values so that can look like anything 
dctUnitConvert = {'msdwswrf':1}

startYear = 1979
endyear = startYear + 4

#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39

rootDir = "/project/6006412/goldford/ECOSPACE/DATA/"
#inputsSubDir = 'CDF/1980' # inserted in loop - GO 20210917
saveTemplate = "ERA5_NEMOgrid_{}{}{}.asc"
saveDir = "/project/6006412/goldford/ECOSPACE/DATA/ASC/"

ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"

for var in lstVars:
    for iyear in range(startYear, endyear):
        
        # input format: one file per year with daily data
        inputsSubDir = 'light_daily/'
        inputTemplate = dctFileTemplates[var]
        filename = inputTemplate.format(iyear) #fixed GO 20210917
        fullpath = os.path.join(rootDir, inputsSubDir,filename)
        print(fullpath)
        dayofyear = 0

        if(os.path.isfile(fullpath)):
          
          ds = nc.Dataset(fullpath)
          vars = ds.variables[var]

          for imon in range(1,13):
            
            print("month")
            print(imon)

            ndays = monthrange(iyear,imon)
            
            # monthly and daily values empty arrays
            monthlyMean = np.empty([vars.shape[1],vars.shape[2]])
            dailymeans = np.empty([ndays[1],vars.shape[1],vars.shape[2]])
            
            bnewMonth = True 
            for iday in  range(1,ndays[1]+1):
                
                dayofyear += 1
                
                m = buildSortableString(imon,2)
                d = buildSortableString(iday,2)
              
                #print(vars[dayofyear-1])

                if(bnewMonth):
                  dailymeans = np.empty([ndays[1],vars.shape[1],vars.shape[2]])
                  bnewMonth = False
                    
                  # put the daily values into array of days for the month
                  dailymeans[iday-1] = vars[dayofyear-1]

        
            #Get the unit conversion multiplier
            UnitConvert = dctUnitConvert[var]
            
            #get the mean across the days for this month
            monthlyMean = dailymeans.mean(axis=0)
            monthlyMean = monthlyMean * UnitConvert
            print(monthlyMean)

            #Save the monthly means to an .asc file
            savepath = os.path.join(saveDir, var) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            
            savefn = saveTemplate.format(var,iyear,m)
            saveASCFile(os.path.join(savepath, savefn) ,monthlyMean)
        else:
          print("Failed to open netCDF file  " + fullpath)            



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
##Make sure the numpy array slicing/indexing and means used above worked
##Do it the old fashion way... by hand
##unbelievable slow
#n = 10
#means = np.empty([vars.shape[2],vars.shape[3]])
#for x in range(vars.shape[2]):     
#    for y in range(vars.shape[3]):
#        for idepth in range(0,n):
#            means[x,y] +=vars[0,idepth,x,y]

#        means[x,y]= means[x,y]/n

#np.savetxt("Z:\\Projects\\Salish Sea\\NEMO\\Data\\Temp_Mean_10m.txt",means)
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    

