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
lstVars = ['votemper', 'vosaline','vovecrtz','vozocrtx','vomecrty'];

#The key in the dictionary need to match the variable names in lstVars
dctFileTemplates = {'votemper':'SalishSea1500-RUN200_1d_grid_T_{}{}{}-{}{}{}.nc','vosaline':'SalishSea1500-RUN200_1d_grid_T_{}{}{}-{}{}{}.nc','vovecrtz':'SalishSea1500-RUN200_1d_grid_W_{}{}{}-{}{}{}.nc'
                   ,'vozocrtx':'SalishSea1500-RUN200_1d_grid_U_{}{}{}-{}{}{}.nc','vomecrty':'SalishSea1500-RUN200_1d_grid_V_{}{}{}-{}{}{}.nc' }
#Unit Conversions
#'vovecrtz':1000 upwelling/downwelling (vertical velocity) is in meters per second but his makes it really hard parameterize the EwE UI
#so arbitrarily scale them bigger 
#This doesn't matter because we create a response function on the values so that can look like anything 
dctUnitConvert = {'votemper':1,'vosaline':1,'vovecrtz':1000,'vozocrtx':100,'vomecrty':100}


startYear = 1979
endyear = startYear + 1

#Number of depth strata to use
#Depth indexes are zero based
#So 10 will be the 9 index
nDepthIndexes = 10 #= 9.001736   10.003407
#nDepthIndexes = 23 #= 31.101034   39.11886
#nDepthIndexes = 26 #= 86.96747   109.73707 

#array offsets for clipping EwE grid out of NEMO grid
#Calculated by Greig in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39

rootDir = "Z:\\Projects\\Salish Sea\\NEMO\\Data\\"
inputsSubDir = 'Daily-1979\\RUN200-1979-daily'
saveTemplate = "SalishSea1500-RUN200_{}_{}_{}.asc"

ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"

for var in lstVars:
    for iyear in range(startYear, endyear):
        for imon in range(1,13):
            ndays = monthrange(iyear,imon)
            #means = np.empty([vars.shape[2],vars.shape[3]])
            #dailymeans = np.empty([ndays[1],vars.shape[2],vars.shape[3]])
            bnewMonth = True 
            for iday in  range(1,ndays[1]+1):
                m = buildSortableString(imon,2)
                d = buildSortableString(iday,2)
                inputTemplate = dctFileTemplates[var]
                filename = inputTemplate.format(startYear,m,d,startYear,m,d)
                fullpath = os.path.join(rootDir, inputsSubDir,filename)

                print(fullpath)

                if(os.path.isfile(fullpath)):
                    ds = nc.Dataset(fullpath)
                    #print(ds)
                    vars = ds.variables[var]

                    #For looking at the depth/layers
                    #print(ds.variables["deptht"][:])
                    #print(ds.variables["deptht_bounds"][:])
                    #print(vars[:,0:nDepthIndexes,:])

                    if(bnewMonth):
                        dailymeans = np.empty([ndays[1],vars.shape[2],vars.shape[3]])
                        bnewMonth = False
                    #Slice the array on the depth indexes 1-nDepths
                    #Then get the mean of the X and Y across the nDepths (10 for debugging) depth indexes sliced out
                    
                    mean = vars[:,0:nDepthIndexes,:].mean(axis=1)                     

                    #put the daily means into an array of days for the month
                    dailymeans[iday-1] = mean[0]   
                else:
                   print("Failed to open netCDF file  " + fullpath)
        
            #Get the unit conversion multiplier
            UnitConvert = dctUnitConvert[var]
            #get the mean across the days for this month
            monthlyMean = dailymeans.mean(axis=0)
            monthlyMean = monthlyMean * UnitConvert

            #Save the monthly means to an .asc file
            savepath = os.path.join(rootDir, var) 
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)
            
            savefn = saveTemplate.format(var,iyear,m)
            saveASCFile(os.path.join(savepath, savefn) ,monthlyMean)


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

    

