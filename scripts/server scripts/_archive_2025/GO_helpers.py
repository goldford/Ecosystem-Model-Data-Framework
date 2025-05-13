import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import os
from calendar import monthrange
import math


def is_leap_year(year):
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

def buildSortableString(number, nZeros): 
    newstring = str(number)
    while (len(newstring) < nZeros) :
        tmpstr = "0" + newstring
        newstring = tmpstr
    return newstring   


#def saveASCFile(filename, data, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe,sigdigfmt,ASCheader):
#
#    trimmedData = []
#    i =  bottomleft_row_ewe 
#    while i < upperleft_row_ewe:
#        trimmedData.append(data[i,upperleft_col_ewe:])
#        i += 1
#    trimmedData.reverse()
#
#    np.savetxt(filename,trimmedData,fmt=sigdigfmt, delimiter=" ", comments='',header=ASCheader)

# GO2021-12-22 fmt= has big effect on ASC file size, added argument for this
def saveASCFile(filename, data, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe, sigdigfmt, ASCheader, dfPlumeMask=None, dfLandMask=None):
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
    
    if not (dfPlumeMask is None):
      dfEwEGrid = dfEwEGrid.mask(dfPlumeMask)
      
    if not (dfLandMask is None):
      dfEwEGrid = dfEwEGrid.mask(dfLandMask)

    dfEwEGrid.fillna(0.0,inplace=True)
    
    np.savetxt(filename,dfEwEGrid.to_numpy(),fmt=sigdigfmt, delimiter=" ", comments='',header=ASCheader) # GO 20211222 - change from %0.5f

#Open a Dataframe
def getDataFrame(fullPath,NaN):
    if os.path.exists(fullPath):
        nas = [NaN]
        df = pd.read_table(fullPath, skiprows=6, header=None, delim_whitespace=True,na_values=nas)
        return df
