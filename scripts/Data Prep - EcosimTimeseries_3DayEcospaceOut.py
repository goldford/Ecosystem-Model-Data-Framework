# Created March 28 2025 by G Oldford
# Purpose: process forcings from Ecospace 3day model for forcings of a 1 yr Ecospace and Ecosim model
# Input: NetCDF created from Ecospace 3day out
# Output: ASC forcings
#
# Notes:
#  -

import os
import pandas as pd
import netCDF4 as nc
import numpy as np
import xarray as xr
from helpers import saveASCFile

in_p = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//"
in_f = "FULLKEY_Scv51_5-PAR_PI_AllPPTemp_Wind_1978-2018.nc"

out_p = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//forcing//"
out_sub_p1 = "ECOSPACE3Day_out_monthlyavg_ECOSIM//"
out_sub_p2 = "ECOSPACE3Day_out_monthlyavg_ECOSPACE//"

# takes a good amount of time to do averaging
# does one field at a time, one nc avg'd file per field
# turn switch off if already done
avg_fields_switch = True

##########################################
########## LOAD ECOSPACE OUT NC ##########
ecospace_nc = os.path.join(in_p, in_f)
ds = xr.open_dataset(ecospace_nc)



########## MASKS ##########
# load mask for land, plume etc.
dep = ds['depth'].values
print(dep)

# monthly average
# List of data variables to process
variables = list(ds.data_vars)


# Process one variable at a time
if avg_fields_switch:
    for var in variables:
        print(f"Processing variable: {var}")

        if "PZ" in var:

            # Output file name
            output_file = "monthly_averages_SCv51_" + var + ".nc"

            # Select one variable at a time
            ds_var = ds[[var]]

            # Compute monthly averages while keeping years separate
            ds_var_monthly = ds_var.resample(time="1M").mean()

            # Save to file (appending each variable to avoid large memory use)
            mode = "w" if var == variables[0] else "a"  # "w" for first variable, "a" (append) for others
            ds_var_monthly.to_netcdf(output_file, mode=mode)

            print(ds_var_monthly)
            # Clean up memory
            del ds_var, ds_var_monthly

            print(f"Finished processing {var}")

# Close dataset
ds.close()




# convert to TS for Ecosim
# convert to TS for ECOSPACE

print("done")
exit()




# revisited this 2024-05-23 but not completed

# OLD CODE
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#Local Functions

#Open a Dataframe
def getDataFrame(fullPath,NaN):
    if os.path.exists(fullPath):
        nas = [NaN]
        df = pd.read_table(fullPath, skiprows=6, header=None, delim_whitespace=True,na_values=nas)
        return df
 
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#InputPath = "C:\\Sync\\PSF\\EwE\\Georgia Strait 2021\\LTL_model\\DATA\\SalishSea1500-RUN202\\NEMO_out_ASC\\MixingTurboZ"
#outputPath = "C:\\Sync\\PSF\\EwE\\Georgia Strait 2021\\LTL_model\\DATA\\SalishSea1500-RUN202\\NEMO_out_ASC\\TimeSeries_Files"

inpath = "C://Users//Greig//Documents//github//Ecosystem-Model-Data-Framework//data//forcing//NEMO216_ECOSPACE_in_daily_vars_20240523//"
outpath = "C://Users//Greig//Documents//github//Ecosystem-Model-Data-Framework//data//forcing//NEMO216_ECOSIM_in_daily_20240523//"


outfilename = "SalishSea1500-RUN{}_{}_Daily_Timeseries.csv"

minDepth = 0.01

DepthFilename = "W:\\Sync\\PSF\\EwE\\Georgia Strait 2021\\LTL_model\\DATA\\NEMO_1.5km_Depth.asc"

lstFiles = os.listdir(InputPath)
timeindexes=range(len(lstFiles))# number of files
varname = 'MixingTurboZ'
lstCols = [varname];

dfOut = pd.DataFrame(index=timeindexes, columns=lstCols)

dfDepth = getDataFrame(DepthFilename,"-9999.00000000")
dfDepthMask = dfDepth < minDepth
iTs = 0      
for ascFile in lstFiles:
    if ascFile.endswith(".asc"): 
        fullpath = os.path.join(InputPath, ascFile)
        df = getDataFrame(fullpath,"-9999.00000000")
        df = df.mask(dfDepthMask)
        mean =  df.stack().mean()
        n= df.stack().count()
        dfOut._set_value(iTs,varname, mean)
        iTs+=1
        #print(n)
               
outfile = os.path.join(outputPath, outfilename.format(varname))
# dfOut.to_csv(path_or_buf=outfile,sep=',',index=False)   
dfOut.to_csv(path_or_buf=outfile,sep=',')   
print(outfile)
    
#print("Done")   
               
 