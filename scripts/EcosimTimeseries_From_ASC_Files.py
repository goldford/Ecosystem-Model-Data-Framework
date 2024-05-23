import os 
import pandas as pd
                 
 
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#Local Functions

#Open a Dataframe
def getDataFrame(fullPath,NaN):
    if os.path.exists(fullPath):
        nas = [NaN]
        df = pd.read_table(fullPath, skiprows=6, header=None, delim_whitespace=True,na_values=nas)
        return df
 
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


InputPath = "W:\\Sync\\PSF\\EwE\\Georgia Strait 2021\\LTL_model\\DATA\\SalishSea1500-RUN202\\NEMO_out_ASC\\MixingTurboZ"
outputPath = "W:\\Sync\\PSF\\EwE\\Georgia Strait 2021\\LTL_model\\DATA\\SalishSea1500-RUN202\\NEMO_out_ASC\\TimeSeries_Files"

outfilename = "SalishSea1500-RUN202_{}_Timeseries.csv"

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
        fullpath = os.path.join(InputPath,ascFile)       
        df = getDataFrame(fullpath,"-9999.00000000")
        df = df.mask(dfDepthMask)
        mean =  df.stack().mean()
        n= df.stack().count()
        dfOut._set_value(iTs,varname, mean)
        iTs+=1
        #print(n)
               
outfile = os.path.join(outputPath,outfilename.format(varname))
# dfOut.to_csv(path_or_buf=outfile,sep=',',index=False)   
dfOut.to_csv(path_or_buf=outfile,sep=',')   
print(outfile)
    
#print("Done")   
               
 