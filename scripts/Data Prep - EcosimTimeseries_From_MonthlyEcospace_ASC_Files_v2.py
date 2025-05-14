
import os
import pandas as pd

NEMO_RUN = 216
inpath = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//forcing//NEMO_prepped_as_ASC_monthly//"
outpath = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//forcing//ECOSIM_in_NEMO_monthly//"
maskpath = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//basemap//"
outfilename_template = "SalishSea1500-RUN{}_{}_MonthlyAvg_Timeseries.csv"
climatology_outfilename_template = "SalishSea1500-RUN{}_{}_MonthlyClimatology.csv"

minDepth = 0.01
asc_suffix_template = "_{}_"

# === FUNCTION DEFINITIONS ===
def getDataFrame(fullPath, NaN):
    return pd.read_table(fullPath, skiprows=6, header=None, sep=r'\s+', na_values=[NaN])

# === LOAD MASK ===
ecospacegrid_f = os.path.join(maskpath, "ecospacedepthgrid.asc")
dfDepth = getDataFrame(ecospacegrid_f, "-9999.00000000")
dfDepthMask = dfDepth < minDepth

# === LOOP THROUGH VARIABLES ===
varnames = [d for d in os.listdir(inpath) if os.path.isdir(os.path.join(inpath, d))]

for varname in varnames:
    print(f"Processing variable: {varname}")

    var_folder = os.path.join(inpath, varname)
    asc_suffix = asc_suffix_template.format(varname)
    lstFiles = sorted([f for f in os.listdir(var_folder) if f.lower().endswith(".asc")])

    dates = []
    means = []

    for ascFile in lstFiles:
        try:
            # Extract date info from filename
            date_part = ascFile.split(asc_suffix)[1].replace(".asc", "")
            year, month = map(int, date_part.split("_"))
            date = pd.Timestamp(year=year, month=month, day=15)
        except Exception as e:
            print(f"Skipping file {ascFile}: Could not parse date - {e}")
            continue

        fullpath = os.path.join(var_folder, ascFile)
        df = getDataFrame(fullpath, "-9999.00000000")
        df = df.mask(dfDepthMask)
        mean_val = df.stack().mean()

        dates.append(date)
        means.append(mean_val)

    # Create and save CSV
    dfOut = pd.DataFrame({varname: means}, index=pd.to_datetime(dates))
    dfOut.index.name = "Date"
    outfile = os.path.join(outpath, outfilename_template.format(NEMO_RUN, varname))
    dfOut.to_csv(outfile)
    print(f"Saved: {outfile}")

    # === Calculate and write monthly climatology ===
    dfOut['Month'] = dfOut.index.month
    climatology = dfOut.groupby('Month')[varname].mean().to_frame()
    climatology.index.name = "Month"

    clim_outfile = os.path.join(outpath, climatology_outfilename_template.format(NEMO_RUN, varname))
    climatology.to_csv(clim_outfile)
    print(f"Saved monthly climatology: {clim_outfile}")


exit()

# revisited this 2024-05-23 but not completed
NEMO_RUN = 216
inpath = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//forcing//NEMO_prepped_as_ASC_monthly"
outpath = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//forcing//ECOSIM_in_NEMO_monthly"
maskpath = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//basemap"
outfilename = "SalishSea1500-RUN{}_{}_Daily_Timeseries.csv"
minDepth = 0.01
varname = 'TempVertMean4m'
# file name example
# SalishSea1500-RUN216_TempVertMean4m_1979_01.asc

#Open a Dataframe
def getDataFrame(fullPath,NaN):
    if os.path.exists(fullPath):
        nas = [NaN]
        df = pd.read_table(fullPath, skiprows=6, header=None, delim_whitespace=True,na_values=nas)
        return df

# Land mask
ecospacegrid_f = os.path.join(maskpath, "ecospacedepthgrid.asc")
df_dep = getDataFrame(ecospacegrid_f,"-9999.00000000")
dfLandMask = df_dep == 0 #mask is where elev =0


InputPath = os.path.join(inpath, varname)
lstFiles = os.listdir(InputPath)
timeindexes=range(len(lstFiles))# number of files

lstCols = [varname];

dfOut = pd.DataFrame(index=timeindexes, columns=lstCols)
dfDepth = getDataFrame(ecospacegrid_f,"-9999.00000000")
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
               
outfile = os.path.join(outpath, outfilename.format(varname))
# dfOut.to_csv(path_or_buf=outfile,sep=',',index=False)   
dfOut.to_csv(path_or_buf=outfile,sep=',')   
print(outfile)
    

 