# Created by G Oldford
# July 24, 2024
# Purpose:
#   - visualise forcings by reading NEMO and RDRS output files prepped for Ecospace (3 day blocks)
#     to help  with calibrating bloom timing in ECOSPACE model
#     'region 2' in Ecospace outputs, aka central SoG corresponding to Suchy et al., 2021, approximately
#
# Source:
#
# Input:
#
# Output:
#
import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta



# paths
path_ecospace_map = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap"
file_ecospace_map = "Ecospace_grid_20210208_rowscols.csv"

path_NEMO_ASCs_root = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/forcing"
# after processing to NC they are too large for github, so moved here
path_NEMO_NCs_root = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/NEMO forcings"
path_RDRS_ASCs_root = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings/Wind_RDRS/Ecospace/"
path_RDRS_NCs_root = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings/Wind_RDRS/Ecospace/"

subpath_NEMO_ASCs_PAR = "/ECOSPACE_in_3day_PAR3_Sal4m_2003-2018_20240527/PAR-VarZ-VarK"
subpath_NEMO_ASCs_vars = "/ECOSPACE_in_3day_vars_2003-2018_20240527/ECOSPACE_in_3day_vars"
subfolder_NEMO_PAR = "/PAR-VarZ-VarK"
subfolder_NEMO_PARxMixing = "/RUN216_PARxMixing"
subfolder_NEMO_temp0to10m = "/vartemp1_C_0-10mAvg"
subfolder_NEMO_mixing = "/varmixing_m"

subfolders = {
    "PAR":         path_NEMO_ASCs_root + "ECOSPACE_in_3day_PAR3_Sal4m_2003-2018_20240527/PAR-VarZ-VarK",
    "PARxMixing":  path_NEMO_ASCs_root + "ECOSPACE_in_3day_PAR3_Sal4m_2003-2018_20240527/RUN216_PARxMixing",
    "MixingZ":     path_NEMO_ASCs_root + "ECOSPACE_in_3day_vars_2003-2018_20240527/ECOSPACE_in_3day_vars/varmixing_m",
    "Wind_Stress_10m_RDRS":     path_RDRS_ASCs_root + "stress",
    "Wind_Speed_10m_RDRS":     path_RDRS_ASCs_root + "speed",
    "Temp_0to10m": "ECOSPACE_in_3day_vars_2003-2018_20240527/ECOSPACE_in_3day_vars/vartemp1_C_0-10mAvg",
    "Salt_0to4m":"ECOSPACE_in_3day_vars_2003-2018_20240527/ECOSPACE_in_3day_vars/varsalt2_PSU_0-4m"
}

ASC_file_fmts = {"PAR": "PAR-VarZ-VarK_{}_{}.asc", # month, doy
                 "PARxMixing": "RUN216_PARxMixing_{}_{}.asc",
                 "MixingZ": "varmixing_m_{}_{}.asc",
                 "Wind_Stress_10m_RDRS": "RDRS_windstress10m_{}_{}.asc",
                 "Wind_Speed_10m_RDRS": "RDRS_windspeed10m_{}_{}.asc",
                 "Temp_0to10m": "vartemp1_C_0-10mAvg_{}_{}.asc",
                 "Salt_0to4m": "varsalt2_PSU_0-4m_{}_{}.asc"
                 }

path_regions_asc = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap"
file_regions_asc = "ecospace_regions_3day.asc"

# ECOSPACE MODEL domain map with lats / lons
# OID,lat,lon,depth,NEMO_col,NEMO_row,EWE_row,EWE_col
file_path = os.path.join(path_ecospace_map, file_ecospace_map)
if os.path.exists(file_path):
    ecospace_map_df = pd.read_csv(file_path)
    metadata = {
        "Column Names": ecospace_map_df.columns.tolist(),
        "Number of Rows": len(ecospace_map_df)
    }
    lats = ecospace_map_df['lat']
    lons = ecospace_map_df['lon']
else:
    print("Error: File not found.")

# ECOSPACE regions file prep
file_path = os.path.join(path_regions_asc, file_regions_asc)
skiprows = 6
with open(file_path) as f:
    dat_regions = np.loadtxt(f, skiprows=skiprows)

rows = dat_regions.shape[0]
cols = dat_regions.shape[1]
print(rows, cols)

process_asc = False
if process_asc:
    # Loop over subfolders and files
    for var, subfolder in subfolders.items():

        ####################################
        # if var != 'Wind_Speed_10m_RDRS':
        #     continue
        ####################################

        timestamps = []
        file_paths = []
        folder_path = os.path.join(path_NEMO_ASCs_root, subfolder)
        print(folder_path)
        for file_name in os.listdir(folder_path):
            if file_name.endswith('.asc'):
                # Extract year and doy from file name
                parts = file_name.split('_')
                year = int(parts[-2])
                doy = int(parts[-1].split('.')[0])

                # Create timestamp
                timestamp = pd.Timestamp(year, 1, 1) + pd.Timedelta(days=doy - 1)
                timestamps.append(timestamp)

                # Store full file path
                file_paths.append(os.path.join(folder_path, file_name))

        # Sort timestamps and corresponding file paths
        timestamps, file_paths = zip(*sorted(zip(timestamps, file_paths)))

        # Create an empty xarray dataset
        ds = xr.Dataset(
            coords={
                'time': list(timestamps),
                'row': range(1, rows + 1),
                'col': range(1, cols + 1)
            },
            attrs={'description': 'dataset of monthly ASC files'}
        )

        # Add lat, lon, depth, EWE_col, EWE_row as data variables
        ds['lat'] = (('row', 'col'), np.full((rows, cols), np.nan))
        ds['lon'] = (('row', 'col'), np.full((rows, cols), np.nan))
        ds['depth'] = (('row', 'col'), np.full((rows, cols), np.nan))
        ds['EWE_col'] = (('row', 'col'), np.full((rows, cols), np.nan))
        ds['EWE_row'] = (('row', 'col'), np.full((rows, cols), np.nan))

        # Populate the dataset with lat, lon, depth, NEMO_col, NEMO_row from the CSV file
        for _, row in ecospace_map_df.iterrows():
            ewe_row = int(row['EWE_row']) - 1
            ewe_col = int(row['EWE_col']) - 1
            ds['lat'][ewe_row, ewe_col] = row['lat']
            ds['lon'][ewe_row, ewe_col] = row['lon']
            ds['depth'][ewe_row, ewe_col] = row['depth']
            ds['EWE_col'][ewe_row, ewe_col] = row['EWE_col']
            ds['EWE_row'][ewe_row, ewe_col] = row['EWE_row']

        ds[var] = (('time', 'row', 'col'), np.full((len(timestamps), rows, cols), np.nan))

        # Read and load data from each ASC file into the dataset
        for t, file_path in enumerate(file_paths):
            parts = file_path.split('_')
            year = int(parts[-2])
            doy = int(parts[-1].split('.')[0])
            with open(file_path) as f:
                data = np.loadtxt(f, skiprows=skiprows)
                ds[var][t-1, :, :] = data

        # Save the dataset to a file
        if (var != "Wind_Stress_10m_RDRS") & (var != "Wind_Speed_10m_RDRS"):
            ds.to_netcdf(os.path.join(path_NEMO_NCs_root, f'{var}.nc'))
            print(f"Saved {var}.nc to {path_NEMO_NCs_root}")
        else:
            ds.to_netcdf(os.path.join(path_RDRS_NCs_root, f'{var}.nc'))
            print(f"Saved {var}.nc to {path_RDRS_NCs_root}")

doy_suchy = [100, 68, 50, 83, 115, 115, 100, 100, 92, 88, 92, 92, 55, 77]

dates_suchy = []
for i, doy in enumerate(doy_suchy):
    year = 2003 + i
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    dates_suchy.append(date)

# lag 1: minus 3 days to account for time it takes for phyto to respond
dates_suchy_lag1 = []
for i, doy in enumerate(doy_suchy):
    year = 2003 + i
    date = datetime(year, 1, 1) + timedelta(days=doy - 1 - 3)
    dates_suchy_lag1.append(date)

print(dates_suchy)
print(dates_suchy_lag1)

# List of variables
vars = ["PAR", "PARxMixing", "MixingZ", "Wind_Stress_10m_RDRS", "Wind_Speed_10m_RDRS", "Temp_0to10m", "Salt_0to4m"]

# Create a figure with subplots
fig, axes = plt.subplots(len(vars), 1, figsize=(12, len(vars) * 3), sharex=False)

# Initialize lists to store the extracted values - for scatter plot investigation
par_values = []
mixing_values = []
wind_stress_values = []
wind_speed_values = []
par_values_lag1 = []
mixing_values_lag1 = []
years = []

# Function to find the closest date in a dataset - match suchy dates to ecospace
def find_closest_date(ds, target_date):
    closest_date = ds['time'].sel(time=target_date, method='nearest')
    return closest_date

alt_Mixing = False # flag to use just the cells that actually vary in terms of MixingZ

# Loop over the variables and create a subplot for each
for i, var in enumerate(vars):
    # Load the dataset
    if (var != 'Wind_Stress_10m_RDRS') & (var != 'Wind_Speed_10m_RDRS'):
        ds = xr.open_dataset(os.path.join(path_NEMO_NCs_root, f'{var}.nc'))
    else:
        ds = xr.open_dataset(os.path.join(path_RDRS_NCs_root, f'{var}.nc'))

    # Compute the average for the region where dat_regions == 2
    if (var == 'MixingZ') & alt_Mixing:
        # mask = ((dat_regions == 2) & ((ds[var] > 10.1)))
        mask = (dat_regions == 2)
        ds[var] = np.log(ds[var]+0.1)
        var_avg = ds[var].where(mask).mean(dim=['row', 'col'])
    else:
        mask = (dat_regions == 2)
        var_avg = ds[var].where(mask).mean(dim=['row', 'col'])


    # find the matching data corresponding to bloom
    if var == 'MixingZ' or var == 'PAR' or var =='Wind_Stress_10m_RDRS' or var =='Wind_Speed_10m_RDRS':

        # Extract values for each date in dates_suchy
        for date in dates_suchy:
            closest_date = find_closest_date(ds, date)
            value = var_avg.sel(time=closest_date).values

            if var == 'PAR':
                par_values.append(value)
            elif var == 'MixingZ':
                mixing_values.append(value)
            elif var == 'Wind_Stress_10m_RDRS':
                wind_stress_values.append(value)
            elif var == 'Wind_Speed_10m_RDRS':
                wind_speed_values.append(value)

        # Lag 1  - Extract values for each date in dates_suchy (using time delta)
        for date in dates_suchy:
            closest_date = find_closest_date(ds, date - pd.Timedelta(days=3))
            value_lag1 = var_avg.sel(time=closest_date).values

            if var == 'PAR':
                par_values_lag1.append(value_lag1)
            elif var == 'MixingZ':
                mixing_values_lag1.append(value_lag1)

    # Plot the time series
    axes[i].plot(ds['time'], var_avg, label=var, linewidth=0.6)
    axes[i].set_ylabel(var)
    axes[i].legend()

    # Add vertical dashed lines at each date in dates_suchy
    for date in dates_suchy:
        axes[i].axvline(date, color='black', linestyle='--', linewidth=1)

    for date in dates_suchy:
        axes[i].axvline(datetime(date.year, 1, 1), color='blue', linestyle='--', linewidth=0.6)

    # Add vertical dashed lines at each date in dates_suchy
    # for date in dates_suchy_lag1:
    #     axes[i].axvline(date, color='grey', linestyle='--', linewidth=0.6)

    # Calculate medians and suchy threshold
    var_avg_df = var_avg.to_dataframe().reset_index()
    var_avg_df['year'] = var_avg_df['time'].dt.year
    yearly_medians = var_avg_df.groupby('year')[var].median()
    yearly_threshold = yearly_medians * 1.05

    for year, threshold in yearly_threshold.items():
        start_of_year = pd.Timestamp(f'{year}-01-01')
        end_of_year = pd.Timestamp(f'{year}-12-31')
        axes[i].hlines(threshold, start_of_year, end_of_year, colors='red', linestyle='-', linewidth=1)


axes[-1].set_xlabel('Time')
plt.suptitle('Time Series of Variables')
if alt_Mixing:
    plt.savefig('..//figs//' + 'forcing_asc_ts_altMixing.png', bbox_inches='tight', dpi=300)
else:
    plt.savefig('..//figs//' + 'forcing_asc_ts.png', bbox_inches='tight', dpi=300)
plt.show()


# Extract values for each date in dates_suchy
for date in dates_suchy:
    years.append(date.year)

#################################################
################## SCATTER PLOTS ################
# 1/ scatter plot - PAR vs Mixing on day of bloom
plt.figure(figsize=(10, 6))
plt.scatter(doy_suchy, par_values, c='blue', marker='o')
# Annotate each point with the year
for i, date in enumerate(dates_suchy):
    plt.annotate(str(date.year) + ", " + str(date.timetuple().tm_yday), (doy_suchy[i], par_values[i]),
                 textcoords="offset points", xytext=(5,5), ha='center', fontsize=8)
plt.xlabel('Bloom Day')
plt.ylabel('PAR')
plt.title('Scatter Plot of Bloom Day vs PAR')
plt.grid(True)
plt.savefig('..//figs//' + 'scatter_PAR_vs_BloomDay.png', bbox_inches='tight', dpi=300)
plt.show()

# 2/ scatter plot - Mixing vs Bloom Day on day of bloom
plt.figure(figsize=(10, 6))
plt.scatter(doy_suchy, mixing_values, c='blue', marker='o')
# Annotate each point with the year
for i, date in enumerate(dates_suchy):
    plt.annotate(str(date.year) + ", " + str(date.timetuple().tm_yday), (doy_suchy[i], mixing_values[i]),
                 textcoords="offset points", xytext=(5,5), ha='center', fontsize=8)
plt.xlabel('Bloom Day')
plt.ylabel('Mixing Z')
plt.title('Scatter Plot of Bloom Day vs Mixing')
plt.grid(True)
plt.savefig('..//figs//' + 'scatter_Mixing_vs_BloomDay.png', bbox_inches='tight', dpi=300)
plt.show()

# 3/ scatter plot - Bloom Day vs Mixing on day of bloom - -3 days lag
doy_suchy_lag = [doy - 3 for doy in doy_suchy]

plt.figure(figsize=(10, 6))
plt.scatter(doy_suchy_lag, mixing_values_lag1, c='blue', marker='o')
# Annotate each point with the year
for i, date in enumerate(dates_suchy):
    plt.annotate(str(date.year) + ", " + str(date.timetuple().tm_yday - 3), (doy_suchy_lag[i], mixing_values_lag1[i]),
                 textcoords="offset points", xytext=(5,5), ha='center', fontsize=8)
plt.xlabel('Bloom Day')
plt.ylabel('Mixing Z')
plt.title('Scatter Plot of Bloom Day vs Mixing - 3 Day Lag')
plt.grid(True)
plt.savefig('..//figs//' + 'scatter_Mixing_vs_BloomDay_3DayLag.png', bbox_inches='tight', dpi=300)
plt.show()

# 4/ scatter plot - Bloom Day vs Mixing on day of bloom - -3 days lag
plt.figure(figsize=(10, 6))
plt.scatter(doy_suchy_lag, par_values_lag1, c='blue', marker='o')
# Annotate each point with the year
for i, date in enumerate(dates_suchy):
    plt.annotate(str(date.year) + ", " + str(date.timetuple().tm_yday - 3), (doy_suchy_lag[i], par_values_lag1[i]),
                 textcoords="offset points", xytext=(5, 5), ha='center', fontsize=8)
plt.xlabel('Bloom Day')
plt.ylabel('PAR')
plt.title('Scatter Plot of Bloom Day vs PAR - 3 Day Lag')
plt.grid(True)
plt.savefig('..//figs//' + 'scatter_PAR_vs_BloomDay_3DayLag.png', bbox_inches='tight', dpi=300)
plt.show()


# 5/ scatter plot - PAR vs Mixing on day of bloom
plt.figure(figsize=(10, 6))
plt.scatter(mixing_values, par_values, c='blue', marker='o')
# Annotate each point with the year
for i, date in enumerate(dates_suchy):
    plt.annotate(str(date.year) + ", " + str(date.timetuple().tm_yday), (mixing_values[i], par_values[i]),
                 textcoords="offset points", xytext=(5,5), ha='center', fontsize=8)
plt.xlabel('Mixing')
plt.ylabel('PAR')
plt.title('Scatter Plot of Mixing vs PAR')
plt.grid(True)
plt.savefig('..//figs//' + 'scatter_PAR_vs_Mixing.png', bbox_inches='tight', dpi=300)
plt.show()


# /6 scatter plot - lag - par vs mixing on a 3 day lag from bloom
plt.figure(figsize=(10, 6))
plt.scatter(par_values_lag1, mixing_values_lag1, c='blue', marker='o')
# Annotate each point with the year
for i, date in enumerate(dates_suchy_lag1):
    plt.annotate(str(date.year) + ", " + str(date.timetuple().tm_yday), (par_values_lag1[i], mixing_values_lag1[i]),
                 textcoords="offset points", xytext=(5,5), ha='center', fontsize=8)
    i += 1
plt.xlabel('PAR Values')
plt.ylabel('Mixing Values')
plt.title('Scatter Plot of PAR vs Mixing Values - 3Day Lag')
plt.grid(True)
plt.savefig('..//figs//' + 'scatter_PAR_vs_Mixing_3DayLag.png', bbox_inches='tight', dpi=300)
plt.show()


# /7 scatter plot - wind stress versus bloom day suchy
plt.figure(figsize=(10, 6))
plt.scatter(doy_suchy, wind_stress_values, c='blue', marker='o')
# Annotate each point with the year
for i, date in enumerate(dates_suchy):
    plt.annotate(str(date.year) + ", " + str(date.timetuple().tm_yday), (doy_suchy[i], wind_stress_values[i]),
                 textcoords="offset points", xytext=(5,5), ha='center', fontsize=8)
plt.xlabel('Bloom Day')
plt.ylabel('Wind Stress')
plt.title('Scatter Plot of Bloom Day vs Wind Stress')
plt.grid(True)
plt.savefig('..//figs//' + 'scatter_WindStress_vs_BloomDay.png', bbox_inches='tight', dpi=300)
plt.show()

###################