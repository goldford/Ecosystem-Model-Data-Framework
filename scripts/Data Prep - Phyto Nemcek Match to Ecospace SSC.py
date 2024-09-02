# Created: July 11, 2024
# Author: G Oldford
# Purpose: Link Ecospace and SSCast outputs to Nemcek et al 2023 phyto obs
#
# Source of Data:
#  - SalishSeaCast phyto data downloaded from url below using another script
#    https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1moV19-05.html
# - Nemcek, N., Hennekes, M., Sastri, A., & Perry, R. I. (2023). Seasonal and spatial dynamics of the phytoplankton
#      community in the Salish Sea, 2015–2019. Progress in Oceanography, 217, 103108.
#      https://doi.org/10.1016/j.pocean.2023.103108
# File name examples:
#  - SScast: SalishSeaCast_biology_2008.nc
#  - Nemcek et al: Nemcek_Supp_Data.csv
#  - Scv7-PARMixingNut90Temp_2003-2018.nc
# Metadata:
#  Years Covered:
#    SSCast: 2008 - 2018
#    Nemcek et al. 2023: 2015 - 2019
#    Ecospace - 2003 - 2018 (limited to years with adequate bloom timing data analysed elsewhere)
#
# Notes:
#   -

import os
import pandas as pd
import netCDF4 as nc
from helpers import load_yaml, get_sscast_grd_idx
import numpy as np
import os
import requests
from matplotlib.path import Path as Path
# match points to do (blind) using both model and SSC to compare to obs
from helpers import find_nearest_point, find_nearest_point_from1D
from helpers import find_closest_date, find_closest_depth, read_sdomains, get_sscast_data

# Paths to the files
path_SSC = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/SalishSeaCast_BioFields_2008-2018/ORIGINAL"
path_Ecospace = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/ECOSPACE_OUT"
path_Nemcek = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/Phytoplankton Salish Sea Nemcek2023 2015-2019/MODIFIED"
path_evalout = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
path_ecospace_map = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap"

file_SSC_mo = "SalishSeaCast_biology_2008.nc"
# https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSnBathymetryV21-08.nc?bathymetry[(0.0):1:(897.0)][(0.0):1:(397.0)],latitude[(0.0):1:(897.0)][(0.0):1:(397.0)],longitude[(0.0):1:(897.0)][(0.0):1:(397.0)]
file_SSC_grd = "ubcSSnBathymetryV21-08_a29d_efc9_4047.nc"
file_Ecospace = "Scv7-PARMixingNut90Temp_2003-2018.nc"
file_Nemcek = "Nemcek_Supp_Data.csv"
file_Nemcek_matched = "Nemcek_matched_to_model_out.csv"
file_ecospace_map = "Ecospace_grid_20210208_rowscols.csv"

# Function to get metadata from a CSV file
# Function to get metadata from a CSV file
def get_csv_metadata(file_path):
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        metadata = {
            "Column Names": df.columns.tolist(),
            "Number of Rows": len(df)
        }
        return metadata, df
    else:
        return {"Error": "File not found."}, None


def get_ecospace_times(file_path):
    if os.path.exists(file_path):
        try:
            with nc.Dataset(file_path, 'r') as dataset:
                metadata = {
                    "Dimensions": {dim: len(dataset.dimensions[dim]) for dim in dataset.dimensions},
                    "Variables": {var: dataset.variables[var].dimensions for var in dataset.variables}
                }
                time_var = dataset.variables['time']
                time_units = time_var.units
                time_calendar = time_var.calendar if 'calendar' in time_var.ncattrs() else 'standard'
                times = nc.num2date(time_var[:], units=time_units, calendar=time_calendar)

            return metadata, times
        except OSError as e:
            return {"Error": str(e)}, None, None
    else:
        return {"Error": "File not found."}, None, None


# Paths to files
ssc_file_path = os.path.join(path_SSC, file_SSC_mo)
ecospace_file_path = os.path.join(path_Ecospace, file_Ecospace)
nemcek_file_path = os.path.join(path_Nemcek, file_Nemcek)

# Get metadata from files
ssc_metadata, ssc_times, _, _, _ = get_sscast_data(ssc_file_path)
ecospace_metadata, ecospace_times = get_ecospace_times(ecospace_file_path)
nemcek_metadata, nemcek_df = get_csv_metadata(nemcek_file_path)

# Print the metadata
print("SalishSeaCast NetCDF Metadata:")
print(ssc_metadata)
print("\nEcospace NetCDF Metadata:")
print(ecospace_metadata)
print("\nNemcek CSV Metadata:")
print(nemcek_metadata)

# Inspect time
print("\nSalishSeaCast Times (first 10 entries):")
print(ssc_times[:10])

print("\nEcospace Times (first 10 entries):")
print(ecospace_times[:10])

if nemcek_df is not None:
    print("\nNemcek 'Date.Time' Column Sample:")
    print(nemcek_df['Date.Time'].head())
    print("\nNemcek 'Date.Time' Column Data Type:")
    print(nemcek_df['Date.Time'].dtype)

##############################################################
####### FIND obs 'province' AND MATCH TO ECOSPACE OUT ########
##############################################################
# domain_p = "C:/Users/Greig/Documents/github/HOTTSea_v1/serverside/pypkg/config"
domain_p = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
# domain_f = "CTD_analysis_domain_config_template-SalishSea.yml"
domain_f = "analysis_domains_jarnikova.yml"
domain_fullp = os.path.join(domain_p, domain_f)

sdomains = read_sdomains(domain_fullp)

# get subdomain for each point
# for index, row in nemcek_df.iterrows():
#     pos = row['Lat'], row['Lon']
#     sd_obs = ""
#
#     for sd in sdomains.keys():
#         if Path(sdomains[sd]).contains_point(pos):
#             sd_obs = sd

# Add a new column for subdomain
nemcek_df['sdomain'] = ""

# get ECOSPACE MODEL domain map with lats / lons
ecospace_coords_md, ecospace_coords_df = get_csv_metadata(os.path.join(path_ecospace_map, file_ecospace_map))

print("Ecospace coords metadata: ")
print(ecospace_coords_md)
model_lats = ecospace_coords_df['lat']
model_lons = ecospace_coords_df['lon']
model_rows = ecospace_coords_df['EWE_row']
model_cols = ecospace_coords_df['EWE_col']
model_deps = ecospace_coords_df['depth']
#0.01 deg lat = ~1.1 km
dx = 0.01
mask_ecospace = model_deps != 0
mod_lats_nearest = []; mod_lons_nearest = []
mod_rows_nearest = []; mod_cols_nearest = []
nemcek_df['ecospace_closest_lat'] = 0.0
nemcek_df['ecospace_closest_lon'] = 0.0

# Convert 'Date.Time' to datetime objects
nemcek_df['Date.Time'] = pd.to_datetime(nemcek_df['Date.Time'])
nemcek_df['closest_ecospace_time'] = pd.to_datetime('1900-04-23 23:42') # just initialising

nemcek_df['PP1-DIA'] = 0.0
nemcek_df['PP2-NAN'] = 0.0
nemcek_df['PP3-PIC'] = 0.0
nemcek_df['PZ1-CIL'] = 0.0
nemcek_df['PZ2-DIN'] = 0.0

nemcek_df['ssc-DIA'] = 0.0
nemcek_df['ssc-FLA'] = 0.0
nemcek_df['ssc-CIL'] = 0.0

# open sample file to get grid
ssc_grid_metadata, ssc_lats, ssc_lons, ssc_gridY, ssc_gridX, ssc_grid_bathy = get_sscast_grd_idx(os.path.join(path_SSC, file_SSC_grd))
print("SSCast grid info: ")
print(ssc_grid_metadata)
mask_ssc = ssc_grid_bathy != 0
# for index, row in nemcek_df.iterrows():

# Get subdomain and ecospace lat lon (closest) for each point
with nc.Dataset(ecospace_file_path, 'r') as dataset:
    for index, row in nemcek_df.iterrows():
        # subdomain
        pos = (row['Lat'], row['Lon'])
        sd_obs = ""
        for sd in sdomains.keys():
            if Path(sdomains[sd]).contains_point(pos):
                sd_obs = sd
                break
        nemcek_df.at[index, 'sdomain'] = sd_obs

        # ecospace lat lon (closest)
        lat = row['Lat']
        lon = row['Lon']
        dx = 0.01
        p = find_nearest_point_from1D(lon, lat, model_lons, model_lats, mask_ecospace, dx)
        nearest_idx = p[0]
        if not nearest_idx is np.nan:
            nemcek_df.at[index, 'ecospace_closest_lat'] = model_lats[nearest_idx]
            nemcek_df.at[index, 'ecospace_closest_lon'] = model_lons[nearest_idx]
            row_idx = model_rows[nearest_idx]
            col_idx = model_cols[nearest_idx]
        else:
            nemcek_df.at[index, 'ecospace_closest_lat'] = -999
            nemcek_df.at[index, 'ecospace_closest_lon'] = -999
            row_idx = -999
            col_idx = -999

        # ecospace times (closest)
        obs_time = row['Date.Time']
        obs_depth = row['Pressure']
        closest_time = find_closest_date(ecospace_times, obs_time)
        nemcek_df.at[index, 'closest_ecospace_time'] = closest_time

        # grab ecospace model outputs
        time_idx = np.where(ecospace_times == closest_time)[0][0]

        # Retrieve PP1-DIA value
        if row_idx != -999 and col_idx != -999:

            pp1_dia_value = dataset.variables['PP1-DIA'][time_idx, row_idx, col_idx]
            pp2_nan_value = dataset.variables['PP2-NAN'][time_idx, row_idx, col_idx]
            pp3_pic_value = dataset.variables['PP3-PIC'][time_idx, row_idx, col_idx]
            pz1_cil_value = dataset.variables['PZ1-CIL'][time_idx, row_idx, col_idx]
            pz2_din_value = dataset.variables['PZ2-DIN'][time_idx, row_idx, col_idx]
            nemcek_df.at[index, 'PP1-DIA'] = pp1_dia_value
            nemcek_df.at[index, 'PP2-NAN'] = pp2_nan_value
            nemcek_df.at[index, 'PP3-PIC'] = pp3_pic_value
            nemcek_df.at[index, 'PZ1-CIL'] = pz1_cil_value
            nemcek_df.at[index, 'PZ2-DIN'] = pz2_din_value

        # find the closest SSCast grid cell
        # 0.01 deg lat = ~1.1 km so 1/2 to get max dist from cell corner to centre for 500 m cell wid
        dx = 0.01 * 0.5
        # returns the grid indices
        p = find_nearest_point(lon, lat, ssc_lons, ssc_lats, mask_ssc, dx)
        ssc_possible_depths = [0.5000003, 1.5000031, 2.5000114, 3.5000305]
        ssc_depth = find_closest_depth(obs_depth, ssc_possible_depths)
        # download corresponding SSCast files for each obs
        filename_dt = obs_time.strftime('%Y-%m-%d')
        filename = f"SalishSeaCast_biology_v1905_{filename_dt}.nc"

        download_ssc = False
        if download_ssc:
            print("Downloading")
            # Generate the URL
            ssc_col = p[0]
            ssc_row = p[1]
            formatted_dt = obs_time.strftime('%Y-%m-%dT%H:00')
            print(formatted_dt)

            # Path to the NetCDF file
            # file_path = 'C:/temp/SalishSeaCast_biology_v1905_2017-04-09T23.nc'
            #
            #
            # # Function to inspect the contents of the NetCDF file
            # def inspect_nc_file(file_path):
            #     if os.path.exists(file_path):
            #         try:
            #             with nc.Dataset(file_path, 'r') as dataset:
            #                 # Print global attributes
            #                 print("Global Attributes:")
            #                 for attr in dataset.ncattrs():
            #                     print(f"{attr}: {dataset.getncattr(attr)}")
            #
            #                 # Print dimensions
            #                 print("\nDimensions:")
            #                 for dim in dataset.dimensions.values():
            #                     print(dim)
            #
            #                 # Print variables
            #                 print("\nVariables:")
            #                 for var in dataset.variables.values():
            #                     print(var)
            #
            #                 # Check if there's any data
            #                 data_present = False
            #                 for var_name, var in dataset.variables.items():
            #                     if var.size > 0:
            #                         data_present = True
            #                         break
            #
            #                 if data_present:
            #                     print("\nData is present in the file.")
            #                 else:
            #                     print("\nNo data is present in the file.")
            #
            #         except OSError as e:
            #             print(f"Error opening file: {e}")
            #     else:
            #         print("File not found.")
            #
            # # Inspect the NetCDF file
            # inspect_nc_file(file_path)
            #
            # exit()
            griddap_url = ( # hourly fields from v19-05 (v21-11 has no ciliates)
                    "https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1hV19-05.nc?"                    
                    "ciliates[(" + formatted_dt + ":00Z):1:(" + formatted_dt + ":00Z)][(" + str(ssc_depth) + "):1:(" + str(ssc_depth) + ")][(" +
                    str(ssc_row) + "):1:(" + str(ssc_row) + ")][(" + str(ssc_col) + "):1:(" + str(ssc_col) + ")],"
                    "diatoms[(" + formatted_dt + ":00Z):1:(" + formatted_dt + ":00Z)][(" + str(ssc_depth) + "):1:(" + str(ssc_depth) + ")][(" +
                    str(ssc_row) + "):1:(" + str(ssc_row) + ")][(" + str(ssc_col) + "):1:(" + str(ssc_col) + ")],"
                    "flagellates[(" + formatted_dt + ":00Z):1:(" + formatted_dt + ":00Z)][(" + str(ssc_depth) + "):1:(" + str(ssc_depth) + ")][(" +
                    str(ssc_row) + "):1:(" + str(ssc_row) + ")][(" + str(ssc_col) + "):1:(" + str(ssc_col) + ")]"
            )
            print(griddap_url)

            # Download the file
            response = requests.get(griddap_url)
            if response.status_code == 200:
                with open(filename, 'wb') as f:
                    f.write(response.content)
                print(f"Downloaded: {filename}")
            else:
                print(f"Failed to download data for {formatted_dt}. HTTP Status code: {response.status_code}")

        #open file
        path_ssc_da = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/SalishSeaCast_BioFields_2008-2018/ORIGINAL"
        ssc_hrly_f = os.path.join(path_ssc_da, filename)
        ssc_meta, ssc_times, ssc_ciliates, ssc_flagellates, ssc_diatoms = get_sscast_data(ssc_hrly_f)

        if 'Error' in ssc_meta:
            print(f"Error encountered: {ssc_meta['Error']}")
        else:
            nemcek_df.at[index, 'ssc-DIA'] = ssc_diatoms[0, 0, 0, 0]
            nemcek_df.at[index, 'ssc-FLA'] = ssc_flagellates[0, 0, 0, 0]
            nemcek_df.at[index, 'ssc-CIL'] = ssc_ciliates[0, 0, 0, 0]

# Write the DataFrame to a CSV file
output_file_path = os.path.join(path_evalout, file_Nemcek_matched)
nemcek_df.to_csv(output_file_path, index=False)

# Verify the updated dataframe
print(nemcek_df.head())
print(len(nemcek_df))
print("SGS Obs Count")
print(nemcek_df[(nemcek_df['sdomain'] == 'SGS') & (nemcek_df['ecospace_closest_lat'] != -999) & (nemcek_df['ecospace_closest_lat'] != 0)].shape[0])
print("SGN Obs Count")
print(nemcek_df[(nemcek_df['sdomain'] == 'SGN') & (nemcek_df['ecospace_closest_lat'] != -999) & (nemcek_df['ecospace_closest_lat'] != 0)].shape[0])
print("SGI Obs Count")
print(nemcek_df[(nemcek_df['sdomain'] == 'SGI') & (nemcek_df['ecospace_closest_lat'] != -999) & (nemcek_df['ecospace_closest_lat'] != 0)].shape[0])
print(nemcek_df['PP1-DIA'])
print(nemcek_df[(nemcek_df['PP1-DIA'] != 0.0) & (nemcek_df['ecospace_closest_lat'] != -999) & (nemcek_df['ecospace_closest_lat'] != 0)].shape[0])

# summary statistics, histograms by month by region
