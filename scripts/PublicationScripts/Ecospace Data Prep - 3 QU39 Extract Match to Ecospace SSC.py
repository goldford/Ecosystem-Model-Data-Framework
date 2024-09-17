# Created by G Oldford
# July-Sep 2024
# Purpose: use previously prepped data from del bel belluz and extract from
#          ecospace and ssc from same location
# Source:
#    Del Bel Belluz, J. (2024). Protistan plankton time series from the northern Salish Sea and central coast,
#    British Columbia, Canada—Ocean Biodiversity Information System (a62c37c4-6fbe-4a93-8c87-c67dda36436c)
#    [Dataset]. https://obis.org/dataset/071a38b3-0d30-47ba-85aa-d98313460075
# Input: QU39 phyto data, SalishSeaCast outputs, Ecospace output
#
# Output: table of QU39 obs matched to salishseacast and ecospace outputs
#
#
from helpers import (get_sscast_grd_idx, find_nearest_point, get_sscast_data,
                     find_closest_depth, find_nearest_point_from1D, find_closest_date)
import os
import netCDF4 as nc
import pandas as pd
import numpy as np
import requests
from datetime import datetime

do_ecospace_matching = True
download_ssc = False
match_ssc = True

pathfile_QU39_in = r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\Phyto Concent del Bel Belluz 2024 2016 - 2019\MODIFIED\QU39_joined.csv'
pathfile_QU39_out_p = 'C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//28. Phytoplankton//Phyto Concent del Bel Belluz 2024 2016 - 2019//MODIFIED//'
matched_QU39_out_f = 'QU39_joined_matchtoEcospace'
path_SSC = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/SalishSeaCast_BioFields_2008-2018/ORIGINAL"
path_SSC_dwnld = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/SalishSeaCast_BioFields_2008-2018/ORIGINAL/matched_QU39_dwnld/"
path_ecospace_out = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/ECOSPACE_OUT"
path_ecospace_map = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap"

# file_Ecospace_out = "Scv7-PARMixingNut90Temp_2003-2018.nc" # better fits?
# file_Ecospace_out = "Scv39-PARenv_PI_Temp_Wind_2003-2018.nc" # better fits?
file_ecospace_out = "Scv51_4_2-PAR_PI_AllPPTemp_Wind_2000-2018.nc" # same as key run, just shorter
ecospace_code = "SC51_4_2"

file_ecospace_map = "Ecospace_grid_20210208_rowscols.csv"
file_SSC_bathy = "ubcSSnBathymetryV21-08_a29d_efc9_4047.nc"
file_SSC_meshm_3D = 'ubcSSn3DMeshMaskV17-02.nc'
file_SSC_meshm_2D = 'ubcSSn2DMeshMaskV17-02.nc'

# get lat lon QU39
QU39_lat = 50.0307 # hardcoding
QU39_lon = -125.0992

#####################################################
################ LOAD OBSERVATIONS ##################
QU39_df = pd.read_csv(pathfile_QU39_in, sep=',', engine='python')
# QU39_df = pd.read_csv(pathfile_QU39_out, sep=',', engine='python') # temp switch
# print(QU39_df)
print(QU39_df.head()) #7387 rows, 53 cols
print(QU39_df.columns)

skip = False
if not skip:
    # Combine the two columns into a single datetime string
    # eg minimumDepthInMeters, maximumDepthInMeters, season, eventTime, year, month, day, 'anomaly_fr_mean', 'anomaly_fr_median'
    # kingdom', 'phylum', 'class', 'order', 'family', 'genus'
    QU39_df['PP1-DIA'] = -999.9; QU39_df['PP2-NAN'] = -999.9
    QU39_df['PP3-PIC'] = -999.9; QU39_df['PZ1-CIL'] = -999.9
    QU39_df['PZ2-DIN'] = -999.9
    QU39_df['ecospace_closest_lat'] = 0.0; QU39_df['ecospace_closest_lon'] = 0.0
    QU39_df['ecospace_closest_gridY'] = 0.0; QU39_df['ecospace_closest_gridX'] = 0.0
    QU39_df['ecospace_closest_dep'] = 0.0;
    QU39_df['closest_ecospace_time'] = pd.to_datetime('1900-04-23 23:42') # initialising

# fix time date issues for conversions and matching
QU39_df['eventDate'] = QU39_df['eventDate'].astype(str)
QU39_df['eventTime'] = QU39_df['eventTime'].astype(str)
QU39_df['DateTime'] = QU39_df['eventDate'] + ' ' + QU39_df['eventTime']
QU39_df['DateTime'] = pd.to_datetime(QU39_df['DateTime'], format='%Y-%m-%d %H:%M:%SZ')
QU39_df['ssc-DIA'] = -999.9; QU39_df['ssc-FLA'] = -999.9; QU39_df['ssc-CIL'] = -999.9
QU39_df['ssc_closest_lat'] = 0.0;
QU39_df['ssc_closest_lon'] = 0.0
QU39_df['ssc_closest_gridY'] = 0.0;
QU39_df['ssc_closest_gridX'] = 0.0
QU39_df['closest_ssc_time'] = pd.to_datetime('1900-04-23 23:42') #
QU39_df['closest_ssc_dep'] = 0.0;

#####################################################
############# GET ECOSPACE MATCHING  ################
# get Ecospace grid info
ecospace_coords_df = pd.read_csv(os.path.join(path_ecospace_map, file_ecospace_map), sep=',', engine='python')
ecospace_lats = ecospace_coords_df['lat']
ecospace_lons = ecospace_coords_df['lon']
ecospace_rows = ecospace_coords_df['EWE_row']
ecospace_cols = ecospace_coords_df['EWE_col']
ecospace_deps = ecospace_coords_df['depth']
#0.01 deg lat = ~1.1 km
dx = 0.01
mask_ecospace = ecospace_deps != 0
mod_lats_nearest = []; mod_lons_nearest = []
mod_rows_nearest = []; mod_cols_nearest = []

# Get ecospace lat lon (closest) for each QU39 point
if do_ecospace_matching:
    with nc.Dataset(os.path.join(path_ecospace_out, file_ecospace_out), 'r') as dataset:

        time_var = dataset.variables['time']
        time_units = time_var.units
        time_calendar = time_var.calendar if 'calendar' in time_var.ncattrs() else 'standard'
        ecospace_times = nc.num2date(time_var[:], units=time_units, calendar=time_calendar)

        for index, row in QU39_df.iterrows():
            # pos = (QU39_lat['Lat'], row['Lon'])

            # ecospace lat lon (closest)
            dx = 0.01
            p = find_nearest_point_from1D(QU39_lon, QU39_lat, ecospace_lons, ecospace_lats, mask_ecospace, dx)
            nearest_idx = p[0]

            if not nearest_idx is np.nan:
                QU39_df.at[index, 'ecospace_closest_lat'] = ecospace_lats[nearest_idx]
                QU39_df.at[index, 'ecospace_closest_lon'] = ecospace_lons[nearest_idx]
                QU39_df.at[index, 'ecospace_closest_gridY'] = ecospace_rows[nearest_idx]
                QU39_df.at[index, 'ecospace_closest_gridX'] = ecospace_cols[nearest_idx]
                row_idx = ecospace_rows[nearest_idx]
                col_idx = ecospace_cols[nearest_idx]
            else:
                QU39_df.at[index, 'ecospace_closest_lat'] = -999
                QU39_df.at[index, 'ecospace_closest_lon'] = -999
                QU39_df.at[index, 'ecospace_closest_gridY'] = -999
                QU39_df.at[index, 'ecospace_closest_gridX'] = -999
                row_idx = -999
                col_idx = -999

            # ecospace times (closest)
            obs_time = row['DateTime']

            if obs_time.year <= 2018:

                closest_time = find_closest_date(ecospace_times, obs_time)
                QU39_df.at[index, 'closest_ecospace_time'] = closest_time

                # grab ecospace model outputs
                time_idx = np.where(ecospace_times == closest_time)[0][0]

                # Retrieve PP1-DIA value
                if row_idx != -999 and col_idx != -999:

                    pp1_dia_value = dataset.variables['PP1-DIA'][time_idx, row_idx, col_idx]
                    pp2_nan_value = dataset.variables['PP2-NAN'][time_idx, row_idx, col_idx]
                    pp3_pic_value = dataset.variables['PP3-PIC'][time_idx, row_idx, col_idx]
                    pz1_cil_value = dataset.variables['PZ1-CIL'][time_idx, row_idx, col_idx]
                    pz2_din_value = dataset.variables['PZ2-DIN'][time_idx, row_idx, col_idx]
                    QU39_df.at[index, 'PP1-DIA'] = pp1_dia_value
                    QU39_df.at[index, 'PP2-NAN'] = pp2_nan_value
                    QU39_df.at[index, 'PP3-PIC'] = pp3_pic_value
                    QU39_df.at[index, 'PZ1-CIL'] = pz1_cil_value
                    QU39_df.at[index, 'PZ2-DIN'] = pz2_din_value
    QU39_df.to_csv(os.path.join(pathfile_QU39_out_p, matched_QU39_out_f + "_" + ecospace_code + ".csv"), index=False)

else:
    pd.read_csv(os.path.join(pathfile_QU39_out_p, matched_QU39_out_f + "_" + ecospace_code + ".csv"), sep=',', engine='python')
    # necessary?
    QU39_df['DateTime'] = pd.to_datetime(QU39_df['DateTime'], format='%Y-%m-%d %H:%M:%SZ')

#####################################################
########### GET SSC MATCHING GRID CELL ##############
# get 2D grid idx for SalishSeaCast from the bathy file
ssc_grid_metadata, ssc_lats, ssc_lons, ssc_gridY, ssc_gridX, ssc_grid_bathy = get_sscast_grd_idx(os.path.join(path_SSC, file_SSC_bathy))
mask_ssc = ssc_grid_bathy != 0

# find the closest SSCast grid cell, p
# 0.01 deg lat = ~1.1 km so 1/2 to get max dist from cell corner to centre for 500 m cell wid
dx = 0.01 * 0.5
# returns the grid indices (returns p [col, row, dist])
p = find_nearest_point(QU39_lon, QU39_lat, ssc_lons, ssc_lats, mask_ssc, dx)
ssc_col = p[0]
ssc_row = p[1]
ssc_lat = ssc_lats[p[1], p[0]]
ssc_lon = ssc_lons[p[1], p[0]]

print(ssc_lat, ssc_lon)
print(ssc_row, ssc_col)

##### GET SSC MATCHING GRID DEPTH #######
# get gdept_0 for grid cell from the 3D Mesh mask
with nc.Dataset(os.path.join(path_SSC, file_SSC_meshm_3D)) as mesh:
    tmask = mesh.variables['tmask'][:]  # 0's where depth lev exceeds
    e3t0 = mesh.variables['e3t_0'][:]  # 'widths' or 'weights' for each depth lev
    gdept_0 = mesh.variables['gdept_0'][:]  # depth under t grid points (1, 40, 898, 398)

ssc_possible_depths = gdept_0[:,:,p[1],p[0]]

###### DOWNLOAD SSC FILES
# just get unique dates and loop for downloading

if download_ssc:
    unique_df = QU39_df.drop_duplicates(subset=['DateTime', 'minimumDepthInMeters', 'maximumDepthInMeters'])

    # loop over obs, get time and depth, download corresponding ssc file
    for index, row in unique_df.iterrows():

        obs_depth = (row['minimumDepthInMeters'] + row['maximumDepthInMeters']) / 2
        ssc_depth = find_closest_depth(obs_depth, ssc_possible_depths[0][:])
        obs_time = row['DateTime']

        ssc_filename_dt = obs_time.strftime('%Y-%m-%d')
        ssc_filename = path_SSC_dwnld + f"SalishSeaCast_biology_v1905_{ssc_filename_dt}-{obs_depth}.nc"

        print("Downloading")
        # Generate the URL
        formatted_dt = obs_time.strftime('%Y-%m-%dT%H:00')
        print(formatted_dt)

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
            with open(ssc_filename, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded: {ssc_filename}")
        else:
            print(f"Failed to download data for {formatted_dt}. HTTP Status code: {response.status_code}")

    #open file
    ssc_meta, ssc_times, ssc_ciliates, ssc_flagellates, ssc_diatoms = get_sscast_data(ssc_filename)

    if 'Error' in ssc_meta:
        print(f"Error encountered: {ssc_meta['Error']}")
    else:
        print('Successfully found file' + ssc_filename)


if match_ssc:
    for index, row in QU39_df.iterrows():
        obs_depth = (row['minimumDepthInMeters'] + row['maximumDepthInMeters']) / 2
        ssc_filename_dt = row['DateTime'].strftime('%Y-%m-%d')
        ssc_filename = path_SSC_dwnld + f"SalishSeaCast_biology_v1905_{ssc_filename_dt}-{obs_depth}.nc"
        ssc_meta, ssc_times, ssc_ciliates, ssc_flagellates, ssc_diatoms = get_sscast_data(ssc_filename)
        if 'Error' in ssc_meta:
            print(f"Error encountered: {ssc_meta['Error']}")
        else:
            print('Successfully found file' + ssc_filename)

        QU39_df.at[index, 'ssc-DIA'] = ssc_diatoms[0, 0, 0, 0]
        QU39_df.at[index, 'ssc-FLA'] = ssc_flagellates[0, 0, 0, 0]
        QU39_df.at[index, 'ssc-CIL'] = ssc_ciliates[0, 0, 0, 0]
        QU39_df.at[index, 'ssc_closest_gridY'] = ssc_row
        QU39_df.at[index, 'ssc_closest_gridX'] = ssc_col
        QU39_df.at[index, 'ssc_closest_lat'] = ssc_lat
        QU39_df.at[index, 'ssc_closest_lon'] = ssc_lon
        date_converted = datetime(ssc_times[0].year, ssc_times[0].month,
                                  ssc_times[0].day, ssc_times[0].hour,
                                  ssc_times[0].minute, ssc_times[0].second)
        QU39_df.at[index, 'closest_ssc_time'] = date_converted #shape?

    QU39_df.to_csv(os.path.join(pathfile_QU39_out_p, matched_QU39_out_f + "_SSC_" + ecospace_code + ".csv"), index=False)
    print("wrote " + matched_QU39_out_f + "_SSC_" + ecospace_code + ".csv")
else:
    pd.read_csv(os.path.join(pathfile_QU39_out_p, matched_QU39_out_f + "_SSC_" + ecospace_code + ".csv"), sep=',', engine='python')
    # necessary?
    QU39_df['DateTime'] = pd.to_datetime(QU39_df['DateTime'], format='%Y-%m-%d %H:%M:%SZ')



