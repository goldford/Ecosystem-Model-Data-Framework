'''
============================================================================================
Script:    Ecospace Data Prep - Match Zooplankton Obs to Ecospace Outputs
Author:    G. Oldford (based on prior Nemcek match script)
Created:   May 2025

Purpose:
   Match zooplankton observations from survey CSV to Ecospace model outputs

Inputs:
   - Ecospace output NetCDF file (e.g., Scv88_1-All_Groups_20250506_1978-2018.nc)
   - Zooplankton observation CSV (e.g., Zooplankton_B_C_gm2_EWEMODELGRP_Wide.csv)
   - Ecospace grid cell mapping (Ecospace_grid_20210208_rowscols.csv)

Outputs:
   - CSV file combining zooplankton observations with matched Ecospace model values
============================================================================================
'''

import os
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime
from helpers import (
    find_nearest_point_from1D,
    find_closest_date
)

# -------------------------------
# Configuration
# -------------------------------
file_Ecospace = "Scv114_2-All_Groups_20250602_1978-2018.nc"
ecospace_code = 'Scv114_2'
file_Zoop = "Zooplankton_B_C_gm2_EWEMODELGRP_Wide.csv"
file_ecospace_map = "Ecospace_grid_20210208_rowscols.csv"

path_Ecospace = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/ECOSPACE_OUT"
path_Zoop = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/4. Zooplankton/Zoop_Perryetal_2021/MODIFIED"
path_evalout = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
path_ecospace_map = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap"

# -------------------------------
# Utility Functions
# -------------------------------

def get_csv_metadata(file_path):
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        return {"Column Names": df.columns.tolist(), "Number of Rows": len(df)}, df
    return {"Error": "File not found."}, None


def get_ecospace_times(file_path):
    if not os.path.exists(file_path):
        return {"Error": "File not found."}, None
    try:
        with nc.Dataset(file_path, 'r') as ds:
            time_var = ds.variables['time']
            units = time_var.units
            calendar = time_var.calendar if 'calendar' in time_var.ncattrs() else 'standard'
            times = nc.num2date(time_var[:], units=units, calendar=calendar)
            return {
                "Dimensions": {dim: len(ds.dimensions[dim]) for dim in ds.dimensions},
                "Variables": list(ds.variables.keys())
            }, times
    except OSError as e:
        return {"Error": str(e)}, None


# -------------------------------
# Load Input Data
# -------------------------------
ecospace_file_path = os.path.join(path_Ecospace, file_Ecospace)
zoop_file_path = os.path.join(path_Zoop, file_Zoop)
ecospace_map_path = os.path.join(path_ecospace_map, file_ecospace_map)

zoop_df = pd.read_csv(zoop_file_path)
ecospace_coords_df = pd.read_csv(ecospace_map_path)
model_lats = ecospace_coords_df['lat']
model_lons = ecospace_coords_df['lon']
model_rows = ecospace_coords_df['EWE_row']
model_cols = ecospace_coords_df['EWE_col']
model_deps = ecospace_coords_df['depth']
mask_ecospace = model_deps != 0

_, ecospace_times = get_ecospace_times(ecospace_file_path)

# -------------------------------
# Setup Output Columns
# -------------------------------
model_zoop_vars = ["ZC1-EUP", "ZC2-AMP", "ZC3-DEC", "ZC4-CLG", "ZC5-CSM",
                   "ZS2-CTH", "ZS3-CHA", "ZS4-LAR", "misc", "ZF1-ICT", "ZS1-JEL"]

for var in model_zoop_vars:
    zoop_df[f"EWE-{var}"] = np.nan

zoop_df['Date'] = pd.to_datetime(dict(year=zoop_df['Year'], month=zoop_df['Month'], day=15))
zoop_df['ecospace_closest_lat'] = np.nan
zoop_df['ecospace_closest_lon'] = np.nan
zoop_df['closest_ecospace_time'] = pd.NaT


# -------------------------------
# Main Matching Loop
# -------------------------------
with nc.Dataset(ecospace_file_path, 'r') as eco_ds:
    for index, row in zoop_df.iterrows():
        print(index)
        lat, lon, obs_time = row['Latitude.W.'], row['Longitude.N.'], row['Date']

        # Find nearest Ecospace cell
        eco_idx = find_nearest_point_from1D(lon, lat, model_lons, model_lats, mask_ecospace, dx=0.01)

        if eco_idx[0] is np.nan:
            continue

        zoop_df.at[index, 'ecospace_closest_lat'] = model_lats[eco_idx[0]]
        zoop_df.at[index, 'ecospace_closest_lon'] = model_lons[eco_idx[0]]
        row_idx, col_idx = model_rows[eco_idx[0]], model_cols[eco_idx[0]]

        # Find closest model time
        closest_time = find_closest_date(ecospace_times, obs_time)
        time_idx = np.where(ecospace_times == closest_time)[0][0]
        zoop_df.at[index, 'closest_ecospace_time'] = closest_time

        # Extract model data
        for var in model_zoop_vars:
            try:
                val = eco_ds.variables[var][time_idx, row_idx, col_idx]
                zoop_df.at[index, f"EWE-{var}"] = np.round(val, 3)
            except:
                continue

# -------------------------------
# Export Results
# -------------------------------
output_file = f"Zooplankton_matched_to_model_out_{ecospace_code}.csv"
zoop_df.to_csv(os.path.join(path_evalout, output_file), index=False)
print(f"Saved matched output to: {output_file}")
