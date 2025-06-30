import os
import numpy as np
import pandas as pd
from helpers import getDataFrame, buildSortableString

"""
Create Ecosim forcing time series
by: G Oldford, 2023-2025
----------------------------------------------------------
It can be time consuming to prep both sim and space files at once. 
This script is used in case Ecospace files are prepped
but the Ecosim files have not yet been prepped. 

Inputs:
- Ecospace ASC files

Notes:
"""

# --------------------------------------
# CONFIGURATION
# --------------------------------------
startyear = 1980
endyear = 2018
timestep = "3day"
out_code = "var"
# wind already done
timestep = "3day"
out_code = "var"
anomalies = True
var_keys = [
    #"vartemp_at150m",
    "vartemp_avg150toBot"#,
    # "vartemp_avg150toBot"#,
    # "RDRS_windstress10m",
    # "varsalt1_PSU_0-10mAvg",
    # "vartemp1_C_0-10mAvg",
    # "PAR-VarZ-VarK"
]

NEMO_run = "216"
sigdig = 1  # number of decimal places to round
num_digits_ecospace_asc_month = 2

# --------------------------------------
# FILE PATHS
# --------------------------------------
forcing_root = "C:/Users/Greig/Documents/GitHub/Ecosystem-Model-Data-Framework/data/forcing"
forcing_root_RDRS = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings"


# Mask files (must be 151x93 trimmed ASC format used by saveASCFile)
land_mask_file = "C:/Users/Greig/Documents/GitHub/Ecosystem-Model-Data-Framework/data/basemap/ecospacedepthgrid.asc"
plume_mask_file = "C:/Users/Greig/Documents/GitHub/Ecosystem-Model-Data-Framework/data/basemap/Fraser_Plume_Region.asc"


plume_mask = getDataFrame(plume_mask_file, "-9999.00000000") == 1
land_mask = getDataFrame(land_mask_file, "-9999.00000000") == 0
combined_mask = plume_mask | land_mask

from helpers import buildSortableString

# --------------------------------------
# CUSTOM DATA READER
# --------------------------------------
def getDataFrame_custom(f_n, skiprows=6):
    with open(f_n) as f:
        data = np.loadtxt(f, skiprows=skiprows)

    # Mask invalid values
    data[(data == -9999.0) | (data == 0.0)] = np.nan

    # Fix bottom-left corner issue
    data[140:, :15] = np.nan

    return data

# --------------------------------------
# LOOP OVER VARIABLES
# --------------------------------------
for var_key in var_keys:
    print(f"Processing variable: {var_key}")

    if anomalies:
        asc_dir = f"{forcing_root}/ECOSPACE_in_3day_vars_1980-2018/{var_key}/anomalies/"
    else:
        asc_dir = f"{forcing_root}/ECOSPACE_in_3day_vars_1980-2018/{var_key}/"

    if var_key == "PAR-VarZ-VarK":
        asc_dir = f"{forcing_root}/ECOSPACE_in_3day_PAR3_Sal4m_1980-2018/{var_key}/"

    if var_key == "RDRS_windstress10m":
        asc_dir = f"{forcing_root_RDRS}/Wind_RDRS/Ecospace/stress_/"

    if var_key == "vartemp_at150m" or var_key == "vartemp_avg150toBot":
        num_digits_ecospace_asc_month = 3

    if anomalies:
        var_key = var_key + "_anom"

    records = []
    for year in range(startyear, endyear + 1):
        print(year)
        is_leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
        nsteps = 120

        for iday in range(1, nsteps + 1):
            middle_day = (iday - 1) * 3 + 2
            if iday >= 120:
                middle_day += 2 if is_leap else 1

            if var_key == "RDRS_windstress10m": # correcting for a file naming glitch!
                fname = f"{var_key}__{year}_{buildSortableString(middle_day, num_digits_ecospace_asc_month)}.asc"
            else:
                fname = f"{var_key}_{year}_{buildSortableString(middle_day, num_digits_ecospace_asc_month)}.asc"

            fpath = os.path.join(asc_dir, fname)

            if not os.path.exists(fpath):
                if middle_day >= 360:
                    middle_day = 363 # some vars have inconsistently named final day-of-year
                    if var_key == "RDRS_windstress10m":  # correcting for a file naming glitch!
                        fname = f"{var_key}__{year}_{buildSortableString(middle_day, num_digits_ecospace_asc_month)}.asc"
                    else:
                        fname = f"{var_key}_{year}_{buildSortableString(middle_day, num_digits_ecospace_asc_month)}.asc"
                    fpath = os.path.join(asc_dir, fname)
                    if not os.path.exists(fpath): # some files have the final time step named 363, 362, etc
                        print(f"Skipping missing file: {fpath}")
                        continue
                else:
                    print(f"Skipping missing file: {fpath}")
                    continue

            data = getDataFrame_custom(fpath)

            if data.shape != combined_mask.shape:
                raise ValueError(f"Mismatch in shape between data ({data.shape}) and mask ({combined_mask.shape})")

            data = np.ma.masked_array(data, mask=combined_mask)
            mean_val = round(np.nanmean(data), sigdig)
            records.append([year, middle_day, mean_val])

    output_csv = f"{forcing_root}/ECOSIM_in_3day_vars_1980-2018_fromASC_202506/ECOSIM_in_NEMO_{var_key}_1980-2018.csv"

    df = pd.DataFrame(records, columns=["year", "dayofyear", f"{out_code}{var_key}"])
    df["threeday_yrmo"] = range(1, len(df) + 1)
    df["threeday_yr"] = df.index + 1
    df = df[["threeday_yrmo", "threeday_yr", "year", "dayofyear", f"{out_code}{var_key}"]]
    df.to_csv(output_csv, index=False)
    print(f"Saved CSV for {var_key} to:\n{output_csv}")

    # CREATE DUMMY YEARS 1978 and 1979 for spin-up
    df = pd.read_csv(output_csv)
    df_1980 = df[df["year"] == 1980].copy()

    # Create 1979 and 1978 versions by copying and changing the year
    df_1979 = df_1980.copy()
    df_1979["year"] = 1979
    df_1979["threeday_yrmo"] = df_1979["threeday_yrmo"] - len(df_1980)
    df_1979["threeday_yr"] = df_1979["threeday_yr"] - len(df_1980)

    df_1978 = df_1980.copy()
    df_1978["year"] = 1978
    df_1978["threeday_yrmo"] = df_1979["threeday_yrmo"] - len(df_1980)
    df_1978["threeday_yr"] = df_1979["threeday_yr"] - len(df_1980)

    # Concatenate in new order: 1978, 1979, original
    df_extended = pd.concat([df_1978, df_1979, df], ignore_index=True)
    df_extended.reset_index(drop=True, inplace=True)

    # Recalculate new sequential IDs
    df_extended["threeday_yrmo"] = range(1, len(df_extended) + 1)
    df_extended["threeday_yr"] = df_extended["threeday_yrmo"]

    # Replace first 12 time steps with a single climatological mean over 1980-2018
    data_col = f"{out_code}{var_key}"
    mean_val = round(df_extended[(df_extended["year"] >= 1980) & (df_extended["year"] <= 2018)][data_col].mean(),2)
    df_extended.loc[:11, data_col] = mean_val  # Overwrite first 12 rows with same value

    # Save updated CSV
    df_extended.to_csv(output_csv, index=False)
    print(f"Prepended 1978 and 1979 using 1980 data and applied uniform climatological mean to first 12 steps for {var_key}")
