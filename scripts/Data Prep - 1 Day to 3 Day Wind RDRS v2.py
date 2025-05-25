# Data Prep - 1 Day to 3 Day Wind RDRS
# Author: G Oldford
# Date: Aug 7-8, 2024 - updated May 2025
#
# Purpose:
#   Processes daily RDRS v2.1 wind data (10m u/v components) into 3-day maximum
#   wind speed and wind stress fields for use as forcing inputs to Ecospace and
#   Ecosim models.
#
# Key Features:
#   - Converts daily wind vectors to wind speed (sqrt(u^2 + v^2)).
#   - Aggregates data into 3-day blocks and extracts the maximum per block.
#   - Applies spatial clipping to Ecospace domain, with optional plume masking.
#   - Offers conditional override: user may define a constant value to be used
#     for all non-target months (e.g., retain wind data only for winter-spring).
#   - Saves Ecospace forcing files as `.asc` (wind speed and stress).
#   - Computes spatial averages and saves Ecosim time series as `.csv`.
#
# User Parameters:
#   - startyear, endyear         → Year range to process
#   - active_months              → Months where real wind data is preserved
#   - override_value             → Value used in non-active months
#   - output_suffix              → Used to tag output filenames for clarity
#   - apply_plume_mask           → Whether to apply plume region mask (True/False)
#
# Outputs:
#   - Ecospace `.asc` files for wind speed and stress (clipped & optionally masked)
#   - Ecosim `.csv` files containing 3-day-interval time series (expanded to 12x for monthly format)
#
# Notes:
#   - The script assumes 3-day blocks throughout the year, with a special case to lump the final days.
#   - Ecosim CSVs are written with a repeated row-per-month structure to match expected input formats.
#   - Designed for offline processing with pre-averaged daily RDRS inputs.
# ---

import os
import numpy as np
import netCDF4 as nc
import regex as re
from helpers import is_leap_year, saveASCFile, getDataFrame, buildSortableString
import csv
from datetime import datetime, timedelta

# === Paths and Configuration ===
wind_daily_p = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings/Wind_RDRS/RDRS21_NEMOgrid_wind_daily"
wind_daily_f = "RDRS21_NEMOgrid_wind_daily_{}.nc"
ecospace_out_p = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings/Wind_RDRS/Ecospace"
ecospace_out_speed_subdir = "speed"
ecospace_out_stress_subdir = "stress"
ecosim_out_p = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings/Wind_RDRS/Ecosim"
plume_dir = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//basemap//"
ecospacegrid_dir = plume_dir

scaling = "sqrt"
startyear = 1980
endyear = 2018
timestep = "3day"
apply_plume_mask = True  # Set to True if plume mask should be applied
lump_final = True  # Whether to lump final 1-2 days into last 3-day block
override_value = 1.0  # Value to use when overriding wind values
active_months = [1, 2, 3, 4, 12]  # Months to keep original wind data
output_suffix = "winter_spring"  # Label to apply to output files
ecospace_out_speed_subdir = ecospace_out_speed_subdir + "_" + output_suffix
ecospace_out_stress_subdir = ecospace_out_stress_subdir + "_" + output_suffix

# === ASC Parameters ===
saveTemplate = "{}_{}_{}.asc"
ASCheader = "ncols 93\nnrows 151\nxllcorner -123.80739\nyllcorner 48.29471\ncellsize 0.013497347\nNODATA_value 0.0"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102

# === Masks ===
dfPlumeMask = getDataFrame(os.path.join(plume_dir, "Fraser_Plume_Region.asc"), "-9999.00000000") == 1 if apply_plume_mask else None
ecospacegrid_f = os.path.join(ecospacegrid_dir, "ecospacedepthgrid.asc")
df_dep = getDataFrame(ecospacegrid_f, "-9999.00000000")
dfLandMask = df_dep == 0  # mask is where elev =0

# === Processing Loop ===
speed_ecosim_all = []
stress_ecosim_all = []
for iyear in range(startyear, endyear + 1):
    print(iyear)
    ds_var = nc.Dataset(os.path.join(wind_daily_p, wind_daily_f.format(iyear)))
    speed = np.sqrt(ds_var['u10m'][:] ** 2 + ds_var['v10m'][:] ** 2)
    leap = is_leap_year(iyear)
    days_in_year = 366 if leap else 365
    num_blocks = days_in_year // 3
    if not leap:
        num_blocks += 1

    j = 0
    for iday in range(1, num_blocks + 1):
        day_strt = (iday - 1) + (j * 2)
        day_end = day_strt + 2
        middle_day = day_strt + 2

        if lump_final and iday == 120:
            middle_day += 1 if not leap else 2
            varday1 = speed[day_strt:]
        else:
            if not leap and iday == num_blocks:
                day_end = day_strt + 1
                middle_day = day_strt + 1
            varday1 = speed[day_strt:day_end, :, :]

        speed_3day = np.ma.max(varday1, axis=0)
        stress_3day = np.square(np.clip(speed_3day, a_min=None, a_max=np.sqrt(np.finfo(speed_3day.dtype).max)))

        middle_date = datetime(iyear, 1, 1) + timedelta(days=middle_day - 1)
        if middle_date.month not in active_months:
            speed_3day[:, :] = override_value
            stress_3day[:, :] = override_value

        for var, name, fmt, subdir in [
            (speed_3day, f"RDRS_windspeed10m_{output_suffix}", '%0.2f', ecospace_out_speed_subdir),
            (stress_3day, f"RDRS_windstress10m_{output_suffix}", '%0.1f', ecospace_out_stress_subdir)]:
            path = os.path.join(ecospace_out_p, subdir)
            os.makedirs(path, exist_ok=True)
            fname = saveTemplate.format(name, iyear, buildSortableString(middle_day, 2))
            saveASCFile(os.path.join(path, fname), var, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe, fmt, ASCheader, dfPlumeMask, dfLandMask)

        for arr, store, fmt in [
            (speed_3day, speed_ecosim_all, '%0.2f'),
            (stress_3day, stress_ecosim_all, '%0.1f')]:
            val = np.mean(arr[arr != 0])
            precision = int(re.search(r'\.(\d+)f', fmt).group(1))
            store.append([iyear, middle_day, np.round(val, precision)])

        if lump_final and iday == 120:
            break
        j += 1


# === Write Ecosim CSV Files ===
for label, data, fmt in [
    (f"windspeed10m_{output_suffix}", speed_ecosim_all, '%0.2f'),
    (f"windstress10m_{output_suffix}", stress_ecosim_all, '%0.1f')]:
    out = [[i + 1] + row for i, row in enumerate([[idx + 1] + r for idx, r in enumerate(data) for _ in range(12)])]
    with open(os.path.join(ecosim_out_p, f"RDRS_{label}_{timestep}_{startyear}-{endyear}.csv"), 'w', newline='') as f:
        csv.writer(f).writerows([['threeday_yrmo', 'threeday_yr', 'year', 'dayofyear', f"RDRS_{label}"]] + out)

print('done the time series, on to the climatology')

# === Compute and Save Climatology ===
print("Computing climatology from all 3-day blocks...")

# Re-initialize arrays for stacking
all_speed_blocks = []
all_stress_blocks = []

for iyear in range(startyear, endyear + 1):
    ds_var = nc.Dataset(os.path.join(wind_daily_p, wind_daily_f.format(iyear)))
    speed = np.sqrt(ds_var['u10m'][:] ** 2 + ds_var['v10m'][:] ** 2)
    leap = is_leap_year(iyear)
    days_in_year = 366 if leap else 365
    num_blocks = days_in_year // 3 + (1 if not leap else 0)

    j = 0
    for iday in range(1, num_blocks + 1):
        day_strt = (iday - 1) + (j * 2)
        day_end = day_strt + 2
        middle_day = day_strt + 2

        if lump_final and iday == 120:
            middle_day += 1 if not leap else 2
            varday1 = speed[day_strt:]
        else:
            if not leap and iday == num_blocks:
                day_end = day_strt + 1
                middle_day = day_strt + 1
            varday1 = speed[day_strt:day_end, :, :]

        speed_3day = np.ma.max(varday1, axis=0)
        stress_3day = np.square(speed_3day)

        middle_date = datetime(iyear, 1, 1) + timedelta(days=middle_day - 1)
        if middle_date.month not in active_months:
            speed_3day[:, :] = override_value
            stress_3day[:, :] = override_value

        all_speed_blocks.append(speed_3day[np.newaxis, ...])
        all_stress_blocks.append(stress_3day[np.newaxis, ...])
        if lump_final and iday == 120:
            break
        j += 1

# Stack all and compute mean across time
clim_speed = np.mean(np.ma.concatenate(all_speed_blocks, axis=0), axis=0)
clim_stress = np.mean(np.ma.concatenate(all_stress_blocks, axis=0), axis=0)

# Save ASC files
for var, name, fmt, subdir in [
    (clim_speed, f"RDRS_CLIM_windspeed10m_{output_suffix}", '%0.2f', ecospace_out_speed_subdir),
    (clim_stress, f"RDRS_CLIM_windstress10m_{output_suffix}", '%0.1f', ecospace_out_stress_subdir)]:
    path = os.path.join(ecospace_out_p, subdir)
    os.makedirs(path, exist_ok=True)
    fname = f"{name}_climatology.asc"
    saveASCFile(os.path.join(path, fname), var, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe, fmt, ASCheader, dfPlumeMask, dfLandMask)

print("Climatology saved.")

print('DONE')
