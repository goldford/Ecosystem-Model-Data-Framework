import os
import numpy as np
import pandas as pd
from helpers import getDataFrame, buildSortableString

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
var_keys = [
    "varsalt1_PSU_0-10mAvg",
    "vartemp1_C_0-10mAvg",
    "PAR-VarZ-VarK"
]

NEMO_run = "216"
sigdig = 1  # number of decimal places to round
num_digits_ecospace_asc_month = 2

# --------------------------------------
# FILE PATHS
# --------------------------------------
forcing_root = "C:/Users/Greig/Documents/GitHub/Ecosystem-Model-Data-Framework/data/forcing"


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

    asc_dir = f"{forcing_root}/ECOSPACE_in_3day_vars_1980-2018/{var_key}/"
    if var_key == "PAR-VarZ-VarK":
        asc_dir = f"{forcing_root}/ECOSPACE_in_3day_PAR3_Sal4m_1980-2018/{var_key}/"

    output_csv = f"{forcing_root}/ECOSIM_in_3day_vars_1980-2018_fromASC_202506/ECOSIM_in_NEMO_{var_key}_1980-2018.csv"

    records = []
    for year in range(startyear, endyear + 1):
        print(year)
        is_leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
        nsteps = 120

        for iday in range(1, nsteps + 1):
            middle_day = (iday - 1) * 3 + 2
            if iday >= 120:
                middle_day += 2 if is_leap else 1

            fname = f"{var_key}_{year}_{buildSortableString(middle_day, num_digits_ecospace_asc_month)}.asc"
            fpath = os.path.join(asc_dir, fname)

            if not os.path.exists(fpath):
                print(f"Skipping missing file: {fpath}")
                continue

            data = getDataFrame_custom(fpath)

            if data.shape != combined_mask.shape:
                raise ValueError(f"Mismatch in shape between data ({data.shape}) and mask ({combined_mask.shape})")

            data = np.ma.masked_array(data, mask=combined_mask)
            mean_val = round(np.nanmean(data), sigdig)
            records.append([year, middle_day, mean_val])

    df = pd.DataFrame(records, columns=["year", "dayofyear", f"{out_code}{var_key}"])
    df["threeday_yrmo"] = range(1, len(df) + 1)
    df["threeday_yr"] = df.index + 1
    df = df[["threeday_yrmo", "threeday_yr", "year", "dayofyear", f"{out_code}{var_key}"]]
    df.to_csv(output_csv, index=False)
    print(f"Saved CSV for {var_key} to:\n{output_csv}")

