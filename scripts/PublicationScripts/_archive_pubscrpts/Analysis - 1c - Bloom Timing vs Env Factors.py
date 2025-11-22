# Stacked Bloom Forcings - Visualization Script
# ------------------------------------------------------------
# Purpose:
#   Generate stacked time series plots of key environmental forcings (PAR, MixingZ, Wind Stress, Temp)
#   for the central Strait of Georgia (Region 2 in Ecospace).
#   Each year's data is plotted in a separate row, with bloom dates (from Suchy et al., 2021) overlaid.
#   The visualization helps identify which factors may limit or enable bloom initiation.
#
# Input:
#   - NetCDF files (see archived Visuals - ASC Forcings by Region.py)
#     containing 3-day averaged time series for each variable:
#       - PAR.nc (Photosynthetically Active Radiation)
#       - MixingZ.nc (vertical mixing)
#       - Wind_Stress_10m_RDRS.nc (surface wind stress)
#       - Temp_0to10m.nc (0–10 m average temperature)
#       - etc
#   - A region map file (ecospace_regions_3day.asc) to identify Region 2 cells.
#   - A predefined list of bloom dates and DOYs from Suchy et al. (2003–2016).
#
# Output:
#   - A multi-panel figure (stacked per year and per variable) showing each variable’s time series.
#   - Red dashed vertical lines indicate the bloom dates.
#   - Output figure saved to '../figs/stacked_timeseries_bloom_overlay.png'
#
# Author: G. Oldford
# Date: Initial July 2024, Refactor May 2025
# Refactored script for stacked time series visualization with bloom overlays
# Focus: Stacked plots for each year with PAR, Mixing, Wind Stress, and Temp

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")
from datetime import datetime, timedelta
import matplotlib.dates as mdates

# --- User Configuration ---
NEMO_NC_ROOT = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/NEMO forcings/NC_3day/"
RDRS_NC_ROOT = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings/Wind_RDRS/Ecospace/NC_3day/"
REGIONS_ASC = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap/ecospace_regions_3day.asc"
NEMO_ASCs_ROOT = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/forcing/"
FIG_OUT_PATH = "..//..//figs//"

BLOOM_SOURCE = "suchy" # suchy or C09

REGION_ID = 2  # Central SoG (this method can reduce size of NC's)

variables = {
    "PAR": os.path.join(NEMO_NC_ROOT, "PAR" + f"_region{REGION_ID}.nc"),
    # "MixingZ": os.path.join(NEMO_NC_ROOT, "MixingZ" + f"_region{REGION_ID}.nc"),
    "Wind_Stress_10m_RDRS": os.path.join(RDRS_NC_ROOT, "Wind_Stress_10m_RDRS" + f"_region{REGION_ID}.nc"),
    "Temp_0to10m": os.path.join(NEMO_NC_ROOT, "Temp_0to10m" + f"_region{REGION_ID}.nc"),
    #"Temp_30to40m": os.path.join(NEMO_NC_ROOT, "Temp_30to40m" + f"_region{REGION_ID}.nc")
}

skiprows = 6 # for ASC header


# bloom DOYs and corresponding dates
if BLOOM_SOURCE == "suchy":
    doy_blooms = [100, 68, 50, 83, 115, 115, 100, 100, 92, 88, 92, 92, 55, 77] #suchy satellite dates, 2003 - 2016
    doy_bloom_strt = 2003
elif BLOOM_SOURCE == "C09":
    doy_blooms = [94, 78, 81, 82, 87, 76, 75, 87, 99, 76,
                 81, 78, 77, 55, 86, 86, 67, 87, 77, 104,
                 66, 70, 92, 86, 81, 62, 88, 100, 90, 97,
                 104, 103, 98, 88, 86, 77, 88, 104, 77]
    doy_bloom_strt = 1980
else:
    print("error with bloom reference data")
    exit()

bloom_dates = [datetime(doy_bloom_strt + i, 1, 1) + timedelta(days=d - 1) for i, d in enumerate(doy_blooms)]

# --- Load region mask ---
with open(REGIONS_ASC) as f:
    region_mask = np.loadtxt(f, skiprows=skiprows)

# --- Load datasets ---
datasets = {var: xr.open_dataset(path) for var, path in variables.items()}

# --- Create one figure per year ---
for year_idx, bloom_date in enumerate(bloom_dates):
    year = bloom_date.year
    print(year)
    start = datetime(year, 1, 1)
    end = datetime(year, 12, 31)

    fig, axes = plt.subplots(len(variables), 1, figsize=(12, 8), sharex=True)

    for var_idx, (var, ds) in enumerate(datasets.items()):
        mask = (region_mask == REGION_ID) & (ds[var] > 0)
        var_avg = ds[var].where(mask).mean(dim=["row", "col"])

        yearly_data = var_avg.sel(time=slice(start, end))
        times = yearly_data["time"].values
        doy = [pd.Timestamp(t).dayofyear for t in times]

        ax = axes[var_idx]
        ax.plot(doy, yearly_data, label=var, color='tab:blue', linewidth=0.6)
        ax.set_xticks(range(0, 366, 20))
        ax.set_xticklabels([str(d) for d in range(0, 366, 20)], fontsize=8)

        ax.axvline(bloom_date.timetuple().tm_yday, color='red', linestyle='--', linewidth=1)

        try:
            bloom_value = var_avg.sel(time=bloom_date, method='nearest').values.item()
            ax.axhline(bloom_value, color='red', linestyle='--', linewidth=0.6)
            ax.text(bloom_date.timetuple().tm_yday, bloom_value, f"{bloom_value:.1f}",
                    color='red', fontsize=8, va='bottom', ha='left')
        except Exception:
            pass

        try:
            prior_window = var_avg.sel(time=slice(bloom_date - timedelta(days=6), bloom_date - timedelta(days=1)))
            if prior_window.size > 0:
                mean_prior = prior_window.mean().values.item()
                ax.axhline(mean_prior, color='black', linestyle='--', linewidth=0.6)
                ax.text((bloom_date - timedelta(days=20)).timetuple().tm_yday, mean_prior, f"{mean_prior:.1f}", color='black', fontsize=8,
                        va='bottom', ha='left')

        except Exception:
            pass

        ax.set_ylabel(f"{var}", fontsize=9)
        ax.grid(True, linewidth=0.2)

    axes[-1].set_xlabel("Day of Year", fontsize=10)
    fig.suptitle(f"{year} Forcings with Bloom Overlay", fontsize=14)
    # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(FIG_OUT_PATH, f"forcing_by_year_{year}_vs_{BLOOM_SOURCE}.png"), dpi=300)
    plt.close(fig)

print("DONE")
