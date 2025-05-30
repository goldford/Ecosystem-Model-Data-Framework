# Modified version of Analysis - 1c - Bloom Timing vs Env Factors.py
# to produce multi-year environmental variable visualizations
# with climatological bands (like script 4a)

import os
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# --- User Configuration ---
NEMO_NC_ROOT = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/NEMO forcings/NC_3day/"
RDRS_NC_ROOT = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings/Wind_RDRS/Ecospace/NC_3day/"
REGIONS_ASC = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap/ecospace_regions_3day.asc"
FIG_OUT_PATH = "..//..//figs//"

REGION_ID = 2  # Central SoG
YEARS_TO_PLOT = list(range(2015, 2019))  # individual lines
CLIM_YEARS = list(range(1980, 2018))     # shaded bands

variables = {
    "Temp_0to10m": os.path.join(NEMO_NC_ROOT, f"Temp_0to10m_region{REGION_ID}.nc"),
    "Temp_30to40m": os.path.join(NEMO_NC_ROOT, f"Temp_30to40m_region{REGION_ID}.nc")
}

skiprows = 6
with open(REGIONS_ASC) as f:
    region_mask = np.loadtxt(f, skiprows=skiprows)

datasets = {var: xr.open_dataset(path) for var, path in variables.items()}

# --- Generate plot per variable ---
for var, ds in datasets.items():
    mask = (region_mask == REGION_ID) & (ds[var] > 0)
    var_avg = ds[var].where(mask).mean(dim=["row", "col"])
    df = var_avg.to_dataframe(name="value").reset_index()
    df["DoY"] = pd.to_datetime(df["time"]).dt.dayofyear
    df["Year"] = pd.to_datetime(df["time"]).dt.year

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))

    # Climatology
    df_clim = df[df["Year"].isin(CLIM_YEARS)]
    clim_stats = df_clim.groupby("DoY")["value"].agg(["mean", "quantile"])
    q10 = df_clim.groupby("DoY")["value"].quantile(0.10)
    q90 = df_clim.groupby("DoY")["value"].quantile(0.90)
    q25 = df_clim.groupby("DoY")["value"].quantile(0.25)
    q75 = df_clim.groupby("DoY")["value"].quantile(0.75)

    ax.fill_between(q10.index, q10, q90, color="0.6", alpha=0.2)
    ax.fill_between(q25.index, q25, q75, color="0.4", alpha=0.2)
    ax.plot(clim_stats.index, clim_stats["mean"], color="0.2", linewidth=2, label="Climatology Mean")

    # Plot selected years
    for yr in YEARS_TO_PLOT:
        df_yr = df[df["Year"] == yr]
        ax.plot(df_yr["DoY"], df_yr["value"], label=str(yr), linewidth=1)

    ax.set_title(f"{var} - Region {REGION_ID}")
    ax.set_xlabel("Day of Year")
    ax.set_ylabel(f"{var} (avg)")
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_OUT_PATH, f"{var}_envtimeseries_region{REGION_ID}_years{min(YEARS_TO_PLOT)}-{max(YEARS_TO_PLOT)}.png"))
    plt.show()
    plt.close()

print("DONE")
