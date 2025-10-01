"""
 Script for Exporting PP B Anoms as 1D TS from Ecospace Outputs
----------------------------------------------------------
Steps

Outputs:
-
-
"""

import xarray as xr
import numpy as np
import pandas as pd
import os

# Set paths and configuration
SCENARIO = 'LTL_SC85_1'
ECOSPACE_CODE = "Scv85_1-All_Groups_20250506"
NETCDF_FILE = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//" \
              + ECOSPACE_CODE +"_1978-2018.nc"
MASK_PATH = "..//..//data/evaluation//suchy_ecospace_mask.nc"
PHYTO_GROUPS = ["PP3-PIC"]
PP_GRP_TAG = "PP3-PIC"
OUTPUT_CSV = "..//..//data//forcing//phytoplankton_anomalies_" + SCENARIO + "_" + PP_GRP_TAG + ".csv"

# Load datasets
ds = xr.open_dataset(NETCDF_FILE)
mask = xr.open_dataset(MASK_PATH)['mask']

# Subset to years of interest
ds = ds.sel(time=slice("1980-01-01", "2018-12-31"))

# Apply spatial mask and calculate domain average for each timestep
def compute_masked_mean(ds, mask, groups):
    masked = sum(ds[var].where(mask) for var in groups if var in ds)
    return masked.mean(dim=["row", "col"], skipna=True)

domain_mean_ts = compute_masked_mean(ds, mask, PHYTO_GROUPS).to_dataframe(name="biomass").reset_index()
domain_mean_ts['time'] = pd.to_datetime(domain_mean_ts['time'])

# Step 1: Compute monthly climatology over all years
domain_mean_ts['month'] = domain_mean_ts['time'].dt.month
monthly_clim = domain_mean_ts.groupby('month')['biomass'].mean()

# Step 2: Normalize climatology so mean = 1 → seasonal shape
seasonal_shape = monthly_clim / monthly_clim.mean()

# Step 3: Calculate interannual anomalies for each month
domain_mean_ts['month'] = domain_mean_ts['time'].dt.month
domain_mean_ts['year'] = domain_mean_ts['time'].dt.year
monthly_avg = domain_mean_ts.groupby([domain_mean_ts['time'].dt.to_period('M')]).mean(numeric_only=True)
monthly_avg.index = monthly_avg.index.to_timestamp()
monthly_avg['month'] = monthly_avg.index.month

# Merge monthly climatology into full monthly average dataframe
monthly_avg['climatology'] = monthly_avg['month'].map(monthly_clim)
monthly_avg['anomaly'] = monthly_avg['biomass'] - monthly_avg['climatology']

# Step 4: Normalize anomalies to 0–2 scale centered at 1 - EwE requires this
zmin = monthly_avg['anomaly'].min()
zmax = monthly_avg['anomaly'].max()
renorm = 2 * (monthly_avg['anomaly'] - zmin) / (zmax - zmin)
z0_to = 2 * (0 - zmin) / (zmax - zmin)
monthly_avg['interannual_multiplier'] = (renorm - (z0_to - 1)).clip(0, 2)

# Step 5: Combine seasonal shape with interannual multiplier
monthly_avg['seasonal_shape'] = monthly_avg['month'].map(seasonal_shape)
monthly_avg['final_forcing_multiplier'] = monthly_avg['seasonal_shape'] * monthly_avg['interannual_multiplier']

# Save to CSV
monthly_avg[['biomass', 'climatology', 'anomaly', 'interannual_multiplier', 'seasonal_shape', 'final_forcing_multiplier']].to_csv(OUTPUT_CSV)

print("done!")