"""
Script for Analyzing Biomass Anomalies of Ecospace Groups
---------------------------------------------------------------------------
This script processes Ecospace output netCDF files to extract seasonal mean biomass
anomalies (log-transformed, standardized) for a selection of groups
across a specified time period (1980–2018) and season
(default: Spring). The anomalies are computed relative to a climatology baseline.

The script compares interannual anomaly patterns with the North Pacific Gyre Oscillation
(NPGO) index via Pearson correlation and cross-correlation (CCF) analysis, optionally
accounting for lagged effects.

Outputs:
- Time series plots of biomass anomalies vs NPGO
- Correlation and cross-correlation statistics per group
"""

import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import pearsonr
from statsmodels.tsa.stattools import ccf
from helpers import read_sdomains
from matplotlib.path import Path
import yaml


# === CONFIGURATION ===
SCENARIO = "SC95_1"
FILENAME = "Scv95_1-All_Groups_20250506_1978-2018.nc"
NC_PATH = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//"
NPGO_PATH = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//evaluation//npgo.csv"

# === Subdomain Masking Configuration ===
USE_SUBDOMAIN_MASK = False
SUBDOMAIN_NAME = "SGS"  # Choose from keys in analysis_domains_jarnikova.yml
DOMAIN_FILE = "analysis_domains_jarnikova.yml"
DOMAIN_PATH = "C://Users//Greig//Documents//github//Ecosystem-Model-Data-Framework//data//evaluation"


target_groups = ['NK1-COH', 'NK2-CHI']
# target_groups = ['ZC1-EUP', 'ZC2-AMP', 'ZC3-DEC', 'ZC4-CLG', 'ZC5-CSM']
# target_groups = ['PZ1-CIL', 'PZ2-DIN', 'PZ3-HNF', 'PP1-DIA', 'PP2-NAN', 'PP3-PIC']
season_months = {
    "Spring": [3, 4, 5],
    "Summer": [6, 7, 8],
    "Spring_Summer": [3, 4, 5, 6, 7, 8],
    "Fall": [9, 10, 11],
    "Winter": [12, 1, 2],
    "Winter_Spring": [12, 1, 2, 3, 4, 5],
    "Spring_Summer_Fall": [3, 4, 5, 6, 7, 8, 9, 10, 11],
    "Year_Round": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
}
SEASON_USE_ECOL = "Year_Round" # choose one from season_months
SEASON_USE_NPGO = "Spring" # choose one from season_months
APPLY_LAG_PLOT = 3 # lag NPGO (forward?)
APPLY_LAG_STATS = 3
log_offset = 1e-3  # To safely log-transform

analysis_start_year = 1980
analysis_end_year = 2018
climatology_years = list(range(1980, 2018))


# === Read NPGO monthly data and compute annual average ===
df_npgo = pd.read_csv(NPGO_PATH)
df_npgo.columns = ["date", "npgo"]
df_npgo["date"] = pd.to_datetime(df_npgo["date"])
df_npgo = df_npgo[df_npgo["npgo"] > -99]
df_npgo["year"] = df_npgo["date"].dt.year
df_npgo["month"] = df_npgo["date"].dt.month

# Limit to certain months?
df_npgo_seasonal = df_npgo[df_npgo["month"].isin(season_months[SEASON_USE_NPGO])]
df_npgo_annual = df_npgo_seasonal.groupby("year")["npgo"].mean()
# we take more uears than needed from TS for lag analysis
df_npgo_annual = df_npgo_annual[(df_npgo_annual.index >= analysis_start_year-10) & (df_npgo_annual.index <= analysis_end_year+10)]


# === Utility ===
def safe_log(x):
    return xr.where(x > 0, np.log(x), np.nan)


# === Load Ecospace data ===
ds = xr.open_dataset(os.path.join(NC_PATH, FILENAME))

# === Generate spatial mask if enabled ===
if USE_SUBDOMAIN_MASK:
    domain_fp = os.path.join(DOMAIN_PATH, DOMAIN_FILE)
    sdomains = read_sdomains(domain_fp)

    if SUBDOMAIN_NAME not in sdomains:
        raise ValueError(f"Subdomain '{SUBDOMAIN_NAME}' not found in {DOMAIN_FILE}")

    polygon = Path(sdomains[SUBDOMAIN_NAME])
    lat = ds['lat'].values
    lon = ds['lon'].values

    # Flattened grid for point-in-polygon test
    points = np.vstack((lat.flatten(), lon.flatten())).T
    in_poly = polygon.contains_points(points).reshape(lat.shape)

    # Mask: inside polygon AND depth > 0
    mask = (ds['depth'].values > 0) & in_poly
else:
    mask = (ds['depth'].values > 0)

# anomaly calc
results = []
for group in target_groups:
    print(group)
    da = ds[group].where(mask)
    da_seasonal = da.sel(time=da['time.month'].isin(season_months[SEASON_USE_ECOL]))
    da_log = safe_log(da_seasonal)
    da_mean = da_log.mean(dim=('row', 'col'), skipna=True)
    series = da_mean.groupby('time.year').mean().to_pandas()
    series = series[(series.index >= analysis_start_year) & (series.index <= analysis_end_year)]
    clim = series.loc[series.index.isin(climatology_years)]
    clim_mean = clim.mean()
    anomaly = (series - clim.mean()) / clim.std()
    df = pd.DataFrame({
        'Year': series.index,
        'Group': group,
        'Anomaly': anomaly.values
    })
    results.append(df)

# === Combine and plot ===
anomaly_df = pd.concat(results)

# === Create panel of plots ===
n_groups = len(target_groups)
fig, axs = plt.subplots(n_groups, 1, figsize=(8, 3 * n_groups), sharex=True)
if n_groups == 1:
    axs = [axs]  # make axs iterable

letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']
df_npgo_annual_plot = df_npgo_annual.shift(APPLY_LAG_PLOT)
df_npgo_annual_plot = df_npgo_annual_plot[(df_npgo_annual_plot.index >= analysis_start_year) & (df_npgo_annual_plot.index <= analysis_end_year)]


for i, group in enumerate(target_groups):
    sub = anomaly_df[anomaly_df['Group'] == group].set_index('Year')
    colors = ['blue' if x > 0 else 'red' for x in sub['Anomaly']]
    axs[i].bar(sub.index, sub['Anomaly'], color=colors, label='Model Anomaly')
    axs[i].plot(df_npgo_annual_plot.index, df_npgo_annual_plot.values, color='black', label='NPGO Index')
    axs[i].axhline(0, color='grey', linestyle='--')
    axs[i].text(0.05, 0.95, f'{letters[i]} {group}', transform=axs[i].transAxes, fontsize=12, verticalalignment='top')
    # axs[i].set_title(f'{letters[i]} {group}', loc='left', fontsize=12)
    axs[i].set_ylabel('Standardised Annual Anomaly')
    axs[i].legend(loc='upper right')

axs[-1].set_xlabel('Year')
plt.tight_layout()
plt.show()

# Lag NPGO? optional based on CCF
df_npgo_annual_lagged = df_npgo_annual.shift(APPLY_LAG_STATS)

for group in target_groups:
    print(group)
    sub = anomaly_df[anomaly_df['Group'] == group].set_index('Year')

    # === Correlation and CCF ===
    common_years = sub.index.intersection(df_npgo_annual_lagged.index)
    model_vals = sub.loc[common_years, 'Anomaly']
    npgo_vals = df_npgo_annual_lagged.loc[common_years]

    if len(common_years) > 2:
        r, p = pearsonr(model_vals, npgo_vals)
        ccf_vals = ccf(npgo_vals.fillna(0), model_vals.fillna(0), adjusted=False)[:5]
        print(f"{group} Pearson r = {r:.2f}, p = {p:.3f}")
        print(f"{group} Max CCF = {np.max(np.abs(ccf_vals)):.2f} at lag {np.argmax(np.abs(ccf_vals))}")
    else:
        print(f"{group} Not enough data for correlation analysis")

print('done')