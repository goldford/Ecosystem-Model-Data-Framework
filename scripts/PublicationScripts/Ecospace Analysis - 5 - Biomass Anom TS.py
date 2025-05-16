# =============================================================
# Script: Compare Ecospace Anomalies with NPGO and Zooplankton
# Author: G. Oldford
# Description:
#   - Loads and computes normalized biomass anomalies from Ecospace model output
#   - Compares model anomalies to NPGO climate index and observational zooplankton anomalies
#   - Computes Pearson correlation and CCF statistics for each comparison
#   - Computes annual anomalies from seasonal spring/summer means, using 1996–2010 climatology
#   - Now supports winter and fall anomaly output as well
#   - Computes matching seasonal anomalies for NPGO data
#   - Supports cached processing to improve performance
#
#   To do / Notes:
#   - the anomaly fits to perry et al data (2021) downloaded from pacea is poor
#   - will try to do point wise model obs matching then, probably disregard this
# =============================================================
import pandas as pd
import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from statsmodels.tsa.stattools import ccf
import pickle

# === User toggles ===
log_transform = False

# === Paths and constants ===
SCENARIO = "SC88_1"
FILENM = "Scv88_1-All_Groups_20250506_1978-2018.nc"
NC_PATH = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//"
netcdf_path = os.path.join(NC_PATH, FILENM)
npgo_path = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//evaluation//npgo.csv"
zoop_path = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//evaluation//zooplankton_sog_anomalies.csv"
analysis_start_year = 1980
analysis_end_year = 2018
climatology_years = list(range(1996, 2011))

if log_transform:
    lg_tag = "LogTrns"
else:
    lg_tag = "NoLogTrns"
npgo_cache_path = f"cached_anomalies_for_npgo_{SCENARIO}_{lg_tag}.pkl"
zoop_cache_path = f"cached_anomalies_for_zoop_{SCENARIO}_{lg_tag}.pkl"

season_months = {
    "spring": [3, 4, 5],
    "summer": [6, 7, 8],
    "fall": [9, 10, 11],
    "winter": [12, 1, 2]
}
zoop_mapping = {
    "amphipods_gammarid": "ZC2-AMP",
    "benthic_larvae": "ZC3-DEC",
    "euphausiids": "ZC1-EUP",
    "fish": "ZF1-ICT",
    "ctenophora": "ZS2-CTH",
    "chaetognatha": "ZS3-CHA",
    "larvacea": "ZS4-LAR",
    "calanoid_copepods_large": "ZC4-CLG",
    "non_calanoid_copeopods": "ZC5-CSM",
    "scyphozoa": "ZS1-JEL"
}

# === Utility function ===
def safe_log(x):
    return xr.where(x > 0, np.log(x), np.nan)


# === Function: Compute anomalies for comparison with NPGO ===
def compute_anomalies_for_npgo(netcdf_path, climatology_years, season_months, start_year, end_year, log_transform=True):
    print('computing anomalies for npgo')
    ds = xr.open_dataset(netcdf_path)
    mask = ds['depth'] > 0
    exclude_vars = {'lat', 'lon', 'depth', 'EWE_col', 'EWE_row', 'NEMO_col', 'NEMO_row'}
    group_names = [var for var in ds.data_vars if var not in exclude_vars]

    seasonal_means = {}
    for season, months in season_months.items():
        print(season)
        seasonal_means[season] = {}
        for group in group_names:
            print(group)
            da = ds[group].where(mask)
            da = da.sel(time=da['time.month'].isin(months))
            if log_transform:
                da = safe_log(da)
            spatial_mean = da.mean(dim=("row", "col"), skipna=True)
            seasonal_mean = spatial_mean.groupby('time.year').mean()
            seasonal_means[season][group] = seasonal_mean

This is where I think the 'total' zoop B should be done differently. Should not be totalling the anomalies!
Should total the B first - either log B or B, doesn't matter'

    spring_df = pd.DataFrame({g: seasonal_means['spring'][g].to_pandas() for g in group_names})
    summer_df = pd.DataFrame({g: seasonal_means['summer'][g].to_pandas() for g in group_names})
    fall_df = pd.DataFrame({g: seasonal_means['fall'][g].to_pandas() for g in group_names})
    winter_df = pd.DataFrame({g: seasonal_means['winter'][g].to_pandas() for g in group_names})

    combined_df = pd.concat([spring_df, summer_df, fall_df, winter_df], axis=1, keys=['spring', 'summer', 'fall', 'winter'])
    annual_df = combined_df.groupby(level=1, axis=1).mean()
    annual_df = annual_df[(annual_df.index >= start_year) & (annual_df.index <= end_year)]
    clim = annual_df.loc[annual_df.index.isin(climatology_years)]
    anomaly_df = (annual_df - clim.mean()) / clim.std()

    return {
        'spring': seasonal_means['spring'],
        'summer': seasonal_means['summer'],
        'fall': seasonal_means['fall'],
        'winter': seasonal_means['winter'],
        'annual': anomaly_df
    }


# === Function: Compute anomalies for comparison with zooplankton ===
def compute_anomalies_for_zoop(netcdf_path, climatology_years, season_months, start_year, end_year, log_transform=True):
    print('computing zoop anom comparison')
    ds = xr.open_dataset(netcdf_path)
    mask = ds['depth'] > 0
    exclude_vars = {'lat', 'lon', 'depth', 'EWE_col', 'EWE_row', 'NEMO_col', 'NEMO_row'}
    group_names = [var for var in ds.data_vars if var not in exclude_vars]

    seasonal_means = {}
    for season in ['spring', 'summer', 'fall', 'winter']:
        print(season)
        seasonal_means[season] = {}
        months = season_months[season]
        for group in group_names:
            print(group)
            da = ds[group].where(mask)
            da = da.sel(time=da['time.month'].isin(months))
            if log_transform:
                da = safe_log(da)
            spatial_mean = da.mean(dim=("row", "col"), skipna=True)
            seasonal_mean = spatial_mean.groupby('time.year').mean()
            seasonal_means[season][group] = seasonal_mean

    spring_df = pd.DataFrame({g: seasonal_means['spring'][g].to_pandas() for g in group_names})
    summer_df = pd.DataFrame({g: seasonal_means['summer'][g].to_pandas() for g in group_names})
    fall_df = pd.DataFrame({g: seasonal_means['fall'][g].to_pandas() for g in group_names})
    winter_df = pd.DataFrame({g: seasonal_means['winter'][g].to_pandas() for g in group_names})

    combined_df = pd.concat([spring_df, summer_df, fall_df, winter_df], axis=1, keys=['spring', 'summer', 'fall', 'winter'])
    annual_df = combined_df.groupby(level=1, axis=1).mean()
    annual_df = annual_df[(annual_df.index >= start_year) & (annual_df.index <= end_year)]
    clim = annual_df.loc[annual_df.index.isin(climatology_years)]
    anomaly_df = (annual_df - clim.mean()) / clim.std()

    # Compute total biomass by summing mapped groups
    mapped_groups = list(zoop_mapping.values())
    total_biomass = annual_df[mapped_groups].sum(axis=1)
    total_clim = total_biomass.loc[total_biomass.index.isin(climatology_years)]
    total_anomaly = (total_biomass - total_clim.mean()) / total_clim.std()
    anomaly_df['ZOO-TOTAL'] = total_anomaly

    return {
        'spring': seasonal_means['spring'],
        'summer': seasonal_means['summer'],
        'fall': seasonal_means['fall'],
        'winter': seasonal_means['winter'],
        'annual': anomaly_df
    }


# === Check for cached anomalies ===
if os.path.exists(npgo_cache_path) and os.path.exists(zoop_cache_path):
    with open(npgo_cache_path, 'rb') as f1:
        npgo_anoms = pickle.load(f1)
    with open(zoop_cache_path, 'rb') as f2:
        zoop_anoms = pickle.load(f2)
    print(f"Loaded cached anomalies: {npgo_cache_path}, {zoop_cache_path}")
else:
    npgo_anoms = compute_anomalies_for_npgo(netcdf_path, climatology_years, season_months, analysis_start_year, analysis_end_year, log_transform)
    zoop_anoms = compute_anomalies_for_zoop(netcdf_path, climatology_years, season_months, analysis_start_year, analysis_end_year, log_transform)
    with open(npgo_cache_path, 'wb') as f1:
        pickle.dump(npgo_anoms, f1)
    with open(zoop_cache_path, 'wb') as f2:
        pickle.dump(zoop_anoms, f2)
    print(f"Saved computed anomalies to: {npgo_cache_path}, {zoop_cache_path}")

npgo_anomalies = npgo_anoms['annual']
zoop_anomalies = zoop_anoms['annual']


# === Load and process NPGO data ===
df = pd.read_csv(npgo_path)
df.columns = ["date", "npgo"]
df["date"] = pd.to_datetime(df["date"])
df = df[df["npgo"] > -99]
df["year"] = df["date"].dt.year
df["month"] = df["date"].dt.month
ss_df = df[df["month"].isin(season_months["spring"] + season_months["summer"])]
npgo_index = ss_df.groupby("year")["npgo"].mean()
npgo_index = npgo_index[npgo_index.index.isin(npgo_anomalies.index)]

# === NPGO Panel Plot ===
years = npgo_anomalies.index
n_groups = len(npgo_anomalies.columns)
n_cols = 4
n_rows = int(np.ceil(n_groups / n_cols))
fig, axs = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows), sharex=True)
axs = axs.flatten()

for i, group in enumerate(npgo_anomalies.columns):
    anomalies = npgo_anomalies[group]
    common_years = anomalies.index.intersection(npgo_index.index)
    r, p = pearsonr(anomalies.loc[common_years], npgo_index.loc[common_years])
    ccf_vals = ccf(npgo_index.fillna(0), anomalies.fillna(0), adjusted=False)[:5]
    max_lag = int(np.argmax(np.abs(ccf_vals)))
    max_val = ccf_vals[max_lag]
    colors = ["blue" if val >= 0 else "red" for val in anomalies]

    axs[i].bar(years, anomalies, color=colors)
    axs[i].plot(npgo_index.index, npgo_index.values, color='black', linestyle='-', linewidth=1.5)
    axs[i].axhline(0, color='grey', linestyle='--', linewidth=0.8)
    axs[i].set_title(f"{group}\n(r={r:.2f}, p={p:.2f}, lag={max_lag}, ccf={max_val:.2f})")
    axs[i].set_ylim([-3, 3])

for j in range(n_groups, len(axs)):
    axs[j].axis('off')

fig.suptitle("Model Anomalies with NPGO Overlay", fontsize=16)
fig.supxlabel("Year")
fig.supylabel("Normalized Anomaly")
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()

# === Load and process zooplankton data ===
zoop_df = pd.read_csv(zoop_path)
zoop_df = zoop_df[(zoop_df['year'] >= 1996) & (zoop_df['year'] <= analysis_end_year)]
zoop_df.set_index("year", inplace=True)
zoop_df = zoop_df[zoop_df.index.isin(zoop_anomalies.index)]

zoop_mapping2 = {
    "amphipods_gammarid": "ZC2-AMP",
    "benthic_larvae": "ZC3-DEC",
    "euphausiids": "ZC1-EUP",
    "fish": "ZF1-ICT",
    "ctenophora": "ZS2-CTH",
    "chaetognatha": "ZS3-CHA",
    "larvacea": "ZS4-LAR",
    "calanoid_copepods_large": "ZC4-CLG",
    "non_calanoid_copeopods": "ZC5-CSM",
    "scyphozoa": "ZS1-JEL",
    "total_biomass": "ZOO-TOTAL"
}


zoop_mapping2 = {
    "amphipods_hyperiid": "ZC2-AMP",
    "benthic_larvae": "ZC3-DEC",
    "euphausiids": "ZC1-EUP",
    "fish": "ZF1-ICT",
    "calanoid_copepods_large": "ZC4-CLG",
    "total_biomass": "ZOO-TOTAL"
}


zoop_anoms = {model_group: zoop_df[obs_col] for obs_col, model_group in zoop_mapping2.items() if
              obs_col in zoop_df.columns}
match_groups = [g for g in zoop_anoms if g in zoop_anomalies.columns]
n_groups = len(match_groups)
n_cols = 4
n_rows = int(np.ceil(n_groups / n_cols))
fig, axs = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows), sharex=True)
axs = axs.flatten()

for i, group in enumerate(match_groups):
    model = zoop_anomalies[group]
    obs = zoop_anoms[group].reindex(model.index)
    valid_idx = model.notna() & obs.notna()
    if valid_idx.sum() > 1:
        r, p = pearsonr(model[valid_idx], obs[valid_idx])
    else:
        r, p = np.nan, np.nan
    ccf_vals = ccf(obs.fillna(0), model.fillna(0), adjusted=False)[:5]
    max_lag = int(np.argmax(np.abs(ccf_vals)))
    max_val = ccf_vals[max_lag]
    colors = ["blue" if val >= 0 else "red" for val in model]

    axs[i].bar(model.index, model, color=colors)
    axs[i].plot(obs.index, obs.values, color='black', linestyle='-', linewidth=1.5)
    axs[i].axhline(0, color='grey', linestyle='--', linewidth=0.8)
    axs[i].set_title(f"{group}\n(r={r:.2f}, p={p:.2f}, lag={max_lag}, ccf={max_val:.2f})")
    axs[i].set_ylim([-3, 3])

for j in range(n_groups, len(axs)):
    axs[j].axis('off')

fig.suptitle("Model Anomalies with Zooplankton Overlay", fontsize=16)
fig.supxlabel("Year")
fig.supylabel("Normalized Anomaly")
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()

print("Done")
