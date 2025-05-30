import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from statsmodels.tsa.stattools import ccf
import xarray as xr
import os

# === USER CONFIGURATION ===
data_path = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
input_file = "Zooplankton_matched_to_model_out_Scv114_1.csv"  # used when use_raw_ecospace = False
ecospace_run = "Scv114_1"
raw_ecospace_nc = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//Scv114_1-All_Groups_20250523_1978-2018.nc"
npgo_file = "npgo.csv"
output_fig = f"../../figs/{ecospace_run}_anom_panel_Zoop_vs_Model_with_NPGO.png"
total_fig = f"../../figs/{ecospace_run}_anom_TOTAL_stacked_Model_vs_Obs.png"
# groups = ['ZC1-EUP', 'ZC2-AMP', 'ZC3-DEC', 'ZC4-CLG', 'ZC5-CSM']
groups = ['ZC1-EUP', 'ZC3-DEC', 'ZC4-CLG'] # just EUP, DEC, and CLG?
include_total = True
use_raw_ecospace = True  # if True, use raw Ecospace NetCDF, otherwise use pre-matched CSV
season_order = ['Winter', 'Spring', 'Summer', 'Fall']
log_offset = 1e-6
lag_years = 3

# === LOAD NPGO ===
df_npgo = pd.read_csv(f"{data_path}/{npgo_file}")
df_npgo.columns = ["date", "npgo"]
df_npgo["date"] = pd.to_datetime(df_npgo["date"])
df_npgo = df_npgo[df_npgo["npgo"] > -99]
df_npgo["year"] = df_npgo["date"].dt.year
df_npgo["month"] = df_npgo["date"].dt.month
spring_months = [3, 4, 5]
df_npgo_spring = df_npgo[df_npgo["month"].isin(spring_months)]
df_npgo_annual = df_npgo_spring.groupby("year")["npgo"].mean()
df_npgo_lagged = df_npgo_annual.shift(lag_years)

# === LOAD OBSERVATION DATA ===
df = pd.read_csv(f"{data_path}/{input_file}")
for col in groups:
    nonzero_vals = df[col][df[col] > 0]
    if not nonzero_vals.empty:
        min_nonzero = nonzero_vals.min()
        df[col] = df[col].apply(lambda x: np.random.uniform(0, 0.5 * min_nonzero) if x == 0 else x)

df_station = df.groupby(['Station', 'Year', 'Season'])[groups + [f"EWE-{g}" for g in groups]].mean().reset_index()
for col in groups:
    df_station[f'log10_{col}'] = np.log10(df_station[col] + log_offset)

results_obs = []
for group in groups:
    col = f'log10_{group}'
    df_temp = df_station.groupby('Year')[col].mean().dropna()
    df_temp = df_temp[(df_temp.index >= 1980) & (df_temp.index <= 2018)]
    clim = df_temp[df_temp.index.isin(range(1980, 2018))]
    anomaly = (df_temp - clim.mean()) / clim.std()
    results_obs.append(pd.DataFrame({'Year': df_temp.index, 'Group': group, 'Anomaly': anomaly.values, 'Source': 'Observed'}))

if include_total:
    df_station['TOT'] = df_station[groups].sum(axis=1)
    df_station['log10_TOT'] = np.log10(df_station['TOT'] + log_offset)
    df_temp = df_station.groupby('Year')['log10_TOT'].mean().dropna()
    df_temp = df_temp[(df_temp.index >= 1980) & (df_temp.index <= 2018)]
    clim = df_temp[df_temp.index.isin(range(1980, 2018))]
    anomaly = (df_temp - clim.mean()) / clim.std()
    results_obs.append(pd.DataFrame({'Year': df_temp.index, 'Group': 'TOTAL', 'Anomaly': anomaly.values, 'Source': 'Observed'}))

# === LOAD MODEL DATA FROM RAW ECOSPACE ===
ds = xr.open_dataset(raw_ecospace_nc)
results_mod = []
for group in groups:
    da = ds[group]
    da_mean = da.mean(dim=('row', 'col'), skipna=True)
    series = da_mean.groupby('time.year').mean().to_pandas()
    series = series[(series.index >= 1980) & (series.index <= 2018)]
    clim = series.loc[series.index.isin(range(1980, 2018))]
    anomaly = (series - clim.mean()) / clim.std()
    results_mod.append(pd.DataFrame({'Year': series.index, 'Group': group, 'Anomaly': anomaly.values, 'Source': 'Model'}))

if include_total:
    da_total = sum([ds[g] for g in groups])
    da_mean = da_total.mean(dim=('row', 'col'), skipna=True)
    series = da_mean.groupby('time.year').mean().to_pandas()
    series = series[(series.index >= 1980) & (series.index <= 2018)]
    clim = series.loc[series.index.isin(range(1980, 2018))]
    anomaly = (series - clim.mean()) / clim.std()
    results_mod.append(pd.DataFrame({'Year': series.index, 'Group': 'TOTAL', 'Anomaly': anomaly.values, 'Source': 'Model'}))

# === Combine and Plot ===
anom_df = pd.concat(results_obs + results_mod)
all_groups = groups + (['TOTAL'] if include_total else [])

fig, axs = plt.subplots(len(all_groups), 2, figsize=(14, 3.5 * len(all_groups)), sharex=True)
if len(all_groups) == 1:
    axs = [axs]

letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)']

for i, group in enumerate(all_groups):
    for j, source in enumerate(["Model", "Observed"]):
        sub = anom_df[(anom_df["Group"] == group) & (anom_df["Source"] == source)].set_index("Year")
        npgo = df_npgo_lagged[sub.index]
        colors = ['blue' if x > 0 else 'red' for x in sub['Anomaly']]
        axs[i][j].bar(sub.index, sub['Anomaly'], color=colors)
        axs[i][j].plot(npgo.index, npgo.values, color='black', label='NPGO (lagged)')
        axs[i][j].set_xlim(left=df_station['Year'].min())
        axs[i][j].axhline(0, color='grey', linestyle='--')
        axs[i][j].set_ylabel('Standardised Anomaly')
        axs[i][j].legend(loc='lower left')
        axs[i][j].text(0.05, 0.95, f'{letters[i]} {group} - {source}', transform=axs[i][j].transAxes, fontsize=12, verticalalignment='top')


plt.tight_layout()
plt.savefig(output_fig)
plt.show()
print(f"Saved panel plot to {output_fig}")

# === Plot Total Only as Stacked Panels if include_total ===
if include_total:
    total_obs = anom_df[(anom_df['Group'] == 'TOTAL') & (anom_df['Source'] == 'Observed')].set_index('Year')
    total_mod = anom_df[(anom_df['Group'] == 'TOTAL') & (anom_df['Source'] == 'Model')].set_index('Year')
    npgo = df_npgo_lagged[total_obs.index]

    fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    axs[0].bar(total_obs.index, total_obs['Anomaly'], color=['blue' if x > 0 else 'red' for x in total_obs['Anomaly']])
    axs[0].plot(npgo.index, npgo.values, color='black', label='NPGO (lagged)')
    axs[0].axhline(0, color='grey', linestyle='--')
    axs[0].set_ylabel('Observed Anomaly')
    axs[0].set_title('(a) Total Zooplankton - Observed')
    axs[0].legend(loc='lower left')

    axs[1].bar(total_mod.index, total_mod['Anomaly'], color=['blue' if x > 0 else 'red' for x in total_mod['Anomaly']])
    axs[1].plot(npgo.index, npgo.values, color='black', label='NPGO (lagged)')
    axs[1].axhline(0, color='grey', linestyle='--')
    axs[1].set_ylabel('Modelled Anomaly')
    axs[1].set_title('(b) Total Zooplankton - Model')
    axs[1].legend(loc='lower left')
    axs[1].set_xlabel('Year')

    for ax in axs:
        ax.set_xlim(left=df_station['Year'].min())

    plt.tight_layout()
    plt.savefig(total_fig)
    plt.show()
    print(f"Saved total panel plot to {total_fig}")