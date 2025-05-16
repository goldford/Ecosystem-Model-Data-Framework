import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import pearsonr
from statsmodels.tsa.stattools import ccf


# === CONFIGURATION ===
SCENARIO = "SC88_1"
FILENAME = "Scv88_1-All_Groups_20250506_1978-2018.nc"
NC_PATH = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//"
NPGO_PATH = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//evaluation//npgo.csv"

target_groups = ['NK1-COH', 'NK2-CHI']
season_months = {
    "Spring": [3, 4, 5],
    "Summer": [6, 7, 8],
    "Fall": [9, 10, 11],
    "Spring_Summer_Fall": [3, 4, 5, 6, 7, 8, 9, 10, 11],
    "Year_Round": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
}
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

# Limit to spring and summer months
df_npgo_seasonal = df_npgo[df_npgo["month"].isin(season_months["Year_Round"])]
df_npgo_annual = df_npgo_seasonal.groupby("year")["npgo"].mean()
df_npgo_annual = df_npgo_annual[(df_npgo_annual.index >= analysis_start_year-2) & (df_npgo_annual.index <= analysis_end_year+2)]

# Lag NPGO by 2 years
df_npgo_annual_lagged = df_npgo_annual.shift(-1)

# === Utility ===
def safe_log(x):
    return xr.where(x > 0, np.log(x), np.nan)


# === Load Ecospace data ===
ds = xr.open_dataset(os.path.join(NC_PATH, FILENAME))
mask = ds['depth'] > 0
results = []

for group in target_groups:
    print(group)
    da = ds[group].where(mask)
    da_seasonal = da.sel(time=da['time.month'].isin(season_months["Spring_Summer_Fall"]))
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

for group in target_groups:
    print(group)
    plt.figure(figsize=(10, 4))
    sub = anomaly_df[anomaly_df['Group'] == group].set_index('Year')
    colors = ['blue' if x > 0 else 'red' for x in sub['Anomaly']]
    plt.bar(sub.index, sub['Anomaly'], color=colors, label='Model Anomaly')
    plt.plot(df_npgo_annual_lagged.index, df_npgo_annual_lagged.values, color='black', label='NPGO Index')
    plt.axhline(0, color='grey', linestyle='--')
    plt.title(f'{group} Spring+Summer Biomass Anomalies vs NPGO')
    plt.xlabel('Year')
    plt.ylabel('Anomaly (log scale)')
    plt.legend()
    plt.tight_layout()
    plt.show()

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