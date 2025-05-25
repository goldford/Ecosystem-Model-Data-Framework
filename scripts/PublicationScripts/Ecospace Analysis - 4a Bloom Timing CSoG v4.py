"""
Analyzing Bloom Timing in Ecospace Outputs
by: G Oldford, 2023-2025
----------------------------------------------------------
Compares bloom timing from Ecospace model outputs against:
1. Satellite-derived bloom dates (Suchy et al. 2022)
2. Allen 1D model outputs at S3 (Collins et al. 2009; Allen & Latournelle)
3. SSC 3D model outputs (SalishSeaCast v201905)

Outputs:
- CSV files of bloom timing
- Multiple comparative figures

Notes:
- to do: if recompute bloom timing is False, and user changes the Vars to analyse - does it work?
"""

# ====== Imports ======
import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path
from datetime import datetime, timedelta
from sklearn.metrics import mean_squared_error, mean_absolute_error


from helpers import (
    read_sdomains, find_nearest_point, find_bloom_doy,
    buildSortableString
)


# ================================
# Configuration
# ================================
ECOSPACE_OUT_PATH = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//"
DOMAIN_CONFIG_PATH = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
DOMAIN_FILE = "analysis_domains_suchy.yml"
SAT_MASK_PATH = os.path.join(DOMAIN_CONFIG_PATH, '..//..//data//evaluation//suchy_ecospace_mask.nc')

RECOMPUTE_BLOOM_TIMING_SAT = True  # Set to True to force recomputation as needed, saves time
RECOMPUTE_BLOOM_TIMING_C09 = True

SCENARIO = 'FULLKEY_SC88_2'
ECOSPACE_CODE = "Scv88_2-All_Groups_20250506"
FILENM_STRT_YR = 1978
FILENM_END_YR = 2018
START_YEAR = 1980 # analysis years (exclude spinup?)
END_YEAR = 2018

# if multiple vars listed here, it will sum across them when computing anomalies, bloom timing etc
#OK with run 96 this is first time model fit has been okay with all pp groups
VARIABLES_TO_ANALYZE_SAT = ["PP1-DIA", "PP2-NAN", "PP3-PIC"]
VARIABLES_TO_ANALYZE_C09 = ["PP1-DIA"]

ANNUAL_AVG_METHOD_SAT = "all" # should average bloom compared against be from all years, or just one year
ANNUAL_AVG_METHOD_C09 = "annual" # annual or all

EXCLUDE_DEC_JAN_SAT = True
EXCLUDE_DEC_JAN_C09 = False

# Bloom detection method
LOG_TRANSFORM = True
MEAN_OR_MEDIAN = "median"
THRESHOLD_FACTOR = 1.05
SUB_THRESHOLD_FACTOR = 0.7

MIN_Y_TICK = 38
CREATE_MASKS = False  # Set to True to regenerate masks
DO_NUTRIENTS = False


# ================================
# Utility Functions
# ================================

def convert_doy_to_date(doys, start_year):
    return [datetime(start_year + i, 1, 1) + timedelta(days=doy - 1) for i, doy in enumerate(doys)]

def create_bloom_df(years, doys, label="avg"):
    dates = convert_doy_to_date(doys, years[0])
    return pd.DataFrame({
        'Year': years,
        'Bloom Date': dates,
        'Day of Year': doys,
        'Bloom Early Late': label if isinstance(label, list) else [label]*len(years)
    })

def classify_bloom_by_std(doys, mean_val, std_val):
    results = []
    for doy in doys:
        if doy <= (mean_val - std_val):
            results.append("early")
        elif doy >= (mean_val + std_val):
            results.append("late")
        elif doy <= (mean_val + std_val - 1) and doy >= (mean_val - std_val + 1):
            results.append("avg")
        else:
            results.append("cusp")
    # exploring 'cusp' method
    # for doy in doys:
    #     if doy + 4 <= (mean_val - std_val):
    #         results.append("early")
    #     elif doy - 4 >= (mean_val + std_val):
    #         results.append("late")
    #     elif doy + 4 <= (mean_val + std_val - 1) and doy - 4 >= (mean_val - std_val + 1):
    #         results.append("avg")
    #     else:
    #         results.append("cusp")
    return results

# ================================
# Mask Generation
# ================================

def generate_2d_mask(ds, domain_file_path, region_name='SGC2', depth_var='depth'):
    lat = ds['lat'].values
    lon = ds['lon'].values
    dep = ds[depth_var].values

    sdomains = read_sdomains(domain_file_path)
    polygon_path = Path(sdomains[region_name])
    points = np.vstack((lat.flatten(), lon.flatten())).T
    region_mask = polygon_path.contains_points(points).reshape(lat.shape)
    depth_mask = dep > 0
    combined_mask = region_mask & depth_mask

    mask_ds = xr.Dataset({
        'mask': (('row', 'col'), combined_mask)
    }, coords={'lat': (('row', 'col'), lat), 'lon': (('row', 'col'), lon)})

    out_path = os.path.join(DOMAIN_CONFIG_PATH, '..//..//data//evaluation//suchy_ecospace_mask.nc')
    mask_ds.to_netcdf(out_path)
    print(f"Saved 2D mask to {out_path}")
    return mask_ds

def generate_1d_mask(ds, lat_target, lon_target, dx=0.01):
    lats = ds['lat'].values
    lons = ds['lon'].values
    deps = ds['depth'].values
    valid_mask = deps != 0
    return find_nearest_point(lon_target, lat_target, lons, lats, valid_mask, dx)

# ================================
# Bloom Observation Datasets
# ================================

def load_observation_bloom_dfs():
    # Suchy satellite data
    doy_suchy = [100, 68, 50, 83, 115, 115, 100, 100, 92, 88, 92, 92, 55, 77]
    years_suchy = list(range(2003, 2017))
    suchy_df = create_bloom_df(years_suchy, doy_suchy, ["avg", "avg", "early", "avg",
                                                        "late", "late", "avg", "avg",
                                                        "avg", "avg", "avg", "avg",
                                                        "early", "avg"])

    # Gower satellite data
    doy_gower = [83, 68, 76, 85, 71, 56, 84, 49, 71, 58, 90, 97, 91, 69, 94, 66, 86]
    years_gower = list(range(2000, 2017))
    gower_df = create_bloom_df(years_gower, doy_gower, "na")
    gower_df.columns = [col + "_Gower" if col != 'Year' else 'Year_Gower' for col in gower_df.columns]

    # Allen 1D model data
    doy_allen = [94, 78, 81, 82, 87, 76, 75, 87, 99, 76,
                 81, 78, 77, 55, 86, 86, 67, 87, 77, 104,
                 66, 70, 92, 86, 81, 62, 88, 100, 90, 97,
                 104, 103, 98, 88, 86, 77, 88, 104, 77]
    years_allen = list(range(1980, 2019))
    mean_allen = np.mean(doy_allen)
    std_allen = np.std(doy_allen)
    labels_allen = classify_bloom_by_std(doy_allen, mean_allen, std_allen)
    allen_df = create_bloom_df(years_allen, doy_allen, labels_allen)
    allen_df.columns = [col + "_Allen" if col != 'Year' else 'Year_Allen' for col in allen_df.columns]

    return suchy_df, gower_df, allen_df

# ================================
# Model Bloom Detection
# ================================

def load_ecospace_dataset():
    fname = f"{ECOSPACE_CODE}_{FILENM_STRT_YR}-{FILENM_END_YR}.nc"
    path = os.path.join(ECOSPACE_OUT_PATH, fname)
    return xr.open_dataset(path)

def compute_bloom_timing(ds, var_name, mask=None, row=None, col=None,
                         bloom_early=68, bloom_late=108,
                         yr_strt=1980, yr_end=2018,
                         exclude_dec_jan=False,
                         avg_method="annual"):

    years = range(yr_strt, yr_end + 1)
    values, timestamps = [], []
    print("computing bloom timing")

    for year in years:
        print(year)
        yearly_ds = ds.sel(time=str(year))

        if mask is not None:
            yearly_ds = yearly_ds.where(mask)
        elif row is not None and col is not None:
            yearly_ds = yearly_ds.isel(row=row, col=col)
        else:
            raise ValueError("Must provide either a 2D mask or specific row/col indices.")

        if isinstance(var_name, list):
            total = None
            for g in var_name:
                if g in yearly_ds:
                    total = yearly_ds[g] if total is None else total + yearly_ds[g]
                else:
                    print(f"Warning: {g} not found in dataset for year {year}")
            yearly_var = total
        else:
            yearly_var = yearly_ds[var_name]

        if LOG_TRANSFORM:
            yearly_var = np.log(yearly_var + 0.01)
        else:
            yearly_var = yearly_ds[var_name]
            if LOG_TRANSFORM:
                yearly_var = np.log(yearly_var + 0.01)

        for ts in yearly_ds.time:
            val = yearly_var.sel(time=ts)
            value = np.nanmedian(val) if MEAN_OR_MEDIAN == "median" else np.nanmean(val)
            values.append(value)
            timestamps.append(pd.Timestamp(ts.values))

    df = pd.DataFrame({
        'Year': pd.to_datetime(timestamps).year,
        'Date': pd.to_datetime(timestamps),
        'Value': values
    })

    bloom_dates, bloom_doys, bloom_categories = find_bloom_doy(
        df, 'Value',
        thrshld_fctr=THRESHOLD_FACTOR,
        sub_thrshld_fctr=SUB_THRESHOLD_FACTOR,
        average_from=avg_method,
        mean_or_median=MEAN_OR_MEDIAN,
        exclude_juntojan=exclude_dec_jan,
        bloom_early=bloom_early,
        bloom_late=bloom_late
    )

    return pd.DataFrame({
        'Year': years,
        'Bloom Date': bloom_dates,
        'Day of Year': bloom_doys,
        'Bloom Early Late': bloom_categories
    })


# ================================
# Plotting Functions
# ================================
def standardize_columns(df):
    for col in df.columns:
        if "Year" in col and col != "Year" and col != "Day of Year":
            df = df.rename(columns={col: "Year"})
            break  # Stop after the first match
    for col in df.columns:
        if "Day of Year" in col and col != "Year":
            df = df.rename(columns={col: "Day of Year_Obs"})
            break  # Stop after the first match
    return df

def plot_bloom_comparison(df_model, df_obs, label_model="Ecospace",
                          label_obs="Suchy", title="Bloom Timing Comparison",
                          filename="bloom_timing_plot.png"):

    df_obs = standardize_columns(df_obs)

    df_merged = df_model.merge(df_obs, on="Year", suffixes=("_Model", "_Obs"))

    if label_obs == 'Satellite':
        fig_w = 7; fig_h = 4
    else:
        fig_w = 10; fig_h = 5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    for col in df_merged.columns:
        if "Day of Year" in col and col != "Day of Year_Obs":
            df_merged = df_merged.rename(columns={col: "Day of Year_Model"})
            break  # Stop after the first match
    ax.errorbar(df_merged['Year'], df_merged['Day of Year_Model'],
                yerr=1.5, fmt='o', markersize=4,
                label=label_model, color='black')
    ax.errorbar(df_merged['Year'], df_merged['Day of Year_Obs'],
                yerr=4, fmt='s', markersize=4,
                label=label_obs, color='blue')

    mean = np.nanmean(df_merged['Day of Year_Obs'])
    std = np.nanstd(df_merged['Day of Year_Obs'])

    ax.axhline(y=mean+std, linestyle='--', color='grey', label='') # late threshold
    ax.axhline(y=mean, linestyle='-', color='grey', label='')      # average
    ax.axhline(y=mean-std, linestyle='--', color='grey', label='') # early threshold

    # ax.fill_between([x_min-0.5, x_max+0.5], 68, 108, color='lightgrey', alpha=0.3)

    ax.set_xlabel("Year")
    ax.set_ylabel("Day of Year")
    ax.set_title(title)
    ax.set_ylim([MIN_Y_TICK, df_merged[["Day of Year_Model", "Day of Year_Obs"]].max().max() + 10])
    ax.legend()
    plt.tight_layout()
    plt.savefig('..//..//figs//' + filename)
    plt.show()
    plt.close()

# ================================
# Evaluation Statistics
# ================================

def willmott1981(obs, mod):
    mod = np.asarray(mod)
    obs = np.asarray(obs)
    num = np.nansum((mod - obs) ** 2)
    obs_mean = np.nanmean(obs)
    dM = np.abs(mod - obs_mean)
    dO = np.abs(obs - obs_mean)
    den = np.nansum((dM + dO) ** 2)
    if den == 0:
        return np.nan
    else:
        return np.max([0, 1 - num / den])

def evaluate_model(obs, mod, label=""):
    obs = np.asarray(obs)
    mod = np.asarray(mod)
    valid = ~np.isnan(obs) & ~np.isnan(mod)
    obs_valid = obs[valid]
    mod_valid = mod[valid]

    rmse = np.sqrt(mean_squared_error(obs_valid, mod_valid))
    mse = mean_squared_error(obs_valid, mod_valid)
    mae = mean_absolute_error(obs_valid, mod_valid)
    bias = np.mean(mod_valid - obs_valid)
    r = np.corrcoef(obs_valid, mod_valid)[0, 1] if len(obs_valid) > 1 else np.nan
    skill = willmott1981(obs_valid, mod_valid)
    std_obs = np.std(obs_valid)
    std_mod = np.std(mod_valid)

    return {
        "Label": label,
        "RMSE": rmse,
        "MAE": mae,
        "MSE": mse,
        "Bias": bias,
        "R": r,
        "Willmott Skill": skill,
        "Obs StdDev": std_obs,
        "Model StdDev": std_mod
    }

def evaluate_bloom_categories(df_obs, df_mod, col_obs='Bloom Early Late', col_mod='Bloom Early Late'):
    df_obs = standardize_columns(df_obs)

    for col in df_obs.columns:
        if "Bloom Early Late" in col and col != "Bloom Early Late":
            df_obs = df_obs.rename(columns={col: "Bloom Early Late"})
            break  # Stop after the first match

    df = df_obs.merge(df_mod, on='Year', suffixes=('_obs', '_mod'))
    agree_count = (df[f'{col_obs}_obs'] == df[f'{col_mod}_mod']).sum()
    total = len(df)
    return agree_count, total

def evaluate_overlap_by_timing(df_obs, df_mod, obs_col='Day of Year', mod_col='Day of Year'):

    df_obs = standardize_columns(df_obs)
    for col in df_obs.columns: #standardise
        if "Day of Year" in col and col != "Day of Year":
            df_obs = df_obs.rename(columns={col: "Day of Year"})
            break  # Stop after the first match


    df = df_obs.merge(df_mod, on='Year', suffixes=('_obs', '_mod'))
    # Define bounds
    obs_low = df[f'{obs_col}_obs'] - 4
    obs_high = df[f'{obs_col}_obs'] + 4
    mod_low = df[f'{mod_col}_mod'] - 1.5
    mod_high = df[f'{mod_col}_mod'] + 1.5
    overlap = (mod_high >= obs_low) & (mod_low <= obs_high)
    return overlap.sum(), len(overlap)

# ================================
# Nutrient Proxy from Biomass
# ================================

def compute_nutrient_concentration(ds, vars_list, year=2005, include_only=None, exclude=None, mask=None):
    # Constants
    N_free_init = 41.0   # g N m^-2 (96%)
    N_bound_init = 1.54  # g N m^-2 (4%)
    total_N = N_free_init + N_bound_init

    # Molar mass ratios (Redfield C:N is 106:16)
    C_to_N_ratio = 106 / 16

    # Apply inclusion/exclusion filters
    if include_only:
        vars_list = [var for var in vars_list if var in include_only]
    if exclude:
        vars_list = [var for var in vars_list if var not in exclude]

    # Subset one year
    year_ds = ds.sel(time=str(year))

    # Compute total biomass over time (assumed carbon units)
    total_biomass = None
    for var in vars_list:
        if var in year_ds:
            if total_biomass is None:
                total_biomass = year_ds[var]
            else:
                total_biomass += year_ds[var]

    if total_biomass is None:
        raise ValueError("None of the specified biomass variables found in dataset.")

    # Apply spatial mask if provided
    if mask is not None:
        total_biomass = total_biomass.where(mask)

    # Mean over space (depth integrated already assumed)
    biomass_time_series = total_biomass.mean(dim=["row", "col"], skipna=True)

    # Convert carbon to nitrogen (g C to g N)
    N_bound = biomass_time_series / C_to_N_ratio

    # Estimate remaining free N = total_N - N_bound
    N_free = total_N - N_bound
    N_free = N_free.where(N_free >= 0, 0)  # avoid negative nutrient concentrations

    # Prepare DataFrame
    df = pd.DataFrame({
        "Date": pd.to_datetime(year_ds.time.values),
        "Day of Year": pd.to_datetime(year_ds.time.values).dayofyear,
        "Year": year,
        "Total Biomass (C)": biomass_time_series.values,
        "N Bound (g m^-2)": N_bound.values,
        "N Free (g m^-2)": N_free.values
    })
    return df

def plot_nutrient_concentration(df_nutrient_all, ds=None, include_only=None, mask=None):
    plt.figure(figsize=(10, 5))

    # Climatological stats
    grouped = df_nutrient_all.groupby("Day of Year")
    mean_nutrient = grouped["N Free (g m^-2)"].mean()
    q10 = grouped["N Free (g m^-2)"].quantile(0.10)
    q20 = grouped["N Free (g m^-2)"].quantile(0.20)
    q30 = grouped["N Free (g m^-2)"].quantile(0.30)
    q40 = grouped["N Free (g m^-2)"].quantile(0.40)
    q60 = grouped["N Free (g m^-2)"].quantile(0.60)
    q70 = grouped["N Free (g m^-2)"].quantile(0.70)
    q80 = grouped["N Free (g m^-2)"].quantile(0.80)
    q90 = grouped["N Free (g m^-2)"].quantile(0.90)

    # Shaded quantile bands
    plt.fill_between(mean_nutrient.index, q10, q90, color="0.6", alpha=0.2, label="")
    plt.fill_between(mean_nutrient.index, q20, q80, color="0.5", alpha=0.2, label="")
    plt.fill_between(mean_nutrient.index, q30, q70, color="0.4", alpha=0.1, label="")
    plt.fill_between(mean_nutrient.index, q40, q60, color="0.4", alpha=0.1, label="")
    plt.plot(mean_nutrient.index, mean_nutrient, color="0.3", linewidth=2, label="Nutrients - Clim. Avg.")

    plt.xlabel("Day of Year")
    plt.ylabel("Nutrient Concentration (g N m$^{-2}$)")
    plt.title("Climatological Nutrient Concentration vs Total Biomass (All Years) -" + SCENARIO)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    # Add secondary axis for biomass
    if ds is not None and include_only is not None:
        ax = plt.gca()
        ax2 = ax.twinx()
        colors = plt.cm.tab10.colors

        for i, var in enumerate(include_only):
            biomasses = []
            for yr in df_nutrient_all['Year'].unique():
                try:
                    da = ds[var].sel(time=str(yr))
                    if mask is not None:
                        da = da.where(mask)
                    ts = da.mean(dim=["row", "col"], skipna=True)
                    biomasses.append(pd.DataFrame({
                        'Day of Year': pd.to_datetime(da['time'].values).dayofyear,
                        'Biomass': ts.values,
                        'Year': yr
                    }))
                except KeyError:
                    continue

            df_biomass = pd.concat(biomasses, ignore_index=True)
            grouped_biomass = df_biomass.groupby("Day of Year")
            mean_b = grouped_biomass['Biomass'].mean()
            q10_b = grouped_biomass['Biomass'].quantile(0.10)
            q90_b = grouped_biomass['Biomass'].quantile(0.90)

            ax2.fill_between(mean_b.index, q10_b, q90_b, color=colors[i % len(colors)], alpha=0.15)
            ax2.plot(mean_b.index, mean_b, color=colors[i % len(colors)], label=f"{var} Biomass")

        ax2.set_ylabel("Biomass (g C m$^{-2}$)")
        ax2.legend(loc="upper right")

    plt.show()
    plt.savefig('..//..//figs//' + "nutrient_concentration_climatology_" + SCENARIO + ".png")
    plt.close()

# ================================
# Entry Point
# ================================

def main():
    # Load observation datasets
    suchy_df, gower_df, allen_df = load_observation_bloom_dfs()

    # Load Ecospace model outputs
    ds = load_ecospace_dataset()

    # Generate masks if required
    if CREATE_MASKS:
        generate_2d_mask(ds, os.path.join(DOMAIN_CONFIG_PATH, DOMAIN_FILE))
        lat_1d, lon_1d = 49.1887166, -123.4966697
        idx_1d = generate_1d_mask(ds, lat_1d, lon_1d)
        print("1D index (row, col, depth):", idx_1d)

    mask_ds = xr.open_dataset(SAT_MASK_PATH)

    bloom_csv_path_suchy = f"..//..//data//evaluation//ecospace_bloom_timing_SSoG_{SCENARIO}.csv"
    bloom_csv_path_allen = f"..//..//data//evaluation//ecospace_bloom_timing_CSoG_{SCENARIO}.csv"

    var_name_C09 = VARIABLES_TO_ANALYZE_C09

    if RECOMPUTE_BLOOM_TIMING_C09:
        # Hardcoded Allen 1D location
        row_allen, col_allen = 52, 100

        bloom_df_allen = compute_bloom_timing(
            ds, var_name_C09, mask=None, row=col_allen, col=row_allen,
            bloom_early=allen_df['Day of Year_Allen'].mean() - allen_df['Day of Year_Allen'].std(),
            bloom_late=allen_df['Day of Year_Allen'].mean() + allen_df['Day of Year_Allen'].std(),
            yr_strt=1980, yr_end=2018,
            exclude_dec_jan=EXCLUDE_DEC_JAN_C09,
            avg_method=ANNUAL_AVG_METHOD_C09
        )
        bloom_df_allen.to_csv(bloom_csv_path_allen, index=False)
        print(f"Saved computed bloom timing for Allen 1D to {bloom_csv_path_allen}")

    else:
        bloom_df_allen = pd.read_csv(bloom_csv_path_allen)
        print(f"Loaded bloom timing from cached files: {bloom_csv_path_allen}")

    var_name_SAT = VARIABLES_TO_ANALYZE_SAT

    if RECOMPUTE_BLOOM_TIMING_SAT:
        bloom_df_suchy = compute_bloom_timing(
            ds, var_name_SAT, mask=mask_ds['mask'],
            bloom_early=68, bloom_late=108,
            yr_strt=2003, yr_end=2016,
            exclude_dec_jan=EXCLUDE_DEC_JAN_SAT,
            avg_method=ANNUAL_AVG_METHOD_SAT
        )
        bloom_df_suchy.to_csv(bloom_csv_path_suchy, index=False)
        print(f"Saved computed bloom timing for SSoG (Suchy) to {bloom_csv_path_suchy}")

    else:
        bloom_df_suchy = pd.read_csv(bloom_csv_path_suchy)
        print(f"Loaded bloom timing from cached files: {bloom_csv_path_suchy}")



    print('bloom doy mean, C09:')
    print(allen_df['Day of Year_Allen'].mean())
    print('bloom doy early and late thresholds, CO9:')
    print(allen_df['Day of Year_Allen'].mean() - allen_df['Day of Year_Allen'].std())
    print(allen_df['Day of Year_Allen'].mean() + allen_df['Day of Year_Allen'].std())

    # Print summaries
    print("Loaded bloom timing observation data:")
    print("Satellite (Suchy et al):", suchy_df.shape)
    print("Satellite (Gower et al):", gower_df.shape)
    print("C09:", allen_df.shape)
    print("Ecospace (Satellite):", bloom_df_suchy.shape)
    print("Ecospace (C09):", bloom_df_allen.shape)

    # Plot comparison with Allen 1D model
    plot_bloom_comparison(
        bloom_df_allen, allen_df,
        label_model="Ecospace Model", label_obs="C09 1D Model",
        title=f"Diatom Bloom Timing: {SCENARIO} vs C09",
        filename=f"ecospace_vs_C09_{SCENARIO}.png")

    # Plot comparison with Suchy data
    plot_bloom_comparison(
        bloom_df_suchy, suchy_df,
        label_model="Ecospace Model", label_obs="Satellite",
        title=f"Diatom Bloom Timing: {SCENARIO} vs Satellite",
        filename=f"ecospace_vs_satell_{SCENARIO}.png")

    # Compute and print statistics
    stats_suchy = evaluate_model(
        suchy_df["Day of Year"], bloom_df_suchy["Day of Year"], label="Suchy")
    stats_allen = evaluate_model(
        allen_df["Day of Year_Allen"], bloom_df_allen["Day of Year"], label="Allen")

    print("\nEvaluation Statistics:")
    for stat in [stats_suchy, stats_allen]:
        print(
            f"{stat['Label']}: R = {stat['R']:.3f}, "
            f"RMSE = {stat['RMSE']:.2f}, "
            f"MAE = {stat['MAE']:.2f}, "
            f"Bias = {stat['Bias']:.2f}, "
            f"Obs σ = {stat['Obs StdDev']:.2f}, "
            f"Model σ = {stat['Model StdDev']:.2f}, "
            f"Willmott = {stat['Willmott Skill']:.3f}")

    # Additional categorical comparisons
    agree_cat_suchy, total_cat_suchy = evaluate_bloom_categories(suchy_df, bloom_df_suchy)
    agree_cat_allen, total_cat_allen = evaluate_bloom_categories(allen_df, bloom_df_allen,
                                                                 col_obs='Bloom Early Late')
    print("\nCategorical Agreement:")
    print(f"Suchy: {agree_cat_suchy}/{total_cat_suchy} years agree in category")
    print(f"Allen: {agree_cat_allen}/{total_cat_allen} years agree in category")

    # Overlap by timing comparison
    overlap_suchy, n_suchy = evaluate_overlap_by_timing(suchy_df, bloom_df_suchy)
    overlap_allen, n_allen = evaluate_overlap_by_timing(allen_df, bloom_df_allen,
                                                        obs_col='Day of Year')
    if DO_NUTRIENTS:
        print("\nOverlap of Bloom Timing Windows:")
        print(f"Suchy: {overlap_suchy}/{n_suchy} years overlap within timing window")
        print(f"Allen: {overlap_allen}/{n_allen} years overlap within timing window")

        # Example: Nutrient diagnostic for one year
        all_biomass_vars = [
            "NK1-COH", "NK2-CHI", "NK3-FOR",
            "ZF1-ICT", "ZC1-EUP", "ZC2-AMP",
            "ZC3-DEC", "ZC4-CLG", "ZC5-CSM",
            "ZS1-JEL", "ZS2-CTH", "ZS3-CHA",
            "ZS4-LAR", "PZ1-CIL","PZ2-DIN",
            "PZ3-HNF","PP1-DIA", "PP2-NAN",
            "PP3-PIC", "BAC-BA1"
        ]

        include_only = ["PP1-DIA", "PP2-NAN", "PP3-PIC"]  # example filter
        years_to_plot = list(range(1980, 2019))

        nutrient_dfs = []
        for yr in years_to_plot:
            print(yr)
            df = compute_nutrient_concentration(ds, all_biomass_vars, year=yr, include_only=include_only,
                                                mask=mask_ds['mask'])
            nutrient_dfs.append(df)

        df_nutrient_all = pd.concat(nutrient_dfs, ignore_index=True)
        plot_nutrient_concentration(df_nutrient_all, ds=ds, include_only=include_only, mask=mask_ds['mask'])
        print("Saved nutrient climatology plot for 1980–2018.")

if __name__ == "__main__":
    main()




