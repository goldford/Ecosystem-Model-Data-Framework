"""
Ecosim Bloom Timing Evaluation Script
Author: G Oldford, 2025

Description:
This script evaluates bloom timing using **Ecosim model outputs** (1D CSV time series).
It compares Ecosim-derived bloom timing against observational datasets:
- Satellite-derived bloom dates (Suchy et al. 2022)
- Allen 1D model outputs (Collins et al. 2009)

Outputs:
- CSV files summarising bloom timing
- Comparative plots
- Evaluation statistics (RMSE, MAE, Bias, R, Willmott skill)

Note: Designed to be independent from the Ecospace bloom evaluation script.
      Some shared or semi-shared scripts should be combined and put in helpers.
"""

# -------------------------------------------
# Imports
# -------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from datetime import datetime, timedelta
from sklearn.metrics import mean_squared_error, mean_absolute_error
import os
import ecosim_config # eval config file

# -------------------------------------------
# Configuration
# -------------------------------------------

SCENARIO = ecosim_config.SCENARIO
ECOSIM_CSV_PATH = ecosim_config.ECOSIM_F_PREPPED_SINGLERUN
STATS_OUT_PATH = ecosim_config.OUTPUT_DIR_EVAL
FIGS_OUT_PATH = ecosim_config.OUTPUT_DIR_FIGS

# User-specified biomass columns for Satellite and C09
BIOMASS_COLS_SATELLITE = ecosim_config.BIOMASS_COLS_SATELLITE
BIOMASS_COLS_C09 = ecosim_config.BIOMASS_COLS_C09

TOTAL_BIOMASS_COL_SATELLITE = ecosim_config.TOTAL_BIOMASS_COL_SATELLITE
TOTAL_BIOMASS_COL_C09 = ecosim_config.TOTAL_BIOMASS_COL_C09

# Bloom detection parameters
THRESHOLD_FACTOR = ecosim_config.THRESHOLD_FACTOR
SUB_THRESHOLD_FACTOR = ecosim_config.SUB_THRESHOLD_FACTOR
LOG_TRANSFORM = ecosim_config.LOG_TRANSFORM
MEAN_OR_MEDIAN = ecosim_config.MEAN_OR_MEDIAN


# -------------------------------------------
# Load observation datasets
# -------------------------------------------

def load_observation_bloom_dfs():
    doy_satellite = [100, 68, 50, 83, 115, 115, 100, 100, 92, 88, 92, 92, 55, 77]
    years_satellite = list(range(2003, 2017))
    satellite_df = pd.DataFrame({
        'Year': years_satellite,
        'Day of Year': doy_satellite
    })
    doy_C09 = [94, 78, 81, 82, 87, 76, 75, 87, 99, 76,
               81, 78, 77, 55, 86, 86, 67, 87, 77, 104,
               66, 70, 92, 86, 81, 62, 88, 100, 90, 97,
               104, 103, 98, 88, 86, 77, 88, 104, 77]
    years_C09 = list(range(1980, 2019))
    C09_df = pd.DataFrame({
        'Year': years_C09,
        'Day of Year': doy_C09
    })

    return satellite_df, C09_df


# -------------------------------------------
# Load Ecosim model outputs and compute total biomass
# -------------------------------------------

def load_ecosim_dataset():
    df = pd.read_csv(ECOSIM_CSV_PATH)
    df['Date'] = pd.to_datetime(df['date'])
    df['Year'] = df['Date'].dt.year
    df['Day of Year'] = df['Date'].dt.dayofyear

    # Sum specified biomass columns for Satellite and C09
    df[TOTAL_BIOMASS_COL_SATELLITE] = df[BIOMASS_COLS_SATELLITE].sum(axis=1, skipna=True)
    df[TOTAL_BIOMASS_COL_C09] = df[BIOMASS_COLS_C09].sum(axis=1, skipna=True)

    return df


# -------------------------------------------
# Bloom detection
# -------------------------------------------

def find_bloom_doy(df, biomass_col, threshold_factor=1.05):
    bloom_dates = []
    bloom_doys = []

    for year, group in df.groupby('Year'):
        ts = group.copy()
        if LOG_TRANSFORM:
            ts[biomass_col] = np.log(ts[biomass_col] + 0.01)

        base = ts[biomass_col].median() if MEAN_OR_MEDIAN == "median" else ts[biomass_col].mean()
        threshold = base * threshold_factor

        bloom = ts[ts[biomass_col] > threshold]
        if not bloom.empty:
            bloom_date = bloom.iloc[0]['Date']
            bloom_dates.append(bloom_date)
            bloom_doys.append(bloom_date.timetuple().tm_yday)
        else:
            bloom_dates.append(None)
            bloom_doys.append(None)

    return pd.DataFrame({
        'Year': df['Year'].unique(),
        'Bloom Date': bloom_dates,
        'Day of Year': bloom_doys
    })

#  Align years
def align_years(df_model, df_obs):
    return df_model[df_model['Year'].isin(df_obs['Year'])]

# ------------------------------------------
# -
#  Plotting
# -------------------------------------------
def plot_bloom_comparison(df_model, df_obs, label_model="Ecosim", label_obs="Observation", filename="bloom_timing_ecosim.png"):

    df_merged = df_model.merge(df_obs, on="Year", suffixes=("_Model", "_Obs"))

    plt.figure(figsize=(10, 5))

    # Plot model with error bars and line
    plt.errorbar(df_merged['Year'], df_merged['Day of Year_Model'], yerr=1.5,
                 fmt='o', color='black', label=label_model,
                 markersize=3, capsize=3)
    plt.plot(df_merged['Year'], df_merged['Day of Year_Model'], '-', color='black',
             markersize=0)

    # Plot observations with error bars and line
    plt.errorbar(df_merged['Year'], df_merged['Day of Year_Obs'], yerr=4,
                 fmt='s', color='blue', label=label_obs,
                 markersize=3, capsize=3)
    plt.plot(df_merged['Year'], df_merged['Day of Year_Obs'], '-', color='blue',
             markersize=0)

    plt.xlabel("Year")
    plt.ylabel("Day of Year")
    plt.title(f"Bloom Timing Comparison: {label_model} vs {label_obs}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(STATS_OUT_PATH, filename))
    plt.show()


# -------------------------------------------
# Evaluation
# -------------------------------------------

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
        return max(0, 1 - num / den)


def evaluate_model(obs, mod):
    obs = np.asarray(obs)
    mod = np.asarray(mod)
    valid = ~np.isnan(obs) & ~np.isnan(mod)

    obs_valid = obs[valid]
    mod_valid = mod[valid]

    mse = mean_squared_error(obs_valid, mod_valid)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(obs_valid, mod_valid)
    bias = np.mean(mod_valid - obs_valid)
    r = np.corrcoef(obs_valid, mod_valid)[0, 1] if len(obs_valid) > 1 else np.nan
    skill = willmott1981(obs_valid, mod_valid)
    std_obs = np.std(obs_valid)
    std_mod = np.std(mod_valid)

    return {
        "MSE": round(mse, 3),
        "RMSE": round(rmse, 3),
        "MAE": round(mae, 3),
        "Bias": round(bias, 3),
        "R": round(r, 3) if not np.isnan(r) else np.nan,
        "Willmott Skill": round(skill, 3) if not np.isnan(skill) else np.nan,
        "Obs StdDev": round(std_obs, 3),
        "Model StdDev": round(std_mod, 3)
    }


def export_evaluation_stats(stats_list, out_path, scenario):
    df_stats = pd.DataFrame(stats_list)
    outfile = os.path.join(STATS_OUT_PATH, f"ecosim_bloom_eval_stats_{scenario}.csv")
    df_stats.to_csv(outfile, index=False)
    print("Exported bloom timing stats to: " + outfile)


# -------------------------------------------
# Main script
# -------------------------------------------

def run_bloom_eval():
    satellite_df, C09_df = load_observation_bloom_dfs()
    ecosim_df = load_ecosim_dataset()

    # Bloom detection for Satellite biomass columns
    bloom_df_satellite = find_bloom_doy(ecosim_df, biomass_col=TOTAL_BIOMASS_COL_SATELLITE, threshold_factor=THRESHOLD_FACTOR)
    bloom_df_satellite.to_csv(os.path.join(STATS_OUT_PATH, f"ecosim_bloom_timing_satellite_{SCENARIO}.csv"), index=False)

    # Bloom detection for C09 biomass columns
    bloom_df_C09 = find_bloom_doy(ecosim_df, biomass_col=TOTAL_BIOMASS_COL_C09, threshold_factor=THRESHOLD_FACTOR)
    bloom_df_C09.to_csv(os.path.join(STATS_OUT_PATH, f"ecosim_bloom_timing_C09_{SCENARIO}.csv"), index=False)

    # Align model to observation years
    bloom_df_satellite_aligned = bloom_df_satellite[bloom_df_satellite['Year'].isin(satellite_df['Year'])]
    bloom_df_C09_aligned = bloom_df_C09[bloom_df_C09['Year'].isin(C09_df['Year'])]

    # Evaluation
    stats_suchy = evaluate_model(satellite_df['Day of Year'], bloom_df_satellite_aligned['Day of Year'])
    stats_allen = evaluate_model(C09_df['Day of Year'], bloom_df_C09_aligned['Day of Year'])
    export_evaluation_stats([stats_suchy, stats_allen], STATS_OUT_PATH, SCENARIO)

    print("Evaluation Statistics vs Suchy:", stats_suchy)
    print("Evaluation Statistics vs Allen:", stats_allen)

    # Plotting
    plot_bloom_comparison(bloom_df_satellite_aligned, satellite_df, label_model="Ecosim", label_obs="Suchy",
                          filename=f"ecosim_{SCENARIO}vs_satellite.png")
    plot_bloom_comparison(bloom_df_C09_aligned, C09_df, label_model="Ecosim", label_obs="Allen",
                          filename=f"ecosim_{SCENARIO}vs_C09.png")


if __name__ == "__main__":
    run_bloom_eval()

