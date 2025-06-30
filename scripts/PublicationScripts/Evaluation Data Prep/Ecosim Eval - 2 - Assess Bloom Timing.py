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
"""

# ====== Imports ======
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from datetime import datetime, timedelta
from sklearn.metrics import mean_squared_error, mean_absolute_error
import os

# ====== Configuration ======
SCENARIO = "SC123"
ECOSIM_CSV_PATH = f"..//..//..//data/evaluation/ecosim_{SCENARIO}_onerun_B_dates_seasons.csv"
STATS_OUT_PATH = "..//..//..//data//evaluation//"
FIGS_OUT_PATH = "..//..//..//figs//"


# Bloom detection parameters
THRESHOLD_FACTOR = 1.05
SUB_THRESHOLD_FACTOR = 0.7
LOG_TRANSFORM = True
MEAN_OR_MEDIAN = "median"

# ====== Load observation datasets ======
def load_observation_bloom_dfs():
    doy_suchy = [100, 68, 50, 83, 115, 115, 100, 100, 92, 88, 92, 92, 55, 77]
    years_suchy = list(range(2003, 2017))
    suchy_df = pd.DataFrame({
        'Year': years_suchy,
        'Day of Year': doy_suchy
    })
    doy_allen = [94, 78, 81, 82, 87, 76, 75, 87, 99, 76,
                 81, 78, 77, 55, 86, 86, 67, 87, 77, 104,
                 66, 70, 92, 86, 81, 62, 88, 100, 90, 97,
                 104, 103, 98, 88, 86, 77, 88, 104, 77]
    years_allen = list(range(1980, 2019))
    allen_df = pd.DataFrame({
        'Year': years_allen,
        'Day of Year': doy_allen
    })

    return suchy_df, allen_df

# ====== Load Ecosim model outputs ======
def load_ecosim_dataset():
    df = pd.read_csv(ECOSIM_CSV_PATH)
    df['Date'] = pd.to_datetime(df['date'])
    df['Year'] = df['Date'].dt.year
    df['Day of Year'] = df['Date'].dt.dayofyear
    return df

# ====== Bloom detection ======
def find_bloom_doy(df, biomass_col='Biomass', threshold_factor=1.05):
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

# ====== Align years ======
def align_years(df_model, df_obs):
    return df_model[df_model['Year'].isin(df_obs['Year'])]

# ====== Plotting ======
def plot_bloom_comparison(df_model, df_obs, label_model="Ecosim", label_obs="Observation", filename="bloom_timing_ecosim.png"):
    df_merged = df_model.merge(df_obs, on="Year", suffixes=("_Model", "_Obs"))

    plt.figure(figsize=(10, 5))
    plt.plot(df_merged['Year'], df_merged['Day of Year_Model'], 'o-', label=label_model)
    plt.plot(df_merged['Year'], df_merged['Day of Year_Obs'], 's-', label=label_obs)

    plt.xlabel("Year")
    plt.ylabel("Day of Year")
    plt.title(f"Bloom Timing Comparison: {label_model} vs {label_obs}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIGS_OUT_PATH, filename))
    plt.show()

# ====== Evaluation ======
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
    df_stats.to_csv(os.path.join(STATS_OUT_PATH, f"ecosim_bloom_eval_stats_{scenario}.csv"), index=False)


# ====== Main script ======
def main():
    suchy_df, allen_df = load_observation_bloom_dfs()
    ecosim_df = load_ecosim_dataset()

    # Bloom detection
    bloom_df = find_bloom_doy(ecosim_df, biomass_col='17', threshold_factor=THRESHOLD_FACTOR)
    bloom_df.to_csv(os.path.join(STATS_OUT_PATH, f"ecosim_bloom_timing_{SCENARIO}.csv"), index=False)

    # Align model to observation years
    bloom_df_suchy_aligned = align_years(bloom_df, suchy_df)
    bloom_df_allen_aligned = align_years(bloom_df, allen_df)

    # Evaluation
    stats_suchy = evaluate_model(suchy_df['Day of Year'], bloom_df_suchy_aligned['Day of Year'])
    stats_allen = evaluate_model(allen_df['Day of Year'], bloom_df_allen_aligned['Day of Year'])

    print("Evaluation Statistics vs Suchy:", stats_suchy)
    print("Evaluation Statistics vs Allen:", stats_allen)

    # Plotting
    plot_bloom_comparison(bloom_df_suchy_aligned, suchy_df, label_model="Ecosim", label_obs="Suchy",
                          filename=f"ecosim_{SCENARIO}vs_satellite.png")
    plot_bloom_comparison(bloom_df_allen_aligned, allen_df, label_model="Ecosim", label_obs="Allen",
                          filename=f"ecosim_{SCENARIO}vs_C09.png")


if __name__ == "__main__":
    main()

