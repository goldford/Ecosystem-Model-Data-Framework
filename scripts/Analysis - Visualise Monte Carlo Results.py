# Created Mar 3, 2025 by G Oldford
# Purpose: Review and visualise results from Ecosim Monte Carlo

# Data In:
#  1/ monte carlo outputs from ecosim
#  2/ reference time series
#       eg,
#       smolt-to-age 2 M from CWT
#       SoG M from Nelson et al 2024
#
# Output:
#  1/ plots

# C:\Users\Greig\Documents\EwE output\GeorgiaStrait2023_v97_v6_7_0_18858_64b\mc_EcosimScen97\
# mc_output_trial0006\mortality_monthly.csv
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, mean_absolute_error

# Define scenario directories
scenarios = {
    "SC104 - NO PP Anom": "SC104",
    "SC105 - WITH PP Anom": "SC105"
}

# Define paths
root_dir_mc = "C://Users//Greig//Documents//EwE output//"
root_dir_cwt_chin = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//1. Salmon//Chinook Survival Freshwater 2022//MODIFIED//"
root_dir_cwt_coho = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//1. Salmon//All Species Survival Exploitation//"
root_dir_nelson = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//1. Salmon//Coho_Mort_Seal_Nelson_2024\ORIGINAL//"

mc_root_base = {
    "SC104": os.path.join(root_dir_mc, "GeorgiaStrait2023_v103_v6_7_0_18858_64b//mc_EcosimScen104"),
    "SC105": os.path.join(root_dir_mc, "GeorgiaStrait2023_v103_v6_7_0_18858_64b//mc_EcosimScen105")
}

# Define column groups for each species/subtype
column_groups = {
    "Chinook H": ["21", "22", "23"],
    "Chinook WO": ["27", "28", "29"],
    "Chinook WS": ["34"],
    "Coho H": ["39"],
    "Coho W": ["43"]
}


# Load reference datasets
def load_reference_data():

    refs = {}

    # Chinook WO
    f = "salish_M_ocntp_sog_1971-2017.csv"
    df = pd.read_csv(os.path.join(root_dir_cwt_chin, f))
    df['year'] = df['brood_year'] + 1
    refs["Chinook WO"] = df[['year', 'mean_M_oc_sog']].rename(columns={'mean_M_oc_sog': 'obs'})
    refs["Chinook H"] = df[['year', 'mean_M_oc_sog']].rename(columns={'mean_M_oc_sog': 'obs'})


    # Chinook WS
    f = "salish_M_strmtp_sog_1971-2017.csv"
    df = pd.read_csv(os.path.join(root_dir_cwt_chin, f))
    df['year'] = df['brood_year'] + 1
    refs["Chinook WS"] = df[['year', 'mean_M_st_sog']].rename(columns={'mean_M_st_sog': 'obs'})
    refs["Chinook H"] = df[['year', 'mean_M_st_sog']].rename(columns={'mean_M_st_sog': 'obs'})


    # Coho
    f = "coho_M_EWE_TS Jan 2025.csv"
    df = pd.read_csv(os.path.join(root_dir_cwt_coho, f))
    df['year'] = df['year'] + 1
    refs["Coho W"] = df[['year', 'AVERAGE']].rename(columns={'AVERAGE': 'obs'})
    refs["Coho H"] = df[['year', 'AVERAGE']].rename(columns={'AVERAGE': 'obs'})

    return refs


# Load model outputs for all scenarios
def load_mc_data():
    all_data = {scen: {} for scen in scenarios.values()}
    for scen_label, scen_key in scenarios.items():
        mc_root = mc_root_base[scen_key]
        folders = [f for f in os.listdir(mc_root) if f.startswith("mc_output_trial")]
        for group, cols in column_groups.items():
            yearly_data = {}
            for folder in folders:
                path = os.path.join(mc_root, folder, "mortality_annual.csv")
                if os.path.exists(path):
                    df = pd.read_csv(path)
                    if "year\\group" in df.columns:
                        df["year\\group"] = df["year\\group"].astype(int)
                        df[group] = df[cols].sum(axis=1)
                        for year, val in zip(df["year\\group"], df[group]):
                            yearly_data.setdefault(year, []).append(val)
            all_data[scen_key][group] = yearly_data
    return all_data


# Plot comparison
def plot_comparison(all_data, refs):
    for group in column_groups.keys():
        plt.figure(figsize=(12, 6))

        # Plot reference
        ref_df = refs[group]
        plt.scatter(ref_df['year'], ref_df['obs'], color='black', label='Observed (CWT)', zorder=3)

        # Calculate and store fit statistics
        stats = []

        # Plot each scenario
        for label, scen_key in scenarios.items():
            years = sorted(all_data[scen_key][group].keys())
            matrix = np.array([all_data[scen_key][group][y] for y in years])
            percentiles = np.percentile(matrix, [5, 25, 50, 75, 95], axis=1)
            median = percentiles[2]

            plt.fill_between(years, percentiles[0], percentiles[4], alpha=0.2, label=f"{label} 5-95%")
            plt.plot(years, median, label=f"{label} median")

            # Calculate fit stats for overlapping years
            ref_merge = ref_df[ref_df['year'].isin(years)]
            model_vals = [np.nanmedian(all_data[scen_key][group][y]) for y in ref_merge['year']]
            # rmse = np.sqrt(mean_squared_error(ref_merge['obs'], model_vals))
            # mae = mean_absolute_error(ref_merge['obs'], model_vals)
            # stats.append(f"{label}: RMSE={rmse:.2f}, MAE={mae:.2f}")

        plt.title(f"Mortality Comparison - {group}")
        plt.xlabel("Year")
        plt.ylabel("Instantaneous Mortality Rate (M)")
        plt.legend()
        plt.grid(True)
        for i, stat in enumerate(stats):
            plt.text(0.01, 0.95 - 0.05*i, stat, transform=plt.gca().transAxes, fontsize=9, verticalalignment='top')

        plt.xlim(1970, 2018)
        plt.tight_layout()
        plt.show()


# Run all steps
reference_data = load_reference_data()
monte_carlo_data = load_mc_data()
plot_comparison(monte_carlo_data, reference_data)
