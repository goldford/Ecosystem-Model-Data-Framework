# Created Mar 3, 2025 by G Oldford
# Purpose: Review and visualise results from Ecosim Monte Carlo
# C:\Users\Greig\Documents\EwE output\GeorgiaStrait2023_v97_v6_7_0_18858_64b\mc_EcosimScen97\
# mc_output_trial0006\mortality_monthly.csv

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# need a purge script since MC generates so much shite
# purge everything without 'annual' in name
# purge everything with 'discards' in name


# mortality time series (for ref, fitting)
# Chinook - ocean and stream type, Freshwater et al., 2022
mort_path = "C:\\Users\\Greig\\Sync\\6. SSMSP Model\\Model Greig\\Data\\1. Salmon\\Chinook Survival Freshwater 2022\\MODIFIED\\"

# OCEAN TYPE
mort_f = "salish_M_ocntp_sog_1971-2017.csv"
m_p = os.path.join(mort_path, mort_f)
print(m_p)
if os.path.exists(m_p):
    print("found")
    df_mort_ref_chin_wo_sogfr = pd.read_csv(m_p)
    # since it's brood year, add one year to M
    df_mort_ref_chin_wo_sogfr['year'] = df_mort_ref_chin_wo_sogfr['brood_year'] + 1
    # mean_M_oc
    M_ref_col = 'mean_M_oc_sog'

# # STREAM TYPE
# mort_f = "salish_M_strmtp_sog_1971-2017.csv"
# m_p = os.path.join(mort_path, mort_f)
# print(m_p)
# if os.path.exists(m_p):
#     print("found")
#     df_mort_ref_chin_ws_sogfr = pd.read_csv(m_p)
#     # since it's brood year, add one year to M
#     df_mort_ref_chin_ws_sogfr['year'] = df_mort_ref_chin_ws_sogfr['brood_year'] + 1
#     # mean_M_oc
#     M_ref_col = 'mean_M_st_sog'

# Define root folder
root_folder = "C:\\Users\\Greig\\Documents\\EwE output\\GeorgiaStrait2023_v99_v6_7_0_18858_64b\\mc_EcosimScen99\\"
# C:\Users\Greig\Documents\EwE output\GeorgiaStrait2023_v99_v6_7_0_18858_64b\mc_EcosimScen97


# Define columns of interest and their labels
column_groups = {
    "Smolt to Age 2 M - Chinook WO": ["27", "28", "29"],  # Summing these columns
}
# # Define columns of interest and their labels
# column_groups = {
#     "Smolt to Age 2 M - Chinook WS": ["34"],  # Summing these columns
# }


# Find all Monte Carlo output folders
mc_folders = [f for f in os.listdir(root_folder) if f.startswith("mc_output_trial")]

# Initialize a dictionary to store aggregated data
aggregated_data = {}

for label, columns in column_groups.items():
    aggregated_data[label] = {}  # Store data for each group

    for folder in mc_folders:
        file_path = os.path.join(root_folder, folder, "mortality_annual.csv")

        if os.path.exists(file_path):
            # Read the file
            df = pd.read_csv(file_path)

            # Ensure 'year/group' exists
            if "year\group" in df.columns:
                df["year\group"] = df["year\group"].astype(int)  # Convert to integer years

                # Sum the selected columns
                df[label] = df[columns].sum(axis=1)

                # Store results in the dictionary
                for year, value in zip(df["year\group"], df[label]):
                    if year not in aggregated_data[label]:
                        aggregated_data[label][year] = []
                    aggregated_data[label][year].append(value)

# Create a fan plot
for label, yearly_data in aggregated_data.items():
    years = sorted(yearly_data.keys())
    data_matrix = np.array([yearly_data[year] for year in years])

    # Compute percentiles
    percentiles = np.percentile(data_matrix, [5, 25, 50, 75, 95], axis=1)

    # Plot
    plt.figure(figsize=(10, 5))
    plt.fill_between(years, percentiles[0], percentiles[4], color="lightblue", alpha=0.3, label="5th-95th Percentile")
    plt.fill_between(years, percentiles[1], percentiles[3], color="blue", alpha=0.5, label="25th-75th Percentile")
    plt.plot(years, percentiles[2], color="black", linestyle="--", label="Median")

    # Overlay reference data as points
    plt.scatter(df_mort_ref_chin_wo_sogfr['year'],
                df_mort_ref_chin_wo_sogfr[M_ref_col],
                color="red", label="M ref", zorder=3)

    plt.xlabel("Year")
    plt.ylabel("Summed Mortality")
    plt.title(f"Fan Plot for {label}")
    plt.legend()
    plt.grid(True)
    plt.show()