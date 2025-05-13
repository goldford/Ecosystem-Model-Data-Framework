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

# need a purge script since MC generates so much shite
# purge everything without 'annual' in name
# purge everything with 'discards' in name


# mortality time series (for ref, fitting)
# Chinook - ocean and stream type, Freshwater et al., 2022

# SCENARIO = "SC100 - NO PP Anom"
SCENARIO = "SC103 - WITH PP Anom"
# mc_root_p = "C:\\Users\\Greig\\Documents\\EwE output\\GeorgiaStrait2023_v100_v6_7_0_18858_64b\\mc_EcosimScen100\\"
mc_root_p = "C:\\Users\\Greig\\Documents\\EwE output\\GeorgiaStrait2023_v103_v6_7_0_18858_64b\\mc_EcosimScen103\\"
mort_chin_cwt_p = "C:\\Users\\Greig\\Sync\\6. SSMSP Model\\Model Greig\\Data\\1. Salmon\\Chinook Survival Freshwater 2022\\MODIFIED\\"
mort_coho_nelson_p = "C:\\Users\\Greig\\Sync\\6. SSMSP Model\\Model Greig\\Data\\1. Salmon\\Coho_Mort_Seal_Nelson_2024\ORIGINAL\\"
mort_coho_cwt_p = "C:\\Users\\Greig\\Sync\\6. SSMSP Model\\Model Greig\\Data\\1. Salmon\\All Species Survival Exploitation\\"
# C:\Users\Greig\Documents\EwE output\GeorgiaStrait2023_v99_v6_7_0_18858_64b\mc_EcosimScen97


# OCEAN TYPE CHINOOK
# reference CWT data
mort_f = "salish_M_ocntp_sog_1971-2017.csv"
m_p = os.path.join(mort_chin_cwt_p, mort_f)
print(m_p)
if os.path.exists(m_p):
    print("found")
    df_mort_ref_chin_wo_sogfr = pd.read_csv(m_p)
    # since it's brood year, add one year to M
    df_mort_ref_chin_wo_sogfr['year'] = df_mort_ref_chin_wo_sogfr['brood_year'] + 1
    # mean_M_oc
    M_ref_col_ch_wo = 'mean_M_oc_sog'

# STREAM TYPE CHINOOK
# reference CWT data
mort_f = "salish_M_strmtp_sog_1971-2017.csv"
m_p = os.path.join(mort_chin_cwt_p, mort_f)
print(m_p)
if os.path.exists(m_p):
    print("found")
    df_mort_ref_chin_ws_sogfr = pd.read_csv(m_p)
    # since it's brood year, add one year to M
    df_mort_ref_chin_ws_sogfr['year'] = df_mort_ref_chin_ws_sogfr['brood_year'] + 1
    # mean_M_oc
    M_ref_col_ch_ws = 'mean_M_st_sog'

#############################
# COHO
# reference M from ***CWT***
M_cwt_coho_f = "coho_M_EWE_TS Jan 2025.csv"
m_p = os.path.join(mort_coho_cwt_p, M_cwt_coho_f)
print(m_p)
if os.path.exists(m_p):
    print("found " + m_p)
    df_m_ref_coho_cwt = pd.read_csv(m_p)
    # since it's brood year, add one year to M
    df_m_ref_coho_cwt['year'] = df_m_ref_coho_cwt['year'] + 1
    # mean_M_oc
    M_ref_col_co_cwt = 'AVERAGE'

# reference M from ***Nelson***
M_nelson_f = "nelson_2024_taba4.csv"
m_p = os.path.join(mort_coho_nelson_p, M_nelson_f)
print(m_p)
if os.path.exists(m_p):
    print("found " + m_p)
    df_mort_ref_coho_nelson = pd.read_csv(m_p)
    # since it's brood year, add one year to M
    df_mort_ref_coho_nelson['year'] = df_mort_ref_coho_nelson['year'] + 1
    M_ref_col_co_nels = 'm_seals_coho'


###############################################
# Define columns of interest and their labels
# Chinook WO
column_groups = {
    "Smolt to Age 2 M - Chinook WO": ["27", "28", "29"],  # Summing these columns
    "Smolt to Age 2 M - Chinook WS": ["34"],
    "Smolt to Age 2 M - Coho W": ["43"]
}


# Find all Monte Carlo output folders
mc_folders = [f for f in os.listdir(mc_root_p) if f.startswith("mc_output_trial")]

# Initialize a dictionary to store aggregated data
aggregated_data = {}

for label, columns in column_groups.items():
    aggregated_data[label] = {}  # Store data for each group
    aggregated_data[label + " Seal M"] = {}

    for folder in mc_folders:
        file_path = os.path.join(mc_root_p, folder, "mortality_annual.csv")

        if os.path.exists(file_path):

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

        # repeated loop as above to grab predation files to look at M from specific preds
        file_path = os.path.join(mc_root_p, folder, "predation_coho2-w-emar_annual.csv")
        pred_col = 'Harbour Seals'

        if os.path.exists(file_path):

            df_pred = pd.read_csv(file_path)
             # Store data for each group
            # Ensure 'year/group' exists
            if "year\group" in df_pred.columns:
                df_pred["year\group"] = df_pred["year\group"].astype(int)  # Convert to integer years

                # # Sum the selected columns
                # df_pred[label] = df_pred[columns].sum(axis=1)
                # df_pred[label + " Seal M"] = df_pred[pred_col]

                # Store results in the dictionary
                for year, value in zip(df_pred["year\group"], df_pred[pred_col]):
                    if year not in aggregated_data[label + " Seal M"]:
                        aggregated_data[label + " Seal M"][year] = []
                    aggregated_data[label + " Seal M"][year].append(value)


######## CHIN #########
# fan plot
for label, yearly_data in aggregated_data.items():

    # Overlay reference data as points
    if label == "Smolt to Age 2 M - Chinook WO":
        df_ref_cwt = df_mort_ref_chin_wo_sogfr
        M_ref_col = M_ref_col_ch_wo
        min_M_norm = 3.1
        max_M_norm = 4.2
    elif label == "Smolt to Age 2 M - Chinook WS":
        df_ref_cwt = df_mort_ref_chin_ws_sogfr
        M_ref_col = M_ref_col_ch_ws
        min_M_norm = 2.3
        max_M_norm = 3
    elif label == "Smolt to Age 2 M - Coho W":
        df_ref_cwt = df_m_ref_coho_cwt
        M_ref_col = M_ref_col_co_cwt
        min_M_norm = 2.3
        max_M_norm = 3
    else:
        print("skipping label with seal M for now (see below)")
        continue # jump out if it's seal M - instead we match below for each other label


    years = sorted(yearly_data.keys())
    data_matrix = np.array([yearly_data[year] for year in years])

    # Compute percentiles
    percentiles = np.percentile(data_matrix, [5, 25, 50, 75, 95], axis=1)

    # Plot total M
    plt.figure(figsize=(10, 7))
    plt.fill_between(years, percentiles[0]*0.85, percentiles[4]*1.15, color="lightblue", alpha=0.3, label="")
    plt.fill_between(years, percentiles[1]*0.95, percentiles[3]*1.05, color="blue", alpha=0.5, label="")
    plt.plot(years, percentiles[2], color="black", linestyle="--", label="Total M - Model")



    # plot seal M
    label_sealM = label + " Seal M"
    yearly_sealM_data = aggregated_data[label_sealM]
    years_sealM = sorted(yearly_sealM_data.keys())
    data_matrix_sealM = np.array([yearly_sealM_data[year] for year in years_sealM])
    percentiles_sealM = np.percentile(data_matrix_sealM, [5, 25, 50, 75, 95], axis=1)
    # Plot total M
    plt.fill_between(years_sealM, percentiles_sealM[0]*0.85, percentiles_sealM[4]*1.15, color="lightblue", alpha=0.3, label="")
    plt.fill_between(years_sealM, percentiles_sealM[1]*0.95, percentiles_sealM[3]*1.05, color="blue", alpha=0.5, label="")
    plt.plot(years_sealM, percentiles_sealM[2], color="blue", linestyle="--", label="M Seals - Model")

    # plot CWT reference M
    plt.scatter(df_ref_cwt['year'],
                df_ref_cwt[M_ref_col],
                color="red", label="Smolt to Age 2 M - CWT", zorder=3)

    plt.axhspan(min_M_norm, max_M_norm, color='gray', alpha=0.3, label="Normal?")

    # if coho, compare to nelson predictions
    if label == "Smolt to Age 2 M - Coho W":
        plt.plot(df_mort_ref_coho_nelson['year'], df_mort_ref_coho_nelson[M_ref_col_co_nels], color="red", linestyle="--", label="Seal M - Nelson 2024")
        plt.fill_between(df_mort_ref_coho_nelson['year'], df_mort_ref_coho_nelson['m_low_ci'], df_mort_ref_coho_nelson['m_hi_ci'], color="red", alpha=0.3,
                         label="5th-95th Percentile")

    plt.xlabel("Year")
    plt.ylabel("Mortality (/yr)")
    plt.title(f"Fan Plot for {label} - " + SCENARIO)
    plt.legend()
    plt.grid(True)
    plt.xlim(1970, 2018)
    plt.show()

    # FOR PPT - just plot our est versus Nelson's
    if label == "Smolt to Age 2 M - Coho W":
        # plot seal M
        label_sealM = label + " Seal M"
        yearly_sealM_data = aggregated_data[label_sealM]
        years_sealM = sorted(yearly_sealM_data.keys())
        data_matrix_sealM = np.array([yearly_sealM_data[year] for year in years_sealM])
        percentiles_sealM = np.percentile(data_matrix_sealM, [5, 25, 50, 75, 95], axis=1)

        plt.figure(figsize=(10, 7))
        # Plot total M
        plt.fill_between(years_sealM, percentiles_sealM[0]*0.5, percentiles_sealM[4]*2, color="lightblue", alpha=0.3,
                         label="")
        plt.fill_between(years_sealM, percentiles_sealM[1]*0.8, percentiles_sealM[3]*1.2, color="blue", alpha=0.5, label="")
        plt.plot(years_sealM, percentiles_sealM[2], color="blue", linestyle="--", label="M from Seals - Model")

        plt.plot(df_mort_ref_coho_nelson['year'], df_mort_ref_coho_nelson[M_ref_col_co_nels], color="red", linestyle="--", label="M from Seals - Nelson et al. 2024")
        plt.fill_between(df_mort_ref_coho_nelson['year'], df_mort_ref_coho_nelson['m_low_ci'], df_mort_ref_coho_nelson['m_hi_ci'], color="red", alpha=0.3,
                         label="")

        #plot CWT reference M
        plt.scatter(df_ref_cwt['year'],
                    df_ref_cwt[M_ref_col],
                    color="red", label="M - CWT", zorder=3)

        plt.axhspan(2.3, 3, color='gray', alpha=0.3, label="Normal?")

        plt.xlabel("Year")
        plt.ylabel("Mortality (/yr)")
        plt.title(f"Fan Plot for {label} - " + SCENARIO)
        plt.xlim(1970, 2018)
        plt.legend()
        plt.grid(True)
        plt.show()

