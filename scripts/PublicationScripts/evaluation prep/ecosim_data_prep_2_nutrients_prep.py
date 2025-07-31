"""
Preparation of nutrient data to use in evaluation
Author: G Oldford, 2025

Description:
Prepares nutrient data from two sources
- Citizen Science Data 2015 - 2023 from SoG (https://soggy2.zoology.ubc.ca/geonetwork/srv/eng/catalog.search#/metadata/dc72987d-2eb0-4b36-a6bf-62d04272bf29)
- IOS Bottle Data 1930 - present (https://catalogue.cioos.ca/dataset/ca-cioos_58418aab-008a-4ed9-be71-2f7679b6cbf7)

Input:
- CSV (from XLSX) pre-filtered to only include EwE domain area (used QGIS to do this) and filtering empty records

Outputs:
- combined CSV containing Nitrate and Nitrite

Note:
"""
import pandas as pd
import numpy as np
import os
from datetime import datetime
from matplotlib.path import Path
import matplotlib.pyplot as plt
from helpers import read_sdomains
import matplotlib
matplotlib.use('TkAgg')

# -------------------------------
# Configuration
# -------------------------------

BASE_PATH = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/29. Oceanographic Atmospheric/"
FILE_CSOP = "CitSci_Nutrients_2015-2023_justEwEDomain.csv"
PATH_CSOP = "CSOP Nutrients/MODIFIED"
FILE_IOS = "IOS_BOT_1978to2019_JustEwEArea.csv"
PATH_IOS = "IOS Rosette Bottle Data/MODIFIED"
OUTPUT_PATH = "..//..//..//data/evaluation"
OUTPUT_PATH_FIGS = "..//..//..//figs"


# -------------------------------
# Load data
# -------------------------------

# Define full paths
csop_path = os.path.join(BASE_PATH, PATH_CSOP, FILE_CSOP)
ios_path = os.path.join(BASE_PATH, PATH_IOS, FILE_IOS)

# Read CSV files
csop_df = pd.read_csv(csop_path)
ios_df = pd.read_csv(ios_path)

# -------------------------------
# Process CSOP data
# -------------------------------

# Select needed columns and rename for consistency
csop_df_sub = csop_df[["date", "latitude", "longitude", "depth", "no3"]].copy()
csop_df_sub.rename(columns={"no3": "nitrogen"}, inplace=True)
csop_df_sub["dataset"] = "csop"
# Adjust CSOP depth if 0 to 0.1
csop_df_sub.loc[csop_df_sub["depth"] == 0, "depth"] = 0.1

# -------------------------------
# Process IOS data
# -------------------------------

# Clean IOS latitude/longitude if needed (remove leading quotes)
ios_df["latitude"] = ios_df["latitude"].astype(str).str.replace("'", "").astype(float)
ios_df["longitude"] = ios_df["longitude"].astype(str).str.replace("'", "").astype(float)

# Convert IOS time to date only
if "time" in ios_df.columns:
    ios_df["date"] = pd.to_datetime(ios_df["time"]).dt.date

# Select needed columns and rename for consistency
ios_df_sub = ios_df[["date", "latitude", "longitude", "depth", "NTRZAAZ1"]].copy()
ios_df_sub.rename(columns={"NTRZAAZ1": "nitrogen"}, inplace=True)
ios_df_sub["dataset"] = "ios"

# -------------------------------
# Combine datasets
# -------------------------------

combined_df = pd.concat([csop_df_sub, ios_df_sub], ignore_index=True)

# -------------------------------
# Add month and season columns
# -------------------------------

combined_df['date'] = pd.to_datetime(combined_df['date'])
combined_df['month'] = combined_df['date'].dt.month

# Define seasons
# matches choices by McEwan et al 2023
season_map = {1: 'Winter', 2: 'Winter', 3: 'Spring', 4: 'Spring',
              5: 'Spring', 6: 'Summer', 7: 'Summer', 8: 'Summer',
              9: 'Fall', 10: 'Fall', 11: 'Winter', 12: 'Winter'}
combined_df['season'] = combined_df['month'].map(season_map)

# -------------------------------
# Export to CSV
# -------------------------------

output_file = os.path.join(OUTPUT_PATH, "nutrients_ios_csop_combined_sampled.csv")
combined_df.to_csv(output_file, index=False)

print(f"Combined dataset saved to: {output_file}")

# -------------------------------
# Generate seasonal panel plots
# -------------------------------

# [Binning depths by 2.5m]
combined_df['depth_bin'] = (combined_df['depth'] / 2.5).round() * 2.5

# Generate seasonal panel plots with mean lines and error bars
seasons = ['Winter', 'Spring', 'Summer', 'Fall']
fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)
axes = axes.flatten()

for ax, season in zip(axes, seasons):
    season_df = combined_df[combined_df['season'] == season]
    for dataset in ['ios', 'csop']:
        subset = season_df[season_df['dataset'] == dataset]
        grouped = subset.groupby('depth_bin').agg(
            mean_nitrogen=('nitrogen', 'mean'),
            std_nitrogen=('nitrogen', 'std')
        ).reset_index().sort_values('depth_bin')
        ax.errorbar(grouped['mean_nitrogen'], grouped['depth_bin'],
                    xerr=grouped['std_nitrogen'], fmt='-o', label=f"{dataset}")
    ax.set_title(season)
    ax.set_xlabel('Nitrogen (umol/L)')
    ax.set_ylabel('Depth (m)')
    ax.legend()
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

# -------------------------------
# Generate bi-weekly averaged nitrogen time series
# -------------------------------

combined_df['date'] = pd.to_datetime(combined_df['date'])
combined_df['julian_day'] = combined_df['date'].dt.dayofyear

# Define bi-weekly bins (14-day intervals)
combined_df['biweek'] = ((combined_df['julian_day'] - 1) // 14 + 1) * 14

# Define shallow vs deep threshold (e.g. 10 m)
shallow_threshold = 10

# Group by biweek and depth_bin, calculate mean nitrogen per bin
biweek_depth_grouped = combined_df.groupby(['biweek','depth_bin']).agg(mean_nitrogen=('nitrogen','mean')).reset_index()

# Determine if both shallow and deep bins are present in each biweek
biweek_depth_grouped['is_shallow'] = biweek_depth_grouped['depth_bin'] <= shallow_threshold
biweek_depth_grouped['is_deep'] = biweek_depth_grouped['depth_bin'] > shallow_threshold

# Prepare list to store valid biweeks
valid_biweeks = []
for biweek, group in biweek_depth_grouped.groupby('biweek'):
    has_shallow = group['is_shallow'].any()
    has_deep = group['is_deep'].any()
    if has_shallow and has_deep:
        valid_biweeks.append(biweek)

# Filter to valid biweeks only
biweek_depth_grouped_valid = biweek_depth_grouped[biweek_depth_grouped['biweek'].isin(valid_biweeks)]

# Calculate overall MEAN treating each depth_bin equally for valid biweeks only
biweek_overall_mean = biweek_depth_grouped_valid.groupby('biweek').agg(overall_mean_nitrogen=('mean_nitrogen','mean')).reset_index()

# Plot time series
plt.figure(figsize=(10,6))
plt.plot(biweek_overall_mean['biweek'], biweek_overall_mean['overall_mean_nitrogen'], '-o')
plt.xlabel('Julian Day of Year (bi-weekly bins)')
plt.ylabel('Mean Nitrogen (umol/L)')
plt.title('Bi-weekly Mean Nitrogen')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_PATH_FIGS, 'biweekly_nitrogen_mean_filtered.png'))
plt.show()


# Print overall mean nitrogen concentration for filtered dataset
overall_mean_filtered = biweek_overall_mean['overall_mean_nitrogen'].mean()
print(f"Mean nitrogen concentration across all filtered biweekly blocks: {overall_mean_filtered:.2f} umol/L")

# MEDIAN VERSION
# Calculate overall mean treating each depth_bin equally for valid biweeks only
biweek_overall_median = biweek_depth_grouped_valid.groupby('biweek').agg(overall_median_nitrogen=('mean_nitrogen','median')).reset_index()

# Plot time series
plt.figure(figsize=(10,6))
plt.plot(biweek_overall_median['biweek'], biweek_overall_median['overall_median_nitrogen'], '-o')
plt.xlabel('Julian Day of Year (bi-weekly bins)')
plt.ylabel('Median Nitrogen (umol/L)')
plt.title('Bi-weekly Median Nitrogen')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_PATH_FIGS, 'biweekly_nitrogen_median_filtered.png'))
plt.show()

overall_median_filtered = biweek_overall_median['overall_median_nitrogen'].median()
print(f"Median nitrogen concentration across all filtered biweekly blocks: {overall_median_filtered:.2f} umol/L")