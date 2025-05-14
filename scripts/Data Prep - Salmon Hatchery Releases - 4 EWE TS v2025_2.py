# File: HatcheryRel_Ecosim_TS_mo_2025.csv
# Purpose: Monthly hatchery release time series for input to Ecosim model
# Author: G. Oldford
# Last Modified: May 13, 2025
#
# Description:
# This file contains time series data of hatchery-origin juvenile salmon releases
# by species, formatted for Ecosim's monthly timestep input. Biomass (MT/km²) and
# total number released (N) are included for five major salmon groups: Chinook,
# Coho, Chum, Pink, and Sockeye.
#
# Data Sources:
# - EPAD Hatchery release data (DFO/SEP) aggregated by species and release date.
# - Release timing and numbers from 'HatcheryRel_TS_ForNextstep_2025.csv'.
# - Biomass estimated using species-specific weights and normalized to study area.
#
# Notes:
# - Biomass values (MT/km²) are calculated per month per species.
# - Time series uses calendar month resolution from Jan 1950 to Dec 2020.
# - Months without hatchery releases are filled with zeros.
# - Header rows ('Title', 'Weight', 'Type', 'Timestep') follow EwE format conventions.
#   * Type = -1 (biomass forcing)
#   * Weight = 1 (uniform)
# - TIMESTEP is a linear monthly counter starting from 1 in Jan 1950.
#
# Columns:
# YEAR                - Calendar year
# TIMESTEP            - Ecosim timestep (1 = Jan 1950, 2 = Feb 1950, ..., 852 = Dec 2020)
# CHIN_H_MT_KM2       - Chinook hatchery biomass (MT/km²)
# COHO_H_MT_KM2       - Coho hatchery biomass (MT/km²)
# CHUM_H_MT_KM2       - Chum hatchery biomass (MT/km²)
# SOCKEYE_H_MT_KM2    - Sockeye hatchery biomass (MT/km²)
# PINK_H_MT_KM2       - Pink hatchery biomass (MT/km²)
# CHIN_H_N            - Chinook hatchery abundance (number released)
# COHO_H_N            - Coho hatchery abundance (number released)
# CHUM_H_N            - Chum hatchery abundance (number released)
# SOCKEYE_H_N         - Sockeye hatchery abundance (number released)
# PINK_H_N            - Pink hatchery abundance (number released)
#
# Units:
# - Biomass: metric tonnes per square kilometer (MT/km²)
# - Abundance: number of fish
#

import pandas as pd
import numpy as np
from decimal import Decimal, ROUND_UP
import copy

# --- Parameters ---
study_area = 11274
start_year = 1950
end_year = 2020

# --- Read and Prepare Raw Data ---
localpath_in = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//1. Salmon//All Species Hatchery Releases//EPADHatcherReleasesGST"
releases_df = pd.read_csv(localpath_in + "//MODIFIED//HatcheryRel_TS_ForNextstep_2025.csv")

# Convert to datetime, extract year and month
releases_df['RELEASE_DATE'] = pd.to_datetime(releases_df['release_avg_date'])
releases_df['YEAR'] = releases_df['RELEASE_DATE'].dt.year
releases_df['MONTH'] = releases_df['RELEASE_DATE'].dt.month
releases_df['EWE_TIMESTEP'] = (releases_df['YEAR'] - start_year) * 12 + releases_df['MONTH']

# --- Fix biomass estimates ---
coho_weight = 0.020  # kg
chin_weight = 0.0062 # kg

releases_df.loc[releases_df['SPECIES_NAME'] == 'Chinook', 'BIOMASS_MT2'] = releases_df['TOTRELEASE_NO'] * chin_weight * 0.001
releases_df.loc[releases_df['SPECIES_NAME'] == 'Coho', 'BIOMASS_MT2'] = releases_df['TOTRELEASE_NO'] * coho_weight * 0.001
releases_df.loc[~releases_df['SPECIES_NAME'].isin(['Coho', 'Chinook']), 'BIOMASS_MT2'] = releases_df['BIOMASS_MT']

# Normalize biomass by study area
releases_df['BIOMASS_MT2'] = releases_df['BIOMASS_MT2'] / study_area
releases_df['BIOMASS_MT']  = releases_df['BIOMASS_MT'] / study_area

# Round values to 5 decimal places
releases_df['BIOMASS_MT2'] = releases_df['BIOMASS_MT2'].apply(lambda x: Decimal(str(x)).quantize(Decimal('.00001'), rounding=ROUND_UP))
releases_df['BIOMASS_MT']  = releases_df['BIOMASS_MT'].apply(lambda x: Decimal(str(x)).quantize(Decimal('.00001'), rounding=ROUND_UP))

# --- Aggregate to species level (across all areas) ---
species_groupcodes = pd.DataFrame({
    'EWE_GROUP_CODE': ['CHINOOK-H-1', 'COHO-H-1', 'CHUM-H', 'PINK-H', 'SOCKEYE-H'],
    'SPECIES': ['CHINOOK', 'COHO', 'CHUM', 'PINK', 'SOCKEYE']
})

releases_df = releases_df.drop(columns=['release_avg_date', 'FINAL_LAT', 'FINAL_LON',
                                        'ROW_EWE', 'COL_EWE', 'SOURCE_ID', 'RELEASE_DATE'])

releases_df = releases_df.groupby(['SPECIES_NAME', 'EWE_TIMESTEP', 'YEAR', 'MONTH']).sum().reset_index()
releases_df['EWE_GROUP_CODE'] = releases_df['SPECIES_NAME']
releases_df = pd.merge(releases_df, species_groupcodes, left_on='EWE_GROUP_CODE', right_on='EWE_GROUP_CODE', how='left')
releases_df = releases_df.drop(columns=['SPECIES_NAME'])
releases_df['SPECIES'] = releases_df['EWE_GROUP_CODE'].str.upper()

# --- Create full monthly time series (fill missing months with zero) ---
timesteps = pd.DataFrame({'EWE_TIMESTEP': range(1, (end_year - start_year + 1) * 12 + 1)})

species_names = species_groupcodes['SPECIES'].unique()
all_combos = pd.MultiIndex.from_product([timesteps['EWE_TIMESTEP'], species_names], names=['EWE_TIMESTEP', 'SPECIES'])
full_df = pd.DataFrame(index=all_combos).reset_index()
full_df['YEAR'] = ((full_df['EWE_TIMESTEP'] - 1) // 12) + start_year

releases_df = pd.merge(full_df, releases_df, on=['EWE_TIMESTEP', 'YEAR', 'SPECIES'], how='left')
releases_df = releases_df.fillna(0)

# --- Pivot for Ecosim ---
df_ecosim = releases_df[['EWE_TIMESTEP', 'YEAR', 'SPECIES', 'TOTRELEASE_NO', 'BIOMASS_MT2']]
df_ecosim = df_ecosim.rename(columns={'EWE_TIMESTEP': 'TIMESTEP', 'SPECIES': 'TITLE'})

# Fix species names to uppercase to match expected keys
df_ecosim['TITLE'] = df_ecosim['TITLE'].str.upper()

pivot = df_ecosim.pivot_table(index=['TIMESTEP', 'YEAR'],
                               columns='TITLE',
                               values=['TOTRELEASE_NO', 'BIOMASS_MT2'],
                               aggfunc='sum',
                               fill_value=0).reset_index()

# --- Flatten columns ---
pivot.columns = ['_'.join(col).strip('_') for col in pivot.columns.values]

# --- Rename to match Ecosim expectation ---
col_map = {
    'BIOMASS_MT2_CHINOOK': 'CHIN_H_MT_KM2',
    'BIOMASS_MT2_COHO': 'COHO_H_MT_KM2',
    'BIOMASS_MT2_CHUM': 'CHUM_H_MT_KM2',
    'BIOMASS_MT2_PINK': 'PINK_H_MT_KM2',
    'BIOMASS_MT2_SOCKEYE': 'SOCKEYE_H_MT_KM2',
    'TOTRELEASE_NO_CHINOOK': 'CHIN_H_N',
    'TOTRELEASE_NO_COHO': 'COHO_H_N',
    'TOTRELEASE_NO_CHUM': 'CHUM_H_N',
    'TOTRELEASE_NO_PINK': 'PINK_H_N',
    'TOTRELEASE_NO_SOCKEYE': 'SOCKEYE_H_N'
}
pivot = pivot.rename(columns=col_map)

# --- Final column order ---
final_cols = ['YEAR', 'TIMESTEP'] + list(col_map.values())
pivot = pivot[final_cols]

# --- EwE header rows ---
header_lines = [
    "Title," + ",".join(final_cols[2:]),
    "Weight," + ",".join(['1'] * (len(final_cols) - 2)),
    "Type," + ",".join(['-1'] * (len(final_cols) - 2)),
    "Timestep," + ",".join(['Interval'] * (len(final_cols) - 2))
]

# --- Write final file ---
output_file = '../data/forcing/HatcheryRel_Ecosim_TS_mo_2025.csv'

with open(output_file, 'w', newline='') as f:
    for line in header_lines:
        f.write(line + '\n')
    pivot.to_csv(f, index=False, header=False)

print("done")