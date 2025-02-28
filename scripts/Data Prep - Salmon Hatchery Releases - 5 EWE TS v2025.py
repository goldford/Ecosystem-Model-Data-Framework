### Step 5 - Modify Ecosim TS for Model
#  G Oldford (2025)
# Purpose:
#    - Open the TS from previous, and adjust for directly pasting into EwE
#    - Adjusting to make TS for MS groups relative to base year
#    - Also adjusting based on B ratios to N, with details in the tech report
#
# Data In:
#   - HatcheryRel_Ecosim_TS_yr_2025.csv
#   - HatcheryRel_Ecosim_TS_mo_2025.csv
#
# Data Out:
#  -
#
# Notes:
# - CHECK THE COL NAMES IN INPUT FOR TRUNCATION

import pandas as pd
import numpy as np
import datetime
from dateutil.parser import parse
from decimal import *

pd.set_option('display.max_columns', None)

# 1. open time series CSV
# 2. create new fields for each species (coho, chin, sockeye, pink, chum).
#    there are two fields for each species
#    eg, [SPECIES]_H_MT_KM2 becomes [SPECIES]_H_MT_KM2_adj
#    eg, [SPECIES]_H_MT_KM2 becomes [SPECIES]_H_MT_KM2_adj
# 3. skip header (4 lines)
# 4. For the _first_ year (ie, row) only, set baseline 'fake' hatchery N,B
#    for columns suffixed _N, values for first year should be, eg, [SPECIES]_H_N = N_1950_W_[species] * 0.2
#    for columns suffixed _MT_KM2, values for first year should be, eg, x = [SPECIES]_H_N * B_1950_W_[species] / N_1950_W_[species]
# 5. For all years, except first one

# # juve B's and N's for MS based on _FRESH_ stanza for ratio below
# # source: see tech report
# B_1950_W_coho   = 0.055
# B_1950_W_chin_o = 0.0000627
# B_1950_W_pink   = 0.0034
# B_1950_W_sock   = 0.066
# B_1950_W_chum   = 0.0032
#
# N_1950_W_coho   = 15500000
# N_1950_W_chin_o = 35000000
# N_1950_W_pink   = 227000000
# N_1950_W_sock   = 66000000
# N_1950_W_chum   = 82000000

# Define biomass (B) and abundance (N) values for 1950
B_1950_W = {
    "COHO": 0.055,
    "CHIN": 0.0000627,
    "PINK": 0.0034,
    "SOCKEYE": 0.066,
    "CHUM": 0.0032
}

N_1950_W = {
    "COHO": 15500000,
    "CHIN": 35000000,
    "PINK": 227000000,
    "SOCKEYE": 66000000,
    "CHUM": 82000000
}

yr_strt = 1950
yr_end = 2020

df = pd.read_csv("../data/forcing/HatcheryRel_Ecosim_TS_yr_2025.csv", skiprows=[1,2,3])
print(df.columns)
# Index(['Title', 'YEAR', 'TIMESTEP', 'CHIN_H_MT_KM2', 'COHO_H_MT_KM2',
#        'CHUM_H_MT_KM2', 'SOCKEYE_H_MT_KM2', 'PINK_H_MT_KM2', 'CHIN_H_N',
#        'COHO_H_N', 'CHUM_H_N', 'SOCKEYE_H_N', 'PINK_H_N'],
#       dtype='object')
# Ensure the dataframe has a year column
if "YEAR_" not in df.columns:
    raise ValueError("CSV file must contain a 'year' column.")

# Create new adjusted columns
species_list = ["COHO", "CHIN", "PINK", "SOCKEYE", "CHUM"]

for species in species_list:
    df[f"{species}_H_N_adj"] = np.nan
    df[f"{species}_H_MT_KM2_adj"] = np.nan

# Process first year
first_year_idx = df[df["YEAR_"] == yr_strt].index[0]

for species in species_list:
    # Adjust _N column
    df.at[first_year_idx, f"{species}_H_N_adj"] = N_1950_W[species] * 0.2

    # Adjust _MT_KM2 column
    x = df.at[first_year_idx, f"{species}_H_N_adj"] * B_1950_W[species] / N_1950_W[species]
    df.at[first_year_idx, f"{species}_H_MT_KM2_adj"] = x

# Process each species
for species in species_list:
    N_col = f"{species}_H_N_"  # Original N column
    N_adj_col = f"{species}_H_N_adj"  # Adjusted column

    if N_col in df.columns:

        # Get the first year's value
        first_year_value = df.loc[df.index[0], N_adj_col]
        print(first_year_value)

        # Apply the logic: set to the real hatchery number unless hatchery prduction falls below 1/5 natural
        df[N_adj_col] = np.where(df[N_col] < first_year_value, first_year_value, df[N_col])

    B_col = f"{species}_H_MT_"
    B_adj_col = f"{species}_H_MT_KM2_adj"  # Adjusted column

    # logic - establish original ratio and apply to get hatchery B
    # reasoning in tech report, based on B being time integrated in model params but B being a point in time in H data
    B_N_ratio = B_1950_W[species] / N_1950_W[species]
    print(B_N_ratio)

    if B_col in df.columns:

        # Apply the logic: set to hatchery B to number of hatchery (adjusted) * weight in mt km2 of each fish
        df[B_adj_col] = df[N_adj_col] * B_N_ratio

print(df)

df.to_csv('../scratch/temp_yr_adj.csv', index=True)
df.to_csv('../data/forcing/HatcheryRel_Ecosim_TS_yr_2025_adj.csv', index=True)







# # params
# start_year = 1950
# end_year = 2020
#
# aggregate_time = "year" # month or year
# aggregate_all_areas = "yes" # yes means aspatial
# aggregate_to_level = "species" # otherwise will use codes in EWE_GROUP_CODE
# species_groupcodes = pd.DataFrame(data = {'EWE_GROUP_CODE': ['CHINOOK-H-1','COHO-H-1',
#                                                              'CHUM-H', 'PINK-H',
#                                                              'SOCKEYE-H'],
#                                           'SPECIES':['CHINOOK','COHO', 'CHUM', 'PINK', 'SOCKEYE']})
#
# # locations table from the SSMSP SOGDC (may have more lats / lons added than source at RMIS)
# localpath_in = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/1. Salmon/All Species Hatchery Releases/EPADHatcherReleasesGST"
# releases_df = pd.read_csv(localpath_in + "/MODIFIED/HatcheryRel_Ecosim_TS_2025.csv")
#
# print(releases_df)


# # cross check for one year
# print(releases_df.loc[releases_df['YEAR']==1980].groupby(['EWE_GROUP_CODE','YEAR']).sum().reset_index())
# # chinook in 1980 should be 6930550, coho in 1980 should be 4552663 - GO 01152025
#
# # round to 5 decimal places due to issues with floating point data
# # storage causing rounding to not work so
# # using decimal library https://stackoverflow.com/questions/56820/round-doesnt-seem-to-be-rounding-properly
#
# releases_df['BIOMASS_MT2'] = releases_df['BIOMASS_MT2'].apply(lambda x: Decimal(str(x)).quantize(Decimal('.00001'), rounding=ROUND_UP))
# releases_df['BIOMASS_MT'] = releases_df['BIOMASS_MT'].apply(lambda x: Decimal(str(x)).quantize(Decimal('.00001'), rounding=ROUND_UP))
#
# # add dummy variable containing all timesteps
# dummy = pd.Series(range(1,((end_year - start_year)*12)))
# dummy_df = (dummy.to_frame())
# dummy_df['EWE_TIMESTEP'] = dummy_df[0]
# dummy_df['EWE_GROUP_CODE'] = "DUMMY"
# dummy_df['YEAR'] = (dummy_df['EWE_TIMESTEP'] // 12)+start_year
#
# dummy_df = dummy_df[['EWE_GROUP_CODE','EWE_TIMESTEP','YEAR']]
#
# #releases_df = releases_df.append(dummy_df, ignore_index = True) # deprecated - 2025
# releases_df = pd.concat([releases_df, dummy_df], ignore_index=True)
# print(releases_df)
#
#
# ################ For Ecosim ########################
# releasesEcosim = releases_df[['EWE_TIMESTEP', 'BIOMASS_MT2', 'TOTRELEASE_NO', 'EWE_GROUP_CODE', 'YEAR']]
# releasesEcosim = releasesEcosim.fillna(0)
#
# # sum by EWE_GROUP_CODE and timestep
# releasesEcosim = releasesEcosim.rename(columns={'EWE_TIMESTEP': 'TIMESTEP','EWE_GROUP_CODE': 'TITLE'})
#
# # for timestep = monthly
# releasesEcosim_gp_mo = releasesEcosim.groupby(['TIMESTEP', 'TITLE', 'YEAR']).sum().reset_index()
#
#
# # pivot wide
# releasesEcosim_wide_mo = releasesEcosim_gp_mo.pivot_table(
#         values=['BIOMASS_MT2', 'TOTRELEASE_NO'],
#         index=['TIMESTEP', 'YEAR'],
#         columns='TITLE',
#         aggfunc=np.sum).reset_index()
#
# # reset the multilevel index via hack
# releasesEcosim_wide_mo['CHIN_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Chinook')].astype(float)
# releasesEcosim_wide_mo['COHO_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Coho')].astype(float)
# releasesEcosim_wide_mo['CHUM_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Chum')].astype(float)
# releasesEcosim_wide_mo['SOCKEYE_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Sockeye')].astype(float)
# releasesEcosim_wide_mo['PINK_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Pink')].astype(float)
#
# releasesEcosim_wide_mo['CHIN_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Chinook')].astype(float)
# releasesEcosim_wide_mo['COHO_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Coho')].astype(float)
# releasesEcosim_wide_mo['CHUM_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Chum')].astype(float)
# releasesEcosim_wide_mo['SOCKEYE_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Sockeye')].astype(float)
# releasesEcosim_wide_mo['PINK_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Pink')].astype(float)
#
# releasesEcosim_wide_mo['TIMESTEP'] = releasesEcosim_wide_mo[('TIMESTEP', '')].astype(float)
# releasesEcosim_wide_mo['YEAR'] = releasesEcosim_wide_mo[('YEAR', '')].astype(float)
# releasesEcosim_wide_mo = releasesEcosim_wide_mo[['YEAR',
#                                                  'TIMESTEP',
#                                                  'CHIN_H_MT',
#                                                  'COHO_H_MT',
#                                                  'CHUM_H_MT',
#                                                  'SOCKEYE_H_MT',
#                                                  'PINK_H_MT',
#                                                  'CHIN_H_N',
#                                                  'COHO_H_N',
#                                                  'CHUM_H_N',
#                                                  'SOCKEYE_H_N',
#                                                  'PINK_H_N'
#                                                  ]]
# releasesEcosim_wide_mo.columns = [f"{x}_{y}" for x, y in releasesEcosim_wide_mo.columns.to_flat_index()]
#
#
#
# # releasesEcosim_wide = releasesEcosim_wide.drop(columns=[('BIOMASS_MT2',   'DUMMY')])
#
# # fill NaNs with zeros (required by ecosim)
# releasesEcosim_wide_mo = releasesEcosim_wide_mo.fillna(0)
#
# #releasesEcosim_wide_mo = pd.DataFrame(releasesEcosim_wide_mo.to_records())
# print(releasesEcosim_wide_mo)
# #print(releasesEcosim_wide_mo.columns)
#
# # use average monthly for annual time series
# # releasesEcosim_wide_yr = releasesEcosim_wide_mo.groupby(['YEAR_']).mean().reset_index()
# # releasesEcosim_wide_yr = releasesEcosim_wide_yr[['YEAR_',
# #                                                  'CHIN_H_MT_',
# #                                                  'COHO_H_MT_',
# #                                                  'CHUM_H_MT_',
# #                                                  'SOCKEYE_H_MT_',
# #                                                  'PINK_H_MT_',
# #                                                  'CHIN_H_N_',
# #                                                  'COHO_H_N_',
# #                                                  'CHUM_H_N_',
# #                                                  'SOCKEYE_H_N_',
# #                                                  'PINK_H_N_'
# #                                                  ]]
#
# releasesEcosim_wide_yr = releasesEcosim_wide_mo.groupby(['YEAR_']).agg({
#     'CHIN_H_MT_': 'mean',
#     'COHO_H_MT_': 'mean',
#     'CHUM_H_MT_': 'mean',
#     'SOCKEYE_H_MT_': 'mean',
#     'PINK_H_MT_': 'mean',
#     'CHIN_H_N_': 'sum',
#     'COHO_H_N_': 'sum',
#     'CHUM_H_N_': 'sum',
#     'SOCKEYE_H_N_': 'sum',
#     'PINK_H_N_': 'sum'
# }).reset_index()
#
# # #if aggregate_time == "year":
# # # releasesEcosim_wide = releasesEcosim_wide.drop(columns="('TIMESTEP', '')", axis=1)
# # #     releasesEcosim_wide['Chinook'] = releasesEcosim_wide["('BIOMASS_MT2', 'Chinook')"].astype(float)
# # #     releasesEcosim_wide['Coho'] = releasesEcosim_wide["('BIOMASS_MT2', 'Coho')"].astype(float)
# # #     releasesEcosim_wide = releasesEcosim_wide.groupby("('YEAR', '')").mean().reset_index()
# #
# # write to temp file
# releasesEcosim_wide_yr.to_csv('../scratch/temp_yr.csv', index=True)
# releasesEcosim_wide_mo.to_csv('../scratch/temp_mo.csv', index=True)
#
# # this repeats same avg value each month, for silly workaround
# repeated_yr_avg = pd.merge(releasesEcosim_wide_mo, releasesEcosim_wide_yr, on=['YEAR_'], how='left')
# repeated_yr_avg = repeated_yr_avg[['YEAR_',
#                                    'TIMESTEP_',
#                                    'CHIN_H_MT__y',
#                                    'COHO_H_MT__y',
#                                    'CHUM_H_MT__y',
#                                    'SOCKEYE_H_MT__y',
#                                    'PINK_H_MT__y',
#                                    'CHIN_H_N__y',
#                                    'COHO_H_N__y',
#                                    'CHUM_H_N__y',
#                                    'SOCKEYE_H_N__y',
#                                    'PINK_H_N__y'
#                                    ]]
#
# repeated_yr_avg = repeated_yr_avg.rename(columns={'TIMESTEP_': 'TIMESTEP',
#                                                   'YEAR_': 'YEAR',
#                                                   'CHIN_H_MT__y': 'CHIN_H_MT_KM2',
#                                                   'COHO_H_MT__y': 'COHO_H_MT_KM2',
#                                                   'CHUM_H_MT__y': 'CHUM_H_MT_KM2',
#                                                   'PINK_H_MT__y': 'PINK_H_MT_KM2',
#                                                   'SOCKEYE_H_MT__y': 'SOCKEYE_H_MT_KM2',
#                                                   'CHIN_H_N__y': 'CHIN_H_N',
#                                                   'COHO_H_N__y': 'COHO_H_N',
#                                                   'CHUM_H_N__y': 'CHUM_H_N',
#                                                   'PINK_H_N__y': 'PINK_H_N',
#                                                   'SOCKEYE_H_N__y': 'SOCKEYE_H_N'
#                                                   })
# repeated_yr_avg.to_csv('../scratch/temp_yr_rep.csv', index=True)
#
# # ===================================
# # open temp file and insert header
# # ===================================
#
# # Title	Combined_GST_FR_Escape_RelB_NuSEDS	Chin_Hatch_RelB_CW	Chin_1stYrM_1_CW	Chin_1stYrM_2_CW	Chin_C_Rel_CW
# # Weight	1	1	1	1	1
# # Pool Code	14	18	16	15	14
# # Type	0	0	5	5	61
# # 1979	11.26655002	3.84	3.449022245	3.449022245	0.35
# # 1980	11.07767237	6.93	3.021428984	3.021428984	0.371
# # 1981	11.23108247	8.75	3.354206073	3.354206073	0.2533
#
# # codes for 'type'
# # relative biomass = 0
# # absolute biomass = 1
# # biomass forcing = -1
# # fishing mortality = 4
# # relative fishing mortality = 104
# # total mortality = 5
# # constant total mortality = -5 (forcing?)
# # catches = 6
# # catches forcing = -6
# # relative catches = 61
# # average weight = 7
#
# import copy
#
# f = open('../scratch/temp_yr.csv', "r")
# contents = f.readlines()
# f.close()
#
# line1 = contents[0].split(',')
# line1[0] = 'Title'
#
# line2 = copy.deepcopy(line1)
# line2[0] = 'Weight'
# i = 0
# for line in line2:
#     if i > 0:
#         if i == (len(line2) - 1):
#             line2[i] = '1\n'
#         else:
#             line2[i] = 1
#     i += 1
#
# line3 = copy.deepcopy(line1)
# line3[0] = 'Type'
# i = 0
# for line in line3:
#     if i > 0:
#         if i == (len(line3) - 1):
#             line3[i] = '-1\n'
#         else:
#             line3[i] = -1
#     i += 1
#
# line4 = copy.deepcopy(line1)
# line4[0] = 'Timestep'
# i = 0
# for line in line4:
#     if i > 0:
#         if i == (len(line4) - 1):
#             line4[i] = 'Interval\n'
#         else:
#             line4[i] = 'Interval'
#     i += 1
#
# s = ""
# contents.insert(1, ','.join(str(line) for line in line1))
# contents.insert(2, ','.join(str(line) for line in line2))
# contents.insert(3, ','.join(str(line) for line in line3))
# contents.insert(4, ','.join(str(line) for line in line4))
#
# i = 0
# with open('../scratch/HatcheryRel_Ecosim_TS_1.csv', 'w') as a_writer:
#     for line in contents:
#         if i > 0:
#             a_writer.writelines(line)
#         i += 1
#
# f = open('../scratch/temp_yr_rep.csv', "r")
# contents = f.readlines()
# f.close()
#
# line1 = contents[0].split(',')
# line1[0] = 'Title'
#
# line2 = copy.deepcopy(line1)
# line2[0] = 'Weight'
# i = 0
# for line in line2:
#     if i > 0:
#         if i == (len(line2) - 1):
#             line2[i] = '1\n'
#         else:
#             line2[i] = 1
#     i += 1
#
# line3 = copy.deepcopy(line1)
# line3[0] = 'Type'
# i = 0
# for line in line3:
#     if i > 0:
#         if i == (len(line3) - 1):
#             line3[i] = '-1\n'
#         else:
#             line3[i] = -1
#     i += 1
#
# line4 = copy.deepcopy(line1)
# line4[0] = 'Timestep'
# i = 0
# for line in line4:
#     if i > 0:
#         if i == (len(line4) - 1):
#             line4[i] = 'Interval\n'
#         else:
#             line4[i] = 'Interval'
#     i += 1
#
# s = ""
# contents.insert(1, ','.join(str(line) for line in line1))
# contents.insert(2, ','.join(str(line) for line in line2))
# contents.insert(3, ','.join(str(line) for line in line3))
# contents.insert(4, ','.join(str(line) for line in line4))
#
# i = 0
# with open('../scratch/HatcheryRel_Ecosim_TS.csv', 'w') as a_writer:
#     for line in contents:
#         if i > 0:
#             a_writer.writelines(line)
#         i += 1
#
# line1 = contents[0].split(',')
# line1[0] = 'Title'
#
# line2 = copy.deepcopy(line1)
# line2[0] = 'Weight'
# i = 0
# for line in line2:
#     if i > 0:
#         if i == (len(line2) - 1):
#             line2[i] = '1\n'
#         else:
#             line2[i] = 1
#     i += 1
#
# line3 = copy.deepcopy(line1)
# line3[0] = 'Type'
# i = 0
# for line in line3:
#     if i > 0:
#         if i == (len(line3) - 1):
#             line3[i] = '-1\n'
#         else:
#             line3[i] = -1
#     i += 1
#
# line4 = copy.deepcopy(line1)
# line4[0] = 'Timestep'
# i = 0
# for line in line4:
#     if i > 0:
#         if i == (len(line4) - 1):
#             line4[i] = 'Interval\n'
#         else:
#             line4[i] = 'Interval'
#     i += 1
#
# s = ""
# contents.insert(1, ','.join(str(line) for line in line1))
# contents.insert(2, ','.join(str(line) for line in line2))
# contents.insert(3, ','.join(str(line) for line in line3))
# contents.insert(4, ','.join(str(line) for line in line4))
#
# i = 0
# with open('../data/forcing/HatcheryRel_Ecosim_TS_2025.csv', 'w') as a_writer:
#     for line in contents:
#         if i > 0:
#             a_writer.writelines(line)
#         i += 1