### Export Custom Time Series file for Ecosim / Ecospace
#  G Oldford (2021; last mod: Feb, 2025)
# Purpose: Export a hatchery forcing or time series file to .csv's or ASCII's for EWE
#    - export the monthly timestep forcing file that Ecosim expects
#    - export the monthly timestep spatial forcing file that Ecospace expects
#
# Data In:
#   - HatcheryRel_TS_ForNextStep.csv - EPAD data from Carl (DFO / SEP)- from 'step 3'
#
# Data Out:
#  - HatcheryRel_Ecosim_TS_Mo_2025.csv # monthly forcing time series (B, N hatchery)
#  - HatcheryRel_Ecosim_TS_yr_2025.csv # annual forcing time series (B, N hatchery)
#
# Notes:
# - EPAD data from Carl Walters and RMIS locations data from SOGDC
# - Apr 12, 2021 - the average weight or the weight field is off! not sure why,
#                 more so with coho so I just went to EPAD data and got avg n
#                 non-zero weight for 1980 - 1990 from Puntledge
# - Apr 12, 2021 - the annual out should be average of monthly b_mt released! Not sum over yr
# - Jan 15, 2025 - issue encountered: need spacing of monthly not annual TS for Ecosim forcing
# - Feb 24, 2025 - fixed issues, exporting numbers of hatchery fish as well as B's

import pandas as pd
import numpy as np
import datetime
from dateutil.parser import parse
from decimal import *

pd.set_option('display.max_columns', None)

# params
study_area = 11274 # used to calculate biomass density (mt / km^2) - updated 2025
start_year = 1950
end_year = 2020

aggregate_time = "year" # month or year
aggregate_all_areas = "yes" # yes means aspatial
aggregate_to_level = "species" # otherwise will use codes in EWE_GROUP_CODE
species_groupcodes = pd.DataFrame(data = {'EWE_GROUP_CODE': ['CHINOOK-H-1','COHO-H-1',
                                                             'CHUM-H', 'PINK-H',
                                                             'SOCKEYE-H'],
                                          'SPECIES':['CHINOOK','COHO', 'CHUM', 'PINK', 'SOCKEYE']})

# locations table from the SSMSP SOGDC (may have more lats / lons added than source at RMIS)
localpath_in = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/1. Salmon/All Species Hatchery Releases/EPADHatcherReleasesGST"
releases_df = pd.read_csv(localpath_in + "/MODIFIED/HatcheryRel_TS_ForNextstep_2025.csv")


# fix date
releases_df['RELEASE_DATE'] = releases_df['release_avg_date'].astype('datetime64[ns]')
releases_df['YEAR'] = pd.DatetimeIndex(releases_df['RELEASE_DATE']).year
releases_df['MONTH'] = pd.DatetimeIndex(releases_df['RELEASE_DATE']).month
releases_df['EWE_TIMESTEP'] = releases_df['MONTH'] + ((releases_df['YEAR'] - start_year) * 12)
#print(releases_df['BIOMASS_MT'].sum().round())

#print(releases_df.columns)

# Apr 12 2022 - fix mean weight of releases
# something is wrong with avg_weight field. Coho are consistently released at
# 20 g but cross-checks indicate weight from spreadsheet far too low.
# (from EPAD Puntledge river releases, hatchery-reared, tossing avg_weight = 0)
coho_weight = 0.020 # kg
chin_weight = 0.0062 # kg
releases_df.loc[releases_df['SPECIES_NAME'] == 'Chinook', 'BIOMASS_MT2'] = releases_df['TOTRELEASE_NO'] * chin_weight * 0.001
releases_df.loc[releases_df['SPECIES_NAME'] == 'Coho', 'BIOMASS_MT2'] = releases_df['TOTRELEASE_NO'] * coho_weight * 0.001
releases_df.loc[(releases_df['SPECIES_NAME'] != 'Coho') &
                (releases_df['SPECIES_NAME'] != 'Chinook'), 'BIOMASS_MT2'] = releases_df['BIOMASS_MT']
print(releases_df['BIOMASS_MT2'].sum().round())
print(releases_df.head())

if aggregate_all_areas == "yes":
    releases_df2 = releases_df.drop(
        ['release_avg_date', 'FINAL_LAT', 'FINAL_LON', 'ROW_EWE', 'COL_EWE', 'SOURCE_ID', 'RELEASE_DATE'], axis=1)
    releases_df2 = releases_df2.groupby(['EWE_GROUP_CODE', 'SPECIES_NAME', 'EWE_TIMESTEP', 'YEAR',
                                         'MONTH']).agg('sum').reset_index()
    releases_df = releases_df2

if aggregate_to_level == "species":
    releases_df2 = releases_df.drop(['EWE_GROUP_CODE'], axis=1)
    releases_df2 = releases_df2.groupby(['SPECIES_NAME', 'EWE_TIMESTEP', 'YEAR',
                                         'MONTH']).agg('sum').reset_index()
    releases_df2['EWE_GROUP_CODE'] = releases_df2['SPECIES_NAME']
    releases_df2 = pd.merge(releases_df2, species_groupcodes, on=['EWE_GROUP_CODE'], how='left')

    releases_df = releases_df2.drop(['SPECIES_NAME'], axis=1)

releases_df['BIOMASS_MT2'] = releases_df['BIOMASS_MT2'] / study_area
releases_df['BIOMASS_MT'] = releases_df['BIOMASS_MT'] / study_area
print(releases_df)


# cross check for one year
print(releases_df.loc[releases_df['YEAR']==1980].groupby(['EWE_GROUP_CODE','YEAR']).sum().reset_index())
# chinook in 1980 should be 6930550, coho in 1980 should be 4552663 - GO 01152025

# round to 5 decimal places due to issues with floating point data
# storage causing rounding to not work so
# using decimal library https://stackoverflow.com/questions/56820/round-doesnt-seem-to-be-rounding-properly

releases_df['BIOMASS_MT2'] = releases_df['BIOMASS_MT2'].apply(lambda x: Decimal(str(x)).quantize(Decimal('.00001'), rounding=ROUND_UP))
releases_df['BIOMASS_MT'] = releases_df['BIOMASS_MT'].apply(lambda x: Decimal(str(x)).quantize(Decimal('.00001'), rounding=ROUND_UP))

# add dummy variable containing all timesteps
dummy = pd.Series(range(1,((end_year - start_year)*12)))
dummy_df = (dummy.to_frame())
dummy_df['EWE_TIMESTEP'] = dummy_df[0]
dummy_df['EWE_GROUP_CODE'] = "DUMMY"
dummy_df['YEAR'] = (dummy_df['EWE_TIMESTEP'] // 12)+start_year

dummy_df = dummy_df[['EWE_GROUP_CODE','EWE_TIMESTEP','YEAR']]

#releases_df = releases_df.append(dummy_df, ignore_index = True) # deprecated - 2025
releases_df = pd.concat([releases_df, dummy_df], ignore_index=True)
print(releases_df)


################ For Ecosim ########################
releasesEcosim = releases_df[['EWE_TIMESTEP', 'BIOMASS_MT2', 'TOTRELEASE_NO', 'EWE_GROUP_CODE', 'YEAR']]
releasesEcosim = releasesEcosim.fillna(0)

# sum by EWE_GROUP_CODE and timestep
releasesEcosim = releasesEcosim.rename(columns={'EWE_TIMESTEP': 'TIMESTEP','EWE_GROUP_CODE': 'TITLE'})

# for timestep = monthly
releasesEcosim_gp_mo = releasesEcosim.groupby(['TIMESTEP', 'TITLE', 'YEAR']).sum().reset_index()


# pivot wide
releasesEcosim_wide_mo = releasesEcosim_gp_mo.pivot_table(
        values=['BIOMASS_MT2', 'TOTRELEASE_NO'],
        index=['TIMESTEP', 'YEAR'],
        columns='TITLE',
        aggfunc=np.sum).reset_index()

# reset the multilevel index via hack
releasesEcosim_wide_mo['CHIN_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Chinook')].astype(float)
releasesEcosim_wide_mo['COHO_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Coho')].astype(float)
releasesEcosim_wide_mo['CHUM_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Chum')].astype(float)
releasesEcosim_wide_mo['SOCKEYE_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Sockeye')].astype(float)
releasesEcosim_wide_mo['PINK_H_MT'] = releasesEcosim_wide_mo[('BIOMASS_MT2', 'Pink')].astype(float)

releasesEcosim_wide_mo['CHIN_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Chinook')].astype(float)
releasesEcosim_wide_mo['COHO_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Coho')].astype(float)
releasesEcosim_wide_mo['CHUM_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Chum')].astype(float)
releasesEcosim_wide_mo['SOCKEYE_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Sockeye')].astype(float)
releasesEcosim_wide_mo['PINK_H_N'] = releasesEcosim_wide_mo[('TOTRELEASE_NO', 'Pink')].astype(float)

releasesEcosim_wide_mo['TIMESTEP'] = releasesEcosim_wide_mo[('TIMESTEP', '')].astype(float)
releasesEcosim_wide_mo['YEAR'] = releasesEcosim_wide_mo[('YEAR', '')].astype(float)
releasesEcosim_wide_mo = releasesEcosim_wide_mo[['YEAR',
                                                 'TIMESTEP',
                                                 'CHIN_H_MT',
                                                 'COHO_H_MT',
                                                 'CHUM_H_MT',
                                                 'SOCKEYE_H_MT',
                                                 'PINK_H_MT',
                                                 'CHIN_H_N',
                                                 'COHO_H_N',
                                                 'CHUM_H_N',
                                                 'SOCKEYE_H_N',
                                                 'PINK_H_N'
                                                 ]]
releasesEcosim_wide_mo.columns = [f"{x}_{y}" for x, y in releasesEcosim_wide_mo.columns.to_flat_index()]



# releasesEcosim_wide = releasesEcosim_wide.drop(columns=[('BIOMASS_MT2',   'DUMMY')])

# fill NaNs with zeros (required by ecosim)
releasesEcosim_wide_mo = releasesEcosim_wide_mo.fillna(0)

#releasesEcosim_wide_mo = pd.DataFrame(releasesEcosim_wide_mo.to_records())
print(releasesEcosim_wide_mo)
#print(releasesEcosim_wide_mo.columns)

releasesEcosim_wide_yr = releasesEcosim_wide_mo.groupby(['YEAR_']).agg({
    'CHIN_H_MT_': 'mean',
    'COHO_H_MT_': 'mean',
    'CHUM_H_MT_': 'mean',
    'SOCKEYE_H_MT_': 'mean',
    'PINK_H_MT_': 'mean',
    'CHIN_H_N_': 'sum',
    'COHO_H_N_': 'sum',
    'CHUM_H_N_': 'sum',
    'SOCKEYE_H_N_': 'sum',
    'PINK_H_N_': 'sum'
}).reset_index()

# #if aggregate_time == "year":
# # releasesEcosim_wide = releasesEcosim_wide.drop(columns="('TIMESTEP', '')", axis=1)
# #     releasesEcosim_wide['Chinook'] = releasesEcosim_wide["('BIOMASS_MT2', 'Chinook')"].astype(float)
# #     releasesEcosim_wide['Coho'] = releasesEcosim_wide["('BIOMASS_MT2', 'Coho')"].astype(float)
# #     releasesEcosim_wide = releasesEcosim_wide.groupby("('YEAR', '')").mean().reset_index()
#
# write to temp file
releasesEcosim_wide_yr.to_csv('../scratch/temp_yr.csv', index=True)
releasesEcosim_wide_mo.to_csv('../scratch/temp_mo.csv', index=True)

# this repeats same avg value each month, for silly workaround
repeated_yr_avg = pd.merge(releasesEcosim_wide_mo, releasesEcosim_wide_yr, on=['YEAR_'], how='left')
repeated_yr_avg = repeated_yr_avg[['YEAR_',
                                   'TIMESTEP_',
                                   'CHIN_H_MT__y',
                                   'COHO_H_MT__y',
                                   'CHUM_H_MT__y',
                                   'SOCKEYE_H_MT__y',
                                   'PINK_H_MT__y',
                                   'CHIN_H_N__y',
                                   'COHO_H_N__y',
                                   'CHUM_H_N__y',
                                   'SOCKEYE_H_N__y',
                                   'PINK_H_N__y'
                                   ]]

repeated_yr_avg = repeated_yr_avg.rename(columns={'TIMESTEP_': 'TIMESTEP',
                                                  'YEAR_': 'YEAR',
                                                  'CHIN_H_MT__y': 'CHIN_H_MT_KM2',
                                                  'COHO_H_MT__y': 'COHO_H_MT_KM2',
                                                  'CHUM_H_MT__y': 'CHUM_H_MT_KM2',
                                                  'PINK_H_MT__y': 'PINK_H_MT_KM2',
                                                  'SOCKEYE_H_MT__y': 'SOCKEYE_H_MT_KM2',
                                                  'CHIN_H_N__y': 'CHIN_H_N',
                                                  'COHO_H_N__y': 'COHO_H_N',
                                                  'CHUM_H_N__y': 'CHUM_H_N',
                                                  'PINK_H_N__y': 'PINK_H_N',
                                                  'SOCKEYE_H_N__y': 'SOCKEYE_H_N'
                                                  })
repeated_yr_avg.to_csv('../scratch/temp_yr_rep.csv', index=True)

# ===================================
# INSERT HEADER
# ===================================

# Title	Combined_GST_FR_Escape_RelB_NuSEDS	Chin_Hatch_RelB_CW	Chin_1stYrM_1_CW	Chin_1stYrM_2_CW	Chin_C_Rel_CW
# Weight	1	1	1	1	1
# Pool Code	14	18	16	15	14
# Type	0	0	5	5	61
# 1979	11.26655002	3.84	3.449022245	3.449022245	0.35
# 1980	11.07767237	6.93	3.021428984	3.021428984	0.371
# 1981	11.23108247	8.75	3.354206073	3.354206073	0.2533

# codes for 'type'
# relative biomass = 0
# absolute biomass = 1
# biomass forcing = -1
# fishing mortality = 4
# relative fishing mortality = 104
# total mortality = 5
# constant total mortality = -5 (forcing?)
# catches = 6
# catches forcing = -6
# relative catches = 61
# average weight = 7

import copy

#####################################################
# insert header to version with annual time increment
f = open('../scratch/temp_yr.csv', "r")
contents = f.readlines()
f.close()

line1 = contents[0].split(',')
line1[0] = 'Title'

line2 = copy.deepcopy(line1)
line2[0] = 'Weight'
i = 0
for line in line2:
    if i > 0:
        if i == (len(line2) - 1):
            line2[i] = '1\n'
        else:
            line2[i] = 1
    i += 1

line3 = copy.deepcopy(line1)
line3[0] = 'Type'
i = 0
for line in line3:
    if i > 0:
        if i == (len(line3) - 1):
            line3[i] = '-1\n'
        else:
            line3[i] = -1
    i += 1

line4 = copy.deepcopy(line1)
line4[0] = 'Timestep'
i = 0
for line in line4:
    if i > 0:
        if i == (len(line4) - 1):
            line4[i] = 'Interval\n'
        else:
            line4[i] = 'Interval'
    i += 1

s = ""
contents.insert(1, ','.join(str(line) for line in line1))
contents.insert(2, ','.join(str(line) for line in line2))
contents.insert(3, ','.join(str(line) for line in line3))
contents.insert(4, ','.join(str(line) for line in line4))

i = 0
with open('../data/forcing/HatcheryRel_Ecosim_TS_yr_2025.csv', 'w') as a_writer:
    for line in contents:
        if i > 0:
            a_writer.writelines(line)
        i += 1
a_writer.close()

#####################################################
# insert header to version with monthly time increment
f = open('../scratch/temp_yr_rep.csv', "r")
contents = f.readlines()
f.close()

line1 = contents[0].split(',')
line1[0] = 'Title'

line2 = copy.deepcopy(line1)
line2[0] = 'Weight'
i = 0
for line in line2:
    if i > 0:
        if i == (len(line2) - 1):
            line2[i] = '1\n'
        else:
            line2[i] = 1
    i += 1

line3 = copy.deepcopy(line1)
line3[0] = 'Type'
i = 0
for line in line3:
    if i > 0:
        if i == (len(line3) - 1):
            line3[i] = '-1\n'
        else:
            line3[i] = -1
    i += 1

line4 = copy.deepcopy(line1)
line4[0] = 'Timestep'
i = 0
for line in line4:
    if i > 0:
        if i == (len(line4) - 1):
            line4[i] = 'Interval\n'
        else:
            line4[i] = 'Interval'
    i += 1

s = ""
contents.insert(1, ','.join(str(line) for line in line1))
contents.insert(2, ','.join(str(line) for line in line2))
contents.insert(3, ','.join(str(line) for line in line3))
contents.insert(4, ','.join(str(line) for line in line4))

i = 0
with open('../data/forcing/HatcheryRel_Ecosim_TS_mo_2025.csv', 'w') as a_writer:
    for line in contents:
        if i > 0:
            a_writer.writelines(line)
        i += 1
a_writer.close()

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
# with open('../data/forcing/HatcheryRel_Ecosim_TS_Mo_2025.csv', 'w') as a_writer:
#     for line in contents:
#         if i > 0:
#             a_writer.writelines(line)
#         i += 1