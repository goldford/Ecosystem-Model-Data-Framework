
"""
ecosim config.py
Author: G Oldford, 2025

Description:
This script stores config params that the ecosim eval scripts have in common (single run)

Notes:
    - currently only set up for 3 day time step
"""

import pandas as pd # would rather not do imports here

# -------------------------------------------
# General settings
# -------------------------------------------
SCENARIO = "SC126"

# ===== Paths =====
YEAR_START_FULLRUN = 1978
YEAR_END_FULLRUN = 2018
TIMESTEP_DAYS = 3
MAX_TIMESTEPS = 120 # per year (will lump >120 with 120 if 3-day)

# Base directories
ECOPATH_F_NM = "ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v14_BigPC"
# ECOPATH_F_NM = "ECOSPACE_LTL_Carb_3day_ewe6_7_19295_v13_Jun30CopyVCCW_Laptop"
ECOSIM_F_RAW_SINGLERUN = "biomass_monthly.csv"
ECOSIM_F_RAW_HEADERN = 14
ECOSIM_RAW_DIR = f"C://Users//Greig//Documents//EwE output//{ECOPATH_F_NM}//ecosim_Ecosim_{SCENARIO}//{ECOSIM_F_RAW_SINGLERUN}"
OUTPUT_DIR_EVAL = "..//..//data//evaluation//"
OUTPUT_DIR_FIGS = "..//..//figs//"
ECOSIM_F_PREPPED_SINGLERUN = f"{OUTPUT_DIR_EVAL}//ecosim_{SCENARIO}_onerun_B_dates_seasons.csv"
NUTRIENTS_F_PREPPED = f"{OUTPUT_DIR_EVAL}//nutrients_ios_csop_combined_sampled.csv"


# -------------------------------------------
# evaluation 1 - rel PP multiplier calc
#  (corrects for problem in Ecosim)
# -------------------------------------------

START_SPINUP = '1978-01-01'
END_SPINUP = '1995-12-31'
START_FULL = '1996-01-01'
END_FULL = '2018-12-31'
GROUPS_relPP = [14, 15, 16, 17, 18, 19, 20]
# must be updated when Ecopath changed!
# (note Ecosim out time = 0 does not include Ecopath base values (it's one time-step in))
GROUPS_ECOPATH_B = {
    1: 0.0046,
    2: 0.00018,
    3: 0.28,
    4: 0.03,
    5: 0.6,
    6: 0.52,
    7: 0.09,
    8: 0.64,
    9: 0.5,
    10: 0.003,
    11: 0.03,
    12: 0.26,
    13: 0.04,
    14: 1.1,
    15: 0.87,
    16: 0.46,
    17: 2.31, # PP1
    18: 1.5,  # PP2
    19: 0.35, # PP3
    20: 0.17
}


# -------------------------------------------
# evaluation 2 - seasonal average B
# -------------------------------------------

START_DATA_PHYT_SEAS = '2015-01-01'
END_DATE_PHYT_SEAS = '2018-12-31'


# -------------------------------------------
# evaluation 3 - nutrients and annual pattern
# -------------------------------------------

# C_TO_N_RATIO = 106 / 16  # molar Redfield ratio
N_FREE_AVG_INIT = 18 # see ecosim_data_prep_2_nutrients.py, depth average N, umol /L
P_FREE_INIT = 0.55 # must be updated when Ecosim change is made!
N_BOUND_GROUPS = [1, 2, 3, 4, 5,
                  6, 7, 8, 9, 10,
                  11, 12, 13, 14, 15,
                  16, 17, 18, 19, 20]

TOT_B_INIT = 0
for group_col in N_BOUND_GROUPS:
    INIT_B = GROUPS_ECOPATH_B[group_col]
    TOT_B_INIT += INIT_B

N_B_INIT = N_FREE_AVG_INIT *  (1 - P_FREE_INIT) / P_FREE_INIT # inferred absolute amount bound in biomass in nitrogen units (see text)

# if there is seasonal 'nutrient loading' set up as penalty for flux reduced in spring summer
# Account (crudely) for seasonally reduced upward nutrient flux to the surface layer
# (e.g., stratified spring/summer) by scaling *model-derived free nutrient* (or, equivalently,
# model-derived nutrient drawdown) before climatology and comparison to observations
USE_N_MULT = True # MAKE SURE FALSE IF ALSO NO LOADING FUNC APPLIED IN ECOSIM
N_MULT_TYPE = "monthly" # seasonal, monthly
SEASONAL_N_FLUX_MULT = {
    "Winter": 1.0,
    "Spring": 0.8,
    "Summer": 0.6,
    "Fall": 0.9
}

MONTHLY_N_FLUX_MULT = {
    1: 1.0, 2: 1.0, 3: 0.92, 4: 0.82, 5: 0.67, 6: 0.5,
    7: 0.42, 8: 0.45, 9: 0.74, 10: 0.94, 11: 1.0, 12: 1.0
}

SEASON_MAP = {
    1:  "Winter", 2: "Winter", 12: "Winter",
    3:  "Spring", 4: "Spring", 5:  "Spring",
    6:  "Summer", 7: "Summer", 8:  "Summer",
    9:  "Fall",   10: "Fall",   11: "Winter",
}

print(f"Nitrogen bound at initialisation: {N_B_INIT:.3f}")




# # N_FREE_MIN = 7 # the lowest the depth int N falls
# # N_FREE_DRAWDOWN = (N_FREE_AVG - N_FREE_MIN)
#
# PP_COLUMNS = ['17', '18', '19'] # for calc of N
#
# # FIX SO IT'S READING FROM FIRST TIME STEP
#
# # let's move all this to nutrients script
# # df = pd.read_csv(ECOSIM_F_PREPPED_SINGLERUN, skiprows=0)
# INITIAL_BIOMASS = 8.5  # g C m^-2, for initial phytoplankton biomass only
#
# MAX_BIOMASS = 15.2 # under bloom conditions
# # New input assumptions
# INITIAL_NUTRIENT_CONC_FREE = 18  # umol N L^-1 annual avg
# MIN_N_CONC = 7        # μmol N L^-1, observed minimum concentration
#
# # Calculation of 'Free' Nutrients (THIS SHOULD MATCH ECOSIM SCENARIO)
# B_MAX_INCREASE_PROP = (MAX_BIOMASS - INITIAL_BIOMASS) / INITIAL_BIOMASS
# N_MAX_DECREASE_PROP = 1 - (MIN_N_CONC / INITIAL_NUTRIENT_CONC_FREE)
# N_EQ_TO_B_EQ = B_MAX_INCREASE_PROP / N_MAX_DECREASE_PROP
#
# # Calculate total nitrogen in the system
# # T = N_eq + P_eq
# # Let P_eq = 1 unit (arbitrary), then N_eq = N_eq_to_P_eq * P_eq
# P_eq = 1  # define as 1 unit for proportional calculation
# N_eq = N_EQ_TO_B_EQ * P_eq
# T = N_eq + P_eq
# # Calculate the proportion that is 'free' (i.e. as dissolved nutrients)
# PROP_N_FREE = N_eq / T

# -------------------------------------------
# evaluation 4 - bloom timing
# -------------------------------------------

START_FULL_BLM = '1997-01-01' # 1996 is off because of imprecise spinup data prep (bloom is always day 20
END_FULL_BLM = '2018-12-31'

# Biomass columns for bloom detection (output from model is named by model group number)
BIOMASS_COLS_SATELLITE = ['17', '18', '19'] # ['17']
BIOMASS_COLS_C09 = ['17']

# adds a col before doing eval (sum across groups above, optionally)
TOTAL_BIOMASS_COL_SATELLITE = "Biomass_Total_Satellite"
TOTAL_BIOMASS_COL_C09 = "Biomass_Total_C09"

# Thresholds and bloom detection params
THRESHOLD_FACTOR = 1.05
SUB_THRESHOLD_FACTOR = 0.7
LOG_TRANSFORM = True
MEAN_OR_MEDIAN = "median"

# -------------------------------------------
# evaluation 5 - zooplankton eval
# -------------------------------------------

# Mapping from obs column names to numeric column indices in the model output
Z_GROUP_MAP = {
    'ZF1-ICT': 4,
    'ZC1-EUP': 5,
    'ZC2-AMP': 6,
    'ZC3-DEC': 7,
    'ZC4-CLG': 8,
    'ZC5-CSM': 9,
    'ZS1-JEL': 10,
    'ZS2-CTH': 11,
    'ZS3-CHA': 12,
    'ZS4-LAR': 14,
}
