
"""
ecosim config.py
Author: G Oldford, 2025

Description:
This script stores config params that the ecosim eval scripts have in common (single run)

Notes:
    - currently only set up for 3 day time step
"""

# -------------------------------------------
# General settings
# -------------------------------------------
SCENARIO = "SC123"

# ===== Paths =====
YEAR_START_FULLRUN = 1978
YEAR_END_FULLRUN = 2018
TIMESTEP_DAYS = 3
MAX_TIMESTEPS = 120 # per year (will lump >120 with 120 if 3-day)

# Base directories
ECOPATH_F_NM = "ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v13_BigPC"
ECOSIM_F_RAW_SINGLERUN = "biomass_monthly.csv"
ECOSIM_F_RAW_HEADERN = 14
ECOSIM_RAW_DIR = f"C://Users//Greig//Documents//EwE output//{ECOPATH_F_NM}//ecosim_Ecosim_{SCENARIO}//{ECOSIM_F_RAW_SINGLERUN}"
OUTPUT_DIR_EVAL = "..//..//data//evaluation//"
OUTPUT_DIR_FIGS = "..//..//figs//"
ECOSIM_F_PREPPED_SINGLERUN = f"{OUTPUT_DIR_EVAL}//ecosim_{SCENARIO}_onerun_B_dates_seasons.csv"

# -------------------------------------------
# evaluation 1 - seasonal average B
# -------------------------------------------

START_DATA_PHYT_SEAS = '2015-01-01'
END_DATE_PHYT_SEAS = '2018-12-31'

# -------------------------------------------
# evaluation 2 - bloom timing
# -------------------------------------------

# Biomass columns for bloom detection (output from model is named by model group number)
BIOMASS_COLS_SATELLITE = ['17', '18']
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
# evaluation 3 - nutrients and annual pattern
# -------------------------------------------

C_TO_N_RATIO = 106 / 16  # molar Redfield ratio
N_FREE_AVG = 18 # see ecosim_data_prep_2_nutrients.py
N_FREE_MIN = 7
N_FREE_DRAWDOWN = (N_FREE_AVG - N_FREE_MIN)
N_BOUND_PROP = N_FREE_DRAWDOWN / N_FREE_AVG
PP_COLUMNS = [17, 18, 19]

# FIX SO IT'S READING FROM FIRST TIME STEP
INIT_TOT_PP_B = 3.7  # g C m^-2, for initial phytoplankton biomass only