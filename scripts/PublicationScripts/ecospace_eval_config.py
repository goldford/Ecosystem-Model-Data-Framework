
"""
ecospace_eval_config.py
Author: G Oldford, 2025

Description:
This script stores config params that the ecospace eval scripts have in common (single run)

Notes:
    - currently only set up for 3 day time step
"""

import pandas as pd # would rather not do imports here
import os
from typing import Dict, List

# -------------------------------------------
# General settings
# -------------------------------------------
ECOSPACE_SC = "SC138_w2"
ECOSPACE_SC_FULL = "SC138_w2"
ECOPATH_F_NM = "ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v15_BigPC"
ECOSPACE_RAW_DIR = f"C://Users//Greig//Documents//EwE output//{ECOPATH_F_NM}//Ecospace_{ECOSPACE_SC_FULL}//asc//"

ECOSPACE_RN_STR_YR = 1978
ECOSPACE_RN_END_YR = 2018
ECOSPACE_RN_STR_MO = 1
ECOSPACE_RN_STR_DA = 2
ECOSPACE_RN_END_MO = 12
ECOSPACE_RN_END_DA = 30

NC_FILENAME = ECOSPACE_SC + "_" + str(ECOSPACE_RN_STR_YR) + "-" + str(ECOSPACE_RN_END_YR) + ".nc"

# NC_PATH_OUT = "..//..//data//ecospace_out//"
NC_PATH_OUT = r"C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/ECOSPACE_OUT/"
NEMO_EWE_CSV = "..//..//data//basemap//Ecospace_grid_20210208_rowscols.csv"

DO_NC_CROSSCHECK = False # can be memory intesnive to crosscheck NC with ASC after (they always match anyway)

EWE_ROWS = 151
EWE_COLS = 93
SKIPROWS = 6 # header

EVALOUT_P = "C://Users//Greig/Documents//github//Ecosystem-Model-Data-Framework//data//evaluation//"

# -------------------------------------------
# Prep 2 - mcewan seasonal eval
# -------------------------------------------

MW_SHOW_PLTS = False # causes halt of pipeline - optional
MW_SEASON_DEF: Dict[str, List[int]] = {
    "Winter": [11,12, 1, 2],
    "Spring": [3, 4, 5],
    "Summer": [6, 7, 8],
    # "Fall":   [9, 10],  # not used for McEwan comparison
}
MW_GROUP_MAP = {
    "Diatoms": "PP1-DIA",
    "Nano":    "PP2-NAN",
    "Other":   "PP3-PIC",  # adjust if your "Other" aggregates differently
}
MW_GROUP_COLORS = {
    "Diatoms": "#1f77b4",
    "Nano":    "#ff7f0e",
    "Other":   "#2ca02c",
}
MW_DOMAIN_FP = "..//..//data/evaluation//analysis_domains_mcewan.yml"
MW_STATS_OUT = "..//..//data/evaluation"
MW_FIGS_OUT  = "..//..//figs"

MW_DOMAINS = {
    "NSoG": "NSoG",   # Northern Strait of Georgia
    "CSoG": "CSoG",   # Central/Southern Strait of Georgia (rename if needed)
}
MW_OBS_MCEWAN = {
    "NSoG": {
        'Diatoms': {'Winter': 0.77, 'Spring': 6.00, 'Summer': 1.16},
        'Nano':    {'Winter': 1.38, 'Spring': 1.85, 'Summer': 1.90},
        'Other':   {'Winter': 0.12, 'Spring': 0.26, 'Summer': 0.34},
    },
    "CSoG": {
        'Diatoms': {'Winter': 0.65, 'Spring': 8.06, 'Summer': 1.34},
        'Nano':    {'Winter': 1.57, 'Spring': 1.58, 'Summer': 1.53},
        'Other':   {'Winter': 0.12, 'Spring': 0.19, 'Summer': 0.74},
    }
}
# Spatial reduction across the mask: "mean" or "median"
MW_SPATIAL_REDUCTION = "mean"
MW_START_DATE = '2015-01-01'
MW_END_DATE   = '2018-12-31'




# -------------------------------------------
# Prep 2 - nemcek dataset eval
# -------------------------------------------

NM_SHOW_PLTS = False # causes halt of pipeline - optional

SSC_MO_F = "SalishSeaCast_biology_2008.nc"
SSC_GRD_F = "ubcSSnBathymetryV21-08_a29d_efc9_4047.nc"
NEMCEK_F = "Nemcek_Supp_Data.csv"
NEMCEK_MATCHED_F = f"Nemcek_matched_to_ecospace_{ECOSPACE_SC}.csv"
ECOSPACE_MAP_F = "Ecospace_grid_20210208_rowscols.csv"

# File paths
SSC_P = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/SalishSeaCast_BioFields_2008-2018/ORIGINAL"
ECOSPACE_P = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/ECOSPACE_OUT"
NEMCEK_P = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/Phytoplankton Salish Sea Nemcek2023 2015-2019/MODIFIED"

ECOSPACE_MAP_P = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap"
DOMAIN_P = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
DOMAIN_F = "analysis_domains_jarnikova.yml"

NEMCEK_MATCHED_PF = os.path.join(EVALOUT_P, NEMCEK_MATCHED_F)

NM_LOG_TRANSFORM_OBS = False     # Apply log transform to observational values
NM_LOG_TRANSFORM_MODEL = False  # Apply log transform to model values only for anomaly calc
NM_EPSILON = 1e-4        # Small constant to avoid log(0)
NM_DO_EXPLOR_PLOTS = False # exploratory histograms, data summary
NM_APPLY_NEMCEK_SEASON = False


# -------------------------------------------
# Prep 3 - qu39 dataset eval
# -------------------------------------------

QU39_SHOW_PLTS = False # causes halt of pipeline - optional
QU39_DO_ECOSPACE_MATCHING = True       # Perform Ecospace matching
QU39_DOWNLOAD_SSC = False              # Download SSC ERDDAP files
QU39_MATCH_SSC = True                  # Match existing SSC files

QU39_IN_P = r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\Phyto Concent del Bel Belluz 2024 2016 - 2019\MODIFIED\QU39_joined.csv'

QU39_OUT_P = r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\Phyto Concent del Bel Belluz 2024 2016 - 2019\MODIFIED'
QU39_MATCHED_OUT_F = 'QU39_joined_matchtoEcospace' # gets tagged w scen later

SSC_DDWNLD_P = os.path.join(SSC_P, 'matched_QU39_dwnld')
SSC_BATHY_F = "ubcSSnBathymetryV21-08_a29d_efc9_4047.nc"
SSC_MESH3D_F = "ubcSSn3DMeshMaskV17-02.nc"
SSC_MESHM_2D_F = "ubcSSn2DMeshMaskV17-02.nc"

# ----- QU39 sampling site) -----
QU39_LAT = 50.0307
QU39_lON = -125.0992

QU39_INIT_DT_PLACEHOLDER = pd.to_datetime('1900-04-23 23:42')
QU39_FILL_VAL = -999.9

QU39_SSC_ECOSPACE_PF = (
    r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton'
    r'\Phyto Concent del Bel Belluz 2024 2016 - 2019\MODIFIED'
    r'\QU39_joined_matchtoEcospace_SSC_' + ECOSPACE_SC + '.csv'
)

QU39_INCLUDE_SSC = True
QU39_FIGS_OUT_P = os.path.normpath("..//..//figs//") + os.sep
QU39_STATS_OUT_P = os.path.normpath("..//..//data//evaluation//") + os.sep
QU39_FILL_VALUE = -999

QU39_DO_HIST = False


# -------------------------------------------
# Prep 4 - B Timing
# -------------------------------------------
BT_DOMAIN_CONFIG_PATH = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
BT_DOMAIN_FILE = "analysis_domains_suchy.yml"
BT_SAT_MASK_PF = r'../../data/evaluation/suchy_ecospace_mask.nc'
BT_RECOMPUTE_BLOOM_TIMING_SAT = True  # Set to True to force recomputation as needed, saves time
BT_RECOMPUTE_BLOOM_TIMING_C09 = True
BT_START_YEAR = 1980 # analysis years (exclude spinup?)
BT_END_YEAR = 2018

BT_CSV_SUCHY_PF = f"..//..//data//evaluation//ecospace_bloom_timing_SSoG_{ECOSPACE_SC}.csv"
BT_CSV_ALLEN_PF = f"..//..//data//evaluation//ecospace_bloom_timing_CSoG_{ECOSPACE_SC}.csv"

# if multiple vars listed here, it will sum across them when computing anomalies, bloom timing etc!
#OK with run 96 this is first time model fit has been okay with all pp groups
# VARIABLES_TO_ANALYZE_SAT = ["PP1-DIA"] #for 88_2 - pp1-dai only, annual, keep janfeb
BT_VARS_TO_ANALYZE_SAT = ["PP1-DIA", "PP2-NAN", "PP3-PIC"]
BT_VARS_TO_ANALYZE_C09 = ["PP1-DIA"]
BT_ANNUAL_AVG_METHOD_SAT = "annual" # should average bloom compared against be from all years, or just one year
BT_ANNUAL_AVG_METHOD_C09 = "all" # annual or all
# Bloom detection method
BT_LOG_TRANSFORM = True
BT_MEAN_OR_MEDIAN = "median"
BT_THRESHOLD_FACTOR = 1.05
BT_SUB_THRESHOLD_FACTOR = 0.7

# use mask c09 instead of pnt?
BT_CREATE_MASKS = False  # Set to True to regenerate masks
BT_USE_SAT_MASK_CO9 = False
BT_USEC09_MASK_FOR_SAT = False

BT_MASK_REGNM = "SGC3" # SGC2 is more accurate to map in suchy pub, but in 2005 clouds mask northern portion, affecting accuracy so SGC3 trims north
BT_DO_NUTRIENTS = False # another script does this now (#9?)
# OVERRIDE_REDFIELD = True # added by GO to help eval 2025-06-03

BT_EXCLUDE_DEC_JAN_SAT = False # not hooked up ?
BT_EXCLUDE_DEC_JAN_C09 = False

BT_MIN_Y_TICK = 38







# # ===== Paths =====
# YEAR_START_FULLRUN = 1978
# YEAR_END_FULLRUN = 2018
# TIMESTEP_DAYS = 3
# MAX_TIMESTEPS = 120 # per year (will lump >120 with 120 if 3-day)


#
#
# # Base directories
# ECOPATH_F_NM = "ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v14_BigPC_CW"
# # ECOPATH_F_NM = "ECOSPACE_LTL_Carb_3day_ewe6_7_19295_v13_Jun30CopyVCCW_Laptop"
# ECOSIM_F_RAW_SINGLERUN = "biomass_monthly.csv"
# ECOSIM_F_RAW_HEADERN = 14
# ECOSIM_RAW_DIR = f"C://Users//Greig//Documents//EwE output//{ECOPATH_F_NM}//ecosim_Ecosim_{SCENARIO}//{ECOSIM_F_RAW_SINGLERUN}"
# OUTPUT_DIR_EVAL = "..//..//data//evaluation//"
# OUTPUT_DIR_FIGS = "..//..//figs//"
# ECOSIM_F_PREPPED_SINGLERUN = f"{OUTPUT_DIR_EVAL}//ecosim_{SCENARIO}_onerun_B_dates_seasons.csv"
# NUTRIENTS_F_PREPPED = f"{OUTPUT_DIR_EVAL}//nutrients_ios_csop_combined_sampled.csv"
#
#
# # -------------------------------------------
# # evaluation 1 - rel PP multiplier calc
# #  (corrects for problem in Ecosim)
# # -------------------------------------------
#
# START_SPINUP = '1978-01-01'
# END_SPINUP = '1995-12-31'
# START_FULL = '1996-01-01'
# END_FULL = '2018-12-31'
# GROUPS_relPP = [1, 2, 3, 4, 5,
#                 6, 7, 8, 9, 10,
#                 11, 12, 13, 14, 15,
#                 16, 17, 18, 19, 20]
# # must be updated when Ecopath changed!
# # (note Ecosim out time = 0 does not include Ecopath base values (it's one time-step in))
# GROUPS_ECOPATH_B = {
#     1: 0.0046,
#     2: 0.00018,
#     3: 0.28,
#     4: 0.03,
#     5: 0.6,
#     6: 0.52,
#     7: 0.09,
#     8: 0.64,
#     9: 0.5,
#     10: 0.003,
#     11: 0.03,
#     12: 0.26,
#     13: 0.04,
#     14: 1.1,
#     15: 0.87,
#     16: 0.46,
#     17: 2.31, # PP1
#     18: 1.5,  # PP2
#     19: 0.35, # PP3
#     20: 0.17,
#     21: 4.3, #DE1
#     22: 10.1 #DE2
# }
#
#
#
# # -------------------------------------------
# # evaluation 2 - seasonal average B
# # -------------------------------------------
#
# START_DATA_PHYT_SEAS = '2015-01-01'
# END_DATE_PHYT_SEAS = '2018-12-31'
#
#
#
# # -------------------------------------------
# # evaluation 3 - nutrients and annual pattern
# # -------------------------------------------
#
# # C_TO_N_RATIO = 106 / 16  # molar Redfield ratio
# N_FREE_AVG_INIT = 18 # see ecosim_data_prep_2_nutrients.py, depth average N, umol /L
# P_FREE_INIT = 0.78# must be updated when Ecosim change is made!
# # N_BOUND_GROUPS = [1, 2, 3, 4, 5,
# #                   6, 7, 8, 9, 10,
# #                   11, 12, 13, 14, 15,
# #                   16, 17, 18, 19, 20]
# N_BOUND_GROUPS = [17, 18, 19]
#
# OBS_AVG_TYPE = "mean" # mean, median
#
# TOT_B_INIT = 0
# for group_col in N_BOUND_GROUPS:
#     INIT_B = GROUPS_ECOPATH_B[group_col]
#     TOT_B_INIT += INIT_B
#
# N_B_INIT = N_FREE_AVG_INIT *  (1 - P_FREE_INIT) / P_FREE_INIT # inferred absolute amount bound in biomass in nitrogen units (see text)
#
# # if there is seasonal 'nutrient loading' set up as penalty for flux reduced in spring summer
# # Account (crudely) for seasonally reduced upward nutrient flux to the surface layer
# # (e.g., stratified spring/summer) by scaling *model-derived free nutrient* (or, equivalently,
# # model-derived nutrient drawdown) before climatology and comparison to observations
# USE_N_MULT = True # MAKE SURE FALSE IF ALSO NO LOADING FUNC APPLIED IN ECOSIM
# # EWE_NUTR_LOADING_FILE = "Forcings Data Prep//Nutrient_loading_using_RDRS.csv"
# EWE_NUTR_LOADING_FILE = "..//..///data//forcing//ECOSIM_in_3day_vars_1980-2018_fromASC_202506//ECOSIM_in_NEMO_varmixing_m_stdfilter_1980-2018.csv"
# N_MULT_TYPE = "3day" # seasonal, monthly, 3day
#
# SEASONAL_N_FLUX_MULT = {
#     "Winter": 1.0,
#     "Spring": 0.8,
#     "Summer": 0.6,
#     "Fall": 0.9
# }
#
#
# # nutrient loading function values
# # Scenarios prio to 129
# # MONTHLY_N_FLUX_MULT = {
# #     1: 1.0, 2: 1.0, 3: 0.92, 4: 0.82, 5: 0.67, 6: 0.5,
# #     7: 0.42, 8: 0.45, 9: 0.74, 10: 0.94, 11: 1.0, 12: 1.0
# # }
#
# # scenario 129
# MONTHLY_N_FLUX_MULT = {
#     1: 1.37, 2: 1.36, 3: 1.28, 4: 1.0, 5: 0.65, 6: 0.55,
#     7: 0.53, 8: 0.47, 9: 0.75, 10: 0.95, 11: 1.4, 12: 1.55
# }
#
#
# # MONTHLY_N_FLUX_MULT = {
# #     1: 1.0, 2: 1.0, 3: 1, 4: 1, 5: 1, 6: 1,
# #     7: 1, 8: 1, 9: 1, 10: 1, 11: 1.0, 12: 1.0
# # }
#
# SEASON_MAP = {
#     1:  "Winter", 2: "Winter", 12: "Winter",
#     3:  "Spring", 4: "Spring", 5:  "Spring",
#     6:  "Summer", 7: "Summer", 8:  "Summer",
#     9:  "Fall",   10: "Fall",   11: "Winter",
# }
#
# ECOSIM_F_W_NUTRIENTS =  f"{OUTPUT_DIR_EVAL}//ecosim_{SCENARIO}_onerun_B_dates_seasons_nutrients.csv"
# print(f"Nitrogen bound at initialisation: {N_B_INIT:.3f}")
#
#
#
# # -------------------------------------------
# # evaluation 4 - bloom timing
# # -------------------------------------------
#
# START_FULL_BLM = '1997-01-01' # 1996 is off because of imprecise spinup data prep (bloom is always day 20
# END_FULL_BLM = '2018-12-31'
#
# # Biomass columns for bloom detection (output from model is named by model group number)
# BIOMASS_COLS_SATELLITE = ['17', '18', '19']
# BIOMASS_COLS_C09 = ['17']
#
# # adds a col before doing eval (sum across groups above, optionally)
# TOTAL_BIOMASS_COL_SATELLITE = "Biomass_Total_Satellite"
# TOTAL_BIOMASS_COL_C09 = "Biomass_Total_C09"
#
# # Thresholds and bloom detection params
# THRESHOLD_FACTOR = 1.05
# SUB_THRESHOLD_FACTOR = 0.78
# LOG_TRANSFORM = True
# MEAN_OR_MEDIAN = "median"
# NUTRIENT_DRAWDOWN_FRAC = 0.6 # for nutrient definition of bloom
#
#
# # -------------------------------------------
# # evaluation 5 - zooplankton eval
# # -------------------------------------------
#
# Z_F_SEAS = "Zoopl_SofG_1996-2018_df_summary.csv" # this is output by long R script
# Z_F_TOWLEV = "Zooplankton_B_C_gm2_EWEMODELGRP_Wide_NEMO3daymatch.csv" # this is output by short one
# Z_P_PREPPED = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/4. Zooplankton/Zoop_Perryetal_2021/MODIFIED"
#
# # Mapping from obs column names to numeric column indices in the model output
# Z_GROUP_MAP = {
#     'ZF1-ICT': 4,
#     'ZC1-EUP': 5,
#     'ZC2-AMP': 6,
#     'ZC3-DEC': 7,
#     'ZC4-CLG': 8,
#     'ZC5-CSM': 9,
#     'ZS1-JEL': 10,
#     'ZS2-CTH': 11,
#     'ZS3-CHA': 12,
#     'ZS4-LAR': 14,
# }
#
# # USER SETTING: choose one season and plot type ('bar' or 'line')
# ZP_SEASON_CHOICE = 'Spring'  # e.g., 'Winter','Spring','Summer','Fall'
# ZP_PLOT_TYPE = 'bar'     # 'bar' or 'line'
# # USER SETTING: toggle friendly labels
# ZP_USE_FRIENDLY_LABELS = True  # set False to use cryptic codes
# ZP_FRIENDLY_MAP_ZC = {
#     'ZF1-ICT': 'Ichthyo',
#     'ZC1-EUP': 'Euphausiids',
#     'ZC2-AMP': 'Amphipods',
#     'ZC3-DEC': 'Decapods',
#     'ZC4-CLG': 'Lg. Calanoid Copepods',
#     'ZC5-CSM': 'Other Copepods',
#     'ZS1-JEL': 'Jellyfish',
#     'ZS2-CTH': 'Ctenophores',
#     'ZS3-CHA': 'Chaetognaths',
#     'ZS4-LAR': 'Larvaceans',
#     'ZF1-ICH': 'Ichthyoplankton',
#     'misc': 'Other',
#     'Total': 'Total'
# }
# ZP_LOG_TRANSFORM = True
# ZP_YEAR_START = 2000
# ZP_YEAR_END   = 2018
#
# ZP_SHOW_CNTS = False