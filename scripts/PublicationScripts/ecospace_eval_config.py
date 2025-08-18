
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
ECOSPACE_SC = "SC141_1"
ECOSPACE_SC_FULL = "SC141_1"
ECOPATH_F_NM = "ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v16_BigPC"
ECOSPACE_RAW_DIR = f"C://Users//Greig//Documents//EwE output//{ECOPATH_F_NM}//Ecospace_{ECOSPACE_SC_FULL}//asc//"

FIGS_P = "../../figs"
EVALOUT_P = "../../data/evaluation"

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


# -------------------------------------------
# Prep 2 - mcewan seasonal eval
# -------------------------------------------

MW_SHOW_PLTS = True # causes halt of pipeline - optional
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
MW_STATS_OUT = EVALOUT_P
MW_FIGS_OUT  = FIGS_P

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

NM_SHOW_PLTS = True # causes halt of pipeline - optional

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
BT_ANNUAL_AVG_METHOD_SAT = "annual" # should average bloom compared against be from all years, or just one year
BT_VARS_TO_ANALYZE_C09 = ["PP1-DIA"]
BT_ANNUAL_AVG_METHOD_C09 = "annual" # annual or all (this is always annual for C09 PCT MAX method)

# Bloom detection method
BT_LOG_TRANSFORM_SAT = True
BT_MEAN_OR_MEDIAN_SAT = "median"
BT_THRESHOLD_FACTOR_SAT = 1.05
BT_SUB_THRESHOLD_FACTOR_SAT = 0.7

BT_LOG_TRANSFORM_C09 = False
BT_USE_PCT_MAX_C09 = False # new method added 2025-9-12 GO
BT_PCT_MAX_C09 = 0.9
BT_PCT_MAX_WINDOW_DAYS_C09 = 6

# use mask c09 instead of pnt?
BT_CREATE_MASKS = True  # Set to True to regenerate masks
BT_USE_SAT_MASK_CO9 = False
BT_USEC09_MASK_FOR_SAT = False

BT_MASK_REGNM = "SGC3" # SGC2 is more accurate to map in suchy pub, but in 2005 clouds mask northern portion, affecting accuracy so SGC3 trims north
BT_DO_NUTRIENTS = False # another script does this now (#9?)
# OVERRIDE_REDFIELD = True # added by GO to help eval 2025-06-03

BT_EXCLUDE_DEC_JAN_SAT = False # not hooked up ?
BT_EXCLUDE_DEC_JAN_C09 = False

BT_MIN_Y_TICK = 38


# -------------------------------------------
# evaluation 5 - zooplankton eval
# -------------------------------------------

Z_F_SEAS = "Zoopl_SofG_1996-2018_df_summary.csv" # this is output by long R script
Z_F_TOWLEV = "Zooplankton_B_C_gm2_EWEMODELGRP_Wide_NEMO3daymatch.csv" # this is output by short one
Z_P_PREPPED = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/4. Zooplankton/Zoop_Perryetal_2021/MODIFIED"




# -------------------------------------------
# evaluation 6 - taylor plots
# -------------------------------------------
TY_STATS_FPAT_SPC_FMT = "ecospace_bloom_timing_stats_*.csv"
TY_STATS_FPAT_SIM_FMT = "ecosim_bloom_eval_stats_*.csv"

TY_STATS_F_SPC_FMT = r"^ecospace_bloom_timing_stats_(.+)\.csv$"
TY_STATS_F_SIM_FMT = r"^ecosim_bloom_eval_stats_(.+)\.csv$"


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