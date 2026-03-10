
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

SCENARIO = "SC214"
ECOPATH_F_NM = "LTL_Carb_3day_ewe6_7_19295_v18_BigPC_ECOSPACEPARAMZ"

# ===== Paths =====
YEAR_START_FULLRUN = 1978
YEAR_END_FULLRUN = 2018
TIMESTEP_DAYS = 3
MAX_TIMESTEPS = 120 # per year (will lump >120 with 120 if 3-day)


# ECOPATH_F_NM = "ECOSPACE_LTL_Carb_3day_ewe6_7_19295_v13_Jun30CopyVCCW_Laptop"
ECOSIM_F_RAW_SINGLERUN = "biomass_monthly.csv"
ECOSIM_F_RAW_HEADERN = 14
ECOSIM_RAW_DIR = f"C://Users//Greig//Documents//EwE output//{ECOPATH_F_NM}//ecosim_Ecosim_{SCENARIO}//{ECOSIM_F_RAW_SINGLERUN}"
OUTPUT_DIR_EVAL = "..//..//data//evaluation//"
OUTPUT_DIR_FIGS = "..//..//figs//"
ECOSIM_F_PREPPED_SINGLERUN = f"{OUTPUT_DIR_EVAL}//ecosim_{SCENARIO}_onerun_B_dates_seasons.csv"
NUTRIENTS_F_PREPPED = f"{OUTPUT_DIR_EVAL}//nutrients_ios_csop_combined_sampled.csv"

MATCH_MCEWAN_SEAS = False # lumps Feb in with Winter, lumps June in with Spring

# -------------------------------------------
# evaluation 1 - rel PP multiplier calc
#  (corrects for problem in Ecosim)
# -------------------------------------------

START_SPINUP = '1978-01-01'
END_SPINUP = '1995-12-31'
START_FULL = '1996-01-01'
END_FULL = '2018-12-31'
GROUPS_relPP = [1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                11, 12, 13, 14, 15,
                16, 17, 18, 19, 20]

# must be updated when Ecopath changed!
# (note Ecosim out time = 0 does not include Ecopath base values (it's one time-step in))
# added same PP mults
GROUPS_ECOPATH_B = {
    1: 0.0046,  # NK1-COH
    2: 0.00018, # NK2-CHI
    3: 0.28,    # NK3-FOR
    4: 0.03,    # ZF1-ICT
    5: 0.6,     # ZC1-EUP
    6: 0.52,    # ZC2-AMP
    7: 0.09,    # ZC3-DEC
    8: 0.64,    # ZC4-CLG
    9: 0.5,     # ZC5-CSM
    10: 0.003,  # ZS1-JEL
    11: 0.03,   # ZS2-CTG
    12: 0.26,   # ZS3-CHA
    13: 0.04,   # ZS4-LAR
    14: 1.1,    # PZ1-CIL
    15: 0.87,   # PZ2-DIN
    16: 0.46,   # PZ3-HNF
    17: 2.31, #* 2.5,   # PP1-DIA
    18: 1.5, #* 1.6,    # PP2-NAN
    19: 0.35, #* 1.5,   # PP3-PIC
    20: 0.17,   # BA1-BAC
    21: 4.3,    # DE2-DOC
    22: 10.1    # DE1-POC
}


# -------------------------------------------
# evaluation 2 - seasonal average B
# -------------------------------------------

START_DATA_PHYT_SEAS = '2015-01-01'
END_DATE_PHYT_SEAS = '2018-12-31'


# -------------------------------------------
# evaluation 3 - nutrients and annual pattern
# -------------------------------------------

P_FREE_INIT = 0.59# must be updated when Ecosim change is made!

# C_TO_N_RATIO = 106 / 16  # molar Redfield ratio
# NOT USING THESE UNITS BELOW! - Feb 2026
N_FREE_AVG_INIT = 18 # see ecosim_data_prep_2_nutrients.py, depth average 18 N, umol /L

N_BOUND_GROUPS = [1, 2, 3, 4, 5,
                  6, 7, 8, 9, 10,
                  11, 12, 13, 14, 15,
                  16, 17, 18, 19, 20, 22]
# N_BOUND_GROUPS = [17, 18, 19]

OBS_AVG_TYPE = "mean" # mean, median

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
# EWE_NUTR_LOADING_FILE = "Forcings Data Prep//Nutrient_loading_using_RDRS.csv"
EWE_NUTR_LOADING_FILE = "..//..///data//forcing//ECOSIM_in_3day_vars_1980-2018_fromASC_202506//ECOSIM_in_NEMO_varmixing_m_stdfilter_1980-2018.csv"
N_MULT_TYPE = "3day" # seasonal, monthly, 3day

SEASONAL_N_FLUX_MULT = {
    "Winter": 1.0,
    "Spring": 0.8,
    "Summer": 0.6,
    "Fall": 0.9
}


# nutrient loading function values
# Scenarios prio to 129
# MONTHLY_N_FLUX_MULT = {
#     1: 1.0, 2: 1.0, 3: 0.92, 4: 0.82, 5: 0.67, 6: 0.5,
#     7: 0.42, 8: 0.45, 9: 0.74, 10: 0.94, 11: 1.0, 12: 1.0
# }

# scenario 129
MONTHLY_N_FLUX_MULT = {
    1: 1.37, 2: 1.36, 3: 1.28, 4: 1.0, 5: 0.65, 6: 0.55,
    7: 0.53, 8: 0.47, 9: 0.75, 10: 0.95, 11: 1.4, 12: 1.55
}


# MONTHLY_N_FLUX_MULT = {
#     1: 1.0, 2: 1.0, 3: 1, 4: 1, 5: 1, 6: 1,
#     7: 1, 8: 1, 9: 1, 10: 1, 11: 1.0, 12: 1.0
# }

SEASON_MAP = {
    1:  "Winter", 2: "Winter", 12: "Winter",
    3:  "Spring", 4: "Spring", 5:  "Spring",
    6:  "Summer", 7: "Summer", 8:  "Summer",
    9:  "Fall",   10: "Fall",   11: "Winter",
}

ECOSIM_F_W_NUTRIENTS =  f"{OUTPUT_DIR_EVAL}//ecosim_{SCENARIO}_onerun_B_dates_seasons_nutrients.csv"
print(f"Nitrogen bound at initialisation: {N_B_INIT:.3f}")




# -------------------------------------------------------------------------
# Nutrient units harmonization (2026-02)
# -------------------------------------------------------------------------
# NOTE:
#   In the older nutrient workflow we mixed:
#     - biomass in areal carbon units (g C m-2)
#     - "free nutrient" as a concentration (umol/L)
#   That makes the implied budget inconsistent.
#
# Going forward, we evaluate nutrients as an *areal inventory* (g N m-2) over a
# specified surface layer thickness (e.g., 0–20 m), matching the Ecospace eval.
#
# Unit conversion reminders:
#   1 umol/L  == 1 mmol/m3
#   inventory (mmol/m2) = conc (mmol/m3) * thickness (m)
#   g N m-2   = mmol/m2 * 0.014   (since 1 mol N = 14 g)
# -------------------------------------------------------------------------

# Layer used for nutrient inventory (should match Ecospace NU_ZMIN_M/NU_ZMAX_M)
N_EVAL_ZMIN_M = 0.1
N_EVAL_ZMAX_M = 20.0
N_EVAL_LAYER_THICKNESS_M = float(N_EVAL_ZMAX_M - N_EVAL_ZMIN_M)

MMOL_TO_GN = 14.0 / 1000.0  # 0.014 g N per mmol N

def umolL_to_gNm2(conc_umol_L: float, thickness_m: float = N_EVAL_LAYER_THICKNESS_M) -> float:
    """Convert a concentration in umol/L to an areal inventory in g N m-2 over `thickness_m`."""
    # umol/L == mmol/m3
    return float(conc_umol_L) * float(thickness_m) * MMOL_TO_GN

# Convert initial *concentration* to initial *inventory*
N_FREE_INIT_GNM2 = 3.4

# C -> N multipliers (g N per g C). Keep aligned with Ecospace defaults.
N_C_TO_N_LIVING = 0.176
N_C_TO_N_DOM = 0.15
N_C_TO_N_POM = 0.07

# Group-specific overrides by numeric group id (Ecosim output columns are group numbers)
N_C_TO_N_BY_GROUP = {
    22: N_C_TO_N_POM,  # DE1-POC
    # 21: N_C_TO_N_DOM,  # DE2-DOC (if/when included)
}

def group_c_to_n(group_id: int) -> float:
    return float(N_C_TO_N_BY_GROUP.get(int(group_id), N_C_TO_N_LIVING))

# Initial N bound in biomass pools (g N m-2), based on Ecopath baseline biomasses
N_BOUND_INIT_GNM2 = sum(
    float(GROUPS_ECOPATH_B[g]) * group_c_to_n(g)
    for g in N_BOUND_GROUPS
    if g in GROUPS_ECOPATH_B
)

# Total inventory implied by the initialization (g N m-2)
N_TOTAL_INIT_GNM2 = float(N_FREE_INIT_GNM2 + N_BOUND_INIT_GNM2)

# How to anchor the initial budget when evaluating model time series
N_BOUND_INIT_MODE = "ecopath"          # {"ecopath", "t0_series"}
N_FREE_INIT_MODE = "config"            # {"config", "t0_preserve_total"}

# How to apply the upward-flux multiplier (if USE_N_MULT True)
#  - "scale_free":  N_free_used = N_free * mult
#  - "scale_total": N_free_used = (N_TOTAL_INIT_GNM2 * mult) - N_bound(t)
N_FLUX_APPLY_MODE = "scale_free"

# Plotting toggle for ecosim_eval_3_assess_nutrients.py
N_SHOW_PLOT = True
N_PLOT_INCLUDE_OBS = True   # set True if we also want an obs overlay in the Ecosim-only plot

print(f"[units] N_FREE_INIT_GNM2={N_FREE_INIT_GNM2:.3f} gN/m2  "
      f"N_BOUND_INIT_GNM2={N_BOUND_INIT_GNM2:.3f} gN/m2  "
      f"N_TOTAL_INIT_GNM2={N_TOTAL_INIT_GNM2:.3f} gN/m2")


# -------------------------------------------
# evaluation 4 - bloom timing
# -------------------------------------------

START_FULL_BLM = '1980-01-01' # 1996 is off because of imprecise spinup data prep (bloom is always day 20
END_FULL_BLM = '2018-12-31'

# Biomass columns for bloom detection (output from model is named by model group number)
BIOMASS_COLS_SATELLITE = ['17', '18', '19']
BIOMASS_COLS_C09 = ['17', '18', '19']

# adds a col before doing eval (sum across groups above, optionally)
TOTAL_BIOMASS_COL_SATELLITE = "Biomass_Total_Satellite"
TOTAL_BIOMASS_COL_C09 = "Biomass_Total_C09"

# Thresholds and bloom detection params
C09_USE_PCT_MAX = False # alt is to use the Sat methods see below
C09_PCT_MAX = 0.9
C09_LOG_TRNSFRM = False
C09_PCT_MAX_WINDOW_DAYS = 6
C09_USE_ANNUALORALL = "annual" # a future swithc; for threshold, use just that year's data to determine bloom?
C09_MEAN_OR_MEDIAN = "NOT SET UP" # a future switch (uses mean for now, or logmean)

SAT_LOG_TRNSFRM = False
SAT_THRESHOLD_FACTOR = 1.05
SAT_SUB_THRESHOLD_FACTOR = 0.7
SAT_MEAN_OR_MEDIAN = "median"
SAT_USE_ANNUALORALL = "annual" # for threshold, use just that year's data to determine bloom?
NUTRIENT_DRAWDOWN_FRAC = 0.6 # experimental - for nutrient definition of bloom (unused)


# -------------------------------------------
# evaluation 5 - zooplankton eval
# -------------------------------------------

Z_F_SEAS = "Zoopl_SofG_1996-2018_df_summary.csv" # this is output by long R script
Z_F_TOWLEV = "Zooplankton_B_C_gm2_EWEMODELGRP_Wide_NEMO3daymatch.csv" # this is output by short one
Z_P_PREPPED = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/4. Zooplankton/Zoop_Perryetal_2021/MODIFIED"

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
    'ZS4-LAR': 13,
}

# USER SETTING: choose one season and plot type ('bar' or 'line')
ZP_SEASON_CHOICE = 'Spring'  # e.g., 'Winter','Spring','Summer','Fall'
ZP_PLOT_TYPE = 'bar'     # 'bar' or 'line'
# USER SETTING: toggle friendly labels
ZP_USE_FRIENDLY_LABELS = True  # set False to use cryptic codes
ZP_FRIENDLY_MAP_ZC = {
    'ZF1-ICT': 'Ichthyo',
    'ZC1-EUP': 'Euphausiids',
    'ZC2-AMP': 'Amphipods',
    'ZC3-DEC': 'Decapods',
    'ZC4-CLG': 'Lg. Calanoid Copepods',
    'ZC5-CSM': 'Other Copepods',
    'ZS1-JEL': 'Jellyfish',
    'ZS2-CTH': 'Ctenophores',
    'ZS3-CHA': 'Chaetognaths',
    'ZS4-LAR': 'Larvaceans',
    'ZF1-ICH': 'Ichthyoplankton',
    'misc': 'Other',
    'Total': 'Total'
}
ZP_LOG_TRANSFORM = True
ZP_YEAR_START = 2000
ZP_YEAR_END   = 2018

ZP_FULLRN_START = 1980
ZP_FULLRN_END = 2018

ZP_SHOW_CNTS = True

# -----------------------------------------------------------------------------
# Tow filtering & anomaly stuff

#  context:
# - Ecosim anomalies are computed from the *paired tow table* (obs matched to model).
# - If we filter which tows are eligible (e.g., “deep/complete tows only”), we are
#   also changing the climatology used for z-scoring (mean/std) and the sample counts.
#
#  Other consideration
# - Sampling is temporally uneven (many more tows in recent years). If we pool all
#   tows to compute the climatology, then years with more tows (often recent years)
#   have more influence on the mean/std. That will tend to pull recent anomalies
#   closer to 0 by construction.
#
# These options allow us to:
# (A) Turn the tow filter ON/OFF so Ecosim can match Ecospace’s “include all matched tows”
# (B) Choose whether the climatology is tow-weighted vs year-weighted
# -----------------------------------------------------------------------------

# Tow eligibility filter applied BEFORE matching (only affects Ecosim workflow).
# "deep_or_complete": keep tows with (tow_prop >= min_prop) OR (start_depth >= min_start_depth_m)
# "none": keep all tows (closest to Ecospace workflow)
ZP_TOW_FILTER_MODE = "none"   # {"deep_or_complete", "none"}
ZP_TOW_MIN_PROP = 0.7
ZP_TOW_MIN_START_DEPTH_M = 150.0

# "tows": pool all tows in the season+year window (years with more tows weigh more)
# "yearly_mean": compute per-year means first; climatology is across years (each year weighs equally)
ZP_ANOM_CLIM_MODE = "tows"               # {"tows", "yearly_mean"}
