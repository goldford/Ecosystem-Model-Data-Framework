
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
ECOSPACE_SC = "SC202"
ECOSPACE_SC_FULL = "SC202"
ECOPATH_F_NM = "LTL_Carb_3day_ewe6_7_19295_v17_BigPC_ECOSPACEPARAMZ"
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
# 1 b - B before and after spin-up
# -------------------------------------------

ES_GROUPS_RELB = {
    "NK1-COH", "NK2-CHI", "NK3-FOR",
    "ZF1-ICT",
    "ZC1-EUP", "ZC2-AMP", "ZC3-DEC", "ZC4-CLG", "ZC5-CSM",
    "ZS1-JEL", "ZS2-CTH", "ZS3-CHA", "ZS4-LAR",
    "PZ1-CIL", "PZ2-DIN", "PZ3-HNF",
    "PP1-DIA", "PP2-NAN", "PP3-PIC",
    "BA1-BAC",
}

ES_GROUPS_ECOPATH_B = {
    "NK1-COH": 0.0046,
    "NK2-CHI": 0.00018,
    "NK3-FOR": 0.28,
    "ZF1-ICT": 0.03,
    "ZC1-EUP": 0.6,
    "ZC2-AMP": 0.52,
    "ZC3-DEC": 0.09,
    "ZC4-CLG": 0.64,
    "ZC5-CSM": 0.5,
    "ZS1-JEL": 0.003,
    "ZS2-CTH": 0.03,
    "ZS3-CHA": 0.26,
    "ZS4-LAR": 0.04,
    "PZ1-CIL": 1.1,
    "PZ2-DIN": 0.87,
    "PZ3-HNF": 0.46,
    "PP1-DIA": 2.31,
    "PP2-NAN": 1.50,
    "PP3-PIC": 0.35,
    "BA1-BAC": 0.17
}

ES_START_SPINUP = "1978-01-02"
ES_END_SPINUP =   "1979-12-31"



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

# use mask c09 instead of pnt?
BT_CREATE_MASKS = True  # Set to True to regenerate masks

# 20260102 - we may look at just applying this to that one year - or year by year would be ideal
# SGC2 is more accurate to map in suchy pub, but in 2005 clouds mask northern portion, affecting accuracy so SGC3 trims north
# results somewhat sensitive
BT_MASK_REGNM = "SGC4"
BT_USEC09_MASK_FOR_SAT = False # overrides above!
BT_USE_SAT_MASK_CO9 = False
BT_C09_ROW = 100 # overridden by above!
BT_C09_COL = 52

BT_CSV_SUCHY_PF = f"..//..//data//evaluation//ecospace_bloom_timing_SSoG_{ECOSPACE_SC}.csv"
BT_CSV_ALLEN_PF = f"..//..//data//evaluation//ecospace_bloom_timing_C09_{ECOSPACE_SC}.csv"

# if multiple vars listed here, it will sum across them when computing anomalies, bloom timing etc!
#OK with run 96 this is first time model fit has been okay with all pp groups
# VARIABLES_TO_ANALYZE_SAT = ["PP1-DIA"] #for 88_2 - pp1-dai only, annual, keep janfeb
BT_VARS_TO_ANALYZE_SAT = ["PP1-DIA", "PP2-NAN", "PP3-PIC"]
BT_ANNUAL_AVG_METHOD_SAT = "annual" # 'annual' or 'all' - should average bloom compared against be from all years, or just one year
BT_VARS_TO_ANALYZE_C09 = ["PP1-DIA"]
BT_ANNUAL_AVG_METHOD_C09 = "annual" # annual or all (this is always annual for C09 PCT MAX method)

# Bloom detection method
BT_LOG_TRANSFORM_SAT = False #
BT_MEAN_OR_MEDIAN_SAT = "median" # 'median" or "mean"
BT_THRESHOLD_FACTOR_SAT = 1.05
BT_SUB_THRESHOLD_FACTOR_SAT = 0.7

BT_EXCLUDE_DEC_JAN_SAT = False # working
BT_EXCLUDE_DEC_JAN_C09 = False

BT_LOG_TRANSFORM_C09 = True
BT_USE_PCT_MAX_C09 = False # new method added 2025-9-12 GO
BT_PCT_MAX_C09 = 0.9
BT_PCT_MAX_WINDOW_DAYS_C09 = 6



BT_DO_NUTRIENTS = False # another script does this now (#9?)
# OVERRIDE_REDFIELD = True # added by GO to help eval 2025-06-03


BT_MIN_Y_TICK = 38



# -------------------------------------------
# evaluation 5 - zooplankton eval
# -------------------------------------------

Z_F_SEAS = "Zoopl_SofG_1996-2018_df_summary.csv" # this is output by long R script
Z_F_TOWLEV = "Zooplankton_B_C_gm2_EWEMODELGRP_Wide_NEMO3daymatch.csv" # this is output by short one
Z_P_PREPPED = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/4. Zooplankton/Zoop_Perryetal_2021/MODIFIED"
Z_F_MATCH = f"Zooplankton_matched_to_model_out_{ECOSPACE_SC}.csv" # ecospace matched to zoop data

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
ZP_LOG_TRANSFORM = False
ZP_SHOW_CNTS = True             # annotate n_tows
ZP_ANOM_ALL_SEASONS = True     # True => loop Winter/Spring/Summer/Fall
ZP_MAKE_SCATTER = True          # make scatter figs too
ZP_SCATTER_LOG10 = False         # tow-level scatter uses log10 axes

ZP_YEAR_START = 2000
ZP_YEAR_END   = 2018

ZP_FULLRN_START = 1980
ZP_FULLRN_END = 2018

ZP_SHOW_CNTS = True

# -------------------------------------------------------------------------
# Anomaly climatology weighting (sampling-effort imbalance)
#
# Problem: tow sampling effort is uneven across years (often far more tows in
# recent years). If the climatology (mean/std) is computed by pooling all tows,
# heavily sampled years dominate the baseline and their anomalies tend to be
# closer to 0 by construction. If the climatology is computed from annual means
# with equal-year weighting, sparsely sampled years (e.g., 1–2 tows) can exert
# disproportionate influence because their annual mean is noisy.
#
# These options control how the climatological mean/std are estimated:
#   - "tow": use all tows pooled (effort-weighted; dominated by sampled years)
#   - "year_equal": use annual means; each year weight=1 (year-balanced)
#   - "year_weighted": use annual means with year weights ~ (n_tows**power),
#       optionally capped, and with a minimum tows-per-year threshold.
# -------------------------------------------------------------------------
ZP_ANOM_CLIM_MODE = "year_weighted"      # or "tow" / "year_equal"
ZP_ANOM_CLIM_WEIGHT_POWER = 0.5         # sqrt(n) compromise
ZP_ANOM_CLIM_WEIGHT_CAP = 25            # optional cap
ZP_ANOM_CLIM_MIN_TOWS_PER_YEAR = 3      # years with <3 tows won’t define baseline

# -------------------------------------------------------------------------
# Optional tow-eligibility filter (to mirror Ecosim workflow)
#   - "deep_or_complete": keep tows with (tow_prop >= min_prop) OR (start_depth >= min_start_depth_m)
#   - "none": keep all tows
# NOTE: if you change these, set RECOMPUTE_MATCH=True! or delete the old matched CSV.
# -------------------------------------------------------------------------
ZP_TOW_FILTER_MODE = "none"   # {"deep_or_complete", "none"}
ZP_TOW_MIN_PROP = 0.7
ZP_TOW_MIN_START_DEPTH_M = 150.0
ZP_TOW_MAX_BOTTOM_DEPTH_M = None   # optional additional cutoff



# -------------------------------------------
# evaluation 6 - nutrients (obs vs inferred model free N)
# -------------------------------------------

# --- Inputs ---
NU_F_PREPPED = f"{EVALOUT_P}//nutrients_ios_csop_combined_sampled.csv"  # must contain date,ewe_row,ewe_col,depth,nitrogen(umol/L)

# Years for the plotted climatology (exclusive end, like range())
NU_PLT_YR_ST = 1980
NU_PLT_YR_EN = 2019  # through 2018

# --- 6a: observation depth integration + binning ---
# Restrict obs years (END exclusive)
NU_OBS_YR_ST = 2015
NU_OBS_YR_EN = 2019  # through 2018

# Pooling behavior
NU_POOL_BY_CELL = True  # keep row/col in the biweekly bins (recommended)

# Depth integration window (m)
# NOTE: 6a warns that ZMIN=0.0 can cause all casts to fail coverage checks.
NU_ZMIN_M = 0.1
NU_ZMAX_M = 20.0

# Coverage logic
NU_REQUIRE_FULL_COVERAGE = True

# Endpoint padding tolerances (m): if a cast is within tolerance of endpoint, pad using nearest measured conc.
NU_SURFACE_PAD_TOL_M = 0.11
NU_DEEP_PAD_TOL_M = 0.0

# If beyond padding tolerance, allow extrapolation (constant endpoint value)
NU_DEPTH_EXTRAPOLATE = False

# Allow a cast with only one depth (constant profile) — usually False
NU_ALLOW_SINGLE_DEPTH = False

# Which cast metric becomes the pooled observational column:
#   - "integral" => mmol/m^2 inventory over [zmin,zmax] (then also write gN/m^2)
#   - "average"  => mmol/L depth-average over [zmin,zmax]
NU_OBS_VALUE_MODE = "integral"

# Pool casts within each (year,biweekly,row,col) bin
NU_OBS_BIN_AVG_TYPE = "mean"  # or "median"

# Optional cast-quality filters
NU_MIN_DEPTHS_PER_CAST = None
NU_MIN_MAXDEPTH_M = None

# Match tolerance (days) for nearest model time
NU_TIME_TOL_DAYS = 7


# --- 6b: plotting + aggregation behavior ---
# Interactive plotting toggle
NU_SHOW_PLOT = True

# Force climatology x-axis to include biweekly bins 1..26 so missing late-year bins show as gaps
NU_FORCE_FULL_BIWEEK_AXIS = True
NU_BIWEEK_MAX = 26

# If True: only aggregate rows where BOTH obs and model exist for that (cell,time)
# If False: aggregate obs and model independently (recommended when obs winter coverage is sparse)
NU_REQUIRE_BOTH_FOR_AGG = True

# Center statistic used for the climatology overlay
# (6b reads OBS_AVG_TYPE and MODEL_AVG_TYPE)
OBS_AVG_TYPE = "mean"    # "mean" or "median"
MODEL_AVG_TYPE = "mean"  # usually mean


# --- 6b: C -> N conversions to compute model bound N ---
# Default C->N multiplier for living pools (Redfield-ish; in mass units): gN = gC * 0.176
NU_C_TO_N_LIVING = 0.176

# Optional detritus-specific multipliers (mass-based)
NU_C_TO_N_DOM = 0.15   # if you later include DE2-DOC (DOM)
NU_C_TO_N_POM = 0.07   # if you later include DE1-POC (POM)

# Per-group override map (if/when you include detrital pools)
NU_C_TO_N_BY_GROUP = {
    "DE2-DOC": NU_C_TO_N_DOM,
    "DE1-POC": NU_C_TO_N_POM,
}


def _nu_c_to_n(group_code: str) -> float:
    """Return gN per gC for a given group."""
    return float(NU_C_TO_N_BY_GROUP.get(group_code, NU_C_TO_N_LIVING))


# --- Which model groups define N_bound(t)? ---
# Use an explicit list to avoid surprises from set ordering.
# This should match the *columns you expect to exist* in the 6a matched table.
# (Ecospace often does NOT export detrital groups; keep them excluded for now.)
NU_MODEL_GROUPS = [
    "NK1-COH", "NK2-CHI", "NK3-FOR",
    "ZF1-ICT",
    "ZC1-EUP", "ZC2-AMP", "ZC3-DEC", "ZC4-CLG", "ZC5-CSM",
    "ZS1-JEL", "ZS2-CTH", "ZS3-CHA", "ZS4-LAR",
    "PZ1-CIL", "PZ2-DIN", "PZ3-HNF",
    "PP1-DIA", "PP2-NAN", "PP3-PIC",
    "BA1-BAC",
    # "DE2-DOC",  # uncomment if/when Ecospace exports it
    # "DE1-POC",
]

# Tell 6a exactly which groups to extract (prevents it defaulting to sorted(ES_GROUPS_RELB))
NU_EXTRACT_GRPS = list(NU_MODEL_GROUPS)


# --- Initial total N inventory (g N m-2) for the depth window used in 6a ---
# Set this from your established conversion (e.g., your 18 umol/L -> 3.5 gN/m2 calculation).
# IMPORTANT: keep this consistent with [NU_ZMIN_M, NU_ZMAX_M].
NU_FREE_INIT_GNM2 = 3.5

# Choose either:
# (A) compute initial bound N from Ecopath initial biomasses of the selected model groups, OR
# (B) hard-code the value you computed elsewhere (e.g., 1.56 gN/m2 for non-detritals).

# (A) computed
NU_BOUND_INIT_GNM2 = sum(
    float(ES_GROUPS_ECOPATH_B[g]) * _nu_c_to_n(g)
    for g in NU_MODEL_GROUPS
    if g in ES_GROUPS_ECOPATH_B
)

# (B) override (uncomment to force)
# NU_BOUND_INIT_GNM2 = 1.56

# Convenience: implied initial fraction free (sanity check)
NU_P_FREE_INIT = NU_FREE_INIT_GNM2 / (NU_FREE_INIT_GNM2 + NU_BOUND_INIT_GNM2)

# Total initial areal N inventory for the inference (gN/m2)
NU_TOTAL_INIT_GNM2 = NU_FREE_INIT_GNM2 + NU_BOUND_INIT_GNM2



# -------------------------------------------
# evaluation 7 - taylor plots
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