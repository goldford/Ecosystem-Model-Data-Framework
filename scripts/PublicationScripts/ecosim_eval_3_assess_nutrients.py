"""
Script: ecosim_eval_3_assess_nutrients.py
Purpose: Analyze nutrient concentrations from 1D Ecosim outputs
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import ecosim_eval_config as cfg
import matplotlib
matplotlib.use('TkAgg')

# ====== Config ======

SCENARIO = cfg.SCENARIO
ECOSIM_F_PREPPED_SINGLERUN = cfg.ECOSIM_F_PREPPED_SINGLERUN
OUTPUT_DIR_EVAL = cfg.OUTPUT_DIR_EVAL
OUTPUT_DIR_FIGS = cfg.OUTPUT_DIR_FIGS
N_DATA_REF = cfg.NUTRIENTS_F_PREPPED
N_BOUND_GROUPS = cfg.N_BOUND_GROUPS
N_FREE_AVG_INIT = cfg.N_FREE_AVG_INIT
N_B_INIT = cfg.N_B_INIT # inferred the absolute B in N units (see text and config file)
B_INIT_TOT = cfg.TOT_B_INIT
USE_N_MULT = cfg.USE_N_MULT
N_MULT_TYPE = cfg.N_MULT_TYPE
SEASON_MAP = cfg.SEASON_MAP
SEASONAL_N_FLUX_MULT = cfg.SEASONAL_N_FLUX_MULT
MONTHLY_N_FLUX_MULT = cfg.MONTHLY_N_FLUX_MULT
EWE_NUTR_LOADING_FILE = cfg.EWE_NUTR_LOADING_FILE
OBS_AVG_TYPE = cfg.OBS_AVG_TYPE

# exclude any 'spin-up' in Ecosim or things get thrown off
START_FULL = cfg.START_FULL
END_FULL = cfg.END_FULL




# New input assumptions

# INITIAL_NUTRIENT_CONC_FREE = cfg.INITIAL_NUTRIENT_CONC_FREE  # umol / L average annual
# MIN_N_CONC = cfg.MIN_N_CONC # observed min typical umol / L
# N_EQ_TO_B_EQ = cfg.N_EQ_TO_B_EQ
# PROP_N_FREE_INIT = cfg.PROP_N_FREE

# # Define a function to calculate % nutrient decline given % biomass increase
# def nutrient_decline(biomass_increase_fraction):
#     """
#     Calculate % decline in nutrient concentration
#     given a biomass increase (as a fraction).
#
#     biomass_increase_fraction: float
#         e.g. 0.2 for +20% increase
#     """
#     decline_fraction = biomass_increase_fraction / N_EQ_TO_B_EQ
#     return decline_fraction * 100

def derive_season_from_month(month: int) -> str:
    """Return season label for month using SEASON_MAP; default to 'Winter' if missing."""
    return SEASON_MAP.get(month, 'Winter')


def run_nutrient_eval():

    # ----------------------------------------------
    # Open data tables, fix dates
    # ----------------------------------------------

    # model output
    # filter to year!
    mdf = pd.read_csv(ECOSIM_F_PREPPED_SINGLERUN)
    mdf["date"] = pd.to_datetime(mdf["date"])
    mdf["year"] = mdf["date"].dt.year
    mdf["month"] = mdf["date"].dt.month
    mdf["day_of_year"] = mdf["date"].dt.dayofyear
    mdf["biweekly"] = ((mdf["day_of_year"] - 1) // 14 + 1)

    if 'season' not in mdf.columns:
        mdf['season'] = mdf['date'].dt.month.apply(derive_season_from_month)

    if N_MULT_TYPE == "monthly":
        mdf["multiplier"] = mdf["month"].map(MONTHLY_N_FLUX_MULT)
    elif N_MULT_TYPE == "seasonal":
        mdf["multiplier"] = mdf['season'].map(SEASONAL_N_FLUX_MULT).fillna(1.0)
    elif N_MULT_TYPE == "3day":
        forc = pd.read_csv(EWE_NUTR_LOADING_FILE)
        forc = forc.rename(columns={
            "dayofyear": "day_of_year",
            "stdfilter_varmixing_m": "multiplier"
        })
        mdf = mdf.merge(
            forc[["year", "day_of_year", "multiplier"]],
            on=["year", "day_of_year"],
            how="left"
        )
        mdf["multiplier"] = mdf["multiplier"].fillna(1.0)

    else:
        mdf["multiplier"] = 1.0



    # ----------------------------------------------
    # Compute 'bound' N
    # ----------------------------------------------

    N_BOUND_GROUPS_STR = []
    for group_col in N_BOUND_GROUPS:
        print(group_col)
        N_BOUND_GROUPS_STR.append(str(group_col))

    mdf["Total_Biomass_C"] = mdf[N_BOUND_GROUPS_STR].sum(axis=1)
    mdf["Total_Biomass_N"] = mdf["Total_Biomass_C"]  * N_B_INIT / B_INIT_TOT
    mdf["N_Free"] = N_B_INIT + N_FREE_AVG_INIT - mdf["Total_Biomass_N"]
    mdf["N_Free_Adjusted"] = mdf["N_Free"] * mdf["multiplier"] # in case of a seasonal or monthly ecosim 'nutrient loading'

    if USE_N_MULT:
        field_mdf = "N_Free_Adjusted"
    else:
        field_mdf = "N_Free"


    # ----------------------------------------------
    # Save file
    # ----------------------------------------------

    out_csv = cfg.ECOSIM_F_W_NUTRIENTS
    mdf.to_csv(out_csv, index=False)
    print(f"Nutrient analysis saved to {out_csv}")

    # exclude any ecosim 'spin-up':
    mdf = mdf[(mdf['date'] >= START_FULL) & (mdf['date'] <= END_FULL)]


    # ----------------------------------------------
    # observational data
    # ----------------------------------------------

    obs_df = pd.read_csv(N_DATA_REF)
    obs_df["date"] = pd.to_datetime(obs_df["date"])
    obs_df["year"] = obs_df["date"].dt.year
    obs_df["day_of_year"] = obs_df["date"].dt.dayofyear

    # Define biweekly bins
    obs_df["biweekly"] = ((obs_df["day_of_year"] - 1) // 14 + 1)
    # Adjust depth = 0 to 0.1 if not already done
    obs_df.loc[obs_df["depth"] == 0, "depth"] = 0.1
    # Depth integrate per year + biweekly
    if OBS_AVG_TYPE == "mean":
        depth_integrated = obs_df.groupby(["year", "biweekly"]).agg(
            avg_nitrogen=("nitrogen", "mean")
        ).reset_index()
    elif OBS_AVG_TYPE == "median":
        depth_integrated = obs_df.groupby(["year", "biweekly"]).agg(
            avg_nitrogen=("nitrogen", "median")
        ).reset_index()

    # Compute observed climatology (mean and quantiles) per biweekly across years
    biweek_grouped_obs = depth_integrated.groupby("biweekly")
    if OBS_AVG_TYPE == "mean":
        avg_obs = biweek_grouped_obs["avg_nitrogen"].mean()
    elif OBS_AVG_TYPE == "median":
        avg_obs = biweek_grouped_obs["avg_nitrogen"].median()

    q10_obs = biweek_grouped_obs["avg_nitrogen"].quantile(0.1)
    q90_obs = biweek_grouped_obs["avg_nitrogen"].quantile(0.9)


    # clima
    grouped_model = mdf.groupby("biweekly")
    mean_model = grouped_model[field_mdf].mean()
    q10_model = grouped_model[field_mdf].quantile(0.1)
    q90_model = grouped_model[field_mdf].quantile(0.9)

    # ----------------------------------------------
    # Plotting
    # ----------------------------------------------

    # ==== Model nutrient climatology (assuming you have these already computed) ====

    # ==== Plotting ====
    # plt.figure(figsize=(12, 6))

    # Plot model
    plt.fill_between(mean_model.index, q10_model, q90_model, alpha=0.2, color='lightblue', label="Model 10–90%")
    plt.plot(mean_model.index, mean_model, color='blue', label="Model Avg")

    # Plot observed
    plt.fill_between(avg_obs.index, q10_obs, q90_obs, alpha=0.2, color='lightgreen', label="Obs 10–90%")
    plt.plot(avg_obs.index, avg_obs, color='green', label="Obs Avg")

    month_start_doy = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    month_ticks = [((d - 1) // 14 + 1) for d in month_start_doy]
    month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    ax = plt.gca()
    ax.set_xticks(month_ticks)
    ax.set_xticklabels(month_labels)
    ax.set_xlim(1, mean_model.index.max())
    ax.set_xlabel("Month")

    plt.ylabel("Nitrogen (g N m$^{-2}$ or umol/L)")
    plt.title(f"Nutrient Biweekly Climatology – Ecosim {SCENARIO}")
    plt.legend()
    plt.grid()
    plt.tight_layout()

    plt.savefig(f"{OUTPUT_DIR_FIGS}//nutrient_climatology_overlay_ecosim{SCENARIO}.png")
    plt.show()
    plt.close()



if __name__ == "__main__":
    run_nutrient_eval()
