"""
Script: ecosim_eval_3_assess_nutrients.py
Purpose: Analyze nutrient concentrations from 1D Ecosim outputs
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import ecosim_eval_config
import matplotlib
matplotlib.use('TkAgg')

# ====== Config ======

SCENARIO = ecosim_eval_config.SCENARIO
ECOSIM_F_PREPPED_SINGLERUN = ecosim_eval_config.ECOSIM_F_PREPPED_SINGLERUN
OUTPUT_DIR_EVAL = ecosim_eval_config.OUTPUT_DIR_EVAL
OUTPUT_DIR_FIGS = ecosim_eval_config.OUTPUT_DIR_FIGS
N_DATA_REF = ecosim_eval_config.NUTRIENTS_F_PREPPED
N_BOUND_GROUPS = ecosim_eval_config.N_BOUND_GROUPS
N_FREE_AVG_INIT = ecosim_eval_config.N_FREE_AVG_INIT
N_B_INIT = ecosim_eval_config.N_B_INIT # inferred the absolute B in N units (see text and config file)
B_INIT_TOT = ecosim_eval_config.TOT_B_INIT
USE_N_MULT = ecosim_eval_config.USE_N_MULT
N_MULT_TYPE = ecosim_eval_config.N_MULT_TYPE
SEASON_MAP = ecosim_eval_config.SEASON_MAP
SEASONAL_N_FLUX_MULT = ecosim_eval_config.SEASONAL_N_FLUX_MULT
MONTHLY_N_FLUX_MULT = ecosim_eval_config.MONTHLY_N_FLUX_MULT

# exclude any 'spin-up' in Ecosim or things get thrown off
START_FULL = ecosim_eval_config.START_FULL
END_FULL = ecosim_eval_config.END_FULL


# New input assumptions

# INITIAL_NUTRIENT_CONC_FREE = ecosim_eval_config.INITIAL_NUTRIENT_CONC_FREE  # umol / L average annual
# MIN_N_CONC = ecosim_eval_config.MIN_N_CONC # observed min typical umol / L
# N_EQ_TO_B_EQ = ecosim_eval_config.N_EQ_TO_B_EQ
# PROP_N_FREE_INIT = ecosim_eval_config.PROP_N_FREE

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
    # filter to year! todo
    mdf = pd.read_csv(ECOSIM_F_PREPPED_SINGLERUN)
    mdf["date"] = pd.to_datetime(mdf["date"])

    # exclude any ecosim 'spin-up':
    mdf = mdf[(mdf['date'] >= START_FULL) & (mdf['date'] <= END_FULL)]

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
    else:
        mdf["multiplier"] = 1.0

        # observational data
    obs_df = pd.read_csv(N_DATA_REF)
    obs_df["date"] = pd.to_datetime(obs_df["date"])
    obs_df["year"] = obs_df["date"].dt.year
    obs_df["day_of_year"] = obs_df["date"].dt.dayofyear

    # Define biweekly bins
    obs_df["biweekly"] = ((obs_df["day_of_year"] - 1) // 14 + 1)
    # Adjust depth = 0 to 0.1 if not already done
    obs_df.loc[obs_df["depth"] == 0, "depth"] = 0.1
    # Depth integrate per year + biweekly
    depth_integrated = obs_df.groupby(["year", "biweekly"]).agg(
        mean_nitrogen=("nitrogen", "mean")
    ).reset_index()

    # Compute observed climatology (mean and quantiles) per biweekly across years
    biweek_grouped_obs = depth_integrated.groupby("biweekly")
    mean_obs = biweek_grouped_obs["mean_nitrogen"].mean()
    q10_obs = biweek_grouped_obs["mean_nitrogen"].quantile(0.1)
    q90_obs = biweek_grouped_obs["mean_nitrogen"].quantile(0.9)


    # ----------------------------------------------
    # Compute 'bound' N
    # ----------------------------------------------

    N_BOUND_GROUPS_STR = []
    for group_col in N_BOUND_GROUPS:
        N_BOUND_GROUPS_STR.append(str(N_BOUND_GROUPS[group_col-1]))

    mdf["Total_Biomass_C"] = mdf[N_BOUND_GROUPS_STR].sum(axis=1)
    mdf["Total_Biomass_N"] = mdf["Total_Biomass_C"]  * N_B_INIT / B_INIT_TOT
    mdf["N_Free"] = N_B_INIT + N_FREE_AVG_INIT - mdf["Total_Biomass_N"]
    mdf["N_Free_Adjusted"] = mdf["N_Free"] * mdf["multiplier"] # in case of a seasonal or monthly ecosim 'nutrient loading'

    if USE_N_MULT:
        field_mdf = "N_Free_Adjusted"
    else:
        field_mdf = "N_Free"

    # clima
    grouped_model = mdf.groupby("biweekly")
    mean_model = grouped_model[field_mdf].mean()
    q10_model = grouped_model[field_mdf].quantile(0.1)
    q90_model = grouped_model[field_mdf].quantile(0.9)

    # ----------------------------------------------
    # Save file
    # ----------------------------------------------

    out_csv = os.path.join(OUTPUT_DIR_EVAL, f"ecosim_{SCENARIO}_nutrients_revised.csv")
    mdf.to_csv(out_csv, index=False)
    print(f"Nutrient analysis saved to {out_csv}")

    # ----------------------------------------------
    # Plotting
    # ----------------------------------------------

    # ==== Model nutrient climatology (assuming you have these already computed) ====

    # ==== Plotting ====
    # plt.figure(figsize=(12, 6))

    # Plot model
    plt.fill_between(mean_model.index, q10_model, q90_model, alpha=0.2, color='lightblue', label="Model 10–90%")
    plt.plot(mean_model.index, mean_model, color='blue', label="Model Mean")

    # Plot observed
    plt.fill_between(mean_obs.index, q10_obs, q90_obs, alpha=0.2, color='lightgreen', label="Obs 10–90%")
    plt.plot(mean_obs.index, mean_obs, color='green', label="Obs Mean")

    plt.xlabel("Biweekly period")
    plt.ylabel("Nitrogen (g N m$^{-2}$ or umol/L)")
    plt.title(f"Nutrient Biweekly Climatology – {SCENARIO}")
    plt.legend()
    plt.grid()
    plt.tight_layout()

    plt.savefig(f"{OUTPUT_DIR_FIGS}//nutrient_climatology_overlay_{SCENARIO}.png")
    plt.show()
    plt.close()



if __name__ == "__main__":
    run_nutrient_eval()
