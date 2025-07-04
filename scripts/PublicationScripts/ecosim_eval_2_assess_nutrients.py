"""
Script: ecosim_eval_3_assess_nutrients.py
Purpose: Analyze nutrient concentrations from 1D Ecosim outputs
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import ecosim_config
import matplotlib
matplotlib.use('TkAgg')

# ====== Config ======
SCENARIO = ecosim_config.SCENARIO
OUTPUT_DIR_EVAL = ecosim_config.OUTPUT_DIR_EVAL
OUTPUT_DIR_FIGS = ecosim_config.OUTPUT_DIR_FIGS
ECOSIM_F_PREPPED_SINGLERUN = ecosim_config.ECOSIM_F_PREPPED_SINGLERUN

C_TO_N_RATIO = ecosim_config.C_TO_N_RATIO  # molar Redfield ratio
N_FREE_AVG = ecosim_config.N_FREE_AVG # ecopath base estimated N mmol / m-2
N_BOUND_INIT = ecosim_config.N_FREE_DRAWDOWN
INIT_TOT_PP_B = ecosim_config.INIT_TOT_PP_B # g C m^-2, for initial phytoplankton biomass only
N_BOUND_GROUPS = ecosim_config.N_BOUND_GROUPS
PP_COLUMNS = ecosim_config.PP_COLUMNS
N_DATA_REF = ecosim_config.NUTRIENTS_F_PREPPED # nutrients file to compare


def run_nutrient_eval():
    # Load the prepped Ecosim data
    df = pd.read_csv(ECOSIM_F_PREPPED_SINGLERUN)

    # Compute total phytoplankton biomass
    # Adjust this depending on your actual column names for PP groups in the CSV

    df["Total_Biomass_C"] = df[PP_COLUMNS].sum(axis=1)

    # Calculate N bound
    ratio_bound_per_C = N_BOUND_INIT / INIT_TOT_PP_B
    df["N_Bound_gm2"] = df["Total_Biomass_C"] * ratio_bound_per_C

    # Calculate N free
    total_N = N_FREE_AVG + N_BOUND_INIT
    df["N_Free_gm2"] = total_N - df["N_Bound_gm2"]
    df["N_Free_gm2"] = df["N_Free_gm2"].clip(lower=0)

    # Save results
    out_csv = f"{OUTPUT_DIR_EVAL}//ecosim_{SCENARIO}_nutrients.csv"
    df.to_csv(out_csv, index=False)
    print(f"Nutrient analysis saved to {out_csv}")

    # ------------------------------------------------------------
    # PLOTS
    # ------------------------------------------------------------

    # ==== Load observed data ====
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

    # model
    df["date"] = pd.to_datetime(df["date"])

    # Create 'day of year' field
    df["day_of_year"] = df["date"].dt.dayofyear

    # Create 'biweekly' field
    # This will number each biweekly period from 1 to ~26
    # Formula: day_of_year divided by 14 days per biweekly period, then ceil to get group number
    df["biweekly"] = ((df["day_of_year"] - 1) // 14 + 1).astype(int)

    # Verify by printing unique biweekly periods
    print(df["biweekly"].unique())

    # Group by day of year
    # grouped = df.groupby("day_of_year")
    # mean_nutrient = grouped["N_Free_gm2"].mean()
    # q10 = grouped["N_Free_gm2"].quantile(0.1)
    # q90 = grouped["N_Free_gm2"].quantile(0.9)

    df["biweekly"] = ((df["day_of_year"] - 1) // 14 + 1)
    grouped_model = df.groupby("biweekly")
    mean_model = grouped_model["N_Free_gm2"].mean()
    q10_model = grouped_model["N_Free_gm2"].quantile(0.1)
    q90_model = grouped_model["N_Free_gm2"].quantile(0.9)

    # Plot model
    plt.figure(figsize=(10, 5))
    # plt.fill_between(mean_nutrient.index, q10, q90, alpha=0.2, color='lightblue', label="10-90% range")
    # plt.plot(mean_nutrient.index, mean_nutrient, color='blue', label="Mean N Free")

    plt.fill_between(mean_model.index, q10_model, q90_model, alpha=0.2, color='lightblue', label="Model 10–90%")
    plt.plot(mean_model.index, mean_model, color='blue', label="Model Mean")

    # Plot observed
    plt.fill_between(mean_obs.index, q10_obs, q90_obs, alpha=0.2, color='lightgreen', label="Obs 10–90%")
    plt.plot(mean_obs.index, mean_obs, color='green', label="Obs Mean")

    plt.xlabel("Day of Year")
    plt.ylabel("N Free (g N m$^{-2}$)")
    plt.title(f"Nutrient Drawdown – {SCENARIO}")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR_FIGS}//nutrients_{SCENARIO}.png")
    plt.show()
    plt.close()


if __name__ == "__main__":
    run_nutrient_eval()
