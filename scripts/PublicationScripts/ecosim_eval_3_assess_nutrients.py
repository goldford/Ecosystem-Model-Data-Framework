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
N_FREE_AVG = ecosim_config.N_FREE_AVG
N_BOUND_INIT = ecosim_config.N_FREE_DRAWDOWN
INIT_TOT_PP_B = ecosim_config.INIT_TOT_PP_B # g C m^-2, for initial phytoplankton biomass only
PP_COLUMNS = ecosim_config.PP_COLUMNS

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

    # Optional plot: N Free over time
    plt.figure(figsize=(10, 5))
    plt.plot(pd.to_datetime(df["date"]), df["N_Free_gm2"], label="N Free")
    plt.xlabel("Date")
    plt.ylabel("N Free (g N m^-2)")
    plt.title(f"Nutrient Time Series – {SCENARIO}")
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig(f"{OUTPUT_DIR_FIGS}//nutrient_timeseries_{SCENARIO}.png")
    plt.close()


if __name__ == "__main__":
    run_nutrient_eval()
