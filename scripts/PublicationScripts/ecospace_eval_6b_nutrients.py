"""
Analyzing Nutrients from Ecospace Outputs
by: G Oldford 2025
----------------------------------------------------------
Visualises nutrient pattern from Ecospace

Inputs:
- NetCDF of Ecospace outputs

- needs full revamp so instead of just doing nutrient climatology or masked nutrient pattern extract,
and to instead do some matching

Notes:
"""


import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import ecospace_eval_config as cfg

# ====== Configuration ======
SCENARIO           = cfg.ECOSPACE_SC
ECOSPACE_OUT_PATH  = cfg.NC_PATH_OUT
ECOSPACE_CODE      = cfg.ECOSPACE_SC
FILENM_STRT_YR     = cfg.ECOSPACE_RN_STR_YR
FILENM_END_YR      = cfg.ECOSPACE_RN_END_YR
ECOSPACE_OUT_F     = f"{ECOSPACE_CODE}_{FILENM_STRT_YR}-{FILENM_END_YR}.nc"
# DOMAIN_CONFIG_PATH = cfg.DOMAIN_P
# DOMAIN_FILE        = cfg.BT_DOMAIN_FILE
DOMAIN_MASK_PF     = cfg.BT_SAT_MASK_PF
OUTPUT_FIG_DIR     = cfg.FIGS_P

YEARS_TO_PLOT = list(range(cfg.NU_PLT_YR_ST, cfg.NU_PLT_YR_EN))

# C_TO_N_RATIO = cfg.C_TO_N_RATIO  # molar Redfield ratio # not used
N_B_INIT = cfg.N_B_INIT
N_TO_B_RATIO_INIT = N_B_INIT / cfg.TOT_B_INIT
N_FREE_INIT  = cfg.NU_FREE_AVG_INIT

INCLUDE_GRPS = cfg.NU_INCLUDE_GRPS
NU_PLT_GRPS = cfg.NU_PLT_GRPS

# INCLUDE_ONLY = ["PP1-DIA", "PP2-NAN", "PP3-PIC",
#                 "PZ1-CIL", "PZ2-DIN"]
# PROP_PP_MIXO = [1, 1, 1, 0.25, 0.25]  # Proportion of each group that is autotrophic

# ====== Load Data ======
def load_ecospace_dataset():
    return xr.open_dataset(os.path.join(ECOSPACE_OUT_PATH, ECOSPACE_OUT_F))

def load_mask():
    return xr.open_dataset(DOMAIN_MASK_PF)['mask']

# ====== Compute Nutrient Concentration ======
def compute_nutrient_concentration(ds, year, include_only, mask):

    # total_N = N_FREE_INIT + N_B_INIT
    # ratio_bound_per_C = N_B_INIT / INIT_TOT_PP_B

    year_ds = ds.sel(time=str(year))

    total_biomass = None
    biomass_by_group = {}
    for var in include_only:
        if var in year_ds:
            da = year_ds[var]
            if mask is not None:
                da = da.where(mask)
            mean_b = da.mean(dim=["row", "col"], skipna=True)
            biomass_by_group[var] = mean_b
            total_biomass = mean_b if total_biomass is None else total_biomass + mean_b

    N_bound = total_biomass * N_TO_B_RATIO_INIT
    N_free = N_B_INIT + N_FREE_INIT - N_bound
    N_free = N_free.where(N_free >= 0, 0)

    df = pd.DataFrame({
        "Date": pd.to_datetime(year_ds.time.values),
        "Day of Year": pd.to_datetime(year_ds.time.values).dayofyear,
        "Year": year,
        "Total Biomass (C)": total_biomass.values,
        "N Bound (g m^-2)": N_bound.values,
        "N Free (g m^-2)": N_free.values
    })

    for var, b in biomass_by_group.items():
        df[var] = b.values

    return df

# ====== Plot Nutrient Climatology ======
def plot_nutrient_concentration(df_nutrient_all, scenario):
    fig, ax1 = plt.subplots(figsize=(10, 5))

    grouped = df_nutrient_all.groupby("Day of Year")
    mean_nutrient = grouped["N Free (g m^-2)"].mean()
    q10 = grouped["N Free (g m^-2)"].quantile(0.1)
    q90 = grouped["N Free (g m^-2)"].quantile(0.9)

    ax1.fill_between(mean_nutrient.index, q10, q90, alpha=0.2, color='lightblue', label="10-90% range")
    ax1.plot(mean_nutrient.index, mean_nutrient, color='blue', label="Mean N Free")

    ax1.set_xlabel("Day of Year")
    ax1.set_ylabel("N Free (g N m$^{-2}$)", color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.set_title(f"Nutrient Drawdown Climatology – {scenario}")
    ax1.grid(True)

    ax2 = ax1.twinx()
    colors = plt.cm.tab10.colors
    for i, var in enumerate(NU_PLT_GRPS):
        if var in df_nutrient_all:
            grouped_b = df_nutrient_all.groupby("Day of Year")[var]
            mean_b = grouped_b.mean()
            q10_b = grouped_b.quantile(0.1)
            q90_b = grouped_b.quantile(0.9)
            ax2.fill_between(mean_b.index, q10_b, q90_b, alpha=0.1, color=colors[i % len(colors)])
            ax2.plot(mean_b.index, mean_b, label=var, color=colors[i % len(colors)])

    ax2.set_ylabel("Biomass (g C m$^{-2}$)", color='black')
    ax2.tick_params(axis='y', labelcolor='black')
    ax2.legend(loc='upper right')

    ax1.legend(loc='upper left')
    fig.tight_layout()
    plt.savefig(os.path.join(OUTPUT_FIG_DIR, f"nutrient_climatology_{scenario}.png"))
    plt.show()

# ====== Main Entry ======
def main():
    ds = load_ecospace_dataset()
    mask = load_mask()
    nutrient_dfs = []

    for yr in YEARS_TO_PLOT:
        print(f"Processing year: {yr}")
        df = compute_nutrient_concentration(ds, yr, INCLUDE_GRPS, mask)
        nutrient_dfs.append(df)

    df_nutrient_all = pd.concat(nutrient_dfs, ignore_index=True)

    # Save to CSV
    output_csv_path = os.path.join(OUTPUT_FIG_DIR, f"nutrient_concentration_{SCENARIO}.csv")
    df_nutrient_all.to_csv(output_csv_path, index=False)
    print(f"Nutrient concentration data saved to: {output_csv_path}")

    plot_nutrient_concentration(df_nutrient_all, scenario=SCENARIO)
    print("Done.")

if __name__ == "__main__":
    main()
