"""
Analyzing Nutrients from Ecospace Outputs
by: G Oldford 2025
----------------------------------------------------------
Visualises nutrient pattern from Ecospace

Inputs:
- NetCDF of Ecospace outputs

Notes:
"""


import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

# ====== Configuration ======
SCENARIO = 'FULLKEY_SC117_8'
ECOSPACE_OUT_F = "Scv117_8-All_Groups_20250602_1978-2018.nc"
ECOSPACE_OUT_PATH = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//"
DOMAIN_MASK_PATH = "C://Users//Greig//Documents//github//Ecosystem-Model-Data-Framework//data//evaluation//suchy_ecospace_mask.nc"
OUTPUT_FIG_DIR = "..//..//figs//"

YEARS_TO_PLOT = list(range(1980, 2019))

C_TO_N_RATIO = 106 / 16  # molar Redfield ratio
N_FREE_INIT = 18.75
N_BOUND_INIT = 6.25
INIT_TOT_PP_B = 3.7  # g C m^-2, for initial phytoplankton biomass only

ALL_BIOMASS_VARS = [
    "NK1-COH", "NK2-CHI", "NK3-FOR",
    "ZF1-ICT", "ZC1-EUP", "ZC2-AMP",
    "ZC3-DEC", "ZC4-CLG", "ZC5-CSM",
    "ZS1-JEL", "ZS2-CTH", "ZS3-CHA",
    "ZS4-LAR", "PZ1-CIL", "PZ2-DIN",
    "PZ3-HNF", "PP1-DIA", "PP2-NAN",
    "PP3-PIC", "BAC-BA1"
]

INCLUDE_ONLY = ["PP1-DIA", "PP2-NAN", "PP3-PIC",
                "PZ1-CIL", "PZ2-DIN"]
PROP_PP_MIXO = [1, 1, 1, 0.25, 0.25]  # Proportion of each group that is autotrophic

# ====== Load Data ======
def load_ecospace_dataset():
    return xr.open_dataset(os.path.join(ECOSPACE_OUT_PATH, ECOSPACE_OUT_F))

def load_mask():
    return xr.open_dataset(DOMAIN_MASK_PATH)['mask']

# ====== Compute Nutrient Concentration ======
def compute_nutrient_concentration(ds, vars_list, year, include_only, mixo_props, mask):
    total_N = N_FREE_INIT + N_BOUND_INIT
    ratio_bound_per_C = N_BOUND_INIT / INIT_TOT_PP_B

    year_ds = ds.sel(time=str(year))

    total_biomass = None
    biomass_by_group = {}
    for var, prop in zip(include_only, mixo_props):
        if var in year_ds:
            da = year_ds[var] * prop
            if mask is not None:
                da = da.where(mask)
            mean_b = da.mean(dim=["row", "col"], skipna=True)
            biomass_by_group[var] = mean_b
            total_biomass = mean_b if total_biomass is None else total_biomass + mean_b

    N_bound = total_biomass * ratio_bound_per_C
    N_free = total_N - N_bound
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
    for i, var in enumerate(INCLUDE_ONLY):
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
        df = compute_nutrient_concentration(ds, ALL_BIOMASS_VARS, yr, INCLUDE_ONLY, PROP_PP_MIXO, mask)
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
