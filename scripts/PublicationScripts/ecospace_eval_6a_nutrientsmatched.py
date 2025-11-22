"""
Match 2-D Ecospace outputs to nutrient observations that already have ewe_row/ewe_col.


"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import ecospace_eval_config as cfg

# ---- Config (edit paths/names) ----

TIME_TOL_DAYS   = 7      # None → always nearest; else drop if farther than tol
DERIVE_N_FREE   = True   # set False to skip

NUTRIENTS_CSV = cfg.NU_F_PREPPED # from prep1_nutrients.py

SCENARIO           = cfg.ECOSPACE_SC
ECOSPACE_OUT_PATH  = cfg.NC_PATH_OUT
ECOSPACE_CODE      = cfg.ECOSPACE_SC
FILENM_STRT_YR     = cfg.ECOSPACE_RN_STR_YR
FILENM_END_YR      = cfg.ECOSPACE_RN_END_YR
ECOSPACE_NC      = rf"{ECOSPACE_OUT_PATH}/{ECOSPACE_CODE}_{FILENM_STRT_YR}-{FILENM_END_YR}.nc"

OUTPUT_DIR_EVAL = cfg.EVALOUT_P
OUTPUT_DIR_FIGS = cfg.FIGS_P

N_BOUND_GROUPS = cfg.NU_INCLUDE_GRPS
N_FREE_AVG_INIT = cfg.NU_FREE_AVG_INIT
N_B_INIT = cfg.N_B_INIT # inferred the absolute B in N units (see text and config file)

TOT_B_INIT = cfg.TOT_B_INIT
# USE_N_MULT = cfg.USE_N_MULT
# N_MULT_TYPE = cfg.N_MULT_TYPE

# SEASON_MAP = cfg.SEASON_MAP
# SEASONAL_N_FLUX_MULT = cfg.SEASONAL_N_FLUX_MULT
# MONTHLY_N_FLUX_MULT = cfg.MONTHLY_N_FLUX_MULT
# EWE_NUTR_LOADING_FILE = cfg.EWE_NUTR_LOADING_FILE
# OBS_AVG_TYPE = cfg.OBS_AVG_TYPE

OUT_MATCHED_CSV_1 = rf"{OUTPUT_DIR_EVAL}/ecospace_{SCENARIO}_nutrients_model_matched.csv"
OUT_MATCHED_CSV_2 = rf"{OUTPUT_DIR_EVAL}/ecospace_{SCENARIO}_nutrients_biweek_model_matched.csv"

# ---- Helpers ----
def nearest_time_join(obs_df: pd.DataFrame, model_times: pd.DatetimeIndex, tol_days: int | None):
    """Add the nearest model time to each obs row (within tolerance)."""
    # Build tiny frame of model times for asof-merge
    tdf = pd.DataFrame({"model_time": model_times}).sort_values("model_time")
    obs = obs_df.copy()
    obs["date"] = pd.to_datetime(obs["date"])  # already done in prep, but safe
    obs = obs.sort_values("date")
    merged = pd.merge_asof(
        obs, tdf, left_on="date", right_on="model_time", direction="nearest"
    )
    if tol_days is not None:
        too_far = (merged["model_time"] - merged["date"]).abs() > pd.Timedelta(days=tol_days)
        merged.loc[too_far, "model_time"] = pd.NaT
    return merged

def extract_model_at_points(ds: xr.Dataset, df: pd.DataFrame, vars_to_get: list[str]) -> pd.DataFrame:
    """Loop over rows (simple/clear) to pull model[var][time,row,col]."""
    out_rows = []
    # map model time to integer index for speed
    time_indexer = pd.Index(pd.to_datetime(ds["time"].values))
    i = 0
    for r in df.itertuples(index=False):
        print(r)
        if pd.isna(r.model_time):
            out = {f"{v}": np.nan for v in vars_to_get}
        else:
            t_idx = int(time_indexer.get_indexer([pd.to_datetime(r.model_time)], method="nearest")[0])
            rr, cc = int(r.ewe_row), int(r.ewe_col)
            out = {}
            for v in vars_to_get:
                if v in ds:
                    # guard on bounds
                    if 0 <= rr < ds.dims["row"] and 0 <= cc < ds.dims["col"]:
                        out[f"{v}"] = float(ds[v].isel(time=t_idx, row=rr, col=cc).values)
                    else:
                        out[f"{v}"] = np.nan
                else:
                    out[f"{v}"] = np.nan
        out_rows.append(out)
        #if i > 20: break
        i += 1
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(out_rows)], axis=1)

def derive_n_free(rowvals: pd.Series) -> float:
    # Sum biomass C across selected groups and convert using fixed ratio
    # total_c = np.nansum([rowvals.get(f"{v}", np.nan) for v in N_BOUND_GROUPS])
    total_c = rowvals["Total_Biomass_C"]
    total_n = total_c * (N_B_INIT / TOT_B_INIT)

    n_free  = N_B_INIT + N_FREE_AVG_INIT - total_n
    return max(0.0, n_free)

# ---- Main ----
def run():
    # 1) Observations (already prepped: month/season etc.)  :contentReference[oaicite:3]{index=3}
    obs = pd.read_csv(NUTRIENTS_CSV)
    obs["date"] = pd.to_datetime(obs["date"])

    # 2) Open model once
    ds = xr.open_dataset(ECOSPACE_NC)
    # Ensure datetime64 time
    if not np.issubdtype(ds["time"].dtype, np.datetime64):
        ds = ds.assign_coords(time=pd.to_datetime(ds["time"].values))
    model_times = pd.DatetimeIndex(pd.to_datetime(ds["time"].values))

    # 3) Nearest-time join (with tolerance)
    obs2 = nearest_time_join(obs, model_times, tol_days=TIME_TOL_DAYS)

    # 4) Pull model fields at [time,row,col]
    matched = extract_model_at_points(ds, obs2, N_BOUND_GROUPS)

    # Ensure types are clean
    matched["date"] = pd.to_datetime(matched["date"])
    matched["ewe_row"] = matched["ewe_row"].astype("Int64")
    matched["ewe_col"] = matched["ewe_col"].astype("Int64")

    matched["Total_Biomass_C"] = matched[N_BOUND_GROUPS].sum(axis=1)
    matched["model_N_free"] = matched.apply(derive_n_free, axis=1)

    # ---------- 2) Depth-integrate observations per (date,cell)
    # Floor depth==0 to 0.1 (defensive; prep likely did this already)
    matched.loc[matched["depth"] == 0, "depth"] = 0.1

    matched.to_csv(OUT_MATCHED_CSV_1, index=False)
    print(f"Saved matched table → {OUT_MATCHED_CSV_1}")


    # Choose statistic
    OBS_AVG_TYPE = "mean"  # or "median"

    if OBS_AVG_TYPE == "mean":
        obs_daily_cell = (
            matched.groupby(["date", "ewe_row", "ewe_col"], dropna=False)
            .agg(nitrogen_obs_int=("nitrogen", "mean"))
            .reset_index()
        )
    else:
        obs_daily_cell = (
            matched.groupby(["date", "ewe_row", "ewe_col"], dropna=False)
            .agg(nitrogen_obs_int=("nitrogen", "median"))
            .reset_index()
        )

    # ---------- 3) Make biweekly bins by (year,biweekly,cell) for obs
    obs_daily_cell["year"] = obs_daily_cell["date"].dt.year
    obs_daily_cell["doy"] = obs_daily_cell["date"].dt.dayofyear
    obs_daily_cell["biweekly"] = ((obs_daily_cell["doy"] - 1) // 14 + 1)

    obs_biweek_cell = (
        obs_daily_cell.groupby(["year", "biweekly", "ewe_row", "ewe_col"], dropna=False)
        .agg(avg_nitrogen_obs=("nitrogen_obs_int", "mean"))
        .reset_index()
    )

    # ---------- 3b) Build the model biweekly-by-cell table
    # If kept per-row model fields (e.g., PP1-DIA, ... and model_N_free) in `matched`:
    # First, collapse to (date,cell) in case an obs day had multiple depths mapping to same model time

    model_daily_cell = (
        matched.dropna(subset=["model_time"])
        .assign(year=lambda d: d["date"].dt.year,
                doy=lambda d: d["date"].dt.dayofyear,
                biweekly=lambda d: ((d["doy"] - 1) // 14 + 1))
        .groupby(["year", "biweekly", "ewe_row", "ewe_col"], dropna=False)
    )

    # pick what to compare:
    MODEL_COMPARE_COLS = ["model_N_free"]
    # MODEL_COMPARE_COLS = ["model_N_free"] + [v for v in matched.columns if
    #                                          v in {"PP1-DIA", "PP2-NAN", "PP3-PIC", "PZ1-CIL", "PZ2-DIN"}]

    model_biweek_cell = (
        model_daily_cell.agg({c: "mean" for c in MODEL_COMPARE_COLS})
        .reset_index()
    )

    # ---------- 4) Join obs↔model by (year,biweekly,row,col) and compute residuals
    paired = (
        obs_biweek_cell.merge(model_biweek_cell,
                              on=["year", "biweekly", "ewe_row", "ewe_col"],
                              how="left",
                              validate="m:1")
    )

    if "model_N_free" in paired.columns:
        paired["residual_N_free"] = paired["avg_nitrogen_obs"] - paired["model_N_free"]

    # Save biweekly-by-cell table ready for plotting/skill
    paired.to_csv(OUT_MATCHED_CSV_2, index=False)
    print(f"Saved biweekly matched table → {OUT_MATCHED_CSV_2}")

if __name__ == "__main__":
    run()
