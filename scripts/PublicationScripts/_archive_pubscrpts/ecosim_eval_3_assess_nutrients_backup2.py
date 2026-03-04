"""
Script: ecosim_eval_3_assess_nutrients.py (UNITS FIXED)
Purpose: Infer a 1D Ecosim "free dissolved N" inventory time series (g N m-2)
         from Ecosim biomass outputs (g C m-2), and optionally plot a biweekly
         climatology.

Why this exists
---------------
Ecosim does not explicitly track dissolved nutrients. We approximate a conserved
surface-layer N inventory:

    N_free(t) = N_free_init + N_bound_init - N_bound(t)

where:
- N_free_init is the initial dissolved inventory (g N m-2) over the evaluation layer.
- N_bound_init is initial N bound in biomass pools (g N m-2), from Ecopath baseline.
- N_bound(t) is time-varying N bound in those pools, from Ecosim biomass (g C m-2)
  converted to g N m-2 using C->N multipliers.

Flux multiplier
---------------
If a seasonal/3-day "nutrient loading" multiplier is provided (e.g., to mimic
reduced upward flux under stratification), we can apply it in two ways:

- scale_free:  N_free_used = N_free * mult
- scale_total: N_free_used = (N_total_init * mult) - N_bound(t)

The "scale_total" option is usually what we want when the multiplier represents
a change in supply (total available inventory for the layer), not a direct scaling
of the dissolved pool after drawdown.

Inputs
------
- ECOSIM_F_PREPPED_SINGLERUN: must contain at least:
    date, and biomass columns as strings matching group ids (e.g. "17", "18", ...)

- NUTRIENTS_F_PREPPED (optional obs overlay):
    date, depth (m), nitrogen (umol/L) [and optionally ewe_row/ewe_col]

Outputs
-------
- ECOSIM_F_W_NUTRIENTS: CSV with:
    model_N_bound_gN_m2, model_N_free_gN_m2, model_N_free_used_gN_m2, multiplier, etc.

Created by: G Oldford
Unit harmonization update: 2026-02
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd

import ecosim_eval_config as cfg  # type: ignore


import matplotlib
if not bool(getattr(cfg, "N_SHOW_PLOT", False)):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------

SCENARIO = cfg.SCENARIO
ECOSIM_F_PREPPED_SINGLERUN = cfg.ECOSIM_F_PREPPED_SINGLERUN
OUTPUT_DIR_EVAL = cfg.OUTPUT_DIR_EVAL
OUTPUT_DIR_FIGS = cfg.OUTPUT_DIR_FIGS

# Optional obs overlay
N_DATA_REF = getattr(cfg, "NUTRIENTS_F_PREPPED", None)
OBS_AVG_TYPE = getattr(cfg, "OBS_AVG_TYPE", "mean")  # mean|median

# Groups that "bind" N (numeric IDs; Ecosim output columns are strings)
N_BOUND_GROUPS = list(getattr(cfg, "N_BOUND_GROUPS", []))

# Initial inventory terms (g N m-2). These are added in ecosim_eval_config_units_v2.py.
N_FREE_INIT_GNM2 = float(getattr(cfg, "N_FREE_INIT_GNM2", np.nan))
N_BOUND_INIT_GNM2 = float(getattr(cfg, "N_BOUND_INIT_GNM2", np.nan))
N_TOTAL_INIT_GNM2 = float(getattr(cfg, "N_TOTAL_INIT_GNM2", np.nan))

N_BOUND_INIT_MODE = str(getattr(cfg, "N_BOUND_INIT_MODE", "ecopath")).lower()
N_FREE_INIT_MODE = str(getattr(cfg, "N_FREE_INIT_MODE", "config")).lower()

# Layer definition for obs integration
ZMIN = float(getattr(cfg, "N_EVAL_ZMIN_M", 0.1))
ZMAX = float(getattr(cfg, "N_EVAL_ZMAX_M", 20.0))
MMOL_TO_GN = float(getattr(cfg, "MMOL_TO_GN", 14.0 / 1000.0))

# C->N multipliers
N_C_TO_N_LIVING = float(getattr(cfg, "N_C_TO_N_LIVING", 0.176))
N_C_TO_N_BY_GROUP: Dict[int, float] = dict(getattr(cfg, "N_C_TO_N_BY_GROUP", {}) or {})

# Flux multiplier controls
USE_N_MULT = bool(getattr(cfg, "USE_N_MULT", False))
N_MULT_TYPE = str(getattr(cfg, "N_MULT_TYPE", "none")).lower()
SEASON_MAP = getattr(cfg, "SEASON_MAP", {})
SEASONAL_N_FLUX_MULT = getattr(cfg, "SEASONAL_N_FLUX_MULT", {})
MONTHLY_N_FLUX_MULT = getattr(cfg, "MONTHLY_N_FLUX_MULT", {})
EWE_NUTR_LOADING_FILE = getattr(cfg, "EWE_NUTR_LOADING_FILE", None)
N_FLUX_APPLY_MODE = str(getattr(cfg, "N_FLUX_APPLY_MODE", "scale_total")).lower()

# Exclude any spin-up window
START_FULL = getattr(cfg, "START_FULL", None)
END_FULL = getattr(cfg, "END_FULL", None)

# Output
OUT_CSV = getattr(cfg, "ECOSIM_F_W_NUTRIENTS", os.path.join(OUTPUT_DIR_EVAL, f"ecosim_{SCENARIO}_nutrients.csv"))

# Plot toggles
N_SHOW_PLOT = bool(getattr(cfg, "N_SHOW_PLOT", False))
N_PLOT_INCLUDE_OBS = bool(getattr(cfg, "N_PLOT_INCLUDE_OBS", False))


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def derive_season_from_month(month: int) -> str:
    return SEASON_MAP.get(month, "Winter")


def group_c_to_n(group_id: int) -> float:
    """Return gN per gC for a numeric Ecosim group id."""
    return float(N_C_TO_N_BY_GROUP.get(int(group_id), N_C_TO_N_LIVING))


def _compute_multiplier(mdf: pd.DataFrame) -> pd.Series:
    if not USE_N_MULT:
        return pd.Series(1.0, index=mdf.index)

    if N_MULT_TYPE == "monthly":
        return mdf["month"].map(MONTHLY_N_FLUX_MULT).astype(float).fillna(1.0)

    if N_MULT_TYPE == "seasonal":
        return mdf["season"].map(SEASONAL_N_FLUX_MULT).astype(float).fillna(1.0)

    if N_MULT_TYPE == "3day":
        if not EWE_NUTR_LOADING_FILE:
            raise ValueError("N_MULT_TYPE='3day' but EWE_NUTR_LOADING_FILE is not set in config.")
        forc = pd.read_csv(EWE_NUTR_LOADING_FILE)
        forc = forc.rename(columns={
            "dayofyear": "day_of_year",
            "stdfilter_varmixing_m": "multiplier",
        })
        dfm = mdf.merge(
            forc[["year", "day_of_year", "multiplier"]],
            on=["year", "day_of_year"],
            how="left",
        )
        return pd.to_numeric(dfm["multiplier"], errors="coerce").fillna(1.0)

    # default (unknown type)
    return pd.Series(1.0, index=mdf.index)


def _t0_bound_from_series(df: pd.DataFrame) -> float:
    """Bound N at the first available timestep in the *series* (1D: earliest date)."""
    d = df.dropna(subset=["date", "model_N_bound_gN_m2"]).sort_values("date")
    if d.empty:
        raise ValueError("No non-NaN bound values to compute t0.")
    return float(d.iloc[0]["model_N_bound_gN_m2"])


def _umolL_profile_to_gNm2(depth_m: np.ndarray, conc_umol_L: np.ndarray, *, zmin: float, zmax: float) -> float:
    """
    Integrate a single vertical profile of concentration (umol/L == mmol/m3) over zmin..zmax
    and return g N m-2.
    """
    # Clean/sort
    ok = np.isfinite(depth_m) & np.isfinite(conc_umol_L)
    depth = depth_m[ok].astype(float)
    conc = conc_umol_L[ok].astype(float)

    if depth.size < 2:
        return np.nan

    # enforce positive depths
    depth = np.clip(depth, 0.0, None)
    order = np.argsort(depth)
    depth = depth[order]
    conc = conc[order]

    # Restrict to zmin..zmax with endpoint interpolation
    if depth.max() < zmin or depth.min() > zmax:
        return np.nan

    # Build an integration grid including endpoints
    z_nodes = depth[(depth >= zmin) & (depth <= zmax)]
    if z_nodes.size == 0:
        # no interior nodes; rely on endpoints
        z_nodes = np.array([], dtype=float)

    # Add endpoints explicitly (interpolate within available depth range)
    z_grid = np.unique(np.concatenate([[zmin], z_nodes, [zmax]]))
    # interpolate conc onto z_grid (no extrapolation beyond min/max)
    zmin_clip = max(zmin, depth.min())
    zmax_clip = min(zmax, depth.max())
    if zmin_clip > zmax_clip:
        return np.nan
    z_grid = z_grid[(z_grid >= zmin_clip) & (z_grid <= zmax_clip)]
    if z_grid.size < 2:
        return np.nan

    conc_grid = np.interp(z_grid, depth, conc)  # conc in mmol/m3

    mmol_m2 = float(np.trapz(conc_grid, z_grid))  # mmol/m2
    return mmol_m2 * MMOL_TO_GN  # g N m-2


def _compute_obs_climatology(obs_csv: str) -> pd.DataFrame:
    """
    Compute a biweekly climatology of observed dissolved N inventories (g N m-2)
    from a long-format obs table with columns:
        date, depth, nitrogen(umol/L), [optional ewe_row,ewe_col]
    """
    obs = pd.read_csv(obs_csv)
    obs["date"] = pd.to_datetime(obs["date"], errors="coerce")
    obs = obs.dropna(subset=["date"])

    # Standardize column names
    if "nitrogen" not in obs.columns:
        raise ValueError(f"Obs file missing 'nitrogen' column: {obs_csv}")

    if "depth" not in obs.columns:
        raise ValueError(f"Obs file missing 'depth' column: {obs_csv}")

    obs["depth"] = pd.to_numeric(obs["depth"], errors="coerce")
    obs["nitrogen"] = pd.to_numeric(obs["nitrogen"], errors="coerce")
    obs.loc[obs["depth"] == 0, "depth"] = 0.1

    obs["year"] = obs["date"].dt.year
    obs["day_of_year"] = obs["date"].dt.dayofyear
    obs["biweekly"] = ((obs["day_of_year"] - 1) // 14 + 1)

    # profile id: (date, cell) if present; else (date) only
    prof_keys = ["date"]
    if "ewe_row" in obs.columns and "ewe_col" in obs.columns:
        prof_keys += ["ewe_row", "ewe_col"]

    # integrate each profile
    prof = []
    for key, d in obs.groupby(prof_keys):
        gNm2 = _umolL_profile_to_gNm2(
            depth_m=d["depth"].to_numpy(),
            conc_umol_L=d["nitrogen"].to_numpy(),
            zmin=ZMIN,
            zmax=ZMAX,
        )
        if not np.isfinite(gNm2):
            continue
        row = dict(zip(prof_keys, key if isinstance(key, tuple) else (key,)))
        row.update({
            "year": int(d["year"].iloc[0]),
            "biweekly": int(d["biweekly"].iloc[0]),
            "obs_gN_m2": float(gNm2),
        })
        prof.append(row)

    if not prof:
        return pd.DataFrame()

    prof_df = pd.DataFrame(prof)

    # Reduce within (year, biweekly) across profiles
    reducer = np.mean if OBS_AVG_TYPE == "mean" else np.median
    yrbw = prof_df.groupby(["year", "biweekly"], as_index=False).agg(avg_obs=("obs_gN_m2", reducer))

    # Climatology across years
    g = yrbw.groupby("biweekly")["avg_obs"]
    clima = pd.DataFrame(index=sorted(yrbw["biweekly"].unique()))
    clima.index.name = "biweekly"
    clima["obs_center"] = g.mean() if OBS_AVG_TYPE == "mean" else g.median()
    clima["obs_q10"] = g.quantile(0.1)
    clima["obs_q90"] = g.quantile(0.9)
    return clima


def _compute_model_climatology(mdf: pd.DataFrame, value_col: str) -> pd.DataFrame:
    # Reduce within (year, biweekly) first so each year counts equally
    yrbw = mdf.groupby(["year", "biweekly"], as_index=False).agg(avg_model=(value_col, "mean"))
    g = yrbw.groupby("biweekly")["avg_model"]

    clima = pd.DataFrame(index=sorted(yrbw["biweekly"].unique()))
    clima.index.name = "biweekly"
    clima["model_center"] = g.mean()
    clima["model_q10"] = g.quantile(0.1)
    clima["model_q90"] = g.quantile(0.9)
    return clima


def _plot_overlay(obs_clima: pd.DataFrame, model_clima: pd.DataFrame, *, out_png: str) -> None:
    plt.figure(figsize=(10, 5))

    # Model
    x = model_clima.index.to_numpy()
    plt.fill_between(x, model_clima["model_q10"], model_clima["model_q90"], alpha=0.2, label="Model 10–90%")
    plt.plot(x, model_clima["model_center"], label="Model center")

    # Obs
    if obs_clima is not None and not obs_clima.empty:
        xo = obs_clima.index.to_numpy()
        plt.fill_between(xo, obs_clima["obs_q10"], obs_clima["obs_q90"], alpha=0.2, label="Obs 10–90%")
        plt.plot(xo, obs_clima["obs_center"], label="Obs center")

    month_start_doy = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    month_ticks = [((d - 1) // 14 + 1) for d in month_start_doy]
    month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    ax = plt.gca()
    ax.set_xticks(month_ticks)
    ax.set_xticklabels(month_labels)
    ax.set_xlim(1, max(26, max(month_ticks)))
    ax.set_xlabel("Month")

    plt.ylabel("Dissolved N inventory (g N m$^{-2}$)")
    plt.title(f"Biweekly nutrient climatology – Ecosim {SCENARIO}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(out_png, dpi=200)
    if N_SHOW_PLOT:
        plt.show()
    plt.close()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def run_nutrient_eval() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns
    -------
    mdf : pd.DataFrame
        Full model table with inferred nutrient series.
    model_clima : pd.DataFrame
        Biweekly climatology of model_N_free_used_gN_m2.
    """

    # ------------------
    # Load model output
    # ------------------
    mdf = pd.read_csv(ECOSIM_F_PREPPED_SINGLERUN)
    mdf["date"] = pd.to_datetime(mdf["date"], errors="coerce")
    mdf = mdf.dropna(subset=["date"]).copy()

    mdf["year"] = mdf["date"].dt.year
    mdf["month"] = mdf["date"].dt.month
    mdf["day_of_year"] = mdf["date"].dt.dayofyear
    mdf["biweekly"] = ((mdf["day_of_year"] - 1) // 14 + 1)

    if "season" not in mdf.columns:
        mdf["season"] = mdf["month"].apply(derive_season_from_month)

    mdf["multiplier"] = _compute_multiplier(mdf).astype(float)

    # ------------------
    # Compute bound N(t)
    # ------------------
    group_cols = [str(g) for g in N_BOUND_GROUPS if str(g) in mdf.columns]
    if not group_cols:
        raise ValueError("No N_BOUND_GROUPS columns found in model output. Check N_BOUND_GROUPS vs CSV columns.")

    # Convert each group biomass to N and sum
    bound_terms = []
    for g in N_BOUND_GROUPS:
        col = str(g)
        if col not in mdf.columns:
            continue
        mult = group_c_to_n(g)
        x = pd.to_numeric(mdf[col], errors="coerce").clip(lower=0)
        bound_terms.append(x * mult)

    mdf["model_N_bound_gN_m2"] = np.sum(bound_terms, axis=0)

    # ------------------
    # Choose initialization anchor
    # ------------------
    bound_init_ecopath = N_BOUND_INIT_GNM2
    bound_init_used = bound_init_ecopath

    if N_BOUND_INIT_MODE == "t0_series":
        bound_init_used = _t0_bound_from_series(mdf)

    if N_FREE_INIT_MODE == "t0_preserve_total":
        # preserve total inventory implied by config (free + ecopath bound)
        total_ref = float(N_FREE_INIT_GNM2 + bound_init_ecopath)
        bound_t0 = _t0_bound_from_series(mdf)
        free_init_used = max(0.0, total_ref - bound_t0)
    else:
        free_init_used = N_FREE_INIT_GNM2

    mdf["model_N_free_gN_m2"] = (free_init_used + bound_init_used) - mdf["model_N_bound_gN_m2"]
    mdf.loc[mdf["model_N_free_gN_m2"] < 0, "model_N_free_gN_m2"] = 0.0

    # ------------------
    # Apply multiplier (optional)
    # ------------------
    if USE_N_MULT:
        if N_FLUX_APPLY_MODE == "scale_free":
            mdf["model_N_free_used_gN_m2"] = mdf["model_N_free_gN_m2"] * mdf["multiplier"]
        elif N_FLUX_APPLY_MODE == "scale_total":
            mdf["model_N_free_used_gN_m2"] = (N_TOTAL_INIT_GNM2 * mdf["multiplier"]) - mdf["model_N_bound_gN_m2"]
        else:
            raise ValueError("N_FLUX_APPLY_MODE must be 'scale_free' or 'scale_total'")
        mdf.loc[mdf["model_N_free_used_gN_m2"] < 0, "model_N_free_used_gN_m2"] = 0.0
        value_col = "model_N_free_used_gN_m2"
    else:
        value_col = "model_N_free_gN_m2"

    # ------------------
    # Save full series
    # ------------------
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    mdf.to_csv(OUT_CSV, index=False)
    print(f"[ecosim nutrients] Saved → {OUT_CSV}")

    # Exclude spin-up period for climatology
    if START_FULL is not None:
        mdf = mdf[mdf["date"] >= pd.to_datetime(START_FULL)]
    if END_FULL is not None:
        mdf = mdf[mdf["date"] <= pd.to_datetime(END_FULL)]

    # ------------------
    # Climatology
    # ------------------
    model_clima = _compute_model_climatology(mdf, value_col=value_col)

    # ------------------
    # Plot (optional)
    # ------------------
    if N_SHOW_PLOT:
        obs_clima = None
        if N_PLOT_INCLUDE_OBS and N_DATA_REF:
            obs_clima = _compute_obs_climatology(N_DATA_REF)

        out_png = os.path.join(OUTPUT_DIR_FIGS, f"ecosim_nutrient_clima_overlay_{SCENARIO}_gNm2.png")
        os.makedirs(os.path.dirname(out_png), exist_ok=True)
        _plot_overlay(obs_clima, model_clima, out_png=out_png)
        print(f"[ecosim nutrients] Plot saved → {out_png}")

    return mdf, model_clima


if __name__ == "__main__":
    run_nutrient_eval()
