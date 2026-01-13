"""ecospace_eval_6b_nutrients.py (REDO)

Goal
----
Infer an Ecospace "free dissolved N" time series (g N m-2) from model biomass pools
that are exported in carbon units (g C m-2), and overlay a biweekly climatology
against depth-integrated nutrient observations produced by ecospace_eval_6a_nutrientsmatched.py.

Key idea
--------
Ecospace tracks nutrients implicitly. We approximate a conserved, vertically-integrated
N budget over the evaluation layer (e.g., 0–20 m):

    N_free(t) = N_free_init + N_bound_init - N_bound(t)

where:
- N_free_init is the initial dissolved inventory (g N m-2) for the layer
  (e.g., from 18 umol/L climatological average -> user-provided 3.5 g N m-2).
- N_bound_init is initial N bound in model biomass pools included in the sum (g N m-2)
  (user-provided 1.56 g N m-2, OR computed from Ecopath B_init + C->N factors).
- N_bound(t) is the time-varying N bound in those biomass pools, obtained by converting
  each model group B (g C m-2) to g N m-2 using group-specific C->N multipliers.

Detritus
--------
If detrital pools (DOM/POM) are not exported by Ecospace, they cannot be included in
N_bound(t). For now this script defaults to using only non-detrital groups exported
in the matched table.

If detrital columns exist in your NetCDF / matched table later (e.g., DE1-POC / DE2-DOC),
add them to the conversion map and they will automatically be included.

"""

from __future__ import annotations

import os
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
import ecospace_eval_config as cfg


# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------

SCENARIO = cfg.ECOSPACE_SC
OUTPUT_DIR_EVAL = cfg.EVALOUT_P
OUTPUT_DIR_FIGS = cfg.FIGS_P

# Produced by ecospace_eval_6a_nutrientsmatched.py
PAIRED_BIWEEK_CSV = os.path.join(
    OUTPUT_DIR_EVAL,
    f"ecospace_{SCENARIO}_nutrients_biweek_model_matched.csv",
)

# Optional year filter
NU_PLT_YR_ST = getattr(cfg, "NU_PLT_YR_ST", None)
NU_PLT_YR_EN = getattr(cfg, "NU_PLT_YR_EN", None)  # exclusive

# Obs statistic within (year, biweekly) across cells
OBS_AVG_TYPE = getattr(cfg, "OBS_AVG_TYPE", getattr(cfg, "NU_OBS_AVG_TYPE", "mean"))  # mean|median
MODEL_AVG_TYPE = getattr(cfg, "MODEL_AVG_TYPE", "mean")  # mean|median

# ---- Conversion factors (g N per g C) ----
C_TO_N_LIVING = float(getattr(cfg, "NU_C_TO_N_LIVING", 0.176))
C_TO_N_DOM = float(getattr(cfg, "NU_C_TO_N_DOM", 0.15))
C_TO_N_POM = float(getattr(cfg, "NU_C_TO_N_POM", 0.07))

# Optional explicit overrides per group name, e.g.:
# NU_GROUP_C_TO_N = {"BA1-BAC": 0.20}
NU_GROUP_C_TO_N: Dict[str, float] = dict(getattr(cfg, "NU_GROUP_C_TO_N", {}))

# Candidate detrital column names if they exist
DOM_COL_CANDIDATES = list(getattr(cfg, "NU_DOM_COLS", ["DE2-DOC", "DE1-DOC", "DE1-DOM", "DE2-DOM"]))
POM_COL_CANDIDATES = list(getattr(cfg, "NU_POM_COLS", ["DE1-POC", "DE2-POC", "DE1-POM", "DE2-POM"]))

# ---- Initial inventory terms (g N m-2) ----
# Prefer explicit values if you add them to ecospace_eval_config.py.
# (If you only have a concentration in umol/L, convert to an areal inventory in gN/m2 first.)
NU_FREE_INIT_GNM2 = float(getattr(cfg, "NU_FREE_INIT_GNM2", 3.5))

# If provided, use directly. Otherwise compute from cfg.ES_GROUPS_ECOPATH_B.
NU_BOUND_INIT_GNM2 = getattr(cfg, "NU_BOUND_INIT_GNM2", None)

# Plot display toggle (interactive). If running headless/HPC, set False.
NU_SHOW_PLOT = bool(getattr(cfg, "NU_SHOW_PLOT", False))

# If True, force climatology index to include all biweekly bins 1..NU_BIWEEK_MAX
NU_FORCE_FULL_BIWEEK_AXIS = bool(getattr(cfg, "NU_FORCE_FULL_BIWEEK_AXIS", True))
NU_BIWEEK_MAX = int(getattr(cfg, "NU_BIWEEK_MAX", 26))

# Aggregation behavior:
# - True (default): keep only rows where BOTH obs and model exist for that (cell,time)
# - False: allow model-only rows into model climatology and obs-only rows into obs climatology
NU_REQUIRE_BOTH_FOR_AGG = bool(getattr(cfg, "NU_REQUIRE_BOTH_FOR_AGG", True))


# Which model groups to treat as biomass pools contributing to N_bound(t).
# If not set, we use the intersection of matched-table columns with:
#   cfg.ES_GROUPS_RELB (preferred) or cfg.NU_INCLUDE_GRPS (fallback)
NU_MODEL_GROUPS = list(getattr(cfg, "NU_MODEL_GROUPS", []))

# When True, include detrital columns if present in the matched table.
NU_INCLUDE_DETRITUS_IF_AVAILABLE = bool(getattr(cfg, "NU_INCLUDE_DETRITUS_IF_AVAILABLE", True))


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _pick_obs_gN_col(df: pd.DataFrame) -> str:
    """Return the best available obs column in g N m-2."""
    if "avg_nitrogen_gN_m2_obs" in df.columns:
        return "avg_nitrogen_gN_m2_obs"

    # fallback: compute from mmol/m^2 inventory if present
    if "avg_nitrogen_mmol_m2_obs" in df.columns:
        df["avg_nitrogen_gN_m2_obs"] = df["avg_nitrogen_mmol_m2_obs"] * (14.0067 / 1000.0)
        return "avg_nitrogen_gN_m2_obs"

    # older pipelines sometimes used a unitless name
    if "avg_nitrogen_obs" in df.columns:
        return "avg_nitrogen_obs"  # WARNING: unknown units

    raise ValueError(
        "Could not find an obs column to use. Expected 'avg_nitrogen_gN_m2_obs' (preferred) "
        "or 'avg_nitrogen_mmol_m2_obs'."
    )


def _infer_group_cols(df: pd.DataFrame) -> list[str]:
    """Choose model biomass columns to include in N_bound(t)."""

    if NU_MODEL_GROUPS:
        cols = [c for c in NU_MODEL_GROUPS if c in df.columns]
        if not cols:
            raise ValueError(
                f"NU_MODEL_GROUPS was provided but none were found in the matched table. "
                f"First few columns: {list(df.columns[:30])}"
            )
        return cols

    base_candidates: Iterable[str]
    if hasattr(cfg, "ES_GROUPS_RELB"):
        base_candidates = list(getattr(cfg, "ES_GROUPS_RELB"))
    else:
        base_candidates = list(getattr(cfg, "NU_INCLUDE_GRPS", []))

    cols = [c for c in base_candidates if c in df.columns]

    # Optionally include detritus (if present in matched table)
    if NU_INCLUDE_DETRITUS_IF_AVAILABLE:
        for c in DOM_COL_CANDIDATES + POM_COL_CANDIDATES:
            if c in df.columns and c not in cols:
                cols.append(c)

    if not cols:
        raise ValueError(
            "Could not infer any model biomass group columns from the matched table. "
            "Consider setting cfg.NU_MODEL_GROUPS explicitly."
        )

    return cols


def _build_c_to_n_map(group_cols: list[str]) -> Dict[str, float]:
    """Group-specific gN/gC conversion map."""
    c2n = {g: C_TO_N_LIVING for g in group_cols}

    for g in group_cols:
        if g in DOM_COL_CANDIDATES:
            c2n[g] = C_TO_N_DOM
        if g in POM_COL_CANDIDATES:
            c2n[g] = C_TO_N_POM

    # Explicit overrides win
    for g, v in NU_GROUP_C_TO_N.items():
        if g in c2n:
            c2n[g] = float(v)

    return c2n


def _compute_bound_N(df: pd.DataFrame, group_cols: list[str], c2n: Dict[str, float]) -> pd.Series:
    """Compute N_bound(t) from biomass pools (g C m-2 -> g N m-2)."""
    d = df[group_cols].apply(pd.to_numeric, errors="coerce")

    # Defensively clip negative biomass
    d = d.clip(lower=0)

    # Weighted sum with per-column conversion
    out = pd.Series(0.0, index=df.index, dtype=float)
    for g in group_cols:
        out = out + d[g].fillna(0.0) * float(c2n[g])

    # If all groups are NaN in a row, mark bound N as NaN
    all_nan = d.isna().all(axis=1)
    out.loc[all_nan] = np.nan

    return out


def _compute_bound_init_gN_m2(group_cols: list[str], c2n: Dict[str, float]) -> float:
    """Initial bound N for the same pools used in N_bound(t)."""
    if NU_BOUND_INIT_GNM2 is not None:
        return float(NU_BOUND_INIT_GNM2)

    # Compute from Ecopath B_init if available
    if hasattr(cfg, "ES_GROUPS_ECOPATH_B"):
        b0 = dict(getattr(cfg, "ES_GROUPS_ECOPATH_B"))
        s = 0.0
        missing = []
        for g in group_cols:
            if g in b0:
                s += float(b0[g]) * float(c2n[g])
            else:
                missing.append(g)

        # If detrital columns are in group_cols but not in ES_GROUPS_ECOPATH_B, we ignore them here
        # unless you provide NU_BOUND_INIT_GNM2 explicitly.
        if missing:
            print(
                "[6b] NOTE: Some group columns are not in cfg.ES_GROUPS_ECOPATH_B and were ignored "
                f"when computing bound_init: {missing}"
            )

        return float(s)

    raise ValueError(
        "NU_BOUND_INIT_GNM2 was not provided and cfg.ES_GROUPS_ECOPATH_B is not available. "
        "Add cfg.NU_BOUND_INIT_GNM2 (g N m-2) or provide ES_GROUPS_ECOPATH_B."
    )


def _load_paired_table(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"Paired biweekly table not found: {csv_path}\n"
            "Run ecospace_eval_6a_nutrientsmatched.py first."
        )

    df = pd.read_csv(csv_path)

    # Minimal required structure
    for c in ("year", "biweekly"):
        if c not in df.columns:
            raise ValueError(f"Paired table missing required column: {c!r}")

    df["year"] = df["year"].astype(int)
    df["biweekly"] = df["biweekly"].astype(int)

    if NU_PLT_YR_ST is not None:
        df = df[df["year"] >= int(NU_PLT_YR_ST)]
    if NU_PLT_YR_EN is not None:
        df = df[df["year"] < int(NU_PLT_YR_EN)]

    return df


def _aggregate_year_biweek(df: pd.DataFrame, obs_col: str, model_col: str) -> pd.DataFrame:
    """Collapse from (year, biweekly, cell) to (year, biweekly).

    If NU_REQUIRE_BOTH_FOR_AGG=True, rows missing either obs or model are dropped.
    If False, we keep obs and model aggregation *independently* by allowing NaNs.

    Output always contains columns:
      year, biweekly, avg_obs, avg_model, n_cells_obs, n_cells_model
    """

    if OBS_AVG_TYPE not in {"mean", "median"}:
        raise ValueError(f"OBS_AVG_TYPE must be 'mean' or 'median', got: {OBS_AVG_TYPE}")
    if MODEL_AVG_TYPE not in {"mean", "median"}:
        raise ValueError(f"MODEL_AVG_TYPE must be 'mean' or 'median', got: {MODEL_AVG_TYPE}")

    obs_func = "mean" if OBS_AVG_TYPE == "mean" else "median"
    mod_func = "mean" if MODEL_AVG_TYPE == "mean" else "median"

    d = df.copy()

    if NU_REQUIRE_BOTH_FOR_AGG:
        d = d.dropna(subset=[obs_col, model_col], how="any")

        out = (
            d.groupby(["year", "biweekly"], as_index=False)
            .agg(
                avg_obs=(obs_col, obs_func),
                avg_model=(model_col, mod_func),
                n_cells_obs=(obs_col, "count"),
                n_cells_model=(model_col, "count"),
            )
            .sort_values(["year", "biweekly"])
            .reset_index(drop=True)
        )
        return out

    # Independent aggregation (do not require both)
    g = d.groupby(["year", "biweekly"], as_index=False)
    out = (
        g.agg(
            avg_obs=(obs_col, obs_func),
            avg_model=(model_col, mod_func),
            n_cells_obs=(obs_col, "count"),
            n_cells_model=(model_col, "count"),
        )
        .sort_values(["year", "biweekly"])
        .reset_index(drop=True)
    )

    # NOTE: with independent aggregation, avg_obs can be NaN where n_cells_obs=0, etc.
    return out


def _compute_climatology(year_biweek: pd.DataFrame) -> pd.DataFrame:
    g = year_biweek.groupby("biweekly")

    clima = pd.DataFrame({"biweekly": sorted(year_biweek["biweekly"].unique())}).set_index("biweekly")

    clima["obs_center"] = g["avg_obs"].mean() if OBS_AVG_TYPE == "mean" else g["avg_obs"].median()
    clima["obs_q10"] = g["avg_obs"].quantile(0.1)
    clima["obs_q90"] = g["avg_obs"].quantile(0.9)

    clima["model_center"] = g["avg_model"].mean() if MODEL_AVG_TYPE == "mean" else g["avg_model"].median()
    clima["model_q10"] = g["avg_model"].quantile(0.1)
    clima["model_q90"] = g["avg_model"].quantile(0.9)

    # context
    clima["n_years_obs"] = g["avg_obs"].count()
    clima["n_years_model"] = g["avg_model"].count()

    # Optionally force 1..NU_BIWEEK_MAX so late-year bins show up as gaps instead of disappearing
    if NU_FORCE_FULL_BIWEEK_AXIS:
        full = pd.Index(range(1, NU_BIWEEK_MAX + 1), name="biweekly")
        clima = clima.reindex(full)

    return clima


def _plot_overlay(clima: pd.DataFrame, *, scenario: str) -> str:
    clima = clima.copy()
    clima.index = clima.index.astype(int)

    x = clima.index.to_numpy()

    plt.figure(figsize=(10, 5))

    # ---- Model ----
    m_fill = np.isfinite(clima["model_q10"].to_numpy()) & np.isfinite(clima["model_q90"].to_numpy())
    m_line = np.isfinite(clima["model_center"].to_numpy())
    if m_fill.any():
        plt.fill_between(
            x,
            clima["model_q10"].to_numpy(),
            clima["model_q90"].to_numpy(),
            where=m_fill,
            interpolate=True,
            alpha=0.2,
            label="Model 10–90%",
        )
    if m_line.any():
        plt.plot(x[m_line], clima["model_center"].to_numpy()[m_line], label="Model center")

    # ---- Obs ----
    o_fill = np.isfinite(clima["obs_q10"].to_numpy()) & np.isfinite(clima["obs_q90"].to_numpy())
    o_line = np.isfinite(clima["obs_center"].to_numpy())
    if o_fill.any():
        plt.fill_between(
            x,
            clima["obs_q10"].to_numpy(),
            clima["obs_q90"].to_numpy(),
            where=o_fill,
            interpolate=True,
            alpha=0.2,
            label="Obs 10–90%",
        )
    if o_line.any():
        plt.plot(x[o_line], clima["obs_center"].to_numpy()[o_line], label="Obs center")

    # Month ticks (as before)
    month_start_doy = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    month_ticks = [((d - 1) // 14 + 1) for d in month_start_doy]
    month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    ax = plt.gca()
    ax.set_xticks(month_ticks)
    ax.set_xticklabels(month_labels)

    # Keep the axis extending through the last biweekly bin
    ax.set_xlim(1, max(NU_BIWEEK_MAX, max(month_ticks)))
    ax.set_xlabel("Month")

    plt.ylabel("N (g N m$^{-2}$)")
    plt.title(f"Biweekly nutrient climatology (obs vs inferred model free N) – Ecospace {scenario}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    if NU_SHOW_PLOT:
        plt.show()

    _ensure_dir(OUTPUT_DIR_FIGS)
    out_png = os.path.join(OUTPUT_DIR_FIGS, f"nutrient_climatology_overlay_ecospace{scenario}_freeNredo.png")
    plt.savefig(out_png, dpi=200)


    plt.close()

    return out_png


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def run_nutrient_overlay() -> None:
    """Pipeline entrypoint (keeps the same name used by ecospace_evaluation_master_pipeline.py)."""

    df = _load_paired_table(PAIRED_BIWEEK_CSV)

    obs_col = _pick_obs_gN_col(df)
    group_cols = _infer_group_cols(df)
    c2n = _build_c_to_n_map(group_cols)

    bound_init = _compute_bound_init_gN_m2(group_cols, c2n)

    # Compute time-varying bound N and inferred free N
    df = df.copy()
    df["model_N_bound_gN_m2"] = _compute_bound_N(df, group_cols, c2n)
    df["model_N_free_gN_m2"] = NU_FREE_INIT_GNM2 + bound_init - df["model_N_bound_gN_m2"]

    # Guard against tiny negatives from rounding
    df.loc[df["model_N_free_gN_m2"] < 0, "model_N_free_gN_m2"] = 0.0

    print("[6b] Using obs column:", obs_col)
    print("[6b] Number of biomass pools contributing to bound N:", len(group_cols))
    print("[6b] C->N factors: living=", C_TO_N_LIVING, "DOM=", C_TO_N_DOM, "POM=", C_TO_N_POM)
    print("[6b] N_free_init (gN/m2)=", NU_FREE_INIT_GNM2)
    print("[6b] N_bound_init (gN/m2)=", bound_init)

    # Quick diagnostics for missing bins
    if "model_time" in df.columns:
        n_no_time = int(df["model_time"].isna().sum())
        print(f"[6b] Rows with no matched model_time: {n_no_time:,}")
    n_nan_bound = int(df["model_N_bound_gN_m2"].isna().sum())
    print(f"[6b] Rows with model_N_bound NaN (all group cols NaN): {n_nan_bound:,}")


    # Collapse from (year, biweekly, cell) to (year, biweekly)
    year_biweek = _aggregate_year_biweek(df, obs_col=obs_col, model_col="model_N_free_gN_m2")

    # Save year-biweekly series (debug)
    _ensure_dir(OUTPUT_DIR_EVAL)
    out_series_csv = os.path.join(OUTPUT_DIR_EVAL, f"ecospace_{SCENARIO}_nutrients_year_biweekly_freeNredo.csv")
    year_biweek.to_csv(out_series_csv, index=False)

    # Climatology
    clima = _compute_climatology(year_biweek)

    out_clima_csv = os.path.join(OUTPUT_DIR_EVAL, f"ecospace_{SCENARIO}_nutrients_climatology_biweekly_freeNredo.csv")
    clima.reset_index().to_csv(out_clima_csv, index=False)

    out_png = _plot_overlay(clima, scenario=SCENARIO)

    print(f"Loaded paired table: {PAIRED_BIWEEK_CSV}")
    print(f"Saved year-biweekly series → {out_series_csv}")
    print(f"Saved biweekly climatology → {out_clima_csv}")
    print(f"Saved overlay plot → {out_png}")


if __name__ == "__main__":
    run_nutrient_overlay()
