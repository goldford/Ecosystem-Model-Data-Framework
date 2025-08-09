"""
Script:   Match Ecosim data to zooplankton observations and generate seasonal
          and anomaly plots
Author:   G. Oldford (2025) – **refactored for modular plotting**

Purpose
-------
* Match Ecosim outputs to zooplankton observations (wide‑form).
* Compute seasonal means and annual anomalies.
* Produce two figure products:
  1. Seasonal biomass comparison (Obs vs Model)
  2. Annual anomaly time‑series (bar or line)
* Export comparison tables and skill metrics.

Public API
----------
The module still exposes **`run_zoop_eval()`** with **no arguments**, so
`ecosim_evaluation_master_pipeline.py` can keep importing & calling it
unchanged.
"""

from __future__ import annotations

# ── Imports ────────────────────────────────────────────────────────────────
import os
import logging
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import ecosim_eval_config as cfg

# Use same backend as original script to keep behaviour identical when run
# interactively but avoid GUI locks when executed headless in a pipeline.
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")

# ── Configuration pulled from ecosim_eval_config ──────────────────────────
SCENARIO           = cfg.SCENARIO
ECOSIM_CSV         = cfg.ECOSIM_F_PREPPED_SINGLERUN
OBS_CSV_SEAS       = os.path.join(cfg.Z_P_PREPPED, cfg.Z_F_SEAS)
OBS_CSV_TOWLEV     = os.path.join(cfg.Z_P_PREPPED, cfg.Z_F_TOWLEV)
OUTPUT_DIR_STATS   = cfg.OUTPUT_DIR_EVAL
OUTPUT_DIR_FIGS    = cfg.OUTPUT_DIR_FIGS

GROUP_MAP          = cfg.Z_GROUP_MAP                # {code: column‑id}
ZC_GROUPS          = [g for g in GROUP_MAP if g.startswith("ZC")]  # five Crustacean ‑type
ZS_GROUPS          = [g for g in GROUP_MAP if g.startswith(("ZS2","ZS3","ZS4"))] # four soft type but only 3 match obs
ZF_GROUPS          = [g for g in GROUP_MAP if g.startswith("ZF")] # one fish type (could be combined with ZS! (no ict in obs tho)

SEASONS            = ["Winter", "Spring", "Summer", "Fall"]
NCOLS              = 3
TIME_TOL           = pd.Timedelta(days=cfg.TIMESTEP_DAYS)

SEASON_CHOICE      = cfg.ZP_SEASON_CHOICE           # season for anomaly panels
PLOT_TYPE          = cfg.ZP_PLOT_TYPE               # 'bar' or 'line'
USE_FRIENDLY_LABELS= cfg.ZP_USE_FRIENDLY_LABELS
FRIENDLY_MAP       = cfg.ZP_FRIENDLY_MAP_ZC
ZP_YEAR_START      = cfg.ZP_YEAR_START
ZP_YEAR_END        = cfg.ZP_YEAR_END
ZP_FULLRN_START    = cfg.ZP_FULLRN_START
ZP_FULLRN_END      = cfg.ZP_FULLRN_END
ZP_LOG_TRANSFORM   = cfg.ZP_LOG_TRANSFORM
SHW_CNTS           = cfg.ZP_SHOW_CNTS

# Titles for panels
TITLE_MAP = {g: g for g in ZC_GROUPS}; TITLE_MAP["Total"] = "Total"

# ── Helper: skill metrics ─────────────────────────────────────────────────

def compute_stats(df: pd.DataFrame, obs_col: str, mod_col: str, *, log_or_anom: bool = False) -> dict:
    """Return classic skill metrics for two columns in *df* (NaN‑robust)."""
    o_ser, m_ser = df[obs_col], df[mod_col]
    mask         = o_ser.notna() & m_ser.notna()
    o, m         = o_ser[mask].values, m_ser[mask].values
    N            = len(o)
    if N == 0:
        return dict(N=0, MB=np.nan, MAE=np.nan, RMSE=np.nan, NRMSE=np.nan,
                    r=np.nan, R2=np.nan, MAPE=np.nan, NSE=np.nan, WSS=np.nan)
    mb   = np.nanmean(m - o)
    mae  = np.nanmean(np.abs(m - o))
    rmse = np.sqrt(np.nanmean((m - o) ** 2))
    nrmse= rmse if log_or_anom else rmse / np.nanmean(o) if np.nanmean(o) != 0 else np.nan
    r    = np.corrcoef(m, o)[0, 1] if N > 1 else np.nan
    r2   = r ** 2 if not np.isnan(r) else np.nan
    mape = np.nanmean(np.abs((m - o) / o)) * 100 if np.all(o != 0) else np.nan
    denom_nse = np.nansum((o - np.nanmean(o)) ** 2)
    nse  = 1 - np.nansum((m - o) ** 2) / denom_nse if denom_nse != 0 else np.nan
    denom_wss= np.nansum((np.abs(m - np.nanmean(o)) + np.abs(o - np.nanmean(o))) ** 2)
    wss  = 1 - np.nansum((m - o) ** 2) / denom_wss if denom_wss != 0 else np.nan
    return dict(N=N, MB=mb, MAE=mae, RMSE=rmse, NRMSE=nrmse,
                r=r, R2=r2, MAPE=mape, NSE=nse, WSS=wss)

# ── Plotting helpers ──────────────────────────────────────────────────────
def plot_seasonal_comparison(comp: pd.DataFrame, PLOT_GROUPS=ZC_GROUPS, label="ZC") -> None:
    """Barplots of seasonal mean biomass for each group (Obs vs Model)."""
    logger.info("Generating seasonal comparison figure …")
    fig, axes = plt.subplots(int(np.ceil(len(PLOT_GROUPS) / NCOLS)), NCOLS,
                             sharey=False,
                             figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS) / NCOLS)))
    axes = axes.flatten(); x = np.arange(len(SEASONS)); w = 0.35
    comp_idx = comp.set_index(["group", "season"])
    for i, (grp, ax) in enumerate(zip(PLOT_GROUPS, axes)):
        obs_vals = [comp_idx.loc[(grp, s), "obs_biomass"] if (grp, s) in comp_idx.index else np.nan for s in SEASONS]
        mod_vals = [comp_idx.loc[(grp, s), "model_biomass"] if (grp, s) in comp_idx.index else np.nan for s in SEASONS]
        ax.bar(x - w / 2, obs_vals, w, label="Obs", zorder=2); ax.bar(x + w / 2, mod_vals, w, label="Model", zorder=2)
        ax.grid(True, which='major', axis='both', linestyle='-', alpha=0.5, zorder=1)
        ax.set_xticks(x); ax.set_xticklabels(SEASONS); ax.set_title(TITLE_MAP.get(grp, grp))
        if i % NCOLS == 0: ax.set_ylabel("Biomass (g m⁻² DW)")
        if i == 0: ax.legend()
    # for ax in axes[len(PLOT_GROUPS):]: ax.axis("off")

    fig.tight_layout(); os.makedirs(OUTPUT_DIR_FIGS, exist_ok=True)
    plt.show()
    fname = Path(OUTPUT_DIR_FIGS) / f"ecosim_zoop_obs_vs_model_{SCENARIO}_{label}.png"; fig.savefig(fname, dpi=300); plt.close(fig)
    logger.info("Seasonal comparison figure saved to %s", fname)


# def plot_anomaly_panel(annual: pd.DataFrame, PLOT_GROUPS=ZC_GROUPS, label="ZC") -> None:
#     """Panel of anomaly time‑series for each group in *annual*."""
#     logger.info("Generating anomaly panel (%s, %s)…", SEASON_CHOICE, PLOT_TYPE)
#     years = sorted(annual["year"].unique()); x = np.arange(len(years)); w = 0.35
#     fig, axes = plt.subplots(int(np.ceil(len(PLOT_GROUPS) / NCOLS)), NCOLS,
#                              sharey=True,
#                              figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS) / NCOLS)))
#     axes = axes.flatten(); annual_idx = annual.set_index(["group", "year"])
#     for i, (grp, ax) in enumerate(zip(PLOT_GROUPS, axes)):
#         y_obs = [annual_idx.loc[(grp, yr), "obs_anom"] if (grp, yr) in annual_idx.index else np.nan for yr in years]
#         y_mod = [annual_idx.loc[(grp, yr), "model_anom"] if (grp, yr) in annual_idx.index else np.nan for yr in years]
#         if PLOT_TYPE == "bar":
#             ax.bar(x - w / 2, y_obs, w, label="Obs"); ax.bar(x + w / 2, y_mod, w, label="Model")
#         else:
#             ax.plot(years, y_obs, "-o", label="Obs"); ax.plot(years, y_mod, "-o", label="Model")
#         ax.set_xticks(x); ax.set_xticklabels(years, rotation=45)
#         ax.set_title(f"{FRIENDLY_MAP.get(grp, grp)} ({grp})" if USE_FRIENDLY_LABELS else grp)
#         if i % NCOLS == 0: ax.set_ylabel("Anomaly (z‑score)");
#         if i == 0: ax.legend()
#     for ax in axes[len(PLOT_GROUPS):]: ax.axis("off")
#     fig.tight_layout(); os.makedirs(OUTPUT_DIR_FIGS, exist_ok=True)
#
#     plt.show()
#     fname = Path(OUTPUT_DIR_FIGS) / f"zoop_anomalies_{SCENARIO}_{label}.png"; fig.savefig(fname, dpi=300); plt.close(fig)
#     logger.info("Anomaly panel saved to %s", fname)

def plot_anomaly_panel_paired(annual: pd.DataFrame,
                       PLOT_GROUPS=ZC_GROUPS,
                       label="ZC",
                       # NEW ↓
                       counts: pd.DataFrame | None = None,
                       season: str | None = None,
                       show_counts: bool = True,
                       count_yoffset: float = 0.05,   # as fraction of y-range
                       count_fontsize: int = 9,
                       count_color: str = "black") -> None:
    """Panel of anomaly time-series for each group in *annual*."""
    logger.info("Generating anomaly panel (%s, %s)…", SEASON_CHOICE, PLOT_TYPE)
    years = sorted(annual["year"].unique()); x = np.arange(len(years)); w = 0.35

    # NEW: build a dict of year -> n_tows for the chosen season
    counts_map = {}
    if show_counts and counts is not None:
        c = counts.copy()
        # Normalize column names if needed
        rename_map = {}
        if "Year" in c.columns: rename_map["Year"] = "year"
        if "Season" in c.columns: rename_map["Season"] = "season"
        if "n_tows" not in c.columns:
            # allow alternate names if you used something else
            for alt in ("n", "count", "n_rows", "N"):
                if alt in c.columns: rename_map[alt] = "n_tows"; break
        if rename_map: c = c.rename(columns=rename_map)

        if "season" in c.columns and season:
            c["season"] = c["season"].astype(str).str.capitalize()
            c = c[c["season"] == season.capitalize()]
        if "year" not in c.columns or "n_tows" not in c.columns:
            logger.warning("Counts table missing required columns; skipping annotations.")
        else:
            counts_map = dict(zip(c["year"], c["n_tows"]))

    fig, axes = plt.subplots(int(np.ceil(len(PLOT_GROUPS) / NCOLS)), NCOLS,
                             sharey=True,
                             figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS) / NCOLS)))
    axes = axes.flatten(); annual_idx = annual.set_index(["group", "year"])

    for i, (grp, ax) in enumerate(zip(PLOT_GROUPS, axes)):
        y_obs = [annual_idx.loc[(grp, yr), "obs_anom"] if (grp, yr) in annual_idx.index else np.nan for yr in years]
        y_mod = [annual_idx.loc[(grp, yr), "model_anom"] if (grp, yr) in annual_idx.index else np.nan for yr in years]

        ax.grid(True, which='major', axis='both', linestyle='-', alpha=0.5, zorder=1)

        if PLOT_TYPE == "bar":
            ax.bar(x - w / 2, y_obs, w, label="Obs", zorder=2)
            ax.bar(x + w / 2, y_mod, w, label="Model", zorder=2)
        else:
            ax.plot(years, y_obs, "-o", label="Obs", zorder=2)
            ax.plot(years, y_mod, "-o", label="Model", zorder=2)

        ax.set_xticks(x); ax.set_xticklabels(years, rotation=45)
        ax.set_title(f"{FRIENDLY_MAP.get(grp, grp)} ({grp})" if USE_FRIENDLY_LABELS else grp)
        if i % NCOLS == 0: ax.set_ylabel("Anomaly (z-score)")
        if i == 0: ax.legend()

        # NEW: add headroom for labels and annotate counts centered at each year
        if counts_map:
            ax.margins(y=0.15)  # give some vertical padding for text
            ymin, ymax = ax.get_ylim()
            yrange = ymax - ymin if ymax > ymin else 1.0
            for j, yr in enumerate(years):
                n = counts_map.get(yr)
                if n is None:
                    continue
                # place label above the taller of obs/model for that year
                vals = []
                if j < len(y_obs): vals.append(y_obs[j])
                if j < len(y_mod): vals.append(y_mod[j])
                # If both NaN, skip; else use nanmax
                if len(vals) == 0 or all(pd.isna(v) for v in vals):
                    continue
                y_top = np.nanmax(vals)
                # If both series are negative, this still places text just above the more negative bar/point
                y_text = y_top + yrange * count_yoffset
                ax.text(j, y_text, f"{int(n)}", ha="center", va="bottom",
                        fontsize=count_fontsize, color=count_color)

    for ax in axes[len(PLOT_GROUPS):]:
        ax.axis("off")

    fig.tight_layout(); os.makedirs(OUTPUT_DIR_FIGS, exist_ok=True)
    plt.show()
    fname = Path(OUTPUT_DIR_FIGS) / f"zoop_anomalies_ecosimvobs_{SCENARIO}_{label}.png"
    fig.savefig(fname, dpi=300); plt.close(fig)
    logger.info("Anomaly panel saved to %s", fname)

def plot_anomaly_panel_single(model_seas_anom: pd.DataFrame,
                       PLOT_GROUPS=ZC_GROUPS,
                       label="ZC",
                       season="Spring") -> None:

    """Panel of anomaly time-series for each group in model"""
    logger.info("Generating anomaly panel (%s, %s)…", season, PLOT_TYPE)
    years = sorted(model_seas_anom["year"].unique()); x = np.arange(len(years)); w = 0.8

    fig, axes = plt.subplots(int(np.ceil(len(PLOT_GROUPS) / NCOLS)), NCOLS,
                             sharey=True,
                             figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS) / NCOLS)))
    axes = axes.flatten(); annual_idx = model_seas_anom.set_index(["group", "year"])

    for i, (grp, ax) in enumerate(zip(PLOT_GROUPS, axes)):
        y_mod = [annual_idx.loc[(grp, yr), "model_anom"] if (grp, yr) in annual_idx.index else np.nan for yr in years]
        ax.grid(True, which='major', axis='both', linestyle='-', alpha=0.5, zorder=1)
        if PLOT_TYPE == "bar":
            ax.bar(x, y_mod, w, label="Model", zorder=2)
        else:
            ax.plot(years, y_mod, "-o", label="Model")

        ax.set_xticks(x)
        ax.set_xticklabels([str(yr) if yr % 5 == 0 else "" for yr in years], rotation=45)

        ax.set_title(f"{FRIENDLY_MAP.get(grp, grp)} ({grp})" if USE_FRIENDLY_LABELS else grp)
        if i % NCOLS == 0: ax.set_ylabel("Anomaly (z-score)")
        if i == 0: ax.legend()

    for ax in axes[len(PLOT_GROUPS):]:
        ax.axis("off")

    fig.tight_layout(); os.makedirs(OUTPUT_DIR_FIGS, exist_ok=True)
    plt.show()
    fname = Path(OUTPUT_DIR_FIGS) / f"zoop_anomalies_ecosim_fullrun_{SCENARIO}_{label}.png"
    fig.savefig(fname, dpi=300); plt.close(fig)
    logger.info("Anomaly panel saved to %s", fname)

def plot_seasonal_boxplots(dist: pd.DataFrame, PLOT_GROUPS=ZC_GROUPS, label="ZC") -> None:
    """
    Side-by-side Obs vs Model boxplots for each season and zoop group.
    *dist* must have columns: year, season, group, obs_biomass, model_biomass.
    """
    logger.info("Generating seasonal **box-plot** figure …")
    fig, axes = plt.subplots(int(np.ceil(len(PLOT_GROUPS) / NCOLS)), NCOLS,
                             figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS) / NCOLS)),
                             sharey=False)
    axes = axes.flatten()
    for i, (grp, ax) in enumerate(zip(PLOT_GROUPS, axes)):
        sub = dist[dist["group"] == grp]
        boxes, positions = [], []
        for j, seas in enumerate(SEASONS):
            s = sub[sub["season"] == seas]
            boxes.extend([s["obs_biomass"].values,
                          s["model_biomass"].values])
            positions.extend([j - 0.18, j + 0.18])   # little gap

        bp = ax.boxplot(boxes, positions=positions,
                        widths=0.30, showfliers=False,
                        patch_artist=True, zorder=2)
        # colour alternating boxes
        for k, b in enumerate(bp["boxes"]):
            b.set_facecolor("#4C72B0" if k % 2 == 0 else "#55A868")

        ax.grid(True, which='major', axis='both', linestyle='-', alpha=0.5, zorder=1)
        ax.set_xticks(range(len(SEASONS)))
        ax.set_xticklabels(SEASONS)
        ax.set_title(TITLE_MAP.get(grp, grp))
        if i % NCOLS == 0:
            ax.set_ylabel("Biomass (g C m⁻²)")
        if i == 0:
            ax.legend([Patch(color="#4C72B0"), Patch(color="#55A868")],
                      ["Obs", "Model"], loc="upper left")

    # for ax in axes[len(PLOT_GROUPS):]:
    #     ax.axis("off")
    fig.tight_layout()

    os.makedirs(OUTPUT_DIR_FIGS, exist_ok=True)
    fname = Path(OUTPUT_DIR_FIGS) / f"zoop_obs_vs_model_box_{SCENARIO}_{label}.png"
    fig.savefig(fname, dpi=300)
    plt.show()
    logger.info("Seasonal boxplot figure saved to %s", fname)



# ───────────────────────── I / O  ─────────────────────────────────────────
def load_observed_seasonal(csv: str, *, zp_groups: list[str],
                           plot_groups: list[str]) -> pd.DataFrame:
    """Return long-format seasonal obs table (cols: season, group, obs_biomass)."""
    obs = (
        pd.read_csv(csv)
          .assign(season=lambda d: d["season"].str.capitalize())
          .rename(columns={"modelgroup": "group", "adj_C_g_m2": "obs_biomass"})
    )

    obs_seas = (obs.groupby(["season", "group"], as_index=False)
                    ["obs_biomass"].sum())

    obs_wide = (obs_seas
                .pivot(index="season", columns="group", values="obs_biomass")
                .reset_index())
    obs_wide["Total"] = obs_wide[zp_groups].sum(axis=1)

    return obs_wide.melt(id_vars=["season"],
                         value_vars=plot_groups,
                         var_name="group",
                         value_name="obs_biomass")


def load_model_seasonal(csv: str, *, group_map: dict[str, int],
                        plot_groups: list[str],
                        year_range: tuple[int, int]) -> pd.DataFrame:
    """Return long-format seasonal means from Ecosim output (season, group, model_biomass)."""
    y0, y1 = year_range
    df = (pd.read_csv(csv, parse_dates=["date"])
            .assign(season=lambda d: d["season"].str.capitalize(),
                    year   =lambda d: d["date"].dt.year)
            .query(" @y0 <= year <= @y1 ")
            .rename(columns={str(v): k for k, v in group_map.items()}))

    df["Total"] = df[[g for g in group_map if g.startswith("ZC")]].sum(axis=1)

    long = pd.melt(df, id_vars=["season"],
                   value_vars=plot_groups,
                   var_name="group",
                   value_name="model_biomass")

    return (long.groupby(["season", "group"], as_index=False)
                ["model_biomass"].mean())


def merge_seasonal(obs_long: pd.DataFrame,
                   mod_seasonal: pd.DataFrame) -> pd.DataFrame:
    return obs_long.merge(mod_seasonal, on=["season", "group"], how="inner")


# ── I / O  (distribution-preserving) ──────────────────────────────────────
def load_observed_seasonal_dist(csv: str, *, zp_groups: list[str],
                                plot_groups: list[str]) -> pd.DataFrame:
    """
    Long-format table with one row per *year × season × group* and
    column 'obs_biomass'.
    """
    obs = (pd.read_csv(csv)
             .assign(season=lambda d: d["season"].str.capitalize())
             .rename(columns={"modelgroup": "group",
                              "adj_C_g_m2": "obs_biomass"}))

    # collapse species -> group, but *retain year*
    obs = (obs
           .groupby(["year", "season", "group"], as_index=False)
           ["obs_biomass"].sum())

    # add total
    obs_total = (obs[obs["group"].isin(zp_groups)]
                 .groupby(["year", "season"], as_index=False)["obs_biomass"]
                 .sum()
                 .assign(group="Total"))
    obs = pd.concat([obs, obs_total], ignore_index=True)

    return obs[obs["group"].isin(plot_groups)]


def load_model_seasonal_dist(csv: str, *, group_map: dict,
                             plot_groups: list[str],
                             year_range: tuple[int, int]) -> pd.DataFrame:
    """
    Long-format table with one row per *year × season × group* and
    column 'model_biomass'.
    """
    y0, y1 = year_range
    df = (pd.read_csv(csv, parse_dates=["date"])
            .assign(season=lambda d: d["season"].str.capitalize(),
                    year=lambda d: d["date"].dt.year)
            .query("@y0 <= year <= @y1")
            .rename(columns={str(v): k for k, v in group_map.items()}))

    df["Total"] = df[[g for g in group_map if g.startswith("ZC")]].sum(axis=1)

    long = pd.melt(df,
                   id_vars=["year", "season"],
                   value_vars=plot_groups,
                   var_name="group",
                   value_name="model_biomass")

    # average across monthly timesteps → one value per year×season×group
    return (long.groupby(["year", "season", "group"], as_index=False)
                 ["model_biomass"].mean())


# ── OBS loader from tow-level file ────────────────────────────────────────
def load_observed_seasonal_dist_from_tows(
    csv: str,
    *,
    plot_groups: list[str],
    zp_groups: list[str],            # groups that form 'Total'
    years: tuple[int, int] | None = None,
    region: str | None = None,
    agg: str | None = "yearly_mean"  # or None to use each tow directly
) -> pd.DataFrame:
    """
    Returns long-format obs with columns: year, season, group, obs_biomass.

    - If agg == 'yearly_mean': average across all tows within each (year, season, group)
      so each year contributes one value to the box.
    - If agg is None: keep tow-level values (one row per tow × season × group).
    """
    df = pd.read_csv(csv)

    # Standardize names present in your screenshot
    rename = {}
    # if "Year" in df.columns:   rename["Year"]   = "year"
    if "Season" in df.columns: rename["Season"] = "season"
    if "Date" in df.columns and "year" not in df.columns:
        df["year"] = pd.to_datetime(df["Date"]).dt.year
    df = df.rename(columns=rename)

    # Clean season labels to Title case (e.g., "Fall")
    if "season" not in df.columns or "year" not in df.columns:
        raise ValueError("Need 'season' and 'year' columns in the tow-level file.")
    df["season"] = df["season"].astype(str).str.strip().str.title()

    # Optional filters
    if years is not None:
        y0, y1 = years
        df = df.query("@y0 <= year <= @y1")
    if region is not None and "Region" in df.columns:
        df = df[df["Region"] == region]

    # Build Total (sum across the core zooplankton groups per tow)
    missing = [g for g in zp_groups if g not in df.columns]
    if missing:
        raise ValueError(f"These expected group columns are missing in obs file: {missing}")
    df["Total"] = df[zp_groups].sum(axis=1)

    # Keep only the groups we plan to plot
    keep_cols = ["year", "season"] + [g for g in plot_groups if g in df.columns]
    obs = df[keep_cols].copy()

    # Long form
    long = obs.melt(id_vars=["year", "season"],
                    var_name="group",
                    value_name="obs_biomass")

    # Aggregate choice
    if agg == "yearly_mean":
        # Each year contributes one number per season × group
        out = (long.groupby(["year", "season", "group"], as_index=False)
                    ["obs_biomass"].mean())
    elif agg is None:
        # Keep tow-level values (more rows; boxes reflect sampling intensity)
        out = long
    else:
        raise ValueError("agg must be 'yearly_mean' or None")

    return out


# ─────────────────────── tow-level prep  ──────────────────────────────────
def prepare_tow_obs(csv_tow: str) -> pd.DataFrame:
    """Return tow-level obs filtered to deep/complete tows, with date & season."""
    df = pd.read_csv(csv_tow)
    df["tow_depth_range"] = df["Tow_start_depth.m."].abs() - df["Tow_end_depth.m."].abs()
    df["tow_prop"] = df["tow_depth_range"] / df["Bottom_depth.m."]
    df = df[(df["tow_prop"] >= 0.7) | (df["Tow_start_depth.m."] >= 150)].copy()

    if "date" in df:
        df["date"] = pd.to_datetime(df["date"])
    else:
        df["date"] = pd.to_datetime(df[["Year", "Month", "Day"]])

    if "season" not in df:
        df["season"] = df["Month"].map({
            12: "Winter", 1: "Winter", 2: "Winter",
             3: "Spring", 4: "Spring", 5: "Spring",
             6: "Summer", 7: "Summer", 8: "Summer",
             9: "Fall", 10: "Fall", 11: "Fall"})
    df["season"] = df["season"].str.capitalize()
    return df


def prepare_model_for_match(csv: str, *, group_map: dict[str, int],
                            zp_groups: list[str]) -> pd.DataFrame:
    df = pd.read_csv(csv, parse_dates=["date"])
    df["season"] = df["season"].str.capitalize()
    for k, v in group_map.items():
        df = df.rename(columns={str(v): f"EWE-{k}"})
    keep = ["date", "season"] + [f"EWE-{g}" for g in zp_groups]
    return df[keep]


def match_obs_model(obs: pd.DataFrame, mod: pd.DataFrame,
                    *, tol: pd.Timedelta) -> pd.DataFrame:
    """Nearest-day merge; rows with NA model values are dropped."""
    obs = obs.sort_values("date")
    mod = mod.sort_values("date")
    matched = pd.merge_asof(obs, mod, on="date",
                            direction="nearest", tolerance=tol)
    # remove any unmatched rows (all ZC groups NA)
    drop_cols = [c for c in matched if c.startswith("EWE-")]

    matched = matched.dropna(subset=drop_cols)

    # --- unify the season column ----------------------------------------
    if "season_x" in matched.columns:
        matched["season"] = matched["season_x"].combine_first(matched["season_y"])
        matched = matched.drop(columns=[c for c in ["season_x", "season_y"] if c in matched])

    return matched


# ─────────────────────────  anomalies ────────────────────────────
def build_paired_long(matched: pd.DataFrame, *,
                      zp_groups: list[str]) -> pd.DataFrame:
    """Return long table with obs and model biomass for each group/Total."""
    recs: list[pd.DataFrame] = []
    for g in zp_groups + ["Total"]:
        if g == "Total":
            matched["Total"] = matched[zp_groups].sum(axis=1)
            obs_c, mod_c = "Total", None
        else:
            obs_c, mod_c = g, f"EWE-{g}"

        o = matched[["date", "season", "Index", obs_c]]\
                .rename(columns={obs_c: "obs_biomass"})
        if mod_c:
            m = matched[["date", "season", "Index", mod_c]]\
                    .rename(columns={mod_c: "model_biomass"})
        else:
            matched["model_biomass"] = matched[[f"EWE-{x}"
                                                for x in zp_groups]].sum(axis=1)
            m = matched[["date", "season", "Index", "model_biomass"]]

        recs.append(o.merge(m, on=["date", "season", "Index"])
                      .assign(group=g))
    paired = pd.concat(recs, ignore_index=True)
    paired["year"] = paired["date"].dt.year
    return paired


def build_model_long(mod_df: pd.DataFrame, *,
                      zp_groups: list[str]) -> pd.DataFrame:
    """Return long table with model biomass for each group/Total."""
    recs: list[pd.DataFrame] = []
    for g in zp_groups + ["Total"]:
        if g == "Total":
            mod_df["Total"] = mod_df[[f"EWE-{x}" for x in zp_groups]].sum(axis=1)
            # mod_df["Total"] = mod_df[zp_groups].sum(axis=1)
            mod_c = "Total"
        else:
            mod_c = f"EWE-{g}"

        if mod_c:
            m = mod_df[["date", "season", mod_c]]\
                    .rename(columns={mod_c: "model_biomass"})
        else:
            mod_df["model_biomass"] = mod_df[[f"EWE-{x}" for x in zp_groups]].sum(axis=1)
            m = mod_df[["date", "season", "model_biomass"]]

        recs.append(m.assign(group=g))
    mod_df_long = pd.concat(recs, ignore_index=True)
    mod_df_long["year"] = mod_df_long["date"].dt.year
    return mod_df_long


def compute_anomalies_paired(paired: pd.DataFrame, *, season: str,
                      year_range: tuple[int, int], log_transform: bool
                      ) -> pd.DataFrame:
    y0, y1 = year_range

    # keep only the chosen season *and* year window
    df = (paired
          .query("season == @season and @y0 <= year <= @y1")
          .copy())

    if log_transform:
        df["obs_biomass"]   = np.log(df["obs_biomass"]   + 1e-6)
        df["model_biomass"] = np.log(df["model_biomass"] + 1e-6)

    # climatology on the same subset
    clim = (df.groupby("group", as_index=False)
              .agg(mean_obs  = ("obs_biomass","mean"),
                   std_obs   = ("obs_biomass","std"),
                   mean_mod  = ("model_biomass","mean"),
                   std_mod   = ("model_biomass","std")))

    df = df.merge(clim, on="group", how="left")
    df["obs_anom"]   = (df["obs_biomass"]   - df["mean_obs"]) / df["std_obs"]
    df["model_anom"] = (df["model_biomass"] - df["mean_mod"]) / df["std_mod"]

    annual = (df.groupby(["year","group"], as_index=False)
                .agg(obs_anom=("obs_anom","mean"),
                     model_anom=("model_anom","mean")))
    return annual


def compute_anomalies_model(mod_df: pd.DataFrame, *, season: str,
                      year_range: tuple[int, int], log_transform: bool
                      ) -> pd.DataFrame:
    y0, y1 = year_range

    # keep only the chosen season *and* year window
    mod_df['year'] = mod_df["date"].dt.year

    df = (mod_df
          .query("season == @season and @y0 <= year <= @y1")
          .copy())

    if log_transform:
        df["model_biomass"] = np.log(df["model_biomass"] + 1e-6)

    # climatology on the same subset
    clim = (df.groupby("group", as_index=False)
              .agg( mean_mod  = ("model_biomass","mean"),
                   std_mod   = ("model_biomass","std")))

    df = df.merge(clim, on="group", how="left")
    df["model_anom"] = (df["model_biomass"] - df["mean_mod"]) / df["std_mod"]

    annual = (df.groupby(["year","group"], as_index=False)
                .agg(model_anom=("model_anom","mean")))
    return annual


# ───────────────────────── metrics helper  ────────────────────────────────
def skill_metrics_by_group(annual: pd.DataFrame, *,
                           groups: list[str], log_or_anom: bool) -> pd.DataFrame:
    rows = []
    for g in groups:
        d = annual[annual["group"] == g]
        if len(d) > 1:
            s = compute_stats(d, "obs_anom", "model_anom", log_or_anom=log_or_anom)
            s["group"] = g
            rows.append(s)
    return pd.DataFrame(rows)


def tow_counts_by_season_year(csv_path: str, out_csv: str | None = None) -> pd.DataFrame:
    """
    Count tow-level rows by (season, year) from OBS_CSV_TOWLEV.
    Tries to auto-detect column names like 'season'/'Season', 'year'/'Year', or derive from month.
    """
    df = pd.read_csv(csv_path)

    # --- Season column
    if "season" in df.columns:
        df["season"] = df["season"].astype(str).str.capitalize()
        season_col = "season"
    elif "Season" in df.columns:
        df["Season"] = df["Season"].astype(str).str.capitalize()
        season_col = "Season"
    else:
        # Derive season from month if needed
        if "Month" in df.columns:
            month_to_season = {
                12: "Winter", 1: "Winter", 2: "Winter",
                3: "Spring", 4: "Spring", 5: "Spring",
                6: "Summer", 7: "Summer", 8: "Summer",
                9: "Fall", 10: "Fall", 11: "Fall",
            }
            df["season"] = df["Month"].map(month_to_season)
            season_col = "season"
        else:
            raise ValueError("No season column and no Month column to derive from.")

    # --- Year column
    if "year" in df.columns:
        year_col = "year"
    elif "Year" in df.columns:
        year_col = "Year"
    else:
        # Try to parse from a date column if present
        date_col = "date" if "date" in df.columns else ("Date" if "Date" in df.columns else None)
        if not date_col:
            raise ValueError("No year/date column found.")
        df[date_col] = pd.to_datetime(df[date_col])
        df["year"] = df[date_col].dt.year
        year_col = "year"

    # --- Count rows
    counts = (
        df.groupby([year_col, season_col])
          .size()
          .rename("n_tows")
          .reset_index()
          .sort_values([year_col, season_col])
    )

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        counts.to_csv(out_csv, index=False)

    return counts


# ───────────────────────── orchestrator  ────────────────────────────────
def run_zoop_eval() -> None:
    """
    Match Ecosim outputs to zooplankton observations, plot seasonal & anomaly
    panels, and write comparison tables + skill metrics.  Keeps the same public
    contract (no arguments, side-effects only).
    """

    passes = [
        ("ZC",      ZC_GROUPS + ["Total"], ZC_GROUPS),
        ("ZSZF",  ZS_GROUPS + ["Total"], ZS_GROUPS)
    ]

    for label, plot_groups, groups in passes:
        PLOT_GROUPS = plot_groups # awkward but this includes 'total' for plotting
        GROUPS = groups

        # note the 'zc_groups' variable throughout the f(n)'s should be changed to just 'groups' or something

        # 1. Seasonal comparison ------------------------------------------------
        obs_long = load_observed_seasonal(
            OBS_CSV_SEAS,
            zp_groups=GROUPS,
            plot_groups=PLOT_GROUPS
        )

        mod_seasonal = load_model_seasonal(
            ECOSIM_CSV,
            group_map=GROUP_MAP,
            plot_groups=PLOT_GROUPS,
            year_range=(ZP_YEAR_START, ZP_YEAR_END)
        )

        comp = merge_seasonal(obs_long, mod_seasonal)
        comp.to_csv(Path(OUTPUT_DIR_STATS) /
                    f"ecosim_zoop_obs_vs_model_{SCENARIO}_{label}.csv", index=False)
        plot_seasonal_comparison(comp, PLOT_GROUPS, label)


        # --- OPTIONAL box-plot figure ----------------------------------------
        obs_dist = load_observed_seasonal_dist_from_tows(
            OBS_CSV_TOWLEV,
            plot_groups=PLOT_GROUPS,
            zp_groups=GROUPS,
            years=(ZP_YEAR_START, ZP_YEAR_END),  # or None
            agg="yearly_mean"  # or None for per-tow boxes
        )

        # 2) Load model distribution the way we set up earlier
        mod_dist = load_model_seasonal_dist(
            ECOSIM_CSV, group_map=GROUP_MAP,
            plot_groups=PLOT_GROUPS,
            year_range=(ZP_YEAR_START, ZP_YEAR_END)
        )

        # 3) Merge and plot
        dist_long = obs_dist.merge(mod_dist, on=["year", "season", "group"], how="inner")
        plot_seasonal_boxplots(dist_long, PLOT_GROUPS, label)  # <- the helper we added earlier


        # 2. Anomaly comparison -------------------------------------------------
        obs_tow   = prepare_tow_obs(OBS_CSV_TOWLEV)
        counts = tow_counts_by_season_year(OBS_CSV_TOWLEV, out_csv=f"{OUTPUT_DIR_STATS}/zoop_tow_counts_by_season_year.csv")
        print(counts)
        mod_prepped = prepare_model_for_match(ECOSIM_CSV,
                                            group_map=GROUP_MAP,
                                            zp_groups=GROUPS)
        matched   = match_obs_model(obs_tow, mod_prepped, tol=TIME_TOL)
        paired    = build_paired_long(matched, zp_groups=GROUPS)


        for SEAS in SEASONS:
            print(SEAS)
            matched_seas_anom    = compute_anomalies_paired(
                paired,
                season=SEAS,
                year_range=(ZP_YEAR_START, ZP_YEAR_END),
                log_transform=ZP_LOG_TRANSFORM
            )
            matched_seas_anom.to_csv(Path(OUTPUT_DIR_STATS) /
                          f"ecosim_zoop_anomalies_{SCENARIO}_{SEAS}_{label}.csv", index=False)

            plot_anomaly_panel_paired(matched_seas_anom, PLOT_GROUPS, label=label,
                               counts=counts, season=SEAS,
                               show_counts=SHW_CNTS, count_yoffset=0.06)

            # 3. Skill metrics ------------------------------------------------------
            model_seas_FULL    = build_model_long(mod_prepped, zp_groups=GROUPS)

            model_seas_anom_FULL = compute_anomalies_model(
                model_seas_FULL,
                season=SEAS,
                year_range=(ZP_FULLRN_START, ZP_FULLRN_END),
                log_transform=ZP_LOG_TRANSFORM
            )

            plot_anomaly_panel_single(model_seas_anom_FULL,
                                      PLOT_GROUPS, label=label,
                                      season=SEAS)



            # 3. Skill metrics ------------------------------------------------------
            stats_df = skill_metrics_by_group(
                matched_seas_anom,
                groups=PLOT_GROUPS,
                log_or_anom=ZP_LOG_TRANSFORM
            )
            print(stats_df)
            stats_df.to_csv(Path(OUTPUT_DIR_STATS) /
                            f"ecosim_zoop_performance_{SCENARIO}_{label}_{SEAS}.csv",
                            index=False)

# ── Run as script ──────────────────────────────────────────────────────
if __name__ == "__main__":
    run_zoop_eval()





