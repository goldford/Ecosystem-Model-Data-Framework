"""
Ecospace Zooplankton Evaluation — Combined 5b/5c/5d (aligned to 4a pattern)
by: G. Oldford

Process
--------------
1) Match zooplankton observations to Ecospace outputs (former 5b)
2) Visualize & compute stats (former 5c)
3) Make anomaly comparisons incl. NPGO overlays (former 5d)

----------
- Configure paths and options in `ecospace_eval_config.py` (mirrors 4a style)
- Then either:
    from ecospace_eval_5_combined_zoop import run_zoop_eval
    run_zoop_eval(recompute_match=True)

Outputs
-------
- Matched CSV:  Zooplankton_matched_to_model_out_{ECOSPACE_CODE}.csv
- Skill CSV:    model_skill_stats_log10_{ECOSPACE_CODE}.csv
- Anom figs:    {ECOSPACE_CODE}_anom_panel_Zoop_vs_Model_with_NPGO.png, etc.

Notes
-----
- Keeps logic and naming close to  4a script for consistency
- Uses cfg.* for paths, run switches, years, scenario codes, etc.
"""


from __future__ import annotations
import os
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from dataclasses import field
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import seaborn as sns

import netCDF4 as nc
from datetime import datetime

# Local helpers & config (same approach as 4a)
from helpers import (
    find_nearest_point_from1D,
    find_closest_date,
)
import ecospace_eval_config as cfg

# ============================================================
# Configuration adaptor (pulls needed bits from cfg like 4a)
# ============================================================
@dataclass
class Eval5Config:
    # INPUTS
    ecospace_code: str = cfg.ECOSPACE_SC
    ecospace_nc_name: str = f"{cfg.ECOSPACE_SC}_{cfg.ECOSPACE_RN_STR_YR}-{cfg.ECOSPACE_RN_END_YR}.nc"

    # PATHS
    NC_PATH_OUT: str = cfg.NC_PATH_OUT
    EVALOUT_P: str = cfg.EVALOUT_P
    FIGSOUT_P: str = cfg.FIGS_P
    BASEMAP_P: str = cfg.ECOSPACE_MAP_P

    # Zoop & mapping input files (set these in cfg or keep defaults here)
    ZOOP_P: str = getattr(cfg, "Z_P_PREPPED", "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/4. Zooplankton/Zoop_Perryetal_2021/MODIFIED")
    ZOOP_CSV: str = getattr(cfg, "Z_F_TOWLEV", "Zooplankton_B_C_gm2_EWEMODELGRP_Wide.csv")
    GRID_MAP_CSV: str = getattr(cfg, "ECOSPACE_GRID_RC_CSV", "Ecospace_grid_20210208_rowscols.csv")
    ZOOP_CSV_MATCH: str = getattr(cfg, "Z_F_MATCH", f"Zooplankton_matched_to_model_out_{ecospace_code}.csv")

    # Optional external series
    NPGO_CSV: str = getattr(cfg, "NPGO_CSV", "npgo.csv")

    # YEARS
    START_YEAR: int = getattr(cfg, "ZP_START_YEAR", 1980)
    END_YEAR: int = getattr(cfg, "ZP_END_YEAR", 2018)

    # GROUPS (consistent w/ 5b/5c/5d)
    ZOOP_GROUPS: List[str] = field(default_factory=lambda: ["ZC1-EUP", "ZC2-AMP", "ZC3-DEC", "ZC4-CLG", "ZC5-CSM",
            "ZS2-CTH", "ZS3-CHA", "ZS4-LAR", "ZF1-ICT"])


    # RUN SWITCHES
    RECOMPUTE_MATCH: bool = getattr(cfg, "ZP_RECOMPUTE_MATCH", True)
    MAKE_VIZ: bool = getattr(cfg, "ZP_MAKE_VIZ", True)
    MAKE_ANOM: bool = getattr(cfg, "ZP_MAKE_ANOM", True)

   # climatology weighting for anomalies
    # - ------------------------
    # Options:
    #   "tow"        : current behavior (clim from all tows pooled)
    #   "year_equal" : clim from annual means, each year weight=1
    #   "year_weighted" : clim from annual means, year weight ~ (n_tows**power), optional cap
    ANOM_CLIM_MODE: str = getattr(cfg, "ZP_ANOM_CLIM_MODE", "tow")
    # For year_weighted:
    ANOM_CLIM_WEIGHT_POWER: float = float(getattr(cfg, "ZP_ANOM_CLIM_WEIGHT_POWER", 0.5))  # 0.5=sqrt(n)
    ANOM_CLIM_WEIGHT_CAP: int | None = getattr(cfg, "ZP_ANOM_CLIM_WEIGHT_CAP", None)  # e.g., 20 or 30
    ANOM_CLIM_MIN_TOWS_PER_YEAR: int = int(getattr(cfg, "ZP_ANOM_CLIM_MIN_TOWS_PER_YEAR", 3))
    ANOM_CLIM_EPS_STD: float = float(getattr(cfg, "ZP_ANOM_CLIM_EPS_STD", 1e-12))  # Numerical safety

    # ANOM/FIG OPTIONS
    SEASON_ORDER: List[str] = field(default_factory=lambda: ["Winter", "Spring", "Summer", "Fall"])
    LOG_OFFSET: float = 1e-6
    LAG_YEARS: int = getattr(cfg, "NPGO_LAG_YEARS", 0)
    LOG_TRANSFORM: bool = getattr(cfg, "ZP_LOG_TRANSFORM", True)

    ADD_FIG_SUPTITLE: bool = getattr(cfg, "ZP_ADD_FIG_SUPTITLE", True)
    PREFIX_AX_TITLES_WITH_ECOSPACE: bool = getattr(cfg, "ZP_PREFIX_AX_TITLES_WITH_ECOSPACE", True)

    # Match ecosim-style anomaly/scatter knobs where possible
    ANOM_SEASON_CHOICE: str = getattr(cfg, "ZP_SEASON_CHOICE", "Spring")
    ANOM_PLOT_TYPE: str = getattr(cfg, "ZP_PLOT_TYPE", "bar")  # "bar" or "line"
    SHOW_COUNTS: bool = getattr(cfg, "ZP_SHOW_CNTS", True)
    ANOM_ALL_SEASONS: bool = getattr(cfg, "ZP_ANOM_ALL_SEASONS", False)
    MAKE_SCATTER: bool = getattr(cfg, "ZP_MAKE_SCATTER", True)
    SCATTER_LOG10: bool = getattr(cfg, "ZP_SCATTER_LOG10", True)

    # Tow filtering (Ecosim parity)
    TOW_FILTER_MODE: str = getattr(cfg, "ZP_TOW_FILTER_MODE", "none")
    TOW_MIN_PROP: float = float(getattr(cfg, "ZP_TOW_MIN_PROP", 0.7))
    TOW_MIN_START_DEPTH_M: float = float(getattr(cfg, "ZP_TOW_MIN_START_DEPTH_M", 150.0))
    TOW_MAX_BOTTOM_DEPTH_M: Optional[float] = getattr(cfg, "ZP_TOW_MAX_BOTTOM_DEPTH_M", None)

# ============================================================
# Stage 1 — Match Zooplankton obs to Ecospace (ex-5b)
# ============================================================
cfg5 = Eval5Config()
def _load_ecospace_times(nc_path: str) -> Tuple[Dict, np.ndarray]:
    with nc.Dataset(nc_path, 'r') as ds:
        time_var = ds.variables['time']
        units = time_var.units
        calendar = getattr(time_var, 'calendar', 'standard')

        # Get times (may be Python datetimes or cftime objects)
        times = nc.num2date(time_var[:], units=units, calendar=calendar)

        # Coerce cftime -> Python datetime when possible
        def to_py_datetime(t):
            # cftime objects have year/month/day/hour/minute/second attributes
            if hasattr(t, 'year') and not isinstance(t, datetime):
                # Handle fractional seconds if present
                sec = float(getattr(t, 'second', 0))
                micro = int(round((sec - int(sec)) * 1_000_000))
                return datetime(t.year, t.month, t.day, getattr(t, 'hour', 0),
                                getattr(t, 'minute', 0), int(sec), micro)
            return t  # already a datetime

        times = np.array([to_py_datetime(t) for t in np.asarray(times).ravel()])

        return {
            "Dimensions": {dim: len(ds.dimensions[dim]) for dim in ds.dimensions},
            "Variables": list(ds.variables.keys())
        }, times


def apply_tow_filters(
    df: pd.DataFrame,
    *,
    filter_mode: str,
    min_prop: float,
    min_start_depth_m: float,
    max_bottom_depth_m: Optional[float] = None,
) -> pd.DataFrame:
    """
    Mirror Ecosim's optional tow-eligibility filter.

    Keeps tows if:
      (tow_prop >= min_prop) OR (Tow_start_depth.m. >= min_start_depth_m)
    where tow_prop = (|start|-|end|) / Bottom_depth.m.

    If required columns are missing, returns df unchanged.
    """
    mode = str(filter_mode).lower().strip()
    if mode in ("none", "all", "", "no_filter", "unfiltered"):
        return df

    req = {"Tow_start_depth.m.", "Tow_end_depth.m.", "Bottom_depth.m."}
    missing = [c for c in req if c not in df.columns]
    if missing:
        print(f"[5b][WARN] Tow filter requested ({filter_mode}) but missing columns: {missing}. Skipping filter.")
        return df

    d = df.copy()
    d["tow_depth_range"] = d["Tow_start_depth.m."].abs() - d["Tow_end_depth.m."].abs()
    d["tow_prop"] = d["tow_depth_range"] / d["Bottom_depth.m."]

    keep = (d["tow_prop"] >= float(min_prop)) | (d["Tow_start_depth.m."] >= float(min_start_depth_m))
    if max_bottom_depth_m is not None:
        keep = keep & (d["Bottom_depth.m."] <= float(max_bottom_depth_m))

    before = len(d)
    d = d.loc[keep].copy()
    after = len(d)
    print(f"[5b] Tow filter '{filter_mode}': kept {after}/{before} rows ({after/before:.1%}).")
    return d

def match_zoop_to_ecospace(cfg5: Eval5Config) -> str:
    """Produce matched CSV if needed; return path to file."""
    ecospace_nc = os.path.join(cfg5.NC_PATH_OUT, cfg5.ecospace_nc_name)
    zoop_csv = os.path.join(cfg5.ZOOP_P, cfg5.ZOOP_CSV)
    grid_csv = os.path.join(cfg5.BASEMAP_P, cfg5.GRID_MAP_CSV)
    zoop_ecospace_match_csv = os.path

    out_csv = os.path.join(cfg5.EVALOUT_P, cfg5.ZOOP_CSV_MATCH)
    if (not cfg5.RECOMPUTE_MATCH) and os.path.exists(out_csv):
        print(f"[5b] Using existing matched CSV: {out_csv}")
        return out_csv

    print("[5b] Matching zooplankton obs to Ecospace…")
    zoop_df = pd.read_csv(zoop_csv)

    n_raw = len(zoop_df)
    zoop_df = apply_tow_filters(
        zoop_df,
        filter_mode=cfg5.TOW_FILTER_MODE,
        min_prop=cfg5.TOW_MIN_PROP,
        min_start_depth_m=cfg5.TOW_MIN_START_DEPTH_M,
        max_bottom_depth_m=cfg5.TOW_MAX_BOTTOM_DEPTH_M,
    )
    print(f"[5b] Obs rows: raw={n_raw}, after_tow_filter={len(zoop_df)}")
    
    grid_df = pd.read_csv(grid_csv)

    model_lats = grid_df['lat']
    model_lons = grid_df['lon']
    model_rows = grid_df['EWE_row']
    model_cols = grid_df['EWE_col']
    model_deps = grid_df['depth']
    mask_ecospace = model_deps != 0

    _, ecospace_times = _load_ecospace_times(ecospace_nc)

    # Prepare columns
    for var in cfg5.ZOOP_GROUPS:
        zoop_df[f"EWE-{var}"] = np.nan

    # force a Date col (mid‑month) if not present
    if 'Date' not in zoop_df.columns:
        zoop_df['Date'] = pd.to_datetime(dict(year=zoop_df['Year'], month=zoop_df['Month'], day=15))

    zoop_df['ecospace_closest_lat'] = np.nan
    zoop_df['ecospace_closest_lon'] = np.nan
    zoop_df['closest_ecospace_time'] = pd.NaT

    with nc.Dataset(ecospace_nc, 'r') as eco_ds:
        for idx, row in zoop_df.iterrows():
            lat, lon, obs_time = row['Latitude.W.'], row['Longitude.N.'], row['Date']
            eco_idx = find_nearest_point_from1D(lon, lat, model_lons, model_lats, mask_ecospace, dx=0.01)
            if eco_idx[0] is np.nan:
                continue

            zoop_df.at[idx, 'ecospace_closest_lat'] = model_lats[eco_idx[0]]
            zoop_df.at[idx, 'ecospace_closest_lon'] = model_lons[eco_idx[0]]
            r, c = model_rows[eco_idx[0]], model_cols[eco_idx[0]]

            closest_time = find_closest_date(ecospace_times, obs_time)
            t_idx = int(np.where(ecospace_times == closest_time)[0][0])
            zoop_df.at[idx, 'closest_ecospace_time'] = closest_time

            for var in cfg5.ZOOP_GROUPS:
                if var in eco_ds.variables:
                    try:
                        val = eco_ds.variables[var][t_idx, r, c]
                        zoop_df.at[idx, f"EWE-{var}"] = np.round(val, 3)
                    except Exception:
                        pass

    zoop_df.to_csv(out_csv, index=False)
    print(f"[5b] Saved: {out_csv}")
    return out_csv


# ============================================================
# Stage 2 — Viz & skill stats (ex-5c)
# ============================================================

def _perry_zero_replacement(df: pd.DataFrame, cols: List[str]) -> None:
    for col in cols:
        nonzero = df[col][df[col] > 0]
        if not nonzero.empty:
            mn = nonzero.min()
            df[col] = df[col].apply(lambda x: np.random.uniform(0, 0.5 * mn) if x == 0 else x)


def _compute_skill_statistics(obs: np.ndarray, mod: np.ndarray) -> Dict[str, float]:
    obs = np.asarray(obs)
    mod = np.asarray(mod)
    msk = ~np.isnan(obs) & ~np.isnan(mod)
    obs, mod = obs[msk], mod[msk]
    if len(obs) == 0:
        return {k: np.nan for k in ["obs_std","mod_std","R","RMSE","Bias","WSS"]}

    obs_std = np.std(obs)
    mod_std = np.std(mod)
    R = np.corrcoef(obs, mod)[0, 1] if len(obs) > 1 else np.nan
    RMSE = float(np.sqrt(np.mean((obs - mod) ** 2)))
    Bias = float(np.mean(mod - obs))
    WSS = 1 - (np.sum((obs - mod)**2) / np.sum((np.abs(obs - np.mean(obs)) + np.abs(mod - np.mean(obs)))**2))
    return {"obs_std": obs_std, "mod_std": mod_std, "R": R, "RMSE": RMSE, "Bias": Bias, "WSS": WSS}


def panel_seasonal_boxplots(
    df_station: pd.DataFrame,
    *,
    groups: list[str],
    cfg5,
    extra_total: str | None = None,
    suffix: str = "",
) -> None:
    """
    Make a multi-panel seasonal Obs vs Model boxplot figure.

    - df_station: output of the station-level aggregation in visualize_and_stats
    - groups: list of *observation* group column names (e.g. ["ZC1-EUP", "ZC2-AMP"])
    - extra_total: optional extra column to append (e.g. "TOT_ZC", "TOT_SOFT")
    - suffix: string used in the output filename (e.g. "ZC", "soft")
    """

    if not groups:
        return

    # groups to actually plot (optionally add total column)
    plot_groups = list(groups)
    if extra_total and extra_total in df_station.columns:
        plot_groups.append(extra_total)

    records = []
    for g in plot_groups:
        obs_col = g

        # figure out the matching model column
        if g == extra_total:
            # built these in visualize_and_stats:
            #   "TOT_ZC"   -> "EWE-TOT_ZC"
            #   "TOT_SOFT" -> "EWE-TOT_SOFT"
            mod_col = f"EWE-{extra_total}"
        else:
            mod_col = f"EWE-{g}"

        if obs_col not in df_station.columns or mod_col not in df_station.columns:
            continue

        sub = df_station[["Season", obs_col, mod_col]].copy()

        if cfg5.LOG_TRANSFORM:
            sub["Obs"] = np.log10(sub[obs_col] + cfg5.LOG_OFFSET)
            sub["Model"] = np.log10(sub[mod_col] + cfg5.LOG_OFFSET)
        else:
            sub["Obs"] = sub[obs_col]
            sub["Model"] = sub[mod_col]

        long = sub.melt(
            id_vars="Season",
            value_vars=["Obs", "Model"],
            var_name="Source",
            value_name="plot_value",
        )
        long["group"] = g
        records.append(long)


    if not records:
        return

    plot_df = pd.concat(records, ignore_index=True)

    # order seasons consistently
    plot_df["Season"] = pd.Categorical(
        plot_df["Season"],
        categories=cfg5.SEASON_ORDER,
        ordered=True,
    )

    unique_groups = list(dict.fromkeys(plot_df["group"]))  # preserve order
    n_groups = len(unique_groups)
    ncols = 3
    nrows = int(np.ceil(n_groups / ncols))

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5 * ncols, 4 * nrows),
        sharey=False,
    )
    axes = np.atleast_1d(axes).flatten()

    for i, g in enumerate(unique_groups):
        ax = axes[i]
        dd = plot_df[plot_df["group"] == g]

        sns.boxplot(
            data=dd,
            x="Season",
            y="plot_value",
            showfliers=False,  # match Ecosim behaviour
            hue="Source",
            ax=ax,
        )

        # nicer labels
        if g == "TOT_ZC":
            title = "Total crustaceans (ZC)"
        elif g == "TOT_SOFT":
            title = "Total soft-bodied/other"
        else:
            title = g

        ax.set_title(title)
        ax.set_xlabel("Season")
        ylabel = (
            "log10(g C m$^{-2}$ + offset)"
            if cfg5.LOG_TRANSFORM
            else "Biomass (g C m$^{-2}$)"
        )
        ax.set_ylabel(ylabel)

        if i == 0:
            ax.legend(title="")
        else:
            # remove duplicate legends
            if ax.get_legend() is not None:
                ax.get_legend().remove()

        ax.grid(True, alpha=0.3)

    # turn off any unused axes
    for j in range(n_groups, len(axes)):
        axes[j].axis("off")

    if cfg5.ADD_FIG_SUPTITLE:

        # pick a sensible year range
        if "Year" in df_station.columns and df_station["Year"].notna().any():
            yr0 = int(df_station["Year"].min())
            yr1 = int(df_station["Year"].max())
        else:
            yr0, yr1 = cfg5.START_YEAR, cfg5.END_YEAR

        log_state = "log10" if cfg5.LOG_TRANSFORM else "raw"
        st = fig.suptitle(
            _fig_title(
                cfg5,
                kind="Seasonal biomass (Obs vs Model)",
                years=(yr0, yr1),
                suffix=suffix,
                log_state=log_state,
            ),
            y=0.995,  # <-- was 1.02
            fontsize=14,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.95])  # <-- a touch more room than 0.96
    else:
        st = None
        fig.tight_layout()

    out_plt = os.path.join(
        cfg5.FIGSOUT_P,
        f"ecospace_{cfg5.ecospace_code}_seasB_modobs_panel_{suffix}.png",
    )

    fig.savefig(
        out_plt,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.2,
        bbox_extra_artists=([st] if st is not None else None),
    )

    plt.show()
    plt.close(fig)
    print(f"[5c] Saved seasonal panel boxplots: {out_plt}")


def plot_seasonal_boxpanels(
    df_station: pd.DataFrame,
    groups: list[str],
    *,
    suffix: str,
    title_prefix: str = ""
) -> None:
    """
    Panel of Obs vs Model seasonal boxplots for a set of groups.

    df_station must contain, for each group G in `groups`,
    columns:
      - G         (observed biomass)
      - 'EWE-G'   (model biomass)
      - 'Season'  (categorical)
    """
    import math

    ncols = 3
    nrows = int(math.ceil(len(groups) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5 * ncols, 4 * nrows),
        sharey=False,
    )
    axes = np.atleast_1d(axes).flatten()

    for i, group in enumerate(groups):
        ax = axes[i]

        model_col = f"EWE-{group}"
        if group not in df_station.columns or model_col not in df_station.columns:
            ax.axis("off")
            continue

        tmp = df_station[["Season", group, model_col]].copy()
        tmp = tmp.melt(
            id_vars="Season",
            value_vars=[group, model_col],
            var_name="Source",
            value_name="Biomass",
        )

        # Consistent season ordering
        tmp["Season"] = pd.Categorical(
            tmp["Season"],
            categories=cfg5.SEASON_ORDER,
            ordered=True,
        )

        tmp["Source"] = tmp["Source"].map({group: "Obs", model_col: "Model"})

        sns.boxplot(
            data=tmp,
            x="Season",
            y="Biomass",
            hue="Source",
            ax=ax,
            showfliers=False,
        )
        # ax.set_yscale("log")
        ax.grid(True, which="major", axis="both", linestyle="-", alpha=0.4)

        # nicer titles for the totals
        display_name = {
            "TOT_ZC": "Total (ZC groups)",
            "TOT_SOFT": "Total (soft groups)",
        }.get(group, group)
        ax.set_title(display_name)

        ax.set_xlabel("")
        if i % ncols == 0:
            ax.set_ylabel("Biomass (g C m$^{-2}$)")
        else:
            ax.set_ylabel("")

        # shared legend only on first panel
        if i == 0 and ax.get_legend() is not None:
            ax.legend()
        else:
            if ax.get_legend() is not None:
                ax.get_legend().remove()

    # hide any unused axes
    for j in range(len(groups), len(axes)):
        axes[j].axis("off")

    if title_prefix:
        fig.suptitle(title_prefix, y=1.02)

    if cfg5.ADD_FIG_SUPTITLE:
        # pick a sensible year range
        if "Year" in df_station.columns and df_station["Year"].notna().any():
            yr0 = int(df_station["Year"].min())
            yr1 = int(df_station["Year"].max())
        else:
            yr0, yr1 = cfg5.START_YEAR, cfg5.END_YEAR

        log_state = "log10" if cfg5.LOG_TRANSFORM else "raw"
        fig.suptitle(
            _fig_title(
                cfg5,
                kind="Seasonal biomass (Obs vs Model)",
                years=(yr0, yr1),
                suffix=suffix,
                log_state=log_state,
            ),
            y=1.02,
            fontsize=14,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.96])  # leave room for suptitle
    else:
        fig.tight_layout()

    os.makedirs(cfg5.FIGSOUT_P, exist_ok=True)
    out_plt = os.path.join(
        cfg5.FIGSOUT_P,
        f"ecospace_{cfg5.ecospace_code}_seasB_modobs_panel_{suffix}.png",
    )
    fig.savefig(out_plt, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[5d] Saved seasonal boxplot panels → {out_plt}")



def visualize_and_stats(matched_csv: str, cfg5: Eval5Config) -> str:
    print("[5c] Visualizing & computing skill stats…")
    df = pd.read_csv(matched_csv)

    obs_cols = [c for c in cfg5.ZOOP_GROUPS if c in df.columns]
    model_cols = [f"EWE-{c}" for c in obs_cols]

    _perry_zero_replacement(df, obs_cols)

    # Aggregate by closest_ecospace_time + Index
    df_agg = df.groupby(['closest_ecospace_time', 'Index'])[obs_cols + model_cols].mean().reset_index()

    # Totals
    df_agg['TOT'] = df_agg[obs_cols].sum(axis=1)
    df_agg['EWE-TOT'] = df_agg[model_cols].sum(axis=1)

    # Station-level aggregation
    df_station = df.groupby(['Station', 'Year', 'Season', 'Region'])[obs_cols + model_cols].mean().reset_index()
    df_station['TOT'] = df_station[obs_cols].sum(axis=1)
    df_station['EWE-TOT'] = df_station[model_cols].sum(axis=1)

    # --- NEW: ZC vs soft-bodied/group splits -----------------------------
    crust_cols = [g for g in obs_cols if g.startswith("ZC")]
    soft_cols = [g for g in obs_cols if not g.startswith("ZC")]

    if crust_cols:
        df_station["TOT_ZC"] = df_station[crust_cols].sum(axis=1)
        df_station["EWE-TOT_ZC"] = df_station[[f"EWE-{g}" for g in crust_cols]].sum(axis=1)

    if soft_cols:
        df_station["TOT_SOFT"] = df_station[soft_cols].sum(axis=1)
        df_station["EWE-TOT_SOFT"] = df_station[[f"EWE-{g}" for g in soft_cols]].sum(axis=1)

    # Values used for skill stats (log or raw, controlled by cfg5.LOG_TRANSFORM)
    for col in obs_cols + ["TOT"]:
        mcol = f"EWE-{col}" if col != "TOT" else "EWE-TOT"

        if cfg5.LOG_TRANSFORM:
            df_station[f"skill_obs_{col}"] = np.log10(df_station[col] + cfg5.LOG_OFFSET)
            df_station[f"skill_mod_{col}"] = np.log10(df_station[mcol] + cfg5.LOG_OFFSET)
        else:
            df_station[f"skill_obs_{col}"] = df_station[col]
            df_station[f"skill_mod_{col}"] = df_station[mcol]

    # (Optional) Perry-style climatology on the same transformed scale
    # Only needed  for station-level anomalies for something later.
    if cfg5.LOG_TRANSFORM:
        value_cols_obs = [f"skill_obs_{c}" for c in obs_cols + ["TOT"]]
        value_cols_mod = [f"skill_mod_{c}" for c in obs_cols + ["TOT"]]

        seasonal_means = (
            df_station
            .groupby("Season")[value_cols_obs + value_cols_mod]
            .mean()
            .reset_index()
        )
        annual_means = seasonal_means.mean(numeric_only=True)

        df_clim_anom = df_station.copy()
        for col in obs_cols + ["TOT"]:
            ocol = f"skill_obs_{col}"
            mcol = f"skill_mod_{col}"
            df_clim_anom[f"anom_annmean_{ocol}"] = df_clim_anom[ocol] - annual_means[ocol]
            df_clim_anom[f"anom_annmean_{mcol}"] = df_clim_anom[mcol] - annual_means[mcol]
    else:
        df_clim_anom = df_station.copy()  # placeholder, not used further


    # Compute skill
    rows = []
    for col in obs_cols + ["TOT"]:
        ovals = df_station[f"skill_obs_{col}"].values
        mvals = df_station[f"skill_mod_{col}"].values
        s = _compute_skill_statistics(ovals, mvals)
        s["Group"] = col
        rows.append(s)
    skill_df = pd.DataFrame(rows).sort_values("WSS", ascending=False)

    suffix = "log10" if cfg5.LOG_TRANSFORM else "raw"
    out_skill = os.path.join(
        cfg5.EVALOUT_P,
        f"model_skill_stats_{suffix}_{cfg5.ecospace_code}.csv",
    )
    skill_df.to_csv(out_skill, index=False)
    print(f"[5c] Saved skill stats: {out_skill}")


    # figures
    if cfg5.MAKE_VIZ:
        # Total scatter by Region
        fig, ax = plt.subplots(figsize=(6, 6))

        if cfg5.LOG_TRANSFORM:
            x = df_station["skill_obs_TOT"]
            y = df_station["skill_mod_TOT"]
            xlab = "log10(Observed + offset)"
            ylab = "log10(Modelled + offset)"
        else:
            x = df_station["TOT"]
            y = df_station["EWE-TOT"]
            xlab = "Observed (g C m$^{-2}$)"
            ylab = "Modelled (g C m$^{-2}$)"

        sns.scatterplot(
            x=x,
            y=y,
            hue=df_station["Region"],
            alpha=0.7,
            ax=ax,
        )

        # 1:1 line + sensible limits
        lo = np.nanmin([x.min(), y.min()])
        hi = np.nanmax([x.max(), y.max()])
        pad = 0.05 * (hi - lo if hi > lo else 1.0)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k--", lw=1)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)

        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(
            f"Ecospace {cfg5.ecospace_code}- Model vs. Obs Total Biomass by Station"
            + (" (log10)" if cfg5.LOG_TRANSFORM else "")
        )
        plt.tight_layout()
        out_plt = os.path.join(
            cfg5.FIGSOUT_P,
            f"ecospace_{cfg5.ecospace_code}_scatter_modobs_total.png",
        )
        plt.savefig(out_plt)
        plt.show()
        plt.close()

        # NEW: multi-panel seasonal boxplots (log scale), same aesthetics
        # NEW: multi-panel seasonal boxplots (log10), split into
        #      (a) crustacean ZC groups + Total ZC
        #      (b) soft-bodied/other groups + Total soft
        if crust_cols:
            panel_seasonal_boxplots(
                df_station,
                groups=crust_cols,
                cfg5=cfg5,
                extra_total="TOT_ZC",
                suffix="ZC",
            )

        if soft_cols:
            panel_seasonal_boxplots(
                df_station,
                groups=soft_cols,
                cfg5=cfg5,
                extra_total="TOT_SOFT",
                suffix="soft",
            )

    return out_skill


# ============================================================
# Stage 3 — Anomaly comparison & NPGO (ex-5d)
# ============================================================
# ============================================================
# Stage 3 helpers — reshape matched CSV and compute anomalies
# (aligned with Ecosim 4a pipeline)
# ============================================================

def _finalize_and_save(fig, outpath: str, *, cfg5, st=None):
    """
    Standard figure finalization to avoid suptitle clipping.
    Save first, then show/close.
    """
    if st is not None:
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    else:
        fig.tight_layout()

    fig.savefig(
        outpath,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.2,
        bbox_extra_artists=([st] if st is not None else None),
    )
    plt.show()
    plt.close(fig)


def _fig_title(cfg5, *, kind: str, years: tuple[int, int] | None = None,
               season: str | None = None, suffix: str | None = None,
               label: str | None = None, plot_type: str | None = None,
               log_state: str | None = None) -> str:
    bits = [f"Ecospace {cfg5.ecospace_code}", kind]
    if years:
        bits.append(f"{years[0]}–{years[1]}")
    if season:
        bits.append(season)
    if plot_type:
        bits.append(plot_type)
    if log_state:
        bits.append(log_state)
    if suffix:
        bits.append(suffix)
    if label:
        bits.append(label)
    return " | ".join(bits)


def build_paired_long_ecospace(df_match: pd.DataFrame,
                               *,
                               groups: List[str]) -> pd.DataFrame:
    """
    Turn the matched Ecospace CSV into a long table with
    obs + model biomass for each group and for Total.

    Expected columns in df_match:
      - 'Index'
      - 'Season' (string; Winter/Spring/Summer/Fall)
      - 'Date' (or 'date')
      - one column per obs group, e.g. 'ZC1-EUP'
      - one column per model group, e.g. 'EWE-ZC1-EUP'
    """
    df = df_match.copy()

    # Normalise date/season column names
    if "date" in df.columns:
        df["date"] = pd.to_datetime(df["date"])
    elif "Date" in df.columns:
        df["date"] = pd.to_datetime(df["Date"])
    else:
        raise ValueError("Matched CSV needs a 'Date' or 'date' column.")

    if "season" in df.columns:
        df["season"] = df["season"].astype(str).str.capitalize()
    elif "Season" in df.columns:
        df["season"] = df["Season"].astype(str).str.capitalize()
    else:
        raise ValueError("Matched CSV needs a 'Season' or 'season' column.")

    if "Index" not in df.columns:
        raise ValueError("Matched CSV must have an 'Index' column (tow ID).")

    recs: list[pd.DataFrame] = []

    # Groups + Total
    for g in groups + ["Total"]:
        if g == "Total":
            # Sum obs + model across all zoop groups
            df["Total"] = df[groups].sum(axis=1)
            df["EWE-Total"] = df[[f"EWE-{x}" for x in groups]].sum(axis=1)
            obs_col, mod_col = "Total", "EWE-Total"
        else:
            obs_col, mod_col = g, f"EWE-{g}"
            if obs_col not in df.columns or mod_col not in df.columns:
                # Skip groups that didn't make it into the matched file
                continue

        o = (
            df[["date", "season", "Index", obs_col]]
            .rename(columns={obs_col: "obs_biomass"})
        )
        m = (
            df[["date", "season", "Index", mod_col]]
            .rename(columns={mod_col: "model_biomass"})
        )

        recs.append(
            o.merge(m, on=["date", "season", "Index"])
             .assign(group=g)
        )

    paired = pd.concat(recs, ignore_index=True)
    paired["year"] = paired["date"].dt.year
    return paired


def _wmean(x: np.ndarray, w: np.ndarray) -> float:
    wsum = float(np.sum(w))
    return float(np.sum(w * x) / wsum) if wsum > 0 else np.nan


def _wstd_unbiased(x: np.ndarray, w: np.ndarray) -> float:
    """
    Unbiased-ish weighted std using the effective-dof correction:
      Var = sum(w*(x-mu)^2) / (sum(w) - sum(w^2)/sum(w))
    Works well for "reliability"/effort weights too (and avoids the cap dominating).
    """
    x = np.asarray(x, dtype=float)
    w = np.asarray(w, dtype=float)

    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x, w = x[m], w[m]
    if x.size < 2:
        return np.nan

    mu = _wmean(x, w)
    sw = float(np.sum(w))
    sw2 = float(np.sum(w * w))
    denom = sw - (sw2 / sw) if sw > 0 else 0.0
    if denom <= 0:
        return np.nan

    var = float(np.sum(w * (x - mu) ** 2) / denom)
    return float(np.sqrt(var))


def compute_anomalies_paired(
    paired: pd.DataFrame,
    *,
    season: Optional[str],
    year_range: tuple[int, int],
    log_transform: bool,
    # --- NEW ---
    clim_mode: str = "tow",                 # "tow" | "year_equal" | "year_weighted"
    log_offset: float = 1e-6,
    year_weight_power: float = 0.5,         # used only for year_weighted
    year_weight_cap: Optional[int] = None,  # used only for year_weighted
    min_tows_per_year: int = 3,
    eps_std: float = 1e-12,
) -> pd.DataFrame:
    """
    Returns: columns ['year','group','obs_anom','model_anom'].

    clim_mode:
      - "tow":       climatology from all tows pooled (your current behavior)
      - "year_equal": climatology from annual means; each year gets weight=1
      - "year_weighted": climatology from annual means; year weight ~ (n_tows**power) with optional cap
    """
    y0, y1 = year_range

    df = paired.copy()
    df = df.query("@y0 <= year <= @y1").copy()
    if season is not None:
        df = df.query("season == @season").copy()

    season_label = season if season is not None else "All"
    if df.empty:
        raise ValueError(f"No data in paired table for season={season_label} and years {y0}-{y1}")

    # Transform on the same scale you want anomalies on
    if log_transform:
        df["obs_biomass"] = np.log(df["obs_biomass"] + log_offset)
        df["model_biomass"] = np.log(df["model_biomass"] + log_offset)

    clim_mode = str(clim_mode).strip().lower()

    # ------------------------------------------------------------------
    # Mode A: tow-weighted climatology (CURRENT behavior)
    # ------------------------------------------------------------------
    print("computing anomalies")
    if clim_mode == "tow":
        clim = (
            df.groupby("group", as_index=False)
              .agg(
                  mean_obs=("obs_biomass", "mean"),
                  std_obs=("obs_biomass", "std"),
                  mean_mod=("model_biomass", "mean"),
                  std_mod=("model_biomass", "std"),
              )
        )

        # Guard against zero/degenerate std
        clim.loc[clim["std_obs"].abs() < eps_std, "std_obs"] = np.nan
        clim.loc[clim["std_mod"].abs() < eps_std, "std_mod"] = np.nan

        df = df.merge(clim, on="group", how="left")
        df["obs_anom"] = (df["obs_biomass"] - df["mean_obs"]) / df["std_obs"]
        df["model_anom"] = (df["model_biomass"] - df["mean_mod"]) / df["std_mod"]

        annual = (
            df.groupby(["year", "group"], as_index=False)
              .agg(
                  obs_anom=("obs_anom", "mean"),
                  model_anom=("model_anom", "mean"),
              )
        )
        return annual

    # ------------------------------------------------------------------
    # Mode B/C: year-based climatology (climatology computed from annual means)
    # ------------------------------------------------------------------
    # Annual means first (this is what you actually plot/analyze)
    if "Index" in df.columns:
        n_tows = df.groupby(["year", "group"])["Index"].nunique()
    else:
        n_tows = df.groupby(["year", "group"]).size()

    annual_means = (
        df.groupby(["year", "group"], as_index=False)
          .agg(
              obs_mean=("obs_biomass", "mean"),
              mod_mean=("model_biomass", "mean"),
          )
    )
    annual_means["n_tows"] = n_tows.reset_index(drop=True).values

    # Use only "reasonably estimated" years to define the climatology baseline
    clim_base = annual_means[annual_means["n_tows"] >= int(min_tows_per_year)].copy()
    if clim_base.empty:
        raise ValueError(
            f"No years have >= {min_tows_per_year} tows for season={season_label}. "
            "Lower min_tows_per_year or use clim_mode='tow'."
        )

    rows = []
    for grp, gdf in clim_base.groupby("group"):
        xo = gdf["obs_mean"].to_numpy(dtype=float)
        xm = gdf["mod_mean"].to_numpy(dtype=float)

        if clim_mode == "year_equal":
            mean_obs = float(np.nanmean(xo))
            std_obs = float(np.nanstd(xo, ddof=1)) if np.isfinite(xo).sum() > 1 else np.nan
            mean_mod = float(np.nanmean(xm))
            std_mod = float(np.nanstd(xm, ddof=1)) if np.isfinite(xm).sum() > 1 else np.nan

        elif clim_mode == "year_weighted":
            ny = gdf["n_tows"].to_numpy(dtype=float)
            if year_weight_cap is not None:
                ny = np.minimum(ny, float(year_weight_cap))
            w = np.power(ny, float(year_weight_power))

            mean_obs = _wmean(xo, w)
            std_obs = _wstd_unbiased(xo, w)
            mean_mod = _wmean(xm, w)
            std_mod = _wstd_unbiased(xm, w)

        else:
            raise ValueError(
                f"Unknown clim_mode='{clim_mode}'. Use one of: "
                "'tow', 'year_equal', 'year_weighted'."
            )

        if np.isfinite(std_obs) and abs(std_obs) < eps_std:
            std_obs = np.nan
        if np.isfinite(std_mod) and abs(std_mod) < eps_std:
            std_mod = np.nan

        rows.append(
            dict(
                group=grp,
                mean_obs=mean_obs,
                std_obs=std_obs,
                mean_mod=mean_mod,
                std_mod=std_mod,
            )
        )

    clim = pd.DataFrame(rows)

    out = annual_means.merge(clim, on="group", how="left")
    out["obs_anom"] = (out["obs_mean"] - out["mean_obs"]) / out["std_obs"]
    out["model_anom"] = (out["mod_mean"] - out["mean_mod"]) / out["std_mod"]

    return out[["year", "group", "obs_anom", "model_anom"]]


# ───────────────────────── metrics helper (ecosim-style) ─────────────────────
def compute_stats(df: pd.DataFrame, obs_col: str, mod_col: str, *, log_or_anom: bool = False) -> Dict[str, float]:
    """Return classic skill metrics for two columns in *df* (NaN-robust)."""
    o_ser, m_ser = df[obs_col], df[mod_col]
    mask = o_ser.notna() & m_ser.notna()
    o, m = o_ser[mask].values, m_ser[mask].values
    N = len(o)
    if N == 0:
        return dict(N=0, MB=np.nan, MAE=np.nan, RMSE=np.nan, NRMSE=np.nan,
                    r=np.nan, R2=np.nan, MAPE=np.nan, NSE=np.nan, WSS=np.nan)

    mb = float(np.nanmean(m - o))
    mae = float(np.nanmean(np.abs(m - o)))
    rmse = float(np.sqrt(np.nanmean((m - o) ** 2)))

    # For anomalies/logs, NRMSE isn't meaningful as RMSE / mean(obs)
    if log_or_anom:
        nrmse = rmse
    else:
        denom = float(np.nanmean(o))
        nrmse = rmse / denom if denom != 0 else np.nan

    r = float(np.corrcoef(m, o)[0, 1]) if N > 1 else np.nan
    r2 = float(r ** 2) if not np.isnan(r) else np.nan
    mape = float(np.nanmean(np.abs((m - o) / o)) * 100) if np.all(o != 0) else np.nan

    denom_nse = float(np.nansum((o - np.nanmean(o)) ** 2))
    nse = float(1 - np.nansum((m - o) ** 2) / denom_nse) if denom_nse != 0 else np.nan

    denom_wss = float(np.nansum((np.abs(m - np.nanmean(o)) + np.abs(o - np.nanmean(o))) ** 2))
    wss = float(1 - np.nansum((m - o) ** 2) / denom_wss) if denom_wss != 0 else np.nan

    return dict(N=N, MB=mb, MAE=mae, RMSE=rmse, NRMSE=nrmse,
                r=r, R2=r2, MAPE=mape, NSE=nse, WSS=wss)


def skill_metrics_by_group(annual: pd.DataFrame, *, groups: List[str], log_or_anom: bool) -> pd.DataFrame:
    """Compute skill metrics per group using annual columns obs_anom/model_anom."""
    rows: list[dict] = []
    for g in groups:
        d = annual[annual["group"] == g]
        if len(d) > 1:
            s = compute_stats(d, "obs_anom", "model_anom", log_or_anom=log_or_anom)
            s["Group"] = g
            rows.append(s)
    return pd.DataFrame(rows)


def counts_by_season_year_from_paired(
    paired: pd.DataFrame,
    *,
    season: Optional[str],
    group_for_counts: str = "Total",
    unique_index: bool = True,
    valid_pairs_only: bool = True,
) -> pd.DataFrame:
    """
    Compute a *single* tow/sample count per (year, season) suitable for annotating
    the anomaly panels.

    By default, counts are based on the 'Total' group and unique Index values.
    This tends to match the "how many tows went into that year/season mean"
    interpretation, and avoids double-counting if the matched table has duplicates.
    """
    df = paired.copy()
    if season is not None:
        df = df.query("season == @season").copy()

    season_label = season if season is not None else "All"

    if group_for_counts is not None and "group" in df.columns:
        df = df[df["group"] == group_for_counts]

    if valid_pairs_only:
        df = df.dropna(subset=["obs_biomass", "model_biomass"])

    if df.empty:
        return pd.DataFrame({"year": [], "season": [], "n_tows": []})

    if unique_index and "Index" in df.columns:
        c = df.groupby("year")["Index"].nunique()
    else:
        c = df.groupby("year").size()

    out = c.rename("n_tows").reset_index()
    out["season"] = season_label
    return out


def plot_anomaly_panel_paired(
    annual: pd.DataFrame,
    *,
    cfg5: Eval5Config,
    plot_groups: List[str],
    label: str,
    season: str,
    plot_type: str = "bar",
    counts: Optional[pd.DataFrame] = None,
    show_counts: bool = True,
    npgo_ann: Optional[pd.Series] = None,
    count_yoffset: float = 0.05,
    count_fontsize: int = 9,
) -> str:
    """
    Ecosim-style anomaly panel: one subplot per group, with Obs vs Model shown
    together (bars or lines) and optional per-year tow counts.
    """
    years = sorted(annual["year"].unique())
    if not years:
        raise ValueError("No years found in annual anomalies table.")

    # year -> n_tows map (for labels)
    counts_map: Dict[int, float] = {}
    if show_counts and counts is not None and not counts.empty:
        c = counts.copy()
        # normalize columns
        if "Year" in c.columns and "year" not in c.columns:
            c = c.rename(columns={"Year": "year"})
        if "Season" in c.columns and "season" not in c.columns:
            c = c.rename(columns={"Season": "season"})
        if "n_tows" not in c.columns:
            for alt in ("n", "count", "n_rows", "N"):
                if alt in c.columns:
                    c = c.rename(columns={alt: "n_tows"})
                    break

        if "season" in c.columns:
            c["season"] = c["season"].astype(str).str.capitalize()
            c = c[c["season"] == season.capitalize()]

        if {"year", "n_tows"}.issubset(c.columns):
            counts_map = dict(zip(c["year"].astype(int), c["n_tows"]))
        else:
            print("[WARN] Counts table missing required columns; skipping annotations.")

    ncols = 3
    nrows = int(np.ceil(len(plot_groups) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        sharey=True,
        figsize=(5 * ncols, 4 * nrows),
    )
    axes = np.atleast_1d(axes).flatten()

    annual_idx = annual.set_index(["group", "year"])
    x = np.arange(len(years))
    w = 0.35

    for i, grp in enumerate(plot_groups):
        ax = axes[i]
        y_obs = [
            annual_idx.loc[(grp, yr), "obs_anom"] if (grp, yr) in annual_idx.index else np.nan
            for yr in years
        ]
        y_mod = [
            annual_idx.loc[(grp, yr), "model_anom"] if (grp, yr) in annual_idx.index else np.nan
            for yr in years
        ]

        ax.grid(True, which="major", axis="both", linestyle="-", alpha=0.5, zorder=1)

        if plot_type.lower() == "bar":
            ax.bar(x - w / 2, y_obs, w, label="Obs", zorder=2)
            ax.bar(x + w / 2, y_mod, w, label="Model", zorder=2)
        else:
            ax.plot(years, y_obs, "-o", label="Obs", zorder=2)
            ax.plot(years, y_mod, "-o", label="Model", zorder=2)

        # Optional NPGO overlay (same axis, like the current ecospace plot)
        if npgo_ann is not None:
            x_npgo = [yr for yr in years if yr in npgo_ann.index]
            y_npgo = [npgo_ann.loc[yr] for yr in x_npgo]
            if len(x_npgo) > 0:
                ax.plot(
                    [years.index(yr) for yr in x_npgo],
                    y_npgo,
                    "k-",
                    label=f"NPGO (lag {cfg5.LAG_YEARS}y)",
                )

        ax.set_xticks(x)
        ax.set_xticklabels(years, rotation=45)
        ax.set_title(grp)

        if i % ncols == 0:
            ax.set_ylabel("Anomaly (z-score)")

        if i == 0:
            ax.legend()

        # Count annotations
        if counts_map:
            ax.margins(y=0.15)
            ymin, ymax = ax.get_ylim()
            yrange = ymax - ymin if ymax > ymin else 1.0
            for j, yr in enumerate(years):
                n = counts_map.get(int(yr))
                if n is None:
                    continue
                vals = []
                if j < len(y_obs):
                    vals.append(y_obs[j])
                if j < len(y_mod):
                    vals.append(y_mod[j])
                if len(vals) == 0 or all(pd.isna(v) for v in vals):
                    continue
                y_top = np.nanmax(vals)
                y_text = y_top + yrange * count_yoffset
                ax.text(
                    j,
                    y_text,
                    f"{int(n)}",
                    ha="center",
                    va="bottom",
                    fontsize=count_fontsize,
                )

    # turn off unused axes
    for j in range(len(plot_groups), len(axes)):
        axes[j].axis("off")

    # ---- Figure title (ONE suptitle only; informative) ----
    log_state = "anom z-score"
    yr0, yr1 = int(years[0]), int(years[-1])

    st = None
    if getattr(cfg5, "ADD_FIG_SUPTITLE", True):
        st = fig.suptitle(
            _fig_title(
                cfg5,
                kind="Obs vs Model annual anomalies",
                years=(yr0, yr1),
                season=season,
                label=label,
                plot_type=plot_type,
                log_state=log_state,
            ),
            y=0.995,  # keep inside canvas
            fontsize=14,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.95])  # reserve space for title
    else:
        fig.tight_layout()

    # ---- Save (save before show; tight bbox includes suptitle) ----
    os.makedirs(cfg5.FIGSOUT_P, exist_ok=True)
    out_panel = os.path.join(
        cfg5.FIGSOUT_P,
        f"ecospace_{cfg5.ecospace_code}_zoop_anomalies_panel_{season}_{label}.png",
    )

    fig.savefig(
        out_panel,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.2,
        bbox_extra_artists=([st] if st is not None else None),
    )
    plt.show()
    plt.close(fig)

    print(f"[INFO] Saved anomaly panel: {out_panel}")
    return out_panel


def plot_scatter_matched_tows(
    paired_long: pd.DataFrame,
    *,
    cfg5: Eval5Config,
    groups: List[str],
    label: str,
    season: str,  # <-- ADD THIS (so it can appear in the suptitle)
    outname: str = "scatter_tows",
    log10: bool = True,
) -> str:
    """
    Scatter of tow-level matched values (obs vs model), one subplot per group.
    Mirrors ecosim_eval_5_zoop.py behaviour.
    """
    d = paired_long.copy()
    d = d.dropna(subset=["group", "obs_biomass", "model_biomass"])
    d = d[d["group"].isin(groups)]

    if d.empty:
        raise ValueError("No paired tow-level values available for scatter plot.")

    if log10:
        eps = 1e-12
        d = d[(d["obs_biomass"] > 0) & (d["model_biomass"] > 0)]
        d["obs_x"] = np.log10(d["obs_biomass"] + eps)
        d["mod_y"] = np.log10(d["model_biomass"] + eps)
        xlab = "Obs (log10 g C m$^{-2}$)"
        ylab = "Model (log10 g C m$^{-2}$)"
    else:
        d["obs_x"] = d["obs_biomass"]
        d["mod_y"] = d["model_biomass"]
        xlab = "Obs (g C m$^{-2}$)"
        ylab = "Model (g C m$^{-2}$)"

    ncols = 3
    nrows = int(np.ceil(len(groups) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for i, grp in enumerate(groups):
        ax = axes[i]
        g = d[d["group"] == grp]
        ax.scatter(g["obs_x"], g["mod_y"], s=12, alpha=0.6)

        if len(g) > 0:
            lo = min(g["obs_x"].min(), g["mod_y"].min())
            hi = max(g["obs_x"].max(), g["mod_y"].max())
            ax.plot([lo, hi], [lo, hi], "--")

            if len(g) > 1:
                r = np.corrcoef(g["obs_x"], g["mod_y"])[0, 1]
                ax.set_title(f"{grp}  r={r:.2f}")
            else:
                ax.set_title(grp)

        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.grid(True, alpha=0.3)

    for j in range(len(groups), len(axes)):
        axes[j].axis("off")

    # ---- Figure suptitle (informative + consistent) ----
    st = None
    if getattr(cfg5, "ADD_FIG_SUPTITLE", True):
        log_state = "log10" if log10 else "raw"

        st = fig.suptitle(
            _fig_title(
                cfg5,
                kind="Tow-level scatter (Obs vs Model)",
                season=season,           # <-- ADD THIS
                plot_type="scatter",
                log_state=log_state,
                suffix=label,            # keep label as suffix for consistency
            ),
            y=0.995,
            fontsize=14,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    else:
        fig.tight_layout()

    # ---- Save ----
    os.makedirs(cfg5.FIGSOUT_P, exist_ok=True)
    outpath = os.path.join(
        cfg5.FIGSOUT_P,
        f"ecospace_{cfg5.ecospace_code}_zoop_{outname}_{label}.png",
    )
    fig.savefig(
        outpath,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.2,
        bbox_extra_artists=([st] if st is not None else None),
    )
    plt.show()
    plt.close(fig)
    print(f"[INFO] Saved scatter plot: {outpath}")
    return outpath



def plot_scatter_anomalies(
    annual: pd.DataFrame,
    *,
    cfg5: Eval5Config,
    groups: List[str],
    label: str,
    season: str,  # <-- ADD THIS
    outname: str = "scatter_anoms",
) -> str:
    """
    Scatter of annual anomalies (obs vs model), one subplot per group.
    Mirrors ecosim_eval_5_zoop.py behaviour.
    """
    d = annual.dropna(subset=["group", "obs_anom", "model_anom"]).copy()
    d = d[d["group"].isin(groups)]

    if d.empty:
        raise ValueError("No annual anomalies available for anomaly scatter plot.")

    ncols = 3
    nrows = int(np.ceil(len(groups) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for i, grp in enumerate(groups):
        ax = axes[i]
        g = d[d["group"] == grp]
        ax.scatter(g["obs_anom"], g["model_anom"], s=25, alpha=0.7)

        if len(g) > 0:
            lo = min(g["obs_anom"].min(), g["model_anom"].min())
            hi = max(g["obs_anom"].max(), g["model_anom"].max())
            ax.plot([lo, hi], [lo, hi], "--")

            if len(g) > 1:
                r = np.corrcoef(g["obs_anom"], g["model_anom"])[0, 1]
                ax.set_title(f"{grp}  r={r:.2f}")
            else:
                ax.set_title(grp)

        ax.set_xlabel("Obs anomaly (z-score)")
        ax.set_ylabel("Model anomaly (z-score)")
        ax.grid(True, alpha=0.3)

    for j in range(len(groups), len(axes)):
        axes[j].axis("off")

    # ---- Figure suptitle (informative + consistent) ----
    st = None
    if getattr(cfg5, "ADD_FIG_SUPTITLE", True):
        # prefer actual range in the data if present
        if "year" in d.columns and d["year"].notna().any():
            yr0 = int(d["year"].min())
            yr1 = int(d["year"].max())
        elif "Year" in d.columns and d["Year"].notna().any():
            yr0 = int(d["Year"].min())
            yr1 = int(d["Year"].max())
        else:
            # fall back to config
            yr0 = getattr(cfg5, "START_YEAR", None)
            yr1 = getattr(cfg5, "END_YEAR", None)

        years = (yr0, yr1) if (yr0 is not None and yr1 is not None) else None

        st = fig.suptitle(
            _fig_title(
                cfg5,
                kind="Anomaly scatter (Obs vs Model)",
                years=years,
                season=season,          # <-- ADD THIS
                plot_type="scatter",
                log_state="anom z-score",
                suffix=label,
            ),
            y=0.995,
            fontsize=14,
        )
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    else:
        fig.tight_layout()

    # ---- Save (save before show; include suptitle) ----
    os.makedirs(cfg5.FIGSOUT_P, exist_ok=True)
    outpath = os.path.join(
        cfg5.FIGSOUT_P,
        f"ecospace_{cfg5.ecospace_code}_zoop_{outname}_{label}.png",
    )

    fig.savefig(
        outpath,
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.2,
        bbox_extra_artists=([st] if st is not None else None),
    )
    plt.show()
    plt.close(fig)
    print(f"[INFO] Saved anomaly scatter plot: {outpath}")
    return outpath


def anomaly_comparisons(
    cfg5: Eval5Config,
    groups: Optional[List[str]] = None,
    *,
    include_total: bool = True,
    seasons: Optional[List[str]] = None,
) -> Tuple[str, Optional[str]]:
    """
    Create ecosim-style model-vs-obs anomaly panels for Ecospace matched outputs.

    To mirror the Ecosim workflow, this runs *two passes* when possible:
      1) Crustaceans (ZC*), with Total = sum(ZC*)
      2) Soft-bodied (ZS2–ZS4), with Total = sum(ZS2–ZS4)

    """
    print("[5d] Building anomaly comparisons (Model vs Obs)…")

    # --- load matched table
    matched_csv = os.path.join(cfg5.EVALOUT_P, cfg5.ZOOP_CSV_MATCH)
    if not os.path.exists(matched_csv):
        raise FileNotFoundError(
            f"Matched CSV not found: {matched_csv}. Run 5b first or set RECOMPUTE_MATCH=True."
        )

    df_match = pd.read_csv(matched_csv)

    # Determine which groups are actually present in the matched file
    candidate = groups or cfg5.ZOOP_GROUPS
    groups_present = [g for g in candidate if g in df_match.columns and f"EWE-{g}" in df_match.columns]
    if not groups_present:
        raise ValueError(
            "No zoop group columns found in matched CSV. "
            "Expected columns like 'ZC1-EUP' and 'EWE-ZC1-EUP'."
        )

    # Split to match Ecosim passes (same definitions as ecosim_eval_5_zoop.py)
    crust_groups = [g for g in groups_present if str(g).startswith("ZC")]
    zs_groups = [g for g in groups_present if str(g).startswith(("ZS2", "ZS3", "ZS4"))]

    passes: list[tuple[str, list[str]]] = []
    if crust_groups:
        passes.append(("ZC", crust_groups))
    if zs_groups:
        passes.append(("ZS", zs_groups))

    # Optional: warn if you have groups present that won't be included in either pass (e.g., ZF1-ICT)
    extra_groups = [g for g in groups_present if g not in set(crust_groups + zs_groups)]
    if extra_groups:
        print(f"[WARN] Groups present but excluded from anomaly passes to match Ecosim definitions: {extra_groups}")

    # --- optional external NPGO series (spring mean only, lagged)
    npgo_ann = None
    if cfg5.NPGO_CSV and os.path.exists(cfg5.NPGO_CSV):
        try:
            npgo = pd.read_csv(cfg5.NPGO_CSV)
            npgo["date"] = pd.to_datetime(npgo["date"])
            npgo["year"] = npgo["date"].dt.year
            npgo["month"] = npgo["date"].dt.month
            spring = npgo[npgo["month"].isin([3, 4, 5])]
            npgo_ann = (
                spring.groupby("year")["npgo"]
                .mean()
                .shift(cfg5.LAG_YEARS)
                .dropna()
            )
        except Exception as e:
            print(f"[WARN] Failed to load NPGO series from {cfg5.NPGO_CSV}: {e}")
            npgo_ann = None

    # --- seasons to run
    if seasons is not None:
        seasons_to_run = seasons
    else:
        seasons_to_run = cfg5.SEASON_ORDER if cfg5.ANOM_ALL_SEASONS else [cfg5.ANOM_SEASON_CHOICE]

    out_panel_last: Optional[str] = None
    out_total_last: Optional[str] = None

    for pass_label, pass_groups in passes:
        # build paired long table (tow-level) for this pass (Total is computed from pass_groups)
        paired = build_paired_long_ecospace(df_match, groups=pass_groups)

        plot_groups = list(pass_groups) + (["Total"] if include_total else [])
        if "Total" in paired["group"].unique() and "Total" not in plot_groups and include_total:
            plot_groups.append("Total")

        seasons_plus_all = ["All"] + list(seasons_to_run)

        for season in seasons_plus_all:
            season_filter = None if str(season).lower() == "all" else season
            # Annual anomalies
            annual = compute_anomalies_paired(
                paired,
                season=season_filter,
                year_range=(cfg5.START_YEAR, cfg5.END_YEAR),
                log_transform=cfg5.LOG_TRANSFORM,
                clim_mode=cfg5.ANOM_CLIM_MODE,
                log_offset=cfg5.LOG_OFFSET,
                year_weight_power=cfg5.ANOM_CLIM_WEIGHT_POWER,
                year_weight_cap=cfg5.ANOM_CLIM_WEIGHT_CAP,
                min_tows_per_year=cfg5.ANOM_CLIM_MIN_TOWS_PER_YEAR,
                eps_std=cfg5.ANOM_CLIM_EPS_STD,
            )

            # Counts for labels (use Total by default)
            counts = counts_by_season_year_from_paired(
                paired,
                season=season_filter,
                group_for_counts="Total" if "Total" in paired["group"].unique() else None,
                unique_index=True,
                valid_pairs_only=True,
            )

            # Panel plot (Obs vs Model together)
            out_panel_last = plot_anomaly_panel_paired(
                annual,
                cfg5=cfg5,
                plot_groups=plot_groups,
                label=f"zoop_{pass_label}",
                season=season,
                plot_type=cfg5.ANOM_PLOT_TYPE,
                counts=counts,
                show_counts=cfg5.SHOW_COUNTS,
                npgo_ann=npgo_ann,
            )

            # Optional standalone Total plot (matches old behavior)
            if include_total and "Total" in plot_groups:
                years = sorted(annual["year"].unique().tolist())
                tot = annual[annual["group"] == "Total"].set_index("year").reindex(years)

                fig, ax = plt.subplots(figsize=(12, 4))
                ax.axhline(0, lw=1, linestyle="--", alpha=0.6)
                ax.bar(np.arange(len(years)) - 0.18, tot["obs_anom"].values, width=0.36, label="Observed")
                ax.bar(np.arange(len(years)) + 0.18, tot["model_anom"].values, width=0.36, label="Model")
                ax.set_xticks(np.arange(len(years)))
                ax.set_xticklabels([str(y) for y in years], rotation=45, ha="right")
                ax.set_title(f"Ecospace vs Obs annual anomalies — Total ({pass_label}) — {season}")
                ax.set_ylabel("Standardized anomaly")

                # annotate counts once (for Total)
                if cfg5.SHOW_COUNTS and not counts.empty:
                    cmap = dict(zip(counts["year"].astype(int), counts["n_tows"]))
                    ax.margins(y=0.15)
                    ymin, ymax = ax.get_ylim()
                    yrange = ymax - ymin if ymax > ymin else 1.0
                    for j, yr in enumerate(years):
                        n = cmap.get(int(yr))
                        if n is None:
                            continue
                        y_top = np.nanmax([tot["obs_anom"].values[j], tot["model_anom"].values[j]])
                        ax.text(j, y_top + yrange * 0.05, f"{int(n)}", ha="center", va="bottom", fontsize=9)

                ax.legend()
                fig.tight_layout()
                os.makedirs(cfg5.FIGSOUT_P, exist_ok=True)
                out_total_last = os.path.join(
                    cfg5.FIGSOUT_P,
                    f"ecospace_{cfg5.ecospace_code}_zoop_anomalies_total_{pass_label}_{season}.png",
                )
                plt.show()
                fig.savefig(out_total_last, dpi=300)
                plt.close(fig)
                print(f"[INFO] Saved total anomaly plot: {out_total_last}")

            # Skill metric summary (print + CSV), like Ecosim
            try:
                stats_df = skill_metrics_by_group(
                    annual,
                    groups=plot_groups,
                    log_or_anom=True,  # anomalies are dimensionless
                ).sort_values("WSS", ascending=False)
            except Exception as e:
                print(f"[WARN] Skill metrics failed for {pass_label}/{season}: {e}")
                stats_df = pd.DataFrame()

            if not stats_df.empty:
                print(f"\n[{pass_label}] {season} skill metrics (anomalies):")
                print(stats_df.to_string(index=False))

                out_stats = os.path.join(
                    cfg5.EVALOUT_P,
                    f"ecospace_zoop_performance_{cfg5.ecospace_code}_{pass_label}_{season}.csv",
                )
                stats_df.to_csv(out_stats, index=False)
                print(f"[INFO] Saved skill metrics: {out_stats}")

            # Scatter plots (tow-level + anomaly scatter)
            if cfg5.MAKE_SCATTER:
                paired_seas = paired.copy() if season_filter is None else paired.query("season == @season").copy()
                plot_scatter_matched_tows(
                    paired_seas,
                    cfg5=cfg5,
                    groups=plot_groups,
                    label=f"zoop_{pass_label}",
                    outname=f"scatter_tows_{pass_label}_{season}",
                    log10=cfg5.SCATTER_LOG10,
                    season=season,
                )
                plot_scatter_anomalies(
                    annual,
                    cfg5=cfg5,
                    groups=plot_groups,
                    label=f"zoop_{pass_label}",
                    outname=f"scatter_anoms_{pass_label}_{season}",
                    season=season,
                )

    return out_panel_last or "", out_total_last


def run_zoop_eval(
    recompute_match: Optional[bool] = None,
    make_viz: Optional[bool] = None,
    make_anom: Optional[bool] = None,
    groups: Optional[List[str]] = None,
) -> None:
    cfg5 = Eval5Config()
    if recompute_match is not None:
        cfg5.RECOMPUTE_MATCH = recompute_match
    if make_viz is not None:
        cfg5.MAKE_VIZ = make_viz
    if make_anom is not None:
        cfg5.MAKE_ANOM = make_anom

    # Stage 1: matching
    matched = match_zoop_to_ecospace(cfg5)

    # Stage 2: viz & stats
    skill_csv = visualize_and_stats(matched, cfg5)

    # Stage 3: anomaly panels
    if cfg5.MAKE_ANOM:
        anomaly_comparisons(cfg5, groups=groups)

    print("[5] DONE.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--recompute-match", action="store_true",
                        help="Force recomputation of obs–model match")
    parser.add_argument("--no-viz", action="store_true",
                        help="Skip making figures")
    parser.add_argument("--no-anom", action="store_true",
                        help="Skip anomaly calculations")
    parser.add_argument("--groups", nargs="*", default=None,
                        help="Subset of zoop groups to include (e.g. ZC4-CLG ZC5-CSM)")

    args = parser.parse_args()

    run_zoop_eval(
        recompute_match=args.recompute_match if args.recompute_match else False,
        make_viz=not args.no_viz,
        make_anom=not args.no_anom,
        groups=args.groups,
    )