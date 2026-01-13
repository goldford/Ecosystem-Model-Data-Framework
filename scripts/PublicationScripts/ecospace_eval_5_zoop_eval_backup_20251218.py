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

    # ANOM/FIG OPTIONS
    SEASON_ORDER: List[str] = field(default_factory=lambda: ["Winter", "Spring", "Summer", "Fall"])
    LOG_OFFSET: float = 1e-6
    LAG_YEARS: int = getattr(cfg, "NPGO_LAG_YEARS", 0)
    LOG_TRANSFORM: bool = getattr(cfg, "ZP_LOG_TRANSFORM", True)


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
            # we built these in visualize_and_stats:
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

    fig.tight_layout()
    out_plt = os.path.join(
        cfg5.FIGSOUT_P,
        f"ecospace_{cfg5.ecospace_code}_seasB_modobs_panel_{suffix}.png",
    )
    fig.savefig(out_plt, dpi=300)
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
    # Only needed if you want these station-level anomalies for something later.
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
            "Total Biomass by Station"
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

def build_paired_long_ecospace(df_match: pd.DataFrame,
                               *,
                               groups: List[str]) -> pd.DataFrame:
    """
    Turn the matched Ecospace CSV into a long table with
    obs + model biomass for each group and for Total.

    Expected columns in df_match (which you have):
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


def compute_anomalies_paired(
    paired: pd.DataFrame,
    *,
    season: str,
    year_range: tuple[int, int],
    log_transform: bool,
) -> pd.DataFrame:
    """
    Ecosim-style anomalies:

      1. Filter to chosen season + year window.
      2. (Optionally) log-transform obs + model.
      3. For each group:
         - compute mean + std over all years (climatology)
         - anomaly = (value - mean) / std
      4. Average anomalies within each year (if >1 tow/year).

    Returns: columns ['year','group','obs_anom','model_anom'].
    """
    y0, y1 = year_range

    df = (
        paired
        .query("season == @season and @y0 <= year <= @y1")
        .copy()
    )

    if df.empty:
        raise ValueError(f"No data in paired table for season={season} and years {y0}-{y1}")

    if log_transform:
        # Natural log is fine here; using log10 would only rescale,
        # but if you prefer to match Perry-style, you can switch to np.log10.
        df["obs_biomass"] = np.log(df["obs_biomass"] + 1e-6)
        df["model_biomass"] = np.log(df["model_biomass"] + 1e-6)

    clim = (
        df.groupby("group", as_index=False)
          .agg(
              mean_obs=("obs_biomass", "mean"),
              std_obs=("obs_biomass", "std"),
              mean_mod=("model_biomass", "mean"),
              std_mod=("model_biomass", "std"),
          )
    )

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

def anomaly_comparisons(
    cfg5: Eval5Config,
    groups: Optional[List[str]] = None,
    include_total: bool = True,
) -> Tuple[str, Optional[str]]:
    """
    Stage 3 — anomaly comparisons (Ecospace vs Obs), aligned with Ecosim:

    - anomalies computed from the *matched tow-level pairs*
      using compute_anomalies_paired()
    - same z-score style for both obs and model
    - seasonal filter + year range from cfg
    - NPGO spring (MAM) overlay with optional lag
    """
    print("[5d] Building anomaly comparisons (Model vs Obs) …")

    # 1) Decide which groups to use
    if groups is None:
        groups = cfg5.ZOOP_GROUPS[:]  # base zoop groups

    # 2) Load NPGO (Spring mean, with lag)
    npgo_path = os.path.join(cfg5.EVALOUT_P, cfg5.NPGO_CSV)
    df_npgo = pd.read_csv(npgo_path)
    df_npgo.columns = ["date", "npgo"]
    df_npgo["date"] = pd.to_datetime(df_npgo["date"])
    df_npgo = df_npgo[df_npgo["npgo"] > -99]

    df_npgo["year"] = df_npgo["date"].dt.year
    df_npgo["month"] = df_npgo["date"].dt.month

    spring = df_npgo[df_npgo["month"].isin([3, 4, 5])]
    npgo_ann = spring.groupby("year")["npgo"].mean().shift(cfg5.LAG_YEARS)

    # 3) Load matched CSV (Stage 1 output) and build "paired" long table
    matched_csv = os.path.join(
        cfg5.EVALOUT_P,
        f"Zooplankton_matched_to_model_out_{cfg5.ecospace_code}.csv",
    )
    df_match = pd.read_csv(matched_csv)

    paired = build_paired_long_ecospace(df_match, groups=groups)

    # 4) Compute anomalies using the same logic as Ecosim
    #    Here I use Spring, but you can point to a config knob if you prefer.
    season_choice = "Spring"   # or getattr(cfg, "ZP_SEASON_CHOICE", "Spring")
    year_range = (cfg5.START_YEAR, cfg5.END_YEAR)

    annual = compute_anomalies_paired(
        paired,
        season=season_choice,
        year_range=year_range,
        log_transform=cfg5.LOG_TRANSFORM,
    )

    # 5) Build a tidy table similar to your current anom_df:
    #    columns: Year, Group, Anomaly, Source ('Observed'/'Model')
    annual["Year"] = annual["year"]
    annual["Group"] = annual["group"]

    obs_df = (
        annual[["Year", "Group", "obs_anom"]]
        .rename(columns={"obs_anom": "Anomaly"})
        .assign(Source="Observed")
    )
    mod_df = (
        annual[["Year", "Group", "model_anom"]]
        .rename(columns={"model_anom": "Anomaly"})
        .assign(Source="Model")
    )

    anom_df = pd.concat([obs_df, mod_df], ignore_index=True)

    # 6) Decide which groups to plot in panels
    all_groups = groups[:]
    if include_total:
        if "Total" not in all_groups:
            all_groups.append("Total")

    # 7) Panel figure: model vs obs anomalies + NPGO overlay
    fig, axs = plt.subplots(len(all_groups), 2,
                            figsize=(14, 3.5 * len(all_groups)),
                            sharex=True)
    if len(all_groups) == 1:
        axs = [axs]  # make it iterable


    for i, g in enumerate(all_groups):
        for j, source in enumerate(["Model", "Observed"]):
            sub = (
                anom_df[(anom_df["Group"] == g) &
                        (anom_df["Source"] == source)]
                .set_index("Year")
                .sort_index()
            )
            if sub.empty:
                continue

            # align NPGO with the available years
            npgo = npgo_ann[sub.index]

            colors = ["blue" if x > 0 else "red" for x in sub["Anomaly"]]
            axs[i][j].bar(sub.index, sub["Anomaly"], color=colors)
            axs[i][j].plot(npgo.index, npgo.values,
                           color="black", label="NPGO (lagged)")

            axs[i][j].set_xlim(left=year_range[0])
            axs[i][j].axhline(0, color="grey", linestyle="--")
            axs[i][j].set_ylabel("Standardised Anomaly")
            axs[i][j].legend(loc="lower left")

            letter = f"({chr(97 + i)})"  # 97 = 'a', so i=0 -> '(a)', i=1 -> '(b)', etc.

            axs[i][j].text(
                0.05,
                0.95,
                f"{letter} {g} - {source}",
                transform=axs[i][j].transAxes,
                fontsize=12,
                va="top",
            )

    out_panel = os.path.join(
        cfg5.FIGSOUT_P,
        f"ecospace_{cfg5.ecospace_code}_anom_panel_Zoop_vs_Model_with_NPGO.png",
    )
    plt.tight_layout()
    plt.savefig(out_panel)
    plt.show()
    print(f"[5d] Saved panel fig: {out_panel}")

    # 8) Total-only stacked figure, using the *same* anomalies
    out_total = None
    if include_total and "Total" in all_groups:
        tot_obs = (
            anom_df[(anom_df["Group"] == "Total") &
                    (anom_df["Source"] == "Observed")]
            .set_index("Year")
            .sort_index()
        )
        tot_mod = (
            anom_df[(anom_df["Group"] == "Total") &
                    (anom_df["Source"] == "Model")]
            .set_index("Year")
            .sort_index()
        )
        if not tot_obs.empty and not tot_mod.empty:
            npgo = npgo_ann[tot_obs.index]

            fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

            axs[0].bar(
                tot_obs.index,
                tot_obs["Anomaly"],
                color=["blue" if x > 0 else "red"
                       for x in tot_obs["Anomaly"]],
            )
            axs[0].plot(npgo.index, npgo.values,
                        color="black", label="NPGO (lagged)")
            axs[0].axhline(0, color="grey", linestyle="--")
            axs[0].set_ylabel("Observed Anomaly")
            axs[0].set_title("(a) Total Zooplankton - Observed")
            axs[0].legend(loc="lower left")

            axs[1].bar(
                tot_mod.index,
                tot_mod["Anomaly"],
                color=["blue" if x > 0 else "red"
                       for x in tot_mod["Anomaly"]],
            )
            axs[1].plot(npgo.index, npgo.values,
                        color="black", label="NPGO (lagged)")
            axs[1].axhline(0, color="grey", linestyle="--")
            axs[1].set_ylabel("Modelled Anomaly")
            axs[1].set_title("(b) Total Zooplankton - Model")
            axs[1].legend(loc="lower left")
            axs[1].set_xlabel("Year")
            axs[1].set_xlim(left=year_range[0])

            plt.tight_layout()
            out_total = os.path.join(
                cfg5.FIGSOUT_P,
                f"ecospace_{cfg5.ecospace_code}_anom_TOTAL_stacked_Model_vs_Obs.png",
            )
            plt.savefig(out_total)
            plt.show()
            print(f"[5d] Saved total fig: {out_total}")

    return out_panel, out_total

# def anomaly_comparisons(cfg5: Eval5Config, groups: Optional[List[str]] = None, include_total: bool = True) -> Tuple[str, Optional[str]]:
#     print("[5d] Building anomaly comparisons (Model vs Obs) …")
#     groups = groups or [g for g in cfg5.ZOOP_GROUPS if g.startswith('ZC') or g.startswith('ZS')][:5]
#
#     # Load NPGO
#     npgo_path = os.path.join(cfg5.EVALOUT_P, cfg5.NPGO_CSV)
#     df_npgo = pd.read_csv(npgo_path)
#     df_npgo.columns = ["date", "npgo"]
#     df_npgo["date"] = pd.to_datetime(df_npgo["date"])
#     df_npgo = df_npgo[df_npgo["npgo"] > -99]
#     df_npgo["year"] = df_npgo["date"].dt.year
#     df_npgo["month"] = df_npgo["date"].dt.month
#     spring = df_npgo[df_npgo["month"].isin([3,4,5])]
#     npgo_ann = spring.groupby("year")["npgo"].mean().shift(cfg5.LAG_YEARS)
#
#     # Load matched CSV
#     matched_csv = os.path.join(cfg5.EVALOUT_P, cfg5.ZOOP_CSV_MATCH)
#     df = pd.read_csv(matched_csv)
#
#     # Zero replacement for obs
#     for col in groups:
#         nonzero = df[col][df[col] > 0]
#         if not nonzero.empty:
#             mn = nonzero.min()
#             df[col] = df[col].apply(lambda x: np.random.uniform(0, 0.5 * mn) if x == 0 else x)
#
#     # Station aggregation
#     df_station = df.groupby(['Station', 'Year', 'Season'])[groups + [f"EWE-{g}" for g in groups]].mean().reset_index()
#     for col in groups:
#         df_station[f'log10_{col}'] = np.log10(df_station[col] + cfg5.LOG_OFFSET)
#
#     # OBS anomalies (annual mean of station means, z‑scored relative to 1980‑2017)
#     results_obs = []
#     for g in groups:
#         ser = df_station.groupby('Year')[f'log10_{g}'].mean().dropna()
#         ser = ser[(ser.index >= cfg5.START_YEAR) & (ser.index <= cfg5.END_YEAR)]
#         clim = ser.loc[ser.index.isin(range(cfg5.START_YEAR, cfg5.END_YEAR))]
#         anom = (ser - clim.mean()) / clim.std()
#         results_obs.append(pd.DataFrame({'Year': ser.index, 'Group': g, 'Anomaly': anom.values, 'Source': 'Observed'}))
#
#     if include_total:
#         df_station['TOT'] = df_station[groups].sum(axis=1)
#         df_station['log10_TOT'] = np.log10(df_station['TOT'] + cfg5.LOG_OFFSET)
#         ser = df_station.groupby('Year')['log10_TOT'].mean().dropna()
#         ser = ser[(ser.index >= cfg5.START_YEAR) & (ser.index <= cfg5.END_YEAR)]
#         clim = ser.loc[ser.index.isin(range(cfg5.START_YEAR, cfg5.END_YEAR))]
#         anom = (ser - clim.mean()) / clim.std()
#         results_obs.append(pd.DataFrame({'Year': ser.index, 'Group': 'TOTAL', 'Anomaly': anom.values, 'Source': 'Observed'}))
#
#     # MODEL anomalies direct from raw Ecospace
#     ds = xr.open_dataset(os.path.join(cfg5.NC_PATH_OUT, cfg5.ecospace_nc_name))
#     results_mod = []
#     for g in groups:
#         da = ds[g]
#         da_mean = da.mean(dim=('row', 'col'), skipna=True)
#         series = da_mean.groupby('time.year').mean().to_pandas()
#         series = series[(series.index >= cfg5.START_YEAR) & (series.index <= cfg5.END_YEAR)]
#         clim = series.loc[series.index.isin(range(cfg5.START_YEAR, cfg5.END_YEAR))]
#         anom = (series - clim.mean()) / clim.std()
#         results_mod.append(pd.DataFrame({'Year': series.index, 'Group': g, 'Anomaly': anom.values, 'Source': 'Model'}))
#
#     if include_total:
#         da_total = sum([ds[g] for g in groups])
#         da_mean = da_total.mean(dim=('row', 'col'), skipna=True)
#         series = da_mean.groupby('time.year').mean().to_pandas()
#         series = series[(series.index >= cfg5.START_YEAR) & (series.index <= cfg5.END_YEAR)]
#         clim = series.loc[series.index.isin(range(cfg5.START_YEAR, cfg5.END_YEAR))]
#         anom = (series - clim.mean()) / clim.std()
#         results_mod.append(pd.DataFrame({'Year': series.index, 'Group': 'TOTAL', 'Anomaly': anom.values, 'Source': 'Model'}))
#
#     anom_df = pd.concat(results_obs + results_mod, ignore_index=True)
#
#     all_groups = groups + (['TOTAL'] if include_total else [])
#     fig, axs = plt.subplots(len(all_groups), 2, figsize=(14, 3.5 * len(all_groups)), sharex=True)
#     if len(all_groups) == 1:
#         axs = [axs]
#
#     letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)']
#     for i, g in enumerate(all_groups):
#         for j, source in enumerate(["Model", "Observed"]):
#             sub = anom_df[(anom_df["Group"] == g) & (anom_df["Source"] == source)].set_index("Year")
#             npgo = npgo_ann[sub.index]
#             colors = ['blue' if x > 0 else 'red' for x in sub['Anomaly']]
#             axs[i][j].bar(sub.index, sub['Anomaly'], color=colors)
#             axs[i][j].plot(npgo.index, npgo.values, color='black', label='NPGO (lagged)')
#             axs[i][j].set_xlim(left=cfg5.START_YEAR)
#             axs[i][j].axhline(0, color='grey', linestyle='--')
#             axs[i][j].set_ylabel('Standardised Anomaly')
#             axs[i][j].legend(loc='lower left')
#             axs[i][j].text(0.05, 0.95, f'{letters[i]} {g} - {source}', transform=axs[i][j].transAxes, fontsize=12, va='top')
#
#     out_panel = os.path.join(cfg5.FIGSOUT_P, f"ecospace_{cfg5.ecospace_code}_anom_panel_Zoop_vs_Model_with_NPGO.png")
#     plt.tight_layout(); plt.savefig(out_panel); plt.show()
#     print(f"[5d] Saved panel fig: {out_panel}")
#
#     out_total = None
#     if include_total:
#         tot_obs = anom_df[(anom_df['Group'] == 'TOTAL') & (anom_df['Source'] == 'Observed')].set_index('Year')
#         tot_mod = anom_df[(anom_df['Group'] == 'TOTAL') & (anom_df['Source'] == 'Model')].set_index('Year')
#         npgo = npgo_ann[tot_obs.index]
#         fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
#         axs[0].bar(tot_obs.index, tot_obs['Anomaly'], color=['blue' if x>0 else 'red' for x in tot_obs['Anomaly']])
#         axs[0].plot(npgo.index, npgo.values, color='black', label='NPGO (lagged)')
#         axs[0].axhline(0, color='grey', linestyle='--'); axs[0].set_ylabel('Observed Anomaly'); axs[0].set_title('(a) Total Zooplankton - Observed'); axs[0].legend(loc='lower left')
#         axs[1].bar(tot_mod.index, tot_mod['Anomaly'], color=['blue' if x>0 else 'red' for x in tot_mod['Anomaly']])
#         axs[1].plot(npgo.index, npgo.values, color='black', label='NPGO (lagged)')
#         axs[1].axhline(0, color='grey', linestyle='--'); axs[1].set_ylabel('Modelled Anomaly'); axs[1].set_title('(b) Total Zooplankton - Model'); axs[1].legend(loc='lower left'); axs[1].set_xlabel('Year')
#         axs[1].set_xlim(left=cfg5.START_YEAR)
#         plt.tight_layout()
#         out_total = os.path.join(cfg5.FIGSOUT_P, f"ecospace_{cfg5.ecospace_code}_anom_TOTAL_stacked_Model_vs_Obs.png")
#         plt.savefig(out_total); plt.show()
#         print(f"[5d] Saved total fig: {out_total}")
#
#     return out_panel, out_total
#

# ============================================================
# Orchestrator — callable from other scripts (like 4a)
# ============================================================

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