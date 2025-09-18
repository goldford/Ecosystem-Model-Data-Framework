"""
Ecospace Zooplankton Evaluation — Combined 5b/5c/5d (aligned to 4a pattern)
by: G. Oldford + ChatGPT helper, 2025-08-11

What this does
--------------
1) Match zooplankton observations to Ecospace outputs (former 5b)
2) Visualize & compute stats (former 5c)
3) Make anomaly comparisons incl. NPGO overlays (former 5d)

How to use
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
- Keeps logic and naming close to your 4a script for consistency
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
    BASEMAP_P: str = cfg.ECOSPACE_MAP_P

    # Zoop & mapping input files (set these in cfg or keep defaults here)
    ZOOP_P: str = getattr(cfg, "Z_P_PREPPED", "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/4. Zooplankton/Zoop_Perryetal_2021/MODIFIED")
    ZOOP_CSV: str = getattr(cfg, "Z_F_TOWLEV", "Zooplankton_B_C_gm2_EWEMODELGRP_Wide.csv")
    GRID_MAP_CSV: str = getattr(cfg, "ECOSPACE_GRID_RC_CSV", "Ecospace_grid_20210208_rowscols.csv")

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


# ============================================================
# Stage 1 — Match Zooplankton obs to Ecospace (ex-5b)
# ============================================================

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

    out_csv = os.path.join(cfg5.EVALOUT_P, f"Zooplankton_matched_to_model_out_{cfg5.ecospace_code}.csv")
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

    # Log10 for skill stats
    for col in obs_cols + ['TOT']:
        df_station[f'log10_{col}'] = np.log10(df_station[col] + cfg5.LOG_OFFSET)
        mcol = f"EWE-{col}" if col != 'TOT' else 'EWE-TOT'
        df_station[f'log10_{mcol}'] = np.log10(df_station[mcol] + cfg5.LOG_OFFSET)

    # Seasonal mean then annual mean of seasonal means (Perry-style)
    seasonal_means = df_station.groupby('Season')[
        [f'log10_{c}' for c in obs_cols + ['TOT']] + [f'log10_EWE-{c}' for c in obs_cols] + ['log10_EWE-TOT']
    ].mean().reset_index()
    annual_means = seasonal_means.mean(numeric_only=True)

    df_clim_anom = df_station.copy()
    for col in obs_cols + ['TOT']:
        ocol = f'log10_{col}'
        mcol = f'log10_EWE-{col}' if col != 'TOT' else 'log10_EWE-TOT'
        df_clim_anom[f'anom_annmean_{ocol}'] = df_clim_anom[ocol] - annual_means[ocol]
        df_clim_anom[f'anom_annmean_{mcol}'] = df_clim_anom[mcol] - annual_means[mcol]

    # Compute skill
    rows = []
    for col in obs_cols + ['TOT']:
        olog = df_station[f'log10_{col}']
        mcol = f"EWE-{col}" if col != 'TOT' else 'EWE-TOT'
        mlog = df_station[f'log10_{mcol}']
        s = _compute_skill_statistics(olog.values, mlog.values)
        s['Group'] = col
        rows.append(s)
    skill_df = pd.DataFrame(rows).sort_values('WSS', ascending=False)

    out_skill = os.path.join(cfg5.EVALOUT_P, f"model_skill_stats_log10_{cfg5.ecospace_code}.csv")
    skill_df.to_csv(out_skill, index=False)
    print(f"[5c] Saved skill stats: {out_skill}")

    # Optional figures (kept concise)
    if cfg5.MAKE_VIZ:
        # Total scatter by Region
        fig, ax = plt.subplots(figsize=(6,6))
        sns.scatterplot(x=df_station['log10_TOT'], y=df_station['log10_EWE-TOT'], hue=df_station['Region'], alpha=0.7, ax=ax)
        ax.plot([-3,3],[-3,3],'k--',lw=1); ax.set_xlim(-3,3); ax.set_ylim(-3,3)
        ax.set_xlabel("log10(Observed + 1e-6)"); ax.set_ylabel("log10(Modelled + 1e-6)"); ax.set_title("Total Biomass by Station (log10)")
        plt.tight_layout(); plt.show()

        # Per-group seasonal boxplots (log scale)
        for group in obs_cols + ['TOT']:
            mgroup = f"EWE-{group}" if group != 'TOT' else 'EWE-TOT'
            dfp = df_station[['Season', group, mgroup]].copy().melt(id_vars='Season', var_name='Source', value_name='Biomass')
            dfp['Source'] = dfp['Source'].map({group:'Observed', mgroup:'Modelled'})
            dfp['Season'] = pd.Categorical(dfp['Season'], categories=cfg5.SEASON_ORDER, ordered=True)
            plt.figure(figsize=(8,6)); sns.boxplot(data=dfp, x='Season', y='Biomass', hue='Source')
            plt.yscale('log'); plt.ylabel('Log-Biomass (g C m⁻²)'); plt.xlabel('Season'); plt.title(f'Seasonal Biomass — {group}')
            plt.tight_layout(); plt.show()

    return out_skill


# ============================================================
# Stage 3 — Anomaly comparison & NPGO (ex-5d)
# ============================================================

def anomaly_comparisons(cfg5: Eval5Config, groups: Optional[List[str]] = None, include_total: bool = True) -> Tuple[str, Optional[str]]:
    print("[5d] Building anomaly comparisons (Model vs Obs) …")
    groups = groups or [g for g in cfg5.ZOOP_GROUPS if g.startswith('ZC') or g.startswith('ZS')][:5]

    # Load NPGO
    npgo_path = os.path.join(cfg5.EVALOUT_P, cfg5.NPGO_CSV)
    df_npgo = pd.read_csv(npgo_path)
    df_npgo.columns = ["date", "npgo"]
    df_npgo["date"] = pd.to_datetime(df_npgo["date"])
    df_npgo = df_npgo[df_npgo["npgo"] > -99]
    df_npgo["year"] = df_npgo["date"].dt.year
    df_npgo["month"] = df_npgo["date"].dt.month
    spring = df_npgo[df_npgo["month"].isin([3,4,5])]
    npgo_ann = spring.groupby("year")["npgo"].mean().shift(cfg5.LAG_YEARS)

    # Load matched CSV
    matched_csv = os.path.join(cfg5.EVALOUT_P, f"Zooplankton_matched_to_model_out_{cfg5.ecospace_code}.csv")
    df = pd.read_csv(matched_csv)

    # Zero replacement for obs
    for col in groups:
        nonzero = df[col][df[col] > 0]
        if not nonzero.empty:
            mn = nonzero.min()
            df[col] = df[col].apply(lambda x: np.random.uniform(0, 0.5 * mn) if x == 0 else x)

    # Station aggregation
    df_station = df.groupby(['Station', 'Year', 'Season'])[groups + [f"EWE-{g}" for g in groups]].mean().reset_index()
    for col in groups:
        df_station[f'log10_{col}'] = np.log10(df_station[col] + cfg5.LOG_OFFSET)

    # OBS anomalies (annual mean of station means, z‑scored relative to 1980‑2017)
    results_obs = []
    for g in groups:
        ser = df_station.groupby('Year')[f'log10_{g}'].mean().dropna()
        ser = ser[(ser.index >= cfg5.START_YEAR) & (ser.index <= cfg5.END_YEAR)]
        clim = ser.loc[ser.index.isin(range(cfg5.START_YEAR, cfg5.END_YEAR))]
        anom = (ser - clim.mean()) / clim.std()
        results_obs.append(pd.DataFrame({'Year': ser.index, 'Group': g, 'Anomaly': anom.values, 'Source': 'Observed'}))

    if include_total:
        df_station['TOT'] = df_station[groups].sum(axis=1)
        df_station['log10_TOT'] = np.log10(df_station['TOT'] + cfg5.LOG_OFFSET)
        ser = df_station.groupby('Year')['log10_TOT'].mean().dropna()
        ser = ser[(ser.index >= cfg5.START_YEAR) & (ser.index <= cfg5.END_YEAR)]
        clim = ser.loc[ser.index.isin(range(cfg5.START_YEAR, cfg5.END_YEAR))]
        anom = (ser - clim.mean()) / clim.std()
        results_obs.append(pd.DataFrame({'Year': ser.index, 'Group': 'TOTAL', 'Anomaly': anom.values, 'Source': 'Observed'}))

    # MODEL anomalies direct from raw Ecospace
    ds = xr.open_dataset(os.path.join(cfg5.NC_PATH_OUT, cfg5.ecospace_nc_name))
    results_mod = []
    for g in groups:
        da = ds[g]
        da_mean = da.mean(dim=('row', 'col'), skipna=True)
        series = da_mean.groupby('time.year').mean().to_pandas()
        series = series[(series.index >= cfg5.START_YEAR) & (series.index <= cfg5.END_YEAR)]
        clim = series.loc[series.index.isin(range(cfg5.START_YEAR, cfg5.END_YEAR))]
        anom = (series - clim.mean()) / clim.std()
        results_mod.append(pd.DataFrame({'Year': series.index, 'Group': g, 'Anomaly': anom.values, 'Source': 'Model'}))

    if include_total:
        da_total = sum([ds[g] for g in groups])
        da_mean = da_total.mean(dim=('row', 'col'), skipna=True)
        series = da_mean.groupby('time.year').mean().to_pandas()
        series = series[(series.index >= cfg5.START_YEAR) & (series.index <= cfg5.END_YEAR)]
        clim = series.loc[series.index.isin(range(cfg5.START_YEAR, cfg5.END_YEAR))]
        anom = (series - clim.mean()) / clim.std()
        results_mod.append(pd.DataFrame({'Year': series.index, 'Group': 'TOTAL', 'Anomaly': anom.values, 'Source': 'Model'}))

    anom_df = pd.concat(results_obs + results_mod, ignore_index=True)

    all_groups = groups + (['TOTAL'] if include_total else [])
    fig, axs = plt.subplots(len(all_groups), 2, figsize=(14, 3.5 * len(all_groups)), sharex=True)
    if len(all_groups) == 1:
        axs = [axs]

    letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)']
    for i, g in enumerate(all_groups):
        for j, source in enumerate(["Model", "Observed"]):
            sub = anom_df[(anom_df["Group"] == g) & (anom_df["Source"] == source)].set_index("Year")
            npgo = npgo_ann[sub.index]
            colors = ['blue' if x > 0 else 'red' for x in sub['Anomaly']]
            axs[i][j].bar(sub.index, sub['Anomaly'], color=colors)
            axs[i][j].plot(npgo.index, npgo.values, color='black', label='NPGO (lagged)')
            axs[i][j].set_xlim(left=cfg5.START_YEAR)
            axs[i][j].axhline(0, color='grey', linestyle='--')
            axs[i][j].set_ylabel('Standardised Anomaly')
            axs[i][j].legend(loc='lower left')
            axs[i][j].text(0.05, 0.95, f'{letters[i]} {g} - {source}', transform=axs[i][j].transAxes, fontsize=12, va='top')

    out_panel = os.path.join(cfg5.EVALOUT_P, f"{cfg5.ecospace_code}_anom_panel_Zoop_vs_Model_with_NPGO.png")
    plt.tight_layout(); plt.savefig(out_panel); plt.show()
    print(f"[5d] Saved panel fig: {out_panel}")

    out_total = None
    if include_total:
        tot_obs = anom_df[(anom_df['Group'] == 'TOTAL') & (anom_df['Source'] == 'Observed')].set_index('Year')
        tot_mod = anom_df[(anom_df['Group'] == 'TOTAL') & (anom_df['Source'] == 'Model')].set_index('Year')
        npgo = npgo_ann[tot_obs.index]
        fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        axs[0].bar(tot_obs.index, tot_obs['Anomaly'], color=['blue' if x>0 else 'red' for x in tot_obs['Anomaly']])
        axs[0].plot(npgo.index, npgo.values, color='black', label='NPGO (lagged)')
        axs[0].axhline(0, color='grey', linestyle='--'); axs[0].set_ylabel('Observed Anomaly'); axs[0].set_title('(a) Total Zooplankton - Observed'); axs[0].legend(loc='lower left')
        axs[1].bar(tot_mod.index, tot_mod['Anomaly'], color=['blue' if x>0 else 'red' for x in tot_mod['Anomaly']])
        axs[1].plot(npgo.index, npgo.values, color='black', label='NPGO (lagged)')
        axs[1].axhline(0, color='grey', linestyle='--'); axs[1].set_ylabel('Modelled Anomaly'); axs[1].set_title('(b) Total Zooplankton - Model'); axs[1].legend(loc='lower left'); axs[1].set_xlabel('Year')
        axs[1].set_xlim(left=cfg5.START_YEAR)
        plt.tight_layout()
        out_total = os.path.join(cfg5.EVALOUT_P, f"{cfg5.ecospace_code}_anom_TOTAL_stacked_Model_vs_Obs.png")
        plt.savefig(out_total); plt.show()
        print(f"[5d] Saved total fig: {out_total}")

    return out_panel, out_total


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
    run_zoop_eval()
