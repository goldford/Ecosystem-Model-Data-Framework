"""
Script:   Match ecosim data to zoop obs and run seasonal barplots
By: G Oldford, 2025

Purpose:
   - Match Ecosim outputs to Zoop data (wide-form obs with group columns)
   - Compute per-timestep averages of obs and model, then seasonal means
   - Generate seasonal barplots comparing model vs obs for key zooplankton groups
Inputs:
   - Wide-form CSV of zooplankton observations (date + one column per group)
   - Ecosim single-run CSV (with date and season columns)
Outputs:
   - CSV of matched wide table
   - CSV of matched long table
   - PNG of seasonal barplots (6 panels: ZC1–ZC5 + Total)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ecosim_eval_config as cfg
import matplotlib
matplotlib.use('TkAgg')

# --- Configuration ---
SCENARIO = cfg.SCENARIO
ECOSIM_CSV = cfg.ECOSIM_F_PREPPED_SINGLERUN
OBS_CSV = os.path.join(cfg.Z_P_PREPPED, cfg.Z_F_PREPPED)
OUTPUT_DIR = cfg.OUTPUT_DIR_EVAL

# Full list of zoop groups in config (wide-form obs columns)
GROUP_MAP = cfg.Z_GROUP_MAP
# Select only the five C-type groups for plotting, plus a total
ZC_GROUPS = [g for g in GROUP_MAP.keys() if g.startswith('ZC')]
PLOT_GROUPS = ZC_GROUPS + ['Total']
# Seasons in order
SEASONS = ['Winter', 'Spring', 'Summer', 'Fall']
NCOLS = 3
# Tolerance for matching dates (in days)
TIME_TOL = pd.Timedelta(days=cfg.TIMESTEP_DAYS)

# Short labels for panels (strip suffixes)
LABEL_MAP = {grp: grp.split('-')[0] for grp in ZC_GROUPS}
LABEL_MAP['Total'] = 'Total'

# Anomaly analysis window (years) - define in config
default_years = getattr(cfg, 'Z_ANOM_YR_RG', (1980, 2018))
ANOM_START, ANOM_END = default_years

# Labels
title_map = {grp: grp for grp in ZC_GROUPS}
title_map['Total'] = 'Total'


def run_zoop_eval():
    # Load obs and model
    obs = pd.read_csv(OBS_CSV)
    obs['date'] = pd.to_datetime(obs['date']) if 'date' in obs.columns else pd.to_datetime(obs[['Year','Month','Day']])

    mod = pd.read_csv(ECOSIM_CSV, parse_dates=['date'])
    mod['date'] = pd.to_datetime(mod['date'])

    # Prepare model columns
    if 'season' in mod.columns:
        mod_df = mod[['date','season'] + list(map(str,GROUP_MAP.values()))]
    else:
        mod_df = mod.rename(columns={'Season':'season'})[['date','season'] + list(map(str,GROUP_MAP.values()))]
    mod_df = mod_df.rename(columns={str(v): f"EWE-{k}" for k,v in GROUP_MAP.items()})

    # Match by nearest date
    obs = obs.sort_values('date')
    mod_df = mod_df.sort_values('date')
    matched = pd.merge_asof(obs, mod_df, on='date', direction='nearest', tolerance=TIME_TOL)
    matched = matched.dropna(subset=[f"EWE-{g}" for g in GROUP_MAP])

    # Save matched tables
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    wide_csv = os.path.join(OUTPUT_DIR, f"ecosim_{SCENARIO}_zoop_matched_wide.csv")
    matched.to_csv(wide_csv, index=False)
    print(f"Saved wide-match: {wide_csv}")

    # Compute totals
    matched['Total'] = matched[ZC_GROUPS].sum(axis=1)
    matched['EWE-Total'] = matched[[f"EWE-{g}" for g in ZC_GROUPS]].sum(axis=1)

    # Build long format
    long_dfs = []
    for g in PLOT_GROUPS:
        long_dfs.append(
            pd.DataFrame({
                'date': matched['date'],
                'season': matched['season'],
                'group': g,
                'obs_biomass': matched[g],
                'model_biomass': matched[f"EWE-{g}"]
            })
        )
    paired = pd.concat(long_dfs, ignore_index=True)
    long_csv = os.path.join(OUTPUT_DIR, f"ecosim_{SCENARIO}_zoop_matched_long.csv")
    paired.to_csv(long_csv, index=False)
    print(f"Saved long-match: {long_csv}")

    # Seasonal barplots
    seasonal = (
        paired.groupby(['season','group'], as_index=False)
              .agg(obs_biomass=('obs_biomass','mean'), model_biomass=('model_biomass','mean'))
    )
    fig, axes = plt.subplots(int(np.ceil(len(PLOT_GROUPS)/NCOLS)), NCOLS, sharey=True, figsize=(5*NCOLS,4*np.ceil(len(PLOT_GROUPS)/NCOLS)))
    axes = axes.flatten()
    x = np.arange(len(SEASONS)); w=0.35
    for i,(ax,g) in enumerate(zip(axes, PLOT_GROUPS)):
        d = seasonal[seasonal.group==g]
        obs_vals = [d.loc[d.season==s,'obs_biomass'].values[0] if s in d.season.values else np.nan for s in SEASONS]
        mod_vals = [d.loc[d.season==s,'model_biomass'].values[0] if s in d.season.values else np.nan for s in SEASONS]
        ax.bar(x-w/2, obs_vals, w, label='Obs'); ax.bar(x+w/2, mod_vals, w, label='Model')
        ax.set_xticks(x); ax.set_xticklabels(SEASONS); ax.set_title(title_map[g])
        if i%NCOLS==0: ax.set_ylabel('Biomass')
        if g==PLOT_GROUPS[0]: ax.legend()
    for ax in axes[len(PLOT_GROUPS):]: ax.axis('off')
    fig.tight_layout()
    out1 = os.path.join(OUTPUT_DIR, f"ecosim_{SCENARIO}_zoop_seasonal_barplots.png")
    fig.savefig(out1, dpi=300)
    print(f"Saved seasonal plots: {out1}")

    # Compute seasonal anomalies per year
    paired['year'] = paired['date'].dt.year
    hist = paired[(paired.year >= ANOM_START) & (paired.year <= ANOM_END)]
    anoms = []
    for g in PLOT_GROUPS:
        df = hist[hist.group == g].copy()
        clim = df.groupby('season').agg(obs_mean=('obs_biomass', 'mean'), obs_sd=('obs_biomass', 'std'),
                                        mod_mean=('model_biomass', 'mean'),
                                        mod_sd=('model_biomass', 'std')).reset_index()
        df = df.merge(clim, on='season')
        df['obs_anom'] = (df.obs_biomass - df.obs_mean) / df.obs_sd
        df['mod_anom'] = (df.model_biomass - df.mod_mean) / df.mod_sd
        anoms.append(df[['year', 'season', 'group', 'obs_anom', 'mod_anom']])
    anom_df = pd.concat(anoms, ignore_index=True)
    anom_csv = os.path.join(OUTPUT_DIR, f"ecosim_{SCENARIO}_zoop_anomalies.csv")
    anom_df.to_csv(anom_csv, index=False)
    print(f"Saved anomalies: {anom_csv}")



    # Plot anomaly barplots by year
    fig2, axes2 = plt.subplots(int(np.ceil(len(PLOT_GROUPS) / NCOLS)), NCOLS, sharey=True,
                               figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS) / NCOLS)))
    axes2 = axes2.flatten()
    for i, (ax, g) in enumerate(zip(axes2, PLOT_GROUPS)):
        # Select only numeric anomaly columns for yearly mean
        sub = anom_df[anom_df.group == g][['year', 'obs_anom', 'mod_anom']]
        df = sub.groupby('year', as_index=False).mean()
        ax.bar(df.year - 0.2, df.obs_anom, 0.4, label='Obs')
        ax.bar(df.year + 0.2, df.mod_anom, 0.4, label='Model')
        ax.set_title(title_map[g])
        ax.set_xticks(df.year)
        ax.set_xticklabels(df.year, rotation=90)
        if i % NCOLS == 0:
            ax.set_ylabel('Anom')
        if g == PLOT_GROUPS[0]:
            ax.legend()
    # Turn off any unused axes
    for ax in axes2[len(PLOT_GROUPS):]:
        ax.axis('off')
    fig2.tight_layout()
    out2 = os.path.join(OUTPUT_DIR, f"ecosim_{SCENARIO}_zoop_anomaly_barplots.png")
    fig2.savefig(out2, dpi=300)
    print(f"Saved anomaly plots: {out2}")
    fig2.tight_layout()
    fig2.show()
OUTPUT_DIR SHOULD BE FIGS
    out2 = os.path.join(OUTPUT_DIR, f"ecosim_{SCENARIO}_zoop_anomaly_barplots.png")
    fig2.savefig(out2, dpi=300)
    print(f"Saved anomaly plots: {out2}")


if __name__=='__main__':
    run_zoop_eval()