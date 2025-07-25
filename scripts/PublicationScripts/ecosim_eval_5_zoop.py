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
OBS_CSV_SEAS = os.path.join(cfg.Z_P_PREPPED, cfg.Z_F_SEAS)
OBS_CSV_TOWLEV = os.path.join(cfg.Z_P_PREPPED, cfg.Z_F_TOWLEV)
OUTPUT_DIR_STATS = cfg.OUTPUT_DIR_EVAL
OUTPUT_DIR_FIGS = cfg.OUTPUT_DIR_FIGS


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

# Labels
title_map = {grp: grp for grp in ZC_GROUPS}
title_map['Total'] = 'Total'

# USER SETTING: choose one season and plot type ('bar' or 'line')
SEASON_CHOICE = cfg.ZP_SEASON_CHOICE  # e.g., 'Winter','Spring','Summer','Fall'
PLOT_TYPE = cfg.ZP_PLOT_TYPE  # 'bar' or 'line'

# USER SETTING: toggle friendly labels
USE_FRIENDLY_LABELS = cfg.ZP_USE_FRIENDLY_LABELS  # set False to use cryptic codes
FRIENDLY_MAP = cfg.ZP_FRIENDLY_MAP_ZC
# USER SETTING: specify anomaly year range
ANOM_YEAR_START = cfg.ZP_ANOM_YEAR_START
ANOM_YEAR_END = cfg.ZP_ANOM_YEAR_END


def compute_stats(df, obs_col, mod_col):
    """
    Compute a suite of performance stats, handling NaNs.
    """
    # Extract and mask NaNs
    o_series = df[obs_col]
    m_series = df[mod_col]
    mask = o_series.notna() & m_series.notna()
    o = o_series[mask].values
    m = m_series[mask].values
    N = len(o)
    if N == 0:
        # No data to compute
        return dict(MB=np.nan, MAE=np.nan, RMSE=np.nan, NRMSE=np.nan,
                    r=np.nan, R2=np.nan, MAPE=np.nan, NSE=np.nan, WSS=np.nan)
    # Core metrics
    mb   = np.nanmean(m - o)
    mae  = np.nanmean(np.abs(m - o))
    rmse = np.sqrt(np.nanmean((m - o)**2))
    nrmse= rmse / np.nanmean(o) if np.nanmean(o) != 0 else np.nan
    # Correlation
    r    = np.corrcoef(m, o)[0,1] if N > 1 else np.nan
    r2   = r**2 if not np.isnan(r) else np.nan
    # Percent error
    mape = np.nanmean(np.abs((m - o) / o)) * 100 if np.all(o != 0) else np.nan
    # Efficiency
    denom_nse = np.nansum((o - np.nanmean(o))**2)
    nse = 1 - np.nansum((m - o)**2) / denom_nse if denom_nse != 0 else np.nan
    # Willmott's skill score
    denom_wss = np.nansum((np.abs(m - np.nanmean(o)) + np.abs(o - np.nanmean(o)))**2)
    wss = 1 - np.nansum((m - o)**2) / denom_wss if denom_wss != 0 else np.nan
    return {
        'N': N,
        'MB': mb, 'MAE': mae, 'RMSE': rmse, 'NRMSE': nrmse,
        'r': r, 'R2': r2, 'MAPE': mape, 'NSE': nse, 'WSS': wss
    }



def run_zoop_eval():

    # # -----------------------------------------
    # # data prep
    # # -----------------------------------------
    #
    # obs = pd.read_csv(OBS_CSV_SEAS)
    # obs['date'] = pd.to_datetime(obs['date']) if 'date' in obs.columns else pd.to_datetime(obs[['Year','Month','Day']])
    #
    # mod = pd.read_csv(ECOSIM_CSV, parse_dates=['date'])
    # mod['date'] = pd.to_datetime(mod['date'])
    #
    # # Prepare model columns
    # if 'season' in mod.columns:
    #     mod_df = mod[['date','season'] + list(map(str,GROUP_MAP.values()))]
    # else:
    #     mod_df = mod.rename(columns={'Season':'season'})[['date','season'] + list(map(str,GROUP_MAP.values()))]
    # mod_df = mod_df.rename(columns={str(v): f"EWE-{k}" for k,v in GROUP_MAP.items()})
    #
    # # -----------------------------------------
    # # Matching ecosim to zoop obs
    # # -----------------------------------------
    #
    # # Match by nearest date
    # obs = obs.sort_values('date')
    # mod_df = mod_df.sort_values('date')
    # matched = pd.merge_asof(obs, mod_df, on='date', direction='nearest', tolerance=TIME_TOL)
    # matched = matched.dropna(subset=[f"EWE-{g}" for g in GROUP_MAP])
    #
    # # Save matched tables
    # os.makedirs(OUTPUT_DIR_STATS, exist_ok=True)
    # wide_csv = os.path.join(OUTPUT_DIR_STATS, f"ecosim_{SCENARIO}_zoop_matched_wide.csv")
    # matched.to_csv(wide_csv, index=False)
    # print(f"Saved wide-match: {wide_csv}")
    #
    # # Compute totals
    # matched['Total'] = matched[ZC_GROUPS].sum(axis=1)
    # matched['EWE-Total'] = matched[[f"EWE-{g}" for g in ZC_GROUPS]].sum(axis=1)
    #
    # # Build long format
    # long_dfs = []
    # for g in PLOT_GROUPS:
    #     long_dfs.append(
    #         pd.DataFrame({
    #             'date': matched['date'],
    #             'season': matched['season'],
    #             'group': g,
    #             'obs_biomass': matched[g],
    #             'model_biomass': matched[f"EWE-{g}"]
    #         })
    #     )
    # paired = pd.concat(long_dfs, ignore_index=True)
    # long_csv = os.path.join(OUTPUT_DIR_STATS, f"ecosim_{SCENARIO}_zoop_matched_long.csv")
    # paired.to_csv(long_csv, index=False)
    # print(f"Saved long-match: {long_csv}")
    #
    #
    # # -----------------------------------------
    # # stratify by year→season before climatolog
    # # -----------------------------------------
    #
    # # add calendar year
    # paired['year'] = paired['date'].dt.year
    #
    # # 1) Compute mean biomass per (year, season, group)
    # year_seas = (
    #     paired
    #     .groupby(['year', 'season', 'group'], as_index=False)
    #     .agg(
    #         obs_mean=('obs_biomass', 'mean'),
    #         model_mean=('model_biomass', 'mean'),
    #         n_tows=('obs_biomass', 'count')
    #     )
    # )
    #
    # # 2) Compute climatological seasonal mean across years
    # climatology = (
    #     year_seas
    #     .groupby(['season', 'group'], as_index=False)
    #     .agg(
    #         obs_clim=('obs_mean', 'mean'),
    #         model_clim=('model_mean', 'mean')
    #     )
    # )
    #
    # # 3) Prepare panel settings
    # SEASONS = ['Winter', 'Spring', 'Summer', 'Fall']
    # PLOT_GROUPS2 = sorted(climatology['group'].unique())  # ensure a consistent order
    # NCOLS = 2
    #
    # # 4) Generate seasonal barplots in panels
    # fig, axes = plt.subplots(
    #     int(np.ceil(len(PLOT_GROUPS2) / NCOLS)),
    #     NCOLS,
    #     sharey=True,
    #     figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS2) / NCOLS))
    # )
    #
    # axes = axes.flatten()
    # x = np.arange(len(SEASONS))
    # w = 0.35
    #
    # for i, (ax, g) in enumerate(zip(axes, PLOT_GROUPS2)):
    #     d = climatology[climatology['group'] == g]
    #     obs_vals = [
    #         d.loc[d['season'] == s, 'obs_clim'].values[0] if s in d['season'].values else np.nan
    #         for s in SEASONS
    #     ]
    #     mod_vals = [
    #         d.loc[d['season'] == s, 'model_clim'].values[0] if s in d['season'].values else np.nan
    #         for s in SEASONS
    #     ]
    #     ax.bar(x - w / 2, obs_vals, w, label='Obs')
    #     ax.bar(x + w / 2, mod_vals, w, label='Model')
    #     ax.set_xticks(x)
    #     ax.set_xticklabels(SEASONS)
    #     ax.set_title(g)
    #     if i % NCOLS == 0:
    #         ax.set_ylabel('Biomass (g C m⁻²)')
    #     if i == 0:
    #         ax.legend()
    #
    # # Turn off any unused panels
    # for ax in axes[len(PLOT_GROUPS2):]:
    #     ax.axis('off')
    #
    # fig.tight_layout()
    # out = os.path.join(OUTPUT_DIR_FIGS, f"ecosim_{SCENARIO}_zoop_seasonal_barplots.png")
    # fig.savefig(out, dpi=300)
    # print(f"Saved seasonal plots: {out}")
    # plt.show()



    # === Load Observed Seasonal Summary ===
    obs = pd.read_csv(OBS_CSV_SEAS)
    obs['season'] = obs['season'].str.capitalize()
    obs = obs.rename(columns={'modelgroup': 'group',
                              'adj_C_g_m2': 'obs_biomass'})
    # Collapse species-level entries to one value per season-group
    obs_grouped = (
        obs
        .groupby(['season', 'group'], as_index=False)
        ['obs_biomass']
        .sum()
    )
    # Pivot to wide to compute Total
    obs_wide = obs_grouped.pivot(index='season', columns='group', values='obs_biomass').reset_index()
    obs_wide['Total'] = obs_wide[ZC_GROUPS].sum(axis=1)
    # Melt back to long format
    obs_long = obs_wide.melt(id_vars=['season'],
                             value_vars=PLOT_GROUPS,
                             var_name='group',
                             value_name='obs_biomass')

    # === Load Model Output & Compute Seasonal Means ===
    mod = pd.read_csv(ECOSIM_CSV, parse_dates=['date'])
    mod['season'] = mod['season'].str.capitalize()
    # Rename numeric columns to group codes
    mod = mod.rename(columns={str(v): k for k, v in GROUP_MAP.items()})
    # Compute total biomass across groups
    mod['Total'] = mod[ZC_GROUPS].sum(axis=1)
    # Melt to long format for consistency
    long_mod = pd.melt(mod,
                       id_vars=['season'],
                       value_vars=PLOT_GROUPS,
                       var_name='group',
                       value_name='model_biomass')
    # Compute mean per season×group
    seasonal_model = (
        long_mod
        .groupby(['season', 'group'], as_index=False)
        ['model_biomass']
        .mean()
    )

    # === Merge Observations and Model Summaries ===
    comp = pd.merge(obs_long, seasonal_model,
                    on=['season', 'group'],
                    how='inner')

    # === Export Comparison Table ===
    os.makedirs(OUTPUT_DIR_STATS, exist_ok=True)
    out_csv = os.path.join(OUTPUT_DIR_STATS,
                           f"zoop_obs_vs_model_{SCENARIO}.csv")
    comp.to_csv(out_csv, index=False)
    print(f"Comparison written to: {out_csv}")

    # === Plot Obs vs Model Seasonal Barplots ===
    fig, axes = plt.subplots(int(np.ceil(len(PLOT_GROUPS) / NCOLS)), NCOLS,
                             sharey=True,
                             figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS) / NCOLS)))
    axes = axes.flatten()
    x = np.arange(len(SEASONS));
    w = 0.35
    for i, (grp, ax) in enumerate(zip(PLOT_GROUPS, axes)):
        d = comp[comp.group == grp]
        obs_vals = [d.loc[d.season == s, 'obs_biomass'].values[0] if s in d.season.values else np.nan for s in SEASONS]
        mod_vals = [d.loc[d.season == s, 'model_biomass'].values[0] if s in d.season.values else np.nan for s in
                    SEASONS]
        ax.bar(x - w / 2, obs_vals, w, label='Obs')
        ax.bar(x + w / 2, mod_vals, w, label='Model')
        ax.set_xticks(x)
        ax.set_xticklabels(SEASONS)
        ax.set_title(title_map.get(grp, grp))
        if i % NCOLS == 0:
            ax.set_ylabel('Biomass (g/m² DW)')
        if grp == PLOT_GROUPS[0]:
            ax.legend()
    # Disable unused axes
    for ax in axes[len(PLOT_GROUPS):]:
        ax.axis('off')
    fig.tight_layout()
    fig.show()
    os.makedirs(OUTPUT_DIR_FIGS, exist_ok=True)
    fig.savefig(os.path.join(OUTPUT_DIR_FIGS,
                             f"zoop_obs_vs_model_{SCENARIO}.png"), dpi=300)



    # # -------------------------------------------
    # # Compute seasonal anomalies per year
    # # -------------------------------------------
    # paired = comp # quick fix
    # paired['year'] = paired['date'].dt.year
    # hist = paired[(paired.year >= ANOM_START) & (paired.year <= ANOM_END)]
    # anoms = []
    # for g in PLOT_GROUPS:
    #     df = hist[hist.group == g].copy()
    #     clim = df.groupby('season').agg(obs_mean=('obs_biomass', 'mean'), obs_sd=('obs_biomass', 'std'),
    #                                     mod_mean=('model_biomass', 'mean'),
    #                                     mod_sd=('model_biomass', 'std')).reset_index()
    #     df = df.merge(clim, on='season')
    #     df['obs_anom'] = (df.obs_biomass - df.obs_mean) / df.obs_sd
    #     df['mod_anom'] = (df.model_biomass - df.mod_mean) / df.mod_sd
    #     anoms.append(df[['year', 'season', 'group', 'obs_anom', 'mod_anom']])
    # anom_df = pd.concat(anoms, ignore_index=True)
    # anom_csv = os.path.join(OUTPUT_DIR_STATS, f"ecosim_{SCENARIO}_zoop_anomalies.csv")
    # anom_df.to_csv(anom_csv, index=False)
    # print(f"Saved anomalies: {anom_csv}")
    #
    # # -------------------------------------------
    # # Plot anomaly barplots by year
    # # -------------------------------------------
    #
    # fig2, axes2 = plt.subplots(int(np.ceil(len(PLOT_GROUPS) / NCOLS)), NCOLS, sharey=True,
    #                            figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS) / NCOLS)))
    # axes2 = axes2.flatten()
    # for i, (ax, g) in enumerate(zip(axes2, PLOT_GROUPS)):
    #     # Select only numeric anomaly columns for yearly mean
    #     sub = anom_df[anom_df.group == g][['year', 'obs_anom', 'mod_anom']]
    #     df = sub.groupby('year', as_index=False).mean()
    #     ax.bar(df.year - 0.2, df.obs_anom, 0.4, label='Obs')
    #     ax.bar(df.year + 0.2, df.mod_anom, 0.4, label='Model')
    #     ax.set_title(title_map[g])
    #     ax.set_xticks(df.year)
    #     ax.set_xticklabels(df.year, rotation=90)
    #     if i % NCOLS == 0:
    #         ax.set_ylabel('Anom')
    #     if g == PLOT_GROUPS[0]:
    #         ax.legend()
    # # Turn off any unused axes
    # for ax in axes2[len(PLOT_GROUPS):]:
    #     ax.axis('off')
    #
    # fig2.tight_layout()
    # out2 = os.path.join(OUTPUT_DIR_FIGS, f"ecosim_{SCENARIO}_zoop_anomaly_barplots.png")
    # fig2.savefig(out2, dpi=300)
    # print(f"Saved anomaly plots: {out2}")
    # plt.show()

    # === Anomaly Section (single-season, toggle plot type) ===

    # Prepare tow-level obs for anomalies
    obs_raw = pd.read_csv(os.path.join(cfg.Z_P_PREPPED, OBS_CSV_TOWLEV))
    obs_raw['tow_depth_range'] = obs_raw['Tow_start_depth.m.'].abs() - obs_raw['Tow_end_depth.m.'].abs()
    obs_raw['tow_prop'] = obs_raw['tow_depth_range'] / obs_raw['Bottom_depth.m.']
    obs_raw = obs_raw[(obs_raw['tow_prop'] >= 0.7) | (obs_raw['Tow_start_depth.m.'] >= 150)].copy()
    if 'date' in obs_raw:
        obs_raw['date'] = pd.to_datetime(obs_raw['date'])
    else:
        obs_raw['date'] = pd.to_datetime(obs_raw[['Year', 'Month', 'Day']])
    if 'season' in obs_raw:
        obs_raw['season'] = obs_raw['season'].str.capitalize()
    else:
        obs_raw['season'] = obs_raw['Month'].map(
            {12: 'Winter', 1: 'Winter', 2: 'Winter', 3: 'Spring', 4: 'Spring', 5: 'Spring', 6: 'Summer', 7: 'Summer',
             8: 'Summer', 9: 'Fall', 10: 'Fall', 11: 'Fall'})

    # Load and prepare model for matching
    mod_df = pd.read_csv(ECOSIM_CSV, parse_dates=['date'])
    mod_df['season'] = mod_df['season'].str.capitalize()
    for k, v in GROUP_MAP.items():
        mod_df = mod_df.rename(columns={str(v): f"EWE-{k}"})
    mod_cols = ['date', 'season'] + [f"EWE-{g}" for g in ZC_GROUPS]
    mod_df = mod_df[mod_cols]

    # Merge_asof by date
    obs_raw = obs_raw.sort_values('date')
    mod_df = mod_df.sort_values('date')
    matched = pd.merge_asof(obs_raw, mod_df, on='date', direction='nearest', tolerance=TIME_TOL)
    matched = matched.dropna(subset=[f"EWE-{g}" for g in ZC_GROUPS])
    matched['season'] = matched['Season'] if 'Season' in matched else matched['season']

    # Build paired long-format
    records = []
    for g in ZC_GROUPS + ['Total']:
        if g == 'Total':
            matched['Total'] = matched[ZC_GROUPS].sum(axis=1)
            obs_c, mod_c = 'Total', None
        else:
            obs_c, mod_c = g, f"EWE-{g}"
        o = matched[['date', 'season', 'Index', obs_c]].rename(columns={obs_c: 'obs_biomass'})
        if mod_c:
            m = matched[['date', 'season', 'Index', mod_c]].rename(columns={mod_c: 'model_biomass'})
        else:
            matched['model_biomass'] = matched[[f"EWE-{x}" for x in ZC_GROUPS]].sum(axis=1)
            m = matched[['date', 'season', 'Index', 'model_biomass']]
        df_pair = o.merge(m, on=['date', 'season', 'Index'])
        df_pair['group'] = g
        records.append(df_pair)
    paired = pd.concat(records, ignore_index=True)
    paired['year'] = paired['date'].dt.year

    # Filter to selected season
    paired_season = paired[paired['season'] == SEASON_CHOICE]



    # Compute anomalies
    # Use the selected year range instead of default config
    start, end = ANOM_YEAR_START, ANOM_YEAR_END
    hist = paired_season[(paired_season['year'] >= start) & (paired_season['year'] <= end)]
    clim = hist.groupby(['group'], as_index=False).agg(
        mean_obs=('obs_biomass', 'mean'), std_obs=('obs_biomass', 'std'),
        mean_mod=('model_biomass', 'mean'), std_mod=('model_biomass', 'std')
    )
    paired_season = paired_season.merge(clim, on='group', how='left')
    paired_season['obs_anom'] = (paired_season['obs_biomass'] - paired_season['mean_obs']) / paired_season['std_obs']
    paired_season['mod_anom'] = (paired_season['model_biomass'] - paired_season['mean_mod']) / paired_season['std_mod']

    # Aggregate annual anomaly for that season, restricting to selected range
    paired_filtered = paired_season[
        (paired_season['year'] >= ANOM_YEAR_START) & (paired_season['year'] <= ANOM_YEAR_END)]
    annual = paired_filtered.groupby(['year', 'group'], as_index=False).agg(
        obs_anom=('obs_anom', 'mean'),
        model_anom=('mod_anom', 'mean')
    )
    years = sorted(annual['year'].unique())

    # Plot for each group
    grid = plt.subplots(int(np.ceil(len(PLOT_GROUPS) / NCOLS)), NCOLS, sharey=True,
                        figsize=(5 * NCOLS, 4 * np.ceil(len(PLOT_GROUPS) / NCOLS)))
    axes = grid[1].flatten()
    x = np.arange(len(years));
    w = 0.35
    for i, (g, ax) in enumerate(zip(PLOT_GROUPS, axes)):
        d = annual[annual['group'] == g]
        y_obs = [d.loc[d.year == yr, 'obs_anom'].values[0] if yr in d.year.values else np.nan for yr in years]
        y_mod = [d.loc[d.year == yr, 'model_anom'].values[0] if yr in d.year.values else np.nan for yr in years]
        if PLOT_TYPE == 'bar':
            ax.bar(x - w / 2, y_obs, w, label='Obs')
            ax.bar(x + w / 2, y_mod, w, label='Model')
        else:
            ax.plot(years, y_obs, '-o', label='Obs')
            ax.plot(years, y_mod, '-o', label='Model')
        ax.set_xticks(x);
        ax.set_xticklabels(years, rotation=45)
        # inside your plot loop, replace the ax_title logic with:
        if USE_FRIENDLY_LABELS:
            # e.g. “Euphausiids (ZC1‑EUP)”
            ax_title = f"{FRIENDLY_MAP.get(g, g)} ({g})"
        else:
            ax_title = g
        ax.set_title(ax_title)
        ax.set_title(ax_title)
        if i % NCOLS == 0: ax.set_ylabel('Anomaly (z-score)')
        if g == PLOT_GROUPS[0]: ax.legend()

    for ax in axes[len(PLOT_GROUPS):]: ax.axis('off')
    plt.tight_layout()
    plt.show()

    plt.savefig(os.path.join(OUTPUT_DIR_FIGS, f"zoop_anomalies_{SCENARIO}.png"), dpi=300)


    # EVAL stats
    stats_list = []
    for grp in PLOT_GROUPS:
        df_grp = paired_filtered[paired_filtered['group'] == grp]
        if len(df_grp) > 1:
            stats = compute_stats(df_grp, 'obs_anom', 'mod_anom')
            stats['group'] = grp
            stats_list.append(stats)

    stats_df = pd.DataFrame(stats_list)
    print(f"Performance metrics by group (seasonal anomalies - {SEASON_CHOICE}):")
    print(stats_df.to_string(index=False))

    # Write metrics to CSV
    out_stats_csv = os.path.join(OUTPUT_DIR_STATS, f"zoop_performance_{SCENARIO}_{SEASON_CHOICE}.csv")
    stats_df.to_csv(out_stats_csv, index=False)
    print(f"Metrics written to: {out_stats_csv}")

if __name__=='__main__':
    run_zoop_eval()