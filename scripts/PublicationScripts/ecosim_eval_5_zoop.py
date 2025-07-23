"""
Script:   Match ecosim data to zoop obs and run analysis to compare
By: G Oldford, 2025
Purpose:
   - Match Ecosim outputs to Zoop data (previously processed so model groups as cols)
   - Run ecosim analyses of anoms of z groups for calib and tuning
Input:
  - CSV of zooplankton groups processed in previous script
  - Mainly from Mackas and Perry publications (2021)

Output:
 -

"""

import pandas as pd
import os
import ecosim_eval_config as cfg

SCENARIO = cfg.SCENARIO
# 1D model output and obs tables
ECOSIM_CSV = cfg.ECOSIM_F_PREPPED_SINGLERUN
OBS_CSV = os.path.join(cfg.Z_P_PREPPED, cfg.Z_F_PREPPED)
OUTPUT_DIR = cfg.OUTPUT_DIR_EVAL

# mapping from obs groups to model numeric columns
GROUP_MAP = cfg.Z_GROUP_MAP
# tolerance for matching dates (in days)
TIME_TOL = pd.Timedelta(days=cfg.TIMESTEP_DAYS)


def run_zoop_eval():
    # load observations
    obs = pd.read_csv(OBS_CSV)
    # ensure date column
    if 'date' in obs.columns:
        obs['date'] = pd.to_datetime(obs['date'])
    else:
        obs['date'] = pd.to_datetime(obs[['Year', 'Month', 'Day']])

    # load model output
    mod = pd.read_csv(ECOSIM_CSV, parse_dates=['date'])
    mod['date'] = pd.to_datetime(mod['date'])

    # sort for merge_asof
    obs = obs.sort_values('date')
    mod = mod.sort_values('date')

    # prepare model columns: rename numeric cols to EWE-<group>
    mod_ewe = mod[['date'] + [str(v) for v in GROUP_MAP.values()]].copy()
    rename_map = {str(v): f"EWE-{k}" for k, v in GROUP_MAP.items()}
    mod_ewe = mod_ewe.rename(columns=rename_map)

    # merge on nearest date
    matched = pd.merge_asof(
        left=obs,
        right=mod_ewe,
        on='date',
        direction='nearest',
        tolerance=TIME_TOL
    )

    # drop rows without a match
    matched = matched.dropna(subset=[f"EWE-{g}" for g in GROUP_MAP.keys()])

    # save wide matched table
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    wide_out = os.path.join(OUTPUT_DIR, f"ecosim_{SCENARIO}_zoop_matched_wide.csv")
    matched.to_csv(wide_out, index=False)
    print(f"Wide match saved to: {wide_out}")

    # build long paired table
    long_rows = []
    for grp in GROUP_MAP.keys():
        df = pd.DataFrame({
            'date': matched['date'],
            'group': grp,
            'obs_biomass': matched[grp],
            'model_biomass': matched[f"EWE-{grp}"]
        })
        long_rows.append(df)
    paired_long = pd.concat(long_rows, ignore_index=True)

    # save paired long table
    long_out = os.path.join(OUTPUT_DIR, f"ecosim_{SCENARIO}_zoop_obs_model_matched.csv")
    paired_long.to_csv(long_out, index=False)
    print(f"Long-paired table saved to: {long_out}")


if __name__ == '__main__':
    run_zoop_eval()

