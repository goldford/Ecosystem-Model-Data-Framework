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
import ecosim_eval_config

# ====== Config ======
SCENARIO = ecosim_eval_config.SCENARIO
# 1D Ecosim (wide-form) output CSV with numeric columns for zooplankton groups
ecosim_csv = ecosim_eval_config.ECOSIM_F_PREPPED_SINGLERUN
# Semi-prepared zooplankton observations CSV
oobs_csv = ecosim_eval_config.ZOOP_F_PREPPED  # define this in config
# Where to write matched table
OUTPUT_DIR_EVAL = ecosim_eval_config.OUTPUT_DIR_EVAL

# Mapping from obs column names to numeric column indices in the model output
GROUP_MAP = ecosim_eval_config.Z_GROUP_MAP

def run_zoop_matching():
    """
    Load obs and model 1D outputs, pivot to long form, match on date + group,
    and save paired table for downstream analysis.
    """
    # ---- Load observations ----
    obs = pd.read_csv(oobs_csv)
    # If date column is not already present, construct it
    if 'date' not in obs.columns:
        obs['date'] = pd.to_datetime(obs[['Year', 'Month', 'Day']])
    else:
        obs['date'] = pd.to_datetime(obs['date'])

    # Pivot obs to long form
    obs_long = (
        obs
        .melt(
            id_vars=['date'],
            value_vars=list(GROUP_MAP.keys()),
            var_name='group',
            value_name='obs_biomass'
        )
        .dropna(subset=['obs_biomass'])
    )

    # ---- Load model output ----
    mod = pd.read_csv(ecosim_csv, parse_dates=['date'])

    # Build long-form from numeric columns via mapping
    mod_list = []
    for grp, col_idx in GROUP_MAP.items():
        col_name = str(col_idx)
        if col_name not in mod.columns:
            raise KeyError(f"Column '{col_name}' not found in model output.")
        df = mod[['date', col_name]].rename(columns={col_name: 'model_biomass'})
        df['group'] = grp
        mod_list.append(df)
    mod_long = pd.concat(mod_list, ignore_index=True)

    # ---- Match observations to model ----
    paired = pd.merge(
        obs_long,
        mod_long,
        on=['date', 'group'],
        how='left'
    )

    # ---- Save results ----
    os.makedirs(OUTPUT_DIR_EVAL, exist_ok=True)
    out_file = os.path.join(
        OUTPUT_DIR_EVAL,
        f"ecosim_{SCENARIO}_zoop_obs_model_matched.csv"
    )
    paired.to_csv(out_file, index=False)
    print(f"Matched zooplankton table saved to: {out_file}")


if __name__ == '__main__':
    run_zoop_matching()
