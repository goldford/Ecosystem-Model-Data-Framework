from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

STATS_GLOB = "ecospace_bloom_timing_stats_*.csv"
AGREE_GLOB = "ecospace_bloom_timing_agreement_*.csv"


def _scenario_from_stats_name(name: str) -> str:
    """Extract <SCENARIO> from 'ecospace_bloom_timing_stats_<SCENARIO>...'."""
    # Prefer the common '..._stats_<SCENARIO>_bt_' naming
    m = re.match(r"ecospace_bloom_timing_stats_(.+?)_bt_", name)
    if m:
        return m.group(1)
    # Fallback: between '..._stats_' and '.csv'
    m = re.match(r"ecospace_bloom_timing_stats_(.+)\.csv$", name)
    return m.group(1) if m else "unknown"


def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan


def _find_run_dirs(root: Path) -> list[Path]:
    runs = root / "runs"
    if not runs.exists():
        return []
    return [p for p in runs.iterdir() if p.is_dir()]


def _read_stats(stats_path: Path, run_tag: str) -> pd.DataFrame:
    df = pd.read_csv(stats_path)
    df["run_tag"] = run_tag
    df["scenario"] = _scenario_from_stats_name(stats_path.name)
    return df


def _read_agree(agree_path: Path, run_tag: str) -> pd.DataFrame:
    df = pd.read_csv(agree_path)
    df["run_tag"] = run_tag
    # scenario is usually parallel to stats; infer if possible
    # expected: ecospace_bloom_timing_agreement_<SCENARIO>_bt_...
    m = re.match(r"ecospace_bloom_timing_agreement_(.+?)_bt_", agree_path.name)
    df["scenario"] = m.group(1) if m else "unknown"
    return df

def build_leaderboard(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (wide_metrics_df, stats_long_df, agree_long_df)."""
    run_dirs = _find_run_dirs(root)
    stats_long = []
    agree_long = []

    for rd in run_dirs:
        run_tag = rd.name
        # stats
        stats_files = sorted(rd.glob(STATS_GLOB))
        for sf in stats_files:
            try:
                stats_long.append(_read_stats(sf, run_tag))
            except Exception as e:
                print(f"WARN: failed reading stats {sf}: {e}")
        # agreement
        agree_files = sorted(rd.glob(AGREE_GLOB))
        for af in agree_files:
            try:
                agree_long.append(_read_agree(af, run_tag))
            except Exception as e:
                print(f"WARN: failed reading agreement {af}: {e}")

    stats_long_df = pd.concat(stats_long, ignore_index=True) if stats_long else pd.DataFrame()
    agree_long_df = pd.concat(agree_long, ignore_index=True) if agree_long else pd.DataFrame()

    if stats_long_df.empty:
        return pd.DataFrame(), stats_long_df, agree_long_df

    # Pivot stats to wide: one row per run_tag+scenario, with sat_* and C09_* columns
    metric_cols = [c for c in stats_long_df.columns if c not in ("run_tag", "scenario", "Label")]
    wide = stats_long_df.pivot_table(
        index=["run_tag", "scenario"],
        columns="Label",
        values=metric_cols,
        aggfunc="first",
    )
    # Flatten MultiIndex columns: (<metric>, <label>) -> <label>_<metric>
    wide.columns = [f"{lbl}_{met}" for (met, lbl) in wide.columns.to_list()]
    wide = wide.reset_index()

    # Add agreement metrics (optional)
    if not agree_long_df.empty:
        agree_wide = agree_long_df.pivot_table(
            index=["run_tag", "scenario"],
            columns=["Label", "Type"],
            values=["Count", "Total", "Proportion"],
            aggfunc="first",
        )
        agree_wide.columns = [f"{lbl}_{typ}_{val}" for (val, lbl, typ) in agree_wide.columns.to_list()]
        agree_wide = agree_wide.reset_index()
        wide = wide.merge(agree_wide, on=["run_tag", "scenario"], how="left")

    # Attach runlog status if present
    runlog_pf = root / "runlog.csv"
    if runlog_pf.exists():
        try:
            rl = pd.read_csv(runlog_pf)
            # expected columns: run_tag, status, ...
            if "run_tag" in rl.columns:
                wide = wide.merge(rl, on="run_tag", how="left")
        except Exception as e:
            print(f"WARN: failed reading runlog.csv: {e}")

    return wide, stats_long_df, agree_long_df


def add_score(df: pd.DataFrame, score: str) -> tuple[pd.DataFrame, bool]:
    """Add 'score' column. Returns (df, higher_is_better)."""
    d = df.copy()

    def col(name: str) -> pd.Series:
        return pd.to_numeric(d.get(name), errors="coerce")

    score = score.lower().strip()

    # Lower-is-better
    if score in ("rmse_mean", "mean_rmse"):
        d["score"] = np.nanmean(np.vstack([col("sat_RMSE"), col("C09_RMSE")]), axis=0)
        hib = False
    elif score in ("rmse_sum", "sum_rmse"):
        d["score"] = col("sat_RMSE") + col("C09_RMSE")
        hib = False
    elif score == "sat_rmse":
        d["score"] = col("sat_RMSE")
        hib = False
    elif score == "c09_rmse":
        d["score"] = col("C09_RMSE")
        hib = False
    elif score in ("mae_mean", "mean_mae"):
        d["score"] = np.nanmean(np.vstack([col("sat_MAE"), col("C09_MAE")]), axis=0)
        hib = False

    # Higher-is-better
    elif score in ("willmott_mean", "mean_willmott"):
        d["score"] = np.nanmean(np.vstack([col("sat_Willmott Skill"), col("C09_Willmott Skill")]), axis=0)
        hib = True
    elif score in ("r_mean", "mean_r"):
        d["score"] = np.nanmean(np.vstack([col("sat_R"), col("C09_R")]), axis=0)
        hib = True
    elif score == "sat_willmott":
        d["score"] = col("sat_Willmott Skill")
        hib = True
    elif score == "c09_willmott":
        d["score"] = col("C09_Willmott Skill")
        hib = True
    elif score == "sat_agree":
        # Uses agreement file proportion if present
        d["score"] = col("sat_Categorical Agreement_Proportion")
        hib = True
    else:
        raise ValueError(
            "Unknown --score. Try one of: rmse_mean (default), rmse_sum, sat_rmse, c09_rmse, "
            "mae_mean, willmott_mean, r_mean, sat_willmott, c09_willmott, sat_agree"
        )

    return d, hib


def main() -> None:
    ap = argparse.ArgumentParser(description="Create a leaderboard from bt_sweep_4a outputs")
    ap.add_argument("--root", required=True, help="Sweep root folder (contains runs/)")
    ap.add_argument("--score", default="rmse_mean", help="Scoring rule (see script for options)")
    ap.add_argument("--scenario", default=None, help="Optional: filter to a specific scenario")
    ap.add_argument("--top", type=int, default=20, help="How many top runs to write as a separate CSV")
    args = ap.parse_args()

    root = Path(args.root).expanduser().resolve()
    if not root.exists():
        raise SystemExit(f"Root not found: {root}")

    wide, stats_long, agree_long = build_leaderboard(root)
    if wide.empty:
        raise SystemExit(f"No stats files found under: {root / 'runs'}")

    if args.scenario:
        wide = wide[wide["scenario"] == args.scenario].copy()
        if wide.empty:
            raise SystemExit(f"No runs match scenario={args.scenario}")

    wide_scored, hib = add_score(wide, args.score)
    wide_scored = wide_scored.sort_values("score", ascending=not hib).reset_index(drop=True)
    wide_scored.insert(0, "rank", np.arange(1, len(wide_scored) + 1))

    out_all = root / "leaderboard.csv"
    out_top = root / f"leaderboard_top{args.top}.csv"

    wide_scored.to_csv(out_all, index=False)
    wide_scored.head(args.top).to_csv(out_top, index=False)

    print(f"Wrote: {out_all}")
    print(f"Wrote: {out_top}")


if __name__ == "__main__":
    main()
