"""
Visualise Ecosim Monte Carlo (MC) outputs vs. EwE-style time series reference file.

This script generalises the earlier "mortality-only" visualiser to handle multiple
time-series types (e.g., TotalMortality, FishingMortality, BiomassAbs, BiomassRel,
BiomassForcing, Catches, CatchesRel, Landings), using an EwE time-series CSV
structured like:

  Title, <Series1>, <Series2>, ...
  Weight, ...
  Pool code, ...
  Pool code 2, ...
  Type, ...
  1950, ...
  1951, ...
  ...

Where:
- "Pool code" typically references either a group code (e.g., 6 for Harbour seals)
  or a fleet code (for catch/landings series), depending on Type.
- "Pool code 2" is optional and is used by some fleet-by-group series (e.g., gear
  catches where pool code = fleet, pool code 2 = group).

Key features:
- Config-driven scenarios + output locations.
- Auto-parses the reference time-series CSV into long format.
- Maps reference Type -> MC output annual CSV (configurable filename candidates).
- Supports composites (multistanza / summed groups) for both reference and model.
- Relative series (BiomassRel, CatchesRel) are normalised prior to comparison.

Outputs:
- One PNG per selected reference series (or composite series).
- Optional CSV summary of fit stats per series x scenario.

"""

from __future__ import annotations

import os
import re
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# =========================
# CONFIG
# =========================

# ---- Paths (override with env vars if you like) ----
EWE_OUTPUT_DIR = Path(os.environ.get("EWE_OUTPUT_DIR", r"C:\Users\Greig\Documents\EwE output"))
MODEL_RUN_DIR = os.environ.get("EWE_MODEL_RUN_DIR", "GeorgiaStrait2025_DietProjectFork_v103_v6_7_0_18858_64b")

# Reference time-series file (uploaded example: TimeSeries_2025_v59.csv)
TIME_SERIES_FILE = Path(os.environ.get("EWE_TIME_SERIES_FILE", r"C:\Users\Greig\Documents\GitHub\ecospace_seal_diets\EwE\time_series\TimeSeries_2025_v59.csv"))

# ---- Scenarios ----
@dataclass(frozen=True)
class Scenario:
    scen_id: str
    label: str
    mc_subdir: str

SCENARIOS: List[Scenario] = [
    Scenario("mc_EcosimScenCalib6", "mc_EcosimScenCalib6 - NO PP Anom",   "mc_EcosimScenCalib6"),
    Scenario("mc_EcosimScenCalib6", "mc_EcosimScenCalib6 - WITH PP Anom", "mc_EcosimScenCalib6"),
]


# ---- Time-series Types ----
# Types in your uploaded file include:
# ['BiomassAbs','BiomassForcing','BiomassRel','Catches','CatchesRel','FishingMortality','Landings','TotalMortality']
RELATIVE_TYPES = {"BiomassRel", "CatchesRel"}
ABSOLUTE_TYPES = {"BiomassAbs", "BiomassForcing", "Catches", "Landings", "FishingMortality", "TotalMortality"}

# For your "main priority is absolute series", set INCLUDE_TYPES to ABSOLUTE_TYPES.
# You can add RELATIVE_TYPES once the core is working.
INCLUDE_TYPES: Optional[set[str]] = set(ABSOLUTE_TYPES)  # set to None to include all types
EXCLUDE_TYPES: set[str] = set()  # e.g., {"BiomassForcing"} if you don't want to plot those

# Optional: limit to a subset of reference series columns by *name*.
# If empty, all matching INCLUDE_TYPES are plotted.
INCLUDE_SERIES_NAMES: List[str] = []  # e.g. ["WCT_mid_B", "SRKW", "Humpback"]
EXCLUDE_SERIES_NAMES: List[str] = []  # e.g. ["Herring_C_Gear1_Cleary2019"]

# Optional: include series by regex (matched against the column name)
INCLUDE_SERIES_REGEX: List[str] = []  # e.g. [r"^Seal_", r"Hake_"]
EXCLUDE_SERIES_REGEX: List[str] = []

# Optional: include only series that reference particular pool codes (group codes) or pool_code2 values.
# Leave empty to disable.
INCLUDE_POOL_CODES: List[int] = []     # e.g. [1,2,3,5,6,8]
INCLUDE_POOL_CODE2: List[int] = []     # e.g. [44,49] (for fleet-by-group catches)

# Optional: shift the reference years for particular series (e.g. brood_year -> year)
SERIES_YEAR_OFFSETS: Dict[str, int] = {
    # "Chinook_WO_M": +1,
}


# Optional: rename series for plotting / filenames
SERIES_ALIASES: Dict[str, str] = {
    # "WCT_mid_B": "Orca-WCT (mid) biomass",
    # "SRKW": "Orca-Resident (SRKW)",
}

# ---- Composite / multistanza series ----
# Create new series by summing existing reference series columns.
# Example (uncomment and adjust):
# COMPOSITES = {
#     "HAKE_C_total": {"members": ["Hake_C", "Hake_C_old"], "op": "sum"},
# }
COMPOSITES: Dict[str, Dict] = {}

# ---- MC output file mapping ----
@dataclass(frozen=True)
class OutputSpec:
    type_name: str
    filename_candidates: Sequence[str]

# You may need to adjust these to match the filenames produced by your MC export.
OUTPUT_SPECS: Dict[str, OutputSpec] = {
    "TotalMortality":   OutputSpec("TotalMortality",   ["mortality_annual.csv", "total_mortality_annual.csv"]),
    "FishingMortality": OutputSpec("FishingMortality", ["fishing_mortality_annual.csv", "F_annual.csv", "f_annual.csv"]),
    "BiomassAbs":       OutputSpec("BiomassAbs",       ["biomass_annual.csv", "Biomass_annual.csv"]),
    "BiomassRel":       OutputSpec("BiomassRel",       ["biomass_annual.csv", "Biomass_annual.csv"]),
    "BiomassForcing":   OutputSpec("BiomassForcing",   ["biomass_annual.csv", "Biomass_annual.csv"]),
    "Catches":          OutputSpec("Catches",          ["catch_annual.csv", "catches_annual.csv", "Catch_annual.csv"]),
    "CatchesRel":       OutputSpec("CatchesRel",       ["catch_annual.csv", "catches_annual.csv", "Catch_annual.csv"]),
    "Landings":         OutputSpec("Landings",         ["landings_annual.csv", "Landings_annual.csv", "landing_annual.csv"]),
}

# ---- Relative-series normalisation ----
# We normalise BOTH observed and model to a baseline, so they are comparable even if one is absolute.
REL_NORM_ENABLED = True
REL_BASELINE_MODE = "overlap_mean"   # "overlap_mean" or "first_year"
REL_BASELINE_YEARS: Optional[Tuple[int, int]] = None  # e.g. (1980, 2010) to force baseline window

# ---- Plot / outputs ----
TRIAL_PREFIX = "mc_output_trial"
YEAR_COL_CANDIDATES = ["year\\group", "year/group", "year\\pool", "year/pool", "year_group", "year", "Year"]

OUT_DIR = Path(os.environ.get("MC_PLOT_DIR", str(EWE_OUTPUT_DIR / "_mc_timeseries_plots")))
SAVE_PLOTS = True
SHOW_PLOTS = True
DPI = 200

YEAR_MIN: Optional[int] = 1970
YEAR_MAX: Optional[int] = 2019

QUANTILES = [5, 25, 50, 75, 95]

COMPUTE_FIT_STATS = True
FIT_STATS_FILENAME = "fit_stats_timeseries.csv"

LOG_LEVEL = os.environ.get("MC_VIZ_LOG_LEVEL", "INFO").upper()


# =========================
# Utilities
# =========================

def setup_logging() -> None:
    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL, logging.INFO),
        format="%(levelname)s | %(message)s",
    )

def safe_mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def sanitize_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", s).strip("_")

def find_year_column(columns: Sequence[str]) -> Optional[str]:
    cols = list(columns)
    for c in YEAR_COL_CANDIDATES:
        if c in cols:
            return c
    for c in cols:
        if "year" in str(c).lower():
            return c
    return None

def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.nanmean((a - b) ** 2)))

def mae(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.nanmean(np.abs(a - b)))


# =========================
# Reference time-series parsing
# =========================

@dataclass(frozen=True)
class RefSeriesMeta:
    name: str
    type_name: str
    pool_code: int
    pool_code2: Optional[int]
    weight: float
    display_name: str

def parse_ewe_timeseries_wide(path: Path) -> Tuple[Dict[str, RefSeriesMeta], pd.DataFrame]:
    """
    Returns:
      meta: dict series_name -> RefSeriesMeta
      ref_long: DataFrame with columns [series, type, pool_code, pool_code2, weight, year, obs]
    """
    wide = pd.read_csv(path)
    if "Title" not in wide.columns:
        raise ValueError(f"Expected a 'Title' column in {path}")

    series_cols = [c for c in wide.columns if c != "Title"]

    # metadata rows
    meta_rows = wide[wide["Title"].isin(["Weight", "Pool code", "Pool code 2", "Type"])].copy()
    if meta_rows.shape[0] < 4:
        raise ValueError("Missing one or more of required metadata rows: Weight, Pool code, Pool code 2, Type")

    meta_rows = meta_rows.set_index("Title")

    weights = pd.to_numeric(meta_rows.loc["Weight", series_cols], errors="coerce").fillna(0.0)
    pool = pd.to_numeric(meta_rows.loc["Pool code", series_cols], errors="coerce")
    pool2 = pd.to_numeric(meta_rows.loc["Pool code 2", series_cols], errors="coerce")
    types = meta_rows.loc["Type", series_cols].astype(str)

    # year rows
    year_mask = wide["Title"].astype(str).str.fullmatch(r"\d{4}")
    year_rows = wide[year_mask].copy()
    if year_rows.empty:
        raise ValueError("No year rows found (expected Title values like 1950, 1951, ...).")

    years = pd.to_numeric(year_rows["Title"], errors="coerce").astype(int)

    meta: Dict[str, RefSeriesMeta] = {}
    long_parts: List[pd.DataFrame] = []

    for col in series_cols:
        t = str(types[col]).strip()
        if INCLUDE_TYPES is not None and t not in INCLUDE_TYPES:
            continue
        if t in EXCLUDE_TYPES:
            continue

        if INCLUDE_SERIES_NAMES and col not in INCLUDE_SERIES_NAMES:
            continue
        if col in EXCLUDE_SERIES_NAMES:
            continue

        if INCLUDE_SERIES_REGEX and not any(re.search(p, col) for p in INCLUDE_SERIES_REGEX):
            continue
        if EXCLUDE_SERIES_REGEX and any(re.search(p, col) for p in EXCLUDE_SERIES_REGEX):
            continue

        pc_tmp = int(pool[col])
        pc2_tmp = None if pd.isna(pool2[col]) else int(pool2[col])

        if INCLUDE_POOL_CODES and pc_tmp not in INCLUDE_POOL_CODES:
            # note: for series where pool_code is a fleet and pool_code2 is a group, you may want to filter on INCLUDE_POOL_CODE2 instead
            continue
        if INCLUDE_POOL_CODE2 and (pc2_tmp is None or pc2_tmp not in INCLUDE_POOL_CODE2):
            continue

        pc = pc_tmp
        pc2 = pc2_tmp
        w = float(weights[col]) if not pd.isna(weights[col]) else 0.0
        disp = SERIES_ALIASES.get(col, col)

        meta[col] = RefSeriesMeta(
            name=col,
            type_name=t,
            pool_code=pc,
            pool_code2=pc2,
            weight=w,
            display_name=disp,
        )

        obs = pd.to_numeric(year_rows[col], errors="coerce")
        yoff = int(SERIES_YEAR_OFFSETS.get(col, 0))
        years_use = (years + yoff).values
        tmp = pd.DataFrame({
            "series": col,
            "type": t,
            "pool_code": pc,
            "pool_code2": pc2,
            "weight": w,
            "year": years_use,
            "obs": obs.values,
        }).dropna(subset=["obs"])

        long_parts.append(tmp)

    if not long_parts:
        return meta, pd.DataFrame(columns=["series", "type", "pool_code", "pool_code2", "weight", "year", "obs"])

    ref_long = pd.concat(long_parts, ignore_index=True)
    ref_long["year"] = ref_long["year"].astype(int)
    return meta, ref_long


def apply_reference_composites(meta: Dict[str, RefSeriesMeta], ref_long: pd.DataFrame) -> Tuple[Dict[str, RefSeriesMeta], pd.DataFrame]:
    """
    Add composite reference series (sum of member series).
    COMPOSITES format:
      { "NewName": {"members": ["SeriesA","SeriesB"], "op": "sum"} }
    """
    if not COMPOSITES:
        return meta, ref_long

    ref = ref_long.copy()
    meta_out = dict(meta)

    for new_name, spec in COMPOSITES.items():
        members = spec.get("members", [])
        op = spec.get("op", "sum")
        if not members:
            continue

        missing = [m for m in members if m not in meta]
        if missing:
            logging.warning(f"Composite {new_name}: missing members {missing}; skipping.")
            continue

        # require same type for now (simplifies mapping to output files)
        type_set = {meta[m].type_name for m in members}
        if len(type_set) != 1:
            logging.warning(f"Composite {new_name}: members have multiple types {type_set}; skipping.")
            continue
        t = next(iter(type_set))

        sub = ref[ref["series"].isin(members)].copy()
        if sub.empty:
            continue

        if op != "sum":
            logging.warning(f"Composite {new_name}: unsupported op={op!r}; only 'sum' is implemented; skipping.")
            continue

        comp = sub.groupby("year", as_index=False)["obs"].sum()
        comp["series"] = new_name
        comp["type"] = t
        comp["pool_code"] = -1
        comp["pool_code2"] = None
        comp["weight"] = float(np.mean([meta[m].weight for m in members]))

        ref = pd.concat([ref, comp], ignore_index=True)

        meta_out[new_name] = RefSeriesMeta(
            name=new_name,
            type_name=t,
            pool_code=-1,
            pool_code2=None,
            weight=float(np.mean([meta[m].weight for m in members])),
            display_name=SERIES_ALIASES.get(new_name, new_name),
        )

    return meta_out, ref


# =========================
# MC output extraction
# =========================

@dataclass(frozen=True)
class Target:
    series: str
    display_name: str
    type_name: str
    pool_code: int
    pool_code2: Optional[int]
    weight: float
    is_composite: bool = False
    members: Optional[List[str]] = None  # for composite targets only


def build_targets(meta: Dict[str, RefSeriesMeta]) -> Dict[str, Target]:
    targets: Dict[str, Target] = {}
    for name, m in meta.items():
        targets[name] = Target(
            series=name,
            display_name=m.display_name,
            type_name=m.type_name,
            pool_code=m.pool_code,
            pool_code2=m.pool_code2,
            weight=m.weight,
        )
    # Mark composites if defined
    for new_name, spec in COMPOSITES.items():
        if new_name in targets:
            targets[new_name] = Target(
                series=new_name,
                display_name=SERIES_ALIASES.get(new_name, new_name),
                type_name=targets[new_name].type_name,
                pool_code=-1,
                pool_code2=None,
                weight=targets[new_name].weight,
                is_composite=True,
                members=list(spec.get("members", [])),
            )
    return targets


def resolve_output_file(trial_dir: Path, type_name: str) -> Optional[Path]:
    spec = OUTPUT_SPECS.get(type_name)
    if spec is None:
        return None
    for fname in spec.filename_candidates:
        p = trial_dir / fname
        if p.exists():
            return p
    return None


def select_model_columns(df: pd.DataFrame, pool_code: int, pool_code2: Optional[int]) -> Optional[List[str]]:
    """
    Heuristic selection of the appropriate output column(s) for a reference series.
    - If pool_code2 is None: expects group-level outputs -> column named like '12'
    - If pool_code2 is not None: expects a combined key -> column containing both codes
      (tries 'pool\pool2', 'pool2\pool', underscores, hyphens, etc.)

    Returns list of column names to sum (usually length 1), or None if not found.
    """
    cols = [str(c) for c in df.columns]

    # Exclude year col from matching attempts (we'll find year separately)
    year_col = find_year_column(cols)
    cols_no_year = [c for c in cols if c != year_col]

    if pool_code2 is None:
        # exact match
        if str(pool_code) in cols_no_year:
            return [str(pool_code)]
        # numeric-ish match (e.g., "12.0")
        for c in cols_no_year:
            try:
                if int(float(c)) == int(pool_code):
                    return [c]
            except Exception:
                pass
        # boundary match inside strings
        token = str(pool_code)
        patt = re.compile(rf"(?<!\d){re.escape(token)}(?!\d)")
        hits = [c for c in cols_no_year if patt.search(c)]
        if len(hits) == 1:
            return hits
        if len(hits) > 1:
            # prefer shortest / simplest
            hits_sorted = sorted(hits, key=len)
            return [hits_sorted[0]]
        return None

    # pool_code2 present: try common concatenations
    a, b = str(pool_code), str(pool_code2)
    candidates = [
        f"{a}\\{b}", f"{b}\\{a}",
        f"{a}/{b}", f"{b}/{a}",
        f"{a}_{b}", f"{b}_{a}",
        f"{a}-{b}", f"{b}-{a}",
        f"{a} {b}", f"{b} {a}",
    ]
    for cand in candidates:
        if cand in cols_no_year:
            return [cand]

    # fallback: columns containing both tokens as standalone numbers
    patt_a = re.compile(rf"(?<!\d){re.escape(a)}(?!\d)")
    patt_b = re.compile(rf"(?<!\d){re.escape(b)}(?!\d)")
    hits = [c for c in cols_no_year if patt_a.search(c) and patt_b.search(c)]
    if len(hits) == 1:
        return hits
    if len(hits) > 1:
        hits_sorted = sorted(hits, key=len)
        return [hits_sorted[0]]
    return None


def load_mc_long(trials_needed_types: Sequence[str], targets: Dict[str, Target]) -> pd.DataFrame:
    """
    Returns long DF: [scenario_id, scenario_label, trial, series, type, year, value]
    Only loads the output files corresponding to types in trials_needed_types,
    and only extracts the columns needed for the targets.
    """
    rows: List[pd.DataFrame] = []

    # group targets by type for efficiency
    targets_by_type: Dict[str, List[Target]] = {}
    for t in targets.values():
        targets_by_type.setdefault(t.type_name, []).append(t)

    for scen in SCENARIOS:
        scen_dir = EWE_OUTPUT_DIR / MODEL_RUN_DIR / scen.mc_subdir
        if not scen_dir.exists():
            logging.warning(f"Scenario directory not found: {scen_dir}")
            continue

        trial_dirs = sorted([p for p in scen_dir.iterdir() if p.is_dir() and p.name.startswith(TRIAL_PREFIX)])
        logging.info(f"{scen.scen_id}: found {len(trial_dirs)} trial folders in {scen_dir}")

        for trial_dir in trial_dirs:
            # cache output dfs per type per trial
            trial_cache: Dict[str, pd.DataFrame] = {}

            for type_name in trials_needed_types:
                if type_name not in targets_by_type:
                    continue

                out_path = resolve_output_file(trial_dir, type_name)
                if out_path is None:
                    logging.debug(f"{scen.scen_id} | {trial_dir.name}: no output file found for type {type_name}")
                    continue

                try:
                    df_out = pd.read_csv(out_path, skiprows=14)
                except Exception as e:
                    logging.warning(f"Failed reading {out_path}: {e}")
                    continue

                # force string column labels for matching
                df_out.columns = [str(c) for c in df_out.columns]
                year_col = find_year_column(df_out.columns)
                if year_col is None:
                    logging.warning(f"No year column found in {out_path}")
                    continue

                trial_cache[type_name] = df_out

                years = pd.to_numeric(df_out[year_col], errors="coerce").astype("Int64")

                # extract each target of this type
                for tgt in targets_by_type[type_name]:
                    if tgt.is_composite:
                        continue  # composites built later
                    cols = select_model_columns(df_out, tgt.pool_code, tgt.pool_code2)
                    if not cols:
                        logging.debug(
                            f"{scen.scen_id} | {trial_dir.name} | {tgt.series}: no matching output col for "
                            f"pool_code={tgt.pool_code}, pool_code2={tgt.pool_code2}"
                        )
                        continue

                    vals = df_out[cols].apply(pd.to_numeric, errors="coerce").sum(axis=1)
                    tmp = pd.DataFrame({
                        "scenario_id": scen.scen_id,
                        "scenario_label": scen.label,
                        "trial": trial_dir.name,
                        "series": tgt.series,
                        "type": tgt.type_name,
                        "year": years,
                        "value": vals,
                    }).dropna(subset=["year", "value"])

                    tmp["year"] = tmp["year"].astype(int)
                    rows.append(tmp)

    if not rows:
        return pd.DataFrame(columns=["scenario_id", "scenario_label", "trial", "series", "type", "year", "value"])
    return pd.concat(rows, ignore_index=True)


def apply_mc_composites(mc_long: pd.DataFrame, targets: Dict[str, Target]) -> pd.DataFrame:
    """
    Create composite series by summing member series per (scenario, trial, year).
    """
    if not COMPOSITES:
        return mc_long

    out = mc_long.copy()

    for comp_name, spec in COMPOSITES.items():
        members = spec.get("members", [])
        if not members:
            continue
        if comp_name not in targets or not targets[comp_name].is_composite:
            continue

        # require same type
        t = targets[comp_name].type_name
        sub = out[(out["type"] == t) & (out["series"].isin(members))].copy()
        if sub.empty:
            logging.warning(f"Composite {comp_name}: no MC data found for members; skipping.")
            continue

        comp = (
            sub.groupby(["scenario_id", "scenario_label", "trial", "year"], as_index=False)["value"]
            .sum()
        )
        comp["series"] = comp_name
        comp["type"] = t
        out = pd.concat([out, comp], ignore_index=True)

    return out


# =========================
# Relative normalisation + summaries
# =========================

def get_baseline_years(ref_years: np.ndarray, model_years: np.ndarray) -> np.ndarray:
    """
    Determine baseline years to normalise relative series.
    """
    if REL_BASELINE_YEARS is not None:
        y0, y1 = REL_BASELINE_YEARS
        years = np.arange(y0, y1 + 1)
        return years

    overlap = np.intersect1d(ref_years, model_years)
    if overlap.size == 0:
        # fallback to model years
        return model_years

    if REL_BASELINE_MODE == "first_year":
        return np.array([int(overlap.min())])

    # default: overlap mean
    return overlap


def normalise_relative(
    ref_df: pd.DataFrame,
    mc_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """
    Normalise both ref and model (per trial) to baseline mean (or first year).
    Returns (ref_norm, mc_norm, baseline_info)
    """
    ref = ref_df.copy()
    mc = mc_df.copy()

    ref_years = ref["year"].to_numpy(dtype=int)
    model_years = mc["year"].to_numpy(dtype=int)

    baseline_years = get_baseline_years(ref_years, model_years)
    baseline_info = {"baseline_mode": REL_BASELINE_MODE, "baseline_years_min": int(baseline_years.min()), "baseline_years_max": int(baseline_years.max())}

    # observed baseline
    ref_base = ref[ref["year"].isin(baseline_years)]["obs"]
    ref_scale = float(ref_base.mean()) if ref_base.notna().any() else np.nan

    # per-trial baseline for model
    mc["value_norm"] = np.nan
    for trial, sub in mc.groupby("trial"):
        base = sub[sub["year"].isin(baseline_years)]["value"]
        scale = float(base.mean()) if base.notna().any() else np.nan
        if np.isfinite(scale) and scale != 0:
            mc.loc[sub.index, "value_norm"] = sub["value"] / scale

    if np.isfinite(ref_scale) and ref_scale != 0:
        ref["obs_norm"] = ref["obs"] / ref_scale
    else:
        ref["obs_norm"] = np.nan

    return ref, mc, baseline_info


def mc_quantiles_by_year(mc_df: pd.DataFrame, use_norm: bool) -> pd.DataFrame:
    """
    Summarise MC distribution by year using nanpercentile across trials.
    mc_df contains columns: [year, trial, value] (and possibly value_norm)
    """
    val_col = "value_norm" if use_norm else "value"
    sub = mc_df[["year", "trial", val_col]].dropna(subset=[val_col]).copy()
    if sub.empty:
        return pd.DataFrame(columns=["year"] + [f"q{q}" for q in QUANTILES])

    pivot = sub.pivot_table(index="year", columns="trial", values=val_col, aggfunc="first").sort_index()
    arr = pivot.to_numpy(dtype=float)

    q = np.nanpercentile(arr, QUANTILES, axis=1)
    out = pd.DataFrame({"year": pivot.index.values})
    for i, qq in enumerate(QUANTILES):
        out[f"q{qq}"] = q[i, :]
    return out


# =========================
# Plotting + stats
# =========================

def plot_one_series(
    series: str,
    targets: Dict[str, Target],
    ref_long: pd.DataFrame,
    mc_long: pd.DataFrame
) -> List[Dict]:
    """
    Plot a single reference series across scenarios.
    Returns fit-stat rows for this series (possibly empty).
    """
    tgt = targets[series]
    tname = tgt.type_name
    is_relative = (tname in RELATIVE_TYPES) and REL_NORM_ENABLED

    ref_df = ref_long[ref_long["series"] == series][["year", "obs"]].sort_values("year").copy()
    if ref_df.empty:
        logging.warning(f"{series}: no reference data; skipping plot.")
        return []

    fit_rows: List[Dict] = []

    plt.figure(figsize=(12, 6))

    # Plot observed (possibly normalised)
    # We'll compute normalisation per scenario after we have model years; to keep plots consistent across scenarios,
    # we normalise using the union of model years across all scenarios/trials available for this series.
    mc_all = mc_long[(mc_long["series"] == series) & (mc_long["type"] == tname)].copy()

    if mc_all.empty:
        logging.warning(f"{series}: no MC model data found for type={tname}; skipping plot.")
        return []

    ref_plot = ref_df.copy()
    mc_plot_all = mc_all.copy()
    baseline_info = {}

    if is_relative:
        ref_plot, mc_plot_all, baseline_info = normalise_relative(ref_plot, mc_plot_all)
        plt.scatter(ref_plot["year"], ref_plot["obs_norm"], color="black", label="Observed (normalised)", zorder=3)
        y_label = f"{tname} (normalised)"
    else:
        plt.scatter(ref_plot["year"], ref_plot["obs"], color="black", label="Observed", zorder=3)
        y_label = tname

    # Now per scenario envelopes
    for scen in SCENARIOS:
        mc_s = mc_plot_all[mc_plot_all["scenario_id"] == scen.scen_id].copy()
        if mc_s.empty:
            continue

        summ = mc_quantiles_by_year(mc_s, use_norm=is_relative)
        if summ.empty:
            continue

        years = summ["year"].values
        q5, q25, q50, q75, q95 = (summ["q5"].values, summ["q25"].values, summ["q50"].values, summ["q75"].values, summ["q95"].values)

        plt.fill_between(years, q5, q95, alpha=0.15, label=f"{scen.label} 5–95%")
        plt.fill_between(years, q25, q75, alpha=0.25, label=f"{scen.label} 25–75%")
        plt.plot(years, q50, linewidth=2, label=f"{scen.label} median")

        # Fit stats on overlapping years
        if COMPUTE_FIT_STATS:
            ref_col = "obs_norm" if is_relative else "obs"
            mod_col = "q50"

            merged = pd.merge(
                ref_plot[["year", ref_col]].rename(columns={ref_col: "obs_use"}),
                summ[["year", mod_col]].rename(columns={mod_col: "mod_use"}),
                on="year",
                how="inner",
            ).dropna()

            if not merged.empty:
                obs = merged["obs_use"].to_numpy(dtype=float)
                mod = merged["mod_use"].to_numpy(dtype=float)
                row = {
                    "series": series,
                    "display_name": tgt.display_name,
                    "type": tname,
                    "is_relative": bool(is_relative),
                    "weight": tgt.weight,
                    "scenario_id": scen.scen_id,
                    "scenario_label": scen.label,
                    "n_years": int(len(merged)),
                    "rmse": rmse(obs, mod),
                    "mae": mae(obs, mod),
                    "bias_model_minus_obs": float(np.nanmean(mod - obs)),
                    "year_min": int(merged["year"].min()),
                    "year_max": int(merged["year"].max()),
                }
                row.update(baseline_info)
                fit_rows.append(row)

    title = f"{tgt.display_name}  |  {tname}"
    plt.title(title)
    plt.xlabel("Year")
    plt.ylabel(y_label)
    plt.grid(True)
    plt.legend(loc="best", fontsize=9)
    if YEAR_MIN is not None and YEAR_MAX is not None:
        plt.xlim(YEAR_MIN, YEAR_MAX)
    plt.tight_layout()

    if SAVE_PLOTS:
        safe_mkdir(OUT_DIR)
        fname = OUT_DIR / f"mc_{sanitize_filename(tname)}_{sanitize_filename(tgt.display_name)}.png"
        plt.savefig(fname, dpi=DPI)
        logging.info(f"Saved: {fname}")

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()

    return fit_rows


def main() -> None:
    setup_logging()

    if not TIME_SERIES_FILE.exists():
        logging.error(f"TIME_SERIES_FILE not found: {TIME_SERIES_FILE}")
        return

    logging.info(f"TimeSeries file: {TIME_SERIES_FILE}")
    logging.info(f"EWE_OUTPUT_DIR:  {EWE_OUTPUT_DIR}")
    logging.info(f"MODEL_RUN_DIR:   {MODEL_RUN_DIR}")
    logging.info(f"OUT_DIR:         {OUT_DIR}")

    meta, ref_long = parse_ewe_timeseries_wide(TIME_SERIES_FILE)
    meta, ref_long = apply_reference_composites(meta, ref_long)
    targets = build_targets(meta)

    if ref_long.empty:
        logging.error("No reference series loaded after filters. Check INCLUDE_TYPES / INCLUDE_SERIES_NAMES.")
        return

    # Determine which types we need to load from MC outputs
    needed_types = sorted({t.type_name for t in targets.values() if t.type_name in OUTPUT_SPECS})
    if not needed_types:
        logging.error("No target types matched OUTPUT_SPECS. Update OUTPUT_SPECS mapping.")
        return

    mc_long = load_mc_long(needed_types, targets)
    if mc_long.empty:
        logging.error("No MC data loaded. Check scenario paths and output filename mapping in OUTPUT_SPECS.")
        return

    mc_long = apply_mc_composites(mc_long, targets)

    # Plot each series
    fit_all: List[Dict] = []
    for series in targets.keys():
        # only plot series that exist in ref_long (some targets might be missing if composite skipped)
        if series not in set(ref_long["series"].unique()):
            continue
        fit_all.extend(plot_one_series(series, targets, ref_long, mc_long))

    # Save fit stats
    if COMPUTE_FIT_STATS and fit_all:
        safe_mkdir(OUT_DIR)
        fit_df = pd.DataFrame(fit_all).sort_values(["type", "series", "scenario_id"])
        out_path = OUT_DIR / FIT_STATS_FILENAME
        fit_df.to_csv(out_path, index=False)
        logging.info(f"Saved fit stats: {out_path}")


if __name__ == "__main__":
    main()
