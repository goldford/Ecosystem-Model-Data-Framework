#!/usr/bin/env python3
"""
-------------------
Standalone converter from a herring spawn CSV (with Ecospace grid indices)
to Ecospace-ready ASC rasters. Grinnell et al. 2023;
https://open.canada.ca/data/en/dataset/d892511c-d851-4f85-a0ec-708bc05d2810

Created by G Oldford Sep 2025

Purpose: create a series of ASC files that use the raw data on spawn locations to convert
         to a 'habitat capacity' layer for EwE Ecospace to ingest.

Inputs:
1) Time series:  Pacific_SoG_herring_spawn_index_data_2024_EN_taggedEwERowsCols.csv

Outputs:
1) Time series:   HERRING_SPAWN_YYYY_MM.asc  (per-month density; sums to 1 over water)
2) All-years:     HERRING_SPAWN_BASE_ALLYEARS.asc  (across all months & years; sums to 1)
3) Climatology:   HERRING_SPAWN_CLIM_01..12.asc  (per-month across-years; sums to 1)

“Probability per month across all years” = these CLIM_* maps.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from helpers import saveASCFile
from helpers_ewe_asc_netcdf import (
    read_asc_header_and_dims,
    load_mask_asc,
    normalize_to_one,
    _normalize_grid,
    month_iter
)

# ---------------------------
# defaults
# ---------------------------

DEFAULT_SPAWN_P      = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//")
DEFAULT_SPAWN_FP     = Path(DEFAULT_SPAWN_P, "Pacific_SoG_herring_spawn_index_data_2024_EN_taggedEwERowsCols.csv")
DEF_OUT_P            = Path(DEFAULT_SPAWN_P, "SPAWN_DATA_ASC")

# example ASC used only for header/dimensions
DEF_TEMPLATE_ASC_P   = Path("..//data//forcing//ECOSPACE_in_climayr_1980-2018_PAR3_Sal4m_20240523//PAR-VarZ-VarK//")
DEF_TEMPLATE_ASC_FP  = Path(DEF_TEMPLATE_ASC_P, "PAR-VarZ-VarK_1980-2018.asc")

# optional 0=land ASC (nonzero = water)
DEF_LAND_P           = Path("..//data//")
DEF_LAND_F           = "ecospacedepthgrid.asc"
DEF_LAND_PF          = Path(DEF_LAND_P, DEF_LAND_F)

# spawn/value/time parsing
DEF_BINARY_MODE      = False    # set True to make each month 0/1 before normalization
DEF_NORMALISE        = True   # have a binary map or a 'probability' map (normalise ensures probabilities, ie map sums to 1)
DEF_NORM_METHOD      = "max"   # "sum" (across map) or "max" (maximum val across map within each month)
DEF_BIN_THRESHOLD    = 0.0     # > threshold counts as presence
DEF_SPWN_COL         = "RelativeEggB"
DEF_DATE_COL         = "StartDate"
DEF_YEAR             = None  # keep; not used by default (process all years)

# export/normalization knobs

DEF_FMT              = "%0.2f"
DEF_INDEX_BASE_EWE   = 1        # EwE 1-indexed
DEF_BASEOUTNAME      = "HERRING_SPAWN"
DEF_MIN_FILL         = 0.05     # your previous min fill
DEF_MIN_FLOOR        = DEF_MIN_FILL
DEF_LOG_NORM         = False
DEF_LOG_EPS          = 1e-6


# ---------------------------
# Utils
# ---------------------------

def infer_year_month(df: pd.DataFrame, date_col: str | None):
    if {"Year","Month"}.issubset(df.columns):
        out = df.copy()
        out["Year"] = pd.to_numeric(out["Year"], errors="coerce").astype("Int64")
        out["Month"] = pd.to_numeric(out["Month"], errors="coerce").astype("Int64")
        return out.dropna(subset=["Year","Month"])
    if date_col is None or date_col not in df.columns:
        raise KeyError("Provide Year/Month columns or a valid --date_col.")
    out = df.copy()
    dt = pd.to_datetime(out[date_col], errors="coerce")
    out["Year"] = dt.dt.year.astype("Int64")
    out["Month"] = dt.dt.month.astype("Int64")
    return out.dropna(subset=["Year","Month"])

def choose_value_col(df: pd.DataFrame, value_col: str | None):
    if value_col and value_col in df.columns:
        return value_col
    # common fallbacks used in previous scripts
    for c in ["RelativeEggB", "SpawnIntensity", "value", "intensity"]:
        if c in df.columns:
            return c
    raise KeyError("Could not find a spawn intensity column. Use --value_col.")


# ---------------------------
# Main
# ---------------------------

def csv_to_ascs(
    spawn_csv: Path,
    out_dir: Path,
    example_asc: Path,
    land_mask_asc: Path | None = None,   # optional 0=land mask in ASC form
    value_col: str | None = None,
    date_col: str | None = "StartDate",
    index_base: int = 1,
    fmt: str = "%0.2f",
    min_floor: float = 0.0,
    log_norm: bool = False,
    log_eps: float = 1e-6,
    quiet: bool = False,
    basename: str = "HERRING_SPAWN",
    norm: bool = False,
    norm_method: str = "max",
    binary_mode: bool = True,
    bin_threshold: float = 0.0
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not quiet:
        print(f"[INFO] Reading CSV: {spawn_csv}")
    df = pd.read_csv(spawn_csv)

    for req in ["ecospace_row", "ecospace_col"]:
        if req not in df.columns:
            raise KeyError(f"CSV missing required column: {req}")

    val_col = choose_value_col(df, value_col)
    df["_val"] = pd.to_numeric(df[val_col], errors="coerce").fillna(0.0)

    df = infer_year_month(df, date_col=date_col)
    df["ecospace_row"] = pd.to_numeric(df["ecospace_row"], errors="coerce").astype("Int64")
    df["ecospace_col"] = pd.to_numeric(df["ecospace_col"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["ecospace_row", "ecospace_col"])

    header, nrows, ncols = read_asc_header_and_dims(example_asc)
    include = np.ones((nrows, ncols), dtype=bool)

    # Optional water mask from land-depth ASC (0 = land)
    if land_mask_asc and Path(land_mask_asc).exists():
        include = load_mask_asc(land_mask_asc, nrows, ncols)

    # Accumulators
    all_years_accum = np.zeros((nrows, ncols), dtype=float)
    # For climatology (month across years)
    clim_accum = {m: np.zeros((nrows, ncols), dtype=float) for m in month_iter()}

    years = sorted(int(y) for y in df["Year"].dropna().unique())
    if not years:
        raise ValueError("No valid Year values found after parsing.")

    for yy in years:
        df_y = df[df["Year"] == yy]
        for mm in month_iter():
            sub = df_y[df_y["Month"] == mm]
            grid = np.zeros((nrows, ncols), dtype=float)
            if not sub.empty:
                rows = (sub["ecospace_row"].to_numpy() - index_base).astype(int)
                cols = (sub["ecospace_col"].to_numpy() - index_base).astype(int)
                vals = sub["_val"].to_numpy(float)
                ok = (rows >= 0) & (rows < nrows) & (cols >= 0) & (cols < ncols)
                rows, cols, vals = rows[ok], cols[ok], vals[ok]
                np.add.at(grid, (rows, cols), vals)

            if binary_mode:
                grid = (grid > bin_threshold).astype(float)

            # optional log compression before normalization
            if log_norm:
                pos = grid > 0
                grid[pos] = np.log(grid[pos] + log_eps)

            # normalize to 1 over water; apply floor
            if norm:
                out_grid = _normalize_grid(grid, include, norm_method, min_floor=min_floor)
            else:
                out_grid = grid.copy()

            # write year-month ASC
            out_name = f"{basename}_{yy}_{mm:02d}.asc"
            saveASCFile(str(out_dir / out_name),
                        np.flipud(out_grid),
                        bottomleft_row_ewe=0,
                        upperleft_row_ewe=nrows,
                        upperleft_col_ewe=0,
                        ASCheader=header,
                        sigdigfmt=fmt,
                        dfPlumeMask=None,
                        dfLandMask=None)

            if not quiet:
                s = out_grid[include].sum()
                print(f"[INFO] Wrote {out_name} (sum over water={s:.6f})")

            # accumulate for BASE and CLIM
            all_years_accum += out_grid
            clim_accum[mm] += out_grid

    # -------- All-years BASE (normalize once to 1 over water) ----------
    if norm:
        base = _normalize_grid(all_years_accum, include, norm_method, min_floor=min_floor)
    else:
        base = all_years_accum.copy()

    saveASCFile(str(out_dir / f"{basename}_BASE_ALLYEARS.asc"),
                np.flipud(base),
                bottomleft_row_ewe=0,
                upperleft_row_ewe=nrows,
                upperleft_col_ewe=0,
                ASCheader=header,
                sigdigfmt=fmt,
                dfPlumeMask=None,
                dfLandMask=None)

    if not quiet:
        print("[INFO] Wrote BASE_ALLYEARS")

    # -------- Monthly CLIMATOLOGY (across years) -----------------------
    # For each month, normalize its across-years sum to 1 over water
    for mm in month_iter():
        if norm:
            # clim = normalize_to_one(clim_accum[mm], include, min_floor=min_floor)
            clim = _normalize_grid(clim_accum[mm], include, norm_method, min_floor=min_floor)
        else:
            clim = clim_accum[mm].copy()

        saveASCFile(str(out_dir / f"{basename}_CLIM_{mm:02d}.asc"),
                    np.flipud(clim),
                    bottomleft_row_ewe=0,
                    upperleft_row_ewe=nrows,
                    upperleft_col_ewe=0,
                    ASCheader=header,
                    sigdigfmt=fmt,
                    dfPlumeMask=None,
                    dfLandMask=None)


    if not quiet:
        print("[INFO] Wrote monthly CLIM_01..12")

def main(argv=None):
    ap = argparse.ArgumentParser(description="CSV→ASC herring spawn generator (time series + all-years + monthly climatology).")
    ap.add_argument("--spawn_csv",      type=Path, default=DEFAULT_SPAWN_FP)
    ap.add_argument("--out_dir",        type=Path, default=DEF_OUT_P)
    ap.add_argument("--example_asc",    type=Path, default=DEF_TEMPLATE_ASC_FP)
    ap.add_argument("--land_mask_asc",  type=Path, default=DEF_LAND_PF)
    ap.add_argument("--value_col",      type=str,  default=DEF_SPWN_COL)
    ap.add_argument("--date_col",       type=str,  default=DEF_DATE_COL)
    ap.add_argument("--index_base",     type=int,  default=DEF_INDEX_BASE_EWE, choices=[0,1])
    ap.add_argument("--fmt",            type=str,  default=DEF_FMT)
    ap.add_argument("--min_floor",      type=float,default=DEF_MIN_FLOOR)
    ap.add_argument("--log_norm",       action="store_true", default=DEF_LOG_NORM)
    ap.add_argument("--log_eps",        type=float,default=DEF_LOG_EPS)
    ap.add_argument("--basename",       type=str,  default=DEF_BASEOUTNAME)
    ap.add_argument("--norm_method", type=str, choices=["sum", "max"], default=DEF_NORM_METHOD)
    ap.add_argument("--normalise_map",    default=DEF_NORMALISE)
    ap.add_argument("--binary_mode",    default=DEF_BINARY_MODE)
    ap.add_argument("--bin_threshold", type=float, default=DEF_BIN_THRESHOLD)
    args = ap.parse_args(argv)


    csv_to_ascs(
        spawn_csv=args.spawn_csv,
        out_dir=args.out_dir,
        example_asc=args.example_asc,
        land_mask_asc=args.land_mask_asc,
        value_col=args.value_col,
        date_col=args.date_col,
        index_base=args.index_base,
        fmt=args.fmt,
        min_floor=args.min_floor,
        log_norm=args.log_norm,
        log_eps=args.log_eps,
        basename=args.basename,
        norm=args.normalise_map,
        norm_method=args.norm_method,
        binary_mode=args.binary_mode,
        bin_threshold=args.bin_threshold
    )

if __name__ == "__main__":
    main()
