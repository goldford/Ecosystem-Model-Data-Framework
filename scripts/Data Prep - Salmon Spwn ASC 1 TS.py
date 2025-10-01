"""
G Oldford Sep 2025

Data Prep – Salmon Monthly ASC from Migration Probabilities (Wilson & Peacock, 2025; 10.1139/cjfas-2024-0213)
Creates 12 monthly Ecospace ASC layers for a salmon species (“CM”, type “rt”)
from a migration-probabilities table, with coastal/offshore emphasis
overwritten as (1 - monthly_probability) in the indicated cells.

Logic
-----
For each month m:
  grid[water] = p_m
  grid[(chum_rivers==1 or offshore_areas==1) & water] = 1 - p_m
  (land remains 0 via writer’s dfLandMask semantics, as in your pipeline)

Inputs
------
- migration CSV with monthly probabilities (supports long or wide format)
- template ASC: provides header/nrows/ncols & orientation conventions
- NC bundle: contains static masks (depth, chum_rivers, offshore_areas)
  (depth>0 is treated as water)

Outputs
-------
12 ASC files:  CHUM_CM_rt_01.asc ... CHUM_CM_rt_12.asc

Dependencies
------------
- pandas, numpy, netCDF4
- helpers already used:
    - saveASCFile(path, array2d_flipped, ..., ASCheader=header, dfLandMask=..., dfPlumeMask=...)
    - getDataFrame(path, nodata_str)  [used only to build land mask for writer semantics]
"""


from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset

# existing utilities
from helpers import saveASCFile, getDataFrame
from helpers_ewe_asc_netcdf import smooth_gaussian_masked

# ------------------------------
# DEFAULTS
# ------------------------------

DEF_INPUT_MIGPROB_CSV = "..//..//ecospace_seal_diets//Data//migration_probabilities_normalisedforEwE.csv"
DEF_EXAMP_ASC = "..//data//basemap//ecospacedepthgrid.asc"
DEF_NC_PATH = "..//..//ecospace_seal_diets//Data//basemaps//ecospace_bundled_asc.nc" # multi ASC stored in this NC

DEF_SPECIES = "CHUM"
DEF_SPECIES_CODE = "CM"
DEF_TYPE_CODE = "rt" # rt for returning, oe for ocean entrant
DEF_BASENAME = DEF_SPECIES + "_" + DEF_SPECIES_CODE + "_" + DEF_TYPE_CODE

DEF_OUT_DIR = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//1. Salmon//All Species Run Timing Wilson Peacock 2025//MODIFIED//ASC for Adult Run Timing//",DEF_SPECIES +" SPWN//")

DEF_FMT = "%0.3f"
DEF_FLOOR = 0.01
DEF_DEP_VAR = "depth"

DEF_RIVERS_VAR = DEF_SPECIES.lower() + "_" + "rivers"
DEF_RIVERS_VAR = "chum" + "_" + "rivers" # TEMP - only have chum right now

DEF_RIV_MOS = "8-12" # restrict application of river mouth to habitat cap to certain months?
DEF_OFFSHORE_VAR = "offshore_areas_chum" # TO DO - only have chum r now
DEF_OFFSHORE_MOS = "1-7"

DEF_GAUSS_SIGMA = 2.0 # cells for gaussian smoothing

# ------------------------------
# ASC header util
# ------------------------------

def read_asc_header_and_dims(example_asc: Path) -> tuple[str, int, int]:
    """Read first 6 lines of an ASC to get header text and (nrows, ncols)."""
    with open(example_asc, "r") as f:
        header_lines = [next(f).rstrip("\n") for _ in range(6)]
    header = "\n".join(header_lines)

    def _find_int(key: str) -> int:
        for line in header_lines:
            if line.lower().startswith(key):
                return int(line.split()[-1])
        raise ValueError(f"Could not find {key} in ASC header.")

    ncols = _find_int("ncols")
    nrows = _find_int("nrows")
    return header, nrows, ncols


# add near other imports
def parse_months(s: str | None) -> set[int]:
    """
    Parse month spec into a set of ints 1..12.
    Examples:
      "all"                     -> {1..12}
      "none" or ""              -> set()
      "1,2,3" or "01,02,03"     -> {1,2,3}
      "3-6"                     -> {3,4,5,6}
      "1,3-5,12"                -> {1,3,4,5,12}
    """
    if not s or s.lower() == "none":
        return set()
    if s.lower() == "all":
        return set(range(1,13))
    out = set()
    for part in s.split(","):
        part = part.strip()
        if "-" in part:
            a,b = part.split("-",1)
            a, b = int(a), int(b)
            out.update(range(min(a,b), max(a,b)+1))
        else:
            out.add(int(part))
    # clamp to 1..12 just in case
    return {m for m in out if 1 <= m <= 12}

# ------------------------------
# CSV monthly probabilities
# ------------------------------

def extract_monthly_probs(df: pd.DataFrame,
                          species: str = "CM",
                          typ: str = "rt") -> dict[int, float]:
    """
    Return {1..12: p} for given species/type.
    Accepts long (month/probability) or wide (Jan..Dec, M01..M12, prob_01..prob_12).
    """
    df2 = df.copy()
    lower_cols = {c.lower(): c for c in df2.columns}
    if "species" not in lower_cols or "type" not in lower_cols:
        raise KeyError("CSV must include 'species' and 'type' columns (any case).")

    df2 = df2[df2[lower_cols["species"]] == species]
    df2 = df2[df2[lower_cols["type"]] == typ]
    if df2.empty:
        raise ValueError(f"No rows for species='{species}' and type='{typ}'")

    # Long form?
    month_col = next((c for c in df2.columns if c.lower() == "month"), None)
    prob_col = next((c for c in df2.columns if c.lower() in ("probability","prob")), None)
    if month_col and prob_col:
        df2 = df2.dropna(subset=[month_col, prob_col]).copy()
        df2[month_col] = pd.to_numeric(df2[month_col], errors="coerce").astype(int)
        out = {}
        for m in range(1, 13):
            sub = df2[df2[month_col] == m]
            if sub.empty:
                raise ValueError(f"No probability for month {m} in long-form CSV.")
            out[m] = float(sub.iloc[0][prob_col])
        return out

    # Wide form
    month_aliases = {
        **{f"prob_{i:02d}": i for i in range(1,13)},
        **{f"m{i:02d}": i for i in range(1,13)},
        **{name.lower(): i+1 for i,name in enumerate(
            ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
        )},
        **{str(i): i for i in range(1,13)}
    }
    mapping = {}
    for c in df2.columns:
        low = c.strip().lower()
        if low in month_aliases:
            mapping[month_aliases[low]] = c
        elif low.startswith("prob_") and low.split("prob_")[-1].isdigit():
            m = int(low.split("prob_")[-1])
            if 1 <= m <= 12: mapping[m] = c
        elif low.startswith("prob") and low[4:].isdigit():
            m = int(low[4:])
            if 1 <= m <= 12: mapping[m] = c

    if len(mapping) < 12:
        raise ValueError("Could not find 12 monthly probability columns in wide-form CSV.")

    row = df2.iloc[0]
    return {m: float(row[col]) for m, col in mapping.items()}


# ------------------------------
# Read masks from NetCDF
# ------------------------------

def load_masks_from_nc(nc_path: Path,
                       depth_var: str = "depth",
                       rivers_var: str = "chum_rivers",
                       offshore_var: str = "offshore_areas") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns (water_mask, rivers_mask, offshore_mask) from the NC bundle.
    - water_mask: depth > 0
    - rivers_mask: (var == 1)
    - offshore_mask: (var == 1)
    Assumes (y,x) layout for these 2D statics.
    """
    with Dataset(nc_path, "r") as ds:
        if depth_var not in ds.variables:
            raise KeyError(f"'{depth_var}' not found in {nc_path}")
        depth = ds.variables[depth_var][:]
        # squeeze and ensure 2D
        depth2d = np.array(depth).squeeze()
        if depth2d.ndim != 2:
            raise ValueError(f"{depth_var} must be 2D; got shape {depth2d.shape}")
        water_mask = depth2d > 0

        def _fetch_mask(name: str) -> np.ndarray:
            if name not in ds.variables:
                print(ds.variables)
                raise KeyError(f"'{name}' not found in {nc_path}")
            arr = np.array(ds.variables[name][:]).squeeze()
            if arr.ndim != 2:
                raise ValueError(f"{name} must be 2D; got shape {arr.shape}")
            return arr == 1

        rivers_mask = _fetch_mask(rivers_var)
        offshore_mask = _fetch_mask(offshore_var)

    return water_mask, rivers_mask, offshore_mask


# ------------------------------
# Main builder
# ------------------------------

def build_salmon_monthlies_from_nc(migration_csv: Path,
                                 example_asc: Path,
                                 nc_path: Path,
                                 out_dir: Path,
                                 species: str = "CM",
                                 typ: str = "rt",
                                 basename: str = "CHUM_CM_rt",
                                 fmt: str = "%0.3f",
                                 floor: float = 0.01,
                                 depth_var: str = "depth",
                                 rivers_var: str = "chum_rivers",
                                 river_months: set[int] | None = None,
                                 offshore_var: str = "offshore_areas",
                                 offshore_months: set[int] | None = None,
                                 gauss_sigma: int = None,
                                 verbose: bool = True):
    # Probabilities
    if verbose:
        print(f"[INFO] Reading migration table: {migration_csv}")
    df = pd.read_csv(migration_csv)
    probs = extract_monthly_probs(df, species=species, typ=typ)
    if verbose:
        print("[INFO] Monthly probabilities:", {k: round(v,5) for k,v in sorted(probs.items())})

    # ASC header / dims
    header, nrows, ncols = read_asc_header_and_dims(example_asc)
    if verbose:
        print(f"[INFO] Grid dims from template: {nrows} x {ncols}")

    # Land mask for writer semantics (land==0 in  usual depth ASC);
    # we won’t use it for logic, only to keep writer behavior.
    # If prefer not to read any ASC at all here, can pass None to dfLandMask,
    # but keeping it preserves  pipeline’s exact handling.
    df_land_for_writer = getDataFrame(str(example_asc), "-9999.00000000")
    if df_land_for_writer.shape != (nrows, ncols):
        raise ValueError("Template ASC shape mismatch with its own header (unexpected).")
    land_writer_mask = (df_land_for_writer == 0)

    # NC masks
    water_mask, rivers_mask, offshore_mask = load_masks_from_nc(
        nc_path, depth_var=depth_var, rivers_var=rivers_var, offshore_var=offshore_var
    )
    if water_mask.shape != (nrows, ncols):
        raise ValueError(f"NC mask shape {water_mask.shape} != template {(nrows, ncols)}")
    if rivers_mask.shape != (nrows, ncols) or offshore_mask.shape != (nrows, ncols):
        raise ValueError("NC river/offshore mask shape mismatch with template.")

    overwrite_mask = (rivers_mask | offshore_mask) & water_mask

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build/write per month
    for m in range(1, 13):
        p = float(probs[m])
        grid = np.zeros((nrows, ncols), dtype=float)

        grid[water_mask] = p

        # month-aware overwrite:
        rivers_active = (m in river_months) if river_months is not None else True
        offshore_active = (m in offshore_months) if offshore_months is not None else True
        month_overwrite_mask = ((offshore_mask if offshore_active else False) |
                                (rivers_mask if rivers_active else False)) & water_mask

        grid[month_overwrite_mask] = 1.0 - p

        # Optional smoothing
        if gauss_sigma is not None and gauss_sigma > 0:
            grid = smooth_gaussian_masked(grid, water_mask.astype(float), gauss_sigma)

        np.clip(grid, 0.0, 1.0, out=grid)

        # enforce a minimum on water cells (land remains 0)
        if floor is not None and floor > 0:
            grid[water_mask] = np.maximum(grid[water_mask], floor)

        out_path = out_dir / f"{basename}_{m:02d}.asc"
        if verbose:
            print(f"[INFO] Writing {out_path}  (p={p:.3f} ; overwrite to {1.0-p:.3f})")

        # Flip vertically if that’s what saveASCFile expects (matching prior workflow)
        saveASCFile(
            str(out_path),
            np.flipud(grid),
            bottomleft_row_ewe=0,
            upperleft_row_ewe=nrows,
            upperleft_col_ewe=0,
            sigdigfmt=fmt,
            ASCheader=header,
            dfPlumeMask=None,
            dfLandMask=land_writer_mask
        )

    if verbose:
        print("[INFO] Done creating monthly salmon layers from NC masks.")



# ------------------------------
# CLI
# ------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description="Create 12 monthly salmon ASC maps using masks from a NetCDF bundle.")
    ap.add_argument("--migration_csv", type=Path, default=DEF_INPUT_MIGPROB_CSV, help="CSV with monthly probabilities (must include species/type).")
    ap.add_argument("--example_asc", type=Path, default=DEF_EXAMP_ASC, help="Template ASC for header/nrows/ncols/orientation.")
    ap.add_argument("--nc_path", type=Path, default=DEF_NC_PATH, help="NetCDF bundle with static masks (depth, rivers, offshore_areas).")
    ap.add_argument("--out_dir", type=Path, default=DEF_OUT_DIR, help="Output directory for monthly ASCs.")
    ap.add_argument("--species", type=str, default=DEF_SPECIES_CODE)
    ap.add_argument("--type", dest="typ", type=str, default=DEF_TYPE_CODE)
    ap.add_argument("--basename", type=str, default=DEF_BASENAME)
    ap.add_argument("--fmt", type=str, default=DEF_FMT)
    ap.add_argument("--depth_var", type=str, default=DEF_DEP_VAR)
    ap.add_argument("--rivers_var", type=str, default=DEF_RIVERS_VAR)
    ap.add_argument("--river_months", type=str, default=DEF_RIV_MOS,
                    help='Months to apply river overwrite. Examples: "all", "none", "1,2,3", "3-6", "1,3-5,12". Default "all".')
    ap.add_argument("--offshore_var", type=str, default=DEF_OFFSHORE_VAR)
    ap.add_argument("--offshore_months", type=str, default=DEF_OFFSHORE_MOS,
                    help='Months to apply offshore area overwrite. Examples: "all", "none", "1,2,3", "3-6", "1,3-5,12". Default "all".')
    ap.add_argument("--floor", type=float, default=DEF_FLOOR, help="Minimum value to enforce on water cells (e.g., 0.001).")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--gauss_sigma", type=float, default=DEF_GAUSS_SIGMA,
                    help="Std dev (cells) for Gaussian smoothing, e.g., 1.0 or 2.0.")
    args = ap.parse_args(argv)

    river_months = parse_months(args.river_months)
    offshore_months = parse_months(args.offshore_months)

    build_salmon_monthlies_from_nc(
        migration_csv=args.migration_csv,
        example_asc=args.example_asc,
        nc_path=args.nc_path,
        out_dir=args.out_dir,
        species=args.species,
        typ=args.typ,
        basename=args.basename,
        fmt=args.fmt,
        floor=args.floor,
        depth_var=args.depth_var,
        rivers_var=args.rivers_var,
        river_months=river_months,
        offshore_var=args.offshore_var,
        offshore_months=offshore_months,
        gauss_sigma=args.gauss_sigma,
        verbose=(not args.quiet),
    )



if __name__ == "__main__":
    sys.exit(main())