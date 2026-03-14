#!/usr/bin/env python3
"""
Data Prep – Herring Habitat Capacity (Weighted Sum)
---------------------------------------------------
Combine monthly herring layers with static habitat layers (depth, offshore, etc.)
to produce Ecospace-ready ASC rasters to use as habitat capacity layers.

Static layers are read from a NetCDF bundle (DEF_NC_PATH), typically produced by
ASC→NC converter. By default we use:
  - DEF_DEPTH_VAR            = "depth"                   # used to derive water mask (depth > 0) and as a weighted layer
  - DEF_OFFSHORE_VAR         = "offshore_areas_herring"  # 0/1 or 0..1, treated as a weighted layer
  - DEF_RIVERS_VAR           = None                      # optional additional layer (weighted if set)

can easily extend with more weighted NC variables via DEF_EXTRA_NC_VARS.

Outputs
-------
- {BASE}_{YYYY}_{MM}.asc     : monthly weighted-sum habitat capacity
- {BASE}_CLIM_{MM}.asc       : monthly climatology of the combined layer (across years)
- {BASE}_BASE_ALLYEARS.asc   : all-years sum of the combined layer

Notes
-----
- Water/include mask comes from NC depth > 0 (fallback to land ASC if needed).
- Fraser inland mask ASC (==1) is applied to zero out fake inland river cells.
- Each weighted layer is scaled to [0..1] by its max BEFORE weighting, unless noted.
- Final per-month map is normalized by DEF_OUTPUT_NORM_METHOD ("max" by default),
  then DEF_MIN_FLOOR is applied on non-zero water cells and the chosen normalization
  is re-enforced.

Dependencies
------------
- numpy, pandas, netCDF4
- (optional) helpers.saveASCFile if available on PYTHONPATH; otherwise a plain writer is used.
"""

from __future__ import annotations
from pathlib import Path
import re
import numpy as np
from netCDF4 import Dataset, num2date
from helpers import saveASCFile
import argparse
from scipy.ndimage import gaussian_filter

from helpers_ewe_asc_netcdf import (
    read_asc_header_and_dims,
    collect_herring_from_nc,
    depth_pref_from_depth,
    read_asc_array,
    parse_months_str,
    norm_array,
    scale_to_unit_by_max,
    gaussian_blur_masked,
    read_nc_var2d
)


# ===========================
# PARAMS (edit these)
# ===========================
# Inputs
DEF_IN_ASC_DIR          = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//SPAWN_DATA_ASC")
DEF_TEMPLATE_ASC_FP     = Path("..//data//forcing//ECOSPACE_in_climayr_1980-2018_PAR3_Sal4m_20240523//PAR-VarZ-VarK//PAR-VarZ-VarK_1980-2018.asc")

# NetCDF bundle (static layers)
DEF_NC_PATH             = Path("..//..//ecospace_seal_diets//Data//basemaps//ecospace_bundled_asc.nc")
DEF_DEPTH_VAR           = "depth"                    # 2D, meters; >0 => water
DEF_RIVERS_VAR          = None                       # e.g., "chum_rivers"; if set, treated as weighted layer
DEF_OFFSHORE_VAR        = "offshore_areas_herring"   # 0/1 (or 0..1) weight layer

# NEW: Herring from NC options
DEF_HERRING_FROM_NC     = True                       # True => read time series from NC instead of ASC files
DEF_HERRING_TS_VAR      = "herring_spawn_intensity"  # 3-D (time, y, x)
DEF_YEAR_COORD_VAR      = "year"                     # optional 1-D (time,)
DEF_MONTH_COORD_VAR     = "month"                    # optional 1-D (time,)
DEF_TIME_COORD_VAR      = "time"                     # fallback if year/month not present

# Optional ASC masks
# DEF_FRASER_MASK_ASC_FP  = Path("..//data//fraser_inland_mask.asc")
DEF_FRASER_MASK_VAR     = "fraser_inland_mask"
# DEF_LAND_ASC_FP         = Path("..//data//ecospacedepthgrid.asc")     # fallback include mask if NC depth missing (0=land)
DEF_LAND_MASK_VAR     = "depth"
DEF_SPAWN_MASK_VAR    = "no_spawn_shallows_herring"

# Herring file pattern (year-month only) [ASC fallback]
DEF_HERRING_PATTERN     = r"^HERRING_SPAWN_(\d{4})_(\d{2})\.asc$"

# Output
DEF_OUT_DIR             = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//HAB_CAP_ASC//")
DEF_BASE_OUTNAME        = "HERRING_HABCAP"
DEF_FMT                 = "%0.3f"

# Weights
DEF_W_HERRING           = 1.0
DEF_W_DEPTH             = 0.0
DEF_W_OFFSHORE          = 1.0

# Extra weighted NC variables (extend later)
# Each dict: {"var": "nc_var_name", "weight": 0.5, "apply_months": None}
DEF_EXTRA_NC_VARS       = []  # e.g., [{"var": "kelp_canopy", "weight": 0.7, "apply_months": [5,6,7,8]}]

# Herring transform before blending
DEF_HERRING_BINARY      = True
DEF_HERRING_BIN_THR     = 0.0
DEF_HERRING_NORM_METHOD = "sum"     # "sum" | "max" | "none"

# Final output normalization (after weighted sum & masks)
DEF_OUTPUT_NORM_METHOD  = "max"     # "sum" | "max" | "none"
                                    # - (sum across max and then normalise each cell to that val
                                    #   or max across map and normalise each cell rel to that val)
DEF_MIN_FLOOR           = 0.001

# === Applicability by month (1..12); None = all months ===
# adding some overlap so fish have time to move
DEF_OFFSHORE_APPLY_MONTHS = [4,5,6,7,8,9,10]         # e.g., [6,7,8,9] to apply only Jun–Sep
DEF_DEPTH_APPLY_MONTHS    = [10,11,12,1,2,3,4]         # e.g., [1,2,3,4,5] if depth matters only in winter/spring

# === Depth preference as a range (in meters, using NC 'depth') ===
# Modes:
#   "scale"   -> scale depth to [0..1] by max (deeper => larger)
#   "invert"  -> 1 - scale (shallower => larger)
#   "range_binary"  -> 1 inside [MIN, MAX], 0 outside
#   "range_trapezoid"-> 0 outside [MIN-SOFTEN, MAX+SOFTEN], ramps to 1 at [MIN..MAX]
DEF_DEPTH_MODE        = "range_trapezoid"   # or "range_binary", "scale", "invert"
DEF_DEPTH_MIN_M       = 0.5
DEF_DEPTH_MAX_M       = 10.0
DEF_DEPTH_SOFTEN_M    = 5.0                # only used for trapezoid; set 0 for hard edges

# Write monthly climatology & all-years base of the combined layer
DEF_WRITE_CLIM_AND_BASE = True

# Try to use helper writer (keeps orientation consistent); else plain writer
DEF_TRY_HELPER_WRITER   = True

# Gaussian blur sigma (cells). 0 disables.
DEF_GAUSS_SIGMA         = 2.0


# ===========================
# main
# ===========================
def main(argv=None):
    # ===== NEW: argparse (defaults from DEF_* above) =====
    ap = argparse.ArgumentParser(description="Weighted habitat capacity blend for herring (reads statics from NC; herring from NC or ASC).")
    ap.add_argument("--in_asc_dir", type=Path, default=DEF_IN_ASC_DIR)
    ap.add_argument("--template_asc", type=Path, default=DEF_TEMPLATE_ASC_FP)
    ap.add_argument("--nc_path", type=Path, default=DEF_NC_PATH)
    ap.add_argument("--depth_var", type=str, default=DEF_DEPTH_VAR)
    ap.add_argument("--rivers_var", type=str, default=DEF_RIVERS_VAR)
    ap.add_argument("--offshore_var", type=str, default=DEF_OFFSHORE_VAR)
    ap.add_argument("--fraser_mask_var", type=str, default=DEF_FRASER_MASK_VAR,
                    help="NC variable name of Fraser inland mask (1=mask out)")
    ap.add_argument("--land_mask_var", type=str, default=DEF_LAND_MASK_VAR,
                    help="NC variable name of land/water mask (0=land, nonzero=water). If not set, use depth>0")
    ap.add_argument("--herring_spawn_mask_var", type=str, default=DEF_SPAWN_MASK_VAR,
                    help="NC variable name of areas that are shallow where not much herring spawning is recorded")
    ap.add_argument("--herring_from_nc", action="store_true" if DEF_HERRING_FROM_NC else "store_false")
    ap.add_argument("--enable_herring_from_nc", dest="herring_from_nc", action="store_true")
    ap.add_argument("--disable_herring_from_nc", dest="herring_from_nc", action="store_false")
    ap.set_defaults(herring_from_nc=DEF_HERRING_FROM_NC)
    ap.add_argument("--herring_nc_var", type=str, default=DEF_HERRING_TS_VAR)
    ap.add_argument("--herring_year_var", type=str, default=DEF_YEAR_COORD_VAR)
    ap.add_argument("--herring_month_var", type=str, default=DEF_MONTH_COORD_VAR)
    ap.add_argument("--herring_time_var", type=str, default=DEF_TIME_COORD_VAR)

    # ap.add_argument("--fraser_mask", type=Path, default=DEF_FRASER_MASK_ASC_FP)
    # ap.add_argument("--land_asc", type=Path, default=DEF_LAND_ASC_FP)

    ap.add_argument("--out_dir", type=Path, default=DEF_OUT_DIR)
    ap.add_argument("--basename", type=str, default=DEF_BASE_OUTNAME)
    ap.add_argument("--fmt", type=str, default=DEF_FMT)

    ap.add_argument("--w_herring", type=float, default=DEF_W_HERRING)
    ap.add_argument("--w_depth",   type=float, default=DEF_W_DEPTH)
    ap.add_argument("--w_offshore",type=float, default=DEF_W_OFFSHORE)

    ap.add_argument("--extra_nc_vars", type=str, default="", help='Comma list of "var:weight" (e.g., "kelp:0.6,eelgrass:0.4")')

    ap.add_argument("--herring_binary", action="store_true" if DEF_HERRING_BINARY else "store_false")
    ap.add_argument("--enable_herring_binary", dest="herring_binary", action="store_true")
    ap.add_argument("--disable_herring_binary", dest="herring_binary", action="store_false")
    ap.set_defaults(herring_binary=DEF_HERRING_BINARY)
    ap.add_argument("--herring_bin_thr", type=float, default=DEF_HERRING_BIN_THR)
    ap.add_argument("--herring_norm_method", type=str, choices=["sum","max","none"], default=DEF_HERRING_NORM_METHOD)

    ap.add_argument("--output_norm_method", type=str, choices=["sum","max","none"], default=DEF_OUTPUT_NORM_METHOD)
    ap.add_argument("--min_floor", type=float, default=DEF_MIN_FLOOR)

    ap.add_argument("--offshore_months", type=str, default=",".join(map(str, DEF_OFFSHORE_APPLY_MONTHS)) if DEF_OFFSHORE_APPLY_MONTHS else "all")
    ap.add_argument("--depth_months",    type=str, default=",".join(map(str, DEF_DEPTH_APPLY_MONTHS)) if DEF_DEPTH_APPLY_MONTHS else "all")

    ap.add_argument("--depth_mode",  type=str, choices=["scale","invert","range_binary","range_trapezoid"], default=DEF_DEPTH_MODE)
    ap.add_argument("--depth_min",   type=float, default=DEF_DEPTH_MIN_M)
    ap.add_argument("--depth_max",   type=float, default=DEF_DEPTH_MAX_M)
    ap.add_argument("--depth_soften",type=float, default=DEF_DEPTH_SOFTEN_M)

    ap.add_argument("--write_clim_and_base", action="store_true" if DEF_WRITE_CLIM_AND_BASE else "store_false")
    ap.add_argument("--enable_clim_base", dest="write_clim_and_base", action="store_true")
    ap.add_argument("--disable_clim_base", dest="write_clim_and_base", action="store_false")
    ap.set_defaults(write_clim_and_base=DEF_WRITE_CLIM_AND_BASE)

    ap.add_argument("--gauss_sigma", type=float, default=DEF_GAUSS_SIGMA, help="Gaussian blur sigma (cells), 0=off")

    args = ap.parse_args(argv)
    # ===== end argparse =====

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Header/dims from example ASC
    header, nrows, ncols = read_asc_header_and_dims(args.template_asc)

    # Load NC statics
    depth_nc     = read_nc_var2d(args.nc_path, args.depth_var,     (nrows, ncols))
    offshore_nc  = read_nc_var2d(args.nc_path, args.offshore_var,  (nrows, ncols))
    rivers_nc    = read_nc_var2d(args.nc_path, args.rivers_var,    (nrows, ncols)) if args.rivers_var else None

    # land/water include mask
    landmask_nc = read_nc_var2d(args.nc_path, args.land_mask_var, (nrows, ncols)) if args.land_mask_var else None
    if landmask_nc is not None:
        # treat 0 as land, nonzero as water; also handle NaNs
        include = np.nan_to_num(landmask_nc) != 0
    elif depth_nc is not None:
        include = depth_nc > 0
    else:
        raise RuntimeError("No water mask available: set --land_mask_var in NC or provide --depth_var")

    # Fraser mask (1 => mask out)
    fraser_nc = read_nc_var2d(args.nc_path, args.fraser_mask_var, (nrows, ncols)) if args.fraser_mask_var else None
    fraser_block = (np.nan_to_num(fraser_nc) >= 1) if fraser_nc is not None else np.zeros((nrows, ncols), dtype=bool)

    herring_nospawn_nc = read_nc_var2d(args.nc_path, args.herring_spawn_mask_var, (nrows, ncols)) if args.herring_spawn_mask_var else None
    herring_nospawn_block = (np.nan_to_num(herring_nospawn_nc) >= 1) if herring_nospawn_nc is not None else np.zeros((nrows, ncols), dtype=bool)

    # CHECK
    # with Dataset(DEF_NC_PATH, "r") as ds:
    #     v = ds.variables["herring_spawn_intensity"]  # or args.herring_nc_var
    #     print("dtype:", v.dtype, "shape:", v.shape)
    #     print("_FillValue:", getattr(v, "_FillValue", None), "missing_value:", getattr(v, "missing_value", None))
    #     ds.set_auto_maskandscale(True)
    #     arr = v[:]
    #     print("MaskedArray?", isinstance(arr, np.ma.MaskedArray))
    #     if isinstance(arr, np.ma.MaskedArray):
    #         print("masked %:", (arr.mask.sum() / arr.size) * 100, "%")
    #         print("min/max (filled 0):", np.nanmin(arr.filled(0.0)), np.nanmax(arr.filled(0.0)))
    #     else:
    #         print("min/max:", np.nanmin(arr), np.nanmax(arr))


    # Static weighted layers (scaled 0..1 over water)
    if depth_nc is not None:
        depth_pref = depth_pref_from_depth(
            depth_m=depth_nc,
            include=include,
            mode=args.depth_mode,
            dmin=args.depth_min,
            dmax=args.depth_max,
            soften=args.depth_soften,
        )
    else:
        depth_pref = np.zeros((nrows, ncols), float)

    offshore = scale_to_unit_by_max(offshore_nc, include) if offshore_nc is not None else np.zeros((nrows, ncols), float)
    rivers   = scale_to_unit_by_max(rivers_nc, include) if rivers_nc is not None else None

    # Extra NC vars (simple parser: "var:weight,var2:weight2")
    extra_vars = []
    if args.extra_nc_vars.strip():
        for tok in args.extra_nc_vars.split(","):
            tok = tok.strip()
            if not tok:
                continue
            if ":" in tok:
                v, w = tok.split(":")
                v = v.strip(); w = float(w.strip())
            else:
                v, w = tok, 1.0
            arr = read_nc_var2d(args.nc_path, v, (nrows, ncols))
            if arr is None:
                continue
            arr = scale_to_unit_by_max(arr, include)
            extra_vars.append({"var": v, "array": arr, "weight": w, "apply_months": None})

    # Month gating
    off_months = parse_months_str(args.offshore_months)
    dep_months = parse_months_str(args.depth_months)

    # Collect herring inputs: from NC or from ASC files
    herring_items = None
    if args.herring_from_nc:
        try:
            herring_items = collect_herring_from_nc(
                args.nc_path, args.herring_nc_var, nrows, ncols,
                args.herring_year_var, args.herring_month_var, args.herring_time_var
            )
        except Exception as e:
            print(f"[WARN] herring_from_nc failed: {e}. Falling back to ASC files.")
            herring_items = None

    if herring_items is None:
        rx = re.compile(DEF_HERRING_PATTERN)
        files = []
        for fp in Path(args.in_asc_dir).glob("HERRING_SPAWN_*.asc"):
            m = rx.match(fp.name)
            if m:
                files.append((int(m.group(1)), int(m.group(2)), fp))
        files.sort()
        if not files:
            raise FileNotFoundError(f"No herring ASCs matching {DEF_HERRING_PATTERN} in {args.in_asc_dir}")
        # convert to items [(yy,mm, array)]
        herring_items = []
        for (yy, mm, fp) in files:
            H = np.loadtxt(fp, skiprows=6)
            herring_items.append((yy, mm, H))

    # Accumulators for climatology/base of combined layer
    clim_accum = {m: np.zeros((nrows, ncols), dtype=float) for m in range(1, 13)}
    base_accum = np.zeros((nrows, ncols), dtype=float)

    # Process month-by-month
    for (yy, mm, H) in herring_items:
        if args.herring_binary:
            H = (H > args.herring_bin_thr).astype(float)

        H = norm_array(H, include, args.herring_norm_method, min_floor=0.0)

        use_off = (off_months is None) or (mm in off_months)
        use_dep = (dep_months is None) or (mm in dep_months)

        off_contrib = (args.w_offshore * offshore) if use_off else 0.0
        dep_contrib = (args.w_depth    * depth_pref) if use_dep else 0.0

        combined = args.w_herring * H + dep_contrib + off_contrib

        if rivers is not None:
            combined += rivers

        for ex in extra_vars:
            # (extend later with per-layer months if   like)
            combined += ex["weight"] * ex["array"]

        # Masks
        combined[~include] = 0.0
        combined[fraser_block] = 0.0
        combined[herring_nospawn_block] = 0.0

        # ===== NEW: optional Gaussian blur BEFORE the final normalization =====
        if args.gauss_sigma and args.gauss_sigma > 0:
            combined = gaussian_blur_masked(combined, include, args.gauss_sigma)

        # Final norm + floor
        combined = norm_array(combined, include, args.output_norm_method, min_floor=args.min_floor)

        # Save
        out_fp = args.out_dir / f"{args.basename}_{yy}_{mm:02d}.asc"
        saveASCFile(
            str(out_fp),
            np.flipud(combined),
            bottomleft_row_ewe=0,
            upperleft_row_ewe=nrows,
            upperleft_col_ewe=0,
            sigdigfmt=args.fmt,
            ASCheader=header,
            dfPlumeMask=None,
            dfLandMask=None
        )
        s = combined[include].sum()
        print(f"[INFO] Wrote {out_fp} (sum over water={s:.6f})")

        clim_accum[mm] += combined
        base_accum     += combined

    if args.write_clim_and_base:
        for mm in range(1, 13):
            clim = norm_array(clim_accum[mm], include, args.output_norm_method, min_floor=args.min_floor)
            out_name = args.out_dir / f"{args.basename}_CLIM_{mm:02d}.asc"
            saveASCFile(
                str(out_name),
                np.flipud(clim),
                bottomleft_row_ewe=0,
                upperleft_row_ewe=nrows,
                upperleft_col_ewe=0,
                sigdigfmt=args.fmt,
                ASCheader=header,
                dfPlumeMask=None,
                dfLandMask=None
            )
            s = clim[include].sum()
            print(f"[INFO] Wrote {out_name} (sum over water={s:.6f})")

        base = norm_array(base_accum, include, args.output_norm_method, min_floor=args.min_floor)
        out_name = f"{args.basename}_BASE_ALLYEARS.asc"
        saveASCFile(
            str(args.out_dir / out_name),
            np.flipud(base),
            bottomleft_row_ewe=0,
            upperleft_row_ewe=nrows,
            upperleft_col_ewe=0,
            sigdigfmt=args.fmt,
            ASCheader=header,
            dfPlumeMask=None,
            dfLandMask=None
        )
        s = base[include].sum()
        print(f"[INFO] Wrote {out_name} (sum over water={s:.6f})")

if __name__ == "__main__":
    main()



#
#
#
#
# from __future__ import annotations
# import argparse
# from pathlib import Path
# import sys
# import numpy as np
# import pandas as pd
# from netCDF4 import Dataset
# from scipy.ndimage import gaussian_filter
#
# # Import helpers from existing module
# # (helpers.py needs to be reachable – same folder or PYTHONPATH)
# from helpers import (
#     saveASCFile, getDataFrame
# )
#
# from helpers_ewe_asc_netcdf import (
#     parse_depth_bands,
#     read_asc_header_and_dims,
#     build_inclusion_mask,
#     month_iter,
#     resolve_value_column,
#     parse_year_month,
#     clean_numeric_column,
#     load_masks_from_nc,
#     parse_months,
#     compose_monthly_weights,
#     make_depth_weight,
#     build_output_name,
#     smooth_gaussian_masked,
#     apply_min_floor,
#     _scale_to_unit_interval
# )
#
# # ------------------------------
# # paths and default params
# # ------------------------------
# DEFAULT_SPAWN_P    = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//")
# DEFAULT_SPAWN_FP  = Path(DEFAULT_SPAWN_P, "Pacific_SoG_herring_spawn_index_data_2024_EN_taggedEwERowsCols.csv")
# DEF_OUT_P = Path(DEFAULT_SPAWN_P, "ECOSPACE_ASC")
# DEF_TEMPLATE_ASC_P = Path("..//data//forcing//ECOSPACE_in_climayr_1980-2018_PAR3_Sal4m_20240523//PAR-VarZ-VarK//")
# DEF_TEMPLATE_ASC_FP = Path(DEF_TEMPLATE_ASC_P, "PAR-VarZ-VarK_1980-2018.asc")
# # DEF_PLUME_P = Path("..//data//")
# # DEF_PLUME_FP = Path(DEF_PLUME_P, "Fraser_Plume_Region.asc")
# DEF_PLUME_FP = None
# DEF_LAND_P = "..//data//"
#
# DEF_LAND_F = "ecospacedepthgrid.asc"
# DEF_LAND_PF = Path("..//data//", DEF_LAND_F)
#
# DEF_MONTHLY_ACROSS_ALL = True  # determines whether it's monthly average hab cap based on all years in spawn data
#
# DEF_SPWN_COL = "RelativeEggB"
# DEF_DATE_COL = "StartDate"
# DEF_YEAR = None # optional to restrict to just one year (default to None)
# DEF_FMT = "%0.2f"
#
# DEF_INDEX_BASE_EWE = 1 # EwE is 1-indexed
# DEF_BASEOUTNAME = "HERRING_SPAWN"
# DEF_LOG_NORM_BASEASC = False
# DEF_BINARYMAPS_TF = True
# DEF_MIN_FILL = 0.05
#
# --- New defaults for NC masks + depth preference + smoothing ---
# DEF_NC_PATH = "..//..//ecospace_seal_diets//Data//basemaps//ecospace_bundled_asc.nc" # multi ASC stored in this NC
# DEF_DEPTH_VAR = "depth"
# DEF_RIVERS_VAR = None                                # e.g., "chum_rivers" (or exact var)
# DEF_OFFSHORE_VAR = "offshore_areas_herring"                              # e.g., "offshore_areas"
#
# # Month gating as CLI-friendly strings; parsed later via parse_months()
# DEF_RIVERS_MONTHS_STR = "none"                       # "none" | "all" | "9-12" | "1,2,3,6-8"
# DEF_OFFSHORE_MONTHS_STR = "6-10"
#
# DEF_RIVERS_WEIGHT = 1.5                              # not used
# DEF_OFFSHORE_WEIGHT = 1.2
# DEF_COMBINE_MODE = "multiply"                        # "multiply" | "overwrite"
# DEF_GAUSS_SIGMA = 1.0
#
# # Keep one place of truth for floors:
# DEF_MIN_FLOOR = DEF_MIN_FILL                         # use same floor for monthly + base
#
# # Depth preference
# DEF_DEPTH_PREF_ON = False
# DEF_DEPTH_BANDS_STR = "1-10"                           # e.g., "0-15,25-45"
# DEF_DEPTH_PEAK_W = 5.0
# DEF_DEPTH_FLOOR_W = 1.0
#
#
#
#
# # ------------------------------
# # Presence-Probability mode (weighted + smoothed)
# # ------------------------------
# import numpy as np
# import pandas as pd
# from pathlib import Path
#
# def write_monthly_presence_probability(
#     df: pd.DataFrame,
#     example_asc: Path,
#     out_dir: Path,
#     output_basename: str,
#     index_base: int,
#     date_col: str | None,
#     plume_mask_path: Path | None,
#     land_mask_path: Path | None,
#     years_denominator: str | int,
#     min_fill: float,
#     fmt: str,
#     verbose: bool,
#     # --- NEW: equal-weight factor blend controls ---
#     combine_strategy: str = "mean",        # "mean" (default) or "geomean"
#     include_depth: bool = True,
#     include_offshore: bool = True,
#     # thresholding for *data* presence (defines the data-driven factor)
#     presence_threshold: float = 0.0,
#     presence_quantile: float | None = None,
#     # NC/masks
#     nc_path: Path | None = None,
#     depth_var: str = "depth",
#     rivers_var: str | None = None,         # not used as a separate factor here can add later)
#     offshore_var: str | None = None,
#     offshore_months: set[int] = frozenset(),
#     # depth preference (factor)
#     depth_pref_on: bool = False,
#     depth_bands: list[tuple[float, float]] | None = None,
#     # If used weights before, we now map depth to 0–1 directly, ignoring weights
#     # smoothing
#     gauss_sigma: float = 0.0,
# ):
#     """
#     Build 12 presence-probability maps by *equally blending* multiple 0–1 factors:
#
#       combined[m] = mean( factors[m] )
#
#     Factors (0–1):
#       1) Data presence frequency (fraction of years with presence in month m)
#       2) Depth suitability (0–1 membership for configured bands) [optional]
#       3) Offshore suitability (0/1, month-gated) [optional]
#
#     No renormalization (they are probabilities/suitabilities). Apply a min floor at the end.
#     """
#     out_dir = Path(out_dir)
#     out_dir.mkdir(parents=True, exist_ok=True)
#
#     ASCheader, nrows, ncols = read_asc_header_and_dims(example_asc)
#
#     # ASC masks -> inclusion
#     dfPlumeMask = None
#     if plume_mask_path is not None and Path(plume_mask_path).exists():
#         dfPlumeMask = getDataFrame(str(plume_mask_path), "-9999.00000000") == 1
#
#     dfLandMask = None
#     if land_mask_path is not None and Path(land_mask_path).exists():
#         dfLandMask = getDataFrame(str(land_mask_path), "-9999.00000000") == 0
#
#     inclusion = build_inclusion_mask(nrows, ncols, dfPlumeMask, dfLandMask)
#     if inclusion.shape != (nrows, ncols):
#         raise ValueError("inclusion mask shape mismatch")
#
#     # Optional NC bundle
#     depth_2d = None
#     water_mask = None
#     offshore_mask = None
#     if nc_path is not None:
#         depth_2d, water_mask, _rivers_mask, offshore_mask = load_masks_from_nc(
#             Path(nc_path),
#             depth_var=depth_var,
#             rivers_var=rivers_var,
#             offshore_var=offshore_var
#         )
#         if depth_2d.shape != (nrows, ncols):
#             raise ValueError(f"NC depth shape {depth_2d.shape} != ASC {(nrows, ncols)}")
#
#     if water_mask is None:
#         water_mask = inclusion.copy()
#
#     # -------- Factor 2: Depth suitability (0–1) --------
#     # Make a *normalized membership* map (0 outside bands, 1 at band peaks).
#     depth_suit = None
#     if include_depth and depth_pref_on and depth_2d is not None and (depth_bands or []):
#         # Build membership as max across bands with triangular membership per band.
#         depth_suit = np.zeros((nrows, ncols), dtype=float)
#         for (z0, z1) in (depth_bands or []):
#             if z1 <= z0:
#                 continue
#             mid = 0.5 * (z0 + z1)
#             halfw = max(0.5 * (z1 - z0), 1e-6)
#             tri = np.clip(1.0 - np.abs(depth_2d - mid) / halfw, 0.0, 1.0)  # 0..1
#             depth_suit = np.maximum(depth_suit, tri)
#         # zero on land; already 0..1
#         depth_suit = np.where(water_mask, depth_suit, 0.0)
#         if verbose:
#             mn, mx = depth_suit[inclusion].min(initial=0.0), depth_suit[inclusion].max(initial=0.0)
#             print(f"[PPROB] Depth suitability range (on included): {mn:.3g}..{mx:.3g}")
#     # else: depth_suit stays None (factor excluded)
#
#     # -------- Factor 3: Offshore suitability (0/1, month-gated) --------
#     # Build per-month maps later because of gating.
#
#     # Distinct years
#     all_years = sorted([int(y) for y in df["Year"].dropna().unique()])
#
#     # Helper to compute data presence frequency for one month (factor 1)
#     def data_presence_frequency_for_month(mm: int) -> np.ndarray:
#         df_m = df[df["Month"] == mm]
#         years_with_month = sorted([int(y) for y in df_m["Year"].dropna().unique()])
#         bin_slices = []
#         for yy in years_with_month:
#             sub = df_m[df_m["Year"] == yy]
#             grid = np.zeros((nrows, ncols), dtype=float)
#             if not sub.empty:
#                 # choose value column
#                 val_candidates = [c for c in ["_val", "intensity", "value", "spawn_intensity", "spawn_value"] if c in sub.columns]
#                 if not val_candidates:
#                     raise KeyError("No value column found for data factor")
#                 vcol = val_candidates[0]
#                 rows = (sub["ecospace_row"].to_numpy() - index_base).astype(int)
#                 cols = (sub["ecospace_col"].to_numpy() - index_base).astype(int)
#                 vals = sub[vcol].to_numpy(float)
#                 ok = (rows >= 0) & (rows < nrows) & (cols >= 0) & (cols < ncols)
#                 rows, cols, vals = rows[ok], cols[ok], vals[ok]
#                 np.add.at(grid, (rows, cols), vals)
#
#             # Apply inclusion and decide presence with threshold/quantile
#             grid[~inclusion] = 0.0
#             thr = presence_threshold
#             if presence_quantile is not None:
#                 pos = grid[inclusion & (grid > 0)]
#                 if pos.size:
#                     thr = float(np.quantile(pos, presence_quantile))
#             bin_slices.append((grid > thr).astype(float))
#
#         # Denominator
#         if isinstance(years_denominator, int):
#             D = int(years_denominator)
#         elif isinstance(years_denominator, str):
#             key = years_denominator.lower()
#             if key == "monthly":
#                 D = len(years_with_month)
#             elif key == "global":
#                 D = len(all_years)
#             else:
#                 raise ValueError('years_denominator must be "monthly", "global" or an integer')
#         else:
#             raise ValueError('years_denominator must be "monthly", "global" or an integer')
#
#         if not bin_slices or D == 0:
#             return np.zeros((nrows, ncols), dtype=float)
#         count_present = np.sum(bin_slices, axis=0)
#         P = count_present / float(D)  # 0..1
#         P[~inclusion] = 0.0
#         return P
#
#     for mm in month_iter():
#         # Factor 1: data-driven presence frequency
#         P_data = data_presence_frequency_for_month(mm)
#
#         # Factor 2: depth suitability (same for all months if provided)
#         F_depth = depth_suit if depth_suit is not None else None
#
#         # Factor 3: offshore suitability (month-gated)
#         F_off = None
#         if include_offshore and offshore_mask is not None:
#             if (mm in offshore_months) if offshore_months else False:
#                 # active ⇒ 1 on offshore cells, 0 elsewhere
#                 F_off = np.where(offshore_mask & water_mask, 1.0, 0.0).astype(float)
#             else:
#                 # inactive ⇒ neutral 0.5 everywhere on water
#                 F_off = np.where(water_mask, 0.5, 0.0).astype(float)
#
#         # Collect active factors
#         factors = [P_data]
#         if F_depth is not None:
#             factors.append(F_depth)
#         if F_off is not None:
#             factors.append(F_off)
#
#         # Combine (equal weights)
#         if combine_strategy == "geomean":
#             # geometric mean (avoid zeros by adding tiny epsilon on included cells)
#             eps = 1e-12
#             prod = np.ones((nrows, ncols), dtype=float)
#             for f in factors:
#                 prod *= np.where(inclusion, np.clip(f, 0.0, 1.0) + eps, 0.0)
#             combined = np.where(inclusion, prod ** (1.0 / len(factors)) - eps, 0.0)
#         else:
#             # arithmetic mean
#             stacked = np.stack([np.clip(f, 0.0, 1.0) for f in factors], axis=0)
#             combined = stacked.mean(axis=0)
#             combined[~inclusion] = 0.0
#
#         # Optional smoothing (edge-aware)
#         if gauss_sigma and gauss_sigma > 0:
#             edge_mask = water_mask.astype(float) if water_mask is not None else inclusion.astype(float)
#             combined = smooth_gaussian_masked(combined, edge_mask, gauss_sigma)
#
#         # Apply floor (no renormalization)
#         if min_fill > 0:
#             combined = np.clip(combined, 0.0, 1.0)
#             combined = apply_min_floor(combined, inclusion, min_fill, renormalize=False)
#             combined = np.minimum(combined, 1.0)
#
#         # Write map
#         out_name = f"{output_basename}_PPROB_{mm:02d}.asc"
#         out_path = out_dir / out_name
#         if verbose:
#             facs = ["data"] + (["depth"] if F_depth is not None else []) + (["offshore"] if F_off is not None else [])
#             print(f"[PPROB] Month {mm:02d}: factors={facs}, combine={combine_strategy} -> {out_name}")
#         saveASCFile(
#             str(out_path),
#             np.flipud(combined),
#             bottomleft_row_ewe=0,
#             upperleft_row_ewe=nrows,
#             upperleft_col_ewe=0,
#             sigdigfmt=fmt,
#             ASCheader=ASCheader,
#             dfPlumeMask=dfPlumeMask,
#             dfLandMask=dfLandMask
#         )
#
#
#
# # ------------------------------
# # Core conversion
# # ------------------------------
#
#
# def csv_to_monthly_ascs(
#     spawn_csv: Path,
#     out_dir: Path,
#     example_asc: Path,
#     value_col: str | None = None,
#     date_col: str | None = "StartDate",
#     year: int | None = None,
#     fmt: str = "%0.1f",
#     plume_mask_path: Path | None = None,
#     land_mask_path: Path | None = None,
#     index_base: int = 1,
#     output_basename: str = "HERRING_SPAWN",
#     log_norm: bool = False,
#     log_eps: float = 1e-6,
#     binary: bool = False,
#     binary_threshold: float = 0.0,
#     min_fill: float = 0.0,     # <- keep existing name; used as the 'floor'
#     presence_prob: bool = False,
#     years_denominator: str | int = "monthly",
#     verbose: bool = True,
#     # --- NEW: NC mask bundle & weighting controls ---
#     nc_path: Path | None = None,
#     depth_var: str = "depth",
#     rivers_var: str | None = None,
#     offshore_var: str | None = None,
#     rivers_months: set[int] = frozenset(),
#     offshore_months: set[int] = frozenset(),
#     rivers_weight: float = 1.5,
#     offshore_weight: float = 1.2,
#     combine_mode: str = "multiply",  # "multiply" (recommended) or "overwrite"
#     gauss_sigma: float = 0.0,
#     # --- NEW: depth preference controls ---
#     depth_pref_on: bool = False,
#     depth_bands: list[tuple[float, float]] | None = None,
#     depth_peak_w: float = 2.0,
#     depth_floor_w: float = 1.0,
# ) -> None:
#     """
#     Convert a tagged herring spawn CSV into 12 monthly ASC grids, and also writes a single **all-years base layer** that sums across all years/months and normalizes to sum=1 over included cells.
#
#     Normalization step (per-month):
#     --------------------------------
#     After accumulating cell values, we apply masks (land/plume) in-memory
#     and normalize so that the *sum over included cells* equals 1.0.
#     If the monthly sum is 0, we leave the grid at zeros (no normalization).Optionally apply
#     NC-provided mask layers (rivers/offshore), month-gated, plus a depth-preference
#     weighting and optional Gaussian smoothing before normalization. Also writes an
#     all-years base layer (summed across months/years then normalized to sum=1
#     over included cells)
#
#     Parameters
#     ----------
#     spawn_csv : Path
#         CSV with at least: ecospace_row, ecospace_col, intensity column,
#         and either (Year, Month) or a date column (e.g., StartDate).
#     out_dir : Path
#         Where to write the ASC files.
#     example_asc : Path
#         Example ASC used to read header & grid dims (nrows/ncols).
#     value_col : str, optional
#         Column in the CSV to aggregate onto the grid.
#         If None, we try to auto-detect common names.
#     date_col : str, optional
#         If Year/Month absent, parse these from this date column.
#     year : int, optional
#         If provided, filter to this calendar year; otherwise use the year(s)
#         present in the CSV (writing a set of 12 files for each year found).
#     fmt : str
#         Numeric format string for ASC writing (e.g., "%0.1f", "%0.2f").
#     plume_mask_path : Path, optional
#         ASC file path for a plume mask (1 inside; masked to 0). Optional.
#     land_mask_path : Path, optional
#         ASC file path for a land mask (0=land). Optional.
#     index_base : int
#         1 if ecospace_row/col are 1-based; 0 if zero-based.
#     output_basename : str
#         Prefix for output file names.
#     verbose : bool
#         Print progress.
#     """
#     if verbose:
#         print(f"[INFO] Reading: {spawn_csv}")
#     df = pd.read_csv(spawn_csv)
#
#     # Guard: required row/col fields
#     required = ["ecospace_row", "ecospace_col"]
#     missing = [c for c in required if c not in df.columns]
#     if missing:
#         raise KeyError(f"CSV missing required columns: {missing}")
#
#     # Resolve value and year/month
#     val_col = resolve_value_column(df, value_col)
#     df2 = parse_year_month(df, date_col=date_col)
#
#     # Clean/convert value column to numeric once up-front
#     df2["_val"] = clean_numeric_column(df2, val_col, verbose=verbose)
#
#     if presence_prob:
#         write_monthly_presence_probability(
#             df=df2,
#             example_asc=example_asc,
#             out_dir=out_dir,
#             output_basename=output_basename,
#             index_base=index_base,
#             date_col=date_col,
#             plume_mask_path=plume_mask_path,
#             land_mask_path=land_mask_path,
#             years_denominator=years_denominator,
#             min_fill=min_fill,
#             fmt=fmt,
#             verbose=verbose,
#             # NEW equal-weight blend:
#             combine_strategy="mean",  # or "geomean"
#             include_depth=depth_pref_on,  # tie to toggle
#             include_offshore=(offshore_var is not None),
#             presence_threshold=binary_threshold,  # keep CLI threshold
#             presence_quantile=None,  # or set to e.g. 0.10
#             nc_path=nc_path,
#             depth_var=depth_var,
#             rivers_var=rivers_var,
#             offshore_var=offshore_var,
#             offshore_months=offshore_months,
#             depth_pref_on=depth_pref_on,
#             depth_bands=depth_bands,
#             gauss_sigma=gauss_sigma,
#         )
#         if verbose:
#             print("[INFO] Done (presence-probability equal-blend mode).")
#         return
#
#     # Filter by year if requested
#     if year is not None:
#         df2 = df2[df2["Year"] == year]
#         if df2.empty:
#             raise ValueError(f"No records found for --year {year} in {spawn_csv}")
#
#     # Load ASC header/dims
#     header, nrows, ncols = read_asc_header_and_dims(example_asc)
#     if verbose:
#         print(f"[INFO] Grid dims from example ASC: nrows={nrows}, ncols={ncols}")
#
#     # Optional masks (as DataFrames) for saveASCFile and normalization
#     dfPlumeMask = None
#     if plume_mask_path is not None and Path(plume_mask_path).exists():
#         dfPlumeMask = getDataFrame(str(plume_mask_path), "-9999.00000000") == 1
#
#     dfLandMask = None
#     if land_mask_path is not None and Path(land_mask_path).exists():
#         # land mask: (depth==0) in many configs; we mask to 0.0
#         dfLandMask = getDataFrame(str(land_mask_path), "-9999.00000000") == 0
#
#     # Prepare inclusion mask for normalization (True = include in sum)
#     inclusion = build_inclusion_mask(nrows, ncols, dfPlumeMask, dfLandMask)
#
#     # --- NC masks: depth/water/rivers/offshore ---
#     depth_2d = None
#     water_mask = None
#     rivers_mask = None
#     offshore_mask = None
#
#     if nc_path is not None:
#         depth_2d, water_mask, rivers_mask, offshore_mask = load_masks_from_nc(
#             Path(nc_path),
#             depth_var=depth_var,
#             rivers_var=rivers_var,
#             offshore_var=offshore_var
#         )
#         # shape checks vs ASC dims
#         if depth_2d.shape != (nrows, ncols):
#             raise ValueError(f"NC depth shape {depth_2d.shape} != ASC dims {(nrows, ncols)}")
#
#     # If still no water_mask, derive it from inclusion (or land mask)
#     if water_mask is None:
#         # inclusion “cells that count”; this keeps land out
#         water_mask = inclusion.copy()
#
#     # 0..1 depth membership (max across bands); None if disabled or missing depth
#     depth_suit = None
#     if depth_pref_on and depth_2d is not None and (depth_bands or []):
#         # reuse triangular helper: peak=1, floor=0 → a proper membership
#         depth_suit = make_depth_weight(depth_2d,
#                                        bands=depth_bands,
#                                        peak_weight=1.0,
#                                        floor_weight=0.0)
#         depth_suit = np.where(water_mask, depth_suit, 0.0)
#
#     # Accumulator for all-years base map (post-mask, pre-normalization)
#     base_accum = np.zeros((nrows, ncols), dtype=float)
#     # Ensure output dir
#     out_dir = Path(out_dir)
#     out_dir.mkdir(parents=True, exist_ok=True)
#
#     # Make files for each year present (or the specified year)
#     years = sorted([int(y) for y in df2["Year"].dropna().unique()])
#     if year is not None and year not in years:
#         years = [year]  # should not happen due to filter above
#
#     for yy in years:
#
#         df_y = df2[df2["Year"] == yy].copy()
#         for mm in month_iter():
#
#             df_m = df_y[df_y["Month"] == mm]
#
#             # --- (A) Build raw data grid for this (year,month) ---
#             raw = np.zeros((nrows, ncols), dtype=float)
#             if not df_m.empty:
#                 rows = (df_m["ecospace_row"].to_numpy() - index_base).astype(int)
#                 cols = (df_m["ecospace_col"].to_numpy() - index_base).astype(int)
#                 vals = df_m["_val"].to_numpy(float)
#                 ok = (rows >= 0) & (rows < nrows) & (cols >= 0) & (cols < ncols)
#                 rows, cols, vals = rows[ok], cols[ok], vals[ok]
#                 np.add.at(raw, (rows, cols), vals)
#             # respect inclusion (land/plume) early
#             raw[~inclusion] = 0.0
#
#             # --- (B) FACTOR 1: data suitability 0..1 (optionally log-compressed) ---
#             data_suit = _scale_to_unit_interval(
#                 raw, inclusion,
#                 log_norm=log_norm, log_eps=log_eps,
#                 quantile=None  # or set e.g. 0.99 if outliers dominate
#             )
#
#             # --- (C) FACTOR 2: depth suitability (0..1), precomputed ---
#             F_depth = depth_suit  # may be None
#
#             # --- (D) FACTOR 3: offshore suitability (month-gated) ---
#             F_off = None
#             if offshore_mask is not None:
#                 # active months: 1 on offshore cells; inactive: neutral 0.5 on water
#                 off_active = (mm in offshore_months) if offshore_months else False
#                 if off_active:
#                     F_off = np.where(offshore_mask & water_mask, 1.0, 0.0).astype(float)
#                 else:
#                     F_off = np.where(water_mask, 0.5, 0.0).astype(float)
#
#             # --- (E) Combine factors equally (arithmetic mean by default) ---
#             factors = [data_suit]
#             if F_depth is not None:
#                 factors.append(F_depth)
#             if F_off is not None:
#                 factors.append(F_off)
#
#             combined = np.stack([np.clip(f, 0.0, 1.0) for f in factors], axis=0).mean(axis=0)
#             combined[~inclusion] = 0.0
#
#             # --- (F) Optional Gaussian smoothing (edge-aware) ---
#             if gauss_sigma and gauss_sigma > 0:
#                 edge_mask = water_mask.astype(float) if water_mask is not None else inclusion.astype(float)
#                 combined = smooth_gaussian_masked(combined, edge_mask, gauss_sigma)
#
#             # Keep a pre-norm copy for the BASE layer accumulation
#             for_base = np.where(inclusion, combined, 0.0)
#             base_accum += for_base
#
#             # --- (G) Output: binary vs normalized density ---
#             if binary:
#                 # threshold to binary, then (optionally) floor without renormalizing
#                 out_grid = (combined > binary_threshold).astype(float)
#                 if min_fill and min_fill > 0:
#                     out_grid = apply_min_floor(out_grid, inclusion, min_fill, renormalize=False)
#                     out_grid = np.minimum(out_grid, 1.0)
#                 if verbose:
#                     print(f"[INFO] {yy}-{mm:02d}: factor-mean → binary (thr={binary_threshold}).")
#             else:
#                 # sum-to-1 density while guaranteeing the floor
#                 out_grid = np.clip(combined, 0.0, 1.0)
#                 out_grid = apply_min_floor(out_grid, inclusion, min_fill,
#                                            renormalize=True, target_sum=1.0)
#                 if verbose:
#                     s = out_grid[inclusion].sum()
#                     print(f"[INFO] {yy}-{mm:02d}: factor-mean → normalized (sum={s:.6f}).")
#
#             # --- (H) Write ASC ---
#             out_name = build_output_name(output_basename, yy, mm)
#             out_path = out_dir / out_name
#             if verbose:
#                 print(f"[INFO] Writing: {out_path}")
#             saveASCFile(
#                 str(out_path),
#                 np.flipud(out_grid),
#                 bottomleft_row_ewe=0,
#                 upperleft_row_ewe=nrows,
#                 upperleft_col_ewe=0,
#                 sigdigfmt=fmt,
#                 ASCheader=header,
#                 dfPlumeMask=dfPlumeMask,
#                 dfLandMask=dfLandMask
#             )
#
#     # ---- Write all-years base layer -----------------------------------------
#     if binary:
#         base_grid = (base_accum > binary_threshold).astype(float)
#         if verbose:
#             print(f"[INFO] BASE: binary presence-any map (threshold {binary_threshold}).")
#     else:
#         base_work = base_accum.copy()
#         # optional log-normalization for base layer as well
#         if log_norm:
#             pos = base_work > 0
#             if pos.any():
#                 base_work[pos] = np.log(base_work[pos] + log_eps)
#                 base_work[~pos] = 0.0
#         base_total = base_work.sum()
#         if base_total > 0:
#             base_grid = base_work / base_total
#             if verbose:
#                 print(f"[INFO] BASE: normalized all-years sum to 1.0 (pre-norm sum was {base_total:.6g}).")
#         else:
#             base_grid = base_work
#             if verbose:
#                 print("[INFO] BASE: total is zero; writing all-zero base map.")
#
#     if min_fill and min_fill > 0:
#         base_grid[inclusion] = np.maximum(base_grid[inclusion], min_fill)
#         base_grid = np.clip(base_grid, 0.0, 1.0)
#
#     base_out_name = f"{output_basename}_BASE.asc"
#     base_out_path = out_dir / base_out_name
#     if verbose:
#         print(f"[INFO] Writing base layer: {base_out_path}")
#     grid_to_write = np.flipud(base_grid)  # counteracts saveASCFile’s internal reverse
#     saveASCFile(str(base_out_path),
#                 grid_to_write,
#                 bottomleft_row_ewe=0,
#                 upperleft_row_ewe=nrows,
#                 upperleft_col_ewe=0,
#                 sigdigfmt=fmt,
#                 ASCheader=header,
#                 dfPlumeMask=dfPlumeMask,
#                 dfLandMask=dfLandMask)
#     # -------------------------------------------------------------------------
#
#     if verbose:
#         print("[INFO] Done.")
#
#
# # ------------------------------
# # CLI
# # ------------------------------
#
# def main(argv=None):
#     p = argparse.ArgumentParser(description="Convert tagged herring spawn CSV into monthly Ecospace ASC rasters (per-month normalized to sum=1 over included cells).")
#     p.add_argument("--spawn_csv", type=Path, default=DEFAULT_SPAWN_FP, help="Path to tagged CSV with ecospace_row/col and intensity")
#     p.add_argument("--out_dir", type=Path, default=DEF_OUT_P, help="Directory to write ASC files")
#     p.add_argument("--example_asc", type=Path, default=DEF_TEMPLATE_ASC_FP, help="ASC to copy header/dimensions from (e.g., an existing Ecospace input)")
#     p.add_argument("--value_col", type=str, default=DEF_SPWN_COL, help="Column holding the intensity to map (default: RelativeEggB)")
#     p.add_argument("--date_col", type=str, default=DEF_DATE_COL, help="Date column to derive Month/Year if needed (default: StartDate)")
#     p.add_argument("--year", type=int, default=DEF_YEAR,
#                    help="Restrict to a single year (default: None)")
#     p.add_argument("--fmt", type=str, default=DEF_FMT, help="Numeric format string for ASC values (default: %0.1f)")
#     p.add_argument("--plume_mask", type=Path, default=DEF_PLUME_FP, help="Optional: plume mask ASC (1=plume region)")
#     p.add_argument("--land_mask", type=Path, default=DEF_LAND_PF, help="Optional: land/depth ASC to mask (0=land)")
#     p.add_argument("--index_base", type=int, default=DEF_INDEX_BASE_EWE, choices=[0,1], help="Indexing base for ecospace_row/col (default: 1)")
#     p.add_argument("--basename", type=str, default=DEF_BASEOUTNAME, help="Output file basename/prefix")
#     p.add_argument("--log_norm", action="store_true",
#                    help="Apply log transform before normalization (compress extremes)")
#     p.add_argument("--log_eps", type=float, default=1e-6, help="Epsilon added inside log(v+eps)")
#
#     # Binary toggle (default from DEF_BINARYMAPS_TF)
#     p.add_argument("--binary", dest="binary", action="store_true",
#                    help="Use presence/absence instead of intensity (v>threshold -> 1)")
#
#     p.add_argument("--no-binary", dest="binary", action="store_false",
#                    help="Disable binary mode; use intensity weights.")
#
#     p.set_defaults(binary=DEF_BINARYMAPS_TF)
#     p.add_argument("--binary_threshold", type=float, default=0.0, help="Threshold for presence (default: > 0)")
#
#     p.add_argument("--min_fill", type=float, default=DEF_MIN_FILL, help="Floor for values in map")
#
#     # 'Presence-probability' toggle (default from DEF_MONTHLY_ACROSS_ALL)
#     # poor name - means the output is just twelve maps with most popular spawning areas across all years
#     p.add_argument("--presence_prob", dest="presence_prob", action="store_true",
#                    help="Write 12 monthly probability-of-presence maps across all years (bypasses per-year outputs).")
#     p.add_argument("--no-presence_prob", dest="presence_prob", action="store_false",
#                    help="Disable presence-probability mode; generate per-year/month maps.")
#     p.set_defaults(presence_prob=DEF_MONTHLY_ACROSS_ALL)
#
#     p.add_argument("--years_denominator", type=str, default="monthly",
#                    help='Denominator for probability: "monthly" (default), "global", or an integer (e.g., 40).')
#     p.add_argument("--quiet", action="store_true", help="Suppress progress messages")
#     g = p.add_argument_group("Mask & habitat weighting (optional)")
#
#     g.add_argument("--nc_path", type=str, default=DEF_NC_PATH,
#                    help="Path to NC bundle that contains static 2D layers (depth, rivers, offshore).")
#
#     g.add_argument("--depth_var", type=str, default=DEF_DEPTH_VAR,
#                    help="Depth variable name in --nc_path (2D, meters; >0 = water). Default: depth")
#
#     g.add_argument("--rivers_var", type=str, default=DEF_RIVERS_VAR,
#                    help="Rivers mask variable name in --nc_path (2D, 1=river-mouth cells).")
#
#     g.add_argument("--offshore_var", type=str, default=DEF_OFFSHORE_VAR,
#                    help="Offshore mask variable name in --nc_path (2D, 1=offshore cells).")
#
#     g.add_argument("--rivers_months", type=str, default=DEF_RIVERS_MONTHS_STR,
#                    help="Months to apply rivers mask: 'all' | 'none' | '9-12' | '1,2,3,6-8'. Default: none")
#
#     g.add_argument("--offshore_months", type=str, default=DEF_OFFSHORE_MONTHS_STR,
#                    help="Months to apply offshore mask: 'all' | 'none' | '1-6', etc. Default: none")
#
#     g.add_argument("--rivers_weight", type=float, default=DEF_RIVERS_WEIGHT,
#                    help="Multiplicative weight where rivers mask is 1 (when active). Default: 1.5")
#
#     g.add_argument("--offshore_weight", type=float, default=DEF_OFFSHORE_WEIGHT,
#                    help="Multiplicative weight where offshore mask is 1 (when active). Default: 1.2")
#
#     g.add_argument("--combine_mode", type=str, choices=["multiply", "overwrite"],
#                    default=DEF_COMBINE_MODE,
#                    help=("How to combine masks with the month grid. "
#                          "'multiply' gently boosts habitat (recommended for herring). "
#                          "'overwrite' forces a salmon-style 1-p replacement."))
#
#     g.add_argument("--gauss_sigma", type=float, default=DEF_GAUSS_SIGMA,
#                    help="Gaussian blur sigma (in cells). 0 disables smoothing. Default: 0")
#
#     g.add_argument("--min_floor", type=float, default=DEF_MIN_FLOOR,
#                    help="Minimum non-zero floor applied to water cells after normalisation. Default: 0.001")
#
#     d = p.add_argument_group("Depth preference (optional)")
#
#     # Depth preference toggle (default from DEF_DEPTH_PREF_ON)
#     d.add_argument("--depth_pref_on", dest="depth_pref_on", action="store_true",
#                    help="Enable depth preference weighting.")
#     d.add_argument("--no-depth_pref_on", dest="depth_pref_on", action="store_false",
#                    help="Disable depth preference weighting.")
#     p.set_defaults(depth_pref_on=DEF_DEPTH_PREF_ON)
#
#     d.add_argument("--depth_bands", type=str, default=DEF_DEPTH_BANDS_STR,
#                    help="Preferred depth band(s) like '0-15' or '0-15,25-45' (meters).")
#
#     d.add_argument("--depth_peak_w", type=float, default=DEF_DEPTH_PEAK_W,
#                    help="Peak multiplicative weight at the center of each band. Default: 2.0")
#
#     d.add_argument("--depth_floor_w", type=float, default=DEF_DEPTH_FLOOR_W,
#                    help="Minimum multiplicative weight (outside bands). Default: 1.0")
#
#
#     args = p.parse_args(argv)  # parse ONCE
#     # post-process compact strings
#     args.rivers_months = parse_months(args.rivers_months)
#     args.offshore_months = parse_months(args.offshore_months)
#     args.depth_bands = parse_depth_bands(args.depth_bands)
#
#     csv_to_monthly_ascs(
#         spawn_csv=args.spawn_csv,
#         out_dir=args.out_dir,
#         example_asc=args.example_asc,
#         value_col=args.value_col,
#         date_col=args.date_col,
#         year=args.year,
#         fmt=args.fmt,
#         plume_mask_path=args.plume_mask,
#         land_mask_path=args.land_mask,
#         index_base=args.index_base,
#         output_basename=args.basename,
#         verbose=(not args.quiet),
#         log_norm=args.log_norm,
#         log_eps=args.log_eps,
#         binary=args.binary,
#         binary_threshold=args.binary_threshold,
#         min_fill=args.min_floor,          # CLI --min_floor -> function min_fill
#         presence_prob=args.presence_prob,
#         years_denominator=args.years_denominator,
#
#         # NC & weighting
#         nc_path=args.nc_path,
#         depth_var=args.depth_var,
#         rivers_var=args.rivers_var,
#         offshore_var=args.offshore_var,
#         rivers_months=args.rivers_months,
#         offshore_months=args.offshore_months,
#         rivers_weight=args.rivers_weight,
#         offshore_weight=args.offshore_weight,
#         combine_mode=args.combine_mode,
#         gauss_sigma=args.gauss_sigma,
#
#         # depth preference
#         depth_pref_on=args.depth_pref_on,
#         depth_bands=args.depth_bands,
#         depth_peak_w=args.depth_peak_w,
#         depth_floor_w=args.depth_floor_w,
#     )
#
#
# if __name__ == "__main__":
#     sys.exit(main())