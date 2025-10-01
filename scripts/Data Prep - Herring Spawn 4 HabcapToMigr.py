#!/usr/bin/env python
"""
Rescale 12 monthly HERRING_HABCAP climatology layers by per-cell annual max,
write back to NC as new static vars, and export each month to ASC using your
existing helper functions & header template.

Depends on:
  - helpers_ewe_asc_netcdf.py  (ASC/NC helpers)
  - helpers.py                  (general utils)

Example:
  python rescale_herring_habcap_clim_by_cellmax.py \
    --nc_in "..//..//ecospace_seal_diets//Data//basemaps//ecospace_bundled_asc.nc" \
    --nc_out SAME \
    --prefix "herring_habcap_climatology_m" \
    --out_suffix "_scaled" \
    --asc_template "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//some_template.asc" \
    --asc_out_dir "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//MIGR_ASC//" \
    --overwrite
"""

from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import numpy as np
from netCDF4 import Dataset

# Reuse your helpers (ASC header I/O, simple ASC writer, etc.)
from helpers_ewe_asc_netcdf import (
    read_asc_header_and_dims   # header_text, nrows, ncols
)
from helpers import (           # general helpers / your style
    saveASCFile
)

# ------------------------------
# Defaults (your style at top)
# ------------------------------
DEF_NC_IN      = Path("..//..//ecospace_seal_diets//Data//basemaps//ecospace_bundled_asc.nc")
DEF_NC_OUT     = "SAME"   # "SAME" => in-place update
DEF_PREFIX     = "herring_habcap_climatology_m"
DEF_MONTHS     = "1-12"
DEF_OUT_PREFIX = "herring_migration_"
DEF_OUT_SUFFIX = "_scaled"
DEF_OVERWRITE  = True

DEF_FLOOR      = 0.01
DEF_MAX        = 4.0

DEF_BINARY     = True
DEF_BIN_THRES  = 0.15

# ASC export
DEF_ASC_OUT_DIR      = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//MIGR_ASC//")
DEF_TEMPLATE_ASC_P   = Path("..//data//forcing//ECOSPACE_in_climayr_1980-2018_PAR3_Sal4m_20240523//PAR-VarZ-VarK//")
DEF_ASC_TEMPLATE     = Path(DEF_TEMPLATE_ASC_P, "PAR-VarZ-VarK_1980-2018.asc") # REQUIRED to export; use an existing ASC with desired header
DEF_ASC_FMT          = "%0.3f"       # matches your usual formatting

# ------------------------------
# Minimal new logic
# ------------------------------
def parse_months(spec: str) -> list[int]:
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            a = int(a); b = int(b)
            out.extend(range(min(a,b), max(a,b)+1))
        else:
            out.append(int(part))
    out = sorted(set([m for m in out if 1 <= m <= 12]))
    return out


def rescale_by_cellmax(stack: np.ndarray) -> np.ndarray:
    """
    stack: (12, y, x) float array. Returns same shape.
    Per-cell: divide by annual max; if max<=0 -> set to 0 to avoid div-by-zero.
    """
    annual_max = np.nanmax(stack, axis=0)  # ignore NaNs
    denom = np.where(annual_max > 0, annual_max, 1.0)
    out = stack / denom
    out[:, annual_max <= 0] = 0.0
    return out


def rescale_by_monthly_mapmax(stack: np.ndarray, target_max: float) -> np.ndarray:
    """
    Scale each month by its own map-wide maximum so that
    max(month) == target_max (e.g., 1.0, 3.0, 4.0).
    stack: (12, y, x) float array (NaNs allowed). Returns same shape.
    """
    out = np.empty_like(stack, dtype=float)
    for i in range(stack.shape[0]):
        month = stack[i]
        mmax = np.nanmax(month)
        if not np.isfinite(mmax) or mmax <= 0:
            # no information or all non-positive: keep zeros
            out[i] = 0.0
        else:
            scale = target_max / mmax
            out[i] = month * scale
    return out


def apply_floor_ceiling(arr: np.ndarray, floor: float | None, ceiling: float | None, mode: str) -> np.ndarray:
    """
    Clamp values to [floor, ceiling] with control over how the floor is applied.
    - mode='valid': apply floor to every unmasked/finite cell
    - mode='nonzero': apply floor only where arr > 0
    - mode='none': skip floor
    Masked/NaN cells are untouched.
    """
    out = arr.copy()
    valid = np.isfinite(out)

    if ceiling is not None:
        out[valid] = np.minimum(out[valid], ceiling)

    if floor is not None and mode != "none":
        if mode == "valid":
            out[valid] = np.maximum(out[valid], floor)
        elif mode == "nonzero":
            nz = valid & (out > 0)
            out[nz] = np.maximum(out[nz], floor)
        else:
            raise ValueError(f"Unknown floor_mode: {mode}")

    return out

# ------------------------------
# I/O helpers that only wrap netCDF read/write (no duplication of your helpers)
# ------------------------------
def read_month_2d(ds: Dataset, vname: str) -> np.ndarray:
    if vname not in ds.variables:
        raise KeyError(f"Missing var: {vname}")
    v = ds.variables[vname]
    arr = np.array(v[:], dtype=float).squeeze()
    if arr.ndim != 2:
        raise ValueError(f"{vname} not 2D. Got shape {arr.shape}")
    return arr


def write_month_2d(ds, like_var: str, out_name: str, arr: np.ndarray, overwrite: bool):
    """
    Safe writer:
    - If var exists and overwrite=False -> raise (so you notice).
    - If var exists and overwrite=True  -> update data/attrs in place (no recreate).
    - If var missing                    -> create with compression.
    """
    import numpy as np

    # sanity
    if like_var not in ds.variables:
        raise KeyError(f"Reference var not found: {like_var}")
    vref = ds.variables[like_var]
    dims = vref.dimensions
    dtype = vref.dtype
    fv = getattr(vref, "_FillValue", np.nan if np.issubdtype(dtype, np.floating) else -9999)

    if out_name in ds.variables:
        if not overwrite:
            raise ValueError(f"Var exists: {out_name}. Use --overwrite to replace.")
        # ---- Update-in-place path (no recreate) ----
        vout = ds.variables[out_name]
        # If shape mismatch, you can either error or try to proceed; error is safer:
        if vout.shape != (len(ds.dimensions[dims[0]]), len(ds.dimensions[dims[1]])):
            raise ValueError(f"Existing var {out_name} has shape {vout.shape}, "
                             f"expected {(len(ds.dimensions[dims[0]]), len(ds.dimensions[dims[1]]))}. "
                             f"Delete it manually or change --out_suffix.")
        # Write data and refresh a couple attrs for clarity
        vout[:, :] = np.nan_to_num(arr, nan=fv)
        vout.long_name = f"{getattr(vref, 'long_name', out_name)} (rescaled by per-cell annual max)"
        if hasattr(vref, "units"):
            vout.units = getattr(vref, "units")
        if hasattr(vref, "grid_mapping"):
            vout.grid_mapping = getattr(vref, "grid_mapping")
        # done
        return

    # ---- Create fresh variable path ----
    vout = ds.createVariable(out_name, dtype, dims, zlib=True, complevel=4, fill_value=fv)
    for attr in ("units", "grid_mapping"):
        if hasattr(vref, attr):
            setattr(vout, attr, getattr(vref, attr))
    vout.long_name = f"{getattr(vref, 'long_name', out_name)} (rescaled by per-cell annual max)"
    vout[:, :] = np.nan_to_num(arr, nan=fv)


def binarize_array(arr: np.ndarray, threshold: float, nan_to_zero: bool = False) -> np.ndarray:
    """
    Return a binary array where:
      - values >= threshold -> 1
      - values  < threshold -> 0
    NaN handling:
      - if nan_to_zero is False (default), leave NaNs as NaN (preserves land masks)
      - if True, set NaNs to 0
    """
    out = np.zeros_like(arr, dtype=float)

    # Identify finite & NaN
    finite = np.isfinite(arr)
    out[finite] = (arr[finite] >= threshold).astype(float)

    if nan_to_zero:
        out[~finite] = 0.0
    else:
        # keep NaNs as NaN (don’t overwrite)
        out[~finite] = np.nan

    # out *= 10 # experiment w/ increasing migration values

    return out

# ------------------------------
# Main
# ------------------------------
def main():
    ap = argparse.ArgumentParser(description="Rescale HERRING_HABCAP climatology by per-cell annual max; write NC + export ASC.")
    ap.add_argument("--nc_in", default=str(DEF_NC_IN))
    ap.add_argument("--nc_out", default=DEF_NC_OUT, help='Output NC or "SAME" for in-place')
    ap.add_argument("--prefix", default=DEF_PREFIX)
    ap.add_argument("--months", default=DEF_MONTHS)

    ap.add_argument("--floor", type=float, default=DEF_FLOOR,
                    help="Minimum value to enforce (e.g., 0.001).")
    ap.add_argument("--target_max", type=float,
                    default=DEF_MAX,
                    help="Scale each month so its map-wide maximum equals this value (e.g., 1, 3, 4)."
    )
    ap.add_argument("--floor_mode", choices=["none", "nonzero", "valid"], default="none",
                    help="'valid' applies floor to all unmasked cells; 'nonzero' only to cells > 0.")

    ap.add_argument("--out_prefix", default=DEF_OUT_PREFIX)
    ap.add_argument("--out_suffix", default=DEF_OUT_SUFFIX)
    ap.add_argument("--overwrite", action="store_true", default=DEF_OVERWRITE)

    # ASC export
    ap.add_argument("--asc_out_dir", default=str(DEF_ASC_OUT_DIR))
    ap.add_argument("--asc_template", default=DEF_ASC_TEMPLATE, help="Path to an existing ASC with the desired header")
    ap.add_argument("--asc_fmt", default=DEF_ASC_FMT)

    ap.add_argument("--binary", action="store_true", default=DEF_BINARY,
                    help="If set, convert outputs to binary 0/1 using --binary_threshold.")
    ap.add_argument("--binary_threshold", type=float, default=DEF_BIN_THRES,
                    help="Cells >= threshold -> 1, below -> 0 (used only with --binary).")
    ap.add_argument("--binary_nan_to_zero", action="store_true", default=True,
                    help="When --binary, convert NaNs to 0 instead of leaving them as NaN.")

    args = ap.parse_args()
    months = parse_months(args.months)

    nc_in = Path(args.nc_in)
    if not nc_in.exists():
        ap.error(f"NC not found: {nc_in}")

    # Handle nc_out (copy-on-write if different target)
    if args.nc_out == "SAME":
        nc_out = nc_in
    else:
        nc_out = Path(args.nc_out)
        if nc_out.resolve() != nc_in.resolve():
            shutil.copy2(nc_in, nc_out)

    # Open writable
    with Dataset(nc_out, "r+") as ds:
        # Load monthly layers
        monthly = []
        like_var = None
        for m in months:
            vname = f"{args.prefix}{m:02d}"
            arr = read_month_2d(ds, vname)
            if like_var is None:
                like_var = vname
            monthly.append(arr)

        # Stack -> rescale
        stack = np.stack(monthly, axis=0)  # (12, y, x)
        rescaled = rescale_by_monthly_mapmax(stack, args.target_max)

        # Write back to NC (one 2D static var per month)
        for i, m in enumerate(months):

            src = f"{args.prefix}{m:02d}"
            out_name = f"{args.out_prefix}{m:02d}{args.out_suffix}"

            # Optional floor logic (only if you kept it)
            # arr = apply_floor_ceiling(arr, args.floor, None, args.floor_mode)

            arr = rescaled[i]

            # Optional floor/ceiling if you still want it (left commented as before)
            # arr = apply_floor_ceiling(arr, args.floor, None, args.floor_mode)

            # Optional binary conversion
            if args.binary:
                arr = binarize_array(arr, threshold=args.binary_threshold, nan_to_zero=args.binary_nan_to_zero)

            write_month_2d(ds, like_var=src, out_name=out_name, arr=arr, overwrite=args.overwrite)

            print(f"[OK] wrote {out_name} to {nc_out.name}")

    # ASC export (optional, but you asked for it)
    asc_out_dir = Path(args.asc_out_dir) if args.asc_out_dir else None
    if asc_out_dir:
        if not args.asc_template:
            raise ValueError("--asc_template is required to export ASC (for header consistency).")
        header_text, nrows, ncols = read_asc_header_and_dims(Path(args.asc_template))  # your helper
        asc_out_dir.mkdir(parents=True, exist_ok=True)

        # Reopen (read) to avoid keeping a handle open
        with Dataset(nc_out, "r") as ds:
            for m in months:
                # was: vname = f"{args.prefix}{m:02d}{args.out_suffix}"
                # Should match how you *wrote* them (out_prefix + month + out_suffix)
                vname = f"{args.out_prefix}{m:02d}{args.out_suffix}"
                if vname not in ds.variables:
                    raise KeyError(f"Expected output var missing: {vname}")
                arr = np.array(ds.variables[vname][:], dtype=float)
                # If shapes differ, still write—your helper header forces target nrows/ncols.
                # But better to alert loudly so you can fix ASAP.
                print(f"[WARN] {vname} shape {arr.shape} != {(nrows, ncols)} (template). Writing anyway.")
                out_path = asc_out_dir / f"HERRING_MIGR_CLIM_{m:02d}_scaled.asc"

                # save_asc_plain(out_path, arr, header_text, fmt=args.asc_fmt)  # your helper writer
                saveASCFile(
                    str(out_path),
                    np.flipud(arr),
                    bottomleft_row_ewe=0,
                    upperleft_row_ewe=nrows,
                    upperleft_col_ewe=0,
                    sigdigfmt=args.asc_fmt,
                    ASCheader=header_text,
                    dfPlumeMask=None,
                    dfLandMask=None
                )

                print(f"[OK] exported {out_path}")

    print("[INFO] Done.")

if __name__ == "__main__":
    main()
