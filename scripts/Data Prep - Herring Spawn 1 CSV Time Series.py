"""
Match Pacific herring spawning locations to Ecospace (row, col)

Created by G Oldford Oct 2025

Purpose
-------
Given:
  1) a CSV of herring spawn locations (with columns 'Latitude' and 'Longitude'), and
  2) an Ecospace grid CSV containing at least row/col, lat/lon, and depth,

this script finds, for each spawn record, the *nearest* Ecospace grid cell
with depth > 0 and appends:
    - ecospace_row
    - ecospace_col
    - ewe_lat
    - ewe_lon
    - ewe_depth
    - distance_km  (great-circle distance to matched grid point)

Data In:
---------
herring CSV - https://open.canada.ca/data/en/dataset/d892511c-d851-4f85-a0ec-708bc05d2810
- modified by GO to add column and compute 'Spawn Intensity' (L x W x EggDens)
ecospace table of rows cols and _real_ lats lons (rotation makes lats lons unreliable)

Usage (example)
---------------
python match_herring_spawn_to_ewe_rowcol.py \
  --spawn_csv "/path/to/Pacific_herring_spawn_index_data_2024_EN_GO.csv" \
  --grid_csv  "/path/to/basemaps/ecospace_regions_coords_rowscols_cfg1.csv" \
  --output_csv "/path/to/output/Pacific_herring_spawn_tagged_cfg1.csv" \
  --max_km1 1.5 --max_km2 3.0

Notes
-----
- We match to the *closest* grid point with depth > 0.
- Optional --max_km lets you flag improbable matches (they'll still match to
  the nearest ocean cell but will be flagged via a boolean column 'over_max_km').

"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd

# ============================
# CONFIG: Default parameters
# ============================
DEFAULT_SPAWN_P    = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//")
DEFAULT_SPAWN_CSV  = Path(DEFAULT_SPAWN_P, "Pacific_herring_spawn_index_data_2024_EN_GO.csv")
DEFAULT_GRID_CSV   = Path("..//data//basemap//Ecospace_grid_20210208_rowscols.csv")
DEFAULT_OUTPUT_CSV = Path(DEFAULT_SPAWN_P, "Pacific_SoG_herring_spawn_index_data_2024_EN_taggedEwERowsCols.csv")
DEFAULT_MAX_KM     = 5.0  # purely for flagging 'over_max_km', not a hard cutoff
DEFAULT_RETAIN_OVER_MAX = False

def _log(msg: str) -> None:
    print(f"[INFO] {msg}", flush=True)


def haversine_km(lat1, lon1, lat2, lon2):
    """
    Great-circle distance (km) using the Haversine formula.
    Works with scalar-scalar, scalar-array, array-scalar, or array-array.
    """
    R = 6371.0  # Earth radius in km

    # Ensure numpy arrays (broadcastable), then convert to radians
    lat1 = np.radians(np.asarray(lat1, dtype=float))
    lon1 = np.radians(np.asarray(lon1, dtype=float))
    lat2 = np.radians(np.asarray(lat2, dtype=float))
    lon2 = np.radians(np.asarray(lon2, dtype=float))

    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    return 2.0 * R * np.arcsin(np.sqrt(a))

def _resolve_col(df: pd.DataFrame, candidates) -> str | None:
    """Return the first matching column name from 'candidates' present in df, else None."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _normalize_grid_columns(df_grid: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to standard names: row, col, lat, lon, depth (where possible)."""
    rename_map = {}

    # latitude / longitude
    lat_col = _resolve_col(df_grid, ["lat", "latitude", "Latitude", "LAT", "Lat"])
    lon_col = _resolve_col(df_grid, ["lon", "longitude", "Longitude", "LON", "Lon"])
    if lat_col: rename_map[lat_col] = "lat"
    if lon_col: rename_map[lon_col] = "lon"

    # row / col
    row_col = _resolve_col(df_grid, ["row", "Row", "EWE_row", "ROW"])
    col_col = _resolve_col(df_grid, ["col", "Col", "EWE_col", "COL"])
    if row_col: rename_map[row_col] = "row"
    if col_col: rename_map[col_col] = "col"

    # depth
    depth_col = _resolve_col(df_grid, ["depth", "Depth", "DEPTH", "depth_m", "Depth_m", "bathymetry", "Bathy", "BATHY"])
    if depth_col: rename_map[depth_col] = "depth"

    df_grid = df_grid.rename(columns=rename_map)

    required = {"row", "col", "lat", "lon"}
    missing = required.difference(df_grid.columns)
    if missing:
        raise ValueError(f"Grid CSV must contain columns {sorted(required)} (case-insensitive variants allowed). Missing: {sorted(missing)}")

    # if depth is missing, assume ocean everywhere (depth>0). Otherwise coerce to numeric.
    if "depth" not in df_grid.columns:
        df_grid["depth"] = np.inf  # treat as ocean
    else:
        df_grid["depth"] = pd.to_numeric(df_grid["depth"], errors="coerce").fillna(0.0)

    # ensure numeric types
    df_grid["lat"] = pd.to_numeric(df_grid["lat"], errors="coerce")
    df_grid["lon"] = pd.to_numeric(df_grid["lon"], errors="coerce")
    df_grid["row"] = pd.to_numeric(df_grid["row"], errors="coerce").astype("Int64")
    df_grid["col"] = pd.to_numeric(df_grid["col"], errors="coerce").astype("Int64")

    # drop invalid coordinate rows
    df_grid = df_grid.dropna(subset=["lat", "lon"]).reset_index(drop=True)
    return df_grid


def _prepare_spawn(df_spawn: pd.DataFrame) -> pd.DataFrame:
    """Ensure spawn CSV has usable Latitude/Longitude."""
    # tolerate various casings
    if "Latitude" not in df_spawn.columns or "Longitude" not in df_spawn.columns:
        # attempt to find alternates
        lat_alt = [c for c in df_spawn.columns if c.lower() == "latitude"]
        lon_alt = [c for c in df_spawn.columns if c.lower() == "longitude"]
        if lat_alt: df_spawn = df_spawn.rename(columns={lat_alt[0]: "Latitude"})
        if lon_alt: df_spawn = df_spawn.rename(columns={lon_alt[0]: "Longitude"})

    if "Latitude" not in df_spawn.columns or "Longitude" not in df_spawn.columns:
        raise ValueError("Spawn CSV must contain 'Latitude' and 'Longitude' columns.")

    df_spawn["Latitude"] = pd.to_numeric(df_spawn["Latitude"], errors="coerce")
    df_spawn["Longitude"] = pd.to_numeric(df_spawn["Longitude"], errors="coerce")
    n0 = len(df_spawn)
    df_spawn = df_spawn.dropna(subset=["Latitude", "Longitude"]).reset_index(drop=True)
    if len(df_spawn) < n0:
        _log(f"Dropped {n0 - len(df_spawn)} rows with missing/invalid coordinates.")
    return df_spawn


def tag_herring_to_rowscols(spawn_fp: Path, grid_fp: Path, output_fp: Path, max_km_flag: float, retain_over_max: bool) -> None:
    """Core routine: read inputs, nearest ocean-cell match (depth>0), write output."""
    _log(f"Reading spawn CSV: {spawn_fp}")
    df_spawn = pd.read_csv(spawn_fp)
    df_spawn = _prepare_spawn(df_spawn)

    _log(f"Reading grid CSV: {grid_fp}")
    df_grid = pd.read_csv(grid_fp)
    df_grid = _normalize_grid_columns(df_grid)

    # Only consider *ocean* cells (depth > 0); drop depth<=0 and missing
    df_ocean = df_grid[df_grid["depth"] > 0].copy()
    if df_ocean.empty:
        raise ValueError("After filtering depth>0, no grid cells remain. Check 'depth' column mapping.")

    _log(f"Matching {len(df_spawn)} rows to nearest (row,col) among {len(df_ocean)} ocean cells...")

    # Vectorized nearest neighbor via brute-force (adequate for typical grid sizes).
    # If performance becomes an issue, we can swap to a KD-tree in projected meters.
    ocean_lat = df_ocean["lat"].to_numpy(dtype=float)
    ocean_lon = df_ocean["lon"].to_numpy(dtype=float)

    matches = []
    for _, r in df_spawn.iterrows():
        lat = float(r["Latitude"])
        lon = float(r["Longitude"])
        dists = haversine_km(lat, lon, ocean_lat, ocean_lon)
        idx = int(np.argmin(dists))
        best = df_ocean.iloc[idx]

        matches.append({
            "ecospace_row": int(best["row"]),
            "ecospace_col": int(best["col"]),
            "ewe_lat": float(best["lat"]),
            "ewe_lon": float(best["lon"]),
            "ewe_depth": float(best["depth"]),
            "distance_km": float(dists[idx]),
            "over_max_km": bool(dists[idx] > max_km_flag) if max_km_flag is not None else False,
        })

    df_out = pd.concat([df_spawn, pd.DataFrame(matches)], axis=1)

    # Drop rows over max_km unless retain_over_max
    if max_km_flag is not None and not retain_over_max:
        n_before = len(df_out)
        df_out = df_out[df_out["distance_km"] <= float(max_km_flag)].reset_index(drop=True)
        _log(f"Filtered out {n_before - len(df_out)} rows with distance_km > {max_km_flag}.")

    # Summary
    n_total = len(df_spawn)
    n_out = len(df_out)
    n_flagged = int((df_out["over_max_km"] == True).sum()) if "over_max_km" in df_out.columns else 0
    _log(f"Summary: input={n_total}, kept={n_out}, flagged_over_max_km={n_flagged}")

    output_fp.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(output_fp, index=False)
    _log(f"Saved tagged CSV: {output_fp}")



def main(argv=None):
    p = argparse.ArgumentParser(description="Match herring spawning locations to Ecospace (row,col) with depth>0.")
    p.add_argument("--spawn_csv", type=Path, default=DEFAULT_SPAWN_CSV, help="Path to herring spawn CSV (Latitude,Longitude)")
    p.add_argument("--grid_csv", type=Path, default=DEFAULT_GRID_CSV, help="Path to Ecospace grid CSV (row,col,lat,lon,depth)")
    p.add_argument("--output_csv", type=Path, default=DEFAULT_OUTPUT_CSV, help="Where to write the tagged CSV")
    p.add_argument("--max_km", type=float, default=DEFAULT_MAX_KM, help="Flag rows with match distance > this (adds 'over_max_km')")
    p.add_argument("--retain_over_max", default=DEFAULT_RETAIN_OVER_MAX, action="store_true",
                   help="Keep rows farther than --max_km (they will be flagged but not dropped).")
    args = p.parse_args(argv)

    tag_herring_to_rowscols(
        spawn_fp=args.spawn_csv,
        grid_fp=args.grid_csv,
        output_fp=args.output_csv,
        max_km_flag=args.max_km,
        retain_over_max=args.retain_over_max,
    )


if __name__ == "__main__":
    sys.exit(main())