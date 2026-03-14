
# helpers for moving between netcdf and ASC for Ecospace


from __future__ import annotations

from pathlib import Path
import re
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd
from datetime import datetime
from netCDF4 import Dataset, date2num, num2date
from scipy.ndimage import gaussian_filter


# Lazy import to avoid hard dependency when just reading/writing ASC.
def _nc4():
    from netCDF4 import Dataset, date2num, num2date
    return Dataset, date2num, num2date

# ----------------------------------------------------
# ESRI ASCII I/O
# ----------------------------------------------------

def _parse_asc_header(lines: List[str]) -> Dict[str, float]:
    """Parse the first 6 header lines of ESRI ASCII grid.
    Supports xllcorner/xllcenter and yllcorner/yllcenter.
    """
    header = {}
    for i in range(6):
        key_val = lines[i].strip().split()
        if len(key_val) >= 2:
            key = key_val[0].lower()
            try:
                val = float(key_val[1])
            except ValueError:
                # Occasionally headers include ints or unexpected tokens; try safe cast
                try:
                    val = float(int(key_val[1]))
                except Exception:
                    raise
            header[key] = val
    # Normalize keys
    if "ncols" not in header or "nrows" not in header:
        raise ValueError("ASC header must include ncols and nrows.")
    if "xllcorner" not in header and "xllcenter" not in header:
        raise ValueError("ASC header must include xllcorner or xllcenter.")
    if "yllcorner" not in header and "yllcenter" not in header:
        raise ValueError("ASC header must include yllcorner or yllcenter.")
    if "cellsize" not in header:
        raise ValueError("ASC header must include cellsize.")
    if "nodata_value" not in header:
        # default (some files omit)
        header["nodata_value"] = -9999.0
    return header


def read_esri_ascii(path: str | Path) -> Dict[str, Any]:
    """Read an ESRI ASCII grid file into numpy array with header metadata.
    Returns dict with keys: 'header', 'data' (2-D array, float64 with NaN for NODATA).
    """
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    header = _parse_asc_header(lines)
    ncols = int(header["ncols"]); nrows = int(header["nrows"])
    nodata = header.get("nodata_value", -9999.0)
    # Try fast path: one line per row after header
    try:
        data = np.loadtxt(path, skiprows=6)
        if data.shape != (nrows, ncols):
            raise ValueError("Unexpected shape from loadtxt; falling back to manual parse")
    except Exception:
        data_vals = " ".join(l.strip() for l in lines[6:]).split()
        if len(data_vals) != ncols * nrows:
            raise ValueError(f"ASC appears malformed: expected {nrows*ncols} numbers, found {len(data_vals)}")
        data = np.array(list(map(float, data_vals)), dtype=float).reshape(nrows, ncols)
    data = data.astype(float)
    data[data == nodata] = np.nan
    return {"header": header, "data": data}


def write_esri_ascii(path: str | Path, data: np.ndarray, header_template: Dict[str, float]) -> None:
    """Write a numpy array to ESRI ASCII using a provided header (must include ncols, nrows, cellsize, origin, nodata_value)."""
    path = Path(path)
    nrows, ncols = data.shape
    # Ensure header matches data shape
    hdr = dict(header_template)  # copy
    hdr["nrows"] = float(nrows)
    hdr["ncols"] = float(ncols)
    nodata = hdr.get("nodata_value", -9999.0)
    # Build header text
    def kv(k, v): return f"{k} {int(v) if k in ('ncols','nrows') else v}"
    keys_out = []
    # Keep xllcorner/xllcenter choice as provided
    if "xllcorner" in hdr:
        keys_out.extend([
            kv("ncols", hdr["ncols"]),
            kv("nrows", hdr["nrows"]),
            kv("xllcorner", hdr["xllcorner"]),
            kv("yllcorner", hdr["yllcorner"]),
            kv("cellsize", hdr["cellsize"]),
            kv("NODATA_value", nodata),
        ])
    else:
        keys_out.extend([
            kv("ncols", hdr["ncols"]),
            kv("nrows", hdr["nrows"]),
            kv("xllcenter", hdr["xllcenter"]),
            kv("yllcenter", hdr["yllcenter"]),
            kv("cellsize", hdr["cellsize"]),
            kv("NODATA_value", nodata),
        ])
    # Replace NaN with nodata for output
    out = np.array(data, copy=True)
    out[np.isnan(out)] = nodata
    with path.open("w", encoding="utf-8") as f:
        f.write("\n".join(keys_out) + "\n")
        np.savetxt(f, out, fmt="%.6f")

# ----------------------------------------------------
# Grid utilities
# ----------------------------------------------------

def build_latlon_from_rowscols_csv(csv_path: str | Path,
                                   row_col_names: Tuple[str, str] = ("EWE_row","EWE_col"),
                                   lat_lon_names: Tuple[str, str] = ("lat","lon"),
                                   y_first: bool = True) -> Tuple[np.ndarray, np.ndarray, Tuple[int,int], Dict[str, Any]]:
    """Build 2-D lat/lon arrays from a CSV that has row, col, lat, lon.
    Returns (lat2d, lon2d, (nrows, ncols), meta).
    If y_first is True, arrays are shaped (nrows, ncols) with row as y and col as x.
    """
    df = pd.read_csv(csv_path)
    rname, cname = row_col_names
    latname, lonname = lat_lon_names
    # Attempt to auto-detect column names if not present
    cols_lower = {c.lower(): c for c in df.columns}
    for want, default in [(rname, "row"), (cname, "col"), (latname, "lat"), (lonname, "lon")]:
        if want not in df.columns:
            cand = cols_lower.get(want.lower()) or cols_lower.get((want+"itude").lower())
            if cand:
                if want == rname: rname = cand
                elif want == cname: cname = cand
                elif want == latname: latname = cand
                elif want == lonname: lonname = cand
    # Detect 1-based row/col and normalize to 0-based for array placement
    rmin, cmin = int(df[rname].min()), int(df[cname].min())
    one_based = (rmin == 1 and cmin == 1)
    rr = df[rname].astype(int).to_numpy() - (1 if one_based else 0)
    cc = df[cname].astype(int).to_numpy() - (1 if one_based else 0)
    nrows = max(rr)+1; ncols = max(cc)+1
    lat2d = np.full((nrows, ncols), np.nan, dtype=float)
    lon2d = np.full((nrows, ncols), np.nan, dtype=float)
    lat2d[rr, cc] = df[latname].to_numpy(dtype=float)
    lon2d[rr, cc] = df[lonname].to_numpy(dtype=float)
    meta = {"one_based": one_based, "row_name": rname, "col_name": cname, "lat_name": latname, "lon_name": lonname}
    return lat2d, lon2d, (nrows, ncols), meta

# ----------------------------------------------------
# Time parsing
# ----------------------------------------------------

_DATE_PATTERNS = [
    r"(?P<year>\d{4})[-_](?P<month>\d{1,2})",   # 2005-4 or 2005_04
    r"(?P<year>\d{4})(?P<month>\d{2})",         # 200504
    r"(?P<year>\d{4})",                         # year only (month defaults to 1)
]

def guess_time_from_filename(name: str) -> Optional[datetime]:
    """Extract a datetime from a filename using common patterns. Returns None if not found."""
    base = Path(name).stem
    for pat in _DATE_PATTERNS:
        m = re.search(pat, base)
        if m:
            year = int(m.group("year"))
            month = int(m.groupdict().get("month") or 1)
            day = 1
            try:
                return datetime(year, month, day)
            except Exception:
                return None
    return None

# ----------------------------------------------------
# NetCDF writing/reading
# ----------------------------------------------------

def write_netcdf(path_out: str | Path,
                 lat2d: np.ndarray, lon2d: np.ndarray,
                 static_vars: Dict[str, Tuple[np.ndarray, Dict[str, Any]]] | None,
                 time_series_vars: List[Tuple[str, List[datetime], np.ndarray, Dict[str, Any]]] | None,
                 global_attrs: Dict[str, Any] | None = None,
                 calendar: str = "proleptic_gregorian") -> None:
    """Write a CF-1.8 NetCDF with coordinates and variables.

    Parameters
    ----------
    path_out : output NetCDF path
    lat2d, lon2d : 2-D arrays of shape (y, x)
    static_vars : dict name -> (2-D array [y,x], var_attrs)
    time_series_vars : list of tuples (var_name, times, data_stack, var_attrs),
        where data_stack has shape (T, y, x)
    global_attrs : dict of attributes for the file
    calendar : calendar for the 'time' variable
    """
    Dataset, date2num, _ = _nc4()
    path_out = Path(path_out)
    ny, nx = lat2d.shape
    ds = Dataset(path_out, "w", format="NETCDF4")
    try:
        # Dimensions
        ds.createDimension("y", ny)
        ds.createDimension("x", nx)
        has_time = bool(time_series_vars)
        if has_time:
            ds.createDimension("time", None)
            time_var = ds.createVariable("time", "f8", ("time",))
            time_var.units = "days since 1900-01-01 00:00:00"
            time_var.calendar = calendar
            time_var.long_name = "time"

        # Coordinates
        latv = ds.createVariable("lat", "f8", ("y" ,"x"), zlib=True, complevel=4, fill_value=np.nan)
        lonv = ds.createVariable("lon", "f8", ("y" ,"x"), zlib=True, complevel=4, fill_value=np.nan)
        latv[:, :] = lat2d
        lonv[:, :] = lon2d
        latv.standard_name = "latitude"; latv.units = "degrees_north"
        lonv.standard_name = "longitude"; lonv.units = "degrees_east"

        # Static variables
        if static_vars:
            for name, (arr, attrs) in static_vars.items():
                if arr.shape != (ny, nx):
                    raise ValueError(f"Static var '{name}' shape {arr.shape} != {(ny ,nx)}")
                v = ds.createVariable(name, "f4", ("y" ,"x"), zlib=True, complevel=4, fill_value=np.nan)
                v[:, :] = arr
                v.coordinates = "lat lon"
                if attrs:
                    for k, val in attrs.items():
                        setattr(v, k, val)

        # Time series variables
        if has_time and time_series_vars:
            all_times = []
            for (vname, times, data_stack, attrs) in time_series_vars:
                all_times.extend(times)
            uniq_times = sorted(set(all_times))
            # Write time coordinate
            ds_times = date2num(uniq_times, units=time_var.units, calendar=time_var.calendar)
            time_var[:] = ds_times
            # Index map for quick placement
            time_index = {t: i for i, t in enumerate(uniq_times)}
            # Write each series aligned to global time
            for (vname, times, data_stack, attrs) in time_series_vars:
                if data_stack.shape[1:] != (ny, nx):
                    raise ValueError(f"Var '{vname}' stack spatial shape {data_stack.shape[1:]} != {(ny ,nx)}")
                T = len(uniq_times)
                v = ds.createVariable(vname, "f4", ("time" ,"y" ,"x"), zlib=True, complevel=4, fill_value=np.nan)
                v.coordinates = "lat lon"
                v.long_name = attrs.get("long_name", vname.replace("_" ," "))
                if "units" in attrs: v.units = attrs["units"]
                # Initialize with NaN
                v[:, :, :] = np.nan
                # Insert slices
                for t, arr in zip(times, data_stack):
                    idx = time_index[t]
                    v[idx, :, :] = arr
        # Global attributes
        ds.Conventions = "CF-1.8"
        if global_attrs:
            for k, v in global_attrs.items():
                setattr(ds, k, v)
    finally:
        ds.close()

def read_netcdf(path: str | Path) -> Dict[str, Any]:
    """Read NetCDF metadata into a dict for quick inspection."""
    Dataset, _, _ = _nc4()
    path = Path(path)
    out = {"path": str(path), "dims": {}, "vars": {}, "global_attrs": {}}
    with Dataset(path, "r") as ds:
        # Dims
        for k, dim in ds.dimensions.items():
            out["dims"][k] = len(dim) if not dim.isunlimited() else "unlimited"
        # Vars
        for k, v in ds.variables.items():
            attrs = {a: getattr(v, a) for a in v.ncattrs()}
            out["vars"][k] = {"dtype": str(v.dtype), "shape": v.shape, "attrs": attrs}
        # Globals
        for a in ds.ncattrs():
            out["global_attrs"][a] = getattr(ds, a)
    return out

def read_nc_var2d(nc_path: Path, varname: str, shape_hint: tuple[int,int]) -> np.ndarray | None:
    if varname is None:
        return None
    with Dataset(nc_path, "r") as ds:
        if varname not in ds.variables:
            return None
        arr = ds.variables[varname][:]
    arr = np.asarray(arr).squeeze()
    if arr.ndim != 2:
        raise ValueError(f"NC var '{varname}' is not 2D after squeeze (ndim={arr.ndim}).")
    if arr.shape != shape_hint:
        raise ValueError(f"NC var '{varname}' shape {arr.shape} != {shape_hint}.")
    return arr


# ----------------------------------------------------
# Collection helpers
# ----------------------------------------------------

def collect_time_series_from_folder(folder: str | Path, pattern: str, var_name: str, units: str = "1") -> Tuple[
    str, List[datetime], np.ndarray, Dict[str, Any]]:
    """Collect time-varying .asc files from a folder matching a regex pattern,
    parse time from filenames, and stack into (T, y, x).
    Returns (var_name, times_sorted, stack, attrs).
    """
    folder = Path(folder)
    files = [p for p in folder.iterdir() if p.is_file() and re.search(pattern, p.name)]
    records = []
    header_template = None
    for p in files:
        t = guess_time_from_filename(p.name)
        if t is None:
            continue
        asc = read_esri_ascii(p)
        if header_template is None:
            header_template = asc["header"]
        records.append((t, asc["data"]))
    if not records:
        raise ValueError(f"No time-stamped ASC files matched pattern: {pattern}")
    # Ensure consistent shapes
    shapes = {r[1].shape for r in records}
    if len(shapes) != 1:
        raise ValueError(f"Inconsistent shapes among ASC files: {shapes}")
    times = sorted([r[0] for r in records])
    stack = np.stack([r[1] for r in sorted(records, key=lambda x: x[0])], axis=0)  # (T, y, x)
    attrs = {"long_name": var_name.replace("_", " "), "units": units}
    return var_name, times, stack, attrs


def collect_static_from_list(files: List[str | Path], names: Optional[List[str]] = None,
                             units: Optional[List[str]] = None, long_names: Optional[List[str]] = None) -> Dict[
    str, Tuple[np.ndarray, Dict[str, Any]]]:
    """Collect static ASC files into a dict suitable for write_netcdf()."""
    out = {}
    for i, f in enumerate(files):
        asc = read_esri_ascii(f)
        arr = asc["data"]
        name = names[i] if names and i < len(names) else Path(f).stem.lower()
        attrs = {}
        if long_names and i < len(long_names): attrs["long_name"] = long_names[i]
        if units and i < len(units): attrs["units"] = units[i]
        out[name] = (arr, attrs)
    return out

def month_iter():
    for m in range(1, 13):
        yield m


# ----------------------------------------------------
# output
# ----------------------------------------------------

def export_var_to_asc_series(nc_path: str | Path, var_name: str, template_asc: str | Path, out_dir: str | Path,
                             filename_fmt: str = "{var}_{YYYY}-{MM}.asc") -> List[Path]:
    """Export a (time,y,x) variable to a series of ASC files using the header from template_asc.

    filename_fmt supports tokens: {var}, {YYYY}, {MM}, {DD}, {index}
    """
    out_paths = []
    out_dir = Path(out_dir);
    out_dir.mkdir(parents=True, exist_ok=True)
    # Load template header
    tmpl = read_esri_ascii(template_asc)["header"]
    with Dataset(nc_path, "r") as ds:
        if var_name not in ds.variables:
            raise KeyError(f"Variable '{var_name}' not found in {nc_path}")
        v = ds.variables[var_name]
        if v.ndim != 3 or v.shape[0] != len(ds.dimensions["time"]):
            raise ValueError(f"Variable '{var_name}' must be (time,y,x)")
        times_num = ds.variables["time"][:]
        units = ds.variables["time"].units
        calendar = getattr(ds.variables["time"], "calendar", "proleptic_gregorian")
        times = num2date(times_num, units=units, calendar=calendar)
        for i, t in enumerate(times):
            arr = v[i, :, :].filled(np.nan) if hasattr(v[i, :, :], "filled") else v[i, :, :]
            YYYY = f"{t.year:04d}";
            MM = f"{t.month:02d}";
            DD = f"{t.day:02d}"
            name = filename_fmt.format(var=var_name, YYYY=YYYY, MM=MM, DD=DD, index=i)
            path = out_dir / name
            write_esri_ascii(path, np.array(arr), tmpl)
            out_paths.append(path)
    return out_paths


def read_asc_header_and_dims(example_asc: Path) -> tuple[str, int, int]:
    """
    Read the first 6 lines of an ASC template EwE Ecospace file to
    capture the header and
    parse ncols/nrows. Returns (header_text, nrows, ncols).
    """
    with open(example_asc, 'r') as f:
        header_lines = [next(f).rstrip("\n") for _ in range(6)]
    header = "\n".join(header_lines)

    # Extract ncols / nrows
    def _find_int(key: str) -> int:
        for line in header_lines:
            if line.lower().startswith(key):
                # Split by whitespace and take the last token
                parts = line.split()
                return int(parts[-1])
        raise ValueError(f"Could not find {key} in ASC header:\n{header}")

    ncols = _find_int("ncols")
    nrows = _find_int("nrows")

    return header, nrows, ncols

def build_output_name(base_name: str, year: int, month: int) -> str:
    return f"{base_name}_{year}_{month:02d}.asc"

# ----------------------------------------------------
# Masks
# ----------------------------------------------------

def build_inclusion_mask(nrows: int, ncols: int,
                         dfPlumeMask: pd.DataFrame | None,
                         dfLandMask: pd.DataFrame | None) -> np.ndarray:
    """Return a boolean mask of cells that should *participate* in normalization.
    True = included; False = excluded (land or outside allowed region).

    We interpret:
      - dfLandMask == True where land -> EXCLUDE
      - dfPlumeMask == True where inside plume region -> INCLUDE
    """
    include = np.ones((nrows, ncols), dtype=bool)

    if dfLandMask is not None:
        land_bool = dfLandMask.values if hasattr(dfLandMask, "values") else np.asarray(dfLandMask)
        if land_bool.shape == include.shape:
            include &= ~land_bool
        else:
            print("[WARN] Land mask shape mismatch; skipping land mask in normalization.")

    if dfPlumeMask is not None:
        plume_bool = dfPlumeMask.values if hasattr(dfPlumeMask, "values") else np.asarray(dfPlumeMask)
        if plume_bool.shape == include.shape:
            include &= plume_bool
        else:
            print("[WARN] Plume mask shape mismatch; skipping plume mask in normalization.")

    return include



# ---------------------------
# IO utilities
# ---------------------------

def read_asc_array(path: Path, nrows: int, ncols: int) -> np.ndarray:
    arr = np.loadtxt(path, skiprows=6)
    if arr.shape != (nrows, ncols):
        raise ValueError(f"ASC shape {arr.shape} != {(nrows, ncols)} for {path}")
    return arr

def load_mask_asc(mask_path: Path, nrows: int, ncols: int, nodata=None):
    """Return boolean water mask from an ASC where 0 = land."""
    arr = np.loadtxt(mask_path, skiprows=6)
    if arr.shape != (nrows, ncols):
        raise ValueError(f"Mask shape {arr.shape} != {(nrows, ncols)}")
    if nodata is not None:
        arr = np.where(arr == nodata, 0, arr)
    water = (arr != 0).astype(bool)
    return water


def load_masks_from_nc(nc_path: Path,
                       depth_var: str = "depth",
                       rivers_var: str | None = None,
                       offshore_var: str | None = None):
    with Dataset(nc_path, "r") as ds:
        depth = np.array(ds.variables[depth_var][:]).squeeze()
        if depth.ndim != 2:
            raise ValueError(f"{depth_var} must be 2D; got {depth.shape}")
        water_mask = depth > 0

        def _load_bool(name):
            if name is None: return None
            if name not in ds.variables:
                raise KeyError(f"'{name}' not found in {nc_path}")
            arr = np.array(ds.variables[name][:]).squeeze()
            if arr.ndim != 2:
                raise ValueError(f"{name} must be 2D; got {arr.shape}")
            return (arr == 1)

        rivers_mask = _load_bool(rivers_var)
        offshore_mask = _load_bool(offshore_var)

    return depth, water_mask, rivers_mask, offshore_mask


def resolve_value_column(df: pd.DataFrame, explicit: str | None) -> str:
    """
    Pick a value/intensity column.
    If 'explicit' provided, use it (and error if missing).
    Otherwise try a few common names.
    """
    if explicit is not None:
        if explicit in df.columns:
            return explicit
        else:
            raise KeyError(f"Requested value column '{explicit}' not in CSV. Columns: {list(df.columns)}")

    candidates = [
        "SpawnIntensity", "spawn_intensity", "Intensity", "intensity",
        "Eggs", "EggDensity", "Spawn Index", "Spawn_Index"
    ]
    for c in candidates:
        if c in df.columns:
            return c

    # Fallback: if no candidates, error
    raise KeyError("Could not locate a spawn intensity/value column. "
                   "Use --value_col to specify the desired column.")



# ===== read herring time series from NC (time,y,x) and year/month coords =====
def collect_herring_from_nc(nc_path: Path,
                            var_name: str,
                            nrows: int,
                            ncols: int,
                            year_var: str | None,
                            month_var: str | None,
                            time_var: str | None):
    with Dataset(nc_path, "r") as ds:
        # 1) Turn OFF auto mask/scale so zeros stay zeros (no masked→NaN surprise)
        ds.set_auto_maskandscale(False)

        if var_name not in ds.variables:
            raise KeyError(f"NC missing var '{var_name}'")

        v = ds.variables[var_name]

        data = np.asarray(v[:], dtype=float)  # plain ndarray, may contain NaNs
        data = np.nan_to_num(data, nan=0.0)  # <-- key line: NaN → 0.0

        if data.ndim != 3:
            raise ValueError(f"NC var '{var_name}' must be 3D (time,y,x); got {data.ndim}D")
        t, y, x = data.shape
        if (y, x) != (nrows, ncols):
            raise ValueError(f"'{var_name}' shape {(y,x)} != {(nrows,ncols)}")

        # ---- DEBUG: tell us if it *looks* empty/masked ----
        fv = getattr(v, "_FillValue", None)
        mv = getattr(v, "missing_value", None)
        print(f"[DEBUG] {var_name}: shape={data.shape}, dtype={v.dtype}, "
              f"_FillValue={fv}, missing_value={mv}, "
              f"min={np.nanmin(data):.6g}, max={np.nanmax(data):.6g}, "
              f"nan%={(np.isnan(data).sum()/data.size)*100:.2f}%")

        # 2) If you actually *want* masking applied but without NaNs, do this instead:
        # ds.set_auto_maskandscale(True)
        # data = v[:]
        # if isinstance(data, np.ma.MaskedArray):
        #     print(f"[DEBUG] masked%={(data.mask.sum()/data.size)*100:.2f}%")
        #     data = data.filled(0.0)  # <- choose 0.0 so "no spawn" stays 0, not NaN

        # years/months
        yy = None; mm = None
        if year_var and (year_var in ds.variables):
            yy = np.asarray(ds.variables[year_var][:], dtype=int).ravel()
        if month_var and (month_var in ds.variables):
            mm = np.asarray(ds.variables[month_var][:], dtype=int).ravel()

        # Fallback via 'time' coord
        if (yy is None or mm is None) and time_var and (time_var in ds.variables) and (num2date is not None):
            tv = ds.variables[time_var]
            try:
                dts = num2date(tv[:], units=getattr(tv, "units", None))
                yy = np.array([int(getattr(d, "year", 0)) for d in dts], dtype=int)
                mm = np.array([int(getattr(d, "month", 0)) for d in dts], dtype=int)
            except Exception:
                pass

        if yy is None or mm is None or len(yy) != t or len(mm) != t:
            print("[WARN] Could not decode year/month from NC; fabricating sequence.")
            mm = np.array([(i % 12) + 1 for i in range(t)], dtype=int)
            yy = np.zeros((t,), dtype=int)

        items = []
        for i in range(t):
            # ensure plain float array (no masked dtype sneaking through)
            items.append((int(yy[i]), int(mm[i]), data[i, :, :].astype(float)))
        return items


def save_asc_plain(path: Path, grid: np.ndarray, header: str, fmt: str):
    with open(path, "w", newline="\n") as f:
        f.write(header.rstrip("\r\n") + "\n")
        np.savetxt(f, grid, fmt=fmt, newline="\n")

# ----------------------------------------------------
# parsing related
# ----------------------------------------------------

def parse_year_month(df: pd.DataFrame, date_col: str | None) -> pd.DataFrame:
    """
    Ensure df has 'Year' and 'Month' integer columns.
    If a 'Year' column already exists and looks numeric, keep it.
    If 'Month' exists, keep it.
    Otherwise we parse from 'date_col' (ISO-like) into Year/Month.
    """
    out = df.copy()

    # 1) Year
    if "Year" in out.columns:
        # coerce to numeric and keep
        out["Year"] = pd.to_numeric(out["Year"], errors="coerce").astype("Int64")
    else:
        if date_col is None:
            raise KeyError("No 'Year' column and no --date_col provided to derive Year/Month.")
        # parse Year from date_col
        out[date_col] = pd.to_datetime(out[date_col], errors="coerce", utc=False)
        if out[date_col].isna().all():
            raise ValueError(f"Could not parse any dates from column '{date_col}'.")
        out["Year"] = out[date_col].dt.year

    # 2) Month
    if "Month" in out.columns:
        out["Month"] = pd.to_numeric(out["Month"], errors="coerce").astype("Int64")
    else:
        if date_col is None and "StartDate" in out.columns:
            date_col = "StartDate"
        if date_col is None:
            raise KeyError("No 'Month' column and no --date_col provided to derive Month.")
        if not np.issubdtype(out[date_col].dtype, np.datetime64):
            out[date_col] = pd.to_datetime(out[date_col], errors="coerce", utc=False)
        out["Month"] = out[date_col].dt.month

    return out


def parse_depth_bands(s: str | None) -> list[tuple[float, float]]:
    """Parse '0-15,25-45' -> [(0.0,15.0),(25.0,45.0)]."""
    if not s:
        return []
    bands = []
    for part in s.split(","):
        part = part.strip()
        if "-" not in part:
            raise ValueError(f"Depth band '{part}' must look like 'a-b'")
        a, b = map(float, part.split("-", 1))
        bands.append((min(a, b), max(a, b)))
    return bands


def parse_months(s: str | None) -> set[int]:
    if not s or s.lower() == "none":
        return set()
    if s.lower() == "all":
        return set(range(1, 13))
    out = set()
    for part in s.split(","):
        part = part.strip()
        if "-" in part:
            a, b = map(int, part.split("-", 1))
            out.update(range(min(a, b), max(a, b) + 1))
        else:
            out.add(int(part))
    return {m for m in out if 1 <= m <= 12}

# ===== NEW: parse "1,2,6-9" → [1,2,6,7,8,9] =====
def parse_months_str(s: str | None):
    if not s or s.lower() in ("none", "all"):
        return None if not s or s.lower() == "all" else []
    parts = []
    for token in s.split(","):
        token = token.strip()
        if "-" in token:
            a, b = token.split("-")
            parts.extend(range(int(a), int(b) + 1))
        elif token:
            parts.append(int(token))
    return sorted(set([m for m in parts if 1 <= m <= 12]))

# ----------------------------------------------------
# cleaning
# ----------------------------------------------------

def clean_numeric_column(df: pd.DataFrame, col: str, verbose: bool = True) -> pd.Series:
    """
    Coerce a column to float, handling common non-numeric tokens gracefully.
    - Strips whitespace and thousands separators (",").
    - Converts Excel-style error tokens (e.g., '#VALUE!', '#DIV/0!', '#N/A', '#REF!', '#NAME?') to NaN.
    - Converts parenthesized negatives like '(123)' to '-123'.
    - Coerces to numeric; NaNs are set to 0.0.
    """
    s = df[col]
    if pd.api.types.is_numeric_dtype(s):
        clean = s.astype(float)
        n_bad = 0
    else:
        s_str = s.astype(str).str.strip()
        s_str = s_str.str.replace(',', '', regex=False)                    # remove thousands separators
        s_str = s_str.str.replace(r'^\((.*)\)$', r'-\1', regex=True)       # (123) -> -123
        bad_tokens = {
            '#VALUE!': np.nan, '#DIV/0!': np.nan, '#N/A': np.nan,
            '#REF!': np.nan, '#NAME?': np.nan, 'nan': np.nan,
            'NaN': np.nan, 'None': np.nan, '': np.nan, 'N/A': np.nan,
        }
        s_str = s_str.replace(bad_tokens)
        clean = pd.to_numeric(s_str, errors='coerce')
        n_bad = int(clean.isna().sum())

    clean = clean.replace([np.inf, -np.inf], np.nan)
    if n_bad > 0 and verbose:
        print(f"[WARN] Column '{col}': {n_bad} non-numeric/invalid values coerced to 0.0")
    return clean.fillna(0.0).astype(float)


# ----------------------------------------------------
# weights
# ----------------------------------------------------

# Instead of “overwrite to 1−p” (salmon did that by design),
# use multiplicative boosts for herring (gentler, plays nicely with normalisation
def compose_monthly_weights(nrows, ncols,
                            water_mask: np.ndarray,
                            rivers_mask: np.ndarray | None,
                            offshore_mask: np.ndarray | None,
                            rivers_active: bool,
                            offshore_active: bool,
                            rivers_w: float = 1.5,
                            offshore_w: float = 1.2):
    W = np.ones((nrows, ncols), float)
    W[~water_mask] = 0.0  # keep land as 0 to neutralize any later multiply
    if rivers_mask is not None and rivers_active:
        W[rivers_mask & water_mask] *= rivers_w
    if offshore_mask is not None and offshore_active:
        W[offshore_mask & water_mask] *= offshore_w
    return W


# A smooth triangular band (peak weight at center of the band, tapering to 0 at edges) works
# well and avoids discontinuities:
def make_depth_weight(depth_2d: np.ndarray,
                      bands: list[tuple[float, float]] = [(5.0, 40.0)],
                      peak_weight: float = 2.0,
                      floor_weight: float = 1.0):
    """
    Returns a weight grid >= floor_weight. For each band [zmin, zmax],
    apply a triangular preference peaking midway; combine by max across bands.
    """
    h, w = depth_2d.shape
    out = np.full((h, w), floor_weight, float)
    for (z0, z1) in bands:
        if z1 <= z0: continue
        mid = 0.5 * (z0 + z1)
        halfw = 0.5 * (z1 - z0)
        # triangular shape: 0 at edges, 1 at mid; then scale to peak_weight
        tri = np.clip(1.0 - np.abs(depth_2d - mid) / max(halfw, 1e-6), 0.0, 1.0)
        band_w = floor_weight + (peak_weight - floor_weight) * tri
        out = np.maximum(out, band_w)
    return out


def depth_pref_from_depth(depth_m: np.ndarray,
                           include: np.ndarray,
                           mode: str,
                           dmin: float,
                           dmax: float,
                           soften: float) -> np.ndarray:
    """
    Build a 0..1 preference from raw depth (meters, positive).
    """
    pref = np.zeros_like(depth_m, dtype=float)
    d = depth_m.astype(float)
    pref[~include] = 0.0

    if mode == "scale":
        m = d[include].max()
        if m > 0:
            pref[include] = d[include] / m

    elif mode == "invert":
        m = d[include].max()
        if m > 0:
            tmp = d / m
            pref[include] = 1.0 - tmp[include]

    elif mode == "range_binary":
        inside = include & (d >= dmin) & (d <= dmax)
        pref[inside] = 1.0

    elif mode == "range_trapezoid":
        s = max(float(soften), 0.0)
        core = (d >= dmin) & (d <= dmax)
        left  = (d >= dmin - s) & (d < dmin) if s > 0 else np.zeros_like(core)
        right = (d >  dmax) & (d <= dmax + s) if s > 0 else np.zeros_like(core)
        pref[core & include] = 1.0
        if s > 0:
            pref[left & include] = (d[left & include] - (dmin - s)) / s
            pref[right & include] = ((dmax + s) - d[right & include]) / s
        pref[pref < 0] = 0.0
        pref[pref > 1] = 1.0

    else:
        raise ValueError(f"Unknown DEF_DEPTH_MODE: {mode}")

    return pref



# ------------------------------
# Gaussian filter
# ------------------------------

def smooth_gaussian_masked(arr: np.ndarray, mask: np.ndarray, sigma: float) -> np.ndarray:
    """
    Gaussian blur only within mask. Values near land are renormalized so they
    don't get artificially depressed by zeros outside the mask.
    """
    if sigma is None or sigma <= 0:
        return arr
    if gaussian_filter is None:
        raise RuntimeError("scipy not available; use --box_radius instead or install scipy")

    arr_masked = arr * mask
    blurred_num = gaussian_filter(arr_masked, sigma=sigma, mode="nearest")
    blurred_den = gaussian_filter(mask.astype(float), sigma=sigma, mode="nearest")
    # Avoid divide by zero; keep original where denominator tiny
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(blurred_den > 1e-9, blurred_num / blurred_den, arr)
    return out


def gaussian_blur_masked(arr: np.ndarray, mask: np.ndarray, sigma: float) -> np.ndarray:
    if sigma is None or sigma <= 0:
        return arr
    if gaussian_filter is None:
        print("[WARN] scipy not available; skipping gaussian blur.")
        return arr
    # blur value*mask and mask, then divide
    val = gaussian_filter(arr * mask, sigma=sigma, mode="nearest")
    den = gaussian_filter(mask.astype(float), sigma=sigma, mode="nearest")
    out = np.zeros_like(arr)
    good = den > 1e-12
    out[good] = val[good] / den[good]
    return out

# ------------------------------
# math, vals, floor, normalise etc
# ------------------------------

def norm_array(a: np.ndarray, include: np.ndarray, method: str, min_floor: float) -> np.ndarray:
    out = a.astype(float).copy()
    out[~include] = 0.0
    if method == "sum":
        s = out[include].sum()
        if s > 0:
            out[include] /= s
    elif method == "max":
        m = out[include].max()
        if m > 0:
            out[include] /= m
    elif method == "none":
        pass
    else:
        raise ValueError(f"Unknown norm method: {method}")
    if min_floor > 0:
        nz = include & (out > 0)
        out[nz] = np.maximum(out[nz], min_floor)
        if method == "sum":
            s = out[include].sum()
            if s > 0:
                out[include] /= s
        elif method == "max":
            m = out[include].max()
            if m > 0:
                out[include] /= m
    return out

def scale_to_unit_by_max(a: np.ndarray, include: np.ndarray, invert: bool = False) -> np.ndarray:
    out = a.astype(float).copy()
    out[~include] = 0.0
    m = out[include].max()
    if m > 0:
        out[include] /= m
    if invert:
        out[include] = 1.0 - out[include]
    return out


def apply_min_floor(grid: np.ndarray,
                    inclusion: np.ndarray,
                    min_floor: float,
                    renormalize: bool = False,
                    target_sum: float = 1.0,
                    eps: float = 1e-12) -> np.ndarray:
    """
    Enforce a minimum value on all INCLUDED cells. Optionally renormalize to a target sum
    (e.g., 1.0). Guarantees no included cell < min_floor. If renormalization is requested
    and the requested min_floor is infeasible (min_floor * N > target_sum), it auto-clips
    the floor to target_sum / N.

    - grid, inclusion: same shape; inclusion=True means 'water/included' cell
    - min_floor: floor to enforce on inclusion cells
    - renormalize=False is recommended for 'presence probability' outputs
    """
    G = grid.copy()
    inc = inclusion.astype(bool)

    if min_floor <= 0:
        # nothing to do
        return G

    # Feasibility check (only matters if we renormalize)
    if renormalize:
        N = int(inc.sum())
        if N > 0 and min_floor * N > target_sum:
            # Auto-clip the floor to the maximum feasible
            min_floor = target_sum / N

    # Raise all included cells below floor (including zeros)
    low_mask = inc & (G < min_floor)
    G[low_mask] = min_floor

    if not renormalize:
        return G

    # Keep the sum constraint while honoring the floor
    s = G[inc].sum()
    if s <= eps:
        return G

    excess = s - target_sum  # >0 means we need to subtract from some cells
    if abs(excess) <= eps:
        return G

    # Only redistribute among cells strictly above floor
    hi_mask = inc & (G > min_floor + eps)
    if not np.any(hi_mask):
        # No headroom to redistribute; fallback to uniform feasible vector
        G[inc] = target_sum / float(inc.sum())
        return G

    headroom = G[hi_mask] - min_floor
    cap_sum = headroom.sum()
    if cap_sum <= eps:
        G[inc] = target_sum / float(inc.sum())
        return G

    # Proportional redistribution
    G[hi_mask] -= excess * (headroom / cap_sum)

    # Final clean-up and tiny rescale to hit target exactly
    G[inc] = np.maximum(G[inc], min_floor)
    s2 = G[inc].sum()
    if s2 > eps:
        G[inc] *= (target_sum / s2)
    return G


def _scale_to_unit_interval(arr: np.ndarray,
                            mask: np.ndarray,
                            log_norm: bool = False,
                            log_eps: float = 1e-6,
                            quantile: float | None = None) -> np.ndarray:
    """
    Return a 0..1 suitability from a non-negative 'arr' using either max-scaling or
    quantile scaling over masked (included) positives. Optional log compression first.
    """
    A = arr.copy()
    A[~mask] = 0.0

    if log_norm:
        pos = A > 0
        if np.any(pos):
            A[pos] = np.log(A[pos] + log_eps)
            # shift to start at zero over included cells
            mn = A[mask].min(initial=0.0)
            if mn < 0:
                A[mask] -= mn
        else:
            return np.zeros_like(A)

    posvals = A[mask & (A > 0)]
    if posvals.size == 0:
        return np.zeros_like(A)

    if quantile is not None:
        scale = float(np.quantile(posvals, quantile))
    else:
        scale = float(posvals.max())

    if scale <= 0:
        return np.zeros_like(A)

    S = np.clip(A / scale, 0.0, 1.0)
    S[~mask] = 0.0
    return S


def normalize_to_one(grid: np.ndarray,
                     include: np.ndarray,
                     min_floor: float = 0.0) -> np.ndarray:
    out = grid.astype(float).copy()
    out[~include] = 0.0
    s = out[include].sum()
    if s > 0:
        out[include] = out[include] / s
    # Apply floor on water cells, then renormalize to exactly 1
    if min_floor > 0:
        out[include] = np.maximum(out[include], min_floor)
        s2 = out[include].sum()
        if s2 > 0:
            out[include] = out[include] / s2
    return out


def _normalize_grid(grid: np.ndarray, include: np.ndarray, method: str, min_floor: float) -> np.ndarray:
    out = grid.astype(float).copy()
    out[~include] = 0.0

    if method == "sum":
        s = out[include].sum()
        if s > 0:
            out[include] /= s
    elif method == "max":
        m = out[include].max()
        if m > 0:
            out[include] /= m
    else:
        raise ValueError(f"Unknown norm method: {method}")

    # Impose floor on any non-zero water cell, then re-enforce the chosen normalization
    if min_floor > 0:
        mask = include & (out > 0)
        out[mask] = np.maximum(out[mask], min_floor)
        if method == "sum":
            s = out[include].sum()
            if s > 0:
                out[include] /= s
        else:  # "max"
            m = out[include].max()
            if m > 0:
                out[include] /= m

    return out