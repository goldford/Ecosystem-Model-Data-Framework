"""
G Oldford Sep 2025

Purpose:
- Load ASC files into a single NetCDF for storage

Process:
- Reads grid CSV for lat/lon
- Loads static ASC layers
- Loads time-varying seal intensity ASCs (expects YYYY-M in filename)
- Writes a single NetCDF bundle with CF-1.8 coords
- Prints a quick summary of the resulting file


"""

from pathlib import Path
from glob import glob
from helpers_ewe_asc_netcdf import (
    build_latlon_from_rowscols_csv,
    collect_static_from_list,
    collect_time_series_from_folder,
    write_netcdf,
    read_netcdf,
)

PATH1 = Path("..//data//basemap//")
PATH2 = Path("..//..//ecospace_seal_diets//Data//basemaps//")
PATH3 = Path("..//..//ecospace_seal_diets//Data//seal_haulouts//model_out//")
PATH4 = Path("..//..//ecospace_seal_diets//Data//basemaps//")
PATH5 = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//SPAWN_DATA_ASC//")
PATH6 = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//HAB_CAP_ASC//")


GRID_CSV = PATH1 / "Ecospace_grid_20210208_rowscols.csv"

STATIC_FILES = [
    PATH1 / "Fraser_Plume_Region.asc",
    PATH1 / "ecospacedepthgrid.asc",
    PATH2 / "offshore_areas_map.asc",
    PATH2 / "offshore_areas_map_chum.asc",
    PATH2 / "offshore_areas_map_herring.asc",
    PATH2 / "chum_rivers.asc",
    PATH2 / "fraser_inland_mask.asc",
    PATH2 / "fraser_boundary_nospawn_areas_herring.asc"
]
STATIC_NAMES = ["fraser_plume_region",
                "depth", "offshore_areas",
                "offshore_areas_chum",
                "offshore_areas_herring",
                "chum_rivers",
                "fraser_inland_mask",
                "no_spawn_shallows_herring"
                ]

STATIC_UNITS = ["1", "m", "1", "1"]
STATIC_LN = ["Fraser Plume Region",
             "Ecospace depth",
             "Offshore areas mask",
             "Offshore areas Chum mask",
             "Offshore areas herring mask",
             "Chum rivers mask",
             "Fraser inland mask",
             "shallows with limited herring spawn"
             ]

# --- Add HERRING BASE (all-years, all-months) as a static layer ---
HERRING_BASE = PATH5 / "HERRING_SPAWN_BASE_ALLYEARS.asc"
if HERRING_BASE.exists():
    STATIC_FILES.append(HERRING_BASE)
    STATIC_NAMES.append("herring_spawn_base_all_years")
    STATIC_UNITS.append("1")
    STATIC_LN.append("Herring spawn (across all years & months)")
else:
    print("[WARN] Missing HERRING_SPAWN_BASE_ALLYEARS.asc in", PATH5)

# --- Add Herring Habcap
HERRING_HABCAP = PATH6 / "HERRING_HABCAP_BASE_ALLYEARS.asc"
if HERRING_HABCAP.exists():
    STATIC_FILES.append(HERRING_HABCAP)
    STATIC_NAMES.append("herring_habcap_base_all_years")
    STATIC_UNITS.append("1")
    STATIC_LN.append("Herring hab cap (across all years & months)")
else:
    print("[WARN] Missing HERRING_HABCAP_BASE_ALLYEARS.asc in", PATH6)

# --- Add HERRING SPAWN monthly climatology (CLIM_01..12) as static layers ---
clim_files = sorted(glob(str(PATH5 / "HERRING_SPAWN_CLIM_*.asc")))
if clim_files:
    for fp in clim_files:
        fn = Path(fp).name
        # Expect name like HERRING_SPAWN_CLIM_01.asc
        mm = fn.split("_")[-1].split(".")[0]  # "01"
        STATIC_FILES.append(Path(fp))
        STATIC_NAMES.append(f"herring_spawn_climatology_m{mm}")
        STATIC_UNITS.append("1")
        STATIC_LN.append(f"Herring spawn monthly climatology (m={mm})")
else:
    print("[WARN] No CLIM files found in", PATH5)

# --- Add HERRING HABCAP monthly climatology (CLIM_01..12) as static layers ---
clim_files = sorted(glob(str(PATH6 / "HERRING_HABCAP_CLIM_*.asc")))
if clim_files:
    for fp in clim_files:
        fn = Path(fp).name
        # Expect name like HERRING_HABCAP_CLIM_01.asc
        mm = fn.split("_")[-1].split(".")[0]  # "01"
        STATIC_FILES.append(Path(fp))
        STATIC_NAMES.append(f"herring_habcap_climatology_m{mm}")
        STATIC_UNITS.append("1")
        STATIC_LN.append(f"Herring habcap monthly climatology (m={mm})")
else:
    print("[WARN] No CLIM files found in", PATH6)

# Define multiple time-series sources here. Add/remove blocks as needed.
TS_CONFIGS = [
    {
        "folder": PATH3,
        "pattern": r"sealforagingintens.*\.asc$",     # regex to match files
        "var_name": "seal_foraging_intensity",        # NetCDF variable name
        "units": "1",
        "long_name": "Seal foraging intensity"
    },
    {
        "folder": PATH5,
        "pattern": r"HERRING_SPAWN_\d{4}_\d{2}\.asc$",  # YYYY_MM only
        "var_name": "herring_spawn_intensity",
        "units": "1",
        "long_name": "Herring spawn intensity (per year-month)"
    },
    {
        "folder": PATH6,
        "pattern": r"HERRING_HABCAP_\d{4}_\d{2}\.asc$",  # YYYY_MM only
        "var_name": "herring_habcap",
        "units": "1",
        "long_name": "Herring habcap (per year-month)"
    }
]

OUT_NC = PATH4 / "ecospace_bundled_asc.nc"
# ---------- End of config ----------

def main():

    # 1) Build coords
    lat2d, lon2d, shape, meta = build_latlon_from_rowscols_csv(GRID_CSV)
    print(f"[INFO] Grid shape (y,x): {shape}, one_based={meta['one_based']}")

    # 2) Static layers
    static_layers = collect_static_from_list(
        files=STATIC_FILES,
        names=STATIC_NAMES,
        units=STATIC_UNITS,
        long_names=STATIC_LN,
    )
    print(f"[INFO] Loaded {len(static_layers)} static layers.")

    # 3) Collect multiple time-series
    ts_vars = []
    for cfg in TS_CONFIGS:
        folder = Path(cfg["folder"])
        pattern = cfg["pattern"]
        var_name = cfg["var_name"]
        units = cfg.get("units", "1")
        try:
            ts_tuple = collect_time_series_from_folder(folder=folder, pattern=pattern, var_name=var_name, units=units)
            # Optionally override long_name
            ln = cfg.get("long_name")
            if ln:
                # ts_tuple = (var_name, times, stack, attrs)
                ts_tuple = (ts_tuple[0], ts_tuple[1], ts_tuple[2], {**ts_tuple[3], "long_name": ln})
            ts_vars.append(ts_tuple)
            print(f"[INFO] Collected '{var_name}' from {folder} -> {len(ts_tuple[1])} time steps.")
        except Exception as e:
            print(f"[WARN] Skipping series '{var_name}' from {folder}: {e}")

    if not ts_vars:
        print("[WARN] No time series collected; writing only static layers.")

    # 4) Write NetCDF
    write_netcdf(
        path_out=OUT_NC,
        lat2d=lat2d, lon2d=lon2d,
        static_vars=static_layers,
        time_series_vars=ts_vars if ts_vars else None,
        global_attrs={
            "title": "EwE Ecospace bundle (multiple time-series)",
            "source": "Greig uploads",
            "history": "Created by build_ecospace_nc_multits.py"
        }
    )

    # 5) Inspect
    meta = read_netcdf(OUT_NC)
    print("[INFO] Wrote:", meta["path"])
    print("[INFO] Dims:", meta["dims"])
    print("[INFO] Vars:", list(meta["vars"].keys()))

if __name__ == "__main__":
    main()