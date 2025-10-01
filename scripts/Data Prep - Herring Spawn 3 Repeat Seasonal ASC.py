"""
Repeat January Map by Year (ASC)

Create by G Oldford Sep 2025

Purpose
-------
Given a single January ASC map (e.g., HERRING_SPAWN_PPROB_01.asc), create
year-specific copies for a specified range:
    HERRING_SPAWN_PPROB_01_1950.asc
    HERRING_SPAWN_PPROB_01_1951.asc
    ...
    HERRING_SPAWN_PPROB_01_2024.asc

Usage
-----
python repeat_january_map_by_year.py \
    --jan_file /path/to/HERRING_SPAWN_PPROB_01.asc \
    --start_year 1950 --end_year 2024 \
    --out_dir /path/to/output

Optional:
- You can override the base name used for outputs with --basename.
  If not provided, the source file's stem is used (without extension).
- Use --overwrite to replace existing outputs.
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

# ------------------------------
# paths and default params
# ------------------------------
DEFAULT_SPAWN_P    = Path("C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Herring_Spawning_2023//MODIFIED//")
DEF_OUT_P = Path(DEFAULT_SPAWN_P, "ECOSPACE_ASC")
DEF_FILE_PATTERN = 'HERRING_SPAWN_PPROB_{mm}.asc'
DEF_BASE_NAME = 'HERRING_SPAWN_PPROB'
DEF_START_YEAR = 1950
DEF_END_YEAR = 2024
DEF_ORDER = "year_month"
DEF_MONTHS = '1-12'
DEF_YR_SUBDIRS = False

def parse_months(spec: str) -> list[int]:
    months = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^(\d{1,2})-(\d{1,2})$", part)
        if m:
            a = int(m.group(1)); b = int(m.group(2))
            if a > b:
                raise ValueError(f"Invalid months range: {part}")
            months.extend(range(a, b + 1))
        else:
            months.append(int(part))
    months = sorted(set(months))
    for m in months:
        if m < 1 or m > 12:
            raise ValueError(f"Month out of range: {m}")
    return months


def resolve_month_file(input_dir: Path, pattern: str, month: int) -> Path:
    if "{mm}" not in pattern:
        raise ValueError("Pattern must include '{mm}', e.g., 'HERRING_SPAWN_PPROB_{mm}.asc'")
    return (input_dir / pattern.replace("{mm}", f"{month:02d}")).resolve()


def make_out_path(base_dir: Path, basename: str, year: int, month: int, order: str, year_subdirs: bool, suffix: str) -> Path:
    fname = (
        f"{basename}_{month:02d}_{year}{suffix}" if order == "month_year"
        else f"{basename}_{year}_{month:02d}{suffix}"
    )
    return (base_dir / str(year) / fname) if year_subdirs else (base_dir / fname)


def main():
    ap = argparse.ArgumentParser(description="Copy the 12 monthly ASC files across a year range with simple renaming.")
    ap.add_argument("--input_dir", default=DEF_OUT_P, type=Path, help="Directory containing the 12 monthly inputs.")
    ap.add_argument("--pattern", default=DEF_FILE_PATTERN, type=str, help="Filename pattern with {mm}, e.g., 'HERRING_SPAWN_PPROB_{mm}.asc'")
    ap.add_argument("--basename", default=DEF_BASE_NAME, type=str, help="Base name for output files (e.g., 'HERRING_SPAWN_PPROB')")
    ap.add_argument("--start_year", default=DEF_START_YEAR, type=int, help="Start year (inclusive).")
    ap.add_argument("--end_year", default=DEF_END_YEAR, type=int, help="End year (inclusive).")
    ap.add_argument("--out_dir", default=DEF_OUT_P, type=Path, help="Output directory.")
    ap.add_argument("--order", default=DEF_ORDER, type=str, choices=["month_year", "year_month"], help="Output filename order.")
    ap.add_argument("--months", default=DEF_MONTHS, type=str, help="Months to include, e.g., '1-12' or '3,4,5' or '1,2,7-9'.")
    ap.add_argument("--year_subdirs", action="store_true", default=DEF_YR_SUBDIRS, help="Create per-year subfolders.")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    args = ap.parse_args()

    if args.end_year < args.start_year:
        ap.error("end_year must be >= start_year")

    months = parse_months(args.months)

    input_dir = args.input_dir.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Determine suffix from first month file (keeps '.asc' or any extension)
    first_src = resolve_month_file(input_dir, args.pattern, months[0])
    if not first_src.exists():
        ap.error(f"Missing monthly file for month={months[0]:02d}: {first_src}")
    suffix = first_src.suffix

    # Verify all monthly sources exist
    for m in months:
        src = resolve_month_file(input_dir, args.pattern, m)
        if not src.exists():
            ap.error(f"Missing monthly file for month={m:02d}: {src}")

    total = 0
    for year in range(args.start_year, args.end_year + 1):
        for m in months:
            src = resolve_month_file(input_dir, args.pattern, m)
            dst = make_out_path(out_dir, args.basename, year, m, args.order, args.year_subdirs, suffix)
            dst.parent.mkdir(parents=True, exist_ok=True)

            if dst.exists() and not args.overwrite:
                print(f"[SKIP] Exists: {dst}")
                continue

            shutil.copy2(src, dst)
            print(f"[OK] {dst}")
            total += 1

    print(f"[INFO] Wrote {total} file(s) to {out_dir}")


if __name__ == "__main__":
    main()