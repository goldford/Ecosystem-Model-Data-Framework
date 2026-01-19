"""bt_sweep_4a.py

Run a parameter sweep over bloom-timing settings and re-run the Ecospace
bloom-timing evaluation script (4a) for each combo.

How it works
------------
- Uses environment variables to override config values at import time.
- Spawns a fresh Python process per run (so imports are clean).
- Writes each run's outputs into its own folder, so nothing is overwritten.
- Uses a shared cache folder for:
  - the 2D Suchy/Sat mask (cached per region)
  - the C09 bloom timing CSV (computed once and reused)

Files expected to be in the same directory / import path:
- ecospace_eval_config_sweep.py
- ecospace_eval_4a_bloomt_CSoG_sweep.py
- helpers.py (your existing helpers module)

Usage (from the folder that can already run script 4a):
    python bt_sweep_4a.py

Optional environment variables:
- BT_SWEEP_ROOT: output root folder (default: ..//..//data//evaluation//bt_sweep_4a)
- BT_SWEEP_LIMIT: int; if set, run only the first N combos (handy for smoke tests)
"""

from __future__ import annotations

import itertools
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def _env_bool(v: str | None, default: bool = False) -> bool:
    if v is None:
        return default
    return v.strip().lower() in ("1", "true", "t", "yes", "y")


def _run_tag(
    *,
    region: str,
    vars_to_analyze: tuple[str, ...],
    annual_avg_method: str,
    log_transform: bool,
    mean_or_median: str,
    exclude_dec_jan: bool,
) -> str:
    vars_tag = "PP1" if vars_to_analyze == ("PP1-DIA",) else "PP123"
    avg_tag = "ann" if annual_avg_method == "annual" else "all"
    log_tag = "log" if log_transform else "lin"
    mm_tag = "med" if mean_or_median == "median" else "mean"
    ex_tag = "exdj1" if exclude_dec_jan else "exdj0"
    return f"bt_{region}_{vars_tag}_{avg_tag}_{log_tag}_{mm_tag}_{ex_tag}"


def main() -> int:
    # --- sweep grid ---
    regions = ("SGC2", "SGC3", "SGC4")
    vars_options = (("PP1-DIA", "PP2-NAN", "PP3-PIC"), ("PP1-DIA",))
    annual_avg_methods = ("annual", "all")
    log_transform_opts = (False, True)
    mean_or_median_opts = ("median", "mean")
    exclude_dec_jan_opts = (False, True)
    combos = list(
        itertools.product(
            regions,
            vars_options,
            annual_avg_methods,
            log_transform_opts,
            mean_or_median_opts,
            exclude_dec_jan_opts,
        )
    )

    # Pull scenario code from config (used for shared cache filenames)
    try:
        import ecospace_eval_config_sweep as cfg  # expected to be alongside this script
        scenario = getattr(cfg, 'ECOSPACE_SC', 'SC204')
    except Exception:
        scenario = 'SC204'

    # --- output dirs ---
    sweep_root = Path(
        os.getenv("BT_SWEEP_ROOT", f"..//..//data//evaluation//bt_sweep_4a_{scenario}")
    ).resolve()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    sweep_root = sweep_root / stamp

    cache_dir = sweep_root / "_cache"
    runs_dir = sweep_root / "runs"
    cache_dir.mkdir(parents=True, exist_ok=True)
    runs_dir.mkdir(parents=True, exist_ok=True)

    limit_str = os.getenv("BT_SWEEP_LIMIT")
    limit = int(limit_str) if (limit_str and limit_str.isdigit()) else None

    # Lightweight run log
    runlog_path = sweep_root / "runlog.csv"
    runlog_path.write_text(
        "run_tag,region,vars,avg_method,log_transform,mean_or_median,exclude_dec_jan,status\n",
        encoding="utf-8",
    )

    print(f"Sweep root: {sweep_root}")
    print(f"Combos: {len(combos)}")
    if limit is not None:
        print(f"BT_SWEEP_LIMIT={limit} (smoke-test mode)")

    for i, (region, vars_to_analyze, avg_method, log_tf, mm, exdj) in enumerate(combos, start=1):
        if limit is not None and i > limit:
            break

        tag = _run_tag(
            region=region,
            vars_to_analyze=vars_to_analyze,
            annual_avg_method=avg_method,
            log_transform=log_tf,
            mean_or_median=mm,
            exclude_dec_jan=exdj,
        )

        run_dir = runs_dir / tag
        figs_dir = run_dir / "figs"
        run_dir.mkdir(parents=True, exist_ok=True)
        figs_dir.mkdir(parents=True, exist_ok=True)

        # shared caches
        c09_cache_csv = cache_dir / f"ecospace_bloom_timing_C09_{scenario}.csv"
        mask_cache_nc = cache_dir / f"suchy_ecospace_mask_{region}.nc"

        env = os.environ.copy()
        env.update(
            {
                # output layout
                "BT_RUN_TAG": tag,
                "BT_OUTDIR": str(run_dir),
                "BT_FIGS_OUTDIR": str(figs_dir),
                "BT_CACHE_OUTDIR": str(cache_dir),

                # suppress interactive windows
                "BT_SHOW_PLOTS": "0",
                "MPLBACKEND": "Agg",

                # sweep knobs
                "BT_MASK_REGNM": region,
                "BT_VARS_TO_ANALYZE_SAT": ",".join(vars_to_analyze),
                "BT_ANNUAL_AVG_METHOD_SAT": avg_method,
                "BT_LOG_TRANSFORM_SAT": "1" if log_tf else "0",
                "BT_MEAN_OR_MEDIAN_SAT": mm,
                "BT_EXCLUDE_DEC_JAN_SAT": "1" if exdj else "0",

                # recompute/caching controls
                "BT_RECOMPUTE_BLOOM_TIMING_SAT": "1",  # sat depends on sweep settings
                "BT_RECOMPUTE_BLOOM_TIMING_C09": "0" if c09_cache_csv.exists() else "1",
                "BT_CREATE_MASKS": "0" if mask_cache_nc.exists() else "1",
            }
        )

        cmd = [
            sys.executable,
            "-c",
            "import ecospace_eval_4a_bloomt_CSoG_sweep as m; m.run_bt_eval()",
        ]

        print(f"[{i}/{len(combos)}] {tag}")
        try:
            res = subprocess.run(
                cmd,
                env=env,
                cwd=str(Path.cwd()),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            (run_dir / "stdout.txt").write_text(res.stdout, encoding="utf-8")
            status = "ok" if res.returncode == 0 else f"fail_{res.returncode}"
        except Exception as e:
            status = f"exception_{type(e).__name__}"
            (run_dir / "stdout.txt").write_text(str(e), encoding="utf-8")

        with runlog_path.open("a", encoding="utf-8") as f:
            f.write(
                f"{tag},{region},\"{';'.join(vars_to_analyze)}\",{avg_method},{int(log_tf)},{mm},{int(exdj)},{status}\n"
            )

        if status != "ok":
            print(f"  -> {status} (see {run_dir / 'stdout.txt'})")

    print(f"Done. Run log: {runlog_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
