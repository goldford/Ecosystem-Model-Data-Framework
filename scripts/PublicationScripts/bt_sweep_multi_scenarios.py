"""
bt_sweep_multi_scenarios.py

Run the bloom-timing parameter sweep (bt_sweep_4a.py) for multiple Ecospace
scenarios *sequentially*, so you can "set it and forget it" (as long as your PC
stays awake).

Prereqs:
  - Make sure your active config module supports scenario override via env vars:
      ECOSPACE_SC, ECOSPACE_SC_FULL
    (Use ecospace_eval_config_sweep_envscenario.py as a drop-in replacement for
     ecospace_eval_config_sweep.py.)

Typical usage:
  python bt_sweep_multi_scenarios.py --scenarios SC205 SC206 SC207

Optional:
  python bt_sweep_multi_scenarios.py --scenarios SC205 SC206 SC207 --leaderboard
  python bt_sweep_multi_scenarios.py --scenarios SC205 SC206 SC207 --base-root "C:\path\to\outputs"
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

SWEEP_ROOT_RE = re.compile(r"^Sweep root:\s*(.+)\s*$")


def _stream_subprocess(cmd: List[str], env: dict, cwd: Path) -> Tuple[int, str]:
    """
    Run a subprocess, stream stdout live, and also capture it (so we can parse Sweep root).
    Returns: (returncode, full_stdout)
    """
    proc = subprocess.Popen(
        cmd,
        env=env,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert proc.stdout is not None
    out_lines: List[str] = []
    for line in proc.stdout:
        print(line, end="")  # live stream (nice in PyCharm)
        out_lines.append(line)
    proc.wait()
    return proc.returncode, "".join(out_lines)


def _parse_sweep_root(stdout: str) -> Optional[str]:
    for line in stdout.splitlines():
        m = SWEEP_ROOT_RE.match(line.strip())
        if m:
            return m.group(1)
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--scenarios",
        nargs="+",
        default=["SC205", "SC206", "SC207"],
        help="Scenario codes to run sequentially, e.g. SC205 SC206 SC207",
    )
    ap.add_argument(
        "--sweep-script",
        default="bt_sweep_4a.py",
        help="Path to the sweep script (default: bt_sweep_4a.py).",
    )
    ap.add_argument(
        "--base-root",
        default=None,
        help="Optional base output folder. If set, each scenario uses <base-root>/bt_sweep_4a_<SCENARIO>/",
    )
    ap.add_argument(
        "--leaderboard",
        action="store_true",
        help="Also run bt_results_leaderboard.py after each scenario sweep completes.",
    )
    ap.add_argument(
        "--leaderboard-script",
        default="bt_results_leaderboard.py",
        help="Leaderboard builder script path (default: bt_results_leaderboard.py).",
    )
    ap.add_argument(
        "--continue-on-error",
        action="store_true",
        help="If a scenario sweep fails, continue to the next scenario.",
    )
    args = ap.parse_args()

    cwd = Path.cwd()

    sweep_script = Path(args.sweep_script)
    if not sweep_script.is_file():
        print(f"ERROR: sweep script not found: {sweep_script}")
        return 2

    if args.leaderboard:
        lb_script = Path(args.leaderboard_script)
        if not lb_script.is_file():
            print(f"ERROR: leaderboard requested but script not found: {lb_script}")
            return 2

    summary_rows: List[str] = []

    for s in args.scenarios:
        scenario = s.strip()
        if not scenario:
            continue

        env = os.environ.copy()

        # --- key: scenario override so you don't edit config between runs ---
        env["ECOSPACE_SC"] = scenario
        env["ECOSPACE_SC_FULL"] = scenario

        # Ensure distinct output root per scenario (bt_sweep_4a.py will add its own timestamp folder)
        if args.base_root:
            base_root = Path(args.base_root).expanduser()
            env["BT_SWEEP_ROOT"] = str((base_root / f"bt_sweep_4a_{scenario}").resolve())
        else:
            # Explicit default to keep things predictable from any working dir
            env["BT_SWEEP_ROOT"] = str((cwd / f"..//..//data//evaluation//bt_sweep_4a_{scenario}").resolve())

        # Nice-to-have: immediate flushing in PyCharm console
        env.setdefault("PYTHONUNBUFFERED", "1")

        cmd = [sys.executable, str(sweep_script)]
        print("\n" + "=" * 88)
        print(f"RUNNING SWEEP FOR {scenario}")
        print("=" * 88)

        rc, stdout = _stream_subprocess(cmd, env=env, cwd=cwd)
        sweep_root = _parse_sweep_root(stdout) or ""

        status = "ok" if rc == 0 else f"fail_{rc}"
        summary_rows.append(f"{scenario},{status},{sweep_root}")

        if rc != 0:
            print(f"\nSweep for {scenario} FAILED with return code {rc}.")
            if not args.continue_on_error:
                break
            continue

        if args.leaderboard and sweep_root:
            print("\n" + "-" * 88)
            print(f"BUILDING LEADERBOARD for {scenario}")
            print("-" * 88)

            lb_cmd = [sys.executable, str(args.leaderboard_script), "--root", sweep_root]
            lb_rc, _ = _stream_subprocess(lb_cmd, env=os.environ.copy(), cwd=cwd)
            if lb_rc != 0:
                print(f"Leaderboard step failed for {scenario} (rc={lb_rc}). Continuing.")

    print("\n" + "=" * 88)
    print("MULTI-SCENARIO SWEEP SUMMARY")
    print("=" * 88)
    print("scenario,status,sweep_root")
    for row in summary_rows:
        print(row)

    # Also write a small summary CSV for convenience
    try:
        out_csv = cwd / "bt_multi_sweep_summary.csv"
        out_csv.write_text("scenario,status,sweep_root\n" + "\n".join(summary_rows) + "\n", encoding="utf-8")
        print(f"\nWrote summary CSV: {out_csv}")
    except Exception as e:
        print(f"\nCould not write summary CSV: {e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
