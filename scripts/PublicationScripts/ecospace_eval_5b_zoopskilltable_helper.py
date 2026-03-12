from __future__ import annotations

from pathlib import Path
import pandas as pd
import re


SEASON_ORDER = ["Winter", "Spring", "Summer", "Fall", "All"]


def _parse_filename(path: Path) -> tuple[str, str, str]:
    """
    Expect names like:
        ecospace_zoop_performance_SC215_ZC_Spring.csv
        ecospace_zoop_performance_SC215_ZS_Winter.csv

    Returns
    -------
    scenario, pass_name, season
    """
    m = re.match(
        r"ecospace_zoop_performance_(?P<scenario>.+?)_(?P<pass>ZC|ZS)_(?P<season>All|Winter|Spring|Summer|Fall)\.csv$",
        path.name,
    )
    if not m:
        raise ValueError(f"Filename does not match expected pattern: {path.name}")
    return m.group("scenario"), m.group("pass"), m.group("season")


def _season_sort(df: pd.DataFrame, season_col: str = "season") -> pd.DataFrame:
    order_map = {s: i for i, s in enumerate(SEASON_ORDER)}
    return (
        df.assign(_season_order=df[season_col].map(order_map).fillna(999))
          .sort_values(["_season_order"] + [c for c in df.columns if c != season_col and c != "_season_order"])
          .drop(columns="_season_order")
          .reset_index(drop=True)
    )


def build_manuscript_skill_table(
    scenario: str,
    in_dir: str | Path,
    *,
    out_csv: str | Path | None = None,
    include_all: bool = False,
    total_name_in_source: str = "Total",
) -> pd.DataFrame:
    """
    Build manuscript table:
        season, N, r, RMSE, MAE, WSS

    Uses the TOTAL row from each seasonal ZC file only.
    """
    in_dir = Path(in_dir)
    files = sorted(in_dir.glob(f"ecospace_zoop_performance_{scenario}_ZC_*.csv"))

    rows = []
    for f in files:
        _, _, season = _parse_filename(f)
        if (not include_all) and season == "All":
            continue

        df = pd.read_csv(f)
        sub = df[df["Group"].astype(str) == total_name_in_source].copy()
        if sub.empty:
            continue

        sub.insert(0, "season", season)
        sub = sub[["season", "N", "r", "RMSE", "MAE", "WSS"]]
        rows.append(sub)

    if not rows:
        raise ValueError(f"No matching ZC seasonal files with Group == '{total_name_in_source}' found for scenario {scenario}.")

    out = pd.concat(rows, ignore_index=True)
    out = _season_sort(out, "season")

    if out_csv is None:
        out_csv = in_dir / f"ecospace_forpub_zoop_skill_{scenario}_ZC.csv"

    out.to_csv(out_csv, index=False)
    print(f"Saved manuscript table: {out_csv}")
    return out


def build_supplement_skill_table(
    scenario: str,
    in_dir: str | Path,
    *,
    out_csv: str | Path | None = None,
    include_all: bool = False,
    passes: tuple[str, ...] = ("ZC", "ZS"),
    total_name_in_source: str = "Total",
    total_name_fmt: str = "{pass_name}-TOTAL",
) -> pd.DataFrame:
    """
    Build supplement table:
        season, group, N, r, RMSE, MAE, WSS

    Stacks ZC and ZS files and keeps individual groups.
    Renames source 'Total' rows to e.g. ZC-TOTAL / ZS-TOTAL.
    """
    in_dir = Path(in_dir)

    rows = []
    for pass_name in passes:
        files = sorted(in_dir.glob(f"ecospace_zoop_performance_{scenario}_{pass_name}_*.csv"))
        for f in files:
            _, _, season = _parse_filename(f)
            if (not include_all) and season == "All":
                continue

            df = pd.read_csv(f).copy()
            df.insert(0, "season", season)
            df["group"] = df["Group"].astype(str)
            df.loc[df["group"] == total_name_in_source, "group"] = total_name_fmt.format(pass_name=pass_name)
            df = df[["season", "group", "N", "r", "RMSE", "MAE", "WSS"]]
            rows.append(df)

    if not rows:
        raise ValueError(f"No matching ZC/ZS seasonal files found for scenario {scenario}.")

    out = pd.concat(rows, ignore_index=True)

    # Keep seasons in the desired order; within season, put totals first, then alphabetical
    season_order = {s: i for i, s in enumerate(SEASON_ORDER)}
    out["_season_order"] = out["season"].map(season_order).fillna(999)
    out["_group_order"] = out["group"].str.contains("TOTAL", na=False).map({True: 0, False: 1})
    out = (
        out.sort_values(["_season_order", "_group_order", "group"])
           .drop(columns=["_season_order", "_group_order"])
           .reset_index(drop=True)
    )

    if out_csv is None:
        out_csv = in_dir / f"ecospace_forpub_zoop_skill_{scenario}_ALLGROUPS.csv"

    out.to_csv(out_csv, index=False)
    print(f"Saved supplement table: {out_csv}")
    return out


def export_forpub_zoop_skill_tables(
    scenario: str,
    in_dir: str | Path,
    *,
    include_all: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Convenience wrapper:
      1) manuscript table from ZC total rows
      2) supplement table from ZC + ZS seasonal files
    """
    ms = build_manuscript_skill_table(
        scenario=scenario,
        in_dir=in_dir,
        include_all=include_all,
    )
    supp = build_supplement_skill_table(
        scenario=scenario,
        in_dir=in_dir,
        include_all=include_all,
    )
    return ms, supp


if __name__ == "__main__":
    # Example:
    #   python zoop_skill_table_helper.py
    scenario = "SC215"
    in_dir = Path("..//..//data//evaluation//")
    export_forpub_zoop_skill_tables(scenario=scenario, in_dir=in_dir, include_all=False)
