#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_CLOCKS = [
    "Horvath",
    "Hannum",
    "PhenoAge",
    "GrimAgeV1",
    "GrimAgeV2",
    "DunedinPACE",
    "DunedinPACE_20k",
    "DunedinPACE_original_lowmem",
]


def read_csv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}")
    return pd.read_csv(path)


def summarize_series(
    values: pd.Series,
    dataset: str,
    clock: str,
    group: str = "overall",
) -> dict:
    numeric = pd.to_numeric(values, errors="coerce")

    n_total = len(numeric)
    n_nonmissing = int(numeric.notna().sum())
    n_missing = int(numeric.isna().sum())

    if n_nonmissing == 0:
        return {
            "dataset": dataset,
            "group": group,
            "clock": clock,
            "n_total": n_total,
            "n_nonmissing": 0,
            "n_missing": n_missing,
            "mean": None,
            "sd": None,
            "median": None,
            "min": None,
            "max": None,
            "q1": None,
            "q3": None,
        }

    return {
        "dataset": dataset,
        "group": group,
        "clock": clock,
        "n_total": n_total,
        "n_nonmissing": n_nonmissing,
        "n_missing": n_missing,
        "mean": float(numeric.mean()),
        "sd": float(numeric.std(ddof=1)) if n_nonmissing >= 2 else 0.0,
        "median": float(numeric.median()),
        "min": float(numeric.min()),
        "max": float(numeric.max()),
        "q1": float(numeric.quantile(0.25)),
        "q3": float(numeric.quantile(0.75)),
    }


def summarize_clock_csv(
    path: Path,
    dataset: str,
    clocks: list[str],
) -> list[dict]:
    df = read_csv(path)

    rows = []
    for clock in clocks:
        if clock not in df.columns:
            continue
        rows.append(
            summarize_series(
                df[clock],
                dataset=dataset,
                clock=clock,
                group="overall",
            )
        )

    return rows


def summarize_split_files(
    clock_compare_dir: Path,
    clocks: list[str],
) -> list[dict]:
    rows = []

    if not clock_compare_dir.is_dir():
        return rows

    for clock in clocks:
        # Do not look for files for alternative DunedinPACE columns.
        # clock_compare_outputs usually contains DunedinPACE.training.csv only.
        file_clock = clock
        if clock in {"DunedinPACE_20k", "DunedinPACE_original_lowmem"}:
            continue

        for split in ["training", "testing"]:
            path = clock_compare_dir / f"{file_clock}.{split}.csv"
            if not path.is_file():
                continue

            df = pd.read_csv(path)

            for dataset in ["EPICv2", "MSA"]:
                if dataset not in df.columns:
                    continue

                rows.append(
                    summarize_series(
                        df[dataset],
                        dataset=dataset,
                        clock=clock,
                        group=split,
                    )
                )

    return rows


def make_overall_wide(summary_long: pd.DataFrame) -> pd.DataFrame:
    if summary_long.empty:
        return pd.DataFrame()

    rows = []

    for clock in summary_long["clock"].drop_duplicates():
        sub = summary_long[summary_long["clock"] == clock]

        row = {"clock": clock}

        for dataset in ["EPICv2", "MSA"]:
            ds = sub[sub["dataset"] == dataset]
            if ds.empty:
                continue

            r = ds.iloc[0]
            row[f"{dataset}_n"] = r["n_nonmissing"]
            row[f"{dataset}_mean"] = r["mean"]
            row[f"{dataset}_sd"] = r["sd"]
            row[f"{dataset}_median"] = r["median"]

        if "EPICv2_mean" in row and "MSA_mean" in row:
            row["MSA_minus_EPICv2_mean"] = (
                row["MSA_mean"] - row["EPICv2_mean"]
            )

        if "EPICv2_sd" in row and "MSA_sd" in row:
            row["MSA_minus_EPICv2_sd"] = (
                row["MSA_sd"] - row["EPICv2_sd"]
            )

        rows.append(row)

    return pd.DataFrame(rows)


def make_split_wide(summary_long: pd.DataFrame) -> pd.DataFrame:
    if summary_long.empty:
        return pd.DataFrame()

    rows = []

    for group in summary_long["group"].drop_duplicates():
        for clock in summary_long["clock"].drop_duplicates():
            sub = summary_long[
                (summary_long["group"] == group)
                & (summary_long["clock"] == clock)
            ]

            if sub.empty:
                continue

            row = {
                "group": group,
                "clock": clock,
            }

            for dataset in ["EPICv2", "MSA"]:
                ds = sub[sub["dataset"] == dataset]
                if ds.empty:
                    continue

                r = ds.iloc[0]
                row[f"{dataset}_n"] = r["n_nonmissing"]
                row[f"{dataset}_mean"] = r["mean"]
                row[f"{dataset}_sd"] = r["sd"]
                row[f"{dataset}_median"] = r["median"]

            if "EPICv2_mean" in row and "MSA_mean" in row:
                row["MSA_minus_EPICv2_mean"] = (
                    row["MSA_mean"] - row["EPICv2_mean"]
                )

            rows.append(row)

    return pd.DataFrame(rows)


def round_numeric(df: pd.DataFrame, digits: int = 6) -> pd.DataFrame:
    if df.empty:
        return df

    out = df.copy()
    numeric_cols = out.select_dtypes(include=["number"]).columns
    out[numeric_cols] = out[numeric_cols].round(digits)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize mean and SD of epigenetic clock values for EPICv2 and MSA."
    )
    parser.add_argument(
        "--epicv2-clock",
        type=Path,
        default=Path("/work/postprocess/EPICv2_biolearn_default_grimagev1/clock.csv"),
    )
    parser.add_argument(
        "--msa-clock",
        type=Path,
        default=Path("/work/postprocess/MSA_biolearn_default_grimagev1/clock.csv"),
    )
    parser.add_argument(
        "--clock-compare-dir",
        type=Path,
        default=Path("/work/clock_compare_outputs_grimagev1"),
        help="Directory containing Clock.training.csv and Clock.testing.csv files.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("/work/postprocess/clock_summary_stats"),
    )
    parser.add_argument(
        "--clocks",
        nargs="+",
        default=DEFAULT_CLOCKS,
        help="Clock columns to summarize.",
    )
    parser.add_argument(
        "--digits",
        type=int,
        default=6,
    )

    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    overall_rows = []
    overall_rows.extend(
        summarize_clock_csv(
            args.epicv2_clock,
            dataset="EPICv2",
            clocks=args.clocks,
        )
    )
    overall_rows.extend(
        summarize_clock_csv(
            args.msa_clock,
            dataset="MSA",
            clocks=args.clocks,
        )
    )

    overall_long = pd.DataFrame(overall_rows)
    overall_wide = make_overall_wide(overall_long)

    split_rows = summarize_split_files(
        args.clock_compare_dir,
        clocks=args.clocks,
    )
    split_long = pd.DataFrame(split_rows)
    split_wide = make_split_wide(split_long)

    overall_long = round_numeric(overall_long, args.digits)
    overall_wide = round_numeric(overall_wide, args.digits)
    split_long = round_numeric(split_long, args.digits)
    split_wide = round_numeric(split_wide, args.digits)

    overall_long_path = args.out_dir / "overall_clock_mean_sd_long.csv"
    overall_wide_path = args.out_dir / "overall_clock_mean_sd_wide.csv"
    split_long_path = args.out_dir / "split_clock_mean_sd_long.csv"
    split_wide_path = args.out_dir / "split_clock_mean_sd_wide.csv"
    xlsx_path = args.out_dir / "clock_mean_sd_summary.xlsx"

    overall_long.to_csv(overall_long_path, index=False)
    overall_wide.to_csv(overall_wide_path, index=False)

    if not split_long.empty:
        split_long.to_csv(split_long_path, index=False)
        split_wide.to_csv(split_wide_path, index=False)

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        overall_wide.to_excel(writer, sheet_name="Overall_Wide", index=False)
        overall_long.to_excel(writer, sheet_name="Overall_Long", index=False)

        if not split_wide.empty:
            split_wide.to_excel(writer, sheet_name="Split_Wide", index=False)
            split_long.to_excel(writer, sheet_name="Split_Long", index=False)

    print("[Done] Summary files written to:")
    print(f"  {overall_long_path}")
    print(f"  {overall_wide_path}")
    if not split_long.empty:
        print(f"  {split_long_path}")
        print(f"  {split_wide_path}")
    print(f"  {xlsx_path}")

    print()
    print("[Overall summary]")
    if not overall_wide.empty:
        print(overall_wide.to_string(index=False))


if __name__ == "__main__":
    main()
