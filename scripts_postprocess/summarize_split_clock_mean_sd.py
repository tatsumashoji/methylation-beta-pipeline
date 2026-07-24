#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_CLOCKS = [
    "Horvath",
    "Hannum",
    "PhenoAge",
    "GrimAgeV1",
    "GrimAgeV2",
    "DunedinPACE",
]


def choose_epic_col(df: pd.DataFrame) -> str | None:
    if "EPICv2" in df.columns:
        return "EPICv2"
    if "EPIC" in df.columns:
        return "EPIC"
    return None


def summarize_series(values: pd.Series) -> dict:
    s = pd.to_numeric(values, errors="coerce")

    return {
        "n_total": len(s),
        "n_nonmissing": int(s.notna().sum()),
        "n_missing": int(s.isna().sum()),
        "mean": float(s.mean()) if s.notna().sum() else None,
        "sd": float(s.std(ddof=1)) if s.notna().sum() >= 2 else None,
        "median": float(s.median()) if s.notna().sum() else None,
        "min": float(s.min()) if s.notna().sum() else None,
        "max": float(s.max()) if s.notna().sum() else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--clock-compare-dir",
        type=Path,
        default=Path("/work/clock_compare_outputs_grimagev1"),
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
    )
    parser.add_argument("--digits", type=int, default=6)

    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    long_rows = []
    wide_rows = []
    input_check_rows = []

    for clock in args.clocks:
        for split in ["training", "testing"]:
            path = args.clock_compare_dir / f"{clock}.{split}.csv"

            if not path.is_file():
                input_check_rows.append({
                    "clock": clock,
                    "split": split,
                    "file": str(path),
                    "status": "missing",
                    "n_rows": None,
                    "paired_nonmissing": None,
                })
                continue

            df = pd.read_csv(path)
            epic_col = choose_epic_col(df)

            if epic_col is None or "MSA" not in df.columns:
                input_check_rows.append({
                    "clock": clock,
                    "split": split,
                    "file": str(path),
                    "status": "missing_EPIC_or_MSA_column",
                    "columns": ";".join(df.columns),
                    "n_rows": len(df),
                    "paired_nonmissing": None,
                })
                continue

            epic = pd.to_numeric(df[epic_col], errors="coerce")
            msa = pd.to_numeric(df["MSA"], errors="coerce")
            paired_nonmissing = int((epic.notna() & msa.notna()).sum())

            input_check_rows.append({
                "clock": clock,
                "split": split,
                "file": str(path),
                "status": "ok",
                "n_rows": len(df),
                "epic_column": epic_col,
                "EPICv2_nonmissing": int(epic.notna().sum()),
                "MSA_nonmissing": int(msa.notna().sum()),
                "paired_nonmissing": paired_nonmissing,
            })

            for dataset, col in [("EPICv2", epic_col), ("MSA", "MSA")]:
                stats = summarize_series(df[col])
                row = {
                    "clock": clock,
                    "split": split,
                    "dataset": dataset,
                }
                row.update(stats)
                long_rows.append(row)

            epic_stats = summarize_series(df[epic_col])
            msa_stats = summarize_series(df["MSA"])

            wide_row = {
                "clock": clock,
                "split": split,
                "EPICv2_n": epic_stats["n_nonmissing"],
                "EPICv2_mean": epic_stats["mean"],
                "EPICv2_sd": epic_stats["sd"],
                "EPICv2_median": epic_stats["median"],
                "MSA_n": msa_stats["n_nonmissing"],
                "MSA_mean": msa_stats["mean"],
                "MSA_sd": msa_stats["sd"],
                "MSA_median": msa_stats["median"],
                "paired_nonmissing": paired_nonmissing,
            }

            if epic_stats["mean"] is not None and msa_stats["mean"] is not None:
                wide_row["MSA_minus_EPICv2_mean"] = (
                    msa_stats["mean"] - epic_stats["mean"]
                )

            if epic_stats["sd"] is not None and msa_stats["sd"] is not None:
                wide_row["MSA_minus_EPICv2_sd"] = (
                    msa_stats["sd"] - epic_stats["sd"]
                )

            wide_rows.append(wide_row)

    long_df = pd.DataFrame(long_rows)
    wide_df = pd.DataFrame(wide_rows)
    check_df = pd.DataFrame(input_check_rows)

    for df in [long_df, wide_df]:
        if not df.empty:
            numeric_cols = df.select_dtypes(include=["number"]).columns
            df[numeric_cols] = df[numeric_cols].round(args.digits)

    long_path = args.out_dir / "split_clock_mean_sd_long.csv"
    wide_path = args.out_dir / "split_clock_mean_sd_wide.csv"
    check_path = args.out_dir / "split_clock_input_check.csv"
    xlsx_path = args.out_dir / "split_clock_mean_sd_summary.xlsx"

    long_df.to_csv(long_path, index=False)
    wide_df.to_csv(wide_path, index=False)
    check_df.to_csv(check_path, index=False)

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        wide_df.to_excel(writer, sheet_name="Mean_SD_Wide", index=False)
        long_df.to_excel(writer, sheet_name="Mean_SD_Long", index=False)
        check_df.to_excel(writer, sheet_name="Input_Check", index=False)

    print("[Done] Wrote:")
    print(f"  {wide_path}")
    print(f"  {long_path}")
    print(f"  {check_path}")
    print(f"  {xlsx_path}")
    print()
    print(wide_df.to_string(index=False))


if __name__ == "__main__":
    main()
