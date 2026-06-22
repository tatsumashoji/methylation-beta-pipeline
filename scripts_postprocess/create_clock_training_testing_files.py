#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Create per-clock training/testing CSV files by merging EPICv2/MSA clock.csv outputs
with training/testing template CSV files.

Expected inputs:
  - postprocess/EPICv2/clock.csv
  - postprocess/MSA/clock.csv
  - 251104_DNAmEpiclockAge_MSA_vs_EPICv2_training.csv
  - 251104_DNAmEpiclockAge_MSA_vs_EPICv2_testing.csv

Outputs:
  - <out_dir>/<ClockName>.training.csv
  - <out_dir>/<ClockName>.testing.csv

For training:
  columns = Name, Gender, Chronological_Age, MSA, EPICv2

For testing:
  columns = first 5 columns of testing template + EPIC + MSA
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd


DEFAULT_CLOCKS = [
    "Horvath",
    "Hannum",
    "PhenoAge",
    "GrimAgeV2",
    "DunedinPACE",
    "EpiclockAge",
]


def normalize_sample_id(x: object) -> str:
    """Normalize sample IDs for fallback matching.

    Exact matching is always tried first. This normalization is used only as a fallback.
    It removes common tissue suffixes such as _Blood/_Cheek and trims whitespace.
    """
    s = str(x).strip()
    s = re.sub(r"\.idat(?:\.gz)?$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"_(Red|Grn)$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"_(Blood|Cheek)$", "", s, flags=re.IGNORECASE)
    return s


def load_clock_csv(path: Path) -> pd.DataFrame:
    """Load clock.csv robustly and return DataFrame indexed by sample ID."""
    df = pd.read_csv(path)

    if df.shape[1] < 2:
        raise ValueError(f"{path} has too few columns: {df.shape}")

    first = str(df.columns[0])
    sample_col_candidates = [
        "Sample_ID", "SampleID", "sample_id", "sample", "Sample",
        "Name", "ID", "id", "Unnamed: 0", ""
    ]

    if first in sample_col_candidates or first.startswith("Unnamed"):
        df = df.set_index(df.columns[0])
    else:
        # If first column is not obviously a sample column but contains non-numeric values,
        # treat it as sample ID. Otherwise keep as is and require explicit sample-like index.
        first_vals = df.iloc[:, 0].astype(str)
        if first_vals.str.contains(r"[A-Za-z_]", regex=True).any():
            df = df.set_index(df.columns[0])
        else:
            raise ValueError(
                f"Cannot infer sample ID column in {path}. "
                f"First columns: {df.columns[:5].tolist()}"
            )

    df.index = df.index.astype(str).str.strip()
    df = df.loc[~df.index.isna(), :]

    # Convert clock columns to numeric where possible.
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if df.index.duplicated().any():
        dup = df.index[df.index.duplicated()].unique().tolist()[:10]
        raise ValueError(f"Duplicated sample IDs in {path}: {dup}")

    return df


def build_lookup(clock_df: pd.DataFrame, clock_name: str) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Return exact and normalized lookup maps for one clock."""
    if clock_name not in clock_df.columns:
        raise KeyError(
            f"Clock '{clock_name}' was not found. Available columns: {clock_df.columns.tolist()}"
        )

    exact = clock_df[clock_name].to_dict()

    # Build normalized map only for unique normalized keys.
    tmp = pd.DataFrame({
        "sample": clock_df.index.astype(str),
        "norm": [normalize_sample_id(x) for x in clock_df.index],
        "value": clock_df[clock_name].values,
    })

    counts = tmp["norm"].value_counts()
    unique_norms = set(counts[counts == 1].index)
    norm_map = tmp.loc[tmp["norm"].isin(unique_norms)].set_index("norm")["value"].to_dict()

    return exact, norm_map


def get_clock_value(sample_id: object, exact: Dict[str, float], norm_map: Dict[str, float]) -> float:
    """Get value by exact sample ID, then normalized sample ID, then sample_Blood fallback."""
    sid = str(sample_id).strip()

    if sid in exact:
        return exact[sid]

    # Common fallback: template uses YS_246 while clock index may be YS_246_Blood.
    sid_blood = f"{sid}_Blood"
    if sid_blood in exact:
        return exact[sid_blood]

    sid_cheek = f"{sid}_Cheek"
    if sid_cheek in exact:
        return exact[sid_cheek]

    nsid = normalize_sample_id(sid)
    if nsid in norm_map:
        return norm_map[nsid]

    return np.nan


def make_training_file(
    template: pd.DataFrame,
    epic_clock: pd.DataFrame,
    msa_clock: pd.DataFrame,
    clock_name: str,
    out_path: Path,
) -> pd.DataFrame:
    required = ["Name", "Gender", "Chronological_Age"]
    missing = [c for c in required if c not in template.columns]
    if missing:
        raise KeyError(f"Training template is missing required columns: {missing}")

    epic_exact, epic_norm = build_lookup(epic_clock, clock_name)
    msa_exact, msa_norm = build_lookup(msa_clock, clock_name)

    out = template.loc[:, required].copy()
    out["MSA"] = [get_clock_value(x, msa_exact, msa_norm) for x in out["Name"]]
    out["EPICv2"] = [get_clock_value(x, epic_exact, epic_norm) for x in out["Name"]]

    out.to_csv(out_path, index=False)
    return out


def make_testing_file(
    template: pd.DataFrame,
    epic_clock: pd.DataFrame,
    msa_clock: pd.DataFrame,
    clock_name: str,
    out_path: Path,
) -> pd.DataFrame:
    if "Name" not in template.columns:
        raise KeyError("Testing template is missing required column: Name")

    if template.shape[1] < 5:
        raise ValueError("Testing template must have at least 5 columns.")

    first5 = template.columns[:5].tolist()

    epic_exact, epic_norm = build_lookup(epic_clock, clock_name)
    msa_exact, msa_norm = build_lookup(msa_clock, clock_name)

    out = template.loc[:, first5].copy()
    out["EPIC"] = [get_clock_value(x, epic_exact, epic_norm) for x in out["Name"]]
    out["MSA"] = [get_clock_value(x, msa_exact, msa_norm) for x in out["Name"]]

    out.to_csv(out_path, index=False)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--epicv2-clock", default="postprocess/EPICv2/clock.csv")
    ap.add_argument("--msa-clock", default="postprocess/MSA/clock.csv")
    ap.add_argument("--training-template", default="251104_DNAmEpiclockAge_MSA_vs_EPICv2_training.csv")
    ap.add_argument("--testing-template", default="251104_DNAmEpiclockAge_MSA_vs_EPICv2_testing.csv")
    ap.add_argument("--out-dir", default="clock_compare_outputs")
    ap.add_argument("--clocks", nargs="+", default=DEFAULT_CLOCKS)
    ap.add_argument("--strict", action="store_true", help="fail if any EPIC/MSA clock values are missing")
    args = ap.parse_args()

    epic_path = Path(args.epicv2_clock)
    msa_path = Path(args.msa_clock)
    train_path = Path(args.training_template)
    test_path = Path(args.testing_template)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[Load] EPICv2 clock: {epic_path}")
    epic_clock = load_clock_csv(epic_path)
    print(f"  shape: {epic_clock.shape}")

    print(f"[Load] MSA clock: {msa_path}")
    msa_clock = load_clock_csv(msa_path)
    print(f"  shape: {msa_clock.shape}")

    print(f"[Load] training template: {train_path}")
    train = pd.read_csv(train_path)
    print(f"  shape: {train.shape}, columns: {train.columns.tolist()}")

    print(f"[Load] testing template: {test_path}")
    test = pd.read_csv(test_path)
    print(f"  shape: {test.shape}, columns: {test.columns.tolist()}")

    summary_rows = []

    for clock in args.clocks:
        print(f"\n[Clock] {clock}")

        train_out_path = out_dir / f"{clock}.training.csv"
        test_out_path = out_dir / f"{clock}.testing.csv"

        train_out = make_training_file(train, epic_clock, msa_clock, clock, train_out_path)
        test_out = make_testing_file(test, epic_clock, msa_clock, clock, test_out_path)

        train_missing_msa = int(train_out["MSA"].isna().sum())
        train_missing_epic = int(train_out["EPICv2"].isna().sum())
        test_missing_msa = int(test_out["MSA"].isna().sum())
        test_missing_epic = int(test_out["EPIC"].isna().sum())

        summary_rows.append({
            "clock": clock,
            "training_file": str(train_out_path),
            "testing_file": str(test_out_path),
            "training_n": train_out.shape[0],
            "testing_n": test_out.shape[0],
            "training_missing_MSA": train_missing_msa,
            "training_missing_EPICv2": train_missing_epic,
            "testing_missing_MSA": test_missing_msa,
            "testing_missing_EPIC": test_missing_epic,
        })

        print(f"  wrote: {train_out_path}")
        print(f"  wrote: {test_out_path}")
        print(
            f"  missing: training MSA={train_missing_msa}, training EPICv2={train_missing_epic}, "
            f"testing MSA={test_missing_msa}, testing EPIC={test_missing_epic}"
        )

        if args.strict and any([train_missing_msa, train_missing_epic, test_missing_msa, test_missing_epic]):
            raise RuntimeError(f"Missing values were found for {clock}. See summary/missing report.")

    summary = pd.DataFrame(summary_rows)
    summary_path = out_dir / "summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"\n[Done] summary: {summary_path}")

    missing_total = int(
        summary[
            ["training_missing_MSA", "training_missing_EPICv2", "testing_missing_MSA", "testing_missing_EPIC"]
        ].sum().sum()
    )
    if missing_total > 0:
        print(f"[Warning] Total missing clock values: {missing_total}")
        print("          Use summary.csv to check which clock/source has missing samples.")
        print("          Add --strict if you want the script to fail on missing values.")


if __name__ == "__main__":
    main()
