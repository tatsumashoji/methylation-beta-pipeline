#!/usr/bin/env python3
"""
Merge batch-wise SeSAMe collapsed_beta_matrix.txt files into one sample-row beta matrix,
and rename sample IDs using an EPICv2_blood.txt / MSA_blood.txt sample sheet.

Input batch matrix default orientation:
    rows    = CpGs / probe IDs
    columns = IDAT prefixes / samples

Output:
    rows    = samples
    columns = CpGs / probe IDs

Example:
    python scripts_postprocess/prepare_collapsed_beta.py \
      --batch-root results/EPICv2_batches \
      --sample-sheet EPICv2_blood.txt \
      --out-dir postprocess/EPICv2
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


def log(msg: str) -> None:
    print(msg, flush=True)


def normalize_id(value: object) -> str:
    """Normalize a sample/path/IDAT filename to an Illumina IDAT prefix."""
    if pd.isna(value):
        return ""
    s = str(value).strip().strip('"').strip("'")
    if not s:
        return ""
    s = os.path.basename(s)
    if s.startswith("._"):
        s = s[2:]

    # Remove common compression and IDAT suffixes.
    for suffix in [
        "_Red.idat.gz", "_Grn.idat.gz",
        "_Red.idat", "_Grn.idat",
        ".idat.gz", ".idat",
        ".gz",
    ]:
        if s.endswith(suffix):
            s = s[: -len(suffix)]
            break
    return s


def read_sample_sheet(path: Path) -> Tuple[Dict[str, str], pd.DataFrame]:
    """
    Build mapping from IDAT prefix / sample name to first column sample name.

    Expected columns are like:
        Sample_Name1, Sample_Name2, Grn, Red

    The output sample name is always the first column, as requested.
    """
    if not path.exists():
        raise FileNotFoundError(f"sample sheet does not exist: {path}")

    # sep=None lets pandas infer tab/comma. EPICv2_blood.txt is tab-separated.
    ss = pd.read_csv(path, sep=None, engine="python", dtype=str)
    if ss.shape[1] < 1:
        raise ValueError(f"sample sheet has no columns: {path}")

    first_col = ss.columns[0]
    mapping: Dict[str, str] = {}

    def put(key: object, value: object):
        k = normalize_id(key)
        v = str(value).strip() if not pd.isna(value) else ""
        if not k or not v:
            return
        if k in mapping and mapping[k] != v:
            raise ValueError(
                f"Duplicate mapping with different values for key={k}: "
                f"{mapping[k]} vs {v}"
            )
        mapping[k] = v

    for _, row in ss.iterrows():
        target = row[first_col]

        # map sample names as well as IDAT paths
        for col in ss.columns:
            put(row[col], target)

        # explicitly map basename prefixes from Red/Grn columns when present
        for col in ["Red", "Grn", "red", "grn"]:
            if col in ss.columns:
                put(row[col], target)

    return mapping, ss


def rename_sample(sample_id: object, mapping: Dict[str, str]) -> Tuple[str, str, bool]:
    original = str(sample_id)
    key = normalize_id(original)
    if key in mapping:
        return mapping[key], key, True

    # Some SeSAMe names may already be sample names or contain path-like prefixes.
    if original in mapping:
        return mapping[original], original, True

    return original, key, False


def find_batch_files(batch_root: Path, filename: str, start_batch: int = 1) -> List[Path]:
    if not batch_root.exists():
        raise FileNotFoundError(f"batch root does not exist: {batch_root}")

    files = []
    for p in sorted(batch_root.glob("batch_*")):
        if not p.is_dir():
            continue
        m = re.search(r"batch_(\d+)$", p.name)
        if m and int(m.group(1)) < start_batch:
            continue
        f = p / filename
        if f.exists():
            files.append(f)
        else:
            log(f"[WARN] missing {filename}: {f}")

    if not files:
        raise FileNotFoundError(
            f"No {filename} files found under {batch_root}/batch_*"
        )
    return files


def cpg_like_ratio(values: Iterable[object], max_n: int = 1000) -> float:
    vals = list(values)[:max_n]
    if not vals:
        return 0.0
    n = 0
    for x in vals:
        s = str(x)
        if s.startswith("cg") or s.startswith("ch") or s.startswith("rs"):
            n += 1
    return n / len(vals)


def orient_matrix(df: pd.DataFrame, orientation: str) -> pd.DataFrame:
    """
    Return sample-row matrix.

    cpg_rows:   input rows CpGs, columns samples -> transpose
    sample_rows: input already rows samples -> as-is
    auto: decide from CpG-like row/column labels
    """
    if orientation == "cpg_rows":
        return df.T
    if orientation == "sample_rows":
        return df
    if orientation != "auto":
        raise ValueError(f"unknown orientation: {orientation}")

    row_ratio = cpg_like_ratio(df.index)
    col_ratio = cpg_like_ratio(df.columns)
    if row_ratio >= col_ratio:
        return df.T
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-root", required=True, help="e.g. results/EPICv2_batches")
    ap.add_argument("--sample-sheet", required=True, help="EPICv2_blood.txt or MSA_blood.txt")
    ap.add_argument("--out-dir", required=True, help="output directory")
    ap.add_argument("--batch-file-name", default="collapsed_beta_matrix.txt")
    ap.add_argument("--output-name", default="collapsed_beta_matrix.txt")
    ap.add_argument("--orientation", choices=["cpg_rows", "sample_rows", "auto"], default="cpg_rows")
    ap.add_argument("--start-batch", type=int, default=1)
    ap.add_argument("--dtype", default="float32")
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--float-format", default="%.6g")
    ap.add_argument("--strict-renaming", action="store_true",
                    help="fail if any sample could not be mapped to sample sheet")
    ap.add_argument("--drop-unmapped", action="store_true",
                    help="drop samples that could not be mapped to the sample sheet")
    ap.add_argument("--allow-duplicate-samples", action="store_true")
    args = ap.parse_args()

    batch_root = Path(args.batch_root)
    sample_sheet = Path(args.sample_sheet)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / args.output_name
    report_path = out_dir / "sample_renaming_report.csv"
    summary_path = out_dir / "merge_summary.txt"

    mapping, ss = read_sample_sheet(sample_sheet)
    files = find_batch_files(batch_root, args.batch_file_name, args.start_batch)

    log(f"[Input] batch files: {len(files)}")
    log(f"[Input] sample sheet rows: {ss.shape[0]}")
    log(f"[Input] mapping keys: {len(mapping)}")
    log(f"[Output] beta matrix: {out_path}")

    if out_path.exists():
        out_path.unlink()
    if report_path.exists():
        report_path.unlink()

    reference_cols = None
    total_samples = 0
    n_unmapped = 0
    seen_samples = set()
    report_rows = []

    for i, f in enumerate(files, 1):
        log(f"[{i}/{len(files)}] Reading {f}")
        # Do not pass dtype=args.dtype here.
        # The first column contains CpG/probe IDs such as cg00000029,
        # and some files may have non-numeric index/header fields.
        df = pd.read_csv(
            f,
            sep=args.sep,
            index_col=0,
            na_values=["", "NA", "NaN", "nan"],
            low_memory=False,
        )

        # Convert beta-value cells only to numeric after the index has been set.
        # Non-numeric cells become NaN.
        df = df.apply(pd.to_numeric, errors="coerce").astype(args.dtype)

        sdf = orient_matrix(df, args.orientation)
        del df

        new_index = []
        mapped_flags = []

        for old in sdf.index:
            new, key, mapped = rename_sample(old, mapping)
            if not mapped:
                n_unmapped += 1

            report_rows.append({
                "old_sample_id": str(old),
                "normalized_key": key,
                "new_sample_id": new,
                "mapped": mapped,
                "source_file": str(f),
            })

            new_index.append(new)
            mapped_flags.append(mapped)

        sdf.index = new_index
        sdf.index.name = "Sample_ID"

        if args.drop_unmapped:
            keep_mask = np.array(mapped_flags, dtype=bool)
            n_drop_batch = int((~keep_mask).sum())
            if n_drop_batch > 0:
                log(f"  dropping {n_drop_batch} unmapped samples in {f}")
            sdf = sdf.iloc[keep_mask, :]
            if sdf.shape[0] == 0:
                log(f"  no mapped samples remained in {f}; skipping this batch")
                continue

        if not args.allow_duplicate_samples:
            dup = pd.Index(sdf.index)[pd.Index(sdf.index).duplicated()].unique().tolist()
            dup_global = [x for x in sdf.index if x in seen_samples]
            if dup or dup_global:
                raise ValueError(
                    "Duplicate sample names after renaming. "
                    f"within_batch={dup[:10]}, already_seen={dup_global[:10]}"
                )

        for x in sdf.index:
            seen_samples.add(x)

        if reference_cols is None:
            reference_cols = pd.Index(sdf.columns.astype(str))
            sdf.columns = reference_cols
            sdf.to_csv(out_path, sep=args.sep, mode="w", header=True, index=True,
                       float_format=args.float_format)
        else:
            sdf.columns = sdf.columns.astype(str)
            cur_cols = pd.Index(sdf.columns)
            if not cur_cols.equals(reference_cols):
                if set(cur_cols) == set(reference_cols):
                    log("[WARN] column order differs; reordering to first batch columns")
                    sdf = sdf.loc[:, reference_cols]
                else:
                    missing = list(reference_cols.difference(cur_cols))[:10]
                    extra = list(cur_cols.difference(reference_cols))[:10]
                    raise ValueError(
                        f"CpG columns differ in {f}. "
                        f"missing first10={missing}; extra first10={extra}"
                    )

            sdf.to_csv(out_path, sep=args.sep, mode="a", header=False, index=True,
                       float_format=args.float_format)

        total_samples += sdf.shape[0]
        log(f"  wrote {sdf.shape[0]} samples x {sdf.shape[1]} CpGs")

    report = pd.DataFrame(report_rows)
    report.to_csv(report_path, index=False)

    if args.strict_renaming and (not args.drop_unmapped) and n_unmapped > 0:
        raise RuntimeError(
            f"{n_unmapped} samples were not mapped to sample sheet. "
            f"See {report_path}"
        )

    with open(summary_path, "w") as fh:
        fh.write(f"batch_root\t{batch_root}\n")
        fh.write(f"sample_sheet\t{sample_sheet}\n")
        fh.write(f"batch_files\t{len(files)}\n")
        fh.write(f"samples\t{total_samples}\n")
        fh.write(f"cpgs\t{len(reference_cols) if reference_cols is not None else 0}\n")
        fh.write(f"unmapped_samples\t{n_unmapped}\n")
        fh.write(f"output\t{out_path}\n")
        fh.write(f"renaming_report\t{report_path}\n")

    log("[Done]")
    log(f"  samples: {total_samples}")
    log(f"  CpGs: {len(reference_cols) if reference_cols is not None else 0}")
    log(f"  unmapped samples: {n_unmapped}")
    log(f"  matrix: {out_path}")
    log(f"  report: {report_path}")


if __name__ == "__main__":
    main()
