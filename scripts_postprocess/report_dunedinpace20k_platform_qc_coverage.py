#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path
from typing import Iterable

import pandas as pd
from biolearn.model_gallery import ModelGallery
from biolearn.util import get_data_file


CG_RE = re.compile(r"cg\d{8}", re.IGNORECASE)
NA_VALUES = {"", "NA", "NaN", "nan", "NAN", "None", "none", "null", "NULL"}


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def normalize_cpg(x) -> str | None:
    if x is None:
        return None
    m = CG_RE.search(str(x))
    if not m:
        return None
    return m.group(0).lower()


def load_dunedinpace_sets() -> dict[str, set[str]]:
    gallery = ModelGallery()
    model = gallery.get("DunedinPACE", imputation_method="none")

    core = {normalize_cpg(x) for x in model.methylation_sites()}
    core = {x for x in core if x}

    gold = pd.read_csv(
        get_data_file("DunedinPACE_Gold_Means.csv"),
        index_col=0,
    )
    gold_set = {normalize_cpg(x) for x in gold.index.astype(str)}
    gold_set = {x for x in gold_set if x}

    background = gold_set - core
    union = core | gold_set

    return {
        "core_173": core,
        "gold_20k": gold_set,
        "background_only_19827": background,
        "union_requested": union,
    }


def scan_manifest(manifest_path: Path, requested: set[str]) -> set[str]:
    found = set()

    with open_text(manifest_path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            for m in CG_RE.finditer(line):
                cg = m.group(0).lower()
                if cg in requested:
                    found.add(cg)

            if i % 1_000_000 == 0:
                print(
                    f"[Manifest] scanned {i:,} lines; "
                    f"found {len(found):,}/{len(requested):,}",
                    flush=True,
                )

    return found


def read_probe_list(path: Path) -> set[str]:
    result = set()

    with open_text(path, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            first = line.strip().split("\t")[0].split(",")[0]
            cg = normalize_cpg(first)
            if cg:
                result.add(cg)

    return result


def find_existing(batch_dir: Path, names: list[str]) -> Path | None:
    for name in names:
        p = batch_dir / name
        if p.is_file():
            return p
    return None


def collect_pre_post_and_qc(batch_root: Path) -> tuple[set[str], set[str], set[str], int]:
    pre_union = set()
    post_union = set()
    qc_affected = set()
    batch_count = 0

    batch_dirs = sorted(batch_root.glob("batch_*"))
    if not batch_dirs:
        raise FileNotFoundError(f"No batch_* directories under {batch_root}")

    for batch_dir in batch_dirs:
        batch_count += 1

        pre_path = find_existing(
            batch_dir,
            ["pre_qc_collapsed_probe_ids.txt.gz", "pre_qc_probe_ids.txt.gz"],
        )
        post_path = find_existing(
            batch_dir,
            ["post_qc_collapsed_probe_ids.txt.gz", "post_qc_probe_ids.txt.gz"],
        )
        qc_path = find_existing(
            batch_dir,
            ["qc_probe_failures_by_sample.tsv.gz", "qc_probe_failures.tsv.gz"],
        )

        if pre_path is None:
            raise FileNotFoundError(f"Pre-QC probe list not found in {batch_dir}")
        if post_path is None:
            raise FileNotFoundError(f"Post-QC probe list not found in {batch_dir}")

        pre_set = read_probe_list(pre_path)
        post_set = read_probe_list(post_path)

        pre_union.update(pre_set)
        post_union.update(post_set)

        if qc_path is not None:
            with open_text(qc_path, "rt") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                if reader.fieldnames:
                    probe_col = None
                    for cand in ["ProbeID", "probe_id", "CpG", "cpg_id"]:
                        if cand in reader.fieldnames:
                            probe_col = cand
                            break
                    if probe_col is not None:
                        for row in reader:
                            cg = normalize_cpg(row.get(probe_col))
                            if cg:
                                qc_affected.add(cg)

        print(
            f"[Batch] {batch_dir.name}: "
            f"pre={len(pre_set):,}, post={len(post_set):,}",
            flush=True,
        )

    return pre_union, post_union, qc_affected, batch_count


def infer_orientation(header: list[str], first_row: list[str]) -> str:
    header_cg_count = sum(
        1 for x in header[1:min(len(header), 2000)]
        if normalize_cpg(x)
    )
    first_col_is_cg = normalize_cpg(first_row[0]) is not None if first_row else False

    if header_cg_count > 10:
        return "sample_rows"
    if first_col_is_cg:
        return "cpg_rows"

    raise ValueError(
        "Could not infer final matrix orientation. "
        f"Header head: {header[:5]}, first row head: {first_row[:5]}"
    )


def is_missing(x: str) -> bool:
    return str(x).strip() in NA_VALUES


def read_final_matrix_presence(final_beta: Path, requested: set[str], sep: str = "\t") -> dict:
    matrix_features = set()
    requested_features = set()
    nonmissing_counter = Counter()
    n_samples = 0

    with open_text(final_beta, "rt") as fh:
        reader = csv.reader(fh, delimiter=sep)
        header = next(reader)
        first_row = next(reader, None)

        if first_row is None:
            raise ValueError(f"Empty final beta matrix: {final_beta}")

        orientation = infer_orientation(header, first_row)

        if orientation == "sample_rows":
            requested_col_index = {}

            for idx, col in enumerate(header):
                cg = normalize_cpg(col)
                if cg:
                    matrix_features.add(cg)
                    if cg in requested:
                        requested_col_index[idx] = cg
                        requested_features.add(cg)

            def process_sample(row: list[str]):
                nonlocal n_samples
                n_samples += 1
                for idx, cg in requested_col_index.items():
                    if idx < len(row) and not is_missing(row[idx]):
                        nonmissing_counter[cg] += 1

            process_sample(first_row)

            for row in reader:
                process_sample(row)
                if n_samples % 25 == 0:
                    print(f"[Final] processed samples: {n_samples:,}", flush=True)

        else:
            n_samples = max(len(header) - 1, 0)

            def process_cpg(row: list[str]):
                if not row:
                    return
                cg = normalize_cpg(row[0])
                if not cg:
                    return

                matrix_features.add(cg)

                if cg not in requested:
                    return

                requested_features.add(cg)
                values = row[1:]
                nonmissing_counter[cg] += sum(
                    1 for v in values if not is_missing(v)
                )

            process_cpg(first_row)

            for i, row in enumerate(reader, start=2):
                process_cpg(row)
                if i % 100_000 == 0:
                    print(f"[Final] processed CpG rows: {i:,}", flush=True)

    final_any = {
        cg for cg in requested_features
        if nonmissing_counter.get(cg, 0) > 0
    }
    final_all = {
        cg for cg in requested_features
        if n_samples > 0 and nonmissing_counter.get(cg, 0) == n_samples
    }

    return {
        "orientation": orientation,
        "n_samples": n_samples,
        "matrix_features": matrix_features,
        "requested_features": requested_features,
        "final_any": final_any,
        "final_all": final_all,
    }


def write_list(path: Path, values: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        for v in sorted(values):
            fh.write(f"{v}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, choices=["EPICv2", "MSA"])
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--batch-root", required=True, type=Path)
    parser.add_argument("--final-beta", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--sep", default="\t")

    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    probe_sets = load_dunedinpace_sets()
    requested_union = probe_sets["union_requested"]

    print(f"[DunedinPACE20kQC] dataset: {args.dataset}")
    for k, v in probe_sets.items():
        print(f"[DunedinPACE20kQC] {k}: {len(v):,}")

    print("[DunedinPACE20kQC] scanning manifest...")
    manifest_present = scan_manifest(args.manifest, requested_union)

    print("[DunedinPACE20kQC] collecting pre/post QC probe lists...")
    pre_qc, post_qc, qc_affected, batch_count = collect_pre_post_and_qc(args.batch_root)

    print("[DunedinPACE20kQC] reading final matrix...")
    final = read_final_matrix_presence(
        args.final_beta,
        requested_union,
        sep=args.sep,
    )

    rows = []

    for set_name, req in probe_sets.items():
        row = {
            "dataset": args.dataset,
            "probe_set": set_name,
            "requested_n": len(req),

            "manifest_present_n": len(req & manifest_present),
            "manifest_missing_n": len(req - manifest_present),

            "pre_qc_present_n": len(req & pre_qc),
            "pre_qc_missing_n": len(req - pre_qc),

            "post_qc_present_n": len(req & post_qc),
            "post_qc_missing_n": len(req - post_qc),

            "qc_affected_at_least_one_sample_n": len(req & qc_affected),
            "qc_affected_pct_of_requested": (
                100 * len(req & qc_affected) / len(req) if len(req) else 0
            ),

            "lost_pre_to_post_qc_n": len((req & pre_qc) - post_qc),
            "lost_pre_to_final_any_n": len((req & pre_qc) - final["final_any"]),

            "final_matrix_feature_present_n": len(req & final["matrix_features"]),
            "final_value_present_any_sample_n": len(req & final["final_any"]),
            "final_value_present_all_samples_n": len(req & final["final_all"]),
            "final_value_missing_all_samples_n": len(req - final["final_any"]),
            "final_value_missing_at_least_one_sample_n": len(req - final["final_all"]),

            "batch_count": batch_count,
            "final_beta_orientation": final["orientation"],
            "final_beta_sample_n": final["n_samples"],
        }

        rows.append(row)

        prefix = f"{args.dataset}.{set_name}"

        write_list(args.out_dir / f"{prefix}.present_in_manifest.txt.gz", req & manifest_present)
        write_list(args.out_dir / f"{prefix}.missing_from_manifest.txt.gz", req - manifest_present)
        write_list(args.out_dir / f"{prefix}.present_pre_qc.txt.gz", req & pre_qc)
        write_list(args.out_dir / f"{prefix}.present_post_qc.txt.gz", req & post_qc)
        write_list(args.out_dir / f"{prefix}.qc_affected.txt.gz", req & qc_affected)
        write_list(args.out_dir / f"{prefix}.present_final_any_sample.txt.gz", req & final["final_any"])
        write_list(args.out_dir / f"{prefix}.present_final_all_samples.txt.gz", req & final["final_all"])
        write_list(args.out_dir / f"{prefix}.lost_pre_to_post_qc.txt.gz", (req & pre_qc) - post_qc)
        write_list(args.out_dir / f"{prefix}.lost_pre_to_final_any.txt.gz", (req & pre_qc) - final["final_any"])

    summary_path = args.out_dir / f"{args.dataset}.DunedinPACE20k_platform_qc_coverage.csv"
    pd.DataFrame(rows).to_csv(summary_path, index=False)

    metadata = {
        "dataset": args.dataset,
        "manifest": str(args.manifest),
        "batch_root": str(args.batch_root),
        "final_beta": str(args.final_beta),
        "out_dir": str(args.out_dir),
        "note": "DunedinPACE 20k platform/QC/final-matrix coverage report.",
    }

    with open(args.out_dir / f"{args.dataset}.DunedinPACE20k_platform_qc_coverage.metadata.json", "w") as fh:
        json.dump(metadata, fh, indent=2)

    print("[DunedinPACE20kQC] Wrote:", summary_path)


if __name__ == "__main__":
    main()
