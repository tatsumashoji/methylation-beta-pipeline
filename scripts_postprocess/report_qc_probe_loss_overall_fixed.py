#!/usr/bin/env python3
"""Create corrected overall QC probe-loss and clock-overlap reports.

This script recomputes the overall EPICv2/MSA QC report from per-batch
qc_probe_failures_by_sample.tsv.gz files and the integrated collapsed beta
matrix. It fixes two common problems in older reports:

1. sample-rows collapsed beta matrices were accidentally interpreted as
   CpG-rows matrices, leading to integrated_matrix_cpg_n being equal to
   sample count.
2. QC loss is derived from per-sample pre/post QC beta comparisons rather
   than from probe names alone.

Definitions:
- qc_affected: CpG failed QC in at least one sample.
- qc_completely_lost: CpG failed QC in all matched samples in the dataset.
- qc_partially_affected: CpG failed QC in at least one but not all samples.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Set, Tuple


CLOCK_MODELS = {
    "Horvath": "Horvathv1",
    "Hannum": "Hannum",
    "PhenoAge": "PhenoAge",
    "GrimAgeV1": "GrimAgeV1",
    "GrimAgeV2": "GrimAgeV2",
    "DunedinPACE": "DunedinPACE",
}


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def parse_int(value, default: int = 0) -> int:
    if value is None:
        return default
    s = str(value).strip()
    if not s or s.lower() == "nan":
        return default
    try:
        return int(float(s))
    except Exception:
        return default


def parse_bool(value) -> bool:
    if value is None:
        return False
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def looks_like_cpg(value: str) -> bool:
    value = str(value).strip()
    return value.startswith("cg") or value.startswith("ch") or value.startswith("rs")


def read_probe_list(path: Path) -> Set[str]:
    probes: Set[str] = set()
    with open_text(path, "rt") as fh:
        for line in fh:
            value = line.strip().split("\t")[0].split(",")[0]
            if not value:
                continue
            if value.lower() in {"probeid", "probe_id", "cpg", "cpg_id"}:
                continue
            probes.add(value)
    return probes


def read_final_beta_probe_ids(path: Path, sep: str = "\t") -> Tuple[Set[str], str, int]:
    """Return CpG IDs, orientation, and sample count from collapsed beta matrix.

    Supports both:
    - CpG rows x sample columns: first column is ProbeID/CpG.
    - Sample rows x CpG columns: first column is Sample_ID and header columns are CpGs.
    """
    with open_text(path, "rt") as fh:
        header_line = fh.readline().rstrip("\n\r")
        if not header_line:
            return set(), "empty", 0
        header = header_line.split(sep)
        second_line = fh.readline().rstrip("\n\r")
        second = second_line.split(sep) if second_line else []

    header_cpg_count = sum(1 for x in header[1:] if looks_like_cpg(x))
    first_data_is_cpg = bool(second) and looks_like_cpg(second[0])

    if header_cpg_count >= max(10, int(0.5 * max(1, len(header[1:])))):
        # sample rows x CpG columns
        sample_count = 0
        with open_text(path, "rt") as fh:
            next(fh, None)
            for line in fh:
                if line.strip():
                    sample_count += 1
        return set(x.strip() for x in header[1:] if x.strip()), "sample_rows", sample_count

    if first_data_is_cpg:
        probes: Set[str] = set()
        sample_count = max(0, len(header) - 1)
        with open_text(path, "rt") as fh:
            next(fh, None)
            for line in fh:
                if not line.strip():
                    continue
                probe = line.split(sep, 1)[0].strip()
                if probe:
                    probes.add(probe)
        return probes, "cpg_rows", sample_count

    # Fallback: count first column as rows but mark unknown.
    probes = set()
    with open_text(path, "rt") as fh:
        next(fh, None)
        for line in fh:
            if line.strip():
                probes.add(line.split(sep, 1)[0].strip())
    return probes, "unknown", max(0, len(header) - 1)


def load_clock_sites() -> Dict[str, Set[str]]:
    from biolearn.model_gallery import ModelGallery

    gallery = ModelGallery()
    out: Dict[str, Set[str]] = {}
    for clock, model_name in CLOCK_MODELS.items():
        model = gallery.get(model_name, imputation_method="none")
        if not hasattr(model, "methylation_sites"):
            raise RuntimeError(f"Biolearn model {model_name} lacks methylation_sites().")
        out[clock] = {str(x).strip() for x in model.methylation_sites() if str(x).strip()}
    return out


def find_column(fieldnames: Iterable[str], candidates: Iterable[str]) -> str | None:
    lower_map = {str(f).strip().lower(): f for f in fieldnames}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def read_batch_summary(path: Path) -> dict:
    if not path.is_file():
        return {}
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            return row
    return {}


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_tsv_gz(path: Path, rows_iter, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows_iter:
            writer.writerow(row)


def process_by_sample_qc(
    qc_path: Path,
    sample_ids: Set[str],
    failed_sample_counter: Counter,
    affected_batch_counter: Counter,
    probe_to_samples: dict[str, Set[str]],
) -> None:
    with open_text(qc_path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return
        sample_col = find_column(reader.fieldnames, ["SampleID", "sample_id", "sample", "Sample_ID"])
        probe_col = find_column(reader.fieldnames, ["ProbeID", "probe_id", "probe", "CpG", "cpg_id"])
        if sample_col is None or probe_col is None:
            raise ValueError(f"{qc_path} lacks SampleID/ProbeID columns. Columns={reader.fieldnames}")
        batch_seen_probes: Set[str] = set()
        for row in reader:
            sample = str(row.get(sample_col, "")).strip()
            probe = str(row.get(probe_col, "")).strip()
            if not sample or not probe:
                continue
            sample_ids.add(sample)
            probe_to_samples[probe].add(sample)
            failed_sample_counter[probe] += 1
            batch_seen_probes.add(probe)
        for probe in batch_seen_probes:
            affected_batch_counter[probe] += 1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--batch-root", type=Path, required=True)
    ap.add_argument("--final-beta", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--sep", default="\t")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    batch_dirs = sorted(args.batch_root.glob("batch_*"))
    if not batch_dirs:
        raise FileNotFoundError(f"No batch directories found under {args.batch_root}")

    pre_qc_union: Set[str] = set()
    post_qc_union: Set[str] = set()
    sample_ids: Set[str] = set()
    failed_sample_counter: Counter = Counter()
    affected_batch_counter: Counter = Counter()
    probe_to_samples: dict[str, Set[str]] = defaultdict(set)
    total_samples_from_summary = 0

    for idx, batch_dir in enumerate(batch_dirs, start=1):
        print(f"[overall-qc] {args.dataset}: reading {batch_dir.name} ({idx}/{len(batch_dirs)})", flush=True)
        pre_path = batch_dir / "pre_qc_collapsed_probe_ids.txt.gz"
        post_path = batch_dir / "post_qc_collapsed_probe_ids.txt.gz"
        qc_by_sample = batch_dir / "qc_probe_failures_by_sample.tsv.gz"
        qc_agg = batch_dir / "qc_probe_failures.tsv.gz"
        summary_path = batch_dir / "qc_batch_summary.csv"

        if not pre_path.is_file():
            raise FileNotFoundError(f"Missing {pre_path}")
        pre_qc_union.update(read_probe_list(pre_path))
        if post_path.is_file():
            post_qc_union.update(read_probe_list(post_path))

        summary = read_batch_summary(summary_path)
        total_samples_from_summary += parse_int(summary.get("n_samples") or summary.get("sample_count"), default=0)

        if qc_by_sample.is_file():
            process_by_sample_qc(qc_by_sample, sample_ids, failed_sample_counter, affected_batch_counter, probe_to_samples)
        elif qc_agg.is_file():
            # Fallback for older reports: no sample-level info.
            with open_text(qc_agg, "rt") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                if reader.fieldnames is None:
                    continue
                probe_col = find_column(reader.fieldnames, ["ProbeID", "probe_id", "probe", "CpG", "cpg_id"])
                failed_col = find_column(reader.fieldnames, ["qc_failed_samples", "failed_sample_count"])
                if probe_col is None:
                    raise ValueError(f"{qc_agg} lacks probe column")
                for row in reader:
                    probe = str(row.get(probe_col, "")).strip()
                    if not probe:
                        continue
                    n = parse_int(row.get(failed_col), default=1) if failed_col else 1
                    failed_sample_counter[probe] += n
                    affected_batch_counter[probe] += 1
        else:
            raise FileNotFoundError(f"Missing QC file in {batch_dir}")

    n_samples = len(sample_ids) if sample_ids else total_samples_from_summary
    final_probe_ids, final_orientation, final_sample_count = read_final_beta_probe_ids(args.final_beta, sep=args.sep)
    if final_sample_count and not n_samples:
        n_samples = final_sample_count

    affected_cpgs = set(failed_sample_counter.keys())
    completely_lost_cpgs = {p for p in affected_cpgs if failed_sample_counter[p] >= n_samples and n_samples > 0}
    partially_affected_cpgs = affected_cpgs - completely_lost_cpgs
    post_qc_estimated = pre_qc_union - completely_lost_cpgs

    clock_sites = load_clock_sites()
    clock_names = list(CLOCK_MODELS.keys())

    def required_flags(probe: str) -> dict:
        required = []
        row = {}
        for clock in clock_names:
            hit = probe in clock_sites[clock]
            row[f"required_by_{clock}"] = hit
            if hit:
                required.append(clock)
        row["required_clock_count"] = len(required)
        row["required_clocks"] = ";".join(required)
        return row

    detail_fields = [
        "ProbeID", "qc_failed_samples", "affected_batch_count", "qc_failed_fraction",
        "qc_completely_lost", "qc_partially_affected", "present_before_qc",
        "present_after_qc_estimated", "present_in_integrated_matrix",
    ] + [f"required_by_{c}" for c in clock_names] + ["required_clock_count", "required_clocks"]

    def detail_rows(probes: Iterable[str]):
        for probe in sorted(probes):
            row = {
                "ProbeID": probe,
                "qc_failed_samples": failed_sample_counter.get(probe, 0),
                "affected_batch_count": affected_batch_counter.get(probe, 0),
                "qc_failed_fraction": (failed_sample_counter.get(probe, 0) / n_samples if n_samples else 0.0),
                "qc_completely_lost": probe in completely_lost_cpgs,
                "qc_partially_affected": probe in partially_affected_cpgs,
                "present_before_qc": probe in pre_qc_union,
                "present_after_qc_estimated": probe in post_qc_estimated,
                "present_in_integrated_matrix": probe in final_probe_ids,
            }
            row.update(required_flags(probe))
            yield row

    write_tsv_gz(args.out_dir / "qc_affected_cpgs.tsv.gz", detail_rows(affected_cpgs), detail_fields)
    write_tsv_gz(args.out_dir / "qc_completely_lost_cpgs.tsv.gz", detail_rows(completely_lost_cpgs), detail_fields)

    clock_summary_rows = []
    clock_detail_fields = [
        "clock", "biolearn_model", "ProbeID", "present_before_qc", "qc_affected",
        "qc_completely_lost", "qc_partially_affected", "qc_failed_samples",
        "affected_batch_count", "present_after_qc_estimated", "present_in_integrated_matrix",
    ]

    def clock_detail_rows():
        for clock, sites in clock_sites.items():
            for probe in sorted(sites):
                yield {
                    "clock": clock,
                    "biolearn_model": CLOCK_MODELS[clock],
                    "ProbeID": probe,
                    "present_before_qc": probe in pre_qc_union,
                    "qc_affected": probe in affected_cpgs,
                    "qc_completely_lost": probe in completely_lost_cpgs,
                    "qc_partially_affected": probe in partially_affected_cpgs,
                    "qc_failed_samples": failed_sample_counter.get(probe, 0),
                    "affected_batch_count": affected_batch_counter.get(probe, 0),
                    "present_after_qc_estimated": probe in post_qc_estimated,
                    "present_in_integrated_matrix": probe in final_probe_ids,
                }

    for clock, sites in clock_sites.items():
        required_n = len(sites)
        present_before = sites & pre_qc_union
        present_after = sites & post_qc_estimated
        in_final = sites & final_probe_ids
        qc_aff = sites & affected_cpgs
        qc_lost = sites & completely_lost_cpgs
        qc_partial = sites & partially_affected_cpgs
        clock_summary_rows.append({
            "dataset": args.dataset,
            "clock": clock,
            "biolearn_model": CLOCK_MODELS[clock],
            "required_cpg_n": required_n,
            "present_before_qc_n": len(present_before),
            "missing_before_qc_n": required_n - len(present_before),
            "qc_affected_cpg_n": len(qc_aff),
            "qc_affected_pct_of_required": 100.0 * len(qc_aff) / required_n if required_n else 0.0,
            "qc_completely_lost_cpg_n": len(qc_lost),
            "qc_completely_lost_pct_of_required": 100.0 * len(qc_lost) / required_n if required_n else 0.0,
            "qc_partially_affected_cpg_n": len(qc_partial),
            "present_after_qc_n": len(present_after),
            "missing_after_qc_n": required_n - len(present_after),
            "present_in_integrated_matrix_n": len(in_final),
            "absent_from_integrated_matrix_n": required_n - len(in_final),
        })

    write_csv(
        args.out_dir / "clock_qc_overlap_summary.csv",
        clock_summary_rows,
        [
            "dataset", "clock", "biolearn_model", "required_cpg_n",
            "present_before_qc_n", "missing_before_qc_n", "qc_affected_cpg_n",
            "qc_affected_pct_of_required", "qc_completely_lost_cpg_n",
            "qc_completely_lost_pct_of_required", "qc_partially_affected_cpg_n",
            "present_after_qc_n", "missing_after_qc_n", "present_in_integrated_matrix_n",
            "absent_from_integrated_matrix_n",
        ],
    )
    write_tsv_gz(args.out_dir / "clock_qc_overlap_cpgs.tsv.gz", clock_detail_rows(), clock_detail_fields)

    summary_rows = [
        {"metric": "dataset", "value": args.dataset},
        {"metric": "batch_count", "value": len(batch_dirs)},
        {"metric": "sample_count", "value": n_samples},
        {"metric": "sample_count_from_batch_summaries", "value": total_samples_from_summary},
        {"metric": "final_beta_orientation", "value": final_orientation},
        {"metric": "final_beta_sample_n", "value": final_sample_count},
        {"metric": "pre_qc_cpg_n", "value": len(pre_qc_union)},
        {"metric": "post_qc_cpg_n_estimated", "value": len(post_qc_estimated)},
        {"metric": "integrated_matrix_cpg_n", "value": len(final_probe_ids)},
        {"metric": "qc_affected_cpg_n", "value": len(affected_cpgs)},
        {"metric": "qc_completely_lost_cpg_n", "value": len(completely_lost_cpgs)},
        {"metric": "qc_partially_affected_cpg_n", "value": len(partially_affected_cpgs)},
    ]
    write_csv(args.out_dir / "qc_probe_loss_summary.csv", summary_rows, ["metric", "value"])

    metadata = {
        "dataset": args.dataset,
        "batch_root": str(args.batch_root),
        "final_beta": str(args.final_beta),
        "out_dir": str(args.out_dir),
        "definition": {
            "qc_affected": "CpG failed QC in at least one sample.",
            "qc_completely_lost": "CpG failed QC in all matched samples.",
            "present_after_qc_n": "required CpGs present before QC and not completely lost by QC.",
        },
    }
    with open(args.out_dir / "qc_report_metadata.json", "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2)

    print(f"[Done] dataset: {args.dataset}", flush=True)
    print(f"[Done] output: {args.out_dir}", flush=True)
    print(f"[Done] final beta orientation: {final_orientation}", flush=True)
    print(f"[Done] pre-QC CpGs: {len(pre_qc_union):,}", flush=True)
    print(f"[Done] QC-affected CpGs: {len(affected_cpgs):,}", flush=True)
    print(f"[Done] QC-completely-lost CpGs: {len(completely_lost_cpgs):,}", flush=True)


if __name__ == "__main__":
    main()
