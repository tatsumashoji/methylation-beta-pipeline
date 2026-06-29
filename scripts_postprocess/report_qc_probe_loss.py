#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Set


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


def parse_bool(value) -> bool:
    if value is None:
        return False
    s = str(value).strip().lower()
    return s in {"true", "t", "1", "yes", "y"}


def parse_int(value, default: int = 0) -> int:
    if value is None:
        return default
    s = str(value).strip()
    if s == "" or s.lower() == "nan":
        return default
    try:
        return int(float(s))
    except Exception:
        return default


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


def read_final_beta_probe_ids(path: Path, sep: str = "\t") -> Set[str]:
    probes: Set[str] = set()
    with open_text(path, "rt") as fh:
        header = fh.readline()
        for line in fh:
            if not line.strip():
                continue
            probe = line.split(sep, 1)[0].strip()
            if probe:
                probes.add(probe)
    return probes


def load_clock_sites() -> Dict[str, Set[str]]:
    from biolearn.model_gallery import ModelGallery

    gallery = ModelGallery()
    result: Dict[str, Set[str]] = {}

    for clock_name, model_name in CLOCK_MODELS.items():
        model = gallery.get(model_name, imputation_method="none")
        if not hasattr(model, "methylation_sites"):
            raise RuntimeError(
                f"Biolearn model {model_name} does not expose methylation_sites()."
            )
        result[clock_name] = {
            str(x).strip()
            for x in model.methylation_sites()
            if str(x).strip()
        }

    return result


def find_column(fieldnames: Iterable[str], candidates: Iterable[str]) -> str | None:
    field_set = {f.lower(): f for f in fieldnames}
    for cand in candidates:
        if cand.lower() in field_set:
            return field_set[cand.lower()]
    return None


def read_batch_summary(path: Path) -> dict:
    if not path.is_file():
        return {}

    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            return row

    return {}


def process_qc_failures(
    qc_path: Path,
    affected_cpgs: Set[str],
    completely_lost_cpgs: Set[str],
    partially_affected_cpgs: Set[str],
    failed_sample_counter: Counter,
    affected_batch_counter: Counter,
) -> None:
    with open_text(qc_path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return

        probe_col = find_column(
            reader.fieldnames,
            ["ProbeID", "probe_id", "probe", "CpG", "cpg_id"],
        )
        failed_samples_col = find_column(
            reader.fieldnames,
            [
                "qc_failed_samples",
                "failed_sample_count",
                "failed_samples",
                "n_failed_samples",
            ],
        )
        completely_lost_col = find_column(
            reader.fieldnames,
            [
                "qc_completely_lost",
                "completely_lost",
                "dropped_from_batch_collapsed",
            ],
        )
        partially_affected_col = find_column(
            reader.fieldnames,
            [
                "qc_partially_affected",
                "partially_affected",
            ],
        )

        if probe_col is None:
            raise ValueError(f"{qc_path} has no probe ID column. Columns: {reader.fieldnames}")

        for row in reader:
            probe = str(row.get(probe_col, "")).strip()
            if not probe:
                continue

            affected_cpgs.add(probe)
            affected_batch_counter[probe] += 1

            failed_n = parse_int(row.get(failed_samples_col), default=1) if failed_samples_col else 1
            failed_sample_counter[probe] += failed_n

            if completely_lost_col and parse_bool(row.get(completely_lost_col)):
                completely_lost_cpgs.add(probe)

            if partially_affected_col and parse_bool(row.get(partially_affected_col)):
                partially_affected_cpgs.add(probe)


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--batch-root", required=True, type=Path)
    parser.add_argument("--final-beta", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--sep", default="\t")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    batch_dirs = sorted(args.batch_root.glob("batch_*"))
    if not batch_dirs:
        raise FileNotFoundError(f"No batch directories found under {args.batch_root}")

    pre_qc_union: Set[str] = set()
    post_qc_union: Set[str] = set()
    affected_cpgs: Set[str] = set()
    completely_lost_cpgs_from_batches: Set[str] = set()
    partially_affected_cpgs: Set[str] = set()
    failed_sample_counter: Counter = Counter()
    affected_batch_counter: Counter = Counter()

    total_samples = 0
    batch_count = 0

    for batch_dir in batch_dirs:
        batch_count += 1

        pre_path = batch_dir / "pre_qc_collapsed_probe_ids.txt.gz"
        post_path = batch_dir / "post_qc_collapsed_probe_ids.txt.gz"
        qc_path = batch_dir / "qc_probe_failures.tsv.gz"
        summary_path = batch_dir / "qc_batch_summary.csv"

        if not pre_path.is_file():
            raise FileNotFoundError(f"Missing {pre_path}")
        if not post_path.is_file():
            raise FileNotFoundError(f"Missing {post_path}")
        if not qc_path.is_file():
            raise FileNotFoundError(f"Missing {qc_path}")

        pre_qc_union.update(read_probe_list(pre_path))
        post_qc_union.update(read_probe_list(post_path))

        summary = read_batch_summary(summary_path)
        total_samples += parse_int(
            summary.get("n_samples") or summary.get("sample_count"),
            default=0,
        )

        process_qc_failures(
            qc_path=qc_path,
            affected_cpgs=affected_cpgs,
            completely_lost_cpgs=completely_lost_cpgs_from_batches,
            partially_affected_cpgs=partially_affected_cpgs,
            failed_sample_counter=failed_sample_counter,
            affected_batch_counter=affected_batch_counter,
        )

    final_probe_ids = read_final_beta_probe_ids(args.final_beta, sep=args.sep)

    # Completely lost at integrated matrix level:
    # present before QC, affected by QC, and absent after QC / absent from final matrix.
    completely_lost_cpgs = {
        p for p in affected_cpgs
        if p in pre_qc_union and p not in post_qc_union
    }

    completely_lost_cpgs |= {
        p for p in affected_cpgs
        if p in pre_qc_union and p not in final_probe_ids
    }

    partially_affected_cpgs |= affected_cpgs - completely_lost_cpgs

    clock_sites = load_clock_sites()
    all_clock_names = list(CLOCK_MODELS.keys())

    def required_by_flags(probe: str) -> dict:
        flags = {}
        required_clocks = []
        for clock in all_clock_names:
            hit = probe in clock_sites[clock]
            flags[f"required_by_{clock}"] = hit
            if hit:
                required_clocks.append(clock)
        flags["required_clock_count"] = len(required_clocks)
        flags["required_clocks"] = ";".join(required_clocks)
        return flags

    base_detail_fields = [
        "ProbeID",
        "qc_failed_samples",
        "affected_batch_count",
        "qc_completely_lost",
        "qc_partially_affected",
        "present_before_qc",
        "present_after_qc",
        "present_in_integrated_matrix",
    ]
    clock_flag_fields = []
    for clock in all_clock_names:
        clock_flag_fields.append(f"required_by_{clock}")
    detail_fields = base_detail_fields + clock_flag_fields + [
        "required_clock_count",
        "required_clocks",
    ]

    def detail_rows(probes: Iterable[str]):
        for probe in sorted(probes):
            row = {
                "ProbeID": probe,
                "qc_failed_samples": failed_sample_counter.get(probe, 0),
                "affected_batch_count": affected_batch_counter.get(probe, 0),
                "qc_completely_lost": probe in completely_lost_cpgs,
                "qc_partially_affected": probe in partially_affected_cpgs,
                "present_before_qc": probe in pre_qc_union,
                "present_after_qc": probe in post_qc_union,
                "present_in_integrated_matrix": probe in final_probe_ids,
            }
            row.update(required_by_flags(probe))
            yield row

    write_tsv_gz(
        args.out_dir / "qc_affected_cpgs.tsv.gz",
        detail_rows(affected_cpgs),
        detail_fields,
    )

    write_tsv_gz(
        args.out_dir / "qc_completely_lost_cpgs.tsv.gz",
        detail_rows(completely_lost_cpgs),
        detail_fields,
    )

    clock_summary_rows = []
    clock_detail_fields = [
        "clock",
        "biolearn_model",
        "ProbeID",
        "present_before_qc",
        "qc_affected",
        "qc_completely_lost",
        "qc_partially_affected",
        "qc_failed_samples",
        "affected_batch_count",
        "present_after_qc",
        "present_in_integrated_matrix",
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
                    "present_after_qc": probe in post_qc_union,
                    "present_in_integrated_matrix": probe in final_probe_ids,
                }

    for clock, sites in clock_sites.items():
        required_n = len(sites)
        qc_affected = sites & affected_cpgs
        qc_lost = sites & completely_lost_cpgs
        qc_partial = sites & partially_affected_cpgs

        row = {
            "dataset": args.dataset,
            "clock": clock,
            "biolearn_model": CLOCK_MODELS[clock],
            "required_cpg_n": required_n,
            "present_before_qc_n": len(sites & pre_qc_union),
            "missing_before_qc_n": len(sites - pre_qc_union),
            "qc_affected_cpg_n": len(qc_affected),
            "qc_affected_pct_of_required": (
                100.0 * len(qc_affected) / required_n if required_n else 0.0
            ),
            "qc_completely_lost_cpg_n": len(qc_lost),
            "qc_completely_lost_pct_of_required": (
                100.0 * len(qc_lost) / required_n if required_n else 0.0
            ),
            "qc_partially_affected_cpg_n": len(qc_partial),
            "present_after_qc_n": len(sites & post_qc_union),
            "missing_after_qc_n": len(sites - post_qc_union),
            "present_in_integrated_matrix_n": len(sites & final_probe_ids),
            "absent_from_integrated_matrix_n": len(sites - final_probe_ids),
        }
        clock_summary_rows.append(row)

    write_csv(
        args.out_dir / "clock_qc_overlap_summary.csv",
        clock_summary_rows,
        [
            "dataset",
            "clock",
            "biolearn_model",
            "required_cpg_n",
            "present_before_qc_n",
            "missing_before_qc_n",
            "qc_affected_cpg_n",
            "qc_affected_pct_of_required",
            "qc_completely_lost_cpg_n",
            "qc_completely_lost_pct_of_required",
            "qc_partially_affected_cpg_n",
            "present_after_qc_n",
            "missing_after_qc_n",
            "present_in_integrated_matrix_n",
            "absent_from_integrated_matrix_n",
        ],
    )

    write_tsv_gz(
        args.out_dir / "clock_qc_overlap_cpgs.tsv.gz",
        clock_detail_rows(),
        clock_detail_fields,
    )

    pipeline_summary_rows = [
        {"metric": "dataset", "value": args.dataset},
        {"metric": "batch_count", "value": batch_count},
        {"metric": "sample_count_from_batch_summaries", "value": total_samples},
        {"metric": "pre_qc_cpg_n", "value": len(pre_qc_union)},
        {"metric": "post_qc_cpg_n", "value": len(post_qc_union)},
        {"metric": "integrated_matrix_cpg_n", "value": len(final_probe_ids)},
        {"metric": "qc_affected_cpg_n", "value": len(affected_cpgs)},
        {"metric": "qc_completely_lost_cpg_n", "value": len(completely_lost_cpgs)},
        {"metric": "qc_partially_affected_cpg_n", "value": len(partially_affected_cpgs)},
    ]

    write_csv(
        args.out_dir / "qc_probe_loss_summary.csv",
        pipeline_summary_rows,
        ["metric", "value"],
    )

    metadata = {
        "dataset": args.dataset,
        "batch_root": str(args.batch_root),
        "final_beta": str(args.final_beta),
        "out_dir": str(args.out_dir),
        "batch_count": batch_count,
        "note": "Memory-efficient streaming QC probe-loss report.",
    }

    with open(args.out_dir / "qc_report_metadata.json", "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2)

    print(f"[Done] dataset: {args.dataset}")
    print(f"[Done] output: {args.out_dir}")
    print(f"[Done] QC-affected CpGs: {len(affected_cpgs):,}")
    print(f"[Done] QC-completely-lost CpGs: {len(completely_lost_cpgs):,}")


if __name__ == "__main__":
    main()
