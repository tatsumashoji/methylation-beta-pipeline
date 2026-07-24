#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

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


def canonical_id(x) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    s = re.sub(r"([_\-\.\s]*)blood$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\.0$", "", s)
    s = s.upper()
    s = re.sub(r"[^A-Z0-9]", "", s)
    return s


def read_any_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def choose_sample_col(df: pd.DataFrame) -> str:
    candidates = [
        "Name",
        "Sample_ID",
        "SampleID",
        "sample_id",
        "sample",
        "ID",
        "Name_Original",
    ]
    for col in candidates:
        if col in df.columns:
            return col

    for col in df.columns:
        if not str(col).startswith("Unnamed"):
            return col

    return df.columns[0]


def load_split_ids(path: Path) -> set[str]:
    df = read_any_csv(path)
    col = choose_sample_col(df)
    return {
        str(x).strip()
        for x in df[col].dropna()
        if str(x).strip() and str(x).strip().lower() != col.lower()
    }


def load_dunedinpace_sets() -> dict[str, set[str]]:
    gallery = ModelGallery()
    model = gallery.get("DunedinPACE", imputation_method="none")

    core = {normalize_cpg(x) for x in model.methylation_sites()}
    core = {x for x in core if x}

    gold = pd.read_csv(get_data_file("DunedinPACE_Gold_Means.csv"), index_col=0)
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
    found: set[str] = set()

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


def raw_id_from_idat_path(path_value: str) -> str:
    name = Path(str(path_value)).name
    name = re.sub(r"_(Grn|Red)\.idat(\.gz)?$", "", name, flags=re.IGNORECASE)
    return name


def load_raw_to_final_sample_map(sample_sheet: Path) -> dict[str, str]:
    df = read_any_csv(sample_sheet)
    raw_to_final: dict[str, str] = {}

    final_candidates = [
        "Sample_Name2",
        "Sample_Name1",
        "Sample_ID",
        "SampleID",
        "Name",
        "sample_id",
    ]

    final_col = None
    for c in final_candidates:
        if c in df.columns:
            final_col = c
            break

    if final_col is None:
        final_col = df.columns[0]

    path_cols = [c for c in df.columns if c.lower() in {"grn", "red"}]
    if not path_cols:
        path_cols = [
            c for c in df.columns
            if df[c].astype(str).str.contains(r"\.idat", case=False, regex=True).any()
        ]

    for _, row in df.iterrows():
        final_id = str(row.get(final_col, "")).strip()
        if not final_id:
            continue

        for c in path_cols:
            raw = raw_id_from_idat_path(row.get(c, ""))
            if raw:
                raw_to_final[raw] = final_id

    return raw_to_final


def load_optional_map(path: Path | None, final_ids: set[str]) -> dict[str, str]:
    if path is None or not path.is_file():
        return {}

    df = read_any_csv(path)
    final_canon = {canonical_id(x): x for x in final_ids}

    result: dict[str, str] = {}

    for source_col in df.columns:
        for target_col in df.columns:
            if source_col == target_col:
                continue

            for _, row in df.iterrows():
                source = str(row.get(source_col, "")).strip()
                target = str(row.get(target_col, "")).strip()

                if not source or not target:
                    continue

                target_final = None
                if target in final_ids:
                    target_final = target
                elif canonical_id(target) in final_canon:
                    target_final = final_canon[canonical_id(target)]

                if target_final:
                    result[source] = target_final

    return result


def build_meta_mapping(meta: Path | None, split_ids: set[str], final_ids: set[str]) -> dict[str, str]:
    if meta is None or not meta.is_file():
        return {}

    df = read_any_csv(meta)
    final_canon = {canonical_id(x): x for x in final_ids}
    split_canon = {canonical_id(x): x for x in split_ids}

    result: dict[str, str] = {}

    for source_col in df.columns:
        for target_col in df.columns:
            if source_col == target_col:
                continue

            for _, row in df.iterrows():
                source_raw = str(row.get(source_col, "")).strip()
                target_raw = str(row.get(target_col, "")).strip()

                source_can = canonical_id(source_raw)
                target_can = canonical_id(target_raw)

                if source_can in split_canon and target_can in final_canon:
                    result[split_canon[source_can]] = final_canon[target_can]

                if target_can in split_canon and source_can in final_canon:
                    result[split_canon[target_can]] = final_canon[source_can]

    return result


def infer_orientation(header: list[str], first_row: list[str]) -> str:
    header_cg_count = sum(1 for x in header[1:min(len(header), 2000)] if normalize_cpg(x))
    first_col_is_cg = normalize_cpg(first_row[0]) is not None if first_row else False

    if header_cg_count > 10:
        return "sample_rows"
    if first_col_is_cg:
        return "cpg_rows"

    raise ValueError(
        "Could not infer final matrix orientation. "
        f"Header: {header[:5]}, first row: {first_row[:5]}"
    )


def is_missing(x: str) -> bool:
    return str(x).strip() in NA_VALUES


def get_final_matrix_samples_and_cpgs(final_beta: Path, sep: str = "\t") -> tuple[str, list[str], set[str]]:
    with open_text(final_beta, "rt") as fh:
        reader = csv.reader(fh, delimiter=sep)
        header = next(reader)
        first_row = next(reader)

        orientation = infer_orientation(header, first_row)

        if orientation == "sample_rows":
            sample_ids = [first_row[0]]
            cpgs = {normalize_cpg(x) for x in header[1:]}
            cpgs = {x for x in cpgs if x}

            for row in reader:
                if row:
                    sample_ids.append(row[0])

            return orientation, sample_ids, cpgs

        sample_ids = header[1:]
        cpgs = set()
        cg = normalize_cpg(first_row[0])
        if cg:
            cpgs.add(cg)

        for row in reader:
            if row:
                cg = normalize_cpg(row[0])
                if cg:
                    cpgs.add(cg)

        return orientation, sample_ids, cpgs


def resolve_split_to_final_ids(
    split_ids: set[str],
    final_sample_ids: list[str],
    meta: Path | None,
    sample_map: Path | None,
) -> tuple[dict[str, str], pd.DataFrame]:
    final_set = set(final_sample_ids)
    final_canon = {canonical_id(x): x for x in final_sample_ids}

    explicit_map = load_optional_map(sample_map, final_set)
    meta_map = build_meta_mapping(meta, split_ids, final_set)

    rows = []
    resolved: dict[str, str] = {}

    for sid in sorted(split_ids):
        final_id = None
        method = None

        if sid in final_set:
            final_id = sid
            method = "exact_final_id"
        elif canonical_id(sid) in final_canon:
            final_id = final_canon[canonical_id(sid)]
            method = "canonical_final_id"
        elif sid in explicit_map:
            final_id = explicit_map[sid]
            method = "explicit_sample_map"
        elif sid in meta_map:
            final_id = meta_map[sid]
            method = "meta_map"
        else:
            sid_can = canonical_id(sid)
            for k, v in explicit_map.items():
                if canonical_id(k) == sid_can:
                    final_id = v
                    method = "explicit_sample_map_canonical"
                    break

        matched = final_id is not None

        if matched:
            resolved[sid] = final_id

        rows.append({
            "split_sample_id": sid,
            "final_sample_id": final_id if final_id else "",
            "matched": matched,
            "method": method if method else "unmatched",
        })

    return resolved, pd.DataFrame(rows)


def read_final_matrix_split_presence(
    final_beta: Path,
    requested: set[str],
    split_to_final_ids: dict[str, set[str]],
    sep: str = "\t",
) -> dict:
    counters = {split: Counter() for split in split_to_final_ids}
    observed_samples = {split: set() for split in split_to_final_ids}
    orientation = None

    with open_text(final_beta, "rt") as fh:
        reader = csv.reader(fh, delimiter=sep)
        header = next(reader)
        first_row = next(reader)

        orientation = infer_orientation(header, first_row)

        if orientation == "sample_rows":
            selected_cols = {}

            for idx, col in enumerate(header):
                cg = normalize_cpg(col)
                if cg and cg in requested:
                    selected_cols[idx] = cg

            def process_row(row: list[str]):
                if not row:
                    return
                sample_id = row[0]

                for split, final_ids in split_to_final_ids.items():
                    if sample_id not in final_ids:
                        continue

                    observed_samples[split].add(sample_id)

                    for idx, cg in selected_cols.items():
                        if idx < len(row) and not is_missing(row[idx]):
                            counters[split][cg] += 1

            process_row(first_row)

            for row in reader:
                process_row(row)

        else:
            sample_cols = header[1:]
            split_col_indices = {
                split: [
                    i + 1
                    for i, sample in enumerate(sample_cols)
                    if sample in final_ids
                ]
                for split, final_ids in split_to_final_ids.items()
            }

            for split, idxs in split_col_indices.items():
                observed_samples[split] = {
                    sample_cols[i - 1] for i in idxs
                }

            def process_cpg_row(row: list[str]):
                if not row:
                    return
                cg = normalize_cpg(row[0])
                if cg not in requested:
                    return

                for split, idxs in split_col_indices.items():
                    for idx in idxs:
                        if idx < len(row) and not is_missing(row[idx]):
                            counters[split][cg] += 1

            process_cpg_row(first_row)

            for row in reader:
                process_cpg_row(row)

    return {
        "orientation": orientation,
        "counters": counters,
        "observed_samples": observed_samples,
    }


def read_qc_affected_by_split(
    batch_root: Path,
    raw_to_final: dict[str, str],
    split_final_ids: dict[str, set[str]],
) -> dict[str, set[str]]:
    affected = {split: set() for split in split_final_ids}

    for batch_dir in sorted(batch_root.glob("batch_*")):
        qc_path = batch_dir / "qc_probe_failures_by_sample.tsv.gz"

        if not qc_path.is_file():
            qc_path = batch_dir / "qc_probe_failures.tsv.gz"

        if not qc_path.is_file():
            continue

        with open_text(qc_path, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")

            if reader.fieldnames is None:
                continue

            sample_col = None
            for c in ["SampleID", "Sample_ID", "sample_id", "sample"]:
                if c in reader.fieldnames:
                    sample_col = c
                    break

            probe_col = None
            for c in ["ProbeID", "probe_id", "CpG", "cpg_id"]:
                if c in reader.fieldnames:
                    probe_col = c
                    break

            if probe_col is None:
                continue

            for row in reader:
                probe = normalize_cpg(row.get(probe_col))
                if not probe:
                    continue

                final_id = None

                if sample_col is not None:
                    raw_sample = str(row.get(sample_col, "")).strip()
                    final_id = raw_to_final.get(raw_sample)

                    if final_id is None:
                        final_id = raw_sample

                for split, final_ids in split_final_ids.items():
                    if final_id in final_ids:
                        affected[split].add(probe)

    return affected


def write_list(path: Path, values: set[str]) -> None:
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
    parser.add_argument("--sample-sheet", required=True, type=Path)
    parser.add_argument("--training-samples", required=True, type=Path)
    parser.add_argument("--testing-samples", required=True, type=Path)
    parser.add_argument("--meta", type=Path, default=None)
    parser.add_argument("--sample-map", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--sep", default="\t")

    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    probe_sets = load_dunedinpace_sets()
    requested = probe_sets["union_requested"]

    print(f"[SplitCoverage] dataset: {args.dataset}")
    print(f"[SplitCoverage] requested union: {len(requested):,}")

    print("[SplitCoverage] scanning manifest...")
    manifest_present = scan_manifest(args.manifest, requested)

    orientation, final_sample_ids, final_matrix_cpgs = get_final_matrix_samples_and_cpgs(
        args.final_beta,
        sep=args.sep,
    )

    print(f"[SplitCoverage] final matrix orientation: {orientation}")
    print(f"[SplitCoverage] final samples: {len(final_sample_ids):,}")
    print(f"[SplitCoverage] final CpG features: {len(final_matrix_cpgs):,}")

    training_ids = load_split_ids(args.training_samples)
    testing_ids = load_split_ids(args.testing_samples)

    train_map, train_map_df = resolve_split_to_final_ids(
        training_ids,
        final_sample_ids,
        args.meta,
        args.sample_map,
    )

    test_map, test_map_df = resolve_split_to_final_ids(
        testing_ids,
        final_sample_ids,
        args.meta,
        args.sample_map,
    )

    train_map_df.insert(0, "dataset", args.dataset)
    train_map_df.insert(1, "split", "training")
    test_map_df.insert(0, "dataset", args.dataset)
    test_map_df.insert(1, "split", "testing")

    sample_map_df = pd.concat([train_map_df, test_map_df], ignore_index=True)
    sample_map_df.to_csv(
        args.out_dir / f"{args.dataset}.split_sample_mapping.csv",
        index=False,
    )

    split_final_ids = {
        "training": set(train_map.values()),
        "testing": set(test_map.values()),
    }

    print(
        f"[SplitCoverage] matched training samples: "
        f"{len(split_final_ids['training'])}/{len(training_ids)}"
    )
    print(
        f"[SplitCoverage] matched testing samples: "
        f"{len(split_final_ids['testing'])}/{len(testing_ids)}"
    )

    final_presence = read_final_matrix_split_presence(
        args.final_beta,
        requested,
        split_final_ids,
        sep=args.sep,
    )

    raw_to_final = load_raw_to_final_sample_map(args.sample_sheet)
    qc_affected = read_qc_affected_by_split(
        args.batch_root,
        raw_to_final,
        split_final_ids,
    )

    rows = []

    for split in ["training", "testing"]:
        n_split_samples = len(split_final_ids[split])
        counter = final_presence["counters"][split]
        observed_final_samples = final_presence["observed_samples"][split]

        for set_name, req in probe_sets.items():
            final_any = {cg for cg in req if counter.get(cg, 0) > 0}
            final_all = {
                cg for cg in req
                if n_split_samples > 0 and counter.get(cg, 0) == n_split_samples
            }
            affected = req & qc_affected[split]

            row = {
                "dataset": args.dataset,
                "split": split,
                "probe_set": set_name,
                "requested_n": len(req),
                "manifest_present_n": len(req & manifest_present),
                "manifest_missing_n": len(req - manifest_present),
                "final_matrix_feature_present_n": len(req & final_matrix_cpgs),
                "final_value_present_any_split_sample_n": len(final_any),
                "final_value_present_all_split_samples_n": len(final_all),
                "final_value_missing_all_split_samples_n": len(req - final_any),
                "final_value_missing_at_least_one_split_sample_n": len(req - final_all),
                "qc_affected_at_least_one_split_sample_n": len(affected),
                "qc_affected_pct_of_requested": (
                    100 * len(affected) / len(req) if len(req) else 0
                ),
                "expected_split_sample_n": len(training_ids) if split == "training" else len(testing_ids),
                "matched_final_sample_n": n_split_samples,
                "observed_final_sample_n": len(observed_final_samples),
                "final_beta_orientation": final_presence["orientation"],
            }

            rows.append(row)

            prefix = f"{args.dataset}.{split}.{set_name}"
            write_list(args.out_dir / f"{prefix}.present_in_manifest.txt.gz", req & manifest_present)
            write_list(args.out_dir / f"{prefix}.present_final_any_split_sample.txt.gz", final_any)
            write_list(args.out_dir / f"{prefix}.present_final_all_split_samples.txt.gz", final_all)
            write_list(args.out_dir / f"{prefix}.qc_affected.txt.gz", affected)

    summary = pd.DataFrame(rows)

    csv_path = args.out_dir / f"{args.dataset}.DunedinPACE20k_split_coverage.csv"
    summary.to_csv(csv_path, index=False)

    metadata = {
        "dataset": args.dataset,
        "manifest": str(args.manifest),
        "batch_root": str(args.batch_root),
        "final_beta": str(args.final_beta),
        "training_samples": str(args.training_samples),
        "testing_samples": str(args.testing_samples),
        "sample_sheet": str(args.sample_sheet),
        "meta": str(args.meta) if args.meta else None,
        "sample_map": str(args.sample_map) if args.sample_map else None,
    }

    with open(args.out_dir / f"{args.dataset}.DunedinPACE20k_split_coverage.metadata.json", "w") as fh:
        json.dump(metadata, fh, indent=2)

    print("[SplitCoverage] wrote:", csv_path)


if __name__ == "__main__":
    main()
