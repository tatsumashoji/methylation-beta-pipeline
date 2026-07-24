#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

CLOCK_MODELS = {
    "Horvath": "Horvathv1",
    "Hannum": "Hannum",
    "PhenoAge": "PhenoAge",
    "GrimAgeV1": "GrimAgeV1",
    "GrimAgeV2": "GrimAgeV2",
    "DunedinPACE": "DunedinPACE",
}

SPLIT_ID_COLUMNS = [
    "Name",
    "Name_Original",
    "SampleID",
    "Sample_ID",
    "sample_id",
    "Sample",
    "ID",
    "id",
    "Unnamed: 0",
]

RAW_BASENAME_RE = re.compile(r"([0-9]{8,}_R[0-9]{2}C[0-9]{2})")
ID_LIKE_RE = re.compile(r"^(YS_|KAZU_|[A-Za-z]*[0-9]{8,}_R[0-9]{2}C[0-9]{2}|[0-9]{8,}_R[0-9]{2}C[0-9]{2})")
CPG_RE = re.compile(r"^(cg|ch|rs|ctl|CTL|SNP)")


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def sniff_delimiter(path: Path) -> str:
    with open_text(path, "rt") as fh:
        sample = fh.read(8192)
    if not sample:
        return ","
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        return dialect.delimiter
    except Exception:
        # Most pipeline files are tab-delimited. Most templates are comma-delimited.
        first = sample.splitlines()[0] if sample.splitlines() else ""
        return "\t" if "\t" in first and "," not in first else ","


def read_table(path: Path) -> Tuple[List[str], List[dict], str]:
    delimiter = sniff_delimiter(path)
    with open_text(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            return [], [], delimiter
        return list(reader.fieldnames), [dict(r) for r in reader], delimiter


def first_existing_column(fieldnames: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    lookup = {str(f).strip().lower(): str(f) for f in fieldnames}
    for c in candidates:
        if str(c).lower() in lookup:
            return lookup[str(c).lower()]
    return None


def clean_value(value: object) -> str:
    if value is None:
        return ""
    s = str(value).strip()
    if s.lower() in {"", "na", "nan", "none", "null"}:
        return ""
    return s


def raw_basename(value: str) -> str:
    value = clean_value(value)
    if not value:
        return ""
    name = Path(value).name
    name = re.sub(r"_(Grn|Red)\.idat(\.gz)?$", "", name)
    m = RAW_BASENAME_RE.search(name)
    return m.group(1) if m else name


def is_id_like(value: str, known_ids: Set[str]) -> bool:
    value = clean_value(value)
    if not value:
        return False
    if value in known_ids:
        return True
    if RAW_BASENAME_RE.search(value):
        return True
    if ID_LIKE_RE.search(value):
        return True
    # KAZU_01 etc and YS_002 are already covered; keep this conservative.
    return False


class UnionFind:
    def __init__(self) -> None:
        self.parent: Dict[str, str] = {}

    def add(self, x: str) -> None:
        x = clean_value(x)
        if x and x not in self.parent:
            self.parent[x] = x

    def find(self, x: str) -> str:
        self.add(x)
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        a = clean_value(a)
        b = clean_value(b)
        if not a or not b:
            return
        ra = self.find(a)
        rb = self.find(b)
        if ra != rb:
            # Stable-ish parent choice.
            self.parent[rb] = ra

    def groups(self) -> Dict[str, Set[str]]:
        out: Dict[str, Set[str]] = defaultdict(set)
        for x in list(self.parent):
            out[self.find(x)].add(x)
        return dict(out)

    def aliases(self, x: str) -> Set[str]:
        x = clean_value(x)
        if not x:
            return set()
        if x not in self.parent:
            return {x}
        return self.groups().get(self.find(x), {x})


def read_split_ids(path: Path, explicit_column: Optional[str] = None) -> List[str]:
    fieldnames, rows, _ = read_table(path)
    if not fieldnames:
        ids: List[str] = []
        with open_text(path, "rt") as fh:
            for line in fh:
                x = clean_value(line.strip().split("\t")[0].split(",")[0])
                if x and x.lower() not in {"sampleid", "sample_id", "sample", "id", "name"}:
                    ids.append(x)
        return ids

    if explicit_column:
        if explicit_column not in fieldnames:
            raise ValueError(f"Column {explicit_column!r} not found in {path}. Columns: {fieldnames}")
        col = explicit_column
    else:
        col = first_existing_column(fieldnames, SPLIT_ID_COLUMNS) or fieldnames[0]

    ids = []
    for row in rows:
        v = clean_value(row.get(col))
        if v:
            ids.append(v)
    return ids


def read_probe_list(path: Path) -> Set[str]:
    probes: Set[str] = set()
    with open_text(path, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            v = clean_value(line.strip().split("\t")[0].split(",")[0])
            if v and v.lower() not in {"probeid", "probe_id", "cpg", "cpg_id", "sample_id", "sampleid"}:
                probes.add(v)
    return probes


def read_final_beta_ids(path: Path, sep: str = "\t") -> Tuple[Set[str], Set[str], str]:
    """Return (sample_ids, probe_ids, orientation).

    orientation is sample_rows if rows are samples and columns are CpGs;
    cpg_rows if rows are CpGs and columns are samples.
    """
    with open_text(path, "rt") as fh:
        header_line = fh.readline().rstrip("\n\r")
        if not header_line:
            raise ValueError(f"Empty final beta matrix: {path}")
        header = header_line.split(sep)
        second_line = fh.readline().rstrip("\n\r")
        second = second_line.split(sep) if second_line else []

    if len(header) < 2:
        raise ValueError(f"Final beta matrix has too few columns: {path}")

    header_cpg_count = sum(1 for x in header[1: min(len(header), 2000)] if CPG_RE.search(x))
    first_data_value = second[0] if second else ""
    first_data_is_cpg = bool(CPG_RE.search(first_data_value))

    if header_cpg_count >= max(1, min(100, len(header[1: min(len(header), 2000)]) // 2)) and not first_data_is_cpg:
        orientation = "sample_rows"
        probe_ids = set(header[1:])
        sample_ids: Set[str] = set()
        with open_text(path, "rt") as fh:
            fh.readline()
            for line in fh:
                if not line.strip():
                    continue
                sample_ids.add(clean_value(line.split(sep, 1)[0]))
        sample_ids.discard("")
        return sample_ids, probe_ids, orientation

    orientation = "cpg_rows"
    sample_ids = {clean_value(x) for x in header[1:] if clean_value(x)}
    probe_ids: Set[str] = set()
    with open_text(path, "rt") as fh:
        fh.readline()
        for line in fh:
            if not line.strip():
                continue
            probe_ids.add(clean_value(line.split(sep, 1)[0]))
    probe_ids.discard("")
    return sample_ids, probe_ids, orientation


def add_sample_sheet_aliases(path: Path, uf: UnionFind, known_ids: Set[str]) -> Set[str]:
    processed_raw_ids: Set[str] = set()
    if not path.is_file():
        return processed_raw_ids
    fieldnames, rows, _ = read_table(path)
    for row in rows:
        ids: Set[str] = set()
        for col in fieldnames:
            value = clean_value(row.get(col))
            if not value:
                continue
            if col.lower() in {"grn", "red", "green", "grn_idat", "red_idat"} or value.endswith(".idat") or ".idat" in value:
                rb = raw_basename(value)
                if rb:
                    ids.add(rb)
                    processed_raw_ids.add(rb)
            elif is_id_like(value, known_ids):
                ids.add(value)
        ids = {x for x in ids if x}
        for x in ids:
            uf.add(x)
        ids_list = sorted(ids)
        for x in ids_list[1:]:
            uf.union(ids_list[0], x)
    return processed_raw_ids


def add_meta_aliases(path: Optional[Path], uf: UnionFind, known_ids: Set[str]) -> None:
    if path is None or not path.is_file():
        return
    fieldnames, rows, _ = read_table(path)
    for row in rows:
        ids: Set[str] = set()
        for col in fieldnames:
            value = clean_value(row.get(col))
            if not value:
                continue
            rb = raw_basename(value) if (".idat" in value or RAW_BASENAME_RE.search(value)) else ""
            if rb and is_id_like(rb, known_ids):
                ids.add(rb)
            if is_id_like(value, known_ids):
                ids.add(value)
        ids = {x for x in ids if x}
        if len(ids) >= 2:
            ids_list = sorted(ids)
            for x in ids_list:
                uf.add(x)
            for x in ids_list[1:]:
                uf.union(ids_list[0], x)


def add_sample_map_aliases(path: Optional[Path], uf: UnionFind, known_ids: Set[str]) -> None:
    if path is None or not path.is_file():
        return
    fieldnames, rows, _ = read_table(path)
    if not fieldnames:
        return
    for row in rows:
        ids = {clean_value(row.get(c)) for c in fieldnames}
        ids = {raw_basename(x) if ".idat" in x else x for x in ids if x}
        ids = {x for x in ids if is_id_like(x, known_ids) or x in known_ids}
        if len(ids) >= 2:
            ids_list = sorted(ids)
            for x in ids_list:
                uf.add(x)
            for x in ids_list[1:]:
                uf.union(ids_list[0], x)


def load_qc_processed_sample_ids(batch_root: Path) -> Set[str]:
    ids: Set[str] = set()
    # Prefer qc_summary.txt if available.
    for path in sorted(batch_root.glob("batch_*/qc_summary.txt")):
        try:
            fieldnames, rows, _ = read_table(path)
        except Exception:
            continue
        if not fieldnames or not rows:
            continue
        col = first_existing_column(fieldnames, ["SampleID", "SampleName", "sample_id", "Sample"])
        if col is None:
            continue
        for row in rows:
            v = clean_value(row.get(col))
            if v:
                ids.add(v)
    # Also include raw IDs that appear in failure files.
    for path in sorted(batch_root.glob("batch_*/qc_probe_failures_by_sample.tsv.gz")):
        with open_text(path, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames:
                continue
            sample_col = first_existing_column(reader.fieldnames, ["SampleID", "sample_id", "Sample", "SampleName"])
            if sample_col is None:
                continue
            for row in reader:
                v = clean_value(row.get(sample_col))
                if v:
                    ids.add(v)
    return ids


def load_clock_sites() -> Dict[str, Set[str]]:
    from biolearn.model_gallery import ModelGallery
    gallery = ModelGallery()
    out: Dict[str, Set[str]] = {}
    for label, model_name in CLOCK_MODELS.items():
        model = gallery.get(model_name, imputation_method="none")
        if not hasattr(model, "methylation_sites"):
            raise RuntimeError(f"{model_name} does not expose methylation_sites()")
        out[label] = {str(x).strip() for x in model.methylation_sites() if str(x).strip()}
    return out


def write_csv(path: Path, rows: Iterable[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_tsv_gz(path: Path, rows: Iterable[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--batch-root", required=True, type=Path)
    ap.add_argument("--final-beta", required=True, type=Path)
    ap.add_argument("--sample-sheet", required=True, type=Path)
    ap.add_argument("--training-samples", required=True, type=Path)
    ap.add_argument("--validation-samples", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--meta", type=Path, default=None)
    ap.add_argument("--sample-map", type=Path, default=None,
                    help="Optional explicit alias table connecting training/testing IDs to final or raw IDs.")
    ap.add_argument("--training-id-column", default=None)
    ap.add_argument("--validation-id-column", default=None)
    ap.add_argument("--sep", default="\t")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    training_ids = set(read_split_ids(args.training_samples, args.training_id_column))
    validation_ids = set(read_split_ids(args.validation_samples, args.validation_id_column))
    split_defs = {"training": training_ids, "validation": validation_ids}

    final_sample_ids, final_probe_ids, final_orientation = read_final_beta_ids(args.final_beta, args.sep)
    processed_sample_ids = load_qc_processed_sample_ids(args.batch_root)

    known_ids = set(training_ids) | set(validation_ids) | set(final_sample_ids) | set(processed_sample_ids)

    uf = UnionFind()
    for x in known_ids:
        uf.add(x)

    sample_sheet_raw_ids = add_sample_sheet_aliases(args.sample_sheet, uf, known_ids)
    known_ids |= sample_sheet_raw_ids
    for x in sample_sheet_raw_ids:
        uf.add(x)

    # Re-run with expanded known IDs to capture row aliases in name_table.
    add_sample_sheet_aliases(args.sample_sheet, uf, known_ids)
    add_meta_aliases(args.meta, uf, known_ids)
    add_sample_map_aliases(args.sample_map, uf, known_ids)

    processed_raw_ids = processed_sample_ids | sample_sheet_raw_ids

    # Map processed raw IDs to splits through alias groups.
    split_processed_ids: Dict[str, Set[str]] = {"training": set(), "validation": set()}
    split_alias_matches: Dict[str, Dict[str, Set[str]]] = {"training": defaultdict(set), "validation": defaultdict(set)}

    for raw_id in processed_raw_ids:
        aliases = uf.aliases(raw_id)
        for split, ids in split_defs.items():
            hits = aliases & ids
            if hits:
                split_processed_ids[split].add(raw_id)
                for h in hits:
                    split_alias_matches[split][h].add(raw_id)

    # In case split IDs directly appear as processed IDs.
    for split, ids in split_defs.items():
        for sid in ids:
            if sid in processed_raw_ids:
                split_processed_ids[split].add(sid)
                split_alias_matches[split][sid].add(sid)

    pre_qc_union: Set[str] = set()
    # Do not use post_qc names as proof of availability; names can remain even when values are NA.
    failed_sample_sets: Dict[str, Dict[str, Set[str]]] = {
        "training": defaultdict(set),
        "validation": defaultdict(set),
    }
    affected_batch_count: Dict[str, Counter] = {"training": Counter(), "validation": Counter()}

    for batch_dir in sorted(args.batch_root.glob("batch_*")):
        pre_path = batch_dir / "pre_qc_collapsed_probe_ids.txt.gz"
        q_path = batch_dir / "qc_probe_failures_by_sample.tsv.gz"
        if not pre_path.is_file():
            raise FileNotFoundError(pre_path)
        if not q_path.is_file():
            raise FileNotFoundError(q_path)

        pre_qc_union.update(read_probe_list(pre_path))
        touched_in_batch: Dict[str, Set[str]] = {"training": set(), "validation": set()}

        with open_text(q_path, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None:
                continue
            sample_col = first_existing_column(reader.fieldnames, ["SampleID", "sample_id", "Sample", "SampleName"])
            probe_col = first_existing_column(reader.fieldnames, ["ProbeID", "probe_id", "CpG", "cpg_id"])
            if sample_col is None or probe_col is None:
                raise ValueError(f"{q_path} must contain SampleID and ProbeID columns. Columns: {reader.fieldnames}")

            for row in reader:
                raw_id = clean_value(row.get(sample_col))
                probe = clean_value(row.get(probe_col))
                if not raw_id or not probe:
                    continue
                aliases = uf.aliases(raw_id)
                for split, ids in split_defs.items():
                    if aliases & ids or raw_id in split_processed_ids[split]:
                        failed_sample_sets[split][probe].add(raw_id)
                        touched_in_batch[split].add(probe)

        for split in split_defs:
            for probe in touched_in_batch[split]:
                affected_batch_count[split][probe] += 1

    clock_sites = load_clock_sites()
    clock_names = list(CLOCK_MODELS.keys())

    def required_flags(probe: str) -> dict:
        clocks = []
        row = {}
        for c in clock_names:
            hit = probe in clock_sites[c]
            row[f"required_by_{c}"] = hit
            if hit:
                clocks.append(c)
        row["required_clock_count"] = len(clocks)
        row["required_clocks"] = ";".join(clocks)
        return row

    detail_fields = [
        "dataset", "split", "ProbeID", "qc_failed_samples", "n_samples_in_split",
        "qc_failed_fraction", "affected_batch_count", "qc_completely_lost_in_split",
        "qc_partially_affected_in_split", "present_before_qc", "present_in_integrated_matrix",
    ] + [f"required_by_{c}" for c in clock_names] + ["required_clock_count", "required_clocks"]

    affected_rows: List[dict] = []
    summary_rows: List[dict] = []
    clock_summary_rows: List[dict] = []
    clock_detail_rows: List[dict] = []

    for split in ["training", "validation"]:
        expected = len(split_defs[split])
        matched = len(split_processed_ids[split])
        affected = set(failed_sample_sets[split].keys())
        completely_lost = {p for p in affected if matched > 0 and len(failed_sample_sets[split][p]) >= matched}
        partially = affected - completely_lost
        post_qc_cpg_n_estimated = max(0, len(pre_qc_union) - len(completely_lost))

        for probe in sorted(affected):
            failed_n = len(failed_sample_sets[split][probe])
            row = {
                "dataset": args.dataset,
                "split": split,
                "ProbeID": probe,
                "qc_failed_samples": failed_n,
                "n_samples_in_split": matched,
                "qc_failed_fraction": failed_n / matched if matched else 0,
                "affected_batch_count": affected_batch_count[split][probe],
                "qc_completely_lost_in_split": probe in completely_lost,
                "qc_partially_affected_in_split": probe in partially,
                "present_before_qc": probe in pre_qc_union,
                "present_in_integrated_matrix": probe in final_probe_ids,
            }
            row.update(required_flags(probe))
            affected_rows.append(row)

        summary_rows.append({
            "dataset": args.dataset,
            "split": split,
            "expected_sample_n_from_split_file": expected,
            "matched_sample_n_in_qc_outputs": matched,
            "qc_affected_cpg_n": len(affected),
            "qc_completely_lost_cpg_n": len(completely_lost),
            "qc_partially_affected_cpg_n": len(partially),
            "pre_qc_cpg_n": len(pre_qc_union),
            "post_qc_cpg_n_estimated_for_split": post_qc_cpg_n_estimated,
            "integrated_matrix_cpg_n": len(final_probe_ids),
            "final_beta_orientation": final_orientation,
        })

        for clock, sites in clock_sites.items():
            req = len(sites)
            affected_clock = sites & affected
            lost_clock = sites & completely_lost
            partial_clock = sites & partially
            clock_summary_rows.append({
                "dataset": args.dataset,
                "split": split,
                "clock": clock,
                "biolearn_model": CLOCK_MODELS[clock],
                "required_cpg_n": req,
                "qc_affected_cpg_n": len(affected_clock),
                "qc_affected_pct_of_required": 100 * len(affected_clock) / req if req else 0,
                "qc_completely_lost_cpg_n": len(lost_clock),
                "qc_completely_lost_pct_of_required": 100 * len(lost_clock) / req if req else 0,
                "qc_partially_affected_cpg_n": len(partial_clock),
                "present_before_qc_n": len(sites & pre_qc_union),
                "present_in_integrated_matrix_n": len(sites & final_probe_ids),
            })
            for probe in sorted(affected_clock):
                failed_n = len(failed_sample_sets[split][probe])
                clock_detail_rows.append({
                    "dataset": args.dataset,
                    "split": split,
                    "clock": clock,
                    "biolearn_model": CLOCK_MODELS[clock],
                    "ProbeID": probe,
                    "qc_failed_samples": failed_n,
                    "n_samples_in_split": matched,
                    "qc_failed_fraction": failed_n / matched if matched else 0,
                    "qc_completely_lost_in_split": probe in completely_lost,
                    "present_in_integrated_matrix": probe in final_probe_ids,
                })

    write_csv(args.out_dir / "split_qc_probe_loss_summary.csv", summary_rows, list(summary_rows[0].keys()) if summary_rows else [])
    write_tsv_gz(args.out_dir / "split_qc_affected_cpgs.tsv.gz", affected_rows, detail_fields)
    write_csv(args.out_dir / "split_clock_qc_overlap_summary.csv", clock_summary_rows, list(clock_summary_rows[0].keys()) if clock_summary_rows else [])
    write_tsv_gz(
        args.out_dir / "split_clock_qc_overlap_cpgs.tsv.gz",
        clock_detail_rows,
        list(clock_detail_rows[0].keys()) if clock_detail_rows else [
            "dataset", "split", "clock", "biolearn_model", "ProbeID", "qc_failed_samples",
            "n_samples_in_split", "qc_failed_fraction", "qc_completely_lost_in_split", "present_in_integrated_matrix",
        ],
    )

    unmatched_rows: List[dict] = []
    for split, ids in split_defs.items():
        matched_split_ids = set(split_alias_matches[split].keys())
        for sid in sorted(ids - matched_split_ids):
            unmatched_rows.append({"dataset": args.dataset, "split": split, "unmatched_sample_id": sid})
    write_csv(args.out_dir / "split_unmatched_samples.csv", unmatched_rows, ["dataset", "split", "unmatched_sample_id"])

    # Mapping diagnostics: useful for checking YS/KAZU/IDAT mapping.
    mapping_rows: List[dict] = []
    for split, matches in split_alias_matches.items():
        for split_id, raws in sorted(matches.items()):
            for raw in sorted(raws):
                mapping_rows.append({
                    "dataset": args.dataset,
                    "split": split,
                    "split_sample_id": split_id,
                    "qc_raw_sample_id": raw,
                    "aliases": ";".join(sorted(uf.aliases(raw))),
                })
    write_csv(args.out_dir / "split_sample_id_mapping.csv", mapping_rows, ["dataset", "split", "split_sample_id", "qc_raw_sample_id", "aliases"])

    metadata = {
        "dataset": args.dataset,
        "batch_root": str(args.batch_root),
        "final_beta": str(args.final_beta),
        "final_beta_orientation": final_orientation,
        "final_sample_n": len(final_sample_ids),
        "final_cpg_n": len(final_probe_ids),
        "sample_sheet": str(args.sample_sheet),
        "meta": str(args.meta) if args.meta else None,
        "sample_map": str(args.sample_map) if args.sample_map else None,
        "training_samples": str(args.training_samples),
        "validation_samples": str(args.validation_samples),
        "training_expected_n": len(training_ids),
        "validation_expected_n": len(validation_ids),
        "training_matched_n": len(split_processed_ids["training"]),
        "validation_matched_n": len(split_processed_ids["validation"]),
        "notes": [
            "Final beta orientation is auto-detected, so sample_rows matrices are handled correctly.",
            "QC affected CpGs are counted from qc_probe_failures_by_sample.tsv.gz, i.e. sample-level pre_qc_beta non-NA and post_qc_beta NA records.",
            "post_qc_cpg_n_estimated_for_split = pre_qc_cpg_n - qc_completely_lost_cpg_n.",
        ],
    }
    with open(args.out_dir / "split_qc_report_metadata.json", "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2, ensure_ascii=False)

    print(f"[Done] {args.dataset}: {args.out_dir}")
    print(f"[Done] final beta orientation: {final_orientation}")
    print(f"[Done] final CpGs: {len(final_probe_ids):,}; final samples: {len(final_sample_ids):,}")
    for row in summary_rows:
        print(
            f"[Done] {row['dataset']} {row['split']}: "
            f"matched n={row['matched_sample_n_in_qc_outputs']}, "
            f"affected CpGs={row['qc_affected_cpg_n']}, "
            f"completely lost CpGs={row['qc_completely_lost_cpg_n']}"
        )


if __name__ == "__main__":
    main()
