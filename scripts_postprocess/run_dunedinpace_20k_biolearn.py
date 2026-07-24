#!/usr/bin/env python3
"""Recompute DunedinPACE using measured 20,000 gold/background probes.

This memory-efficient version is designed for very wide collapsed beta matrices
such as Sample_ID x ~900k CpG columns.  It avoids pandas.read_csv() on the full
wide matrix and instead streams the file line-by-line, extracting only the
DunedinPACE 20k gold/background probes plus the 173 core CpGs.

Supported input matrix orientations:
  1. sample rows x CpG columns: Sample_ID, cg000..., cg000...
  2. CpG rows x sample columns: ProbeID, sample1, sample2, ...
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from biolearn.data_library import GeoData
from biolearn.model_gallery import ModelGallery
from biolearn.util import get_data_file

CPG_PREFIXES = ("cg", "ch")


def is_cpg_like(value: object) -> bool:
    s = str(value).strip()
    return s.startswith(CPG_PREFIXES)


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def infer_delimiter(path: Path) -> str:
    with open_text(path, "rt") as fh:
        first = fh.readline()
    return "\t" if "\t" in first else ","


def read_header_and_second_line(path: Path, sep: str) -> tuple[list[str], list[str]]:
    with open_text(path, "rt") as fh:
        header = fh.readline().rstrip("\n\r").split(sep)
        second = fh.readline().rstrip("\n\r").split(sep)
    return header, second


def infer_matrix_orientation(path: Path, sep: str) -> tuple[str, str]:
    header, second = read_header_and_second_line(path, sep)
    first_col = header[0]
    header_cpg_count = sum(is_cpg_like(x) for x in header[1:min(len(header), 2000)])
    first_value = second[0] if second else ""

    if header_cpg_count >= 5:
        return "sample_rows", first_col
    if is_cpg_like(first_value) or first_col.lower() in {"probeid", "probe_id", "cpg", "cpg_id"}:
        return "cpg_rows", first_col

    raise RuntimeError(
        f"Could not infer matrix orientation for {path}. "
        f"first_col={first_col!r}, first_value={first_value!r}, "
        f"cpg-like header count={header_cpg_count}."
    )


def get_dunedinpace_probe_sets() -> tuple[set[str], set[str], set[str]]:
    gallery = ModelGallery()
    model = gallery.get("DunedinPACE")
    core = {str(x).strip() for x in model.methylation_sites() if str(x).strip()}

    gold = pd.read_csv(get_data_file("DunedinPACE_Gold_Means.csv"), index_col=0)
    gold_set = {str(x).strip() for x in gold.index if str(x).strip()}

    return core, gold_set, core | gold_set


def parse_float(value: str) -> float:
    value = value.strip()
    if value == "" or value.upper() == "NA" or value.lower() == "nan":
        return np.nan
    try:
        return float(value)
    except Exception:
        return np.nan


def load_sample_rows_subset_streaming(
    path: Path,
    sep: str,
    wanted_probes: set[str],
    sample_id_col: str,
    progress_every: int,
) -> tuple[pd.DataFrame, dict]:
    """Read Sample_ID x CpG-column matrix without pandas wide parsing."""
    with open_text(path, "rt") as fh:
        header_line = fh.readline()
        if not header_line:
            raise RuntimeError(f"Empty input file: {path}")
        header = header_line.rstrip("\n\r").split(sep)

        if sample_id_col not in header:
            sample_idx = 0
            sample_id_col = header[0]
        else:
            sample_idx = header.index(sample_id_col)

        # Identify requested columns from the header.  Keep input order to avoid
        # unnecessary sorting overhead during parsing.
        selected_indices: list[int] = []
        selected_probes: list[str] = []
        matrix_cpg_n = 0
        for i, col in enumerate(header):
            if i == sample_idx:
                continue
            if is_cpg_like(col):
                matrix_cpg_n += 1
            if col in wanted_probes:
                selected_indices.append(i)
                selected_probes.append(col)

        if not selected_indices:
            raise RuntimeError("None of the requested DunedinPACE probes are present in the matrix header.")

        print(
            f"[DunedinPACE20k] sample-row matrix columns: {len(header):,}; "
            f"CpG-like columns: {matrix_cpg_n:,}; selected requested probes: {len(selected_probes):,}",
            flush=True,
        )

        sample_ids: list[str] = []
        row_arrays: list[np.ndarray] = []
        max_idx = max(max(selected_indices), sample_idx)

        for line_no, line in enumerate(fh, start=2):
            if not line.strip():
                continue

            parts = line.rstrip("\n\r").split(sep)
            if len(parts) <= max_idx:
                raise RuntimeError(
                    f"Line {line_no} has {len(parts):,} fields, but expected at least {max_idx + 1:,}."
                )

            sample_id = str(parts[sample_idx]).strip()
            sample_ids.append(sample_id)

            values = np.empty(len(selected_indices), dtype=np.float64)
            for j, idx in enumerate(selected_indices):
                values[j] = parse_float(parts[idx])
            row_arrays.append(values)

            if progress_every > 0 and len(sample_ids) % progress_every == 0:
                print(f"[DunedinPACE20k] parsed samples: {len(sample_ids):,}", flush=True)

    if not row_arrays:
        raise RuntimeError("No samples were read from the matrix.")

    mat_sample_by_probe = np.vstack(row_arrays)
    # Biolearn expects CpG rows x sample columns.
    mat_probe_by_sample = mat_sample_by_probe.T

    dnam = pd.DataFrame(
        mat_probe_by_sample,
        index=pd.Index(selected_probes, name="CpG"),
        columns=pd.Index(sample_ids, name="Sample_ID"),
    )
    dnam.index = dnam.index.astype(str)
    dnam.columns = dnam.columns.astype(str)

    info = {
        "orientation": "sample_rows",
        "sample_id_col": sample_id_col,
        "matrix_cpg_n": matrix_cpg_n,
        "selected_measured_probe_n": len(selected_probes),
        "sample_n": len(sample_ids),
    }
    return dnam, info


def load_cpg_rows_subset(
    path: Path,
    sep: str,
    wanted_probes: set[str],
    probe_id_col: str,
    chunksize: int,
) -> tuple[pd.DataFrame, dict]:
    frames: list[pd.DataFrame] = []
    total_cpg_rows = 0

    for i, chunk in enumerate(pd.read_csv(path, sep=sep, chunksize=chunksize, low_memory=False), start=1):
        chunk[probe_id_col] = chunk[probe_id_col].astype(str)
        total_cpg_rows += int(chunk[probe_id_col].map(is_cpg_like).sum())
        sub = chunk[chunk[probe_id_col].isin(wanted_probes)].copy()
        if sub.empty:
            continue
        sub = sub.set_index(probe_id_col)
        sub.index = sub.index.astype(str)
        sub = sub.apply(pd.to_numeric, errors="coerce")
        frames.append(sub)
        print(f"[DunedinPACE20k] cpg-row chunks processed: {i:,}; selected so far: {sum(f.shape[0] for f in frames):,}", flush=True)

    if not frames:
        raise RuntimeError("None of the requested DunedinPACE probes are present as matrix rows.")

    dnam = pd.concat(frames, axis=0)
    dnam = dnam[~dnam.index.duplicated(keep="first")]
    dnam.columns = dnam.columns.astype(str)

    info = {
        "orientation": "cpg_rows",
        "sample_id_col": None,
        "matrix_cpg_n": total_cpg_rows,
        "selected_measured_probe_n": dnam.shape[0],
        "sample_n": dnam.shape[1],
    }
    return dnam, info


def load_dnam_subset(
    path: Path,
    wanted_probes: set[str],
    sep: str | None,
    chunksize: int,
    progress_every: int,
) -> tuple[pd.DataFrame, dict]:
    path = Path(path)
    if sep is None or sep == "auto":
        sep = infer_delimiter(path)

    orientation, id_col = infer_matrix_orientation(path, sep)
    print(f"[DunedinPACE20k] inferred matrix orientation: {orientation}", flush=True)

    if orientation == "sample_rows":
        return load_sample_rows_subset_streaming(path, sep, wanted_probes, id_col, progress_every)
    return load_cpg_rows_subset(path, sep, wanted_probes, id_col, chunksize)


def read_meta(meta_path: Path | None, sample_index: Sequence[str]) -> pd.DataFrame:
    sample_index = [str(x) for x in sample_index]
    if meta_path is None or not Path(meta_path).is_file():
        return pd.DataFrame(index=sample_index)

    meta = pd.read_csv(meta_path, low_memory=False)
    candidates = [
        "Sample_ID", "SampleID", "sample_id", "sample", "Sample", "ID",
        "Name", "Sample_Name", "Sample_Name1", "Sample_Name2",
    ]
    id_col = next((c for c in candidates if c in meta.columns), meta.columns[0])
    meta[id_col] = meta[id_col].astype(str)
    meta = meta.set_index(id_col)
    meta.index = meta.index.astype(str)
    meta = meta.reindex(sample_index)
    return meta


def coerce_prediction(pred, sample_index: Sequence[str]) -> pd.Series:
    sample_index = pd.Index([str(x) for x in sample_index])

    if isinstance(pred, pd.Series):
        s = pred.copy()
    elif isinstance(pred, pd.DataFrame):
        preferred = ["DunedinPACE", "Predicted", "predicted", "Prediction", "score"]
        col = next((c for c in preferred if c in pred.columns), None)
        if col is None:
            numeric_cols = pred.select_dtypes(include=[np.number]).columns.tolist()
            if not numeric_cols:
                raise RuntimeError(f"Could not find a numeric prediction column in {list(pred.columns)}")
            col = numeric_cols[0]
        s = pred[col].copy()
    else:
        raise RuntimeError(f"Unsupported prediction type: {type(pred)}")

    s.index = s.index.astype(str)
    return pd.to_numeric(s.reindex(sample_index), errors="coerce")


def write_list(path: Path, values: Iterable[str], gz: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            for v in sorted(set(map(str, values))):
                fh.write(v + "\n")
    else:
        with open(path, "w", encoding="utf-8") as fh:
            for v in sorted(set(map(str, values))):
                fh.write(v + "\n")


def update_clock_csv(clock_csv: Path, pred_df: pd.DataFrame, replace_existing: bool) -> Path:
    clock_csv = Path(clock_csv)
    if not clock_csv.is_file():
        raise FileNotFoundError(f"clock.csv not found: {clock_csv}")

    clock = pd.read_csv(clock_csv, low_memory=False)
    id_candidates = ["Sample_ID", "SampleID", "sample_id", "sample", "Sample", "ID", "Name"]
    id_col = next((c for c in id_candidates if c in clock.columns), clock.columns[0])

    clock[id_col] = clock[id_col].astype(str)
    pred = pred_df.rename(columns={"Sample_ID": id_col})[[id_col, "DunedinPACE_20k"]].copy()
    pred[id_col] = pred[id_col].astype(str)

    out = clock.merge(pred, on=id_col, how="left")
    if replace_existing:
        out["DunedinPACE_original_lowmem"] = out["DunedinPACE"] if "DunedinPACE" in out.columns else np.nan
        out["DunedinPACE"] = out["DunedinPACE_20k"]

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = clock_csv.with_suffix(clock_csv.suffix + f".bak_before_dunedinpace20k_{ts}")
    shutil.copy2(clock_csv, backup)
    out.to_csv(clock_csv, index=False)
    return backup



def patch_biolearn_dunedinpace_writable_arrays():
    """Patch Biolearn DunedinPACE normalization for pandas/numpy read-only arrays.

    Some pandas/numpy combinations return read-only arrays from DataFrame.values.
    Biolearn's DunedinPACE quantile normalization performs in-place assignment.
    This patch replaces that function with an equivalent implementation that
    explicitly copies input data into writable float64 arrays before assignment.
    """
    import numpy as _np
    from scipy.stats import rankdata as _rankdata
    import biolearn.dunedin_pace as _dp

    def _quantile_normalize_using_target_writable(data, target_values):
        data = _np.array(data, dtype=_np.float64, copy=True)
        sorted_target = _np.sort(
            _np.array(target_values, dtype=_np.float64, copy=True)
        )

        n_target = len(sorted_target)
        if n_target == 0:
            raise ValueError("DunedinPACE target_values is empty.")

        for j in range(data.shape[1]):
            column_data = data[:, j]
            ranks = _rankdata(column_data)
            rank_floor_values = _np.floor(ranks).astype(int)

            # Biolearn's original implementation uses 1-based rank positions.
            rank_floor_values = _np.clip(rank_floor_values, 1, n_target)

            has_decimal_above_0_4 = (ranks - rank_floor_values) > 0.4

            idx_low = rank_floor_values - 1
            idx_high = _np.clip(rank_floor_values, 0, n_target - 1)

            normalized = sorted_target[idx_low].copy()

            if _np.any(has_decimal_above_0_4):
                normalized[has_decimal_above_0_4] = 0.5 * (
                    sorted_target[idx_low[has_decimal_above_0_4]]
                    + sorted_target[idx_high[has_decimal_above_0_4]]
                )

            data[:, j] = normalized

        return data

    _dp.quantile_normalize_using_target = _quantile_normalize_using_target_writable

    # Be explicit: the function object used by dunedin_pace_normalization
    # resolves globals from the biolearn.dunedin_pace module dictionary.
    if hasattr(_dp, "dunedin_pace_normalization"):
        _dp.dunedin_pace_normalization.__globals__[
            "quantile_normalize_using_target"
        ] = _quantile_normalize_using_target_writable

    print(
        "[DunedinPACE20k] Patched Biolearn DunedinPACE normalization "
        "to use writable NumPy arrays."
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path, help="Collapsed beta matrix")
    parser.add_argument("--meta", type=Path, default=None, help="Optional metadata CSV")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--clock-csv", type=Path, default=None, help="Existing Biolearn clock.csv to update")
    parser.add_argument("--replace-existing", action="store_true", help="Replace DunedinPACE in clock.csv with DunedinPACE_20k")
    parser.add_argument("--imputation-method", default="default", help="default, none, averaging, dunedin, sesame_450k, etc.")
    parser.add_argument("--sep", default="auto")
    parser.add_argument("--chunksize", type=int, default=50000)
    parser.add_argument("--progress-every", type=int, default=25)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    core, gold, wanted = get_dunedinpace_probe_sets()
    print(f"[DunedinPACE20k] core CpGs from model.methylation_sites(): {len(core):,}", flush=True)
    print(f"[DunedinPACE20k] gold/background probes from DunedinPACE_Gold_Means.csv: {len(gold):,}", flush=True)
    print(f"[DunedinPACE20k] union requested probes: {len(wanted):,}", flush=True)

    dnam, matrix_info = load_dnam_subset(args.input, wanted, args.sep, args.chunksize, args.progress_every)
    sample_index = list(dnam.columns.astype(str))
    meta = read_meta(args.meta, sample_index)

    measured = set(dnam.index.astype(str))
    core_present = core & measured
    gold_present = gold & measured
    union_present = wanted & measured

    print(f"[DunedinPACE20k] samples: {len(sample_index):,}", flush=True)
    print(f"[DunedinPACE20k] measured requested probes: {len(union_present):,} / {len(wanted):,}", flush=True)
    print(f"[DunedinPACE20k] measured core CpGs: {len(core_present):,} / {len(core):,}", flush=True)
    print(f"[DunedinPACE20k] measured gold/background probes: {len(gold_present):,} / {len(gold):,}", flush=True)

    # Make writable float64 copy because DunedinPACE normalization may modify arrays internally.
    dnam = dnam.astype("float64", copy=True)

    data = GeoData(metadata=meta, dnam=dnam)
    gallery = ModelGallery()
    imputation = None if args.imputation_method == "default" else args.imputation_method
    model = gallery.get("DunedinPACE", imputation_method=imputation)

    print("[DunedinPACE20k] running Biolearn prediction...", flush=True)
    patch_biolearn_dunedinpace_writable_arrays()
    pred = model.predict(data)
    s = coerce_prediction(pred, sample_index)

    pred_df = pd.DataFrame({"Sample_ID": sample_index, "DunedinPACE_20k": s.values})
    pred_path = args.out_dir / "DunedinPACE_20k.predictions.csv"
    pred_df.to_csv(pred_path, index=False)

    coverage_rows = [
        {"metric": "input", "value": str(args.input)},
        {"metric": "matrix_orientation", "value": matrix_info["orientation"]},
        {"metric": "matrix_cpg_n", "value": matrix_info["matrix_cpg_n"]},
        {"metric": "sample_n", "value": len(sample_index)},
        {"metric": "dunedinpace_core_cpg_n", "value": len(core)},
        {"metric": "dunedinpace_gold_background_probe_n", "value": len(gold)},
        {"metric": "dunedinpace_requested_union_probe_n", "value": len(wanted)},
        {"metric": "measured_core_cpg_n", "value": len(core_present)},
        {"metric": "missing_core_cpg_n", "value": len(core - measured)},
        {"metric": "measured_gold_background_probe_n", "value": len(gold_present)},
        {"metric": "missing_gold_background_probe_n", "value": len(gold - measured)},
        {"metric": "measured_requested_union_probe_n", "value": len(union_present)},
        {"metric": "missing_requested_union_probe_n", "value": len(wanted - measured)},
        {"metric": "imputation_method", "value": args.imputation_method},
        {"metric": "clock_csv_updated", "value": bool(args.clock_csv)},
        {"metric": "replace_existing_DunedinPACE", "value": bool(args.replace_existing)},
    ]
    coverage_path = args.out_dir / "DunedinPACE_20k.probe_coverage.csv"
    pd.DataFrame(coverage_rows).to_csv(coverage_path, index=False)

    write_list(args.out_dir / "DunedinPACE_20k.missing_core_cpgs.txt", core - measured, gz=False)
    write_list(args.out_dir / "DunedinPACE_20k.missing_gold_probes.txt.gz", gold - measured, gz=True)

    metadata = {
        "script": "run_dunedinpace_20k_biolearn.py",
        "input": str(args.input),
        "out_dir": str(args.out_dir),
        "biolearn_model": "DunedinPACE",
        "imputation_method": args.imputation_method,
        "matrix_info": matrix_info,
        "core_cpg_n": len(core),
        "gold_background_probe_n": len(gold),
        "requested_union_probe_n": len(wanted),
        "measured_core_cpg_n": len(core_present),
        "measured_gold_background_probe_n": len(gold_present),
        "measured_requested_union_probe_n": len(union_present),
    }
    with open(args.out_dir / "DunedinPACE_20k.metadata.json", "w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2)

    if args.clock_csv is not None:
        backup = update_clock_csv(args.clock_csv, pred_df, replace_existing=args.replace_existing)
        print(f"[DunedinPACE20k] updated clock.csv: {args.clock_csv}", flush=True)
        print(f"[DunedinPACE20k] backup: {backup}", flush=True)

    print(f"[DunedinPACE20k] wrote: {pred_path}", flush=True)
    print(f"[DunedinPACE20k] wrote: {coverage_path}", flush=True)
    print("[DunedinPACE20k] done", flush=True)


if __name__ == "__main__":
    main()
