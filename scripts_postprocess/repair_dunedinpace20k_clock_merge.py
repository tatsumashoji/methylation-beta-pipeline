#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


ID_CANDIDATES = [
    "Sample_ID",
    "SampleID",
    "sample_id",
    "sample",
    "Sample",
    "Name",
    "Name_Original",
    "ID",
    "Unnamed: 0",
]

PRED_CANDIDATES = [
    "DunedinPACE_20k",
    "DunedinPACE",
    "Predicted",
    "prediction",
    "predicted",
    "score",
]


def normalize_id_basic(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def canonical_id(x) -> str:
    """Canonicalize sample IDs for safer matching.

    Examples:
      YS_438_Blood -> YS438
      YS_438       -> YS438
      KAZU_01      -> KAZU01
      KAZU-01      -> KAZU01
    """
    s = normalize_id_basic(x)
    if not s:
        return ""

    s = s.strip()

    # Remove common blood suffixes.
    s = re.sub(r"([_\-\.\s]*)blood$", "", s, flags=re.IGNORECASE)

    # Remove trailing .0 from Excel-like numeric IDs.
    s = re.sub(r"\.0$", "", s)

    # Uppercase and remove non-alphanumeric characters.
    s = s.upper()
    s = re.sub(r"[^A-Z0-9]", "", s)

    return s


def choose_prediction_column(pred: pd.DataFrame) -> str:
    for col in PRED_CANDIDATES:
        if col in pred.columns:
            s = pd.to_numeric(pred[col], errors="coerce")
            if s.notna().sum() > 0:
                return col

    numeric_cols = []
    for col in pred.columns:
        s = pd.to_numeric(pred[col], errors="coerce")
        if s.notna().sum() > 0:
            numeric_cols.append((col, int(s.notna().sum())))

    if len(numeric_cols) == 1:
        return numeric_cols[0][0]

    raise ValueError(
        "Could not determine prediction column. "
        f"Numeric candidates: {numeric_cols}. Columns: {list(pred.columns)}"
    )


def candidate_id_columns(df: pd.DataFrame) -> list[str]:
    cols = [c for c in ID_CANDIDATES if c in df.columns]

    for col in df.columns:
        if col in cols:
            continue

        if df[col].dtype == object:
            vals = df[col].dropna().astype(str)
            if len(vals) == 0:
                continue

            # Avoid obvious numeric-value columns.
            numeric_like = pd.to_numeric(vals.head(50), errors="coerce").notna().mean()
            if numeric_like < 0.5:
                cols.append(col)

    return cols


def best_id_overlap(clock: pd.DataFrame, pred: pd.DataFrame):
    """Find best matching ID columns and normalization mode."""
    best = None

    modes = [
        ("exact", normalize_id_basic),
        ("canonical_strip_blood", canonical_id),
    ]

    for ccol in candidate_id_columns(clock):
        for pcol in candidate_id_columns(pred):
            for mode_name, norm_func in modes:
                cset = {norm_func(x) for x in clock[ccol].dropna()}
                pset = {norm_func(x) for x in pred[pcol].dropna()}
                cset.discard("")
                pset.discard("")

                overlap = len(cset & pset)

                if best is None or overlap > best["overlap"]:
                    best = {
                        "overlap": overlap,
                        "clock_col": ccol,
                        "pred_col": pcol,
                        "mode": mode_name,
                        "norm_func": norm_func,
                    }

    return best


def build_prediction_map(
    pred: pd.DataFrame,
    pred_id_col: str,
    pred_value_col: str,
    norm_func,
) -> dict[str, float]:
    tmp = pred[[pred_id_col, pred_value_col]].copy()
    tmp["_merge_id"] = tmp[pred_id_col].map(norm_func)
    tmp["_pred_value"] = pd.to_numeric(tmp[pred_value_col], errors="coerce")

    tmp = tmp[(tmp["_merge_id"] != "") & tmp["_pred_value"].notna()].copy()

    # If duplicate IDs exist, keep the first non-NA value.
    # Also report duplicates for diagnostics.
    duplicated = tmp["_merge_id"].duplicated().sum()
    if duplicated > 0:
        print(f"[Repair] WARNING: duplicate normalized prediction IDs: {duplicated}")

    tmp = tmp.drop_duplicates("_merge_id", keep="first")

    return tmp.set_index("_merge_id")["_pred_value"].to_dict()


def write_unmatched_report(
    out_dir: Path,
    clock: pd.DataFrame,
    pred: pd.DataFrame,
    ccol: str | None,
    pcol: str | None,
    norm_func,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    if ccol is None or pcol is None:
        return

    clock_ids_raw = clock[ccol].map(normalize_id_basic)
    pred_ids_raw = pred[pcol].map(normalize_id_basic)

    clock_ids_norm = clock[ccol].map(norm_func)
    pred_ids_norm = pred[pcol].map(norm_func)

    clock_norm_set = {x for x in clock_ids_norm if x}
    pred_norm_set = {x for x in pred_ids_norm if x}

    missing_norm = sorted(clock_norm_set - pred_norm_set)
    extra_norm = sorted(pred_norm_set - clock_norm_set)

    pd.DataFrame({"clock_normalized_id_not_in_predictions": missing_norm}).to_csv(
        out_dir / "DunedinPACE_20k.clock_ids_missing_in_predictions.csv",
        index=False,
    )

    pd.DataFrame({"prediction_normalized_id_not_in_clock": extra_norm}).to_csv(
        out_dir / "DunedinPACE_20k.prediction_ids_extra.csv",
        index=False,
    )

    # Detailed raw-ID table for manual inspection.
    clock_detail = pd.DataFrame(
        {
            "clock_raw_id": clock_ids_raw,
            "clock_normalized_id": clock_ids_norm,
        }
    )
    clock_detail["matched_in_predictions"] = clock_detail["clock_normalized_id"].isin(pred_norm_set)
    clock_detail.to_csv(
        out_dir / "DunedinPACE_20k.clock_id_match_detail.csv",
        index=False,
    )

    pred_detail = pd.DataFrame(
        {
            "prediction_raw_id": pred_ids_raw,
            "prediction_normalized_id": pred_ids_norm,
        }
    )
    pred_detail["matched_in_clock"] = pred_detail["prediction_normalized_id"].isin(clock_norm_set)
    pred_detail.to_csv(
        out_dir / "DunedinPACE_20k.prediction_id_match_detail.csv",
        index=False,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clock-csv", required=True, type=Path)
    parser.add_argument("--pred-csv", required=True, type=Path)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--replace-existing", action="store_true")
    parser.add_argument(
        "--allow-row-order",
        action="store_true",
        help=(
            "Use row order only when ID overlap is insufficient AND row counts match. "
            "Do not use this when clock.csv and predictions.csv have different row counts."
        ),
    )
    args = parser.parse_args()

    if not args.clock_csv.is_file():
        raise FileNotFoundError(f"clock.csv not found: {args.clock_csv}")

    if not args.pred_csv.is_file():
        raise FileNotFoundError(f"prediction file not found: {args.pred_csv}")

    clock = pd.read_csv(args.clock_csv)
    pred = pd.read_csv(args.pred_csv)

    pred_col = choose_prediction_column(pred)
    pred_values = pd.to_numeric(pred[pred_col], errors="coerce")

    print("[Repair] clock:", args.clock_csv)
    print("[Repair] pred:", args.pred_csv)
    print("[Repair] clock rows:", len(clock))
    print("[Repair] pred rows:", len(pred))
    print("[Repair] prediction column:", pred_col)
    print("[Repair] pred nonNA:", int(pred_values.notna().sum()))
    print("[Repair] pred mean:", pred_values.mean())
    print("[Repair] pred sd:", pred_values.std())

    if pred_values.notna().sum() == 0:
        raise SystemExit("Prediction column is all NA; cannot repair clock.csv.")

    out = args.out or args.clock_csv

    if "DunedinPACE" in clock.columns and "DunedinPACE_original_lowmem" not in clock.columns:
        clock["DunedinPACE_original_lowmem"] = clock["DunedinPACE"]

    best = best_id_overlap(clock, pred)

    if best is None:
        raise SystemExit("No candidate ID columns found.")

    overlap = best["overlap"]
    ccol = best["clock_col"]
    pcol = best["pred_col"]
    mode = best["mode"]
    norm_func = best["norm_func"]

    print(f"[Repair] best ID overlap: {overlap} using clock.{ccol} vs pred.{pcol} ({mode})")

    new_values = pd.Series(
        np.full(len(clock), np.nan, dtype=float),
        index=clock.index,
        dtype="float64",
    )

    used_method = None

    # EPICv2 may have clock rows=495 and predictions=496.
    # Therefore, if overlap covers most clock rows, ID merge is valid.
    threshold = max(1, int(len(clock) * 0.80))

    if overlap >= threshold:
        pred_map = build_prediction_map(
            pred=pred,
            pred_id_col=pcol,
            pred_value_col=pred_col,
            norm_func=norm_func,
        )

        new_values = clock[ccol].map(norm_func).map(pred_map)
        used_method = f"id_merge:{ccol}<->{pcol}:{mode}"

    elif args.allow_row_order and len(clock) == len(pred):
        print("[Repair] using row-order replacement")
        new_values = pred_values.reset_index(drop=True)
        used_method = "row_order"

    else:
        write_unmatched_report(
            args.clock_csv.parent,
            clock,
            pred,
            ccol,
            pcol,
            norm_func,
        )
        raise SystemExit(
            "Could not safely merge predictions into clock.csv. "
            f"Best overlap was {overlap}/{len(clock)} using {ccol} vs {pcol} ({mode}). "
            "Diagnostic files were written to the clock output directory."
        )

    clock["DunedinPACE_20k"] = pd.to_numeric(new_values, errors="coerce")

    if args.replace_existing:
        clock["DunedinPACE"] = clock["DunedinPACE_20k"]

    s = pd.to_numeric(clock["DunedinPACE_20k"], errors="coerce")

    print("[Repair] method:", used_method)
    print("[Repair] output nonNA:", int(s.notna().sum()))
    print("[Repair] output mean:", s.mean())
    print("[Repair] output sd:", s.std())

    if s.notna().sum() == 0:
        write_unmatched_report(
            args.clock_csv.parent,
            clock,
            pred,
            ccol,
            pcol,
            norm_func,
        )
        raise SystemExit("Repaired DunedinPACE_20k is all NA. Aborting.")

    # Backup only after successful repair is prepared.
    if out == args.clock_csv:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup = args.clock_csv.with_name(
            args.clock_csv.name + f".bak_before_successful_dunedinpace20k_repair_{ts}"
        )
        args.clock_csv.rename(backup)
        print("[Repair] backup:", backup)

    clock.to_csv(out, index=False)
    print("[Repair] wrote:", out)


if __name__ == "__main__":
    main()
