#!/usr/bin/env python3

import argparse
import hashlib
import subprocess
import re
from pathlib import Path

import numpy as np
import pandas as pd


CLOCK_CANDIDATES = {
    "Horvath": ["Horvathv1", "Horvath"],
    "Hannum": ["Hannum"],
    "PhenoAge": ["PhenoAge"],
    "GrimAgeV2": ["GrimAgeV2"],
    "DunedinPACE": ["DunedinPACE"],
    "GrimAgeV1": ["GrimAgeV1"],
}



# ---------------------------------------------------------------------
# Low-memory CpG subsetting utilities
# ---------------------------------------------------------------------

CPG_ID_PATTERN = re.compile(r"^(cg|ch)[A-Za-z0-9_:.-]+$", re.IGNORECASE)


def _looks_like_cpg(x):
    if not isinstance(x, str):
        return False
    return bool(CPG_ID_PATTERN.match(x))


def _collect_cpg_ids(obj, out, seen=None, depth=0, max_depth=7):
    """
    Recursively collect CpG-like IDs from Biolearn model objects.
    This is intentionally generic because Biolearn model classes store
    coefficient CpGs in different attributes depending on model type.
    """
    if seen is None:
        seen = set()

    if obj is None or depth > max_depth:
        return out

    oid = id(obj)
    if oid in seen:
        return out
    seen.add(oid)

    if isinstance(obj, str):
        if _looks_like_cpg(obj):
            out.add(obj)
        return out

    if isinstance(obj, pd.Index):
        for x in obj.astype(str):
            if _looks_like_cpg(x):
                out.add(x)
        return out

    if isinstance(obj, pd.Series):
        _collect_cpg_ids(obj.index, out, seen, depth + 1, max_depth)
        if obj.dtype == "object":
            for x in obj.dropna().astype(str).head(200000):
                if _looks_like_cpg(x):
                    out.add(x)
        return out

    if isinstance(obj, pd.DataFrame):
        _collect_cpg_ids(obj.index, out, seen, depth + 1, max_depth)
        _collect_cpg_ids(pd.Index(obj.columns.astype(str)), out, seen, depth + 1, max_depth)

        # Also scan object columns because some model definitions store CpGs
        # as values in columns such as "CpG", "cpg", "marker", etc.
        for col in obj.columns:
            if str(col).lower() in {"cpg", "cpgs", "cpgid", "cpg_id", "probe", "probeid", "name"}:
                for x in obj[col].dropna().astype(str):
                    if _looks_like_cpg(x):
                        out.add(x)
        return out

    if isinstance(obj, dict):
        for k, v in obj.items():
            _collect_cpg_ids(k, out, seen, depth + 1, max_depth)
            _collect_cpg_ids(v, out, seen, depth + 1, max_depth)
        return out

    if isinstance(obj, (list, tuple, set)):
        for x in obj:
            _collect_cpg_ids(x, out, seen, depth + 1, max_depth)
        return out

    if isinstance(obj, np.ndarray):
        if obj.dtype.kind in {"U", "S", "O"}:
            flat = obj.ravel()
            for x in flat[:200000]:
                if isinstance(x, bytes):
                    x = x.decode(errors="ignore")
                else:
                    x = str(x)
                if _looks_like_cpg(x):
                    out.add(x)
        return out

    if hasattr(obj, "__dict__"):
        _collect_cpg_ids(vars(obj), out, seen, depth + 1, max_depth)

    return out


def get_required_cpgs_for_biolearn_models(clocks, imputation_method):
    """
    Low-memory CpG extraction.

    Do NOT recursively scan the whole Biolearn model object.
    Some Biolearn models or imputation wrappers contain large reference objects,
    and recursively scanning them can cause OOM.

    Instead, use the Biolearn model API:
      model.methylation_sites()
    """
    from biolearn.model_gallery import ModelGallery

    gallery = ModelGallery()
    all_cpgs = set()

    for label in clocks:
        candidates = CLOCK_CANDIDATES[label]
        model = None
        selected_model_name = None
        last_error = None

        for model_name in candidates:
            try:
                if imputation_method == "default":
                    model = gallery.get(model_name)
                else:
                    model = gallery.get(model_name, imputation_method=imputation_method)

                selected_model_name = model_name
                print(
                    f"[LowMem] Loaded model for CpG list: {label} / {model_name}",
                    flush=True,
                )
                break

            except Exception as e:
                last_error = e
                print(
                    f"[LowMem][WARN] Could not load {label} candidate {model_name}: "
                    f"{type(e).__name__}: {e}",
                    flush=True,
                )

        if model is None:
            raise RuntimeError(
                f"Could not load any model for {label}. Last error: {last_error}"
            )

        if not hasattr(model, "methylation_sites"):
            raise RuntimeError(
                f"Biolearn model {label}/{selected_model_name} does not have methylation_sites()."
            )

        sites = model.methylation_sites()

        if sites is None:
            raise RuntimeError(
                f"Biolearn model {label}/{selected_model_name} returned None from methylation_sites()."
            )

        cpgs = []
        for x in sites:
            x = str(x)
            if re.match(r"^(cg|ch)[A-Za-z0-9_:.-]+$", x, flags=re.IGNORECASE):
                cpgs.append(x)

        cpgs = sorted(set(cpgs))

        print(
            f"[LowMem] {label}: required CpGs from methylation_sites() = {len(cpgs)}",
            flush=True,
        )

        if len(cpgs) == 0:
            raise RuntimeError(
                f"No CpG IDs detected from methylation_sites() for {label}/{selected_model_name}."
            )

        # Safety guard: if this becomes huge, something is wrong.
        if len(cpgs) > 100000:
            raise RuntimeError(
                f"{label}/{selected_model_name} returned too many CpGs: {len(cpgs)}. "
                "This is likely not a clock-specific CpG list and may cause OOM."
            )

        all_cpgs.update(cpgs)

    print(f"[LowMem] Union of required Biolearn CpGs = {len(all_cpgs)}", flush=True)

    if len(all_cpgs) == 0:
        raise RuntimeError("No required CpGs were detected from Biolearn models.")

    return all_cpgs


def load_beta_matrix_lowmem(path, required_cpgs, sep="	", dtype="float32"):
    """
    Load only CpG columns required by Biolearn clocks.

    pandas.read_csv(usecols=...) can still use too much memory when the
    original file has ~1M columns. Therefore, first create a slim TSV using
    the Linux cut command, then read that slim TSV with pandas.
    """
    path = Path(path)
    print(f"[Load beta lowmem] {path}", flush=True)

    header = pd.read_csv(path, sep=sep, nrows=0)
    columns = list(header.columns)

    if len(columns) < 2:
        raise RuntimeError(f"Input matrix seems invalid: {path}")

    sample_col = columns[0]
    beta_cols = columns[1:]

    required_cpgs = set(required_cpgs)

    if len(required_cpgs) == 0:
        raise RuntimeError(
            "Could not detect required CpGs from Biolearn models. "
            "Refusing to load the full EPICv2 matrix to avoid OOM."
        )

    selected_cpgs = [c for c in beta_cols if c in required_cpgs]
    missing_required = sorted(required_cpgs - set(beta_cols))

    print(f"[LowMem] matrix CpGs in file = {len(beta_cols)}", flush=True)
    print(f"[LowMem] required CpGs detected = {len(required_cpgs)}", flush=True)
    print(f"[LowMem] required CpGs present in file = {len(selected_cpgs)}", flush=True)
    print(f"[LowMem] required CpGs absent from file = {len(missing_required)}", flush=True)

    if len(selected_cpgs) == 0:
        raise RuntimeError("None of the required Biolearn CpGs were found in the beta matrix.")

    # Save absent CpGs for QC.
    absent_path = path.parent / "biolearn_required_cpgs_absent.txt"
    with open(absent_path, "w") as f:
        for cpg in missing_required:
            f.write(str(cpg) + "\n")

    print(f"[LowMem] absent required CpG list saved: {absent_path}", flush=True)

    # Build 1-based column indices for cut.
    col_to_pos = {c: i for i, c in enumerate(columns, start=1)}
    field_numbers = [1] + [col_to_pos[c] for c in selected_cpgs]
    fields_arg = ",".join(str(i) for i in field_numbers)

    # Deterministic slim file name based on selected CpGs.
    digest = hashlib.md5(("\n".join(selected_cpgs)).encode()).hexdigest()[:12]
    slim_path = path.parent / f"{path.stem}.biolearn_required_{digest}.tsv"

    if slim_path.exists() and slim_path.stat().st_size > 0:
        print(f"[LowMem] Reusing existing slim matrix: {slim_path}", flush=True)
    else:
        print(f"[LowMem] Creating slim matrix with cut: {slim_path}", flush=True)
        print(f"[LowMem] fields to keep = {len(field_numbers)}", flush=True)

        with open(slim_path, "wb") as fout:
            result = subprocess.run(
                ["cut", "-f", fields_arg, str(path)],
                stdout=fout,
                stderr=subprocess.PIPE,
                check=False,
            )

        if result.returncode != 0:
            raise RuntimeError(
                "cut command failed:\n"
                + result.stderr.decode(errors="replace")
            )

        if not slim_path.exists() or slim_path.stat().st_size == 0:
            raise RuntimeError(f"Failed to create slim matrix: {slim_path}")

    print(f"[LowMem] Reading slim matrix with pandas: {slim_path}", flush=True)

    df = pd.read_csv(
        slim_path,
        sep=sep,
        index_col=0,
        na_values=["", "NA", "NaN", "nan"],
        low_memory=False,
    )

    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    df = df.apply(pd.to_numeric, errors="coerce").astype(dtype)
    df = df.clip(0.0, 1.0)

    print(f"[Matrix lowmem] samples={df.shape[0]}, CpGs loaded={df.shape[1]}", flush=True)
    print(f"[Matrix lowmem] missing cells={int(df.isna().sum().sum())}", flush=True)

    return df


def strip_sample_suffix(x):
    return re.sub(r"_(Blood|blood|BLOOD|Cheek|cheek|CHEEK)$", "", str(x))


def sample_suffix_priority(x):
    x = str(x).lower()

    if x.endswith("_blood"):
        return 0
    if not x.endswith(("_blood", "_cheek")):
        return 1
    if x.endswith("_cheek"):
        return 2
    return 3


def load_beta_matrix(path, sep="\t", dtype="float32"):
    print(f"[Load beta] {path}")

    df = pd.read_csv(
        path,
        sep=sep,
        index_col=0,
        na_values=["", "NA", "NaN", "nan"],
        low_memory=False,
    )

    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    df = df.apply(pd.to_numeric, errors="coerce").astype(dtype)
    df = df.clip(0.0, 1.0)

    print(f"[Matrix] samples={df.shape[0]}, CpGs={df.shape[1]}")
    print(f"[Matrix] missing cells={int(df.isna().sum().sum())}")

    return df


def detect_sample_col(df):
    for c in ["Sample_ID", "Name", "検体名", "Unnamed: 0", "sample", "Sample"]:
        if c in df.columns:
            return c
    return df.columns[0]


def detect_sex_col(df):
    for c in ["sex", "Sex", "Gender", "gender", "生物学的性別", "Sex__c"]:
        if c in df.columns:
            return c
    raise ValueError("meta.csv must contain a sex column")


def detect_age_col(df):
    for c in ["ChronologicalAge", "Chronological_Age", "age", "Age", "ChronologicalAge__c"]:
        if c in df.columns:
            return c
    raise ValueError("meta.csv must contain ChronologicalAge column")


def normalize_sex(x):
    """
    Biolearn metadata sex coding used here:
      0 = female
      1 = male
    """
    if pd.isna(x):
        return np.nan

    if isinstance(x, str):
        z = x.strip().lower()

        if z in {"0", "female", "f", "woman", "女性", "女"}:
            return 0.0
        if z in {"1", "male", "m", "man", "男性", "男"}:
            return 1.0

        try:
            x = float(z)
        except Exception:
            return np.nan

    try:
        v = float(x)
        if v == 0:
            return 0.0
        if v == 1:
            return 1.0
        return 1.0 if v >= 0.5 else 0.0
    except Exception:
        return np.nan


def load_meta_map(path):
    df = pd.read_csv(path)

    sample_col = detect_sample_col(df)
    sex_col = detect_sex_col(df)
    age_col = detect_age_col(df)

    df[sample_col] = df[sample_col].astype(str).map(strip_sample_suffix)
    df[sex_col] = df[sex_col].map(normalize_sex)
    df[age_col] = pd.to_numeric(df[age_col], errors="coerce")

    df = df.dropna(subset=[sample_col]).drop_duplicates(subset=[sample_col], keep="last")

    sex_map = dict(zip(df[sample_col], df[sex_col]))
    age_map = dict(zip(df[sample_col], df[age_col]))

    print(f"[Metadata] meta loaded: {len(df)} samples")
    print(f"[Metadata] columns: sample={sample_col}, sex={sex_col}, age={age_col}")
    print(f"[Metadata] sex counts: {pd.Series(list(sex_map.values())).value_counts(dropna=False).to_dict()}")
    print(f"[Metadata] age missing in meta: {int(pd.Series(list(age_map.values())).isna().sum())}")

    return sex_map, age_map


def build_metadata(sample_ids, sex_map, age_map):
    sample_ids = pd.Index(sample_ids.astype(str), name="Sample_ID")
    clean = sample_ids.map(strip_sample_suffix)

    meta = pd.DataFrame(index=sample_ids)
    meta["Sample_ID_clean"] = clean
    meta["sex"] = [sex_map.get(x, np.nan) for x in clean]
    meta["age"] = [age_map.get(x, np.nan) for x in clean]

    print(f"[Metadata] samples={meta.shape[0]}")
    print(f"[Metadata] missing sex={int(meta['sex'].isna().sum())}")
    print(f"[Metadata] missing age={int(meta['age'].isna().sum())}")

    if meta["sex"].isna().any():
        print("[Metadata][WARN] first missing sex:", meta.index[meta["sex"].isna()].tolist()[:30])

    if meta["age"].isna().any():
        print("[Metadata][WARN] first missing age:", meta.index[meta["age"].isna()].tolist()[:30])

    return meta


def make_dnam(df_samples_by_cpg):
    """
    Biolearn expects CpGs x samples.
    Make a writable copy to avoid read-only numpy issues in DunedinPACE.
    """
    dnam0 = df_samples_by_cpg.T

    arr = np.array(
        dnam0.to_numpy(dtype=np.float64, copy=True),
        dtype=np.float64,
        copy=True,
        order="C",
    )
    arr.setflags(write=True)

    return pd.DataFrame(
        arr,
        index=dnam0.index.astype(str),
        columns=dnam0.columns.astype(str),
    )


def patch_dunedin_pace():
    try:
        import biolearn.dunedin_pace as dp
    except Exception as e:
        print(f"[Patch][WARN] cannot import biolearn.dunedin_pace: {e}")
        return

    if getattr(dp.quantile_normalize_using_target, "_writable_patch_applied", False):
        return

    original = dp.quantile_normalize_using_target

    def patched(data, target_values):
        arr = np.array(data, dtype=np.float64, copy=True, order="C")
        arr.setflags(write=True)
        return original(arr, target_values)

    patched._writable_patch_applied = True
    dp.quantile_normalize_using_target = patched

    print("[Patch] Applied writable-copy patch for DunedinPACE.")


def coerce_prediction(pred, sample_index, label):
    preferred = {
        "GrimAgeV1": ["DNAmGrimAge", "GrimAge", "Predicted"],
        "GrimAgeV2": ["DNAmGrimAge2", "DNAmGrimAge", "GrimAgeV2", "GrimAge", "Predicted"],
        "DunedinPACE": ["DunedinPACE", "Predicted"],
        "Horvath": ["Predicted"],
        "Hannum": ["Predicted"],
        "PhenoAge": ["Predicted"],
        "SexEstimation": ["predicted_sex", "sex", "Predicted", "X", "Y"],
    }

    if isinstance(pred, pd.Series):
        s = pred.copy()

    elif isinstance(pred, pd.DataFrame):
        chosen = None

        for c in preferred.get(label, ["Predicted"]):
            if c in pred.columns:
                chosen = c
                break

        if chosen is not None:
            print(f"[Biolearn] {label}: using output column '{chosen}'")
            s = pred[chosen].copy()

        elif pred.shape[1] == 1:
            chosen = pred.columns[0]
            print(f"[Biolearn] {label}: using only output column '{chosen}'")
            s = pred.iloc[:, 0].copy()

        else:
            raise ValueError(
                f"{label}: multiple output columns but no preferred column found. "
                f"Available columns: {pred.columns.tolist()}"
            )

    else:
        s = pd.Series(np.asarray(pred).reshape(-1), index=sample_index)

    s.index = s.index.astype(str)
    sample_index = pd.Index(sample_index.astype(str))

    if set(sample_index).issubset(set(s.index)):
        s = s.reindex(sample_index)
    elif len(s) == len(sample_index):
        s = pd.Series(s.to_numpy(), index=sample_index)
    else:
        raise ValueError(
            f"{label}: prediction could not be aligned. "
            f"pred_len={len(s)}, sample_len={len(sample_index)}"
        )

    if label == "SexEstimation":
        return s

    return pd.to_numeric(s, errors="coerce")


def fill_missing_sex(df_x, meta):
    if not meta["sex"].isna().any():
        print("[Metadata] sex is complete.")
        return meta

    print("[Metadata] filling missing sex by Biolearn SexEstimation.")

    from biolearn.data_library import GeoData
    from biolearn.model_gallery import ModelGallery

    data = GeoData(metadata=meta.copy(), dnam=make_dnam(df_x))
    pred = ModelGallery().get("SexEstimation").predict(data)

    sex_pred = coerce_prediction(pred, df_x.index.astype(str), "SexEstimation")
    sex_est = sex_pred.map(normalize_sex)

    missing = meta["sex"].isna()
    failed = sex_est[missing & sex_est.isna()].index.tolist()

    if failed:
        raise RuntimeError(f"SexEstimation failed for samples: {failed[:20]}")

    meta = meta.copy()
    meta.loc[missing, "sex"] = sex_est.loc[missing].astype(float)

    print(meta["sex"].value_counts(dropna=False).to_string())
    return meta


def get_model(gallery, model_name, imputation_method):
    if imputation_method == "default":
        return gallery.get(model_name)
    return gallery.get(model_name, imputation_method=imputation_method)


def predict_clock(gallery, GeoData, df_x, meta, label, candidates, imputation_method):
    last_error = None

    for model_name in candidates:
        try:
            print(f"[Biolearn] {label}: model={model_name}, imputation={imputation_method}")

            model = get_model(gallery, model_name, imputation_method)

            data = GeoData(
                metadata=meta.copy(),
                dnam=make_dnam(df_x),
            )

            pred = model.predict(data)
            s = coerce_prediction(pred, df_x.index.astype(str), label)

            print(f"[Biolearn] {label}: OK")
            return s

        except Exception as e:
            last_error = e
            print(f"[Biolearn][WARN] {label} {model_name} failed: {type(e).__name__}: {e}")

    raise RuntimeError(f"{label}: all candidates failed. Last error: {last_error}")


def run_clocks(df_x, meta, clocks, imputation_method, strict):
    patch_dunedin_pace()

    from biolearn.data_library import GeoData
    from biolearn.model_gallery import ModelGallery

    gallery = ModelGallery()

    out = pd.DataFrame(index=df_x.index.astype(str))
    errors = []

    for label in clocks:
        try:
            candidates = CLOCK_CANDIDATES[label]

            if label in {"GrimAgeV1", "GrimAgeV2"}:
                valid = meta[["sex", "age"]].notna().all(axis=1)
                n_invalid = int((~valid).sum())

                if n_invalid:
                    bad = meta.index[~valid].tolist()[:30]
                    print(f"[Biolearn][WARN] {label}: {n_invalid} samples missing sex/age; set NaN. {bad}")

                if int(valid.sum()) == 0:
                    raise RuntimeError(f"No samples with complete sex+age for {label}.")

                pred_valid = predict_clock(
                    gallery,
                    GeoData,
                    df_x.loc[valid, :],
                    meta.loc[valid, ["sex", "age"]].copy(),
                    label,
                    candidates,
                    imputation_method,
                )

                s = pd.Series(np.nan, index=df_x.index.astype(str), dtype=float)
                s.loc[pred_valid.index.astype(str)] = pred_valid.values
                out[label] = s

            else:
                out[label] = predict_clock(
                    gallery,
                    GeoData,
                    df_x,
                    meta[["sex", "age"]].copy(),
                    label,
                    candidates,
                    imputation_method,
                )

        except Exception as e:
            msg = f"{label}: {type(e).__name__}: {e}"
            errors.append(msg)
            print("[Biolearn][ERROR]", msg)

            if strict:
                raise

            out[label] = np.nan

    if errors:
        print("[Biolearn] errors/warnings:")
        for e in errors:
            print("  -", e)

    return out


def deduplicate_preferring_blood(clock_df, out_dir):
    df = clock_df.copy()

    original = pd.Index(df.index.astype(str), name="Original_Sample_ID")
    clean = pd.Index([strip_sample_suffix(x) for x in original], name="Sample_ID")
    priority = [sample_suffix_priority(x) for x in original]

    df["__original_sample_id__"] = original
    df["__clean_sample_id__"] = clean
    df["__suffix_priority__"] = priority
    df["__original_order__"] = range(df.shape[0])

    dup = df["__clean_sample_id__"].duplicated(keep=False)

    report = df.loc[
        dup,
        ["__clean_sample_id__", "__original_sample_id__", "__suffix_priority__", "__original_order__"]
    ].copy()

    if report.shape[0] > 0:
        report = report.rename(columns={
            "__clean_sample_id__": "clean_sample_id",
            "__original_sample_id__": "original_sample_id",
            "__suffix_priority__": "suffix_priority",
            "__original_order__": "original_order",
        })

        report_path = Path(out_dir) / "duplicate_sample_resolution_report.csv"
        report.to_csv(report_path, index=False)

        print(f"[Output] duplicate resolution report: {report_path}")

    df = df.sort_values(
        by=["__clean_sample_id__", "__suffix_priority__", "__original_order__"],
        ascending=[True, True, True],
        kind="mergesort",
    )

    before = df.shape[0]
    df = df.drop_duplicates(subset=["__clean_sample_id__"], keep="first")
    after = df.shape[0]

    if before != after:
        print(f"[Output] dropped {before - after} duplicated rows, preferring _Blood.")

    df.index = pd.Index(df["__clean_sample_id__"].astype(str), name="Sample_ID")

    df = df.drop(
        columns=[
            "__original_sample_id__",
            "__clean_sample_id__",
            "__suffix_priority__",
            "__original_order__",
        ]
    )

    return df


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--input", required=True, help="collapsed_beta_matrix.txt; samples x CpGs; may contain NaN")
    ap.add_argument("--meta", required=True, help="meta.csv with sample ID, sex, ChronologicalAge")
    ap.add_argument("--out-dir", required=True)

    ap.add_argument("--clock-output", default="clock.csv")
    ap.add_argument("--dtype", default="float32")
    ap.add_argument("--sep", default="\t")

    ap.add_argument(
        "--biolearn-imputation-method",
        default="default",
        choices=["default", "none", "averaging", "dunedin", "sesame_450k"],
    )

    ap.add_argument(
        "--clocks",
        nargs="+",
        default=["Horvath", "Hannum", "PhenoAge", "GrimAgeV2", "DunedinPACE", "GrimAgeV1"],
        choices=list(CLOCK_CANDIDATES.keys()),
    )

    ap.add_argument("--strict-clocks", action="store_true")
    ap.add_argument("--allow-missing-age", action="store_true")
    ap.add_argument("--allow-missing-sex", action="store_true")
    ap.add_argument("--fill-missing-sex", action="store_true")

    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    required_cpgs = get_required_cpgs_for_biolearn_models(
        clocks=args.clocks,
        imputation_method=args.biolearn_imputation_method,
    )

    df_x = load_beta_matrix_lowmem(
        Path(args.input),
        required_cpgs=required_cpgs,
        sep=args.sep,
        dtype=args.dtype,
    )

    sex_map, age_map = load_meta_map(Path(args.meta))
    meta = build_metadata(df_x.index, sex_map, age_map)

    if meta["age"].isna().any() and not args.allow_missing_age:
        missing = meta.index[meta["age"].isna()].tolist()
        raise RuntimeError(
            f"ChronologicalAge missing for {len(missing)} samples: {missing[:30]}. "
            "Use --allow-missing-age to continue."
        )

    if meta["sex"].isna().any():
        missing = meta.index[meta["sex"].isna()].tolist()

        if args.fill_missing_sex:
            try:
                meta = fill_missing_sex(df_x, meta)
            except Exception as e:
                if args.allow_missing_sex:
                    print(
                        "[Metadata][WARN] SexEstimation failed, but continuing because "
                        "--allow-missing-sex was specified."
                    )
                    print(f"[Metadata][WARN] SexEstimation error: {type(e).__name__}: {e}")
                    print(f"[Metadata][WARN] Missing sex samples will be NaN: {missing[:30]}")
                else:
                    raise

        elif args.allow_missing_sex:
            print(
                f"[Metadata][WARN] Sex missing for {len(missing)} samples, "
                "but continuing because --allow-missing-sex was specified."
            )
            print(f"[Metadata][WARN] Missing sex samples will be NaN: {missing[:30]}")

        else:
            raise RuntimeError(
                f"Sex missing for {len(missing)} samples: {missing[:30]}. "
                "Add these samples to meta.csv, use --fill-missing-sex, "
                "or use --allow-missing-sex."
            )

    meta.to_csv(out_dir / "clock_metadata_used.csv")
    print(f"[Metadata] saved: {out_dir / 'clock_metadata_used.csv'}")

    print(f"[Biolearn] imputation method: {args.biolearn_imputation_method}")

    clock_df = run_clocks(
        df_x=df_x,
        meta=meta[["sex", "age"]].copy(),
        clocks=args.clocks,
        imputation_method=args.biolearn_imputation_method,
        strict=args.strict_clocks,
    )

    clock_df = deduplicate_preferring_blood(clock_df, out_dir)
    clock_df.index.name = "Sample_ID"

    out_path = out_dir / args.clock_output
    clock_df.to_csv(out_path)

    summary = pd.DataFrame({
        "metric": [
            "n_samples_input",
            "n_cpgs_input",
            "n_missing_cells_input",
            "biolearn_imputation_method",
            "clocks",
        ],
        "value": [
            df_x.shape[0],
            df_x.shape[1],
            int(df_x.isna().sum().sum()),
            args.biolearn_imputation_method,
            ",".join(args.clocks),
        ],
    })

    summary.to_csv(out_dir / "biolearn_default_imputation_summary.csv", index=False)

    print("[Done]")
    print(f"  clock output: {out_path}")
    print(clock_df.head())
    print("\nNA counts:")
    print(clock_df.isna().sum())


if __name__ == "__main__":
    main()
