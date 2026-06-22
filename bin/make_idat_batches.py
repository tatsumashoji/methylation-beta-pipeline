#!/usr/bin/env python3

from pathlib import Path
import argparse
import os
import shutil
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch-size", type=int, default=25)
    ap.add_argument("--min-size", type=int, default=1_000_000)
    args = ap.parse_args()

    src = Path(args.src)
    out = Path(args.out)

    if not src.exists():
        raise SystemExit(f"ERROR: source directory does not exist: {src}")

    if out.exists():
        shutil.rmtree(out)

    out.mkdir(parents=True, exist_ok=True)

    all_files = [p for p in src.rglob("*") if p.is_file()]
    reds = []

    for p in all_files:
        if p.name.startswith("._"):
            continue
        if not p.name.endswith("_Red.idat"):
            continue
        if p.stat().st_size < args.min_size:
            print(f"SKIP small Red file: {p} size={p.stat().st_size}")
            continue
        reds.append(p)

    reds = sorted(reds)

    print(f"Found Red IDAT files: {len(reds)}")

    if not reds:
        print("First 30 files:")
        for p in all_files[:30]:
            print(p)
        sys.exit(1)

    pairs = []

    for red in reds:
        grn = red.with_name(red.name.replace("_Red.idat", "_Grn.idat"))

        if not grn.exists():
            print(f"WARNING: missing Grn for {red}")
            continue

        if grn.name.startswith("._"):
            continue

        if grn.stat().st_size < args.min_size:
            print(f"WARNING: small Grn file: {grn} size={grn.stat().st_size}")
            continue

        pairs.append((red, grn))

    print(f"Complete Red/Grn pairs: {len(pairs)}")

    if not pairs:
        sys.exit("ERROR: no complete IDAT pairs")

    for i in range(0, len(pairs), args.batch_size):
        bdir = out / f"batch_{i // args.batch_size + 1:03d}"
        bdir.mkdir(parents=True, exist_ok=True)

        for red, grn in pairs[i:i + args.batch_size]:
            for f in (red, grn):
                dest = bdir / f.name
                try:
                    os.link(f, dest)
                except OSError:
                    shutil.copy2(f, dest)

    for bdir in sorted(out.glob("batch_*")):
        n_red = len(list(bdir.glob("*_Red.idat")))
        n_grn = len(list(bdir.glob("*_Grn.idat")))
        print(f"{bdir}: Red={n_red}, Grn={n_grn}")


if __name__ == "__main__":
    main()
