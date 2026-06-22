#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="/opt/conda/envs/methylation-beta-pipeline/bin/python"
RSCRIPT_BIN="/opt/conda/envs/methylation-beta-pipeline/bin/Rscript"

timestamp=$(date +"%Y%m%d_%H%M%S")

echo "Replacing python/python3/Rscript calls in shell scripts..."
echo "Python:  ${PYTHON_BIN}"
echo "Rscript: ${RSCRIPT_BIN}"
echo

python3 - <<PY
from pathlib import Path
import re
import shutil
import subprocess

PYTHON_BIN = "${PYTHON_BIN}"
RSCRIPT_BIN = "${RSCRIPT_BIN}"
timestamp = "${timestamp}"

target_dirs = [Path("bin")]

patterns = [
    (re.compile(r'(?<![A-Za-z0-9_./-])python3(?![A-Za-z0-9_.-])'), PYTHON_BIN),
    (re.compile(r'(?<![A-Za-z0-9_./-])python(?![A-Za-z0-9_.-])'), PYTHON_BIN),
    (re.compile(r'(?<![A-Za-z0-9_./-])Rscript(?![A-Za-z0-9_.-])'), RSCRIPT_BIN),
]

def is_text_file(path: Path) -> bool:
    try:
        data = path.read_bytes()
        if b"\\0" in data:
            return False
        data.decode("utf-8")
        return True
    except Exception:
        return False

files = []
for d in target_dirs:
    if not d.exists():
        continue
    for p in d.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix == ".sh" or p.stat().st_mode & 0o111:
            files.append(p)

patched = []

for path in files:
    if not is_text_file(path):
        continue

    text = path.read_text(encoding="utf-8")
    new_text = text

    for pattern, replacement in patterns:
        new_text = pattern.sub(replacement, new_text)

    if new_text != text:
        backup = path.with_name(path.name + f".bak_interpreter_{timestamp}")
        shutil.copy2(path, backup)
        path.write_text(new_text, encoding="utf-8")
        patched.append((path, backup))

if patched:
    for path, backup in patched:
        print(f"Patched: {path}")
        print(f"Backup:  {backup}")
else:
    print("No files were changed.")
PY

echo
echo "Done."
echo
echo "Check remaining commands with:"
echo "  grep -R \"python\\|python3\\|Rscript\" bin"
