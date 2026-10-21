#!/usr/bin/env python3
"""Compare compiled tables, bit for bit, with every upstream data file."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from netCDF4 import Dataset
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("executable", type=Path)
p.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2] / "TUV-x")
p.add_argument("--manifest", type=Path, default=Path(__file__).resolve().parents[1] / "data/manifest.json")
a = p.parse_args()
manifest = json.loads(a.manifest.read_text())
compiled = subprocess.check_output([str(a.executable.resolve())])
offset = 0
for entry in manifest["files"]:
    path = a.source / entry["path"]
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == entry["sha256"], path
    if path.suffix == ".nc":
        with Dataset(path) as nc:
            nc.set_auto_maskandscale(False)
            raw = b"".join(np.asarray(v[:], dtype=np.float64).tobytes() for v in nc.variables.values())
    assert compiled[offset:offset+len(raw)] == raw, path
    offset += len(raw)
assert offset == len(compiled)
print(f"{len(manifest['files'])} datasets: {offset} compiled bytes verified exactly")
