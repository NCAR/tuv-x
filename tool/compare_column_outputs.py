#!/usr/bin/env python3
# Copyright (C) 2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
"""Compare complete column_sweep outputs from two builds at roundoff scale.

Each radiation wavelength is normalized to at least unit incident flux; each
rate is normalized to its maximum magnitude over the heights of that column.
This permits roundoff at cancellation boundaries without an absolute tolerance
that would hide differences in small reaction rates. Requires NumPy.
"""
import argparse
from pathlib import Path

import numpy as np


def compare(reference_path, actual_path, ulps):
    with Path(reference_path).open("rb") as reference_file, Path(actual_path).open("rb") as actual_file:
        header = np.fromfile(reference_file, dtype=np.int32, count=6)
        actual_header = np.fromfile(actual_file, dtype=np.int32, count=6)
        if header.size != 6 or not np.array_equal(header, actual_header) or np.any(header < 0):
            raise ValueError("Missing or different output dimensions")
        columns, heights, wavelengths, photo, heat, dose = map(int, header)
        fields = [(name, wavelengths) for name in ("fdr", "fdn", "fup", "edr", "edn", "eup")]
        fields += [("photolysis", photo), ("heating", heat), ("dose", dose)]
        epsilon = np.finfo(np.float64).eps
        maxima = dict.fromkeys((name for name, _ in fields), 0.0)
        values = 0
        for column in range(columns):
            for name, width in fields:
                count = heights * width
                reference = np.fromfile(reference_file, dtype=np.float64, count=count)
                actual = np.fromfile(actual_file, dtype=np.float64, count=count)
                if reference.size != count or actual.size != count:
                    raise ValueError(f"Truncated {name} in column {column}")
                if not np.all(np.isfinite(reference)) or not np.all(np.isfinite(actual)):
                    raise ValueError(f"Nonfinite {name} in column {column}")
                if count == 0:
                    continue
                reference = reference.reshape((heights, width), order="F")
                actual = actual.reshape((heights, width), order="F")
                scale = np.maximum(np.max(np.abs(reference), axis=0), np.max(np.abs(actual), axis=0))
                scale = np.maximum(scale, 1.0 if name in {"fdr", "fdn", "fup", "edr", "edn", "eup"}
                                   else np.finfo(np.float64).tiny)
                error = np.abs(reference - actual) / scale / epsilon
                maximum = float(np.max(error))
                maxima[name] = max(maxima[name], maximum)
                if maximum > ulps:
                    raise ValueError(f"{name} column {column}: {maximum:.2f} eps exceeds {ulps} eps")
                values += count
        if reference_file.read(1) or actual_file.read(1):
            raise ValueError("Unexpected data after the last column")
    print(f"Compared {values:,} values across {columns} columns")
    for name, maximum in maxima.items():
        print(f"  {name}: {maximum:.2f} eps maximum normalized difference")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference")
    parser.add_argument("actual")
    parser.add_argument("--eps", type=float, default=16.0,
                        help="roundoff bound for full-model comparisons (default: 16)")
    args = parser.parse_args()
    compare(args.reference, args.actual, args.eps)
