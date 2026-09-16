#!/usr/bin/env python3
"""Build a panel-of-normals baseline (mean/SD per arm) from healthy-control
samples, exactly analogous to the "Z-score per chromosome arm ... relative
to healthy female controls" step in the paper's Methods.

Usage:
    build_control_baseline.py \
        --counts control1.arm_counts.tsv control2.arm_counts.tsv ... \
        --out baseline.tsv
"""
import argparse

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--counts", nargs="+", required=True,
                     help="arm_counts.tsv files produced by count_reads_by_arm.py, "
                          "one per healthy control sample")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if len(args.counts) < 10:
        print(f"WARNING: only {len(args.counts)} control samples given; "
              "the paper's approach and general practice for panel-of-normals "
              "Z-scoring calls for a reasonably sized reference set "
              "(dozens, ideally) to get stable per-arm mean/SD estimates.",
              file=__import__("sys").stderr)

    frames = []
    for path in args.counts:
        df = pd.read_csv(path, sep="\t")
        df = df.set_index("arm")["normalized_fraction"]
        df.name = path
        frames.append(df)

    matrix = pd.concat(frames, axis=1)  # arms x samples
    excluded = pd.read_csv(args.counts[0], sep="\t").set_index("arm")["excluded"]

    baseline = pd.DataFrame({
        "mean": matrix.mean(axis=1),
        "sd": matrix.std(axis=1, ddof=1),
        "n_controls": matrix.count(axis=1),
        "excluded": excluded,
    })

    zero_sd = baseline.index[(baseline["sd"] == 0) & (~baseline["excluded"].astype(bool))]
    if len(zero_sd):
        print(f"WARNING: zero SD in included arms {list(zero_sd)}; "
              "Z-scores for these arms will be undefined (NaN/inf) downstream.",
              file=__import__("sys").stderr)

    baseline.to_csv(args.out, sep="\t", index_label="arm")
    print(f"Wrote baseline for {len(baseline)} arms from {len(frames)} controls to {args.out}")


if __name__ == "__main__":
    main()
