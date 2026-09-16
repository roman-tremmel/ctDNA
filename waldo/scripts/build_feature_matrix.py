#!/usr/bin/env python3
"""Assemble a samples x chromosome-arm-Z-score feature matrix from several
call_aneuploidy.py outputs, for training or applying the SVM classifier.

Usage:
    build_feature_matrix.py \
        --scores sampleA.arm_scores.tsv sampleB.arm_scores.tsv ... \
        --names sampleA sampleB ... \
        --out feature_matrix.tsv
"""
import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", nargs="+", required=True,
                     help="call_aneuploidy.py output TSVs")
    ap.add_argument("--names", nargs="+", required=True,
                     help="sample name per --scores file, same order")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if len(args.scores) != len(args.names):
        raise SystemExit("--scores and --names must have the same length")

    rows = {}
    for path, name in zip(args.scores, args.names):
        df = pd.read_csv(path, sep="\t").set_index("arm")
        rows[name] = df["z_score"]

    matrix = pd.DataFrame(rows).T  # samples x arms
    matrix.to_csv(args.out, sep="\t", index_label="sample")
    print(f"Wrote {matrix.shape[0]} samples x {matrix.shape[1]} arms to {args.out}")


if __name__ == "__main__":
    main()
