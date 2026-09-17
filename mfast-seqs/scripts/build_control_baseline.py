#!/usr/bin/env python3
"""Build a panel-of-normals baseline from healthy-control samples.

Reproduces two levels of standardization from the paper's Methods, not
just the first:

  1. Per-arm mean/SD of the normalized LINE-1 read fraction across the
     control panel ("a Z-score per chromosome arm was calculated
     relative to healthy female controls").
  2. The mean/SD, ACROSS THE SAME CONTROL PANEL, of the sum of squared
     per-arm Z-scores itself ("the resulting Z-scores were squared and
     summed and compared to genome-wide squared and summed values of
     the controls").

Level 2 matters: a sum of ~39 squared Z-scores has an expected value of
about 39 under the null (chi-square-like, not zero), so comparing that
raw sum directly to a cutoff like 5.0 would flag essentially every
sample - including the healthy controls themselves - as "high". The
paper's final "aneuploidy score...indicates the number of standard
deviations the sample is deviating from the healthy controls" only
makes sense, and only clusters near 0 for genuinely healthy samples, as
a second-level Z-score of the raw sum against the controls' own
distribution of that same sum. compute_aneuploidy_score.py applies this
second level using the raw_score_mean/raw_score_sd this script computes.

Each control's own raw sum is computed leave-one-out (its per-arm
Z-scores use the mean/SD of the OTHER controls, not including itself),
to avoid the downward bias of scoring a sample against a baseline that
already includes it.

Usage:
    build_control_baseline.py \
        --counts control1.arm_counts.tsv control2.arm_counts.tsv ... \
        --out baseline.tsv
"""
import argparse
import sys

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", nargs="+", required=True,
                     help="arm_counts.tsv files produced by count_reads_by_arm.py, "
                          "one per healthy control sample")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if len(args.counts) < 10:
        print(f"WARNING: only {len(args.counts)} control samples given; "
              "the paper's approach and general practice for panel-of-normals "
              "Z-scoring calls for a reasonably sized reference set "
              "(dozens, ideally) to get stable per-arm mean/SD estimates "
              "and a stable raw-score mean/SD for the second-level "
              "normalization below.", file=sys.stderr)

    frames = []
    for path in args.counts:
        df = pd.read_csv(path, sep="\t")
        df = df.set_index("arm")["normalized_fraction"]
        df.name = path
        frames.append(df)

    matrix = pd.concat(frames, axis=1)  # arms x samples
    excluded = pd.read_csv(args.counts[0], sep="\t").set_index("arm")["excluded"]
    included_arms = excluded.index[~excluded.astype(bool)]

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
              file=sys.stderr)

    # Level 2: each control's own raw sum-of-squared-Z, leave-one-out.
    raw_scores = []
    for col in matrix.columns:
        others = matrix.drop(columns=col)
        loo_mean = others.mean(axis=1)
        loo_sd = others.std(axis=1, ddof=1)
        z = (matrix.loc[included_arms, col] - loo_mean[included_arms]) / loo_sd[included_arms]
        raw_scores.append((z ** 2).sum())
    raw_scores = pd.Series(raw_scores, index=matrix.columns)

    raw_score_mean = raw_scores.mean()
    raw_score_sd = raw_scores.std(ddof=1)
    print(f"Control panel's own (leave-one-out) raw scores: "
          f"mean={raw_score_mean:.2f}, sd={raw_score_sd:.2f}, "
          f"range=[{raw_scores.min():.2f}, {raw_scores.max():.2f}]", file=sys.stderr)

    baseline.to_csv(args.out, sep="\t", index_label="arm")
    with open(args.out, "a") as fh:
        fh.write(f"# raw_score_mean\t{raw_score_mean:.6f}\n")
        fh.write(f"# raw_score_sd\t{raw_score_sd:.6f}\n")

    print(f"Wrote baseline for {len(baseline)} arms from {len(frames)} controls to {args.out}")


if __name__ == "__main__":
    main()
