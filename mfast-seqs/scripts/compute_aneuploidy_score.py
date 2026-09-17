#!/usr/bin/env python3
"""Compute the genome-wide aneuploidy score for one sample.

Reproduces the final steps of the paper's Methods:
  1. Per included chromosome arm, compute a Z-score of the sample's
     normalized LINE-1 read fraction relative to the healthy-control
     panel-of-normals (mean, SD).
  2. Square and sum the Z-scores over all included arms (acrocentric
     short arms 13p/14p/15p/21p/22p and chrY are excluded for lack of
     LINE-1 elements) - this raw sum is chi-square-like and has an
     expected value of roughly the number of included arms (~39) even
     for a genuinely healthy sample, NOT 0.
  3. Standardize that raw sum a SECOND time, against the mean/SD of the
     same raw sum computed (leave-one-out) across the control panel
     itself (build_control_baseline.py's raw_score_mean/raw_score_sd) -
     this is what the paper means by comparing the summed Z-scores "to
     genome-wide squared and summed values of the controls", and is
     what makes the final score cluster near 0 for healthy samples and
     meaningful against a cutoff like 5 ("the number of standard
     deviations the sample is deviating from the healthy controls").
     Skipping this step (comparing the raw chi-square-like sum directly
     to the cutoff) would flag essentially every sample, including
     healthy controls, as "high".

Usage:
    compute_aneuploidy_score.py \
        --counts sample.arm_counts.tsv \
        --baseline baseline.tsv \
        --cutoff 5.0 \
        --out sample.aneuploidy.tsv
"""
import argparse
import sys

import numpy as np
import pandas as pd


def read_baseline_meta(baseline_path):
    """Pull the '# key\tvalue' summary lines build_control_baseline.py
    appends after the per-arm table (raw_score_mean, raw_score_sd)."""
    meta = {}
    with open(baseline_path) as fh:
        for line in fh:
            if line.startswith("#"):
                key, _, value = line[1:].strip().partition("\t")
                meta[key.strip()] = float(value)
    missing = {"raw_score_mean", "raw_score_sd"} - meta.keys()
    if missing:
        raise SystemExit(
            f"ERROR: baseline file '{baseline_path}' is missing {missing} - "
            "it looks like it was built with an older version of "
            "build_control_baseline.py. Rebuild it."
        )
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sample = pd.read_csv(args.counts, sep="\t").set_index("arm")
    baseline = pd.read_csv(args.baseline, sep="\t", comment="#").set_index("arm")
    baseline_meta = read_baseline_meta(args.baseline)

    merged = sample.join(baseline, how="inner", rsuffix="_baseline")
    if len(merged) != len(baseline):
        missing = set(baseline.index) - set(sample.index)
        print(f"WARNING: sample is missing arms present in baseline: {missing}",
              file=sys.stderr)

    included = merged[merged["excluded"].astype(bool) == False]  # noqa: E712

    z = (included["normalized_fraction"] - included["mean"]) / included["sd"]
    included = included.assign(z_score=z, z_squared=z ** 2)

    raw_score = included["z_squared"].sum()
    score = (raw_score - baseline_meta["raw_score_mean"]) / baseline_meta["raw_score_sd"]
    is_high = score >= args.cutoff

    out_df = included[["raw_count", "total_reads", "normalized_fraction",
                        "mean", "sd", "z_score", "z_squared"]].sort_values(
        "z_squared", ascending=False
    )
    out_df.to_csv(args.out, sep="\t", index_label="arm")

    with open(args.out, "a") as fh:
        fh.write(f"\n# raw_chi2_sum\t{raw_score:.4f}\n")
        fh.write(f"# control_raw_score_mean\t{baseline_meta['raw_score_mean']:.4f}\n")
        fh.write(f"# control_raw_score_sd\t{baseline_meta['raw_score_sd']:.4f}\n")
        fh.write(f"# aneuploidy_score\t{score:.4f}\n")
        fh.write(f"# high (>= {args.cutoff})\t{is_high}\n")

    print(f"Aneuploidy score: {score:.3f}  (raw chi2 sum={raw_score:.1f}, "
          f"cutoff={args.cutoff}, high={is_high})")
    print("Top contributing arms:")
    print(out_df.head(5).to_string())


if __name__ == "__main__":
    main()
