#!/usr/bin/env python3
"""Compute the genome-wide aneuploidy score for one sample.

Reproduces the final steps of the paper's Methods:
  1. Per included chromosome arm, compute a Z-score of the sample's
     normalized LINE-1 read fraction relative to the healthy-control
     panel-of-normals (mean, SD).
  2. Square and sum the Z-scores over all included arms (acrocentric
     short arms 13p/14p/15p/21p/22p and chrY are excluded for lack of
     LINE-1 elements).
  3. Report the resulting aneuploidy score, i.e. how many standard
     deviations the sample's genome-wide profile deviates from healthy
     controls, and classify it against a cutoff (default 5, as used in
     Verschoor et al. 2023 / Belic et al. 2016).

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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sample = pd.read_csv(args.counts, sep="\t").set_index("arm")
    baseline = pd.read_csv(args.baseline, sep="\t").set_index("arm")

    merged = sample.join(baseline, how="inner", rsuffix="_baseline")
    if len(merged) != len(baseline):
        missing = set(baseline.index) - set(sample.index)
        print(f"WARNING: sample is missing arms present in baseline: {missing}",
              file=sys.stderr)

    included = merged[merged["excluded"].astype(bool) == False]  # noqa: E712

    z = (included["normalized_fraction"] - included["mean"]) / included["sd"]
    included = included.assign(z_score=z, z_squared=z ** 2)

    score = included["z_squared"].sum()
    is_high = score >= args.cutoff

    out_df = included[["raw_count", "total_reads", "normalized_fraction",
                        "mean", "sd", "z_score", "z_squared"]].sort_values(
        "z_squared", ascending=False
    )
    out_df.to_csv(args.out, sep="\t", index_label="arm")

    with open(args.out, "a") as fh:
        fh.write(f"\n# aneuploidy_score\t{score:.4f}\n")
        fh.write(f"# high (>= {args.cutoff})\t{is_high}\n")

    print(f"Aneuploidy score: {score:.3f}  (cutoff={args.cutoff}, high={is_high})")
    print("Top contributing arms:")
    print(out_df.head(5).to_string())


if __name__ == "__main__":
    main()
