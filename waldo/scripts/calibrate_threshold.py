#!/usr/bin/env python3
"""Derive a per-arm, data-calibrated significance threshold from a panel
of presumed-euploid samples, following Douville et al. 2018's approach:
"A statistically significant gain of a single chromosome arm was defined
as one whose Z-score was >4 above the maximum Z-score observed in the
677 normal WBC samples. Similarly, a ... loss ... <-4 below the minimum
Z-score observed."

Runs the same within-sample cluster test as call_aneuploidy.py on every
sample in the panel (ideally distinct from, or at least not dominated
by, the samples used to build the clusters themselves - the paper
reused its 677-WBC panel for both), then for each arm sets:

    gain_threshold = max(Z over the panel) + margin
    loss_threshold = min(Z over the panel) - margin

Usage:
    calibrate_threshold.py \
        --counts control1.window_counts.tsv control2.window_counts.tsv ... \
        --clusters clusters.tsv \
        --windows-bed resources/windows.500kb.hg38.bed \
        --arms-bed resources/chrom_arms.hg38.bed \
        --margin 4.0 \
        --out threshold_table.tsv
"""
import argparse
import sys

import pandas as pd

from call_aneuploidy import map_windows_to_arms, robust_cluster_moments, scaled_reads


def arm_z_scores(counts_path, members_by_window, window_arm, excluded_arms, outlier_p):
    scaled, _ = scaled_reads(counts_path)
    moments = robust_cluster_moments(scaled, members_by_window, outlier_p)

    rows = []
    for window_id, (mu, var, _n) in moments.items():
        arm = window_arm.get(window_id)
        if arm is None or arm in excluded_arms:
            continue
        r = scaled.get(window_id)
        if pd.isna(r):
            continue
        rows.append((arm, r, mu, var))
    df = pd.DataFrame(rows, columns=["arm", "R", "mu", "var"])
    if df.empty:
        return pd.Series(dtype=float)

    agg = df.groupby("arm").agg(actual=("R", "sum"), expected=("mu", "sum"),
                                 variance=("var", "sum"))
    return (agg["actual"] - agg["expected"]) / agg["variance"].pow(0.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", nargs="+", required=True,
                     help="window_counts.tsv for a panel of presumed-euploid samples")
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--windows-bed", required=True)
    ap.add_argument("--arms-bed", required=True)
    ap.add_argument("--excluded-arms", nargs="*", default=[])
    ap.add_argument("--outlier-p-threshold", type=float, default=0.01)
    ap.add_argument("--margin", type=float, default=4.0,
                     help="Douville et al. 2018 use a margin of 4 (in Z-score units) "
                          "beyond the panel's observed extremes")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    window_arm = map_windows_to_arms(args.windows_bed, args.arms_bed)
    members_by_window = (
        pd.read_csv(args.clusters, sep="\t").groupby("window_id")["member_id"]
        .apply(list).to_dict()
    )

    all_z = []
    for i, path in enumerate(args.counts):
        z = arm_z_scores(path, members_by_window, window_arm, args.excluded_arms,
                          args.outlier_p_threshold)
        z.name = path
        all_z.append(z)
        print(f"  scored {i + 1}/{len(args.counts)} reference samples", file=sys.stderr)

    z_matrix = pd.concat(all_z, axis=1)  # arms x reference samples
    thresholds = pd.DataFrame({
        "n_reference_samples": z_matrix.count(axis=1),
        "observed_min_z": z_matrix.min(axis=1),
        "observed_max_z": z_matrix.max(axis=1),
        "gain_threshold": z_matrix.max(axis=1) + args.margin,
        "loss_threshold": z_matrix.min(axis=1) - args.margin,
    })
    thresholds.to_csv(args.out, sep="\t", index_label="arm")
    print(f"Wrote per-arm thresholds (margin={args.margin}) to {args.out}")

    if (thresholds["n_reference_samples"] < 20).any():
        print("WARNING: fewer than 20 reference samples for some arms - the "
              "paper calibrated this from 677 normal WBC samples. With a small "
              "panel, 'max/min observed + margin' is a poor estimate of the true "
              "tail and will likely under- or over-call. Use more reference "
              "samples before trusting these thresholds.", file=sys.stderr)


if __name__ == "__main__":
    main()
