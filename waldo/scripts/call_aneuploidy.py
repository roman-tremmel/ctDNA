#!/usr/bin/env python3
"""WALDO-style within-sample chromosome-arm aneuploidy test.

For a test sample, each window's expected read count is predicted from
its OWN within-sample cluster total (test_sample's own reads summed over
all windows in that window's cluster) scaled by the window's reference
share-of-cluster ratio (mean/SD learned in build_window_clusters.py from
euploid samples). Because both the "expectation" and the "actual" values
come from the same sample's own read counts (only the ratio comes from
the reference panel), a global depth/batch shift in the test sample
cancels out - this is the "within-sample" comparison Douville et al.
2018 contrast with between-sample normalized-count comparisons (which
mfast-seqs/ uses, and which they note is more sensitive to batch effects).

Per chromosome-arm test statistic (an original approximation of the
paper's "sum the means and variances of all the clusters represented on
that chromosome arm" - see waldo/README.md for caveats):

    A_actual(arm)   = sum of observed reads over the arm's windows
    A_expected(arm) = sum over the arm's windows of
                       ref_mean_ratio_w * (test sample's total reads in w's cluster)
    Var(arm)        = sum over the arm's windows of
                       (ref_sd_ratio_w * (test sample's total reads in w's cluster))^2
    Z(arm)          = (A_actual(arm) - A_expected(arm)) / sqrt(Var(arm))

Requires `bedtools` on PATH (to map windows to chromosome arms).
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def map_windows_to_arms(windows_bed, arms_bed, min_overlap_frac=0.5):
    """Return a DataFrame window_id -> arm, keeping only windows with a
    single arm covering >= min_overlap_frac of the window (drops the
    handful of windows straddling a centromere)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = Path(tmpdir) / "window_arm.bed"
        with open(out_path, "w") as out:
            subprocess.run(
                ["bedtools", "intersect", "-wa", "-wb", "-f", str(min_overlap_frac),
                 "-a", windows_bed, "-b", arms_bed],
                check=True, stdout=out,
            )
        rows = []
        with open(out_path) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                window_id = fields[3]
                arm = fields[7]
                rows.append((window_id, arm))
    df = pd.DataFrame(rows, columns=["window_id", "arm"]).drop_duplicates("window_id")
    return df.set_index("window_id")["arm"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="test sample window_counts.tsv")
    ap.add_argument("--baseline", required=True, help="window_baseline.tsv from build_window_clusters.py")
    ap.add_argument("--windows-bed", required=True)
    ap.add_argument("--arms-bed", required=True)
    ap.add_argument("--excluded-arms", nargs="*", default=[])
    ap.add_argument("--z-cutoff", type=float, default=3.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    counts = pd.read_csv(args.counts, sep="\t").set_index("window_id")
    baseline = pd.read_csv(args.baseline, sep="\t").set_index("window_id")
    window_arm = map_windows_to_arms(args.windows_bed, args.arms_bed)
    all_arms = sorted(set(window_arm.unique()) - set(args.excluded_arms))

    df = counts.join(baseline[["cluster_id", "ref_mean_ratio", "ref_sd_ratio"]], how="inner")
    df = df.join(window_arm, how="inner")
    df = df[df["cluster_id"] != -1]
    df = df[~df["arm"].isin(args.excluded_arms)]
    df = df.dropna(subset=["ref_mean_ratio", "ref_sd_ratio"])
    # A window whose cluster has zero reference variance for its share-of-
    # cluster ratio (most often a singleton cluster, where the ratio is
    # trivially always 1.0) carries no usable information for this test and
    # is dropped. If dropping leaves an arm with no windows at all, that arm
    # is reported as "no_data" below rather than silently disappearing.
    df = df[df["ref_sd_ratio"] > 0]

    cluster_totals = df.groupby("cluster_id")["raw_count"].transform("sum")
    df["expected_reads"] = df["ref_mean_ratio"] * cluster_totals
    df["expected_sd"] = df["ref_sd_ratio"] * cluster_totals

    arm_stats = df.groupby("arm").apply(
        lambda g: pd.Series({
            "n_windows": len(g),
            "actual_reads": g["raw_count"].sum(),
            "expected_reads": g["expected_reads"].sum(),
            "expected_variance": (g["expected_sd"] ** 2).sum(),
        }),
        include_groups=False,
    )
    # Reindex against every arm the windows BED maps to (minus excluded
    # arms) so an arm that lost all of its windows upstream is reported
    # explicitly as having no usable data, instead of just vanishing.
    arm_stats = arm_stats.reindex(all_arms)
    arm_stats["n_windows"] = arm_stats["n_windows"].fillna(0)
    arm_stats["z_score"] = (
        (arm_stats["actual_reads"] - arm_stats["expected_reads"])
        / np.sqrt(arm_stats["expected_variance"])
    )
    arm_stats["call"] = np.select(
        [arm_stats["n_windows"] == 0,
         arm_stats["z_score"] >= args.z_cutoff,
         arm_stats["z_score"] <= -args.z_cutoff],
        ["no_data", "gain", "loss"],
        default="normal",
    )
    arm_stats = arm_stats.sort_values(
        "z_score", key=lambda s: s.abs(), ascending=False, na_position="last"
    )
    arm_stats.to_csv(args.out, sep="\t", index_label="arm")

    n_no_data = (arm_stats["call"] == "no_data").sum()
    n_aberrant = arm_stats["call"].isin(["gain", "loss"]).sum()
    print(f"{n_aberrant} / {len(arm_stats)} arms called gained/lost "
          f"(|Z| >= {args.z_cutoff})")
    if n_no_data:
        print(f"WARNING: {n_no_data} arm(s) had no usable windows left after "
              "cluster/baseline filtering and are reported as 'no_data' - see "
              "waldo/README.md (singleton clusters, low reference sample count).",
              file=sys.stderr)
    print(arm_stats.head(10).to_string())


if __name__ == "__main__":
    main()
