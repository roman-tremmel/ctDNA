#!/usr/bin/env python3
"""WALDO's within-sample chromosome-arm aneuploidy test, following
Douville et al. 2018 (PNAS 115(8):1871-1876) SI Appendix, SI Materials
and Methods ("Identifying chromosome arm gains or losses in a test
sample"):

  1. Scale the test sample's own window read counts genome-wide
     (subtract the sample's own mean, divide by its own SD) to get
     "scaled reads" R_i per window.
  2. For each window i, using its reference-learned cluster C_i (see
     build_window_clusters.py - windows on OTHER chromosomes whose read
     depth is statistically indistinguishable from window i across a
     reference panel), estimate the cluster's mean mu_i and variance
     sigma_i^2 FROM THIS SAME TEST SAMPLE's own R values at the windows
     in C_i, by maximum likelihood with iterative outlier removal
     (any member with two-sided p < --outlier-p-threshold under the
     current Normal(mu_i, sigma_i^2) is dropped, then mu_i/sigma_i^2 are
     re-estimated; repeated until stable).
  3. Since a sum of normal random variables is normal, for a chromosome
     arm made up of windows I: sum(R_i) ~ N(sum(mu_i), sum(sigma_i^2)).
     Z = (sum(R_i) - sum(mu_i)) / sqrt(sum(sigma_i^2)); a positive Z > a
     threshold is a gain, a negative Z < -threshold is a loss.

NOTE on an ambiguity in the source text: the SI describes mu_i/sigma_i^2
as being estimated "from the genomic intervals in each ... cluster" and
illustrates this (SI Appendix Fig. S8B-C) separately for one normal WBC
sample and one aneuploid tumor sample - which is how we read it as a
per-test-sample (within-sample) computation, consistent with the
"within-sample" premise of the method (and with waldo/ being an
independent, good-faith reimplementation rather than the original
code - see waldo/README.md).

Use calibrate_threshold.py to derive a data-driven --z-cutoff (the paper
sets it 4 standard deviations beyond the most extreme value seen in a
large euploid reference panel) instead of trusting a default value.
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def map_windows_to_arms(windows_bed, arms_bed, min_overlap_frac=0.5):
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
                rows.append((fields[3], fields[7]))
    df = pd.DataFrame(rows, columns=["window_id", "arm"]).drop_duplicates("window_id")
    return df.set_index("window_id")["arm"]


def scaled_reads(counts_path):
    df = pd.read_csv(counts_path, sep="\t").set_index("window_id")
    r = df["raw_count"].astype(float)
    scaled = (r - r.mean()) / r.std(ddof=1)
    return scaled, df


def robust_cluster_moments(scaled, members_by_window, outlier_p, max_iter=20, min_members=5):
    """For each window, iteratively trim outlier cluster members (in THIS
    sample's own scaled-read values) and return the converged (mean,
    variance, n_members_used) - the within-sample mu_i, sigma_i^2."""
    from scipy import stats

    moments = {}
    for window_id, members in members_by_window.items():
        values = scaled.reindex(members).dropna().to_numpy()
        for _ in range(max_iter):
            if len(values) < min_members:
                break
            mu, sd = values.mean(), values.std(ddof=1)
            if sd == 0:
                break
            z = (values - mu) / sd
            p = 2 * stats.norm.sf(np.abs(z))
            keep = p >= outlier_p
            if keep.all():
                break
            values = values[keep]
        if len(values) >= min_members and values.std(ddof=1) > 0:
            moments[window_id] = (values.mean(), values.var(ddof=1), len(values))
    return moments


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="test sample window_counts.tsv")
    ap.add_argument("--clusters", required=True,
                     help="window_id/member_id membership TSV from build_window_clusters.py")
    ap.add_argument("--windows-bed", required=True)
    ap.add_argument("--arms-bed", required=True)
    ap.add_argument("--excluded-arms", nargs="*", default=[])
    ap.add_argument("--outlier-p-threshold", type=float, default=0.01,
                     help="two-sided p-value below which a cluster member is "
                          "trimmed as an outlier when estimating mu_i/sigma_i^2")
    ap.add_argument("--z-cutoff", type=float, default=None,
                     help="flat |Z| cutoff for gain/loss calls. Prefer "
                          "--threshold-table from calibrate_threshold.py for a "
                          "per-arm, data-calibrated cutoff.")
    ap.add_argument("--threshold-table",
                     help="per-arm threshold TSV from calibrate_threshold.py "
                          "(columns: arm, gain_threshold, loss_threshold)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.z_cutoff is None and not args.threshold_table:
        ap.error("provide --z-cutoff or --threshold-table")

    scaled, counts_df = scaled_reads(args.counts)
    window_arm = map_windows_to_arms(args.windows_bed, args.arms_bed)
    all_arms = sorted(set(window_arm.unique()) - set(args.excluded_arms))

    members_by_window = (
        pd.read_csv(args.clusters, sep="\t").groupby("window_id")["member_id"]
        .apply(list).to_dict()
    )
    moments = robust_cluster_moments(scaled, members_by_window, args.outlier_p_threshold)

    rows = []
    for window_id, (mu, var, n_members) in moments.items():
        arm = window_arm.get(window_id)
        if arm is None or arm in args.excluded_arms:
            continue
        r = scaled.get(window_id)
        if pd.isna(r):
            continue
        rows.append((window_id, arm, r, mu, var, n_members))
    win_df = pd.DataFrame(rows, columns=["window_id", "arm", "R", "mu", "var", "n_members"])

    arm_stats = win_df.groupby("arm").agg(
        n_windows=("window_id", "size"),
        actual_scaled=("R", "sum"),
        expected_scaled=("mu", "sum"),
        expected_variance=("var", "sum"),
    )
    arm_stats = arm_stats.reindex(all_arms)
    arm_stats["n_windows"] = arm_stats["n_windows"].fillna(0)
    arm_stats["z_score"] = (
        (arm_stats["actual_scaled"] - arm_stats["expected_scaled"])
        / np.sqrt(arm_stats["expected_variance"])
    )

    if args.threshold_table:
        thresh = pd.read_csv(args.threshold_table, sep="\t").set_index("arm")
        arm_stats = arm_stats.join(thresh[["gain_threshold", "loss_threshold"]])
    else:
        arm_stats["gain_threshold"] = args.z_cutoff
        arm_stats["loss_threshold"] = -args.z_cutoff

    arm_stats["call"] = np.select(
        [arm_stats["n_windows"] == 0,
         arm_stats["z_score"] >= arm_stats["gain_threshold"],
         arm_stats["z_score"] <= arm_stats["loss_threshold"]],
        ["no_data", "gain", "loss"],
        default="normal",
    )
    arm_stats = arm_stats.sort_values(
        "z_score", key=lambda s: s.abs(), ascending=False, na_position="last"
    )
    arm_stats.to_csv(args.out, sep="\t", index_label="arm")

    n_no_data = (arm_stats["call"] == "no_data").sum()
    n_aberrant = arm_stats["call"].isin(["gain", "loss"]).sum()
    print(f"{n_aberrant} / {len(arm_stats)} arms called gained/lost")
    if n_no_data:
        print(f"WARNING: {n_no_data} arm(s) had no windows with a usable cluster "
              "(too few surviving members after outlier trimming) and are "
              "reported as 'no_data'.", file=sys.stderr)
    print(arm_stats.head(10).to_string())


if __name__ == "__main__":
    main()
