#!/usr/bin/env python3
"""Learn window clusters and a within-sample reference baseline from a
panel of euploid reference samples.

This is an original, good-faith reimplementation of the *idea* behind
Douville et al. 2018 (PNAS 115(8):1871-1876) - "group the 500-kb intervals
with similar read depths into clusters" - not a byte-exact reproduction:
the paper's precise clustering algorithm and statistical model are in its
SI Appendix, which is not available here. See waldo/README.md for the
exact approximation this script makes and why.

Approach:
  1. From N euploid reference samples, build a windows x samples matrix of
     LINE-1 read fractions (normalized to each sample's total usable reads).
  2. Cluster windows by how correlated their fractions are across the
     reference panel (hierarchical clustering on 1 - Pearson correlation).
     Windows that track together are assumed to share local amplification
     efficiency (GC content, primer accessibility, etc.) - the same
     rationale the paper gives for why some genomic intervals "track
     together."
  3. Within EACH reference sample, compute each window's share of its own
     cluster's total reads (a within-sample ratio, immune to between-run
     depth/batch differences). Store the mean/SD of that ratio per window
     across the reference panel - this is the baseline a test sample's
     own within-sample ratios get compared against in call_aneuploidy.py.

Usage:
    build_window_clusters.py \
        --counts control1.window_counts.tsv control2.window_counts.tsv ... \
        --n-clusters 40 \
        --out-clusters clusters.tsv \
        --out-baseline window_baseline.tsv
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def warn_about_arm_pure_clusters(clusters_out, arms_bed):
    """A cluster whose windows all map to the SAME chromosome arm cannot
    provide any evidence about that arm's aneuploidy status: in
    call_aneuploidy.py, a cluster's expected-vs-actual reads always sum to
    exactly zero net deviation across the windows that make it up (the
    reference ratios sum to 1 by construction), so if a cluster IS one
    arm, that arm's "expected" total is forced to equal its own "actual"
    total. Clusters need to span multiple arms to be informative.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        windows_path = Path(tmpdir) / "windows.bed"
        clusters_out.reset_index()[["chrom", "start", "end", "window_id"]].to_csv(
            windows_path, sep="\t", header=False, index=False
        )
        out_path = Path(tmpdir) / "window_arm.bed"
        with open(out_path, "w") as out:
            subprocess.run(
                ["bedtools", "intersect", "-wa", "-wb", "-f", "0.5",
                 "-a", str(windows_path), "-b", arms_bed],
                check=True, stdout=out,
            )
        rows = []
        with open(out_path) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                rows.append((fields[3], fields[7]))
    window_arm = pd.DataFrame(rows, columns=["window_id", "arm"]).drop_duplicates("window_id")
    merged = clusters_out.reset_index().merge(window_arm, on="window_id", how="inner")
    merged = merged[merged["cluster_id"] != -1]
    n_arms_per_cluster = merged.groupby("cluster_id")["arm"].nunique()
    pure_clusters = n_arms_per_cluster[n_arms_per_cluster == 1].index
    if len(pure_clusters):
        print(f"WARNING: {len(pure_clusters)} cluster(s) are fully contained within a "
              f"single chromosome arm (cluster_id(s): {list(pure_clusters)}) and will "
              "carry no evidence about that arm's aneuploidy status. Consider using "
              "fewer clusters (larger, more likely to span multiple arms) or a larger "
              "reference panel.", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", nargs="+", required=True,
                     help="window_counts.tsv files from count_reads_by_window.py, "
                          "one per euploid reference sample")
    ap.add_argument("--n-clusters", type=int, default=40)
    ap.add_argument("--min-mean-fraction", type=float, default=1e-6,
                     help="windows with a lower mean normalized fraction across "
                          "references are treated as too sparsely covered by "
                          "LINE-1 reads (e.g. centromeric gaps) and excluded")
    ap.add_argument("--out-clusters", required=True)
    ap.add_argument("--out-baseline", required=True)
    ap.add_argument("--arms-bed",
                     help="optional: chrom_arms BED, used only to warn about clusters "
                          "that end up fully contained within a single chromosome arm "
                          "(such a cluster is structurally blind to that arm's "
                          "aneuploidy status - see waldo/README.md)")
    args = ap.parse_args()

    frames = []
    meta = None
    for path in args.counts:
        df = pd.read_csv(path, sep="\t")
        if meta is None:
            meta = df[["window_id", "chrom", "start", "end"]].set_index("window_id")
        frames.append(df.set_index("window_id")["normalized_fraction"].rename(path))

    matrix = pd.concat(frames, axis=1)  # windows x samples
    n_samples = matrix.shape[1]
    if n_samples < 5:
        print(f"WARNING: only {n_samples} reference samples given; correlation-based "
              "window clustering is unreliable with this few replicates. The paper's "
              "own description implies a much larger reference panel used once to "
              "define clusters. Treat cluster assignments here as provisional.",
              file=sys.stderr)

    mean_fraction = matrix.mean(axis=1)
    included_windows = mean_fraction[mean_fraction >= args.min_mean_fraction].index
    excluded_windows = mean_fraction.index.difference(included_windows)
    print(f"{len(excluded_windows)} / {len(mean_fraction)} windows excluded "
          f"(mean normalized fraction < {args.min_mean_fraction}, likely no "
          "LINE-1 coverage - gaps/centromeres/telomeres)", file=sys.stderr)

    sub = matrix.loc[included_windows]

    # Correlate windows' fraction profiles across the reference panel.
    corr = sub.T.corr()
    corr = corr.fillna(0.0)
    dist = 1.0 - corr.values
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2  # enforce exact symmetry against float noise
    condensed = squareform(dist, checks=False)

    n_clusters = min(args.n_clusters, len(included_windows))
    Z = linkage(condensed, method="average")
    cluster_ids = fcluster(Z, t=n_clusters, criterion="maxclust")

    clusters = pd.Series(cluster_ids, index=included_windows, name="cluster_id")

    # Within-sample ratio: each window's share of its own cluster's total
    # reads, computed independently in each reference sample.
    ratio_frames = []
    for col in sub.columns:
        sample_fracs = sub[col]
        cluster_totals = sample_fracs.groupby(clusters).transform("sum")
        ratio = (sample_fracs / cluster_totals.replace(0, np.nan)).rename(col)
        ratio_frames.append(ratio)
    ratios = pd.concat(ratio_frames, axis=1)

    baseline = pd.DataFrame({
        "cluster_id": clusters,
        "ref_mean_ratio": ratios.mean(axis=1),
        "ref_sd_ratio": ratios.std(axis=1, ddof=1),
        "n_ref_samples": ratios.count(axis=1),
    })
    baseline = baseline.join(meta)

    clusters_out = meta.loc[included_windows].copy()
    clusters_out["cluster_id"] = clusters
    excl_out = meta.loc[excluded_windows].copy()
    excl_out["cluster_id"] = -1
    clusters_out = pd.concat([clusters_out, excl_out]).sort_index()

    if args.arms_bed:
        warn_about_arm_pure_clusters(clusters_out, args.arms_bed)

    clusters_out.to_csv(args.out_clusters, sep="\t", index_label="window_id")
    baseline.to_csv(args.out_baseline, sep="\t", index_label="window_id")

    print(f"Wrote {n_clusters} clusters over {len(included_windows)} windows to "
          f"{args.out_clusters}", file=sys.stderr)
    print(f"Wrote per-window baseline ratios to {args.out_baseline}", file=sys.stderr)


if __name__ == "__main__":
    main()
