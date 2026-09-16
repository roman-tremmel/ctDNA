#!/usr/bin/env python3
"""Learn per-window cluster membership from a euploid reference panel,
following the algorithm described in Douville et al. 2018 (PNAS
115(8):1871-1876) SI Appendix, SI Materials and Methods ("Sample
alignment and genomic interval grouping"):

  For each 500-kb genomic interval i, compare its normalized read counts
  across the reference panel to every other interval i' that lies on a
  DIFFERENT chromosome. If the two intervals' read-count vectors are not
  significantly different in mean (paired t-test, p > --mean-p-threshold)
  AND not significantly different in variance (F-test, p >
  --var-p-threshold), interval i' is added to interval i's cluster C_i.

This yields one (typically large, ~200-window in the paper) cluster PER
WINDOW, not a fixed number of mutually-exclusive clusters - a window's
cluster membership is a personalized "these windows behave like me"
neighbor list. Restricting candidates to other chromosomes is
deliberate (per the SI text: "compared ... to ... all other genomic
intervals i' that occurred on the remaining 21 autosomal chromosomes")
and is what guarantees a chromosome arm's own windows can never make up
an entire cluster - see call_aneuploidy.py, which uses these clusters
for a within-sample test that would otherwise be blind to an arm
affecting its own cluster.

With only a handful of reference samples (paper uses 7), the t/F-tests
have low power, so "not significantly different" is a lenient bar and
clusters end up large - this matches the paper's own reported ~200
windows/cluster and is expected behavior, not a bug.

Output: a long-format TSV `window_id  member_id` (one row per
membership edge), consumed by call_aneuploidy.py.
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy import stats


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", nargs="+", required=True,
                     help="window_counts.tsv files from count_reads_by_window.py, "
                          "one per euploid reference sample (paper found 7 sufficient)")
    ap.add_argument("--mean-p-threshold", type=float, default=0.05,
                     help="paired t-test p-value above which two windows' means "
                          "are considered indistinguishable")
    ap.add_argument("--var-p-threshold", type=float, default=0.05,
                     help="F-test p-value above which two windows' variances "
                          "are considered indistinguishable")
    ap.add_argument("--min-mean-fraction", type=float, default=1e-6,
                     help="windows with a lower mean normalized fraction across "
                          "references are treated as too sparsely covered by "
                          "LINE-1 reads (e.g. centromeric gaps) and excluded")
    ap.add_argument("--out", required=True, help="output long-format membership TSV")
    args = ap.parse_args()

    frames = []
    meta = None
    for path in args.counts:
        df = pd.read_csv(path, sep="\t")
        if meta is None:
            meta = df[["window_id", "chrom"]].set_index("window_id")
        frames.append(df.set_index("window_id")["normalized_fraction"].rename(path))

    matrix = pd.concat(frames, axis=1)  # windows x samples
    n_samples = matrix.shape[1]
    if n_samples < 7:
        print(f"WARNING: only {n_samples} reference samples given; the paper found "
              "7 to be sufficient for stable clustering (more did not help). Fewer "
              "than that gives the t/F-tests very low power, which will make "
              "clusters unusually large (almost everything looks 'not "
              "significantly different').", file=sys.stderr)

    mean_fraction = matrix.mean(axis=1)
    included = mean_fraction[mean_fraction >= args.min_mean_fraction].index
    excluded = mean_fraction.index.difference(included)
    print(f"{len(excluded)} / {len(mean_fraction)} windows excluded "
          f"(mean normalized fraction < {args.min_mean_fraction})", file=sys.stderr)

    sub = matrix.loc[included]
    chrom = meta.loc[included, "chrom"]
    X = sub.to_numpy()  # (n_windows, n_samples)
    window_ids = sub.index.to_numpy()
    chroms = chrom.to_numpy()
    n, m = X.shape
    df_t = m - 1

    mean = X.mean(axis=1)
    var = X.var(axis=1, ddof=1)

    with open(args.out, "w") as out:
        out.write("window_id\tmember_id\n")
        for i in range(n):
            diffs = X - X[i][None, :]  # (n, m): paired differences vs window i
            mean_diff = diffs.mean(axis=1)
            sd_diff = diffs.std(axis=1, ddof=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                t_stat = mean_diff / (sd_diff / np.sqrt(m))
            p_mean = np.where(
                sd_diff > 0,
                2 * stats.t.sf(np.abs(t_stat), df_t),
                np.where(mean_diff == 0, 1.0, 0.0),
            )

            with np.errstate(divide="ignore", invalid="ignore"):
                F = var / var[i]
            lower_tail = stats.f.cdf(F, df_t, df_t)
            p_var = np.where(
                var[i] > 0,
                2 * np.minimum(lower_tail, 1 - lower_tail),
                np.where(var == 0, 1.0, 0.0),
            )

            candidate = (chroms != chroms[i]) & (p_mean > args.mean_p_threshold) \
                & (p_var > args.var_p_threshold)
            candidate[i] = False

            members = window_ids[candidate]
            for member_id in members:
                out.write(f"{window_ids[i]}\t{member_id}\n")

            if (i + 1) % 500 == 0 or i + 1 == n:
                print(f"  ...processed {i + 1}/{n} windows", file=sys.stderr)

    sizes = pd.read_csv(args.out, sep="\t").groupby("window_id").size()
    empty = len(included) - len(sizes)
    print(f"Wrote cluster memberships for {len(included)} windows to {args.out}",
          file=sys.stderr)
    if len(sizes):
        print(f"Cluster size: mean={sizes.mean():.1f}, median={sizes.median():.0f}, "
              f"min={sizes.min()}, max={sizes.max()}", file=sys.stderr)
    if empty:
        print(f"WARNING: {empty} window(s) got an empty cluster (no other-chromosome "
              "window passed both equivalence tests) and will be excluded from "
              "call_aneuploidy.py's arm test.", file=sys.stderr)


if __name__ == "__main__":
    main()
