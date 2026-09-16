#!/usr/bin/env python3
"""Count LINE-1-derived reads per chromosome arm from an aligned BAM.

Pipeline (mirrors the read-counting step of mFAST-SeqS):
  1. Keep confidently, uniquely mapped primary reads (MAPQ >= min_mapq,
     not unmapped/secondary/supplementary).
  2. Restrict to reads overlapping an annotated LINE-1 element (the
     amplicon's on-target loci).
  3. Assign each surviving read to the chromosome arm it falls in.
  4. Report the raw count per arm and the count normalized to the total
     number of usable (post-filter) reads in the library, i.e. the same
     "read counts per chromosome arm, normalized to total library size"
     quantity described in the paper's Methods.

Requires `samtools` and `bedtools` on PATH.

Usage:
    count_reads_by_arm.py --bam sample.sorted.bam \
        --line1-bed resources/line1_elements.hg38.bed \
        --arms-bed resources/chrom_arms.hg38.bed \
        --min-mapq 30 \
        --out sample.arm_counts.tsv
"""
import argparse
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


def run(cmd, **kwargs):
    return subprocess.run(cmd, check=True, **kwargs)


def load_arms(arms_bed):
    """Return list of (chrom, start, end, arm, excluded) sorted for scanning."""
    arms = []
    with open(arms_bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, start, end, arm, excluded = line.rstrip("\n").split("\t")
            arms.append((chrom, int(start), int(end), arm, excluded == "1"))
    return arms


def count_reads(bam, line1_bed, arms_bed, min_mapq, tmpdir):
    tmpdir = Path(tmpdir)
    filtered_bam = tmpdir / "filtered.bam"
    reads_bed = tmpdir / "reads.bed"
    ontarget_bed = tmpdir / "ontarget.bed"
    assigned_bed = tmpdir / "assigned.bed"

    # 1. Drop unmapped (0x4), secondary (0x100), supplementary (0x800) reads
    #    and enforce a minimum mapping quality.
    with open(filtered_bam, "wb") as out:
        run(["samtools", "view", "-b", "-q", str(min_mapq), "-F", "0x904", bam],
            stdout=out)

    total_reads = int(
        run(["samtools", "view", "-c", str(filtered_bam)],
            stdout=subprocess.PIPE, text=True).stdout.strip()
    )

    if total_reads == 0:
        return total_reads, Counter()

    with open(reads_bed, "w") as out:
        run(["bedtools", "bamtobed", "-i", str(filtered_bam)], stdout=out)

    # 2. Keep only reads overlapping an annotated LINE-1 element.
    with open(ontarget_bed, "w") as out:
        run(["bedtools", "intersect", "-u", "-a", str(reads_bed), "-b", line1_bed],
            stdout=out)

    # 3. Assign each on-target read to the chromosome arm it overlaps.
    with open(assigned_bed, "w") as out:
        run(["bedtools", "intersect", "-wa", "-wb", "-a", str(ontarget_bed), "-b", arms_bed],
            stdout=out)

    arm_counts = Counter()
    seen_reads = set()
    with open(assigned_bed) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            # reads.bed (bamtobed) has 6 columns: chrom,start,end,name,score,strand
            read_key = (fields[0], fields[1], fields[2], fields[3])
            if read_key in seen_reads:
                # a read spanning an arm boundary could match both arms;
                # count it once, against the first (5'-most) arm it hits.
                continue
            seen_reads.add(read_key)
            arm = fields[9]  # arms_bed col 4 -> offset 6+3
            arm_counts[arm] += 1

    return total_reads, arm_counts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--line1-bed", required=True)
    ap.add_argument("--arms-bed", required=True)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    arms = load_arms(args.arms_bed)
    excluded_by_arm = {arm: excluded for _, _, _, arm, excluded in arms}

    with tempfile.TemporaryDirectory() as tmpdir:
        total_reads, arm_counts = count_reads(
            args.bam, args.line1_bed, args.arms_bed, args.min_mapq, tmpdir
        )

    with open(args.out, "w") as out:
        out.write("arm\texcluded\traw_count\ttotal_reads\tnormalized_fraction\n")
        for arm in sorted(excluded_by_arm):
            raw = arm_counts.get(arm, 0)
            norm = raw / total_reads if total_reads else 0.0
            out.write(
                f"{arm}\t{int(excluded_by_arm[arm])}\t{raw}\t{total_reads}\t{norm:.10f}\n"
            )

    print(f"Total usable reads: {total_reads}", file=sys.stderr)
    print(f"Wrote per-arm counts to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
