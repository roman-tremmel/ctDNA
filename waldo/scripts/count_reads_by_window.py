#!/usr/bin/env python3
"""Count LINE-1-derived reads per 500kb genomic window from an aligned BAM.

Same on-target filtering logic as mfast-seqs/scripts/count_reads_by_arm.py
(confidently/uniquely mapped reads overlapping an annotated LINE-1 locus),
but reads are assigned to fixed-size genomic windows instead of whole
chromosome arms - the finer-grained unit WALDO clusters and tests.

Requires `samtools` and `bedtools` on PATH.
"""
import argparse
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


def run(cmd, **kwargs):
    return subprocess.run(cmd, check=True, **kwargs)


def load_windows(windows_bed):
    windows = []
    with open(windows_bed) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, start, end, window_id = line.rstrip("\n").split("\t")
            windows.append((chrom, int(start), int(end), window_id))
    return windows


def count_reads(bam, line1_bed, windows_bed, min_mapq, tmpdir):
    tmpdir = Path(tmpdir)
    filtered_bam = tmpdir / "filtered.bam"
    reads_bed = tmpdir / "reads.bed"
    ontarget_bed = tmpdir / "ontarget.bed"
    assigned_bed = tmpdir / "assigned.bed"

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

    with open(ontarget_bed, "w") as out:
        run(["bedtools", "intersect", "-u", "-a", str(reads_bed), "-b", line1_bed],
            stdout=out)

    with open(assigned_bed, "w") as out:
        run(["bedtools", "intersect", "-wa", "-wb", "-a", str(ontarget_bed), "-b", windows_bed],
            stdout=out)

    window_counts = Counter()
    seen_reads = set()
    with open(assigned_bed) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            read_key = (fields[0], fields[1], fields[2], fields[3])
            if read_key in seen_reads:
                continue
            seen_reads.add(read_key)
            window_id = fields[9]  # windows.bed col 4 -> offset 6+3
            window_counts[window_id] += 1

    return total_reads, window_counts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--line1-bed", required=True)
    ap.add_argument("--windows-bed", required=True)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    windows = load_windows(args.windows_bed)

    with tempfile.TemporaryDirectory() as tmpdir:
        total_reads, window_counts = count_reads(
            args.bam, args.line1_bed, args.windows_bed, args.min_mapq, tmpdir
        )

    with open(args.out, "w") as out:
        out.write("window_id\tchrom\tstart\tend\traw_count\ttotal_reads\tnormalized_fraction\n")
        for chrom, start, end, window_id in windows:
            raw = window_counts.get(window_id, 0)
            norm = raw / total_reads if total_reads else 0.0
            out.write(f"{window_id}\t{chrom}\t{start}\t{end}\t{raw}\t{total_reads}\t{norm:.12f}\n")

    print(f"Total usable reads: {total_reads}", file=sys.stderr)
    print(f"Wrote per-window counts to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
