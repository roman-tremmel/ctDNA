#!/usr/bin/env python3
"""Build fixed-size genomic windows across GRCh38, the unit WALDO groups
into read-depth "clusters" (Douville et al. 2018, PNAS 115(8):1871-1876
describe 4,361 windows of 500 kb genome-wide).

Output columns (BED, 0-based half-open): chrom start end window_id

Uses the same hardcoded GRCh38 chromosome lengths as
resources/build_chrom_arms_bed.py (duplicated here to keep this pipeline
independent of ../mfast-seqs/).
"""
import argparse

CHROM_LENGTHS = {
    "chr1":  248956422, "chr2":  242193529, "chr3":  198295559,
    "chr4":  190214555, "chr5":  181538259, "chr6":  170805979,
    "chr7":  159345973, "chr8":  145138636, "chr9":  138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16":  90338345, "chr17":  83257441, "chr18":  80373285,
    "chr19":  58617616, "chr20":  64444167, "chr21":  46709983,
    "chr22":  50818468, "chrX":  156040895, "chrY":   57227415,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", default="resources/windows.500kb.hg38.bed")
    ap.add_argument("--window-size", type=int, default=500_000)
    args = ap.parse_args()

    n = 0
    with open(args.output, "w") as fh:
        for chrom, length in CHROM_LENGTHS.items():
            start = 0
            while start < length:
                end = min(start + args.window_size, length)
                fh.write(f"{chrom}\t{start}\t{end}\twin_{chrom}_{start}\n")
                start = end
                n += 1

    print(f"Wrote {n} windows ({args.window_size} bp) to {args.output}")


if __name__ == "__main__":
    main()
