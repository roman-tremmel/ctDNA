#!/usr/bin/env python3
"""Build a chromosome-arm BED file for GRCh38 (hg38).

Coordinates are the centromere gap boundaries and primary-assembly contig
lengths published by UCSC/GRC for GRCh38. They are hardcoded here for
convenience so the pipeline works out of the box, but they should be
cross-checked against the authoritative UCSC "gap"/"cytoBand" tables
(https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/) before any
diagnostic or publication use.

Output columns (BED, 0-based half-open):
    chrom  start  end  arm  excluded

`excluded` is "1" for arms that mFAST-SeqS drops because they carry too
few LINE-1 elements to be quantified reliably: the short arms of the
acrocentric chromosomes (13, 14, 15, 21, 22) and chromosome Y (see
Belic et al. 2016 / Verschoor et al. 2023, npj Breast Cancer 9:61).
"""
import argparse

# (chrom, length, centromere_start, centromere_end)
# Lengths: GRCh38 primary assembly. Centromere gap: UCSC hg38 "gap" track
# (approximate midpoints of the centromeric gap region).
CHROM_INFO = {
    "chr1":  (248956422, 122026460, 125184587),
    "chr2":  (242193529,  92188145,  94090557),
    "chr3":  (198295559,  90772459,  93655574),
    "chr4":  (190214555,  49712061,  51743951),
    "chr5":  (181538259,  46485901,  50059807),
    "chr6":  (170805979,  58553888,  59829934),
    "chr7":  (159345973,  58169653,  60828234),
    "chr8":  (145138636,  44033744,  45877265),
    "chr9":  (138394717,  43236167,  45518558),
    "chr10": (133797422,  39686682,  41593521),
    "chr11": (135086622,  51078349,  54425074),
    "chr12": (133275309,  34769407,  37185252),
    "chr13": (114364328,  16000000,  18051248),
    "chr14": (107043718,  16000000,  18173523),
    "chr15": (101991189,  17083673,  19725254),
    "chr16": ( 90338345,  36311158,  38265669),
    "chr17": ( 83257441,  22813679,  26616164),
    "chr18": ( 80373285,  15460899,  20861206),
    "chr19": ( 58617616,  24498980,  27190874),
    "chr20": ( 64444167,  26436232,  28780655),
    "chr21": ( 46709983,  10864561,  12915808),
    "chr22": ( 50818468,  12954788,  15026459),
    "chrX":  (156040895,  58605579,  62412542),
    "chrY":  ( 57227415,  10316944,  10544039),
}

# Short arms with too few LINE-1 elements to quantify (acrocentric chromosomes),
# plus chrY, are excluded from the aneuploidy score (Verschoor et al. 2023).
EXCLUDED_ARMS = {
    "chr13p", "chr14p", "chr15p", "chr21p", "chr22p",
    "chrYp", "chrYq",
}


def build_rows():
    rows = []
    for chrom, (length, cen_start, cen_end) in CHROM_INFO.items():
        p_arm = f"{chrom}p"
        q_arm = f"{chrom}q"
        rows.append((chrom, 0, cen_start, p_arm, int(p_arm in EXCLUDED_ARMS)))
        rows.append((chrom, cen_end, length, q_arm, int(q_arm in EXCLUDED_ARMS)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", default="resources/chrom_arms.hg38.bed")
    args = ap.parse_args()

    rows = build_rows()
    with open(args.output, "w") as fh:
        for chrom, start, end, arm, excluded in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{arm}\t{excluded}\n")
    print(f"Wrote {len(rows)} arms to {args.output}")


if __name__ == "__main__":
    main()
