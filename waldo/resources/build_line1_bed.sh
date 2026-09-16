#!/usr/bin/env bash
# Build a BED file of LINE-1 (L1) repeat element loci for GRCh38, which the
# pipeline uses to decide which aligned reads are "on-target" mFAST-SeqS
# amplicon reads.
#
# mFAST-SeqS amplifies the ~500,000 LINE-1 copies scattered across the
# genome with a single primer pair; after alignment, reads are assigned to
# the LINE-1 locus (and therefore chromosome arm) they best match. The
# UCSC RepeatMasker track is the standard source for L1 element
# coordinates (this is also what tools such as WISECONDOR/ichorCNA-style
# repeat masking use, just filtered to the L1 family instead of masked out).
#
# This script is not run automatically (it downloads ~500 MB from UCSC) -
# run it once per reference genome build and point config.yaml at the
# resulting BED file.
set -euo pipefail

OUTDIR="$(dirname "$0")"
OUT="${1:-$OUTDIR/line1_elements.hg38.bed}"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "Downloading UCSC RepeatMasker track (hg38)..." >&2
curl -fsSL \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz" \
  -o "$TMP/rmsk.txt.gz"

echo "Filtering for LINE/L1 repeat family and writing BED..." >&2
# rmsk.txt columns (1-based): 6=genoName 7=genoStart 8=genoEnd 10=strand
# 11=repName 12=repClass (e.g. "LINE") 13=repFamily (e.g. "L1", "L2", "CR1")
# - filtering must be on repFamily ($13), not repClass ($12): repClass groups
# all LINEs together, repFamily is what actually says "L1" vs. other LINEs.
zcat "$TMP/rmsk.txt.gz" \
  | awk -F'\t' 'BEGIN{OFS="\t"} $13=="L1" {print $6, $7, $8, $11, 0, $10}' \
  | sort -k1,1 -k2,2n \
  > "$OUT"

echo "Wrote $(wc -l < "$OUT") L1 element intervals to $OUT" >&2
