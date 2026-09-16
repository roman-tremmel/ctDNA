#!/usr/bin/env bash
# Strip the heterogeneity spacer (Fadrosh et al. 2014 dual-indexing scheme)
# and the LINE-1 locus-specific primer from raw mFAST-SeqS reads.
#
# Reads look like:  [0-7 random bases][fixed L1 primer][biological insert]
# cutadapt's non-anchored 5' adapter mode (-g) finds the primer wherever it
# sits and removes it plus everything preceding it, which is exactly the
# standard trick for removing heterogeneity spacers in amplicon protocols.
#
# Usage: trim_line1_reads.sh <in.fastq[.gz]> <out.trimmed.fastq.gz> <primer_seq> [max_error_rate] [min_len]
set -euo pipefail

IN_FASTQ="$1"
OUT_FASTQ="$2"
PRIMER="$3"
MAX_ERR="${4:-0.15}"
MIN_LEN="${5:-20}"

if [[ -z "$PRIMER" || "$PRIMER" == "REPLACE_WITH_REAL_L1_PRIMER_SEQUENCE" ]]; then
  echo "ERROR: no real LINE-1 primer sequence configured (config.yaml primer.forward_seq)" >&2
  exit 1
fi

cutadapt \
  -g "$PRIMER" \
  -e "$MAX_ERR" \
  --no-indels \
  --discard-untrimmed \
  -m "$MIN_LEN" \
  -o "$OUT_FASTQ" \
  "$IN_FASTQ" \
  1>&2
