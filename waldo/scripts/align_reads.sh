#!/usr/bin/env bash
# Align trimmed single-end mFAST-SeqS reads with bwa mem and produce a
# coordinate-sorted, indexed BAM.
#
# Usage: align_reads.sh <trimmed.fastq.gz> <reference.fa> <out.bam> [threads]
set -euo pipefail

FASTQ="$1"
REF="$2"
OUT_BAM="$3"
THREADS="${4:-4}"

if [[ ! -f "${REF}.bwt" ]]; then
  echo "ERROR: reference is not bwa-indexed. Run: bwa index ${REF}" >&2
  exit 1
fi

bwa mem -t "$THREADS" "$REF" "$FASTQ" 2>/dev/null \
  | samtools sort -@ "$THREADS" -o "$OUT_BAM" -
samtools index "$OUT_BAM"
