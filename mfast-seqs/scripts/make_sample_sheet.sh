#!/usr/bin/env bash
# Build a sample sheet (sample_name<TAB>fastq_path), one row per FASTQ,
# for use with run_pipeline_array.sbatch.
#
# Assumes standard Illumina bcl2fastq single-end naming:
#   <SampleName>_S<N>_L<LLL>_R1_001.fastq.gz
# e.g. V91-01_S1_L001_R1_001.fastq.gz -> sample name "V91-01"
#
# NOTE: if a sample was split across multiple lanes (L001, L002, ...),
# this produces one row per lane file with the SAME sample name - decide
# whether you want to concatenate those lane FASTQs into one file per
# sample first (`cat L001 L002 > merged.fastq.gz`), since run_pipeline.sh
# processes one FASTQ per sample-name/output-directory as-is.
#
# Usage: make_sample_sheet.sh <fastq_dir> [glob] > samples.tsv
#   e.g. make_sample_sheet.sh /project/cfDNA_Schroth/test_samples/Fastq
set -euo pipefail

FASTQ_DIR="$1"
GLOB="${2:-*_R1_*.fastq.gz}"

shopt -s nullglob
for fq in "$FASTQ_DIR"/$GLOB; do
  base=$(basename "$fq")
  sample=$(echo "$base" | sed -E 's/_S[0-9]+_L[0-9]+_R1_[0-9]+\.fastq\.gz$//')
  printf '%s\t%s\n' "$sample" "$fq"
done
