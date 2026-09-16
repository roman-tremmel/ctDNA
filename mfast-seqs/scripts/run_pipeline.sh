#!/usr/bin/env bash
# End-to-end mFAST-SeqS-style aneuploidy pipeline for a single sample:
#   raw fastq -> trim spacer+primer -> align -> count reads per arm
#
# The resulting <sample>.arm_counts.tsv is the input to both
# build_control_baseline.py (for healthy controls) and
# compute_aneuploidy_score.py (for a case sample, against a baseline
# built from controls processed the same way).
#
# Usage: run_pipeline.sh <sample_name> <raw.fastq.gz> <outdir> <config.yaml>
set -euo pipefail

SAMPLE="$1"
FASTQ="$2"
OUTDIR="$3"
CONFIG="$4"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUTDIR"

# --- read config.yaml without adding a YAML dependency ---
get_cfg() { python3 -c "
import yaml, sys
with open('$CONFIG') as fh:
    cfg = yaml.safe_load(fh)
node = cfg
for key in '$1'.split('.'):
    node = node[key]
print(node)
"
}

REF_FASTA=$(get_cfg reference.fasta)
LINE1_BED=$(get_cfg reference.line1_bed)
ARMS_BED=$(get_cfg reference.chrom_arms_bed)
PRIMER=$(get_cfg primer.forward_seq)
MAX_ERR=$(get_cfg primer.max_error_rate)
MIN_LEN=$(get_cfg alignment.min_read_length)
THREADS=$(get_cfg alignment.threads)
MIN_MAPQ=$(get_cfg alignment.min_mapq)
MIN_TOTAL_READS=$(get_cfg qc.min_total_reads)

TRIMMED="$OUTDIR/${SAMPLE}.trimmed.fastq.gz"
BAM="$OUTDIR/${SAMPLE}.sorted.bam"
COUNTS="$OUTDIR/${SAMPLE}.arm_counts.tsv"

echo "[1/3] Trimming spacer + LINE-1 primer..." >&2
"$SCRIPT_DIR/trim_line1_reads.sh" "$FASTQ" "$TRIMMED" "$PRIMER" "$MAX_ERR" "$MIN_LEN"

echo "[2/3] Aligning to reference..." >&2
"$SCRIPT_DIR/align_reads.sh" "$TRIMMED" "$REF_FASTA" "$BAM" "$THREADS"

echo "[3/3] Counting LINE-1 reads per chromosome arm..." >&2
python3 "$SCRIPT_DIR/count_reads_by_arm.py" \
  --bam "$BAM" \
  --line1-bed "$LINE1_BED" \
  --arms-bed "$ARMS_BED" \
  --min-mapq "$MIN_MAPQ" \
  --out "$COUNTS"

TOTAL_READS=$(awk -F'\t' 'NR==2{print $4}' "$COUNTS")
if [[ -n "$TOTAL_READS" ]] && (( TOTAL_READS < MIN_TOTAL_READS )); then
  echo "WARNING: sample $SAMPLE has only $TOTAL_READS usable reads (< $MIN_TOTAL_READS QC threshold)" >&2
fi

echo "Done: $COUNTS" >&2
