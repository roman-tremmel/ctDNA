#!/usr/bin/env bash
#SBATCH --job-name=waldo
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# SLURM array job: one task per sample, running run_pipeline.sh (trim ->
# align -> per-window counting). Works for both reference/control and
# case samples - clustering, calibration, and scoring happen afterwards,
# once every array task has finished.
#
# Usage:
#   mkdir -p logs
#   scripts/make_sample_sheet.sh /project/cfDNA_Schroth/test_samples/Fastq > samples.tsv
#   sbatch --array=1-$(wc -l < samples.tsv) \
#     scripts/run_pipeline_array.sbatch samples.tsv results config/config.yaml
set -euo pipefail

SAMPLE_SHEET="$1"
OUTDIR_ROOT="$2"
CONFIG="$3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Activate the conda env if it isn't already active in this shell - adjust
# the path/env name to match your Miniforge install (see ../../README.md
# "Installing without root / admin rights: conda or mamba").
if ! command -v bwa >/dev/null 2>&1; then
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
  conda activate ctdna-aneuploidy
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_SHEET")
if [[ -z "$LINE" ]]; then
  echo "ERROR: no line $SLURM_ARRAY_TASK_ID in $SAMPLE_SHEET (check --array range)" >&2
  exit 1
fi
SAMPLE=$(cut -f1 <<< "$LINE")
FASTQ=$(cut -f2 <<< "$LINE")

echo "[task $SLURM_ARRAY_TASK_ID] sample=$SAMPLE fastq=$FASTQ"
mkdir -p "$OUTDIR_ROOT"

"$SCRIPT_DIR/run_pipeline.sh" "$SAMPLE" "$FASTQ" "$OUTDIR_ROOT/$SAMPLE" "$CONFIG"
