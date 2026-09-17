#!/usr/bin/env bash
#SBATCH --job-name=waldo-score
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#
# SLURM array job: one task per row in samples.tsv, calling per-arm
# gains/losses for that sample's window_counts.tsv (produced by
# run_pipeline_array.sbatch) against cluster membership + calibrated
# thresholds (build_window_clusters.py / calibrate_threshold.py). Runs
# over EVERY row, reference samples included - scoring one against its
# own calibration panel is a reasonable sanity check, not wasted work.
# If you'd rather only score cases, submit with a narrower --array range
# that excludes the reference rows' line numbers in samples.tsv.
#
# samples.tsv's optional 3rd column (group label, see
# run_pipeline_array.sbatch) is honored the same way here, so this reads
# the correct <group>_<sample> output directory that step 4 produced.
#
# --export=ALL,CONDA_ROOT=...,CONDA_ENV_PATH=... and --output/--error
# pointed at your external $RUNDIR are required for the same reasons as
# run_pipeline_array.sbatch (see that script, or the README "How to
# run", for why - not optional). Full example:
#
#   sbatch --array=1-$(wc -l < "$RUNDIR/samples.tsv") \
#     --export=ALL,CONDA_ROOT="$CONDA_ROOT",CONDA_ENV_PATH="$CONDA_ENV_PATH" \
#     --output="$RUNDIR/logs/%x_%A_%a.out" --error="$RUNDIR/logs/%x_%A_%a.err" \
#     "$REPO/waldo/scripts/score_samples_array.sbatch" \
#     "$RUNDIR/samples.tsv" "$RUNDIR/results" "$RUNDIR/results/clusters.tsv" \
#     "$RUNDIR/results/thresholds.tsv" "$REPO/waldo/resources/windows.500kb.hg38.bed" \
#     "$REPO/waldo/resources/chrom_arms.hg38.bed" "$REPO/waldo/scripts"
set -eo pipefail
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-ctdna-aneuploidy}"
CONDA_ENV_PATH="${CONDA_ENV_PATH:-$CONDA_ROOT/envs/$CONDA_ENV_NAME}"

SAMPLE_SHEET="$1"
RESULTS_DIR="$2"
CLUSTERS="$3"
THRESHOLDS="$4"
WINDOWS_BED="$5"
ARMS_BED="$6"
SCRIPT_DIR="$7"
EXCLUDED_ARMS=(chr13p chr14p chr15p chr21p chr22p chrYp chrYq)

if [[ -z "$SCRIPT_DIR" || ! -f "$SCRIPT_DIR/call_aneuploidy.py" ]]; then
  echo "ERROR: 7th argument must be the waldo/scripts directory (got '$SCRIPT_DIR')" >&2
  echo "e.g. sbatch ... score_samples_array.sbatch samples.tsv results clusters.tsv thresholds.tsv windows.bed arms.bed \"\$REPO/waldo/scripts\"" >&2
  exit 1
fi

# Same activation approach as run_pipeline_array.sbatch - see that
# script's comments for why (no `set -u`, mamba.sh preferred, full-path
# activation instead of by-name).
source "$CONDA_ROOT/etc/profile.d/conda.sh"
if [[ -f "$CONDA_ROOT/etc/profile.d/mamba.sh" ]]; then
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$CONDA_ROOT}"
  source "$CONDA_ROOT/etc/profile.d/mamba.sh"
  mamba activate "$CONDA_ENV_PATH"
else
  conda activate "$CONDA_ENV_PATH"
fi

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_SHEET")
if [[ -z "$LINE" ]]; then
  echo "ERROR: no line $SLURM_ARRAY_TASK_ID in $SAMPLE_SHEET (check --array range)" >&2
  exit 1
fi
SAMPLE=$(cut -f1 <<< "$LINE")
GROUP=$(cut -f3 <<< "$LINE")   # empty string if samples.tsv has no 3rd column

if [[ -n "$GROUP" ]]; then
  SAMPLE_ID="${GROUP}_${SAMPLE}"
else
  SAMPLE_ID="$SAMPLE"
fi

echo "[task $SLURM_ARRAY_TASK_ID] scoring $SAMPLE_ID"
python3 "$SCRIPT_DIR/call_aneuploidy.py" \
  --counts "$RESULTS_DIR/$SAMPLE_ID/$SAMPLE_ID.window_counts.tsv" \
  --clusters "$CLUSTERS" \
  --windows-bed "$WINDOWS_BED" \
  --arms-bed "$ARMS_BED" \
  --excluded-arms "${EXCLUDED_ARMS[@]}" \
  --threshold-table "$THRESHOLDS" \
  --out "$RESULTS_DIR/$SAMPLE_ID.arm_scores.tsv"
