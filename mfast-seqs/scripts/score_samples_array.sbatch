#!/usr/bin/env bash
#SBATCH --job-name=mfast-seqs-score
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#
# SLURM array job: one task per row in samples.tsv, scoring that sample's
# arm_counts.tsv (produced by run_pipeline_array.sbatch) against a
# baseline built with build_control_baseline.py. Runs over EVERY row,
# controls included - scoring a control against its own baseline is a
# reasonable sanity check (expect a low score), not wasted work. If you'd
# rather only score cases, submit with a narrower --array range that
# excludes the control rows' line numbers in samples.tsv.
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
#     "$REPO/mfast-seqs/scripts/score_samples_array.sbatch" \
#     "$RUNDIR/samples.tsv" "$RUNDIR/results" "$RUNDIR/results/baseline.tsv" \
#     5.0 "$REPO/mfast-seqs/scripts"
set -eo pipefail
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-ctdna-aneuploidy}"
CONDA_ENV_PATH="${CONDA_ENV_PATH:-$CONDA_ROOT/envs/$CONDA_ENV_NAME}"

SAMPLE_SHEET="$1"
RESULTS_DIR="$2"
BASELINE="$3"
CUTOFF="$4"
SCRIPT_DIR="$5"
if [[ -z "$SCRIPT_DIR" || ! -f "$SCRIPT_DIR/compute_aneuploidy_score.py" ]]; then
  echo "ERROR: 5th argument must be the mfast-seqs/scripts directory (got '$SCRIPT_DIR')" >&2
  echo "e.g. sbatch ... score_samples_array.sbatch samples.tsv results baseline.tsv 5.0 \"\$REPO/mfast-seqs/scripts\"" >&2
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
python3 "$SCRIPT_DIR/compute_aneuploidy_score.py" \
  --counts "$RESULTS_DIR/$SAMPLE_ID/$SAMPLE_ID.arm_counts.tsv" \
  --baseline "$BASELINE" \
  --cutoff "$CUTOFF" \
  --out "$RESULTS_DIR/$SAMPLE_ID.aneuploidy.tsv"
