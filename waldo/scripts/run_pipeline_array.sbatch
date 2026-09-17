#!/usr/bin/env bash
#SBATCH --job-name=waldo
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#
# SLURM array job: one task per sample, running run_pipeline.sh (trim ->
# align -> per-window counting). Works for both reference/control and
# case samples - clustering, calibration, and scoring happen afterwards,
# once every array task has finished.
#
# No default --output/--error here on purpose: this script lives in the
# git clone, but logs must not land inside it (see ../../README.md
# "Nothing is ever written inside the git clone") - always pass explicit
# --output/--error pointing at your external $RUNDIR, e.g.:
#
#   mkdir -p "$RUNDIR"/logs
#   "$REPO/waldo/scripts/make_sample_sheet.sh" /path/to/fastq_dir > "$RUNDIR/samples.tsv"
#   sbatch --array=1-$(wc -l < "$RUNDIR/samples.tsv") \
#     --output="$RUNDIR/logs/%x_%A_%a.out" --error="$RUNDIR/logs/%x_%A_%a.err" \
#     "$REPO/waldo/scripts/run_pipeline_array.sbatch" \
#     "$RUNDIR/samples.tsv" "$RUNDIR/results" "$RUNDIR/config.yaml"
#
# See ../README.md "How to run" for the full sequence.
#
# Conda root/env are overridable without editing this file, e.g.:
#   sbatch --export=ALL,CONDA_ROOT=/other/path,CONDA_ENV_NAME=my-env ...
#   sbatch --export=ALL,CONDA_ENV_PATH=/full/path/to/envs/my-env ...
set -eo pipefail
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-ctdna-aneuploidy}"
# Activate by full PATH, not by name: name-based activation depends on
# conda finding the env in its envs_dirs registry, which a batch job's
# shell doesn't always have set up the same way an interactive login
# shell does (e.g. a .condarc envs_dirs entry that isn't sourced here).
# A full path sidesteps that lookup entirely. Override directly with
# CONDA_ENV_PATH if your env doesn't live under $CONDA_ROOT/envs/.
CONDA_ENV_PATH="${CONDA_ENV_PATH:-$CONDA_ROOT/envs/$CONDA_ENV_NAME}"

SAMPLE_SHEET="$1"
OUTDIR_ROOT="$2"
CONFIG="$3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Activate the conda env. NOTE: deliberately no `set -u` above - conda's
# own activation scripts reference unset variables internally and abort
# under strict mode, which is why "source conda.sh; conda activate ..."
# can silently fail to activate anything in a batch job even though the
# exact same commands work fine in an interactive shell.
source "$CONDA_ROOT/etc/profile.d/conda.sh"
if [[ -f "$CONDA_ROOT/etc/profile.d/mamba.sh" ]]; then
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$CONDA_ROOT}"
  source "$CONDA_ROOT/etc/profile.d/mamba.sh"
  mamba activate "$CONDA_ENV_PATH"
else
  conda activate "$CONDA_ENV_PATH"
fi

if ! command -v bwa >/dev/null 2>&1; then
  echo "ERROR: bwa still not on PATH after activating '$CONDA_ENV_PATH'." >&2
  echo "Check that path exists, or override with CONDA_ENV_PATH (see ../README.md 'How to run')." >&2
  exit 1
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
