#!/usr/bin/env bash
#SBATCH --job-name=ctdna-step
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#
# Generic single-task SLURM job: activates the conda env (identically to
# run_pipeline_array.sbatch) and then runs whatever command you pass as
# arguments. For the non-alignment steps - build_control_baseline.py,
# compute_aneuploidy_score.py - so they run on a compute node like
# everything else, instead of directly on a login node.
#
# --export=ALL,CONDA_ROOT=...,CONDA_ENV_PATH=... is required for the
# same reason as in run_pipeline_array.sbatch (see that script, or the
# README "How to run", for why). --output/--error should point at your
# external $RUNDIR, same as for the array job. Example:
#
#   sbatch --export=ALL,CONDA_ROOT="$CONDA_ROOT",CONDA_ENV_PATH="$CONDA_ENV_PATH" \
#     --output="$RUNDIR/logs/%x_%j.out" --error="$RUNDIR/logs/%x_%j.err" \
#     "$REPO/mfast-seqs/scripts/run_step.sbatch" \
#     python3 "$REPO/mfast-seqs/scripts/build_control_baseline.py" \
#     --counts "$RUNDIR"/results/control_*/control_*.arm_counts.tsv \
#     --out "$RUNDIR/results/baseline.tsv"
#
# See ../README.md "How to run" for the full sequence.
set -eo pipefail
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-ctdna-aneuploidy}"
CONDA_ENV_PATH="${CONDA_ENV_PATH:-$CONDA_ROOT/envs/$CONDA_ENV_NAME}"

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

if [[ $# -eq 0 ]]; then
  echo "ERROR: no command given - pass one after this script, e.g.:" >&2
  echo "  sbatch ... run_step.sbatch python3 some_script.py --args" >&2
  exit 1
fi

echo "[job $SLURM_JOB_ID] running: $*"
exec "$@"
