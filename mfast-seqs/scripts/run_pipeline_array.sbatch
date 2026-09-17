#!/usr/bin/env bash
#SBATCH --job-name=mfast-seqs
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#
# SLURM array job: one task per sample, running run_pipeline.sh (trim ->
# align -> per-arm counting). Works for both control and case samples -
# what you do with the resulting <sample>.arm_counts.tsv (feed it to
# build_control_baseline.py or compute_aneuploidy_score.py) happens
# afterwards, once every array task has finished.
#
# samples.tsv has an OPTIONAL 3rd column: a group label you add by hand
# (e.g. "control" for healthy controls, "case" or blank for everything
# else). When present, it's prefixed onto the sample name (and therefore
# the output directory/file names) as "<group>_<sample>", so a plain
# shell glob like results/control_*/control_*.arm_counts.tsv picks up
# every control sample for build_control_baseline.py without needing
# your original sample IDs to follow any particular naming convention.
#
# No default --output/--error here on purpose: this script lives in the
# git clone, but logs must not land inside it (see ../../README.md
# "Nothing is ever written inside the git clone") - always pass explicit
# --output/--error pointing at your external $RUNDIR. And --export=ALL
# is NOT optional: SLURM does not reliably hand a submitting shell's
# exported variables ($CONDA_ROOT, $CONDA_ENV_PATH) to the job on every cluster's
# configuration - without it, this script silently falls back to its
# hardcoded default ($HOME/miniforge3) and fails if your install is
# anywhere else. Full example:
#
#   mkdir -p "$RUNDIR"/logs
#   "$REPO/mfast-seqs/scripts/make_sample_sheet.sh" /path/to/fastq_dir > "$RUNDIR/samples.tsv"
#   sbatch --array=1-$(wc -l < "$RUNDIR/samples.tsv") \
#     --export=ALL,CONDA_ROOT="$CONDA_ROOT",CONDA_ENV_PATH="$CONDA_ENV_PATH" \
#     --output="$RUNDIR/logs/%x_%A_%a.out" --error="$RUNDIR/logs/%x_%A_%a.err" \
#     "$REPO/mfast-seqs/scripts/run_pipeline_array.sbatch" \
#     "$RUNDIR/samples.tsv" "$RUNDIR/results" "$RUNDIR/config.yaml" \
#     "$REPO/mfast-seqs/scripts"
#
# The 4th argument (this script's own scripts/ directory) is required,
# not optional: SLURM copies a submitted script into a per-job spool
# directory (/var/spool/slurmd/...) and runs THAT copy, so this script
# can't reliably find run_pipeline.sh next to itself via $0/BASH_SOURCE -
# that would resolve to the spool directory, not your git clone.
#
# See ../README.md "How to run" for the full sequence.
#
# If your conda env itself isn't at the default $CONDA_ROOT/envs/<name>
# subpath, override that directly instead - add to the same --export list:
#   ,CONDA_ENV_NAME=my-env             (still under $CONDA_ROOT/envs/)
#   ,CONDA_ENV_PATH=/full/path/to/envs/my-env   (anywhere at all)
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
SCRIPT_DIR="$4"
if [[ -z "$SCRIPT_DIR" || ! -x "$SCRIPT_DIR/run_pipeline.sh" ]]; then
  echo "ERROR: 4th argument must be the mfast-seqs/scripts directory (got '$SCRIPT_DIR')" >&2
  echo "e.g. sbatch ... run_pipeline_array.sbatch samples.tsv results config.yaml \"\$REPO/mfast-seqs/scripts\"" >&2
  exit 1
fi

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
GROUP=$(cut -f3 <<< "$LINE")   # empty string if samples.tsv has no 3rd column

if [[ -n "$GROUP" ]]; then
  SAMPLE_ID="${GROUP}_${SAMPLE}"
else
  SAMPLE_ID="$SAMPLE"
fi

echo "[task $SLURM_ARRAY_TASK_ID] sample=$SAMPLE group=${GROUP:-<none>} -> $SAMPLE_ID fastq=$FASTQ"
mkdir -p "$OUTDIR_ROOT"

"$SCRIPT_DIR/run_pipeline.sh" "$SAMPLE_ID" "$FASTQ" "$OUTDIR_ROOT/$SAMPLE_ID" "$CONFIG"
