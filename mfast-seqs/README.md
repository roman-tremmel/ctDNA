# LINE-1 Aneuploidy Score Pipeline

A from-scratch, FASTQ-to-score reimplementation of the read-processing
side of **mFAST-SeqS** (modified Fast Aneuploidy Screening
Test-Sequencing), the assay used in Verschoor et al. 2023
(*npj Breast Cancer* 9:61, [10.1038/s41523-023-00563-w](https://doi.org/10.1038/s41523-023-00563-w))
and originally described by Belic et al. 2016.

## Background

mFAST-SeqS estimates a tumor-agnostic ctDNA "aneuploidy score" from
plasma cfDNA without any prior knowledge of tumor mutations. A single
primer pair amplifies the ~500,000 LINE-1 (L1) retrotransposon copies
scattered across the genome. Because most solid tumors carry
chromosome-arm-level copy-number gains/losses (aneuploidy), a tumor
sample's plasma will show a skewed distribution of LINE-1 read
counts across chromosome arms relative to a healthy-donor baseline -
even though no single mutation is being tracked.

Wet-lab steps (PCR1 with the LINE-1 primer, PCR2 to add a random
heterogeneity spacer and Illumina indices, sequencing) are **not**
part of this repo; this pipeline starts from the resulting raw FASTQ
files and reproduces the paper's dry-lab analysis:

1. **Trim** the heterogeneity spacer (Fadrosh et al. 2014 dual-indexing
   scheme) and the LINE-1 primer off each read.
2. **Align** the remaining (genomic, non-repetitive-primer) sequence to
   the reference genome.
3. **Count** reads per chromosome arm, restricted to reads that overlap
   an annotated LINE-1 locus, normalized to total library size.
4. **Z-score** each arm's normalized count against a healthy-control
   panel-of-normals (mean/SD per arm) - see "A cohort of healthy
   controls is required" below.
5. **Sum the squared Z-scores** over all included arms (excluding the
   acrocentric short arms 13p/14p/15p/21p/22p and chrY, which carry too
   few LINE-1 elements) to get the genome-wide **aneuploidy score**.
6. Classify against a cutoff (5.0, as used in the paper).

## A cohort of healthy controls is required

The Z-score in step 4 above needs a per-arm mean/SD baseline computed
from a cohort of healthy control samples (the paper: *"a Z-score per
chromosome arm was calculated relative to healthy female controls"*) -
this is a **one-time** cost, not a per-sample one: build the baseline
once from your control cohort (`build_control_baseline.py`, step 5 in
"How to run" below), then score any number of individual case samples
against it independently. Rules of thumb: aim for >=20-30 controls for
a stable per-arm mean/SD, and process them in the same sequencing
run/conditions as your case samples where possible - this "between
sample" comparison is the main way this method is more exposed to
batch effects than `../waldo/`.

## Requirements

System tools: `bwa`, `samtools`, `bedtools`, `cutadapt`. Python 3 with
`numpy`, `pandas`, `scipy`, `pyyaml`. The easiest way to get all of
these without root/admin rights is the top-level `../environment.yml`
conda/mamba environment (covers both pipelines) - see "How to run"
below. `requirements.txt` here covers just the Python side if you
already have the system tools some other way.

## Testing

`tests/test_pipeline_synthetic.py` builds a small synthetic genome and
LINE-1-like loci in memory, simulates reads for several euploid
"controls" and one sample with an engineered chromosome-arm gain/loss,
and runs the full pipeline (cutadapt -> bwa -> samtools -> bedtools ->
scoring) end to end to confirm the code is wired together correctly:

```bash
python3 tests/test_pipeline_synthetic.py
```

This validates pipeline *mechanics*, not the real hg38/LINE-1 biology —
it's a toy genome, not a clinical validation. It's self-contained
(writes to a temp directory) and needs no setup from "How to run" below.

## Key methodological differences from the published assay

- The LINE-1 forward-primer core (`ACACAGGGAGGGGAACAT`) is extracted
  directly from Verschoor et al. 2023 Supplementary Table 1 (see
  `resources/primers_supplementary_table1.tsv`); the index/adapter
  primer portions are also transcribed there for reference, though this
  pipeline only needs the LINE-1-specific core for trimming.
- Chromosome-arm boundaries are hardcoded approximate GRCh38 centromere
  coordinates - verify against UCSC's `cytoBand`/`gap` tables before any
  diagnostic use.
- The paper used SPSS for downstream statistics (Z-scores, Cox models
  for survival); this repo only reproduces the FASTQ -> aneuploidy-score
  computation, not the clinical/survival analysis in the paper.

## How to run

**Nothing below ever writes inside the git clone.** Everything
generated (environment, downloaded resources, your `config.yaml`,
results, logs) goes into an external run directory, so `git pull`
always stays conflict-free no matter what you've run. Set these two
variables once per shell session (adjust to your actual paths):

```bash
export REPO=/path/to/ctDNA        # your git clone (top level, not mfast-seqs/)
export RUNDIR=/path/to/your/rundir   # anywhere OUTSIDE the git clone
mkdir -p "$RUNDIR"/{resources,results,logs}
```

**1. Environment (once per machine, no root needed):**

```bash
export CONDA_ROOT="$HOME/miniforge3"   # note this value - step 4's sbatch calls need it too

curl -fsSL -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p "$CONDA_ROOT"
source "$CONDA_ROOT/etc/profile.d/conda.sh"
mamba env create -f "$REPO/environment.yml"
conda activate ctdna-aneuploidy

# Find the env's ACTUAL full path and note it too - you need this exact
# value for every sbatch call in step 4. Don't assume it's under
# $CONDA_ROOT/envs/: env managers (especially a standalone `mamba`
# binary, as opposed to one bundled with conda/miniconda) often default
# to a DIFFERENT root prefix, e.g. ~/.local/share/mamba/envs/.
conda env list | grep ctdna-aneuploidy
export CONDA_ENV_PATH="$(conda env list | awk '/ctdna-aneuploidy/ {print $NF}')"
echo "CONDA_ENV_PATH=$CONDA_ENV_PATH"   # sanity-check this isn't empty
```

If Miniforge is already installed on your cluster, set `CONDA_ROOT` to
wherever it actually lives instead
(check with `conda info --base`, or `ls -d ~/miniforge3 ~/miniconda3
~/mambaforge ~/anaconda3 2>/dev/null` if `conda` isn't on your PATH by
default) - `$HOME/miniforge3` is only this guide's own install location
from the command above, not a requirement.

If your cluster has no outbound internet access at all, build this
environment on a machine that does, pack it with `conda-pack`, and
transfer the resulting tarball - see the top-level `../README.md`'s
"Running on an air-gapped / offline HPC cluster" section for the exact
commands.

**2. Reference resources (once per machine, into `$RUNDIR`):**

```bash
# LINE-1 element BED - downloads UCSC RepeatMasker (~500MB), needs internet
"$REPO/mfast-seqs/resources/build_line1_bed.sh" "$RUNDIR/resources/line1_elements.hg38.bed"
```

For the genome: reuse an existing site-wide BWA index if your cluster
has one (e.g. an iGenomes GATK GRCh38 bundle) - nothing to build or
copy, just reference its path in step 3. Otherwise:

```bash
cp /path/to/genome.fa "$RUNDIR/resources/genome.fa"
bwa index "$RUNDIR/resources/genome.fa"
```

`$REPO/mfast-seqs/resources/chrom_arms.hg38.bed` already ships with the
repo (static, deterministic) - nothing to build for it.

**3. Your `config.yaml` (once, outside the git clone):**

```bash
cp "$REPO/mfast-seqs/config/config.example.yaml" "$RUNDIR/config.yaml"
```

Edit `$RUNDIR/config.yaml` (never the tracked `config.example.yaml`):
- `reference.fasta` -> your BWA index prefix (site-wide, or the one from step 2)
- `reference.chrom_arms_bed` -> `$REPO/mfast-seqs/resources/chrom_arms.hg38.bed`
- `reference.line1_bed` -> `$RUNDIR/resources/line1_elements.hg38.bed`
- everything else has sensible defaults already filled in

**4. Run every sample (controls and cases alike):**

```bash
"$REPO/mfast-seqs/scripts/make_sample_sheet.sh" /path/to/your/fastq_dir > "$RUNDIR/samples.tsv"
sbatch --array=1-$(wc -l < "$RUNDIR/samples.tsv") \
  --export=ALL,CONDA_ROOT="$CONDA_ROOT",CONDA_ENV_PATH="$CONDA_ENV_PATH" \
  --output="$RUNDIR/logs/%x_%A_%a.out" --error="$RUNDIR/logs/%x_%A_%a.err" \
  "$REPO/mfast-seqs/scripts/run_pipeline_array.sbatch" \
  "$RUNDIR/samples.tsv" "$RUNDIR/results" "$RUNDIR/config.yaml" \
  "$REPO/mfast-seqs/scripts"
```

(Or without SLURM, one sample at a time:
`"$REPO/mfast-seqs/scripts/run_pipeline.sh" <sample> <fastq> "$RUNDIR/results/<sample>" "$RUNDIR/config.yaml"`.)

**Neither `--export=...` nor the final `"$REPO/mfast-seqs/scripts"`
argument is optional** - both are load-bearing, not just-in-case:

- `CONDA_ROOT`/`CONDA_ENV_PATH`: SLURM does not reliably hand a
  submitting shell's exported variables to the job on every cluster's
  configuration, so without passing them explicitly via `--export`,
  `run_pipeline_array.sbatch` falls back to hardcoded defaults
  (`$HOME/miniforge3` and `$CONDA_ROOT/envs/ctdna-aneuploidy`) and
  fails - either with `.../conda.sh: No such file or directory` (wrong
  `CONDA_ROOT`) or `libmamba ... Cannot activate, prefix does not exist`
  (env isn't actually under `$CONDA_ROOT/envs/`, e.g. because a
  standalone `mamba` binary, separate from a `conda`/`miniconda3`
  install, defaults to its own root prefix at `~/.local/share/mamba`).
  Step 1 above has you determine and export the real `$CONDA_ENV_PATH`
  for exactly this reason - use that value here, don't skip it because
  the flag "looks like" it should be optional.
- The `scripts` argument: SLURM copies a submitted script into a
  per-job spool directory (`/var/spool/slurmd/...`) and runs that copy,
  so `run_pipeline_array.sbatch` can't reliably find `run_pipeline.sh`
  next to itself at runtime - it fails with
  `.../run_pipeline.sh: No such file or directory` without it.

This produces `$RUNDIR/results/<sample>/<sample>.arm_counts.tsv` for
every sample.

**5. Build the healthy-control baseline (once, from control samples only):**

```bash
python3 "$REPO/mfast-seqs/scripts/build_control_baseline.py" \
  --counts "$RUNDIR"/results/control_*/control_*.arm_counts.tsv \
  --out "$RUNDIR/results/baseline.tsv"
```

**6. Score each case sample against that baseline:**

```bash
python3 "$REPO/mfast-seqs/scripts/compute_aneuploidy_score.py" \
  --counts "$RUNDIR/results/<sample>/<sample>.arm_counts.tsv" \
  --baseline "$RUNDIR/results/baseline.tsv" \
  --cutoff 5.0 \
  --out "$RUNDIR/results/<sample>.aneuploidy.tsv"
```
