# WALDO-Style Within-Sample Aneuploidy Pipeline

A from-scratch reimplementation of WALDO (Within-Sample AneupLoidy
DetectiOn), described in Douville et al. 2018, *PNAS* 115(8):1871-1876
([10.1073/pnas.1717846115](https://doi.org/10.1073/pnas.1717846115)),
following the algorithm as specified in that paper's **SI Appendix, SI
Materials and Methods**.

This is a **separate, independent pipeline** from `../mfast-seqs/` -
it shares the same underlying wet-lab chemistry (FAST-SeqS: a single
primer pair amplifying LINE-1 elements genome-wide, Kinde et al. 2012)
but implements a fundamentally different analysis strategy. See "How
this differs from mfast-seqs/" below. The SI Appendix's own "Comparison
with a prior technique" section describes essentially the same
arm-level, between-sample Z-score approach as `mfast-seqs/` as the
baseline WALDO was benchmarked against.

## Fidelity to the published method

This reimplements the algorithm **as written in the SI Materials and
Methods**, not a loose approximation:

- **Per-window clusters, not k-means/hierarchical clusters.** For every
  500kb window *i*, its cluster C_i is every OTHER window *i'* on a
  DIFFERENT chromosome whose reference-panel read-count vector is not
  significantly different from *i*'s in either mean (paired t-test,
  p > 0.05) or variance (F-test, p > 0.05). This is a personalized
  neighbor list per window (the paper reports ~200 members per cluster
  on real data), not a partition into a fixed number of mutually
  exclusive groups. Restricting candidates to other chromosomes is
  deliberate in the source text and is what makes the resulting test
  structurally unable to be blind to any single arm's own aneuploidy
  (an early hierarchical-clustering approximation we tried before
  reading the SI didn't have this guarantee - see "History" below).
- **Within-sample cluster statistics.** For a given test sample, each
  window's cluster mean/variance (mu_i, sigma_i^2) are estimated from
  *that same sample's own scaled read values* at its cluster's member
  windows, with iterative outlier trimming (members with two-sided
  p < 0.01 under the running Normal(mu_i, sigma_i^2) estimate are
  dropped, repeated to convergence). Only cluster *membership* comes
  from the reference panel; the actual numbers being tested are the
  sample's own.
- **Sum-of-normals arm test.** Since a chromosome arm's windows'
  scaled reads sum to a random variable that's (approximately) the sum
  of independent normals, `sum(R_i) ~ N(sum(mu_i), sum(sigma_i^2))`,
  giving a closed-form per-arm Z-score.
- **Empirically-calibrated significance threshold**
  (`scripts/calibrate_threshold.py`), not a fixed Z cutoff: the paper
  sets an arm's gain/loss threshold at 4 SD beyond the most extreme
  Z-score observed when this same test is run on a large euploid
  reference panel.

## An acknowledged ambiguity

The SI text describes mu_i/sigma_i^2 as estimated "from the genomic
intervals in each ... cluster" without stating explicitly whether that
estimation uses the reference panel or the sample being tested. We read
it as within-sample (computed from the test sample's own values),
because SI Appendix Figure S8B-C is described as illustrating cluster
normality separately for "a normal WBC sample" and "an aneuploid
Primary Tumor sample" - i.e., as a per-sample property - and because
within-sample estimation is what makes the method's own name
("Within-Sample AneupLoidy DetectiOn") and its stated batch-effect
motivation coherent. If you have access to the paper's actual source
code or a clearer reading of the SI, this is the place to correct.

## How this differs from mfast-seqs/

| | `mfast-seqs/` | `waldo/` (this pipeline) |
|---|---|---|
| Unit of analysis | Whole chromosome arm | 500kb windows, personalized per-window clusters |
| Normalization | *Between-sample*: arm's read fraction vs. a panel-of-normals mean/SD | *Within-sample*: a window's reads vs. the SAME sample's own reads in its own cluster |
| Cluster membership source | N/A | Learned once from a reference panel (paired t-test + F-test equivalence, cross-chromosome only) |
| Significance threshold | Fixed literature cutoff (5.0) | Empirically calibrated from a reference panel's own Z-score range + margin |
| Extra layer | None | Optional SVM classifier over per-arm Z-scores, for samples with many weak (individually non-significant) arm perturbations |

If you just need the simple, well-validated arm-level score from
Verschoor et al. 2023 / Belic et al. 2016, use `../mfast-seqs/`. This
pipeline trades simplicity for (claimed) better batch-effect robustness
and sensitivity at low tumor fractions.

## Pipeline

1. **Trim + align** (identical wet-lab chemistry to mfast-seqs/, scripts
   duplicated here to keep the two pipelines independent):
   `trim_line1_reads.sh`, `align_reads.sh`. The SI specifies Bowtie2
   against GRCh37 with degenerate-base molecular barcoding (UMIs) for
   exact template counting; this pipeline uses bwa mem against hg38
   without UMI deduplication (documented simplification, see "Not
   implemented" below) - align/count logic just needs the arm/window
   BEDs regenerated for whichever build you target.
2. **Count LINE-1 reads per 500kb window**
   (`scripts/count_reads_by_window.py`), restricted to reads overlapping
   an annotated LINE-1 locus.
3. **Learn per-window cluster membership** from a euploid reference
   panel (`scripts/build_window_clusters.py`, paired t-test + F-test
   equivalence, cross-chromosome only). The paper found 7 reference
   samples sufficient; using more doesn't hurt but increases the
   tests' statistical power, which *shrinks* clusters (fewer windows
   look "indistinguishable" from a well-powered test) - this is
   expected, not a bug.
4. **Within-sample chromosome-arm test**
   (`scripts/call_aneuploidy.py`): computes each window's within-sample
   cluster mean/variance (with outlier trimming) and the resulting
   per-arm Z-score.
5. **Calibrate a significance threshold** (`scripts/calibrate_threshold.py`):
   runs step 4 on a euploid reference panel to get each arm's observed
   Z-score range, then sets `gain_threshold = max + margin`,
   `loss_threshold = min - margin` (default margin 4.0, per the paper).
6. **Optional genome-wide SVM classifier**
   (`scripts/train_svm_classifier.py`, `scripts/classify_sample.py`):
   the paper found that at very low neoplastic fractions, no single arm
   reaches its calibrated threshold but many arms show small, consistent
   deviations; an SVM over the vector of per-arm Z-scores (the paper
   uses 39 - all arms except the 5 acrocentric short arms) can pick up
   on that combined pattern. Per the paper, only run this on samples
   where no single arm was already called significant in step 5. This
   ships as a trainable component (no pretrained model, since the
   paper's own model/weights are not published) - train it on your own
   labeled samples (real and/or synthetic).

## Not implemented

- **Molecular barcoding / UMI deduplication.** The paper uses degenerate
  primer bases as UMIs so each original template molecule is counted
  once regardless of PCR duplication. This pipeline counts aligned
  reads directly (as `mfast-seqs/` does), which is simpler but more
  exposed to PCR duplication bias. Adding UMI support would mean
  extending `count_reads_by_window.py` to parse a UMI from the read
  (or its original untrimmed prefix) and deduplicate before counting.
- **Amplicon-size-based reference-sample matching.** The paper matches
  each test sample to reference samples with a similar input-DNA
  fragment-size distribution (computed from amplicon sizes) before
  clustering, to control for fragment-size-dependent amplification
  bias - the same simplification `mfast-seqs/` documents. Doing this
  properly needs per-read fragment/insert-size information our
  single-end counting doesn't currently extract.
- **SNP-based allelic imbalance, sample-identity fingerprinting, somatic
  mutation calling and MSI detection.** All four use the ~26,220 common
  polymorphisms within the amplified LINE-1 loci (24,720 SNPs + 1,500
  indels per the SI). None are implemented here - they need a validated
  SNP-within-LINE-1 site list for your reference build and molecular
  barcoding for the mutation-calling steps, neither of which were
  available while building this.
- **SVM read-depth calibration.** The paper corrects raw SVM scores for
  a systematic read-depth-dependent inflation (low-depth samples score
  artificially high), fit from downsampling experiments on their own
  63-sample WBC cohort. `train_svm_classifier.py`/`classify_sample.py`
  don't apply this correction - if your samples vary a lot in read
  depth, be aware low-depth samples may score higher than warranted.

## Known limitations and caveats

- **Real-scale reference panels needed for calibration.** With few
  reference samples, `calibrate_threshold.py`'s "observed max/min +
  margin" is a poor estimate of the true tail (the paper used 677
  normal WBC samples). `tests/test_waldo_synthetic.py` demonstrates
  this directly at toy scale and checks Z-score direction/magnitude
  rather than the resulting gain/loss call for exactly this reason.
  Use as large a reference panel as you reasonably can.
- **Cluster size scales with genome coverage.** The paper's ~200
  windows/cluster comes from having 4,361 genome-wide windows to search
  among; a small reference genome or a heavily filtered window set will
  yield much smaller (or empty) clusters, and windows with too few
  surviving cluster members after outlier trimming are excluded and
  reported as "no_data" for that arm.
- **The mu_i/sigma_i^2 "within-sample vs. reference-panel" ambiguity**
  noted above - flag if you find contrary evidence.
- Chromosome-arm boundaries are hardcoded approximate GRCh38 centromere
  coordinates (same caveat as mfast-seqs/) - verify against UCSC before
  diagnostic use. The paper itself used GRCh37.

## Testing

```bash
python3 tests/test_waldo_synthetic.py
```

Builds a small synthetic genome with windows/clusters/arms, simulates a
euploid reference panel plus one sample with an engineered arm gain and
an arm loss (with realistic Poisson read-sampling noise and per-window
amplification-efficiency differences that cut across chromosomes, so
cross-chromosome clustering has real structure to find), and runs the
full pipeline end to end - trim, align, count, cluster, calibrate,
score. Validates pipeline *mechanics* on a toy-scale genome; not a
validation of the real biology, of numerical agreement with the
published WALDO tool, or of the calibrated threshold being reliable at
this scale (see "Known limitations" above). It's self-contained (writes
to a temp directory) and needs no setup from "How to run" below.

## How to run

**Nothing below ever writes inside the git clone.** Everything
generated (environment, downloaded resources, your `config.yaml`,
results, logs) goes into an external run directory, so `git pull`
always stays conflict-free no matter what you've run. Set these two
variables once per shell session (adjust to your actual paths):

```bash
export REPO=/path/to/ctDNA        # your git clone (top level, not waldo/)
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
commands. (Same environment as `../mfast-seqs/` - one env covers both
pipelines.)

**2. Reference resources (once per machine, into `$RUNDIR`):**

```bash
# LINE-1 element BED - downloads UCSC RepeatMasker (~500MB), needs internet
"$REPO/waldo/resources/build_line1_bed.sh" "$RUNDIR/resources/line1_elements.hg38.bed"
```

For the genome: reuse an existing site-wide BWA index if your cluster
has one (e.g. an iGenomes GATK GRCh38 bundle) - nothing to build or
copy, just reference its path in step 3. Otherwise:

```bash
cp /path/to/genome.fa "$RUNDIR/resources/genome.fa"
bwa index "$RUNDIR/resources/genome.fa"
```

`$REPO/waldo/resources/chrom_arms.hg38.bed` and
`$REPO/waldo/resources/windows.500kb.hg38.bed` already ship with the
repo (static, deterministic) - nothing to build for either.

**3. Your `config.yaml` (once, outside the git clone):**

```bash
cp "$REPO/waldo/config/config.example.yaml" "$RUNDIR/config.yaml"
```

Edit `$RUNDIR/config.yaml` (never the tracked `config.example.yaml`):
- `reference.fasta` -> your BWA index prefix (site-wide, or the one from step 2)
- `reference.chrom_arms_bed` -> `$REPO/waldo/resources/chrom_arms.hg38.bed`
- `reference.windows_bed` -> `$REPO/waldo/resources/windows.500kb.hg38.bed`
- `reference.line1_bed` -> `$RUNDIR/resources/line1_elements.hg38.bed`
- everything else has sensible defaults already filled in

**4. Run every sample (reference/control and case samples alike):**

```bash
"$REPO/waldo/scripts/make_sample_sheet.sh" /path/to/your/fastq_dir > "$RUNDIR/samples.tsv"
sbatch --array=1-$(wc -l < "$RUNDIR/samples.tsv") \
  --export=ALL,CONDA_ROOT="$CONDA_ROOT",CONDA_ENV_PATH="$CONDA_ENV_PATH" \
  --output="$RUNDIR/logs/%x_%A_%a.out" --error="$RUNDIR/logs/%x_%A_%a.err" \
  "$REPO/waldo/scripts/run_pipeline_array.sbatch" \
  "$RUNDIR/samples.tsv" "$RUNDIR/results" "$RUNDIR/config.yaml" \
  "$REPO/waldo/scripts"
```

(Or without SLURM, one sample at a time:
`"$REPO/waldo/scripts/run_pipeline.sh" <sample> <fastq> "$RUNDIR/results/<sample>" "$RUNDIR/config.yaml"`.)

**Neither `--export=...` nor the final `"$REPO/waldo/scripts"`
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

This produces `$RUNDIR/results/<sample>/<sample>.window_counts.tsv` for
every sample.

**5. Learn cluster membership from your euploid reference samples**
(7 sufficient per the paper; the SAME panel is reused for calibration
in step 6, matching the paper's own approach):

Nothing in `run_pipeline.sh`/`make_sample_sheet.sh` knows which of your
sample IDs are euploid reference samples vs. cases - maintain your own
plain-text list (names must match what you used as `<sample_name>` in
step 4) and build the `--counts` argument from it:

```bash
cat > "$RUNDIR/controls.txt" <<'EOF'
V91-02
V91-07
V91-11
EOF

python3 "$REPO/waldo/scripts/build_window_clusters.py" \
  --counts $("$REPO/waldo/scripts/counts_for_samples.sh" \
               "$RUNDIR/results" .window_counts.tsv "$RUNDIR/controls.txt") \
  --mean-p-threshold 0.05 --var-p-threshold 0.05 \
  --out "$RUNDIR/results/clusters.tsv"
```

(If your sample IDs share a consistent, greppable prefix instead, a
plain shell glob works just as well - `counts_for_samples.sh` is for
when they don't.)

**6. Calibrate per-arm significance thresholds from the same panel:**

```bash
python3 "$REPO/waldo/scripts/calibrate_threshold.py" \
  --counts $("$REPO/waldo/scripts/counts_for_samples.sh" \
               "$RUNDIR/results" .window_counts.tsv "$RUNDIR/controls.txt") \
  --clusters "$RUNDIR/results/clusters.tsv" \
  --windows-bed "$REPO/waldo/resources/windows.500kb.hg38.bed" \
  --arms-bed "$REPO/waldo/resources/chrom_arms.hg38.bed" \
  --excluded-arms chr13p chr14p chr15p chr21p chr22p chrYp chrYq \
  --margin 4.0 \
  --out "$RUNDIR/results/thresholds.tsv"
```

**7. Call per-arm gains/losses for a case sample:**

```bash
python3 "$REPO/waldo/scripts/call_aneuploidy.py" \
  --counts "$RUNDIR/results/<sample>/<sample>.window_counts.tsv" \
  --clusters "$RUNDIR/results/clusters.tsv" \
  --windows-bed "$REPO/waldo/resources/windows.500kb.hg38.bed" \
  --arms-bed "$REPO/waldo/resources/chrom_arms.hg38.bed" \
  --excluded-arms chr13p chr14p chr15p chr21p chr22p chrYp chrYq \
  --threshold-table "$RUNDIR/results/thresholds.tsv" \
  --out "$RUNDIR/results/<sample>.arm_scores.tsv"
```

**8. Optional: train/apply the SVM genome-wide classifier**, once you
have labeled samples (per the paper, only meaningful for samples where
step 7 found no single significant arm):

```bash
python3 "$REPO/waldo/scripts/build_feature_matrix.py" \
  --scores "$RUNDIR/results/sampleA.arm_scores.tsv" "$RUNDIR/results/sampleB.arm_scores.tsv" \
  --names sampleA sampleB \
  --out "$RUNDIR/results/feature_matrix.tsv"

python3 "$REPO/waldo/scripts/train_svm_classifier.py" \
  --features "$RUNDIR/results/feature_matrix.tsv" --labels "$RUNDIR/labels.tsv" \
  --out "$RUNDIR/results/model.joblib"

python3 "$REPO/waldo/scripts/classify_sample.py" \
  --scores "$RUNDIR/results/<sample>.arm_scores.tsv" --model "$RUNDIR/results/model.joblib"
```
