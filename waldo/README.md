# WALDO-Style Within-Sample Aneuploidy Pipeline

An original, from-scratch reimplementation of the *ideas* behind WALDO
(Within-Sample AneupLoidy DetectiOn), described in Douville et al. 2018,
*PNAS* 115(8):1871-1876 ([10.1073/pnas.1717846115](https://doi.org/10.1073/pnas.1717846115)).

This is a **separate, independent pipeline** from `../mfast-seqs/` -
it shares the same underlying wet-lab chemistry (FAST-SeqS: a single
primer pair amplifying LINE-1 elements genome-wide, Kinde et al. 2012)
but implements a different, more elaborate analysis strategy. See
"How this differs from mfast-seqs/" below.

## IMPORTANT: this is a good-faith reconstruction, not the published algorithm

The PNAS paper's exact statistical model, clustering algorithm, and SVM
training procedure are described only at a conceptual level in the main
text; the precise equations are in its *SI Appendix*, which was not
available when this was built. This pipeline implements an **original
statistical test that captures the paper's stated design principles**
(within-sample, cluster-based comparison instead of between-sample
normalized-count comparison), but it is **not a byte-exact reproduction**
of WALDO and will not numerically match the published tool. Treat it as
a methodologically-inspired starting point to adapt and validate against
your own data, not a validated clinical assay. See "Known limitations
and caveats" below before trusting its output.

## How this differs from mfast-seqs/

Both pipelines start from the same amplicon sequencing data (LINE-1
reads) and end with a genome-wide aneuploidy signal, but the statistics
differ fundamentally:

| | `mfast-seqs/` | `waldo/` (this pipeline) |
|---|---|---|
| Unit of analysis | Whole chromosome arm | 500kb windows, grouped into clusters |
| Normalization | *Between-sample*: arm's read fraction vs. a panel-of-normals mean/SD | *Within-sample*: a window's reads vs. the SAME sample's own reads in windows that historically track with it |
| Why it matters | Simple, but sensitive to batch effects (sequencing run, reagent lot, input amount) shifting the whole panel-of-normals comparison | More robust to batch effects, because the comparison never needs read counts to be numerically comparable across runs/samples |
| Output | One genome-wide score (sum of squared arm Z-scores) | Per-arm gain/loss calls + an optional SVM genome-wide classifier for low-tumor-fraction samples |
| Complexity | Low - a handful of scripts, one clear formula | High - clustering, per-window baselines, more assumptions and tunable parameters |

If you just need the simple, well-validated arm-level score from
Verschoor et al. 2023 / Belic et al. 2016, use `../mfast-seqs/`. This
pipeline is for exploring whether the finer-grained, within-sample
approach gives you better sensitivity (e.g. for low ctDNA-fraction
plasma samples), at the cost of more moving parts.

## Pipeline

1. **Trim + align** (identical wet-lab chemistry to mfast-seqs/, scripts
   duplicated here to keep the two pipelines independent):
   `trim_line1_reads.sh`, `align_reads.sh`.
2. **Count LINE-1 reads per 500kb window**
   (`scripts/count_reads_by_window.py`), restricted to reads overlapping
   an annotated LINE-1 locus, same on-target filtering as mfast-seqs/.
3. **Learn window clusters + a within-sample baseline** from a panel of
   euploid reference samples (`scripts/build_window_clusters.py`):
   windows whose normalized read fractions correlate across the
   reference panel are grouped into clusters (hierarchical clustering on
   1 - Pearson correlation). Within each reference sample, each window's
   share of its own cluster's total reads is computed; the mean/SD of
   that share across the reference panel is the baseline.
4. **Within-sample chromosome-arm test** (`scripts/call_aneuploidy.py`):
   for a test sample, a window's *expected* reads are predicted from the
   window's reference share-of-cluster ratio times the *test sample's
   own* cluster total (not an absolute reference count) - so a uniform
   depth/batch shift in the test sample cancels out. Expected
   means/variances are summed over a chromosome arm's windows (sum of
   normals is normal) to get a per-arm Z-score and gain/loss call.
5. **Optional genome-wide SVM classifier**
   (`scripts/train_svm_classifier.py`, `scripts/classify_sample.py`):
   the paper found that at very low neoplastic fractions, no single arm
   reaches significance but many arms show small consistent deviations;
   an SVM over the vector of per-arm Z-scores can pick up on that
   combined pattern. This ships as a trainable component (no pretrained
   model, since the paper's is not published) - train it on your own
   labeled samples (real and/or synthetic).

## Not implemented

The paper describes several additional capabilities built on ~26,220
common SNPs within the amplified LINE-1 loci, none of which are
implemented here (they need resources - a validated SNP panel within
hg38 LINE-1 loci, molecular barcoding info for exact template counting -
that were not available while building this):

- Per-arm allelic imbalance (can catch copy-neutral LOH that read depth
  alone cannot).
- Sample-identity fingerprinting / relatedness checks (matching a plasma
  sample to its primary tumor, or detecting replicate/mislabeled samples).
- Somatic mutation calling and microsatellite-instability detection
  within the LINE amplicons (requires a matched normal sample).

If you need these, you'd have to source or build an hg38 SNP-within-LINE-1
site list and extend `count_reads_by_window.py`-style logic into a
pileup/genotyping step (e.g. with `pysam`).

## Setup

Same reference genome as `../mfast-seqs/` works here; only the
window/cluster resources are specific to this pipeline.

```bash
bwa index resources/genome.fa
python3 resources/build_chrom_arms_bed.py -o resources/chrom_arms.hg38.bed
python3 resources/build_genome_windows.py -o resources/windows.500kb.hg38.bed
resources/build_line1_bed.sh resources/line1_elements.hg38.bed
```

Fill in your real LINE-1 primer sequence in `config/config.yaml`
(`primer.forward_seq`).

See `../mfast-seqs/README.md`'s "Running on an air-gapped / offline HPC
cluster" section - the same applies here: only the genome/LINE-1-BED
setup steps need internet, the pipeline itself runs entirely offline
once those resources exist.

## Running

Per-sample counting (same pattern as mfast-seqs/):

```bash
scripts/run_pipeline.sh <sample_name> <sample.fastq.gz> results/<sample_name> config/config.yaml
```

Build clusters + baseline from a panel of euploid reference samples
(more references = more reliable correlation-based clustering and
per-window SD estimates - the paper's own within-sample premise still
needs *some* reference panel to learn which windows track together):

```bash
python3 scripts/build_window_clusters.py \
  --counts results/control_*/control_*.window_counts.tsv \
  --n-clusters 40 \
  --arms-bed resources/chrom_arms.hg38.bed \
  --out-clusters results/clusters.tsv \
  --out-baseline results/window_baseline.tsv
```

Call per-arm gains/losses for a case sample:

```bash
python3 scripts/call_aneuploidy.py \
  --counts results/<sample_name>/<sample_name>.window_counts.tsv \
  --baseline results/window_baseline.tsv \
  --windows-bed resources/windows.500kb.hg38.bed \
  --arms-bed resources/chrom_arms.hg38.bed \
  --excluded-arms chr13p chr14p chr15p chr21p chr22p chrYp chrYq \
  --z-cutoff 3.0 \
  --out results/<sample_name>.arm_scores.tsv
```

Optional: train and apply the SVM genome-wide classifier once you have
labeled samples (arm-score TSVs from `call_aneuploidy.py` + known
euploid/aneuploid status):

```bash
python3 scripts/build_feature_matrix.py \
  --scores results/sampleA.arm_scores.tsv results/sampleB.arm_scores.tsv ... \
  --names sampleA sampleB ... \
  --out results/feature_matrix.tsv

python3 scripts/train_svm_classifier.py \
  --features results/feature_matrix.tsv --labels labels.tsv \
  --out results/model.joblib

python3 scripts/classify_sample.py \
  --scores results/<sample_name>.arm_scores.tsv --model results/model.joblib
```

## Known limitations and caveats

- **Per-window variance is treated as independent** when summed to the
  arm level. In reality, windows sharing a cluster are compositionally
  linked (their within-sample ratios sum to 1), so this is an
  approximation, not an exact propagation of a joint distribution.
- **A cluster that ends up fully contained within one chromosome arm is
  structurally blind to that arm's aneuploidy status** (its
  expected-vs-actual reads are forced to sum to zero deviation by
  construction). `build_window_clusters.py --arms-bed ...` warns when
  this happens; prefer fewer, larger clusters or a bigger reference
  panel if you see this warning often.
- **Small reference panels give unreliable per-window SD estimates**,
  which can inflate Z-scores and false-positive rates well beyond what
  the nominal Z-score cutoff would suggest (seen directly in
  `tests/test_waldo_synthetic.py`, where a purely euploid held-out
  sample sometimes crosses `arm_z_cutoff`). **Don't trust the default
  cutoff at face value - calibrate `arm_z_cutoff` empirically** against
  your own held-out euploid samples (leave-one-out from your reference
  panel) before using it for calls.
- The paper matches each test sample to reference samples with similar
  input-DNA fragment-size distributions to control for
  fragment-size-dependent amplification bias; this pipeline does not
  implement that matching/correction step.
- Chromosome-arm boundaries are hardcoded approximate GRCh38 centromere
  coordinates (same caveat as mfast-seqs/) - verify against UCSC before
  diagnostic use.

## Testing

```bash
python3 tests/test_waldo_synthetic.py
```

Builds a small synthetic genome with windows/clusters/arms, simulates a
euploid reference panel plus one sample with an engineered arm gain and
an arm loss (both with realistic Poisson read-sampling noise and
per-window amplification-efficiency differences that cut across
chromosome arms, so the test genuinely exercises cross-arm clustering),
and runs the full pipeline end to end. This validates pipeline
*mechanics*; it is not a validation of the real biology, of numerical
agreement with the published WALDO tool, or of the default
`arm_z_cutoff` being well-calibrated for real data (see "Known
limitations" above - the test itself tolerates a couple of
false-positive calls on the euploid holdout for exactly this reason).
