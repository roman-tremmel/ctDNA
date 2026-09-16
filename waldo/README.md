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

## Setup

Same reference genome as `../mfast-seqs/` works here; only the
window/cluster resources are specific to this pipeline.

```bash
bwa index resources/genome.fa
python3 resources/build_chrom_arms_bed.py -o resources/chrom_arms.hg38.bed
python3 resources/build_genome_windows.py -o resources/windows.500kb.hg38.bed
resources/build_line1_bed.sh resources/line1_elements.hg38.bed
```

`config/config.yaml`'s `primer.forward_seq` is pre-filled with the
LINE-1-specific core extracted from Verschoor et al. 2023 Supplementary
Table 1 (`resources/primers_supplementary_table1.tsv`) - the same assay
family WALDO's own paper cites (Kinde et al. 2012 FAST-SeqS), but not
confirmed identical to WALDO's own primer, which wasn't available while
building this. Replace it if you have WALDO's actual primer sequence.

See `../mfast-seqs/README.md`'s "Running on an air-gapped / offline HPC
cluster" section - the same applies here: only the genome/LINE-1-BED
setup steps need internet, the pipeline itself runs entirely offline
once those resources exist.

## Running

Per-sample counting (same pattern as mfast-seqs/):

```bash
scripts/run_pipeline.sh <sample_name> <sample.fastq.gz> results/<sample_name> config/config.yaml
```

Build cluster membership from a panel of euploid reference samples (the
paper found 7 sufficient; this is the SAME panel used for threshold
calibration below, matching the paper's own reuse of its 677-WBC panel
for both purposes):

```bash
python3 scripts/build_window_clusters.py \
  --counts results/control_*/control_*.window_counts.tsv \
  --mean-p-threshold 0.05 --var-p-threshold 0.05 \
  --out results/clusters.tsv
```

Calibrate per-arm significance thresholds from the same reference panel:

```bash
python3 scripts/calibrate_threshold.py \
  --counts results/control_*/control_*.window_counts.tsv \
  --clusters results/clusters.tsv \
  --windows-bed resources/windows.500kb.hg38.bed \
  --arms-bed resources/chrom_arms.hg38.bed \
  --excluded-arms chr13p chr14p chr15p chr21p chr22p chrYp chrYq \
  --margin 4.0 \
  --out results/thresholds.tsv
```

Call per-arm gains/losses for a case sample:

```bash
python3 scripts/call_aneuploidy.py \
  --counts results/<sample_name>/<sample_name>.window_counts.tsv \
  --clusters results/clusters.tsv \
  --windows-bed resources/windows.500kb.hg38.bed \
  --arms-bed resources/chrom_arms.hg38.bed \
  --excluded-arms chr13p chr14p chr15p chr21p chr22p chrYp chrYq \
  --threshold-table results/thresholds.tsv \
  --out results/<sample_name>.arm_scores.tsv
```

Optional: train and apply the SVM genome-wide classifier once you have
labeled samples (per the paper, only meaningful for samples where step
above found no single significant arm):

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
this scale (see "Known limitations" above).
