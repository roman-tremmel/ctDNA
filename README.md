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
   panel-of-normals (mean/SD per arm).
5. **Sum the squared Z-scores** over all included arms (excluding the
   acrocentric short arms 13p/14p/15p/21p/22p and chrY, which carry too
   few LINE-1 elements) to get the genome-wide **aneuploidy score**.
6. Classify against a cutoff (5.0, as used in the paper).

## Requirements

System tools (install via your package manager or conda/mamba):
`bwa`, `samtools`, `bedtools`, `cutadapt`.

Python 3 packages: `numpy`, `pandas`, `scipy`, `pyyaml`
(`pip install -r requirements.txt`).

## Setup

1. Get a reference genome and bwa-index it:
   ```bash
   bwa index resources/genome.fa
   ```
2. Build the chromosome-arm BED (GRCh38 coordinates included, regenerate
   with `resources/build_chrom_arms_bed.py` if needed, or adapt for
   another genome build):
   ```bash
   python3 resources/build_chrom_arms_bed.py -o resources/chrom_arms.hg38.bed
   ```
3. Build the LINE-1 element BED from the UCSC RepeatMasker track (one-time,
   ~500 MB download):
   ```bash
   resources/build_line1_bed.sh resources/line1_elements.hg38.bed
   ```
4. Copy `config/config.yaml`, point it at your reference/BED files, and
   fill in your assay's real LINE-1 primer sequence
   (`primer.forward_seq`) — this is assay-specific and intentionally
   left as a placeholder here.

## Running

Per-sample: raw FASTQ -> trimmed -> aligned -> per-arm counts:

```bash
scripts/run_pipeline.sh <sample_name> <sample.fastq.gz> results/<sample_name> config/config.yaml
```

This produces `results/<sample_name>/<sample_name>.arm_counts.tsv`.

Build a healthy-control baseline from several such control samples
(more controls = more stable per-arm mean/SD; the paper's assay QC target
was >=90,000 usable reads/sample):

```bash
python3 scripts/build_control_baseline.py \
  --counts results/control_*/control_*.arm_counts.tsv \
  --out results/baseline.tsv
```

Compute the aneuploidy score for a case sample against that baseline:

```bash
python3 scripts/compute_aneuploidy_score.py \
  --counts results/<sample_name>/<sample_name>.arm_counts.tsv \
  --baseline results/baseline.tsv \
  --cutoff 5.0 \
  --out results/<sample_name>.aneuploidy.tsv
```

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
it's a toy genome, not a clinical validation.

## Key methodological differences from the published assay

- The published assay's exact LINE-1 primer sequence and error
  tolerances (Supplementary Table 1) are not reproduced here; you must
  supply your own validated primer.
- Chromosome-arm boundaries are hardcoded approximate GRCh38 centromere
  coordinates - verify against UCSC's `cytoBand`/`gap` tables before any
  diagnostic use.
- The paper used SPSS for downstream statistics (Z-scores, Cox models
  for survival); this repo only reproduces the FASTQ -> aneuploidy-score
  computation, not the clinical/survival analysis in the paper.
