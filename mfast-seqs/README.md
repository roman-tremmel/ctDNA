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
4. Copy `config/config.yaml` and point it at your reference/BED files.
   `primer.forward_seq` is pre-filled with the LINE-1-specific primer
   core extracted from Verschoor et al. 2023 Supplementary Table 1
   ("Primers for mFAST-SeqS") - see
   `resources/primers_supplementary_table1.tsv` for the full table and
   derivation. If your own assay uses a different primer, replace it.

## Running on an air-gapped / offline HPC cluster

Once `resources/genome.fa` (+ bwa index) and `resources/line1_elements.*.bed`
exist, **none of the pipeline scripts make any network calls** —
`trim_line1_reads.sh`, `align_reads.sh`, `count_reads_by_arm.py`,
`build_control_baseline.py` and `compute_aneuploidy_score.py` only read
local files and call locally installed tools (bwa/samtools/bedtools/
cutadapt). The only two steps in this repo that need internet access are:

1. Downloading the reference genome.
2. `resources/build_line1_bed.sh`, which fetches the UCSC RepeatMasker
   track (`hgdownload.soe.ucsc.edu`).

If your cluster has no outbound internet access, run those two steps on
any machine that does (your laptop, a login node with general internet,
etc.), then copy the results over:

```bash
# on a machine WITH internet access:
bwa index genome.fa
resources/build_line1_bed.sh line1_elements.hg38.bed

# transfer to the cluster:
rsync -avP genome.fa* line1_elements.hg38.bed \
  cluster:/path/to/ctDNA/resources/
```

`resources/chrom_arms.hg38.bed` is already committed to this repo, so it
does not need to be regenerated or transferred separately.

For the software itself, most HPC clusters already provide `bwa`,
`samtools`, `bedtools`, and `cutadapt` as environment modules
(`module load ...`); if not, build the `environment.yml` conda
environment once on a machine with internet and either recreate it on
the cluster from an internal conda channel/mirror, or ship it with
`conda-pack` (creates a relocatable, self-contained tarball that needs no
internet to unpack and activate).

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
