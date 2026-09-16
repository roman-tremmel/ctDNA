# ctDNA Aneuploidy Pipelines (LINE-1 amplicon sequencing)

Two independent, from-scratch FASTQ-to-result pipelines for estimating
tumor-agnostic ctDNA/aneuploidy signal from LINE-1 amplicon sequencing
data (FAST-SeqS chemistry, Kinde et al. 2012: a single primer pair
amplifies retrotransposon LINE-1 loci scattered across the genome from
very little input cfDNA). Each pipeline reimplements a different
published analysis strategy for turning those reads into an aneuploidy
call. They are unrelated code bases - use whichever fits your samples,
or run both and compare.

## [`mfast-seqs/`](mfast-seqs/README.md)

Reimplements **mFAST-SeqS** as used in Verschoor et al. 2023
(*npj Breast Cancer* 9:61) and Belic et al. 2016: per-chromosome-arm
read counts, normalized to library size, Z-scored against a
healthy-control panel-of-normals, summed as squared Z-scores into one
genome-wide **aneuploidy score** (cutoff 5.0 in the paper). Simple,
well-validated in the literature, but the between-sample panel-of-normals
comparison is more sensitive to batch effects.

## [`waldo/`](waldo/README.md)

An original reimplementation of the *ideas* behind **WALDO**
(Douville et al. 2018, *PNAS* 115(8):1871-1876): genome divided into
500kb windows, windows grouped into clusters that behave similarly
across a reference panel, and a **within-sample** test - a window's
expected reads are predicted from the *same test sample's* own reads in
its cluster, so a global depth/batch shift cancels out. Also includes an
optional SVM layer for classifying genome-wide aneuploidy status from
per-arm Z-scores (useful for low-tumor-fraction samples where no single
arm reaches significance). More complex, more tunable, and - since the
paper's exact statistical model is only in an SI Appendix that wasn't
available while building this - an approximation of the published
method rather than a byte-exact reproduction. See `waldo/README.md`'s
"Known limitations" before trusting its output.

## Method comparison

| | mfast-seqs/ | waldo/ |
|---|---|---|
| Resolution | Chromosome arm | 500kb windows -> clusters -> arm |
| Normalization | Between-sample vs. panel-of-normals | Within-sample vs. own cluster totals |
| Batch-effect robustness | Lower | Higher (by design) |
| Statistics | Sum of squared arm Z-scores -> 1 score | Per-arm Z-test + optional SVM classifier |
| Extra signals | None | (Not implemented here) allelic imbalance, sample fingerprinting, somatic mutations/MSI - see `waldo/README.md` |
| Validated sensitivity | Score >=5 tracks VAF >=5-10% (pre-screening / moderate-to-high ctDNA) | Down to ~1% neoplastic fraction at 99% specificity (SVM, per the paper) |
| Complexity to run/validate | Low | High |

## Requirements

Both pipelines need: `bwa`, `samtools`, `bedtools`, `cutadapt` on PATH,
and Python 3 with `numpy`, `pandas`, `scipy`, `pyyaml` (`waldo/` also
needs `scikit-learn` and `joblib` for its optional SVM layer). See each
subdirectory's `requirements.txt`/`environment.yml`.

## Running on an air-gapped / offline HPC cluster

Neither pipeline's analysis scripts make any network calls - once a
reference genome (+ bwa index) and the LINE-1 element BED exist locally,
everything runs offline. Internet is only needed for two one-time setup
steps (getting the reference genome, and building the LINE-1 BED from
UCSC's RepeatMasker track). See `mfast-seqs/README.md`'s "Running on an
air-gapped / offline HPC cluster" section for the exact commands to run
elsewhere and transfer over.

## Repository layout

```
mfast-seqs/   arm-level, between-sample panel-of-normals pipeline
waldo/        window/cluster-based, within-sample pipeline + SVM layer
```

Each has its own `config/`, `scripts/`, `resources/`, `tests/`, and
`README.md` - they can be copied/deployed independently of one another.
