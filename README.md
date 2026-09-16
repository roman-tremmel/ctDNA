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

Reimplements **WALDO** (Douville et al. 2018, *PNAS* 115(8):1871-1876)
following its SI Appendix, SI Materials and Methods: genome divided into
500kb windows, each window's cluster = other windows on *different*
chromosomes that are statistically indistinguishable from it (paired
t-test + F-test) across a reference panel, and a **within-sample** test
- each window's cluster mean/variance are estimated from the *same test
sample's own* values at its cluster members (with iterative outlier
trimming), so a global depth/batch shift cancels out and no arm's own
windows can ever make up its whole comparison group. A significance
threshold is calibrated empirically from a reference panel (max/min
observed Z + margin, per the paper) rather than assumed. Also includes
an optional SVM layer for classifying genome-wide aneuploidy status from
per-arm Z-scores (useful for low-tumor-fraction samples where no single
arm reaches significance). More complex, more tunable, and - since some
of the SI's phrasing is genuinely ambiguous on one point (see
`waldo/README.md`'s "An acknowledged ambiguity") - not guaranteed to be
byte-exact, but built directly from the published methods text rather
than approximated from the main-text summary alone. See
`waldo/README.md`'s "Known limitations" before trusting its output.

## Method comparison

| | mfast-seqs/ | waldo/ |
|---|---|---|
| Resolution | Chromosome arm | 500kb windows -> per-window clusters -> arm |
| Normalization | Between-sample vs. panel-of-normals | Within-sample: test sample's own values at its clusters' member windows |
| Cluster membership | N/A | Reference-panel paired t-test + F-test, restricted to other chromosomes |
| Batch-effect robustness | Lower | Higher (by design) |
| Significance threshold | Fixed literature cutoff (5.0) | Empirically calibrated per arm from a reference panel (max/min Z + margin) |
| Statistics | Sum of squared arm Z-scores -> 1 score | Per-arm Z-test + optional SVM classifier |
| Extra signals | None | (Not implemented here) allelic imbalance, sample fingerprinting, somatic mutations/MSI - see `waldo/README.md` |
| Validated sensitivity | Score >=5 tracks VAF >=5-10% (pre-screening / moderate-to-high ctDNA) | Down to ~1% neoplastic fraction at 99% specificity (SVM, per the paper) |
| Complexity to run/validate | Low | High |

## Requirements

Both pipelines need: `bwa`, `samtools`, `bedtools`, `cutadapt` on PATH,
and Python 3 with `numpy`, `pandas`, `scipy`, `pyyaml` (`waldo/` also
needs `scikit-learn` and `joblib` for its optional SVM layer). The
top-level `environment.yml` covers the union of both (one env for
everything); each subdirectory also has its own `requirements.txt`/
`environment.yml` if you only want one pipeline.

## Installing without root / admin rights: conda or mamba

None of this needs system package managers or sudo. On an HPC cluster
where you can't `apt-get`/`brew` install anything, a user-space
[Miniforge](https://github.com/conda-forge/miniforge) (conda + mamba,
preconfigured for conda-forge/bioconda) installs entirely under your
home directory:

```bash
curl -fsSL -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"

mamba env create -f environment.yml
conda activate ctdna-aneuploidy
```

That's the entire install - `bwa`, `samtools`, `bedtools`, `cutadapt`
and all Python dependencies for both pipelines land in one environment
in under two minutes, no root needed at any step.

## Running on an air-gapped / offline HPC cluster

Neither pipeline's analysis scripts make any network calls - once a
reference genome (+ bwa index), the LINE-1 element BED, and the conda
environment exist locally, everything runs offline. Three things need
internet at some point, all one-time and all doable on a *different*
machine than the offline cluster if needed:

1. **The reference genome** (see `mfast-seqs/README.md`'s "Running on
   an air-gapped / offline HPC cluster" section for exact commands).
2. **The LINE-1 element BED**, built from UCSC's RepeatMasker track
   (`resources/build_line1_bed.sh` in either pipeline).
3. **The conda environment itself**, if even `mamba env create` can't
   reach conda-forge/bioconda from the cluster. Build it once on a
   machine with internet, then ship it as a self-contained archive with
   [`conda-pack`](https://conda.github.io/conda-pack/) - no conda
   install, no internet, and no root needed to unpack and use it on the
   offline machine:

   ```bash
   # on a machine WITH internet access, after `mamba env create -f environment.yml`:
   mamba install -n base -y -c conda-forge conda-pack
   conda pack -n ctdna-aneuploidy -o ctdna-aneuploidy-env.tar.gz

   # transfer ctdna-aneuploidy-env.tar.gz to the offline cluster, then there:
   mkdir -p ~/ctdna-env && tar -xzf ctdna-aneuploidy-env.tar.gz -C ~/ctdna-env
   source ~/ctdna-env/bin/activate
   conda-unpack   # fixes up hardcoded paths for the new location - run once
   # bwa, samtools, bedtools, cutadapt, python + all deps now just work
   ```

   This was built and verified end-to-end in this repo: both pipelines'
   synthetic tests pass when run from a conda-pack archive re-extracted
   to a completely different path, with no conda and no internet
   involved at that point at all.

## Repository layout

```
mfast-seqs/   arm-level, between-sample panel-of-normals pipeline
waldo/        window/cluster-based, within-sample pipeline + SVM layer
```

Each has its own `config/`, `scripts/`, `resources/`, `tests/`, and
`README.md` - they can be copied/deployed independently of one another.
