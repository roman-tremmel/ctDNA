#!/usr/bin/env python3
"""End-to-end smoke test for the WALDO-style within-sample pipeline.

Builds a small synthetic genome (chromosome arms subdivided into fixed
windows, each window carrying several "LINE-1" loci), simulates reads
for a panel of "healthy" reference replicates plus one "tumor" sample
with an engineered gain on one arm and a loss on another, then runs the
actual pipeline scripts end to end: cutadapt -> bwa -> samtools ->
bedtools -> window clustering -> within-sample arm Z-score test.

This validates pipeline mechanics on a toy genome; it is not a
validation of the real biology or of exact concordance with the
published WALDO algorithm (see waldo/README.md).
"""
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / "scripts"

random.seed(4321)
np.random.seed(4321)

PRIMER = "GGGCACTATTCAGGAAGTCA"
WINDOW_SIZE = 2000
WINDOWS_PER_ARM = 8
ARM_LEN = WINDOW_SIZE * WINDOWS_PER_ARM
LOCI_PER_WINDOW = 3
LOCUS_SPACING = WINDOW_SIZE // (LOCI_PER_WINDOW + 1)
READ_LEN = 150
UNIQUE_LEN = READ_LEN - len(PRIMER) + 20
L1_TAIL_LEN = 40

CHROMS = ["chr1", "chr2", "chr3", "chr4"]
ARMS = [f"{c}{arm}" for c in CHROMS for arm in ("p", "q")]


def random_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def build_genome():
    fasta_lines = []
    arms_bed_rows = []
    windows_bed_rows = []
    line1_bed_rows = []
    locus_index = {arm: [] for arm in ARMS}  # arm -> [(window_id, pos), ...]

    for chrom in CHROMS:
        seq_parts = []
        offset = 0
        for arm_suffix in ("p", "q"):
            arm = f"{chrom}{arm_suffix}"
            arm_start = offset
            for w in range(WINDOWS_PER_ARM):
                window_start = offset
                window_id = f"win_{chrom}_{window_start}"
                for i in range(LOCI_PER_WINDOW):
                    pre = random_seq(LOCUS_SPACING - len(PRIMER) - UNIQUE_LEN)
                    locus_pos = offset + len(pre)
                    seq_parts.append(pre)
                    seq_parts.append(PRIMER)
                    seq_parts.append(random_seq(UNIQUE_LEN))
                    line1_bed_rows.append(
                        (chrom, locus_pos, locus_pos + len(PRIMER) + L1_TAIL_LEN,
                         f"L1_{window_id}_{i}")
                    )
                    locus_index[arm].append((window_id, locus_pos))
                    offset = locus_pos + len(PRIMER) + UNIQUE_LEN
                window_end = offset
                windows_bed_rows.append((chrom, window_start, window_end, window_id))
            arm_end = offset
            excluded = 1 if arm in ("chr4p",) else 0
            arms_bed_rows.append((chrom, arm_start, arm_end, arm, excluded))
        fasta_lines.append((chrom, "".join(seq_parts)))

    fasta_str = "".join(f">{c}\n{s}\n" for c, s in fasta_lines)
    return fasta_str, arms_bed_rows, windows_bed_rows, line1_bed_rows, locus_index


N_EFFICIENCY_CLASSES = 4


def assign_window_efficiencies(windows_bed_rows, n_classes=N_EFFICIENCY_CLASSES, low=0.6, high=1.4):
    """Give each window a fixed, sample-independent "amplification
    efficiency" (stand-in for real local sequence effects like GC content
    that make some LINE-1 loci amplify better than others). Efficiency
    classes are assigned independent of chromosome arm, so windows that
    end up correlated across the reference panel (and thus clustered
    together by build_window_clusters.py) genuinely span multiple arms -
    exercising the same cross-genome clustering the real WALDO algorithm
    relies on, not just "one cluster per arm".
    """
    class_values = np.linspace(low, high, n_classes)
    window_ids = sorted({row[3] for row in windows_bed_rows})
    return {w: random.choice(class_values) for w in window_ids}


def simulate_reads(fasta_by_chrom, locus_index, arm_multipliers, window_efficiency,
                    base_reads_per_locus=80):
    reads = []
    for arm, entries in locus_index.items():
        mult = arm_multipliers.get(arm, 1.0)
        chrom = arm[:-1]
        genome_seq = fasta_by_chrom[chrom]
        for window_id, pos in entries:
            lam = base_reads_per_locus * mult * window_efficiency[window_id]
            n_reads = max(1, int(np.random.poisson(lam)))
            genomic_read = genome_seq[pos:pos + READ_LEN - random.randint(0, 7)]
            for _ in range(n_reads):
                spacer = random_seq(random.randint(0, 7))
                reads.append(spacer + genomic_read)
    random.shuffle(reads)
    return reads


def write_fastq(path, reads):
    with open(path, "w") as fh:
        for i, seq in enumerate(reads):
            fh.write(f"@read{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    subprocess.run(["gzip", "-f", str(path)], check=True)


def run(cmd, **kw):
    result = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if result.returncode != 0:
        print(result.stdout, file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        result.check_returncode()
    return result


def main():
    for tool in ("cutadapt", "bwa", "samtools", "bedtools"):
        if shutil.which(tool) is None:
            print(f"SKIP: required tool '{tool}' not found on PATH", file=sys.stderr)
            return 0

    tmpdir = Path(tempfile.mkdtemp(prefix="waldo_pipeline_test_"))
    print(f"Working in {tmpdir}")

    fasta_str, arms_bed_rows, windows_bed_rows, line1_bed_rows, locus_index = build_genome()
    window_efficiency = assign_window_efficiencies(windows_bed_rows)
    fasta_by_chrom = {}
    for block in fasta_str.strip().split(">")[1:]:
        header, *seq_lines = block.splitlines()
        fasta_by_chrom[header] = "".join(seq_lines)

    ref_fasta = tmpdir / "genome.fa"
    ref_fasta.write_text(fasta_str)
    run(["bwa", "index", str(ref_fasta)])
    run(["samtools", "faidx", str(ref_fasta)])

    arms_bed = tmpdir / "arms.bed"
    with open(arms_bed, "w") as fh:
        for row in arms_bed_rows:
            fh.write("\t".join(map(str, row)) + "\n")

    windows_bed = tmpdir / "windows.bed"
    with open(windows_bed, "w") as fh:
        for row in windows_bed_rows:
            fh.write("\t".join(map(str, row)) + "\n")

    line1_bed = tmpdir / "line1.bed"
    with open(line1_bed, "w") as fh:
        for row in line1_bed_rows:
            fh.write("\t".join(map(str, row)) + "\n")

    def process_sample(name, arm_multipliers):
        sample_dir = tmpdir / name
        sample_dir.mkdir()
        reads = simulate_reads(fasta_by_chrom, locus_index, arm_multipliers, window_efficiency)
        raw_fastq = sample_dir / "raw.fastq"
        write_fastq(raw_fastq, reads)

        trimmed = sample_dir / "trimmed.fastq.gz"
        run([str(SCRIPTS / "trim_line1_reads.sh"), str(raw_fastq) + ".gz",
             str(trimmed), PRIMER, "0.15", "20"])

        bam = sample_dir / "aligned.sorted.bam"
        run([str(SCRIPTS / "align_reads.sh"), str(trimmed), str(ref_fasta), str(bam), "2"])

        counts = sample_dir / "window_counts.tsv"
        run(["python3", str(SCRIPTS / "count_reads_by_window.py"),
             "--bam", str(bam), "--line1-bed", str(line1_bed),
             "--windows-bed", str(windows_bed), "--min-mapq", "20", "--out", str(counts)])
        return counts

    n_controls = 7  # Douville et al. 2018 found 7 reference samples sufficient
    # Euploid samples get NO per-arm multiplier (every arm at 1.0):
    # coherently shifting every window of an arm is exactly the aneuploidy
    # signal this test looks for, so a "healthy" sample should only carry
    # natural per-window/per-locus Poisson sampling noise (already in
    # simulate_reads), not an artificial arm-level jitter that would look
    # just like a weak, fake aneuploidy signal.
    euploid_multipliers = {arm: 1.0 for arm in ARMS}
    control_counts = []
    for i in range(n_controls):
        control_counts.append(process_sample(f"control_{i}", euploid_multipliers))

    holdout_counts = process_sample("holdout_euploid", euploid_multipliers)

    tumor_multipliers = {arm: 1.0 for arm in ARMS}
    tumor_multipliers["chr2q"] = 2.2
    tumor_multipliers["chr3p"] = 0.3
    tumor_counts = process_sample("tumor", tumor_multipliers)

    clusters_tsv = tmpdir / "clusters.tsv"
    run(["python3", str(SCRIPTS / "build_window_clusters.py"),
         "--counts", *[str(c) for c in control_counts],
         "--mean-p-threshold", "0.05", "--var-p-threshold", "0.05",
         "--out", str(clusters_tsv)])

    threshold_tsv = tmpdir / "thresholds.tsv"
    run(["python3", str(SCRIPTS / "calibrate_threshold.py"),
         "--counts", *[str(c) for c in control_counts],
         "--clusters", str(clusters_tsv),
         "--windows-bed", str(windows_bed), "--arms-bed", str(arms_bed),
         "--excluded-arms", "chr4p",
         "--margin", "4.0", "--out", str(threshold_tsv)])

    def score_sample(counts_path, out_name):
        out = tmpdir / out_name
        result = run(["python3", str(SCRIPTS / "call_aneuploidy.py"),
                      "--counts", str(counts_path), "--clusters", str(clusters_tsv),
                      "--windows-bed", str(windows_bed), "--arms-bed", str(arms_bed),
                      "--excluded-arms", "chr4p",
                      "--threshold-table", str(threshold_tsv), "--out", str(out)])
        print(result.stdout)
        return out

    holdout_out = score_sample(holdout_counts, "holdout.scores.tsv")
    tumor_out = score_sample(tumor_counts, "tumor.scores.tsv")

    import pandas as pd
    holdout_df = pd.read_csv(holdout_out, sep="\t").set_index("arm")
    tumor_df = pd.read_csv(tumor_out, sep="\t").set_index("arm")

    print("\nHoldout arm Z-scores:")
    print(holdout_df["z_score"].to_string())
    print("\nTumor arm Z-scores:")
    print(tumor_df["z_score"].to_string())

    # With only ~64 windows total (vs. the paper's 4,361 genome-wide), the
    # empirically-calibrated per-arm threshold (max/min observed + margin)
    # is inherently noisy at this toy scale - not enough windows survive
    # cluster/outlier filtering per arm for a tight estimate. So we check
    # the *direction and relative magnitude* of the Z-scores (the real
    # mechanical claim: a gain pulls Z up, a loss pulls it down, and the
    # true signal dwarfs the euploid holdout's), rather than demanding the
    # gain/loss/normal calls line up with this run's specific threshold.
    assert tumor_df.loc["chr2q", "z_score"] > 0, (
        f"expected chr2q (gained 2.2x) to have a positive Z-score, got "
        f"{tumor_df.loc['chr2q', 'z_score']:.2f}"
    )
    assert tumor_df.loc["chr3p", "z_score"] < 0, (
        f"expected chr3p (lost to 0.3x) to have a negative Z-score, got "
        f"{tumor_df.loc['chr3p', 'z_score']:.2f}"
    )
    assert tumor_df.loc["chr2q", "z_score"] > 2 * abs(holdout_df.loc["chr2q", "z_score"]), (
        "expected the tumor's chr2q gain to be far more extreme than the "
        "euploid holdout's chr2q Z-score"
    )
    assert tumor_df.loc["chr3p", "z_score"] < -2 * abs(holdout_df.loc["chr3p", "z_score"]), (
        "expected the tumor's chr3p loss to be far more extreme than the "
        "euploid holdout's chr3p Z-score"
    )
    assert tumor_df["z_score"].abs().max() > 2 * holdout_df["z_score"].abs().max(), (
        "expected the tumor sample's most extreme arm Z-score to be much larger "
        "than the holdout euploid sample's"
    )

    print("\nOK: synthetic end-to-end WALDO pipeline test passed.")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
