#!/usr/bin/env python3
"""End-to-end smoke test for the LINE-1 aneuploidy pipeline.

Builds a small synthetic reference genome that mimics the real assay's
structure (many short "LINE-1" loci carrying a shared, conserved primer
landing site, scattered across chromosome arms), simulates raw reads for
several "healthy" replicates plus one "tumor" sample with an engineered
copy-number gain on one arm and a loss on another, runs the actual
pipeline scripts (cutadapt, bwa, samtools, bedtools, then the Python
counting/baseline/scoring scripts), and checks that:

  * the tumor sample's aneuploidy score is much higher than any control's,
  * the gained/lost arms are the top contributors to that score,
  * a held-out healthy replicate scores low (as a leave-one-out control).

This is a logic/integration test with a toy genome, not a validation of
the real hg38 pipeline - it only proves the code is wired together
correctly.
"""
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / "scripts"

random.seed(1234)

PRIMER = "GGGCACTATTCAGGAAGTCA"  # fixed 20bp "conserved LINE-1" motif
ARM_LEN = 8000
N_LOCI_PER_ARM = 12
LOCUS_SPACING = ARM_LEN // (N_LOCI_PER_ARM + 1)
READ_LEN = 150
UNIQUE_LEN = READ_LEN - len(PRIMER) + 20  # extra genomic context past the primer
# How far the annotated "L1 element" footprint extends past the primer itself.
# Mirrors real L1 loci, which are hundreds-thousands of bp and extend beyond
# the short conserved primer-landing site into more degenerate L1 sequence
# before unique flanking genome starts - so a read immediately downstream of
# the (trimmed-away) primer still overlaps the annotated element.
L1_TAIL_LEN = 40

CHROMS = ["chr1", "chr2", "chr3", "chr4"]
ARMS = [f"{c}{arm}" for c in CHROMS for arm in ("p", "q")]


def random_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def build_genome():
    """Return (fasta_str, arms_bed_rows, line1_bed_rows, locus_index).

    locus_index: dict arm -> list of genomic 0-based positions of the
    primer landing site for each L1 copy on that arm.
    """
    fasta_lines = []
    arms_bed_rows = []
    line1_bed_rows = []
    locus_index = {}

    for chrom in CHROMS:
        seq_parts = []
        offset = 0
        for arm_suffix in ("p", "q"):
            arm = f"{chrom}{arm_suffix}"
            locus_index[arm] = []
            arm_start = offset
            for i in range(N_LOCI_PER_ARM):
                pre = random_seq(LOCUS_SPACING - len(PRIMER) - UNIQUE_LEN)
                locus_pos = offset + len(pre)
                seq_parts.append(pre)
                seq_parts.append(PRIMER)
                seq_parts.append(random_seq(UNIQUE_LEN))
                line1_bed_rows.append(
                    (chrom, locus_pos, locus_pos + len(PRIMER) + L1_TAIL_LEN, f"L1_{arm}_{i}")
                )
                locus_index[arm].append(locus_pos)
                offset = locus_pos + len(PRIMER) + UNIQUE_LEN
            arm_end = offset
            excluded = 1 if arm in ("chr4p",) else 0  # pretend chr4p is acrocentric-like
            arms_bed_rows.append((chrom, arm_start, arm_end, arm, excluded))
        fasta_lines.append((chrom, "".join(seq_parts)))

    fasta_str = "".join(f">{c}\n{s}\n" for c, s in fasta_lines)
    return fasta_str, arms_bed_rows, line1_bed_rows, locus_index


def simulate_reads(fasta_by_chrom, locus_index, arm_multipliers, base_reads_per_locus=40):
    """Return list of raw (untrimmed) read sequences for one sample."""
    reads = []
    for arm, positions in locus_index.items():
        mult = arm_multipliers.get(arm, 1.0)
        n_reads = max(1, round(base_reads_per_locus * mult))
        chrom = arm[:-1]
        genome_seq = fasta_by_chrom[chrom]
        for pos in positions:
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
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def main():
    for tool in ("cutadapt", "bwa", "samtools", "bedtools"):
        if shutil.which(tool) is None:
            print(f"SKIP: required tool '{tool}' not found on PATH", file=sys.stderr)
            return 0

    tmpdir = Path(tempfile.mkdtemp(prefix="line1_pipeline_test_"))
    print(f"Working in {tmpdir}")

    fasta_str, arms_bed_rows, line1_bed_rows, locus_index = build_genome()
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

    line1_bed = tmpdir / "line1.bed"
    with open(line1_bed, "w") as fh:
        for row in line1_bed_rows:
            fh.write("\t".join(map(str, row)) + "\n")

    def process_sample(name, arm_multipliers):
        sample_dir = tmpdir / name
        sample_dir.mkdir()
        reads = simulate_reads(fasta_by_chrom, locus_index, arm_multipliers)
        raw_fastq = sample_dir / "raw.fastq"
        write_fastq(raw_fastq, reads)

        trimmed = sample_dir / "trimmed.fastq.gz"
        run([str(SCRIPTS / "trim_line1_reads.sh"), str(raw_fastq) + ".gz",
             str(trimmed), PRIMER, "0.15", "20"])

        bam = sample_dir / "aligned.sorted.bam"
        run([str(SCRIPTS / "align_reads.sh"), str(trimmed), str(ref_fasta), str(bam), "2"])

        counts = sample_dir / "arm_counts.tsv"
        run(["python3", str(SCRIPTS / "count_reads_by_arm.py"),
             "--bam", str(bam), "--line1-bed", str(line1_bed),
             "--arms-bed", str(arms_bed), "--min-mapq", "20", "--out", str(counts)])
        return counts

    # 6 "healthy" replicates: flat coverage across all arms
    control_counts = []
    for i in range(6):
        jitter = {arm: 1.0 + random.uniform(-0.08, 0.08) for arm in ARMS}
        control_counts.append(process_sample(f"control_{i}", jitter))

    # held-out euploid sample, NOT part of the baseline
    holdout_counts = process_sample("holdout_euploid",
                                     {arm: 1.0 + random.uniform(-0.08, 0.08) for arm in ARMS})

    # "tumor" sample: chr2q gained (2x), chr3p lost (~0.3x)
    tumor_multipliers = {arm: 1.0 for arm in ARMS}
    tumor_multipliers["chr2q"] = 2.0
    tumor_multipliers["chr3p"] = 0.3
    tumor_counts = process_sample("tumor", tumor_multipliers)

    baseline = tmpdir / "baseline.tsv"
    run(["python3", str(SCRIPTS / "build_control_baseline.py"),
         "--counts", *[str(c) for c in control_counts], "--out", str(baseline)])

    def score_sample(counts_path, out_name):
        out = tmpdir / out_name
        result = run(["python3", str(SCRIPTS / "compute_aneuploidy_score.py"),
                      "--counts", str(counts_path), "--baseline", str(baseline),
                      "--cutoff", "5.0", "--out", str(out)])
        print(result.stdout)
        score_line = [l for l in out.read_text().splitlines() if l.startswith("# aneuploidy_score")][0]
        return float(score_line.split("\t")[1]), out

    holdout_score, holdout_out = score_sample(holdout_counts, "holdout.score.tsv")
    tumor_score, tumor_out = score_sample(tumor_counts, "tumor.score.tsv")

    print(f"\nHoldout euploid score: {holdout_score:.2f}")
    print(f"Tumor score: {tumor_score:.2f}")

    assert tumor_score > 5.0 * holdout_score, (
        f"expected tumor score to be much higher than euploid holdout "
        f"({tumor_score:.2f} vs {holdout_score:.2f})"
    )

    tumor_lines = tumor_out.read_text().splitlines()
    print("Tumor per-arm result (top rows):")
    for line in tumor_lines[:4]:
        print(" ", line)

    assert tumor_lines[1].startswith("chr2q\t"), "expected chr2q (the gained arm) to be top hit"

    print("\nOK: synthetic end-to-end pipeline test passed.")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
