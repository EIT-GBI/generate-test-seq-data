#!/usr/bin/env python3
"""Write the samplesheets, with published URLs baked in.

Paths point at config.BASE_URL, so a pipeline can use a sheet with no
rewriting:

    nextflow run . --samplesheet <BASE_URL>/data/samplesheets/illumina/single.csv \
                   --reference_dir <BASE_URL>/data/references

Every sheet exists to exercise something specific; see the table in README.md.
"""

from pathlib import Path

import config

BASE = config.BASE_URL
REFS = f"{BASE}/data/references"


def illumina_reads(organism, sample, lane=None):
    stem = f"{sample}_{lane}" if lane else sample
    base = f"{BASE}/data/illumina/{organism}/fastq/{stem}"
    return f"{base}_R1.fastq.gz", f"{base}_R2.fastq.gz"


def long_reads(platform, organism, sample, barcode=None):
    stem = f"{sample}.{barcode}" if barcode else sample
    return f"{BASE}/data/{platform}/{organism}/ubam/{stem}.unaligned.bam"


def write(path: Path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [",".join(header)] + [",".join(row) for row in rows]
    path.write_text("\n".join(lines) + "\n")
    print(f"  {path}  ({len(rows)} rows)")


def illumina_sheets(out: Path):
    header = ["sample", "R1", "R2", "reference"]
    ref = "ecoli/genome.fasta"

    write(out / "single.csv", header,
          [["SAMPLE1", *illumina_reads("ecoli", "SAMPLE1"), ref]])

    write(out / "two_samples.csv", header,
          [[s, *illumina_reads("ecoli", s), ref] for s in config.SAMPLES])

    # Two samples on two different genomes in one run: this is what exercises
    # every `combine(..., by: 0)` keyed on the reference name.
    write(out / "multi_reference.csv", header, [
        ["ECOLI1", *illumina_reads("ecoli", "SAMPLE1"), "ecoli/genome.fasta"],
        ["HUMAN1", *illumina_reads("human", "SAMPLE1"), "human/genome.fasta"],
    ])

    # A reference with no index beside it, so BWA_INDEX and SAMTOOLS_FAIDX
    # have to run. Only the first target region is present in this fasta, so
    # roughly a third of the reads map -- that is expected.
    write(out / "bare_reference.csv", header,
          [["SAMPLE1", *illumina_reads("ecoli", "SAMPLE1"),
            "ecoli-no-index/genome.fasta"]])

    # One sample, several files: the case that collides on output naming
    write(out / "multiple_lanes_per_sample.csv", header,
          [["SAMPLE1", *illumina_reads("ecoli", "SAMPLE1", lane), ref]
           for lane in ("L001", "L002")])

    write(out / "with_platform.csv", header + ["platform"],
          [[s, *illumina_reads("ecoli", s), ref, "ILLUMINA"]
           for s in config.SAMPLES])

    # Deliberately broken: no R2 column, and a file that does not exist. The
    # pipeline should fail fast and say why.
    write(out / "malformed.csv", ["sample", "R1", "reference"],
          [["SAMPLE1", illumina_reads("ecoli", "SAMPLE1")[0], ref],
           ["SAMPLE2", f"{BASE}/data/illumina/ecoli/fastq/NOT_A_FILE_R1.fastq.gz", ref]])


def long_sheets(platform: str, out: Path):
    header = ["sample", "reads", "reference"]
    ref = "ecoli/genome.fasta"

    write(out / "single.csv", header,
          [["SAMPLE1", long_reads(platform, "ecoli", "SAMPLE1"), ref]])

    write(out / "two_samples.csv", header,
          [[s, long_reads(platform, "ecoli", s), ref] for s in config.SAMPLES])

    write(out / "multi_reference.csv", header, [
        ["ECOLI1", long_reads(platform, "ecoli", "SAMPLE1"), "ecoli/genome.fasta"],
        ["HUMAN1", long_reads(platform, "human", "SAMPLE1"), "human/genome.fasta"],
    ])

    write(out / "bare_reference.csv", header,
          [["SAMPLE1", long_reads(platform, "ecoli", "SAMPLE1"),
            "ecoli-no-index/genome.fasta"]])

    write(out / "multiple_barcodes_per_sample.csv", header,
          [["SAMPLE1", long_reads(platform, "ecoli", "SAMPLE1", barcode), ref]
           for barcode in config.BARCODES])

    write(out / "with_platform.csv", header + ["platform"],
          [[s, long_reads(platform, "ecoli", s), ref,
            "ONT" if platform == "ont" else "PACBIO"] for s in config.SAMPLES])

    write(out / "malformed.csv", ["sample", "reference"],
          [["SAMPLE1", ref],
           ["SAMPLE2", ref]])

    # FASTQ rather than unaligned BAM, for the minimap2/chopper entry point
    write(out / "fastq_input.csv", header,
          [[s, f"{BASE}/data/{platform}/ecoli/fastq/{s}.fastq.gz", ref]
           for s in config.SAMPLES])


def main() -> None:
    root = Path(config.DATA_DIR) / "samplesheets"
    print("illumina:")
    illumina_sheets(root / "illumina")
    for platform in ("ont", "pacbio"):
        print(f"{platform}:")
        long_sheets(platform, root / platform)


if __name__ == "__main__":
    main()
