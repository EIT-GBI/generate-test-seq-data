#!/usr/bin/env python3
"""Split one ONT sample across two barcodes.

Real runs multiplex a sample over several barcodes and hand the pipeline one
file per barcode. The samplesheet then carries several rows with the same
sample name, which is the case that tends to break output naming -- so there
needs to be data for it.
"""

import subprocess
from pathlib import Path

import config
from util import gzip_text


def main() -> None:
    sample = config.SAMPLES[0]
    for platform, organism in [(p, o) for p in ("ont", "pacbio")
                               for o in config.ORGANISMS]:
        source = Path(config.CACHE_DIR) / f"{platform}.{organism}.{sample}.sam"
        if not source.exists():
            raise SystemExit(f"missing {source}; run `make reads` first")

        header, records = [], []
        for line in source.read_text().splitlines():
            (header if line.startswith("@") else records).append(line)

        out_dir = Path(config.DATA_DIR) / platform / organism / "ubam"
        out_dir.mkdir(parents=True, exist_ok=True)
        for index, barcode in enumerate(config.BARCODES):
            part = records[index::len(config.BARCODES)]
            sam = Path(config.CACHE_DIR) / f"{platform}.{organism}.{sample}.{barcode}.sam"
            sam.write_text("\n".join(header + part) + "\n")
            bam = out_dir / f"{sample}.{barcode}.unaligned.bam"
            subprocess.run(["samtools", "view", "-b", "--no-PG", "-o", str(bam),
                            str(sam)], check=True)
            print(f"{platform}/{organism}: {barcode} {len(part)} reads -> {bam}")




def split_illumina() -> None:
    """Split one Illumina sample over two lanes.

    Illumina has no barcodes in this sense, but it has lanes, and they reach
    the samplesheet the same way: several rows sharing one sample name.
    """
    import gzip

    sample = config.SAMPLES[0]
    for organism in config.ORGANISMS:
        fastq_dir = Path(config.DATA_DIR) / "illumina" / organism / "fastq"
        for read in ("R1", "R2"):
            source = fastq_dir / f"{sample}_{read}.fastq.gz"
            with gzip.open(source, "rt") as handle:
                lines = handle.read().splitlines()
            records = [lines[i:i + 4] for i in range(0, len(lines), 4)]
            for index, lane in enumerate(["L001", "L002"]):
                part = records[index::2]
                target = fastq_dir / f"{sample}_{lane}_{read}.fastq.gz"
                with gzip_text(target) as handle:
                    handle.write("\n".join("\n".join(r) for r in part) + "\n")
        print(f"{organism}: {sample} split across L001/L002")


if __name__ == "__main__":
    main()
    split_illumina()
