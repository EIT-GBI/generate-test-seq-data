#!/usr/bin/env python3
"""Prove the generated data is fit for the tests that will use it.

This is not a unit test of the generators -- it is a check that the *data* has
the properties the pipelines' tests are going to assert:

  * every planted variant is recoverable by a plain bcftools call
  * no extra variants are called
  * the regions promised to be empty really have no reads
  * fastp has adapter to find
  * the MM/ML methylation calls round-trip through alignment at the levels
    truth_meth.bed claims

Run with `make validate`. A non-zero exit means the dataset is not safe to
publish.
"""

import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import config

SCRATCH = Path(config.CACHE_DIR) / "validate"
failures: list[str] = []


def check(condition: bool, message: str) -> None:
    print(("  ok   " if condition else "  FAIL ") + message)
    if not condition:
        failures.append(message)


def run(args, **kwargs):
    return subprocess.run(args, capture_output=True, text=True, check=True, **kwargs)


def normalise(vcf: Path, fasta: Path, out: Path) -> Path:
    run(["bcftools", "norm", "-f", str(fasta), "-m", "-any", "-Oz", "-o", str(out),
         str(vcf)])
    run(["bcftools", "index", "-f", str(out)])
    return out


def variant_recall(organism: str, platform: str = "illumina") -> None:
    fasta = Path(config.REFERENCE_DIR) / organism / "genome.fasta"
    truth = Path(config.REFERENCE_DIR) / organism / "truth.vcf"
    truth_gz = SCRATCH / f"{organism}.truth.vcf.gz"
    run(["bcftools", "view", "-Oz", "-o", str(truth_gz), str(truth)])
    run(["bcftools", "index", "-f", str(truth_gz)])

    for sample in config.SAMPLES:
        carried = SCRATCH / f"{organism}.{sample}.truth.vcf.gz"
        run(["bcftools", "view", "-s", sample, "-c1", "-Oz", "-o", str(carried),
             str(truth_gz)])
        run(["bcftools", "index", "-f", str(carried)])

        called = Path(config.DATA_DIR) / platform / organism / "vcf" / f"{sample}.vcf.gz"
        left = normalise(carried, fasta, SCRATCH / f"{organism}.{sample}.t.vcf.gz")
        right = normalise(called, fasta,
                          SCRATCH / f"{organism}.{sample}.{platform}.c.vcf.gz")

        out = SCRATCH / f"isec.{organism}.{sample}.{platform}"
        run(["bcftools", "isec", "-p", str(out), str(left), str(right)])
        counts = []
        for name in ("0000.vcf", "0001.vcf", "0002.vcf"):
            text = (out / name).read_text().splitlines()
            counts.append(sum(1 for line in text if not line.startswith("#")))
        missed, false_positive, recovered = counts
        check(missed == 0,
              f"{organism}/{platform}/{sample}: all {recovered} planted "
              f"variants recovered" + (f" ({missed} missed)" if missed else ""))
        check(false_positive == 0,
              f"{organism}/{platform}/{sample}: no false-positive calls"
              + (f" ({false_positive} extra)" if false_positive else ""))


def empty_region(organism: str) -> None:
    contig, start, end = config.ORGANISMS[organism]["empty_region"]
    for platform in ("illumina", "ont", "pacbio"):
        bam = Path(config.DATA_DIR) / platform / organism / "bam" / f"{config.SAMPLES[0]}.bam"
        if not bam.exists():
            continue
        result = run(["samtools", "view", "-c", str(bam), f"{contig}:{start + 1}-{end}"])
        check(int(result.stdout.strip()) == 0,
              f"{organism}/{platform}: {contig}:{start}-{end} has no reads, as promised")


def adapter_present(organism: str) -> None:
    fastq = Path(config.DATA_DIR) / "illumina" / organism / "fastq" / f"{config.SAMPLES[0]}_R1.fastq.gz"
    text = subprocess.run(["gzip", "-dc", str(fastq)], capture_output=True,
                          text=True, check=True).stdout.splitlines()
    sequences = text[1::4]
    hits = sum(1 for seq in sequences if "AGATCGGAAGAGC" in seq)
    fraction = hits / max(1, len(sequences))
    check(0.02 <= fraction <= 0.20,
          f"{organism}: {fraction:.1%} of R1 reads carry adapter (want 2-20%)")


CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_mm(mm: str, ml: str, seq: str):
    """Yield (read_offset, probability) for each 5mC call on a forward read."""
    body = mm.split(":", 2)[2].rstrip(";")
    if "," not in body:
        return
    deltas = [int(x) for x in body.split(",")[1:]]
    probs = [int(x) for x in ml.split(",")[1:]]
    c_positions = [i for i, base in enumerate(seq) if base == "C"]
    index = 0
    for delta, prob in zip(deltas, probs):
        index += delta
        if index >= len(c_positions):
            return
        yield c_positions[index], prob
        index += 1


def methylation_round_trip(organism: str) -> None:
    """Re-derive per-site methylation from the aligned BAM and compare.

    Only forward-strand primary alignments are used: their SEQ is the original
    read, so MM offsets apply directly and no reverse-complement bookkeeping is
    needed to check the simulation is sound.
    """
    bam = Path(config.DATA_DIR) / "ont" / organism / "bam" / f"{config.SAMPLES[0]}.bam"
    if not bam.exists():
        return
    truth = {}
    for line in (Path(config.REFERENCE_DIR) / organism / "truth_meth.bed").read_text().splitlines():
        fields = line.split("\t")
        truth[(fields[0], int(fields[1]))] = (float(fields[4]), fields[6])

    observed = defaultdict(lambda: [0, 0])
    stream = run(["samtools", "view", "-F", "0x914", str(bam)]).stdout
    for line in stream.splitlines():
        fields = line.split("\t")
        flag, contig, pos, cigar, seq = int(fields[1]), fields[2], int(fields[3]) - 1, fields[5], fields[9]
        if flag & 0x10:
            continue
        tags = {f[:5]: f for f in fields[11:] if f.startswith(("MM:Z:", "ML:B:"))}
        if "MM:Z:" not in tags or "ML:B:" not in tags:
            continue

        # read offset -> reference position, from the CIGAR
        mapping, read_i, ref_i = {}, 0, pos
        for count, op in CIGAR.findall(cigar):
            count = int(count)
            if op in "M=X":
                for k in range(count):
                    mapping[read_i + k] = ref_i + k
                read_i += count
                ref_i += count
            elif op in "IS":
                read_i += count
            elif op in "DN":
                ref_i += count

        for offset, prob in parse_mm(tags["MM:Z:"], tags["ML:B:"], seq):
            ref_pos = mapping.get(offset)
            if ref_pos is None:
                continue
            bucket = observed[(contig, ref_pos)]
            bucket[1] += 1
            if prob >= 128:
                bucket[0] += 1

    groups = defaultdict(list)
    for site, (modified, total) in observed.items():
        if total < 5 or site not in truth:
            continue
        _level, label = truth[site]
        groups[label].append(modified / total)

    for label, expected in (("methylated", config.METH["methylated_level"]),
                            ("unmethylated", config.METH["unmethylated_level"])):
        values = groups.get(label, [])
        if not values:
            check(False, f"{organism}: no {label} CpG sites recovered from the ONT BAM")
            continue
        mean = sum(values) / len(values)
        check(abs(mean - expected) < 0.10,
              f"{organism}: {len(values)} {label} CpGs average "
              f"{mean:.2f} modified (truth {expected:.2f})")


def samplesheet_urls() -> None:
    """Every URL in every samplesheet must name a file that exists here.

    raw.githubusercontent.com serves the repo tree verbatim, so a URL is just
    BASE_URL + the repo-relative path. Stripping BASE_URL back off must land on
    a real file -- otherwise the sheets are full of 404s that only surface
    inside a pipeline run.
    """
    sheets = sorted((Path(config.DATA_DIR) / "samplesheets").rglob("*.csv"))
    checked, broken = 0, []
    for sheet in sheets:
        rows = sheet.read_text().splitlines()
        for row in rows[1:]:
            for field in row.split(","):
                if not field.startswith("http"):
                    continue
                checked += 1
                if not field.startswith(config.BASE_URL + "/"):
                    broken.append(f"{sheet.name}: not under BASE_URL: {field}")
                    continue
                local = Path(field[len(config.BASE_URL) + 1:])
                # malformed.csv deliberately points at a file that is not there
                if not local.exists() and "NOT_A_FILE" not in field:
                    broken.append(f"{sheet.name}: no such file: {local}")

    check(not broken,
          f"{len(sheets)} samplesheets: all {checked} URLs resolve to real files"
          + (f" ({len(broken)} broken)" if broken else ""))
    for problem in broken[:10]:
        print(f"         {problem}")


def publishable_tree() -> None:
    """Nothing that must ship is excluded by .gitignore, and vice versa."""
    ignored = Path(".gitignore").read_text().split()
    check("data/" not in ignored and "data" not in ignored,
          ".gitignore does not exclude data/ -- it is the payload and must be committed")
    check(".cache/" in ignored,
          ".gitignore excludes .cache/ -- downloads and intermediates must not ship")


def main() -> None:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    for organism in config.ORGANISMS:
        print(f"\n{organism}")
        # Both platforms carry the same planted variants, so both must
        # recover them. PacBio is the weaker case: long reads over short
        # windows thin out at the edges.
        for platform in ("illumina", "pacbio"):
            variant_recall(organism, platform)
        empty_region(organism)
        adapter_present(organism)
        methylation_round_trip(organism)

    print("\npublishing")
    samplesheet_urls()
    publishable_tree()

    print()
    if failures:
        print(f"{len(failures)} check(s) failed")
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
