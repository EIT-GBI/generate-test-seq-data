#!/usr/bin/env python3
"""Simulate reads from the planted haplotypes.

One seeded generator for all three platforms, so the same variants are
recoverable from Illumina, ONT and PacBio data and a test can assert the same
truth.vcf whichever pipeline produced the calls.

The error models are deliberately simple -- enough to make trimming, alignment
and variant calling do real work, not enough to benchmark a caller:

  illumina  substitutions only, quality falling towards the 3' end, plus a
            real adapter on the tail of short fragments so FASTP_TRIM has
            something to find
  ont       substitutions, insertions and deletions, with deletions favoured
            inside homopolymers, and 5mC calls on CpGs as MM/ML tags
  pacbio    HiFi-like: rare substitutions, high quality throughout

Outputs FASTQ for every platform, and for the long-read platforms a SAM stream
on stdout that the Makefile turns into an unaligned BAM.
"""

import argparse
import gzip
import math
import random
import sys
from pathlib import Path

import config
from util import gzip_text

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
# Illumina TruSeq read-1 adapter: what fastp looks for by default
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
BASES = "ACGT"


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


class HapMap:
    """Translate haplotype coordinates back to reference coordinates."""

    def __init__(self, path: Path):
        self.marks: dict[str, list[tuple[int, int]]] = {}
        if not path.exists():
            return
        for line in path.read_text().splitlines():
            contig, ref_pos, shift = line.split("\t")
            self.marks.setdefault(contig, []).append((int(ref_pos), int(shift)))
        for marks in self.marks.values():
            marks.sort()

    def to_reference(self, contig: str, hap_pos: int) -> int:
        shift = 0
        for ref_pos, cumulative in self.marks.get(contig, []):
            # The mark sits at ref_pos, i.e. hap coordinate ref_pos + previous
            if ref_pos + shift <= hap_pos:
                shift = cumulative
            else:
                break
        return hap_pos - shift


def phred(quality: int) -> str:
    return chr(33 + max(2, min(60, quality)))


def homopolymer_ahead(seq: str, pos: int) -> int:
    run, base = 1, seq[pos]
    while pos + run < len(seq) and seq[pos + run] == base:
        run += 1
    return run


def mutate(rng, template, error_rate, qual_mean, indels, quality_slope=0.0):
    """Apply the error model. Returns (sequence, qualities, source positions).

    `sources[i]` is the offset within `template` that output base i came from,
    or None where the simulator inserted a base.
    """
    out, quals, sources = [], [], []
    index = 0
    while index < len(template):
        base = template[index]
        decayed = qual_mean - quality_slope * (len(out) / max(1, len(template)))
        quality = int(rng.gauss(decayed, 3))
        roll = rng.random()

        if roll < error_rate:
            kind = "sub"
            if indels:
                run = homopolymer_ahead(template, index)
                # Homopolymer length errors are the characteristic nanopore
                # failure, so bias towards them as the run gets longer
                weights = [0.35, 0.25, 0.40 + 0.08 * min(run, 6)]
                kind = rng.choices(["sub", "ins", "del"], weights=weights)[0]

            if kind == "sub":
                out.append(rng.choice([b for b in BASES if b != base]))
                quals.append(max(2, quality - 12))
                sources.append(index)
                index += 1
            elif kind == "ins":
                out.append(rng.choice(BASES))
                quals.append(max(2, quality - 10))
                sources.append(None)
            else:
                index += 1                      # deletion: emit nothing
        else:
            out.append(base)
            quals.append(quality)
            sources.append(index)
            index += 1

    return "".join(out), quals, sources


def methylation_tags(rng, seq, sources, hap_map, contig, hap_start,
                     template_len, reverse, meth_levels):
    """Build MM/ML for the 5mC calls on this read.

    A CpG is methylated on both strands, so a reverse-strand read carries the
    call on the C that pairs with the reference G -- which is why the site key
    shifts by one when `reverse` is set.
    """
    deltas, probabilities, skipped = [], [], 0

    for index, base in enumerate(seq):
        if base != "C":
            continue
        source = sources[index]
        if source is None:
            skipped += 1
            continue
        # On a reverse read `source` indexes the revcomp'd template, so walk
        # it back to a forward haplotype coordinate
        offset = (template_len - 1 - source) if reverse else source
        hap_pos = hap_start + offset
        ref_pos = hap_map.to_reference(contig, hap_pos)
        site = ref_pos - 1 if reverse else ref_pos

        level = meth_levels.get((contig, site))
        if level is None:
            skipped += 1
            continue

        if rng.random() < level:
            probability = rng.uniform(0.85, 0.99)
        else:
            probability = rng.uniform(0.01, 0.15)
        deltas.append(skipped)
        probabilities.append(min(255, int(probability * 256)))
        skipped = 0

    if not deltas:
        return None
    mm = "MM:Z:C+m?," + ",".join(str(d) for d in deltas) + ";"
    ml = "ML:B:C," + ",".join(str(p) for p in probabilities)
    return mm, ml


def load_meth_levels(organism):
    path = Path(config.REFERENCE_DIR) / organism / "truth_meth.bed"
    levels = {}
    if not path.exists():
        return levels
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        levels[(fields[0], int(fields[1]))] = float(fields[4])
    return levels


def simulate_illumina(rng, hap, regions, settings, sample, reads_per_region):
    """Yield (name, seq1, qual1, seq2, qual2) tuples."""
    length = settings["read_length"]
    counter = 0
    for contig, start, end in regions:
        seq = hap[contig]
        for _ in range(reads_per_region[(contig, start, end)]):
            counter += 1
            if rng.random() < settings.get("short_fragment_fraction", 0.0):
                fragment = rng.randint(60, length - 5)   # reads into the adapter
            else:
                fragment = int(rng.gauss(settings["insert_mean"],
                                         settings["insert_sd"]))
                fragment = max(length - 40, min(fragment, 900))
            position = rng.randrange(start, max(start + 1, end - fragment))
            template = seq[position:position + fragment]
            if len(template) < 50:
                continue

            # Short fragments read through into the adapter, which is what
            # gives FASTP_TRIM real work to do
            read1 = template[:length]
            read2 = revcomp(template)[:length]
            if fragment < length:
                read1 = (template + ADAPTER)[:length]
                read2 = (revcomp(template) + ADAPTER)[:length]

            slope = settings["qual_mean"] * 0.35
            seq1, q1, _ = mutate(rng, read1, settings["error_rate"],
                                 settings["qual_mean"], indels=False,
                                 quality_slope=slope)
            seq2, q2, _ = mutate(rng, read2, settings["error_rate"],
                                 settings["qual_mean"], indels=False,
                                 quality_slope=slope)
            name = f"{sample}:{contig}:{counter:07d}"
            yield name, seq1, q1, seq2, q2


def simulate_long(rng, hap, regions, settings, sample, reads_per_region,
                  hap_map, meth_levels, platform):
    for contig, start, end in regions:
        seq = hap[contig]
        counter = 0
        for _ in range(reads_per_region[(contig, start, end)]):
            counter += 1
            length = int(rng.gauss(settings["length_mean"], settings["length_sd"]))
            length = max(settings["length_min"], min(length, end - start))
            # Start reads before the window too, otherwise a read this long
            # covers the window edges far more thinly than its middle and
            # variants near the edge get missed. reads_needed() compensates
            # for the reads that land outside.
            lower = max(0, start - length + 1)
            upper = min(end, len(seq) - length)
            position = rng.randrange(lower, max(lower + 1, upper))
            template = seq[position:position + length]
            if len(template) < settings["length_min"]:
                continue

            reverse = rng.random() < 0.5
            oriented = revcomp(template) if reverse else template
            read, quals, sources = mutate(rng, oriented, settings["error_rate"],
                                          settings["qual_mean"],
                                          indels=(platform == "ont"))

            tags = None
            if platform == "ont" and meth_levels:
                tags = methylation_tags(rng, read, sources, hap_map, contig,
                                        position, len(template), reverse,
                                        meth_levels)
            name = f"{sample}:{contig}:{counter:06d}"
            yield name, read, quals, tags


def reads_needed(regions, settings, platform, ploidy):
    """Read counts per region that hit the configured coverage."""
    per_region = {}
    for contig, start, end in regions:
        span = end - start
        if platform == "illumina":
            total = span * settings["coverage"] / (2 * settings["read_length"])
        else:
            # Long reads are sampled from a window widened by one read length
            # at the start, so scale the count to keep in-window coverage at
            # the configured depth.
            widened = span + settings["length_mean"]
            total = widened * settings["coverage"] / settings["length_mean"]
        per_region[(contig, start, end)] = max(1, int(math.ceil(total / ploidy)))
    return per_region


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", required=True, choices=list(config.ORGANISMS))
    parser.add_argument("--platform", required=True, choices=list(config.PLATFORMS))
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sam", type=Path, help="also write an unaligned SAM here")
    args = parser.parse_args()

    spec = config.ORGANISMS[args.organism]
    settings = config.PLATFORMS[args.platform]
    ploidy = spec["ploidy"]
    meth_levels = (load_meth_levels(args.organism)
                   if args.platform == "ont" else {})

    args.outdir.mkdir(parents=True, exist_ok=True)
    hap_dir = Path(config.CACHE_DIR) / "haplotypes"

    if args.platform == "illumina":
        r1_path = args.outdir / f"{args.sample}_R1.fastq.gz"
        r2_path = args.outdir / f"{args.sample}_R2.fastq.gz"
        handles = (gzip_text(r1_path), gzip_text(r2_path))
    else:
        fq_path = args.outdir / f"{args.sample}.fastq.gz"
        handles = (gzip_text(fq_path),)

    sam_lines = []
    written = 0

    for hap_index in range(ploidy):
        name = f"{args.organism}.{args.sample}.hap{hap_index}"
        hap = read_fasta(hap_dir / f"{name}.fasta")
        hap_map = HapMap(hap_dir / f"{name}.map")
        regions = []
        for line in (hap_dir / f"{name}.regions").read_text().splitlines():
            contig, start, end = line.split("\t")
            regions.append((contig, int(start), int(end)))
        per_region = reads_needed(regions, settings, args.platform, ploidy)
        rng = random.Random(f"{config.SEED}:{name}:{args.platform}")

        if args.platform == "illumina":
            for read_name, seq1, q1, seq2, q2 in simulate_illumina(
                    rng, hap, regions, settings, args.sample, per_region):
                tag = f"{read_name}/h{hap_index}"
                handles[0].write(f"@{tag} 1:N:0:1\n{seq1}\n+\n"
                                 f"{''.join(phred(q) for q in q1)}\n")
                handles[1].write(f"@{tag} 2:N:0:1\n{seq2}\n+\n"
                                 f"{''.join(phred(q) for q in q2)}\n")
                written += 1
        else:
            for read_name, seq, quals, tags in simulate_long(
                    rng, hap, regions, settings, args.sample, per_region,
                    hap_map, meth_levels, args.platform):
                tag = f"{read_name}/h{hap_index}"
                quality = "".join(phred(q) for q in quals)
                handles[0].write(f"@{tag}\n{seq}\n+\n{quality}\n")
                if args.sam:
                    fields = [tag, "4", "*", "0", "0", "*", "*", "0", "0",
                              seq, quality]
                    if tags:
                        fields.extend(tags)
                    sam_lines.append("\t".join(fields))
                written += 1

    # Contamination: random sequence that maps nowhere. Appended last so it
    # does not disturb the seeded stream that produced the real reads.
    fraction = settings.get("unmappable_fraction", 0.0)
    if fraction:
        noise = random.Random(f"{config.SEED}:{args.organism}:{args.sample}:noise")
        count = max(1, int(written * fraction))
        length = settings.get("read_length") or settings["length_mean"]
        for index in range(count):
            name = f"{args.sample}:unmappable:{index:05d}"
            if args.platform == "illumina":
                for handle, read in zip(handles, (1, 2)):
                    seq = "".join(noise.choice(BASES) for _ in range(length))
                    quality = "".join(phred(int(noise.gauss(30, 3)))
                                      for _ in range(length))
                    handle.write(f"@{name} {read}:N:0:1\n{seq}\n+\n{quality}\n")
            else:
                seq = "".join(noise.choice(BASES) for _ in range(length))
                quality = "".join(phred(int(noise.gauss(settings["qual_mean"], 3)))
                                  for _ in range(length))
                handles[0].write(f"@{name}\n{seq}\n+\n{quality}\n")
                if args.sam:
                    sam_lines.append("\t".join([name, "4", "*", "0", "0", "*",
                                                 "*", "0", "0", seq, quality]))
        written += count

    for handle in handles:
        handle.close()

    if args.sam:
        header = ["@HD\tVN:1.6\tSO:unsorted",
                  f"@RG\tID:{args.sample}\tSM:{args.sample}\t"
                  f"PL:{'ONT' if args.platform == 'ont' else 'PACBIO'}"]
        args.sam.parent.mkdir(parents=True, exist_ok=True)
        args.sam.write_text("\n".join(header + sam_lines) + "\n")

    print(f"{args.organism}/{args.platform}/{args.sample}: {written} reads",
          file=sys.stderr)


if __name__ == "__main__":
    main()
