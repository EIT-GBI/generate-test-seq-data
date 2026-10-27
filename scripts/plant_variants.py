#!/usr/bin/env python3
"""Plant a known set of variants into each reference and write the truth.

Outputs per organism:
  data/references/<org>/truth.vcf   multi-sample VCF: exactly what was planted
  data/references/<org>/regions.bed the windows reads are simulated from
  .cache/haplotypes/<org>.<sample>.hap<n>.fasta  what reads are drawn from

The haplotypes are intermediates and are not committed. truth.vcf is the point
of the whole exercise: it lets a test assert that a caller recovered a specific
variant, rather than that an md5 did not change.
"""

import argparse
import random
from pathlib import Path

import config

BASES = "ACGT"


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


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def homopolymer_run(seq: str, pos: int) -> int:
    """Length of the homopolymer run that `pos` sits in."""
    base, left, right = seq[pos], pos, pos
    while left > 0 and seq[left - 1] == base:
        left -= 1
    while right + 1 < len(seq) and seq[right + 1] == base:
        right += 1
    return right - left + 1


def pick_positions(rng, seqs, regions, count, spacing, avoid_homopolymer):
    """Deterministically scatter `count` positions across the target regions."""
    chosen, used = [], []
    attempts = 0
    while len(chosen) < count and attempts < count * 400:
        attempts += 1
        contig, start, end = regions[len(chosen) % len(regions)]
        margin = min(config.VARIANT_MARGIN, (end - start) // 3)
        pos = rng.randrange(start + margin, end - margin)
        seq = seqs[contig]
        if seq[pos] not in BASES:
            continue
        if avoid_homopolymer and homopolymer_run(seq, pos) > 3:
            continue        # indels in homopolymers are exactly what callers
                            # disagree about; keep the truth unambiguous
        if any(c == contig and abs(p - pos) < spacing for c, p in used):
            continue
        used.append((contig, pos))
        chosen.append((contig, pos))
    if len(chosen) < count:
        raise SystemExit(f"could only place {len(chosen)}/{count} positions")
    return sorted(chosen)


def make_variants(rng, seqs, spec):
    """Build the variant list: (contig, pos0, ref, alt, kind)."""
    variants = []

    for contig, pos in pick_positions(rng, seqs, spec["regions"], spec["n_snps"],
                                      spacing=400, avoid_homopolymer=False):
        ref = seqs[contig][pos]
        alt = rng.choice([b for b in BASES if b != ref])
        variants.append((contig, pos, ref, alt, "snp"))

    for contig, pos in pick_positions(rng, seqs, spec["regions"], spec["n_indels"],
                                      spacing=1200, avoid_homopolymer=True):
        seq = seqs[contig]
        if rng.random() < 0.5:
            length = rng.choice([1, 2, 3, 6])          # deletion
            ref = seq[pos:pos + length + 1]
            alt = seq[pos]
            kind = f"del{length}"
        else:
            length = rng.choice([1, 2, 3, 5])          # insertion
            ref = seq[pos]
            alt = seq[pos] + "".join(rng.choice(BASES) for _ in range(length))
            kind = f"ins{length}"
        if "N" in ref:
            continue
        variants.append((contig, pos, ref, alt, kind))

    # Drop any variant whose REF span overlaps another: overlapping records are
    # legal VCF but make "did the caller find it" ambiguous
    variants.sort(key=lambda v: (v[0], v[1]))
    kept = []
    for variant in variants:
        if kept and kept[-1][0] == variant[0] and \
           variant[1] <= kept[-1][1] + len(kept[-1][2]):
            continue
        kept.append(variant)
    return kept


def assign_genotypes(rng, variants, ploidy):
    """Give every sample a genotype, guaranteeing each has private variants."""
    genotypes = {}
    for index, variant in enumerate(variants):
        if ploidy == 1:
            # Rotate so SAMPLE1-only, SAMPLE2-only and shared all occur
            pattern = index % 3
            calls = {"SAMPLE1": ["1", "0", "1"][pattern],
                     "SAMPLE2": ["0", "1", "1"][pattern]}
        else:
            pattern = index % 4
            calls = {
                "SAMPLE1": ["0/1", "0/0", "1/1", "0/1"][pattern],
                "SAMPLE2": ["0/0", "0/1", "1/1", "1/1"][pattern],
            }
        genotypes[variant] = calls
    return genotypes


def build_haplotype(seqs, variants, genotypes, sample, hap_index, ploidy):
    """Apply the variants this haplotype carries, and map region coordinates."""
    out, offsets = {}, {}
    for contig, seq in seqs.items():
        applied = []
        for variant in variants:
            if variant[0] != contig:
                continue
            call = genotypes[variant][sample]
            alleles = [call] if ploidy == 1 else call.split("/")
            if alleles[hap_index] == "1":
                applied.append(variant)

        pieces, cursor, shift, marks = [], 0, 0, []
        for contig_, pos, ref, alt, _kind in applied:
            pieces.append(seq[cursor:pos])
            pieces.append(alt)
            cursor = pos + len(ref)
            shift += len(alt) - len(ref)
            marks.append((pos, shift))
        pieces.append(seq[cursor:])
        out[contig] = "".join(pieces)
        offsets[contig] = marks
    return out, offsets


def shift_for(marks, pos):
    """Cumulative length change introduced before reference position `pos`."""
    total = 0
    for mark_pos, cumulative in marks:
        if mark_pos < pos:
            total = cumulative
        else:
            break
    return total


def write_vcf(path, organism, seqs, variants, genotypes, ploidy):
    lines = [
        "##fileformat=VCFv4.2",
        f"##source=nf-testdata plant_variants.py (organism={organism}, seed={config.SEED})",
        '##INFO=<ID=KIND,Number=1,Type=String,Description="Planted variant class">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]
    for contig in seqs:
        lines.append(f"##contig=<ID={contig},length={len(seqs[contig])}>")
    lines.append("#" + "\t".join(
        ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        + config.SAMPLES))

    for contig, pos, ref, alt, kind in variants:
        calls = genotypes[(contig, pos, ref, alt, kind)]
        lines.append("\t".join([
            contig, str(pos + 1), ".", ref, alt, ".", "PASS",
            f"KIND={kind}", "GT",
        ] + [calls[s] for s in config.SAMPLES]))

    path.write_text("\n".join(lines) + "\n")


def write_regions_bed(path, spec):
    rows = [f"{c}\t{s}\t{e}\ttarget" for c, s, e in spec["regions"]]
    contig, start, end = spec["empty_region"]
    rows.append(f"{contig}\t{start}\t{end}\tno_reads_expected")
    path.write_text("\n".join(rows) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", required=True, choices=list(config.ORGANISMS))
    args = parser.parse_args()

    organism = args.organism
    spec = config.ORGANISMS[organism]
    ref_dir = Path(config.REFERENCE_DIR) / organism
    seqs = read_fasta(ref_dir / "genome.fasta")

    rng = random.Random(f"{config.SEED}:{organism}:variants")
    variants = make_variants(rng, seqs, spec)
    genotypes = assign_genotypes(rng, variants, spec["ploidy"])

    write_vcf(ref_dir / "truth.vcf", organism, seqs, variants, genotypes, spec["ploidy"])
    write_regions_bed(ref_dir / "regions.bed", spec)

    hap_dir = Path(config.CACHE_DIR) / "haplotypes"
    hap_dir.mkdir(parents=True, exist_ok=True)
    for sample in config.SAMPLES:
        for hap_index in range(spec["ploidy"]):
            hap, offsets = build_haplotype(seqs, variants, genotypes, sample,
                                           hap_index, spec["ploidy"])
            name = f"{organism}.{sample}.hap{hap_index}"
            with (hap_dir / f"{name}.fasta").open("w") as handle:
                for contig, seq in hap.items():
                    handle.write(f">{contig}\n{wrap(seq)}\n")
            with (hap_dir / f"{name}.regions").open("w") as handle:
                for contig, start, end in spec["regions"]:
                    marks = offsets[contig]
                    handle.write(f"{contig}\t{start + shift_for(marks, start)}"
                                 f"\t{end + shift_for(marks, end)}\n")
            # Indels move haplotype coordinates relative to the reference. The
            # read simulator needs to walk that back to tag methylation at real
            # reference positions, so record every step.
            with (hap_dir / f"{name}.map").open("w") as handle:
                for contig, marks in offsets.items():
                    for ref_pos, cumulative in marks:
                        handle.write(f"{contig}\t{ref_pos}\t{cumulative}\n")

    counts = {}
    for *_rest, kind in variants:
        counts[kind] = counts.get(kind, 0) + 1
    print(f"{organism}: {len(variants)} variants planted {counts}")


if __name__ == "__main__":
    main()
