#!/usr/bin/env python3
"""Decide which CpGs are methylated, and record it.

Writes data/references/<org>/truth_meth.bed: every CpG inside a target region,
with the methylation level the ONT reads will be simulated at. MODKIT_PILEUP
should recover these fractions within sampling error, which is a far better
assertion than an md5 of a bedMethyl file.

Note for the E. coli set: real K-12 methylation is Dam/Dcm (GATC adenine and
CCWGG cytosine), not CpG. The CpG model here is a convenience so the same
simulator and the same assertions work for both organisms.
"""

import argparse
import random
from pathlib import Path

import config
from plant_variants import read_fasta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", required=True, choices=list(config.ORGANISMS))
    args = parser.parse_args()

    spec = config.ORGANISMS[args.organism]
    ref_dir = Path(config.REFERENCE_DIR) / args.organism
    seqs = read_fasta(ref_dir / "genome.fasta")
    rng = random.Random(f"{config.SEED}:{args.organism}:methylation")

    rows, methylated = [], 0
    for contig, start, end in spec["regions"]:
        seq = seqs[contig]
        for pos in range(start, min(end, len(seq)) - 1):
            if seq[pos:pos + 2] != "CG":
                continue
            is_methylated = rng.random() < config.METH["fraction_methylated_sites"]
            level = (config.METH["methylated_level"] if is_methylated
                     else config.METH["unmethylated_level"])
            methylated += is_methylated
            rows.append(f"{contig}\t{pos}\t{pos + 1}\t5mC\t{level:.2f}\t+\t"
                        f"{'methylated' if is_methylated else 'unmethylated'}")

    out = ref_dir / "truth_meth.bed"
    out.write_text("\n".join(rows) + "\n")
    print(f"{args.organism}: {len(rows)} CpG sites, {methylated} methylated "
          f"-> {out}")


if __name__ == "__main__":
    main()
