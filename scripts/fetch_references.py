#!/usr/bin/env python3
"""Fetch the two source genomes the test data is built from.

  ecoli : E. coli K-12 MDS42 (GenBank AP012306.1), the complete 3.98 Mb genome.
  human : a 100 kb slice of each hg38 chromosome, plus the complete chrM.

The human slices keep real chromosome names (chr1, chr2, ... chrM) so that
contig-sorting behaviour is exercised, but their coordinates are local to the
slice. slices.tsv records where each one came from.

Downloads are cached, so re-running is cheap and offline-friendly.
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

NCBI = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        "?db=nuccore&id=AP012306.1&rettype=fasta&retmode=text")
UCSC = "https://api.genome.ucsc.edu/getData/sequence"

HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}
SLICE_LEN = 100_000
MAX_N_FRACTION = 0.01
STEP = 2_000_000


def get(url: str, cache: Path, name: str) -> bytes:
    """Fetch a URL, caching the body under `cache/name`."""
    dest = cache / name
    if dest.exists():
        return dest.read_bytes()
    # curl rather than urllib: it uses the system trust store, which saves
    # depending on certifi just to talk to NCBI and UCSC
    for attempt in range(4):
        done = subprocess.run(["curl", "-sSfL", "--max-time", "120", url],
                              capture_output=True)
        if done.returncode == 0:
            body = done.stdout
            break
        if attempt == 3:
            sys.exit(f"failed to fetch {url}: {done.stderr.decode().strip()}")
        time.sleep(2 ** attempt)
    dest.write_bytes(body)
    return body


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def fetch_ecoli(outdir: Path, cache: Path) -> None:
    raw = get(NCBI, cache, "AP012306.1.fasta").decode()
    seq = "".join(raw.splitlines()[1:]).upper()
    assert len(seq) > 3_900_000, f"MDS42 looks truncated: {len(seq)} bp"

    target = outdir / "ecoli"
    target.mkdir(parents=True, exist_ok=True)
    (target / "genome.fasta").write_text(
        f">MDS42 Escherichia coli K-12 MDS42, complete genome (AP012306.1)\n{wrap(seq)}\n")
    print(f"ecoli: MDS42 {len(seq):,} bp")


def fetch_slice(chrom: str, start: int, cache: Path) -> str:
    end = start + SLICE_LEN
    url = f"{UCSC}?genome=hg38;chrom={chrom};start={start};end={end}"
    body = json.loads(get(url, cache, f"hg38.{chrom}.{start}.json"))
    return body["dna"].upper()


def fetch_human(outdir: Path, cache: Path) -> None:
    records, provenance = [], []

    for chrom, length in HG38.items():
        start = int(length * 0.40)
        for _ in range(20):
            seq = fetch_slice(chrom, start, cache)
            if seq.count("N") / len(seq) <= MAX_N_FRACTION:
                break
            start += STEP                              # walked into a gap; step past it
        else:
            sys.exit(f"{chrom}: no gap-free {SLICE_LEN} bp window found")
        records.append((chrom, seq))
        provenance.append((chrom, start, start + SLICE_LEN))
        print(f"human: {chrom} hg38:{start:,}-{start + SLICE_LEN:,}")

    # chrM is small enough to take whole, so it is the one genuinely complete
    # human sequence in the set
    mito = fetch_slice_exact("chrM", 0, HG38_M, cache)
    records.append(("chrM", mito))
    provenance.append(("chrM", 0, HG38_M))
    print(f"human: chrM complete {len(mito):,} bp")

    target = outdir / "human"
    target.mkdir(parents=True, exist_ok=True)
    with (target / "genome.fasta").open("w") as handle:
        for chrom, seq in records:
            handle.write(f">{chrom}\n{wrap(seq)}\n")
    with (target / "slices.tsv").open("w") as handle:
        handle.write("contig\tsource\thg38_start\thg38_end\n")
        for chrom, start, end in provenance:
            handle.write(f"{chrom}\thg38:{chrom}\t{start}\t{end}\n")


HG38_M = 16569


def fetch_slice_exact(chrom: str, start: int, end: int, cache: Path) -> str:
    url = f"{UCSC}?genome=hg38;chrom={chrom};start={start};end={end}"
    return json.loads(get(url, cache, f"hg38.{chrom}.{start}.json"))["dna"].upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("data/references"))
    parser.add_argument("--cache", type=Path, default=Path(".cache/downloads"))
    args = parser.parse_args()

    args.cache.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    fetch_ecoli(args.outdir, args.cache)
    fetch_human(args.outdir, args.cache)


if __name__ == "__main__":
    main()
