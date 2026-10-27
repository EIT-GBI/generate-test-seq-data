# nf-testdata

Minimal, synthetic-but-real test data for the `nf-dnaseq`, `nf-dnaseq-ont` and
`nf-dnaseq-long` pipelines and the `nf-mod-*` module repos.

**You are on the `scripts` branch**, which holds the generators and this
documentation. The data lives on `main`.

Everything here is **regenerable from fixed seeds** (`./build.sh`) and, more
importantly, **its truth is known**: the variants and methylation in these
reads were planted deliberately and are written down. A test can assert that a
caller found a specific variant, rather than that an md5 did not change.

- 96.5 MB total across 129 files; largest single file 13.3 MB. No Git LFS.
- Served over `raw.githubusercontent.com`, one file at a time. Nothing clones it.
- **Pin a tag, never a branch.** A branch means regenerating this data silently
  changes what every pipeline's CI downloads.

## Using it

```groovy
params {
    testdata_base_path = 'https://raw.githubusercontent.com/eit-gbi/nf-testdata/v1.0.0'
    reference_dir      = "${params.testdata_base_path}/data/references"
    samplesheet        = "${params.testdata_base_path}/data/samplesheets/illumina/two_samples.csv"
}
```

The samplesheets already contain absolute URLs, so they work unmodified. The
`reference` column is a path relative to `reference_dir`, matching how
`modules/utils/references.nf` resolves it.

## Layout

```
data/
├── references/
│   ├── ecoli/
│   │   ├── genome.fasta                 MDS42, complete 3.98 Mb
│   │   ├── genome.fasta.fai
│   │   ├── genome.fasta.{amb,ann,bwt,pac,sa}      bwa index
│   │   ├── genome.map-ont.mmi
│   │   ├── genome.map-hifi.mmi
│   │   ├── regions.bed
│   │   ├── truth.vcf
│   │   └── truth_meth.bed
│   ├── ecoli-no-index/
│   │   └── genome.fasta                 bare: no .fai, no index
│   └── human/                           as ecoli, plus:
│       └── slices.tsv                   hg38 provenance of each contig
│
├── illumina/{ecoli,human}/
│   ├── fastq/
│   │   ├── SAMPLE1_R1.fastq.gz          SAMPLE1_R2.fastq.gz
│   │   ├── SAMPLE2_R1.fastq.gz          SAMPLE2_R2.fastq.gz
│   │   ├── SAMPLE1_L001_R1.fastq.gz     SAMPLE1_L001_R2.fastq.gz
│   │   └── SAMPLE1_L002_R1.fastq.gz     SAMPLE1_L002_R2.fastq.gz
│   ├── bam/
│   │   ├── SAMPLE1.bam  SAMPLE1.bam.bai
│   │   └── SAMPLE2.bam  SAMPLE2.bam.bai
│   └── vcf/
│       ├── SAMPLE1.vcf.gz  SAMPLE1.vcf.gz.csi
│       └── SAMPLE2.vcf.gz  SAMPLE2.vcf.gz.csi
│
├── ont/{ecoli,human}/
│   ├── fastq/   SAMPLE1.fastq.gz  SAMPLE2.fastq.gz
│   ├── ubam/    SAMPLE1.unaligned.bam           SAMPLE2.unaligned.bam
│   │            SAMPLE1.barcode01.unaligned.bam SAMPLE1.barcode02.unaligned.bam
│   └── bam/     SAMPLE1.bam + .bai              SAMPLE2.bam + .bai
│
├── pacbio/{ecoli,human}/                as ont, plus:
│   └── vcf/     SAMPLE1.vcf.gz + .csi   SAMPLE2.vcf.gz + .csi
│
└── samplesheets/{illumina,ont,pacbio}/
```

ONT has no `vcf/`: neither long-read pipeline calls small variants on ONT data
— it goes to `SAMTOOLS_CONSENSUS` and `MODKIT_PILEUP` instead.

Aligned BAMs and pre-called VCFs are shipped so that a module test for, say,
`BCFTOOLS_CONSENSUS` does not have to run an aligner first. That is most of the
value of this repo over a `conf/test.config`.

**minimap2 indexes are preset-specific.** A `.mmi` records the `-k`/`-w` it was
built with and *overrides* whatever the alignment preset would have used, so
aligning HiFi reads against a `map-ont` index silently changes the parameters.
Hence two files per genome. Use the one matching your preset, or pass the plain
fasta and let minimap2 index on the fly (~2 s at these genome sizes).

## Suggested tests, file by file

Module names below are from `nf-dnaseq`, `nf-dnaseq-ont` and `nf-dnaseq-long`.

### references/

| File | Test to write |
| --- | --- |
| `ecoli/genome.fasta` | `SAMTOOLS_FAIDX`: output `.fai` has 1 contig, length 3976195. `human/genome.fasta` gives 25 contigs with chrM at 16569 |
| `genome.fasta.fai` | `faidxExists()` returns true, so `SAMTOOLS_FAIDX` is **not** scheduled — assert via `workflow.trace.tasks()` |
| `genome.fasta.{amb,ann,bwt,pac,sa}` | `bwaIndexExists()` true branch: `BWA_MEM` runs with no `BWA_INDEX` task in the trace |
| `genome.map-ont.mmi` | `MINIMAP2_ALIGN` / `DORADO_ALIGNER` accept a prebuilt index and produce the same mapped count as aligning against the fasta |
| `genome.map-hifi.mmi` | same for `PARABRICKS_MINIMAP2` (stub) and the HiFi path of `MINIMAP2_ALIGN` |
| `regions.bed` | feed the `target` rows to a coverage assertion; the `no_reads_expected` row is the empty-coverage case for `BEDTOOLS_BIGWIG` and `BCFTOOLS_CONSENSUS` |
| `truth.vcf` | the recall test: `bcftools norm -m -any` both sides, then `isec` — expect 20 (ecoli) / 28 (human) recovered per sample, 0 missed, 0 false positives |
| `truth_meth.bed` | join `MODKIT_PILEUP` output on it; `methylated` sites ≈ 0.90, `unmethylated` ≈ 0.05 |
| `ecoli-no-index/genome.fasta` | the opposite branch: `BWA_INDEX` **and** `SAMTOOLS_FAIDX` both appear in the trace, and the run still succeeds |
| `human/slices.tsv` | documentation only — no test |

### illumina/

| File | Test to write |
| --- | --- |
| `fastq/SAMPLE1_R{1,2}.fastq.gz` | `FASTP_TRIM`: `fastp.json` reports a non-zero `adapter_cutting` count (~7% of reads carry TruSeq adapter) and the trimmed output is shorter than the input. `FASTQC_FASTQC`: 2 html + 2 zip per sample |
| `fastq/SAMPLE2_R{1,2}.fastq.gz` | the two-sample case. SAMPLE1 and SAMPLE2 carry **different genotypes**, so their VCFs must differ — a good guard against a channel that silently reuses one sample's data |
| `fastq/SAMPLE1_L00{1,2}_R{1,2}.fastq.gz` | `PREPARE_SAMPLESHEET` groups both lanes under one sample; then assert the pipeline either merges them or fails loudly, rather than overwriting `SAMPLE1.sorted.bam` twice |
| `bam/SAMPLE{1,2}.bam` | `SAMTOOLS_FLAGSTAT`: mapped ≈ **98%**, not 100%. `BEDTOOLS_BIGWIG`: a `.bw` is produced and is non-empty. `BCFTOOLS_CALL` without re-running `BWA_MEM`. `PARABRICKS_DEEPVARIANT` / `MUTECTCALLER` shape checks under `-stub-run` |
| `bam/SAMPLE{1,2}.bam.bai` | `SAMTOOLS_INDEX` idempotence: indexing an already-indexed BAM gives the same `.bai`; also lets modules take `tuple(meta, bam, bai)` directly |
| `vcf/SAMPLE{1,2}.vcf.gz` | input for `BCFTOOLS_VCF`, `BCFTOOLS_CSV` and `BCFTOOLS_CONSENSUS` without running a caller. Consensus over the `no_reads_expected` window should be all `N` |
| `vcf/SAMPLE{1,2}.vcf.gz.csi` | `BCFTOOLS_CONSENSUS` needs the index present; a test that drops it should fail cleanly |

### ont/

| File | Test to write |
| --- | --- |
| `fastq/SAMPLE{1,2}.fastq.gz` | `CHOPPER_FILTER`: with `--min_read_length 3000`, output read count drops but stays > 0 (mean length is 4 kb). `NANOPLOT_NANOPLOT`: report generated, median length in the expected band. Drives `fastq_input.csv` |
| `ubam/SAMPLE{1,2}.unaligned.bam` | `DORADO_TRIM` and `DORADO_ALIGNER` (stub for the real binary). `SAMTOOLS_FASTQ`: round-trips to the same read count. **Critically**: assert `MM`/`ML` tags survive — `samtools view` the output and check the tags are still there, since `samtools fastq` drops them without `-T` |
| `ubam/SAMPLE1.barcode0{1,2}.unaligned.bam` | the multi-file-per-sample case: two rows, one sample name. Combined read count equals `SAMPLE1.unaligned.bam` |
| `bam/SAMPLE{1,2}.bam` | `MODKIT_PILEUP` against `truth_meth.bed` (the main methylation test). `MODKIT_MOTIF_SEARCH` finds CpG as an enriched motif. `SAMTOOLS_FILTER` with `--min_read_quality` drops reads monotonically. `SAMTOOLS_CONSENSUS` with `min_depth 10` emits `N` outside the target regions |
| `bam/SAMPLE{1,2}.bam.bai` | as illumina |

### pacbio/

| File | Test to write |
| --- | --- |
| `fastq/SAMPLE{1,2}.fastq.gz` | `MINIMAP2_ALIGN` with `-x map-hifi`; `CHOPPER_FILTER` at a 2 kb floor keeps most reads (HiFi reads here average 5 kb) |
| `ubam/SAMPLE{1,2}.unaligned.bam` | `SAMTOOLS_FASTQ`, and `PARABRICKS_MINIMAP2` wiring under `-stub-run` |
| `ubam/SAMPLE1.barcode0{1,2}.unaligned.bam` | multi-file-per-sample, long-read form |
| `bam/SAMPLE{1,2}.bam` | `SAMTOOLS_FLAGSTAT`; `BCFTOOLS_CALL` on HiFi — recall against `truth.vcf` is 100%, the same as Illumina, and `./build.sh validate` enforces it |
| `vcf/SAMPLE{1,2}.vcf.gz` | cross-platform check: after `bcftools norm`, HiFi and Illumina calls on the same sample recover the identical planted set |

### Modules with no data here

`DORADO_BASECALLER` (needs POD5 and a GPU model) and the `PARABRICKS_*` modules
(need a GPU) are **stub-only**. Test them with `-stub-run`: assert the output
shapes and that the right processes were or were not scheduled. The data in this
repo cannot exercise them for real, and no amount of test data will.

## What is true about the data

| Property | Where it is written down | Verified by |
| --- | --- | --- |
| Planted SNPs and indels, per sample | `references/<org>/truth.vcf` | `./build.sh validate` |
| 5mC level at every CpG | `references/<org>/truth_meth.bed` | `./build.sh validate` |
| Windows reads were drawn from | `references/<org>/regions.bed` | — |
| A window with **no reads at all** | `regions.bed`, `no_reads_expected` | `./build.sh validate` |
| hg38 provenance of each contig | `references/human/slices.tsv` | — |

Current measured properties:

- **Variant recall is 100%, with zero false positives**, for both organisms,
  both samples and **both Illumina and PacBio**, from a plain
  `bcftools mpileup | call`. The same variants were planted once and are
  recoverable from either platform, so a cross-platform agreement test is
  meaningful.
- **Normalise before comparing.** Indels in the truth and in a caller's output
  are the same variant written differently. `bcftools norm -f <fasta> -m -any`
  on both sides turns 4 spurious misses into 0. Any test comparing to
  `truth.vcf` must do this.
- Illumina mapping rate is **~98%**, not 100%: ~2% of reads are random
  sequence that maps nowhere, so `SAMTOOLS_FLAGSTAT` has something to assert.
- **~7% of R1 reads carry TruSeq adapter**, from fragments shorter than the
  read, so `FASTP_TRIM` has real work to do.

## Methylation

The ONT uBAMs carry base-modification calls in the standard SAM `MM`/`ML`
tags, the same way dorado emits them — so `MODKIT_PILEUP` and
`MODKIT_MOTIF_SEARCH` have real input without anyone needing a GPU or a
basecalling model.

How it is built, and what you can assert:

1. `methylation_truth.py` walks every CpG inside the target regions and assigns
   each one to a **methylated** set (simulated at 90% modified) or an
   **unmethylated** set (5%), about 35/65. That decision is seeded, and written
   to `truth_meth.bed` with the level and the class for every site.
2. When a read is simulated, each CpG it covers gets a per-read modification
   call drawn against that site's level, encoded as an `ML` probability. A CpG
   is methylated on both strands, so a reverse-strand read carries the call on
   the C that pairs with the reference G — the site key shifts by one.
3. `build.sh aligned` pushes the uBAM through minimap2 with `-y`, so the tags
   survive into `data/ont/<org>/bam/`.

So the assertion a pipeline test should make is **statistical, not an md5**:
run `modkit pileup`, join its output to `truth_meth.bed`, and check that sites
labelled `methylated` come out near 0.90 and sites labelled `unmethylated` near
0.05. That survives a modkit version bump; a checksum does not.

`./build.sh validate` proves the chain end to end without modkit: it re-derives
per-site methylation straight from the aligned BAM, walking the CIGAR to map
read offsets to reference positions, and currently gets **0.89 at sites truth
calls 0.90, and 0.05 at sites truth calls 0.05**. If that check ever fails, the
`MM`/`ML` encoding is wrong and no modkit assertion downstream can be trusted.

## Samplesheets

| Sheet | What it exercises |
| --- | --- |
| `single.csv` | the simplest path |
| `two_samples.csv` | the ordinary case |
| `multi_reference.csv` | two samples on two genomes in one run — every `combine(..., by: 0)` keyed on reference name |
| `bare_reference.csv` | a reference with no index: `BWA_INDEX` and `SAMTOOLS_FAIDX` must run |
| `multiple_lanes_per_sample.csv` (illumina) | one sample, several files — the output-naming collision case |
| `multiple_barcodes_per_sample.csv` (ont, pacbio) | same, in long-read form |
| `with_platform.csv` | the optional `platform` column reaching `@RG PL` |
| `fastq_input.csv` (ont, pacbio) | FASTQ rather than uBAM entry point |
| `malformed.csv` | missing column and a nonexistent file — must fail fast and clearly |

`bare_reference.csv` points at `ecoli-no-index/`, which contains only the first
target region, so roughly a third of the reads map. That is expected: the sheet
exists to exercise index building, not alignment quality.

## Organisms

**ecoli** — *E. coli* K-12 MDS42, the complete 3.98 Mb genome, unmodified.
Reads cover three 30 kb windows; the rest of the genome is at zero depth, which
is what exercises the empty-coverage paths in `BEDTOOLS_BIGWIG` and
`BCFTOOLS_CONSENSUS`.

PacBio reads here average **5 kb**, not the 15-20 kb of real HiFi. Reads much
longer than that cover the edges of a 30 kb target window too thinly, and
variants near the edge become unrecoverable on PacBio while staying trivial on
Illumina — a difference that is an artefact of the test data, not of the
pipelines. Planted variants are also kept `VARIANT_MARGIN` (5.5 kb) clear of
window edges for the same reason.

**human** — 100 kb cut from each hg38 chromosome, plus the complete chrM.
Contigs keep real names (`chr1`, `chr2`, ... `chrM`) so contig-sorting
behaviour is exercised, but coordinates are local to each slice; `slices.tsv`
records the source. Reads cover only chr20, chr21 and chrM, so the other 22
contigs are deliberately empty. chr11 landed on pericentromeric satellite
(~27% unique 25-mers) and is kept as a deliberately hard contig.

**chrM** is the mitochondrial genome: a small circular chromosome, 16,569 bp in
hg38, present at hundreds to thousands of copies per cell rather than two. It
is included whole — the only complete, genuinely real human sequence in the set
— because it is the natural short contig for testing. It makes a run finish
fast, it is where high depth is realistic, and a reference with one 16 kb
contig alongside 100 kb ones is exactly the shape that catches contig-ordering
and chrom.sizes bugs.

## Regenerating

```sh
./build.sh                 # everything: references, truth, reads, alignments, sheets, checks
./build.sh reads aligned   # just those stages
./build.sh validate        # just re-run the correctness checks
./build.sh clean           # remove generated data, keep the download cache
```

Stages run in the order given, and the dependency order is
`references → truth → reads → aligned → sheets → validate → manifest`. There is
no incremental build: generation is ~20 s end to end, so it is always a full
rebuild.

Needs `python3`, `curl`, `samtools`, `bcftools`, `bwa`, `minimap2`. Downloads
from NCBI and the UCSC REST API are cached under `.cache/`, so re-running is
offline and fast.

All randomness derives from `SEED` in `scripts/config.py` plus a per-stream
label, so adding a platform does not shift the reads of ones already generated.

Regeneration is **byte-for-byte reproducible**: a full `./build.sh clean &&
./build.sh` produces identical checksums for all 129 files. Two things were
needed for that, and both are easy to undo by accident — gzip output is written
with `mtime=0` (`scripts/util.py`), and the `##bcftools_*` header lines, which
carry the running date, are stripped from the VCFs
(`scripts/build_alignments.sh`). `MANIFEST.md` is the record; if it changes when
nothing else did, one of those two has regressed.

## Publishing

The data is consumed over `raw.githubusercontent.com`, which is a plain file
server in front of the repo:

```
https://raw.githubusercontent.com/<owner>/<repo>/<ref>/<path>
                                                  ^ branch, tag or commit SHA
```

No clone, no history transfer, no LFS bandwidth. Nextflow's `file('https://...')`
stages these natively, so a samplesheet of URLs just works.

**The URL path is the repo path, verbatim.** There is nothing to configure and
no upload step: `data/` appears in every URL purely because `data/` is a
directory at the root of this repo.

```
main:  data/references/ecoli/genome.fasta
URL:   https://raw.githubusercontent.com/eit-gbi/nf-testdata/v1.0.0/data/references/ecoli/genome.fasta
       └─────────── BASE_URL (scripts/config.py) ──────────────────┘└──── repo path ────┘
```

`data/` sits at the root **of the `main` branch**, which is what those URLs
resolve against. On this branch it is gitignored, because it is built here and
published there — see below. `MANIFEST.md` ships alongside it and lists every
repo-relative path, which makes it the index of what is served, plus the
generator commit and seed the dataset was built from.

`./build.sh validate` checks this end of it too: that every URL in every
samplesheet starts with `BASE_URL` and resolves to a file that exists, and that
no `data/` file has been accidentally tracked on the `scripts` branch.

### Two branches

```
scripts   the generators and this README
main      the data, and nothing else
```

Each branch ignores what the other tracks, so both stay clean:

| | `scripts/.gitignore` | `main/.gitignore` |
| --- | --- | --- |
| ignores | `data/`, `MANIFEST.md` | `scripts/`, `build.sh`, `structure.md` |

To publish a new dataset:

```sh
./build.sh                       # on scripts; writes data/ (ignored here)
git switch main
git add -A
git commit -m "Dataset v1.0.1"
git tag v1.0.1
git push origin main v1.0.1
```

Set the default branch to `main` on GitHub so people land on the data.

**Switching to `scripts` deletes `data/`**, because it is tracked on `main` and
not here. That is normal git behaviour and it does not matter: `./build.sh`
rebuilds it in 20 seconds, byte for byte identical. Treat `data/` as
disposable. If the rebuild ever annoys you, check the branches out side by side
instead:

```sh
git worktree add ../nf-testdata-main main
```

Before tagging, bump `TAG` in `scripts/config.py` and re-run
`./build.sh sheets` — the tag is baked into every samplesheet URL, so the two
have to agree.

### Size and history

Each dataset version that changes bytes costs close to full size in git
history — gzipped FASTQ and BAM do not delta-compress. Measured here: the first
commit of 96.6 MB produced an 80 MB `.git`; changing `SEED` and regenerating
took it to 128 MB. So roughly **48 MB per version**, and about 19 versions
before GitHub's 1 GB threshold.

Two things keep that in hand. A no-op rebuild costs nothing, because the output
is byte-for-byte deterministic and git stores no new objects. And the
references do not churn — the 36 MB of fasta, bwa index and `.mmi` only change
if the genome does.

The thing that actually causes growth is development, not releases. Tuning
coverage or read lengths means regenerating repeatedly, so only commit to
`main` when you are tagging a version.

**If `main` does get too big**, the simplest fix is to throw its history away
and keep only the current data:

```sh
git switch main
git checkout --orphan fresh
git commit -m "Dataset v2.0.0 (earlier history dropped)"
git branch -M fresh main
git push --force origin main
```

Old tags still point at the old commits and keep their blobs alive, so delete
the ones you no longer support — that is the real decision, since a tag is what
pipelines pin to:

```sh
for tag in $(git tag --sort=-creatordate | tail -n +6); do   # keep newest 5
  git tag -d "$tag"
  git push origin --delete "$tag"
done
git reflog expire --expire=now --all
git gc --prune=now --aggressive
```

Check what is pinned first: `grep -rn "nf-testdata/v" ../nf-dnaseq* ../nf-mod-*`.
GitHub reclaims lazily, so the size shown on the repo page may lag;
`git count-objects -vH` locally is the honest number.

## Deliberate omissions

**No POD5.** `DORADO_BASECALLER` needs raw signal, which cannot be simulated
credibly, and it also needs a basecalling model (hundreds of MB, GPU-only). It
is therefore stub-only in CI regardless, exactly like the Parabricks modules.
If a real smoke test is ever wanted, subset ~100 reads from a real R10.4.1 run
with `pod5 subset` and add them under `data/ont/pod5/`; everything from
`DORADO_TRIM` onward is already driven by the uBAMs here.

**No `bedmethyl/` directory.** That would be modkit's *output*, not truth.
`references/<org>/truth_meth.bed` says what modkit should find, which is the
better thing to assert against. Generating a reference bedMethyl would require
modkit at generation time and would only ever re-state the truth file.

**E. coli CpG methylation is a convenience.** Real K-12 methylation is Dam/Dcm
(GATC adenine, CCWGG cytosine), not CpG. The CpG model is used for both
organisms so one simulator and one set of assertions covers both.

## Provenance

- *E. coli* K-12 MDS42: GenBank `AP012306.1`, public domain.
- hg38 slices: UCSC Genome Browser REST API; GRCh38 is freely redistributable.
- All reads are simulated. No real sequencing data from any individual is
  included.
