# Layout

What `./build.sh` actually produces. 129 files, 92.5 MB.

```
nf-testdata/
├── README.md                     what each file is and what is TRUE about it
├── MANIFEST.md                   every file, size and sha256 (the determinism tripwire)
├── build.sh                      seeded regeneration of everything
├── scripts/
└── data/
    ├── references/               <- params.reference_dir for all three pipelines
    │   ├── ecoli/                E. coli K-12 MDS42, complete 3.98 Mb (AP012306.1)
    │   │   ├── genome.fasta
    │   │   ├── genome.fasta.fai
    │   │   ├── genome.fasta.{amb,ann,bwt,pac,sa}    bwa index
    │   │   ├── genome.map-ont.mmi                   minimap2, ONT preset
    │   │   ├── genome.map-hifi.mmi                  minimap2, HiFi preset
    │   │   ├── truth.vcf          planted SNPs + indels, genotyped per sample
    │   │   ├── truth_meth.bed     every CpG in the target regions + its 5mC level
    │   │   └── regions.bed        target windows + one promised to have no reads
    │   ├── ecoli-no-index/
    │   │   └── genome.fasta       BARE: forces BWA_INDEX + SAMTOOLS_FAIDX
    │   └── human/                 100 kb of each hg38 chromosome + complete chrM
    │       ├── genome.fasta  .fai  bwa index  .map-ont.mmi  .map-hifi.mmi
    │       ├── truth.vcf  truth_meth.bed  regions.bed
    │       └── slices.tsv         where in hg38 each contig was cut from
    ├── illumina/
    │   ├── ecoli/{fastq,bam,vcf}
    │   └── human/{fastq,bam,vcf}
    ├── ont/
    │   ├── ecoli/{fastq,ubam,bam}   uBAM carries MM/ML 5mC tags
    │   └── human/{fastq,ubam,bam}
    ├── pacbio/
    │   ├── ecoli/{fastq,ubam,bam,vcf}
    │   └── human/{fastq,ubam,bam,vcf}
    └── samplesheets/
        ├── illumina/  single, two_samples, multi_reference, bare_reference,
        │              multiple_lanes_per_sample, with_platform, malformed
        ├── ont/       single, two_samples, multi_reference, bare_reference,
        │              multiple_barcodes_per_sample, with_platform, fastq_input, malformed
        └── pacbio/    (same as ont)
```

## Differences from the first sketch

- **`human_chr20/` became `human/`** — a slice of every chromosome, not one.
  Multi-contig is the point: it exercises contig sorting and the per-reference
  channel joins that a single contig never would.
- **No `ont/pod5/`.** Raw signal cannot be simulated credibly, and
  `DORADO_BASECALLER` also needs a GPU-only model, so it is stub-only in CI
  regardless. Everything from `DORADO_TRIM` onward runs off the uBAMs.
- **No `bedmethyl/`.** That is modkit's *output*, not truth. `truth_meth.bed`
  says what modkit should find, which is the better assertion target and needs
  no modkit at generation time.
- **`multiple_barcodes_per_sample` is long-read only.** Illumina has lanes, not
  barcodes, so the equivalent sheet is `multiple_lanes_per_sample.csv`.
- **`bare_reference.csv` added**, to drive `ecoli-no-index/`.
- **`fastq_input.csv` added** for the long-read pipelines, which can enter at
  minimap2/chopper with FASTQ rather than at dorado with uBAM.
