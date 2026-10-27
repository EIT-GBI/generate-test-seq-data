# nf-testdata

Test data for the `nf-dnaseq`, `nf-dnaseq-ont` and `nf-dnaseq-long` pipelines
and the `nf-mod-*` module repos.

**This branch holds only the data.** The generators that produce it, and the
full documentation of what is in here and what is true about it, live on the
[`scripts` branch](../../tree/scripts).

Every file is served directly over `raw.githubusercontent.com`, at its own path
in this repo:

```
https://raw.githubusercontent.com/OWNER/REPO/TAG/data/references/ecoli/genome.fasta
```

Pin a **tag**, never a branch — `main` is force-pushed on every release, and a
branch ref would silently change what your CI downloads.

```groovy
params {
    testdata_base_path = 'https://raw.githubusercontent.com/OWNER/REPO/TAG'
    reference_dir      = "${params.testdata_base_path}/data/references"
    samplesheet        = "${params.testdata_base_path}/data/samplesheets/illumina/two_samples.csv"
}
```

`MANIFEST.md` lists every file with its size and sha256, and records the
generator commit and seed this dataset was built from.

## Contents

| Path | What |
| --- | --- |
| `data/references/` | E. coli MDS42 and a 25-contig human reference, with indexes, plus `truth.vcf` and `truth_meth.bed` |
| `data/illumina/` | paired FASTQ, aligned BAM, pre-called VCF |
| `data/ont/` | FASTQ, unaligned BAM with MM/ML methylation tags, aligned BAM |
| `data/pacbio/` | FASTQ, unaligned BAM, aligned BAM, pre-called VCF |
| `data/samplesheets/` | 23 sheets, each exercising one thing |

## History

Each release is a **parentless commit** on `main`, tagged. Old versions are
independent, so they can be pruned by deleting their tags without rewriting
anything. See the `scripts` branch README for the retention recipe.
