#!/usr/bin/env bash
# Align the simulated reads and pre-call variants.
#
# Module-level tests should not have to run an aligner to get a BAM, so the
# aligned BAMs and a bcftools VCF are shipped alongside the reads.
set -euo pipefail

REF=data/references

# bcftools stamps the running date and the invoking command into the VCF
# header, so an otherwise identical rebuild gets a different checksum. Strip
# those lines; everything they record is already in this script.
call_vcf() {
  local fasta=$1 ploidy=$2 bam=$3 out=$4
  bcftools mpileup -Ou -f "$fasta" -a AD,DP "$bam" 2>/dev/null \
    | bcftools call --ploidy "$ploidy" -mv -Ov 2>/dev/null \
    | grep -vE '^##(bcftools_|samtools)' \
    | bgzip > "$out"
  bcftools index -f "$out"
}
ORGANISMS=(ecoli human)
SAMPLES=(SAMPLE1 SAMPLE2)

for org in "${ORGANISMS[@]}"; do
  fasta="$REF/$org/genome.fasta"
  ploidy=1
  [[ "$org" == "human" ]] && ploidy=2

  for sample in "${SAMPLES[@]}"; do
    # --- Illumina -----------------------------------------------------------
    mkdir -p "data/illumina/$org/bam" "data/illumina/$org/vcf"
    bam="data/illumina/$org/bam/$sample.bam"
    bwa mem -t 4 -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" "$fasta" \
        "data/illumina/$org/fastq/${sample}_R1.fastq.gz" \
        "data/illumina/$org/fastq/${sample}_R2.fastq.gz" 2>/dev/null \
      | samtools sort --no-PG -o "$bam" -
    samtools index "$bam"

    call_vcf "$fasta" "$ploidy" "$bam" "data/illumina/$org/vcf/$sample.vcf.gz"

    # --- ONT: carry the MM/ML methylation tags through the aligner ----------
    mkdir -p "data/ont/$org/bam"
    ont_bam="data/ont/$org/bam/$sample.bam"
    samtools fastq -T MM,ML "data/ont/$org/ubam/$sample.unaligned.bam" 2>/dev/null \
      | minimap2 -ax map-ont -y -t 4 "$fasta" - 2>/dev/null \
      | samtools sort --no-PG -o "$ont_bam" -
    samtools index "$ont_bam"

    # --- PacBio HiFi --------------------------------------------------------
    mkdir -p "data/pacbio/$org/bam" "data/pacbio/$org/vcf"
    pb_bam="data/pacbio/$org/bam/$sample.bam"
    samtools fastq "data/pacbio/$org/ubam/$sample.unaligned.bam" 2>/dev/null \
      | minimap2 -ax map-hifi -t 4 "$fasta" - 2>/dev/null \
      | samtools sort --no-PG -o "$pb_bam" -
    samtools index "$pb_bam"

    call_vcf "$fasta" "$ploidy" "$pb_bam" "data/pacbio/$org/vcf/$sample.vcf.gz"

    echo "aligned: $org/$sample"
  done
done
