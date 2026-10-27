#!/usr/bin/env bash
# Regenerate the whole test dataset from fixed seeds.
#
#   ./build.sh                 every stage, in order
#   ./build.sh reads aligned   just those stages
#   ./build.sh clean           remove generated data (keeps the download cache)
#
# Stages, in dependency order:
#   references  fetch the two source genomes, build bwa/faidx/minimap2 indexes
#   truth       plant variants and methylation, write the truth files
#   reads       simulate FASTQ and unaligned BAM for all three platforms
#   aligned     align those reads, producing BAM and pre-called VCF
#   sheets      write the samplesheets
#   validate    check the data has the properties the tests will assert
#   manifest    write MANIFEST.md and enforce the size budget
#
# Needs: python3, curl, samtools, bcftools, bwa, minimap2
set -euo pipefail
cd "$(dirname "$0")"

PY=${PY:-python3}
REF=data/references
ORGANISMS=(ecoli human)
SAMPLES=(SAMPLE1 SAMPLE2)
# minimap2 indexes are preset-specific: the index records -k/-w and overrides
# whatever the alignment preset would have used, so ONT and HiFi need separate
# files. Aligning against the wrong one silently changes the parameters.
PRESETS=(map-ont map-hifi)

say() { printf '\n=== %s ===\n' "$1"; }

require() {
  local missing=0
  for tool in "$@"; do
    command -v "$tool" >/dev/null || { echo "missing required tool: $tool" >&2; missing=1; }
  done
  [[ $missing -eq 0 ]] || exit 1
}

stage_references() {
  say references
  $PY scripts/fetch_references.py --outdir "$REF" --cache .cache/downloads

  # A bare copy of one E. coli target region, with no index beside it. This is
  # what exercises the BWA_INDEX / SAMTOOLS_FAIDX branches: the pipeline only
  # builds an index when one is missing.
  mkdir -p "$REF/ecoli-no-index"
  samtools faidx "$REF/ecoli/genome.fasta" MDS42:500001-530000 \
    | sed '1s/.*/>MDS42_region1/' > "$REF/ecoli-no-index/genome.fasta"

  for org in "${ORGANISMS[@]}"; do
    samtools faidx "$REF/$org/genome.fasta"
    bwa index "$REF/$org/genome.fasta" 2>/dev/null
    for preset in "${PRESETS[@]}"; do
      minimap2 -x "$preset" -d "$REF/$org/genome.$preset.mmi" \
        "$REF/$org/genome.fasta" 2>/dev/null
    done
    echo "indexed $org"
  done
}

stage_truth() {
  say truth
  for org in "${ORGANISMS[@]}"; do
    $PY scripts/plant_variants.py --organism "$org"
    $PY scripts/methylation_truth.py --organism "$org"
  done
}

stage_reads() {
  say reads
  for org in "${ORGANISMS[@]}"; do
    for sample in "${SAMPLES[@]}"; do
      $PY scripts/simulate_reads.py --organism "$org" --platform illumina \
        --sample "$sample" --outdir "data/illumina/$org/fastq"
      for plat in ont pacbio; do
        mkdir -p "data/$plat/$org/fastq" "data/$plat/$org/ubam"
        $PY scripts/simulate_reads.py --organism "$org" --platform "$plat" \
          --sample "$sample" --outdir "data/$plat/$org/fastq" \
          --sam ".cache/$plat.$org.$sample.sam"
        samtools view -b --no-PG \
          -o "data/$plat/$org/ubam/$sample.unaligned.bam" \
          ".cache/$plat.$org.$sample.sam"
      done
    done
  done
  # One sample split across two barcodes / two lanes, for the samplesheets
  # that put several files under a single sample name.
  $PY scripts/split_barcodes.py
}

stage_aligned() { say aligned; bash scripts/build_alignments.sh; }
stage_sheets()  { say sheets;  $PY scripts/make_samplesheets.py; }
stage_validate() { say validate; $PY scripts/validate.py; }
stage_manifest() { say manifest; $PY scripts/make_manifest.py; }

stage_clean() {
  say clean
  rm -rf data/illumina data/ont data/pacbio data/samplesheets data/references
  rm -rf .cache/haplotypes .cache/validate
  rm -f .cache/*.sam
  echo "generated data removed (.cache/downloads kept)"
}

ALL_STAGES=(references truth reads aligned sheets validate manifest)
stages=("${@:-}")
[[ -z "${stages[0]:-}" ]] && stages=("${ALL_STAGES[@]}")

require python3 curl samtools bcftools bwa minimap2
mkdir -p .cache/downloads

for stage in "${stages[@]}"; do
  if ! declare -F "stage_$stage" >/dev/null; then
    echo "unknown stage: $stage" >&2
    echo "available: ${ALL_STAGES[*]} clean" >&2
    exit 1
  fi
  "stage_$stage"
done

echo
echo "done: ${stages[*]}"
