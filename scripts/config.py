"""Single source of truth for what the test dataset contains.

Every generator imports this, so changing a region or a sample here changes it
everywhere -- reads, truth files, samplesheets and README tables stay in step.
"""

REFERENCE_DIR = "data/references"
DATA_DIR = "data"
CACHE_DIR = ".cache"

# Master seed. Every generator derives its own stream from this plus a label,
# so adding a platform does not shift the reads of the ones already generated.
SEED = 20260101

ORGANISMS = {
    "ecoli": {
        "ploidy": 1,
        # Windows that reads are simulated from. Everything outside them is
        # deliberately at zero depth: that is what exercises the empty-coverage
        # paths in BEDTOOLS_BIGWIG and BCFTOOLS_CONSENSUS.
        "regions": [
            ("MDS42", 500_000, 530_000),
            ("MDS42", 1_500_000, 1_530_000),
            ("MDS42", 2_750_000, 2_780_000),
        ],
        # Named in regions.bed as a region we promise has no reads at all
        "empty_region": ("MDS42", 3_500_000, 3_550_000),
        "n_snps": 24,
        "n_indels": 6,
    },
    "human": {
        "ploidy": 2,
        "regions": [
            ("chr20", 0, 30_000),
            ("chr21", 0, 30_000),
            ("chrM", 0, 16_569),
        ],
        # chr11 is pericentromeric satellite: ~27% unique 25-mers. Kept as a
        # deliberately hard contig, and never given reads.
        "empty_region": ("chr11", 0, 100_000),
        "n_snps": 30,
        "n_indels": 8,
    },
}

SAMPLES = ["SAMPLE1", "SAMPLE2"]

# Keep planted variants this far from the edge of a target window. Long reads
# are sampled from a window widened by one read length, but a window starting
# at coordinate 0 has nothing to widen into, so its first few kb are thin. Must
# exceed the longest platform's mean read length.
VARIANT_MARGIN = 5_500

# Read simulation. Coverage is per target region; everything else is untouched.
PLATFORMS = {
    "illumina": {
        "read_length": 150,
        "insert_mean": 350,
        "insert_sd": 50,
        "coverage": 25,
        "error_rate": 0.002,
        "qual_mean": 36,
        # Some fragments are shorter than the read, so the read runs off the
        # end into the adapter. Without these fastp has nothing to trim and
        # FASTP_TRIM is untested.
        "short_fragment_fraction": 0.08,
        # A little unmappable contamination, so SAMTOOLS_FLAGSTAT has a
        # mapping rate to assert rather than a flat 100%
        "unmappable_fraction": 0.02,
    },
    "ont": {
        "length_mean": 4_000,
        "length_sd": 2_500,
        "length_min": 500,
        "coverage": 12,
        "error_rate": 0.05,
        "qual_mean": 16,
        "unmappable_fraction": 0.02,
    },
    "pacbio": {
        # 5 kb rather than a realistic 15-20 kb: with 30 kb target windows,
        # reads much longer than this cover the window edges too thinly and
        # variants there become unrecoverable on PacBio but not Illumina.
        "length_mean": 5_000,
        "length_sd": 1_200,
        "length_min": 2_000,
        "coverage": 12,
        "error_rate": 0.001,
        "qual_mean": 33,
        "unmappable_fraction": 0.02,
    },
}

# ONT only: CpGs are split into a methylated and an unmethylated set so the
# modkit pileup has something to recover. truth_meth.bed records which is which.
METH = {
    "fraction_methylated_sites": 0.35,
    "methylated_level": 0.90,
    "unmethylated_level": 0.05,
}

# ONT barcodes, for the "one sample across several barcodes" samplesheet
BARCODES = ["barcode01", "barcode02"]

# Where the data is served from once published. Samplesheets are written with
# these URLs baked in, so they work with no rewriting. Pin a TAG, never a
# branch: a branch means regenerating the data silently changes what every
# pipeline's CI downloads.
REPO = "eit-gbi/nf-testdata"
TAG = "v1.0.0"
BASE_URL = f"https://raw.githubusercontent.com/{REPO}/{TAG}"
