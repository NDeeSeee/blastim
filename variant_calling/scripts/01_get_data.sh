#!/usr/bin/env bash
set -euo pipefail

# ── 01_get_data.sh ──────────────────────────────────────────────────
# Download B. subtilis 168 reference genome from NCBI,
# simulate reads with wgsim (~50x coverage),
# align with bwa mem, sort and index the BAM.
#
# This is an INSTRUCTOR PREP script — run once before class.
#
# Usage:
#   bash scripts/01_get_data.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

INPUT_DIR="${REPO_ROOT}/data/input"
TRUTH_DIR="${REPO_ROOT}/data/truth"
mkdir -p "$INPUT_DIR" "$TRUTH_DIR"

REF="${INPUT_DIR}/reference.fasta"
BAM="${INPUT_DIR}/reads_sorted.bam"

# ── wgsim parameters ──
SEED=42
N_READS=700000    # ~50x coverage for ~4.2 Mb genome (700k pairs * 150bp * 2 / 4.2M)
READ_LEN=150
MUT_RATE=0.001    # base substitution rate
INDEL_FRAC=0.15   # fraction of mutations that are indels
ERROR_RATE=0.005  # sequencing error rate

# ── Check tools ──
for tool in wget samtools bwa wgsim; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool not found. Activate the ngs-variants environment."
        exit 1
    fi
done

echo "======================================"
echo " Data preparation (instructor)"
echo "======================================"
echo ""

START=$(date +%s)

# ── Step 1: Download reference genome ──
if [ -f "$REF" ]; then
    echo "==> Reference already exists: $REF"
else
    echo "==> Downloading B. subtilis 168 reference genome..."
    wget -q -O "${REF}.gz" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz"
    gunzip "${REF}.gz"
    echo "    Downloaded: $REF"
fi

# ── Step 2: Index reference for bwa ──
echo "==> Indexing reference for bwa..."
bwa index "$REF"

# ── Step 3: Simulate reads with wgsim ──
R1="${INPUT_DIR}/sim_R1.fq"
R2="${INPUT_DIR}/sim_R2.fq"

echo "==> Simulating reads with wgsim (seed=$SEED, N=$N_READS, error=$ERROR_RATE)..."
wgsim \
    -S "$SEED" \
    -N "$N_READS" \
    -1 "$READ_LEN" \
    -2 "$READ_LEN" \
    -r "$MUT_RATE" \
    -R "$INDEL_FRAC" \
    -e "$ERROR_RATE" \
    "$REF" \
    "$R1" "$R2" \
    > "${TRUTH_DIR}/wgsim.mutations.txt"

echo "    Reads: $R1, $R2"
echo "    Truth set: ${TRUTH_DIR}/wgsim.mutations.txt"

# ── Step 4: Align reads ──
echo "==> Aligning reads with bwa mem..."
bwa mem -t 4 -R '@RG\tID:sim\tSM:bsub168\tPL:ILLUMINA' "$REF" "$R1" "$R2" \
    | samtools sort -@ 4 -o "$BAM"

# ── Step 5: Index BAM ──
echo "==> Indexing BAM..."
samtools index "$BAM"

# ── Step 6: Clean up FASTQ ──
echo "==> Removing temporary FASTQ files..."
rm -f "$R1" "$R2"

# ── Step 7: Quick stats ──
echo ""
echo "==> Alignment summary:"
samtools flagstat "$BAM"

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> Data preparation finished in ${ELAPSED}s"
echo ""
echo "Key files:"
echo "  Reference:   $REF"
echo "  Sorted BAM:  $BAM"
echo "  BAM index:   ${BAM}.bai"
echo "  Truth set:   ${TRUTH_DIR}/wgsim.mutations.txt"
echo ""
echo "Next step: cp workflow/config.sh.example workflow/config.sh"
echo "           bash scripts/10_prep_reference.sh"
