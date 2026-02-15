#!/usr/bin/env bash
set -euo pipefail

# ── 20_call_freebayes.sh ────────────────────────────────────────────
# Call variants with freebayes.
# Produces raw.vcf with rich INFO fields (AO, RO, SRP, SAP, etc.).
#
# Usage:
#   bash scripts/20_call_freebayes.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# ── Load config ──
CONFIG="${REPO_ROOT}/workflow/config.sh"
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config not found at $CONFIG"
    echo "       cp workflow/config.sh.example workflow/config.sh"
    exit 1
fi
source "$CONFIG"

# ── Resolve paths relative to repo root ──
REFERENCE="${REPO_ROOT}/${REFERENCE}"
BAM_FILE="${REPO_ROOT}/${BAM_FILE}"
OUTDIR="${REPO_ROOT}/${OUTDIR}"

# ── Check tools ──
if ! command -v freebayes &>/dev/null; then
    echo "ERROR: freebayes not found. Activate the ngs-variants environment."
    exit 1
fi

# ── Check inputs ──
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference not found at $REFERENCE"
    exit 1
fi

if [ ! -f "${REFERENCE}.fai" ]; then
    echo "ERROR: Reference index (.fai) not found."
    echo "       Run scripts/10_prep_reference.sh first."
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    echo "ERROR: BAM file not found at $BAM_FILE"
    exit 1
fi

VARIANTS_DIR="${OUTDIR}/variants"
mkdir -p "$VARIANTS_DIR"
RAW_VCF="${VARIANTS_DIR}/raw.vcf"

echo "======================================"
echo " Step 2: Variant calling (freebayes)"
echo "======================================"
echo ""
echo "  Reference: $REFERENCE"
echo "  BAM:       $BAM_FILE"
echo "  Ploidy:    ${PLOIDY}"
echo "  Output:    $RAW_VCF"
echo ""

START=$(date +%s)

freebayes \
    --fasta-reference "$REFERENCE" \
    --ploidy "${PLOIDY}" \
    --bam "$BAM_FILE" \
    > "$RAW_VCF"

END=$(date +%s)
ELAPSED=$(( END - START ))

# ── Quick summary ──
TOTAL=$(grep -cv "^#" "$RAW_VCF" || true)
TOTAL=${TOTAL:-0}
SNPS=$( (grep -v "^#" "$RAW_VCF" || true) | awk -F'\t' 'length($4)==1 && length($5)==1' | wc -l | tr -d ' ')
INDELS=$(( TOTAL - SNPS ))

echo ""
echo "==> freebayes finished in ${ELAPSED}s"
echo ""
echo "  Total raw variants: $TOTAL"
echo "  SNPs:               $SNPS"
echo "  Indels:             $INDELS"
echo ""
echo "  Output: $RAW_VCF"
echo ""
echo "Next step: bash scripts/30_compress_index_vcf.sh"
