#!/usr/bin/env bash
set -euo pipefail

# ── 30_compress_index_vcf.sh ────────────────────────────────────────
# Compress VCF with bgzip and index with tabix.
#
# Teaching moment: bgzip (block-gzip) is NOT the same as gzip.
# Regular gzip does not support random access; bgzip does, and
# tabix requires bgzip-compressed files.
#
# Usage:
#   bash scripts/30_compress_index_vcf.sh
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
OUTDIR="${REPO_ROOT}/${OUTDIR}"
VARIANTS_DIR="${OUTDIR}/variants"
RAW_VCF="${VARIANTS_DIR}/raw.vcf"

# ── Check tools ──
for tool in bgzip tabix; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool not found. Activate the ngs-variants environment."
        exit 1
    fi
done

# ── Check input ──
if [ ! -f "$RAW_VCF" ]; then
    echo "ERROR: raw.vcf not found at $RAW_VCF"
    echo "       Run scripts/20_call_freebayes.sh first."
    exit 1
fi

echo "======================================"
echo " Step 3: Compress + index VCF"
echo "======================================"
echo ""
echo "  Input:  $RAW_VCF"
echo ""

START=$(date +%s)

# ── bgzip (replaces raw.vcf with raw.vcf.gz) ──
echo "==> Compressing with bgzip..."
bgzip -f "$RAW_VCF"
echo "    Created: ${RAW_VCF}.gz"

# ── tabix ──
echo "==> Indexing with tabix..."
tabix -p vcf "${RAW_VCF}.gz"
echo "    Created: ${RAW_VCF}.gz.tbi"

END=$(date +%s)
ELAPSED=$(( END - START ))

# ── Compare sizes ──
GZ_SIZE=$(du -h "${RAW_VCF}.gz" | cut -f1)

echo ""
echo "==> Compress + index finished in ${ELAPSED}s"
echo ""
echo "  Compressed VCF: ${RAW_VCF}.gz ($GZ_SIZE)"
echo "  Tabix index:    ${RAW_VCF}.gz.tbi"
echo ""
echo "  Note: bgzip != gzip. bgzip creates block-compressed files"
echo "        that support random access via tabix."
echo ""
echo "Next step: bash scripts/40_filter_bcftools.sh"
