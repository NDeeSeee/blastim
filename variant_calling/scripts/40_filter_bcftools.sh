#!/usr/bin/env bash
set -euo pipefail

# ── 40_filter_bcftools.sh ───────────────────────────────────────────
# Filter variants with bcftools using QUAL and DP thresholds.
# Produces filtered.vcf.gz (only PASS variants).
#
# Usage:
#   bash scripts/40_filter_bcftools.sh
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
RAW_GZ="${VARIANTS_DIR}/raw.vcf.gz"
FILTERED_GZ="${VARIANTS_DIR}/filtered.vcf.gz"

# ── Check tools ──
if ! command -v bcftools &>/dev/null; then
    echo "ERROR: bcftools not found. Activate the ngs-variants environment."
    exit 1
fi

# ── Check input ──
if [ ! -f "$RAW_GZ" ]; then
    echo "ERROR: raw.vcf.gz not found at $RAW_GZ"
    echo "       Run scripts/30_compress_index_vcf.sh first."
    exit 1
fi

echo "======================================"
echo " Step 4: Filter variants (bcftools)"
echo "======================================"
echo ""
echo "  Input:    $RAW_GZ"
echo "  Output:   $FILTERED_GZ"
echo "  Filters:  QUAL >= ${MIN_QUAL}, DP >= ${MIN_DP}, DP <= ${MAX_DP}"
echo ""

START=$(date +%s)

BEFORE=$(bcftools view -H "$RAW_GZ" | wc -l | tr -d ' ')

bcftools view -i \
    "QUAL >= ${MIN_QUAL} && INFO/DP >= ${MIN_DP} && INFO/DP <= ${MAX_DP}" \
    "$RAW_GZ" \
    -Oz -o "$FILTERED_GZ"

tabix -p vcf "$FILTERED_GZ"

AFTER=$(bcftools view -H "$FILTERED_GZ" | wc -l | tr -d ' ')

END=$(date +%s)
ELAPSED=$(( END - START ))

PASS_RATE=$(awk "BEGIN{if($BEFORE>0) printf \"%.1f\", ($AFTER/$BEFORE)*100; else print \"0\"}")

echo ""
echo "==> Filtering finished in ${ELAPSED}s"
echo ""
echo "  Before filter: $BEFORE variants"
echo "  After filter:  $AFTER variants"
echo "  Pass rate:     ${PASS_RATE}%"
echo ""
echo "  Output: $FILTERED_GZ"
echo ""
echo "Next step: bash scripts/50_stats_bcftools.sh"
