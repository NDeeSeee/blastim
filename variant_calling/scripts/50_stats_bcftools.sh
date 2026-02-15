#!/usr/bin/env bash
set -euo pipefail

# ── 50_stats_bcftools.sh ────────────────────────────────────────────
# Generate variant statistics with bcftools stats.
# Parses SNP/indel counts and Ti/Tv ratio.
#
# Usage:
#   bash scripts/50_stats_bcftools.sh
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

echo "======================================"
echo " Step 5: Variant statistics"
echo "======================================"
echo ""

START=$(date +%s)

# ── Stats for raw VCF ──
if [ -f "$RAW_GZ" ]; then
    echo "── Raw variants ──"
    RAW_STATS="${VARIANTS_DIR}/raw.stats.txt"
    bcftools stats "$RAW_GZ" > "$RAW_STATS"

    RAW_SNPS=$(grep "^SN" "$RAW_STATS" | grep "number of SNPs:" | awk '{print $NF}' || true)
    RAW_INDELS=$(grep "^SN" "$RAW_STATS" | grep "number of indels:" | awk '{print $NF}' || true)
    RAW_RECORDS=$(grep "^SN" "$RAW_STATS" | grep "number of records:" | awk '{print $NF}' || true)
    RAW_TITV=$(grep "^TSTV" "$RAW_STATS" | awk '{print $5}' || true)

    echo "  Records:  ${RAW_RECORDS:-0}"
    echo "  SNPs:     ${RAW_SNPS:-0}"
    echo "  Indels:   ${RAW_INDELS:-0}"
    echo "  Ti/Tv:    ${RAW_TITV:-N/A}"
    echo "  Stats:    $RAW_STATS"
    echo ""
fi

# ── Stats for filtered VCF ──
if [ -f "$FILTERED_GZ" ]; then
    echo "── Filtered variants ──"
    FILT_STATS="${VARIANTS_DIR}/filtered.stats.txt"
    bcftools stats "$FILTERED_GZ" > "$FILT_STATS"

    FILT_SNPS=$(grep "^SN" "$FILT_STATS" | grep "number of SNPs:" | awk '{print $NF}' || true)
    FILT_INDELS=$(grep "^SN" "$FILT_STATS" | grep "number of indels:" | awk '{print $NF}' || true)
    FILT_RECORDS=$(grep "^SN" "$FILT_STATS" | grep "number of records:" | awk '{print $NF}' || true)
    FILT_TITV=$(grep "^TSTV" "$FILT_STATS" | awk '{print $5}' || true)

    echo "  Records:  ${FILT_RECORDS:-0}"
    echo "  SNPs:     ${FILT_SNPS:-0}"
    echo "  Indels:   ${FILT_INDELS:-0}"
    echo "  Ti/Tv:    ${FILT_TITV:-N/A}"
    echo "  Stats:    $FILT_STATS"
    echo ""
fi

# ── Summary table ──
echo "── Summary ──"
echo ""
printf "  %-20s %-10s %-10s\n" "" "Raw" "Filtered"
printf "  %-20s %-10s %-10s\n" "────────────────────" "──────────" "──────────"
printf "  %-20s %-10s %-10s\n" "Total records" "${RAW_RECORDS:-?}" "${FILT_RECORDS:-?}"
printf "  %-20s %-10s %-10s\n" "SNPs" "${RAW_SNPS:-?}" "${FILT_SNPS:-?}"
printf "  %-20s %-10s %-10s\n" "Indels" "${RAW_INDELS:-?}" "${FILT_INDELS:-?}"
printf "  %-20s %-10s %-10s\n" "Ti/Tv" "${RAW_TITV:-?}" "${FILT_TITV:-?}"

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> Statistics finished in ${ELAPSED}s"
echo ""
echo "Next step: bash scripts/60_pick_3_variants.sh"
