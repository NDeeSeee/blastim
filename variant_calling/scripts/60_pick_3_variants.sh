#!/usr/bin/env bash
set -euo pipefail

# ── 60_pick_3_variants.sh ───────────────────────────────────────────
# Pick 3 example variants for class discussion:
#   1. Confident SNP (high QUAL, high DP)
#   2. Indel (confident)
#   3. Low-quality variant (borderline QUAL or low DP)
#
# Usage:
#   bash scripts/60_pick_3_variants.sh
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
echo " Step 6: Pick 3 example variants"
echo "======================================"
echo ""

START=$(date +%s)

EXAMPLES="${VARIANTS_DIR}/example_variants.txt"

{
    echo "# Example variants for class discussion"
    echo "# Generated: $(date)"
    echo "#"
    echo "# ── 1. Confident SNP (high QUAL, high DP) ──"
    echo ""
} > "$EXAMPLES"

# 1. Confident SNP: length(REF)==1, length(ALT)==1, QUAL>100, DP>20
CONFIDENT_SNP=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AO\t%INFO/RO\n' "$RAW_GZ" \
    | awk -F'\t' 'length($3)==1 && length($4)==1 && $5>100 && $6>20' \
    | sort -t$'\t' -k5 -rn \
    | head -1)

if [ -n "$CONFIDENT_SNP" ]; then
    echo "$CONFIDENT_SNP" >> "$EXAMPLES"
    POS=$(echo "$CONFIDENT_SNP" | cut -f2)
    echo "  1. Confident SNP at position $POS"
    echo "     $(echo "$CONFIDENT_SNP" | awk -F'\t' '{printf "REF=%s ALT=%s QUAL=%s DP=%s AO=%s RO=%s\n", $3, $4, $5, $6, $7, $8}')"
else
    echo "  1. [No confident SNP found]"
fi

{
    echo ""
    echo "# ── 2. Confident indel ──"
    echo ""
} >> "$EXAMPLES"

# 2. Confident indel: length(REF)!=length(ALT), QUAL>50, DP>10
CONFIDENT_INDEL=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AO\t%INFO/RO\n' "$RAW_GZ" \
    | awk -F'\t' '(length($3)!=1 || length($4)!=1) && (length($3)!=length($4)) && $5>50 && $6>10' \
    | sort -t$'\t' -k5 -rn \
    | head -1)

if [ -n "$CONFIDENT_INDEL" ]; then
    echo "$CONFIDENT_INDEL" >> "$EXAMPLES"
    POS=$(echo "$CONFIDENT_INDEL" | cut -f2)
    echo "  2. Confident indel at position $POS"
    echo "     $(echo "$CONFIDENT_INDEL" | awk -F'\t' '{printf "REF=%s ALT=%s QUAL=%s DP=%s AO=%s RO=%s\n", $3, $4, $5, $6, $7, $8}')"
else
    echo "  2. [No confident indel found]"
fi

{
    echo ""
    echo "# ── 3. Low-quality variant (borderline) ──"
    echo ""
} >> "$EXAMPLES"

# 3. Low-quality: QUAL between 5 and 25, or DP < 8
LOW_QUAL=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AO\t%INFO/RO\n' "$RAW_GZ" \
    | awk -F'\t' '($5>=5 && $5<=25) || ($6>0 && $6<8)' \
    | head -1)

if [ -n "$LOW_QUAL" ]; then
    echo "$LOW_QUAL" >> "$EXAMPLES"
    POS=$(echo "$LOW_QUAL" | cut -f2)
    echo "  3. Low-quality variant at position $POS"
    echo "     $(echo "$LOW_QUAL" | awk -F'\t' '{printf "REF=%s ALT=%s QUAL=%s DP=%s AO=%s RO=%s\n", $3, $4, $5, $6, $7, $8}')"
else
    echo "  3. [No low-quality variant found]"
fi

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> Picked 3 example variants in ${ELAPSED}s"
echo "    Output: $EXAMPLES"
echo ""
echo "Discussion questions:"
echo "  - Why is variant #1 confident? Look at QUAL, DP, AO/RO ratio."
echo "  - How does the indel differ from the SNP in evidence quality?"
echo "  - Would you trust variant #3? What filter thresholds would catch it?"
echo ""
echo "Workflow complete! Explore results in: $VARIANTS_DIR"
