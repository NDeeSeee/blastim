#!/usr/bin/env bash
set -euo pipefail

# ── run_all.sh ──────────────────────────────────────────────────────
# Run the complete variant calling workflow (steps 10→60).
#
# Usage:
#   bash workflow/run_all.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# ── Check config ──
CONFIG="${REPO_ROOT}/workflow/config.sh"
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file not found."
    echo "       cp workflow/config.sh.example workflow/config.sh"
    exit 1
fi

# Set a timestamped OUTDIR (respects pre-existing export if set).
export OUTDIR="${OUTDIR:-outputs/run_$(date +%Y%m%d_%H%M%S)}"

echo "╔════════════════════════════════════════════════╗"
echo "║  Variant Calling Workflow                      ║"
echo "╚════════════════════════════════════════════════╝"
echo ""
echo "Config:  $CONFIG"
echo "Output:  $OUTDIR"
echo "Start:   $(date)"
echo ""

TOTAL_START=$(date +%s)

# ── Step 1: Prep reference ──
echo "━━━ Running Step 10: Prep reference ━━━"
bash "${REPO_ROOT}/scripts/10_prep_reference.sh"
echo ""

# ── Step 1b: Validate BAM ──
echo "━━━ Running Step 11: Validate BAM ━━━"
bash "${REPO_ROOT}/scripts/11_prep_bam.sh"
echo ""

# ── Step 2: Call variants ──
echo "━━━ Running Step 20: freebayes ━━━"
bash "${REPO_ROOT}/scripts/20_call_freebayes.sh"
echo ""

# ── Step 3: Compress + index ──
echo "━━━ Running Step 30: bgzip + tabix ━━━"
bash "${REPO_ROOT}/scripts/30_compress_index_vcf.sh"
echo ""

# ── Step 4: Filter ──
echo "━━━ Running Step 40: bcftools filter ━━━"
bash "${REPO_ROOT}/scripts/40_filter_bcftools.sh"
echo ""

# ── Step 5: Statistics ──
echo "━━━ Running Step 50: bcftools stats ━━━"
bash "${REPO_ROOT}/scripts/50_stats_bcftools.sh"
echo ""

# ── Step 6: Pick examples ──
echo "━━━ Running Step 60: Pick 3 examples ━━━"
bash "${REPO_ROOT}/scripts/60_pick_3_variants.sh"
echo ""

TOTAL_END=$(date +%s)
TOTAL_ELAPSED=$(( TOTAL_END - TOTAL_START ))
MINUTES=$(( TOTAL_ELAPSED / 60 ))
SECONDS_R=$(( TOTAL_ELAPSED % 60 ))

echo "╔════════════════════════════════════════════════╗"
echo "║  Workflow complete!                            ║"
printf "║  Total time: %-34s║\n" "${MINUTES}m ${SECONDS_R}s"
echo "╚════════════════════════════════════════════════╝"
echo ""
echo "Output directory: ${REPO_ROOT}/${OUTDIR}"
echo ""
echo "Key files to explore:"
echo "  ${OUTDIR}/variants/raw.vcf.gz          - All called variants"
echo "  ${OUTDIR}/variants/filtered.vcf.gz     - Filtered variants"
echo "  ${OUTDIR}/variants/raw.stats.txt       - Raw variant statistics"
echo "  ${OUTDIR}/variants/filtered.stats.txt  - Filtered variant statistics"
echo "  ${OUTDIR}/variants/example_variants.txt - 3 example variants"
