#!/usr/bin/env bash
set -euo pipefail

# ── 11_prep_bam.sh ──────────────────────────────────────────────────
# Validate and inspect the BAM file before variant calling.
# Runs quickcheck, flagstat, idxstats, and depth summary.
#
# Usage:
#   bash scripts/11_prep_bam.sh
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
BAM_FILE="${REPO_ROOT}/${BAM_FILE}"

# ── Check tools ──
if ! command -v samtools &>/dev/null; then
    echo "ERROR: samtools not found. Activate the ngs-variants environment."
    exit 1
fi

# ── Check input ──
if [ ! -f "$BAM_FILE" ]; then
    echo "ERROR: BAM file not found at $BAM_FILE"
    echo "       Run scripts/01_get_data.sh first."
    exit 1
fi

echo "======================================"
echo " Step 1b: Validate BAM file"
echo "======================================"
echo ""
echo "  BAM: $BAM_FILE"
echo ""

START=$(date +%s)

# ── quickcheck ──
echo "── samtools quickcheck ──"
if samtools quickcheck "$BAM_FILE"; then
    echo "  [OK] BAM file passes quickcheck"
else
    echo "  [FAIL] BAM file is corrupt or truncated"
    exit 1
fi

# ── Check index ──
echo ""
if [ -f "${BAM_FILE}.bai" ] || [ -f "${BAM_FILE%.bam}.bai" ]; then
    echo "  [OK] BAM index found"
else
    echo "  [WARN] BAM index not found, creating..."
    samtools index "$BAM_FILE"
fi

# ── flagstat ──
echo ""
echo "── samtools flagstat ──"
samtools flagstat "$BAM_FILE"

# ── idxstats ──
echo ""
echo "── samtools idxstats ──"
samtools idxstats "$BAM_FILE" | column -t

# ── depth summary ──
echo ""
echo "── Coverage summary (samtools depth) ──"
MEAN_DEPTH=$(samtools depth -a "$BAM_FILE" | awk '{s+=$3; n++} END{if(n>0) printf "%.1f", s/n; else print "0"}')
echo "  Mean depth: ${MEAN_DEPTH}x"

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> BAM validation finished in ${ELAPSED}s"
echo ""
echo "Next step: bash scripts/20_call_freebayes.sh"
