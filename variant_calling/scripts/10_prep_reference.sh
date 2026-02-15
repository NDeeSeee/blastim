#!/usr/bin/env bash
set -euo pipefail

# ── 10_prep_reference.sh ────────────────────────────────────────────
# Index the reference genome with samtools faidx.
# Creates .fai index needed by freebayes.
#
# Usage:
#   bash scripts/10_prep_reference.sh
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

# ── Check tools ──
if ! command -v samtools &>/dev/null; then
    echo "ERROR: samtools not found. Activate the ngs-variants environment."
    exit 1
fi

# ── Check input ──
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference not found at $REFERENCE"
    echo "       Run scripts/01_get_data.sh first."
    exit 1
fi

echo "======================================"
echo " Step 1: Index reference (samtools faidx)"
echo "======================================"
echo ""
echo "  Reference: $REFERENCE"
echo ""

START=$(date +%s)

samtools faidx "$REFERENCE"

END=$(date +%s)
ELAPSED=$(( END - START ))

echo "==> Reference indexed in ${ELAPSED}s"
echo "    Index: ${REFERENCE}.fai"
echo ""

# ── Show reference info ──
echo "Reference contigs:"
awk '{printf "  %-30s %s bp\n", $1, $2}' "${REFERENCE}.fai"

TOTAL_BP=$(awk '{s+=$2} END{print s}' "${REFERENCE}.fai")
echo ""
echo "  Total: $TOTAL_BP bp"
echo ""
echo "Next step: bash scripts/11_prep_bam.sh"
