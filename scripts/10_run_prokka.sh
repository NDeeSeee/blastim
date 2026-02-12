#!/usr/bin/env bash
set -euo pipefail

# ── 10_run_prokka.sh ────────────────────────────────────────────────
# Structural annotation with Prokka.
# Predicts CDS, rRNA, tRNA and produces GFF, FAA, FFN, TSV, etc.
#
# Usage:
#   bash scripts/10_run_prokka.sh
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
ASSEMBLY="${REPO_ROOT}/${ASSEMBLY}"
OUTDIR="${REPO_ROOT}/${OUTDIR}"

# ── Check tools ──
if ! command -v prokka &>/dev/null; then
    echo "ERROR: prokka not found. Activate the ngs-annotation environment."
    exit 1
fi

# ── Check input ──
if [ ! -f "$ASSEMBLY" ]; then
    echo "ERROR: Assembly not found at $ASSEMBLY"
    exit 1
fi

PROKKA_OUT="${OUTDIR}/prokka"
mkdir -p "$OUTDIR"

echo "======================================"
echo " Step 1: Structural annotation (Prokka)"
echo "======================================"
echo ""
echo "  Assembly:  $ASSEMBLY"
echo "  Output:    $PROKKA_OUT"
echo "  CPUs:      ${CPUS}"
echo "  Kingdom:   ${PROKKA_KINGDOM}"
echo ""

START=$(date +%s)

prokka \
    --outdir "$PROKKA_OUT" \
    --prefix "${PROKKA_PREFIX}" \
    --kingdom "${PROKKA_KINGDOM}" \
    --cpus "${CPUS}" \
    --force \
    --compliant \
    "$ASSEMBLY"

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> Prokka finished in ${ELAPSED}s"
echo ""
echo "Key output files:"
echo "  GFF:  ${PROKKA_OUT}/${PROKKA_PREFIX}.gff"
echo "  FAA:  ${PROKKA_OUT}/${PROKKA_PREFIX}.faa  (protein sequences)"
echo "  FFN:  ${PROKKA_OUT}/${PROKKA_PREFIX}.ffn  (nucleotide CDS)"
echo "  TSV:  ${PROKKA_OUT}/${PROKKA_PREFIX}.tsv  (feature table)"
echo "  TXT:  ${PROKKA_OUT}/${PROKKA_PREFIX}.txt  (summary stats)"
echo ""

if command -v seqkit &>/dev/null; then
    echo "Quick stats on predicted proteins:"
    seqkit stats "${PROKKA_OUT}/${PROKKA_PREFIX}.faa"
fi

echo ""
echo "Next step: bash scripts/20_run_diamond.sh"
