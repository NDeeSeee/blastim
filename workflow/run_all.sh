#!/usr/bin/env bash
set -euo pipefail

# ── run_all.sh ──────────────────────────────────────────────────────
# Run the complete annotation workflow (steps 10→50).
#
# Usage:
#   bash workflow/run_all.sh
#   bash workflow/run_all.sh --max-proteins 200   # Plan B for hmmscan
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

# Set a timestamped OUTDIR and export it so all child scripts share it.
export OUTDIR="outputs/run_$(date +%Y%m%d_%H%M%S)"

# Forward all arguments to hmmscan step
HMMSCAN_ARGS=("$@")

echo "╔════════════════════════════════════════════════╗"
echo "║  Prokaryotic Genome Annotation Workflow        ║"
echo "╚════════════════════════════════════════════════╝"
echo ""
echo "Config:  $CONFIG"
echo "Output:  $OUTDIR"
echo "Start:   $(date)"
echo ""

TOTAL_START=$(date +%s)

# ── Step 1: Prokka ──
echo "━━━ Running Step 10: Prokka ━━━"
bash "${REPO_ROOT}/scripts/10_run_prokka.sh"
echo ""

# ── Step 2: DIAMOND ──
echo "━━━ Running Step 20: DIAMOND ━━━"
bash "${REPO_ROOT}/scripts/20_run_diamond.sh"
echo ""

# ── Step 3: hmmscan ──
echo "━━━ Running Step 30: hmmscan ━━━"
bash "${REPO_ROOT}/scripts/30_run_hmmscan.sh" "${HMMSCAN_ARGS[@]+"${HMMSCAN_ARGS[@]}"}"
echo ""

# ── Step 4: Sanity checks ──
echo "━━━ Running Step 40: Sanity checks ━━━"
bash "${REPO_ROOT}/scripts/40_sanity_checks.sh"
echo ""

# ── Step 5: Pick 3 genes ──
echo "━━━ Running Step 50: Pick 3 genes ━━━"
bash "${REPO_ROOT}/scripts/50_pick_3_genes.sh"
echo ""

TOTAL_END=$(date +%s)
TOTAL_ELAPSED=$(( TOTAL_END - TOTAL_START ))
MINUTES=$(( TOTAL_ELAPSED / 60 ))
SECONDS_R=$(( TOTAL_ELAPSED % 60 ))

source "$CONFIG"

echo "╔════════════════════════════════════════════════╗"
echo "║  Workflow complete!                            ║"
echo "║  Total time: ${MINUTES}m ${SECONDS_R}s"
echo "╚════════════════════════════════════════════════╝"
echo ""
echo "Output directory: ${REPO_ROOT}/${OUTDIR}"
echo ""
echo "Key files to explore:"
echo "  ${OUTDIR}/prokka/${PROKKA_PREFIX}.gff         - Annotated genome (GFF3)"
echo "  ${OUTDIR}/prokka/${PROKKA_PREFIX}.faa         - Predicted proteins (FASTA)"
echo "  ${OUTDIR}/diamond/blastp_results.tsv   - Homology search results"
echo "  ${OUTDIR}/hmmer/hmmscan_domtblout.txt  - Domain search results"
