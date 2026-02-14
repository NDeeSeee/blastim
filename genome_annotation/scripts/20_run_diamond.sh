#!/usr/bin/env bash
set -euo pipefail

# ── 20_run_diamond.sh ───────────────────────────────────────────────
# Functional annotation by homology using DIAMOND blastp.
# Searches predicted proteins against the teaching protein database.
#
# Usage:
#   bash scripts/20_run_diamond.sh
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

# ── Resolve paths ──
OUTDIR="${REPO_ROOT}/${OUTDIR}"
DIAMOND_DB="${REPO_ROOT}/${DIAMOND_DB}"
PROKKA_OUT="${OUTDIR}/prokka"
FAA="${PROKKA_OUT}/${PROKKA_PREFIX}.faa"

# ── Check tools ──
if ! command -v diamond &>/dev/null; then
    echo "ERROR: diamond not found. Activate the ngs-annotation environment."
    exit 1
fi

# ── Check inputs ──
if [ ! -f "$FAA" ]; then
    echo "ERROR: Prokka protein file not found at $FAA"
    echo "       Run scripts/10_run_prokka.sh first."
    exit 1
fi

if [ ! -f "${DIAMOND_DB}.dmnd" ]; then
    echo "ERROR: DIAMOND database not found at ${DIAMOND_DB}.dmnd"
    echo "       Run scripts/02_make_diamond_db.sh first."
    exit 1
fi

DIAMOND_OUT="${OUTDIR}/diamond"
mkdir -p "$DIAMOND_OUT"
RESULT="${DIAMOND_OUT}/blastp_results.tsv"

echo "======================================"
echo " Step 2: Homology search (DIAMOND)"
echo "======================================"
echo ""
echo "  Query:     $FAA"
echo "  Database:  ${DIAMOND_DB}.dmnd"
echo "  Output:    $RESULT"
echo "  CPUs:      ${CPUS}"
echo ""

START=$(date +%s)

diamond blastp \
    --db "$DIAMOND_DB" \
    --query "$FAA" \
    --out "$RESULT" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --evalue 1e-5 \
    --max-target-seqs 5 \
    --threads "${CPUS}" \
    --very-sensitive

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> DIAMOND finished in ${ELAPSED}s"
echo ""

# ── Quick summary ──
TOTAL_QUERIES=$(grep -c "^>" "$FAA" || echo 0)
QUERIES_WITH_HITS=$(cut -f1 "$RESULT" | sort -u | wc -l | tr -d ' ')
echo "  Total query proteins:  $TOTAL_QUERIES"
echo "  Proteins with hits:    $QUERIES_WITH_HITS"
echo "  Hit rate:              $(awk "BEGIN{printf \"%.1f\", ($QUERIES_WITH_HITS/$TOTAL_QUERIES)*100}")%"
echo ""
echo "  Output columns: qseqid sseqid pident length mismatch gapopen"
echo "                   qstart qend sstart send evalue bitscore stitle"
echo ""
echo "  First 5 hits:"
head -5 "$RESULT" | column -t -s $'\t'

echo ""
echo "Next step: bash scripts/30_run_hmmscan.sh"
