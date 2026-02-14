#!/usr/bin/env bash
set -euo pipefail

# ── 02_make_diamond_db.sh ───────────────────────────────────────────
# Build a DIAMOND database from the teaching protein FASTA.
#
# Usage:
#   bash scripts/02_make_diamond_db.sh
# ────────────────────────────────────────────────────────────────────

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${REPO_ROOT}/data/db"
PROT_FILE="${DB_DIR}/teaching_proteins.fasta"
DMND_FILE="${DB_DIR}/teaching_proteins.dmnd"

echo "==> Building DIAMOND database"

# ── Check prerequisites ──
if ! command -v diamond &>/dev/null; then
    echo "ERROR: diamond not found. Activate the ngs-annotation environment first."
    exit 1
fi

if [ ! -f "$PROT_FILE" ]; then
    echo "ERROR: Protein FASTA not found at $PROT_FILE"
    echo "       Run scripts/01_get_data.sh first."
    exit 1
fi

echo "    Input:  $PROT_FILE ($(du -h "$PROT_FILE" | cut -f1))"

# ── Build ──
if [ -f "$DMND_FILE" ]; then
    echo "    DIAMOND db already exists: $DMND_FILE"
    echo "    Delete it first if you want to rebuild."
else
    diamond makedb \
        --in "$PROT_FILE" \
        --db "${DMND_FILE%.dmnd}" \
        --threads "$(nproc 2>/dev/null || echo 4)"

    echo "    Output: $DMND_FILE ($(du -h "$DMND_FILE" | cut -f1))"
fi

# ── Validate ──
echo ""
echo "    Validating database..."
diamond dbinfo --db "${DMND_FILE%.dmnd}" 2>&1 | head -5
echo ""
echo "==> DIAMOND database ready."
echo "Next step: bash scripts/03_get_pfam.sh --subset"
