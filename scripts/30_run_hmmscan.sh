#!/usr/bin/env bash
set -euo pipefail

# ── 30_run_hmmscan.sh ───────────────────────────────────────────────
# Functional annotation by conserved domains using hmmscan (HMMER).
# Searches predicted proteins against Pfam (subset or full).
#
# Plan B options for slow runs:
#   --subset-only     Use only the subset HMM database (default)
#   --max-proteins N  Only scan the first N proteins (e.g. 200)
#
# Usage:
#   bash scripts/30_run_hmmscan.sh
#   bash scripts/30_run_hmmscan.sh --max-proteins 200
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

# ── Parse command-line overrides ──
CLI_MAX_PROTEINS=""
for arg in "$@"; do
    case "$arg" in
        --subset-only) HMM_DB="${REPO_ROOT}/data/db/pfam_subset.hmm" ;;
        --max-proteins)
            # next arg will be the number; handled below
            ;;
        [0-9]*)
            CLI_MAX_PROTEINS="$arg"
            ;;
        *)
            if [[ "$arg" == --max-proteins=* ]]; then
                CLI_MAX_PROTEINS="${arg#*=}"
            fi
            ;;
    esac
done

# Apply CLI override if given
if [ -n "$CLI_MAX_PROTEINS" ]; then
    MAX_PROTEINS="$CLI_MAX_PROTEINS"
fi

# ── Resolve paths ──
OUTDIR="${REPO_ROOT}/${OUTDIR}"
HMM_DB="${REPO_ROOT}/${HMM_DB}"
PROKKA_OUT="${OUTDIR}/prokka"
FAA="${PROKKA_OUT}/${PROKKA_PREFIX}.faa"

# ── Check tools ──
if ! command -v hmmscan &>/dev/null; then
    echo "ERROR: hmmscan not found. Activate the ngs-annotation environment."
    exit 1
fi

# ── Check inputs ──
if [ ! -f "$FAA" ]; then
    echo "ERROR: Prokka protein file not found at $FAA"
    echo "       Run scripts/10_run_prokka.sh first."
    exit 1
fi

if [ ! -f "$HMM_DB" ]; then
    echo "ERROR: HMM database not found at $HMM_DB"
    echo "       Run scripts/03_get_pfam.sh first."
    exit 1
fi

if [ ! -f "${HMM_DB}.h3m" ]; then
    echo "ERROR: HMM database not pressed (missing .h3m)."
    echo "       Run: hmmpress $HMM_DB"
    exit 1
fi

HMMER_OUT="${OUTDIR}/hmmer"
mkdir -p "$HMMER_OUT"
DOMTBLOUT="${HMMER_OUT}/hmmscan_domtblout.txt"
HMMSCAN_LOG="${HMMER_OUT}/hmmscan.log"

# ── Prepare query ──
QUERY_FAA="$FAA"
if [ "${MAX_PROTEINS:-0}" -gt 0 ]; then
    echo "  Plan B: limiting to first ${MAX_PROTEINS} proteins"
    QUERY_FAA="${HMMER_OUT}/query_subset.faa"
    if command -v seqkit &>/dev/null; then
        seqkit head -n "${MAX_PROTEINS}" "$FAA" > "$QUERY_FAA"
    else
        # Fallback: awk-based head for FASTA
        awk -v n="${MAX_PROTEINS}" '/^>/{c++; if(c>n) exit} {print}' "$FAA" > "$QUERY_FAA"
    fi
fi

NQUERY=$(grep -c "^>" "$QUERY_FAA" || echo 0)

echo "======================================"
echo " Step 3: Domain annotation (hmmscan)"
echo "======================================"
echo ""
echo "  Query:     $QUERY_FAA ($NQUERY proteins)"
echo "  Database:  $HMM_DB"
echo "  Output:    $DOMTBLOUT"
echo "  CPUs:      ${CPUS}"
if [ "${MAX_PROTEINS:-0}" -gt 0 ]; then
    echo "  Mode:      subset ($MAX_PROTEINS proteins)"
else
    echo "  Mode:      all proteins"
fi
echo ""

START=$(date +%s)

hmmscan \
    --cpu "${CPUS}" \
    --domtblout "$DOMTBLOUT" \
    --noali \
    -E 1e-5 \
    --domE 1e-3 \
    -o "$HMMSCAN_LOG" \
    "$HMM_DB" \
    "$QUERY_FAA"

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> hmmscan finished in ${ELAPSED}s"
echo ""

# ── Quick summary ──
NHITS=$(grep -cv "^#" "$DOMTBLOUT" || echo 0)
NPROTEINS_HIT=$(grep -v "^#" "$DOMTBLOUT" | awk '{print $4}' | sort -u | wc -l | tr -d ' ')
NDOMAINS=$(grep -v "^#" "$DOMTBLOUT" | awk '{print $1}' | sort -u | wc -l | tr -d ' ')

echo "  Total domain hits:       $NHITS"
echo "  Proteins with domains:   $NPROTEINS_HIT / $NQUERY"
echo "  Unique domains found:    $NDOMAINS"
echo ""
echo "  Top 10 most frequent domains:"
grep -v "^#" "$DOMTBLOUT" | awk '{print $1}' | sort | uniq -c | sort -rn | head -10
echo ""
echo "Next step: bash scripts/40_sanity_checks.sh"
