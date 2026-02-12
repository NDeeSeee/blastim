#!/usr/bin/env bash
set -euo pipefail

# ── 03_get_pfam.sh ──────────────────────────────────────────────────
# Download Pfam-A HMMs and prepare for hmmscan.
#
# Modes:
#   --subset  (default) Download full Pfam-A, then extract ~200 common
#             bacterial domain HMMs for fast in-class runs (~2-5 min).
#   --full    Download and hmmpress the full Pfam-A (~20k profiles).
#             hmmscan will be slower (~15-30 min on 8 CPUs for a 4 Mb genome).
#
# Usage:
#   bash scripts/03_get_pfam.sh              # subset mode (default)
#   bash scripts/03_get_pfam.sh --subset     # same as above
#   bash scripts/03_get_pfam.sh --full       # full Pfam-A
# ────────────────────────────────────────────────────────────────────

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${REPO_ROOT}/data/db"
MODE="subset"

for arg in "$@"; do
    case "$arg" in
        --subset) MODE="subset" ;;
        --full)   MODE="full" ;;
        *)        echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

mkdir -p "$DB_DIR"

# ── Check tools ──
for tool in hmmpress hmmfetch wget; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool not found. Activate the ngs-annotation environment."
        exit 1
    fi
done

# ── Download Pfam-A ──
PFAM_URL="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_FULL="${DB_DIR}/Pfam-A.hmm"
PFAM_GZ="${DB_DIR}/Pfam-A.hmm.gz"

echo "==> Pfam HMM database preparation (mode: ${MODE})"

if [ ! -f "$PFAM_FULL" ]; then
    if [ ! -f "$PFAM_GZ" ]; then
        echo "    Downloading Pfam-A.hmm.gz from EBI (~300 MB)..."
        wget -q --show-progress -O "$PFAM_GZ" "$PFAM_URL"
    fi
    echo "    Decompressing..."
    gunzip -kf "$PFAM_GZ"
    rm -f "$PFAM_GZ"
    echo "    Full Pfam-A: $(du -h "$PFAM_FULL" | cut -f1)"
else
    echo "    Full Pfam-A already exists: $(du -h "$PFAM_FULL" | cut -f1)"
fi

if [ "$MODE" = "full" ]; then
    # ── Full Pfam mode ──
    echo ""
    echo "    Pressing full Pfam-A database (this may take a few minutes)..."
    if [ ! -f "${PFAM_FULL}.h3m" ]; then
        hmmpress "$PFAM_FULL"
    else
        echo "    Already pressed."
    fi
    echo "==> Full Pfam-A ready for hmmscan."
    echo "    Set HMM_DB=\"data/db/Pfam-A.hmm\" in config.sh"
    exit 0
fi

# ── Subset mode ──
SUBSET_FILE="${DB_DIR}/pfam_subset.hmm"
echo ""
echo "    Creating subset of ~200 common bacterial domain HMMs..."

# List of common bacterial Pfam domains (well-studied, commonly found)
# These cover: metabolism, transport, regulation, cell wall, stress, etc.
DOMAIN_LIST=(
    PF00005 PF00009 PF00012 PF00013 PF00015 PF00023 PF00027 PF00028
    PF00043 PF00044 PF00069 PF00070 PF00072 PF00076 PF00077 PF00078
    PF00080 PF00082 PF00085 PF00091 PF00092 PF00096 PF00106 PF00107
    PF00108 PF00109 PF00111 PF00115 PF00117 PF00118 PF00119 PF00120
    PF00121 PF00122 PF00126 PF00132 PF00133 PF00137 PF00144 PF00146
    PF00149 PF00150 PF00152 PF00155 PF00156 PF00158 PF00162 PF00164
    PF00166 PF00168 PF00170 PF00171 PF00173 PF00174 PF00175 PF00176
    PF00177 PF00180 PF00181 PF00183 PF00185 PF00186 PF00196 PF00198
    PF00202 PF00204 PF00205 PF00206 PF00210 PF00211 PF00213 PF00216
    PF00224 PF00227 PF00230 PF00231 PF00232 PF00233 PF00237 PF00238
    PF00239 PF00244 PF00248 PF00252 PF00253 PF00254 PF00266 PF00268
    PF00270 PF00271 PF00275 PF00276 PF00281 PF00282 PF00283 PF00285
    PF00288 PF00289 PF00290 PF00291 PF00293 PF00294 PF00300 PF00306
    PF00308 PF00309 PF00310 PF00312 PF00313 PF00318 PF00319 PF00326
    PF00327 PF00328 PF00330 PF00333 PF00334 PF00338 PF00344 PF00347
    PF00348 PF00353 PF00355 PF00357 PF00358 PF00361 PF00364 PF00365
    PF00366 PF00370 PF00378 PF00380 PF00381 PF00382 PF00383 PF00384
    PF00389 PF00390 PF00391 PF00392 PF00393 PF00394 PF00396 PF00398
    PF00400 PF00403 PF00406 PF00410 PF00411 PF00440 PF00441 PF00462
    PF00464 PF00465 PF00467 PF00480 PF00481 PF00484 PF00486 PF00488
    PF00496 PF00497 PF00501 PF00503 PF00507 PF00512 PF00515 PF00516
    PF00517 PF00519 PF00520 PF00525 PF00528 PF00534 PF00535 PF00537
    PF00542 PF00543 PF00550 PF00551 PF00561 PF00571 PF00572 PF00575
    PF00578 PF00581 PF00583 PF00585 PF00588 PF00589 PF00590 PF00593
    PF00595 PF00596 PF00623 PF00624 PF00627 PF00628 PF00633 PF00636
)

# Extract subset using awk (handles versioned accessions like PF00005.30)
echo "    Extracting ${#DOMAIN_LIST[@]} domains from full Pfam-A..."

# Build a lookup set from the domain list
ACC_FILE="${DB_DIR}/.pfam_subset_accessions.txt"
printf '%s\n' "${DOMAIN_LIST[@]}" > "$ACC_FILE"

# Extract matching HMM records directly from the flat file.
# Pfam accessions are versioned (PF00005.30); we match the base accession.
awk '
    BEGIN {
        while ((getline acc < ARGV[1]) > 0) wanted[acc] = 1
        delete ARGV[1]
    }
    /^HMMER/ { rec = $0 "\n"; capturing = 0; next }
    /^ACC / {
        base = $2
        sub(/\.[0-9]+$/, "", base)
        capturing = (base in wanted)
        rec = rec $0 "\n"
        next
    }
    /^\/\// {
        if (capturing) print rec $0
        rec = ""
        capturing = 0
        next
    }
    { rec = rec $0 "\n" }
' "$ACC_FILE" "$PFAM_FULL" > "$SUBSET_FILE"

rm -f "$ACC_FILE"

# Count how many we actually got
NHITS=$(grep -c "^ACC " "$SUBSET_FILE" 2>/dev/null || true)
NHITS="${NHITS:-0}"
echo "    Extracted $NHITS / ${#DOMAIN_LIST[@]} HMMs"

if [ "$NHITS" -eq 0 ]; then
    echo "WARNING: No HMMs extracted. Check Pfam-A format."
    echo "         Debug: head -20 $PFAM_FULL"
    exit 1
fi

# Press the subset
echo "    Pressing subset database..."
hmmpress -f "$SUBSET_FILE"

echo ""
echo "==> Subset Pfam database ready: $SUBSET_FILE ($NHITS profiles)"
echo "    Set HMM_DB=\"data/db/pfam_subset.hmm\" in config.sh"
echo ""
echo "Next step: Copy config and run the workflow!"
