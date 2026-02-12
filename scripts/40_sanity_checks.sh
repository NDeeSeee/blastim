#!/usr/bin/env bash
set -euo pipefail

# ── 40_sanity_checks.sh ────────────────────────────────────────────
# Sanity checks on annotation outputs.
# Counts CDS, checks protein length distribution, summarises hits.
#
# Usage:
#   bash scripts/40_sanity_checks.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# ── Load config ──
CONFIG="${REPO_ROOT}/workflow/config.sh"
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config not found at $CONFIG"
    exit 1
fi
source "$CONFIG"

# ── Resolve paths ──
OUTDIR="${REPO_ROOT}/${OUTDIR}"
PROKKA_OUT="${OUTDIR}/prokka"
GFF="${PROKKA_OUT}/${PROKKA_PREFIX}.gff"
FAA="${PROKKA_OUT}/${PROKKA_PREFIX}.faa"
TSV="${PROKKA_OUT}/${PROKKA_PREFIX}.tsv"
TXT="${PROKKA_OUT}/${PROKKA_PREFIX}.txt"
DIAMOND_RESULT="${OUTDIR}/diamond/blastp_results.tsv"
DOMTBLOUT="${OUTDIR}/hmmer/hmmscan_domtblout.txt"

echo "======================================"
echo " Step 4: Sanity Checks"
echo "======================================"
echo ""

# ── 1. Prokka summary ──
echo "── 1. Prokka summary ──"
if [ -f "$TXT" ]; then
    cat "$TXT"
else
    echo "  WARNING: $TXT not found"
fi
echo ""

# ── 2. CDS count from GFF ──
echo "── 2. Feature counts from GFF ──"
if [ -f "$GFF" ]; then
    echo "  Feature type counts:"
    grep -v "^#" "$GFF" | awk -F'\t' '$3!=""{print $3}' | sort | uniq -c | sort -rn | head -15
else
    echo "  WARNING: $GFF not found"
fi
echo ""

# ── 3. Protein length distribution ──
echo "── 3. Protein length distribution ──"
if [ -f "$FAA" ] && command -v seqkit &>/dev/null; then
    echo "  Overall stats:"
    seqkit stats "$FAA"
    echo ""
    echo "  Length distribution (amino acids):"
    seqkit fx2tab -l "$FAA" | awk '{print $NF}' | sort -n | awk '
    BEGIN { n=0 }
    {
        vals[n++] = $1
        sum += $1
    }
    END {
        printf "    Min:    %d\n", vals[0]
        printf "    Q1:     %d\n", vals[int(n*0.25)]
        printf "    Median: %d\n", vals[int(n*0.5)]
        printf "    Q3:     %d\n", vals[int(n*0.75)]
        printf "    Max:    %d\n", vals[n-1]
        printf "    Mean:   %.1f\n", sum/n
        printf "    N:      %d\n", n
    }'

    echo ""
    echo "  Short proteins (<100 aa):"
    NSHORT=$(seqkit fx2tab -l "$FAA" | awk '$NF < 100' | wc -l | tr -d ' ')
    NTOTAL=$(grep -c "^>" "$FAA" || echo 0)
    echo "    ${NSHORT} / ${NTOTAL} proteins ($(awk "BEGIN{printf \"%.1f\", ($NSHORT/$NTOTAL)*100}")%)"
elif [ -f "$FAA" ]; then
    NTOTAL=$(grep -c "^>" "$FAA" || echo 0)
    echo "  Total proteins: $NTOTAL"
    echo "  (Install seqkit for detailed length stats)"
fi
echo ""

# ── 4. DIAMOND hit summary ──
echo "── 4. DIAMOND homology summary ──"
if [ -f "$DIAMOND_RESULT" ]; then
    TOTAL_HITS=$(wc -l < "$DIAMOND_RESULT" | tr -d ' ')
    UNIQUE_QUERIES=$(cut -f1 "$DIAMOND_RESULT" | sort -u | wc -l | tr -d ' ')
    NTOTAL=$(grep -c "^>" "$FAA" || echo 0)

    echo "  Total hits:             $TOTAL_HITS"
    echo "  Proteins with hits:     $UNIQUE_QUERIES / $NTOTAL"
    echo "  Proteins without hits:  $(( NTOTAL - UNIQUE_QUERIES )) (hypothetical / novel)"
    echo ""

    echo "  Identity distribution of best hits:"
    # Take best hit per query (first occurrence)
    sort -k1,1 -k12,12rn "$DIAMOND_RESULT" | sort -k1,1 -u | awk -F'\t' '{print $3}' | awk '
    BEGIN { b[0]="0-30%"; b[1]="30-50%"; b[2]="50-70%"; b[3]="70-90%"; b[4]="90-100%" }
    {
        if ($1 < 30) c[0]++
        else if ($1 < 50) c[1]++
        else if ($1 < 70) c[2]++
        else if ($1 < 90) c[3]++
        else c[4]++
    }
    END {
        for (i=0; i<=4; i++) printf "    %-10s %d\n", b[i], c[i]+0
    }'
else
    echo "  WARNING: DIAMOND results not found"
fi
echo ""

# ── 5. hmmscan domain summary ──
echo "── 5. hmmscan domain summary ──"
if [ -f "$DOMTBLOUT" ]; then
    NHITS=$(grep -cv "^#" "$DOMTBLOUT" || echo 0)
    NPROTEINS_HIT=$(grep -v "^#" "$DOMTBLOUT" | awk '{print $4}' | sort -u | wc -l | tr -d ' ')
    NDOMAINS=$(grep -v "^#" "$DOMTBLOUT" | awk '{print $1}' | sort -u | wc -l | tr -d ' ')

    echo "  Total domain hits:       $NHITS"
    echo "  Proteins with domains:   $NPROTEINS_HIT"
    echo "  Unique domains found:    $NDOMAINS"
    echo ""
    echo "  Top 15 most frequent domains:"
    grep -v "^#" "$DOMTBLOUT" | awk '{print $1, $2}' | sort | uniq -c | sort -rn | head -15 | \
        awk '{printf "    %4d  %-20s %s\n", $1, $2, $3}'
else
    echo "  WARNING: hmmscan results not found"
fi
echo ""

# ── 6. Evidence ladder summary ──
echo "── 6. Evidence ladder summary ──"
if [ -f "$FAA" ] && [ -f "$DIAMOND_RESULT" ] && [ -f "$DOMTBLOUT" ]; then
    NTOTAL=$(grep -c "^>" "$FAA")
    N_DIAMOND=$(cut -f1 "$DIAMOND_RESULT" | sort -u | wc -l | tr -d ' ')
    N_HMMER=$(grep -v "^#" "$DOMTBLOUT" | awk '{print $4}' | sort -u | wc -l | tr -d ' ')

    # Proteins with both homology AND domain evidence
    N_BOTH=$(comm -12 \
        <(cut -f1 "$DIAMOND_RESULT" | sort -u) \
        <(grep -v "^#" "$DOMTBLOUT" | awk '{print $4}' | sort -u) \
        | wc -l | tr -d ' ')

    # Proteins with no evidence at all
    N_NONE=$(comm -23 \
        <(grep "^>" "$FAA" | sed 's/^>//' | awk '{print $1}' | sort) \
        <(cat <(cut -f1 "$DIAMOND_RESULT") <(grep -v "^#" "$DOMTBLOUT" | awk '{print $4}') | sort -u) \
        | wc -l | tr -d ' ')

    echo "  Total predicted proteins:           $NTOTAL"
    echo "  With homology (DIAMOND) only:       $(( N_DIAMOND - N_BOTH ))"
    echo "  With domain (hmmscan) only:         $(( N_HMMER - N_BOTH ))"
    echo "  With BOTH homology + domain:        $N_BOTH  (strongest evidence)"
    echo "  With NO functional annotation:      $N_NONE  (hypothetical proteins)"
fi

echo ""
echo "======================================"
echo " Sanity checks complete."
echo "======================================"
echo ""
echo "Next step: bash scripts/50_pick_3_genes.sh"
