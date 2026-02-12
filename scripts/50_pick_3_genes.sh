#!/usr/bin/env bash
set -euo pipefail

# ── 50_pick_3_genes.sh ─────────────────────────────────────────────
# Pick 3 example genes/proteins with their annotations and hits
# for class discussion.
#
# Picks:
#   1. A well-annotated protein (has both DIAMOND and hmmscan hits)
#   2. A protein with only homology (DIAMOND hit, no domain)
#   3. A hypothetical protein (no hits at all)
#
# Usage:
#   bash scripts/50_pick_3_genes.sh
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
FAA="${PROKKA_OUT}/${PROKKA_PREFIX}.faa"
GFF="${PROKKA_OUT}/${PROKKA_PREFIX}.gff"
DIAMOND_RESULT="${OUTDIR}/diamond/blastp_results.tsv"
DOMTBLOUT="${OUTDIR}/hmmer/hmmscan_domtblout.txt"

echo "======================================"
echo " Step 5: Pick 3 genes for discussion"
echo "======================================"
echo ""

# ── Check files exist ──
for f in "$FAA" "$GFF" "$DIAMOND_RESULT" "$DOMTBLOUT"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        echo "       Run all previous steps first."
        exit 1
    fi
done

# Build gene lists
TMPDIR_WORK=$(mktemp -d)
trap "rm -rf $TMPDIR_WORK" EXIT

# All protein IDs
grep "^>" "$FAA" | sed 's/^>//' | awk '{print $1}' | sort > "$TMPDIR_WORK/all.txt"

# Proteins with DIAMOND hits
cut -f1 "$DIAMOND_RESULT" | sort -u > "$TMPDIR_WORK/diamond.txt"

# Proteins with hmmscan hits
grep -v "^#" "$DOMTBLOUT" | awk '{print $4}' | sort -u > "$TMPDIR_WORK/hmmer.txt"

# Category 1: Both DIAMOND + hmmscan
comm -12 "$TMPDIR_WORK/diamond.txt" "$TMPDIR_WORK/hmmer.txt" > "$TMPDIR_WORK/both.txt"

# Category 2: DIAMOND only (no domain)
comm -23 "$TMPDIR_WORK/diamond.txt" "$TMPDIR_WORK/hmmer.txt" > "$TMPDIR_WORK/diamond_only.txt"

# Category 3: No hits at all
cat "$TMPDIR_WORK/diamond.txt" "$TMPDIR_WORK/hmmer.txt" | sort -u > "$TMPDIR_WORK/any_hit.txt"
comm -23 "$TMPDIR_WORK/all.txt" "$TMPDIR_WORK/any_hit.txt" > "$TMPDIR_WORK/no_hits.txt"

# ── Helper function to show a gene ──
show_gene() {
    local gene_id="$1"
    local category="$2"

    echo "────────────────────────────────────────"
    echo "GENE: $gene_id"
    echo "Category: $category"
    echo "────────────────────────────────────────"

    # GFF entry
    echo ""
    echo "  GFF annotation:"
    grep -w "$gene_id" "$GFF" | grep -v "^#" | head -3 | while IFS= read -r line; do
        echo "    $line"
    done

    # Protein sequence (first 3 lines)
    echo ""
    echo "  Protein sequence (first 80 aa):"
    if command -v seqkit &>/dev/null; then
        seqkit grep -p "$gene_id" "$FAA" | seqkit subseq -r 1:80 | tail -n +2 | while IFS= read -r line; do
            echo "    $line"
        done
    else
        awk -v id="$gene_id" '
            /^>/ { if (index($0, id)) {p=1; print "    "$0; next} else p=0 }
            p { printf "    %.80s\n", $0; p=0 }
        ' "$FAA"
    fi

    # DIAMOND hits
    echo ""
    echo "  DIAMOND hits (top 3):"
    DHITS=$(grep -w "^${gene_id}" "$DIAMOND_RESULT" | head -3)
    if [ -n "$DHITS" ]; then
        echo "$DHITS" | while IFS=$'\t' read -r qid sid pident alen mm go qs qe ss se eval bits desc; do
            echo "    Hit: $sid"
            echo "         Identity: ${pident}%  E-value: $eval  Score: $bits"
            echo "         Description: $desc"
        done
    else
        echo "    (no DIAMOND hits)"
    fi

    # hmmscan domains
    echo ""
    echo "  Pfam domains:"
    HHITS=$(grep -v "^#" "$DOMTBLOUT" | awk -v id="$gene_id" '$4==id' | head -5)
    if [ -n "$HHITS" ]; then
        echo "$HHITS" | while read -r line; do
            domain=$(echo "$line" | awk '{print $1}')
            acc=$(echo "$line" | awk '{print $2}')
            evalue=$(echo "$line" | awk '{print $13}')
            desc=$(echo "$line" | awk '{for(i=23;i<=NF;i++) printf "%s ", $i; print ""}')
            echo "    Domain: $domain ($acc)"
            echo "            E-value: $evalue"
            echo "            Description: $desc"
        done
    else
        echo "    (no Pfam domains detected)"
    fi
    echo ""
}

# ── Pick and display genes ──

echo "Selecting representative genes for discussion..."
echo ""

# Gene 1: Well-annotated (both hits)
if [ -s "$TMPDIR_WORK/both.txt" ]; then
    GENE1=$(shuf -n1 "$TMPDIR_WORK/both.txt" 2>/dev/null || head -1 "$TMPDIR_WORK/both.txt")
    show_gene "$GENE1" "WELL-ANNOTATED (homology + domain evidence)"
else
    echo "  No proteins found with both DIAMOND and hmmscan hits."
fi

# Gene 2: Homology only
if [ -s "$TMPDIR_WORK/diamond_only.txt" ]; then
    GENE2=$(shuf -n1 "$TMPDIR_WORK/diamond_only.txt" 2>/dev/null || head -1 "$TMPDIR_WORK/diamond_only.txt")
    show_gene "$GENE2" "HOMOLOGY ONLY (DIAMOND hit, no Pfam domain)"
else
    echo "  No proteins found with DIAMOND-only hits."
fi

# Gene 3: Hypothetical
if [ -s "$TMPDIR_WORK/no_hits.txt" ]; then
    GENE3=$(shuf -n1 "$TMPDIR_WORK/no_hits.txt" 2>/dev/null || head -1 "$TMPDIR_WORK/no_hits.txt")
    show_gene "$GENE3" "HYPOTHETICAL (no functional annotation)"
else
    echo "  No hypothetical proteins found (all have hits)."
fi

echo "======================================"
echo " Discussion prompts for class:"
echo "======================================"
echo ""
echo "  1. For the well-annotated gene: What can you learn from"
echo "     combining homology AND domain evidence? How confident"
echo "     are you in this annotation?"
echo ""
echo "  2. For the homology-only gene: Why might it lack a Pfam"
echo "     domain? Is the homology hit enough to annotate it?"
echo ""
echo "  3. For the hypothetical protein: What approaches could you"
echo "     use to learn more about this protein? (structure prediction,"
echo "     synteny, expression data, etc.)"
echo ""
echo "======================================"
