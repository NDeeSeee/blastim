#!/usr/bin/env bash
set -euo pipefail

# ── 01_get_data.sh ──────────────────────────────────────────────────
# Download a small bacterial genome assembly and a small protein set
# for the teaching DIAMOND database.
#
# Uses Bacillus subtilis 168 (~4.2 Mb) as the teaching genome.
# Uses a subset of Swiss-Prot bacterial proteins as the teaching db.
#
# Usage:
#   bash scripts/01_get_data.sh              # download both
#   bash scripts/01_get_data.sh --offline    # skip download, check existing
# ────────────────────────────────────────────────────────────────────

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_DIR="${REPO_ROOT}/data/input"
DB_DIR="${REPO_ROOT}/data/db"
OFFLINE=false

for arg in "$@"; do
    case "$arg" in
        --offline) OFFLINE=true ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

mkdir -p "$INPUT_DIR" "$DB_DIR"

# ── 1. Teaching genome: B. subtilis 168 from NCBI ──
ASSEMBLY_FILE="${INPUT_DIR}/assembly.fasta"
# NCBI RefSeq assembly for B. subtilis subsp. subtilis str. 168
GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz"

echo "==> Step 1: Teaching genome (Bacillus subtilis 168, ~4.2 Mb)"
if [ -f "$ASSEMBLY_FILE" ]; then
    echo "    Assembly already exists: $ASSEMBLY_FILE"
    echo "    Size: $(du -h "$ASSEMBLY_FILE" | cut -f1)"
else
    if $OFFLINE; then
        echo "ERROR: $ASSEMBLY_FILE not found and --offline mode is set."
        echo "       Place the assembly FASTA in data/input/assembly.fasta"
        exit 1
    fi
    echo "    Downloading from NCBI..."
    wget -q --show-progress -O "${ASSEMBLY_FILE}.gz" "$GENOME_URL"
    echo "    Decompressing..."
    gunzip -f "${ASSEMBLY_FILE}.gz"
    echo "    Done: $(du -h "$ASSEMBLY_FILE" | cut -f1)"
fi

# Quick validation
if command -v seqkit &>/dev/null; then
    echo "    Sequences: $(seqkit stats -T "$ASSEMBLY_FILE" | tail -1 | cut -f4)"
fi

# ── 2. Teaching protein database (Swiss-Prot bacteria subset) ──
PROT_FILE="${DB_DIR}/teaching_proteins.fasta"
# Swiss-Prot reviewed bacterial proteins (small, ~100k sequences)
SPROT_URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

echo ""
echo "==> Step 2: Teaching protein database"
if [ -f "$PROT_FILE" ]; then
    echo "    Protein file already exists: $PROT_FILE"
    echo "    Size: $(du -h "$PROT_FILE" | cut -f1)"
else
    if $OFFLINE; then
        echo "ERROR: $PROT_FILE not found and --offline mode is set."
        echo "       Place teaching_proteins.fasta in data/db/"
        exit 1
    fi

    SPROT_GZ="${DB_DIR}/uniprot_sprot.fasta.gz"
    if [ ! -f "$SPROT_GZ" ]; then
        echo "    Downloading Swiss-Prot (this is ~90 MB, one-time download)..."
        wget -q --show-progress -O "$SPROT_GZ" "$SPROT_URL"
    fi

    echo "    Extracting first 100,000 proteins for teaching database..."
    # Take the first 100k proteins from Swiss-Prot for a manageable teaching db
    zcat "$SPROT_GZ" \
        | awk 'BEGIN{n=0} /^>/{n++; if(n>100000) exit} {print}' \
        > "$PROT_FILE"

    echo "    Done: $(du -h "$PROT_FILE" | cut -f1)"
fi

if command -v seqkit &>/dev/null && [ -f "$PROT_FILE" ]; then
    echo "    Proteins: $(seqkit stats -T "$PROT_FILE" | tail -1 | cut -f4)"
fi

echo ""
echo "==> Data preparation complete."
echo "    Assembly:  $ASSEMBLY_FILE"
echo "    Proteins:  $PROT_FILE"
echo ""
echo "Next step: bash scripts/02_make_diamond_db.sh"
