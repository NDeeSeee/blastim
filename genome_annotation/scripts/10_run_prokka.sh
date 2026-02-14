#!/usr/bin/env bash
set -euo pipefail

# ── 10_run_prokka.sh ────────────────────────────────────────────────
# Structural annotation with Prokka.
# Predicts CDS, rRNA, tRNA and produces GFF, FAA, FFN, TSV, etc.
##gff-version 3
##sequence-region gnl|Prokka|OHEIKNII_1 1 4215606
gnl|Prokka|OHEIKNII_1   prokka  gene    410     1750    .       +       .       ID=OHEIKNII_00001_gene;Name=dnaA_1;gene=dnaA_1;locus_tag=OHEIKNII_00001
gnl|Prokka|OHEIKNII_1   Prodigal:002006 CDS     410     1750    .       +       0       ID=OHEIKNII_00001;Parent=OHEIKNII_00001_gene;Name=dnaA_1;db_xref=COG:COG0593;gene=dnaA_1;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P05648;locus_tag=OHEIKNII_00001;product=Chromosomal replication initiator protein DnaA;protein_id=gnl|Prokka|OHEIKNII_00001
gnl|Prokka|OHEIKNII_1   prokka  gene    1939    3075    .       +       .       ID=OHEIKNII_00002_gene;Name=dnaN;gene=dnaN;locus_tag=OHEIKNII_00002
gnl|Prokka|OHEIKNII_1   Prodigal:002006 CDS     1939    3075    .       +       0       ID=OHEIKNII_00002;Parent=OHEIKNII_00002_gene;Name=dnaN;db_xref=COG:COG0592;gene=dnaN;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P05649;locus_tag=OHEIKNII_00002;product=Beta sliding clamp;protein_id=gnl|Prokka|OHEIKNII_00002
gnl|Prokka|OHEIKNII_1   prokka  gene    3206    3421    .       +       .       ID=OHEIKNII_00003_gene;locus_tag=OHEIKNII_00003
gnl|Prokka|OHEIKNII_1   Prodigal:002006 CDS     3206    3421    .       +       0       ID=OHEIKNII_00003;Parent=OHEIKNII_00003_gene;inference=ab initio prediction:Prodigal:002006;locus_tag=OHEIKNII_00003;product=hypothetical protein;protein_id=gnl|Prokka|OHEIKNII_00003
gnl|Prokka|OHEIKNII_1   prokka  gene    3437    4549    .       +       .       ID=OHEIKNII_00004_gene;Name=recF;gene=recF;locus_tag=OHEIKNII_00004
gnl|Prokka|OHEIKNII_1   Prodigal:002006 CDS     3437    4549    .       +       0       ID=OHEIKNII_00004;Parent=OHEIKNII_00004_gene;Name=recF;db_xref=COG:COG1195;gene=recF;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:Q8RDL3;locus_tag=OHEIKNII_00004;product=DNA replication and repair protein RecF;protein_id=gnl|Prokka|OHEIKNII_00004
gnl|Prokka|OHEIKNII_1   prokka  gene    4567    4812    .       +       .       ID=OHEIKNII_00005_gene;locus_tag=OHEIKNII_00005
gnl|Prokka|OHEIKNII_1   Prodigal:002006 CDS     4567    4812    .       +       0       ID=OHEIKNII_00005;Parent=OHEIKNII_00005_gene;inference=ab initio prediction:Prodigal:002006;locus_tag=OHEIKNII_00005;product=hypothetical protein;protein_id=gnl|Prokka|OHEIKNII_00005
gnl|Prokka|OHEIKNII_1   prokka  gene    4861    6783    .       +       .       ID=OHEIKNII_00006_gene;Name=gyrB;gene=gyrB;locus_tag=OHEIKNII_00006
gnl|Prokka|OHEIKNII_1   Prodigal:002006 CDS     4861    6783    .       +       0       ID=OHEIKNII_00006;Parent=OHEIKNII_00006_gene;eC_number=5.6.2.2;Name=gyrB;db_xref=COG:COG0187;gene=gyrB;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:Q839Z1;locus_tag=OHEIKNII_00006;product=DNA gyrase subunit B;protein_id=gnl|Prokka|OHEIKNII_00006
gnl|Prokka|OHEIKNII_1   prokka  gene    6994    9459    .       +       .       ID=OHEIKNII_00007_gene;Name=gyrA;gene=gyrA;locus_tag=OHEIKNII_00007
gnl|Prokka|OHEIKNII_1   Prodigal:002006 CDS     6994    9459    .       +       0       ID=OHEIKNII_00007;Parent=OHEIKNII_00007_gene;eC_number=5.6.2.2;Name=gyrA;db_xref=COG:COG0188;gene=gyrA;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P05653;locus_tag=OHEIKNII_00007;product=DNA gyrase subunit A;protein_id=gnl|Prokka|OHEIKNII_00007
gnl|Prokka|OHEIKNII_1   prokka  gene    9815    11361   .       +       .       ID=OHEIKNII_00008_gene;locus_tag=OHEIKNII_00008
gnl|Prokka|OHEIKNII_1   barrnap:0.9     rRNA    9815    11361   0       +       .       ID=OHEIKNII_00008;Parent=OHEIKNII_00008_gene;locus_tag=OHEIKNII_00008;product=16S ribosomal RNA
gnl|Prokka|OHEIKNII_1   prokka  gene    11464   11540   .       +       .       ID=OHEIKNII_00009_gene;locus_tag=OHEIKNII_00009
gnl|Prokka|OHEIKNII_1   Aragorn:001002  tRNA    11464   11540   .       +       .       ID=OHEIKNII_00009;Parent=OHEIKNII_00009_gene;inference=COORDINATES:profile:Aragorn:001002;locus_tag=OHEIKNII_00009;product=tRNA-Ile(gat)
gnl|Prokka|OHEIKNII_1   prokka  gene    11552   11627   .       +       .       ID=OHEIKNII_00010_gene;locus_tag=OHEIKNII_00010
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
