# Prokaryotic Genome Annotation on HPC

A hands-on 2-hour practical for NGS courses. Students learn to annotate a bacterial genome using standard bioinformatics tools on an HPC system.

## What you will learn

1. **Structural annotation** -- predict protein-coding genes (CDS), rRNA, tRNA using Prokka. Understand GFF/FAA/TSV output formats.
2. **Functional annotation by homology** -- search predicted proteins against a reference database with DIAMOND blastp.
3. **Functional annotation by domains** -- identify conserved protein domains with hmmscan (HMMER) against Pfam.
4. **The evidence ladder** -- combine ORF prediction, homology, and domain evidence to assess annotation confidence.

## Prerequisites

- SSH access to the HPC
- Basic Linux command-line skills (`cd`, `ls`, `less`, `head`)
- The `ngs-annotation` conda environment (your instructor will set this up)

## Quickstart

```bash
# 1. Clone the repository
git clone git@github.com:NDeeSeee/blastim.git
cd blastim

# 2. Create and activate the environment
bash env/micromamba_create.sh
micromamba activate ngs-annotation   # or: conda activate ngs-annotation

# 3. Copy and review the config
cp workflow/config.sh.example workflow/config.sh
# Edit workflow/config.sh if needed (e.g. change CPUs)

# 4. Verify setup
bash scripts/00_check_system.sh

# 5. Run the complete workflow
bash workflow/run_all.sh
```

Or step-by-step:

```bash
bash scripts/10_run_prokka.sh       # Structural annotation
bash scripts/20_run_diamond.sh      # Homology search
bash scripts/30_run_hmmscan.sh      # Domain search
bash scripts/40_sanity_checks.sh    # Quality checks
bash scripts/50_pick_3_genes.sh     # Discuss example genes
```

## Understanding the results

After the workflow completes, your output directory (`outputs/run_YYYYMMDD_HHMMSS/`) contains:

| Directory | Key files | What they are |
|-----------|-----------|---------------|
| `prokka/` | `ANNOT.gff` | Annotated genome in GFF3 format |
| | `ANNOT.faa` | Predicted protein sequences (FASTA) |
| | `ANNOT.ffn` | Nucleotide CDS sequences |
| | `ANNOT.tsv` | Tab-separated feature table |
| | `ANNOT.txt` | Summary statistics |
| `diamond/` | `blastp_results.tsv` | DIAMOND blastp hits (BLAST tabular format) |
| `hmmer/` | `hmmscan_domtblout.txt` | Pfam domain hits per protein |

### How to interpret

- **GFF file**: Each line is a genomic feature. Look at column 3 (type) and column 9 (attributes, including product name).
- **DIAMOND TSV**: Standard BLAST tabular. Key columns: query ID, subject ID, % identity, E-value, bit score, subject description.
- **domtblout**: Each line is a domain hit. Key columns: domain name, accession, query protein, E-value.

## Repository structure

```
blastim/
├── README.md              ← You are here
├── INSTRUCTOR.md          ← Instructor setup guide
├── LICENSE
├── Makefile               ← make env / make all / make clean
├── env/
│   ├── environment.yml    ← Conda environment spec
│   └── micromamba_create.sh
├── scripts/
│   ├── 00_check_system.sh ← Verify HPC environment
│   ├── 01_get_data.sh     ← Download genome + proteins
│   ├── 02_make_diamond_db.sh
│   ├── 03_get_pfam.sh     ← Download & subset Pfam
│   ├── 10_run_prokka.sh   ← Structural annotation
│   ├── 20_run_diamond.sh  ← Homology search
│   ├── 30_run_hmmscan.sh  ← Domain search
│   ├── 40_sanity_checks.sh
│   └── 50_pick_3_genes.sh ← Discussion examples
├── workflow/
│   ├── config.sh.example  ← Copy to config.sh
│   └── run_all.sh         ← One-command workflow
├── data/
│   ├── input/             ← Assembly FASTA
│   └── db/                ← DIAMOND + Pfam databases
├── outputs/               ← Created at runtime
└── docs/
    ├── cheatsheet.md      ← All commands copy/paste
    ├── expected_outputs.md
    └── troubleshooting.md
```

## Organism

This course uses **Bacillus subtilis subsp. subtilis str. 168** (~4.2 Mb genome, ~4,200 predicted proteins). It is a well-studied model organism, so most predictions will have homology and domain evidence -- ideal for teaching.

## License

MIT. See [LICENSE](LICENSE).

Databases (Swiss-Prot, Pfam) are downloaded from their official sources and subject to their own licenses. They are not included in this repository.
