# Command Cheatsheet

All commands to copy/paste during the class. Run from the repo root.

## Setup

```bash
# Activate environment
micromamba activate ngs-annotation
# or: conda activate ngs-annotation

# Copy config
cp workflow/config.sh.example workflow/config.sh

# Check system
bash scripts/00_check_system.sh
```

## Run the workflow

```bash
# All at once
bash workflow/run_all.sh

# Or step by step
bash scripts/10_run_prokka.sh
bash scripts/20_run_diamond.sh
bash scripts/30_run_hmmscan.sh
bash scripts/40_sanity_checks.sh
bash scripts/50_pick_3_genes.sh
```

## Explore Prokka outputs

```bash
# Find your output directory
OUTDIR=$(ls -td outputs/run_* | head -1)

# View summary
cat $OUTDIR/prokka/ANNOT.txt

# Count features by type
grep -v "^#" $OUTDIR/prokka/ANNOT.gff | awk -F'\t' '{print $3}' | sort | uniq -c | sort -rn

# Look at the first 5 CDS entries in GFF
grep "CDS" $OUTDIR/prokka/ANNOT.gff | head -5

# View feature table
head -20 $OUTDIR/prokka/ANNOT.tsv

# Protein stats
seqkit stats $OUTDIR/prokka/ANNOT.faa

# Look at first 3 protein sequences
seqkit head -n 3 $OUTDIR/prokka/ANNOT.faa

# Find a specific protein by keyword
grep "ribosomal" $OUTDIR/prokka/ANNOT.tsv
```

## Explore DIAMOND results

```bash
# View first 10 hits
head -10 $OUTDIR/diamond/blastp_results.tsv | column -t

# Count proteins with hits
cut -f1 $OUTDIR/diamond/blastp_results.tsv | sort -u | wc -l

# Best hit per protein (highest bitscore)
sort -k1,1 -k12,12rn $OUTDIR/diamond/blastp_results.tsv | sort -k1,1 -u | head -10 | column -t

# Distribution of percent identity (best hits)
sort -k1,1 -k12,12rn $OUTDIR/diamond/blastp_results.tsv | sort -k1,1 -u | \
  awk -F'\t' '{printf "%d\n", $3}' | sort -n | uniq -c | sort -k2 -n

# Search for a specific protein
grep "ANNOT_00001" $OUTDIR/diamond/blastp_results.tsv | column -t

# Find hits to a keyword (e.g., "kinase")
grep -i "kinase" $OUTDIR/diamond/blastp_results.tsv | head -5 | column -t
```

## Explore hmmscan results

```bash
# View first 10 domain hits (skip comment lines)
grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | head -10

# Count total domain hits
grep -cv "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt

# Top 20 most common domains
grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $1}' | sort | uniq -c | sort -rn | head -20

# Proteins with the most domains
grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $4}' | sort | uniq -c | sort -rn | head -10

# Find hits for a specific protein
grep "ANNOT_00042" $OUTDIR/hmmer/hmmscan_domtblout.txt

# Find a specific domain (e.g., ABC transporter)
grep -i "ABC" $OUTDIR/hmmer/hmmscan_domtblout.txt | head -5
```

## Combined analysis

```bash
# Proteins with DIAMOND hits but no Pfam domains
comm -23 \
  <(cut -f1 $OUTDIR/diamond/blastp_results.tsv | sort -u) \
  <(grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $4}' | sort -u) \
  | head -10

# Proteins with Pfam domains but no DIAMOND hits
comm -13 \
  <(cut -f1 $OUTDIR/diamond/blastp_results.tsv | sort -u) \
  <(grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $4}' | sort -u) \
  | head -10

# Proteins with NO annotation at all
comm -23 \
  <(grep "^>" $OUTDIR/prokka/ANNOT.faa | sed 's/>//' | awk '{print $1}' | sort) \
  <(cat <(cut -f1 $OUTDIR/diamond/blastp_results.tsv) \
        <(grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $4}') | sort -u) \
  | wc -l
```

## Useful seqkit commands

```bash
# Length distribution of proteins
seqkit fx2tab -l $OUTDIR/prokka/ANNOT.faa | awk '{print $NF}' | sort -n | tail -5

# Find the longest protein
seqkit fx2tab -l $OUTDIR/prokka/ANNOT.faa | sort -t$'\t' -k2 -rn | head -1

# Find the shortest protein
seqkit fx2tab -l $OUTDIR/prokka/ANNOT.faa | sort -t$'\t' -k2 -n | head -1

# Extract a specific protein by ID
seqkit grep -p "ANNOT_00042" $OUTDIR/prokka/ANNOT.faa

# Count proteins shorter than 100 aa
seqkit fx2tab -l $OUTDIR/prokka/ANNOT.faa | awk '$NF < 100' | wc -l
```
