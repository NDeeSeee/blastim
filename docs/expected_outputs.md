# Expected Outputs

What files should appear and what "good" looks like after a successful run with B. subtilis 168.

## Output directory structure

After `bash workflow/run_all.sh`:

```
outputs/run_YYYYMMDD_HHMMSS/
├── prokka/
│   ├── ANNOT.gff          ~12-15 MB   Annotated genome (GFF3)
│   ├── ANNOT.faa          ~1.5 MB     Predicted proteins (FASTA)
│   ├── ANNOT.ffn          ~3.5 MB     Nucleotide CDS (FASTA)
│   ├── ANNOT.fna          ~4.2 MB     Input contigs (renamed)
│   ├── ANNOT.fsa          ~4.2 MB     Submission-ready contigs
│   ├── ANNOT.gbk          ~15 MB      GenBank format
│   ├── ANNOT.log          varies      Prokka log
│   ├── ANNOT.sqn          ~15 MB      Sequin format
│   ├── ANNOT.tbl          ~2 MB       Feature table (Sequin)
│   ├── ANNOT.tsv          ~500 KB     Tab-separated features
│   └── ANNOT.txt          <1 KB       Summary statistics
├── diamond/
│   └── blastp_results.tsv ~5-20 MB    DIAMOND blastp hits
└── hmmer/
    ├── hmmscan_domtblout.txt  ~1-5 MB Domain hits table
    ├── hmmscan.log            varies  hmmscan full output
    └── query_subset.faa       (only if --max-proteins was used)
```

## Expected counts (B. subtilis 168, approximate)

### Prokka

| Feature type | Expected count |
|-------------|----------------|
| CDS | 4,200 - 4,400 |
| tRNA | 85 - 90 |
| rRNA | 10 - 30 |
| tmRNA | 1 |
| Total genes | ~4,300 - 4,500 |

### DIAMOND

| Metric | Expected |
|--------|----------|
| Total hits (all per-query, up to 5) | 15,000 - 20,000 |
| Unique proteins with hits | 3,000 - 4,000 |
| Hit rate | 70% - 90% |
| Median % identity (best hits) | 40% - 60% |

### hmmscan (subset, ~200 HMMs)

| Metric | Expected |
|--------|----------|
| Total domain hits | 1,000 - 3,000 |
| Proteins with domains | 800 - 2,000 |
| Unique domains found | 80 - 180 |

### hmmscan (full Pfam-A)

| Metric | Expected |
|--------|----------|
| Total domain hits | 4,000 - 6,000 |
| Proteins with domains | 2,500 - 3,500 |
| Unique domains found | 800 - 1,500 |

### Evidence ladder

| Category | Expected (subset Pfam) |
|----------|----------------------|
| Both DIAMOND + domain | 800 - 2,000 |
| DIAMOND only | 1,000 - 2,500 |
| Domain only | 100 - 500 |
| No annotation (hypothetical) | 300 - 1,000 |

## What "good" looks like

- Prokka runs without errors and finds >4,000 CDS.
- DIAMOND finds hits for 70%+ of proteins.
- hmmscan (subset) identifies domains in 20-50% of proteins.
- The sanity check script shows a plausible protein length distribution (median ~250-350 aa).
- The evidence ladder shows a gradient: most proteins have at least some evidence.
- Very few proteins (< 5%) should be extremely short (< 50 aa) -- if there are many, Prokka may have used overly permissive ORF calling.

## Red flags

- **0 CDS predicted**: Assembly file might be empty or in wrong format.
- **No DIAMOND hits at all**: Database may not have been built correctly. Check `diamond dbinfo`.
- **hmmscan exits immediately**: HMM database not pressed. Run `hmmpress`.
- **Hit rate < 50%**: Database may be too small or wrong organism group. Check Swiss-Prot subset.
- **Prokka reports "0 contigs"**: Input file may have Windows line endings. Fix with `dos2unix`.
