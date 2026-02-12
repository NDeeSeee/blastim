# Database directory

This directory stores databases used during the course.
Files here are **not tracked by git** (too large).

## Contents after instructor prep

| File | Built by | Description |
|------|----------|-------------|
| `teaching_proteins.fasta` | `01_get_data.sh` | Small protein set (~100k seqs) |
| `teaching_proteins.dmnd` | `02_make_diamond_db.sh` | DIAMOND index of above |
| `pfam_subset.hmm` | `03_get_pfam.sh --subset` | ~200 HMM profiles for fast hmmscan |
| `pfam_subset.hmm.h3{m,i,f,p}` | `03_get_pfam.sh --subset` | hmmpress index files |
| `Pfam-A.hmm` (optional) | `03_get_pfam.sh --full` | Full Pfam-A (~20k profiles) |

## How to populate

Run the instructor prep scripts from the repo root:

```bash
bash scripts/01_get_data.sh
bash scripts/02_make_diamond_db.sh
bash scripts/03_get_pfam.sh --subset   # fast class mode
# bash scripts/03_get_pfam.sh --full   # full Pfam (optional, slower)
```
