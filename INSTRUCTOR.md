# Instructor Guide

How to prepare the HPC, stage data, and run this course.

## Overview

| Phase | What | When |
|-------|------|------|
| **Prep** | Install env, download data/dbs | 1-2 days before class |
| **Validate** | Smoke test the full workflow | Day before class |
| **Class** | Students run steps 10-50 | During the 2h session |
| **Reset** | Clean outputs for next cohort | After class |

## 1. Disk requirements

| Component | Size |
|-----------|------|
| Conda environment | ~2 GB |
| Assembly FASTA | ~4 MB |
| Teaching proteins (FASTA) | ~200 MB |
| DIAMOND database (.dmnd) | ~50 MB |
| Pfam-A full (optional) | ~1.5 GB |
| Pfam subset (~200 HMMs) | ~5 MB |
| Per-student outputs | ~50 MB |
| **Total (subset mode)** | **~2.5 GB** |
| **Total (full Pfam mode)** | **~4 GB** |

## 2. Environment setup

```bash
# Option A: micromamba (recommended)
bash env/micromamba_create.sh
micromamba activate ngs-annotation

# Option B: conda/mamba
conda env create -f env/environment.yml
conda activate ngs-annotation
```

Verify the environment:

```bash
bash scripts/00_check_system.sh
```

All items should show `[OK]` except data/db warnings (those come next).

## 3. Stage data and databases

Run from the repo root:

```bash
# Download B. subtilis genome + Swiss-Prot protein subset
bash scripts/01_get_data.sh

# Build DIAMOND database
bash scripts/02_make_diamond_db.sh

# Download Pfam and create subset for fast in-class runs
bash scripts/03_get_pfam.sh --subset
```

Or all at once:

```bash
make prep-data
```

This requires internet access. After completion, all data is in `data/` and the repo can run offline.

### Offline class mode

If the HPC has no internet during class:

1. Run all prep scripts on a machine with internet.
2. Copy the entire `data/` directory to the HPC:
   ```bash
   rsync -avP data/ hpc:/path/to/blastim/data/
   ```
3. Students only need the repo + the pre-staged `data/` directory.

### Optional: Full Pfam

For advanced students or if time permits:

```bash
bash scripts/03_get_pfam.sh --full
```

Then edit `workflow/config.sh`:
```bash
HMM_DB="data/db/Pfam-A.hmm"
```

Warning: hmmscan with full Pfam takes 15-30 min on 8 CPUs for this genome.

## 4. Smoke test

Run the complete workflow once to verify everything works:

```bash
cp workflow/config.sh.example workflow/config.sh
bash workflow/run_all.sh
```

### Expected timings (8 CPUs, subset Pfam)

| Step | Expected time |
|------|---------------|
| Prokka | 3-8 min |
| DIAMOND | 1-2 min |
| hmmscan (subset) | 2-5 min |
| hmmscan (full Pfam) | 15-30 min |
| Sanity checks | <1 min |
| **Total (subset)** | **~10-20 min** |

### Expected output counts (approximate, B. subtilis 168)

| Metric | Expected range |
|--------|---------------|
| Predicted CDS | 4,200-4,400 |
| Predicted tRNA | 80-90 |
| Predicted rRNA | 10-30 |
| Proteins with DIAMOND hits | 3,000-4,000 |
| Proteins with Pfam domains (subset) | 1,500-2,500 |
| Hypothetical proteins | 200-800 |

If numbers differ significantly, check tool versions and database contents.

## 5. Student workflow

Students should:

1. Copy the config file:
   ```bash
   cp workflow/config.sh.example workflow/config.sh
   ```

2. Run either the full workflow:
   ```bash
   bash workflow/run_all.sh
   ```

   Or step-by-step (recommended for teaching):
   ```bash
   bash scripts/10_run_prokka.sh
   bash scripts/20_run_diamond.sh
   bash scripts/30_run_hmmscan.sh
   bash scripts/40_sanity_checks.sh
   bash scripts/50_pick_3_genes.sh
   ```

3. Explore outputs using commands from `docs/cheatsheet.md`.

## 6. Plan B: hmmscan too slow

If hmmscan takes too long during class:

### Option A: Limit proteins
```bash
bash scripts/30_run_hmmscan.sh --max-proteins 200
```
This scans only the first 200 proteins. Runs in ~1 min.

### Option B: Use subset HMMs only
Ensure `config.sh` has:
```bash
HMM_DB="data/db/pfam_subset.hmm"
```

### Option C: Skip hmmscan entirely
Students can still learn from Prokka + DIAMOND results. Run steps 10, 20, 40, 50 only.

## 7. Reset between cohorts

```bash
make clean
```

This removes all `outputs/run_*` directories but preserves databases.

For a full reset (including databases):
```bash
make clean-all
```

## 8. Suggested class timeline (2 hours)

| Time | Activity |
|------|----------|
| 0:00-0:10 | Introduction: what is genome annotation? Evidence ladder concept |
| 0:10-0:15 | SSH in, activate environment, copy config |
| 0:15-0:25 | Run Prokka; while waiting: explain GFF format |
| 0:25-0:35 | Explore Prokka outputs (GFF, FAA, TSV) |
| 0:35-0:45 | Run DIAMOND; explain BLAST tabular format |
| 0:45-0:55 | Explore DIAMOND results; discuss identity/E-value |
| 0:55-1:10 | Run hmmscan; explain HMM profiles and Pfam |
| 1:10-1:20 | Explore hmmscan results; top domains |
| 1:20-1:35 | Run sanity checks; discuss evidence ladder |
| 1:35-1:50 | Pick 3 genes exercise; class discussion |
| 1:50-2:00 | Wrap-up: real-world annotation, limitations, further resources |

## 9. Licensing notes

- **Swiss-Prot**: Downloaded from UniProt (CC BY 4.0). Not committed to the repo.
- **Pfam**: Downloaded from EBI/InterPro. Not committed to the repo.
- **B. subtilis 168 genome**: NCBI RefSeq, public domain.
- **Course materials**: MIT license.

Do not commit full databases to the git repository. They are gitignored.

## 10. Troubleshooting

See `docs/troubleshooting.md` for common issues and fixes.
