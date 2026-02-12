# Troubleshooting

Common issues and fixes.

## Environment

### "command not found: prokka" (or any tool)

**Cause**: The conda environment is not activated.

**Fix**:
```bash
micromamba activate ngs-annotation
# or
conda activate ngs-annotation
```

If the environment does not exist:
```bash
bash env/micromamba_create.sh
```

### "conda: command not found"

**Cause**: No package manager installed.

**Fix**: Install micromamba:
```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Then restart your shell and try again.

### Environment creation fails with "solving environment" forever

**Cause**: Conda solver is slow.

**Fix**: Use micromamba (preferred) or mamba instead of conda:
```bash
# Install mamba into base
conda install -n base -c conda-forge mamba
mamba env create -f env/environment.yml
```

### "PackagesNotFoundError" during env creation

**Cause**: Channel ordering issue or outdated channel data.

**Fix**:
```bash
# Update channel data
conda update --all -n base
# Or try with explicit channels
micromamba create -n ngs-annotation -c conda-forge -c bioconda -c defaults prokka seqkit diamond hmmer csvkit pigz wget curl
```

## Data and databases

### "Assembly not found"

**Fix**: Run the data download script:
```bash
bash scripts/01_get_data.sh
```

Or check that the assembly is at `data/input/assembly.fasta`.

### "DIAMOND database not found"

**Fix**:
```bash
bash scripts/01_get_data.sh    # get proteins first
bash scripts/02_make_diamond_db.sh
```

### "HMM database not pressed (missing .h3m)"

**Cause**: `hmmpress` was not run or failed.

**Fix**:
```bash
hmmpress -f data/db/pfam_subset.hmm
```

### wget fails (no internet)

**Cause**: HPC has no outbound internet.

**Fix**: Download on a machine with internet and transfer:
```bash
# On machine with internet
bash scripts/01_get_data.sh
bash scripts/02_make_diamond_db.sh
bash scripts/03_get_pfam.sh --subset

# Transfer to HPC
rsync -avP data/ hpc:/path/to/blastim/data/
```

### Pfam download fails ("404 Not Found")

**Cause**: Pfam URL may have changed (EBI reorganizes periodically).

**Fix**: Check the current URL at https://www.ebi.ac.uk/interpro/download/Pfam/ and update the URL in `scripts/03_get_pfam.sh`.

## Running scripts

### "Config not found"

**Fix**:
```bash
cp workflow/config.sh.example workflow/config.sh
```

### "Permission denied" when running scripts

**Fix**:
```bash
chmod +x scripts/*.sh workflow/*.sh
```

Or run with `bash`:
```bash
bash scripts/10_run_prokka.sh
```

### Prokka: "Could not run prodigal"

**Cause**: Prokka dependency issue.

**Fix**:
```bash
# Check prodigal is available
which prodigal
prodigal -v

# If missing, install it
micromamba install -n ngs-annotation -c bioconda prodigal
```

### Prokka: "tbl2asn is missing or too old"

**Cause**: Common Prokka issue with newer systems.

**Fix**: This is a warning, not an error. Prokka still produces GFF/FAA/FFN correctly. The .gbk and .sqn files may be affected but are not needed for this course.

### DIAMOND: "Database was made with incompatible version"

**Cause**: DIAMOND version mismatch between db build and search.

**Fix**: Rebuild the database:
```bash
rm data/db/teaching_proteins.dmnd
bash scripts/02_make_diamond_db.sh
```

### hmmscan runs forever

**Cause**: Using full Pfam-A (~20k profiles) instead of the subset.

**Fix**: Switch to subset mode:
```bash
# Edit workflow/config.sh
HMM_DB="data/db/pfam_subset.hmm"

# Or limit proteins
bash scripts/30_run_hmmscan.sh --max-proteins 200
```

### hmmscan: "No such file or directory" for .h3m

**Cause**: The HMM database was not pressed.

**Fix**:
```bash
hmmpress -f data/db/pfam_subset.hmm
```

### "Disk quota exceeded"

**Cause**: HPC disk quota reached.

**Fix**:
```bash
# Check disk usage
du -sh data/ outputs/

# Clean old outputs
make clean

# Check quota
quota -s  # or: df -h .
```

## Output issues

### 0 DIAMOND hits

**Possible causes**:
1. Database is empty or corrupt: `diamond dbinfo --db data/db/teaching_proteins`
2. Query file is empty: `wc -l outputs/run_*/prokka/ANNOT.faa`
3. E-value threshold too strict (unlikely with default 1e-5)

### 0 hmmscan hits

**Possible causes**:
1. HMM database has 0 profiles: `grep -c "^ACC" data/db/pfam_subset.hmm`
2. Database not pressed: check for `.h3m` files
3. E-value too strict (unlikely with defaults)

### Sanity checks show unexpected numbers

Compare your results with `docs/expected_outputs.md`. If numbers are off by more than 50%, check:
1. Is the correct assembly being used?
2. Are databases populated?
3. Are tool versions roughly matching the environment.yml specs?

## General tips

- Always activate the environment before running scripts.
- Run `bash scripts/00_check_system.sh` to diagnose most issues.
- Check the script output messages -- they indicate what to run next.
- If a step fails, fix the issue and re-run that step only (not the whole workflow).
