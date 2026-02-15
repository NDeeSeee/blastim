# ── Makefile for Variant Calling Course ────────────────────────────
# Usage:
#   make env          Create conda/micromamba environment
#   make check        Run system checks
#   make prep-data    Download reference + simulate reads + align
#   make call         Run variant calling (freebayes)
#   make filter       Filter variants (bcftools)
#   make stats        Compute variant statistics
#   make all          Run complete workflow (steps 10→60)
#   make clean        Remove outputs (keep input data)
#   make clean-all    Remove outputs AND input data
#   make homework-check  Self-check HW2 submission
# ──────────────────────────────────────────────────────────────────

SHELL := /bin/bash

.PHONY: env prep-data call filter stats all clean clean-all check homework-check help

help:
	@echo "Available targets:"
	@echo "  make env             Create conda/micromamba environment"
	@echo "  make check           Run system checks"
	@echo "  make prep-data       Download reference + simulate reads + align"
	@echo "  make call            Run freebayes variant calling (steps 10-30)"
	@echo "  make filter          Filter variants with bcftools (step 40)"
	@echo "  make stats           Compute variant statistics (step 50)"
	@echo "  make all             Run complete workflow (steps 10-60)"
	@echo "  make clean           Remove outputs (keep input data)"
	@echo "  make clean-all       Remove outputs AND input data"
	@echo "  make homework-check  Self-check HW2 submission"

env:
	bash env/micromamba_create.sh

check:
	bash scripts/00_check_system.sh

prep-data:
	bash scripts/01_get_data.sh

call:
	bash scripts/10_prep_reference.sh
	bash scripts/11_prep_bam.sh
	bash scripts/20_call_freebayes.sh
	bash scripts/30_compress_index_vcf.sh

filter:
	bash scripts/40_filter_bcftools.sh

stats:
	bash scripts/50_stats_bcftools.sh
	bash scripts/60_pick_3_variants.sh

all:
	bash workflow/run_all.sh

clean:
	rm -rf outputs/run_*
	@echo "Cleaned outputs. Input data in data/input/ preserved."

homework-check:
	@test -n "$$GITHUB_USER" || (echo "ERROR: Set GITHUB_USER first: export GITHUB_USER=<your_github>"; exit 1)
	python scripts/91_homework_check.py --user $$GITHUB_USER

clean-all: clean
	rm -f data/input/reference.fasta data/input/reference.fasta.*
	rm -f data/input/reads_sorted.bam data/input/reads_sorted.bam.bai
	rm -f data/input/sim_R1.fq data/input/sim_R2.fq
	rm -f data/truth/wgsim.mutations.txt
	@echo "Cleaned outputs and input data."
