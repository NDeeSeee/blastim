# ── Makefile for Prokaryotic Genome Annotation Course ──────────────
# Usage:
#   make env          Create conda/micromamba environment
#   make prep-data    Download genome + proteins + build databases
#   make prokka       Run Prokka structural annotation
#   make diamond      Run DIAMOND homology search
#   make pfam         Run hmmscan domain search
#   make all          Run complete workflow (steps 10→50)
#   make clean        Remove outputs (keep databases)
#   make clean-all    Remove outputs AND databases
#   make check        Run system checks
#   make homework-check  Self-check HW1 submission
# ──────────────────────────────────────────────────────────────────

SHELL := /bin/bash

.PHONY: env prep-data prokka diamond pfam all clean clean-all check homework-check help

help:
	@echo "Available targets:"
	@echo "  make env          Create conda/micromamba environment"
	@echo "  make check        Run system checks"
	@echo "  make prep-data    Download genome + proteins + build all databases"
	@echo "  make prokka       Run Prokka (step 10)"
	@echo "  make diamond      Run DIAMOND (step 20)"
	@echo "  make pfam         Run hmmscan (step 30)"
	@echo "  make all          Run complete workflow (steps 10-50)"
	@echo "  make clean            Remove outputs (keep databases)"
	@echo "  make clean-all        Remove outputs AND databases"
	@echo "  make homework-check   Self-check HW1 submission"

env:
	bash env/micromamba_create.sh

check:
	bash scripts/00_check_system.sh

prep-data:
	bash scripts/01_get_data.sh
	bash scripts/02_make_diamond_db.sh
	bash scripts/03_get_pfam.sh --subset

prokka:
	bash scripts/10_run_prokka.sh

diamond:
	bash scripts/20_run_diamond.sh

pfam:
	bash scripts/30_run_hmmscan.sh

all:
	bash workflow/run_all.sh

clean:
	rm -rf outputs/run_*
	@echo "Cleaned outputs. Databases in data/db/ preserved."

homework-check:
	@test -n "$$GITHUB_USER" || (echo "ERROR: Set GITHUB_USER first: export GITHUB_USER=<your_github>"; exit 1)
	python scripts/91_homework_check.py --user $$GITHUB_USER

clean-all: clean
	rm -f data/db/*.dmnd data/db/*.h3* data/db/*.hmm data/db/*.fasta data/db/*.ssi
	rm -f data/input/assembly.fasta
	@echo "Cleaned outputs and databases."
