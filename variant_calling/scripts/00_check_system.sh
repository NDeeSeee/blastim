#!/usr/bin/env bash
set -euo pipefail

# ── 00_check_system.sh ──────────────────────────────────────────────
# Verify that the HPC environment has everything we need.
# Run this BEFORE the class to catch issues early.
# ────────────────────────────────────────────────────────────────────

PASS=0
FAIL=0

ok()   { echo "  [OK]   $1"; PASS=$(( PASS + 1 )); }
fail() { echo "  [FAIL] $1"; FAIL=$(( FAIL + 1 )); }
warn() { echo "  [WARN] $1"; }

echo "======================================"
echo " System check for NGS Variants class"
echo "======================================"
echo ""

# ── OS ──
echo "── Operating system ──"
if [[ "$(uname -s)" == "Linux" ]]; then
    ok "Linux detected ($(uname -r))"
else
    warn "Expected Linux, got $(uname -s)"
fi

if [[ "$(uname -m)" == "x86_64" ]]; then
    ok "x86_64 architecture"
else
    warn "Expected x86_64, got $(uname -m)"
fi

# ── CPU ──
echo ""
echo "── CPUs ──"
NCPU=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)
if (( NCPU >= 4 )); then
    ok "$NCPU CPUs available (need >= 4)"
else
    warn "$NCPU CPUs available (recommend >= 4)"
fi

# ── Disk space ──
echo ""
echo "── Disk space ──"
AVAIL_KB=$(df -Pk . | awk 'NR==2{print $4}')
AVAIL_GB=$(( AVAIL_KB / 1048576 ))
if (( AVAIL_GB >= 2 )); then
    ok "${AVAIL_GB} GB free in current directory (need >= 2 GB)"
else
    fail "${AVAIL_GB} GB free in current directory (need >= 2 GB)"
fi

# ── Required tools ──
echo ""
echo "── Required tools ──"
for tool in samtools bcftools freebayes bwa bgzip tabix; do
    if command -v "$tool" &>/dev/null; then
        ver=$("$tool" --version 2>&1 | head -1 || echo "unknown")
        ok "$tool  ($ver)"
    else
        fail "$tool not found in PATH"
    fi
done

# ── Optional tools ──
echo ""
echo "── Optional tools ──"
for tool in seqkit python3 rg; do
    if command -v "$tool" &>/dev/null; then
        ok "$tool available"
    else
        warn "$tool not found (optional)"
    fi
done

# ── Conda / micromamba ──
echo ""
echo "── Package manager ──"
if command -v micromamba &>/dev/null; then
    ok "micromamba available"
elif command -v mamba &>/dev/null; then
    ok "mamba available"
elif command -v conda &>/dev/null; then
    ok "conda available"
else
    fail "No conda/mamba/micromamba found"
fi

# ── Data files ──
echo ""
echo "── Data files ──"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [ -f "$REPO_ROOT/data/input/reference.fasta" ]; then
    SIZE=$(du -h "$REPO_ROOT/data/input/reference.fasta" | cut -f1)
    ok "Reference genome found ($SIZE)"
else
    warn "data/input/reference.fasta not found (run 01_get_data.sh)"
fi

if [ -f "$REPO_ROOT/data/input/reads_sorted.bam" ]; then
    SIZE=$(du -h "$REPO_ROOT/data/input/reads_sorted.bam" | cut -f1)
    ok "Sorted BAM found ($SIZE)"
else
    warn "data/input/reads_sorted.bam not found (run 01_get_data.sh)"
fi

if [ -f "$REPO_ROOT/data/input/reads_sorted.bam.bai" ]; then
    ok "BAM index found"
else
    warn "BAM index not found (run 01_get_data.sh)"
fi

if [ -f "$REPO_ROOT/data/input/reference.fasta.fai" ]; then
    ok "Reference index (.fai) found"
else
    warn "Reference .fai not found (run 10_prep_reference.sh)"
fi

# ── Summary ──
echo ""
echo "======================================"
echo " Results:  $PASS passed,  $FAIL failed"
echo "======================================"

if (( FAIL > 0 )); then
    echo ""
    echo "Fix the FAIL items above before class."
    exit 1
else
    echo ""
    echo "All checks passed. Ready for class!"
fi
