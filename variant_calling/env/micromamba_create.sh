#!/usr/bin/env bash
set -euo pipefail

# Create the ngs-variants environment using micromamba (fast) or conda (fallback).
# Usage: bash env/micromamba_create.sh

ENV_NAME="ngs-variants"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
ENV_FILE="${REPO_ROOT}/env/environment.yml"

if ! [ -f "$ENV_FILE" ]; then
    echo "ERROR: Cannot find $ENV_FILE"
    exit 1
fi

# ── Try micromamba first, then mamba, then conda ──
if command -v micromamba &>/dev/null; then
    SOLVER=micromamba
elif command -v mamba &>/dev/null; then
    SOLVER=mamba
elif command -v conda &>/dev/null; then
    SOLVER=conda
else
    echo "ERROR: None of micromamba / mamba / conda found in PATH."
    echo ""
    echo "Install micromamba (recommended):"
    echo '  "${SHELL}" <(curl -L micro.mamba.pm/install.sh)'
    exit 1
fi

echo "==> Using solver: $SOLVER"
echo "==> Creating environment '$ENV_NAME' from $ENV_FILE ..."

if [ "$SOLVER" = "micromamba" ]; then
    micromamba create -n "$ENV_NAME" -f "$ENV_FILE" -y
    echo ""
    echo "Activate with:  micromamba activate $ENV_NAME"
else
    $SOLVER env create -n "$ENV_NAME" -f "$ENV_FILE" -y
    echo ""
    echo "Activate with:  conda activate $ENV_NAME"
fi

echo "==> Done."
