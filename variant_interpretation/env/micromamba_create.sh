#!/usr/bin/env bash
set -euo pipefail

# ── micromamba_create.sh ─────────────────────────────────────────────
# Создание conda-окружения ngs-interpret через micromamba.
#
# Usage:
#   bash env/micromamba_create.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/environment.yml"
ENV_NAME="ngs-interpret"

echo "==> Создание окружения: ${ENV_NAME}"
echo "    Спецификация: ${ENV_FILE}"

if micromamba env list 2>/dev/null | grep -q "^${ENV_NAME}"; then
    echo "    Окружение '${ENV_NAME}' уже существует."
    echo "    Для пересоздания: micromamba env remove -n ${ENV_NAME} && bash env/micromamba_create.sh"
    exit 0
fi

if command -v micromamba &>/dev/null; then
    micromamba env create -n "${ENV_NAME}" -f "${ENV_FILE}" --yes
elif command -v mamba &>/dev/null; then
    mamba env create -n "${ENV_NAME}" -f "${ENV_FILE}"
elif command -v conda &>/dev/null; then
    conda env create -n "${ENV_NAME}" -f "${ENV_FILE}"
else
    echo "ERROR: micromamba/mamba/conda не найден."
    echo "       Установите micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
    exit 1
fi

echo ""
echo "==> Окружение '${ENV_NAME}' создано."
echo "    Активация:"
echo "      micromamba activate ${ENV_NAME}"
echo "    или: conda activate ${ENV_NAME}"
