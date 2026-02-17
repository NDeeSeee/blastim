#!/usr/bin/env bash
set -euo pipefail

# ── 02_get_vep_cache.sh ───────────────────────────────────────────────
# Загрузка VEP cache для режима offline.
# ТОЛЬКО для преподавателя / администратора HPC-сервера.
# Студентам этот шаг не нужен — по умолчанию используется VEP REST API.
#
# Требования:
#   ~15 ГБ дискового пространства для GRCh38 Ensembl v110
#
# Usage:
#   bash scripts/02_get_vep_cache.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG="${REPO_ROOT}/workflow/config.sh"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: workflow/config.sh не найден."
    echo "       Выполните: cp workflow/config.sh.example workflow/config.sh"
    echo "       Установите VEP_MODE=offline и VEP_CACHE_DIR=<путь>"
    exit 1
fi

source "$CONFIG"

VEP_ASSEMBLY="${VEP_ASSEMBLY:-GRCh38}"
VEP_VERSION="${VEP_VERSION:-110}"
CACHE_DIR="${VEP_CACHE_DIR:-/shared/vep_cache}"

echo "==> Загрузка VEP cache"
echo "    Ансамбль: ${VEP_ASSEMBLY}, версия Ensembl: ${VEP_VERSION}"
echo "    Путь назначения: ${CACHE_DIR}"
echo ""
echo "    ВНИМАНИЕ: Загрузка занимает ~15 ГБ и 20-40 минут."
echo "    Прогресс отображается ниже."
echo ""

mkdir -p "$CACHE_DIR"

START=$(date +%s)

# Загрузка VEP FASTA (CDS reference)
FASTA_URL="https://ftp.ensembl.org/pub/release-${VEP_VERSION}/fasta/homo_sapiens/dna/Homo_sapiens.${VEP_ASSEMBLY}.dna.toplevel.fa.gz"
FASTA_OUT="${CACHE_DIR}/Homo_sapiens.${VEP_ASSEMBLY}.dna.toplevel.fa.gz"

echo "==> [1/2] Загрузка FASTA (может занять >10 минут)..."
if [ -f "$FASTA_OUT" ]; then
    echo "    Уже существует: ${FASTA_OUT} — пропускаем."
else
    curl -L --progress-bar -o "$FASTA_OUT" "$FASTA_URL"
    echo "    ✓ FASTA загружена"
fi

# Загрузка VEP cache
CACHE_URL="https://ftp.ensembl.org/pub/release-${VEP_VERSION}/variation/indexed_vep_cache/homo_sapiens_vep_${VEP_VERSION}_${VEP_ASSEMBLY}.tar.gz"
CACHE_TAR="${CACHE_DIR}/homo_sapiens_vep_${VEP_VERSION}_${VEP_ASSEMBLY}.tar.gz"

echo ""
echo "==> [2/2] Загрузка VEP cache (~15 ГБ)..."
if [ -d "${CACHE_DIR}/homo_sapiens" ]; then
    echo "    Cache уже распакован: ${CACHE_DIR}/homo_sapiens/ — пропускаем."
else
    if [ ! -f "$CACHE_TAR" ]; then
        curl -L --progress-bar -o "$CACHE_TAR" "$CACHE_URL"
    fi
    echo "    Распаковка..."
    tar xzf "$CACHE_TAR" -C "$CACHE_DIR"
    echo "    ✓ Cache распакован"
fi

# Индексация FASTA
if [ ! -f "${FASTA_OUT%.gz}" ]; then
    echo ""
    echo "==> Распаковка и индексация FASTA..."
    bgzip -d -k "$FASTA_OUT" || gunzip -k "$FASTA_OUT"
    samtools faidx "${FASTA_OUT%.gz}"
    echo "    ✓ FASTA индексирована"
fi

END=$(date +%s)
ELAPSED=$(( END - START ))
MINUTES=$(( ELAPSED / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo "==> VEP cache готов (${MINUTES}m ${SECS}s)"
echo "    Путь: ${CACHE_DIR}"
echo ""
echo "    Обновите workflow/config.sh:"
echo "      VEP_MODE=\"offline\""
echo "      VEP_CACHE_DIR=\"${CACHE_DIR}\""
echo ""
echo "    Следующий шаг: bash scripts/10_run_vep.sh"
