#!/usr/bin/env bash
set -euo pipefail

# ── 00_check_system.sh ────────────────────────────────────────────────
# Проверка системы перед запуском практикума по интерпретации вариантов.
#
# Usage:
#   bash scripts/00_check_system.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

echo "╔══════════════════════════════════════════════════════╗"
echo "║  Проверка системы: Интерпретация вариантов           ║"
echo "╚══════════════════════════════════════════════════════╝"
echo ""

PASS=0
FAIL=0

check_cmd() {
    local cmd="$1"
    local desc="${2:-$cmd}"
    if command -v "$cmd" &>/dev/null; then
        local ver
        ver=$( "$cmd" --version 2>&1 | head -1 || echo "ok" )
        echo "  ✓  ${desc}: ${ver}"
        PASS=$(( PASS + 1 ))
    else
        echo "  ✗  ${desc}: НЕ НАЙДЕН"
        FAIL=$(( FAIL + 1 ))
    fi
}

check_file() {
    local path="$1"
    local desc="$2"
    if [ -f "${REPO_ROOT}/${path}" ]; then
        echo "  ✓  ${desc}: ${path}"
        PASS=$(( PASS + 1 ))
    else
        echo "  ✗  ${desc}: файл не найден (${path})"
        FAIL=$(( FAIL + 1 ))
    fi
}

echo "── Инструменты ─────────────────────────────────────────"
check_cmd "python3"           "Python 3"
check_cmd "vep"               "Ensembl VEP"
check_cmd "bcftools"          "bcftools"
check_cmd "samtools"          "samtools"
check_cmd "bgzip"             "bgzip (htslib)"
check_cmd "tabix"             "tabix (htslib)"
check_cmd "curl"              "curl"
check_cmd "jq"                "jq"

echo ""
echo "── Python-пакеты ────────────────────────────────────────"
for pkg in pandas matplotlib; do
    if python3 -c "import ${pkg}" 2>/dev/null; then
        ver=$(python3 -c "import ${pkg}; print(${pkg}.__version__)" 2>/dev/null || echo "ok")
        echo "  ✓  python: ${pkg} ${ver}"
        PASS=$(( PASS + 1 ))
    else
        echo "  ✗  python: ${pkg} НЕ НАЙДЕН (pip install ${pkg})"
        FAIL=$(( FAIL + 1 ))
    fi
done

echo ""
echo "── Входные данные ───────────────────────────────────────"
check_file "data/input/demo_variants.vcf"     "demo VCF (SNV/indel)"
check_file "data/input/demo_sv.vcf"           "demo SV VCF"
check_file "data/input/demo_variants_vep.tsv" "предобработанный VEP TSV"
check_file "workflow/config.sh"               "config.sh"

echo ""
echo "── Конфигурация ─────────────────────────────────────────"
CONFIG="${REPO_ROOT}/workflow/config.sh"
if [ -f "$CONFIG" ]; then
    source "$CONFIG"
    echo "  VEP_MODE=${VEP_MODE:-не задан}"
    if [ "${VEP_MODE:-rest}" = "offline" ]; then
        if [ -d "${VEP_CACHE_DIR:-/shared/vep_cache}" ]; then
            echo "  ✓  VEP cache: ${VEP_CACHE_DIR}"
            PASS=$(( PASS + 1 ))
        else
            echo "  ✗  VEP cache: не найден (${VEP_CACHE_DIR:-/shared/vep_cache})"
            echo "       Запустите: make vep-cache"
            FAIL=$(( FAIL + 1 ))
        fi
    else
        echo "  (режим REST: cache не нужен)"
    fi
else
    echo "  WARN: workflow/config.sh не найден."
    echo "        Выполните: cp workflow/config.sh.example workflow/config.sh"
fi

echo ""
echo "────────────────────────────────────────────────────────"
echo "  Итог: ${PASS} ОК, ${FAIL} ошибок"
echo ""

if [ "$FAIL" -gt 0 ]; then
    echo "WARN: Исправьте ошибки перед запуском пайплайна."
    echo "      Документация: docs/troubleshooting.md"
    exit 1
else
    echo "Система готова. Следующий шаг: bash scripts/10_run_vep.sh"
fi
