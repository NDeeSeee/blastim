#!/usr/bin/env bash
set -euo pipefail

# ── 01_get_data.sh ────────────────────────────────────────────────────
# Проверка входных данных практикума.
# Демо-данные уже включены в репозиторий (data/input/).
# Скрипт проверяет целостность файлов.
#
# Usage:
#   bash scripts/01_get_data.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${REPO_ROOT}/data/input"

echo "==> Проверка входных данных..."
echo "    Директория: ${INPUT_DIR}"
echo ""

check_vcf() {
    local file="$1"
    local desc="$2"
    if [ -f "$file" ]; then
        local lines
        lines=$(grep -cv "^#" "$file" 2>/dev/null || echo 0)
        echo "  ✓  ${desc}: ${lines} вариантов ($(basename "$file"))"
    else
        echo "  ✗  ${desc}: файл не найден: $(basename "$file")"
        echo "       Убедитесь, что вы склонировали репозиторий полностью:"
        echo "       git clone https://github.com/<org>/blastim_ngs"
        exit 1
    fi
}

check_vcf "${INPUT_DIR}/demo_variants.vcf"     "SNV/indel демо VCF"
check_vcf "${INPUT_DIR}/demo_sv.vcf"           "SV демо VCF"

# Проверка предобработанного VEP
VEP_TSV="${INPUT_DIR}/demo_variants_vep.tsv"
if [ -f "$VEP_TSV" ]; then
    lines=$(grep -cv "^#" "$VEP_TSV" 2>/dev/null || echo 0)
    echo "  ✓  Предобработанный VEP TSV: ${lines} записей"
else
    echo "  ✗  demo_variants_vep.tsv: не найден"
    echo "     Запустите: bash scripts/10_run_vep.sh"
fi

echo ""
echo "==> Все входные данные на месте."
echo "    Данные: ${INPUT_DIR}/"
echo ""
echo "    Следующий шаг: bash scripts/10_run_vep.sh"
echo "    Или для полного пайплайна: bash workflow/run_all.sh"
