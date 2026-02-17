#!/usr/bin/env bash
set -euo pipefail

# ── 10_run_vep.sh ─────────────────────────────────────────────────────
# Шаг 10: Аннотация вариантов с помощью Ensembl VEP.
#
# Два режима (задаётся в workflow/config.sh):
#   VEP_MODE=rest     — VEP REST API (нужен интернет, не нужен cache)
#   VEP_MODE=offline  — локальный VEP с cache (~15 ГБ, нужен заранее)
#
# Если VEP не установлен или нет интернета, скрипт автоматически
# использует предобработанный файл: data/input/demo_variants_vep.tsv
#
# Usage:
#   bash scripts/10_run_vep.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG="${REPO_ROOT}/workflow/config.sh"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: workflow/config.sh не найден."
    echo "       cp workflow/config.sh.example workflow/config.sh"
    exit 1
fi

source "$CONFIG"

# Разрешаем относительные пути
DEMO_VCF="${REPO_ROOT}/${DEMO_VCF}"
PRECOMPUTED_VEP="${REPO_ROOT}/${PRECOMPUTED_VEP}"
OUTDIR="${REPO_ROOT}/${OUTDIR}"

mkdir -p "${OUTDIR}/vep"
VEP_OUT="${OUTDIR}/vep/vep_output.tsv"
VEP_MODE="${VEP_MODE:-rest}"

echo "==> Шаг 10: VEP аннотация"
echo "    VCF: ${DEMO_VCF}"
echo "    Режим: ${VEP_MODE}"
echo "    Выход: ${VEP_OUT}"
echo ""

START=$(date +%s)

# ── Вспомогательная функция ──────────────────────────────────────────
copy_precomputed() {
    echo "    Используем предобработанный файл: ${PRECOMPUTED_VEP}"
    cp "$PRECOMPUTED_VEP" "$VEP_OUT"
    echo "    ✓ Скопировано: ${VEP_OUT}"
}

# ── Режим REST ───────────────────────────────────────────────────────
run_vep_rest() {
    echo "==> Режим REST: отправка вариантов в Ensembl API..."
    python3 "${SCRIPT_DIR}/_vep_rest_helper.py" \
        --vcf  "$DEMO_VCF" \
        --out  "$VEP_OUT" \
        --assembly "${VEP_ASSEMBLY:-GRCh38}"
}

# ── Режим offline ────────────────────────────────────────────────────
run_vep_offline() {
    if [ ! -d "${VEP_CACHE_DIR:-/shared/vep_cache}" ]; then
        echo "  WARN: VEP cache не найден: ${VEP_CACHE_DIR}"
        echo "        Переключаемся на предобработанный файл."
        copy_precomputed
        return
    fi

    echo "==> Режим offline: локальный VEP с cache..."
    vep \
        --input_file  "$DEMO_VCF" \
        --output_file "$VEP_OUT" \
        --format      vcf \
        --tab \
        --cache \
        --offline \
        --dir_cache   "${VEP_CACHE_DIR}" \
        --assembly    "${VEP_ASSEMBLY:-GRCh38}" \
        --symbol \
        --canonical \
        --biotype \
        --hgvs \
        --af_gnomade \
        --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,Extra" \
        --force_overwrite \
        --fork "${CPUS:-4}"
}

# ── Выбор режима ─────────────────────────────────────────────────────
if [ "${VEP_MODE}" = "offline" ]; then
    run_vep_offline
elif [ "${VEP_MODE}" = "rest" ]; then
    if command -v python3 &>/dev/null && python3 -c "import urllib.request" 2>/dev/null; then
        run_vep_rest || {
            echo "  WARN: REST API недоступен. Используем предобработанный файл."
            copy_precomputed
        }
    else
        copy_precomputed
    fi
else
    echo "  WARN: Неизвестный VEP_MODE='${VEP_MODE}'. Используем предобработанный файл."
    copy_precomputed
fi

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> VEP аннотация завершена за ${ELAPSED}с"
echo "    Результат: ${VEP_OUT}"
echo ""
echo "    Следующий шаг: bash scripts/20_prioritize_variants.sh"
