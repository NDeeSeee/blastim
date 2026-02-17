#!/usr/bin/env bash
set -euo pipefail

# ── run_all.sh ───────────────────────────────────────────────────────
# Запуск полного пайплайна интерпретации вариантов (шаги 10→30).
#
# Usage:
#   bash workflow/run_all.sh
# ────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG="${REPO_ROOT}/workflow/config.sh"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Конфигурационный файл не найден."
    echo "       Выполните: cp workflow/config.sh.example workflow/config.sh"
    exit 1
fi

export OUTDIR="outputs/run_$(date +%Y%m%d_%H%M%S)"

echo "╔══════════════════════════════════════════════════════╗"
echo "║  Практикум: Интерпретация геномных вариантов         ║"
echo "║  Организм: Homo sapiens GRCh38                       ║"
echo "╚══════════════════════════════════════════════════════╝"
echo ""
echo "  Конфигурация:  $CONFIG"
echo "  Выходная папка: $OUTDIR"
echo "  Начало:        $(date)"
echo ""

TOTAL_START=$(date +%s)

echo "━━━ Шаг 10: 10_run_vep.sh ━━━"
bash "${REPO_ROOT}/scripts/10_run_vep.sh"
echo ""

echo "━━━ Шаг 20: 20_prioritize_variants.sh ━━━"
bash "${REPO_ROOT}/scripts/20_prioritize_variants.sh"
echo ""

echo "━━━ Шаг 21: 21_make_report.py ━━━"
source "$CONFIG"
ABS_OUTDIR="${REPO_ROOT}/${OUTDIR}"
python3 "${REPO_ROOT}/scripts/21_make_report.py" \
    --tsv "${ABS_OUTDIR}/vep/vep_output.tsv" \
    --vcf "${REPO_ROOT}/${DEMO_VCF}" \
    --out "${ABS_OUTDIR}/report.md"
echo "    Отчёт: ${ABS_OUTDIR}/report.md"
echo ""

echo "━━━ Шаг 30: 30_sv_cnv_intro.sh ━━━"
bash "${REPO_ROOT}/scripts/30_sv_cnv_intro.sh"
echo ""

TOTAL_END=$(date +%s)
ELAPSED=$(( TOTAL_END - TOTAL_START ))
MINUTES=$(( ELAPSED / 60 ))
SECS=$(( ELAPSED % 60 ))

echo "╔══════════════════════════════════════════════════════╗"
echo "║  Пайплайн завершён!                                  ║"
echo "║  Общее время: ${MINUTES}m ${SECS}s                               ║"
echo "╚══════════════════════════════════════════════════════╝"
echo ""
echo "Ключевые файлы:"
echo "  ${OUTDIR}/vep/vep_output.tsv     — Аннотация VEP"
echo "  ${OUTDIR}/prioritized.tsv        — Приоритизированные варианты"
echo "  ${OUTDIR}/report.md              — Итоговый отчёт"
echo ""
echo "Следующий шаг: изучите docs/cheatsheet.md"
