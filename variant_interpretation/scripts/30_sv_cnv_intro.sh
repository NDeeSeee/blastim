#!/usr/bin/env bash
set -euo pipefail

# ── 30_sv_cnv_intro.sh ────────────────────────────────────────────────
# Шаг 30: Мини-демо структурных вариантов (SV) и копийных чисел (CNV).
#
# Демонстрирует:
#   - Формат SV в VCF (SVTYPE, END, SVLEN, символьные ALT)
#   - bcftools stats для SV VCF
#   - Базовую интерпретацию DEL/DUP/INV
#   - ASCII-визуализацию распределения размеров SV
#
# Usage:
#   bash scripts/30_sv_cnv_intro.sh
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

SV_VCF="${REPO_ROOT}/${DEMO_SV_VCF}"
OUTDIR="${REPO_ROOT}/${OUTDIR}"
mkdir -p "${OUTDIR}/sv"

echo "==> Шаг 30: SV/CNV мини-демо"
echo "    SV VCF: ${SV_VCF}"
echo ""

START=$(date +%s)

# ── 30.1: Содержимое SV VCF ─────────────────────────────────────────
echo "━━━ 30.1 Структурные варианты в VCF ━━━"
echo ""
echo "Заголовки INFO для SV:"
grep "^##INFO" "$SV_VCF" | grep -E "SVTYPE|SVLEN|END|CIPOS|CIEND" || true
echo ""
echo "Варианты (без заголовка):"
grep -v "^#" "$SV_VCF" | \
    awk 'BEGIN{OFS="\t"; print "ID","CHROM","POS","SVTYPE","SVLEN","END","FILTER"}
         {
           svtype="-"; svlen="-"; end_pos="-"
           n=split($8,info,";")
           for(i=1;i<=n;i++){
             split(info[i],kv,"=")
             if(kv[1]=="SVTYPE") svtype=kv[2]
             if(kv[1]=="SVLEN")  svlen=kv[2]
             if(kv[1]=="END")    end_pos=kv[2]
           }
           print $3,$1,$2,svtype,svlen,end_pos,$7
         }' | column -t

echo ""

# ── 30.2: Подсчёт SV по типу ────────────────────────────────────────
echo "━━━ 30.2 Статистика по типам SV ━━━"
echo ""
grep -v "^#" "$SV_VCF" | \
    awk '{
        n=split($8,info,";")
        for(i=1;i<=n;i++){
            split(info[i],kv,"=")
            if(kv[1]=="SVTYPE") print kv[2]
        }
    }' | sort | uniq -c | sort -rn | \
    awk '{printf "  %-10s %d\n", $2, $1}'

echo ""

# ── 30.3: Интерпретация SV ──────────────────────────────────────────
echo "━━━ 30.3 Клиническая интерпретация SV ━━━"
echo ""
python3 - <<'PYTHON'
import re

SV_TYPES = {
    "DEL": ("Делеция",     "Потеря генетического материала. Гомозиготная DEL = loss-of-function."),
    "DUP": ("Дупликация",  "Увеличение копийного числа. Может нарушить дозозависимые гены."),
    "INV": ("Инверсия",    "Переворот сегмента ДНК. Может разрушить гены на границах."),
    "BND": ("Транслокация","Перестройка между хромосомами. Часто онкогенна (e.g. BCR-ABL)."),
    "INS": ("Инсерция",    "Вставка последовательности. MEI (мобильные элементы) — частая причина."),
}

for svtype, (name, desc) in SV_TYPES.items():
    print(f"  {svtype:4s} ({name})")
    print(f"       {desc}")
    print()
PYTHON

# ── 30.4: Размеры SV ────────────────────────────────────────────────
echo "━━━ 30.4 Размеры структурных вариантов ━━━"
echo ""
grep -v "^#" "$SV_VCF" | \
    awk '{
        svlen=0
        n=split($8,info,";")
        for(i=1;i<=n;i++){
            split(info[i],kv,"=")
            if(kv[1]=="SVLEN"){
                val=kv[2]+0
                if(val<0) val=-val
                svlen=val
            }
        }
        if(svlen>0) printf "%s\t%d\n", $3, svlen
    }' | sort -k2 -n | \
    awk '{
        if($2<1000)      cat="<1kb"
        else if($2<10000) cat="1-10kb"
        else if($2<100000) cat="10-100kb"
        else              cat=">100kb"
        printf "  %-12s %s  (%d bp)\n", $1, cat, $2
    }'

echo ""

# ── 30.5: Сохраняем аннотированный SV TSV ───────────────────────────
SV_OUT="${OUTDIR}/sv/sv_annotated.tsv"
echo "━━━ 30.5 Сохранение результатов ━━━"
{
    echo -e "ID\tCHROM\tPOS\tSVTYPE\tSVLEN\tEND\tFILTER\tGENE_REGION\tINTERPRETATION"
    grep -v "^#" "$SV_VCF" | \
        awk '{
            svtype="-"; svlen="-"; end_pos="-"; gene="-"
            n=split($8,info,";")
            for(i=1;i<=n;i++){
                split(info[i],kv,"=")
                if(kv[1]=="SVTYPE") svtype=kv[2]
                if(kv[1]=="SVLEN")  svlen=kv[2]
                if(kv[1]=="END")    end_pos=kv[2]
                if(kv[1]=="GENE")   gene=kv[2]
            }
            interp="требует анализа"
            if(svtype=="DEL") interp="потеря материала"
            if(svtype=="DUP") interp="увеличение копийного числа"
            if(svtype=="INV") interp="инверсия"
            print $3"\t"$1"\t"$2"\t"svtype"\t"svlen"\t"end_pos"\t"$7"\t"gene"\t"interp
        }'
} > "$SV_OUT"
echo "    Сохранено: ${SV_OUT}"

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> SV/CNV демо завершено за ${ELAPSED}с"
echo ""
echo "    Следующий шаг: bash scripts/40_modern_context.sh"
echo "    Документация:  docs/vcf_consequence_glossary.md"
