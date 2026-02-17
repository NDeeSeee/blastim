#!/usr/bin/env bash
set -euo pipefail

# ── 20_prioritize_variants.sh ─────────────────────────────────────────
# Шаг 20: Приоритизация вариантов по клинической значимости.
#
# Алгоритм:
#   1. Читаем VEP TSV (из шага 10)
#   2. Фильтруем по MIN_DP и AF_RARE_THRESHOLD
#   3. Считаем приоритет: IMPACT (HIGH+3, MODERATE+2, LOW+1) +
#                         редкость (gnomAD_AF<0.001: +2, <0.01: +1) +
#                         ClinVar (Pathogenic: +3, Likely_pathogenic: +2)
#   4. Сохраняем приоритизированный TSV с рангом
#
# Usage:
#   bash scripts/20_prioritize_variants.sh
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

OUTDIR="${REPO_ROOT}/${OUTDIR}"
VEP_TSV="${OUTDIR}/vep/vep_output.tsv"

# Проверяем наличие VEP output
if [ ! -f "$VEP_TSV" ]; then
    echo "ERROR: VEP TSV не найден: ${VEP_TSV}"
    echo "       Сначала выполните: bash scripts/10_run_vep.sh"
    exit 1
fi

PRIORITIZED="${OUTDIR}/prioritized.tsv"

echo "==> Шаг 20: Приоритизация вариантов"
echo "    VEP TSV: ${VEP_TSV}"
echo "    AF_RARE_THRESHOLD: ${AF_RARE_THRESHOLD:-0.01}"
echo "    MIN_DP: ${MIN_DP:-10}"
echo ""

START=$(date +%s)

python3 - <<PYTHON
import sys
import re

VEP_TSV = "${VEP_TSV}"
OUT     = "${PRIORITIZED}"
AF_THR  = float("${AF_RARE_THRESHOLD:-0.01}")
MIN_DP  = int("${MIN_DP:-10}")

IMPACT_SCORE = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
CLINVAR_SCORE = {
    "Pathogenic": 3,
    "Likely_pathogenic": 2,
    "Pathogenic/Likely_pathogenic": 3,
}

def parse_extra(extra: str) -> dict:
    d = {}
    for token in extra.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            d[k.strip()] = v.strip()
    return d

rows = []
with open(VEP_TSV) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 14:
            continue
        (uploaded, location, allele, gene, feature, ftype, consequence,
         cdna, cds, prot, aa, codons, existing, extra_str) = parts[:14]

        extra = parse_extra(extra_str)
        impact = extra.get("IMPACT", "MODIFIER")
        symbol = extra.get("SYMBOL", gene)
        canonical = extra.get("CANONICAL", "-")
        hgvsc = extra.get("HGVSc", "-")
        hgvsp = extra.get("HGVSp", "-")
        gnomad_af_str = extra.get("gnomADe_AF", "-")
        clinvar = extra.get("ClinVar_CLNSIG", "-")

        # Только canonical транскрипты для сводной таблицы
        if canonical not in ("YES", "1"):
            continue

        # gnomAD AF
        try:
            gnomad_af = float(gnomad_af_str)
        except ValueError:
            gnomad_af = None

        # Scoring
        score = IMPACT_SCORE.get(impact, 0)
        if gnomad_af is not None:
            if gnomad_af < 0.001:
                score += 2
            elif gnomad_af < AF_THR:
                score += 1
        for key, pts in CLINVAR_SCORE.items():
            if key.lower() in clinvar.lower():
                score += pts
                break

        rows.append({
            "score": score,
            "variant": uploaded,
            "location": location,
            "gene": symbol,
            "consequence": consequence,
            "impact": impact,
            "hgvsc": hgvsc,
            "hgvsp": hgvsp,
            "gnomad_af": gnomad_af_str,
            "clinvar": clinvar,
            "existing": existing,
        })

# Сортируем по убыванию score, затем по IMPACT
IMPACT_ORDER = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
rows.sort(key=lambda r: (-r["score"],
                          IMPACT_ORDER.index(r["impact"])
                          if r["impact"] in IMPACT_ORDER else 99))

with open(OUT, "w") as f:
    header = (
        "rank\tscore\tvariant\tlocation\tgene\tconsequence\timpact\t"
        "hgvsc\thgvsp\tgnomad_af\tclinvar\texisting_id\n"
    )
    f.write(header)
    for rank, r in enumerate(rows, 1):
        f.write(
            f"{rank}\t{r['score']}\t{r['variant']}\t{r['location']}\t"
            f"{r['gene']}\t{r['consequence']}\t{r['impact']}\t"
            f"{r['hgvsc']}\t{r['hgvsp']}\t{r['gnomad_af']}\t"
            f"{r['clinvar']}\t{r['existing']}\n"
        )

print(f"    Записано {len(rows)} вариантов → {OUT}")
if rows:
    print(f"    Топ-3 варианта:")
    for r in rows[:3]:
        print(f"      #{r['score']:2d} {r['gene']:10s}  {r['impact']:8s}  {r['consequence']}")
PYTHON

END=$(date +%s)
ELAPSED=$(( END - START ))

echo ""
echo "==> Приоритизация завершена за ${ELAPSED}с"
echo "    Файл: ${PRIORITIZED}"
echo ""
echo "    Следующий шаг: python3 scripts/21_make_report.py (или make prioritize)"
