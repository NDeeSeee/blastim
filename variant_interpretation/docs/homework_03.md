# ДЗ3: Интерпретация геномных вариантов

## Цель

Освоить практическую интерпретацию вариантов: от VEP-аннотации до вердикта по ACMG-критериям.

## Сроки

| Событие | Дата |
|---------|------|
| Выдача задания | После Лекции 3 |
| Дедлайн сдачи | +7 дней |
| Разбор результатов | +10 дней (Лекция 4) |

---

## Задание

### Шаг 1: Получите свой список вариантов

```bash
# Из директории variant_interpretation/
export GITHUB_USER=<ваш_github>
python3 scripts/90_assign_variants.py \
    --vep_tsv data/input/demo_variants_vep.tsv \
    --sv_vcf  data/input/demo_sv.vcf \
    --out     homework/submissions/${GITHUB_USER}/variants_${GITHUB_USER}.txt
```

Вам будут назначены **8 вариантов** с детерминированным выбором на основе вашего username:
- 2× HIGH-impact (патогенные / вероятно патогенные)
- 2× MODERATE-impact (VUS — missense)
- 1× LOW-impact (синонимичный или доброкачественный)
- 1× MODIFIER (интронный / UTR)
- 1× indel / frameshift
- 1× структурный вариант (SV)

### Шаг 2: Заполните шаблон

```bash
mkdir -p homework/submissions/${GITHUB_USER}
cp homework/TEMPLATE_hw3.md \
   homework/submissions/${GITHUB_USER}/hw3_${GITHUB_USER}.md
```

Для каждого из 8 вариантов заполните:

1. **Геномную позицию** — хромосома, координата GRCh38, REF/ALT
2. **Аннотацию VEP** — ген, consequence, IMPACT, HGVSc/p, gnomAD_AF, ClinVar
3. **Качество секвенирования** — QUAL, DP, AF, FILTER из VCF
4. **Приоритетный скор** — по алгоритму из шага 20
5. **Вердикт** — один из 5 классов ACMG + обоснование (≥3 предложения)
6. **Предлагаемый FILTER** — bcftools выражение для фильтрации этого варианта

### Шаг 3: Самопроверка

```bash
make homework-check GITHUB_USER=${GITHUB_USER}
```

Проверяет:
- Наличие файлов `variants_<user>.txt` и `hw3_<user>.md`
- Обязательные ключевые слова: QUAL, DP, VEP, IMPACT, gnomAD, SVTYPE
- Незаполненные плейсхолдеры `<...>`

### Шаг 4: Создайте PR

```bash
git checkout -b hw3-${GITHUB_USER}
git add homework/submissions/${GITHUB_USER}/
git commit -m "HW3: variant defense [${GITHUB_USER}]"
git push origin hw3-${GITHUB_USER}
```

Откройте Pull Request на GitHub с заголовком:
`HW3: Variant Defense — <GITHUB_USER>`

---

## Критерии оценки

| Критерий | Баллы |
|----------|-------|
| Геномные позиции верно заполнены (8/8) | 2.0 |
| VEP-аннотация: ген, consequence, IMPACT (8/8) | 2.0 |
| HGVSc/HGVSp корректны (6/8 минимум) | 1.0 |
| gnomAD_AF и ClinVar указаны (8/8) | 1.0 |
| Вердикт ACMG выбран с обоснованием ≥3 предложений | 2.0 |
| Предлагаемый FILTER: корректное bcftools выражение | 1.0 |
| SV-вариант: тип, размер, гены, клиническое значение | 1.0 |
| **Итого** | **10.0** |
| Рефлексия (необязательно) | +0.5 |

---

## Подсказки

### Где найти данные для каждого варианта

```bash
# VCF с качеством (QUAL, DP, AF, FILTER):
cat data/input/demo_variants.vcf | grep -v "^#"

# VEP-аннотация (consequence, IMPACT, HGVSc/p, gnomAD_AF):
cat data/input/demo_variants_vep.tsv | grep -v "^##"

# SV данные:
cat data/input/demo_sv.vcf | grep -v "^#"
```

### Алгоритм расчёта скора

```
score = 0
if IMPACT == "HIGH":     score += 3
elif IMPACT == "MODERATE": score += 2
elif IMPACT == "LOW":    score += 1

if gnomAD_AF < 0.001:   score += 2
elif gnomAD_AF < 0.01:  score += 1

if ClinVar == "Pathogenic": score += 3
elif ClinVar == "Likely_pathogenic": score += 2
```

### ACMG критерии быстрая шпаргалка

| Класс | Когда выбирать |
|-------|---------------|
| **Патогенный** | HIGH impact + ClinVar Pathogenic + редкий (gnomAD < 0.001%) |
| **Вероятно патогенный** | HIGH/MODERATE + косвенные доказательства патогенности |
| **VUS** | MODERATE impact + Uncertain_significance в ClinVar |
| **Вероятно доброкачественный** | Умеренная частота в gnomAD (0.01–0.05) + LOW impact |
| **Доброкачественный** | BA1: gnomAD AF > 5%; BS3: функциональные данные норма |

---

## Типичные ошибки

1. **Не различать IMPACT и consequence** — IMPACT это категория (HIGH/MODERATE), consequence это конкретный тип (missense_variant, frameshift_variant)
2. **Игнорировать DP < 10** — низкое покрытие снижает доверие к варианту
3. **Путать ClinVar и ACMG вердикт** — ClinVar = база данных, ACMG = ваша классификация
4. **SV без указания размера** — SVLEN обязателен для интерпретации
5. **Плейсхолдеры `<...>` не заменены** — проверка сообщит об ошибке

---

## Дополнительные ресурсы

- `docs/vcf_consequence_glossary.md` — список всех VEP consequences
- `docs/cheatsheet.md` — bcftools и VCF команды
- `docs/advanced_reading.md` — статьи по теме
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) — поиск по rs ID
- [gnomAD browser](https://gnomad.broadinstitute.org/) — частоты популяций
- [ACMG guidelines 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)
