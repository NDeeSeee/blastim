# Практикум 3: Интерпретация геномных вариантов

**Лекция 3** курса по биоинформатике NGS | Организм: *Homo sapiens* GRCh38

## Цели

По завершению практикума студент умеет:

1. Запустить VEP-аннотацию SNV/indel (REST и offline режимы)
2. Интерпретировать поля VEP TSV: IMPACT, HGVSc/p, gnomAD_AF, ClinVar
3. Приоритизировать варианты по клинической значимости (scoring)
4. Читать SV VCF: SVTYPE, SVLEN, END, символьные ALT
5. Понимать иерархию клинических баз данных (ClinVar, gnomAD, OMIM)
6. Применять базовые ACMG/AMP критерии для вердикта по варианту

---

## Быстрый старт

```bash
# 1. Создать conda-окружение (один раз)
make env

# 2. Скопировать конфигурацию
cp workflow/config.sh.example workflow/config.sh

# 3. Проверить систему
make check

# 4. Запустить полный пайплайн (шаги 10–30)
make all

# 5. Просмотреть результаты
ls outputs/run_*/
cat outputs/run_*/report.md
```

> **Примечание:** VEP работает в двух режимах. По умолчанию используется REST API (нужен интернет). Если интернет недоступен, скрипт автоматически переключится на предобработанный файл `data/input/demo_variants_vep.tsv`.

---

## Пайплайн

```
data/input/demo_variants.vcf   ──►  [10] VEP аннотация  ──►  vep_output.tsv
                                                                      │
                                     [20] Приоритизация  ◄────────────┘
                                             │
                                     prioritized.tsv
                                             │
                                     [21] Отчёт  ──►  report.md

data/input/demo_sv.vcf         ──►  [30] SV/CNV демо  ──►  sv_annotated.tsv
```

| Шаг | Скрипт | Что делает |
|-----|--------|-----------|
| 10 | `scripts/10_run_vep.sh` | VEP аннотация вариантов |
| 20 | `scripts/20_prioritize_variants.sh` | Приоритизация по IMPACT/AF/ClinVar |
| 21 | `scripts/21_make_report.py` | Генерация markdown-отчёта |
| 30 | `scripts/30_sv_cnv_intro.sh` | Разбор структурных вариантов |
| 40 | `scripts/40_modern_context.sh` | Обзор современных инструментов (необяз.) |

---

## Демо-данные

### SNV/indel (`data/input/demo_variants.vcf`)

20 вариантов в ключевых онкогенах (GRCh38):

| Ген | Вариант | rsID | Significance |
|-----|---------|------|-------------|
| TP53 | R175H | rs28934574 | Pathogenic |
| TP53 | R248W | rs28934578 | Pathogenic |
| TP53 | P72R | rs1042522 | Benign |
| BRAF | V600E | rs113488022 | Pathogenic |
| EGFR | T790M | rs121913459 | Pathogenic |
| EGFR | L858R | rs121913375 | Pathogenic |
| BRCA1 | 5382insC | rs80357906 | Pathogenic |
| PIK3CA | E545K | rs104886003 | Pathogenic |
| KRAS | G12D | rs121913529 | Pathogenic |
| … | … | … | … |

### SV (`data/input/demo_sv.vcf`)

5 структурных вариантов: DEL TP53, DUP BRCA1, INV BRAF, DUP EGFR, DEL chr1.

---

## Ожидаемые результаты

```
outputs/run_YYYYMMDD_HHMMSS/
├── vep/
│   └── vep_output.tsv         # 20 строк данных VEP
├── sv/
│   └── sv_annotated.tsv       # 5 SV с интерпретацией
├── prioritized.tsv            # Варианты с приоритетным скором
└── report.md                  # Итоговый markdown-отчёт
```

Подробнее: `docs/expected_outputs.md`

---

## Структура модуля

```
variant_interpretation/
├── data/
│   ├── input/
│   │   ├── demo_variants.vcf         # 20 SNV/indel (GRCh38)
│   │   ├── demo_sv.vcf               # 5 структурных вариантов
│   │   └── demo_variants_vep.tsv     # Предобработанная VEP аннотация
│   └── db/                           # Дополнительные базы (опционально)
├── env/
│   ├── environment.yml               # ngs-interpret conda environment
│   └── micromamba_create.sh
├── scripts/
│   ├── 00_check_system.sh            # Проверка зависимостей
│   ├── 01_get_data.sh                # Проверка входных данных
│   ├── 02_get_vep_cache.sh           # Загрузка VEP cache (преподаватель)
│   ├── 10_run_vep.sh                 # VEP аннотация (REST или offline)
│   ├── 20_prioritize_variants.sh     # Приоритизация вариантов
│   ├── 21_make_report.py             # Генерация отчёта
│   ├── 30_sv_cnv_intro.sh            # SV/CNV демо
│   ├── 40_modern_context.sh          # Современные инструменты
│   ├── 90_assign_variants.py         # Назначение вариантов студентам
│   ├── 91_homework_check.py          # Проверка ДЗ
│   ├── _vep_rest_helper.py           # VEP REST API (вызывается из 10)
│   └── _cnv_plot.py                  # ASCII SV визуализация
├── workflow/
│   ├── config.sh.example             # Шаблон конфигурации
│   └── run_all.sh                    # Оркестратор полного пайплайна
├── docs/
│   ├── cheatsheet.md                 # Шпаргалка bcftools/VEP
│   ├── vcf_consequence_glossary.md   # Глоссарий VEP consequences
│   ├── expected_outputs.md           # Ожидаемые результаты
│   ├── troubleshooting.md            # Устранение проблем
│   ├── advanced_reading.md           # Рекомендуемые статьи
│   └── homework_03.md                # Описание ДЗ3
├── homework/
│   ├── TEMPLATE_hw3.md               # Шаблон для защиты вариантов
│   └── submissions/                  # Работы студентов (PR)
├── outputs/                          # Создаётся при запуске
├── Makefile
├── CONTRIBUTING.md
└── INSTRUCTOR.md
```

---

## Окружение

```yaml
name: ngs-interpret
dependencies:
  - python=3.12
  - ensembl-vep=110.0
  - bcftools=1.20
  - samtools=1.20
  - tabix
  - pandas=2.2
  - matplotlib=3.9
```

---

## Домашнее задание (ДЗ3)

**Формат:** "Защита вариантов" — интерпретация 8 персонально назначенных вариантов.

```bash
export GITHUB_USER=<ваш_github>

# Получить список вариантов
python3 scripts/90_assign_variants.py \
    --vep_tsv data/input/demo_variants_vep.tsv \
    --sv_vcf  data/input/demo_sv.vcf \
    --out     homework/submissions/${GITHUB_USER}/variants_${GITHUB_USER}.txt

# Заполнить шаблон
cp homework/TEMPLATE_hw3.md \
   homework/submissions/${GITHUB_USER}/hw3_${GITHUB_USER}.md

# Проверить
make homework-check GITHUB_USER=${GITHUB_USER}
```

Подробнее: `docs/homework_03.md`

---

## Документация

| Файл | Описание |
|------|----------|
| `docs/cheatsheet.md` | bcftools, VEP, python — команды для быстрой работы |
| `docs/vcf_consequence_glossary.md` | Все VEP consequences с описаниями |
| `docs/expected_outputs.md` | Что должен содержать каждый выходной файл |
| `docs/troubleshooting.md` | Решения типичных проблем |
| `docs/advanced_reading.md` | Ключевые статьи и онлайн-ресурсы |
| `docs/homework_03.md` | Условие и критерии ДЗ3 |
| `INSTRUCTOR.md` | Подготовка окружения, план занятия, Q&A |
