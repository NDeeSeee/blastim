# Поиск геномных вариантов (Variant Calling)

Практикум по поиску и анализу геномных вариантов (SNP/InDel) у *Bacillus subtilis* 168 с использованием freebayes, bcftools и samtools на HPC.

## Быстрый старт

```bash
# 1. Создайте окружение
bash env/micromamba_create.sh
micromamba activate ngs-variants

# 2. Скопируйте конфигурацию
cp workflow/config.sh.example workflow/config.sh

# 3. Проверьте систему
bash scripts/00_check_system.sh

# 4. Запустите полный рабочий процесс
bash workflow/run_all.sh

# Или пошагово
bash scripts/10_prep_reference.sh
bash scripts/11_prep_bam.sh
bash scripts/20_call_freebayes.sh
bash scripts/30_compress_index_vcf.sh
bash scripts/40_filter_bcftools.sh
bash scripts/50_stats_bcftools.sh
bash scripts/60_pick_3_variants.sh
```

## Что делает каждый шаг

| Шаг | Скрипт | Что делает |
|-----|--------|------------|
| 10 | `10_prep_reference.sh` | Индексирование референсного генома (`samtools faidx`) |
| 11 | `11_prep_bam.sh` | Валидация BAM-файла (`flagstat`, `idxstats`, `depth`) |
| 20 | `20_call_freebayes.sh` | Поиск вариантов с freebayes → `raw.vcf` |
| 30 | `30_compress_index_vcf.sh` | Сжатие bgzip + индекс tabix |
| 40 | `40_filter_bcftools.sh` | Фильтрация по QUAL, DP → `filtered.vcf.gz` |
| 50 | `50_stats_bcftools.sh` | Статистика: SNP/indel, Ti/Tv |
| 60 | `60_pick_3_variants.sh` | Выбор 3 примеров для обсуждения |

## Ожидаемые результаты (приблизительно)

| Метрика | Значение |
|---------|----------|
| Всего вариантов (raw) | 3 000 - 5 000 |
| SNP | 2 500 - 4 000 |
| Indel | 300 - 800 |
| Ti/Tv | 2.0 - 2.5 |
| Проходят фильтры | 70% - 85% |

Подробнее: `docs/expected_outputs.md`

## Структура практикума

```
variant_calling/
├── README.md              ← Вы здесь
├── INSTRUCTOR.md          ← Руководство преподавателя
├── CONTRIBUTING.md        ← Как сдать ДЗ (fork/PR)
├── Makefile               ← make all / make clean / ...
├── env/
│   ├── environment.yml    ← Conda-окружение ngs-variants
│   └── micromamba_create.sh
├── workflow/
│   ├── config.sh.example  ← Шаблон конфигурации
│   └── run_all.sh         ← Запуск всего пайплайна
├── scripts/
│   ├── 00_check_system.sh
│   ├── 01_get_data.sh     ← Подготовка данных (преподаватель)
│   ├── 10_prep_reference.sh
│   ├── 11_prep_bam.sh
│   ├── 20_call_freebayes.sh
│   ├── 30_compress_index_vcf.sh
│   ├── 40_filter_bcftools.sh
│   ├── 50_stats_bcftools.sh
│   ├── 60_pick_3_variants.sh
│   ├── 90_homework_assign_variants.py
│   └── 91_homework_check.py
├── data/
│   ├── input/             ← Референс + BAM (создаются 01_get_data.sh)
│   └── truth/             ← Истинные мутации wgsim
├── homework/
│   ├── TEMPLATE_variant_defense.md
│   └── submissions/
└── docs/
    ├── glossary.md        ← Термины, концепции, инструменты
    ├── vcf_glossary.md    ← Справочник полей VCF
    ├── cheatsheet.md      ← Шпаргалка по командам
    ├── homework_02.md     ← Задание ДЗ2
    ├── expected_outputs.md
    ├── troubleshooting.md
    └── advanced_reading.md
```

## Предварительные требования

- SSH-доступ к HPC (или Linux/macOS с conda)
- Conda / Micromamba
- ~2 ГБ свободного места на диске
- Подготовленные данные (преподаватель запускает `01_get_data.sh` заранее)

## Домашнее задание

**HW2: "Variant Defense"** — оценка 5 персональных вариантов.

Подробности: `docs/homework_02.md`
