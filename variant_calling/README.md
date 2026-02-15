# Поиск геномных вариантов (Variant Calling)

**Продолжительность:** ~1.5-2 часа

Практикум по поиску и анализу геномных вариантов (SNP/InDel) у *Bacillus subtilis* 168 с использованием freebayes, bcftools и samtools на HPC.

## Чему вы научитесь

1. **Variant calling** — поиск SNP и indel с помощью freebayes из BAM-файла
2. **VCF format** — понимание формата VCF (variant call format) и метрик качества (QUAL, DP, AO/RO)
3. **Фильтрация вариантов** — отбор надёжных вариантов по порогам bcftools
4. **Оценка качества** — интерпретация метрик strand bias (SAP), mapping quality (MQM), Ti/Tv ratio

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

### Как интерпретировать результаты

- **raw.vcf.gz**: Все найденные варианты до фильтрации. Содержит метрики качества для каждого варианта.
- **filtered.vcf.gz**: Варианты, прошедшие пороги QUAL и DP. Рекомендуется для анализа.
- **stats.txt**: Итоговая статистика (SNP/indel counts, Ti/Tv ratio).

**Ключевые метрики VCF** (подробно в `docs/vcf_glossary.md`):
- **QUAL** — Phred-scaled quality score (>= 30 рекомендуется)
- **DP** — Depth (глубина покрытия; >= 20 хорошо)
- **AO/RO** — ALT / REF observations (для гаплоида AO/DP ~100%)
- **SAP** — Strand bias (< 20 хорошо; > 20 подозрительно)
- **MQM** — Mapping quality (>= 20 хорошо)

## Что дальше?

После успешного завершения рабочего процесса:

1. **Проверьте статистику:** `bash scripts/50_stats_bcftools.sh` покажет итоговые цифры
2. **Изучите примеры:** `bash scripts/60_pick_3_variants.sh` выберет 3 варианта для обсуждения
3. **Домашнее задание:** См. `docs/homework_02.md` для HW2 "Variant Defense" (ручная проверка 5 вариантов)

**Полезные команды:**
```bash
make all             # Запустить весь workflow
make clean           # Удалить outputs/ (сохраняет data/)
make homework-check  # Проверить формат ДЗ
```

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

- **SSH-доступ к HPC** (вычислительный узел, не login node) или Linux/macOS
- **Conda / Micromamba** для управления окружением
- **Ресурсы:** ~2 ГБ дискового пространства, 4-8 CPU cores (рекомендуется)
- **Данные:** Референс + BAM подготавливаются преподавателем (`01_get_data.sh`)

> **Для преподавателей:** См. `../shared/HPC_SETUP.md` для детальных инструкций по настройке HPC-окружения и `INSTRUCTOR.md` для подготовки данных.

## Организм

**Bacillus subtilis subsp. subtilis str. 168** (геном ~4.2 Мб). Используются симулированные прочтения (wgsim) с ~50x покрытием и ~3000-5000 внесёнными мутациями (SNP + indel). Это позволяет контролировать истинные варианты и оценивать precision/recall.

## Лицензия

MIT. См. [LICENSE](../LICENSE).
