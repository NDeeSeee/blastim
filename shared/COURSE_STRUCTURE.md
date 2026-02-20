# Структура курса: NGS-анализ на HPC

Этот курс состоит из трёх независимых практикумов, которые преподаются студентам программы по Next-Generation Sequencing.

## Обзор практикумов

| Практикум | Организм | Объём | Ключевые темы | Директория |
|-----------|---------|-------|---------------|------------|
| **1. Аннотация прокариотического генома** | *B. subtilis* 168 | 2 часа | Структурная аннотация (Prokka), функциональная аннотация (DIAMOND), поиск доменов (HMMER/Pfam), "лестница доказательств" | `genome_annotation/` |
| **2. Поиск геномных вариантов** | *B. subtilis* 168 | 2 часа | Выравнивание (BWA), поиск вариантов (freebayes), фильтрация (bcftools), оценка качества вариантов | `variant_calling/` |
| **3. Интерпретация геномных вариантов** | *H. sapiens* GRCh38 | 2 часа | VEP-аннотация, приоритизация (IMPACT/gnomAD/ClinVar), SV/CNV, ACMG критерии | `variant_interpretation/` |

## Организмы

### Практикумы 1–2: *Bacillus subtilis* subsp. subtilis str. 168 (GenoBase)
- Размер генома: ~4.2 Мб
- GC-содержание: ~43%
- Число предсказанных белков: ~4200
- Данные загружаются автоматически при первом запуске (`01_get_data.sh`)

### Практикум 3: *Homo sapiens* (GRCh38 / hg38)
- Демо-данные: 20 вариантов в TP53, BRCA1, BRAF, EGFR
- 5 структурных вариантов (DEL, DUP, INV)
- Данные включены в репозиторий (`data/input/`), загрузки не требуется

## Целевая аудитория

- **genome_annotation**: Структурная и функциональная аннотация прокариотических геномов
- **variant_calling**: Поиск вариантов в популяционной и медицинской геномике
- **variant_interpretation**: Клиническая интерпретация вариантов, базы данных ClinVar/gnomAD, ACMG критерии

## Структура репозитория

```
blastim_ngs/
├── README.md                    ← Вы здесь
├── LICENSE
├── .gitignore
├── shared/                      ← Общие ресурсы (этот файл)
│   ├── README.md
│   ├── COURSE_STRUCTURE.md
│   ├── HPC_SETUP.md
│   └── TROUBLESHOOTING_COMMON.md
├── genome_annotation/           ← Практикум 1: Аннотация генома
│   ├── README.md
│   ├── INSTRUCTOR.md
│   ├── CONTRIBUTING.md
│   ├── Makefile
│   ├── env/                     # ngs-annotation conda env
│   ├── scripts/                 # 00_check_system.sh … 91_homework_check.py
│   ├── workflow/                # config.sh.example, run_all.sh
│   ├── data/                    # input/, db/ (загружаются)
│   ├── homework/                # TEMPLATE_summary.md, submissions/
│   └── docs/                    # cheatsheet, glossary, hw01, troubleshooting
├── variant_calling/             ← Практикум 2: Variant Calling
│   ├── README.md
│   ├── INSTRUCTOR.md
│   ├── CONTRIBUTING.md
│   ├── Makefile
│   ├── env/                     # ngs-variants conda env
│   ├── scripts/                 # 00_check_system.sh … 91_homework_check.py
│   ├── workflow/                # config.sh.example, run_all.sh
│   ├── data/                    # input/, truth/
│   ├── homework/                # TEMPLATE_variant_defense.md, submissions/
│   └── docs/                    # cheatsheet, glossary, hw02, troubleshooting
└── variant_interpretation/      ← Практикум 3: Интерпретация вариантов
    ├── README.md
    ├── INSTRUCTOR.md
    ├── CONTRIBUTING.md
    ├── Makefile
    ├── env/                     # ngs-interpret conda env
    ├── scripts/                 # 00_check_system.sh … 91_homework_check.py
    ├── workflow/                # config.sh.example, run_all.sh
    ├── data/                    # input/ (включены в репо), db/
    ├── homework/                # TEMPLATE_hw3.md, submissions/
    └── docs/                    # cheatsheet, glossary, hw03, troubleshooting
```

## Как начать

### Студентам

1. Выберите интересующий вас практикум:
   ```bash
   cd genome_annotation       # Практикум 1
   cd variant_calling         # Практикум 2
   cd variant_interpretation  # Практикум 3
   ```

2. Скопируйте конфигурацию:
   ```bash
   cp workflow/config.sh.example workflow/config.sh
   ```

3. Проверьте систему и следуйте **README.md** выбранного практикума:
   ```bash
   make check
   make all
   ```

### Преподавателям

1. Прочитайте **shared/HPC_SETUP.md** для подготовки HPC-системы
2. Для каждого практикума см. **INSTRUCTOR.md** в соответствующей папке
3. Консультируйтесь с **shared/TROUBLESHOOTING_COMMON.md** при проблемах

## Домашние задания

- **ДЗ1 (genome_annotation)**: Защита 3 генов с использованием "лестницы доказательств"
  - Условие: `genome_annotation/docs/homework_01.md`
  - Шаблон: `genome_annotation/homework/TEMPLATE_summary.md`
  - Сдача: fork/PR, см. `genome_annotation/CONTRIBUTING.md`

- **ДЗ2 (variant_calling)**: Защита 5 вариантов с анализом качества
  - Условие: `variant_calling/docs/homework_02.md`
  - Шаблон: `variant_calling/homework/TEMPLATE_variant_defense.md`
  - Сдача: fork/PR, см. `variant_calling/CONTRIBUTING.md`

- **ДЗ3 (variant_interpretation)**: Защита 8 вариантов с ACMG-вердиктами
  - Условие: `variant_interpretation/docs/homework_03.md`
  - Шаблон: `variant_interpretation/homework/TEMPLATE_hw3.md`
  - Сдача: fork/PR, см. `variant_interpretation/CONTRIBUTING.md`

## Требования к системе

| Компонент | Практикумы 1–2 | Практикум 3 |
|-----------|----------------|-------------|
| CPU | 4+ cores | 4+ cores |
| RAM | 8 ГБ | 8 ГБ |
| Диск | ~5 ГБ/модуль | ~1 ГБ (без VEP cache) / ~16 ГБ (с cache) |
| Сеть | Нужен (загрузка данных) | Опционально (VEP REST или precomputed) |
| Conda env | `ngs-annotation`, `ngs-variants` | `ngs-interpret` |

Детали см. в **shared/HPC_SETUP.md**.

## Лицензия

MIT. См. [LICENSE](../LICENSE).

Данные и базы (Swiss-Prot, Pfam, B. subtilis 168, gnomAD, ClinVar) загружаются из официальных источников и подчиняются их собственным лицензиям.
