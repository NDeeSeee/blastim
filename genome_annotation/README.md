# Аннотация прокариотического генома на HPC

Двухчасовой практикум для курсов NGS. Студенты учатся аннотировать бактериальный геном с помощью стандартных биоинформатических инструментов на HPC-системе.

## Чему вы научитесь

1. **Структурная аннотация** -- предсказание белок-кодирующих генов (CDS), rRNA, tRNA с помощью Prokka. Понимание форматов GFF/FAA/TSV.
2. **Функциональная аннотация по гомологии** -- поиск предсказанных белков в референсной базе данных с помощью DIAMOND blastp.
3. **Функциональная аннотация по доменам** -- определение консервативных белковых доменов с помощью hmmscan (HMMER) по базе Pfam.
4. **Лестница доказательств** -- объединение предсказания ORF, гомологии и доменной аннотации для оценки достоверности.

## Предварительные требования

- SSH-доступ к HPC
- Базовые навыки командной строки Linux (`cd`, `ls`, `less`, `head`)
- Conda-окружение `ngs-annotation` (его настроит преподаватель)

## Быстрый старт

```bash
# 1. Клонируйте репозиторий и перейдите в практикум
git clone git@github.com:NDeeSeee/blastim_ngs.git
cd blastim_ngs/genome_annotation

# 2. Создайте и активируйте окружение
bash env/micromamba_create.sh
micromamba activate ngs-annotation   # или: conda activate ngs-annotation

# 3. Скопируйте и проверьте конфигурацию
cp workflow/config.sh.example workflow/config.sh
# Отредактируйте workflow/config.sh при необходимости (например, измените CPUs)

# 4. Проверьте настройку
bash scripts/00_check_system.sh

# 5. Запустите полный рабочий процесс
bash workflow/run_all.sh
```

Или пошагово:

```bash
bash scripts/10_run_prokka.sh       # Структурная аннотация
bash scripts/20_run_diamond.sh      # Поиск по гомологии
bash scripts/30_run_hmmscan.sh      # Поиск доменов
bash scripts/40_sanity_checks.sh    # Проверка качества
bash scripts/50_pick_3_genes.sh     # Обсуждение примеров генов
```

## Понимание результатов

После завершения рабочего процесса в директории результатов (`outputs/run_YYYYMMDD_HHMMSS/`) находятся:

| Директория | Ключевые файлы | Описание |
|-----------|-----------|---------------|
| `prokka/` | `ANNOT.gff` | Аннотированный геном в формате GFF3 |
| | `ANNOT.faa` | Предсказанные белковые последовательности (FASTA) |
| | `ANNOT.ffn` | Нуклеотидные последовательности CDS |
| | `ANNOT.tsv` | Таблица признаков (tab-separated) |
| | `ANNOT.txt` | Сводная статистика |
| `diamond/` | `blastp_results.tsv` | Результаты DIAMOND blastp (табличный формат BLAST) |
| `hmmer/` | `hmmscan_domtblout.txt` | Найденные домены Pfam для каждого белка |

### Как интерпретировать

- **GFF-файл**: Каждая строка -- геномный признак (feature). Смотрите колонку 3 (тип) и колонку 9 (атрибуты, включая название продукта).
- **DIAMOND TSV**: Стандартный табличный формат BLAST. Ключевые колонки: ID запроса, ID субъекта, % идентичности, E-value, bit score, описание субъекта.
- **domtblout**: Каждая строка -- найденный домен. Ключевые колонки: имя домена, accession, белок-запрос, E-value.

## Структура практикума

```
genome_annotation/
├── README.md              ← Вы здесь
├── INSTRUCTOR.md          ← Руководство преподавателя
├── Makefile               ← make env / make all / make clean
├── env/
│   ├── environment.yml    ← Спецификация conda-окружения
│   └── micromamba_create.sh
├── scripts/
│   ├── 00_check_system.sh ← Проверка окружения HPC
│   ├── 01_get_data.sh     ← Загрузка генома + белков
│   ├── 02_make_diamond_db.sh
│   ├── 03_get_pfam.sh     ← Загрузка и подмножество Pfam
│   ├── 10_run_prokka.sh   ← Структурная аннотация
│   ├── 20_run_diamond.sh  ← Поиск по гомологии
│   ├── 30_run_hmmscan.sh  ← Поиск доменов
│   ├── 40_sanity_checks.sh
│   └── 50_pick_3_genes.sh ← Примеры для обсуждения
├── workflow/
│   ├── config.sh.example  ← Скопируйте в config.sh
│   └── run_all.sh         ← Запуск всего одной командой
├── data/
│   ├── input/             ← FASTA-файл сборки
│   └── db/                ← Базы данных DIAMOND + Pfam
├── outputs/               ← Создаётся при запуске
├── homework/
│   ├── TEMPLATE_summary.md  ← Шаблон для Части 1
│   ├── TEMPLATE_case.md     ← Шаблон для Части 2
│   └── submissions/         ← Папки студентов (через fork/PR)
└── docs/
    ├── cheatsheet.md        ← Все команды для копирования
    ├── expected_outputs.md
    ├── troubleshooting.md
    ├── advanced_reading.md  ← Углублённое чтение (статьи, туториалы)
    └── homework_01.md       ← Домашнее задание 1
```

## Организм

В этом курсе используется **Bacillus subtilis subsp. subtilis str. 168** (геном ~4.2 Мб, ~4200 предсказанных белков). Это хорошо изученный модельный организм, поэтому большинство предсказаний будут иметь подтверждение гомологией и доменами -- идеально для обучения.

## Лицензия

MIT. См. [LICENSE](../LICENSE).

Базы данных (Swiss-Prot, Pfam) загружаются из официальных источников и подчиняются собственным лицензиям. Они не включены в этот репозиторий.
