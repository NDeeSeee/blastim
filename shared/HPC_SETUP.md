# Подготовка HPC-системы для курса

Руководство для преподавателя: как подготовить HPC (или локальную систему) для проведения практикумов.

## Быстрая подготовка (20–30 минут)

```bash
# 1. Клонируйте репозиторий
git clone git@github.com:NDeeSeee/blastim_ngs.git
cd blastim_ngs

# 2. Создайте все три conda-окружения
bash genome_annotation/env/micromamba_create.sh     # ngs-annotation
bash variant_calling/env/micromamba_create.sh        # ngs-variants
bash variant_interpretation/env/micromamba_create.sh # ngs-interpret

# 3. Подготовьте данные (Практикумы 1–2 загружают данные; Практикум 3 — нет)
micromamba activate ngs-annotation
bash genome_annotation/scripts/01_get_data.sh
bash genome_annotation/scripts/02_make_diamond_db.sh
bash genome_annotation/scripts/03_get_pfam.sh

micromamba activate ngs-variants
bash variant_calling/scripts/01_get_data.sh

micromamba activate ngs-interpret
bash variant_interpretation/scripts/01_get_data.sh  # только проверка, всё уже в репо

# 4. (Опционально) Загрузите VEP cache для Практикума 3 в offline-режиме
# ВНИМАНИЕ: ~15 ГБ, 30–60 минут!
# cp variant_interpretation/workflow/config.sh.example variant_interpretation/workflow/config.sh
# Установите VEP_MODE=offline и VEP_CACHE_DIR=/shared/vep_cache
# bash variant_interpretation/scripts/02_get_vep_cache.sh

# 5. Проверьте, что всё работает
micromamba activate ngs-annotation  && bash genome_annotation/scripts/00_check_system.sh
micromamba activate ngs-variants    && bash variant_calling/scripts/00_check_system.sh
micromamba activate ngs-interpret   && bash variant_interpretation/scripts/00_check_system.sh
```

## Требования к системе

### Минимальные требования

| Параметр | Практикумы 1–2 | Практикум 3 |
|----------|----------------|-------------|
| **Процессоры** | 4+ cores | 4+ cores |
| **Оперативная память** | 8 ГБ | 8 ГБ |
| **Диск** | ~5 ГБ / модуль | ~1 ГБ (без VEP cache) / ~16 ГБ (с cache) |
| **ОС** | Linux / macOS / WSL2 | Linux / macOS / WSL2 |
| **Сеть** | Нужен (загрузка данных) | Опционально (VEP REST или precomputed) |

### Программное обеспечение

- **conda** или **micromamba** (рекомендуется micromamba)
- **git** (для клонирования репозитория)
- **SSH** (если используете удалённый HPC)

## Подробная подготовка

### 1. Установка Micromamba (рекомендуется)

```bash
# Linux / macOS
curl -L https://micro.mamba.pm/install.sh | bash
# Или скачайте с https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

# После установки добавьте в ~/.bashrc или ~/.zshrc:
export PATH="${MAMBA_ROOT_PREFIX}/bin:$PATH"
```

### 2. Клонирование и структура

```bash
git clone git@github.com:NDeeSeee/blastim_ngs.git
cd blastim_ngs
ls -la
# Вы должны увидеть: genome_annotation/, variant_calling/, shared/, README.md, LICENSE
```

### 3. Установка окружений

#### Вариант A: Micromamba (рекомендуется)

```bash
# Практикум 1
bash genome_annotation/env/micromamba_create.sh
# → окружение ngs-annotation

# Практикум 2
bash variant_calling/env/micromamba_create.sh
# → окружение ngs-variants

# Практикум 3
bash variant_interpretation/env/micromamba_create.sh
# → окружение ngs-interpret
```

#### Вариант B: Conda/Mamba

```bash
cd genome_annotation
conda env create -f env/environment.yml
conda activate ngs-annotation

cd ../variant_calling
conda env create -f env/environment.yml
conda activate ngs-variants

cd ../variant_interpretation
conda env create -f env/environment.yml
conda activate ngs-interpret
```

### 4. Подготовка данных

#### Для genome_annotation (этот шаг нужен один раз)

```bash
micromamba activate ngs-annotation
cd genome_annotation

# Загрузка референса B. subtilis 168 и подготовка баз
bash scripts/01_get_data.sh        # ~2 минуты, ~200 МБ

# Создание DIAMOND индекса
bash scripts/02_make_diamond_db.sh # ~1 минута

# Загрузка Pfam (подмножество)
bash scripts/03_get_pfam.sh        # ~1 минута

# Проверка
ls -lh data/input/
ls -lh data/db/
```

**Где это сохраняется:**
- `data/input/assembly.fasta` — B. subtilis 168 геном
- `data/db/swissprot.dmnd` — DIAMOND индекс Swiss-Prot
- `data/db/pfam_subset/` — Pfam подмножество (~200 семейств)

#### Для variant_calling (этот шаг нужен один раз)

```bash
micromamba activate ngs-variants
cd variant_calling

# Загрузка референса + симуляция прочтений + выравнивание (slow!)
bash scripts/01_get_data.sh        # ~10-15 минут

# Проверка
ls -lh data/input/
ls -lh data/truth/
```

**Где это сохраняется:**
- `data/input/reference.fasta` — B. subtilis 168 геном
- `data/input/reads_*.fastq` — Симулированные прочтения (50x coverage)
- `data/input/reads_sorted.bam` — Выровненные прочтения
- `data/truth/variants_truth.txt` — Истинные варианты (из wgsim)

#### Для variant_interpretation (данные уже в репозитории)

```bash
micromamba activate ngs-interpret
cd variant_interpretation
cp workflow/config.sh.example workflow/config.sh

# Данные уже включены в репо — только проверка
bash scripts/01_get_data.sh        # <1 секунды

# Проверка
ls -lh data/input/
```

**Что уже в репозитории:**
- `data/input/demo_variants.vcf` — 20 SNV/indel вариантов (GRCh38)
- `data/input/demo_sv.vcf` — 5 структурных вариантов
- `data/input/demo_variants_vep.tsv` — Предобработанная VEP аннотация

**Опционально: VEP cache для offline-режима (~15 ГБ)**
```bash
# Установите в workflow/config.sh:
# VEP_MODE="offline"
# VEP_CACHE_DIR="/shared/vep_cache"   # общая папка на сервере
bash scripts/02_get_vep_cache.sh      # ~30–60 минут, ~15 ГБ
```

Без VEP cache скрипты автоматически используют предобработанный файл.

### 5. Проверка системы

```bash
micromamba activate ngs-annotation
cd genome_annotation && bash scripts/00_check_system.sh

micromamba activate ngs-variants
cd ../variant_calling && bash scripts/00_check_system.sh

micromamba activate ngs-interpret
cd ../variant_interpretation && bash scripts/00_check_system.sh
```

Все пункты должны показывать `[OK]` / `✓`.

## Управление пространством

### Размер загруженных данных

| Компонент | Размер | Примечание |
|-----------|--------|-----------|
| genome_annotation окружение | ~1.2 ГБ | Prokka, DIAMOND, HMMER, Python |
| variant_calling окружение | ~1 ГБ | samtools, bcftools, freebayes, bwa |
| variant_interpretation окружение | ~800 МБ | VEP, bcftools, samtools, pandas |
| B. subtilis 168 геном | ~4 МБ | Скачивается в модули 1 и 2 |
| Swiss-Prot DIAMOND индекс | ~50 МБ | Только для genome_annotation |
| Pfam подмножество | ~5 МБ | Только для genome_annotation |
| Симулированные прочтения (BAM) | ~300 МБ | Только для variant_calling |
| Демо VCF и VEP TSV | ~50 КБ | Включены в репо (variant_interpretation) |
| VEP cache GRCh38 v110 (опционально) | ~15 ГБ | Только для variant_interpretation offline |
| **Итого (без VEP cache)** | **~3.5 ГБ** | На систему |
| **Итого (с VEP cache)** | **~18.5 ГБ** | На систему |

### Очистка

```bash
# Удалить результаты (outputs) но сохранить данные и окружения
make -C genome_annotation clean
make -C variant_calling clean
make -C variant_interpretation clean

# Удалить всё включая скачанные данные
make -C genome_annotation clean-all
make -C variant_calling clean-all
make -C variant_interpretation clean-all

# Удалить conda-окружения
micromamba env remove -n ngs-annotation
micromamba env remove -n ngs-variants
micromamba env remove -n ngs-interpret
```

## Тестовый запуск перед занятием

Перед каждым занятием рекомендуется проверить, что всё работает:

```bash
# Практикум 1: genome_annotation
micromamba activate ngs-annotation
cd genome_annotation
cp workflow/config.sh.example workflow/config.sh
bash scripts/10_run_prokka.sh      # ~2 минуты

# Практикум 2: variant_calling
micromamba activate ngs-variants
cd ../variant_calling
cp workflow/config.sh.example workflow/config.sh
bash scripts/10_prep_reference.sh  # ~10 секунд
bash scripts/11_prep_bam.sh        # ~10 секунд

# Практикум 3: variant_interpretation
micromamba activate ngs-interpret
cd ../variant_interpretation
cp workflow/config.sh.example workflow/config.sh
bash workflow/run_all.sh           # ~20 секунд (использует precomputed VEP)
```

## Проблемы и решения

### Проблема: "command not found: conda / micromamba"

**Решение:** Проверьте, что Micromamba/Conda установлены и добавлены в PATH:
```bash
which micromamba
# Если пусто:
echo 'export PATH="${MAMBA_ROOT_PREFIX}/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Проблема: "Permission denied" при запуске скриптов

**Решение:** Убедитесь, что скрипты исполняемые:
```bash
chmod +x genome_annotation/scripts/*.sh
chmod +x genome_annotation/env/*.sh
chmod +x variant_calling/scripts/*.sh
chmod +x variant_calling/env/*.sh
chmod +x variant_interpretation/scripts/*.sh
chmod +x variant_interpretation/env/*.sh
```

### Проблема: Не хватает места на диске

**Решение:** Проверьте свободное место и удалите ненужные данные:
```bash
df -h                          # Проверить свободное место
make -C genome_annotation clean-all
make -C variant_calling clean-all
rm -rf ~/miniconda3/pkgs/*      # Очистить кэш conda
```

### Проблема: "Network error" при загрузке данных

**Решение:** Повторите попытку позже или вручную загрузите:
```bash
# Для genome_annotation
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
# Распакуйте и скопируйте в data/input/assembly.fasta
```

Консультируйтесь с **shared/TROUBLESHOOTING_COMMON.md** для других проблем.

## Часто задаваемые вопросы

### Q: Можно ли запустить все три практикума на одной системе?

**A:** Да, все три полностью независимы. Каждый использует собственное conda-окружение (`ngs-annotation`, `ngs-variants`, `ngs-interpret`) и свою директорию данных.

### Q: Сколько времени занимает подготовка?

**A:** ~30 минут для практикумов 1–2 при первом запуске. Практикум 3 не требует загрузки данных — демо-файлы уже в репозитории.

### Q: Нужен ли интернет во время занятия?

**A:** Для практикумов 1–2: нет (данные загружаются заранее). Для практикума 3: нет — VEP REST API опционален, при недоступности скрипты автоматически используют предобработанный файл `data/input/demo_variants_vep.tsv`.

### Q: Можно ли использовать стандартный conda вместо micromamba?

**A:** Да, просто замените `micromamba` на `conda`. Файлы `env/environment.yml` совместимы с обоими.

### Q: Для практикума 3 нужен VEP cache?

**A:** Нет. VEP cache (~15 ГБ) нужен только для режима `VEP_MODE=offline`. По умолчанию используется REST API или предобработанный файл. Студенты могут пройти весь практикум без VEP cache.
