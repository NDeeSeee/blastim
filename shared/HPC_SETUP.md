# Подготовка HPC-системы для курса

Руководство для преподавателя: как подготовить HPC (или локальную систему) для проведения практикумов.

## Быстрая подготовка (15 минут)

```bash
# 1. Клонируйте репозиторий
git clone git@github.com:NDeeSeee/blastim_ngs.git
cd blastim_ngs

# 2. Создайте оба conda-окружения
bash genome_annotation/env/micromamba_create.sh
bash variant_calling/env/micromamba_create.sh

# 3. Подготовьте данные для обоих модулей (один раз)
micromamba activate ngs-annotation
bash genome_annotation/scripts/01_get_data.sh
bash genome_annotation/scripts/02_make_diamond_db.sh
bash genome_annotation/scripts/03_get_pfam.sh

micromamba activate ngs-variants
bash variant_calling/scripts/01_get_data.sh

# 4. Проверьте, что всё работает
bash genome_annotation/scripts/00_check_system.sh
bash variant_calling/scripts/00_check_system.sh
```

## Требования к системе

### Минимальные требования

| Параметр | Значение |
|----------|----------|
| **Процессоры** | 4+ cores (рекомендуется 8+) |
| **Оперативная память** | 8 ГБ (рекомендуется 16+) |
| **Дисковое пространство** | 10 ГБ свободного |
| **ОС** | Linux / macOS / Windows WSL2 |
| **Сеть** | Интернет для загрузки данных |

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
# genome_annotation
bash genome_annotation/env/micromamba_create.sh
# Это создаст окружение ngs-annotation

# variant_calling
bash variant_calling/env/micromamba_create.sh
# Это создаст окружение ngs-variants
```

#### Вариант B: Conda/Mamba

```bash
cd genome_annotation
conda env create -f env/environment.yml
conda activate ngs-annotation

cd ../variant_calling
conda env create -f env/environment.yml
conda activate ngs-variants
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

### 5. Проверка системы

```bash
cd genome_annotation
bash scripts/00_check_system.sh

cd ../variant_calling
bash scripts/00_check_system.sh
```

Все пункты должны показывать `[OK]`.

## Управление пространством

### Размер загруженных данных

| Компонент | Размер | Примечание |
|-----------|--------|-----------|
| genome_annotation окружение | ~1.2 ГБ | Prokka, DIAMOND, HMMER, Python |
| variant_calling окружение | ~1 ГБ | samtools, bcftools, freebayes, bwa |
| B. subtilis 168 геном | ~4 МБ | Скачивается в оба модуля |
| Swiss-Prot DIAMOND индекс | ~50 МБ | Только для genome_annotation |
| Pfam подмножество | ~5 МБ | Только для genome_annotation |
| Симулированные прочтения (BAM) | ~300 МБ | Только для variant_calling |
| **Итого** | **~2.5-3 ГБ** | На систему |

### Очистка

```bash
# Удалить результаты (outputs) но сохранить данные и окружения
make -C genome_annotation clean
make -C variant_calling clean

# Удалить всё включая скачанные данные
make -C genome_annotation clean-all
make -C variant_calling clean-all

# Удалить conda-окружения
micromamba env remove -n ngs-annotation
micromamba env remove -n ngs-variants
```

## Тестовый запуск перед занятием

Перед каждым занятием рекомендуется проверить, что всё работает:

```bash
# Для genome_annotation
micromamba activate ngs-annotation
cd genome_annotation
cp workflow/config.sh.example workflow/config.sh
bash scripts/10_run_prokka.sh      # Проверяет базовую работу (~2 минуты)

# Для variant_calling
micromamba activate ngs-variants
cd ../variant_calling
cp workflow/config.sh.example workflow/config.sh
bash scripts/10_prep_reference.sh  # Быстрая проверка (~10 секунд)
bash scripts/11_prep_bam.sh        # Проверка данных (~10 секунд)
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

### Q: Можно ли запустить оба практикума на одной системе?

**A:** Да, оба полностью независимы. Каждый использует собственное окружение и данные.

### Q: Сколько времени занимает подготовка?

**A:** ~30 минут при первом запуске (включая загрузку данных). На последующих запусках 01_get_data.sh перепроверяет данные (~2 минуты).

### Q: Нужна ли Internet вся время?

**A:** Нет, только для подготовки (01_get_data.sh). После этого все практикумы работают offline.

### Q: Можно ли использовать стандартный conda вместо micromamba?

**A:** Да, просто используйте `conda` вместо `micromamba activate`. Обе работают одинаково.
