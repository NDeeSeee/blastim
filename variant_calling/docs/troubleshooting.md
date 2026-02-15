# Устранение неполадок

Частые проблемы и их решения.

## Окружение

### «command not found: freebayes» (или любой инструмент)

**Причина**: Conda-окружение не активировано.

**Решение**:
```bash
micromamba activate ngs-variants
# или
conda activate ngs-variants
```

Если окружение не существует:
```bash
bash env/micromamba_create.sh
```

### «conda: command not found»

**Причина**: Менеджер пакетов не установлен.

**Решение**: Установите micromamba:
```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

Затем перезапустите шелл и попробуйте снова.

### Создание окружения зависает на «solving environment»

**Причина**: Solver conda работает медленно.

**Решение**: Используйте micromamba (предпочтительно) или mamba вместо conda:
```bash
# Установить mamba в базовое окружение
conda install -n base -c conda-forge mamba
mamba env create -f env/environment.yml
```

## Данные

### «Reference not found»

**Решение**: Запустите скрипт подготовки данных:
```bash
bash scripts/01_get_data.sh
```

### «BAM file not found»

**Решение**: Также запустите `01_get_data.sh` — он загружает геном, симулирует прочтения и создаёт BAM.

### «BAM index not found»

**Причина**: Индекс `.bai` отсутствует или устарел.

**Решение**:
```bash
samtools index data/input/reads_sorted.bam
```

### wget не работает (нет интернета)

**Причина**: На HPC нет исходящего доступа в интернет.

**Решение**: Загрузите на машине с интернетом и перенесите:
```bash
# На машине с интернетом
bash scripts/01_get_data.sh

# Перенос на HPC
rsync -avP data/ hpc:/path/to/blastim_ngs/variant_calling/data/
```

## Запуск скриптов

### «Config not found»

**Решение**:
```bash
cp workflow/config.sh.example workflow/config.sh
```

### «Permission denied» при запуске скриптов

**Решение**: Запускайте через `bash`:
```bash
bash scripts/10_prep_reference.sh
```

### freebayes зависает или работает очень долго

**Возможные причины**:
1. Большой геном или BAM: freebayes работает на одном ядре
2. Слишком высокое покрытие в отдельных регионах

**Решение**:
```bash
# Вариант A: ограничить по региону
freebayes --fasta-reference data/input/reference.fasta \
    --ploidy 1 --bam data/input/reads_sorted.bam \
    --region "$(head -1 data/input/reference.fasta.fai | cut -f1):1-500000" \
    > outputs/run_latest/variants/raw.vcf

# Вариант B: попробовать параллельный запуск (freebayes-parallel)
# Требует скрипт fasta_generate_regions.py из freebayes
```

### freebayes: «Could not open reference»

**Причина**: Отсутствует `.fai` индекс.

**Решение**:
```bash
samtools faidx data/input/reference.fasta
```

### Неправильная плоидность (ploidy)

**Симптомы**: Для бактерий генотипы `0/1` вместо `0` или `1`.

**Причина**: Плоидность не установлена (по умолчанию diploid).

**Решение**: В `config.sh` должно быть `PLOIDY=1`. Перезапустите freebayes.

### Пустой VCF (0 вариантов)

**Возможные причины**:
1. BAM-файл пустой или не содержит выравненных прочтений:
   ```bash
   samtools flagstat data/input/reads_sorted.bam
   ```
2. Референсный геном не совпадает с тем, на который выравнивали:
   ```bash
   samtools idxstats data/input/reads_sorted.bam
   ```
3. Все прочтения имеют MAPQ=0 (маппинг на множественные позиции)

### bgzip vs gzip

**Проблема**: `tabix` не работает с файлом, сжатым обычным `gzip`.

**Причина**: `bgzip` создаёт block-compressed формат с произвольным доступом. Обычный `gzip` не поддерживает индексирование.

**Решение**: Пересжать файл:
```bash
# Если файл уже сжат gzip:
gunzip raw.vcf.gz
bgzip raw.vcf
tabix -p vcf raw.vcf.gz
```

**Как проверить**: `htslib` файл vs `gzip` файл:
```bash
file raw.vcf.gz
# Должно быть: "raw.vcf.gz: gzip compressed data" (оба выглядят одинаково)
# Но tabix откроет только bgzip-файлы
```

### «Disk quota exceeded»

**Причина**: Исчерпана квота диска на HPC.

**Решение**:
```bash
# Проверить использование диска
du -sh data/ outputs/

# Очистить старые результаты
make clean

# Проверить квоту
quota -s  # или: df -h .
```

## Проблемы с результатами

### Слишком мало вариантов (< 1000)

**Возможные причины**:
1. Низкое покрытие BAM: проверьте `samtools depth`
2. Неправильный референсный геном
3. Ошибка при симуляции прочтений (проверьте параметры wgsim)

### Слишком много вариантов (> 10000)

**Возможные причины**:
1. Покрытие слишком высокое — увеличьте `MAX_DP` в `config.sh`
2. Высокая частота мутаций при симуляции (проверьте `-r` в 01_get_data.sh)

### Ti/Tv ratio слишком низкий (< 1.5) или слишком высокий (> 3.0)

Ожидаемый Ti/Tv для бактериальных данных: ~2.0-2.5.
- Слишком низкий может указывать на артефакты или ошибки в данных.
- Слишком высокий может быть связан с малым числом вариантов (стохастический эффект).

## Общие советы

- Всегда активируйте окружение перед запуском скриптов.
- Запустите `bash scripts/00_check_system.sh` для диагностики большинства проблем.
- Читайте сообщения скриптов — они указывают, что запускать дальше.
- Если шаг завершился с ошибкой, исправьте проблему и перезапустите только этот шаг (не весь рабочий процесс).
