# Устранение неполадок

Частые проблемы и их решения.

## Окружение

### «command not found: prokka» (или любой инструмент)

**Причина**: Conda-окружение не активировано.

**Решение**:
```bash
micromamba activate ngs-annotation
# или
conda activate ngs-annotation
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

### «PackagesNotFoundError» при создании окружения

**Причина**: Проблема с порядком каналов или устаревшие данные каналов.

**Решение**:
```bash
# Обновить данные каналов
conda update --all -n base
# Или попробовать с явным указанием каналов
micromamba create -n ngs-annotation -c conda-forge -c bioconda -c defaults prokka seqkit diamond hmmer csvkit pigz wget curl
```

## Данные и базы данных

### «Assembly not found»

**Решение**: Запустите скрипт загрузки данных:
```bash
bash scripts/01_get_data.sh
```

Или проверьте, что сборка находится в `data/input/assembly.fasta`.

### «DIAMOND database not found»

**Решение**:
```bash
bash scripts/01_get_data.sh    # сначала загрузить белки
bash scripts/02_make_diamond_db.sh
```

### «HMM database not pressed (missing .h3m)»

**Причина**: `hmmpress` не был запущен или завершился с ошибкой.

**Решение**:
```bash
hmmpress -f data/db/pfam_subset.hmm
```

### wget не работает (нет интернета)

**Причина**: На HPC нет исходящего доступа в интернет.

**Решение**: Загрузите на машине с интернетом и перенесите:
```bash
# На машине с интернетом
bash scripts/01_get_data.sh
bash scripts/02_make_diamond_db.sh
bash scripts/03_get_pfam.sh --subset

# Перенос на HPC
rsync -avP data/ hpc:/path/to/blastim_ngs/genome_annotation/data/
```

### Загрузка Pfam не удалась («404 Not Found»)

**Причина**: URL Pfam мог измениться (EBI периодически реорганизует файлы).

**Решение**: Проверьте актуальный URL на https://www.ebi.ac.uk/interpro/download/Pfam/ и обновите URL в `scripts/03_get_pfam.sh`.

## Запуск скриптов

### «Config not found»

**Решение**:
```bash
cp workflow/config.sh.example workflow/config.sh
```

### «Permission denied» при запуске скриптов

**Решение**:
```bash
chmod +x scripts/*.sh workflow/*.sh
```

Или запускайте через `bash`:
```bash
bash scripts/10_run_prokka.sh
```

### Prokka: «Could not run prodigal»

**Причина**: Проблема с зависимостью Prokka.

**Решение**:
```bash
# Проверить наличие prodigal
which prodigal
prodigal -v

# Если отсутствует, установить
micromamba install -n ngs-annotation -c bioconda prodigal
```

### Prokka: «tbl2asn is missing or too old»

**Причина**: Распространённая проблема Prokka на новых системах.

**Решение**: Это предупреждение, а не ошибка. Prokka всё равно корректно создаёт GFF/FAA/FFN. Файлы .gbk и .sqn могут быть затронуты, но для данного курса они не нужны.

### DIAMOND: «Database was made with incompatible version»

**Причина**: Несоответствие версий DIAMOND при построении базы и поиске.

**Решение**: Пересоберите базу данных:
```bash
rm data/db/teaching_proteins.dmnd
bash scripts/02_make_diamond_db.sh
```

### hmmscan работает бесконечно

**Причина**: Используется полная Pfam-A (~20 тыс. профилей) вместо подмножества.

**Решение**: Переключитесь на режим подмножества:
```bash
# Отредактируйте workflow/config.sh
HMM_DB="data/db/pfam_subset.hmm"

# Или ограничьте число белков
bash scripts/30_run_hmmscan.sh --max-proteins 200
```

### hmmscan: «No such file or directory» для .h3m

**Причина**: База HMM не проиндексирована.

**Решение**:
```bash
hmmpress -f data/db/pfam_subset.hmm
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

### 0 попаданий DIAMOND

**Возможные причины**:
1. База данных пустая или повреждена: `diamond dbinfo --db data/db/teaching_proteins`
2. Файл запросов пустой: `wc -l outputs/run_*/prokka/ANNOT.faa`
3. Порог E-value слишком строгий (маловероятно при значении по умолчанию 1e-5)

### 0 попаданий hmmscan

**Возможные причины**:
1. В базе HMM 0 профилей: `grep -c "^ACC" data/db/pfam_subset.hmm`
2. База не проиндексирована: проверьте наличие файлов `.h3m`
3. E-value слишком строгий (маловероятно при настройках по умолчанию)

### Проверки качества показывают неожиданные числа

Сравните ваши результаты с `docs/expected_outputs.md`. Если значения отличаются более чем на 50%, проверьте:
1. Используется ли правильный файл сборки?
2. Заполнены ли базы данных?
3. Соответствуют ли версии инструментов спецификациям в environment.yml?

## Общие советы

- Всегда активируйте окружение перед запуском скриптов.
- Запустите `bash scripts/00_check_system.sh` для диагностики большинства проблем.
- Читайте сообщения скриптов -- они указывают, что запускать дальше.
- Если шаг завершился с ошибкой, исправьте проблему и перезапустите только этот шаг (не весь рабочий процесс).
