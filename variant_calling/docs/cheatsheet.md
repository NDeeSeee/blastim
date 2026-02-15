# Шпаргалка по командам

Все команды для копирования и вставки во время занятия. Выполняйте из директории практикума (`variant_calling/`).

## Настройка

```bash
# Активация окружения
micromamba activate ngs-variants
# или: conda activate ngs-variants

# Копирование конфигурации
cp workflow/config.sh.example workflow/config.sh

# Проверка системы
bash scripts/00_check_system.sh
```

## Подготовка данных

```bash
# Загрузка референса B. subtilis 168 и симуляция прочтений с помощью wgsim
# ВАЖНО: этот шаг выполняет преподаватель один раз перед практикумом
bash scripts/01_get_data.sh

# Проверка что данные скачаны и выровнены
ls -lh data/input/reference.fasta
ls -lh data/input/reads_sorted.bam
ls -lh data/truth/
```

## Запуск рабочего процесса

```bash
# Всё сразу
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

## Изучение BAM-файла

```bash
# Проверка целостности
samtools quickcheck data/input/reads_sorted.bam && echo "OK"

# Статистика выравнивания
samtools flagstat data/input/reads_sorted.bam

# Статистика по хромосомам/контигам
samtools idxstats data/input/reads_sorted.bam

# Средняя глубина покрытия
samtools depth -a data/input/reads_sorted.bam | awk '{s+=$3; n++} END{printf "%.1fx\n", s/n}'

# Глубина в конкретном регионе
samtools depth data/input/reads_sorted.bam -r "NC_000964.3:1000-2000"
```

## Изучение VCF-файлов

```bash
# Найти директорию результатов
OUTDIR=$(ls -td outputs/run_* | head -1)

# Посмотреть заголовок VCF (мета-информация)
bcftools view -h $OUTDIR/variants/raw.vcf.gz | head -30

# Первые 10 вариантов (без заголовка)
bcftools view -H $OUTDIR/variants/raw.vcf.gz | head -10

# Подсчитать общее число вариантов
bcftools view -H $OUTDIR/variants/raw.vcf.gz | wc -l

# Только SNP
bcftools view -v snps $OUTDIR/variants/filtered.vcf.gz | bcftools view -H | wc -l

# Только indels
bcftools view -v indels $OUTDIR/variants/filtered.vcf.gz | bcftools view -H | wc -l
```

## bcftools query — извлечение данных

```bash
# Позиция + REF + ALT + QUAL
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n' $OUTDIR/variants/filtered.vcf.gz | head -10

# Позиция + QUAL + DP + AO + RO
bcftools query -f '%CHROM\t%POS\t%QUAL\t%INFO/DP\t%INFO/AO\t%INFO/RO\n' $OUTDIR/variants/filtered.vcf.gz | head -10

# Варианты с QUAL > 1000 (самые надёжные)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\n' \
    -i 'QUAL>1000' $OUTDIR/variants/filtered.vcf.gz | head -10

# Варианты с низким покрытием (подозрительные)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\n' \
    -i 'INFO/DP<10' $OUTDIR/variants/raw.vcf.gz | head -10

# Подсчёт вариантов по типу (SNP vs indel)
bcftools view -v snps $OUTDIR/variants/filtered.vcf.gz | bcftools view -H | wc -l
bcftools view -v indels $OUTDIR/variants/filtered.vcf.gz | bcftools view -H | wc -l
```

## bcftools filter — примеры фильтров

```bash
# Оставить варианты с QUAL >= 50 и DP >= 10
bcftools view -i 'QUAL>=50 && INFO/DP>=10' $OUTDIR/variants/raw.vcf.gz | bcftools view -H | wc -l

# Оставить только SNP с высоким QUAL
bcftools view -v snps -i 'QUAL>=100' $OUTDIR/variants/raw.vcf.gz | bcftools view -H | head -5

# Пометить (soft-filter) варианты с низким QUAL
bcftools filter -s "LowQual" -e 'QUAL<20' $OUTDIR/variants/raw.vcf.gz | bcftools view -H | head -5
```

## bcftools stats — статистика

```bash
# Полная статистика
bcftools stats $OUTDIR/variants/filtered.vcf.gz | grep "^SN"

# Ti/Tv ratio
bcftools stats $OUTDIR/variants/filtered.vcf.gz | grep "^TSTV"

# Распределение QUAL
bcftools stats $OUTDIR/variants/filtered.vcf.gz | grep "^QUAL" | head -20
```

## Исследование конкретного варианта

```bash
# Посмотреть вариант в конкретной позиции
bcftools view $OUTDIR/variants/raw.vcf.gz NC_000964.3:12345-12345

# Посмотреть покрытие в этой позиции
samtools depth data/input/reads_sorted.bam -r "NC_000964.3:12345-12345"

# Посмотреть выравнивание в этой позиции (текстовый pileup)
samtools tview data/input/reads_sorted.bam data/input/reference.fasta -p NC_000964.3:12345
# (навигация: стрелки, q = выход)
```

## Сравнение raw vs filtered

```bash
# Сколько вариантов отфильтровано
RAW=$(bcftools view -H $OUTDIR/variants/raw.vcf.gz | wc -l)
FILT=$(bcftools view -H $OUTDIR/variants/filtered.vcf.gz | wc -l)
echo "Raw: $RAW  Filtered: $FILT  Removed: $((RAW - FILT))"

# Распределение QUAL у отфильтрованных вариантов
bcftools view -H $OUTDIR/variants/raw.vcf.gz | awk -F'\t' '{printf "%d\n", $6}' | sort -n | uniq -c | tail -20
```
