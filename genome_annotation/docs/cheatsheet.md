# Шпаргалка по командам

Все команды для копирования и вставки во время занятия. Выполняйте из директории практикума (`genome_annotation/`).

## Настройка

```bash
# Активация окружения
micromamba activate ngs-annotation
# или: conda activate ngs-annotation

# Копирование конфигурации
cp workflow/config.sh.example workflow/config.sh

# Проверка системы
bash scripts/00_check_system.sh
```

## Запуск рабочего процесса

```bash
# Всё сразу
bash workflow/run_all.sh

# Или пошагово
bash scripts/10_run_prokka.sh
bash scripts/20_run_diamond.sh
bash scripts/30_run_hmmscan.sh
bash scripts/40_sanity_checks.sh
bash scripts/50_pick_3_genes.sh
```

## Изучение результатов Prokka

```bash
# Найти директорию результатов
OUTDIR=$(ls -td outputs/run_* | head -1)

# Посмотреть сводку
cat $OUTDIR/prokka/ANNOT.txt

# Подсчитать признаки по типам
grep -v "^#" $OUTDIR/prokka/ANNOT.gff | awk -F'\t' '{print $3}' | sort | uniq -c | sort -rn

# Посмотреть первые 5 записей CDS в GFF
grep "CDS" $OUTDIR/prokka/ANNOT.gff | head -5

# Посмотреть таблицу признаков
head -20 $OUTDIR/prokka/ANNOT.tsv

# Статистика по белкам
seqkit stats $OUTDIR/prokka/ANNOT.faa

# Посмотреть первые 3 белковые последовательности
seqkit head -n 3 $OUTDIR/prokka/ANNOT.faa

# Найти конкретный белок по ключевому слову
grep "ribosomal" $OUTDIR/prokka/ANNOT.tsv
```

## Изучение результатов DIAMOND

```bash
# Посмотреть первые 10 попаданий
head -10 $OUTDIR/diamond/blastp_results.tsv | column -t

# Подсчитать белки с попаданиями
cut -f1 $OUTDIR/diamond/blastp_results.tsv | sort -u | wc -l

# Лучшее попадание для каждого белка (наибольший bitscore)
sort -k1,1 -k12,12rn $OUTDIR/diamond/blastp_results.tsv | sort -k1,1 -u | head -10 | column -t

# Распределение процента идентичности (лучшие попадания)
sort -k1,1 -k12,12rn $OUTDIR/diamond/blastp_results.tsv | sort -k1,1 -u | \
  awk -F'\t' '{printf "%d\n", $3}' | sort -n | uniq -c | sort -k2 -n

# Поиск конкретного белка
grep "ANNOT_00001" $OUTDIR/diamond/blastp_results.tsv | column -t

# Найти попадания по ключевому слову (например, «kinase»)
grep -i "kinase" $OUTDIR/diamond/blastp_results.tsv | head -5 | column -t
```

## Изучение результатов hmmscan

```bash
# Посмотреть первые 10 найденных доменов (пропустить строки комментариев)
grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | head -10

# Подсчитать общее число найденных доменов
grep -cv "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt

# Топ-20 наиболее частых доменов
grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $1}' | sort | uniq -c | sort -rn | head -20

# Белки с наибольшим числом доменов
grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $4}' | sort | uniq -c | sort -rn | head -10

# Найти попадания для конкретного белка
grep "ANNOT_00042" $OUTDIR/hmmer/hmmscan_domtblout.txt

# Найти конкретный домен (например, ABC-транспортер)
grep -i "ABC" $OUTDIR/hmmer/hmmscan_domtblout.txt | head -5
```

## Комбинированный анализ

```bash
# Белки с попаданиями DIAMOND, но без доменов Pfam
comm -23 \
  <(cut -f1 $OUTDIR/diamond/blastp_results.tsv | sort -u) \
  <(grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $4}' | sort -u) \
  | head -10

# Белки с доменами Pfam, но без попаданий DIAMOND
comm -13 \
  <(cut -f1 $OUTDIR/diamond/blastp_results.tsv | sort -u) \
  <(grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $4}' | sort -u) \
  | head -10

# Белки БЕЗ какой-либо аннотации
comm -23 \
  <(grep "^>" $OUTDIR/prokka/ANNOT.faa | sed 's/>//' | awk '{print $1}' | sort) \
  <(cat <(cut -f1 $OUTDIR/diamond/blastp_results.tsv) \
        <(grep -v "^#" $OUTDIR/hmmer/hmmscan_domtblout.txt | awk '{print $4}') | sort -u) \
  | wc -l
```

## Полезные команды seqkit

```bash
# Распределение длин белков
seqkit fx2tab -l $OUTDIR/prokka/ANNOT.faa | awk '{print $NF}' | sort -n | tail -5

# Найти самый длинный белок
seqkit fx2tab -l $OUTDIR/prokka/ANNOT.faa | sort -t$'\t' -k2 -rn | head -1

# Найти самый короткий белок
seqkit fx2tab -l $OUTDIR/prokka/ANNOT.faa | sort -t$'\t' -k2 -n | head -1

# Извлечь конкретный белок по ID
seqkit grep -p "ANNOT_00042" $OUTDIR/prokka/ANNOT.faa

# Подсчитать белки короче 100 аминокислот
seqkit fx2tab -l $OUTDIR/prokka/ANNOT.faa | awk '$NF < 100' | wc -l
```
