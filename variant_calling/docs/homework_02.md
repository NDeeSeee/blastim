# Домашнее задание 2 (~ 60 минут): Variant Defense

**Организм:** *Bacillus subtilis* 168 (модельный штамм)

## Цель

После запуска пайплайна (freebayes → bcftools filter → stats) нужно вынести вердикт по 5 персональным вариантам: настоящий вариант или артефакт.

Это задача, которую реально выполняют в лаборатории при ручной проверке результатов variant calling.

---

## Как сдаём (Git-практика)

1. Сделайте **fork** репозитория преподавателя к себе в GitHub.
2. Клонируйте свой fork, создайте ветку:
   ```bash
   git checkout -b hw2-<github_username>
   ```
3. Запустите назначение персональных вариантов:
   ```bash
   export GITHUB_USER=<ваш_github_username>
   python scripts/90_homework_assign_variants.py \
     --vcf outputs/<ваш_run_папка>/variants/raw.vcf.gz \
     --out homework/submissions/$GITHUB_USER/variants_$GITHUB_USER.txt
   ```
4. Скопируйте шаблон в свою папку:
   ```bash
   mkdir -p homework/submissions/$GITHUB_USER
   cp homework/TEMPLATE_variant_defense.md \
      homework/submissions/$GITHUB_USER/defense_$GITHUB_USER.md
   ```
5. Заполните файл, затем:
   ```bash
   git add homework/submissions/$GITHUB_USER/
   git commit -m "HW2: $GITHUB_USER"
   git push origin hw2-$GITHUB_USER
   ```
6. Откройте **Pull Request** в основной репозиторий преподавателя.

---

## Индивидуализация

Каждому студенту назначается персональный набор из 5 вариантов (детерминированно по GitHub username). Скрипт `90_homework_assign_variants.py` выбирает сбалансированный набор:
- 1-2 надёжных SNP (высокий QUAL, высокий DP)
- 1 indel
- 1-2 пограничных/низкокачественных варианта

---

## Задание: Variant Defense

**Файл для сдачи:** `homework/submissions/<user>/defense_<user>.md`

Для каждого из 5 назначенных вариантов заполните:

### 1. Позиция и основные данные
- CHROM:POS, REF → ALT
- QUAL, DP, AO, RO
- Тип: SNP / insertion / deletion

### 2. Доказательства ЗА (вариант настоящий)
- Высокий QUAL?
- Достаточная глубина (DP)?
- Баланс прочтений ALT по цепям (SAF/SAR)?
- AO/DP (частота ALT аллеля)?

### 3. Доказательства ПРОТИВ (вариант = артефакт)
- Низкий QUAL?
- Strand bias (SAP > 20)?
- Position bias (RPP)?
- Низкий mapping quality (MQM)?
- Вблизи indel (proximity to other variants)?

### 4. Вердикт
- [ ] **ACCEPT** — вариант настоящий, высокая уверенность
- [ ] **ACCEPT with caution** — скорее настоящий, но есть сомнения
- [ ] **REJECT** — скорее артефакт, не прошёл бы ручную проверку
- [ ] **UNCLEAR** — недостаточно данных для решения

### 5. Предложение фильтра
Какое правило `bcftools view -i '...'` отделило бы этот вариант от артефактов?

---

## Как извлечь данные для анализа

```bash
# Полная строка VCF для конкретной позиции
bcftools view outputs/<run>/variants/raw.vcf.gz NC_000964.3:12345-12345

# Ключевые поля
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AO\t%INFO/RO\t%INFO/SAF\t%INFO/SAR\t%INFO/SRP\t%INFO/SAP\t%INFO/MQM\n' \
    outputs/<run>/variants/raw.vcf.gz \
    -r NC_000964.3:12345-12345

# Покрытие в позиции
samtools depth data/input/reads_sorted.bam -r NC_000964.3:12345-12345
```

Справочник полей VCF: `docs/vcf_glossary.md`

---

## Самопроверка (опционально)

```bash
python scripts/91_homework_check.py --user $GITHUB_USER
```

Или через Make:

```bash
make homework-check
```

---

## Критерии оценивания (5 баллов)

| Балл | Критерий |
|------|----------|
| 1 | Корректное извлечение данных из VCF (QUAL, DP, AO/RO, strand info) |
| 2 | Логичная аргументация ЗА и ПРОТИВ для каждого варианта |
| 3 | Вердикт обоснован и согласуется с доказательствами |
| 4 | Предложение фильтра конкретное и применимое |
| 5 | Понятный и структурированный отчёт |
