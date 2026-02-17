# Инструкция для преподавателя: Практикум 3

## Обзор

| Параметр | Значение |
|----------|---------|
| Лекция | 3 из N |
| Продолжительность | 90–120 минут |
| Размер группы | 15–30 человек |
| Уровень | Начинающий–средний |
| Предварительные знания | Практикумы 1 (genome_annotation) и 2 (variant_calling) |

---

## Подготовка (за 1–3 дня до занятия)

### 1. Проверьте сервер

```bash
bash scripts/00_check_system.sh
```

Убедитесь:
- Python ≥ 3.10, pandas, matplotlib установлены
- bcftools, samtools, tabix доступны
- Интернет: `curl https://rest.ensembl.org` возвращает ответ (для VEP REST)

### 2. Создайте conda-окружение на сервере

```bash
cd variant_interpretation/
make env
conda activate ngs-interpret   # или micromamba activate ngs-interpret
```

### 3. (Опционально) Подготовьте VEP cache для offline-режима

> Требуется ~15 ГБ дискового пространства и 30–60 минут загрузки.
> Без этого шага практикум работает через REST API или предобработанный файл.

```bash
# В workflow/config.sh:
VEP_MODE="offline"
VEP_CACHE_DIR="/shared/vep_cache"   # укажите общую директорию

# Загрузить (только один раз на сервер):
bash scripts/02_get_vep_cache.sh
```

### 4. Проверьте входные данные

```bash
bash scripts/01_get_data.sh
```

Все файлы включены в репозиторий — дополнительных загрузок не требуется.

### 5. Сделайте тестовый прогон

```bash
cp workflow/config.sh.example workflow/config.sh
make all
```

Убедитесь, что `outputs/run_*/report.md` создан.

---

## Структура занятия (90 мин)

| Время | Тема | Слайды |
|-------|------|--------|
| 0–10 | Обзор практикума, цели, данные | 1–3 |
| 10–25 | Формат VCF: QUAL, DP, AF, INFO | 4–8 |
| 25–40 | VEP: что такое consequence, IMPACT, HGVS | 9–14 |
| 40–55 | Приоритизация: scoring, ClinVar, gnomAD | 15–20 |
| 55–70 | SV: SVTYPE, DEL/DUP/INV, размеры | 21–25 |
| 70–80 | Современные инструменты (AlphaMissense, SpliceAI) | 26–28 |
| 80–90 | ДЗ3: условия, вопросы | — |

---

## Ход практикума (шаг за шагом)

### Шаг 1: Ориентация (5 мин)
```bash
make check
bash scripts/01_get_data.sh
```
Показываем структуру директории. Объясняем, что демо-данные уже внутри репо.

### Шаг 2: Изучаем VCF (10 мин)
```bash
grep -v "^#" data/input/demo_variants.vcf | head -5 | column -t
grep "^##INFO" data/input/demo_variants.vcf
```
Разбираем поля QUAL, DP, AF, gnomAD_AF, ClinVar.

### Шаг 3: Запускаем VEP (5 мин)
```bash
make vep
```
Пока VEP работает (или сразу получает предобработанный файл), объясняем что происходит.

### Шаг 4: Изучаем VEP output (10 мин)
```bash
grep -v "^#" outputs/run_*/vep/vep_output.tsv | head -3
grep "IMPACT=HIGH" outputs/run_*/vep/vep_output.tsv
```
Разбираем: consequence, IMPACT, HGVSc/p, gnomADe_AF, ClinVar_CLNSIG.

### Шаг 5: Приоритизация (10 мин)
```bash
make prioritize
column -t outputs/run_*/prioritized.tsv | head -10
```
Объясняем алгоритм скоринга.

### Шаг 6: Отчёт (5 мин)
```bash
cat outputs/run_*/report.md
```

### Шаг 7: SV демо (10 мин)
```bash
make sv-demo
```

### Шаг 8: Современные инструменты (10 мин)
```bash
make modern
```
Обсуждаем AlphaMissense, SpliceAI, gnomAD v4.

---

## Назначение вариантов для ДЗ3

```bash
# После занятия студенты выполняют сами:
export GITHUB_USER=<username>
python3 scripts/90_assign_variants.py \
    --vep_tsv data/input/demo_variants_vep.tsv \
    --sv_vcf  data/input/demo_sv.vcf \
    --out homework/submissions/${GITHUB_USER}/variants_${GITHUB_USER}.txt
```

Алгоритм детерминирован: один и тот же username → одни и те же варианты.
20 SNV/indel + 5 SV → уникальные наборы по 8 вариантов на студента.

---

## Разбор ДЗ3 (Лекция 4)

Типичные ошибки студентов:
1. **Не различают IMPACT и consequence** — особо объяснить
2. **Игнорируют DP < 10** — низкое покрытие = ненадёжный вызов
3. **Путают ClinVar источник с ACMG классом** — ClinVar это база, ACMG — их вывод
4. **Для SV не указывают размер** — SVLEN обязателен
5. **Плейсхолдеры `<...>` не заменены** — проверка `make homework-check` ловит это

---

## Часто задаваемые вопросы

**Q: VEP REST API слишком медленный или недоступен.**
A: Скрипт автоматически использует `data/input/demo_variants_vep.tsv`. Практикум продолжается.

**Q: Ошибка при `make env` — conda не найден.**
A: Установите micromamba. Инструкция: mamba.readthedocs.io

**Q: Сколько времени занимает VEP REST?**
A: 30–120 с для 20 вариантов. Если >3 мин — проблема с сетью, переключитесь на precomputed.

**Q: Почему `prioritized.tsv` содержит только 12–18 строк вместо 20?**
A: Фильтрация по `CANONICAL=YES`: только один транскрипт на вариант.

**Q: Студент получил другой список вариантов, чем сосед.**
A: Правильно! Каждый username → уникальный набор вариантов (SHA256-seed).

---

## Сброс и повтор

```bash
# Удалить все outputs (данные сохраняются)
make clean

# Полный сброс
make clean-all
cp workflow/config.sh.example workflow/config.sh

# Перезапустить
make all
```

---

## Дисковые требования

| Компонент | Размер |
|-----------|--------|
| Репозиторий (данные практикума) | ~500 KB |
| conda-окружение ngs-interpret | ~800 MB |
| VEP cache GRCh38 v110 (опционально) | ~15 GB |
| Outputs одного прогона | ~200 KB |
