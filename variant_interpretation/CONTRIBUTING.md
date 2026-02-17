# Руководство по сдаче домашнего задания (ДЗ3)

## Форк и ветка

1. **Форкните** репозиторий на GitHub (кнопка Fork)
2. **Склонируйте** свой форк:
   ```bash
   git clone https://github.com/<ваш_github>/blastim_ngs
   cd blastim_ngs/variant_interpretation
   ```
3. **Создайте ветку** `hw3-<ваш_github>`:
   ```bash
   git checkout -b hw3-<ваш_github>
   ```

## Подготовка файлов

```bash
export GITHUB_USER=<ваш_github>

# Получить назначенные варианты
python3 scripts/90_assign_variants.py \
    --vep_tsv data/input/demo_variants_vep.tsv \
    --sv_vcf  data/input/demo_sv.vcf \
    --out     homework/submissions/${GITHUB_USER}/variants_${GITHUB_USER}.txt

# Создать отчёт из шаблона
mkdir -p homework/submissions/${GITHUB_USER}
cp homework/TEMPLATE_hw3.md \
   homework/submissions/${GITHUB_USER}/hw3_${GITHUB_USER}.md

# Заполнить hw3_<user>.md (см. docs/homework_03.md)
# ...

# Проверить оформление
make homework-check GITHUB_USER=${GITHUB_USER}
```

## Структура сдаваемых файлов

```
homework/submissions/<GITHUB_USER>/
├── variants_<GITHUB_USER>.txt    # Назначенные варианты (генерируется скриптом)
└── hw3_<GITHUB_USER>.md          # Защита вариантов (заполняется вручную)
```

## Коммит и PR

```bash
git add homework/submissions/${GITHUB_USER}/
git commit -m "HW3: variant defense [${GITHUB_USER}]"
git push origin hw3-${GITHUB_USER}
```

Откройте PR на GitHub:
- **Base:** `main` (в оригинальном репозитории)
- **Head:** `hw3-<ваш_github>` (ваш форк)
- **Title:** `HW3: Variant Defense — <GITHUB_USER>`

## Checklist перед отправкой

- [ ] `variants_<user>.txt` создан (`90_assign_variants.py`)
- [ ] `hw3_<user>.md` создан из шаблона
- [ ] Все 8 вариантов заполнены (геномная позиция, VEP, QUAL/DP, вердикт)
- [ ] Вердикт для каждого варианта: выбран `[x]` и написано обоснование ≥3 предложений
- [ ] Нет незаполненных плейсхолдеров `<...>`
- [ ] `make homework-check GITHUB_USER=<user>` проходит без ошибок
- [ ] PR создан с правильным названием

## Формат обоснования вердикта

Каждый вердикт должен содержать:
1. **Тип изменения** — что происходит на молекулярном уровне (missense, frameshift, DEL...)
2. **Доказательства патогенности** — IMPACT, ClinVar, gnomAD_AF, QUAL/DP
3. **ACMG критерии** — какие критерии применимы и почему

Пример хорошего обоснования:
> rs28934574 (TP53 R175H) — миссенс-вариант с IMPACT=HIGH. В ClinVar классифицирован как
> Pathogenic; gnomAD_AF=0.000003 (практически отсутствует в популяции). QUAL=60, DP=45 —
> высококачественный вызов. По ACMG: PS1 (тот же аминокислотный вариант известен как
> патогенный), PM2 (крайне редкий в gnomAD), PP3 (предсказатели патогенности согласны).
> Вердикт: **Патогенный**.

## Политика

- **Самостоятельная работа.** Обсуждать подходы можно, но текст обоснований — ваш.
- **Плагиат** (одинаковые обоснования у разных студентов) → ноль баллов обоим.
- **Не изменяйте** файлы вне `homework/submissions/<ваш_github>/`.
- Вопросы → Issues на GitHub с тегом `[HW3]`.
