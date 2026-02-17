# Устранение проблем: Практикум по интерпретации вариантов

## Ошибка: VEP не найден

**Симптом:**
```
command not found: vep
```

**Решение:**
```bash
# Проверьте, активировано ли conda-окружение
conda activate ngs-interpret   # или micromamba activate ngs-interpret

# Если окружение не создано:
make env
```

---

## Ошибка: VEP REST API недоступен

**Симптом:**
```
WARN: REST API недоступен. Используем предобработанный файл.
```
или
```
urllib.error.URLError: <urlopen error [Errno 11001] getaddrinfo failed>
```

**Причина:** Нет доступа к интернету или rest.ensembl.org недоступен.

**Решение:** Скрипт автоматически переключится на предобработанный файл
`data/input/demo_variants_vep.tsv`. Практикум продолжается без изменений.

Для переключения в режим offline с cache:
```bash
# В workflow/config.sh:
VEP_MODE="offline"
VEP_CACHE_DIR="/shared/vep_cache"
```

---

## Ошибка: VEP cache не найден (offline режим)

**Симптом:**
```
WARN: VEP cache не найден: /shared/vep_cache
```

**Решение:**
```bash
# Загрузить cache (только для администратора/преподавателя):
bash scripts/02_get_vep_cache.sh

# Или переключиться обратно в REST режим:
# В workflow/config.sh:
VEP_MODE="rest"
```

---

## Ошибка: workflow/config.sh не найден

**Симптом:**
```
ERROR: workflow/config.sh не найден.
```

**Решение:**
```bash
cp workflow/config.sh.example workflow/config.sh
# Затем отредактируйте под свою среду
```

---

## Ошибка: Python pandas/matplotlib не найден

**Симптом:**
```
ModuleNotFoundError: No module named 'pandas'
```

**Решение:**
```bash
# Убедитесь, что conda-окружение активно
conda activate ngs-interpret

# Если нужно установить вручную:
pip install pandas matplotlib
```

---

## Ошибка: VEP TSV пустой или не содержит данных

**Симптом:** `prioritized.tsv` пуст или `report.md` содержит "0 вариантов"

**Проверка:**
```bash
wc -l outputs/run_*/vep/vep_output.tsv
grep -c "^#" outputs/run_*/vep/vep_output.tsv   # заголовочные строки
grep -vc "^#" outputs/run_*/vep/vep_output.tsv  # строки данных
```

**Причины:**
1. VEP завершился с ошибкой → проверьте лог
2. VCF пустой → `grep -vc "^#" data/input/demo_variants.vcf`
3. Используется неверный TSV → убедитесь, что путь правильный

**Решение:**
```bash
# Скопировать предобработанный файл вручную
mkdir -p outputs/run_latest/vep
cp data/input/demo_variants_vep.tsv outputs/run_latest/vep/vep_output.tsv
```

---

## Ошибка: OUTDIR не задан или outputs/ пуст

**Симптом:**
```
ls -td outputs/run_* 2>/dev/null | head -1
(нет вывода)
```

**Причина:** Пайплайн ещё не запускался или выходные файлы удалены.

**Решение:**
```bash
# Создать OUTDIR вручную и запустить шаги
export OUTDIR="outputs/run_manual"
bash scripts/10_run_vep.sh
```

---

## Ошибка: 21_make_report.py не находит VEP TSV

**Симптом:**
```
ERROR: VEP TSV not found: outputs/.../vep/vep_output.tsv
```

**Решение:**
```bash
# Убедитесь, что шаг 10 выполнен
bash scripts/10_run_vep.sh

# Или укажите путь явно
python3 scripts/21_make_report.py \
    --tsv data/input/demo_variants_vep.tsv \
    --vcf data/input/demo_variants.vcf \
    --out outputs/report_manual.md
```

---

## Ошибка: homework-check не проходит

**Симптом:**
```
[FAIL] variants list: file not found
```

**Причина:** Файл варiantов ещё не сгенерирован.

**Решение:**
```bash
export GITHUB_USER=<ваш_github>
python3 scripts/90_assign_variants.py \
    --vep_tsv data/input/demo_variants_vep.tsv \
    --sv_vcf  data/input/demo_sv.vcf \
    --out     homework/submissions/${GITHUB_USER}/variants_${GITHUB_USER}.txt
```

---

## Makefile: `make` завершается с ошибкой

**Симптом:**
```
Makefile:XX: *** missing separator.  Stop.
```

**Причина:** Makefile требует символы табуляции (Tab), не пробелы.

**Проверка:**
```bash
cat -A Makefile | grep "^\^I"  # ^I = Tab
```

**Решение:** При редактировании Makefile убедитесь, что ваш редактор вставляет Tab, а не пробелы.

---

## Общая диагностика

```bash
# Проверка всей системы
make check

# Просмотр последней ошибки
bash scripts/10_run_vep.sh 2>&1 | tail -30

# Проверка синтаксиса bash-скриптов
bash -n scripts/10_run_vep.sh && echo "OK"
bash -n scripts/20_prioritize_variants.sh && echo "OK"
```
