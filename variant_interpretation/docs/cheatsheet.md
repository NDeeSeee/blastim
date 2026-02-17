# Шпаргалка: Интерпретация вариантов

## VCF — быстрый просмотр

```bash
# Варианты без заголовка (только данные)
grep -v "^#" data/input/demo_variants.vcf

# Заголовки INFO (что означают поля)
grep "^##INFO" data/input/demo_variants.vcf

# Только конкретные поля через awk
grep -v "^#" demo.vcf | awk '{print $1,$2,$3,$4,$5,$6}' | column -t

# Подсчёт вариантов
grep -vc "^#" data/input/demo_variants.vcf
```

---

## bcftools — просмотр и фильтрация

```bash
# Показать все варианты
bcftools view data/input/demo_variants.vcf

# Статистика
bcftools stats data/input/demo_variants.vcf | grep "^SN"

# Фильтрация по QUAL и INFO/DP
bcftools view -i 'QUAL>30 && INFO/DP>10' demo.vcf

# Фильтрация по INFO/gnomAD_AF (редкие варианты)
bcftools view -i 'INFO/gnomAD_AF<0.01 || INFO/gnomAD_AF="."' demo.vcf

# Извлечение конкретных полей в TSV
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/gnomAD_AF\n' demo.vcf

# Фильтрация по нескольким критериям
bcftools view -i 'QUAL>20 && INFO/DP>=10 && INFO/AF>=0.2' demo.vcf

# Выбор вариантов конкретного типа (SNP)
bcftools view --type snps demo.vcf

# Выбор вариантов в регионе
bcftools view demo.vcf chr17:7665000-7680000
```

---

## VEP TSV — просмотр

```bash
# Показать заголовок
grep "^#" data/input/demo_variants_vep.tsv | tail -1

# Данные без заголовков (## и #)
grep -v "^#" data/input/demo_variants_vep.tsv

# Только HIGH-impact варианты (поле Extra содержит IMPACT=HIGH)
grep "IMPACT=HIGH" data/input/demo_variants_vep.tsv

# Варианты с ClinVar Pathogenic
grep "ClinVar_CLNSIG=Pathogenic" data/input/demo_variants_vep.tsv

# Подсчёт по последствию (consequence)
grep -v "^#" data/input/demo_variants_vep.tsv | \
    awk '{print $7}' | sort | uniq -c | sort -rn

# Только canonical транскрипты
grep "CANONICAL=YES" data/input/demo_variants_vep.tsv
```

---

## Python — работа с VEP TSV и VCF

```python
import pandas as pd

# Читаем VEP TSV (пропускаем строки ##)
df = pd.read_csv("data/input/demo_variants_vep.tsv",
                 sep="\t", comment="#",
                 header=0)
# Переименовываем первую колонку (убираем #)
df.columns = [c.lstrip("#") for c in df.columns]

# Расширяем поле Extra
def parse_extra(extra: str) -> dict:
    d = {}
    for token in str(extra).split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            d[k] = v
    return d

extras = df["Extra"].apply(parse_extra).apply(pd.Series)
df = pd.concat([df, extras], axis=1)

# Фильтрация HIGH-impact
high = df[df["IMPACT"] == "HIGH"]

# Топ по gnomAD_AF
df["gnomAD_num"] = pd.to_numeric(df["gnomADe_AF"], errors="coerce")
rare = df[df["gnomAD_num"] < 0.01].sort_values("gnomAD_num")
```

---

## SV VCF — структурные варианты

```bash
# Показать SV варианты
grep -v "^#" data/input/demo_sv.vcf | \
    awk '{print $3,$1,$2,$4,$5}' | column -t

# Извлечь тип SV
grep -v "^#" data/input/demo_sv.vcf | \
    grep -oP "SVTYPE=[^;]+" | sort | uniq -c

# Размер структурных вариантов
grep -v "^#" data/input/demo_sv.vcf | \
    grep -oP "SVLEN=-?\d+" | \
    awk -F= '{print ($2<0 ? -$2 : $2)" bp"}'

# Только PASS SV
grep -v "^#" data/input/demo_sv.vcf | awk '$7=="PASS"'
```

---

## Приоритизация — скрипт вручную

```bash
# Запустить шаг приоритизации
bash scripts/20_prioritize_variants.sh

# Просмотреть результат
column -t outputs/run_*/prioritized.tsv | head -20

# Топ-5 вариантов
sort -k2 -rn outputs/run_*/prioritized.tsv | head -6 | column -t
```

---

## VEP через REST API (вручную, curl)

```bash
# Аннотировать один вариант
curl -X POST "https://rest.ensembl.org/vep/human/region" \
  -H "Content-Type: application/json" \
  -H "Accept: application/json" \
  -d '{"regions":["17 7674220 7674220 C/T +"]}' | jq .

# С параметрами: symbol, canonical, hgvs
curl -X POST \
  "https://rest.ensembl.org/vep/human/region?symbol=1&canonical=1&hgvs=1&af_gnomade=1" \
  -H "Content-Type: application/json" \
  -d '{"regions":["17 7674220 7674220 C/T +"]}' | \
  jq '.[0].transcript_consequences[] | select(.canonical==1) | {gene_symbol, consequence_terms, impact, hgvsc, hgvsp}'
```

---

## HGVSc / HGVSp — расшифровка

| HGVS | Расшифровка |
|------|-------------|
| `c.524G>A` | кодирующая позиция 524: G→A |
| `c.5266dupC` | дупликация C после позиции 5266 |
| `c.1100delC` | делеция C в позиции 1100 |
| `p.Arg175His` | аргинин в позиции 175 → гистидин |
| `p.Lys745Glu` | лизин 745 → глутаминовая кислота |
| `p.Val600Glu` | валин 600 → глутаминовая кислота (BRAF V600E) |
| `p.Glu55*` | глутамат 55 → стоп-кодон (нонсенс) |
| `p.Thr790Met` | треонин 790 → метионин (EGFR T790M) |

---

## Полезные ссылки

| Ресурс | URL |
|--------|-----|
| ClinVar | https://www.ncbi.nlm.nih.gov/clinvar/ |
| gnomAD | https://gnomad.broadinstitute.org/ |
| Ensembl VEP | https://www.ensembl.org/vep |
| VEP REST API | https://rest.ensembl.org/documentation/info/vep_region_post |
| OMIM | https://www.omim.org/ |
| UniProt | https://www.uniprot.org/ |
| Franklin (ACMG) | https://franklin.genoox.com/ |
