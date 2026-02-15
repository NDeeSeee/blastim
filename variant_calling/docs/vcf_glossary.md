# Справочник полей VCF

Описание полей формата VCF (Variant Call Format) с фокусом на freebayes.

## Фиксированные поля (обязательные, 8 столбцов)

| Поле | Описание | Пример |
|------|----------|--------|
| **CHROM** | Хромосома / контиг | `NC_000964.3` |
| **POS** | Позиция (1-based) | `12345` |
| **ID** | Идентификатор варианта (часто `.`) | `.` |
| **REF** | Референсный аллель | `A` |
| **ALT** | Альтернативный аллель | `G` |
| **QUAL** | Оценка качества (Phred-scaled) | `1234.56` |
| **FILTER** | Статус фильтра (`.` или `PASS` = прошёл) | `.` |
| **INFO** | Дополнительная информация (ключ=значение) | `DP=45;AO=40;...` |

### QUAL — Phred-scaled quality score

QUAL = -10 * log10(P(вариант ложный)).

| QUAL | Вероятность ошибки | Надёжность |
|------|-------------------|------------|
| 10 | 10% | Низкая |
| 20 | 1% | Средняя |
| 30 | 0.1% | Хорошая |
| 100 | 10^-10 | Высокая |

## Поля INFO (freebayes-специфичные)

### Ключевые поля

| Поле | Описание | Интерпретация |
|------|----------|---------------|
| **DP** | Total read depth | Общая глубина покрытия в позиции |
| **RO** | Reference Observation count | Число прочтений, подтверждающих REF |
| **AO** | Alternate Observation count | Число прочтений, подтверждающих ALT |
| **AF** | Allele Frequency (estimated) | Оценка частоты ALT аллеля (AO/DP) |
| **TYPE** | Тип варианта | `snp`, `ins`, `del`, `mnp`, `complex` |

### Поля поддержки прочтений (strand/position bias)

| Поле | Описание |
|------|----------|
| **SRF** | Number of reference observations on the forward strand |
| **SRR** | Number of reference observations on the reverse strand |
| **SAF** | Number of alternate observations on the forward strand |
| **SAR** | Number of alternate observations on the reverse strand |
| **SRP** | Strand balance probability for REF (Phred-scaled) |
| **SAP** | Strand balance probability for ALT (Phred-scaled) |
| **RPP** | Read Placement Probability (position bias for ALT) |
| **RPL** | Reads Placed Left (ALT reads placed to the left of variant) |
| **RPR** | Reads Placed Right (ALT reads placed to the right of variant) |

### Интерпретация SRP/SAP

Если SAP > 20 (т.е. P < 0.01), прочтения с ALT аллелем сильно смещены на одну цепь — возможный артефакт секвенирования.

### Поля качества маппинга

| Поле | Описание |
|------|----------|
| **MQM** | Mean mapping quality of observed alternate alleles |
| **MQMR** | Mean mapping quality of observed reference alleles |
| **PAIRED** | Proportion of ALT observations from properly paired reads |
| **PAIREDR** | Proportion of REF observations from properly paired reads |

## Поля FORMAT (per-sample, столбец 9+)

| Поле | Описание | Пример |
|------|----------|--------|
| **GT** | Genotype | `0/1` (het), `1/1` (hom ALT), `0/0` (hom REF) |
| **DP** | Sample depth | `45` |
| **AD** | Allelic depths (REF, ALT) | `5,40` |
| **RO** | Reference allele observation count | `5` |
| **QR** | Sum of quality of reference observations | `180` |
| **AO** | Alternate allele observation count | `40` |
| **QA** | Sum of quality of alternate observations | `1440` |
| **GL** | Genotype likelihoods (log10-scaled) | `-120.5,0,-2.3` |

### Интерпретация GT для гаплоидных организмов (ploidy=1)

| GT | Значение |
|----|----------|
| `0` | Референс |
| `1` | Альтернативный аллель |

Для диплоидов (ploidy=2):

| GT | Значение |
|----|----------|
| `0/0` | Гомозигота по референсу |
| `0/1` | Гетерозигота |
| `1/1` | Гомозигота по ALT |

## Пример строки VCF

```
NC_000964.3  12345  .  A  G  1234.56  .  AB=0;AO=40;DP=45;RO=5;TYPE=snp;...  GT:DP:AD:RO:QR:AO:QA:GL  1:45:5,40:5:180:40:1440:-120.5,0
```

Разбор:
- **Позиция**: 12345 на хромосоме NC_000964.3
- **Замена**: A → G (SNP)
- **Качество**: QUAL=1234.56 (очень надёжный)
- **Покрытие**: DP=45 (45 прочтений)
- **Поддержка**: 40 прочтений за ALT, 5 за REF → AF ≈ 0.89
- **Генотип**: `1` (гаплоид, ALT аллель)
