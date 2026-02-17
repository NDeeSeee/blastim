# Рекомендуемые материалы

## Ключевые статьи

### Инструменты аннотации

1. **McLaren et al., 2016** — Ensembl Variant Effect Predictor
   *Genome Biology* 17, 122
   doi: 10.1186/s13059-016-0974-4
   → Основная статья VEP. Читать в первую очередь.

2. **Wang et al., 2010** — ANNOVAR: functional annotation of genetic variants
   *Nucleic Acids Research* 38(16): e164
   doi: 10.1093/nar/gkq603

3. **Cingolani et al., 2012** — SnpEff: Variants effect prediction
   *Fly* 6(2): 80–92
   doi: 10.4161/fly.19695

### Классификация вариантов

4. **Richards et al., 2015** — ACMG/AMP Standards and Guidelines for variant classification
   *Genetics in Medicine* 17(5): 405–424
   doi: 10.1038/gim.2015.30
   → Стандарт клинической классификации. Обязательно.

5. **Harrison et al., 2019** — Using ClinVar as a resource for variant interpretation
   *Current Protocols in Human Genetics*
   doi: 10.1002/cphg.93

### Популяционные частоты

6. **Karczewski et al., 2020** — gnomAD v2.1.1: The mutational constraint spectrum
   *Nature* 581: 434–443
   doi: 10.1038/s41586-020-2308-7
   → Описание gnomAD v2. Ключевая база данных популяционных частот.

7. **Chen et al., 2024** — gnomAD v4.0: 730,947 genomes
   *Nature* 634: 1–10
   doi: 10.1038/s41586-024-07450-7

### Предсказание патогенности

8. **Rentzsch et al., 2019** — CADD: predicting the deleteriousness of variants
   *Nucleic Acids Research* 47: D886–D894
   doi: 10.1093/nar/gky1016

9. **Cheng et al., 2023** — Accurate proteome-wide missense variant effect prediction
   *Science* 381: eadg7492 (AlphaMissense)
   doi: 10.1126/science.adg7492
   → Революционный deep learning метод от DeepMind.

10. **Jaganathan et al., 2019** — Predicting splicing from primary sequence with deep learning
    *Cell* 176: 535–548 (SpliceAI)
    doi: 10.1016/j.cell.2018.12.015

### Структурные варианты

11. **Chen et al., 2016** — Manta: rapid detection of SVs
    *Bioinformatics* 32(8): 1220–1222
    doi: 10.1093/bioinformatics/btv710

12. **Riggs et al., 2020** — Technical standards for CNV analysis
    *Genetics in Medicine* 22: 245–257
    doi: 10.1038/s41436-019-0686-8

---

## Обзорные статьи (рекомендуется для общего понимания)

13. **Biesecker & Green, 2014** — Diagnostic clinical genome and exome sequencing
    *New England Journal of Medicine* 370: 2418–2425
    doi: 10.1056/NEJMra1312543

14. **Claussnitzer et al., 2020** — A brief history of human disease genetics
    *Nature* 577: 179–189
    doi: 10.1038/s41586-019-1879-7

15. **Turnbull et al., 2018** — 100,000 Genomes Project: rare disease and cancer
    *American Journal of Human Genetics* 103: 684–693
    doi: 10.1016/j.ajhg.2018.08.016

---

## Онлайн-ресурсы

| Ресурс | Описание | URL |
|--------|----------|-----|
| VEP Documentation | Полная документация по VEP флагам | ensembl.org/vep |
| gnomAD Browser | Интерактивный просмотр популяционных частот | gnomad.broadinstitute.org |
| ClinVar | Поиск клинических интерпретаций | ncbi.nlm.nih.gov/clinvar |
| ClinGen | Клиническая значимость генов | clinicalgenome.org |
| Franklin (Genoox) | AI-ассистент ACMG классификации | franklin.genoox.com |
| VariantValidator | Валидация HGVS нотации | variantvalidator.org |
| MutPred2 | Предсказание патогенности missense | mutpred.mutdb.org |
| OMIM | Менделирующие заболевания | omim.org |
| Orphanet | Редкие болезни | orpha.net |

---

## Курсы и обучение

| Курс | Платформа |
|------|-----------|
| Bioinformatics for Biologists (Wellcome Connecting Science) | Coursera |
| Genomic Data Science Specialization | Coursera / JHU |
| Introduction to Clinical Genomics | edX |
| EMBL-EBI Variant Analysis courses | embl-ebi.ac.uk/training |

---

## Базы данных — краткая справка

| База | Тип данных | Релиз | Главная ценность |
|------|-----------|-------|-----------------|
| **gnomAD v4** | ~730k геномов | 2024 | Популяционные частоты |
| **ClinVar** | Клинические | ежемесячно | Курируемые интерпретации |
| **dbSNP** | Все варианты | Обновляется | rs ID идентификация |
| **COSMIC** | Соматические | квартально | Рак-ассоциированные мутации |
| **ClinGen** | Гены + варианты | Обновляется | Доказательная база патогенности |
| **HGMD** | Патогенные | квартально | Наследственные болезни |
| **LOVD** | Ген-специфичные | Обновляется | Специфика по гену/болезни |
