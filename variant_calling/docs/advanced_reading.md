# Углублённое чтение: поиск геномных вариантов

Подборка сильных источников (статьи, туториалы), которые дают глубину: что такое variant calling, как устроены алгоритмы, какие артефакты бывают, и как мыслить про качество вариантов.

Источники сгруппированы по смысловым слоям — можно пройти как мини-курс для самоподготовки.

---

## A) Большая картина: что такое variant calling, где ломается

1. **Обзор: от прочтений до вариантов**
   Nielsen et al., 2011 — классический обзор подходов к variant calling на NGS-данных.
   [Nature Reviews Genetics 12](https://doi.org/10.1038/nrg2986)

2. **GATK Best Practices** — индустриальный стандарт для пайплайнов variant calling.
   - [Обзор](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
   - Van der Auwera & O'Connor, 2020 — книга "Genomics in the Cloud" (O'Reilly).

> **Рекомендация:** начать с Nielsen + GATK Best Practices для понимания, зачем фильтруем и что именно может пойти не так.

---

## B) Выравнивание: от FASTQ до BAM

3. **BWA** (Li & Durbin, 2009) — стандарт для выравнивания коротких прочтений.
   [Bioinformatics 25(14)](https://doi.org/10.1093/bioinformatics/btp324)

4. **SAM/BAM формат** (Li et al., 2009) — спецификация формата SAM, описание флагов и тегов.
   [Bioinformatics 25(16)](https://doi.org/10.1093/bioinformatics/btp352)

---

## C) Variant calling: алгоритмы

5. **freebayes** (Garrison & Marth, 2012) — Bayesian haplotype-based variant caller.
   Основа нашего пайплайна. Объясняет философию "haplotype-based" подхода.
   [arXiv:1207.3907](https://arxiv.org/abs/1207.3907)

6. **bcftools mpileup/call** (Danecek et al., 2021) — альтернативный подход на основе pileup.
   [GigaScience 10(2)](https://doi.org/10.1093/gigascience/giab008)

7. **Сравнение caller'ов** — полезно для понимания, почему разные инструменты дают разные результаты.
   Supernat et al., 2018.
   [PeerJ 6:e5737](https://doi.org/10.7717/peerj.5737)

> **Рекомендация:** прочитать Garrison & Marth для понимания freebayes, затем Danecek et al. для альтернативного подхода.

---

## D) Формат VCF и инструменты

8. **VCF спецификация** (Danecek et al., 2011) — формальное описание формата VCF.
   [Bioinformatics 27(15)](https://doi.org/10.1093/bioinformatics/btr330)

9. **BCFtools** (Li, 2011) — основной инструментарий для работы с VCF/BCF.
   [Bioinformatics 27(21)](https://doi.org/10.1093/bioinformatics/btr509)

10. **htslib** — библиотека для работы с форматами HTS (BAM, CRAM, VCF, BCF).
    [github.com/samtools/htslib](https://github.com/samtools/htslib)

---

## E) Фильтрация и качество вариантов

11. **Что такое QUAL, GQ и как их интерпретировать**
    GATK Technical Documentation:
    - [Understanding QUAL](https://gatk.broadinstitute.org/hc/en-us/articles/360035531572)
    - [Hard-filtering germline variants](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471)

12. **Strand bias, position bias и другие артефакты**
    Guo et al., 2012 — систематический обзор артефактов при variant calling.
    [BMC Genomics 13:666](https://doi.org/10.1186/1471-2164-13-666)

---

## F) Variant calling в бактериальной геномике

13. **Snippy** — популярный пайплайн для бактериального variant calling (обёртка над freebayes).
    [github.com/tseemann/snippy](https://github.com/tseemann/snippy)

14. **Bacterial WGS analysis guide** (Kwong et al., 2015) — обзор WGS-анализа в медицинской микробиологии.
    [Clinical Microbiology Reviews 28(3)](https://doi.org/10.1128/CMR.00056-14)

---

## G) Практические туториалы

15. **Galaxy Training: Variant Calling** — пошаговый туториал с визуализацией.
    [training.galaxyproject.org](https://training.galaxyproject.org/training-material/topics/variant-analysis/)

16. **Data Carpentry: Variant Calling Workflow** — отличный учебный материал для начинающих.
    [datacarpentry.org](https://datacarpentry.org/wrangling-genomics/)

---

## Минимальный маршрут на 3-5 вечеров

| Вечер | Что читать | Зачем |
|-------|-----------|-------|
| 1 | Nielsen review + GATK Best Practices | Общая картина и стандарты |
| 2 | Li & Durbin (BWA) + SAM spec | Понимание выравнивания |
| 3 | Garrison & Marth (freebayes) | Алгоритм variant calling |
| 4 | Danecek (VCF spec) + BCFtools | Формат данных и инструменты |
| 5 | GATK hard-filtering + Galaxy tutorial | Фильтрация и практика |
