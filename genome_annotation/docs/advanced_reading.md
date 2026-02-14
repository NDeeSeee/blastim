# Углублённое чтение: аннотация прокариотического генома

Подборка сильных источников (статьи, книги, туториалы), которые дают глубину: что такое структурная и функциональная аннотация, почему поиск ORF сложнее "найти старт/стоп", как правильно мыслить про гомологию и домены, и как это превращается в воспроизводимую практику на HPC.

Источники сгруппированы по смысловым слоям — можно пройти как мини-курс для самоподготовки.

---

## A) Большая картина: что такое аннотация, где ломается, как её проверять

1. **Обзор: structural + functional annotation + QC, ре-аннотация и будущее**
   Ejigu & Jung, 2020 (review, open access). Очень хороший "каркас мышления" для всей темы.
   [Biology 9(9):295](https://doi.org/10.3390/biology9090295)

2. **Как делает "индустриальный стандарт" NCBI (PGAP): процесс и evidence**
   - [Общий обзор PGAP](https://www.ncbi.nlm.nih.gov/refseq/annotation_prok/) — что это и зачем
   - [Описание процесса](https://www.ncbi.nlm.nih.gov/refseq/annotation_prok/process/) — как комбинируют ORF + HMM + curated proteins
   - [Evidence types](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/) — HMMs, BlastRules, domain architectures. Идеальная философия аннотации: как строится "лестница доказательств"
   - Tatusova et al., 2016 — академическая версия описания PGAP.
     [Nucleic Acids Res. 44(D1)](https://doi.org/10.1093/nar/gkw569)
   - Haft et al., 2018 — про качество и эволюцию подхода RefSeq/аннотации прокариот.
     [Nucleic Acids Res. 46(D1)](https://doi.org/10.1093/nar/gkx1068)

> **Рекомендация:** прочитать Ejigu + PGAP evidence и выстроить понимание как "evidence ladder": ORF -> гомология -> домены -> curated rules -> sanity checks.

---

## B) Структурная аннотация: поиск белок-кодирующих генов

3. **Prodigal** (Hyatt et al., 2010) — классика прокариотического gene calling. Объясняет, почему алгоритм предсказания генов — это не тривиальный поиск стартов/стопов.
   [BMC Bioinformatics 11:119](https://doi.org/10.1186/1471-2105-11-119)

4. **Сравнение инструментов предсказания генов** (DiMonaco et al., 2022) — полезно для понимания ограничений отдельных инструментов и новых ML-подходов.
   [Bioinformatics 38(5)](https://doi.org/10.1093/bioinformatics/btab827)

---

## C) Функциональная аннотация по гомологии: BLAST-мышление

Статистика важнее названий. "Best hit" != функция.

5. **BLAST** (Altschul et al., 1990) — "конституция" локального выравнивания и смысл E-value.
   [J. Mol. Biol. 215(3)](https://doi.org/10.1016/S0022-2836(05)80360-2)

6. **DIAMOND** (Buchfink et al., 2015) — современный быстрый BLAST-like инструмент для белков, практический must-know.
   [Nature Methods 12](https://doi.org/10.1038/nmeth.3176)

> **Рекомендация:** разобрать триаду identity + coverage + E-value, и отдельно проработать проблему паралогов и мультидоменных белков.

---

## D) Консервативные домены: почему HMM — это "более честная гомология"

7. **HMMER3** (Eddy, 2011) — что такое profile HMM и как они стали быстрыми. Очень сильная статья.
   [PLOS Comput. Biol. 7(10)](https://doi.org/10.1371/journal.pcbi.1002195)

8. **Pfam** (Mistry et al., 2021) — что именно за домены, как устроена база, почему это curated ресурс.
   [Nucleic Acids Res. 49(D1)](https://doi.org/10.1093/nar/gkaa913)
   Исторический фундамент: Finn et al., 2014.
   [Nucleic Acids Res. 42(D1)](https://doi.org/10.1093/nar/gkt1223)

---

## E) Практика: пошаговые туториалы

9. **Prokka** (Seemann, 2014) — коротко и по делу, с ориентиром по времени выполнения.
    [Bioinformatics 30(14)](https://doi.org/10.1093/bioinformatics/btu153)

10. **Prokka GitHub** — актуальные параметры, заметки и важная ремарка про Bakta как next-gen.
    [github.com/tseemann/prokka](https://github.com/tseemann/prokka)

11. **Metagenomics workshop: Prokka walkthrough** — учебный текст с понятным объяснением, что происходит внутри.
    [metagenomics-workshop.readthedocs.io](https://metagenomics-workshop.readthedocs.io/en/2014-11-uppsala/functional-annotation/prokka.html)

12. **Galaxy Training: Prokka tutorial** — педагогически хорошо сделанный материал.
    [training.galaxyproject.org](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/annotation-with-prokka/tutorial.html)

---

## F) Современная замена Prokka

13. **Bakta** (Schwengers et al., 2021) — современный аннотатор для бактерий, плазмид и MAGs. Более стандартизованный output по сравнению с Prokka.
    [Microb. Genom. 7(11)](https://doi.org/10.1099/mgen.0.000685)
    GitHub: [github.com/oschwengers/bakta](https://github.com/oschwengers/bakta)
    Bakta Web (2025) — веб-интерфейс с низким барьером входа: [Nucleic Acids Res. (2025)](https://doi.org/10.1093/nar/gkaf335)

> **Рекомендация:** Prokka — как основной инструмент (легко объяснить, быстро запустить). Bakta — как "что часто используют сейчас, если нужен более стандартизованный output".

---

## G) Богатая функциональная аннотация (GO/EC/paths)

Для углублённого изучения — не обязательно запускать на занятии, но полезно понимать.

14. **InterProScan 5** (Jones et al., 2014) — "тяжёлая артиллерия" доменных и сигнатурных баз.
    [Bioinformatics 30(9)](https://doi.org/10.1093/bioinformatics/btu031)

15. **eggNOG-mapper v2** (Cantalapiedra et al., 2021) — мост между "гомология" и "ортология -> функция", объясняет, почему ортология важна.
    [Mol. Biol. Evol. 38(12)](https://doi.org/10.1093/molbev/msab293)
    Оригинальная версия: Huerta-Cepas et al., 2017.
    [Mol. Biol. Evol. 34(8)](https://doi.org/10.1093/molbev/msx148)

---

## Минимальный маршрут на 3-5 вечеров

| Вечер | Что читать | Зачем |
|-------|-----------|-------|
| 1 | Ejigu review + PGAP evidence | Понимание системы доказательств |
| 2 | Prodigal paper | Понимание gene calling |
| 3 | BLAST + DIAMOND | Гомология и статистика |
| 4 | HMMER + Pfam | Домены как уровень выше |
| 5 | Prokka tutorial + запуск самому | Практика |
