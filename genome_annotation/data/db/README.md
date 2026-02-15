# Директория баз данных

В этой директории хранятся базы данных, используемые на занятии.
Файлы здесь **не отслеживаются git** (слишком большие).

## Содержимое после подготовки преподавателем

| Файл                          | Создаётся скриптом          | Описание                                   |
|-------------------------------|-----------------------------|--------------------------------------------|
| `teaching_proteins.fasta`     | `01_get_data.sh`            | Небольшой набор белков (~100 тыс. последовательностей) |
| `teaching_proteins.dmnd`      | `02_make_diamond_db.sh`     | Индекс DIAMOND для набора белков           |
| `pfam_subset.hmm`             | `03_get_pfam.sh --subset`   | ~200 HMM-профилей для быстрого hmmscan     |
| `pfam_subset.hmm.h3{m,i,f,p}` | `03_get_pfam.sh --subset`  | Индексные файлы hmmpress                   |
| `Pfam-A.hmm` (опционально)   | `03_get_pfam.sh --full`     | Полная Pfam-A (~20 тыс. профилей)          |

## Как заполнить

Запустите скрипты подготовки из директории практикума (`genome_annotation/`):

```bash
bash scripts/01_get_data.sh
bash scripts/02_make_diamond_db.sh
bash scripts/03_get_pfam.sh --subset   # быстрый режим для занятия
# bash scripts/03_get_pfam.sh --full   # полная Pfam (опционально, медленнее)
```
