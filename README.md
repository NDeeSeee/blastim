# Курс NGS-биоинформатики на HPC

Набор практикумов по анализу данных секвенирования нового поколения (NGS) на высокопроизводительных вычислительных системах (HPC).

## Структура курса

| Лекция | Тема | Статус |
|--------|------|--------|
| [genome_annotation/](genome_annotation/) | Аннотация прокариотического генома | Готово |
| [variant_calling/](variant_calling/) | Поиск геномных вариантов (SNP/InDel) | Готово |
| [variant_interpretation/](variant_interpretation/) | Интерпретация вариантов (VEP, ClinVar, SV) | Готово |

## Структура репозитория

```
blastim_ngs/
├── README.md                  ← Вы здесь
├── LICENSE
├── shared/                    ← Общие ресурсы курса
│   ├── README.md
│   ├── COURSE_STRUCTURE.md    ← Обзор курса и практикумов
│   ├── HPC_SETUP.md           ← Настройка HPC (для преподавателей)
│   └── TROUBLESHOOTING_COMMON.md ← Общие проблемы и решения
├── genome_annotation/         ← Практикум: аннотация генома
│   ├── README.md
│   ├── INSTRUCTOR.md
│   ├── CONTRIBUTING.md
│   ├── Makefile
│   ├── env/
│   ├── scripts/
│   ├── workflow/
│   ├── data/
│   ├── homework/
│   └── docs/
├── variant_calling/           ← Практикум: поиск вариантов
│   ├── README.md
│   ├── INSTRUCTOR.md
│   ├── CONTRIBUTING.md
│   ├── Makefile
│   ├── env/
│   ├── scripts/
│   ├── workflow/
│   ├── data/
│   ├── homework/
│   └── docs/
└── variant_interpretation/    ← Практикум: интерпретация вариантов
    ├── README.md
    ├── INSTRUCTOR.md
    ├── CONTRIBUTING.md
    ├── Makefile
    ├── env/
    ├── scripts/
    ├── workflow/
    ├── data/
    ├── homework/
    └── docs/
```

## Быстрый старт

```bash
# Клонируйте репозиторий
git clone git@github.com:NDeeSeee/blastim_ngs.git
cd blastim_ngs

# Перейдите в нужный практикум
cd genome_annotation

# Следуйте инструкциям в README.md выбранного практикума
```

## Предварительные требования

- SSH-доступ к HPC
- Базовые навыки командной строки Linux
- Conda / Micromamba для управления окружениями

Подробные требования описаны в README каждого практикума.

## Для преподавателей и администраторов

- **[shared/HPC_SETUP.md](shared/HPC_SETUP.md)** — Подробное руководство по настройке HPC-системы
- **[shared/COURSE_STRUCTURE.md](shared/COURSE_STRUCTURE.md)** — Обзор структуры курса
- **[shared/TROUBLESHOOTING_COMMON.md](shared/TROUBLESHOOTING_COMMON.md)** — Решения частых проблем

Для каждого практикума см. файл `INSTRUCTOR.md` в соответствующей папке.

## Лицензия

MIT. См. [LICENSE](LICENSE).
