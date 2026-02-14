# Курс NGS-биоинформатики на HPC

Набор практикумов по анализу данных секвенирования нового поколения (NGS) на высокопроизводительных вычислительных системах (HPC).

## Структура курса

| Лекция | Тема | Статус |
|--------|------|--------|
| [genome_annotation/](genome_annotation/) | Аннотация прокариотического генома | Готово |
| [variant_calling/](variant_calling/) | Поиск геномных вариантов (SNP/InDel) | В разработке |

## Структура репозитория

```
blastim_ngs/
├── README.md                  ← Вы здесь
├── LICENSE
├── shared/                    ← Общие ресурсы курса
├── genome_annotation/         ← Практикум: аннотация генома
│   ├── README.md
│   ├── INSTRUCTOR.md
│   ├── Makefile
│   ├── env/
│   ├── scripts/
│   ├── workflow/
│   ├── data/
│   └── docs/
└── variant_calling/           ← Практикум: поиск вариантов
    └── README.md
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

## Лицензия

MIT. См. [LICENSE](LICENSE).
