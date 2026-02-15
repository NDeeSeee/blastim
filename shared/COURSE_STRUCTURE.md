# Структура курса: NGS-анализ на HPC

Этот курс состоит из двух независимых практикумов, которые преподаются отдельно студентам программы по Next-Generation Sequencing.

## Обзор практикумов

| Практикум | Язык | Объём | Ключевые темы | Версия |
|-----------|------|-------|---------------|--------|
| **1. Аннотация прокариотического генома** | Russian | 2 часа | Структурная аннотация (Prokka), функциональная аннотация по гомологии (DIAMOND), поиск доменов (HMMER/Pfam), "лестница доказательств" | `genome_annotation/` |
| **2. Поиск геномных вариантов (Variant Calling)** | Russian | 2 часа | Выравнивание (BWA), поиск вариантов (freebayes), фильтрация (bcftools), оценка качества вариантов, контроль качества | `variant_calling/` |

## Организм

Оба практикума используют **Bacillus subtilis subsp. subtilis str. 168** (GenoBase reference):
- Размер генома: ~4.2 Мб
- GC-содержание: ~43%
- Число предсказанных белков: ~4200
- Статус: хорошо изучен, идеален для обучения

Генетический материал (референсный геном, базы данных Swiss-Prot/Pfam) загружаются каждым модулем автоматически при первом запуске (`01_get_data.sh`).

## Целевая аудитория

- **genome_annotation**: Студенты, изучающие структурную и функциональную аннотацию геномов
- **variant_calling**: Студенты, изучающие вариационный анализ в популяционной генетике и медицинской геномике

## Структура репозитория

```
blastim_ngs/
├── README.md                    ← Вы здесь
├── LICENSE
├── .gitignore
├── shared/                      ← Общие ресурсы
│   ├── README.md
│   ├── COURSE_STRUCTURE.md      ← Этот файл
│   ├── HPC_SETUP.md
│   └── TROUBLESHOOTING_COMMON.md
├── genome_annotation/           ← Практикум 1: Аннотация
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
└── variant_calling/             ← Практикум 2: Variant Calling
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

## Как начать

### Студентам

1. Выберите интересующий вас практикум:
   ```bash
   cd genome_annotation          # Или: cd variant_calling
   ```

2. Следуйте инструкциям в соответствующем **README.md**

3. Ознакомьтесь с материалами в **docs/** папке

### Преподавателям

1. Прочитайте **shared/HPC_SETUP.md** для подготовки HPC-системы

2. Для каждого практикума см. **INSTRUCTOR.md** в соответствующей папке

3. Консультируйтесь с **shared/TROUBLESHOOTING_COMMON.md** при возникновении проблем

## Домашние задания

- **HW1 (genome_annotation)**: Защита 3 генов с использованием "лестницы доказательств"
  - Читай: `genome_annotation/docs/homework_01.md`
  - Сдача: через fork/PR, см. `genome_annotation/CONTRIBUTING.md`

- **HW2 (variant_calling)**: Защита 5 вариантов с анализом качества
  - Читай: `variant_calling/docs/homework_02.md`
  - Сдача: через fork/PR, см. `variant_calling/CONTRIBUTING.md`

## Требования к системе

Оба практикума требуют:
- SSH-доступ к HPC или локальная система Linux/macOS
- Conda / Micromamba
- ~2-4 ГБ свободного места на диске (на модуль)
- Базовые знания Linux командной строки

Детали см. в **shared/HPC_SETUP.md**.

## Лицензия

MIT. См. [LICENSE](../LICENSE).

Данные и базы (Swiss-Prot, Pfam, B. subtilis 168 геном) загружаются из официальных источников и подчиняются их собственным лицензиям.
