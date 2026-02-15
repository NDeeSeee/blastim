# Как сдать домашнее задание (fork / PR)

Краткая инструкция по отправке ДЗ через Pull Request.

## 1. Fork и клонирование

```bash
# Сделайте fork через GitHub UI (кнопка "Fork" на странице репозитория)

# Клонируйте свой fork
git clone git@github.com:<ваш_username>/blastim_ngs.git
cd blastim_ngs/variant_calling
```

## 2. Создание ветки

Формат имени ветки: `hw<номер>-<github_username>`

```bash
git checkout -b hw2-myusername
```

## 3. Выполнение задания

Следуйте инструкциям в `docs/homework_02.md`. Все файлы размещайте в своей папке:

```
homework/submissions/<github_username>/
├── variants_<username>.txt       ← сгенерирован скриптом
└── defense_<username>.md         ← заполняется студентом
```

## 4. Коммит и пуш

```bash
git add homework/submissions/<github_username>/
git commit -m "HW2: <github_username>"
git push origin hw2-<github_username>
```

Коммитьте **только свою папку** в `homework/submissions/`. Не изменяйте скрипты, шаблоны и другие файлы репозитория.

## 5. Создание Pull Request

1. Откройте свой fork на GitHub.
2. Нажмите **"Compare & pull request"**.
3. Заголовок PR: `HW2: <github_username>`
4. Убедитесь, что PR направлен в основной репозиторий преподавателя (`main` ветка).
5. Нажмите **"Create pull request"**.

## Самопроверка перед отправкой

```bash
python scripts/91_homework_check.py --user $GITHUB_USER
```

Скрипт проверит наличие файлов и базовую структуру.

## Частые ошибки

- **Не создали ветку** — коммиты попали в `main`. Создайте ветку и перенесите коммит: `git checkout -b hw2-username && git push origin hw2-username`.
- **Изменили чужие файлы** — в PR должны быть только файлы из вашей папки `homework/submissions/<username>/`.
- **Забыли файл вариантов** — сначала запустите `90_homework_assign_variants.py`, потом заполняйте шаблон.
