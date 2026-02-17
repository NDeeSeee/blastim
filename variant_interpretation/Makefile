# ── Makefile для практикума: Интерпретация вариантов ─────────────────
SHELL := /bin/bash

.PHONY: env check get-data vep-cache vep prioritize sv-demo modern all clean clean-all homework-check help

help:
	@echo "Доступные цели:"
	@echo "  make env              Создать conda-окружение ngs-interpret"
	@echo "  make check            Проверить систему и зависимости"
	@echo "  make get-data         Проверить входные данные"
	@echo "  make vep-cache        Загрузить VEP cache (только для преподавателя)"
	@echo "  make vep              Запустить VEP аннотацию (шаг 10)"
	@echo "  make prioritize       Приоритизировать варианты (шаги 20-21)"
	@echo "  make sv-demo          SV/CNV мини-демо (шаг 30)"
	@echo "  make modern           Контекст современных инструментов (шаг 40)"
	@echo "  make all              Запустить полный пайплайн (шаги 10-30)"
	@echo "  make clean            Удалить outputs/ (сохраняет данные)"
	@echo "  make clean-all        Удалить outputs/ и VEP cache"
	@echo "  make homework-check   Проверить оформление ДЗ3"

env:
	bash env/micromamba_create.sh

check:
	bash scripts/00_check_system.sh

get-data:
	bash scripts/01_get_data.sh

vep-cache:
	bash scripts/02_get_vep_cache.sh

vep:
	bash scripts/10_run_vep.sh

prioritize:
	bash scripts/20_prioritize_variants.sh
	python3 scripts/21_make_report.py \
		--tsv "$$(ls -td outputs/run_* 2>/dev/null | head -1)/vep/vep_output.tsv" \
		--vcf data/input/demo_variants.vcf \
		--out "$$(ls -td outputs/run_* 2>/dev/null | head -1)/report.md"

sv-demo:
	bash scripts/30_sv_cnv_intro.sh

modern:
	bash scripts/40_modern_context.sh

all:
	bash workflow/run_all.sh

clean:
	rm -rf outputs/run_*
	@echo "Outputs удалены. Данные в data/ сохранены."

clean-all: clean
	@echo "Для удаления VEP cache выполните: rm -rf $${VEP_CACHE_DIR}"
	@echo "Путь VEP_CACHE_DIR задаётся в workflow/config.sh"

homework-check:
	@test -n "$$GITHUB_USER" || (echo "ERROR: export GITHUB_USER=<ваш_github>"; exit 1)
	python3 scripts/91_homework_check.py --user $$GITHUB_USER
