#!/usr/bin/env python3
"""Generate a markdown report from VEP output and prioritized variants.

Usage:
    python3 scripts/21_make_report.py \
        --tsv outputs/<run>/vep/vep_output.tsv \
        --vcf data/input/demo_variants.vcf \
        --out outputs/<run>/report.md
"""

import argparse
from collections import Counter
from datetime import date
from pathlib import Path


# -- helpers ------------------------------------------------------------------

def parse_vep_tsv(tsv_path: Path) -> list[dict]:
    """Parse VEP TSV output.

    Skips lines starting with '##'.  The first line starting with a single
    '#' is treated as the column-header row (leading '#' stripped).
    Returns a list of dicts keyed by column name; the Extra field is expanded
    into individual key=value pairs.
    """
    records = []
    header = None
    with tsv_path.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("##"):
                continue
            if header is None and line.startswith("#"):
                header = line.lstrip("#").split("\t")
                continue
            if header is None:
                continue
            parts = line.split("\t")
            row = dict(zip(header, parts))
            # Expand semi-colon-delimited Extra column into individual keys
            extra = row.get("Extra", "")
            for item in extra.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    row[k] = v
            records.append(row)
    return records


def parse_vcf(vcf_path: Path) -> dict[str, dict]:
    """Parse a VCF file; return dict keyed by 'CHROM_POS_REF/ALT'.

    Extracts QUAL and INFO fields: DP, AF, gnomAD_AF, ClinVar.
    """
    variants: dict[str, dict] = {}
    with vcf_path.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _vid, ref, alt, qual = (
                parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
            )
            info_str = parts[7]
            key = f"{chrom}_{pos}_{ref}/{alt}"
            info: dict[str, str] = {}
            for item in info_str.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
            variants[key] = {
                "QUAL": qual,
                "DP": info.get("DP", "."),
                "AF": info.get("AF", "."),
                "gnomAD_AF": info.get("gnomAD_AF", "."),
                "ClinVar": info.get("ClinVar", "."),
            }
    return variants


def _to_float(value: str, default: float = 1.0) -> float:
    """Convert string to float; return *default* on failure."""
    try:
        return float(value)
    except (ValueError, TypeError):
        return default


def fmt_af(value: str) -> str:
    """Format an allele-frequency string for markdown display."""
    try:
        f = float(value)
        if f == 0:
            return "0"
        if f < 0.001:
            return f"{f:.2e}"
        return f"{f:.4f}"
    except (ValueError, TypeError):
        return value or "."


def get_impact(record: dict) -> str:
    """Return IMPACT field from a parsed VEP record, defaulting to MODIFIER."""
    return record.get("IMPACT", "MODIFIER").upper()


# -- report sections ----------------------------------------------------------

def build_summary(records: list[dict]) -> str:
    """Return the markdown summary section with impact and gene tables."""
    total = len(records)
    impact_counts: Counter = Counter(get_impact(r) for r in records)
    gene_counts: Counter = Counter(
        r.get("SYMBOL", r.get("Gene", ".")) for r in records
    )

    lines = [
        "## Сводка",
        "",
        f"- **Всего аннотированных вариантов:** {total}",
        "",
        "### По уровню воздействия (IMPACT)",
        "",
        "| IMPACT | Количество |",
        "|--------|-----------|",
    ]
    for impact in ("HIGH", "MODERATE", "LOW", "MODIFIER"):
        lines.append(f"| {impact} | {impact_counts.get(impact, 0)} |")

    lines += [
        "",
        "### По генам (топ-10)",
        "",
        "| Ген | Количество вариантов |",
        "|-----|---------------------|",
    ]
    for gene, cnt in gene_counts.most_common(10):
        lines.append(f"| {gene} | {cnt} |")

    return "\n".join(lines)


def build_top_variants(records: list[dict], vcf_info: dict[str, dict]) -> str:
    """Return top-5 HIGH-impact, rare variants as a markdown table."""
    high_rare = [
        r for r in records
        if get_impact(r) == "HIGH"
        and _to_float(r.get("gnomADe_AF", r.get("gnomAD_AF", "1"))) < 0.01
    ]
    high_rare.sort(
        key=lambda r: _to_float(r.get("gnomADe_AF", r.get("gnomAD_AF", "1")))
    )
    top5 = high_rare[:5]

    lines = [
        "## Топ-5 приоритетных вариантов",
        "",
        "| ВАРИАНТ | ГЕН | CONSEQUENCE | HGVSp | gnomAD_AF | ClinVar |",
        "|---------|-----|-------------|-------|-----------|---------|",
    ]
    for r in top5:
        variant_id = r.get("Uploaded_variation", ".")
        gene = r.get("SYMBOL", r.get("Gene", "."))
        consequence = r.get("Consequence", ".")
        hgvsp = r.get("HGVSp", ".")
        af_raw = r.get("gnomADe_AF", r.get("gnomAD_AF", "."))
        vcf_entry = vcf_info.get(variant_id, {})
        clinvar = vcf_entry.get("ClinVar", r.get("ClinVar_CLNSIG", "."))
        lines.append(
            f"| {variant_id} | {gene} | {consequence} "
            f"| {hgvsp} | {fmt_af(af_raw)} | {clinvar} |"
        )

    if not top5:
        lines.append(
            "| — | Нет HIGH-impact вариантов с gnomAD_AF < 0.01 | | | | |"
        )

    return "\n".join(lines)


def build_per_consequence(records: list[dict]) -> str:
    """Return one example variant per unique consequence category."""
    by_consequence: dict[str, dict] = {}
    for r in records:
        csq = r.get("Consequence", "unknown").split(",")[0]
        if csq not in by_consequence:
            by_consequence[csq] = r

    lines = [
        "## Пример варианта для каждой категории последствий",
        "",
        "| CONSEQUENCE | ВАРИАНТ | ГЕН | IMPACT | HGVSp / HGVSc |",
        "|-------------|---------|-----|--------|---------------|",
    ]
    for csq in sorted(by_consequence):
        r = by_consequence[csq]
        variant_id = r.get("Uploaded_variation", ".")
        gene = r.get("SYMBOL", r.get("Gene", "."))
        impact = get_impact(r)
        hgvs = r.get("HGVSp", r.get("HGVSc", "."))
        lines.append(f"| {csq} | {variant_id} | {gene} | {impact} | {hgvs} |")

    return "\n".join(lines)


def build_sv_section() -> str:
    """Return a placeholder section for SV/CNV results."""
    return (
        "## Структурные варианты (SV / CNV)\n"
        "\n"
        "Полный анализ структурных вариантов и CNV выполняется скриптом "
        "`scripts/30_sv_cnv_intro.sh`.\n"
        "Результаты доступны в директории `outputs/<run>/sv/`.\n"
        "\n"
        "Краткое описание:\n"
        "- Запуск: `bash scripts/30_sv_cnv_intro.sh`\n"
        "- Форматы вывода: VCF (структурные варианты), BED (CNV регионы)\n"
        "- Смотрите также: `data/input/demo_sv.vcf` для примеров SV записей\n"
    )


def build_report(records: list[dict], vcf_info: dict[str, dict]) -> str:
    """Assemble the full markdown report string."""
    today = date.today().isoformat()
    sections = [
        "# Отчёт: Интерпретация вариантов",
        "",
        f"_Дата генерации: {today}_",
        "",
        build_summary(records),
        "",
        build_top_variants(records, vcf_info),
        "",
        build_per_consequence(records),
        "",
        build_sv_section(),
        "",
        "---",
        f"_Отчёт создан автоматически скриптом `21_make_report.py` · {today}_",
    ]
    return "\n".join(sections) + "\n"


# -- entry point --------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--tsv",
        required=True,
        help="VEP output TSV file (e.g. outputs/<run>/vep/vep_output.tsv)",
    )
    ap.add_argument(
        "--vcf",
        required=True,
        help="Original VCF with QUAL/DP/AF/ClinVar INFO (e.g. data/input/demo_variants.vcf)",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output markdown report path (e.g. outputs/<run>/report.md)",
    )
    args = ap.parse_args()

    tsv_path = Path(args.tsv)
    vcf_path = Path(args.vcf)
    out_path = Path(args.out)

    if not tsv_path.exists():
        raise SystemExit(f"ERROR: VEP TSV not found: {tsv_path}")
    if not vcf_path.exists():
        raise SystemExit(f"ERROR: VCF not found: {vcf_path}")

    print(f"==> Parsing VEP TSV: {tsv_path}")
    records = parse_vep_tsv(tsv_path)
    print(f"    {len(records)} annotation record(s) loaded")

    print(f"==> Parsing VCF: {vcf_path}")
    vcf_info = parse_vcf(vcf_path)
    print(f"    {len(vcf_info)} variant(s) loaded")

    print("==> Building report ...")
    report_text = build_report(records, vcf_info)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report_text, encoding="utf-8")
    print(f"==> Report written: {out_path} ({out_path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
