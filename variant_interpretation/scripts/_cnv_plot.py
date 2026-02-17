#!/usr/bin/env python3
"""
_cnv_plot.py — ASCII и matplotlib визуализация SV/CNV для демо.

Используется из 30_sv_cnv_intro.sh (необязательно).
Генерирует простую ASCII-карту SV на хромосоме.

Usage:
    python3 _cnv_plot.py --sv_vcf data/input/demo_sv.vcf \
        --out outputs/run_*/sv/sv_map.txt
"""

import argparse
import sys
from pathlib import Path


def parse_sv_vcf(vcf_path: Path) -> list[dict]:
    svs = []
    with vcf_path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            info: dict[str, str] = {}
            for item in parts[7].split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
            svlen = abs(int(info.get("SVLEN", "0") or "0"))
            svs.append({
                "id":     parts[2],
                "chrom":  parts[0],
                "pos":    int(parts[1]),
                "svtype": info.get("SVTYPE", "?"),
                "svlen":  svlen,
                "end":    int(info.get("END", str(int(parts[1]) + svlen))),
                "gene":   info.get("GENE", "."),
                "filter": parts[6],
            })
    return svs


def ascii_map(svs: list[dict], width: int = 60) -> str:
    """Draw a simple ASCII bar chart of SV positions per chromosome."""
    if not svs:
        return "Нет SV данных.\n"

    # Group by chrom
    by_chrom: dict[str, list[dict]] = {}
    for sv in svs:
        by_chrom.setdefault(sv["chrom"], []).append(sv)

    SVTYPE_CHAR = {"DEL": "▼", "DUP": "▲", "INV": "↔", "BND": "⊗", "INS": "▶"}

    lines = []
    lines.append(f"{'─' * (width + 20)}")
    lines.append(f"  Карта структурных вариантов (ASCII)")
    lines.append(f"  DEL=▼  DUP=▲  INV=↔  BND=⊗  INS=▶")
    lines.append(f"{'─' * (width + 20)}")

    for chrom, chrom_svs in sorted(by_chrom.items()):
        max_end = max(sv["end"] for sv in chrom_svs)
        chrom_len = max_end * 1.1  # 10% padding

        lines.append(f"\n  {chrom}  (0 — {max_end:,} bp)")
        lines.append(f"  {'─' * width}")

        for sv in sorted(chrom_svs, key=lambda x: x["pos"]):
            char = SVTYPE_CHAR.get(sv["svtype"], "?")
            # Scale position to bar width
            start_x = int(sv["pos"] / chrom_len * width)
            end_x = int(sv["end"] / chrom_len * width)
            bar_len = max(1, end_x - start_x)

            bar = " " * start_x + char * bar_len
            flag = " [LowQual]" if sv["filter"] not in ("PASS", ".") else ""
            lines.append(
                f"  {bar:<{width}}  {sv['id']} {sv['svtype']} "
                f"{sv['svlen']:,}bp {sv['gene']}{flag}"
            )

        lines.append(f"  {'─' * width}")

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description="ASCII SV map generator")
    parser.add_argument("--sv_vcf", required=True, help="SV VCF file")
    parser.add_argument("--out", default="-", help="Output file (- for stdout)")
    args = parser.parse_args()

    svs = parse_sv_vcf(Path(args.sv_vcf))
    chart = ascii_map(svs)

    if args.out == "-":
        print(chart)
    else:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(chart)
        print(f"==> SV map saved: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
