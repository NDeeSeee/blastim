#!/usr/bin/env python3
"""Assign a personalized set of 8 variants for HW3.

Deterministically selects variants based on the student's GitHub username.
Tries to include a mix:
  - 2 HIGH impact variants (pathogenic/likely_pathogenic)
  - 2 MODERATE impact (missense VUS)
  - 1 LOW impact (synonymous or benign)
  - 1 intronic/UTR
  - 1 indel/frameshift
  - 1 SV record

Usage:
    export GITHUB_USER=<your_github>
    python3 scripts/90_assign_variants.py \
        --vep_tsv data/input/demo_variants_vep.tsv \
        --sv_vcf data/input/demo_sv.vcf \
        --out homework/submissions/$GITHUB_USER/variants_$GITHUB_USER.txt
"""

import argparse
import hashlib
import os
import random
from pathlib import Path


# -- parsers ------------------------------------------------------------------

def parse_vep_tsv(tsv_path: Path) -> list[dict]:
    """Parse VEP TSV output into a list of dicts.

    Skips '##' header comment lines.  The first '#'-prefixed line is the
    column header.  The Extra field is expanded into individual key=value
    pairs so IMPACT, SYMBOL, ClinVar_CLNSIG etc. are top-level keys.
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
            extra = row.get("Extra", "")
            for item in extra.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    row[k] = v
            records.append(row)
    return records


def parse_sv_vcf(sv_vcf_path: Path) -> list[str]:
    """Parse SV VCF; return list of formatted SV strings.

    Each entry: 'CHROM:POS\tSVTYPE\tSVLEN\tGENE'
    """
    svs = []
    if not sv_vcf_path.exists():
        return svs
    with sv_vcf_path.open(encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos = parts[0], parts[1]
            info_str = parts[7]
            info: dict[str, str] = {}
            for item in info_str.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
            svtype = info.get("SVTYPE", "SV")
            svlen = info.get("SVLEN", ".")
            gene = info.get("GENE", ".")
            svs.append(f"{chrom}:{pos}\t{svtype}\t{svlen}\t{gene}")
    return svs


# -- categorisation -----------------------------------------------------------

def categorise_records(records: list[dict]) -> dict[str, list[dict]]:
    """Sort VEP records into variant-class buckets.

    Buckets:
        high_pathogenic   – HIGH impact + Pathogenic or Likely_pathogenic
        moderate_vus      – MODERATE impact + Uncertain_significance
        low_benign        – LOW impact
        modifier_intronic – MODIFIER impact (intronic, UTR, upstream …)
        indel_frameshift  – contains frameshift or inframe in consequence
    """
    cats: dict[str, list[dict]] = {
        "high_pathogenic": [],
        "moderate_vus": [],
        "low_benign": [],
        "modifier_intronic": [],
        "indel_frameshift": [],
    }
    for r in records:
        impact = r.get("IMPACT", "").upper()
        clinvar = r.get("ClinVar_CLNSIG", "")
        consequence = r.get("Consequence", "").lower()

        if "frameshift" in consequence or "inframe" in consequence:
            cats["indel_frameshift"].append(r)
        elif impact == "HIGH" and any(
            s in clinvar for s in ("Pathogenic", "Likely_pathogenic")
        ):
            cats["high_pathogenic"].append(r)
        elif impact == "MODERATE" and "Uncertain" in clinvar:
            cats["moderate_vus"].append(r)
        elif impact == "LOW":
            cats["low_benign"].append(r)
        elif impact == "MODIFIER":
            cats["modifier_intronic"].append(r)

    return cats


def record_to_line(r: dict) -> str:
    """Format a VEP record as a tab-separated output line."""
    variant_id = r.get("Uploaded_variation", ".")
    gene = r.get("SYMBOL", r.get("Gene", "."))
    consequence = r.get("Consequence", ".")
    impact = r.get("IMPACT", ".")
    return f"{variant_id}\t{gene}\t{consequence}\t{impact}"


# -- deterministic selection --------------------------------------------------

def seeded_pick(items: list, rng: random.Random, n: int) -> list:
    """Return up to *n* items from *items* using a deterministic shuffle."""
    shuffled = list(items)
    rng.shuffle(shuffled)
    return shuffled[:n]


def build_variant_list(
    cats: dict[str, list[dict]], rng: random.Random, n: int
) -> list[str]:
    """Build the mixed candidate list and pick *n* SNV/indel lines."""
    # Target distribution (adjusted if not enough candidates in a bucket)
    targets = [
        ("high_pathogenic", 2),
        ("moderate_vus", 2),
        ("low_benign", 1),
        ("modifier_intronic", 1),
        ("indel_frameshift", 1),
    ]

    picked_lines: list[str] = []
    seen_ids: set[str] = set()

    for cat_name, want in targets:
        candidates = cats.get(cat_name, [])
        chosen = seeded_pick(candidates, rng, want)
        for r in chosen:
            vid = r.get("Uploaded_variation", ".")
            if vid not in seen_ids:
                seen_ids.add(vid)
                picked_lines.append(record_to_line(r))

    # If we are short, fill from remaining records in any category
    if len(picked_lines) < n - 1:  # -1 reserved for SV
        all_remaining = [
            r for bucket in cats.values()
            for r in bucket
            if r.get("Uploaded_variation", ".") not in seen_ids
        ]
        extras = seeded_pick(all_remaining, rng, (n - 1) - len(picked_lines))
        for r in extras:
            vid = r.get("Uploaded_variation", ".")
            if vid not in seen_ids:
                seen_ids.add(vid)
                picked_lines.append(record_to_line(r))

    return picked_lines[: n - 1]  # keep slot for SV


# -- main ---------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--vep_tsv",
        required=True,
        help="VEP-annotated TSV (e.g. data/input/demo_variants_vep.tsv)",
    )
    ap.add_argument(
        "--sv_vcf",
        required=True,
        help="SV VCF file (e.g. data/input/demo_sv.vcf)",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output variant list (e.g. homework/submissions/<user>/variants_<user>.txt)",
    )
    ap.add_argument(
        "--n",
        type=int,
        default=8,
        help="Total number of variants to assign, including 1 SV (default: 8)",
    )
    args = ap.parse_args()

    user = os.environ.get("GITHUB_USER", "").strip()
    if not user:
        raise SystemExit(
            "ERROR: Set GITHUB_USER first:\n  export GITHUB_USER=<your_github>"
        )

    vep_tsv = Path(args.vep_tsv)
    sv_vcf = Path(args.sv_vcf)
    out_path = Path(args.out)

    if not vep_tsv.exists():
        raise SystemExit(f"ERROR: VEP TSV not found: {vep_tsv}")

    # Deterministic RNG seeded by username
    seed = int(hashlib.sha256(user.encode()).hexdigest()[:16], 16)
    rng = random.Random(seed)

    print(f"==> Parsing VEP TSV: {vep_tsv}")
    records = parse_vep_tsv(vep_tsv)
    print(f"    {len(records)} record(s) loaded")

    cats = categorise_records(records)
    for cat_name, bucket in cats.items():
        print(f"    {cat_name}: {len(bucket)} candidate(s)")

    snv_lines = build_variant_list(cats, rng, args.n)

    print(f"==> Parsing SV VCF: {sv_vcf}")
    svs = parse_sv_vcf(sv_vcf)
    if svs:
        chosen_sv = seeded_pick(svs, rng, 1)[0]
        sv_line = f"{chosen_sv}\tSV"
    else:
        sv_line = "NO_SV_FOUND\t.\t.\t.\tSV"

    all_lines = snv_lines + [sv_line]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(all_lines) + "\n", encoding="utf-8")

    print(f"\n==> Assigned {len(all_lines)} variants for {user}")
    print(f"    Output: {out_path}")
    for i, line in enumerate(all_lines, 1):
        print(f"    {i}. {line}")


if __name__ == "__main__":
    main()
