#!/usr/bin/env python3
"""Assign a personalized set of variants for HW2.

Deterministically picks 5 variant positions based on the student's
GitHub username. Ensures a mix of categories:
  - 1-2 confident SNPs (high QUAL, high DP)
  - 1 indel
  - 1-2 borderline/low-quality variants (low QUAL or low DP)

Usage:
    export GITHUB_USER=<your_github>
    python scripts/90_homework_assign_variants.py \
        --vcf outputs/run_latest/variants/raw.vcf.gz \
        --out homework/submissions/$GITHUB_USER/variants_$GITHUB_USER.txt
"""

import argparse
import gzip
import hashlib
import os
import random
from pathlib import Path


def parse_vcf_records(vcf_path: Path) -> list[dict]:
    """Parse VCF records into a list of dicts with key fields."""
    records = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open

    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, _, ref, alt, qual, filt, info = parts[:8]

            try:
                qual_f = float(qual) if qual != "." else 0.0
            except ValueError:
                qual_f = 0.0

            # Parse DP from INFO
            dp = 0
            for field in info.split(";"):
                if field.startswith("DP="):
                    try:
                        dp = int(field[3:])
                    except ValueError:
                        pass
                    break

            # Handle multi-allelic ALT (e.g., "A,G") — use first allele
            first_alt = alt.split(",")[0]

            records.append(
                {
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "qual": qual_f,
                    "dp": dp,
                    "is_snp": len(ref) == 1 and len(first_alt) == 1,
                    "is_indel": len(ref) != len(first_alt),
                    "line": line.rstrip("\n"),
                }
            )
    return records


def categorize(records: list[dict]) -> dict[str, list[dict]]:
    """Split records into categories for balanced assignment."""
    cats: dict[str, list[dict]] = {
        "snp_high": [],
        "snp_low": [],
        "indel": [],
        "borderline": [],
    }

    for r in records:
        if r["is_snp"] and r["qual"] >= 100 and r["dp"] >= 20:
            cats["snp_high"].append(r)
        elif r["is_snp"] and (r["qual"] < 30 or r["dp"] < 10):
            cats["snp_low"].append(r)
        elif r["is_indel"] and r["qual"] >= 30:
            cats["indel"].append(r)
        elif r["qual"] < 30 or r["dp"] < 8:
            cats["borderline"].append(r)

    return cats


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Assign personalized variant positions for HW2."
    )
    ap.add_argument(
        "--vcf",
        required=True,
        help="Path to VCF file (raw.vcf.gz)",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output file (e.g. homework/submissions/<user>/variants_<user>.txt)",
    )
    ap.add_argument(
        "--n", type=int, default=5, help="Number of variants to assign (default: 5)"
    )
    args = ap.parse_args()

    user = os.environ.get("GITHUB_USER", "").strip()
    if not user:
        raise SystemExit(
            "ERROR: Set GITHUB_USER first:\n  export GITHUB_USER=<your_github>"
        )

    vcf_path = Path(args.vcf)
    if not vcf_path.exists():
        raise SystemExit(f"ERROR: VCF file not found: {vcf_path}")

    print(f"==> Parsing VCF: {vcf_path}")
    records = parse_vcf_records(vcf_path)
    if not records:
        raise SystemExit(f"ERROR: No variant records found in {vcf_path}")
    print(f"    Total records: {len(records)}")

    cats = categorize(records)
    print(
        f"    Categories: snp_high={len(cats['snp_high'])}, "
        f"indel={len(cats['indel'])}, "
        f"snp_low={len(cats['snp_low'])}, "
        f"borderline={len(cats['borderline'])}"
    )

    # Deterministic shuffle seeded by username
    seed = int(hashlib.sha256(user.encode()).hexdigest()[:16], 16)
    rng = random.Random(seed)

    for cat_list in cats.values():
        rng.shuffle(cat_list)

    # Pick balanced set: 2 snp_high, 1 indel, 1 snp_low/borderline, 1 borderline
    picked: list[dict] = []

    # 2 confident SNPs
    picked.extend(cats["snp_high"][: min(2, len(cats["snp_high"]))])
    # 1 indel
    if cats["indel"]:
        picked.append(cats["indel"][0])
    # 1 low-quality SNP
    if cats["snp_low"]:
        picked.append(cats["snp_low"][0])
    elif cats["borderline"]:
        picked.append(cats["borderline"][0])
    # 1 borderline
    remaining_border = [r for r in cats["borderline"] if r not in picked]
    if remaining_border:
        picked.append(remaining_border[0])

    # If we don't have enough, fill from any category
    if len(picked) < args.n:
        all_remaining = [r for r in records if r not in picked]
        rng.shuffle(all_remaining)
        picked.extend(all_remaining[: args.n - len(picked)])

    picked = picked[: args.n]

    # Sort by position for readability
    picked.sort(key=lambda r: r["pos"])

    # Write output
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    for r in picked:
        cat = "SNP_HIGH" if r["is_snp"] and r["qual"] >= 100 else \
              "SNP_LOW" if r["is_snp"] and r["qual"] < 30 else \
              "INDEL" if r["is_indel"] else \
              "BORDERLINE"
        lines.append(
            f"{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t"
            f"QUAL={r['qual']:.1f}\tDP={r['dp']}\t{cat}"
        )

    header = (
        f"# Variant assignment for {user} (HW2)\n"
        f"# Generated from: {vcf_path}\n"
        f"# CHROM\tPOS\tREF\tALT\tQUAL\tDP\tCATEGORY\n"
    )
    out_path.write_text(header + "\n".join(lines) + "\n")

    print(f"\n==> Assigned {len(picked)} variants for {user}")
    print(f"    Output: {out_path}")
    for i, r in enumerate(picked, 1):
        print(f"    {i}. {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} QUAL={r['qual']:.0f} DP={r['dp']}")


if __name__ == "__main__":
    main()
