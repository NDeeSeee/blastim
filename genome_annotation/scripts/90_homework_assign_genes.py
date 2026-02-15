#!/usr/bin/env python3
"""Assign a personalized set of genes for HW1.

Deterministically picks N gene IDs based on the student's GitHub username.
Prefers genes that have both DIAMOND and Pfam hits (so there is material
to analyse).

Usage:
    export GITHUB_USER=<your_github>
    python scripts/90_homework_assign_genes.py \
        --run_dir outputs/run_latest \
        --out homework/submissions/$GITHUB_USER/genes_$GITHUB_USER.txt
"""

import argparse
import hashlib
import os
import random
import re
from pathlib import Path


def parse_gff_cds_ids(gff_path: Path) -> list[str]:
    """Extract protein IDs from CDS features in a GFF file."""
    ids = []
    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            match = re.search(r"\bID=([^;]+)", parts[8])
            if match:
                ids.append(match.group(1))
    return ids


def ids_present_in_file(ids: list[str], path: Path) -> set[str]:
    """Return the subset of *ids* that appear anywhere in *path*."""
    if not path.exists():
        return set()
    text = path.read_text(errors="ignore")
    return {i for i in ids if i in text}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--run_dir",
        required=True,
        help="Path to the pipeline run directory (e.g. outputs/run_latest)",
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output file for the gene list (e.g. homework/submissions/<user>/genes_<user>.txt)",
    )
    ap.add_argument("--n", type=int, default=12, help="Number of genes to assign (default: 12)")
    args = ap.parse_args()

    user = os.environ.get("GITHUB_USER", "").strip()
    if not user:
        raise SystemExit(
            "ERROR: Set GITHUB_USER first:\n  export GITHUB_USER=<your_github>"
        )

    run_dir = Path(args.run_dir)
    gff = run_dir / "prokka" / "ANNOT.gff"
    diamond = run_dir / "diamond" / "blastp_results.tsv"
    pfam = run_dir / "hmmer" / "hmmscan_domtblout.txt"

    if not gff.exists():
        raise SystemExit(f"ERROR: GFF not found: {gff}")

    all_ids = parse_gff_cds_ids(gff)
    if not all_ids:
        raise SystemExit(f"ERROR: No CDS IDs found in {gff}")

    # Prefer genes present in both DIAMOND and Pfam results
    in_diamond = ids_present_in_file(all_ids, diamond)
    in_pfam = ids_present_in_file(all_ids, pfam)
    preferred = sorted(in_diamond & in_pfam)

    if len(preferred) < args.n:
        # Fall back to any gene with at least one type of evidence
        preferred = sorted(in_diamond | in_pfam) or all_ids

    # Deterministic shuffle seeded by username
    seed = int(hashlib.sha256(user.encode()).hexdigest()[:16], 16)
    rng = random.Random(seed)
    rng.shuffle(preferred)

    picked = preferred[: args.n]

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(picked) + "\n")

    print(f"==> Assigned {len(picked)} genes for {user}")
    print(f"    Output: {out_path}")
    for i, gid in enumerate(picked, 1):
        print(f"    {i}. {gid}")


if __name__ == "__main__":
    main()
