#!/usr/bin/env python3
"""Quick self-check for HW2 submissions.

Verifies that the required files exist and contain expected keywords.

Usage:
    python scripts/91_homework_check.py --user <GITHUB_USER>
"""

import argparse
import sys
from pathlib import Path

REQUIRED_KEYWORDS_DEFENSE = ["QUAL", "DP", "FILTER"]


def check_file(path: Path, label: str) -> bool:
    if not path.exists():
        print(f"[FAIL] {label}: file not found -> {path}")
        return False
    if path.stat().st_size == 0:
        print(f"[FAIL] {label}: file is empty -> {path}")
        return False
    print(f"[ OK ] {label}: exists ({path.stat().st_size} bytes)")
    return True


def check_keywords(path: Path, label: str, keywords: list[str]) -> None:
    if not path.exists():
        return
    text = path.read_text(errors="ignore")
    for kw in keywords:
        if kw.lower() not in text.lower():
            print(f"[WARN] {label}: missing keyword '{kw}'")


def check_variant_count(path: Path, label: str, expected: int = 5) -> None:
    if not path.exists():
        return
    lines = [
        line
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if len(lines) < expected:
        print(
            f"[WARN] {label}: found {len(lines)} variant lines, expected {expected}"
        )
    else:
        print(f"[ OK ] {label}: {len(lines)} variant lines")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--user", required=True, help="GitHub username")
    ap.add_argument(
        "--root",
        default="homework/submissions",
        help="Root directory for submissions",
    )
    args = ap.parse_args()

    base = Path(args.root) / args.user
    variants = base / f"variants_{args.user}.txt"
    defense = base / f"defense_{args.user}.md"

    print(f"Checking HW2 for: {args.user}")
    print(f"Directory: {base}")
    print()

    ok = True
    ok &= check_file(variants, "variants list")
    ok &= check_file(defense, "variant defense")

    print()
    if ok:
        check_variant_count(variants, "variants list")
        check_keywords(defense, "variant defense", REQUIRED_KEYWORDS_DEFENSE)
        print()
        print(
            "[OK] All files present. Warnings above (if any) are hints, not errors."
        )
    else:
        print("[FAIL] Some required files are missing.")
        sys.exit(1)


if __name__ == "__main__":
    main()
