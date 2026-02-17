#!/usr/bin/env python3
"""Quick self-check for HW3 submissions.

Verifies that the required files exist and contain expected keywords.

Usage:
    python3 scripts/91_homework_check.py --user <GITHUB_USER>
"""

import argparse
import sys
from pathlib import Path

REQUIRED_KEYWORDS = ["QUAL", "DP", "VEP", "IMPACT", "gnomAD", "SVTYPE"]


def check_file(path: Path, label: str) -> bool:
    """Return True if *path* exists and is non-empty; print OK/FAIL."""
    if not path.exists():
        print(f"[FAIL] {label}: file not found -> {path}")
        return False
    if path.stat().st_size == 0:
        print(f"[FAIL] {label}: file is empty -> {path}")
        return False
    print(f"[ OK ] {label}: exists ({path.stat().st_size} bytes)")
    return True


def check_keywords(path: Path, label: str) -> None:
    """Warn about any expected keywords missing from *path*."""
    if not path.exists():
        return
    text = path.read_text(errors="ignore")
    for kw in REQUIRED_KEYWORDS:
        if kw not in text:
            print(f"[WARN] {label}: missing keyword '{kw}'")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--user", required=True, help="GitHub username")
    ap.add_argument(
        "--root",
        default="homework/submissions",
        help="Root directory for submissions (default: homework/submissions)",
    )
    args = ap.parse_args()

    base = Path(args.root) / args.user
    variants_file = base / f"variants_{args.user}.txt"
    report_file = base / f"hw3_{args.user}.md"

    print(f"Checking HW3 for: {args.user}")
    print(f"Directory: {base}")
    print()

    ok = True
    ok &= check_file(variants_file, "variants list")
    ok &= check_file(report_file, "HW3 report")

    print()
    if ok:
        check_keywords(report_file, "HW3 report")
        print()
        print("[ OK ] All files present. Warnings above (if any) are hints, not errors.")
    else:
        print("[FAIL] Some required files are missing.")
        sys.exit(1)


if __name__ == "__main__":
    main()
