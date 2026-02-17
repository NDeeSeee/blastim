#!/usr/bin/env python3
"""
_vep_rest_helper.py — Вспомогательный скрипт для VEP REST API.

Читает демо VCF, отправляет варианты в Ensembl VEP REST API пакетами,
записывает результат в TSV (тот же формат, что и VEP --tab).

Usage:
    python3 _vep_rest_helper.py --vcf data/input/demo_variants.vcf \
        --out outputs/run_*/vep/vep_output.tsv --assembly GRCh38
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.request

# VEP REST API endpoint
VEP_REST_URL = "https://rest.ensembl.org/vep/human/region"
BATCH_SIZE = 20  # max 200 per POST, keep small for reliability

HEADER_COLS = (
    "Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t"
    "Consequence\tcDNA_position\tCDS_position\tProtein_position\t"
    "Amino_acids\tCodons\tExisting_variation\tExtra"
)


def parse_vcf(vcf_path: str) -> list[dict]:
    """Parses a VCF and returns list of variant dicts."""
    variants = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, vid, ref, alt = parts[:5]
            variants.append({"chrom": chrom, "pos": int(pos), "id": vid,
                              "ref": ref, "alt": alt})
    return variants


def variant_to_rest_format(v: dict) -> str:
    """Converts variant dict to VEP REST region format: chr pos end allele +"""
    chrom = v["chrom"].lstrip("chr")
    start = v["pos"]
    end = start + len(v["ref"]) - 1
    return f"{chrom} {start} {end} {v['ref']}/{v['alt']} +"


def query_vep_rest(regions: list[str], assembly: str = "GRCh38",
                   retries: int = 3) -> list[dict]:
    """Sends a batch POST to VEP REST API and returns parsed JSON."""
    url = f"{VEP_REST_URL}?assembly={assembly}&canonical=1&hgvs=1&af_gnomade=1&symbol=1"
    payload = json.dumps({"regions": regions}).encode("utf-8")
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=payload, headers=headers)
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            print(f"  WARN: REST API ошибка (попытка {attempt}/{retries}): {e}",
                  file=sys.stderr)
            if attempt < retries:
                time.sleep(5 * attempt)
    raise RuntimeError("VEP REST API недоступен после нескольких попыток.")


def consequence_to_tsv(variant: dict, rest_response: list[dict]) -> list[str]:
    """Formats VEP REST response into TSV rows matching --tab output."""
    rows = []
    vid = variant["id"] if variant["id"] != "." else (
        f"{variant['chrom']}_{variant['pos']}_{variant['ref']}/{variant['alt']}")
    location = f"{variant['chrom']}:{variant['pos']}"

    for result in rest_response:
        for tc in result.get("transcript_consequences", [result]):
            # Extract fields
            gene = tc.get("gene_id", "-")
            feature = tc.get("transcript_id", "-")
            ftype = "Transcript" if "transcript_id" in tc else "-"
            csq = ",".join(tc.get("consequence_terms", ["-"]))
            cdna = tc.get("cdna_start", "-")
            cds = tc.get("cds_start", "-")
            prot = tc.get("protein_start", "-")
            aa = tc.get("amino_acids", "-")
            codons = tc.get("codons", "-")
            existing = tc.get("colocated_variants", [{}])
            existing_str = ",".join(
                v.get("id", "") for v in existing) if existing else "-"

            # Extra field
            impact = tc.get("impact", "-")
            symbol = tc.get("gene_symbol", "-")
            canonical = "YES" if tc.get("canonical") else "-"
            biotype = tc.get("biotype", "-")
            hgvsc = tc.get("hgvsc", "-")
            hgvsp = tc.get("hgvsp", "-")
            gnomad = tc.get("gnomade_af", "-")

            extra = (f"IMPACT={impact};STRAND=1;SYMBOL={symbol};"
                     f"CANONICAL={canonical};BIOTYPE={biotype};"
                     f"HGVSc={hgvsc};HGVSp={hgvsp};gnomADe_AF={gnomad}")

            rows.append(
                f"{vid}\t{location}\t{variant['alt']}\t{gene}\t{feature}\t"
                f"{ftype}\t{csq}\t{cdna}\t{cds}\t{prot}\t{aa}\t{codons}\t"
                f"{existing_str}\t{extra}"
            )
    return rows if rows else [
        f"{vid}\t{location}\t{variant['alt']}\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-"
    ]


def main():
    parser = argparse.ArgumentParser(description="VEP REST API helper")
    parser.add_argument("--vcf",      required=True, help="Input VCF file")
    parser.add_argument("--out",      required=True, help="Output TSV file")
    parser.add_argument("--assembly", default="GRCh38", help="Genome assembly")
    args = parser.parse_args()

    print(f"==> VEP REST: читаем {args.vcf} ...", file=sys.stderr)
    variants = parse_vcf(args.vcf)
    print(f"    Найдено вариантов: {len(variants)}", file=sys.stderr)

    all_rows: list[str] = []

    for i in range(0, len(variants), BATCH_SIZE):
        batch = variants[i:i + BATCH_SIZE]
        regions = [variant_to_rest_format(v) for v in batch]
        print(f"    Запрос {i + 1}-{i + len(batch)} / {len(variants)} ...",
              file=sys.stderr)
        try:
            results = query_vep_rest(regions, args.assembly)
            for variant, result in zip(batch, results
                                       if isinstance(results, list) else [results]):
                all_rows.extend(consequence_to_tsv(variant, [result]))
        except RuntimeError as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            sys.exit(1)
        # Соблюдаем rate limit (15 req/s для REST)
        time.sleep(0.1)

    with open(args.out, "w") as f:
        f.write(f"## VEP REST API output — assembly {args.assembly}\n")
        f.write(f"#Uploaded_variation\tLocation\tAllele\tGene\tFeature\t"
                f"Feature_type\tConsequence\tcDNA_position\tCDS_position\t"
                f"Protein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n")
        for row in all_rows:
            f.write(row + "\n")

    print(f"==> Записано {len(all_rows)} строк → {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
