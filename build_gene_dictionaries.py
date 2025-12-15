#!/usr/bin/env python3
"""
Build two gene dictionary JSON artifacts for DNABY:

1) gene_dictionary_pc_v1.json
   - Canonical coordinate system for sparse vectors (COO/CSR)
   - Source: GENCODE human GTF (protein-coding genes only)
   - Keys: Ensembl Gene IDs (ENSG...), version suffix stripped by default

2) gene_dictionary_clinvar_v1.json
   - Evidence overlay / clinical attention set (NOT a coordinate system)
   - Source: ClinVar tab_delimited (prefers gene_summary, falls back to variant_summary)
   - Keys: NCBI Gene IDs (Entrez GeneID), because ClinVar exports are most stable there

Usage:
  python3 build_gene_dictionaries.py
  python3 build_gene_dictionaries.py --out-dir ./artifacts

Notes:
- These two dictionaries intentionally use different ID systems because they come from
  different authoritative sources. You should map between them (Ensembl <-> Entrez)
  later in your DB layer (HGNC/Ensembl BioMart/HGNC mapping file).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Dict, Iterable, Optional, Tuple

import requests


# -----------------------------
# Defaults / Sources
# -----------------------------

DEFAULT_GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_49/gencode.v49.annotation.gtf.gz"
)

CLINVAR_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"
CLINVAR_CANDIDATES = [
    "gene_summary.txt.gz",      # preferred (smaller + direct)
    "variant_summary.txt.gz",   # fallback (large)
]


# -----------------------------
# Types / Helpers
# -----------------------------

@dataclass(frozen=True)
class SourceInfo:
    url: str
    fetched_at_utc: str
    http_etag: Optional[str]
    http_last_modified: Optional[str]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def http_head_exists(url: str, timeout: int = 30) -> bool:
    r = requests.head(url, timeout=timeout, allow_redirects=True)
    return r.status_code == 200


def http_download_bytes(url: str, timeout: int = 180) -> Tuple[bytes, SourceInfo]:
    r = requests.get(url, stream=True, timeout=timeout)
    r.raise_for_status()

    source = SourceInfo(
        url=url,
        fetched_at_utc=utc_now_iso(),
        http_etag=r.headers.get("ETag"),
        http_last_modified=r.headers.get("Last-Modified"),
    )

    buf = io.BytesIO()
    for chunk in r.iter_content(chunk_size=1024 * 1024):
        if chunk:
            buf.write(chunk)

    return buf.getvalue(), source


def write_json(path: str, payload: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    print(f"Wrote {path}")


# -----------------------------
# 1) GENCODE protein-coding dict
# -----------------------------

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_gtf_attributes(attr_field: str) -> Dict[str, str]:
    return {k: v for k, v in ATTR_RE.findall(attr_field)}


def strip_ensembl_version(ensembl_id: str) -> str:
    return ensembl_id.split(".", 1)[0]


def build_gencode_pc_dictionary(
    *,
    gtf_url: str,
    version: str,
    genome_build: str,
    keep_gene_version_suffix: bool,
) -> dict:
    # Stream GTF.GZ without saving to disk
    r = requests.get(gtf_url, stream=True, timeout=180)
    r.raise_for_status()

    source = SourceInfo(
        url=gtf_url,
        fetched_at_utc=utc_now_iso(),
        http_etag=r.headers.get("ETag"),
        http_last_modified=r.headers.get("Last-Modified"),
    )

    r.raw.decode_content = True
    gz = gzip.GzipFile(fileobj=r.raw, mode="rb")
    text = io.TextIOWrapper(gz, encoding="utf-8", newline="")

    gene_name_by_id: Dict[str, str] = {}

    for line in text:
        if not line or line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue

        feature = parts[2]
        if feature != "gene":
            continue

        attrs = parse_gtf_attributes(parts[8])

        gene_type = attrs.get("gene_type") or attrs.get("gene_biotype")
        if gene_type != "protein_coding":
            continue

        gene_id = attrs.get("gene_id")
        if not gene_id:
            continue

        if not keep_gene_version_suffix:
            gene_id = strip_ensembl_version(gene_id)

        gene_name = attrs.get("gene_name", "")

        # first win
        gene_name_by_id.setdefault(gene_id, gene_name)

    gene_ids_sorted = sorted(gene_name_by_id.keys())
    genes_index = {gid: i for i, gid in enumerate(gene_ids_sorted)}

    return {
        "version": version,
        "source": {
            "url": source.url,
            "fetched_at_utc": source.fetched_at_utc,
            "http_etag": source.http_etag,
            "http_last_modified": source.http_last_modified,
        },
        "genome_build": genome_build,
        "gene_id_system": "Ensembl",
        "gene_class": "protein_coding",
        "keep_gene_version_suffix": keep_gene_version_suffix,
        "gene_count": len(genes_index),
        "genes": genes_index,
        "gene_name_by_id": {gid: gene_name_by_id.get(gid, "") for gid in gene_ids_sorted},
    }


# -----------------------------
# 2) ClinVar-derived dict
# -----------------------------

def open_gz_tsv_bytes(gz_bytes: bytes) -> Iterable[Dict[str, str]]:
    with gzip.GzipFile(fileobj=io.BytesIO(gz_bytes), mode="rb") as gz:
        text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
        reader = csv.DictReader(text, delimiter="\t")
        for row in reader:
            yield { (k or "").strip(): (v or "").strip() for k, v in row.items() }


def pick_column(headers: Iterable[str], wanted: Iterable[str]) -> Optional[str]:
    header_set = {h.strip(): h for h in headers if h}
    for w in wanted:
        if w in header_set:
            return header_set[w]
    lower_map = {h.lower(): h for h in header_set.keys()}
    for w in wanted:
        if w.lower() in lower_map:
            return lower_map[w.lower()]
    return None


def extract_geneid_symbol_from_variant_row(row: Dict[str, str]) -> Iterable[Tuple[str, Optional[str]]]:
    keys = list(row.keys())
    gene_id_col = pick_column(keys, ["GeneID", "Gene ID", "EntrezGeneID", "Entrez Gene ID"])
    gene_symbol_col = pick_column(keys, ["GeneSymbol", "Gene Symbol", "Symbol", "Gene(s)"])

    gene_id_val = row.get(gene_id_col, "") if gene_id_col else ""
    gene_symbol_val = row.get(gene_symbol_col, "") if gene_symbol_col else ""

    # Simple case: one integer id
    if gene_id_val and re.fullmatch(r"[0-9]+", gene_id_val):
        yield (gene_id_val, gene_symbol_val or None)
        return

    # Split lists
    id_tokens = [t for t in re.split(r"[;,| ]+", gene_id_val) if t]
    sym_tokens = [t for t in re.split(r"[;,| ]+", gene_symbol_val) if t]

    id_ints = [t for t in id_tokens if re.fullmatch(r"[0-9]+", t)]
    if id_ints:
        for i, gid in enumerate(id_ints):
            sym = sym_tokens[i] if i < len(sym_tokens) else None
            yield (gid, sym)
        return

    # Combined tokens like "7157:TP53"
    combined = []
    combined.extend(re.split(r"[;|]+", gene_id_val))
    combined.extend(re.split(r"[;|]+", gene_symbol_val))
    for token in combined:
        token = token.strip()
        if not token:
            continue
        m = re.match(r"(?P<gid>[0-9]+)\s*[:/]\s*(?P<sym>[A-Za-z0-9_.-]+)", token)
        if m:
            yield (m.group("gid"), m.group("sym"))


def build_clinvar_dictionary(*, version: str) -> dict:
    chosen_url = None
    response = None
    last_error: Optional[Exception] = None

    # Prefer the smaller gene_summary, but fall back to variant_summary if needed.
    for name in CLINVAR_CANDIDATES:
        url = CLINVAR_BASE + name
        try:
            resp = requests.get(url, stream=True, timeout=180)
            resp.raise_for_status()
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            continue
        chosen_url = url
        response = resp
        break

    if not chosen_url or response is None:
        raise RuntimeError(f"Could not fetch ClinVar exports at {CLINVAR_BASE} ({CLINVAR_CANDIDATES})") from last_error

    with response as r:
        source = SourceInfo(
            url=chosen_url,
            fetched_at_utc=utc_now_iso(),
            http_etag=r.headers.get("ETag"),
            http_last_modified=r.headers.get("Last-Modified"),
        )

        # Stream + decompress directly to avoid buffering the multi-GB variant_summary.gz in memory.
        r.raw.decode_content = True
        gz = gzip.GzipFile(fileobj=r.raw, mode="rb")
        text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
        reader = csv.DictReader(text, delimiter="\t")

        try:
            first_row = next(reader)
        except StopIteration:
            raise RuntimeError(f"{chosen_url} appears to be empty")

        headers = reader.fieldnames or list(first_row.keys())
        is_gene_summary = chosen_url.endswith("gene_summary.txt.gz")

        gene_symbol_by_id: Dict[str, str] = {}

        def record_symbol(gid: str, sym: Optional[str]) -> None:
            sym_val = sym or ""
            existing = gene_symbol_by_id.get(gid)
            if existing is None or (not existing and sym_val):
                gene_symbol_by_id[gid] = sym_val

        if is_gene_summary:
            gene_id_col = pick_column(headers, ["GeneID", "Gene ID", "EntrezGeneID", "Entrez Gene ID"])
            gene_symbol_col = pick_column(headers, ["GeneSymbol", "Gene Symbol", "Symbol"])
            if not gene_id_col:
                raise RuntimeError(f"ClinVar gene_summary missing GeneID column. Headers: {headers}")

            def process(row: Dict[str, str]) -> None:
                gid = (row.get(gene_id_col, "") or "").strip()
                if not gid or not re.fullmatch(r"[0-9]+", gid):
                    return
                sym = (row.get(gene_symbol_col, "") or "").strip() if gene_symbol_col else ""
                record_symbol(gid, sym)

            process(first_row)
            for row in reader:
                process(row)

        else:
            # variant_summary
            for gid, sym in extract_geneid_symbol_from_variant_row(first_row):
                if gid:
                    record_symbol(gid, sym)

            for row in reader:
                for gid, sym in extract_geneid_symbol_from_variant_row(row):
                    if gid:
                        record_symbol(gid, sym)

    gene_ids_sorted = sorted(gene_symbol_by_id.keys(), key=lambda x: int(x))
    genes_index = {gid: i for i, gid in enumerate(gene_ids_sorted)}

    return {
        "version": version,
        "source": {
            "url": source.url,
            "fetched_at_utc": source.fetched_at_utc,
            "http_etag": source.http_etag,
            "http_last_modified": source.http_last_modified,
        },
        "genome_build": None,  # ClinVar gene list is not build-specific; variants are.
        "gene_id_system": "NCBI_GeneID",
        "gene_class": "clinvar_exposed_loci",
        "gene_count": len(genes_index),
        "genes": genes_index,
        "gene_symbol_by_id": {gid: gene_symbol_by_id.get(gid, "") for gid in gene_ids_sorted},
    }


# -----------------------------
# Main / CLI
# -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=".", help="Output directory for JSON files")
    ap.add_argument("--gencode-gtf-url", default=DEFAULT_GENCODE_GTF_URL, help="GENCODE GTF.GZ URL")
    ap.add_argument("--genome-build", default="GRCh38", help="Genome build label for GENCODE dictionary metadata")

    ap.add_argument("--pc-version", default="pc_v1", help="Version tag for protein-coding dictionary")
    ap.add_argument("--clinvar-version", default="clinvar_v1", help="Version tag for ClinVar dictionary")

    ap.add_argument("--keep-gene-version-suffix", action="store_true",
                    help="Keep Ensembl gene_id version suffix (ENSG... .N). Default strips it.")

    args = ap.parse_args()

    out_pc = os.path.join(args.out_dir, "gene_dictionary_pc_v1.json")
    out_clinvar = os.path.join(args.out_dir, "gene_dictionary_clinvar_v1.json")

    if os.path.exists(out_pc):
        print(f"Skipping existing {out_pc} (GENCODE dictionary is treated as frozen).")
    else: 
        pc = build_gencode_pc_dictionary(
            gtf_url=args.gencode_gtf_url,
            version=args.pc_version,
            genome_build=args.genome_build,
            keep_gene_version_suffix=args.keep_gene_version_suffix,
        )
        write_json(out_pc, pc)

    clinvar = build_clinvar_dictionary(version=args.clinvar_version)
    write_json(out_clinvar, clinvar)

    print("\nDone.")
    print(f"- {out_pc}")
    print(f"- {out_clinvar}")


if __name__ == "__main__":
    main()
