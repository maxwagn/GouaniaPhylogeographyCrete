#!/usr/bin/env python3
"""
make_gene_name_table_zfish_fallback.py  (Python 3.7+)


EXAMPLE:
python make_gene_name_table_zfish_fallback.py \
        background_genes.txt \
        --out background_genes_UPDATE.tsv \
        --ncbi-chunk 150 \
        --ncbi-sleep 0.2 \
        --homology-sleep 0.2


Fixes:
- Chunks NCBI ESummary requests to avoid HTTP 414 (Request-URI Too Long).

What it does
------------
1) Reads a gene list file (one entry per line) containing mixed symbols + ENSGWIG IDs.
2) Extracts ENSGWIG IDs (if present anywhere in a line).
3) Batch-queries Ensembl REST /lookup/id for Gouania ENSGWIG IDs to get:
   - display_name, description, biotype
4) If Ensembl display_name is missing or "novel gene", parses external IDs from Ensembl description:
   - [Source:NCBI gene;Acc:<GeneID>]
   and queries NCBI to get NCBI symbol (often LOC* for non-model species).
5) If NCBI symbol is LOC* (or missing) for a gene that also lacks a good Ensembl symbol,
   tries to get a zebrafish ortholog symbol via Ensembl homology (best single ortholog) and zebrafish lookup.
6) Writes a TSV including ALL original input entries (including non-ENS symbols).

Short_id priority
-----------------
1) Ensembl display_name if informative (not missing/"novel gene")
2) NCBI symbol if present and NOT LOC*
3) Zebrafish ortholog symbol if available (used when NCBI is LOC* or missing and Ensembl name is missing/novel)
4) ENSGWIG ID

Output columns
--------------
input_gene, ens_id, short_id, short_id_source,
ensembl_display_name, ensembl_description, ensembl_biotype,
external_source, external_acc, ncbi_symbol,
zfish_gene_id, zfish_symbol, zfish_description,
orthology_type, orthology_perc_id
"""

import re
import json
import time
import argparse
import urllib.request
import urllib.parse
from typing import Dict, List, Any, Optional, Tuple

import pandas as pd

# Ensembl REST endpoints
ENSEMBL_LOOKUP_POST = "https://rest.ensembl.org/lookup/id"
HOMOLOGY_URL_TMPL = (
    "https://rest.ensembl.org/homology/id/{species}/{gene_id}"
    "?target_species=danio_rerio;type=orthologues;sequence=none"
)

# Regexes
ENS_GOU_RE = re.compile(r"ENSGWIG\d+", re.IGNORECASE)
ENS_ZF_RE = re.compile(r"ENSDARG\d+", re.IGNORECASE)
ACC_RE = re.compile(r"Source:([^;\]]+);Acc:([A-Za-z0-9\-]+)")

ORTHO_PRIORITY = {
    "ortholog_one2one": 0,
    "ortholog_one2many": 1,
    "ortholog_many2many": 2
}

def clean_field(x: str) -> str:
    return (x or "").replace("\t", " ").replace("\n", " ").strip()

def chunked(lst: List[str], n: int) -> List[List[str]]:
    return [lst[i:i+n] for i in range(0, len(lst), n)]

def http_post_json(url: str, payload: Dict[str, Any], timeout: int = 90) -> Any:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=data,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "make_gene_name_table_zfish_fallback/1.1",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode("utf-8")
    return json.loads(body)

def http_get_json(url: str, timeout: int = 90) -> Any:
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "make_gene_name_table_zfish_fallback/1.1",
        },
        method="GET",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode("utf-8")
    return json.loads(body)

def parse_external_acc(desc: str) -> Tuple[str, str]:
    if not isinstance(desc, str):
        return ("", "")
    m = ACC_RE.search(desc)
    if not m:
        return ("", "")
    return (m.group(1).strip(), m.group(2).strip())

def ensembl_good_symbol(display_name: str) -> bool:
    if not display_name:
        return False
    dn = display_name.strip()
    if not dn:
        return False
    if dn.lower() == "novel gene":
        return False
    return True

def ncbi_geneid_to_symbol(
    gene_ids: List[str],
    chunk_size: int = 200,
    sleep: float = 0.15,
    timeout: int = 90
) -> Dict[str, str]:
    """
    Map numeric NCBI GeneIDs to gene symbols using NCBI ESummary JSON.
    Uses chunking to avoid HTTP 414 (Request-URI Too Long).
    Returns dict {gene_id: symbol}.
    """
    gene_ids = [g for g in gene_ids if isinstance(g, str) and g.isdigit()]
    gene_ids = sorted(set(gene_ids))
    if not gene_ids:
        return {}

    out: Dict[str, str] = {}
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    for i in range(0, len(gene_ids), chunk_size):
        chunk = gene_ids[i:i + chunk_size]
        url = base + "?db=gene&id=" + ",".join(chunk) + "&retmode=json"
        with urllib.request.urlopen(url, timeout=timeout) as r:
            data = json.loads(r.read().decode("utf-8"))

        uids = data.get("result", {}).get("uids", [])
        for gid in uids:
            rec = data["result"].get(gid, {})
            out[gid] = (rec.get("name", "") or "").strip()

        time.sleep(sleep)

    return out

def pick_best_ortholog(homologies: List[Dict[str, Any]]) -> Tuple[str, str, float]:
    """
    Returns (zfish_gene_id, orthology_type, perc_id)
    """
    candidates = []
    for h in homologies:
        otype = h.get("type", "") or ""
        target = h.get("target", {}) or {}
        z_id = target.get("id", "") or ""
        if not z_id or not ENS_ZF_RE.match(z_id):
            continue

        perc = h.get("perc_id", None)
        if perc is None:
            perc = target.get("perc_id", None)
        try:
            perc_f = float(perc) if perc is not None else -1.0
        except Exception:
            perc_f = -1.0

        prio = ORTHO_PRIORITY.get(otype, 99)
        candidates.append((prio, -perc_f, z_id, otype, perc_f))

    if not candidates:
        return ("", "", -1.0)

    candidates.sort()
    _, _, z_id, otype, perc_f = candidates[0]
    return (z_id, otype, perc_f)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="Gene list file (one entry per line; mixed symbols + ENSGWIG IDs)")
    ap.add_argument("--out", default="gene_name_table.tsv", help="Output TSV filename")
    ap.add_argument("--species", default="gouania_willdenowi", help="Ensembl species name (default: gouania_willdenowi)")

    ap.add_argument("--lookup-chunk", type=int, default=800, help="IDs per Ensembl /lookup/id batch (default: 800)")
    ap.add_argument("--sleep", type=float, default=0.25, help="Sleep between Ensembl batch calls (default: 0.25s)")
    ap.add_argument("--homology-sleep", type=float, default=0.08, help="Sleep between homology calls (default: 0.08s)")

    ap.add_argument("--ncbi-chunk", type=int, default=200, help="NCBI GeneIDs per request (default: 200)")
    ap.add_argument("--ncbi-sleep", type=float, default=0.15, help="Sleep between NCBI requests (default: 0.15s)")
    args = ap.parse_args()

    # Read input
    with open(args.input, "r", encoding="utf-8") as f:
        raw_lines = [ln.rstrip("\n") for ln in f]
    input_genes = [ln.strip() for ln in raw_lines if ln.strip()]

    # Base rows and ENSGWIG IDs
    base_rows = []
    ens_ids = []
    for g in input_genes:
        m = ENS_GOU_RE.search(g)
        ens_id = m.group(0).upper() if m else ""
        if ens_id:
            ens_ids.append(ens_id)
        base_rows.append({"input_gene": g, "ens_id": ens_id})

    ens_ids = sorted(set([x for x in ens_ids if x]))
    print(f"[OK] Found {len(ens_ids)} unique ENSGWIG IDs in input.")

    # Batch lookup for Gouania ENS IDs
    gou_lookup: Dict[str, Any] = {}
    for ch in chunked(ens_ids, args.lookup_chunk):
        resp = http_post_json(ENSEMBL_LOOKUP_POST, {"ids": ch})
        gou_lookup.update(resp)
        time.sleep(args.sleep)

    # Collect NCBI GeneIDs from descriptions for genes lacking good Ensembl symbol
    ncbi_ids = []
    for gid in ens_ids:
        rec = gou_lookup.get(gid)
        if not isinstance(rec, dict):
            continue
        dn = (rec.get("display_name", "") or "").strip()
        desc = rec.get("description", "") or ""
        if not ensembl_good_symbol(dn):
            src, acc = parse_external_acc(desc)
            if src == "NCBI gene" and acc.isdigit():
                ncbi_ids.append(acc)

    ncbi_map = ncbi_geneid_to_symbol(
        ncbi_ids,
        chunk_size=args.ncbi_chunk,
        sleep=args.ncbi_sleep
    )
    print(f"[OK] Pulled NCBI symbols for {len(ncbi_map)} GeneIDs (may include LOC*).")

    # Determine which genes need zebrafish fallback:
    # Ensembl symbol missing/novel AND (NCBI symbol missing OR LOC*)
    need_zfish = []
    gou_external: Dict[str, Tuple[str, str, str]] = {}  # gid -> (src, acc, ncbi_symbol)

    for gid in ens_ids:
        rec = gou_lookup.get(gid)
        if not isinstance(rec, dict):
            gou_external[gid] = ("", "", "")
            need_zfish.append(gid)
            continue

        dn = (rec.get("display_name", "") or "").strip()
        desc = rec.get("description", "") or ""
        src, acc = parse_external_acc(desc)

        ncbi_symbol = ""
        if src == "NCBI gene" and acc.isdigit():
            ncbi_symbol = (ncbi_map.get(acc, "") or "").strip()

        gou_external[gid] = (src, acc, ncbi_symbol)

        if not ensembl_good_symbol(dn):
            if (not ncbi_symbol) or ncbi_symbol.upper().startswith("LOC"):
                need_zfish.append(gid)

    need_zfish = sorted(set(need_zfish))
    print(f"[OK] Will try zebrafish ortholog fallback for {len(need_zfish)} genes.")

    # Homology lookup per gene (only those that need fallback)
    gou_to_zf: Dict[str, Dict[str, Any]] = {}
    zf_ids_set = set()

    for gid in need_zfish:
        url = HOMOLOGY_URL_TMPL.format(
            species=args.species,
            gene_id=urllib.parse.quote(gid)
        )
        try:
            hjson = http_get_json(url)
        except Exception:
            gou_to_zf[gid] = {"zfish_gene_id": "", "orthology_type": "", "perc_id": -1.0}
            time.sleep(args.homology_sleep)
            continue

        homologies = []
        try:
            data = hjson.get("data", [])
            if data and isinstance(data, list):
                homologies = data[0].get("homologies", []) or []
        except Exception:
            homologies = []

        z_id, otype, perc = pick_best_ortholog(homologies)
        gou_to_zf[gid] = {"zfish_gene_id": z_id, "orthology_type": otype, "perc_id": perc}

        if z_id:
            zf_ids_set.add(z_id)

        time.sleep(args.homology_sleep)

    zf_ids = sorted(zf_ids_set)
    print(f"[OK] Found zebrafish ortholog IDs for {len(zf_ids)} unique zebrafish genes.")

    # Batch lookup zebrafish IDs for symbol/description
    zf_lookup: Dict[str, Any] = {}
    for ch in chunked(zf_ids, args.lookup_chunk):
        resp = http_post_json(ENSEMBL_LOOKUP_POST, {"ids": ch})
        zf_lookup.update(resp)
        time.sleep(args.sleep)

    # Build output rows for ALL inputs
    out_rows = []
    for r in base_rows:
        inp = r["input_gene"]
        gid = r["ens_id"]

        if not gid:
            out_rows.append({
                "input_gene": clean_field(inp),
                "ens_id": "",
                "short_id": clean_field(inp),
                "short_id_source": "input",
                "ensembl_display_name": "",
                "ensembl_description": "",
                "ensembl_biotype": "",
                "external_source": "",
                "external_acc": "",
                "ncbi_symbol": "",
                "zfish_gene_id": "",
                "zfish_symbol": "",
                "zfish_description": "",
                "orthology_type": "",
                "orthology_perc_id": "",
            })
            continue

        rec = gou_lookup.get(gid)
        if not isinstance(rec, dict):
            out_rows.append({
                "input_gene": clean_field(inp),
                "ens_id": gid,
                "short_id": gid,
                "short_id_source": "ens_id",
                "ensembl_display_name": "",
                "ensembl_description": "",
                "ensembl_biotype": "",
                "external_source": "",
                "external_acc": "",
                "ncbi_symbol": "",
                "zfish_gene_id": "",
                "zfish_symbol": "",
                "zfish_description": "",
                "orthology_type": "",
                "orthology_perc_id": "",
            })
            continue

        g_dn = (rec.get("display_name", "") or "").strip()
        g_desc = rec.get("description", "") or ""
        g_bio = rec.get("biotype", "") or ""

        ext_src, ext_acc, ncbi_symbol = gou_external.get(gid, ("", "", ""))
        ncbi_symbol = ncbi_symbol or ""

        zinfo = gou_to_zf.get(gid, {})
        z_id = zinfo.get("zfish_gene_id", "") or ""
        otype = zinfo.get("orthology_type", "") or ""
        perc = zinfo.get("perc_id", -1.0)
        try:
            perc_str = f"{float(perc):.2f}" if perc is not None and float(perc) >= 0 else ""
        except Exception:
            perc_str = ""

        z_dn = ""
        z_desc = ""
        if z_id and isinstance(zf_lookup.get(z_id), dict):
            z_dn = (zf_lookup[z_id].get("display_name", "") or "").strip()
            z_desc = zf_lookup[z_id].get("description", "") or ""

        # Decide short_id
        if ensembl_good_symbol(g_dn):
            short_id = g_dn
            short_src = "ensembl"
        else:
            if ncbi_symbol and (not ncbi_symbol.upper().startswith("LOC")):
                short_id = ncbi_symbol
                short_src = "ncbi"
            else:
                if z_dn:
                    short_id = z_dn
                    short_src = "zebrafish"
                else:
                    short_id = gid
                    short_src = "ens_id"

        out_rows.append({
            "input_gene": clean_field(inp),
            "ens_id": gid,
            "short_id": clean_field(short_id),
            "short_id_source": short_src,
            "ensembl_display_name": clean_field(g_dn),
            "ensembl_description": clean_field(g_desc),
            "ensembl_biotype": clean_field(g_bio),
            "external_source": clean_field(ext_src),
            "external_acc": clean_field(ext_acc),
            "ncbi_symbol": clean_field(ncbi_symbol),
            "zfish_gene_id": clean_field(z_id),
            "zfish_symbol": clean_field(z_dn),
            "zfish_description": clean_field(z_desc),
            "orthology_type": clean_field(otype),
            "orthology_perc_id": clean_field(perc_str),
        })

    df = pd.DataFrame(out_rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] Wrote {len(df)} rows to {args.out}")
    print(df.head(25).to_string(index=False))

if __name__ == "__main__":
    main()

