"""Gene roles (UniProt keywords) and population negatives (gnomAD v4 API).

Gene role
    UniProt keywords KW-0656 "Proto-oncogene" -> oncogene (gas pedal),
    KW-0043 "Tumor suppressor" -> tumour suppressor (brakes). CC BY 4.0.
    Coarse: some genes carry both or neither.

Negatives
    Missense variants seen in at least ``MIN_AF`` of gnomAD v4 exomes+genomes.
    A variant this common in the general population is very unlikely to be a
    somatic cancer driver. These are *population negatives*, not expert labels,
    and Koel always reports results with and without them.
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request

import pandas as pd

from .download import RAW, USER_AGENT

MIN_AF = 1e-3
GNOMAD_CACHE = RAW / "gnomad"
ROLE_CACHE = RAW / "uniprot" / "gene_roles.tsv"

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
       "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
       "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V"}
HGVSP = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def _get(url: str, data: bytes | None = None) -> dict:
    headers = {"User-Agent": USER_AGENT}
    if data is not None:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < 4:
                time.sleep(10 * (attempt + 1))
                continue
            raise
    raise RuntimeError("unreachable")


def gene_roles(accessions: list[str]) -> pd.DataFrame:
    if ROLE_CACHE.exists():
        cached = pd.read_csv(ROLE_CACHE, sep="\t")
        if set(accessions) <= set(cached["uniprot_id"]):
            return cached
    rows = []
    for i in range(0, len(accessions), 50):
        batch = accessions[i:i + 50]
        q = " OR ".join(f"accession:{a}" for a in batch)
        url = ("https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
            {"query": q, "fields": "accession,keyword", "format": "json", "size": 500}))
        for r in _get(url).get("results", []):
            kws = {k["name"] for k in r.get("keywords", [])}
            og, tsg = "Proto-oncogene" in kws, "Tumor suppressor" in kws
            role = "both" if og and tsg else "oncogene" if og else "tumour_suppressor" if tsg else "unknown"
            rows.append({"uniprot_id": r["primaryAccession"], "gene_role": role})
    df = pd.DataFrame(rows)
    ROLE_CACHE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(ROLE_CACHE, sep="\t", index=False)
    return df


GNOMAD_QUERY = """
query($g: String!) {
  gene(gene_symbol: $g, reference_genome: GRCh38) {
    canonical_transcript_id
    variants(dataset: gnomad_r4) {
      variant_id consequence hgvsp transcript_id
      exome { ac an } genome { ac an }
    }
  }
}"""


def gnomad_missense(gene: str) -> pd.DataFrame:
    GNOMAD_CACHE.mkdir(parents=True, exist_ok=True)
    path = GNOMAD_CACHE / f"{gene}.json"
    if path.exists():
        payload = json.loads(path.read_text())
    else:
        body = json.dumps({"query": GNOMAD_QUERY, "variables": {"g": gene}}).encode()
        payload = _get("https://gnomad.broadinstitute.org/api", body)
        path.write_text(json.dumps(payload))
        time.sleep(6)  # stay well inside gnomAD's public API rate limit
    g = (payload.get("data") or {}).get("gene") or {}
    rows = []
    for v in g.get("variants") or []:
        if v["consequence"] != "missense_variant" or not v.get("hgvsp"):
            continue
        m = HGVSP.match(v["hgvsp"])
        if not m or m.group(1) not in AA3 or m.group(3) not in AA3:
            continue
        ac = sum((v.get(k) or {}).get("ac", 0) for k in ("exome", "genome"))
        an = sum((v.get(k) or {}).get("an", 0) for k in ("exome", "genome"))
        if an == 0:
            continue
        rows.append({"gene": gene,
                     "variant": f"{AA3[m.group(1)]}{m.group(2)}{AA3[m.group(3)]}",
                     "gnomad_variant_id": v["variant_id"], "gnomad_af": ac / an})
    df = pd.DataFrame(rows, columns=["gene", "variant", "gnomad_variant_id", "gnomad_af"])
    return df.sort_values("gnomad_af", ascending=False).drop_duplicates(["gene", "variant"])


def gnomad_for_genes(genes: list[str]) -> pd.DataFrame:
    out = []
    for i, gene in enumerate(genes, 1):
        out.append(gnomad_missense(gene))
        if i % 10 == 0:
            print(f"  gnomAD: {i}/{len(genes)} genes")
    return pd.concat(out, ignore_index=True)


INTOGEN = RAW / "intogen" / "2024-06-18_IntOGen-Drivers" / "Compendium_Cancer_Genes.tsv"
INTOGEN_ROLE = {"Act": "oncogene", "LoF": "tumour_suppressor", "ambiguous": "both"}


def intogen_roles() -> pd.Series:
    """Gene symbol -> role from IntOGen's driver compendium (CC0).

    A gene appears once per cohort; take its most frequent role across cohorts.
    """
    c = pd.read_csv(INTOGEN, sep="\t", usecols=["SYMBOL", "ROLE"])
    r = c.groupby("SYMBOL")["ROLE"].agg(lambda s: s.value_counts().index[0])
    return r.map(INTOGEN_ROLE)


ESM1B_ZIP = RAW / "esm1b" / "ALL_hum_isoforms_ESM1b_LLR.zip"


def load_esm1b(accessions: set[str]) -> pd.DataFrame:
    """Precomputed ESM1b log-likelihood ratios for canonical isoforms (CC BY-NC 4.0).

    Each CSV is a matrix: rows = alternative amino acid, columns = "REF POS".
    Lower LLR means more damaging.
    """
    import zipfile

    rows = []
    with zipfile.ZipFile(ESM1B_ZIP) as z:
        names = {n.rsplit("/", 1)[-1].removesuffix("_LLR.csv"): n
                 for n in z.namelist() if n.endswith("_LLR.csv")}
        for acc in sorted(accessions):
            if acc not in names:
                continue
            with z.open(names[acc]) as f:
                m = pd.read_csv(f, index_col=0)
            long = m.stack().reset_index()
            long.columns = ["alt", "refpos", "esm1b_llr"]
            ref_pos = long["refpos"].str.split(" ", expand=True)
            long["variant"] = ref_pos[0] + ref_pos[1] + long["alt"]
            long = long[ref_pos[0] != long["alt"]]
            long["uniprot_id"] = acc
            rows.append(long[["uniprot_id", "variant", "esm1b_llr"]])
    return pd.concat(rows, ignore_index=True)


CLINVAR = RAW / "clinvar" / "variant_summary.txt.gz"
CLINVAR_PROTEIN = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
ONCO_POS = {"Oncogenic", "Likely oncogenic"}
ONCO_NEG = {"Benign", "Likely benign"}


def clinvar_oncogenicity(genes: set[str]) -> pd.DataFrame:
    """ClinVar somatic Oncogenicity classifications for missense variants (public domain).

    Oncogenic / Likely oncogenic -> 1; Benign / Likely benign -> 0; others dropped.
    Uses GRCh38 rows only so each variant is counted once.
    """
    cols = ["Type", "Name", "GeneSymbol", "Assembly", "Oncogenicity",
            "OncogenicityLastEvaluated", "ReviewStatusOncogenicity", "VariationID"]
    df = pd.read_csv(CLINVAR, sep="\t", usecols=cols, dtype=str, low_memory=False)
    df = df[(df["Assembly"] == "GRCh38") & (df["Type"] == "single nucleotide variant")
            & df["GeneSymbol"].isin(genes)
            & df["Oncogenicity"].isin(ONCO_POS | ONCO_NEG)]
    p = df["Name"].str.extract(CLINVAR_PROTEIN)
    df = df[p[0].isin(AA3) & p[2].isin(AA3)].copy()
    p = p.loc[df.index]
    df["variant"] = p[0].map(AA3) + p[1] + p[2].map(AA3)
    df["clinvar_oncogenic"] = df["Oncogenicity"].isin(ONCO_POS).astype(int)
    df = df.rename(columns={"GeneSymbol": "gene", "VariationID": "clinvar_variation_id",
                            "OncogenicityLastEvaluated": "clinvar_last_evaluated",
                            "ReviewStatusOncogenicity": "clinvar_review_status",
                            "Oncogenicity": "clinvar_oncogenicity"})
    keep = ["gene", "variant", "clinvar_oncogenic", "clinvar_oncogenicity",
            "clinvar_last_evaluated", "clinvar_review_status", "clinvar_variation_id"]
    df = df[keep]
    # Drop variants with conflicting ClinVar calls across records.
    n = df.groupby(["gene", "variant"])["clinvar_oncogenic"].transform("nunique")
    return df[n == 1].drop_duplicates(["gene", "variant"])
