"""Build Koel's first labelled table: CIViC missense variants joined to AlphaMissense.

Output: data/processed/civic_missense_alphamissense.parquet (+ a CSV preview).

Labels (v0, deliberately conservative and explicit):
  * ``civic_oncogenic``: an accepted CIViC "Oncogenic" evidence item that Supports
    oncogenicity -> positive; one that Does Not Support -> negative.
  * ``civic_actionable``: an accepted Predictive/Diagnostic/Prognostic item at level
    A or B that Supports its claim. A *proxy* for a functional driver, not proof.
  * ``label``: 1 = CIViC oncogenic; 0 = CIViC not oncogenic, or a common gnomAD
    missense variant (population negative, see enrich.py).
Gene role comes from UniProt keywords (enrich.py).
"""

from __future__ import annotations

import gzip
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from .download import RAW, ROOT, USER_AGENT

PROCESSED = ROOT / "data" / "processed"
MISSENSE = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")
UNIPROT_CACHE = RAW / "uniprot" / "gene_to_uniprot.tsv"


def load_civic_missense() -> pd.DataFrame:
    v = pd.read_csv(RAW / "civic/VariantSummaries.tsv", sep="\t", low_memory=False)
    v = v[v["feature_type"] == "Gene"].copy()
    v["variant"] = v["variant"].astype(str).str.strip()
    v = v[v["variant"].str.match(MISSENSE)]
    return v[["variant_id", "gene", "variant", "single_variant_molecular_profile_id",
              "representative_transcript", "variant_civic_url"]]


def civic_labels() -> pd.DataFrame:
    e = pd.read_csv(RAW / "civic/ClinicalEvidenceSummaries.tsv", sep="\t", low_memory=False)
    e = e[e["evidence_status"] == "accepted"]
    onc = e[e["evidence_type"] == "Oncogenic"]
    act = e[e["evidence_type"].isin(["Predictive", "Diagnostic", "Prognostic"])
            & e["evidence_level"].isin(["A", "B"])
            & (e["evidence_direction"] == "Supports")]

    def agg(df: pd.DataFrame, name: str) -> pd.DataFrame:
        g = df.groupby("molecular_profile_id")
        return pd.DataFrame({
            f"{name}_n": g.size(),
            f"{name}_evidence_ids": g["evidence_id"].apply(lambda s: ";".join(map(str, sorted(s)))),
        })

    sup = agg(onc[onc["evidence_direction"] == "Supports"], "onc_supports")
    neg = agg(onc[onc["evidence_direction"] == "Does Not Support"], "onc_against")
    a = agg(act, "actionable")
    return sup.join(neg, how="outer").join(a, how="outer")


def gene_to_uniprot(genes: list[str]) -> dict[str, str]:
    cache: dict[str, str] = {}
    if UNIPROT_CACHE.exists():
        c = pd.read_csv(UNIPROT_CACHE, sep="\t", dtype=str).fillna("")
        cache = dict(zip(c["gene"], c["uniprot_id"]))
    missing = [g for g in genes if g not in cache]
    for i, gene in enumerate(missing, 1):
        q = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
        url = ("https://rest.uniprot.org/uniprotkb/search?"
               + urllib.parse.urlencode({"query": q, "fields": "accession", "format": "json", "size": 1}))
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=30) as r:
            res = json.load(r).get("results", [])
        cache[gene] = res[0]["primaryAccession"] if res else ""
        if i % 25 == 0:
            print(f"  uniprot: {i}/{len(missing)}")
        time.sleep(0.1)
    UNIPROT_CACHE.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(sorted(cache.items()), columns=["gene", "uniprot_id"]).to_csv(
        UNIPROT_CACHE, sep="\t", index=False)
    return cache


def load_alphamissense(uniprot_ids: set[str]) -> pd.DataFrame:
    path = RAW / "alphamissense/AlphaMissense_aa_substitutions.tsv.gz"
    chunks = []
    with gzip.open(path, "rt") as f:
        reader = pd.read_csv(f, sep="\t", comment="#", chunksize=2_000_000)
        for ch in reader:
            chunks.append(ch[ch["uniprot_id"].isin(uniprot_ids)])
    return pd.concat(chunks, ignore_index=True)


def build() -> Path:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    v = load_civic_missense()
    print(f"CIViC missense variants: {len(v)} across {v['gene'].nunique()} genes")

    labels = civic_labels()
    v = v.merge(labels, left_on="single_variant_molecular_profile_id",
                right_index=True, how="left")
    for col in ("onc_supports_n", "onc_against_n", "actionable_n"):
        v[col] = v[col].fillna(0).astype(int)
    v["civic_oncogenic"] = pd.NA
    v.loc[v["onc_supports_n"] > 0, "civic_oncogenic"] = 1
    v.loc[(v["onc_against_n"] > 0) & (v["onc_supports_n"] == 0), "civic_oncogenic"] = 0
    v["civic_actionable"] = (v["actionable_n"] > 0).astype(int)
    v["source"] = "civic"
    civic_genes = sorted(v["gene"].unique())

    # ClinVar somatic oncogenicity, in CIViC genes and IntOGen driver genes.
    from .enrich import clinvar_oncogenicity, intogen_roles
    cv = clinvar_oncogenicity(set(civic_genes) | set(intogen_roles().index))
    print(f"ClinVar oncogenicity missense: {len(cv)} "
          f"({cv['clinvar_oncogenic'].sum()} oncogenic) across {cv['gene'].nunique()} genes")
    v = v.merge(cv, on=["gene", "variant"], how="outer")
    v["source"] = v["source"].fillna("clinvar")
    for col in ("onc_supports_n", "onc_against_n", "actionable_n", "civic_actionable"):
        v[col] = v[col].fillna(0).astype(int)

    mapping = gene_to_uniprot(sorted(v["gene"].unique()))
    v["uniprot_id"] = v["gene"].map(mapping).fillna("")

    # Population negatives: common gnomAD missense in the CIViC genes.
    from .enrich import MIN_AF, gene_roles, gnomad_for_genes
    gn = gnomad_for_genes(civic_genes)
    common = gn[gn["gnomad_af"] >= MIN_AF].copy()
    civic_keys = set(zip(v["gene"], v["variant"]))  # CIViC + ClinVar
    common = common[[k not in civic_keys for k in zip(common["gene"], common["variant"])]]
    common["uniprot_id"] = common["gene"].map(mapping).fillna("")
    common["source"] = "gnomad_common"
    v = v.merge(gn[["gene", "variant", "gnomad_af"]], on=["gene", "variant"], how="left")
    v = pd.concat([v, common], ignore_index=True)
    print(f"gnomAD common missense negatives (AF >= {MIN_AF}): {len(common)}")

    # Unified label: 1 = oncogenic, 0 = not oncogenic, NA = unknown.
    v["label"] = pd.NA
    v["label_source"] = pd.NA
    pos = v["civic_oncogenic"] == 1
    v.loc[pos, ["label", "label_source"]] = [1, "civic_oncogenic"]
    neg = v["civic_oncogenic"] == 0
    v.loc[neg, ["label", "label_source"]] = [0, "civic_not_oncogenic"]
    # ClinVar fills in where CIViC has no oncogenicity call; disagreements are flagged.
    cvl = v["clinvar_oncogenic"]
    fill = v["label"].isna() & cvl.notna()
    v.loc[fill & (cvl == 1), ["label", "label_source"]] = [1, "clinvar_oncogenic"]
    v.loc[fill & (cvl == 0), ["label", "label_source"]] = [0, "clinvar_benign"]
    v["label_conflict"] = v["civic_oncogenic"].notna() & cvl.notna() & (
        v["civic_oncogenic"].astype("Float64") != cvl.astype("Float64"))
    gneg = v["source"] == "gnomad_common"
    v.loc[gneg, ["label", "label_source"]] = [0, "gnomad_common"]

    roles = gene_roles(sorted(set(filter(None, mapping.values()))))
    v = v.merge(roles, on="uniprot_id", how="left")
    v["gene_role_uniprot"] = v["gene_role"].fillna("unknown")
    # Primary role: IntOGen (CC0); fall back to UniProt keywords.
    from .enrich import intogen_roles
    v["gene_role_intogen"] = v["gene"].map(intogen_roles())
    v["gene_role"] = v["gene_role_intogen"].fillna(v["gene_role_uniprot"])
    v["gene_role_source"] = v["gene_role_intogen"].notna().map({True: "intogen", False: "uniprot"})

    am = load_alphamissense(set(filter(None, mapping.values())))
    out = v.merge(am, left_on=["uniprot_id", "variant"],
                  right_on=["uniprot_id", "protein_variant"], how="left")
    out = out.drop(columns=["protein_variant"])

    from .enrich import load_esm1b
    esm = load_esm1b(set(filter(None, mapping.values())))
    out = out.merge(esm, on=["uniprot_id", "variant"], how="left")
    print(f"Variants with an ESM1b score: {out['esm1b_llr'].notna().sum()}/{len(out)}")

    path = PROCESSED / "civic_missense_alphamissense.parquet"
    out.to_parquet(path, index=False)
    out.head(200).to_csv(PROCESSED / "civic_missense_alphamissense.preview.csv", index=False)

    matched = out["am_pathogenicity"].notna()
    print(f"UniProt mapped genes: {sum(1 for x in mapping.values() if x)}/{len(mapping)}")
    print(f"Variants with an AlphaMissense score: {matched.sum()}/{len(out)}")
    print(out.groupby(["gene_role", "label_source"], dropna=False).size().to_string())
    print(f"CIViC/ClinVar label conflicts: {int(out['label_conflict'].sum())}")
    print(f"Wrote {path.relative_to(ROOT)}")
    return path
