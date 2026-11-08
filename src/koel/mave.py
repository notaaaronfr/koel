"""Lab-measured variant effects (deep mutational scans / saturation genome editing).

Datasets come from MaveDB and are CC0. For each assay Koel:
  1. converts every entry to a single amino-acid change (e.g. R175H);
  2. uses synonymous variants as the "normal" reference and nonsense (stop)
     variants as the "broken" reference, so assay direction is learned from the
     data rather than assumed;
  3. labels a missense variant
       * 1 (functionally damaging) if it lies past the midpoint between the
         synonymous and nonsense medians AND outside the synonymous 95% range;
       * 0 (functionally normal) if it lies inside the synonymous 5-95% range;
       * NA otherwise (intermediate).
These are *functional* labels, not oncogenicity calls, and are reported separately.
"""

from __future__ import annotations

import io
import json
import re
import urllib.request

import numpy as np
import pandas as pd

from .download import RAW, USER_AGENT
from .enrich import AA3

MAVE_DIR = RAW / "mavedb"

ASSAYS = [
    {"gene": "TP53", "urn": "urn:mavedb:00001213-a-1", "protein_col": "HGVS(protein)",
     "label": "TP53 DMS in HCT116 colorectal cancer cells"},
    {"gene": "PTEN", "urn": "urn:mavedb:00000054-a-1", "protein_col": "hgvs_pro",
     "label": "PTEN lipid phosphatase activity (Mighell 2018)"},
    {"gene": "BRCA1", "urn": "urn:mavedb:00000097-0-2", "protein_col": None,
     "transcript": "NM_007294.3", "label": "BRCA1 saturation genome editing (Findlay 2018)"},
]

PRO = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|Ter|\*)\)?$")
CODON = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_bases):
    for j, b in enumerate(_bases):
        for k, c in enumerate(_bases):
            CODON[a + b + c] = _aas[16 * i + 4 * j + k]
ONE_TO_THREE = {v: k for k, v in AA3.items()}


def _fetch(url: str, path) -> str:
    if path.exists():
        return path.read_text()
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:
        text = r.read().decode()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return text


def _scores(urn: str) -> pd.DataFrame:
    safe = urn.replace(":", "_")
    meta = json.loads(_fetch(f"https://api.mavedb.org/api/v1/score-sets/{urn}", MAVE_DIR / f"{safe}.json"))
    lic = (meta.get("license") or {}).get("shortName")
    if lic != "CC0":
        raise ValueError(f"{urn} licence is {lic}, expected CC0")
    text = _fetch(f"https://api.mavedb.org/api/v1/score-sets/{urn}/scores", MAVE_DIR / f"{safe}.csv")
    return pd.read_csv(io.StringIO(text))


def _cds(transcript: str) -> str:
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
           f"&id={transcript}&rettype=fasta_cds_na&retmode=text")
    fasta = _fetch(url, MAVE_DIR / f"{transcript}.cds.fa")
    return "".join(line.strip() for line in fasta.splitlines() if not line.startswith(">")).upper()


def _protein_from_nt(df: pd.DataFrame, transcript: str) -> pd.Series:
    cds = _cds(transcript)
    nt = re.compile(r"c\.(\d+)([ACGT])>([ACGT])$")
    out = []
    for h in df["hgvs_nt"].astype(str):
        m = nt.search(h)
        if not m:
            out.append(None)
            continue
        pos, ref, alt = int(m.group(1)), m.group(2), m.group(3)
        if pos > len(cds) or cds[pos - 1] != ref:
            out.append(None)
            continue
        ci = (pos - 1) // 3
        codon = cds[3 * ci:3 * ci + 3]
        mut = codon[:(pos - 1) % 3] + alt + codon[(pos - 1) % 3 + 1:]
        a, b = CODON[codon], CODON[mut]
        if a == "*":
            out.append(None)
        elif b == "*":
            out.append(f"p.{ONE_TO_THREE[a]}{ci + 1}Ter")
        elif a == b:
            out.append(f"p.{ONE_TO_THREE[a]}{ci + 1}=")
        else:
            out.append(f"p.{ONE_TO_THREE[a]}{ci + 1}{ONE_TO_THREE[b]}")
    return pd.Series(out, index=df.index)


def _classify(pro: pd.Series) -> pd.DataFrame:
    rows = []
    for h in pro:
        if not isinstance(h, str):
            rows.append((None, None))
            continue
        m = PRO.search(h.split(":")[-1])
        if not m or m.group(1) not in AA3:
            rows.append((None, None))
            continue
        ref, pos, alt = m.group(1), m.group(2), m.group(3)
        if alt == "=":
            rows.append(("synonymous", None))
        elif alt in ("Ter", "*"):
            rows.append(("nonsense", None))
        elif alt in AA3:
            rows.append(("missense", f"{AA3[ref]}{pos}{AA3[alt]}"))
        else:
            rows.append((None, None))
    return pd.DataFrame(rows, columns=["kind", "variant"], index=pro.index)


def assay_labels(assay: dict) -> tuple[pd.DataFrame, dict]:
    df = _scores(assay["urn"])
    pro = (_protein_from_nt(df, assay["transcript"]) if assay["protein_col"] is None
           else df[assay["protein_col"]])
    df = pd.concat([df[["score"]], _classify(pro)], axis=1).dropna(subset=["score", "kind"])
    syn, non = df[df.kind == "synonymous"]["score"], df[df.kind == "nonsense"]["score"]
    mis = df[df.kind == "missense"].copy()
    n_med = non.median()
    if len(syn) >= 20:
        method = "synonymous/nonsense reference"
        s_med = syn.median()
        sign = 1.0 if n_med > s_med else -1.0       # direction of "broken"
        zm = sign * (mis["score"] - s_med)           # larger = more broken
        lo, hi = np.percentile(sign * (syn - s_med), [5, 95])
        mid = abs(n_med - s_med) / 2
        damaging = (zm >= mid) & (zm > hi)
        normal = (zm >= lo) & (zm <= hi)
    else:
        # No synonymous controls: two-component Gaussian mixture on missense scores;
        # the component nearer the nonsense median is "broken".
        from sklearn.mixture import GaussianMixture
        method = "2-component mixture (no synonymous controls)"
        x = mis[["score"]].to_numpy()
        gm = GaussianMixture(2, random_state=0).fit(x)
        broken = int(np.argmin(np.abs(gm.means_.ravel() - n_med)))
        post = gm.predict_proba(x)[:, broken]
        s_med = float(gm.means_.ravel()[1 - broken])
        sign = 1.0 if n_med > s_med else -1.0
        zm = pd.Series(sign * (mis["score"].to_numpy() - s_med), index=mis.index)
        damaging = pd.Series(post >= 0.95, index=mis.index)
        normal = pd.Series(post <= 0.05, index=mis.index)
    mis["mave_label"] = pd.NA
    mis.loc[damaging, "mave_label"] = 1
    mis.loc[normal, "mave_label"] = 0
    mis["mave_score_signed"] = zm
    mis = (mis.groupby("variant", as_index=False)
              .agg(mave_score_signed=("mave_score_signed", "mean"),
                   mave_label=("mave_label", lambda s: s.dropna().iloc[0] if s.dropna().nunique() == 1 else pd.NA)))
    mis["gene"], mis["mave_assay"] = assay["gene"], assay["urn"]
    stats = {"gene": assay["gene"], "urn": assay["urn"], "assay": assay["label"],
             "n_synonymous": len(syn), "n_nonsense": len(non), "n_missense": len(mis),
             "n_damaging": int((mis["mave_label"] == 1).sum()),
             "n_normal": int((mis["mave_label"] == 0).sum()),
             "broken_direction": "higher" if sign > 0 else "lower", "method": method}
    return mis, stats


def all_assays() -> tuple[pd.DataFrame, pd.DataFrame]:
    frames, stats = [], []
    for a in ASSAYS:
        m, s = assay_labels(a)
        frames.append(m)
        stats.append(s)
    return pd.concat(frames, ignore_index=True), pd.DataFrame(stats)
