"""First look: how well does AlphaMissense separate labels, per gene class?

Prints AUROC with a bootstrap 95% CI for each gene role, under two label sets:
  * ``expert``: CIViC labels only (small, but no population-frequency leakage).
  * ``expert+population``: adds common gnomAD missense as negatives.

Caveat built into the output: AlphaMissense was trained partly on the idea that
variants common in human and primate populations are benign, so population
negatives flatter it. The expert-only numbers are the honest ones.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from .table import PROCESSED

ROLES = ["oncogene", "tumour_suppressor", "both", "unknown", "all"]


def _boot_auc(y: np.ndarray, s: np.ndarray, n: int = 2000, seed: int = 7) -> tuple[float, float, float]:
    auc = roc_auc_score(y, s)
    rng = np.random.default_rng(seed)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    boots = []
    for _ in range(n):
        idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        boots.append(roc_auc_score(y[idx], s[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return auc, lo, hi


def summarise() -> pd.DataFrame:
    d = pd.read_parquet(PROCESSED / "civic_missense_alphamissense.parquet")
    d = d[d["am_pathogenicity"].notna() & d["label"].notna()].copy()
    d["label"] = d["label"].astype(int)
    rows = []
    for labelset, frame in {
        "expert": d[d["label_source"] != "gnomad_common"],
        "expert+population": d,
    }.items():
        for role in ROLES:
            sub = frame if role == "all" else frame[frame["gene_role"] == role]
            n_pos, n_neg = int((sub["label"] == 1).sum()), int((sub["label"] == 0).sum())
            row = {"labels": labelset, "gene_role": role, "n_pos": n_pos, "n_neg": n_neg,
                   "auroc": np.nan, "ci_low": np.nan, "ci_high": np.nan}
            if n_pos >= 5 and n_neg >= 5:
                row["auroc"], row["ci_low"], row["ci_high"] = _boot_auc(
                    sub["label"].to_numpy(), sub["am_pathogenicity"].to_numpy())
            rows.append(row)
    res = pd.DataFrame(rows)
    out = PROCESSED / "first_look_auroc.csv"
    res.to_csv(out, index=False)
    with pd.option_context("display.float_format", "{:.3f}".format):
        print(res.to_string(index=False))
    print("\nAUROC shown only where both classes have >= 5 variants.")
    print("Population negatives flatter AlphaMissense (trained with population data as benign).")
    return res
