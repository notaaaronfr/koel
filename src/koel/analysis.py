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
# Predictor name -> (column, sign). Higher signed score = more likely damaging.
PREDICTORS = {"AlphaMissense": ("am_pathogenicity", 1.0), "ESM1b": ("esm1b_llr", -1.0)}


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
    # Compare predictors on the same variants: keep rows scored by all of them.
    cols = [c for c, _ in PREDICTORS.values()]
    d = d[d[cols].notna().all(axis=1)].copy()
    rows = []
    act = d[d["label"].isna() & (d["civic_actionable"] == 1)].copy()
    act["label"] = 1
    d = d[d["label"].notna()].copy()
    d["label"] = d["label"].astype(int)
    expert = d[d["label_source"] != "gnomad_common"]
    for labelset, frame in {
        "expert (CIViC+ClinVar)": expert,
        "expert+actionable": pd.concat([expert, act]),
        "expert+population": d,
    }.items():
        for role in ROLES:
            sub = frame if role == "all" else frame[frame["gene_role"] == role]
            n_pos, n_neg = int((sub["label"] == 1).sum()), int((sub["label"] == 0).sum())
            for name, (col, sign) in PREDICTORS.items():
                row = {"labels": labelset, "gene_role": role, "predictor": name,
                       "n_pos": n_pos, "n_neg": n_neg,
                       "auroc": np.nan, "ci_low": np.nan, "ci_high": np.nan}
                if n_pos >= 5 and n_neg >= 5:
                    row["auroc"], row["ci_low"], row["ci_high"] = _boot_auc(
                        sub["label"].to_numpy(), sign * sub[col].to_numpy())
                rows.append(row)
    res = pd.DataFrame(rows)
    out = PROCESSED / "first_look_auroc.csv"
    res.to_csv(out, index=False)
    with pd.option_context("display.float_format", "{:.3f}".format):
        print(res.to_string(index=False))
    print("\nAUROC shown only where both classes have >= 5 variants.")
    print("Population negatives flatter AlphaMissense (trained with population data as benign).")
    return res


def lab_tested() -> pd.DataFrame:
    """Compare predictors with lab-measured function (MaveDB, CC0).

    Brakes (TP53, PTEN, BRCA1): damaging vs normal.
    Gas pedals (KRAS, BRAF, EGFR, ERBB2, MAP2K1): three questions:
      * activating vs everything else (the oncogenicity question);
      * activating vs neutral (does the tool notice switched-on variants at all?);
      * activating vs loss (can it tell switched ON from switched OFF?).
    """
    from scipy.stats import spearmanr

    from .enrich import load_esm1b
    from .mave import all_assays
    from .table import gene_to_uniprot, load_alphamissense

    mave, stats = all_assays()
    mapping = gene_to_uniprot(sorted(mave["gene"].unique()))
    mave["uniprot_id"] = mave["gene"].map(mapping)
    acc = set(mave["uniprot_id"])
    am = load_alphamissense(acc).rename(columns={"protein_variant": "variant"})
    d = (mave.merge(am, on=["uniprot_id", "variant"], how="inner")
             .merge(load_esm1b(acc), on=["uniprot_id", "variant"], how="inner"))
    d.to_parquet(PROCESSED / "lab_tested.parquet", index=False)

    for gene, var in [("TP53", "R175H"), ("PTEN", "C124S"), ("BRCA1", "C61G"),
                      ("KRAS", "G12D"), ("BRAF", "V600E"), ("EGFR", "L858R"),
                      ("ERBB2", "L755S"), ("MAP2K1", "Q56P")]:
        r = d[(d.gene == gene) & (d.variant == var)]
        print(f"check {gene} {var}: lab effect = {r['mave_effect'].tolist()}")

    comparisons = {
        "tumour_suppressor": [("damaging vs normal", ["damaging"], ["normal"])],
        "oncogene": [("activating vs rest", ["activating"], ["neutral", "loss"]),
                     ("activating vs neutral", ["activating"], ["neutral"]),
                     ("activating vs loss", ["activating"], ["loss"])],
    }
    rows = []
    groups = [(g, sub) for g, sub in d.groupby("gene")] + \
             [(f"ALL {role}", sub) for role, sub in d.groupby("gene_role")]
    for gene, g in groups:
        role = g["gene_role"].iloc[0]
        for comp, pos_e, neg_e in comparisons[role]:
            sub = g[g["mave_effect"].isin(pos_e + neg_e)]
            y = sub["mave_effect"].isin(pos_e).astype(int).to_numpy()
            if y.sum() < 5 or (y == 0).sum() < 5:
                continue
            for name, (col, sign) in PREDICTORS.items():
                auc, lo, hi = _boot_auc(y, sign * sub[col].to_numpy(), n=500)
                rho = (spearmanr(sign * g[col], g["mave_score_signed"]).statistic
                       if role == "tumour_suppressor" else np.nan)
                rows.append({"gene": gene, "gene_role": role, "comparison": comp,
                             "predictor": name, "n_pos": int(y.sum()),
                             "n_neg": int((y == 0).sum()), "auroc": auc,
                             "ci_low": lo, "ci_high": hi, "spearman_vs_lab": rho})
    res = pd.DataFrame(rows)
    res.to_csv(PROCESSED / "lab_tested_auroc.csv", index=False)
    stats.to_csv(PROCESSED / "lab_tested_assays.csv", index=False)
    with pd.option_context("display.float_format", "{:.3f}".format, "display.width", 200):
        print(stats[["gene", "assay", "n_missense", "n_damaging", "n_normal"]
                    + (["n_loss"] if "n_loss" in stats else [])].to_string(index=False))
        print(res.to_string(index=False))
    print("\nLabels are functional, not oncogenicity calls. AUROC 0.5 = no better than chance.")
    print("'activating vs loss' below 0.5 means the tool scores switched-OFF variants as MORE"
          " harmful than switched-ON ones.")
    return res


AM_RELEASE = pd.Timestamp("2023-09-19")
# Thresholds: AlphaMissense developer "likely pathogenic" (0.564) and the stricter
# ClinGen-calibrated cut-off reported for germline PP3 (0.792); ESM1b LLR -7.5
# (Brandes et al. 2023). Germline calibrations, used here only as reference points.
THRESHOLDS = [("AlphaMissense", "am_pathogenicity", ">=", 0.564),
              ("AlphaMissense", "am_pathogenicity", ">=", 0.792),
              ("ESM1b", "esm1b_llr", "<=", -7.5)]


def _wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def time_split() -> pd.DataFrame:
    """Held-out test: oncogenic calls first evaluated after AlphaMissense's release.

    Test positives: ClinVar Oncogenic/Likely oncogenic, last evaluated on or after
    19 Sep 2023, with assertion criteria, and NOT already in CIViC (so the tools
    cannot have been tuned on them). Reported per gene class:
      * detection rate at fixed thresholds (Wilson 95% CI);
      * false-flag rate on common gnomAD variants in the same genes;
      * AUROC against those gnomAD negatives (secondary; flatters AlphaMissense).
    """
    d = pd.read_parquet(PROCESSED / "civic_missense_alphamissense.parquet")
    d = d[d["am_pathogenicity"].notna() & d["esm1b_llr"].notna()].copy()
    when = pd.to_datetime(d["clinvar_last_evaluated"], errors="coerce", format="mixed")
    test = d[(d["clinvar_oncogenic"] == 1) & (when >= AM_RELEASE)
             & (d["source"] == "clinvar")
             & d["clinvar_review_status"].str.startswith("criteria provided", na=False)]
    neg = d[d["label_source"] == "gnomad_common"]
    test.to_parquet(PROCESSED / "time_split_test.parquet", index=False)
    print(f"Held-out oncogenic variants: {len(test)} "
          f"({test['gene'].nunique()} genes); gnomAD negatives: {len(neg)}")

    rows = []
    for role in ["oncogene", "tumour_suppressor", "all"]:
        pos_r = test if role == "all" else test[test["gene_role"] == role]
        neg_r = neg if role == "all" else neg[neg["gene_role"] == role]
        for name, col, op, thr in THRESHOLDS:
            hit_p = pos_r[col] >= thr if op == ">=" else pos_r[col] <= thr
            hit_n = neg_r[col] >= thr if op == ">=" else neg_r[col] <= thr
            k, n = int(hit_p.sum()), len(pos_r)
            lo, hi = _wilson(k, n)
            fk, fn = int(hit_n.sum()), len(neg_r)
            rows.append({"gene_role": role, "predictor": name, "threshold": f"{op}{thr}",
                         "n_test": n, "detected": k, "detection_rate": k / n if n else np.nan,
                         "det_ci_low": lo, "det_ci_high": hi,
                         "n_gnomad": fn, "false_flag_rate": fk / fn if fn else np.nan})
        if len(pos_r) >= 5 and len(neg_r) >= 5:
            y = np.r_[np.ones(len(pos_r)), np.zeros(len(neg_r))].astype(int)
            for name, (col, sign) in PREDICTORS.items():
                s = sign * np.r_[pos_r[col].to_numpy(), neg_r[col].to_numpy()]
                auc, lo, hi = _boot_auc(y, s)
                rows.append({"gene_role": role, "predictor": name, "threshold": "AUROC vs gnomAD",
                             "n_test": len(pos_r), "detection_rate": auc,
                             "det_ci_low": lo, "det_ci_high": hi, "n_gnomad": len(neg_r)})
    res = pd.DataFrame(rows)
    res.to_csv(PROCESSED / "time_split_results.csv", index=False)
    with pd.option_context("display.float_format", "{:.3f}".format, "display.width", 200):
        print(res.to_string(index=False))
    print("\ndetection_rate = share of held-out oncogenic variants the tool flags.")
    print("false_flag_rate = share of common (presumed harmless) gnomAD variants it flags.")
    return res


def _gene_boot(frames: dict, stat, n: int = 2000, seed: int = 3) -> tuple[float, float]:
    """95% CI resampling whole genes (clusters), not variants."""
    rng = np.random.default_rng(seed)
    keys = list(frames)
    out = []
    for _ in range(n):
        d = pd.concat([frames[keys[i]] for i in rng.integers(0, len(keys), len(keys))])
        v = stat(d)
        if v is not None and not np.isnan(v):
            out.append(v)
    lo, hi = np.percentile(out, [2.5, 97.5])
    return float(lo), float(hi)


def robustness() -> pd.DataFrame:
    """Gene-level (cluster) bootstrap for the headline claims.

    Variant-level intervals treat mutations in the same gene as independent and look
    sharper than the data supports (point raised by S.-G. Lee). Here whole genes are
    resampled. Run after ``timesplit`` and ``lab``.
    """
    rows = []
    t = pd.read_parquet(PROCESSED / "time_split_test.parquet")
    t = t[t["gene_role"].isin(["oncogene", "tumour_suppressor"])].copy()
    t["hit"] = t["am_pathogenicity"] >= 0.792
    rate = {r: t.loc[t["gene_role"] == r, "hit"].mean() for r in ["oncogene", "tumour_suppressor"]}
    by = {r: t[t["gene_role"] == r].groupby("gene")["hit"].agg(["sum", "count"])
          for r in ["oncogene", "tumour_suppressor"]}
    rng = np.random.default_rng(3)
    diffs = []
    for _ in range(5000):
        s = {r: g.iloc[rng.integers(0, len(g), len(g))].sum() for r, g in by.items()}
        diffs.append(s["oncogene"]["sum"] / s["oncogene"]["count"]
                     - s["tumour_suppressor"]["sum"] / s["tumour_suppressor"]["count"])
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    rows.append({"claim": "time split: detection difference at AM>=0.792 (oncogene - TSG)",
                 "estimate": rate["oncogene"] - rate["tumour_suppressor"], "gene_ci_low": lo,
                 "gene_ci_high": hi, "n_genes": f"{len(by['oncogene'])} vs {len(by['tumour_suppressor'])}"})

    lab = pd.read_parquet(PROCESSED / "lab_tested.parquet")
    for claim, role, pos, neg in [
        ("lab: AUROC activating vs loss (oncogenes)", "oncogene", ["activating"], ["loss"]),
        ("lab: AUROC activating vs rest (oncogenes)", "oncogene", ["activating"], ["neutral", "loss"]),
        ("lab: AUROC damaging vs normal (TSGs)", "tumour_suppressor", ["damaging"], ["normal"]),
    ]:
        d0 = lab[(lab["gene_role"] == role) & lab["mave_effect"].isin(pos + neg)]
        frames = dict(tuple(d0.groupby("gene")))

        def stat(d, pos=pos):
            y = d["mave_effect"].isin(pos).astype(int)
            return roc_auc_score(y, d["am_pathogenicity"]) if y.nunique() == 2 else None

        lo, hi = _gene_boot(frames, stat, n=1000)
        rows.append({"claim": claim, "estimate": stat(d0), "gene_ci_low": lo,
                     "gene_ci_high": hi, "n_genes": str(len(frames))})
    res = pd.DataFrame(rows)
    res.to_csv(PROCESSED / "robustness_gene_level.csv", index=False)
    with pd.option_context("display.float_format", "{:.3f}".format, "display.width", 200):
        print(res.to_string(index=False))
    print("\nWhole genes resampled. With 3-5 genes per lab group, read the intervals, not the point estimates.")
    return res
