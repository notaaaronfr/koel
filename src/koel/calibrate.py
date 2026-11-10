"""Trust table: how much evidence is a predictor score worth, per gene type?

ClinGen/CGC/VICC somatic oncogenicity SOP (Horak et al. 2022) lets a computational
prediction count as evidence (OP1, +1 point; SBP1, -1 point). It uses one fixed
strength for every gene. Koel asks: what strength does the data actually support,
separately for oncogenes (gas pedals) and tumour suppressors (brakes)?

Method (after Tavtigian 2018 / Pejaver 2022, v0):
  * Local likelihood ratio at score t: among variants scoring near t (a window holding
    ~10% of labelled variants), LR = P(near t | oncogenic) / P(near t | benign).
    Local, not "at or past t": a cumulative LR lets very high scores lend evidence to
    middling ones, which the first version of this module got wrong.
  * Bootstrap (2,000 resamples) of whole genes, not single variants. Variants in one
    gene are not independent, and a few genes carry many of the labels, so resampling
    variants gives intervals that are too narrow. A threshold earns a strength only if
    the one-sided 95% lower bound of the local LR clears the cut-off at t AND at every
    higher score (benign side: upper bound below the reciprocal at t and every lower
    score). The variant-level threshold is kept in ``threshold_variant_level`` so the
    cost of the honest interval is visible.
  * Cut-offs (Tavtigian 2018 OddsPath): supporting 2.08, moderate 4.33, strong 18.7,
    very strong 350. Points as in the somatic SOP: 1, 2, 4, 8.
  * The threshold reported is the loosest one that reaches the strength.

Caveats printed with every run:
  * With only 3-5 genes (the lab sets) a gene bootstrap is coarse: read those rows as
    "what survives if any one gene is dropped or doubled", not as a precise interval.
  * Windowed LRs are a simplification of Pejaver's kernel-smoothed local posteriors.
  * Negatives are mostly common gnomAD variants, which AlphaMissense saw as "benign"
    in training, so LRs are optimistic, more so on the benign side.
  * Research use only. This is not a clinical calibration.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .table import PROCESSED

STRENGTHS = [("supporting", 2.08, 1), ("moderate", 4.33, 2), ("strong", 18.7, 4),
             ("very_strong", 350.0, 8)]
# name -> (column, sign). sign * column: higher = more damaging.
PREDICTORS = {"AlphaMissense": ("am_pathogenicity", 1.0), "ESM1b": ("esm1b_llr", -1.0)}
N_BOOT = 2000
GRID = 200
WINDOW_SHARE = 0.10  # each local estimate uses ~10% of the labelled variants
MIN_PER_CLASS = 30


def _windows(pos: np.ndarray, neg: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Grid of thresholds and a score window around each (signed scale).

    Each window spans the pooled variants ranked within +/- ``half`` of the grid
    point, so every local estimate rests on roughly the same amount of data.
    """
    pooled = np.sort(np.r_[pos, neg])
    n = len(pooled)
    half = max(15, int(WINDOW_SHARE * n / 2))
    idx = np.unique(np.linspace(0, n - 1, GRID).round().astype(int))
    lo = pooled[np.clip(idx - half, 0, n - 1)]
    hi = pooled[np.clip(idx + half, 0, n - 1)]
    return pooled[idx], lo, hi


def _local_lr(pos: np.ndarray, neg: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """LR = P(score in window | oncogenic) / P(score in window | benign)."""
    ps, ns = np.sort(pos), np.sort(neg)
    eps = 0.5  # continuity correction so empty windows don't give infinite LRs
    p = np.searchsorted(ps, hi, "right") - np.searchsorted(ps, lo, "left") + eps
    q = np.searchsorted(ns, hi, "right") - np.searchsorted(ns, lo, "left") + eps
    return (p / (len(ps) + eps)) / (q / (len(ns) + eps))


def _boot_bounds(pos, neg, lo, hi, seed=11):
    """Variant-level bootstrap: every variant treated as independent."""
    rng = np.random.default_rng(seed)
    b = np.array([_local_lr(rng.choice(pos, len(pos)), rng.choice(neg, len(neg)), lo, hi)
                  for _ in range(N_BOOT)])
    # One-sided 95%: lower bound for oncogenic evidence, upper bound for benign.
    return np.percentile(b, 5, axis=0), np.percentile(b, 95, axis=0)


def _window_counts(x: np.ndarray, genes: np.ndarray, names: np.ndarray, lo, hi):
    """Per gene: how many of its variants fall in each window, and how many it has."""
    counts = np.zeros((len(names), len(lo)))
    totals = np.zeros(len(names))
    for i, g in enumerate(names):
        s = np.sort(x[genes == g])
        counts[i] = np.searchsorted(s, hi, "right") - np.searchsorted(s, lo, "left")
        totals[i] = len(s)
    return counts, totals


def _gene_boot_bounds(pos, neg, pos_genes, neg_genes, lo, hi, seed=11):
    """Gene-level (cluster) bootstrap: resample whole genes with replacement."""
    names = np.unique(np.r_[pos_genes, neg_genes])
    pc, pn = _window_counts(pos, pos_genes, names, lo, hi)
    qc, qn = _window_counts(neg, neg_genes, names, lo, hi)
    rng = np.random.default_rng(seed)
    g = len(names)
    m = rng.multinomial(g, np.full(g, 1 / g), size=N_BOOT).astype(float)  # gene multiplicities
    eps = 0.5
    b = ((m @ pc + eps) / (m @ pn + eps)[:, None]) / ((m @ qc + eps) / (m @ qn + eps)[:, None])
    return np.percentile(b, 5, axis=0), np.percentile(b, 95, axis=0)


def _loosest(lb, ub, direction: str, cut: float):
    """Index of the loosest grid point that earns ``cut``, or None.

    A threshold earns a strength only if every score beyond it does too: running minimum
    of the lower bound from the top, running maximum of the upper bound from the bottom.
    """
    if direction == "oncogenic":
        ok = np.where(np.minimum.accumulate(lb[::-1])[::-1] >= cut)[0]
        return ok[0] if len(ok) else None
    ok = np.where(np.maximum.accumulate(ub) <= 1 / cut)[0]
    return ok[-1] if len(ok) else None


def _fmt(t: float, direction: str, sign: float) -> str:
    op = ">=" if (direction == "oncogenic") == (sign > 0) else "<="
    return f"{op}{sign * t:.3f}"


def trust_rows(pos: np.ndarray, neg: np.ndarray, sign: float, labels: str, role: str,
               predictor: str, pos_genes: np.ndarray | None = None,
               neg_genes: np.ndarray | None = None) -> list[dict]:
    """Thresholds per strength. With gene names, bounds come from resampling whole genes."""
    grid, lo, hi = _windows(pos, neg)
    lr = _local_lr(pos, neg, lo, hi)
    v_lb, v_ub = _boot_bounds(pos, neg, lo, hi)
    by_gene = pos_genes is not None and neg_genes is not None
    lb, ub = _gene_boot_bounds(pos, neg, pos_genes, neg_genes, lo, hi) if by_gene else (v_lb, v_ub)
    base = {"labels": labels, "gene_role": role, "predictor": predictor,
            "n_oncogenic": len(pos), "n_benign": len(neg),
            "resampling": "genes" if by_gene else "variants",
            "n_genes": len(np.unique(np.r_[pos_genes, neg_genes])) if by_gene else np.nan}
    rows = []
    for direction in ("oncogenic", "benign"):
        for name, cut, pts in STRENGTHS:
            j = _loosest(v_lb, v_ub, direction, cut)
            row = {**base, "direction": direction, "strength": name,
                   "points": pts if direction == "oncogenic" else -pts,
                   "threshold_variant_level":
                       "not reached" if j is None else _fmt(grid[j], direction, sign)}
            i = _loosest(lb, ub, direction, cut)
            if i is None:
                rows.append({**row, "threshold": "not reached"})
                continue
            t = grid[i]
            if direction == "oncogenic":
                reach_p, reach_n, bound = (pos >= t).mean(), (neg >= t).mean(), lb[i]
            else:
                reach_p, reach_n, bound = (pos <= t).mean(), (neg <= t).mean(), ub[i]
            rows.append({**row, "threshold": _fmt(t, direction, sign), "local_lr": lr[i],
                         "lr_bound": bound, "share_oncogenic_reached": float(reach_p),
                         "share_benign_reached": float(reach_n)})
    return rows


def _expert_population() -> pd.DataFrame:
    d = pd.read_parquet(PROCESSED / "civic_missense_alphamissense.parquet")
    cols = [c for c, _ in PREDICTORS.values()]
    d = d[d[cols].notna().all(axis=1) & d["label"].notna()].copy()
    d["label"] = d["label"].astype(int)
    return d


def _table(d: pd.DataFrame, labels: str, roles: list[str], label_col: str = "label") -> list[dict]:
    rows = []
    for role in roles:
        sub = d if role == "all" else d[d["gene_role"] == role]
        for name, (col, sign) in PREDICTORS.items():
            is_pos, is_neg = sub[label_col] == 1, sub[label_col] == 0
            pos = sign * sub.loc[is_pos, col].to_numpy()
            neg = sign * sub.loc[is_neg, col].to_numpy()
            if len(pos) < MIN_PER_CLASS or len(neg) < MIN_PER_CLASS:
                continue
            rows += trust_rows(pos, neg, sign, labels, role, name,
                               sub.loc[is_pos, "gene"].to_numpy(),
                               sub.loc[is_neg, "gene"].to_numpy())
    return rows


def validate(d: pd.DataFrame, trained: pd.DataFrame) -> pd.DataFrame:
    """Check thresholds learnt without the time-split positives on those positives.

    For each strength, the held-out LR is measured inside its own score band (from its
    threshold up to the next strength's threshold), so it is a local check too.
    """
    test = pd.read_parquet(PROCESSED / "time_split_test.parquet")
    keys = set(zip(test["gene"], test["variant"]))
    is_test = np.array([k in keys for k in zip(d["gene"], d["variant"])])
    cuts = {s: c for s, c, _ in STRENGTHS}
    t = trained[(trained["direction"] == "oncogenic") & (trained["threshold"] != "not reached")]
    out = []
    for (role, pred), g in t.groupby(["gene_role", "predictor"], sort=False):
        col, sign = PREDICTORS[pred]
        thr = [sign * float(x[2:]) for x in g["threshold"]] + [np.inf]
        pos = sign * d.loc[is_test & (d["label"] == 1) & (d["gene_role"] == role), col].to_numpy()
        neg = sign * d.loc[(d["label"] == 0) & (d["gene_role"] == role), col].to_numpy()
        for (_, r), lo, hi in zip(g.iterrows(), thr, thr[1:]):
            p = ((pos >= lo) & (pos < hi)).sum()
            n = ((neg >= lo) & (neg < hi)).sum()
            lr = ((p + 0.5) / (len(pos) + 0.5)) / ((n + 0.5) / (len(neg) + 0.5))
            out.append({"gene_role": role, "predictor": pred, "strength": r["strength"],
                        "threshold": r["threshold"], "n_heldout": len(pos),
                        "heldout_in_band": int(p), "benign_in_band": int(n),
                        "band_lr_heldout": lr, "holds": bool(lr >= cuts[r["strength"]])})
    return pd.DataFrame(out)


def build() -> pd.DataFrame:
    d = _expert_population()
    roles = ["oncogene", "tumour_suppressor", "all"]
    rows = _table(d, "expert+population", roles)

    # Lab-measured function: brakes (damaging vs normal) and gas pedals (activating vs
    # neutral-or-loss, since a switched-off oncogene is not oncogenic).
    lab_path = PROCESSED / "lab_tested.parquet"
    if lab_path.exists():
        lab = pd.read_parquet(lab_path)
        lab = lab[lab["mave_label"].notna()].copy()
        lab["label"] = lab["mave_label"].astype(int)
        if "gene_role" not in lab:
            lab["gene_role"] = "tumour_suppressor"
        rows += _table(lab[lab["gene_role"] == "tumour_suppressor"],
                       "lab function (TP53, PTEN, BRCA1)", ["tumour_suppressor"])
        rows += _table(lab[lab["gene_role"] == "oncogene"],
                       "lab activity (KRAS, BRAF, EGFR, ERBB2, MAP2K1)", ["oncogene"])
    res = pd.DataFrame(rows)
    res.to_csv(PROCESSED / "trust_table.csv", index=False)

    # Held-out check: re-learn thresholds without the post-release positives.
    test = pd.read_parquet(PROCESSED / "time_split_test.parquet")
    keys = set(zip(test["gene"], test["variant"]))
    train = d[[k not in keys for k in zip(d["gene"], d["variant"])]]
    trained = pd.DataFrame(_table(train, "train (pre-release positives)", roles[:2]))
    val = validate(d, trained)
    val.to_csv(PROCESSED / "trust_table_validation.csv", index=False)

    show = res[res["threshold"] != "not reached"][
        ["labels", "gene_role", "n_genes", "predictor", "direction", "strength", "points",
         "threshold", "threshold_variant_level", "local_lr", "lr_bound", "share_oncogenic_reached", "share_benign_reached"]]
    with pd.option_context("display.float_format", "{:.2f}".format, "display.width", 220,
                           "display.max_rows", 200):
        print(show.to_string(index=False))
        print("\nNot reached:")
        nr = res[res["threshold"] == "not reached"]
        print(nr.groupby(["labels", "gene_role", "predictor", "direction"])["strength"]
                .agg(", ".join).to_string())
        print("\nHeld-out check (thresholds learnt without post-release positives):")
        print(val.to_string(index=False))
    print("\nBounds resample whole genes; threshold_variant_level is what resampling single"
          " variants would have given.")
    print("\nlocal_lr = LR among variants scoring near the threshold; lr_bound = its one-sided"
          " 95% bound, which must clear the cut-off at the threshold and every score beyond.")
    print("Negatives are mostly common gnomAD variants: LRs are optimistic. Research use only.")
    return res
