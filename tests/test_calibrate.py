"""Calibration tests on synthetic data (no downloads needed)."""

import numpy as np

from koel.calibrate import (STRENGTHS, _boot_bounds, _gene_boot_bounds, _local_lr, _windows,
                            trust_rows)


def _synthetic(n=400, seed=0):
    rng = np.random.default_rng(seed)
    pos = np.clip(rng.normal(0.85, 0.12, n), 0, 1)
    neg = np.clip(rng.normal(0.25, 0.15, n), 0, 1)
    return pos, neg


def test_local_lr_is_high_for_high_scores_and_low_for_low_scores():
    pos, neg = _synthetic()
    grid, lo, hi = _windows(pos, neg)
    lr = _local_lr(pos, neg, lo, hi)
    assert lr[grid > 0.9].min() > 2.08
    assert lr[grid < 0.2].max() < 1 / 2.08


def test_thresholds_are_ordered_by_strength():
    pos, neg = _synthetic()
    rows = [r for r in trust_rows(pos, neg, 1.0, "synthetic", "x", "P")
            if r["direction"] == "oncogenic" and r["threshold"] != "not reached"]
    thr = [float(r["threshold"][2:]) for r in rows]
    assert len(thr) >= 2
    assert thr == sorted(thr)


def test_no_signal_earns_no_evidence():
    rng = np.random.default_rng(1)
    pos, neg = rng.uniform(size=300), rng.uniform(size=300)
    rows = trust_rows(pos, neg, 1.0, "noise", "x", "P")
    assert all(r["threshold"] == "not reached" for r in rows)


def test_negative_sign_reports_original_units():
    pos, neg = _synthetic()
    rows = trust_rows(pos, neg, -1.0, "s", "x", "P")
    onc = [r for r in rows if r["direction"] == "oncogenic" and r["threshold"] != "not reached"]
    assert onc and onc[0]["threshold"].startswith("<=-")


def test_strength_cutoffs_match_tavtigian():
    assert [c for _, c, _ in STRENGTHS] == [2.08, 4.33, 18.7, 350.0]


def test_gene_bootstrap_matches_variant_bootstrap_when_every_variant_is_its_own_gene():
    pos, neg = _synthetic(200)
    _, lo, hi = _windows(pos, neg)
    v_lb, v_ub = _boot_bounds(pos, neg, lo, hi)
    genes_p, genes_n = np.arange(200), np.arange(200, 400)
    g_lb, g_ub = _gene_boot_bounds(pos, neg, genes_p, genes_n, lo, hi)
    for gene, variant in ((g_lb, v_lb), (g_ub, v_ub)):
        gap = np.abs(np.log(gene) - np.log(variant))
        assert np.median(gap) < 0.1  # same answer, up to bootstrap noise
        assert np.percentile(gap, 95) < 0.3


def test_gene_bootstrap_is_wider_when_one_gene_carries_the_signal():
    rng = np.random.default_rng(3)
    # Ten genes; only gene 0 has high-scoring positives. The rest are noise.
    pos = np.r_[rng.uniform(0.9, 1.0, 200), rng.uniform(0, 1, 180)]
    pos_genes = np.r_[np.zeros(200, int), np.repeat(np.arange(1, 10), 20)]
    neg = rng.uniform(0, 1, 400)
    neg_genes = np.repeat(np.arange(10), 40)
    by_variant = trust_rows(pos, neg, 1.0, "s", "x", "P")
    by_gene = trust_rows(pos, neg, 1.0, "s", "x", "P", pos_genes, neg_genes)
    reached = lambda rows: [r["strength"] for r in rows
                            if r["direction"] == "oncogenic" and r["threshold"] != "not reached"]
    assert "supporting" in reached(by_variant)
    assert reached(by_gene) == []
    assert by_gene[0]["threshold_variant_level"] == by_variant[0]["threshold"]
