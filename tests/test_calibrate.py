"""Calibration tests on synthetic data (no downloads needed)."""

import numpy as np

from koel.calibrate import STRENGTHS, _local_lr, _windows, trust_rows


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
