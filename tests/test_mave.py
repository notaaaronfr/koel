"""Oncogene lab-label logic on a tiny synthetic score set (no downloads)."""

import pandas as pd

from koel import mave


def test_onco_labels_use_author_calls(monkeypatch):
    fake = pd.DataFrame({
        "hgvs_pro": ["KRAS:p.Gly12Asp", "KRAS:p.Gly12Asp", "KRAS:p.Ala11Val",
                     "KRAS:p.Gly10Arg", "KRAS:p.Leu6Pro", "KRAS:p.Ala11=", "KRAS:p.Ala11Ter"],
        "score": [4.0, 4.2, 1.0, 0.3, 2.0, 1.0, 0.2],
        "variant_classification": ["high", "high", "wt-like", "low", "high", "wt-like", "low"],
    })
    # Leu6Pro measured twice with conflicting calls -> dropped.
    fake = pd.concat([fake, pd.DataFrame({"hgvs_pro": ["KRAS:p.Leu6Pro"], "score": [1.0],
                                          "variant_classification": ["wt-like"]})])
    monkeypatch.setattr(mave, "_scores", lambda urn: fake.reset_index(drop=True))
    out, stats = mave.onco_assay_labels({"gene": "KRAS", "urn": "x", "label": "test"})
    got = dict(zip(out["variant"], out["mave_effect"]))
    assert got == {"G12D": "activating", "A11V": "neutral", "G10R": "loss"}
    assert dict(zip(out["variant"], out["mave_label"])) == {"G12D": 1, "A11V": 0, "G10R": 0}
    assert stats["n_damaging"] == 1 and stats["n_loss"] == 1
