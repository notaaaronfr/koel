"""Export Koel's headline results as one JSON file for the website.

``uv run koel export <path>`` writes a small, self-describing JSON so the site never
hand-copies numbers. Run it after ``koel all``.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pandas as pd

from .table import PROCESSED


def _num(x, nd=3):
    return None if pd.isna(x) else round(float(x), nd)


def build(path: str | Path) -> Path:
    trust = pd.read_csv(PROCESSED / "trust_table.csv")
    lab = pd.read_csv(PROCESSED / "lab_tested_auroc.csv")
    labstats = pd.read_csv(PROCESSED / "lab_tested_assays.csv")
    ts = pd.read_csv(PROCESSED / "time_split_results.csv")
    val = pd.read_csv(PROCESSED / "trust_table_validation.csv")

    t = trust[trust["threshold"] != "not reached"].copy()
    trust_rows = [{
        "labels": r.labels, "geneRole": r.gene_role, "predictor": r.predictor,
        "direction": r.direction, "strength": r.strength, "points": int(r.points),
        "threshold": r.threshold, "lrBound": _num(r.lr_bound, 2),
        "shareOncogenicReached": _num(r.share_oncogenic_reached),
        "shareBenignReached": _num(r.share_benign_reached),
    } for r in t.itertuples()]
    not_reached = [{
        "labels": r.labels, "geneRole": r.gene_role, "predictor": r.predictor,
        "direction": r.direction, "strength": r.strength,
    } for r in trust[trust["threshold"] == "not reached"].itertuples()]

    lab_rows = [{
        "gene": r.gene, "geneRole": r.gene_role, "comparison": r.comparison,
        "predictor": r.predictor, "nPos": int(r.n_pos), "nNeg": int(r.n_neg),
        "auroc": _num(r.auroc), "ciLow": _num(r.ci_low), "ciHigh": _num(r.ci_high),
    } for r in lab.itertuples()]
    assays = [{
        "gene": r.gene, "assay": r.assay, "urn": r.urn, "geneRole": r.gene_role,
        "nMissense": int(r.n_missense),
    } for r in labstats.itertuples()]

    ts_rows = [{
        "geneRole": r.gene_role, "predictor": r.predictor, "threshold": r.threshold,
        "nTest": int(r.n_test), "rate": _num(r.detection_rate),
        "ciLow": _num(r.det_ci_low), "ciHigh": _num(r.det_ci_high),
        "falseFlagRate": _num(r.false_flag_rate),
    } for r in ts.itertuples()]

    out = {
        "generated": date.today().isoformat(),
        "repo": "https://github.com/notaaaronfr/koel",
        "note": "Research use only. Derived from AlphaMissense (CC BY-NC-SA 4.0) and "
                "ESM1b (CC BY-NC 4.0) scores; non-commercial use.",
        "trust": trust_rows, "trustNotReached": not_reached,
        "trustValidation": val.where(val.notna(), None).to_dict("records"),
        "lab": lab_rows, "labAssays": assays, "timeSplit": ts_rows,
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    print(f"Wrote {path}")
    return path
