"""README figure: how much evidence an AlphaMissense score is worth, per gene type.

Writes docs/figures/trust_curve_{light,dark}.png. Needs the optional ``figures``
dependency group: ``uv sync --group figures``.
"""

from __future__ import annotations

import numpy as np

from .calibrate import STRENGTHS, _expert_population, _local_lr, _windows
from .download import ROOT

OUT = ROOT / "docs" / "figures"
ROLES = [("oncogene", "Oncogenes (gas pedal)"), ("tumour_suppressor", "Tumour suppressors (brakes)")]
THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#e4e3df",
              "series": ["#2a78d6", "#eb6834"]},
    "dark": {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "grid": "#3a3936",
             "series": ["#3987e5", "#d95926"]},
}


def _curves():
    d = _expert_population()
    out = {}
    for role, _ in ROLES:
        sub = d[d["gene_role"] == role]
        pos = sub.loc[sub["label"] == 1, "am_pathogenicity"].to_numpy()
        neg = sub.loc[sub["label"] == 0, "am_pathogenicity"].to_numpy()
        grid, lo, hi = _windows(pos, neg)
        lr = _local_lr(pos, neg, lo, hi)
        # Light smoothing for display only (the table uses the unsmoothed values).
        k = 9
        sm = np.convolve(np.log(lr), np.ones(k) / k, mode="same")
        sm[: k // 2], sm[-(k // 2):] = np.log(lr[: k // 2]), np.log(lr[-(k // 2):])
        reach = _reach(role)
        out[role] = (grid, np.exp(sm), reach)
    return out


def _reach(role: str) -> dict[str, float]:
    """Share of oncogenic variants reaching each strength, exactly as in the trust table."""
    import pandas as pd

    from .table import PROCESSED
    t = pd.read_csv(PROCESSED / "trust_table.csv")
    t = t[(t["labels"] == "expert+population") & (t["gene_role"] == role)
          & (t["predictor"] == "AlphaMissense") & (t["direction"] == "oncogenic")
          & (t["threshold"] != "not reached")]
    return {r.strength: float(r.share_oncogenic_reached) for r in t.itertuples()
            if r.strength in ("supporting", "moderate")}


def _thresholds(role: str) -> dict[str, float]:
    import pandas as pd

    from .table import PROCESSED
    t = pd.read_csv(PROCESSED / "trust_table.csv")
    t = t[(t["labels"] == "expert+population") & (t["gene_role"] == role)
          & (t["predictor"] == "AlphaMissense") & (t["direction"] == "oncogenic")
          & (t["threshold"] != "not reached")]
    return {r.strength: float(r.threshold[2:]) for r in t.itertuples()
            if r.strength in ("supporting", "moderate")}


def build() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    curves = _curves()
    OUT.mkdir(parents=True, exist_ok=True)
    for mode, th in THEMES.items():
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                             "text.color": th["ink"], "axes.labelcolor": th["ink2"],
                             "xtick.color": th["ink2"], "ytick.color": th["ink2"]})
        fig, (ax, bx) = plt.subplots(1, 2, figsize=(12, 4.6), dpi=160,
                                     gridspec_kw={"width_ratios": [1.7, 1]})
        fig.patch.set_facecolor(th["surface"])
        for a in (ax, bx):
            a.set_facecolor(th["surface"])
            for s in ("top", "right"):
                a.spines[s].set_visible(False)
            for s in ("left", "bottom"):
                a.spines[s].set_color(th["grid"])

        # Left: local likelihood ratio along the score.
        for (cut_name, cut, _) in STRENGTHS[:3]:
            ax.axhline(cut, color=th["grid"], lw=1, ls=(0, (4, 3)), zorder=1)
            ax.text(0.105, cut * 1.08, cut_name, color=th["ink2"], fontsize=9, va="bottom")
        ax.axhline(1, color=th["grid"], lw=1, zorder=1)
        for (role, label), color in zip(ROLES, th["series"]):
            g, lr, _ = curves[role]
            keep = (g >= 0.1)
            keep &= g <= 0.995  # last windows run off the top of the scale
            ax.plot(g[keep], lr[keep], color=color, lw=2, solid_capstyle="round", zorder=3)
        ax.text(0.86, 11, "Tumour suppressors", color=th["ink"], fontsize=10, ha="right")
        ax.text(0.97, 2.4, "Oncogenes", color=th["ink"], fontsize=10, ha="right")
        ax.set_yscale("log")
        ax.set_xlim(0.1, 1)
        ax.set_ylim(0.03, 100)
        ax.set_xlabel("AlphaMissense score")
        ax.set_ylabel("Evidence for 'oncogenic' (local likelihood ratio)")
        ax.set_title("Expert labels: point estimates", loc="left", fontsize=13,
                     fontweight="bold", color=th["ink"])
        ax.grid(False)

        # Right: share of real oncogenic variants that reach each strength.
        levels = ["supporting", "moderate"]
        y = np.arange(len(levels))
        h = 0.36
        for i, ((role, label), color) in enumerate(zip(ROLES, th["series"])):
            vals = [curves[role][2].get(s, 0) for s in levels]
            bars = bx.barh(y + (i - 0.5) * (h + 0.04), vals, height=h, color=color, zorder=3)
            for b, v in zip(bars, vals):
                bx.text(v + 0.015, b.get_y() + h / 2, f"{v:.0%}", va="center",
                        color=th["ink"], fontsize=10)
        bx.set_yticks(y, [f"{s.capitalize()}\n(+{p} pt)" for s, p in zip(levels, (1, 2))])
        bx.invert_yaxis()
        bx.set_xlim(0, 1.12)
        bx.set_xticks([0, 0.5, 1], ["0%", "50%", "100%"])
        bx.set_title("Mutations earning the points", loc="left", fontsize=13,
                     fontweight="bold", color=th["ink"])
        handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in th["series"]]
        bx.legend(handles, [lab for _, lab in ROLES], frameon=False, loc="upper left",
                  bbox_to_anchor=(-0.02, -0.12), ncol=1, fontsize=9, labelcolor=th["ink2"])

        fig.text(0.01, 0.01, "Koel v0 · expert labels (CIViC, ClinVar) + common gnomAD negatives · "
                 "bars use bounds that resample whole genes · research use only", color=th["ink2"], fontsize=8)
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        path = OUT / f"trust_curve_{mode}.png"
        fig.savefig(path, facecolor=th["surface"])
        plt.close(fig)
        print(f"Wrote {path.relative_to(ROOT)}")
