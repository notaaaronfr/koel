# Koel

**An open trust checker for cancer mutation scores.**

Cancer tumour tests find many mutations, and most are labelled *uncertain*: nobody
knows whether they drive the cancer. Computer predictors such as AlphaMissense give
each mutation a "how harmful" score, but they are less reliable for mutations that
switch a gene **on** (oncogenes, the gas pedal) than for mutations that switch a gene
**off** (tumour suppressors, the brakes).

Koel measures that gap openly and turns it into the evidence strengths that variant
curators actually use (ClinGen/CGC/VICC somatic oncogenicity rules, criterion OP1),
separately for each gene class and each predictor.

> **Research use only.** Koel does not interpret any individual patient's results and
> must never be used to make diagnosis or treatment decisions.

## What Koel will provide

1. **Engine:** a reproducible benchmark that re-runs end to end with one command.
2. **Trust tables:** per-predictor, per-gene-class calibrated thresholds and
   evidence strengths, downloadable as files other tools (PCGR, VEP, cBioPortal) can load.
3. **Platform:** a website with gene cards, mutation pages, a predictor leaderboard
   and an API. Proven evidence and hypotheses are always shown separately.

## Data policy

Only openly licensed data is used (for example CIViC, CC0). AlphaMissense
predictions are CC BY-NC-SA 4.0 (non-commercial, share-alike), so any Koel output
that contains them is non-commercial and carries the same licence. Non-commercial resources (OncoKB, COSMIC,
CADD, REVEL) are never bundled. No controlled-access or patient-identifiable data is
ever stored here. See `src/koel/sources.py` for every source and its licence.

## Prior work

Koel builds on, and does not claim the discovery of, the predictor gap reported in
Lee, *An openly licensed benchmark and per-gene calibration map for missense
pathogenicity predictors on activating cancer drivers*, bioRxiv (23 Jul 2026).

## Quick start

```bash
uv sync
uv run koel download --list      # show data sources, sizes and licences
uv run koel download civic       # small files
uv run koel download alphamissense   # about 1.2 GB
uv run koel table                # build the labelled table
uv run koel analyse              # AUROC per gene type
uv run koel lab                  # compare with lab-measured function
uv run koel timesplit            # held-out test on calls made after AlphaMissense
uv run koel trust                # trust table (evidence strength per gene type)
```

## First trust table (v0, research use only)

`uv run koel trust` writes `data/processed/trust_table.csv`. For each predictor and
gene type it gives the score needed for each ClinGen evidence strength, using local
likelihood ratios with a bootstrap 95% bound (method in `src/koel/calibrate.py`).

AlphaMissense, oncogenic side, expert labels plus common gnomAD negatives:

| Strength (points) | Oncogenes (gas pedal) | Tumour suppressors (brakes) |
|---|---|---|
| Supporting (+1) | ≥ 0.885 (77% of oncogenic variants reach it) | ≥ 0.830 (83%) |
| Moderate (+2) | ≥ 0.982 (54%) | ≥ 0.943 (73%) |

In oncogenes, the same strength of evidence needs a higher score, and fewer real
oncogenic variants reach it. These thresholds, re-learnt without variants classified
after AlphaMissense's release, held on those held-out variants at the supporting and
moderate levels. Caveats: the negatives are mostly common population variants (which
make predictors look better), some curators may have seen AlphaMissense scores, and
labels are a few hundred variants per gene type. Not for clinical use.
