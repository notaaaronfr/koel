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
```

## Dedication

Named for **Kokila**, whom her family also calls **Jenny**. *Kokila* means the koel,
the songbird whose call announces a new season.
