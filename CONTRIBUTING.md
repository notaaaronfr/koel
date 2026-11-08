# Contributing to Koel

Thank you for helping. Koel is a research tool, so correctness and openness come first.

## Ways to help

- **Data:** point us to openly licensed (CC0 / CC BY / public domain) lab datasets,
  especially deep mutational scans of **oncogenes**, or expert-confirmed benign somatic
  variants. Open an issue with the "Data source" template.
- **Methods:** calibration, statistics, new predictors (only if their licence allows).
- **Domain review:** variant curators and pathologists, tell us where the trust table
  does not match practice.
- **Code:** bug fixes, tests, speed-ups.

## Ground rules

1. **Open data only.** No controlled-access, patient-identifiable or non-commercial-only
   data may be committed. Every new source goes in `src/koel/sources.py` with its licence.
2. **Never commit data files.** `data/` is git-ignored; downloads are reproducible.
3. **Report uncertainty.** New results come with confidence intervals and caveats.
4. **Research use only.** Nothing may suggest Koel interprets an individual patient.

## Development

```bash
uv sync --group dev
uv run pytest
```

Tests run on small synthetic data and need no downloads. Please add a test with
any change to `calibrate.py` or `analysis.py`.

## Pull requests

Keep them focused, explain *why*, and note any change to published numbers.
