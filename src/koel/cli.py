"""Koel command line: ``uv run koel <command>``."""

from __future__ import annotations

import argparse
import sys

from . import download
from .sources import PENDING, SOURCES, by_group


def _cmd_download(args: argparse.Namespace) -> int:
    if args.list or not args.targets:
        print(f"{'key':28} {'size':10} {'licence':10} redistributable")
        for s in SOURCES:
            print(f"{s.key:28} {s.approx_size:10} {s.licence:10} {s.redistributable}")
        print("\nStill to choose:")
        for k, v in PENDING.items():
            print(f"  {k}: {v}")
        return 0
    for target in args.targets:
        chosen = SOURCES if target == "all" else by_group(target)
        if not chosen:
            print(f"Unknown source or group: {target}", file=sys.stderr)
            return 2
        for s in chosen:
            download.fetch(s, force=args.force)
    return 0


def _cmd_all(args: argparse.Namespace) -> int:
    from . import analysis, calibrate, table
    for title, step in [("Building the labelled table", table.build),
                        ("Accuracy per gene type", analysis.summarise),
                        ("Lab-measured function", analysis.lab_tested),
                        ("Held-out time split", analysis.time_split),
                        ("Trust table", calibrate.build)]:
        print(f"\n=== {title} ===")
        step()
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="koel", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    d = sub.add_parser("download", help="download open data sources")
    d.add_argument("targets", nargs="*", help="source key, group (civic, alphamissense) or 'all'")
    d.add_argument("--list", action="store_true", help="list sources and licences")
    d.add_argument("--force", action="store_true", help="re-download even if present")
    d.set_defaults(func=_cmd_download)

    t = sub.add_parser("table", help="build the labelled CIViC x AlphaMissense table")
    t.set_defaults(func=lambda a: (__import__("koel.table", fromlist=["build"]).build(), 0)[1])

    a = sub.add_parser("analyse", help="first per-gene-class AlphaMissense AUROC")
    a.set_defaults(func=lambda a: (__import__("koel.analysis", fromlist=["summarise"]).summarise(), 0)[1])

    lt = sub.add_parser("lab", help="compare predictors with lab-measured function (TP53, PTEN, BRCA1)")
    lt.set_defaults(func=lambda a: (__import__("koel.analysis", fromlist=["lab_tested"]).lab_tested(), 0)[1])

    ts = sub.add_parser("timesplit", help="held-out test on oncogenic calls made after AlphaMissense's release")
    ts.set_defaults(func=lambda a: (__import__("koel.analysis", fromlist=["time_split"]).time_split(), 0)[1])

    tt = sub.add_parser("trust", help="trust table: evidence strength per predictor and gene type")
    tt.set_defaults(func=lambda a: (__import__("koel.calibrate", fromlist=["build"]).build(), 0)[1])

    fg = sub.add_parser("figures", help="redraw the README chart (needs: uv sync --group figures)")
    fg.set_defaults(func=lambda a: (__import__("koel.figures", fromlist=["build"]).build(), 0)[1])

    al = sub.add_parser("all", help="build the table and run every analysis in order")
    al.set_defaults(func=_cmd_all)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
