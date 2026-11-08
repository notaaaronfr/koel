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

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
