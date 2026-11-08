"""Download Koel's open data sources with resume support and a provenance log."""

from __future__ import annotations

import hashlib
import json
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from .sources import Source

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
MANIFEST = RAW / "manifest.json"

USER_AGENT = "koel-research-downloader/0.1"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_manifest() -> dict:
    if MANIFEST.exists():
        return json.loads(MANIFEST.read_text())
    return {}


def _save_manifest(m: dict) -> None:
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(m, indent=2, sort_keys=True))


def fetch(source: Source, force: bool = False) -> Path:
    dest = RAW / source.filename
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        print(f"[skip] {source.key}: already at {dest.relative_to(ROOT)}")
        return dest

    part = dest.with_suffix(dest.suffix + ".part")
    start = part.stat().st_size if part.exists() else 0
    req = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    if start:
        req.add_header("Range", f"bytes={start}-")

    print(f"[get ] {source.key} ({source.approx_size}) <- {source.url}")
    with urllib.request.urlopen(req, timeout=60) as resp:
        mode = "ab" if start and resp.status == 206 else "wb"
        with part.open(mode) as out:
            shutil.copyfileobj(resp, out, length=1 << 20)
    part.rename(dest)

    m = _load_manifest()
    m[source.key] = {
        "url": source.url,
        "file": str(dest.relative_to(ROOT)),
        "licence": source.licence,
        "bytes": dest.stat().st_size,
        "sha256": _sha256(dest),
        "downloaded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    _save_manifest(m)
    print(f"[done] {source.key}: {dest.stat().st_size / 1e6:.1f} MB")
    return dest
