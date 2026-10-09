"""Era stamp for analysis runs.

Every harness output must embed this stamp so a number can always be traced
to the exact code state and inputs that produced it. Read-only: runs git
purely to *read* the repository state.

Usage:
    from era_stamp import era_stamp
    stamp = era_stamp(input_paths=["../datasets/faner2019/soybean_curves.csv"])
    json.dump({"stamp": stamp, "results": ...}, out)
"""
from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def _git(*args: str) -> str:
    out = subprocess.run(
        ["git", "-C", str(REPO_ROOT), *args],
        capture_output=True, text=True, check=True,
    )
    return out.stdout.strip()


def _is_git_checkout() -> bool:
    out = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "--is-inside-work-tree"],
                         capture_output=True, text=True, check=False)
    return out.returncode == 0 and out.stdout.strip() == "true"


def _package_source() -> dict:
    """The public reproducibility package is not a git checkout; it carries the
    commit it was cut from in BUILD_RECORD.json at its root."""
    record = REPO_ROOT / "BUILD_RECORD.json"
    if record.exists():
        import json
        built = json.loads(record.read_text(encoding="utf-8"))
        return {"repo_head": f"public package cut from {built.get('source_commit')}",
                "repo_branch": built.get("source_branch"), "worktree_dirty": None}
    raise RuntimeError("not a git checkout and no BUILD_RECORD.json: the era cannot be stamped")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def era_stamp(input_paths: list[str | Path] | None = None) -> dict:
    """Collect the era stamp. Raises if git state cannot be read."""
    if _is_git_checkout():
        source = {"repo_head": _git("rev-parse", "HEAD"),
                  "repo_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
                  "worktree_dirty": _git("status", "--porcelain") != ""}
    else:
        source = _package_source()
    try:
        import numpy
        numpy_version = numpy.__version__
    except ImportError:
        numpy_version = None
    inputs = {}
    for p in input_paths or []:
        p = Path(p)
        if not p.is_absolute():
            p = (Path(__file__).resolve().parent / p).resolve()
        inputs[str(p.relative_to(REPO_ROOT)) if p.is_relative_to(REPO_ROOT)
               else str(p)] = _sha256(p)
    return {
        **source,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.split()[0],
        "numpy": numpy_version,
        "platform": platform.platform(),
        "input_sha256": inputs,
        "claim_note": "era-bound analysis output; regenerate after any change "
                      "to the model source; not a qualification artifact",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(era_stamp(), indent=2))
