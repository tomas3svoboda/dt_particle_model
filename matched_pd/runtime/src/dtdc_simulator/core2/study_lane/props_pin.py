"""Era-binding for the study-lane surrogate pack.

A surrogate is a *fit to a specific implementation*.  If the pinned property
authority is re-keyed — a coefficient corrected, a tolerance changed, a solver
branch restructured (D1 is exactly such an event) — every coefficient in the
pack becomes a fit to a function that no longer exists.  The pack is therefore
bound to the byte content of the pinned property surfaces it was built against,
and the binding is asserted at LOAD time, not at build time: a stale pack
refuses to activate rather than silently answering with last era's physics.

The digest is taken over the source bytes exactly as stored (``rb``), so it is
insensitive to nothing at all — comments included.  That is deliberate.  A
comment-only edit to a byte-pinned surface is itself a pin event, and re-running
``scripts/build_study_surrogates.py`` is cheap.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

__all__ = [
    "PINNED_PROPS_FILES",
    "props_pin_digests",
    "props_source_root",
]

#: Every pinned property surface the surrogate pack is fitted against, relative
#: to the ``dtdc_simulator/core2`` package directory.
PINNED_PROPS_FILES: tuple[str, ...] = (
    "props/__init__.py",
    "props/binary_gas.py",
    "props/critical_volume.py",
    "props/hexane.py",
    "props/sorption.py",
    "props/state.py",
    "props/water.py",
)


def props_source_root() -> Path:
    """Directory holding ``props/`` (the installed ``core2`` package)."""

    return Path(__file__).resolve().parent.parent


def props_pin_digests() -> dict[str, str]:
    """``{relative path: 'sha256:...'}`` over the pinned property sources."""

    root = props_source_root()
    digests: dict[str, str] = {}
    for relative in PINNED_PROPS_FILES:
        payload = (root / relative).read_bytes()
        digests[relative] = "sha256:" + hashlib.sha256(payload).hexdigest()
    return digests
