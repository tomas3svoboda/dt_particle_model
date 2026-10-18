"""The study-lane surrogate pack: schema, digest, load-time era assertion.

The pack is a plain JSON document under ``docs/evidence/``.  It carries, for
every surrogated call:

* the DECLARED DOMAIN (an interval, or a box) outside which the surrogate is
  never evaluated — the facade falls back to the pinned authority instead;
* the Chebyshev degree and coefficients;
* the MEASURED verification certificate — sup-norm and maximum relative error
  over an independent grid at least ten times denser than the fit grid, with
  the argument at which the maximum occurred;
* an ``adopted`` flag.  A surrogate whose measured error misses the pack's
  declared acceptance threshold is serialized anyway (the honest negative is
  evidence) but is NOT installed.

Nothing in this module evaluates physics.  Loading a pack asserts the props
era-binding and refuses on mismatch; installing it is a separate, explicitly
opted-in act (see ``activation`` and ``surrogate_props``).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dtdc_simulator.core2.study_lane.props_pin import props_pin_digests

__all__ = [
    "DEFAULT_PACK_PATH",
    "PACK_RECORD_ID",
    "PackEraMismatch",
    "SurrogatePack",
    "canonical_pack_digest",
    "default_pack_path",
]

#: Re-keyed AGAIN 2026-09-03 by the D3 props fail-open ceremony (commit
#: 0992d5d, record ``docs/GT_PS2_D3_PROPS_FAIL_OPEN_CEREMONY_2026-09-03.md``),
#: which re-keyed ``props/water.py`` and ``props/hexane.py`` a second time on
#: the same calendar day -- hence the ``B`` suffix: the plain ``2026-09-03``
#: name belongs to the D2-era pack and both are history.  The era-binding below
#: FIRED on that pack (68 pure ``PackEraMismatch``), exactly as designed, and
#: the pack was rebuilt rather than re-pointed.  Superseded packs stay on disk.
#: previous value 2026-09-03 (D2 era): "GT-PS-2-STUDY-SURROGATE-PACK-2026-09-03-01"
#: previous value 2026-09-02: "GT-PS-2-STUDY-SURROGATE-PACK-2026-09-02-01"
PACK_RECORD_ID = "GT-PS-2-STUDY-SURROGATE-PACK-2026-09-03B-01"

#: Location relative to the repository root.
#: previous value 2026-09-03 (D2 era): "docs/evidence/GT_PS2_STUDY_SURROGATE_PACK_2026-09-03.json"
#: previous value 2026-09-02: "docs/evidence/GT_PS2_STUDY_SURROGATE_PACK_2026-09-02.json"
DEFAULT_PACK_PATH = "docs/evidence/GT_PS2_STUDY_SURROGATE_PACK_2026-10-09.json"


class PackEraMismatch(RuntimeError):
    """The pack was fitted against a different pinned-props era."""


def default_pack_path() -> Path:
    """Absolute path to the shipped pack, derived from this file's location."""

    # src/dtdc_simulator/core2/study_lane/pack.py -> repository root
    root = Path(__file__).resolve().parents[4]
    return root / DEFAULT_PACK_PATH


def canonical_pack_digest(payload: dict[str, Any]) -> str:
    """``sha256`` over the pack with its own ``pack_digest`` field removed."""

    reduced = {key: value for key, value in payload.items() if key != "pack_digest"}
    encoded = json.dumps(reduced, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class SurrogatePack:
    """A loaded, era-verified surrogate pack."""

    path: Path
    payload: dict[str, Any]
    pack_digest: str

    @property
    def record(self) -> str:
        return str(self.payload["record"])

    @property
    def surrogates(self) -> dict[str, Any]:
        return self.payload["surrogates"]

    @property
    def acceptance(self) -> dict[str, Any]:
        return self.payload["acceptance"]

    def adopted(self) -> dict[str, Any]:
        """Only the surrogates whose measured certificate passed acceptance."""

        return {
            name: entry
            for name, entry in self.surrogates.items()
            if bool(entry.get("adopted"))
        }

    def declaration(self) -> dict[str, Any]:
        """The self-declaration every study-lane evidence record must carry."""

        adopted = self.adopted()
        return {
            "lane": "study",
            "surrogate_pack_record": self.record,
            "surrogate_pack_digest": self.pack_digest,
            "surrogate_pack_path": self.payload.get("pack_path", DEFAULT_PACK_PATH),
            "certificate_reference": self.payload["certificate_reference"],
            "certificate_bounds": self.payload["certificate_bounds"],
            "props_pin_digests": self.payload["props_pin_digests"],
            "adopted_surrogates": sorted(adopted),
            "declined_surrogates": sorted(set(self.surrogates) - set(adopted)),
        }


def load_pack(path: Path | str | None = None, *, verify_era: bool = True) -> SurrogatePack:
    """Read a pack, verify its digest and (by default) its props era-binding."""

    resolved = Path(path) if path is not None else default_pack_path()
    payload = json.loads(resolved.read_text(encoding="utf-8"))

    recomputed = canonical_pack_digest(payload)
    stored = payload.get("pack_digest")
    if stored != recomputed:
        raise PackEraMismatch(
            f"surrogate pack {resolved} is internally inconsistent: "
            f"stored digest {stored} != recomputed {recomputed}"
        )

    if verify_era:
        expected = payload["props_pin_digests"]
        actual = props_pin_digests()
        if expected != actual:
            differing = sorted(
                name
                for name in set(expected) | set(actual)
                if expected.get(name) != actual.get(name)
            )
            raise PackEraMismatch(
                "surrogate pack was fitted against a different pinned-props era; "
                f"re-run scripts/build_study_surrogates.py.  Changed surfaces: {differing}"
            )

    return SurrogatePack(path=resolved, payload=payload, pack_digest=recomputed)
