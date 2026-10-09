"""The shared bounded exact-argument memo on a packet caloric envelope inversion.

WHAT THIS IS.  One reusable copy of the construction that already wraps the
tag-1 combined caloric inversion
(``high_loading_packet_adapter._exact_envelope_memo``, speed work 2026-09-02,
following the owner-approved exact-memoization precedent in
``particle/_exact_memo.py``).  Owner ruling 2026-09-06 ("I agree with all 4
points.  Feel free to apply them and proceed", item 4 of the science lane's
four cheap items): memoize the dry-sorbate caloric envelope inversion exactly
the way the tag-1 inversion is memoized, bit-identical by construction.

WHY IT IS A SEPARATE MODULE.  ``high_loading_packet_adapter.py`` is a
BYTE-PINNED source module (``scripts/rs2_pin_drift_sweep.py`` lists it among
the pinned files), so its own decorator could not be generalized in place, and
its private ``_envelope_memo_key`` hard-codes the tag-1 keyword names -- it
would silently BYPASS on a tag-2 or tag-3 argument set.  This module therefore
carries the SAME construction parameterized by the variant's own keyword order,
and imports the two behaviour constants (entry limit, disable variable) FROM
the tag-1 module so they can never drift apart and so no new numeric literal is
introduced.  ``tests/test_core2_envelope_inversion_memo_parity.py`` asserts
that the generic key, driven with the tag-1 keyword order, is EQUAL to the
tag-1 module's own private key for the same adapter and arguments.

Rules of the construction (unchanged from the ruled precedent):

- the key is the EXACT binary64 image of every float argument, plus the
  complete set of context pins the result depends on (adapter class,
  configuration digest, energy datum, fixed pressure bits, the component datum
  adapter's definition digest, the sphere/material parameter repr -- which is
  what ``sphere._material_wet_params`` reads -- and both live caloric datum
  signatures), so a monkeypatched datum MISSES instead of being served a stale
  value;
- an argument that is not an exact ``float`` bypasses the cache entirely
  (packing would merge ``1`` with ``1.0``);
- REFUSALS ARE NEVER CACHED: every typed refusal, including
  ``DrySorbateReWettingRefusal``, ``DrySorbateExhaustedRefusal`` and
  ``DrySorbateCaloricBracketError`` (objection O28 of the audit refuter,
  ledger ``GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md``), propagates uncached and
  re-raises identically on the next call;
- the cache is bounded and evicts the oldest entry, so a long march stays
  memory-safe;
- ``wrapped.__wrapped__`` keeps the uncached original for bit-identity
  equivalence tests, and ``cache_info()`` exposes the SAME field names the
  march driver already reports for tag-1.

No gate, validator, or decode verdict is cached here; the owner declined that
class of caching.  Every issue path still runs its full validation.
"""

from __future__ import annotations

import functools
import os

from . import high_loading_packet_adapter as tag1
from .exact_cache import float_bits_key_flat
from .props import hexane, water

#: Bound FROM the tag-1 module, not restated, so the two can never disagree.
_ENVELOPE_MEMO_ENTRY_LIMIT = tag1._ENVELOPE_MEMO_ENTRY_LIMIT  # noqa: SLF001
#: The same escape hatch variable the tag-1 memo publishes: exactly "1" runs
#: the uncached original.
ENVELOPE_MEMO_DISABLE_ENVIRONMENT_VARIABLE = (
    tag1.ENVELOPE_MEMO_DISABLE_ENVIRONMENT_VARIABLE
)
_ENVELOPE_MEMO_ENABLED = os.environ.get(ENVELOPE_MEMO_DISABLE_ENVIRONMENT_VARIABLE) != "1"
_ENVELOPE_MEMO_MISS = object()


def envelope_memo_key(adapter, keywords: dict, argument_order: tuple[str, ...]) -> tuple | None:
    """Return the exact memo key, or ``None`` to bypass the memo entirely."""

    if len(keywords) != len(argument_order):
        return None
    try:
        arguments = tuple(keywords[name] for name in argument_order)
    except KeyError:
        return None
    if any(type(value) is not float for value in arguments):
        return None
    pressure_pa = adapter.pressure_pa
    if type(pressure_pa) is not float:
        return None
    try:
        return (
            type(adapter).__module__,
            type(adapter).__qualname__,
            adapter.configuration_digest,
            adapter.energy_datum_id,
            adapter.component_datum_adapter.definition_digest,
            repr(adapter.sphere_params),
            hexane.caloric_datum_signature(),
            water.caloric_datum_signature(),
            float_bits_key_flat(pressure_pa, *arguments),
        )
    except Exception:
        # A context pin that cannot even be read is never a cache hit; the
        # authority itself decides what to do with the malformed adapter.
        return None


def exact_envelope_memo(*argument_order: str):
    """Return the decorator that wraps one pure caloric inversion in the memo.

    ``argument_order`` is the variant's own keyword-only argument order; it
    fixes the packing order of the flat binary64 image and nothing else.
    """

    def decorate(function):
        cache: dict = {}
        statistics = {
            "hits": 0,
            "misses": 0,
            "bypasses": 0,
            "evictions": 0,
            "uncached_refusals": 0,
        }

        @functools.wraps(function)
        def wrapped(self, **keywords):
            if not _ENVELOPE_MEMO_ENABLED:
                statistics["bypasses"] += 1
                return function(self, **keywords)
            key = envelope_memo_key(self, keywords, argument_order)
            if key is None:
                statistics["bypasses"] += 1
                return function(self, **keywords)
            stored = cache.get(key, _ENVELOPE_MEMO_MISS)
            if stored is not _ENVELOPE_MEMO_MISS:
                statistics["hits"] += 1
                return stored
            try:
                value = function(self, **keywords)
            except BaseException:
                statistics["uncached_refusals"] += 1
                raise
            statistics["misses"] += 1
            if len(cache) >= _ENVELOPE_MEMO_ENTRY_LIMIT:
                del cache[next(iter(cache))]
                statistics["evictions"] += 1
            cache[key] = value
            return value

        def cache_info() -> dict:
            return {
                **statistics,
                "currsize": len(cache),
                "maxsize": _ENVELOPE_MEMO_ENTRY_LIMIT,
                "enabled": _ENVELOPE_MEMO_ENABLED,
            }

        def cache_clear() -> None:
            cache.clear()
            for name in statistics:
                statistics[name] = 0

        wrapped.cache_info = cache_info
        wrapped.cache_clear = cache_clear
        wrapped.envelope_memo_cache = cache
        wrapped.envelope_memo_argument_order = argument_order
        return wrapped

    return decorate


__all__ = [
    "ENVELOPE_MEMO_DISABLE_ENVIRONMENT_VARIABLE",
    "envelope_memo_key",
    "exact_envelope_memo",
]
