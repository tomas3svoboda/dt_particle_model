"""Exact-argument memoization for pure thermodynamic state functions.

Owner-approved 2026-08-30 (GT_PS2_EXACT_MEMOIZATION_RULING): the decorated
functions are pure maps from exact binary64/frozen-dataclass arguments —
PLUS the declared ambient datum anchors — to exact results, so serving a
repeat from the cache returns the SAME BITS the computation would produce;
precision and robustness are untouched by construction, and the
accompanying tests assert bit-identity against the uncached originals.

Measured motivation (one march interval, 2026-08-30): 85,720 calls to
``wet_core.specific_energy`` carried only 954 unique argument sets (98.9 %
exact repeats); the four decorated functions dominate ~130 s of a ~180 s
interval.

Rules of the construction:
- keys are the exact argument tuples (kwargs sorted by name) EXTENDED by
  the ``ambient`` tuple — the declared module-level datum anchors the
  function's result also depends on (e.g. the hexane common-datum
  references, which the datum-covariance suites legitimately monkeypatch);
  an undeclared ambient dependency is a defect, not a tuning choice;
- an unhashable argument silently bypasses the cache (never fails);
- exceptions are never cached — a refusing input re-runs and re-raises
  identically;
- the cache is bounded; at the bound, new entries are simply not stored
  (no eviction reordering, no nondeterminism);
- ``function.__wrapped__`` keeps the uncached original for equivalence
  tests, and ``function.exact_memo_cache`` exposes the cache for tests.
"""

from __future__ import annotations

import functools
from typing import Callable

_MISS = object()
_CACHE_ENTRY_LIMIT = 2_000_000


def exact_memo(*, ambient: Callable[[], tuple] | None = None):
    def decorate(function):
        cache: dict = {}

        @functools.wraps(function)
        def wrapped(*args, **kwargs):
            try:
                key = (
                    args,
                    tuple(sorted(kwargs.items())) if kwargs else (),
                    ambient() if ambient is not None else (),
                )
                hit = cache.get(key, _MISS)
            except TypeError:
                return function(*args, **kwargs)
            if hit is not _MISS:
                return hit
            value = function(*args, **kwargs)
            if len(cache) < _CACHE_ENTRY_LIMIT:
                cache[key] = value
            return value

        wrapped.exact_memo_cache = cache
        return wrapped

    return decorate
