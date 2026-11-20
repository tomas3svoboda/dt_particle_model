"""Kernel-envelope DRY_FACE_FLUX_MEMO: the pinned dry face flux memoized on exact bits within a solve.

The lead's Class-A numerics of 2026-11-02 (record BIRTH_COST_LEVERS_RECORD), measured on a tower-born
particle's first plant seconds (the record's section 2): inside one domain-aware Newton solve, 72 % of
the calls of the dry face flux (``cut_transport._dry_flux``, one of the sixteen pinned S1-parity files,
not edited) repeat an earlier call of the same solve with bit-identical arguments, because a dense
finite-difference Jacobian column moves one unknown and leaves every face whose two cells it does not
touch exactly as it was.  Armed, the flux reached through ``core.cut._dry_flux`` is wrapped, at the
innermost level, by an exact memo:

* the scope is one Newton solve (``particle_numerical_repairs.domain_aware_newton`` opens a frame on
  entry and drops it on exit; nested solves get their own frame); outside every solve the flux runs
  its former statements;
* the key is the exact IEEE-754 bits of every argument (floats packed with their type, enums by name,
  frozen dataclasses field by field) together with the identity of every in-memory substitution the
  harness scopes make inside the flux's call tree (the proof allowance's point evaluator, the segment
  classifier, the activity box, the interface-composition solve, the boundary validator, the cached
  equilibrium, the sorption activity and retained laws, the equilibrium states, the stress-multiplier
  table) and the K1 suppression depth; an argument of any other type makes the call unkeyable (run as
  before, counted);
* a hit returns the identical result (the very object when its whole graph is immutable: frozen
  dataclasses, tuples, numbers, strings and enums; a deep copy otherwise) and replays the proof
  allowance's accounting the
  call made (points started, memoized points, the last path and point), exactly as the repeated call
  would book them: the repeated call evaluates the same points, every one already in the allowance's
  own point memo; a hit whose points would cross the allowance's point clause, or any hit while the
  declared wall-clock clause is armed, runs the flux instead (the clause then acts as before);
* a call that raises is not memoized; the memo holds at most FACE_FLUX_MEMO_BOUND entries per frame
  (least recently used evicted, counted).

Unchanged: every flux value, every proof verdict and point count, every gate, tolerance, budget and
allowance.  Off (the knob absent), nothing is installed and nothing is kept.  The shadow audit
(``SHADOW_AUDIT``, set only by tests and measurement scripts, never by a knob) runs the flux beside
every hit and raises FaceFluxMemoDisagreement on any differing bit of the result or of the replayed
accounting.

physically_qualifying: false.
"""

from __future__ import annotations

from collections import OrderedDict
from contextlib import contextmanager
import copy
import dataclasses
import enum
import struct

import numpy as np

import particle_kernel_envelope as envelope

BIRTH_COST_LEVERS_RECORD = "docs/GT_PS2_DTDC_KERNEL_BIRTH_COST_LEVERS_2026-11-02.md"
# The bound per solve frame: above the largest frame measured (the record's section 4: a born job's
# birth solve keeps 183 residual evaluations of at most 6 faces), so no measured solve evicts.
FACE_FLUX_MEMO_BOUND = 4096  # WAS absent until 2026-11-02 (no memo)
# Test and measurement hook only (never a knob): run the flux beside every hit and compare.
SHADOW_AUDIT = False

_PACK = struct.Struct("<d").pack
_STATE = {"memo": None}


def _immutable(value, seen):
    """True when the whole object graph is immutable (frozen dataclasses with no attribute beyond
    their fields, tuples and frozensets of such, numbers, strings, bytes, enums, numpy scalars)."""
    kind = type(value)
    if value is None or kind in (bool, int, float, str, bytes, complex):
        return True
    if isinstance(value, (enum.Enum, np.generic)):
        return True
    marker = id(value)
    if marker in seen:
        return True
    if kind in (tuple, frozenset):
        seen.add(marker)
        return all(_immutable(item, seen) for item in value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        if not kind.__dataclass_params__.frozen:
            return False
        names = {field.name for field in dataclasses.fields(value)}
        if hasattr(value, "__dict__") and not set(vars(value)) <= names:
            return False
        seen.add(marker)
        return all(_immutable(getattr(value, name), seen) for name in names)
    return False


class FaceFluxMemoDisagreement(AssertionError):
    """The shadow audit found a hit whose result or accounting differs from the flux's own."""


class _Unkeyable(Exception):
    pass


def _canon(value, cache):
    if value is None or value is True or value is False:
        return value
    kind = type(value)
    if kind is float:
        return ("f", _PACK(value))
    if kind is int:
        return ("i", value)
    if kind is str:
        return ("s", value)
    if isinstance(value, enum.Enum):
        return ("e", kind.__qualname__, value.name)
    if isinstance(value, np.floating):
        return ("nf", kind.__name__, _PACK(float(value)))
    if kind is tuple:
        return ("t",) + tuple(_canon(item, cache) for item in value)
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        if not kind.__dataclass_params__.frozen:
            raise _Unkeyable(kind.__qualname__)
        held = cache.get(id(value))
        if held is not None and held[0] is value:
            return held[1]
        canonical = ("dc", kind.__module__, kind.__qualname__) + tuple(
            (field.name, _canon(getattr(value, field.name), cache))
            for field in dataclasses.fields(value)
        )
        cache[id(value)] = (value, canonical)
        return canonical
    raise _Unkeyable(kind.__qualname__)


class FaceFluxMemo:
    """The installed memo: its frames (one per open Newton solve) and its counters."""

    def __init__(self, cut, acfp, substitutions, suppression):
        self.cut = cut
        self.acfp = acfp
        self.substitutions = substitutions
        self.suppression = suppression
        self.frames = []
        self.canonical = {}
        evaluator = acfp._PathActivities._evaluate_coordinates
        # the native point evaluator (no proof allowance installed at the memo's installation)
        self.native_evaluator = None if hasattr(evaluator, "proof_allowance") else evaluator
        self.counts = {
            "calls": 0,
            "calls_outside_solves": 0,
            "hits": 0,
            "misses": 0,
            "unkeyable": 0,
            "raised_not_memoized": 0,
            "hits_declined_at_point_clause": 0,
            "hits_declined_wall_clause": 0,
            "evictions": 0,
            "frames": 0,
            "replayed_points": 0,
            "shadow_audited_hits": 0,
            "mutable_results_copied": 0,
        }

    # -------------------------------------------------------------- frames

    @contextmanager
    def frame(self):
        self.frames.append(OrderedDict())
        self.counts["frames"] += 1
        try:
            yield
        finally:
            self.frames.pop()
            if not self.frames:
                self.canonical.clear()

    # -------------------------------------------------------------- the wrapper

    def _environment(self):
        evaluator = self.acfp._PathActivities._evaluate_coordinates
        return evaluator, tuple(get() for get in self.substitutions)

    @staticmethod
    def _allowance_snapshot(allowance):
        if allowance is None:
            return None
        audit = allowance.audit
        return (
            audit["points_started"],
            audit.get("memoized_points"),
            audit.get("last_path"),
            audit.get("last_temperature_k"),
            audit.get("last_y_hexane"),
        )

    def wrap(self, original):
        memo = self

        def memoized_dry_flux(*args, **kwargs):
            memo.counts["calls"] += 1
            if not memo.frames:
                memo.counts["calls_outside_solves"] += 1
                return original(*args, **kwargs)
            evaluator, held = memo._environment()
            allowance = getattr(evaluator, "proof_allowance", None)
            if allowance is None and evaluator is not memo.native_evaluator:
                # an evaluator this module cannot account for: the flux runs as before
                memo.counts["unkeyable"] += 1
                return original(*args, **kwargs)
            try:
                key = (
                    tuple(id(obj) for obj in held),
                    memo.suppression["depth"],
                    id(evaluator),
                    _canon(args, memo.canonical),
                    tuple(sorted((name, _canon(v, memo.canonical)) for name, v in kwargs.items())),
                )
            except _Unkeyable:
                memo.counts["unkeyable"] += 1
                return original(*args, **kwargs)
            frame = memo.frames[-1]
            entry = frame.get(key)
            if (
                entry is not None
                and entry[0] is evaluator
                and all(a is b for a, b in zip(entry[1], held))
            ):
                _, _, result, points, final, _ = entry
                replay = memo._replayable(allowance, points)
                if replay:
                    if SHADOW_AUDIT:
                        memo._shadow(original, args, kwargs, allowance, result, points, final)
                    memo._replay(allowance, points, final)
                    frame.move_to_end(key)
                    memo.counts["hits"] += 1
                    memo.counts["replayed_points"] += points
                    return result if entry[5] else copy.deepcopy(result)
            before = memo._allowance_snapshot(allowance)
            try:
                result = original(*args, **kwargs)
            except BaseException:
                memo.counts["raised_not_memoized"] += 1
                raise
            memo.counts["misses"] += 1
            after = memo._allowance_snapshot(allowance)
            points = 0 if allowance is None else after[0] - before[0]
            final = None if allowance is None else after[2:]
            # the entry holds the evaluator and the substitutions themselves, so an identity
            # matched on a hit is never a recycled id
            shared = _immutable(result, set())
            if not shared:
                memo.counts["mutable_results_copied"] += 1
            frame[key] = (
                evaluator,
                held,
                result if shared else copy.deepcopy(result),
                points,
                final,
                shared,
            )
            if len(frame) > FACE_FLUX_MEMO_BOUND:
                frame.popitem(last=False)
                memo.counts["evictions"] += 1
            return result

        memoized_dry_flux.face_flux_memo = memo
        memoized_dry_flux.original = original
        return memoized_dry_flux

    def _replayable(self, allowance, points):
        if allowance is None or points == 0:
            return True
        if allowance.wall_clause:
            self.counts["hits_declined_wall_clause"] += 1
            return False
        # the repeated call would refuse at its (k + 1)-th point if points_started + k reached the
        # point clause for some k < points; then the flux itself runs and refuses as before
        if allowance.audit["points_started"] + points > allowance.point_budget:
            self.counts["hits_declined_at_point_clause"] += 1
            return False
        return True

    @staticmethod
    def _replay(allowance, points, final):
        if allowance is None or points == 0:
            return
        audit = allowance.audit
        audit["points_started"] += points
        if allowance.memoize_identical_points:
            # every point of the repeated path is already in this allowance's point memo
            audit["memoized_points"] += points
        audit["last_path"], audit["last_temperature_k"], audit["last_y_hexane"] = final

    def _shadow(self, original, args, kwargs, allowance, result, points, final):
        self.counts["shadow_audited_hits"] += 1
        before = None if allowance is None else dict(allowance.audit)
        fresh = original(*args, **kwargs)
        if not _same_bits(fresh, result):
            raise FaceFluxMemoDisagreement("a memo hit's flux differs from the flux's own")
        if allowance is not None:
            audit = allowance.audit
            expected = dict(before)
            if points:
                expected["points_started"] += points
                if allowance.memoize_identical_points:
                    expected["memoized_points"] += points
                expected["last_path"], expected["last_temperature_k"], expected["last_y_hexane"] = (
                    final
                )
            if not _same_bits(expected, dict(audit)):
                raise FaceFluxMemoDisagreement(
                    "a memo hit's replayed proof accounting differs from the flux's own"
                )
            # undo the shadow call's booking; the replay books it once
            audit.clear()
            audit.update(before)

    def stamp(self):
        return {
            **self.counts,
            "bound_per_frame": FACE_FLUX_MEMO_BOUND,
            "key": "exact bits of every argument with the identity of the harness's in-memory "
            "substitutions in the flux's call tree and the K1 suppression depth",
            "scope": "one domain-aware Newton solve",
            "cumulative": "since the memo's installation",
            "record": BIRTH_COST_LEVERS_RECORD,
            "physically_qualifying": False,
        }


def _same_bits(a, b):
    try:
        return repr(_canon_any(a)) == repr(_canon_any(b))
    except RecursionError:
        return False


def _canon_any(value):
    if isinstance(value, float):
        return ("f", type(value).__name__, _PACK(float(value)))
    if isinstance(value, (list, tuple)):
        return (type(value).__name__,) + tuple(_canon_any(v) for v in value)
    if isinstance(value, dict):
        return ("d",) + tuple((repr(k), _canon_any(v)) for k, v in value.items())
    if isinstance(value, np.ndarray):
        return ("nd", value.dtype.str, value.shape, value.tobytes())
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return ("dc", type(value).__qualname__) + tuple(
            (f.name, _canon_any(getattr(value, f.name))) for f in dataclasses.fields(value)
        )
    return ("r", repr(value))


def active():
    """The installed memo, or None."""
    return _STATE["memo"]


@contextmanager
def installed(environ=None):
    """Install the memo around ``core.cut._dry_flux`` for the scope when DRY_FACE_FLUX_MEMO is armed
    (nothing otherwise); an enclosing installation is kept as it is.  Enter it before any harness
    scope captures the flux (a run's or a job's start)."""
    if not envelope.armed(envelope.DRY_FACE_FLUX_MEMO, environ) or _STATE["memo"] is not None:
        yield _STATE["memo"]
        return
    import experiment_core2_radial_fast_wetting as wetting
    import native_water_path_bounds as path_bounds  # noqa: F401  (the activity box's owner)
    import particle_numerical_repairs as numerics

    core = wetting.dry.core
    cut, acfp = core.cut, core.cut.acfp
    substitutions = (
        lambda: acfp._active_set_on_segment,
        lambda: acfp._activity_box,
        lambda: cut._solve_interface_composition_with_audit,
        lambda: cut._validate_dry_boundary,
        lambda: core.cp._evaluate_equilibrium_cached,
        lambda: core.sp.water_activity,
        lambda: core.sp.water_retained,
        lambda: wetting.dry.ct._equilibrium_states,
        lambda: core.efr._ALLOWED_STRESS_MULTIPLIERS,
    )
    memo = FaceFluxMemo(cut, acfp, substitutions, wetting._CAP_SUPPRESS)
    original = cut._dry_flux
    wrapper = memo.wrap(original)
    cut._dry_flux = wrapper
    _STATE["memo"] = memo
    numerics._FACE_FLUX_MEMO = memo
    try:
        yield memo
    finally:
        numerics._FACE_FLUX_MEMO = None
        _STATE["memo"] = None
        if cut._dry_flux is wrapper:
            cut._dry_flux = original


def stamp(environ=None):
    """The installed memo's counters (None when DRY_FACE_FLUX_MEMO is not armed or none installed)."""
    memo = _STATE["memo"]
    if memo is None or not envelope.armed(envelope.DRY_FACE_FLUX_MEMO, environ):
        return None
    return memo.stamp()
