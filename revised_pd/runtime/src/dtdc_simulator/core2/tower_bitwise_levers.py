"""Bitwise-identical work-avoidance levers of the engineering tower march.

WHAT THIS IS.  The shared machinery of the tower-march speed build of
2026-10-31 (record ``GT_PS2_DTDC_TOWER_MARCH_BITWISE_LEVERS_BUILD_2026-10-31.md``;
the levers were measured as recorded variant trees in
``GT_PS2_DTDC_TOWER_MARCH_LEVERS_MEASURED_2026-10-31.md``).  Every lever here
returns the bits the former path returns and refuses exactly where the former
path refused; none of them moves a law, a gate, a tolerance, a budget, a
declared value or a codec exactness check.  They avoid repeating work whose
result is already known:

- an exact, bounded memo of the codec's encoder (``encode_state``) for the
  tag-2 and tag-3 packet adapters (lever "encode memo");
- an exact, bounded memo of the native zero-gauge component-datum adapter's
  definition digest, read at every tag-2/tag-3 issue, decode and envelope-memo
  key (lever "datum digest memo");
- one complete-payload validation per issued tag-2/tag-3 payload instead of
  two, with the former refusal order restored on every refusing path (lever
  "single validation", :func:`encode_and_decode_validated`);
- the tag-2/tag-3 caloric inversion's bisection replayed from a certified
  bracket (lever "replayed bisection", :func:`replayed_bisection`), with the
  sampled energies checked monotone on every replay and the former bisection
  run on any failure;
- the bounded exact memos of the pre-desolventizer level-0 saturated state
  (``pd_k_dome_bypass_tray_integration._level0_packet_at_saturation``): the
  seed energy and the advanced candidate, shared across packets and solves;
- the boiling root's value-only residual (owner ruling of 2026-10-31, record
  ``GT_PS2_DTDC_TOWER_BOILING_ROOT_VALUE_ONLY_RESIDUAL_BUILD_2026-10-31.md``):
  inside the level-0 boiling root's false-position and bisection loops the
  residual is evaluated from the saturated seed's energy without building the
  candidate packet, and the accepted end state is then built and validated
  exactly as the former path built it.  This is the one lever here that is
  NOT a pure refactor: an interior root iterate is no longer required to be a
  codec-admissible state (the ruling).  Its bit-identity on the measured
  marches is checked by the march comparisons, and a sampled typed audit
  (``DTDC_TOWER_BOILING_ROOT_AUDIT_EVERY_N``, :func:`boiling_root_audit_due`,
  :func:`record_boiling_root_audit`) re-runs the former full-candidate
  residual, memos bypassed, at the two final ends of every N-th root and
  refuses on any difference.

WHAT IS NOT HERE.  The byte-pinned source modules (the RS-2 pin sweep's rows
and the sixteen S1-parity files) are not edited by this build, so the levers
that live inside them (the tag-1 halves of the validation and bisection
levers, the encoder and digest memos inside their own modules, the sorption
quadrature batch in ``particle/wet_core.py`` and the Law-2 cell solve reuse in
``cell_engineering_feasibility.py``) are not built; the record says why.

THE RULES OF EVERY MEMO HERE (the envelope-memo precedent,
``envelope_inversion_memo.py``): the key is the exact binary64 image of every
float argument plus every context pin the result depends on, an object the
result depends on by identity is kept in the entry and served only to that same
object, an argument that is not of an exact keyable type bypasses the memo, a
refusal is never stored (it propagates and re-raises from a fresh evaluation),
the size is bounded and the oldest entry evicted, and the hits, misses,
bypasses, evictions and uncached refusals are counted and stamped by the march
(:func:`lever_statistics`).  No decode verdict is cached (the owner declined
that class).

THE SWITCH.  ``DTDC_TOWER_BITWISE_LEVERS_DISABLE`` set to exactly ``"1"`` (read
once at import) runs every former path, for A/B timing only.  Absent, the
levers are the shipped path.  It also runs the former boiling root (every
candidate built); asking for the boiling-root audit with the switch set is
refused, because there is no value-only root to audit.

physically_qualifying is False: a speed build, not a physics claim.
"""

from __future__ import annotations

import os
import struct
from collections.abc import Callable, Hashable
from typing import Any

from . import component_energy_datum_adapter as component_datum
from . import state_codec
from .props import hexane
from .props import water

physically_qualifying = False

#: Exactly "1" runs every former path (A/B timing only); anything else, or
#: absent, runs the levers.  Read once at import, like the envelope memo's
#: own disable variable.
TOWER_BITWISE_LEVERS_DISABLE_ENVIRONMENT_VARIABLE = "DTDC_TOWER_BITWISE_LEVERS_DISABLE"
LEVERS_ENABLED = os.environ.get(TOWER_BITWISE_LEVERS_DISABLE_ENVIRONMENT_VARIABLE) != "1"

#: The boiling root's sampled audit (owner ruling of 2026-10-31,
#: GT_PS2_DTDC_TOWER_BOILING_ROOT_VALUE_ONLY_RESIDUAL_BUILD_2026-10-31.md): a
#: positive decimal integer N audits every N-th value-only root; absent, no
#: audit.  Read once at import, like the disable switch; any other value is
#: refused (typed) at the first boiling root, not clamped or ignored.
BOILING_ROOT_AUDIT_EVERY_N_ENVIRONMENT_VARIABLE = "DTDC_TOWER_BOILING_ROOT_AUDIT_EVERY_N"
BOILING_ROOT_AUDIT_EVERY_N_RAW = os.environ.get(BOILING_ROOT_AUDIT_EVERY_N_ENVIRONMENT_VARIABLE)


class BoilingRootAuditKnobRefused(ValueError):
    """The audit knob is not a positive decimal integer, or is set with the levers disabled."""


class BoilingRootAuditDisagreement(RuntimeError):
    """The value-only boiling root and the former full-candidate residual differ."""


# Memo sizes are entry COUNTS (memory bounds), not numerical gates: a bound
# changes how often a value is recomputed, never which bits are returned.  The
# saturated-state, seed-energy and encode bounds are the 16,384 entries the
# measured variants used (GT_PS2_DTDC_TOWER_MARCH_LEVERS_MEASURED_2026-10-31.md,
# variants V1A, V1VA and V3E); the digest memo holds one entry per live
# component-datum configuration, of which a march has one.
#: Entry bound of the encode memo (GT_PS2_DTDC_TOWER_MARCH_BITWISE_LEVERS_BUILD_2026-10-31.md).
ENCODE_MEMO_ENTRY_LIMIT = 16384
#: Entry bound of the datum-digest memo (GT_PS2_DTDC_TOWER_MARCH_BITWISE_LEVERS_BUILD_2026-10-31.md).
DATUM_DIGEST_MEMO_ENTRY_LIMIT = 64
#: Entry bound of each level-0 saturation memo (GT_PS2_DTDC_TOWER_MARCH_BITWISE_LEVERS_BUILD_2026-10-31.md).
SATURATION_MEMO_ENTRY_LIMIT = 16384

# The replayed bisection's two declared numerical settings, measured as variant
# V2B in GT_PS2_DTDC_TOWER_MARCH_LEVERS_MEASURED_2026-10-31.md section 4: the
# certified bracket is narrowed by Illinois false position to this width before
# the replay, and a shipped midpoint at least this margin outside the certified
# bracket is decided without evaluating it.  Neither is a gate on the result:
# the returned temperature is the former bisection's, and every replay whose
# sampled energies are not monotone runs the former bisection instead.
#: Certified-bracket width, K (GT_PS2_DTDC_TOWER_MARCH_BITWISE_LEVERS_BUILD_2026-10-31.md).
REPLAY_CERTIFIED_BRACKET_WIDTH_K = 1.0e-8
#: Decision margin outside the certified bracket, K (GT_PS2_DTDC_TOWER_MARCH_BITWISE_LEVERS_BUILD_2026-10-31.md).
REPLAY_DECISION_MARGIN_K = 1.0e-9

_MISS = object()


class BoundedExactMemo:
    """One bounded exact-key memo with the envelope-memo accounting.

    ``get_or_compute(key, compute, anchors)`` returns the stored value for an
    exact ``key`` whose stored anchors are the SAME objects (``is``) as
    ``anchors``; otherwise it calls ``compute()``, stores a returned value and
    returns it.  A ``None`` key bypasses the memo.  An exception from
    ``compute`` is never stored.
    """

    __slots__ = ("name", "limit", "_entries", "_statistics")

    def __init__(self, name: str, limit: int) -> None:
        self.name = name
        self.limit = limit
        self._entries: dict = {}
        self._statistics = {
            "hits": 0,
            "misses": 0,
            "bypasses": 0,
            "evictions": 0,
            "uncached_refusals": 0,
            "anchor_mismatches": 0,
        }

    def get_or_compute(
        self,
        key: Hashable | None,
        compute: Callable[[], Any],
        anchors: tuple = (),
    ) -> Any:
        statistics = self._statistics
        if key is None:
            statistics["bypasses"] += 1
            return compute()
        stored = self._entries.get(key, _MISS)
        if stored is not _MISS:
            stored_anchors, value = stored
            if len(stored_anchors) == len(anchors) and all(
                left is right for left, right in zip(stored_anchors, anchors)
            ):
                statistics["hits"] += 1
                return value
            statistics["anchor_mismatches"] += 1
            del self._entries[key]
        try:
            value = compute()
        except BaseException:
            statistics["uncached_refusals"] += 1
            raise
        statistics["misses"] += 1
        if len(self._entries) >= self.limit:
            del self._entries[next(iter(self._entries))]
            statistics["evictions"] += 1
        self._entries[key] = (anchors, value)
        return value

    def statistics(self) -> dict:
        return {
            **self._statistics,
            "currsize": len(self._entries),
            "maxsize": self.limit,
            "enabled": LEVERS_ENABLED,
        }

    def clear(self) -> None:
        self._entries.clear()
        for name in self._statistics:
            self._statistics[name] = 0


ENCODE_MEMO = BoundedExactMemo("state_codec_encode_tag2_tag3", ENCODE_MEMO_ENTRY_LIMIT)
DATUM_DIGEST_MEMO = BoundedExactMemo(
    "native_zero_gauge_datum_definition_digest", DATUM_DIGEST_MEMO_ENTRY_LIMIT
)
SATURATION_SEED_ENERGY_MEMO = BoundedExactMemo(
    "pd_level0_saturation_seed_energy", SATURATION_MEMO_ENTRY_LIMIT
)
SATURATION_STATE_MEMO = BoundedExactMemo("pd_level0_saturated_state", SATURATION_MEMO_ENTRY_LIMIT)
#: The value-only seed energies of the boiling root's interior iterates
#: (GT_PS2_DTDC_TOWER_BOILING_ROOT_VALUE_ONLY_RESIDUAL_BUILD_2026-10-31.md).  Kept
#: apart from the issued-seed memo above, so that a value computed without the
#: tag-2 issue is never served where the former path issues the seed.
SATURATION_SEED_ENERGY_VALUE_ONLY_MEMO = BoundedExactMemo(
    "pd_level0_saturation_seed_energy_value_only", SATURATION_MEMO_ENTRY_LIMIT
)


def float_bits(*values: float) -> tuple[bytes, ...] | None:
    """The exact binary64 images, or ``None`` when any value is not an exact float."""

    if not all(type(value) is float for value in values):
        return None
    return tuple(struct.pack("<d", value) for value in values)


def caloric_context() -> tuple:
    """The live context pins a packet state depends on beyond its codec and arguments.

    Both caloric datum signatures (a monkeypatched datum misses instead of
    being served a stale value, as in the envelope memo) and the exactness
    knob the tag-1/tag-2 validations read at call time.
    """

    from . import high_loading_packet_adapter as tag1

    return (
        hexane.caloric_datum_signature(),
        water.caloric_datum_signature(),
        tag1.exactness_mass_primary_armed(),
    )


# -- the encode memo ----------------------------------------------------------


def _codec_value_key(value: object) -> tuple | None:
    kind = type(value)
    if kind is float:
        return ("f", struct.pack("<d", value))
    if kind is str:
        return ("s", value)
    if kind is int:
        return ("i", value)
    if kind is bytes:
        return ("b", value)
    return None


def encode_state(schema: state_codec.StateSchema, values: tuple) -> bytes:
    """``state_codec.encode_state`` through the exact encode memo.

    The encoder is a deterministic function of the (frozen) schema object and
    the value tuple; the key is the schema's identity (the schema kept in the
    entry) and the exact image of every value (a float by its binary64 bits,
    so ``0.0`` and ``-0.0`` differ).  A value of any other type than float,
    str, int or bytes bypasses the memo.
    """

    if not LEVERS_ENABLED:
        return state_codec.encode_state(schema, values)
    key: tuple | None = None
    if type(values) is tuple:
        parts = []
        for value in values:
            part = _codec_value_key(value)
            if part is None:
                parts = None
                break
            parts.append(part)
        if parts is not None:
            key = (id(schema), tuple(parts))
    return ENCODE_MEMO.get_or_compute(
        key, lambda: state_codec.encode_state(schema, values), anchors=(schema,)
    )


# -- the datum-digest memo ----------------------------------------------------


def component_datum_definition_digest(adapter: object) -> str:
    """``adapter.definition_digest``, computed once per exact configuration.

    Only the native zero-gauge adapter is memoized: its digest is a SHA-256
    over the digest domain, the native datum identity, its four component
    identities, constant zero sentinels and offsets, and the two LIVE caloric
    datum signatures, so the key carries the domain, the datum identity, the
    four identities and both signatures read at this call.  Every other
    adapter type reads its own property.
    """

    if (
        not LEVERS_ENABLED
        or type(adapter) is not component_datum.NativeZeroGaugeComponentEnergyDatumAdapter
    ):
        return adapter.definition_digest  # type: ignore[attr-defined]
    key = (
        component_datum.NATIVE_ZERO_GAUGE_DATUM_ADAPTER_DIGEST_DOMAIN,
        adapter.energy_datum_id,
        adapter.component_ids,
        hexane.caloric_datum_signature(),
        water.caloric_datum_signature(),
    )
    return DATUM_DIGEST_MEMO.get_or_compute(key, lambda: adapter.definition_digest)


# -- one validation per issued payload ----------------------------------------

_ISSUE_STATISTICS = {"issues": 0, "refusal_rechecks": 0}


def encode_and_decode_validated(
    *,
    validate: Callable[[], None],
    encode: Callable[[], bytes],
    decode: Callable[[bytes], Any],
) -> Any:
    """Issue a payload with one complete-payload validation instead of two.

    The former issue validated the complete state, encoded it, and decoded
    the payload, whose decode validated the same values again with the same
    function.  Here the first validation runs only when the encode or the
    decode raises, BEFORE the exception propagates: a state the validation
    refuses therefore raises the validation's refusal (the former first
    refusal) whichever step noticed it, and a valid state that the encoder or
    the decoder refuses raises that step's own refusal, as before.  An
    accepted issue returns the decode of the same payload bytes as before.
    """

    if not LEVERS_ENABLED:
        validate()
        return decode(encode())
    _ISSUE_STATISTICS["issues"] += 1
    try:
        payload = encode()
    except Exception:
        _ISSUE_STATISTICS["refusal_rechecks"] += 1
        validate()
        raise
    try:
        return decode(payload)
    except Exception:
        _ISSUE_STATISTICS["refusal_rechecks"] += 1
        validate()
        raise


# -- the replayed bisection ---------------------------------------------------

_BISECTION_STATISTICS = {
    "replays": 0,
    "evaluations": 0,
    "decided_without_evaluation": 0,
    "fallbacks_bracket_not_certified": 0,
    "fallbacks_false_position_raised": 0,
    "fallbacks_replay_raised": 0,
    "fallbacks_not_monotone": 0,
    "largest_non_monotone_separation_k": 0.0,
}


def replayed_bisection(
    energy_at: Callable[[float], float],
    bisection: Callable[[Callable[[float], bool]], float],
    *,
    lo: float,
    hi: float,
    energy_lo: float,
    energy_hi: float,
    target: float,
) -> float:
    """The former caloric bisection's result with fewer energy evaluations.

    ``bisection(decide)`` is the adapter's own former bisection loop with the
    comparison ``energy_at(middle) < target`` abstracted as ``decide(middle)``;
    its result depends on those decisions only.  The replay first narrows a
    certified bracket ``[a, b]`` (``energy(a) < target <= energy(b)``, both
    evaluated) by Illinois false position to
    :data:`REPLAY_CERTIFIED_BRACKET_WIDTH_K`, then runs the same loop deciding
    a midpoint at or below ``a - REPLAY_DECISION_MARGIN_K`` as below and one at
    or above ``b + REPLAY_DECISION_MARGIN_K`` as above without evaluating it,
    and every other midpoint exactly as the former loop does.

    GUARD (fail closed): every energy the replay evaluated, the two bracket
    ends included, is checked non-decreasing in temperature on every replay;
    when it is not, or when the bracket ends are not certified, or when any
    evaluation raises, the former bisection runs from the start and its result
    (or its refusal) is returned.  The skipped decisions rest on the energy
    being monotone at separations of the margin and more; the check verifies
    monotonicity on every sampled point, not at the skipped midpoints.
    """

    if not LEVERS_ENABLED:
        return bisection(lambda middle: energy_at(middle) < target)
    statistics = _BISECTION_STATISTICS
    statistics["replays"] += 1

    def shipped() -> float:
        return bisection(lambda middle: energy_at(middle) < target)

    a, b = lo, hi
    fa, fb = energy_lo - target, energy_hi - target
    if not (fa < 0.0 <= fb):
        statistics["fallbacks_bracket_not_certified"] += 1
        return shipped()
    samples = [(lo, energy_lo), (hi, energy_hi)]
    side = 0
    try:
        for _ in range(100):
            if b - a <= REPLAY_CERTIFIED_BRACKET_WIDTH_K:
                break
            point = b - fb * (b - a) / (fb - fa)
            if not a < point < b:
                point = 0.5 * (a + b)
                if point in (a, b):
                    break
            statistics["evaluations"] += 1
            value = energy_at(point)
            samples.append((point, value))
            if value < target:
                a, fa = point, value - target
                if side == -1:
                    fb *= 0.5
                side = -1
            else:
                b, fb = point, value - target
                if side == 1:
                    fa *= 0.5
                side = 1
    except Exception:
        statistics["fallbacks_false_position_raised"] += 1
        return shipped()
    below_bound = a - REPLAY_DECISION_MARGIN_K
    above_bound = b + REPLAY_DECISION_MARGIN_K

    def decide(middle: float) -> bool:
        if middle <= below_bound:
            statistics["decided_without_evaluation"] += 1
            return True
        if middle >= above_bound:
            statistics["decided_without_evaluation"] += 1
            return False
        statistics["evaluations"] += 1
        value = energy_at(middle)
        samples.append((middle, value))
        return value < target

    try:
        result = bisection(decide)
    except Exception:
        statistics["fallbacks_replay_raised"] += 1
        return shipped()
    separation = _largest_non_monotone_separation_k(samples)
    if separation is not None:
        statistics["fallbacks_not_monotone"] += 1
        if separation > statistics["largest_non_monotone_separation_k"]:
            statistics["largest_non_monotone_separation_k"] = separation
        return shipped()
    return result


def _largest_non_monotone_separation_k(samples: list[tuple[float, float]]) -> float | None:
    """``None`` when the sampled energies never decrease with temperature.

    Otherwise the largest temperature separation of an adjacent sampled pair
    whose energy decreases (or differs at one temperature), a diagnostic of
    the scale at which the computed energy is not monotone.
    """

    ordered = sorted(samples, key=lambda item: item[0])
    largest = None
    for (t_left, e_left), (t_right, e_right) in zip(ordered, ordered[1:]):
        if (t_right == t_left and not e_right == e_left) or (
            t_right != t_left and not e_right >= e_left
        ):
            separation = t_right - t_left
            if largest is None or separation > largest:
                largest = separation
    return largest


# -- the boiling root's value-only residual and its sampled audit --------------
#
# Owner ruling of 2026-10-31 (record
# GT_PS2_DTDC_TOWER_BOILING_ROOT_VALUE_ONLY_RESIDUAL_BUILD_2026-10-31.md): the
# level-0 boiling root may evaluate its residual from the seed energy without
# building every candidate state; the accepted end state is still built and
# validated as before; a sampled typed audit of root equality is allowed.  The
# residual itself lives at its call site
# (``pd_k_dome_bypass_tray_integration``); here are its knob and its counters.
# WAS (chain 256): every interior iterate built through the saturation memos.

_BOILING_ROOT_STATISTICS = {
    "value_only_evaluations": 0,
    "value_only_roots": 0,
    "end_states_built_after_the_root": 0,
    "audited_roots": 0,
    "audited_ends": 0,
    "audit_disagreements": 0,
}


def boiling_root_audit_every_n() -> int | None:
    """The audit knob's value: ``None`` when absent, else a positive integer (refused otherwise)."""

    raw = BOILING_ROOT_AUDIT_EVERY_N_RAW
    if raw is None:
        return None
    if not LEVERS_ENABLED:
        raise BoilingRootAuditKnobRefused(
            f"{BOILING_ROOT_AUDIT_EVERY_N_ENVIRONMENT_VARIABLE} is set while "
            f"{TOWER_BITWISE_LEVERS_DISABLE_ENVIRONMENT_VARIABLE} runs the former boiling "
            "root, which has no value-only root to audit - refuse, do not ignore"
        )
    if not (raw.isascii() and raw.isdigit()) or raw.startswith("0"):
        raise BoilingRootAuditKnobRefused(
            f"{BOILING_ROOT_AUDIT_EVERY_N_ENVIRONMENT_VARIABLE} must be a positive decimal "
            f"integer without leading zeros, got {raw!r} - refuse, do not clamp"
        )
    return int(raw)


def value_only_root_armed() -> bool:
    """True when the boiling root runs its value-only residual (the levers' shipped path).

    Validates the audit knob first, so a malformed knob refuses at the first
    boiling root whichever path runs.
    """

    boiling_root_audit_every_n()
    return LEVERS_ENABLED


def note_value_only_evaluation() -> None:
    _BOILING_ROOT_STATISTICS["value_only_evaluations"] += 1


def note_end_state_built_after_the_root() -> None:
    _BOILING_ROOT_STATISTICS["end_states_built_after_the_root"] += 1


def boiling_root_audit_due() -> bool:
    """Count one completed value-only root; True when it is the N-th since the last audit."""

    every = boiling_root_audit_every_n()
    _BOILING_ROOT_STATISTICS["value_only_roots"] += 1
    return every is not None and _BOILING_ROOT_STATISTICS["value_only_roots"] % every == 0


def record_boiling_root_audit(where: str, ends: list[dict]) -> None:
    """Compare the value-only root's ends with the former full-candidate residual's.

    Each end is a dict with ``label``, ``removed_kg``, ``value_used`` (the
    residual value the root compared), ``value_former`` (the former residual
    at the same removal, memos bypassed), ``payload_used`` and
    ``payload_former`` (the encoded bytes of the accepted state and of the
    former candidate; ``None`` for an end whose state the root did not keep)
    and ``former_refusal`` (the former path's refusal there, or ``None``).
    Any refusal, any value differing in its binary64 bits, or any accepted
    state differing in a single encoded byte raises
    :class:`BoilingRootAuditDisagreement`.
    """

    _BOILING_ROOT_STATISTICS["audited_roots"] += 1
    _BOILING_ROOT_STATISTICS["audited_ends"] += len(ends)
    findings = []
    for end in ends:
        if end["former_refusal"] is not None:
            findings.append(
                f"{end['label']} end (removal {end['removed_kg']!r} kg): the former "
                f"full-candidate residual refuses: {end['former_refusal']}"
            )
            continue
        if float_bits(end["value_used"]) != float_bits(end["value_former"]):
            findings.append(
                f"{end['label']} end (removal {end['removed_kg']!r} kg): residual "
                f"{end['value_used']!r} (value-only) against {end['value_former']!r} (former)"
            )
        if end["payload_used"] is not None and end["payload_used"] != end["payload_former"]:
            findings.append(
                f"{end['label']} end (removal {end['removed_kg']!r} kg): the accepted "
                "state's encoded bytes differ from the former candidate's"
            )
    if findings:
        _BOILING_ROOT_STATISTICS["audit_disagreements"] += 1
        raise BoilingRootAuditDisagreement(
            f"{where}: the value-only boiling root disagrees with the former path: "
            + "; ".join(findings)
            + " - refuse, do not accept"
        )


def boiling_root_statistics() -> dict:
    return {
        **_BOILING_ROOT_STATISTICS,
        "audit_every_n_environment_variable": BOILING_ROOT_AUDIT_EVERY_N_ENVIRONMENT_VARIABLE,
        "audit_every_n_raw": BOILING_ROOT_AUDIT_EVERY_N_RAW,
    }


# -- the census the march stamps ----------------------------------------------


def lever_statistics() -> dict:
    """Every lever's counters, for the march record (diagnostic only)."""

    return {
        "enabled": LEVERS_ENABLED,
        "disable_environment_variable": TOWER_BITWISE_LEVERS_DISABLE_ENVIRONMENT_VARIABLE,
        "memos": {
            memo.name: memo.statistics()
            for memo in (
                SATURATION_SEED_ENERGY_MEMO,
                SATURATION_STATE_MEMO,
                SATURATION_SEED_ENERGY_VALUE_ONLY_MEMO,
                ENCODE_MEMO,
                DATUM_DIGEST_MEMO,
            )
        },
        "single_validation": dict(_ISSUE_STATISTICS),
        "replayed_bisection": dict(_BISECTION_STATISTICS),
        "boiling_root_value_only": boiling_root_statistics(),
    }


def clear_lever_state() -> None:
    """Empty every memo and zero every counter (tests only)."""

    for memo in (
        SATURATION_SEED_ENERGY_MEMO,
        SATURATION_STATE_MEMO,
        SATURATION_SEED_ENERGY_VALUE_ONLY_MEMO,
        ENCODE_MEMO,
        DATUM_DIGEST_MEMO,
    ):
        memo.clear()
    for table in (_ISSUE_STATISTICS, _BISECTION_STATISTICS, _BOILING_ROOT_STATISTICS):
        for name in table:
            table[name] = 0.0 if name.endswith("_k") else 0


__all__ = (
    "BOILING_ROOT_AUDIT_EVERY_N_ENVIRONMENT_VARIABLE",
    "BOILING_ROOT_AUDIT_EVERY_N_RAW",
    "DATUM_DIGEST_MEMO",
    "ENCODE_MEMO",
    "LEVERS_ENABLED",
    "REPLAY_CERTIFIED_BRACKET_WIDTH_K",
    "REPLAY_DECISION_MARGIN_K",
    "SATURATION_SEED_ENERGY_MEMO",
    "SATURATION_SEED_ENERGY_VALUE_ONLY_MEMO",
    "SATURATION_STATE_MEMO",
    "TOWER_BITWISE_LEVERS_DISABLE_ENVIRONMENT_VARIABLE",
    "BoilingRootAuditDisagreement",
    "BoilingRootAuditKnobRefused",
    "BoundedExactMemo",
    "boiling_root_audit_due",
    "boiling_root_audit_every_n",
    "boiling_root_statistics",
    "caloric_context",
    "clear_lever_state",
    "component_datum_definition_digest",
    "encode_and_decode_validated",
    "encode_state",
    "float_bits",
    "lever_statistics",
    "note_end_state_built_after_the_root",
    "note_value_only_evaluation",
    "record_boiling_root_audit",
    "replayed_bisection",
    "value_only_root_armed",
)
