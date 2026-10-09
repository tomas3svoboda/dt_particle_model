r"""Condensed-solvent-aware reduced-film boundary variant (BC-3, design D4).

AUTHORITY.  The ratified M2 design's D4 seam
(``docs/GT_PS2_RR2_BIRTH_CONTINUATION_DESIGN_2026-08-21.md``: the
reduced-film condensed-hexane variant, ``y_dew(T_s)`` per trial), as
RELOCATED by the BC-2 epoch measurement
(``docs/GT_PS2_BC2_EPOCH_MEASUREMENT_2026-08-22.md`` section 2): the epoch
clears the FLAKE surface decisively, but the reduced-film closure's own
sample-temperature sweep still evaluates the gas at the pinned BOUNDARY
composition across the film span whose cold end is the newborn dry piece,
and the open gas-only band refuses any sample at/below the hexane dew
locus.  Physically such a sample is not an error: the local gas there is
in local equilibrium with condensed solvent - the mobile-liquid branch
the ratified design names.

THE CONSTRUCTION.  A sweep sample strictly ABOVE the dew locus evaluates
byte-identically to the legacy closure.  A sample at/below the locus has
its composition become the saturation-locus value ``y_dew(T_sample)`` -
the exact frozen chart's own hexane-saturation composition at the sample
temperature, side-refined to the equilibrium gate's INCLUSIVE admission
(``a_h in (0, 1]``, Amendment 18) - and it validates against that
construction instead of the open gas-only band.  Nothing is clamped: the
construction is the chart's own phase boundary, every condensed sample
carries typed evidence into the caller-owned journal, and no tolerance
anywhere changes.

THE A1 COMMIT CHANNEL (BC-4a, ratified 2026-08-22).  Under A1
(``docs/GT_PS2_BC3_CONDENSED_SURFACE_BIRTH_2026-08-22.md``, Status block)
the typed mobile-liquid refusal on the ACCEPTED assembly becomes the
LAWFUL COMMIT PATH when and only when the boundary carries a
:class:`CondensedRootCommitChannel`: the accepted condensed root's typed
admission is RECORDED into the channel instead of refusing, so the birth
commit can book the condensate into the attached-hexane film channel
(the intensive twin of
``surface_active_set.ExternalSurfaceInventories.attached_hexane_kg_per_kg_dry``,
built beside in ``cut_birth_condensed_film``).  The channel is default-off
(``None``): a condensed root WITHOUT the channel refuses byte-identically
to the pre-A1 gate, and every non-birth caller keeps exactly today's
behavior.

TYPE SURFACE.  :class:`CondensedSolventReducedFilmPoreBoundary` is a
frozen dataclass subclass of the frozen
:class:`coupled_transport.ReducedFilmPoreBoundary`, so every existing
``isinstance`` admission (the birth assembly's boundary gate, the birth
integrator's boundary gate, ``prepare_reduced_film_surface_state``) admits
it unchanged.  Its equality is honest - a variant never compares equal to
a plain boundary - so the radial transaction host's boundary-identity
gate consults :func:`is_condensed_solvent_promotion_of` explicitly.  The
one behavioral dispatch lives at the D4 seam itself
(``cut_transport._dry_face_fluxes``), default-off: a plain boundary takes
exactly the legacy path.

Every claim flag reachable from this module is False.  Nothing here makes
the particle model physically qualifying, plant predictive, or production
wired.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, fields

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct

#: Hard cap on admissible-side ulp refinement steps for the locus root.
DEW_LOCUS_SIDE_REFINEMENT_LIMIT = 64

#: Default cap on RETAINED detailed condensed-sample records per journal.
#: Counters and extrema keep accumulating past the cap; only the detailed
#: per-sample records stop being appended.  Record-only: nothing numeric
#: reads the journal.
DEFAULT_JOURNAL_RECORD_CAPACITY = 512


class CondensedSolventTopologyError(ValueError):
    """The condensed-solvent construction is unlawful here; fail closed."""


@dataclass(frozen=True, slots=True)
class CondensedSolventLocusConstruction:
    """The exact chart's hexane dew-locus composition at one ``(T, P)``.

    ``y_dew`` is the frozen gas-only chart's own hexane-saturation upper
    composition bound at this temperature and pressure, refined
    monotonically downward by at most
    :data:`DEW_LOCUS_SIDE_REFINEMENT_LIMIT` ulps until the equilibrium
    gate's inclusive admission ``a_h <= 1`` holds.  The construction is
    exact-form, never a pinned float: re-derive it from the chart, do not
    copy the number.
    """

    temperature_k: float
    pressure_pa: float
    y_dew: float
    hexane_activity: float
    water_activity: float
    side_refinement_steps: int
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.temperature_k,
            self.pressure_pa,
            self.y_dew,
            self.hexane_activity,
            self.water_activity,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("dew-locus construction values must be finite")
        if not 0.0 < self.y_dew < 1.0:
            raise ValueError("dew-locus composition must lie strictly inside (0,1)")
        if not 0.0 < self.hexane_activity <= 1.0:
            raise ValueError(
                "dew-locus construction must satisfy the inclusive equilibrium "
                "admission a_h in (0,1]"
            )
        if (
            isinstance(self.side_refinement_steps, bool)
            or type(self.side_refinement_steps) is not int
            or self.side_refinement_steps < 0
        ):
            raise ValueError("side refinement steps must be a nonnegative exact int")


def dew_locus_composition(
    temperature_k: float,
    pressure_pa: float,
    pore: cp.CoupledPoreParams,
) -> CondensedSolventLocusConstruction:
    """Construct ``y_dew(T)`` from the frozen chart's own saturation bound.

    Fail-closed refusals:

    * pure hexane vapor is not saturated at this ``(T, P)`` (the locus does
      not bind - the measured ``T_max`` bracketing limit);
    * the Luikov lower authority binds before hexane saturation (the
      condensed construction would leave the lawful water-side band);
    * the admissible-side refinement does not reach ``a_h <= 1`` inside the
      declared ulp budget.
    """

    if type(pore) is not cp.CoupledPoreParams:
        raise TypeError("dew-locus construction requires exact CoupledPoreParams")
    _, hexane_at_one = cp.binary_gas_component_activities(
        temperature_k, pressure_pa, 1.0, pore
    )
    if not hexane_at_one > 1.0:
        raise CondensedSolventTopologyError(
            "the hexane saturation locus does not bind at this state: pure "
            f"hexane vapor has a_h(y=1) = {hexane_at_one!r} <= 1 at "
            f"T = {temperature_k!r} K, P = {pressure_pa!r} Pa; there is no "
            "dew composition to construct"
        )
    interval = cp.gas_only_composition_interval(temperature_k, pressure_pa, pore)
    if interval.upper_constraint is not cp.GasOnlyCompositionConstraint.HEXANE_SATURATION:
        raise CondensedSolventTopologyError(
            "the gas-only upper composition bound at "
            f"T = {temperature_k!r} K is {interval.upper_constraint.value!r}, "
            "not hexane saturation: the condensed-solvent construction would "
            "leave the lawful water-side band; refuse, do not clamp"
        )
    composition = interval.upper_y_hexane
    refinement_steps = 0
    water_activity = math.nan
    hexane_activity = math.nan
    for _ in range(DEW_LOCUS_SIDE_REFINEMENT_LIMIT + 1):
        water_activity, hexane_activity = cp.binary_gas_component_activities(
            temperature_k, pressure_pa, composition, pore
        )
        if hexane_activity <= 1.0:
            break
        composition = math.nextafter(composition, 0.0)
        refinement_steps += 1
    else:
        raise CondensedSolventTopologyError(
            "dew-locus admissible-side refinement did not reach a_h <= 1 "
            f"inside {DEW_LOCUS_SIDE_REFINEMENT_LIMIT} ulps at "
            f"T = {temperature_k!r} K (last a_h = {hexane_activity!r})"
        )
    if not composition > interval.lower_y_hexane:
        raise CondensedSolventTopologyError(
            "dew-locus composition collapsed onto the lower gas-only bound; "
            "the condensed construction has no lawful interior representative"
        )
    return CondensedSolventLocusConstruction(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        y_dew=composition,
        hexane_activity=hexane_activity,
        water_activity=water_activity,
        side_refinement_steps=refinement_steps,
    )


@dataclass(frozen=True, slots=True)
class CondensedSolventSampleEvidence:
    """Typed evidence for ONE sweep sample that took the condensed branch.

    The frozen seam evaluates each sweep sample at TWO temperatures: the
    sample temperature itself (the open-band validation) and the symmetric
    face midpoint ``0.5*(T_cell + T_sample)`` (the face-state equilibrium
    inside the flux).  The construction is therefore taken at
    ``evaluation_temperature_k = min`` of the two - the coldest temperature
    the seam evaluates this sample at - which at the sweep's binding cold
    end (``T_sample == T_cell``) is EXACTLY the ratified D4 reading
    ``y_dew(T_sample)``.
    """

    temperature_k: float
    face_temperature_k: float
    evaluation_temperature_k: float
    pressure_pa: float
    boundary_y_hexane: float
    sample_hexane_activity: float
    face_hexane_activity: float
    trigger_hexane_activity: float
    constructed_y_hexane: float
    constructed_hexane_activity: float
    constructed_water_activity: float
    side_refinement_steps: int
    physically_qualifying: bool = field(default=False, init=False)

    def journal_record(self) -> dict:
        return {
            "temperature_k": self.temperature_k,
            "face_temperature_k": self.face_temperature_k,
            "evaluation_temperature_k": self.evaluation_temperature_k,
            "pressure_pa": self.pressure_pa,
            "boundary_y_hexane": self.boundary_y_hexane,
            "sample_hexane_activity": self.sample_hexane_activity,
            "face_hexane_activity": self.face_hexane_activity,
            "trigger_hexane_activity": self.trigger_hexane_activity,
            "constructed_y_hexane": self.constructed_y_hexane,
            "constructed_hexane_activity": self.constructed_hexane_activity,
            "constructed_water_activity": self.constructed_water_activity,
            "side_refinement_steps": self.side_refinement_steps,
            "physically_qualifying": False,
        }


class CondensedSolventSweepJournal:
    """Append-only, record-only evidence carrier for one closure consumer.

    The journal NEVER influences a numeric result: it is excluded from the
    boundary variant's equality and repr, so replay determinism and the
    radial transaction digests are untouched.  Counters and extrema always
    accumulate; detailed per-sample records stop at ``record_capacity``.
    """

    __slots__ = (
        "record_capacity",
        "_records",
        "_above_locus_samples",
        "_condensed_samples",
        "_condensed_root_evaluations",
        "_coldest_condensed_root_temperature_k",
        "_elided_condensed_records",
        "_coldest_condensed_temperature_k",
        "_warmest_condensed_temperature_k",
        "_minimum_constructed_y_hexane",
        "_maximum_trigger_hexane_activity",
    )

    def __init__(self, record_capacity: int = DEFAULT_JOURNAL_RECORD_CAPACITY) -> None:
        if (
            isinstance(record_capacity, bool)
            or type(record_capacity) is not int
            or record_capacity < 1
        ):
            raise ValueError("journal record capacity must be a positive exact int")
        self.record_capacity = record_capacity
        self._records: list[CondensedSolventSampleEvidence] = []
        self._above_locus_samples = 0
        self._condensed_samples = 0
        self._condensed_root_evaluations = 0
        self._coldest_condensed_root_temperature_k = math.inf
        self._elided_condensed_records = 0
        self._coldest_condensed_temperature_k = math.inf
        self._warmest_condensed_temperature_k = -math.inf
        self._minimum_constructed_y_hexane = math.inf
        self._maximum_trigger_hexane_activity = -math.inf

    def record_above_locus(self) -> None:
        self._above_locus_samples += 1

    def record_condensed_root(self, evaluation_temperature_k: float) -> None:
        self._condensed_root_evaluations += 1
        self._coldest_condensed_root_temperature_k = min(
            self._coldest_condensed_root_temperature_k, evaluation_temperature_k
        )

    def record_condensed(self, evidence: CondensedSolventSampleEvidence) -> None:
        if type(evidence) is not CondensedSolventSampleEvidence:
            raise TypeError("journal accepts exact CondensedSolventSampleEvidence")
        self._condensed_samples += 1
        # The extrema track the CONSTRUCTION temperature (the coldest
        # temperature the seam evaluated the sample at), which is what the
        # span-below-the-locus measurement reads.
        self._coldest_condensed_temperature_k = min(
            self._coldest_condensed_temperature_k, evidence.evaluation_temperature_k
        )
        self._warmest_condensed_temperature_k = max(
            self._warmest_condensed_temperature_k, evidence.evaluation_temperature_k
        )
        self._minimum_constructed_y_hexane = min(
            self._minimum_constructed_y_hexane, evidence.constructed_y_hexane
        )
        self._maximum_trigger_hexane_activity = max(
            self._maximum_trigger_hexane_activity, evidence.trigger_hexane_activity
        )
        if len(self._records) < self.record_capacity:
            self._records.append(evidence)
        else:
            self._elided_condensed_records += 1

    @property
    def above_locus_sample_count(self) -> int:
        return self._above_locus_samples

    @property
    def condensed_sample_count(self) -> int:
        return self._condensed_samples

    @property
    def total_sample_evaluations(self) -> int:
        return self._above_locus_samples + self._condensed_samples

    @property
    def condensed_root_evaluations(self) -> int:
        return self._condensed_root_evaluations

    @property
    def coldest_condensed_root_temperature_k(self) -> float:
        return self._coldest_condensed_root_temperature_k

    @property
    def condensed_records(self) -> tuple[CondensedSolventSampleEvidence, ...]:
        return tuple(self._records)

    @property
    def elided_condensed_records(self) -> int:
        return self._elided_condensed_records

    @property
    def coldest_condensed_temperature_k(self) -> float:
        return self._coldest_condensed_temperature_k

    @property
    def warmest_condensed_temperature_k(self) -> float:
        return self._warmest_condensed_temperature_k

    @property
    def minimum_constructed_y_hexane(self) -> float:
        return self._minimum_constructed_y_hexane

    @property
    def maximum_trigger_hexane_activity(self) -> float:
        return self._maximum_trigger_hexane_activity

    def journal_record(self) -> dict:
        return {
            "journal": "condensed-solvent-reduced-film-sweep",
            "above_locus_sample_count": self._above_locus_samples,
            "condensed_sample_count": self._condensed_samples,
            "condensed_root_evaluations": self._condensed_root_evaluations,
            "coldest_condensed_root_temperature_k": (
                None
                if self._condensed_root_evaluations == 0
                else self._coldest_condensed_root_temperature_k
            ),
            "elided_condensed_records": self._elided_condensed_records,
            "coldest_condensed_temperature_k": (
                None
                if self._condensed_samples == 0
                else self._coldest_condensed_temperature_k
            ),
            "warmest_condensed_temperature_k": (
                None
                if self._condensed_samples == 0
                else self._warmest_condensed_temperature_k
            ),
            "minimum_constructed_y_hexane": (
                None
                if self._condensed_samples == 0
                else self._minimum_constructed_y_hexane
            ),
            "maximum_trigger_hexane_activity": (
                None
                if self._condensed_samples == 0
                else self._maximum_trigger_hexane_activity
            ),
            "condensed_records": [
                record.journal_record() for record in self._records
            ],
            "physically_qualifying": False,
        }


@dataclass(frozen=True, slots=True)
class CondensedRootCommitAdmission:
    """Typed evidence that ONE accepted surface root took the condensed branch.

    Recorded by :func:`require_root_above_locus` on the ACCEPTED
    (threshold-enforced) assembly when and only when the boundary carries a
    :class:`CondensedRootCommitChannel` (A1's lawful commit path).  The
    construction values are the same exact-chart dew-locus construction the
    trial composition path uses (:func:`condensed_root_equilibrium_composition`
    with identical arguments), so the admission and the committed surface
    equilibrium are one construction, never two estimates.
    """

    surface_temperature_k: float
    face_temperature_k: float
    evaluation_temperature_k: float
    cell_temperature_k: float
    pressure_pa: float
    boundary_y_hexane: float
    sample_hexane_activity: float
    face_hexane_activity: float
    constructed_y_hexane: float
    constructed_hexane_activity: float
    constructed_water_activity: float
    side_refinement_steps: int
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.surface_temperature_k,
            self.face_temperature_k,
            self.evaluation_temperature_k,
            self.cell_temperature_k,
            self.pressure_pa,
            self.boundary_y_hexane,
            self.sample_hexane_activity,
            self.face_hexane_activity,
            self.constructed_y_hexane,
            self.constructed_hexane_activity,
            self.constructed_water_activity,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("condensed-root admission values must be finite")
        if not 0.0 < self.constructed_y_hexane < 1.0:
            raise ValueError(
                "condensed-root admission composition must lie strictly inside (0,1)"
            )
        if not 0.0 < self.constructed_hexane_activity <= 1.0:
            raise ValueError(
                "condensed-root admission must satisfy the inclusive equilibrium "
                "admission a_h in (0,1]"
            )
        if not max(self.sample_hexane_activity, self.face_hexane_activity) >= 1.0:
            raise ValueError(
                "a condensed-root admission requires the root's own condensed "
                "classification (a_h >= 1 at the sample or its face midpoint)"
            )
        if (
            isinstance(self.side_refinement_steps, bool)
            or type(self.side_refinement_steps) is not int
            or self.side_refinement_steps < 0
        ):
            raise ValueError("side refinement steps must be a nonnegative exact int")

    def journal_record(self) -> dict:
        return {
            "admission": "condensed-root-commit",
            "surface_temperature_k": self.surface_temperature_k,
            "face_temperature_k": self.face_temperature_k,
            "evaluation_temperature_k": self.evaluation_temperature_k,
            "cell_temperature_k": self.cell_temperature_k,
            "pressure_pa": self.pressure_pa,
            "boundary_y_hexane": self.boundary_y_hexane,
            "sample_hexane_activity": self.sample_hexane_activity,
            "face_hexane_activity": self.face_hexane_activity,
            "constructed_y_hexane": self.constructed_y_hexane,
            "constructed_hexane_activity": self.constructed_hexane_activity,
            "constructed_water_activity": self.constructed_water_activity,
            "side_refinement_steps": self.side_refinement_steps,
            "physically_qualifying": False,
        }


class CondensedRootCommitChannel:
    """Caller-owned single-admission carrier for A1's lawful condensed commit.

    One channel books exactly ONE accepted condensed root: a second recording
    is a typed refusal, which forces the fresh-channel-per-solve discipline
    (exactly the sweep journal's own per-evaluation discipline).  Like the
    journal, the channel never influences a numeric result and is excluded
    from the boundary variant's equality and repr, so replay determinism and
    the radial transaction digests are untouched.
    """

    __slots__ = ("_admission",)

    def __init__(self) -> None:
        self._admission: CondensedRootCommitAdmission | None = None

    @property
    def admission(self) -> CondensedRootCommitAdmission | None:
        return self._admission

    @property
    def admitted(self) -> bool:
        return self._admission is not None

    def record_admission(self, admission: CondensedRootCommitAdmission) -> None:
        if type(admission) is not CondensedRootCommitAdmission:
            raise TypeError(
                "the commit channel accepts exact CondensedRootCommitAdmission"
            )
        if self._admission is not None:
            raise CondensedSolventTopologyError(
                "the condensed-root commit channel already carries an admission: "
                "one channel books exactly one accepted root; construct a fresh "
                "channel per solve, never reuse one across solves"
            )
        self._admission = admission

    def journal_record(self) -> dict:
        return {
            "channel": "condensed-root-commit",
            "admitted": self._admission is not None,
            "admission": (
                None if self._admission is None else self._admission.journal_record()
            ),
            "physically_qualifying": False,
        }


@dataclass(frozen=True)
class CondensedSolventReducedFilmPoreBoundary(ct.ReducedFilmPoreBoundary):
    """D4's condensed-solvent-aware reduced-film boundary variant.

    Identical bulk state and film pair as the plain boundary it promotes;
    the ONLY difference is the lawful evaluation of sweep samples at/below
    the hexane dew locus, dispatched at the D4 seam in
    ``cut_transport._dry_face_fluxes``.  The journal is caller-owned,
    record-only, and excluded from equality/repr so replay determinism is
    untouched.  Equality stays honest: a variant never equals a plain
    boundary - the transaction host consults
    :func:`is_condensed_solvent_promotion_of` for the lawful promotion.
    """

    sweep_journal: CondensedSolventSweepJournal = field(
        compare=False, repr=False, kw_only=True
    )
    #: BC-4a (A1): default-off commit channel.  ``None`` keeps the pre-A1
    #: refusal byte-identical; a channel makes the accepted condensed root
    #: the lawful commit path.  Excluded from equality and repr exactly like
    #: the journal, so promotion identity and replay digests are untouched.
    condensed_root_commit_channel: CondensedRootCommitChannel | None = field(
        default=None, compare=False, repr=False, kw_only=True
    )

    def __post_init__(self) -> None:
        super().__post_init__()
        if type(self.sweep_journal) is not CondensedSolventSweepJournal:
            raise TypeError(
                "the condensed-solvent boundary requires an exact "
                "CondensedSolventSweepJournal"
            )
        if self.condensed_root_commit_channel is not None and (
            type(self.condensed_root_commit_channel) is not CondensedRootCommitChannel
        ):
            raise TypeError(
                "the condensed-solvent boundary's commit channel must be an "
                "exact CondensedRootCommitChannel or None"
            )


def promote_reduced_film_boundary(
    boundary: ct.ReducedFilmPoreBoundary,
    *,
    sweep_journal: CondensedSolventSweepJournal,
    condensed_root_commit_channel: CondensedRootCommitChannel | None = None,
) -> CondensedSolventReducedFilmPoreBoundary:
    """Promote one exact plain reduced-film boundary to the condensed variant.

    Every base field is copied unchanged (the promotion IS the same film
    boundary); only the condensed sweep evaluation and its evidence journal
    are added.  Exact-type gated: a Dirichlet oracle has no condensed
    variant, and a variant is never re-promoted.  The optional
    ``condensed_root_commit_channel`` (BC-4a, A1) arms the lawful commit path
    for the accepted condensed root; omitted, the promotion behaves
    byte-identically to the pre-A1 variant.
    """

    if type(boundary) is not ct.ReducedFilmPoreBoundary:
        raise TypeError(
            "promotion requires an exact plain ReducedFilmPoreBoundary; got "
            f"{type(boundary).__name__}"
        )
    base_values = {
        base_field.name: getattr(boundary, base_field.name)
        for base_field in fields(ct.ReducedFilmPoreBoundary)
    }
    return CondensedSolventReducedFilmPoreBoundary(
        **base_values,
        sweep_journal=sweep_journal,
        condensed_root_commit_channel=condensed_root_commit_channel,
    )


def is_condensed_solvent_promotion_of(actual: object, expected: object) -> bool:
    """Return whether ``actual`` is the lawful condensed promotion of ``expected``.

    True exactly when ``expected`` is an exact plain reduced-film boundary,
    ``actual`` is the exact condensed variant, and every base field is
    equal.  This is the ONE comparison the radial transaction host uses to
    admit a promoted boundary through its identity gate; everything else
    keeps honest class-exact dataclass equality.
    """

    if type(expected) is not ct.ReducedFilmPoreBoundary:
        return False
    if type(actual) is not CondensedSolventReducedFilmPoreBoundary:
        return False
    return all(
        getattr(actual, base_field.name) == getattr(expected, base_field.name)
        for base_field in fields(ct.ReducedFilmPoreBoundary)
    )


def _sample_locus_classification(
    boundary: CondensedSolventReducedFilmPoreBoundary,
    surface_temperature_k: float,
    cell_temperature_k: float,
    pressure_pa: float,
    pore: cp.CoupledPoreParams,
) -> tuple[bool, float, float, float]:
    """Classify one sweep sample against the seam's own two evaluations.

    The frozen seam evaluates a surface sweep sample at the boundary
    composition TWICE: the open-band validation at the sample temperature
    (strict, ``a_h < 1``) and the face-state equilibrium at the exact
    symmetric midpoint ``0.5*(T_cell + T_sample)`` (inclusive,
    ``a_h <= 1``).  The condensed branch engages exactly when either
    admission would refuse on the hexane-saturation side.  Returns
    ``(condensed, face_temperature, a_h_sample, a_h_face)``.
    """

    # EXACTLY _dry_flux's symmetric endpoint expression, same operand order.
    face_temperature_k = 0.5 * (cell_temperature_k + surface_temperature_k)
    _, sample_hexane_activity = cp.binary_gas_component_activities(
        surface_temperature_k, pressure_pa, boundary.y_hexane, pore
    )
    _, face_hexane_activity = cp.binary_gas_component_activities(
        face_temperature_k, pressure_pa, boundary.y_hexane, pore
    )
    condensed = sample_hexane_activity >= 1.0 or face_hexane_activity > 1.0
    return condensed, face_temperature_k, sample_hexane_activity, face_hexane_activity


def resolve_condensed_surface_sample(
    boundary: CondensedSolventReducedFilmPoreBoundary,
    *,
    surface_temperature_k: float,
    cell_temperature_k: float,
    pressure_pa: float,
    pore: cp.CoupledPoreParams,
) -> CondensedSolventSampleEvidence | None:
    """Classify one sweep sample against the dew locus at the boundary ``y``.

    Returns ``None`` when both of the seam's own evaluations of this sample
    at the boundary composition are strictly above the locus (the caller
    then evaluates the EXACT legacy path), or the typed condensed-sample
    evidence with the constructed ``y_dew`` taken at the COLDEST temperature
    the seam evaluates this sample at, ``min(T_sample, T_face)`` - at the
    sweep's binding cold end (``T_sample == T_cell``, where the BC-2
    measurement pinned the refusal) that is exactly the ratified D4 reading
    ``y_dew(T_sample)``.  A condensed sample is certified against the frozen
    equilibrium gate's own inclusive admission at that same temperature
    before it is returned; certification failing is a typed fail-closed
    refusal, never a clamp.
    """

    if type(boundary) is not CondensedSolventReducedFilmPoreBoundary:
        raise TypeError(
            "condensed sample resolution requires the exact condensed-solvent "
            "boundary variant"
        )
    (
        condensed,
        face_temperature_k,
        sample_hexane_activity,
        face_hexane_activity,
    ) = _sample_locus_classification(
        boundary, surface_temperature_k, cell_temperature_k, pressure_pa, pore
    )
    if not condensed:
        boundary.sweep_journal.record_above_locus()
        return None
    evaluation_temperature_k = min(surface_temperature_k, face_temperature_k)
    trigger_hexane_activity = max(sample_hexane_activity, face_hexane_activity)
    construction = dew_locus_composition(evaluation_temperature_k, pressure_pa, pore)
    if not construction.y_dew <= boundary.y_hexane:
        raise CondensedSolventTopologyError(
            "dew-locus construction lost activity/composition monotonicity: "
            f"y_dew = {construction.y_dew!r} above the boundary composition "
            f"{boundary.y_hexane!r} while a_h = {trigger_hexane_activity!r} "
            "triggered the condensed branch"
        )
    try:
        cp.evaluate_equilibrium(
            evaluation_temperature_k,
            pressure_pa,
            construction.y_dew,
            pore,
        )
    except cp.CoupledPoreTopologyError as exc:
        raise CondensedSolventTopologyError(
            "the condensed-solvent construction is not admitted by the frozen "
            f"equilibrium gate at T = {evaluation_temperature_k!r} K, "
            f"y_dew = {construction.y_dew!r}: {exc}"
        ) from exc
    evidence = CondensedSolventSampleEvidence(
        temperature_k=surface_temperature_k,
        face_temperature_k=face_temperature_k,
        evaluation_temperature_k=evaluation_temperature_k,
        pressure_pa=pressure_pa,
        boundary_y_hexane=boundary.y_hexane,
        sample_hexane_activity=sample_hexane_activity,
        face_hexane_activity=face_hexane_activity,
        trigger_hexane_activity=trigger_hexane_activity,
        constructed_y_hexane=construction.y_dew,
        constructed_hexane_activity=construction.hexane_activity,
        constructed_water_activity=construction.water_activity,
        side_refinement_steps=construction.side_refinement_steps,
    )
    boundary.sweep_journal.record_condensed(evidence)
    return evidence


def condensed_root_equilibrium_composition(
    boundary: CondensedSolventReducedFilmPoreBoundary,
    *,
    surface_temperature_k: float,
    cell_temperature_k: float,
    pressure_pa: float,
    pore: cp.CoupledPoreParams,
) -> float:
    """Return the lawful surface-equilibrium composition for a converged root.

    An above-locus root returns the boundary composition UNCHANGED (the
    legacy value, bit for bit).  A condensed root - lawful only on
    un-thresholded TRIAL assemblies, per D4's ``y_dew(T_s)`` per-trial
    construction - returns the constructed ``y_dew`` at the coldest
    temperature the seam evaluates the root at, and records the condensed
    root in the journal.  The ACCEPTED (threshold-enforced) assembly never
    reaches this branch: :func:`require_root_above_locus` refuses first.
    """

    if type(boundary) is not CondensedSolventReducedFilmPoreBoundary:
        raise TypeError(
            "root composition resolution requires the exact condensed-solvent "
            "boundary variant"
        )
    condensed, face_temperature_k, _, _ = _sample_locus_classification(
        boundary, surface_temperature_k, cell_temperature_k, pressure_pa, pore
    )
    if not condensed:
        return boundary.y_hexane
    evaluation_temperature_k = min(surface_temperature_k, face_temperature_k)
    construction = dew_locus_composition(evaluation_temperature_k, pressure_pa, pore)
    boundary.sweep_journal.record_condensed_root(evaluation_temperature_k)
    return construction.y_dew


def require_root_above_locus(
    boundary: CondensedSolventReducedFilmPoreBoundary,
    *,
    surface_temperature_k: float,
    cell_temperature_k: float,
    pressure_pa: float,
    pore: cp.CoupledPoreParams,
) -> None:
    """Gate a converged surface ROOT on the condensed branch: refuse or admit.

    The condensed construction admits sub-locus sweep SAMPLES so the scalar
    closure can bracket lawfully; a converged surface root whose own
    evaluation took the condensed branch means the surface itself is on the
    mobile-liquid branch.  Without a commit channel that is outside this
    variant's admission and refuses typed, byte-identically to the pre-A1
    gate - same predicate as the sample branch, fail closed, never a clamp.

    WITH a :class:`CondensedRootCommitChannel` on the boundary (BC-4a, under
    A1 RATIFIED 2026-08-22), the refusal IS the lawful commit path: the
    accepted condensed root's typed admission - the same exact-chart
    ``y_dew`` construction the trial composition path uses - is recorded
    into the channel so the birth commit can book the condensate into the
    attached-hexane film channel.  The channel books exactly one root;
    everything else about the gate is unchanged.
    """

    if type(boundary) is not CondensedSolventReducedFilmPoreBoundary:
        raise TypeError(
            "root admission requires the exact condensed-solvent boundary variant"
        )
    (
        condensed,
        face_temperature_k,
        sample_hexane_activity,
        face_hexane_activity,
    ) = _sample_locus_classification(
        boundary, surface_temperature_k, cell_temperature_k, pressure_pa, pore
    )
    if not condensed:
        return
    channel = boundary.condensed_root_commit_channel
    if channel is None:
        raise CondensedSolventTopologyError(
            "the reduced-film surface root landed at/below the hexane dew "
            f"locus: a_h(T_s = {surface_temperature_k!r} K, y_b = "
            f"{boundary.y_hexane!r}) = {sample_hexane_activity!r}, a_h at the "
            f"face midpoint {face_temperature_k!r} K = "
            f"{face_hexane_activity!r}.  The surface is on the mobile-liquid "
            "branch; the condensed-solvent variant admits sub-locus sweep "
            "SAMPLES only.  Activate the mobile-liquid surface topology, do "
            "not clamp"
        )
    # A1's lawful commit: the SAME construction the trial composition path
    # evaluates (condensed_root_equilibrium_composition with identical
    # arguments), recorded as the typed admission the film channel books.
    evaluation_temperature_k = min(surface_temperature_k, face_temperature_k)
    construction = dew_locus_composition(evaluation_temperature_k, pressure_pa, pore)
    channel.record_admission(
        CondensedRootCommitAdmission(
            surface_temperature_k=surface_temperature_k,
            face_temperature_k=face_temperature_k,
            evaluation_temperature_k=evaluation_temperature_k,
            cell_temperature_k=cell_temperature_k,
            pressure_pa=pressure_pa,
            boundary_y_hexane=boundary.y_hexane,
            sample_hexane_activity=sample_hexane_activity,
            face_hexane_activity=face_hexane_activity,
            constructed_y_hexane=construction.y_dew,
            constructed_hexane_activity=construction.hexane_activity,
            constructed_water_activity=construction.water_activity,
            side_refinement_steps=construction.side_refinement_steps,
        )
    )


__all__ = [
    "DEW_LOCUS_SIDE_REFINEMENT_LIMIT",
    "DEFAULT_JOURNAL_RECORD_CAPACITY",
    "CondensedRootCommitAdmission",
    "CondensedRootCommitChannel",
    "CondensedSolventLocusConstruction",
    "CondensedSolventReducedFilmPoreBoundary",
    "CondensedSolventSampleEvidence",
    "CondensedSolventSweepJournal",
    "CondensedSolventTopologyError",
    "condensed_root_equilibrium_composition",
    "dew_locus_composition",
    "is_condensed_solvent_promotion_of",
    "promote_reduced_film_boundary",
    "require_root_above_locus",
    "resolve_condensed_surface_sample",
]
