"""PART-02 opener: the cell-zero extinction-localizer composition.

AUTHORITY.  ``docs/GT_PS2_PART01_RECONVENED_CEREMONY_PACKET_2026-08-22.md``
(RULED): C-5 routes "the cell-zero extinction-localizer composition (the
arrival terminus's own next step)" into PART-02, and the roadmap's PART-02
criterion (``docs/GT_PS2_TRAY_PARTICLE_PROCESS_MODEL_QUALIFICATION_ROADMAP_
2026-08-15.md``) demands that the canonical trajectory reach "an OUTCOME
rather than a software refusal" through "activation, wet lift, birth, every
face, extinction and dry continuation".  BC-7's certified N=2 arrival
(``docs/GT_PS2_BC7_FACE_EVENT_MARCH_AND_CONTINUED_TRAVERSE_2026-08-22.md``
section 4.1) committed into cell zero and halted TYPED on the very next step
with ``CutEventMacrostepError`` code
``cell_zero_extinction_localizer_required`` - the honest window terminus.
This module closes that terminus BESIDE the march entry, exactly as BC-7's
``qsc_hexane_face_event_march`` closed the interior-face seam: nothing in the
certified extinction chain is edited or re-implemented.

THE COMPOSITION.  Cell zero has no interior inner face; its terminus belongs
to the separately certified extinction family, and every stage of that family
already exists:

  * the decreasing-target ladder with monotone times comes from
    ``cut_extinction_time.solve_positive_core_target`` (each rung solves the
    UNCHANGED positive-core ALE/RH equations at a fixed ``z_* > 0`` with the
    step duration as the exchanged unknown) chained by
    ``cut_extinction_time.continuation_seed``;
  * the finite-extinction certificate is
    ``cut_extinction_time.certify_finite_extinction_time`` - decreasing
    targets, strictly increasing times, contracting tail increments, a
    sufficiently small last target, and a geometric tail bound inside the
    declared tolerance, or NO certificate;
  * the exact ``z=0`` projection and the post-extinction fully-dry remainder
    are committed ATOMICALLY by
    ``cut_extinction_orchestrator.advance_numerical_extinction_macrostep``,
    which internally routes ``cut_extinction_time.
    project_certified_extinction`` and - when the requested duration extends
    beyond the certified extinction time - the existing fully-dry surface
    ``coupled_transport.advance_fully_dry_backward_euler`` on the
    label-preserving ``to_uniform_fully_dry_transport_state`` handoff.

What this layer adds is exactly the routing the march never had: the ladder
SCHEDULE (declared, typed, geometric), the first rung's seed built from the
last ACCEPTED step's own primitives and committed front rate (the identical
accepted-history discipline BC-7's bracket used - no rejected candidate, no
clamp, no manufactured event-time authority), the projection seed mapped
structurally from the deepest accepted rung, and the typed outcome payload.

THE DRY CONTINUATION.  The certified surface for a fully-dry particle's
remaining evolution is ``coupled_transport.advance_fully_dry_backward_euler``
reached ONLY through the extinction orchestrator's
``EXTINCTION_AND_FULLY_DRY_REMAINDER`` branch, guarded by the
label-preserving handoff (a nonuniform residual-oil endpoint refuses typed
with ``LABEL_PRESERVING_HANDOFF_REQUIRED`` rather than homogenising).  This
module never calls the dry surface directly and adds no law to it.

THE FILM COMPANION DOES NOT ARISE.  The BC-7 crossing into cell zero was
measured film-free (minimum surface margin above T* = +0.999 K; the A5 latch
cleared long before the clean march).  The extinction family defines no
film-carrying chart, so this composition ASSERTS the film-free entry typed:
a caller holding a film companion is refused with
``film_companion_present_at_extinction`` - never silently dropped.

Typed everything; exact rollback: no entry point mutates an input, and every
refusal from the certified layers (``PositiveCoreTargetStepError``,
``ExtinctionTimeLimitError``, ``ExtinctionMacrostepError``) propagates
untouched with its own rollback identity.  No tolerance is widened anywhere:
the finite-limit refinement contract is ``cut_extinction_time.
ExtinctionTimeRefinementControls``'s own defaults unless the caller declares
a stricter one, and the residual/ledger/condition contracts are the state's
own.  Every claim flag reachable from this module is False; nothing here
makes the particle model physically qualifying, plant predictive, or
production wired.
"""

from __future__ import annotations

import math
from dataclasses import InitVar, dataclass, field
from enum import Enum
from typing import ClassVar

from . import qsc_hexane_post_birth_stepper as pbs
from .particle import coupled_transport as ct
from .particle import cut_event_orchestrator as ceo
from .particle import cut_extinction_event as cee
from .particle import cut_extinction_orchestrator as ceo_ext
from .particle import cut_extinction_projection as cep
from .particle import cut_extinction_time as cet
from .particle import cut_face_departure_integrator as di
from .particle import cut_integrator as ci
from .particle import cut_transport as cut

QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_ID = (
    "qsc-hexane-extinction-localizer-march-v1"
)
QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_REVISION = 1
QSC_HEXANE_EXTINCTION_LOCALIZER_ORCHESTRATION_ID = (
    "cell-zero-extinction-localizer-ladder-certificate-projection-v1"
)

#: Seed-scale floor for Stefan-flux scales, mirroring the march runners' own.
_STEFAN_SCALE_FLOOR_MOL_M2_S = 1.0e-3

#: OWNER-RULED (P-1, 2026-08-22, option (a); verbatim "P-1 I approve your
#: recommended option (a)"): THE ENVIRONMENT-CONDITIONED RESIDUAL BOUND FOR
#: THE PART-02 EXTINCTION CERTIFICATE, AND FOR NOTHING ELSE.  The frozen
#: ``2e-11`` positive-target residual contract stays untouched for every
#: other consumer.  Sized from the fresh 2026-08-22 floor re-measurement at
#: HEAD (``docs/GT_PS2_PART02_EXTINCTION_FLOOR_PROOF_RECORD_2026-08-22.md``):
#: the certified positive-target residual floor at the N=2 true extinction
#: terminus grows superlinearly (~sweep^2) with the single-solve sweep -
#: 2.062e-11 at 9 %, 1.064e-10 at 20 %, 1.170e-9 at 62.5 % - and the
#: certificate's own finite-limit contract (final target fraction <= 1e-6, a
#: ~100 % sweep) rides that floor toward its full-sweep plateau; the deepest
#: scanned rung (96.9 % sweep) measured
#: PART02_MEASURED_DEEP_LADDER_RESIDUAL_FLOOR below.  The declared bound is
#: that measured full-ladder maximum with an explicit margin factor (a
#: resolution bound over a BLAS/ULP-fragile cancellation band per BC-7
#: section 4.4, never a bit-pin): bound = measured maximum rounded up to one
#: significant figure after a 1.5x environment-drift margin
#: (3.007e-9 * 1.5 = 4.51e-9 -> 5e-9; the s^2 sweep law projects the 100 %
#: plateau at ~3.2-3.5e-9, also inside).  Consumed ONLY through
#: :func:`part02_environment_conditioned_resolution` by this composition's
#: certificate issuance path.
PART02_MEASURED_DEEP_LADDER_RESIDUAL_FLOOR = 3.007e-9
PART02_EXTINCTION_CERTIFICATE_ENVIRONMENT_CONDITIONED_RESIDUAL_BOUND = 5.0e-9

#: The frozen positive-target residual contract this declaration conditions.
PART02_FROZEN_POSITIVE_TARGET_RESIDUAL_CONTRACT = 2.0e-11

#: OWNER-RULED (P-2, 2026-08-23, verbatim "P-2 approved"; raised by P-3,
#: 2026-08-23, verbatim "Agreed totally, we need one more measured budget
#: raise, even if it means 10 h of compute. Proceed now, be as efficient as
#: possible."): THE DECLARED EVALUATION BUDGET FOR THE SAME SINGLE PART-02
#: CERTIFICATE PATH, AND FOR NOTHING ELSE.  The frozen 1200-evaluation
#: nonlinear budget of the certified N=2 trajectory state stays untouched
#: for every other consumer.  Sized from the measured deep-rung growth law
#: (floor record section 6): converged needs ~3.6e3 solver evaluations at
#: 87.5 % sweep and ~5.3e3 at 96.9 % (both xtol-terminated inside the ruled
#: residual resolution), with the 99.2 % rung exhausting the P-2 budget of
#: 6000 - the growth extrapolates the worst-rung need to ~9e3, and the
#: ruled 16000 carries ~1.8x margin over that extrapolation (a
#: measured-need bound, never a bit-pin; the P-3 stop rule forbids any
#: further raise without a further ruling).  Consumed ONLY through
#: :func:`part02_environment_conditioned_resolution`.
PART02_MEASURED_DEEP_RUNG_EVALUATION_NEED = 5300
PART02_EXTINCTION_CERTIFICATE_DECLARED_EVALUATION_BUDGET = 16000

#: The frozen nonlinear evaluation budget this declaration conditions (the
#: certified N=2 trajectory state's own controls value).
PART02_FROZEN_NONLINEAR_EVALUATION_BUDGET = 1200

_PART02_RULING = (
    "P-1 owner ruling 2026-08-22, option (a), verbatim 'P-1 I approve your "
    "recommended option (a)': PART-02 typed outcome at demonstration tier "
    "via an S7c-style proven-floor record; the extinction certificate "
    "carries a declared environment-conditioned residual bound sized from "
    "fresh measurement; the frozen 2e-11 contract stays untouched for every "
    "other consumer. P-2 owner ruling 2026-08-23, verbatim 'P-2 approved': "
    "the same declaration carries a declared evaluation budget for the same "
    "single certificate path (S7c escalated-retry pattern); the frozen "
    "1200-evaluation budget stays untouched for every other consumer. P-3 "
    "owner ruling 2026-08-23, verbatim 'Agreed totally, we need one more "
    "measured budget raise, even if it means 10 h of compute. Proceed now, "
    "be as efficient as possible.': the declared budget is raised to 16000 "
    "(~1.8x the extrapolated worst-rung need per the measured growth law "
    "3.6e3 -> 5.3e3 -> >6e3), with the stop rule that exhaustion at any "
    "rung stops typed with no further raise absent a further ruling"
)
_PART02_MEASUREMENT = (
    "fresh floor re-measurement 2026-08-22 at HEAD, N=2 true extinction "
    "terminus (front_z 2.9161622614389213e-09, t 14.721579691506928 s), "
    "plus the 2026-08-23 budget-floor probes (probe B: worst-case deep rung "
    "converges at ~3.6e3 evaluations, residual 1.798e-9): "
    "docs/GT_PS2_PART02_EXTINCTION_FLOOR_PROOF_RECORD_2026-08-22.md"
)


def part02_environment_conditioned_resolution():
    """Return THE single ruled declaration for the PART-02 certificate.

    Every value is the module-level declared constant; the certified layer
    refuses the declaration against any state whose frozen contract is not
    exactly ``PART02_FROZEN_POSITIVE_TARGET_RESIDUAL_CONTRACT``.
    """

    from .particle import cut_extinction_time as _cet

    return _cet.EnvironmentConditionedResidualResolution(
        declared_bound=(
            PART02_EXTINCTION_CERTIFICATE_ENVIRONMENT_CONDITIONED_RESIDUAL_BOUND
        ),
        frozen_contract_value=PART02_FROZEN_POSITIVE_TARGET_RESIDUAL_CONTRACT,
        ruling=_PART02_RULING,
        measurement=_PART02_MEASUREMENT,
        declared_maximum_function_evaluations=(
            PART02_EXTINCTION_CERTIFICATE_DECLARED_EVALUATION_BUDGET
        ),
        frozen_budget_value=PART02_FROZEN_NONLINEAR_EVALUATION_BUDGET,
    )


class ExtinctionLocalizerRefusalCode(Enum):
    NONFINITE_INPUT = "nonfinite_input"
    BOUNDARY_NOT_SUPPORTED = "boundary_not_supported"
    NOT_CELL_ZERO_POSITIVE_CORE = "not_cell_zero_positive_core"
    TERMINUS_NOT_TYPED = "terminus_not_typed"
    TERMINUS_STATE_MISMATCH = "terminus_state_mismatch"
    HISTORY_NOT_TYPED = "accepted_history_not_typed"
    HISTORY_STATE_MISMATCH = "accepted_history_state_mismatch"
    HISTORY_TANGENT_NOT_RECEDING = "accepted_history_tangent_not_receding"
    FILM_COMPANION_PRESENT = "film_companion_present_at_extinction"
    LADDER_NOT_TYPED = "ladder_controls_not_typed"
    REFINEMENT_NOT_TYPED = "refinement_controls_not_typed"
    LADDER_EXHAUSTED = "ladder_exhausted_without_finite_limit_certificate"
    LOCALIZATION_NOT_TYPED = "localization_payload_not_typed"
    PROJECTION_SEED_RANK_MISMATCH = "projection_seed_rank_mismatch"


class ExtinctionLocalizerError(ValueError):
    """Typed refusal for the extinction-localizer composition; exact rollback.

    No entry point in this module ever mutates its inputs, so a caller that
    catches this has lost nothing: the state it passed in is the state it
    still holds.
    """

    def __init__(
        self,
        code: ExtinctionLocalizerRefusalCode,
        message: str,
        *,
        state: object | None = None,
        evidence: dict | None = None,
    ) -> None:
        super().__init__(f"{code.value}: {message}")
        self.code = code
        self.state = state
        self.evidence = dict(evidence) if evidence else {}


_LOCALIZER_SEAL = object()


class _FalseLocalizerClaims:
    extinction_localizer_composition_only: ClassVar[bool] = True
    calls_certified_extinction_family_never_reimplements: ClassVar[bool] = True
    extinction_entry_film_free_asserted: ClassVar[bool] = True
    imported_by_production: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False


#: The two accepted step families whose committed motion may seed the ladder.
AcceptedExtinctionHistoryStep = ceo.AcceptedFaceHistoryStep


@dataclass(frozen=True)
class ExtinctionLadderControls:
    """Declared decreasing-target schedule for the localizer ladder.

    The schedule is geometry only: head fractions ease the first rungs away
    from the live front, then the tail contracts geometrically by
    ``tail_contraction`` per level until either the certified finite-limit
    contract is satisfied or ``maximum_levels`` is exhausted (typed refusal;
    never a widened tolerance).  Every rung is still the certified
    positive-core target solve with its own unchanged acceptance contracts.

    The default head is DENSE near the terminus - the certified extinction
    suite's own proven continuation pattern (0.99 -> 0.95 -> ...): the first
    live N=2 attempt measured that a sparse head (first target 0.75, a 25 %
    of-core sweep in one solve) exhausts the state's own nonlinear budget at
    9.8e-3 scaled residual, while march-scale first jumps converge.  The
    finite-limit contract constrains only the TAIL (last increments), so a
    dense head is lawful and adds evidence, never weakens it.
    """

    head_fractions: tuple[float, ...] = (
        0.99,
        0.97,
        0.94,
        0.90,
        0.85,
        0.78,
        0.70,
        0.60,
        0.50,
        0.375,
        0.25,
        0.125,
    )
    tail_contraction: float = 0.25
    maximum_levels: int = 40
    duration_bounds_s: tuple[float, float] = (1.0e-9, 1.0e4)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.head_fractions or not all(
            math.isfinite(value) and 0.0 < value < 1.0
            for value in self.head_fractions
        ):
            raise ValueError("ladder head fractions must lie strictly in (0,1)")
        if any(
            new >= old
            for old, new in zip(self.head_fractions, self.head_fractions[1:])
        ):
            raise ValueError("ladder head fractions must decrease strictly")
        if not math.isfinite(self.tail_contraction) or not (
            0.0 < self.tail_contraction < 1.0
        ):
            raise ValueError("ladder tail contraction must lie strictly in (0,1)")
        if (
            isinstance(self.maximum_levels, bool)
            or not isinstance(self.maximum_levels, int)
            or self.maximum_levels < len(self.head_fractions) + 1
        ):
            raise ValueError(
                "ladder maximum levels must be an integer covering the head "
                "plus at least one tail rung"
            )
        lower, upper = self.duration_bounds_s
        if not all(
            math.isfinite(value) and value > 0.0 for value in (lower, upper)
        ) or not lower < upper:
            raise ValueError(
                "ladder duration bounds must be positive, finite, and ordered"
            )


@dataclass(frozen=True, kw_only=True)
class CellZeroExtinctionLocalization(_FalseLocalizerClaims):
    """The certified finite-limit localization at one cell-zero terminus.

    ``steps`` are the live positive-core target roots exactly as the
    certified solver returned them (uncommitted; the shared ``before`` is the
    terminus state itself), and ``certificate`` is the certified finite-limit
    payload over exactly those roots.  Nothing is committed yet: commitment
    belongs to :func:`commit_cell_zero_extinction_outcome`.
    """

    before: ci.CutIntegratorState
    boundary: ct.PoreBoundary
    steps: tuple[cet.PositiveCoreTargetStep, ...]
    certificate: cet.ExtinctionTimeLimitCertificate
    ladder: ExtinctionLadderControls
    refinement: cet.ExtinctionTimeRefinementControls
    history_front_rate_z_s: float
    certification_attempts: int
    film_free_entry_asserted: bool
    _seal: InitVar[object] = None

    schema_id: ClassVar[str] = QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_ID
    schema_revision: ClassVar[int] = (
        QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_REVISION
    )
    orchestration_id: ClassVar[str] = (
        QSC_HEXANE_EXTINCTION_LOCALIZER_ORCHESTRATION_ID
    )

    def __post_init__(self, _seal: object) -> None:
        if _seal is not _LOCALIZER_SEAL:
            raise TypeError(
                "cell-zero extinction localizations are composition-sealed"
            )
        if self.certificate.before is not self.before:
            raise ValueError(
                "the localization certificate must identify the exact terminus"
            )
        if len(self.steps) != len(self.certificate.steps) or any(
            live is not carried
            for live, carried in zip(self.steps, self.certificate.steps)
        ):
            raise ValueError(
                "the localization roots must be identically the certificate's"
            )
        if self.film_free_entry_asserted is not True:
            raise ValueError(
                "a localization without the typed film-free entry assertion "
                "is not this composition's payload"
            )

    def journal_record(self) -> dict:
        certificate = self.certificate
        declared = self.refinement.environment_conditioned_resolution
        return {
            "schema_id": QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_ID,
            "record": "cell-zero-extinction-localization",
            # Both residual tiers, explicitly labeled, so no consumer can
            # mistake the environment-conditioned tier for the frozen one.
            "frozen_positive_target_residual_contract": (
                self.before.controls.nonlinear_residual_tolerance
            ),
            "environment_conditioned_residual_bound": (
                None if declared is None else declared.declared_bound
            ),
            "environment_conditioned_residual_ruling": (
                None if declared is None else declared.ruling
            ),
            "environment_conditioned_residual_measurement": (
                None if declared is None else declared.measurement
            ),
            "frozen_nonlinear_evaluation_budget": (
                self.before.controls.maximum_function_evaluations
            ),
            "declared_evaluation_budget": (
                None
                if declared is None
                else declared.declared_maximum_function_evaluations
            ),
            "front_z_at_terminus": self.before.transport.geometry.front.z,
            "history_front_rate_z_s": self.history_front_rate_z_s,
            "ladder_levels": len(self.steps),
            "certification_attempts": self.certification_attempts,
            "targets_z": [step.ledger.target_z for step in self.steps],
            "target_fractions": [
                step.ledger.target_fraction_of_initial_z for step in self.steps
            ],
            "event_durations_s": [
                step.ledger.event_duration_s for step in self.steps
            ],
            "maximum_scaled_residuals": [
                step.ledger.maximum_scaled_residual for step in self.steps
            ],
            "condition_proxies": [
                step.ledger.condition_proxy for step in self.steps
            ],
            "nonlinear_evaluations": [
                step.ledger.nonlinear_evaluations for step in self.steps
            ],
            "extrapolated_event_duration_s": (
                certificate.extrapolated_event_duration_s
            ),
            "last_positive_target_duration_s": (
                certificate.last_positive_target_duration_s
            ),
            "estimated_tail_bound_s": certificate.estimated_tail_bound_s,
            "maximum_observed_tail_increment_ratio": (
                certificate.maximum_observed_tail_increment_ratio
            ),
            "film_free_entry_asserted": self.film_free_entry_asserted,
            "physically_qualifying": False,
        }


@dataclass(frozen=True, kw_only=True)
class CellZeroExtinctionOutcome(_FalseLocalizerClaims):
    """THE TYPED OUTCOME: extinction certified and the after-state committed.

    ``macrostep`` is the certified extinction orchestrator's own atomic
    payload.  On the ``EXTINCTION_AND_FULLY_DRY_REMAINDER`` branch the
    committed ``after`` is a fully-dry transport state advanced by the
    existing fully-dry surface for exactly the post-extinction remainder; on
    the ``EXACT_EXTINCTION_ENDPOINT`` branch it is the exact ``z=0``
    committed endpoint.  Either way the trajectory ends in an OUTCOME, not a
    refusal - which is precisely the PART-02 criterion at this tier.
    """

    localization: CellZeroExtinctionLocalization
    macrostep: ceo_ext.ExtinctionMacrostep
    requested_dt_s: float
    _seal: InitVar[object] = None

    schema_id: ClassVar[str] = QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_ID
    schema_revision: ClassVar[int] = (
        QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_REVISION
    )
    orchestration_id: ClassVar[str] = (
        QSC_HEXANE_EXTINCTION_LOCALIZER_ORCHESTRATION_ID
    )

    def __post_init__(self, _seal: object) -> None:
        if _seal is not _LOCALIZER_SEAL:
            raise TypeError(
                "cell-zero extinction outcomes are composition-sealed"
            )
        if self.macrostep.before is not self.localization.before:
            raise ValueError(
                "the committed macrostep must begin at the localized terminus"
            )
        if self.macrostep.certificate is not self.localization.certificate:
            raise ValueError(
                "the committed macrostep must carry the localization's own "
                "certificate"
            )

    @property
    def branch(self) -> ceo_ext.ExtinctionMacrostepBranch:
        return self.macrostep.ledger.branch

    @property
    def after(self):
        return self.macrostep.after

    @property
    def dry_continuation_committed(self) -> bool:
        return self.macrostep.fully_dry_step is not None

    def journal_record(self) -> dict:
        ledger = self.macrostep.ledger
        record = {
            "schema_id": QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_ID,
            "record": "cell-zero-extinction-outcome",
            "branch": ledger.branch.value,
            "requested_dt_s": self.requested_dt_s,
            "numerical_extinction_duration_s": (
                ledger.numerical_extinction_duration_s
            ),
            "estimated_geometric_tail_s": ledger.estimated_geometric_tail_s,
            "fully_dry_remainder_s": ledger.fully_dry_remainder_s,
            "dry_continuation_committed": self.dry_continuation_committed,
            "duration_closure_error_s": ledger.duration_closure_error_s,
            "absolute_time_closure_error_s": (
                ledger.absolute_time_closure_error_s
            ),
            "event_ordering": list(ledger.event_ordering),
            "maximum_stage_scaled_residual": (
                ledger.maximum_stage_scaled_residual
            ),
            "maximum_stage_step_ledger_residual": (
                ledger.maximum_stage_step_ledger_residual
            ),
            "maximum_independent_cumulative_ledger_residual": (
                ledger.maximum_independent_cumulative_ledger_residual
            ),
            "positive_target_localization_root_count": (
                ledger.positive_target_localization_root_count
            ),
            "total_nonlinear_evaluations": ledger.total_nonlinear_evaluations,
            "after_time_s": self.after.time_s,
            "localization": self.localization.journal_record(),
            "physically_qualifying": False,
        }
        return record


def extinction_history_from_step(taken: object) -> (
    AcceptedExtinctionHistoryStep | None
):
    """Extract the ladder-seeding history one accepted march step provides.

    Identical extraction discipline to BC-7's
    ``face_event_history_from_step``: the committed macrostep's accepted
    same-cell or face-departure step is lawful accepted history; a stall
    routes nothing.
    """

    macrostep = getattr(taken, "macrostep", None)
    if macrostep is None:
        return None
    same_cell = getattr(macrostep, "same_cell_step", None)
    if isinstance(same_cell, ci.CutIntegratorStep):
        return same_cell
    departure = getattr(macrostep, "face_departure_step", None)
    if isinstance(departure, di.FaceDepartureIntegratorStep):
        return departure
    return None


def _require_supported_boundary(boundary: object) -> None:
    if not isinstance(
        boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.BOUNDARY_NOT_SUPPORTED,
            "the extinction localizer requires a supported pore boundary; "
            f"got {type(boundary).__name__}",
        )


def _require_cell_zero_state(state: object) -> ci.CutIntegratorState:
    if not isinstance(state, ci.CutIntegratorState):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.NONFINITE_INPUT,
            "the extinction localizer consumes a cut integrator state; got "
            f"{type(state).__name__}",
        )
    geometry = state.transport.geometry
    if (
        geometry.front.regime != "partial"
        or state.transport.layout.cut_cell_index != 0
        or not 0.0 < geometry.front.z < 1.0
    ):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.NOT_CELL_ZERO_POSITIVE_CORE,
            "extinction localisation begins only on the strict cell-zero "
            "positive-core chart",
            state=state,
        )
    return state


def _require_typed_terminus(
    terminus: object,
    state: ci.CutIntegratorState,
) -> ceo.CutEventMacrostepError:
    if not isinstance(terminus, ceo.CutEventMacrostepError) or (
        terminus.code is not ceo.CutMacrostepFailureCode.EXTINCTION_REQUIRED
    ):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.TERMINUS_NOT_TYPED,
            "the localizer routes exactly the orchestrator's own typed "
            "cell_zero_extinction_localizer_required terminus; got "
            f"{type(terminus).__name__}"
            + (
                f" with code {terminus.code.value!r}"
                if isinstance(terminus, ceo.CutEventMacrostepError)
                else ""
            ),
        )
    if terminus.rollback_state is not state:
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.TERMINUS_STATE_MISMATCH,
            "the typed terminus must carry exactly the marched state as its "
            "rollback (identity, not equality); a terminus from another "
            "state would localize the wrong trajectory",
            state=state,
        )
    return terminus


def _require_film_free(film_carrier: object, state: ci.CutIntegratorState) -> None:
    if film_carrier is None:
        return
    raise ExtinctionLocalizerError(
        ExtinctionLocalizerRefusalCode.FILM_COMPANION_PRESENT,
        "the certified crossing into cell zero was film-free and the "
        "extinction family defines no film-carrying chart; a caller holding "
        f"a {type(film_carrier).__name__} companion is refused - the "
        "companion is never silently dropped",
        state=state,
        evidence={
            "film_kg_per_kg_dry": getattr(
                film_carrier, "attached_hexane_kg_per_kg_dry", None
            )
        },
    )


def _history_front_rate(
    history: object,
    state: ci.CutIntegratorState,
) -> float:
    """Return the accepted history's own committed front rate (z/s).

    Accepted-trajectory evidence only: the committed front motion of the last
    ACCEPTED step ending in exactly the terminus state.  A non-receding rate
    refuses typed - a ladder seeded on a stalled or advancing tangent has no
    honest first duration.
    """

    if not isinstance(
        history,
        (ci.CutIntegratorStep, di.FaceDepartureIntegratorStep),
    ):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.HISTORY_NOT_TYPED,
            "the ladder seed consumes an accepted same-cell or "
            f"face-departure step; got {type(history).__name__} (a stall "
            "step has no receding tangent - the terminus after a stall "
            "cannot seed the ladder from it)",
        )
    if history.after is not state:
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.HISTORY_STATE_MISMATCH,
            "the accepted history must end in exactly the state being "
            "localized (identity, not equality)",
            state=state,
        )
    if not getattr(history.ledger, "accepted", False):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.HISTORY_NOT_TYPED,
            "the ladder seed consumes only an ACCEPTED history step; this "
            "ledger is not accepted",
            state=state,
        )
    # The identical committed-motion extraction the certified accepted-history
    # bracket layer performs (cut_event_orchestrator.
    # _accepted_history_face_event_payload): a same-cell step starts on its
    # strict cut state, a face-departure step starts on its exact-face
    # committed state; both end in exactly ``state``.
    if isinstance(history, ci.CutIntegratorStep):
        before_z = history.before.transport.geometry.front.z
    else:
        before_z = history.before.geometry.front.z
    after_z = state.transport.geometry.front.z
    duration = history.ledger.dt_s
    rate = math.fsum((after_z, -before_z)) / duration
    if not math.isfinite(rate) or rate >= 0.0:
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.HISTORY_TANGENT_NOT_RECEDING,
            "the accepted history's committed front rate is not strictly "
            f"receding ({rate!r} /s); the extinction ladder has no honest "
            "first-duration forecast from it",
            state=state,
            evidence={"front_rate_z_s": rate, "duration_s": duration},
        )
    return rate


def _first_target_seed(
    state: ci.CutIntegratorState,
    target_z: float,
    front_rate_z_s: float,
    duration_bounds_s: tuple[float, float],
) -> cet.PositiveCoreTargetSeed:
    """Seed the first rung from the state's OWN primitives and accepted rate.

    The candidate is the terminus state's differential primitives with the
    front placed exactly at the target (``cut_transport.
    candidate_from_state`` - the march runners' own accepted-state seed
    discipline); the duration forecast is the swept front distance over the
    accepted rate.  No rejected candidate, no clamp, no field tangent.
    """

    current_z = state.transport.geometry.front.z
    forecast = math.fsum((current_z, -target_z)) / (-front_rate_z_s)
    lower, upper = duration_bounds_s
    if not lower < forecast < upper:
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.NONFINITE_INPUT,
            "the accepted-rate first-rung forecast left the declared ladder "
            f"duration bounds: {forecast!r} s outside ({lower!r}, {upper!r})",
            state=state,
            evidence={"forecast_s": forecast},
        )
    candidate = cut.candidate_from_state(
        state.transport,
        dry_total_stefan_fluxes_mol_m2_s=(
            state.last_total_stefan_fluxes_mol_m2_s
        ),
        interface_temperature_k=state.last_interface_temperature_k,
        front_z=target_z,
    )
    scales = tuple(
        max(abs(value), _STEFAN_SCALE_FLOOR_MOL_M2_S)
        for value in state.last_total_stefan_fluxes_mol_m2_s
    )
    return cet.PositiveCoreTargetSeed(
        candidate,
        forecast,
        scales,
        (
            "cell-zero extinction ladder first rung: terminus-state "
            "primitives with the front at the declared target, duration "
            "forecast from the accepted history's committed front rate"
        ),
    )


def localize_cell_zero_extinction(
    state: ci.CutIntegratorState,
    boundary: ct.PoreBoundary,
    terminus: ceo.CutEventMacrostepError,
    *,
    history: AcceptedExtinctionHistoryStep,
    film_carrier: pbs.PostBirthFilmCarrier | None = None,
    ladder: ExtinctionLadderControls | None = None,
    refinement: cet.ExtinctionTimeRefinementControls | None = None,
    environment_conditioned_resolution: (
        cet.EnvironmentConditionedResidualResolution | None
    ) = None,
    on_rung=None,
) -> CellZeroExtinctionLocalization:
    """Run the certified decreasing-target ladder to a finite-limit certificate.

    Composition only: every rung is ``cut_extinction_time.
    solve_positive_core_target`` on the UNCHANGED equations, chained by
    ``cut_extinction_time.continuation_seed``; certification is
    ``cut_extinction_time.certify_finite_extinction_time`` under the module's
    own default refinement contract (never widened here; a caller may pass a
    different typed contract explicitly).  The ladder extends geometrically
    until the certificate exists or ``maximum_levels`` is exhausted (typed
    refusal carrying the last certification error).  Certified-layer refusals
    propagate untouched with exact rollback.  ``on_rung`` is an optional
    journaling hook called as ``on_rung(level_index, step)`` after each
    ACCEPTED rung - observation only, it receives the certified payload and
    its return value is ignored.
    """

    state = _require_cell_zero_state(state)
    _require_supported_boundary(boundary)
    _require_typed_terminus(terminus, state)
    _require_film_free(film_carrier, state)
    front_rate = _history_front_rate(history, state)
    if ladder is None:
        ladder = ExtinctionLadderControls()
    elif not isinstance(ladder, ExtinctionLadderControls):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.LADDER_NOT_TYPED,
            "the ladder schedule must be the typed ExtinctionLadderControls; "
            f"got {type(ladder).__name__}",
        )
    if environment_conditioned_resolution is not None and not isinstance(
        environment_conditioned_resolution,
        cet.EnvironmentConditionedResidualResolution,
    ):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.REFINEMENT_NOT_TYPED,
            "the environment-conditioned resolution must be the certified "
            "module's typed EnvironmentConditionedResidualResolution; got "
            f"{type(environment_conditioned_resolution).__name__}",
        )
    if refinement is None:
        refinement = cet.ExtinctionTimeRefinementControls(
            duration_bounds_s=ladder.duration_bounds_s,
            environment_conditioned_resolution=(
                environment_conditioned_resolution
            ),
        )
    elif not isinstance(refinement, cet.ExtinctionTimeRefinementControls):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.REFINEMENT_NOT_TYPED,
            "the finite-limit contract must be the certified module's typed "
            "ExtinctionTimeRefinementControls; got "
            f"{type(refinement).__name__}",
        )
    elif (
        environment_conditioned_resolution is not None
        and refinement.environment_conditioned_resolution
        != environment_conditioned_resolution
    ):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.REFINEMENT_NOT_TYPED,
            "the caller's refinement contract and its declared "
            "environment-conditioned resolution disagree; one declaration "
            "governs one certificate",
        )
    effective_resolution = refinement.environment_conditioned_resolution
    if on_rung is not None and not callable(on_rung):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.NONFINITE_INPUT,
            "the rung journaling hook must be callable or None",
        )

    initial_z = state.transport.geometry.front.z
    steps: list[cet.PositiveCoreTargetStep] = []
    fractions = list(ladder.head_fractions)
    certificate: cet.ExtinctionTimeLimitCertificate | None = None
    certification_attempts = 0
    last_certification_error: Exception | None = None
    level = 0
    while level < ladder.maximum_levels:
        if level < len(fractions):
            fraction = fractions[level]
        else:
            fraction = fractions[-1] * ladder.tail_contraction
            fractions.append(fraction)
        target = fraction * initial_z
        if steps:
            seed = cet.continuation_seed(
                steps[-1],
                target,
                label=(
                    "cell-zero extinction ladder rung "
                    f"{level + 1}: target fraction {fraction!r}"
                ),
            )
        else:
            seed = _first_target_seed(
                state,
                target,
                front_rate,
                ladder.duration_bounds_s,
            )
        steps.append(
            cet.solve_positive_core_target(
                state,
                target,
                boundary,
                seed,
                duration_bounds_s=ladder.duration_bounds_s,
                environment_conditioned_resolution=effective_resolution,
            )
        )
        if on_rung is not None:
            on_rung(level, steps[-1])
        level += 1
        if (
            len(steps) >= refinement.minimum_levels
            and fraction <= refinement.maximum_final_target_fraction
        ):
            certification_attempts += 1
            try:
                certificate = cet.certify_finite_extinction_time(
                    tuple(steps),
                    refinement,
                )
            except cet.ExtinctionTimeLimitError as exc:
                last_certification_error = exc
                continue
            break
    if certificate is None:
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.LADDER_EXHAUSTED,
            "the declared ladder was exhausted at "
            f"{len(steps)} accepted rungs without a certified finite "
            "extinction-time limit; the last certification refusal is "
            "chained (exact rollback: nothing was committed)",
            state=state,
            evidence={
                "accepted_rungs": len(steps),
                "certification_attempts": certification_attempts,
                "last_certification_error": (
                    None
                    if last_certification_error is None
                    else str(last_certification_error)[:500]
                ),
                "last_target_fraction": (
                    steps[-1].ledger.target_fraction_of_initial_z
                    if steps
                    else None
                ),
                "last_event_duration_s": (
                    steps[-1].ledger.event_duration_s if steps else None
                ),
            },
        ) from last_certification_error
    return CellZeroExtinctionLocalization(
        before=state,
        boundary=boundary,
        steps=tuple(steps),
        certificate=certificate,
        ladder=ladder,
        refinement=refinement,
        history_front_rate_z_s=front_rate,
        certification_attempts=certification_attempts,
        film_free_entry_asserted=True,
        _seal=_LOCALIZER_SEAL,
    )


def extinction_projection_seed_from_localization(
    localization: CellZeroExtinctionLocalization,
) -> cep.ExtinctionProjectionSeed:
    """Map the deepest accepted rung onto the exact-dry projection rank.

    Structural mapping only, on the certified layouts' own indices: the cut
    chart's dry pieces at a cell-zero cut are exactly the master cells
    (``cut_extinction_event.layout_for_extinction`` requires
    ``dry_cell_indices == range(n)``), so the rung's dry temperatures and
    compositions carry over cell by cell; the endpoint's positive-area faces
    are the master faces ``1..n``, so the rung's interface-face flux (whose
    face area vanishes at the endpoint) is dropped and the remaining face
    fluxes carry over face by face.  The event time is EXACTLY the certified
    extrapolated duration - ``project_certified_extinction`` refuses anything
    else.
    """

    if not isinstance(localization, CellZeroExtinctionLocalization):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.LOCALIZATION_NOT_TYPED,
            "the projection seed maps only this composition's sealed "
            f"localization; got {type(localization).__name__}",
        )
    before = localization.before
    layout = cee.layout_for_extinction(before.transport)
    deepest = localization.steps[-1].assembly.candidate
    n = layout.dry_cell_count
    if len(deepest.dry_temperatures_k) != n or len(
        deepest.dry_total_stefan_fluxes_mol_m2_s
    ) != n + 1:
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.PROJECTION_SEED_RANK_MISMATCH,
            "the deepest rung's dry rank does not map structurally onto the "
            "extinction endpoint rank "
            f"(dry pieces {len(deepest.dry_temperatures_k)}, fluxes "
            f"{len(deepest.dry_total_stefan_fluxes_mol_m2_s)}, endpoint "
            f"cells {n})",
            state=before,
        )
    candidate = cee.ExtinctionUnknowns(
        dry_temperatures_k=deepest.dry_temperatures_k,
        dry_y_hexane=deepest.dry_y_hexane,
        positive_area_total_stefan_fluxes_mol_m2_s=(
            deepest.dry_total_stefan_fluxes_mol_m2_s[1:]
        ),
        event_time_s=localization.certificate.extrapolated_event_duration_s,
    )
    scales = tuple(
        max(abs(value), _STEFAN_SCALE_FLOOR_MOL_M2_S)
        for value in candidate.positive_area_total_stefan_fluxes_mol_m2_s
    )
    return cep.ExtinctionProjectionSeed(
        candidate,
        scales,
        (
            "cell-zero extinction projection seed: the deepest accepted "
            "positive-target rung mapped structurally onto the endpoint "
            "rank at the certified extrapolated duration"
        ),
    )


def commit_cell_zero_extinction_outcome(
    localization: CellZeroExtinctionLocalization,
    dt_s: float,
    *,
    projection_seed: cep.ExtinctionProjectionSeed | None = None,
    boundary_discontinuity_times_s: tuple[float, ...] = (),
) -> CellZeroExtinctionOutcome:
    """Atomically commit the certified extinction and any dry remainder.

    Pure hand-off to ``cut_extinction_orchestrator.
    advance_numerical_extinction_macrostep`` - the certified transaction that
    validates the certificate against its live roots, projects the exact
    ``z=0`` endpoint, routes any remaining requested duration through the
    existing fully-dry surface on the label-preserving handoff, and enforces
    the macro conservation telescope.  Its typed refusals
    (``ExtinctionMacrostepError`` with exact rollback) propagate untouched.
    """

    if not isinstance(localization, CellZeroExtinctionLocalization):
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.LOCALIZATION_NOT_TYPED,
            "outcome commitment consumes only this composition's sealed "
            f"localization; got {type(localization).__name__}",
        )
    if type(dt_s) is not float or not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ExtinctionLocalizerError(
            ExtinctionLocalizerRefusalCode.NONFINITE_INPUT,
            "the outcome duration must be a strictly positive finite binary64",
        )
    if projection_seed is None:
        projection_seed = extinction_projection_seed_from_localization(
            localization
        )
    macrostep = ceo_ext.advance_numerical_extinction_macrostep(
        localization.before,
        dt_s,
        localization.boundary,
        localization.certificate,
        projection_seed,
        boundary_discontinuity_times_s=boundary_discontinuity_times_s,
    )
    return CellZeroExtinctionOutcome(
        localization=localization,
        macrostep=macrostep,
        requested_dt_s=dt_s,
        _seal=_LOCALIZER_SEAL,
    )


def resolve_cell_zero_extinction_terminus(
    state: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.PoreBoundary,
    terminus: ceo.CutEventMacrostepError,
    *,
    history: AcceptedExtinctionHistoryStep,
    film_carrier: pbs.PostBirthFilmCarrier | None = None,
    ladder: ExtinctionLadderControls | None = None,
    refinement: cet.ExtinctionTimeRefinementControls | None = None,
    environment_conditioned_resolution: (
        cet.EnvironmentConditionedResidualResolution | None
    ) = None,
) -> CellZeroExtinctionOutcome:
    """The composed terminus resolution: localize, certify, project, commit.

    The one-call mirror of BC-7's composed march step for the extinction
    terminus: on the orchestrator's typed
    ``cell_zero_extinction_localizer_required`` refusal, run the certified
    ladder to its finite-limit certificate and atomically commit the
    extinction (plus the fully-dry remainder when ``dt_s`` extends beyond the
    certified duration).  Every certified-layer refusal propagates untouched;
    this layer adds only its own typed gates and never mutates an input.
    """

    localization = localize_cell_zero_extinction(
        state,
        boundary,
        terminus,
        history=history,
        film_carrier=film_carrier,
        ladder=ladder,
        refinement=refinement,
        environment_conditioned_resolution=environment_conditioned_resolution,
    )
    return commit_cell_zero_extinction_outcome(localization, dt_s)


__all__ = [
    "PART02_EXTINCTION_CERTIFICATE_DECLARED_EVALUATION_BUDGET",
    "PART02_EXTINCTION_CERTIFICATE_ENVIRONMENT_CONDITIONED_RESIDUAL_BOUND",
    "PART02_FROZEN_NONLINEAR_EVALUATION_BUDGET",
    "PART02_FROZEN_POSITIVE_TARGET_RESIDUAL_CONTRACT",
    "PART02_MEASURED_DEEP_LADDER_RESIDUAL_FLOOR",
    "PART02_MEASURED_DEEP_RUNG_EVALUATION_NEED",
    "QSC_HEXANE_EXTINCTION_LOCALIZER_ORCHESTRATION_ID",
    "QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_ID",
    "QSC_HEXANE_EXTINCTION_LOCALIZER_SCHEMA_REVISION",
    "AcceptedExtinctionHistoryStep",
    "CellZeroExtinctionLocalization",
    "CellZeroExtinctionOutcome",
    "ExtinctionLadderControls",
    "ExtinctionLocalizerError",
    "ExtinctionLocalizerRefusalCode",
    "commit_cell_zero_extinction_outcome",
    "extinction_history_from_step",
    "extinction_projection_seed_from_localization",
    "localize_cell_zero_extinction",
    "part02_environment_conditioned_resolution",
    "resolve_cell_zero_extinction_terminus",
]
