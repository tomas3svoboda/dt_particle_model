r"""Atomic finite-extinction handoff for the Gate-1g cut solver.

The positive-core chart owns extinction-time localisation.  Once
:mod:`cut_extinction_time` has accepted a *numerical* finite ``z -> 0+`` limit
hypothesis under its observed geometric-contraction criteria, the separate
exact-dry projection owns the topology change.  The estimated geometric tail
is not a rigorous physical event-time/error interval: future contraction is
not independently bounded.  Consequently every payload in this module
remains ``physically_qualifying=False``.  If a requested macrostep extends
beyond the numerical projection duration, the remaining time is advanced by
the existing fully-dry equations.  All stages remain provisional until their
local and cumulative component/common-datum-energy ledgers and an independent
macro telescope pass.

No positive-radius cutoff, epsilon core, endpoint Rankine--Hugoniot equation,
or failure-selected contact is introduced here.  The existing fully-dry
predictor carries one residual-oil label, so a nonuniform labelled endpoint
may be committed exactly at extinction but cannot yet enter a post-extinction
remainder.  That unsupported handoff fails closed instead of homogenising the
material labels.

This module does not solve multiple interior master-face events.  Such a
sequence needs a face-to-face chart regular at both zero-volume endpoints.
Using an arbitrary positive-overlap departure solely to re-enter the ordinary
cut chart would make the result depend on a numerical split, so the one-face
orchestrator continues to reject that route explicitly.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum

from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_extinction_projection as projection
from dtdc_simulator.core2.particle import cut_extinction_time as extinction_time
from dtdc_simulator.core2.particle import cut_integrator as ci


class ExtinctionMacrostepBranch(str, Enum):
    """Accepted topology path through one requested duration."""

    EXACT_EXTINCTION_ENDPOINT = "exact_extinction_endpoint"
    EXTINCTION_AND_FULLY_DRY_REMAINDER = (
        "extinction_and_fully_dry_remainder"
    )


class ExtinctionMacrostepFailureCode(str, Enum):
    """Fail-closed reason owned by the extinction transaction."""

    BOUNDARY_DISCONTINUITY = "boundary_discontinuity_inside_macrostep"
    INVALID_LIMIT_CERTIFICATE = "invalid_finite_extinction_limit_certificate"
    EXTINCTION_AFTER_MACROSTEP = "extinction_after_macrostep"
    EXTINCTION_ORDERING_UNCERTIFIED = "extinction_ordering_uncertified"
    EXTINCTION_PROJECTION_FAILED = "exact_extinction_projection_failed"
    RETAINED_WATER_APPLICABILITY_GUARD_FAILED = (
        "wet_retained_water_applicability_guard_failed"
    )
    LABEL_PRESERVING_HANDOFF_REQUIRED = "label_preserving_fully_dry_handoff_required"
    FULLY_DRY_REMAINDER_FAILED = "fully_dry_remainder_failed"
    TOPOLOGY_HANDOFF_FAILED = "extinction_topology_handoff_failed"
    MACRO_LEDGER_FAILED = "extinction_macrostep_ledger_failed"
    INTERNAL_FAILURE = "unclassified_extinction_macrostep_failure"


class ExtinctionMacrostepError(RuntimeError):
    """Rejected transaction carrying the exact original rollback state."""

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState,
        *,
        code: ExtinctionMacrostepFailureCode,
        provisional_extinction_duration_s: float | None = None,
        nonlinear_evaluations: int = 0,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.code = code
        self.provisional_extinction_duration_s = (
            provisional_extinction_duration_s
        )
        self.nonlinear_evaluations = nonlinear_evaluations


@dataclass(frozen=True)
class ExtinctionMacrostepLedger:
    """Chronology, duration, work, and conservation for one transaction."""

    branch: ExtinctionMacrostepBranch
    requested_dt_s: float
    last_positive_target_duration_s: float
    numerical_extinction_duration_s: float
    estimated_geometric_tail_s: float
    fully_dry_remainder_s: float
    duration_closure_error_s: float
    absolute_time_closure_error_s: float
    accepted_balance_interval_count: int
    positive_target_localization_root_count: int
    positive_target_localization_nonlinear_evaluations: int
    accepted_nonlinear_evaluations: int
    total_nonlinear_evaluations: int
    event_ordering: tuple[str, ...]
    maximum_stage_scaled_residual: float
    maximum_stage_step_ledger_residual: float
    maximum_stage_cumulative_ledger_residual: float
    maximum_final_cumulative_ledger_residual: float
    maximum_independent_cumulative_ledger_residual: float
    water_step_residual_mol: float
    hexane_step_residual_mol: float
    energy_step_residual_j: float
    normalized_water_step_residual: float
    normalized_hexane_step_residual: float
    normalized_energy_step_residual: float
    water_cumulative_residual_mol: float
    hexane_cumulative_residual_mol: float
    energy_cumulative_residual_j: float
    normalized_water_cumulative_residual: float
    normalized_hexane_cumulative_residual: float
    normalized_energy_cumulative_residual: float
    boundary_water_out_mol: float
    boundary_hexane_out_mol: float
    boundary_energy_out_j: float
    accepted: bool
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def maximum_step_ledger_residual(self) -> float:
        return max(
            self.normalized_water_step_residual,
            self.normalized_hexane_step_residual,
            self.normalized_energy_step_residual,
        )

    @property
    def maximum_cumulative_ledger_residual(self) -> float:
        return max(
            self.normalized_water_cumulative_residual,
            self.normalized_hexane_cumulative_residual,
            self.normalized_energy_cumulative_residual,
        )


@dataclass(frozen=True)
class ExtinctionMacrostep:
    """One atomically accepted extinction endpoint or endpoint plus remainder."""

    before: ci.CutIntegratorState
    after: projection.ExtinctionCommittedState | ct.FullyDryTransportState
    boundary: ct.PoreBoundary
    certificate: extinction_time.ExtinctionTimeLimitCertificate
    projection_step: projection.ExtinctionProjectionStep
    fully_dry_step: ct.FullyDryTransportStep | None
    ledger: ExtinctionMacrostepLedger
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        has_remainder = self.fully_dry_step is not None
        expected = (
            self.ledger.branch
            is ExtinctionMacrostepBranch.EXTINCTION_AND_FULLY_DRY_REMAINDER
        )
        if has_remainder != expected:
            raise ValueError("extinction branch and remainder payload disagree")
        if not self.ledger.accepted:
            raise ValueError("an accepted extinction transaction needs its ledger")
        if self.projection_step.before is not self.before:
            raise ValueError("extinction projection changed the transaction origin")
        if has_remainder:
            assert self.fully_dry_step is not None
            if self.fully_dry_step.after is not self.after:
                raise ValueError("fully-dry remainder is not the retained endpoint")
        elif self.projection_step.after is not self.after:
            raise ValueError("exact extinction is not the retained endpoint")


class _AtomicRejection(RuntimeError):
    def __init__(self, code: ExtinctionMacrostepFailureCode, message: str) -> None:
        super().__init__(message)
        self.code = code


def advance_numerical_extinction_macrostep(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.PoreBoundary,
    certificate: extinction_time.ExtinctionTimeLimitCertificate,
    projection_seed: projection.ExtinctionProjectionSeed,
    *,
    boundary_discontinuity_times_s: Sequence[float] = (),
) -> ExtinctionMacrostep:
    """Commit a nonqualifying numerical extinction and any dry remainder.

    ``certificate`` is the existing numerical finite-limit payload; its name
    does not make its estimated geometric tail a rigorous event-time error
    bound.  It and every positive-target root must identify the exact immutable
    ``before`` state and the same constant boundary.  A known boundary change
    strictly inside the requested interval rejects before any projection work;
    callers must split exactly at that forcing timestamp and build a numerical
    limit hypothesis for the applicable constant-boundary interval.
    """

    if not isinstance(before, ci.CutIntegratorState):
        raise TypeError("extinction transaction requires CutIntegratorState")
    provisional_duration: float | None = None
    attempted_evaluations = 0
    try:
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("extinction transaction needs one supported pore boundary")
        if not math.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("extinction macrostep duration must be positive")
        ci.certify_enabled_wet_retained_water_state_applicability(
            before.controls,
            before.transport,
        )
        _reject_internal_boundary_discontinuity(
            before.transport.time_s,
            dt_s,
            boundary_discontinuity_times_s,
        )
        _validate_certificate(before, boundary, certificate)
        if not isinstance(projection_seed, projection.ExtinctionProjectionSeed):
            raise _AtomicRejection(
                ExtinctionMacrostepFailureCode.EXTINCTION_PROJECTION_FAILED,
                "exact z=0 projection needs its fixed-time endpoint seed payload",
            )
        provisional_duration = certificate.extrapolated_event_duration_s
        lower_duration = certificate.last_positive_target_duration_s
        event_duration = certificate.extrapolated_event_duration_s

        if dt_s < lower_duration:
            raise _AtomicRejection(
                ExtinctionMacrostepFailureCode.EXTINCTION_AFTER_MACROSTEP,
                "the last accepted positive-core target lies after the requested "
                "macrostep endpoint",
            )
        if dt_s < event_duration:
            raise _AtomicRejection(
                ExtinctionMacrostepFailureCode.EXTINCTION_ORDERING_UNCERTIFIED,
                "the requested endpoint lies inside the unresolved positive-core "
                "geometric tail; refine/split rather than assign event ordering",
            )

        try:
            endpoint = extinction_time.project_certified_extinction(
                certificate,
                projection_seed,
            )
        except (
            extinction_time.ExtinctionTimeLimitError,
            projection.ExtinctionProjectionStepError,
        ) as exc:
            attempted_evaluations += int(
                getattr(exc, "nonlinear_evaluations", 0)
            )
            rollback = getattr(exc, "rollback_state", before)
            if rollback is not before:
                raise _AtomicRejection(
                    ExtinctionMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                    "extinction projection lost the original rollback object",
                ) from exc
            if getattr(
                exc,
                "retained_water_applicability_guard_failed",
                False,
            ):
                raise _AtomicRejection(
                    ExtinctionMacrostepFailureCode.RETAINED_WATER_APPLICABILITY_GUARD_FAILED,
                    "extinction projection left the selected smooth retained-water scope",
                ) from exc
            raise _AtomicRejection(
                ExtinctionMacrostepFailureCode.EXTINCTION_PROJECTION_FAILED,
                f"exact z=0 projection did not certify: {exc}",
            ) from exc

        attempted_evaluations += int(endpoint.ledger.nonlinear_evaluations)
        _validate_projection_handoff(
            before,
            boundary,
            certificate,
            endpoint,
        )
        remainder = math.fsum((dt_s, -event_duration))
        if remainder < 0.0:
            raise _AtomicRejection(
                ExtinctionMacrostepFailureCode.EXTINCTION_ORDERING_UNCERTIFIED,
                "numerical extinction duration left a negative floating-point "
                "remainder",
            )

        if remainder == 0.0:
            branch = ExtinctionMacrostepBranch.EXACT_EXTINCTION_ENDPOINT
            dry_step = None
            after = endpoint.after
            stage_ledgers = (endpoint.ledger,)
            ordering = (
                "positive_core_numerical_finite_limit_hypothesis",
                "fixed_time_exact_z_zero_projection",
            )
        else:
            try:
                fully_dry_before = endpoint.after.to_uniform_fully_dry_transport_state()
            except projection.ExtinctionProjectionCompatibilityError as exc:
                raise _AtomicRejection(
                    ExtinctionMacrostepFailureCode.LABEL_PRESERVING_HANDOFF_REQUIRED,
                    "the extinction endpoint has nonuniform residual-oil labels; "
                    "post-extinction continuation needs a label-preserving fully-dry "
                    "solver and must not homogenise them",
                ) from exc
            try:
                dry_step = ct.advance_fully_dry_backward_euler(
                    fully_dry_before,
                    remainder,
                    boundary,
                )
            except ct.CoupledTransportStepError as exc:
                attempted_evaluations += int(exc.nonlinear_evaluations)
                if exc.rollback_state is not fully_dry_before:
                    raise _AtomicRejection(
                        ExtinctionMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                        "fully-dry rejection lost the provisional dry endpoint",
                    ) from exc
                raise _AtomicRejection(
                    ExtinctionMacrostepFailureCode.FULLY_DRY_REMAINDER_FAILED,
                    f"post-extinction remainder did not certify: {exc}",
                ) from exc
            attempted_evaluations += int(dry_step.ledger.nonlinear_evaluations)
            if dry_step.before is not fully_dry_before:
                raise _AtomicRejection(
                    ExtinctionMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                    "fully-dry continuation changed its exact handoff state",
                )
            branch = (
                ExtinctionMacrostepBranch.EXTINCTION_AND_FULLY_DRY_REMAINDER
            )
            after = dry_step.after
            stage_ledgers = (endpoint.ledger, dry_step.ledger)
            ordering = (
                "positive_core_numerical_finite_limit_hypothesis",
                "fixed_time_exact_z_zero_projection",
                "fully_dry_requested_remainder",
            )

        ledger = _compose_ledger(
            before,
            after,
            dt_s,
            certificate,
            remainder,
            branch,
            stage_ledgers,
            ordering,
        )
        ci.certify_enabled_wet_retained_water_state_applicability(
            before.controls,
            before.transport,
        )
        return ExtinctionMacrostep(
            before,
            after,
            boundary,
            certificate,
            endpoint,
            dry_step,
            ledger,
        )
    except ExtinctionMacrostepError:
        raise
    except _AtomicRejection as exc:
        raise ExtinctionMacrostepError(
            f"atomic extinction macrostep rejected with exact rollback: {exc}",
            before,
            code=exc.code,
            provisional_extinction_duration_s=provisional_duration,
            nonlinear_evaluations=attempted_evaluations,
        ) from exc
    except ci.WetRetainedWaterApplicabilityError as exc:
        raise ExtinctionMacrostepError(
            f"atomic extinction macrostep rejected with exact rollback: {exc}",
            before,
            code=(
                ExtinctionMacrostepFailureCode.RETAINED_WATER_APPLICABILITY_GUARD_FAILED
            ),
            provisional_extinction_duration_s=provisional_duration,
            nonlinear_evaluations=attempted_evaluations,
        ) from exc
    except Exception as exc:
        raise ExtinctionMacrostepError(
            f"atomic extinction macrostep rejected with exact rollback: {exc}",
            before,
            code=ExtinctionMacrostepFailureCode.INTERNAL_FAILURE,
            provisional_extinction_duration_s=provisional_duration,
            nonlinear_evaluations=attempted_evaluations,
        ) from exc


def _reject_internal_boundary_discontinuity(
    start_time_s: float,
    dt_s: float,
    discontinuities_s: Sequence[float],
) -> None:
    times = tuple(float(value) for value in discontinuities_s)
    if not all(math.isfinite(value) and value >= 0.0 for value in times):
        raise ValueError("boundary discontinuity times must be finite and non-negative")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ValueError("boundary discontinuity times must increase strictly")
    endpoint = math.fsum((start_time_s, dt_s))
    for value in times:
        if start_time_s < value < endpoint:
            raise _AtomicRejection(
                ExtinctionMacrostepFailureCode.BOUNDARY_DISCONTINUITY,
                "a forcing discontinuity lies strictly inside the extinction "
                "macrostep; split at its exact timestamp",
            )


def _validate_certificate(
    before: ci.CutIntegratorState,
    boundary: ct.PoreBoundary,
    certificate: extinction_time.ExtinctionTimeLimitCertificate,
) -> None:
    if not isinstance(
        certificate,
        extinction_time.ExtinctionTimeLimitCertificate,
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.INVALID_LIMIT_CERTIFICATE,
            "extinction routing requires the positive-core finite-limit payload",
        )
    if certificate.before is not before or certificate.boundary != boundary:
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.INVALID_LIMIT_CERTIFICATE,
            "numerical finite-limit payload identifies another state or boundary",
        )
    geometry = before.transport.geometry
    if (
        geometry.front.regime != "partial"
        or before.transport.layout.cut_cell_index != 0
        or not 0.0 < geometry.front.z < 1.0
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.INVALID_LIMIT_CERTIFICATE,
            "finite extinction begins only on the strict cell-zero positive-core chart",
        )
    roots = tuple(certificate.steps)
    if not roots:
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.INVALID_LIMIT_CERTIFICATE,
            "numerical finite-limit payload contains no positive-core roots",
        )
    if any(
        not isinstance(step, extinction_time.PositiveCoreTargetStep)
        for step in roots
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.INVALID_LIMIT_CERTIFICATE,
            "numerical finite-limit payload contains a non-live target root",
        )
    try:
        reproduced = extinction_time.certify_finite_extinction_time(
            roots,
            certificate.controls,
        )
    except (
        TypeError,
        ValueError,
        extinction_time.ExtinctionTimeLimitError,
    ) as exc:
        if getattr(
            exc,
            "retained_water_applicability_guard_failed",
            False,
        ):
            raise _AtomicRejection(
                ExtinctionMacrostepFailureCode.RETAINED_WATER_APPLICABILITY_GUARD_FAILED,
                "positive-core extinction evidence left the selected smooth "
                "retained-water scope",
            ) from exc
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.INVALID_LIMIT_CERTIFICATE,
            f"positive-core finite-limit evidence did not reproduce: {exc}",
        ) from exc
    if (
        reproduced.before is not certificate.before
        or reproduced.boundary != certificate.boundary
        or len(reproduced.steps) != len(certificate.steps)
        or any(
            live is not carried
            for live, carried in zip(reproduced.steps, certificate.steps)
        )
        or reproduced.extrapolated_event_duration_s
        != certificate.extrapolated_event_duration_s
        or reproduced.last_positive_target_duration_s
        != certificate.last_positive_target_duration_s
        or reproduced.estimated_tail_bound_s
        != certificate.estimated_tail_bound_s
        or reproduced.maximum_observed_tail_increment_ratio
        != certificate.maximum_observed_tail_increment_ratio
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.INVALID_LIMIT_CERTIFICATE,
            "numerical finite-limit payload does not reproduce from its live roots",
        )
    last = certificate.last_positive_target_duration_s
    event = certificate.extrapolated_event_duration_s
    tail = certificate.estimated_tail_bound_s
    if not all(math.isfinite(value) for value in (last, event, tail)):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.INVALID_LIMIT_CERTIFICATE,
            "numerical finite-limit chronology is non-finite",
        )
    tail_closure = math.fsum((last, tail, -event))
    tail_roundoff = 4.0 * max(
        math.ulp(max(abs(last), 1.0)),
        math.ulp(max(abs(event), 1.0)),
        math.ulp(max(abs(tail), 1.0)),
    )
    if (
        roots[-1].ledger.event_duration_s != last
        or last <= 0.0
        or event < last
        or tail < 0.0
        or abs(tail_closure) > tail_roundoff
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.INVALID_LIMIT_CERTIFICATE,
            "numerical finite-limit duration and estimated tail do not telescope to "
            "floating-point roundoff",
        )


def _validate_projection_handoff(
    before: ci.CutIntegratorState,
    boundary: ct.PoreBoundary,
    certificate: extinction_time.ExtinctionTimeLimitCertificate,
    step: projection.ExtinctionProjectionStep,
) -> None:
    duration = certificate.extrapolated_event_duration_s
    if (
        step.before is not before
        or step.boundary != boundary
        or step.ledger.event_duration_s != duration
        or not step.ledger.accepted
        or step.after.time_s != math.fsum((before.transport.time_s, duration))
        or step.after.grid is not before.transport.geometry.master_grid
        or step.after.accepted_steps != before.accepted_steps + 1
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
            "exact projection payload does not match the numerical limit duration",
        )
    if (
        step.assembly.endpoint_geometry.front.z != 0.0
        or step.assembly.endpoint_geometry.front.radius_m != 0.0
        or step.assembly.endpoint_geometry.front.regime != "fully_dry"
        or step.assembly.endpoint_geometry.cut_cell_index is not None
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
            "extinction projection did not land on exact fully-dry geometry",
        )


def _compensated_vector_total(
    high: Sequence[float],
    compensation: Sequence[float],
) -> float:
    if len(high) != len(compensation):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "extinction cumulative vector and its compensation do not align",
        )
    return math.fsum(
        math.fsum((float(value), float(correction)))
        for value, correction in zip(high, compensation)
    )


def _compose_ledger(
    before: ci.CutIntegratorState,
    after: projection.ExtinctionCommittedState | ct.FullyDryTransportState,
    requested_dt_s: float,
    certificate: extinction_time.ExtinctionTimeLimitCertificate,
    remainder_s: float,
    branch: ExtinctionMacrostepBranch,
    stage_ledgers: Sequence[object],
    ordering: tuple[str, ...],
) -> ExtinctionMacrostepLedger:
    ledgers = tuple(stage_ledgers)
    if not ledgers or not all(
        (
            bool(stage.accepted)
            if hasattr(stage, "accepted")
            else isinstance(stage, ct.FullyDryStepLedger)
        )
        for stage in ledgers
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "only accepted balance stages may enter the extinction telescope",
        )
    expected_stage_count = (
        1
        if branch is ExtinctionMacrostepBranch.EXACT_EXTINCTION_ENDPOINT
        else 2
    )
    if len(ledgers) != expected_stage_count:
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "extinction branch and provisional balance-stage count disagree",
        )
    dry_ledger = ledgers[1] if len(ledgers) == 2 else None

    def stage_scaled(stage: object) -> float:
        if hasattr(stage, "maximum_scaled_residual"):
            return float(stage.maximum_scaled_residual)
        return max(
            float(stage.max_scaled_component_residual),
            float(stage.max_scaled_energy_residual),
        )

    def stage_step(stage: object) -> float:
        if hasattr(stage, "maximum_step_ledger_residual"):
            return float(stage.maximum_step_ledger_residual)
        return float(stage.max_normalized_accepted_step_ledger)

    def stage_cumulative(stage: object) -> float:
        if hasattr(stage, "maximum_cumulative_ledger_residual"):
            return float(stage.maximum_cumulative_ledger_residual)
        return float(stage.max_normalized_cumulative_ledger)

    water_step = math.fsum(
        float(stage.water_step_residual_mol) for stage in ledgers
    )
    hexane_step = math.fsum(
        float(stage.hexane_step_residual_mol) for stage in ledgers
    )
    energy_step = math.fsum(float(stage.energy_step_residual_j) for stage in ledgers)
    boundary_water = math.fsum(
        float(stage.boundary_water_out_mol) for stage in ledgers
    )
    boundary_hexane = math.fsum(
        float(stage.boundary_hexane_out_mol) for stage in ledgers
    )
    boundary_energy = math.fsum(
        float(stage.boundary_energy_out_j) for stage in ledgers
    )
    water_scale = max(
        before.reference_inventory.total_water_mol,
        abs(boundary_water),
        1.0e-300,
    )
    hexane_scale = max(
        before.reference_inventory.total_hexane_mol,
        abs(boundary_hexane),
        1.0e-300,
    )
    energy_scale = max(
        before.reference_capacity_energy_scale_j,
        abs(boundary_energy),
        1.0e-300,
    )
    normalized_step = (
        abs(water_step) / water_scale,
        abs(hexane_step) / hexane_scale,
        abs(energy_step) / energy_scale,
    )
    maximum_stage_scaled = max(stage_scaled(stage) for stage in ledgers)
    maximum_stage_step = max(stage_step(stage) for stage in ledgers)
    maximum_stage_cumulative = max(
        stage_cumulative(stage) for stage in ledgers
    )
    maximum_final_cumulative = stage_cumulative(ledgers[-1])

    if dry_ledger is None:
        assert isinstance(after, projection.ExtinctionCommittedState) or hasattr(
            after,
            "cumulative_material_water_change_cell_mol",
        )
        water_material = _compensated_vector_total(
            after.cumulative_material_water_change_cell_mol,
            after.water_change_compensation_cell_mol,
        )
        hexane_material = _compensated_vector_total(
            after.cumulative_material_hexane_change_cell_mol,
            after.hexane_change_compensation_cell_mol,
        )
        energy_material = _compensated_vector_total(
            after.cumulative_material_energy_change_cell_j,
            after.energy_change_compensation_cell_j,
        )
        cumulative_boundary_water = (
            after.corrected_cumulative_boundary_water_out_mol
        )
        cumulative_boundary_hexane = (
            after.corrected_cumulative_boundary_hexane_out_mol
        )
        cumulative_boundary_energy = (
            after.corrected_cumulative_boundary_energy_out_j
        )
    else:
        assert isinstance(after, ct.FullyDryTransportState)
        water_material = after.cumulative_water_inventory_change_mol
        hexane_material = after.cumulative_hexane_inventory_change_mol
        energy_material = after.cumulative_energy_inventory_change_j
        cumulative_boundary_water = after.cumulative_boundary_water_out_mol
        cumulative_boundary_hexane = after.cumulative_boundary_hexane_out_mol
        cumulative_boundary_energy = after.cumulative_boundary_energy_out_j

    water_cumulative = math.fsum(
        (water_material, cumulative_boundary_water)
    )
    hexane_cumulative = math.fsum(
        (hexane_material, cumulative_boundary_hexane)
    )
    energy_cumulative = math.fsum(
        (energy_material, cumulative_boundary_energy)
    )
    cumulative_water_scale = max(
        before.reference_inventory.total_water_mol,
        abs(water_material),
        abs(cumulative_boundary_water),
        after.cumulative_absolute_water_transfer_mol,
        1.0e-300,
    )
    cumulative_hexane_scale = max(
        before.reference_inventory.total_hexane_mol,
        abs(hexane_material),
        abs(cumulative_boundary_hexane),
        after.cumulative_absolute_hexane_transfer_mol,
        1.0e-300,
    )
    cumulative_energy_scale = max(
        before.reference_capacity_energy_scale_j,
        abs(energy_material),
        abs(cumulative_boundary_energy),
        after.cumulative_absolute_energy_transfer_j,
        1.0e-300,
    )
    normalized_cumulative = (
        abs(water_cumulative) / cumulative_water_scale,
        abs(hexane_cumulative) / cumulative_hexane_scale,
        abs(energy_cumulative) / cumulative_energy_scale,
    )
    maximum_independent_cumulative = max(normalized_cumulative)

    event_duration = certificate.extrapolated_event_duration_s
    duration_closure = math.fsum(
        (event_duration, remainder_s, -requested_dt_s)
    )
    target_time = math.fsum((before.transport.time_s, requested_dt_s))
    after_time = after.time_s
    absolute_time_closure = math.fsum((after_time, -target_time))
    time_roundoff = 8.0 * max(
        math.ulp(max(abs(before.transport.time_s), 1.0)),
        math.ulp(max(abs(target_time), 1.0)),
        math.ulp(max(abs(requested_dt_s), 1.0)),
    )
    if (
        abs(duration_closure) > time_roundoff
        or abs(absolute_time_closure) > time_roundoff
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
            "extinction and remainder do not close the exact requested duration",
        )

    tolerance = min(before.controls.ledger_tolerance, 1.0e-10)
    # P-1 ruling 2026-08-22 option (a): a certificate carrying the typed
    # environment-conditioned residual declaration conditions ONLY the
    # nonlinear residual-magnitude clause of ITS OWN stages; the certified
    # per-stage solvers keep their own internal gates, and every other
    # certificate (declaration None) takes the bit-identical frozen path.
    declared = getattr(
        certificate.controls,
        "environment_conditioned_resolution",
        None,
    )
    if declared is None:
        nonlinear_tolerance = before.controls.nonlinear_residual_tolerance
    else:
        if (
            declared.frozen_contract_value
            != before.controls.nonlinear_residual_tolerance
        ):
            raise _AtomicRejection(
                ExtinctionMacrostepFailureCode.MACRO_LEDGER_FAILED,
                "the certificate's environment-conditioned declaration does "
                "not condition this state's frozen nonlinear contract",
            )
        nonlinear_tolerance = declared.declared_bound
    if maximum_stage_scaled > nonlinear_tolerance:
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "an extinction stage exceeded the unchanged nonlinear contract"
            if declared is None
            else (
                "an extinction stage exceeded even the declared "
                "environment-conditioned resolution"
            ),
        )
    if (
        maximum_stage_step > tolerance
        or maximum_stage_cumulative > tolerance
        or max(normalized_step) > tolerance
        or maximum_independent_cumulative > tolerance
    ):
        raise _AtomicRejection(
            ExtinctionMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "stage, macrostep, or cumulative conservation exceeded 1e-10",
        )

    localization_evaluations = sum(
        int(step.ledger.nonlinear_evaluations) for step in certificate.steps
    )
    accepted_evaluations = sum(
        int(stage.nonlinear_evaluations) for stage in ledgers
    )
    return ExtinctionMacrostepLedger(
        branch=branch,
        requested_dt_s=requested_dt_s,
        last_positive_target_duration_s=(
            certificate.last_positive_target_duration_s
        ),
        numerical_extinction_duration_s=event_duration,
        estimated_geometric_tail_s=certificate.estimated_tail_bound_s,
        fully_dry_remainder_s=remainder_s,
        duration_closure_error_s=duration_closure,
        absolute_time_closure_error_s=absolute_time_closure,
        accepted_balance_interval_count=len(ledgers),
        positive_target_localization_root_count=len(certificate.steps),
        positive_target_localization_nonlinear_evaluations=(
            localization_evaluations
        ),
        accepted_nonlinear_evaluations=accepted_evaluations,
        total_nonlinear_evaluations=(
            localization_evaluations + accepted_evaluations
        ),
        event_ordering=ordering,
        maximum_stage_scaled_residual=maximum_stage_scaled,
        maximum_stage_step_ledger_residual=maximum_stage_step,
        maximum_stage_cumulative_ledger_residual=maximum_stage_cumulative,
        maximum_final_cumulative_ledger_residual=maximum_final_cumulative,
        maximum_independent_cumulative_ledger_residual=(
            maximum_independent_cumulative
        ),
        water_step_residual_mol=water_step,
        hexane_step_residual_mol=hexane_step,
        energy_step_residual_j=energy_step,
        normalized_water_step_residual=normalized_step[0],
        normalized_hexane_step_residual=normalized_step[1],
        normalized_energy_step_residual=normalized_step[2],
        water_cumulative_residual_mol=water_cumulative,
        hexane_cumulative_residual_mol=hexane_cumulative,
        energy_cumulative_residual_j=energy_cumulative,
        normalized_water_cumulative_residual=normalized_cumulative[0],
        normalized_hexane_cumulative_residual=normalized_cumulative[1],
        normalized_energy_cumulative_residual=normalized_cumulative[2],
        boundary_water_out_mol=boundary_water,
        boundary_hexane_out_mol=boundary_hexane,
        boundary_energy_out_j=boundary_energy,
        accepted=True,
    )


__all__ = [
    "ExtinctionMacrostep",
    "ExtinctionMacrostepBranch",
    "ExtinctionMacrostepError",
    "ExtinctionMacrostepFailureCode",
    "ExtinctionMacrostepLedger",
    "advance_numerical_extinction_macrostep",
]
