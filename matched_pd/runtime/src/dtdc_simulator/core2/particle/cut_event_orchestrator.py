r"""Atomic one-interior-face macrosteps for the Gate-1g cut solver.

This module composes the existing fixed-rank solvers without changing their
equations.  A requested backward-Euler macrostep first attempts the ordinary
same-cell chart.  If that chart has no accepted root, an interior master-face
arrival may be localized and the unconsumed duration solved from the exact
face by the coalescence-regular departure chart.

The two event equations form a triangular composite system: the exact-arrival
root depends only on the original state, while the departure root depends on
that arrival state.  Both roots are provisional until the complete requested
duration and an independently recomputed macro ledger pass.  A later failure
therefore rolls back to the exact original :class:`CutIntegratorState`; an
exact-face state is never exposed as a committed partial result.

Only one interior face may be crossed.  Cell-zero extinction, a second face
inside the remainder, and a boundary discontinuity strictly inside the
requested interval are rejected explicitly.  Event ordering has two distinct
routes: an independent absolute-error certificate, or explicitly
nonqualifying ``open_interval_numerical_ordering`` evidence.  The latter is
available only when the unknown-time arrival is solved directly on the exact
open chart ``(0, requested_dt)`` and a fresh physical-Jacobian Newton
correction remains representably inside that chart.  It carries no interval
error bound and makes no physical-qualification claim.

This deliberately narrow numerical interpretation follows the one-sided
event-point/consistency principles of Lopez and Maset, J. Sci. Comput. 100:3
(2024), DOI 10.1007/s10915-024-02546-w.  Unlike the sign brackets refined by
SUNDIALS IDA and PETSc TS, the direct coupled solve does not enclose a true
root; consequently its evidence must never be relabelled a rigorous event-time
certificate.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
import dataclasses
from dataclasses import dataclass, field
from enum import Enum

from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_face_event as face
from dtdc_simulator.core2.particle import cut_face_departure_integrator as di
from dtdc_simulator.core2.particle import cut_face_event_integrator as fi
from dtdc_simulator.core2.particle import cut_face_to_face_integrator as ffi
from dtdc_simulator.core2.particle import cut_face_tangent as tangent
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_retained_cap as wrc


# This is a nonlinear-coordinate scale, not a flux floor or constitutive
# clamp.  It matches the independent physical-column scale used by the exact
# face-arrival solver.
ACCEPTED_HISTORY_STEFAN_SCALE_FLOOR_MOL_M2_S = 1.0e-3


AcceptedFaceHistoryStep = ci.CutIntegratorStep | di.FaceDepartureIntegratorStep


@dataclass(frozen=True)
class AcceptedHistoryFaceEventBracket:
    """Deterministic face-routing input derived from one accepted motion.

    The word ``bracket`` is deliberately kinematic: the current accepted front
    lies above its inner face while the constant accepted-history tangent lies
    below that face at the requested endpoint.  It is not an interval proof of
    the nonlinear event time and therefore carries no physical-qualification
    claim.  The exact arrival solve and its separate ordering evidence remain
    mandatory before topology can change.

    Every seed primitive is reconstructed from ``source_step`` without
    clipping.  The same-cell equations receive the exact current primitives
    with only their excluded current-front endpoint represented at the neutral
    temporal-continuation point of the unchanged chart.  No exact-face seed is
    representable in this forecast: only after same-cell rejection does the
    production wrapper project the current accepted state onto exact-face
    rank.  Rejected nonlinear iterates are never seed authority.
    """

    source_step: AcceptedFaceHistoryStep
    before: ci.CutIntegratorState
    requested_dt_s: float
    source_dt_s: float
    source_start_front_z: float
    source_end_front_z: float
    inferred_front_velocity_z_s: float
    inner_face_z: float
    predicted_endpoint_front_z: float
    predicted_event_duration_s: float
    same_cell_seed: ci.CutStepSeed
    face_event_controls: fi.FaceEventControls
    provenance: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(
            self.source_step,
            (ci.CutIntegratorStep, di.FaceDepartureIntegratorStep),
        ):
            raise TypeError(
                "accepted-history face routing needs an accepted same-cell or "
                "face-departure step"
            )
        if not isinstance(self.before, ci.CutIntegratorState):
            raise TypeError("accepted-history face routing needs a strict cut state")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("accepted-history face routing needs provenance")
        expected = _accepted_history_face_event_payload(
            self.source_step,
            self.requested_dt_s,
        )
        actual = (
            self.before,
            self.source_dt_s,
            self.source_start_front_z,
            self.source_end_front_z,
            self.inferred_front_velocity_z_s,
            self.inner_face_z,
            self.predicted_endpoint_front_z,
            self.predicted_event_duration_s,
            self.same_cell_seed,
            self.face_event_controls,
        )
        if actual != expected:
            raise ValueError(
                "accepted-history face-routing payload was changed after its "
                "deterministic reconstruction"
            )


def bracket_inner_face_from_accepted_history(
    source_step: AcceptedFaceHistoryStep,
    requested_dt_s: float,
    *,
    provenance: str = (
        "accepted trajectory tangent; no rejected-root, clamp, or manufactured "
        "event-time authority"
    ),
) -> AcceptedHistoryFaceEventBracket | None:
    """Return one deterministic predicted crossing, or ``None`` if unbracketed.

    A non-crossing accepted tangent is ordinary and returns ``None``.  Invalid,
    disconnected, non-conservative, or non-draining history fails closed.
    """

    try:
        payload = _accepted_history_face_event_payload(
            source_step,
            requested_dt_s,
        )
    except _NoAcceptedHistoryFaceBracket:
        return None
    (
        before,
        source_dt_s,
        source_start_front_z,
        source_end_front_z,
        velocity,
        inner_face_z,
        predicted_endpoint,
        predicted_event_duration,
        same_cell_seed,
        event_controls,
    ) = payload
    return AcceptedHistoryFaceEventBracket(
        source_step=source_step,
        before=before,
        requested_dt_s=requested_dt_s,
        source_dt_s=source_dt_s,
        source_start_front_z=source_start_front_z,
        source_end_front_z=source_end_front_z,
        inferred_front_velocity_z_s=velocity,
        inner_face_z=inner_face_z,
        predicted_endpoint_front_z=predicted_endpoint,
        predicted_event_duration_s=predicted_event_duration,
        same_cell_seed=same_cell_seed,
        face_event_controls=event_controls,
        provenance=provenance,
    )


def advance_bracketed_interior_face_macrostep(
    bracket: AcceptedHistoryFaceEventBracket,
    boundary: ct.PoreBoundary,
    *,
    boundary_discontinuity_times_s: Sequence[float] = (),
) -> "CutEventMacrostep":
    """Atomically try same-cell then solve one accepted-history face route."""

    if not isinstance(bracket, AcceptedHistoryFaceEventBracket):
        raise TypeError("bracketed face macrostep needs accepted-history evidence")
    # Reconstruct at use time as well.  This catches low-level mutation or a
    # forged instance that bypassed the frozen dataclass initializer.
    expected = _accepted_history_face_event_payload(
        bracket.source_step,
        bracket.requested_dt_s,
    )
    actual = (
        bracket.before,
        bracket.source_dt_s,
        bracket.source_start_front_z,
        bracket.source_end_front_z,
        bracket.inferred_front_velocity_z_s,
        bracket.inner_face_z,
        bracket.predicted_endpoint_front_z,
        bracket.predicted_event_duration_s,
        bracket.same_cell_seed,
        bracket.face_event_controls,
    )
    if actual != expected:
        raise ValueError("accepted-history face bracket failed use-time reconstruction")
    return advance_one_interior_face_macrostep(
        bracket.before,
        bracket.requested_dt_s,
        boundary,
        bracket.same_cell_seed,
        face_event_controls=bracket.face_event_controls,
        face_event_seed_factory=(
            lambda: face_event_seed_from_accepted_history_bracket(bracket)
        ),
        boundary_discontinuity_times_s=boundary_discontinuity_times_s,
    )


class _NoAcceptedHistoryFaceBracket(RuntimeError):
    """The accepted tangent does not cross the adjacent inner face."""


class AcceptedHistoryFaceEndpointUncertified(RuntimeError):
    """A tangent hit at the macro endpoint has no rigorous ordering proof."""


def _accepted_history_face_event_payload(
    source_step: AcceptedFaceHistoryStep,
    requested_dt_s: float,
) -> tuple[
    ci.CutIntegratorState,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    ci.CutStepSeed,
    fi.FaceEventControls,
]:
    """Reconstruct the complete physics-neutral routing payload."""

    if not isinstance(
        source_step,
        (ci.CutIntegratorStep, di.FaceDepartureIntegratorStep),
    ):
        raise TypeError(
            "face-routing history must be an accepted same-cell or departure step"
        )
    if not math.isfinite(requested_dt_s) or requested_dt_s <= 0.0:
        raise ValueError("face-routing requested duration must be positive and finite")
    if not getattr(source_step.ledger, "accepted", False):
        raise ValueError("face-routing history must have an accepted ledger")
    before = source_step.after
    if not isinstance(before, ci.CutIntegratorState):
        raise TypeError("face-routing history must end in one strict cut state")
    if source_step.ledger.maximum_step_ledger_residual > (
        before.controls.ledger_tolerance
    ) or source_step.ledger.maximum_cumulative_ledger_residual > (
        before.controls.ledger_tolerance
    ):
        raise ValueError("face-routing history exceeds the conservation contract")
    if source_step.ledger.maximum_scaled_residual > (
        before.controls.nonlinear_residual_tolerance
    ):
        raise ValueError("face-routing history exceeds the nonlinear contract")
    current_duals = (
        before.transport.effective_wet_retained_water_capacity_duals_over_rt
    )
    if any(value != 0.0 for value in current_duals):
        raise ValueError(
            "face-routing exact-arrival/departure charts do not yet carry the "
            "wet retained-cap graph dual; active history fails closed"
        )

    if isinstance(source_step, ci.CutIntegratorStep):
        if source_step.before.accepted_steps + 1 != before.accepted_steps:
            raise ValueError("same-cell face-routing history lost its accepted counter")
        source_dt_s = source_step.ledger.dt_s
        source_start_front_z = (
            source_step.before.transport.geometry.front.z
        )
    else:
        if source_step.before.accepted_steps + 1 != before.accepted_steps:
            raise ValueError("departure face-routing history lost its accepted counter")
        if source_step.before.geometry.master_grid is not (
            before.transport.geometry.master_grid
        ):
            raise ValueError("departure face-routing history changed its master grid")
        source_dt_s = source_step.ledger.dt_s
        source_start_front_z = source_step.before.geometry.front.z
    if not math.isfinite(source_dt_s) or source_dt_s <= 0.0:
        raise ValueError("face-routing history needs a positive accepted duration")

    source_end_front_z = before.transport.geometry.front.z
    velocity = math.fsum((source_end_front_z, -source_start_front_z)) / source_dt_s
    if not math.isfinite(velocity) or velocity >= 0.0:
        raise ValueError("face-routing history must be strict accepted primary drainage")
    cut_index = before.transport.layout.cut_cell_index
    if cut_index is None or cut_index <= 0:
        raise ValueError("accepted-history routing needs one interior cut face")
    grid = before.transport.geometry.master_grid
    inner_face_z = (grid.faces[cut_index] / grid.R) ** 3
    if not inner_face_z < source_end_front_z:
        raise ValueError("current strict cut state is not above its inner face")
    predicted_endpoint = math.fsum(
        (source_end_front_z, velocity * requested_dt_s)
    )
    if predicted_endpoint == inner_face_z:
        raise AcceptedHistoryFaceEndpointUncertified(
            "accepted front tangent lands exactly on the macro endpoint face, "
            "but production has no rigorous zero-width event-time certificate; "
            "refine the requested interval instead of falling through"
        )
    if not predicted_endpoint < inner_face_z:
        raise _NoAcceptedHistoryFaceBracket
    predicted_event_duration = (
        math.fsum((inner_face_z, -source_end_front_z)) / velocity
    )
    if not 0.0 < predicted_event_duration < requested_dt_s:
        raise RuntimeError("accepted-history face kinematics lost strict ordering")

    flux_scales = tuple(
        max(
            abs(value),
            ACCEPTED_HISTORY_STEFAN_SCALE_FLOOR_MOL_M2_S,
        )
        for value in before.last_total_stefan_fluxes_mol_m2_s
    )
    same_cell_seed = _accepted_history_same_cell_seed(
        source_step,
        before,
        inner_face_z=inner_face_z,
        inferred_front_velocity_z_s=velocity,
        predicted_event_duration_s=predicted_event_duration,
        flux_scales_mol_m2_s=flux_scales,
    )
    face_event_controls = fi.FaceEventControls((0.0, requested_dt_s))
    return (
        before,
        source_dt_s,
        source_start_front_z,
        source_end_front_z,
        velocity,
        inner_face_z,
        predicted_endpoint,
        predicted_event_duration,
        same_cell_seed,
        face_event_controls,
    )


def _accepted_history_same_cell_seed(
    source_step: AcceptedFaceHistoryStep,
    before: ci.CutIntegratorState,
    *,
    inner_face_z: float,
    inferred_front_velocity_z_s: float,
    predicted_event_duration_s: float,
    flux_scales_mol_m2_s: tuple[float, ...],
) -> ci.CutStepSeed:
    """Return one history-continuous strict-interior same-cell seed.

    The accepted tangent is backtracked in time, never componentwise clipped.
    This construction performs no nonlinear solve and consumes no rejected-root
    state.  Its first trial is halfway from the current front to the predicted
    face; further exact dyadic backtracking is used only when another primitive
    leaves its declared seed chart.
    """

    last_error: Exception | None = None
    for backtrack_exponent in range(48):
        seed_duration_s = math.ldexp(
            predicted_event_duration_s,
            -(backtrack_exponent + 1),
        )
        try:
            candidate = _accepted_history_same_cell_candidate(
                source_step,
                before,
                seed_duration_s=seed_duration_s,
                inferred_front_velocity_z_s=inferred_front_velocity_z_s,
            )
            if not (
                inner_face_z
                < candidate.front_z
                < before.transport.geometry.front.z
            ):
                raise ValueError(
                    "accepted-history same-cell seed is not strictly inside "
                    "the current material-cell chart"
                )
            seed = ci.CutStepSeed(
                candidate,
                flux_scales_mol_m2_s,
                (
                    "accepted-history temporal-continuation same-cell seed; "
                    f"half-face-horizon dyadic backtrack exponent "
                    f"{backtrack_exponent}; exact physical increments, no "
                    "clipping and no rejected-root authority"
                ),
            )
            ci._validate_seed(before, seed)  # noqa: SLF001
        except ValueError as exc:
            last_error = exc
            continue
        return seed
    raise ValueError(
        "accepted-history temporal-continuation seed could not enter the "
        "strict same-cell chart after deterministic dyadic backtracking"
    ) from last_error


def _accepted_history_same_cell_candidate(
    source_step: AcceptedFaceHistoryStep,
    before: ci.CutIntegratorState,
    *,
    seed_duration_s: float,
    inferred_front_velocity_z_s: float,
) -> cut.CutTransportUnknowns:
    """Evaluate the accepted differential tangent at one positive seed time."""

    current = before.transport
    if isinstance(source_step, ci.CutIntegratorStep):
        old = source_step.before.transport
        source_dt_s = source_step.ledger.dt_s
        old_duals = old.effective_wet_retained_water_capacity_duals_over_rt
        current_duals = (
            current.effective_wet_retained_water_capacity_duals_over_rt
        )
        if any(value != 0.0 for value in (*old_duals, *current_duals)):
            raise ValueError(
                "accepted-history same-cell tangent cannot backtrack an active "
                "wet retained-cap dual"
            )

        def extrapolate(
            old_values: Sequence[float],
            new_values: Sequence[float],
        ) -> tuple[float, ...]:
            ratio = seed_duration_s / source_dt_s
            return tuple(
                math.fsum((new_value, ratio * (new_value - old_value)))
                for old_value, new_value in zip(old_values, new_values, strict=True)
            )

        def extrapolate_wet_graph() -> tuple[float, ...]:
            if before.controls.wet_water_bounds[1] != (
                current.config.wet.luikov.W_cap
            ):
                return extrapolate(
                    old.wet_retained_water_loadings,
                    current.wet_retained_water_loadings,
                )
            wet_chart = wrc.WetRetainedCapSemismoothChart(
                before.controls.wet_water_bounds[0],
                current.config.wet.luikov,
            )
            ratio = seed_duration_s / source_dt_s
            result: list[float] = []
            for old_loading, new_loading in zip(
                old.wet_retained_water_loadings,
                current.wet_retained_water_loadings,
                strict=True,
            ):
                old_coordinate = wet_chart.encode(old_loading, 0.0)
                new_coordinate = wet_chart.encode(new_loading, 0.0)
                predicted_coordinate = math.fsum(
                    (
                        new_coordinate,
                        ratio * (new_coordinate - old_coordinate),
                    )
                )
                point = wet_chart.decode(predicted_coordinate)
                if point.capacity_dual_over_rt != 0.0:
                    raise ValueError(
                        "accepted-history wet graph seed reaches retained-cap "
                        "contact"
                    )
                result.append(point.retained_water_loading)
            return tuple(result)

        wet_temperature = extrapolate(
            old.wet_temperatures_k,
            current.wet_temperatures_k,
        )
        wet_water = extrapolate_wet_graph()
        dry_temperature = extrapolate(
            old.dry_temperatures_k,
            current.dry_temperatures_k,
        )
        dry_y = extrapolate(
            old.dry_y_hexane,
            current.dry_y_hexane,
        )
        interface_temperature = math.fsum(
            (
                before.last_interface_temperature_k,
                seed_duration_s
                / source_dt_s
                * (
                    before.last_interface_temperature_k
                    - source_step.before.last_interface_temperature_k
                ),
            )
        )
    else:
        # The exact-face source has one fewer dry primitive than the strict
        # departure result.  Current accepted rank-safe primitives are the
        # unique zero-state-increment continuation seed; only the accepted
        # front tangent has a rank-preserving differential history.
        wet_temperature = current.wet_temperatures_k
        wet_water = current.wet_retained_water_loadings
        dry_temperature = current.dry_temperatures_k
        dry_y = current.dry_y_hexane
        interface_temperature = before.last_interface_temperature_k

    return cut.CutTransportUnknowns(
        wet_temperatures_k=tuple(wet_temperature),
        wet_retained_water_loadings=tuple(wet_water),
        dry_temperatures_k=tuple(dry_temperature),
        dry_y_hexane=tuple(dry_y),
        dry_total_stefan_fluxes_mol_m2_s=(
            before.last_total_stefan_fluxes_mol_m2_s
        ),
        front_z=math.fsum(
            (
                current.geometry.front.z,
                inferred_front_velocity_z_s * seed_duration_s,
            )
        ),
        interface_temperature_k=interface_temperature,
        wet_retained_water_capacity_duals_over_rt=(0.0,) * len(wet_water),
    )


def face_event_seed_from_accepted_history_bracket(
    bracket: AcceptedHistoryFaceEventBracket,
) -> fi.FaceEventSeed:
    """Project the current accepted state onto exact-face rank for seeding only.

    This helper is intentionally separate from bracket construction.  The
    production wrapper calls it exactly once, and only after the unchanged
    same-cell equations reject with exact rollback.  The accepted tangent
    supplies event time only; no wet, dry, or interface field is extrapolated.
    """

    if not isinstance(bracket, AcceptedHistoryFaceEventBracket):
        raise TypeError("face-event seed construction needs a frozen bracket")
    current = bracket.before.transport
    if any(
        value != 0.0
        for value in current.effective_wet_retained_water_capacity_duals_over_rt
    ):
        raise ValueError(
            "current accepted face seed cannot discard a wet retained-cap dual"
        )
    candidate = face.FaceArrivalUnknowns(
        wet_temperatures_k=current.wet_temperatures_k[:-1],
        wet_retained_water_loadings=(
            current.wet_retained_water_loadings[:-1]
        ),
        dry_temperatures_k=current.dry_temperatures_k,
        dry_y_hexane=current.dry_y_hexane,
        dry_total_stefan_fluxes_mol_m2_s=(
            bracket.before.last_total_stefan_fluxes_mol_m2_s
        ),
        event_time_s=bracket.predicted_event_duration_s,
        interface_temperature_k=(bracket.before.last_interface_temperature_k),
    )
    layout = face.layout_for_arrival(current)
    if len(candidate.vector()) != layout.unknown_count:
        raise ValueError("current-state face seed lost exact arrival rank")
    return fi.FaceEventSeed(
        candidate,
        bracket.same_cell_seed.stefan_flux_scales_mol_m2_s,
        (
            "current accepted state projected onto exact-face rank after "
            "same-cell rejection; accepted tangent supplies event time only; "
            "no field extrapolation, clipping, or rejected-root authority"
        ),
    )


class CutMacrostepBranch(str, Enum):
    """Accepted topology path through one requested duration."""

    SAME_CELL = "same_cell"
    FACE_ARRIVAL_ENDPOINT = "face_arrival_endpoint"
    FACE_ARRIVAL_AND_DEPARTURE = "face_arrival_and_departure"


class EventOrderingMode(str, Enum):
    """Authority used to place an exact-face event inside a macrostep."""

    RIGOROUS_INTERVAL_CERTIFICATE = "rigorous_interval_certificate"
    OPEN_INTERVAL_NUMERICAL_ORDERING = "open_interval_numerical_ordering"


class CutMacrostepFailureCode(str, Enum):
    """Fail-closed reason owned by the atomic orchestration layer."""

    BOUNDARY_DISCONTINUITY = "boundary_discontinuity_inside_macrostep"
    SAME_CELL_FAILED = "same_cell_solve_failed"
    EXTINCTION_REQUIRED = "cell_zero_extinction_localizer_required"
    MISSING_FACE_EVENT_INPUT = "missing_face_event_input"
    FACE_ARRIVAL_FAILED = "exact_face_arrival_failed"
    MISSING_EVENT_TIME_CERTIFICATE = "missing_event_time_certificate"
    INVALID_EVENT_TIME_CERTIFICATE = "invalid_event_time_certificate"
    EVENT_TIME_ORDERING_UNCERTIFIED = "event_time_ordering_uncertified"
    FACE_EVENT_AFTER_MACROSTEP = "face_event_not_strictly_before_macrostep_end"
    FACE_TANGENT_FAILED = "exact_face_tangent_failed"
    ADDITIONAL_FACE_EVENT_REQUIRED = "additional_face_event_required"
    EXTINCTION_REMAINDER_REQUIRED = "extinction_inside_departure_remainder"
    FACE_DEPARTURE_FAILED = "face_departure_remainder_failed"
    RETAINED_WATER_APPLICABILITY_GUARD_FAILED = (
        "wet_retained_water_applicability_guard_failed"
    )
    TOPOLOGY_HANDOFF_FAILED = "topology_handoff_failed"
    MACRO_LEDGER_FAILED = "macrostep_ledger_failed"
    INTERNAL_FAILURE = "unclassified_atomic_macrostep_failure"


class CutEventMacrostepError(RuntimeError):
    """Rejected atomic macrostep carrying the exact original rollback state."""

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState,
        *,
        code: CutMacrostepFailureCode,
        same_cell_error: ci.CutIntegratorStepError | None = None,
        provisional_event_duration_s: float | None = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.code = code
        self.same_cell_error = same_cell_error
        self.provisional_event_duration_s = provisional_event_duration_s


@dataclass(frozen=True)
class FaceEventTimeCertificate:
    """Independent interval certificate for one computed arrival duration.

    ``absolute_error_bound_s`` is an error bound, not an event deadband.  The
    orchestrator accepts a post-event remainder only when the entire certified
    interval lies strictly before the requested endpoint.  A zero bound is
    reserved for exact/manufactured authorities.
    """

    before_time_s: float
    event_duration_s: float
    absolute_error_bound_s: float
    arrival_face_index: int
    label: str
    mode: EventOrderingMode = field(
        default=EventOrderingMode.RIGOROUS_INTERVAL_CERTIFICATE,
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.before_time_s,
            self.event_duration_s,
            self.absolute_error_bound_s,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("event-time certificate values must be finite")
        if self.before_time_s < 0.0:
            raise ValueError("certified before-time cannot be negative")
        if self.event_duration_s <= 0.0:
            raise ValueError("certified event duration must be positive")
        if self.absolute_error_bound_s < 0.0:
            raise ValueError("event-time error bound cannot be negative")
        if (
            isinstance(self.arrival_face_index, bool)
            or not isinstance(self.arrival_face_index, int)
            or self.arrival_face_index <= 0
        ):
            raise ValueError("an interior arrival certificate needs face index > 0")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("event-time certificate needs provenance")

    @property
    def lower_duration_s(self) -> float:
        return math.fsum((self.event_duration_s, -self.absolute_error_bound_s))

    @property
    def upper_duration_s(self) -> float:
        return math.fsum((self.event_duration_s, self.absolute_error_bound_s))


FaceEventTimeCertifier = Callable[[fi.FaceEventStep], FaceEventTimeCertificate]
FaceEventSeedFactory = Callable[[], fi.FaceEventSeed]


@dataclass(frozen=True)
class OpenIntervalNumericalOrderingEvidence:
    r"""Nonqualifying evidence that a computed event is numerically interior.

    ``newton_guard_*`` is the symmetric hull generated by the magnitude of one
    independently assembled physical Newton time correction.  It is a
    representability/near-overlap guard only—not an error interval.  One
    floating-point neighbour must remain between that hull and each endpoint.
    """

    before_time_s: float
    requested_dt_s: float
    event_duration_s: float
    newton_corrected_event_duration_s: float
    linearized_event_time_correction_s: float
    newton_guard_lower_duration_s: float
    newton_guard_upper_duration_s: float
    arrival_face_index: int
    arrival_radius_change_m: float
    interface_swept_volume_rate_m3_s: float
    linearization_audit: fi.FaceEventTimeLinearizationAudit
    label: str
    event_time_absolute_error_bound_s: None = field(default=None, init=False)
    rigorous_interval_error_bound_available: bool = field(default=False, init=False)
    mode: EventOrderingMode = field(
        default=EventOrderingMode.OPEN_INTERVAL_NUMERICAL_ORDERING,
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.before_time_s,
            self.requested_dt_s,
            self.event_duration_s,
            self.newton_corrected_event_duration_s,
            self.linearized_event_time_correction_s,
            self.newton_guard_lower_duration_s,
            self.newton_guard_upper_duration_s,
            self.arrival_radius_change_m,
            self.interface_swept_volume_rate_m3_s,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("open-interval ordering evidence must be finite")
        if self.before_time_s < 0.0 or self.requested_dt_s <= 0.0:
            raise ValueError("open-interval ordering needs positive forward time")
        if not isinstance(
            self.linearization_audit,
            fi.FaceEventTimeLinearizationAudit,
        ):
            raise TypeError("open-interval ordering needs the physical-Jacobian audit")
        if not self.linearization_audit.jacobian_certificate.full_rank:
            raise ValueError("open-interval ordering lost full nonlinear rank")
        if self.linearization_audit.ordering_certified:
            raise ValueError(
                "the local linearization audit must not claim rigorous ordering"
            )
        if self.linearized_event_time_correction_s != (
            self.linearization_audit.linearized_event_time_correction_s
        ):
            raise ValueError("ordering evidence changed the audited Newton correction")
        if self.arrival_radius_change_m >= 0.0 or (
            self.interface_swept_volume_rate_m3_s >= 0.0
        ):
            raise ValueError("open-interval event is not primary drainage")
        if (
            isinstance(self.arrival_face_index, bool)
            or not isinstance(self.arrival_face_index, int)
            or self.arrival_face_index <= 0
        ):
            raise ValueError("open-interval ordering needs an interior face")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("open-interval ordering evidence needs provenance")

        corrected = math.fsum(
            (self.event_duration_s, self.linearized_event_time_correction_s)
        )
        radius = abs(self.linearized_event_time_correction_s)
        expected_lower = math.fsum((self.event_duration_s, -radius))
        expected_upper = math.fsum((self.event_duration_s, radius))
        if corrected != self.newton_corrected_event_duration_s or (
            expected_lower != self.newton_guard_lower_duration_s
            or expected_upper != self.newton_guard_upper_duration_s
        ):
            raise ValueError("open-interval Newton time accounting is inconsistent")
        if not (
            0.0 < self.event_duration_s < self.requested_dt_s
            and 0.0
            < self.newton_corrected_event_duration_s
            < self.requested_dt_s
        ):
            raise ValueError(
                "candidate and Newton-corrected event times must be strictly interior"
            )
        representable_lower = math.nextafter(
            self.newton_guard_lower_duration_s,
            -math.inf,
        )
        representable_upper = math.nextafter(
            self.newton_guard_upper_duration_s,
            math.inf,
        )
        if not 0.0 < representable_lower or not (
            representable_upper < self.requested_dt_s
        ):
            raise ValueError(
                "Newton correction guard is not representably separated from an "
                "open-interval endpoint"
            )

    @property
    def condition_number_2(self) -> float:
        return self.linearization_audit.jacobian_certificate.condition_number_2

    @property
    def maximum_scaled_residual(self) -> float:
        return max(abs(value) for value in self.linearization_audit.scaled_residuals)

    @property
    def newton_corrected_maximum_scaled_residual(self) -> float:
        return self.linearization_audit.newton_corrected_maximum_scaled_residual


EventOrderingEvidence = (
    FaceEventTimeCertificate | OpenIntervalNumericalOrderingEvidence
)


@dataclass(frozen=True)
class CutEventMacrostepLedger:
    """Independent stage telescope and atomic-commit diagnostics."""

    branch: CutMacrostepBranch
    requested_dt_s: float
    event_duration_s: float
    event_time_error_bound_s: float | None
    departure_remainder_s: float
    duration_closure_error_s: float
    absolute_time_closure_error_s: float
    accepted_balance_interval_count: int
    accepted_nonlinear_evaluations: int
    total_attempted_nonlinear_evaluations: int
    rejected_same_cell_nonlinear_evaluations: int
    uncommitted_continuation_root_count: int
    uncommitted_continuation_nonlinear_evaluations: int
    event_ordering_mode: EventOrderingMode | None
    event_ordering: tuple[str, ...]
    maximum_stage_scaled_residual: float
    maximum_stage_step_ledger_residual: float
    maximum_final_cumulative_ledger_residual: float
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
class CutEventMacrostep:
    """One atomically accepted same-cell or one-face composite macrostep."""

    before: ci.CutIntegratorState
    after: ci.CutIntegratorState | fi.FaceArrivalCommittedState
    boundary: ct.PoreBoundary
    same_cell_step: ci.CutIntegratorStep | None
    face_event_step: fi.FaceEventStep | None
    face_departure_step: di.FaceDepartureIntegratorStep | None
    ledger: CutEventMacrostepLedger
    event_ordering_evidence: EventOrderingEvidence | None = None
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        stages = (
            self.same_cell_step is not None,
            self.face_event_step is not None,
            self.face_departure_step is not None,
        )
        expected = {
            CutMacrostepBranch.SAME_CELL: (True, False, False),
            CutMacrostepBranch.FACE_ARRIVAL_ENDPOINT: (False, True, False),
            CutMacrostepBranch.FACE_ARRIVAL_AND_DEPARTURE: (False, True, True),
        }[self.ledger.branch]
        if stages != expected:
            raise ValueError("macrostep branch and accepted stage payloads disagree")
        if not self.ledger.accepted:
            raise ValueError("an accepted macrostep requires an accepted ledger")
        if self.ledger.branch is CutMacrostepBranch.SAME_CELL:
            if self.event_ordering_evidence is not None:
                raise ValueError("same-cell macrostep cannot carry event evidence")
        elif self.event_ordering_evidence is None:
            raise ValueError("an event macrostep must retain its ordering evidence")
        elif self.ledger.event_ordering_mode is not (
            self.event_ordering_evidence.mode
        ):
            raise ValueError("event evidence and macro ledger use different authorities")
        elif (
            self.face_event_step is None
            or self.event_ordering_evidence.before_time_s
            != self.before.transport.time_s
            or self.event_ordering_evidence.event_duration_s
            != self.face_event_step.ledger.event_time_s
            or self.event_ordering_evidence.arrival_face_index
            != self.face_event_step.after.arrival_face_index
        ):
            raise ValueError("event evidence does not identify the retained face solve")
        if isinstance(
            self.event_ordering_evidence,
            OpenIntervalNumericalOrderingEvidence,
        ):
            if self.ledger.event_time_error_bound_s is not None:
                raise ValueError("numerical ordering must not populate an error bound")
            if self.event_ordering_evidence.requested_dt_s != (
                self.ledger.requested_dt_s
            ):
                raise ValueError("numerical ordering used a different requested interval")


class _AtomicRejection(RuntimeError):
    def __init__(self, code: CutMacrostepFailureCode, message: str) -> None:
        super().__init__(message)
        self.code = code


def _certify_same_cell_forecast_disproof(
    before: ci.CutIntegratorState,
    step: ci.CutIntegratorStep,
) -> None:
    """Certify that a corrected same-cell root invalidates only the forecast."""

    if not isinstance(step, ci.CutIntegratorStep):
        raise _AtomicRejection(
            CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
            "same-cell solver returned the wrong accepted payload type",
        )
    if step.before is not before:
        raise _AtomicRejection(
            CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
            "same-cell root lost the exact immutable input state",
        )
    old = before.transport
    new = step.after.transport
    cut_index = old.layout.cut_cell_index
    if cut_index is None or new.layout.cut_cell_index != cut_index:
        raise _AtomicRejection(
            CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
            "same-cell root changed the fixed material-cell chart",
        )
    grid = old.geometry.master_grid
    inner_face_z = (grid.faces[cut_index] / grid.R) ** 3
    old_front_z = old.geometry.front.z
    new_front_z = new.geometry.front.z
    if not inner_face_z < new_front_z < old_front_z:
        raise _AtomicRejection(
            CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
            "same-cell root did not retain a strict monotone primary-drainage "
            "sweep above the adjacent inner material face",
        )
    ledger = step.ledger
    if not ledger.accepted or (
        ledger.maximum_scaled_residual
        > before.controls.nonlinear_residual_tolerance
    ):
        raise _AtomicRejection(
            CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "same-cell forecast disproof lacks an accepted nonlinear root",
        )
    if (
        ledger.maximum_step_ledger_residual > before.controls.ledger_tolerance
        or ledger.maximum_cumulative_ledger_residual
        > before.controls.ledger_tolerance
    ):
        raise _AtomicRejection(
            CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "same-cell forecast disproof exceeds the conservation contract",
        )
    ci.certify_enabled_wet_retained_water_state_applicability(
        before.controls,
        new,
    )


def advance_one_interior_face_macrostep(
    before: ci.CutIntegratorState,
    dt_s: float,
    boundary: ct.PoreBoundary,
    same_cell_seed: ci.CutStepSeed,
    *,
    face_event_controls: fi.FaceEventControls | None = None,
    face_event_seed: fi.FaceEventSeed | None = None,
    face_event_seed_factory: FaceEventSeedFactory | None = None,
    event_time_certifier: FaceEventTimeCertifier | None = None,
    boundary_discontinuity_times_s: Sequence[float] = (),
) -> CutEventMacrostep:
    """Advance one requested duration with at most one interior-face event.

    A failed same-cell solve is diagnostic work only.  If the exact-arrival
    route is used, the arrival and requested-remainder roots are provisional
    until the final macro ledger passes.  No accepted adaptive subdivision is
    performed here.
    """

    if not isinstance(before, ci.CutIntegratorState):
        raise TypeError("atomic cut macrostep requires CutIntegratorState")
    same_cell_error: ci.CutIntegratorStepError | None = None
    provisional_event_duration: float | None = None
    try:
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("atomic cut macrostep requires one supported pore boundary")
        if not math.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("atomic cut macrostep duration must be positive and finite")
        if face_event_seed is not None and face_event_seed_factory is not None:
            raise ValueError(
                "supply either a materialized face-event seed or one lazy "
                "factory, never both"
            )
        if face_event_seed_factory is not None and not callable(
            face_event_seed_factory
        ):
            raise TypeError("face-event seed factory must be callable")
        ci.certify_enabled_wet_retained_water_state_applicability(
            before.controls,
            before.transport,
        )
        _reject_internal_boundary_discontinuity(
            before.transport.time_s,
            dt_s,
            boundary_discontinuity_times_s,
        )

        try:
            same_cell = ci.advance_same_cell_backward_euler(
                before,
                dt_s,
                boundary,
                same_cell_seed,
            )
        except ci.CutIntegratorStepError as exc:
            same_cell_error = exc
            if exc.retained_water_applicability_guard_failed:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.RETAINED_WATER_APPLICABILITY_GUARD_FAILED,
                    "same-cell root left the selected smooth retained-water scope",
                ) from exc
        else:
            if same_cell.before is not before:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                    "same-cell solver did not preserve the exact before object",
                )
            _certify_same_cell_forecast_disproof(before, same_cell)
            ledger = _compose_macro_ledger(
                before,
                same_cell.after,
                dt_s,
                CutMacrostepBranch.SAME_CELL,
                (same_cell.ledger,),
                event_duration_s=0.0,
                event_time_error_bound_s=0.0,
                departure_remainder_s=0.0,
                failed_same_cell_evaluations=0,
                uncommitted_continuation_roots=(),
                event_ordering_mode=None,
                event_ordering=("same_cell",),
            )
            ci.certify_enabled_wet_retained_water_state_applicability(
                before.controls,
                same_cell.after.transport,
            )
            return CutEventMacrostep(
                before,
                same_cell.after,
                boundary,
                same_cell,
                None,
                None,
                ledger,
            )

        if same_cell_error.rollback_state is not before:
            raise _AtomicRejection(
                CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                "same-cell rejection did not preserve the exact rollback object",
            )
        cut_index = before.transport.layout.cut_cell_index
        if cut_index == 0:
            raise _AtomicRejection(
                CutMacrostepFailureCode.EXTINCTION_REQUIRED,
                "same-cell cell-zero failure requires the separately certified "
                "extinction-time limit and exact-dry projection",
            )
        if cut_index is None or cut_index < 0:
            raise _AtomicRejection(
                CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                "atomic interior-face routing requires one strict cut cell",
            )
        if face_event_controls is None or (
            face_event_seed is None and face_event_seed_factory is None
        ):
            raise _AtomicRejection(
                CutMacrostepFailureCode.MISSING_FACE_EVENT_INPUT,
                "same-cell failure needs explicit face-event controls and seed",
            )
        if face_event_seed is None:
            assert face_event_seed_factory is not None
            try:
                face_event_seed = face_event_seed_factory()
            except Exception as exc:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.FACE_ARRIVAL_FAILED,
                    f"lazy current-state face seed did not certify: {exc}",
                ) from exc
            if not isinstance(face_event_seed, fi.FaceEventSeed):
                raise _AtomicRejection(
                    CutMacrostepFailureCode.FACE_ARRIVAL_FAILED,
                    "lazy face-event seed factory returned the wrong payload type",
                )

        try:
            face_event = fi.solve_inner_face_arrival(
                before,
                face_event_controls,
                boundary,
                face_event_seed,
            )
        except fi.FaceEventIntegratorStepError as exc:
            if exc.rollback_state is not before:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                    "face-arrival rejection lost the original rollback object",
                ) from exc
            if exc.retained_water_applicability_guard_failed:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.RETAINED_WATER_APPLICABILITY_GUARD_FAILED,
                    "face-arrival root left the selected smooth retained-water scope",
                ) from exc
            raise _AtomicRejection(
                CutMacrostepFailureCode.FACE_ARRIVAL_FAILED,
                f"exact inner-face arrival did not certify: {exc}",
            ) from exc
        provisional_event_duration = face_event.ledger.event_time_s
        if face_event.before is not before or (
            face_event.after.arrival_face_index != cut_index
        ):
            raise _AtomicRejection(
                CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                "face-arrival payload does not match the original cut face",
            )

        ordering_evidence = _event_time_ordering_evidence(
            before,
            dt_s,
            face_event,
            face_event_controls,
            event_time_certifier,
        )
        tau = face_event.ledger.event_time_s
        if isinstance(ordering_evidence, FaceEventTimeCertificate):
            certificate = ordering_evidence
            if certificate.lower_duration_s <= 0.0:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.EVENT_TIME_ORDERING_UNCERTIFIED,
                    "event-time interval does not lie strictly after macrostep start",
                )

            # Exact endpoint acceptance needs an exact/manufactured zero-width
            # certificate.  Numerical open-interval evidence can never take
            # this branch.
            if certificate.absolute_error_bound_s == 0.0 and tau == dt_s:
                ledger = _compose_macro_ledger(
                    before,
                    face_event.after,
                    dt_s,
                    CutMacrostepBranch.FACE_ARRIVAL_ENDPOINT,
                    (face_event.ledger,),
                    event_duration_s=tau,
                    event_time_error_bound_s=0.0,
                    departure_remainder_s=0.0,
                    failed_same_cell_evaluations=(
                        same_cell_error.nonlinear_evaluations
                        if same_cell_error is not None
                        else 0
                    ),
                    uncommitted_continuation_roots=(),
                    event_ordering_mode=certificate.mode,
                    event_ordering=(
                        "same_cell_rejected",
                        "rigorous_interval_certificate",
                        "exact_face_arrival",
                    ),
                )
                ci.certify_enabled_wet_retained_water_state_applicability(
                    before.controls,
                    face_event.after,
                    capacity_duals_known_zero=True,
                )
                return CutEventMacrostep(
                    before,
                    face_event.after,
                    boundary,
                    None,
                    face_event,
                    None,
                    ledger,
                    certificate,
                )

            if certificate.lower_duration_s >= dt_s:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.FACE_EVENT_AFTER_MACROSTEP,
                    "certified face arrival is not strictly before macrostep end",
                )
            if certificate.upper_duration_s >= dt_s:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.EVENT_TIME_ORDERING_UNCERTIFIED,
                    "event-time interval overlaps macrostep end; refine rather than "
                    "collapse a physical remainder",
                )
            event_time_error_bound_s: float | None = (
                certificate.absolute_error_bound_s
            )
        else:
            event_time_error_bound_s = None
        remainder = math.fsum((dt_s, -tau))
        if remainder <= 0.0:
            raise _AtomicRejection(
                CutMacrostepFailureCode.EVENT_TIME_ORDERING_UNCERTIFIED,
                "pre-end face arrival produced no representable remainder",
            )

        try:
            tangent_candidate = tangent.candidate_from_committed_face(face_event.after)
            tangent_assembly = tangent.assemble_face_tangent(
                face_event.after,
                tangent_candidate,
            )
        except tangent.FaceTangentStepError as exc:
            raise _AtomicRejection(
                CutMacrostepFailureCode.FACE_TANGENT_FAILED,
                f"exact-face tangent did not certify: {exc}",
            ) from exc

        face_radius = face_event.after.geometry.front.radius_m
        inner_radius = face_event.after.geometry.master_grid.faces[cut_index - 1]
        predicted_radius = math.fsum(
            (
                face_radius,
                tangent_assembly.candidate.front_speed_m_s * remainder,
            )
        )
        if predicted_radius <= inner_radius:
            code = (
                CutMacrostepFailureCode.EXTINCTION_REMAINDER_REQUIRED
                if cut_index == 1
                else CutMacrostepFailureCode.ADDITIONAL_FACE_EVENT_REQUIRED
            )
            raise _AtomicRejection(
                code,
                "exact tangent predicts that the unconsumed remainder reaches "
                "another topology event; a face-to-face/extinction solver is required",
            )
        if predicted_radius >= face_radius:
            raise _AtomicRejection(
                CutMacrostepFailureCode.FACE_TANGENT_FAILED,
                "exact tangent is not a strict primary-drainage departure",
            )

        try:
            overlap = tangent.continue_to_positive_sweep(
                tangent_assembly,
                remainder,
                boundary,
            )
        except tangent.FaceTangentStepError as exc:
            raise _AtomicRejection(
                CutMacrostepFailureCode.FACE_TANGENT_FAILED,
                f"finite exact-tangent departure seed failed: {exc}",
            ) from exc
        departure_seed = _departure_seed_from_overlap(
            overlap,
            face_event_seed,
        )
        try:
            departure_step = di.solve_face_departure(
                face_event.after,
                before.controls,
                remainder,
                boundary,
                departure_seed,
            )
        except di.FaceDepartureIntegratorStepError as exc:
            if exc.rollback_state is not face_event.after:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                    "departure rejection lost the provisional exact-face rollback",
                ) from exc
            if exc.retained_water_applicability_guard_failed:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.RETAINED_WATER_APPLICABILITY_GUARD_FAILED,
                    "face-departure root left the selected smooth retained-water scope",
                ) from exc
            raise _AtomicRejection(
                CutMacrostepFailureCode.FACE_DEPARTURE_FAILED,
                f"exact-face departure remainder did not certify: {exc}",
            ) from exc

        if departure_step.before is not face_event.after or (
            departure_step.after.transport.layout.cut_cell_index != cut_index - 1
        ):
            raise _AtomicRejection(
                CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                "departure did not enter exactly the adjacent strict cut cell",
            )
        auxiliary = tuple(
            getattr(departure_step.ledger, "auxiliary_stages", ())
        )
        if any(getattr(stage, "committed", True) for stage in auxiliary):
            raise _AtomicRejection(
                CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
                "a departure continuation root was incorrectly marked committed",
            )
        ledger = _compose_macro_ledger(
            before,
            departure_step.after,
            dt_s,
            CutMacrostepBranch.FACE_ARRIVAL_AND_DEPARTURE,
            (face_event.ledger, departure_step.ledger),
            event_duration_s=tau,
            event_time_error_bound_s=event_time_error_bound_s,
            departure_remainder_s=remainder,
            failed_same_cell_evaluations=(
                same_cell_error.nonlinear_evaluations
                if same_cell_error is not None
                else 0
            ),
            uncommitted_continuation_roots=auxiliary,
            event_ordering_mode=ordering_evidence.mode,
            event_ordering=(
                "same_cell_rejected",
                ordering_evidence.mode.value,
                "exact_face_arrival",
                "exact_face_departure_remainder",
            ),
        )
        ci.certify_enabled_wet_retained_water_state_applicability(
            before.controls,
            departure_step.after.transport,
        )
        return CutEventMacrostep(
            before,
            departure_step.after,
            boundary,
            None,
            face_event,
            departure_step,
            ledger,
            ordering_evidence,
        )
    except CutEventMacrostepError:
        raise
    except _AtomicRejection as exc:
        raise CutEventMacrostepError(
            f"atomic one-face macrostep rejected with exact rollback: {exc}",
            before,
            code=exc.code,
            same_cell_error=same_cell_error,
            provisional_event_duration_s=provisional_event_duration,
        ) from exc
    except ci.WetRetainedWaterApplicabilityError as exc:
        raise CutEventMacrostepError(
            f"atomic one-face macrostep rejected with exact rollback: {exc}",
            before,
            code=(
                CutMacrostepFailureCode.RETAINED_WATER_APPLICABILITY_GUARD_FAILED
            ),
            same_cell_error=same_cell_error,
            provisional_event_duration_s=provisional_event_duration,
        ) from exc
    except Exception as exc:
        raise CutEventMacrostepError(
            f"atomic one-face macrostep rejected with exact rollback: {exc}",
            before,
            code=CutMacrostepFailureCode.INTERNAL_FAILURE,
            same_cell_error=same_cell_error,
            provisional_event_duration_s=provisional_event_duration,
        ) from exc


def _reject_internal_boundary_discontinuity(
    start_time_s: float,
    dt_s: float,
    discontinuity_times_s: Sequence[float],
) -> None:
    times = tuple(float(value) for value in discontinuity_times_s)
    if not all(math.isfinite(value) and value >= 0.0 for value in times):
        raise ValueError("boundary discontinuity times must be finite and non-negative")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ValueError("boundary discontinuity times must increase strictly")
    end_time = math.fsum((start_time_s, dt_s))
    inside = tuple(value for value in times if start_time_s < value < end_time)
    if inside:
        raise _AtomicRejection(
            CutMacrostepFailureCode.BOUNDARY_DISCONTINUITY,
            "one constant-boundary BE macrostep cannot straddle discontinuities "
            f"at {inside!r}",
        )


def _event_time_ordering_evidence(
    before: ci.CutIntegratorState,
    requested_dt_s: float,
    face_event: fi.FaceEventStep,
    controls: fi.FaceEventControls,
    certifier: FaceEventTimeCertifier | None,
) -> EventOrderingEvidence:
    if certifier is None:
        if isinstance(controls, fi.FaceEventControls) and (
            controls.event_time_bounds_s[0] == 0.0
        ):
            if controls.event_time_bounds_s[1] != requested_dt_s:
                raise _AtomicRejection(
                    CutMacrostepFailureCode.EVENT_TIME_ORDERING_UNCERTIFIED,
                    "open-interval numerical ordering requires the exact direct "
                    "event chart (0, requested_dt)",
                )
            return _open_interval_numerical_ordering(
                before,
                requested_dt_s,
                face_event,
            )
        raise _AtomicRejection(
            CutMacrostepFailureCode.MISSING_EVENT_TIME_CERTIFICATE,
            "event routing needs either an independent error certificate or the "
            "exact nonqualifying open chart (0, requested_dt)",
        )
    try:
        certificate = certifier(face_event)
    except Exception as exc:
        raise _AtomicRejection(
            CutMacrostepFailureCode.INVALID_EVENT_TIME_CERTIFICATE,
            f"event-time certifier failed: {exc}",
        ) from exc
    if not isinstance(certificate, FaceEventTimeCertificate):
        raise _AtomicRejection(
            CutMacrostepFailureCode.INVALID_EVENT_TIME_CERTIFICATE,
            "event-time certifier returned the wrong payload type",
        )
    if (
        certificate.before_time_s != before.transport.time_s
        or certificate.event_duration_s != face_event.ledger.event_time_s
        or certificate.arrival_face_index != face_event.after.arrival_face_index
    ):
        raise _AtomicRejection(
            CutMacrostepFailureCode.INVALID_EVENT_TIME_CERTIFICATE,
            "event-time certificate does not identify the solved before/event pair",
        )
    return certificate


def _open_interval_numerical_ordering(
    before: ci.CutIntegratorState,
    requested_dt_s: float,
    face_event: fi.FaceEventStep,
) -> OpenIntervalNumericalOrderingEvidence:
    """Build local numerical ordering evidence without inventing an error bound."""

    try:
        audit = fi.audit_face_event_time_linearization(face_event)
        tau = face_event.ledger.event_time_s
        correction = audit.linearized_event_time_correction_s
        corrected = math.fsum((tau, correction))
        correction_radius = abs(correction)
        guard_lower = math.fsum((tau, -correction_radius))
        guard_upper = math.fsum((tau, correction_radius))
        evidence = OpenIntervalNumericalOrderingEvidence(
            before_time_s=before.transport.time_s,
            requested_dt_s=requested_dt_s,
            event_duration_s=tau,
            newton_corrected_event_duration_s=corrected,
            linearized_event_time_correction_s=correction,
            newton_guard_lower_duration_s=guard_lower,
            newton_guard_upper_duration_s=guard_upper,
            arrival_face_index=face_event.after.arrival_face_index,
            arrival_radius_change_m=math.fsum(
                (
                    face_event.after.geometry.front.radius_m,
                    -before.transport.geometry.front.radius_m,
                )
            ),
            interface_swept_volume_rate_m3_s=(
                face_event.assembly.swept_geometry.interface_swept_volume_rate_m3_s
            ),
            linearization_audit=audit,
            label=(
                "direct one-sided exact-face solve on (0, requested_dt); "
                "physical-Newton interiority evidence only, not an error interval"
            ),
        )
        nonlinear_tolerance = before.controls.nonlinear_residual_tolerance
        if evidence.maximum_scaled_residual > nonlinear_tolerance:
            raise ValueError(
                "candidate or Newton-corrected residual exceeds the unchanged "
                "nonlinear acceptance contract"
            )
        if evidence.newton_corrected_maximum_scaled_residual > nonlinear_tolerance:
            # GT-PS-2-R5S7: the corrected-point check is meaningful only when
            # the Newton correction is resolvable.  The correction is solved
            # from the audit's finite-difference Jacobian, whose per-column
            # central probe steps the audit records; a correction smaller
            # than EVERY probe step that built the Jacobian is a displacement
            # the linear model cannot distinguish from the certified center,
            # so its re-evaluated residual measures the assembly's off-root
            # evaluation floor, not the root.  In that sub-resolution case
            # the center residual (checked above against the UNCHANGED
            # contract) and the rank/condition certificate govern, and the
            # branch is disclosed on the evidence label with the corrected
            # residual recorded but not certified.  No tolerance, contract,
            # audit quantity, or dataclass changes; the comparison uses only
            # quantities the audit already records.
            correction_resolvable = any(
                abs(scale * component) >= probe
                for scale, component, probe in zip(
                    audit.physical_variable_scales,
                    audit.dimensionless_newton_correction,
                    audit.physical_jacobian_perturbations,
                )
            )
            if correction_resolvable:
                raise ValueError(
                    "candidate or Newton-corrected residual exceeds the "
                    "unchanged nonlinear acceptance contract"
                )
            evidence = dataclasses.replace(
                evidence,
                label=(
                    evidence.label
                    + "; R5S7 sub-probe-resolution Newton correction: center "
                    "certificate governs (corrected residual "
                    f"{evidence.newton_corrected_maximum_scaled_residual:.6e}"
                    " recorded, not certified)"
                ),
            )
        ledger_tolerance = before.controls.ledger_tolerance
        if (
            face_event.ledger.maximum_step_ledger_residual > ledger_tolerance
            or face_event.ledger.maximum_cumulative_ledger_residual
            > ledger_tolerance
        ):
            raise ValueError(
                "open-interval event exceeds the unchanged conservation contract"
            )
        return evidence
    except _AtomicRejection:
        raise
    except Exception as exc:
        raise _AtomicRejection(
            CutMacrostepFailureCode.EVENT_TIME_ORDERING_UNCERTIFIED,
            f"open-interval numerical ordering failed closed: {exc}",
        ) from exc


def _departure_seed_from_overlap(
    overlap: tangent.FaceTangentOverlap,
    event_seed: fi.FaceEventSeed,
) -> di.FaceDepartureSeed:
    candidate = overlap.candidate
    event_scales = event_seed.stefan_flux_scales_mol_m2_s
    if len(candidate.dry_total_stefan_fluxes_mol_m2_s) != len(event_scales) + 1:
        raise _AtomicRejection(
            CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
            "departure added an unexpected number of dry faces",
        )
    newborn_scale = max(
        abs(candidate.dry_total_stefan_fluxes_mol_m2_s[0]),
        event_scales[0],
    )
    return di.FaceDepartureSeed(
        candidate,
        (newborn_scale, *event_scales),
        "atomic macrostep exact-face tangent / positive-sweep seed",
        overlap.newborn_trace_gradients,
    )


def _compose_macro_ledger(
    before: ci.CutIntegratorState,
    after: ci.CutIntegratorState | fi.FaceArrivalCommittedState,
    requested_dt_s: float,
    branch: CutMacrostepBranch,
    stage_ledgers: Sequence[object],
    *,
    event_duration_s: float,
    event_time_error_bound_s: float | None,
    departure_remainder_s: float,
    failed_same_cell_evaluations: int,
    uncommitted_continuation_roots: Sequence[object],
    event_ordering_mode: EventOrderingMode | None,
    event_ordering: tuple[str, ...],
) -> CutEventMacrostepLedger:
    ledgers = tuple(stage_ledgers)
    if not ledgers or not all(getattr(ledger, "accepted", False) for ledger in ledgers):
        raise _AtomicRejection(
            CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "only accepted physical balance stages may enter a macro ledger",
        )
    water_step = math.fsum(
        float(ledger.water_step_residual_mol) for ledger in ledgers
    )
    hexane_step = math.fsum(
        float(ledger.hexane_step_residual_mol) for ledger in ledgers
    )
    energy_step = math.fsum(float(ledger.energy_step_residual_j) for ledger in ledgers)
    boundary_water = math.fsum(float(ledger.boundary_water_out_mol) for ledger in ledgers)
    boundary_hexane = math.fsum(
        float(ledger.boundary_hexane_out_mol) for ledger in ledgers
    )
    boundary_energy = math.fsum(float(ledger.boundary_energy_out_j) for ledger in ledgers)

    old_inventory = ci.inventory_snapshot(before.transport)
    if isinstance(after, ci.CutIntegratorState):
        new_inventory = ci.inventory_snapshot(after.transport)
        # The macro telescope needs state inventories, not another local
        # heat-capacity stencil.  Each accepted stage has already certified
        # its phase-safe current capacity scale.  Reprobing a final gas-only
        # state at fixed composition can cross a phase boundary even though
        # the state and its conserved energy are valid.  Use the immutable
        # cumulative reference scale here; this is at least as strict whenever
        # it is smaller than a current-state diagnostic scale.
        after_capacity = after.reference_capacity_energy_scale_j
        after_time_s = after.transport.time_s
        after_accepted_steps = after.accepted_steps
    else:
        new_inventory = after.current_inventory
        after_capacity = after.reference_capacity_energy_scale_j
        after_time_s = after.time_s
        after_accepted_steps = after.accepted_steps
    water_scale = max(
        before.reference_inventory.total_water_mol,
        old_inventory.total_water_mol,
        new_inventory.total_water_mol,
        abs(boundary_water),
        1.0e-300,
    )
    hexane_scale = max(
        before.reference_inventory.total_hexane_mol,
        old_inventory.total_hexane_mol,
        new_inventory.total_hexane_mol,
        abs(boundary_hexane),
        1.0e-300,
    )
    energy_scale = max(
        before.reference_capacity_energy_scale_j,
        after_capacity,
        abs(boundary_energy),
        1.0e-300,
    )
    normalized_water = abs(water_step) / water_scale
    normalized_hexane = abs(hexane_step) / hexane_scale
    normalized_energy = abs(energy_step) / energy_scale

    final = ledgers[-1]
    maximum_stage_scaled = max(
        float(ledger.maximum_scaled_residual) for ledger in ledgers
    )
    maximum_stage_step = max(
        float(ledger.maximum_step_ledger_residual) for ledger in ledgers
    )
    maximum_final_cumulative = float(final.maximum_cumulative_ledger_residual)
    stage_evaluations = sum(int(ledger.nonlinear_evaluations) for ledger in ledgers)
    auxiliary_evaluations = sum(
        int(getattr(root, "nonlinear_evaluations", 0))
        for root in uncommitted_continuation_roots
    )
    if auxiliary_evaluations < 0 or auxiliary_evaluations > stage_evaluations:
        raise _AtomicRejection(
            CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "uncommitted continuation work is inconsistent with stage work",
        )
    accepted_evaluations = stage_evaluations - auxiliary_evaluations
    total_attempted_evaluations = (
        failed_same_cell_evaluations + stage_evaluations
    )
    balance_intervals = len(ledgers)
    if after_accepted_steps - before.accepted_steps != balance_intervals:
        raise _AtomicRejection(
            CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "accepted-step counter includes a missing or auxiliary balance interval",
        )
    if any(getattr(root, "committed", True) for root in uncommitted_continuation_roots):
        raise _AtomicRejection(
            CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "an auxiliary continuation root entered the committed ledger",
        )

    duration_closure = math.fsum(
        (event_duration_s, departure_remainder_s, -requested_dt_s)
        if event_duration_s > 0.0
        else (requested_dt_s, -requested_dt_s)
    )
    target_time_s = math.fsum((before.transport.time_s, requested_dt_s))
    absolute_time_closure = math.fsum((after_time_s, -target_time_s))
    time_roundoff = 8.0 * max(
        math.ulp(max(abs(before.transport.time_s), 1.0)),
        math.ulp(max(abs(target_time_s), 1.0)),
        math.ulp(max(abs(requested_dt_s), 1.0)),
    )
    if abs(duration_closure) > time_roundoff or abs(absolute_time_closure) > time_roundoff:
        raise _AtomicRejection(
            CutMacrostepFailureCode.TOPOLOGY_HANDOFF_FAILED,
            "event/remainder durations do not close the requested macrostep to "
            "floating-point roundoff",
        )

    ledger_tolerance = before.controls.ledger_tolerance
    nonlinear_tolerance = before.controls.nonlinear_residual_tolerance
    if maximum_stage_scaled > nonlinear_tolerance:
        raise _AtomicRejection(
            CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "an accepted stage exceeded the unchanged nonlinear residual contract",
        )
    if (
        maximum_stage_step > ledger_tolerance
        or maximum_final_cumulative > ledger_tolerance
        or max(normalized_water, normalized_hexane, normalized_energy)
        > ledger_tolerance
    ):
        raise _AtomicRejection(
            CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "stage, macrostep, or cumulative conservation exceeded 1e-10",
        )
    if event_ordering_mode is EventOrderingMode.OPEN_INTERVAL_NUMERICAL_ORDERING:
        if event_time_error_bound_s is not None:
            raise _AtomicRejection(
                CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
                "open-interval numerical ordering cannot carry an error bound",
            )
    elif event_ordering_mode is EventOrderingMode.RIGOROUS_INTERVAL_CERTIFICATE:
        if event_time_error_bound_s is None:
            raise _AtomicRejection(
                CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
                "rigorous interval ordering must retain its error bound",
            )
    elif event_duration_s > 0.0:
        raise _AtomicRejection(
            CutMacrostepFailureCode.MACRO_LEDGER_FAILED,
            "an event macrostep has no declared ordering authority",
        )

    return CutEventMacrostepLedger(
        branch=branch,
        requested_dt_s=requested_dt_s,
        event_duration_s=event_duration_s,
        event_time_error_bound_s=event_time_error_bound_s,
        departure_remainder_s=departure_remainder_s,
        duration_closure_error_s=duration_closure,
        absolute_time_closure_error_s=absolute_time_closure,
        accepted_balance_interval_count=balance_intervals,
        accepted_nonlinear_evaluations=accepted_evaluations,
        total_attempted_nonlinear_evaluations=total_attempted_evaluations,
        rejected_same_cell_nonlinear_evaluations=failed_same_cell_evaluations,
        uncommitted_continuation_root_count=len(uncommitted_continuation_roots),
        uncommitted_continuation_nonlinear_evaluations=auxiliary_evaluations,
        event_ordering_mode=event_ordering_mode,
        event_ordering=event_ordering,
        maximum_stage_scaled_residual=maximum_stage_scaled,
        maximum_stage_step_ledger_residual=maximum_stage_step,
        maximum_final_cumulative_ledger_residual=maximum_final_cumulative,
        water_step_residual_mol=water_step,
        hexane_step_residual_mol=hexane_step,
        energy_step_residual_j=energy_step,
        normalized_water_step_residual=normalized_water,
        normalized_hexane_step_residual=normalized_hexane,
        normalized_energy_step_residual=normalized_energy,
        water_cumulative_residual_mol=float(final.water_cumulative_residual_mol),
        hexane_cumulative_residual_mol=float(final.hexane_cumulative_residual_mol),
        energy_cumulative_residual_j=float(final.energy_cumulative_residual_j),
        normalized_water_cumulative_residual=float(
            final.normalized_water_cumulative_residual
        ),
        normalized_hexane_cumulative_residual=float(
            final.normalized_hexane_cumulative_residual
        ),
        normalized_energy_cumulative_residual=float(
            final.normalized_energy_cumulative_residual
        ),
        boundary_water_out_mol=boundary_water,
        boundary_hexane_out_mol=boundary_hexane,
        boundary_energy_out_j=boundary_energy,
        accepted=True,
    )


class ExactFaceEventSequenceError(RuntimeError):
    """Atomic sequence rejection carrying the original exact-face state."""

    def __init__(
        self,
        message: str,
        rollback_state: fi.FaceArrivalCommittedState,
        *,
        failed_stage_index: int,
        provisional_stages: Sequence[ffi.FaceToFaceStep],
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.failed_stage_index = failed_stage_index
        self.provisional_stages = tuple(provisional_stages)


@dataclass(frozen=True)
class ExactFaceEventSequenceLedger:
    """Independent atomic ledger for consecutive exact-face events."""

    stage_count: int
    total_duration_s: float
    duration_closure_error_s: float
    boundary_step_count: int
    nonlinear_evaluations: int
    maximum_stage_scaled_residual: float
    maximum_stage_step_ledger_residual: float
    maximum_final_cumulative_ledger_residual: float
    normalized_water_sequence_residual: float
    normalized_hexane_sequence_residual: float
    normalized_energy_sequence_residual: float
    normalized_water_raw_final_cumulative_residual: float
    normalized_hexane_raw_final_cumulative_residual: float
    normalized_energy_raw_final_cumulative_residual: float
    accepted: bool
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def maximum_sequence_ledger_residual(self) -> float:
        return max(
            self.normalized_water_sequence_residual,
            self.normalized_hexane_sequence_residual,
            self.normalized_energy_sequence_residual,
        )

    @property
    def maximum_raw_final_cumulative_residual(self) -> float:
        return max(
            self.normalized_water_raw_final_cumulative_residual,
            self.normalized_hexane_raw_final_cumulative_residual,
            self.normalized_energy_raw_final_cumulative_residual,
        )


@dataclass(frozen=True)
class ExactFaceEventSequence:
    """Atomically accepted adjacent-face events with aligned boundary steps."""

    before: fi.FaceArrivalCommittedState
    after: fi.FaceArrivalCommittedState
    stages: tuple[ffi.FaceToFaceStep, ...]
    boundaries: tuple[ct.PoreBoundary, ...]
    ledger: ExactFaceEventSequenceLedger
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.stages or len(self.stages) != len(self.boundaries):
            raise ValueError("exact-face sequence lost its stage/boundary rank")
        if self.stages[0].before is not self.before:
            raise ValueError("exact-face sequence lost its immutable start")
        if self.stages[-1].after is not self.after:
            raise ValueError("exact-face sequence lost its final state")
        if any(
            left.after is not right.before
            for left, right in zip(self.stages, self.stages[1:])
        ):
            raise ValueError("exact-face sequence contains a disconnected handoff")
        if not self.ledger.accepted:
            raise ValueError("exact-face sequence requires an accepted atomic ledger")


def advance_exact_face_event_sequence(
    before: fi.FaceArrivalCommittedState,
    solver_controls: ci.CutSolverControls,
    event_controls: Sequence[fi.FaceEventControls],
    boundaries: Sequence[ct.PoreBoundary],
    seeds: Sequence[ffi.FaceToFaceSeed],
    *,
    boundary_discontinuity_times_s: Sequence[float] = (),
) -> ExactFaceEventSequence:
    """Solve consecutive adjacent-face events and commit them atomically.

    Each supplied boundary is constant over its corresponding unknown-duration
    event.  A change between consecutive entries is aligned with the exact
    shared face state and is therefore representable without smearing.  A
    declared discontinuity strictly inside a solved event rejects the complete
    sequence; no event-time root before that rejection is committed.
    """

    provisional: list[ffi.FaceToFaceStep] = []
    stage_index = 0
    try:
        if not isinstance(before, fi.FaceArrivalCommittedState):
            raise TypeError("exact-face sequence requires a committed face state")
        controls = tuple(event_controls)
        boundary_values = tuple(boundaries)
        seed_values = tuple(seeds)
        if not controls or not (
            len(controls) == len(boundary_values) == len(seed_values)
        ):
            raise ValueError(
                "exact-face sequence requires equally ranked non-empty inputs"
            )
        if not all(
            isinstance(
                boundary,
                (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
            )
            for boundary in boundary_values
        ):
            raise TypeError(
                "exact-face sequence boundaries must be supported pore boundaries"
            )
        discontinuities = _validated_discontinuity_times(
            boundary_discontinuity_times_s
        )
        current = before
        for stage_index, (event_control, boundary, seed) in enumerate(
            zip(controls, boundary_values, seed_values)
        ):
            try:
                step = ffi.solve_adjacent_face_event(
                    current,
                    solver_controls,
                    event_control,
                    boundary,
                    seed,
                )
            except ffi.FaceToFaceIntegratorStepError as exc:
                if exc.rollback_state is not current:
                    raise RuntimeError(
                        "face-to-face rejection lost the current exact-face state"
                    ) from exc
                raise
            if step.before is not current:
                raise RuntimeError("face-to-face stage lost its exact handoff object")
            if step.after.arrival_face_index != current.arrival_face_index - 1:
                raise RuntimeError("face-to-face stage skipped or repeated a face")
            if (
                step.ledger.maximum_scaled_residual
                > solver_controls.nonlinear_residual_tolerance
                or step.ledger.maximum_step_ledger_residual
                > solver_controls.ledger_tolerance
                or step.ledger.maximum_cumulative_ledger_residual
                > solver_controls.ledger_tolerance
            ):
                raise RuntimeError("face-to-face stage exceeded an unchanged gate")
            provisional.append(step)
            _reject_discontinuity_inside_exact_event(
                current.time_s,
                step.after.time_s,
                discontinuities,
            )
            current = step.after
        ledger = _compose_exact_face_sequence_ledger(
            before,
            current,
            provisional,
            boundary_values,
            solver_controls,
        )
        return ExactFaceEventSequence(
            before,
            current,
            tuple(provisional),
            boundary_values,
            ledger,
        )
    except ExactFaceEventSequenceError:
        raise
    except Exception as exc:
        if not isinstance(before, fi.FaceArrivalCommittedState):
            raise
        raise ExactFaceEventSequenceError(
            f"atomic exact-face sequence rejected with exact rollback: {exc}",
            before,
            failed_stage_index=stage_index,
            provisional_stages=provisional,
        ) from exc


def _validated_discontinuity_times(values: Sequence[float]) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if not all(math.isfinite(value) and value >= 0.0 for value in result):
        raise ValueError("boundary discontinuity times must be finite and non-negative")
    if any(right <= left for left, right in zip(result, result[1:])):
        raise ValueError("boundary discontinuity times must increase strictly")
    return result


def _reject_discontinuity_inside_exact_event(
    start_time_s: float,
    end_time_s: float,
    discontinuities: Sequence[float],
) -> None:
    if any(start_time_s < value < end_time_s for value in discontinuities):
        raise ValueError(
            "boundary discontinuity lies strictly inside an unknown-duration "
            "exact-face event; split at the physical step"
        )


def _compose_exact_face_sequence_ledger(
    before: fi.FaceArrivalCommittedState,
    after: fi.FaceArrivalCommittedState,
    stages: Sequence[ffi.FaceToFaceStep],
    boundaries: Sequence[ct.PoreBoundary],
    controls: ci.CutSolverControls,
) -> ExactFaceEventSequenceLedger:
    ledgers = tuple(stage.ledger for stage in stages)
    durations = tuple(ledger.event_time_s for ledger in ledgers)
    total_duration = math.fsum(durations)
    duration_closure = math.fsum((after.time_s, -before.time_s, -total_duration))
    roundoff = 8.0 * max(
        math.ulp(max(abs(before.time_s), 1.0)),
        math.ulp(max(abs(after.time_s), 1.0)),
        math.ulp(max(abs(total_duration), 1.0)),
    )
    if abs(duration_closure) > roundoff:
        raise RuntimeError("exact-face sequence duration failed roundoff closure")
    if after.accepted_steps - before.accepted_steps != len(stages):
        raise RuntimeError("exact-face sequence accepted-step counter lost a stage")
    water_sequence = math.fsum(
        ledger.water_step_residual_mol for ledger in ledgers
    )
    hexane_sequence = math.fsum(
        ledger.hexane_step_residual_mol for ledger in ledgers
    )
    energy_sequence = math.fsum(
        ledger.energy_step_residual_j for ledger in ledgers
    )
    boundary_water = math.fsum(
        ledger.boundary_water_out_mol for ledger in ledgers
    )
    boundary_hexane = math.fsum(
        ledger.boundary_hexane_out_mol for ledger in ledgers
    )
    boundary_energy = math.fsum(
        ledger.boundary_energy_out_j for ledger in ledgers
    )
    water_scale = max(
        before.reference_inventory.total_water_mol,
        before.current_inventory.total_water_mol,
        after.current_inventory.total_water_mol,
        abs(boundary_water),
        1.0e-300,
    )
    hexane_scale = max(
        before.reference_inventory.total_hexane_mol,
        before.current_inventory.total_hexane_mol,
        after.current_inventory.total_hexane_mol,
        abs(boundary_hexane),
        1.0e-300,
    )
    energy_scale = max(
        before.reference_capacity_energy_scale_j,
        ffi._capacity_scale(stages[-1].assembly),  # noqa: SLF001
        abs(boundary_energy),
        1.0e-300,
    )
    sequence_normalized = (
        abs(water_sequence) / water_scale,
        abs(hexane_sequence) / hexane_scale,
        abs(energy_sequence) / energy_scale,
    )
    raw_water = math.fsum(
        (
            *after.cumulative_material_water_change_cell_mol,
            *after.water_change_compensation_cell_mol,
            after.cumulative_boundary_water_out_mol,
            after.boundary_water_compensation_mol,
        )
    )
    raw_hexane = math.fsum(
        (
            *after.cumulative_material_hexane_change_cell_mol,
            *after.hexane_change_compensation_cell_mol,
            after.cumulative_boundary_hexane_out_mol,
            after.boundary_hexane_compensation_mol,
        )
    )
    raw_energy = math.fsum(
        (
            *after.cumulative_material_energy_change_cell_j,
            *after.energy_change_compensation_cell_j,
            after.cumulative_boundary_energy_out_j,
            after.boundary_energy_compensation_j,
        )
    )
    corrected_water_out = after.corrected_cumulative_boundary_water_out_mol
    corrected_hexane_out = after.corrected_cumulative_boundary_hexane_out_mol
    raw_water_scale = max(
        after.reference_inventory.total_water_mol,
        after.current_inventory.total_water_mol,
        abs(corrected_water_out),
        after.cumulative_absolute_water_transfer_mol,
        1.0e-300,
    )
    raw_hexane_scale = max(
        after.reference_inventory.total_hexane_mol,
        after.current_inventory.total_hexane_mol,
        abs(corrected_hexane_out),
        after.cumulative_absolute_hexane_transfer_mol,
        1.0e-300,
    )
    raw_energy_scale = max(
        after.reference_capacity_energy_scale_j,
        ffi._capacity_scale(stages[-1].assembly),  # noqa: SLF001
        after.cumulative_absolute_energy_transfer_j,
        1.0e-300,
    )
    raw_normalized = (
        abs(raw_water) / raw_water_scale,
        abs(raw_hexane) / raw_hexane_scale,
        abs(raw_energy) / raw_energy_scale,
    )
    maximum_stage_scaled = max(
        ledger.maximum_scaled_residual for ledger in ledgers
    )
    maximum_stage_step = max(
        ledger.maximum_step_ledger_residual for ledger in ledgers
    )
    maximum_final_cumulative = ledgers[-1].maximum_cumulative_ledger_residual
    if maximum_stage_scaled > controls.nonlinear_residual_tolerance:
        raise RuntimeError("exact-face sequence contains an uncertified root")
    if max(
        maximum_stage_step,
        maximum_final_cumulative,
        *sequence_normalized,
        *raw_normalized,
    ) > controls.ledger_tolerance:
        raise RuntimeError(
            "exact-face sequence stage, telescope, or raw cumulative ledger "
            "exceeded the unchanged 1e-10 gate"
        )
    boundary_step_count = sum(
        left != right for left, right in zip(boundaries, boundaries[1:])
    )
    return ExactFaceEventSequenceLedger(
        stage_count=len(stages),
        total_duration_s=total_duration,
        duration_closure_error_s=duration_closure,
        boundary_step_count=boundary_step_count,
        nonlinear_evaluations=sum(
            ledger.nonlinear_evaluations for ledger in ledgers
        ),
        maximum_stage_scaled_residual=maximum_stage_scaled,
        maximum_stage_step_ledger_residual=maximum_stage_step,
        maximum_final_cumulative_ledger_residual=maximum_final_cumulative,
        normalized_water_sequence_residual=sequence_normalized[0],
        normalized_hexane_sequence_residual=sequence_normalized[1],
        normalized_energy_sequence_residual=sequence_normalized[2],
        normalized_water_raw_final_cumulative_residual=raw_normalized[0],
        normalized_hexane_raw_final_cumulative_residual=raw_normalized[1],
        normalized_energy_raw_final_cumulative_residual=raw_normalized[2],
        accepted=True,
    )


__all__ = [
    "CutEventMacrostep",
    "CutEventMacrostepError",
    "CutEventMacrostepLedger",
    "CutMacrostepBranch",
    "CutMacrostepFailureCode",
    "AcceptedHistoryFaceEventBracket",
    "EventOrderingEvidence",
    "EventOrderingMode",
    "ExactFaceEventSequence",
    "ExactFaceEventSequenceError",
    "ExactFaceEventSequenceLedger",
    "FaceEventTimeCertificate",
    "FaceEventTimeCertifier",
    "OpenIntervalNumericalOrderingEvidence",
    "advance_exact_face_event_sequence",
    "advance_bracketed_interior_face_macrostep",
    "advance_one_interior_face_macrostep",
    "bracket_inner_face_from_accepted_history",
    "face_event_seed_from_accepted_history_bracket",
]
