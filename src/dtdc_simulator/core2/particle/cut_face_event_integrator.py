r"""Safeguarded nonlinear solve for the exact inner-face Gate-1g event.

The event time is an unknown and the final front is the exact inner master
face.  Open physical charts are represented by one-to-one logit transforms;
no variable is clipped and no epsilon phase is created.  Small systems use a
dense finite-difference oracle, while larger systems use a conservative
bordered-banded structural pattern with ``scipy.optimize.least_squares``.

An accepted result commits an exact-face payload.  It is intentionally *not*
a :class:`cut_transport.CutTransportState`: the next finite recession needs a
separate zero-old-dry-subpiece face-departure chart.  This module records all
material fields and signed compensated ledgers needed by that future chart,
but does not invent an epsilon restart.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field, replace
from typing import Sequence

import numpy as np
from scipy import optimize, sparse

from dtdc_simulator.core2.particle import conditioning
from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_face_event as face
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp


class FaceDepartureChartRequired(RuntimeError):
    """An exact face state cannot be used as a strict partial-cut state."""


class FaceEventIntegratorStepError(RuntimeError):
    """Rejected exact-face solve carrying the immutable rollback state."""

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState,
        *,
        nonlinear_evaluations: int,
        rejected_trial_evaluations: int,
        maximum_scaled_residual: float = math.inf,
        condition_proxy: float = math.inf,
        last_candidate: face.FaceArrivalUnknowns | None = None,
        scaled_residuals: tuple[float, ...] = (),
        residual_blocks: face.FaceArrivalResidualBlocks | None = None,
        minimum_fractional_distance_to_bound: float = 0.0,
        optimizer_status: int | None = None,
        optimizer_success: bool | None = None,
        optimizer_message: str | None = None,
        optimizer_nfev: int | None = None,
        optimizer_njev: int | None = None,
        optimizer_cost: float | None = None,
        optimizer_optimality: float | None = None,
        optimizer_active_mask: tuple[int, ...] = (),
        retained_water_applicability_guard_failed: bool = False,
        retained_water_applicability_audit: (
            ci.WetRetainedWaterApplicabilityAudit | None
        ) = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.nonlinear_evaluations = nonlinear_evaluations
        self.rejected_trial_evaluations = rejected_trial_evaluations
        self.maximum_scaled_residual = maximum_scaled_residual
        self.condition_proxy = condition_proxy
        self.last_candidate = last_candidate
        self.scaled_residuals = scaled_residuals
        self.residual_blocks = residual_blocks
        self.minimum_fractional_distance_to_bound = (
            minimum_fractional_distance_to_bound
        )
        self.optimizer_status = optimizer_status
        self.optimizer_success = optimizer_success
        self.optimizer_message = optimizer_message
        self.optimizer_nfev = optimizer_nfev
        self.optimizer_njev = optimizer_njev
        self.optimizer_cost = optimizer_cost
        self.optimizer_optimality = optimizer_optimality
        self.optimizer_active_mask = optimizer_active_mask
        self.retained_water_applicability_guard_failed = (
            retained_water_applicability_guard_failed
        )
        self.retained_water_applicability_audit = (
            retained_water_applicability_audit
        )


class FaceEventTimeLinearizationAuditError(RuntimeError):
    """The accepted face root has no valid local physical-Jacobian audit."""


@dataclass(frozen=True)
class FaceEventControls:
    """The open chart needed when event duration is an unknown.

    The declared bounds may start at exactly zero.  The logit map still
    represents only the strict interior, so ``(0, upper)`` introduces no
    positive epsilon time and can never return the initial endpoint.
    """

    event_time_bounds_s: tuple[float, float]
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        lower, upper = self.event_time_bounds_s
        if not all(math.isfinite(value) for value in (lower, upper)):
            raise ValueError("event-time bounds must be finite")
        if not 0.0 <= lower < upper:
            raise ValueError(
                "event-time bounds need a non-negative lower endpoint and "
                "positive ordered width"
            )


@dataclass(frozen=True)
class FaceEventSeed:
    """One explicit event seed and physical scales for dry Stefan fluxes."""

    candidate: face.FaceArrivalUnknowns
    stefan_flux_scales_mol_m2_s: tuple[float, ...]
    label: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.candidate, face.FaceArrivalUnknowns):
            raise TypeError("face-event seed requires FaceArrivalUnknowns")
        if not self.stefan_flux_scales_mol_m2_s or not all(
            math.isfinite(value) and value > 0.0
            for value in self.stefan_flux_scales_mol_m2_s
        ):
            raise ValueError("every event Stefan face needs a positive scale")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("face-event seed needs a provenance label")


@dataclass(frozen=True)
class FaceEventResidualScales:
    """Fixed component-inventory and caloric-capacity row scales."""

    wet_water_mol_s: tuple[float, ...]
    wet_energy_w: tuple[float, ...]
    dry_water_mol_s: tuple[float, ...]
    dry_hexane_mol_s: tuple[float, ...]
    dry_energy_w: tuple[float, ...]
    rh_water_mol_s: float
    rh_hexane_mol_s: float
    rh_energy_w: float
    dry_energy_capacity_stencil_audits: tuple[
        ci.DryEnergyCapacityStencilAudit,
        ...,
    ]
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.dry_energy_capacity_stencil_audits, tuple):
            raise TypeError("face-event dry-capacity audits must be an immutable tuple")
        if len(self.dry_energy_capacity_stencil_audits) < len(self.dry_energy_w):
            raise ValueError(
                "every face-event dry-energy scale needs an Amendment-19 "
                "stencil audit"
            )
        if not all(
            isinstance(audit, ci.DryEnergyCapacityStencilAudit)
            for audit in self.dry_energy_capacity_stencil_audits
        ):
            raise TypeError("face-event dry-capacity audit payload has the wrong type")

    def __setstate__(self, state: dict[str, object]) -> None:
        """Restore literal-name identity for byte-stable restart payloads."""

        for name, value in state.items():
            object.__setattr__(self, name, value)
        audits = self.dry_energy_capacity_stencil_audits
        if isinstance(audits, tuple) and all(
            isinstance(audit, ci.DryEnergyCapacityStencilAudit)
            for audit in audits
        ):
            object.__setattr__(
                self,
                "dry_energy_capacity_stencil_audits",
                tuple(
                    replace(
                        audit,
                        primitive_domain_kind=sys.intern(
                            audit.primitive_domain_kind
                        ),
                        stencil_kind=sys.intern(audit.stencil_kind),
                    )
                    for audit in audits
                ),
            )
        self.__post_init__()

    @property
    def vector(self) -> tuple[float, ...]:
        result = (
            *self.wet_water_mol_s,
            *self.wet_energy_w,
            *self.dry_water_mol_s,
            *self.dry_hexane_mol_s,
            *self.dry_energy_w,
            self.rh_water_mol_s,
            self.rh_hexane_mol_s,
            self.rh_energy_w,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in result):
            raise RuntimeError("every face-event residual scale must be positive")
        return result


@dataclass(frozen=True)
class FaceEventLedger:
    """Nonlinear, conditioning, signed step, and cumulative diagnostics."""

    event_time_s: float
    nonlinear_evaluations: int
    rejected_trial_evaluations: int
    nonlinear_message: str
    residual_scales: FaceEventResidualScales
    scaled_residuals: tuple[float, ...]
    maximum_scaled_residual: float
    condition_proxy: float
    minimum_fractional_distance_to_bound: float
    used_sparse_jacobian: bool
    scalability_class: str
    jacobian_structural_color_count: int
    jacobian_structural_nnz: int
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
class FaceArrivalCommittedState:
    """Exact-face handoff retaining everything a departure chart will need."""

    geometry: cg.CutGeometry
    config: cut.CutTransportConfig
    solver_controls: ci.CutSolverControls
    time_s: float
    arrival_face_index: int
    wet_temperatures_k: tuple[float, ...]
    wet_retained_water_loadings: tuple[float, ...]
    dry_temperatures_k: tuple[float, ...]
    dry_y_hexane: tuple[float, ...]
    historical_hexane_loadings: tuple[float, ...]
    oil_fraction_labels: tuple[float, ...]
    interface_temperature_k: float
    dry_total_stefan_fluxes_mol_m2_s: tuple[float, ...]
    current_inventory: ci.MaterialInventorySnapshot
    reference_inventory: ci.MaterialInventorySnapshot
    reference_capacity_energy_scale_j: float
    cumulative_boundary_water_out_mol: float
    cumulative_boundary_hexane_out_mol: float
    cumulative_boundary_energy_out_j: float
    boundary_water_compensation_mol: float
    boundary_hexane_compensation_mol: float
    boundary_energy_compensation_j: float
    cumulative_absolute_water_transfer_mol: float
    cumulative_absolute_hexane_transfer_mol: float
    cumulative_absolute_energy_transfer_j: float
    cumulative_material_water_change_cell_mol: tuple[float, ...]
    cumulative_material_hexane_change_cell_mol: tuple[float, ...]
    cumulative_material_energy_change_cell_j: tuple[float, ...]
    water_change_compensation_cell_mol: tuple[float, ...]
    hexane_change_compensation_cell_mol: tuple[float, ...]
    energy_change_compensation_cell_j: tuple[float, ...]
    accepted_steps: int
    cumulative_nonlinear_evaluations: int
    transition_authority: str = (
        "next recession requires the zero-old-dry-subpiece face-departure chart"
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.solver_controls, ci.CutSolverControls):
            raise TypeError("committed face state must retain CutSolverControls")
        geometry = self.geometry
        if (
            geometry.cut_cell_index is not None
            or not geometry.front.at_master_face
            or geometry.front.master_face_index != self.arrival_face_index
        ):
            raise ValueError("committed face state must have exact no-cut geometry")
        if not math.isfinite(self.time_s) or self.time_s < 0.0:
            raise ValueError("committed event time must be finite and non-negative")
        n = geometry.master_grid.n
        if len(self.wet_temperatures_k) != self.arrival_face_index:
            raise ValueError("committed wet fields do not end at the arrival face")
        if len(self.dry_temperatures_k) != n - self.arrival_face_index:
            raise ValueError("committed dry fields do not start at the arrival face")
        for values in (
            self.historical_hexane_loadings,
            self.oil_fraction_labels,
            self.cumulative_material_water_change_cell_mol,
            self.cumulative_material_hexane_change_cell_mol,
            self.cumulative_material_energy_change_cell_j,
            self.water_change_compensation_cell_mol,
            self.hexane_change_compensation_cell_mol,
            self.energy_change_compensation_cell_j,
        ):
            if len(values) != n:
                raise ValueError("committed material arrays must align with the master grid")

    @property
    def requires_face_departure_chart(self) -> bool:
        return True

    @property
    def is_strict_cut_transport_state(self) -> bool:
        return False

    def to_strict_cut_transport_state(self) -> cut.CutTransportState:
        raise FaceDepartureChartRequired(
            f"{self.transition_authority}; an epsilon-z restart is forbidden"
        )

    @property
    def corrected_cumulative_boundary_water_out_mol(self) -> float:
        return math.fsum(
            (
                self.cumulative_boundary_water_out_mol,
                self.boundary_water_compensation_mol,
            )
        )

    @property
    def corrected_cumulative_boundary_hexane_out_mol(self) -> float:
        return math.fsum(
            (
                self.cumulative_boundary_hexane_out_mol,
                self.boundary_hexane_compensation_mol,
            )
        )

    @property
    def corrected_cumulative_boundary_energy_out_j(self) -> float:
        return math.fsum(
            (
                self.cumulative_boundary_energy_out_j,
                self.boundary_energy_compensation_j,
            )
        )


@dataclass(frozen=True)
class FaceEventStep:
    """One certified exact inner-face arrival and its handoff payload."""

    before: ci.CutIntegratorState
    after: FaceArrivalCommittedState
    assembly: face.FaceArrivalAssembly
    boundary: ct.PoreBoundary
    controls: FaceEventControls
    seed: FaceEventSeed
    ledger: FaceEventLedger
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceEventTimeLinearizationAudit:
    r"""Local physical-Jacobian evidence, deliberately not a time bound.

    The dimensionless Newton correction solves

    ``J_star delta_x_star = -R_star``

    at the accepted exact-face candidate, where ``J_star`` is certified by a
    fresh central finite difference in physical variables.  This establishes
    local rank and reports the linearized event-time correction without
    assigning it the semantics of an error interval.

    A rigorous event-ordering certificate additionally needs an enclosure of
    the nonlinear Jacobian over a neighborhood (for example, a validated
    interval or a proven Jacobian-Lipschitz bound).  The current constitutive
    stack exposes no such authority.  Therefore the two ordering fields below
    are intentionally fixed to their fail-closed values.
    """

    physical_variable_scale_provenance: str
    physical_variable_scale_labels: tuple[str, ...]
    physical_variable_scales: tuple[float, ...]
    physical_jacobian_perturbations: tuple[float, ...]
    jacobian_certificate: conditioning.UnitInvariantJacobianCertificate
    scaled_residuals: tuple[float, ...]
    dimensionless_newton_correction: tuple[float, ...]
    physical_newton_correction: tuple[float, ...]
    newton_corrected_scaled_residuals: tuple[float, ...]
    newton_corrected_maximum_scaled_residual: float
    event_time_variable_index: int
    linearized_event_time_correction_s: float
    maximum_dimensionless_newton_correction: float
    event_time_absolute_error_bound_s: None = field(default=None, init=False)
    nonlinear_jacobian_enclosure_available: bool = field(default=False, init=False)
    ordering_certified: bool = field(default=False, init=False)
    missing_ordering_authority: str = field(
        default=(
            "no validated neighborhood Jacobian-Lipschitz or interval enclosure "
            "is available for the coupled face-arrival residual"
        ),
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class _CoordinateChart:
    wet_temperature: tuple[float, float]
    wet_water: tuple[float, float]
    dry_temperature: tuple[float, float]
    dry_y_hexane: tuple[float, float] | None
    dry_config: ct.FullyDryTransportConfig
    dry_pores: tuple[cp.CoupledPoreParams, ...]
    event_time: tuple[float, float]
    interface_temperature: tuple[float, float]


@dataclass(frozen=True)
class _EventSparsity:
    matrix: sparse.csr_matrix
    color_count: int

    @property
    def nnz(self) -> int:
        return int(self.matrix.nnz)


_FACE_EVENT_STEFAN_SCALE_FLOOR_MOL_M2_S = 1.0e-3
_FACE_EVENT_JACOBIAN_RANK_RELATIVE_TOLERANCE = 1.0e-13
_FACE_EVENT_PHYSICAL_PROBE_RELATIVE_STEP = 1.0e-5
_FACE_EVENT_PHYSICAL_PROBE_MAX_HALVINGS = 16


def solve_inner_face_arrival(
    before: ci.CutIntegratorState,
    controls: FaceEventControls,
    boundary: ct.PoreBoundary,
    seed: FaceEventSeed,
) -> FaceEventStep:
    """Solve and certify one exact inner-master-face event or roll back."""

    evaluations = 0
    rejected_trials = 0
    maximum_scaled = math.inf
    condition_proxy = math.inf
    last_candidate = None
    last_scaled: tuple[float, ...] = ()
    last_blocks = None
    minimum_bound_distance = 0.0
    optimizer_status = None
    optimizer_success = None
    optimizer_message = None
    optimizer_nfev = None
    optimizer_njev = None
    optimizer_cost = None
    optimizer_optimality = None
    optimizer_active_mask: tuple[int, ...] = ()
    try:
        if not isinstance(before, ci.CutIntegratorState):
            raise TypeError("before must be a CutIntegratorState")
        if not isinstance(controls, FaceEventControls):
            raise TypeError("controls must be FaceEventControls")
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("boundary must be a supported pore boundary")
        layout = face.layout_for_arrival(before.transport)
        _validate_event_chart(before, controls)
        _validate_seed(before, controls, seed, layout)
        chart = _coordinate_chart(before, controls)
        initial_coordinates = _encode_candidate(
            seed.candidate,
            chart,
            seed.stefan_flux_scales_mol_m2_s,
        )
        initial = face.assemble_face_arrival(
            before.transport,
            seed.candidate,
            boundary,
            enforce_reduced_film_thresholds=False,
        )
        evaluations += 1
        scales = _residual_scales(before.transport, initial)
        scale_vector = np.asarray(scales.vector, dtype=float)
        use_sparse = layout.unknown_count > before.controls.maximum_dense_unknowns
        structural = _event_sparsity(layout) if use_sparse else None

        def scaled_residual(coordinates: np.ndarray) -> np.ndarray:
            nonlocal evaluations, rejected_trials
            candidate = _decode_candidate(
                coordinates,
                layout,
                chart,
                seed.stefan_flux_scales_mol_m2_s,
            )
            try:
                assembly = face.assemble_face_arrival(
                    before.transport,
                    candidate,
                    boundary,
                    enforce_reduced_film_thresholds=False,
                )
            except face.FaceArrivalStepError:
                rejected_trials += 1
                raise
            evaluations += 1
            return np.asarray(
                _datum_covariant_residual_vector(assembly.residuals),
                dtype=float,
            ) / scale_vector

        lower, upper = ci._coordinate_solver_bounds(  # noqa: SLF001
            layout,
            before.controls.maximum_logit_magnitude,
        )
        solution = optimize.least_squares(
            scaled_residual,
            initial_coordinates,
            method="dogbox" if structural is not None else "trf",
            bounds=(lower, upper),
            ftol=before.controls.nonlinear_step_tolerance,
            xtol=before.controls.nonlinear_step_tolerance,
            gtol=before.controls.nonlinear_step_tolerance,
            max_nfev=before.controls.maximum_function_evaluations,
            x_scale="jac",
            jac_sparsity=None if structural is None else structural.matrix,
            tr_options=(
                {"atol": 1.0e-14, "btol": 1.0e-14, "maxiter": 1000}
                if structural is not None
                else None
            ),
        )
        optimizer_status = int(solution.status)
        optimizer_success = bool(solution.success)
        optimizer_message = str(solution.message)
        optimizer_nfev = int(solution.nfev)
        optimizer_njev = None if solution.njev is None else int(solution.njev)
        optimizer_cost = float(solution.cost)
        optimizer_optimality = float(solution.optimality)
        optimizer_active_mask = tuple(int(value) for value in solution.active_mask)
        candidate = _decode_candidate(
            solution.x,
            layout,
            chart,
            seed.stefan_flux_scales_mol_m2_s,
        )
        converged = face.assemble_face_arrival(before.transport, candidate, boundary)
        evaluations += 1
        covariant = _datum_covariant_residual_vector(converged.residuals)
        signed_scaled = tuple(
            residual / scale for residual, scale in zip(covariant, scales.vector)
        )
        absolute_scaled = tuple(abs(value) for value in signed_scaled)
        last_candidate = candidate
        last_scaled = signed_scaled
        last_blocks = converged.residuals
        fractions = _bounded_fractions(candidate, chart)
        minimum_bound_distance = min(
            min(value, 1.0 - value) for value in fractions
        )
        maximum_scaled = max(absolute_scaled)
        condition_proxy = ci._jacobian_condition_proxy(solution.jac)  # noqa: SLF001
        if not solution.success:
            raise RuntimeError(
                f"face-event nonlinear root did not certify success: {solution.message}"
            )
        if np.any(solution.active_mask != 0):
            raise RuntimeError(
                "face-event solution touched a logit guard; open-chart root failed"
            )
        if maximum_scaled > before.controls.nonlinear_residual_tolerance:
            raise RuntimeError(
                "face-event regional/RH residual contract failed: "
                f"{maximum_scaled:.3e} > "
                f"{before.controls.nonlinear_residual_tolerance:.3e}"
            )
        if not math.isfinite(condition_proxy) or (
            condition_proxy > before.controls.maximum_condition_proxy
        ):
            raise RuntimeError(
                "face-event Jacobian condition proxy is uncertified: "
                f"{condition_proxy:.3e}"
            )
        step = _accept_event(
            before,
            controls,
            boundary,
            seed,
            converged,
            scales,
            signed_scaled,
            evaluations,
            rejected_trials,
            str(solution.message),
            condition_proxy,
            chart,
            structural,
        )
        if (
            step.ledger.maximum_step_ledger_residual
            > before.controls.ledger_tolerance
            or step.ledger.maximum_cumulative_ledger_residual
            > before.controls.ledger_tolerance
        ):
            raise RuntimeError(
                "face-event accepted/cumulative ledger exceeded the 1e-10 contract"
            )
        ci.certify_explicit_wet_retained_water_applicability(
            before.controls,
            before.transport.wet_retained_water_loadings,
            before.transport.effective_wet_retained_water_capacity_duals_over_rt,
            before.transport.config.wet.luikov,
        )
        ci.certify_explicit_wet_retained_water_applicability(
            before.controls,
            step.after.wet_retained_water_loadings,
            (0.0,) * len(step.after.wet_retained_water_loadings),
            step.after.config.wet.luikov,
        )
        return step
    except FaceEventIntegratorStepError:
        raise
    except Exception as exc:
        if not isinstance(before, ci.CutIntegratorState):
            raise
        applicability_error = (
            exc
            if isinstance(exc, ci.WetRetainedWaterApplicabilityError)
            else None
        )
        raise FaceEventIntegratorStepError(
            f"exact face-event solve rejected with exact rollback: {exc}",
            before,
            nonlinear_evaluations=evaluations,
            rejected_trial_evaluations=rejected_trials,
            maximum_scaled_residual=maximum_scaled,
            condition_proxy=condition_proxy,
            last_candidate=last_candidate,
            scaled_residuals=last_scaled,
            residual_blocks=last_blocks,
            minimum_fractional_distance_to_bound=minimum_bound_distance,
            optimizer_status=optimizer_status,
            optimizer_success=optimizer_success,
            optimizer_message=optimizer_message,
            optimizer_nfev=optimizer_nfev,
            optimizer_njev=optimizer_njev,
            optimizer_cost=optimizer_cost,
            optimizer_optimality=optimizer_optimality,
            optimizer_active_mask=optimizer_active_mask,
            retained_water_applicability_guard_failed=(
                applicability_error is not None
            ),
            retained_water_applicability_audit=(
                None if applicability_error is None else applicability_error.audit
            ),
        ) from exc


def audit_face_event_time_linearization(
    step: FaceEventStep,
) -> FaceEventTimeLinearizationAudit:
    r"""Independently audit the accepted root's local physical linearization.

    This function reassembles central probes in the physical unknown vector;
    it does not reuse the solver-coordinate Jacobian.  The resulting rank,
    condition number, and Newton correction are unit-invariant under the
    declared ``D_R``/``D_X`` contract.

    The result is intentionally not convertible into an event-time ordering
    certificate: a point Jacobian cannot bound the nonlinear correction tail.
    """

    try:
        if not isinstance(step, FaceEventStep):
            raise TypeError("face-event time audit requires a FaceEventStep")
        if not step.ledger.accepted:
            raise ValueError("face-event time audit requires an accepted step")
        if step.assembly.before is not step.before.transport:
            raise ValueError("face-event audit lost its immutable before state")
        layout = step.assembly.layout
        if not layout.is_square:
            raise ValueError("face-event audit requires a square arrival rank")

        labels, variable_scales = _face_event_physical_variable_scale_contract(
            step
        )
        center = np.asarray(step.assembly.candidate.vector(), dtype=float)
        residual_scales = np.asarray(step.ledger.residual_scales.vector, dtype=float)
        center_residual = _face_event_physical_residual(step, center)
        scaled_residual = center_residual / residual_scales
        physical_jacobian, perturbations = _face_event_physical_jacobian(
            step,
            center,
            np.asarray(variable_scales, dtype=float),
        )
        certificate = conditioning.certify_unit_invariant_jacobian(
            physical_jacobian,
            residual_scales,
            variable_scales,
            rank_relative_tolerance=(
                _FACE_EVENT_JACOBIAN_RANK_RELATIVE_TOLERANCE
            ),
        )
        if (
            certificate.condition_number_2
            > step.before.controls.maximum_condition_proxy
        ):
            raise conditioning.JacobianCertificateError(
                "face-event physical Jacobian exceeds the declared condition "
                f"ceiling: {certificate.condition_number_2:.3e} > "
                f"{step.before.controls.maximum_condition_proxy:.3e}"
            )
        dimensionless_jacobian = np.asarray(
            certificate.dimensionless_jacobian,
            dtype=float,
        )
        try:
            normalized_correction = np.linalg.solve(
                dimensionless_jacobian,
                -scaled_residual,
            )
        except np.linalg.LinAlgError as exc:
            raise conditioning.JacobianCertificateError(
                "face-event certified Jacobian could not produce a Newton correction"
            ) from exc
        physical_correction = np.asarray(variable_scales) * normalized_correction
        if not np.all(np.isfinite(physical_correction)):
            raise conditioning.JacobianCertificateError(
                "face-event physical Newton correction is non-finite"
            )
        corrected_residual = _face_event_physical_residual(
            step,
            center + physical_correction,
        )
        corrected_scaled_residual = corrected_residual / residual_scales
        event_time_index = (
            2 * layout.wet_piece_count
            + 2 * layout.dry_piece_count
            + layout.dry_face_count
        )
        return FaceEventTimeLinearizationAudit(
            physical_variable_scale_provenance=(
                "declared feasibility widths; unit mole-fraction scale for dry y; "
                "accepted-before Stefan history with a fixed 1e-3 mol m-2 s-1 "
                "floor; declared event-time chart width"
            ),
            physical_variable_scale_labels=labels,
            physical_variable_scales=variable_scales,
            physical_jacobian_perturbations=perturbations,
            jacobian_certificate=certificate,
            scaled_residuals=tuple(float(value) for value in scaled_residual),
            dimensionless_newton_correction=tuple(
                float(value) for value in normalized_correction
            ),
            physical_newton_correction=tuple(
                float(value) for value in physical_correction
            ),
            newton_corrected_scaled_residuals=tuple(
                float(value) for value in corrected_scaled_residual
            ),
            newton_corrected_maximum_scaled_residual=float(
                np.linalg.norm(corrected_scaled_residual, ord=np.inf)
            ),
            event_time_variable_index=event_time_index,
            linearized_event_time_correction_s=float(
                physical_correction[event_time_index]
            ),
            maximum_dimensionless_newton_correction=float(
                np.linalg.norm(normalized_correction, ord=np.inf)
            ),
        )
    except FaceEventTimeLinearizationAuditError:
        raise
    except Exception as exc:
        raise FaceEventTimeLinearizationAuditError(
            f"face-event time linearization audit failed closed: {exc}"
        ) from exc


def _face_event_physical_variable_scale_contract(
    step: FaceEventStep,
) -> tuple[tuple[str, ...], tuple[float, ...]]:
    """Return root-independent physical column scales for one accepted solve."""

    layout = step.assembly.layout
    chart = _coordinate_chart(step.before, step.controls)
    wet_temperature_width = chart.wet_temperature[1] - chart.wet_temperature[0]
    wet_water_width = chart.wet_water[1] - chart.wet_water[0]
    dry_temperature_width = chart.dry_temperature[1] - chart.dry_temperature[0]
    event_time_width = chart.event_time[1] - chart.event_time[0]
    interface_temperature_width = (
        chart.interface_temperature[1] - chart.interface_temperature[0]
    )
    if len(step.before.last_total_stefan_fluxes_mol_m2_s) != (
        layout.dry_face_count
    ):
        raise ValueError("before-state Stefan history lost face-event rank")
    stefan_scales = tuple(
        max(abs(value), _FACE_EVENT_STEFAN_SCALE_FLOOR_MOL_M2_S)
        for value in step.before.last_total_stefan_fluxes_mol_m2_s
    )
    labels = (
        *(
            f"wet_temperature[{index}]"
            for index in range(layout.wet_piece_count)
        ),
        *(
            f"wet_retained_water[{index}]"
            for index in range(layout.wet_piece_count)
        ),
        *(
            f"dry_temperature[{index}]"
            for index in range(layout.dry_piece_count)
        ),
        *(f"dry_y_hexane[{index}]" for index in range(layout.dry_piece_count)),
        *(
            f"total_stefan_flux[{index}]"
            for index in range(layout.dry_face_count)
        ),
        "event_time",
        "interface_temperature",
    )
    scales = (
        *((wet_temperature_width,) * layout.wet_piece_count),
        *((wet_water_width,) * layout.wet_piece_count),
        *((dry_temperature_width,) * layout.dry_piece_count),
        *((1.0,) * layout.dry_piece_count),
        *stefan_scales,
        event_time_width,
        interface_temperature_width,
    )
    if len(labels) != layout.unknown_count or len(scales) != layout.unknown_count:
        raise RuntimeError("face-event physical scale contract lost nonlinear rank")
    if not all(math.isfinite(value) and value > 0.0 for value in scales):
        raise ValueError("face-event physical scales must be positive and finite")
    return tuple(labels), tuple(float(value) for value in scales)


def _face_event_physical_residual(
    step: FaceEventStep,
    vector: np.ndarray,
) -> np.ndarray:
    layout = step.assembly.layout
    candidate = face.FaceArrivalUnknowns.from_vector(layout, vector)
    _encode_candidate(
        candidate,
        _coordinate_chart(step.before, step.controls),
        step.seed.stefan_flux_scales_mol_m2_s,
    )
    assembly = face.assemble_face_arrival(
        step.before.transport,
        candidate,
        step.boundary,
        enforce_reduced_film_thresholds=False,
    )
    return np.asarray(
        _datum_covariant_residual_vector(assembly.residuals),
        dtype=float,
    )


def _face_event_physical_jacobian(
    step: FaceEventStep,
    center: np.ndarray,
    variable_scales: np.ndarray,
) -> tuple[np.ndarray, tuple[float, ...]]:
    layout = step.assembly.layout
    matrix = np.empty((layout.residual_count, layout.unknown_count), dtype=float)
    perturbations = np.empty(layout.unknown_count, dtype=float)
    for column, characteristic in enumerate(variable_scales):
        delta = _FACE_EVENT_PHYSICAL_PROBE_RELATIVE_STEP * characteristic
        for _ in range(_FACE_EVENT_PHYSICAL_PROBE_MAX_HALVINGS):
            upper = center.copy()
            lower = center.copy()
            upper[column] += delta
            lower[column] -= delta
            try:
                upper_residual = _face_event_physical_residual(step, upper)
                lower_residual = _face_event_physical_residual(step, lower)
            except (ValueError, RuntimeError):
                delta *= 0.5
                continue
            matrix[:, column] = (
                upper_residual - lower_residual
            ) / (2.0 * delta)
            perturbations[column] = delta
            break
        else:
            raise ValueError(
                "face-event physical Jacobian column has no phase-safe central probe"
            )
    return matrix, tuple(float(value) for value in perturbations)


def _coordinate_chart(
    before: ci.CutIntegratorState,
    controls: FaceEventControls,
) -> _CoordinateChart:
    layout = face.layout_for_arrival(before.transport)
    dry_config = before.transport.config.dry
    return _CoordinateChart(
        wet_temperature=before.controls.wet_temperature_bounds_k,
        wet_water=before.controls.wet_water_bounds,
        dry_temperature=dry_config.conditioned_temperature_domain.solver_bounds_k,
        dry_y_hexane=(
            dry_config.primitive_band.y_hexane_bounds
            if isinstance(dry_config.primitive_band, ct.JointPrimitiveBand)
            else None
        ),
        dry_config=dry_config,
        dry_pores=tuple(
            replace(
                dry_config.pore,
                w_o=before.transport.oil_fraction_labels[cell],
            )
            for cell in layout.dry_cell_indices
        ),
        event_time=controls.event_time_bounds_s,
        interface_temperature=before.controls.interface_temperature_bounds_k,
    )


def _encode_candidate(
    candidate: face.FaceArrivalUnknowns,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
) -> np.ndarray:
    values = (
        *(
            ci._encode_open(value, chart.wet_temperature)  # noqa: SLF001
            for value in candidate.wet_temperatures_k
        ),
        *(
            ci._encode_open(value, chart.wet_water)  # noqa: SLF001
            for value in candidate.wet_retained_water_loadings
        ),
        *(
            ci._encode_open(value, chart.dry_temperature)  # noqa: SLF001
            for value in candidate.dry_temperatures_k
        ),
        *(
            _encode_dry_y(temperature, value, pore, chart)
            for temperature, value, pore in zip(
                candidate.dry_temperatures_k,
                candidate.dry_y_hexane,
                chart.dry_pores,
            )
        ),
        *(
            value / scale
            for value, scale in zip(
                candidate.dry_total_stefan_fluxes_mol_m2_s,
                stefan_scales,
            )
        ),
        ci._encode_open(candidate.event_time_s, chart.event_time),  # noqa: SLF001
        ci._encode_open(  # noqa: SLF001
            candidate.interface_temperature_k,
            chart.interface_temperature,
        ),
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("encoded face-event coordinates must be finite")
    return np.asarray(values, dtype=float)


def _decode_candidate(
    coordinates: Sequence[float],
    layout: face.FaceArrivalLayout,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
) -> face.FaceArrivalUnknowns:
    values = tuple(float(value) for value in coordinates)
    if len(values) != layout.unknown_count or not all(
        math.isfinite(value) for value in values
    ):
        raise ValueError("face-event coordinate vector is non-finite or wrong-rank")
    nw = layout.wet_piece_count
    nd = layout.dry_piece_count
    nf = layout.dry_face_count
    cursor = 0

    def take(count: int) -> tuple[float, ...]:
        nonlocal cursor
        result = values[cursor : cursor + count]
        cursor += count
        return result

    wet_t = tuple(
        ci._decode_open(value, chart.wet_temperature)  # noqa: SLF001
        for value in take(nw)
    )
    wet_w = tuple(
        ci._decode_open(value, chart.wet_water)  # noqa: SLF001
        for value in take(nw)
    )
    dry_t = tuple(
        ci._decode_open(value, chart.dry_temperature)  # noqa: SLF001
        for value in take(nd)
    )
    dry_y_coordinates = take(nd)
    dry_y = tuple(
        _decode_dry_y(temperature, value, pore, chart)
        for temperature, value, pore in zip(
            dry_t,
            dry_y_coordinates,
            chart.dry_pores,
        )
    )
    nt = tuple(value * scale for value, scale in zip(take(nf), stefan_scales))
    event_time = ci._decode_open(take(1)[0], chart.event_time)  # noqa: SLF001
    interface_t = ci._decode_open(  # noqa: SLF001
        take(1)[0],
        chart.interface_temperature,
    )
    return face.FaceArrivalUnknowns(
        wet_temperatures_k=wet_t,
        wet_retained_water_loadings=wet_w,
        dry_temperatures_k=dry_t,
        dry_y_hexane=dry_y,
        dry_total_stefan_fluxes_mol_m2_s=nt,
        event_time_s=event_time,
        interface_temperature_k=interface_t,
    )


def _encode_dry_y(
    temperature_k: float,
    y_hexane: float,
    pore: cp.CoupledPoreParams,
    chart: _CoordinateChart,
) -> float:
    if chart.dry_y_hexane is not None:
        return ci._encode_open(y_hexane, chart.dry_y_hexane)  # noqa: SLF001
    return cp.encode_gas_only_y(
        temperature_k,
        chart.dry_config.pressure_pa,
        y_hexane,
        pore,
    )


def _decode_dry_y(
    temperature_k: float,
    coordinate: float,
    pore: cp.CoupledPoreParams,
    chart: _CoordinateChart,
) -> float:
    if chart.dry_y_hexane is not None:
        return ci._decode_open(coordinate, chart.dry_y_hexane)  # noqa: SLF001
    return cp.decode_gas_only_y(
        temperature_k,
        chart.dry_config.pressure_pa,
        coordinate,
        pore,
    )


def _validate_event_chart(
    before: ci.CutIntegratorState,
    controls: FaceEventControls,
) -> None:
    ci.certify_explicit_wet_retained_water_applicability(
        before.controls,
        before.transport.wet_retained_water_loadings,
        before.transport.effective_wet_retained_water_capacity_duals_over_rt,
        before.transport.config.wet.luikov,
    )
    layout = face.layout_for_arrival(before.transport)
    k = layout.arrival_face_index
    for temperature in (
        before.controls.interface_temperature_bounds_k[0],
        0.5 * math.fsum(before.controls.interface_temperature_bounds_k),
        before.controls.interface_temperature_bounds_k[1],
    ):
        face.evaluate_face_interface_state(
            temperature,
            before.transport.config,
            before.transport.historical_hexane_loadings[k],
            before.transport.oil_fraction_labels[k],
            before.transport.oil_fraction_labels[k],
        )
    if controls.event_time_bounds_s[0] < 0.0:
        raise ValueError("face-event time chart cannot start before zero")


def _validate_seed(
    before: ci.CutIntegratorState,
    controls: FaceEventControls,
    seed: FaceEventSeed,
    layout: face.FaceArrivalLayout,
) -> None:
    if len(seed.candidate.vector()) != layout.unknown_count:
        raise ValueError("face-event seed does not match the event rank")
    if len(seed.stefan_flux_scales_mol_m2_s) != layout.dry_face_count:
        raise ValueError("one face-event Stefan scale is required per dry face")
    _encode_candidate(
        seed.candidate,
        _coordinate_chart(before, controls),
        seed.stefan_flux_scales_mol_m2_s,
    )


def _bounded_fractions(
    candidate: face.FaceArrivalUnknowns,
    chart: _CoordinateChart,
) -> tuple[float, ...]:
    pairs = [
        *((value, chart.wet_temperature) for value in candidate.wet_temperatures_k),
        *((value, chart.wet_water) for value in candidate.wet_retained_water_loadings),
        *((value, chart.dry_temperature) for value in candidate.dry_temperatures_k),
        (candidate.event_time_s, chart.event_time),
        (candidate.interface_temperature_k, chart.interface_temperature),
    ]
    result = [
        (value - bounds[0]) / (bounds[1] - bounds[0])
        for value, bounds in pairs
    ]
    if chart.dry_y_hexane is not None:
        result.extend(
            (value - chart.dry_y_hexane[0])
            / (chart.dry_y_hexane[1] - chart.dry_y_hexane[0])
            for value in candidate.dry_y_hexane
        )
    else:
        activity_at_ref = sp.water_activity(
            chart.dry_config.pore.luikov.W_ref,
            chart.dry_config.pore.luikov,
        )
        for temperature, composition, pore in zip(
            candidate.dry_temperatures_k,
            candidate.dry_y_hexane,
            chart.dry_pores,
        ):
            interval = cp.gas_only_composition_interval(
                temperature,
                chart.dry_config.pressure_pa,
                pore,
            )
            state = cp.evaluate_equilibrium(
                temperature,
                chart.dry_config.pressure_pa,
                composition,
                pore,
            )
            result.append(
                (composition - interval.lower_y_hexane)
                / (interval.upper_y_hexane - interval.lower_y_hexane)
            )
            for gap in (
                1.0 - state.water_activity,
                1.0 - state.hexane_activity,
                state.water_activity - activity_at_ref,
            ):
                result.append(gap / (1.0 + gap))
    result_tuple = tuple(result)
    if not all(0.0 < value < 1.0 for value in result_tuple):
        raise RuntimeError("accepted face-event primitive touched a chart boundary")
    return result_tuple


def _datum_covariant_residual_vector(
    residuals: face.FaceArrivalResidualBlocks,
) -> tuple[float, ...]:
    dry_energy = tuple(
        math.fsum((energy, -hx.H_REF * hexane))
        for energy, hexane in zip(
            residuals.dry_energy_w,
            residuals.dry_hexane_mol_s,
        )
    )
    rh_energy = math.fsum(
        (residuals.rh_energy_w, -hx.H_REF * residuals.rh_hexane_mol_s)
    )
    return (
        *residuals.wet_water_mol_s,
        *residuals.wet_energy_w,
        *residuals.dry_water_mol_s,
        *residuals.dry_hexane_mol_s,
        *dry_energy,
        residuals.rh_water_mol_s,
        residuals.rh_hexane_mol_s,
        rh_energy,
    )


def _residual_scales(
    before: cut.CutTransportState,
    initial: face.FaceArrivalAssembly,
) -> FaceEventResidualScales:
    tau = initial.candidate.event_time_s
    grid = before.geometry.master_grid
    layout = initial.layout
    k = layout.arrival_face_index
    floor = 1.0e-300
    covariant = _datum_covariant_residual_vector(initial.residuals)
    old_wet = initial.old_wet_cells
    old_dry = initial.old_dry_cells
    new_wet = initial.wet_cells
    new_dry = initial.dry_cells

    wet_water = tuple(
        max(
            grid.volumes[index]
            * old_cell.retained_water_concentration_mol_m3
            / tau,
            grid.volumes[index]
            * new_cell.retained_water_concentration_mol_m3
            / tau,
            abs(initial.residuals.wet_water_mol_s[position]),
            floor,
        )
        for position, (index, old_cell, new_cell) in enumerate(
            zip(layout.wet_cell_indices, old_wet, new_wet)
        )
    )
    wet_energy = tuple(
        max(
            _wet_capacity_energy(
                grid.volumes[index],
                old_cell,
                before,
            )
            / tau,
            _wet_capacity_energy(
                grid.volumes[index],
                new_cell,
                before,
            )
            / tau,
            abs(initial.residuals.wet_energy_w[position]),
            floor,
        )
        for position, (index, old_cell, new_cell) in enumerate(
            zip(layout.wet_cell_indices, old_wet, new_wet)
        )
    )

    old_cut = before.geometry.cells[k]
    old_wet_cut = old_wet[-1]
    old_dry_cut = old_dry[0]
    first_new = new_dry[0]
    old_wet_hexane = (
        before.config.wet.wet.rho_dm_p
        * before.historical_hexane_loadings[k]
        / hx.M
    )

    def first_component_scale(
        old_wet_concentration: float,
        old_dry_concentration: float,
        new_concentration: float,
        residual: float,
    ) -> float:
        return max(
            old_cut.wet_volume_m3 * old_wet_concentration / tau,
            old_cut.dry_volume_m3 * old_dry_concentration / tau,
            grid.volumes[k] * new_concentration / tau,
            abs(residual),
            floor,
        )

    dry_water_values = [
        first_component_scale(
            old_wet_cut.retained_water_concentration_mol_m3,
            old_dry_cut.total_water_concentration_mol_m3,
            first_new.total_water_concentration_mol_m3,
            initial.residuals.dry_water_mol_s[0],
        )
    ]
    dry_hexane_values = [
        first_component_scale(
            old_wet_hexane,
            old_dry_cut.total_hexane_concentration_mol_m3,
            first_new.total_hexane_concentration_mol_m3,
            initial.residuals.dry_hexane_mol_s[0],
        )
    ]
    first_old_wet_capacity = _wet_capacity_energy(
        old_cut.wet_volume_m3,
        old_wet_cut,
        before,
    )
    first_old_dry_capacity, first_old_dry_audit = (
        _dry_capacity_energy_with_audit(
            old_cut.dry_volume_m3,
            old_dry_cut,
            before,
            k,
        )
    )
    first_old_capacity = math.fsum(
        (first_old_wet_capacity, first_old_dry_capacity)
    )
    first_new_capacity, first_new_audit = _dry_capacity_energy_with_audit(
        grid.volumes[k],
        first_new,
        before,
        k,
    )
    dry_capacity_audits = [first_old_dry_audit, first_new_audit]
    dry_energy_values = [
        max(
            first_old_capacity / tau,
            first_new_capacity / tau,
            abs(
                math.fsum(
                    (
                        initial.residuals.dry_energy_w[0],
                        -hx.H_REF * initial.residuals.dry_hexane_mol_s[0],
                    )
                )
            ),
            floor,
        )
    ]
    for position, cell_index in enumerate(layout.dry_cell_indices[1:], start=1):
        old_cell = old_dry[position]
        new_cell = new_dry[position]
        volume = grid.volumes[cell_index]
        dry_water_values.append(
            max(
                volume * old_cell.total_water_concentration_mol_m3 / tau,
                volume * new_cell.total_water_concentration_mol_m3 / tau,
                abs(initial.residuals.dry_water_mol_s[position]),
                floor,
            )
        )
        dry_hexane_values.append(
            max(
                volume * old_cell.total_hexane_concentration_mol_m3 / tau,
                volume * new_cell.total_hexane_concentration_mol_m3 / tau,
                abs(initial.residuals.dry_hexane_mol_s[position]),
                floor,
            )
        )
        old_capacity, old_audit = _dry_capacity_energy_with_audit(
            volume,
            old_cell,
            before,
            cell_index,
        )
        new_capacity, new_audit = _dry_capacity_energy_with_audit(
            volume,
            new_cell,
            before,
            cell_index,
        )
        dry_capacity_audits.extend((old_audit, new_audit))
        dry_energy_values.append(
            max(
                old_capacity / tau,
                new_capacity / tau,
                abs(
                    math.fsum(
                        (
                            initial.residuals.dry_energy_w[position],
                            -hx.H_REF
                            * initial.residuals.dry_hexane_mol_s[position],
                        )
                    )
                ),
                floor,
            )
        )
    return FaceEventResidualScales(
        wet_water_mol_s=wet_water,
        wet_energy_w=wet_energy,
        dry_water_mol_s=tuple(dry_water_values),
        dry_hexane_mol_s=tuple(dry_hexane_values),
        dry_energy_w=tuple(dry_energy_values),
        rh_water_mol_s=max(
            old_cut.wet_volume_m3
            * old_wet_cut.retained_water_concentration_mol_m3
            / tau,
            old_cut.dry_volume_m3
            * old_dry_cut.total_water_concentration_mol_m3
            / tau,
            abs(initial.residuals.rh_water_mol_s),
            floor,
        ),
        rh_hexane_mol_s=max(
            old_cut.wet_volume_m3 * old_wet_hexane / tau,
            old_cut.dry_volume_m3
            * old_dry_cut.total_hexane_concentration_mol_m3
            / tau,
            abs(initial.residuals.rh_hexane_mol_s),
            floor,
        ),
        rh_energy_w=max(
            first_old_capacity / tau,
            abs(covariant[-1]),
            floor,
        ),
        dry_energy_capacity_stencil_audits=tuple(dry_capacity_audits),
    )


def _wet_capacity_energy(
    volume_m3: float,
    cell,
    state: cut.CutTransportState,
) -> float:
    capacity = wet_core.energy_capacity(
        cell.temperature_k,
        cell.historical_hexane_loading,
        replace(
            state.config.wet.wet,
            X_water=cell.retained_water_loading,
            w_o=cell.oil_fraction_label,
        ),
    )
    result = volume_m3 * capacity * max(abs(cell.temperature_k), 1.0)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("wet event capacity scale must be positive")
    return result


def _dry_capacity_energy(
    volume_m3: float,
    cell: cp.EquilibriumPoreState,
    state: cut.CutTransportState,
    cell_index: int,
) -> float:
    """Compatibility scalar view of the audited event capacity scale."""

    return _dry_capacity_energy_with_audit(
        volume_m3,
        cell,
        state,
        cell_index,
    )[0]


def _dry_capacity_energy_with_audit(
    volume_m3: float,
    cell: cp.EquilibriumPoreState,
    state: cut.CutTransportState,
    cell_index: int,
) -> tuple[float, ci.DryEnergyCapacityStencilAudit]:
    capacity, audit = ci._dry_energy_capacity_with_audit(  # noqa: SLF001
        cell.temperature_k,
        cell.y_hexane,
        replace(state.config.dry.pore, w_o=state.oil_fraction_labels[cell_index]),
        state.config,
    )
    result = volume_m3 * capacity * max(abs(cell.temperature_k), 1.0)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("dry event capacity scale must be positive")
    return result, audit


def _event_capacity_scale(assembly: face.FaceArrivalAssembly) -> float:
    state = assembly.before
    grid = state.geometry.master_grid
    values = [
        _wet_capacity_energy(grid.volumes[index], cell, state)
        for index, cell in zip(assembly.layout.wet_cell_indices, assembly.wet_cells)
    ]
    values.extend(
        _dry_capacity_energy(grid.volumes[index], cell, state, index)
        for index, cell in zip(assembly.layout.dry_cell_indices, assembly.dry_cells)
    )
    result = math.fsum(values)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("event whole-particle capacity scale must be positive")
    return result


def _event_inventory(
    assembly: face.FaceArrivalAssembly,
) -> ci.MaterialInventorySnapshot:
    grid = assembly.event_geometry.master_grid
    water = [0.0] * grid.n
    hexane = [0.0] * grid.n
    energy = [0.0] * grid.n
    for index, cell in zip(assembly.layout.wet_cell_indices, assembly.wet_cells):
        volume = grid.volumes[index]
        water[index] = volume * cell.retained_water_concentration_mol_m3
        hexane[index] = (
            volume
            * assembly.before.config.wet.wet.rho_dm_p
            * assembly.before.historical_hexane_loadings[index]
            / hx.M
        )
        energy[index] = volume * cell.energy_density_j_m3
    for index, cell in zip(assembly.layout.dry_cell_indices, assembly.dry_cells):
        volume = grid.volumes[index]
        water[index] = volume * cell.total_water_concentration_mol_m3
        hexane[index] = volume * cell.total_hexane_concentration_mol_m3
        energy[index] = volume * cell.energy_density_j_m3
    return ci.MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def _material_changes(
    assembly: face.FaceArrivalAssembly,
) -> ci.MaterialInventorySnapshot:
    n = assembly.event_geometry.master_grid.n
    water = [0.0] * n
    hexane = [0.0] * n
    energy = [0.0] * n
    for position, index in enumerate(assembly.layout.wet_cell_indices):
        water[index] = assembly.inventory_changes.wet_water_mol[position]
        energy[index] = assembly.inventory_changes.wet_energy_j[position]
    for position, index in enumerate(assembly.layout.dry_cell_indices):
        water[index] = assembly.inventory_changes.dry_water_mol[position]
        hexane[index] = assembly.inventory_changes.dry_hexane_mol[position]
        energy[index] = assembly.inventory_changes.dry_energy_j[position]
    return ci.MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def _accept_event(
    before: ci.CutIntegratorState,
    controls: FaceEventControls,
    boundary: ct.PoreBoundary,
    seed: FaceEventSeed,
    assembly: face.FaceArrivalAssembly,
    scales: FaceEventResidualScales,
    scaled_residuals: tuple[float, ...],
    evaluations: int,
    rejected_trials: int,
    nonlinear_message: str,
    condition_proxy: float,
    chart: _CoordinateChart,
    structural: _EventSparsity | None,
) -> FaceEventStep:
    tau = assembly.candidate.event_time_s
    if not math.isfinite(before.transport.time_s + tau):
        raise ValueError("accepted event time overflowed")
    current_inventory = _event_inventory(assembly)
    changes = _material_changes(assembly)
    area = assembly.event_geometry.master_grid.areas[-1]
    surface = assembly.dry_face_fluxes[-1]
    water_out = tau * area * surface.component.conserved_water_flux_mol_m2_s
    hexane_out = tau * area * surface.component.conserved_hexane_flux_mol_m2_s
    energy_out = tau * area * surface.energy.total_energy_flux_w_m2
    water_step = math.fsum((*changes.water_cell_mol, water_out))
    hexane_step = math.fsum((*changes.hexane_cell_mol, hexane_out))
    energy_step = math.fsum((*changes.energy_cell_j, energy_out))
    old_inventory = ci.inventory_snapshot(before.transport)
    after_capacity = _event_capacity_scale(assembly)
    step_water_scale = max(
        old_inventory.total_water_mol,
        current_inventory.total_water_mol,
        abs(water_out),
        1.0e-300,
    )
    step_hexane_scale = max(
        old_inventory.total_hexane_mol,
        current_inventory.total_hexane_mol,
        abs(hexane_out),
        1.0e-300,
    )
    step_energy_scale = max(
        ci.capacity_energy_scale(before.transport),
        after_capacity,
        1.0e-300,
    )

    boundary_water, boundary_water_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_water_out_mol,
        before.boundary_water_compensation_mol,
        water_out,
    )
    boundary_hexane, boundary_hexane_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_hexane_out_mol,
        before.boundary_hexane_compensation_mol,
        hexane_out,
    )
    boundary_energy, boundary_energy_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_energy_out_j,
        before.boundary_energy_compensation_j,
        energy_out,
    )
    cumulative_water, water_comp = ci._compensated_add_vectors(  # noqa: SLF001
        before.cumulative_material_water_change_cell_mol,
        before.water_change_compensation_cell_mol,
        changes.water_cell_mol,
    )
    cumulative_hexane, hexane_comp = ci._compensated_add_vectors(  # noqa: SLF001
        before.cumulative_material_hexane_change_cell_mol,
        before.hexane_change_compensation_cell_mol,
        changes.hexane_cell_mol,
    )
    cumulative_energy, energy_comp = ci._compensated_add_vectors(  # noqa: SLF001
        before.cumulative_material_energy_change_cell_j,
        before.energy_change_compensation_cell_j,
        changes.energy_cell_j,
    )
    corrected_boundary_water = math.fsum((boundary_water, boundary_water_comp))
    corrected_boundary_hexane = math.fsum(
        (boundary_hexane, boundary_hexane_comp)
    )
    cumulative_water_residual = math.fsum(
        (
            *cumulative_water,
            *water_comp,
            boundary_water,
            boundary_water_comp,
        )
    )
    cumulative_hexane_residual = math.fsum(
        (
            *cumulative_hexane,
            *hexane_comp,
            boundary_hexane,
            boundary_hexane_comp,
        )
    )
    cumulative_energy_residual = math.fsum(
        (
            *cumulative_energy,
            *energy_comp,
            boundary_energy,
            boundary_energy_comp,
        )
    )
    cumulative_water_scale = max(
        before.reference_inventory.total_water_mol,
        current_inventory.total_water_mol,
        abs(corrected_boundary_water),
        before.cumulative_absolute_water_transfer_mol + abs(water_out),
        1.0e-300,
    )
    cumulative_hexane_scale = max(
        before.reference_inventory.total_hexane_mol,
        current_inventory.total_hexane_mol,
        abs(corrected_boundary_hexane),
        before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out),
        1.0e-300,
    )
    cumulative_energy_scale = max(
        before.reference_capacity_energy_scale_j,
        after_capacity,
        before.cumulative_absolute_energy_transfer_j + abs(energy_out),
        1.0e-300,
    )
    candidate = assembly.candidate
    after = FaceArrivalCommittedState(
        geometry=assembly.event_geometry,
        config=before.transport.config,
        solver_controls=before.controls,
        time_s=before.transport.time_s + tau,
        arrival_face_index=assembly.layout.arrival_face_index,
        wet_temperatures_k=candidate.wet_temperatures_k,
        wet_retained_water_loadings=candidate.wet_retained_water_loadings,
        dry_temperatures_k=candidate.dry_temperatures_k,
        dry_y_hexane=candidate.dry_y_hexane,
        historical_hexane_loadings=before.transport.historical_hexane_loadings,
        oil_fraction_labels=before.transport.oil_fraction_labels,
        interface_temperature_k=candidate.interface_temperature_k,
        dry_total_stefan_fluxes_mol_m2_s=(
            candidate.dry_total_stefan_fluxes_mol_m2_s
        ),
        current_inventory=current_inventory,
        reference_inventory=before.reference_inventory,
        reference_capacity_energy_scale_j=before.reference_capacity_energy_scale_j,
        cumulative_boundary_water_out_mol=boundary_water,
        cumulative_boundary_hexane_out_mol=boundary_hexane,
        cumulative_boundary_energy_out_j=boundary_energy,
        boundary_water_compensation_mol=boundary_water_comp,
        boundary_hexane_compensation_mol=boundary_hexane_comp,
        boundary_energy_compensation_j=boundary_energy_comp,
        cumulative_absolute_water_transfer_mol=(
            before.cumulative_absolute_water_transfer_mol + abs(water_out)
        ),
        cumulative_absolute_hexane_transfer_mol=(
            before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out)
        ),
        cumulative_absolute_energy_transfer_j=(
            before.cumulative_absolute_energy_transfer_j + abs(energy_out)
        ),
        cumulative_material_water_change_cell_mol=cumulative_water,
        cumulative_material_hexane_change_cell_mol=cumulative_hexane,
        cumulative_material_energy_change_cell_j=cumulative_energy,
        water_change_compensation_cell_mol=water_comp,
        hexane_change_compensation_cell_mol=hexane_comp,
        energy_change_compensation_cell_j=energy_comp,
        accepted_steps=before.accepted_steps + 1,
        cumulative_nonlinear_evaluations=(
            before.cumulative_nonlinear_evaluations + evaluations
        ),
    )
    fractions = _bounded_fractions(candidate, chart)
    ledger = FaceEventLedger(
        event_time_s=tau,
        nonlinear_evaluations=evaluations,
        rejected_trial_evaluations=rejected_trials,
        nonlinear_message=nonlinear_message,
        residual_scales=scales,
        scaled_residuals=scaled_residuals,
        maximum_scaled_residual=max(abs(value) for value in scaled_residuals),
        condition_proxy=condition_proxy,
        minimum_fractional_distance_to_bound=min(
            min(value, 1.0 - value) for value in fractions
        ),
        used_sparse_jacobian=structural is not None,
        scalability_class=(
            "bordered-banded grouped finite-difference sparse event Jacobian; "
            "runtime qualification remains open"
            if structural is not None
            else "certified small-N dense finite-difference event oracle"
        ),
        jacobian_structural_color_count=(
            assembly.layout.unknown_count
            if structural is None
            else structural.color_count
        ),
        jacobian_structural_nnz=(
            assembly.layout.unknown_count**2
            if structural is None
            else structural.nnz
        ),
        water_step_residual_mol=water_step,
        hexane_step_residual_mol=hexane_step,
        energy_step_residual_j=energy_step,
        normalized_water_step_residual=abs(water_step) / step_water_scale,
        normalized_hexane_step_residual=abs(hexane_step) / step_hexane_scale,
        normalized_energy_step_residual=abs(energy_step) / step_energy_scale,
        water_cumulative_residual_mol=cumulative_water_residual,
        hexane_cumulative_residual_mol=cumulative_hexane_residual,
        energy_cumulative_residual_j=cumulative_energy_residual,
        normalized_water_cumulative_residual=(
            abs(cumulative_water_residual) / cumulative_water_scale
        ),
        normalized_hexane_cumulative_residual=(
            abs(cumulative_hexane_residual) / cumulative_hexane_scale
        ),
        normalized_energy_cumulative_residual=(
            abs(cumulative_energy_residual) / cumulative_energy_scale
        ),
        boundary_water_out_mol=water_out,
        boundary_hexane_out_mol=hexane_out,
        boundary_energy_out_j=energy_out,
        accepted=True,
    )
    return FaceEventStep(before, after, assembly, boundary, controls, seed, ledger)


def _event_sparsity(layout: face.FaceArrivalLayout) -> _EventSparsity:
    unknown_spans = _spans(layout.unknown_blocks)
    residual_spans = _spans(layout.residual_blocks)
    unknown = {name: (start, size) for name, start, size in unknown_spans}
    residual = {name: (start, size) for name, start, size in residual_spans}
    rows = [set() for _ in range(layout.residual_count)]

    def column(name: str, index: int = 0) -> int:
        start, size = unknown[name]
        if not 0 <= index < size:
            raise IndexError(name)
        return start + index

    for equation_name in ("wet_water", "wet_energy"):
        row_start, _ = residual[equation_name]
        for piece in range(layout.wet_piece_count):
            for neighbour in range(
                max(0, piece - 1),
                min(layout.wet_piece_count, piece + 2),
            ):
                rows[row_start + piece].add(column("wet_temperature", neighbour))
                rows[row_start + piece].add(
                    column("wet_retained_water", neighbour)
                )
    for equation_name in ("dry_water", "dry_hexane", "dry_energy"):
        row_start, _ = residual[equation_name]
        for piece in range(layout.dry_piece_count):
            row = rows[row_start + piece]
            for neighbour in range(
                max(0, piece - 1),
                min(layout.dry_piece_count, piece + 2),
            ):
                row.add(column("dry_temperature", neighbour))
                row.add(column("dry_y_hexane", neighbour))
            row.add(column("dry_total_stefan_flux", piece))
            row.add(column("dry_total_stefan_flux", piece + 1))
            if piece == 0:
                row.add(column("wet_temperature", layout.wet_piece_count - 1))
                row.add(
                    column("wet_retained_water", layout.wet_piece_count - 1)
                )
    interface_columns = (
        column("wet_temperature", layout.wet_piece_count - 1),
        column("wet_retained_water", layout.wet_piece_count - 1),
        column("dry_temperature", 0),
        column("dry_y_hexane", 0),
        column("dry_total_stefan_flux", 0),
    )
    for name in ("rh_water", "rh_hexane", "rh_energy"):
        rows[residual[name][0]].update(interface_columns)
    borders = (column("event_time"), column("interface_temperature"))
    for row in rows:
        row.update(borders)
    row_columns = tuple(tuple(sorted(row)) for row in rows)
    indptr = [0]
    indices: list[int] = []
    for row in row_columns:
        indices.extend(row)
        indptr.append(len(indices))
    matrix = sparse.csr_matrix(
        (
            np.ones(len(indices), dtype=bool),
            np.asarray(indices, dtype=np.int64),
            np.asarray(indptr, dtype=np.int64),
        ),
        shape=(layout.residual_count, layout.unknown_count),
    )
    colors = _greedy_colors(layout.unknown_count, row_columns)
    return _EventSparsity(matrix, max(colors, default=-1) + 1)


def _spans(
    blocks: tuple[tuple[str, int], ...],
) -> tuple[tuple[str, int, int], ...]:
    start = 0
    result = []
    for name, size in blocks:
        result.append((name, start, size))
        start += size
    return tuple(result)


def _greedy_colors(
    count: int,
    rows: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    neighbours = [set() for _ in range(count)]
    for columns in rows:
        for position, column in enumerate(columns):
            neighbours[column].update(columns[:position])
            neighbours[column].update(columns[position + 1 :])
    order = sorted(range(count), key=lambda value: (-len(neighbours[value]), value))
    colors = [-1] * count
    for column in order:
        forbidden = {
            colors[neighbour]
            for neighbour in neighbours[column]
            if colors[neighbour] >= 0
        }
        color = 0
        while color in forbidden:
            color += 1
        colors[column] = color
    return tuple(colors)


__all__ = [
    "FaceArrivalCommittedState",
    "FaceDepartureChartRequired",
    "FaceEventControls",
    "FaceEventIntegratorStepError",
    "FaceEventLedger",
    "FaceEventResidualScales",
    "FaceEventSeed",
    "FaceEventStep",
    "FaceEventTimeLinearizationAudit",
    "FaceEventTimeLinearizationAuditError",
    "audit_face_event_time_linearization",
    "solve_inner_face_arrival",
]
