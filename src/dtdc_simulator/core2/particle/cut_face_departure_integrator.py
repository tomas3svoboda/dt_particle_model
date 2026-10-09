r"""Safeguarded fixed-dt nonlinear integrator for inner-face departure.

The input is the exact-face :class:`FaceArrivalCommittedState` produced by the
safeguarded face-arrival solver.  It is adapted only through
``adapt_face_arrival_committed_state``; no epsilon strict-cut restart exists.
The fixed-dt corrector uses a coalescence-regular chart: existing material
cells are represented by rates, newborn traces by their one-sided gradients,
and the moving geometry by the exact positive newborn-volume rate ``Q``.
The ordinary full ALE/RH residual is still assembled and accepted unchanged.

Gas-only barriers are feasibility guards only; they are never inverted to
recover a newborn gradient.  The direct ``g_T`` and ``g_y`` coordinates remain
the authoritative traces passed to assembly.  The newborn composition is
checked at its storage temperature, moving ``T_Gamma`` face, and first
symmetric internal face.  Existing dry compositions are checked at every
storage/face context used by transport.

The requested duration is solved directly first.  Only a genuine basin miss
activates uncommitted same-old-state auxiliary durations; those numerical
roots must satisfy the unchanged ordinary residual and phase contracts, while
only the requested-duration root is conditioning-certified and reaches the
commit machinery.  The fallback start is a non-physical representability/
basin heuristic and never changes the residual or accepted duration.  The
same safeguarded direct Newton uses dense central differences below the
declared dense-size limit and structurally verified colored central
differences with sparse LU above it.  The bordered-radial nonlinear Jacobian
has O(N) storage and a bounded color count; the mandatory final physical
``J_star`` SVD remains dense and is reported honestly as the remaining
O(N^2)-memory/O(N^3)-factorization seam.  There is no penalty residual,
clipping, tolerance relaxation, contact/stall fallback, or physics-changing
predictor.

On acceptance the exact-face commit's original reference inventories,
capacity scale, signed Neumaier accumulators, absolute-transfer histories, and
counters are advanced in place semantically (through a new immutable state),
never reinitialized.  The output is an ordinary strict-partial
:class:`cut_integrator.CutIntegratorState` ready for same-cell integration.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Sequence

import numpy as np
from scipy import sparse as scipy_sparse
from scipy.sparse import linalg as sparse_linalg

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import conditioning
from dtdc_simulator.core2.particle import cut_face_departure as departure
from dtdc_simulator.core2.particle import cut_face_event_integrator as face_solver
from dtdc_simulator.core2.particle import cut_face_tangent as face_tangent
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_sparsity
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.props import hexane as hx


class FaceDepartureIntegratorStepError(RuntimeError):
    """Rejected departure solve with exact committed-state rollback."""

    def __init__(
        self,
        message: str,
        rollback_state: face_solver.FaceArrivalCommittedState,
        *,
        nonlinear_evaluations: int,
        rejected_trial_evaluations: int,
        maximum_scaled_residual: float = math.inf,
        condition_proxy: float = math.inf,
        raw_solver_condition_proxy: float = math.inf,
        last_candidate: departure.FaceDepartureUnknowns | None = None,
        scaled_residuals: tuple[float, ...] = (),
        residual_blocks: cut.CutResidualBlocks | None = None,
        minimum_fractional_distance_to_bound: float = 0.0,
        nonlinear_method: str = "not-started",
        retained_water_applicability_guard_failed: bool = False,
        retained_water_applicability_audit: (
            ci.WetRetainedWaterApplicabilityAudit | None
        ) = None,
        direct_newton_iteration_audits: tuple[
            "_DirectNewtonIterationAudit", ...
        ] = (),
        cauchy_globalization_accepted_steps: int = 0,
        cauchy_line_search_trials: int = 0,
        cauchy_residual_evaluations: int = 0,
        dogleg_globalization_accepted_steps: int = 0,
        dogleg_jacobian_linearizations: int = 0,
        dogleg_trial_attempts: int = 0,
        dogleg_residual_evaluations: int = 0,
        direct_dogleg_iteration_audits: tuple[
            "_DirectDoglegIterationAudit", ...
        ] = (),
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.nonlinear_evaluations = nonlinear_evaluations
        self.rejected_trial_evaluations = rejected_trial_evaluations
        self.maximum_scaled_residual = maximum_scaled_residual
        self.condition_proxy = condition_proxy
        self.raw_solver_condition_proxy = raw_solver_condition_proxy
        self.last_candidate = last_candidate
        self.scaled_residuals = scaled_residuals
        self.residual_blocks = residual_blocks
        self.minimum_fractional_distance_to_bound = (
            minimum_fractional_distance_to_bound
        )
        self.nonlinear_method = nonlinear_method
        self.retained_water_applicability_guard_failed = (
            retained_water_applicability_guard_failed
        )
        self.retained_water_applicability_audit = (
            retained_water_applicability_audit
        )
        self.direct_newton_iteration_audits = direct_newton_iteration_audits
        self.cauchy_globalization_accepted_steps = (
            cauchy_globalization_accepted_steps
        )
        self.cauchy_line_search_trials = cauchy_line_search_trials
        self.cauchy_residual_evaluations = cauchy_residual_evaluations
        self.dogleg_globalization_accepted_steps = (
            dogleg_globalization_accepted_steps
        )
        self.dogleg_jacobian_linearizations = dogleg_jacobian_linearizations
        self.dogleg_trial_attempts = dogleg_trial_attempts
        self.dogleg_residual_evaluations = dogleg_residual_evaluations
        self.direct_dogleg_iteration_audits = direct_dogleg_iteration_audits


class _EvaluationBudgetError(Exception):
    """A direct stage exhausted its declared residual-evaluation budget."""


class _DepartureConditioningError(Exception):
    """A converged requested-duration root failed its fixed certificate."""


class _DirectStageBasinError(Exception):
    """A valid direct stage could not take a safeguarded Newton step."""


@dataclass(frozen=True)
class FaceDepartureSeed:
    """One explicit transformed-coordinate seed and Stefan-flux scales."""

    candidate: departure.FaceDepartureUnknowns
    stefan_flux_scales_mol_m2_s: tuple[float, ...]
    label: str
    newborn_trace_gradients: departure.NewbornDryTraceGradients | None = None
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.candidate, departure.FaceDepartureUnknowns):
            raise TypeError("departure seed requires FaceDepartureUnknowns")
        if not self.stefan_flux_scales_mol_m2_s or not all(
            math.isfinite(value) and value > 0.0
            for value in self.stefan_flux_scales_mol_m2_s
        ):
            raise ValueError("every departure Stefan face needs a positive scale")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("departure seed needs non-empty provenance")
        if self.newborn_trace_gradients is not None and not isinstance(
            self.newborn_trace_gradients,
            departure.NewbornDryTraceGradients,
        ):
            raise TypeError("departure seed trace gradients have the wrong type")


@dataclass(frozen=True)
class FaceDepartureIntegratorLedger:
    """Nonlinear, conditioning, conservation, and handoff diagnostics."""

    dt_s: float
    nonlinear_method: str
    nonlinear_evaluations: int
    rejected_trial_evaluations: int
    nonlinear_message: str
    residual_scales: ci.CutResidualScales
    scaled_residuals: tuple[float, ...]
    maximum_scaled_residual: float
    # ``condition_proxy`` is the unit-invariant physical certificate.  The
    # raw solver-coordinate number is retained only as a numerical diagnostic.
    condition_proxy: float
    raw_solver_condition_proxy: float
    physical_variable_scale_provenance: str
    physical_variable_scale_labels: tuple[str, ...]
    physical_variable_scales: tuple[float, ...]
    physical_jacobian_perturbations: tuple[float, ...]
    jacobian_certificate: conditioning.UnitInvariantJacobianCertificate
    # This certificate belongs only to the accepted requested-duration root.
    # Same-old-state continuation roots never enter this ledger or its
    # conservation update; their separately recorded diagnostics below make no
    # accepted-state conditioning claim.
    accepted_state_conditioning_certified: bool = field(default=True, init=False)
    conditioning_certificate_scope: str = field(
        default="accepted requested-duration root only",
        init=False,
    )
    minimum_fractional_distance_to_bound: float
    minimum_front_z_distance_to_chart_boundary: float
    used_sparse_jacobian: bool
    scalability_class: str
    jacobian_structural_color_count: int
    jacobian_structural_nnz: int
    jacobian_structural_storage_bytes: int
    jacobian_dense_equivalent_bytes: int
    continuation_fallback_used: bool
    direct_requested_attempt_evaluations: int
    direct_requested_failure_message: str | None
    direct_requested_failure_maximum_scaled_residual: float | None
    auxiliary_stages: tuple["FaceDepartureAuxiliaryStage", ...]
    cauchy_globalization_used: bool
    cauchy_globalization_accepted_steps: int
    cauchy_line_search_trials: int
    cauchy_residual_evaluations: int
    dogleg_globalization_used: bool
    dogleg_globalization_accepted_steps: int
    dogleg_jacobian_linearizations: int
    dogleg_trial_attempts: int
    dogleg_residual_evaluations: int
    newborn_trace_gradients: departure.NewbornDryTraceGradients
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
    source_accepted_steps: int
    source_cumulative_nonlinear_evaluations: int
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
class FaceDepartureIntegratorStep:
    """One certified face departure and strict-cut handoff."""

    before: face_solver.FaceArrivalCommittedState
    after: ci.CutIntegratorState
    assembly: departure.FaceDepartureAssembly
    boundary: ct.PoreBoundary
    controls: ci.CutSolverControls
    seed: FaceDepartureSeed
    ledger: FaceDepartureIntegratorLedger
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceDepartureAuxiliaryStage:
    """Residual/phase-safe same-old-state root deliberately not committed.

    Auxiliary roots only seed the next solve from the identical old state.
    They receive neither an accepted-state condition certificate nor any
    trajectory/conservation commit; all of their residual evaluations remain
    charged to the public solve budget.
    """

    dt_s: float
    nonlinear_evaluations: int
    rejected_trial_evaluations: int
    maximum_scaled_residual: float
    condition_proxy: float | None
    raw_solver_condition_proxy: float | None
    cauchy_globalization_accepted_steps: int
    cauchy_line_search_trials: int
    cauchy_residual_evaluations: int
    dogleg_globalization_accepted_steps: int
    dogleg_jacobian_linearizations: int
    dogleg_trial_attempts: int
    dogleg_residual_evaluations: int
    committed: bool = field(default=False, init=False)
    entered_accepted_trajectory: bool = field(default=False, init=False)
    entered_conservation_ledger: bool = field(default=False, init=False)
    accepted_state_conditioning_certified: bool = field(default=False, init=False)
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceDepartureSparsityAudit:
    """Independent columnwise verification of the colored departure pattern."""

    unknown_count: int
    structural_nnz: int
    color_count: int
    residual_evaluations: int
    maximum_relative_derivative_outside_pattern: float
    maximum_relative_colored_difference_inside_pattern: float
    independent_condition_proxy: float
    accepted_condition_proxy: float
    relative_condition_proxy_difference: float
    passed: bool
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class _CoordinateChart:
    dt_s: float
    wet_temperature: tuple[float, float]
    wet_water: tuple[float, float]
    dry_temperature: tuple[float, float]
    dry_y_hexane: tuple[float, float] | None
    front_z: tuple[float, float]
    interface_temperature: tuple[float, float]
    reference_wet_temperatures_k: tuple[float, ...]
    reference_wet_water_loadings: tuple[float, ...]
    reference_outer_dry_temperatures_k: tuple[float, ...]
    reference_outer_dry_y_hexane: tuple[float, ...]
    dry_pores: tuple[cp.CoupledPoreParams, ...]
    face_conditioned_chart: ci._CoordinateChart  # noqa: SLF001
    newborn_temperature_gradient_bounds_k_m: tuple[float, float]
    newborn_composition_gradient_bounds_m_inv: tuple[float, float]
    outer_y_rate_bounds_s_inv: tuple[tuple[float, float], ...]
    q_scale_m3_s: float
    q_ratio_bounds: tuple[float, float]


@dataclass(frozen=True)
class _DecodedTrial:
    candidate: departure.FaceDepartureUnknowns
    newborn_trace_gradients: departure.NewbornDryTraceGradients


@dataclass(frozen=True)
class _DirectStageSolution:
    assembly: departure.FaceDepartureAssembly
    newborn_trace_gradients: departure.NewbornDryTraceGradients
    scales: ci.CutResidualScales
    signed_scaled_residuals: tuple[float, ...]
    evaluations: int
    rejected_trials: int
    message: str
    raw_condition_proxy: float
    condition_certificate: conditioning.UnitInvariantJacobianCertificate | None
    physical_variable_scale_labels: tuple[str, ...]
    physical_variable_scales: tuple[float, ...]
    physical_jacobian_perturbations: tuple[float, ...]
    chart: _CoordinateChart
    structural: cut_sparsity.CutJacobianSparsity | None
    cauchy_globalization_accepted_steps: int
    cauchy_line_search_trials: int
    cauchy_residual_evaluations: int
    dogleg_globalization_accepted_steps: int
    dogleg_jacobian_linearizations: int
    dogleg_trial_attempts: int
    dogleg_residual_evaluations: int


@dataclass(frozen=True)
class _DirectNewtonTrialAudit:
    """One behavior-neutral observation of a Newton line-search trial."""

    fraction: float
    rejection_reason: str | None
    residual_l2_norm: float | None
    maximum_scaled_residual: float | None
    strictly_decreased_l2_merit: bool


@dataclass(frozen=True)
class _DirectNewtonIterationAudit:
    """Behavior-neutral diagnostics for one departure Newton iteration."""

    iteration: int
    current_l2_norm: float
    current_maximum_scaled_residual: float
    jacobian_route: str
    minimum_coordinate_perturbation: float
    maximum_coordinate_perturbation: float
    minimum_relative_coordinate_perturbation: float
    maximum_relative_coordinate_perturbation: float
    predicted_merit_directional_derivative: float
    ideal_newton_directional_derivative: float
    scaled_correction_l2_norm: float
    maximum_open_bound_fraction: float
    line_search_trials: tuple[_DirectNewtonTrialAudit, ...]
    accepted_fraction: float | None
    accepted_route: str | None
    cauchy_gradient_l2_norm: float | None
    cauchy_numerator: float | None
    cauchy_denominator: float | None
    cauchy_alpha: float | None
    cauchy_scaled_correction_l2_norm: float | None
    cauchy_maximum_open_bound_fraction: float | None
    cauchy_line_search_trials: tuple[_DirectNewtonTrialAudit, ...]
    cauchy_accepted_fraction: float | None
    equations_changed: bool = field(default=False, init=False)


@dataclass(frozen=True)
class _DirectDoglegTrialAudit:
    """One phase-safe trust-radius trial in the diagnostic rescue."""

    trial: int
    trust_radius_before: float
    dogleg_branch: str
    scaled_step_l2_norm: float
    predicted_merit_reduction: float
    rejection_reason: str | None
    residual_l2_norm: float | None
    maximum_scaled_residual: float | None
    actual_merit_reduction: float | None
    trust_ratio: float | None
    accepted: bool
    trust_radius_after: float
    equations_changed: bool = field(default=False, init=False)


@dataclass(frozen=True)
class _DirectDoglegIterationAudit:
    """One Jacobian linearization in the diagnostic Powell-dogleg phase."""

    iteration: int
    current_l2_norm: float
    current_maximum_scaled_residual: float
    scaled_newton_correction_l2_norm: float
    scaled_cauchy_correction_l2_norm: float
    cauchy_gradient_l2_norm: float
    cauchy_numerator: float
    cauchy_denominator: float
    cauchy_alpha: float
    trust_radius_on_entry: float
    trials: tuple[_DirectDoglegTrialAudit, ...]
    accepted: bool
    trust_radius_on_exit: float
    equations_changed: bool = field(default=False, init=False)


@dataclass
class _DirectStageProgress:
    evaluations: int = 0
    rejected_trials: int = 0
    maximum_scaled_residual: float = math.inf
    condition_proxy: float = math.inf
    raw_condition_proxy: float = math.inf
    candidate: departure.FaceDepartureUnknowns | None = None
    signed_scaled_residuals: tuple[float, ...] = ()
    residual_blocks: cut.CutResidualBlocks | None = None
    current_candidate: departure.FaceDepartureUnknowns | None = None
    current_newborn_trace_gradients: (
        departure.NewbornDryTraceGradients | None
    ) = None
    current_solver_coordinates: tuple[float, ...] = ()
    direct_newton_iteration_audits: list[_DirectNewtonIterationAudit] = field(
        default_factory=list
    )
    cauchy_globalization_accepted_steps: int = 0
    cauchy_line_search_trials: int = 0
    cauchy_residual_evaluations: int = 0
    dogleg_globalization_accepted_steps: int = 0
    dogleg_jacobian_linearizations: int = 0
    dogleg_trial_attempts: int = 0
    dogleg_residual_evaluations: int = 0
    direct_dogleg_iteration_audits: list[_DirectDoglegIterationAudit] = field(
        default_factory=list
    )


@dataclass(frozen=True)
class _PhysicalScaleReference:
    tangent_temperature_gradient_k_m: float
    tangent_composition_gradient_m_inv: float
    tangent_stefan_fluxes_mol_m2_s: tuple[float, ...]
    tangent_q_m3_s: float
    tangent_interface_temperature_k: float
    newborn_temperature_gradient_scale_k_m: float
    newborn_composition_gradient_scale_m_inv: float
    stefan_flux_scales_mol_m2_s: tuple[float, ...]
    q_scale_m3_s: float


_STEFAN_SCALE_FLOOR_MOL_M2_S = 1.0e-3
_GRADIENT_SCALE_FLOOR = 1.0
_Q_SCALE_FLOOR_M3_S = 1.0e-18
_JACOBIAN_RANK_RELATIVE_TOLERANCE = 1.0e-13
_FALLBACK_REPRESENTABILITY_FLOOR_S = 1.0e-6
_SPARSITY_RELATIVE_ZERO_TOLERANCE = 1.0e-8


def _scaled_cauchy_correction(
    scaled_jacobian: scipy_sparse.spmatrix | np.ndarray,
    scaled_residual: np.ndarray,
) -> tuple[np.ndarray, float, float, float, float]:
    """Return the exact linearized least-squares Cauchy correction.

    This is a globalization direction only.  It uses the already assembled
    scaled residual/Jacobian and changes neither equations nor tolerances.
    """

    residual = np.asarray(scaled_residual, dtype=float)
    if residual.ndim != 1 or not np.all(np.isfinite(residual)):
        raise _DirectStageBasinError(
            "departure scaled Cauchy residual is non-finite or wrong-rank"
        )
    gradient = np.asarray(scaled_jacobian.T @ residual, dtype=float).reshape(-1)
    if gradient.shape != residual.shape or not np.all(np.isfinite(gradient)):
        raise _DirectStageBasinError(
            "departure scaled Cauchy gradient is non-finite or wrong-rank"
        )
    jacobian_gradient = np.asarray(
        scaled_jacobian @ gradient,
        dtype=float,
    ).reshape(-1)
    if jacobian_gradient.shape != residual.shape or not np.all(
        np.isfinite(jacobian_gradient)
    ):
        raise _DirectStageBasinError(
            "departure scaled Cauchy Jacobian action is non-finite or wrong-rank"
        )
    numerator = float(np.dot(gradient, gradient))
    denominator = float(np.dot(jacobian_gradient, jacobian_gradient))
    if (
        not math.isfinite(numerator)
        or not math.isfinite(denominator)
        or numerator <= 0.0
        or denominator <= 0.0
    ):
        raise _DirectStageBasinError(
            "departure scaled Cauchy direction has no finite positive curvature"
        )
    alpha = numerator / denominator
    correction = -alpha * gradient
    if not math.isfinite(alpha) or alpha <= 0.0 or not np.all(
        np.isfinite(correction)
    ):
        raise _DirectStageBasinError(
            "departure scaled Cauchy correction is non-finite"
        )
    return (
        correction,
        float(np.linalg.norm(gradient)),
        numerator,
        denominator,
        alpha,
    )


def _scaled_newton_correction(
    scaled_jacobian: scipy_sparse.spmatrix | np.ndarray,
    scaled_residual: np.ndarray,
) -> np.ndarray:
    """Solve the existing scaled Newton system without globalization policy."""

    try:
        correction = (
            sparse_linalg.splu(scaled_jacobian.tocsc()).solve(-scaled_residual)
            if scipy_sparse.issparse(scaled_jacobian)
            else np.linalg.solve(scaled_jacobian, -scaled_residual)
        )
    except (np.linalg.LinAlgError, RuntimeError) as exc:
        raise _DirectStageBasinError(
            "departure direct Newton Jacobian is singular"
        ) from exc
    result = np.asarray(correction, dtype=float)
    if result.shape != scaled_residual.shape or not np.all(np.isfinite(result)):
        raise _DirectStageBasinError(
            "departure direct Newton correction is non-finite or wrong-rank"
        )
    return result


def _scaled_powell_dogleg_correction(
    scaled_newton_correction: np.ndarray,
    scaled_cauchy_correction: np.ndarray,
    trust_radius: float,
) -> tuple[np.ndarray, str]:
    """Return the classical Powell dogleg point at one scaled trust radius."""

    newton = np.asarray(scaled_newton_correction, dtype=float)
    cauchy = np.asarray(scaled_cauchy_correction, dtype=float)
    if (
        newton.ndim != 1
        or cauchy.shape != newton.shape
        or not np.all(np.isfinite(newton))
        or not np.all(np.isfinite(cauchy))
        or not math.isfinite(trust_radius)
        or trust_radius <= 0.0
    ):
        raise _DirectStageBasinError(
            "departure scaled dogleg inputs are non-finite or invalid"
        )
    newton_norm = float(np.linalg.norm(newton))
    cauchy_norm = float(np.linalg.norm(cauchy))
    if not math.isfinite(newton_norm) or not math.isfinite(cauchy_norm):
        raise _DirectStageBasinError(
            "departure scaled dogleg correction norm is non-finite"
        )
    if newton_norm <= trust_radius:
        return newton.copy(), "newton_inside_trust_region"
    if cauchy_norm >= trust_radius:
        if cauchy_norm <= 0.0:
            raise _DirectStageBasinError(
                "departure scaled dogleg has a zero Cauchy vertex"
            )
        return (
            (trust_radius / cauchy_norm) * cauchy,
            "scaled_gradient_boundary",
        )
    difference = newton - cauchy
    quadratic = float(np.dot(difference, difference))
    linear = 2.0 * float(np.dot(cauchy, difference))
    constant = float(np.dot(cauchy, cauchy)) - trust_radius * trust_radius
    discriminant = linear * linear - 4.0 * quadratic * constant
    if (
        not math.isfinite(quadratic)
        or quadratic <= 0.0
        or not math.isfinite(discriminant)
        or discriminant < 0.0
    ):
        raise _DirectStageBasinError(
            "departure scaled dogleg segment intersection is invalid"
        )
    tau = (-linear + math.sqrt(discriminant)) / (2.0 * quadratic)
    if not math.isfinite(tau) or not 0.0 < tau < 1.0:
        raise _DirectStageBasinError(
            "departure scaled dogleg segment parameter left (0, 1)"
        )
    result = cauchy + tau * difference
    if not np.all(np.isfinite(result)):
        raise _DirectStageBasinError(
            "departure scaled dogleg correction is non-finite"
        )
    return result, "cauchy_to_newton_segment"


def _updated_dogleg_trust_radius(
    trust_radius: float,
    scaled_step_l2_norm: float,
    trust_ratio: float,
) -> float:
    """Apply the fixed Powell-dogleg radius update thresholds."""

    values = (trust_radius, scaled_step_l2_norm, trust_ratio)
    if (
        not all(math.isfinite(value) for value in values)
        or trust_radius <= 0.0
        or scaled_step_l2_norm < 0.0
    ):
        raise _DirectStageBasinError(
            "departure diagnostic dogleg radius update is invalid"
        )
    if trust_ratio < 0.25:
        return 0.25 * trust_radius
    if trust_ratio > 0.75 and scaled_step_l2_norm >= 0.9 * trust_radius:
        return 2.0 * trust_radius
    return trust_radius


def _exact_tangent_scale_reference(
    before: face_solver.FaceArrivalCommittedState,
) -> _PhysicalScaleReference:
    """Derive all seam scales from the unique exact-face tangent, not a seed."""

    tangent_candidate = face_tangent.candidate_from_committed_face(before)
    tangent_assembly = face_tangent.assemble_face_tangent(
        before,
        tangent_candidate,
    )
    k = before.arrival_face_index - 1
    pore = replace(before.config.dry.pore, w_o=before.oil_fraction_labels[k])
    derivative = cp.isothermal_potential_composition_derivatives(
        tangent_candidate.interface_temperature_k,
        before.config.dry.pressure_pa,
        tangent_assembly.interface.y_hexane,
        pore,
    ).gas_exchange_potential_derivative_per_y_hexane
    composition_gradient = (
        tangent_candidate.dry_gas_force_gradient_m_inv / derivative
    )
    raw_stefan = (
        tangent_assembly.direct_interface_flux.recovered_total_stefan_flux_mol_m2_s,
        tangent_candidate.outer_total_stefan_flux_mol_m2_s,
        *before.dry_total_stefan_fluxes_mol_m2_s[1:],
    )
    stefan_scales = tuple(
        max(abs(value), _STEFAN_SCALE_FLOOR_MOL_M2_S)
        for value in raw_stefan
    )
    q_scale = max(
        tangent_assembly.ledger.interface_area_m2
        * abs(tangent_candidate.front_speed_m_s),
        _Q_SCALE_FLOOR_M3_S,
    )
    return _PhysicalScaleReference(
        tangent_temperature_gradient_k_m=(
            tangent_candidate.dry_temperature_gradient_k_m
        ),
        tangent_composition_gradient_m_inv=composition_gradient,
        tangent_stefan_fluxes_mol_m2_s=raw_stefan,
        tangent_q_m3_s=(
            tangent_assembly.ledger.interface_area_m2
            * abs(tangent_candidate.front_speed_m_s)
        ),
        tangent_interface_temperature_k=tangent_candidate.interface_temperature_k,
        newborn_temperature_gradient_scale_k_m=max(
            abs(tangent_candidate.dry_temperature_gradient_k_m),
            _GRADIENT_SCALE_FLOOR,
        ),
        newborn_composition_gradient_scale_m_inv=max(
            abs(composition_gradient),
            _GRADIENT_SCALE_FLOOR,
        ),
        stefan_flux_scales_mol_m2_s=stefan_scales,
        q_scale_m3_s=q_scale,
    )


def solve_face_departure(
    before: face_solver.FaceArrivalCommittedState,
    controls: ci.CutSolverControls,
    dt_s: float,
    boundary: ct.PoreBoundary,
    seed: FaceDepartureSeed,
) -> FaceDepartureIntegratorStep:
    """Try the requested root directly, then use uncommitted fallback stages."""

    evaluations = 0
    rejected_trials = 0
    maximum_scaled = math.inf
    condition_proxy = math.inf
    raw_solver_condition_proxy = math.inf
    last_candidate = None
    last_scaled: tuple[float, ...] = ()
    last_blocks = None
    minimum_bound_distance = 0.0
    nonlinear_method = "safeguarded-direct-newton"
    continuation_fallback_used = False
    direct_requested_attempt_evaluations = 0
    direct_requested_failure_message: str | None = None
    direct_requested_failure_maximum_scaled_residual: float | None = None
    cauchy_globalization_accepted_steps = 0
    cauchy_line_search_trials = 0
    cauchy_residual_evaluations = 0
    dogleg_globalization_accepted_steps = 0
    dogleg_jacobian_linearizations = 0
    dogleg_trial_attempts = 0
    dogleg_residual_evaluations = 0
    audit_progress: _DirectStageProgress | None = None
    try:
        if not isinstance(before, face_solver.FaceArrivalCommittedState):
            raise TypeError("before must be FaceArrivalCommittedState")
        if not isinstance(controls, ci.CutSolverControls):
            raise TypeError("controls must be CutSolverControls")
        if controls != before.solver_controls:
            raise ValueError(
                "face departure must preserve the exact arrival solver controls"
            )
        ci.certify_explicit_wet_retained_water_applicability(
            controls,
            before.wet_retained_water_loadings,
            (0.0,) * len(before.wet_retained_water_loadings),
            before.config.wet.luikov,
        )
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("boundary must be a supported pore boundary")
        if not isinstance(seed, FaceDepartureSeed):
            raise TypeError("seed must be FaceDepartureSeed")
        if not math.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("departure dt must be positive and finite")
        state = departure.adapt_face_arrival_committed_state(before)
        layout = departure.layout_for_departure(state)
        structural = (
            _departure_jacobian_sparsity(layout)
            if layout.unknown_count > controls.maximum_dense_unknowns
            else None
        )
        scale_reference = _exact_tangent_scale_reference(before)
        if len(scale_reference.stefan_flux_scales_mol_m2_s) != layout.dry_face_count:
            raise RuntimeError("exact-tangent Stefan scale rank is inconsistent")
        target_physical = _physical_departure_vector(
            seed.candidate,
            _seed_trace_gradients(seed, state),
            state,
            dt_s,
        )
        auxiliary: list[FaceDepartureAuxiliaryStage] = []
        final: _DirectStageSolution | None = None
        direct_seed = FaceDepartureSeed(
            seed.candidate,
            scale_reference.stefan_flux_scales_mol_m2_s,
            f"{seed.label}; direct requested dt={dt_s:.17g} s",
            _seed_trace_gradients(seed, state),
        )
        direct_progress = _DirectStageProgress()
        audit_progress = direct_progress
        try:
            final = _solve_direct_stage(
                state,
                controls,
                dt_s,
                boundary,
                direct_seed,
                scale_reference,
                controls.maximum_function_evaluations,
                direct_progress,
                certify=True,
                structural=structural,
            )
        except (
            _EvaluationBudgetError,
            _DepartureConditioningError,
            _DirectStageBasinError,
            conditioning.JacobianCertificateError,
            departure.FaceDepartureTopologyError,
            departure.FaceDepartureStepError,
        ) as direct_error:
            evaluations += direct_progress.evaluations
            cauchy_globalization_accepted_steps += (
                direct_progress.cauchy_globalization_accepted_steps
            )
            cauchy_line_search_trials += direct_progress.cauchy_line_search_trials
            cauchy_residual_evaluations += (
                direct_progress.cauchy_residual_evaluations
            )
            dogleg_globalization_accepted_steps += (
                direct_progress.dogleg_globalization_accepted_steps
            )
            dogleg_jacobian_linearizations += (
                direct_progress.dogleg_jacobian_linearizations
            )
            dogleg_trial_attempts += direct_progress.dogleg_trial_attempts
            dogleg_residual_evaluations += (
                direct_progress.dogleg_residual_evaluations
            )
            direct_requested_attempt_evaluations = direct_progress.evaluations
            rejected_trials += direct_progress.rejected_trials
            maximum_scaled = direct_progress.maximum_scaled_residual
            condition_proxy = direct_progress.condition_proxy
            raw_solver_condition_proxy = direct_progress.raw_condition_proxy
            last_candidate = direct_progress.candidate
            last_scaled = direct_progress.signed_scaled_residuals
            last_blocks = direct_progress.residual_blocks
            direct_requested_failure_message = str(direct_error)
            direct_requested_failure_maximum_scaled_residual = (
                direct_progress.maximum_scaled_residual
            )
            if isinstance(
                direct_error,
                (
                    _EvaluationBudgetError,
                    _DepartureConditioningError,
                    conditioning.JacobianCertificateError,
                ),
            ):
                raise
            continuation_fallback_used = True
        else:
            evaluations = final.evaluations
            direct_requested_attempt_evaluations = final.evaluations
            cauchy_globalization_accepted_steps = (
                final.cauchy_globalization_accepted_steps
            )
            cauchy_line_search_trials = final.cauchy_line_search_trials
            cauchy_residual_evaluations = final.cauchy_residual_evaluations
            dogleg_globalization_accepted_steps = (
                final.dogleg_globalization_accepted_steps
            )
            dogleg_jacobian_linearizations = (
                final.dogleg_jacobian_linearizations
            )
            dogleg_trial_attempts = final.dogleg_trial_attempts
            dogleg_residual_evaluations = final.dogleg_residual_evaluations

        if continuation_fallback_used:
            physical_seed = target_physical
            durations = _auxiliary_departure_durations(dt_s)
            for stage_index, stage_dt in enumerate(durations):
                decoded_seed = _decode_physical_departure_vector(
                    physical_seed,
                    state,
                    layout,
                    stage_dt,
                )
                stage_seed = FaceDepartureSeed(
                    decoded_seed.candidate,
                    scale_reference.stefan_flux_scales_mol_m2_s,
                    f"{seed.label}; same-old-state stage dt={stage_dt:.17g} s",
                    decoded_seed.newborn_trace_gradients,
                )
                progress = _DirectStageProgress()
                audit_progress = progress
                is_requested_stage = stage_index + 1 == len(durations)
                try:
                    stage = _solve_direct_stage(
                        state,
                        controls,
                        stage_dt,
                        boundary,
                        stage_seed,
                        scale_reference,
                        controls.maximum_function_evaluations - evaluations,
                        progress,
                        certify=is_requested_stage,
                        structural=structural,
                    )
                except Exception:
                    evaluations += progress.evaluations
                    cauchy_globalization_accepted_steps += (
                        progress.cauchy_globalization_accepted_steps
                    )
                    cauchy_line_search_trials += (
                        progress.cauchy_line_search_trials
                    )
                    cauchy_residual_evaluations += (
                        progress.cauchy_residual_evaluations
                    )
                    dogleg_globalization_accepted_steps += (
                        progress.dogleg_globalization_accepted_steps
                    )
                    dogleg_jacobian_linearizations += (
                        progress.dogleg_jacobian_linearizations
                    )
                    dogleg_trial_attempts += progress.dogleg_trial_attempts
                    dogleg_residual_evaluations += (
                        progress.dogleg_residual_evaluations
                    )
                    rejected_trials += progress.rejected_trials
                    maximum_scaled = progress.maximum_scaled_residual
                    condition_proxy = progress.condition_proxy
                    raw_solver_condition_proxy = progress.raw_condition_proxy
                    last_candidate = progress.candidate
                    last_scaled = progress.signed_scaled_residuals
                    last_blocks = progress.residual_blocks
                    raise
                evaluations += stage.evaluations
                cauchy_globalization_accepted_steps += (
                    stage.cauchy_globalization_accepted_steps
                )
                cauchy_line_search_trials += stage.cauchy_line_search_trials
                cauchy_residual_evaluations += stage.cauchy_residual_evaluations
                dogleg_globalization_accepted_steps += (
                    stage.dogleg_globalization_accepted_steps
                )
                dogleg_jacobian_linearizations += (
                    stage.dogleg_jacobian_linearizations
                )
                dogleg_trial_attempts += stage.dogleg_trial_attempts
                dogleg_residual_evaluations += stage.dogleg_residual_evaluations
                rejected_trials += stage.rejected_trials
                maximum_scaled = max(
                    abs(value) for value in stage.signed_scaled_residuals
                )
                condition_proxy = (
                    math.nan
                    if stage.condition_certificate is None
                    else stage.condition_certificate.condition_number_2
                )
                raw_solver_condition_proxy = stage.raw_condition_proxy
                last_candidate = stage.assembly.candidate
                last_scaled = stage.signed_scaled_residuals
                last_blocks = stage.assembly.residuals
                physical_seed = _physical_departure_vector(
                    stage.assembly.candidate,
                    stage.newborn_trace_gradients,
                    state,
                    stage_dt,
                )
                if is_requested_stage:
                    final = replace(
                        stage,
                        message=(
                            "direct physical Newton root after same-old-state "
                            "continuation fallback; unchanged ordinary residual"
                        ),
                    )
                else:
                    auxiliary.append(
                        FaceDepartureAuxiliaryStage(
                            dt_s=stage_dt,
                            nonlinear_evaluations=stage.evaluations,
                            rejected_trial_evaluations=stage.rejected_trials,
                            maximum_scaled_residual=maximum_scaled,
                            condition_proxy=None,
                            raw_solver_condition_proxy=None,
                            cauchy_globalization_accepted_steps=(
                                stage.cauchy_globalization_accepted_steps
                            ),
                            cauchy_line_search_trials=(
                                stage.cauchy_line_search_trials
                            ),
                            cauchy_residual_evaluations=(
                                stage.cauchy_residual_evaluations
                            ),
                            dogleg_globalization_accepted_steps=(
                                stage.dogleg_globalization_accepted_steps
                            ),
                            dogleg_jacobian_linearizations=(
                                stage.dogleg_jacobian_linearizations
                            ),
                            dogleg_trial_attempts=stage.dogleg_trial_attempts,
                            dogleg_residual_evaluations=(
                                stage.dogleg_residual_evaluations
                            ),
                        )
                    )
        if final is None or final.assembly.dt_s != dt_s:
            raise RuntimeError("departure continuation lost its requested-dt root")
        if final.condition_certificate is None:
            raise RuntimeError("requested departure root lacks its fixed certificate")
        maximum_scaled = max(abs(value) for value in final.signed_scaled_residuals)
        condition_proxy = final.condition_certificate.condition_number_2
        raw_solver_condition_proxy = final.raw_condition_proxy
        last_candidate = final.assembly.candidate
        last_scaled = final.signed_scaled_residuals
        last_blocks = final.assembly.residuals
        fractions = _bounded_fractions(
            final.assembly.candidate,
            state,
            final.chart,
        )
        minimum_bound_distance = min(
            min(value, 1.0 - value) for value in fractions
        )
        if dogleg_globalization_accepted_steps:
            nonlinear_method = (
                "safeguarded-direct-newton-cauchy-with-diagnostic-"
                "powell-dogleg-rescue"
            )
        elif cauchy_globalization_accepted_steps:
            nonlinear_method = (
                "safeguarded-direct-newton-with-scaled-cauchy-globalization"
            )
        step = _accept_departure(
            before,
            controls,
            boundary,
            seed,
            final.assembly,
            final.scales,
            final.signed_scaled_residuals,
            evaluations,
            rejected_trials,
            final.message,
            final.condition_certificate.condition_number_2,
            final.raw_condition_proxy,
            final.physical_variable_scale_labels,
            final.physical_variable_scales,
            final.physical_jacobian_perturbations,
            final.condition_certificate,
            final.chart,
            final.structural,
            nonlinear_method,
            tuple(auxiliary),
            final.newborn_trace_gradients,
            continuation_fallback_used,
            direct_requested_attempt_evaluations,
            direct_requested_failure_message,
            direct_requested_failure_maximum_scaled_residual,
            cauchy_globalization_accepted_steps,
            cauchy_line_search_trials,
            cauchy_residual_evaluations,
            dogleg_globalization_accepted_steps,
            dogleg_jacobian_linearizations,
            dogleg_trial_attempts,
            dogleg_residual_evaluations,
        )
        if (
            step.ledger.maximum_step_ledger_residual
            > controls.ledger_tolerance
            or step.ledger.maximum_cumulative_ledger_residual
            > controls.ledger_tolerance
        ):
            raise RuntimeError(
                "face-departure accepted/cumulative ledger exceeded 1e-10"
            )
        ci.certify_explicit_wet_retained_water_applicability(
            controls,
            before.wet_retained_water_loadings,
            (0.0,) * len(before.wet_retained_water_loadings),
            before.config.wet.luikov,
        )
        ci.certify_explicit_wet_retained_water_applicability(
            controls,
            step.after.transport.wet_retained_water_loadings,
            step.after.transport.effective_wet_retained_water_capacity_duals_over_rt,
            step.after.transport.config.wet.luikov,
        )
        return step
    except FaceDepartureIntegratorStepError:
        raise
    except Exception as exc:
        if not isinstance(before, face_solver.FaceArrivalCommittedState):
            raise
        applicability_error = (
            exc
            if isinstance(exc, ci.WetRetainedWaterApplicabilityError)
            else None
        )
        raise FaceDepartureIntegratorStepError(
            f"face-departure solve rejected with exact rollback: {exc}",
            before,
            nonlinear_evaluations=evaluations,
            rejected_trial_evaluations=rejected_trials,
            maximum_scaled_residual=maximum_scaled,
            condition_proxy=condition_proxy,
            raw_solver_condition_proxy=raw_solver_condition_proxy,
            last_candidate=last_candidate,
            scaled_residuals=last_scaled,
            residual_blocks=last_blocks,
            minimum_fractional_distance_to_bound=minimum_bound_distance,
            nonlinear_method=nonlinear_method,
            retained_water_applicability_guard_failed=(
                applicability_error is not None
            ),
            retained_water_applicability_audit=(
                None if applicability_error is None else applicability_error.audit
            ),
            direct_newton_iteration_audits=(
                tuple(audit_progress.direct_newton_iteration_audits)
                if audit_progress is not None
                else ()
            ),
            cauchy_globalization_accepted_steps=(
                cauchy_globalization_accepted_steps
            ),
            cauchy_line_search_trials=cauchy_line_search_trials,
            cauchy_residual_evaluations=cauchy_residual_evaluations,
            dogleg_globalization_accepted_steps=(
                dogleg_globalization_accepted_steps
            ),
            dogleg_jacobian_linearizations=dogleg_jacobian_linearizations,
            dogleg_trial_attempts=dogleg_trial_attempts,
            dogleg_residual_evaluations=dogleg_residual_evaluations,
            direct_dogleg_iteration_audits=(
                tuple(audit_progress.direct_dogleg_iteration_audits)
                if audit_progress is not None
                else ()
            ),
        ) from exc


def audit_face_departure_sparsity(
    step: FaceDepartureIntegratorStep,
) -> FaceDepartureSparsityAudit:
    """Compare colored groups with independent one-column central probes.

    This is an explicit qualification audit, not part of the committed solve;
    its residual evaluations are reported separately and never hidden in the
    step ledger.  Both routes use the accepted ordinary ALE/RH residual, the
    accepted fixed row scales, and the same direct physical chart.
    """

    if not isinstance(step, FaceDepartureIntegratorStep):
        raise TypeError("departure sparsity audit requires an accepted step")
    before = departure.adapt_face_arrival_committed_state(step.before)
    layout = departure.layout_for_departure(before)
    structural = _departure_jacobian_sparsity(layout)
    reference = _exact_tangent_scale_reference(step.before)
    chart = _coordinate_chart(
        before,
        step.controls,
        step.assembly.candidate,
        step.ledger.newborn_trace_gradients,
        step.assembly.dt_s,
        step.boundary,
        reference.q_scale_m3_s,
    )
    coordinates = _encode_candidate(
        step.assembly.candidate,
        before,
        chart,
        reference.stefan_flux_scales_mol_m2_s,
        step.ledger.newborn_trace_gradients,
    )
    lower, upper = _coordinate_solver_bounds(
        layout,
        chart,
        step.controls.maximum_logit_magnitude,
    )
    coordinate_scales = np.asarray(
        _solver_coordinate_characteristic_scales(before, chart, reference),
        dtype=float,
    )
    residual_scale_vector = np.asarray(
        step.ledger.residual_scales.vector,
        dtype=float,
    )
    evaluations = 0

    def residual(
        trial_coordinates: np.ndarray,
    ) -> tuple[np.ndarray, _DecodedTrial, departure.FaceDepartureAssembly]:
        nonlocal evaluations
        decoded = _decode_trial(
            trial_coordinates,
            before,
            layout,
            chart,
            reference.stefan_flux_scales_mol_m2_s,
        )
        assembly = departure.assemble_face_departure(
            before,
            decoded.candidate,
            step.assembly.dt_s,
            step.boundary,
            newborn_trace_gradients=decoded.newborn_trace_gradients,
            enforce_reduced_film_thresholds=False,
        )
        evaluations += 1
        values = np.asarray(
            departure.datum_covariant_residual_vector(assembly.residuals),
            dtype=float,
        ) / residual_scale_vector
        return values, decoded, assembly

    independent, independent_perturbations = _direct_coordinate_jacobian(
        coordinates,
        lower,
        upper,
        coordinate_scales,
        residual,
        layout,
        None,
    )
    colored, _ = _direct_coordinate_jacobian(
        coordinates,
        lower,
        upper,
        coordinate_scales,
        residual,
        layout,
        structural,
    )
    outside = _certify_departure_sparsity_pattern(independent, structural)
    colored_dense = np.asarray(colored.toarray(), dtype=float)
    row_scales = np.maximum(np.max(np.abs(independent), axis=1), 1.0e-300)
    mask = structural.jac_sparsity.toarray()
    relative_difference = np.where(
        mask,
        np.abs(independent - colored_dense) / row_scales[:, None],
        0.0,
    )
    maximum_colored_difference = float(np.max(relative_difference))
    if maximum_colored_difference > _SPARSITY_RELATIVE_ZERO_TOLERANCE:
        raise RuntimeError(
            "colored departure Jacobian disagrees with independent columns: "
            f"{maximum_colored_difference:.3e}"
        )
    independent_certificate, _, _, _ = _certify_affine_direct_jacobian(
        before,
        independent,
        independent_perturbations,
        step.ledger.residual_scales,
        chart,
        reference,
    )
    independent_condition = independent_certificate.condition_number_2
    accepted_condition = step.ledger.condition_proxy
    relative_condition_difference = abs(
        independent_condition - accepted_condition
    ) / max(independent_condition, accepted_condition)
    return FaceDepartureSparsityAudit(
        unknown_count=layout.unknown_count,
        structural_nnz=structural.nnz,
        color_count=structural.color_count,
        residual_evaluations=evaluations,
        maximum_relative_derivative_outside_pattern=outside,
        maximum_relative_colored_difference_inside_pattern=(
            maximum_colored_difference
        ),
        independent_condition_proxy=independent_condition,
        accepted_condition_proxy=accepted_condition,
        relative_condition_proxy_difference=relative_condition_difference,
        passed=True,
    )


def _auxiliary_departure_durations(
    requested_dt_s: float,
    *,
    representability_floor_s: float = _FALLBACK_REPRESENTABILITY_FLOOR_S,
) -> tuple[float, ...]:
    """Geometric basin continuation; the floor has no physical meaning."""

    if not math.isfinite(requested_dt_s) or requested_dt_s <= 0.0:
        raise ValueError("requested departure duration must be positive and finite")
    if not math.isfinite(representability_floor_s) or representability_floor_s <= 0.0:
        raise ValueError("fallback representability floor must be positive and finite")
    if requested_dt_s <= representability_floor_s:
        return (requested_dt_s,)
    start = min(
        requested_dt_s,
        max(representability_floor_s, 1.0e-3 * requested_dt_s),
    )
    durations = [start]
    while 10.0 * durations[-1] < requested_dt_s:
        durations.append(10.0 * durations[-1])
    if durations[-1] != requested_dt_s:
        durations.append(requested_dt_s)
    result = tuple(durations)
    if not all(
        left < right for left, right in zip(result, result[1:])
    ) or result[-1] != requested_dt_s:
        raise ValueError("same-old-state departure durations are unresolved")
    return result


def _departure_jacobian_sparsity(
    layout: departure.FaceDepartureLayout,
) -> cut_sparsity.CutJacobianSparsity:
    """Map the affine departure chart to the exact cut residual structure.

    The direct departure blocks have exactly the cut ordering.  Within the
    dry blocks, ``g_T`` and ``g_y`` replace only the newborn cell primitives;
    ``Q`` is the affine replacement for the dense ``front_z`` border.  The
    direct outer composition rates do not carry the triangular temperature
    dependence of the ordinary face-conditioned logit chart, so that optional
    expansion is deliberately disabled.  ``T_Gamma`` and ``Q`` remain dense
    border columns.
    """

    if not isinstance(layout, departure.FaceDepartureLayout):
        raise TypeError("departure sparsity requires FaceDepartureLayout")
    transport_layout = cut.CutTransportLayout(
        wet_cell_indices=layout.wet_cell_indices,
        dry_cell_indices=layout.dry_cell_indices,
        cut_cell_index=layout.cut_cell_index,
    )
    if (
        transport_layout.unknown_blocks != layout.unknown_blocks
        or transport_layout.residual_blocks != layout.residual_blocks
    ):
        raise RuntimeError("departure/cut structural block order is inconsistent")
    return cut_sparsity.build_cut_jacobian_sparsity(
        transport_layout,
        face_conditioned_composition=False,
    )


def _seed_trace_gradients(
    seed: FaceDepartureSeed,
    before: departure.FaceDepartureState,
) -> departure.NewbornDryTraceGradients:
    if seed.newborn_trace_gradients is not None:
        return seed.newborn_trace_gradients
    geometry = cg.partition_master_grid(
        before.geometry.master_grid,
        z=seed.candidate.front_z,
    )
    layout = cut.layout_for_geometry(geometry)
    distance = _newborn_center_distance(geometry, layout.cut_cell_index)
    interface = cut.evaluate_interface_state(
        seed.candidate.interface_temperature_k,
        before.config,
        before.historical_hexane_loadings[layout.cut_cell_index],
        before.oil_fraction_labels[layout.cut_cell_index],
    )
    return departure.NewbornDryTraceGradients(
        (
            seed.candidate.dry_temperatures_k[0]
            - seed.candidate.interface_temperature_k
        )
        / distance,
        (seed.candidate.dry_y_hexane[0] - interface.y_hexane) / distance,
    )


def _tangent_reference_physical_vector(
    layout: departure.FaceDepartureLayout,
    reference: _PhysicalScaleReference,
) -> tuple[float, ...]:
    values = (
        *((0.0,) * layout.wet_piece_count),
        *((0.0,) * layout.wet_piece_count),
        reference.tangent_temperature_gradient_k_m,
        *((0.0,) * (layout.dry_piece_count - 1)),
        reference.tangent_composition_gradient_m_inv,
        *((0.0,) * (layout.dry_piece_count - 1)),
        *reference.tangent_stefan_fluxes_mol_m2_s,
        reference.tangent_q_m3_s,
        reference.tangent_interface_temperature_k,
    )
    if len(values) != layout.unknown_count:
        raise RuntimeError("exact-tangent physical reference lost block rank")
    return values


def _solve_direct_stage(
    before: departure.FaceDepartureState,
    controls: ci.CutSolverControls,
    dt_s: float,
    boundary: ct.PoreBoundary,
    seed: FaceDepartureSeed,
    scale_reference: _PhysicalScaleReference,
    evaluation_budget: int,
    progress: _DirectStageProgress,
    *,
    certify: bool,
    structural: cut_sparsity.CutJacobianSparsity | None,
    diagnostic_scaled_cauchy_globalization: bool = False,
    diagnostic_powell_dogleg_rescue: bool = False,
) -> _DirectStageSolution:
    if not isinstance(diagnostic_scaled_cauchy_globalization, bool):
        raise TypeError("diagnostic Cauchy switch must be boolean")
    if not isinstance(diagnostic_powell_dogleg_rescue, bool):
        raise TypeError("diagnostic dogleg switch must be boolean")
    if evaluation_budget < 1:
        raise _EvaluationBudgetError(
            "departure solve exhausted its cumulative function-evaluation budget"
        )
    layout = departure.layout_for_departure(before)
    reference_decoded = _decode_physical_departure_vector(
        _tangent_reference_physical_vector(layout, scale_reference),
        before,
        layout,
        dt_s,
    )
    chart = _coordinate_chart(
        before,
        controls,
        reference_decoded.candidate,
        reference_decoded.newborn_trace_gradients,
        dt_s,
        boundary,
        scale_reference.q_scale_m3_s,
    )
    _validate_chart(before, controls, chart)
    _validate_seed(seed, before, layout, chart)
    coordinates = _encode_candidate(
        seed.candidate,
        before,
        chart,
        seed.stefan_flux_scales_mol_m2_s,
        seed.newborn_trace_gradients,
    )
    lower, upper = _coordinate_solver_bounds(
        layout,
        chart,
        controls.maximum_logit_magnitude,
    )
    coordinate_scales = np.asarray(
        _solver_coordinate_characteristic_scales(
            before,
            chart,
            scale_reference,
        ),
        dtype=float,
    )
    initial = _decode_trial(
        coordinates,
        before,
        layout,
        chart,
        seed.stefan_flux_scales_mol_m2_s,
    )
    assembly = departure.assemble_face_departure(
        before,
        initial.candidate,
        dt_s,
        boundary,
        newborn_trace_gradients=initial.newborn_trace_gradients,
        enforce_reduced_film_thresholds=False,
    )
    evaluations = 1
    rejected_trials = 0
    progress.evaluations = evaluations
    progress.candidate = assembly.candidate
    progress.residual_blocks = assembly.residuals
    progress.current_candidate = assembly.candidate
    progress.current_newborn_trace_gradients = initial.newborn_trace_gradients
    progress.current_solver_coordinates = tuple(float(value) for value in coordinates)
    if evaluations >= evaluation_budget:
        raise _EvaluationBudgetError(
            "departure solve exhausted its cumulative function-evaluation budget"
        )
    reference_assembly = departure.assemble_face_departure(
        before,
        reference_decoded.candidate,
        dt_s,
        boundary,
        newborn_trace_gradients=reference_decoded.newborn_trace_gradients,
        enforce_reduced_film_thresholds=False,
    )
    evaluations += 1
    scales = _residual_scales(before, reference_assembly, dt_s)
    scale_vector = np.asarray(scales.vector, dtype=float)
    progress.evaluations = evaluations
    progress.candidate = assembly.candidate
    progress.residual_blocks = assembly.residuals

    def residual(
        trial_coordinates: np.ndarray,
    ) -> tuple[np.ndarray, _DecodedTrial, departure.FaceDepartureAssembly]:
        nonlocal evaluations, rejected_trials
        if evaluations >= evaluation_budget:
            raise _EvaluationBudgetError(
                "departure solve exhausted its cumulative function-evaluation budget"
            )
        evaluations += 1
        progress.evaluations = evaluations
        try:
            decoded = _decode_trial(
                trial_coordinates,
                before,
                layout,
                chart,
                seed.stefan_flux_scales_mol_m2_s,
            )
            trial_assembly = departure.assemble_face_departure(
                before,
                decoded.candidate,
                dt_s,
                boundary,
                newborn_trace_gradients=decoded.newborn_trace_gradients,
                enforce_reduced_film_thresholds=False,
            )
        except (ValueError, departure.FaceDepartureStepError):
            rejected_trials += 1
            progress.rejected_trials = rejected_trials
            raise
        values = np.asarray(
            departure.datum_covariant_residual_vector(trial_assembly.residuals),
            dtype=float,
        ) / scale_vector
        progress.candidate = decoded.candidate
        progress.residual_blocks = trial_assembly.residuals
        progress.signed_scaled_residuals = tuple(float(value) for value in values)
        progress.maximum_scaled_residual = float(np.max(np.abs(values)))
        return values, decoded, trial_assembly

    current = np.asarray(
        departure.datum_covariant_residual_vector(assembly.residuals),
        dtype=float,
    ) / scale_vector
    progress.signed_scaled_residuals = tuple(float(value) for value in current)
    progress.maximum_scaled_residual = float(np.max(np.abs(current)))
    decoded = initial
    raw_condition = math.inf
    for _iteration in range(20):
        if float(np.max(np.abs(current))) <= controls.nonlinear_residual_tolerance:
            break
        iteration_current_norm = float(np.linalg.norm(current))
        iteration_current_maximum = float(np.max(np.abs(current)))
        jacobian, jacobian_perturbations = _direct_coordinate_jacobian(
            coordinates,
            lower,
            upper,
            coordinate_scales,
            residual,
            layout,
            structural,
        )
        scaled_jacobian = (
            jacobian * coordinate_scales[None, :]
            if structural is None
            else jacobian.multiply(coordinate_scales[None, :])
        )
        scaled_correction = _scaled_newton_correction(
            scaled_jacobian,
            current,
        )
        correction = coordinate_scales * scaled_correction
        jacobian_action = scaled_jacobian @ scaled_correction
        predicted_directional_derivative = float(
            np.dot(current, jacobian_action)
        )
        ideal_directional_derivative = -float(np.dot(current, current))
        bound_fractions = np.full(coordinates.shape, math.inf, dtype=float)
        positive = correction > 0.0
        negative = correction < 0.0
        bound_fractions[positive] = (
            (upper[positive] - coordinates[positive]) / correction[positive]
        )
        bound_fractions[negative] = (
            (coordinates[negative] - lower[negative]) / -correction[negative]
        )
        maximum_open_bound_fraction = float(np.min(bound_fractions))
        accepted = False
        accepted_fraction: float | None = None
        accepted_route: str | None = None
        line_search_trials: list[_DirectNewtonTrialAudit] = []
        for backtrack in range(24):
            fraction = 0.5**backtrack
            trial_coordinates = coordinates + fraction * correction
            if np.any(trial_coordinates <= lower) or np.any(trial_coordinates >= upper):
                rejected_trials += 1
                progress.rejected_trials = rejected_trials
                line_search_trials.append(
                    _DirectNewtonTrialAudit(
                        fraction=fraction,
                        rejection_reason="outside_open_coordinate_bounds",
                        residual_l2_norm=None,
                        maximum_scaled_residual=None,
                        strictly_decreased_l2_merit=False,
                    )
                )
                continue
            try:
                trial, trial_decoded, trial_assembly = residual(trial_coordinates)
            except (ValueError, departure.FaceDepartureStepError) as exc:
                line_search_trials.append(
                    _DirectNewtonTrialAudit(
                        fraction=fraction,
                        rejection_reason=f"{type(exc).__name__}: {exc}",
                        residual_l2_norm=None,
                        maximum_scaled_residual=None,
                        strictly_decreased_l2_merit=False,
                    )
                )
                continue
            trial_norm = float(np.linalg.norm(trial))
            decreased = trial_norm < iteration_current_norm
            line_search_trials.append(
                _DirectNewtonTrialAudit(
                    fraction=fraction,
                    rejection_reason=(None if decreased else "l2_merit_not_decreased"),
                    residual_l2_norm=trial_norm,
                    maximum_scaled_residual=float(np.max(np.abs(trial))),
                    strictly_decreased_l2_merit=decreased,
                )
            )
            if decreased:
                coordinates = trial_coordinates
                current = trial
                decoded = trial_decoded
                assembly = trial_assembly
                progress.current_candidate = trial_assembly.candidate
                progress.current_newborn_trace_gradients = (
                    trial_decoded.newborn_trace_gradients
                )
                progress.current_solver_coordinates = tuple(
                    float(value) for value in trial_coordinates
                )
                accepted = True
                accepted_fraction = fraction
                accepted_route = "newton"
                break
        cauchy_gradient_l2_norm: float | None = None
        cauchy_numerator: float | None = None
        cauchy_denominator: float | None = None
        cauchy_alpha: float | None = None
        cauchy_scaled_correction_l2_norm: float | None = None
        cauchy_maximum_open_bound_fraction: float | None = None
        cauchy_line_search_trials: list[_DirectNewtonTrialAudit] = []
        cauchy_accepted_fraction: float | None = None
        cauchy_enabled = (
            diagnostic_scaled_cauchy_globalization
            or diagnostic_powell_dogleg_rescue
        )
        if not accepted and cauchy_enabled:
            (
                scaled_cauchy_correction,
                cauchy_gradient_l2_norm,
                cauchy_numerator,
                cauchy_denominator,
                cauchy_alpha,
            ) = _scaled_cauchy_correction(scaled_jacobian, current)
            cauchy_scaled_correction_l2_norm = float(
                np.linalg.norm(scaled_cauchy_correction)
            )
            cauchy_correction = coordinate_scales * scaled_cauchy_correction
            cauchy_bound_fractions = np.full(
                coordinates.shape,
                math.inf,
                dtype=float,
            )
            cauchy_positive = cauchy_correction > 0.0
            cauchy_negative = cauchy_correction < 0.0
            cauchy_bound_fractions[cauchy_positive] = (
                (upper[cauchy_positive] - coordinates[cauchy_positive])
                / cauchy_correction[cauchy_positive]
            )
            cauchy_bound_fractions[cauchy_negative] = (
                (coordinates[cauchy_negative] - lower[cauchy_negative])
                / -cauchy_correction[cauchy_negative]
            )
            cauchy_maximum_open_bound_fraction = float(
                np.min(cauchy_bound_fractions)
            )
            for backtrack in range(24):
                fraction = 0.5**backtrack
                progress.cauchy_line_search_trials += 1
                trial_coordinates = coordinates + fraction * cauchy_correction
                if np.any(trial_coordinates <= lower) or np.any(
                    trial_coordinates >= upper
                ):
                    rejected_trials += 1
                    progress.rejected_trials = rejected_trials
                    cauchy_line_search_trials.append(
                        _DirectNewtonTrialAudit(
                            fraction=fraction,
                            rejection_reason="outside_open_coordinate_bounds",
                            residual_l2_norm=None,
                            maximum_scaled_residual=None,
                            strictly_decreased_l2_merit=False,
                        )
                    )
                    continue
                evaluations_before_trial = evaluations
                try:
                    trial, trial_decoded, trial_assembly = residual(
                        trial_coordinates
                    )
                except (ValueError, departure.FaceDepartureStepError) as exc:
                    cauchy_line_search_trials.append(
                        _DirectNewtonTrialAudit(
                            fraction=fraction,
                            rejection_reason=f"{type(exc).__name__}: {exc}",
                            residual_l2_norm=None,
                            maximum_scaled_residual=None,
                            strictly_decreased_l2_merit=False,
                        )
                    )
                    continue
                finally:
                    progress.cauchy_residual_evaluations += (
                        evaluations - evaluations_before_trial
                    )
                trial_norm = float(np.linalg.norm(trial))
                decreased = trial_norm < iteration_current_norm
                cauchy_line_search_trials.append(
                    _DirectNewtonTrialAudit(
                        fraction=fraction,
                        rejection_reason=(
                            None if decreased else "l2_merit_not_decreased"
                        ),
                        residual_l2_norm=trial_norm,
                        maximum_scaled_residual=float(np.max(np.abs(trial))),
                        strictly_decreased_l2_merit=decreased,
                    )
                )
                if decreased:
                    coordinates = trial_coordinates
                    current = trial
                    decoded = trial_decoded
                    assembly = trial_assembly
                    progress.current_candidate = trial_assembly.candidate
                    progress.current_newborn_trace_gradients = (
                        trial_decoded.newborn_trace_gradients
                    )
                    progress.current_solver_coordinates = tuple(
                        float(value) for value in trial_coordinates
                    )
                    progress.cauchy_globalization_accepted_steps += 1
                    accepted = True
                    accepted_route = "scaled_cauchy"
                    cauchy_accepted_fraction = fraction
                    break
        relative_perturbations = jacobian_perturbations / coordinate_scales
        progress.direct_newton_iteration_audits.append(
            _DirectNewtonIterationAudit(
                iteration=_iteration,
                current_l2_norm=iteration_current_norm,
                current_maximum_scaled_residual=iteration_current_maximum,
                jacobian_route=(
                    "independent_dense_columns"
                    if structural is None
                    else "structurally_colored_sparse_columns"
                ),
                minimum_coordinate_perturbation=float(
                    np.min(jacobian_perturbations)
                ),
                maximum_coordinate_perturbation=float(
                    np.max(jacobian_perturbations)
                ),
                minimum_relative_coordinate_perturbation=float(
                    np.min(relative_perturbations)
                ),
                maximum_relative_coordinate_perturbation=float(
                    np.max(relative_perturbations)
                ),
                predicted_merit_directional_derivative=(
                    predicted_directional_derivative
                ),
                ideal_newton_directional_derivative=(
                    ideal_directional_derivative
                ),
                scaled_correction_l2_norm=float(
                    np.linalg.norm(scaled_correction)
                ),
                maximum_open_bound_fraction=maximum_open_bound_fraction,
                line_search_trials=tuple(line_search_trials),
                accepted_fraction=accepted_fraction,
                accepted_route=accepted_route,
                cauchy_gradient_l2_norm=cauchy_gradient_l2_norm,
                cauchy_numerator=cauchy_numerator,
                cauchy_denominator=cauchy_denominator,
                cauchy_alpha=cauchy_alpha,
                cauchy_scaled_correction_l2_norm=(
                    cauchy_scaled_correction_l2_norm
                ),
                cauchy_maximum_open_bound_fraction=(
                    cauchy_maximum_open_bound_fraction
                ),
                cauchy_line_search_trials=tuple(cauchy_line_search_trials),
                cauchy_accepted_fraction=cauchy_accepted_fraction,
            )
        )
        if not accepted:
            if cauchy_enabled:
                message = (
                    "departure direct Newton and diagnostic scaled Cauchy "
                    "corrections have no residual-decreasing phase-safe step"
                )
            else:
                message = (
                    "departure direct Newton correction has no "
                    "residual-decreasing phase-safe step"
                )
            raise _DirectStageBasinError(message)
    dogleg_trust_radius: float | None = None
    if (
        diagnostic_powell_dogleg_rescue
        and float(np.max(np.abs(current)))
        > controls.nonlinear_residual_tolerance
    ):
        for dogleg_iteration in range(20):
            if (
                float(np.max(np.abs(current)))
                <= controls.nonlinear_residual_tolerance
            ):
                break
            iteration_current_norm = float(np.linalg.norm(current))
            iteration_current_maximum = float(np.max(np.abs(current)))
            jacobian, _ = _direct_coordinate_jacobian(
                coordinates,
                lower,
                upper,
                coordinate_scales,
                residual,
                layout,
                structural,
            )
            progress.dogleg_jacobian_linearizations += 1
            scaled_jacobian = (
                jacobian * coordinate_scales[None, :]
                if structural is None
                else jacobian.multiply(coordinate_scales[None, :])
            )
            scaled_newton = _scaled_newton_correction(
                scaled_jacobian,
                current,
            )
            (
                scaled_cauchy,
                cauchy_gradient_norm,
                cauchy_numerator,
                cauchy_denominator,
                cauchy_alpha,
            ) = _scaled_cauchy_correction(scaled_jacobian, current)
            cauchy_norm = float(np.linalg.norm(scaled_cauchy))
            if dogleg_trust_radius is None:
                dogleg_trust_radius = cauchy_norm
            if (
                not math.isfinite(dogleg_trust_radius)
                or dogleg_trust_radius <= 0.0
            ):
                raise _DirectStageBasinError(
                    "departure diagnostic dogleg trust radius is invalid"
                )
            trust_radius_on_entry = dogleg_trust_radius
            dogleg_trials: list[_DirectDoglegTrialAudit] = []
            dogleg_accepted = False
            for dogleg_trial in range(24):
                progress.dogleg_trial_attempts += 1
                radius_before = dogleg_trust_radius
                scaled_dogleg, dogleg_branch = (
                    _scaled_powell_dogleg_correction(
                        scaled_newton,
                        scaled_cauchy,
                        radius_before,
                    )
                )
                step_norm = float(np.linalg.norm(scaled_dogleg))
                linearized_residual = np.asarray(
                    current + scaled_jacobian @ scaled_dogleg,
                    dtype=float,
                ).reshape(-1)
                predicted_reduction = 0.5 * (
                    iteration_current_norm * iteration_current_norm
                    - float(np.dot(linearized_residual, linearized_residual))
                )
                if (
                    not math.isfinite(predicted_reduction)
                    or predicted_reduction <= 0.0
                ):
                    raise _DirectStageBasinError(
                        "departure diagnostic dogleg has no positive finite "
                        "predicted merit reduction"
                    )
                trial_coordinates = (
                    coordinates + coordinate_scales * scaled_dogleg
                )
                if np.any(trial_coordinates <= lower) or np.any(
                    trial_coordinates >= upper
                ):
                    rejected_trials += 1
                    progress.rejected_trials = rejected_trials
                    dogleg_trust_radius = 0.25 * radius_before
                    dogleg_trials.append(
                        _DirectDoglegTrialAudit(
                            trial=dogleg_trial,
                            trust_radius_before=radius_before,
                            dogleg_branch=dogleg_branch,
                            scaled_step_l2_norm=step_norm,
                            predicted_merit_reduction=predicted_reduction,
                            rejection_reason=(
                                "outside_open_coordinate_bounds"
                            ),
                            residual_l2_norm=None,
                            maximum_scaled_residual=None,
                            actual_merit_reduction=None,
                            trust_ratio=None,
                            accepted=False,
                            trust_radius_after=dogleg_trust_radius,
                        )
                    )
                    continue
                evaluations_before_trial = evaluations
                try:
                    trial, trial_decoded, trial_assembly = residual(
                        trial_coordinates
                    )
                except (ValueError, departure.FaceDepartureStepError) as exc:
                    dogleg_trust_radius = 0.25 * radius_before
                    dogleg_trials.append(
                        _DirectDoglegTrialAudit(
                            trial=dogleg_trial,
                            trust_radius_before=radius_before,
                            dogleg_branch=dogleg_branch,
                            scaled_step_l2_norm=step_norm,
                            predicted_merit_reduction=predicted_reduction,
                            rejection_reason=f"{type(exc).__name__}: {exc}",
                            residual_l2_norm=None,
                            maximum_scaled_residual=None,
                            actual_merit_reduction=None,
                            trust_ratio=None,
                            accepted=False,
                            trust_radius_after=dogleg_trust_radius,
                        )
                    )
                    continue
                finally:
                    progress.dogleg_residual_evaluations += (
                        evaluations - evaluations_before_trial
                    )
                trial_norm = float(np.linalg.norm(trial))
                trial_maximum = float(np.max(np.abs(trial)))
                actual_reduction = 0.5 * (
                    iteration_current_norm * iteration_current_norm
                    - trial_norm * trial_norm
                )
                trust_ratio = actual_reduction / predicted_reduction
                if (
                    math.isfinite(actual_reduction)
                    and math.isfinite(trust_ratio)
                    and actual_reduction > 0.0
                    and trust_ratio >= 0.1
                ):
                    dogleg_trust_radius = _updated_dogleg_trust_radius(
                        radius_before,
                        step_norm,
                        trust_ratio,
                    )
                    dogleg_trials.append(
                        _DirectDoglegTrialAudit(
                            trial=dogleg_trial,
                            trust_radius_before=radius_before,
                            dogleg_branch=dogleg_branch,
                            scaled_step_l2_norm=step_norm,
                            predicted_merit_reduction=predicted_reduction,
                            rejection_reason=None,
                            residual_l2_norm=trial_norm,
                            maximum_scaled_residual=trial_maximum,
                            actual_merit_reduction=actual_reduction,
                            trust_ratio=trust_ratio,
                            accepted=True,
                            trust_radius_after=dogleg_trust_radius,
                        )
                    )
                    coordinates = trial_coordinates
                    current = trial
                    decoded = trial_decoded
                    assembly = trial_assembly
                    progress.current_candidate = trial_assembly.candidate
                    progress.current_newborn_trace_gradients = (
                        trial_decoded.newborn_trace_gradients
                    )
                    progress.current_solver_coordinates = tuple(
                        float(value) for value in trial_coordinates
                    )
                    progress.dogleg_globalization_accepted_steps += 1
                    dogleg_accepted = True
                    break
                dogleg_trust_radius = 0.25 * radius_before
                dogleg_trials.append(
                    _DirectDoglegTrialAudit(
                        trial=dogleg_trial,
                        trust_radius_before=radius_before,
                        dogleg_branch=dogleg_branch,
                        scaled_step_l2_norm=step_norm,
                        predicted_merit_reduction=predicted_reduction,
                        rejection_reason=(
                            "l2_merit_not_decreased"
                            if not math.isfinite(actual_reduction)
                            or actual_reduction <= 0.0
                            else "trust_ratio_below_acceptance_floor"
                        ),
                        residual_l2_norm=trial_norm,
                        maximum_scaled_residual=trial_maximum,
                        actual_merit_reduction=actual_reduction,
                        trust_ratio=trust_ratio,
                        accepted=False,
                        trust_radius_after=dogleg_trust_radius,
                    )
                )
            progress.direct_dogleg_iteration_audits.append(
                _DirectDoglegIterationAudit(
                    iteration=dogleg_iteration,
                    current_l2_norm=iteration_current_norm,
                    current_maximum_scaled_residual=(
                        iteration_current_maximum
                    ),
                    scaled_newton_correction_l2_norm=float(
                        np.linalg.norm(scaled_newton)
                    ),
                    scaled_cauchy_correction_l2_norm=cauchy_norm,
                    cauchy_gradient_l2_norm=cauchy_gradient_norm,
                    cauchy_numerator=cauchy_numerator,
                    cauchy_denominator=cauchy_denominator,
                    cauchy_alpha=cauchy_alpha,
                    trust_radius_on_entry=trust_radius_on_entry,
                    trials=tuple(dogleg_trials),
                    accepted=dogleg_accepted,
                    trust_radius_on_exit=dogleg_trust_radius,
                )
            )
            if not dogleg_accepted:
                raise _DirectStageBasinError(
                    "departure diagnostic Powell dogleg has no accepted "
                    "phase-safe trust-region step"
                )
    maximum = float(np.max(np.abs(current)))
    if maximum > controls.nonlinear_residual_tolerance:
        raise _DirectStageBasinError(
            "face-departure regional/RH residual contract failed: "
            f"{maximum:.3e} > {controls.nonlinear_residual_tolerance:.3e}"
        )
    signed_scaled = tuple(float(value) for value in current)
    certificate = None
    labels: tuple[str, ...] = ()
    physical_scales: tuple[float, ...] = ()
    perturbations: tuple[float, ...] = ()
    if certify:
        final_jacobian, coordinate_perturbations = _direct_coordinate_jacobian(
            coordinates,
            lower,
            upper,
            coordinate_scales,
            residual,
            layout,
            structural,
        )
        raw_condition = ci._jacobian_condition_proxy(  # noqa: SLF001
            final_jacobian
        )
        certificate, labels, physical_scales, perturbations = (
            _certify_affine_direct_jacobian(
                before,
                final_jacobian,
                coordinate_perturbations,
                scales,
                chart,
                scale_reference,
            )
        )
        progress.evaluations = evaluations
        progress.condition_proxy = certificate.condition_number_2
        progress.raw_condition_proxy = raw_condition
        progress.candidate = assembly.candidate
        progress.residual_blocks = assembly.residuals
        progress.signed_scaled_residuals = signed_scaled
        progress.maximum_scaled_residual = maximum
        if certificate.condition_number_2 > controls.maximum_condition_proxy:
            raise _DepartureConditioningError(
                "face-departure unit-invariant physical Jacobian is uncertified: "
                f"{certificate.condition_number_2:.3e}"
            )
    if certify and isinstance(boundary, ct.ReducedFilmPoreBoundary):
        surface_flux = assembly.dry_face_fluxes[-1]
        if surface_flux.surface_film_audit is None:
            raise RuntimeError("reduced-film departure lost its surface audit")
        enforced_audit = ct.enforce_reduced_film_surface_audit(
            surface_flux.surface_film_audit
        )
        assembly = replace(
            assembly,
            dry_face_fluxes=(
                *assembly.dry_face_fluxes[:-1],
                replace(surface_flux, surface_film_audit=enforced_audit),
            ),
        )
    progress.evaluations = evaluations
    progress.condition_proxy = (
        math.nan if certificate is None else certificate.condition_number_2
    )
    progress.raw_condition_proxy = raw_condition
    progress.candidate = assembly.candidate
    progress.residual_blocks = assembly.residuals
    progress.signed_scaled_residuals = signed_scaled
    progress.maximum_scaled_residual = maximum
    return _DirectStageSolution(
        assembly=assembly,
        newborn_trace_gradients=decoded.newborn_trace_gradients,
        scales=scales,
        signed_scaled_residuals=signed_scaled,
        evaluations=evaluations,
        rejected_trials=rejected_trials,
        message=(
            "sparse colored direct physical Newton root after diagnostic "
            "Powell-dogleg rescue; unchanged ordinary residual"
            if structural is not None
            and progress.dogleg_globalization_accepted_steps
            else "direct physical Newton root after diagnostic Powell-dogleg "
            "rescue; unchanged ordinary residual"
            if progress.dogleg_globalization_accepted_steps
            else
            "sparse colored direct physical Newton root with scaled-Cauchy "
            "globalization; unchanged ordinary residual"
            if structural is not None
            and progress.cauchy_globalization_accepted_steps
            else "direct physical Newton root with scaled-Cauchy globalization; "
            "unchanged ordinary residual"
            if progress.cauchy_globalization_accepted_steps
            else "sparse colored direct physical Newton root; unchanged ordinary residual"
            if structural is not None
            else "direct physical Newton root; unchanged ordinary residual"
        ),
        raw_condition_proxy=raw_condition,
        condition_certificate=certificate,
        physical_variable_scale_labels=labels,
        physical_variable_scales=physical_scales,
        physical_jacobian_perturbations=perturbations,
        chart=chart,
        structural=structural,
        cauchy_globalization_accepted_steps=(
            progress.cauchy_globalization_accepted_steps
        ),
        cauchy_line_search_trials=progress.cauchy_line_search_trials,
        cauchy_residual_evaluations=progress.cauchy_residual_evaluations,
        dogleg_globalization_accepted_steps=(
            progress.dogleg_globalization_accepted_steps
        ),
        dogleg_jacobian_linearizations=(
            progress.dogleg_jacobian_linearizations
        ),
        dogleg_trial_attempts=progress.dogleg_trial_attempts,
        dogleg_residual_evaluations=progress.dogleg_residual_evaluations,
    )


def _direct_coordinate_jacobian(
    coordinates: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    coordinate_scales: np.ndarray,
    residual,
    layout: departure.FaceDepartureLayout,
    structural: cut_sparsity.CutJacobianSparsity | None,
) -> tuple[scipy_sparse.csr_matrix | np.ndarray, np.ndarray]:
    if structural is not None:
        return _colored_direct_coordinate_jacobian(
            coordinates,
            lower,
            upper,
            coordinate_scales,
            residual,
            layout,
            structural,
        )
    matrix = np.empty((layout.residual_count, layout.unknown_count), dtype=float)
    perturbations = np.empty(layout.unknown_count, dtype=float)
    for column, characteristic in enumerate(coordinate_scales):
        room = min(
            coordinates[column] - lower[column],
            upper[column] - coordinates[column],
        )
        delta = min(1.0e-5 * characteristic, 0.2 * room)
        if not math.isfinite(delta) or delta <= 0.0:
            raise RuntimeError("departure Jacobian has no open central probe")
        for _ in range(14):
            upper_probe = coordinates.copy()
            lower_probe = coordinates.copy()
            upper_probe[column] += delta
            lower_probe[column] -= delta
            try:
                upper_residual, _, _ = residual(upper_probe)
                lower_residual, _, _ = residual(lower_probe)
            except (ValueError, departure.FaceDepartureStepError):
                delta *= 0.5
                continue
            matrix[:, column] = (
                upper_residual - lower_residual
            ) / (2.0 * delta)
            perturbations[column] = delta
            break
        else:
            raise RuntimeError(
                "departure Jacobian column has no phase-safe central probe"
            )
    return matrix, perturbations


def _colored_direct_coordinate_jacobian(
    coordinates: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    coordinate_scales: np.ndarray,
    residual,
    layout: departure.FaceDepartureLayout,
    structural: cut_sparsity.CutJacobianSparsity,
) -> tuple[scipy_sparse.csr_matrix, np.ndarray]:
    """Central-difference one structurally orthogonal column group at a time."""

    if structural.shape != (
        layout.residual_count,
        layout.unknown_count,
    ):
        raise RuntimeError("departure colored Jacobian has the wrong structure")
    column_rows: list[list[int]] = [
        [] for _ in range(layout.unknown_count)
    ]
    for row, columns in enumerate(structural.row_columns):
        for column in columns:
            column_rows[column].append(row)
    if any(not rows for rows in column_rows):
        raise RuntimeError("departure sparsity left an unknown column unsupported")
    matrix = scipy_sparse.lil_matrix(
        (layout.residual_count, layout.unknown_count),
        dtype=float,
    )
    perturbations = np.empty(layout.unknown_count, dtype=float)
    for columns in structural.column_groups:
        deltas = np.asarray(
            tuple(
                min(
                    1.0e-5 * coordinate_scales[column],
                    0.2
                    * min(
                        coordinates[column] - lower[column],
                        upper[column] - coordinates[column],
                    ),
                )
                for column in columns
            ),
            dtype=float,
        )
        if not np.all(np.isfinite(deltas)) or np.any(deltas <= 0.0):
            raise _DirectStageBasinError(
                "departure colored Jacobian has no open central probe"
            )
        for _ in range(14):
            upper_probe = coordinates.copy()
            lower_probe = coordinates.copy()
            for position, column in enumerate(columns):
                upper_probe[column] += deltas[position]
                lower_probe[column] -= deltas[position]
            try:
                upper_residual, _, _ = residual(upper_probe)
                lower_residual, _, _ = residual(lower_probe)
            except (ValueError, departure.FaceDepartureStepError):
                deltas *= 0.5
                continue
            difference = upper_residual - lower_residual
            for position, column in enumerate(columns):
                rows = column_rows[column]
                matrix[rows, column] = (
                    difference[rows] / (2.0 * deltas[position])
                )[:, None]
                perturbations[column] = deltas[position]
            break
        else:
            raise _DirectStageBasinError(
                "departure colored Jacobian group has no phase-safe central probe"
            )
    return matrix.tocsr(), perturbations


def _certify_departure_sparsity_pattern(
    independently_differenced_jacobian: np.ndarray,
    structural: cut_sparsity.CutJacobianSparsity,
) -> float:
    """Fail if an independently perturbed column escapes the declared support."""

    matrix = np.asarray(independently_differenced_jacobian, dtype=float)
    if matrix.shape != structural.shape or not np.all(np.isfinite(matrix)):
        raise ValueError("departure sparsity audit Jacobian is invalid")
    row_scales = np.maximum(np.max(np.abs(matrix), axis=1), 1.0e-300)
    normalized = np.abs(matrix) / row_scales[:, None]
    outside = np.where(structural.jac_sparsity.toarray(), 0.0, normalized)
    maximum = float(np.max(outside))
    if maximum > _SPARSITY_RELATIVE_ZERO_TOLERANCE:
        row, column = np.unravel_index(np.argmax(outside), outside.shape)
        raise RuntimeError(
            "departure Jacobian sparsity omitted a resolved derivative: "
            f"row={row}, column={column}, relative={maximum:.3e}"
        )
    return maximum


def _coordinate_chart(
    before: departure.FaceDepartureState,
    controls: ci.CutSolverControls,
    seed: departure.FaceDepartureUnknowns,
    seed_trace_gradients: departure.NewbornDryTraceGradients | None = None,
    dt_s: float = 1.0e-3,
    surface_boundary: ct.PoreBoundary | None = None,
    q_scale_m3_s: float | None = None,
) -> _CoordinateChart:
    """Build direct physical coordinates with explicit state-derived bounds."""

    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError("departure chart dt must be positive and finite")
    if surface_boundary is not None and not isinstance(
        surface_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("departure surface boundary must be a supported pore boundary")
    grid = before.geometry.master_grid
    k = before.departure_face_index
    inner_z = (grid.faces[k - 1] / grid.R) ** 3
    face_z = before.geometry.front.z
    if not inner_z < face_z:
        raise departure.FaceDepartureTopologyError(
            "departure chart has no positive z width"
        )
    if len(seed.dry_temperatures_k) != grid.n - k + 1 or len(
        seed.dry_y_hexane
    ) != grid.n - k + 1:
        raise ValueError("departure chart seed has the wrong dry-piece rank")
    dry_temperature_band = (
        before.config.dry.conditioned_temperature_domain.solver_bounds_k
    )
    numerical_y_bounds = (
        None
        if isinstance(before.config.dry.primitive_band, ct.GasOnlyPrimitiveDomain)
        else before.config.dry.primitive_band.y_hexane_bounds
    )
    dry_pores = tuple(
        replace(before.config.dry.pore, w_o=before.oil_fraction_labels[index])
        for index in range(k - 1, grid.n)
    )
    conditioned_chart = ci._CoordinateChart(  # noqa: SLF001
        wet_temperature=controls.wet_temperature_bounds_k,
        wet_water=controls.wet_water_bounds,
        # This event chart currently fails closed on a nonzero wet cap dual;
        # it must not silently reinterpret an active graph coordinate.
        wet_retained_cap=None,
        dry_temperature=dry_temperature_band,
        dry_y_hexane=numerical_y_bounds,
        dry_config=before.config.dry,
        dry_pores=dry_pores,
        surface_boundary=surface_boundary,
        front_z=(inner_z, face_z),
        interface_temperature=controls.interface_temperature_bounds_k,
        moving_interface_composition_force_authority=(
            before.config.moving_interface_composition_force_authority
        ),
    )
    seed_geometry = cg.partition_master_grid(grid, z=seed.front_z)
    seed_layout = cut.layout_for_geometry(seed_geometry)
    seed_distance = _newborn_center_distance(
        seed_geometry,
        seed_layout.cut_cell_index,
    )
    seed_gradient = (
        seed_trace_gradients.temperature_gradient_k_m
        if seed_trace_gradients is not None
        else (
            seed.dry_temperatures_k[0] - seed.interface_temperature_k
        )
        / seed_distance
    )
    if seed_trace_gradients is not None and math.fsum(
        (seed.interface_temperature_k, seed_gradient * seed_distance)
    ) != seed.dry_temperatures_k[0]:
        raise departure.FaceDepartureTopologyError(
            "departure seed temperature is inconsistent with its trace gradient"
        )
    temperature_gradient_bounds = _newborn_temperature_gradient_bounds(
        seed.interface_temperature_k,
        seed_distance,
        dry_temperature_band,
    )
    q_scale = (
        grid.total_volume * math.fsum((face_z, -seed.front_z)) / dt_s
        if q_scale_m3_s is None
        else q_scale_m3_s
    )
    q_max = grid.volumes[k - 1] / dt_s
    if not 0.0 < q_scale < q_max:
        raise departure.FaceDepartureTopologyError(
            "departure seed must have a positive sub-cell newborn volume rate"
        )
    interface = cut.evaluate_interface_state(
        seed.interface_temperature_k,
        before.config,
        before.historical_hexane_loadings[k - 1],
        before.oil_fraction_labels[k - 1],
    )
    seed_y_gradient = (
        seed_trace_gradients.composition_gradient_m_inv
        if seed_trace_gradients is not None
        else (seed.dry_y_hexane[0] - interface.y_hexane) / seed_distance
    )
    if seed_trace_gradients is not None and math.fsum(
        (interface.y_hexane, seed_y_gradient * seed_distance)
    ) != seed.dry_y_hexane[0]:
        raise departure.FaceDepartureTopologyError(
            "departure seed composition is inconsistent with its trace gradient"
        )
    newborn_contexts = _newborn_composition_context_temperatures(
        seed.dry_temperatures_k,
        seed.interface_temperature_k,
        before.config.moving_interface_composition_force_authority,
    )
    newborn_y_bounds = _newborn_conditioned_interval(
        newborn_contexts,
        seed.interface_temperature_k,
        before.config.dry.pressure_pa,
        dry_pores[0],
        numerical_y_bounds,
        before.config.moving_interface_composition_force_authority,
    )
    composition_gradient_bounds = (
        (newborn_y_bounds[0] - interface.y_hexane) / seed_distance,
        (newborn_y_bounds[1] - interface.y_hexane) / seed_distance,
    )
    if not composition_gradient_bounds[0] < seed_y_gradient < composition_gradient_bounds[1]:
        raise departure.FaceDepartureTopologyError(
            "departure seed composition gradient is outside its three-context domain"
        )
    outer_rate_bounds: list[tuple[float, float]] = []
    for position, reference_y in enumerate(before.dry_y_hexane, start=1):
        contexts = ci._dry_composition_context_temperatures(  # noqa: SLF001
            position,
            seed.dry_temperatures_k,
            seed.interface_temperature_k,
            conditioned_chart,
        )
        if numerical_y_bounds is None:
            y_bounds = ci._face_conditioned_composition_interval(  # noqa: SLF001
                contexts,
                before.config.dry.pressure_pa,
                dry_pores[position],
            )
        else:
            y_bounds = numerical_y_bounds
        strict_y = _strict_bounds(y_bounds)
        outer_rate_bounds.append(
            (
                (strict_y[0] - reference_y) / dt_s,
                (strict_y[1] - reference_y) / dt_s,
            )
        )
    return _CoordinateChart(
        dt_s=dt_s,
        wet_temperature=controls.wet_temperature_bounds_k,
        wet_water=controls.wet_water_bounds,
        dry_temperature=dry_temperature_band,
        dry_y_hexane=numerical_y_bounds,
        front_z=(inner_z, face_z),
        interface_temperature=controls.interface_temperature_bounds_k,
        reference_wet_temperatures_k=before.wet_temperatures_k,
        reference_wet_water_loadings=before.wet_retained_water_loadings,
        reference_outer_dry_temperatures_k=before.dry_temperatures_k,
        reference_outer_dry_y_hexane=before.dry_y_hexane,
        dry_pores=dry_pores,
        face_conditioned_chart=conditioned_chart,
        newborn_temperature_gradient_bounds_k_m=temperature_gradient_bounds,
        newborn_composition_gradient_bounds_m_inv=composition_gradient_bounds,
        outer_y_rate_bounds_s_inv=tuple(outer_rate_bounds),
        q_scale_m3_s=q_scale,
        q_ratio_bounds=(0.0, q_max / q_scale),
    )


def _validate_chart(
    before: departure.FaceDepartureState,
    controls: ci.CutSolverControls,
    chart: _CoordinateChart,
) -> None:
    k = before.departure_face_index - 1
    for temperature in (
        controls.interface_temperature_bounds_k[0],
        0.5 * math.fsum(controls.interface_temperature_bounds_k),
        controls.interface_temperature_bounds_k[1],
    ):
        cut.evaluate_interface_state(
            temperature,
            before.config,
            before.historical_hexane_loadings[k],
            before.oil_fraction_labels[k],
        )
    if not chart.front_z[0] < chart.front_z[1]:
        raise ValueError("departure front chart must remain open and ordered")
    if not math.isfinite(chart.q_scale_m3_s) or chart.q_scale_m3_s <= 0.0:
        raise ValueError("departure Q scale must be positive")
    for bounds in (
        chart.newborn_temperature_gradient_bounds_k_m,
        chart.newborn_composition_gradient_bounds_m_inv,
        *chart.outer_y_rate_bounds_s_inv,
    ):
        if not all(math.isfinite(value) for value in bounds) or not bounds[0] < bounds[1]:
            raise ValueError("departure direct-coordinate bounds are unresolved")


def _validate_seed(
    seed: FaceDepartureSeed,
    before: departure.FaceDepartureState,
    layout: departure.FaceDepartureLayout,
    chart: _CoordinateChart,
) -> None:
    if not isinstance(seed, FaceDepartureSeed):
        raise TypeError("seed must be FaceDepartureSeed")
    if len(seed.candidate.vector()) != layout.unknown_count:
        raise ValueError("departure seed does not match the nonlinear rank")
    if len(seed.stefan_flux_scales_mol_m2_s) != layout.dry_face_count:
        raise ValueError("one departure Stefan scale is required per dry face")
    _encode_candidate(
        seed.candidate,
        before,
        chart,
        seed.stefan_flux_scales_mol_m2_s,
        seed.newborn_trace_gradients,
    )


def _validate_direct_candidate_contexts(
    candidate: departure.FaceDepartureUnknowns,
    before: departure.FaceDepartureState,
    chart: _CoordinateChart,
) -> None:
    """Fail closed on every primitive and every temperature used by fluxes."""

    if not all(ci._inside(value, chart.wet_temperature) for value in candidate.wet_temperatures_k):  # noqa: SLF001
        raise departure.FaceDepartureTopologyError(
            "wet temperature left its direct-coordinate authority"
        )
    if not all(ci._inside(value, chart.wet_water) for value in candidate.wet_retained_water_loadings):  # noqa: SLF001
        raise departure.FaceDepartureTopologyError(
            "wet water left its direct-coordinate authority"
        )
    if not all(ci._inside(value, chart.dry_temperature) for value in candidate.dry_temperatures_k):  # noqa: SLF001
        raise departure.FaceDepartureTopologyError(
            "dry temperature left its conditioned authority"
        )
    if not ci._inside(candidate.interface_temperature_k, chart.interface_temperature):  # noqa: SLF001
        raise departure.FaceDepartureTopologyError(
            "interface temperature left its direct-coordinate authority"
        )
    newborn_contexts = _newborn_composition_context_temperatures(
        candidate.dry_temperatures_k,
        candidate.interface_temperature_k,
        before.config.moving_interface_composition_force_authority,
    )
    _newborn_conditioned_coordinate_and_derivative(
        newborn_contexts,
        candidate.interface_temperature_k,
        candidate.dry_y_hexane[0],
        before.config.dry.pressure_pa,
        chart.dry_pores[0],
        chart.dry_y_hexane,
        before.config.moving_interface_composition_force_authority,
    )
    for position, composition in enumerate(candidate.dry_y_hexane[1:], start=1):
        contexts = ci._dry_composition_context_temperatures(  # noqa: SLF001
            position,
            candidate.dry_temperatures_k,
            candidate.interface_temperature_k,
            chart.face_conditioned_chart,
        )
        _outer_conditioned_coordinate_and_derivative(
            contexts,
            composition,
            chart.face_conditioned_chart,
            position,
        )


def _encode_candidate(
    candidate: departure.FaceDepartureUnknowns,
    before: departure.FaceDepartureState,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
    newborn_trace_gradients: departure.NewbornDryTraceGradients | None = None,
) -> np.ndarray:
    """Encode direct rates/gradients, normalized fluxes/Q, and direct TΓ."""

    geometry = cg.partition_master_grid(
        before.geometry.master_grid,
        z=candidate.front_z,
    )
    layout = cut.layout_for_geometry(geometry)
    newborn_distance = _newborn_center_distance(geometry, layout.cut_cell_index)
    newborn_temperature_gradient = (
        newborn_trace_gradients.temperature_gradient_k_m
        if newborn_trace_gradients is not None
        else (
            candidate.dry_temperatures_k[0]
            - candidate.interface_temperature_k
        )
        / newborn_distance
    )
    if newborn_trace_gradients is not None and math.fsum(
        (
            candidate.interface_temperature_k,
            newborn_temperature_gradient * newborn_distance,
        )
    ) != candidate.dry_temperatures_k[0]:
        raise departure.FaceDepartureTopologyError(
            "encoded newborn temperature lost its trace gradient"
        )
    interface = cut.evaluate_interface_state(
        candidate.interface_temperature_k,
        before.config,
        before.historical_hexane_loadings[layout.cut_cell_index],
        before.oil_fraction_labels[layout.cut_cell_index],
    )
    newborn_y_gradient = (
        newborn_trace_gradients.composition_gradient_m_inv
        if newborn_trace_gradients is not None
        else (
            candidate.dry_y_hexane[0] - interface.y_hexane
        )
        / newborn_distance
    )
    if newborn_trace_gradients is not None and math.fsum(
        (interface.y_hexane, newborn_y_gradient * newborn_distance)
    ) != candidate.dry_y_hexane[0]:
        raise departure.FaceDepartureTopologyError(
            "encoded newborn composition lost its trace gradient"
        )
    _validate_direct_candidate_contexts(candidate, before, chart)
    outer_y_rates = tuple(
        (value - reference) / chart.dt_s
        for value, reference in zip(
            candidate.dry_y_hexane[1:],
            chart.reference_outer_dry_y_hexane,
        )
    )
    q_m3_s = (
        before.geometry.master_grid.total_volume
        * math.fsum((chart.front_z[1], -candidate.front_z))
        / chart.dt_s
    )
    values = (
        *(
            (value - reference) / chart.dt_s
            for value, reference in zip(
                candidate.wet_temperatures_k,
                chart.reference_wet_temperatures_k,
            )
        ),
        *(
            (value - reference) / chart.dt_s
            for value, reference in zip(
                candidate.wet_retained_water_loadings,
                chart.reference_wet_water_loadings,
            )
        ),
        newborn_temperature_gradient,
        *(
            (value - reference) / chart.dt_s
            for value, reference in zip(
                candidate.dry_temperatures_k[1:],
                chart.reference_outer_dry_temperatures_k,
            )
        ),
        newborn_y_gradient,
        *outer_y_rates,
        *(
            value / scale
            for value, scale in zip(
                candidate.dry_total_stefan_fluxes_mol_m2_s,
                stefan_scales,
            )
        ),
        q_m3_s / chart.q_scale_m3_s,
        candidate.interface_temperature_k,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("encoded departure coordinates must remain finite")
    return np.asarray(values, dtype=float)


def _decode_trial(
    coordinates: Sequence[float],
    before: departure.FaceDepartureState,
    layout: departure.FaceDepartureLayout,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
) -> _DecodedTrial:
    values = tuple(float(value) for value in coordinates)
    if len(values) != layout.unknown_count or not all(
        math.isfinite(value) for value in values
    ):
        raise ValueError("departure coordinate vector is non-finite or wrong-rank")
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
        math.fsum((reference, chart.dt_s * rate))
        for reference, rate in zip(
            chart.reference_wet_temperatures_k,
            take(nw),
        )
    )
    wet_w = tuple(
        math.fsum((reference, chart.dt_s * rate))
        for reference, rate in zip(
            chart.reference_wet_water_loadings,
            take(nw),
        )
    )
    newborn_gradient = take(1)[0]
    outer_temperature_rates = take(nd - 1)
    newborn_y_gradient = take(1)[0]
    outer_y_rates = take(nd - 1)
    nt = tuple(value * scale for value, scale in zip(take(nf), stefan_scales))
    q_ratio = take(1)[0]
    interface_t = take(1)[0]
    if not ci._inside(interface_t, chart.interface_temperature):  # noqa: SLF001
        raise departure.FaceDepartureTopologyError(
            "departure interface temperature left its open authority"
        )
    if not ci._inside(q_ratio, chart.q_ratio_bounds):  # noqa: SLF001
        raise departure.FaceDepartureTopologyError(
            "departure Q ratio left its open geometric interval"
        )
    q_m3_s = chart.q_scale_m3_s * q_ratio
    front_z = math.fsum(
        (
            chart.front_z[1],
            -q_m3_s * chart.dt_s / before.geometry.master_grid.total_volume,
        )
    )
    geometry = cg.partition_master_grid(before.geometry.master_grid, z=front_z)
    cut_layout = cut.layout_for_geometry(geometry)
    newborn_distance = _newborn_center_distance(
        geometry,
        cut_layout.cut_cell_index,
    )
    newborn_temperature = math.fsum(
        (interface_t, newborn_gradient * newborn_distance)
    )
    dry_t = (
        newborn_temperature,
        *(
            math.fsum((reference, chart.dt_s * rate))
            for reference, rate in zip(
                chart.reference_outer_dry_temperatures_k,
                outer_temperature_rates,
            )
        ),
    )
    interface = cut.evaluate_interface_state(
        interface_t,
        before.config,
        before.historical_hexane_loadings[cut_layout.cut_cell_index],
        before.oil_fraction_labels[cut_layout.cut_cell_index],
    )
    newborn_y = math.fsum(
        (interface.y_hexane, newborn_y_gradient * newborn_distance)
    )
    outer_y = tuple(
        math.fsum((reference, chart.dt_s * rate))
        for reference, rate in zip(
            chart.reference_outer_dry_y_hexane,
            outer_y_rates,
        )
    )
    candidate = departure.FaceDepartureUnknowns(
        wet_temperatures_k=wet_t,
        wet_retained_water_loadings=wet_w,
        dry_temperatures_k=dry_t,
        dry_y_hexane=(newborn_y, *outer_y),
        dry_total_stefan_fluxes_mol_m2_s=nt,
        front_z=front_z,
        interface_temperature_k=interface_t,
    )
    decoded = _DecodedTrial(
        candidate,
        departure.NewbornDryTraceGradients(
            newborn_gradient,
            newborn_y_gradient,
        ),
    )
    _validate_direct_candidate_contexts(candidate, before, chart)
    return decoded


def _coordinate_solver_bounds(
    layout: departure.FaceDepartureLayout,
    chart: _CoordinateChart,
    maximum_logit: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return direct physical guards in documented nonlinear block order."""

    del maximum_logit

    wet_t_bounds = tuple(
        _rate_bounds(reference, chart.wet_temperature, chart.dt_s)
        for reference in chart.reference_wet_temperatures_k
    )
    wet_w_bounds = tuple(
        _rate_bounds(reference, chart.wet_water, chart.dt_s)
        for reference in chart.reference_wet_water_loadings
    )
    gradient_bounds = _strict_bounds(
        chart.newborn_temperature_gradient_bounds_k_m
    )
    outer_t_bounds = tuple(
        _rate_bounds(reference, chart.dry_temperature, chart.dt_s)
        for reference in chart.reference_outer_dry_temperatures_k
    )
    newborn_y_guard = _strict_bounds(
        chart.newborn_composition_gradient_bounds_m_inv
    )
    outer_y_guards = chart.outer_y_rate_bounds_s_inv
    q_guard = _strict_bounds(chart.q_ratio_bounds)
    interface_guard = _strict_bounds(chart.interface_temperature)
    blocks = (
        *wet_t_bounds,
        *wet_w_bounds,
        gradient_bounds,
        *outer_t_bounds,
        newborn_y_guard,
        *outer_y_guards,
        *((-math.inf, math.inf),) * layout.dry_face_count,
        q_guard,
        interface_guard,
    )
    if len(blocks) != layout.unknown_count:
        raise RuntimeError("departure coalescing solver bounds lost block rank")
    return (
        np.asarray(tuple(bounds[0] for bounds in blocks), dtype=float),
        np.asarray(tuple(bounds[1] for bounds in blocks), dtype=float),
    )


def _solver_coordinate_characteristic_scales(
    before: departure.FaceDepartureState,
    chart: _CoordinateChart,
    scale_reference: _PhysicalScaleReference,
) -> tuple[float, ...]:
    """Declared absolute differencing scales in direct coordinates."""

    layout = departure.layout_for_departure(before)
    wet_t_width = chart.wet_temperature[1] - chart.wet_temperature[0]
    wet_w_width = chart.wet_water[1] - chart.wet_water[0]
    dry_t_width = chart.dry_temperature[1] - chart.dry_temperature[0]
    values = (
        *((wet_t_width,) * layout.wet_piece_count),
        *((wet_w_width,) * layout.wet_piece_count),
        scale_reference.newborn_temperature_gradient_scale_k_m,
        *((dry_t_width,) * (layout.dry_piece_count - 1)),
        scale_reference.newborn_composition_gradient_scale_m_inv,
        *((1.0,) * (layout.dry_piece_count - 1)),
        *((1.0,) * layout.dry_face_count),
        1.0,
        5.0,
    )
    if len(values) != layout.unknown_count or not all(
        math.isfinite(value) and value > 0.0 for value in values
    ):
        raise RuntimeError("departure coordinate differencing scales lost rank")
    return values


def _strict_bounds(bounds: tuple[float, float]) -> tuple[float, float]:
    lower, upper = bounds
    result = (math.nextafter(lower, upper), math.nextafter(upper, lower))
    if not result[0] < result[1]:
        raise ValueError("open authority has no representable strict interior")
    return result


def _rate_bounds(
    reference: float,
    primitive_bounds: tuple[float, float],
    dt_s: float,
) -> tuple[float, float]:
    lower, upper = _strict_bounds(primitive_bounds)
    result = ((lower - reference) / dt_s, (upper - reference) / dt_s)
    if not all(math.isfinite(value) for value in result) or not result[0] < result[1]:
        raise ValueError("existing-cell rate bounds are unresolved")
    return result


def _bounded_fractions(
    candidate: departure.FaceDepartureUnknowns,
    before: departure.FaceDepartureState,
    chart: _CoordinateChart,
) -> tuple[float, ...]:
    newborn_contexts = _newborn_composition_context_temperatures(
        candidate.dry_temperatures_k,
        candidate.interface_temperature_k,
        before.config.moving_interface_composition_force_authority,
    )
    newborn_y_bounds = _newborn_conditioned_interval(
        newborn_contexts,
        candidate.interface_temperature_k,
        before.config.dry.pressure_pa,
        chart.dry_pores[0],
        chart.dry_y_hexane,
        before.config.moving_interface_composition_force_authority,
    )
    dry_y_pairs = [(candidate.dry_y_hexane[0], newborn_y_bounds)]
    for position, composition in enumerate(candidate.dry_y_hexane[1:], start=1):
        if chart.dry_y_hexane is None:
            contexts = ci._dry_composition_context_temperatures(  # noqa: SLF001
                position,
                candidate.dry_temperatures_k,
                candidate.interface_temperature_k,
                chart.face_conditioned_chart,
            )
            bounds = ci._face_conditioned_composition_interval(  # noqa: SLF001
                contexts,
                before.config.dry.pressure_pa,
                chart.dry_pores[position],
            )
        else:
            bounds = chart.dry_y_hexane
        dry_y_pairs.append((composition, bounds))
    pairs = (
        *((value, chart.wet_temperature) for value in candidate.wet_temperatures_k),
        *((value, chart.wet_water) for value in candidate.wet_retained_water_loadings),
        *((value, chart.dry_temperature) for value in candidate.dry_temperatures_k),
        *dry_y_pairs,
        (candidate.front_z, chart.front_z),
        (candidate.interface_temperature_k, chart.interface_temperature),
    )
    result = tuple(
        (value - bounds[0]) / (bounds[1] - bounds[0])
        for value, bounds in pairs
    )
    if not all(0.0 < value < 1.0 for value in result):
        raise RuntimeError("accepted departure primitive touched a chart boundary")
    return result


def _newborn_center_distance(geometry: cg.CutGeometry, cell_index: int) -> float:
    distance = math.fsum(
        (
            cut._dry_piece_center(  # noqa: SLF001
                geometry.cells[cell_index],
                geometry.front.radius_m,
            ),
            -geometry.front.radius_m,
        )
    )
    if not math.isfinite(distance) or distance <= 0.0:
        raise departure.FaceDepartureTopologyError(
            "newborn dry-cell center must remain strictly outside the front"
        )
    return distance


def _newborn_temperature_gradient_bounds(
    interface_temperature_k: float,
    center_distance_m: float,
    dry_temperature_bounds_k: tuple[float, float],
) -> tuple[float, float]:
    """Open gradient interval keeping the newborn temperature in authority.

    The moving-interface force is reconstructed at ``T_Gamma`` itself, so
    gas-only topology does not require a nonnegative dry-side temperature
    gradient.  Coupled ``T-y`` paths with ``T_newborn < T_Gamma`` remain
    admissible when their composition moves into the gas-only interior.
    """

    t_lo, t_hi = dry_temperature_bounds_k
    bounds = (
        (t_lo - interface_temperature_k) / center_distance_m,
        (t_hi - interface_temperature_k) / center_distance_m,
    )
    if not all(math.isfinite(value) for value in bounds) or not (
        bounds[0] < bounds[1]
    ):
        raise departure.FaceDepartureTopologyError(
            "newborn dry-temperature gradient chart is unresolved"
        )
    return bounds


def _newborn_composition_context_temperatures(
    dry_temperatures_k: Sequence[float],
    interface_temperature_k: float,
    authority: cut.MovingInterfaceCompositionForceAuthority,
) -> tuple[float, ...]:
    """Return every common-temperature context used by the newborn state.

    At an interior master-face departure there is always an existing dry cell
    outside the newborn one.  The first ordinary dry/dry face therefore uses
    the symmetric temperature ``(T_new + T_outer)/2`` and must participate in
    the newborn composition chart from its first positive-volume trial.  The
    legacy moving face additionally uses ``T_Gamma``.  Amendment 20's actual
    path owns that face-domain certificate and therefore omits only this
    virtual common-temperature context.
    """

    temperatures = tuple(float(value) for value in dry_temperatures_k)
    if len(temperatures) < 2:
        raise departure.FaceDepartureTopologyError(
            "inner-face departure needs an outer dry cell for its first face"
        )
    if not isinstance(authority, cut.MovingInterfaceCompositionForceAuthority):
        raise TypeError("newborn composition contexts require a named authority")
    internal_temperature = 0.5 * math.fsum((temperatures[0], temperatures[1]))
    if authority is cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH:
        result = (temperatures[0], internal_temperature)
    else:
        result = (
            temperatures[0],
            float(interface_temperature_k),
            internal_temperature,
        )
    if not all(math.isfinite(value) and value > 0.0 for value in result):
        raise departure.FaceDepartureTopologyError(
            "newborn composition contexts must be positive and finite"
        )
    return result


def _newborn_conditioned_interval(
    context_temperatures_k: Sequence[float],
    interface_temperature_k: float,
    pressure_pa: float,
    pore: cp.CoupledPoreParams,
    numerical_bounds: tuple[float, float] | None,
    authority: cut.MovingInterfaceCompositionForceAuthority,
) -> tuple[float, float]:
    physical_lower, physical_upper = (
        ci._face_conditioned_composition_interval(  # noqa: SLF001
            context_temperatures_k,
            pressure_pa,
            pore,
        )
    )
    if authority is cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH:
        lower = physical_lower
    else:
        smooth_lower = math.nextafter(
            cp.smooth_luikov_composition_lower_bound(
                interface_temperature_k,
                pressure_pa,
                pore,
            ),
            physical_upper,
        )
        lower = max(physical_lower, smooth_lower)
    upper = physical_upper
    if numerical_bounds is not None:
        lower = max(lower, numerical_bounds[0])
        upper = min(upper, numerical_bounds[1])
    if not lower < upper:
        raise departure.FaceDepartureTopologyError(
            "newborn storage/interface/internal-face composition domains do not overlap"
        )
    return lower, upper


def _newborn_conditioned_coordinate_and_derivative(
    context_temperatures_k: Sequence[float],
    interface_temperature_k: float,
    y_hexane: float,
    pressure_pa: float,
    pore: cp.CoupledPoreParams,
    numerical_bounds: tuple[float, float] | None,
    authority: cut.MovingInterfaceCompositionForceAuthority,
) -> tuple[float, float]:
    r"""Smooth strictly increasing lift for the newborn ``g_y`` unknown.

    The exact gas-only barriers for every common-temperature evaluation
    context are summed, so a change in which context is tight does not switch
    coordinates.  The legacy authority also places its lower open-logit bound
    on the approved smooth Luikov side.  Amendment 20 permits a resolved
    retained-cap crossing and leaves that path topology to the actual-path
    operator.  This numerical chart adds no flux, storage, or phase mechanism.
    """

    lower, upper = _newborn_conditioned_interval(
        context_temperatures_k,
        interface_temperature_k,
        pressure_pa,
        pore,
        numerical_bounds,
        authority,
    )
    if not math.isfinite(y_hexane) or not lower < y_hexane < upper:
        raise departure.FaceDepartureTopologyError(
            "newborn composition is infeasible in an evaluated dry context"
        )
    exact_coordinate, exact_derivative = (
        ci._face_conditioned_coordinate_and_derivative(  # noqa: SLF001
            context_temperatures_k,
            pressure_pa,
            y_hexane,
            pore,
        )
    )
    seam_coordinate = ci._encode_open(y_hexane, (lower, upper))  # noqa: SLF001
    seam_derivative = 1.0 / (y_hexane - lower) + 1.0 / (upper - y_hexane)
    coordinate = math.fsum((exact_coordinate, seam_coordinate))
    derivative = math.fsum((exact_derivative, seam_derivative))
    if not math.isfinite(coordinate) or not math.isfinite(derivative) or derivative <= 0.0:
        raise departure.FaceDepartureTopologyError(
            "newborn smooth composition barrier lost strict monotonicity"
        )
    return coordinate, derivative


def _outer_conditioned_coordinate_and_derivative(
    context_temperatures_k: Sequence[float],
    y_hexane: float,
    chart: ci._CoordinateChart,  # noqa: SLF001
    piece: int,
) -> tuple[float, float]:
    if chart.dry_y_hexane is None:
        return ci._face_conditioned_coordinate_and_derivative(  # noqa: SLF001
            context_temperatures_k,
            chart.dry_config.pressure_pa,
            y_hexane,
            chart.dry_pores[piece],
        )
    coordinate = ci._encode_open(y_hexane, chart.dry_y_hexane)  # noqa: SLF001
    lower, upper = chart.dry_y_hexane
    derivative = 1.0 / (y_hexane - lower) + 1.0 / (upper - y_hexane)
    return coordinate, derivative


def _physical_departure_vector(
    candidate: departure.FaceDepartureUnknowns,
    gradients: departure.NewbornDryTraceGradients,
    before: departure.FaceDepartureState,
    dt_s: float,
) -> tuple[float, ...]:
    """Direct physical variables used only by the conditioning audit."""

    q_m3_s = (
        before.geometry.master_grid.total_volume
        * math.fsum((before.geometry.front.z, -candidate.front_z))
        / dt_s
    )
    values = (
        *((new - old) / dt_s for new, old in zip(
            candidate.wet_temperatures_k,
            before.wet_temperatures_k,
        )),
        *((new - old) / dt_s for new, old in zip(
            candidate.wet_retained_water_loadings,
            before.wet_retained_water_loadings,
        )),
        gradients.temperature_gradient_k_m,
        *((new - old) / dt_s for new, old in zip(
            candidate.dry_temperatures_k[1:],
            before.dry_temperatures_k,
        )),
        gradients.composition_gradient_m_inv,
        *((new - old) / dt_s for new, old in zip(
            candidate.dry_y_hexane[1:],
            before.dry_y_hexane,
        )),
        *candidate.dry_total_stefan_fluxes_mol_m2_s,
        q_m3_s,
        candidate.interface_temperature_k,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("physical departure vector must remain finite")
    return values


def _decode_physical_departure_vector(
    values: Sequence[float],
    before: departure.FaceDepartureState,
    layout: departure.FaceDepartureLayout,
    dt_s: float,
) -> _DecodedTrial:
    """Decode the audit variables without using solver-coordinate scaling."""

    vector = tuple(float(value) for value in values)
    if len(vector) != layout.unknown_count or not all(
        math.isfinite(value) for value in vector
    ):
        raise ValueError("physical departure vector is non-finite or wrong-rank")
    nw = layout.wet_piece_count
    nd = layout.dry_piece_count
    nf = layout.dry_face_count
    cursor = 0

    def take(count: int) -> tuple[float, ...]:
        nonlocal cursor
        result = vector[cursor : cursor + count]
        cursor += count
        return result

    wet_t_rates = take(nw)
    wet_w_rates = take(nw)
    g_t = take(1)[0]
    outer_t_rates = take(nd - 1)
    g_y = take(1)[0]
    outer_y_rates = take(nd - 1)
    nt = take(nf)
    q_m3_s = take(1)[0]
    interface_t = take(1)[0]
    q_max = before.geometry.master_grid.volumes[layout.cut_cell_index] / dt_s
    if not 0.0 < q_m3_s < q_max:
        raise departure.FaceDepartureTopologyError(
            "physical Jacobian probe left the open newborn-volume interval"
        )
    front_z = math.fsum(
        (
            before.geometry.front.z,
            -q_m3_s * dt_s / before.geometry.master_grid.total_volume,
        )
    )
    geometry = cg.partition_master_grid(before.geometry.master_grid, z=front_z)
    cut_layout = cut.layout_for_geometry(geometry)
    distance = _newborn_center_distance(geometry, cut_layout.cut_cell_index)
    wet_t = tuple(
        math.fsum((old, dt_s * rate))
        for old, rate in zip(before.wet_temperatures_k, wet_t_rates)
    )
    wet_w = tuple(
        math.fsum((old, dt_s * rate))
        for old, rate in zip(before.wet_retained_water_loadings, wet_w_rates)
    )
    outer_t = tuple(
        math.fsum((old, dt_s * rate))
        for old, rate in zip(before.dry_temperatures_k, outer_t_rates)
    )
    newborn_t = math.fsum((interface_t, g_t * distance))
    dry_t = (newborn_t, *outer_t)
    interface = cut.evaluate_interface_state(
        interface_t,
        before.config,
        before.historical_hexane_loadings[layout.cut_cell_index],
        before.oil_fraction_labels[layout.cut_cell_index],
    )
    newborn_y = math.fsum((interface.y_hexane, g_y * distance))
    outer_y = tuple(
        math.fsum((old, dt_s * rate))
        for old, rate in zip(before.dry_y_hexane, outer_y_rates)
    )
    return _DecodedTrial(
        departure.FaceDepartureUnknowns(
            wet_temperatures_k=wet_t,
            wet_retained_water_loadings=wet_w,
            dry_temperatures_k=dry_t,
            dry_y_hexane=(newborn_y, *outer_y),
            dry_total_stefan_fluxes_mol_m2_s=nt,
            front_z=front_z,
            interface_temperature_k=interface_t,
        ),
        departure.NewbornDryTraceGradients(g_t, g_y),
    )


def _physical_variable_scale_contract(
    before: departure.FaceDepartureState,
    chart: _CoordinateChart,
    scale_reference: _PhysicalScaleReference,
) -> tuple[tuple[str, ...], tuple[float, ...]]:
    """Declare reproducible ``D_X`` entries for the physical Jacobian.

    Existing-state rates use the previously declared one-second
    characteristic-rate convention: their numerical scales equal the
    authorized feasibility widths per one second.  This is fixed across dt
    refinement.  Converting the time unit transforms both the rate variables
    and these scales, leaving ``J_star`` invariant.  Every seam scale comes
    from the unique exact-face tangent with fixed documented floors.
    ``T_Gamma`` uses the fixed 5 K interface contract.  No entry depends on
    the caller seed, the converged root, or its Jacobian.
    """

    layout = departure.layout_for_departure(before)
    wet_t_width = chart.wet_temperature[1] - chart.wet_temperature[0]
    wet_w_width = chart.wet_water[1] - chart.wet_water[0]
    dry_t_width = chart.dry_temperature[1] - chart.dry_temperature[0]
    labels = (
        *(f"wet_temperature_rate[{index}]" for index in range(layout.wet_piece_count)),
        *(f"wet_water_rate[{index}]" for index in range(layout.wet_piece_count)),
        "newborn_temperature_gradient",
        *(f"outer_dry_temperature_rate[{index}]" for index in range(layout.dry_piece_count - 1)),
        "newborn_composition_gradient",
        *(f"outer_dry_composition_rate[{index}]" for index in range(layout.dry_piece_count - 1)),
        *(f"total_stefan_flux[{index}]" for index in range(layout.dry_face_count)),
        "newborn_volume_rate_Q",
        "interface_temperature",
    )
    scales = (
        *((wet_t_width,) * layout.wet_piece_count),
        *((wet_w_width,) * layout.wet_piece_count),
        scale_reference.newborn_temperature_gradient_scale_k_m,
        *((dry_t_width,) * (layout.dry_piece_count - 1)),
        scale_reference.newborn_composition_gradient_scale_m_inv,
        *((1.0,) * (layout.dry_piece_count - 1)),
        *scale_reference.stefan_flux_scales_mol_m2_s,
        scale_reference.q_scale_m3_s,
        5.0,
    )
    if len(labels) != layout.unknown_count or len(scales) != layout.unknown_count:
        raise RuntimeError("physical Jacobian scale contract lost block rank")
    if not all(math.isfinite(value) and value > 0.0 for value in scales):
        raise ValueError("physical Jacobian scales must be positive and finite")
    return tuple(labels), tuple(scales)


def _certify_affine_direct_jacobian(
    before: departure.FaceDepartureState,
    direct_jacobian: scipy_sparse.spmatrix | np.ndarray,
    coordinate_perturbations: np.ndarray,
    residual_scales: ci.CutResidualScales,
    chart: _CoordinateChart,
    scale_reference: _PhysicalScaleReference,
) -> tuple[
    conditioning.UnitInvariantJacobianCertificate,
    tuple[str, ...],
    tuple[float, ...],
    tuple[float, ...],
]:
    r"""Certify the final physical Jacobian without duplicating its probes.

    The active corrector chart is affine: rates and gradients are already
    physical, while only ``N_t`` and ``Q`` are divided by fixed tangent
    scales.  The solver Jacobian is the derivative of ``R / D_R`` with
    respect to those affine coordinates.  Therefore

    ``J_phys = diag(D_R) J_direct diag(1 / A)``

    where ``A = dx_phys / dx_direct``.  This is an exact change of variables,
    not a condition-number approximation.  The central probes used for the
    final Newton Jacobian are consequently the same declared physical probes
    after multiplication by ``A``; a second dense ``2n`` residual audit would
    be redundant.
    """

    layout = departure.layout_for_departure(before)
    matrix = (
        np.asarray(direct_jacobian.toarray(), dtype=float)
        if scipy_sparse.issparse(direct_jacobian)
        else np.asarray(direct_jacobian, dtype=float)
    )
    coordinate_deltas = np.asarray(coordinate_perturbations, dtype=float)
    expected_shape = (layout.residual_count, layout.unknown_count)
    if matrix.shape != expected_shape or coordinate_deltas.shape != (
        layout.unknown_count,
    ):
        raise conditioning.JacobianCertificateError(
            "direct departure Jacobian certificate lost nonlinear rank"
        )
    labels, variable_scales = _physical_variable_scale_contract(
        before,
        chart,
        scale_reference,
    )
    affine_factors = np.asarray(
        (
            *((1.0,) * layout.wet_piece_count),
            *((1.0,) * layout.wet_piece_count),
            1.0,
            *((1.0,) * (layout.dry_piece_count - 1)),
            1.0,
            *((1.0,) * (layout.dry_piece_count - 1)),
            *scale_reference.stefan_flux_scales_mol_m2_s,
            scale_reference.q_scale_m3_s,
            1.0,
        ),
        dtype=float,
    )
    declared_coordinate_scales = np.asarray(variable_scales) / affine_factors
    solver_coordinate_scales = np.asarray(
        _solver_coordinate_characteristic_scales(
            before,
            chart,
            scale_reference,
        )
    )
    if not np.array_equal(declared_coordinate_scales, solver_coordinate_scales):
        raise conditioning.JacobianCertificateError(
            "direct chart and physical D_X contracts are not the same affine map"
        )
    residual_scale_vector = np.asarray(residual_scales.vector, dtype=float)
    physical_jacobian = (
        matrix * residual_scale_vector[:, None] / affine_factors[None, :]
    )
    physical_perturbations = coordinate_deltas * affine_factors
    if not np.all(np.isfinite(physical_perturbations)) or np.any(
        physical_perturbations <= 0.0
    ):
        raise conditioning.JacobianCertificateError(
            "direct departure Jacobian perturbations are not physical"
        )
    certificate = conditioning.certify_unit_invariant_jacobian(
        physical_jacobian,
        residual_scales.vector,
        variable_scales,
        rank_relative_tolerance=_JACOBIAN_RANK_RELATIVE_TOLERANCE,
    )
    return (
        certificate,
        labels,
        variable_scales,
        tuple(float(value) for value in physical_perturbations),
    )


def _residual_scales(
    before: departure.FaceDepartureState,
    initial: departure.FaceDepartureAssembly,
    dt_s: float,
) -> ci.CutResidualScales:
    layout = initial.layout
    floor = 1.0e-300
    covariant = departure.datum_covariant_residual_vector(initial.residuals)
    wet_water = []
    wet_energy = []
    for position, cell_index in enumerate(layout.wet_cell_indices):
        old_cell = initial.old_wet_cells[position]
        new_cell = initial.wet_cells[position]
        old_volume = before.geometry.cells[cell_index].wet_volume_m3
        new_volume = initial.candidate_geometry.cells[cell_index].wet_volume_m3
        wet_water.append(
            max(
                old_volume * old_cell.retained_water_concentration_mol_m3 / dt_s,
                new_volume * new_cell.retained_water_concentration_mol_m3 / dt_s,
                abs(initial.residuals.wet_water_mol_s[position]),
                floor,
            )
        )
        wet_energy.append(
            max(
                _wet_capacity_energy(old_volume, old_cell, before) / dt_s,
                _wet_capacity_energy(new_volume, new_cell, before) / dt_s,
                abs(initial.residuals.wet_energy_w[position]),
                floor,
            )
        )
    dry_water = []
    dry_hexane = []
    dry_energy = []
    dry_capacity_audits: list[ci.DryEnergyCapacityStencilAudit] = []
    new_dry_capacity_energies: list[float] = []
    for position, cell_index in enumerate(layout.dry_cell_indices):
        new_cell = initial.dry_cells[position]
        new_volume = initial.candidate_geometry.cells[cell_index].dry_volume_m3
        new_capacity, new_capacity_audit = _dry_capacity_energy_with_audit(
            new_volume,
            new_cell,
            before,
            cell_index,
        )
        new_dry_capacity_energies.append(new_capacity)
        if new_capacity_audit is not None:
            dry_capacity_audits.append(new_capacity_audit)
        if position == 0:
            old_water = 0.0
            old_hexane = 0.0
            old_capacity = 0.0
        else:
            old_cell = initial.old_dry_cells[position - 1]
            old_volume = before.geometry.cells[cell_index].dry_volume_m3
            old_water = old_volume * old_cell.total_water_concentration_mol_m3
            old_hexane = old_volume * old_cell.total_hexane_concentration_mol_m3
            old_capacity, old_capacity_audit = _dry_capacity_energy_with_audit(
                old_volume,
                old_cell,
                before,
                cell_index,
            )
            if old_capacity_audit is not None:
                dry_capacity_audits.append(old_capacity_audit)
        dry_water.append(
            max(
                old_water / dt_s,
                new_volume * new_cell.total_water_concentration_mol_m3 / dt_s,
                abs(initial.residuals.dry_water_mol_s[position]),
                floor,
            )
        )
        dry_hexane.append(
            max(
                old_hexane / dt_s,
                new_volume * new_cell.total_hexane_concentration_mol_m3 / dt_s,
                abs(initial.residuals.dry_hexane_mol_s[position]),
                floor,
            )
        )
        dry_energy.append(
            max(
                old_capacity / dt_s,
                new_capacity / dt_s,
                abs(covariant[
                    2 * layout.wet_piece_count
                    + 2 * layout.dry_piece_count
                    + position
                ]),
                floor,
            )
        )
    cut_wet = initial.wet_cells[-1]
    cut_dry = initial.dry_cells[0]
    cut_index = layout.cut_cell_index
    new_wet_volume = initial.candidate_geometry.cells[cut_index].wet_volume_m3
    new_dry_volume = initial.candidate_geometry.cells[cut_index].dry_volume_m3
    return ci.CutResidualScales(
        wet_water_mol_s=tuple(wet_water),
        wet_energy_w=tuple(wet_energy),
        dry_water_mol_s=tuple(dry_water),
        dry_hexane_mol_s=tuple(dry_hexane),
        dry_energy_w=tuple(dry_energy),
        rh_water_mol_s=max(
            new_wet_volume * cut_wet.retained_water_concentration_mol_m3 / dt_s,
            new_dry_volume * cut_dry.total_water_concentration_mol_m3 / dt_s,
            abs(initial.residuals.rh_water_mol_s),
            floor,
        ),
        rh_hexane_mol_s=max(
            new_wet_volume
            * before.config.wet.wet.rho_dm_p
            * before.historical_hexane_loadings[cut_index]
            / hx.M
            / dt_s,
            new_dry_volume * cut_dry.total_hexane_concentration_mol_m3 / dt_s,
            abs(initial.residuals.rh_hexane_mol_s),
            floor,
        ),
        rh_energy_w=max(
            _wet_capacity_energy(new_wet_volume, cut_wet, before) / dt_s,
            new_dry_capacity_energies[0] / dt_s,
            abs(covariant[-1]),
            floor,
        ),
        dry_energy_capacity_stencil_audits=tuple(dry_capacity_audits),
    )


def _wet_capacity_energy(
    volume_m3: float,
    cell,
    state: departure.FaceDepartureState,
) -> float:
    if volume_m3 == 0.0:
        return 0.0
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
        raise ValueError("wet departure capacity scale must be positive")
    return result


def _dry_capacity_energy(
    volume_m3: float,
    cell: cp.EquilibriumPoreState,
    state: departure.FaceDepartureState,
    cell_index: int,
) -> float:
    """Compatibility scalar view of the audited departure capacity scale."""

    return _dry_capacity_energy_with_audit(
        volume_m3,
        cell,
        state,
        cell_index,
    )[0]


def _dry_capacity_energy_with_audit(
    volume_m3: float,
    cell: cp.EquilibriumPoreState,
    state: departure.FaceDepartureState,
    cell_index: int,
) -> tuple[float, ci.DryEnergyCapacityStencilAudit | None]:
    if volume_m3 == 0.0:
        return 0.0, None
    capacity, audit = ci._dry_energy_capacity_with_audit(  # noqa: SLF001
        cell.temperature_k,
        cell.y_hexane,
        replace(state.config.dry.pore, w_o=state.oil_fraction_labels[cell_index]),
        state.config,
    )
    result = volume_m3 * capacity * max(abs(cell.temperature_k), 1.0)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("dry departure capacity scale must be positive")
    return result, audit


def _departure_dry_energy_capacity(
    temperature_k: float,
    y_hexane: float,
    pore: cp.CoupledPoreParams,
    config: cut.CutTransportConfig,
) -> float:
    """Compatibility view delegated to the shared Amendment-19 helper."""

    return ci._dry_energy_capacity(  # noqa: SLF001
        temperature_k,
        y_hexane,
        pore,
        config,
    )


def _face_inventory(
    state: departure.FaceDepartureState,
) -> ci.MaterialInventorySnapshot:
    grid = state.geometry.master_grid
    k = state.departure_face_index
    wet = departure._evaluate_old_wet_cells(state)  # noqa: SLF001
    dry = departure._evaluate_old_dry_cells(state)  # noqa: SLF001
    water = [0.0] * grid.n
    hexane = [0.0] * grid.n
    energy = [0.0] * grid.n
    for index, cell in zip(range(k), wet):
        volume = grid.volumes[index]
        water[index] = volume * cell.retained_water_concentration_mol_m3
        hexane[index] = (
            volume
            * state.config.wet.wet.rho_dm_p
            * state.historical_hexane_loadings[index]
            / hx.M
        )
        energy[index] = volume * cell.energy_density_j_m3
    for index, cell in zip(range(k, grid.n), dry):
        volume = grid.volumes[index]
        water[index] = volume * cell.total_water_concentration_mol_m3
        hexane[index] = volume * cell.total_hexane_concentration_mol_m3
        energy[index] = volume * cell.energy_density_j_m3
    return ci.MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def _face_capacity_scale(state: departure.FaceDepartureState) -> float:
    grid = state.geometry.master_grid
    k = state.departure_face_index
    wet = departure._evaluate_old_wet_cells(state)  # noqa: SLF001
    dry = departure._evaluate_old_dry_cells(state)  # noqa: SLF001
    values = [
        _wet_capacity_energy(grid.volumes[index], cell, state)
        for index, cell in zip(range(k), wet)
    ]
    values.extend(
        _dry_capacity_energy(grid.volumes[index], cell, state, index)
        for index, cell in zip(range(k, grid.n), dry)
    )
    result = math.fsum(values)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("face-state capacity scale must be positive")
    return result


def _material_changes(
    assembly: departure.FaceDepartureAssembly,
) -> ci.MaterialInventorySnapshot:
    n = assembly.candidate_geometry.master_grid.n
    water = [0.0] * n
    hexane = [0.0] * n
    energy = [0.0] * n
    for position, index in enumerate(assembly.layout.wet_cell_indices):
        water[index] = math.fsum(
            (water[index], assembly.inventory_changes.wet_water_mol[position])
        )
        hexane[index] = math.fsum(
            (hexane[index], assembly.inventory_changes.wet_hexane_mol[position])
        )
        energy[index] = math.fsum(
            (energy[index], assembly.inventory_changes.wet_energy_j[position])
        )
    for position, index in enumerate(assembly.layout.dry_cell_indices):
        water[index] = math.fsum(
            (water[index], assembly.inventory_changes.dry_water_mol[position])
        )
        hexane[index] = math.fsum(
            (hexane[index], assembly.inventory_changes.dry_hexane_mol[position])
        )
        energy[index] = math.fsum(
            (energy[index], assembly.inventory_changes.dry_energy_j[position])
        )
    return ci.MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def _departure_new_inventory(
    assembly: departure.FaceDepartureAssembly,
) -> ci.MaterialInventorySnapshot:
    n = assembly.candidate_geometry.master_grid.n
    water = [0.0] * n
    hexane = [0.0] * n
    energy = [0.0] * n
    for position, index in enumerate(assembly.layout.wet_cell_indices):
        volume = assembly.candidate_geometry.cells[index].wet_volume_m3
        cell = assembly.wet_cells[position]
        water[index] += volume * cell.retained_water_concentration_mol_m3
        hexane[index] += (
            volume
            * assembly.before.config.wet.wet.rho_dm_p
            * assembly.before.historical_hexane_loadings[index]
            / hx.M
        )
        energy[index] += volume * cell.energy_density_j_m3
    for position, index in enumerate(assembly.layout.dry_cell_indices):
        volume = assembly.candidate_geometry.cells[index].dry_volume_m3
        cell = assembly.dry_cells[position]
        water[index] += volume * cell.total_water_concentration_mol_m3
        hexane[index] += volume * cell.total_hexane_concentration_mol_m3
        energy[index] += volume * cell.energy_density_j_m3
    return ci.MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def _departure_new_capacity_scale(
    assembly: departure.FaceDepartureAssembly,
) -> float:
    values = []
    for position, index in enumerate(assembly.layout.wet_cell_indices):
        values.append(
            _wet_capacity_energy(
                assembly.candidate_geometry.cells[index].wet_volume_m3,
                assembly.wet_cells[position],
                assembly.before,
            )
        )
    for position, index in enumerate(assembly.layout.dry_cell_indices):
        values.append(
            _dry_capacity_energy(
                assembly.candidate_geometry.cells[index].dry_volume_m3,
                assembly.dry_cells[position],
                assembly.before,
                index,
            )
        )
    result = math.fsum(values)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("accepted departure capacity scale must be positive")
    return result


def _accept_departure(
    before: face_solver.FaceArrivalCommittedState,
    controls: ci.CutSolverControls,
    boundary: ct.PoreBoundary,
    seed: FaceDepartureSeed,
    assembly: departure.FaceDepartureAssembly,
    scales: ci.CutResidualScales,
    scaled_residuals: tuple[float, ...],
    evaluations: int,
    rejected_trials: int,
    nonlinear_message: str,
    condition_proxy: float,
    raw_solver_condition_proxy: float,
    physical_variable_scale_labels: tuple[str, ...],
    physical_variable_scales: tuple[float, ...],
    physical_jacobian_perturbations: tuple[float, ...],
    jacobian_certificate: conditioning.UnitInvariantJacobianCertificate,
    chart: _CoordinateChart,
    structural: cut_sparsity.CutJacobianSparsity | None,
    nonlinear_method: str,
    auxiliary_stages: tuple[FaceDepartureAuxiliaryStage, ...],
    newborn_trace_gradients: departure.NewbornDryTraceGradients,
    continuation_fallback_used: bool,
    direct_requested_attempt_evaluations: int,
    direct_requested_failure_message: str | None,
    direct_requested_failure_maximum_scaled_residual: float | None,
    cauchy_globalization_accepted_steps: int,
    cauchy_line_search_trials: int,
    cauchy_residual_evaluations: int,
    dogleg_globalization_accepted_steps: int,
    dogleg_jacobian_linearizations: int,
    dogleg_trial_attempts: int,
    dogleg_residual_evaluations: int,
) -> FaceDepartureIntegratorStep:
    candidate = assembly.candidate
    if not math.isfinite(before.time_s + assembly.dt_s):
        raise ValueError("accepted departure time overflowed")
    after_transport = cut.CutTransportState(
        geometry=assembly.candidate_geometry,
        config=before.config,
        time_s=before.time_s + assembly.dt_s,
        wet_temperatures_k=candidate.wet_temperatures_k,
        wet_retained_water_loadings=candidate.wet_retained_water_loadings,
        dry_temperatures_k=candidate.dry_temperatures_k,
        dry_y_hexane=candidate.dry_y_hexane,
        historical_hexane_loadings=before.historical_hexane_loadings,
        oil_fraction_labels=before.oil_fraction_labels,
    )
    old_inventory = _face_inventory(assembly.before)
    if old_inventory != before.current_inventory:
        raise RuntimeError(
            "committed face inventory differs from the departure adapter state"
        )
    new_inventory = _departure_new_inventory(assembly)
    changes = _material_changes(assembly)
    area = after_transport.geometry.master_grid.areas[-1]
    surface = assembly.dry_face_fluxes[-1]
    water_out = (
        assembly.dt_s
        * area
        * surface.component.conserved_water_flux_mol_m2_s
    )
    hexane_out = (
        assembly.dt_s
        * area
        * surface.component.conserved_hexane_flux_mol_m2_s
    )
    energy_out = assembly.dt_s * area * surface.energy.total_energy_flux_w_m2
    water_step = math.fsum((*changes.water_cell_mol, water_out))
    hexane_step = math.fsum((*changes.hexane_cell_mol, hexane_out))
    energy_step = math.fsum((*changes.energy_cell_j, energy_out))
    old_capacity = _face_capacity_scale(assembly.before)
    new_capacity = _departure_new_capacity_scale(assembly)
    step_water_scale = max(
        old_inventory.total_water_mol,
        new_inventory.total_water_mol,
        abs(water_out),
        1.0e-300,
    )
    step_hexane_scale = max(
        old_inventory.total_hexane_mol,
        new_inventory.total_hexane_mol,
        abs(hexane_out),
        1.0e-300,
    )
    step_energy_scale = max(old_capacity, new_capacity, 1.0e-300)

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
    corrected_boundary_energy = math.fsum(
        (boundary_energy, boundary_energy_comp)
    )
    cumulative_water_residual = math.fsum(
        (*cumulative_water, *water_comp, corrected_boundary_water)
    )
    cumulative_hexane_residual = math.fsum(
        (*cumulative_hexane, *hexane_comp, corrected_boundary_hexane)
    )
    cumulative_energy_residual = math.fsum(
        (*cumulative_energy, *energy_comp, corrected_boundary_energy)
    )
    cumulative_water_scale = max(
        before.reference_inventory.total_water_mol,
        new_inventory.total_water_mol,
        abs(corrected_boundary_water),
        before.cumulative_absolute_water_transfer_mol + abs(water_out),
        1.0e-300,
    )
    cumulative_hexane_scale = max(
        before.reference_inventory.total_hexane_mol,
        new_inventory.total_hexane_mol,
        abs(corrected_boundary_hexane),
        before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out),
        1.0e-300,
    )
    cumulative_energy_scale = max(
        before.reference_capacity_energy_scale_j,
        new_capacity,
        before.cumulative_absolute_energy_transfer_j + abs(energy_out),
        1.0e-300,
    )
    after = ci.CutIntegratorState(
        transport=after_transport,
        controls=controls,
        reference_inventory=before.reference_inventory,
        reference_capacity_energy_scale_j=before.reference_capacity_energy_scale_j,
        last_interface_temperature_k=candidate.interface_temperature_k,
        last_total_stefan_fluxes_mol_m2_s=(
            candidate.dry_total_stefan_fluxes_mol_m2_s
        ),
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
    fractions = _bounded_fractions(candidate, assembly.before, chart)
    front_distances = (
        candidate.front_z - chart.front_z[0],
        chart.front_z[1] - candidate.front_z,
    )
    ledger = FaceDepartureIntegratorLedger(
        dt_s=assembly.dt_s,
        nonlinear_method=nonlinear_method,
        nonlinear_evaluations=evaluations,
        rejected_trial_evaluations=rejected_trials,
        nonlinear_message=nonlinear_message,
        residual_scales=scales,
        scaled_residuals=scaled_residuals,
        maximum_scaled_residual=max(abs(value) for value in scaled_residuals),
        condition_proxy=condition_proxy,
        raw_solver_condition_proxy=raw_solver_condition_proxy,
        physical_variable_scale_provenance=(
            "unique exact-face tangent; fixed 1 s existing-rate widths; "
            "|g_T|/|g_y|/N_t/Q documented floors; T_Gamma=5 K"
        ),
        physical_variable_scale_labels=physical_variable_scale_labels,
        physical_variable_scales=physical_variable_scales,
        physical_jacobian_perturbations=physical_jacobian_perturbations,
        jacobian_certificate=jacobian_certificate,
        minimum_fractional_distance_to_bound=min(
            min(value, 1.0 - value) for value in fractions
        ),
        minimum_front_z_distance_to_chart_boundary=min(front_distances),
        used_sparse_jacobian=structural is not None,
        scalability_class=(
            "bordered-banded colored central differences with sparse LU; "
            "O(N) structural storage, dense raw-condition/J* SVD diagnostics"
            if structural is not None
            else "dense finite-difference safeguarded direct-Newton correctness "
            "kernel"
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
        jacobian_structural_storage_bytes=(
            8 * assembly.layout.unknown_count**2
            if structural is None
            else (
                8 * structural.nnz
                + structural.jac_sparsity.indices.nbytes
                + structural.jac_sparsity.indptr.nbytes
            )
        ),
        jacobian_dense_equivalent_bytes=(
            8 * assembly.layout.unknown_count**2
        ),
        continuation_fallback_used=continuation_fallback_used,
        direct_requested_attempt_evaluations=(
            direct_requested_attempt_evaluations
        ),
        direct_requested_failure_message=direct_requested_failure_message,
        direct_requested_failure_maximum_scaled_residual=(
            direct_requested_failure_maximum_scaled_residual
        ),
        auxiliary_stages=auxiliary_stages,
        cauchy_globalization_used=(
            cauchy_globalization_accepted_steps > 0
        ),
        cauchy_globalization_accepted_steps=(
            cauchy_globalization_accepted_steps
        ),
        cauchy_line_search_trials=cauchy_line_search_trials,
        cauchy_residual_evaluations=cauchy_residual_evaluations,
        dogleg_globalization_used=(
            dogleg_globalization_accepted_steps > 0
        ),
        dogleg_globalization_accepted_steps=(
            dogleg_globalization_accepted_steps
        ),
        dogleg_jacobian_linearizations=dogleg_jacobian_linearizations,
        dogleg_trial_attempts=dogleg_trial_attempts,
        dogleg_residual_evaluations=dogleg_residual_evaluations,
        newborn_trace_gradients=newborn_trace_gradients,
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
        source_accepted_steps=before.accepted_steps,
        source_cumulative_nonlinear_evaluations=(
            before.cumulative_nonlinear_evaluations
        ),
        accepted=True,
    )
    return FaceDepartureIntegratorStep(
        before,
        after,
        assembly,
        boundary,
        controls,
        seed,
        ledger,
    )


__all__ = [
    "FaceDepartureAuxiliaryStage",
    "FaceDepartureIntegratorLedger",
    "FaceDepartureIntegratorStep",
    "FaceDepartureIntegratorStepError",
    "FaceDepartureSparsityAudit",
    "FaceDepartureSeed",
    "audit_face_departure_sparsity",
    "solve_face_departure",
]
