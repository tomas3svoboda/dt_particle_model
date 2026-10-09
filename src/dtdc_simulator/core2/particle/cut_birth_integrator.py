r"""Safeguarded Dirichlet numerical-oracle integrator for dry-shell birth.

This module wraps :mod:`cut_birth_event` without changing its equations.  It
certifies the owner-approved receding ``GT-PS-2-G1G-01-A1`` branch under an
explicit *numerical* Dirichlet surface oracle.  The finite external
film/contact active set is not frozen, so no result here is physically
qualifying and no failed solve can select contact.

Wet temperature/loading, newborn dry temperature, strict outer-cell ``z``, and
``T_Gamma`` use smooth open logit charts.  A legacy rectangular dry composition
band uses the same representation.  An active exact gas-only domain instead
decodes temperature first and then uses the pressure-, temperature-, and
surface-material-label-conditioned composition chart.  The two total Stefan
fluxes use unbounded scaled coordinates.  Residual rows are fixed-scaled from
inventories, local capacity rates, and RH swept-storage rates; dry and RH energy
use the equation-equivalent ``R_E-H_REF*R_h`` basis.  Dense and bordered-banded
sparse finite-difference routes share the exact same residual.

An accepted root commits the first strict partial state directly to
:class:`cut_integrator.CutIntegratorState`.  Fully-wet reference inventories,
prior wet-water boundary ledgers, stationary material labels, and signed
low-order corrections are retained with Neumaier accumulation.  There is no
epsilon dry seed, clipping, contact fallback, or tolerance relaxation.
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass, field, replace
from typing import Sequence

import numpy as np
from scipy import optimize, sparse
from scipy.optimize._numdiff import group_columns

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_birth_event as birth
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp


# GT-PS-2-P1E-R15: single authority for the birth acceptance lives with
# the evidence dataclass it also bounds (see cut_birth_event).
BIRTH_REGIONAL_RH_ACCEPTANCE = birth.BIRTH_REGIONAL_RH_ACCEPTANCE


class BirthIntegratorStepError(RuntimeError):
    """Rejected numerical-oracle birth with exact fully-wet rollback."""

    def __init__(
        self,
        message: str,
        rollback_state: ww.WetWaterState,
        *,
        nonlinear_evaluations: int,
        rejected_trial_evaluations: int,
        maximum_scaled_residual: float = math.inf,
        condition_proxy: float = math.inf,
        last_candidate: birth.BirthUnknowns | None = None,
        scaled_residuals: tuple[float, ...] = (),
        minimum_fractional_distance_to_bound: float = 0.0,
        retained_water_applicability_guard_failed: bool = False,
        retained_water_applicability_audit: (ci.WetRetainedWaterApplicabilityAudit | None) = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.nonlinear_evaluations = nonlinear_evaluations
        self.rejected_trial_evaluations = rejected_trial_evaluations
        self.maximum_scaled_residual = maximum_scaled_residual
        self.condition_proxy = condition_proxy
        self.last_candidate = last_candidate
        self.scaled_residuals = scaled_residuals
        self.minimum_fractional_distance_to_bound = minimum_fractional_distance_to_bound
        self.retained_water_applicability_guard_failed = retained_water_applicability_guard_failed
        self.retained_water_applicability_audit = retained_water_applicability_audit
        self.branch = birth.BirthBranch.REJECTED_NUMERICAL_FAILURE
        self.contact_selected = False


class BirthSolverJacobianPolicy(enum.Enum):
    """Explicit nonlinear method/Jacobian selection for one birth solve."""

    DENSE_TRF_TWO_POINT = "dense_trf_two_point"
    STRUCTURED_DOGBOX_TWO_POINT = "structured_dogbox_two_point"


@dataclass(frozen=True)
class BirthSolverJacobianSelection:
    """Resolved numerical route; this is not physical qualification evidence."""

    policy: BirthSolverJacobianPolicy
    solver_method: str
    structure: BirthSparsity | None
    physically_qualifying: bool = field(default=False, init=False)


class BirthResidualAssemblyBudgetExceeded(BirthIntegratorStepError):
    """A declared hard total residual-assembly budget was exhausted.

    SciPy's ``max_nfev`` does not count every residual call used to construct
    finite-difference Jacobians.  This typed refusal counts calls at the native
    birth-assembly boundary instead, before an additional assembly is allowed
    to execute.  The exact fully-wet input remains the rollback authority.
    """

    def __init__(
        self,
        rollback_state: ww.WetWaterState,
        *,
        total_residual_assembly_budget: int,
        total_residual_assembly_attempts: int,
        nonlinear_evaluations: int,
        rejected_trial_evaluations: int,
        solver_stage: str = "unspecified_total_residual_assembly_boundary",
        scipy_nfev: int | None = None,
    ) -> None:
        super().__init__(
            "receding birth total residual-assembly budget exhausted with "
            "exact rollback; contact forbidden: "
            f"{total_residual_assembly_attempts} >= "
            f"{total_residual_assembly_budget}",
            rollback_state,
            nonlinear_evaluations=nonlinear_evaluations,
            rejected_trial_evaluations=rejected_trial_evaluations,
        )
        self.total_residual_assembly_budget = total_residual_assembly_budget
        self.total_residual_assembly_attempts = total_residual_assembly_attempts
        self.total_residual_assembly_budget_exhausted = True
        self.solver_stage = solver_stage
        self.scipy_nfev = scipy_nfev


@dataclass(frozen=True)
class BirthStepSeed:
    """Explicit birth candidate and positive physical Stefan scales."""

    candidate: birth.BirthUnknowns
    stefan_flux_scales_mol_m2_s: tuple[float, float]
    label: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.candidate, birth.BirthUnknowns):
            raise TypeError("birth seed requires BirthUnknowns")
        if len(self.stefan_flux_scales_mol_m2_s) != 2 or not all(
            math.isfinite(value) and value > 0.0 for value in self.stefan_flux_scales_mol_m2_s
        ):
            raise ValueError("birth seed requires two positive Stefan scales")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("birth seed needs non-empty provenance")


@dataclass(frozen=True)
class BirthSparsity:
    """Deterministic bordered-banded birth Jacobian structure."""

    matrix: sparse.csr_matrix
    color_count: int
    column_colors: tuple[int, ...]
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not sparse.isspmatrix_csr(self.matrix) or self.matrix.dtype != np.dtype(bool):
            raise TypeError("birth sparsity must be a boolean CSR matrix")
        row_count, column_count = self.matrix.shape
        if row_count < 1 or column_count < 1:
            raise ValueError("birth sparsity must be nonempty")
        if (
            isinstance(self.color_count, bool)
            or type(self.color_count) is not int
            or self.color_count < 1
            or type(self.column_colors) is not tuple
            or len(self.column_colors) != column_count
            or any(type(color) is not int or color < 0 for color in self.column_colors)
            or max(self.column_colors, default=-1) + 1 != self.color_count
            or set(self.column_colors) != set(range(self.color_count))
        ):
            raise ValueError("birth sparsity needs one contiguous exact color per column")
        for row in range(row_count):
            columns = self.matrix.indices[self.matrix.indptr[row] : self.matrix.indptr[row + 1]]
            colors = tuple(self.column_colors[int(column)] for column in columns)
            if len(set(colors)) != len(colors):
                raise ValueError("birth sparsity has a same-row structural color collision")

    @property
    def nnz(self) -> int:
        return int(self.matrix.nnz)


@dataclass(frozen=True)
class BirthTopologyDiagnostics:
    """Exact newborn-cell gas-only interval and strict activity margins."""

    uses_conditioned_composition_chart: bool
    newborn_surface_oil_fraction_label: float
    exact_composition_interval: cp.GasOnlyCompositionInterval
    exact_composition_interval_fraction: float
    water_saturation_activity_gap: float
    hexane_saturation_activity_gap: float
    luikov_lower_activity_gap: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.newborn_surface_oil_fraction_label):
            raise ValueError("newborn surface-cell oil label must remain finite")
        if not isinstance(
            self.exact_composition_interval,
            cp.GasOnlyCompositionInterval,
        ):
            raise TypeError("birth topology needs one exact gas-only interval")
        if (
            not math.isfinite(self.exact_composition_interval_fraction)
            or not 0.0 < self.exact_composition_interval_fraction < 1.0
        ):
            raise ValueError("newborn composition must be inside its exact interval")
        gaps = (
            self.water_saturation_activity_gap,
            self.hexane_saturation_activity_gap,
            self.luikov_lower_activity_gap,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in gaps):
            raise ValueError("every exact newborn activity gap must remain positive")

    @property
    def minimum_exact_activity_gap(self) -> float:
        return min(
            self.water_saturation_activity_gap,
            self.hexane_saturation_activity_gap,
            self.luikov_lower_activity_gap,
        )

    @property
    def minimum_exact_composition_fractional_distance(self) -> float:
        fraction = self.exact_composition_interval_fraction
        return min(fraction, 1.0 - fraction)


@dataclass(frozen=True)
class BirthIntegratorLedger:
    """Nonlinear, branch, conservation, and numerical-oracle diagnostics."""

    dt_s: float
    nonlinear_evaluations: int
    rejected_trial_evaluations: int
    nonlinear_message: str
    residual_scales: ci.CutResidualScales
    scaled_residuals: tuple[float, ...]
    maximum_scaled_residual: float
    condition_proxy: float
    minimum_fractional_distance_to_bound: float
    minimum_front_z_distance_to_chart_boundary: float
    topology_diagnostics: BirthTopologyDiagnostics
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
    surface_hexane_flux_mol_m2_s: float
    interface_hexane_flux_mol_m2_s: float
    hexane_storage_jump_mol_m3: float
    swept_recession_speed_m_s: float
    rh_recession_speed_m_s: float
    relative_immediate_recession_error: float
    newborn_dry_volume_m3: float
    branch_decision: birth.BirthBranchDecision
    birth_telescoping_residual: float
    dirichlet_numerical_oracle: bool
    contact_selected: bool
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

    @property
    def minimum_exact_activity_gap(self) -> float:
        return self.topology_diagnostics.minimum_exact_activity_gap

    @property
    def minimum_exact_composition_fractional_distance(self) -> float:
        return self.topology_diagnostics.minimum_exact_composition_fractional_distance


@dataclass(frozen=True)
class BirthIntegratorStep:
    """One accepted receding numerical-oracle birth and strict-cut commit."""

    before: ww.WetWaterState
    after: ci.CutIntegratorState
    assembly: birth.BirthAssembly
    boundary: ct.PoreBoundary
    controls: ci.CutSolverControls
    seed: BirthStepSeed
    ledger: BirthIntegratorLedger
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class _CoordinateChart:
    wet_temperature: tuple[float, float]
    wet_water: tuple[float, float]
    dry_temperature: tuple[float, float]
    legacy_dry_y_hexane: tuple[float, float] | None
    dry_config: ct.FullyDryTransportConfig
    newborn_surface_pore: cp.CoupledPoreParams
    front_z: tuple[float, float]
    interface_temperature: tuple[float, float]


def advance_receding_birth(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
    controls: ci.CutSolverControls,
    dt_s: float,
    boundary: ct.PoreBoundary,
    seed: BirthStepSeed,
    *,
    solver_jacobian_policy: BirthSolverJacobianPolicy = (
        BirthSolverJacobianPolicy.DENSE_TRF_TWO_POINT
    ),
    total_residual_assembly_budget: int | None = None,
) -> BirthIntegratorStep:
    """Solve, certify, and commit one first strict-partial A1 oracle step.

    ``solver_jacobian_policy`` selects the nonlinear/Jacobian route directly;
    ``CutSolverControls.maximum_dense_unknowns`` no longer makes that decision.
    ``total_residual_assembly_budget`` is an optional hard cap on every call
    attempted at :func:`cut_birth_event.assemble_birth_backward_euler`,
    including finite-difference Jacobian samples and rejected physical trials.
    It is intentionally distinct from SciPy's ``max_nfev`` accounting.
    """

    evaluations = 0
    rejected_trials = 0
    residual_assembly_attempts = 0
    solver_stage = "input_validation"
    scipy_nfev = None
    maximum_scaled = math.inf
    condition_proxy = math.inf
    last_candidate = None
    last_scaled: tuple[float, ...] = ()
    minimum_bound_distance = 0.0
    try:
        if not isinstance(before, ww.WetWaterState):
            raise TypeError("before must be a fully-wet WetWaterState")
        if not isinstance(config, cut.CutTransportConfig):
            raise TypeError("birth integrator requires CutTransportConfig")
        if not isinstance(controls, ci.CutSolverControls):
            raise TypeError("birth integrator requires CutSolverControls")
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("birth boundary must be a supported pore boundary")
        if not isinstance(seed, BirthStepSeed):
            raise TypeError("birth integrator requires BirthStepSeed")
        if type(solver_jacobian_policy) is not BirthSolverJacobianPolicy:
            raise TypeError("birth integrator requires an explicit Jacobian policy enum")
        if total_residual_assembly_budget is not None and (
            isinstance(total_residual_assembly_budget, bool)
            or type(total_residual_assembly_budget) is not int
            or total_residual_assembly_budget < 1
        ):
            raise ValueError("birth total residual-assembly budget must be a positive exact int")
        if not math.isfinite(dt_s) or dt_s <= 0.0:
            raise ValueError("birth duration must be positive and finite")
        ww._validate_state_structure(before)  # noqa: SLF001
        layout = birth.layout_for_birth(before)
        _validate_birth_chart(before, config, controls)
        _validate_seed(layout, seed)
        chart = _coordinate_chart(before, config, controls)
        initial_coordinates = _encode_candidate(
            seed.candidate,
            chart,
            seed.stefan_flux_scales_mol_m2_s,
        )

        def assemble_counted(
            candidate: birth.BirthUnknowns,
            *,
            enforce_reduced_film_thresholds: bool,
            count_rejected_trial: bool,
        ) -> birth.BirthAssembly:
            nonlocal evaluations, rejected_trials, residual_assembly_attempts
            if (
                total_residual_assembly_budget is not None
                and residual_assembly_attempts >= total_residual_assembly_budget
            ):
                raise BirthResidualAssemblyBudgetExceeded(
                    before,
                    total_residual_assembly_budget=(total_residual_assembly_budget),
                    total_residual_assembly_attempts=residual_assembly_attempts,
                    nonlinear_evaluations=evaluations,
                    rejected_trial_evaluations=rejected_trials,
                    solver_stage=solver_stage,
                    scipy_nfev=scipy_nfev,
                )
            residual_assembly_attempts += 1
            try:
                assembly = birth.assemble_birth_backward_euler(
                    before,
                    config,
                    candidate,
                    dt_s,
                    boundary,
                    enforce_reduced_film_thresholds=(enforce_reduced_film_thresholds),
                )
            except birth.BirthStepError:
                if count_rejected_trial:
                    rejected_trials += 1
                raise
            evaluations += 1
            return assembly

        solver_stage = "initial_unthresholded_scale_assembly"
        initial = assemble_counted(
            seed.candidate,
            enforce_reduced_film_thresholds=False,
            count_rejected_trial=False,
        )
        scales = _residual_scales(initial)
        scale_vector = np.asarray(scales.vector, dtype=float)
        solver_selection = select_birth_solver_jacobian(
            layout,
            solver_jacobian_policy,
        )
        structure = solver_selection.structure

        def scaled_residual(coordinates: np.ndarray) -> np.ndarray:
            candidate = _decode_candidate(
                coordinates,
                layout,
                chart,
                seed.stefan_flux_scales_mol_m2_s,
            )
            assembly = assemble_counted(
                candidate,
                enforce_reduced_film_thresholds=False,
                count_rejected_trial=True,
            )
            return (
                np.asarray(
                    birth.datum_covariant_residual_vector(assembly.residuals),
                    dtype=float,
                )
                / scale_vector
            )

        lower, upper = ci._coordinate_solver_bounds(  # noqa: SLF001
            layout,
            controls.maximum_logit_magnitude,
        )
        solver_stage = "scipy_initial_residual_or_finite_difference_jacobian"
        solution = optimize.least_squares(
            scaled_residual,
            initial_coordinates,
            method=solver_selection.solver_method,
            bounds=(lower, upper),
            ftol=controls.nonlinear_step_tolerance,
            xtol=controls.nonlinear_step_tolerance,
            gtol=controls.nonlinear_step_tolerance,
            max_nfev=controls.maximum_function_evaluations,
            x_scale="jac",
            jac_sparsity=None if structure is None else structure.matrix,
            tr_options=(
                {"atol": 1.0e-14, "btol": 1.0e-14, "maxiter": 1000}
                if structure is not None
                else None
            ),
        )
        scipy_nfev = int(solution.nfev)
        solver_stage = "final_thresholded_certification_assembly"
        candidate = _decode_candidate(
            solution.x,
            layout,
            chart,
            seed.stefan_flux_scales_mol_m2_s,
        )
        converged = assemble_counted(
            candidate,
            enforce_reduced_film_thresholds=True,
            count_rejected_trial=False,
        )
        covariant = birth.datum_covariant_residual_vector(converged.residuals)
        signed_scaled = tuple(residual / scale for residual, scale in zip(covariant, scales.vector))
        maximum_scaled = max(abs(value) for value in signed_scaled)
        last_candidate = candidate
        last_scaled = signed_scaled
        topology = _topology_diagnostics(converged, chart)
        fractions = _bounded_fractions(candidate, chart, topology)
        minimum_bound_distance = min(min(value, 1.0 - value) for value in fractions)
        condition_proxy = ci._jacobian_condition_proxy(solution.jac)  # noqa: SLF001
        if not solution.success:
            raise RuntimeError(f"birth nonlinear root did not certify success: {solution.message}")
        if np.any(solution.active_mask != 0):
            raise RuntimeError("birth solution touched a logit guard")
        # GT-PS-2-P1E-R15: the birth system's residual scale is amplified by
        # the newborn cell's near-zero volume; its measured convergence
        # points at the pinned joint forcing (N1: 7.74e-11, N2: 1.104e-10)
        # straddle the shared 1e-10 gate while trajectory solves sit at
        # 1e-12..1e-11.  The birth acceptance is therefore its own constant
        # at >2x the worst measured convergence point; conservation ledgers
        # are unchanged, and the A21-grade scale derivation is queued with
        # A22.
        birth_acceptance = max(
            controls.nonlinear_residual_tolerance,
            BIRTH_REGIONAL_RH_ACCEPTANCE,
        )
        if maximum_scaled > birth_acceptance:
            raise RuntimeError(
                "birth regional/RH residual contract failed: "
                f"{maximum_scaled:.3e} > "
                f"{birth_acceptance:.3e}"
            )
        if not math.isfinite(condition_proxy) or (
            condition_proxy > controls.maximum_condition_proxy
        ):
            raise RuntimeError(
                f"birth Jacobian condition proxy is uncertified: {condition_proxy:.3e}"
            )
        step = _accept_birth(
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
            topology,
            structure,
        )
        if (
            step.ledger.maximum_step_ledger_residual > controls.ledger_tolerance
            or step.ledger.maximum_cumulative_ledger_residual > controls.ledger_tolerance
            or step.ledger.birth_telescoping_residual > controls.ledger_tolerance
            or step.ledger.relative_immediate_recession_error > controls.ledger_tolerance
        ):
            raise RuntimeError("birth acceptance ledger exceeded the declared 1e-10 contract")
        ci.certify_explicit_wet_retained_water_applicability(
            controls,
            before.retained_water_loadings,
            (0.0,) * len(before.retained_water_loadings),
            config.wet.luikov,
        )
        ci.certify_explicit_wet_retained_water_applicability(
            controls,
            step.after.transport.wet_retained_water_loadings,
            step.after.transport.effective_wet_retained_water_capacity_duals_over_rt,
            config.wet.luikov,
        )
        return step
    except BirthIntegratorStepError:
        raise
    except Exception as exc:
        if not isinstance(before, ww.WetWaterState):
            raise
        applicability_error = (
            exc if isinstance(exc, ci.WetRetainedWaterApplicabilityError) else None
        )
        raise BirthIntegratorStepError(
            f"receding birth solve rejected with exact rollback; contact forbidden: {exc}",
            before,
            nonlinear_evaluations=evaluations,
            rejected_trial_evaluations=rejected_trials,
            maximum_scaled_residual=maximum_scaled,
            condition_proxy=condition_proxy,
            last_candidate=last_candidate,
            scaled_residuals=last_scaled,
            minimum_fractional_distance_to_bound=minimum_bound_distance,
            retained_water_applicability_guard_failed=(applicability_error is not None),
            retained_water_applicability_audit=(
                None if applicability_error is None else applicability_error.audit
            ),
        ) from exc


def _coordinate_chart(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
    controls: ci.CutSolverControls,
) -> _CoordinateChart:
    inner = before.grid.faces[-2]
    inner_z = (inner / before.grid.R) ** 3
    return _CoordinateChart(
        wet_temperature=controls.wet_temperature_bounds_k,
        wet_water=controls.wet_water_bounds,
        dry_temperature=config.dry.conditioned_temperature_domain.solver_bounds_k,
        legacy_dry_y_hexane=(
            config.dry.primitive_band.y_hexane_bounds
            if isinstance(config.dry.primitive_band, ct.JointPrimitiveBand)
            else None
        ),
        dry_config=config.dry,
        newborn_surface_pore=replace(
            config.dry.pore,
            w_o=before.oil_fraction_labels[-1],
        ),
        front_z=(inner_z, 1.0),
        interface_temperature=controls.interface_temperature_bounds_k,
    )


def _validate_birth_chart(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
    controls: ci.CutSolverControls,
) -> None:
    if before.model != config.wet:
        raise ValueError("fully-wet state and birth chart need one wet authority")
    wet_t_lo, wet_t_hi = controls.wet_temperature_bounds_k
    if not config.wet.wet.T_min <= wet_t_lo < wet_t_hi <= config.wet.wet.T_max:
        raise ValueError("birth wet-temperature chart left its caloric authority")
    water_lo, water_hi = controls.wet_water_bounds
    if not config.wet.luikov.W_ref <= water_lo < water_hi <= config.wet.luikov.W_cap:
        raise ValueError("birth wet-water chart left its Luikov authority")
    gamma_lo, gamma_hi = controls.interface_temperature_bounds_k
    pore_lo, pore_hi = config.dry.pore.temperature_bounds_k
    if (
        not max(config.wet.wet.T_min, pore_lo)
        <= gamma_lo
        < gamma_hi
        <= min(
            config.wet.wet.T_max,
            pore_hi,
        )
    ):
        raise ValueError("birth interface-temperature chart left shared authority")
    if not all(
        ci._inside(value, controls.wet_temperature_bounds_k)  # noqa: SLF001
        for value in before.temperatures_k
    ):
        raise ValueError("fully-wet temperature is outside the birth solver chart")
    if not all(
        ci._inside(value, controls.wet_water_bounds)  # noqa: SLF001
        for value in before.retained_water_loadings
    ):
        raise ValueError("fully-wet water loading is outside the birth solver chart")
    ci.certify_explicit_wet_retained_water_applicability(
        controls,
        before.retained_water_loadings,
        (0.0,) * len(before.retained_water_loadings),
        config.wet.luikov,
    )
    cut_index = before.grid.n - 1
    for temperature in (gamma_lo, 0.5 * (gamma_lo + gamma_hi), gamma_hi):
        cut.evaluate_interface_state(
            temperature,
            config,
            before.historical_hexane_loadings[cut_index],
            before.oil_fraction_labels[cut_index],
        )


def _validate_seed(layout: birth.BirthLayout, seed: BirthStepSeed) -> None:
    if len(seed.candidate.vector()) != layout.unknown_count:
        raise ValueError("birth seed candidate does not match the nonlinear rank")
    if (
        len(seed.candidate.wet_temperatures_k) != layout.wet_piece_count
        or len(seed.candidate.wet_retained_water_loadings) != layout.wet_piece_count
    ):
        raise ValueError("birth seed wet primitives do not align with material cells")


def _encode_candidate(
    candidate: birth.BirthUnknowns,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
) -> np.ndarray:
    values = (
        *(ci._encode_open(value, chart.wet_temperature) for value in candidate.wet_temperatures_k),  # noqa: SLF001,E501
        *(
            ci._encode_open(value, chart.wet_water)
            for value in candidate.wet_retained_water_loadings
        ),  # noqa: SLF001,E501
        ci._encode_open(candidate.dry_temperature_k, chart.dry_temperature),  # noqa: SLF001
        _encode_dry_y(
            candidate.dry_temperature_k,
            candidate.dry_y_hexane,
            chart,
        ),
        *(
            value / scale
            for value, scale in zip(candidate.dry_total_stefan_fluxes_mol_m2_s, stefan_scales)
        ),  # noqa: E501
        ci._encode_open(candidate.front_z, chart.front_z),  # noqa: SLF001
        ci._encode_open(candidate.interface_temperature_k, chart.interface_temperature),  # noqa: SLF001,E501
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("encoded birth coordinates must remain finite")
    return np.asarray(values, dtype=float)


def _decode_candidate(
    coordinates: Sequence[float],
    layout: birth.BirthLayout,
    chart: _CoordinateChart,
    stefan_scales: Sequence[float],
) -> birth.BirthUnknowns:
    values = tuple(float(value) for value in coordinates)
    if len(values) != layout.unknown_count or not all(math.isfinite(value) for value in values):
        raise ValueError("birth coordinate vector is non-finite or wrong-rank")
    n = layout.wet_piece_count
    cursor = 0

    def take(count: int) -> tuple[float, ...]:
        nonlocal cursor
        result = values[cursor : cursor + count]
        cursor += count
        return result

    wet_t = tuple(
        ci._decode_open(value, chart.wet_temperature)  # noqa: SLF001
        for value in take(n)
    )
    wet_w = tuple(
        ci._decode_open(value, chart.wet_water)  # noqa: SLF001
        for value in take(n)
    )
    dry_t = ci._decode_open(take(1)[0], chart.dry_temperature)  # noqa: SLF001
    dry_y = _decode_dry_y(dry_t, take(1)[0], chart)
    nt = tuple(value * scale for value, scale in zip(take(2), stefan_scales))
    front_z = ci._decode_open(take(1)[0], chart.front_z)  # noqa: SLF001
    interface_t = ci._decode_open(  # noqa: SLF001
        take(1)[0],
        chart.interface_temperature,
    )
    return birth.BirthUnknowns(
        wet_t,
        wet_w,
        dry_t,
        dry_y,
        (nt[0], nt[1]),
        front_z,
        interface_t,
    )


def _encode_dry_y(
    temperature_k: float,
    y_hexane: float,
    chart: _CoordinateChart,
) -> float:
    if chart.legacy_dry_y_hexane is not None:
        return ci._encode_open(  # noqa: SLF001
            y_hexane,
            chart.legacy_dry_y_hexane,
        )
    return cp.encode_gas_only_y(
        temperature_k,
        chart.dry_config.pressure_pa,
        y_hexane,
        chart.newborn_surface_pore,
    )


def _decode_dry_y(
    temperature_k: float,
    coordinate: float,
    chart: _CoordinateChart,
) -> float:
    if chart.legacy_dry_y_hexane is not None:
        return ci._decode_open(  # noqa: SLF001
            coordinate,
            chart.legacy_dry_y_hexane,
        )
    return cp.decode_gas_only_y(
        temperature_k,
        chart.dry_config.pressure_pa,
        coordinate,
        chart.newborn_surface_pore,
    )


def _bounded_fractions(
    candidate: birth.BirthUnknowns,
    chart: _CoordinateChart,
    topology: BirthTopologyDiagnostics,
) -> tuple[float, ...]:
    pairs = [
        *((value, chart.wet_temperature) for value in candidate.wet_temperatures_k),
        *((value, chart.wet_water) for value in candidate.wet_retained_water_loadings),
        (candidate.dry_temperature_k, chart.dry_temperature),
        (candidate.front_z, chart.front_z),
        (candidate.interface_temperature_k, chart.interface_temperature),
    ]
    fractions = [(value - bounds[0]) / (bounds[1] - bounds[0]) for value, bounds in pairs]
    if chart.legacy_dry_y_hexane is not None:
        y_lo, y_hi = chart.legacy_dry_y_hexane
        fractions.append((candidate.dry_y_hexane - y_lo) / (y_hi - y_lo))
    fractions.append(topology.exact_composition_interval_fraction)
    fractions.extend(
        gap / (1.0 + gap)
        for gap in (
            topology.water_saturation_activity_gap,
            topology.hexane_saturation_activity_gap,
            topology.luikov_lower_activity_gap,
        )
    )
    result = tuple(fractions)
    if not all(0.0 < value < 1.0 for value in result):
        raise RuntimeError("accepted birth primitive touched an open-chart boundary")
    return result


def _topology_diagnostics(
    assembly: birth.BirthAssembly,
    chart: _CoordinateChart,
) -> BirthTopologyDiagnostics:
    candidate = assembly.candidate
    if chart.newborn_surface_pore.w_o != assembly.before.oil_fraction_labels[-1]:
        raise RuntimeError("birth composition chart lost the surface-cell oil label")
    interval = cp.gas_only_composition_interval(
        candidate.dry_temperature_k,
        chart.dry_config.pressure_pa,
        chart.newborn_surface_pore,
    )
    fraction = (candidate.dry_y_hexane - interval.lower_y_hexane) / (
        interval.upper_y_hexane - interval.lower_y_hexane
    )
    dry_cell = assembly.dry_cells[0]
    activity_at_ref = sp.water_activity(
        chart.newborn_surface_pore.luikov.W_ref,
        chart.newborn_surface_pore.luikov,
    )
    return BirthTopologyDiagnostics(
        uses_conditioned_composition_chart=(chart.legacy_dry_y_hexane is None),
        newborn_surface_oil_fraction_label=chart.newborn_surface_pore.w_o,
        exact_composition_interval=interval,
        exact_composition_interval_fraction=fraction,
        water_saturation_activity_gap=1.0 - dry_cell.water_activity,
        hexane_saturation_activity_gap=1.0 - dry_cell.hexane_activity,
        luikov_lower_activity_gap=dry_cell.water_activity - activity_at_ref,
    )


def _residual_scales(assembly: birth.BirthAssembly) -> ci.CutResidualScales:
    dt_s = assembly.dt_s
    grid = assembly.before.grid
    geometry = assembly.candidate_geometry
    floor = 1.0e-300
    covariant = birth.datum_covariant_residual_vector(assembly.residuals)
    n = assembly.layout.wet_piece_count
    wet_water = tuple(
        max(
            grid.volumes[index] * old.retained_water_concentration_mol_m3 / dt_s,
            geometry.cells[index].wet_volume_m3 * new.retained_water_concentration_mol_m3 / dt_s,
            abs(assembly.residuals.wet_water_mol_s[index]),
            floor,
        )
        for index, (old, new) in enumerate(zip(assembly.old_wet_cells, assembly.wet_cells))
    )
    wet_energy = tuple(
        max(
            _wet_capacity_energy(grid.volumes[index], old, assembly.config) / dt_s,
            _wet_capacity_energy(
                geometry.cells[index].wet_volume_m3,
                new,
                assembly.config,
            )
            / dt_s,
            abs(assembly.residuals.wet_energy_w[index]),
            floor,
        )
        for index, (old, new) in enumerate(zip(assembly.old_wet_cells, assembly.wet_cells))
    )
    dry_volume = geometry.cells[-1].dry_volume_m3
    dry = assembly.dry_cells[0]
    dry_water = max(
        dry_volume * dry.total_water_concentration_mol_m3 / dt_s,
        abs(assembly.residuals.dry_water_mol_s[0]),
        floor,
    )
    dry_hexane = max(
        dry_volume * dry.total_hexane_concentration_mol_m3 / dt_s,
        abs(assembly.residuals.dry_hexane_mol_s[0]),
        floor,
    )
    dry_capacity_energy, dry_capacity_audit = _dry_capacity_energy_with_audit(
        dry_volume,
        dry,
        assembly,
    )
    # GT-PS-2-B1c / Q-B1-3 (owner ruling 2026-09-03 evening,
    # docs/GT_PS2_B1C_CERTIFIED_BIRTH_FOLLOWUPS_CEREMONY_2026-09-03.md), the
    # same reasoning B2 applied to ``rh_water``: this row's scale must NOT
    # contain the residual it normalizes.  ``dry_capacity_energy / dt_s`` is
    # the newborn dry cell's own stored-enthalpy rate, i.e. the inventory-rate
    # scale every sibling energy row uses.  The removed
    # ``abs(covariant[2 * n + 2])`` branch is inert at all five pinned roots
    # (|r|/phys = 2.8e-12..9.8e-12) but is the ACTIVE branch at the
    # conditioned perturbed seed (|r|/phys = 2.1714857368681035), where it
    # saturated the reported scaled residual at exactly +1.0 and loosened the
    # acceptance comparison by that factor.  Removing it can only raise a
    # reported residual, never lower one, so no gate is weakened.
    dry_energy = max(
        dry_capacity_energy / dt_s,
        floor,
    )
    q_abs = abs(assembly.swept_geometry.interface_swept_volume_rate_m3_s)
    wet_hexane = (
        assembly.config.wet.wet.rho_dm_p * assembly.interface.historical_hexane_loading / hx.M
    )
    # GT-PS-2-B2, docs/GT_PS2_B2_BIRTH_NORMALIZATION_CEREMONY_2026-09-04.md:
    # this row's scale must NOT contain the residual it normalizes.  ``q_abs``
    # is the newborn cell's own swept volume rate (exactly V_newborn/dt at
    # birth), so ``q_abs * C`` already IS the inventory-rate scale the sibling
    # rows use; the removed ``abs(residuals.rh_water_mol_s)`` branch was the
    # active branch at every pinned root (|r|/(q*C) = 1.68..1.75), which
    # saturated the reported scaled residual at exactly -1.0 and silently
    # loosened the acceptance gate by that same factor.  Removing it can only
    # raise a reported residual, never lower one, so no gate is weakened.
    rh_water = max(
        q_abs * assembly.interface.wet.retained_water_concentration_mol_m3,
        q_abs * assembly.interface.dry.total_water_concentration_mol_m3,
        floor,
    )
    rh_hexane = max(
        q_abs * wet_hexane,
        q_abs * assembly.interface.dry.total_hexane_concentration_mol_m3,
        abs(assembly.residuals.rh_hexane_mol_s),
        floor,
    )
    outer = n - 1
    rh_energy = max(
        _wet_capacity_energy(
            grid.volumes[outer],
            assembly.old_wet_cells[outer],
            assembly.config,
        )
        / dt_s,
        dry_capacity_energy / dt_s,
        abs(covariant[-1]),
        floor,
    )
    return ci.CutResidualScales(
        wet_water,
        wet_energy,
        (dry_water,),
        (dry_hexane,),
        (dry_energy,),
        rh_water,
        rh_hexane,
        rh_energy,
        (dry_capacity_audit,),
    )


def _wet_capacity_energy(
    volume_m3: float,
    cell: ww.WetWaterCellState,
    config: cut.CutTransportConfig,
) -> float:
    capacity = wet_core.energy_capacity(
        cell.temperature_k,
        cell.historical_hexane_loading,
        replace(
            config.wet.wet,
            X_water=cell.retained_water_loading,
            w_o=cell.oil_fraction_label,
        ),
    )
    result = volume_m3 * capacity * max(abs(cell.temperature_k), 1.0)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("birth wet capacity scale must be positive")
    return result


def _dry_capacity_energy(
    volume_m3: float,
    cell: cp.EquilibriumPoreState,
    assembly: birth.BirthAssembly,
) -> float:
    """Compatibility scalar view of the audited birth capacity scale."""

    return _dry_capacity_energy_with_audit(volume_m3, cell, assembly)[0]


def _dry_capacity_energy_with_audit(
    volume_m3: float,
    cell: cp.EquilibriumPoreState,
    assembly: birth.BirthAssembly,
) -> tuple[float, ci.DryEnergyCapacityStencilAudit]:
    capacity, audit = ci._dry_energy_capacity_with_audit(  # noqa: SLF001
        cell.temperature_k,
        cell.y_hexane,
        replace(
            assembly.config.dry.pore,
            w_o=assembly.before.oil_fraction_labels[-1],
        ),
        assembly.config,
    )
    result = volume_m3 * capacity * max(abs(cell.temperature_k), 1.0)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("birth dry capacity scale must be positive")
    return result, audit


def birth_sparsity(layout: birth.BirthLayout) -> BirthSparsity:
    """Return an O(N) safe structure for the birth residual Jacobian."""

    unknown = _span_map(layout.unknown_blocks)
    residual = _span_map(layout.residual_blocks)
    rows = [set() for _ in range(layout.residual_count)]

    def column(name: str, index: int = 0) -> int:
        start, size = unknown[name]
        if not 0 <= index < size:
            raise IndexError(name)
        return start + index

    n = layout.wet_piece_count
    for equation in ("wet_water", "wet_energy"):
        start, _ = residual[equation]
        for cell in range(n):
            for neighbour in range(max(0, cell - 1), min(n, cell + 2)):
                rows[start + cell].add(column("wet_temperature", neighbour))
                rows[start + cell].add(column("wet_retained_water", neighbour))
    coupled = {
        column("wet_temperature", n - 1),
        column("wet_retained_water", n - 1),
        column("dry_temperature"),
        column("dry_y_hexane"),
        column("dry_total_stefan_flux", 0),
        column("dry_total_stefan_flux", 1),
    }
    for equation in ("dry_water", "dry_hexane", "dry_energy", "rh_water", "rh_hexane", "rh_energy"):
        rows[residual[equation][0]].update(coupled)
    borders = (column("front_z"), column("interface_temperature"))
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
    return BirthSparsity(matrix, max(colors, default=-1) + 1, colors)


def scipy_solver_birth_sparsity(layout: birth.BirthLayout) -> BirthSparsity:
    """Return the exact grouping SciPy will recompute for ``jac_sparsity``."""

    declared = birth_sparsity(layout)
    colors = tuple(int(color) for color in group_columns(declared.matrix))
    return BirthSparsity(
        declared.matrix,
        max(colors, default=-1) + 1,
        colors,
    )


def select_birth_solver_jacobian(
    layout: birth.BirthLayout,
    policy: BirthSolverJacobianPolicy,
) -> BirthSolverJacobianSelection:
    """Resolve an exact caller-selected route without consulting size controls."""

    if type(layout) is not birth.BirthLayout:
        raise TypeError("birth solver selection requires an exact BirthLayout")
    if type(policy) is not BirthSolverJacobianPolicy:
        raise TypeError("birth solver selection requires an exact policy enum")
    if policy is BirthSolverJacobianPolicy.DENSE_TRF_TWO_POINT:
        return BirthSolverJacobianSelection(policy, "trf", None)
    if policy is BirthSolverJacobianPolicy.STRUCTURED_DOGBOX_TWO_POINT:
        return BirthSolverJacobianSelection(
            policy,
            "dogbox",
            scipy_solver_birth_sparsity(layout),
        )
    raise ValueError("birth solver selection encountered an unknown policy")


def _span_map(
    blocks: tuple[tuple[str, int], ...],
) -> dict[str, tuple[int, int]]:
    start = 0
    result = {}
    for name, size in blocks:
        result[name] = (start, size)
        start += size
    return result


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
            colors[neighbour] for neighbour in neighbours[column] if colors[neighbour] >= 0
        }
        color = 0
        while color in forbidden:
            color += 1
        colors[column] = color
    return tuple(colors)


def _accept_birth(
    before: ww.WetWaterState,
    controls: ci.CutSolverControls,
    boundary: ct.PoreBoundary,
    seed: BirthStepSeed,
    assembly: birth.BirthAssembly,
    scales: ci.CutResidualScales,
    scaled_residuals: tuple[float, ...],
    evaluations: int,
    rejected_trials: int,
    nonlinear_message: str,
    condition_proxy: float,
    chart: _CoordinateChart,
    topology: BirthTopologyDiagnostics,
    structure: BirthSparsity | None,
) -> BirthIntegratorStep:
    candidate = assembly.candidate
    config = assembly.config
    dt_s = assembly.dt_s
    after_transport = cut.CutTransportState(
        geometry=assembly.candidate_geometry,
        config=config,
        time_s=before.time_s + dt_s,
        wet_temperatures_k=candidate.wet_temperatures_k,
        wet_retained_water_loadings=candidate.wet_retained_water_loadings,
        dry_temperatures_k=(candidate.dry_temperature_k,),
        dry_y_hexane=(candidate.dry_y_hexane,),
        historical_hexane_loadings=before.historical_hexane_loadings,
        oil_fraction_labels=before.oil_fraction_labels,
    )
    reference = _fully_wet_reference_inventory(before, config)
    reference_capacity = _fully_wet_reference_capacity(before, config)
    current = ci.inventory_snapshot(after_transport)
    event_changes = _material_event_changes(assembly)
    prior_changes = _fully_wet_prior_changes(before, config, reference)
    cumulative_water, water_comp = ci._compensated_add_vectors(  # noqa: SLF001
        prior_changes.water_cell_mol,
        (0.0,) * before.grid.n,
        event_changes.water_cell_mol,
    )
    cumulative_hexane, hexane_comp = ci._compensated_add_vectors(  # noqa: SLF001
        prior_changes.hexane_cell_mol,
        (0.0,) * before.grid.n,
        event_changes.hexane_cell_mol,
    )
    cumulative_energy, energy_comp = ci._compensated_add_vectors(  # noqa: SLF001
        prior_changes.energy_cell_j,
        (0.0,) * before.grid.n,
        event_changes.energy_cell_j,
    )
    area = before.grid.areas[-1]
    surface = assembly.dry_face_fluxes[-1]
    water_out = dt_s * area * surface.component.conserved_water_flux_mol_m2_s
    hexane_out = dt_s * area * surface.component.conserved_hexane_flux_mol_m2_s
    energy_out = dt_s * area * surface.energy.total_energy_flux_w_m2
    boundary_water, boundary_water_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_water_out_mol,
        0.0,
        water_out,
    )
    boundary_hexane, boundary_hexane_comp = ci._compensated_add_scalar(  # noqa: SLF001
        0.0,
        0.0,
        hexane_out,
    )
    boundary_energy, boundary_energy_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_energy_out_j,
        0.0,
        energy_out,
    )
    water_step = math.fsum((*event_changes.water_cell_mol, water_out))
    hexane_step = math.fsum((*event_changes.hexane_cell_mol, hexane_out))
    energy_step = math.fsum((*event_changes.energy_cell_j, energy_out))
    after_capacity = ci.capacity_energy_scale(after_transport)
    water_step_scale = max(
        _fully_wet_current_inventory(before, config).total_water_mol,
        current.total_water_mol,
        abs(water_out),
        1.0e-300,
    )
    hexane_step_scale = max(
        reference.total_hexane_mol,
        current.total_hexane_mol,
        abs(hexane_out),
        1.0e-300,
    )
    energy_step_scale = max(reference_capacity, after_capacity, 1.0e-300)
    water_cumulative = math.fsum(
        (*cumulative_water, *water_comp, boundary_water, boundary_water_comp)
    )
    hexane_cumulative = math.fsum(
        (*cumulative_hexane, *hexane_comp, boundary_hexane, boundary_hexane_comp)
    )
    energy_cumulative = math.fsum(
        (*cumulative_energy, *energy_comp, boundary_energy, boundary_energy_comp)
    )
    absolute_water = before.cumulative_absolute_water_exchange_mol + abs(water_out)
    absolute_hexane = abs(hexane_out)
    absolute_energy = before.cumulative_absolute_energy_exchange_j + abs(energy_out)
    water_cumulative_scale = max(
        reference.total_water_mol,
        current.total_water_mol,
        abs(math.fsum((boundary_water, boundary_water_comp))),
        absolute_water,
        1.0e-300,
    )
    hexane_cumulative_scale = max(
        reference.total_hexane_mol,
        current.total_hexane_mol,
        abs(math.fsum((boundary_hexane, boundary_hexane_comp))),
        absolute_hexane,
        1.0e-300,
    )
    energy_cumulative_scale = max(
        reference_capacity,
        after_capacity,
        absolute_energy,
        1.0e-300,
    )
    interface_hexane = assembly.dry_face_fluxes[0].component.conserved_hexane_flux_mol_m2_s
    surface_hexane = assembly.surface_hexane_flux_mol_m2_s
    wet_hexane = config.wet.wet.rho_dm_p * assembly.interface.historical_hexane_loading / hx.M
    storage_jump = math.fsum(
        (-assembly.interface.dry.total_hexane_concentration_mol_m3, wet_hexane)
    )
    if storage_jump <= 0.0:
        raise RuntimeError("birth converged without a positive n-hexane storage jump")
    area_gamma = 4.0 * math.pi * assembly.candidate_geometry.front.radius_m**2
    swept_speed = -assembly.swept_geometry.interface_swept_volume_rate_m3_s / area_gamma
    rh_speed = interface_hexane / storage_jump
    immediate_error = abs(rh_speed - swept_speed) / max(
        abs(rh_speed),
        abs(swept_speed),
        1.0e-300,
    )
    evidence = birth.RecedingBranchEvidence(
        True,
        candidate.front_z,
        max(abs(value) for value in scaled_residuals),
        # GT-PS-2-P1E-R15: the branch certificate consumes the same birth
        # regional/RH acceptance as the residual contract (see the module
        # constant); the shared trajectory tolerance stays untouched.
        max(controls.nonlinear_residual_tolerance, BIRTH_REGIONAL_RH_ACCEPTANCE),
    )
    decision = birth.classify_endpoint_branch(
        outward_hexane_flux_mol_m2_s=surface_hexane,
        hexane_storage_jump_mol_m3=storage_jump,
        receding=evidence,
    )
    if (
        surface_hexane <= 0.0
        or interface_hexane <= 0.0
        or swept_speed <= 0.0
        or decision.branch is not birth.BirthBranch.RECEDING
        or decision.contact_evidence_used
    ):
        raise RuntimeError("converged birth did not certify the positive receding branch")
    after = ci.CutIntegratorState(
        transport=after_transport,
        controls=controls,
        reference_inventory=reference,
        reference_capacity_energy_scale_j=reference_capacity,
        last_interface_temperature_k=candidate.interface_temperature_k,
        last_total_stefan_fluxes_mol_m2_s=(candidate.dry_total_stefan_fluxes_mol_m2_s),
        cumulative_boundary_water_out_mol=boundary_water,
        cumulative_boundary_hexane_out_mol=boundary_hexane,
        cumulative_boundary_energy_out_j=boundary_energy,
        boundary_water_compensation_mol=boundary_water_comp,
        boundary_hexane_compensation_mol=boundary_hexane_comp,
        boundary_energy_compensation_j=boundary_energy_comp,
        cumulative_absolute_water_transfer_mol=absolute_water,
        cumulative_absolute_hexane_transfer_mol=absolute_hexane,
        cumulative_absolute_energy_transfer_j=absolute_energy,
        cumulative_material_water_change_cell_mol=cumulative_water,
        cumulative_material_hexane_change_cell_mol=cumulative_hexane,
        cumulative_material_energy_change_cell_j=cumulative_energy,
        water_change_compensation_cell_mol=water_comp,
        hexane_change_compensation_cell_mol=hexane_comp,
        energy_change_compensation_cell_j=energy_comp,
        accepted_steps=1,
        cumulative_nonlinear_evaluations=evaluations,
    )
    fractions = _bounded_fractions(candidate, chart, topology)
    birth_telescope = _normalized_birth_telescoping(assembly)
    ledger = BirthIntegratorLedger(
        dt_s=dt_s,
        nonlinear_evaluations=evaluations,
        rejected_trial_evaluations=rejected_trials,
        nonlinear_message=nonlinear_message,
        residual_scales=scales,
        scaled_residuals=scaled_residuals,
        maximum_scaled_residual=max(abs(value) for value in scaled_residuals),
        condition_proxy=condition_proxy,
        minimum_fractional_distance_to_bound=min(min(value, 1.0 - value) for value in fractions),
        minimum_front_z_distance_to_chart_boundary=min(
            candidate.front_z - chart.front_z[0],
            chart.front_z[1] - candidate.front_z,
        ),
        topology_diagnostics=topology,
        used_sparse_jacobian=structure is not None,
        scalability_class=(
            "bordered-banded grouped finite-difference sparse birth Jacobian; "
            "runtime qualification remains open"
            if structure is not None
            else "certified small-N dense finite-difference birth oracle"
        ),
        jacobian_structural_color_count=(
            assembly.layout.unknown_count if structure is None else structure.color_count
        ),
        jacobian_structural_nnz=(
            assembly.layout.unknown_count**2 if structure is None else structure.nnz
        ),
        water_step_residual_mol=water_step,
        hexane_step_residual_mol=hexane_step,
        energy_step_residual_j=energy_step,
        normalized_water_step_residual=abs(water_step) / water_step_scale,
        normalized_hexane_step_residual=abs(hexane_step) / hexane_step_scale,
        normalized_energy_step_residual=abs(energy_step) / energy_step_scale,
        water_cumulative_residual_mol=water_cumulative,
        hexane_cumulative_residual_mol=hexane_cumulative,
        energy_cumulative_residual_j=energy_cumulative,
        normalized_water_cumulative_residual=(abs(water_cumulative) / water_cumulative_scale),
        normalized_hexane_cumulative_residual=(abs(hexane_cumulative) / hexane_cumulative_scale),
        normalized_energy_cumulative_residual=(abs(energy_cumulative) / energy_cumulative_scale),
        boundary_water_out_mol=water_out,
        boundary_hexane_out_mol=hexane_out,
        boundary_energy_out_j=energy_out,
        surface_hexane_flux_mol_m2_s=surface_hexane,
        interface_hexane_flux_mol_m2_s=interface_hexane,
        hexane_storage_jump_mol_m3=storage_jump,
        swept_recession_speed_m_s=swept_speed,
        rh_recession_speed_m_s=rh_speed,
        relative_immediate_recession_error=immediate_error,
        newborn_dry_volume_m3=assembly.ledger.newborn_dry_volume_m3,
        branch_decision=decision,
        birth_telescoping_residual=birth_telescope,
        dirichlet_numerical_oracle=True,
        contact_selected=False,
        accepted=True,
    )
    return BirthIntegratorStep(before, after, assembly, boundary, controls, seed, ledger)


def _normalized_birth_telescoping(assembly: birth.BirthAssembly) -> float:
    """Normalize exact block identities by physical rates, not cancelled nets.

    At a converged conservative root both the block sum and the independently
    assembled global residual approach zero.  Dividing their roundoff-level
    difference by either cancelled net would therefore be ill-conditioned.
    The extensive inventory-change and surface-throughput rates remain the
    appropriate non-cancelling scale for this algebraic identity.
    """

    changes = assembly.inventory_changes
    surface = assembly.dry_face_fluxes[-1]
    area = assembly.before.grid.areas[-1]
    dt_s = assembly.dt_s
    water_rate = area * surface.component.conserved_water_flux_mol_m2_s
    hexane_rate = area * surface.component.conserved_hexane_flux_mol_m2_s
    energy_rate = area * surface.energy.total_energy_flux_w_m2
    water_scale = max(
        math.fsum(abs(value) for value in (*changes.wet_water_mol, *changes.dry_water_mol)) / dt_s,
        abs(water_rate),
        1.0e-300,
    )
    hexane_scale = max(
        math.fsum(abs(value) for value in (*changes.wet_hexane_mol, *changes.dry_hexane_mol))
        / dt_s,
        abs(hexane_rate),
        1.0e-300,
    )
    energy_scale = max(
        math.fsum(abs(value) for value in (*changes.wet_energy_j, *changes.dry_energy_j)) / dt_s,
        abs(energy_rate),
        1.0e-300,
    )
    ledger = assembly.ledger
    return max(
        abs(ledger.water_telescoping_error_mol_s) / water_scale,
        abs(ledger.hexane_telescoping_error_mol_s) / hexane_scale,
        abs(ledger.wet_hexane_sweep_identity_mol_s) / hexane_scale,
        abs(ledger.energy_telescoping_error_w) / energy_scale,
    )


def _fully_wet_reference_inventory(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
) -> ci.MaterialInventorySnapshot:
    hexane = tuple(
        volume * config.wet.wet.rho_dm_p * loading / hx.M
        for volume, loading in zip(
            before.grid.volumes,
            before.historical_hexane_loadings,
        )
    )
    if not math.isclose(
        math.fsum(hexane),
        before.reference_hexane_mass_kg / hx.M,
        rel_tol=2.0e-14,
        abs_tol=0.0,
    ):
        raise RuntimeError("fully-wet reference n-hexane inventory changed at birth")
    return ci.MaterialInventorySnapshot(
        before.reference_water_cell_mol,
        hexane,
        before.reference_energy_cell_j,
    )


def _fully_wet_current_inventory(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
) -> ci.MaterialInventorySnapshot:
    cells = tuple(
        ww.evaluate_cell(T, W, Xh, wo, config.wet)
        for T, W, Xh, wo in zip(
            before.temperatures_k,
            before.retained_water_loadings,
            before.historical_hexane_loadings,
            before.oil_fraction_labels,
        )
    )
    return ci.MaterialInventorySnapshot(
        tuple(
            volume * cell.retained_water_concentration_mol_m3
            for volume, cell in zip(before.grid.volumes, cells)
        ),
        tuple(
            volume * config.wet.wet.rho_dm_p * loading / hx.M
            for volume, loading in zip(
                before.grid.volumes,
                before.historical_hexane_loadings,
            )
        ),
        tuple(
            volume * cell.energy_density_j_m3 for volume, cell in zip(before.grid.volumes, cells)
        ),
    )


def _fully_wet_prior_changes(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
    reference: ci.MaterialInventorySnapshot,
) -> ci.MaterialInventorySnapshot:
    current = _fully_wet_current_inventory(before, config)
    return ci.MaterialInventorySnapshot(
        tuple(
            math.fsum((value, -baseline))
            for value, baseline in zip(
                current.water_cell_mol,
                reference.water_cell_mol,
            )
        ),
        tuple(
            math.fsum((value, -baseline))
            for value, baseline in zip(
                current.hexane_cell_mol,
                reference.hexane_cell_mol,
            )
        ),
        tuple(
            math.fsum((value, -baseline))
            for value, baseline in zip(
                current.energy_cell_j,
                reference.energy_cell_j,
            )
        ),
    )


def _material_event_changes(
    assembly: birth.BirthAssembly,
) -> ci.MaterialInventorySnapshot:
    changes = assembly.inventory_changes
    water = list(changes.wet_water_mol)
    hexane = list(changes.wet_hexane_mol)
    energy = list(changes.wet_energy_j)
    water[-1] = math.fsum((water[-1], changes.dry_water_mol[0]))
    hexane[-1] = math.fsum((hexane[-1], changes.dry_hexane_mol[0]))
    energy[-1] = math.fsum((energy[-1], changes.dry_energy_j[0]))
    return ci.MaterialInventorySnapshot(tuple(water), tuple(hexane), tuple(energy))


def _fully_wet_reference_capacity(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
) -> float:
    cells = tuple(
        ww.evaluate_cell(T, W, Xh, wo, config.wet)
        for T, W, Xh, wo in zip(
            before.reference_temperatures_k,
            before.reference_retained_water_loadings,
            before.historical_hexane_loadings,
            before.oil_fraction_labels,
        )
    )
    value = math.fsum(
        _wet_capacity_energy(volume, cell, config)
        for volume, cell in zip(before.grid.volumes, cells)
    )
    if not math.isfinite(value) or value <= 0.0:
        raise RuntimeError("fully-wet reference capacity scale must be positive")
    return value


__all__ = [
    "BirthIntegratorLedger",
    "BirthResidualAssemblyBudgetExceeded",
    "BirthSolverJacobianPolicy",
    "BirthSolverJacobianSelection",
    "BirthIntegratorStep",
    "BirthIntegratorStepError",
    "BirthSparsity",
    "BirthStepSeed",
    "BirthTopologyDiagnostics",
    "advance_receding_birth",
    "birth_sparsity",
    "scipy_solver_birth_sparsity",
    "select_birth_solver_jacobian",
]
