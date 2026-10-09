r"""Fixed-time fully-dry projection at the exact ``z=0`` seam.

The exact endpoint audit in :mod:`cut_extinction_event` proves that an
unknown extinction time leaves one scalar degree of freedom.  This module
does not fill that degree with an invented endpoint condition.  It accepts a
*supplied* positive duration, removes it from the nonlinear coordinates, and
solves the resulting square ``3*N`` conservative projection.

Consequently an accepted result certifies only the fixed-time projection and
its component/common-datum-energy ledgers.  It is not evidence that the
supplied duration is the physical extinction time.  A production trajectory
must obtain that duration by integration and root localisation on the strict
positive-core chart before calling this projection.

All coefficient and Dirichlet-boundary results remain numerical oracles and
report ``physically_qualifying=False``.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field, replace
from typing import Sequence

import numpy as np
from scipy import optimize, sparse

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_extinction_event as extinction
from dtdc_simulator.core2.particle import cut_integrator as ci
from dtdc_simulator.core2.props import sorption as sp


class ExtinctionProjectionCompatibilityError(RuntimeError):
    """The preserved endpoint labels cannot use a reduced legacy dry state."""


class ExtinctionProjectionStepError(RuntimeError):
    """Rejected fixed-time projection carrying the exact rollback object."""

    def __init__(
        self,
        message: str,
        rollback_state: ci.CutIntegratorState,
        *,
        nonlinear_evaluations: int,
        maximum_scaled_residual: float = math.inf,
        condition_proxy: float = math.inf,
        last_candidate: extinction.ExtinctionUnknowns | None = None,
        retained_water_applicability_guard_failed: bool = False,
        retained_water_applicability_audit: (
            ci.WetRetainedWaterApplicabilityAudit | None
        ) = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.nonlinear_evaluations = nonlinear_evaluations
        self.maximum_scaled_residual = maximum_scaled_residual
        self.condition_proxy = condition_proxy
        self.last_candidate = last_candidate
        self.retained_water_applicability_guard_failed = (
            retained_water_applicability_guard_failed
        )
        self.retained_water_applicability_audit = (
            retained_water_applicability_audit
        )


@dataclass(frozen=True)
class ExtinctionProjectionSeed:
    """One fixed-time endpoint seed and explicit positive Stefan scales."""

    candidate: extinction.ExtinctionUnknowns
    stefan_flux_scales_mol_m2_s: tuple[float, ...]
    label: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.candidate, extinction.ExtinctionUnknowns):
            raise TypeError("projection seed requires ExtinctionUnknowns")
        if not self.stefan_flux_scales_mol_m2_s or not all(
            math.isfinite(value) and value > 0.0
            for value in self.stefan_flux_scales_mol_m2_s
        ):
            raise ValueError("every positive-area Stefan face needs a positive scale")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("projection seed needs a non-empty provenance label")


@dataclass(frozen=True)
class ExtinctionProjectionResidualScales:
    """Fixed local component and datum-independent capacity scales."""

    dry_water_mol_s: tuple[float, ...]
    dry_hexane_mol_s: tuple[float, ...]
    dry_energy_w: tuple[float, ...]
    dry_energy_capacity_stencil_audits: tuple[
        ci.DryEnergyCapacityStencilAudit,
        ...,
    ]
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.dry_energy_capacity_stencil_audits, tuple):
            raise TypeError(
                "extinction dry-capacity audits must be an immutable tuple"
            )
        if len(self.dry_energy_capacity_stencil_audits) != len(self.dry_energy_w):
            raise ValueError(
                "every extinction dry-energy scale needs exactly one "
                "Amendment-19 stencil audit"
            )
        if not all(
            isinstance(audit, ci.DryEnergyCapacityStencilAudit)
            for audit in self.dry_energy_capacity_stencil_audits
        ):
            raise TypeError("extinction dry-capacity audit payload has the wrong type")

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
        values = (
            *self.dry_water_mol_s,
            *self.dry_hexane_mol_s,
            *self.dry_energy_w,
        )
        if not values or not all(
            math.isfinite(value) and value > 0.0 for value in values
        ):
            raise RuntimeError("every extinction-projection row scale must be positive")
        return values


@dataclass(frozen=True)
class ExtinctionProjectionTopologyDiagnostics:
    """Exact per-cell gas-only intervals and strict activity margins."""

    uses_conditioned_composition_chart: bool
    oil_fraction_labels: tuple[float, ...]
    exact_composition_intervals: tuple[cp.GasOnlyCompositionInterval, ...]
    exact_composition_interval_fractions: tuple[float, ...]
    water_saturation_activity_gaps: tuple[float, ...]
    hexane_saturation_activity_gaps: tuple[float, ...]
    luikov_lower_activity_gaps: tuple[float, ...]
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        count = len(self.oil_fraction_labels)
        if count < 1 or not all(
            len(values) == count
            for values in (
                self.exact_composition_intervals,
                self.exact_composition_interval_fractions,
                self.water_saturation_activity_gaps,
                self.hexane_saturation_activity_gaps,
                self.luikov_lower_activity_gaps,
            )
        ):
            raise ValueError("one exact extinction topology diagnostic is required per cell")
        if not all(math.isfinite(value) for value in self.oil_fraction_labels):
            raise ValueError("extinction oil labels must remain finite")
        if not all(
            isinstance(interval, cp.GasOnlyCompositionInterval)
            for interval in self.exact_composition_intervals
        ):
            raise TypeError("extinction topology diagnostics need exact gas-only intervals")
        if not all(
            math.isfinite(value) and 0.0 < value < 1.0
            for value in self.exact_composition_interval_fractions
        ):
            raise ValueError("every endpoint composition must be inside its exact interval")
        if not all(
            math.isfinite(value) and value > 0.0
            for values in (
                self.water_saturation_activity_gaps,
                self.hexane_saturation_activity_gaps,
                self.luikov_lower_activity_gaps,
            )
            for value in values
        ):
            raise ValueError("every exact endpoint activity gap must remain strictly positive")

    @property
    def minimum_exact_activity_gap(self) -> float:
        return min(
            *self.water_saturation_activity_gaps,
            *self.hexane_saturation_activity_gaps,
            *self.luikov_lower_activity_gaps,
        )

    @property
    def minimum_exact_composition_fractional_distance(self) -> float:
        return min(
            min(fraction, 1.0 - fraction)
            for fraction in self.exact_composition_interval_fractions
        )


@dataclass(frozen=True)
class ExtinctionProjectionLedger:
    """Nonlinear, structural, step, and cumulative conservation evidence."""

    event_duration_s: float
    nonlinear_evaluations: int
    nonlinear_message: str
    residual_scales: ExtinctionProjectionResidualScales
    scaled_residuals: tuple[float, ...]
    maximum_scaled_residual: float
    condition_proxy: float
    minimum_fractional_distance_to_bound: float
    topology_diagnostics: ExtinctionProjectionTopologyDiagnostics
    used_sparse_jacobian: bool
    jacobian_structural_nnz: int
    jacobian_structural_color_count: int
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
    event_time_was_localized_here: bool = field(default=False, init=False)
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
        return (
            self.topology_diagnostics.minimum_exact_composition_fractional_distance
        )


@dataclass(frozen=True)
class ExtinctionCommittedState:
    """Exact dry endpoint with every material label and compensated history.

    ``historical_hexane_loadings`` remains an audit label after dry-out.  The
    legacy fully-dry predictor can be used only when every residual-oil label
    is exactly the single value in its configuration; otherwise conversion
    fails rather than silently homogenising material.
    """

    grid: object
    config: ct.FullyDryTransportConfig
    controls: ci.CutSolverControls
    time_s: float
    temperatures_k: tuple[float, ...]
    y_hexane: tuple[float, ...]
    historical_hexane_loadings: tuple[float, ...]
    oil_fraction_labels: tuple[float, ...]
    reference_inventory: ci.MaterialInventorySnapshot
    reference_capacity_energy_scale_j: float
    last_positive_area_total_stefan_fluxes_mol_m2_s: tuple[float, ...]
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
    transition_authority: str = field(
        default=(
            "fixed-time exact-z=0 projection after separately supplied positive-core "
            "event duration; no endpoint RH equation and no epsilon core"
        ),
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        n = self.grid.n
        for values in (
            self.temperatures_k,
            self.y_hexane,
            self.historical_hexane_loadings,
            self.oil_fraction_labels,
            self.last_positive_area_total_stefan_fluxes_mol_m2_s,
            self.cumulative_material_water_change_cell_mol,
            self.cumulative_material_hexane_change_cell_mol,
            self.cumulative_material_energy_change_cell_j,
            self.water_change_compensation_cell_mol,
            self.hexane_change_compensation_cell_mol,
            self.energy_change_compensation_cell_j,
        ):
            if len(values) != n:
                raise ValueError("extinction committed arrays must align with the grid")
            if not all(math.isfinite(value) for value in values):
                raise ValueError("extinction committed arrays must remain finite")
        if not math.isfinite(self.time_s) or self.time_s <= 0.0:
            raise ValueError("extinction committed time must be positive and finite")
        if self.accepted_steps < 1 or self.cumulative_nonlinear_evaluations < 1:
            raise ValueError("extinction commit must include accepted work counters")

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

    def to_uniform_fully_dry_transport_state(self) -> ct.FullyDryTransportState:
        """Return the existing dry predictor state only without label loss."""

        configured = self.config.pore.w_o
        if any(value != configured for value in self.oil_fraction_labels):
            raise ExtinctionProjectionCompatibilityError(
                "the existing fully-dry predictor has one oil label; a nonuniform "
                "endpoint must not be homogenised"
            )

        def corrected(high: Sequence[float], compensation: Sequence[float]) -> float:
            return math.fsum(
                math.fsum((value, correction))
                for value, correction in zip(high, compensation)
            )

        return ct.FullyDryTransportState(
            grid=self.grid,
            config=self.config,
            time_s=self.time_s,
            temperatures_k=self.temperatures_k,
            y_hexane=self.y_hexane,
            reference_water_mol=self.reference_inventory.total_water_mol,
            reference_hexane_mol=self.reference_inventory.total_hexane_mol,
            reference_energy_j=self.reference_inventory.total_energy_j,
            cumulative_boundary_water_out_mol=(
                self.corrected_cumulative_boundary_water_out_mol
            ),
            cumulative_boundary_hexane_out_mol=(
                self.corrected_cumulative_boundary_hexane_out_mol
            ),
            cumulative_boundary_energy_out_j=(
                self.corrected_cumulative_boundary_energy_out_j
            ),
            cumulative_water_inventory_change_mol=corrected(
                self.cumulative_material_water_change_cell_mol,
                self.water_change_compensation_cell_mol,
            ),
            cumulative_hexane_inventory_change_mol=corrected(
                self.cumulative_material_hexane_change_cell_mol,
                self.hexane_change_compensation_cell_mol,
            ),
            cumulative_energy_inventory_change_j=corrected(
                self.cumulative_material_energy_change_cell_j,
                self.energy_change_compensation_cell_j,
            ),
            cumulative_absolute_water_transfer_mol=(
                self.cumulative_absolute_water_transfer_mol
            ),
            cumulative_absolute_hexane_transfer_mol=(
                self.cumulative_absolute_hexane_transfer_mol
            ),
            cumulative_absolute_energy_transfer_j=(
                self.cumulative_absolute_energy_transfer_j
            ),
        )


@dataclass(frozen=True)
class ExtinctionProjectionStep:
    """One accepted fixed-time endpoint projection and exact dry payload."""

    before: ci.CutIntegratorState
    after: ExtinctionCommittedState
    assembly: extinction.ExtinctionAssembly
    boundary: ct.PoreBoundary
    seed: ExtinctionProjectionSeed
    ledger: ExtinctionProjectionLedger
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class _ProjectionChart:
    dry_temperature: tuple[float, float]
    legacy_dry_y_hexane: tuple[float, float] | None
    dry_config: ct.FullyDryTransportConfig
    dry_pores: tuple[cp.CoupledPoreParams, ...]


@dataclass(frozen=True)
class _ProjectionSparsity:
    matrix: sparse.csr_matrix
    column_colors: tuple[int, ...]

    @property
    def nnz(self) -> int:
        return int(self.matrix.nnz)

    @property
    def color_count(self) -> int:
        return max(self.column_colors, default=-1) + 1


def solve_fixed_time_extinction_projection(
    before: ci.CutIntegratorState,
    event_duration_s: float,
    boundary: ct.PoreBoundary,
    seed: ExtinctionProjectionSeed,
    *,
    environment_conditioned_resolution: object | None = None,
) -> ExtinctionProjectionStep:
    """Solve a supplied-duration endpoint projection; never localize its time.

    ``environment_conditioned_resolution`` defaults to ``None`` - the frozen
    residual contract, bit-identical to the pre-declaration path.  When the
    P-1-ruled typed declaration (``cut_extinction_time.
    EnvironmentConditionedResidualResolution``) arrives through the certified
    certificate hand-off, ONLY the residual-magnitude acceptance clause reads
    the declared bound; success, logit-guard, condition-proxy, ledger, and
    retained-water gates are untouched.
    """

    evaluations = 0
    maximum_scaled = math.inf
    condition_proxy = math.inf
    last_candidate = None
    try:
        if not isinstance(before, ci.CutIntegratorState):
            raise TypeError("before must be CutIntegratorState")
        if environment_conditioned_resolution is None:
            residual_limit = before.controls.nonlinear_residual_tolerance
            evaluation_budget = before.controls.maximum_function_evaluations
        else:
            # Runtime import: cut_extinction_time imports this module at top
            # level, so the typed declaration class is resolved lazily here.
            from dtdc_simulator.core2.particle import (
                cut_extinction_time as _extinction_time,
            )

            residual_limit = _extinction_time._residual_limit_in_force(  # noqa: SLF001
                before,
                environment_conditioned_resolution,
            )
            evaluation_budget = _extinction_time._evaluation_budget_in_force(  # noqa: SLF001
                before,
                environment_conditioned_resolution,
            )
        if not math.isfinite(event_duration_s) or event_duration_s <= 0.0:
            raise ValueError("fixed extinction duration must be positive and finite")
        if not isinstance(
            boundary,
            (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
        ):
            raise TypeError("boundary must be a supported pore boundary")
        ci.certify_explicit_wet_retained_water_applicability(
            before.controls,
            before.transport.wet_retained_water_loadings,
            before.transport.effective_wet_retained_water_capacity_duals_over_rt,
            before.transport.config.wet.luikov,
        )
        layout = extinction.layout_for_extinction(before.transport)
        chart = _coordinate_chart(before)
        _validate_seed(event_duration_s, seed, layout, chart)
        initial_coordinates = _encode_candidate(
            seed.candidate,
            chart,
            seed.stefan_flux_scales_mol_m2_s,
        )
        initial = extinction.assemble_extinction_endpoint(
            before.transport,
            seed.candidate,
            boundary,
            enforce_reduced_film_thresholds=False,
        )
        evaluations += 1
        scales = _residual_scales(before, initial, event_duration_s)
        scale_vector = np.asarray(scales.vector, dtype=float)
        use_sparse = layout.fixed_time_unknown_count > before.controls.maximum_dense_unknowns
        structural = _projection_sparsity(layout.dry_cell_count) if use_sparse else None

        def scaled_residual(coordinates: np.ndarray) -> np.ndarray:
            nonlocal evaluations
            candidate = _decode_candidate(
                coordinates,
                chart,
                seed.stefan_flux_scales_mol_m2_s,
                event_duration_s,
            )
            assembly = extinction.assemble_extinction_endpoint(
                before.transport,
                candidate,
                boundary,
                enforce_reduced_film_thresholds=False,
            )
            evaluations += 1
            return np.asarray(
                assembly.residuals.datum_covariant_vector,
                dtype=float,
            ) / scale_vector

        n = layout.dry_cell_count
        guard = before.controls.maximum_logit_magnitude
        lower = np.asarray((*([-guard] * (2 * n)), *([-np.inf] * n)))
        upper = np.asarray((*([guard] * (2 * n)), *([np.inf] * n)))
        solution = optimize.least_squares(
            scaled_residual,
            initial_coordinates,
            method="dogbox" if structural is not None else "trf",
            bounds=(lower, upper),
            ftol=before.controls.nonlinear_step_tolerance,
            xtol=before.controls.nonlinear_step_tolerance,
            gtol=before.controls.nonlinear_step_tolerance,
            max_nfev=evaluation_budget,
            x_scale="jac",
            jac_sparsity=None if structural is None else structural.matrix,
            tr_options=(
                {"atol": 1.0e-14, "btol": 1.0e-14, "maxiter": 1000}
                if structural is not None
                else None
            ),
        )
        candidate = _decode_candidate(
            solution.x,
            chart,
            seed.stefan_flux_scales_mol_m2_s,
            event_duration_s,
        )
        converged = extinction.assemble_extinction_endpoint(
            before.transport,
            candidate,
            boundary,
        )
        evaluations += 1
        signed_scaled = tuple(
            residual / scale
            for residual, scale in zip(
                converged.residuals.datum_covariant_vector,
                scales.vector,
            )
        )
        maximum_scaled = max(abs(value) for value in signed_scaled)
        condition_proxy = ci._jacobian_condition_proxy(solution.jac)  # noqa: SLF001
        last_candidate = candidate
        topology_diagnostics = _topology_diagnostics(
            before,
            candidate,
            converged.dry_cells,
            chart,
        )
        fractions = _bounded_fractions(candidate, chart, topology_diagnostics)
        minimum_fraction = min(min(value, 1.0 - value) for value in fractions)
        if not solution.success:
            raise RuntimeError(
                f"fixed-time endpoint root did not certify success: {solution.message}"
            )
        if np.any(solution.active_mask != 0):
            raise RuntimeError("fixed-time endpoint touched the numerical logit guard")
        if maximum_scaled > residual_limit:
            if environment_conditioned_resolution is None:
                raise RuntimeError(
                    "fixed-time endpoint residual contract failed: "
                    f"{maximum_scaled:.3e} > "
                    f"{before.controls.nonlinear_residual_tolerance:.3e}"
                )
            raise RuntimeError(
                "fixed-time endpoint residual exceeds even the declared "
                "environment-conditioned resolution: "
                f"{maximum_scaled:.3e} > declared bound "
                f"{residual_limit:.3e} (frozen contract "
                f"{before.controls.nonlinear_residual_tolerance:.3e})"
            )
        if not math.isfinite(condition_proxy) or (
            condition_proxy > before.controls.maximum_condition_proxy
        ):
            raise RuntimeError(
                "fixed-time endpoint Jacobian condition proxy is uncertified: "
                f"{condition_proxy:.3e}"
            )
        step = _accept_projection(
            before,
            boundary,
            seed,
            converged,
            scales,
            signed_scaled,
            evaluations,
            str(solution.message),
            condition_proxy,
            minimum_fraction,
            topology_diagnostics,
            structural,
        )
        if (
            step.ledger.maximum_step_ledger_residual
            > before.controls.ledger_tolerance
            or step.ledger.maximum_cumulative_ledger_residual
            > before.controls.ledger_tolerance
        ):
            raise RuntimeError(
                "fixed-time endpoint accepted/cumulative ledger exceeded 1e-10"
            )
        ci.certify_explicit_wet_retained_water_applicability(
            before.controls,
            before.transport.wet_retained_water_loadings,
            before.transport.effective_wet_retained_water_capacity_duals_over_rt,
            before.transport.config.wet.luikov,
        )
        return step
    except ExtinctionProjectionStepError:
        raise
    except Exception as exc:
        if not isinstance(before, ci.CutIntegratorState):
            raise
        applicability_error = (
            exc
            if isinstance(exc, ci.WetRetainedWaterApplicabilityError)
            else None
        )
        raise ExtinctionProjectionStepError(
            f"fixed-time extinction projection rejected with exact rollback: {exc}",
            before,
            nonlinear_evaluations=evaluations,
            maximum_scaled_residual=maximum_scaled,
            condition_proxy=condition_proxy,
            last_candidate=last_candidate,
            retained_water_applicability_guard_failed=(
                applicability_error is not None
            ),
            retained_water_applicability_audit=(
                None if applicability_error is None else applicability_error.audit
            ),
        ) from exc


def _validate_seed(
    event_duration_s: float,
    seed: ExtinctionProjectionSeed,
    layout: extinction.ExtinctionLayout,
    chart: _ProjectionChart,
) -> None:
    if seed.candidate.event_time_s != event_duration_s:
        raise ValueError(
            "the seed duration must equal the separately supplied fixed event duration"
        )
    if len(seed.candidate.vector()) != layout.unknown_count:
        raise ValueError("projection seed does not match the endpoint rank")
    if len(seed.stefan_flux_scales_mol_m2_s) != layout.positive_area_face_count:
        raise ValueError("one projection Stefan scale is required per positive-area face")
    _encode_candidate(
        seed.candidate,
        chart,
        seed.stefan_flux_scales_mol_m2_s,
    )


def _coordinate_chart(before: ci.CutIntegratorState) -> _ProjectionChart:
    dry_config = before.transport.config.dry
    primitive_domain = dry_config.primitive_band
    return _ProjectionChart(
        dry_temperature=dry_config.conditioned_temperature_domain.solver_bounds_k,
        legacy_dry_y_hexane=(
            primitive_domain.y_hexane_bounds
            if isinstance(primitive_domain, ct.JointPrimitiveBand)
            else None
        ),
        dry_config=dry_config,
        dry_pores=tuple(
            replace(dry_config.pore, w_o=oil_label)
            for oil_label in before.transport.oil_fraction_labels
        ),
    )


def _encode_candidate(
    candidate: extinction.ExtinctionUnknowns,
    chart: _ProjectionChart,
    stefan_scales: Sequence[float],
) -> np.ndarray:
    if not (
        len(candidate.dry_temperatures_k)
        == len(candidate.dry_y_hexane)
        == len(chart.dry_pores)
    ):
        raise ValueError("projection dry primitives must align with material labels")
    values = (
        *(
            ci._encode_open(value, chart.dry_temperature)  # noqa: SLF001
            for value in candidate.dry_temperatures_k
        ),
        *(
            _encode_dry_y(temperature, composition, pore, chart)
            for temperature, composition, pore in zip(
                candidate.dry_temperatures_k,
                candidate.dry_y_hexane,
                chart.dry_pores,
            )
        ),
        *(
            value / scale
            for value, scale in zip(
                candidate.positive_area_total_stefan_fluxes_mol_m2_s,
                stefan_scales,
            )
        ),
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("encoded projection coordinates must remain finite")
    return np.asarray(values, dtype=float)


def _decode_candidate(
    coordinates: Sequence[float],
    chart: _ProjectionChart,
    stefan_scales: Sequence[float],
    event_duration_s: float,
) -> extinction.ExtinctionUnknowns:
    values = tuple(float(value) for value in coordinates)
    n = len(chart.dry_pores)
    if len(values) != 3 * n or not all(math.isfinite(value) for value in values):
        raise ValueError("projection coordinate vector is non-finite or wrong-rank")
    temperatures = tuple(
        ci._decode_open(value, chart.dry_temperature)  # noqa: SLF001
        for value in values[:n]
    )
    return extinction.ExtinctionUnknowns(
        dry_temperatures_k=temperatures,
        dry_y_hexane=tuple(
            _decode_dry_y(temperature, coordinate, pore, chart)
            for temperature, coordinate, pore in zip(
                temperatures,
                values[n : 2 * n],
                chart.dry_pores,
            )
        ),
        positive_area_total_stefan_fluxes_mol_m2_s=tuple(
            value * scale
            for value, scale in zip(values[2 * n :], stefan_scales)
        ),
        event_time_s=event_duration_s,
    )


def _encode_dry_y(
    temperature_k: float,
    y_hexane: float,
    pore: cp.CoupledPoreParams,
    chart: _ProjectionChart,
) -> float:
    if chart.legacy_dry_y_hexane is not None:
        return ci._encode_open(y_hexane, chart.legacy_dry_y_hexane)  # noqa: SLF001
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
    chart: _ProjectionChart,
) -> float:
    if chart.legacy_dry_y_hexane is not None:
        return ci._decode_open(coordinate, chart.legacy_dry_y_hexane)  # noqa: SLF001
    return cp.decode_gas_only_y(
        temperature_k,
        chart.dry_config.pressure_pa,
        coordinate,
        pore,
    )


def _bounded_fractions(
    candidate: extinction.ExtinctionUnknowns,
    chart: _ProjectionChart,
    topology: ExtinctionProjectionTopologyDiagnostics,
) -> tuple[float, ...]:
    pairs = tuple(
        (value, chart.dry_temperature) for value in candidate.dry_temperatures_k
    )
    fractions = [
        (value - bounds[0]) / (bounds[1] - bounds[0])
        for value, bounds in pairs
    ]
    if chart.legacy_dry_y_hexane is not None:
        fractions.extend(
            (value - chart.legacy_dry_y_hexane[0])
            / (chart.legacy_dry_y_hexane[1] - chart.legacy_dry_y_hexane[0])
            for value in candidate.dry_y_hexane
        )
    fractions.extend(topology.exact_composition_interval_fractions)
    for gaps in (
        topology.water_saturation_activity_gaps,
        topology.hexane_saturation_activity_gaps,
        topology.luikov_lower_activity_gaps,
    ):
        fractions.extend(gap / (1.0 + gap) for gap in gaps)
    fractions_tuple = tuple(fractions)
    if not all(0.0 < value < 1.0 for value in fractions_tuple):
        raise RuntimeError("accepted projection primitive touched a chart boundary")
    return fractions_tuple


def _topology_diagnostics(
    before: ci.CutIntegratorState,
    candidate: extinction.ExtinctionUnknowns,
    dry_cells: Sequence[cp.EquilibriumPoreState],
    chart: _ProjectionChart,
) -> ExtinctionProjectionTopologyDiagnostics:
    intervals = tuple(
        cp.gas_only_composition_interval(
            temperature,
            chart.dry_config.pressure_pa,
            pore,
        )
        for temperature, pore in zip(
            candidate.dry_temperatures_k,
            chart.dry_pores,
        )
    )
    fractions = tuple(
        (composition - interval.lower_y_hexane)
        / (interval.upper_y_hexane - interval.lower_y_hexane)
        for composition, interval in zip(candidate.dry_y_hexane, intervals)
    )
    return ExtinctionProjectionTopologyDiagnostics(
        uses_conditioned_composition_chart=(
            chart.legacy_dry_y_hexane is None
        ),
        oil_fraction_labels=before.transport.oil_fraction_labels,
        exact_composition_intervals=intervals,
        exact_composition_interval_fractions=fractions,
        water_saturation_activity_gaps=tuple(
            1.0 - cell.water_activity for cell in dry_cells
        ),
        hexane_saturation_activity_gaps=tuple(
            1.0 - cell.hexane_activity for cell in dry_cells
        ),
        luikov_lower_activity_gaps=tuple(
            cell.water_activity - sp.water_activity(pore.luikov.W_ref, pore.luikov)
            for cell, pore in zip(dry_cells, chart.dry_pores)
        ),
    )


def _residual_scales(
    before: ci.CutIntegratorState,
    initial: extinction.ExtinctionAssembly,
    event_duration_s: float,
) -> ExtinctionProjectionResidualScales:
    grid = before.transport.geometry.master_grid
    water_rates = tuple(
        area * flux.component.conserved_water_flux_mol_m2_s
        for area, flux in zip(grid.areas, initial.dry_face_fluxes)
    )
    hexane_rates = tuple(
        area * flux.component.conserved_hexane_flux_mol_m2_s
        for area, flux in zip(grid.areas, initial.dry_face_fluxes)
    )

    def component_scales(
        old: Sequence[float],
        new: Sequence[float],
        rates: Sequence[float],
    ) -> tuple[float, ...]:
        return tuple(
            max(
                abs(old[index]) / event_duration_s,
                abs(new[index]) / event_duration_s,
                abs(rates[index]),
                abs(rates[index + 1]),
                1.0e-300,
            )
            for index in range(grid.n)
        )

    energy_scales = []
    dry_capacity_audits = []
    for index, cell in enumerate(initial.dry_cells):
        pore = cp.CoupledPoreParams(
            **{
                **before.transport.config.dry.pore.__dict__,
                "w_o": before.transport.oil_fraction_labels[index],
            }
        )
        capacity, audit = ci._dry_energy_capacity_with_audit(  # noqa: SLF001
            cell.temperature_k,
            cell.y_hexane,
            pore,
            before.transport.config,
        )
        dry_capacity_audits.append(audit)
        energy_scales.append(
            grid.volumes[index]
            * capacity
            * max(abs(cell.temperature_k), 1.0)
            / event_duration_s
        )
    return ExtinctionProjectionResidualScales(
        dry_water_mol_s=component_scales(
            initial.old_inventory.water_cell_mol,
            initial.current_inventory.water_cell_mol,
            water_rates,
        ),
        dry_hexane_mol_s=component_scales(
            initial.old_inventory.hexane_cell_mol,
            initial.current_inventory.hexane_cell_mol,
            hexane_rates,
        ),
        dry_energy_w=tuple(energy_scales),
        dry_energy_capacity_stencil_audits=tuple(dry_capacity_audits),
    )


def _accept_projection(
    before: ci.CutIntegratorState,
    boundary: ct.PoreBoundary,
    seed: ExtinctionProjectionSeed,
    assembly: extinction.ExtinctionAssembly,
    scales: ExtinctionProjectionResidualScales,
    signed_scaled: tuple[float, ...],
    evaluations: int,
    nonlinear_message: str,
    condition_proxy: float,
    minimum_fraction: float,
    topology_diagnostics: ExtinctionProjectionTopologyDiagnostics,
    structural: _ProjectionSparsity | None,
) -> ExtinctionProjectionStep:
    duration = assembly.candidate.event_time_s
    changes = assembly.inventory_changes
    surface_area = assembly.endpoint_geometry.master_grid.areas[-1]
    surface = assembly.dry_face_fluxes[-1]
    water_out = (
        duration
        * surface_area
        * surface.component.conserved_water_flux_mol_m2_s
    )
    hexane_out = (
        duration
        * surface_area
        * surface.component.conserved_hexane_flux_mol_m2_s
    )
    energy_out = duration * surface_area * surface.energy.total_energy_flux_w_m2
    water_step = math.fsum((*changes.water_cell_mol, water_out))
    hexane_step = math.fsum((*changes.hexane_cell_mol, hexane_out))
    energy_step = math.fsum((*changes.energy_cell_j, energy_out))
    water_scale = max(
        assembly.old_inventory.total_water_mol,
        assembly.current_inventory.total_water_mol,
        abs(water_out),
        1.0e-300,
    )
    hexane_scale = max(
        assembly.old_inventory.total_hexane_mol,
        assembly.current_inventory.total_hexane_mol,
        abs(hexane_out),
        1.0e-300,
    )
    endpoint_capacity = math.fsum(scales.dry_energy_w) * duration
    energy_scale = max(
        before.reference_capacity_energy_scale_j,
        endpoint_capacity,
        1.0e-300,
    )

    water_high, water_boundary_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_water_out_mol,
        before.boundary_water_compensation_mol,
        water_out,
    )
    hexane_high, hexane_boundary_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_hexane_out_mol,
        before.boundary_hexane_compensation_mol,
        hexane_out,
    )
    energy_high, energy_boundary_comp = ci._compensated_add_scalar(  # noqa: SLF001
        before.cumulative_boundary_energy_out_j,
        before.boundary_energy_compensation_j,
        energy_out,
    )
    water_changes, water_comp = ci._compensated_add_vectors(  # noqa: SLF001
        before.cumulative_material_water_change_cell_mol,
        before.water_change_compensation_cell_mol,
        changes.water_cell_mol,
    )
    hexane_changes, hexane_comp = ci._compensated_add_vectors(  # noqa: SLF001
        before.cumulative_material_hexane_change_cell_mol,
        before.hexane_change_compensation_cell_mol,
        changes.hexane_cell_mol,
    )
    energy_changes, energy_comp = ci._compensated_add_vectors(  # noqa: SLF001
        before.cumulative_material_energy_change_cell_j,
        before.energy_change_compensation_cell_j,
        changes.energy_cell_j,
    )

    def corrected(values: Sequence[float], compensation: Sequence[float]):
        return tuple(
            math.fsum((value, correction))
            for value, correction in zip(values, compensation)
        )

    cumulative_water_out = math.fsum((water_high, water_boundary_comp))
    cumulative_hexane_out = math.fsum((hexane_high, hexane_boundary_comp))
    cumulative_energy_out = math.fsum((energy_high, energy_boundary_comp))
    cumulative_water = math.fsum(
        (*corrected(water_changes, water_comp), cumulative_water_out)
    )
    cumulative_hexane = math.fsum(
        (*corrected(hexane_changes, hexane_comp), cumulative_hexane_out)
    )
    cumulative_energy = math.fsum(
        (*corrected(energy_changes, energy_comp), cumulative_energy_out)
    )
    cumulative_water_scale = max(
        before.reference_inventory.total_water_mol,
        assembly.current_inventory.total_water_mol,
        abs(cumulative_water_out),
        before.cumulative_absolute_water_transfer_mol + abs(water_out),
        1.0e-300,
    )
    cumulative_hexane_scale = max(
        before.reference_inventory.total_hexane_mol,
        assembly.current_inventory.total_hexane_mol,
        abs(cumulative_hexane_out),
        before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out),
        1.0e-300,
    )
    cumulative_energy_scale = max(
        before.reference_capacity_energy_scale_j,
        endpoint_capacity,
        before.cumulative_absolute_energy_transfer_j + abs(energy_out),
        1.0e-300,
    )
    after = ExtinctionCommittedState(
        grid=assembly.endpoint_geometry.master_grid,
        config=before.transport.config.dry,
        controls=before.controls,
        time_s=before.transport.time_s + duration,
        temperatures_k=assembly.candidate.dry_temperatures_k,
        y_hexane=assembly.candidate.dry_y_hexane,
        historical_hexane_loadings=before.transport.historical_hexane_loadings,
        oil_fraction_labels=before.transport.oil_fraction_labels,
        reference_inventory=before.reference_inventory,
        reference_capacity_energy_scale_j=before.reference_capacity_energy_scale_j,
        last_positive_area_total_stefan_fluxes_mol_m2_s=(
            assembly.candidate.positive_area_total_stefan_fluxes_mol_m2_s
        ),
        cumulative_boundary_water_out_mol=water_high,
        cumulative_boundary_hexane_out_mol=hexane_high,
        cumulative_boundary_energy_out_j=energy_high,
        boundary_water_compensation_mol=water_boundary_comp,
        boundary_hexane_compensation_mol=hexane_boundary_comp,
        boundary_energy_compensation_j=energy_boundary_comp,
        cumulative_absolute_water_transfer_mol=(
            before.cumulative_absolute_water_transfer_mol + abs(water_out)
        ),
        cumulative_absolute_hexane_transfer_mol=(
            before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out)
        ),
        cumulative_absolute_energy_transfer_j=(
            before.cumulative_absolute_energy_transfer_j + abs(energy_out)
        ),
        cumulative_material_water_change_cell_mol=water_changes,
        cumulative_material_hexane_change_cell_mol=hexane_changes,
        cumulative_material_energy_change_cell_j=energy_changes,
        water_change_compensation_cell_mol=water_comp,
        hexane_change_compensation_cell_mol=hexane_comp,
        energy_change_compensation_cell_j=energy_comp,
        accepted_steps=before.accepted_steps + 1,
        cumulative_nonlinear_evaluations=(
            before.cumulative_nonlinear_evaluations + evaluations
        ),
    )
    n = assembly.layout.dry_cell_count
    ledger = ExtinctionProjectionLedger(
        event_duration_s=duration,
        nonlinear_evaluations=evaluations,
        nonlinear_message=nonlinear_message,
        residual_scales=scales,
        scaled_residuals=signed_scaled,
        maximum_scaled_residual=max(abs(value) for value in signed_scaled),
        condition_proxy=condition_proxy,
        minimum_fractional_distance_to_bound=minimum_fraction,
        topology_diagnostics=topology_diagnostics,
        used_sparse_jacobian=structural is not None,
        jacobian_structural_nnz=(
            (3 * n) ** 2 if structural is None else structural.nnz
        ),
        jacobian_structural_color_count=(
            3 * n if structural is None else structural.color_count
        ),
        water_step_residual_mol=water_step,
        hexane_step_residual_mol=hexane_step,
        energy_step_residual_j=energy_step,
        normalized_water_step_residual=abs(water_step) / water_scale,
        normalized_hexane_step_residual=abs(hexane_step) / hexane_scale,
        normalized_energy_step_residual=abs(energy_step) / energy_scale,
        water_cumulative_residual_mol=cumulative_water,
        hexane_cumulative_residual_mol=cumulative_hexane,
        energy_cumulative_residual_j=cumulative_energy,
        normalized_water_cumulative_residual=(
            abs(cumulative_water) / cumulative_water_scale
        ),
        normalized_hexane_cumulative_residual=(
            abs(cumulative_hexane) / cumulative_hexane_scale
        ),
        normalized_energy_cumulative_residual=(
            abs(cumulative_energy) / cumulative_energy_scale
        ),
        boundary_water_out_mol=water_out,
        boundary_hexane_out_mol=hexane_out,
        boundary_energy_out_j=energy_out,
        accepted=True,
    )
    return ExtinctionProjectionStep(before, after, assembly, boundary, seed, ledger)


def _projection_sparsity(cell_count: int) -> _ProjectionSparsity:
    """Return the bordered-free nearest-neighbour ``3*N`` endpoint pattern."""

    if cell_count < 1:
        raise ValueError("projection sparsity needs at least one cell")
    n = cell_count
    row_columns: list[tuple[int, ...]] = []
    for equation in range(3):
        for cell in range(n):
            columns: set[int] = set()
            for neighbor in range(max(0, cell - 1), min(n, cell + 2)):
                columns.add(neighbor)
                columns.add(n + neighbor)
            if cell > 0:
                columns.add(2 * n + cell - 1)
            columns.add(2 * n + cell)
            row_columns.append(tuple(sorted(columns)))

    # Both construction and coloring stay O(N): every row has at most eight
    # columns, and the column-intersection graph therefore has bounded degree.
    indptr = np.empty(3 * n + 1, dtype=np.int64)
    indptr[0] = 0
    for row, columns in enumerate(row_columns, start=1):
        indptr[row] = indptr[row - 1] + len(columns)
    indices = np.fromiter(
        (column for columns in row_columns for column in columns),
        dtype=np.int64,
        count=int(indptr[-1]),
    )
    matrix = sparse.csr_matrix(
        (np.ones(len(indices), dtype=bool), indices, indptr),
        shape=(3 * n, 3 * n),
    )
    colors = _bounded_degree_column_coloring(3 * n, row_columns)
    return _ProjectionSparsity(matrix, colors)


def _bounded_degree_column_coloring(
    column_count: int,
    row_columns: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    """Greedily color a bounded-stencil column graph in linear work."""

    neighbours = [set() for _ in range(column_count)]
    for columns in row_columns:
        for position, column in enumerate(columns):
            neighbours[column].update(columns[:position])
            neighbours[column].update(columns[position + 1 :])
    colors = [-1] * column_count
    for column, adjacent in enumerate(neighbours):
        forbidden = {
            colors[other] for other in adjacent if other < column
        }
        color = 0
        while color in forbidden:
            color += 1
        colors[column] = color
    return tuple(colors)


__all__ = [
    "ExtinctionCommittedState",
    "ExtinctionProjectionCompatibilityError",
    "ExtinctionProjectionLedger",
    "ExtinctionProjectionResidualScales",
    "ExtinctionProjectionSeed",
    "ExtinctionProjectionStep",
    "ExtinctionProjectionStepError",
    "ExtinctionProjectionTopologyDiagnostics",
    "solve_fixed_time_extinction_projection",
]
