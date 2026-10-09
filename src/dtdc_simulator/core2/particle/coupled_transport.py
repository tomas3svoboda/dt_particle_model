r"""Bounded numerical transport slice for the fully dry stationary sphere.

This module closes one deliberately narrow part of Gate 1g.  Each material
shell owns two primitives, temperature and pore-gas n-hexane mole fraction.
At prescribed bed pressure the summed component equation recovers the total
Stefan flux algebraically; the remaining n-hexane equation and the
common-datum energy equation form a square ``2*N`` backward-Euler system.

The transport branch is explicitly ``COMMON_FACE_T_NEGLIGIBLE_SORET``.  At
each face both compositions are evaluated at the same symmetric face
temperature before the frozen Maxwell--Stefan and retained-water laws are
called.  Thus temperature gradients drive Fourier heat conduction, but no
unidentified Soret term is smuggled into either mass-transfer potential.

An optional finite prescribed-pressure transition evaluates old storage at
``P_n`` and the complete backward-Euler candidate, interface fluxes, and
algebraic Stefan recovery at ``P_(n+1)``.  It is a finite conservative step,
never an instantaneous state projection.  Only pressure may differ between
the two immutable configurations, and the target configuration is committed
only after every nonlinear and conservation check succeeds.

All coefficient values and the jointly admissible primitive band come from
the caller.  There is no water reservoir, clipping, default transport point,
front, ALE term, tolerance relaxation, or claim of physical Gate-1g
qualification.  The finite-pressure branch remains a numerical foundation
until a finite-film and bed-hydraulic pressure authority is connected.  The
immutable input state is the exact rollback state whenever a step fails.
"""

from __future__ import annotations

import math
import operator
from collections.abc import Callable
from dataclasses import dataclass, replace as dataclass_replace
from enum import Enum
from functools import lru_cache
from typing import Sequence

import numpy as np
from scipy import optimize
from scipy.special import expit

from dtdc_simulator.core2.particle import bed_film as bf
from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import external_film_reduced as efr
from dtdc_simulator.core2.particle import transport_coefficients as tc
from dtdc_simulator.core2.particle.grid import SphericalGrid
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import water as wa


class DryTransportBranch(str, Enum):
    """Named approximation used by this bounded numerical slice."""

    COMMON_FACE_T_NEGLIGIBLE_SORET = "common_face_temperature_negligible_soret"


class CoupledTransportStepError(RuntimeError):
    """A rejected step carrying the exact immutable rollback state."""

    def __init__(
        self,
        message: str,
        rollback_state: "FullyDryTransportState",
        *,
        nonlinear_evaluations: int,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.nonlinear_evaluations = nonlinear_evaluations


@dataclass(frozen=True)
class ThermalConductivitySelection:
    """One explicit dry-particle conductivity case, with provenance."""

    value_w_m_k: float
    label: str

    def __post_init__(self) -> None:
        if not math.isfinite(self.value_w_m_k) or self.value_w_m_k <= 0.0:
            raise ValueError("thermal conductivity must be positive and finite")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("thermal conductivity needs a non-empty provenance label")

    @property
    def physically_qualifying(self) -> bool:
        return False


@dataclass(frozen=True)
class JointPrimitiveBand:
    """Caller-declared rectangular subset of one gas-only topology.

    This is a nonlinear-solver feasibility band, not a clip.  Every trial is
    represented by smooth bounded coordinates and an inadmissible evaluation
    rejects the whole step.  The provenance label must say why the band is
    appropriate for the particular campaign.
    """

    temperature_bounds_k: tuple[float, float]
    y_hexane_bounds: tuple[float, float]
    label: str

    def __post_init__(self) -> None:
        t_lo, t_hi = self.temperature_bounds_k
        y_lo, y_hi = self.y_hexane_bounds
        if not all(math.isfinite(value) for value in (t_lo, t_hi, y_lo, y_hi)):
            raise ValueError("primitive-band bounds must be finite")
        if not 0.0 < t_lo < t_hi:
            raise ValueError("temperature bounds must be positive and ordered")
        if not 0.0 < y_lo < y_hi < 1.0:
            raise ValueError("n-hexane mole-fraction bounds must lie inside (0,1)")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("primitive band needs a non-empty provenance label")

    @property
    def uses_exact_phase_conditioned_composition(self) -> bool:
        return False


@dataclass(frozen=True)
class GasOnlyPrimitiveDomain:
    """Requested temperature authority with exact local composition charts.

    Unlike :class:`JointPrimitiveBand`, this domain does not declare one
    rectangular composition interval.  At every decoded temperature the
    admissible ``y_h`` interval comes from the frozen pressure-conditioned
    gas-only phase/Luikov inequalities in :mod:`coupled_pore`.
    """

    temperature_bounds_k: tuple[float, float]
    label: str

    def __post_init__(self) -> None:
        lower, upper = self.temperature_bounds_k
        if not all(math.isfinite(value) for value in (lower, upper)):
            raise ValueError("gas-only temperature bounds must be finite")
        if not 0.0 < lower < upper:
            raise ValueError("gas-only temperature bounds must be positive and ordered")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("gas-only primitive domain needs a non-empty provenance label")

    @property
    def uses_exact_phase_conditioned_composition(self) -> bool:
        return True


@dataclass(frozen=True)
class PressureConditionedTemperatureDomain:
    """Connected representable gas-only temperature interval at fixed pressure."""

    requested_bounds_k: tuple[float, float]
    solver_bounds_k: tuple[float, float]
    pressure_pa: float
    lower_is_phase_boundary: bool
    lower_boundary_bracket_k: tuple[float, float]
    physically_qualifying: bool = False


@dataclass(frozen=True)
class FullyDryTransportConfig:
    """Explicit coefficients and numerical contract for the predictive slice."""

    pressure_pa: float
    binary_diffusivity: tc.EffectiveBinaryDiffusivitySelection
    retained_water_mobility: cp.MobilitySelection
    thermal_conductivity: ThermalConductivitySelection
    primitive_band: JointPrimitiveBand | GasOnlyPrimitiveDomain
    pore: cp.CoupledPoreParams = cp.CoupledPoreParams()
    nonlinear_residual_tolerance: float = 1.0e-10
    ledger_tolerance: float = 1.0e-10
    nonlinear_step_tolerance: float = 2.0e-11
    maximum_function_evaluations: int = 400

    def __post_init__(self) -> None:
        if not math.isfinite(self.pressure_pa) or self.pressure_pa <= 0.0:
            raise ValueError("pressure must be positive and finite")
        if not isinstance(self.binary_diffusivity, tc.EffectiveBinaryDiffusivitySelection):
            raise TypeError("binary diffusivity must be an explicit selection")
        if not isinstance(self.retained_water_mobility, cp.MobilitySelection):
            raise TypeError("retained-water mobility must be an explicit selection")
        if not isinstance(self.thermal_conductivity, ThermalConductivitySelection):
            raise TypeError("thermal conductivity must be an explicit selection")
        if not isinstance(
            self.primitive_band,
            (JointPrimitiveBand, GasOnlyPrimitiveDomain),
        ):
            raise TypeError(
                "a rectangular numerical band or exact gas-only primitive domain is required"
            )
        tolerances = (
            self.nonlinear_residual_tolerance,
            self.ledger_tolerance,
            self.nonlinear_step_tolerance,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in tolerances):
            raise ValueError("solver tolerances must be positive and finite")
        if self.nonlinear_residual_tolerance > 1.0e-10:
            raise ValueError("accepted nonlinear residual tolerance cannot exceed 1e-10")
        if self.ledger_tolerance > 1.0e-10:
            raise ValueError("accepted ledger tolerance cannot exceed 1e-10")
        try:
            max_evaluations = operator.index(self.maximum_function_evaluations)
        except TypeError as exc:
            raise ValueError("maximum function evaluations must be an integer") from exc
        if isinstance(self.maximum_function_evaluations, bool) or max_evaluations < 1:
            raise ValueError("maximum function evaluations must be positive")
        _validate_declared_band(self)

    @property
    def branch(self) -> DryTransportBranch:
        return DryTransportBranch.COMMON_FACE_T_NEGLIGIBLE_SORET

    @property
    def physically_qualifying(self) -> bool:
        return False

    @property
    def conditioned_temperature_domain(self) -> PressureConditionedTemperatureDomain:
        return _pressure_conditioned_temperature_domain(self)


@dataclass(frozen=True)
class DirichletPoreBoundary:
    """Numerical Dirichlet oracle for external temperature and composition.

    It is not the frozen external-film active set and is consequently never
    physical Gate-1g evidence.
    """

    temperature_k: float
    y_hexane: float
    label: str
    physically_qualifying: bool = False

    def __post_init__(self) -> None:
        if not math.isfinite(self.temperature_k) or self.temperature_k <= 0.0:
            raise ValueError("boundary temperature must be positive and finite")
        if not math.isfinite(self.y_hexane) or not 0.0 < self.y_hexane < 1.0:
            raise ValueError("boundary n-hexane mole fraction must lie in (0,1)")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("boundary needs a non-empty provenance label")


@dataclass(frozen=True)
class ReducedFilmPoreBoundary:
    """Bounded fast-mass/exact-A finite-heat engineering boundary.

    ``temperature_k`` and ``y_hexane`` are bulk-bed values.  The surface
    composition uses the audited fast-film asymptote ``y_s=y_b``.  The
    surface temperature is resolved by the scalar solid/film heat closure
    using the exact component-flux Ackermann factor.

    The mass-fraction tolerance is the caller-owned production guard.  The
    Ackermann departure threshold is retained only to report how badly a
    counterfactual ``A=1`` approximation would perform.  Production ``A`` uses
    exact PHY-053 component partial-enthalpy secants between surface and bulk
    at this fixed ``(P,y,k_wh)``.  The two constant heat capacities are kept
    only as the explicitly named legacy engineering-regression counterfactual.
    Dynamic DTDC validity-band and production Ackermann-domain failures remain
    fail-closed.
    """

    temperature_k: float
    y_hexane: float
    correlated_pair: bf.CorrelatedBedFilmPair
    stress_multiplier: float
    liquid_saturation: float
    forcing_timescale_s: float
    absolute_mass_fraction_departure_threshold: float
    absolute_ackermann_factor_departure_threshold: float
    legacy_engineering_regression_water_heat_capacity_j_kg_k: float
    legacy_engineering_regression_hexane_heat_capacity_j_kg_k: float
    binary_gas_interaction_k_wh: float
    label: str
    fidelity_band: efr.FidelityBand = efr.FidelityBand.NOMINAL
    drainage_regime: efr.DrainageRegime = efr.DrainageRegime.PRIMARY_DRAINAGE
    physically_qualifying: bool = False

    def __post_init__(self) -> None:
        if not math.isfinite(self.temperature_k) or self.temperature_k <= 0.0:
            raise ValueError("reduced-film bulk temperature must be positive and finite")
        if not math.isfinite(self.y_hexane) or not 0.0 < self.y_hexane < 1.0:
            raise ValueError("reduced-film bulk n-hexane mole fraction must lie in (0,1)")
        if not isinstance(self.correlated_pair, bf.CorrelatedBedFilmPair):
            raise TypeError("reduced boundary requires one CorrelatedBedFilmPair")
        if (
            isinstance(self.stress_multiplier, bool)
            or not math.isfinite(self.stress_multiplier)
            or float(self.stress_multiplier) not in efr.ALLOWED_STRESS_MULTIPLIERS
        ):
            raise ValueError("reduced boundary common stress multiplier must be {0.5,0.6264,1,2}")
        if not math.isfinite(self.liquid_saturation) or not 0.0 <= self.liquid_saturation <= 0.8:
            raise ValueError("reduced boundary liquid saturation must lie in [0,0.8]")
        if not math.isfinite(self.forcing_timescale_s) or self.forcing_timescale_s <= 0.0:
            raise ValueError("reduced boundary forcing timescale must be positive")
        thresholds = (
            self.absolute_mass_fraction_departure_threshold,
            self.absolute_ackermann_factor_departure_threshold,
        )
        if not all(math.isfinite(value) and 0.0 < value <= 1.0 for value in thresholds):
            raise ValueError("reduced boundary audit thresholds must lie in (0,1]")
        capacities = (
            self.legacy_engineering_regression_water_heat_capacity_j_kg_k,
            self.legacy_engineering_regression_hexane_heat_capacity_j_kg_k,
        )
        if not all(math.isfinite(value) and value > 0.0 for value in capacities):
            raise ValueError("legacy engineering-regression heat capacities must be positive")
        if not bg.K_WH_MIN <= self.binary_gas_interaction_k_wh <= bg.K_WH_MAX:
            raise ValueError("reduced boundary k_wh is outside the frozen PHY-053 interval")
        if self.fidelity_band not in (efr.FidelityBand.NOMINAL, efr.FidelityBand.BOUNDED_ANOMALY):
            raise ValueError("reduced boundary fidelity band must be explicit")
        if self.drainage_regime is not efr.DrainageRegime.PRIMARY_DRAINAGE:
            raise ValueError("reduced boundary is restricted to primary drainage")
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("reduced boundary needs a non-empty provenance label")


PoreBoundary = DirichletPoreBoundary | ReducedFilmPoreBoundary


@dataclass(frozen=True)
class ReducedFilmSurfaceState:
    """A=1 reference geometry used to bracket the exact surface closure."""

    boundary: ReducedFilmPoreBoundary
    coefficients: efr.ReducedFilmCoefficients
    cell_temperature_k: float
    surface_temperature_k: float
    surface_y_hexane: float
    surface_hexane_mass_fraction: float
    outer_half_cell_distance_m: float
    outer_half_cell_heat_resistance_m2_k_w: float
    external_heat_resistance_m2_k_w: float
    series_heat_transfer_w_m2_k: float
    production_conductive_heat_flux_w_m2: float
    solid_temperature_gradient_k_m: float


@dataclass(frozen=True)
class ReducedFilmSurfaceFluxEvaluation:
    """Component fluxes and caller payload evaluated at one surface trial.

    ``non_fourier_reduced_heat_flux_w_m2`` is an explicitly selected internal
    reciprocal/heat-of-transport contribution.  It defaults to exact zero so
    every historical caller retains the original Fourier/Ackermann scalar
    closure bit for bit.  A nonzero value is legal only when the scalar solver
    is also given an admissible surface-temperature bracket; this prevents a
    caller from appending an unclosed heat term after the surface root.
    """

    component: cp.ComponentMolarFluxes
    payload: object | None = None
    non_fourier_reduced_heat_flux_w_m2: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.component, cp.ComponentMolarFluxes):
            raise TypeError("surface closure requires ComponentMolarFluxes")
        if not math.isfinite(self.non_fourier_reduced_heat_flux_w_m2):
            raise ValueError("surface non-Fourier reduced heat must be finite")

    @property
    def component_mass_fluxes(self) -> efr.ComponentMassFluxes:
        return efr.ComponentMassFluxes(
            water_kg_m2_s=(self.component.conserved_water_flux_mol_m2_s * wa.M),
            hexane_kg_m2_s=(self.component.conserved_hexane_flux_mol_m2_s * hx.M),
        )


@dataclass(frozen=True)
class ReducedFilmSurfaceClosure:
    """Exact scalar Ackermann/solid/film surface closure and final flux."""

    reference: ReducedFilmSurfaceState
    coefficients: efr.ReducedFilmCoefficients
    surface_temperature_k: float
    surface_y_hexane: float
    surface_hexane_mass_fraction: float
    exact_partial_enthalpy_secant_heat_capacities: efr.PartialEnthalpySecantHeatCapacities
    ackermann: efr.AckermannHeatTransfer
    surface_flux_evaluation: ReducedFilmSurfaceFluxEvaluation
    external_ackermann_heat_resistance_m2_k_w: float
    series_heat_transfer_w_m2_k: float
    conductive_heat_flux_outward_w_m2: float
    fourier_heat_flux_outward_w_m2: float
    reciprocal_heat_of_transport_flux_outward_w_m2: float
    reduced_heat_flux_outward_w_m2: float
    external_ackermann_heat_flux_outward_w_m2: float
    local_reduced_heat_closure_residual_w_m2: float
    solid_temperature_gradient_k_m: float
    heat_closure_residual_w_m2: float
    normalized_heat_closure_residual: float
    heat_closure_roundoff_bound_w_m2: float
    temperature_bracket_k: tuple[float, float]
    scalar_function_evaluations: int
    monotone_unique_root_diagnostic_passed: bool
    exact_equal_temperature_branch: bool
    exact_fourier_series_identity_applicable: bool
    used_generalized_non_fourier_surface_closure: bool


@dataclass(frozen=True)
class SurfaceFilmAudit:
    """Exact-A heat evidence plus the one remaining fast-mass reduction."""

    common_provenance: str
    stress_multiplier: float
    heat_transfer_w_m2_k: float
    binary_mass_transfer_m_s: float
    coefficient_covariance_preserved: bool
    bulk_temperature_k: float
    surface_temperature_k: float
    pressure_pa: float
    binary_gas_interaction_k_wh: float
    bulk_y_hexane: float
    assumed_surface_y_hexane: float
    assumed_surface_hexane_mass_fraction: float
    outer_half_cell_heat_resistance_m2_k_w: float
    external_heat_resistance_m2_k_w: float
    series_heat_transfer_w_m2_k: float
    production_conductive_heat_flux_w_m2: float
    particle_fourier_heat_flux_w_m2: float
    particle_reciprocal_heat_of_transport_flux_w_m2: float
    particle_reduced_heat_flux_w_m2: float
    external_ackermann_heat_flux_w_m2: float
    local_reduced_heat_closure_residual_w_m2: float
    full_ackermann_conductive_heat_flux_w_m2: float
    full_ackermann_surface_temperature_k: float
    linear_heat_absolute_error_w_m2: float
    unity_ackermann_conductive_heat_flux_w_m2: float
    unity_ackermann_surface_temperature_k: float
    unity_ackermann_series_heat_transfer_w_m2_k: float
    mass_departure: efr.SuppliedFluxCompositionDeparture
    ackermann: efr.AckermannHeatTransfer
    exact_water_partial_enthalpy_secant_heat_capacity_j_kg_k: float
    exact_hexane_partial_enthalpy_secant_heat_capacity_j_kg_k: float
    exact_secant_equal_temperature_derivative_limit_used: bool
    exact_secant_derivative_stencil_half_width_k: float | None
    legacy_engineering_regression_water_heat_capacity_j_kg_k: float
    legacy_engineering_regression_hexane_heat_capacity_j_kg_k: float
    legacy_engineering_regression_ackermann_factor: float
    legacy_engineering_regression_ackermann_relative_error: float
    legacy_engineering_regression_series_heat_transfer_w_m2_k: float
    legacy_engineering_regression_series_conductance_relative_error: float
    legacy_engineering_regression_series_conductance_relative_error_threshold: float
    legacy_engineering_regression_series_conductance_guard_passed: bool
    ackermann_absolute_factor_departure: float
    caller_ackermann_factor_departure_threshold: float
    water_gas_partial_enthalpy_j_mol: float
    hexane_gas_partial_enthalpy_j_mol: float
    production_common_datum_energy_outward_w_m2: float
    particle_inventory_energy_source_w_m2: float
    bed_inventory_energy_source_w_m2: float
    closed_pair_energy_residual_w_m2: float
    fast_mass_guard_passed: bool
    linear_heat_guard_passed: bool
    validity_band_checked: bool
    reduction_thresholds_enforced: bool
    heat_closure_residual_w_m2: float
    normalized_heat_closure_residual: float
    heat_closure_roundoff_bound_w_m2: float
    surface_temperature_bracket_k: tuple[float, float]
    scalar_closure_function_evaluations: int
    exact_equal_temperature_branch: bool
    monotone_unique_root_diagnostic_passed: bool
    exact_fourier_series_identity_applicable: bool
    used_generalized_non_fourier_surface_closure: bool
    production_used_exact_ackermann_factor: bool = True
    production_used_exact_phy053_partial_enthalpy_secant_heat_capacities: bool = True
    legacy_engineering_regression_heat_capacities_are_counterfactual_only: bool = True
    unity_ackermann_is_counterfactual_only: bool = True
    production_used_exact_unity_ackermann_factor: bool = False
    production_used_fast_mass_surface_equals_bulk: bool = True
    conserved_water_crosses_external_film_as_gas_once: bool = True


def _hexane_mass_fraction_from_mole_fraction(y_hexane: float) -> float:
    hexane_mass = y_hexane * hx.M
    water_mass = (1.0 - y_hexane) * wa.M
    return hexane_mass / (hexane_mass + water_mass)


def prepare_reduced_film_surface_state(
    boundary: ReducedFilmPoreBoundary,
    *,
    cell_temperature_k: float,
    pressure_pa: float,
    particle_radius_m: float,
    outer_half_cell_distance_m: float,
    particle_thermal_conductivity_w_m_k: float,
    binary_gas_interaction_k_wh: float,
) -> ReducedFilmSurfaceState:
    """Build the A=1 reference geometry for the exact scalar closure."""

    if not isinstance(boundary, ReducedFilmPoreBoundary):
        raise TypeError("a ReducedFilmPoreBoundary is required")
    values = (
        cell_temperature_k,
        pressure_pa,
        particle_radius_m,
        outer_half_cell_distance_m,
        particle_thermal_conductivity_w_m_k,
        binary_gas_interaction_k_wh,
    )
    if not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError("reduced-film surface geometry and state must be positive")
    if binary_gas_interaction_k_wh != boundary.binary_gas_interaction_k_wh:
        raise efr.ReducedFilmValidityError(
            "reduced-film Ackermann caloric k_wh differs from the pore PHY-053 state"
        )

    scaled_heat = float(boundary.stress_multiplier) * boundary.correlated_pair.heat_transfer_w_m2_k
    if not math.isfinite(scaled_heat) or scaled_heat <= 0.0:
        raise efr.ReducedFilmValidityError("reduced boundary requires positive h_Q")
    solid_resistance = outer_half_cell_distance_m / particle_thermal_conductivity_w_m_k
    film_resistance = 1.0 / scaled_heat
    series_coefficient = 1.0 / (solid_resistance + film_resistance)
    conductive_outward = series_coefficient * (cell_temperature_k - boundary.temperature_k)
    surface_temperature = math.fsum((cell_temperature_k, -conductive_outward * solid_resistance))
    solid_gradient = (surface_temperature - cell_temperature_k) / outer_half_cell_distance_m

    operating_point = efr.ReducedFilmOperatingPoint(
        bulk_temperature_k=boundary.temperature_k,
        interface_temperature_k=surface_temperature,
        pressure_pa=pressure_pa,
        particle_radius_m=particle_radius_m,
        liquid_saturation=boundary.liquid_saturation,
        forcing_timescale_s=boundary.forcing_timescale_s,
        drainage_regime=boundary.drainage_regime,
        fidelity_band=boundary.fidelity_band,
    )
    coefficients = efr.adapt_correlated_bed_film(
        boundary.correlated_pair,
        operating_point,
        stress_multiplier=float(boundary.stress_multiplier),
    )
    if coefficients.heat_transfer_w_m2_k != scaled_heat:
        raise RuntimeError("series heat boundary separated the correlated coefficient pair")
    return ReducedFilmSurfaceState(
        boundary=boundary,
        coefficients=coefficients,
        cell_temperature_k=cell_temperature_k,
        surface_temperature_k=surface_temperature,
        surface_y_hexane=boundary.y_hexane,
        surface_hexane_mass_fraction=_hexane_mass_fraction_from_mole_fraction(boundary.y_hexane),
        outer_half_cell_distance_m=outer_half_cell_distance_m,
        outer_half_cell_heat_resistance_m2_k_w=solid_resistance,
        external_heat_resistance_m2_k_w=film_resistance,
        series_heat_transfer_w_m2_k=series_coefficient,
        production_conductive_heat_flux_w_m2=conductive_outward,
        solid_temperature_gradient_k_m=solid_gradient,
    )


def _reduced_film_coefficients_at_surface_temperature(
    reference: ReducedFilmSurfaceState,
    surface_temperature_k: float,
) -> efr.ReducedFilmCoefficients:
    boundary = reference.boundary
    operating_point = efr.ReducedFilmOperatingPoint(
        bulk_temperature_k=boundary.temperature_k,
        interface_temperature_k=surface_temperature_k,
        pressure_pa=reference.coefficients.operating_point.pressure_pa,
        particle_radius_m=reference.coefficients.operating_point.particle_radius_m,
        liquid_saturation=boundary.liquid_saturation,
        forcing_timescale_s=boundary.forcing_timescale_s,
        drainage_regime=boundary.drainage_regime,
        fidelity_band=boundary.fidelity_band,
    )
    return efr.adapt_correlated_bed_film(
        boundary.correlated_pair,
        operating_point,
        stress_multiplier=float(boundary.stress_multiplier),
    )


def solve_reduced_film_surface_closure(
    reference: ReducedFilmSurfaceState,
    component_flux_evaluator: Callable[[float, float], ReducedFilmSurfaceFluxEvaluation],
    *,
    admissible_surface_temperature_bounds_k: tuple[float, float] | None = None,
) -> ReducedFilmSurfaceClosure:
    r"""Close ``T_s``, component fluxes, Ackermann ``A``, and heat exactly.

    The historical zero-cross-effect scalar residual is

    ``F(T_s)=k_p*(T_cell-T_s)/delta-h_Q*A(N_w,N_h)*(T_s-T_b)``.

    Positive finite resistances and ``A`` bracket a root between the cell and
    bulk temperatures for either heat-flow sign.  Five deterministic bracket
    samples must be strictly decreasing before safeguarded Brent iteration;
    topology, Ackermann-domain, nonfinite, and nonunique diagnostics all fail
    closed.  The component-flux evaluation at the returned root is retained
    so mass and energy use the identical surface state.

    A caller that selects a non-Fourier reciprocal heat contribution must
    supply its complete admissible surface-temperature interval.  That branch
    instead closes

    ``F(T_s)=q_F(T_s)+q_D(T_s)-h_Q*A(N_w,N_h)*(T_s-T_b)``.

    It scans the caller-owned interval deterministically, requires a strictly
    decreasing sampled residual and one bracketed root, and fails closed on a
    missing/nonunique root.  No external ``q_D`` is invented.  Omitting the
    interval retains the original branch, operations, bracket, and exact
    equal-temperature shortcut; any nonzero ``q_D`` returned there is an
    error rather than an after-root correction.
    """

    if not isinstance(reference, ReducedFilmSurfaceState):
        raise TypeError("a ReducedFilmSurfaceState reference is required")
    if not callable(component_flux_evaluator):
        raise TypeError("surface component-flux evaluator must be callable")

    generalized_non_fourier_closure = admissible_surface_temperature_bounds_k is not None
    if generalized_non_fourier_closure:
        bounds = admissible_surface_temperature_bounds_k
        assert bounds is not None
        if not isinstance(bounds, tuple) or len(bounds) != 2:
            raise TypeError("surface-temperature admissible bounds must be a two-value tuple")
        admissible_lower, admissible_upper = bounds
        if (
            not all(
                math.isfinite(value) and value > 0.0
                for value in (admissible_lower, admissible_upper)
            )
            or not admissible_lower < admissible_upper
        ):
            raise ValueError(
                "surface-temperature admissible bounds must be finite, positive, "
                "and strictly ordered"
            )
        if not (
            admissible_lower < reference.cell_temperature_k < admissible_upper
            and admissible_lower < reference.boundary.temperature_k < admissible_upper
        ):
            raise efr.ReducedFilmValidityError(
                "generalized surface closure bounds must strictly contain both "
                "cell and bulk temperatures"
            )

    boundary = reference.boundary
    cell_temperature = reference.cell_temperature_k
    bulk_temperature = boundary.temperature_k
    solid_resistance = reference.outer_half_cell_heat_resistance_m2_k_w
    evaluations: dict[
        float,
        tuple[
            float,
            efr.ReducedFilmCoefficients,
            efr.PartialEnthalpySecantHeatCapacities,
            efr.AckermannHeatTransfer,
            ReducedFilmSurfaceFluxEvaluation,
            float,
            float,
            float,
            float,
        ],
    ] = {}

    def evaluate(surface_temperature: float):
        cached = evaluations.get(surface_temperature)
        if cached is not None:
            return cached
        if not math.isfinite(surface_temperature):
            raise efr.ReducedFilmValidityError(
                "surface scalar closure produced a non-finite temperature"
            )
        conductive_outward = (cell_temperature - surface_temperature) / solid_resistance
        flux_evaluation = component_flux_evaluator(
            surface_temperature,
            conductive_outward,
        )
        if not isinstance(flux_evaluation, ReducedFilmSurfaceFluxEvaluation):
            raise TypeError("surface component-flux evaluator returned the wrong type")
        coefficients = _reduced_film_coefficients_at_surface_temperature(
            reference,
            surface_temperature,
        )
        exact_capacities = efr.phy053_partial_enthalpy_secant_heat_capacities(
            surface_temperature_k=surface_temperature,
            bulk_temperature_k=bulk_temperature,
            pressure_pa=coefficients.operating_point.pressure_pa,
            y_hexane=reference.surface_y_hexane,
            k_wh=boundary.binary_gas_interaction_k_wh,
        )
        ackermann = efr.ackermann_heat_transfer(
            coefficients,
            fluxes=flux_evaluation.component_mass_fluxes,
            water_film_heat_capacity_j_kg_k=exact_capacities.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=exact_capacities.hexane_j_kg_k,
        )
        film_outward = (
            coefficients.heat_transfer_w_m2_k
            * ackermann.factor
            * (surface_temperature - bulk_temperature)
        )
        reciprocal_heat = flux_evaluation.non_fourier_reduced_heat_flux_w_m2
        if reciprocal_heat == 0.0:
            reduced_heat = conductive_outward
            # Keep the established subtraction exactly unchanged at q_D=0.
            residual = conductive_outward - film_outward
        else:
            if not generalized_non_fourier_closure:
                raise efr.ReducedFilmValidityError(
                    "a nonzero surface reciprocal heat requires explicit "
                    "admissible-temperature bounds and coupled root closure"
                )
            reduced_heat = math.fsum((conductive_outward, reciprocal_heat))
            residual = math.fsum((reduced_heat, -film_outward))
        if not all(
            math.isfinite(value)
            for value in (
                conductive_outward,
                reciprocal_heat,
                reduced_heat,
                film_outward,
                residual,
            )
        ):
            raise efr.ReducedFilmValidityError(
                "surface heat closure produced a non-finite residual"
            )
        result = (
            residual,
            coefficients,
            exact_capacities,
            ackermann,
            flux_evaluation,
            conductive_outward,
            reciprocal_heat,
            reduced_heat,
            film_outward,
        )
        evaluations[surface_temperature] = result
        return result

    if not generalized_non_fourier_closure:
        lower = min(cell_temperature, bulk_temperature)
        upper = max(cell_temperature, bulk_temperature)
        if lower == upper:
            root_temperature = lower
            exact_equal_temperature_branch = True
            monotone_unique = True
        else:
            sample_temperatures = tuple(lower + (upper - lower) * index / 4.0 for index in range(5))
            sample_residuals = tuple(
                evaluate(temperature)[0] for temperature in sample_temperatures
            )
            if sample_residuals[0] < 0.0 or sample_residuals[-1] > 0.0:
                raise efr.ReducedFilmValidityError(
                    "surface heat closure lost its positive-resistance bracket"
                )
            monotone_unique = all(
                left > right for left, right in zip(sample_residuals, sample_residuals[1:])
            )
            if not monotone_unique:
                raise efr.ReducedFilmValidityError(
                    "surface heat closure failed its monotone unique-root diagnostic"
                )
            exact_nodes = tuple(
                temperature
                for temperature, residual in zip(
                    sample_temperatures,
                    sample_residuals,
                )
                if residual == 0.0
            )
            if exact_nodes:
                root_temperature = exact_nodes[0]
            else:
                root_temperature = float(
                    optimize.brentq(
                        lambda temperature: evaluate(float(temperature))[0],
                        lower,
                        upper,
                        xtol=np.nextafter(0.0, 1.0),
                        rtol=4.0 * np.finfo(float).eps,
                        maxiter=64,
                        disp=True,
                    )
                )
            exact_equal_temperature_branch = False
    else:
        assert admissible_surface_temperature_bounds_k is not None
        lower, upper = admissible_surface_temperature_bounds_k
        # The reciprocal term can displace the root outside the old
        # min/max(cell, bulk) interval and can produce a nonzero residual when
        # cell and bulk temperatures are equal.  Scan the complete caller-
        # authorized open property interval at a frozen deterministic rank.
        sample_temperatures = tuple(lower + (upper - lower) * index / 16.0 for index in range(17))
        sample_residuals = tuple(evaluate(temperature)[0] for temperature in sample_temperatures)
        if sample_residuals[0] < 0.0 or sample_residuals[-1] > 0.0:
            raise efr.ReducedFilmValidityError(
                "generalized reduced-heat surface closure has no root in the "
                "caller-supplied admissible temperature interval"
            )
        monotone_unique = all(
            left > right for left, right in zip(sample_residuals, sample_residuals[1:])
        )
        if not monotone_unique:
            raise efr.ReducedFilmValidityError(
                "generalized reduced-heat surface closure failed its monotone "
                "unique-root diagnostic"
            )
        exact_nodes = tuple(
            temperature
            for temperature, residual in zip(
                sample_temperatures,
                sample_residuals,
            )
            if residual == 0.0
        )
        if exact_nodes:
            root_temperature = exact_nodes[0]
        else:
            root_temperature = float(
                optimize.brentq(
                    lambda temperature: evaluate(float(temperature))[0],
                    lower,
                    upper,
                    xtol=np.nextafter(0.0, 1.0),
                    rtol=4.0 * np.finfo(float).eps,
                    maxiter=64,
                    disp=True,
                )
            )
        exact_equal_temperature_branch = False

    (
        closure_residual,
        coefficients,
        exact_capacities,
        ackermann,
        final_flux_evaluation,
        conductive_outward,
        reciprocal_heat,
        reduced_heat,
        film_outward,
    ) = evaluate(root_temperature)
    film_scale = (
        coefficients.heat_transfer_w_m2_k
        * ackermann.factor
        * abs(cell_temperature - bulk_temperature)
    )
    physical_heat_scale = max(
        abs(conductive_outward),
        abs(reciprocal_heat),
        abs(reduced_heat),
        abs(film_outward),
        film_scale,
        np.finfo(float).tiny,
    )
    normalized_residual = abs(closure_residual) / physical_heat_scale
    temperature_roundoff = max(
        math.ulp(cell_temperature),
        math.ulp(bulk_temperature),
        math.ulp(root_temperature),
    )
    heat_roundoff_bound = (
        8.0
        * (1.0 / solid_resistance + coefficients.heat_transfer_w_m2_k * ackermann.factor)
        * temperature_roundoff
    )
    if abs(closure_residual) > max(
        2.0e-12 * physical_heat_scale,
        heat_roundoff_bound,
    ):
        raise efr.ReducedFilmValidityError(
            "surface scalar heat closure did not converge to its tight residual contract"
        )
    external_resistance = 1.0 / (coefficients.heat_transfer_w_m2_k * ackermann.factor)
    series_coefficient = 1.0 / (solid_resistance + external_resistance)
    series_heat = series_coefficient * (cell_temperature - bulk_temperature)
    series_identity_residual = math.fsum((conductive_outward, -series_heat))
    series_scale = max(
        abs(conductive_outward),
        abs(series_heat),
        np.finfo(float).tiny,
    )
    if not generalized_non_fourier_closure and abs(series_identity_residual) > max(
        2.0e-12 * series_scale,
        2.0 * heat_roundoff_bound,
    ):
        raise efr.ReducedFilmValidityError(
            "surface closure failed the exact Ackermann series identity"
        )
    return ReducedFilmSurfaceClosure(
        reference=reference,
        coefficients=coefficients,
        surface_temperature_k=root_temperature,
        surface_y_hexane=reference.surface_y_hexane,
        surface_hexane_mass_fraction=reference.surface_hexane_mass_fraction,
        exact_partial_enthalpy_secant_heat_capacities=exact_capacities,
        ackermann=ackermann,
        surface_flux_evaluation=final_flux_evaluation,
        external_ackermann_heat_resistance_m2_k_w=external_resistance,
        series_heat_transfer_w_m2_k=series_coefficient,
        conductive_heat_flux_outward_w_m2=conductive_outward,
        fourier_heat_flux_outward_w_m2=conductive_outward,
        reciprocal_heat_of_transport_flux_outward_w_m2=reciprocal_heat,
        reduced_heat_flux_outward_w_m2=reduced_heat,
        external_ackermann_heat_flux_outward_w_m2=film_outward,
        local_reduced_heat_closure_residual_w_m2=closure_residual,
        solid_temperature_gradient_k_m=(
            (root_temperature - cell_temperature) / reference.outer_half_cell_distance_m
        ),
        heat_closure_residual_w_m2=closure_residual,
        normalized_heat_closure_residual=normalized_residual,
        heat_closure_roundoff_bound_w_m2=heat_roundoff_bound,
        temperature_bracket_k=(lower, upper),
        scalar_function_evaluations=len(evaluations),
        monotone_unique_root_diagnostic_passed=monotone_unique,
        exact_equal_temperature_branch=exact_equal_temperature_branch,
        exact_fourier_series_identity_applicable=(not generalized_non_fourier_closure),
        used_generalized_non_fourier_surface_closure=(generalized_non_fourier_closure),
    )


def audit_reduced_film_surface_flux(
    surface: ReducedFilmSurfaceClosure,
    *,
    conserved_water_flux_mol_m2_s: float,
    conserved_hexane_flux_mol_m2_s: float,
    water_gas_partial_enthalpy_j_mol: float,
    hexane_gas_partial_enthalpy_j_mol: float,
    enforce_reduction_thresholds: bool = True,
) -> SurfaceFilmAudit:
    """Audit exact-A heat transfer and the remaining fast-mass reduction.

    Algebraic, profile, and operating-band failures always fail closed.  A
    nonlinear solver may defer the caller-owned fast-mass threshold while
    evaluating uncommitted iterates; the converged candidate must enforce it
    before commitment.  Production always uses exact PHY-053 secant
    capacities.  The fixed 0.005 legacy-``cp`` series-conductance comparison
    and former unity-A shortcut are reported counterfactuals only and cannot
    reject an exact-production root.
    """

    if not isinstance(surface, ReducedFilmSurfaceClosure):
        raise TypeError("a resolved ReducedFilmSurfaceClosure is required")
    if not isinstance(enforce_reduction_thresholds, bool):
        raise TypeError("reduction-threshold enforcement flag must be boolean")
    if not all(
        math.isfinite(value)
        for value in (
            conserved_water_flux_mol_m2_s,
            conserved_hexane_flux_mol_m2_s,
            water_gas_partial_enthalpy_j_mol,
            hexane_gas_partial_enthalpy_j_mol,
        )
    ):
        raise ValueError("surface component fluxes must be finite")
    component = surface.surface_flux_evaluation.component
    if (
        conserved_water_flux_mol_m2_s != component.conserved_water_flux_mol_m2_s
        or conserved_hexane_flux_mol_m2_s != component.conserved_hexane_flux_mol_m2_s
    ):
        raise ValueError("surface audit fluxes must be the exact scalar-closure fluxes")
    reference = surface.reference
    boundary = reference.boundary
    water_mass_flux = conserved_water_flux_mol_m2_s * wa.M
    hexane_mass_flux = conserved_hexane_flux_mol_m2_s * hx.M
    external_mass_fluxes = efr.ComponentMassFluxes(
        water_kg_m2_s=water_mass_flux,
        hexane_kg_m2_s=hexane_mass_flux,
    )
    departure = efr.supplied_flux_composition_departure(
        surface.coefficients,
        bulk_hexane_mass_fraction=surface.surface_hexane_mass_fraction,
        supplied_fluxes=external_mass_fluxes,
        acceptance_threshold_mass_fraction=(boundary.absolute_mass_fraction_departure_threshold),
    )
    ackermann = surface.ackermann
    if external_mass_fluxes != surface.surface_flux_evaluation.component_mass_fluxes:
        raise RuntimeError("surface closure changed component mass-flux identity")
    exact_capacities = surface.exact_partial_enthalpy_secant_heat_capacities
    legacy_ackermann = efr.ackermann_heat_transfer(
        surface.coefficients,
        fluxes=external_mass_fluxes,
        water_film_heat_capacity_j_kg_k=(
            boundary.legacy_engineering_regression_water_heat_capacity_j_kg_k
        ),
        hexane_film_heat_capacity_j_kg_k=(
            boundary.legacy_engineering_regression_hexane_heat_capacity_j_kg_k
        ),
        enforce_production_domain=False,
    )
    legacy_external_heat_resistance = 1.0 / (
        surface.coefficients.heat_transfer_w_m2_k * legacy_ackermann.factor
    )
    legacy_series_conductance = 1.0 / (
        reference.outer_half_cell_heat_resistance_m2_k_w + legacy_external_heat_resistance
    )
    ackermann_relative_error = abs(legacy_ackermann.factor - ackermann.factor) / abs(
        ackermann.factor
    )
    series_conductance_relative_error = abs(
        legacy_series_conductance - surface.series_heat_transfer_w_m2_k
    ) / abs(surface.series_heat_transfer_w_m2_k)
    legacy_cp_series_guard_passed = (
        series_conductance_relative_error <= efr.LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT
    )
    ackermann_departure = abs(ackermann.factor - 1.0)
    fast_mass_guard_passed = departure.within_caller_threshold
    linear_heat_guard_passed = (
        ackermann_departure <= boundary.absolute_ackermann_factor_departure_threshold
    )
    # The boundary energy uses the one locally closed reduced-heat flux.  At
    # q_D=0 this is bit-identical to the historical conductive value.  For a
    # reciprocal counterfactual, q_D has already participated in the scalar
    # root and must not be appended again downstream.
    production_heat = surface.reduced_heat_flux_outward_w_m2
    full_heat = production_heat
    full_surface_temperature = surface.surface_temperature_k
    unity_heat = reference.production_conductive_heat_flux_w_m2
    water_enthalpy_outward = conserved_water_flux_mol_m2_s * water_gas_partial_enthalpy_j_mol
    hexane_enthalpy_outward = conserved_hexane_flux_mol_m2_s * hexane_gas_partial_enthalpy_j_mol
    production_energy = math.fsum(
        (
            production_heat,
            water_enthalpy_outward,
            hexane_enthalpy_outward,
        )
    )
    particle_energy_source = -production_energy
    bed_energy_source = production_energy
    audit = SurfaceFilmAudit(
        common_provenance=surface.coefficients.common_provenance,
        stress_multiplier=surface.coefficients.stress_multiplier,
        heat_transfer_w_m2_k=surface.coefficients.heat_transfer_w_m2_k,
        binary_mass_transfer_m_s=surface.coefficients.binary_mass_transfer_m_s,
        coefficient_covariance_preserved=(surface.coefficients.coefficient_covariance_preserved),
        bulk_temperature_k=boundary.temperature_k,
        surface_temperature_k=surface.surface_temperature_k,
        pressure_pa=surface.coefficients.operating_point.pressure_pa,
        binary_gas_interaction_k_wh=boundary.binary_gas_interaction_k_wh,
        bulk_y_hexane=boundary.y_hexane,
        assumed_surface_y_hexane=surface.surface_y_hexane,
        assumed_surface_hexane_mass_fraction=(surface.surface_hexane_mass_fraction),
        outer_half_cell_heat_resistance_m2_k_w=(reference.outer_half_cell_heat_resistance_m2_k_w),
        external_heat_resistance_m2_k_w=(surface.external_ackermann_heat_resistance_m2_k_w),
        series_heat_transfer_w_m2_k=surface.series_heat_transfer_w_m2_k,
        production_conductive_heat_flux_w_m2=production_heat,
        particle_fourier_heat_flux_w_m2=surface.fourier_heat_flux_outward_w_m2,
        particle_reciprocal_heat_of_transport_flux_w_m2=(
            surface.reciprocal_heat_of_transport_flux_outward_w_m2
        ),
        particle_reduced_heat_flux_w_m2=surface.reduced_heat_flux_outward_w_m2,
        external_ackermann_heat_flux_w_m2=(surface.external_ackermann_heat_flux_outward_w_m2),
        local_reduced_heat_closure_residual_w_m2=(surface.local_reduced_heat_closure_residual_w_m2),
        full_ackermann_conductive_heat_flux_w_m2=full_heat,
        full_ackermann_surface_temperature_k=full_surface_temperature,
        linear_heat_absolute_error_w_m2=abs(full_heat - unity_heat),
        unity_ackermann_conductive_heat_flux_w_m2=unity_heat,
        unity_ackermann_surface_temperature_k=(reference.surface_temperature_k),
        unity_ackermann_series_heat_transfer_w_m2_k=(reference.series_heat_transfer_w_m2_k),
        mass_departure=departure,
        ackermann=ackermann,
        exact_water_partial_enthalpy_secant_heat_capacity_j_kg_k=(exact_capacities.water_j_kg_k),
        exact_hexane_partial_enthalpy_secant_heat_capacity_j_kg_k=(exact_capacities.hexane_j_kg_k),
        exact_secant_equal_temperature_derivative_limit_used=(
            exact_capacities.equal_temperature_derivative_limit_used
        ),
        exact_secant_derivative_stencil_half_width_k=(
            exact_capacities.derivative_stencil_half_width_k
        ),
        legacy_engineering_regression_water_heat_capacity_j_kg_k=(
            boundary.legacy_engineering_regression_water_heat_capacity_j_kg_k
        ),
        legacy_engineering_regression_hexane_heat_capacity_j_kg_k=(
            boundary.legacy_engineering_regression_hexane_heat_capacity_j_kg_k
        ),
        legacy_engineering_regression_ackermann_factor=legacy_ackermann.factor,
        legacy_engineering_regression_ackermann_relative_error=(ackermann_relative_error),
        legacy_engineering_regression_series_heat_transfer_w_m2_k=(legacy_series_conductance),
        legacy_engineering_regression_series_conductance_relative_error=(
            series_conductance_relative_error
        ),
        legacy_engineering_regression_series_conductance_relative_error_threshold=(
            efr.LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT
        ),
        legacy_engineering_regression_series_conductance_guard_passed=(
            legacy_cp_series_guard_passed
        ),
        ackermann_absolute_factor_departure=ackermann_departure,
        caller_ackermann_factor_departure_threshold=(
            boundary.absolute_ackermann_factor_departure_threshold
        ),
        water_gas_partial_enthalpy_j_mol=water_gas_partial_enthalpy_j_mol,
        hexane_gas_partial_enthalpy_j_mol=hexane_gas_partial_enthalpy_j_mol,
        production_common_datum_energy_outward_w_m2=production_energy,
        particle_inventory_energy_source_w_m2=particle_energy_source,
        bed_inventory_energy_source_w_m2=bed_energy_source,
        closed_pair_energy_residual_w_m2=math.fsum((particle_energy_source, bed_energy_source)),
        fast_mass_guard_passed=fast_mass_guard_passed,
        linear_heat_guard_passed=linear_heat_guard_passed,
        validity_band_checked=surface.coefficients.validity_band_checked,
        reduction_thresholds_enforced=False,
        heat_closure_residual_w_m2=surface.heat_closure_residual_w_m2,
        normalized_heat_closure_residual=(surface.normalized_heat_closure_residual),
        heat_closure_roundoff_bound_w_m2=(surface.heat_closure_roundoff_bound_w_m2),
        surface_temperature_bracket_k=surface.temperature_bracket_k,
        scalar_closure_function_evaluations=(surface.scalar_function_evaluations),
        exact_equal_temperature_branch=surface.exact_equal_temperature_branch,
        monotone_unique_root_diagnostic_passed=(surface.monotone_unique_root_diagnostic_passed),
        exact_fourier_series_identity_applicable=(surface.exact_fourier_series_identity_applicable),
        used_generalized_non_fourier_surface_closure=(
            surface.used_generalized_non_fourier_surface_closure
        ),
    )
    return enforce_reduced_film_surface_audit(audit) if enforce_reduction_thresholds else audit


def enforce_reduced_film_surface_audit(
    audit: SurfaceFilmAudit,
) -> SurfaceFilmAudit:
    """Certify the production fast-mass reduction before state commitment.

    The reported 0.005 legacy-``cp`` comparison is diagnostic only: exact
    PHY-053 secant capacities are already production, so failure of that
    counterfactual cannot reject an otherwise valid exact-production root.
    """

    if not isinstance(audit, SurfaceFilmAudit):
        raise TypeError("a SurfaceFilmAudit is required")
    if not audit.fast_mass_guard_passed:
        raise efr.ReducedFilmValidityError(
            "fast-mass surface-equals-bulk departure exceeded the caller threshold"
        )
    if audit.reduction_thresholds_enforced:
        return audit
    return dataclass_replace(audit, reduction_thresholds_enforced=True)


def reduced_film_surface_energy_flux(
    surface: ReducedFilmSurfaceClosure,
    component_fluxes: cp.ComponentMolarFluxes,
    *,
    water_gas_partial_enthalpy_j_mol: float,
    hexane_gas_partial_enthalpy_j_mol: float,
) -> cp.ComponentEnergyFlux:
    """Collapse the zero-storage surface into one common-datum gas flux.

    Retained water reaching the external boundary crosses the gas film as
    water vapor.  Its conserved molar rate is therefore multiplied by the
    surface gas partial enthalpy exactly once; no retained-water enthalpy term
    is added at this boundary.  Internal pore faces retain the unchanged
    G1G-01 gas/retained split.
    """

    if not isinstance(surface, ReducedFilmSurfaceClosure):
        raise TypeError("a resolved ReducedFilmSurfaceClosure is required")
    if component_fluxes != surface.surface_flux_evaluation.component:
        raise ValueError("surface energy must reuse the exact closure component flux")
    values = (
        water_gas_partial_enthalpy_j_mol,
        hexane_gas_partial_enthalpy_j_mol,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("surface gas partial enthalpies must be finite")
    water_energy = component_fluxes.conserved_water_flux_mol_m2_s * water_gas_partial_enthalpy_j_mol
    hexane_energy = (
        component_fluxes.conserved_hexane_flux_mol_m2_s * hexane_gas_partial_enthalpy_j_mol
    )
    total = math.fsum(
        (
            surface.reduced_heat_flux_outward_w_m2,
            water_energy,
            hexane_energy,
        )
    )
    return cp.ComponentEnergyFlux(
        conductive_heat_flux_w_m2=(surface.reduced_heat_flux_outward_w_m2),
        gas_water_enthalpy_flux_w_m2=water_energy,
        gas_hexane_enthalpy_flux_w_m2=hexane_energy,
        retained_water_enthalpy_flux_w_m2=0.0,
        total_energy_flux_w_m2=total,
    )


def reduced_film_surface_energy_and_audit(
    surface: ReducedFilmSurfaceClosure,
    *,
    water_gas_partial_enthalpy_j_mol: float,
    hexane_gas_partial_enthalpy_j_mol: float,
    enforce_fast_mass_threshold: bool = True,
) -> tuple[cp.ComponentEnergyFlux, SurfaceFilmAudit]:
    """Return one exact-A surface energy flux and its matching audit."""

    component = surface.surface_flux_evaluation.component
    audit = audit_reduced_film_surface_flux(
        surface,
        conserved_water_flux_mol_m2_s=(component.conserved_water_flux_mol_m2_s),
        conserved_hexane_flux_mol_m2_s=(component.conserved_hexane_flux_mol_m2_s),
        water_gas_partial_enthalpy_j_mol=water_gas_partial_enthalpy_j_mol,
        hexane_gas_partial_enthalpy_j_mol=hexane_gas_partial_enthalpy_j_mol,
        enforce_reduction_thresholds=enforce_fast_mass_threshold,
    )
    energy = reduced_film_surface_energy_flux(
        surface,
        component,
        water_gas_partial_enthalpy_j_mol=water_gas_partial_enthalpy_j_mol,
        hexane_gas_partial_enthalpy_j_mol=(hexane_gas_partial_enthalpy_j_mol),
    )
    if energy.total_energy_flux_w_m2 != (audit.production_common_datum_energy_outward_w_m2):
        raise RuntimeError("surface energy and audit lost their common-datum identity")
    return energy, audit


@dataclass(frozen=True)
class FullyDryTransportState:
    """Immutable stationary-shell state and cumulative conservation ledger."""

    grid: SphericalGrid
    config: FullyDryTransportConfig
    time_s: float
    temperatures_k: tuple[float, ...]
    y_hexane: tuple[float, ...]
    reference_water_mol: float
    reference_hexane_mol: float
    reference_energy_j: float
    cumulative_boundary_water_out_mol: float = 0.0
    cumulative_boundary_hexane_out_mol: float = 0.0
    cumulative_boundary_energy_out_j: float = 0.0
    cumulative_water_inventory_change_mol: float = 0.0
    cumulative_hexane_inventory_change_mol: float = 0.0
    cumulative_energy_inventory_change_j: float = 0.0
    cumulative_absolute_water_transfer_mol: float = 0.0
    cumulative_absolute_hexane_transfer_mol: float = 0.0
    cumulative_absolute_energy_transfer_j: float = 0.0
    physically_qualifying: bool = False

    def __post_init__(self) -> None:
        if not math.isfinite(self.time_s) or self.time_s < 0.0:
            raise ValueError("state time must be finite and non-negative")
        if len(self.temperatures_k) != self.grid.n or len(self.y_hexane) != self.grid.n:
            raise ValueError("state primitives must align with the stationary grid")
        for temperature, y_hexane in zip(self.temperatures_k, self.y_hexane):
            _validate_primitive_inside_band(temperature, y_hexane, self.config)
        values = (
            self.reference_water_mol,
            self.reference_hexane_mol,
            self.reference_energy_j,
            self.cumulative_boundary_water_out_mol,
            self.cumulative_boundary_hexane_out_mol,
            self.cumulative_boundary_energy_out_j,
            self.cumulative_water_inventory_change_mol,
            self.cumulative_hexane_inventory_change_mol,
            self.cumulative_energy_inventory_change_j,
            self.cumulative_absolute_water_transfer_mol,
            self.cumulative_absolute_hexane_transfer_mol,
            self.cumulative_absolute_energy_transfer_j,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("state conservation ledger must remain finite")
        if self.reference_water_mol <= 0.0 or self.reference_hexane_mol <= 0.0:
            raise ValueError("reference component inventories must be positive")
        if any(
            value < 0.0
            for value in (
                self.cumulative_absolute_water_transfer_mol,
                self.cumulative_absolute_hexane_transfer_mol,
                self.cumulative_absolute_energy_transfer_j,
            )
        ):
            raise ValueError("cumulative absolute transfers cannot be negative")


@dataclass(frozen=True)
class FullyDryStepLedger:
    """Nonlinear, local, accepted-step, and cumulative diagnostics."""

    dt_s: float
    pressure_before_pa: float
    pressure_after_pa: float
    finite_pressure_transition: bool
    branch: DryTransportBranch
    nonlinear_evaluations: int
    nonlinear_message: str
    used_conduction_predictor: bool
    component_residual_scales_mol_s: tuple[float, ...]
    energy_residual_scales_w: tuple[float, ...]
    residual_scale_component_mol_s: float
    residual_scale_energy_w: float
    max_scaled_component_residual: float
    max_scaled_energy_residual: float
    component_balance: cp.ComponentBalanceResiduals
    stefan_recovery: cp.TotalStefanFluxRecovery
    energy_residuals_w: tuple[float, ...]
    water_step_residual_mol: float
    hexane_step_residual_mol: float
    energy_step_residual_j: float
    normalized_water_step_residual: float
    normalized_hexane_step_residual: float
    normalized_energy_step_residual: float
    normalized_water_cumulative_residual: float
    normalized_hexane_cumulative_residual: float
    normalized_energy_cumulative_residual: float
    boundary_water_out_mol: float
    boundary_hexane_out_mol: float
    boundary_energy_out_j: float
    surface_film_audit: SurfaceFilmAudit | None
    maximum_entropy_production_w_m3_k: float
    minimum_total_entropy_production_w_m3_k: float
    local_total_entropy_production_nonnegative: bool
    used_source_domain_extrapolation: bool
    physically_qualifying: bool = False

    @property
    def max_normalized_accepted_step_ledger(self) -> float:
        return max(
            self.normalized_water_step_residual,
            self.normalized_hexane_step_residual,
            self.normalized_energy_step_residual,
        )

    @property
    def max_normalized_cumulative_ledger(self) -> float:
        return max(
            self.normalized_water_cumulative_residual,
            self.normalized_hexane_cumulative_residual,
            self.normalized_energy_cumulative_residual,
        )

    @property
    def scaled_component_residuals(self) -> tuple[float, ...]:
        """Absolute independent-component residual in each physical row scale."""
        return tuple(
            abs(residual) / scale
            for residual, scale in zip(
                self.component_balance.hexane_residuals_mol_s,
                self.component_residual_scales_mol_s,
            )
        )

    @property
    def scaled_energy_residuals(self) -> tuple[float, ...]:
        """Absolute energy residual in each local capacity/transport scale."""
        return tuple(
            abs(residual) / scale
            for residual, scale in zip(
                self.energy_residuals_w,
                self.energy_residual_scales_w,
            )
        )


@dataclass(frozen=True)
class FullyDryTransportStep:
    """One accepted backward-Euler transition."""

    before: FullyDryTransportState
    after: FullyDryTransportState
    boundary: PoreBoundary
    component_face_fluxes: tuple[cp.ComponentMolarFluxes, ...]
    energy_face_fluxes: tuple[cp.ComponentEnergyFlux, ...]
    surface_film_audit: SurfaceFilmAudit | None
    ledger: FullyDryStepLedger
    physically_qualifying: bool = False


@dataclass(frozen=True)
class _Inventories:
    water_mol: tuple[float, ...]
    hexane_mol: tuple[float, ...]
    energy_j: tuple[float, ...]

    @property
    def total_water_mol(self) -> float:
        return math.fsum(self.water_mol)

    @property
    def total_hexane_mol(self) -> float:
        return math.fsum(self.hexane_mol)

    @property
    def total_energy_j(self) -> float:
        return math.fsum(self.energy_j)


@dataclass(frozen=True)
class _FaceConstitutiveData:
    independent_water_fluxes_mol_m2_s: tuple[float, ...]
    retained_water_fluxes_mol_m2_s: tuple[float, ...]
    conductive_heat_fluxes_w_m2: tuple[float, ...]
    left_face_states: tuple[cp.EquilibriumPoreState, ...]
    right_face_states: tuple[cp.EquilibriumPoreState, ...]
    centered_face_states: tuple[cp.EquilibriumPoreState, ...]
    minimum_total_entropy_production_w_m3_k: float
    maximum_total_entropy_production_w_m3_k: float
    reduced_surface_state: ReducedFilmSurfaceState | None


@dataclass(frozen=True)
class _ResidualEvaluation:
    states: tuple[cp.EquilibriumPoreState, ...]
    inventories: _Inventories
    faces: _FaceConstitutiveData
    stefan: cp.TotalStefanFluxRecovery
    component_fluxes: tuple[cp.ComponentMolarFluxes, ...]
    energy_fluxes: tuple[cp.ComponentEnergyFlux, ...]
    component_balance: cp.ComponentBalanceResiduals
    energy_residuals_w: tuple[float, ...]
    surface_film_audit: SurfaceFilmAudit | None


@dataclass(frozen=True)
class _ReducedFullyDrySurfacePayload:
    faces: _FaceConstitutiveData
    stefan: cp.TotalStefanFluxRecovery
    component_fluxes: tuple[cp.ComponentMolarFluxes, ...]


def initialize_fully_dry_state(
    grid: SphericalGrid,
    config: FullyDryTransportConfig,
    temperatures_k: float | Sequence[float],
    y_hexane: float | Sequence[float],
    *,
    time_s: float = 0.0,
) -> FullyDryTransportState:
    """Create a dry stationary-shell state from local-equilibrium storage."""
    temperatures = _expand_cell_values(temperatures_k, grid.n, "temperatures")
    compositions = _expand_cell_values(y_hexane, grid.n, "n-hexane mole fractions")
    states = _equilibrium_states(temperatures, compositions, config)
    inventories = _inventories(grid, states)
    return FullyDryTransportState(
        grid=grid,
        config=config,
        time_s=time_s,
        temperatures_k=temperatures,
        y_hexane=compositions,
        reference_water_mol=inventories.total_water_mol,
        reference_hexane_mol=inventories.total_hexane_mol,
        reference_energy_j=inventories.total_energy_j,
    )


def pressure_only_target_config(
    config: FullyDryTransportConfig,
    pressure_pa: float,
) -> FullyDryTransportConfig:
    """Return an explicit pressure-only target configuration.

    Construction re-runs the complete pressure-conditioned topology checks.
    It does not transform or project a state at the new pressure.
    """

    if not isinstance(config, FullyDryTransportConfig):
        raise TypeError("pressure target requires a fully dry transport config")
    return dataclass_replace(config, pressure_pa=pressure_pa)


def validate_pressure_only_target_config(
    before: FullyDryTransportConfig,
    target: FullyDryTransportConfig,
) -> None:
    """Reject every target-configuration change except prescribed pressure."""

    if not isinstance(before, FullyDryTransportConfig) or not isinstance(
        target,
        FullyDryTransportConfig,
    ):
        raise TypeError("pressure transition requires fully dry transport configs")
    expected = dataclass_replace(before, pressure_pa=target.pressure_pa)
    if target != expected:
        raise ValueError(
            "finite-pressure target may change pressure only; coefficients, "
            "topology charts, pore authorities, and solver contracts must be identical"
        )


def advance_fully_dry_backward_euler(
    before: FullyDryTransportState,
    dt_s: float,
    boundary: PoreBoundary,
    *,
    target_config: FullyDryTransportConfig | None = None,
) -> FullyDryTransportStep:
    """Advance one predictive ``2*N`` step or reject with exact rollback.

    The function never mutates ``before``.  All nonlinear failures, topology
    excursions, and ledger failures raise :class:`CoupledTransportStepError`
    whose ``rollback_state`` is that exact object.
    """
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError("time step must be positive and finite")
    config = before.config if target_config is None else target_config
    try:
        if not isinstance(boundary, (DirichletPoreBoundary, ReducedFilmPoreBoundary)):
            raise TypeError("boundary must be DirichletPoreBoundary or ReducedFilmPoreBoundary")
        validate_pressure_only_target_config(before.config, config)
        _validate_primitive_inside_band(
            boundary.temperature_k,
            boundary.y_hexane,
            config,
            what="boundary",
        )
    except (TypeError, ValueError, cp.CoupledPoreTopologyError) as exc:
        raise CoupledTransportStepError(
            f"finite-pressure target rejected with exact rollback: {exc}",
            before,
            nonlinear_evaluations=0,
        ) from exc
    old_states = _equilibrium_states(
        before.temperatures_k,
        before.y_hexane,
        before.config,
    )
    old_inventories = _inventories(before.grid, old_states)

    evaluation_count = 0

    def evaluate(temperatures: Sequence[float], compositions: Sequence[float]):
        nonlocal evaluation_count
        evaluation_count += 1
        return _evaluate_residual(
            before.grid,
            old_inventories,
            temperatures,
            compositions,
            dt_s,
            boundary,
            config,
        )

    try:
        initial_compositions = _pressure_target_composition_seed(
            before,
            before.temperatures_k,
            config,
        )
        initial = evaluate(before.temperatures_k, initial_compositions)
        component_scales = _component_residual_scales(
            before.grid,
            old_inventories,
            initial,
            dt_s,
        )
        energy_scales = _energy_residual_scales(
            before.grid,
            before.temperatures_k,
            initial_compositions,
            initial,
            dt_s,
            config,
        )
    except (ValueError, cp.CoupledPoreTopologyError) as exc:
        raise CoupledTransportStepError(
            f"fully dry residual scaling left its declared topology: {exc}",
            before,
            nonlinear_evaluations=evaluation_count,
        ) from exc

    exact_invariant = _is_exact_invariant(initial)
    used_conduction_predictor = False
    if exact_invariant:
        converged = initial
        nonlinear_message = "exact no-drive invariant"
    else:
        try:
            predicted_temperatures, used_conduction_predictor = _conduction_predictor(
                before,
                dt_s,
                boundary,
                config,
            )
        except (ValueError, cp.CoupledPoreTopologyError) as exc:
            raise CoupledTransportStepError(
                f"fully dry conduction predictor left its declared topology: {exc}",
                before,
                nonlinear_evaluations=evaluation_count,
            ) from exc
        try:
            predicted_compositions = _pressure_target_composition_seed(
                before,
                predicted_temperatures,
                config,
            )
            initial_coordinates = _encode_coordinates(
                predicted_temperatures,
                predicted_compositions,
                config,
            )
        except (ValueError, cp.CoupledPoreTopologyError) as exc:
            raise CoupledTransportStepError(
                f"fully dry pressure-target seed left its declared topology: {exc}",
                before,
                nonlinear_evaluations=evaluation_count,
            ) from exc

        def scaled_residual(coordinates: np.ndarray) -> np.ndarray:
            temperatures, compositions = _decode_coordinates(
                coordinates,
                before.grid.n,
                config,
            )
            try:
                candidate = evaluate(temperatures, compositions)
            except (ValueError, cp.CoupledPoreTopologyError) as exc:
                raise _TrialTopologyError(str(exc)) from exc
            return np.asarray(
                tuple(
                    value / scale
                    for value, scale in zip(
                        candidate.component_balance.hexane_residuals_mol_s,
                        component_scales,
                    )
                )
                + tuple(
                    value / scale
                    for value, scale in zip(
                        candidate.energy_residuals_w,
                        energy_scales,
                    )
                ),
                dtype=float,
            )

        try:
            solution = optimize.root(
                scaled_residual,
                initial_coordinates,
                method="hybr",
                options={
                    "xtol": config.nonlinear_step_tolerance,
                    "maxfev": config.maximum_function_evaluations,
                },
            )
            temperatures, compositions = _decode_coordinates(
                solution.x,
                before.grid.n,
                config,
            )
            converged = evaluate(temperatures, compositions)
        except (_TrialTopologyError, ValueError, cp.CoupledPoreTopologyError) as exc:
            raise CoupledTransportStepError(
                f"fully dry coupled trial left its declared topology: {exc}",
                before,
                nonlinear_evaluations=evaluation_count,
            ) from exc
        maximum_scaled_component = max(
            abs(value) / scale
            for value, scale in zip(
                converged.component_balance.hexane_residuals_mol_s,
                component_scales,
            )
        )
        maximum_scaled_energy = max(
            abs(value) / scale for value, scale in zip(converged.energy_residuals_w, energy_scales)
        )
        if (
            not solution.success
            or maximum_scaled_component > config.nonlinear_residual_tolerance
            or maximum_scaled_energy > config.nonlinear_residual_tolerance
        ):
            raise CoupledTransportStepError(
                "fully dry coupled solve rejected: "
                f"success={solution.success}, component={maximum_scaled_component:.3e}, "
                f"energy={maximum_scaled_energy:.3e}, message={solution.message}",
                before,
                nonlinear_evaluations=evaluation_count,
            )
        nonlinear_message = str(solution.message)
    return _accept_step(
        before,
        boundary,
        dt_s,
        old_inventories,
        converged,
        component_scales,
        energy_scales,
        evaluation_count,
        nonlinear_message,
        used_conduction_predictor,
        config,
    )


class _TrialTopologyError(RuntimeError):
    """Internal marker converting an inadmissible trial into step rollback."""


def _recover_and_compose_fully_dry_fluxes(
    grid: SphericalGrid,
    old: _Inventories,
    new: _Inventories,
    faces: _FaceConstitutiveData,
    dt_s: float,
) -> tuple[
    cp.TotalStefanFluxRecovery,
    tuple[cp.ComponentMolarFluxes, ...],
]:
    stefan = cp.recover_total_stefan_flux(
        grid,
        old.water_mol,
        new.water_mol,
        old.hexane_mol,
        new.hexane_mol,
        faces.retained_water_fluxes_mol_m2_s,
        dt_s,
    )
    component_fluxes: list[cp.ComponentMolarFluxes] = []
    for (
        left_face_state,
        right_face_state,
        centered_face_state,
        total_stefan,
        independent_water,
        retained_water,
    ) in zip(
        faces.left_face_states,
        faces.right_face_states,
        faces.centered_face_states,
        stefan.total_stefan_fluxes_mol_m2_s,
        faces.independent_water_fluxes_mol_m2_s,
        faces.retained_water_fluxes_mol_m2_s,
    ):
        advective_state = (
            left_face_state
            if total_stefan > 0.0
            else right_face_state
            if total_stefan < 0.0
            else centered_face_state
        )
        component_fluxes.append(
            cp.compose_component_fluxes(
                advective_state.y_water,
                advective_state.y_hexane,
                total_stefan,
                independent_water,
                retained_water,
            )
        )
    return stefan, tuple(component_fluxes)


def _evaluate_residual(
    grid: SphericalGrid,
    old: _Inventories,
    temperatures_k: Sequence[float],
    y_hexane: Sequence[float],
    dt_s: float,
    boundary: PoreBoundary,
    config: FullyDryTransportConfig,
) -> _ResidualEvaluation:
    states = _equilibrium_states(temperatures_k, y_hexane, config)
    new = _inventories(grid, states)
    reduced_surface_closure: ReducedFilmSurfaceClosure | None = None
    reduced_surface_energy: cp.ComponentEnergyFlux | None = None
    surface_film_audit: SurfaceFilmAudit | None = None
    if isinstance(boundary, DirichletPoreBoundary):
        faces = _face_constitutive_data(grid, states, boundary, config)
        stefan, component_fluxes = _recover_and_compose_fully_dry_fluxes(
            grid,
            old,
            new,
            faces,
            dt_s,
        )
    else:
        reference = prepare_reduced_film_surface_state(
            boundary,
            cell_temperature_k=states[-1].temperature_k,
            pressure_pa=config.pressure_pa,
            particle_radius_m=grid.R,
            outer_half_cell_distance_m=grid.R - grid.centers[-1],
            particle_thermal_conductivity_w_m_k=(config.thermal_conductivity.value_w_m_k),
            binary_gas_interaction_k_wh=config.pore.k_wh,
        )

        def evaluate_surface_flux(
            surface_temperature_k: float,
            conductive_outward_w_m2: float,
        ) -> ReducedFilmSurfaceFluxEvaluation:
            trial_faces = _face_constitutive_data(
                grid,
                states,
                boundary,
                config,
                reduced_surface_reference=reference,
                reduced_surface_temperature_k=surface_temperature_k,
                reduced_conductive_heat_flux_w_m2=conductive_outward_w_m2,
            )
            trial_stefan, trial_components = _recover_and_compose_fully_dry_fluxes(
                grid,
                old,
                new,
                trial_faces,
                dt_s,
            )
            return ReducedFilmSurfaceFluxEvaluation(
                component=trial_components[-1],
                payload=_ReducedFullyDrySurfacePayload(
                    faces=trial_faces,
                    stefan=trial_stefan,
                    component_fluxes=trial_components,
                ),
            )

        reduced_surface_closure = solve_reduced_film_surface_closure(
            reference,
            evaluate_surface_flux,
        )
        payload = reduced_surface_closure.surface_flux_evaluation.payload
        if not isinstance(payload, _ReducedFullyDrySurfacePayload):
            raise RuntimeError("fully dry surface closure lost its final flux payload")
        faces = payload.faces
        stefan = payload.stefan
        component_fluxes = payload.component_fluxes
        reduced_surface_thermo = cp.evaluate_equilibrium(
            reduced_surface_closure.surface_temperature_k,
            config.pressure_pa,
            reduced_surface_closure.surface_y_hexane,
            config.pore,
        )
        reduced_surface_energy, surface_film_audit = reduced_film_surface_energy_and_audit(
            reduced_surface_closure,
            water_gas_partial_enthalpy_j_mol=(
                reduced_surface_thermo.water_gas_partial_enthalpy_j_mol
            ),
            hexane_gas_partial_enthalpy_j_mol=(
                reduced_surface_thermo.hexane_gas_partial_enthalpy_j_mol
            ),
            enforce_fast_mass_threshold=False,
        )
    component_balance = cp.component_balance_residuals(
        grid,
        old.water_mol,
        new.water_mol,
        old.hexane_mol,
        new.hexane_mol,
        component_fluxes,
        dt_s,
    )
    energy_fluxes_list: list[cp.ComponentEnergyFlux] = []
    for face_index, (
        component_flux,
        conductive_heat_flux,
        left_state,
        right_state,
        center_state,
    ) in enumerate(
        zip(
            component_fluxes,
            faces.conductive_heat_fluxes_w_m2,
            faces.left_face_states,
            faces.right_face_states,
            faces.centered_face_states,
        )
    ):
        if (
            reduced_surface_closure is not None
            and reduced_surface_energy is not None
            and face_index == len(component_fluxes) - 1
        ):
            if component_flux != reduced_surface_closure.surface_flux_evaluation.component:
                raise RuntimeError("fully dry surface closure lost component identity")
            energy_fluxes_list.append(reduced_surface_energy)
            continue
        water_enthalpy_state = _upwind_state(
            component_flux.gas_water_flux_mol_m2_s,
            left_state,
            right_state,
            center_state,
        )
        hexane_enthalpy_state = _upwind_state(
            component_flux.gas_hexane_flux_mol_m2_s,
            left_state,
            right_state,
            center_state,
        )
        retained_enthalpy_state = _upwind_state(
            component_flux.retained_water_flux_mol_m2_s,
            left_state,
            right_state,
            center_state,
        )
        energy_fluxes_list.append(
            cp.component_energy_flux(
                component_flux,
                conductive_heat_flux,
                water_enthalpy_state.water_gas_partial_enthalpy_j_mol,
                hexane_enthalpy_state.hexane_gas_partial_enthalpy_j_mol,
                retained_enthalpy_state.retained_water_enthalpy_j_mol,
            )
        )
    energy_fluxes = tuple(energy_fluxes_list)
    energy_rates = tuple(
        area * flux.total_energy_flux_w_m2 for area, flux in zip(grid.areas, energy_fluxes)
    )
    energy_residuals = tuple(
        math.fsum(
            (
                (new.energy_j[cell] - old.energy_j[cell]) / dt_s,
                energy_rates[cell + 1],
                -energy_rates[cell],
            )
        )
        for cell in range(grid.n)
    )
    return _ResidualEvaluation(
        states=states,
        inventories=new,
        faces=faces,
        stefan=stefan,
        component_fluxes=component_fluxes,
        energy_fluxes=energy_fluxes,
        component_balance=component_balance,
        energy_residuals_w=energy_residuals,
        surface_film_audit=surface_film_audit,
    )


def _face_constitutive_data(
    grid: SphericalGrid,
    states: Sequence[cp.EquilibriumPoreState],
    boundary: PoreBoundary,
    config: FullyDryTransportConfig,
    *,
    reduced_surface_reference: ReducedFilmSurfaceState | None = None,
    reduced_surface_temperature_k: float | None = None,
    reduced_conductive_heat_flux_w_m2: float | None = None,
) -> _FaceConstitutiveData:
    independent = [0.0]
    retained = [0.0]
    conductive = [0.0]
    left_face_states = [states[0]]
    right_face_states = [states[0]]
    centered_face_states = [states[0]]
    total_entropy_productions = [0.0]
    reduced_surface_state: ReducedFilmSurfaceState | None = None

    def append_face(
        left: cp.EquilibriumPoreState,
        right_temperature_k: float,
        right_y_hexane: float,
        distance_m: float,
        *,
        boundary_face: bool,
        conductive_heat_flux_override_w_m2: float | None = None,
    ) -> None:
        face_temperature = 0.5 * (left.temperature_k + right_temperature_k)
        left_force_state = cp.evaluate_equilibrium(
            face_temperature,
            config.pressure_pa,
            left.y_hexane,
            config.pore,
        )
        right_force_state = cp.evaluate_equilibrium(
            face_temperature,
            config.pressure_pa,
            right_y_hexane,
            config.pore,
        )
        binary_face = tc.symmetric_binary_face_state(
            left_force_state.binary_gas.molar_density_mol_m3,
            left_force_state.y_water,
            left_force_state.y_hexane,
            right_force_state.binary_gas.molar_density_mol_m3,
            right_force_state.y_water,
            right_force_state.y_hexane,
        )
        evaluated_mobility = tc.evaluate_binary_mobility_interval(
            config.binary_diffusivity.interval,
            binary_face,
        ).select(config.binary_diffusivity.fraction)
        gas_flux = cp.independent_maxwell_stefan_flux(
            left_force_state,
            right_force_state,
            distance_m,
            evaluated_mobility.as_coupled_pore_selection(),
        )
        water_flux = cp.retained_water_flux(
            left_force_state,
            right_force_state,
            distance_m,
            config.retained_water_mobility,
        )
        face_y_hexane = right_y_hexane if boundary_face else 0.5 * (left.y_hexane + right_y_hexane)
        face_state = cp.evaluate_equilibrium(
            face_temperature,
            config.pressure_pa,
            face_y_hexane,
            config.pore,
        )
        independent.append(gas_flux.water_diffusive_flux_mol_m2_s)
        retained.append(water_flux.flux_mol_m2_s)
        temperature_gradient = (right_temperature_k - left.temperature_k) / distance_m
        conductive_heat_flux = (
            -config.thermal_conductivity.value_w_m_k * temperature_gradient
            if conductive_heat_flux_override_w_m2 is None
            else conductive_heat_flux_override_w_m2
        )
        conductive.append(conductive_heat_flux)
        left_face_states.append(left_force_state)
        right_face_states.append(right_force_state)
        centered_face_states.append(face_state)
        conduction_entropy = (
            config.thermal_conductivity.value_w_m_k
            * temperature_gradient
            * temperature_gradient
            / (left.temperature_k * right_temperature_k)
        )
        total_entropy_productions.append(
            math.fsum(
                (
                    gas_flux.entropy_production_w_m3_k,
                    water_flux.entropy_production_w_m3_k,
                    conduction_entropy,
                )
            )
        )

    for left, right, distance in zip(
        states[:-1],
        states[1:],
        (
            right_center - left_center
            for left_center, right_center in zip(grid.centers, grid.centers[1:])
        ),
    ):
        append_face(
            left,
            right.temperature_k,
            right.y_hexane,
            distance,
            boundary_face=False,
        )
    if isinstance(boundary, DirichletPoreBoundary):
        # Preserve the original numerical-oracle path exactly.
        append_face(
            states[-1],
            boundary.temperature_k,
            boundary.y_hexane,
            grid.R - grid.centers[-1],
            boundary_face=True,
        )
    else:
        if (
            reduced_surface_reference is None
            or reduced_surface_temperature_k is None
            or reduced_conductive_heat_flux_w_m2 is None
        ):
            raise RuntimeError("reduced production faces require the exact scalar surface closure")
        reduced_surface_state = reduced_surface_reference
        surface_temperature = reduced_surface_temperature_k
        conductive_heat = reduced_conductive_heat_flux_w_m2
        _validate_primitive_inside_band(
            surface_temperature,
            reduced_surface_state.surface_y_hexane,
            config,
            what="reduced-film particle surface",
        )
        append_face(
            states[-1],
            surface_temperature,
            reduced_surface_state.surface_y_hexane,
            grid.R - grid.centers[-1],
            boundary_face=True,
            conductive_heat_flux_override_w_m2=conductive_heat,
        )
    return _FaceConstitutiveData(
        independent_water_fluxes_mol_m2_s=tuple(independent),
        retained_water_fluxes_mol_m2_s=tuple(retained),
        conductive_heat_fluxes_w_m2=tuple(conductive),
        left_face_states=tuple(left_face_states),
        right_face_states=tuple(right_face_states),
        centered_face_states=tuple(centered_face_states),
        minimum_total_entropy_production_w_m3_k=min(total_entropy_productions),
        maximum_total_entropy_production_w_m3_k=max(total_entropy_productions),
        reduced_surface_state=reduced_surface_state,
    )


def _accept_step(
    before: FullyDryTransportState,
    boundary: PoreBoundary,
    dt_s: float,
    old: _Inventories,
    converged: _ResidualEvaluation,
    component_scales: Sequence[float],
    energy_scales: Sequence[float],
    nonlinear_evaluations: int,
    nonlinear_message: str,
    used_conduction_predictor: bool,
    target_config: FullyDryTransportConfig,
) -> FullyDryTransportStep:
    if converged.surface_film_audit is not None:
        try:
            strict_surface_audit = enforce_reduced_film_surface_audit(converged.surface_film_audit)
        except efr.ReducedFilmValidityError as exc:
            raise CoupledTransportStepError(
                f"fully dry converged root failed its reduced-film acceptance guard: {exc}",
                before,
                nonlinear_evaluations=nonlinear_evaluations,
            ) from exc
        converged = dataclass_replace(
            converged,
            surface_film_audit=strict_surface_audit,
        )
    area = before.grid.areas[-1]
    boundary_components = converged.component_fluxes[-1]
    boundary_energy_flux = converged.energy_fluxes[-1].total_energy_flux_w_m2
    water_out = dt_s * area * boundary_components.conserved_water_flux_mol_m2_s
    hexane_out = dt_s * area * boundary_components.conserved_hexane_flux_mol_m2_s
    energy_out = dt_s * area * boundary_energy_flux
    # Difference shell by shell before summation.  Subtracting two global
    # inventories would discard the tiny physical transfer in this stiff
    # sorption problem before the conservative telescope is even audited.
    water_change = math.fsum(
        new - previous for new, previous in zip(converged.inventories.water_mol, old.water_mol)
    )
    hexane_change = math.fsum(
        new - previous for new, previous in zip(converged.inventories.hexane_mol, old.hexane_mol)
    )
    energy_change = math.fsum(
        new - previous for new, previous in zip(converged.inventories.energy_j, old.energy_j)
    )
    water_step_residual = math.fsum((water_change, water_out))
    hexane_step_residual = math.fsum((hexane_change, hexane_out))
    energy_step_residual = math.fsum((energy_change, energy_out))
    normalized_water_step = _normalized_residual(
        water_step_residual,
        old.total_water_mol,
        converged.inventories.total_water_mol,
        water_out,
    )
    normalized_hexane_step = _normalized_residual(
        hexane_step_residual,
        old.total_hexane_mol,
        converged.inventories.total_hexane_mol,
        hexane_out,
    )
    normalized_energy_step = _normalized_residual(
        energy_step_residual,
        old.total_energy_j,
        converged.inventories.total_energy_j,
        energy_out,
    )

    cumulative_water_out = before.cumulative_boundary_water_out_mol + water_out
    cumulative_hexane_out = before.cumulative_boundary_hexane_out_mol + hexane_out
    cumulative_energy_out = before.cumulative_boundary_energy_out_j + energy_out
    cumulative_water_change = math.fsum(
        (before.cumulative_water_inventory_change_mol, water_change)
    )
    cumulative_hexane_change = math.fsum(
        (before.cumulative_hexane_inventory_change_mol, hexane_change)
    )
    cumulative_energy_change = math.fsum(
        (before.cumulative_energy_inventory_change_j, energy_change)
    )
    absolute_water = before.cumulative_absolute_water_transfer_mol + abs(water_out)
    absolute_hexane = before.cumulative_absolute_hexane_transfer_mol + abs(hexane_out)
    absolute_energy = before.cumulative_absolute_energy_transfer_j + abs(energy_out)
    water_cumulative_residual = math.fsum((cumulative_water_change, cumulative_water_out))
    hexane_cumulative_residual = math.fsum((cumulative_hexane_change, cumulative_hexane_out))
    energy_cumulative_residual = math.fsum((cumulative_energy_change, cumulative_energy_out))
    normalized_water_cumulative = _normalized_residual(
        water_cumulative_residual,
        before.reference_water_mol,
        converged.inventories.total_water_mol,
        cumulative_water_change,
        cumulative_water_out,
        absolute_water,
    )
    normalized_hexane_cumulative = _normalized_residual(
        hexane_cumulative_residual,
        before.reference_hexane_mol,
        converged.inventories.total_hexane_mol,
        cumulative_hexane_change,
        cumulative_hexane_out,
        absolute_hexane,
    )
    normalized_energy_cumulative = _normalized_residual(
        energy_cumulative_residual,
        before.reference_energy_j,
        converged.inventories.total_energy_j,
        cumulative_energy_change,
        cumulative_energy_out,
        absolute_energy,
    )
    max_scaled_component = max(
        abs(value) / scale
        for value, scale in zip(
            converged.component_balance.hexane_residuals_mol_s,
            component_scales,
        )
    )
    max_scaled_energy = max(
        abs(value) / scale for value, scale in zip(converged.energy_residuals_w, energy_scales)
    )
    after = FullyDryTransportState(
        grid=before.grid,
        config=target_config,
        time_s=before.time_s + dt_s,
        temperatures_k=tuple(state.temperature_k for state in converged.states),
        y_hexane=tuple(state.y_hexane for state in converged.states),
        reference_water_mol=before.reference_water_mol,
        reference_hexane_mol=before.reference_hexane_mol,
        reference_energy_j=before.reference_energy_j,
        cumulative_boundary_water_out_mol=cumulative_water_out,
        cumulative_boundary_hexane_out_mol=cumulative_hexane_out,
        cumulative_boundary_energy_out_j=cumulative_energy_out,
        cumulative_water_inventory_change_mol=cumulative_water_change,
        cumulative_hexane_inventory_change_mol=cumulative_hexane_change,
        cumulative_energy_inventory_change_j=cumulative_energy_change,
        cumulative_absolute_water_transfer_mol=absolute_water,
        cumulative_absolute_hexane_transfer_mol=absolute_hexane,
        cumulative_absolute_energy_transfer_j=absolute_energy,
    )
    ledger = FullyDryStepLedger(
        dt_s=dt_s,
        pressure_before_pa=before.config.pressure_pa,
        pressure_after_pa=target_config.pressure_pa,
        finite_pressure_transition=(target_config.pressure_pa != before.config.pressure_pa),
        branch=target_config.branch,
        nonlinear_evaluations=nonlinear_evaluations,
        nonlinear_message=nonlinear_message,
        used_conduction_predictor=used_conduction_predictor,
        component_residual_scales_mol_s=tuple(component_scales),
        energy_residual_scales_w=tuple(energy_scales),
        residual_scale_component_mol_s=max(component_scales),
        residual_scale_energy_w=max(energy_scales),
        max_scaled_component_residual=max_scaled_component,
        max_scaled_energy_residual=max_scaled_energy,
        component_balance=converged.component_balance,
        stefan_recovery=converged.stefan,
        energy_residuals_w=converged.energy_residuals_w,
        water_step_residual_mol=water_step_residual,
        hexane_step_residual_mol=hexane_step_residual,
        energy_step_residual_j=energy_step_residual,
        normalized_water_step_residual=normalized_water_step,
        normalized_hexane_step_residual=normalized_hexane_step,
        normalized_energy_step_residual=normalized_energy_step,
        normalized_water_cumulative_residual=normalized_water_cumulative,
        normalized_hexane_cumulative_residual=normalized_hexane_cumulative,
        normalized_energy_cumulative_residual=normalized_energy_cumulative,
        boundary_water_out_mol=water_out,
        boundary_hexane_out_mol=hexane_out,
        boundary_energy_out_j=energy_out,
        surface_film_audit=converged.surface_film_audit,
        maximum_entropy_production_w_m3_k=(converged.faces.maximum_total_entropy_production_w_m3_k),
        minimum_total_entropy_production_w_m3_k=(
            converged.faces.minimum_total_entropy_production_w_m3_k
        ),
        local_total_entropy_production_nonnegative=(
            converged.faces.minimum_total_entropy_production_w_m3_k >= 0.0
        ),
        used_source_domain_extrapolation=any(
            state.has_source_domain_extrapolation
            for state in (
                *converged.states,
                *converged.faces.left_face_states,
                *converged.faces.right_face_states,
                *converged.faces.centered_face_states,
            )
        ),
    )
    if (
        ledger.max_scaled_component_residual > target_config.nonlinear_residual_tolerance
        or ledger.max_scaled_energy_residual > target_config.nonlinear_residual_tolerance
        or ledger.max_normalized_accepted_step_ledger > target_config.ledger_tolerance
        or ledger.max_normalized_cumulative_ledger > target_config.ledger_tolerance
        or ledger.stefan_recovery.normalized_telescope_residual > target_config.ledger_tolerance
        or not ledger.local_total_entropy_production_nonnegative
    ):
        raise CoupledTransportStepError(
            "fully dry coupled step failed its acceptance ledger: "
            f"nonlinear={max(max_scaled_component, max_scaled_energy):.3e}, "
            f"step={ledger.max_normalized_accepted_step_ledger:.3e}, "
            f"cumulative={ledger.max_normalized_cumulative_ledger:.3e}",
            before,
            nonlinear_evaluations=nonlinear_evaluations,
        )
    return FullyDryTransportStep(
        before=before,
        after=after,
        boundary=boundary,
        component_face_fluxes=converged.component_fluxes,
        energy_face_fluxes=converged.energy_fluxes,
        surface_film_audit=converged.surface_film_audit,
        ledger=ledger,
    )


def _upwind_state(
    outward_flux: float,
    left: cp.EquilibriumPoreState,
    right: cp.EquilibriumPoreState,
    centered: cp.EquilibriumPoreState,
) -> cp.EquilibriumPoreState:
    """Return the donor state for one signed component flux."""
    if outward_flux > 0.0:
        return left
    if outward_flux < 0.0:
        return right
    return centered


def _conduction_predictor(
    before: FullyDryTransportState,
    dt_s: float,
    boundary: PoreBoundary,
    target_config: FullyDryTransportConfig,
) -> tuple[tuple[float, ...], bool]:
    """Return a conservative heat-only BE predictor for the nonlinear solve.

    This state is never installed or accepted on its own.  It supplies the
    thermal penetration scale that a finite-difference Jacobian cannot infer
    reliably from a perfectly uniform state when the component residual is
    initially exactly zero.  The final solve still closes both components and
    common-datum energy using the full nonlinear constitutive model.
    """
    grid = before.grid
    t_lo, t_hi = _solver_temperature_bounds(target_config)
    capacities: list[float] = []
    for temperature, y_hexane in zip(before.temperatures_k, before.y_hexane):
        try:
            capacity = _equilibrium_energy_capacity(
                temperature,
                y_hexane,
                before.config,
            )
        except ValueError:
            return before.temperatures_k, False
        capacities.append(capacity)

    matrix = np.zeros((grid.n, grid.n), dtype=float)
    right_hand_side = np.zeros(grid.n, dtype=float)
    conductivity = target_config.thermal_conductivity.value_w_m_k
    for cell in range(grid.n):
        storage_rate = capacities[cell] * grid.volumes[cell] / dt_s
        matrix[cell, cell] += storage_rate
        right_hand_side[cell] += storage_rate * before.temperatures_k[cell]
        if cell > 0:
            left_conductance = (
                grid.areas[cell] * conductivity / (grid.centers[cell] - grid.centers[cell - 1])
            )
            matrix[cell, cell] += left_conductance
            matrix[cell, cell - 1] -= left_conductance
        if cell < grid.n - 1:
            right_conductance = (
                grid.areas[cell + 1] * conductivity / (grid.centers[cell + 1] - grid.centers[cell])
            )
            matrix[cell, cell] += right_conductance
            matrix[cell, cell + 1] -= right_conductance
        else:
            if isinstance(boundary, DirichletPoreBoundary):
                # Preserve the original Dirichlet predictor arithmetic.
                boundary_conductance = grid.areas[-1] * conductivity / (grid.R - grid.centers[-1])
            else:
                surface = prepare_reduced_film_surface_state(
                    boundary,
                    cell_temperature_k=before.temperatures_k[-1],
                    pressure_pa=target_config.pressure_pa,
                    particle_radius_m=grid.R,
                    outer_half_cell_distance_m=grid.R - grid.centers[-1],
                    particle_thermal_conductivity_w_m_k=conductivity,
                    binary_gas_interaction_k_wh=target_config.pore.k_wh,
                )
                boundary_conductance = grid.areas[-1] * surface.series_heat_transfer_w_m2_k
            matrix[cell, cell] += boundary_conductance
            right_hand_side[cell] += boundary_conductance * boundary.temperature_k
    try:
        predicted = np.linalg.solve(matrix, right_hand_side)
    except np.linalg.LinAlgError:
        return before.temperatures_k, False
    if not np.all(np.isfinite(predicted)):
        return before.temperatures_k, False
    predicted_temperatures = tuple(float(value) for value in predicted)
    if not all(t_lo < value < t_hi for value in predicted_temperatures):
        return before.temperatures_k, False
    return predicted_temperatures, True


def _equilibrium_energy_capacity(
    temperature_k: float,
    y_hexane: float,
    config: FullyDryTransportConfig,
) -> float:
    """Return ``|d(U-H_ref*C_h)/dT|_y`` for scaling/prediction only.

    This fixed-``y`` derivative is a coordinate storage derivative, not a
    closed-cell heat capacity: local-equilibrium sorbate inventories also
    change along it, so its sign need not be positive over the complete
    gas-only chart.  Only its magnitude is used for a residual scale and an
    uncommitted heat-only predictor.  Accepted energy remains the unchanged
    simultaneous component/common-datum balance.
    """
    t_lo, t_hi = _solver_temperature_bounds(config)
    difference = min(
        1.0e-3,
        0.25 * (temperature_k - t_lo),
        0.25 * (t_hi - temperature_k),
    )
    if difference <= math.sqrt(np.finfo(float).eps) * temperature_k:
        raise ValueError("temperature is too close to the predictor stencil boundary")
    hotter = cp.evaluate_equilibrium(
        temperature_k + difference,
        config.pressure_pa,
        y_hexane,
        config.pore,
    )
    colder = cp.evaluate_equilibrium(
        temperature_k - difference,
        config.pressure_pa,
        y_hexane,
        config.pore,
    )
    # A shift of the hexane caloric datum adds ``H_REF*C_h`` to the stored
    # energy.  Remove that coordinate contribution before forming a numerical
    # capacity so the solver scale is a property of the thermal state, not of
    # an arbitrary enthalpy zero.  This is the same invertible energy-row basis
    # used by the cut solver.
    hotter_reduced = math.fsum(
        (
            hotter.energy_density_j_m3,
            -hx.H_REF * hotter.total_hexane_concentration_mol_m3,
        )
    )
    colder_reduced = math.fsum(
        (
            colder.energy_density_j_m3,
            -hx.H_REF * colder.total_hexane_concentration_mol_m3,
        )
    )
    storage_derivative = (hotter_reduced - colder_reduced) / (2.0 * difference)
    if not math.isfinite(storage_derivative) or storage_derivative == 0.0:
        raise ValueError("equilibrium energy storage derivative must be finite and nonzero")
    return abs(storage_derivative)


@lru_cache(maxsize=256, typed=True)
def _pressure_conditioned_temperature_domain(
    config: FullyDryTransportConfig,
) -> PressureConditionedTemperatureDomain:
    requested = config.primitive_band.temperature_bounds_k
    if isinstance(config.primitive_band, JointPrimitiveBand):
        return PressureConditionedTemperatureDomain(
            requested,
            requested,
            config.pressure_pa,
            False,
            (requested[0], requested[0]),
        )

    lower_requested, upper_requested = requested

    def feasible(temperature_k: float) -> bool:
        try:
            cp.gas_only_composition_interval(
                temperature_k,
                config.pressure_pa,
                config.pore,
            )
        except cp.CoupledPoreTopologyError:
            return False
        return True

    samples = tuple(
        lower_requested + index * (upper_requested - lower_requested) / 64.0 for index in range(65)
    )
    statuses = tuple(feasible(temperature) for temperature in samples)
    if not any(statuses):
        raise ValueError(
            "the requested temperature interval contains no gas-only state at "
            f"P={config.pressure_pa} Pa"
        )
    first = statuses.index(True)
    if any(not status for status in statuses[first:]):
        raise ValueError(
            "the requested gas-only temperature domain is disconnected at the configured pressure"
        )
    if first == 0:
        solver_lower = lower_requested
        bracket = (lower_requested, lower_requested)
        phase_boundary = False
    else:
        infeasible_lower = samples[first - 1]
        feasible_upper = samples[first]
        for _ in range(100):
            midpoint = infeasible_lower + 0.5 * (feasible_upper - infeasible_lower)
            if midpoint in (infeasible_lower, feasible_upper):
                break
            if feasible(midpoint):
                feasible_upper = midpoint
            else:
                infeasible_lower = midpoint
        solver_lower = feasible_upper
        bracket = (infeasible_lower, feasible_upper)
        phase_boundary = True
    if not solver_lower < upper_requested or not feasible(upper_requested):
        raise ValueError(
            "the requested upper temperature does not belong to the connected gas-only domain"
        )
    return PressureConditionedTemperatureDomain(
        requested,
        (solver_lower, upper_requested),
        config.pressure_pa,
        phase_boundary,
        bracket,
    )


def _solver_temperature_bounds(
    config: FullyDryTransportConfig,
) -> tuple[float, float]:
    return config.conditioned_temperature_domain.solver_bounds_k


def _validate_declared_band(config: FullyDryTransportConfig) -> None:
    """Fail closed when the requested solver domain crosses another topology."""

    if isinstance(config.primitive_band, GasOnlyPrimitiveDomain):
        t_lo, t_hi = _solver_temperature_bounds(config)
        for fraction in (0.01, 0.25, 0.5, 0.75, 0.99):
            temperature = t_lo + fraction * (t_hi - t_lo)
            for coordinate in (-20.0, 0.0, 20.0):
                composition = cp.decode_gas_only_y(
                    temperature,
                    config.pressure_pa,
                    coordinate,
                    config.pore,
                )
                cp.evaluate_equilibrium(
                    temperature,
                    config.pressure_pa,
                    composition,
                    config.pore,
                )
        return

    t_lo, t_hi = config.primitive_band.temperature_bounds_k
    y_lo, y_hi = config.primitive_band.y_hexane_bounds
    for t_fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        temperature = t_lo + t_fraction * (t_hi - t_lo)
        for y_fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
            y_hexane = y_lo + y_fraction * (y_hi - y_lo)
            try:
                cp.evaluate_equilibrium(
                    temperature,
                    config.pressure_pa,
                    y_hexane,
                    config.pore,
                )
            except (ValueError, cp.CoupledPoreTopologyError) as exc:
                raise ValueError(
                    "declared joint primitive band crosses an unauthorized topology "
                    f"near T={temperature}, y_h={y_hexane}"
                ) from exc


def _validate_primitive_inside_band(
    temperature_k: float,
    y_hexane: float,
    config: FullyDryTransportConfig,
    *,
    what: str = "state",
) -> None:
    if not math.isfinite(temperature_k) or not math.isfinite(y_hexane):
        raise ValueError(f"{what} primitives must be finite")
    t_lo, t_hi = _solver_temperature_bounds(config)
    if not t_lo < temperature_k < t_hi:
        raise ValueError(
            f"{what} temperature must lie strictly inside its pressure-conditioned solver domain"
        )
    if isinstance(config.primitive_band, GasOnlyPrimitiveDomain):
        try:
            cp.encode_gas_only_y(
                temperature_k,
                config.pressure_pa,
                y_hexane,
                config.pore,
            )
        except cp.CoupledPoreTopologyError as exc:
            raise ValueError(f"{what} composition is outside the exact gas-only topology") from exc
        return
    y_lo, y_hi = config.primitive_band.y_hexane_bounds
    if not y_lo < y_hexane < y_hi:
        raise ValueError(f"{what} composition must lie strictly inside the declared numerical band")


def _equilibrium_states(
    temperatures_k: Sequence[float],
    y_hexane: Sequence[float],
    config: FullyDryTransportConfig,
) -> tuple[cp.EquilibriumPoreState, ...]:
    if len(temperatures_k) != len(y_hexane):
        raise ValueError("temperature and composition arrays must align")
    states = []
    for temperature, composition in zip(temperatures_k, y_hexane):
        _validate_primitive_inside_band(temperature, composition, config)
        states.append(
            cp.evaluate_equilibrium(
                temperature,
                config.pressure_pa,
                composition,
                config.pore,
            )
        )
    return tuple(states)


def _pressure_target_composition_seed(
    before: FullyDryTransportState,
    target_temperatures_k: Sequence[float],
    target_config: FullyDryTransportConfig,
) -> tuple[float, ...]:
    """Return an uncommitted target-chart seed, never an accepted remap.

    A pressure change can move an otherwise valid old composition outside the
    target gas-only interval.  For the exact gas-only chart, preserve each old
    state's dimensionless barrier coordinate and decode that coordinate at the
    target ``(T, P)``.  This only initializes the nonlinear iteration; storage,
    fluxes, and the committed state still come from the simultaneous BE root.
    """

    temperatures = tuple(target_temperatures_k)
    if len(temperatures) != before.grid.n:
        raise ValueError("pressure-target seed temperatures must align with cells")
    if target_config.pressure_pa == before.config.pressure_pa:
        return before.y_hexane
    if not isinstance(target_config.primitive_band, GasOnlyPrimitiveDomain):
        return before.y_hexane
    if not isinstance(before.config.primitive_band, GasOnlyPrimitiveDomain):
        raise ValueError(
            "a pressure transition needing an exact gas-only target chart also "
            "requires an exact gas-only old chart; do not project a rectangle"
        )
    old_coordinates = tuple(
        cp.encode_gas_only_y(
            temperature,
            before.config.pressure_pa,
            composition,
            before.config.pore,
        )
        for temperature, composition in zip(
            before.temperatures_k,
            before.y_hexane,
        )
    )
    return tuple(
        cp.decode_gas_only_y(
            temperature,
            target_config.pressure_pa,
            coordinate,
            target_config.pore,
        )
        for temperature, coordinate in zip(temperatures, old_coordinates)
    )


def _inventories(
    grid: SphericalGrid,
    states: Sequence[cp.EquilibriumPoreState],
) -> _Inventories:
    if len(states) != grid.n:
        raise ValueError("equilibrium states must align with the grid")
    return _Inventories(
        water_mol=tuple(
            volume * state.total_water_concentration_mol_m3
            for volume, state in zip(grid.volumes, states)
        ),
        hexane_mol=tuple(
            volume * state.total_hexane_concentration_mol_m3
            for volume, state in zip(grid.volumes, states)
        ),
        energy_j=tuple(
            volume * state.energy_density_j_m3 for volume, state in zip(grid.volumes, states)
        ),
    )


def _expand_cell_values(
    value: float | Sequence[float],
    count: int,
    name: str,
) -> tuple[float, ...]:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        values = (float(value),) * count
    else:
        values = tuple(value)
    if len(values) != count:
        raise ValueError(f"{name} must align with all grid cells")
    if not all(isinstance(item, (int, float)) and math.isfinite(item) for item in values):
        raise ValueError(f"{name} must be finite real values")
    return tuple(float(item) for item in values)


def _encode_coordinates(
    temperatures_k: Sequence[float],
    y_hexane: Sequence[float],
    config: FullyDryTransportConfig,
) -> np.ndarray:
    t_lo, t_hi = _solver_temperature_bounds(config)

    def encode(value: float, lower: float, upper: float) -> float:
        fraction = (value - lower) / (upper - lower)
        if not 0.0 < fraction < 1.0:
            raise ValueError("cannot encode a primitive on or outside its open band")
        return math.log(fraction) - math.log1p(-fraction)

    temperatures = tuple(temperatures_k)
    compositions = tuple(y_hexane)
    if len(temperatures) != len(compositions):
        raise ValueError("temperature and composition coordinates must align")
    temperature_coordinates = tuple(encode(value, t_lo, t_hi) for value in temperatures)
    if isinstance(config.primitive_band, GasOnlyPrimitiveDomain):
        composition_coordinates = tuple(
            cp.encode_gas_only_y(
                temperature,
                config.pressure_pa,
                composition,
                config.pore,
            )
            for temperature, composition in zip(temperatures, compositions)
        )
    else:
        y_lo, y_hi = config.primitive_band.y_hexane_bounds
        composition_coordinates = tuple(encode(value, y_lo, y_hi) for value in compositions)
    return np.asarray(
        (*temperature_coordinates, *composition_coordinates),
        dtype=float,
    )


def _decode_coordinates(
    coordinates: Sequence[float],
    count: int,
    config: FullyDryTransportConfig,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if len(coordinates) != 2 * count:
        raise ValueError("coupled coordinate vector must have exactly 2*N entries")
    if not all(math.isfinite(float(value)) for value in coordinates):
        raise ValueError("coupled coordinates must remain finite")
    t_lo, t_hi = _solver_temperature_bounds(config)
    temperatures = tuple(
        t_lo + (t_hi - t_lo) * float(expit(value)) for value in coordinates[:count]
    )
    if isinstance(config.primitive_band, GasOnlyPrimitiveDomain):
        compositions = tuple(
            cp.decode_gas_only_y(
                temperature,
                config.pressure_pa,
                float(value),
                config.pore,
            )
            for temperature, value in zip(temperatures, coordinates[count:])
        )
    else:
        y_lo, y_hi = config.primitive_band.y_hexane_bounds
        compositions = tuple(
            y_lo + (y_hi - y_lo) * float(expit(value)) for value in coordinates[count:]
        )
    return temperatures, compositions


def _component_residual_scales(
    grid: SphericalGrid,
    old: _Inventories,
    evaluation: _ResidualEvaluation,
    dt_s: float,
) -> tuple[float, ...]:
    """Return fixed row-local hexane scales for one nonlinear solve.

    Each scale is physical: the old cell inventory rate and its two adjacent
    extensive face rates.  Using the whole-particle inventory for every row
    would increasingly underweight the smaller shells as the mesh is refined.
    The old/initial values keep the nonlinear scaling fixed and smooth.
    """
    rates = tuple(
        area * flux.conserved_hexane_flux_mol_m2_s
        for area, flux in zip(grid.areas, evaluation.component_fluxes)
    )
    return tuple(
        max(
            old.hexane_mol[cell] / dt_s,
            abs(rates[cell]),
            abs(rates[cell + 1]),
        )
        for cell in range(grid.n)
    )


def _energy_residual_scales(
    grid: SphericalGrid,
    temperatures_k: Sequence[float],
    y_hexane: Sequence[float],
    evaluation: _ResidualEvaluation,
    dt_s: float,
    config: FullyDryTransportConfig,
) -> tuple[float, ...]:
    """Return fixed row-local, capacity-based common-energy scales.

    The storage term uses ``V (dU/dT) T / dt`` rather than ``|U|/dt`` so an
    arbitrary component caloric-datum shift cannot change the thermal scale.
    Adjacent extensive energy rates retain the actual local transport scale.
    """
    rates = tuple(
        area * flux.total_energy_flux_w_m2
        for area, flux in zip(grid.areas, evaluation.energy_fluxes)
    )
    capacities = tuple(
        _equilibrium_energy_capacity(temperature, composition, config)
        for temperature, composition in zip(temperatures_k, y_hexane)
    )
    return tuple(
        max(
            grid.volumes[cell] * capacities[cell] * max(abs(temperatures_k[cell]), 1.0) / dt_s,
            abs(rates[cell]),
            abs(rates[cell + 1]),
        )
        for cell in range(grid.n)
    )


def _is_exact_invariant(evaluation: _ResidualEvaluation) -> bool:
    return all(
        value == 0.0 for value in evaluation.component_balance.hexane_residuals_mol_s
    ) and all(value == 0.0 for value in evaluation.energy_residuals_w)


def _normalized_residual(residual: float, *terms: float) -> float:
    scale = max(*(abs(value) for value in terms), np.finfo(float).tiny)
    return abs(residual) / scale
