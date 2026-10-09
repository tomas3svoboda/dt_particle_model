r"""Non-qualifying retained-water/energy predictor for the wet particle.

This module is the stationary-domain foundation needed before the Gate-1g
moving-front solve is allowed to carry wet-side water.  Its deliberately
limited topology is

* one liquid-filled wet region with no pore gas;
* immutable material labels ``X_h,w^a`` and ``w_o``;
* exactly zero total relative n-hexane flux; and
* two evolving cell primitives, temperature and retained-water loading.

Retained water follows the frozen modified-Luikov inverse on its qualified
positive-moisture branch.  The wet energy is not re-derived here: replacing
only ``WetCoreParams.X_water`` and ``WetCoreParams.w_o`` cell by cell is
algebraically identical to the G1G-01 wet-core state function,

``e = h_dm + W_w h_w,l + X_h,w^a u_h,l - B_h(a_h=1)``.

The transport law is a caller-selected positive mobility multiplying the
dimensionless retained chemical-potential gradient.  Both traces at a face
are evaluated at one symmetric face temperature.  Hence the pure-liquid
standard-state term cancels and this predictor adds no unapproved Soret law.
The energy flux is Fourier heat plus retained-water molar flux times an
upwind, sign-aware retained-water enthalpy, exactly once.

The outer Dirichlet state is a numerical stress boundary, not a physical bed
closure.  Consequently every state, flux, ledger, and result in this module
reports ``physically_qualifying=False``.  No state is clipped.  A failed step
raises :class:`WetWaterStepError` and carries the exact untouched input state.
Nonzero forcing below the double-precision resolution needed to certify the
strict ledgers is likewise rejected; it is never converted to a zero-drive
state by a numerical deadband.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass, replace
from functools import lru_cache
from typing import Sequence

import numpy as np
from scipy.optimize import least_squares

from dtdc_simulator.core2.exact_cache import float_bits_key
from dtdc_simulator.core2.particle import (
    coupled_pore,
    nonisothermal_potential,
    wet_retained_cap,
    wet_core,
)
from dtdc_simulator.core2.particle.grid import SphericalGrid
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa


REFERENCE_FUGACITY_PA = 1.0


class WetWaterStepError(RuntimeError):
    """A predictor step failed without changing its input state."""

    def __init__(
        self,
        message: str,
        rollback_state: "WetWaterState",
        nonlinear_evaluations: int = 0,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state
        self.nonlinear_evaluations = nonlinear_evaluations


@dataclass(frozen=True)
class WetWaterModel:
    """Frozen state authorities plus an explicit retained-water mobility."""

    retained_water_mobility: coupled_pore.MobilitySelection
    wet: wet_core.WetCoreParams = wet_core.WetCoreParams()
    luikov: sp.LuikovParams = sp.LuikovParams()
    ledger_tolerance: float = 1.0e-10
    nonlinear_tolerance: float = 2.0e-11
    max_nonlinear_evaluations: int = 2500

    def __post_init__(self) -> None:
        # Reuse the wet-core authority's complete validation rather than
        # maintaining a second, potentially divergent parameter checklist.
        wet_core._validate_params(self.wet)
        if not isinstance(
            self.retained_water_mobility, coupled_pore.MobilitySelection
        ):
            raise TypeError("retained-water mobility must be an explicit selection")
        if not 0.0 < self.ledger_tolerance <= 1.0e-10:
            raise ValueError("ledger tolerance must be positive and no looser than 1e-10")
        if not 0.0 < self.nonlinear_tolerance < self.ledger_tolerance:
            raise ValueError("nonlinear tolerance must be tighter than the ledger")
        try:
            evaluation_limit = operator.index(self.max_nonlinear_evaluations)
        except TypeError as exc:
            raise ValueError("max nonlinear evaluations must be an integer") from exc
        if isinstance(self.max_nonlinear_evaluations, bool) or evaluation_limit < 1:
            raise ValueError("max nonlinear evaluations must be a positive integer")
        continued = isinstance(self.luikov, sp.ContinuedPositiveLuikovParams)
        if not (
            0.0 <= self.luikov.W_ref < self.luikov.W_cap
            if continued
            else 0.0 < self.luikov.W_ref < self.luikov.W_cap
        ):
            raise ValueError("Luikov retained-water interval is empty")


@dataclass(frozen=True)
class WetWaterDirichletBoundary:
    """Numerical outer-face trace; this is not a physical bed boundary."""

    temperature_k: float
    retained_water_loading: float
    label: str = "numerical Dirichlet stress trace"
    physically_qualifying: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label.strip():
            raise ValueError("Dirichlet stress boundary needs a non-empty label")
        if self.physically_qualifying:
            raise ValueError("the numerical Dirichlet boundary cannot be qualifying")


@dataclass(frozen=True)
class WetWaterCellState:
    """Complete local wet state, with no pore-gas degree of freedom."""

    temperature_k: float
    retained_water_loading: float
    historical_hexane_loading: float
    oil_fraction_label: float
    retained_water_activity: float
    retained_water_fugacity_pa: float
    retained_water_equilibrium_potential: float
    retained_water_potential: float
    retained_water_capacity_dual_over_rt: float
    retained_water_active_set: wet_retained_cap.WetRetainedWaterBranch
    retained_water_concentration_mol_m3: float
    energy_density_j_m3: float
    retained_water_enthalpy_j_mol: float
    total_hexane_concentration_kg_m3: float
    pore_gas_volume_fraction: float = 0.0
    total_relative_hexane_flux_kg_m2_s: float = 0.0
    physically_qualifying: bool = False


@dataclass(frozen=True)
class WetWaterFaceFlux:
    """Outward-positive water and energy flux at one spherical face."""

    retained_water_flux_mol_m2_s: float
    conductive_heat_flux_w_m2: float
    retained_water_energy_flux_w_m2: float
    total_energy_flux_w_m2: float
    dimensionless_force_gradient_m_inv: float
    entropy_production_w_m3_k: float
    enthalpy_donor: str
    total_relative_hexane_flux_kg_m2_s: float = 0.0
    used_explicit_nonisothermal_potential: bool = False
    retained_water_chemical_force_gradient_m_inv: float = 0.0
    retained_water_thermal_force_gradient_m_inv: float = 0.0
    retained_water_thermal_force_factor: float = 0.0
    log_temperature_gradient_m_inv: float = 0.0
    fourier_heat_flux_w_m2: float = 0.0
    reciprocal_heat_of_transport_flux_w_m2: float = 0.0
    reduced_heat_flux_w_m2: float = 0.0
    conduction_entropy_production_w_m3_k: float = 0.0
    discrete_conjugate_temperature_k: float | None = None
    physically_qualifying: bool = False

    @property
    def total_entropy_production_w_m3_k(self) -> float:
        return math.fsum(
            (
                self.entropy_production_w_m3_k,
                self.conduction_entropy_production_w_m3_k,
            )
        )


@dataclass(frozen=True)
class WetWaterState:
    """Immutable predictor state and cumulative conservation references."""

    grid: SphericalGrid
    model: WetWaterModel
    time_s: float
    temperatures_k: tuple[float, ...]
    retained_water_loadings: tuple[float, ...]
    historical_hexane_loadings: tuple[float, ...]
    oil_fraction_labels: tuple[float, ...]
    reference_temperatures_k: tuple[float, ...]
    reference_retained_water_loadings: tuple[float, ...]
    reference_water_cell_mol: tuple[float, ...]
    reference_energy_cell_j: tuple[float, ...]
    reference_hexane_mass_kg: float
    cumulative_boundary_water_out_mol: float = 0.0
    cumulative_boundary_energy_out_j: float = 0.0
    cumulative_absolute_water_exchange_mol: float = 0.0
    cumulative_absolute_energy_exchange_j: float = 0.0
    physically_qualifying: bool = False


@dataclass(frozen=True)
class WetWaterLedger:
    """Per-step and cumulative extensive balances."""

    step_water_change_mol: float
    step_water_boundary_out_mol: float
    step_water_residual_mol: float
    step_water_relative_residual: float
    step_energy_change_j: float
    step_energy_boundary_out_j: float
    step_energy_residual_j: float
    step_energy_relative_residual: float
    cumulative_water_change_mol: float
    cumulative_water_boundary_out_mol: float
    cumulative_water_residual_mol: float
    cumulative_water_relative_residual: float
    cumulative_energy_change_j: float
    cumulative_energy_boundary_out_j: float
    cumulative_energy_residual_j: float
    cumulative_energy_relative_residual: float
    max_local_relative_residual: float
    historical_hexane_change_kg: float
    max_total_relative_hexane_flux_kg_m2_s: float
    nonlinear_evaluations: int
    accepted: bool
    physically_qualifying: bool = False

    @property
    def maximum_conservation_relative_residual(self) -> float:
        return max(
            self.step_water_relative_residual,
            self.step_energy_relative_residual,
            self.cumulative_water_relative_residual,
            self.cumulative_energy_relative_residual,
            self.max_local_relative_residual,
        )


@dataclass(frozen=True)
class WetWaterStepResult:
    state: WetWaterState
    cells: tuple[WetWaterCellState, ...]
    face_fluxes: tuple[WetWaterFaceFlux, ...]
    ledger: WetWaterLedger
    physically_qualifying: bool = False


def evaluate_cell(
    temperature_k: float,
    retained_water_loading: float,
    historical_hexane_loading: float,
    oil_fraction_label: float,
    model: WetWaterModel,
    *,
    retained_water_capacity_dual_over_rt: float = 0.0,
) -> WetWaterCellState:
    """Evaluate a wet cell, memoized for bit-exact repeated inputs."""

    return _evaluate_cell_cached(
        temperature_k,
        retained_water_loading,
        historical_hexane_loading,
        oil_fraction_label,
        model,
        retained_water_capacity_dual_over_rt,
        float_bits_key(
            temperature_k,
            retained_water_loading,
            historical_hexane_loading,
            oil_fraction_label,
            retained_water_capacity_dual_over_rt,
        ),
        hx.caloric_datum_signature(),
        wa.caloric_datum_signature(),
    )


@lru_cache(maxsize=8192, typed=True)
def _evaluate_cell_cached(
    temperature_k: float,
    retained_water_loading: float,
    historical_hexane_loading: float,
    oil_fraction_label: float,
    model: WetWaterModel,
    retained_water_capacity_dual_over_rt: float,
    _scalar_key: tuple[bytes, ...],
    _hexane_datum_key: tuple[bytes, ...],
    _water_datum_key: tuple[bytes, ...],
) -> WetWaterCellState:
    """Evaluate one local state using the unmodified wet-core energy identity."""

    values = (
        temperature_k,
        retained_water_loading,
        historical_hexane_loading,
        oil_fraction_label,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("wet-water cell primitives and labels must be finite")
    if not model.wet.T_min <= temperature_k <= model.wet.T_max:
        raise ValueError("wet-water temperature is outside the caloric bracket")
    if historical_hexane_loading < 0.0:
        raise ValueError("historical wet n-hexane loading must be non-negative")

    if isinstance(model.luikov, sp.ContinuedPositiveLuikovParams) and retained_water_loading <= 0.0:
        raise ValueError("continued wet material requires positive water")
    activity = sp.water_activity(retained_water_loading, model.luikov)
    retained_graph = wet_retained_cap.classify_exact_graph_point(
        retained_water_loading,
        retained_water_capacity_dual_over_rt,
        model.luikov,
    )
    cell_params = replace(
        model.wet,
        X_water=retained_water_loading,
        w_o=oil_fraction_label,
    )
    # This call both preserves the frozen common-datum state function and
    # fails closed if the immutable historical hexane label cannot support the
    # saturated, no-pore-gas wet topology at this temperature.
    energy_density = wet_core.energy_density(
        temperature_k, historical_hexane_loading, cell_params
    )
    water_liquid = wa.state_Tp(temperature_k, model.wet.pressure_pa, "liquid")
    fugacity = activity * water_liquid.fugacity
    continued = isinstance(model.luikov, sp.ContinuedPositiveLuikovParams)
    if not math.isfinite(fugacity) or (fugacity <= 0.0 and not continued):
        raise ValueError("retained-water fugacity must be positive and finite")
    equilibrium_potential = (
        sp.water_log_activity(retained_water_loading, model.luikov)
        + math.log(water_liquid.fugacity / REFERENCE_FUGACITY_PA)
        if continued
        else math.log(fugacity / REFERENCE_FUGACITY_PA)
    )
    generalized_potential = math.fsum(
        (equilibrium_potential, retained_water_capacity_dual_over_rt)
    )
    return WetWaterCellState(
        temperature_k=temperature_k,
        retained_water_loading=retained_water_loading,
        historical_hexane_loading=historical_hexane_loading,
        oil_fraction_label=oil_fraction_label,
        retained_water_activity=activity,
        retained_water_fugacity_pa=fugacity,
        retained_water_equilibrium_potential=equilibrium_potential,
        retained_water_potential=generalized_potential,
        retained_water_capacity_dual_over_rt=(
            retained_water_capacity_dual_over_rt
        ),
        retained_water_active_set=retained_graph.branch,
        retained_water_concentration_mol_m3=(
            model.wet.rho_dm_p * retained_water_loading / wa.M
        ),
        energy_density_j_m3=energy_density,
        retained_water_enthalpy_j_mol=water_liquid.h_mass * wa.M,
        total_hexane_concentration_kg_m3=(
            model.wet.rho_dm_p * historical_hexane_loading
        ),
    )


def clear_evaluate_cell_cache() -> None:
    """Drop every cached wet-cell state and its complete model key."""

    _evaluate_cell_cached.cache_clear()


def evaluate_cell_cache_info():
    """Return bounded-cache hit/miss/size diagnostics."""

    return _evaluate_cell_cached.cache_info()


def initialize(
    grid: SphericalGrid,
    temperatures_k: Sequence[float],
    retained_water_loadings: Sequence[float],
    historical_hexane_loadings: Sequence[float],
    model: WetWaterModel,
    *,
    oil_fraction_labels: Sequence[float] | None = None,
    time_s: float = 0.0,
) -> WetWaterState:
    """Create a validated stationary wet state and exact ledger references."""

    temperatures = _tuple_of_length("temperatures", temperatures_k, grid.n)
    water = _tuple_of_length(
        "retained-water loadings", retained_water_loadings, grid.n
    )
    hexane = _tuple_of_length(
        "historical n-hexane loadings", historical_hexane_loadings, grid.n
    )
    if oil_fraction_labels is None:
        oil = (model.wet.w_o,) * grid.n
    else:
        oil = _tuple_of_length("oil labels", oil_fraction_labels, grid.n)
    if not math.isfinite(time_s) or time_s < 0.0:
        raise ValueError("wet-water time must be finite and non-negative")
    cells = tuple(
        evaluate_cell(T, W, Xh, wo, model)
        for T, W, Xh, wo in zip(temperatures, water, hexane, oil)
    )
    reference_water = tuple(
        volume * cell.retained_water_concentration_mol_m3
        for volume, cell in zip(grid.volumes, cells)
    )
    reference_energy = tuple(
        volume * cell.energy_density_j_m3
        for volume, cell in zip(grid.volumes, cells)
    )
    reference_hexane = math.fsum(
        volume * cell.total_hexane_concentration_kg_m3
        for volume, cell in zip(grid.volumes, cells)
    )
    return WetWaterState(
        grid=grid,
        model=model,
        time_s=time_s,
        temperatures_k=temperatures,
        retained_water_loadings=water,
        historical_hexane_loadings=hexane,
        oil_fraction_labels=oil,
        reference_temperatures_k=temperatures,
        reference_retained_water_loadings=water,
        reference_water_cell_mol=reference_water,
        reference_energy_cell_j=reference_energy,
        reference_hexane_mass_kg=reference_hexane,
    )


def step(
    state: WetWaterState,
    dt_s: float,
    boundary: WetWaterDirichletBoundary,
) -> WetWaterStepResult:
    """Advance the coupled ``2N`` stationary-sphere predictor by backward Euler."""

    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise WetWaterStepError("time step must be positive and finite", state)
    try:
        old_cells = _evaluate_state_cells(state)
        _validate_boundary(boundary, state.model)
        next_time_s = state.time_s + dt_s
        if not math.isfinite(next_time_s):
            raise ValueError("wet-water time overflowed during the requested step")
    except Exception as exc:
        raise WetWaterStepError(
            f"invalid wet-water step input; exact rollback: {exc}", state
        ) from exc

    if _is_uniform_invariant(state, boundary):
        fluxes = _face_fluxes(old_cells, state.grid, boundary, state.model)
        new_state = replace(state, time_s=next_time_s)
        try:
            ledger = _build_ledger(
                state,
                new_state,
                old_cells,
                old_cells,
                fluxes,
                dt_s,
                max_local_relative_residual=0.0,
                nonlinear_evaluations=0,
            )
            _require_acceptable_ledger(ledger, state.model.ledger_tolerance)
        except Exception as exc:
            raise WetWaterStepError(
                f"uniform wet-water invariant rejected with exact rollback: {exc}",
                state,
            ) from exc
        return WetWaterStepResult(new_state, old_cells, fluxes, ledger)

    lower, upper = state.model.luikov.W_ref, state.model.luikov.W_cap
    if any(not lower < value < upper for value in state.retained_water_loadings):
        raise WetWaterStepError(
            "a driven Luikov endpoint needs the named retained-water active-set "
            "solver; this smooth predictor does not move or clip it",
            state,
        )

    try:
        coordinates = _encode_primitives(
            state.temperatures_k,
            state.retained_water_loadings,
            state.model,
        )
        residual_scales = _residual_scales(
            state.grid, old_cells, boundary, state.model, dt_s
        )
        evaluations = 0

        def residual_function(candidate: np.ndarray) -> np.ndarray:
            nonlocal evaluations
            evaluations += 1
            temperatures, water = _decode_primitives(candidate, state.model)
            cells = tuple(
                evaluate_cell(T, W, Xh, wo, state.model)
                for T, W, Xh, wo in zip(
                    temperatures,
                    water,
                    state.historical_hexane_loadings,
                    state.oil_fraction_labels,
                )
            )
            fluxes = _face_fluxes(cells, state.grid, boundary, state.model)
            normalized, _, _ = _balance_residuals(
                state.grid,
                old_cells,
                cells,
                fluxes,
                dt_s,
                residual_scales,
                state.model,
            )
            return np.asarray(normalized, dtype=float)

        solution = least_squares(
            residual_function,
            coordinates,
            method="trf",
            ftol=1.0e-13,
            xtol=1.0e-13,
            gtol=1.0e-13,
            max_nfev=state.model.max_nonlinear_evaluations,
            x_scale="jac",
        )
        temperatures, water = _decode_primitives(solution.x, state.model)
        new_cells = tuple(
            evaluate_cell(T, W, Xh, wo, state.model)
            for T, W, Xh, wo in zip(
                temperatures,
                water,
                state.historical_hexane_loadings,
                state.oil_fraction_labels,
            )
        )
        fluxes = _face_fluxes(new_cells, state.grid, boundary, state.model)
        normalized, _, _ = _balance_residuals(
            state.grid,
            old_cells,
            new_cells,
            fluxes,
            dt_s,
            residual_scales,
            state.model,
        )
        local_residual = max(abs(value) for value in normalized)
        if not solution.success or local_residual > state.model.nonlinear_tolerance:
            raise RuntimeError(
                f"2N solve did not meet its residual contract: success={solution.success}, "
                f"max relative residual={local_residual:.3e}, status={solution.status}"
            )

        provisional = replace(
            state,
            time_s=next_time_s,
            temperatures_k=tuple(temperatures),
            retained_water_loadings=tuple(water),
        )
        boundary_water = (
            dt_s
            * state.grid.areas[-1]
            * fluxes[-1].retained_water_flux_mol_m2_s
        )
        boundary_energy = (
            dt_s * state.grid.areas[-1] * fluxes[-1].total_energy_flux_w_m2
        )
        new_state = replace(
            provisional,
            cumulative_boundary_water_out_mol=(
                state.cumulative_boundary_water_out_mol + boundary_water
            ),
            cumulative_boundary_energy_out_j=(
                state.cumulative_boundary_energy_out_j + boundary_energy
            ),
            cumulative_absolute_water_exchange_mol=(
                state.cumulative_absolute_water_exchange_mol + abs(boundary_water)
            ),
            cumulative_absolute_energy_exchange_j=(
                state.cumulative_absolute_energy_exchange_j + abs(boundary_energy)
            ),
        )
        ledger = _build_ledger(
            state,
            new_state,
            old_cells,
            new_cells,
            fluxes,
            dt_s,
            max_local_relative_residual=local_residual,
            nonlinear_evaluations=evaluations,
        )
        _require_acceptable_ledger(ledger, state.model.ledger_tolerance)
        return WetWaterStepResult(new_state, new_cells, fluxes, ledger)
    except WetWaterStepError:
        raise
    except Exception as exc:
        count = locals().get("evaluations", 0)
        raise WetWaterStepError(
            f"wet-water step rejected with exact rollback: {exc}",
            state,
            nonlinear_evaluations=count,
        ) from exc


def _evaluate_state_cells(state: WetWaterState) -> tuple[WetWaterCellState, ...]:
    _validate_state_structure(state)
    return tuple(
        evaluate_cell(T, W, Xh, wo, state.model)
        for T, W, Xh, wo in zip(
            state.temperatures_k,
            state.retained_water_loadings,
            state.historical_hexane_loadings,
            state.oil_fraction_labels,
        )
    )


def _one_face_flux(
    left_temperature: float,
    left_water: float,
    right_temperature: float,
    right_water: float,
    distance_m: float,
    model: WetWaterModel,
    left_name: str,
    right_name: str,
) -> WetWaterFaceFlux:
    if not math.isfinite(distance_m) or distance_m <= 0.0:
        raise ValueError("wet-water face distance must be positive and finite")
    left_activity = sp.water_activity(left_water, model.luikov)
    right_activity = sp.water_activity(right_water, model.luikov)
    # Both chemical-potential traces use this same face temperature.  Their
    # pure-liquid standard fugacity is therefore identical and cancels
    # algebraically.  Cancelling before evaluation is equation-exact and also
    # avoids contaminating a small Luikov force by subtracting two large logs.
    gradient = (
        sp.water_log_activity(right_water, model.luikov)
        - sp.water_log_activity(left_water, model.luikov)
        if isinstance(model.luikov, sp.ContinuedPositiveLuikovParams)
        else math.log(right_activity) - math.log(left_activity)
    ) / distance_m
    water_flux = -model.retained_water_mobility.value_mol_m_s * gradient
    conductive_flux = (
        model.wet.conductivity
        * (left_temperature - right_temperature)
        / distance_m
    )
    if water_flux > 0.0:
        donor_temperature, donor_name = left_temperature, left_name
    elif water_flux < 0.0:
        donor_temperature, donor_name = right_temperature, right_name
    else:
        donor_name = "none (zero retained-water flux)"
        water_energy_flux = 0.0
    if water_flux != 0.0:
        donor_enthalpy = (
            wa.state_Tp(donor_temperature, model.wet.pressure_pa, "liquid").h_mass
            * wa.M
        )
        water_energy_flux = water_flux * donor_enthalpy
    return WetWaterFaceFlux(
        retained_water_flux_mol_m2_s=water_flux,
        conductive_heat_flux_w_m2=conductive_flux,
        retained_water_energy_flux_w_m2=water_energy_flux,
        total_energy_flux_w_m2=conductive_flux + water_energy_flux,
        dimensionless_force_gradient_m_inv=gradient,
        entropy_production_w_m3_k=(
            bg.R
            * model.retained_water_mobility.value_mol_m_s
            * gradient
            * gradient
        ),
        enthalpy_donor=donor_name,
        retained_water_chemical_force_gradient_m_inv=gradient,
        log_temperature_gradient_m_inv=(
            nonisothermal_potential.log_temperature_gradient_m_inv(
                left_temperature,
                right_temperature,
                distance_m,
            )
        ),
        fourier_heat_flux_w_m2=conductive_flux,
        reduced_heat_flux_w_m2=conductive_flux,
        conduction_entropy_production_w_m3_k=(
            model.wet.conductivity
            * ((right_temperature - left_temperature) / distance_m) ** 2
            / (left_temperature * right_temperature)
        ),
        discrete_conjugate_temperature_k=(
            nonisothermal_potential.discrete_conjugate_temperature_k(
                left_temperature,
                right_temperature,
            )
        ),
    )


def _face_fluxes(
    cells: Sequence[WetWaterCellState],
    grid: SphericalGrid,
    boundary: WetWaterDirichletBoundary,
    model: WetWaterModel,
) -> tuple[WetWaterFaceFlux, ...]:
    zero = WetWaterFaceFlux(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "center symmetry")
    faces = [zero]
    for face in range(1, grid.n):
        left, right = cells[face - 1], cells[face]
        faces.append(
            _one_face_flux(
                left.temperature_k,
                left.retained_water_loading,
                right.temperature_k,
                right.retained_water_loading,
                grid.centers[face] - grid.centers[face - 1],
                model,
                f"cell {face - 1}",
                f"cell {face}",
            )
        )
    last = cells[-1]
    faces.append(
        _one_face_flux(
            last.temperature_k,
            last.retained_water_loading,
            boundary.temperature_k,
            boundary.retained_water_loading,
            grid.R - grid.centers[-1],
            model,
            f"cell {grid.n - 1}",
            "Dirichlet boundary",
        )
    )
    return tuple(faces)


def _balance_residuals(
    grid: SphericalGrid,
    old_cells: Sequence[WetWaterCellState],
    new_cells: Sequence[WetWaterCellState],
    fluxes: Sequence[WetWaterFaceFlux],
    dt_s: float,
    residual_scales: tuple[float, float],
    model: WetWaterModel,
) -> tuple[list[float], list[float], list[float]]:
    normalized: list[float] = []
    water_residuals: list[float] = []
    energy_residuals: list[float] = []
    for i, (old, new) in enumerate(zip(old_cells, new_cells)):
        water_accumulation = _water_amount_change(
            grid.volumes[i], old, new, model.wet.rho_dm_p
        )
        water_right = (
            dt_s * grid.areas[i + 1] * fluxes[i + 1].retained_water_flux_mol_m2_s
        )
        water_left = (
            dt_s * grid.areas[i] * fluxes[i].retained_water_flux_mol_m2_s
        )
        water_residual = water_accumulation + water_right - water_left

        energy_accumulation = _energy_amount_change(
            grid.volumes[i], old, new, model.wet.rho_dm_p
        )
        energy_right = (
            dt_s * grid.areas[i + 1] * fluxes[i + 1].total_energy_flux_w_m2
        )
        energy_left = dt_s * grid.areas[i] * fluxes[i].total_energy_flux_w_m2
        energy_residual = energy_accumulation + energy_right - energy_left
        water_residuals.append(water_residual)
        energy_residuals.append(energy_residual)
        normalized.extend(
            (
                water_residual / residual_scales[0],
                energy_residual / residual_scales[1],
            )
        )
    return normalized, water_residuals, energy_residuals


def _water_amount_change(
    volume_m3: float,
    old: WetWaterCellState,
    new: WetWaterCellState,
    rho_dm_p: float,
) -> float:
    """Stable molar accumulation from primitive loadings, without subtraction."""

    return (
        volume_m3
        * rho_dm_p
        / wa.M
        * (new.retained_water_loading - old.retained_water_loading)
    )


def _energy_amount_change(
    volume_m3: float,
    old: WetWaterCellState,
    new: WetWaterCellState,
    rho_dm_p: float,
) -> float:
    """Stable, common-datum-exact wet-energy accumulation.

    Wet energy is affine in retained-water loading with coefficient ``h_w,l``.
    Removing that term before differencing avoids cancellation of two large
    energy densities; adding its exact product difference back preserves the
    frozen wet-core state function and its arbitrary common datum.
    """

    if (
        old.historical_hexane_loading != new.historical_hexane_loading
        or old.oil_fraction_label != new.oil_fraction_label
    ):
        raise ValueError("wet material labels changed inside an energy increment")
    old_h_mass = old.retained_water_enthalpy_j_mol / wa.M
    new_h_mass = new.retained_water_enthalpy_j_mol / wa.M
    if new.temperature_k == old.temperature_k:
        reduced_energy_change = 0.0
    else:
        old_reduced = (
            old.energy_density_j_m3 / rho_dm_p
            - old.retained_water_loading * old_h_mass
        )
        new_reduced = (
            new.energy_density_j_m3 / rho_dm_p
            - new.retained_water_loading * new_h_mass
        )
        reduced_energy_change = new_reduced - old_reduced
    water_energy_change = math.fsum(
        (
            old.retained_water_loading * (new_h_mass - old_h_mass),
            (new.retained_water_loading - old.retained_water_loading) * new_h_mass,
        )
    )
    return volume_m3 * rho_dm_p * math.fsum(
        (reduced_energy_change, water_energy_change)
    )


def _residual_scales(
    grid: SphericalGrid,
    old_cells: Sequence[WetWaterCellState],
    boundary: WetWaterDirichletBoundary,
    model: WetWaterModel,
    dt_s: float,
) -> tuple[float, float]:
    """Fixed extensive scales keep the nonlinear residual smooth at zero flux."""

    initial_fluxes = _face_fluxes(old_cells, grid, boundary, model)
    water_face_amounts = tuple(
        dt_s * area * face.retained_water_flux_mol_m2_s
        for area, face in zip(grid.areas, initial_fluxes)
    )
    energy_face_amounts = tuple(
        dt_s * area * face.total_energy_flux_w_m2
        for area, face in zip(grid.areas, initial_fluxes)
    )
    water_inventory = math.fsum(
        volume * cell.retained_water_concentration_mol_m3
        for volume, cell in zip(grid.volumes, old_cells)
    )
    # A capacity scale is invariant to the arbitrary common caloric datum.
    # Scaling on |U| would make the nonlinear acceptance depend on an EOS
    # reference-state shift even though the balance equations do not.
    thermal_capacity_scale = math.fsum(
        volume
        * wet_core.energy_capacity(
            cell.temperature_k,
            cell.historical_hexane_loading,
            replace(
                model.wet,
                X_water=cell.retained_water_loading,
                w_o=cell.oil_fraction_label,
            ),
        )
        * max(abs(cell.temperature_k), 1.0)
        for volume, cell in zip(grid.volumes, old_cells)
    )
    return (
        max(
            *(abs(value) for value in water_face_amounts),
            1.0e-14 * abs(water_inventory),
            1.0e-30,
        ),
        max(
            *(abs(value) for value in energy_face_amounts),
            1.0e-14 * thermal_capacity_scale,
            1.0e-30,
        ),
    )


def _build_ledger(
    old_state: WetWaterState,
    new_state: WetWaterState,
    old_cells: Sequence[WetWaterCellState],
    new_cells: Sequence[WetWaterCellState],
    fluxes: Sequence[WetWaterFaceFlux],
    dt_s: float,
    *,
    max_local_relative_residual: float,
    nonlinear_evaluations: int,
) -> WetWaterLedger:
    grid = old_state.grid
    water_cell_changes = tuple(
        _water_amount_change(
            volume, old, new, old_state.model.wet.rho_dm_p
        )
        for volume, old, new in zip(grid.volumes, old_cells, new_cells)
    )
    energy_cell_changes = tuple(
        _energy_amount_change(
            volume, old, new, old_state.model.wet.rho_dm_p
        )
        for volume, old, new in zip(grid.volumes, old_cells, new_cells)
    )
    step_water_change = math.fsum(water_cell_changes)
    step_energy_change = math.fsum(energy_cell_changes)
    step_water_out = (
        dt_s * grid.areas[-1] * fluxes[-1].retained_water_flux_mol_m2_s
    )
    step_energy_out = dt_s * grid.areas[-1] * fluxes[-1].total_energy_flux_w_m2
    step_water_residual = step_water_change + step_water_out
    step_energy_residual = step_energy_change + step_energy_out
    step_water_relative = _relative_residual(
        step_water_residual,
        step_water_change,
        step_water_out,
        math.fsum(abs(value) for value in water_cell_changes),
    )
    step_energy_relative = _relative_residual(
        step_energy_residual,
        step_energy_change,
        step_energy_out,
        math.fsum(abs(value) for value in energy_cell_changes),
    )

    reference_cells = tuple(
        evaluate_cell(T, W, Xh, wo, new_state.model)
        for T, W, Xh, wo in zip(
            new_state.reference_temperatures_k,
            new_state.reference_retained_water_loadings,
            new_state.historical_hexane_loadings,
            new_state.oil_fraction_labels,
        )
    )
    reconstructed_reference_water = tuple(
        volume * cell.retained_water_concentration_mol_m3
        for volume, cell in zip(grid.volumes, reference_cells)
    )
    reconstructed_reference_energy = tuple(
        volume * cell.energy_density_j_m3
        for volume, cell in zip(grid.volumes, reference_cells)
    )
    if (
        reconstructed_reference_water != new_state.reference_water_cell_mol
        or reconstructed_reference_energy != new_state.reference_energy_cell_j
    ):
        raise RuntimeError("wet-water restart reference inventories are inconsistent")
    cumulative_water_change = math.fsum(
        _water_amount_change(
            volume,
            reference,
            current,
            new_state.model.wet.rho_dm_p,
        )
        for volume, reference, current in zip(
            grid.volumes, reference_cells, new_cells
        )
    )
    cumulative_energy_change = math.fsum(
        _energy_amount_change(
            volume,
            reference,
            current,
            new_state.model.wet.rho_dm_p,
        )
        for volume, reference, current in zip(
            grid.volumes, reference_cells, new_cells
        )
    )
    cumulative_water_residual = (
        cumulative_water_change + new_state.cumulative_boundary_water_out_mol
    )
    cumulative_energy_residual = (
        cumulative_energy_change + new_state.cumulative_boundary_energy_out_j
    )
    cumulative_water_relative = _relative_residual(
        cumulative_water_residual,
        cumulative_water_change,
        new_state.cumulative_boundary_water_out_mol,
        new_state.cumulative_absolute_water_exchange_mol,
    )
    cumulative_energy_relative = _relative_residual(
        cumulative_energy_residual,
        cumulative_energy_change,
        new_state.cumulative_boundary_energy_out_j,
        new_state.cumulative_absolute_energy_exchange_j,
    )
    current_hexane = math.fsum(
        volume * cell.total_hexane_concentration_kg_m3
        for volume, cell in zip(grid.volumes, new_cells)
    )
    ledger = WetWaterLedger(
        step_water_change_mol=step_water_change,
        step_water_boundary_out_mol=step_water_out,
        step_water_residual_mol=step_water_residual,
        step_water_relative_residual=step_water_relative,
        step_energy_change_j=step_energy_change,
        step_energy_boundary_out_j=step_energy_out,
        step_energy_residual_j=step_energy_residual,
        step_energy_relative_residual=step_energy_relative,
        cumulative_water_change_mol=cumulative_water_change,
        cumulative_water_boundary_out_mol=(
            new_state.cumulative_boundary_water_out_mol
        ),
        cumulative_water_residual_mol=cumulative_water_residual,
        cumulative_water_relative_residual=cumulative_water_relative,
        cumulative_energy_change_j=cumulative_energy_change,
        cumulative_energy_boundary_out_j=(
            new_state.cumulative_boundary_energy_out_j
        ),
        cumulative_energy_residual_j=cumulative_energy_residual,
        cumulative_energy_relative_residual=cumulative_energy_relative,
        max_local_relative_residual=max_local_relative_residual,
        historical_hexane_change_kg=current_hexane - new_state.reference_hexane_mass_kg,
        max_total_relative_hexane_flux_kg_m2_s=max(
            abs(face.total_relative_hexane_flux_kg_m2_s) for face in fluxes
        ),
        nonlinear_evaluations=nonlinear_evaluations,
        accepted=True,
    )
    if (
        ledger.historical_hexane_change_kg != 0.0
        or ledger.max_total_relative_hexane_flux_kg_m2_s != 0.0
    ):
        raise RuntimeError("immutable wet n-hexane ledger changed")
    return ledger


def _relative_residual(residual: float, *scales: float) -> float:
    denominator = max(*(abs(value) for value in scales), 1.0e-30)
    return abs(residual) / denominator


def _require_acceptable_ledger(ledger: WetWaterLedger, tolerance: float) -> None:
    residuals = (
        ledger.step_water_relative_residual,
        ledger.step_energy_relative_residual,
        ledger.cumulative_water_relative_residual,
        ledger.cumulative_energy_relative_residual,
        ledger.max_local_relative_residual,
    )
    if not all(math.isfinite(value) and value >= 0.0 for value in residuals):
        raise RuntimeError("wet-water conservation ledger is non-finite or negative")
    maximum = max(residuals)
    if maximum > tolerance:
        raise RuntimeError(
            "candidate state failed the <=1e-10 conservation contract: "
            f"{maximum:.3e}"
        )


def _encode_primitives(
    temperatures: Sequence[float], water: Sequence[float], model: WetWaterModel
) -> np.ndarray:
    encoded_temperature = [
        _logit_from_interval(value, model.wet.T_min, model.wet.T_max)
        for value in temperatures
    ]
    encoded_water = [
        _logit_from_interval(value, model.luikov.W_ref, model.luikov.W_cap)
        for value in water
    ]
    return np.asarray(encoded_temperature + encoded_water, dtype=float)


def _decode_primitives(
    coordinates: Sequence[float], model: WetWaterModel
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    count = len(coordinates) // 2
    if len(coordinates) != 2 * count:
        raise ValueError("wet-water primitive vector must contain 2N values")
    temperatures = tuple(
        _interval_from_logit(value, model.wet.T_min, model.wet.T_max)
        for value in coordinates[:count]
    )
    water = tuple(
        _interval_from_logit(value, model.luikov.W_ref, model.luikov.W_cap)
        for value in coordinates[count:]
    )
    return temperatures, water


def _logit_from_interval(value: float, lower: float, upper: float) -> float:
    if not lower < value < upper:
        raise ValueError("a driven smooth primitive must lie strictly inside its bounds")
    fraction = (value - lower) / (upper - lower)
    return math.log(fraction / (1.0 - fraction))


def _interval_from_logit(value: float, lower: float, upper: float) -> float:
    if not math.isfinite(value):
        raise ValueError("nonlinear coordinate must be finite")
    if value >= 0.0:
        exponential = math.exp(-value) if value < 745.0 else 0.0
        fraction = 1.0 / (1.0 + exponential)
    else:
        exponential = math.exp(value) if value > -745.0 else 0.0
        fraction = exponential / (1.0 + exponential)
    decoded = lower + (upper - lower) * fraction
    if not lower < decoded < upper:
        raise ValueError(
            "smooth nonlinear coordinate reached an exact active-set endpoint; "
            "reject rather than round or clip"
        )
    return decoded


def _validate_boundary(
    boundary: WetWaterDirichletBoundary, model: WetWaterModel
) -> None:
    if not isinstance(boundary, WetWaterDirichletBoundary):
        raise TypeError("wet-water predictor requires its explicit Dirichlet boundary")
    if not math.isfinite(boundary.temperature_k):
        raise ValueError("Dirichlet temperature must be finite")
    if not model.wet.T_min <= boundary.temperature_k <= model.wet.T_max:
        raise ValueError("Dirichlet temperature is outside the caloric bracket")
    sp.water_activity(boundary.retained_water_loading, model.luikov)
    wa.state_Tp(boundary.temperature_k, model.wet.pressure_pa, "liquid")


def _validate_state_structure(state: WetWaterState) -> None:
    n = state.grid.n
    for name, values in (
        ("temperatures", state.temperatures_k),
        ("retained-water loadings", state.retained_water_loadings),
        ("historical n-hexane labels", state.historical_hexane_loadings),
        ("oil labels", state.oil_fraction_labels),
        ("reference temperatures", state.reference_temperatures_k),
        (
            "reference retained-water loadings",
            state.reference_retained_water_loadings,
        ),
        ("reference water inventory", state.reference_water_cell_mol),
        ("reference energy inventory", state.reference_energy_cell_j),
    ):
        if len(values) != n:
            raise ValueError(f"{name} has length {len(values)}, expected {n}")
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"{name} must contain only finite values")
    scalar_values = (
        state.time_s,
        state.reference_hexane_mass_kg,
        state.cumulative_boundary_water_out_mol,
        state.cumulative_boundary_energy_out_j,
        state.cumulative_absolute_water_exchange_mol,
        state.cumulative_absolute_energy_exchange_j,
    )
    if not all(math.isfinite(value) for value in scalar_values):
        raise ValueError("wet-water restart metadata must be finite")
    if state.time_s < 0.0:
        raise ValueError("wet-water restart time must be non-negative")
    if any(value < 0.0 for value in state.reference_water_cell_mol):
        raise ValueError("reference retained-water inventories must be non-negative")
    if state.reference_hexane_mass_kg < 0.0:
        raise ValueError("reference n-hexane inventory must be non-negative")
    if (
        state.cumulative_absolute_water_exchange_mol < 0.0
        or state.cumulative_absolute_energy_exchange_j < 0.0
    ):
        raise ValueError("cumulative absolute exchanges must be non-negative")
    if state.physically_qualifying:
        raise ValueError("wet-water foundation cannot be physically qualifying")


def _is_uniform_invariant(
    state: WetWaterState, boundary: WetWaterDirichletBoundary
) -> bool:
    return (
        all(value == boundary.temperature_k for value in state.temperatures_k)
        and all(
            value == boundary.retained_water_loading
            for value in state.retained_water_loadings
        )
    )


def _tuple_of_length(
    name: str, values: Sequence[float], expected: int
) -> tuple[float, ...]:
    result = tuple(values)
    if len(result) != expected:
        raise ValueError(f"{name} has length {len(result)}, expected {expected}")
    return result


__all__ = [
    "WetWaterCellState",
    "WetWaterDirichletBoundary",
    "WetWaterFaceFlux",
    "WetWaterLedger",
    "WetWaterModel",
    "WetWaterState",
    "WetWaterStepError",
    "WetWaterStepResult",
    "clear_evaluate_cell_cache",
    "evaluate_cell",
    "evaluate_cell_cache_info",
    "initialize",
    "step",
]
