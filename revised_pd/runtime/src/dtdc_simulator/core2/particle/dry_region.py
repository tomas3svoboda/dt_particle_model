r"""Conservative finite-volume ledgers for a moving spherical dry region.

This module is deliberately constitutive-free.  A caller supplies the physical
component fluxes ``j_k``, heat flux ``q``, and transported partial enthalpies
``h_k`` at the new-time interior faces.  For each face the kernel forms the
exact relative ALE rates

    Q_f   = 4*pi/3 * (r_f,new**3 - r_f,old**3) / dt,
    H_k,f = A_f,new*j_k,f - Q_f*C_k,upwind,
    H_E,f = A_f,new*(q_f + sum_k h_k,f*j_k,f) - Q_f*u_upwind.

Both boundary inputs are already relative extensive rates ``H`` in the
increasing-r direction.  Consequently cell ``i`` has the backward-Euler
residual

    V_i,new*y_i,new - V_i,old*y_i,old + dt*(H_right-H_left),

evaluated in the cancellation-resistant equivalent form

    V_i,old*(y_i,new-y_i,old) + dt*(Q_right-Q_left)*y_i,new
        + dt*(H_right-H_left).

and summing the cells cancels every interior rate exactly.  Keeping the inner
and outer relative rates explicit is essential at the wet/dry interface: the
front jump supplies the former and the particle surface closure supplies the
latter, without either being reinterpreted as a laboratory-frame flux.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from dtdc_simulator.core2.particle.grid import SphericalGrid, spherical_volume_difference


class DryRegionTopologyError(ValueError):
    """The supplied grids do not describe admissible primary drainage."""


@dataclass(frozen=True)
class DryRegionState:
    """Cell state with component-major concentrations and one energy field."""

    components: tuple[str, ...]
    component_concentrations_kg_m3: tuple[tuple[float, ...], ...]
    energy_density_j_m3: tuple[float, ...]

    def __post_init__(self) -> None:
        _validate_components(self.components)
        if len(self.component_concentrations_kg_m3) != len(self.components):
            raise ValueError("dry-region concentration fields must align with components")
        cell_count = len(self.energy_density_j_m3)
        if cell_count < 1:
            raise ValueError("dry-region state needs at least one cell")
        if any(len(field) != cell_count for field in self.component_concentrations_kg_m3):
            raise ValueError("dry-region state fields must have one value per cell")
        concentrations = (value for field in self.component_concentrations_kg_m3 for value in field)
        if not all(math.isfinite(value) and value >= 0.0 for value in concentrations):
            raise ValueError("component concentrations must be finite and non-negative")
        if not all(math.isfinite(value) for value in self.energy_density_j_m3):
            raise ValueError("energy densities must be finite")

    @property
    def n(self) -> int:
        return len(self.energy_density_j_m3)


@dataclass(frozen=True)
class InteriorFluxes:
    """Laboratory-frame physical flux data on all new-time interior faces.

    Component fluxes and partial enthalpies are component-major.  Heat fluxes
    are face-major.  Every flux is outward in the increasing-r direction.
    """

    components: tuple[str, ...]
    component_fluxes_kg_m2_s: tuple[tuple[float, ...], ...]
    partial_enthalpies_j_kg: tuple[tuple[float, ...], ...]
    heat_fluxes_w_m2: tuple[float, ...]

    def __post_init__(self) -> None:
        _validate_components(self.components)
        component_count = len(self.components)
        if (
            len(self.component_fluxes_kg_m2_s) != component_count
            or len(self.partial_enthalpies_j_kg) != component_count
        ):
            raise ValueError("interior flux fields must align with components")
        face_count = len(self.heat_fluxes_w_m2)
        fields = (*self.component_fluxes_kg_m2_s, *self.partial_enthalpies_j_kg)
        if any(len(field) != face_count for field in fields):
            raise ValueError("interior flux fields must have one value per face")
        values = (value for field in fields for value in field)
        if not all(math.isfinite(value) for value in values) or not all(
            math.isfinite(value) for value in self.heat_fluxes_w_m2
        ):
            raise ValueError("interior flux data must be finite")

    @property
    def face_count(self) -> int:
        return len(self.heat_fluxes_w_m2)


@dataclass(frozen=True)
class RelativeBoundaryRates:
    """Prescribed extensive ALE rates ``H`` at one dry-region boundary."""

    components: tuple[str, ...]
    component_rates_kg_s: tuple[float, ...]
    energy_rate_w: float

    def __post_init__(self) -> None:
        _validate_components(self.components)
        if len(self.component_rates_kg_s) != len(self.components):
            raise ValueError("boundary component rates must align with components")
        if not all(
            math.isfinite(value) for value in (*self.component_rates_kg_s, self.energy_rate_w)
        ):
            raise ValueError("boundary relative rates must be finite")


@dataclass(frozen=True)
class RelativeFaceRates:
    """One interior face's relative component and total-energy rates."""

    components: tuple[str, ...]
    component_rates_kg_s: tuple[float, ...]
    energy_rate_w: float
    swept_volume_rate_m3_s: float
    upwind_side: str


@dataclass(frozen=True)
class CellResidualLedger:
    """Backward-Euler cell residuals and the interior rates that formed them."""

    components: tuple[str, ...]
    component_residuals_kg: tuple[tuple[float, ...], ...]
    energy_residuals_j: tuple[float, ...]
    swept_volume_rates_m3_s: tuple[float, ...]
    interior_relative_rates: tuple[RelativeFaceRates, ...]

    @property
    def component_cell_sums_kg(self) -> tuple[float, ...]:
        return tuple(math.fsum(field) for field in self.component_residuals_kg)

    @property
    def energy_cell_sum_j(self) -> float:
        return math.fsum(self.energy_residuals_j)


@dataclass(frozen=True)
class TotalResidualLedger:
    """Whole-annulus residual computed independently of interior face rates."""

    components: tuple[str, ...]
    component_inventory_changes_kg: tuple[float, ...]
    component_boundary_terms_kg: tuple[float, ...]
    component_residuals_kg: tuple[float, ...]
    energy_inventory_change_j: float
    energy_boundary_term_j: float
    energy_residual_j: float


def swept_volume_rates(
    grid_prev: SphericalGrid, grid: SphericalGrid, dt: float
) -> tuple[float, ...]:
    """Return exact swept-volume rates at every corresponding grid face."""
    _validate_grid_pair(grid_prev, grid, dt)
    return tuple(
        spherical_volume_difference(r_old, r_new) / dt
        for r_old, r_new in zip(grid_prev.faces, grid.faces)
    )


def relative_face_rates(
    *,
    components: Sequence[str],
    area_m2: float,
    swept_volume_rate_m3_s: float,
    component_fluxes_kg_m2_s: Sequence[float],
    partial_enthalpies_j_kg: Sequence[float],
    heat_flux_w_m2: float,
    left_concentrations_kg_m3: Sequence[float],
    right_concentrations_kg_m3: Sequence[float],
    left_energy_density_j_m3: float,
    right_energy_density_j_m3: float,
) -> RelativeFaceRates:
    """Form ``A*j-Q*C_upwind`` and ``A*(q+sum(h*j))-Q*u_upwind``.

    If ``Q < 0`` the mesh face moves inward and the relative transport is
    outward, so the left state is upwind.  If ``Q > 0`` the right state is
    upwind.  At ``Q == 0`` the selected state is immaterial; ``left`` is
    reported for deterministic diagnostics.
    """
    component_tuple = tuple(components)
    _validate_components(component_tuple)
    vectors = (
        component_fluxes_kg_m2_s,
        partial_enthalpies_j_kg,
        left_concentrations_kg_m3,
        right_concentrations_kg_m3,
    )
    if any(len(values) != len(component_tuple) for values in vectors):
        raise ValueError("relative face vectors must align with components")
    scalar_values = (
        area_m2,
        swept_volume_rate_m3_s,
        heat_flux_w_m2,
        left_energy_density_j_m3,
        right_energy_density_j_m3,
    )
    vector_values = (value for values in vectors for value in values)
    if not all(math.isfinite(value) for value in scalar_values) or not all(
        math.isfinite(value) for value in vector_values
    ):
        raise ValueError("relative face inputs must be finite")
    if area_m2 < 0.0:
        raise ValueError("face area must be non-negative")
    if any(value < 0.0 for value in left_concentrations_kg_m3) or any(
        value < 0.0 for value in right_concentrations_kg_m3
    ):
        raise ValueError("face concentrations must be non-negative")
    if swept_volume_rate_m3_s <= 0.0:
        upwind_concentrations = left_concentrations_kg_m3
        upwind_energy = left_energy_density_j_m3
        upwind_side = "left"
    else:
        upwind_concentrations = right_concentrations_kg_m3
        upwind_energy = right_energy_density_j_m3
        upwind_side = "right"
    component_rates = tuple(
        area_m2 * flux - swept_volume_rate_m3_s * concentration
        for flux, concentration in zip(component_fluxes_kg_m2_s, upwind_concentrations)
    )
    laboratory_energy_flux = math.fsum(
        (
            heat_flux_w_m2,
            *(
                enthalpy * flux
                for enthalpy, flux in zip(partial_enthalpies_j_kg, component_fluxes_kg_m2_s)
            ),
        )
    )
    energy_rate = area_m2 * laboratory_energy_flux - swept_volume_rate_m3_s * upwind_energy
    return RelativeFaceRates(
        components=component_tuple,
        component_rates_kg_s=component_rates,
        energy_rate_w=energy_rate,
        swept_volume_rate_m3_s=swept_volume_rate_m3_s,
        upwind_side=upwind_side,
    )


def backward_euler_residuals(
    grid_prev: SphericalGrid,
    grid: SphericalGrid,
    state_prev: DryRegionState,
    state: DryRegionState,
    dt: float,
    interior_fluxes: InteriorFluxes,
    inner_boundary: RelativeBoundaryRates,
    outer_boundary: RelativeBoundaryRates,
) -> CellResidualLedger:
    """Assemble conservative cell residuals on a moving dry annulus."""
    swept = swept_volume_rates(grid_prev, grid, dt)
    _validate_ledger_inputs(
        grid_prev,
        state_prev,
        state,
        interior_fluxes,
        inner_boundary,
        outer_boundary,
    )
    internal_rates = tuple(
        relative_face_rates(
            components=state.components,
            area_m2=grid.areas[face],
            swept_volume_rate_m3_s=swept[face],
            component_fluxes_kg_m2_s=tuple(
                field[face - 1] for field in interior_fluxes.component_fluxes_kg_m2_s
            ),
            partial_enthalpies_j_kg=tuple(
                field[face - 1] for field in interior_fluxes.partial_enthalpies_j_kg
            ),
            heat_flux_w_m2=interior_fluxes.heat_fluxes_w_m2[face - 1],
            left_concentrations_kg_m3=tuple(
                field[face - 1] for field in state.component_concentrations_kg_m3
            ),
            right_concentrations_kg_m3=tuple(
                field[face] for field in state.component_concentrations_kg_m3
            ),
            left_energy_density_j_m3=state.energy_density_j_m3[face - 1],
            right_energy_density_j_m3=state.energy_density_j_m3[face],
        )
        for face in range(1, grid.n)
    )
    component_faces = tuple(
        (inner_boundary.component_rates_kg_s[index],)
        + tuple(rate.component_rates_kg_s[index] for rate in internal_rates)
        + (outer_boundary.component_rates_kg_s[index],)
        for index in range(len(state.components))
    )
    energy_faces = (
        (inner_boundary.energy_rate_w,)
        + tuple(rate.energy_rate_w for rate in internal_rates)
        + (outer_boundary.energy_rate_w,)
    )
    volume_changes = _stable_cell_volume_changes(swept, dt)
    component_residuals = tuple(
        tuple(
            math.fsum(
                (
                    grid_prev.volumes[cell] * (current[cell] - previous[cell]),
                    volume_changes[cell] * current[cell],
                    dt * (face_rates[cell + 1] - face_rates[cell]),
                )
            )
            for cell in range(grid.n)
        )
        for previous, current, face_rates in zip(
            state_prev.component_concentrations_kg_m3,
            state.component_concentrations_kg_m3,
            component_faces,
        )
    )
    energy_residuals = tuple(
        math.fsum(
            (
                grid_prev.volumes[cell]
                * (state.energy_density_j_m3[cell] - state_prev.energy_density_j_m3[cell]),
                volume_changes[cell] * state.energy_density_j_m3[cell],
                dt * (energy_faces[cell + 1] - energy_faces[cell]),
            )
        )
        for cell in range(grid.n)
    )
    return CellResidualLedger(
        components=state.components,
        component_residuals_kg=component_residuals,
        energy_residuals_j=energy_residuals,
        swept_volume_rates_m3_s=swept,
        interior_relative_rates=internal_rates,
    )


def telescoping_totals(
    grid_prev: SphericalGrid,
    grid: SphericalGrid,
    state_prev: DryRegionState,
    state: DryRegionState,
    dt: float,
    inner_boundary: RelativeBoundaryRates,
    outer_boundary: RelativeBoundaryRates,
) -> TotalResidualLedger:
    """Compute the whole-region ledger without using any interior face rate."""
    _validate_grid_pair(grid_prev, grid, dt)
    _validate_state_and_boundaries(grid_prev, state_prev, state, inner_boundary, outer_boundary)
    swept = swept_volume_rates(grid_prev, grid, dt)
    volume_changes = _stable_cell_volume_changes(swept, dt)
    inventory_changes = tuple(
        math.fsum(
            math.fsum(
                (
                    old_volume * (new_value - old_value),
                    volume_change * new_value,
                )
            )
            for old_volume, volume_change, old_value, new_value in zip(
                grid_prev.volumes, volume_changes, previous, current
            )
        )
        for previous, current in zip(
            state_prev.component_concentrations_kg_m3,
            state.component_concentrations_kg_m3,
        )
    )
    boundary_terms = tuple(
        dt * (outer - inner)
        for inner, outer in zip(
            inner_boundary.component_rates_kg_s,
            outer_boundary.component_rates_kg_s,
        )
    )
    component_residuals = tuple(
        inventory + boundary for inventory, boundary in zip(inventory_changes, boundary_terms)
    )
    energy_change = math.fsum(
        math.fsum(
            (
                old_volume * (new_value - old_value),
                volume_change * new_value,
            )
        )
        for old_volume, volume_change, old_value, new_value in zip(
            grid_prev.volumes,
            volume_changes,
            state_prev.energy_density_j_m3,
            state.energy_density_j_m3,
        )
    )
    energy_boundary = dt * (outer_boundary.energy_rate_w - inner_boundary.energy_rate_w)
    return TotalResidualLedger(
        components=state.components,
        component_inventory_changes_kg=inventory_changes,
        component_boundary_terms_kg=boundary_terms,
        component_residuals_kg=component_residuals,
        energy_inventory_change_j=energy_change,
        energy_boundary_term_j=energy_boundary,
        energy_residual_j=energy_change + energy_boundary,
    )


def _stable_cell_volume_changes(
    swept_volume_rates_m3_s: Sequence[float], dt: float
) -> tuple[float, ...]:
    """Cell-volume increments from the authoritative face sweeps.

    This is algebraically ``V_new-V_old`` but remains resolved for a microscopic
    front move on an already-large annulus, where subtracting the two cell
    volumes loses the entire increment.
    """
    return tuple(
        dt * math.fsum((right, -left))
        for left, right in zip(swept_volume_rates_m3_s, swept_volume_rates_m3_s[1:])
    )


def _validate_components(components: Sequence[str]) -> None:
    if not components or any(not isinstance(name, str) or not name for name in components):
        raise ValueError("component names must be non-empty strings")
    if len(set(components)) != len(components):
        raise ValueError("component names must be unique")


def _validate_grid_pair(grid_prev: SphericalGrid, grid: SphericalGrid, dt: float) -> None:
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("dry-region time step must be positive and finite")
    if grid_prev.n != grid.n:
        raise DryRegionTopologyError("dry-region ALE step cannot change cell count")
    if grid_prev.R != grid.R:
        raise DryRegionTopologyError("dry-region ALE step cannot change particle radius")
    if grid.inner_radius > grid_prev.inner_radius:
        raise DryRegionTopologyError("primary-drainage wet front cannot move outward")
    _validate_grid_metrics(grid_prev)
    _validate_grid_metrics(grid)


def _validate_grid_metrics(grid: SphericalGrid) -> None:
    if grid.n < 1 or len(grid.faces) != grid.n + 1 or len(grid.areas) != grid.n + 1:
        raise DryRegionTopologyError("dry-region grid arrays are not aligned")
    if len(grid.volumes) != grid.n:
        raise DryRegionTopologyError("dry-region grid volumes are not aligned")
    if grid.inner_radius < 0.0 or grid.inner_radius >= grid.R:
        raise DryRegionTopologyError("dry region must have positive exact volume")
    if grid.faces[-1] != grid.R or any(
        right <= left for left, right in zip(grid.faces, grid.faces[1:])
    ):
        raise DryRegionTopologyError("dry-region grid faces have invalid topology")
    expected_volumes = (
        spherical_volume_difference(left, right) for left, right in zip(grid.faces, grid.faces[1:])
    )
    if any(
        not math.isclose(actual, expected, rel_tol=2.0e-15, abs_tol=0.0)
        for actual, expected in zip(grid.volumes, expected_volumes)
    ):
        raise DryRegionTopologyError("dry-region cell volumes are not exact spherical measures")
    if any(
        not math.isclose(area, 4.0 * math.pi * radius**2, rel_tol=2.0e-15, abs_tol=0.0)
        for area, radius in zip(grid.areas, grid.faces)
    ):
        raise DryRegionTopologyError("dry-region face areas are not exact spherical measures")


def _validate_state_and_boundaries(
    grid_prev: SphericalGrid,
    state_prev: DryRegionState,
    state: DryRegionState,
    inner_boundary: RelativeBoundaryRates,
    outer_boundary: RelativeBoundaryRates,
) -> None:
    if state_prev.n != grid_prev.n or state.n != grid_prev.n:
        raise ValueError("dry-region states must align with both grids")
    objects = (state_prev, state, inner_boundary, outer_boundary)
    if any(obj.components != state.components for obj in objects):
        raise ValueError("dry-region states and boundaries must use identical components")


def _validate_ledger_inputs(
    grid_prev: SphericalGrid,
    state_prev: DryRegionState,
    state: DryRegionState,
    interior_fluxes: InteriorFluxes,
    inner_boundary: RelativeBoundaryRates,
    outer_boundary: RelativeBoundaryRates,
) -> None:
    _validate_state_and_boundaries(grid_prev, state_prev, state, inner_boundary, outer_boundary)
    if interior_fluxes.components != state.components:
        raise ValueError("dry-region fluxes and states must use identical components")
    if interior_fluxes.face_count != grid_prev.n - 1:
        raise ValueError("dry-region fluxes must have one value per interior face")
