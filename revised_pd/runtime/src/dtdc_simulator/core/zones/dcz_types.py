"""Data contracts and acceptance tolerances for the DCZ solver.

This module contains state only.  Governing equations and nonlinear
algorithms belong in :mod:`dtdc_simulator.core.zones.dcz` and its numerical
submodules.  Keeping these immutable contracts separate makes the scientific
API reviewable without coupling it to a particular solution algorithm.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from dtdc_simulator.core import thermo
from dtdc_simulator.core.zones import particle as pt

# Physical fixed-point and accepted-ledger tolerances.
EXIT_TOL_T = 1.0e-2
EXIT_TOL_X1 = 1.0e-4
EXIT_TOL_X2 = 1.0e-6
EXIT_TOL_FLOW = 1.0e-4
BALANCE_TOL_FLOW = 1.0e-9
ENERGY_BALANCE_TOL_W = 1.0
PHASE_EQUILIBRIUM_TOL_K = 1.0e-3
PROFILE_TOL_T = 2.0e-2
PROFILE_TOL_FLOW = 2.0e-3


@dataclass(frozen=True)
class DCZResiduals:
    """Final fixed-point residuals for every state exported by the DCZ map."""

    solid_out_X2: float = math.inf
    solid_out_T: float = math.inf
    solid_out_X1: float = math.inf
    vapor_top_T: float = math.inf
    vapor_top_water_flow: float = math.inf
    vapor_top_hexane_flow: float = math.inf
    vapor_profile_T_max: float = math.inf
    vapor_profile_water_flow_max: float = math.inf
    vapor_profile_hexane_flow_max: float = math.inf
    evaporated_water_flow_max: float = math.inf
    free_water_max: float = math.inf
    water_boundary_balance: float = math.inf
    vapor_energy_balance_max: float = math.inf
    interphase_energy_mismatch_max: float = math.inf
    phase_equilibrium_temperature_max: float = math.inf

    @property
    def maximum_scaled(self) -> float:
        """Largest residual normalized by its physical acceptance tolerance."""
        return max(
            self.solid_out_X2 / EXIT_TOL_X2,
            self.solid_out_T / EXIT_TOL_T,
            self.solid_out_X1 / EXIT_TOL_X1,
            self.vapor_top_T / PROFILE_TOL_T,
            self.vapor_top_water_flow / EXIT_TOL_FLOW,
            self.vapor_top_hexane_flow / EXIT_TOL_FLOW,
            self.vapor_profile_T_max / PROFILE_TOL_T,
            self.vapor_profile_water_flow_max / PROFILE_TOL_FLOW,
            self.vapor_profile_hexane_flow_max / PROFILE_TOL_FLOW,
            self.evaporated_water_flow_max / PROFILE_TOL_FLOW,
            self.free_water_max / EXIT_TOL_X1,
            self.water_boundary_balance / BALANCE_TOL_FLOW,
            self.vapor_energy_balance_max / ENERGY_BALANCE_TOL_W,
            self.interphase_energy_mismatch_max / ENERGY_BALANCE_TOL_W,
            self.phase_equilibrium_temperature_max / PHASE_EQUILIBRIUM_TOL_K,
        )


@dataclass(frozen=True)
class DCZWarmStart:
    """Complete lagged state of the nested DCZ fixed-point iteration."""

    bed_height_m: float
    particles: tuple[pt.ParticleState, ...]
    vapor_T: tuple[float, ...]
    vapor_water_kg_s: tuple[float, ...]
    vapor_hexane_kg_s: tuple[float, ...]
    X1: tuple[float, ...]
    free_water: tuple[float, ...]
    condensed_water_kg_s: tuple[float, ...]
    evaporated_water_kg_s: tuple[float, ...]
    bound_water_mass_rate_w_m3: tuple[float, ...]
    dwpg2_dt: tuple[tuple[float, ...], ...]
    temperature_relaxation: float = 0.5
    hexane_relaxation: float = 0.5
    water_relaxation: float = 0.1
    water_active_set: tuple[tuple[bool, bool, bool], ...] | None = None
    cell_widths_m: tuple[float, ...] = ()


@dataclass(frozen=True)
class DCZConstants:
    """Physical and constitutive constants consumed by a DCZ zone solve."""

    diameter_m: float
    bed_height_m: float
    hM: float
    hQ: float
    aV: float
    D_ax: float
    k_mixL: float
    rho_V: float
    cp_V: float
    alpha_V: float
    alpha_L: float
    particle: pt.ParticleConstants
    dH_vap_water: float
    antoine_water: thermo.AntoineParams
    antoine_hexane: thermo.AntoineParams
    luikov: thermo.LuikovParams
    water_diffusivity: float
    vapor_enthalpy_ref: thermo.VaporEnthalpyRef
    pressure_pa: float = thermo.ATM_PRESSURE_PA


@dataclass(frozen=True)
class VaporState:
    """Vapor composition and temperature at one DCZ face."""

    wV2: float
    T: float


@dataclass(frozen=True)
class DCZCellResult:
    """Accepted physical state and exact-ledger telemetry for one DCZ cell."""

    vapor_top: VaporState
    particle: pt.ParticleState
    X2_bulk: float
    X1_bulk: float
    free_water_bulk: float
    condensed_water_kg_s: float
    vapor_water_kg_s: float
    vapor_hexane_kg_s: float
    evaporated_water_kg_s: float = 0.0
    water_evaporation_energy_w: float = 0.0
    vapor_energy_balance_residual_w: float = 0.0
    particle_energy_accumulation_w: float = 0.0
    particle_sorption_source_w: float = 0.0
    particle_axial_conduction_source_w: float = 0.0
    particle_surface_heat_to_vapor_w: float = 0.0
    particle_energy_balance_residual_w: float = 0.0
    interphase_energy_mismatch_w: float = 0.0
    combined_energy_balance_residual_w: float = 0.0
    wet_phase_constraint_active: bool = False
    wet_phase_water_limited: bool = False
    evaporative_pin_energy_removed_w: float = 0.0
    particle_temperature_clamped: bool = False

    @property
    def vapor_flow_kg_s(self) -> float:
        return self.vapor_water_kg_s + self.vapor_hexane_kg_s


@dataclass(frozen=True)
class DCZZoneResult:
    """Accepted DCZ profile and convergence metadata."""

    cells: tuple[DCZCellResult, ...]
    iterations: int
    converged: bool = True
    residuals: DCZResiduals = DCZResiduals()
    warm_start: DCZWarmStart | None = None
    coupling_tolerance: float | None = None
    cell_widths_m: tuple[float, ...] = ()

    @property
    def vapor_out(self) -> VaporState:
        return self.cells[0].vapor_top

    @property
    def solid_out_X2(self) -> float:
        return self.cells[-1].X2_bulk

    @property
    def solid_out_X1(self) -> float:
        return self.cells[-1].X1_bulk

    @property
    def total_condensed_kg_s(self) -> float:
        return sum(cell.condensed_water_kg_s for cell in self.cells)

    @property
    def vapor_water_out_kg_s(self) -> float:
        return self.cells[0].vapor_water_kg_s

    @property
    def vapor_hexane_out_kg_s(self) -> float:
        return self.cells[0].vapor_hexane_kg_s

    @property
    def vapor_flow_out_kg_s(self) -> float:
        return self.cells[0].vapor_flow_kg_s

    @property
    def total_evaporative_pin_energy_removed_w(self) -> float:
        return sum(cell.evaporative_pin_energy_removed_w for cell in self.cells)

    @property
    def evaporative_pin_cell_count(self) -> int:
        return sum(cell.evaporative_pin_energy_removed_w > 0.0 for cell in self.cells)

    @property
    def total_evaporated_water_kg_s(self) -> float:
        return sum(cell.evaporated_water_kg_s for cell in self.cells)

    @property
    def total_water_evaporation_energy_w(self) -> float:
        return sum(cell.water_evaporation_energy_w for cell in self.cells)

    @property
    def total_vapor_energy_balance_residual_w(self) -> float:
        return sum(cell.vapor_energy_balance_residual_w for cell in self.cells)

    @property
    def maximum_absolute_vapor_energy_balance_residual_w(self) -> float:
        return max(abs(cell.vapor_energy_balance_residual_w) for cell in self.cells)

    @property
    def total_combined_energy_balance_residual_w(self) -> float:
        return sum(cell.combined_energy_balance_residual_w for cell in self.cells)

    @property
    def maximum_absolute_combined_energy_balance_residual_w(self) -> float:
        return max(abs(cell.combined_energy_balance_residual_w) for cell in self.cells)

    @property
    def wet_phase_constraint_cell_count(self) -> int:
        return sum(cell.wet_phase_constraint_active for cell in self.cells)

    @property
    def wet_phase_water_limited_cell_count(self) -> int:
        return sum(cell.wet_phase_water_limited for cell in self.cells)

    @property
    def simultaneous_condensation_evaporation_cell_count(self) -> int:
        return sum(
            cell.condensed_water_kg_s > 1.0e-12
            and cell.evaporated_water_kg_s > 1.0e-12
            for cell in self.cells
        )

    @property
    def particle_temperature_clamp_cell_count(self) -> int:
        return sum(cell.particle_temperature_clamped for cell in self.cells)
