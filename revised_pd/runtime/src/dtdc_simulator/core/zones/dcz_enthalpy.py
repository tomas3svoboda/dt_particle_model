"""Conservative vapor-cell enthalpy equations for the DCZ."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from dtdc_simulator.core import thermo
from dtdc_simulator.core.critical_volume import (
    CriticalFixedPhaseOccupancy,
    CriticalPhaseLedger,
    CriticalVolumeBasis,
    DrySideFrontInventory,
    FrontMobileJumpLedger,
    solve_critical_phase_ledger,
    solve_front_mobile_jump_ledger,
)
from dtdc_simulator.core.zones.dcz_types import DCZConstants
from dtdc_simulator.core.zones import particle as pt


def vapor_specific_enthalpies_j_kg(
    temperature_k: float,
    ref: thermo.VaporEnthalpyRef,
) -> tuple[float, float]:
    """Water and hexane vapor enthalpies on one common normal-boiling datum."""

    return (
        ref.dH_vap_water
        + ref.cp_water_vapor * (temperature_k - ref.T_boil_water),
        ref.dH_vap_hexane
        + ref.cp_hexane_vapor * (temperature_k - ref.T_boil_hexane),
    )


def vapor_enthalpy_j(
    water_kg: float,
    hexane_kg: float,
    temperature_k: float,
    ref: thermo.VaporEnthalpyRef,
) -> float:
    h_water, h_hexane = vapor_specific_enthalpies_j_kg(temperature_k, ref)
    return water_kg * h_water + hexane_kg * h_hexane


def vapor_temperature_from_enthalpy_k(
    water_kg: float,
    hexane_kg: float,
    enthalpy_j: float,
    ref: thermo.VaporEnthalpyRef,
) -> float:
    """Invert the affine binary-vapor enthalpy without a numerical root."""

    capacity = (
        water_kg * ref.cp_water_vapor
        + hexane_kg * ref.cp_hexane_vapor
    )
    if water_kg < 0.0 or hexane_kg < 0.0 or capacity <= 0.0:
        raise ValueError("vapor holdup must be non-negative and nonempty")
    intercept = water_kg * (
        ref.dH_vap_water - ref.cp_water_vapor * ref.T_boil_water
    ) + hexane_kg * (
        ref.dH_vap_hexane - ref.cp_hexane_vapor * ref.T_boil_hexane
    )
    temperature = (enthalpy_j - intercept) / capacity
    if temperature <= 0.0 or not np.isfinite(temperature):
        raise ValueError("vapor enthalpy implies an invalid temperature")
    return float(temperature)


def adiabatic_vapor_mix_temperature_k(
    streams: tuple[tuple[float, float, float], ...],
    ref: thermo.VaporEnthalpyRef,
) -> float:
    """Mix ``(water kg/s, hexane kg/s, T K)`` streams on total enthalpy."""

    if not streams:
        raise ValueError("at least one vapor stream is required")
    water_total = 0.0
    hexane_total = 0.0
    enthalpy_total = 0.0
    for water_kg_s, hexane_kg_s, temperature_k in streams:
        if water_kg_s < 0.0 or hexane_kg_s < 0.0 or temperature_k <= 0.0:
            raise ValueError("vapor streams require non-negative flow and positive temperature")
        water_total += water_kg_s
        hexane_total += hexane_kg_s
        enthalpy_total += vapor_enthalpy_j(
            water_kg_s,
            hexane_kg_s,
            temperature_k,
            ref,
        )
    return vapor_temperature_from_enthalpy_k(
        water_total,
        hexane_total,
        enthalpy_total,
        ref,
    )


def particle_hexane_inventory_density_kg_m3(
    activity: float,
    temperature_k: float,
    c: pt.ParticleConstants,
) -> float:
    """Exact A.22 pore-plus-sorbed/oil hexane inventory density."""

    x2_so, _ = thermo.x2_so_and_slope(
        activity,
        temperature_k,
        c.X3,
        c.gab,
        c.oil,
    )
    return (
        c.alpha_pg * c.rho_pg * activity
        + c.alpha_ps * c.rho_ps * x2_so
    )


@dataclass(frozen=True)
class ParticleConstitutiveState:
    """Intensive particle state recovered from conserved shell variables."""

    activity: float
    temperature_k: float
    saturation_reserve_kg_m3: float


def particle_hexane_saturation_density_kg_m3(
    temperature_k: float,
    c: pt.ParticleConstants,
) -> float:
    """Largest equilibrium A.22 inventory at unit hexane activity."""

    return particle_hexane_inventory_density_kg_m3(1.0, temperature_k, c)


def legacy_geometric_pore_reserve_capacity_kg_m3(
    c: pt.ParticleConstants,
) -> float:
    """Return Coletto's raw unswollen-pore reserve oracle.

    This is deliberately *not* the production moving-front capacity.  It is
    retained as the mandatory ``kappa_v=1`` geometric limiting oracle.
    """

    return c.alpha_pg * max(
        c.rho_hexane_liquid_ref - c.rho_pg,
        0.0,
    )


def particle_wet_critical_phase_ledger(
    temperature_k: float,
    c: pt.ParticleConstants,
) -> CriticalPhaseLedger:
    """Return the approved dry-clean wet-side critical phase ledger.

    At the wet side of the front, retained and mobile n-hexane occupy the
    effective critical correlation volume with liquid-like partial volume.
    Pore gas belongs to the *dry* side of the jump and is therefore zero
    here.  Free water is external to this dry-clean solvent-front closure;
    a future co-located-water formulation must enter explicitly rather than
    silently reducing this volume.
    """

    if temperature_k <= 0.0 or not math.isfinite(temperature_k):
        raise ValueError("critical-ledger temperature must be positive and finite")
    dry_density = c.alpha_ps * c.rho_ps
    retained_loading, _ = thermo.x2_so_and_slope(
        1.0,
        temperature_k,
        c.X3,
        c.gab,
        c.oil,
    )
    liquid_density = thermo.rho_hexane_liquid(temperature_k)
    basis = CriticalVolumeBasis(
        particle_porosity=c.alpha_pg,
        dry_meal_envelope_density_kg_m3=dry_density,
        critical_volume_multiplier=c.critical_volume_multiplier,
    )
    return solve_critical_phase_ledger(
        basis,
        CriticalFixedPhaseOccupancy(
            retained_hexane_loading_kg_per_kg_dry_meal=retained_loading,
            pore_gas_hexane_loading_kg_per_kg_dry_meal=0.0,
            free_water_loading_kg_per_kg_dry_meal=0.0,
            pore_gas_specific_volume_m3_per_kg_dry_meal=0.0,
            liquid_hexane_density_kg_m3=liquid_density,
            # Inactive in the approved dry-clean ledger, but kept explicit
            # by the generic phase-ledger contract.
            liquid_water_density_kg_m3=1000.0,
        ),
    )


def particle_critical_hexane_loading_kg_kg_dry(
    temperature_k: float,
    c: pt.ParticleConstants,
) -> float:
    """Total internal wet-side loading calculated from one canonical ledger."""

    return (
        particle_wet_critical_phase_ledger(temperature_k, c)
        .critical_hexane_loading_kg_per_kg_dry_meal
    )


def particle_front_mobile_jump_ledger(
    temperature_k: float,
    pore_hexane_density_kg_m3: float,
    c: pt.ParticleConstants,
) -> FrontMobileJumpLedger:
    """Exact wet-minus-dry n-hexane jump at the moving interface."""

    if pore_hexane_density_kg_m3 < 0.0 or not math.isfinite(
        pore_hexane_density_kg_m3
    ):
        raise ValueError("dry-side pore-hexane density must be finite and non-negative")
    wet = particle_wet_critical_phase_ledger(temperature_k, c)
    dry_density = c.alpha_ps * c.rho_ps
    retained_loading, _ = thermo.x2_so_and_slope(
        1.0,
        temperature_k,
        c.X3,
        c.gab,
        c.oil,
    )
    dry_side = DrySideFrontInventory.from_pore_gas_hexane_density(
        retained_hexane_loading_kg_per_kg_dry_meal=retained_loading,
        pore_gas_hexane_density_kg_m3=pore_hexane_density_kg_m3,
        pore_gas_specific_volume_m3_per_kg_dry_meal=(
            c.alpha_pg / dry_density
        ),
    )
    return solve_front_mobile_jump_ledger(
        wet_critical_ledger=wet,
        dry_side_inventory=dry_side,
    )


def particle_front_mobile_jump_capacity_kg_m3(
    temperature_k: float,
    pore_hexane_density_kg_m3: float,
    c: pt.ParticleConstants,
) -> float:
    """Physical sharp-front reserve bound from the exact phase ledger."""

    return particle_front_mobile_jump_ledger(
        temperature_k,
        pore_hexane_density_kg_m3,
        c,
    ).mobile_jump_density_kg_m3_particle_envelope


def particle_mobile_reserve_numerical_upper_bound_kg_m3(
    c: pt.ParticleConstants,
) -> float:
    """Loose state-independent solver bound, never a physical capacity.

    The full wet critical concentration bounds the smaller wet-minus-dry
    jump for every admissible dry-side inventory.  The configured reference
    liquid density is used only to size a nonlinear-search box.
    """

    return (
        c.critical_volume_multiplier
        * c.alpha_pg
        * c.rho_hexane_liquid_ref
    )


def particle_activity_from_inventory_density(
    inventory_density_kg_m3: float,
    temperature_k: float,
    c: pt.ParticleConstants,
) -> float:
    """Invert the monotone A.22 storage relation on ``0 <= activity <= 1``."""

    if not math.isfinite(inventory_density_kg_m3):
        raise ValueError("particle hexane inventory must be finite")
    if temperature_k <= 0.0 or not math.isfinite(temperature_k):
        raise ValueError("particle temperature must be positive and finite")
    zero_inventory = particle_hexane_inventory_density_kg_m3(
        0.0,
        temperature_k,
        c,
    )
    saturated_inventory = particle_hexane_saturation_density_kg_m3(
        temperature_k,
        c,
    )
    scale = max(abs(saturated_inventory), 1.0)
    tolerance = 1.0e-12 * scale
    if inventory_density_kg_m3 < zero_inventory - tolerance:
        raise ValueError("particle hexane inventory is below the dry-state value")
    if inventory_density_kg_m3 > saturated_inventory + tolerance:
        raise ValueError("equilibrium inventory exceeds unit activity")
    if inventory_density_kg_m3 <= zero_inventory + tolerance:
        return 0.0
    if inventory_density_kg_m3 >= saturated_inventory - tolerance:
        return 1.0

    def residual(activity: float) -> float:
        return (
            particle_hexane_inventory_density_kg_m3(
                activity,
                temperature_k,
                c,
            )
            - inventory_density_kg_m3
        )

    return float(brentq(residual, 0.0, 1.0, xtol=1.0e-13, rtol=1.0e-13))


def liquid_hexane_specific_enthalpy_j_kg(
    temperature_k: float,
    particle: pt.ParticleConstants,
    ref: thermo.VaporEnthalpyRef,
) -> float:
    """Liquid-hexane sensible enthalpy on the common boiling-point datum."""

    return particle.cp_hexane_liquid * (
        temperature_k - ref.T_boil_hexane
    )


def particle_total_enthalpy_density_from_phase_coordinates_j_m3(
    activity: float,
    temperature_k: float,
    total_water_kg_kg_dry: float,
    pore_hexane_inventory_kg_m3: float,
    saturation_reserve_kg_m3: float,
    particle: pt.ParticleConstants,
    ref: thermo.VaporEnthalpyRef,
    *,
    dry_reference_temperature_k: float = 273.15,
    allow_signed_reserve: bool = False,
) -> float:
    """Particle enthalpy for an explicit pore-hexane inventory coordinate.

    ``saturation_reserve_kg_m3`` is the wet-minus-dry component jump, not a
    second liquid mass added on top of an unchanged dry pore gas.  A partially
    wet shell therefore interpolates the *phase ledger*: the same wet fraction
    that adds the mobile mass jump also converts the corresponding dry-side
    pore-gas inventory from vapor to mobile-liquid enthalpy.  At a full
    reserve this recovers the approved wet state (retained plus mobile liquid,
    with no wet pore gas) exactly; at zero reserve it recovers the dry state.
    """

    if pore_hexane_inventory_kg_m3 < 0.0:
        raise ValueError("pore hexane inventory must be non-negative")
    bound_loading, _ = thermo.x2_so_and_slope(
        activity,
        temperature_k,
        particle.X3,
        particle.gab,
        particle.oil,
    )
    adsorbed_loading = thermo.gab_hexane_content(
        activity,
        temperature_k,
        particle.gab,
    )
    oil_loading = thermo.oil_hexane_content(activity, particle.oil)
    retained_hexane = (
        particle.alpha_ps * particle.rho_ps * bound_loading
    )
    _, h_hexane = vapor_specific_enthalpies_j_kg(temperature_k, ref)
    binding_potential = particle.alpha_ps * particle.rho_ps * (
        pt.sorption_enthalpy_potential_j_kg_dry(
            adsorbed_loading,
            particle,
        )
        + particle.X3 * particle.dH_vap_hexane * oil_loading
    )
    dry_sensible = (
        particle.alpha_ps
        * particle.rho_ps
        * particle.cp_ps
        * (temperature_k - dry_reference_temperature_k)
    )
    oil_sensible = (
        particle.alpha_ps
        * particle.rho_ps
        * particle.X3
        * particle.cp_oil
        * (temperature_k - dry_reference_temperature_k)
    )
    water_sensible = (
        particle.alpha_ps
        * particle.rho_ps
        * total_water_kg_kg_dry
        * particle.cp_water_liquid
        * (temperature_k - ref.T_boil_water)
    )
    if saturation_reserve_kg_m3 < 0.0 and not allow_signed_reserve:
        raise ValueError("hexane saturation reserve must be non-negative")
    front_phase_energy_correction = 0.0
    if saturation_reserve_kg_m3 != 0.0:
        pore_density = pore_hexane_inventory_kg_m3 / particle.alpha_pg
        reserve_capacity = particle_front_mobile_jump_capacity_kg_m3(
            temperature_k,
            pore_density,
            particle,
        )
        if reserve_capacity <= 0.0:
            raise ValueError(
                "nonzero hexane saturation reserve requires a positive "
                "critical front jump"
            )
        if (
            not allow_signed_reserve
            and saturation_reserve_kg_m3 > reserve_capacity + 1.0e-10
        ):
            raise ValueError(
                "hexane saturation reserve exceeds the critical front jump"
            )
        wet_fraction = saturation_reserve_kg_m3 / reserve_capacity
        h_hexane_liquid = liquid_hexane_specific_enthalpy_j_kg(
            temperature_k,
            particle,
            ref,
        )
        front_phase_energy_correction = (
            saturation_reserve_kg_m3 * h_hexane_liquid
            + wet_fraction
            * pore_hexane_inventory_kg_m3
            * (h_hexane_liquid - h_hexane)
        )
    return (
        dry_sensible
        + oil_sensible
        + water_sensible
        + (pore_hexane_inventory_kg_m3 + retained_hexane) * h_hexane
        - binding_potential
        + front_phase_energy_correction
    )


def particle_total_enthalpy_density_j_m3(
    activity: float,
    temperature_k: float,
    total_water_kg_kg_dry: float,
    particle: pt.ParticleConstants,
    ref: thermo.VaporEnthalpyRef,
    *,
    dry_reference_temperature_k: float = 273.15,
    saturation_reserve_kg_m3: float = 0.0,
    allow_signed_reserve: bool = False,
) -> float:
    """Conservative particle enthalpy underlying steady and dynamic DCZ use.

    Coletto A.22 transports pore, adsorbed, and oil-dissolved hexane; A.23
    and A.30 account for sensible energy and differential sorption heat; and
    A.34 carries transferred hexane enthalpy in the vapor.  Combining those
    terms gives a state function: dry/water sensible enthalpy plus the vapor
    enthalpy of all retained hexane, minus the integrated binding potential.
    A diffusive hexane flux therefore carries the same vapor enthalpy on both
    sides of the particle/vapor interface, while the potential supplies the
    finite heat of sorption without a point-derivative approximation.
    """

    return particle_total_enthalpy_density_from_phase_coordinates_j_m3(
        activity,
        temperature_k,
        total_water_kg_kg_dry,
        particle.alpha_pg * particle.rho_pg * activity,
        saturation_reserve_kg_m3,
        particle,
        ref,
        dry_reference_temperature_k=dry_reference_temperature_k,
        allow_signed_reserve=allow_signed_reserve,
    )


def particle_state_from_inventory_enthalpy(
    inventory_density_kg_m3: float,
    enthalpy_density_j_m3: float,
    total_water_kg_kg_dry: float,
    particle: pt.ParticleConstants,
    ref: thermo.VaporEnthalpyRef,
    *,
    dry_reference_temperature_k: float = 273.15,
    minimum_temperature_k: float = 250.0,
    maximum_temperature_k: float = 480.0,
) -> ParticleConstitutiveState:
    """Recover activity, temperature, and saturation reserve conservatively.

    Below saturation the A.22 storage relation is inverted for activity.  At
    saturation activity stays exactly one and excess conserved inventory is
    carried as a pore-volume-bounded liquid-enthalpy reserve.  Temperature is
    the root of the corresponding total-enthalpy state function, so neither
    component mass nor energy is projected away at the saturation switch.
    """

    if inventory_density_kg_m3 < 0.0 or not math.isfinite(
        inventory_density_kg_m3
    ):
        raise ValueError("particle hexane inventory must be non-negative and finite")
    if not math.isfinite(enthalpy_density_j_m3):
        raise ValueError("particle enthalpy must be finite")
    if total_water_kg_kg_dry < 0.0 or not math.isfinite(
        total_water_kg_kg_dry
    ):
        raise ValueError("particle water loading must be non-negative and finite")
    if (
        minimum_temperature_k <= 0.0
        or maximum_temperature_k <= minimum_temperature_k
    ):
        raise ValueError("invalid particle-temperature inversion interval")

    def state_at_temperature(
        temperature_k: float,
    ) -> tuple[float, float, float]:
        saturated_inventory = particle_hexane_saturation_density_kg_m3(
            temperature_k,
            particle,
        )
        if inventory_density_kg_m3 <= saturated_inventory:
            activity = particle_activity_from_inventory_density(
                inventory_density_kg_m3,
                temperature_k,
                particle,
            )
            reserve = 0.0
        else:
            activity = 1.0
            reserve = inventory_density_kg_m3 - saturated_inventory
            reserve_capacity = particle_front_mobile_jump_capacity_kg_m3(
                temperature_k,
                particle.rho_pg,
                particle,
            )
            if reserve > reserve_capacity + 1.0e-10:
                raise ValueError(
                    "hexane saturation reserve exceeds the critical front jump"
                )
        calculated_enthalpy = particle_total_enthalpy_density_j_m3(
            activity,
            temperature_k,
            total_water_kg_kg_dry,
            particle,
            ref,
            dry_reference_temperature_k=dry_reference_temperature_k,
            saturation_reserve_kg_m3=reserve,
        )
        return activity, reserve, calculated_enthalpy

    def residual(temperature_k: float) -> float:
        return state_at_temperature(temperature_k)[2] - enthalpy_density_j_m3

    residual_minimum = residual(minimum_temperature_k)
    residual_maximum = residual(maximum_temperature_k)
    enthalpy_scale = max(abs(enthalpy_density_j_m3), 1.0)
    root_tolerance = 1.0e-11 * enthalpy_scale
    if abs(residual_minimum) <= root_tolerance:
        temperature = minimum_temperature_k
    elif abs(residual_maximum) <= root_tolerance:
        temperature = maximum_temperature_k
    elif residual_minimum * residual_maximum > 0.0:
        raise ValueError(
            "particle inventory/enthalpy state has no temperature root "
            f"within [{minimum_temperature_k:g}, {maximum_temperature_k:g}] K"
        )
    else:
        temperature = float(
            brentq(
                residual,
                minimum_temperature_k,
                maximum_temperature_k,
                xtol=1.0e-10,
                rtol=1.0e-12,
            )
        )
    activity, reserve, _ = state_at_temperature(temperature)
    return ParticleConstitutiveState(
        activity=float(activity),
        temperature_k=float(temperature),
        saturation_reserve_kg_m3=float(reserve),
    )


def evaporation_required_for_wet_temperature_constraint(
    temperature_after_previous_evaporation_k: float,
    saturation_temperature_k: float,
    previous_evaporation_kg_s: float,
    cell_thermal_coefficient_w_k: float,
    phase_change_enthalpy_j_kg: float,
) -> float:
    """Invert the affine cell-energy equation for total evaporation."""
    if cell_thermal_coefficient_w_k <= 0.0:
        raise ValueError("cell_thermal_coefficient_w_k must be positive")
    if phase_change_enthalpy_j_kg <= 0.0:
        raise ValueError("phase_change_enthalpy_j_kg must be positive")
    previous_evaporation_kg_s = max(previous_evaporation_kg_s, 0.0)
    required_energy_w = (
        previous_evaporation_kg_s * phase_change_enthalpy_j_kg
        + cell_thermal_coefficient_w_k
        * (temperature_after_previous_evaporation_k - saturation_temperature_k)
    )
    return max(required_energy_w / phase_change_enthalpy_j_kg, 0.0)


@dataclass(frozen=True)
class VaporCellEnthalpyResult:
    water_out_kg_s: float
    hexane_out_kg_s: float
    temperature_k: float
    condensed_water_kg_s: float
    energy_residual_w: float


def vapor_cell_enthalpy_residual_at_state(
    *,
    water_in_kg_s: float,
    hexane_in_kg_s: float,
    temperature_in_k: float,
    water_out_kg_s: float,
    hexane_out_kg_s: float,
    temperature_out_k: float,
    bound_water_to_solid_kg_s: float,
    free_water_evaporation_kg_s: float,
    condensed_water_kg_s: float,
    external_duty_w: float,
    interphase_heat_to_vapor_w: float,
    particle_surface_temperature_k: float,
    bound_water_surface_temperature_k: float,
    free_water_temperature_k: float,
    c: DCZConstants,
) -> float:
    """First-law residual of an arbitrary accepted DCZ vapor-cell state."""
    ref = c.vapor_enthalpy_ref
    enthalpy_in_w = water_in_kg_s * thermo.vapor_enthalpy_water_basis(
        hexane_in_kg_s / max(water_in_kg_s, 1.0e-12),
        temperature_in_k,
        ref,
    )
    enthalpy_out_w = water_out_kg_s * thermo.vapor_enthalpy_water_basis(
        hexane_out_kg_s / max(water_out_kg_s, 1.0e-12),
        temperature_out_k,
        ref,
    )
    transferred_hexane_enthalpy_w = (
        (hexane_out_kg_s - hexane_in_kg_s)
        * (
            ref.dH_vap_hexane
            + ref.cp_hexane_vapor
            * (particle_surface_temperature_k - ref.T_boil_hexane)
        )
    )
    bound_liquid_enthalpy_w = (
        bound_water_to_solid_kg_s
        * thermo.liquid_water_enthalpy(
            bound_water_surface_temperature_k,
            ref,
            c.particle.cp_water_liquid,
        )
    )
    free_liquid_enthalpy_w = (
        free_water_evaporation_kg_s
        * thermo.liquid_water_enthalpy(
            free_water_temperature_k,
            ref,
            c.particle.cp_water_liquid,
        )
    )
    condensed_liquid_enthalpy_w = (
        condensed_water_kg_s
        * thermo.liquid_water_enthalpy(
            temperature_out_k,
            ref,
            c.particle.cp_water_liquid,
        )
    )
    return enthalpy_out_w - (
        enthalpy_in_w
        + external_duty_w
        + interphase_heat_to_vapor_w
        + transferred_hexane_enthalpy_w
        - bound_liquid_enthalpy_w
        + free_liquid_enthalpy_w
        - condensed_liquid_enthalpy_w
    )


def solve_vapor_cell_temperature_at_state(
    *,
    water_in_kg_s: float,
    hexane_in_kg_s: float,
    temperature_in_k: float,
    water_out_kg_s: float,
    hexane_out_kg_s: float,
    bound_water_to_solid_kg_s: float,
    free_water_evaporation_kg_s: float,
    condensed_water_kg_s: float,
    external_duty_w: float,
    interphase_ua_w_k: float,
    particle_surface_temperature_k: float,
    bound_water_surface_temperature_k: float,
    free_water_temperature_k: float,
    c: DCZConstants,
    interphase_particle_temperature_k: float | None = None,
) -> float:
    """Exact fixed-mass temperature projection for one accepted cell."""
    interphase_temperature = (
        particle_surface_temperature_k
        if interphase_particle_temperature_k is None
        else interphase_particle_temperature_k
    )

    def residual(temperature_out_k: float) -> float:
        return vapor_cell_enthalpy_residual_at_state(
            water_in_kg_s=water_in_kg_s,
            hexane_in_kg_s=hexane_in_kg_s,
            temperature_in_k=temperature_in_k,
            water_out_kg_s=water_out_kg_s,
            hexane_out_kg_s=hexane_out_kg_s,
            temperature_out_k=temperature_out_k,
            bound_water_to_solid_kg_s=bound_water_to_solid_kg_s,
            free_water_evaporation_kg_s=free_water_evaporation_kg_s,
            condensed_water_kg_s=condensed_water_kg_s,
            external_duty_w=external_duty_w,
            interphase_heat_to_vapor_w=interphase_ua_w_k
            * (interphase_temperature - temperature_out_k),
            particle_surface_temperature_k=particle_surface_temperature_k,
            bound_water_surface_temperature_k=bound_water_surface_temperature_k,
            free_water_temperature_k=free_water_temperature_k,
            c=c,
        )

    residual_zero = residual(0.0)
    slope = residual(1.0) - residual_zero
    if not math.isfinite(slope) or abs(slope) < 1.0e-12:
        raise ValueError("DCZ fixed-state enthalpy equation is singular")
    temperature = -residual_zero / slope
    if not math.isfinite(temperature):
        raise ValueError("DCZ fixed-state enthalpy projection is non-finite")
    return temperature


def solve_vapor_cell_enthalpy(
    *,
    water_in_kg_s: float,
    hexane_in_kg_s: float,
    temperature_in_k: float,
    hexane_transfer_kg_s: float,
    bound_water_to_solid_kg_s: float,
    free_water_evaporation_kg_s: float,
    external_duty_w: float,
    interphase_ua_w_k: float,
    particle_surface_temperature_k: float,
    bound_water_surface_temperature_k: float,
    free_water_temperature_k: float,
    c: DCZConstants,
    interphase_particle_temperature_k: float | None = None,
) -> VaporCellEnthalpyResult:
    """Close one DCZ vapor cell on component flow and total enthalpy."""
    if water_in_kg_s <= 0.0:
        raise ValueError("DCZ vapor-cell water inlet flow must be positive")
    if interphase_ua_w_k < 0.0:
        raise ValueError("DCZ vapor-cell interphase UA must be non-negative")

    interphase_temperature = (
        particle_surface_temperature_k
        if interphase_particle_temperature_k is None
        else interphase_particle_temperature_k
    )
    water_in_kg_s = max(water_in_kg_s, 1.0e-12)
    hexane_out_kg_s = max(hexane_in_kg_s + hexane_transfer_kg_s, 0.0)
    hexane_transfer_kg_s = hexane_out_kg_s - hexane_in_kg_s
    bound_water_to_solid_kg_s = min(
        bound_water_to_solid_kg_s,
        water_in_kg_s + free_water_evaporation_kg_s - 1.0e-12,
    )
    water_before_condensation_kg_s = max(
        water_in_kg_s
        - bound_water_to_solid_kg_s
        + free_water_evaporation_kg_s,
        1.0e-12,
    )

    def residual_at(
        temperature_k: float,
        condensed_water_kg_s: float,
    ) -> float:
        water_out_kg_s = water_before_condensation_kg_s - condensed_water_kg_s
        return vapor_cell_enthalpy_residual_at_state(
            water_in_kg_s=water_in_kg_s,
            hexane_in_kg_s=hexane_in_kg_s,
            temperature_in_k=temperature_in_k,
            water_out_kg_s=water_out_kg_s,
            hexane_out_kg_s=hexane_out_kg_s,
            temperature_out_k=temperature_k,
            bound_water_to_solid_kg_s=bound_water_to_solid_kg_s,
            free_water_evaporation_kg_s=free_water_evaporation_kg_s,
            condensed_water_kg_s=condensed_water_kg_s,
            external_duty_w=external_duty_w,
            interphase_heat_to_vapor_w=interphase_ua_w_k
            * (interphase_temperature - temperature_k),
            particle_surface_temperature_k=particle_surface_temperature_k,
            bound_water_surface_temperature_k=bound_water_surface_temperature_k,
            free_water_temperature_k=free_water_temperature_k,
            c=c,
        )

    water_out = water_before_condensation_kg_s
    vapor_heat_capacity_w_k = (
        water_out * c.vapor_enthalpy_ref.cp_water_vapor
        + hexane_out_kg_s * c.vapor_enthalpy_ref.cp_hexane_vapor
    )
    residual_at_zero_k = residual_at(0.0, 0.0)
    temperature_candidate = -residual_at_zero_k / (
        vapor_heat_capacity_w_k + interphase_ua_w_k
    )
    candidate_ratio = hexane_out_kg_s / water_out
    dew_temperature = thermo.dew_point_temperature(
        candidate_ratio,
        c.antoine_water,
        P=c.pressure_pa,
    )
    if temperature_candidate >= dew_temperature:
        residual = residual_at(temperature_candidate, 0.0)
        return VaporCellEnthalpyResult(
            water_out_kg_s=water_out,
            hexane_out_kg_s=hexane_out_kg_s,
            temperature_k=temperature_candidate,
            condensed_water_kg_s=0.0,
            energy_residual_w=residual,
        )

    upper = water_before_condensation_kg_s * 0.999

    def saturated_residual(condensed_water_kg_s: float) -> float:
        remaining_water = water_before_condensation_kg_s - condensed_water_kg_s
        ratio = hexane_out_kg_s / remaining_water
        dew_point = thermo.dew_point_temperature(
            ratio,
            c.antoine_water,
            P=c.pressure_pa,
        )
        return residual_at(dew_point, condensed_water_kg_s)

    residual_lo = saturated_residual(0.0)
    residual_hi = saturated_residual(upper)
    if residual_lo * residual_hi > 0.0:
        raise ValueError(
            "DCZ saturated vapor-cell enthalpy balance has no bounded "
            "condensation root "
            f"(residual_zero_condensation={residual_lo:.6g} W, "
            f"residual_near_total_condensation={residual_hi:.6g} W, "
            f"water_in={water_in_kg_s:.6g} kg/s, "
            f"hexane_in={hexane_in_kg_s:.6g} kg/s, "
            f"temperature_in={temperature_in_k:.6g} K, "
            f"external_duty={external_duty_w:.6g} W, "
            f"interphase_UA={interphase_ua_w_k:.6g} W/K, "
            f"particle_surface_temperature={particle_surface_temperature_k:.6g} K)"
        )
    condensed = brentq(saturated_residual, 0.0, upper)
    water_out = water_before_condensation_kg_s - condensed
    ratio = hexane_out_kg_s / water_out
    temperature_out = thermo.dew_point_temperature(
        ratio,
        c.antoine_water,
        P=c.pressure_pa,
    )
    return VaporCellEnthalpyResult(
        water_out_kg_s=water_out,
        hexane_out_kg_s=hexane_out_kg_s,
        temperature_k=temperature_out,
        condensed_water_kg_s=condensed,
        energy_residual_w=residual_at(temperature_out, condensed),
    )
