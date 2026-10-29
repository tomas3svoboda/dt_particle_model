"""Tray-level PHY-022 film activation (T2b tranche 2).

Composition layer over the whitelisted particle-side film machinery —
imports only, no new constitutive physics
(docs/T2B_FILM_DESIGN_MAP_2026-08-06.md).  Exposes the coupled binary
water/hexane film primitives a tray layer needs:

- the sanctioned Coletto B.7–B.10 coefficient pair behind the
  fail-closed `ReducedFilmOperatingPoint` envelope;
- the rank-one Stefan partition of a LAYER-SUPPLIED total mass flux
  (the film never invents j_t — frozen rank contract);
- the component-capacity Ackermann heat flux with exact PHY-053
  partial-enthalpy secant capacities, operating point rebuilt per call
  (stale-temperature hazard);
- the two frozen limit oracles as testable functions: the
  independent-linear low-flux limit and the Faner pure-hexane
  heat-controlled capacity.

Sign conventions follow the whitelisted modules: component mass fluxes
positive solid→gas (kg m⁻² s⁻¹); heat positive gas→solid.  Area is NOT
applied here (T2c owns the area law; fluxes are per interfacial area).
physically_qualifying is False throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .particle import bed_film as bf
from .particle import external_film_reduced as efr
from .props import binary_gas
from .props import hexane as _hexane
from .props import water as _water


class TrayFilmError(RuntimeError):
    """A tray-film contract was violated; the evaluation is unusable."""


def hexane_mass_fraction_from_mole_fraction(y_hexane: float) -> float:
    """Mole → mass fraction with the frozen molar masses."""
    if not 0.0 <= y_hexane <= 1.0:
        raise TrayFilmError(f"y_hexane {y_hexane!r} outside [0, 1]")
    numerator = y_hexane * _hexane.M
    return numerator / (numerator + (1.0 - y_hexane) * _water.M)


@dataclass(frozen=True)
class TrayFilmHeat:
    """Ackermann-corrected film heat flux at one (T_s, T_b) evaluation."""

    heat_into_solid_w_m2: float
    ackermann_factor: float
    ackermann_beta: float
    water_capacity_j_kg_k: float
    hexane_capacity_j_kg_k: float

    @property
    def physically_qualifying(self) -> bool:
        return False


class TrayLayerFilm:
    """PHY-022 film primitives for one tray layer at the certified point.

    One instance = one hydrodynamic state (gas properties + superficial
    velocity + particle/bed geometry).  All evaluations go through the
    whitelisted fail-closed constructors; envelope violations raise the
    underlying `ReducedFilmValidityError` untouched.
    """

    physically_qualifying = False

    def __init__(
        self,
        *,
        particle_radius_m: float,
        bed_void_fraction: float,
        gas_density_kg_m3: float,
        gas_superficial_velocity_m_s: float,
        gas_dynamic_viscosity_pa_s: float,
        gas_heat_capacity_j_kg_k: float,
        gas_thermal_conductivity_w_m_k: float,
        pressure_pa: float,
        forcing_timescale_s: float,
        stress_multiplier: float = 1.0,
        state_provenance: str = "T2b tray-layer gas state",
    ) -> None:
        self._state = bf.PackedBedFilmState(
            particle_radius_m=particle_radius_m,
            bed_void_fraction=bed_void_fraction,
            gas_density_kg_m3=gas_density_kg_m3,
            gas_superficial_velocity_m_s=gas_superficial_velocity_m_s,
            gas_dynamic_viscosity_pa_s=gas_dynamic_viscosity_pa_s,
            gas_heat_capacity_j_kg_k=gas_heat_capacity_j_kg_k,
            gas_thermal_conductivity_w_m_k=gas_thermal_conductivity_w_m_k,
            bulk_binary_diffusivity_m2_s=efr.COLETTO_NOMINAL_BINARY_DIFFUSIVITY_M2_S,
            state_provenance=state_provenance,
            reynolds_convention_provenance=(
                "Faner (2008) Eq. 4.23 solid-fraction superficial Reynolds "
                "convention, D20 2026-09-23 (was: Coletto voidage-corrected)"
            ),
        )
        self._pair = bf.coletto_b7_b10_pair(self._state)
        self._pressure_pa = pressure_pa
        self._forcing_timescale_s = forcing_timescale_s
        self._stress_multiplier = stress_multiplier

    # -- coefficient adaptation (rebuilt per temperature pair) -----------
    def coefficients(self, *, bulk_temperature_k: float, interface_temperature_k: float):
        operating_point = efr.ReducedFilmOperatingPoint(
            bulk_temperature_k=bulk_temperature_k,
            interface_temperature_k=interface_temperature_k,
            pressure_pa=self._pressure_pa,
            particle_radius_m=self._state.particle_radius_m,
            liquid_saturation=0.0,
            forcing_timescale_s=self._forcing_timescale_s,
            drainage_regime=efr.DrainageRegime.PRIMARY_DRAINAGE,
        )
        return efr.adapt_correlated_bed_film(
            self._pair, operating_point, stress_multiplier=self._stress_multiplier
        )

    # -- rank-one Stefan partition ---------------------------------------
    def partition_total_mass_flux(
        self,
        *,
        total_mass_flux_kg_m2_s: float,
        interface_y_hexane: float,
        bulk_y_hexane: float,
        bulk_temperature_k: float,
        interface_temperature_k: float,
    ) -> efr.ComponentMassFluxes:
        """Partition a LAYER-SUPPLIED total mass flux into components.

        The total flux is the layer balance's degree of freedom (frozen
        rank contract) — this method never invents it.
        """
        film = self.coefficients(
            bulk_temperature_k=bulk_temperature_k,
            interface_temperature_k=interface_temperature_k,
        )
        return efr.component_fluxes_for_supplied_total_mass_flux(
            film,
            total_mass_flux_kg_m2_s=total_mass_flux_kg_m2_s,
            interface_hexane_mass_fraction=hexane_mass_fraction_from_mole_fraction(
                interface_y_hexane
            ),
            bulk_hexane_mass_fraction=hexane_mass_fraction_from_mole_fraction(
                bulk_y_hexane
            ),
        )

    # -- Ackermann heat ---------------------------------------------------
    def ackermann_heat(
        self,
        *,
        fluxes: efr.ComponentMassFluxes,
        bulk_temperature_k: float,
        interface_temperature_k: float,
        y_hexane: float,
        k_wh: float = binary_gas.K_WH_CENTRAL,
    ) -> TrayFilmHeat:
        """Component-capacity Ackermann heat flux (gas → solid positive).

        Operating point and coefficients are rebuilt at THIS temperature
        pair (the whitelisted primitive reads temperatures from the
        operating point, never from arguments).
        """
        film = self.coefficients(
            bulk_temperature_k=bulk_temperature_k,
            interface_temperature_k=interface_temperature_k,
        )
        capacities = efr.phy053_partial_enthalpy_secant_heat_capacities(
            surface_temperature_k=interface_temperature_k,
            bulk_temperature_k=bulk_temperature_k,
            pressure_pa=self._pressure_pa,
            y_hexane=y_hexane,
            k_wh=k_wh,
        )
        result = efr.ackermann_heat_transfer(
            film,
            fluxes=fluxes,
            water_film_heat_capacity_j_kg_k=capacities.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=capacities.hexane_j_kg_k,
        )
        return TrayFilmHeat(
            heat_into_solid_w_m2=result.heat_into_particle_w_m2,
            ackermann_factor=result.factor,
            ackermann_beta=result.beta,
            water_capacity_j_kg_k=capacities.water_j_kg_k,
            hexane_capacity_j_kg_k=capacities.hexane_j_kg_k,
        )


def faner_heat_controlled_capacity_kg_m2_s(
    *,
    heat_transfer_w_m2_k: float,
    gas_heat_capacity_j_kg_k: float,
    gas_temperature_k: float,
    boiling_temperature_k: float,
    latent_heat_j_kg: float,
    dry_shell_resistance: float = 0.0,
) -> float:
    """Faner's analytic pure-solvent heat-controlled evaporation capacity.

    N = (h/c_p)·ln(1 + c_p·ΔT/(ΔH·(1+R_shell))) per unit transfer area —
    the exact elimination of the implicit Ackermann factor (Faner eqs 3,
    5, 6; thin-layer oracle uses R_shell = 0).  Returns 0 when the gas is
    not superheated.  This is the mandated LIMIT oracle, never the
    general binary law.
    """
    if heat_transfer_w_m2_k <= 0.0 or gas_heat_capacity_j_kg_k <= 0.0:
        raise TrayFilmError("h and c_p must be strictly positive")
    if latent_heat_j_kg <= 0.0:
        raise TrayFilmError("latent heat must be strictly positive")
    superheat = gas_temperature_k - boiling_temperature_k
    if superheat <= 0.0:
        return 0.0
    return (heat_transfer_w_m2_k / gas_heat_capacity_j_kg_k) * math.log1p(
        gas_heat_capacity_j_kg_k
        * superheat
        / (latent_heat_j_kg * (1.0 + dry_shell_resistance))
    )
