"""Owner-authorized bounded Phase-1 particle engineering configuration.

This module is the single construction path for ``P1EF`` qualification.  It
does not promote the normative F2/F3 physics gates: every object remains
``physically_qualifying=False`` and ``plant_predictive=False``.  The
construction preserves the frozen conservative water, n-hexane, and
common-datum energy equations while selecting only the narrowly authorized
coefficient manifold, fixed geometry, structural conductivity, and compact
external-film reduction recorded by ``GT-PS-2-P1E-01``.

The coefficient pair is selected once at the declared smooth design state and
then held fixed for a trajectory.  The heat/mass film pair is generated from
one Coletto B.7--B.10 state and receives one common multiplier.  No fitted
state law, coefficient floor, or independent heat/mass corner is introduced.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

from dtdc_simulator.core2.particle import bed_film as bf
from dtdc_simulator.core2.particle import coefficient_manifold as cm
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import external_film_reduced as efr
from dtdc_simulator.core2.particle import nonisothermal_potential as nip
from dtdc_simulator.core2.particle import transport_coefficients as tc
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.props import binary_gas as bg


P1EF_DECISION_ID = "GT-PS-2-P1E-01"
PRODUCTION_RADIUS_M = 0.885e-3
DESIGN_TEMPERATURE_K = 385.65
DESIGN_PRESSURE_PA = 130_000.0
DESIGN_Y_HEXANE = 0.80
NOMINAL_PROCESS_PRESSURE_PA = 101_325.0

NOMINAL_BINARY_DIFFUSIVITY_M2_S = 4.0e-10
NOMINAL_APPARENT_WATER_DIFFUSIVITY_M2_S = math.sqrt(
    cm.APPARENT_WATER_DIFFUSIVITY_MIN_M2_S
    * cm.APPARENT_WATER_DIFFUSIVITY_MAX_M2_S
)

NOMINAL_STRUCTURAL_CONDUCTIVITY_W_M_K = 0.24
STRUCTURAL_CONDUCTIVITY_SENSITIVITY_W_M_K = 0.29
ALLOWED_STRUCTURAL_CONDUCTIVITIES_W_M_K = frozenset(
    (
        NOMINAL_STRUCTURAL_CONDUCTIVITY_W_M_K,
        STRUCTURAL_CONDUCTIVITY_SENSITIVITY_W_M_K,
    )
)

# The declared lower corner (GT-PS-2-OWNER-RULING-FILM-MULTIPLIER-RANGE-2026-09-30-01,
# external_film_reduced.REDUCED_FILM_VALIDITY_LOWER_CORNER_DECLARATION); 0.5 is
# kept so that it still refuses at the validity check below.
COMPACT_FILM_STRESS_MULTIPLIERS = (
    0.5,
    efr.REDUCED_FILM_VALIDITY_LOWER_CORNER_MULTIPLIER,
    1.0,
    2.0,
)
FAST_FILM_MASS_FRACTION_DEPARTURE_LIMIT = 0.01
UNITY_ACKERMANN_FACTOR_DEPARTURE_LIMIT = 0.005
CONDITIONAL_RETAINED_SURFACE_COEFFICIENT_FLOOR_M_S = 0.006
MAXIMUM_RETAINED_SURFACE_RESISTANCE_FRACTION = 0.05
LEGACY_ENGINEERING_REGRESSION_WATER_HEAT_CAPACITY_J_KG_K = 1900.0
LEGACY_ENGINEERING_REGRESSION_HEXANE_HEAT_CAPACITY_J_KG_K = 1650.0


@dataclass(frozen=True)
class EngineeringCoefficientCase:
    """One manifold-generated trajectory-wise constant coefficient pair."""

    identifier: str
    binary_diffusivity_m2_s: float
    apparent_water_diffusivity_m2_s: float
    partition: cm.CoefficientPartitionSelection
    physically_qualifying: bool = field(default=False, init=False)
    plant_predictive: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.identifier, str) or not self.identifier.strip():
            raise ValueError("engineering coefficient case needs an identifier")
        if not isinstance(self.partition, cm.CoefficientPartitionSelection):
            raise TypeError("engineering coefficient case needs a manifold partition")
        request = self.partition.request
        if request.binary_diffusivity.value_m2_s != self.binary_diffusivity_m2_s:
            raise ValueError("binary coefficient differs from its manifold request")
        if (
            request.apparent_water_diffusivity.value_m2_s
            != self.apparent_water_diffusivity_m2_s
        ):
            raise ValueError("apparent-water target differs from its manifold request")


@dataclass(frozen=True)
class EngineeringParticleConfiguration:
    """Matched wet/dry/front configuration for one P1EF campaign case."""

    coefficient_case: EngineeringCoefficientCase
    structural_conductivity_w_m_k: float
    process_pressure_pa: float
    dry: ct.FullyDryTransportConfig
    wet: ww.WetWaterModel
    cut: cut.CutTransportConfig
    radius_m: float = PRODUCTION_RADIUS_M
    decision_id: str = P1EF_DECISION_ID
    conditional_retained_surface_coefficient_floor_m_s: float = (
        CONDITIONAL_RETAINED_SURFACE_COEFFICIENT_FLOOR_M_S
    )
    maximum_retained_surface_resistance_fraction: float = (
        MAXIMUM_RETAINED_SURFACE_RESISTANCE_FRACTION
    )
    retained_surface_equilibrium_limit_is_conditional: bool = field(
        default=True,
        init=False,
    )
    retained_surface_coefficient_is_target_identified: bool = field(
        default=False,
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)
    plant_predictive: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.radius_m != PRODUCTION_RADIUS_M:
            raise ValueError("P1EF production geometry is fixed at 0.885 mm")
        if self.conditional_retained_surface_coefficient_floor_m_s != (
            CONDITIONAL_RETAINED_SURFACE_COEFFICIENT_FLOOR_M_S
        ):
            raise ValueError(
                "P1EF equilibrium retained-surface limit requires the declared "
                "conditional K_w,s floor"
            )
        if self.maximum_retained_surface_resistance_fraction != (
            MAXIMUM_RETAINED_SURFACE_RESISTANCE_FRACTION
        ):
            raise ValueError(
                "P1EF retained-surface reduction threshold must remain 0.05"
            )
        if self.structural_conductivity_w_m_k not in (
            ALLOWED_STRUCTURAL_CONDUCTIVITIES_W_M_K
        ):
            raise ValueError("P1EF conductivity must be the 0.24 or 0.29 source prior")
        if self.dry.thermal_conductivity.value_w_m_k != (
            self.structural_conductivity_w_m_k
        ):
            raise ValueError("dry structural conductivity differs from the declared case")
        if self.wet.wet.conductivity != self.structural_conductivity_w_m_k:
            raise ValueError("wet structural conductivity differs from the declared case")
        if self.dry.pressure_pa != self.process_pressure_pa:
            raise ValueError("dry pressure differs from the declared process pressure")
        if self.wet.wet.pressure_pa != self.process_pressure_pa:
            raise ValueError("wet pressure differs from the declared process pressure")
        if self.cut.dry is not self.dry or self.cut.wet is not self.wet:
            raise ValueError("cut configuration must reuse the exact matched wet/dry objects")


def _fixed_binary_selection(value_m2_s: float) -> tc.EffectiveBinaryDiffusivitySelection:
    return tc.EffectiveBinaryDiffusivityInterval(
        value_m2_s,
        value_m2_s,
        (
            "P1EF Cardarelli-scale pore-effective numerical/model-form coordinate; "
            "not soybean-DT identification"
        ),
    ).select(0.0)


def select_engineering_coefficient_case(
    identifier: str,
    binary_diffusivity_m2_s: float,
    apparent_water_diffusivity_m2_s: float,
) -> EngineeringCoefficientCase:
    """Generate one P1EF coefficient case through the guarded manifold."""

    request = cm.LocalCoefficientManifoldRequest(
        temperature_k=DESIGN_TEMPERATURE_K,
        pressure_pa=DESIGN_PRESSURE_PA,
        y_hexane=DESIGN_Y_HEXANE,
        particle_radius_m=PRODUCTION_RADIUS_M,
        binary_diffusivity=_fixed_binary_selection(binary_diffusivity_m2_s),
        apparent_water_diffusivity=cm.ApparentWaterDiffusivityTarget(
            apparent_water_diffusivity_m2_s,
            (
                "P1EF declared apparent-water numerical/model-form bracket; "
                "not a target-property confidence interval"
            ),
        ),
    )
    partition = cm.select_local_coefficient_partition(request)
    return EngineeringCoefficientCase(
        identifier=identifier,
        binary_diffusivity_m2_s=binary_diffusivity_m2_s,
        apparent_water_diffusivity_m2_s=apparent_water_diffusivity_m2_s,
        partition=partition,
    )


def engineering_coefficient_cases() -> tuple[EngineeringCoefficientCase, ...]:
    """Return the frozen nominal point followed by all rectangle vertices."""

    coordinates = (
        (
            "nominal_log_center",
            NOMINAL_BINARY_DIFFUSIVITY_M2_S,
            NOMINAL_APPARENT_WATER_DIFFUSIVITY_M2_S,
        ),
        (
            "db_min_dapp_min",
            cm.BINARY_DIFFUSIVITY_MIN_M2_S,
            cm.APPARENT_WATER_DIFFUSIVITY_MIN_M2_S,
        ),
        (
            "db_min_dapp_max",
            cm.BINARY_DIFFUSIVITY_MIN_M2_S,
            cm.APPARENT_WATER_DIFFUSIVITY_MAX_M2_S,
        ),
        (
            "db_max_dapp_min",
            cm.BINARY_DIFFUSIVITY_MAX_M2_S,
            cm.APPARENT_WATER_DIFFUSIVITY_MIN_M2_S,
        ),
        (
            "db_max_dapp_max",
            cm.BINARY_DIFFUSIVITY_MAX_M2_S,
            cm.APPARENT_WATER_DIFFUSIVITY_MAX_M2_S,
        ),
    )
    return tuple(select_engineering_coefficient_case(*item) for item in coordinates)


def build_engineering_particle_configuration(
    coefficient_case: EngineeringCoefficientCase,
    *,
    structural_conductivity_w_m_k: float = (
        NOMINAL_STRUCTURAL_CONDUCTIVITY_W_M_K
    ),
    process_pressure_pa: float = NOMINAL_PROCESS_PRESSURE_PA,
    mass_force_mode: nip.NonisothermalMassForceMode = (
        nip.NonisothermalMassForceMode.COMMON_FACE_TEMPERATURE_NEGLIGIBLE_SORET
    ),
    complete_potential_reference_gauge: nip.CompletePotentialReferenceGauge = (
        nip.DEFAULT_COMPLETE_POTENTIAL_REFERENCE_GAUGE
    ),
    binary_thermal_diffusion_factor: nip.BinaryThermalDiffusionFactorSelection = (
        nip.ZERO_BINARY_THERMAL_DIFFUSION_FACTOR
    ),
    retained_water_thermal_force_factor: (
        nip.RetainedWaterThermalForceFactorSelection
    ) = nip.ZERO_RETAINED_WATER_THERMAL_FORCE_FACTOR,
) -> EngineeringParticleConfiguration:
    """Construct matched conservative wet/dry/front authorities for P1EF."""

    if not isinstance(coefficient_case, EngineeringCoefficientCase):
        raise TypeError("P1EF configuration requires an engineering coefficient case")
    if structural_conductivity_w_m_k not in ALLOWED_STRUCTURAL_CONDUCTIVITIES_W_M_K:
        raise ValueError("P1EF conductivity must be exactly 0.24 or 0.29 W/(m K)")
    if not math.isfinite(process_pressure_pa) or not 101_000.0 <= process_pressure_pa <= 170_000.0:
        raise ValueError("P1EF process pressure must lie in [101, 170] kPa")
    if nip.reference_invariant_thermal_factor_mode(mass_force_mode):
        if (
            binary_thermal_diffusion_factor.alpha_t
            not in nip.P1EF_GAS_REPRESENTATIVE_ALPHA_T_VALUES
        ):
            raise ValueError(
                "P1EF gas thermal-factor ablation is frozen to {-0.5,0,+0.5}"
            )
        retained_allowed = (
            *nip.P1EF_RETAINED_REPRESENTATIVE_ALPHA_VALUES,
            *nip.P1EF_RETAINED_ANOMALOUS_STRESS_ALPHA_VALUES,
        )
        if retained_water_thermal_force_factor.alpha_ret not in retained_allowed:
            raise ValueError(
                "P1EF retained thermal-factor ablation is frozen to "
                "representative {-25,0,+25} plus anomalous {-40,+40}"
            )

    partition = coefficient_case.partition
    pore = partition.request.pore
    retained = partition.as_retained_water_mobility_selection()
    dry = ct.FullyDryTransportConfig(
        pressure_pa=process_pressure_pa,
        binary_diffusivity=partition.request.binary_diffusivity,
        retained_water_mobility=retained,
        thermal_conductivity=ct.ThermalConductivitySelection(
            structural_conductivity_w_m_k,
            (
                "P1EF shared wet/dry constant source prior; correlated 0.24/0.29 "
                "W/(m K) structural sensitivity"
            ),
        ),
        primitive_band=ct.GasOnlyPrimitiveDomain(
            (323.15, 433.0),
            (
                "P1EF local property envelope; exact pressure-conditioned phase "
                "topology remains authoritative"
            ),
        ),
        pore=pore,
        nonlinear_residual_tolerance=2.0e-11,
        ledger_tolerance=1.0e-10,
        nonlinear_step_tolerance=1.0e-12,
        maximum_function_evaluations=800,
    )
    wet_params = replace(
        wet_core.WetCoreParams(),
        pressure_pa=process_pressure_pa,
        conductivity=structural_conductivity_w_m_k,
        rho_dm_p=pore.rho_dm_p,
        epsilon_p=pore.eps_g,
        cp_dry_meal=pore.cp_dry_meal,
        cp_oil=pore.cp_oil,
        T_ref_solid=pore.T_ref_solid,
        w_o=pore.w_o,
        gab=pore.gab,
        oil=pore.oil,
    )
    wet = ww.WetWaterModel(
        retained,
        wet=wet_params,
        luikov=pore.luikov,
        ledger_tolerance=1.0e-10,
        nonlinear_tolerance=2.0e-11,
        max_nonlinear_evaluations=2500,
    )
    combined = cut.CutTransportConfig(
        dry,
        wet,
        cut.InterfaceCompositionBracket(
            (0.50, 0.999),
            (
                "P1EF saturated-interface root bracket; phase/fugacity equations "
                "remain unchanged"
            ),
        ),
        mass_force_mode=mass_force_mode,
        complete_potential_reference_gauge=complete_potential_reference_gauge,
        binary_thermal_diffusion_factor=binary_thermal_diffusion_factor,
        retained_water_thermal_force_factor=(
            retained_water_thermal_force_factor
        ),
        moving_interface_composition_force_authority=(
            cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH
        ),
    )
    return EngineeringParticleConfiguration(
        coefficient_case=coefficient_case,
        structural_conductivity_w_m_k=structural_conductivity_w_m_k,
        process_pressure_pa=process_pressure_pa,
        dry=dry,
        wet=wet,
        cut=combined,
    )


def source_correlated_film_pair() -> bf.CorrelatedBedFilmPair:
    """Reproduce the single frozen Coletto B.7--B.10 source-state pair."""

    source_mass_flow_kg_s = 3.9
    source_tower_diameter_m = 6.0
    gas_density_kg_m3 = 0.6
    cross_section_m2 = math.pi * source_tower_diameter_m**2 / 4.0
    superficial_velocity_m_s = source_mass_flow_kg_s / (
        gas_density_kg_m3 * cross_section_m2
    )
    state = bf.PackedBedFilmState(
        particle_radius_m=PRODUCTION_RADIUS_M,
        bed_void_fraction=0.40,
        gas_density_kg_m3=gas_density_kg_m3,
        gas_superficial_velocity_m_s=superficial_velocity_m_s,
        gas_dynamic_viscosity_pa_s=1.329e-5,
        gas_heat_capacity_j_kg_k=1926.0,
        gas_thermal_conductivity_w_m_k=0.02371,
        bulk_binary_diffusivity_m2_s=(
            efr.COLETTO_NOMINAL_BINARY_DIFFUSIVITY_M2_S
        ),
        state_provenance=(
            "Coletto source prior: 3.9 kg/s, 6 m tower, rho=0.6 kg/m3, "
            "void fraction=0.40, cp=1926 J/(kg K), k=0.02371 W/(m K)"
        ),
        reynolds_convention_provenance=(
            "Faner (2008) Eq. 4.23 solid-fraction superficial-velocity "
            "convention, D20 2026-09-23 (was: Coletto B.7--B.10 "
            "voidage-corrected)"
        ),
    )
    return bf.coletto_b7_b10_pair(state)


def make_reduced_film_boundary(
    *,
    bulk_temperature_k: float,
    bulk_y_hexane: float,
    stress_multiplier: float,
    forcing_timescale_s: float,
    label: str,
) -> ct.ReducedFilmPoreBoundary:
    """Construct a preclassified fast-mass/finite-heat P1EF boundary.

    The nominal/anomaly classification is determined from ``D/h_M**2`` before
    a nonlinear solve.  It therefore cannot become an equation-changing
    fallback selected by a rejected trial.
    """

    pair = source_correlated_film_pair()
    if stress_multiplier not in COMPACT_FILM_STRESS_MULTIPLIERS:
        raise ValueError("P1EF film multiplier must be exactly 0.5, 0.6264, 1, or 2")
    if not math.isfinite(forcing_timescale_s) or forcing_timescale_s <= 0.0:
        raise ValueError("P1EF forcing timescale must be positive and finite")
    scaled_mass = stress_multiplier * pair.binary_mass_transfer_m_s
    relaxation_s = pair.state.bulk_binary_diffusivity_m2_s / scaled_mass**2
    ratio = relaxation_s / forcing_timescale_s
    if ratio <= 0.01:
        fidelity = efr.FidelityBand.NOMINAL
    elif ratio <= 0.10:
        fidelity = efr.FidelityBand.BOUNDED_ANOMALY
    else:
        raise efr.ReducedFilmValidityError(
            "film relaxation/forcing ratio exceeds the bounded P1EF anomaly limit"
        )
    return ct.ReducedFilmPoreBoundary(
        temperature_k=bulk_temperature_k,
        y_hexane=bulk_y_hexane,
        correlated_pair=pair,
        stress_multiplier=stress_multiplier,
        liquid_saturation=0.0,
        forcing_timescale_s=forcing_timescale_s,
        absolute_mass_fraction_departure_threshold=(
            FAST_FILM_MASS_FRACTION_DEPARTURE_LIMIT
        ),
        absolute_ackermann_factor_departure_threshold=(
            UNITY_ACKERMANN_FACTOR_DEPARTURE_LIMIT
        ),
        legacy_engineering_regression_water_heat_capacity_j_kg_k=(
            LEGACY_ENGINEERING_REGRESSION_WATER_HEAT_CAPACITY_J_KG_K
        ),
        legacy_engineering_regression_hexane_heat_capacity_j_kg_k=(
            LEGACY_ENGINEERING_REGRESSION_HEXANE_HEAT_CAPACITY_J_KG_K
        ),
        binary_gas_interaction_k_wh=bg.K_WH_CENTRAL,
        label=label,
        fidelity_band=fidelity,
    )


__all__ = [
    "ALLOWED_STRUCTURAL_CONDUCTIVITIES_W_M_K",
    "COMPACT_FILM_STRESS_MULTIPLIERS",
    "CONDITIONAL_RETAINED_SURFACE_COEFFICIENT_FLOOR_M_S",
    "DESIGN_PRESSURE_PA",
    "DESIGN_TEMPERATURE_K",
    "DESIGN_Y_HEXANE",
    "EngineeringCoefficientCase",
    "EngineeringParticleConfiguration",
    "FAST_FILM_MASS_FRACTION_DEPARTURE_LIMIT",
    "LEGACY_ENGINEERING_REGRESSION_HEXANE_HEAT_CAPACITY_J_KG_K",
    "LEGACY_ENGINEERING_REGRESSION_WATER_HEAT_CAPACITY_J_KG_K",
    "MAXIMUM_RETAINED_SURFACE_RESISTANCE_FRACTION",
    "NOMINAL_APPARENT_WATER_DIFFUSIVITY_M2_S",
    "NOMINAL_BINARY_DIFFUSIVITY_M2_S",
    "NOMINAL_PROCESS_PRESSURE_PA",
    "NOMINAL_STRUCTURAL_CONDUCTIVITY_W_M_K",
    "P1EF_DECISION_ID",
    "PRODUCTION_RADIUS_M",
    "STRUCTURAL_CONDUCTIVITY_SENSITIVITY_W_M_K",
    "UNITY_ACKERMANN_FACTOR_DEPARTURE_LIMIT",
    "build_engineering_particle_configuration",
    "engineering_coefficient_cases",
    "make_reduced_film_boundary",
    "select_engineering_coefficient_case",
    "source_correlated_film_pair",
]
