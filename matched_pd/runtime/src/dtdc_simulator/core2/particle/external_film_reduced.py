r"""Bounded reduced external water/n-hexane film for the Phase-1 particle.

This module is deliberately smaller than a general non-isothermal
Maxwell--Stefan film.  Within a declared DTDC operating envelope it uses one
ideal-binary, constant-property mass-coordinate Stefan relation.  The
relation has rank one: it partitions a *supplied/coupled* total mass flux into
water and n-hexane fluxes, but it never invents the remaining total-flux
degree of freedom.

All fluxes are positive from the particle interface toward the bed gas.  The
heat flux ``q_in`` is positive in the opposite direction, from bulk gas to
particle.  Consequently the common-datum energy flux leaving the particle is

``J_E = -q_in + N_w*hbar_w + N_h*hbar_h``.

The component partial molar enthalpies occur exactly once.  No separate
latent-heat term belongs in this assembly.

The reduction is fail-closed.  It accepts only the primary-drainage DTDC band
recorded below, the inseparable Coletto B.7--B.10 coefficient pair, and one
common coefficient stress multiplier.  Invalid logarithms, inadmissible
profiles, and band violations are rejected; none is clipped or epsilon
regularized.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

from dtdc_simulator.core2.particle import bed_film as bf
from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import water as wa


COLETTO_NOMINAL_BINARY_DIFFUSIVITY_M2_S = 1.33e-5
LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT = 0.005
# GT-PS-2-OWNER-RULING-FILM-MULTIPLIER-RANGE-2026-09-30-01 (owner, 2026-09-30):
# the lower corner of the film-multiplier range of the composition screen and of
# item 08 is the smallest multiplier the reduced film's own validity check admits
# (engineering_foundation.make_reduced_film_boundary: relaxation D/h_M**2 over the
# 0.075 s forcing at or below the bounded-anomaly limit 0.10, the limit not moved).
# Measured by bisection of that check between 0.5 and 1.0 on all 19 screened states
# (the three frozen candidates on the five coefficient cases and on conductivity
# 0.29, and item 08's nominal state): the smallest admitted binary64 multiplier is
# 0.6263060583204717 on every state (the ratio is state-independent: one frozen
# Coletto source-state film pair), where the ratio is 0.1 to the last bit; the
# declared corner is that value rounded up to 1e-4, ratio 0.09997000809857089, so
# that a one-ulp host difference in the film pair cannot flip the check.  The
# former corner 0.5 (ratio 0.15690371147557045 on the D20 tree, refused since D20)
# stays admissible to the grid so that it still refuses at the validity check.
REDUCED_FILM_VALIDITY_LOWER_CORNER_MULTIPLIER = 0.6264
REDUCED_FILM_VALIDITY_LOWER_CORNER_DECLARATION = {
    "provenance_class": "DECLARED_ENGINEERING_ASSUMPTION",
    "ruling_id": "GT-PS-2-OWNER-RULING-FILM-MULTIPLIER-RANGE-2026-09-30-01",
    "criterion": (
        "the reduced external film's quasi-steady validity: relaxation D/h_M**2 over "
        "the forcing timescale at or below the bounded-anomaly limit 0.10"
    ),
    "value": 0.6264,
    "measured_smallest_admitted_multiplier": 0.6263060583204717,
    "rounding": "the measured smallest admitted multiplier rounded up to 1e-4",
    "bracket": {
        "ratio_at_former_corner_0p5": 0.15690371147557045,
        "ratio_at_measured_smallest_admitted": 0.1,
        "ratio_at_declared_corner": 0.09997000809857089,
        "ratio_at_nominal_1p0": 0.03922592786889261,
    },
    "rejecting_outcome": (
        "the validity check refusing at the declared corner, or the composition "
        "screen refusing at another corner"
    ),
    "former_corner": 0.5,
    "validity_limit_moved": False,
    "full_film_model_substituted": False,
}
_ALLOWED_STRESS_MULTIPLIERS = frozenset(
    (0.5, REDUCED_FILM_VALIDITY_LOWER_CORNER_MULTIPLIER, 1.0, 2.0)
)
ALLOWED_STRESS_MULTIPLIERS = _ALLOWED_STRESS_MULTIPLIERS


class ReducedFilmValidityError(ValueError):
    """The requested state lies outside the bounded engineering reduction."""


class FilmProfileAdmissibilityError(ValueError):
    """The proposed fluxes do not define an admissible binary film profile."""


class PureComponentBranch(str, Enum):
    WATER = "water"
    N_HEXANE = "n-hexane"


class DrainageRegime(str, Enum):
    PRIMARY_DRAINAGE = "primary_drainage"


class FidelityBand(str, Enum):
    NOMINAL = "nominal_dtdc"
    BOUNDED_ANOMALY = "bounded_dtdc_anomaly"


class FilmComponent(str, Enum):
    WATER = "water"
    N_HEXANE = "n-hexane"


class SuppliedFluxProfileAuditBranch(str, Enum):
    """Numerical representation used to certify a supplied-flux profile."""

    EXACT_UNIFORM_CONVECTIVE = "exact_uniform_convective"
    ADJACENT_UNIFORM_CONVECTIVE_LIMIT = "adjacent_uniform_convective_limit"
    FINITE_STEFAN_LOGARITHM = "finite_stefan_logarithm"


@dataclass(frozen=True)
class EquationRankMetadata:
    """Structural rank, independent of a particular numerical state."""

    component_flux_unknown_count: int
    independent_mass_residual_count: int
    remaining_total_flux_degrees_of_freedom: int
    remaining_closure: str
    total_flux_is_fabricated_by_film: bool = field(default=False, init=False)


IDEAL_BINARY_MASS_RANK = EquationRankMetadata(
    component_flux_unknown_count=2,
    independent_mass_residual_count=1,
    remaining_total_flux_degrees_of_freedom=1,
    remaining_closure="coupled component inventories and common-datum energy balance",
)

PURE_COMPONENT_BRANCH_RANK = EquationRankMetadata(
    component_flux_unknown_count=1,
    independent_mass_residual_count=0,
    remaining_total_flux_degrees_of_freedom=1,
    remaining_closure="caller-supplied total flux from the coupled balance",
)


@dataclass(frozen=True)
class ReducedFilmOperatingPoint:
    """DTDC state against which the compact-film assumptions are checked."""

    bulk_temperature_k: float
    interface_temperature_k: float
    pressure_pa: float
    particle_radius_m: float
    liquid_saturation: float
    forcing_timescale_s: float
    drainage_regime: DrainageRegime
    fidelity_band: FidelityBand = FidelityBand.NOMINAL
    components: tuple[FilmComponent, FilmComponent] = (
        FilmComponent.WATER,
        FilmComponent.N_HEXANE,
    )

    def __post_init__(self) -> None:
        scalar_values = (
            self.bulk_temperature_k,
            self.interface_temperature_k,
            self.pressure_pa,
            self.particle_radius_m,
            self.liquid_saturation,
            self.forcing_timescale_s,
        )
        if not all(math.isfinite(value) for value in scalar_values):
            raise ReducedFilmValidityError("reduced-film operating inputs must be finite")
        if not (
            323.15 <= self.bulk_temperature_k <= 433.0
            and 323.15 <= self.interface_temperature_k <= 433.0
        ):
            raise ReducedFilmValidityError("both film temperatures must lie in [323.15, 433] K")
        if not 101_000.0 <= self.pressure_pa <= 170_000.0:
            raise ReducedFilmValidityError("film pressure must lie in [101, 170] kPa")
        if not 0.5e-3 <= self.particle_radius_m <= 1.5e-3:
            raise ReducedFilmValidityError("particle radius must lie in [0.5, 1.5] mm")
        if not 0.0 <= self.liquid_saturation <= 0.80:
            raise ReducedFilmValidityError("liquid saturation must lie in [0, 0.80]")
        if self.forcing_timescale_s <= 0.0:
            raise ReducedFilmValidityError("forcing timescale must be strictly positive")
        if self.drainage_regime is not DrainageRegime.PRIMARY_DRAINAGE:
            raise ReducedFilmValidityError("the reduced film is restricted to primary drainage")
        if self.fidelity_band not in (FidelityBand.NOMINAL, FidelityBand.BOUNDED_ANOMALY):
            raise ReducedFilmValidityError("an explicit supported fidelity band is required")
        if len(self.components) != 2 or set(self.components) != {
            FilmComponent.WATER,
            FilmComponent.N_HEXANE,
        }:
            raise ReducedFilmValidityError("the reduced film is binary water/n-hexane only")

        film_temperature = 0.5 * (self.bulk_temperature_k + self.interface_temperature_k)
        relative_temperature_jump = (
            abs(self.bulk_temperature_k - self.interface_temperature_k) / film_temperature
        )
        if relative_temperature_jump > 0.20:
            raise ReducedFilmValidityError("relative bulk/interface temperature jump exceeds 0.20")


@dataclass(frozen=True)
class ReducedFilmCoefficients:
    """One jointly scaled and envelope-checked heat/mass coefficient pair."""

    base_pair: bf.CorrelatedBedFilmPair
    operating_point: ReducedFilmOperatingPoint
    stress_multiplier: float
    heat_transfer_w_m2_k: float
    binary_mass_transfer_m_s: float
    film_density_kg_m3: float
    binary_diffusivity_m2_s: float
    film_relaxation_time_s: float
    relaxation_to_forcing_ratio: float
    equation_rank: EquationRankMetadata = field(default=IDEAL_BINARY_MASS_RANK, init=False)
    validity_band_checked: bool = field(default=True, init=False)

    def __post_init__(self) -> None:
        if self.stress_multiplier not in _ALLOWED_STRESS_MULTIPLIERS:
            raise ReducedFilmValidityError(
                "stress multiplier must be exactly one of {0.5, 0.6264, 1, 2}"
            )
        expected_heat = self.stress_multiplier * self.base_pair.heat_transfer_w_m2_k
        expected_mass = self.stress_multiplier * self.base_pair.binary_mass_transfer_m_s
        if self.heat_transfer_w_m2_k != expected_heat:
            raise ReducedFilmValidityError("heat coefficient was separated from the common pair")
        if self.binary_mass_transfer_m_s != expected_mass:
            raise ReducedFilmValidityError("mass coefficient was separated from the common pair")
        if self.film_density_kg_m3 != self.base_pair.state.gas_density_kg_m3:
            raise ReducedFilmValidityError("film density must come from the common pair state")
        if self.binary_diffusivity_m2_s != self.base_pair.state.bulk_binary_diffusivity_m2_s:
            raise ReducedFilmValidityError(
                "binary diffusivity must come from the common pair state"
            )
        expected_relaxation = self.binary_diffusivity_m2_s / self.binary_mass_transfer_m_s**2
        if self.film_relaxation_time_s != expected_relaxation:
            raise ReducedFilmValidityError("film relaxation time is inconsistent with D/h_M^2")
        expected_ratio = expected_relaxation / self.operating_point.forcing_timescale_s
        if self.relaxation_to_forcing_ratio != expected_ratio:
            raise ReducedFilmValidityError("film/forcing timescale ratio is inconsistent")

    @property
    def common_provenance(self) -> str:
        return (
            f"{self.base_pair.common_provenance}; compact ideal-binary mass-coordinate film; "
            f"common lambda={self.stress_multiplier:g}; {self.operating_point.fidelity_band.value}"
        )

    @property
    def coefficient_covariance_preserved(self) -> bool:
        return (
            self.heat_transfer_w_m2_k / self.base_pair.heat_transfer_w_m2_k
            == self.binary_mass_transfer_m_s / self.base_pair.binary_mass_transfer_m_s
            == self.stress_multiplier
        )


@dataclass(frozen=True)
class ComponentMassFluxes:
    """Outward component mass fluxes, kg m^-2 s^-1."""

    water_kg_m2_s: float
    hexane_kg_m2_s: float

    def __post_init__(self) -> None:
        if not all(math.isfinite(value) for value in (self.water_kg_m2_s, self.hexane_kg_m2_s)):
            raise ValueError("component mass fluxes must be finite")

    @property
    def total_kg_m2_s(self) -> float:
        return self.water_kg_m2_s + self.hexane_kg_m2_s


@dataclass(frozen=True)
class StefanMassResidual:
    """One evaluated ideal-binary film equation and its admissibility trace."""

    residual: float
    residual_form: str
    total_mass_flux_kg_m2_s: float
    mass_peclet: float
    implied_bulk_hexane_mass_fraction: float
    log_numerator_kg_m2_s: float | None
    log_denominator_kg_m2_s: float | None
    equation_rank: EquationRankMetadata = field(default=IDEAL_BINARY_MASS_RANK, init=False)


@dataclass(frozen=True)
class ConditionedStefanMassResidual:
    """Solver-facing flux-coordinate residual plus its log-identity audit.

    ``residual_kg_m2_s`` is algebraically root-equivalent to the logarithmic
    film equation wherever that equation is admissible, but retains finite
    sensitivity as total flux crosses zero.  ``logarithmic_audit`` is a
    diagnostic representation of the same single equation, not a second
    independent residual.
    """

    residual_kg_m2_s: float
    total_mass_flux_kg_m2_s: float
    mass_peclet: float
    logarithmic_audit: StefanMassResidual
    equation_rank: EquationRankMetadata = field(default=IDEAL_BINARY_MASS_RANK, init=False)
    logarithmic_audit_adds_equation_rank: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FastFilmCompositionDepartureBound:
    """Caller-thresholded low-flux demand; no transfer rate is inferred."""

    supplied_component: FilmComponent
    supplied_internal_component_flux_kg_m2_s: float
    film_mass_conductance_kg_m2_s: float
    signed_hexane_mass_fraction_departure: float
    absolute_composition_departure_bound: float
    caller_acceptance_threshold_mass_fraction: float
    dimensionless_demand_to_threshold_ratio: float
    required_film_mass_conductance_kg_m2_s: float
    available_to_required_conductance_ratio: float | None
    within_caller_threshold: bool
    internal_flux_was_supplied: bool = field(default=True, init=False)
    film_determined_internal_flux: bool = field(default=False, init=False)


@dataclass(frozen=True)
class SuppliedFluxCompositionDeparture:
    """Exact film-implied surface departure for supplied component fluxes."""

    bulk_hexane_mass_fraction: float
    implied_surface_hexane_mass_fraction: float
    signed_surface_departure: float
    absolute_surface_departure: float
    caller_acceptance_threshold_mass_fraction: float
    within_caller_threshold: bool
    total_mass_flux_kg_m2_s: float
    mass_peclet: float
    conditioned_residual_kg_m2_s: float
    profile_audit_branch: SuppliedFluxProfileAuditBranch
    asymptotic_component_bound: FastFilmCompositionDepartureBound
    dimensionless_total_flux_demand: float
    supplied_fluxes_determine_departure_not_rate: bool = field(default=True, init=False)
    profile_admissibility_audited: bool = field(default=True, init=False)


@dataclass(frozen=True)
class AckermannHeatTransfer:
    """Stable component-capacity Ackermann correction and inward heat flux."""

    component_capacity_flux_w_m2_k: float
    beta: float
    factor: float
    heat_into_particle_w_m2: float


@dataclass(frozen=True)
class PartialEnthalpySecantHeatCapacities:
    """PHY-053 component partial-enthalpy slopes at fixed ``(P, y)``.

    Away from equal endpoint temperatures the values are exact endpoint
    secants of the frozen PHY-053 partial enthalpies.  At equal temperatures a
    deterministic symmetric five-point derivative supplies the continuous
    limit; the fixed dyadic stencil is reported rather than hidden.
    """

    surface_temperature_k: float
    bulk_temperature_k: float
    pressure_pa: float
    y_water: float
    y_hexane: float
    k_wh: float
    water_j_kg_k: float
    hexane_j_kg_k: float
    equal_temperature_derivative_limit_used: bool
    derivative_stencil_half_width_k: float | None
    production_authority: str = field(
        default="PHY-053 exact component partial-enthalpy endpoint secant",
        init=False,
    )


@dataclass(frozen=True)
class OutwardEnergyFlux:
    """Common-datum energy transfer and exact particle/bed pair ledger."""

    ackermann: AckermannHeatTransfer
    water_molar_flux_mol_m2_s: float
    hexane_molar_flux_mol_m2_s: float
    conductive_outward_w_m2: float
    water_enthalpy_outward_w_m2: float
    hexane_enthalpy_outward_w_m2: float
    total_outward_w_m2: float
    particle_inventory_source_w_m2: float
    bed_inventory_source_w_m2: float
    closed_pair_residual_w_m2: float


def adapt_correlated_bed_film(
    pair: bf.CorrelatedBedFilmPair,
    operating_point: ReducedFilmOperatingPoint,
    *,
    stress_multiplier: float = 1.0,
) -> ReducedFilmCoefficients:
    """Adapt an inseparable Coletto pair inside the bounded DTDC envelope."""

    if not isinstance(pair, bf.CorrelatedBedFilmPair):
        raise TypeError("the reduced film requires one CorrelatedBedFilmPair")
    if not isinstance(operating_point, ReducedFilmOperatingPoint):
        raise TypeError("a ReducedFilmOperatingPoint is required")
    if isinstance(stress_multiplier, bool) or not math.isfinite(stress_multiplier):
        raise ReducedFilmValidityError("stress multiplier must be a finite numeric value")
    stress_multiplier = float(stress_multiplier)
    if stress_multiplier not in _ALLOWED_STRESS_MULTIPLIERS:
        raise ReducedFilmValidityError(
            "stress multiplier must be exactly one of {0.5, 0.6264, 1, 2}"
        )

    state = pair.state
    regenerated_pair = bf.coletto_b7_b10_pair(state)
    if (
        pair.reynolds_voidage,
        pair.prandtl,
        pair.schmidt,
        pair.nusselt_voidage,
        pair.heat_transfer_w_m2_k,
        pair.binary_mass_transfer_m_s,
    ) != (
        regenerated_pair.reynolds_voidage,
        regenerated_pair.prandtl,
        regenerated_pair.schmidt,
        regenerated_pair.nusselt_voidage,
        regenerated_pair.heat_transfer_w_m2_k,
        regenerated_pair.binary_mass_transfer_m_s,
    ):
        raise ReducedFilmValidityError(
            "coefficient pair is not the inseparable B.7--B.10 result of its state"
        )
    if state.gas_superficial_velocity_m_s <= 0.0:
        raise ReducedFilmValidityError("the zero-flow correlation limit is outside this film")
    if pair.heat_transfer_w_m2_k <= 0.0 or pair.binary_mass_transfer_m_s <= 0.0:
        raise ReducedFilmValidityError("the reduced film requires a positive coefficient pair")
    if not math.isclose(
        state.particle_radius_m,
        operating_point.particle_radius_m,
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        raise ReducedFilmValidityError("operating and correlated-pair particle radii differ")
    if not 10.0 <= pair.reynolds_voidage <= 300.0:
        raise ReducedFilmValidityError("voidage Reynolds number must lie in [10, 300]")
    if not 0.6 <= pair.prandtl <= 1.3:
        raise ReducedFilmValidityError("Prandtl number must lie in [0.6, 1.3]")
    # The locally frozen Coletto direct-steam state (rho=0.6 kg/m3,
    # mu=1.329e-5 Pa s, D_HW=1.33e-5 m2/s) has Sc=1.6654135.  The bounded
    # engineering band therefore extends to 1.8 so it contains its own source
    # point; this is an explicit envelope correction, not extrapolation by
    # clipping.
    if not 0.6 <= pair.schmidt <= 1.8:
        raise ReducedFilmValidityError("Schmidt number must lie in [0.6, 1.8]")
    if not math.isclose(
        state.bulk_binary_diffusivity_m2_s,
        COLETTO_NOMINAL_BINARY_DIFFUSIVITY_M2_S,
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        raise ReducedFilmValidityError("binary diffusivity must use the Coletto nominal value")

    heat = stress_multiplier * pair.heat_transfer_w_m2_k
    mass = stress_multiplier * pair.binary_mass_transfer_m_s
    relaxation = state.bulk_binary_diffusivity_m2_s / mass**2
    timescale_ratio = relaxation / operating_point.forcing_timescale_s
    maximum_ratio = 0.01 if operating_point.fidelity_band is FidelityBand.NOMINAL else 0.10
    if timescale_ratio > maximum_ratio:
        raise ReducedFilmValidityError(
            "film relaxation is not fast enough relative to the declared forcing timescale"
        )

    return ReducedFilmCoefficients(
        base_pair=pair,
        operating_point=operating_point,
        stress_multiplier=stress_multiplier,
        heat_transfer_w_m2_k=heat,
        binary_mass_transfer_m_s=mass,
        film_density_kg_m3=state.gas_density_kg_m3,
        binary_diffusivity_m2_s=state.bulk_binary_diffusivity_m2_s,
        film_relaxation_time_s=relaxation,
        relaxation_to_forcing_ratio=timescale_ratio,
    )


def _validate_mass_fractions(interface_hexane: float, bulk_hexane: float) -> None:
    if not all(math.isfinite(value) for value in (interface_hexane, bulk_hexane)):
        raise FilmProfileAdmissibilityError("film mass fractions must be finite")
    if not 0.0 <= interface_hexane <= 1.0 or not 0.0 <= bulk_hexane <= 1.0:
        raise FilmProfileAdmissibilityError("film mass fractions must lie in [0,1]")


def _inverse_exprel(value: float) -> float:
    """Return ``value/expm1(value)`` without losing the zero limit."""

    if value == 0.0:
        return 1.0
    if abs(value) < 1.0e-5:
        square = value * value
        return 1.0 - 0.5 * value + square / 12.0 - square * square / 720.0
    if value > 700.0:
        return 0.0
    return value / math.expm1(value)


def _exprel(value: float) -> float:
    """Return ``expm1(value)/value`` without losing the zero limit."""

    if value == 0.0:
        return 1.0
    if abs(value) < 1.0e-5:
        square = value * value
        return 1.0 + 0.5 * value + square / 6.0 + square * value / 24.0
    try:
        result = math.expm1(value) / value
    except OverflowError as exc:
        raise FilmProfileAdmissibilityError("mass-coordinate profile overflowed") from exc
    if not math.isfinite(result):
        raise FilmProfileAdmissibilityError("mass-coordinate profile is non-finite")
    return result


def ideal_binary_mass_coordinate_residual(
    film: ReducedFilmCoefficients,
    *,
    interface_hexane_mass_fraction: float,
    bulk_hexane_mass_fraction: float,
    fluxes: ComponentMassFluxes,
) -> StefanMassResidual:
    r"""Evaluate the one mass-coordinate Stefan residual.

    For nonzero ``j_t = j_h + j_w`` the exact residual is

    ``ln[(j_h-w_h^b*j_t)/(j_h-w_h^I*j_t)] - j_t/(rho_f*h_M)``.

    At exactly zero total flux its continuous, full-rank branch is
    ``j_h-rho_f*h_M*(w_h^I-w_h^b)``.  The latter is Coletto A.22 in the
    low-flux limit, with an exactly counter-diffusing water flux.
    """

    if not isinstance(film, ReducedFilmCoefficients):
        raise TypeError("validated ReducedFilmCoefficients are required")
    if not isinstance(fluxes, ComponentMassFluxes):
        raise TypeError("ComponentMassFluxes are required")
    _validate_mass_fractions(
        interface_hexane_mass_fraction,
        bulk_hexane_mass_fraction,
    )

    total = fluxes.total_kg_m2_s
    conductance = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
    if total == 0.0:
        implied_bulk = interface_hexane_mass_fraction - fluxes.hexane_kg_m2_s / conductance
        if not math.isfinite(implied_bulk) or not 0.0 <= implied_bulk <= 1.0:
            raise FilmProfileAdmissibilityError(
                "zero-total-flux candidate implies a profile outside [0,1]"
            )
        residual = fluxes.hexane_kg_m2_s - conductance * (
            interface_hexane_mass_fraction - bulk_hexane_mass_fraction
        )
        return StefanMassResidual(
            residual=residual,
            residual_form="exact_zero_total_mass_flux",
            total_mass_flux_kg_m2_s=0.0,
            mass_peclet=0.0,
            implied_bulk_hexane_mass_fraction=implied_bulk,
            log_numerator_kg_m2_s=None,
            log_denominator_kg_m2_s=None,
        )

    peclet = total / conductance
    denominator = fluxes.hexane_kg_m2_s - interface_hexane_mass_fraction * total
    numerator = fluxes.hexane_kg_m2_s - bulk_hexane_mass_fraction * total
    if not all(math.isfinite(value) for value in (peclet, numerator, denominator)):
        raise FilmProfileAdmissibilityError("Stefan logarithm inputs must be finite")
    if numerator == 0.0 or denominator == 0.0:
        raise FilmProfileAdmissibilityError(
            "Stefan logarithm arguments must be strictly nonzero; use an exact branch helper"
        )
    if math.copysign(1.0, numerator) != math.copysign(1.0, denominator):
        raise FilmProfileAdmissibilityError(
            "Stefan logarithm arguments must have strictly the same sign"
        )

    implied_bulk = interface_hexane_mass_fraction - (denominator / conductance) * _exprel(peclet)
    if not math.isfinite(implied_bulk) or not 0.0 <= implied_bulk <= 1.0:
        raise FilmProfileAdmissibilityError(
            "candidate implies a mass-fraction profile outside [0,1]"
        )

    relative_increment = (
        (interface_hexane_mass_fraction - bulk_hexane_mass_fraction) * total / denominator
    )
    if math.isfinite(relative_increment) and relative_increment > -1.0:
        log_ratio = math.log1p(relative_increment)
    else:
        log_ratio = math.log(abs(numerator)) - math.log(abs(denominator))
    residual = log_ratio - peclet
    if not math.isfinite(residual):
        raise FilmProfileAdmissibilityError("Stefan residual is non-finite")
    return StefanMassResidual(
        residual=residual,
        residual_form="finite_total_mass_flux_logarithm",
        total_mass_flux_kg_m2_s=total,
        mass_peclet=peclet,
        implied_bulk_hexane_mass_fraction=implied_bulk,
        log_numerator_kg_m2_s=numerator,
        log_denominator_kg_m2_s=denominator,
    )


def component_fluxes_for_supplied_total_mass_flux(
    film: ReducedFilmCoefficients,
    *,
    total_mass_flux_kg_m2_s: float,
    interface_hexane_mass_fraction: float,
    bulk_hexane_mass_fraction: float,
) -> ComponentMassFluxes:
    """Partition a caller-supplied total flux without closing its magnitude."""

    if not isinstance(film, ReducedFilmCoefficients):
        raise TypeError("validated ReducedFilmCoefficients are required")
    if not math.isfinite(total_mass_flux_kg_m2_s):
        raise ValueError("supplied total mass flux must be finite")
    _validate_mass_fractions(
        interface_hexane_mass_fraction,
        bulk_hexane_mass_fraction,
    )
    if (
        interface_hexane_mass_fraction == bulk_hexane_mass_fraction
        and interface_hexane_mass_fraction in (0.0, 1.0)
        and total_mass_flux_kg_m2_s != 0.0
    ):
        component = (
            PureComponentBranch.N_HEXANE
            if interface_hexane_mass_fraction == 1.0
            else PureComponentBranch.WATER
        )
        return pure_component_fluxes_from_supplied_total_mass_flux(
            component=component,
            total_mass_flux_kg_m2_s=total_mass_flux_kg_m2_s,
        )

    conductance = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
    if total_mass_flux_kg_m2_s == 0.0:
        hexane = conductance * (interface_hexane_mass_fraction - bulk_hexane_mass_fraction)
    else:
        peclet = total_mass_flux_kg_m2_s / conductance
        diffusive_weight = conductance * _inverse_exprel(peclet)
        hexane = interface_hexane_mass_fraction * total_mass_flux_kg_m2_s + diffusive_weight * (
            interface_hexane_mass_fraction - bulk_hexane_mass_fraction
        )
    return ComponentMassFluxes(
        water_kg_m2_s=total_mass_flux_kg_m2_s - hexane,
        hexane_kg_m2_s=hexane,
    )


def conditioned_ideal_binary_mass_coordinate_residual(
    film: ReducedFilmCoefficients,
    *,
    interface_hexane_mass_fraction: float,
    bulk_hexane_mass_fraction: float,
    fluxes: ComponentMassFluxes,
) -> ConditionedStefanMassResidual:
    r"""Evaluate the full-rank flux-coordinate form of the film equation.

    With ``K=rho_f*h_M``, ``Pe=j_t/K``, and
    ``phi=Pe/expm1(Pe)``, the conditioned residual is

    ``j_h - w_h^I*j_t - K*(w_h^I-w_h^b)*phi``.

    Its exact ``j_t=0`` value is
    ``j_h-K*(w_h^I-w_h^b)``.  Thus the equation and its Jacobian do not
    disappear near counter-diffusion, unlike the unscaled logarithmic
    identity.  The latter is still evaluated for strict profile/admissibility
    auditing.
    """

    logarithmic_audit = ideal_binary_mass_coordinate_residual(
        film,
        interface_hexane_mass_fraction=interface_hexane_mass_fraction,
        bulk_hexane_mass_fraction=bulk_hexane_mass_fraction,
        fluxes=fluxes,
    )
    total = fluxes.total_kg_m2_s
    conductance = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
    peclet = total / conductance
    residual = (
        fluxes.hexane_kg_m2_s
        - interface_hexane_mass_fraction * total
        - conductance
        * (interface_hexane_mass_fraction - bulk_hexane_mass_fraction)
        * _inverse_exprel(peclet)
    )
    if not math.isfinite(residual):
        raise FilmProfileAdmissibilityError("conditioned Stefan residual is non-finite")
    return ConditionedStefanMassResidual(
        residual_kg_m2_s=residual,
        total_mass_flux_kg_m2_s=total,
        mass_peclet=peclet,
        logarithmic_audit=logarithmic_audit,
    )


def fast_film_composition_departure_bound(
    film: ReducedFilmCoefficients,
    *,
    supplied_component: FilmComponent,
    supplied_internal_component_flux_kg_m2_s: float,
    acceptance_threshold_mass_fraction: float,
) -> FastFilmCompositionDepartureBound:
    r"""Quantify the fast-film composition demand from a supplied flux.

    In the exact zero-total-flux/Coletto-A.22 asymptote,

    ``Delta w_h = j_h/(rho_f*h_M) = -j_w/(rho_f*h_M)``.

    The absolute value is the dimensionless film demand used here as a
    conservative departure bound.  The internal component flux and the
    acceptable departure are both caller inputs.  This helper therefore
    measures whether a proposed external-film campaign has enough
    conductance; it neither predicts an internal rate nor creates a closure.
    """

    if not isinstance(film, ReducedFilmCoefficients):
        raise TypeError("validated ReducedFilmCoefficients are required")
    if not isinstance(supplied_component, FilmComponent):
        raise TypeError("an exact FilmComponent is required")
    if not math.isfinite(supplied_internal_component_flux_kg_m2_s):
        raise ValueError("supplied internal component flux must be finite")
    if (
        not math.isfinite(acceptance_threshold_mass_fraction)
        or not 0.0 < acceptance_threshold_mass_fraction <= 1.0
    ):
        raise ValueError("caller acceptance threshold must lie in (0,1]")

    conductance = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
    orientation = 1.0 if supplied_component is FilmComponent.N_HEXANE else -1.0
    signed_departure = orientation * supplied_internal_component_flux_kg_m2_s / conductance
    absolute_departure = abs(signed_departure)
    demand_ratio = absolute_departure / acceptance_threshold_mass_fraction
    required_conductance = (
        abs(supplied_internal_component_flux_kg_m2_s) / acceptance_threshold_mass_fraction
    )
    available_to_required = (
        None if required_conductance == 0.0 else conductance / required_conductance
    )
    return FastFilmCompositionDepartureBound(
        supplied_component=supplied_component,
        supplied_internal_component_flux_kg_m2_s=(supplied_internal_component_flux_kg_m2_s),
        film_mass_conductance_kg_m2_s=conductance,
        signed_hexane_mass_fraction_departure=signed_departure,
        absolute_composition_departure_bound=absolute_departure,
        caller_acceptance_threshold_mass_fraction=acceptance_threshold_mass_fraction,
        dimensionless_demand_to_threshold_ratio=demand_ratio,
        required_film_mass_conductance_kg_m2_s=required_conductance,
        available_to_required_conductance_ratio=available_to_required,
        within_caller_threshold=absolute_departure <= acceptance_threshold_mass_fraction,
    )


def supplied_flux_composition_departure(
    film: ReducedFilmCoefficients,
    *,
    bulk_hexane_mass_fraction: float,
    supplied_fluxes: ComponentMassFluxes,
    acceptance_threshold_mass_fraction: float,
) -> SuppliedFluxCompositionDeparture:
    r"""Infer the exact film departure demanded by already-supplied fluxes.

    Rearranging the conditioned ideal-binary relation gives

    ``w_I = (j_h + K*phi*w_b)/(j_t + K*phi)``,

    where ``K=rho_f*h_M`` and ``phi=Pe/expm1(Pe)``.  This is an audit, not a
    rate closure: both component fluxes arrive from the coupled particle
    balance.  The full finite-Stefan result guards water/total-flow effects
    that the zero-total asymptotic ``|j_h|/K`` cannot see.
    """

    if not isinstance(film, ReducedFilmCoefficients):
        raise TypeError("validated ReducedFilmCoefficients are required")
    if not isinstance(supplied_fluxes, ComponentMassFluxes):
        raise TypeError("supplied ComponentMassFluxes are required")
    if not math.isfinite(bulk_hexane_mass_fraction) or not (
        0.0 <= bulk_hexane_mass_fraction <= 1.0
    ):
        raise FilmProfileAdmissibilityError("bulk mass fraction must lie in [0,1]")
    if (
        not math.isfinite(acceptance_threshold_mass_fraction)
        or not 0.0 < acceptance_threshold_mass_fraction <= 1.0
    ):
        raise ValueError("caller acceptance threshold must lie in (0,1]")

    conductance = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
    total = supplied_fluxes.total_kg_m2_s
    peclet = total / conductance
    phi = _inverse_exprel(peclet)
    denominator = total + conductance * phi
    numerator = supplied_fluxes.hexane_kg_m2_s + conductance * phi * bulk_hexane_mass_fraction
    if not math.isfinite(denominator) or denominator == 0.0:
        raise FilmProfileAdmissibilityError(
            "supplied fluxes make the film composition inversion singular"
        )
    implied_surface = numerator / denominator
    if not math.isfinite(implied_surface) or not 0.0 <= implied_surface <= 1.0:
        raise FilmProfileAdmissibilityError(
            "supplied component fluxes imply a surface composition outside [0,1]"
        )
    departure = implied_surface - bulk_hexane_mass_fraction
    conditioned_value = (
        supplied_fluxes.hexane_kg_m2_s - implied_surface * total - conductance * departure * phi
    )
    if not math.isfinite(conditioned_value):
        raise FilmProfileAdmissibilityError("conditioned film audit is non-finite")

    # A perfectly uniform convective profile is the analytic zero-diffusion
    # branch of the same equation; its logarithmic representation is 0/0 and
    # must not be epsilon-regularized.  Binary64 inversion can also place the
    # implied surface and bulk values on adjacent representable numbers while
    # one logarithm argument rounds to exact zero.  That is an exact numerical
    # chart boundary, not a physical finite departure.  Admit only that
    # one-ULP adjacency branch; every resolvably nonuniform profile receives
    # the strict logarithmic/profile audit.
    uniform_convective_branch = (
        departure == 0.0 and supplied_fluxes.hexane_kg_m2_s == bulk_hexane_mass_fraction * total
    )
    logarithm_numerator = supplied_fluxes.hexane_kg_m2_s - bulk_hexane_mass_fraction * total
    logarithm_denominator = supplied_fluxes.hexane_kg_m2_s - implied_surface * total
    adjacent_surface_and_bulk = (
        math.nextafter(implied_surface, bulk_hexane_mass_fraction) == bulk_hexane_mass_fraction
        or math.nextafter(bulk_hexane_mass_fraction, implied_surface) == implied_surface
    )
    adjacent_uniform_limit = (
        not uniform_convective_branch
        and adjacent_surface_and_bulk
        and (logarithm_numerator == 0.0 or logarithm_denominator == 0.0)
    )
    if uniform_convective_branch:
        profile_audit_branch = SuppliedFluxProfileAuditBranch.EXACT_UNIFORM_CONVECTIVE
    elif adjacent_uniform_limit:
        profile_audit_branch = SuppliedFluxProfileAuditBranch.ADJACENT_UNIFORM_CONVECTIVE_LIMIT
    else:
        conditioned = conditioned_ideal_binary_mass_coordinate_residual(
            film,
            interface_hexane_mass_fraction=implied_surface,
            bulk_hexane_mass_fraction=bulk_hexane_mass_fraction,
            fluxes=supplied_fluxes,
        )
        conditioned_value = conditioned.residual_kg_m2_s
        profile_audit_branch = SuppliedFluxProfileAuditBranch.FINITE_STEFAN_LOGARITHM

    asymptotic = fast_film_composition_departure_bound(
        film,
        supplied_component=FilmComponent.N_HEXANE,
        supplied_internal_component_flux_kg_m2_s=(supplied_fluxes.hexane_kg_m2_s),
        acceptance_threshold_mass_fraction=acceptance_threshold_mass_fraction,
    )
    absolute_departure = abs(departure)
    return SuppliedFluxCompositionDeparture(
        bulk_hexane_mass_fraction=bulk_hexane_mass_fraction,
        implied_surface_hexane_mass_fraction=implied_surface,
        signed_surface_departure=departure,
        absolute_surface_departure=absolute_departure,
        caller_acceptance_threshold_mass_fraction=acceptance_threshold_mass_fraction,
        within_caller_threshold=(absolute_departure <= acceptance_threshold_mass_fraction),
        total_mass_flux_kg_m2_s=total,
        mass_peclet=peclet,
        conditioned_residual_kg_m2_s=conditioned_value,
        profile_audit_branch=profile_audit_branch,
        asymptotic_component_bound=asymptotic,
        dimensionless_total_flux_demand=abs(peclet),
    )


def pure_component_fluxes_from_supplied_total_mass_flux(
    *,
    component: PureComponentBranch,
    total_mass_flux_kg_m2_s: float,
) -> ComponentMassFluxes:
    """Return an exact pure-component split for an externally supplied rate.

    At a pure endpoint the binary composition equation cannot determine the
    magnitude of the only remaining component flux.  Requiring the total
    flux as an argument makes that missing closure explicit and avoids both
    epsilon compositions and fabricated rates.
    """

    if not isinstance(component, PureComponentBranch):
        raise TypeError("an exact PureComponentBranch is required")
    if not math.isfinite(total_mass_flux_kg_m2_s):
        raise ValueError("supplied pure-component total flux must be finite")
    if component is PureComponentBranch.WATER:
        return ComponentMassFluxes(
            water_kg_m2_s=total_mass_flux_kg_m2_s,
            hexane_kg_m2_s=0.0,
        )
    return ComponentMassFluxes(
        water_kg_m2_s=0.0,
        hexane_kg_m2_s=total_mass_flux_kg_m2_s,
    )


def phy053_partial_enthalpy_secant_heat_capacities(
    *,
    surface_temperature_k: float,
    bulk_temperature_k: float,
    pressure_pa: float,
    y_hexane: float,
    k_wh: float = bg.K_WH_CENTRAL,
) -> PartialEnthalpySecantHeatCapacities:
    r"""Return exact PHY-053 component partial-enthalpy secants, J kg^-1 K^-1.

    Production Ackermann heat transfer uses

    ``cp_i = (hbar_i(T_b,P,y)-hbar_i(T_s,P,y))/(T_b-T_s)``

    at the fixed film pressure and composition.  This is a caloric slope only:
    the common-datum component partial enthalpies remain the sole transported
    enthalpy terms in the energy ledger.

    If ``T_b == T_s`` exactly, a symmetric five-point derivative of the same
    frozen enthalpy functions supplies the deterministic continuous limit.  A
    dyadic ``T*2^-12`` half-width avoids platform-dependent adaptive choices
    while remaining small on the smooth P1EF gas envelope.
    """

    values = (
        surface_temperature_k,
        bulk_temperature_k,
        pressure_pa,
        y_hexane,
        k_wh,
    )
    if not all(math.isfinite(value) for value in values):
        raise ReducedFilmValidityError("partial-enthalpy secant inputs must be finite")
    if surface_temperature_k <= 0.0 or bulk_temperature_k <= 0.0:
        raise ReducedFilmValidityError("partial-enthalpy secant temperatures must be positive")
    if pressure_pa <= 0.0:
        raise ReducedFilmValidityError("partial-enthalpy secant pressure must be positive")
    if not 0.0 <= y_hexane <= 1.0:
        raise ReducedFilmValidityError("partial-enthalpy secant composition must lie in [0,1]")
    if not bg.K_WH_MIN <= k_wh <= bg.K_WH_MAX:
        raise ReducedFilmValidityError("partial-enthalpy secant k_wh is outside PHY-053")

    y_water = 1.0 - y_hexane

    def partial_enthalpies(temperature_k: float) -> tuple[float, float]:
        return (
            bg.water_partial_enthalpy_mass(
                temperature_k,
                pressure_pa,
                y_water,
                y_hexane,
                k_wh=k_wh,
            ),
            bg.hexane_partial_enthalpy_mass(
                temperature_k,
                pressure_pa,
                y_water,
                y_hexane,
                k_wh=k_wh,
            ),
        )

    temperature_difference = bulk_temperature_k - surface_temperature_k
    if temperature_difference != 0.0:
        surface_enthalpies = partial_enthalpies(surface_temperature_k)
        bulk_enthalpies = partial_enthalpies(bulk_temperature_k)
        capacities = tuple(
            (bulk - surface) / temperature_difference
            for surface, bulk in zip(surface_enthalpies, bulk_enthalpies)
        )
        derivative_stencil = None
        equal_temperature_limit = False
    else:
        temperature = surface_temperature_k
        derivative_stencil = math.ldexp(temperature, -12)
        if temperature - 2.0 * derivative_stencil <= 0.0:
            raise ReducedFilmValidityError(
                "equal-temperature partial-enthalpy derivative stencil crossed zero"
            )
        minus_two = partial_enthalpies(temperature - 2.0 * derivative_stencil)
        minus_one = partial_enthalpies(temperature - derivative_stencil)
        plus_one = partial_enthalpies(temperature + derivative_stencil)
        plus_two = partial_enthalpies(temperature + 2.0 * derivative_stencil)
        denominator = 12.0 * derivative_stencil
        capacities = tuple(
            math.fsum((-plus_two_i, 8.0 * plus_one_i, -8.0 * minus_one_i, minus_two_i))
            / denominator
            for minus_two_i, minus_one_i, plus_one_i, plus_two_i in zip(
                minus_two,
                minus_one,
                plus_one,
                plus_two,
            )
        )
        equal_temperature_limit = True

    if not all(math.isfinite(value) and value > 0.0 for value in capacities):
        raise ReducedFilmValidityError(
            "PHY-053 partial-enthalpy secant heat capacities must be positive"
        )
    return PartialEnthalpySecantHeatCapacities(
        surface_temperature_k=surface_temperature_k,
        bulk_temperature_k=bulk_temperature_k,
        pressure_pa=pressure_pa,
        y_water=y_water,
        y_hexane=y_hexane,
        k_wh=k_wh,
        water_j_kg_k=capacities[0],
        hexane_j_kg_k=capacities[1],
        equal_temperature_derivative_limit_used=equal_temperature_limit,
        derivative_stencil_half_width_k=derivative_stencil,
    )


def ackermann_heat_transfer(
    film: ReducedFilmCoefficients,
    *,
    fluxes: ComponentMassFluxes,
    water_film_heat_capacity_j_kg_k: float,
    hexane_film_heat_capacity_j_kg_k: float,
    enforce_production_domain: bool = True,
) -> AckermannHeatTransfer:
    r"""Return ``A=beta/expm1(beta)`` using component capacity fluxes.

    ``beta=(j_w*cp_w,f+j_h*cp_h,f)/h_Q``.  This remains meaningful for
    counter-diffusion, where total mass flux alone would lose the component
    heat-capacity transport.  Production retains the frozen ``|beta|<=1``
    validity guard.  Disabling that guard is reserved for an explicitly
    labelled counterfactual evaluated from an already-resolved production
    root; it must never select or accept a production state.
    """

    if not isinstance(film, ReducedFilmCoefficients):
        raise TypeError("validated ReducedFilmCoefficients are required")
    if not isinstance(fluxes, ComponentMassFluxes):
        raise TypeError("ComponentMassFluxes are required")
    if not isinstance(enforce_production_domain, bool):
        raise TypeError("Ackermann production-domain flag must be boolean")
    capacities = (
        water_film_heat_capacity_j_kg_k,
        hexane_film_heat_capacity_j_kg_k,
    )
    if not all(math.isfinite(value) and value > 0.0 for value in capacities):
        raise ReducedFilmValidityError("component film heat capacities must be positive")
    component_capacity_flux = (
        fluxes.water_kg_m2_s * water_film_heat_capacity_j_kg_k
        + fluxes.hexane_kg_m2_s * hexane_film_heat_capacity_j_kg_k
    )
    beta = component_capacity_flux / film.heat_transfer_w_m2_k
    if not math.isfinite(beta) or (enforce_production_domain and abs(beta) > 1.0):
        raise ReducedFilmValidityError("Ackermann |beta| must not exceed 1")
    factor = _inverse_exprel(beta)
    heat_into_particle = (
        film.heat_transfer_w_m2_k
        * factor
        * (film.operating_point.bulk_temperature_k - film.operating_point.interface_temperature_k)
    )
    return AckermannHeatTransfer(
        component_capacity_flux_w_m2_k=component_capacity_flux,
        beta=beta,
        factor=factor,
        heat_into_particle_w_m2=heat_into_particle,
    )


def common_datum_outward_energy_flux(
    film: ReducedFilmCoefficients,
    *,
    fluxes: ComponentMassFluxes,
    water_film_heat_capacity_j_kg_k: float,
    hexane_film_heat_capacity_j_kg_k: float,
    water_partial_molar_enthalpy_j_mol: float,
    hexane_partial_molar_enthalpy_j_mol: float,
) -> OutwardEnergyFlux:
    """Assemble the outward energy flux and the equal/opposite pair ledger."""

    partial_enthalpies = (
        water_partial_molar_enthalpy_j_mol,
        hexane_partial_molar_enthalpy_j_mol,
    )
    if not all(math.isfinite(value) for value in partial_enthalpies):
        raise ValueError("supplied partial molar enthalpies must be finite")
    ackermann = ackermann_heat_transfer(
        film,
        fluxes=fluxes,
        water_film_heat_capacity_j_kg_k=water_film_heat_capacity_j_kg_k,
        hexane_film_heat_capacity_j_kg_k=hexane_film_heat_capacity_j_kg_k,
    )
    water_molar_flux = fluxes.water_kg_m2_s / wa.M
    hexane_molar_flux = fluxes.hexane_kg_m2_s / hx.M
    conductive_outward = -ackermann.heat_into_particle_w_m2
    water_enthalpy_outward = water_molar_flux * water_partial_molar_enthalpy_j_mol
    hexane_enthalpy_outward = hexane_molar_flux * hexane_partial_molar_enthalpy_j_mol
    total_outward = conductive_outward + water_enthalpy_outward + hexane_enthalpy_outward
    particle_source = -total_outward
    bed_source = total_outward
    return OutwardEnergyFlux(
        ackermann=ackermann,
        water_molar_flux_mol_m2_s=water_molar_flux,
        hexane_molar_flux_mol_m2_s=hexane_molar_flux,
        conductive_outward_w_m2=conductive_outward,
        water_enthalpy_outward_w_m2=water_enthalpy_outward,
        hexane_enthalpy_outward_w_m2=hexane_enthalpy_outward,
        total_outward_w_m2=total_outward,
        particle_inventory_source_w_m2=particle_source,
        bed_inventory_source_w_m2=bed_source,
        closed_pair_residual_w_m2=particle_source + bed_source,
    )


__all__ = [
    "COLETTO_NOMINAL_BINARY_DIFFUSIVITY_M2_S",
    "IDEAL_BINARY_MASS_RANK",
    "LEGACY_CP_SERIES_CONDUCTANCE_RELATIVE_ERROR_LIMIT",
    "PURE_COMPONENT_BRANCH_RANK",
    "AckermannHeatTransfer",
    "ComponentMassFluxes",
    "ConditionedStefanMassResidual",
    "DrainageRegime",
    "EquationRankMetadata",
    "FidelityBand",
    "FilmComponent",
    "FilmProfileAdmissibilityError",
    "FastFilmCompositionDepartureBound",
    "OutwardEnergyFlux",
    "PartialEnthalpySecantHeatCapacities",
    "PureComponentBranch",
    "ReducedFilmCoefficients",
    "ReducedFilmOperatingPoint",
    "ReducedFilmValidityError",
    "StefanMassResidual",
    "SuppliedFluxCompositionDeparture",
    "SuppliedFluxProfileAuditBranch",
    "ackermann_heat_transfer",
    "adapt_correlated_bed_film",
    "common_datum_outward_energy_flux",
    "conditioned_ideal_binary_mass_coordinate_residual",
    "component_fluxes_for_supplied_total_mass_flux",
    "fast_film_composition_departure_bound",
    "ideal_binary_mass_coordinate_residual",
    "phy053_partial_enthalpy_secant_heat_capacities",
    "pure_component_fluxes_from_supplied_total_mass_flux",
    "supplied_flux_composition_departure",
]
