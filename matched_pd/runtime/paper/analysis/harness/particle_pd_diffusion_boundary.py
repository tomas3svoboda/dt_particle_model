"""Explicit diffusion-cell PD boundary, distinct from Faner bed convection.

Source transport anchors: saved NIST dilute-vapour tables at 337.051 K.
Viscosity: Wilke (1950). Conductivity: the arithmetic/harmonic mixture
approximation documented by Cantera MixTransport. Gas diffusivity: Fuller
et al. (1966), equation 4 and ORIGINAL table I (water 12.7, C16.5/H1.98).
These are estimated mixture properties. No H/W mixture measurement is claimed.

The geometry alternatives are assumptions, not empirical confidence bounds:
an unconfined diffusion sphere, or a concentric equal-volume bed cell.
No forced main-steam through-bed convection is imposed in PD. The native
finite-flow component and Ackermann operators retain generated Stefan flow.
"""

from contextlib import contextmanager
from dataclasses import dataclass, field
import math

from dtdc_simulator.core2 import cell_native_caloric as caloric
from dtdc_simulator.core2.particle import bed_film as bed
from dtdc_simulator.core2.particle import external_film_reduced as film
from dtdc_simulator.core2.particle import coupled_transport as coupled
from dtdc_simulator.core2.props import hexane as hx, water as wa

NIST_T_K = 337.051
ANCHOR_TOLERANCE_K = 0.0006  # six-digit NIST output rounding, NOT a property extrapolation range


def mixture_transport(temperature_k, pressure_pa, y_hexane):
    if not math.isfinite(temperature_k) or abs(temperature_k - NIST_T_K) > ANCHOR_TOLERANCE_K:
        raise ValueError("PD property snapshot requires its sourced NIST temperature")
    if not math.isfinite(pressure_pa) or not 101000 <= pressure_pa <= 170000:
        raise ValueError("PD transport pressure outside the native project range")
    if not math.isfinite(y_hexane) or not 0 <= y_hexane <= 1:
        raise ValueError("invalid binary composition")
    fractions = (1 - y_hexane, y_hexane)
    mass = (wa.M, hx.M)
    mu = (11.0233e-6, 7.10619e-6)
    conductivity = (0.0212780, 0.0172559)
    mixed_mu = 0.0
    for i in range(2):
        denominator = 0.0
        for j in range(2):
            phi = (1 + math.sqrt(mu[i] / mu[j]) * (mass[j] / mass[i]) ** 0.25) ** 2 / math.sqrt(
                8 * (1 + mass[i] / mass[j])
            )
            denominator += fractions[j] * phi
        mixed_mu += fractions[i] * mu[i] / denominator
    arithmetic = math.fsum(x * k for x, k in zip(fractions, conductivity))
    harmonic = 1 / math.fsum(x / k for x, k in zip(fractions, conductivity))
    # Equation 4 gives cm^2/s, K, atm, g/mol. Convert ONCE to SI.
    volume_h = 6 * 16.5 + 14 * 1.98
    volume_w = 12.7
    diffusivity = (
        1e-4
        * 1e-3
        * temperature_k**1.75
        * math.sqrt(1 / (wa.M * 1000) + 1 / (hx.M * 1000))
        / ((pressure_pa / 101325) * (volume_w ** (1 / 3) + volume_h ** (1 / 3)) ** 2)
    )
    return dict(
        viscosity_pa_s=mixed_mu,
        conductivity_w_m_k=0.5 * (arithmetic + harmonic),
        binary_diffusivity_m2_s=diffusivity,
        conductivity_arithmetic=arithmetic,
        conductivity_harmonic=harmonic,
    )


@dataclass(frozen=True)
class DiffusivePDPairV1(bed.CorrelatedBedFilmPair):
    correlation_id: str = field(default="PRESCRIBED-PD-SPHERICAL-DIFFUSION-CELL-V1", init=False)
    effective_film_length_m: float = 0.0015
    outer_radius_m: float = math.inf
    geometry_assumption: str = "unconfined_sphere"

    def __post_init__(self):
        super().__post_init__()
        radius = self.state.particle_radius_m
        if self.geometry_assumption not in ("unconfined_sphere", "equal_volume_cell"):
            raise ValueError("unknown PD diffusion geometry")
        if not 0 < self.effective_film_length_m <= radius or not self.outer_radius_m > radius:
            raise ValueError("invalid spherical diffusion-cell radii")
        expected = (
            radius
            if math.isinf(self.outer_radius_m)
            else radius * (1 - radius / self.outer_radius_m)
        )
        if self.effective_film_length_m != expected:
            raise ValueError("diffusion-cell length differs from declared radii")
        if (
            self.heat_transfer_w_m2_k != self.state.gas_thermal_conductivity_w_m_k / expected
            or self.binary_mass_transfer_m_s != self.state.bulk_binary_diffusivity_m2_s / expected
        ):
            raise ValueError("PD diffusion conductances do not match their geometry/properties")


@dataclass(frozen=True)
class DiffusivePDCoefficientsV1(film.ReducedFilmCoefficients):
    validity_band_checked: bool = field(default=False, init=False)
    prescribed_diffusion_geometry_checked: bool = field(default=True, init=False)
    empirical_bed_correlation_applied: bool = field(default=False, init=False)


def pair(source, radius_m=0.0015, *, geometry="unconfined_sphere"):
    from pd_prescribed_bed_bath import PdReference

    reference = PdReference(**source["reference"])
    temperature, y = source["temperature_k"], source["y_hexane"]
    properties = mixture_transport(temperature, reference.pressure_pa, y)
    molar_mass = (1 - y) * wa.M + y * hx.M
    cp = (
        caloric.mixture_heat_capacity(temperature, reference.pressure_pa, y).value_j_mol_k
        / molar_mass
    )
    density = source["gas_density_kg_m3"]
    solid_fraction = reference.bed_dry_density_kg_m3 / 1159.65
    outer = math.inf if geometry == "unconfined_sphere" else radius_m / solid_fraction ** (1 / 3)
    length = radius_m if math.isinf(outer) else radius_m * (1 - radius_m / outer)
    state = bed.PackedBedFilmState(
        radius_m,
        1 - solid_fraction,
        density,
        source["total_flow_kg_s"] / (density * reference.tray_area_m2),
        properties["viscosity_pa_s"],
        cp,
        properties["conductivity_w_m_k"],
        properties["binary_diffusivity_m2_s"],
        "native bulk rho/cp; NIST dilute pure transport; Wilke viscosity; Cantera conductivity approximation; Fuller ORIGINAL table I diffusivity",
        "generated-flow Re is diagnostic only; no forced-through-bed convection correlation",
    )
    reynolds = (
        2
        * radius_m
        * density
        * state.gas_superficial_velocity_m_s
        / (state.gas_dynamic_viscosity_pa_s * solid_fraction)
    )
    pr = cp * state.gas_dynamic_viscosity_pa_s / state.gas_thermal_conductivity_w_m_k
    sc = state.gas_dynamic_viscosity_pa_s / (density * state.bulk_binary_diffusivity_m2_s)
    heat = state.gas_thermal_conductivity_w_m_k / length
    mass = state.bulk_binary_diffusivity_m2_s / length
    return DiffusivePDPairV1(
        state,
        reynolds,
        pr,
        sc,
        2 * radius_m / length,
        heat,
        mass,
        effective_film_length_m=length,
        outer_radius_m=outer,
        geometry_assumption=geometry,
    )


@contextmanager
def installed_pd_diffusion():
    """Dispatch only the distinctly typed new closure; old correlation unchanged."""
    original = film.adapt_correlated_bed_film
    original_prepare = coupled.prepare_reduced_film_surface_state

    def adapt(pair, operating_point, *, stress_multiplier=1.0):
        if not isinstance(pair, DiffusivePDPairV1):
            return original(pair, operating_point, stress_multiplier=stress_multiplier)
        if stress_multiplier != 1.0:
            raise ValueError(
                "PD geometry alternatives must be declared by radius, not a hidden stress factor"
            )
        state = pair.state
        if (
            operating_point.particle_radius_m != state.particle_radius_m
            or not 101000 <= operating_point.pressure_pa <= 170000
        ):
            raise ValueError("PD diffusion geometry/pressure identity differs")
        relaxation = state.bulk_binary_diffusivity_m2_s / pair.binary_mass_transfer_m_s**2
        ratio = relaxation / operating_point.forcing_timescale_s
        maximum = 0.01 if operating_point.fidelity_band is film.FidelityBand.NOMINAL else 0.10
        if ratio > maximum:
            raise film.ReducedFilmValidityError(
                "PD diffusion relaxation does not meet the unchanged quasi-steady criterion"
            )
        return DiffusivePDCoefficientsV1(
            pair,
            operating_point,
            1.0,
            pair.heat_transfer_w_m2_k,
            pair.binary_mass_transfer_m_s,
            state.gas_density_kg_m3,
            state.bulk_binary_diffusivity_m2_s,
            relaxation,
            ratio,
        )

    def prepare(boundary, **kwargs):
        if not isinstance(boundary.correlated_pair, DiffusivePDPairV1):
            return original_prepare(boundary, **kwargs)
        # This boundary freezes molecular properties at the sourced BULK
        # state. The series-conduction predictor is a geometric auxiliary,
        # not the actual interface temperature or a phase certificate.
        # A 49 C body can give a predictor below the old convective oracle's
        # 50 C calibration floor. Do not present that auxiliary as an
        # accepted thermodynamic state or relax the original oracle.
        radius = kwargs["particle_radius_m"]
        distance = kwargs["outer_half_cell_distance_m"]
        conductivity = kwargs["particle_thermal_conductivity_w_m_k"]
        cell_t = kwargs["cell_temperature_k"]
        if not all(math.isfinite(v) and v > 0 for v in (radius, distance, conductivity, cell_t)):
            raise ValueError("invalid PD reference geometry")
        if kwargs["binary_gas_interaction_k_wh"] != boundary.binary_gas_interaction_k_wh:
            raise ValueError("PD native binary caloric identity differs")
        point = film.ReducedFilmOperatingPoint(
            boundary.temperature_k,
            boundary.temperature_k,
            kwargs["pressure_pa"],
            radius,
            boundary.liquid_saturation,
            boundary.forcing_timescale_s,
            boundary.drainage_regime,
            boundary.fidelity_band,
        )
        coefficients = adapt(
            boundary.correlated_pair, point, stress_multiplier=boundary.stress_multiplier
        )
        solid_r = distance / conductivity
        gas_r = 1 / coefficients.heat_transfer_w_m2_k
        series = 1 / (solid_r + gas_r)
        q = series * (cell_t - boundary.temperature_k)
        surface_t = math.fsum((cell_t, -q * solid_r))
        return coupled.ReducedFilmSurfaceState(
            boundary,
            coefficients,
            cell_t,
            surface_t,
            boundary.y_hexane,
            coupled._hexane_mass_fraction_from_mole_fraction(boundary.y_hexane),
            distance,
            solid_r,
            gas_r,
            series,
            q,
            (surface_t - cell_t) / distance,
        )

    film.adapt_correlated_bed_film = adapt
    coupled.prepare_reduced_film_surface_state = prepare
    try:
        yield
    finally:
        film.adapt_correlated_bed_film = original
        coupled.prepare_reduced_film_surface_state = original_prepare
