r"""Gate 1e: PHY-039 moving-front jumps and local fugacity equilibrium.

Outward radial flux is positive and ``[psi]_w^d = psi_d - psi_w``.  The
continuous Rankine--Hugoniot laws are

    [j_k - s_dot C_k]_w^d = 0
    [q + sum_k h_k j_k - s_dot u]_w^d = 0.

The numerical kernel uses their exact finite-step ALE form with
``z=(s/R)^3``:

    Q_Gamma = 4*pi/3 (s_new^3-s_old^3)/dt
    H_k     = A_Gamma j_k - Q_Gamma C_k
    H_E     = A_Gamma (q + sum_k h_k j_k) - Q_Gamma u.

Analytically, ``Q_Gamma = V_R (z_new-z_old)/dt``.  The radius form is evaluated
from the front-fitted grid face and the resulting *single scalar* is shared by
the wet region, jump, and dry region.  This avoids ULP-level disagreement from
independently evaluating the equivalent ``z`` and cube-root expressions.  The
jump residual is ``H_d-H_w``.  This is the same swept-volume convention as gate
1d and avoids approximating ``A*s_dot``.  In particular, a wet-core rate
returned by ``wet_core.residual_and_jacobian`` is already a relative ALE energy
rate; it must not be treated as conductive ``q`` and have ``Q*u`` subtracted a
second time.

Faner's ``s(X)`` relation appears here only as an explicitly named validation
oracle.  It is never called by the RH front advance.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from typing import Callable, Sequence

from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp


class FrontTopologyError(ValueError):
    """The one-front primary-drainage closure is not physically admissible."""


@dataclass(frozen=True)
class FrontCoordinate:
    """Sharp-front volume coordinate with exact endpoint semantics."""

    particle_radius_m: float
    z: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.particle_radius_m) or self.particle_radius_m <= 0.0:
            raise ValueError("particle radius must be positive and finite")
        if not math.isfinite(self.z) or not 0.0 <= self.z <= 1.0:
            raise ValueError("front volume coordinate z must lie in [0, 1]")

    @property
    def radius_m(self) -> float:
        if self.z == 0.0:
            return 0.0
        if self.z == 1.0:
            return self.particle_radius_m
        return self.particle_radius_m * self.z ** (1.0 / 3.0)

    @property
    def area_m2(self) -> float:
        return 4.0 * math.pi * self.radius_m**2

    @property
    def particle_volume_m3(self) -> float:
        return 4.0 / 3.0 * math.pi * self.particle_radius_m**3

    @property
    def regime(self) -> str:
        if self.z == 0.0:
            return "dry"
        if self.z == 1.0:
            return "fully_wet"
        return "partial"

    @classmethod
    def from_radius(cls, particle_radius_m: float, front_radius_m: float) -> "FrontCoordinate":
        if not math.isfinite(front_radius_m) or not 0.0 <= front_radius_m <= particle_radius_m:
            raise ValueError("front radius must lie in [0, R]")
        if front_radius_m == 0.0:
            z = 0.0
        elif front_radius_m == particle_radius_m:
            z = 1.0
        else:
            z = (front_radius_m / particle_radius_m) ** 3
        return cls(particle_radius_m, z)


@dataclass(frozen=True)
class FrontTrace:
    """One side of the interface, using mass-based component fluxes."""

    components: tuple[str, ...]
    concentrations_kg_m3: tuple[float, ...]
    fluxes_kg_m2_s: tuple[float, ...]
    partial_enthalpies_j_kg: tuple[float, ...]
    energy_density_j_m3: float
    heat_flux_w_m2: float

    def __post_init__(self) -> None:
        lengths = (
            len(self.components),
            len(self.concentrations_kg_m3),
            len(self.fluxes_kg_m2_s),
            len(self.partial_enthalpies_j_kg),
        )
        if lengths[0] == 0 or len(set(lengths)) != 1:
            raise ValueError("front trace component arrays must be non-empty and aligned")
        if len(set(self.components)) != len(self.components):
            raise ValueError("front trace component names must be unique")
        values = (
            *self.concentrations_kg_m3,
            *self.fluxes_kg_m2_s,
            *self.partial_enthalpies_j_kg,
            self.energy_density_j_m3,
            self.heat_flux_w_m2,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("front trace values must be finite")
        if any(value < 0.0 for value in self.concentrations_kg_m3):
            raise ValueError("front concentrations must be non-negative")


@dataclass(frozen=True)
class RelativeRates:
    """Component and energy rates relative to one moving interface, SI."""

    components: tuple[str, ...]
    component_rates_kg_s: tuple[float, ...]
    energy_rate_w: float


@dataclass(frozen=True)
class JumpResidual:
    """Dry-minus-wet Rankine--Hugoniot residuals in extensive rates."""

    components: tuple[str, ...]
    component_residuals_kg_s: tuple[float, ...]
    energy_residual_w: float


@dataclass(frozen=True)
class ContinuousHexaneJump:
    """One-component continuous RH solution per unit interface area."""

    front_speed_m_s: float
    dry_hexane_flux_kg_m2_s: float
    mass_residual_kg_m2_s: float
    energy_residual_w_m2: float
    energy_denominator_j_m3: float


@dataclass(frozen=True)
class InterfaceEquilibrium:
    """Local pure-liquid/binary-pore-gas n-hexane fugacity equality."""

    temperature_k: float
    pressure_liquid_pa: float
    pressure_gas_pa: float
    y_water: float
    y_hexane: float
    liquid_fugacity_pa: float
    gas_hexane_fugacity_pa: float
    log_fugacity_residual: float
    phi_hexane: float
    gas_hexane_density_kg_m3: float
    gas_hexane_partial_enthalpy_j_kg: float
    gas_model: str
    outside_project_pressure: bool
    below_cross_direct_evidence: bool


@dataclass(frozen=True)
class EquilibriumDryStorage:
    """Fugacity-matched dry-side pore-gas plus retained n-hexane storage."""

    equilibrium: InterfaceEquilibrium
    pore_gas_hexane_density_kg_m3: float
    retained_loading: float
    total_concentration_kg_m3: float


@dataclass(frozen=True)
class FrontStep:
    """Conservative front move with an exactly located dry-out event."""

    before: FrontCoordinate
    after: FrontCoordinate
    elapsed_wet_s: float
    remaining_dry_s: float
    integrated_surface_area_time_m2_s: float
    transferred_hexane_kg: float
    component_ledger_residual_kg: float
    converted_wet_volume_m3: float
    converted_energy_difference_j: float
    extinguished: bool


def relative_rates(trace: FrontTrace, area_m2: float, swept_volume_rate_m3_s: float) -> RelativeRates:
    """Return ``A*j-Q*C`` and ``A*(q+sum h*j)-Q*u`` for one side."""
    if not math.isfinite(area_m2) or area_m2 < 0.0:
        raise ValueError("interface area must be finite and non-negative")
    if not math.isfinite(swept_volume_rate_m3_s):
        raise ValueError("interface swept-volume rate must be finite")
    component_rates = tuple(
        area_m2 * flux - swept_volume_rate_m3_s * concentration
        for flux, concentration in zip(trace.fluxes_kg_m2_s, trace.concentrations_kg_m3)
    )
    local_energy_flux = trace.heat_flux_w_m2 + sum(
        h * flux for h, flux in zip(trace.partial_enthalpies_j_kg, trace.fluxes_kg_m2_s)
    )
    energy_rate = (
        area_m2 * local_energy_flux
        - swept_volume_rate_m3_s * trace.energy_density_j_m3
    )
    return RelativeRates(trace.components, component_rates, energy_rate)


def jump_residual(
    wet: FrontTrace, dry: FrontTrace, area_m2: float, swept_volume_rate_m3_s: float
) -> JumpResidual:
    """Evaluate exact finite-step RH residuals ``H_d-H_w``."""
    if wet.components != dry.components:
        raise ValueError("wet and dry front traces must use identical component ordering")
    wet_rates = relative_rates(wet, area_m2, swept_volume_rate_m3_s)
    dry_rates = relative_rates(dry, area_m2, swept_volume_rate_m3_s)
    return JumpResidual(
        components=wet.components,
        component_residuals_kg_s=tuple(
            dry_value - wet_value
            for wet_value, dry_value in zip(
                wet_rates.component_rates_kg_s, dry_rates.component_rates_kg_s
            )
        ),
        energy_residual_w=dry_rates.energy_rate_w - wet_rates.energy_rate_w,
    )


def swept_geometry(
    before: FrontCoordinate, after: FrontCoordinate, dt: float
) -> tuple[float, float]:
    """Return backward-Euler ``(A_Gamma, Q_Gamma)`` at the new front.

    ``Q_Gamma`` is evaluated from the boundary radii, exactly as in the wet
    and dry ALE grid operators.  Coupled assemblers that already own those
    grids should call :func:`swept_volume_rate` on the actual shared face once
    and pass that scalar to every regional/interface residual.
    """
    if before.particle_radius_m != after.particle_radius_m:
        raise ValueError("front step cannot change particle radius")
    if after.z > before.z:
        raise FrontTopologyError("primary-drainage front cannot advance outward")
    Q = swept_volume_rate(before.radius_m, after.radius_m, dt)
    return after.area_m2, Q


def swept_volume_rate(
    previous_front_radius_m: float, front_radius_m: float, dt: float
) -> float:
    """Return the exact spherical volume swept per time by one grid face.

    The argument order follows a time step: previous radius, new radius.  A
    receding primary-drainage interface therefore returns a non-positive rate.
    The formula deliberately matches the ALE operators in ``wet_core`` and
    ``dry_region`` operation-for-operation.
    """
    radii = (previous_front_radius_m, front_radius_m)
    if not all(math.isfinite(radius) and radius >= 0.0 for radius in radii):
        raise ValueError("front face radii must be finite and non-negative")
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("front time step must be positive and finite")
    if front_radius_m > previous_front_radius_m:
        raise FrontTopologyError("primary-drainage front cannot advance outward")
    return (
        4.0
        / 3.0
        * math.pi
        * (front_radius_m**3 - previous_front_radius_m**3)
        / dt
    )


def solve_continuous_hexane_jump(
    *,
    wet_concentration_kg_m3: float,
    dry_concentration_kg_m3: float,
    wet_energy_density_j_m3: float,
    dry_energy_density_j_m3: float,
    wet_heat_flux_w_m2: float,
    dry_heat_flux_w_m2: float,
    dry_hexane_enthalpy_j_kg: float,
    fixed_wet_energy_flux_w_m2: float = 0.0,
    fixed_dry_energy_flux_w_m2: float = 0.0,
) -> ContinuousHexaneJump:
    """Solve the one-component continuous mass+energy jumps without a rate law.

    ``fixed_*_energy_flux`` carries any already-closed non-hexane enthalpy flux.
    There is deliberately no latent-heat input: phase/binding duty arises from
    the common-datum ``u`` jump and the transported hexane enthalpy.
    """
    values = (
        wet_concentration_kg_m3,
        dry_concentration_kg_m3,
        wet_energy_density_j_m3,
        dry_energy_density_j_m3,
        wet_heat_flux_w_m2,
        dry_heat_flux_w_m2,
        dry_hexane_enthalpy_j_kg,
        fixed_wet_energy_flux_w_m2,
        fixed_dry_energy_flux_w_m2,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("continuous front-jump inputs must be finite")
    if wet_concentration_kg_m3 < 0.0 or dry_concentration_kg_m3 < 0.0:
        raise ValueError("front concentrations must be non-negative")
    delta_c = wet_concentration_kg_m3 - dry_concentration_kg_m3
    if delta_c <= 0.0:
        raise FrontTopologyError("wet-side historical hexane concentration must exceed dry storage")
    delta_u = dry_energy_density_j_m3 - wet_energy_density_j_m3
    denominator = dry_hexane_enthalpy_j_kg * (
        dry_concentration_kg_m3 - wet_concentration_kg_m3
    ) - delta_u
    scale = max(
        abs(dry_hexane_enthalpy_j_kg * delta_c),
        abs(delta_u),
        1.0,
    )
    if abs(denominator) <= 1.0e-14 * scale:
        raise FrontTopologyError("front energy jump is singular")
    fixed_jump = (
        dry_heat_flux_w_m2
        + fixed_dry_energy_flux_w_m2
        - wet_heat_flux_w_m2
        - fixed_wet_energy_flux_w_m2
    )
    speed = -fixed_jump / denominator
    dry_flux = -speed * delta_c
    if speed > 0.0 or dry_flux < 0.0:
        raise FrontTopologyError(
            "front jump requires condensation/outward front motion on the primary-drainage branch"
        )
    mass_residual = dry_flux - speed * (
        dry_concentration_kg_m3 - wet_concentration_kg_m3
    )
    energy_residual = (
        fixed_jump
        + dry_hexane_enthalpy_j_kg * dry_flux
        - speed * delta_u
    )
    return ContinuousHexaneJump(
        front_speed_m_s=speed,
        dry_hexane_flux_kg_m2_s=dry_flux,
        mass_residual_kg_m2_s=mass_residual,
        energy_residual_w_m2=energy_residual,
        energy_denominator_j_m3=denominator,
    )


def solve_equilibrium_hexane_jump(
    equilibrium: InterfaceEquilibrium, **jump_inputs: float
) -> ContinuousHexaneJump:
    """Solve RH jumps using the partial h_h from the same gas fugacity model."""
    if "dry_hexane_enthalpy_j_kg" in jump_inputs:
        raise ValueError("equilibrium jump enthalpy is supplied by the fugacity/caloric state")
    return solve_continuous_hexane_jump(
        dry_hexane_enthalpy_j_kg=equilibrium.gas_hexane_partial_enthalpy_j_kg,
        **jump_inputs,
    )


def wet_hexane_concentration(rho_dm_p: float, historical_loading: float) -> float:
    """``C_h,w = rho_dm,p X_h,w^a`` using the immutable activation trace."""
    if not all(math.isfinite(value) for value in (rho_dm_p, historical_loading)):
        raise ValueError("wet concentration inputs must be finite")
    if rho_dm_p <= 0.0 or historical_loading < 0.0:
        raise ValueError("wet concentration inputs are outside physical bounds")
    return rho_dm_p * historical_loading


def dry_hexane_concentration(
    eps_g: float, pore_gas_hexane_density_kg_m3: float, rho_dm_p: float, retained_loading: float
) -> float:
    """``C_h,d = eps_g*c_h,g + rho_dm,p*W_h``; pore gas is never omitted."""
    values = (eps_g, pore_gas_hexane_density_kg_m3, rho_dm_p, retained_loading)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("dry concentration inputs must be finite")
    if not 0.0 <= eps_g <= 1.0 or pore_gas_hexane_density_kg_m3 < 0.0:
        raise ValueError("dry pore-gas storage inputs are outside physical bounds")
    if rho_dm_p <= 0.0 or retained_loading < 0.0:
        raise ValueError("dry retained storage inputs are outside physical bounds")
    return eps_g * pore_gas_hexane_density_kg_m3 + rho_dm_p * retained_loading


def log_fugacity_residual(
    T: float,
    pressure_gas_pa: float,
    y_water: float,
    y_hexane: float,
    *,
    pressure_liquid_pa: float | None = None,
    k_wh: float = bg.K_WH_CENTRAL,
    gas_model: str = "virial",
) -> tuple[float, float, float, float, float, float, bool]:
    """Return ``ln(f_h,l^pure)-ln(y_h phi_h P_g)`` and state diagnostics."""
    composition = (y_water, y_hexane)
    if not all(math.isfinite(value) and value >= 0.0 for value in composition):
        raise ValueError("front gas mole fractions must be finite and non-negative")
    if abs(y_water + y_hexane - 1.0) > 1.0e-12:
        raise ValueError("front gas mole fractions must sum to one")
    if y_hexane <= 0.0:
        raise FrontTopologyError("active liquid-hexane front requires positive gas hexane fraction")
    P_l = pressure_gas_pa if pressure_liquid_pa is None else pressure_liquid_pa
    if not all(math.isfinite(value) and value > 0.0 for value in (P_l, pressure_gas_pa)):
        raise ValueError("front liquid and gas pressures must be positive and finite")
    if P_l != pressure_gas_pa:
        raise ValueError(
            "PHY-037 planar nominal requires equal liquid and gas interface pressures"
        )
    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("front equilibrium temperature must be positive and finite")
    liquid = hx.state_Tp(T, P_l, "liquid")
    if gas_model == "virial":
        mixture = bg.state(T, pressure_gas_pa, y_water, y_hexane, k_wh=k_wh)
        ln_phi = mixture.ln_phi_hexane
        phi = mixture.phi_hexane
        below = mixture.below_cross_direct_evidence
        gas_density = mixture.hexane_density_kg_m3
        partial_enthalpy = (
            hx.h_ideal(T) + mixture.partial_residual_enthalpy_hexane_molar
        ) / hx.M
    elif gas_model == "ideal":
        # Mandatory PHY-053 ablation oracle, not the production branch.
        ln_phi, phi, below = 0.0, 1.0, False
        gas_density = y_hexane * pressure_gas_pa * hx.M / (bg.R * T)
        partial_enthalpy = hx.h_ideal(T) / hx.M
    else:
        raise ValueError(f"unknown front gas fugacity model {gas_model!r}")
    ln_gas = math.log(y_hexane * pressure_gas_pa) + ln_phi
    return (
        liquid.ln_fugacity - ln_gas,
        liquid.fugacity,
        math.exp(ln_gas),
        phi,
        gas_density,
        partial_enthalpy,
        below,
    )


def solve_interface_equilibrium(
    pressure_gas_pa: float,
    y_water: float,
    y_hexane: float,
    *,
    pressure_liquid_pa: float | None = None,
    temperature_bounds_k: tuple[float, float] = (310.0, 430.0),
    k_wh: float = bg.K_WH_CENTRAL,
    gas_model: str = "virial",
) -> InterfaceEquilibrium:
    """Solve the active PHY-039 fugacity equality by safeguarded bisection."""
    lo, hi = temperature_bounds_k
    if not (math.isfinite(lo) and math.isfinite(hi) and 0.0 < lo < hi):
        raise ValueError("invalid front-equilibrium temperature bracket")

    def evaluate(T: float):
        return log_fugacity_residual(
            T,
            pressure_gas_pa,
            y_water,
            y_hexane,
            pressure_liquid_pa=pressure_liquid_pa,
            k_wh=k_wh,
            gas_model=gas_model,
        )

    r_lo = evaluate(lo)[0]
    r_hi = evaluate(hi)[0]
    if r_lo == 0.0:
        root = lo
    elif r_hi == 0.0:
        root = hi
    elif r_lo * r_hi > 0.0:
        raise FrontTopologyError("liquid/binary-gas fugacity root is not bracketed")
    else:
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            residual = evaluate(mid)[0]
            if residual == 0.0:
                lo = hi = mid
                break
            if r_lo * residual <= 0.0:
                hi = mid
                r_hi = residual
            else:
                lo = mid
                r_lo = residual
            if hi - lo < 1.0e-10:
                break
        root = 0.5 * (lo + hi)
    residual, f_liquid, f_gas, phi, gas_density, partial_enthalpy, below = evaluate(root)
    P_l = pressure_gas_pa if pressure_liquid_pa is None else pressure_liquid_pa
    return InterfaceEquilibrium(
        temperature_k=root,
        pressure_liquid_pa=P_l,
        pressure_gas_pa=pressure_gas_pa,
        y_water=y_water,
        y_hexane=y_hexane,
        liquid_fugacity_pa=f_liquid,
        gas_hexane_fugacity_pa=f_gas,
        log_fugacity_residual=residual,
        phi_hexane=phi,
        gas_hexane_density_kg_m3=gas_density,
        gas_hexane_partial_enthalpy_j_kg=partial_enthalpy,
        gas_model=gas_model,
        outside_project_pressure=not (
            bg.PROJECT_PRESSURE_MIN_PA
            <= pressure_gas_pa
            <= bg.PROJECT_PRESSURE_MAX_PA
        ),
        below_cross_direct_evidence=below,
    )


def equilibrium_dry_storage(
    equilibrium: InterfaceEquilibrium,
    *,
    eps_g: float,
    rho_dm_p: float,
    w_o: float = sp.W_O_REF,
    gab: sp.GabParams = sp.GabParams(),
    oil: sp.OilIsotherm = sp.OilIsotherm(),
) -> EquilibriumDryStorage:
    """Build C_h,d at the same active fugacity-equilibrium state.

    At a liquid-bearing front ``a_h=1`` by definition.  The pore-gas density
    comes from the selected binary pressure-virial state, not gate 1c's ideal
    Cardarelli ``c/c_sat`` oracle.
    """
    retained = sp.retained_hexane(1.0, equilibrium.temperature_k, w_o, gab, oil)
    total = dry_hexane_concentration(
        eps_g,
        equilibrium.gas_hexane_density_kg_m3,
        rho_dm_p,
        retained,
    )
    return EquilibriumDryStorage(
        equilibrium=equilibrium,
        pore_gas_hexane_density_kg_m3=equilibrium.gas_hexane_density_kg_m3,
        retained_loading=retained,
        total_concentration_kg_m3=total,
    )


def exact_pure_hexane_equilibrium_temperature(
    pressure_pa: float, *, temperature_bounds_k: tuple[float, float] = (310.0, 430.0)
) -> float:
    """Exact Span--Wagner pure-fluid oracle ``Psat(T)=P``."""
    if not math.isfinite(pressure_pa) or pressure_pa <= 0.0:
        raise ValueError("pure-hexane equilibrium pressure must be positive and finite")
    lo, hi = temperature_bounds_k
    r_lo = hx.saturation_pressure(lo) - pressure_pa
    r_hi = hx.saturation_pressure(hi) - pressure_pa
    if r_lo * r_hi > 0.0:
        raise ValueError("pure-hexane saturation root is not bracketed")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        residual = hx.saturation_pressure(mid) - pressure_pa
        if r_lo * residual <= 0.0:
            hi = mid
            r_hi = residual
        else:
            lo = mid
            r_lo = residual
        if hi - lo < 1.0e-10:
            break
    return 0.5 * (lo + hi)


def backward_euler_front_step(
    before: FrontCoordinate,
    dt: float,
    concentration_jump_kg_m3: float,
    outward_flux_at_new: Callable[[FrontCoordinate], float],
) -> FrontCoordinate:
    """Solve the BE component jump for a prescribed outward dry-side flux.

    This isolated integrator is a manufactured-order qualifier.  The coupled
    sphere later solves the same mass, energy, and fugacity residuals together.
    """
    if before.z == 0.0:
        raise FrontTopologyError("dry endpoint has no active front equation")
    if dt <= 0.0 or not math.isfinite(dt) or concentration_jump_kg_m3 <= 0.0:
        raise ValueError("invalid backward-Euler front input")

    def residual(z_new: float) -> float:
        trial = FrontCoordinate(before.particle_radius_m, z_new)
        flux = outward_flux_at_new(trial)
        if not math.isfinite(flux) or flux < 0.0:
            raise FrontTopologyError("primary-drainage front requires non-negative outward flux")
        Q = before.particle_volume_m3 * (z_new - before.z) / dt
        return Q * concentration_jump_kg_m3 + trial.area_m2 * flux

    lo, hi = 0.0, before.z
    r_lo, r_hi = residual(lo), residual(hi)
    if r_lo > 0.0 or r_hi < 0.0:
        raise FrontTopologyError("backward-Euler front residual is not bracketed")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        r_mid = residual(mid)
        if r_mid <= 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-14:
            break
    return FrontCoordinate(before.particle_radius_m, 0.5 * (lo + hi))


def step_constant_flux_to_extinction(
    before: FrontCoordinate,
    dt: float,
    concentration_jump_kg_m3: float,
    outward_flux_kg_m2_s: float,
    wet_minus_dry_energy_density_j_m3: float,
) -> FrontStep:
    """Exact constant-flux Stefan trajectory with a located ``z=0`` event."""
    values = (
        dt,
        concentration_jump_kg_m3,
        outward_flux_kg_m2_s,
        wet_minus_dry_energy_density_j_m3,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("front-step inputs must be finite")
    if dt <= 0.0 or concentration_jump_kg_m3 <= 0.0 or outward_flux_kg_m2_s < 0.0:
        raise ValueError("front-step inputs are outside physical bounds")
    if before.z == 0.0:
        return FrontStep(before, before, 0.0, dt, 0.0, 0.0, 0.0, 0.0, 0.0, False)
    if outward_flux_kg_m2_s == 0.0:
        return FrontStep(before, before, dt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, False)
    speed_magnitude = outward_flux_kg_m2_s / concentration_jump_kg_m3
    extinction_time = before.radius_m / speed_magnitude
    event_tolerance = 8.0 * sys.float_info.epsilon * max(dt, extinction_time)
    if dt >= extinction_time - event_tolerance:
        # Set the event endpoint exactly.  Algebraically equivalent evaluations
        # of the event time can differ by an ulp; subtraction would then leave a
        # microscopic ghost wet core.  Keep the reported substep times
        # non-negative when ``dt`` lies just below the computed event time.
        elapsed = min(dt, extinction_time)
        radius_after = 0.0
    else:
        elapsed = dt
        radius_after = before.radius_m - speed_magnitude * elapsed
    after = FrontCoordinate.from_radius(before.particle_radius_m, radius_after)
    inventory_before = before.particle_volume_m3 * before.z * concentration_jump_kg_m3
    inventory_after = after.particle_volume_m3 * after.z * concentration_jump_kg_m3
    converted_volume = before.particle_volume_m3 * before.z - after.particle_volume_m3 * after.z
    area_time = (
        4.0
        * math.pi
        * (before.radius_m**3 - after.radius_m**3)
        / (3.0 * speed_magnitude)
    )
    # Independent surface-flux integral over s(t)=s0-(j/DeltaC)t.
    transferred = outward_flux_kg_m2_s * area_time
    return FrontStep(
        before=before,
        after=after,
        elapsed_wet_s=elapsed,
        remaining_dry_s=max(0.0, dt - extinction_time),
        integrated_surface_area_time_m2_s=area_time,
        transferred_hexane_kg=transferred,
        component_ledger_residual_kg=inventory_after - inventory_before + transferred,
        converted_wet_volume_m3=converted_volume,
        converted_energy_difference_j=(
            converted_volume * wet_minus_dry_energy_density_j_m3
        ),
        extinguished=after.z == 0.0 and before.z > 0.0,
    )


def validate_film_front_complementarity(
    attached_loading: float, coordinate: FrontCoordinate
) -> None:
    """Enforce ``X_f >= 0``, ``R-s >= 0``, and ``X_f(R-s)=0`` exactly."""
    if not math.isfinite(attached_loading) or attached_loading < 0.0:
        raise ValueError("attached n-hexane loading must be finite and non-negative")
    if attached_loading > 0.0 and coordinate.z != 1.0:
        raise FrontTopologyError("positive attached inventory requires an exactly full wet core")


def faner_radius_oracle(
    total_loading: float,
    equilibrium_loading: float,
    critical_loading: float,
    particle_radius_m: float,
) -> float:
    """Faner validation oracle; never a production front-position equation."""
    _validate_faner_inputs(
        total_loading, equilibrium_loading, critical_loading, particle_radius_m
    )
    if total_loading == equilibrium_loading:
        return 0.0
    if total_loading == critical_loading:
        return particle_radius_m
    return particle_radius_m * (
        (total_loading - equilibrium_loading)
        / (critical_loading - equilibrium_loading)
    ) ** (1.0 / 3.0)


def faner_loading_oracle(
    front_radius_m: float,
    equilibrium_loading: float,
    critical_loading: float,
    particle_radius_m: float,
) -> float:
    """Inverse Faner validation oracle, including exact dry/wet endpoints."""
    _validate_faner_inputs(
        equilibrium_loading, equilibrium_loading, critical_loading, particle_radius_m
    )
    if not math.isfinite(front_radius_m) or not 0.0 <= front_radius_m <= particle_radius_m:
        raise ValueError("Faner oracle radius must lie in [0, R]")
    return equilibrium_loading + (critical_loading - equilibrium_loading) * (
        front_radius_m / particle_radius_m
    ) ** 3


def _validate_faner_inputs(
    total_loading: float,
    equilibrium_loading: float,
    critical_loading: float,
    particle_radius_m: float,
) -> None:
    values: Sequence[float] = (
        total_loading,
        equilibrium_loading,
        critical_loading,
        particle_radius_m,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Faner oracle inputs must be finite")
    if particle_radius_m <= 0.0 or critical_loading <= equilibrium_loading:
        raise ValueError("Faner oracle has invalid radius or loading interval")
    if not equilibrium_loading <= total_loading <= critical_loading:
        raise ValueError("Faner oracle loading lies outside [Xe, Xc]")
