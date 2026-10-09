r"""Dry-side water/n-hexane thermodynamic closure for the coupled sphere.

At fixed temperature and nominal pore pressure, the only gas components are
water and n-hexane.  PHY-053 supplies one pressure-second-virial state for the
EOS, fugacity, residual calorics, and transported partial enthalpy.  The
n-hexane activity and complete dry-side storage are

    a_h = f_h,g(T, P, y) / f_h,l^pure(T, P)
    C_h,d = eps_g c_h,g + rho_dm,p W_h(a_h, T, w_o).

The internal-energy density uses the same component data exactly once:

    u_d = rho_dm,p [h_dm(T,w_o) + X_w h_w,l + U_h,ret]
          + eps_g n_g u_g,mix,

where ``U_h,ret`` already contains the retained-state binding potential and
partial-volume convention.  Consequently callers must not add independent
latent-heat or sorption-heat sources.  ``h_dm`` is PHY-048's
reference-anchored native-composite/oil construction, not a native value plus
a duplicate full oil sensible term.

This module supersedes the ideal ``c/c_sat`` activity used only by gate 1c's
isothermal Cardarelli oracle.  Its fixed ``X_water`` state remains the
explicitly nonqualifying legacy-trajectory ablation; it is not the G1G-01
production retained/pore-water equilibrium.  Production variable total-water
storage and transport live in :mod:`coupled_pore`.  This module deliberately
contains no transport operator.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa

GAB_DIRECT_T_MIN_K = 323.15
GAB_DIRECT_T_MAX_K = 368.15


class DryThermoTopologyError(ValueError):
    """A state is incompatible with the mobile-liquid-free dry branch."""


@dataclass(frozen=True)
class DryThermoParams:
    """Fixed dry-region material and nominal-pressure parameters."""

    pressure_pa: float = 101_325.0
    eps_g: float = 0.141
    rho_dm_p: float = 1159.65
    cp_dry_meal: float = sp.CP_NATIVE_DRY_MEAL
    cp_oil: float = sp.CP_OIL
    T_ref_solid: float = 298.15
    X_water: float = 0.10
    w_o: float = sp.W_O_REF
    gab: sp.GabParams = sp.GabParams()
    oil: sp.OilIsotherm = sp.OilIsotherm()
    k_wh: float = bg.K_WH_CENTRAL


@dataclass(frozen=True)
class DryThermoState:
    """One complete fixed-(T,P,y) dry-side storage and caloric state."""

    temperature_k: float
    pressure_pa: float
    y_water: float
    y_hexane: float
    binary_gas: bg.BinaryVirialState
    hexane_gas_density_kg_m3: float
    gas_hexane_fugacity_pa: float
    pure_liquid_hexane_fugacity_pa: float
    hexane_activity: float
    gas_water_fugacity_pa: float
    pure_liquid_water_fugacity_pa: float
    water_activity: float
    retained_hexane_loading: float
    total_hexane_concentration_kg_m3: float
    dry_meal_heat_capacity_density_j_m3_k: float
    dry_meal_sensible_energy_density_j_m3: float
    retained_water_energy_density_j_m3: float
    retained_hexane_energy_density_j_m3: float
    pore_gas_energy_density_j_m3: float
    energy_density_j_m3: float
    pore_gas_internal_energy_molar_j_mol: float
    transported_hexane_partial_enthalpy_j_kg: float
    outside_project_pressure: bool
    below_cross_direct_evidence: bool
    outside_gab_direct_temperature_evidence: bool

    @property
    def has_source_domain_extrapolation(self) -> bool:
        """Whether any pressure/cross/GAB applicability flag is active."""
        return (
            self.outside_project_pressure
            or self.below_cross_direct_evidence
            or self.outside_gab_direct_temperature_evidence
        )


def hexane_gas_density_from_mole_fraction(
    T: float, y_hexane: float, p: DryThermoParams = DryThermoParams()
) -> float:
    """Map ``y_h`` to physical gas-phase hexane density using PHY-053."""
    _validate_params(p)
    _validate_temperature(T)
    y_h = _validate_mole_fraction(y_hexane)
    return bg.state(T, p.pressure_pa, 1.0 - y_h, y_h, k_wh=p.k_wh).hexane_density_kg_m3


def hexane_mole_fraction_from_gas_density(
    T: float, hexane_gas_density_kg_m3: float, p: DryThermoParams = DryThermoParams()
) -> float:
    """Invert physical ``c_h,g`` to ``y_h`` by a bracketed EOS solve.

    The inversion covers the binary-gas composition interval itself.  The
    subsequent dry-state evaluation separately rejects a supersaturated state;
    neither operation clips its input.
    """
    _validate_params(p)
    _validate_temperature(T)
    c_h = hexane_gas_density_kg_m3
    if not math.isfinite(c_h) or c_h < 0.0:
        raise ValueError("hexane gas density must be finite and non-negative")
    if c_h == 0.0:
        return 0.0
    c_max = hexane_gas_density_from_mole_fraction(T, 1.0, p)
    if c_h > c_max:
        raise ValueError(
            f"hexane gas density {c_h} exceeds the fixed-(T,P) binary maximum "
            f"{c_max}; reject, do not clamp"
        )
    if c_h == c_max:
        return 1.0

    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        c_mid = hexane_gas_density_from_mole_fraction(T, mid, p)
        if c_mid < c_h:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 2.0e-15:
            break
    y_h = 0.5 * (lo + hi)
    residual = hexane_gas_density_from_mole_fraction(T, y_h, p) - c_h
    if abs(residual) > 5.0e-13 * max(c_h, c_max, 1.0):
        raise ValueError("hexane gas-density EOS inversion did not converge")
    return y_h


def evaluate(
    T: float, y_hexane: float, p: DryThermoParams = DryThermoParams()
) -> DryThermoState:
    """Evaluate complete dry storage, energy, and transported hexane enthalpy."""
    _validate_params(p)
    _validate_temperature(T)
    y_h = _validate_mole_fraction(y_hexane)
    y_w = 1.0 - y_h
    mixture = bg.state(T, p.pressure_pa, y_w, y_h, k_wh=p.k_wh)

    liquid_hexane = hx.state_Tp(T, p.pressure_pa, "liquid")
    liquid_fugacity = liquid_hexane.fugacity
    if not math.isfinite(liquid_fugacity) or liquid_fugacity <= 0.0:
        raise ValueError("pure-liquid hexane fugacity must be positive and finite")
    gas_fugacity = mixture.fugacity_hexane_pa
    activity = 0.0 if y_h == 0.0 else gas_fugacity / liquid_fugacity
    if not math.isfinite(activity) or activity < 0.0:
        raise ValueError("hexane fugacity activity must be finite and non-negative")
    # This is a branch condition, not a numerical clamp.  Root-level roundoff
    # at an a_h=1 interface is admitted without changing the evaluated value.
    if activity > 1.0 and math.log(activity) > 1.0e-10:
        raise DryThermoTopologyError(
            f"dry-side hexane activity {activity} exceeds liquid equilibrium; "
            "activate the mobile-liquid branch"
        )

    retained = sp.retained_hexane(activity, T, p.w_o, p.gab, p.oil)
    c_h = mixture.hexane_density_kg_m3
    storage = p.eps_g * c_h + p.rho_dm_p * retained

    dry_meal_cp = sp.composite_dry_meal_heat_capacity(p.w_o, p.cp_dry_meal, p.cp_oil)
    dry_meal_capacity = p.rho_dm_p * dry_meal_cp
    dry_meal_energy = p.rho_dm_p * sp.composite_dry_meal_sensible_energy(
        T, p.T_ref_solid, p.w_o, p.cp_dry_meal, p.cp_oil
    )
    liquid_water = wa.state_Tp(T, p.pressure_pa, "liquid")
    # G1G-01 nominal retained-water bookkeeping assigns zero retained partial
    # volume.  The selected liquid-like retained enthalpy therefore enters the
    # mixture internal-energy ledger directly; subtracting P/rho_l would be
    # the separately named pure-liquid-partial-volume sensitivity.
    water_retained_u = liquid_water.h_mass
    water_liquid_fugacity = liquid_water.fugacity
    gas_water_fugacity = mixture.fugacity_water_pa
    water_activity = gas_water_fugacity / water_liquid_fugacity
    retained_water_energy = p.rho_dm_p * p.X_water * water_retained_u
    retained_hexane_energy = p.rho_dm_p * sp.retained_internal_energy(
        retained,
        T,
        p.w_o,
        p.gab,
        p.oil,
        pressure_pa=p.pressure_pa,
        activity=activity,
    )

    h_ideal_mix = (
        y_w * _water_ideal_enthalpy_molar(T) + y_h * hx.h_ideal(T)
    )
    h_mix = h_ideal_mix + mixture.residual_enthalpy_molar
    u_mix = h_mix - p.pressure_pa / mixture.molar_density_mol_m3
    pore_gas_energy = p.eps_g * mixture.molar_density_mol_m3 * u_mix
    transported_h = (
        hx.h_ideal(T) + mixture.partial_residual_enthalpy_hexane_molar
    ) / hx.M
    energy = (
        dry_meal_energy
        + retained_water_energy
        + retained_hexane_energy
        + pore_gas_energy
    )
    values = (
        gas_fugacity,
        activity,
        gas_water_fugacity,
        water_liquid_fugacity,
        water_activity,
        retained,
        storage,
        dry_meal_capacity,
        dry_meal_energy,
        retained_water_energy,
        retained_hexane_energy,
        u_mix,
        pore_gas_energy,
        transported_h,
        energy,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("dry thermodynamic state produced a non-finite value")

    return DryThermoState(
        temperature_k=T,
        pressure_pa=p.pressure_pa,
        y_water=y_w,
        y_hexane=y_h,
        binary_gas=mixture,
        hexane_gas_density_kg_m3=c_h,
        gas_hexane_fugacity_pa=gas_fugacity,
        pure_liquid_hexane_fugacity_pa=liquid_fugacity,
        hexane_activity=activity,
        gas_water_fugacity_pa=gas_water_fugacity,
        pure_liquid_water_fugacity_pa=water_liquid_fugacity,
        water_activity=water_activity,
        retained_hexane_loading=retained,
        total_hexane_concentration_kg_m3=storage,
        dry_meal_heat_capacity_density_j_m3_k=dry_meal_capacity,
        dry_meal_sensible_energy_density_j_m3=dry_meal_energy,
        retained_water_energy_density_j_m3=retained_water_energy,
        retained_hexane_energy_density_j_m3=retained_hexane_energy,
        pore_gas_energy_density_j_m3=pore_gas_energy,
        energy_density_j_m3=energy,
        pore_gas_internal_energy_molar_j_mol=u_mix,
        transported_hexane_partial_enthalpy_j_kg=transported_h,
        outside_project_pressure=mixture.outside_project_pressure,
        below_cross_direct_evidence=mixture.below_cross_direct_evidence,
        outside_gab_direct_temperature_evidence=not (
            GAB_DIRECT_T_MIN_K <= T <= GAB_DIRECT_T_MAX_K
        ),
    )


def evaluate_from_gas_density(
    T: float,
    hexane_gas_density_kg_m3: float,
    p: DryThermoParams = DryThermoParams(),
) -> DryThermoState:
    """Evaluate the same state with physical ``c_h,g`` as coordinate."""
    y_h = hexane_mole_fraction_from_gas_density(T, hexane_gas_density_kg_m3, p)
    return evaluate(T, y_h, p)


def _water_ideal_enthalpy_molar(T: float) -> float:
    """IAPWS-95 ideal-water enthalpy on its component datum, J/mol.

    The public water authority is mass-based.  Subtracting its analytic
    zero-density residual enthalpy from any finite-density state recovers the
    density-independent ideal contribution without introducing a second fit.
    """
    return bg.water_ideal_enthalpy_molar(T)


def _validate_mole_fraction(y_hexane: float) -> float:
    if not math.isfinite(y_hexane) or not 0.0 <= y_hexane <= 1.0:
        raise ValueError("hexane gas mole fraction must be finite and lie in [0, 1]")
    return y_hexane


def _validate_temperature(T: float) -> None:
    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("dry thermodynamic temperature must be positive and finite")
    if T >= hx.TC:
        raise ValueError(
            "dry hexane activity requires a subcritical pure-liquid reference state"
        )


def _validate_params(p: DryThermoParams) -> None:
    scalars = (
        p.pressure_pa,
        p.eps_g,
        p.rho_dm_p,
        p.cp_dry_meal,
        p.cp_oil,
        p.T_ref_solid,
        p.X_water,
        p.w_o,
        p.k_wh,
        p.gab.Xm,
        p.gab.C0,
        p.gab.dHC_R,
        p.gab.K0,
        p.gab.dHK_R,
        p.oil.A0,
        p.oil.B,
    )
    if not all(math.isfinite(value) for value in scalars):
        raise ValueError("dry thermodynamic parameters must be finite")
    if p.pressure_pa <= 0.0:
        raise ValueError("dry nominal pressure must be positive")
    if not 0.0 < p.eps_g < 1.0:
        raise ValueError("dry gas holdup must lie strictly between zero and one")
    if p.rho_dm_p <= 0.0 or p.cp_dry_meal <= 0.0 or p.cp_oil <= 0.0:
        raise ValueError("dry-meal density, native-meal cp, and oil cp must be positive")
    water_bounds = sp.LuikovParams()
    if not water_bounds.W_ref <= p.X_water <= water_bounds.W_cap:
        raise ValueError(
            "retained-water loading must lie on the qualified positive-moisture "
            f"Luikov branch [{water_bounds.W_ref}, {water_bounds.W_cap}] "
            "inside the evidence cap"
        )
    if not sp.W_O_MIN <= p.w_o <= sp.W_O_MAX:
        raise ValueError("residual-oil fraction is outside the PHY-048 envelope")
    if min(p.gab.Xm, p.gab.C0, p.gab.K0, p.oil.A0, p.oil.B) <= 0.0:
        raise ValueError("dry sorption parameters must be positive")
    if p.gab.dHC_R < 0.0 or p.gab.dHK_R < 0.0:
        raise ValueError("dry sorption temperature scales must be non-negative")
    # Validate the frozen cross-correlation interval even for a pure endpoint.
    if not bg.K_WH_MIN <= p.k_wh <= bg.K_WH_MAX:
        raise ValueError("k_wh is outside the frozen qualification interval")
