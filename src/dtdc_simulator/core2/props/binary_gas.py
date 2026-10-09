r"""PHY-053 water/n-hexane pressure-second-virial mixture authority.

The frozen low-pressure closure is

    B_m = sum_i sum_j y_i y_j B_ij
    ln(phi_i) = P/(R T) [2 sum_j y_j B_ij - B_m].

Pure virials and their analytic temperature derivatives come from the frozen
IAPWS-95 and Span--Wagner authorities.  The water/hexane cross term is the
Plyasunov--Shock implementation of the Tsonopoulos correlation with the frozen
central interaction ``k_wh=0.494``.  The same derivatives supply the residual
mixture and partial molar enthalpies, so fugacity and front energy transport do
not use inconsistent thermodynamics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from dtdc_simulator.core2.exact_cache import float_bits_key
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import water as wa

R = hx.R
K_WH_CENTRAL = 0.494
K_WH_MIN = 0.477
K_WH_MAX = 0.511
DIRECT_EVIDENCE_MIN_T = 363.2
PROJECT_PRESSURE_MIN_PA = 101_000.0
PROJECT_PRESSURE_MAX_PA = 170_000.0

# Source-consistent constants frozen under PHY-053.  These belong to the cross
# correlation and intentionally retain its 507.6 K n-hexane critical constant;
# pure properties continue to use the PHY-049 authority's own constants.
TC_W = 647.096
TC_H = 507.6
VC_W = 55.95e-6  # m3/mol
VC_H = 371.0e-6  # m3/mol
ZC_W = 0.2294
ZC_H = 0.2660
OMEGA_W = 0.3443
OMEGA_H = 0.3010


@dataclass(frozen=True)
class CrossVirial:
    """Water/n-hexane cross second virial and analytic derivative."""

    T: float
    k_wh: float
    Tc12: float
    B_wh: float
    dB_wh_dT: float
    below_direct_evidence: bool
    #: CELL-02a (PHY-057): analytic second temperature derivative.
    d2B_wh_dT2: float = 0.0
    #: CELL-02a2: analytic third temperature derivative.
    d3B_wh_dT3: float = 0.0


@dataclass(frozen=True)
class BinaryVirialState:
    """One fixed-(T,P,y) binary-gas fugacity/caloric state in molar SI."""

    T: float
    pressure_pa: float
    y_water: float
    y_hexane: float
    B_ww: float
    B_wh: float
    B_hh: float
    dB_ww_dT: float
    dB_wh_dT: float
    dB_hh_dT: float
    B_mix: float
    dB_mix_dT: float
    compressibility: float
    molar_density_mol_m3: float
    ln_phi_water: float
    ln_phi_hexane: float
    residual_enthalpy_molar: float
    partial_residual_enthalpy_water_molar: float
    partial_residual_enthalpy_hexane_molar: float
    virial_correction_magnitude: float
    outside_project_pressure: bool
    below_cross_direct_evidence: bool

    @property
    def phi_water(self) -> float:
        return math.exp(self.ln_phi_water)

    @property
    def phi_hexane(self) -> float:
        return math.exp(self.ln_phi_hexane)

    @property
    def fugacity_water_pa(self) -> float:
        return self.y_water * self.phi_water * self.pressure_pa

    @property
    def fugacity_hexane_pa(self) -> float:
        return self.y_hexane * self.phi_hexane * self.pressure_pa

    @property
    def water_density_kg_m3(self) -> float:
        return self.y_water * self.molar_density_mol_m3 * wa.M

    @property
    def hexane_density_kg_m3(self) -> float:
        return self.y_hexane * self.molar_density_mol_m3 * hx.M


def cross_second_virial(T: float, k_wh: float = K_WH_CENTRAL) -> CrossVirial:
    """Return B_wh and analytic dB_wh/dT from frozen combining rules."""
    if T <= 0.0 or not math.isfinite(T):
        raise ValueError("binary-virial temperature must be positive and finite")
    if not K_WH_MIN <= k_wh <= K_WH_MAX:
        raise ValueError(
            f"k_wh={k_wh} outside frozen qualification interval [{K_WH_MIN}, {K_WH_MAX}]"
        )
    Tc12 = math.sqrt(TC_W * TC_H) * (1.0 - k_wh)
    Vc12 = ((VC_W ** (1.0 / 3.0) + VC_H ** (1.0 / 3.0)) / 2.0) ** 3
    Zc12 = 0.5 * (ZC_W + ZC_H)
    omega12 = 0.5 * (OMEGA_W + OMEGA_H)
    Pc12 = Zc12 * R * Tc12 / Vc12
    Tr = T / Tc12
    f0 = 0.1445 - 0.330 / Tr - 0.1385 / Tr**2 - 0.0121 / Tr**3 - 0.000607 / Tr**8
    f1 = 0.0637 + 0.331 / Tr**2 - 0.423 / Tr**3 - 0.008 / Tr**8
    df0_dT = (0.330 / Tr + 2.0 * 0.1385 / Tr**2 + 3.0 * 0.0121 / Tr**3 + 8.0 * 0.000607 / Tr**8) / T
    df1_dT = (-2.0 * 0.331 / Tr**2 + 3.0 * 0.423 / Tr**3 + 8.0 * 0.008 / Tr**8) / T
    scale = R * Tc12 / Pc12
    # CELL-02a: for a term a/Tr^n, d2/dT2 = n (n+1) a / (Tr^n T^2).
    d2f0_dT2 = (
        -(2.0 * 0.330 / Tr + 6.0 * 0.1385 / Tr**2 + 12.0 * 0.0121 / Tr**3 + 72.0 * 0.000607 / Tr**8)
        / T**2
    )
    d2f1_dT2 = (6.0 * 0.331 / Tr**2 - 12.0 * 0.423 / Tr**3 - 72.0 * 0.008 / Tr**8) / T**2
    return CrossVirial(
        T=T,
        k_wh=k_wh,
        Tc12=Tc12,
        B_wh=scale * (f0 + omega12 * f1),
        dB_wh_dT=scale * (df0_dT + omega12 * df1_dT),
        below_direct_evidence=T < DIRECT_EVIDENCE_MIN_T,
        d2B_wh_dT2=scale * (d2f0_dT2 + omega12 * d2f1_dT2),
        d3B_wh_dT3=scale
        * (
            (
                6.0 * 0.330 / Tr
                + 24.0 * 0.1385 / Tr**2
                + 60.0 * 0.0121 / Tr**3
                + 720.0 * 0.000607 / Tr**8
            )
            / T**3
            + omega12
            * (-24.0 * 0.331 / Tr**2 + 60.0 * 0.423 / Tr**3 + 720.0 * 0.008 / Tr**8)
            / T**3
        ),
    )


def state(
    T: float,
    pressure_pa: float,
    y_water: float,
    y_hexane: float,
    *,
    k_wh: float = K_WH_CENTRAL,
) -> BinaryVirialState:
    """Evaluate the frozen closure, memoized for bit-exact repeated inputs."""

    return _state_cached(
        T,
        pressure_pa,
        y_water,
        y_hexane,
        k_wh,
        float_bits_key(T, pressure_pa, y_water, y_hexane, k_wh),
    )


@lru_cache(maxsize=4096, typed=True)
def _state_cached(
    T: float,
    pressure_pa: float,
    y_water: float,
    y_hexane: float,
    k_wh: float,
    _scalar_key: tuple[bytes, ...],
) -> BinaryVirialState:
    """Unmodified binary fugacity/caloric evaluation behind the exact cache."""

    if pressure_pa <= 0.0 or not math.isfinite(pressure_pa):
        raise ValueError("binary-gas pressure must be positive and finite")
    yw, yh = _mole_fractions(y_water, y_hexane)
    cross = cross_second_virial(T, k_wh)
    Bww, Bhh = wa.second_virial(T), hx.second_virial(T)
    dBww, dBhh = wa.second_virial_dT(T), hx.second_virial_dT(T)
    Bwh, dBwh = cross.B_wh, cross.dB_wh_dT
    Bmix = yw * yw * Bww + 2.0 * yw * yh * Bwh + yh * yh * Bhh
    dBmix = yw * yw * dBww + 2.0 * yw * yh * dBwh + yh * yh * dBhh
    A_water = 2.0 * (yw * Bww + yh * Bwh) - Bmix
    A_hexane = 2.0 * (yw * Bwh + yh * Bhh) - Bmix
    dA_water = 2.0 * (yw * dBww + yh * dBwh) - dBmix
    dA_hexane = 2.0 * (yw * dBwh + yh * dBhh) - dBmix
    prefactor = pressure_pa / (R * T)
    hbar_water = pressure_pa * (A_water - T * dA_water)
    hbar_hexane = pressure_pa * (A_hexane - T * dA_hexane)
    h_residual = pressure_pa * (Bmix - T * dBmix)
    compressibility = 1.0 + Bmix * pressure_pa / (R * T)
    if compressibility <= 0.0:
        raise ValueError("binary pressure-virial compressibility is non-positive")
    molar_density = pressure_pa / (compressibility * R * T)
    return BinaryVirialState(
        T=T,
        pressure_pa=pressure_pa,
        y_water=yw,
        y_hexane=yh,
        B_ww=Bww,
        B_wh=Bwh,
        B_hh=Bhh,
        dB_ww_dT=dBww,
        dB_wh_dT=dBwh,
        dB_hh_dT=dBhh,
        B_mix=Bmix,
        dB_mix_dT=dBmix,
        compressibility=compressibility,
        molar_density_mol_m3=molar_density,
        ln_phi_water=prefactor * A_water,
        ln_phi_hexane=prefactor * A_hexane,
        residual_enthalpy_molar=h_residual,
        partial_residual_enthalpy_water_molar=hbar_water,
        partial_residual_enthalpy_hexane_molar=hbar_hexane,
        virial_correction_magnitude=abs(Bmix * pressure_pa / (R * T)),
        outside_project_pressure=not (
            PROJECT_PRESSURE_MIN_PA <= pressure_pa <= PROJECT_PRESSURE_MAX_PA
        ),
        below_cross_direct_evidence=(cross.below_direct_evidence and yw > 0.0 and yh > 0.0),
    )


def clear_state_cache() -> None:
    """Drop every cached binary-gas ``(T,P,y,k_wh)`` state."""

    _state_cached.cache_clear()


def state_cache_info():
    """Return bounded-cache hit/miss/size diagnostics."""

    return _state_cached.cache_info()


def hexane_partial_enthalpy_mass(
    T: float,
    pressure_pa: float,
    y_water: float,
    y_hexane: float,
    *,
    k_wh: float = K_WH_CENTRAL,
) -> float:
    """Dry-pore-gas partial h_h, J/kg hexane, on the PHY-049 datum."""
    mixture = state(T, pressure_pa, y_water, y_hexane, k_wh=k_wh)
    return (hx.h_ideal(T) + mixture.partial_residual_enthalpy_hexane_molar) / hx.M


def water_ideal_enthalpy_molar(T: float) -> float:
    """IAPWS-95 ideal-water enthalpy on its component datum, J/mol.

    IAPWS-95 is mass based.  Removing the analytic zero-density residual
    enthalpy from a finite-density state gives the density-independent ideal
    contribution without introducing another correlation or caloric datum.
    """
    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("water ideal-gas temperature must be positive and finite")
    rho_mass = 1.0
    pure = wa.state(T, rho_mass)
    delta, tau = rho_mass / wa.RHOC, wa.TC / T
    residual = wa.residual(delta, tau)
    residual_enthalpy_mass = wa.R * T * (tau * residual.phir_t + delta * residual.phir_d)
    return (pure.h_mass - residual_enthalpy_mass) * wa.M


def water_partial_enthalpy_mass(
    T: float,
    pressure_pa: float,
    y_water: float,
    y_hexane: float,
    *,
    k_wh: float = K_WH_CENTRAL,
) -> float:
    """Dry-pore-gas partial ``h_w``, J/kg water, on the IAPWS datum."""
    mixture = state(T, pressure_pa, y_water, y_hexane, k_wh=k_wh)
    return (water_ideal_enthalpy_molar(T) + mixture.partial_residual_enthalpy_water_molar) / wa.M


def _mole_fractions(y_water: float, y_hexane: float) -> tuple[float, float]:
    values = (y_water, y_hexane)
    if not all(math.isfinite(value) and value >= 0.0 for value in values):
        raise ValueError("binary-gas mole fractions must be finite and non-negative")
    total = y_water + y_hexane
    if abs(total - 1.0) > 1.0e-12:
        raise ValueError(f"binary-gas mole fractions must sum to one, got {total}")
    if total <= 0.0:
        raise ValueError("binary-gas composition is empty")
    return y_water / total, y_hexane / total
