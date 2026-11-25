"""Native-law caloric closure for the engineering cell kernels (CELL-02b).

Frozen-law caloric quantities on the native component datums with the zero
additive gauge (``QSC10-NATIVE-PROPERTY-PLUS-ZERO-GAUGE-V1``): Span-Wagner
n-hexane, IAPWS-95 water, and the PHY-053 binary pressure-virial mixture
closure, exactly as ruled in PHY-057. This module mints no property number:
every returned value is an evaluation of a frozen law, and the only declared
content is the datum identity string itself, which a caller must state
explicitly rather than receive by default.

SCOPE BOUNDARY. The closure is defined for the water/hexane binary only. The
positive-carrier topology's inert third species has no frozen property law,
so the Law-1 kernel refuses ``NativeGasProperties`` at construction instead
of receiving invented carrier caloric constants; whether a carrier caloric
law should ever exist is an owner (Class-B) decision, not a default.

MODEL FORM (PHY-053 pressure-explicit virial, truncated after B, molar SI):

* partial molar residual enthalpy ``hbar_res_i = p (A_i - T dA_i/dT)`` with
  ``A_i = 2 sum_j y_j B_ij - B_mix`` (the same ``A_i`` as the fugacity path);
* mixture enthalpy ``h = y_w h0_w + y_h h0_h + p (B_mix - T dB_mix/dT)``;
* mixture heat capacity ``cp = y_w cp0_w + y_h cp0_h - p T d2B_mix/dT2``;
* saturation latent heats and their analytic slopes from the frozen
  Helmholtz saturation maps (PHY-057 selection 1).

Every tangent returned here is analytic (PHY-057 selection 2): the
temperature tangent of ``cp`` consumes the CELL-02a2 third-order virial
ladder and the ideal-gas cp slopes, and the latent slopes are the analytic
Clapeyron forms. The water ideal-gas caloric pieces reuse the established
residual-removal convention of ``binary_gas.water_ideal_enthalpy_molar``:
the ideal part is density-independent, so evaluating one finite-density
state and removing the analytic residual leaves the exact ideal value.
All tangents are caged by Richardson finite-difference gates in
``tests/test_core2_cell_native_caloric.py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import ClassVar

from .props import binary_gas
from .props import hexane as hx
from .props import water as wa

#: The one datum identity this closure is defined on (PHY-057 selection 4).
NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID = "QSC10-NATIVE-PROPERTY-PLUS-ZERO-GAUGE-V1"

#: Provenance marker: values here are frozen-law evaluations, not declared
#: engineering assumptions and not fitted parameters.
FROZEN_PROPERTY_LAW_EVALUATION = "FROZEN_PROPERTY_LAW_EVALUATION"

#: Probe density for the water ideal-part extraction; any positive value gives
#: the same ideal result because the analytic residual is removed exactly.
_WATER_IDEAL_PROBE_RHO_KG_M3 = 1.0


class NativeCaloricConfigurationError(ValueError):
    """A native caloric request contradicts the closure's declared contract."""


@dataclass(frozen=True, slots=True, kw_only=True)
class NativeGasProperties:
    """Native-law caloric mode marker for the engineering cell kernels.

    Carries no numbers: molar masses come from the frozen laws themselves and
    every caloric quantity is evaluated on demand. The datum identity is the
    single field so that callers state the native zero-gauge datum explicitly
    and datum-coherence gates compare real strings, never defaults.
    """

    energy_datum_id: str

    provenance_class: ClassVar[str] = FROZEN_PROPERTY_LAW_EVALUATION
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.energy_datum_id) is not str:
            raise NativeCaloricConfigurationError("energy_datum_id must be an exact str")
        if self.energy_datum_id != NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID:
            raise NativeCaloricConfigurationError(
                "the native caloric closure is defined only on "
                f"{NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID!r}; got {self.energy_datum_id!r}"
            )

    @property
    def hexane_molar_mass_kg_mol(self) -> float:
        return hx.M

    @property
    def water_molar_mass_kg_mol(self) -> float:
        return wa.M


@dataclass(frozen=True, slots=True, kw_only=True)
class NativeEnthalpyTangents:
    """A molar enthalpy value and its analytic (T, p, y_hexane) tangents."""

    value_j_mol: float
    d_dT_j_mol_k: float
    d_dp_j_mol_pa: float
    d_dy_hexane_j_mol: float


@dataclass(frozen=True, slots=True, kw_only=True)
class NativeHeatCapacityTangents:
    """A molar cp value and its analytic (T, p, y_hexane) tangents."""

    value_j_mol_k: float
    d_dT_j_mol_k2: float
    d_dp_j_mol_k_pa: float
    d_dy_hexane_j_mol_k: float


@dataclass(frozen=True, slots=True, kw_only=True)
class NativeLatentHeatTangent:
    """A saturation latent heat and its analytic slope along the locus."""

    value_j_mol: float
    d_dT_j_mol_k: float


def _require_state(T: float, pressure_pa: float, y_hexane: float) -> None:
    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("native caloric temperature must be positive and finite")
    if not math.isfinite(pressure_pa) or pressure_pa <= 0.0:
        raise ValueError("native caloric pressure must be positive and finite")
    if not math.isfinite(y_hexane) or not 0.0 <= y_hexane <= 1.0:
        raise ValueError("native caloric hexane mole fraction must lie in [0, 1]")


def _virial_orders(T: float) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    """Return (ww, wh, hh) second-virial ladders at orders 0..3 in d/dT."""

    cross = binary_gas.cross_second_virial(T)
    ww = (
        wa.second_virial(T),
        wa.second_virial_dT(T),
        wa.second_virial_dT2(T),
        wa.second_virial_dT3(T),
    )
    wh = (cross.B_wh, cross.dB_wh_dT, cross.d2B_wh_dT2, cross.d3B_wh_dT3)
    hh = (
        hx.second_virial(T),
        hx.second_virial_dT(T),
        hx.second_virial_dT2(T),
        hx.second_virial_dT3(T),
    )
    return ww, wh, hh


def _water_ideal_cp_and_slope_molar(T: float) -> tuple[float, float]:
    """Exact IAPWS-95 ideal-gas cp0 (J/(mol K)) and dcp0/dT (J/(mol K^2)).

    ``cv0 = cv(T, rho) + R tau^2 phir_tt`` removes the analytic residual from
    one finite-density state (the residual-removal convention), and
    ``dcp0/dT = -2 cv0 / T + R tau^3 phi0_ttt / T`` is the CELL-02a2 ideal
    slope with ``phi0_tt`` eliminated through ``cv0 = -R tau^2 phi0_tt``.
    """

    tau = wa.TC / T
    pure = wa.state(T, _WATER_IDEAL_PROBE_RHO_KG_M3)
    res = wa.residual(_WATER_IDEAL_PROBE_RHO_KG_M3 / wa.RHOC, tau)
    cv0_mass = pure.cv_mass + wa.R * tau * tau * res.phir_tt
    cp0_molar = (cv0_mass + wa.R) * wa.M
    dcp0_molar = (-2.0 * cv0_mass / T + wa.R * tau**3 * wa.phi0_ttt(tau) / T) * wa.M
    return cp0_molar, dcp0_molar


def mixture_enthalpy(T: float, pressure_pa: float, y_hexane: float) -> NativeEnthalpyTangents:
    """Molar mixture enthalpy on the native zero-gauge datum, with tangents.

    The temperature tangent is exactly the mixture cp; the composition
    tangent is the direct derivative, which equals ``hbar_h - hbar_w`` by
    Gibbs-Duhem (the identity is a test gate, not an assumption here).
    """

    _require_state(T, pressure_pa, y_hexane)
    y_h = y_hexane
    y_w = 1.0 - y_h
    ww, wh, hh = _virial_orders(T)
    b_mix = tuple(y_w * y_w * ww[k] + 2.0 * y_w * y_h * wh[k] + y_h * y_h * hh[k] for k in range(3))
    db_mix_dy = tuple(2.0 * (-y_w * ww[k] + (y_w - y_h) * wh[k] + y_h * hh[k]) for k in range(2))
    h0_w = binary_gas.water_ideal_enthalpy_molar(T)
    h0_h = hx.h_ideal(T)
    cp0_w, _ = _water_ideal_cp_and_slope_molar(T)
    cp0_h = hx.cp0(T)
    return NativeEnthalpyTangents(
        value_j_mol=y_w * h0_w + y_h * h0_h + pressure_pa * (b_mix[0] - T * b_mix[1]),
        d_dT_j_mol_k=y_w * cp0_w + y_h * cp0_h - pressure_pa * T * b_mix[2],
        d_dp_j_mol_pa=b_mix[0] - T * b_mix[1],
        d_dy_hexane_j_mol=(h0_h - h0_w) + pressure_pa * (db_mix_dy[0] - T * db_mix_dy[1]),
    )


def _partial_enthalpy(
    T: float, pressure_pa: float, y_hexane: float, *, species: str
) -> NativeEnthalpyTangents:
    _require_state(T, pressure_pa, y_hexane)
    y_h = y_hexane
    y_w = 1.0 - y_h
    ww, wh, hh = _virial_orders(T)
    b_mix = tuple(y_w * y_w * ww[k] + 2.0 * y_w * y_h * wh[k] + y_h * y_h * hh[k] for k in range(3))
    db_mix_dy = tuple(2.0 * (-y_w * ww[k] + (y_w - y_h) * wh[k] + y_h * hh[k]) for k in range(2))
    if species == "hexane":
        a = tuple(2.0 * (y_w * wh[k] + y_h * hh[k]) - b_mix[k] for k in range(3))
        da_dy = tuple(2.0 * (hh[k] - wh[k]) - db_mix_dy[k] for k in range(2))
        h0 = hx.h_ideal(T)
        cp0 = hx.cp0(T)
    else:
        a = tuple(2.0 * (y_w * ww[k] + y_h * wh[k]) - b_mix[k] for k in range(3))
        da_dy = tuple(2.0 * (wh[k] - ww[k]) - db_mix_dy[k] for k in range(2))
        h0 = binary_gas.water_ideal_enthalpy_molar(T)
        cp0, _ = _water_ideal_cp_and_slope_molar(T)
    return NativeEnthalpyTangents(
        value_j_mol=h0 + pressure_pa * (a[0] - T * a[1]),
        d_dT_j_mol_k=cp0 - pressure_pa * T * a[2],
        d_dp_j_mol_pa=a[0] - T * a[1],
        d_dy_hexane_j_mol=pressure_pa * (da_dy[0] - T * da_dy[1]),
    )


def hexane_partial_enthalpy(
    T: float, pressure_pa: float, y_hexane: float
) -> NativeEnthalpyTangents:
    """Partial molar hexane enthalpy in the binary gas, with tangents."""

    return _partial_enthalpy(T, pressure_pa, y_hexane, species="hexane")


def water_partial_enthalpy(T: float, pressure_pa: float, y_hexane: float) -> NativeEnthalpyTangents:
    """Partial molar water enthalpy in the binary gas, with tangents."""

    return _partial_enthalpy(T, pressure_pa, y_hexane, species="water")


def mixture_heat_capacity(
    T: float, pressure_pa: float, y_hexane: float
) -> NativeHeatCapacityTangents:
    """Molar mixture cp with analytic tangents (CELL-02a2 third ladder).

    ``d_dT`` is the Ackermann-path derivative the PHY-057 design-basis
    amendment exists for: it consumes ``d3B_mix/dT3`` and both ideal-gas cp
    slopes, so a dropped term cannot hide behind a finite-difference cage.
    """

    _require_state(T, pressure_pa, y_hexane)
    y_h = y_hexane
    y_w = 1.0 - y_h
    ww, wh, hh = _virial_orders(T)
    b_mix = tuple(y_w * y_w * ww[k] + 2.0 * y_w * y_h * wh[k] + y_h * y_h * hh[k] for k in range(4))
    db_mix_dy2 = 2.0 * (-y_w * ww[2] + (y_w - y_h) * wh[2] + y_h * hh[2])
    cp0_w, dcp0_w = _water_ideal_cp_and_slope_molar(T)
    cp0_h = hx.cp0(T)
    dcp0_h = hx.cp0_dT(T)
    return NativeHeatCapacityTangents(
        value_j_mol_k=y_w * cp0_w + y_h * cp0_h - pressure_pa * T * b_mix[2],
        d_dT_j_mol_k2=y_w * dcp0_w + y_h * dcp0_h - pressure_pa * (b_mix[2] + T * b_mix[3]),
        d_dp_j_mol_k_pa=-T * b_mix[2],
        d_dy_hexane_j_mol_k=(cp0_h - cp0_w) - pressure_pa * T * db_mix_dy2,
    )


#: CELL-02e-pd inverse-caloric domain: inside both frozen laws' validity and
#: covering the Law-2 cell domain (330-435 K) with margin on both sides.
#: Refuse-not-clamp: an iterate leaving this interval is a typed refusal.
INVERSE_CALORIC_TEMPERATURE_DOMAIN_K = (250.0, 500.0)
INVERSE_CALORIC_MAX_ITERATIONS = 64
INVERSE_CALORIC_TOLERANCE_J_MOL = 1.0e-9


def temperature_from_mixture_enthalpy(
    h_j_mol: float,
    pressure_pa: float,
    y_hexane: float,
    *,
    seed_temperature_k: float = 360.0,
) -> float:
    """Invert ``h_mix(T, p, y) = h`` for T by bracketed Newton on the laws.

    Well-posedness: ``dh_mix/dT = cp_mix > 0`` throughout the declared
    domain (measured 35-194 J/(mol K) across compositions), so the map is
    strictly monotone and the root, when the target enthalpy lies inside
    the domain's enthalpy range, exists and is unique. The solver is the
    standard safeguarded (bracketed) Newton: the domain endpoints form a
    guaranteed bracket, the Newton derivative is the analytic,
    Richardson-caged mixture cp - no finite differencing and no fitted
    structure - and any Newton trial that leaves the current bracket is
    replaced by its bisection midpoint, which preserves guaranteed
    convergence on a monotone function while keeping Newton's terminal
    rate. Acceptance is the NaN-safe residual criterion at
    ``INVERSE_CALORIC_TOLERANCE_J_MOL`` (about 1e-11 K in temperature).
    A target outside the domain's enthalpy range, a non-positive cp, or
    budget exhaustion are typed refusals, never clamps.
    """

    if not math.isfinite(h_j_mol):
        raise ValueError("inverse caloric enthalpy must be finite")
    low, high = INVERSE_CALORIC_TEMPERATURE_DOMAIN_K
    if not low <= seed_temperature_k <= high:
        raise ValueError("inverse caloric seed temperature left the declared domain")
    low_residual = mixture_enthalpy(low, pressure_pa, y_hexane).value_j_mol - h_j_mol
    high_residual = mixture_enthalpy(high, pressure_pa, y_hexane).value_j_mol - h_j_mol
    if not low_residual <= INVERSE_CALORIC_TOLERANCE_J_MOL or not (
        high_residual >= -INVERSE_CALORIC_TOLERANCE_J_MOL
    ):
        raise ValueError(
            "inverse caloric target enthalpy is outside the declared domain's "
            "enthalpy range; refused without clamping"
        )
    bracket_low, bracket_high = low, high
    temperature = seed_temperature_k
    for _ in range(INVERSE_CALORIC_MAX_ITERATIONS):
        state = mixture_enthalpy(temperature, pressure_pa, y_hexane)
        residual = state.value_j_mol - h_j_mol
        if abs(residual) <= INVERSE_CALORIC_TOLERANCE_J_MOL:
            return temperature
        if residual > 0.0:
            bracket_high = temperature
        else:
            bracket_low = temperature
        heat_capacity = state.d_dT_j_mol_k
        if not heat_capacity > 0.0:
            raise ValueError(
                "inverse caloric monotonicity guarantee failed: cp_mix is not positive"
            )
        trial = temperature - residual / heat_capacity
        if not bracket_low < trial < bracket_high:
            trial = 0.5 * (bracket_low + bracket_high)
        temperature = trial
    raise ValueError("inverse caloric solve did not converge within the frozen budget")


def hexane_latent_heat(T: float) -> NativeLatentHeatTangent:
    """Span-Wagner saturation latent heat (J/mol) and its analytic slope."""

    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("native caloric temperature must be positive and finite")
    sat = hx.saturation(T)
    if not sat.converged:
        raise ValueError(f"hexane saturation did not converge at T={T}")
    return NativeLatentHeatTangent(value_j_mol=sat.dh_vap, d_dT_j_mol_k=hx.dh_vap_dT(T))


def water_latent_heat(T: float) -> NativeLatentHeatTangent:
    """IAPWS-95 saturation latent heat (J/mol) and its analytic slope."""

    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("native caloric temperature must be positive and finite")
    sat = wa.saturation(T)
    if not sat.converged:
        raise ValueError(f"water saturation did not converge at T={T}")
    return NativeLatentHeatTangent(
        value_j_mol=sat.dh_vap * wa.M, d_dT_j_mol_k=wa.dh_vap_dT(T) * wa.M
    )


__all__ = [
    "INVERSE_CALORIC_MAX_ITERATIONS",
    "INVERSE_CALORIC_TEMPERATURE_DOMAIN_K",
    "INVERSE_CALORIC_TOLERANCE_J_MOL",
    "FROZEN_PROPERTY_LAW_EVALUATION",
    "NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID",
    "NativeCaloricConfigurationError",
    "NativeEnthalpyTangents",
    "NativeGasProperties",
    "NativeHeatCapacityTangents",
    "NativeLatentHeatTangent",
    "hexane_latent_heat",
    "hexane_partial_enthalpy",
    "mixture_enthalpy",
    "mixture_heat_capacity",
    "temperature_from_mixture_enthalpy",
    "water_latent_heat",
    "water_partial_enthalpy",
]
