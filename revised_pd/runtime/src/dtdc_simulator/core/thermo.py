"""Pure thermophysical / sorption / VLLE property functions — Coletto, Bandoni &
Blanco (2022) *J. Food Eng.* 318, 110870, Appendix A (sorption/critical-hexane,
eqs. 3-7 and A.31) and Appendix B (thermophysical correlations, B.1-B.12).

BuildSpec §3 invariant: this module must never import `config/`, `engine/`, or do
file/network I/O — it mirrors `core/model.py`'s purity constraint. Callers (in
practice `config/builder.py`, later `core/model.py`) translate validated cold
config into the plain dataclasses defined here.

This is BuildSpec §14 **M1**: thermo/sorption/properties only — pure, unit
-tested, no solver. It does not implement the zone sub-models (PHZ/FTRZ/DCZ,
§14 M2) or touch `core/model.py`.

PROVENANCE — read before trusting any number fed into these functions: the
soybean parameters in `properties/soybean.yaml` are tagged `[PAPER]`/`[PLACE]`/
`[STD]`/`[DERIVED]`; notably the oil-sorption power-law (`A0`, `B`) and the
heat-of-sorption constants (`sorption_C0`, `sorption_C1`) are cited by Coletto
(2022) to two unpublished PhD theses. Cardarelli (1998) is now held in part —
its §2.4 sorption-energy chapter and nomenclature pages, see
`docs/GT_PS2_CARDARELLI1998_INTAKE_RECORD_2026-09-01.md` — and those pages print
the measured net isosteric heat *data* (Fig 2.23), not a fitted C0/C1 pair, so
the pair remains Coletto's own fit. Faner (2008) is still not available to us.
Their *functional forms* below are implemented exactly as published; the soybean
C0/C1 values are a documented fit to Cardarelli & Crapiste's published caloric
bounds, not a thesis-printed pair.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

R_GAS = 8.314462618  # J/(mol K)
M_HEXANE = 0.08618  # kg/mol, n-hexane (C6H14)
M_WATER = 0.018015  # kg/mol
ATM_PRESSURE_PA = 101325.0


# ---------------------------------------------------------------------------
# Sorption isotherms (§7.4; Cardarelli & Crapiste 1996 eqs. [2]-[4]; Coletto eq. 7)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GabParams:
    """GAB isotherm coefficients. `C`/`K` are exponential (van 't Hoff-type) in
    temperature — NOT linear; `Xm` is temperature-independent (Cardarelli &
    Crapiste 1996, Table 2)."""

    Xm: float  # kg/kg dry solid, monolayer capacity (their Hm)
    C0: float
    dHC_R: float  # K, delta_H_C / R
    K0: float
    dHK_R: float  # K, delta_H_K / R


@dataclass(frozen=True)
class OilIsotherm:
    A0: float
    B: float


COMPOSITE_MEAL_REFERENCE_OIL_MASS_FRACTION = 0.0195
COMPOSITE_MEAL_MIN_OIL_MASS_FRACTION = 0.012
COMPOSITE_MEAL_MAX_OIL_MASS_FRACTION = 0.030


@dataclass(frozen=True)
class CompositeMealHexaneParams:
    """Reference-anchored retained-hexane closure for composite dry meal.

    ``oil_mass_fraction`` is kg residual oil per kg solvent- and water-free
    composite meal, including that oil.  This is deliberately not legacy
    ``X3 = kg oil / kg non-oil solid``.

    Cardarelli and Crapiste's native-soybean GAB is the *total-meal* anchor at
    ``reference_oil_mass_fraction``.  The oil power law supplies only the
    signed composition correction away from that anchor; it is never added
    once more at the reference composition.
    """

    native_gab: GabParams
    oil: OilIsotherm
    reference_oil_mass_fraction: float = COMPOSITE_MEAL_REFERENCE_OIL_MASS_FRACTION
    minimum_oil_mass_fraction: float = COMPOSITE_MEAL_MIN_OIL_MASS_FRACTION
    maximum_oil_mass_fraction: float = COMPOSITE_MEAL_MAX_OIL_MASS_FRACTION

    def __post_init__(self) -> None:
        values = (
            self.reference_oil_mass_fraction,
            self.minimum_oil_mass_fraction,
            self.maximum_oil_mass_fraction,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("composite-meal oil bounds must be finite")
        if not (
            0.0
            <= self.minimum_oil_mass_fraction
            <= self.reference_oil_mass_fraction
            <= self.maximum_oil_mass_fraction
            < 1.0
        ):
            raise ValueError("oil bounds must satisfy 0 <= minimum <= reference <= maximum < 1")


@dataclass(frozen=True)
class ProvisionalHexaneSorptionCaloricParams:
    """Legacy empirical constants for a provisional retained-hexane potential.

    ``dH_vap_hexane + sorption_C0*W**sorption_C1`` is the configured
    differential heat of sorption.  The low-loading continuation uses the
    same ``loading_floor_fraction * native_gab.Xm`` convention as the current
    particle model, but integrates it into a finite state function. This is a
    conservative adapter for the current empirical heat law, not a frozen
    GT-PS-2 binding-energy truth; the native potential remains open under
    PHY-021.
    """

    dH_vap_hexane: float
    sorption_C0: float
    sorption_C1: float
    loading_floor_fraction: float = 0.02

    def __post_init__(self) -> None:
        if (
            not math.isfinite(self.dH_vap_hexane)
            or self.dH_vap_hexane <= 0.0
            or not math.isfinite(self.sorption_C0)
            or self.sorption_C0 < 0.0
            or not math.isfinite(self.sorption_C1)
            or not math.isfinite(self.loading_floor_fraction)
            or self.loading_floor_fraction <= 0.0
        ):
            raise ValueError("invalid retained-hexane caloric parameters")


_GAB_MAX_KA = 0.999
"""Largest `K*a_h` the GAB multilayer term is evaluated at before its `(1 - K*a_h)`
denominator diverges. In normal DT operation `K*a_h` stays well below 1; this bound
only engages for a cold OFF-DESIGN transient (e.g. a weak/under-set sparge driving the
pore gas to hexane saturation, `a_h -> 1`, where `K > 1`)."""


def _gab_clamp_activity(a_h: float, K: float) -> float:
    """Clamp `a_h` so `K*a_h <= _GAB_MAX_KA`, keeping the GAB isotherm finite at its own
    divergence boundary. Graceful degradation (return the value at the largest valid
    activity) rather than raising -- a robustness guard for off-design/transient inputs,
    NOT a physics change: at every calibrated operating point `K*a_h` is far from 1, so
    this never engages there."""
    if K > 0.0 and K * a_h > _GAB_MAX_KA:
        return _GAB_MAX_KA / K
    return a_h


def gab_hexane_content(a_h: float, T: float, p: GabParams) -> float:
    """W2_eq(a_h, T) — hexane adsorbed on the particle solid phase, kg/kg dry
    solid (Cardarelli & Crapiste 1996, eq. [2], with C/K per eqs. [3]-[4]).

    Valid only for `K*a_h < 1` (the GAB model's multilayer term diverges at the
    boundary, same as the water-sorption GAB literature it's adapted from).
    """
    if not 0.0 <= a_h <= 1.0:
        raise ValueError(f"hexane activity a_h must be in [0,1], got {a_h}")
    C = p.C0 * math.exp(p.dHC_R / T)
    K = p.K0 * math.exp(p.dHK_R / T)
    a_h = _gab_clamp_activity(a_h, K)
    return p.Xm * C * K * a_h / ((1.0 - K * a_h) * (1.0 - K * a_h + C * K * a_h))


def gab_hexane_content_and_slope(
    a_h: float,
    T: float,
    p: GabParams,
    *,
    strict_domain: bool = False,
) -> tuple[float, float]:
    """`(W2_eq, dW2_eq/da_h)` at `(a_h, T)` -- the value from
    `gab_hexane_content` plus its EXACT analytic derivative w.r.t. `a_h`.

    Replaces the old central-finite-difference slope (h=1e-5) used throughout
    the DCZ particle scale: the GAB isotherm is a closed-form rational
    function, so its derivative is closed-form too. Verified to agree with the
    former central difference to ~1e-9 relative across the operating range --
    a pure speed/accuracy change (one isotherm evaluation instead of two/three,
    and no finite-difference truncation or roundoff), NOT a physics change.

    With ``strict_domain=False`` this preserves the legacy off-design clamp.
    The GT composite-meal path uses ``strict_domain=True``: it evaluates the
    exact GAB expression for ``K*a_h < 1`` and rejects the divergence instead
    of flattening the property and destroying its inverse/Jacobian."""
    if not 0.0 <= a_h <= 1.0:
        raise ValueError(f"hexane activity a_h must be in [0,1], got {a_h}")
    C = p.C0 * math.exp(p.dHC_R / T)
    K = p.K0 * math.exp(p.dHK_R / T)
    if strict_domain:
        if K * a_h >= 1.0:
            raise ValueError(
                "strict GAB domain requires K(T)*a_h < 1; clamping is not "
                "permitted in the ground-truth property path"
            )
    else:
        a_h = _gab_clamp_activity(a_h, K)
    # W2 = Xm*C*K*a / D, with D = (1 - Ka)(1 - Ka + CKa). Quotient rule on
    # N = Xm*C*K*a, D as above (d1 = 1 - Ka, d2 = 1 - Ka + CKa).
    N = p.Xm * C * K * a_h
    d1 = 1.0 - K * a_h
    d2 = 1.0 - K * a_h + C * K * a_h
    D = d1 * d2
    N_prime = p.Xm * C * K
    d1_prime = -K
    d2_prime = -K + C * K
    D_prime = d1_prime * d2 + d1 * d2_prime
    value = N / D
    slope = (N_prime * D - N * D_prime) / (D * D)
    return value, slope


def gab_hexane_slope(a_h: float, T: float, p: GabParams) -> float:
    """Exact ``dW2_eq/da_h`` without computing the unused loading value."""
    if not 0.0 <= a_h <= 1.0:
        raise ValueError(f"hexane activity a_h must be in [0,1], got {a_h}")
    C = p.C0 * math.exp(p.dHC_R / T)
    K = p.K0 * math.exp(p.dHK_R / T)
    a_h = _gab_clamp_activity(a_h, K)
    x = K * a_h
    denominator = (1.0 - x) * (1.0 - x + C * x)
    # d/da of Xm*C*x/((1-x)*(1+(C-1)*x)), where x=K*a.
    return p.Xm * C * K * (1.0 + (C - 1.0) * x * x) / (denominator * denominator)


def oil_hexane_content_and_slope(a_h: float, p: OilIsotherm) -> tuple[float, float]:
    """`(qo, dqo/da_h)` for the power-law oil isotherm `qo = A0*a_h^B` (eq. 7).
    Analytic companion to `oil_hexane_content`, same rationale as
    `gab_hexane_content_and_slope`. `dqo/da_h = A0*B*a_h^(B-1)` (0 at a_h=0 for
    the usual B>1)."""
    if not 0.0 <= a_h <= 1.0:
        raise ValueError(f"hexane activity a_h must be in [0,1], got {a_h}")
    value = p.A0 * a_h**p.B
    slope = p.A0 * p.B * a_h ** (p.B - 1.0) if a_h > 0.0 else 0.0
    return value, slope


def oil_hexane_slope(a_h: float, p: OilIsotherm) -> float:
    """Exact ``dqo/da_h`` without computing the unused oil loading."""
    if not 0.0 <= a_h <= 1.0:
        raise ValueError(f"hexane activity a_h must be in [0,1], got {a_h}")
    return p.A0 * p.B * a_h ** (p.B - 1.0) if a_h > 0.0 else 0.0


def hexane_activity_from_loading(W2: float, T: float, p: GabParams, a_h_max: float = 1.0) -> float:
    """Inverse of `gab_hexane_content`: the hexane activity `a_h` whose GAB
    equilibrium loading equals `W2` (kg/kg dry solid) at temperature `T`. The
    isotherm is monotone in `a_h` on (0, 1/K), so the root is unique; clamped to
    `(0, min(a_h_max, 1/K))`. Used by the DC hexane desorption driving force
    (`core/dc.py`): the solid's own equilibrium hexane partial pressure is
    `a_h * p_sat(T)`, the "escaping tendency" that collapses at low temperature."""
    if W2 <= 0.0:
        return 0.0
    K = p.K0 * math.exp(p.dHK_R / T)
    upper = min(a_h_max, 1.0 / K) - 1.0e-9
    if upper <= 1.0e-12:
        return 0.0
    if gab_hexane_content(upper, T, p) <= W2:
        return upper
    from scipy.optimize import brentq

    return brentq(lambda a: gab_hexane_content(a, T, p) - W2, 1.0e-12, upper)


def oil_hexane_content(a_h: float, p: OilIsotherm) -> float:
    """qo(a_h) = A0 * a_h^B — hexane absorbed in the oil phase, kg/kg oil (eq. 7)."""
    if not 0.0 <= a_h <= 1.0:
        raise ValueError(f"hexane activity a_h must be in [0,1], got {a_h}")
    return p.A0 * a_h**p.B


def validate_composite_meal_oil_mass_fraction(
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
) -> float:
    """Return an approved mass-property oil fraction, rejecting rather than clamping.

    Passing this check does not by itself energy-qualify the state; the
    equilibrium-derived caloric applicability decision remains open under
    PHY-021.
    """

    if (
        not math.isfinite(oil_mass_fraction)
        or oil_mass_fraction < p.minimum_oil_mass_fraction
        or oil_mass_fraction > p.maximum_oil_mass_fraction
    ):
        raise ValueError(
            "composite dry-meal oil mass fraction must lie in "
            f"[{p.minimum_oil_mass_fraction}, {p.maximum_oil_mass_fraction}], "
            f"got {oil_mass_fraction}"
        )
    return oil_mass_fraction


def _reference_anchored_composite_weights(
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
) -> tuple[float, float]:
    w_o = validate_composite_meal_oil_mass_fraction(oil_mass_fraction, p)
    denominator = 1.0 - p.reference_oil_mass_fraction
    return (
        (1.0 - w_o) / denominator,
        (w_o - p.reference_oil_mass_fraction) / denominator,
    )


def reference_anchored_composite_property(
    native_composite_value: float,
    pure_oil_value: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
) -> float:
    """Mass-consistent property anchored to measured native composite meal.

    If ``P_ref=(1-w_ref)*P_non-oil+w_ref*P_oil``, eliminating the unmeasured
    non-oil property gives ``P(w)=alpha*P_ref+beta*P_oil`` with the same
    ``alpha`` and ``beta`` used by the retained-hexane closure. At
    ``w=w_ref`` the result is exactly ``P_ref``: the native oil contribution
    is not counted twice.
    """

    if not math.isfinite(native_composite_value) or not math.isfinite(pure_oil_value):
        raise ValueError("reference-anchored component properties must be finite")
    alpha, beta = _reference_anchored_composite_weights(oil_mass_fraction, p)
    return alpha * native_composite_value + beta * pure_oil_value


def reference_anchored_composite_property_and_oil_fraction_slope(
    native_composite_value: float,
    pure_oil_value: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
) -> tuple[float, float]:
    """Anchored property and exact composition slope.

    This generic combiner is the frozen PHY-048 caloric contract: a future
    PHY-021 native binding potential can use it without adopting the current
    empirical ``C0/C1`` potential adapter.
    """

    value = reference_anchored_composite_property(
        native_composite_value, pure_oil_value, oil_mass_fraction, p
    )
    slope = (pure_oil_value - native_composite_value) / (1.0 - p.reference_oil_mass_fraction)
    return value, slope


def composite_retained_hexane_content_and_slopes(
    a_h: float,
    T: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
) -> tuple[float, float, float]:
    """Return retained hexane and exact activity/composition derivatives.

    The result is kg n-hexane per kg composite dry meal:
    ``W = alpha*W_native(a,T) + beta*q_o(a)``.

    The native-soybean GAB already includes its reference oil. Consequently
    ``q_o`` is only a signed off-reference correction, not an additional
    storage term. This pressure-agnostic function does not claim that every
    activity up to one is reachable at ``(T,P)``; use
    :func:`validate_pressure_feasible_composite_retention` when qualifying an
    inversion domain.
    """

    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("retained-hexane temperature must be positive and finite")
    alpha, beta = _reference_anchored_composite_weights(oil_mass_fraction, p)
    native, native_da = gab_hexane_content_and_slope(
        a_h,
        T,
        p.native_gab,
        strict_domain=True,
    )
    oil, oil_da = oil_hexane_content_and_slope(a_h, p.oil)
    value = alpha * native + beta * oil
    if value < 0.0:
        raise ValueError(
            "reference-anchored retained-hexane law is negative at this extrapolated state"
        )
    denominator = 1.0 - p.reference_oil_mass_fraction
    return (
        value,
        alpha * native_da + beta * oil_da,
        (oil - native) / denominator,
    )


def composite_retained_hexane_content(
    a_h: float,
    T: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
) -> float:
    """Retained n-hexane, kg/kg composite dry meal."""

    value, _, _ = composite_retained_hexane_content_and_slopes(a_h, T, oil_mass_fraction, p)
    return value


def composite_retained_hexane_activity_slope(
    a_h: float,
    T: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
) -> float:
    """Exact ``dW/da_h`` on the canonical composite-meal basis."""

    _, slope, _ = composite_retained_hexane_content_and_slopes(a_h, T, oil_mass_fraction, p)
    return slope


def composite_retained_hexane_oil_fraction_slope(
    a_h: float,
    T: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
) -> float:
    """Exact ``dW/dw_o`` at fixed activity and temperature."""

    _, _, slope = composite_retained_hexane_content_and_slopes(a_h, T, oil_mass_fraction, p)
    return slope


def minimum_composite_retained_hexane_activity_slope(
    T: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
    *,
    a_h_max: float,
    subdivisions: int = 24,
) -> float:
    """Numerically certify the minimum analytic ``dW/da_h`` on an interval.

    Bounded minimization over short subintervals prevents a coarse activity
    grid from silently missing the shallow negative pocket that invalidated
    the former 1.0 wt% envelope. This qualification helper is not a hot-path
    constitutive evaluation.
    """

    if not math.isfinite(a_h_max) or a_h_max < 0.0 or a_h_max > 1.0 or subdivisions < 1:
        raise ValueError("a_h_max must lie in [0,1] and subdivisions must be positive")
    validate_composite_meal_oil_mass_fraction(oil_mass_fraction, p)
    if a_h_max == 0.0:
        return composite_retained_hexane_activity_slope(0.0, T, oil_mass_fraction, p)

    from scipy.optimize import minimize_scalar

    def slope(activity: float) -> float:
        return composite_retained_hexane_activity_slope(activity, T, oil_mass_fraction, p)

    minimum = min(slope(0.0), slope(a_h_max))
    width = a_h_max / subdivisions
    for index in range(subdivisions):
        lower = index * width
        upper = (index + 1) * width
        result = minimize_scalar(
            slope,
            bounds=(lower, upper),
            method="bounded",
            options={"xatol": max(1.0e-13, width * 1.0e-10)},
        )
        if not result.success:
            raise RuntimeError("retained-hexane monotonicity minimization failed")
        minimum = min(minimum, result.fun)
    return minimum


def hexane_activity_from_composite_retained_content(
    retained_hexane_kg_kg_composite_dry: float,
    T: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
    *,
    a_h_max: float,
    certify_monotonicity: bool = False,
) -> float:
    """Invert the retained-hexane law on a prequalified activity domain.

    Call :func:`validate_pressure_feasible_composite_retention` once when a
    material/temperature/pressure domain is established, then reuse its
    ``a_h_max`` here. The hot inverse performs only a five-point derivative
    sanity check plus a strict bracket. Set ``certify_monotonicity=True`` for
    the expensive segmented minimization when no prior certification exists.
    The inverse never chooses an arbitrary root or clamps an extrapolated
    state.
    """

    if (
        not math.isfinite(retained_hexane_kg_kg_composite_dry)
        or retained_hexane_kg_kg_composite_dry < 0.0
    ):
        raise ValueError("retained-hexane loading must be finite and non-negative")
    if certify_monotonicity:
        minimum_slope = minimum_composite_retained_hexane_activity_slope(
            T, oil_mass_fraction, p, a_h_max=a_h_max
        )
        slopes = (minimum_slope,)
    else:
        if not math.isfinite(a_h_max) or a_h_max < 0.0 or a_h_max > 1.0:
            raise ValueError("a_h_max must lie in [0,1]")
        slopes = tuple(
            composite_retained_hexane_activity_slope(fraction * a_h_max, T, oil_mass_fraction, p)
            for fraction in (0.0, 0.25, 0.5, 0.75, 1.0)
        )
    if min(slopes) <= 0.0:
        raise ValueError(
            "retained-hexane law is not strictly monotone on the requested activity domain"
        )
    if retained_hexane_kg_kg_composite_dry == 0.0:
        return 0.0

    upper_loading = composite_retained_hexane_content(a_h_max, T, oil_mass_fraction, p)
    tolerance = 1.0e-12 * max(1.0, upper_loading)
    if retained_hexane_kg_kg_composite_dry > upper_loading + tolerance:
        raise ValueError("retained-hexane loading exceeds the requested activity domain")
    if retained_hexane_kg_kg_composite_dry >= upper_loading - tolerance:
        return a_h_max

    from scipy.optimize import brentq

    return brentq(
        lambda activity: (
            composite_retained_hexane_content(activity, T, oil_mass_fraction, p)
            - retained_hexane_kg_kg_composite_dry
        ),
        0.0,
        a_h_max,
    )


# ---------------------------------------------------------------------------
# Heat of sorption (eq. A.31)
# ---------------------------------------------------------------------------


def heat_of_sorption(W2: float, dH_lv2: float, sorption_C0: float, sorption_C1: float) -> float:
    """delta_H_s = delta_H_lv2 + C0*W2^C1 — net isosteric heat of sorption plus
    the latent heat of vaporization (eq. A.31)."""
    return dH_lv2 + sorption_C0 * W2**sorption_C1


def _native_hexane_sorption_potential_and_loading_slope(
    native_loading: float,
    p: CompositeMealHexaneParams,
    caloric: ProvisionalHexaneSorptionCaloricParams,
) -> tuple[float, float]:
    """Integral of the native-meal differential heat and its ``dPhi/dW``."""

    if not math.isfinite(native_loading) or native_loading < 0.0:
        raise ValueError("native retained-hexane loading must be finite and non-negative")
    floor = caloric.loading_floor_fraction * p.native_gab.Xm
    exponent = caloric.sorption_C1
    if abs(exponent + 1.0) < 1.0e-12:
        if native_loading <= floor:
            excess_integral = native_loading / floor
        else:
            excess_integral = 1.0 + math.log(native_loading / floor)
    elif native_loading <= floor:
        excess_integral = floor**exponent * native_loading
    else:
        floor_term = floor ** (exponent + 1.0)
        excess_integral = floor_term + (native_loading ** (exponent + 1.0) - floor_term) / (
            exponent + 1.0
        )
    potential = caloric.dH_vap_hexane * native_loading + caloric.sorption_C0 * excess_integral
    differential_heat = (
        caloric.dH_vap_hexane + caloric.sorption_C0 * max(native_loading, floor) ** exponent
    )
    return potential, differential_heat


def provisional_composite_retained_hexane_sorption_potential_and_slopes(
    a_h: float,
    T: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
    caloric: ProvisionalHexaneSorptionCaloricParams,
) -> tuple[float, float, float]:
    """Return finite binding potential and exact ``a_h``/``w_o`` slopes.

    The native total-meal potential is anchored at ``w_ref``. The signed oil
    correction is liquid-like, ``dH_vap_hexane*q_o``. Both use exactly the
    storage closure's ``alpha,beta`` weights, so retained hexane represented
    at vapor enthalpy can subtract this potential once without a duplicate
    latent or oil source.
    """

    if not math.isfinite(T) or T <= 0.0:
        raise ValueError("retained-hexane temperature must be positive and finite")
    alpha, beta = _reference_anchored_composite_weights(oil_mass_fraction, p)
    native, native_da = gab_hexane_content_and_slope(
        a_h,
        T,
        p.native_gab,
        strict_domain=True,
    )
    oil, oil_da = oil_hexane_content_and_slope(a_h, p.oil)
    native_potential, native_dpotential_dloading = (
        _native_hexane_sorption_potential_and_loading_slope(native, p, caloric)
    )
    oil_potential = caloric.dH_vap_hexane * oil
    denominator = 1.0 - p.reference_oil_mass_fraction
    return (
        alpha * native_potential + beta * oil_potential,
        alpha * native_dpotential_dloading * native_da + beta * caloric.dH_vap_hexane * oil_da,
        (oil_potential - native_potential) / denominator,
    )


def provisional_composite_retained_hexane_sorption_potential(
    a_h: float,
    T: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
    caloric: ProvisionalHexaneSorptionCaloricParams,
) -> float:
    """Finite retained-hexane sorption potential, J/kg composite dry meal."""

    potential, _, _ = provisional_composite_retained_hexane_sorption_potential_and_slopes(
        a_h, T, oil_mass_fraction, p, caloric
    )
    return potential


def provisional_composite_retained_hexane_sorption_potential_change(
    entering_state: tuple[float, float, float],
    leaving_state: tuple[float, float, float],
    p: CompositeMealHexaneParams,
    caloric: ProvisionalHexaneSorptionCaloricParams,
) -> float:
    """Exact endpoint potential change for ``(a_h,T,w_o)`` states.

    Using endpoint differences makes adjacent finite-volume increments
    telescope exactly, including simultaneous temperature and oil-composition
    changes.
    """

    entering = provisional_composite_retained_hexane_sorption_potential(*entering_state, p, caloric)
    leaving = provisional_composite_retained_hexane_sorption_potential(*leaving_state, p, caloric)
    return leaving - entering


def provisional_composite_retained_hexane_enthalpy_and_slopes(
    a_h: float,
    T: float,
    oil_mass_fraction: float,
    hexane_vapor_specific_enthalpy_j_kg: float,
    p: CompositeMealHexaneParams,
    caloric: ProvisionalHexaneSorptionCaloricParams,
) -> tuple[float, float, float]:
    """Retained-phase enthalpy contribution on a vapor-enthalpy datum.

    ``h_ret = W*h_hexane,vapor - Phi`` in J/kg composite dry meal. Slopes are
    at fixed ``T`` and fixed supplied vapor specific enthalpy.
    """

    if not math.isfinite(hexane_vapor_specific_enthalpy_j_kg):
        raise ValueError("hexane vapor specific enthalpy must be finite")
    loading, loading_da, loading_dw = composite_retained_hexane_content_and_slopes(
        a_h, T, oil_mass_fraction, p
    )
    potential, potential_da, potential_dw = (
        provisional_composite_retained_hexane_sorption_potential_and_slopes(
            a_h, T, oil_mass_fraction, p, caloric
        )
    )
    return (
        loading * hexane_vapor_specific_enthalpy_j_kg - potential,
        loading_da * hexane_vapor_specific_enthalpy_j_kg - potential_da,
        loading_dw * hexane_vapor_specific_enthalpy_j_kg - potential_dw,
    )


def composite_meal_heat_capacity_and_oil_fraction_slope(
    oil_mass_fraction: float,
    native_composite_cp_j_kg_k: float,
    pure_oil_cp_j_kg_k: float,
    p: CompositeMealHexaneParams,
) -> tuple[float, float]:
    """Reference-anchored dry-meal heat capacity and exact ``dcp/dw_o``.

    ``native_composite_cp_j_kg_k`` is the measured/native meal value (Coletto
    supplement ``CPps``), not the unmeasured non-oil matrix value.
    """

    if (
        not math.isfinite(native_composite_cp_j_kg_k)
        or native_composite_cp_j_kg_k <= 0.0
        or not math.isfinite(pure_oil_cp_j_kg_k)
        or pure_oil_cp_j_kg_k <= 0.0
    ):
        raise ValueError("composite and oil heat capacities must be positive and finite")
    heat_capacity = reference_anchored_composite_property(
        native_composite_cp_j_kg_k,
        pure_oil_cp_j_kg_k,
        oil_mass_fraction,
        p,
    )
    if heat_capacity <= 0.0:
        raise ValueError("reference-anchored composite heat capacity is non-positive")
    slope = (pure_oil_cp_j_kg_k - native_composite_cp_j_kg_k) / (
        1.0 - p.reference_oil_mass_fraction
    )
    return heat_capacity, slope


def composite_meal_sensible_enthalpy_j_kg(
    T: float,
    oil_mass_fraction: float,
    native_composite_cp_j_kg_k: float,
    pure_oil_cp_j_kg_k: float,
    p: CompositeMealHexaneParams,
    *,
    reference_temperature_k: float,
) -> float:
    """Composite dry-meal sensible enthalpy on one explicit datum."""

    if (
        not math.isfinite(T)
        or T <= 0.0
        or not math.isfinite(reference_temperature_k)
        or reference_temperature_k <= 0.0
    ):
        raise ValueError("sensible-enthalpy temperatures must be positive and finite")
    heat_capacity, _ = composite_meal_heat_capacity_and_oil_fraction_slope(
        oil_mass_fraction,
        native_composite_cp_j_kg_k,
        pure_oil_cp_j_kg_k,
        p,
    )
    return heat_capacity * (T - reference_temperature_k)


def temperature_from_composite_meal_sensible_enthalpy(
    sensible_enthalpy_j_kg: float,
    oil_mass_fraction: float,
    native_composite_cp_j_kg_k: float,
    pure_oil_cp_j_kg_k: float,
    p: CompositeMealHexaneParams,
    *,
    reference_temperature_k: float,
) -> float:
    """Exact inverse of :func:`composite_meal_sensible_enthalpy_j_kg`."""

    if not math.isfinite(sensible_enthalpy_j_kg):
        raise ValueError("composite-meal sensible enthalpy must be finite")
    if not math.isfinite(reference_temperature_k) or reference_temperature_k <= 0.0:
        raise ValueError("reference temperature must be positive and finite")
    heat_capacity, _ = composite_meal_heat_capacity_and_oil_fraction_slope(
        oil_mass_fraction,
        native_composite_cp_j_kg_k,
        pure_oil_cp_j_kg_k,
        p,
    )
    temperature = reference_temperature_k + sensible_enthalpy_j_kg / heat_capacity
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("composite-meal sensible enthalpy implies invalid temperature")
    return temperature


# ---------------------------------------------------------------------------
# Critical / equilibrium hexane content (eqs. 4-6)
# ---------------------------------------------------------------------------


def rho_hexane_liquid(T: float) -> float:
    """Daubert & Danner correlation for liquid n-hexane density, kg/m3 (eqs.
    B.11-B.12). `T` in K; 507.6 K is n-hexane's critical temperature."""
    f1 = 1.0 + (1.0 - T / 507.6) ** 0.27537
    return 61.034 / (0.26411**f1)


def rho_hexane_vapor(T: float, P: float = ATM_PRESSURE_PA) -> float:
    """Ideal-gas density of pure hexane vapor at temperature T (K) and
    pressure P (Pa) — same ideal-gas form as eq. B.3."""
    return P * M_HEXANE / (R_GAS * T)


def binary_vapor_hexane_mole_fraction(
    hexane_mass_fraction: float,
) -> float:
    """Convert binary water/hexane vapor mass fraction to mole fraction."""

    if not 0.0 <= hexane_mass_fraction <= 1.0:
        raise ValueError("hexane vapor mass fraction must lie in [0, 1]")
    hexane_moles = hexane_mass_fraction / M_HEXANE
    water_moles = (1.0 - hexane_mass_fraction) / M_WATER
    return hexane_moles / (hexane_moles + water_moles)


def hexane_activity_from_partial_pressure(
    partial_pressure_pa: float,
    temperature_k: float,
    antoine_hexane: AntoineParams,
    *,
    fugacity_coefficient_ratio: float = 1.0,
) -> float:
    """Return ``a_h=phi_h*p_h/(phi_h,sat*P_h,sat)``.

    The fugacity-coefficient ratio is one in the ideal-gas limit. Values
    above one are retained for supersaturated trial states; a phase solver,
    rather than this pure property function, owns admissibility.
    """

    if (
        partial_pressure_pa < 0.0
        or not math.isfinite(partial_pressure_pa)
        or temperature_k <= 0.0
        or not math.isfinite(temperature_k)
        or fugacity_coefficient_ratio <= 0.0
        or not math.isfinite(fugacity_coefficient_ratio)
    ):
        raise ValueError("invalid state for hexane activity")
    return (
        fugacity_coefficient_ratio
        * partial_pressure_pa
        / antoine_pressure_pa(temperature_k, antoine_hexane)
    )


def hexane_activity_from_binary_vapor_mass_fraction(
    hexane_mass_fraction: float,
    temperature_k: float,
    pressure_pa: float,
    antoine_hexane: AntoineParams,
    *,
    fugacity_coefficient_ratio: float = 1.0,
) -> float:
    """Map binary-vapor mass fraction to the isotherm's hexane activity."""

    if pressure_pa <= 0.0 or not math.isfinite(pressure_pa):
        raise ValueError("vapor pressure must be positive and finite")
    mole_fraction = binary_vapor_hexane_mole_fraction(hexane_mass_fraction)
    return hexane_activity_from_partial_pressure(
        mole_fraction * pressure_pa,
        temperature_k,
        antoine_hexane,
        fugacity_coefficient_ratio=fugacity_coefficient_ratio,
    )


def ideal_hexane_partial_density_kg_m3(
    partial_pressure_pa: float,
    temperature_k: float,
) -> float:
    """Ideal-gas hexane partial mass density, independent of diluent."""

    if (
        partial_pressure_pa < 0.0
        or not math.isfinite(partial_pressure_pa)
        or temperature_k <= 0.0
        or not math.isfinite(temperature_k)
    ):
        raise ValueError("invalid state for hexane partial density")
    return partial_pressure_pa * M_HEXANE / (R_GAS * temperature_k)


def ideal_hexane_saturation_density_kg_m3(
    temperature_k: float,
    antoine_hexane: AntoineParams,
) -> float:
    """Ideal-gas hexane density at its saturation pressure."""

    return ideal_hexane_partial_density_kg_m3(
        antoine_pressure_pa(temperature_k, antoine_hexane),
        temperature_k,
    )


def x2_critical(
    alpha_pg: float,
    rho_hexL: float,
    alpha_ps: float,
    rho_ps: float,
    empirical: float | None = None,
) -> float:
    """Legacy scalar critical-loading oracle.

    With ``empirical=None`` this is Coletto's raw unswollen geometric-pore
    expression.  Supplying ``empirical`` reproduces the former fixed-threshold
    implementation for historical tests only.  Ground-truth production code
    must use :mod:`dtdc_simulator.core.critical_volume`, which keeps the
    reference-anchored multiplier fixed while porosity, local liquid density,
    and the explicit phase ledger determine ``X_c(theta)``.
    """
    if empirical is not None:
        return empirical
    return (alpha_pg * rho_hexL) / (alpha_ps * rho_ps)


def ypg2_equilibrium(rho_hexV: float, alpha_ps: float, rho_ps: float) -> float:
    """Ypg2,eq (eq. 6) — dry-basis gas-phase hexane content when the pores hold
    only saturated hexane VAPOR."""
    return rho_hexV * (1.0 - alpha_ps) / (rho_ps * alpha_ps)


def x2_so_and_slope(
    a_h: float,
    T: float,
    X3: float,
    gab: GabParams,
    oil: OilIsotherm,
    h: float = 1.0e-5,
) -> tuple[float, float]:
    """X2,so(a_h) = W2(a_h,T) + X3*qo(a_h) (eq. A.26) and its slope dX2,so/da_h
    (the DCZ particle scale's `Ca`, eq. A.28, treats the local pore-gas mass
    fraction `wpg2` itself as the isotherm's activity variable `a_h` -- the
    same convention `x2_equilibrium`/`x2_critical` already use at the fixed
    point `a_h=1`, generalized here to an arbitrary local value).

    The slope is the EXACT analytic derivative of the two closed-form isotherms
    (`gab_hexane_content_and_slope`/`oil_hexane_content_and_slope`), which
    remain this module's single source of truth for the isotherm algebra. This
    replaced a former central finite difference (h=1e-5); verified to agree
    with it to ~1e-9 relative -- a speed/accuracy change, not a physics one.
    `h` is retained in the signature (unused) for backward compatibility with
    existing callers/tests that may still pass it.
    """
    del h  # retained for signature compatibility; the analytic slope needs no step
    W2, dW2 = gab_hexane_content_and_slope(a_h, T, gab)
    qo, dqo = oil_hexane_content_and_slope(a_h, oil)
    value = W2 + X3 * qo
    slope = dW2 + X3 * dqo
    return value, slope


def x2_so_slope(
    a_h: float,
    T: float,
    X3: float,
    gab: GabParams,
    oil: OilIsotherm,
) -> float:
    """Derivative-only companion to :func:`x2_so_and_slope`."""
    if not 0.0 <= a_h <= 1.0:
        raise ValueError(f"hexane activity a_h must be in [0,1], got {a_h}")
    C = gab.C0 * math.exp(gab.dHC_R / T)
    K = gab.K0 * math.exp(gab.dHK_R / T)
    clamped_a_h = _gab_clamp_activity(a_h, K)
    x = K * clamped_a_h
    denominator = (1.0 - x) * (1.0 - x + C * x)
    gab_slope = gab.Xm * C * K * (1.0 + (C - 1.0) * x * x) / (denominator * denominator)
    oil_slope = oil.A0 * oil.B * a_h ** (oil.B - 1.0) if a_h > 0.0 else 0.0
    return gab_slope + X3 * oil_slope


def x2_equilibrium(
    T: float,
    X3: float,
    gab: GabParams,
    oil: OilIsotherm,
    alpha_pg: float,
    alpha_ps: float,
    rho_ps: float,
    P: float = ATM_PRESSURE_PA,
) -> float:
    """X2,eq(T) (eq. 5) — hexane content when the pores hold only saturated
    hexane vapor, in equilibrium with adsorbed/absorbed hexane at that same
    saturation condition (a_h = 1 by definition of "pores saturated with gas
    hexane", per Coletto §2.3.2)."""
    W2_eq = gab_hexane_content(1.0, T, gab)
    qo_eq = oil_hexane_content(1.0, oil)
    rho_hexV = rho_hexane_vapor(T, P)
    Ypg2_eq = ypg2_equilibrium(rho_hexV, alpha_ps, rho_ps)
    return W2_eq + X3 * qo_eq + Ypg2_eq


# ---------------------------------------------------------------------------
# Thermophysical properties — Appendix B (implement exactly)
# ---------------------------------------------------------------------------


def rho_vip(yV1: float, yV2: float, T: float, P: float = ATM_PRESSURE_PA) -> float:
    """rho_Vip (eq. B.3) — interparticle vapor density, ideal-gas mixture of
    water (mol frac yV1) and hexane (mol frac yV2) at temperature T (K)."""
    return (yV1 * M_WATER + yV2 * M_HEXANE) * P / (R_GAS * T)


def rho_lmix(alpha_L: float, rho_L: float, alpha_V: float, rho_vip_: float) -> float:
    """rho_Lmix (eq. B.1) — cake/inter-particle mixture density in the PHZ."""
    return alpha_L * rho_L + alpha_V * rho_vip_


def rho_l(alpha_ps: float, rho_ps: float, X1: float, X2: float, X3: float) -> float:
    """rho_L (eq. B.2) — descending solid (wet meal) density."""
    return alpha_ps * rho_ps * (1.0 + X1 + X2 + X3)


def cp_l(
    w_water: float,
    w_hexane: float,
    w_oil: float,
    w_ds: float,
    cps: tuple[float, float, float, float],
) -> float:
    """CPL (eq. B.5) — descending-solid specific heat, mass-fraction-weighted
    sum over the 4 components (water, hexane, oil, dry solid), in that order,
    matching `cps = (cp_water_liquid, cp_hexane_liquid, cp_oil, cp_solid)`."""
    weights = (w_water, w_hexane, w_oil, w_ds)
    return sum(w * cp for w, cp in zip(weights, cps))


def mixture_cp_per_kg_dry_solid(
    X1: float,
    X2: float,
    X3: float,
    cp_water_liquid: float,
    cp_hexane_liquid: float,
    cp_oil: float,
    cp_solid: float,
) -> float:
    """Effective heat capacity of the wet solid stream, J/(kg dry solid . K)
    -- `cp_l` (eq. B.5) applied to the solid's own composition (X1/X2/X3 are
    kg-per-kg-dry-solid ratios, converted to the mass fractions `cp_l` itself
    expects). Promoted from `zones/phz.py`'s own private
    `_mixture_cp_per_kg_dry_solid` (that module's only prior caller, which now
    delegates here) so `core/balance.py`'s independent solid-side energy
    checks can share the exact same, already-tested constitutive formula
    without duplicating it -- a balance CHECK duplicating physics ad hoc would
    itself be a maintenance/drift risk, whereas this is pure composition
    -weighting, not a suspect balance term."""
    m_total = 1.0 + X1 + X2 + X3
    w_water, w_hexane, w_oil, w_ds = X1 / m_total, X2 / m_total, X3 / m_total, 1.0 / m_total
    cp_per_kg_wet = cp_l(
        w_water, w_hexane, w_oil, w_ds, (cp_water_liquid, cp_hexane_liquid, cp_oil, cp_solid)
    )
    return cp_per_kg_wet * m_total  # J/(kg WET . K) -> J/(kg DRY solid . K)


def cp_vip(w_water_vapor: float, w_hexane_vapor: float, cps: tuple[float, float]) -> float:
    """CPVip (eq. B.6) — intra-particle vapor specific heat, mass-fraction
    -weighted sum over the 2 vapor components (water, hexane), matching
    `cps = (cp_water_vapor, cp_hexane_vapor)`."""
    weights = (w_water_vapor, w_hexane_vapor)
    return sum(w * cp for w, cp in zip(weights, cps))


def nu_from_reynolds(Re: float, Pr: float) -> float:
    """Nuε = 0.6949 * Reε^0.579 * Pr_V^(1/3) -- Coletto's eq. B.7 EXACTLY as
    printed in the main-paper appendix (Faner 2008, from experimental data on
    oilseed-meal desolventizing). D3 (GROUNDING_MATRIX.md): a prior version
    replaced this with the canonical Ranz-Marshall `2 + 0.6*Re^0.5*Pr^(1/3)`
    on the mistaken premise that B.7 was "unrecoverable" -- it is printed in
    the paper. Restored here to stay Coletto-faithful. Note B.7 carries NO
    `+2` conduction floor: at low Reε it goes to zero, which the paper accepts;
    the packed-bed regime this is used in does not reach that limit."""
    return 0.6949 * Re**0.579 * Pr ** (1.0 / 3.0)


def hq_from_nu(Nu: float, r_P: float, k_V: float, alpha_V: float, alpha_L: float) -> float:
    """hQ from Nu_epsilon = 2*hQ*rP/kV * (alphaV/alphaL) (eq. B.8), solved for hQ."""
    return Nu * k_V * alpha_L / (2.0 * r_P * alpha_V)


def schmidt_number(mu_V: float, rho_V: float, D_HW: float) -> float:
    """Sc_p = muV / (rhoV * D_HW) (eq. B.10)."""
    return mu_V / (rho_V * D_HW)


def hm_from_hq(hQ: float, rho_V: float, cp_V: float, Pr: float, Sc: float) -> float:
    """hM from the Chilton-Colburn analogy (eq. B.9)."""
    return (hQ / (rho_V * cp_V)) * (Pr / Sc) ** (2.0 / 3.0)


# ---------------------------------------------------------------------------
# VLLE dew-point curve (§A.2.3 / §2.3.4) — V-SCAL (superheated) -> V-SAT (dew curve)
#
# The paper describes the physics (vapor enters FTRZ superheated, cools until it
# hits the dew point for its current composition, then evolves along the dew
# curve as water condenses) but does not spell out the solve algorithm or the
# enthalpy datum. Two decisions made here (confirmed with the user):
#
# 1. Dew point = single-condensable-component Raoult's law on WATER only
#    (y_water*P = P_water,sat(T_dew)), not the two-liquid azeotrope co-boiling
#    point — valid because the paper explicitly restricts FTRZ to "hexane content
#    below the azeotropic value," where the only liquid that forms is water and
#    hexane stays entirely in the vapor as an inert carrier. This also matches
#    the paper's own Fig. 8(b), where the reported FTRZ dew temperature dips
#    below 100 C as hexane dilutes the vapor and suppresses water's partial
#    pressure.
# 2. Enthalpy datum: liquid water at its own boiling point and liquid hexane at
#    its own boiling point are the zero-enthalpy references for the respective
#    vapor-phase components (standard humid-gas-type convention). This is a
#    reasonable, documented construction filling a gap the paper leaves
#    implicit -- not a value stated in the paper itself.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AntoineParams:
    A: float
    B: float
    C: float


def antoine_pressure_pa(T: float, p: AntoineParams) -> float:
    """log10(P[bar]) = A - B/(C+T) -> saturation pressure in Pa."""
    return (10.0 ** (p.A - p.B / (p.C + T))) * 1.0e5


def pressure_feasible_hexane_activity_limit(
    temperature_k: float,
    pressure_pa: float,
    antoine_hexane: AntoineParams,
    *,
    fugacity_coefficient_ratio: float = 1.0,
) -> float:
    """Largest single-vapor-phase hexane activity reachable at ``(T,P)``.

    Since ``p_h <= P``,
    ``a_h <= min(1, phi_h/phi_h,sat * P/P_h,sat(T))``. This helper makes the
    physical inversion interval explicit; it does not clamp a requested trial
    activity.
    """

    if (
        not math.isfinite(temperature_k)
        or temperature_k <= 0.0
        or not math.isfinite(pressure_pa)
        or pressure_pa <= 0.0
        or not math.isfinite(fugacity_coefficient_ratio)
        or fugacity_coefficient_ratio <= 0.0
    ):
        raise ValueError("invalid temperature, pressure, or fugacity ratio")
    return min(
        1.0,
        fugacity_coefficient_ratio
        * pressure_pa
        / antoine_pressure_pa(temperature_k, antoine_hexane),
    )


def validate_pressure_feasible_hexane_activity(
    a_h: float,
    temperature_k: float,
    pressure_pa: float,
    antoine_hexane: AntoineParams,
    *,
    fugacity_coefficient_ratio: float = 1.0,
) -> float:
    """Reject an activity outside the reachable single-vapor-phase domain."""

    if not math.isfinite(a_h) or a_h < 0.0:
        raise ValueError("hexane activity must be finite and non-negative")
    limit = pressure_feasible_hexane_activity_limit(
        temperature_k,
        pressure_pa,
        antoine_hexane,
        fugacity_coefficient_ratio=fugacity_coefficient_ratio,
    )
    if a_h > limit:
        raise ValueError(f"hexane activity {a_h} exceeds pressure-feasible limit {limit}")
    return a_h


def validate_pressure_feasible_composite_retention(
    temperature_k: float,
    pressure_pa: float,
    oil_mass_fraction: float,
    p: CompositeMealHexaneParams,
    antoine_hexane: AntoineParams,
    *,
    fugacity_coefficient_ratio: float = 1.0,
    subdivisions: int = 24,
) -> tuple[float, float]:
    """Validate strict invertibility on the pressure-feasible activity interval.

    Returns ``(a_h_max, minimum_dW_da_h)``. A non-monotone closure is rejected,
    but the caller is responsible for supplying a saturation-pressure relation
    valid at ``temperature_k``. The project-wide physical qualification remains
    separate from this local numerical validator.
    """

    a_h_max = pressure_feasible_hexane_activity_limit(
        temperature_k,
        pressure_pa,
        antoine_hexane,
        fugacity_coefficient_ratio=fugacity_coefficient_ratio,
    )
    minimum_slope = minimum_composite_retained_hexane_activity_slope(
        temperature_k,
        oil_mass_fraction,
        p,
        a_h_max=a_h_max,
        subdivisions=subdivisions,
    )
    if minimum_slope <= 0.0:
        raise ValueError(
            "composite retained-hexane law is not strictly monotone over the "
            "pressure-feasible activity domain"
        )
    return a_h_max, minimum_slope


def antoine_saturation_temperature_k(
    pressure_pa: float,
    p: AntoineParams,
) -> float:
    """Invert the Antoine equation at a positive saturation pressure."""

    if pressure_pa <= 0.0 or not math.isfinite(pressure_pa):
        raise ValueError("saturation pressure must be positive and finite")
    denominator = p.A - math.log10(pressure_pa / 1.0e5)
    if denominator == 0.0:
        raise ValueError("Antoine inverse is singular at this pressure")
    temperature = p.B / denominator - p.C
    if temperature <= 0.0 or not math.isfinite(temperature):
        raise ValueError("Antoine inverse returned a non-physical temperature")
    return temperature


def _y_water_mole_fraction(Y_V2: float) -> float:
    """Y_V2 = kg hexane / kg water (water basis, eq. A.2b) -> mole fraction of
    water in the (water + hexane) vapor."""
    mole_ratio_hex_to_water = Y_V2 * M_WATER / M_HEXANE
    return 1.0 / (1.0 + mole_ratio_hex_to_water)


def dew_point_temperature(
    Y_V2: float,
    antoine_water: AntoineParams,
    P: float = ATM_PRESSURE_PA,
    T_bounds: tuple[float, float] = (230.0, 450.0),
) -> float:
    """T_dew(Y_V2) solving y_water*P = P_water,sat(T_dew) (Raoult's law, water
    condensing from a hexane-diluted vapor acting as an inert carrier).

    Closed-form inverse of the Antoine equation, replacing a former `brentq`
    root find (which was called ~10k times per DT solve). Antoine is
    `log10(P[bar]) = A - B/(C+T)`, so with `target` in Pa the dew point is
    `T = B/(A - log10(target/1e5)) - C` directly -- exact, no iteration.
    Clamped to `T_bounds` to preserve the old bracketed solver's range
    (unchanged behaviour wherever the root was already interior, which is
    everywhere this is reached in practice)."""
    lo, hi = T_bounds
    target = _y_water_mole_fraction(Y_V2) * P
    # Guard: a non-physical (negative/zero) water partial pressure -- e.g. from a
    # transient negative Y_V2 mid-iteration -- has no real dew point; return the
    # lower bound rather than crash `log10`. (`brentq` used to raise here; the
    # closed form must degrade gracefully so a bad transient can't kill the solve.)
    if target <= 0.0:
        return lo
    T = antoine_saturation_temperature_k(target, antoine_water)
    return min(max(T, lo), hi)


def water_activity(
    Y_V2: float, T: float, antoine_water: AntoineParams, P: float = ATM_PRESSURE_PA
) -> float:
    """a_w = y_water*P / P_water,sat(T) -- water's own partial pressure
    (Raoult's law, same convention as `dew_point_temperature`) relative to
    its saturation pressure at the LOCAL temperature. Deliberately NOT
    clamped to [0,1]: `a_w >= 1` is mathematically the same condition as
    `T <= dew_point_temperature(Y_V2, ...)` (both say "water can't stay
    vapor here") -- callers route that case to condensation, not this
    isotherm (see `luikov_equilibrium_moisture`'s own docstring)."""
    return _y_water_mole_fraction(Y_V2) * P / antoine_pressure_pa(T, antoine_water)


@dataclass(frozen=True)
class LuikovParams:
    """Modified LUIKOV (1978) desorption isotherm coefficients -- Gianini,
    Luz, Sousa, Jorge & Paraíso (2006), *Ciênc. Tecnol. Aliment.* 26(2):
    408-413, Table 7, fit to soybean meal sampled DIRECTLY from a
    desolventizer/toaster's own outlet (not a generic food-science
    isotherm) across a combined 15-70 °C dataset -- temperature-INDEPENDENT
    by design (the paper's own finding: T has no significant effect on
    equilibrium moisture in that range, confirmed both graphically and
    statistically, Tables 5 vs. 7)."""

    A1: float
    A2: float


LUIKOV_MAX_VALIDATED_UR = 0.799
"""Gianini et al.'s own highest tested water activity (their KCl saturated
-salt-solution data point) -- callers must clamp `a_w` here before calling
`luikov_equilibrium_moisture` (see that function's own docstring): the fitted
curve climbs steeply toward its asymptote `A1` as `a_w -> 1` with zero
supporting data past this point (confirmed: evaluating it unclamped where the
local vapor phase sits close to saturation gives `Xe > 0.5`, a pure
extrapolation artifact, not a real equilibrium). Shared here (not left as a
private per-caller constant) because more than one caller now needs the exact
same bound (`core/zones/dcz.py`'s DCZ moisture balance, `core/dc.py`'s
dryer/cooler air-contact isotherm)."""


def luikov_equilibrium_moisture(a_w: float, p: LuikovParams) -> float:
    """Xe(a_w) = A1 / (1 + A2*ln(1/a_w)) -- equilibrium solid moisture (kg
    water/kg dry solid) at water activity `a_w`. Valid domain is (0, 1),
    matching the cited paper's own tested range (its own UR never exceeded
    ~0.8) -- callers must keep `a_w` there themselves (e.g. `min(a_w, 1.0 -
    eps)`); `a_w >= 1` is a genuinely different regime (condensation, not
    sorption) with its own separate mechanism, not this isotherm extrapolated
    past its own domain.

    EXTRAPOLATION CAVEAT, stated not hidden: the cited paper's own
    temperature range (15-70 °C / 288-343 K) sits below this project's own
    DT-internal operating temperatures (DCZ currently reaches ~380 K+). Its
    finding that temperature barely matters is reassuring but doesn't cover
    that gap -- applying it there is a documented assumption, not a
    validated one.
    """
    return p.A1 / (1.0 + p.A2 * math.log(1.0 / a_w))


# ---------------------------------------------------------------------------
# Soybean-meal air-drying correlations (Luz et al.) -- used by the DC
# (dryer/cooler) stage, `core/dc.py`. NOT for the DT-internal DCZ zone, which
# keeps the Gianini/Luikov desorption isotherm above: these are purpose-built
# for the falling-rate AIR-drying regime (25-100 C, ambient/heated air over a
# hygroscopic bed), the exact regime the DT-internal (superheated-vapor,
# hexane-dominated) sorption model does NOT cover.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LuzDryingParams:
    """Soybean-meal air-drying correlations from Luz, dos Santos Conceição,
    de Matos Jorge, Paraíso & Andrade (2010) *Food Bioprod. Process.* 88:
    90-98 ("Dynamic modeling and control of soybean meal drying in a direct
    rotary dryer"), themselves citing Luz et al. (2006a/c).

    MASS-TRANSFER COEFFICIENT K [1/s] -- Luz eq. (4), the FRACTIONAL-moisture
    form (X_s in kg water/kg dry solid, not %):

        K = (k_a2*T_a + k_b2)*X_s^2 + (k_a1*T_a + k_b1)*X_s + k_c

    Both Luz and Silva et al. (2012) establish soybean-meal drying is
    ENTIRELY falling-rate (internal-diffusion-controlled) -- there is no
    constant-rate period -- so this rate, NOT the air's saturation capacity,
    sets how fast moisture leaves. The temperature terms are near-negligible
    next to `k_c` (K ~ 8.44e-3 /s over the whole band), so evaluating `T_a`
    in K here rather than Luz's own °C shifts K by <0.1% -- immaterial, and
    kept in K for SI consistency with the rest of `core/`.

    EQUILIBRIUM MOISTURE X_e [kg water/kg dry solid] -- Luz eq. (5), a
    modified-Halsey isotherm (Luz et al. 2006a):

        X_e = xe_num / (1 + xe_coef * T_s * ln(1/ur))

    with `T_s` the SOLID temperature (K) and `ur` the local air relative
    humidity (0-1). Unlike the temperature-INDEPENDENT Gianini/Luikov
    isotherm (`LuikovParams`, fit to desolventizer-outlet meal at 15-70 C),
    this one is explicitly temperature-dependent and was regressed against
    the drying process itself -- the appropriate `X_e` for the DC regime."""

    k_a2: float  # 1/(s.K), coeff of T_a in the X_s^2 term        (Luz: -0.33e-11 in 1/(s.°C))
    k_b2: float  # 1/s, constant of the X_s^2 term                (Luz: 4.60e-9)
    k_a1: float  # 1/(s.K), coeff of T_a in the X_s term          (Luz: 7e-8 in 1/(s.°C))
    k_b1: float  # 1/s, constant of the X_s term                  (Luz: 1.42e-5)
    k_c: float  # 1/s, constant term (dominant)                   (Luz: 8.44e-3)
    xe_num: float  # kg/kg dry solid, isotherm numerator          (Luz: 0.834)
    xe_coef: float  # 1/K, isotherm temperature/activity factor   (Luz: 0.036)


def luz_mass_transfer_coefficient(T_a: float, X_s: float, p: LuzDryingParams) -> float:
    """K(T_a, X_s) [1/s] -- Luz eq. (4). `T_a` air temperature (K), `X_s`
    solid moisture (kg water/kg dry solid). Clamped at `X_s >= 0`; the
    quadratic is monotone-safe over the physical moisture range."""
    x = max(X_s, 0.0)
    return (p.k_a2 * T_a + p.k_b2) * x * x + (p.k_a1 * T_a + p.k_b1) * x + p.k_c


def luz_equilibrium_moisture(T_s: float, ur: float, p: LuzDryingParams) -> float:
    """X_e(T_s, ur) [kg water/kg dry solid] -- Luz eq. (5). `T_s` solid
    temperature (K), `ur` local air relative humidity. `ur` is clamped to
    (0, 1) open interval: `ur -> 1` gives the isotherm's own moisture ceiling
    (`xe_num`), `ur -> 0` gives bone-dry equilibrium, both finite."""
    ur_c = min(max(ur, 1.0e-6), 1.0 - 1.0e-6)
    return p.xe_num / (1.0 + p.xe_coef * T_s * math.log(1.0 / ur_c))


@dataclass(frozen=True)
class VaporEnthalpyRef:
    """Zero-enthalpy datum for the water-basis vapor enthalpy curve: liquid
    water at its own bp, liquid hexane at its own bp."""

    dH_vap_water: float
    cp_water_vapor: float
    T_boil_water: float
    dH_vap_hexane: float
    cp_hexane_vapor: float
    T_boil_hexane: float


def vapor_enthalpy_water_basis(Y_V2: float, T_V: float, ref: VaporEnthalpyRef) -> float:
    """Hvbw = f(YV2, TV) (eq. A.5) in the superheated (V-SCAL) regime: a
    mass-fraction-weighted (water-basis) sum of the water and hexane vapor
    specific enthalpies, each referenced to its own liquid at its own bp."""
    H_water = ref.dH_vap_water + ref.cp_water_vapor * (T_V - ref.T_boil_water)
    H_hexane = ref.dH_vap_hexane + ref.cp_hexane_vapor * (T_V - ref.T_boil_hexane)
    return H_water + Y_V2 * H_hexane


def liquid_water_enthalpy(
    T: float,
    ref: VaporEnthalpyRef,
    cp_water_liquid: float,
) -> float:
    """Specific liquid-water enthalpy on ``VaporEnthalpyRef``'s datum."""
    return cp_water_liquid * (T - ref.T_boil_water)


def water_phase_change_enthalpy(
    T: float,
    ref: VaporEnthalpyRef,
    cp_water_liquid: float,
) -> float:
    """Return ``h_water,vapor(T) - h_water,liquid(T)`` in J/kg.

    This is the energy paired with liquid-to-vapor water transfer at one
    temperature.  Evaluating both phase enthalpies on the same datum gives
    the constant-heat-capacity integral of Kirchhoff's relation,
    ``dH_vap(T_boil) + (cp_v - cp_l) * (T - T_boil)``.
    """
    h_vapor = ref.dH_vap_water + ref.cp_water_vapor * (T - ref.T_boil_water)
    h_liquid = liquid_water_enthalpy(T, ref, cp_water_liquid)
    return h_vapor - h_liquid


def temperature_from_vapor_enthalpy(H_vbw: float, Y_V2: float, ref: VaporEnthalpyRef) -> float:
    """Inverse of `vapor_enthalpy_water_basis`: T_V from (Hvbw, YV2). Closed
    -form since the enthalpy is linear in T_V (BuildSpec §7.4 asks for this
    inverse alongside the forward curve)."""
    a = (ref.dH_vap_water - ref.cp_water_vapor * ref.T_boil_water) + Y_V2 * (
        ref.dH_vap_hexane - ref.cp_hexane_vapor * ref.T_boil_hexane
    )
    b = ref.cp_water_vapor + Y_V2 * ref.cp_hexane_vapor
    return (H_vbw - a) / b


def dew_point_enthalpy_water_basis(
    Y_V2: float,
    antoine_water: AntoineParams,
    ref: VaporEnthalpyRef,
    P: float = ATM_PRESSURE_PA,
) -> float:
    """Hvbw = f(YV2) (eq. A.5) once the vapor is on the dew curve (V-SAT): the
    same enthalpy curve evaluated at T = T_dew(YV2), collapsing the remaining
    degree of freedom exactly as §A.2.3 describes."""
    T_dew = dew_point_temperature(Y_V2, antoine_water, P)
    return vapor_enthalpy_water_basis(Y_V2, T_dew, ref)
