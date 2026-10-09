r"""Sorption isotherms and the equilibrium-derived retained-phase caloric
potentials for the resolved particle (Phase 1).

Frozen decisions implemented here (the *selection* of each, per
`release/physics_decisions.yaml`):

- PHY-043 native-soybean total-meal GAB (Cardarelli & Crapiste 1996): retained
  n-hexane equilibrium W_native(a_h, T); no duplicate explicit oil term at the
  native anchor.
- PHY-048 residual-oil composite closure: W_total = alpha_o*W_native +
  beta_o*q_o(a_h), q_o = 0.9635*a^2.7036, over 0.012 <= w_o <= 0.030; the
  same weights give cp_dm = alpha_o*cp_native + beta_o*cp_oil.  The native
  anchor already contains its native oil, so no separate w_o*cp_oil term is
  added.  States outside the envelope are REJECTED, not clamped.
- PHY-020 local-equilibrium storage: retained n-hexane is an algebraic function
  of activity/T/w_o; activity is recovered by monotone inversion (no fitted
  site-rate state).
- PHY-021.hexane equilibrium-derived retained potential: the chemical-potential
  primitive psi, the net isosteric heat q_net, the binding deficit B, and the
  retained enthalpy/internal energy H_ret/U_ret, on ONE Span-Wagner component
  datum (PHY-049) with the PHY-051 liquid-like retained partial volume. No
  independent sorption-heat or latent source.
- PHY-031 / PHY-021.water modified-Luikov retained water: W_w(a_w) with zero
  nominal net excess binding enthalpy (calorically liquid-like), evidence cap.

Prohibitions honored: no activity clamp (the legacy `_gab_clamp_activity` guard
is prohibited here per PHY-021 — clamping destroys the derivative/inverse
identities); reject a_h*K(T) >= 1 and out-of-envelope w_o instead.

Numerics note (Phase 1 = correctness-first): psi and B use adaptive quadrature
over the inverse isotherm and are cross-checked against the Gibbs-Helmholtz
identity `B = (1/M_h) * integral_0^W q_net dxi`. The analytic partial-fraction
GAB potential (PHY-021's "no hot-loop quadrature" runtime kernel) is a Phase-4
optimization; the identities it must satisfy are gated here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from dtdc_simulator.core2.props import hexane as hx

R = 8.31451              # J/(mol K)  (matches the hexane authority)
W_O_REF = 0.0195         # PHY-040/048 native residual-oil mass fraction
W_O_MIN = 0.012          # PHY-048 composition-applicability envelope
W_O_MAX = 0.030
CP_NATIVE_DRY_MEAL = 2317.0  # J/(kg dry composite meal K), native anchor
CP_OIL = 2000.0              # J/(kg oil K), signed off-anchor correction


# ---------------------------------------------------------------------------
# Parameters (introspectable; seeded with the frozen soybean values)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class GabParams:
    """Cardarelli & Crapiste (1996) native-soybean GAB (properties/soybean.yaml)."""

    Xm: float = 5.183e-3     # kg/kg dry, monolayer capacity
    C0: float = 3.117e-3
    dHC_R: float = 2262.0    # K
    K0: float = 9.172e-2
    dHK_R: float = 729.6     # K

    def C(self, T: float) -> float:
        return self.C0 * math.exp(self.dHC_R / T)

    def K(self, T: float) -> float:
        return self.K0 * math.exp(self.dHK_R / T)


@dataclass(frozen=True)
class OilIsotherm:
    """Hexane-in-oil power law q_o = A0*a^B (Coletto suppl. / PHY-048)."""

    A0: float = 0.9635
    B: float = 2.7036

    def q(self, a: float) -> float:
        return self.A0 * a**self.B

    def dq_da(self, a: float) -> float:
        return self.A0 * self.B * a ** (self.B - 1.0)


@dataclass(frozen=True)
class LuikovParams:
    """Modified-Luikov retained-water isotherm (PHY-031)."""

    A1: float = 0.880
    A2: float = 12.184
    W_cap: float = 0.2356709040   # evidence cap W_w(a_w=0.799)
    W_ref: float = 0.023          # qualified lower moisture anchor


@dataclass(frozen=True)
class ContinuedPositiveLuikovParams(LuikovParams):
    """Owner-approved study extrapolation to positive water below A5's anchor.

    Authorization: GT_PS2_PARTICLE_WATER_CONTINUATION_RULING_2026-10-06.md.
    Zero is an activity-map boundary, never an admitted wet material state.
    The original coefficients and upper evidence cap are retained.
    """

    W_ref: float = 0.0

    def __post_init__(self):
        if self.W_ref != 0.0:
            raise ValueError("continued positive-water branch has an open zero lower endpoint")


def water_log_activity(W_w: float, p: LuikovParams) -> float:
    """Positive-water potential without exponential underflow in the study branch."""
    if isinstance(p, ContinuedPositiveLuikovParams):
        if not math.isfinite(W_w) or not 0.0 < W_w <= p.W_cap:
            raise ValueError("continued water potential requires positive loading inside the cap")
        return (1.0 - p.A1 / W_w) / p.A2
    return math.log(water_activity(W_w, p))


# ---------------------------------------------------------------------------
# PHY-043 native GAB isotherm and derivatives
# ---------------------------------------------------------------------------
def gab_native(a: float, T: float, p: GabParams) -> float:
    """W_native(a, T), kg n-hexane/kg dry composite meal (no clamp)."""
    C, K = p.C(T), p.K(T)
    x = K * a
    if x >= 1.0:
        raise ValueError(
            f"GAB singular: K(T)*a_h = {x:.4f} >= 1 (a_h={a}, T={T}); reject, do not clamp"
        )
    return p.Xm * C * x / ((1.0 - x) * (1.0 - x + C * x))


def gab_native_da(a: float, T: float, p: GabParams) -> float:
    """d W_native/d a at fixed T (exact)."""
    C, K = p.C(T), p.K(T)
    x = K * a
    denom = (1.0 - x) * (1.0 - x + C * x)
    # d/da of Xm*C*x/((1-x)(1-x+Cx)), x=Ka  ->  Xm*C*K*(1 + (C-1)x^2)/denom^2
    return p.Xm * C * K * (1.0 + (C - 1.0) * x * x) / (denom * denom)


def gab_native_dT(a: float, T: float, p: GabParams) -> float:
    """d W_native/d T at fixed a (exact, via C(T)/K(T) van't Hoff derivatives)."""
    C, K = p.C(T), p.K(T)
    dC = -C * p.dHC_R / T**2
    dK = -K * p.dHK_R / T**2
    x = K * a
    d1 = 1.0 - x
    d2 = 1.0 - x + C * x
    N = p.Xm * C * x
    D = d1 * d2
    # dN/dT and dD/dT with x=Ka -> dx/dT = a*dK
    dx = a * dK
    dN = p.Xm * (dC * x + C * dx)
    dd1 = -dx
    dd2 = -dx + dC * x + C * dx
    dD = dd1 * d2 + d1 * dd2
    return (dN * D - N * dD) / (D * D)


def gab_native_dTT(a: float, T: float, p: GabParams) -> float:
    """d2 W_native/d T2 at fixed activity (exact).

    The wet-core caloric Jacobian needs this second derivative because the
    saturated retained/mobile partition changes with temperature while total
    shell n-hexane remains fixed.  Keeping it analytic avoids hiding a finite-
    difference derivative inside the energy Newton solve.
    """
    C, K = p.C(T), p.K(T)
    C_T = -C * p.dHC_R / T**2
    C_TT = C * (2.0 * p.dHC_R / T**3 + p.dHC_R**2 / T**4)
    x = K * a
    x_T = -x * p.dHK_R / T**2
    x_TT = x * (2.0 * p.dHK_R / T**3 + p.dHK_R**2 / T**4)

    N = p.Xm * C * x
    N_T = p.Xm * (C_T * x + C * x_T)
    N_TT = p.Xm * (C_TT * x + 2.0 * C_T * x_T + C * x_TT)

    d1 = 1.0 - x
    d1_T = -x_T
    d1_TT = -x_TT
    d2 = 1.0 - x + C * x
    d2_T = -x_T + C_T * x + C * x_T
    d2_TT = -x_TT + C_TT * x + 2.0 * C_T * x_T + C * x_TT
    D = d1 * d2
    D_T = d1_T * d2 + d1 * d2_T
    D_TT = d1_TT * d2 + 2.0 * d1_T * d2_T + d1 * d2_TT

    return (
        N_TT / D
        - 2.0 * N_T * D_T / D**2
        - N * D_TT / D**2
        + 2.0 * N * D_T**2 / D**3
    )


# ---------------------------------------------------------------------------
# PHY-048 composite (oil-dependent) retained hexane
# ---------------------------------------------------------------------------
def _oil_weights(w_o: float) -> tuple[float, float]:
    if not (W_O_MIN - 1e-12 <= w_o <= W_O_MAX + 1e-12):
        raise ValueError(
            f"w_o={w_o} outside the frozen PHY-048 envelope [{W_O_MIN}, {W_O_MAX}]; reject, do not clamp"
        )
    alpha = (1.0 - w_o) / (1.0 - W_O_REF)
    beta = (w_o - W_O_REF) / (1.0 - W_O_REF)
    return alpha, beta


def composite_dry_meal_heat_capacity(
    w_o: float,
    cp_native: float = CP_NATIVE_DRY_MEAL,
    cp_oil: float = CP_OIL,
) -> float:
    r"""Reference-anchored PHY-048 dry-meal heat capacity, J/(kg K).

    ``cp_native`` is the heat capacity of the *whole native composite meal*
    at ``W_O_REF``.  Therefore the mass-consistent construction is

    ``cp_dm = alpha_o*cp_native + beta_o*cp_oil``,

    not ``cp_native + w_o*cp_oil``.  Below the anchor ``beta_o`` is a signed
    correction, exactly as it is in the frozen retention and caloric-potential
    closures.
    """
    if not all(math.isfinite(value) and value > 0.0 for value in (cp_native, cp_oil)):
        raise ValueError("native-meal and oil heat capacities must be positive and finite")
    alpha, beta = _oil_weights(w_o)
    capacity = alpha * cp_native + beta * cp_oil
    if not math.isfinite(capacity) or capacity <= 0.0:
        raise ValueError("composite dry-meal heat capacity must be positive and finite")
    return capacity


def composite_dry_meal_sensible_energy(
    T: float,
    T_ref: float,
    w_o: float,
    cp_native: float = CP_NATIVE_DRY_MEAL,
    cp_oil: float = CP_OIL,
) -> float:
    """PHY-048 sensible state function, J/kg dry composite meal."""
    if not all(math.isfinite(value) for value in (T, T_ref)):
        raise ValueError("dry-meal sensible-energy temperatures must be finite")
    return composite_dry_meal_heat_capacity(w_o, cp_native, cp_oil) * (T - T_ref)


def retained_hexane(a: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm) -> float:
    """W_total(a, T, w_o) — PHY-048 reference-anchored composite retention."""
    alpha, beta = _oil_weights(w_o)
    return alpha * gab_native(a, T, p) + beta * oil.q(a)


def retained_hexane_da(
    a: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm
) -> float:
    alpha, beta = _oil_weights(w_o)
    return alpha * gab_native_da(a, T, p) + beta * oil.dq_da(a)


def retained_hexane_dT(
    a: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm
) -> float:
    alpha, _ = _oil_weights(w_o)
    # oil term q_o(a) has no explicit T dependence (PHY-048: liquid-like, zero excess)
    return alpha * gab_native_dT(a, T, p)


def retained_hexane_dTT(
    a: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm
) -> float:
    """d2 W_total/d T2 at fixed activity (exact)."""
    alpha, _ = _oil_weights(w_o)
    return alpha * gab_native_dTT(a, T, p)


def activity_max(T: float, p: GabParams) -> float:
    """Largest admissible activity from the GAB singularity K(T)*a < 1."""
    return 1.0 / p.K(T)


def activity_from_retained(
    W: float,
    T: float,
    w_o: float,
    p: GabParams,
    oil: OilIsotherm,
    *,
    a_ceiling: float | None = None,
) -> float:
    """Monotone inverse a_h(W, T, w_o) by bracketed bisection (no clamp).

    ``a_ceiling`` optionally bounds the search to a pressure-feasible activity
    (a_h*Psat_h(T) <= P); otherwise the GAB singularity ceiling is used.
    """
    if W < 0.0:
        raise ValueError(f"retained loading must be non-negative, got {W}")
    if W == 0.0:
        return 0.0
    hi = activity_max(T, p) * (1.0 - 1e-9)
    if a_ceiling is not None:
        hi = min(hi, a_ceiling)
    W_hi = retained_hexane(hi, T, w_o, p, oil)
    if W > W_hi:
        raise ValueError(
            f"retained loading {W} exceeds the admissible maximum {W_hi} at "
            f"T={T}, w_o={w_o}; reject, do not clamp"
        )
    lo, hi_a = 0.0, hi
    for _ in range(200):
        mid = 0.5 * (lo + hi_a)
        if retained_hexane(mid, T, w_o, p, oil) < W:
            lo = mid
        else:
            hi_a = mid
        if hi_a - lo < 1e-14 * max(hi_a, 1e-12):
            break
    return 0.5 * (lo + hi_a)


# ---------------------------------------------------------------------------
# PHY-021 retained-hexane caloric potential (equilibrium-derived, one datum)
# ---------------------------------------------------------------------------
def q_net_isosteric(
    a: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm
) -> float:
    """Net excess (isosteric) heat q_net(W,T), J/mol hexane, at the activity a.

    q_net = -R*T^2 * (dW/dT)_a / (a * (dW/da)_T)   (PHY-021 / A2b).
    Finite and positive as W -> 0 (native-GAB zero-loading limit).
    """
    dW_dT = retained_hexane_dT(a, T, w_o, p, oil)
    dW_da = retained_hexane_da(a, T, w_o, p, oil)
    return -R * T * T * dW_dT / (a * dW_da)


def _gauss_legendre_5():
    # nodes/weights on [-1, 1]
    n = (
        -0.906179845938664, -0.538469310105683, 0.0,
        0.538469310105683, 0.906179845938664,
    )
    w = (
        0.236926885056189, 0.478628670499366, 0.568888888888889,
        0.478628670499366, 0.236926885056189,
    )
    return n, w


def _integrate_over_loading(f, W: float, T: float, w_o: float, p: GabParams,
                            oil: OilIsotherm, n_panels: int = 24) -> float:
    """Composite 5-point Gauss-Legendre of f(a(xi)) over loading xi in [0, W]."""
    if W <= 0.0:
        return 0.0
    nodes, weights = _gauss_legendre_5()
    h = W / n_panels
    total = 0.0
    for k in range(n_panels):
        xi0 = k * h
        for xi_n, wt in zip(nodes, weights):
            xi = xi0 + 0.5 * h * (xi_n + 1.0)
            a = activity_from_retained(xi, T, w_o, p, oil)
            total += 0.5 * h * wt * f(a, T, w_o, p, oil)
    return total


def psi_potential(
    W: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm
) -> float:
    """Excess-Gibbs chemical-potential primitive psi(T,W,w_o), J/kg dry meal.

    psi = (R*T/M_h) * integral_0^W ln a_h(xi,T,w_o) dxi.
    """
    def integrand(a, T, w_o, p, oil):
        return math.log(a)

    integral = _integrate_over_loading(integrand, W, T, w_o, p, oil)
    return (R * T / hx.M) * integral


def binding_deficit(
    W: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm
) -> float:
    """Integrated binding deficit B(T,W,w_o), J/kg dry meal.

    B = (1/M_h) * integral_0^W q_net(xi) dxi   (Gibbs-Helmholtz consistent form).
    """
    def integrand(a, T, w_o, p, oil):
        return q_net_isosteric(a, T, w_o, p, oil)

    return (1.0 / hx.M) * _integrate_over_loading(integrand, W, T, w_o, p, oil)


def binding_deficit_closed_form(
    activity: float,
    T: float,
    w_o: float,
    p: GabParams,
    oil: OilIsotherm,
) -> float:
    r"""Exact partial-fraction form of the PHY-021 binding deficit.

    This is algebraically identical to :func:`binding_deficit`, whose
    loading-space quadrature remains the independent qualification oracle.  If
    ``x=K(T)*a`` and

    ``W_native=Xm*C*x/((1-x)*(1+(C-1)*x))``, then changing integration
    variable from loading to activity gives

    ``B = alpha*R/M_h * [dHC_R*Xm*C*x/(1+(C-1)*x) + dHK_R*W_native]``.

    The off-anchor oil term has no explicit temperature dependence and hence
    zero excess binding deficit under PHY-048.  Activity is deliberately the
    coordinate here: the signed off-anchor oil correction can make a
    loading-only inverse non-unique outside the admissible positive-storage
    branch.  This closed form is the runtime kernel; no fitted approximation
    or activity clamp is used.
    """
    if not math.isfinite(activity) or activity < 0.0:
        raise ValueError("hexane activity must be finite and non-negative")
    retained = retained_hexane(activity, T, w_o, p, oil)
    if retained < 0.0:
        raise ValueError("composite retained loading is negative at this activity/state")
    if activity == 0.0:
        return 0.0
    C = p.C(T)
    x = p.K(T) * activity
    native = gab_native(activity, T, p)
    alpha, _ = _oil_weights(w_o)
    return (
        R
        / hx.M
        * alpha
        * (
            p.dHC_R * p.Xm * C * x / (1.0 + (C - 1.0) * x)
            + p.dHK_R * native
        )
    )


def binding_deficit_from_psi(
    W: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm, dT: float = 0.05
) -> float:
    """B via B = T*(d psi/dT)_W - psi (central difference in T).

    Independent route used to gate Gibbs-Helmholtz consistency against
    :func:`binding_deficit`.
    """
    psi = psi_potential(W, T, w_o, p, oil)
    dpsi_dT = (
        psi_potential(W, T + dT, w_o, p, oil) - psi_potential(W, T - dT, w_o, p, oil)
    ) / (2.0 * dT)
    return T * dpsi_dT - psi


def retained_enthalpy(
    W: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm,
    pressure_pa: float = 101325.0,
    *,
    activity: float | None = None,
) -> float:
    """Retained-hexane enthalpy H_ret = W*h_l(T) - B, J/kg dry meal.

    h_l is the Span-Wagner liquid n-hexane enthalpy on the common component
    datum (PHY-049/PHY-044); latent/binding enter once, through states (no
    duplicate source).
    """
    h_l = hx.state_Tp(T, pressure_pa, "liquid").h_mass       # J/kg hexane
    if activity is None:
        deficit = binding_deficit(W, T, w_o, p, oil)
    else:
        expected = retained_hexane(activity, T, w_o, p, oil)
        if not math.isclose(expected, W, rel_tol=2.0e-13, abs_tol=2.0e-15):
            raise ValueError("retained loading and supplied activity are inconsistent")
        deficit = binding_deficit_closed_form(activity, T, w_o, p, oil)
    return W * h_l - deficit


def retained_internal_energy(
    W: float, T: float, w_o: float, p: GabParams, oil: OilIsotherm,
    pressure_pa: float = 101325.0,
    *,
    activity: float | None = None,
) -> float:
    """U_ret = H_ret - P*V_ret, V_ret = W/rho_l(T) (PHY-051), J/kg dry meal."""
    rho_l = hx.state_Tp(T, pressure_pa, "liquid").rho_mass   # kg/m3
    v_ret = W / rho_l                                        # m3/kg dry meal
    return retained_enthalpy(
        W,
        T,
        w_o,
        p,
        oil,
        pressure_pa,
        activity=activity,
    ) - pressure_pa * v_ret


# ---------------------------------------------------------------------------
# PHY-031 retained water (modified Luikov) — calorically liquid-like
# ---------------------------------------------------------------------------
def water_activity(W_w: float, p: LuikovParams) -> float:
    """a_w(W_w) on the qualified positive-moisture Luikov branch.

    Packet A5 explicitly does not authorize a continuation below ``W_ref``:
    the state potential is singular as ``W_w -> 0``.  A genuinely dry
    diagnostic state therefore needs a separately named, non-qualifying
    closure; it must not enter this authority as the convenient point
    ``(W_w, a_w) = (0, 0)``.
    """
    if isinstance(p, ContinuedPositiveLuikovParams):
        # Reuse the separately named analytic map; exact zero is only a chart
        # boundary. Wet-state evaluators independently require positive water.
        from dtdc_simulator.core2.particle.dry_limit_moisture_branch import dry_limit_water_activity

        return dry_limit_water_activity(W_w, LuikovParams(p.A1, p.A2, p.W_cap))
    if not math.isfinite(W_w) or W_w < p.W_ref:
        raise ValueError(
            f"retained-water loading must be finite and >= the qualified "
            f"Luikov lower bound {p.W_ref}; no low-moisture continuation is authorized"
        )
    if W_w > p.W_cap:
        raise ValueError(
            f"retained-water loading {W_w} exceeds the evidence cap {p.W_cap}; "
            "apply the explicit active set, do not clamp"
        )
    return math.exp((1.0 - p.A1 / W_w) / p.A2)


def water_retained(a_w: float, p: LuikovParams) -> float:
    """W_w(a_w) on the qualified positive-moisture Luikov branch."""
    if not math.isfinite(a_w) or not 0.0 < a_w <= 1.0:
        raise ValueError("water activity must be finite and lie in (0, 1]")
    continued_below_anchor = isinstance(p, ContinuedPositiveLuikovParams) and a_w < water_activity(
        LuikovParams().W_ref, LuikovParams(p.A1, p.A2, p.W_cap)
    )
    retained = p.A1 / (
        1.0 + p.A2 * (-math.log(a_w) if continued_below_anchor else math.log(1.0 / a_w))
    )
    lower_rounding_tolerance = 1.0e-12 * max(abs(p.W_ref), 1.0)
    if retained < p.W_ref - lower_rounding_tolerance:
        raise ValueError(
            f"water activity {a_w} implies retained loading {retained} below "
            f"the qualified Luikov lower bound {p.W_ref}; no low-moisture "
            "continuation is authorized"
        )
    # W_cap is a published ten-decimal evidence value; admit only its stated
    # decimal-rounding uncertainty.  The returned value is not clipped.
    cap_rounding_tolerance = 1.0e-10 * max(abs(p.W_cap), 1.0)
    if retained > p.W_cap + cap_rounding_tolerance:
        raise ValueError(
            f"water activity {a_w} implies retained loading {retained} above "
            f"the evidence cap {p.W_cap}; apply the explicit active set, do not clamp"
        )
    return retained


def water_q_net(p: LuikovParams) -> float:
    """Nominal net excess binding enthalpy of retained water = 0 (PHY-021.water).

    The modified-Luikov nominal has no explicit T dependence, so q_net,w = 0:
    retained water is calorically liquid-like. Vaporisation remains the full
    IAPWS-95 h_v - h_l (carried by the state difference, not a source).
    """
    return 0.0
