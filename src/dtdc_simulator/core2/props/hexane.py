r"""PHY-049 — pure n-hexane thermodynamic property authority.

Published Span-Wagner (2003) 12-term technical Helmholtz EOS:

    R. Span, W. Wagner, "Equations of State for Technical Applications. II.
    Results for Nonpolar Fluids," Int. J. Thermophysics 24 (2003) 41-109,
    DOI 10.1023/A:1022310214958.

The residual coefficients are transcribed in
``docs/GT_PS2_PURE_HEXANE_PROPERTY_AUTHORITY.md`` (the owner-approved authority
record). CoolProp/REFPROP are qualification oracles only and are NOT imported.

Property relations (p, s, u, cv, h, cp, w, and the second virial) follow
Span-Wagner (2003) Part I, Table II, verified against the primary source
(``literature_sources/Span_Wagner_Equations_of_State_for_Technical_Applications.pdf``).
The 12-term functional form matches Part I Eq. (6) (all tau exponents are j/8).

Reduced coordinates: ``delta = rho/rho_c``, ``tau = Tc/T``.

    alpha^r(delta, tau) = sum_i n_i * delta^d_i * tau^t_i * exp(-delta^l_i)

with the exponential omitted where l_i = 0. Every measurable property is a
derivative of this one potential, so p, f, rho, h/s/cp and dh_vap are mutually
consistent (PHY-049's stated purpose).

Datum
-----
Absolute enthalpy/entropy use a single common caloric datum (h_ref, s_ref at
T_ref); its numerical value is arbitrary and cancels in every closed-system
balance (PHY-021 datum-invariance rule; verified by tests).

Ideal-gas heat capacity — authority source, validated
-----------------------------------------------------
Saturation pressure, saturated densities, latent heat and fugacity depend ONLY
on the residual part ``alpha^r`` (the ideal-gas contribution cancels between
phases at equal T), so those are exact. Absolute cp/cv/h/s additionally need the
ideal-gas cp0(T); this uses the ``cp0`` correlation Span-Wagner Part II Table I
cites for n-hexane — Jaeschke & Schley (1995), Eq. (1) — analytically integrated
for h0/s0. The complete implementation (residual + ideal + entropy) reproduces
Span-Wagner Part II Table III (cp0, p, cp, h2-h1, s2-s1) to ~1e-4 (PHY-049 gate).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from dtdc_simulator.core2.exact_cache import float_bits_key
from dtdc_simulator.core2.props.state import HelmholtzResult, PhaseState

# --- Span-Wagner (2003) n-hexane constants (see the authority doc) -----------
R = 8.31451  # J/(mol K)  — the paper's gas constant
M = 0.08617536  # kg/mol
TC = 507.82  # K
RHOC = 2705.8778751  # mol/m3
PC = 3.0341e6  # Pa (Span-Wagner reducing pressure; informational)

_D = (1, 1, 1, 2, 3, 7, 2, 5, 1, 4, 3, 4)
_L = (0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3)
_T = (0.25, 1.125, 1.5, 1.375, 0.25, 0.875, 0.625, 1.75, 3.625, 3.625, 14.5, 12.0)
_N = (
    1.0553238,
    -2.6120616,
    0.76613883,
    -0.29770321,
    0.11879908,
    0.00027922861,
    0.4634759,
    0.011433197,
    -0.48256969,
    -0.093750559,
    -0.0067273247,
    -0.0051141584,
)

# --- Ideal-gas cp0(T): Jaeschke & Schley (1995), Eq. (1), n-hexane -----------
# The authority cited by Span-Wagner Part II Table I for the n-hexane alpha^0:
#   M. Jaeschke, P. Schley, "Ideal-Gas Thermodynamic Properties for Natural-Gas
#   Applications," Int. J. Thermophysics 16 (1995) 1381-1392.
# Expanded Aly-Lee form, Eq. (1):
#   cp0/R = B + C[(D/T)/sinh(D/T)]^2 + E[(F/T)/cosh(F/T)]^2 + G[(H/T)/sinh(H/T)]^2
# with the paper's R = 8.31451 J/mol/K (== R above). Table I gives n-hexane with
# seven parameters (no I/J cosh term; p. 1388). Verified against Span-Wagner
# Part II Table III (cp0, cp, h2-h1, s2-s1) to ~1e-4.
_JS_B = 4.0
_JS_C, _JS_D = 11.6977, 182.326  # sinh term
_JS_E, _JS_F = 26.8142, 859.207  # cosh term
_JS_G, _JS_H = 38.6164, 1826.59  # sinh term

# Common caloric datum: h_ref = 0, s_ref = 0 for the ideal gas at T_ref (K),
# p_ref (Pa). This is exactly the Span-Wagner Part II reference state (p. 46:
# "enthalpy and entropy were zero for the ideal gas at T0 = 298.15 K and
# p0 = 0.101325 MPa"). The numerical zero is arbitrary and cancels under
# datum-invariance (tested).
T_REF = 298.15
P_REF = 101325.0
H_REF = 0.0
S_REF = 0.0


# ---------------------------------------------------------------------------
# Residual Helmholtz energy and reduced derivatives
# ---------------------------------------------------------------------------
def residual(delta: float, tau: float) -> HelmholtzResult:
    """alpha^r(delta, tau) and its reduced first/second derivatives.

    Exact-bits cached (2026-08-31 speed mandate): the result depends only
    on (delta, tau) and the frozen Span-Wagner coefficient tables — no
    caloric datum term is baked in (unlike the water residual), so the
    key is the argument bits alone.
    """

    return _residual_cached(delta, tau, float_bits_key(delta, tau))


@lru_cache(maxsize=65536)
def _residual_cached(delta: float, tau: float, _bits: tuple[bytes, ...]) -> HelmholtzResult:
    return _residual_uncached(delta, tau)


def _residual_uncached(delta: float, tau: float) -> HelmholtzResult:
    phir = phir_d = phir_dd = phir_t = phir_tt = phir_dt = phir_dtt = phir_dttt = 0.0
    for n, d, exponent, t in zip(_N, _D, _L, _T):
        dl = delta**exponent if exponent else 0.0
        e = math.exp(-dl) if exponent else 1.0
        base = n * delta**d * tau**t * e
        phir += base
        # d/ddelta of delta^d exp(-delta^l) = delta^(d-1) exp(-delta^l) (d - l delta^l)
        f1 = d - exponent * dl
        phir_d += base / delta * f1
        phir_dd += base / delta**2 * (f1 * (f1 - 1.0) - exponent * exponent * dl)
        phir_t += base * t / tau
        phir_tt += base * t * (t - 1.0) / tau**2
        phir_dt += base / delta * f1 * t / tau
        phir_dtt += base / delta * f1 * t * (t - 1.0) / tau**2
        phir_dttt += base / delta * f1 * t * (t - 1.0) * (t - 2.0) / tau**3
    return HelmholtzResult(
        delta=delta,
        tau=tau,
        phir=phir,
        phir_d=phir_d,
        phir_dd=phir_dd,
        phir_t=phir_t,
        phir_tt=phir_tt,
        phir_dt=phir_dt,
        phir_dtt=phir_dtt,
        phir_dttt=phir_dttt,
    )


# ---------------------------------------------------------------------------
# Ideal-gas caloric integrals (analytic, from the Aly-Lee cp0)
# ---------------------------------------------------------------------------
def cp0(T: float) -> float:
    """Ideal-gas isobaric heat capacity, J/(mol K)  (Jaeschke-Schley Eq. 1)."""
    dT, fT, hT = _JS_D / T, _JS_F / T, _JS_H / T
    return R * (
        _JS_B
        + _JS_C * (dT / math.sinh(dT)) ** 2
        + _JS_E * (fT / math.cosh(fT)) ** 2
        + _JS_G * (hT / math.sinh(hT)) ** 2
    )


def _h0_indef(T: float) -> float:
    """Indefinite integral of cp0 dT, J/mol (up to the datum constant).

    int of a sinh term C[(D/T)/sinh(D/T)]^2 is C*D*coth(D/T); of a cosh term
    E[(F/T)/cosh(F/T)]^2 is -E*F*tanh(F/T).
    """
    return R * (
        _JS_B * T
        + _JS_C * _JS_D / math.tanh(_JS_D / T)
        - _JS_E * _JS_F * math.tanh(_JS_F / T)
        + _JS_G * _JS_H / math.tanh(_JS_H / T)
    )


def _s0_cp_indef(T: float) -> float:
    """Indefinite integral of (cp0/T) dT, J/(mol K) (up to the datum constant).

    Sinh term -> C[(D/T)coth(D/T) - ln sinh(D/T)];
    cosh term -> E[ln cosh(F/T) - (F/T)tanh(F/T)].
    """
    dT, fT, hT = _JS_D / T, _JS_F / T, _JS_H / T
    return R * (
        _JS_B * math.log(T)
        + _JS_C * (dT / math.tanh(dT) - math.log(math.sinh(dT)))
        + _JS_E * (math.log(math.cosh(fT)) - fT * math.tanh(fT))
        + _JS_G * (hT / math.tanh(hT) - math.log(math.sinh(hT)))
    )


def h_ideal(T: float) -> float:
    """Ideal-gas molar enthalpy on the common datum, J/mol."""
    return H_REF + (_h0_indef(T) - _h0_indef(T_REF))


def s_ideal(T: float, p: float) -> float:
    """Ideal-gas molar entropy on the common datum, J/(mol K)."""
    return S_REF + (_s0_cp_indef(T) - _s0_cp_indef(T_REF)) - R * math.log(p / P_REF)


# ---------------------------------------------------------------------------
# Pressure / density root / fugacity
# ---------------------------------------------------------------------------
def pressure(T: float, rho: float) -> float:
    """Pa, from p = rho R T (1 + delta alpha^r_delta)."""
    res = residual(rho / RHOC, TC / T)
    return rho * R * T * res.compressibility


def _dp_drho(T: float, rho: float) -> float:
    """dp/drho at fixed T, Pa/(mol/m3)."""
    delta, tau = rho / RHOC, TC / T
    res = residual(delta, tau)
    # p = rho R T (1 + delta phir_d); dp/drho = RT(1 + 2 delta phir_d + delta^2 phir_dd)
    return R * T * (1.0 + 2.0 * delta * res.phir_d + delta**2 * res.phir_dd)


def density(T: float, p: float, phase: str) -> float:
    """Molar density (mol/m3) root at (T, p) for ``phase`` in {'liquid','vapor'}.

    Newton from a phase-appropriate seed with a pressure-derivative guard.

    Converge-or-raise, the twin of the water surface (D2 ceremony, owner-ruled
    2026-09-03; record
    ``docs/GT_PS2_D2_VAPOUR_DENSITY_REFUSAL_CEREMONY_2026-09-03.md``).  Above
    the vapour spinodal ``p_sp(T)`` the requested vapour branch carries NO
    root.  The twin manifests DIFFERENTLY here than in water: the wider
    ``1.05`` guard push walks the iterate across the two-phase region well
    inside the 100-iteration cap, so in the operating envelope this surface
    never exhausted its cap - it silently returned the LIQUID root for a
    ``'vapor'`` request (measured: 126/126 of the no-root grid points, e.g.
    7843.17 mol/m3 at 275 K / 217 kPa).  Cap exhaustion does occur, but only
    at 180-210 K, below the fluid's operating range.  Both failure modes now
    raise.

    The branch test is exact, not a heuristic: for ``T < TC`` every root on
    the vapour branch lies at or below the vapour spinodal density, which is
    strictly below ``RHOC`` (measured sup 0.829 RHOC at 0.997 Tc; 0.310 RHOC
    over the operating envelope), so a converged root at ``rho >= RHOC`` is
    provably not on the requested branch.  At or above ``TC`` the branch
    label carries no meaning and the test is not applied.

    LIQUID MIRROR (D3 ceremony, owner-ruled 2026-09-03; record
    ``docs/GT_PS2_D3_PROPS_FAIL_OPEN_CEREMONY_2026-09-03.md``), the twin of
    the water surface.  Below the liquid spinodal pressure ``p_sp,liq(T)`` the
    requested LIQUID branch carries no root, the ``df <= 0`` guard walks the
    seed DOWN on the ``0.98`` factor, and the Newton converges on the VAPOUR
    root and returns it for a ``'liquid'`` request.  Measured on a 1,170-point
    grid: 125/125 of the no-root points escaped this way, all at
    ``T >= 472 K`` (e.g. 2.38 mol/m3 at 472 K / 10 kPa), which is exactly the
    hazard ``heteroazeotrope._liquid_fugacity`` documents and guards against
    with the same ``rho > RHOC`` test.  That caller's guard is now
    redundant-but-harmless; this authority no longer relies on it.

    The mirror test is exact in the direction that matters: every genuine
    liquid root lies at or above the liquid spinodal density, and that density
    is strictly ABOVE ``RHOC`` on every subcritical isotherm of this EOS
    (measured inf 1.1222 RHOC at 506 K over 160 isotherms, 2.60 RHOC at 300 K),
    so ``rho <= RHOC`` is provably not the liquid branch and the test cannot
    produce a false refusal.  As on the water side it is not claimed to be a
    complete detector — Span-Wagner carries spurious interior stable segments
    at ``rho/rho_c`` roughly [0.78, 1.43] — but no escape onto one was measured
    on either fluid's grid below 615 K.

    Arithmetic is untouched: every call that converged to a root on the
    requested branch before returns the identical double.
    """
    if phase == "vapor":
        rho = p / (R * T)  # ideal-gas seed
    elif phase == "liquid":
        rho = 8.5 * RHOC / 3.0  # ~ liquid branch seed (~ 7660 mol/m3)
    else:
        raise ValueError(f"phase must be 'liquid' or 'vapor', got {phase!r}")
    for _ in range(100):
        f = pressure(T, rho) - p
        df = _dp_drho(T, rho)
        if df <= 0.0:
            # push off a mechanically unstable guess toward the target branch
            rho *= 1.05 if phase == "vapor" else 0.98
            continue
        step = f / df
        rho_new = rho - step
        if rho_new <= 0.0:
            rho_new = 0.5 * rho
        if abs(rho_new - rho) <= 1e-10 * rho:
            if phase == "vapor" and T < TC and rho_new >= RHOC:
                raise ValueError(
                    f"no hexane vapor-branch root at T={T} K, p={p} Pa: the "
                    f"Newton left the vapour branch and converged on the "
                    f"liquid root rho={rho_new} mol/m3 (>= rho_c={RHOC}); the "
                    f"request is above the vapour spinodal"
                )
            if phase == "liquid" and T < TC and rho_new <= RHOC:
                raise ValueError(
                    f"no hexane liquid-branch root at T={T} K, p={p} Pa: the "
                    f"Newton left the liquid branch and converged on the "
                    f"vapour root rho={rho_new} mol/m3 (<= rho_c={RHOC}); the "
                    f"request is below the liquid spinodal"
                )
            return rho_new
        rho = rho_new
    raise ValueError(
        f"hexane {phase} density did not converge at T={T} K, p={p} Pa "
        f"(100 Newton iterations); the last iterate {rho} mol/m3 is not a root"
    )


def ln_fugacity(T: float, rho: float) -> float:
    """ln(f/Pa) for the pure fluid: ln f = ln(rho R T) + alpha^r + delta alpha^r_delta."""
    res = residual(rho / RHOC, TC / T)
    return math.log(rho * R * T) + res.phir + res.delta * res.phir_d


def second_virial(T: float) -> float:
    """Second thermal virial coefficient B(T), m^3/mol.

    Span-Wagner (2003) Part I, Table II: ``B(T)*rho_c = lim_{delta->0}
    alpha^r_delta``. This is the zero-density limit ``B_hh`` that PHY-053 draws
    from the frozen pure-fluid authority. Note the technical EOS does not fit
    gas-phase B (Part I, Sec. 3.1), so B(T) here is a *derived* — not
    independently fitted — quantity, exactly as PHY-053 specifies.
    """
    # as delta -> 0 only the d=1 terms survive alpha^r_delta; a tiny delta
    # evaluates that limit without special-casing the term list.
    return residual(1.0e-8, TC / T).phir_d / RHOC


def second_virial_dT(T: float) -> float:
    """Analytic dB/dT, m3/(mol K), from the same zero-density Helmholtz limit."""
    tau = TC / T
    return -(tau / T) * residual(1.0e-8, tau).phir_dt / RHOC


def second_virial_dT2(T: float) -> float:
    """Analytic d2B/dT2, m3/(mol K^2), CELL-02a (PHY-057).

    With B = phir_delta(0, tau)/rho_c and dtau/dT = -tau/T,
    d2B/dT2 = (tau/T^2) * (2 phir_dtau + tau phir_dtautau) / rho_c.
    """
    tau = TC / T
    res = residual(1.0e-8, tau)
    return (tau / T**2) * (2.0 * res.phir_dt + tau * res.phir_dtt) / RHOC


def second_virial_dT3(T: float) -> float:
    """Analytic d3B/dT3, m3/(mol K^3), CELL-02a2 (PHY-057 amendment)."""
    tau = TC / T
    res = residual(1.0e-8, tau)
    return (
        -(tau / T**3)
        * (6.0 * res.phir_dt + 6.0 * tau * res.phir_dtt + tau**2 * res.phir_dttt)
        / RHOC
    )


def cp0_dT(T: float) -> float:
    """Analytic dcp0/dT, J/(mol K^2), from the Jaeschke-Schley hyperbolic form."""
    total = 0.0
    for coeff, scale, kind in (
        (_JS_C, _JS_D, "sinh"),
        (_JS_E, _JS_F, "cosh"),
        (_JS_G, _JS_H, "sinh"),
    ):
        x = scale / T
        if kind == "sinh":
            total += (2.0 * coeff * x**2 / T) * (
                (x * math.cosh(x) - math.sinh(x)) / math.sinh(x) ** 3
            )
        else:
            total += (2.0 * coeff * x**2 / T) * (
                (x * math.sinh(x) - math.cosh(x)) / math.cosh(x) ** 3
            )
    return R * total


def _h_molar_from_residual(T: float, rho: float, res: "HelmholtzResult") -> float:
    """h(T, rho) = RT [1 + tau(phi0_t + phir_t) + delta phir_d] without ideal part.

    The saturation latent-heat slope only ever needs vapor-minus-liquid
    differences at one (T), so the ideal-gas tau(phi0_t) contribution cancels
    exactly and is deliberately omitted from both branches here.
    """
    return R * T * (1.0 + res.tau * res.phir_t + res.delta * res.phir_d)


def dh_vap_dT(T: float) -> float:
    """Analytic d(dh_vap)/dT along saturation, J/(mol K), CELL-02a (PHY-057).

    Total derivative of h_vap(T, rho_v(T)) - h_liq(T, rho_l(T)) along the
    saturation locus: per branch dh/dT = dh/dT|rho + dh/drho|T * drho/dT with
    drho/dT = (dpsat/dT - dp/dT|rho) / (dp/drho|T) and the exact Clapeyron
    slope dpsat/dT = dh_vap / (T (1/rho_v - 1/rho_l)).  Every piece is an
    analytic Helmholtz-derivative expression; the ideal-gas caloric part
    cancels between the branches at equal T.
    """
    sat = saturation(T)
    if not sat.converged:
        raise ValueError(f"hexane saturation did not converge at T={T}")
    dpsat_dT = sat.dh_vap / (T * (1.0 / sat.rho_vapor - 1.0 / sat.rho_liquid))
    slope = 0.0
    for rho, sign in ((sat.rho_vapor, 1.0), (sat.rho_liquid, -1.0)):
        delta, tau = rho / RHOC, TC / T
        res = residual(delta, tau)
        # dp/dT|rho and dp/drho|T from the reduced derivatives.
        dp_dT = rho * R * (1.0 + delta * res.phir_d - delta * tau * res.phir_dt)
        dp_drho = _dp_drho(T, rho)
        drho_dT = (dpsat_dT - dp_dT) / dp_drho
        # dh/dT|rho, residual part only (ideal part cancels between branches):
        # d/dT [RT(1 + tau phir_t + delta phir_d)] with tau phir_t = (TC/T) phir_t
        # and T * tau = TC constant:
        #   R(1 + delta phir_d) + R TC d(phir_t)/dT + R T delta d(phir_d)/dT
        # where d(phir_t)/dT|rho = phir_tt * (-tau/T)
        #   and d(phir_d)/dT|rho = phir_dt * (-tau/T).
        dh_dT_rho = R * (
            1.0 + delta * res.phir_d - tau * tau * res.phir_tt - delta * tau * res.phir_dt
        )
        # dh/drho|T = (RT/rho_c)(tau phir_dt + phir_d + delta phir_dd).
        dh_drho = (R * T / RHOC) * (tau * res.phir_dt + res.phir_d + delta * res.phir_dd)
        slope += sign * (dh_dT_rho + dh_drho * drho_dT)
    return slope


# ---------------------------------------------------------------------------
# Saturation (equal pressure + equal fugacity between the two density roots)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Saturation:
    T: float  # K
    p: float  # Pa
    rho_liquid: float  # mol/m3
    rho_vapor: float  # mol/m3
    dh_vap: float  # J/mol  (h_vapor - h_liquid at saturation)
    converged: bool


def _psat_seed(T: float) -> float:
    """Rough Psat(T) seed (Pa) from a hexane Antoine fit; only a Newton seed."""
    # NIST hexane Antoine (bar, K), valid ~286-343 K but adequate as a seed.
    A, B, C = 4.00266, 1171.53, -48.784
    return 10.0 ** (A - B / (T + C)) * 1.0e5


def saturation(T: float) -> Saturation:
    """Coexistence at temperature ``T`` (K), T < TC."""
    if T >= TC:
        raise ValueError(f"no saturation above Tc={TC} K (got {T})")
    p = _psat_seed(T)
    rho_l = density(T, p, "liquid")
    rho_v = density(T, p, "vapor")
    converged = False
    for _ in range(200):
        # residuals: equal pressure and equal fugacity (2x2 Newton on rho_l, rho_v)
        pl, pv = pressure(T, rho_l), pressure(T, rho_v)
        lfl, lfv = ln_fugacity(T, rho_l), ln_fugacity(T, rho_v)
        r1 = pl - pv
        r2 = lfl - lfv
        # D1 phase-2 fix, option (i-a) - the water twin's identical structural
        # edit.  Owner-ruled 2026-09-02 (the Class-B ruling packet G-0
        # demanded); measured in
        # docs/GT_PS2_D1_WATER_SATURATION_CHARACTERIZATION_2026-09-02.md
        # (sections 3 and 5) against
        # docs/GT_PS2_G0_WATER_SATURATION_CONVERGENCE_RECORD_2026-08-23.md.
        # `p` is now assigned from the CURRENT iterate unconditionally, so a
        # non-accepted exit (200-iteration cap, or the det == 0 break below)
        # returns the settled branch-pressure mean instead of the Antoine
        # seed; `converged` becomes pure metadata.  Bit-neutral on every
        # accepted exit by construction (identical pl, pv, identical doubles).
        # Tolerances, the iteration cap and the flag semantics are UNTOUCHED.
        p = 0.5 * (pl + pv)
        if abs(r1) < 1e-11 * max(pl, 1.0) and abs(r2) < 1e-13:
            converged = True
            break
        dpl, dpv = _dp_drho(T, rho_l), _dp_drho(T, rho_v)
        # d ln f/drho = (dp/drho)/(rho R T)   [exact for a pure fluid]
        dlfl = dpl / (rho_l * R * T)
        dlfv = dpv / (rho_v * R * T)
        # J = [[dpl, -dpv], [dlfl, -dlfv]]
        det = dpl * (-dlfv) - (-dpv) * dlfl
        if det == 0.0:
            break
        # solve J [drl; drv] = [r1; r2]  (Newton: x -= J^-1 r)
        drl = (r1 * (-dlfv) - (-dpv) * r2) / det
        drv = (dpl * r2 - dlfl * r1) / det
        rho_l -= drl
        rho_v -= drv
        if rho_l <= 0.0:
            rho_l = 0.5 * (rho_l + drl)
        if rho_v <= 0.0:
            rho_v = 0.5 * (rho_v + drv)
    hl = _residual_h(T, rho_l)
    hv = _residual_h(T, rho_v)
    return Saturation(
        T=T,
        p=p,
        rho_liquid=rho_l,
        rho_vapor=rho_v,
        dh_vap=hv - hl,
        converged=converged,
    )


def _residual_h(T: float, rho: float) -> float:
    """Residual molar enthalpy h - h_ideal, J/mol."""
    res = residual(rho / RHOC, TC / T)
    return R * T * (res.tau * res.phir_t + res.delta * res.phir_d)


def saturation_pressure(T: float) -> float:
    """Psat(T), Pa."""
    return saturation(T).p


def dh_vap(T: float) -> float:
    """Latent heat of vaporisation at temperature T, J/mol (datum-free)."""
    return saturation(T).dh_vap


# ---------------------------------------------------------------------------
# Full single-phase state
# ---------------------------------------------------------------------------
def state(T: float, rho: float) -> PhaseState:
    """Full molar thermodynamic state at (T, rho)."""
    delta, tau = rho / RHOC, TC / T
    res = residual(delta, tau)
    Z = res.compressibility
    p = rho * R * T * Z
    h_res = R * T * (tau * res.phir_t + delta * res.phir_d)
    s_res = R * (tau * res.phir_t - res.phir)
    cv_res = -R * tau**2 * res.phir_tt
    h = h_ideal(T) + h_res
    u = h - p / rho  # u = h - p v (molar)
    # Entropy departure s_res is defined at constant density, so the ideal-gas
    # part must be evaluated at the ideal-gas pressure rho*R*T (NOT the real p);
    # the two differ by R*ln(Z). Verified against Span-Wagner Part II Table III.
    s = s_ideal(T, rho * R * T) + s_res
    cv = (cp0(T) - R) + cv_res
    # cp from the standard Helmholtz relation
    num = (1.0 + delta * res.phir_d - delta * tau * res.phir_dt) ** 2
    den = 1.0 + 2.0 * delta * res.phir_d + delta**2 * res.phir_dd
    cp = cv + R * num / den
    lnf = math.log(rho * R * T) + res.phir + delta * res.phir_d
    return PhaseState(
        T=T,
        rho=rho,
        p=p,
        Z=Z,
        h=h,
        u=u,
        s=s,
        cv=cv,
        cp=cp,
        ln_fugacity=lnf,
        molar_mass=M,
    )


_STATE_TP_CACHE_MAXSIZE = 4096


def caloric_datum_signature() -> tuple[bytes, ...]:
    """Exact cache signature for the mutable common caloric datum constants."""

    return float_bits_key(H_REF, S_REF, T_REF, P_REF)


def state_Tp(T: float, p: float, phase: str) -> PhaseState:
    """Full state at (T, p), memoized only for bit-exact repeated inputs."""

    return _state_Tp_cached(
        T,
        p,
        phase,
        float_bits_key(T, p),
        caloric_datum_signature(),
    )


@lru_cache(maxsize=_STATE_TP_CACHE_MAXSIZE, typed=True)
def _state_Tp_cached(
    T: float,
    p: float,
    phase: str,
    _scalar_key: tuple[bytes, ...],
    _datum_key: tuple[bytes, ...],
) -> PhaseState:
    # The exact keys deliberately participate only in memoization. Evaluation
    # still uses the original values and the unmodified Helmholtz equations.
    return state(T, density(T, p, phase))


def clear_state_tp_cache() -> None:
    """Drop every cached pure-hexane ``(T,p,phase,datum)`` state."""

    _state_Tp_cached.cache_clear()


def state_tp_cache_info():
    """Return bounded-cache hit/miss/size diagnostics."""

    return _state_Tp_cached.cache_info()
