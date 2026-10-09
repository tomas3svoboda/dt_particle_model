r"""PHY-030 — pure water thermodynamic property authority (IAPWS-95).

Reference Helmholtz formulation of Wagner & Pruss (2002), released as IAPWS
R6-95. Implemented self-contained; IAPWS-IF97 / CoolProp / steam tables are
qualification oracles only (the frozen decision makes IAPWS-95 the authority
because IF97 is a piecewise industrial fit awkward for the pure second virial
that PHY-053 needs).

Per-mass (specific) formulation, SI. Dimensionless Helmholtz energy

    a(rho, T)/(R T) = phi0(delta, tau) + phir(delta, tau),
    delta = rho/rho_c,  tau = Tc/T.

Residual part has 56 terms: 7 polynomial, 44 polynomial*exp(-delta^c), 3
Gaussian bell-shaped, 2 non-analytic (critical-region) terms — all included so
the module reproduces the official IAPWS-95 verification table, including the
near-critical point, not only the DT operating range.

One common water component datum (PHY-021 water rule): the IAPWS-95 ideal-gas
constants n1/n2 fix the triple-point reference; its numerical zero is arbitrary
and cancels under datum shifts (tested).

Property relations (p, s, u, cv, h, cp, w, and the second virial) follow the
same Helmholtz-derivative table as the hexane authority — Span-Wagner (2003)
Part I, Table II — which is the general reduced-Helmholtz property table, not
fluid-specific.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from dtdc_simulator.core2.exact_cache import float_bits_key
from dtdc_simulator.core2.props.state import HelmholtzResult, MassPhaseState

# --- IAPWS-95 reducing constants ---------------------------------------------
R = 461.51805  # J/(kg K)  — specific gas constant of water
M = 0.018015268  # kg/mol
TC = 647.096  # K
RHOC = 322.0  # kg/m3
PC = 22.064e6  # Pa
TT = 273.16  # K, triple point
PT = 611.655  # Pa, triple-point pressure

# --- Ideal-gas part (Table 1) ------------------------------------------------
_N0 = (
    -8.3204464837497,
    6.6832105275932,
    3.00632,
    0.012436,
    0.97315,
    1.27950,
    0.96956,
    0.24873,
)
_G0 = (0.0, 0.0, 0.0, 1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105)

# --- Residual part (Table 2) -------------------------------------------------
# Group 1: i=1..7  (no exponential, c=0)
# Group 2: i=8..51 (exp(-delta^c))
_C = (
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    4,
    6,
    6,
    6,
    6,
)
_D = (
    1,
    1,
    1,
    2,
    2,
    3,
    4,
    1,
    1,
    1,
    2,
    2,
    3,
    4,
    4,
    5,
    7,
    9,
    10,
    11,
    13,
    15,
    1,
    2,
    2,
    2,
    3,
    4,
    4,
    4,
    5,
    6,
    6,
    7,
    9,
    9,
    9,
    9,
    9,
    10,
    10,
    12,
    3,
    4,
    4,
    5,
    14,
    3,
    6,
    6,
    6,
)
_T = (
    -0.5,
    0.875,
    1.0,
    0.5,
    0.75,
    0.375,
    1.0,
    4.0,
    6.0,
    12.0,
    1.0,
    5.0,
    4.0,
    2.0,
    13.0,
    9.0,
    3.0,
    4.0,
    11.0,
    4.0,
    13.0,
    1.0,
    7.0,
    1.0,
    9.0,
    10.0,
    10.0,
    3.0,
    7.0,
    10.0,
    10.0,
    6.0,
    10.0,
    10.0,
    1.0,
    2.0,
    3.0,
    4.0,
    8.0,
    6.0,
    9.0,
    8.0,
    16.0,
    22.0,
    23.0,
    23.0,
    10.0,
    50.0,
    44.0,
    46.0,
    50.0,
)
_NR = (
    0.12533547935523e-1,
    0.78957634722828e1,
    -0.87803203303561e1,
    0.31802509345418,
    -0.26145533859358,
    -0.78199751687981e-2,
    0.88089493102134e-2,
    -0.66856572307965,
    0.20433810950965,
    -0.66212605039687e-4,
    -0.19232721156002,
    -0.25709043003438,
    0.16074868486251,
    -0.40092828925807e-1,
    0.39343422603254e-6,
    -0.75941377088144e-5,
    0.56250979351888e-3,
    -0.15608652257135e-4,
    0.11537996422951e-8,
    0.36582165144204e-6,
    -0.13251180074668e-11,
    -0.62639586912454e-9,
    -0.10793600908932,
    0.17611491008752e-1,
    0.22132295167546,
    -0.40247669763528,
    0.58083399985759,
    0.49969146990806e-2,
    -0.31358700712549e-1,
    -0.74315929710341,
    0.47807329915480,
    0.20527940895948e-1,
    -0.13636435110343,
    0.14180634400617e-1,
    0.83326504880713e-2,
    -0.29052336009585e-1,
    0.38615085574206e-1,
    -0.20393486513704e-1,
    -0.16554050063734e-2,
    0.19955571979541e-2,
    0.15870308324157e-3,
    -0.16388568342530e-4,
    0.43613615723811e-1,
    0.34994005463765e-1,
    -0.76788197844621e-1,
    0.22446277332006e-1,
    -0.62689710414685e-4,
    -0.55711118565645e-9,
    -0.19905718354408,
    0.31777497330738,
    -0.11841182425981,
)

# Group 3: i=52..54  Gaussian bell-shaped
_G_D = (3, 3, 3)
_G_T = (0.0, 1.0, 4.0)
_G_N = (-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4)
_G_ALPHA = (20.0, 20.0, 20.0)
_G_BETA = (150.0, 150.0, 250.0)
_G_GAMMA = (1.21, 1.21, 1.25)
_G_EPS = (1.0, 1.0, 1.0)

# Group 4: i=55..56  non-analytic
_NA_A = (3.5, 3.5)
_NA_B = (0.85, 0.95)
_NA_BB = (0.2, 0.2)  # capital B_i
_NA_N = (-0.14874640856724, 0.31806110878444)
_NA_C = (28.0, 32.0)
_NA_D = (700.0, 800.0)
_NA_AA = (0.32, 0.32)  # capital A_i
_NA_BETA = (0.3, 0.3)

# Common datum reference state for reporting (IAPWS uses saturated-liquid triple
# point h=s=0; the n1/n2 constants above encode a fixed reference — its numerical
# zero is arbitrary for closed-system balances, verified by the datum test).


# ---------------------------------------------------------------------------
# Ideal-gas dimensionless Helmholtz energy and derivatives
# ---------------------------------------------------------------------------
def _phi0(delta: float, tau: float) -> tuple[float, float, float]:
    """Return (phi0, phi0_tau, phi0_tautau)."""
    phi = math.log(delta) + _N0[0] + _N0[1] * tau + _N0[2] * math.log(tau)
    phi_t = _N0[1] + _N0[2] / tau
    phi_tt = -_N0[2] / tau**2
    for i in range(3, 8):
        e = math.exp(-_G0[i] * tau)
        one_minus = 1.0 - e
        phi += _N0[i] * math.log(one_minus)
        phi_t += _N0[i] * _G0[i] * e / one_minus
        phi_tt += -_N0[i] * _G0[i] ** 2 * e / one_minus**2
    return phi, phi_t, phi_tt


# ---------------------------------------------------------------------------
# Residual dimensionless Helmholtz energy and derivatives
# ---------------------------------------------------------------------------
def residual(delta: float, tau: float) -> HelmholtzResult:
    """alpha^r(delta, tau) and its reduced first/second derivatives (56 terms).

    Exact-bit memoized (GT_PS2_EXACT_MEMOIZATION_RULING 2026-08-30 standing
    mandate; measured 94.0 % exact-argument repeats per march interval).
    The returned object bakes in the ideal ``phi0*`` terms from ``_N0``,
    so the key carries the same caloric datum signature the ``state_Tp``
    cache already uses.
    """

    return _residual_cached(delta, tau, float_bits_key(delta, tau), caloric_datum_signature())


_RESIDUAL_CACHE_MAXSIZE = 65536


@lru_cache(maxsize=_RESIDUAL_CACHE_MAXSIZE, typed=True)
def _residual_cached(
    delta: float,
    tau: float,
    _scalar_key: tuple[bytes, ...],
    _datum_key: tuple[bytes, ...],
) -> HelmholtzResult:
    # The original values, not decoded key bytes, enter the evaluation.
    return _residual_uncached(delta, tau)


def _residual_uncached(delta: float, tau: float) -> HelmholtzResult:
    """Uncached 56-term evaluation (equivalence-test surface)."""
    phir = phir_d = phir_dd = phir_t = phir_tt = phir_dt = phir_dtt = phir_dttt = 0.0

    # Groups 1 & 2: polynomial and polynomial * exp(-delta^c)
    for n, c, d, t in zip(_NR, _C, _D, _T):
        if c == 0:
            e = 1.0
            f1 = d
            f2 = d * (d - 1.0)
        else:
            dc = delta**c
            e = math.exp(-dc)
            f1 = d - c * dc
            f2 = f1 * (f1 - 1.0) - c * c * dc
        base = n * delta**d * tau**t * e
        phir += base
        phir_d += base / delta * f1
        phir_dd += base / delta**2 * f2
        phir_t += base * t / tau
        phir_tt += base * t * (t - 1.0) / tau**2
        phir_dt += base / delta * f1 * t / tau
        phir_dtt += base / delta * f1 * t * (t - 1.0) / tau**2
        phir_dttt += base / delta * f1 * t * (t - 1.0) * (t - 2.0) / tau**3

    # Group 3: Gaussian bell-shaped terms
    for n, d, t, al, be, ga, ep in zip(_G_N, _G_D, _G_T, _G_ALPHA, _G_BETA, _G_GAMMA, _G_EPS):
        g = math.exp(-al * (delta - ep) ** 2 - be * (tau - ga) ** 2)
        base = n * delta**d * tau**t * g
        gd = d / delta - 2.0 * al * (delta - ep)
        gt = t / tau - 2.0 * be * (tau - ga)
        phir += base
        phir_d += base * gd
        phir_dd += base * (gd * gd - d / delta**2 - 2.0 * al)
        phir_t += base * gt
        phir_tt += base * (gt * gt - t / tau**2 - 2.0 * be)
        phir_dt += base * gd * gt
        phir_dtt += base * gd * (gt * gt - t / tau**2 - 2.0 * be)
        gt_prime = -t / tau**2 - 2.0 * be
        phir_dttt += base * gd * (gt**3 + 3.0 * gt * gt_prime + 2.0 * t / tau**3)

    # Group 4: non-analytic terms
    for n, a, b, BB, C, D, AA, be in zip(
        _NA_N, _NA_A, _NA_B, _NA_BB, _NA_C, _NA_D, _NA_AA, _NA_BETA
    ):
        dm1 = delta - 1.0
        tm1 = tau - 1.0
        # guard the exact ridge delta==1 (never at verification points)
        if dm1 == 0.0:
            dm1 = 1.0e-12
        dm1sq = dm1 * dm1
        theta = (1.0 - tau) + AA * dm1sq ** (1.0 / (2.0 * be))
        Delta = theta * theta + BB * dm1sq**a
        psi = math.exp(-C * dm1sq - D * tm1 * tm1)

        # psi derivatives
        psi_d = -2.0 * C * dm1 * psi
        psi_dd = (2.0 * C * dm1sq - 1.0) * 2.0 * C * psi
        psi_t = -2.0 * D * tm1 * psi
        psi_tt = (2.0 * D * tm1 * tm1 - 1.0) * 2.0 * D * psi
        psi_dt = 4.0 * C * D * dm1 * tm1 * psi

        # Delta derivatives
        dDelta_d = dm1 * (
            AA * theta * (2.0 / be) * dm1sq ** (1.0 / (2.0 * be) - 1.0)
            + 2.0 * BB * a * dm1sq ** (a - 1.0)
        )
        d2Delta_dd = dDelta_d / dm1 + dm1sq * (
            4.0 * BB * a * (a - 1.0) * dm1sq ** (a - 2.0)
            + 2.0 * AA * AA * (1.0 / be) ** 2 * (dm1sq ** (1.0 / (2.0 * be) - 1.0)) ** 2
            + AA * theta * (4.0 / be) * (1.0 / (2.0 * be) - 1.0) * dm1sq ** (1.0 / (2.0 * be) - 2.0)
        )

        # Delta^b derivatives
        Db = Delta**b
        dDb_d = b * Delta ** (b - 1.0) * dDelta_d
        d2Db_dd = b * (
            Delta ** (b - 1.0) * d2Delta_dd + (b - 1.0) * Delta ** (b - 2.0) * dDelta_d * dDelta_d
        )
        dDb_t = -2.0 * theta * b * Delta ** (b - 1.0)
        d2Db_tt = 2.0 * b * Delta ** (b - 1.0) + 4.0 * theta * theta * b * (b - 1.0) * Delta ** (
            b - 2.0
        )
        d2Db_dt = (
            -AA * b * (2.0 / be) * Delta ** (b - 1.0) * dm1 * dm1sq ** (1.0 / (2.0 * be) - 1.0)
            - 2.0 * theta * b * (b - 1.0) * Delta ** (b - 2.0) * dDelta_d
        )

        phir += n * Db * delta * psi
        phir_d += n * (Db * (psi + delta * psi_d) + dDb_d * delta * psi)
        phir_dd += n * (
            Db * (2.0 * psi_d + delta * psi_dd)
            + 2.0 * dDb_d * (psi + delta * psi_d)
            + d2Db_dd * delta * psi
        )
        phir_t += n * delta * (dDb_t * psi + Db * psi_t)
        phir_tt += n * delta * (d2Db_tt * psi + 2.0 * dDb_t * psi_t + Db * psi_tt)
        phir_dt += n * (
            Db * (psi_t + delta * psi_dt)
            + delta * dDb_d * psi_t
            + dDb_t * (psi + delta * psi_d)
            + d2Db_dt * delta * psi
        )
        # CELL-02a: third mixed derivative d3/(ddelta dtau^2) of n Db delta psi.
        psi_dtt = 4.0 * C * D * dm1 * psi * (1.0 - 2.0 * D * tm1 * tm1)
        P_geom = dm1 * dm1sq ** (1.0 / (2.0 * be) - 1.0)
        d3Db_dtt = (
            4.0 * AA * b * (b - 1.0) * (2.0 / be) * theta * Delta ** (b - 2.0) * P_geom
            + 2.0 * b * (b - 1.0) * Delta ** (b - 2.0) * dDelta_d
            + 4.0 * b * (b - 1.0) * (b - 2.0) * theta * theta * Delta ** (b - 3.0) * dDelta_d
        )
        phir_dtt += n * (
            d2Db_tt * (psi + delta * psi_d)
            + 2.0 * dDb_t * (psi_t + delta * psi_dt)
            + Db * (psi_tt + delta * psi_dtt)
            + d3Db_dtt * delta * psi
            + 2.0 * d2Db_dt * delta * psi_t
        )
        # CELL-02a2: fourth mixed derivative d4/(ddelta dtau^3).
        psi_ttt = 4.0 * D * D * tm1 * psi * (3.0 - 2.0 * D * tm1 * tm1)
        psi_dttt = -8.0 * C * D * D * dm1 * tm1 * psi * (3.0 - 2.0 * D * tm1 * tm1)
        d3Db_ttt = -12.0 * theta * b * (b - 1.0) * Delta ** (b - 2.0) - 8.0 * theta**3 * b * (
            b - 1.0
        ) * (b - 2.0) * Delta ** (b - 3.0)
        Q_geom = AA * (2.0 / be)
        A1 = 4.0 * AA * b * (b - 1.0) * (2.0 / be)
        A2 = 2.0 * b * (b - 1.0)
        A3 = 4.0 * b * (b - 1.0) * (b - 2.0)
        d4Db_dttt = (
            A1 * P_geom * (-(Delta ** (b - 2.0)) - 2.0 * theta**2 * (b - 2.0) * Delta ** (b - 3.0))
            + A2
            * (
                -2.0 * theta * (b - 2.0) * Delta ** (b - 3.0) * dDelta_d
                - Q_geom * P_geom * Delta ** (b - 2.0)
            )
            + A3
            * (
                -2.0 * theta * Delta ** (b - 3.0) * dDelta_d
                - 2.0 * theta**3 * (b - 3.0) * Delta ** (b - 4.0) * dDelta_d
                - Q_geom * P_geom * theta**2 * Delta ** (b - 3.0)
            )
        )
        phir_dttt += n * (
            d3Db_ttt * (psi + delta * psi_d)
            + 3.0 * d2Db_tt * (psi_t + delta * psi_dt)
            + 3.0 * dDb_t * (psi_tt + delta * psi_dtt)
            + Db * (psi_ttt + delta * psi_dttt)
            + d4Db_dttt * delta * psi
            + 3.0 * d3Db_dtt * delta * psi_t
            + 3.0 * d2Db_dt * delta * psi_tt
        )

    phi0, phi0_t, phi0_tt = _phi0(delta, tau)
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
        phi0=phi0,
        phi0_t=phi0_t,
        phi0_tt=phi0_tt,
    )


# ---------------------------------------------------------------------------
# Pressure / density / fugacity
# ---------------------------------------------------------------------------
def pressure(T: float, rho: float) -> float:
    """Pa, from p = rho R T (1 + delta phir_delta)."""
    res = residual(rho / RHOC, TC / T)
    return rho * R * T * res.compressibility


def _dp_drho(T: float, rho: float) -> float:
    delta, tau = rho / RHOC, TC / T
    res = residual(delta, tau)
    return R * T * (1.0 + 2.0 * delta * res.phir_d + delta**2 * res.phir_dd)


def density(T: float, p: float, phase: str) -> float:
    """Molar-free specific density (kg/m3) root at (T, p) for the chosen branch.

    Exact-bit memoized (measured 83.2 % repeats); residual-only, datum-free.
    """

    return _density_cached(T, p, phase, float_bits_key(T, p))


@lru_cache(maxsize=16384, typed=True)
def _density_cached(T: float, p: float, phase: str, _scalar_key: tuple[bytes, ...]) -> float:
    return _density_uncached(T, p, phase)


def _density_uncached(T: float, p: float, phase: str) -> float:
    """Uncached density root (equivalence-test surface).

    Converge-or-raise (D2 ceremony, owner-ruled 2026-09-03; record
    ``docs/GT_PS2_D2_VAPOUR_DENSITY_REFUSAL_CEREMONY_2026-09-03.md``).  Above
    the vapour spinodal ``p_sp(T)`` the requested vapour branch carries NO
    root, and the pre-D2 loop answered such a request with a number anyway:
    either the wandering iterate at the 200-iteration cap (measured 98.04
    kg/m3 at 290 K / 250 kPa - neither the vapour root, which does not exist
    there, nor the liquid one at 998.9) or, when the ``df <= 0`` guard walked
    the iterate across the two-phase region in under 200 steps, the LIQUID
    root returned for a ``'vapor'`` request.  Both are silent wrong answers
    and both now raise.

    The branch test is exact, not a heuristic: for ``T < TC`` every root on
    the vapour branch lies at or below the vapour spinodal density, which is
    strictly below ``RHOC`` (measured sup 0.777 RHOC at 0.997 Tc; 0.054 RHOC
    over the operating envelope), so a converged root at ``rho >= RHOC`` is
    provably not on the requested branch.  At or above ``TC`` the branch
    label carries no meaning and the test is not applied.

    LIQUID MIRROR (D3 ceremony, owner-ruled 2026-09-03; record
    ``docs/GT_PS2_D3_PROPS_FAIL_OPEN_CEREMONY_2026-09-03.md``).  D2 left the
    opposite escape open and reported it: below the liquid spinodal pressure
    ``p_sp,liq(T)`` the requested LIQUID branch carries no root either, the
    ``df <= 0`` guard then walks the seed DOWN across the two-phase region on
    the ``0.99`` factor, and the Newton converges on the VAPOUR root, which is
    returned for a ``'liquid'`` request.  Measured: ``_density_uncached(600.0,
    1000.0, 'liquid')`` returned 3.611e-03 kg/m3 — a vapour root, labelled
    liquid, at a machine-precision pressure residual, so no downstream check
    on the residual can catch it.  It now raises.

    The mirror test is exact in the direction that matters: every genuine
    liquid root lies at or above the liquid spinodal density, and that density
    is strictly ABOVE ``RHOC`` on every subcritical isotherm of this EOS
    (measured inf 1.2389 RHOC at 645 K over 160 isotherms, 2.85 RHOC at 275 K),
    so ``rho <= RHOC`` is provably not the liquid branch and the test cannot
    produce a false refusal.  It is deliberately NOT claimed to be a complete
    detector: IAPWS-95 carries spurious interior stable segments in
    ``rho/rho_c`` roughly [0.87, 1.24] deep inside the two-phase region, and a
    single measured escape (615 K, 1.487 MPa) converged onto one at 1.043 RHOC
    and is therefore not caught.  That corner is 145 K above the operating
    envelope; see the ceremony record §7.

    Arithmetic is untouched: every call that converged to a root on the
    requested branch before returns the identical double.
    """
    if phase == "vapor":
        rho = p / (R * T)
    elif phase == "liquid":
        rho = 3.1 * RHOC  # ~ 1000 kg/m3 liquid seed
    else:
        raise ValueError(f"phase must be 'liquid' or 'vapor', got {phase!r}")
    for _ in range(200):
        f = pressure(T, rho) - p
        df = _dp_drho(T, rho)
        if df <= 0.0:
            rho *= 1.02 if phase == "vapor" else 0.99
            continue
        rho_new = rho - f / df
        if rho_new <= 0.0:
            rho_new = 0.5 * rho
        if abs(rho_new - rho) <= 1.0e-9 * rho:
            if phase == "vapor" and T < TC and rho_new >= RHOC:
                raise ValueError(
                    f"no water vapor-branch root at T={T} K, p={p} Pa: the "
                    f"Newton left the vapour branch and converged on the "
                    f"liquid root rho={rho_new} kg/m3 (>= rho_c={RHOC}); the "
                    f"request is above the vapour spinodal"
                )
            if phase == "liquid" and T < TC and rho_new <= RHOC:
                raise ValueError(
                    f"no water liquid-branch root at T={T} K, p={p} Pa: the "
                    f"Newton left the liquid branch and converged on the "
                    f"vapour root rho={rho_new} kg/m3 (<= rho_c={RHOC}); the "
                    f"request is below the liquid spinodal"
                )
            return rho_new
        rho = rho_new
    raise ValueError(
        f"water {phase} density did not converge at T={T} K, p={p} Pa "
        f"(200 Newton iterations); the last iterate {rho} kg/m3 is not a root"
    )


def ln_fugacity(T: float, rho: float) -> float:
    """ln(f) with f in Pa; pure fluid ln f = ln(rho R T) + phir + delta phir_delta."""
    res = residual(rho / RHOC, TC / T)
    return math.log(rho * R * T) + res.phir + res.delta * res.phir_d


def second_virial(T: float) -> float:
    """Second thermal virial coefficient B(T), m^3/mol.

    Span-Wagner (2003) Part I, Table II: ``B(T)*rho_c = lim_{delta->0}
    phir_delta``. IAPWS-95 fits the gas region well, so this is the ``B_ww``
    PHY-053 draws from the frozen pure-water authority. ``rho_c`` is converted
    to a molar basis (IAPWS-95 is specific/per-mass).
    """
    rho_c_molar = RHOC / M  # mol/m3
    return residual(1.0e-8, TC / T).phir_d / rho_c_molar


def second_virial_dT(T: float) -> float:
    """Analytic dB/dT, m3/(mol K), from the IAPWS-95 zero-density limit."""
    tau = TC / T
    rho_c_molar = RHOC / M
    return -(tau / T) * residual(1.0e-8, tau).phir_dt / rho_c_molar


def second_virial_dT2(T: float) -> float:
    """Analytic d2B/dT2, m3/(mol K^2), CELL-02a (PHY-057)."""
    tau = TC / T
    rho_c_molar = RHOC / M
    res = residual(1.0e-8, tau)
    return (tau / T**2) * (2.0 * res.phir_dt + tau * res.phir_dtt) / rho_c_molar


def second_virial_dT3(T: float) -> float:
    """Analytic d3B/dT3, m3/(mol K^3), CELL-02a2 (PHY-057 amendment)."""
    tau = TC / T
    rho_c_molar = RHOC / M
    res = residual(1.0e-8, tau)
    return (
        -(tau / T**3)
        * (6.0 * res.phir_dt + 6.0 * tau * res.phir_dtt + tau**2 * res.phir_dttt)
        / rho_c_molar
    )


def phi0_ttt(tau: float) -> float:
    """Third tau-derivative of the ideal-gas part, CELL-02a2."""
    total = 2.0 * _N0[2] / tau**3
    for i in range(3, 8):
        e = math.exp(-_G0[i] * tau)
        one_minus = 1.0 - e
        total += _N0[i] * _G0[i] ** 3 * e * (1.0 + e) / one_minus**3
    return total


def dh_vap_dT(T: float) -> float:
    """Analytic d(dh_vap)/dT along saturation, J/(kg K), CELL-02a (PHY-057).

    Mirrors the hexane derivation on the mass basis: per branch
    dh/dT = dh/dT|rho + dh/drho|T * drho/dT, with the exact Clapeyron slope
    and the ideal-gas caloric part cancelling between branches at equal T.
    """
    sat = saturation(T)
    if not sat.converged:
        raise ValueError(f"water saturation did not converge at T={T}")
    dpsat_dT = sat.dh_vap / (T * (1.0 / sat.rho_vapor - 1.0 / sat.rho_liquid))
    slope = 0.0
    for rho, sign in ((sat.rho_vapor, 1.0), (sat.rho_liquid, -1.0)):
        delta, tau = rho / RHOC, TC / T
        res = residual(delta, tau)
        dp_dT = rho * R * (1.0 + delta * res.phir_d - delta * tau * res.phir_dt)
        dp_drho = _dp_drho(T, rho)
        drho_dT = (dpsat_dT - dp_dT) / dp_drho
        dh_dT_rho = R * (
            1.0 + delta * res.phir_d - tau * tau * res.phir_tt - delta * tau * res.phir_dt
        )
        dh_drho = (R * T / RHOC) * (tau * res.phir_dt + res.phir_d + delta * res.phir_dd)
        slope += sign * (dh_dT_rho + dh_drho * drho_dT)
    return slope


# ---------------------------------------------------------------------------
# Saturation
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Saturation:
    T: float
    p: float
    rho_liquid: float
    rho_vapor: float
    dh_vap: float  # J/kg
    converged: bool


def _psat_seed(T: float) -> float:
    """IAPWS-97-style saturation seed (Pa); only a Newton seed, not the authority."""
    # Wagner saturation-pressure ansatz on the IAPWS-95 critical constants.
    theta = 1.0 - T / TC
    a = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    lnpr = (
        TC
        / T
        * (
            a[0] * theta
            + a[1] * theta**1.5
            + a[2] * theta**3
            + a[3] * theta**3.5
            + a[4] * theta**4
            + a[5] * theta**7.5
        )
    )
    return PC * math.exp(lnpr)


def _residual_h(T: float, rho: float) -> float:
    """Residual specific enthalpy h - h_ideal, J/kg."""
    res = residual(rho / RHOC, TC / T)
    return R * T * (res.tau * res.phir_t + res.delta * res.phir_d)


def saturation(T: float) -> Saturation:
    """Coexistence at temperature T (K), TT <= T < TC.

    Exact-bit memoized (measured 93.2 % repeats); residual-only — the
    ideal caloric datum cancels in coexistence and in dh_vap.
    """

    return _saturation_cached(T, float_bits_key(T))


@lru_cache(maxsize=8192, typed=True)
def _saturation_cached(T: float, _scalar_key: tuple[bytes, ...]) -> Saturation:
    return _saturation_uncached(T)


def _saturation_uncached(T: float) -> Saturation:
    """Uncached coexistence solve (equivalence-test surface)."""
    if T >= TC:
        raise ValueError(f"no saturation above Tc={TC} K (got {T})")
    p = _psat_seed(T)
    rho_l = density(T, p, "liquid")
    rho_v = density(T, p, "vapor")
    converged = False
    for _ in range(200):
        pl, pv = pressure(T, rho_l), pressure(T, rho_v)
        lfl, lfv = ln_fugacity(T, rho_l), ln_fugacity(T, rho_v)
        r1, r2 = pl - pv, lfl - lfv
        # D1 phase-2 fix, option (i-a).  Owner-ruled 2026-09-02 (the Class-B
        # ruling packet G-0 demanded); measured in
        # docs/GT_PS2_D1_WATER_SATURATION_CHARACTERIZATION_2026-09-02.md
        # (section 5) against
        # docs/GT_PS2_G0_WATER_SATURATION_CONVERGENCE_RECORD_2026-08-23.md.
        # `p` is now assigned from the CURRENT iterate unconditionally, so a
        # non-accepted exit (200-iteration cap, or the det == 0 break below)
        # returns the settled branch-pressure mean instead of the Wagner seed;
        # `converged` becomes pure metadata.  Bit-neutral on every accepted
        # exit by construction (identical pl, pv, identical doubles: measured
        # 0 changed bits over 11,974 converged points).  Tolerances, the
        # iteration cap and the flag semantics are UNTOUCHED.
        p = 0.5 * (pl + pv)
        if abs(r1) < 1e-10 * max(pl, 1.0) and abs(r2) < 1e-13:
            converged = True
            break
        dpl, dpv = _dp_drho(T, rho_l), _dp_drho(T, rho_v)
        dlfl, dlfv = dpl / (rho_l * R * T), dpv / (rho_v * R * T)
        det = dpl * (-dlfv) - (-dpv) * dlfl
        if det == 0.0:
            break
        drl = (r1 * (-dlfv) - (-dpv) * r2) / det
        drv = (dpl * r2 - dlfl * r1) / det
        rho_l = max(rho_l - drl, 0.5 * rho_l)
        rho_v = max(rho_v - drv, 0.5 * rho_v)
    hl, hv = _residual_h(T, rho_l), _residual_h(T, rho_v)
    return Saturation(
        T=T,
        p=p,
        rho_liquid=rho_l,
        rho_vapor=rho_v,
        dh_vap=hv - hl,
        converged=converged,
    )


def saturation_pressure(T: float) -> float:
    return saturation(T).p


def saturation_temperature(p: float) -> float:
    """Invert Psat(T) = p for T (K) by bracketed secant.

    Exact-bit memoized (measured 99.0 % repeats); saturation-only.
    """

    return _saturation_temperature_cached(p, float_bits_key(p))


@lru_cache(maxsize=4096, typed=True)
def _saturation_temperature_cached(p: float, _scalar_key: tuple[bytes, ...]) -> float:
    return _saturation_temperature_uncached(p)


def _saturation_temperature_uncached(p: float) -> float:
    """Uncached bracketed-bisection inversion (equivalence-test surface).

    Domain-guarded (D3 ceremony, owner-ruled 2026-09-03; record
    ``docs/GT_PS2_D3_PROPS_FAIL_OPEN_CEREMONY_2026-09-03.md``).  The loop
    below is a pure bisection on ``[TT, TC - 0.1]`` with no bracket test, so
    before D3 it FAILED OPEN in three ways, all silently:

    * ``p`` below ``Psat(TT)`` — every comparison takes the ``hi = mid``
      branch and the loop returns the bracket floor.  Measured: ``p = 100 Pa``
      returned 273.16000000272 K, whose ``Psat`` is 611.65 Pa, a **relative
      residual of 5.12**; at ``p = 1e-3 Pa`` the residual is 6.1e5;
    * ``p`` above ``Psat(TC - 0.1)`` — mirror image, returns 646.99599999728 K
      for every request up to and beyond ``PC``;
    * ``p`` non-finite — ``NaN`` compares ``False`` against everything, so the
      loop again returns 273.16000000272 K.  A ``NaN`` in, a plausible
      *finite* temperature out, is the worst of the three.

    All three now raise ``ValueError``.  The honest domain is unchanged: the
    bisection arithmetic is untouched and every ``p`` inside
    ``[Psat(TT), Psat(TC - 0.1)]`` returns the identical double.  The two
    endpoint evaluations are served from ``saturation``'s exact-bit cache
    after the first call.
    """
    lo, hi = TT, TC - 0.1
    if not math.isfinite(p):
        raise ValueError(
            f"water saturation_temperature requires a finite pressure, "
            f"got p={p!r} Pa; the bisection cannot bracket a non-finite "
            f"target and would return the bracket floor {lo} K"
        )
    p_lo, p_hi = saturation_pressure(lo), saturation_pressure(hi)
    if p < p_lo or p > p_hi:
        raise ValueError(
            f"no water saturation temperature at p={p} Pa: the inversion "
            f"bracket [{lo} K, {hi} K] carries Psat in [{p_lo} Pa, {p_hi} Pa] "
            f"and the request lies outside it; the bisection would return the "
            f"clamped bracket end, not a root"
        )
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if saturation_pressure(mid) < p:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-8:
            break
    return 0.5 * (lo + hi)


def dh_vap(T: float) -> float:
    """Latent heat of vaporisation at T, J/kg (datum-free)."""
    return saturation(T).dh_vap


# ---------------------------------------------------------------------------
# Full single-phase state
# ---------------------------------------------------------------------------
def state(T: float, rho: float) -> MassPhaseState:
    delta, tau = rho / RHOC, TC / T
    res = residual(delta, tau)
    Z = res.compressibility
    p = rho * R * T * Z
    h = R * T * (1.0 + tau * (res.phi0_t + res.phir_t) + delta * res.phir_d)
    u = R * T * (tau * (res.phi0_t + res.phir_t))
    s = R * (tau * (res.phi0_t + res.phir_t) - (res.phi0 + res.phir))
    cv = -R * tau**2 * (res.phi0_tt + res.phir_tt)
    num = (1.0 + delta * res.phir_d - delta * tau * res.phir_dt) ** 2
    den = 1.0 + 2.0 * delta * res.phir_d + delta**2 * res.phir_dd
    cp = cv + R * num / den
    lnf = math.log(rho * R * T) + res.phir + delta * res.phir_d
    return MassPhaseState(
        T=T,
        rho_mass=rho,
        p=p,
        Z=Z,
        h_mass=h,
        u_mass=u,
        s_mass=s,
        cv_mass=cv,
        cp_mass=cp,
        ln_fugacity=lnf,
        molar_mass=M,
    )


def speed_of_sound(T: float, rho: float) -> float:
    """m/s — used by the IAPWS verification oracle."""
    delta, tau = rho / RHOC, TC / T
    res = residual(delta, tau)
    num = (1.0 + delta * res.phir_d - delta * tau * res.phir_dt) ** 2
    den = tau**2 * (res.phi0_tt + res.phir_tt)
    w2 = R * T * (1.0 + 2.0 * delta * res.phir_d + delta**2 * res.phir_dd - num / den)
    return math.sqrt(w2)


_STATE_TP_CACHE_MAXSIZE = 4096


def caloric_datum_signature() -> tuple[bytes, ...]:
    """Exact signature of the IAPWS ideal terms carrying its caloric datum."""

    return float_bits_key(*_N0)


def state_Tp(T: float, p: float, phase: str) -> MassPhaseState:
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
) -> MassPhaseState:
    # The original values, not decoded or canonicalized key values, enter EOS.
    return state(T, density(T, p, phase))


def clear_state_tp_cache() -> None:
    """Drop every cached pure-water ``(T,p,phase,datum)`` state."""

    _state_Tp_cached.cache_clear()


def state_tp_cache_info():
    """Return bounded-cache hit/miss/size diagnostics."""

    return _state_Tp_cached.cache_info()
