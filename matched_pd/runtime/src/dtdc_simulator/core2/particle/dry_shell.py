r"""PHY-019 dry-shell n-hexane transport + PHY-020 local-equilibrium storage.

Cardarelli (2002) isothermal limit, in conservation form (registry
`PB-SOLVENT-DRY-STORAGE` / `PB-SOLVENT-DRY-FLUX`):

    d/dt[ eps_g*c_g + rho_dm,p*W_h(a_h,T,w_o) ]  =  (1/r^2) d/dr[ r^2 * D_eff * dc_g/dr ]

- pore-gas hexane mass density ``c_g`` (kg/m3) is the primary variable;
- retained n-hexane is local-equilibrium storage W_h (PHY-020), advanced through
  the conservation law with activity recovered by monotone inversion — the
  variable equilibrium capacity stays in STORAGE, never inside the divergence as
  a face mobility (PHY-020 rule);
- flux j = -D_eff * dc_g/dr with the effective (pore) diffusivity D_eff; the gas
  holdup eps_g appears only in storage;
- activity a_h = c_g / c_g^sat(T), c_g^sat = Psat_h(T)*M_h/(R*T) (Span-Wagner Psat).

Scope: isothermal, pure-n-hexane dry shell — exactly the qualified oracle case
PHY-019 freezes. Non-isothermal coupling and the water/binary-Maxwell-Stefan
bracket are later work; nothing here reaches into the wet core (no diffusion
operator crosses the front).

Numerics: conservative cell-centered finite volume (flux-form => exact discrete
conservation by telescoping), implicit backward-Euler, Newton with the analytic
tridiagonal Jacobian (nonlinear only through the storage capacity dC/dc_g).
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass
from functools import lru_cache
from typing import Sequence

from dtdc_simulator.core2.particle.grid import SphericalGrid
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp

R_GAS = 8.31451     # J/(mol K)


@dataclass(frozen=True)
class DryShellParams:
    eps_g: float = 0.141              # gas holdup = particle porosity
    rho_dm_p: float = 1159.65         # kg dry meal / m3 particle
    D_eff: float = 4.0e-10            # m2/s, Cardarelli effective pore diffusivity
    w_o: float = 0.0195               # residual-oil mass fraction
    gab: sp.GabParams = sp.GabParams()
    oil: sp.OilIsotherm = sp.OilIsotherm()


@lru_cache(maxsize=4096)
def cg_saturation(T: float) -> float:
    """Saturated pore-gas hexane mass density c_g^sat(T), kg/m3 (a_h=1).

    Cached: the saturation Newton is pure in T. Isothermal solves hit one entry;
    the non-isothermal model will replace this with the coupled front closure.
    """
    return hx.saturation_pressure(T) * hx.M / (R_GAS * T)


def activity(c_g: float, T: float) -> float:
    """a_h = c_g / c_g^sat(T)."""
    return c_g / cg_saturation(T)


def storage(c_g: float, T: float, p: DryShellParams) -> float:
    """C_h,d = eps_g*c_g + rho_dm,p*W_h(a_h,T,w_o), kg hexane / m3 particle."""
    a = activity(c_g, T)
    W = sp.retained_hexane(a, T, p.w_o, p.gab, p.oil)
    return p.eps_g * c_g + p.rho_dm_p * W


def storage_capacity(c_g: float, T: float, p: DryShellParams) -> float:
    """dC_h,d/dc_g (analytic) = eps_g + rho_dm,p * dW/da * da/dc_g."""
    cg_sat = cg_saturation(T)
    a = c_g / cg_sat
    dW_da = sp.retained_hexane_da(a, T, p.w_o, p.gab, p.oil)
    return p.eps_g + p.rho_dm_p * dW_da / cg_sat


def cg_from_storage(C: float, T: float, p: DryShellParams) -> float:
    """Invert C_h,d -> c_g (monotone bisection)."""
    lo, hi = 0.0, cg_saturation(T) * (1.0 / p.gab.K(T)) / activity(cg_saturation(T), T)
    # cheap safe upper bound: the activity ceiling K(T)*a<1 caps c_g
    hi = cg_saturation(T) / p.gab.K(T) * 0.999999
    if storage(hi, T, p) < C:
        raise ValueError(f"C={C} exceeds admissible dry-shell storage; reject, do not clamp")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if storage(mid, T, p) < C:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-16 * max(hi, 1e-30):
            break
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Diffusion operator (flux form) and backward-Euler step
# ---------------------------------------------------------------------------
def _diffusion_flux_divergence(
    grid: SphericalGrid, c: list[float], D: float, surface,
) -> tuple[list[float], list[float], list[float], list[float]]:
    """Return per-cell (net_flux, sub, diag, sup) for the diffusion operator.

    net_flux[i] = A_R*D*(c[i+1]-c[i])/dr_R - A_L*D*(c[i]-c[i-1])/dr_L
    plus the tridiagonal Jacobian of net_flux wrt c (sub/diag/sup), so the caller
    can assemble the implicit residual/Jacobian. ``surface`` is the r=R boundary:
      ("closed",)                 zero flux
      ("dirichlet", c_surface)    prescribed pore-gas density at r=R
      ("flux", F_surface)         prescribed outward flux at r=R
    """
    n = grid.n
    net = [0.0] * n
    sub = [0.0] * n   # d net[i]/d c[i-1]
    diag = [0.0] * n  # d net[i]/d c[i]
    sup = [0.0] * n   # d net[i]/d c[i+1]
    for i in range(n):
        # left face (i-1/2)
        if i > 0:
            dr_L = grid.centers[i] - grid.centers[i - 1]
            gL = grid.areas[i] * D / dr_L
            net[i] += -gL * (c[i] - c[i - 1])
            diag[i] += -gL
            sub[i] += gL
        # inner face at r=0 has area 0 -> no term (symmetry)
        # right face (i+1/2)
        if i < n - 1:
            dr_R = grid.centers[i + 1] - grid.centers[i]
            gR = grid.areas[i + 1] * D / dr_R
            net[i] += gR * (c[i + 1] - c[i])
            diag[i] += -gR
            sup[i] += gR
        else:
            # outermost right face at r=R
            if surface[0] == "closed":
                pass
            elif surface[0] == "dirichlet":
                c_s = surface[1]
                dr_R = grid.R - grid.centers[i]
                gR = grid.areas[-1] * D / dr_R
                net[i] += gR * (c_s - c[i])
                diag[i] += -gR
            elif surface[0] == "flux":
                # F_surface is outward flux (kg/m2/s); accumulation loses A_R*F
                net[i] += -grid.areas[-1] * surface[1]
            else:
                raise ValueError(f"unknown surface BC {surface[0]!r}")
    return net, sub, diag, sup


def _thomas(sub, diag, sup, rhs):
    """Solve a tridiagonal system in place (Thomas algorithm)."""
    n = len(rhs)
    cp = [0.0] * n
    dp = [0.0] * n
    cp[0] = sup[0] / diag[0]
    dp[0] = rhs[0] / diag[0]
    for i in range(1, n):
        m = diag[i] - sub[i] * cp[i - 1]
        cp[i] = sup[i] / m if i < n - 1 else 0.0
        dp[i] = (rhs[i] - sub[i] * dp[i - 1]) / m
    x = [0.0] * n
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def _validate_step_inputs(
    grid: SphericalGrid,
    c_prev: Sequence[float],
    dt: float,
    T: float,
    p: DryShellParams,
    surface,
    source: Sequence[float] | None,
    tol: float,
    max_iter: int,
) -> tuple[list[float], list[float], tuple]:
    """Validate one dry-shell step without altering any physical state.

    The GAB closure has a genuine singular boundary at ``K(T)*a_h = 1``.
    Newton iterates therefore have to remain in the open interval below that
    boundary; accepting or clipping an invalid iterate would change PHY-020.
    """
    if grid.n < 1:
        raise ValueError("dry-shell step needs at least one grid cell")
    if len(c_prev) != grid.n:
        raise ValueError("c_prev must contain one value per grid cell")
    if source is not None and len(source) != grid.n:
        raise ValueError("source must contain one value per grid cell")
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be finite and positive")
    if not math.isfinite(T) or not (0.0 < T < hx.TC):
        raise ValueError(f"T must be finite and below the n-hexane critical temperature {hx.TC} K")
    if not math.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive")
    try:
        iterations = operator.index(max_iter)
    except TypeError as exc:
        raise ValueError("max_iter must be a positive integer") from exc
    if isinstance(max_iter, bool) or iterations < 1:
        raise ValueError("max_iter must be a positive integer")

    parameter_values = {
        "eps_g": p.eps_g,
        "rho_dm_p": p.rho_dm_p,
        "D_eff": p.D_eff,
        "w_o": p.w_o,
    }
    if not all(math.isfinite(value) for value in parameter_values.values()):
        raise ValueError("dry-shell parameters must be finite")
    if p.eps_g < 0.0 or p.rho_dm_p <= 0.0 or p.D_eff < 0.0:
        raise ValueError("need eps_g >= 0, rho_dm_p > 0, and D_eff >= 0")

    csat = cg_saturation(T)
    K = p.gab.K(T)
    if not math.isfinite(csat) or csat <= 0.0 or not math.isfinite(K) or K <= 0.0:
        raise ValueError("invalid saturation/GAB state for dry-shell step")
    c_ceiling = csat / K

    state = list(c_prev)
    if not all(math.isfinite(value) for value in state):
        raise ValueError("c_prev values must be finite")
    if any(value < 0.0 or value >= c_ceiling for value in state):
        raise ValueError(
            "c_prev must be non-negative and below the GAB singular density; reject, do not clamp"
        )

    sources = [0.0] * grid.n if source is None else list(source)
    if not all(math.isfinite(value) for value in sources):
        raise ValueError("source values must be finite")

    if not isinstance(surface, (tuple, list)) or not surface:
        raise ValueError("surface BC must be a non-empty tuple or list")
    kind = surface[0]
    if kind == "closed":
        if len(surface) != 1:
            raise ValueError("closed surface BC takes no value")
        boundary = ("closed",)
    elif kind in {"dirichlet", "flux"}:
        if len(surface) != 2 or not isinstance(surface[1], (int, float)):
            raise ValueError(f"{kind} surface BC needs one finite numeric value")
        value = float(surface[1])
        if not math.isfinite(value):
            raise ValueError(f"{kind} surface BC value must be finite")
        if kind == "dirichlet" and not (0.0 <= value < c_ceiling):
            raise ValueError(
                "Dirichlet pore-gas density must be non-negative and below the GAB singular density"
            )
        boundary = (kind, value)
    else:
        raise ValueError(f"unknown surface BC {kind!r}")

    # Exercise the constitutive domain once up front (including the residual-oil
    # applicability guard) so bad inputs are rejected before Newton starts.
    for value in state:
        capacity = storage_capacity(value, T, p)
        if not math.isfinite(capacity) or capacity <= 0.0:
            raise ValueError("dry-shell storage capacity must be finite and positive")
    return state, sources, boundary


def _step_residual(
    grid: SphericalGrid,
    c: Sequence[float],
    C_prev: Sequence[float],
    dt: float,
    T: float,
    p: DryShellParams,
    surface,
    source: Sequence[float],
) -> tuple[list[float], list[float], list[float], list[float], float]:
    """Return the backward-Euler residual, flux Jacobian bands, and norm."""
    net, sub, diag, sup = _diffusion_flux_divergence(grid, list(c), p.D_eff, surface)
    residual = [
        grid.volumes[i] * (storage(c[i], T, p) - C_prev[i]) / dt
        - net[i]
        - grid.volumes[i] * source[i]
        for i in range(grid.n)
    ]
    norm = math.sqrt(
        sum(
            (residual[i] / max(grid.volumes[i], 1e-300)) ** 2
            for i in range(grid.n)
        )
        / grid.n
    )
    return residual, sub, diag, sup, norm


def step(
    grid: SphericalGrid,
    c_prev: Sequence[float],
    dt: float,
    T: float,
    p: DryShellParams,
    surface=("closed",),
    source: Sequence[float] | None = None,
    *,
    tol: float = 1e-12,
    max_iter: int = 50,
) -> list[float]:
    """One implicit backward-Euler step; returns the new pore-gas density field.

    Residual per cell:
      G_i = V_i*(storage(c_i) - C_prev_i)/dt - net_flux_i - V_i*S_i = 0.
    """
    c, S, boundary = _validate_step_inputs(
        grid, c_prev, dt, T, p, surface, source, tol, max_iter
    )
    n = grid.n
    C_prev = [storage(value, T, p) for value in c]
    c_ceiling = cg_saturation(T) / p.gab.K(T)
    rnorm = math.inf

    # ``max_iter`` counts Newton updates.  The extra residual evaluation after
    # the last update prevents a converged final iterate from being rejected,
    # while still guaranteeing that an unconverged state is never returned.
    for iteration in range(max_iter + 1):
        G, sub, diag, sup, rnorm = _step_residual(
            grid, c, C_prev, dt, T, p, boundary, S
        )
        if not math.isfinite(rnorm):
            raise RuntimeError("dry-shell Newton produced a non-finite residual")
        if rnorm < tol:
            return c
        if iteration == max_iter:
            break

        Jd = [
            grid.volumes[i] * storage_capacity(c[i], T, p) / dt - diag[i]
            for i in range(n)
        ]
        Jsub = [-sub[i] for i in range(n)]
        Jsup = [-sup[i] for i in range(n)]
        try:
            delta = _thomas(Jsub, Jd, Jsup, G)
        except ZeroDivisionError as exc:
            raise RuntimeError("dry-shell Newton Jacobian is singular") from exc
        if not all(math.isfinite(value) for value in delta):
            raise RuntimeError("dry-shell Newton produced a non-finite update")

        # Fraction-to-boundary damping keeps every iterate inside the exact
        # constitutive domain.  No negative density or GAB singularity is ever
        # clipped into a seemingly valid state.
        alpha = 1.0
        for value, correction in zip(c, delta):
            if correction > 0.0:
                alpha = min(alpha, 0.99 * value / correction)
            elif correction < 0.0:
                alpha = min(alpha, 0.99 * (c_ceiling - value) / -correction)
        if not math.isfinite(alpha) or alpha <= 0.0:
            raise RuntimeError(
                "dry-shell Newton has no admissible positive-density update; reject, do not clamp"
            )

        accepted = False
        for _ in range(60):
            trial = [value - alpha * correction for value, correction in zip(c, delta)]
            if all(0.0 <= value < c_ceiling and math.isfinite(value) for value in trial):
                *_, trial_norm = _step_residual(
                    grid, trial, C_prev, dt, T, p, boundary, S
                )
                if math.isfinite(trial_norm) and (trial_norm < rnorm or trial_norm < tol):
                    c = trial
                    accepted = True
                    break
            alpha *= 0.5
        if not accepted:
            raise RuntimeError(
                "dry-shell Newton line search failed inside the admissible storage domain"
            )

    raise RuntimeError(
        f"dry-shell Newton did not converge after {max_iter} updates "
        f"(residual density norm={rnorm:.6e}, tolerance={tol:.6e})"
    )


def total_hexane(grid: SphericalGrid, c: list[float], T: float, p: DryShellParams) -> float:
    """Total particle n-hexane inventory, kg (integral of storage over volume)."""
    return sum(V * storage(ci, T, p) for V, ci in zip(grid.volumes, c))


def steady_diffusion_solve(
    grid: SphericalGrid, D: float, source: list[float], surface
) -> list[float]:
    """Solve the linear steady diffusion `net_flux_i + V_i*S_i = 0` directly.

    Storage-free: isolates the spatial FV operator for manufactured-solution
    order verification (the diffusion operator is linear in c).
    """
    n = grid.n
    net0, sub, diag, sup = _diffusion_flux_divergence(grid, [0.0] * n, D, surface)
    rhs = [-(net0[i] + grid.volumes[i] * source[i]) for i in range(n)]
    return _thomas(sub, diag, sup, rhs)
