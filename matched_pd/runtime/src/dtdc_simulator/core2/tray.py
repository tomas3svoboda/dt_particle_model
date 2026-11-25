r"""Single quasi-steady equilibrium tray — the thin vertical slice.

This is the smallest instance of the whole-DTDC pattern the new core will use:
one conservative control volume with a **quasi-steady (algebraic) binary gas**
(PHY-004: no gas mass/energy storage — Coletto's steady bed-gas limit), coupled
to a well-mixed solid holdup, closed by local phase equilibrium (PHY-039 spirit)
and one first-law energy balance, and solved as a Newton residual system. It is
deliberately NOT the full model: the resolved particle/front (PHY-013/033/038),
the separate solid/gas energy nodes with interfacial transfer, the wall node
(PHY-008), and the nodal pressure network (PHY-017) all arrive in later phases.

What it demonstrates now, end to end, on the frozen property authorities:
  * both pure-fluid authorities in anger — Span-Wagner n-hexane (PHY-049) and
    IAPWS-95 water (PHY-030): saturation pressures and phase enthalpies;
  * the binary water/n-hexane DT gas with no noncondensable (PHY-011);
  * retained-water equilibrium via the modified-Luikov isotherm (PHY-031) and
    its PHY-021 caloric consequence (retained water is calorically liquid-like,
    zero net excess binding enthalpy — no separate sorption-heat source);
  * latent heat carried ONLY as the authority vapor-minus-liquid enthalpy
    difference (FLUX-ENERGY-INTERFACE: no duplicate latent source);
  * a conserved component + energy ledger that closes to ~1e-10.

Slice simplifications (documented, not frozen claims): one common tray
temperature for solid and exit gas (the full model resolves separate nodes with
an interfacial transfer law); hexane activity a_h=1 while mobile/attached hexane
is present (PHY-046 film-active high-loading regime); ideal partial-pressure
gas mixing (the frozen PHY-053 cross-virial correction is a later layer).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp
from dtdc_simulator.core2.props import water as wa

ATM_PA = 101325.0
M_HEX = hx.M     # kg/mol
M_WAT = wa.M     # kg/mol


# ---------------------------------------------------------------------------
# Introspectable parameter set (registry-friendly; the calibration/sensitivity
# tooling in scripts/ + calibration/ will read these, not hard-coded literals).
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TrayParams:
    """Physical/closure parameters for one equilibrium tray."""

    pressure_pa: float = ATM_PA           # PHY-017 local tray pressure
    cp_dry_meal: float = 2000.0           # J/(kg K), composite dry meal (PHY-030 uncertain const)
    T_ref_solid: float = 298.15           # K, dry-meal sensible-enthalpy datum
    luikov_A1: float = 0.880              # PHY-031 retained-water isotherm
    luikov_A2: float = 12.184
    w_water_cap: float = 0.2356709040     # PHY-031 evidence-cap capacity, kg/kg dry
    w_water_ref: float = 0.023            # PHY-031 qualified lower loading, kg/kg dry

    def as_dict(self) -> dict[str, float]:
        """Flat introspection view for the parameter registry."""
        return {
            "pressure_pa": self.pressure_pa,
            "cp_dry_meal": self.cp_dry_meal,
            "T_ref_solid": self.T_ref_solid,
            "luikov_A1": self.luikov_A1,
            "luikov_A2": self.luikov_A2,
            "w_water_cap": self.w_water_cap,
            "w_water_ref": self.w_water_ref,
        }


@dataclass(frozen=True)
class SolidStream:
    """Descending solid stream (per unit time)."""

    F_dm: float       # kg/s dry composite meal
    X_w: float        # kg retained water / kg dry meal
    X_h: float        # kg n-hexane / kg dry meal (retained + mobile/attached, liquid)
    T: float          # K


@dataclass(frozen=True)
class GasStream:
    """Ascending binary vapor stream (per unit time)."""

    G_w: float        # kg/s water vapor
    G_h: float        # kg/s n-hexane vapor
    T: float          # K


@dataclass(frozen=True)
class TraySolution:
    T: float                # K, common tray temperature
    solid_out: SolidStream
    gas_out: GasStream
    a_water: float          # retained-water activity at exit
    y_hexane: float         # exit gas hexane mole fraction
    converged: bool
    iterations: int
    residual_norm: float
    ledger: "TrayLedger" = field(default=None)  # type: ignore[assignment]

    @property
    def film_regime_valid(self) -> bool:
        """The a_h=1 film closure (PHY-046) only holds while mobile/attached
        n-hexane remains. If the solve drives X_h to zero the particle has dried
        out and the dry-shell diffusion closure (PHY-019, Phase 1) is required;
        the tray fails closed rather than returning an off-regime state."""
        return self.solid_out.X_h > 1.0e-9

    @property
    def physically_qualifying(self) -> bool:
        """This deliberately simplified Phase-0 slice is never a release model."""
        return False


# ---------------------------------------------------------------------------
# Constitutive closures
# ---------------------------------------------------------------------------
def luikov_activity(W_w: float, p: TrayParams) -> float:
    """Retained-water activity a_w(W_w) — inverse modified-Luikov (PHY-031).

    ln a_w = (1 - A1/W_w)/A2.  Both evidence bounds fail closed: this
    verification slice has neither the retained-cap active set nor an approved
    low-moisture continuation, so silently clipping either boundary would make
    the property oracle look more complete than it is.
    """
    return sp.water_activity(
        W_w,
        sp.LuikovParams(
            A1=p.luikov_A1,
            A2=p.luikov_A2,
            W_cap=p.w_water_cap,
            W_ref=p.w_water_ref,
        ),
    )


def _mole_fractions(G_w: float, G_h: float) -> tuple[float, float]:
    n_w, n_h = G_w / M_WAT, G_h / M_HEX
    tot = n_w + n_h
    if tot <= 0.0:
        return 0.5, 0.5
    return n_w / tot, n_h / tot


# --- stream enthalpies (J/kg of the carrier), on the frozen authorities ------
def _h_dry_meal(T: float, p: TrayParams) -> float:
    return p.cp_dry_meal * (T - p.T_ref_solid)


def _h_water_liquid(T: float, p: TrayParams) -> float:
    return wa.state_Tp(T, p.pressure_pa, "liquid").h_mass


def _h_hexane_liquid(T: float, p: TrayParams) -> float:
    return hx.state_Tp(T, p.pressure_pa, "liquid").h_mass


def _h_water_vapor(T: float, p_partial: float) -> float:
    # subsaturated water vapor at its partial pressure (near-ideal at DT range)
    return wa.state_Tp(T, max(p_partial, 1.0), "vapor").h_mass


def _h_hexane_vapor(T: float) -> float:
    # a_h = 1 -> hexane vapor is the saturated vapor at T (PHY-046 film regime)
    sat = hx.saturation(T)
    return hx.state(T, sat.rho_vapor).h_mass


# ---------------------------------------------------------------------------
# Energy tally (independent of the solver — used by the ledger gate)
# ---------------------------------------------------------------------------
def _stream_enthalpy_flow(
    solid: SolidStream, gas: GasStream, p: TrayParams, *, gas_partial_water: float
) -> float:
    """Total enthalpy flow (W) of a (solid + gas) side on the common datum."""
    h_solid = (
        solid.F_dm * _h_dry_meal(solid.T, p)
        + solid.F_dm * solid.X_w * _h_water_liquid(solid.T, p)
        + solid.F_dm * solid.X_h * _h_hexane_liquid(solid.T, p)
    )
    h_gas = (
        gas.G_w * _h_water_vapor(gas.T, gas_partial_water)
        + gas.G_h * _h_hexane_vapor(gas.T)
    )
    return h_solid + h_gas


@dataclass(frozen=True)
class TrayLedger:
    """Extensive conservation residuals, recomputed independently of the solve."""

    dry_meal_kg_s: float
    water_kg_s: float
    hexane_kg_s: float
    energy_w: float

    @property
    def max_mass_residual(self) -> float:
        return max(abs(self.dry_meal_kg_s), abs(self.water_kg_s), abs(self.hexane_kg_s))


# ---------------------------------------------------------------------------
# Residual system and Newton solve
# ---------------------------------------------------------------------------
def _residuals(
    u: list[float],
    solid_in: SolidStream,
    gas_in: GasStream,
    Q_w: float,
    p: TrayParams,
) -> list[float]:
    """5 residuals in the 5 unknowns u = [T, X_h_out, W_w_out, G_h_out, G_w_out]."""
    T, X_h_out, W_w_out, G_h_out, G_w_out = u
    F = solid_in.F_dm

    solid_out = SolidStream(F_dm=F, X_w=W_w_out, X_h=X_h_out, T=T)
    gas_out = GasStream(G_w=G_w_out, G_h=G_h_out, T=T)

    # component mass balances (kg/s): in = out
    r_hex = (F * solid_in.X_h + gas_in.G_h) - (F * X_h_out + G_h_out)
    r_wat = (F * solid_in.X_w + gas_in.G_w) - (F * W_w_out + G_w_out)

    # phase equilibrium at exit T and tray pressure (binary partial pressures)
    y_w, y_h = _mole_fractions(G_w_out, G_h_out)
    p_h = y_h * p.pressure_pa
    p_w = y_w * p.pressure_pa
    a_w = luikov_activity(W_w_out, p)
    r_vle_h = p_h - hx.saturation_pressure(T)           # a_h = 1 (film regime)
    r_vle_w = p_w - a_w * wa.saturation_pressure(T)

    # first-law energy balance (W): H_in + Q = H_out
    h_in = _stream_enthalpy_flow(
        solid_in, gas_in, p, gas_partial_water=_mole_fractions(gas_in.G_w, gas_in.G_h)[0] * p.pressure_pa
    )
    h_out = _stream_enthalpy_flow(solid_out, gas_out, p, gas_partial_water=p_w)
    r_energy = (h_in + Q_w) - h_out

    return [r_hex, r_wat, r_vle_h, r_vle_w, r_energy]


def _numeric_jacobian(f, u: list[float], f0: list[float]) -> list[list[float]]:
    """Central-difference Jacobian. Central differencing plus a step large
    enough to clear the saturation-solver's residual floor keeps the VLE rows'
    temperature column clean (see module tests)."""
    n = len(u)
    # per-variable step scales: T ~ 1e-3 K, loadings ~ 1e-6, flows ~ 1e-6 relative
    scales = [1.0e-3, 1.0e-8, 1.0e-8, 1.0e-9, 1.0e-9]
    J = [[0.0] * n for _ in range(n)]
    for j in range(n):
        du = max(scales[j], 1.0e-8 * abs(u[j]))
        up, um = list(u), list(u)
        up[j] += du
        um[j] -= du
        fp, fm = f(up), f(um)
        for i in range(n):
            J[i][j] = (fp[i] - fm[i]) / (2.0 * du)
    return J


def _solve_linear(J: list[list[float]], b: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting (n=5)."""
    n = len(b)
    A = [row[:] + [b[i]] for i, row in enumerate(J)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(A[r][col]))
        A[col], A[piv] = A[piv], A[col]
        pivot = A[col][col]
        if pivot == 0.0:
            raise ZeroDivisionError("singular Jacobian in tray Newton solve")
        for r in range(n):
            if r == col:
                continue
            factor = A[r][col] / pivot
            for c in range(col, n + 1):
                A[r][c] -= factor * A[col][c]
    return [A[i][n] / A[i][i] for i in range(n)]


def solve_tray(
    solid_in: SolidStream,
    gas_in: GasStream,
    Q_w: float,
    params: TrayParams | None = None,
    *,
    tol: float = 1.0e-11,
    max_iter: int = 100,
    seed: list[float] | None = None,
) -> TraySolution:
    """Solve one quasi-steady equilibrium tray and return the closed solution.

    Newton on the 5-residual system with a numerically differentiated Jacobian
    (the analytic sparse Jacobian is the Phase-4 speed work; here correctness and
    the closed ledger are what matter).
    """
    p = params or TrayParams()
    F = solid_in.F_dm
    if not p.w_water_ref <= solid_in.X_w <= p.w_water_cap:
        raise ValueError(
            "tray inlet retained water must lie on the qualified Luikov branch; "
            "the Phase-0 slice has no low-moisture continuation or cap active set"
        )

    # seed: gas unchanged, solid keeps its load, T at the inlet solid
    # temperature — or a caller-supplied warm start (owner freeze-opening
    # 2026-08-06: optional, backward-compatible; seed=None is
    # bit-identical to the frozen behavior).
    if seed is None:
        u = [solid_in.T, solid_in.X_h, solid_in.X_w, gas_in.G_h, gas_in.G_w]
    else:
        if len(seed) != 5 or not all(math.isfinite(v) for v in seed):
            raise ValueError("warm-start seed must be five finite values")
        u = list(seed)

    def f(uu: list[float]) -> list[float]:
        return _residuals(uu, solid_in, gas_in, Q_w, p)

    converged = False
    it = 0
    rnorm = math.inf
    # T, X_h and gas flows stay strictly positive.  W_w additionally remains
    # inside the evidence-bounded Luikov branch; this slice fails closed rather
    # than pretending to implement the later retained-cap active set.
    positive = (0, 1, 3, 4)
    for it in range(1, max_iter + 1):
        f0 = f(u)
        rnorm = _weighted_norm(f0, F, p)
        if rnorm < tol:
            converged = True
            break
        J = _numeric_jacobian(f, u, f0)
        # Newton step: J s = f0, update u <- u - s
        s = _solve_linear(J, f0)
        # cap the temperature move to keep the stiff VLE linearisation in range
        if abs(s[0]) > 25.0:
            scale = 25.0 / abs(s[0])
            s = [si * scale for si in s]
        # fraction-to-boundary: never step a positive unknown to <= 0
        alpha = 1.0
        for idx in positive:
            if s[idx] > 0.0 and u[idx] - s[idx] <= 0.0:
                alpha = min(alpha, 0.9 * u[idx] / s[idx])
        if s[2] > 0.0 and u[2] - s[2] < p.w_water_ref:
            alpha = min(alpha, 0.9 * (u[2] - p.w_water_ref) / s[2])
        if s[2] < 0.0 and u[2] - s[2] > p.w_water_cap:
            alpha = min(alpha, 0.9 * (p.w_water_cap - u[2]) / (-s[2]))
        u = [u[i] - alpha * s[i] for i in range(5)]

    T, X_h_out, W_w_out, G_h_out, G_w_out = u
    solid_out = SolidStream(F_dm=F, X_w=W_w_out, X_h=X_h_out, T=T)
    gas_out = GasStream(G_w=G_w_out, G_h=G_h_out, T=T)
    y_w, y_h = _mole_fractions(G_w_out, G_h_out)
    ledger = tray_ledger(solid_in, gas_in, solid_out, gas_out, Q_w, p)

    return TraySolution(
        T=T, solid_out=solid_out, gas_out=gas_out,
        a_water=luikov_activity(W_w_out, p), y_hexane=y_h,
        converged=converged, iterations=it, residual_norm=rnorm, ledger=ledger,
    )


def _weighted_norm(r: list[float], F: float, p: TrayParams) -> float:
    """Scale residuals to comparable magnitude for the convergence test."""
    mass_scale = max(F, 1.0e-6)
    press_scale = p.pressure_pa
    energy_scale = max(F * 1.0e5, 1.0)          # ~ F * cp * dT scale
    scaled = [
        r[0] / mass_scale, r[1] / mass_scale,
        r[2] / press_scale, r[3] / press_scale,
        r[4] / energy_scale,
    ]
    return math.sqrt(sum(x * x for x in scaled))


def tray_ledger(
    solid_in: SolidStream,
    gas_in: GasStream,
    solid_out: SolidStream,
    gas_out: GasStream,
    Q_w: float,
    p: TrayParams,
) -> TrayLedger:
    """Recompute extensive in-minus-out for every conserved quantity.

    Independent of the Newton residuals: it re-tallies each stream from the
    property authorities, so a nonzero value flags an inconsistency between the
    closure and the authorities, not merely an unconverged step.
    """
    dry = solid_in.F_dm - solid_out.F_dm
    water = (
        (solid_in.F_dm * solid_in.X_w + gas_in.G_w)
        - (solid_out.F_dm * solid_out.X_w + gas_out.G_w)
    )
    hexane = (
        (solid_in.F_dm * solid_in.X_h + gas_in.G_h)
        - (solid_out.F_dm * solid_out.X_h + gas_out.G_h)
    )
    p_w_in = _mole_fractions(gas_in.G_w, gas_in.G_h)[0] * p.pressure_pa
    p_w_out = _mole_fractions(gas_out.G_w, gas_out.G_h)[0] * p.pressure_pa
    h_in = _stream_enthalpy_flow(solid_in, gas_in, p, gas_partial_water=p_w_in)
    h_out = _stream_enthalpy_flow(solid_out, gas_out, p, gas_partial_water=p_w_out)
    energy = (h_in + Q_w) - h_out
    return TrayLedger(
        dry_meal_kg_s=dry, water_kg_s=water, hexane_kg_s=hexane, energy_w=energy
    )
