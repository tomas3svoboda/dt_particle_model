"""Two-liquid water/n-hexane heteroazeotrope solver (T2b, PHY-007/PHY-053),
generalized to the ACTIVITY-SCALED interface of both species (T3b).

Solves the frozen two-liquid condition at a given pressure:

    y_h * phi_h(T, P, y) * P = a_h * f_h^{L,pure}(T, P)
    y_w * phi_w(T, P, y) * P = a_w * f_w^{L,pure}(T, P),   y_w = 1 - y_h

with pure compressed-liquid fugacities from the frozen Helmholtz
authorities (Poynting implicit in the EOS at (T, P) liquid) and vapor
fugacity coefficients from the frozen PHY-053 cross-virial closure
(`props/binary_gas.py`).  This is the certified production form of the
2026-08-06 diagnostic probe.

`a_h = a_w = 1` is the certified two-liquid heteroazeotrope, bit-exactly.
Lowering either activity is the sorption-suppressed interface of a solid
whose corresponding free liquid phase is exhausted:

- `a_w < 1` (T2b): free water gone, retained water on the PHY-031
  modified-Luikov branch.  Raises T_I above the azeotrope.
- `a_h < 1` (T3b, NEW): free/mobile n-hexane gone — X_h below the
  PHY-046 critical loading X_c(T) — so the surface hexane is the
  PHY-043/048 GAB sorbate and its escaping tendency is a_h(W_h, T, w_o)
  from `props/sorption.activity_from_retained`.  In the limit a_h -> 0
  the hexane equality degenerates (y_h -> 0) and the solve returns the
  pure SORBED-WATER boiling point a_w * f_w^L(T) = P.  That limit is the
  physical content of the DT meal outlet: at the PHY-031 evidence cap
  (a_w = 0.799) it is ~379.5 K at 1 atm, against a measured 378-379 K.

WHY THE BRACKET CEILING IS NOW COMPUTED, NOT DECLARED.  The old constant
`_BRACKET_HIGH_K = 340.5` encoded "strictly below the pure-hexane
boiling point at 1 atm" — correct ONLY at a_h = 1, where the hexane
equality loses its y_h < 1 solution at 343.46 K.  The true ceiling is
the temperature at which a_h * f_h^L(T) = P, which MOVES with a_h
(a_h = 0.2 -> 412 K).  Pinning it at 340.5 therefore capped every
interface temperature 38 K below the benchmark meal outlet as a side
effect of the a_h = 1 assumption.  The ceiling is now solved from the
authorities themselves; at a_h = 1.0 the frozen constants are used
unchanged so the certified pin is bit-identical.

DOMAIN HAZARD, MEASURED AND GUARDED (do not remove).  The frozen hexane
authority's `density(T, P, "liquid")` SILENTLY returns the VAPOUR root
above ~471 K at 1 atm (rho 4678 -> 26.2 kg/m3 between 470 and 472 K)
with no exception, so `_pure_liquid_fugacity` collapses from 1.18 MPa to
99.8 kPa and every downstream number becomes quietly wrong.  Water does
the same above ~590 K.  `_pure_liquid_fugacity` therefore verifies the
root really is the dense one before returning.  This guard cannot fire
anywhere on the frozen a_h = 1 path (ceiling 340.5 K).

Gate context (Packet T2 / pre-flight ask A, owner-blessed 2026-08-06):
the Kemper atmospheric oracle requires |T_az - 335.15 K| <= 2 K and
hexane mass fraction 0.94 +/- 0.02, swept over the frozen k_wh interval,
with the F3 engineering-status disclosure.  `ideal_vapor=True` gives the
mandatory ablation oracle (phi_i = 1).

DISCLOSED CONVENTION BIAS.  The GAB (PHY-043/048) and modified-Luikov
(PHY-031) isotherms are sorption fits on the PRESSURE-ratio convention
a = p/Psat, and are used here on the FUGACITY-ratio convention
a = f/f^{L,pure}.  The two differ by f^L/Psat, which this module exposes
as `activity_convention_ratio(T, P)`: 0.9961 for water and 0.9547 for
n-hexane at 328 K / 1 atm.  The water side has carried this since T2b;
the hexane side inherits it.  Named, not hidden; ruling deferred.

physically_qualifying is False: engineering closure evidence, not plant
validation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from .props import binary_gas
from .props import hexane as _hexane
from .props import water as _water

# Memo size for the pure-liquid fugacity (see `_pure_liquid_fugacity`).
# Bounded so a long march cannot grow it without limit; the same pattern
# `particle/dry_shell.cg_saturation` already uses.
_FUGACITY_CACHE_SIZE = 200_000

# Bracket strictly below the pure-hexane boiling point at 1 atm
# (341.9 K): above it the hexane equality has no y_h < 1 solution.
# FROZEN a_h = 1 PATH ONLY — see `_bracket_for` for a_h < 1.
_BRACKET_LOW_K = 328.0
_BRACKET_HIGH_K = 340.5
_TEMPERATURE_TOLERANCE_K = 1.0e-10
# Bracket ENDPOINT search only — never the converged root (see `_bracket_for`).
_BRACKET_SEARCH_TOLERANCE_K = 1.0e-4
_COMPOSITION_TOLERANCE = 1.0e-14
_MAX_FIXED_POINT_ITERATIONS = 200
_MAX_BISECTIONS = 200

# Hard ceiling for the generalized bracket.  Set by the MEASURED collapse
# of the hexane compressed-liquid root at ~471 K (module docstring), with
# margin.  The `_pure_liquid_fugacity` root guard is the real defence;
# this keeps the bisection from ever probing there.
_LIQUID_REFERENCE_MAX_K = 465.0
# Fractional standoff from the y_h -> 1 singularity of the hexane
# equality, so the fixed point stays inside (0, 1).
_CEILING_STANDOFF = 0.995


class HeteroazeotropeError(RuntimeError):
    """The two-liquid solve failed a contract; the result is unusable."""


@dataclass(frozen=True)
class HeteroazeotropeState:
    """Converged two-liquid co-boiling point."""

    temperature_k: float
    y_hexane: float
    hexane_mass_fraction: float
    pressure_pa: float
    k_wh: float
    ideal_vapor: bool
    water_residual_pa: float
    water_activity: float = 1.0
    hexane_activity: float = 1.0

    @property
    def physically_qualifying(self) -> bool:
        return False

    @property
    def sorption_suppressed_hexane(self) -> bool:
        """True when the hexane equality ran below its pure-liquid pin —
        i.e. the surface carries GAB sorbate, not mobile n-hexane."""
        return self.hexane_activity < 1.0


@lru_cache(maxsize=_FUGACITY_CACHE_SIZE)
def _hexane_liquid_fugacity(T: float, pressure_pa: float) -> float:
    return _liquid_fugacity(_hexane, T, pressure_pa)


@lru_cache(maxsize=_FUGACITY_CACHE_SIZE)
def _water_liquid_fugacity(T: float, pressure_pa: float) -> float:
    return _liquid_fugacity(_water, T, pressure_pa)


def _liquid_fugacity(module, T: float, pressure_pa: float) -> float:
    """Pure compressed-liquid fugacity, with the liquid root VERIFIED.

    The frozen authorities do not raise when the metastable liquid root
    disappears (hexane above ~471 K at 1 atm, water above ~590 K) — they
    quietly hand back the vapour root.  MEASURED: hexane's density at
    1 atm falls 4678 -> 26.2 mol/m3 between 470 and 472 K and its
    fugacity collapses 1.18 MPa -> 99.8 kPa, with no exception.

    The test is `rho > rho_c`: a liquid branch is denser than critical,
    a low-pressure vapour branch is far less dense.  That is a
    structural statement, not a tuned cutoff, and it is cheap — the
    earlier saturated-vapour-density form cost a full saturation Newton
    and was ~half the module's wall time.  Margin is large everywhere we
    solve (hexane 4678 vs rho_c 2706 mol/m3 at the 470 K edge; water
    668 vs 322 kg/m3 at 580 K), and it never fires below the frozen
    a_h = 1 bracket ceiling of 340.5 K.  Note each authority reports in
    its own convention — hexane molar, water mass — and `RHOC` follows
    the same convention, so the comparison is unit-consistent in both.
    """
    rho_liquid = module.density(T, pressure_pa, "liquid")
    if not rho_liquid > module.RHOC:
        raise HeteroazeotropeError(
            f"{module.__name__.rsplit('.', 1)[-1]} compressed-liquid root "
            f"collapsed to the vapour branch at T = {T} K, P = {pressure_pa} "
            f"Pa (rho = {rho_liquid} <= rho_c = {module.RHOC}); the "
            "pure-liquid reference state does not exist here — reject, do not "
            "extrapolate"
        )
    return math.exp(module.ln_fugacity(T, rho_liquid))


def _pure_liquid_fugacity(module, T: float, pressure_pa: float) -> float:
    """Memoized `_liquid_fugacity`.

    PURE in (T, P), so the cache returns the SAME float object and every
    result stays bit-identical — this is a speed change only.  It matters
    a great deal: profiling one solve showed 169 calls costing 98 % of
    the wall time (0.462 s of 0.470 s), roughly half of that spent in the
    saturation solve the liquid-root guard needs.  Marching a tray
    revisits temperatures constantly, so the hit rate is high.
    """
    if module is _hexane:
        return _hexane_liquid_fugacity(T, pressure_pa)
    if module is _water:
        return _water_liquid_fugacity(T, pressure_pa)
    return _liquid_fugacity(module, T, pressure_pa)


def activity_convention_ratio(
    T: float, pressure_pa: float = 101325.0
) -> tuple[float, float]:
    """(water, hexane) f^{L,pure}/Psat at (T, P) — the DISCLOSED bias
    between the sorption fits' pressure-ratio activity convention and
    this module's fugacity-ratio one (module docstring)."""
    return (
        _pure_liquid_fugacity(_water, T, pressure_pa)
        / _water.saturation_pressure(T),
        _pure_liquid_fugacity(_hexane, T, pressure_pa)
        / _hexane.saturation_pressure(T),
    )


def _phi(T, pressure_pa, y_hexane, k_wh, ideal_vapor):
    if ideal_vapor:
        return 1.0, 1.0
    state = binary_gas.state(T, pressure_pa, 1.0 - y_hexane, y_hexane, k_wh=k_wh)
    return state.phi_water, state.phi_hexane


def _y_hexane_from_hexane_equality(
    T, pressure_pa, k_wh, ideal_vapor, hexane_activity=1.0
):
    """Fixed point y_h = a_h*f_h^L / (phi_h(y) * P); None outside (0, 1).

    At `hexane_activity = 1.0` the multiplication is the IEEE identity
    (1.0 * x == x for every finite x), so the frozen path is bit-exact.
    """
    y_h = 0.8
    f_h = _pure_liquid_fugacity(_hexane, T, pressure_pa)
    for _ in range(_MAX_FIXED_POINT_ITERATIONS):
        _, phi_h = _phi(T, pressure_pa, y_h, k_wh, ideal_vapor)
        y_new = hexane_activity * f_h / (phi_h * pressure_pa)
        if not 0.0 < y_new < 1.0:
            return None
        if abs(y_new - y_h) < _COMPOSITION_TOLERANCE:
            return y_new
        y_h = y_new
    raise HeteroazeotropeError("hexane-equality fixed point did not converge")


def _water_residual(
    T, pressure_pa, k_wh, ideal_vapor, water_activity=1.0, hexane_activity=1.0
):
    y_h = _y_hexane_from_hexane_equality(
        T, pressure_pa, k_wh, ideal_vapor, hexane_activity
    )
    if y_h is None:
        return None, None
    f_w = _pure_liquid_fugacity(_water, T, pressure_pa)
    phi_w, _ = _phi(T, pressure_pa, y_h, k_wh, ideal_vapor)
    return (1.0 - y_h) * phi_w * pressure_pa - water_activity * f_w, y_h


def hexane_equality_ceiling_k(
    hexane_activity: float, pressure_pa: float
) -> float:
    """T at which a_h * f_h^L(T) = P — where the hexane equality loses its
    y_h < 1 solution.  This is what `_BRACKET_HIGH_K = 340.5` hard-coded
    for the a_h = 1 case (true value 343.46 K at 1 atm); it MOVES with
    a_h, which is why the constant had to stop being a constant.

    Capped at `_LIQUID_REFERENCE_MAX_K` because below some a_h the
    equality has no ceiling inside the authorities' liquid domain at all.
    """
    lo, hi = _BRACKET_LOW_K, _LIQUID_REFERENCE_MAX_K
    try:
        if hexane_activity * _pure_liquid_fugacity(_hexane, hi, pressure_pa) < pressure_pa:
            return hi
    except HeteroazeotropeError:
        pass
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        try:
            value = hexane_activity * _pure_liquid_fugacity(
                _hexane, mid, pressure_pa
            )
        except HeteroazeotropeError:
            hi = mid
            continue
        if value < pressure_pa:
            lo = mid
        else:
            hi = mid
        if hi - lo < _TEMPERATURE_TOLERANCE_K:
            break
    return 0.5 * (lo + hi)


def _hexane_equality_admissible(
    T, pressure_pa, k_wh, ideal_vapor, hexane_activity
) -> bool:
    """Does the hexane equality have a y_h in (0, 1) at this T?"""
    try:
        return (
            _y_hexane_from_hexane_equality(
                T, pressure_pa, k_wh, ideal_vapor, hexane_activity
            )
            is not None
        )
    except HeteroazeotropeError:
        return False


def _bracket_for(
    hexane_activity: float,
    pressure_pa: float,
    k_wh: float,
    ideal_vapor: bool,
) -> tuple[float, float]:
    """The temperature bracket for the water residual.

    `hexane_activity == 1.0` returns the FROZEN constants unchanged, so
    the certified pin's bisection sequence — and hence its last bit — is
    preserved exactly.

    Otherwise the ceiling is SOLVED, not declared.  Note it is NOT simply
    `hexane_equality_ceiling_k` (where a_h·f_h^L = P): the fixed point
    needs y_h = a_h·f_h^L/(φ_h·P) < 1, and φ_h ≈ 0.99 puts the real edge
    BELOW that.  y_h(T) is monotone increasing (f_h^L rises steeply while
    φ_h barely moves), so the admissible set is an interval and its right
    edge is found by bisection on admissibility itself — no standoff
    constant to tune, and no silent extrapolation past it.
    """
    if hexane_activity == 1.0:
        return _BRACKET_LOW_K, _BRACKET_HIGH_K
    low = _BRACKET_LOW_K
    if not _hexane_equality_admissible(
        low, pressure_pa, k_wh, ideal_vapor, hexane_activity
    ):
        raise HeteroazeotropeError(
            f"hexane equality already inadmissible at the bracket floor "
            f"{low} K (a_h = {hexane_activity}, P = {pressure_pa} Pa)"
        )
    ceiling = _LIQUID_REFERENCE_MAX_K
    if _hexane_equality_admissible(
        ceiling, pressure_pa, k_wh, ideal_vapor, hexane_activity
    ):
        return low, ceiling
    # `_hexane_equality_admissible` is False for BOTH failure modes —
    # y_h leaving (0, 1) and the liquid-reference guard firing — and both
    # are monotone in T, so one bisection finds the edge.  Searched to
    # `_BRACKET_SEARCH_TOLERANCE_K`, not to the solve tolerance: this
    # only has to be a valid endpoint of the right residual sign, and
    # the root inside it is still converged to 1e-10 K by the main
    # bisection.  Resolving the endpoint to 1e-10 K instead doubled the
    # cost of the whole solve for no change in the answer.
    lo_ok, hi_bad = low, ceiling
    for _ in range(200):
        mid = 0.5 * (lo_ok + hi_bad)
        if _hexane_equality_admissible(
            mid, pressure_pa, k_wh, ideal_vapor, hexane_activity
        ):
            lo_ok = mid
        else:
            hi_bad = mid
        if hi_bad - lo_ok < _BRACKET_SEARCH_TOLERANCE_K:
            break
    if not lo_ok > low:
        raise HeteroazeotropeError(
            f"no admissible temperature bracket at a_h = {hexane_activity}, "
            f"P = {pressure_pa} Pa (hexane-equality ceiling {ceiling} K)"
        )
    return low, lo_ok


def solve_heteroazeotrope(
    pressure_pa: float = 101325.0,
    *,
    k_wh: float = binary_gas.K_WH_CENTRAL,
    ideal_vapor: bool = False,
    water_activity: float = 1.0,
    hexane_activity: float = 1.0,
) -> HeteroazeotropeState:
    """Solve the activity-scaled co-boiling point fail-closed.

    `water_activity < 1` generalizes the water equality to
    y_w·φ_w·P = a_w·f_w^L — the sorption-suppressed interface of the
    film layer's below-cap branch (grounding freeze §2.2, our own
    defended extension).

    `hexane_activity < 1` (T3b) does the same on the hexane side, for a
    surface below the PHY-046 critical loading X_c(T) where no mobile
    n-hexane phase exists and the escaping tendency is the PHY-043/048
    GAB sorbate activity.  The temperature bracket is then SOLVED from
    the authorities rather than declared (see `_bracket_for`).

    At a_w = a_h = 1 the solve is bit-identical to the frozen two-liquid
    production form — same bracket constants, same bisection sequence.
    """
    if not 0.0 < water_activity <= 1.0:
        raise HeteroazeotropeError(
            f"water activity {water_activity!r} outside (0, 1]"
        )
    if not 0.0 < hexane_activity <= 1.0:
        raise HeteroazeotropeError(
            f"hexane activity {hexane_activity!r} outside (0, 1]"
        )
    if not binary_gas.PROJECT_PRESSURE_MIN_PA <= pressure_pa <= binary_gas.PROJECT_PRESSURE_MAX_PA:
        raise HeteroazeotropeError(
            f"pressure {pressure_pa!r} Pa outside the frozen project band "
            f"[{binary_gas.PROJECT_PRESSURE_MIN_PA}, {binary_gas.PROJECT_PRESSURE_MAX_PA}]"
        )
    if not binary_gas.K_WH_MIN <= k_wh <= binary_gas.K_WH_MAX:
        raise HeteroazeotropeError(
            f"k_wh {k_wh!r} outside the frozen interval "
            f"[{binary_gas.K_WH_MIN}, {binary_gas.K_WH_MAX}]"
        )

    low, high = _bracket_for(hexane_activity, pressure_pa, k_wh, ideal_vapor)
    residual_low, _ = _water_residual(
        low, pressure_pa, k_wh, ideal_vapor, water_activity, hexane_activity
    )
    residual_high, _ = _water_residual(
        high, pressure_pa, k_wh, ideal_vapor, water_activity, hexane_activity
    )
    if residual_low is None or residual_high is None:
        raise HeteroazeotropeError("bracket endpoint left the (0, 1) composition domain")
    if residual_low * residual_high > 0.0:
        raise HeteroazeotropeError(
            f"no sign change in [{low}, {high}] K at P={pressure_pa} Pa "
            f"(a_w = {water_activity}, a_h = {hexane_activity}); the surface "
            "cannot reach saturation inside the authorities' liquid-reference "
            "domain — this is a SUB-SATURATED (non-boiling) interface and "
            "needs the diffusive closure, not a co-boiling pin"
        )

    for _ in range(_MAX_BISECTIONS):
        mid = 0.5 * (low + high)
        residual_mid, _ = _water_residual(
            mid, pressure_pa, k_wh, ideal_vapor, water_activity, hexane_activity
        )
        if residual_mid is None:
            raise HeteroazeotropeError("interior point left the composition domain")
        if residual_mid == 0.0 or high - low < _TEMPERATURE_TOLERANCE_K:
            break
        if residual_low * residual_mid < 0.0:
            high = mid
        else:
            low, residual_low = mid, residual_mid

    T = 0.5 * (low + high)
    residual, y_h = _water_residual(
        T, pressure_pa, k_wh, ideal_vapor, water_activity, hexane_activity
    )
    if y_h is None:
        raise HeteroazeotropeError("converged point left the composition domain")
    w_h = y_h * _hexane.M / (y_h * _hexane.M + (1.0 - y_h) * _water.M)
    return HeteroazeotropeState(
        temperature_k=T,
        y_hexane=y_h,
        hexane_mass_fraction=w_h,
        pressure_pa=pressure_pa,
        k_wh=k_wh,
        ideal_vapor=ideal_vapor,
        water_residual_pa=residual,
        water_activity=water_activity,
        hexane_activity=hexane_activity,
    )


# ---------------------------------------------------------------------------
# G5 CLOSURE (owner ruling G5-R1, 2026-08-22): the derived appearance-gate
# temperature T_het(P).
# ---------------------------------------------------------------------------
#: Ceiling of the G5-grounded pressure band: 1.10 atm, the top of the
#: literature-bracketed band of the G5 grounding record
#: (``docs/GT_PS2_G5_HETEROAZEOTROPE_PRESSURE_GROUNDING_RECORD_2026-08-22.md``
#: sections 3.2-3.3, where both independent derivations were checked at
#: 111 458 Pa: 337.238 K frozen ideal / 337.117 K Tsonopoulos curve).  Held
#: as the construction-exact product ``1.10 * 101_325.0`` per the
#: environment-fragile pin discipline, never a rounded decimal.  The band
#: FLOOR is the frozen project band floor ``binary_gas.PROJECT_PRESSURE_MIN_PA``
#: (101 000 Pa), which this module does not restate and does not widen: the
#: G5 record's 99 300 Pa plant-dome point (its F2) sits BELOW the frozen
#: floor and stays outside the grounded code domain until the floor itself
#: is re-ruled.
G5_GROUNDED_PRESSURE_MAX_PA = 1.10 * 101_325.0


@lru_cache(maxsize=4096)
def grounded_heteroazeotrope_temperature_k(pressure_pa: float) -> float:
    """T_het(P) from the certified production solver, on the G5-grounded band.

    G5 CLOSED 2026-08-22 (owner ruling G5-R1: "derived T_het instead of
    constant is great, approved"): the provisional Kemper-chart constant
    335.15 K is retired and the appearance/condensation mode-(ii) gate of
    both immiscible-pair appearance laws is THIS derived temperature at the
    layer's own pressure.  Authority:
    ``docs/GT_PS2_G5_HETEROAZEOTROPE_PRESSURE_GROUNDING_RECORD_2026-08-22.md``
    - certified anchor 334.4804 K at 101 325 Pa (record section 3.1), local
    slope ~2.7e-4 K/Pa (section 3.2); both are declared CHECK VALUES with
    bounds in the test suite, never bit-pins.

    CLOSURE-DOMAIN DESIGN CHOICE (record section 5, G5-R1): the ruling's
    off-solver leg - the ideal-immiscible closure for pressures below the
    certified solver's floor - is UNREACHABLE in code because the frozen
    project band floor ``binary_gas.PROJECT_PRESSURE_MIN_PA`` (101 000 Pa)
    coincides with the solver's own floor, so this function carries the
    CERTIFIED SOLVER ONLY and the ideal closure is deliberately omitted
    rather than implemented dormant.  Outside the grounded band
    [``binary_gas.PROJECT_PRESSURE_MIN_PA``, ``G5_GROUNDED_PRESSURE_MAX_PA``]
    the request is REFUSED, never clamped; inside it the certified solver's
    own guards still apply and still fail closed.

    Memoized (pure in ``pressure_pa``, same float object returned) because
    every appearance-law layer evaluation consults the gate and a march
    revisits one pressure constantly; this is a speed change only, exactly
    the `_pure_liquid_fugacity` pattern.
    """

    if not (
        binary_gas.PROJECT_PRESSURE_MIN_PA <= pressure_pa <= G5_GROUNDED_PRESSURE_MAX_PA
    ):
        raise HeteroazeotropeError(
            f"pressure {pressure_pa!r} Pa outside the G5-grounded heteroazeotrope "
            f"band [{binary_gas.PROJECT_PRESSURE_MIN_PA}, {G5_GROUNDED_PRESSURE_MAX_PA}] "
            "(GT_PS2_G5_HETEROAZEOTROPE_PRESSURE_GROUNDING_RECORD_2026-08-22 "
            "sections 3.2-3.3); the gate is underivable here - refuse, never clamp"
        )
    return solve_heteroazeotrope(pressure_pa).temperature_k
