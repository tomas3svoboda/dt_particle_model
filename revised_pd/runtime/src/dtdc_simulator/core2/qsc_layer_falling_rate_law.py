r"""B1 stage-3 falling-rate law — the receding-core arm, as pure algebra.

DESIGN AUTHORITY.  `docs/GT_PS2_B1_PHASE_EXHAUSTION_TRANSITION_DESIGN_PACKET_
2026-08-31.md` sections 3.3 (falling-rate arm) and 1b (method grounding),
signed by the owner in `docs/GT_PS2_B1_DESIGN_SIGNOFF_RECORD_2026-08-31.md`
as decision **D-B1-2**: "falling-rate arm = corrected Faner Eq. (9) receding
core with the banked GAB X_e floor (Cardarelli 1996), at L2; dual-scale
upgrade path reserved; Schluender-tension note attached."

WHAT THIS MODULE IS.  Four pure functions and two frozen records.  No I/O, no
module state, no clock, no ledger.  It owns the ALGEBRA of the falling-rate
arm and nothing else: the classifier (design 3.1), the exhaustion clearance
(3.2) and the tray-level wiring (4) are separate surfaces and are not touched
here.  Every authority value the law needs — the critical loading X_c, the
equilibrium floor X_e, the mass Biot number, the unhindered demand — is a
CALLER-SUPPLIED ARGUMENT.  The caller binds the authorities; this module only
evaluates the law.  The single exception is `equilibrium_floor_loading`, which
is an explicit delegating convenience onto the existing authority and computes
nothing of its own.

THE LAW.

Faner, Perez & Crapiste, J Food Process Eng 42 (2019) e12987, Eq. (7): during
the falling-rate period a sharp wet front recedes into the particle and its
radius follows

    D_i / D_p  =  s / R  =  ( (X - X_e) / (X_c - X_e) ) ** (1/3)

on the mobile-solvent interval between the equilibrium floor X_e and the
critical loading X_c.  KNOWN ERRATUM, carried forward rather than silently
corrected: Eq. (9)'s PRINTED falling-rate inequality is reversed.  The
physical interval is ``X_e < X < X_c`` — see
`docs/archive/GT_PS2_XC_XE_BASIS_AUDIT.md:36` ("the physical interval for the
falling-rate relation is X_e < X < X_c, despite the reversed printed
inequality in Eq. (9)"), and the same interval is what the live twin
`core2/sorption_interface.py` calls `HexaneRegime.RECEDING_FRONT`.  Both
endpoints are REFUSED here, each with its own type:

  * ``X <= X_e``  -> EXHAUSTED_TO_EQUILIBRIUM.  No liquid n-hexane remains
    anywhere; the sorbate/pore-diffusion arm owns that regime
    (`sorption_interface.HexaneRegime.DRY_SORBATE`).
  * ``X >= X_c``  -> ABOVE_CRITICAL.  A bulk liquid phase touches the surface;
    the FILM arm owns that regime (`HexaneRegime.FILM_ACTIVE`), and the B1
    design's own trigger (D-B1-1) is exact attached-phase exhaustion, not a
    fitted threshold.

THE RESISTANCE FACTOR, AND WHERE IT COMES FROM.  Design 3.3 specifies
"gas-side coupling through the existing BC-1 lumped surface law with the
receding-front resistance factor" but prints no closed form.  The physical
basis it points at is the quasi-steady spherical dry-shell conductance already
carried by `core2/sorption_interface.py:313-325`,

    k_shell = D_eff * s / (R * (R - s))     [m/s]

— zero at ``s = 0`` and singular at ``s = R``, i.e. film-controlled exactly
where the film arm takes over.  Putting that shell resistance IN SERIES with
the gas film of conductance ``k_gas`` gives the dimensionless multiplicative
factor on the gas-side conductance directly:

    factor = (1/k_gas) / (1/k_gas + 1/k_shell)
           = 1 / (1 + k_gas/k_shell)
           = 1 / (1 + Bi_m * (1 - f) / f),      f = s/R,

with the MASS BIOT NUMBER

    Bi_m = k_gas * R / D_eff

— the same group, with the same definition, that
`sorption_interface.surface_state_from_gas` already names in prose ("a mass
Biot number ``k_gas R / D_eff`` of order 1e5").  ``Bi_m`` is a DECLARED
caller-supplied argument here: this module does not reach for k_gas, R or
D_eff, so a caller may bind a plant Biot number, a fixture one, or a
sensitivity sweep without this law drifting from its own authority.

Consequences worth stating plainly rather than discovering later:
  * ``factor`` is EXACTLY 1.0 at ``f == 1.0`` (a zero-thickness shell has no
    resistance) and falls monotonically to 0 as the front recedes.
  * At the DT-scale ``Bi_m ~ 1e5`` quoted by the interface authority, the
    factor collapses hard as soon as the front leaves the surface: at
    ``f = 0.99`` it is already ~1e-3.  That is the receding-core mechanism
    doing its job, not a numerical pathology — but it means the arm is
    strongly Bi_m-sensitive and the caller's declared Bi_m is a first-order
    modelling choice, not a detail.

DISCLOSED, per the header of the authority this module delegates to.
`core2/sorption_interface.py` is a DISCLOSED NON-QUALIFYING closure
("physically_qualifying is False throughout"): its regime-B shell term has no
reachable transient oracle and is disclosed unqualified, its regime-B surface
activity ignores Kelvin depression (optimistic by up to ~15 % at a 10 nm
pore), it is isothermal across the shell, and it is monodisperse.  Everything
this module returns inherits that disclosure, so `physically_qualifying` is
False on every record here as well.  Design 3.3's own Schluender-2004 tension
note (receding fronts rejected for capillary-porous WATER drying) is resolved
in the signed justification record, not here.

NOT REPRESENTED, named rather than hidden:
  * no radial PDE — this is the algebraic tier D-B1-2 authorizes, with the
    Cardarelli dual-scale D_eff upgrade path explicitly reserved;
  * no temperature dependence of the factor beyond whatever the caller has
    already put into Bi_m, X_c and X_e;
  * no water arm — D-B1-4 makes the design phase-symmetric but implements
    hexane first;
  * no time integration: `evaluate_falling_rate_demand` factors ONE interval's
    already-computed unhindered demand.  Marching, clearance and conservation
    ledgers belong to the transition module (design 3.2/4).
"""

from __future__ import annotations

import enum
import math
from dataclasses import dataclass
from typing import ClassVar

from . import sorption_interface as _sorption_interface

FALLING_RATE_LAW_SCHEMA_ID = "GT-PS-2-B1-FALLING-RATE-LAW"
FALLING_RATE_LAW_SCHEMA_REVISION = 1
FALLING_RATE_LAW_SOURCE_IDENTITY = "dtdc_simulator.core2.qsc_layer_falling_rate_law"
FALLING_RATE_LAW_DESIGN_AUTHORITY = (
    "docs/GT_PS2_B1_PHASE_EXHAUSTION_TRANSITION_DESIGN_PACKET_2026-08-31.md sections 3.3, 1b"
)
FALLING_RATE_LAW_SIGNOFF_AUTHORITY = (
    "docs/GT_PS2_B1_DESIGN_SIGNOFF_RECORD_2026-08-31.md decision D-B1-2"
)

# COMPARISON ANCHOR ONLY (design 3.1 / owner decision D-B1-1): Faner's
# measured average critical loading for soybean and sunflower meals, carried
# as evidence.  It is NEVER read by the algebra below — X_c is an argument.
FANER_MEASURED_CRITICAL_LOADING_KG_KG = 0.20

# COMPARISON ANCHOR ONLY: the banked 373 K GAB SORBATE floor quoted in the
# `sorption_interface` header (lines 38-40, "X_e = retained_hexane(a_h=1) ...
# 0.0105 at 373 K").  That banked number is the SORBATE-ONLY loading; the live
# `equilibrium_floor_loading` below returns the TOTAL floor, which also carries
# the saturated pore GAS (~7.3 % more at 373 K).  The two are different
# quantities and must not be substituted for one another.
BANKED_373K_SORBATE_FLOOR_KG_KG = 0.0105

_ONE_THIRD = 1.0 / 3.0


class FallingRateRefusal(enum.Enum):
    """Every reason this law refuses, as a stable code."""

    NONFINITE_INPUT = "nonfinite_input"
    NONPHYSICAL_LOADING = "nonphysical_loading"
    DEGENERATE_LOADING_INTERVAL = "degenerate_loading_interval"
    EXHAUSTED_TO_EQUILIBRIUM = "exhausted_to_equilibrium"
    ABOVE_CRITICAL = "above_critical"
    FRONT_FRACTION_OUT_OF_RANGE = "front_fraction_out_of_range"
    NONPOSITIVE_MASS_BIOT = "nonpositive_mass_biot"
    NEGATIVE_DEMAND = "negative_demand"
    EQUILIBRIUM_FLOOR_UNAVAILABLE = "equilibrium_floor_unavailable"
    INCONSISTENT_DEMAND_RECORD = "inconsistent_demand_record"


class FallingRateLawError(ValueError):
    """A falling-rate contract was violated; the state is unusable."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.NONPHYSICAL_LOADING


class FallingRateNonFiniteError(FallingRateLawError):
    """An input was NaN or infinite; the law is defined on finite reals only."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.NONFINITE_INPUT


class FallingRateNonPhysicalLoadingError(FallingRateLawError):
    """A loading was negative; there is no such particle state."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.NONPHYSICAL_LOADING


class FallingRateDegenerateIntervalError(FallingRateLawError):
    """X_c is not strictly above X_e; there is no receding-front interval."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.DEGENERATE_LOADING_INTERVAL


class FallingRateExhaustedToEquilibriumError(FallingRateLawError):
    """X <= X_e: no liquid n-hexane remains; the sorbate arm owns this state."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.EXHAUSTED_TO_EQUILIBRIUM


class FallingRateAboveCriticalError(FallingRateLawError):
    """X >= X_c: bulk liquid touches the surface; the FILM arm owns this state."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.ABOVE_CRITICAL


class FallingRateFrontFractionError(FallingRateLawError):
    """The front fraction is outside the open-below/closed-above unit range."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.FRONT_FRACTION_OUT_OF_RANGE


class FallingRateMassBiotError(FallingRateLawError):
    """The declared mass Biot number is not strictly positive and finite."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.NONPOSITIVE_MASS_BIOT


class FallingRateDemandError(FallingRateLawError):
    """The unhindered demand is negative; this arm releases, it does not absorb."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.NEGATIVE_DEMAND


class FallingRateEquilibriumFloorError(FallingRateLawError):
    """The delegated sorption authority refused to supply X_e for this state."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.EQUILIBRIUM_FLOOR_UNAVAILABLE


class FallingRateRecordError(FallingRateLawError):
    """A demand record does not agree with its own inputs."""

    refusal: ClassVar[FallingRateRefusal] = FallingRateRefusal.INCONSISTENT_DEMAND_RECORD


def _as_finite_float(label: str, value: object) -> float:
    """Accept a real number, refuse bool/other types and non-finite values."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be a real number, not {type(value).__name__}")
    number = float(value)
    if not math.isfinite(number):
        raise FallingRateNonFiniteError(f"{label} must be finite, got {value!r}")
    return number


def _validated_interval(
    loading_kg_kg: object,
    critical_loading_kg_kg: object,
    equilibrium_loading_kg_kg: object,
) -> tuple[float, float, float]:
    """Refuse everything outside the corrected Faner interval ``X_e < X < X_c``."""
    loading = _as_finite_float("mean n-hexane loading X", loading_kg_kg)
    critical = _as_finite_float("critical loading X_c", critical_loading_kg_kg)
    floor = _as_finite_float("equilibrium floor X_e", equilibrium_loading_kg_kg)
    for label, value in (
        ("mean n-hexane loading X", loading),
        ("critical loading X_c", critical),
        ("equilibrium floor X_e", floor),
    ):
        if value < 0.0:
            raise FallingRateNonPhysicalLoadingError(f"{label} must be non-negative, got {value!r}")
    if not critical > floor:
        raise FallingRateDegenerateIntervalError(
            f"degenerate receding-front interval: X_e = {floor!r} is not strictly "
            f"below X_c = {critical!r}; there is no falling-rate arm to evaluate"
        )
    if loading <= floor:
        raise FallingRateExhaustedToEquilibriumError(
            f"X = {loading!r} kg/kg is at or below the equilibrium floor "
            f"X_e = {floor!r}: no liquid n-hexane remains, so the receding-front "
            "law does not apply — the sorbate/pore-diffusion arm owns this state"
        )
    if loading >= critical:
        raise FallingRateAboveCriticalError(
            f"X = {loading!r} kg/kg is at or above the critical loading "
            f"X_c = {critical!r}: a bulk liquid phase still touches the surface, "
            "so the FILM arm owns this state, not the falling-rate arm"
        )
    return loading, critical, floor


def _validated_mass_biot(mass_biot_number: object) -> float:
    biot = _as_finite_float("mass Biot number Bi_m", mass_biot_number)
    if biot <= 0.0:
        raise FallingRateMassBiotError(
            f"the declared mass Biot number Bi_m = k_gas*R/D_eff must be strictly "
            f"positive, got {biot!r}"
        )
    return biot


def _validated_front_fraction(front_fraction: object) -> float:
    fraction = _as_finite_float("front fraction s/R", front_fraction)
    if not 0.0 < fraction <= 1.0:
        raise FallingRateFrontFractionError(
            f"front fraction s/R must lie in (0, 1], got {fraction!r}: at s/R = 0 the "
            "core is exhausted (sorbate arm) and above 1 the front is outside the "
            "particle — refuse, do not clamp"
        )
    return fraction


def _validated_demand(unhindered_demand_kg: object) -> float:
    demand = _as_finite_float("unhindered demand", unhindered_demand_kg)
    if demand < 0.0:
        raise FallingRateDemandError(
            f"unhindered demand must be non-negative, got {demand!r} kg: the "
            "falling-rate arm releases n-hexane, it does not condense it — a "
            "condensing interval is a different regime and must be refused here"
        )
    return demand


# ---------------------------------------------------------------------------
# 1. Equilibrium floor X_e — delegated, never recomputed
# ---------------------------------------------------------------------------
def equilibrium_floor_loading(
    T_k: float, oil_mass_fraction: float, *, pressure_pa: float = 101325.0
) -> float:
    """X_e(T, w_o, P) in kg n-hexane / kg dry meal, from the existing authority.

    This is a THIN DELEGATION onto
    `core2.sorption_interface.equilibrium_loading_kg_kg`, which is the frozen
    definition of the floor:

        X_e = [eps_g * c_g(a*, T) + rho_dm,p * W_h(a*, T, w_o)] / rho_dm,p,
        a* = min(1, P / Psat_h(T))

    i.e. the pore GAS at the pressure-feasible activity ceiling plus the
    PHY-043/048 GAB sorbate at the same activity — the same inventory the
    certified `dry_shell.storage` counts.  Nothing is recomputed here and no
    constant is copied; only the authority's refusals are re-typed into this
    module's `FallingRateEquilibriumFloorError` so a caller can catch one
    exception family.

    F6 (2026-09-02).  ``pressure_pa`` is the BED pressure the floor is
    evaluated at and it defaults to 1 atm, which is bit-identical to the
    pre-F6 behaviour everywhere at or below the n-hexane boiling point (the
    ceiling is then exactly 1.0) — including the T3b 334.5 K anchor.  Above
    it the old implicit ``a_h = 1`` overstated the floor, and this floor is
    the DENOMINATOR of the corrected Faner relation below, so it moved the
    front trajectory and the residual tail.  A caller that knows its layer
    pressure must pass it.

    SIGNATURE NOTE (resolved, not silently chosen).  The authority's signature
    is ``equilibrium_loading_kg_kg(temperature_k, params=None)``; the oil mass
    fraction is not a positional argument there but the
    ``residual_oil_mass_fraction`` field of `ParticleSorptionParams`.  This
    wrapper therefore builds a `ParticleSorptionParams` carrying the requested
    ``oil_mass_fraction`` and the authority's OWN DEFAULTS for every other
    particle constant (radius, porosity, dry-meal density, pore diffusivity,
    GAB parameters, oil isotherm).  A caller who needs non-default particle
    constants must call the authority directly and pass X_e into
    `receding_front_fraction` as an argument — which is the normal path here,
    since the caller binds the authorities.

    DOMAIN.  Enforced by the authority, not duplicated here, so it cannot
    drift:
      * ``w_o`` must lie in the frozen PHY-048 envelope [0.012, 0.030] —
        outside it the composite retention closure REJECTS, it does not clamp;
      * ``T`` must be high enough that the GAB isotherm is non-singular at
        unit activity (``K(T) < 1``, i.e. above ~305.4 K for the frozen
        soybean parameters) and below the n-hexane critical temperature.

    DISCLOSED: `sorption_interface` is a non-qualifying closure
    (``physically_qualifying is False throughout``); this value inherits that
    standing in full.
    """
    temperature_k = _as_finite_float("temperature T", T_k)
    if temperature_k <= 0.0:
        raise FallingRateNonPhysicalLoadingError(
            f"temperature must be strictly positive, got {temperature_k!r} K"
        )
    oil_fraction = _as_finite_float("residual-oil mass fraction w_o", oil_mass_fraction)
    if not 0.0 < oil_fraction < 1.0:
        raise FallingRateNonPhysicalLoadingError(
            f"residual-oil mass fraction must lie strictly in (0, 1), got {oil_fraction!r}"
        )
    bed_pressure_pa = _as_finite_float("bed pressure P", pressure_pa)
    if bed_pressure_pa <= 0.0:
        raise FallingRateNonPhysicalLoadingError(
            f"bed pressure must be strictly positive, got {bed_pressure_pa!r} Pa"
        )
    params = _sorption_interface.ParticleSorptionParams(residual_oil_mass_fraction=oil_fraction)
    try:
        floor = _sorption_interface.equilibrium_loading_kg_kg(
            temperature_k, params, pressure_pa=bed_pressure_pa
        )
    except (
        _sorption_interface.SorptionInterfaceError,
        ValueError,
        ArithmeticError,
    ) as exc:
        raise FallingRateEquilibriumFloorError(
            f"the sorption authority refused the equilibrium floor at "
            f"T = {temperature_k!r} K, w_o = {oil_fraction!r}: {exc}"
        ) from exc
    if not math.isfinite(floor) or floor <= 0.0:
        raise FallingRateEquilibriumFloorError(
            f"the sorption authority returned a non-usable equilibrium floor "
            f"{floor!r} at T = {temperature_k!r} K, w_o = {oil_fraction!r}"
        )
    return floor


# ---------------------------------------------------------------------------
# 2. Receding-front position — corrected Faner Eq. (7)/(9)
# ---------------------------------------------------------------------------
def receding_front_fraction(X: float, X_c: float, X_e: float) -> float:
    """``D_i/D_p = s/R = ((X - X_e)/(X_c - X_e))**(1/3)`` on ``X_e < X < X_c``.

    Faner Eq. (7), evaluated on the CORRECTED Eq. (9) inequality (the printed
    one is reversed; `docs/archive/GT_PS2_XC_XE_BASIS_AUDIT.md:36`).  Both
    endpoints are refused with their own type — `FallingRateExhaustedTo
    EquilibriumError` at or below X_e, `FallingRateAboveCriticalError` at or
    above X_c — because each endpoint hands the state to a DIFFERENT arm, and
    a caller must be able to tell which.

    ``X_c`` and ``X_e`` are ARGUMENTS: the caller binds the authorities.  The
    live twins of this shape are
    `sorption_interface.front_radius_m` / `.classify` (which multiply the same
    cube root by the particle radius) and `props/critical_volume` (which
    supplies X_c(T)); this module deliberately imports neither of their values
    so that a fixture, a plant state and a sensitivity sweep all run the same
    algebra.

    The two differences are formed with `math.fsum`, which for two terms is
    bit-identical to the plain IEEE-754 subtraction the live twin uses — no
    rounding trick, no compensated re-association that would make this
    function disagree with `front_radius_m`.

    DECLARED ULP BEHAVIOUR AT THE UPPER ENDPOINT (measured, not clamped).  The
    open interval is exact in the INPUT, but the cube root compresses the
    output: at exactly ONE ULP below X_c the ratio rounds to 1.0 and so does
    its cube root, so this function can return exactly 1.0 for an X that is
    strictly inside the interval.  Measured width of that band: 1 ULP, i.e.
    ~2.8e-16 relative, on both the literature interval (0.0105, 0.20) and the
    live X_c(T)/X_e(T) pair.  This is the CORRECTLY ROUNDED value of the true
    expression, and the state it describes — front at the surface, zero shell
    resistance, film control — is exactly the limit the film arm carries, so
    the two arms meet continuously there.  It is therefore declared and gated,
    not surgically removed: refusing the band would open a hole below X_c that
    no arm owns, and rounding it away would be the clamp this module forbids.
    """
    loading, critical, floor = _validated_interval(X, X_c, X_e)
    mobile = math.fsum((loading, -floor))
    span = math.fsum((critical, -floor))
    return (mobile / span) ** _ONE_THIRD


# ---------------------------------------------------------------------------
# 3. Receding-front resistance factor on the gas-side conductance
# ---------------------------------------------------------------------------
def falling_rate_flux_factor(front_fraction: float, bi_m: float) -> float:
    """``1 / (1 + Bi_m * (1 - f) / f)`` — the dry shell in series with the film.

    The dimensionless variant of the shell resistance
    ``k_shell = D_eff*s/(R*(R - s))`` carried at
    `core2/sorption_interface.py:313-325`.  With a gas film of conductance
    ``k_gas`` in series,

        factor = (1/k_gas) / (1/k_gas + 1/k_shell) = 1 / (1 + Bi_m*(1 - f)/f)

    where ``f = s/R`` is the front fraction and ``Bi_m = k_gas*R/D_eff`` is
    the mass Biot number — DECLARED by the caller, exactly as the interface
    authority defines it in prose.  Multiply the BC-1 lumped gas-side
    conductance (or, equivalently, the unhindered demand it produces) by this
    number to get the receding-front-hindered value.

    Exactness that is relied upon elsewhere: at ``f == 1.0`` the hindrance term
    is ``Bi_m * 0.0 / 1.0 = 0.0`` and the factor is EXACTLY ``1.0`` for every
    admissible Bi_m — a zero-thickness shell has zero resistance, which is the
    same statement as `shell_conductance_m_s` returning ``+inf`` at ``s = R``.
    It is reached either by calling this function directly at ``f = 1.0`` or
    through the one-ULP band just below X_c where `receding_front_fraction`'s
    cube root rounds to 1.0 (declared in that function's docstring).

    Refusals: ``f`` outside ``(0, 1]`` and any non-positive or non-finite
    ``bi_m``.  ``f = 0`` is the exhausted core (infinite shell resistance) and
    belongs to the sorbate arm; it is refused, not returned as 0.0.
    """
    fraction = _validated_front_fraction(front_fraction)
    biot = _validated_mass_biot(bi_m)
    hindrance = biot * (math.fsum((1.0, -fraction)) / fraction)
    return 1.0 / math.fsum((1.0, hindrance))


# ---------------------------------------------------------------------------
# 4. Composition: one interval's demand, factored
# ---------------------------------------------------------------------------
@dataclass(frozen=True, slots=True, kw_only=True)
class FallingRateLawInputs:
    """Everything the falling-rate arm needs for ONE interval, all declared.

    No default is supplied for any field: each one is an authority the caller
    binds, and a silent default would be this module choosing physics.
    """

    loading_kg_kg: float
    critical_loading_kg_kg: float
    equilibrium_loading_kg_kg: float
    mass_biot_number: float
    unhindered_demand_kg: float

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _validated_interval(
            self.loading_kg_kg,
            self.critical_loading_kg_kg,
            self.equilibrium_loading_kg_kg,
        )
        _validated_mass_biot(self.mass_biot_number)
        _validated_demand(self.unhindered_demand_kg)


@dataclass(frozen=True, slots=True, kw_only=True)
class FallingRateDemand:
    """The factored demand plus the front state and every input echoed."""

    loading_kg_kg: float
    critical_loading_kg_kg: float
    equilibrium_loading_kg_kg: float
    mass_biot_number: float
    unhindered_demand_kg: float
    front_fraction: float
    flux_factor: float
    factored_demand_kg: float

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _validated_interval(
            self.loading_kg_kg,
            self.critical_loading_kg_kg,
            self.equilibrium_loading_kg_kg,
        )
        _validated_mass_biot(self.mass_biot_number)
        _validated_demand(self.unhindered_demand_kg)
        fraction = _validated_front_fraction(self.front_fraction)
        factor = _as_finite_float("flux factor", self.flux_factor)
        factored = _as_finite_float("factored demand", self.factored_demand_kg)
        if not 0.0 < factor <= 1.0:
            raise FallingRateRecordError(f"flux factor must lie in (0, 1], got {factor!r}")
        if fraction >= 1.0 and factor != 1.0:
            raise FallingRateRecordError(
                "a full-wetness front fraction must carry an exactly unit flux factor"
            )
        if factored != self.unhindered_demand_kg * factor:
            raise FallingRateRecordError(
                f"factored demand {factored!r} is not the exact product of the "
                f"unhindered demand {self.unhindered_demand_kg!r} and the flux "
                f"factor {factor!r}"
            )


def evaluate_falling_rate_demand(inputs: FallingRateLawInputs) -> FallingRateDemand:
    """Compose Faner Eq. (7) with the series-resistance factor for one interval.

    Pure algebra: the front fraction from `receding_front_fraction`, the
    hindrance factor from `falling_rate_flux_factor`, and their product with
    the caller's already-computed unhindered demand.  No time integration, no
    ledger, no state — the interval's clock, its conservation ledgers and the
    two-arm split itself belong to the transition module (design 3.2/4).

    Because `receding_front_fraction` refuses ``X >= X_c``, the returned factor
    is below 1 and the factored demand is below the unhindered demand — with
    ONE declared exception at binary64 resolution: inside the one-ULP band
    immediately below X_c the cube root rounds to exactly 1.0, the factor is
    then exactly 1.0, and the demand passes through unfactored.  That is the
    film-controlled limit the two arms share, measured at ~2.8e-16 relative
    width, and it is declared here rather than clamped away (see
    `receding_front_fraction`).
    """
    if type(inputs) is not FallingRateLawInputs:
        raise TypeError(
            f"falling-rate evaluation requires an exact FallingRateLawInputs, "
            f"not {type(inputs).__name__}"
        )
    inputs.__post_init__()
    fraction = receding_front_fraction(
        inputs.loading_kg_kg,
        inputs.critical_loading_kg_kg,
        inputs.equilibrium_loading_kg_kg,
    )
    factor = falling_rate_flux_factor(fraction, inputs.mass_biot_number)
    return FallingRateDemand(
        loading_kg_kg=float(inputs.loading_kg_kg),
        critical_loading_kg_kg=float(inputs.critical_loading_kg_kg),
        equilibrium_loading_kg_kg=float(inputs.equilibrium_loading_kg_kg),
        mass_biot_number=float(inputs.mass_biot_number),
        unhindered_demand_kg=float(inputs.unhindered_demand_kg),
        front_fraction=fraction,
        flux_factor=factor,
        factored_demand_kg=float(inputs.unhindered_demand_kg) * factor,
    )


__all__ = (
    "BANKED_373K_SORBATE_FLOOR_KG_KG",
    "FALLING_RATE_LAW_DESIGN_AUTHORITY",
    "FALLING_RATE_LAW_SCHEMA_ID",
    "FALLING_RATE_LAW_SCHEMA_REVISION",
    "FALLING_RATE_LAW_SIGNOFF_AUTHORITY",
    "FALLING_RATE_LAW_SOURCE_IDENTITY",
    "FANER_MEASURED_CRITICAL_LOADING_KG_KG",
    "FallingRateAboveCriticalError",
    "FallingRateDegenerateIntervalError",
    "FallingRateDemand",
    "FallingRateDemandError",
    "FallingRateEquilibriumFloorError",
    "FallingRateExhaustedToEquilibriumError",
    "FallingRateFrontFractionError",
    "FallingRateLawError",
    "FallingRateLawInputs",
    "FallingRateMassBiotError",
    "FallingRateNonFiniteError",
    "FallingRateNonPhysicalLoadingError",
    "FallingRateRecordError",
    "FallingRateRefusal",
    "equilibrium_floor_loading",
    "evaluate_falling_rate_demand",
    "falling_rate_flux_factor",
    "receding_front_fraction",
)
