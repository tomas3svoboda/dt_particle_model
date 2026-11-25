r"""The dry-limit continuation of the qualified positive-moisture Luikov branch.

Why this module exists
----------------------
The frozen modified-Luikov retained-water isotherm (PHY-031,
:func:`dtdc_simulator.core2.props.sorption.water_activity`) is qualified only on
the positive-moisture branch and refuses everything below its lower anchor:

    ``retained-water loading must be finite and >= the qualified Luikov lower
    bound 0.023; no low-moisture continuation is authorized``

and the feed entry refuses with

    ``feed retained water must lie on the qualified positive-moisture Luikov
    branch [0.023, 0.235670904] inside the evidence cap``

Faner's charge is **water-free**: there is no water in that experiment at all.
Every loading curve of Faner 2019 is therefore refused by the moisture branch,
and the refusal is recorded in
``docs/GT_PS2_PARTICLE_DRYING_VALIDATION_2026-09-03.md`` section 5.2.  The frozen
function's own docstring names the remedy exactly: *"A genuinely dry diagnostic
state therefore needs a separately named, non-qualifying closure; it must not
enter this authority as the convenient point ``(W_w, a_w) = (0, 0)``."*

This module is that separately named closure.  It edits nothing.

What the continuation is, and what is and is not singular
---------------------------------------------------------
The isotherm is ``a_w(W) = exp((1 - A1/W)/A2)``.  Two different objects must be
kept apart.

*   **The activity map is not singular at zero.**  As ``W -> 0+`` the exponent
    goes to ``-inf`` and ``a_w -> 0``.  The limit exists, it is zero, and every
    derivative of ``a_w`` with respect to ``W`` also goes to zero: the map
    extends to ``W = 0`` as a ``C^infinity`` function that is flat to all orders
    (an essential singularity of the *inverse*, not of the map).  Defining
    ``a_w(0) = 0`` is therefore the unique continuous continuation, not a choice.

*   **The state potential is singular.**  The retained-water force the transport
    law uses is ``ln(f_w,ret)``, which behaves as ``(1 - A1/W)/A2`` and diverges
    to ``-inf`` like ``-A1/(A2 W)``.  There is no finite limiting force, which is
    exactly why the frozen branch refuses rather than clamps, and why this
    continuation must not be fed back into the qualified authority.

The resolution is the same as for the single-component pore gas: in a water-free
charge the water component is **absent**, so the retained-water flux
``N_w,ret = -L_w,ret grad(mu_w,ret/(R T))`` is identically zero because there is
no water to carry, not because a divergent force was multiplied by a vanishing
mobility.  The potential is reported as ``-inf`` and typed as unusable.

A measured consequence
----------------------
In IEEE-754 double precision the exponential underflows to *exactly* ``0.0``
below a loading this module computes from the frozen parameters
(:func:`luikov_activity_underflow_loading`; about ``9.6e-5`` kg/kg dry at the
frozen ``A1``, ``A2``).  So a charge at ``1e-5`` or ``1e-6`` kg/kg dry and a
charge at exactly zero are **bit-identical** in activity: the continuation is not
an approximation there, it is the same number.  That is the strongest continuity
statement available and it is a measurement, not an argument.

Scope and honesty
-----------------
``physically_qualifying`` is ``False`` on every object here.  On the qualified
band ``[W_ref, W_cap]`` this module returns bit-identical values to the frozen
authority and is therefore a continuation of it, never a replacement; above the
evidence cap it refuses exactly as the frozen authority does.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache

from dtdc_simulator.core2.props import sorption as sp

#: The water-free charge, exactly.
WATER_FREE_CHARGE_LOADING_KG_KG_DRY = 0.0

#: The activity of the water-free charge, exactly.
DRY_LIMIT_WATER_ACTIVITY = 0.0

#: The loadings the qualification of this module is reported at.
QUALIFICATION_LOADINGS_KG_KG_DRY = (1.0e-4, 1.0e-5, 1.0e-6)


class DryLimitMoistureError(ValueError):
    """Raised when a loading is outside the continued branch."""


class DryLimitPotentialError(RuntimeError):
    """Raised when a caller asks the water-free charge for a water potential."""


@dataclass(frozen=True)
class WaterFreeChargeState:
    """One typed water-free charge: no water component is present.

    ``retained_water_potential_isothermal`` is ``-inf`` by construction and must
    never enter a force difference.  Ask
    :attr:`retained_water_flux_is_identically_zero` instead: it is the physically
    meaningful statement, and it is true by the absence of the component.
    """

    retained_water_loading_kg_kg_dry: float
    water_activity: float
    is_water_free: bool
    luikov: sp.LuikovParams
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def retained_water_potential_isothermal(self) -> float:
        """``ln(f_w,ret)`` of the charge; ``-inf`` when it is water free."""

        if self.is_water_free:
            return float("-inf")
        return (1.0 - self.luikov.A1 / self.retained_water_loading_kg_kg_dry) / self.luikov.A2

    @property
    def retained_water_flux_is_identically_zero(self) -> bool:
        """True for a water-free charge: there is no water to carry."""

        return self.is_water_free

    @property
    def is_inside_qualified_luikov_band(self) -> bool:
        """Whether the frozen positive-moisture authority admits this loading."""

        return (
            self.luikov.W_ref
            <= self.retained_water_loading_kg_kg_dry
            <= self.luikov.W_cap
        )

    @property
    def is_bit_identical_to_water_free(self) -> bool:
        """Whether this charge's activity is the water-free one, bit for bit."""

        return self.water_activity == DRY_LIMIT_WATER_ACTIVITY

    def require_finite_potential(self) -> float:
        """Return the retained-water potential, refusing the water-free charge."""

        if self.is_water_free:
            raise DryLimitPotentialError(
                "the water-free charge has no retained-water potential (-inf); the "
                "retained-water flux is identically zero by the absence of the "
                "component, and no force difference may be formed"
            )
        return self.retained_water_potential_isothermal


@dataclass(frozen=True)
class DryLimitRow:
    """One loading of the dry-limit continuity qualification."""

    retained_water_loading_kg_kg_dry: float
    water_activity: float
    absolute_activity_difference_from_water_free: float
    bit_identical_to_water_free: bool
    underflowed_to_zero: bool
    retained_water_potential_isothermal: float
    inside_qualified_luikov_band: bool
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class DryLimitQualification:
    """The measured approach of the continued branch to the water-free charge."""

    luikov: sp.LuikovParams
    water_free: WaterFreeChargeState
    rows: tuple[DryLimitRow, ...]
    underflow_loading_kg_kg_dry: float
    agrees_bitwise_with_frozen_authority_on_qualified_band: bool
    qualified_band_probe_count: int
    maximum_qualified_band_difference: float
    physically_qualifying: bool = field(default=False, init=False)


@lru_cache(maxsize=32)
def luikov_activity_underflow_loading(luikov: sp.LuikovParams = sp.LuikovParams()) -> float:
    """Largest loading whose Luikov activity is exactly ``0.0`` in double precision.

    Found by bisection on the frozen expression itself, so the value follows the
    frozen parameters rather than being a pinned literal.
    """

    lo, hi = 0.0, luikov.W_ref
    if _raw_activity(hi, luikov) == 0.0:
        raise DryLimitMoistureError(
            "the frozen Luikov activity underflows at its own qualified lower anchor"
        )
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        if _raw_activity(mid, luikov) == 0.0:
            lo = mid
        else:
            hi = mid
    return lo


def _raw_activity(loading_kg_kg_dry: float, luikov: sp.LuikovParams) -> float:
    """``exp((1 - A1/W)/A2)`` with the underflow left to the hardware."""

    if loading_kg_kg_dry <= 0.0:
        return 0.0
    exponent = (1.0 - luikov.A1 / loading_kg_kg_dry) / luikov.A2
    try:
        return math.exp(exponent)
    except OverflowError:  # pragma: no cover - unreachable below the cap
        raise DryLimitMoistureError("Luikov exponent overflowed") from None


def dry_limit_water_activity(
    loading_kg_kg_dry: float, luikov: sp.LuikovParams = sp.LuikovParams()
) -> float:
    """``a_w(W)`` continued down to and including ``W = 0``.

    Bit-identical to :func:`dtdc_simulator.core2.props.sorption.water_activity`
    everywhere that function is qualified.  Refuses above the evidence cap and
    refuses a negative or non-finite loading, exactly as the frozen branch does.
    """

    if not math.isfinite(loading_kg_kg_dry) or loading_kg_kg_dry < 0.0:
        raise DryLimitMoistureError(
            "retained-water loading must be finite and non-negative; the dry-limit "
            "continuation reaches zero, never below it"
        )
    if loading_kg_kg_dry > luikov.W_cap:
        raise DryLimitMoistureError(
            f"retained-water loading {loading_kg_kg_dry} exceeds the evidence cap "
            f"{luikov.W_cap}; apply the explicit active set, do not clamp"
        )
    if loading_kg_kg_dry == 0.0:
        return DRY_LIMIT_WATER_ACTIVITY
    return _raw_activity(loading_kg_kg_dry, luikov)


def dry_limit_water_activity_derivative(
    loading_kg_kg_dry: float, luikov: sp.LuikovParams = sp.LuikovParams()
) -> float:
    """``da_w/dW`` on the continued branch; exactly zero at ``W = 0``.

    The flatness at the origin is the reason the continuation is unique: the map
    and every one of its derivatives vanish there.
    """

    activity = dry_limit_water_activity(loading_kg_kg_dry, luikov)
    if loading_kg_kg_dry == 0.0:
        return 0.0
    return activity * luikov.A1 / (luikov.A2 * loading_kg_kg_dry * loading_kg_kg_dry)


def initialize_water_free_charge(
    luikov: sp.LuikovParams = sp.LuikovParams(),
) -> WaterFreeChargeState:
    """The typed feed at exactly zero retained water."""

    return WaterFreeChargeState(
        retained_water_loading_kg_kg_dry=WATER_FREE_CHARGE_LOADING_KG_KG_DRY,
        water_activity=DRY_LIMIT_WATER_ACTIVITY,
        is_water_free=True,
        luikov=luikov,
    )


def initialize_continued_charge(
    loading_kg_kg_dry: float, luikov: sp.LuikovParams = sp.LuikovParams()
) -> WaterFreeChargeState:
    """A typed charge anywhere on the continued branch, including zero."""

    activity = dry_limit_water_activity(loading_kg_kg_dry, luikov)
    return WaterFreeChargeState(
        retained_water_loading_kg_kg_dry=loading_kg_kg_dry,
        water_activity=activity,
        is_water_free=loading_kg_kg_dry == 0.0,
        luikov=luikov,
    )


def qualify_dry_limit(
    loadings_kg_kg_dry: tuple[float, ...] = QUALIFICATION_LOADINGS_KG_KG_DRY,
    luikov: sp.LuikovParams = sp.LuikovParams(),
    qualified_band_probes: int = 65,
) -> DryLimitQualification:
    """Measure the continuation against the water-free charge and the frozen branch.

    Two independent checks: the approach to zero at the requested small loadings,
    and bitwise agreement with the frozen authority at ``qualified_band_probes``
    points spread over ``[W_ref, W_cap]``.
    """

    water_free = initialize_water_free_charge(luikov)
    rows: list[DryLimitRow] = []
    for loading in loadings_kg_kg_dry:
        charge = initialize_continued_charge(loading, luikov)
        rows.append(
            DryLimitRow(
                retained_water_loading_kg_kg_dry=loading,
                water_activity=charge.water_activity,
                absolute_activity_difference_from_water_free=abs(
                    charge.water_activity - DRY_LIMIT_WATER_ACTIVITY
                ),
                bit_identical_to_water_free=charge.is_bit_identical_to_water_free,
                underflowed_to_zero=charge.water_activity == 0.0 and loading > 0.0,
                retained_water_potential_isothermal=(
                    charge.retained_water_potential_isothermal
                ),
                inside_qualified_luikov_band=charge.is_inside_qualified_luikov_band,
            )
        )

    if qualified_band_probes < 2:
        raise DryLimitMoistureError("the qualified-band comparison needs at least two probes")
    worst = 0.0
    agrees = True
    span = luikov.W_cap - luikov.W_ref
    for index in range(qualified_band_probes):
        probe = luikov.W_ref + span * index / (qualified_band_probes - 1)
        frozen = sp.water_activity(probe, luikov)
        continued = dry_limit_water_activity(probe, luikov)
        if continued != frozen:
            agrees = False
            worst = max(worst, abs(continued - frozen))

    return DryLimitQualification(
        luikov=luikov,
        water_free=water_free,
        rows=tuple(rows),
        underflow_loading_kg_kg_dry=luikov_activity_underflow_loading(luikov),
        agrees_bitwise_with_frozen_authority_on_qualified_band=agrees,
        qualified_band_probe_count=qualified_band_probes,
        maximum_qualified_band_difference=worst,
    )


__all__ = [
    "DRY_LIMIT_WATER_ACTIVITY",
    "QUALIFICATION_LOADINGS_KG_KG_DRY",
    "WATER_FREE_CHARGE_LOADING_KG_KG_DRY",
    "DryLimitMoistureError",
    "DryLimitPotentialError",
    "DryLimitQualification",
    "DryLimitRow",
    "WaterFreeChargeState",
    "dry_limit_water_activity",
    "dry_limit_water_activity_derivative",
    "initialize_continued_charge",
    "initialize_water_free_charge",
    "luikov_activity_underflow_loading",
    "qualify_dry_limit",
]
