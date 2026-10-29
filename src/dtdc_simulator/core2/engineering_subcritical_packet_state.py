"""TAG-2 native particle state: subcritical hexane, film exactly exhausted.

Owner-ruled GO 2026-09-01 (`GT_PS2_MORNING_RULING_RECORD_2026-09-01.md`,
disposition 1) opening B1 stage 3c: the additive second packet variant for
the marching band the tag-1 high-loading codec refuses (`total n-hexane is
subcritical throughout tag-1`).

THE CALORIC LAW IS NOT NEW PHYSICS.  In the subcritical marching band the
particle is a receding wet core: sorption sites saturated (a_h = 1) with the
remainder mobile pore liquid — exactly the state the FROZEN wet-core caloric
(`particle/wet_core.specific_energy`) already describes, and its own
`partition` gate fails closed at the sorption floor (mobile liquid gone).
This module therefore only CONTINUES that frozen law below pore saturation:

    e_tag2(T, X) = wet_core.specific_energy(T, X, wet)      for
    W_eq(a_h=1, T) <= X < X_c(T)

and is EXACTLY continuous with tag-1 at the variant boundary by
construction: tag-1's combined law is the same function evaluated at
X_c(T) plus an attached-liquid term that is exactly zero at the crossing —
no latent impulse can occur at a tag-1 -> tag-2 handoff.

Declared engineering simplifications (engineering tier, stated not hidden):
- the pore-GAS hexane energy in the emptied pore fraction is excluded (the
  tag-1 wet core excludes it identically at X_c, so the seam stays exact);
- the sorbed fraction carries the a_h = 1 saturated binding deficit of the
  frozen law throughout the band (exact while mobile liquid persists, which
  is the band's definition).

Typed endpoints, never clamped: X >= X_c(T) is tag-1 territory
(SUPERCRITICAL refusal), X below the sorption floor is the regime-C dry
continuation (SORPTION_FLOOR refusal — PART-02 machinery, not this module).

physically_qualifying: false.  plant_predictive: false.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import ClassVar

from dtdc_simulator.core2.particle import sphere
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.props import critical_volume as cv
from dtdc_simulator.core2.props import sorption as sp

SUBCRITICAL_STATE_SOURCE_IDENTITY = "dtdc_simulator.core2.engineering_subcritical_packet_state"
SUBCRITICAL_CALORIC_LAW_ID = "FROZEN_WET_CORE_CALORIC_CONTINUED_BELOW_PORE_SATURATION_V1"
#: Bisection contract mirrored from the frozen tag-1 inversions.
_BISECTION_ITERATIONS = 120
_TEMPERATURE_RESOLUTION_K = 1.0e-11
_ENDPOINT_RELATIVE_TOLERANCE = 2.0e-12
_CONSISTENCY_RELATIVE_TOLERANCE = 2.0e-10


class SubcriticalStateError(ValueError):
    """A tag-2 native-state operand is malformed or inconsistent."""


class SubcriticalSupercriticalRefusal(SubcriticalStateError):
    """The loading is at or above pore saturation: tag-1 territory."""


class SubcriticalSorptionFloorRefusal(SubcriticalStateError):
    """The loading is below the a_h=1 sorption floor: regime-C territory."""


def _require_finite(name: str, value: float) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise SubcriticalStateError(f"{name} must be finite exact binary64")


@dataclass(frozen=True, slots=True, kw_only=True)
class SubcriticalParticleState:
    """Extensive subcritical state; attached hexane is ZERO BY TYPE.

    There is deliberately no attached-hexane field: a tag-2 particle that
    re-grows a film (condensation) must hand back to tag-1 through the
    codec boundary, never mutate in place.
    """

    dry_meal_mass_kg: float
    oil_label_mass_kg: float
    retained_water_mass_kg: float
    internal_hexane_mass_kg: float
    total_internal_energy_j: float
    temperature_k: float

    attached_hexane_mass_kg: ClassVar[float] = 0.0
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        for name in (
            "dry_meal_mass_kg",
            "oil_label_mass_kg",
            "retained_water_mass_kg",
            "internal_hexane_mass_kg",
            "total_internal_energy_j",
            "temperature_k",
        ):
            _require_finite(name, getattr(self, name))
        if self.dry_meal_mass_kg <= 0.0:
            raise SubcriticalStateError("subcritical dry-meal mass must be positive")
        if self.internal_hexane_mass_kg <= 0.0:
            raise SubcriticalStateError(
                "subcritical internal hexane must be positive; a hexane-free "
                "state is outside every ruled tag-2 band"
            )
        if self.oil_label_mass_kg < 0.0 or self.retained_water_mass_kg < 0.0:
            raise SubcriticalStateError("subcritical component masses must be nonnegative")
        if self.temperature_k <= 0.0:
            raise SubcriticalStateError("subcritical temperature must be positive")

    @property
    def total_hexane_mass_kg(self) -> float:
        return self.internal_hexane_mass_kg

    @property
    def internal_hexane_loading(self) -> float:
        return self.internal_hexane_mass_kg / self.dry_meal_mass_kg

    @property
    def retained_water_loading(self) -> float:
        return self.retained_water_mass_kg / self.dry_meal_mass_kg

    @property
    def oil_fraction(self) -> float:
        return self.oil_label_mass_kg / self.dry_meal_mass_kg


def _critical_loading(temperature_k: float, wet: wet_core.WetCoreParams) -> float:
    return cv.critical_hexane_loading(
        temperature_k, wet.epsilon_p, wet.rho_dm_p, pressure_pa=wet.pressure_pa
    )


def _sorption_floor(temperature_k: float, wet: wet_core.WetCoreParams) -> float:
    """W_eq(a_h=1, T): the frozen partition's own mobile-liquid floor."""

    return sp.retained_hexane(1.0, temperature_k, wet.w_o, wet.gab, wet.oil)


def subcritical_specific_energy(
    temperature_k: float, internal_hexane_loading: float, wet: wet_core.WetCoreParams
) -> float:
    """The frozen wet-core caloric on the subcritical band, J/kg dry.

    Delegates entirely to ``wet_core.specific_energy`` (which fails closed
    at the sorption floor through its own ``partition`` gate) after refusing
    the supercritical side typed — this function adds a DOMAIN, not a law.
    """

    _require_finite("subcritical temperature", temperature_k)
    _require_finite("internal hexane loading", internal_hexane_loading)
    if internal_hexane_loading >= _critical_loading(temperature_k, wet):
        raise SubcriticalSupercriticalRefusal(
            "internal loading is at or above pore saturation X_c(T); the "
            "film-active tag-1 law owns this state"
        )
    try:
        return wet_core.specific_energy(temperature_k, internal_hexane_loading, wet)
    except ValueError as error:
        if "mobile liquid disappeared" in str(error):
            raise SubcriticalSorptionFloorRefusal(
                "internal loading is below the a_h=1 sorption floor W_eq(T); "
                "the regime-C dry continuation owns this state"
            ) from error
        raise SubcriticalStateError(str(error)) from error


def subcritical_temperature_bracket(
    internal_hexane_loading: float, wet: wet_core.WetCoreParams
) -> tuple[float, float]:
    """The feasible temperature band [T_lo, T_hi] for one subcritical loading.

    X_c(T) and W_eq(a_h=1, T) both fall with temperature, so the band is
    bounded above by the pore-saturation crossing (heating a subcritical
    particle re-saturates its pores — liquid expansion) and below by the
    sorption-floor crossing.  Both endpoints are located by the same
    bisection contract the frozen tag-1 inversion uses.
    """

    _require_finite("internal hexane loading", internal_hexane_loading)
    X = internal_hexane_loading
    if X >= _critical_loading(wet.T_min, wet):
        raise SubcriticalSupercriticalRefusal(
            "loading is supercritical throughout the caloric bracket; tag-1 territory"
        )
    if X < _sorption_floor(wet.T_max, wet):
        raise SubcriticalSorptionFloorRefusal(
            "loading is below the sorption floor throughout the caloric "
            "bracket; regime-C territory"
        )
    if X < _critical_loading(wet.T_max, wet):
        hi = wet.T_max
    else:
        lo_t, hi_t = wet.T_min, wet.T_max
        for _ in range(_BISECTION_ITERATIONS):
            mid = 0.5 * (lo_t + hi_t)
            if _critical_loading(mid, wet) > X:
                lo_t = mid
            else:
                hi_t = mid
            if hi_t - lo_t < _TEMPERATURE_RESOLUTION_K:
                break
        hi = lo_t  # strictly on the subcritical side of the crossing
    if X >= _sorption_floor(wet.T_min, wet):
        lo = wet.T_min
    else:
        lo_t, hi_t = wet.T_min, wet.T_max
        for _ in range(_BISECTION_ITERATIONS):
            mid = 0.5 * (lo_t + hi_t)
            if _sorption_floor(mid, wet) > X:
                lo_t = mid
            else:
                hi_t = mid
            if hi_t - lo_t < _TEMPERATURE_RESOLUTION_K:
                break
        lo = hi_t  # strictly on the feasible side of the floor
    if not lo < hi:
        raise SubcriticalStateError(
            "the subcritical loading has no feasible temperature band inside "
            "the caloric bracket; reject, do not clamp"
        )
    return lo, hi


def _energy_scale(
    dry_meal_mass_kg: float, wet: wet_core.WetCoreParams, temperature_k: float
) -> float:
    cp_dm = sp.composite_dry_meal_heat_capacity(wet.w_o, wet.cp_dry_meal, wet.cp_oil)
    return dry_meal_mass_kg * cp_dm * max(abs(temperature_k), 1.0)


def subcritical_temperature_from_energy(
    total_internal_energy_j: float,
    dry_meal_mass_kg: float,
    internal_hexane_mass_kg: float,
    wet: wet_core.WetCoreParams,
) -> float:
    """Invert the subcritical common-datum energy without clipping.

    Mirror of the frozen ``sphere.high_loading_temperature_from_energy``
    contract (monotone capacity guaranteed by the frozen positive-capacity
    gate; endpoint tolerance; fixed-resolution bisection) on the tag-2 band.
    """

    for name, value in (
        ("subcritical energy", total_internal_energy_j),
        ("subcritical dry mass", dry_meal_mass_kg),
        ("subcritical hexane mass", internal_hexane_mass_kg),
    ):
        _require_finite(name, value)
    if dry_meal_mass_kg <= 0.0 or internal_hexane_mass_kg <= 0.0:
        raise SubcriticalStateError("subcritical caloric masses are outside physical bounds")
    loading = internal_hexane_mass_kg / dry_meal_mass_kg
    lo, hi = subcritical_temperature_bracket(loading, wet)

    def energy_at(temperature_k: float) -> float:
        return dry_meal_mass_kg * subcritical_specific_energy(temperature_k, loading, wet)

    energy_lo = energy_at(lo)
    energy_hi = energy_at(hi)
    if energy_hi <= energy_lo:
        raise SubcriticalStateError("subcritical caloric branch lacks positive capacity")
    endpoint_tol = _ENDPOINT_RELATIVE_TOLERANCE * _energy_scale(dry_meal_mass_kg, wet, hi)
    if (
        total_internal_energy_j < energy_lo - endpoint_tol
        or total_internal_energy_j > energy_hi + endpoint_tol
    ):
        raise SubcriticalStateError(
            "subcritical energy lies outside the feasible caloric bracket; " "reject, do not clamp"
        )
    if abs(total_internal_energy_j - energy_lo) <= endpoint_tol:
        return lo
    if abs(total_internal_energy_j - energy_hi) <= endpoint_tol:
        return hi
    for _ in range(_BISECTION_ITERATIONS):
        mid = 0.5 * (lo + hi)
        if energy_at(mid) < total_internal_energy_j:
            lo = mid
        else:
            hi = mid
        if hi - lo < _TEMPERATURE_RESOLUTION_K:
            break
    return 0.5 * (lo + hi)


def initialize_subcritical_state(
    p: sphere.SphereParams,
    *,
    internal_hexane_loading: float,
    temperature_k: float | None = None,
    total_internal_energy_j: float | None = None,
    retained_water_loading: float | None = None,
    oil_fraction: float | None = None,
) -> SubcriticalParticleState:
    """Initialize one canonical tag-2 native state (mirror of PHY-018's form).

    Temperature, energy, or both may be supplied; when both are present they
    must describe the same material state within the frozen consistency
    tolerance.  Material parameter gates (qualified Luikov moisture branch,
    PHY-048 oil envelope) are the SAME frozen gates tag-1 applies, reached
    through the same constructor.
    """

    if temperature_k is None and total_internal_energy_j is None:
        raise SubcriticalStateError("subcritical state requires temperature or energy")
    _require_finite("internal hexane loading", internal_hexane_loading)
    if internal_hexane_loading <= 0.0:
        raise SubcriticalStateError("subcritical loading must be strictly positive")
    X_water = p.wet_core.X_water if retained_water_loading is None else retained_water_loading
    w_o = p.wet_core.w_o if oil_fraction is None else oil_fraction
    try:
        wet = sphere._material_wet_params(p, X_water, w_o)
    except ValueError as error:
        raise SubcriticalStateError(str(error)) from error
    dry_mass = p.dry_meal_mass_kg
    hexane_mass = dry_mass * internal_hexane_loading
    if temperature_k is None:
        assert total_internal_energy_j is not None
        temperature_k = subcritical_temperature_from_energy(
            total_internal_energy_j, dry_mass, hexane_mass, wet
        )
    _require_finite("subcritical temperature", temperature_k)
    lo, hi = subcritical_temperature_bracket(internal_hexane_loading, wet)
    if not lo <= temperature_k <= hi:
        raise SubcriticalStateError(
            "subcritical temperature is outside the feasible band for this loading"
        )
    calculated_energy = dry_mass * subcritical_specific_energy(
        temperature_k, internal_hexane_loading, wet
    )
    if total_internal_energy_j is None:
        total_internal_energy_j = calculated_energy
    else:
        _require_finite("subcritical energy", total_internal_energy_j)
        tolerance = _CONSISTENCY_RELATIVE_TOLERANCE * _energy_scale(dry_mass, wet, temperature_k)
        if abs(total_internal_energy_j - calculated_energy) > tolerance:
            raise SubcriticalStateError(
                "supplied subcritical temperature and energy are inconsistent"
            )
    return SubcriticalParticleState(
        dry_meal_mass_kg=dry_mass,
        oil_label_mass_kg=dry_mass * w_o,
        retained_water_mass_kg=dry_mass * X_water,
        internal_hexane_mass_kg=hexane_mass,
        total_internal_energy_j=total_internal_energy_j,
        temperature_k=temperature_k,
    )


__all__ = (
    "SUBCRITICAL_CALORIC_LAW_ID",
    "SUBCRITICAL_STATE_SOURCE_IDENTITY",
    "SubcriticalParticleState",
    "SubcriticalSorptionFloorRefusal",
    "SubcriticalStateError",
    "SubcriticalSupercriticalRefusal",
    "initialize_subcritical_state",
    "subcritical_specific_energy",
    "subcritical_temperature_bracket",
    "subcritical_temperature_from_energy",
)
