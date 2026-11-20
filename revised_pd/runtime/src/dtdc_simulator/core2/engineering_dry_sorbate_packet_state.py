"""TAG-3 native particle state: all n-hexane sorbed, no mobile pore liquid.

Owner-ruled 2026-09-06 morning, verbatim "I rule all four as recommended,
proceed." (Q-C8-2B, Q-TAG2-FLOOR, Q-T2A-1, Q-D0-AUDIT-PIN); build
specification `GT_PS2_TAG3_DRY_SORBATE_PACKET_BUILD_SPEC_2026-09-06.md`
sections 1 and 6.

THE CALORIC LAW IS NOT NEW PHYSICS.  The frozen wet-core caloric
(`particle/wet_core.specific_energy`) writes a particle's energy per kg dry
meal as the dry-meal sensible energy, plus retained water, plus ALL hexane at
the liquid internal energy, minus the saturated binding deficit B_h(1, T),
which is the Clausius-Clapeyron integral of the frozen GAB isotherm from
activity 0 to activity 1 (`wet_core.saturated_binding_deficit`).  Its
`partition` gate refuses any loading below W_eq(a_h = 1, T) because the mobile
pore liquid has gone.  Tag-3 is that same law continued below that gate: all
hexane sorbed, at the partial activity a*(X, T) that inverts the same frozen
isotherm, with the SAME integral truncated at a*:

    e_tag3(T, X) = h_dm(T, w_o) + X_w h_w,l(T, P) + X u_h,l(T, P)
                   - B_h(a*(X, T), T),      0 < X <= W_eq(a_h=1, T)

At X == W_eq(1, T) exactly, a* is exactly 1.0 and the deficit is the FROZEN,
memoized `saturated_binding_deficit` object itself, so the tag-2 -> tag-3
seam is exact by construction (bit-for-bit, measured), exactly as the
tag-1 -> tag-2 seam is.

Typed endpoints, never clamped: a loading strictly above W_eq(1, T) has
mobile liquid and is tag-2 territory (`DrySorbateReWettingRefusal`); a
loading of exactly zero is the EXHAUSTED regime and carries no sorbate state
at all (`DrySorbateExhaustedRefusal`).

The pressure-feasible activity ceiling is deliberately NOT a state gate here.
A tag-3 particle may carry a*(X, T) above min(1, P / Psat_h(T)) - a
superheated sorbed state that will desorb fast.  The cap belongs where the
activity meets the gas, in the kernel binding, and is declared and recorded
there.

Declared engineering simplifications (engineering tier, stated not hidden):
- the pore-GAS hexane energy is excluded, exactly as the frozen wet-core law
  excludes it at and above the gate, so the seam stays exact;
- the sorbed fraction carries the truncated binding deficit of the same
  frozen isotherm, which is the definition of the band.

physically_qualifying: false.  plant_predictive: false.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from types import MappingProxyType
from typing import ClassVar

from dtdc_simulator.core2.engineering_subcritical_packet_state import (
    _BISECTION_ITERATIONS,
    _CONSISTENCY_RELATIVE_TOLERANCE,
    _ENDPOINT_RELATIVE_TOLERANCE,
    _TEMPERATURE_RESOLUTION_K,
)
from dtdc_simulator.core2.particle import sphere
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp

DRY_SORBATE_STATE_SOURCE_IDENTITY = "dtdc_simulator.core2.engineering_dry_sorbate_packet_state"
DRY_SORBATE_CALORIC_LAW_ID = "FROZEN_WET_CORE_CALORIC_CONTINUED_BELOW_THE_SORPTION_GATE_V1"

#: Provenance class of every declared value in this module (D0).
DECLARED_ENGINEERING_ASSUMPTION = "DECLARED_ENGINEERING_ASSUMPTION"

#: The ONLY new numeric literal in this module.  Declared below.
_ACTIVITY_RESOLUTION = 2.0**-48

DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "activity_resolution": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "the bisection resolution of the frozen tag-1 envelope inversion "
                    "(48 iterations), transported to the unit activity interval"
                ),
                "value": _ACTIVITY_RESOLUTION,
                "sources": (
                    "src/dtdc_simulator/core2/high_loading_packet_adapter.py:118-120 "
                    "(the measured 48 bisection iterations of the tag-1 combined "
                    "caloric inversion); "
                    "src/dtdc_simulator/core2/engineering_subcritical_packet_state.py"
                    ":51-52 (the tag-2 bisection contract this module imports)"
                ),
                "rejecting_outcome": (
                    "a seam handoff whose activity differs from 1.0 by more than "
                    "this while the loading equals the gate to the last bit"
                ),
            }
        ),
    }
)

#: The 24-panel Gauss-Legendre-5 panel count of the frozen activity
#: quadrature.  Read from ``wet_core._integrate_activity_0_1``'s own default;
#: this module never invents a rule.
_ACTIVITY_QUADRATURE_PANELS = 24


class DrySorbateStateError(ValueError):
    """A tag-3 native-state operand is malformed or inconsistent."""


class DrySorbateReWettingRefusal(DrySorbateStateError):
    """The loading is strictly above W_eq(a_h=1, T): tag-2 territory."""


class DrySorbateExhaustedRefusal(DrySorbateStateError):
    """The loading is not strictly positive: the EXHAUSTED regime."""


def _require_finite(name: str, value: float) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise DrySorbateStateError(f"{name} must be finite exact binary64")


@dataclass(frozen=True, slots=True, kw_only=True)
class DrySorbateParticleState:
    """Extensive dry-sorbate state; attached hexane is ZERO BY TYPE.

    There is deliberately no attached-hexane field and no mobile-liquid
    field: a tag-3 particle that re-wets (condensation, or heating past its
    own gate) must hand back to tag-2 through the codec boundary, never
    mutate in place.
    """

    dry_meal_mass_kg: float
    oil_label_mass_kg: float
    retained_water_mass_kg: float
    sorbed_hexane_mass_kg: float
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
            "sorbed_hexane_mass_kg",
            "total_internal_energy_j",
            "temperature_k",
        ):
            _require_finite(name, getattr(self, name))
        if self.dry_meal_mass_kg <= 0.0:
            raise DrySorbateStateError("dry-sorbate dry-meal mass must be positive")
        if self.sorbed_hexane_mass_kg <= 0.0:
            raise DrySorbateExhaustedRefusal(
                "dry-sorbate sorbed hexane must be positive; a hexane-free state "
                "is the EXHAUSTED regime and carries no sorbate state"
            )
        if self.oil_label_mass_kg < 0.0 or self.retained_water_mass_kg < 0.0:
            raise DrySorbateStateError("dry-sorbate component masses must be nonnegative")
        if self.temperature_k <= 0.0:
            raise DrySorbateStateError("dry-sorbate temperature must be positive")

    @property
    def total_hexane_mass_kg(self) -> float:
        return self.sorbed_hexane_mass_kg

    @property
    def sorbed_hexane_loading(self) -> float:
        return self.sorbed_hexane_mass_kg / self.dry_meal_mass_kg

    @property
    def retained_water_loading(self) -> float:
        return self.retained_water_mass_kg / self.dry_meal_mass_kg

    @property
    def oil_fraction(self) -> float:
        return self.oil_label_mass_kg / self.dry_meal_mass_kg


def sorption_gate_loading(temperature_k: float, wet: wet_core.WetCoreParams) -> float:
    """W_eq(a_h=1, T): the frozen partition's own mobile-liquid gate.

    Identical expression to ``engineering_subcritical_packet_state._sorption_floor``
    on the same ``WetCoreParams``; it is the tag-2 band's FLOOR and the tag-3
    band's CEILING, which is why the seam is one point and not an interval.
    """

    return sp.retained_hexane(1.0, temperature_k, wet.w_o, wet.gab, wet.oil)


def a_star(
    sorbed_hexane_loading: float, temperature_k: float, wet: wet_core.WetCoreParams
) -> float:
    """The partial n-hexane activity that inverts the frozen isotherm at X.

    The root of ``sp.retained_hexane(a, T, w_o, gab, oil) = X`` on (0, 1].
    The isotherm is monotone increasing in a, so a fixed-resolution bisection
    on the unit interval is exact to ``_ACTIVITY_RESOLUTION``.  At the gate
    (``X == W_eq(1, T)`` to the last bit) exactly ``1.0`` is returned with no
    bisection, which is what makes the seam value the FROZEN object.
    """

    _require_finite("dry-sorbate temperature", temperature_k)
    _require_finite("sorbed hexane loading", sorbed_hexane_loading)
    if sorbed_hexane_loading <= 0.0:
        raise DrySorbateExhaustedRefusal(
            "a non-positive sorbed loading is the EXHAUSTED regime; reject, do not clamp"
        )
    gate = sorption_gate_loading(temperature_k, wet)
    if sorbed_hexane_loading > gate:
        raise DrySorbateReWettingRefusal(
            "sorbed loading is strictly above the a_h=1 gate W_eq(T); mobile "
            "liquid is present and the tag-2 (subcritical) law owns this state"
        )
    if sorbed_hexane_loading == gate:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(_BISECTION_ITERATIONS):
        mid = 0.5 * (lo + hi)
        if sp.retained_hexane(mid, temperature_k, wet.w_o, wet.gab, wet.oil) < (
            sorbed_hexane_loading
        ):
            lo = mid
        else:
            hi = mid
        if hi - lo < _ACTIVITY_RESOLUTION:
            break
    return 0.5 * (lo + hi)


def partial_binding_deficit(
    activity: float, temperature_k: float, wet: wet_core.WetCoreParams
) -> float:
    """B_h(a*, T), J/kg dry meal: the frozen deficit truncated at a*.

    ``B = -R T^2 / M * integral_0^{a*} (dW/dT)_a / a da`` - the SAME identity
    ``wet_core.saturated_binding_deficit`` evaluates at a* = 1, with the SAME
    R, the SAME M and the SAME quadrature rule.  The nodes are scaled to
    [0, a*] by substituting ``a = a* u`` into the frozen
    ``wet_core._integrate_activity_0_1``, so the node set and the weights are
    literally the frozen ones and no rule is invented here.

    RULE: at ``activity == 1.0`` exactly this function CALLS the frozen,
    memoized ``wet_core.saturated_binding_deficit`` and never the partial
    integrator, so the seam value is the frozen object itself.
    """

    _require_finite("dry-sorbate activity", activity)
    if not 0.0 < activity <= 1.0:
        raise DrySorbateStateError(
            "the partial binding deficit is defined only for an activity in (0, 1]"
        )
    if activity == 1.0:
        return wet_core.saturated_binding_deficit(temperature_k, wet)
    integral = activity * wet_core._integrate_activity_0_1(
        lambda u: (
            sp.retained_hexane_dT(activity * u, temperature_k, wet.w_o, wet.gab, wet.oil)
            / (activity * u)
        ),
        _ACTIVITY_QUADRATURE_PANELS,
    )
    return -sp.R * temperature_k * temperature_k * integral / hx.M


def dry_sorbate_specific_energy(
    temperature_k: float, sorbed_hexane_loading: float, wet: wet_core.WetCoreParams
) -> float:
    """The frozen wet-core caloric continued below its gate, J/kg dry meal.

    Term for term the frozen ``wet_core.specific_energy`` expression, with the
    saturated binding deficit replaced by the SAME integral truncated at
    a*(X, T).  At X == W_eq(1, T) the two are bit-identical (asserted).
    """

    _require_finite("dry-sorbate temperature", temperature_k)
    _require_finite("sorbed hexane loading", sorbed_hexane_loading)
    activity = a_star(sorbed_hexane_loading, temperature_k, wet)
    try:
        u_h, _ = wet_core._hexane_liquid_u_du(temperature_k, wet.pressure_pa)
        h_w_retained, _ = wet_core._water_retained_h_dh(temperature_k, wet.pressure_pa)
    except ValueError as error:
        raise DrySorbateStateError(str(error)) from error
    return (
        sp.composite_dry_meal_sensible_energy(
            temperature_k, wet.T_ref_solid, wet.w_o, wet.cp_dry_meal, wet.cp_oil
        )
        + wet.X_water * h_w_retained
        + sorbed_hexane_loading * u_h
        - partial_binding_deficit(activity, temperature_k, wet)
    )


def dry_sorbate_temperature_bracket(
    sorbed_hexane_loading: float, wet: wet_core.WetCoreParams
) -> tuple[float, float]:
    """The feasible temperature band [T_lo, T_hi] for one sorbed loading.

    W_eq(a_h=1, T) FALLS with temperature, so a fixed sorbed loading stays
    inside the tag-3 band (X <= W_eq) while the particle is COOL and leaves it
    upward: heating shrinks the retained capacity until the excess re-appears
    as mobile liquid and the tag-2 law owns the state again.  The band is
    therefore ``[T_min, min(T_max, T_gate(X))]`` where ``T_gate`` is the root
    of ``W_eq(1, T) = X`` - the exact complement of the tag-2 module's own
    ``subcritical_temperature_bracket``, whose LOWER endpoint is the same
    crossing.  Both endpoints are located by the tag-2 bisection contract.

    (The build specification's section 1.6 states this inequality the other
    way round; the direction implemented here is the one the frozen isotherm
    and the tag-2 bracket both require, and it is measured in the bites.)
    """

    _require_finite("sorbed hexane loading", sorbed_hexane_loading)
    X = sorbed_hexane_loading
    if X <= 0.0:
        raise DrySorbateExhaustedRefusal(
            "a non-positive sorbed loading is the EXHAUSTED regime; reject, do not clamp"
        )
    if X > sorption_gate_loading(wet.T_min, wet):
        raise DrySorbateReWettingRefusal(
            "sorbed loading is above the a_h=1 gate throughout the caloric "
            "bracket; the tag-2 (subcritical) law owns it everywhere"
        )
    if X <= sorption_gate_loading(wet.T_max, wet):
        hi = wet.T_max
    else:
        lo_t, hi_t = wet.T_min, wet.T_max
        for _ in range(_BISECTION_ITERATIONS):
            mid = 0.5 * (lo_t + hi_t)
            if sorption_gate_loading(mid, wet) >= X:
                lo_t = mid
            else:
                hi_t = mid
            if hi_t - lo_t < _TEMPERATURE_RESOLUTION_K:
                break
        hi = lo_t  # strictly on the sorbed side of the crossing
    lo = wet.T_min
    if not lo < hi:
        raise DrySorbateStateError(
            "the sorbed loading has no feasible temperature band inside the "
            "caloric bracket; reject, do not clamp"
        )
    return lo, hi


def _energy_scale(
    dry_meal_mass_kg: float, wet: wet_core.WetCoreParams, temperature_k: float
) -> float:
    cp_dm = sp.composite_dry_meal_heat_capacity(wet.w_o, wet.cp_dry_meal, wet.cp_oil)
    return dry_meal_mass_kg * cp_dm * max(abs(temperature_k), 1.0)


def dry_sorbate_temperature_from_energy(
    total_internal_energy_j: float,
    dry_meal_mass_kg: float,
    sorbed_hexane_mass_kg: float,
    wet: wet_core.WetCoreParams,
) -> float:
    """Invert the dry-sorbate common-datum energy without clipping.

    Mirror of the tag-2 ``subcritical_temperature_from_energy`` contract
    (endpoint tolerance, fixed-resolution bisection) on the tag-3 band.
    """

    for name, value in (
        ("dry-sorbate energy", total_internal_energy_j),
        ("dry-sorbate dry mass", dry_meal_mass_kg),
        ("dry-sorbate hexane mass", sorbed_hexane_mass_kg),
    ):
        _require_finite(name, value)
    if dry_meal_mass_kg <= 0.0:
        raise DrySorbateStateError("dry-sorbate caloric masses are outside physical bounds")
    if sorbed_hexane_mass_kg <= 0.0:
        raise DrySorbateExhaustedRefusal(
            "a non-positive sorbed hexane mass is the EXHAUSTED regime; reject, do not clamp"
        )
    loading = sorbed_hexane_mass_kg / dry_meal_mass_kg
    lo, hi = dry_sorbate_temperature_bracket(loading, wet)

    def energy_at(temperature_k: float) -> float:
        return dry_meal_mass_kg * dry_sorbate_specific_energy(temperature_k, loading, wet)

    energy_lo = energy_at(lo)
    energy_hi = energy_at(hi)
    if energy_hi <= energy_lo:
        raise DrySorbateStateError("dry-sorbate caloric branch lacks positive capacity")
    endpoint_tol = _ENDPOINT_RELATIVE_TOLERANCE * _energy_scale(dry_meal_mass_kg, wet, hi)
    if (
        total_internal_energy_j < energy_lo - endpoint_tol
        or total_internal_energy_j > energy_hi + endpoint_tol
    ):
        raise DrySorbateStateError(
            "dry-sorbate energy lies outside the feasible caloric bracket; reject, do not clamp"
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


def initialize_dry_sorbate_state(
    p: sphere.SphereParams,
    *,
    sorbed_hexane_loading: float,
    temperature_k: float | None = None,
    total_internal_energy_j: float | None = None,
    retained_water_loading: float | None = None,
    oil_fraction: float | None = None,
) -> DrySorbateParticleState:
    """Initialize one canonical tag-3 native state (mirror of the tag-2 form).

    Temperature, energy, or both may be supplied; when both are present they
    must describe the same material state within the frozen consistency
    tolerance.  Material parameter gates (qualified Luikov moisture branch,
    PHY-048 oil envelope) are the SAME frozen gates tag-1 and tag-2 apply,
    reached through the same constructor.
    """

    if temperature_k is None and total_internal_energy_j is None:
        raise DrySorbateStateError("dry-sorbate state requires temperature or energy")
    _require_finite("sorbed hexane loading", sorbed_hexane_loading)
    if sorbed_hexane_loading <= 0.0:
        raise DrySorbateExhaustedRefusal(
            "dry-sorbate loading must be strictly positive; zero is the EXHAUSTED regime"
        )
    X_water = p.wet_core.X_water if retained_water_loading is None else retained_water_loading
    w_o = p.wet_core.w_o if oil_fraction is None else oil_fraction
    try:
        wet = sphere._material_wet_params(p, X_water, w_o)
    except ValueError as error:
        raise DrySorbateStateError(str(error)) from error
    dry_mass = p.dry_meal_mass_kg
    hexane_mass = dry_mass * sorbed_hexane_loading
    if temperature_k is None:
        assert total_internal_energy_j is not None
        temperature_k = dry_sorbate_temperature_from_energy(
            total_internal_energy_j, dry_mass, hexane_mass, wet
        )
    _require_finite("dry-sorbate temperature", temperature_k)
    lo, hi = dry_sorbate_temperature_bracket(sorbed_hexane_loading, wet)
    if not lo <= temperature_k <= hi:
        raise DrySorbateStateError(
            "dry-sorbate temperature is outside the feasible band for this loading"
        )
    calculated_energy = dry_mass * dry_sorbate_specific_energy(
        temperature_k, sorbed_hexane_loading, wet
    )
    if total_internal_energy_j is None:
        total_internal_energy_j = calculated_energy
    else:
        _require_finite("dry-sorbate energy", total_internal_energy_j)
        tolerance = _CONSISTENCY_RELATIVE_TOLERANCE * _energy_scale(dry_mass, wet, temperature_k)
        if abs(total_internal_energy_j - calculated_energy) > tolerance:
            raise DrySorbateStateError(
                "supplied dry-sorbate temperature and energy are inconsistent"
            )
    return DrySorbateParticleState(
        dry_meal_mass_kg=dry_mass,
        oil_label_mass_kg=dry_mass * w_o,
        retained_water_mass_kg=dry_mass * X_water,
        sorbed_hexane_mass_kg=hexane_mass,
        total_internal_energy_j=total_internal_energy_j,
        temperature_k=temperature_k,
    )


__all__ = (
    "DECLARED_ASSUMPTIONS",
    "DECLARED_ENGINEERING_ASSUMPTION",
    "DRY_SORBATE_CALORIC_LAW_ID",
    "DRY_SORBATE_STATE_SOURCE_IDENTITY",
    "DrySorbateExhaustedRefusal",
    "DrySorbateParticleState",
    "DrySorbateReWettingRefusal",
    "DrySorbateStateError",
    "a_star",
    "dry_sorbate_specific_energy",
    "dry_sorbate_temperature_bracket",
    "dry_sorbate_temperature_from_energy",
    "initialize_dry_sorbate_state",
    "partial_binding_deficit",
    "sorption_gate_loading",
)
