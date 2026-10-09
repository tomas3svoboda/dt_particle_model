r"""Gate 1f activation seam for the resolved representative sphere.

This module intentionally implements only the lossless transition

    high-loading attached-film state -> fully wet radial state at z = 1.

It does not advance a partial wet/dry front.  While attached n-hexane is
present, the internal loading is the *live* ``X_c(T)`` and the excess remains
external.  The material temperature is the algebraic inverse of the combined
common-datum energy of that wet core plus attached liquid n-hexane.

The event locator uses independently prescribed outward n-hexane and total
energy rates.  Consequently the event time and temperature are coupled through
both the mass balance and ``X_c(T)``; no loading threshold is applied after an
uncoupled temperature update.  At the event, :func:`wet_core.activate` stores
the historical loading once on an exact ``z=1`` grid.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from dtdc_simulator.core2.particle import front
from dtdc_simulator.core2.particle._exact_memo import exact_memo
from dtdc_simulator.core2.particle import grid as spherical_grid
from dtdc_simulator.core2.particle import wet_core
from dtdc_simulator.core2.props import critical_volume as cv
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import sorption as sp


class ActivationNotReached(ValueError):
    """The prescribed high-loading trajectory has no activation in the step."""


@dataclass(frozen=True)
class SphereParams:
    """Geometry and constitutive template for one representative sphere."""

    particle_radius_m: float = 0.885e-3
    radial_cells: int = 12
    wet_core: wet_core.WetCoreParams = wet_core.WetCoreParams()

    def __post_init__(self) -> None:
        if not math.isfinite(self.particle_radius_m) or self.particle_radius_m <= 0.0:
            raise ValueError("sphere radius must be positive and finite")
        if isinstance(self.radial_cells, bool) or not isinstance(self.radial_cells, int):
            raise ValueError("sphere radial cell count must be an integer")
        if self.radial_cells < 1:
            raise ValueError("sphere needs at least one radial cell")

    @property
    def particle_volume_m3(self) -> float:
        return 4.0 / 3.0 * math.pi * self.particle_radius_m**3

    @property
    def dry_meal_mass_kg(self) -> float:
        return self.wet_core.rho_dm_p * self.particle_volume_m3

    def grid(self) -> spherical_grid.SphericalGrid:
        return spherical_grid.uniform_grid(self.radial_cells, self.particle_radius_m)


@dataclass(frozen=True)
class HighLoadingState:
    """Extensive film-active state; temperature is an algebraic caloric trace."""

    dry_meal_mass_kg: float
    oil_label_mass_kg: float
    retained_water_mass_kg: float
    internal_hexane_mass_kg: float
    attached_hexane_mass_kg: float
    total_internal_energy_j: float
    temperature_k: float

    def __post_init__(self) -> None:
        values = (
            self.dry_meal_mass_kg,
            self.oil_label_mass_kg,
            self.retained_water_mass_kg,
            self.internal_hexane_mass_kg,
            self.attached_hexane_mass_kg,
            self.total_internal_energy_j,
            self.temperature_k,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("high-loading state values must be finite")
        if self.dry_meal_mass_kg <= 0.0:
            raise ValueError("high-loading dry-meal mass must be positive")
        masses = (
            self.oil_label_mass_kg,
            self.retained_water_mass_kg,
            self.internal_hexane_mass_kg,
            self.attached_hexane_mass_kg,
        )
        if any(value < 0.0 for value in masses):
            raise ValueError("high-loading material inventories must be non-negative")

    @property
    def total_hexane_mass_kg(self) -> float:
        return self.internal_hexane_mass_kg + self.attached_hexane_mass_kg

    @property
    def total_hexane_loading(self) -> float:
        return self.total_hexane_mass_kg / self.dry_meal_mass_kg

    @property
    def retained_water_loading(self) -> float:
        return self.retained_water_mass_kg / self.dry_meal_mass_kg

    @property
    def oil_fraction(self) -> float:
        return self.oil_label_mass_kg / self.dry_meal_mass_kg


@dataclass(frozen=True)
class RadialActivationState:
    """Exact fully wet endpoint produced by the activation state map."""

    grid: spherical_grid.SphericalGrid
    coordinate: front.FrontCoordinate
    wet_inventory: wet_core.WetCoreInventory
    temperatures_k: tuple[float, ...]
    oil_fraction_labels: tuple[float, ...]
    retained_water_loadings: tuple[float, ...]
    dry_meal_mass_kg: float
    oil_label_mass_kg: float
    retained_water_mass_kg: float
    internal_hexane_mass_kg: float
    attached_hexane_mass_kg: float
    total_internal_energy_j: float
    temperature_k: float
    dry_domain: None = None

    def __post_init__(self) -> None:
        if self.coordinate.z != 1.0 or self.coordinate.regime != "fully_wet":
            raise ValueError("radial activation must start at exact z=1")
        if self.coordinate.particle_radius_m != self.grid.R:
            raise ValueError("activation coordinate and radial grid radii differ")
        if self.dry_domain is not None:
            raise ValueError("the exact z=1 activation endpoint has no dry domain")
        arrays = (
            self.wet_inventory.loadings,
            self.temperatures_k,
            self.oil_fraction_labels,
            self.retained_water_loadings,
        )
        if any(len(values) != self.grid.n for values in arrays):
            raise ValueError("radial activation arrays must match the wet grid")
        if any(value != self.temperature_k for value in self.temperatures_k):
            raise ValueError("radial activation temperature must be spatially uniform")
        if self.attached_hexane_mass_kg != 0.0:
            raise ValueError("radial activation cannot retain attached n-hexane")

    @property
    def total_hexane_mass_kg(self) -> float:
        return self.internal_hexane_mass_kg


@dataclass(frozen=True)
class ActivationConservationLedger:
    """Independent high-loading-minus-radial transition residuals."""

    dry_meal_residual_kg: float
    oil_label_residual_kg: float
    retained_water_residual_kg: float
    hexane_residual_kg: float
    energy_residual_j: float
    temperature_residual_k: float
    max_normalized_residual: float


@dataclass(frozen=True)
class ActivationEvent:
    """Located event, exact radial map, and both independent balance audits."""

    initial: HighLoadingState
    high_loading_at_event: HighLoadingState
    radial: RadialActivationState
    transition_ledger: ActivationConservationLedger
    elapsed_high_loading_s: float
    remaining_radial_time_s: float
    outward_hexane_mass_kg: float
    outward_energy_j: float
    high_loading_hexane_balance_residual_kg: float
    high_loading_energy_balance_residual_j: float


def critical_loading(temperature_k: float, p: SphereParams) -> float:
    """Live PHY-050 internal loading used only before/at activation."""
    return cv.critical_hexane_loading(
        temperature_k,
        p.wet_core.epsilon_p,
        p.wet_core.rho_dm_p,
        pressure_pa=p.wet_core.pressure_pa,
    )


@exact_memo(ambient=lambda: (hx.T_REF, hx.H_REF, hx.S_REF))
def combined_high_loading_specific_energy(
    temperature_k: float,
    total_hexane_loading: float,
    wet_params: wet_core.WetCoreParams,
) -> float:
    """Wet core at live ``X_c(T)`` plus attached-liquid energy, J/kg dry.

    Internal and attached n-hexane use the same Span--Wagner liquid state and
    common datum.  Thus moving material between those two locations introduces
    no latent impulse.
    """
    if not math.isfinite(total_hexane_loading) or total_hexane_loading < 0.0:
        raise ValueError("total n-hexane loading must be finite and non-negative")
    X_c = cv.critical_hexane_loading(
        temperature_k,
        wet_params.epsilon_p,
        wet_params.rho_dm_p,
        pressure_pa=wet_params.pressure_pa,
    )
    attached = total_hexane_loading - X_c
    if attached < 0.0:
        raise ValueError(
            "bulk-only feed is subcritical and cannot define a unique radial state; "
            "reject, do not clamp"
        )
    liquid_u = hx.state_Tp(temperature_k, wet_params.pressure_pa, "liquid").u_mass
    return wet_core.specific_energy(temperature_k, X_c, wet_params) + attached * liquid_u


def high_loading_temperature_from_energy(
    total_internal_energy_j: float,
    dry_meal_mass_kg: float,
    total_hexane_mass_kg: float,
    wet_params: wet_core.WetCoreParams,
) -> float:
    """Invert combined high-loading common-datum energy without clipping."""
    values = (total_internal_energy_j, dry_meal_mass_kg, total_hexane_mass_kg)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("high-loading caloric inputs must be finite")
    if dry_meal_mass_kg <= 0.0 or total_hexane_mass_kg < 0.0:
        raise ValueError("high-loading caloric masses are outside physical bounds")
    loading = total_hexane_mass_kg / dry_meal_mass_kg
    lo = _lowest_feasible_temperature(loading, wet_params)
    hi = wet_params.T_max

    def energy_at(temperature_k: float) -> float:
        return dry_meal_mass_kg * combined_high_loading_specific_energy(
            temperature_k, loading, wet_params
        )

    E_lo = energy_at(lo)
    E_hi = energy_at(hi)
    if E_hi <= E_lo:
        raise ValueError("high-loading caloric branch does not have positive capacity")
    energy_scale = _energy_scale(dry_meal_mass_kg, wet_params, hi)
    endpoint_tol = 2.0e-12 * energy_scale
    if (
        total_internal_energy_j < E_lo - endpoint_tol
        or total_internal_energy_j > E_hi + endpoint_tol
    ):
        raise ValueError(
            "high-loading energy lies outside the feasible caloric bracket; reject, do not clamp"
        )
    if abs(total_internal_energy_j - E_lo) <= endpoint_tol:
        return lo
    if abs(total_internal_energy_j - E_hi) <= endpoint_tol:
        return hi
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if energy_at(mid) < total_internal_energy_j:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-11:
            break
    return 0.5 * (lo + hi)


@exact_memo(ambient=lambda: (hx.T_REF, hx.H_REF, hx.S_REF))
def initialize_qualified_feed(
    p: SphereParams,
    *,
    total_hexane_loading: float,
    temperature_k: float | None = None,
    total_internal_energy_j: float | None = None,
    retained_water_loading: float | None = None,
    oil_fraction: float | None = None,
) -> HighLoadingState:
    """Initialize the approved bulk high-loading feed boundary (PHY-018).

    Temperature, energy, or both may be supplied.  When both are present they
    must describe the same material state.  A subcritical bulk-only input is
    rejected because total loading cannot determine a front and dry profile.
    """
    if temperature_k is None and total_internal_energy_j is None:
        raise ValueError("qualified feed requires material temperature or energy")
    if not math.isfinite(total_hexane_loading) or total_hexane_loading < 0.0:
        raise ValueError("feed n-hexane loading must be finite and non-negative")
    X_water = p.wet_core.X_water if retained_water_loading is None else retained_water_loading
    w_o = p.wet_core.w_o if oil_fraction is None else oil_fraction
    wet = _material_wet_params(p, X_water, w_o)
    M_dm = p.dry_meal_mass_kg
    M_h = M_dm * total_hexane_loading
    if temperature_k is None:
        assert total_internal_energy_j is not None
        temperature_k = high_loading_temperature_from_energy(
            total_internal_energy_j, M_dm, M_h, wet
        )
    if not math.isfinite(temperature_k) or not wet.T_min <= temperature_k <= wet.T_max:
        raise ValueError("feed material temperature is outside the caloric bracket")
    X_c = cv.critical_hexane_loading(
        temperature_k, wet.epsilon_p, wet.rho_dm_p, pressure_pa=wet.pressure_pa
    )
    if total_hexane_loading < X_c:
        raise ValueError(
            "bulk-only feed is subcritical and cannot define a unique radial state; "
            "reject, do not clamp"
        )
    calculated_energy = M_dm * combined_high_loading_specific_energy(
        temperature_k, total_hexane_loading, wet
    )
    if total_internal_energy_j is None:
        total_internal_energy_j = calculated_energy
    else:
        if not math.isfinite(total_internal_energy_j):
            raise ValueError("feed material energy must be finite")
        tolerance = 2.0e-10 * _energy_scale(M_dm, wet, temperature_k)
        if abs(total_internal_energy_j - calculated_energy) > tolerance:
            raise ValueError("supplied feed temperature and energy are inconsistent")
    state = HighLoadingState(
        dry_meal_mass_kg=M_dm,
        oil_label_mass_kg=M_dm * w_o,
        retained_water_mass_kg=M_dm * X_water,
        internal_hexane_mass_kg=M_dm * X_c,
        attached_hexane_mass_kg=M_dm * (total_hexane_loading - X_c),
        total_internal_energy_j=total_internal_energy_j,
        temperature_k=temperature_k,
    )
    _validate_high_loading_state(state, p)
    return state


def locate_activation_event(
    state: HighLoadingState,
    step_duration_s: float,
    outward_hexane_rate_kg_s: float,
    outward_total_energy_rate_w: float,
    p: SphereParams,
) -> ActivationEvent:
    """Locate high-loading depletion and return unused time for radial advance.

    Rates are constant over the event substep and use the outward-positive
    convention.  The n-hexane rate must be non-negative (no reverse drainage);
    the net total-energy rate is signed so inward conductive heating is allowed.
    The event solves

    ``M_h,a = M_h,0 - t_a m_dot = M_dm X_c(T_a)``
    ``U_a   = U_0   - t_a E_dot = M_dm e_w(T_a, X_c(T_a))``.
    """
    _validate_high_loading_state(state, p)
    rates = (step_duration_s, outward_hexane_rate_kg_s, outward_total_energy_rate_w)
    if not all(math.isfinite(value) for value in rates):
        raise ValueError("activation step and prescribed rates must be finite")
    if step_duration_s <= 0.0:
        raise ValueError("activation step duration must be positive")
    if outward_hexane_rate_kg_s < 0.0:
        raise ValueError("primary drainage forbids a negative outward n-hexane rate")

    if state.attached_hexane_mass_kg == 0.0:
        event_state = state
        elapsed = 0.0
    elif outward_hexane_rate_kg_s == 0.0:
        event_state, elapsed = _zero_mass_rate_event(
            state, step_duration_s, outward_total_energy_rate_w, p
        )
    else:
        event_state, elapsed = _coupled_rate_event(
            state,
            step_duration_s,
            outward_hexane_rate_kg_s,
            outward_total_energy_rate_w,
            p,
        )
    if not 0.0 <= elapsed <= step_duration_s:
        raise RuntimeError("activation locator returned time outside its step")
    remaining = step_duration_s - elapsed
    endpoint_tolerance = 64.0 * math.ulp(max(step_duration_s, 1.0))
    if remaining <= endpoint_tolerance:
        # A root a few ulps before the macro endpoint is the exact endpoint,
        # not a physical radial substep with an unbounded manufactured rate.
        elapsed = step_duration_s
        remaining = 0.0
    radial, transition = activate_radially(event_state, p)
    lost_hexane = outward_hexane_rate_kg_s * elapsed
    lost_energy = outward_total_energy_rate_w * elapsed
    return ActivationEvent(
        initial=state,
        high_loading_at_event=event_state,
        radial=radial,
        transition_ledger=transition,
        elapsed_high_loading_s=elapsed,
        remaining_radial_time_s=remaining,
        outward_hexane_mass_kg=lost_hexane,
        outward_energy_j=lost_energy,
        high_loading_hexane_balance_residual_kg=(
            event_state.total_hexane_mass_kg - state.total_hexane_mass_kg + lost_hexane
        ),
        high_loading_energy_balance_residual_j=(
            event_state.total_internal_energy_j - state.total_internal_energy_j + lost_energy
        ),
    )


def activate_radially(
    event_state: HighLoadingState, p: SphereParams
) -> tuple[RadialActivationState, ActivationConservationLedger]:
    """Apply the exact, lossless high-loading -> ``z=1`` radial state map."""
    _validate_high_loading_state(event_state, p)
    if event_state.attached_hexane_mass_kg != 0.0:
        raise ValueError("radial activation requires exactly zero attached n-hexane")
    wet = _wet_params_for_state(event_state, p)
    grid = p.grid()
    inventory = wet_core.activate(grid, event_state.temperature_k, wet)
    temperatures = (event_state.temperature_k,) * grid.n
    oil_labels = (event_state.oil_fraction,) * grid.n
    water_labels = (event_state.retained_water_loading,) * grid.n

    # These totals are deliberately re-integrated from the radial geometry and
    # shell labels; they are not copied from the high-loading state.
    M_dm = sum(volume * wet.rho_dm_p for volume in grid.volumes)
    L_o = sum(volume * wet.rho_dm_p * label for volume, label in zip(grid.volumes, oil_labels))
    M_w = sum(
        volume * wet.rho_dm_p * loading for volume, loading in zip(grid.volumes, water_labels)
    )
    M_h = wet_core.total_hexane(grid, inventory, wet)
    U = wet_core.total_energy(grid, temperatures, inventory, wet)
    radial = RadialActivationState(
        grid=grid,
        coordinate=front.FrontCoordinate(p.particle_radius_m, 1.0),
        wet_inventory=inventory,
        temperatures_k=temperatures,
        oil_fraction_labels=oil_labels,
        retained_water_loadings=water_labels,
        dry_meal_mass_kg=M_dm,
        oil_label_mass_kg=L_o,
        retained_water_mass_kg=M_w,
        internal_hexane_mass_kg=M_h,
        attached_hexane_mass_kg=0.0,
        total_internal_energy_j=U,
        temperature_k=event_state.temperature_k,
    )
    ledger = _activation_ledger(event_state, radial, wet)
    return radial, ledger


def _material_wet_params(
    p: SphereParams, retained_water_loading: float, oil_fraction: float
) -> wet_core.WetCoreParams:
    if not all(math.isfinite(value) for value in (retained_water_loading, oil_fraction)):
        raise ValueError("feed water and oil labels must be finite")
    water_bounds = sp.LuikovParams()
    if not water_bounds.W_ref <= retained_water_loading <= water_bounds.W_cap:
        raise ValueError(
            "feed retained water must lie on the qualified positive-moisture "
            f"Luikov branch [{water_bounds.W_ref}, {water_bounds.W_cap}] "
            "inside the evidence cap"
        )
    if not sp.W_O_MIN <= oil_fraction <= sp.W_O_MAX:
        raise ValueError("feed residual-oil fraction is outside the PHY-048 envelope")
    return replace(p.wet_core, X_water=retained_water_loading, w_o=oil_fraction)


def _wet_params_for_state(state: HighLoadingState, p: SphereParams) -> wet_core.WetCoreParams:
    return _material_wet_params(p, state.retained_water_loading, state.oil_fraction)


def _validate_high_loading_state(state: HighLoadingState, p: SphereParams) -> None:
    wet = _wet_params_for_state(state, p)
    mass_scale = max(p.dry_meal_mass_kg, 1.0e-300)
    if abs(state.dry_meal_mass_kg - p.dry_meal_mass_kg) > 2.0e-13 * mass_scale:
        raise ValueError("high-loading dry-meal mass is inconsistent with sphere geometry")
    if not wet.T_min <= state.temperature_k <= wet.T_max:
        raise ValueError("high-loading temperature is outside the caloric bracket")
    X_c = cv.critical_hexane_loading(
        state.temperature_k,
        wet.epsilon_p,
        wet.rho_dm_p,
        pressure_pa=wet.pressure_pa,
    )
    expected_internal = state.dry_meal_mass_kg * X_c
    if abs(state.internal_hexane_mass_kg - expected_internal) > 2.0e-13 * max(
        expected_internal, 1.0e-300
    ):
        raise ValueError("high-loading internal n-hexane is not the live X_c(T)")
    # Reconstruct from the live internal loading plus the independently stored
    # attached loading.  At the exact endpoint this preserves ``X=X_c`` despite
    # the harmless rounding in a multiply-then-divide extensive conversion.
    caloric_loading = X_c + state.attached_hexane_mass_kg / state.dry_meal_mass_kg
    expected_energy = state.dry_meal_mass_kg * combined_high_loading_specific_energy(
        state.temperature_k, caloric_loading, wet
    )
    tolerance = 2.0e-10 * _energy_scale(state.dry_meal_mass_kg, wet, state.temperature_k)
    if abs(state.total_internal_energy_j - expected_energy) > tolerance:
        raise ValueError("high-loading temperature is inconsistent with its conserved energy")


def _lowest_feasible_temperature(total_loading: float, wet: wet_core.WetCoreParams) -> float:
    X_lo = cv.critical_hexane_loading(
        wet.T_min, wet.epsilon_p, wet.rho_dm_p, pressure_pa=wet.pressure_pa
    )
    X_hi = cv.critical_hexane_loading(
        wet.T_max, wet.epsilon_p, wet.rho_dm_p, pressure_pa=wet.pressure_pa
    )
    if total_loading < X_hi:
        raise ValueError(
            "n-hexane loading is subcritical throughout the caloric bracket; "
            "reject, do not clamp"
        )
    if total_loading >= X_lo:
        return wet.T_min
    lo, hi = wet.T_min, wet.T_max
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        X_mid = cv.critical_hexane_loading(
            mid, wet.epsilon_p, wet.rho_dm_p, pressure_pa=wet.pressure_pa
        )
        if X_mid > total_loading:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-11:
            break
    return hi


def _coupled_rate_event(
    state: HighLoadingState,
    duration: float,
    mass_rate: float,
    energy_rate: float,
    p: SphereParams,
) -> tuple[HighLoadingState, float]:
    wet = _wet_params_for_state(state, p)
    M_dm = state.dry_meal_mass_kg

    def time_from_temperature(temperature_k: float) -> float:
        event_mass = M_dm * cv.critical_hexane_loading(
            temperature_k,
            wet.epsilon_p,
            wet.rho_dm_p,
            pressure_pa=wet.pressure_pa,
        )
        return (state.total_hexane_mass_kg - event_mass) / mass_rate

    t_min = time_from_temperature(wet.T_min)
    t_max = time_from_temperature(wet.T_max)
    if t_max < 0.0 or t_min > duration:
        raise ActivationNotReached("activation mass condition is outside the prescribed step")
    T_lo = wet.T_min if t_min >= 0.0 else _temperature_at_time(0.0, time_from_temperature, wet)
    T_hi = (
        wet.T_max
        if t_max <= duration
        else _temperature_at_time(duration, time_from_temperature, wet)
    )
    if T_lo > T_hi:
        raise ActivationNotReached("activation mass condition has an empty time bracket")

    def energy_residual(temperature_k: float) -> float:
        elapsed = time_from_temperature(temperature_k)
        X_c = cv.critical_hexane_loading(
            temperature_k,
            wet.epsilon_p,
            wet.rho_dm_p,
            pressure_pa=wet.pressure_pa,
        )
        event_energy = M_dm * combined_high_loading_specific_energy(temperature_k, X_c, wet)
        return state.total_internal_energy_j - energy_rate * elapsed - event_energy

    temperature = _first_bracketed_root(energy_residual, T_lo, T_hi)
    elapsed = time_from_temperature(temperature)
    if not 0.0 <= elapsed <= duration:
        raise ActivationNotReached("coupled activation root lies outside the prescribed step")
    X_a = cv.critical_hexane_loading(
        temperature, wet.epsilon_p, wet.rho_dm_p, pressure_pa=wet.pressure_pa
    )
    event_energy = M_dm * combined_high_loading_specific_energy(temperature, X_a, wet)
    return (
        HighLoadingState(
            dry_meal_mass_kg=M_dm,
            oil_label_mass_kg=state.oil_label_mass_kg,
            retained_water_mass_kg=state.retained_water_mass_kg,
            internal_hexane_mass_kg=M_dm * X_a,
            attached_hexane_mass_kg=0.0,
            total_internal_energy_j=event_energy,
            temperature_k=temperature,
        ),
        elapsed,
    )


def _zero_mass_rate_event(
    state: HighLoadingState,
    duration: float,
    energy_rate: float,
    p: SphereParams,
) -> tuple[HighLoadingState, float]:
    if energy_rate == 0.0:
        raise ActivationNotReached("zero prescribed losses cannot deplete attached n-hexane")
    wet = _wet_params_for_state(state, p)
    X_total = state.total_hexane_loading
    temperature = _temperature_for_critical_loading(X_total, wet)
    event_energy = state.dry_meal_mass_kg * combined_high_loading_specific_energy(
        temperature, X_total, wet
    )
    elapsed = (state.total_internal_energy_j - event_energy) / energy_rate
    if not 0.0 <= elapsed <= duration:
        raise ActivationNotReached("zero-mass-rate activation lies outside the prescribed step")
    return (
        HighLoadingState(
            dry_meal_mass_kg=state.dry_meal_mass_kg,
            oil_label_mass_kg=state.oil_label_mass_kg,
            retained_water_mass_kg=state.retained_water_mass_kg,
            internal_hexane_mass_kg=state.total_hexane_mass_kg,
            attached_hexane_mass_kg=0.0,
            total_internal_energy_j=event_energy,
            temperature_k=temperature,
        ),
        elapsed,
    )


def _temperature_for_critical_loading(loading: float, wet: wet_core.WetCoreParams) -> float:
    def residual(temperature_k: float) -> float:
        return (
            cv.critical_hexane_loading(
                temperature_k,
                wet.epsilon_p,
                wet.rho_dm_p,
                pressure_pa=wet.pressure_pa,
            )
            - loading
        )

    r_lo, r_hi = residual(wet.T_min), residual(wet.T_max)
    if r_lo == 0.0:
        return wet.T_min
    if r_hi == 0.0:
        return wet.T_max
    if r_lo * r_hi > 0.0:
        raise ActivationNotReached("constant-mass activation temperature is not bracketed")
    lo, hi = wet.T_min, wet.T_max
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-11:
            break
    return 0.5 * (lo + hi)


def _temperature_at_time(target: float, time_function, wet: wet_core.WetCoreParams) -> float:
    lo, hi = wet.T_min, wet.T_max
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if time_function(mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-11:
            break
    return 0.5 * (lo + hi)


def _first_bracketed_root(function, lo: float, hi: float) -> float:
    """Find the earliest smooth event root; fail closed if none is bracketed."""
    if lo == hi:
        if abs(function(lo)) <= 1.0e-12:
            return lo
        raise ActivationNotReached("coupled mass/energy activation root is not bracketed")
    subdivisions = 256
    left = lo
    f_left = function(left)
    if f_left == 0.0:
        return left
    bracket = None
    for index in range(1, subdivisions + 1):
        right = lo + (hi - lo) * index / subdivisions
        f_right = function(right)
        if f_right == 0.0:
            return right
        if f_left * f_right < 0.0:
            bracket = (left, right, f_left)
            break
        left, f_left = right, f_right
    if bracket is None:
        raise ActivationNotReached("coupled mass/energy activation root is not bracketed")
    left, right, f_left = bracket
    for _ in range(120):
        mid = 0.5 * (left + right)
        f_mid = function(mid)
        if f_mid == 0.0:
            return mid
        if f_left * f_mid <= 0.0:
            right = mid
        else:
            left, f_left = mid, f_mid
        if right - left < 1.0e-11:
            break
    return 0.5 * (left + right)


def _activation_ledger(
    before: HighLoadingState,
    after: RadialActivationState,
    wet: wet_core.WetCoreParams,
) -> ActivationConservationLedger:
    residuals = (
        after.dry_meal_mass_kg - before.dry_meal_mass_kg,
        after.oil_label_mass_kg - before.oil_label_mass_kg,
        after.retained_water_mass_kg - before.retained_water_mass_kg,
        after.total_hexane_mass_kg - before.total_hexane_mass_kg,
        after.total_internal_energy_j - before.total_internal_energy_j,
        after.temperature_k - before.temperature_k,
    )
    scales = (
        max(abs(before.dry_meal_mass_kg), 1.0e-300),
        max(abs(before.oil_label_mass_kg), 1.0e-300),
        max(abs(before.retained_water_mass_kg), 1.0e-300),
        max(abs(before.total_hexane_mass_kg), 1.0e-300),
        _energy_scale(before.dry_meal_mass_kg, wet, before.temperature_k),
        max(abs(before.temperature_k), 1.0),
    )
    return ActivationConservationLedger(
        dry_meal_residual_kg=residuals[0],
        oil_label_residual_kg=residuals[1],
        retained_water_residual_kg=residuals[2],
        hexane_residual_kg=residuals[3],
        energy_residual_j=residuals[4],
        temperature_residual_k=residuals[5],
        max_normalized_residual=max(
            abs(residual) / scale for residual, scale in zip(residuals, scales)
        ),
    )


def _energy_scale(
    dry_meal_mass_kg: float, wet: wet_core.WetCoreParams, temperature_k: float
) -> float:
    """Capacity scale independent of arbitrary component caloric datums."""
    cp_dm = sp.composite_dry_meal_heat_capacity(wet.w_o, wet.cp_dry_meal, wet.cp_oil)
    return dry_meal_mass_kg * cp_dm * max(abs(temperature_k), 1.0)
