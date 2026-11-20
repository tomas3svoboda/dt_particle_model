"""Counter-current bound-water equations used by the DCZ."""

from __future__ import annotations

import math

from scipy.optimize import brentq

from dtdc_simulator.core import thermo


def countercurrent_vapor_water_profile(
    solid_moisture_profile: list[float],
    solid_moisture_in: float,
    dry_solid_flow_kg_s: float,
    bottom_vapor_water_kg_s: float,
) -> tuple[list[float], float]:
    """Reconstruct bottom-to-top vapor water from a solid-moisture profile."""
    vapor_water_kg_s = [0.0] * len(solid_moisture_profile)
    running_kg_s = bottom_vapor_water_kg_s
    minimum_kg_s = running_kg_s
    for cell in range(len(solid_moisture_profile) - 1, -1, -1):
        solid_water_in = (
            solid_moisture_in
            if cell == 0
            else solid_moisture_profile[cell - 1]
        )
        running_kg_s -= dry_solid_flow_kg_s * (
            solid_moisture_profile[cell] - solid_water_in
        )
        vapor_water_kg_s[cell] = running_kg_s
        minimum_kg_s = min(minimum_kg_s, running_kg_s)
    return vapor_water_kg_s, minimum_kg_s


def exact_linear_driving_force_step(
    value: float,
    equilibrium: float,
    rate_s_inv: float,
    dt_s: float,
) -> float:
    """Exact update of a local linear-driving-force equation."""
    if rate_s_inv < 0.0 or dt_s < 0.0:
        raise ValueError(
            "linear-driving-force rate and time step must be non-negative"
        )
    return equilibrium + (value - equilibrium) * math.exp(-rate_s_inv * dt_s)


def solve_bound_water_countercurrent_bvp(
    *,
    moisture_in: float,
    bottom_vapor_water_kg_s: float,
    vapor_hexane_top_kg_s: list[float],
    particle_surface_temperature_k: list[float],
    dry_solid_flow_kg_s: float,
    rate_s_inv: float,
    cell_residence_s: float | tuple[float, ...],
    bound_capacity: float,
    luikov: thermo.LuikovParams,
    antoine_water: thermo.AntoineParams,
    pressure_pa: float,
) -> tuple[list[float], list[float]]:
    """Solve the phase-free counter-current bound-water BVP by shooting."""
    if dry_solid_flow_kg_s <= 0.0 or bottom_vapor_water_kg_s <= 0.0:
        raise ValueError("counter-current water BVP requires positive flows")
    if len(vapor_hexane_top_kg_s) != len(particle_surface_temperature_k):
        raise ValueError("counter-current water BVP profiles must have equal length")
    if not vapor_hexane_top_kg_s:
        raise ValueError("counter-current water BVP requires at least one cell")
    residence_times = (
        tuple(cell_residence_s for _ in vapor_hexane_top_kg_s)
        if isinstance(cell_residence_s, (int, float))
        else cell_residence_s
    )
    if len(residence_times) != len(vapor_hexane_top_kg_s) or any(
        not math.isfinite(value) or value <= 0.0 for value in residence_times
    ):
        raise ValueError("counter-current water BVP residence times must be positive")

    def march(outlet_guess: float) -> tuple[list[float], list[float], float]:
        solid_profile: list[float] = []
        vapor_water_profile: list[float] = []
        solid_running = moisture_in
        for hexane_kg_s, surface_temperature_k, residence_time_s in zip(
            vapor_hexane_top_kg_s,
            particle_surface_temperature_k,
            residence_times,
        ):
            water_top_kg_s = bottom_vapor_water_kg_s - dry_solid_flow_kg_s * (
                outlet_guess - solid_running
            )
            if water_top_kg_s <= 0.0:
                return solid_profile, vapor_water_profile, math.nan
            mass_ratio_hexane_to_water = hexane_kg_s / water_top_kg_s
            activity = thermo.water_activity(
                mass_ratio_hexane_to_water,
                surface_temperature_k,
                antoine_water,
                P=pressure_pa,
            )
            activity = min(
                max(activity, 1.0e-9),
                thermo.LUIKOV_MAX_VALIDATED_UR,
            )
            equilibrium = thermo.luikov_equilibrium_moisture(activity, luikov)
            solid_running = exact_linear_driving_force_step(
                solid_running,
                equilibrium,
                rate_s_inv,
                residence_time_s,
            )
            solid_running = min(max(solid_running, 0.0), bound_capacity)
            solid_profile.append(solid_running)
            vapor_water_profile.append(water_top_kg_s)
        return solid_profile, vapor_water_profile, solid_running - outlet_guess

    lower = 0.0
    upper = min(
        bound_capacity,
        moisture_in
        + (bottom_vapor_water_kg_s - 1.0e-10) / dry_solid_flow_kg_s,
    )
    if upper <= lower:
        raise ValueError("counter-current water BVP has no positive outlet interval")
    _, _, residual_lower = march(lower)
    _, _, residual_upper = march(upper)
    if not math.isfinite(residual_lower):
        raise ValueError("counter-current water BVP lower bracket is non-finite")
    if not math.isfinite(residual_upper):
        feasible = lower
        infeasible = upper
        for _ in range(80):
            midpoint = 0.5 * (feasible + infeasible)
            _, _, midpoint_residual = march(midpoint)
            if math.isfinite(midpoint_residual):
                feasible = midpoint
            else:
                infeasible = midpoint
        upper = feasible
        _, _, residual_upper = march(upper)
    if not math.isfinite(residual_upper):
        raise ValueError(
            "counter-current water BVP has no locally feasible upper bracket"
        )
    if residual_lower == 0.0:
        outlet = lower
    elif residual_upper == 0.0:
        outlet = upper
    elif residual_lower * residual_upper > 0.0:
        raise ValueError(
            "counter-current water BVP has no bracketed outlet root "
            f"(lower={residual_lower:.6g}, upper={residual_upper:.6g})"
        )
    else:
        outlet = brentq(
            lambda value: march(value)[2],
            lower,
            upper,
            xtol=1.0e-12,
            rtol=1.0e-12,
        )
    solid_profile, vapor_water_profile, residual = march(outlet)
    if abs(residual) > 1.0e-10:
        raise ValueError(
            f"counter-current water BVP root residual is too large: {residual:.6g}"
        )
    return solid_profile, vapor_water_profile
