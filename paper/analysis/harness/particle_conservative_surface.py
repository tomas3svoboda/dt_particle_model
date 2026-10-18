"""Invert the existing native film law; do not add a surface constitutive law."""

import math

import numpy as np
from scipy.optimize import brentq

import finite_film_core2_study as study

core = study.core


def total_flux_for_hexane_target(
    film, target, surface_mass_fraction, bulk_mass_fraction, water_cp, hexane_cp, beta_limit
):
    """Unique native Stefan total flux matching a supplied hexane component flux.

    The component derivative lies between the two mass fractions. The heat
    beta derivative is consequently between min(cp)/h and max(cp)/h. Those
    bounds bracket every total flux inside the unchanged beta interval.
    """
    if not all(
        math.isfinite(x)
        for x in (
            target,
            surface_mass_fraction,
            bulk_mass_fraction,
            water_cp,
            hexane_cp,
            beta_limit,
        )
    ):
        raise ValueError("finite native film inverse inputs required")
    if min(water_cp, hexane_cp, beta_limit) <= 0:
        raise ValueError("positive component heat capacities and beta band required")
    if surface_mass_fraction == bulk_mass_fraction == 0:
        raise ValueError("hexane-free inverse is rank deficient; no total flux fabricated")
    calls = 0

    def evaluate(total):
        nonlocal calls
        calls += 1
        fluxes = core.efr.component_fluxes_for_supplied_total_mass_flux(
            film,
            total_mass_flux_kg_m2_s=total,
            interface_hexane_mass_fraction=surface_mass_fraction,
            bulk_hexane_mass_fraction=bulk_mass_fraction,
        )
        ack = core.efr.ackermann_heat_transfer(
            film,
            fluxes=fluxes,
            water_film_heat_capacity_j_kg_k=water_cp,
            hexane_film_heat_capacity_j_kg_k=hexane_cp,
            enforce_production_domain=False,
        )
        return fluxes, ack.beta

    zero_components, beta_zero = evaluate(0.0)
    if target == zero_components.hexane_kg_m2_s:
        if abs(beta_zero) > beta_limit:
            raise ValueError("zero-total component target exceeds the unchanged beta band")
        return 0.0, zero_components, beta_zero, calls
    bounds = [
        (b - beta_zero) * film.heat_transfer_w_m2_k / cp
        for b in (-beta_limit, beta_limit)
        for cp in (water_cp, hexane_cp)
    ]
    lower, upper = min(bounds), max(bounds)
    f_lower = evaluate(lower)[0].hexane_kg_m2_s - target
    f_upper = evaluate(upper)[0].hexane_kg_m2_s - target
    if f_lower > 0 or f_upper < 0:
        raise ValueError("component target has no native-film root inside the beta enclosure")
    total = brentq(
        lambda x: evaluate(x)[0].hexane_kg_m2_s - target,
        lower,
        upper,
        xtol=math.ulp(0.0),
        rtol=4 * np.finfo(float).eps,
        maxiter=100,
    )
    fluxes, beta = evaluate(total)
    if abs(beta) > beta_limit:
        raise ValueError("component target requires a film beta outside the unchanged band")
    return total, fluxes, beta, calls


def saturated_water_film_target(
    cell, config, radius, distance, boundary, surface_temperature, target_hexane
):
    """The existing water-trace endpoint plus the native component-flux inverse."""
    pore = config.dry.pore
    y = core.cp.decode_gas_only_y(surface_temperature, config.dry.pressure_pa, -32.0, pore)
    reference = core.ct.prepare_reduced_film_surface_state(
        boundary.film_boundary,
        cell_temperature_k=cell.temperature_k,
        pressure_pa=config.dry.pressure_pa,
        particle_radius_m=radius,
        outer_half_cell_distance_m=distance,
        particle_thermal_conductivity_w_m_k=config.dry.thermal_conductivity.value_w_m_k,
        binary_gas_interaction_k_wh=pore.k_wh,
    )
    film = core.ct._reduced_film_coefficients_at_surface_temperature(reference, surface_temperature)
    cp = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=surface_temperature,
        bulk_temperature_k=boundary.temperature_k,
        pressure_pa=config.dry.pressure_pa,
        y_hexane=y,
        k_wh=pore.k_wh,
    )
    total, components, beta, calls = total_flux_for_hexane_target(
        film,
        target_hexane,
        core._hexane_mass_fraction(y),
        core._hexane_mass_fraction(boundary.y_hexane),
        cp.water_j_kg_k,
        cp.hexane_j_kg_k,
        boundary.accepted_beta_limit,
    )
    return total, components, beta, y, calls


def _reduced_film(cell, config, radius, distance, boundary, surface_temperature):
    """The film coefficients exactly as the native water node assembles them."""
    reference = core.ct.prepare_reduced_film_surface_state(
        boundary.film_boundary,
        cell_temperature_k=cell.temperature_k,
        pressure_pa=config.dry.pressure_pa,
        particle_radius_m=radius,
        outer_half_cell_distance_m=distance,
        particle_thermal_conductivity_w_m_k=config.dry.thermal_conductivity.value_w_m_k,
        binary_gas_interaction_k_wh=config.dry.pore.k_wh,
    )
    return core.ct._reduced_film_coefficients_at_surface_temperature(reference, surface_temperature)


def interface_mass_fraction_for_component_fluxes(film, water_flux, hexane_flux, bulk_mass_fraction):
    """Invert the native Stefan partition for its interface mass fraction.

    With K=rho_f h_M, N=j_w+j_h and phi=Pe/expm1(Pe), the unchanged partition
    j_h = w_I N + K (w_I - w_b) phi is affine in w_I with slope
    N + K phi = K Pe exp(Pe)/expm1(Pe) > 0, so the inverse is unique and exact.
    """
    if not all(math.isfinite(x) for x in (water_flux, hexane_flux, bulk_mass_fraction)):
        raise ValueError("finite component fluxes and bulk composition required")
    conductance = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
    total = math.fsum((water_flux, hexane_flux))
    phi = core.efr._inverse_exprel(total / conductance) if total != 0.0 else 1.0
    slope = total + conductance * phi
    if not math.isfinite(slope) or slope <= 0:
        raise ValueError("native Stefan partition has no positive interface slope")
    fraction = math.fsum((hexane_flux, conductance * bulk_mass_fraction * phi)) / slope
    if not 0.0 < fraction < 1.0:
        raise ValueError("component fluxes need an interface composition outside (0, 1)")
    return total, fraction, slope


def absent_film_surface_for_component_targets(
    cell, config, radius, distance, boundary, surface_temperature, target_water, target_hexane
):
    """Gas-only surface: both component fluxes fixed by the particle inventories.

    The interface mole fraction follows from the exact affine inverse of the
    native partition (a few Newton corrections remove the mole/mass round trip).
    It must lie in the open gas-only chart and the resulting heat beta in the
    unchanged accepted band; otherwise the trial is refused, never clipped.
    """
    pore = config.dry.pore
    pressure = config.dry.pressure_pa
    film = _reduced_film(cell, config, radius, distance, boundary, surface_temperature)
    bulk = core._hexane_mass_fraction(boundary.y_hexane)
    total, fraction, slope = interface_mass_fraction_for_component_fluxes(
        film, target_water, target_hexane, bulk
    )
    mh, mw = core.M_HEXANE, core.M_WATER
    y = (fraction / mh) / (fraction / mh + (1.0 - fraction) / mw)
    corrections = 0
    for _ in range(4):
        components = core.efr.component_fluxes_for_supplied_total_mass_flux(
            film,
            total_mass_flux_kg_m2_s=total,
            interface_hexane_mass_fraction=core._hexane_mass_fraction(y),
            bulk_hexane_mass_fraction=bulk,
        )
        defect = components.hexane_kg_m2_s - target_hexane
        if defect == 0.0:
            break
        derivative = slope * mh * mw / (y * mh + (1.0 - y) * mw) ** 2
        candidate = y - defect / derivative
        if candidate == y:
            break
        y = candidate
        corrections += 1
    lo, hi = core.cp.gas_only_composition_interval(
        surface_temperature, pressure, pore
    ).y_hexane_bounds
    if not lo < y < hi:
        raise ValueError("absent-film component targets need a surface outside the gas-only chart")
    components = core.efr.component_fluxes_for_supplied_total_mass_flux(
        film,
        total_mass_flux_kg_m2_s=total,
        interface_hexane_mass_fraction=core._hexane_mass_fraction(y),
        bulk_hexane_mass_fraction=bulk,
    )
    cp = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=surface_temperature,
        bulk_temperature_k=boundary.temperature_k,
        pressure_pa=pressure,
        y_hexane=y,
        k_wh=pore.k_wh,
    )
    beta = core.efr.ackermann_heat_transfer(
        film,
        fluxes=components,
        water_film_heat_capacity_j_kg_k=cp.water_j_kg_k,
        hexane_film_heat_capacity_j_kg_k=cp.hexane_j_kg_k,
        enforce_production_domain=False,
    ).beta
    if abs(beta) > boundary.accepted_beta_limit:
        raise ValueError("absent-film component targets require a beta outside the unchanged band")
    return total, components, beta, y, corrections


def depletion_event_surface_temperature(
    old_liquid, particle_water_gain, duration, film_target, seed, lower, upper
):
    """Exact-zero film over a given duration: the water closure fixes T_s.

    For fixed cell inventories and duration the external hexane flux is fixed
    (the caller's film_target inverts the native saturated film for it at a
    trial surface temperature), and the surface temperature T solves
    g(T) = L_old - water_gain - d j_w(T) = 0. The root is bracketed outward
    from the caller's fixed seed inside the open solver temperature bounds,
    then the representable temperature with the smallest |g| within the final
    bracket is taken, so the closure is a pure function of the trial state.
    Returns (T, total, components, beta, y, inverse_calls, g(T), quantum), where
    quantum = |g(next float above T) - g(T)| is the closure's representable step.
    """
    if not all(math.isfinite(x) for x in (old_liquid, particle_water_gain, duration, seed)):
        raise ValueError("finite depletion-event inventories required")
    if not old_liquid > 0 or not duration > 0 or not lower < seed < upper:
        raise ValueError("depletion event needs positive liquid, duration and an admissible seed")
    calls = 0
    cache = {}

    def closure(temperature):
        nonlocal calls
        if temperature not in cache:
            total, components, beta, y, inverse = film_target(temperature)
            calls += inverse
            value = math.fsum(
                (old_liquid, -particle_water_gain, -duration * components.water_kg_m2_s)
            )
            cache[temperature] = (value, total, components, beta, y)
        return cache[temperature][0]

    first = closure(seed)
    root = seed
    if first != 0.0:
        other = None
        for step in (1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 1.0, 3.0):
            for trial in (seed - step, seed + step):
                if not lower < trial < upper:
                    continue
                try:
                    value = closure(trial)
                except ValueError:
                    continue
                if (value <= 0) != (first <= 0) or value == 0.0:
                    other = trial
                    break
            if other is not None:
                break
        if other is None:
            raise ValueError("no exact-zero film surface temperature brackets near the seed")
        lo, hi = sorted((seed, other))
        if closure(other) == 0.0:
            root = other
        else:
            root = brentq(
                closure, lo, hi, xtol=math.ulp(0.0), rtol=4 * np.finfo(float).eps, maxiter=100
            )
            candidates = [root]
            down = up = root
            for _ in range(4):
                down = float(np.nextafter(down, -np.inf))
                up = float(np.nextafter(up, np.inf))
                candidates.extend((down, up))
            admissible = []
            for candidate in candidates:
                if lo <= candidate <= hi:
                    try:
                        admissible.append((abs(closure(candidate)), candidate))
                    except ValueError:
                        continue
            root = min(admissible)[1]
    value, total, components, beta, y = cache[root]
    try:
        quantum = abs(closure(float(np.nextafter(root, np.inf))) - value)
    except ValueError:
        quantum = None
    return root, total, components, beta, y, calls, value, quantum
