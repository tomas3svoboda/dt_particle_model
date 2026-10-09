"""Explicit numerical coordinates and pure-bulk film closure error budgets.

Physical residuals and stored states are unchanged. Endpoint closure uses a
prospective 128-epsilon arithmetic-scale budget, not a composition clamp.
"""

import math
import numpy as np
from scipy.optimize import OptimizeResult
from dtdc_simulator.core2.particle import external_film_reduced as efr


def domain_aware_newton(
    fun, x0, *, bounds, max_nfev=400, gtol=1e-11, difference_mode="forward", **unused
):
    """Backtrack undefined physical trials; never replace their residuals.

    Jacobians use finite differences with one-sided retreat when a probe is
    outside the physical domain. Every real residual call counts in nfev.
    """
    x = np.array(x0, dtype=float, copy=True)
    lower, upper = map(np.asarray, bounds)
    work = 0
    rejected = 0
    rejection_causes = {}
    if difference_mode not in ("forward", "central"):
        raise ValueError("unknown numerical derivative mode")

    def evaluate(values):
        nonlocal work, rejected
        if work >= max_nfev:
            raise RuntimeError("domain-aware Newton residual-call budget exhausted")
        work += 1
        try:
            result = np.asarray(fun(values), dtype=float)
            if not np.all(np.isfinite(result)):
                raise ValueError("physical residual is nonfinite")
            return result
        except (ValueError, RuntimeError) as exc:
            if "budget" in str(exc).lower() or "source" in str(exc).lower():
                raise
            rejected += 1
            label = f"{type(exc).__name__}: {exc}"
            rejection_causes[label] = rejection_causes.get(label, 0) + 1
            raise

    def jacobian(values, residual):
        columns = []
        for column in range(values.size):
            h = (
                np.cbrt(np.finfo(float).eps)
                if difference_mode == "central"
                else math.sqrt(np.finfo(float).eps)
            ) * max(1.0, abs(values[column]))
            for retreat in range(16):
                found = False
                samples = []
                for sign in (1.0, -1.0):
                    trial = values.copy()
                    trial[column] += sign * h
                    if not lower[column] < trial[column] < upper[column]:
                        continue
                    try:
                        sample = evaluate(trial)
                    except (ValueError, RuntimeError) as exc:
                        if "budget" in str(exc).lower() or "source" in str(exc).lower():
                            raise
                        continue
                    represented_h = trial[column] - values[column]
                    if represented_h == 0:
                        continue
                    samples.append((represented_h, sample))
                    if difference_mode == "forward":
                        break
                if len(samples) == 2:
                    columns.append(
                        (samples[0][1] - samples[1][1]) / (samples[0][0] - samples[1][0])
                    )
                    found = True
                elif len(samples) == 1 and (difference_mode == "forward" or retreat >= 4):
                    columns.append((samples[0][1] - residual) / samples[0][0])
                    found = True
                if found:
                    break
                h *= 0.5
            else:
                raise RuntimeError("no admissible one-sided physical Jacobian probe")
        return np.column_stack(columns)

    residual = evaluate(x)
    best = (float(np.max(np.abs(residual))), x.copy())
    try:
        for iteration in range(max_nfev):
            jac = jacobian(x, residual)
            maximum = float(np.max(np.abs(residual)))
            if maximum <= gtol:
                return OptimizeResult(
                    x=x,
                    fun=residual,
                    jac=jac,
                    success=True,
                    active_mask=np.zeros(x.size),
                    nfev=work,
                    optimality=float(np.max(np.abs(jac.T @ residual))),
                    message="domain-aware Newton converged",
                    rejected_trial_evaluations=rejected,
                )
            column_scale = 1 / np.maximum(np.linalg.norm(jac, axis=0), np.finfo(float).tiny)
            scaled_jac = jac * column_scale
            old_norm = np.linalg.norm(residual)
            advanced = False
            for damping in (0.0, 1e-6, 1e-4, 0.01, 1.0, 100.0):
                matrix = np.vstack((scaled_jac, math.sqrt(damping) * np.eye(x.size)))
                rhs = np.r_[-residual, np.zeros(x.size)]
                step = np.linalg.lstsq(matrix, rhs, rcond=1e-14)[0] * column_scale
                # Tangent search directions near a bound. This changes a
                # search direction only; no stored coordinate is projected.
                numerical_distance = 32 * np.finfo(float).eps * np.maximum(1.0, np.abs(x))
                blocked = ((step > 0) & ((upper - x) < numerical_distance)) | (
                    (step < 0) & ((x - lower) < numerical_distance)
                )
                if np.any(blocked):
                    allowed = ~blocked
                    step[:] = 0.0
                    step[allowed] = (
                        np.linalg.lstsq(matrix[:, allowed], rhs, rcond=1e-14)[0]
                        * column_scale[allowed]
                    )
                alpha = 1.0
                for index, change in enumerate(step):
                    if change > 0:
                        alpha = min(alpha, 0.99 * (upper[index] - x[index]) / change)
                    elif change < 0:
                        alpha = min(alpha, 0.99 * (lower[index] - x[index]) / change)
                for retreat in range(20):
                    trial = x + alpha * step
                    try:
                        candidate = evaluate(trial)
                    except (ValueError, RuntimeError) as exc:
                        if "budget" in str(exc).lower() or "source" in str(exc).lower():
                            raise
                        alpha *= 0.5
                        continue
                    if np.linalg.norm(candidate) < (1 - 1e-4 * alpha) * old_norm:
                        x, residual = trial, candidate
                        if np.max(np.abs(residual)) < best[0]:
                            best = (float(np.max(np.abs(residual))), x.copy())
                        advanced = True
                        break
                    alpha *= 0.5
                if advanced:
                    break
            if not advanced:
                return OptimizeResult(
                    x=x,
                    fun=residual,
                    jac=jac,
                    success=False,
                    active_mask=np.zeros(x.size),
                    nfev=work,
                    optimality=float(np.max(np.abs(jac.T @ residual))),
                    message="domain-aware Newton has no admissible decreasing step",
                    rejected_trial_evaluations=rejected,
                )
    except Exception as exc:
        exc.optimizer_debug = {
            "best_maximum_scaled_residual": best[0],
            "best_coordinates": best[1].tolist(),
            "residual_calls": work,
            "rejected_trials": rejected,
            "rejection_causes": rejection_causes,
        }
        raise
    raise RuntimeError("domain-aware Newton iteration limit exhausted")


def physical_chart_solve(
    fun, x0, *, decode, encode, physical_bounds, temperature_indices, **kwargs
):
    """Pull back the unchanged native residual into physical primitives.

    The encode callback still enforces every original native solver bound.
    Temperature origins are translations, with unit slopes. The returned
    Jacobian is explicitly in physical coordinates, for a conditioning audit.
    """
    initial = np.asarray(decode(x0), dtype=float)
    origin = np.zeros(initial.size)
    origin[list(temperature_indices)] = initial[list(temperature_indices)]
    lower, upper = (np.asarray(v, dtype=float) for v in physical_bounds)
    native_lower, native_upper = (np.asarray(v) for v in kwargs.pop("bounds"))

    def native(values):
        encoded = np.asarray(encode(values), dtype=float)
        if np.any(encoded <= native_lower) or np.any(encoded >= native_upper):
            raise ValueError("physical trial left the original open native solver bounds")
        return encoded

    def residual(values):
        return fun(native(values + origin))

    try:
        result = domain_aware_newton(
            residual, initial - origin, bounds=(lower - origin, upper - origin), **kwargs
        )
    except Exception as exc:
        debug = getattr(exc, "optimizer_debug", None)
        if debug is not None:
            physical = np.asarray(debug["best_coordinates"]) + origin
            debug["temperature_origins_k"] = origin.tolist()
            debug["best_physical_coordinates"] = physical.tolist()
            debug["best_native_coordinates"] = native(physical).tolist()
        raise
    result.physical_coordinates = result.x + origin
    result.x = native(result.physical_coordinates)
    if result.success:
        # Restore the native coordinate Jacobian for its original condition
        # gate. Physical increments remain the nonlinear search coordinates.
        # A raw matrix in K, kg/kg and seconds is not the native dimensionless
        # barrier matrix, particularly near a vanishing event duration.
        pullback = []
        for column in range(result.x.size):
            h = 1e-5 * max(1.0, abs(result.x[column]))
            for retreat in range(16):
                left = result.x.copy()
                right = result.x.copy()
                left[column] -= h
                right[column] += h
                try:
                    if np.any(left <= native_lower) or np.any(right >= native_upper):
                        raise ValueError("Jacobian pullback probe outside native bounds")
                    change = (np.asarray(decode(right)) - np.asarray(decode(left))) / (
                        right[column] - left[column]
                    )
                    pullback.append(change)
                    break
                except ValueError, RuntimeError:
                    h *= 0.5
            else:
                raise RuntimeError("native Jacobian pullback has no admissible probes")
        result.native_jac = result.jac @ np.column_stack(pullback)
    result.physical_coordinate_audit = {
        "jacobian_coordinates": "physical primitives; translated temperatures",
        "temperature_origins_k": origin.tolist(),
        "native_solver_bounds_preserved": True,
        "physical_residual_changed": False,
        "native_jacobian_available": bool(result.success),
    }
    return result


def centered_temperature_solve(solver, fun, x0, *, temperature_index, **kwargs):
    """Solve an affine, unit-slope temperature increment; return physical x."""
    initial = np.array(x0, dtype=float, copy=True)
    origin = float(initial[temperature_index])
    initial[temperature_index] = 0.0
    lower, upper = (np.array(v, dtype=float, copy=True) for v in kwargs["bounds"])
    lower[temperature_index] -= origin
    upper[temperature_index] -= origin
    kwargs["bounds"] = (lower, upper)

    def residual(coordinates):
        physical = np.array(coordinates, copy=True)
        physical[temperature_index] += origin
        return fun(physical)

    result = solver(residual, initial, **kwargs)
    result.x = np.array(result.x, copy=True)
    result.x[temperature_index] += origin
    # The unit-slope affine map leaves every Jacobian column and its units
    # unchanged; the finite-difference probes are now around increments.
    result.temperature_coordinate_audit = {
        "origin_k": origin,
        "temperature_index": temperature_index,
        "affine_slope": 1.0,
        "physical_residual_changed": False,
    }
    return result


def audit_film_profile(
    film,
    *,
    interface_hexane_mass_fraction,
    bulk_hexane_mass_fraction,
    fluxes,
    closure_residual_budget_kg_m2_s=0.0,
):
    """Keep native auditing, with a bounded pure-bulk endpoint fallback.

    Only the native reconstructed-profile range failure can enter this
    branch. The logarithm signs, monotonic direction, flux-coordinate closure
    and reconstructed endpoint must all pass independent error budgets.
    """
    wi, wb = interface_hexane_mass_fraction, bulk_hexane_mass_fraction
    if not math.isfinite(closure_residual_budget_kg_m2_s) or closure_residual_budget_kg_m2_s < 0:
        raise ValueError("finite nonnegative prospective film closure error budget required")
    try:
        efr.conditioned_ideal_binary_mass_coordinate_residual(
            film, interface_hexane_mass_fraction=wi, bulk_hexane_mass_fraction=wb, fluxes=fluxes
        )
        return None
    except efr.FilmProfileAdmissibilityError as exc:
        if (
            str(exc) != "candidate implies a mass-fraction profile outside [0,1]"
            or wb not in (0.0, 1.0)
            or not 0.0 < wi < 1.0
        ):
            raise
        original = exc
    total = fluxes.total_kg_m2_s
    conductance = film.film_density_kg_m3 * film.binary_mass_transfer_m_s
    pe = total / conductance
    numerator = fluxes.hexane_kg_m2_s - wb * total
    denominator = fluxes.hexane_kg_m2_s - wi * total
    if total == 0.0 or numerator * denominator <= 0.0:
        raise original
    phi = efr._inverse_exprel(pe)
    diffusive = conductance * (wi - wb) * phi
    residual = math.fsum((fluxes.hexane_kg_m2_s, -wi * total, -diffusive))
    flux_scale = math.fsum((abs(fluxes.hexane_kg_m2_s), abs(wi * total), abs(diffusive)))
    flux_budget = 128.0 * np.finfo(float).eps * flux_scale + closure_residual_budget_kg_m2_s
    profile_term = denominator / conductance * efr._exprel(pe)
    implied = wi - profile_term
    endpoint_budget = 128.0 * np.finfo(float).eps * math.fsum((1.0, abs(wi), abs(profile_term)))
    endpoint_budget += closure_residual_budget_kg_m2_s / conductance * abs(efr._exprel(pe))
    # d w_h/dz has sign -denominator throughout the exact exponential
    # profile. Require the direction of the prescribed physical endpoints.
    monotone = denominator * (wi - wb) > 0.0
    if (
        not monotone
        or not math.isfinite(implied)
        or abs(residual) > flux_budget
        or abs(implied - wb) > endpoint_budget
    ):
        raise original
    return {
        "branch": "pure_bulk_endpoint_numerical_closure",
        "bulk_hexane_mass_fraction": wb,
        "unmodified_implied_bulk_hexane_mass_fraction": implied,
        "conditioned_residual_kg_m2_s": residual,
        "flux_error_budget_kg_m2_s": flux_budget,
        "prospective_closure_residual_budget_kg_m2_s": closure_residual_budget_kg_m2_s,
        "closure_error_amplification_to_bulk_fraction": abs(efr._exprel(pe)) / conductance,
        "endpoint_error_budget": endpoint_budget,
        "epsilon_multiplier": 128,
        "profile_monotonicity_checked": True,
        "composition_clamped": False,
    }
