"""Numerical search bounds for a resolved, evaporatively cooled surface.

The old outer-cell temperature is an initial state, not a lower thermodynamic
bound on the new surface. Physical pore, phase and flux checks remain native.
"""

from contextlib import contextmanager
from dataclasses import replace
import math

import numpy as np
import finite_film_core2_study as study


@contextmanager
def joint_surface_temperature_bounds(
    before, seed, *, solver_policy="scipy", solver_budget=400, surface_coordinate_limit=60.0
):
    if not 30.0 <= surface_coordinate_limit <= 70.0:
        raise ValueError("declared surface barrier search limit must lie in [30,70]")
    original = study.least_squares
    physical_lower, physical_upper = (
        before.transport.config.dry.conditioned_temperature_domain.solver_bounds_k
    )
    new_lower = float(np.nextafter(physical_lower, math.inf))
    old_cell_temperature = seed.candidate.dry_temperatures_k[-1]
    audit = {
        "reason": "active-water evaporation permits surface cooling below the old outer-cell temperature",
        "physical_dry_temperature_domain_k": [physical_lower, physical_upper],
        "requested_surface_temperature_lower_k": new_lower,
        "joint_calls_modified": 0,
        "other_solver_calls_unchanged": 0,
        "solver_bounds": [],
        "physical_admission_checks_changed": False,
        "residual_tolerances_changed": False,
        "assembly_budget_changed": False,
    }

    def widened_joint_solver(fun, x0, *args, **kwargs):
        if getattr(fun, "__qualname__", "") != "advance_joint_study.<locals>.residual":
            audit["other_solver_calls_unchanged"] += 1
            return original(fun, x0, *args, **kwargs)
        if args or "bounds" not in kwargs:
            raise RuntimeError(
                "joint study solver invocation changed; review the surface-bound adapter"
            )
        old_lower, old_upper = map(lambda value: np.asarray(value, dtype=float), kwargs["bounds"])
        coordinates = np.asarray(x0, dtype=float)
        if (
            coordinates.ndim != 1
            or coordinates.size <= 2
            or old_lower.shape != coordinates.shape
            or old_upper.shape != coordinates.shape
        ):
            raise RuntimeError(
                "joint study coordinate layout changed; surface-bound amendment refused"
            )
        if old_lower[-2] != old_cell_temperature:
            raise RuntimeError("joint study no longer has the old outer-cell lower bound")
        lower = old_lower.copy()
        lower[-2] = new_lower
        upper = old_upper.copy()
        # This barrier sums several phase/authority logarithms and has scale
        # two. A bound of 30 is NOT a saturation endpoint (measured a_h gap
        # 1.9e-5). Widen only its numerical surface search, retaining all
        # strict native gas-only activity inequalities and final ledgers.
        lower[-1] = -surface_coordinate_limit
        upper[-1] = surface_coordinate_limit
        audit["surface_composition_coordinate_limit"] = surface_coordinate_limit
        audit["original_surface_composition_coordinate_limit"] = float(old_upper[-1])
        if not lower[-2] < old_upper[-2]:
            raise RuntimeError("physical surface temperature interval is empty")
        if np.any(coordinates < lower) or np.any(coordinates > old_upper):
            raise RuntimeError(
                "existing joint seed is outside the declared physical numerical bounds"
            )
        audit["joint_calls_modified"] += 1
        audit["solver_bounds"].append(
            {
                "old_lower": old_lower.tolist(),
                "new_lower": lower.tolist(),
                "old_upper": old_upper.tolist(),
                "new_upper": upper.tolist(),
                "initial_coordinates": coordinates.tolist(),
                "initial_coordinates_valid": True,
                "changed_coordinate_index": int(coordinates.size - 2),
                "changed_coordinate": "resolved_surface_temperature_k",
            }
        )
        kwargs["bounds"] = (lower, upper)
        kwargs["jac"] = "2-point"
        kwargs["diff_step"] = None
        audit["jacobian_policy"] = "native-style automatic two-point increments"
        audit["solver_policy"] = solver_policy
        solver = original
        if solver_policy == "newton":
            solver = study.numerics.domain_aware_newton
            kwargs["max_nfev"] = solver_budget - 3
        if solver_policy == "physical":
            core = study.core
            chart = core.ci._coordinate_chart(before)
            layout = before.transport.layout
            scales = seed.stefan_flux_scales_mol_m2_s
            pore = replace(
                before.transport.config.dry.pore, w_o=before.transport.oil_fraction_labels[-1]
            )
            pressure = before.transport.config.dry.pressure_pa
            nw, nd, nf = layout.wet_piece_count, layout.dry_piece_count, layout.dry_face_count

            def fraction(t, y, p):
                interval = core.cp.gas_only_composition_interval(t, pressure, p)
                lo, hi = interval.y_hexane_bounds
                return (y - lo) / (hi - lo)

            def composition(t, q, p):
                if not 0 < q < 1:
                    raise ValueError("gas interval fraction outside its open domain")
                interval = core.cp.gas_only_composition_interval(t, pressure, p)
                lo, hi = interval.y_hexane_bounds
                return lo + (hi - lo) * q

            def decode(values):
                candidate = core.ci._decode_candidate_for_chart(values[:-2], layout, chart, scales)
                y = core.cp.decode_gas_only_y(values[-2], pressure, values[-1], pore)
                fractions = tuple(
                    fraction(t, y, p)
                    for t, y, p in zip(
                        candidate.dry_temperatures_k, candidate.dry_y_hexane, chart.dry_pores
                    )
                )
                candidate = replace(candidate, dry_y_hexane=fractions)
                return np.r_[candidate.vector(), values[-2], fraction(values[-2], y, pore)]

            def encode(values):
                candidate = core.cut.CutTransportUnknowns.from_vector(layout, values[:-2])
                ys = tuple(
                    composition(t, q, p)
                    for t, q, p in zip(
                        candidate.dry_temperatures_k, candidate.dry_y_hexane, chart.dry_pores
                    )
                )
                candidate = replace(candidate, dry_y_hexane=ys)
                native = core.ci._encode_candidate_for_chart(candidate, chart, scales)
                y = core.cp.encode_gas_only_y(
                    values[-2], pressure, composition(values[-2], values[-1], pore), pore
                )
                return np.r_[native, values[-2], y]

            pairs = (
                [chart.wet_temperature] * nw
                + [chart.wet_water] * nw
                + [chart.dry_temperature] * nd
                + [(0.0, 1.0)] * nd
                + [(-math.inf, math.inf)] * nf
                + [chart.front_z, chart.interface_temperature]
                + [(new_lower, old_upper[-2]), (0.0, 1.0)]
            )
            result = study.numerics.physical_chart_solve(
                fun,
                x0,
                decode=decode,
                encode=encode,
                physical_bounds=tuple(zip(*pairs)),
                temperature_indices=[
                    *range(nw),
                    *range(2 * nw, 2 * nw + nd),
                    len(pairs) - 3,
                    len(pairs) - 2,
                ],
                **{
                    **kwargs,
                    "max_nfev": solver_budget - 3,
                    "gtol": before.controls.nonlinear_residual_tolerance,
                },
            )
            audit["physical_coordinate_audit"] = result.physical_coordinate_audit
            audit["composition_coordinate"] = (
                "linear fraction of native temperature-dependent gas interval"
            )
            audit["optimizer_residual_stopping_threshold"] = (
                before.controls.nonlinear_residual_tolerance
            )
            return result
        result = study.numerics.centered_temperature_solve(
            solver, fun, x0, temperature_index=-2, **kwargs
        )
        audit["temperature_coordinate"] = result.temperature_coordinate_audit
        return result

    study.least_squares = widened_joint_solver
    try:
        yield audit
    finally:
        study.least_squares = original
        audit["original_solver_callable_restored"] = study.least_squares is original
