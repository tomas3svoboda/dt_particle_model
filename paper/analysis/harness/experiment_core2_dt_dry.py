"""Stationary DT march using the same native actual-T/Y face law and water node.

Only numerical assembly is added here: native equilibrium, constitutive fluxes,
FV residuals, Stefan recovery and accepted-step ledgers remain Core2 calls.
The fully dry conversion is allowed only by the native uniform-label conversion.
"""

import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import pickle
import traceback
import numpy as np
import experiment_core2_dt_face as boundary_driver
from dtdc_simulator.core2.particle import cut_extinction_event as extinction
from dtdc_simulator.core2.particle.numerical_work import StepBudgetExceeded

driver = boundary_driver.driver
core = driver.core
ct = core.ct
cp = core.cp
study = driver.study


# Declared search constant of the eliminated depletion event: the largest
# duration shift (s) used to remove the water-closure remainder left by the
# representable surface temperature (scaled energy-row effect ~1e-11 or less).
EVENT_DURATION_SHIFT_S = 1e-8

# Owner ruling GT_PS2_OWNER_RULING_PAPER1_GATE_UNITS_2026-10-17, points 1 and 2.
GATE_UNITS_RULING = "GT_PS2_OWNER_RULING_PAPER1_GATE_UNITS_2026-10-17"
# Declared multiple m of the noise-aware energy-row gate max(tolerance, m * floor_i).
NOISE_AWARE_ENERGY_GATE_MULTIPLE = 2.0
# The original strict native component residual limit (checked on the
# _accept_step ledger below); the cell water and hexane search rows are
# weighted so that their search target is exactly this limit (or the
# tolerance, if tighter), as at runtime 131.
NATIVE_COMPONENT_RESIDUAL_LIMIT = 2e-11


class DrySearchTargetRefused(ValueError):
    """Typed refusal: a dry search target differing from the declared gate.

    Ruling point 1: the dry advance searches to its declared gate tolerance;
    no caller can set a tighter (or looser) search target.
    """


def declared_dry_search_target(tolerance, search_tolerance=None):
    """The dry advance's search target: its declared gate tolerance, nothing else."""
    if not math.isfinite(tolerance) or not tolerance > 0:
        raise ValueError("positive finite declared dry tolerance required")
    if search_tolerance is not None and search_tolerance != tolerance:
        raise DrySearchTargetRefused(
            f"dry search target {search_tolerance!r} differs from its declared gate "
            f"{tolerance!r}; ruled equal ({GATE_UNITS_RULING} point 1)"
        )
    return tolerance


def representability_floor(jacobian_rows, physical_coordinates):
    """sum_j |J_ij| ulp(x_j): the row's arithmetic floor at the iterate x."""
    rows = np.asarray(jacobian_rows, dtype=float)
    ulps = np.spacing(np.abs(np.asarray(physical_coordinates, dtype=float)))
    if rows.ndim != 2 or rows.shape[1] != ulps.size or not np.all(np.isfinite(rows)):
        raise ValueError("finite Jacobian rows over the physical unknowns required")
    return np.abs(rows) @ ulps


class NoiseAwareEnergyRowAcceptance:
    """Ruling point 2 on the dry search rows (the domain-aware Newton's acceptance).

    Every search row passes at the declared target, except the cell energy
    rows, which pass at max(target, m * floor_i), floor_i from the dense
    Jacobian at the iterate. The energy search rows are the gate-scaled rows
    (search scale = gate scale, weight 1), so one floor serves both clauses.
    """

    def __init__(self, target, energy_rows, origins, multiple=NOISE_AWARE_ENERGY_GATE_MULTIPLE):
        if not math.isfinite(multiple) or not multiple > 0:
            raise ValueError("positive finite declared gate multiple required")
        self.target = float(target)
        self.energy_rows = energy_rows
        self.origins = np.asarray(origins, dtype=float)
        self.multiple = float(multiple)

    def gates(self, x, jac):
        jac = np.asarray(jac, dtype=float)
        floor = representability_floor(jac[self.energy_rows], np.asarray(x) + self.origins)
        gates = np.full(jac.shape[0], self.target)
        gates[self.energy_rows] = np.maximum(self.target, self.multiple * floor)
        return gates, floor

    def candidate(self, residual):
        fixed = np.ones(np.asarray(residual).size, dtype=bool)
        fixed[self.energy_rows] = False
        return bool(np.all(np.abs(np.asarray(residual)[fixed]) <= self.target))

    def accepts(self, x, residual, jac):
        gates, _ = self.gates(x, jac)
        return bool(np.all(np.abs(np.asarray(residual)) <= gates))


class FilmDepletionBracketed(Exception):
    """Declared stop: a positive-film trial returned a nonpositive exterior film.

    Deliberately neither ValueError nor RuntimeError, so the bounded Newton
    search does not absorb it as a domain rejection; the dispatcher records
    it as the positive-film refusal and hands the remaining allowance to the
    exact-zero depletion event. Raised only when the caller opts in.
    """


def global_component_search_rows(original, cells, area):
    """Invertible replacement of surface mass rows by whole-particle balances."""
    values = np.asarray(original, dtype=float)
    if values.shape != (3 * cells + 3,) or not math.isfinite(area) or area <= 0:
        raise ValueError("native cell/surface row layout and positive area required")
    transformed = values.copy()
    transformed[-3] = math.fsum((*values[:cells], area * values[-3] / core.wa.M))
    transformed[-2] = math.fsum((*values[cells : 2 * cells], area * values[-2] / core.hx.M))
    return transformed


def advance(
    before,
    config,
    surface,
    carry,
    dt,
    boundary,
    *,
    zero_film_event=False,
    boundary_factory=None,
    surface_phase_override=None,
    surface_liquid_transfer_adapter=None,
    search_tolerance=None,
    physical_advance_assembly_budget=None,
    global_component_mass_error_budgets_kg=None,
    use_native_conduction_predictor=False,
    conservative_liquid_elimination=False,
    conservative_event_elimination=False,
    conservative_absent_elimination=False,
    event_duration_seed=None,
    stop_at_nonpositive_liquid=False,
    event_initial_profile=None,
):
    # The five trailing options are declared numerical search strategies only
    # (2026-10-17 film-depletion event-first build). With all of them at their
    # defaults every line below executes exactly as before.
    # Ruling point 1: the search target is the declared gate tolerance (the
    # former min(5e-11, tolerance) and its tighten-only override are removed).
    root_target = declared_dry_search_target(
        before.config.nonlinear_residual_tolerance, search_tolerance
    )
    if global_component_mass_error_budgets_kg is not None and (
        len(global_component_mass_error_budgets_kg) != 2
        or not all(math.isfinite(x) and x > 0 for x in global_component_mass_error_budgets_kg)
    ):
        raise ValueError("positive finite water/hexane per-advance error budgets required")
    grid = before.grid
    n = grid.n
    pressure = config.dry.pressure_pa
    pore = config.dry.pore
    old = ct._inventories(
        grid, ct._equilibrium_states(before.temperatures_k, before.y_hexane, before.config)
    )

    def fractions(ts, ys):
        return tuple(
            (y - cp.gas_only_composition_interval(t, pressure, pore).y_hexane_bounds[0])
            / (
                cp.gas_only_composition_interval(t, pressure, pore).y_hexane_bounds[1]
                - cp.gas_only_composition_interval(t, pressure, pore).y_hexane_bounds[0]
            )
            for t, y in zip(ts, ys)
        )

    gas_only = carry["external_water_loading_kg_m2"] == 0.0
    contact = surface_phase_override == "contact"
    if surface_phase_override is not None:
        if surface_phase_override not in ("liquid", "absent", "contact"):
            raise ValueError("unknown surface-phase override")
        if surface_phase_override == "absent" and not gas_only:
            raise ValueError("positive incoming film cannot be discarded by a phase override")
        if contact and (not gas_only or surface_liquid_transfer_adapter is None):
            raise ValueError(
                "zero-inventory wet contact needs zero incoming film and a transfer law"
            )
        gas_only = surface_phase_override == "absent"
    if zero_film_event and gas_only:
        raise ValueError("depletion event needs positive incoming liquid")
    if conservative_liquid_elimination and (gas_only or contact or zero_film_event):
        raise ValueError("fixed-time liquid elimination requires the positive-film active set")
    if conservative_event_elimination and (not zero_film_event or conservative_liquid_elimination):
        raise ValueError("eliminated depletion event requires the zero-film event active set")
    if conservative_absent_elimination and (not gas_only or contact or zero_film_event):
        raise ValueError("eliminated absent-film set requires the gas-only active set")
    if stop_at_nonpositive_liquid and not conservative_liquid_elimination:
        raise ValueError("the nonpositive-film stop belongs to the eliminated positive-film set")
    if event_duration_seed is not None and (
        not zero_film_event
        or not math.isfinite(event_duration_seed)
        or not 0 < event_duration_seed <= dt
    ):
        raise ValueError("a depletion duration seed must lie inside the attempted interval")
    reduced = (
        conservative_liquid_elimination
        or conservative_event_elimination
        or conservative_absent_elimination
    )
    attempted_duration = dt
    initial_temperatures = before.temperatures_k
    thermal_predictor_used = False
    if use_native_conduction_predictor:
        initial_temperatures, thermal_predictor_used = ct._conduction_predictor(
            before, dt, boundary.film_boundary, before.config
        )
    initial = np.r_[
        initial_temperatures,
        fractions(initial_temperatures, before.y_hexane),
        np.zeros(n),
        surface.temperature_k,
        surface.beta,
        (
            (dt if event_duration_seed is None else event_duration_seed)
            if zero_film_event
            else (
                0.01
                if contact
                else (
                    fractions([surface.temperature_k], [surface.surface_y_hexane])[0]
                    if gas_only
                    else surface.liquid_water_kg_m2
                )
            )
        ),
    ]
    event_initial_profile_used = None
    if event_initial_profile is not None:
        # Search seed only: a cell profile extrapolated by the caller to the
        # predicted event. An unusable profile keeps the former seed.
        if not zero_film_event:
            raise ValueError("an event initial profile belongs to the depletion event")
        event_initial_profile_used = False
        seed_t = np.asarray(event_initial_profile[0], dtype=float)
        seed_y = np.asarray(event_initial_profile[1], dtype=float)
        low, high = config.dry.conditioned_temperature_domain.solver_bounds_k
        if (
            seed_t.shape == (n,)
            and seed_y.shape == (n,)
            and np.all(np.isfinite(seed_t))
            and np.all(np.isfinite(seed_y))
            and np.all((low < seed_t) & (seed_t < high))
        ):
            try:
                seed_q = np.asarray(fractions(seed_t, seed_y), dtype=float)
            except ValueError:
                seed_q = None
            if seed_q is not None and np.all((0 < seed_q) & (seed_q < 1)):
                initial[:n] = seed_t
                initial[n : 2 * n] = seed_q
                event_initial_profile_used = True
    if reduced:
        initial = np.r_[initial[: 3 * n], surface.temperature_k]
    if conservative_event_elimination:
        initial = np.r_[
            initial[: 3 * n],
            attempted_duration if event_duration_seed is None else event_duration_seed,
        ]
    holder = {"x": np.zeros(3)}
    last = {}
    work = 0

    def evaluate(values):
        nonlocal work, dt, boundary
        if physical_advance_assembly_budget is not None:
            counter = physical_advance_assembly_budget
            if counter["used"] >= counter["limit"]:
                raise StepBudgetExceeded("dry physical-advance residual assembly budget exhausted")
            counter["used"] += 1
        ts = values[:n]
        q = values[n : 2 * n]
        stefan = values[2 * n : 3 * n]
        if not all(0 < v < 1 for v in q):
            raise ValueError("dry composition fraction outside strict phase chart")
        ys = []
        for t, v in zip(ts, q):
            lo, hi = cp.gas_only_composition_interval(t, pressure, pore).y_hexane_bounds
            ys.append(lo + (hi - lo) * v)
        cells = ct._equilibrium_states(ts, ys, before.config)
        new = ct._inventories(grid, cells)
        fluxes = [extinction._zero_center_flux()]
        for k in range(1, n):
            fluxes.append(
                core.cut._dry_flux(
                    ts[k - 1],
                    ys[k - 1],
                    ts[k],
                    ys[k],
                    grid.centers[k] - grid.centers[k - 1],
                    stefan[k - 1],
                    pore,
                    pore,
                    config,
                )
            )
        surface_y = None
        if conservative_liquid_elimination:
            from particle_conservative_surface import saturated_water_film_target

            surf_t = float(values[-1])
            target_h = (
                -math.fsum(b - a for a, b in zip(old.hexane_mol, new.hexane_mol))
                * core.hx.M
                / (grid.areas[-1] * dt)
            )
            total, components, beta, surface_y, inverse_calls = saturated_water_film_target(
                cells[-1], config, grid.R, grid.R - grid.centers[-1], boundary, surf_t, target_h
            )
            liquid = math.fsum(
                (
                    carry["external_water_loading_kg_m2"],
                    -math.fsum(b - a for a, b in zip(old.water_mol, new.water_mol))
                    * core.wa.M
                    / grid.areas[-1],
                    -dt * components.water_kg_m2_s,
                )
            )
            if stop_at_nonpositive_liquid and not liquid > 0:
                raise FilmDepletionBracketed(
                    "positive-film trial returned a nonpositive exterior film; "
                    "declared stop hands the allowance to the depletion event"
                )
            if not 0 < liquid < max(0.5, 4 * surface.liquid_water_kg_m2):
                raise ValueError(
                    "conservative liquid inventory left its original open bounds; solve a depletion event"
                )
            holder["x"] = np.array([surf_t, total, liquid])
            last["surface_component_inverse_evaluations"] = inverse_calls
            third = liquid
        elif conservative_event_elimination:
            from particle_conservative_surface import (
                depletion_event_surface_temperature,
                saturated_water_film_target,
            )

            # Duration is a search coordinate (its original declared bound);
            # the exact-zero water closure fixes the surface temperature,
            # the closure's well-conditioned variable, so W and H close by
            # construction (H through the native film inverse).
            dt = float(values[-1])
            if not 0 < dt <= 100.0:
                raise ValueError("event duration coordinate outside its declared bound")
            boundary = boundary_factory(dt)
            area = grid.areas[-1]
            outer = cells[-1]
            hexane_target = (
                -math.fsum(b - a for a, b in zip(old.hexane_mol, new.hexane_mol))
                * core.hx.M
                / (area * dt)
            )

            def film_target(trial_temperature):
                return saturated_water_film_target(
                    outer,
                    config,
                    grid.R,
                    grid.R - grid.centers[-1],
                    boundary,
                    trial_temperature,
                    hexane_target,
                )

            low, high = config.dry.conditioned_temperature_domain.solver_bounds_k
            surf_t, total, components, beta, surface_y, inverse_calls, closure, quantum = (
                depletion_event_surface_temperature(
                    carry["external_water_loading_kg_m2"],
                    math.fsum(b - a for a, b in zip(old.water_mol, new.water_mol))
                    * core.wa.M
                    / area,
                    dt,
                    film_target,
                    surface.temperature_k,
                    low,
                    high,
                )
            )
            # The saturated surface trace carries the native decode's
            # ulp-level scatter, so a representable T_s leaves a closure
            # remainder; at that fixed T_s the closure is smooth in the
            # duration, and one secant shift of at most EVENT_DURATION_SHIFT_S
            # (an energy-row effect far below the root target) removes it.
            # A larger or non-improving shift is not taken.
            shift = 0.0
            if closure != 0.0:
                water_gain = (
                    math.fsum(b - a for a, b in zip(old.water_mol, new.water_mol))
                    * core.wa.M
                    / area
                )
                hexane_loss = hexane_target * dt

                def closure_at(duration):
                    trial_boundary = boundary_factory(duration)
                    parts = saturated_water_film_target(
                        outer,
                        config,
                        grid.R,
                        grid.R - grid.centers[-1],
                        trial_boundary,
                        surf_t,
                        hexane_loss / duration,
                    )
                    value = math.fsum(
                        (
                            carry["external_water_loading_kg_m2"],
                            -water_gain,
                            -duration * parts[1].water_kg_m2_s,
                        )
                    )
                    return value, trial_boundary, parts

                probe = dt * 1e-7
                slope = (closure_at(dt + probe)[0] - closure) / probe
                if math.isfinite(slope) and slope != 0.0:
                    candidate = -closure / slope
                    if abs(candidate) <= EVENT_DURATION_SHIFT_S and 0 < dt + candidate <= 100.0:
                        value, trial_boundary, parts = closure_at(dt + candidate)
                        if abs(value) < abs(closure):
                            shift, closure = candidate, value
                            dt, boundary = dt + candidate, trial_boundary
                            total, components, beta, surface_y, _ = parts
            liquid = 0.0
            holder["x"] = np.array([surf_t, total, liquid])
            last["surface_component_inverse_evaluations"] = inverse_calls
            last["event_water_closure_kg_m2"] = closure
            last["event_water_closure_quantum_kg_m2"] = quantum
            last["event_duration_shift_s"] = shift
            third = liquid
        elif conservative_absent_elimination:
            from particle_conservative_surface import absent_film_surface_for_component_targets

            surf_t = float(values[-1])
            area = grid.areas[-1]
            total, components, beta, surface_y, corrections = (
                absent_film_surface_for_component_targets(
                    cells[-1],
                    config,
                    grid.R,
                    grid.R - grid.centers[-1],
                    boundary,
                    surf_t,
                    -math.fsum(b - a for a, b in zip(old.water_mol, new.water_mol))
                    * core.wa.M
                    / (area * dt),
                    -math.fsum(b - a for a, b in zip(old.hexane_mol, new.hexane_mol))
                    * core.hx.M
                    / (area * dt),
                )
            )
            liquid = 0.0
            holder["x"] = np.array([surf_t, total, liquid])
            last["surface_composition_corrections"] = corrections
            third = liquid
        else:
            surf_t, beta, third = values[-3:]
            if zero_film_event:
                dt = float(third)
                boundary = boundary_factory(dt)
            liquid = 0.0 if zero_film_event or gas_only or contact else third
            if gas_only:
                if not 0 < third < 1:
                    raise ValueError("gas-only surface fraction outside open phase chart")
                lo, hi = cp.gas_only_composition_interval(surf_t, pressure, pore).y_hexane_bounds
                surface_y = lo + (hi - lo) * third
            holder["x"] = np.array(
                [
                    surf_t,
                    driver.film_at_beta(
                        boundary, config, grid.R, surf_t, beta, surface_y_hexane=surface_y
                    )[0],
                    liquid,
                ]
            )
        fluxes.append(
            driver.node.evaluate_water_surface(
                holder,
                dt,
                carry["external_water_loading_kg_m2"],
                carry["external_water_temperature_k"],
                cells[-1],
                grid.R - grid.centers[-1],
                stefan[-1],
                pore,
                config,
                grid.R,
                boundary,
                surface_y_hexane=surface_y,
                surface_liquid_transfer_adapter=surface_liquid_transfer_adapter,
                surface_liquid_transfer_fraction=float(third) if contact else None,
            )
        )
        changes = extinction.ExtinctionInventoryChanges(
            tuple(b - a for a, b in zip(old.water_mol, new.water_mol)),
            tuple(b - a for a, b in zip(old.hexane_mol, new.hexane_mol)),
            tuple(b - a for a, b in zip(old.energy_j, new.energy_j)),
        )
        blocks = extinction._residual_blocks(grid, changes, fluxes, dt)
        last.update(cells=cells, new=new, fluxes=tuple(fluxes), blocks=blocks, changes=changes)
        work += 1
        original = np.r_[
            blocks.dry_water_mol_s,
            blocks.dry_hexane_mol_s,
            blocks.dry_datum_covariant_energy_w,
            fluxes[-1].surface_film_audit.water_node_residual_kg_m2_s,
            fluxes[-1].surface_film_audit.hexane_node_residual_kg_m2_s,
            fluxes[-1].surface_film_audit.energy_node_residual_w_m2,
        ]
        last["original_residual"] = original
        if reduced:
            return np.r_[original[: 3 * n], original[-1]]
        return (
            original
            if global_component_mass_error_budgets_kg is None
            else global_component_search_rows(original, n, grid.areas[-1])
        )

    tl, tu = config.dry.conditioned_temperature_domain.solver_bounds_k
    origins = np.r_[before.temperatures_k, np.zeros(2 * n), surface.temperature_k, 0.0, 0.0]
    bounds = (
        [tl] * n + [0.0] * n + [-math.inf] * n + [tl, -3.0, 0.0],
        [tu] * n
        + [1.0] * n
        + [math.inf] * n
        + [
            tu,
            3.0,
            (
                100.0
                if zero_film_event
                else 1.0 if gas_only or contact else max(0.5, 4 * surface.liquid_water_kg_m2)
            ),
        ],
    )
    if reduced:
        origins = np.r_[before.temperatures_k, np.zeros(2 * n), surface.temperature_k]
        bounds = (
            [tl] * n + [0.0] * n + [-math.inf] * n + [tl],
            [tu] * n + [1.0] * n + [math.inf] * n + [tu],
        )
    if conservative_event_elimination:
        origins = np.r_[before.temperatures_k, np.zeros(2 * n), 0.0]
        bounds = (
            [tl] * n + [0.0] * n + [-math.inf] * n + [0.0],
            [tu] * n + [1.0] * n + [math.inf] * n + [100.0],
        )
    with driver.liquid_internal_energy_node(
        holder,
        dt,
        carry["external_water_loading_kg_m2"],
        carry["external_water_temperature_k"],
        pressure,
    ):
        raw = evaluate(initial)
        capacities = [
            ct._equilibrium_energy_capacity(t, y, before.config) * v
            for t, y, v in zip(before.temperatures_k, before.y_hexane, grid.volumes)
        ]
        scales = np.r_[
            np.maximum(np.asarray(old.water_mol) / dt, np.abs(raw[:n])),
            np.maximum(np.asarray(old.hexane_mol) / dt, np.abs(raw[n : 2 * n])),
            np.maximum(np.asarray(capacities) / dt, np.abs(raw[2 * n : 3 * n])),
            max(carry["external_water_loading_kg_m2"] / dt, 0.05),
            max(
                old.total_hexane_mol * core.hx.M / (grid.areas[-1] * dt),
                abs(last["fluxes"][-1].surface_film_audit.external_hexane_flux_kg_m2_s),
                1e-12,
            ),
            max(abs(carry["external_water_energy_j"]) / (grid.areas[-1] * dt), 5000.0),
        ]
        search_scales = scales.copy()
        search_weights = np.ones(scales.size)
        # The cell water and hexane rows are searched at the strict native
        # component limit min(2e-11, tolerance) that _accept_step's ledger is
        # held to below: weight = target / that limit (5 at the ruled 1e-10,
        # 2.5 under the former 5e-11 target; the same limit either way).
        component_weight = root_target / min(root_target, NATIVE_COMPONENT_RESIDUAL_LIMIT)
        if reduced:
            search_scales = np.r_[scales[: 3 * n], scales[-1]]
            search_weights = np.ones(3 * n + 1)
            search_weights[: 2 * n] = component_weight
        elif global_component_mass_error_budgets_kg is not None:
            # Accepted component defects are bounded in kg per physical advance,
            # independently of the thermal cutoff. The caller reserves these
            # budgets from the original-feed ledger, before starting this solve.
            search_scales[-3] = global_component_mass_error_budgets_kg[0] / (
                core.wa.M * dt * root_target
            )
            search_scales[-2] = global_component_mass_error_budgets_kg[1] / (
                core.hx.M * dt * root_target
            )
            search_weights[: 2 * n] = component_weight

        def residual(values):
            return evaluate(values + origins) / search_scales * search_weights

        energy_rows = slice(2 * n, 3 * n)
        acceptance = NoiseAwareEnergyRowAcceptance(root_target, energy_rows, origins)
        result = study.numerics.domain_aware_newton(
            residual,
            initial - origins,
            bounds=(np.asarray(bounds[0]) - origins, np.asarray(bounds[1]) - origins),
            max_nfev=1400,
            gtol=root_target,
            acceptance=acceptance,
        )
        search_final = residual(result.x)
        final = last["original_residual"] / scales
        # Ruling point 2: the floor from the dense forward-difference Jacobian
        # the search returns at its iterate (its energy rows are the
        # gate-scaled rows); the declared 1e-10 gate clause keeps its value on
        # every other row.
        search_gates, energy_floor = acceptance.gates(result.x, result.jac)
        final_gates = np.full(final.size, 1e-10)
        final_gates[energy_rows] = np.maximum(1e-10, search_gates[energy_rows])
        energy_gate_stamp = {
            "ruling": GATE_UNITS_RULING,
            "declared_multiple": acceptance.multiple,
            "declared_lower_bound": root_target,
            "floor_definition": "sum_j |J_ij| numpy.spacing(x_j), J the domain-aware Newton's "
            "dense forward-difference Jacobian at the returned iterate, x the physical unknowns",
            "rows": [f"energy[{k}]" for k in range(n)],
            "residual": final[energy_rows].tolist(),
            "floor": energy_floor.tolist(),
            "gate": final_gates[energy_rows].tolist(),
            "absorbed_by_floor": (final_gates[energy_rows] > 1e-10).tolist(),
            "absorbed_rows": int(np.count_nonzero(final_gates[energy_rows] > 1e-10)),
            "maximum_residual_over_gate": float(
                np.max(np.abs(final[energy_rows]) / final_gates[energy_rows])
            ),
            "other_rows_gate": 1e-10,
            "search_target": root_target,
            "maximum_other_row_scaled_residual": float(
                np.max(np.abs(np.r_[final[: 2 * n], final[3 * n :]]))
            ),
        }
        if (
            not result.success
            or np.any(np.abs(search_final) > search_gates)
            or np.any(np.abs(final) > final_gates)
        ):
            refusal = RuntimeError("fully dry joint root refused")
            refusal.optimizer_debug = {
                "message": str(result.message),
                "maximum_scaled_residual": float(max(abs(final))),
                "scaled_residual": final.tolist(),
                "physical_coordinates": (result.x + origins).tolist(),
                "residual_calls": int(result.nfev),
                "noise_aware_energy_row_gate": energy_gate_stamp,
                "line_search": getattr(result, "line_search", None),
            }
            raise refusal
        fluxes = last["fluxes"]
        cells = last["cells"]
        new = last["new"]
        audit = fluxes[-1].surface_film_audit
        if abs(audit.beta) > boundary.accepted_beta_limit:
            raise RuntimeError("resolved dry surface beta exceeds declared band")
        recovered = cp.recover_total_stefan_flux(
            grid,
            old.water_mol,
            new.water_mol,
            old.hexane_mol,
            new.hexane_mol,
            tuple(f.retained_water_flux_mol_m2_s for f in fluxes),
            dt,
        )
        balance = cp.component_balance_residuals(
            grid,
            old.water_mol,
            new.water_mol,
            old.hexane_mol,
            new.hexane_mol,
            tuple(f.component for f in fluxes),
            dt,
        )
        # This native record is used for entropy/source auditing, not flux evaluation.
        surface_state = cp.evaluate_equilibrium(
            audit.temperature_k, pressure, audit.surface_y_hexane, pore
        )
        left = (cells[0], *cells)
        right = (cells[0], *cells[1:], surface_state)
        faces = ct._FaceConstitutiveData(
            tuple(f.independent_water_flux_mol_m2_s for f in fluxes),
            tuple(f.retained_water_flux_mol_m2_s for f in fluxes),
            tuple(f.conductive_heat_flux_w_m2 for f in fluxes),
            left,
            right,
            left,
            min(f.total_entropy_production_w_m3_k for f in fluxes),
            max(f.total_entropy_production_w_m3_k for f in fluxes),
            None,
        )
        evaluated = ct._ResidualEvaluation(
            cells,
            new,
            faces,
            recovered,
            tuple(f.component for f in fluxes),
            tuple(f.energy for f in fluxes),
            balance,
            last["blocks"].dry_energy_w,
            None,
        )
        # Raw energy rows must also close; the solver used the native invertible datum-covariant basis.
        energy_scales = np.asarray(scales[2 * n : 3 * n]) + abs(core.hx.H_REF) * np.asarray(
            scales[n : 2 * n]
        )
        # Ruling point 2 in the native ledger's raw-energy basis: its check
        # |raw_i| / s_i <= tolerance becomes |raw_i| / s_i <= gate_i by passing
        # s_i * gate_i / tolerance (factor 1 on every row the floor does not
        # absorb, so those scales are the former bytes). With the covariant
        # row within gate_i and the hexane row within 2e-11, the raw row
        # (covariant plus H_REF times hexane) is within gate_i by the mediant
        # inequality, so this adds no clause of its own. Stamped below.
        native_energy_gate_factors = search_gates[energy_rows] / root_target
        energy_scales = energy_scales * native_energy_gate_factors
        energy_gate_stamp["native_ledger_energy_scale_factors"] = (
            native_energy_gate_factors.tolist()
        )
        step = ct._accept_step(
            before,
            core.ct.DirichletPoreBoundary(
                audit.temperature_k,
                audit.surface_y_hexane,
                "resolved water surface; native actual-T/Y fluxes assembled by paper driver",
            ),
            dt,
            old,
            evaluated,
            scales[n : 2 * n],
            energy_scales,
            work,
            result.message,
            False,
            before.config,
        )
        if (
            global_component_mass_error_budgets_kg is not None or reduced
        ) and step.ledger.max_scaled_component_residual > min(
            2e-11, before.config.nonlinear_residual_tolerance
        ):
            raise RuntimeError("original strict native component residual limit refused")
        area = grid.areas[-1]
        film_energy = (
            area
            * audit.liquid_water_kg_m2
            * core.wa.state_Tp(audit.temperature_k, pressure, "liquid").u_mass
        )
        film_mass = area * audit.liquid_water_kg_m2
        combined_raw = [
            (new.total_water_mol - old.total_water_mol) * core.wa.M
            + film_mass
            - area * carry["external_water_loading_kg_m2"]
            + dt * area * audit.external_water_flux_kg_m2_s,
            (new.total_hexane_mol - old.total_hexane_mol) * core.hx.M
            + dt * area * audit.external_hexane_flux_kg_m2_s,
            new.total_energy_j
            - old.total_energy_j
            + film_energy
            - carry["external_water_energy_j"]
            + dt * area * audit.external_energy_flux_w_m2,
        ]
        combined_scales = [
            max(
                old.total_water_mol * core.wa.M + area * carry["external_water_loading_kg_m2"],
                1e-300,
            ),
            max(old.total_hexane_mol * core.hx.M, 1e-300),
            max(
                abs(old.total_energy_j) + abs(carry["external_water_energy_j"]),
                sum(capacities),
                1e-300,
            ),
        ]
        normalized = [abs(r) / s for r, s in zip(combined_raw, combined_scales)]
        if max(normalized) > 1e-10:
            raise RuntimeError("particle-plus-water-film dry ledger refused")
        report = {
            "time_s": step.after.time_s,
            "dt_s": dt,
            "accepted": True,
            "maximum_scaled_residual": float(max(abs(final))),
            "native_conduction_predictor_used": thermal_predictor_used,
            "conservative_liquid_mass_elimination": conservative_liquid_elimination,
            "surface_component_inverse_evaluations": last.get(
                "surface_component_inverse_evaluations"
            ),
            "maximum_search_scaled_residual": float(max(abs(search_final))),
            "dry_search_target": root_target,
            "noise_aware_energy_row_gate": energy_gate_stamp,
            # Lead decision (Class A), docs/GT_PS2_PAPER1_LINE_SEARCH_MERIT_BUILD_2026-10-18.md:
            # the merit the line search used and every accepted step's length.
            "line_search": getattr(result, "line_search", None),
            "newton_residual_calls": int(result.nfev),
            "newton_certificate_calls_outside_search_budget": int(
                getattr(result, "certificate_jacobian_calls_outside_search_budget", 0)
            ),
            "newton_message": str(result.message),
            "global_component_mass_error_budgets_kg": global_component_mass_error_budgets_kg,
            "surface_search_basis": (
                "algebraic native-film component elimination"
                if conservative_liquid_elimination
                else (
                    "whole-particle component balances"
                    if global_component_mass_error_budgets_kg is not None
                    else "original surface rows"
                )
            ),
            "ledger": asdict(step.ledger),
            "surface": asdict(audit),
            "combined_normalized_residuals": normalized,
            "mass_force_authority": config.moving_interface_composition_force_authority.value,
            "native_config_branch_is_coefficient_metadata_only": True,
            "external_water_phase": (
                "absent"
                if gas_only
                else "depletion event" if zero_film_event else "positive liquid"
            ),
            "temperature_k": list(step.after.temperatures_k),
            "y_hexane": list(step.after.y_hexane),
            "water_mol": new.total_water_mol,
            "hexane_mol": new.total_hexane_mol,
            "energy_j": new.total_energy_j,
        }
        if (
            conservative_event_elimination
            or conservative_absent_elimination
            or stop_at_nonpositive_liquid
            or event_duration_seed is not None
            or event_initial_profile is not None
        ):
            rows = last["original_residual"]
            report["declared_search_options"] = {
                "conservative_event_elimination": conservative_event_elimination,
                "conservative_absent_elimination": conservative_absent_elimination,
                "event_duration_seed_s": event_duration_seed,
                "stop_at_nonpositive_liquid": stop_at_nonpositive_liquid,
                "attempted_duration_s": attempted_duration,
                "event_initial_profile_used": event_initial_profile_used,
            }
            if conservative_event_elimination or conservative_absent_elimination:
                report["surface_search_basis"] = (
                    "algebraic native-film component elimination ("
                    + ("depletion event" if conservative_event_elimination else "absent film")
                    + ")"
                )
            report["global_component_defects_kg_per_advance"] = {
                "summed_original_rows_W_H": [
                    math.fsum((*rows[:n], area * rows[-3] / core.wa.M)) * core.wa.M * dt,
                    math.fsum((*rows[n : 2 * n], area * rows[-2] / core.hx.M)) * core.hx.M * dt,
                ],
                "particle_plus_film_ledger_W_H": combined_raw[:2],
            }
            if global_component_mass_error_budgets_kg is not None and (
                conservative_event_elimination or conservative_absent_elimination
            ):
                defects = report["global_component_defects_kg_per_advance"]
                if any(
                    abs(defects[key][i]) > global_component_mass_error_budgets_kg[i]
                    for key in defects
                    for i in (0, 1)
                ):
                    raise RuntimeError(
                        "eliminated component closure exceeds the declared per-advance budget"
                    )
            if "event_water_closure_kg_m2" in last:
                report["event_water_closure_kg_m2"] = last["event_water_closure_kg_m2"]
                report["event_water_closure_quantum_kg_m2"] = last[
                    "event_water_closure_quantum_kg_m2"
                ]
                report["event_duration_shift_s"] = last["event_duration_shift_s"]
                report["event_duration_shift_limit_s"] = EVENT_DURATION_SHIFT_S
            if "surface_composition_corrections" in last:
                report["surface_composition_corrections"] = last["surface_composition_corrections"]
        carry = {
            **carry,
            "external_water_loading_kg_m2": audit.liquid_water_kg_m2,
            "external_water_energy_j": film_energy,
            "external_water_temperature_k": audit.temperature_k,
        }
        return step.after, audit, carry, report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--zero-film-event", action="store_true")
    parser.add_argument("--until", type=float)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite evidence")
    checkpoint = pickle.loads(args.resume.read_bytes())
    current = checkpoint["current"]
    if hasattr(current, "to_uniform_fully_dry_transport_state"):
        current = current.to_uniform_fully_dry_transport_state()
    if not isinstance(current, ct.FullyDryTransportState):
        raise TypeError("native fully dry state required")
    old_root_tolerance = current.config.nonlinear_residual_tolerance
    # Owner RULE 2: align only the nonlinear numerical tolerance with the
    # preceding CutIntegratorState's 1e-10; all native conservation gates stay.
    current = replace(current, config=replace(current.config, nonlinear_residual_tolerance=1e-10))
    surface = checkpoint["surface"]
    carry = checkpoint["carry"]
    config = checkpoint["cut_config"]
    bath_path = core.DESTINY_PATH.with_name("destiny_v3_exit105.json")
    bath = json.loads(bath_path.read_text())
    sources = [
        Path(__file__),
        Path(boundary_driver.__file__),
        Path(driver.__file__),
        Path(driver.node.__file__),
        Path(study.__file__),
        Path(study.numerics.__file__),
        Path(driver.path_bounds.__file__),
        bath_path,
        *sorted((driver.REPO_ROOT / "src/dtdc_simulator/core2").rglob("*.py")),
    ]
    hashes = {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    record = {
        "complete_dt_destiny": False,
        "sources_at_entry": hashes,
        "attempts": [],
        "declared_numerical_tolerance": {
            "old_dry_root": old_root_tolerance,
            "new_dry_root": 1e-10,
            "optimizer_stopping_threshold": 1e-10,
            "search_target_ruling": GATE_UNITS_RULING,
            "ledger_tolerance": current.config.ledger_tolerance,
            "authority": "owner RULE 2; measured energy residual floor 1.22e-11, step energy error 4.46e-15; no physical coefficient changed",
        },
    }
    for i in range(args.steps):
        if args.until is not None and current.time_s >= args.until:
            break
        step_dt = args.dt if args.until is None else min(args.dt, args.until - current.time_s)
        try:
            boundary = boundary_driver.bath_boundary(bath, current.grid.R, current.time_s + step_dt)
            with (
                driver.bounded_path_proof(3000000, 120.0),
                driver.path_bounds.monotone_water_path_enclosure([]),
                study.magnitude_storage_scale(),
            ):
                after, audit, next_carry, report = advance(
                    current,
                    config,
                    surface,
                    carry,
                    step_dt,
                    boundary,
                    zero_film_event=args.zero_film_event,
                    boundary_factory=lambda duration: boundary_driver.bath_boundary(
                        bath, current.grid.R, current.time_s + duration
                    ),
                )
            current, surface, carry = after, audit, next_carry
            record["attempts"].append(report)
            print(
                json.dumps(
                    {
                        "step": i + 1,
                        "time_s": current.time_s,
                        "maximum_scaled_residual": report["maximum_scaled_residual"],
                    }
                ),
                flush=True,
            )
            if args.zero_film_event:
                break
        except Exception as exc:
            record["attempts"].append(
                {
                    "accepted": False,
                    "error": str(exc),
                    "traceback": traceback.format_exc(),
                    "optimizer_debug": getattr(exc, "optimizer_debug", None),
                }
            )
            break
    output = args.out.with_suffix(".accepted.pkl")
    output.write_bytes(
        pickle.dumps(
            {
                **checkpoint,
                "current": current,
                "surface": surface,
                "carry": carry,
                "source_sha256": {**checkpoint["source_sha256"], **hashes},
                "result_path": str(args.out),
            },
            protocol=5,
        )
    )
    record["accepted_checkpoint"] = str(output)
    record["sources_unchanged_during_run"] = all(
        hashlib.sha256(p.read_bytes()).hexdigest() == hashes[str(p.resolve())] for p in sources
    )
    record["bath_exit_reached"] = current.time_s >= bath["waypoints"][-1]["t_s"]
    record["stamp"] = driver.era_stamp(sources)
    args.out.write_text(json.dumps(record, indent=2))
    print(
        json.dumps(
            {
                "accepted_steps": sum(a["accepted"] for a in record["attempts"]),
                "last_time_s": current.time_s,
                "error": record["attempts"][-1].get("error"),
            }
        )
    )


if __name__ == "__main__":
    main()
