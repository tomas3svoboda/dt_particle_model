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

driver = boundary_driver.driver
core = driver.core
ct = core.ct
cp = core.cp
study = driver.study


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
):
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
    initial = np.r_[
        before.temperatures_k,
        fractions(before.temperatures_k, before.y_hexane),
        np.zeros(n),
        surface.temperature_k,
        surface.beta,
        (
            dt
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
    holder = {"x": np.zeros(3)}
    last = {}
    work = 0

    def evaluate(values):
        nonlocal work, dt, boundary
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
        surf_t, beta, third = values[-3:]
        if zero_film_event:
            dt = float(third)
            boundary = boundary_factory(dt)
        liquid = 0.0 if zero_film_event or gas_only or contact else third
        surface_y = None
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
        return np.r_[
            blocks.dry_water_mol_s,
            blocks.dry_hexane_mol_s,
            blocks.dry_datum_covariant_energy_w,
            fluxes[-1].surface_film_audit.water_node_residual_kg_m2_s,
            fluxes[-1].surface_film_audit.hexane_node_residual_kg_m2_s,
            fluxes[-1].surface_film_audit.energy_node_residual_w_m2,
        ]

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

        def residual(values):
            return evaluate(values + origins) / scales

        result = study.numerics.domain_aware_newton(
            residual,
            initial - origins,
            bounds=(np.asarray(bounds[0]) - origins, np.asarray(bounds[1]) - origins),
            max_nfev=1400,
            gtol=5e-11,
        )
        final = residual(result.x)
        if not result.success or max(abs(final)) > 1e-10:
            refusal = RuntimeError("fully dry joint root refused")
            refusal.optimizer_debug = {
                "message": str(result.message),
                "maximum_scaled_residual": float(max(abs(final))),
                "scaled_residual": final.tolist(),
                "physical_coordinates": (result.x + origins).tolist(),
                "residual_calls": int(result.nfev),
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
            "optimizer_stopping_threshold": 5e-11,
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
