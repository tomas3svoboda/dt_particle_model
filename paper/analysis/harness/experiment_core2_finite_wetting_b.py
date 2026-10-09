"""Finite liquid-contact sensitivity through native Regime B.

Additional retained-water path in the opened shell only; its wet-core
interface is exactly native. Rates are declared sensitivities, not fitted
or identified properties. Trusted local inputs only, one Python worker.
"""

import argparse
from contextlib import ExitStack, contextmanager
from dataclasses import asdict, fields, replace
import hashlib
import json
import math
from pathlib import Path
import pickle
import sys
import time
import traceback

import audit_core2_dt_external_ledger as independent
import experiment_core2_dt_birth as driver
import experiment_core2_dt_face as face
import experiment_core2_dt_projection as projection
import experiment_core2_radial_fast_wetting as wetting
import positive_water_continuation as continued

core, study = driver.core, driver.study


@contextmanager
def per_step_path_budget():
    """A fresh unchanged native proof budget for each physical advance."""
    active = []

    def reset():
        if active:
            active.pop().__exit__(None, None, None)
        context = driver.bounded_path_proof(3000000, 120.0)
        context.__enter__()
        active.append(context)

    try:
        yield reset
    finally:
        if active:
            active.pop().__exit__(None, None, None)


@contextmanager
def finite_shell_path(factor, conductance):
    """Thread the surface transfer through all existing native event charts."""
    original = driver.node.evaluate_water_surface
    with wetting.liquid_connected_fast_path(
        factor,
        exterior_liquid_present=True,
        contact_conductance=conductance,
        exclude_wet_front=True,
    ) as limiter:

        def evaluate(*args, **kwargs):
            if args[2] <= 0:
                raise ValueError(
                    "B prototype requires positive incoming film; no silent phase change"
                )
            return original(*args, **{**kwargs, "surface_liquid_transfer_adapter": limiter})

        driver.node.evaluate_water_surface = evaluate
        try:
            yield limiter
        finally:
            driver.node.evaluate_water_surface = original


def primitive_row(current, config, carry, surface, kind):
    transport = getattr(current, "transport", current)
    if hasattr(transport, "geometry") and hasattr(transport, "wet_temperatures_k"):
        primitives = asdict(transport)
        totals, temperature = independent.inventory(primitives, config)
        grid = transport.geometry.master_grid
        z = transport.geometry.front.z
        time_s = transport.time_s
    else:
        dry = current.to_uniform_fully_dry_transport_state()
        r = wetting.row(dry, surface, carry)
        return {**r, "kind": kind, "front_z": 0.0, "after_transport": asdict(dry)}
    water, hexane, energy = totals
    water *= core.wa.M
    hexane *= core.hx.M
    film = carry["external_water_mass_kg"]
    mass = config.dry.pore.rho_dm_p * math.fsum(grid.volumes)
    return {
        "kind": kind,
        "time_s": time_s,
        "front_z": z,
        "temperature_k": temperature,
        "particle_water_kg_kg_dry": water / mass,
        "external_water_kg_kg_dry": film / mass,
        "hexane_kg_kg_dry": hexane / mass,
        "moisture_wet_basis": (water + film) / (mass + water + film + hexane),
        "combined_W_H_U": [water + film, hexane, energy + carry["external_water_energy_j"]],
        "surface": asdict(surface),
        "after_transport": primitives,
    }


def advance_partial(current, carry, surface, bath, dt, target=None):
    config = current.transport.config
    radius = current.transport.geometry.master_grid.R
    candidate = core.cut.candidate_from_state(
        current.transport,
        dry_total_stefan_fluxes_mol_m2_s=current.last_total_stefan_fluxes_mol_m2_s,
        interface_temperature_k=current.last_interface_temperature_k,
    )
    candidate = replace(candidate, front_z=target or candidate.front_z * 0.999)
    seed = core.backend._seed(candidate, "finite liquid contact; native moving-front candidate")

    def factory(duration):
        return face.bath_boundary(bath, radius, current.transport.time_s + duration)

    boundary = factory(dt)

    def energy(holder):
        return driver.liquid_internal_energy_node(
            holder,
            lambda: holder["dt"],
            carry["external_water_loading_kg_m2"],
            carry["external_water_temperature_k"],
            config.dry.pressure_pa,
        )

    def decode(ts, beta, actual=None):
        return driver.film_at_beta(actual or boundary, config, radius, ts, beta)[0]

    return driver.node.advance_water_node(
        current,
        dt,
        boundary,
        seed,
        old_loading=carry["external_water_loading_kg_m2"],
        old_temperature=carry["external_water_temperature_k"],
        surface_seed=[surface.temperature_k, surface.beta, surface.liquid_water_kg_m2],
        physical_solver=True,
        residual_budget=2400,
        beta_decoder=decode,
        energy_context=energy,
        target_z=target,
        boundary_factory=factory if target else None,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--factor", type=float, default=1000)
    parser.add_argument("--contact-conductance", type=float, required=True)
    parser.add_argument("--dt", type=float, default=0.2)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--until", type=float, default=510.0)
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    if args.dt <= 0 or args.steps <= 0:
        raise ValueError("positive duration and step count required")
    if args.out.exists() or args.out.with_suffix(".accepted.pkl").exists():
        raise ValueError("refusing to overwrite evidence")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    root = driver.REPO_ROOT.resolve()
    feed_path = (
        root / "paper/analysis/results_2026-10-06/shared_model/dt_feed_massfilm_gated_R1p5.json"
    )
    bath_path = core.DESTINY_PATH.with_name("destiny_v3_exit105.json")
    feed, bath = json.loads(feed_path.read_text()), json.loads(bath_path.read_text())
    before, config, controls, carry = driver.initialize_from_feed(
        feed, coefficient_case="db_max_dapp_max", continued_water=True
    )
    sources = {feed_path.resolve(), bath_path.resolve(), Path(__file__).resolve()}
    sources.update(root.joinpath("src/dtdc_simulator/core2").rglob("*.py"))
    for module in list(sys.modules.values()):
        name = getattr(module, "__file__", None)
        if name and Path(name).suffix == ".py" and Path(name).resolve().is_relative_to(root):
            sources.add(Path(name).resolve())
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)}
    record = {
        "experiment": "finite liquid contact from native zero-old-dry-volume birth through Regime B",
        "physically_qualifying": False,
        "conductance_identified": False,
        "contact_conductance_mol_m2_s": args.contact_conductance,
        "shell_mobility_factor": args.factor,
        "law": "N_add=-ln(a_R/a_L)/(d/L_add+1/K_contact) at exterior; N_add=-L_add*grad(ln a) within shell",
        "wet_front_extra_flux": "exactly zero; native wet-core/interface transport retained",
        "native_and_combined_ledger_tolerance": 1e-10,
        "birth_dt_s": 0.01,
        "dt_requested_s": args.dt,
        "source_sha256_at_entry": hashes,
        "history": [],
        "attempts": [],
        "status": "STARTED",
    }
    current = surface = None
    initial_streams = [[], [], []]
    for t, w, v, h, oil in zip(
        before.temperatures_k,
        before.retained_water_loadings,
        before.grid.volumes,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
    ):
        cell = driver.ww.evaluate_cell(t, w, h, oil, config.wet)
        for stream, density in zip(
            initial_streams,
            (
                cell.retained_water_concentration_mol_m3 * core.wa.M,
                cell.total_hexane_concentration_kg_m3,
                cell.energy_density_j_m3,
            ),
        ):
            stream.append(v * density)
    initial_inventory = [math.fsum(v) for v in initial_streams]
    initial_inventory[0] += carry["external_water_mass_kg"]
    initial_inventory[2] += carry["external_water_energy_j"]
    area = before.grid.areas[-1]
    integrated = [[], [], []]
    scales = [
        abs(initial_inventory[0]),
        abs(initial_inventory[1]),
        carry["import_capacity_scale_j"],
    ]
    record["initial_time_s"] = before.time_s
    record["initial_combined_W_H_U"] = initial_inventory
    record["independent_audit_scales"] = scales
    started = time.perf_counter()
    try:
        with (
            finite_shell_path(args.factor, args.contact_conductance) as limiter,
            per_step_path_budget() as reset_path_budget,
            driver.path_bounds.monotone_water_path_enclosure([]),
            driver.path_bounds.monotone_hexane_path_enclosure([]),
            ExitStack() as numerical_scope,
        ):
            if args.resume:
                checkpoint = pickle.loads(args.resume.read_bytes())
                if checkpoint.get("finite_parameters") != [args.factor, args.contact_conductance]:
                    raise ValueError("finite conductance resume parameters changed")
                for name, expected in checkpoint["source_sha256"].items():
                    if (
                        Path(name).is_relative_to(root / "src/dtdc_simulator/core2")
                        or Path(name).resolve() in (feed_path.resolve(), bath_path.resolve())
                    ) and hashes.get(name) != expected:
                        raise ValueError("parent native source changed")
                current, surface, carry = (checkpoint[k] for k in ("current", "surface", "carry"))
                config = getattr(current, "transport", current).config
                record["parent_sha256"] = hashlib.sha256(args.resume.read_bytes()).hexdigest()
                record["parent_checkpoint"] = str(args.resume)
                record["history"].append(primitive_row(current, config, carry, surface, "resumed"))
                initial_inventory = record["history"][-1]["combined_W_H_U"]
                record["initial_time_s"] = getattr(current, "transport", current).time_s
                record["initial_combined_W_H_U"] = initial_inventory
            scale_installed = False
            for index in range(args.steps):
                reset_path_budget()
                if current is not None and not scale_installed:
                    numerical_scope.enter_context(study.magnitude_storage_scale())
                    scale_installed = True
                previous = current
                if current is None:
                    boundary = face.bath_boundary(bath, before.grid.R, before.time_s + 0.01)
                    step, birth_record = driver.advance_birth(
                        before,
                        config,
                        controls,
                        carry,
                        boundary,
                        0.01,
                        2400,
                        solver_policy="physical",
                    )
                    record["birth"] = birth_record
                    kind = "birth"
                else:
                    if not hasattr(current, "transport") and not getattr(
                        current, "requires_face_departure_chart", False
                    ):
                        record["status"] = "REGIME_B_COMPLETED"
                        break
                    transport = getattr(current, "transport", current)
                    if transport.time_s >= args.until:
                        record["status"] = "BOUNDED_WINDOW_COMPLETED"
                        break
                    z = transport.geometry.front.z
                    if getattr(current, "requires_face_departure_chart", False):
                        step, _, work = face.advance_departure(
                            current, carry, surface, bath, 0.05, tangent_residual_bound=2e-10
                        )
                        record["departure_seed_policy"] = {
                            "tangent_residual_bound": 2e-10,
                            "measured_default_refusal": 1.086e-10,
                            "tangent_used_as_accepted_time_step": False,
                            "accepted_component_energy_tolerance": 1e-10,
                            "basis": "owner RULE 2; seed-only search tolerance, native final solve unchanged",
                        }
                        kind = "face_departure"
                    elif z <= 1e-12:
                        targets = [
                            r
                            for r in record["attempts"]
                            if r.get("kind") == "positive_target" and r.get("accepted")
                        ]
                        if len(targets) < 4:
                            raise ValueError(
                                "four sequential target durations required before extinction projection"
                            )
                        durations = [r["ledger"]["dt_s"] for r in targets[-4:]]
                        ratios = [b / a for a, b in zip(durations, durations[1:])]
                        if (
                            not all(0 < q < 0.82 for q in ratios)
                            or durations[-1] * 0.82 / 0.18 > 0.002
                        ):
                            raise ValueError("prospective 2 ms extinction-time contraction refused")
                        dt = durations[-1] * max(ratios) / (1 - max(ratios))
                        boundary = face.bath_boundary(
                            bath, transport.geometry.master_grid.R, transport.time_s + dt
                        )
                        candidate = projection.extinction.ExtinctionUnknowns(
                            transport.dry_temperatures_k,
                            transport.dry_y_hexane,
                            current.last_total_stefan_fluxes_mol_m2_s[1:],
                            dt,
                        )
                        seed = projection.projection.ExtinctionProjectionSeed(
                            candidate,
                            tuple(
                                max(abs(v), 0.01)
                                for v in candidate.positive_area_total_stefan_fluxes_mol_m2_s
                            ),
                            "finite contact; measured native diminishing-core time sequence",
                        )
                        with projection.water_surface_at_dry_endpoint(
                            current, carry, surface, dt, boundary
                        ):
                            step = projection.projection.solve_fixed_time_extinction_projection(
                                current, dt, boundary, seed
                            )
                        record["projection_time_bound_s"] = durations[-1] * 0.82 / 0.18
                        kind = "extinction_projection"
                    elif transport.geometry.cut_cell_index == 1 and z < 0.18:
                        step, _, work = face.advance_arrival(
                            current, carry, surface, bath, time_guess=0.8
                        )
                        kind = "face_arrival"
                    else:
                        target = z * 0.5 if z < 0.03 else None
                        if z <= 1e-7:
                            declarations = {
                                f.name: getattr(current.controls, f.name)
                                for f in fields(core.ci.CutSolverControls)
                                if f.init
                            }
                            current = replace(
                                current,
                                controls=continued.VanishingCoreStudyControls(
                                    **{
                                        **declarations,
                                        "nonlinear_residual_tolerance": 2e-8,
                                        "maximum_condition_proxy": 1e18,
                                    }
                                ),
                            )
                            record["vanishing_core_numerical_policy"] = (
                                "inherited 2026-10-06 local 2e-8, conditioning 1e18 only below z=1e-7; W/H/U remain 1e-10"
                            )
                        step = advance_partial(
                            current,
                            carry,
                            surface,
                            bath,
                            args.dt if target is None else 0.1,
                            target,
                        )
                        kind = "positive_target" if target else "partial"
                new_surface = getattr(step, "surface_film_audit", None)
                if new_surface is None:
                    new_surface = step.assembly.dry_face_fluxes[-1].surface_film_audit
                combined = driver.combined_ledger(
                    before if previous is None else previous, config, carry, step, new_surface
                )
                if combined["maximum_normalized_residual"] > 1e-10:
                    raise RuntimeError("combined component/energy ledger refused")
                next_carry = {
                    **carry,
                    **combined,
                    "external_water_loading_kg_m2": new_surface.liquid_water_kg_m2,
                }
                next_row = primitive_row(step.after, config, next_carry, new_surface, kind)
                elapsed = next_row["time_s"] - (
                    before.time_s
                    if previous is None
                    else getattr(previous, "transport", previous).time_s
                )
                next_integrated = [s.copy() for s in integrated]
                for stream, rate in zip(
                    next_integrated,
                    (
                        new_surface.external_water_flux_kg_m2_s,
                        new_surface.external_hexane_flux_kg_m2_s,
                        new_surface.external_energy_flux_w_m2,
                    ),
                ):
                    stream.append(elapsed * area * rate)
                residuals = [
                    math.fsum((a, -b, math.fsum(stream)))
                    for a, b, stream in zip(
                        next_row["combined_W_H_U"], initial_inventory, next_integrated
                    )
                ]
                normalized = [abs(r) / s for r, s in zip(residuals, scales)]
                next_row["independent_cumulative_residuals"] = residuals
                next_row["independent_normalized_residuals"] = normalized
                if max(normalized) > 1e-10:
                    raise RuntimeError("independent archived-state cumulative W/H/U refused")
                current, surface, carry = step.after, new_surface, next_carry
                integrated = next_integrated
                record["history"].append(next_row)
                record["attempts"].append(
                    {
                        "index": index,
                        "accepted": True,
                        "kind": kind,
                        "ledger": asdict(step.ledger),
                        "combined": combined,
                        "liquid_transfer": limiter.last_audit,
                    }
                )
                checkpoint = {
                    "current": current,
                    "surface": surface,
                    "carry": carry,
                    "cut_config": config,
                    "source_sha256": hashes,
                    "finite_parameters": [args.factor, args.contact_conductance],
                    "result_path": str(args.out),
                }
                args.out.with_suffix(".progress.pkl").write_bytes(
                    pickle.dumps(checkpoint, protocol=5)
                )
                args.out.with_suffix(".progress.json").write_text(
                    json.dumps(
                        {
                            "status": "RUNNING",
                            "accepted_steps": len(record["history"]),
                            "last": {
                                k: v
                                for k, v in record["history"][-1].items()
                                if k not in ("surface", "after_transport")
                            },
                        },
                        indent=2,
                    )
                )
                print(
                    json.dumps(
                        {
                            "index": index,
                            "kind": kind,
                            "time_s": record["history"][-1]["time_s"],
                            "z": record["history"][-1]["front_z"],
                            "film": carry["external_water_loading_kg_m2"],
                        }
                    ),
                    flush=True,
                )
            else:
                record["status"] = "STEP_BUDGET_REACHED"
    except Exception as exc:
        current = previous if "previous" in locals() else current
        record.update(
            status="REFUSED",
            error=str(exc),
            traceback=traceback.format_exc(),
            optimizer_debug=getattr(exc, "optimizer_debug", None),
        )
    record["sources_unchanged_during_run"] = all(
        hashlib.sha256(p.read_bytes()).hexdigest() == hashes[str(p)] for p in sources
    )
    record["wall_s"] = time.perf_counter() - started
    if record["history"]:
        record["last"] = record["history"][-1]
        if "checkpoint" in locals():
            args.out.with_suffix(".accepted.pkl").write_bytes(pickle.dumps(checkpoint, protocol=5))
    args.out.write_text(json.dumps(record, indent=2))
    print(json.dumps({k: record.get(k) for k in ("status", "error", "wall_s")}), flush=True)


if __name__ == "__main__":
    main()
