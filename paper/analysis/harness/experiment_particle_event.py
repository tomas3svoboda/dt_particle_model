"""Shared native particle event runner for the declared constant Faner bath.

Local trusted native checkpoints only. Exact face topology is committed by
the original Core2 event integrators, with the same conservative surface node.
"""

from dataclasses import asdict, fields, replace
from pathlib import Path
import argparse
import hashlib
import json
import math
import pickle
import traceback

import experiment_core2_faner as faner
import experiment_core2_dt_face as events
import experiment_core2_dt_projection as dry_event

FanerStudyControls = faner.FanerStudyControls  # trusted old CLI pickle identity
core, driver = faner.core, faner.shared


def positive_target(before, config, carry, surface, boundary, target, guess, budget):
    candidate = core.cut.candidate_from_state(
        before.transport,
        dry_total_stefan_fluxes_mol_m2_s=before.last_total_stefan_fluxes_mol_m2_s,
        interface_temperature_k=before.last_interface_temperature_k,
    )
    candidate = replace(candidate, front_z=target)
    seed = core.backend._seed(
        candidate, "native positive-target time solve; no imposed time history"
    )
    lo, hi = core.cp.gas_only_composition_interval(
        surface.temperature_k, config.dry.pressure_pa, config.dry.pore
    ).y_hexane_bounds
    step = driver.node.advance_water_node(
        before,
        guess,
        boundary,
        seed,
        old_loading=0.0,
        old_temperature=surface.temperature_k,
        surface_seed=[
            surface.temperature_k,
            surface.beta,
            (surface.surface_y_hexane - lo) / (hi - lo),
        ],
        physical_solver=True,
        residual_budget=budget,
        gas_only_surface=True,
        target_z=target,
        duration_bounds=(1e-12, 10.0),
        difference_mode="central",
        maximum_restart_trial_residual=0.05,
        beta_decoder=lambda t, b, y: driver.film_at_beta(
            boundary, config, before.transport.geometry.master_grid.R, t, b, surface_y_hexane=y
        )[0],
        energy_context=lambda holder: driver.liquid_internal_energy_node(
            holder, guess, 0.0, surface.temperature_k, config.dry.pressure_pa
        ),
    )
    return (
        step,
        step.surface_film_audit,
        {"positive_target_z": target, "duration_solved_from_native_balances": True},
    )


def face_summary(after):
    inv = after.current_inventory
    grid = after.geometry.master_grid
    dry = after.config.wet.wet.rho_dm_p * sum(grid.volumes)
    temperatures = (*after.wet_temperatures_k, *after.dry_temperatures_k)
    return {
        "time_s": after.time_s,
        "front_z": after.geometry.front.z,
        "X_hexane_kg_kg": inv.total_hexane_mol * core.hx.M / dry,
        "X_water_kg_kg": inv.total_water_mol * core.wa.M / dry,
        "T_k": sum(t * v for t, v in zip(temperatures, grid.volumes)) / sum(grid.volumes),
        "interface_temperature_k": after.interface_temperature_k,
        "hexane_mol": inv.total_hexane_mol,
        "water_mol": inv.total_water_mol,
        "energy_j": inv.total_energy_j,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--resume", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--event", choices=["arrival", "departure", "target", "projection"], required=True
    )
    p.add_argument("--target", type=float)
    p.add_argument("--dt", type=float, default=0.1)
    p.add_argument("--upper", type=float, default=10)
    p.add_argument("--budget", type=int, default=2400)
    p.add_argument("--native-budget", type=int)
    p.add_argument("--wall-budget", type=float, default=180)
    p.add_argument("--regional-tolerance", type=float, default=1e-10)
    p.add_argument("--tangent-tolerance", type=float, default=1e-10)
    p.add_argument("--condition-ceiling", type=float)
    p.add_argument("--event-time-tolerance", type=float, default=0.002)
    args = p.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite evidence")
    parent_path = args.resume.with_name(args.resume.name.replace(".accepted.pkl", ".json"))
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    if not parent["sources_unchanged_during_run"]:
        raise ValueError("parent executable changed during run")
    checkpoint = pickle.loads(args.resume.read_bytes())
    before = checkpoint["before"]
    if args.regional_tolerance > 1e-9:
        if not hasattr(before, "transport") or before.transport.geometry.front.z > 1e-7:
            raise ValueError("extended local resolution requires vanishing wet volume <=1e-7")
        kwargs = {
            f.name: getattr(before.controls, f.name)
            for f in fields(core.ci.CutSolverControls)
            if f.init
        }
        controls = faner.continued.VanishingCoreStudyControls(
            **{**kwargs, "nonlinear_residual_tolerance": args.regional_tolerance}
        )
        before = replace(before, controls=controls)
    if hasattr(before, "transport"):
        before = replace(
            before,
            controls=replace(before.controls, nonlinear_residual_tolerance=args.regional_tolerance),
        )
    else:
        before = replace(
            before,
            solver_controls=replace(
                before.solver_controls, nonlinear_residual_tolerance=args.regional_tolerance
            ),
        )
    if args.condition_ceiling is not None:
        if not hasattr(before, "transport"):
            raise ValueError("condition study applies only to radial state")
        before = replace(
            before,
            controls=replace(before.controls, maximum_condition_proxy=args.condition_ceiling),
        )
    config = checkpoint["config"]
    if args.native_budget is not None:
        if not 1 <= args.native_budget <= 24000:
            raise ValueError("declared native evaluation budget must be in [1,24000]")
        field = "controls" if hasattr(before, "transport") else "solver_controls"
        before = replace(
            before,
            **{
                field: replace(
                    getattr(before, field), maximum_function_evaluations=args.native_budget
                )
            },
        )
    carry = checkpoint["carry"]
    surface = checkpoint["surface"]
    source_paths = [
        Path(__file__),
        Path(faner.__file__),
        Path(events.__file__),
        Path(dry_event.__file__),
        Path(driver.__file__),
        Path(driver.node.__file__),
        Path(driver.study.__file__),
        Path(driver.study.numerics.__file__),
        *sorted((driver.REPO_ROOT / "src/dtdc_simulator/core2").rglob("*.py")),
    ]

    def hashes():
        return {str(x.resolve()): hashlib.sha256(x.read_bytes()).hexdigest() for x in source_paths}

    record = {
        "event": args.event,
        "accepted": False,
        "attempts": [],
        "accepted_parent": {
            "result": str(parent_path),
            "result_sha256": hashlib.sha256(parent_path.read_bytes()).hexdigest(),
            "checkpoint": str(args.resume),
            "checkpoint_sha256": hashlib.sha256(args.resume.read_bytes()).hexdigest(),
        },
        "source_sha256_at_entry": hashes(),
        "numerical_controls": asdict(
            before.controls if hasattr(before, "transport") else before.solver_controls
        ),
        "inputs": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "complete_falling_rate_experiment": False,
    }
    for path, expected in parent["source_sha256_at_exit"].items():
        if str(driver.REPO_ROOT / "src/dtdc_simulator/core2") in path:
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
                raise ValueError("native source changed from parent")
    boundary = checkpoint["boundary"]
    try:
        with (
            faner.declared_faner_film_multiplier(boundary.film_boundary.stress_multiplier),
            faner.pure_bulk_domain_adapter(),
            driver.bounded_path_proof(3000000, args.wall_budget),
            driver.path_bounds.monotone_hexane_path_enclosure([]),
            driver.study.magnitude_storage_scale(),
        ):
            if args.event == "arrival":
                step, surface, audit = events.advance_arrival(
                    before,
                    carry,
                    surface,
                    None,
                    time_guess=args.dt,
                    time_upper=args.upper,
                    budget=args.budget,
                    boundary_factory=lambda _: boundary,
                    gas_only_surface=True,
                    difference_mode="central",
                )
            elif args.event == "target":
                step, surface, audit = positive_target(
                    before, config, carry, surface, boundary, args.target, args.dt, args.budget
                )
            elif args.event == "projection":
                rows = checkpoint.get("target_refinements", [])
                if len(rows) < 4 or before.transport.geometry.front.z > 1e-10:
                    raise ValueError("four actual refinements and wet fraction <=1e-10 required")
                d = [r["dt_s"] for r in rows[-4:]]
                ratios = [b / a for a, b in zip(d, d[1:])]
                if not all(0 < q < 0.82 for q in ratios):
                    raise ValueError("refinement contraction not demonstrated")
                bound = d[-1] * 0.82 / (1 - 0.82)
                if bound > args.event_time_tolerance:
                    raise ValueError("conditional event time bound exceeds declared tolerance")
                duration = d[-1] * max(ratios) / (1 - max(ratios))
                candidate = dry_event.extinction.ExtinctionUnknowns(
                    before.transport.dry_temperatures_k,
                    before.transport.dry_y_hexane,
                    before.last_total_stefan_fluxes_mol_m2_s[1:],
                    duration,
                )
                seed = dry_event.projection.ExtinctionProjectionSeed(
                    candidate,
                    tuple(
                        max(abs(v), 0.01)
                        for v in candidate.positive_area_total_stefan_fluxes_mol_m2_s
                    ),
                    "native dry endpoint from measured actual positive-water time refinement",
                )
                with dry_event.water_surface_at_dry_endpoint(
                    before, carry, surface, duration, boundary
                ) as work:
                    step = dry_event.projection.solve_fixed_time_extinction_projection(
                        before, duration, boundary, seed
                    )
                surface = step.assembly.dry_face_fluxes[-1].surface_film_audit
                audit = {
                    "conditional_time_bound_s": bound,
                    "estimated_remaining_time_s": duration,
                    "last_four_durations_s": d,
                    "ratios": ratios,
                    "native_exact_limit_certificate_claimed": False,
                    "declared_event_time_tolerance_s": args.event_time_tolerance,
                    "surface_work": work,
                }
            else:
                step, surface, audit = events.advance_departure(
                    before,
                    carry,
                    surface,
                    None,
                    args.dt,
                    boundary_factory=lambda _: boundary,
                    tangent_residual_bound=args.tangent_tolerance,
                )
            combined = driver.combined_ledger(before, config, carry, step, surface)
            controls = before.controls if hasattr(before, "transport") else before.solver_controls
            if combined["maximum_normalized_residual"] > controls.ledger_tolerance:
                raise RuntimeError("combined event ledger refused")
            dry = (
                config.wet.wet.rho_dm_p
                * 4
                / 3
                * math.pi
                * (
                    before.transport.geometry.master_grid.R
                    if hasattr(before, "transport")
                    else before.geometry.master_grid.R
                )
                ** 3
            )
            if hasattr(step.after, "transport"):
                after = core.receding_summary(step.after, dry)
            elif hasattr(step.after, "geometry"):
                after = face_summary(step.after)
            else:
                inv = step.assembly.current_inventory
                after = {
                    "time_s": step.after.time_s,
                    "front_z": 0.0,
                    "X_hexane_kg_kg": inv.total_hexane_mol * core.hx.M / dry,
                    "X_water_kg_kg": inv.total_water_mol * core.wa.M / dry,
                    "T_k": sum(
                        t * v for t, v in zip(step.after.temperatures_k, step.after.grid.volumes)
                    )
                    / sum(step.after.grid.volumes),
                    "interface_temperature_k": None,
                    "hexane_mol": inv.total_hexane_mol,
                    "water_mol": inv.total_water_mol,
                    "energy_j": inv.total_energy_j,
                }
            record.update(
                accepted=True,
                last_accepted=after,
                ledger=asdict(step.ledger),
                surface=asdict(surface),
                audit=audit,
                combined_conservation=combined,
            )
            checkpoint.update(
                before=step.after,
                surface=surface,
                carry={**carry, **combined, "external_water_loading_kg_m2": 0.0},
                source_record=str(args.out),
            )
            if args.event == "target":
                checkpoint["target_refinements"] = [
                    *checkpoint.get("target_refinements", []),
                    {
                        "target_z": args.target,
                        "time_s": step.after.transport.time_s,
                        "dt_s": step.ledger.dt_s,
                    },
                ]
            if hasattr(step.after, "transport"):
                old_z = (
                    before.transport.geometry.front.z
                    if hasattr(before, "transport")
                    else before.geometry.front.z
                )
                checkpoint["front_rate"] = (
                    old_z - step.after.transport.geometry.front.z
                ) / step.ledger.dt_s
            args.out.with_suffix(".accepted.pkl").write_bytes(pickle.dumps(checkpoint))
    except Exception as error:
        record.update(
            error=str(error),
            debug=getattr(error, "optimizer_debug", None),
            native_residual_debug=getattr(error, "study_debug", None),
            traceback=traceback.format_exc(),
        )
    record["source_sha256_at_exit"] = hashes()
    record["sources_unchanged_during_run"] = (
        record["source_sha256_at_entry"] == record["source_sha256_at_exit"]
    )
    args.out.write_text(json.dumps(record, indent=2, allow_nan=True), encoding="utf-8")
    print(json.dumps({k: record.get(k) for k in ["accepted", "last_accepted", "error"]}, indent=2))


if __name__ == "__main__":
    main()
