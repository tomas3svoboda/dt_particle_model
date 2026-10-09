"""Bounded initialized Faner radial entry, with active water and exact pure bath.

This is a fresh, declared initial-value experiment, not feed-to-birth evidence.
No historical profiles, bridges, donor classes or endpoint tails are imported.
An unaccepted attempt changes no committed trajectory state.
"""

from dataclasses import asdict, replace
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[4]
HARNESS = ROOT / "paper" / "analysis" / "harness"
sys.path.insert(0, str(HARNESS))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scipy.optimize import brentq, least_squares  # noqa: E402
import finite_film_core2_study as study  # noqa: E402
from bulk_endpoint import (  # noqa: E402
    declared_faner_film_multiplier,
    make_pure_bulk_boundary,
    pure_bulk_domain_adapter,
)
from assembly_budget import assembly_budget  # noqa: E402
from surface_bounds import joint_surface_temperature_bounds  # noqa: E402

core = study.core


def phase_point(water, pressure, pore):
    """Native virial/liquid fugacity equations: a_h=1, a_w=Luikov(W)."""
    target_aw = core.sp.water_activity(water, pore.luikov)

    def equations(x):
        temperature, y = map(float, x)
        gas = core.bg.state(temperature, pressure, 1.0 - y, y, k_wh=pore.k_wh)
        water_liquid = core.wa.state_Tp(temperature, pressure, "liquid")
        hexane_liquid = core.hx.state_Tp(temperature, pressure, "liquid")
        return [
            math.log(gas.fugacity_hexane_pa / hexane_liquid.fugacity),
            math.log(gas.fugacity_water_pa / (target_aw * water_liquid.fugacity)),
        ]

    root = least_squares(
        equations,
        [339.0, 0.96],
        bounds=([330.0, 0.6], [342.0, 1.0 - 1e-12]),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=100,
    )
    if not root.success or max(map(abs, equations(root.x))) > 1e-10:
        raise RuntimeError("active-water native phase initializer did not close")
    return tuple(map(float, root.x))


def active_binary_controls(original, config):
    """One numerical temperature chart for every declared moisture case.

    The upper endpoint is the native binary phase point at the existing
    minimum admitted retained-water loading. All moisture cases use this
    same chart. Phase equations and acceptance checks remain unchanged.
    """
    pressure = config.dry.pressure_pa

    minimum_water = original.wet_water_bounds[0]
    endpoint_temperature, endpoint_y = phase_point(minimum_water, pressure, config.dry.pore)
    upper = math.nextafter(endpoint_temperature, -math.inf)
    controls = replace(
        original,
        interface_temperature_bounds_k=(original.interface_temperature_bounds_k[0], upper),
        wet_temperature_bounds_k=(330.0, 345.0),
    )
    audit = {
        "reason": "shared active-water phase roots at low moisture exceed the old 340 K numerical chart",
        "old_interface_temperature_bounds_k": list(original.interface_temperature_bounds_k),
        "new_interface_temperature_bounds_k": list(controls.interface_temperature_bounds_k),
        "old_wet_temperature_bounds_k": list(original.wet_temperature_bounds_k),
        "new_wet_temperature_bounds_k": list(controls.wet_temperature_bounds_k),
        "pressure_pa": pressure,
        "physical_upper_endpoint": "native active-binary phase root at the unchanged minimum retained-water loading",
        "minimum_retained_water_loading": minimum_water,
        "endpoint_temperature_k": endpoint_temperature,
        "endpoint_y_hexane": endpoint_y,
        "old_controls": asdict(original),
        "new_controls": asdict(controls),
        "binary_phase_checks_changed": False,
        "physical_wet_caloric_bounds_changed": False,
        "acceptance_tolerances_changed": False,
        "work_budgets_changed": False,
        "wet_water_bounds_changed": False,
    }
    old_y_bounds = config.interface_composition.y_hexane_bounds
    new_y_bounds = (old_y_bounds[0], max(old_y_bounds[1], (1.0 + endpoint_y) / 2.0))
    config = replace(
        config,
        interface_composition=replace(config.interface_composition, y_hexane_bounds=new_y_bounds),
    )
    audit["old_interface_composition_root_bracket"] = list(old_y_bounds)
    audit["new_interface_composition_root_bracket"] = list(new_y_bounds)
    return controls, config, audit


def run(args):
    source_paths = [
        Path(__file__),
        Path(__file__).with_name("bulk_endpoint.py"),
        Path(__file__).with_name("assembly_budget.py"),
        Path(__file__).with_name("surface_bounds.py"),
        Path(study.__file__),
        Path(study.numerics.__file__),
        Path(core.__file__),
        Path(core.backend.__file__),
        *(
            [
                ROOT / "paper/analysis/harness/experiment_core2_dt_birth.py",
                ROOT / "paper/analysis/harness/core2_water_film_node.py",
                ROOT / "paper/analysis/harness/native_water_path_bounds.py",
            ]
            if args.conservative_node
            else []
        ),
        *sorted((ROOT / "src/dtdc_simulator/core2").rglob("*.py")),
    ]
    source_entry_hashes = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths
    }
    record = {
        "experiment": "initialized active-binary Faner particle entry",
        "complete_falling_rate_experiment": False,
        "feed_to_birth_proven": False,
        "historical_chart_used": False,
        "bulk_temperature_k": args.bath_temperature,
        "bulk_y_hexane": 1.0,
        "initial_water_loading": args.water,
        "water_is_dynamic": True,
        "same_conservative_dt_surface_node": args.conservative_node,
        "native_assembly_budget_per_step": args.budget,
        "initial_step_dt_s": args.dt,
        "continuation_dt_s": args.continuation_dt,
        "source_sha256_at_entry": source_entry_hashes,
        "attempts": [],
    }
    started = time.perf_counter()
    try:
        installed = core.q.install_engineering_case(
            "nominal_log_center",
            stress_multiplier=1.0,
            structural_conductivity_w_m_k=0.24,
        )
        config = installed.configuration.cut
        wet = replace(config.wet.wet, X_water=args.water)
        config = replace(
            config,
            wet=replace(config.wet, wet=wet),
            moving_interface_composition_force_authority=core.cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH,
        )
        record["shared_moving_interface_force_authority"] = (
            config.moving_interface_composition_force_authority.value
        )
        controls, config, controls_audit = active_binary_controls(core.backend.CONTROLS, config)
        record["active_binary_temperature_chart"] = controls_audit
        gamma_t, gamma_y = phase_point(args.water, core.PRESSURE_PA, config.dry.pore)
        wet_t = gamma_t + args.initial_wet_temperature_rise
        dry_t = gamma_t + args.initial_dry_temperature_rise
        target_water_fugacity = (
            core.sp.water_activity(args.water, config.dry.pore.luikov)
            * core.wa.state_Tp(dry_t, core.PRESSURE_PA, "liquid").fugacity
        )
        dry_y = brentq(
            lambda y: core.bg.state(
                dry_t, core.PRESSURE_PA, 1 - y, y, k_wh=config.dry.pore.k_wh
            ).fugacity_water_pa
            - target_water_fugacity,
            1e-12,
            1 - 1e-12,
            xtol=1e-14,
        )
        core.cp.evaluate_equilibrium(dry_t, core.PRESSURE_PA, dry_y, config.dry.pore)
        grid = core.backend.uniform_grid(args.cells, args.radius)
        geometry = core.backend.cg.partition_master_grid(
            grid, radius_m=args.front_fraction * args.radius
        )
        layout = core.cut.layout_for_geometry(geometry)
        dry_mass = wet.rho_dm_p * 4 / 3 * math.pi * args.radius**3

        # Native saturated-void capacity and sorption determine the wet
        # historical label. Total initial loading is an output, never a target
        # used to force extra solvent into a declared thin-shell wet core.
        history = core.wet_core.activate(grid, wet_t, wet).loadings
        radial = core.cut.initialize_cut_state(
            geometry,
            config,
            [wet_t] * layout.wet_piece_count,
            [args.water] * layout.wet_piece_count,
            [dry_t] * layout.dry_piece_count,
            [dry_y] * layout.dry_piece_count,
            history,
            oil_fraction_labels=[wet.w_o] * args.cells,
        )
        before = core.ci.initialize_integrator_state(
            radial, controls, interface_temperature_k=gamma_t
        )
        record["initialization"] = {
            "front_fraction": args.front_fraction,
            "wet_temperature_k": wet_t,
            "wet_interface_phase_temperature_k": gamma_t,
            "declared_initial_wet_temperature_rise_k": args.initial_wet_temperature_rise,
            "dry_temperature_k": dry_t,
            "phase_y_hexane": gamma_y,
            "dry_y_hexane": dry_y,
            "native_activation_historical_hexane_loadings": list(history),
            "initial_inventory": asdict(core.ci.inventory_snapshot(radial)),
            "summary": core.receding_summary(before, dry_mass),
        }
        pair = core.film_pair_at_radius(args.radius)
        boundary = make_pure_bulk_boundary(
            pair,
            temperature_k=args.bath_temperature,
            forcing_timescale_s=args.forcing_timescale,
            beta_limit=args.beta_limit,
            stress_multiplier=args.film_multiplier,
        )
        record["film_pair"] = asdict(pair)
        record["film_common_multiplier"] = args.film_multiplier
        step = None
        for index in range(args.steps):
            if time.perf_counter() - started > args.wall_budget:
                record["stop"] = "bounded wall budget reached between completed steps"
                break
            step_dt = (
                args.dt if index == 0 or args.continuation_dt is None else args.continuation_dt
            )
            candidate = core.cut.candidate_from_state(
                before.transport,
                dry_total_stefan_fluxes_mol_m2_s=(
                    core.stefan_seed_for_faces(before.transport.layout.dry_face_count)
                    if index == 0
                    else before.last_total_stefan_fluxes_mol_m2_s
                ),
                interface_temperature_k=before.last_interface_temperature_k,
            )
            if index == 0:
                front_guess = candidate.front_z * 0.999
            else:
                front_guess = (
                    candidate.front_z
                    - (
                        step.before.transport.geometry.front.z
                        - step.after.transport.geometry.front.z
                    )
                    * step_dt
                    / step.ledger.dt_s
                )
            candidate = replace(candidate, front_z=front_guess)
            seed = core.backend._seed(candidate, "active binary Faner declared initialized state")
            identity_before = before
            budget_audit = None
            bulk_domain_audit = None
            multiplier_audit = None
            surface_bounds_audit = None
            journal = []
            try:
                with assembly_budget(
                    boundary,
                    limit=args.budget,
                    remaining_wall_seconds=args.wall_budget - (time.perf_counter() - started),
                ) as budget_audit:
                    with declared_faner_film_multiplier(args.film_multiplier) as multiplier_audit:
                        with pure_bulk_domain_adapter() as bulk_domain_audit:
                            with study.magnitude_storage_scale() as journal:
                                with joint_surface_temperature_bounds(
                                    before,
                                    seed,
                                    solver_policy=args.solver,
                                    solver_budget=args.budget,
                                    surface_coordinate_limit=args.surface_coordinate_limit,
                                ) as surface_bounds_audit:
                                    if args.conservative_node:
                                        import experiment_core2_dt_birth as dt_driver

                                        surf_t = (
                                            candidate.dry_temperatures_k[-1]
                                            if index == 0
                                            else step.surface_film_audit.temperature_k
                                        )
                                        surf_y = (
                                            candidate.dry_y_hexane[-1]
                                            if index == 0
                                            else step.surface_film_audit.surface_y_hexane
                                        )
                                        lo, hi = core.cp.gas_only_composition_interval(
                                            surf_t, core.PRESSURE_PA, config.dry.pore
                                        ).y_hexane_bounds
                                        with dt_driver.path_bounds.monotone_hexane_path_enclosure(
                                            journal
                                        ):
                                            step = dt_driver.node.advance_water_node(
                                                before,
                                                step_dt,
                                                boundary,
                                                seed,
                                                old_loading=0.0,
                                                old_temperature=surf_t,
                                                surface_seed=[
                                                    surf_t,
                                                    (
                                                        0.0
                                                        if index == 0
                                                        else step.surface_film_audit.beta
                                                    ),
                                                    (surf_y - lo) / (hi - lo),
                                                ],
                                                physical_solver=True,
                                                residual_budget=args.budget,
                                                gas_only_surface=True,
                                                maximum_restart_trial_residual=0.05,
                                                beta_decoder=lambda t, b, y: dt_driver.film_at_beta(
                                                    boundary,
                                                    config,
                                                    args.radius,
                                                    t,
                                                    b,
                                                    surface_y_hexane=y,
                                                )[0],
                                                energy_context=lambda holder: dt_driver.liquid_internal_energy_node(
                                                    holder, step_dt, 0.0, surf_t, core.PRESSURE_PA
                                                ),
                                            )
                                    else:
                                        step = study.advance_joint_study(
                                            before,
                                            step_dt,
                                            boundary,
                                            seed,
                                            surface_seed=(
                                                None
                                                if index == 0
                                                else (
                                                    step.surface_film_audit.surface_temperature_k,
                                                    step.surface_film_audit.surface_y_hexane,
                                                )
                                            ),
                                        )
                record["attempts"].append(
                    {
                        "index": index,
                        "accepted": True,
                        "after": core.receding_summary(step.after, dry_mass),
                        "ledger": asdict(step.ledger),
                        "film": asdict(step.surface_film_audit),
                        "numerical_storage_scale_audit": journal,
                        "assembly_budget": budget_audit,
                        "bulk_domain_adapter": bulk_domain_audit,
                        "film_multiplier_admission": multiplier_audit,
                        "surface_temperature_search_bounds": surface_bounds_audit,
                    }
                )
                before = step.after
            except Exception as exc:
                record["attempts"].append(
                    {
                        "index": index,
                        "accepted": False,
                        "error": str(exc),
                        "error_type": type(exc).__name__,
                        "classification": getattr(exc, "classification", "NATIVE_SOLVER_REFUSAL"),
                        "assembly_budget": budget_audit,
                        "bulk_domain_adapter": bulk_domain_audit,
                        "film_multiplier_admission": multiplier_audit,
                        "surface_temperature_search_bounds": surface_bounds_audit,
                        "numerical_storage_scale_audit": journal,
                        "before_state_identity_preserved": before is identity_before,
                        "debug": getattr(exc, "study_debug", None),
                        "optimizer_debug": getattr(exc, "optimizer_debug", None),
                        "last_accepted": core.receding_summary(before, dry_mass),
                        "traceback": traceback.format_exc(),
                    }
                )
                break
    except Exception as exc:
        record["setup_refusal"] = str(exc)
        record["traceback"] = traceback.format_exc()
    record["wall_seconds"] = time.perf_counter() - started
    source_exit_hashes = {
        str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths
    }
    record["source_sha256"] = source_exit_hashes
    record["source_sha256_at_exit"] = source_exit_hashes
    record["sources_unchanged_during_run"] = source_entry_hashes == source_exit_hashes
    if not record["sources_unchanged_during_run"]:
        record["provenance_refusal"] = "SOURCE_CHANGED_DURING_RUN"
        record["result_is_reportable"] = False
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--water", type=float, default=0.05)
    parser.add_argument("--radius", type=float, default=0.0008078927451503978)
    parser.add_argument("--cells", type=int, default=2)
    parser.add_argument("--front-fraction", type=float, default=0.999)
    parser.add_argument("--initial-dry-temperature-rise", type=float, default=0.3)
    parser.add_argument("--initial-wet-temperature-rise", type=float, default=0.0)
    parser.add_argument("--bath-temperature", type=float, default=393.15)
    parser.add_argument("--film-multiplier", type=float, default=0.14578543636088787)
    parser.add_argument("--forcing-timescale", type=float, default=29.74922457229942)
    parser.add_argument("--beta-limit", type=float, default=1.0)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--continuation-dt", type=float)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--wall-budget", type=float, default=120)
    parser.add_argument("--budget", type=int, default=400)
    parser.add_argument("--surface-coordinate-limit", type=float, default=60.0)
    parser.add_argument("--conservative-node", action="store_true")
    parser.add_argument("--solver", choices=("scipy", "newton", "physical"), default="scipy")
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite an existing result")
    result = run(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "out": str(args.out),
                "setup_refusal": result.get("setup_refusal"),
                "attempts": [
                    {"accepted": a["accepted"], "error": a.get("error")} for a in result["attempts"]
                ],
            }
        )
    )


if __name__ == "__main__":
    main()
