"""Fresh original-feed particle history under one immutable runtime source set.

Reuse the existing feed, birth, ordinary-cut, exact-arrival, conservative mesh,
core-projection and dry/surface operators. This generates candidate evidence;
an individual completed case does not establish mesh/time convergence.
"""

# Capture the command's directory before legacy analysis imports can change it.
# ruff: noqa: E402

from pathlib import Path

INVOCATION_CWD = Path.cwd()

import argparse
from contextlib import nullcontext
from dataclasses import asdict, replace
import hashlib
import json
import math
import pickle
import time
import traceback
from types import SimpleNamespace

import experiment_core2_dt_feed as feed_driver
import experiment_core2_faner as faner
import faner_measured_trace as faner_trace
import experiment_core2_finite_wetting_b as engine
import experiment_core2_radial_fast_wetting as liquid
import particle_kernel_envelope as envelope
import particle_surface_graded_mesh as graded_mesh
import qualify_sliding_front_history as front
from particle_core_extinction import project_core_endpoint
from particle_history_time_control import HistoryTimeControl, depletion_advance_allowance_factor
from particle_postcore_solver import advance_dry_phases, DryPhaseAdvanceError
from particle_solver_configuration import ParticleModelParameters, PAPER_CHECKPOINT_CODEC
from dtdc_simulator.core2.particle.numerical_work import StepBudgetExceeded

ROOT = faner.REPO_ROOT
# The dry stage's declared maximum physical advances: the runner's ceiling, and the count that divides
# the per-advance original-feed component budgets. 1,500 is the declared default (the ladder of
# GT_PS2_PAPER1_FINAL_LADDER_AND_SENSITIVITY_2026-10-27); --dry-advance-ceiling N states another value
# for a declared run (owner ruling of 2026-10-28, "Ok, raise the ceiling, run half step tests": 4,000
# for the DT half-step pair). The record carries the value in force.
DECLARED_DRY_ADVANCE_CEILING = 1500


def declared_dry_advance_ceiling(args):
    """The dry stage's maximum physical advances: the declared default unless the run states
    --dry-advance-ceiling N, a positive integer; an absent attribute is the default."""
    value = getattr(args, "dry_advance_ceiling", None)
    if value is None:
        return DECLARED_DRY_ADVANCE_CEILING
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError("--dry-advance-ceiling must be a positive integer")
    return value


def hashes():
    return {
        p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for directory in (ROOT / "src", ROOT / "paper/analysis/harness", ROOT / "scripts")
        for p in sorted(directory.rglob("*.py"))
    }


def predicted_film_depletion_duration(points):
    """Linear extrapolation of the last two accepted positive exterior films.

    Returns the duration from the last point to the extrapolated zero, or None
    when there is no decreasing positive pair. Search seeding only.
    """
    if len(points) < 2:
        return None
    (t0, l0), (t1, l1) = points[-2], points[-1]
    if not (t1 > t0 and l0 > 0 and l1 > 0 and l1 < l0):
        return None
    duration = l1 * (t1 - t0) / (l0 - l1)
    return duration if math.isfinite(duration) and duration > 0 else None


def curvature_corrected_depletion_duration(points):
    """Quadratic extrapolation of the last three accepted positive exterior films.

    Declared search option (seeding only): the film's decay decelerates as the
    uptake falls, so the two-point line under-predicts the depletion duration;
    under the ruled imbibition bridge the miss carried the event-first solve
    past its attempted interval (evidence 121). Returns the duration from the
    last point to the first zero of the parabola through the last three
    points, the two-point line when only two exist, and None when the films
    do not decrease or the parabola has no zero after the last point.
    """
    if len(points) < 3:
        return predicted_film_depletion_duration(points)
    (t0, l0), (t1, l1), (t2, l2) = points[-3:]
    if not (t0 < t1 < t2 and l0 > l1 > l2 > 0):
        return predicted_film_depletion_duration(points)
    d01 = (l1 - l0) / (t1 - t0)
    d12 = (l2 - l1) / (t2 - t1)
    curvature = (d12 - d01) / (t2 - t0)
    # Newton form about t2: l(t2 + s) = l2 + b s + curvature s^2, b = slope at t2
    slope = d12 + curvature * (t2 - t1)
    if curvature == 0.0:
        return predicted_film_depletion_duration(points)
    discriminant = slope * slope - 4.0 * curvature * l2
    if not discriminant >= 0.0:
        return None
    root = math.sqrt(discriminant)
    candidates = [(-slope - root) / (2.0 * curvature), (-slope + root) / (2.0 * curvature)]
    positive = [c for c in candidates if math.isfinite(c) and c > 0]
    return min(positive) if positive else None


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, default=front.json_default), encoding="utf-8")


def fresh_dt_birth(args, multiplier, pinned):
    feed_path, birth_path = args.out / "feed.json", args.out / "birth.json"
    controls = SimpleNamespace(
        bath=args.bath,
        out=feed_path,
        radius_m=0.0015,
        cells=args.cells,
        dt_s=0.1 / multiplier,
        steam_dt_s=0.02 / multiplier,
        max_time_s=520.0,
        paper_properties=True,
    )
    feed = feed_driver.run(controls)
    write_json(feed_path, feed)
    pressure_schedule = envelope.armed(envelope.PRESSURE_SCHEDULE)
    before, config, native_controls, carry = engine.driver.initialize_from_feed(
        feed,
        cells=args.cells,
        coefficient_case="db_max_dapp_max",
        continued_water=True,
        property_parameters=ParticleModelParameters(),
        dry_local_residual_tolerance=args.dry_local_tolerance,
        **({"pressure_pa": feed["kernel_envelope_feed_pressure_pa"]} if pressure_schedule else {}),
    )
    bath = json.loads(args.bath.read_text())
    dt = 0.01 / multiplier
    boundary = engine.face.bath_boundary(bath, before.grid.R, before.time_s + dt)
    birth_target = {}
    if pressure_schedule:
        birth_pressure = envelope.bath_pressure_pa(bath, before.time_s + dt)
        if birth_pressure != config.dry.pressure_pa:
            birth_target = {
                "target_config": engine.core.cut.pressure_only_target_config(config, birth_pressure)
            }
            # The newborn interface seed and the interface chart's lower end are the
            # declared native two-liquid temperature (+1e-6 K for the chart) at the
            # birth's own pressure, as the feed stage takes them at its own.
            root = envelope.two_liquid_root_k(birth_pressure)
            carry = {
                **carry,
                "feed_two_liquid_interface_temperature_k": carry[
                    "incoming_two_liquid_interface_temperature_k"
                ],
                "incoming_two_liquid_interface_temperature_k": root,
            }
            native_controls = replace(
                native_controls,
                interface_temperature_bounds_k=(
                    root + 1e-6,
                    native_controls.interface_temperature_bounds_k[1],
                ),
            )
            carry["birth_interface_solver_chart_k"] = native_controls.interface_temperature_bounds_k
    engine.wetting.surface_contact_reset()
    with (
        # The adopted film-covered surface law at the declared contact conductance (the one
        # constant particle_kernel_envelope.K_L_MOL_M2_S).
        engine.finite_shell_path(100, envelope.surface_contact_conductance()),
        engine.driver.bounded_path_proof(3000000, 120) as proof,
        engine.driver.path_bounds.monotone_water_path_enclosure([]),
        engine.driver.path_bounds.monotone_hexane_path_enclosure([]),
    ):
        step, audit = engine.driver.advance_birth(
            before,
            config,
            native_controls,
            carry,
            boundary,
            dt,
            2400,
            solver_policy="physical",
            trace_coordinates=True,
            **birth_target,
            # Ruling packet 6 item 2 (GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md)
            **({"budget_per_unknown": True} if envelope.budgets_per_unknown() is not None else {}),
        )
    surface = step.assembly.dry_face_fluxes[-1].surface_film_audit
    birth_contact = dict(engine.wetting.FINITE_CONTACT_LAST)
    checkpoint = {
        "checkpoint_codec": PAPER_CHECKPOINT_CODEC,
        "current": step.after,
        # Repair A2, GT_PS2_UNIVERSAL_KERNEL_BIRTH_PRESSURE_REPAIR_BUILD_2026-10-22
        "config": step.after.transport.config,
        "carry": {
            **carry,
            **audit["whole_particle_and_external_water"],
            "external_water_loading_kg_m2": surface.liquid_water_kg_m2,
        },
        "surface": surface,
        "front_rate": (1.0 - step.after.transport.geometry.front.z) / dt,
        "feed_path": str(feed_path),
        "bath_path": str(args.bath),
    }
    birth_path.with_suffix(".accepted.pkl").write_bytes(pickle.dumps(checkpoint, protocol=5))
    write_json(
        birth_path,
        {
            "checkpoint_codec": PAPER_CHECKPOINT_CODEC,
            "accepted": True,
            "source_sha256_at_exit": {str(ROOT / p): digest for p, digest in pinned.items()},
            "sources_unchanged_during_run": hashes() == pinned,
            "feed_sha256": hashlib.sha256(feed_path.read_bytes()).hexdigest(),
            "audit": audit,
            "proof": proof,
            "dt_s": dt,
            # The adopted film-covered surface law (report only): its declared row and the last
            # surface-face evaluation of the birth's scope, with the conductance it used.
            "surface_contact_law": {
                "declared_row": envelope.surface_contact_law_declared_row(),
                "K_l_mol_m2_s": birth_contact.get("K_l_mol_m2_s"),
                "last_surface_face": birth_contact or None,
            },
        },
    )
    return birth_path


def fresh_faner_birth(args, multiplier):
    path = args.out / "birth.json"
    # Primary Paper 1 source-mean radius and declared initial moisture. Other
    # physical and solver defaults come from the existing study parser.
    controls = faner.argument_parser().parse_args(
        list(
            map(
                str,
                [
                    "--out",
                    path,
                    "--water",
                    0.025,
                    "--radius",
                    0.000885,
                    "--cells",
                    args.cells,
                    "--feed-dt",
                    0.25 / multiplier,
                    "--birth-dt",
                    0.01 / multiplier,
                    "--coefficient-case",
                    "db_max_dapp_max",
                    "--continued-water",
                    "--paper-properties",
                    "--auto-leading-seed",
                    "--dry-local-tolerance",
                    args.dry_local_tolerance,
                    "--steps",
                    0,
                    *(
                        ["--grouped-birth-search"]
                        if getattr(args, "grouped_birth_search", False)
                        else []
                    ),
                ],
            )
        )
    )
    with engine.driver.bounded_path_proof(3000000, 120) as proof:
        record = faner.run(controls)
    record["original_feed_and_birth_proof_work"] = proof
    record["qualification_birth_seed"] = (
        "leading seed if its original balance/profile checks close, otherwise existing native finite-birth seed; "
        "selected before the one physical solve, uniformly at every mesh/time level"
    )
    write_json(path, record)
    if record.get("error") or not path.with_suffix(".accepted.pkl").exists():
        raise RuntimeError("fresh Faner birth refused: " + str(record.get("error")))
    return path


def apply_ruling_options(args):
    """Process-wide declared defaults of the 2026-10-17 numerical gate ruling."""
    engine.driver.PROOF_WALL_CLOCK_CLAUSE = bool(getattr(args, "proof_wall_clock_clause", False))
    engine.study.numerics.CERTIFICATE_JACOBIAN_OUTSIDE_SEARCH_BUDGET = not bool(
        getattr(args, "certificate_within_search_budget", False)
    )
    import particle_imbibition_bridge as imbibition

    mode = getattr(args, "imbibition_bridge", "ruled")
    w_sat = getattr(args, "saturated_loading", imbibition.W_SAT_NOMINAL)
    if imbibition.active() is None:
        imbibition.install(mode, w_sat)
    return {
        "imbibition_bridge": {
            "mode": imbibition.active()[0],
            "saturated_loading_kg_kg_dry": imbibition.active()[1],
            "declared_rows": imbibition.declared_rows(*imbibition.active()),
        },
        "ruling": "GT_PS2_OWNER_RULING_PAPER1_NUMERICAL_GATE_HYGIENE_2026-10-17",
        "dry_local_residual_tolerance": args.dry_local_tolerance,
        "proof_wall_clock_clause": engine.driver.PROOF_WALL_CLOCK_CLAUSE,
        "certificate_jacobian_outside_search_budget": (
            engine.study.numerics.CERTIFICATE_JACOBIAN_OUTSIDE_SEARCH_BUDGET
        ),
        "gate_units": gate_units_declared_rows(),
    }


def gate_units_declared_rows():
    """Typed rows of owner ruling GT_PS2_OWNER_RULING_PAPER1_GATE_UNITS_2026-10-17."""
    import experiment_core2_dt_dry as dry

    birth = engine.driver
    common = {
        "kind": "DECLARED_ENGINEERING_ASSUMPTION",
        "authority": dry.GATE_UNITS_RULING,
        "physically_qualifying": False,
    }
    return {
        "ruling": dry.GATE_UNITS_RULING,
        "declared_rows": [
            {
                **common,
                "name": "dry_advance_search_target_equals_gate",
                "search_target": "the declared dry tolerance (config nonlinear_residual_tolerance); "
                "a differing search_tolerance is refused (DrySearchTargetRefused)",
                "former": "min(5e-11, tolerance) with a tighten-only override",
                "cell_component_rows": "weighted to the strict native component limit "
                "min(2e-11, tolerance), as before",
            },
            {
                **common,
                "name": "dry_energy_row_noise_aware_gate",
                "gate": "max(tolerance, m * floor_i), floor_i = sum_j |J_ij| ulp(x_j) at the iterate",
                "declared_multiple": dry.NOISE_AWARE_ENERGY_GATE_MULTIPLE,
                "rows": "cell energy rows of the dry advance; water, hexane and surface node rows keep 1e-10",
                "jacobian": "the domain-aware Newton's dense forward-difference Jacobian at the returned iterate",
                "stamped": "dry.history.jsonl native_report.noise_aware_energy_row_gate per row; "
                "case.json noise_aware_energy_gate_summary",
            },
            {
                **common,
                "name": "birth_condition_certificate_on_searched_matrix",
                "gated": "np.linalg.cond of the column-2-norm-scaled physical joint birth Jacobian, DT and Faner",
                "maximum_birth_condition_proxy": birth.MAXIMUM_BIRTH_CONDITION_PROXY,
                "unchanged_maximum_condition_proxy_elsewhere": 1e12,
                "reported_not_gated": "raw-column physical proxy, native pullback condition",
            },
            {
                "kind": "DECLARED_ENGINEERING_ASSUMPTION",
                "name": "dry_line_search_merit_in_gate_units",
                "authority": "lead decision (Class A) under the hygiene principle of "
                "GT_PS2_OWNER_RULING_PAPER1_NUMERICAL_GATE_HYGIENE_2026-10-17",
                "record": "docs/GT_PS2_PAPER1_LINE_SEARCH_MERIT_BUILD_2026-10-18.md",
                "physically_qualifying": False,
                "merit": "the domain-aware Newton's sufficient-decrease test on the gate-normalized "
                "search residual r_i / g_i whenever an acceptance clause is passed (the dry advance): "
                "g_i = max(tolerance, m * floor_i) on the cell energy rows, floor_i from the Newton's "
                "Jacobian at the current iterate, and the fixed gate (the search target) on every other row",
                "declared_norm": engine.study.numerics.GATE_NORMALIZED_MERIT_NORM,
                "norm_choice_evidence": "paper/analysis/results_2026-10-16/solver_qualification/"
                "line_search_merit_150_*.json (both norms measured on the refused DT_N64_L0 state)",
                "former": "the unweighted 2-norm of the search residual (evidence 149: refused the "
                "full step the per-row clause accepts, on energy-row rounding noise)",
                "added_residual_calls": 0,
                "unchanged": "gates, targets, the 1,400-call allowance, the 1e-10 tolerance, the "
                "Armijo constant 1e-4, the damping ladder, the search direction, every law; callers "
                "without an acceptance clause (birth, face, projection, extinction) bit for bit",
                "stamped": "dry.history.jsonl native_report.line_search (merit, norm, per accepted "
                "step damping, step_length, retreats, merit before and after)",
            },
        ],
    }


def update_noise_aware_gate_summary(record, report):
    """Case-level summary of the stamped per-row energy gates of accepted dry advances."""
    stamp = (report or {}).get("noise_aware_energy_row_gate")
    if stamp is None:
        return
    summary = record.setdefault(
        "noise_aware_energy_gate_summary",
        {
            "ruling": stamp["ruling"],
            "declared_multiple": stamp["declared_multiple"],
            "accepted_advances": 0,
            "advances_with_absorbed_rows": 0,
            "absorbed_rows_total": 0,
            "maximum_floor": 0.0,
            "maximum_gate": 0.0,
            "maximum_residual_over_gate": 0.0,
            "maximum_absorbed_residual": 0.0,
            "first_absorbed_time_s": None,
        },
    )
    summary["accepted_advances"] += 1
    summary["absorbed_rows_total"] += stamp["absorbed_rows"]
    if stamp["absorbed_rows"]:
        summary["advances_with_absorbed_rows"] += 1
        if summary["first_absorbed_time_s"] is None:
            summary["first_absorbed_time_s"] = report.get("time_s")
        summary["maximum_absorbed_residual"] = max(
            summary["maximum_absorbed_residual"],
            max(abs(r) for r, a in zip(stamp["residual"], stamp["absorbed_by_floor"]) if a),
        )
    summary["maximum_floor"] = max(summary["maximum_floor"], max(stamp["floor"]))
    summary["maximum_gate"] = max(summary["maximum_gate"], max(stamp["gate"]))
    summary["maximum_residual_over_gate"] = max(
        summary["maximum_residual_over_gate"], stamp["maximum_residual_over_gate"]
    )


def dry_physical_advance_assembly_budget(current):
    """One dry physical advance's shared residual-assembly allowance.

    Knob off: the former fixed 1,400 (the state is not read). Kernel-envelope item 7 armed:
    3 (4 u + 6) with u = 3 N + 1, N the dry grid's cells.
    """
    if not envelope.armed(envelope.DRY_ALLOWANCE_PER_UNKNOWN):
        return {"limit": 1400, "used": 0}
    unknowns = envelope.dry_physical_advance_unknowns(current.grid.n)
    return {"limit": envelope.dry_assembly_allowance(unknowns), "used": 0}


def update_native_factor_ceiling_summary(record, report):
    """Kernel-envelope item 7 (armed only): count the report-only ceiling stamps."""
    if not envelope.armed(envelope.NATIVE_LEDGER_FACTOR_CEILING):
        return
    summary = record.setdefault(
        "kernel_envelope_native_ledger_factor_ceiling",
        {
            "ceiling": envelope.NATIVE_LEDGER_FACTOR_REPORTING_CEILING,
            "report_only": True,
            "accepted_advances": 0,
            "advances_above_ceiling": 0,
            "maximum_factor": 0.0,
            "first_above_time_s": None,
        },
    )
    stamp = (report or {}).get("noise_aware_energy_row_gate") or {}
    factors = stamp.get("native_ledger_energy_scale_factors") or []
    summary["accepted_advances"] += 1
    if factors:
        summary["maximum_factor"] = max(summary["maximum_factor"], max(factors))
    if stamp.get("kernel_envelope_native_ledger_factor_ceiling") is not None:
        summary["advances_above_ceiling"] += 1
        if summary["first_above_time_s"] is None:
            summary["first_above_time_s"] = report.get("time_s")


def run(args):
    time_control = HistoryTimeControl(
        args.study, args.time_level, getattr(args, "dry_time_schedule", "uniform")
    )
    args.out.mkdir(parents=True)
    pinned = hashes()
    multiplier = time_control.multiplier
    record = {
        "status": "RUNNING",
        "study": args.study,
        "base_cells": args.cells,
        "time_level": args.time_level,
        "history_time_control": time_control.describe(),
        "native_dry_local_residual_tolerance": args.dry_local_tolerance,
        "conservative_original_feed_components": args.conservative_components,
        "source_sha256_at_entry": pinned,
        "property_parameters": asdict(ParticleModelParameters()),
        "per_area_interface_conditioning": args.experimental_per_area,
        "experimental_per_area_conditioning": args.experimental_per_area
        and not args.per_area_qualification,
        "owner_adoption_required": args.experimental_per_area and not args.per_area_qualification,
        "conditioning_authority": (
            front.AREA_CHART_AUTHORITY if args.per_area_qualification else None
        ),
        "physically_qualifying": False,
        "convergence_established": False,
        "fresh_original_feed": True,
        "single_runtime_source_generation": True,
        "complete_history": False,
        "C_accepted_steps": 0,
        "numerical_gate_ruling_settings": apply_ruling_options(args),
        # The adopted film-covered surface law and its one declared contact conductance.
        "surface_contact_law": envelope.surface_contact_law_declared_row(),
    }
    if envelope.armed_items():
        record["kernel_envelope_hygiene"] = {
            "record": envelope.RECORD,
            "armed": list(envelope.armed_items()),
            "declared_rows": envelope.declared_rows(),
            "physically_qualifying": False,
        }
    if args.study == "Faner" and faner_trace.armed():
        # Kernel-envelope FANER_MEASURED_VAPOUR_TRACE (armed only): the trace file by path,
        # its sha256, the points and the rule.
        record["faner_measured_vapour_trace"] = faner_trace.declared_row()
    if getattr(args, "grouped_birth_search", False):
        record["declared_birth_search_option"] = "--grouped-birth-search"
    memo_option = (
        {"memoize_identical_points": True} if getattr(args, "memoize_proof_points", False) else {}
    )
    if memo_option:
        record["declared_proof_option"] = "--memoize-proof-points"
    if getattr(args, "front_search_at_acceptance_gate", False):
        record["declared_front_search_option"] = "--front-search-at-acceptance-gate"
    if getattr(args, "conservative_core_endpoint", False):
        record["declared_core_endpoint_option"] = "--conservative-core-endpoint"
    if args.bath:
        record["bath_sha256"] = hashlib.sha256(args.bath.read_bytes()).hexdigest()
        # Kernel-envelope INSTRUMENT_SHARED_LAW (armed only): the bath gas's water activity on
        # the shared law beside the ideal Raoult-Dalton reading, report only.
        instrument = envelope.bath_instrument(json.loads(args.bath.read_text()))
        if instrument is not None:
            record["kernel_envelope_instrument_shared_law"] = instrument
    started = time.perf_counter()
    result_path = args.out / "case.json"
    write_json(result_path, record)
    # Kernel-envelope SURFACE_GRADED_MESH (armed only): the declared face set, stamped.
    face_set = getattr(args, "surface_graded_face_set", None)
    if face_set is not None:
        record["surface_graded_mesh"] = face_set.describe(
            0.0015 if args.study == "DT" else 0.000885
        )
        write_json(result_path, record)
    try:
        with graded_mesh.installed(face_set) if face_set is not None else nullcontext():
            birth = (
                fresh_dt_birth(args, multiplier, pinned)
                if args.study == "DT"
                else fresh_faner_birth(args, multiplier)
            )
        record["birth_result"] = str(birth)
        write_json(result_path, record)
        horizon = 1766.5355008060146 if args.study == "DT" else 240.0
        b_path = args.out / "front.json"
        b_args = SimpleNamespace(
            resume=birth.with_suffix(".accepted.pkl"),
            out=b_path,
            study=args.study,
            bath=args.bath,
            maximum_dt=(0.1 if args.study == "DT" else 0.4) / multiplier,
            minimum_dt=1e-5,
            until=horizon,
            steps=10000,
            core_time_tolerance=0.002 / multiplier,
            experimental_per_area=args.experimental_per_area,
            per_area_qualification=args.per_area_qualification,
        )
        if getattr(args, "memoize_proof_points", False):
            b_args.memoize_proof_points = True
        if getattr(args, "front_search_at_acceptance_gate", False):
            b_args.front_search_at_acceptance_gate = True
        b_result = front.run(b_args)
        if b_result.get("error") or not b_result.get("native_endpoint_inputs_ready"):
            raise RuntimeError(
                "front history has no accepted native endpoint inputs: "
                + str(b_result.get("error"))
            )
        checkpoint = pickle.loads(b_path.with_suffix(".accepted.pkl").read_bytes())
        current, carry, surface = (checkpoint[k] for k in ("current", "carry", "surface"))
        config = current.transport.config
        if args.study == "DT":
            bath = json.loads(args.bath.read_text())

            def boundary_factory(clock):
                return engine.face.bath_boundary(
                    bath,
                    (
                        current.transport.geometry.master_grid.R
                        if hasattr(current, "transport")
                        else current.grid.R
                    ),
                    clock,
                )

        else:

            def boundary_factory(clock):
                # Kernel-envelope FANER_MEASURED_VAPOUR_TRACE (armed only): the trace at the
                # clock; off, the constant boundary itself.
                return faner_trace.bath_boundary(checkpoint["boundary"], clock)

        from contextlib import ExitStack

        with ExitStack() as scope:
            if args.study == "Faner":
                boundary = checkpoint["boundary"]
                scope.enter_context(
                    faner.declared_faner_film_multiplier(boundary.film_boundary.stress_multiplier)
                )
                scope.enter_context(faner.pure_bulk_domain_adapter())
            else:
                # The extra liquid path is present during the native core
                # projection; dry active sets own their own scoped path.
                pass
            scope.enter_context(engine.study.magnitude_storage_scale())
            scope.enter_context(engine.driver.path_bounds.monotone_water_path_enclosure([]))
            scope.enter_context(engine.driver.path_bounds.monotone_hexane_path_enclosure([]))
            if args.study == "DT":
                with (
                    engine.finite_shell_path(100, envelope.surface_contact_conductance()),
                    (
                        engine.film_free_surface_scope()
                        if checkpoint.get("kernel_envelope_film_free", False)
                        else nullcontext()
                    ),
                    engine.driver.bounded_path_proof(3000000, 120, **memo_option) as proof,
                ):
                    endpoint = project_core_endpoint(
                        current,
                        carry,
                        surface,
                        checkpoint["positive_core_targets"],
                        checkpoint["missing_face"],
                        boundary_factory,
                        maximum_remaining_time_s=0.002 / multiplier,
                        **(
                            {"conservative_surface": True}
                            if getattr(args, "conservative_core_endpoint", False)
                            else {}
                        ),
                    )
            else:
                with engine.driver.bounded_path_proof(3000000, 120, **memo_option) as proof:
                    endpoint = project_core_endpoint(
                        current,
                        carry,
                        surface,
                        checkpoint["positive_core_targets"],
                        checkpoint["missing_face"],
                        boundary_factory,
                        maximum_remaining_time_s=0.002 / multiplier,
                        **(
                            {"conservative_surface": True}
                            if getattr(args, "conservative_core_endpoint", False)
                            else {}
                        ),
                    )
            current, carry, surface = endpoint.after, endpoint.carry, endpoint.surface
            write_json(
                args.out / "endpoint.json",
                {
                    "accepted": True,
                    "ledger": asdict(endpoint.native_step.ledger),
                    "time_control": asdict(endpoint.time_control),
                    "proof": proof,
                    "base_face_restoration": endpoint.base_face_restoration,
                    # the residual-contract polish's stamp, only where it acted
                    **(
                        {"core_endpoint_residual_polish": endpoint.residual_polish}
                        if endpoint.residual_polish is not None
                        else {}
                    ),
                    "after_transport": asdict(current),
                    "surface": asdict(surface),
                },
            )
            original_balance = None
            if args.conservative_components:
                from particle_original_feed_balance import prepare_prefix

                original_balance, prefix_report = prepare_prefix(
                    args.out, args.study, args.bath, hashes()
                )
                initial_actual = liquid.row(current, surface, carry)["combined_W_H_U"]
                if (
                    original_balance.clock != current.time_s
                    or tuple(initial_actual) != original_balance.last_actual
                ):
                    raise RuntimeError(
                        "original-feed prefix does not match the exact C activation state"
                    )
                write_json(args.out / "original_feed_prefix_audit.json", prefix_report)
                record["original_feed_prefix_audit_passed"] = True
            c_path = args.out / "dry.history.jsonl"
            with c_path.open("x", encoding="utf-8") as stream:
                origin_row = liquid.row(current, surface, carry)
                origin = origin_row["combined_W_H_U"]
                scales = [origin[0], origin[1], carry["import_capacity_scale_j"]]
                integrated = [[], [], []]
                previous_phase = "liquid" if carry["external_water_loading_kg_m2"] > 0 else "absent"
                film_points = [(current.time_s, carry["external_water_loading_kg_m2"])]
                previous_profile = None
                declared_options = {}
                if getattr(args, "event_first_film_depletion", False):
                    declared_options.update(
                        film_depletion_event_first=True,
                        stop_positive_film_at_nonpositive_liquid=True,
                    )
                if getattr(args, "eliminate_event_absent_components", False):
                    declared_options.update(
                        conservative_event_elimination=True,
                        conservative_absent_elimination=True,
                    )
                if declared_options:
                    record["declared_film_depletion_search_options"] = dict(declared_options)
                if getattr(args, "curvature_corrected_depletion_prediction", False):
                    record["declared_depletion_prediction"] = "curvature_corrected_three_point"
                component_budgets = None
                # The declared maximum dry physical advances. Kernel-envelope BUDGETS_PER_UNKNOWN
                # (ruling packet 6 item 5, GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md;
                # armed only): max(1,500, 2 P), P the planned requests from this clock to the
                # horizon; the same count divides the per-advance component budgets.
                maximum_advances = declared_dry_advance_ceiling(args)
                if maximum_advances != DECLARED_DRY_ADVANCE_CEILING:
                    record["history_time_control"]["dry_physical_advance_limit"] = maximum_advances
                    record["declared_dry_advance_ceiling"] = {
                        "value": maximum_advances,
                        "was": DECLARED_DRY_ADVANCE_CEILING,
                        "ruling": "owner, 2026-10-28: raise the ceiling, run the half-step pair",
                        "physically_qualifying": False,
                    }
                budgets_rule = envelope.budgets_per_unknown()
                if budgets_rule is not None:
                    maximum_advances, ceiling_stamp = envelope.dry_advance_ceiling(
                        time_control, current.time_s, horizon, maximum_advances
                    )
                    record["kernel_envelope_budgets_per_unknown"] = {
                        "declared_row": envelope.budgets_per_unknown_declared_row(budgets_rule),
                        "dry_advance_ceiling": ceiling_stamp,
                    }
                    record["history_time_control"]["dry_physical_advance_limit"] = maximum_advances
                if original_balance is not None:
                    from particle_original_feed_balance import component_error_budgets_per_advance

                    component_budgets = component_error_budgets_per_advance(
                        original_balance, origin, maximum_advances
                    )
                    record["global_component_error_budgets_kg_per_advance"] = component_budgets
                stream.write(
                    json.dumps(
                        {**origin_row, "kind": "initial", "after_transport": asdict(current)},
                        default=front.json_default,
                    )
                    + "\n"
                )
                begin_advance = scope.enter_context(
                    engine.driver.physical_advance_proof_budget(3000000, 120, **memo_option)
                )
                for index in range(maximum_advances):
                    if current.time_s >= horizon:
                        record["complete_history"] = True
                        break
                    proof = begin_advance()
                    assembly_budget = dry_physical_advance_assembly_budget(current)
                    request_start_clock = current.time_s
                    # Film-aware last band (ruling packet 7, 2026-10-27): the caller states the film.
                    requested = time_control.dry_request(
                        current.time_s,
                        horizon,
                        film_present=carry["external_water_loading_kg_m2"] > 0,
                    )
                    advance_options = dict(declared_options)
                    if getattr(args, "event_first_film_depletion", False):
                        advance_options["predicted_depletion_duration_s"] = (
                            curvature_corrected_depletion_duration(film_points)
                            if getattr(args, "curvature_corrected_depletion_prediction", False)
                            else predicted_film_depletion_duration(film_points)
                        )
                        advance_options["previous_accepted_profile"] = previous_profile
                    # The depletion advance's allowance doubled (ruling packet 7, 2026-10-27).
                    assembly_budget["limit"] *= depletion_advance_allowance_factor(
                        carry["external_water_loading_kg_m2"],
                        advance_options.get("predicted_depletion_duration_s"),
                        requested,
                    )
                    # Kernel-envelope DRY_HALF_STEP_FIRST (armed only): a film-bearing advance
                    # abandons its full step at the declared assembly threshold.
                    half_step_first = (
                        envelope.armed(envelope.DRY_HALF_STEP_FIRST)
                        and carry["external_water_loading_kg_m2"] > 0
                    )
                    half_step_stamp = None
                    # The adopted film-covered surface law: report-only stamp of every advance.
                    liquid.surface_contact_reset()
                    # LIQUID_MOBILITY_RULE (study knob, armed only): report-only stamp.
                    mobility_rule = envelope.liquid_mobility_rule()
                    if mobility_rule is not None:
                        liquid.surface_rules_reset()
                    for halving in range(13):
                        if half_step_first:
                            if halving == 0:
                                # BUDGETS_PER_UNKNOWN (ruling packet 6 item 4): max(1,200,
                                # floor(5 A / 12)), A this advance's allowance in force.
                                assembly_budget[envelope.HALF_STEP_ABANDON_KEY] = (
                                    envelope.DRY_HALF_STEP_FIRST_ASSEMBLY_THRESHOLD
                                    if budgets_rule is None
                                    else envelope.dry_half_step_threshold(assembly_budget["limit"])
                                )
                            else:
                                assembly_budget.pop(envelope.HALF_STEP_ABANDON_KEY, None)
                        if args.study == "DT" and envelope.armed(envelope.PRESSURE_SCHEDULE):
                            advance_options["target_pressure_pa"] = envelope.bath_pressure_pa(
                                bath, current.time_s + requested / (2**halving)
                            )
                        try:
                            candidate = advance_dry_phases(
                                current,
                                config,
                                carry,
                                surface,
                                requested / (2**halving),
                                boundary_factory,
                                previous_phase=previous_phase,
                                water_inventory_scale_kg=scales[0],
                                assembly_budget=assembly_budget,
                                global_component_mass_error_budgets_kg=component_budgets,
                                conservative_liquid_elimination=args.conservative_components,
                                **advance_options,
                            )
                        except envelope.DryFullStepAbandoned as abandoned:
                            if not (half_step_first and halving == 0):
                                raise
                            assembly_budget.pop(envelope.HALF_STEP_ABANDON_KEY, None)
                            half_step_stamp = {
                                "abandoned_full_step_s": requested,
                                "assemblies_at_abandonment": abandoned.used,
                                "declared_numerics_threshold": abandoned.threshold,
                                "allowance": assembly_budget["limit"],
                                "film_kg_m2_at_start": carry["external_water_loading_kg_m2"],
                                "record": envelope.QUEUE_RECORD,
                                "physically_qualifying": False,
                            }
                            stream.write(
                                json.dumps(
                                    {
                                        "kind": "refusal",
                                        "time_s": current.time_s,
                                        "halving": halving,
                                        "trials": [],
                                        "kernel_envelope_dry_half_step_first": half_step_stamp,
                                    },
                                    default=front.json_default,
                                )
                                + "\n"
                            )
                            stream.flush()
                            record["kernel_envelope_dry_half_step_first_abandonments"] = (
                                record.get("kernel_envelope_dry_half_step_first_abandonments", 0)
                                + 1
                            )
                            continue
                        except StepBudgetExceeded:
                            raise
                        except DryPhaseAdvanceError as error:
                            stream.write(
                                json.dumps(
                                    {
                                        "kind": "refusal",
                                        "time_s": current.time_s,
                                        "halving": halving,
                                        "trials": error.rejected_trials,
                                    },
                                    default=front.json_default,
                                )
                                + "\n"
                            )
                            stream.flush()
                            if halving == 12:
                                raise
                            continue
                        if half_step_first:
                            assembly_budget.pop(envelope.HALF_STEP_ABANDON_KEY, None)
                        row = liquid.row(candidate.after, candidate.surface, candidate.carry)
                        elapsed = candidate.after.time_s - current.time_s
                        histories = [history.copy() for history in integrated]
                        for history, flux in zip(
                            histories,
                            (
                                candidate.surface.external_water_flux_kg_m2_s,
                                candidate.surface.external_hexane_flux_kg_m2_s,
                                candidate.surface.external_energy_flux_w_m2,
                            ),
                        ):
                            history.append(elapsed * candidate.after.grid.areas[-1] * flux)
                        defects = [
                            math.fsum((a, -b, *history))
                            for a, b, history in zip(row["combined_W_H_U"], origin, histories)
                        ]
                        normalized = [abs(d) / s for d, s in zip(defects, scales)]
                        if max(normalized) > config.dry.ledger_tolerance:
                            raise RuntimeError(
                                "independent C primitive/rate cumulative conservation refused: "
                                + str(normalized)
                            )
                        if original_balance is not None:
                            original_balance.accept(
                                candidate.after.time_s,
                                row["combined_W_H_U"],
                                [
                                    candidate.after.grid.areas[-1] * getattr(candidate.surface, k)
                                    for k in (
                                        "external_water_flux_kg_m2_s",
                                        "external_hexane_flux_kg_m2_s",
                                        "external_energy_flux_w_m2",
                                    )
                                ],
                                "native conservative C " + candidate.phase,
                            )
                        previous_profile = (
                            current.time_s,
                            list(current.temperatures_k),
                            list(current.y_hexane),
                        )
                        current, carry, surface = (
                            candidate.after,
                            candidate.carry,
                            candidate.surface,
                        )
                        if current.config.pressure_pa != config.dry.pressure_pa:
                            # Kernel-envelope pressure schedule: the cut configuration
                            # follows the accepted state's pressure-only target.
                            config = engine.core.cut.pressure_only_target_config(
                                config, current.config.pressure_pa
                            )
                        integrated, previous_phase = histories, candidate.phase
                        film_points.append((current.time_s, carry["external_water_loading_kg_m2"]))
                        row.update(
                            kind="depletion" if candidate.depletion_event else candidate.phase,
                            time_control_request={
                                "start_clock_s": request_start_clock,
                                "requested_dt_s": requested,
                                "accepted_halving": halving,
                                "native_attempt_dt_s": requested / (2**halving),
                                "schedule": time_control.dry_schedule,
                            },
                            after_transport=asdict(current),
                            native_report=candidate.native_report,
                            independent_normalized_residuals=normalized,
                            assembly_work=assembly_budget,
                            rejected_trials=candidate.rejected_trials,
                            proof=dict(proof),
                            independent_original_feed_normalized=(
                                original_balance.rows[-1]["normalized_residuals"]
                                if original_balance is not None
                                else None
                            ),
                        )
                        if half_step_stamp is not None:
                            row["kernel_envelope_dry_half_step_first"] = {
                                **half_step_stamp,
                                "accepted_halving": halving,
                            }
                        row["kernel_envelope_finite_contact"] = (
                            liquid.surface_contact_advance_stamp(
                                candidate.native_report,
                                candidate.liquid_supply_audit,
                            )
                        )
                        if mobility_rule is not None:
                            row["kernel_envelope_surface_rules"] = (
                                liquid.surface_rules_advance_stamp(
                                    mobility_rule,
                                    candidate.native_report,
                                    current,
                                    config,
                                )
                            )
                        stream.write(json.dumps(row, default=front.json_default) + "\n")
                        stream.flush()
                        record["C_accepted_steps"] += 1
                        update_noise_aware_gate_summary(record, candidate.native_report)
                        update_native_factor_ceiling_summary(record, candidate.native_report)
                        record["last_accepted"] = {
                            k: v
                            for k, v in row.items()
                            if k
                            not in (
                                "after_transport",
                                "surface",
                                "native_report",
                                "proof",
                                "rejected_trials",
                            )
                        }
                        checkpoint.update(
                            current=current,
                            config=config,
                            carry=carry,
                            surface=surface,
                            previous_surface_phase=previous_phase,
                            C_origin=origin,
                            C_streams=integrated,
                            checkpoint_codec=(
                                "paper1-complete-original-feed-case-v2-original-feed-ledger"
                                if original_balance is not None
                                else "paper1-complete-original-feed-case-v1"
                            ),
                            source_sha256=pinned,
                        )
                        if original_balance is not None:
                            checkpoint["original_feed_conservation"] = {
                                "codec": "original-feed-component-conservation-v1",
                                "initial": original_balance.initial,
                                "scales": original_balance.scales,
                                "clock": original_balance.clock,
                                "integral_terms": original_balance.streams,
                                "maximum_normalized_W_H_U": original_balance.maximum,
                                "last_actual_W_H_U": original_balance.last_actual,
                                "last_signed_defects_W_H_U": original_balance.last_signed_defects,
                            }
                            record["maximum_original_feed_normalized_W_H_U"] = (
                                original_balance.maximum
                            )
                        temporary = args.out / "accepted.tmp"
                        temporary.write_bytes(pickle.dumps(checkpoint, protocol=5))
                        front.replace_with_retry(
                            temporary,
                            args.out / "accepted.pkl",
                            record.setdefault("checkpoint_replace_retries", []),
                        )
                        write_json(result_path, record)
                        if record["C_accepted_steps"] % 25 == 0 or candidate.depletion_event:
                            print(
                                json.dumps(
                                    {
                                        "C_steps": record["C_accepted_steps"],
                                        "time_s": current.time_s,
                                        "phase": previous_phase,
                                    }
                                ),
                                flush=True,
                            )
                        break
                else:
                    raise RuntimeError("declared C physical-advance budget exhausted")
                if current.time_s == horizon:
                    record["complete_history"] = True
        if record["complete_history"]:
            record["status"] = (
                "COMPLETE_CANDIDATE_PENDING_CONVERGENCE_AND_INDEPENDENT_WHOLE_HISTORY_AUDIT"
            )
    except Exception as error:
        record.update(
            status="REFUSED_WITH_ACCEPTED_EVIDENCE_PRESERVED",
            error=str(error),
            traceback=traceback.format_exc(),
            optimizer_debug=getattr(error, "optimizer_debug", None),
        )
    record["wall_seconds"] = time.perf_counter() - started
    record["source_sha256_at_exit"] = hashes()
    record["sources_unchanged_during_run"] = pinned == record["source_sha256_at_exit"]
    write_json(result_path, record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", choices=("DT", "Faner"), required=True)
    parser.add_argument("--cells", type=int, required=True)
    parser.add_argument("--time-level", type=int, choices=(0, 1, 2), required=True)
    parser.add_argument(
        "--dry-advance-ceiling",
        type=int,
        default=DECLARED_DRY_ADVANCE_CEILING,
        help="the dry stage's declared maximum physical advances (default 1,500, a positive "
        "integer; the same count divides the per-advance component budgets); 4,000 for the DT "
        "half-step pair under the owner's ruling of 2026-10-28",
    )
    parser.add_argument(
        "--dry-time-schedule",
        choices=("uniform", "dt-film-resolved"),
        default="uniform",
        help="explicit DT qualification family: caps 0.3125 s before 515 s, "
        "0.625 s to 570 s, then 5 s, all divided by 2**time-level; "
        "uniform preserves the existing controls",
    )
    parser.add_argument("--bath", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--experimental-per-area", action="store_true")
    parser.add_argument("--per-area-qualification", action="store_true")
    # Owner ruling GT_PS2_OWNER_RULING_PAPER1_NUMERICAL_GATE_HYGIENE_2026-10-17:
    # 1e-10 is the declared default; 2e-11 remains reachable for replay.
    parser.add_argument("--dry-local-tolerance", type=float, default=1e-10)
    parser.add_argument(
        "--imbibition-bridge",
        choices=("ruled", "former-bridge", "plain-cap"),
        default="ruled",
        help="GT_PS2_OWNER_RULING_PAPER1_FILM_IMBIBITION_BRIDGE_2026-10-17: retained water above "
        "the Luikov cap on the particle path (ruled default); former-bridge = C8 tier-1 linear "
        "bridge at 0.25; plain-cap = the former capacity active set (replay)",
    )
    parser.add_argument(
        "--saturated-loading",
        type=float,
        default=0.55,
        help="W_sat, the loading at a_w = 1 (kg/kg dry), ruled sweep [0.25, 0.97]",
    )
    parser.add_argument(
        "--proof-wall-clock-clause",
        action="store_true",
        help="replay option: restore the former 120 s wall-clock clause of the proof allowance",
    )
    parser.add_argument(
        "--certificate-within-search-budget",
        action="store_true",
        help="replay option: count the certificate Jacobian within the residual-call budget",
    )
    parser.add_argument("--conservative-components", action="store_true")
    parser.add_argument(
        "--event-first-film-depletion",
        action="store_true",
        help="declared search option: predicted exact-zero film event first, "
        "and stop the eliminated positive-film search at a nonpositive film",
    )
    parser.add_argument(
        "--curvature-corrected-depletion-prediction",
        action="store_true",
        help="declared search option: seed the event-first depletion with the parabola "
        "through the last three accepted films instead of the last two (seeding only)",
    )
    parser.add_argument(
        "--eliminate-event-absent-components",
        action="store_true",
        help="declared search option: algebraic component elimination for the "
        "depletion-event and absent-film active sets (needs --conservative-components)",
    )
    parser.add_argument(
        "--conservative-core-endpoint",
        action="store_true",
        help="declared option: core-extinction projection with the algebraically eliminated "
        "positive-film surface; every native and combined ledger gate unchanged",
    )
    parser.add_argument(
        "--front-search-at-acceptance-gate",
        action="store_true",
        help="declared search option: the front driver searches to the native acceptance "
        "gate (controls.nonlinear_residual_tolerance) instead of the tighter 1e-11; every "
        "native acceptance check unchanged",
    )
    parser.add_argument(
        "--memoize-proof-points",
        action="store_true",
        help="declared search option: within one proof allowance, reuse bit-identical native "
        "path-coordinate evaluations; every requested point still counts against the "
        "unchanged 3,000,000-point and 120 s clauses",
    )
    parser.add_argument(
        "--grouped-birth-search",
        action="store_true",
        help="declared search option (Faner): grouped-colour Newton search directions for "
        "the joint birth; dense final certificate, budget and cap unchanged",
    )
    parser.add_argument(
        "--graded-interior-cells",
        type=int,
        help="kernel-envelope SURFACE_GRADED_MESH only: the uniform interior cell count m",
    )
    parser.add_argument(
        "--graded-band-cells",
        type=int,
        help="kernel-envelope SURFACE_GRADED_MESH only: the surface band's cell count n",
    )
    parser.add_argument(
        "--graded-band-ratio",
        type=float,
        help="kernel-envelope SURFACE_GRADED_MESH only: the band's geometric ratio q > 1",
    )
    args = parser.parse_args()
    try:
        HistoryTimeControl(args.study, args.time_level, args.dry_time_schedule)
    except ValueError as error:
        parser.error(str(error))
    try:
        envelope.validate_environment()
    except envelope.EnvelopeKnobRefused as error:
        parser.error(str(error))
    if faner_trace.armed() and args.study != "Faner":
        parser.error(
            "DTDC_KERNEL_ENVELOPE_FANER_MEASURED_VAPOUR_TRACE is declared for the Faner study "
            "only (the DT bath is the --bath file)"
        )
    graded_values = (args.graded_interior_cells, args.graded_band_cells, args.graded_band_ratio)
    args.surface_graded_face_set = None
    if envelope.armed(envelope.SURFACE_GRADED_MESH):
        if any(value is None for value in graded_values):
            parser.error(
                "DTDC_KERNEL_ENVELOPE_SURFACE_GRADED_MESH needs --graded-interior-cells, "
                "--graded-band-cells and --graded-band-ratio"
            )
        try:
            args.surface_graded_face_set = graded_mesh.parse(*graded_values)
            graded_mesh.native_grid(
                args.surface_graded_face_set, 0.0015 if args.study == "DT" else 0.000885
            )
        except graded_mesh.SurfaceGradedMeshRefused as error:
            parser.error(str(error))
        if args.surface_graded_face_set.cells != args.cells:
            parser.error("--cells must equal the graded total m + n")
    elif any(value is not None for value in graded_values):
        parser.error(
            "the --graded-* options are refused without DTDC_KERNEL_ENVELOPE_SURFACE_GRADED_MESH"
        )
    if args.grouped_birth_search and args.study != "Faner":
        parser.error("--grouped-birth-search is declared for the Faner birth only")
    if args.eliminate_event_absent_components and not args.conservative_components:
        parser.error("--eliminate-event-absent-components requires --conservative-components")
    args.experimental_per_area = args.experimental_per_area or args.per_area_qualification
    args.out = (INVOCATION_CWD / args.out).resolve()
    if args.bath:
        args.bath = (INVOCATION_CWD / args.bath).resolve()
    if args.out.exists() or args.cells < 2 or (args.study == "DT" and args.bath is None):
        parser.error("use a new output directory, at least two base cells, and an explicit DT bath")
    result = run(args)
    print(
        json.dumps(
            {k: result.get(k) for k in ("status", "complete_history", "error", "wall_seconds")}
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
