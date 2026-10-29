"""Whole-front engineering experiment using conservative finite-volume faces.

Native phase, ALE/RH, component/energy and conditioning gates own every
physical advance. Mesh operations advance no time. Accepted/rejected trials,
primitive fields, old references and the controller persist across restarts.
This candidate needs full-history mesh/time refinement before qualification.
"""

import argparse
import contextlib
import faner_measured_trace as faner_trace
from contextlib import ExitStack
from dataclasses import asdict, fields, is_dataclass, replace
import hashlib
import json
import math
from pathlib import Path
import pickle
import struct
import time
import traceback

import numpy as np
import experiment_core2_faner as faner
import experiment_core2_finite_wetting_b as engine
from face_volume_coalescence import coalesce_arrival_face, slide_arrival_face
from particle_solver_configuration import PAPER_CHECKPOINT_CODEC, ParticleModelParameters
from particle_core_extinction import (
    PositiveCoreTargetEvidence,
    estimate_core_time,
    CoreTimeNotResolved,
)
from dtdc_simulator.core2.particle.numerical_work import StepBudgetExceeded
import particle_kernel_envelope as envelope
import particle_kernel_pd_boundary as pd_boundary
import particle_kernel_stall_materiality as materiality

CODEC = "paper1-conservative-sliding-front-v3-approved-chart-control"
AREA_CHART_AUTHORITY = "docs/GT_PS2_PAPER1_PER_AREA_QUALIFICATION_OWNER_ADOPTION_2026-10-17.md"
ROOT = faner.REPO_ROOT


def json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def source_hashes():
    return {
        p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for directory in (ROOT / "src", ROOT / "paper/analysis/harness")
        for p in sorted(directory.rglob("*.py"))
    }


def load_start(args, hashes):
    parent_path = args.resume.with_name(args.resume.name.replace(".accepted.pkl", ".json"))
    parent = json.loads(parent_path.read_text())
    if not parent.get("sources_unchanged_during_run", False):
        raise ValueError("parent source identity was not preserved")
    parent_codec = parent.get("checkpoint_codec")
    if parent_codec not in (CODEC, PAPER_CHECKPOINT_CODEC):
        raise ValueError("unsupported checkpoint codec; legacy fixtures are separate diagnostics")
    for path, digest in parent.get("source_sha256_at_exit", {}).items():
        relative = path.replace("\\", "/").split("/src/dtdc_simulator/core2/")
        if len(relative) == 2 and hashes["src/dtdc_simulator/core2/" + relative[1]] != digest:
            raise ValueError("parent native physical/numerical source changed: " + path)
    checkpoint = pickle.loads(args.resume.read_bytes())
    if checkpoint.get("checkpoint_codec") != parent_codec:
        raise ValueError("checkpoint and parent codec differ")
    if parent_codec == CODEC:
        if checkpoint["source_sha256"] != hashes:
            raise ValueError("candidate restart must use exactly the same runtime sources")
        if checkpoint["controls"] != [
            args.study,
            args.maximum_dt,
            args.minimum_dt,
            args.core_time_tolerance,
            args.experimental_per_area,
            getattr(args, "per_area_qualification", False),
        ]:
            raise ValueError("restart must preserve study and time-step controls")
        return checkpoint, parent
    before = checkpoint["before"] if args.study == "Faner" else checkpoint["current"]
    if before.transport.config is not checkpoint["config"]:
        raise ValueError("accepted parent lost native state/config identity")
    if "surface" in checkpoint:
        surface = checkpoint["surface"]
        front_rate = checkpoint["front_rate"]
    else:
        birth = next(
            a["audit"] for a in parent["attempts"] if a.get("stage") == "birth" and a["accepted"]
        )
        surface = faner.shared.node.WaterFilmNodeAudit(**birth["surface_water_node"])
        front_rate = (1.0 - before.transport.geometry.front.z) / birth["native_birth_ledger"][
            "dt_s"
        ]
    config = before.transport.config
    params = ParticleModelParameters()
    if (
        config.wet.wet.oil_heat_capacity_law != params.oil_law
        or config.dry.pore.oil_heat_capacity_law != params.oil_law
        or config.dry.binary_diffusivity.temperature_law != params.binary_law
    ):
        raise ValueError("candidate requires properties applied before the original feed")
    return {
        "checkpoint_codec": CODEC,
        "current": before,
        "carry": checkpoint["carry"],
        "surface": surface,
        "boundary": checkpoint.get("boundary"),
        "front_rate": front_rate,
        "missing_face": None,
        "next_dt": args.maximum_dt,
        "source_sha256": hashes,
        "controls": [
            args.study,
            args.maximum_dt,
            args.minimum_dt,
            args.core_time_tolerance,
            args.experimental_per_area,
            getattr(args, "per_area_qualification", False),
        ],
        "positive_core_targets": [],
        "audit_origin": None,
        "audit_streams": [[], [], []],
    }, parent


class CheckpointReplaceRefused(RuntimeError):
    """The atomic checkpoint replace stayed locked through every bounded retry."""


CHECKPOINT_REPLACE_DELAYS_S = (0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)


def replace_with_retry(temporary, target, log=None):
    """Atomically replace target by temporary, retrying a transient file lock.

    Windows refuses os.replace while another process (for example an on-access
    scanner) holds the target without share-delete. The same bytes are retried
    with a bounded back-off (16.95 s in total); every refused attempt is
    appended to log; after the last delay a typed refusal is raised. The
    checkpoint content and the accepted sequence are unchanged.
    """
    for attempt, delay in enumerate((*CHECKPOINT_REPLACE_DELAYS_S, None)):
        try:
            temporary.replace(target)
            return attempt
        except PermissionError as error:
            if log is not None:
                log.append(
                    {
                        "target": str(target),
                        "attempt": attempt,
                        "error": str(error),
                        "next_delay_s": delay,
                    }
                )
            if delay is None:
                raise CheckpointReplaceRefused(
                    f"checkpoint replace stayed locked after {attempt + 1} attempts: {target}"
                ) from error
            time.sleep(delay)


# Kernel-envelope PRESSURE_SCHEDULE armed with a bath (the Class-A front-ledger repair of
# 2026-10-26, FRONT_LEDGER_RECORD; the diagnosis GT_PS2_UNIVERSAL_KERNEL_FIXUP_VERIFICATION_AND_
# PRESSURE_DRIFT_READING_2026-10-22.md section 3): the kernel commits every front state on its
# pressure-only target configuration, so the driver's independent audit reads each state's storage
# on that state's own configuration (strict-ward checked to differ from the run's start
# configuration in its pressure alone), and the solver controls' interface-temperature chart is
# re-derived for every pressure-target advance by its own two declared rules over the advance's
# pressure range (section 4).  With the knob absent, or without a bath, none of this is entered;
# armed at a held pressure, each state's own configuration is the run's start configuration.
FRONT_LEDGER_RECORD = (
    "docs/GT_PS2_UNIVERSAL_KERNEL_FRONT_LEDGER_OWN_CONFIGURATION_REPAIR_2026-10-26.md"
)
FRONT_LEDGER_REFUSAL_PREFIX = "front ledger own-configuration ward refused"
# The declared chart offset above the native two-liquid root, the two birth sites' value
# (experiment_core2_dt_birth.initialize_from_feed, run_particle_qualification_case.fresh_dt_birth).
PRESSURE_CHART_ROOT_OFFSET_K = 1e-6


class FrontLedgerConfigurationRefused(RuntimeError):
    """A committed state's own configuration differs from the run's start configuration in a
    field other than its pressure (the strict ward of the own-configuration reading)."""


def configuration_differences(expected, actual, path="config"):
    """Dotted paths of every field in which two configurations differ (floats by their bits)."""
    if expected is actual:
        return []
    if is_dataclass(expected) and not isinstance(expected, type):
        if type(actual) is not type(expected):
            return [path]
        out = []
        for item in fields(expected):
            out += configuration_differences(
                getattr(expected, item.name), getattr(actual, item.name), f"{path}.{item.name}"
            )
        return out
    if isinstance(expected, (tuple, list)):
        if type(actual) is not type(expected) or len(actual) != len(expected):
            return [path]
        out = []
        for index, (a, b) in enumerate(zip(expected, actual)):
            out += configuration_differences(a, b, f"{path}[{index}]")
        return out
    if isinstance(expected, float):
        if isinstance(actual, float) and struct.pack("<d", expected) == struct.pack("<d", actual):
            return []
        return [path]
    return [] if type(actual) is type(expected) and actual == expected else [path]


def own_ledger_configuration(run_start_config, state):
    """The configuration the independent audit reads a state on: the state's own.

    Strict ward: the state's own configuration must be the run's start configuration with its
    pressure moved by the kernel's own pressure-only target (cut.pressure_only_target_config);
    any other differing field is refused, typed.
    """
    own = state.transport.config
    if own is run_start_config:
        return own
    expected = engine.core.cut.pressure_only_target_config(run_start_config, own.dry.pressure_pa)
    differences = configuration_differences(expected, own)
    if differences:
        raise FrontLedgerConfigurationRefused(
            f"{FRONT_LEDGER_REFUSAL_PREFIX}: the state at {state.transport.time_s!r} s differs "
            f"from the run's start configuration beyond its pressure in {differences[:6]}"
        )
    return own


def pressure_target_chart(controls, before_config, target_config):
    """The solver controls' interface-temperature chart for one pressure-target advance.

    Re-derived by its own two declared rules over the advance's pressure range, so that it holds
    on both configurations the integrator validates it on (the before and the target state): the
    lower end the native two-liquid root + 1e-6 K at the higher of the two pressures (the birth
    sites' rule), the upper end one ulp below the pure n-hexane boiling point at the lower
    (positive_water_continuation.continued_controls).  Nothing else in the controls changes.
    """
    import positive_water_continuation as continuation

    low, high = sorted((before_config, target_config), key=lambda c: c.dry.pressure_pa)
    upper = continuation.continued_controls(controls, low).interface_temperature_bounds_k[1]
    lower = envelope.two_liquid_root_k(high.dry.pressure_pa) + PRESSURE_CHART_ROOT_OFFSET_K
    return replace(controls, interface_temperature_bounds_k=(lower, upper))


def run(args):
    hashes = source_hashes()
    checkpoint, parent = load_start(args, hashes)
    current, carry, surface = (checkpoint[k] for k in ("current", "carry", "surface"))
    config = current.transport.config
    bath = json.loads(args.bath.read_text()) if args.bath else None
    record = {
        "checkpoint_codec": CODEC,
        "candidate_spatial_scheme": True,
        "physically_qualifying": False,
        "convergence_established": False,
        "parent_checkpoint": str(args.resume),
        "parent_sha256": hashlib.sha256(args.resume.read_bytes()).hexdigest(),
        "parent_result_sha256": hashlib.sha256(
            args.resume.with_name(args.resume.name.replace(".accepted.pkl", ".json")).read_bytes()
        ).hexdigest(),
        "source_sha256_at_entry": hashes,
        "property_parameters": asdict(ParticleModelParameters()),
        "inputs": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "original_physical_and_conservation_gates_retained": True,
        "accepted_steps": 0,
        "mesh_operations": 0,
        "refusals": 0,
        "whole_history_complete": False,
        "per_area_interface_conditioning": args.experimental_per_area,
        "experimental_per_area_conditioning": args.experimental_per_area
        and not getattr(args, "per_area_qualification", False),
        "owner_adoption_required": args.experimental_per_area
        and not getattr(args, "per_area_qualification", False),
        "conditioning_authority": (
            AREA_CHART_AUTHORITY if getattr(args, "per_area_qualification", False) else None
        ),
    }
    log_path = args.out.with_suffix(".history.jsonl")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    log = log_path.open("x", encoding="utf-8")

    def append(row):
        log.write(json.dumps(row, default=json_default) + "\n")
        log.flush()

    def save():
        record["source_sha256_at_exit"] = source_hashes()
        record["sources_unchanged_during_run"] = record["source_sha256_at_exit"] == hashes
        args.out.write_text(json.dumps(record, indent=2, default=json_default))
        checkpoint.update(current=current, carry=carry, surface=surface)
        target = args.out.with_suffix(".accepted.pkl")
        temp = target.with_suffix(".tmp")
        temp.write_bytes(pickle.dumps(checkpoint, protocol=5))
        replace_with_retry(temp, target, record.setdefault("checkpoint_replace_retries", []))

    # Declared option: search no tighter than the native acceptance gate. Off keeps
    # the former 1e-11 front search target; acceptance checks are unchanged.
    front_search_tolerance = (
        current.controls.nonlinear_residual_tolerance
        if getattr(args, "front_search_at_acceptance_gate", False)
        else 1e-11
    )
    if getattr(args, "front_search_at_acceptance_gate", False):
        record["declared_front_search_tolerance"] = front_search_tolerance
    # Kernel-envelope knobs (particle_kernel_envelope); all off adds nothing below.
    knobs = envelope.validate_environment()
    # A result-neutral knob (envelope.RESULT_NEUTRAL_ITEMS; memory only) arms no audit here.
    envelope_on = envelope.behaviour_armed(knobs)
    if knobs[envelope.FRONT_RELEASE_SOLVE_MEMORY]:
        record["kernel_envelope_release_solve_memory"] = envelope.release_solve_memory_row()
    # Kernel-envelope CORE_TARGET_SEED (valued; armed only): the positive core target's duration
    # seed from the accepted front history; absent, the target step keeps the driver's seed.
    seed_rule = envelope.core_target_seed()
    # The adopted film-covered surface law (2026-10-25): the surface face under the film carries
    # the contact on all its paths (the finite-shell scope at the declared conductance); each
    # accepted step is stamped, report only.
    record["surface_contact_law"] = envelope.surface_contact_law_declared_row()
    if seed_rule is not None:
        seed_points, seed_bootstrap = envelope.core_target_seed_points(
            checkpoint, current.transport.time_s, current.transport.geometry.front.z, args.resume
        )
        record["kernel_envelope_core_target_seed"] = {
            "declared_row": envelope.core_target_seed_declared_row(seed_rule),
            "bootstrap": seed_bootstrap,
            "attempts": [],
        }
    # Kernel-envelope BUDGETS_PER_UNKNOWN (valued; ruling packet 6, record
    # GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md; armed only): the front's proof-point
    # allowance per advance and the arrival's search budget per unknown (the partial, target and
    # event searches read the knob in experiment_core2_finite_wetting_b.advance_partial).
    budgets_rule = envelope.budgets_per_unknown()
    if budgets_rule is not None:
        record["kernel_envelope_budgets_per_unknown"] = {
            "declared_row": envelope.budgets_per_unknown_declared_row(budgets_rule),
            "advances": 0,
            "advances_with_proof_allowance_above_floor": 0,
            "maximum_proof_point_allowance": None,
            "maximum_unknowns": None,
        }
    film_free = bool(checkpoint.get("kernel_envelope_film_free", False))
    film_points = list(
        checkpoint.get("kernel_envelope_film_points")
        or [(current.transport.time_s, carry["external_water_loading_kg_m2"])]
    )
    if envelope_on:
        record["kernel_envelope"] = {
            "armed": [item for item, on in knobs.items() if on],
            "declared_rows": envelope.declared_rows(),
            "film_free_from_s": checkpoint.get("kernel_envelope_film_free_from_s"),
            "film_depletion_events": [],
            "film_depletion_event_refusals": [],
        }
        if knobs[envelope.STALLED_FRONT]:
            record["kernel_envelope"].update(stalled_steps=[], stall_refusals=[], front_resumed=[])
    # Kernel-envelope STALL_FILM_DEPLETION_EVENT (armed only, with STALLED_FRONT): the held
    # front's film-depletion event and the declined-advance materiality rows.
    stall_event_on = envelope_on and knobs[envelope.STALL_FILM_DEPLETION_EVENT]
    stall_event_tried_at = None
    if stall_event_on:
        record["kernel_envelope"].update(
            stalled_film_depletion_events=[],
            stalled_film_depletion_event_refusals=[],
            stall_materiality={
                "rows": 0,
                "rows_evaluated": 0,
                "cumulative_declined_hexane_kg_kg_dry": checkpoint.get(
                    "kernel_envelope_stall_declined_kg_kg_dry", 0.0
                ),
                "record": envelope.STALL_RECORD,
                "report_only": True,
            },
        )
    stalled_last = bool(checkpoint.get("kernel_envelope_stalled_last", False))
    stall_last_demand = checkpoint.get("kernel_envelope_stall_last_demand")
    start = time.perf_counter()
    try:
        with ExitStack() as scope:
            if args.study == "Faner":
                boundary = checkpoint["boundary"]
                scope.enter_context(
                    faner.declared_faner_film_multiplier(boundary.film_boundary.stress_multiplier)
                )
                scope.enter_context(faner.pure_bulk_domain_adapter())

                def boundary_factory(clock):
                    # Kernel-envelope FANER_MEASURED_VAPOUR_TRACE (armed only): the trace
                    # at the clock; off, the constant boundary itself.
                    return faner_trace.bath_boundary(boundary, clock)

            else:
                scope.enter_context(
                    engine.finite_shell_path(100.0, envelope.surface_contact_conductance())
                )
                radius = current.transport.geometry.master_grid.R

                def boundary_factory(clock):
                    return engine.face.bath_boundary(bath, radius, clock)

            scope.enter_context(engine.study.magnitude_storage_scale())
            scope.enter_context(engine.driver.path_bounds.monotone_water_path_enclosure([]))
            scope.enter_context(engine.driver.path_bounds.monotone_hexane_path_enclosure([]))
            begin_advance = scope.enter_context(
                engine.driver.physical_advance_proof_budget(
                    3000000,
                    120.0,
                    **(
                        {"memoize_identical_points": True}
                        if getattr(args, "memoize_proof_points", False)
                        else {}
                    ),
                )
            )
            # PRESSURE_SCHEDULE with a bath (FRONT_LEDGER_RECORD): the independent audit reads
            # each state on its own configuration and the solver chart follows each advance's
            # pressure range; otherwise the configuration is constant and `config` is read.
            own_configuration_reading = bool(
                envelope_on and knobs[envelope.PRESSURE_SCHEDULE] and bath is not None
            )

            def ledger_configuration(state):
                # the configuration the audit reads `state` on (its own under the knob)
                if not own_configuration_reading:
                    return config
                return own_ledger_configuration(config, state)

            def charted(state, options, duration):
                # the state an advance starts from, its interface-temperature chart re-derived
                # for the advance's pressure range when the advance carries a pressure target
                factory = options.get("target_config_factory")
                target = None if factory is None else factory(duration)
                if target is None:
                    return state
                rederived = pressure_target_chart(state.controls, state.transport.config, target)
                if rederived == state.controls:
                    return state
                chart = record["kernel_envelope"].setdefault(
                    "front_ledger_own_configuration",
                    {"record": FRONT_LEDGER_RECORD, "charts_rederived": 0},
                )
                chart["charts_rederived"] += 1
                return replace(state, controls=rederived)

            first = engine.primitive_row(current, config, carry, surface, "initial")
            if checkpoint["audit_origin"] is None:
                checkpoint["audit_origin"] = {
                    "time_s": first["time_s"],
                    "W_H_U": first["combined_W_H_U"],
                    "scales": [
                        first["combined_W_H_U"][0],
                        first["combined_W_H_U"][1],
                        carry["import_capacity_scale_j"],
                    ],
                }
            record["independent_audit_origin"] = checkpoint["audit_origin"]
            append(first)
            save()
            wall_budget = getattr(args, "wall_budget_s", None)
            for index in range(args.steps):
                transport = current.transport
                if transport.time_s >= args.until:
                    record["requested_horizon_reached"] = True
                    break
                if wall_budget is not None and time.perf_counter() - start > wall_budget:
                    # declared run control (a measurement campaign's wall allowance)
                    record["stop"] = "WALL_BUDGET"
                    break
                cut_index = transport.geometry.cut_cell_index
                if (
                    cut_index == 0
                    and transport.geometry.front.z <= 1e-12
                    and len(checkpoint["positive_core_targets"]) >= 4
                ):
                    try:
                        time_control = estimate_core_time(
                            current, checkpoint["positive_core_targets"], args.core_time_tolerance
                        )
                    except CoreTimeNotResolved:
                        # More strictly accepted targets are needed. This
                        # does not suppress a failed solve or admit a state.
                        pass
                    else:
                        record["native_endpoint_inputs_ready"] = True
                        record["conditional_endpoint_time_control"] = asdict(time_control)
                        break
                if budgets_rule is None:
                    path_work = begin_advance()
                else:
                    unknowns = envelope.front_search_unknowns(transport)
                    path_work = begin_advance(
                        point_budget_form=lambda present, n=unknowns: (
                            envelope.front_proof_point_allowance(present, n)
                        )
                    )
                    stamp = record["kernel_envelope_budgets_per_unknown"]
                    stamp["advances"] += 1
                    if path_work["point_budget"] > path_work["present_point_budget"]:
                        stamp["advances_with_proof_allowance_above_floor"] += 1
                    stamp["maximum_proof_point_allowance"] = max(
                        stamp["maximum_proof_point_allowance"] or 0, path_work["point_budget"]
                    )
                    stamp["maximum_unknowns"] = max(stamp["maximum_unknowns"] or 0, unknowns)
                remaining = args.until - transport.time_s
                dt = min(checkpoint["next_dt"], remaining)

                def envelope_options(transport=transport):
                    if not envelope_on:
                        return {}
                    options = {"envelope_audit": {}}
                    if knobs[envelope.FILM_NODE_CERTIFICATE_SCALED]:
                        options["scaled_condition_cap"] = envelope.FILM_NODE_SCALED_CONDITION_CAP
                    if knobs[envelope.PRESSURE_SCHEDULE] and bath is not None:
                        origin, base = transport.time_s, transport.config

                        def target_config_factory(duration):
                            pressure = envelope.bath_pressure_pa(bath, origin + duration)
                            if pressure == base.dry.pressure_pa:
                                return None
                            return engine.core.cut.pressure_only_target_config(base, pressure)

                        options["target_config_factory"] = target_config_factory
                    return options

                event_duration = (
                    envelope.linear_film_zero_duration(film_points)
                    if knobs[envelope.FRONT_FILM_DEPLETION_EVENT] and not film_free
                    else None
                )
                lower_z = (
                    transport.geometry.master_grid.faces[cut_index]
                    / transport.geometry.master_grid.R
                ) ** 3
                rate = checkpoint["front_rate"]
                estimated = (transport.geometry.front.z - lower_z) / rate if rate > 0 else math.inf
                event_attempt = estimated <= dt

                def held_placeholder(duration, transport=transport, lower_z=lower_z, rate=rate):
                    # the seed's encodable front (the stalled construction never uses it)
                    z = transport.geometry.front.z - rate * duration
                    if not lower_z < z < transport.geometry.front.z:
                        z = 0.5 * (lower_z + transport.geometry.front.z)
                    return z

                def held_front_one(duration, gas_only_flag):
                    # STALL_FILM_DEPLETION_EVENT's materiality: the counterfactual advance
                    # demand from this step's start state (scratch overlay, never adopted)
                    probe_options = envelope_options()
                    return materiality.counterfactual_demand(
                        engine,
                        charted(current, probe_options, duration),
                        carry,
                        surface,
                        bath,
                        duration,
                        initial_front_z=held_placeholder(duration),
                        boundary_factory=boundary_factory,
                        gas_only_surface=gas_only_flag,
                        options=probe_options,
                        search_tolerance=front_search_tolerance,
                    )

                def held_front_demand(duration, gas_only_flag, start=None):
                    if start is None:
                        start = held_front_one(materiality.DEMAND_START_PROBE_S, gas_only_flag)
                    return start, held_front_one(duration, gas_only_flag)

                def held_front_row(start, end, elapsed_s, cumulative, hexane=None):
                    if hexane is None:
                        hexane = engine.primitive_row(
                            current, ledger_configuration(current), carry, surface, "materiality"
                        )["hexane_kg_kg_dry"]
                    return materiality.materiality_row(
                        start,
                        end,
                        elapsed_s,
                        materiality.particle_dry_mass_kg(current.transport),
                        hexane,
                        cumulative,
                    )

                pending_windows = []

                def settle_pending_windows():
                    # the refused held-front events' windows, from this advance's start state
                    while pending_windows:
                        entry, located, start = pending_windows.pop(0)
                        if start is None:
                            start = held_front_one(materiality.DEMAND_START_PROBE_S, False)
                        rungs = [
                            held_front_one(fraction * located, False)
                            for fraction in materiality.WINDOW_LADDER_FRACTIONS
                        ]
                        entry["event_window_materiality"] = materiality.window_row(
                            start,
                            rungs,
                            located,
                            materiality.particle_dry_mass_kg(current.transport),
                            engine.primitive_row(
                                current,
                                ledger_configuration(current),
                                carry,
                                surface,
                                "materiality",
                            )["hexane_kg_kg_dry"],
                        )

                engine.wetting.surface_contact_reset()
                while dt >= args.minimum_dt:
                    committed = current
                    positive_target = cut_index == 0 and event_attempt
                    film_event = (
                        event_duration is not None
                        and event_duration <= dt
                        and not (event_attempt and estimated <= event_duration)
                    )
                    kind = (
                        "film_depletion_event"
                        if film_event
                        else (
                            "positive_target"
                            if positive_target
                            else "arrival" if event_attempt else "partial"
                        )
                    )
                    gas_only = args.study == "Faner" or film_free
                    surface_scope = (
                        engine.film_free_surface_scope() if film_free else contextlib.nullcontext()
                    )
                    envelope_audit = None
                    try:
                        if film_event:
                            options = envelope_options()
                            envelope_audit = options.get("envelope_audit")
                            prediction = transport.geometry.front.z - rate * event_duration
                            if not lower_z < prediction < transport.geometry.front.z:
                                prediction = 0.5 * (lower_z + transport.geometry.front.z)
                            step = engine.advance_partial(
                                charted(current, options, event_duration),
                                carry,
                                surface,
                                bath,
                                event_duration,
                                initial_front_z=prediction,
                                search_tolerance=front_search_tolerance,
                                boundary_factory=boundary_factory,
                                gas_only_surface=False,
                                film_depletion_event=True,
                                **options,
                            )
                            new_surface = step.surface_film_audit
                            if new_surface.liquid_water_kg_m2 != 0.0:
                                raise RuntimeError("film depletion event left a nonzero film")
                            transfer = engine.active_shell_transfer_audit()
                            envelope_audit["d9_transfer_at_event"] = (
                                None if transfer is None else dict(transfer)
                            )
                        elif positive_target:
                            area_audit = {} if args.experimental_per_area else None
                            target = 0.5 * transport.geometry.front.z
                            options = envelope_options()
                            envelope_audit = options.get("envelope_audit")
                            duration_seed = 0.5 * estimated
                            if seed_rule is not None:
                                duration_seed, seed_stamp = envelope.core_target_duration_seed(
                                    seed_rule,
                                    seed_points,
                                    transport.time_s,
                                    transport.geometry.front.z,
                                    duration_seed,
                                )
                                record["kernel_envelope_core_target_seed"]["attempts"].append(
                                    seed_stamp
                                )
                            with surface_scope:
                                step = engine.advance_partial(
                                    charted(current, options, duration_seed),
                                    carry,
                                    surface,
                                    bath,
                                    duration_seed,
                                    target=target,
                                    search_tolerance=front_search_tolerance,
                                    boundary_factory=boundary_factory,
                                    gas_only_surface=gas_only,
                                    experimental_interface_area_chart_audit=area_audit,
                                    **options,
                                )
                            new_surface = step.surface_film_audit
                            if area_audit is not None and getattr(
                                args, "per_area_qualification", False
                            ):
                                area_audit.update(
                                    experimental_counterfactual=False,
                                    owner_adoption_required=False,
                                    qualification_authority=AREA_CHART_AUTHORITY,
                                    scope="owner-approved qualification run; not qualified by this step",
                                )
                        elif event_attempt:
                            with surface_scope:
                                step, new_surface, event_audit = engine.face.advance_arrival(
                                    current,
                                    carry,
                                    surface,
                                    bath,
                                    time_guess=estimated,
                                    time_upper=dt,
                                    budget=2400,
                                    boundary_factory=boundary_factory,
                                    gas_only_surface=gas_only,
                                    search_tolerance=front_search_tolerance,
                                    **(
                                        {"budget_per_unknown": True}
                                        if budgets_rule is not None
                                        else {}
                                    ),
                                )
                        else:
                            prediction = transport.geometry.front.z - rate * dt
                            if not lower_z < prediction < transport.geometry.front.z:
                                prediction = 0.5 * (lower_z + transport.geometry.front.z)
                            options = envelope_options()
                            envelope_audit = options.get("envelope_audit")
                            with surface_scope:
                                try:
                                    step = engine.advance_partial(
                                        charted(current, options, dt),
                                        carry,
                                        surface,
                                        bath,
                                        dt,
                                        initial_front_z=prediction,
                                        search_tolerance=front_search_tolerance,
                                        boundary_factory=boundary_factory,
                                        gas_only_surface=gas_only,
                                        **options,
                                    )
                                except Exception as moving_error:
                                    if isinstance(
                                        moving_error,
                                        (
                                            StepBudgetExceeded,
                                            pd_boundary.DewLineCondensateOnFilmFreeSurface,
                                        ),
                                    ) or not (envelope_on and knobs[envelope.STALLED_FRONT]):
                                        raise
                                    # Kernel-envelope stalled front: no receding root at
                                    # this duration; the front held, same duration, same
                                    # gates. The moving refusal is logged first.
                                    append(
                                        {
                                            "kind": "refusal",
                                            "time_s": transport.time_s,
                                            "dt_s": dt,
                                            "attempted": "partial",
                                            "error": str(moving_error),
                                            "exact_rollback": current is committed,
                                            "followed_by": "stalled_front",
                                        }
                                    )
                                    record["refusals"] += 1
                                    kind = "stalled_front"
                                    if (
                                        stall_event_on
                                        and not film_free
                                        and carry["external_water_loading_kg_m2"] > 0.0
                                        and stall_event_tried_at != transport.time_s
                                    ):
                                        # STALL_FILM_DEPLETION_EVENT: the held front's film
                                        # depletion located as an event first, once per start
                                        # state, seeded at the films' linear zero
                                        stall_event_tried_at = transport.time_s
                                        event_seed = envelope.linear_film_zero_duration(film_points)
                                        event_options = envelope_options()
                                        try:
                                            if event_seed is None:
                                                raise RuntimeError(
                                                    "no decreasing positive film pair to seed "
                                                    "the held front's depletion event"
                                                )
                                            step = engine.advance_partial(
                                                charted(current, event_options, event_seed),
                                                carry,
                                                surface,
                                                bath,
                                                event_seed,
                                                initial_front_z=held_placeholder(event_seed),
                                                search_tolerance=front_search_tolerance,
                                                boundary_factory=boundary_factory,
                                                gas_only_surface=False,
                                                stalled_front=True,
                                                film_depletion_event=True,
                                                **event_options,
                                            )
                                            if step.surface_film_audit.liquid_water_kg_m2 != 0.0:
                                                raise RuntimeError(
                                                    "stalled film depletion event left a "
                                                    "nonzero film"
                                                )
                                        except (
                                            StepBudgetExceeded,
                                            pd_boundary.DewLineCondensateOnFilmFreeSurface,
                                        ):
                                            raise
                                        except Exception as event_error:
                                            stamp = getattr(event_error, "stall_stamp", None)
                                            append(
                                                {
                                                    "kind": "refusal",
                                                    "time_s": transport.time_s,
                                                    "dt_s": event_seed,
                                                    "attempted": "stalled_film_depletion_event",
                                                    "error": str(event_error),
                                                    "exact_rollback": current is committed,
                                                    "followed_by": "stalled_front",
                                                    "stall_stamp": stamp,
                                                }
                                            )
                                            record["refusals"] += 1
                                            entry = {
                                                "time_s": transport.time_s,
                                                "duration_seed_s": event_seed,
                                                "error": str(event_error),
                                                "stall_stamp": stamp,
                                                "event_window_materiality": None,
                                            }
                                            record["kernel_envelope"][
                                                "stalled_film_depletion_event_refusals"
                                            ].append(entry)
                                            located = (stamp or {}).get("event_duration_s")
                                            if located:
                                                # the declined advance over the window the
                                                # refused event spans, measured once this
                                                # advance is over (outside its allowance)
                                                pending_windows.append(
                                                    (
                                                        entry,
                                                        located,
                                                        stall_last_demand if stalled_last else None,
                                                    )
                                                )
                                        else:
                                            kind = "stalled_film_depletion_event"
                                            options = event_options
                                            envelope_audit = options.get("envelope_audit")
                                    if kind == "stalled_front":
                                        options = envelope_options()
                                        envelope_audit = options.get("envelope_audit")
                                        try:
                                            step = engine.advance_partial(
                                                charted(current, options, dt),
                                                carry,
                                                surface,
                                                bath,
                                                dt,
                                                initial_front_z=prediction,
                                                search_tolerance=front_search_tolerance,
                                                boundary_factory=boundary_factory,
                                                gas_only_surface=gas_only,
                                                stalled_front=True,
                                                **options,
                                            )
                                        except Exception as stall_error:
                                            stamp = getattr(stall_error, "stall_stamp", None)
                                            record["kernel_envelope"]["stall_refusals"].append(
                                                {
                                                    "time_s": transport.time_s,
                                                    "dt_s": dt,
                                                    "error": str(stall_error),
                                                    "consistency_scaled_residual": (
                                                        None
                                                        if stamp is None
                                                        else stamp["consistency_scaled_residual"]
                                                    ),
                                                }
                                            )
                                            raise
                                    envelope_audit["stalled_front"]["reason"] = str(moving_error)
                            new_surface = step.surface_film_audit
                        combined = engine.driver.combined_ledger(
                            current,
                            config,
                            carry,
                            step,
                            new_surface,
                            **(
                                {
                                    # an exact-face arrival has no pressure target: it
                                    # keeps the start pressure
                                    "film_pressure_pa": getattr(
                                        step.after, "transport", transport
                                    ).config.dry.pressure_pa
                                }
                                if envelope_on and knobs[envelope.PRESSURE_SCHEDULE]
                                else {}
                            ),
                        )
                        if (
                            combined["maximum_normalized_residual"]
                            > current.controls.ledger_tolerance
                        ):
                            raise RuntimeError(
                                "original particle-plus-film conservation gate refused"
                            )
                        new_carry = {
                            **carry,
                            **combined,
                            "external_water_loading_kg_m2": new_surface.liquid_water_kg_m2,
                        }
                        new_current = step.after
                        mesh_audit = None
                        if event_attempt and not positive_target and not film_event:
                            if checkpoint["missing_face"] is None:
                                new_current, mesh_audit = coalesce_arrival_face(new_current)
                            else:
                                new_current, mesh_audit = slide_arrival_face(
                                    new_current, checkpoint["missing_face"]
                                )
                        row = engine.primitive_row(
                            new_current,
                            ledger_configuration(new_current),
                            new_carry,
                            new_surface,
                            kind,
                        )
                        elapsed = row["time_s"] - transport.time_s
                        if elapsed <= 0:
                            raise RuntimeError("physical advance must have positive elapsed time")
                        row["kernel_envelope_finite_contact"] = (
                            engine.wetting.surface_contact_front_stamp(
                                carry["external_water_loading_kg_m2"],
                                elapsed,
                                new_surface,
                                engine.active_shell_transfer_audit(),
                            )
                        )
                        area = transport.geometry.master_grid.areas[-1]
                        streams = [s.copy() for s in checkpoint["audit_streams"]]
                        for stream, flux in zip(
                            streams,
                            (
                                new_surface.external_water_flux_kg_m2_s,
                                new_surface.external_hexane_flux_kg_m2_s,
                                new_surface.external_energy_flux_w_m2,
                            ),
                        ):
                            stream.append(elapsed * area * flux)
                        origin = checkpoint["audit_origin"]
                        residuals = [
                            math.fsum((a, -b, *s))
                            for a, b, s in zip(row["combined_W_H_U"], origin["W_H_U"], streams)
                        ]
                        normalized = [abs(r) / s for r, s in zip(residuals, origin["scales"])]
                        if max(normalized) > current.controls.ledger_tolerance:
                            raise RuntimeError(
                                "independent primitive/rate cumulative conservation refused: "
                                + str(normalized)
                            )
                        new_z = new_current.transport.geometry.front.z
                        if positive_target and not film_event:
                            next_targets = [
                                *checkpoint["positive_core_targets"],
                                PositiveCoreTargetEvidence(current, new_current, step.ledger.dt_s),
                            ][-4:]
                        else:
                            next_targets = []
                        stalled = kind == "stalled_front" or kind == "stalled_film_depletion_event"
                        next_checkpoint = {
                            **checkpoint,
                            "front_rate": (
                                # a stalled step keeps the last moving rate (seeding only)
                                checkpoint["front_rate"]
                                if stalled
                                else (transport.geometry.front.z - new_z) / elapsed
                            ),
                            "next_dt": min(args.maximum_dt, elapsed * 1.5),
                            "audit_streams": streams,
                            "positive_core_targets": next_targets,
                            "missing_face": (
                                mesh_audit.removed_face_radius_m
                                if mesh_audit
                                else checkpoint["missing_face"]
                            ),
                        }
                        row.update(
                            ledger=asdict(step.ledger),
                            combined_conservation=combined,
                            independent_normalized_residuals=normalized,
                            path_work=dict(path_work),
                            mesh_operation=asdict(mesh_audit) if mesh_audit else None,
                            experimental_interface_area_chart=(
                                area_audit if positive_target else None
                            ),
                        )
                        if envelope_on and knobs[envelope.STALLED_FRONT]:
                            ledger = record["kernel_envelope"]
                            if kind == "stalled_film_depletion_event":
                                ledger["stalled_film_depletion_events"].append(
                                    {
                                        "time_s": row["time_s"],
                                        "duration_s": elapsed,
                                        "duration_seed_s": event_seed,
                                        "front_z": new_z,
                                        "film_after_kg_m2": new_surface.liquid_water_kg_m2,
                                        "consistency_scaled_residual": envelope_audit[
                                            "stalled_front"
                                        ]["consistency_scaled_residual"],
                                        "transfer_at_event": (
                                            None
                                            if engine.active_shell_transfer_audit() is None
                                            else dict(engine.active_shell_transfer_audit())
                                        ),
                                    }
                                )
                            elif stalled:
                                ledger["stalled_steps"].append(
                                    {
                                        "time_s": row["time_s"],
                                        "dt_s": elapsed,
                                        "front_z": new_z,
                                        "consistency_scaled_residual": envelope_audit[
                                            "stalled_front"
                                        ]["consistency_scaled_residual"],
                                        "reason": envelope_audit["stalled_front"]["reason"],
                                    }
                                )
                            elif stalled_last:
                                ledger["front_resumed"].append(
                                    {"time_s": row["time_s"], "kind": kind, "front_z": new_z}
                                )
                            stalled_last = stalled
                            next_checkpoint["kernel_envelope_stalled_last"] = stalled
                        if stall_event_on:
                            # the declined-advance materiality row of a held-front step and
                            # the windows of refused held-front events, measured now that the
                            # advance is accepted (scratch overlay; nothing is adopted)
                            try:
                                settle_pending_windows()
                                summary = record["kernel_envelope"]["stall_materiality"]
                                if stalled:
                                    start_rate, end_rate = held_front_demand(
                                        elapsed,
                                        gas_only and kind == "stalled_front",
                                        start=stall_last_demand if stall_last_demand else None,
                                    )
                                    stall_row = held_front_row(
                                        start_rate,
                                        end_rate,
                                        elapsed,
                                        summary["cumulative_declined_hexane_kg_kg_dry"],
                                        hexane=row["hexane_kg_kg_dry"],
                                    )
                                    summary["rows"] += 1
                                    if stall_row["evaluated"]:
                                        summary["rows_evaluated"] += 1
                                        summary["cumulative_declined_hexane_kg_kg_dry"] = stall_row[
                                            "cumulative_declined_hexane_kg_kg_dry"
                                        ]
                                        summary["ratio_to_particle_hexane"] = stall_row[
                                            "ratio_to_particle_hexane"
                                        ]
                                    envelope_audit["stall_materiality"] = stall_row
                                    stall_last_demand = end_rate
                                else:
                                    stall_last_demand = None
                            except Exception as measurement_error:  # never changes the run
                                record["kernel_envelope"].setdefault(
                                    "stall_materiality_errors", []
                                ).append({"time_s": row["time_s"], "error": str(measurement_error)})
                                stall_last_demand = None
                            next_checkpoint["kernel_envelope_stall_last_demand"] = stall_last_demand
                            next_checkpoint["kernel_envelope_stall_declined_kg_kg_dry"] = record[
                                "kernel_envelope"
                            ]["stall_materiality"]["cumulative_declined_hexane_kg_kg_dry"]
                        if envelope_on and (
                            knobs[envelope.CONDUCTION_BOUNDARY]
                            or knobs[envelope.GAS_DEW_LINE_PROJECTION]
                        ):
                            stamp = pd_boundary.step_stamp(
                                boundary_factory(transport.time_s + elapsed),
                                new_surface.temperature_k,
                                elapsed,
                            )
                            if stamp is not None:
                                envelope_audit = {**(envelope_audit or {}), "pd_boundary": stamp}
                                record["kernel_envelope"].setdefault(
                                    "pd_boundary_summary", pd_boundary.empty_summary()
                                )
                                pd_boundary.add_to_summary(
                                    record["kernel_envelope"]["pd_boundary_summary"],
                                    stamp,
                                    row["time_s"],
                                )
                        if envelope_on:
                            row["kernel_envelope"] = envelope_audit
                            if film_event:
                                film_free = True
                                next_checkpoint["kernel_envelope_film_free"] = True
                                next_checkpoint["kernel_envelope_film_free_from_s"] = row["time_s"]
                                record["kernel_envelope"]["film_free_from_s"] = row["time_s"]
                                record["kernel_envelope"]["film_depletion_events"].append(
                                    {
                                        "time_s": row["time_s"],
                                        "duration_s": elapsed,
                                        "predicted_duration_s": event_duration,
                                        "front_z": new_z,
                                        "audit": envelope_audit,
                                    }
                                )
                            elif kind == "stalled_film_depletion_event":
                                # the held front's film is gone: film-free from here, the
                                # moving (receding) chart tried first on every later step
                                film_free = True
                                next_checkpoint["kernel_envelope_film_free"] = True
                                next_checkpoint["kernel_envelope_film_free_from_s"] = row["time_s"]
                                record["kernel_envelope"]["film_free_from_s"] = row["time_s"]
                        film_points.append((row["time_s"], new_surface.liquid_water_kg_m2))
                        if envelope_on and knobs[envelope.FRONT_FILM_DEPLETION_EVENT]:
                            next_checkpoint["kernel_envelope_film_points"] = film_points[-3:]
                        if seed_rule is not None and not stalled:
                            # CORE_TARGET_SEED: the accepted moving front point, kept for restarts
                            seed_points = [*seed_points, (row["time_s"], new_z)][-3:]
                            next_checkpoint[envelope.CORE_TARGET_SEED_POINTS_KEY] = seed_points
                            if positive_target and not film_event:
                                row["kernel_envelope_core_target_seed"] = seed_stamp
                        current, carry, surface = new_current, new_carry, new_surface
                        checkpoint = next_checkpoint
                        if mesh_audit is not None:
                            record["mesh_operations"] += 1
                        append(row)
                        record["accepted_steps"] += 1
                        record["last_accepted"] = {
                            k: v
                            for k, v in row.items()
                            if k
                            not in (
                                "after_transport",
                                "surface",
                                "ledger",
                                "combined_conservation",
                                "path_work",
                            )
                        }
                        save()
                        if event_attempt or record["accepted_steps"] % 25 == 0:
                            print(
                                json.dumps(
                                    {
                                        "steps": record["accepted_steps"],
                                        "mesh_operations": record["mesh_operations"],
                                        "time_s": row["time_s"],
                                        "front_z": new_z,
                                    }
                                ),
                                flush=True,
                            )
                        break
                    except Exception as error:
                        append(
                            {
                                "kind": "refusal",
                                "time_s": transport.time_s,
                                "dt_s": dt,
                                "attempted": kind,
                                "error": str(error),
                                "exact_rollback": current is committed,
                                "optimizer_debug": getattr(error, "optimizer_debug", None),
                            }
                        )
                        record["refusals"] += 1
                        if current is not committed:
                            raise RuntimeError(
                                "post-commit evidence persistence failed; no retry permitted"
                            ) from error
                        if isinstance(error, StepBudgetExceeded):
                            raise
                        if isinstance(error, pd_boundary.DewLineCondensateOnFilmFreeSurface):
                            raise
                        if film_event:
                            # A refused depletion event falls back to the former order.
                            record["kernel_envelope"]["film_depletion_event_refusals"].append(
                                {"time_s": transport.time_s, "error": str(error)}
                            )
                            event_duration = None
                        elif event_attempt:
                            event_attempt = False
                        else:
                            dt *= 0.5
                else:
                    if stall_event_on:
                        try:
                            settle_pending_windows()
                        except Exception as measurement_error:  # never changes the stop
                            record["kernel_envelope"].setdefault(
                                "stall_materiality_errors", []
                            ).append({"time_s": transport.time_s, "error": str(measurement_error)})
                    raise RuntimeError("declared minimum interval reached with exact rollback")
    except Exception as error:
        record.update(error=str(error), traceback=traceback.format_exc())
    finally:
        record["wall_seconds"] = time.perf_counter() - start
        save()
        log.close()
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--study", choices=("Faner", "DT"), required=True)
    parser.add_argument("--bath", type=Path)
    parser.add_argument("--maximum-dt", type=float, default=0.2)
    parser.add_argument("--minimum-dt", type=float, default=1e-5)
    parser.add_argument("--until", type=float, required=True)
    parser.add_argument("--steps", type=int, default=10000)
    parser.add_argument("--core-time-tolerance", type=float, default=0.002)
    parser.add_argument("--experimental-per-area", action="store_true")
    parser.add_argument("--per-area-qualification", action="store_true")
    args = parser.parse_args()
    args.experimental_per_area = args.experimental_per_area or args.per_area_qualification
    if args.out.exists() or args.out.with_suffix(".accepted.pkl").exists():
        parser.error("refusing to overwrite evidence")
    if not 0 < args.minimum_dt <= args.maximum_dt or args.steps <= 0:
        parser.error("invalid declared time controls")
    if not 0 < args.core_time_tolerance <= 0.002:
        parser.error("core-time tolerance may preserve or tighten the existing 2-ms ceiling")
    args.resume, args.out = args.resume.resolve(), args.out.resolve()
    if args.bath:
        args.bath = args.bath.resolve()
    result = run(args)
    print(
        json.dumps(
            {
                k: result.get(k)
                for k in (
                    "accepted_steps",
                    "mesh_operations",
                    "error",
                    "stop",
                    "requested_horizon_reached",
                )
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
