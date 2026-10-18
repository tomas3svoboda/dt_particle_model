"""Whole-front engineering experiment using conservative finite-volume faces.

Native phase, ALE/RH, component/energy and conditioning gates own every
physical advance. Mesh operations advance no time. Accepted/rejected trials,
primitive fields, old references and the controller persist across restarts.
This candidate needs full-history mesh/time refinement before qualification.
"""

import argparse
from contextlib import ExitStack
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import pickle
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
                    return boundary

            else:
                scope.enter_context(engine.finite_shell_path(100.0, 0.32))
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
            for index in range(args.steps):
                transport = current.transport
                if transport.time_s >= args.until:
                    record["requested_horizon_reached"] = True
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
                path_work = begin_advance()
                remaining = args.until - transport.time_s
                dt = min(checkpoint["next_dt"], remaining)
                lower_z = (
                    transport.geometry.master_grid.faces[cut_index]
                    / transport.geometry.master_grid.R
                ) ** 3
                rate = checkpoint["front_rate"]
                estimated = (transport.geometry.front.z - lower_z) / rate if rate > 0 else math.inf
                event_attempt = estimated <= dt
                while dt >= args.minimum_dt:
                    committed = current
                    positive_target = cut_index == 0 and event_attempt
                    kind = (
                        "positive_target"
                        if positive_target
                        else "arrival" if event_attempt else "partial"
                    )
                    try:
                        if positive_target:
                            area_audit = {} if args.experimental_per_area else None
                            target = 0.5 * transport.geometry.front.z
                            step = engine.advance_partial(
                                current,
                                carry,
                                surface,
                                bath,
                                0.5 * estimated,
                                target=target,
                                search_tolerance=front_search_tolerance,
                                boundary_factory=boundary_factory,
                                gas_only_surface=args.study == "Faner",
                                experimental_interface_area_chart_audit=area_audit,
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
                            step, new_surface, event_audit = engine.face.advance_arrival(
                                current,
                                carry,
                                surface,
                                bath,
                                time_guess=estimated,
                                time_upper=dt,
                                budget=2400,
                                boundary_factory=boundary_factory,
                                gas_only_surface=args.study == "Faner",
                                search_tolerance=front_search_tolerance,
                            )
                        else:
                            prediction = transport.geometry.front.z - rate * dt
                            if not lower_z < prediction < transport.geometry.front.z:
                                prediction = 0.5 * (lower_z + transport.geometry.front.z)
                            step = engine.advance_partial(
                                current,
                                carry,
                                surface,
                                bath,
                                dt,
                                initial_front_z=prediction,
                                search_tolerance=front_search_tolerance,
                                boundary_factory=boundary_factory,
                                gas_only_surface=args.study == "Faner",
                            )
                            new_surface = step.surface_film_audit
                        combined = engine.driver.combined_ledger(
                            current, config, carry, step, new_surface
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
                        if event_attempt and not positive_target:
                            if checkpoint["missing_face"] is None:
                                new_current, mesh_audit = coalesce_arrival_face(new_current)
                            else:
                                new_current, mesh_audit = slide_arrival_face(
                                    new_current, checkpoint["missing_face"]
                                )
                        row = engine.primitive_row(
                            new_current, config, new_carry, new_surface, kind
                        )
                        elapsed = row["time_s"] - transport.time_s
                        if elapsed <= 0:
                            raise RuntimeError("physical advance must have positive elapsed time")
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
                        if positive_target:
                            next_targets = [
                                *checkpoint["positive_core_targets"],
                                PositiveCoreTargetEvidence(current, new_current, step.ledger.dt_s),
                            ][-4:]
                        else:
                            next_targets = []
                        next_checkpoint = {
                            **checkpoint,
                            "front_rate": (transport.geometry.front.z - new_z) / elapsed,
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
                        if event_attempt:
                            event_attempt = False
                        else:
                            dt *= 0.5
                else:
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
