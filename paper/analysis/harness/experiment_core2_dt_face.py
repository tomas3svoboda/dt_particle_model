"""Native exact-face arrival with the carried conservative DT water node.

Trusted local accepted checkpoints only. No epsilon shell or state reset.
"""

from dataclasses import asdict, replace
from contextlib import contextmanager
from pathlib import Path
import argparse
import hashlib
import json
import math
import pickle
import os
import traceback
import numpy as np
import experiment_core2_dt_birth as driver
from dtdc_simulator.core2.particle import cut_face_event as face
from dtdc_simulator.core2.particle import cut_face_event_integrator as fi
from dtdc_simulator.core2.particle import cut_face_departure_integrator as di
from dtdc_simulator.core2.particle import cut_face_tangent as tangent

core = driver.core
study = driver.study


def bath_boundary(bath, radius, clock):
    temperature, y = core.interpolate_bath(bath["waypoints"], clock)
    native = core.reduced_film_boundary(
        temperature_k=temperature,
        y_hexane=y,
        pair=core.film_pair_at_radius(radius),
        forcing_timescale_s=0.075,
        label=f"same prescribed DT bath at exact event t={clock}",
    )
    return study.FiniteFilmStudyBoundary(
        temperature,
        y,
        "native exact-face water-node event",
        film_boundary=native,
        accepted_beta_limit=3.0,
        actual_surface_path=True,
        numerical_cap_seam=True,
    )


def advance_arrival(
    before,
    carry,
    surface,
    bath,
    *,
    time_guess=0.14,
    time_upper=2.0,
    budget=1800,
    boundary_factory=None,
    gas_only_surface=False,
    difference_mode="forward",
):
    layout = face.layout_for_arrival(before.transport)
    controls = fi.FaceEventControls((0.0, time_upper))
    chart = fi._coordinate_chart(before, controls)
    candidate = face.FaceArrivalUnknowns(
        before.transport.wet_temperatures_k[:-1],
        before.transport.wet_retained_water_loadings[:-1],
        before.transport.dry_temperatures_k,
        before.transport.dry_y_hexane,
        before.last_total_stefan_fluxes_mol_m2_s,
        time_guess,
        before.last_interface_temperature_k,
    )
    seed = fi.FaceEventSeed(
        candidate,
        tuple(max(abs(v), 0.01) for v in candidate.dry_total_stefan_fluxes_mol_m2_s),
        "native exact-face candidate from accepted partial state",
    )
    pressure = before.transport.config.dry.pressure_pa
    config = before.transport.config
    radius = before.transport.geometry.master_grid.R
    nw, nd, nf = layout.wet_piece_count, layout.dry_piece_count, layout.dry_face_count
    holder = {
        "x": np.array([surface.temperature_k, 0.0, surface.liquid_water_kg_m2]),
        "dt": time_guess,
    }
    # beta is a numerical coordinate; film_at_beta reconstructs the exact
    # native total mass flux. It is not a new transfer law.
    surface_initial = np.array([surface.temperature_k, surface.beta, surface.liquid_water_kg_m2])
    surface_pore = replace(config.dry.pore, w_o=before.transport.oil_fraction_labels[-1])
    if gas_only_surface:
        if carry["external_water_loading_kg_m2"] != 0.0:
            raise ValueError("gas-only arrival requires absent external water")
        lo, hi = core.cp.gas_only_composition_interval(
            surface.temperature_k, pressure, surface_pore
        ).y_hexane_bounds
        surface_initial[-1] = (surface.surface_y_hexane - lo) / (hi - lo)
    mass_scale = max(carry["external_water_loading_kg_m2"] / time_guess, 0.05)
    heat_scale = max(
        abs(carry["external_water_energy_j"]) / (4 * math.pi * radius**2 * time_guess), 5000.0
    )
    work = 0
    final_boundary = None

    def native_candidate(values):
        return fi._decode_candidate(values, layout, chart, seed.stefan_flux_scales_mol_m2_s)

    def assemble(values):
        nonlocal work, final_boundary
        c = native_candidate(values[:-3])
        holder["dt"] = c.event_time_s
        clock = before.transport.time_s + c.event_time_s
        final_boundary = (
            boundary_factory(clock) if boundary_factory else bath_boundary(bath, radius, clock)
        )
        ts, beta, third = values[-3:]
        liquid = third
        y = None
        if gas_only_surface:
            if not 0 < third < 1:
                raise ValueError("gas-only arrival surface fraction must be interior")
            lo, hi = core.cp.gas_only_composition_interval(
                ts, pressure, surface_pore
            ).y_hexane_bounds
            y = lo + (hi - lo) * third
            holder["y_hexane"] = y
            liquid = 0.0
        holder["x"] = np.array(
            [
                ts,
                driver.film_at_beta(final_boundary, config, radius, ts, beta, surface_y_hexane=y)[
                    0
                ],
                liquid,
            ]
        )
        work += 1
        return face.assemble_face_arrival(
            before.transport, c, final_boundary, enforce_reduced_film_thresholds=False
        )

    def fraction(t, y, p):
        lo, hi = core.cp.gas_only_composition_interval(t, pressure, p).y_hexane_bounds
        return (y - lo) / (hi - lo)

    def decode(values):
        c = native_candidate(values[:-3])
        c = replace(
            c,
            dry_y_hexane=tuple(
                fraction(t, y, p)
                for t, y, p in zip(c.dry_temperatures_k, c.dry_y_hexane, chart.dry_pores)
            ),
        )
        return np.r_[c.vector(), values[-3:]]

    def encode(values):
        c = face.FaceArrivalUnknowns.from_vector(layout, values[:-3])
        ys = []
        for t, q, p in zip(c.dry_temperatures_k, c.dry_y_hexane, chart.dry_pores):
            if not 0 < q < 1:
                raise ValueError("gas fraction outside native open interval")
            lo, hi = core.cp.gas_only_composition_interval(t, pressure, p).y_hexane_bounds
            ys.append(lo + (hi - lo) * q)
        c = replace(c, dry_y_hexane=tuple(ys))
        return np.r_[fi._encode_candidate(c, chart, seed.stefan_flux_scales_mol_m2_s), values[-3:]]

    tl, tu = chart.dry_temperature
    pairs = (
        [chart.wet_temperature] * nw
        + [chart.wet_water] * nw
        + [chart.dry_temperature] * nd
        + [(0.0, 1.0)] * nd
        + [(-math.inf, math.inf)] * nf
        + [chart.event_time, chart.interface_temperature]
        + [
            (tl, tu),
            (-3.0, 3.0),
            (
                0.0,
                1.0 if gas_only_surface else max(0.05, 4 * carry["external_water_loading_kg_m2"]),
            ),
        ]
    )
    lower, upper = core.ci._coordinate_solver_bounds(
        layout, before.controls.maximum_logit_magnitude
    )
    x = np.r_[
        fi._encode_candidate(candidate, chart, seed.stefan_flux_scales_mol_m2_s), surface_initial
    ]
    with (
        study.installed_boundary_adapter(numerical_cap_seam=True),
        driver.liquid_internal_energy_node(
            holder,
            lambda: holder["dt"],
            carry["external_water_loading_kg_m2"],
            carry["external_water_temperature_k"],
            pressure,
        ),
    ):
        initial = assemble(x)
        scales = fi._residual_scales(before.transport, initial)

        def residual(values):
            a = assemble(values)
            s = a.dry_face_fluxes[-1].surface_film_audit
            return np.r_[
                np.asarray(fi._datum_covariant_residual_vector(a.residuals)) / scales.vector,
                s.water_node_residual_kg_m2_s / mass_scale,
                s.hexane_node_residual_kg_m2_s / mass_scale,
                s.energy_node_residual_w_m2 / heat_scale,
            ]

        solution = study.numerics.physical_chart_solve(
            residual,
            x,
            decode=decode,
            encode=encode,
            physical_bounds=tuple(zip(*pairs)),
            temperature_indices=[
                *range(nw),
                *range(2 * nw, 2 * nw + nd),
                len(pairs) - 4,
                len(pairs) - 3,
            ],
            bounds=(np.r_[lower, tl, -3.0, 0.0], np.r_[upper, tu, 3.0, pairs[-1][1]]),
            max_nfev=budget,
            gtol=before.controls.nonlinear_residual_tolerance,
            **({"difference_mode": difference_mode} if gas_only_surface else {}),
        )
        values = residual(solution.x)
        assembly = assemble(solution.x)
        maximum = float(np.max(np.abs(values)))
        condition = float(np.linalg.cond(solution.jac))
        if not solution.success or maximum > before.controls.nonlinear_residual_tolerance:
            raise RuntimeError(f"joint native face root refused: {maximum}")
        if condition > before.controls.maximum_condition_proxy:
            raise RuntimeError("native face conditioning refused")
        surface = assembly.dry_face_fluxes[-1].surface_film_audit
        if assembly.ledger.minimum_entropy_production_w_m3_k < 0:
            raise RuntimeError("native face entropy refused")
        if abs(surface.beta) > 3.0:
            raise RuntimeError("water-node beta band refused")
        step = fi._accept_event(
            before,
            controls,
            final_boundary,
            seed,
            assembly,
            scales,
            tuple(values[:-3]),
            work,
            solution.rejected_trial_evaluations,
            str(solution.message),
            condition,
            chart,
            None,
        )
        if (
            max(
                step.ledger.maximum_step_ledger_residual,
                step.ledger.maximum_cumulative_ledger_residual,
            )
            > before.controls.ledger_tolerance
        ):
            raise RuntimeError("native face conservation refused")
        core.ci.certify_explicit_wet_retained_water_applicability(
            before.controls,
            step.after.wet_retained_water_loadings,
            (0.0,) * nw,
            before.transport.config.wet.luikov,
        )
        return (
            step,
            surface,
            {
                "maximum_scaled_residual": maximum,
                "physical_coordinate_audit": solution.physical_coordinate_audit,
            },
        )


@contextmanager
def resolved_water_boundary(config, radius, carry, surface, dt, boundary, time_supplier=None):
    """Eliminate the same three conservative surface rows inside native solves.

    This is a nested nonlinear solve, not a split time step. Old water mass
    and energy are immutable throughout all body trials and Jacobian probes.
    """
    holder = {"x": np.array([surface.temperature_k, 0.0, surface.liquid_water_kg_m2])}
    guess = np.array([surface.temperature_k, surface.beta, surface.liquid_water_kg_m2])
    gas_only = carry["external_water_loading_kg_m2"] == 0.0
    if gas_only:
        lo, hi = core.cp.gas_only_composition_interval(
            surface.temperature_k, config.dry.pressure_pa, config.dry.pore
        ).y_hexane_bounds
        guess[-1] = (surface.surface_y_hexane - lo) / (hi - lo)
    audit = {
        "surface_residual_calls": 0,
        "surface_solves": 0,
        "maximum_accepted_surface_residual": 0.0,
        "old_film_mass_changed_during_trials": False,
        "surface_tolerance": 1e-11,
    }
    mass_scale = max(carry["external_water_loading_kg_m2"] / dt, 0.05)
    heat_scale = max(abs(carry["external_water_energy_j"]) / (4 * math.pi * radius**2 * dt), 5000.0)
    tl, tu = config.dry.conditioned_temperature_domain.solver_bounds_k
    with driver.liquid_internal_energy_node(
        holder,
        dt if time_supplier is None else time_supplier,
        carry["external_water_loading_kg_m2"],
        carry["external_water_temperature_k"],
        config.dry.pressure_pa,
    ):
        original = core.cut._dry_face_fluxes

        def resolved(*args, **kwargs):
            nonlocal guess
            actual_boundary = args[5]
            if not isinstance(actual_boundary, study.FiniteFilmStudyBoundary):
                return original(*args, **kwargs)
            faces = None

            def residual(values):
                nonlocal faces
                ts, beta, third = values
                liquid, y = third, None
                if gas_only:
                    if not 0 < third < 1:
                        raise ValueError("gas-only surface fraction must be interior")
                    lo, hi = core.cp.gas_only_composition_interval(
                        ts, config.dry.pressure_pa, config.dry.pore
                    ).y_hexane_bounds
                    y = lo + (hi - lo) * third
                    holder["y_hexane"] = y
                    liquid = 0.0
                holder["x"] = np.array(
                    [
                        ts,
                        driver.film_at_beta(
                            actual_boundary, config, radius, ts, beta, surface_y_hexane=y
                        )[0],
                        liquid,
                    ]
                )
                audit["surface_residual_calls"] += 1
                faces = original(*args, **kwargs)
                s = faces[-1].surface_film_audit
                return np.array(
                    [
                        s.water_node_residual_kg_m2_s / mass_scale,
                        s.hexane_node_residual_kg_m2_s / mass_scale,
                        s.energy_node_residual_w_m2 / heat_scale,
                    ]
                )

            solution = study.numerics.centered_temperature_solve(
                study.numerics.domain_aware_newton,
                residual,
                guess,
                temperature_index=0,
                bounds=(
                    [tl, -3.0, 0.0],
                    [
                        tu,
                        3.0,
                        1.0 if gas_only else max(0.05, 4 * carry["external_water_loading_kg_m2"]),
                    ],
                ),
                max_nfev=160,
                gtol=1e-11,
            )
            values = residual(solution.x)
            maximum = float(np.max(np.abs(values)))
            if not solution.success or maximum > 1e-11:
                raise RuntimeError(f"nested conservative water surface refused: {maximum}")
            guess = solution.x.copy()
            audit["surface_solves"] += 1
            audit["maximum_accepted_surface_residual"] = max(
                audit["maximum_accepted_surface_residual"], maximum
            )
            return faces

        core.cut._dry_face_fluxes = resolved
        try:
            yield audit
        finally:
            core.cut._dry_face_fluxes = original


def advance_departure(
    before, carry, surface, bath, dt, *, boundary_factory=None, tangent_residual_bound=1e-10
):
    radius = before.geometry.master_grid.R
    clock = before.time_s + dt
    boundary = boundary_factory(clock) if boundary_factory else bath_boundary(bath, radius, clock)
    old_bound = os.environ.get(tangent.FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE)
    # Existing owner-ruled engineering declaration; no outward root admitted.
    os.environ[tangent.FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE] = str(
        tangent_residual_bound
    )
    try:
        t = tangent.assemble_face_tangent(before, tangent.candidate_from_committed_face(before))
        with (
            study.installed_boundary_adapter(numerical_cap_seam=True),
            resolved_water_boundary(
                before.config, radius, carry, surface, dt, boundary
            ) as surface_work,
        ):
            overlap = tangent.continue_to_positive_sweep(t, dt, boundary)
            seed = di.FaceDepartureSeed(
                overlap.candidate,
                tuple(
                    max(abs(v), 0.01) for v in overlap.candidate.dry_total_stefan_fluxes_mol_m2_s
                ),
                "native tangent as numerical seed only",
                overlap.newborn_trace_gradients,
            )
            step = di.solve_face_departure(before, before.solver_controls, dt, boundary, seed)
        surface = step.assembly.dry_face_fluxes[-1].surface_film_audit
        return (
            step,
            surface,
            {
                "nested_surface": surface_work,
                "tangent": asdict(t.ledger),
                "tangent_residual_bound": tangent_residual_bound,
                "tangent_used_as_accepted_time_step": False,
            },
        )
    finally:
        if old_bound is None:
            os.environ.pop(tangent.FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE, None)
        else:
            os.environ[tangent.FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE] = old_bound


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--bath", type=Path, default=core.DESTINY_PATH.with_name("destiny_v3_exit105.json")
    )
    parser.add_argument("--budget", type=int, default=1800)
    parser.add_argument("--mode", choices=("arrival", "departure"), default="arrival")
    parser.add_argument("--dt", type=float, default=0.05)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite evidence")
    source_paths = [
        Path(__file__).resolve(),
        Path(driver.__file__),
        Path(driver.node.__file__),
        Path(study.__file__),
        Path(study.numerics.__file__),
        Path(driver.path_bounds.__file__),
        args.bath.resolve(),
        *sorted((driver.REPO_ROOT / "src/dtdc_simulator/core2").rglob("*.py")),
    ]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    checkpoint = pickle.loads(args.resume.read_bytes())
    before = checkpoint["current"]
    record = {
        "accepted": False,
        "complete_dt_destiny": False,
        "source_sha256": hashes,
        "resume": str(args.resume),
    }
    for p in source_paths:
        if "src" in p.parts and checkpoint["source_sha256"].get(str(p)) not in (
            None,
            hashes[str(p)],
        ):
            raise ValueError("native source differs from accepted checkpoint")
    try:
        journal = []
        with (
            driver.bounded_path_proof(3000000, 120.0) as proof,
            driver.path_bounds.monotone_water_path_enclosure(journal),
            study.magnitude_storage_scale(),
        ):
            bath = json.loads(args.bath.read_text())
            if args.mode == "arrival":
                step, surface, audit = advance_arrival(
                    before, checkpoint["carry"], checkpoint["surface"], bath, budget=args.budget
                )
                config = before.transport.config
                controls = before.controls
            else:
                step, surface, audit = advance_departure(
                    before, checkpoint["carry"], checkpoint["surface"], bath, args.dt
                )
                config = before.config
                controls = before.solver_controls
        combined = driver.combined_ledger(before, config, checkpoint["carry"], step, surface)
        if combined["maximum_normalized_residual"] > controls.ledger_tolerance:
            raise RuntimeError("combined particle/water-film face ledger refused")
        carry = {
            **checkpoint["carry"],
            **combined,
            "external_water_loading_kg_m2": surface.liquid_water_kg_m2,
        }
        record.update(
            accepted=True,
            status=f"NATIVE_EXACT_FACE_{args.mode.upper()}_ACCEPTED",
            after=asdict(step.after),
            ledger=asdict(step.ledger),
            surface=asdict(surface),
            audit=audit,
            proof=proof,
            combined=combined,
        )
        output = args.out.with_suffix(".accepted.pkl")
        output.write_bytes(
            pickle.dumps(
                {
                    **checkpoint,
                    "current": step.after,
                    "surface": surface,
                    "carry": carry,
                    "source_sha256": {**checkpoint["source_sha256"], **hashes},
                    "result_path": str(args.out),
                },
                protocol=5,
            )
        )
        record["accepted_checkpoint"] = str(output)
    except Exception as exc:
        record.update(
            error=str(exc),
            traceback=traceback.format_exc(),
            optimizer_debug=getattr(exc, "optimizer_debug", None),
            exact_accepted_state_preserved=True,
        )
    record["sources_unchanged_during_run"] = all(
        hashlib.sha256(p.read_bytes()).hexdigest() == hashes[str(p)] for p in source_paths
    )
    if not record["sources_unchanged_during_run"]:
        record.update(accepted=False, status="SOURCE_CHANGED_DURING_RUN")
    record["stamp"] = driver.era_stamp(source_paths)
    args.out.write_text(json.dumps(record, indent=2))
    print(json.dumps({k: record.get(k) for k in ("accepted", "status", "error")}))


if __name__ == "__main__":
    main()
