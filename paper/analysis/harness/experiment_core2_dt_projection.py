"""Conservative native dry-state projection after measured sequential time refinement.

The time estimate is a disclosed numerical truncation, not the native parked
same-input extinction certificate. Its phase endpoint and all balances are native.
"""

from contextlib import contextmanager
from dataclasses import asdict, replace
from pathlib import Path
import argparse
import hashlib
import json
import pickle
import traceback
import numpy as np
import experiment_core2_dt_face as boundary_driver
from dtdc_simulator.core2.particle import cut_extinction_event as extinction
from dtdc_simulator.core2.particle import cut_extinction_projection as projection

driver = boundary_driver.driver
core = driver.core
study = driver.study


@contextmanager
def water_surface_at_dry_endpoint(before, carry, surface, dt, boundary):
    config = before.transport.config
    grid = before.transport.geometry.master_grid
    holder = {"x": np.array([surface.temperature_k, 0.0, surface.liquid_water_kg_m2])}
    guess = np.array([surface.temperature_k, surface.beta, surface.liquid_water_kg_m2])
    gas_only = carry["external_water_loading_kg_m2"] == 0.0
    if gas_only:
        lo, hi = core.cp.gas_only_composition_interval(
            surface.temperature_k, config.dry.pressure_pa, config.dry.pore
        ).y_hexane_bounds
        guess[-1] = (surface.surface_y_hexane - lo) / (hi - lo)
    work = {"surface_solves": 0, "surface_residual_calls": 0, "maximum_surface_residual": 0.0}
    mass_scale = max(carry["external_water_loading_kg_m2"] / dt, 0.05)
    heat_scale = max(abs(carry["external_water_energy_j"]) / (grid.areas[-1] * dt), 5000.0)
    tl, tu = config.dry.conditioned_temperature_domain.solver_bounds_k
    original = extinction._fully_dry_face_fluxes
    with driver.liquid_internal_energy_node(
        holder,
        dt,
        carry["external_water_loading_kg_m2"],
        carry["external_water_temperature_k"],
        config.dry.pressure_pa,
    ):

        def faces(old, candidate, cells, actual_boundary, **kwargs):
            nonlocal guess
            if not isinstance(actual_boundary, study.FiniteFilmStudyBoundary):
                return original(old, candidate, cells, actual_boundary, **kwargs)
            cell = cells[-1]
            pore = replace(config.dry.pore, w_o=before.transport.oil_fraction_labels[-1])
            face = None

            def residual(values):
                nonlocal face
                ts, beta, third = values
                liquid, y = third, None
                if gas_only:
                    if not 0 < third < 1:
                        raise ValueError("gas-only endpoint surface fraction must be interior")
                    lo, hi = core.cp.gas_only_composition_interval(
                        ts, config.dry.pressure_pa, pore
                    ).y_hexane_bounds
                    y = lo + (hi - lo) * third
                    liquid = 0.0
                holder["x"] = np.array(
                    [
                        ts,
                        driver.film_at_beta(
                            actual_boundary, config, grid.R, ts, beta, surface_y_hexane=y
                        )[0],
                        liquid,
                    ]
                )
                face = driver.node.evaluate_water_surface(
                    holder,
                    dt,
                    carry["external_water_loading_kg_m2"],
                    carry["external_water_temperature_k"],
                    cell,
                    grid.R - grid.centers[-1],
                    candidate.positive_area_total_stefan_fluxes_mol_m2_s[-1],
                    pore,
                    config,
                    grid.R,
                    actual_boundary,
                    surface_y_hexane=y,
                )
                s = face.surface_film_audit
                work["surface_residual_calls"] += 1
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
            maximum = float(np.max(np.abs(residual(solution.x))))
            if not solution.success or maximum > 1e-11:
                raise RuntimeError("native endpoint water node refused")
            guess = solution.x.copy()
            work["surface_solves"] += 1
            work["maximum_surface_residual"] = max(work["maximum_surface_residual"], maximum)
            resolved = core.ct.DirichletPoreBoundary(
                cell.temperature_k, cell.y_hexane, "native internal faces only"
            )
            internal = original(old, candidate, cells, resolved, **kwargs)
            return (*internal[:-1], face)

        extinction._fully_dry_face_fluxes = faces
        try:
            yield work
        finally:
            extinction._fully_dry_face_fluxes = original


def time_estimate(before, checkpoint):
    data = json.loads(Path(checkpoint["result_path"]).read_text())
    if not data["sources_unchanged_during_run"]:
        raise ValueError("target evidence provenance failed")
    rows = [r for r in data["attempts"] if r["accepted"]]
    if len(rows) < 4 or rows[-1]["time_s"] != before.transport.time_s:
        raise ValueError("four linked accepted time refinements required")
    if rows[-1]["target_z"] != before.transport.geometry.front.z or rows[-1]["target_z"] > 1e-10:
        raise ValueError("last accepted wet volume is not below declared 1e-10 fraction")
    d = [r["ledger"]["dt_s"] for r in rows[-4:]]
    ratios = [b / a for a, b in zip(d, d[1:])]
    if not all(0 < q < 0.82 for q in ratios):
        raise ValueError("time increments do not meet prospective 0.82 contraction")
    estimated = d[-1] * max(ratios) / (1 - max(ratios))
    bound = d[-1] * 0.82 / (1 - 0.82)
    if bound > 0.002:
        raise ValueError("remaining numerical event-time bound exceeds declared 2 ms")
    return estimated, {
        "method": "sequential native positive-target real-time refinement",
        "last_four_step_durations_s": d,
        "observed_ratios": ratios,
        "assumed_future_ratio_upper": 0.82,
        "estimated_remaining_time_s": estimated,
        "conditional_remaining_time_bound_s": bound,
        "declared_time_tolerance_s": 0.002,
        "last_wet_volume_fraction": rows[-1]["target_z"],
        "native_same_input_limit_certificate_claimed": False,
        "new_physical_law": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite evidence")
    checkpoint = pickle.loads(args.resume.read_bytes())
    before = checkpoint["current"]
    bath_path = core.DESTINY_PATH.with_name("destiny_v3_exit105.json")
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
    record = {"accepted": False, "complete_dt_destiny": False, "sources_at_entry": hashes}
    try:
        dt, estimate = time_estimate(before, checkpoint)
        record["event_time_estimate"] = estimate
        boundary = boundary_driver.bath_boundary(
            json.loads(bath_path.read_text()),
            before.transport.geometry.master_grid.R,
            before.transport.time_s + dt,
        )
        candidate = extinction.ExtinctionUnknowns(
            before.transport.dry_temperatures_k,
            before.transport.dry_y_hexane,
            before.last_total_stefan_fluxes_mol_m2_s[1:],
            dt,
        )
        seed = projection.ExtinctionProjectionSeed(
            candidate,
            tuple(max(abs(v), 0.01) for v in candidate.positive_area_total_stefan_fluxes_mol_m2_s),
            "native dry endpoint from actual diminishing core and measured time refinement",
        )
        with (
            driver.bounded_path_proof(3000000, 120.0) as proof,
            driver.path_bounds.monotone_water_path_enclosure([]),
            study.magnitude_storage_scale(),
            water_surface_at_dry_endpoint(
                before, checkpoint["carry"], checkpoint["surface"], dt, boundary
            ) as work,
        ):
            step = projection.solve_fixed_time_extinction_projection(before, dt, boundary, seed)
        surface = step.assembly.dry_face_fluxes[-1].surface_film_audit
        combined = driver.combined_ledger(
            before, before.transport.config, checkpoint["carry"], step, surface
        )
        if combined["maximum_normalized_residual"] > before.controls.ledger_tolerance:
            raise RuntimeError("combined dry-endpoint ledger refused")
        carry = {
            **checkpoint["carry"],
            **combined,
            "external_water_loading_kg_m2": surface.liquid_water_kg_m2,
        }
        record.update(
            accepted=True,
            status="NATIVE_CONSERVATIVE_DRY_ENDPOINT_ACCEPTED",
            after=asdict(step.after),
            ledger=asdict(step.ledger),
            surface=asdict(surface),
            combined=combined,
            proof=proof,
            surface_work=work,
        )
        output = args.out.with_suffix(".accepted.pkl")
        output.write_bytes(
            pickle.dumps(
                {
                    **checkpoint,
                    "current": step.after,
                    "surface": surface,
                    "carry": carry,
                    "cut_config": before.transport.config,
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
        )
    record["sources_unchanged_during_run"] = all(
        hashlib.sha256(p.read_bytes()).hexdigest() == hashes[str(p.resolve())] for p in sources
    )
    record["stamp"] = driver.era_stamp(sources)
    args.out.write_text(json.dumps(record, indent=2))
    print(json.dumps({k: record.get(k) for k in ("accepted", "status", "error")}))


if __name__ == "__main__":
    main()
