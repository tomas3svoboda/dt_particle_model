"""Use a neighboring completed calculation as a numerical departure guess only.

The physical before-state, contact coefficient, bath and native corrector
remain those of K=0.32. No neighboring accepted material or time is imported.
"""
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / "src"))
sys.path.insert(0, str(root / "paper/analysis/harness"))
import experiment_core2_finite_wetting_b as engine

reference_rows = []
reference_attempts = []
reference_hashes = {}
reference_case = root / "paper/analysis/results_2026-10-08/dt_contact_0p32/completed_case.json"
case = json.loads(reference_case.read_text())
assert case["status"] == "COMPLETE_AND_AUDITED" and case["K_l"] == 0.32
reference_hashes[str(reference_case.relative_to(root))] = hashlib.sha256(reference_case.read_bytes()).hexdigest()
for relative in case["b_parts"]:
    path = root / relative
    data = json.loads(path.read_text())
    assert data["sources_unchanged_during_run"]
    reference_hashes[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    reference_rows.extend(r for r in data["history"] if r["kind"] != "resumed")
    reference_attempts.extend(data["attempts"])
assert float(sys.argv[sys.argv.index("--factor") + 1]) == 10000
assert float(sys.argv[sys.argv.index("--contact-conductance") + 1]) == 0.32
reference_before = next(r for r in reference_rows if r["kind"] == "face_arrival")["after_transport"]
reference_after = next(r for r in reference_rows if r["kind"] == "face_departure")["after_transport"]
reference_ledger = next(r for r in reference_attempts if r["kind"] == "face_departure")["ledger"]
reference_dt = reference_ledger["dt_s"]
g_t = reference_ledger["newborn_trace_gradients"]["temperature_gradient_k_m"]
g_y = reference_ledger["newborn_trace_gradients"]["composition_gradient_m_inv"]
reference_geometry = engine.core.cut.cg.partition_master_grid(
    engine.core.backend.uniform_grid(2, 0.0015), z=reference_after["geometry"]["front"]["z"])
reference_distance = engine.face.di._newborn_center_distance(reference_geometry, 0)
reference_gamma = math.fsum((reference_after["dry_temperatures_k"][0], -g_t * reference_distance))
out_path = Path(sys.argv[sys.argv.index("--out") + 1])
audit_path = out_path.with_suffix(".neighbor_seed_attempts.json")
assert not audit_path.exists()
attempts = []
original = engine.face.advance_departure

def departure(before, carry, surface, bath, dt, **kwargs):
    for duration in (dt, dt / 4, dt / 16, dt / 64):
        boundary = engine.face.bath_boundary(bath, before.geometry.master_grid.R, before.time_s + duration)
        fraction = duration / reference_dt
        z = before.geometry.front.z - fraction * (reference_before["geometry"]["front"]["z"] - reference_after["geometry"]["front"]["z"])
        geometry = engine.core.cut.cg.partition_master_grid(before.geometry.master_grid, z=z)
        distance = engine.face.di._newborn_center_distance(geometry, before.arrival_face_index - 1)
        interface = engine.core.cut.evaluate_interface_state(reference_gamma, before.config,
                    before.historical_hexane_loadings[before.arrival_face_index - 1],
                    before.oil_fraction_labels[before.arrival_face_index - 1])
        def increments(name):
            return tuple(old + fraction * (new - prior) for old, new, prior in
                         zip(getattr(before, name), reference_after[name], reference_before[name]))
        dry_t = tuple(old + fraction * (new - prior) for old, new, prior in
                      zip(before.dry_temperatures_k, reference_after["dry_temperatures_k"][1:], reference_before["dry_temperatures_k"]))
        dry_y = tuple(old + fraction * (new - prior) for old, new, prior in
                      zip(before.dry_y_hexane, reference_after["dry_y_hexane"][1:], reference_before["dry_y_hexane"]))
        # Existing fluxes are guesses only; native RH/component rows solve them.
        fluxes = (before.dry_total_stefan_fluxes_mol_m2_s[0], *before.dry_total_stefan_fluxes_mol_m2_s)
        candidate = engine.face.di.departure.FaceDepartureUnknowns(
            increments("wet_temperatures_k"), increments("wet_retained_water_loadings"),
            (math.fsum((reference_gamma, g_t * distance)), *dry_t),
            (math.fsum((interface.y_hexane, g_y * distance)), *dry_y),
            fluxes, z, reference_gamma)
        gradients = engine.face.di.departure.NewbornDryTraceGradients(g_t, g_y)
        seed = engine.face.di.FaceDepartureSeed(candidate, tuple(max(abs(v), 0.01) for v in fluxes),
                        "neighbor K=0.32 factor1000 root as uncommitted initial guess; actual IVP is K=0.32 factor10000", gradients)
        entry = {"before_time_s": before.time_s, "requested_dt_s": duration, "reference_result_sha256": reference_hashes,
                 "candidate": asdict(candidate), "native_acceptance_changed": False,
                 "neighbor_material_or_time_imported": False}
        attempts.append(entry)
        try:
            with engine.study.installed_boundary_adapter(numerical_cap_seam=True), engine.face.resolved_water_boundary(
                    before.config, before.geometry.master_grid.R, carry, surface, duration, boundary) as work:
                step = engine.face.di.solve_face_departure(before, before.solver_controls, duration, boundary, seed)
        except Exception as error:
            entry.update(accepted=False, error=str(error), maximum_scaled_residual=getattr(error, "maximum_scaled_residual", None),
                         nonlinear_evaluations=getattr(error, "nonlinear_evaluations", None),
                         last_candidate=asdict(error.last_candidate) if getattr(error, "last_candidate", None) else None,
                         scaled_residuals=getattr(error, "scaled_residuals", None),
                         exact_rollback_identity=getattr(error, "rollback_state", before) is before)
            audit_path.write_text(json.dumps(attempts, indent=2))
            print(json.dumps({"neighbor_seed_retry": {k: v for k, v in entry.items() if k not in ("candidate", "last_candidate", "reference_result_sha256")}}), flush=True)
        else:
            entry.update(accepted=True, ledger=asdict(step.ledger))
            audit_path.write_text(json.dumps(attempts, indent=2))
            return step, step.assembly.dry_face_fluxes[-1].surface_film_audit, work
    raise RuntimeError("Neighbor-root departure guesses exhausted without advancing physical state")

engine.face.advance_departure = departure
try:
    engine.main()
finally:
    engine.face.advance_departure = original
