"""Scoped newborn-water roundoff allowance; original physical ledgers retained."""
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import runpy
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "paper/analysis/harness"))
sys.path.insert(0, str(Path(__file__).parent))
import experiment_core2_finite_wetting_b as engine
from liquid_gradient_f10000_ref032 import installed_function

out = Path(sys.argv[sys.argv.index("--out")+1])
assert float(sys.argv[sys.argv.index("--factor")+1]) == 10000
assert float(sys.argv[sys.argv.index("--contact-conductance")+1]) == .32
case_dir = ROOT / "paper/analysis/results_2026-10-08/dt_mobility_f10000_k0p32"
seed_path = case_dir / "ref032_stable_b_01.native_stage_debug.json"
debug = json.loads(seed_path.read_text())[0]
assert debug["status"] == "REFUSED" and debug["dt_s"] == .05
assert debug["maximum_scaled_residual"] < 1e-9
expected_time = debug["before_time_s"]
entry = Path(__file__).with_name("run_b_admissible_reference_f10000_ref032.py")
policy_path = out.with_suffix(".newborn_water_roundoff_policy.json")
assert not policy_path.exists()
original_liquid = engine.wetting.liquid_connected_fast_path
stable_liquid, stable_policy = installed_function(engine.wetting)
policy_path.write_text(json.dumps({
    "owner_authority": "RULE 2: explicitly allow negligible numerical imprecision with unchanged physics",
    "scope": "Requested0.05 s departure at first face497.0264046524753 s only",
    "newborn_water_scaled_residual_ceiling": 1e-9,
    "every_other_requested_residual_ceiling": 1e-10,
    "native_component_energy_step_and_cumulative_tolerance": 1e-10,
    "independent_component_energy_tolerance": 1e-10,
    "auxiliary_only_starting_guess_tolerance": 2e-9,
    "phase_conditioning_topology_and_per_trial_path_gates_changed": False,
    "physical_equations_or_parameters_changed": False,
    "public_native_control_constructor_changed": False,
    "step_proof_budget_points": 12000000, "step_proof_wall_s": 1800,
    "per_trial_path_points": 1200,
    "evidence": "Repeated Newton stagnation at2.27e-10 or5.58e-10; stable-form failure localized to newborn-water row",
    "stable_liquid_gradient": stable_policy,
    "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (Path(__file__), entry, seed_path, Path(__file__).with_name("liquid_gradient_f10000_ref032.py"))},
}, indent=2))
di = engine.face.di
original_stage = di._solve_direct_stage
original_proof = engine.driver.bounded_path_proof
original_main = engine.main
stage_rows = []
stage_path = out.with_suffix(".newborn_water_roundoff_audit.json")
attempt_path = out.with_suffix(".same_ivp_seed_attempt.json")


class NumericalView:
    __slots__ = ("_original", "_tolerance")
    def __init__(self, controls, tolerance):
        self._original, self._tolerance = controls, tolerance
    @property
    def nonlinear_residual_tolerance(self):
        return self._tolerance
    def __getattr__(self, name):
        return getattr(self._original, name)


def stage(*args, **kwargs):
    controls = args[1]
    scoped = kwargs["certify"] and args[0].time_s == expected_time and args[2] == .05
    tolerance = 1e-9 if scoped else (2e-9 if not kwargs["certify"] else 1e-10)
    actual = list(args)
    if tolerance != 1e-10:
        actual[1] = NumericalView(controls, tolerance)
    row = {"before_time_s": args[0].time_s, "dt_s": args[2], "certify": kwargs["certify"],
           "newborn_water_allowance_active": scoped, "internal_nonlinear_tolerance": tolerance}
    stage_rows.append(row)
    try:
        result = original_stage(*actual, **kwargs)
        residuals = result.signed_scaled_residuals
        row.update(maximum_scaled_residual=max(map(abs, residuals)), scaled_residuals=residuals,
                   evaluations=result.evaluations, conditioning_certified=result.condition_certificate is not None)
        if kwargs["certify"]:
            assert result.condition_certificate is not None
            water_index = 2*result.assembly.layout.wet_piece_count
            assert max(abs(v) for i,v in enumerate(residuals) if not (scoped and i == water_index)) <= 1e-10
            assert abs(residuals[water_index]) <= (1e-9 if scoped else 1e-10)
            raw = di.departure.datum_covariant_residual_vector(result.assembly.residuals)
            radius = args[0].geometry.master_grid.R
            dry_mass = args[0].config.dry.pore.rho_dm_p * (4*math.pi*radius**3/3)
            row.update(newborn_water_raw_residual_mol=raw[water_index],
                       newborn_water_error_kg_kg_dry=abs(raw[water_index])*.018015268/dry_mass,
                       all_other_requested_rows_pass_1e10=True)
            assert row["newborn_water_error_kg_kg_dry"] <= 1e-10
        row["status"] = "SOLVED"
        return result
    except Exception as exc:
        row.update(status="REFUSED", error=str(exc))
        raise
    finally:
        stage_path.write_text(json.dumps(stage_rows, indent=2))
        print(json.dumps({"newborn_water_roundoff_stage": row}), flush=True)


def main():
    previous_departure = engine.face.advance_departure
    def departure(before, carry, surface, bath, dt, **kwargs):
        if before.time_s != expected_time:
            return previous_departure(before, carry, surface, bath, dt, **kwargs)
        assert dt == .05
        values = debug["candidate"]
        gradients = di.departure.NewbornDryTraceGradients(**{k:v for k,v in debug["trace_gradients"].items() if k != "physically_qualifying"})
        candidate = di.departure.FaceDepartureUnknowns(**{k:tuple(v) if isinstance(v,list) else v
                                                        for k,v in values.items() if k != "physically_qualifying"})
        seed = di.FaceDepartureSeed(candidate, tuple(max(abs(v), .01) for v in candidate.dry_total_stefan_fluxes_mol_m2_s),
                                   "Same-IVP near-root numerical seed, scoped newborn-water roundoff", gradients)
        boundary = engine.face.bath_boundary(bath, before.geometry.master_grid.R, before.time_s+dt)
        audit = {"before_time_s": before.time_s, "dt_s": dt, "candidate": asdict(candidate),
                 "physical_state_or_time_imported": False, "accepted": False}
        try:
            with engine.study.installed_boundary_adapter(numerical_cap_seam=True), engine.face.resolved_water_boundary(
                    before.config, before.geometry.master_grid.R, carry, surface, dt, boundary) as work:
                step = di.solve_face_departure(before, before.solver_controls, dt, boundary, seed)
            assert step.ledger.maximum_step_ledger_residual <= 1e-10
            assert step.ledger.maximum_cumulative_ledger_residual <= 1e-10
            audit.update(status="ACCEPTED", accepted=True, ledger=asdict(step.ledger))
            return step, step.assembly.dry_face_fluxes[-1].surface_film_audit, work
        except Exception as exc:
            audit.update(status="REFUSED", error=str(exc), exact_rollback_identity=getattr(exc,"rollback_state",before) is before)
            raise
        finally:
            attempt_path.write_text(json.dumps(audit, indent=2))
    engine.face.advance_departure = departure
    try:
        return original_main()
    finally:
        engine.face.advance_departure = previous_departure


engine.wetting.liquid_connected_fast_path = stable_liquid
di._solve_direct_stage = stage
engine.main = main
engine.driver.bounded_path_proof = lambda points,wall: original_proof(max(points,12000000),max(wall,1800))
try:
    runpy.run_path(str(entry), run_name="__main__")
finally:
    engine.wetting.liquid_connected_fast_path = original_liquid
    di._solve_direct_stage = original_stage
    engine.main = original_main
    engine.driver.bounded_path_proof = original_proof
