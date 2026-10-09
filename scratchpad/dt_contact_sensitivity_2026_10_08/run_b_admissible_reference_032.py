"""Initialize numerical chart and row normalization from the admissible seed.

The exact-tangent physical column scales are retained. Only the reference
vector used to initialize a chart and normalize residual rows changes.
All physical trial assemblies and native final acceptance remain unchanged.
"""
import hashlib
import json
from pathlib import Path
import runpy
import sys

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / "src"))
sys.path.insert(0, str(root / "paper/analysis/harness"))
import experiment_core2_finite_wetting_b as engine

di = engine.face.di
original_stage = di._solve_direct_stage
original_vector = di._tangent_reference_physical_vector
active = []
out = Path(sys.argv[sys.argv.index("--out") + 1])
audit_path = out.with_suffix(".reference_initialization_audit.json")
assert not audit_path.exists()
audits = []


def reference_vector(layout, reference):
    if not active:
        return original_vector(layout, reference)
    before, dt_s, seed = active[-1]
    gradients = di._seed_trace_gradients(seed, before)
    result = di._physical_departure_vector(seed.candidate, gradients, before, dt_s)
    assert len(result) == layout.unknown_count
    audits.append({
        "before_time_s": before.time_s,
        "dt_s": dt_s,
        "reference_source": "uncommitted admissible numerical seed",
        "original_tangent_vector": original_vector(layout, reference),
        "admissible_reference_vector": result,
        "tangent_based_physical_column_scales_changed": False,
        "native_residual_scale_function_changed": False,
        "physical_residual_assembly_changed": False,
        "native_final_phase_balance_conditioning_changed": False,
    })
    audit_path.write_text(json.dumps(audits, indent=2))
    print(json.dumps({"admissible_reference_initialization": {
        "before_time_s": before.time_s, "dt_s": dt_s,
        "physical_equations_and_final_gates_changed": False,
    }}), flush=True)
    return result


def stage(*args, **kwargs):
    active.append((args[0], args[2], args[4]))
    try:
        return original_stage(*args, **kwargs)
    finally:
        active.pop()


entry = Path(__file__).with_name("run_b_neighbor_seed_policy_032.py")
policy = out.with_suffix(".reference_policy.json")
assert not policy.exists()
policy.write_text(json.dumps({
    "change": "Numerical chart and residual-row normalization reference uses an admissible seed rather than an inadmissible tangent prediction",
    "basis": "Observed zero-evaluation initialization failures; owner RULE 2 authorizes numerical recovery with unchanged physical equations",
    "physical_law_changed": False,
    "native_final_acceptance_changed": False,
    "tangent_based_physical_column_scales_retained": True,
    "component_energy_tolerance": 1e-10,
    "source_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (Path(__file__), entry, entry.with_name("run_b_neighbor_departure_seed_032.py"),
                                root / "src/dtdc_simulator/core2/particle/cut_face_departure_integrator.py")},
}, indent=2))
di._tangent_reference_physical_vector = reference_vector
di._solve_direct_stage = stage
try:
    runpy.run_path(str(entry), run_name="__main__")
finally:
    di._tangent_reference_physical_vector = original_vector
    di._solve_direct_stage = original_stage
