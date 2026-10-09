"""Give the native proof 600 seconds per step; keep its assembly budget."""
import hashlib
import json
from pathlib import Path
import runpy
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "paper/analysis/harness"))
import experiment_core2_finite_wetting_b as engine

out = Path(sys.argv[sys.argv.index("--out") + 1])
assert float(sys.argv[sys.argv.index("--contact-conductance") + 1]) == 0.32
entry = Path(__file__).with_name("run_b_admissible_reference_032.py")
policy = out.with_suffix(".proof_wall_policy.json")
assert not policy.exists()
policy.write_text(json.dumps({
    "owner_authority": "RULE 2: explicit numerical budget recovery with unchanged physics",
    "old_step_proof_wall_s": 120.0, "new_step_proof_wall_s": 600.0,
    "proof_assembly_budget_unchanged": 3000000,
    "component_energy_acceptance": 1e-10,
    "physical_equations_changed": False, "native_acceptance_tolerances_changed": False,
    "sources_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (Path(__file__), entry)},
}, indent=2))
original = engine.driver.bounded_path_proof

def longer_wall(assembly_budget, wall_s, *args, **kwargs):
    return original(assembly_budget, max(wall_s, 600.0), *args, **kwargs)

engine.driver.bounded_path_proof = longer_wall
try:
    runpy.run_path(str(entry), run_name="__main__")
finally:
    engine.driver.bounded_path_proof = original
