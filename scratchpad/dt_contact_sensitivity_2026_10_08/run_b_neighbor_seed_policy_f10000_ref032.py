"""Reuse the existing DT seed-search tolerance with a neighboring root guess."""
import hashlib
import json
import os
from pathlib import Path
import runpy
import sys

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / "src"))
sys.path.insert(0, str(root / "paper/analysis/harness"))
from dtdc_simulator.core2.particle import cut_face_tangent as tangent

entry = Path(__file__).with_name("run_b_neighbor_departure_seed_f10000_ref032.py")
out = Path(sys.argv[sys.argv.index("--out") + 1])
policy = out.with_suffix(".neighbor_seed_policy.json")
assert not policy.exists()
policy.write_text(json.dumps({
    "seed_only_tangent_residual_bound": 3e-10,
    "bound_basis": "Existing B departure policy and owner RULE 2; native final conditioning and W/H/U acceptance unchanged",
    "component_energy_tolerance": 1e-10,
    "neighboring_result_is_numerical_guess_only": True,
    "accepted_physical_coefficient": 0.32,
    "source_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__), entry)},
}, indent=2))
name = tangent.FACE_TANGENT_RESIDUAL_BOUND_ENVIRONMENT_VARIABLE
prior = os.environ.get(name)
os.environ[name] = "3e-10"
try:
    runpy.run_path(str(entry), run_name="__main__")
finally:
    if prior is None:
        os.environ.pop(name, None)
    else:
        os.environ[name] = prior
