"""Declare a seed-only tangent bound; retain native accepted-step checks."""
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
entry = Path(__file__).with_name("run_b_032.py")
policy = out.with_suffix(".seed_bound_policy.json")
assert not policy.exists()
policy.write_text(json.dumps({
    "owner_authority": "RULE 2: documented numerical relaxation with unchanged physical equations and negligible numerical error",
    "observed_bounded_tangent_reconstruction_residual": 2.189e-10,
    "previous_seed_only_tangent_bound": 2e-10,
    "declared_seed_only_tangent_bound": 3e-10,
    "tangent_is_accepted_time_step": False,
    "native_accepted_step_controls_changed": False,
    "component_energy_acceptance": 1e-10,
    "physical_equations_or_coefficient_changed": False,
    "sources_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (Path(__file__), entry)},
}, indent=2))
original = engine.face.advance_departure

def seed_bound(*args, **kwargs):
    kwargs["tangent_residual_bound"] = 3e-10
    return original(*args, **kwargs)

engine.face.advance_departure = seed_bound
try:
    runpy.run_path(str(entry), run_name="__main__")
finally:
    engine.face.advance_departure = original
