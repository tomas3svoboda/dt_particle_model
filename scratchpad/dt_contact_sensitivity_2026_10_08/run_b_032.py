"""Compose existing same-IVP birth recovery and adaptive native departure.

This adds no physical coefficient, accepted-state substitution or gate change.
"""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import runpy
import sys

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / "src"))
sys.path.insert(0, str(root / "paper/analysis/harness"))
import experiment_core2_finite_wetting_b as engine

out = Path(sys.argv[sys.argv.index("--out") + 1])
coefficient = float(sys.argv[sys.argv.index("--contact-conductance") + 1])
assert coefficient == 0.32
out.parent.mkdir(parents=True, exist_ok=True)
original = engine.face.advance_departure
attempts = []
audit = out.with_suffix(".departure_attempts.json")
assert not audit.exists()


def adaptive_departure(before, carry, surface, bath, dt, **kwargs):
    for duration in (dt, dt / 4, dt / 16, dt / 64, dt / 256, dt / 1024):
        entry = {"before_time_s": before.time_s, "requested_dt_s": duration}
        try:
            result = original(before, carry, surface, bath, duration, **kwargs)
        except Exception as error:
            entry.update(accepted=False, error=str(error),
                         exact_rollback_identity=getattr(error, "rollback_state", before) is before)
            attempts.append(entry)
            audit.write_text(json.dumps(attempts, indent=2))
            print(json.dumps({"departure_retry": entry}), flush=True)
            if "tangent predictor left" not in str(error):
                raise
        else:
            entry.update(accepted=True, ledger=asdict(result[0].ledger),
                         native_acceptance_changed=False)
            attempts.append(entry)
            audit.write_text(json.dumps(attempts, indent=2))
            return result
    raise RuntimeError("All declared native departure durations refused")


entry_script = Path(__file__).with_name("run_b_with_birth_diagnostics.py")
policy = out.with_suffix(".recovery_policy.json")
assert not policy.exists()
policy.write_text(json.dumps({
    "contact_conductance_mol_m2_s": coefficient,
    "physical_law_changed": False,
    "native_acceptance_changed": False,
    "same_IVP_failed_birth_iterate_may_be_reused_as_guess": True,
    "other_case_accepted_material_imported": False,
    "departure_duration_factors": [1, .25, .0625, .015625, .00390625, .0009765625],
    "component_energy_tolerance": 1e-10,
    "manuscript_adoption": False,
    "source_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (Path(__file__), entry_script)},
}, indent=2))
engine.face.advance_departure = adaptive_departure
try:
    runpy.run_path(str(entry_script), run_name="__main__")
finally:
    engine.face.advance_departure = original
