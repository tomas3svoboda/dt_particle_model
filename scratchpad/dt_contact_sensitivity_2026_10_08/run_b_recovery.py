"""Resume the same DT IVP with adaptive native face-departure steps only."""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / "src"))
sys.path.insert(0, str(root / "paper/analysis/harness"))
import experiment_core2_finite_wetting_b as engine

out_path = Path(sys.argv[sys.argv.index("--out") + 1])
audit_path = out_path.with_suffix(".departure_attempts.json")
assert not audit_path.exists()
original = engine.face.advance_departure
attempts = []

def adaptive_departure(before, carry, surface, bath, dt, **kwargs):
    for attempt_dt in (dt, dt / 4, dt / 16, dt / 64, dt / 256, dt / 1024):
        try:
            result = original(before, carry, surface, bath, attempt_dt, **kwargs)
        except Exception as error:
            attempts.append({"before_time_s": before.time_s, "requested_dt_s": attempt_dt,
                             "accepted": False, "error": str(error),
                             "before_state_identity": id(before),
                             "rollback_identity_if_exposed": getattr(error, "rollback_state", before) is before})
            audit_path.write_text(json.dumps(attempts, indent=2))
            print(json.dumps({"departure_retry": attempts[-1]}), flush=True)
            if "tangent predictor left" not in str(error):
                raise
        else:
            attempts.append({"before_time_s": before.time_s, "requested_dt_s": attempt_dt,
                             "accepted": True, "ledger": asdict(result[0].ledger),
                             "native_acceptance_changed": False})
            audit_path.write_text(json.dumps(attempts, indent=2))
            return result
    raise RuntimeError("All declared face-departure predictor steps refused; state unchanged")

engine.face.advance_departure = adaptive_departure
provenance = out_path.with_suffix(".recovery_policy.json")
assert not provenance.exists()
provenance.write_text(json.dumps({
    "physical_law_changed": False,
    "native_acceptance_changed": False,
    "inventory_fraction_arithmetic": "D9 drained mass / requested mass; no clipping",
    "face_departure_step_factors": [1, 0.25, 0.0625, 0.015625, 0.00390625, 0.0009765625],
    "ordinary_dt_from_arguments": True,
    "wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}, indent=2))
try:
    engine.main()
finally:
    engine.face.advance_departure = original
