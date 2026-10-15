"""Reuse the existing DT birth/departure recovery at owner-selected F=100."""

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[4]
SNAPSHOT = REPO  # Public export: native particle import closure is at repository root.
sys.path.insert(0, str(SNAPSHOT / "src"))
sys.path.insert(0, str(SNAPSHOT / "paper/analysis/harness"))

import experiment_core2_finite_wetting_b as engine  # noqa: E402


def main():
    seed_path = None
    if "--birth-seed" in sys.argv:
        position = sys.argv.index("--birth-seed")
        seed_path = Path(sys.argv[position + 1])
        del sys.argv[position : position + 2]
    out = Path(sys.argv[sys.argv.index("--out") + 1])
    factor = float(sys.argv[sys.argv.index("--factor") + 1])
    contact = float(sys.argv[sys.argv.index("--contact-conductance") + 1])
    assert factor == 100 and contact == 0.32
    out.parent.mkdir(parents=True, exist_ok=True)
    declaration = json.loads(out.with_name("properties_" + out.name).read_text())
    physical_law = {k: v for k, v in declaration.items() if k != "source_sha256" and k != "argv"}
    seed = json.loads(seed_path.read_text()) if seed_path else None
    if seed:
        assert seed["shell_mobility_factor"] == factor
        assert seed["K_l_mol_m2_s"] == contact
        assert seed["temperature_property_variant"] == physical_law

    original_birth = engine.driver.advance_birth
    original_departure = engine.face.advance_departure
    attempts = []

    def observed_birth(*args, **kwargs):
        if seed:
            assert args[0].time_s == seed["fully_wet_before_time_s"]
            kwargs["root_seed"] = seed
        try:
            return original_birth(*args, **kwargs)
        except Exception as error:
            debug = getattr(error, "study_debug", None)
            if debug:
                payload = {
                    "K_l_mol_m2_s": contact,
                    "shell_mobility_factor": factor,
                    "temperature_property_variant": physical_law,
                    "fully_wet_before_time_s": args[0].time_s,
                    "study_debug": debug,
                    "error": str(error),
                    "failed_trial_advanced_physical_time": False,
                    "exact_fully_wet_rollback_identity_preserved": error.rollback_state is args[0],
                }
                path = out.with_suffix(".birth_debug.json")
                assert not path.exists()
                path.write_text(
                    json.dumps(
                        payload,
                        indent=2,
                        default=lambda v: v.tolist() if hasattr(v, "tolist") else float(v),
                    )
                )
            raise

    def adaptive_departure(before, carry, surface, bath, dt, **kwargs):
        for duration in (dt, dt / 4, dt / 16, dt / 64, dt / 256, dt / 1024):
            entry = {"before_time_s": before.time_s, "requested_dt_s": duration}
            try:
                result = original_departure(before, carry, surface, bath, duration, **kwargs)
            except Exception as error:
                entry.update(
                    accepted=False,
                    error=str(error),
                    exact_rollback_identity=getattr(error, "rollback_state", before) is before,
                )
                attempts.append(entry)
                out.with_suffix(".departure_attempts.json").write_text(
                    json.dumps(attempts, indent=2)
                )
                if "tangent predictor left" not in str(error):
                    raise
            else:
                entry.update(
                    accepted=True, ledger=asdict(result[0].ledger), native_acceptance_changed=False
                )
                attempts.append(entry)
                out.with_suffix(".departure_attempts.json").write_text(
                    json.dumps(attempts, indent=2)
                )
                return result
        raise RuntimeError("All existing native departure durations refused")

    policy = {
        "factor": factor,
        "contact_conductance": contact,
        "reused_policy": "2026-10-14 nominal run_b_032.py and run_b_with_birth_diagnostics.py",
        "root_seed": str(seed_path) if seed_path else None,
        "root_seed_sha256": (
            hashlib.sha256(seed_path.read_bytes()).hexdigest() if seed_path else None
        ),
        "same_IVP_failed_birth_iterate_only": True,
        "other_case_accepted_material_imported": False,
        "departure_duration_factors": [1, 0.25, 0.0625, 0.015625, 0.00390625, 0.0009765625],
        "native_acceptance_changed": False,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    out.with_suffix(".recovery_policy.json").write_text(json.dumps(policy, indent=2))
    engine.driver.advance_birth = observed_birth
    engine.face.advance_departure = adaptive_departure
    try:
        engine.main()
    finally:
        engine.driver.advance_birth = original_birth
        engine.face.advance_departure = original_departure


if __name__ == "__main__":
    main()
