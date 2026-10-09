"""Actual native stationary dry-particle continuation in the constant Faner bath."""

from pathlib import Path
from dataclasses import replace
import argparse
import hashlib
import json
import pickle
import traceback

import experiment_core2_faner as faner
import experiment_core2_dt_dry as dry_driver

FanerStudyControls = faner.FanerStudyControls
driver, core = faner.shared, faner.core


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--resume", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--horizon", type=float, default=240)
    p.add_argument("--dt", type=float, default=5)
    p.add_argument("--steps", type=int, default=200)
    p.add_argument("--wall-budget", type=float, default=180)
    args = p.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite evidence")
    parent_path = args.resume.with_name(args.resume.name.replace(".accepted.pkl", ".json"))
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    if not parent["sources_unchanged_during_run"]:
        raise ValueError("parent sources changed during run")
    checkpoint = pickle.loads(args.resume.read_bytes())
    before = checkpoint["before"]
    if hasattr(before, "to_uniform_fully_dry_transport_state"):
        before = before.to_uniform_fully_dry_transport_state()
    if not isinstance(before, dry_driver.ct.FullyDryTransportState):
        raise TypeError("native fully dry conversion required")
    before = replace(before, config=replace(before.config, nonlinear_residual_tolerance=1e-10))
    config = checkpoint["config"]
    boundary = checkpoint["boundary"]
    carry, surface = checkpoint["carry"], checkpoint["surface"]
    sources = [
        Path(__file__),
        Path(faner.__file__),
        Path(dry_driver.__file__),
        Path(driver.__file__),
        Path(driver.node.__file__),
        Path(driver.study.__file__),
        Path(driver.study.numerics.__file__),
        *sorted((driver.REPO_ROOT / "src/dtdc_simulator/core2").rglob("*.py")),
    ]

    def hashes():
        return {str(x.resolve()): hashlib.sha256(x.read_bytes()).hexdigest() for x in sources}

    for path, expected in parent["source_sha256_at_exit"].items():
        if str(driver.REPO_ROOT / "src/dtdc_simulator/core2") in path:
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
                raise ValueError("native source changed from parent")
    record = {
        "experiment": "native Faner stationary dry continuation",
        "attempts": [],
        "complete_falling_rate_experiment": False,
        "source_sha256_at_entry": hashes(),
        "accepted_parent": {
            "result": str(parent_path),
            "result_sha256": hashlib.sha256(parent_path.read_bytes()).hexdigest(),
            "checkpoint": str(args.resume),
            "checkpoint_sha256": hashlib.sha256(args.resume.read_bytes()).hexdigest(),
        },
        "last_accepted": parent["last_accepted"],
    }
    dt = args.dt
    try:
        with (
            faner.declared_faner_film_multiplier(boundary.film_boundary.stress_multiplier),
            faner.pure_bulk_domain_adapter(),
            driver.bounded_path_proof(3000000, args.wall_budget),
            driver.path_bounds.monotone_hexane_path_enclosure([]),
            driver.study.magnitude_storage_scale(),
        ):
            for index in range(args.steps):
                if before.time_s >= args.horizon:
                    break
                trial = min(dt, args.horizon - before.time_s)
                while trial >= 1e-5:
                    try:
                        after, new_surface, new_carry, report = dry_driver.advance(
                            before, config, surface, carry, trial, boundary
                        )
                        dry = config.wet.wet.rho_dm_p * sum(before.grid.volumes)
                        summary = {
                            "time_s": after.time_s,
                            "front_z": 0.0,
                            "X_hexane_kg_kg": report["hexane_mol"] * core.hx.M / dry,
                            "X_water_kg_kg": report["water_mol"] * core.wa.M / dry,
                            "T_k": sum(
                                t * v for t, v in zip(after.temperatures_k, after.grid.volumes)
                            )
                            / sum(after.grid.volumes),
                            "interface_temperature_k": None,
                            "hexane_mol": report["hexane_mol"],
                            "water_mol": report["water_mol"],
                            "energy_j": report["energy_j"],
                        }
                        record["attempts"].append({"accepted": True, "after": summary, **report})
                        before, surface, carry = after, new_surface, new_carry
                        record["last_accepted"] = summary
                        checkpoint.update(
                            before=before, surface=surface, carry=carry, source_record=str(args.out)
                        )
                        args.out.with_suffix(".accepted.pkl").write_bytes(pickle.dumps(checkpoint))
                        args.out.write_text(json.dumps(record, indent=2), encoding="utf-8")
                        dt = min(args.dt, trial * 1.5)
                        break
                    except Exception as error:
                        record["attempts"].append(
                            {
                                "accepted": False,
                                "at_s": before.time_s,
                                "dt_s": trial,
                                "error": str(error),
                                "optimizer_debug": getattr(error, "optimizer_debug", None),
                            }
                        )
                        if "computational budget exhausted" in str(error):
                            raise
                        trial *= 0.5
                else:
                    raise RuntimeError("minimum interval reached with no state advanced")
        record["complete_falling_rate_experiment"] = before.time_s >= args.horizon >= 240
    except Exception as error:
        record.update(error=str(error), traceback=traceback.format_exc())
    record["source_sha256_at_exit"] = hashes()
    record["sources_unchanged_during_run"] = (
        record["source_sha256_at_entry"] == record["source_sha256_at_exit"]
    )
    args.out.write_text(json.dumps(record, indent=2, allow_nan=True), encoding="utf-8")
    print(
        json.dumps(
            {
                k: record.get(k)
                for k in ["complete_falling_rate_experiment", "last_accepted", "error"]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
