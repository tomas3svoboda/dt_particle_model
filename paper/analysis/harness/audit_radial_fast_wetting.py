"""Reconstruct fast-wetting study storage and exterior flux ledgers independently."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import pickle

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.props import water as wa


def audit(result_path):
    data = json.loads(result_path.read_text())
    parent_path = Path(data["parent_checkpoint"])
    if hashlib.sha256(parent_path.read_bytes()).hexdigest() != data["parent_sha256"]:
        raise ValueError("parent checkpoint hash changed")
    parent = pickle.loads(parent_path.read_bytes())
    state = parent["current"]
    if hasattr(state, "to_uniform_fully_dry_transport_state"):
        state = state.to_uniform_fully_dry_transport_state()
    grid = state.grid
    params = state.config.pore
    pressure = state.config.pressure_pa
    mass = params.rho_dm_p * math.fsum(grid.volumes)
    area = grid.areas[-1]

    def storage(temperatures, compositions, surface):
        cells = [
            cp.evaluate_equilibrium(t, pressure, y, params)
            for t, y in zip(temperatures, compositions)
        ]
        result = [
            math.fsum(v * getattr(c, field) for v, c in zip(grid.volumes, cells))
            for field in (
                "total_water_concentration_kg_m3",
                "total_hexane_concentration_kg_m3",
                "energy_density_j_m3",
            )
        ]
        film = area * surface["liquid_water_kg_m2"]
        result[0] += film
        result[2] += film * wa.state_Tp(surface["temperature_k"], pressure, "liquid").u_mass
        return result

    initial = storage(state.temperatures_k, state.y_hexane, data["initial"]["surface"])
    scales = [max(abs(v), 1e-300) for v in initial]
    scales[2] = max(scales[2], mass * params.cp_dry_meal)
    streams = [[], [], []]
    maximum = [0.0, 0.0, 0.0]
    checks = []
    previous = state.time_s
    for attempt in data["attempts"]:
        if not attempt["accepted"]:
            continue
        dt = attempt["time_s"] - previous
        if dt <= 0 or abs(dt - attempt["dt_s"]) > 1e-10:
            raise ValueError("accepted chronology or interval inconsistent")
        surface = attempt["surface"]
        after = storage(attempt["temperature_k"], attempt["y_hexane"], surface)
        for stream, field in zip(
            streams,
            (
                "external_water_flux_kg_m2_s",
                "external_hexane_flux_kg_m2_s",
                "external_energy_flux_w_m2",
            ),
        ):
            stream.append(dt * area * surface[field])
        residual = [
            abs(math.fsum((a, -i, *s))) / scale
            for a, i, s, scale in zip(after, initial, streams, scales)
        ]
        maximum = [max(a, b) for a, b in zip(maximum, residual)]
        checks.append({"time_s": attempt["time_s"], "normalized_W_H_U": residual})
        previous = attempt["time_s"]
    return {
        "passed": bool(checks) and max(maximum) <= 1e-10 and data["sources_unchanged_during_run"],
        "declared_reporting_tolerance": 1e-10,
        "maximum_normalized_W_H_U": maximum,
        "physical_error_bounds_per_kg_dry_W_H_U": [r * s / mass for r, s in zip(maximum, scales)],
        "accepted_steps": len(checks),
        "last_time_s": previous,
        "checks": checks,
        "storage_from_native_primitives_not_archived_ledger_deltas": True,
        "result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
        "auditor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError("refusing to overwrite audit evidence")
    result = audit(args.result)
    args.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "checks"}))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
