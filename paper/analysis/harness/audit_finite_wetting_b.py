"""Re-evaluate Regime B inventories from archived native primitives.

Does not use solver ledger deltas or the march's reported inventory rows.
Trusted local accepted pickle supplies native constitutive parameter types.
"""

import argparse
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import pickle

import audit_core2_dt_external_ledger as radial_audit
import experiment_core2_dt_birth as driver

core = driver.core


def audit(path, parts=None, dry_path=None):
    data = json.loads(path.read_text())
    if parts:
        combined_history = []
        predecessor = None
        parameters = None
        first = None
        stable = True
        physical_sources = None
        for part in parts:
            segment = json.loads(part.read_text())
            declared = [segment["shell_mobility_factor"], segment["contact_conductance_mol_m2_s"]]
            native_and_bath = {
                name: sha
                for name, sha in segment["source_sha256_at_entry"].items()
                if "core2" in Path(name).parts
                or Path(name).name
                in ("destiny_v3_exit105.json", "dt_feed_massfilm_gated_R1p5.json")
            }
            if first is None:
                first, parameters = segment, declared
                physical_sources = native_and_bath
                if segment.get("parent_checkpoint"):
                    raise ValueError("full B chain must begin at native birth")
            elif (
                declared != parameters
                or segment["parent_sha256"]
                != hashlib.sha256(predecessor.with_suffix(".accepted.pkl").read_bytes()).hexdigest()
            ):
                raise ValueError("linked finite B inputs changed")
            if native_and_bath != physical_sources:
                raise ValueError(
                    "native physics, prescribed bath or original feed changed between B links"
                )
            for row in segment["history"]:
                if row["kind"] == "resumed":
                    if row["time_s"] != combined_history[-1]["time_s"]:
                        raise ValueError("chain chronology changed")
                    continue
                combined_history.append(row)
            stable = stable and segment["sources_unchanged_during_run"]
            predecessor = part
        if predecessor.resolve() != path.resolve():
            raise ValueError("last chain segment must be the audited result")
        data = {
            **data,
            "history": combined_history,
            "initial_time_s": first["initial_time_s"],
            "sources_unchanged_during_run": stable,
        }
        data.pop("parent_checkpoint", None)
    if dry_path:
        tail = json.loads(dry_path.read_text())
        if (
            tail["parent_sha256"]
            != hashlib.sha256(path.with_suffix(".accepted.pkl").read_bytes()).hexdigest()
        ):
            raise ValueError("B/C checkpoint identity changed")
        if (
            tail["additional_parallel_path_total_to_base_mobility_factor"]
            != data["shell_mobility_factor"]
            or tail["finite_liquid_contact_conductance_mol_m2_s"]
            != data["contact_conductance_mol_m2_s"]
        ):
            raise ValueError("B/C finite physical law changed")
        if not tail["bath_exit_reached"] or not tail["sources_unchanged_during_run"]:
            raise ValueError("complete source-stable C required")
        for name, expected in data["source_sha256_at_entry"].items():
            if "core2" in Path(name).parts or Path(name).name == "destiny_v3_exit105.json":
                if tail["source_sha256_at_entry"].get(name) != expected:
                    raise ValueError("native physics or prescribed bath changed at B/C")
        if tail["initial"]["time_s"] != data["history"][-1]["time_s"]:
            raise ValueError("B/C chronology changed")
        for attempt in tail["attempts"]:
            if not attempt.get("accepted"):
                continue
            data["history"].append(
                {
                    "kind": "dry_continuation",
                    "time_s": attempt["time_s"],
                    "surface": attempt["surface"],
                    "after_transport": {
                        "temperatures_k": attempt["temperature_k"],
                        "y_hexane": attempt["y_hexane"],
                    },
                }
            )
    checkpoint = pickle.loads(path.with_suffix(".accepted.pkl").read_bytes())
    config = checkpoint["cut_config"]
    feed_path = (
        driver.REPO_ROOT
        / "paper/analysis/results_2026-10-06/shared_model/dt_feed_massfilm_gated_R1p5.json"
    )
    before, _, _, carry = driver.initialize_from_feed(
        json.loads(feed_path.read_text()), coefficient_case="db_max_dapp_max", continued_water=True
    )
    grid = before.grid
    streams = [[], [], []]
    for t, w, v, h, oil in zip(
        before.temperatures_k,
        before.retained_water_loadings,
        grid.volumes,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
    ):
        c = driver.ww.evaluate_cell(t, w, h, oil, config.wet)
        for s, rho in zip(
            streams,
            (
                c.retained_water_concentration_mol_m3 * core.wa.M,
                c.total_hexane_concentration_kg_m3,
                c.energy_density_j_m3,
            ),
        ):
            s.append(v * rho)
    initial = [math.fsum(s) for s in streams]
    initial[0] += carry["external_water_mass_kg"]
    initial[2] += carry["external_water_energy_j"]
    scales = [abs(initial[0]), abs(initial[1]), carry["import_capacity_scale_j"]]

    def inventory(row):
        t = row["after_transport"]
        if "wet_temperatures_k" in t:
            values, _ = radial_audit.inventory(t, config)
            values[0] *= core.wa.M
            values[1] *= core.hx.M
        else:
            cells = [
                core.cp.evaluate_equilibrium(
                    temp, config.dry.pressure_pa, y, replace(config.dry.pore, w_o=oil)
                )
                for temp, y, oil in zip(
                    t["temperatures_k"], t["y_hexane"], [config.dry.pore.w_o] * grid.n
                )
            ]
            values = [
                math.fsum(v * getattr(c, key) for v, c in zip(grid.volumes, cells))
                for key in (
                    "total_water_concentration_kg_m3",
                    "total_hexane_concentration_kg_m3",
                    "energy_density_j_m3",
                )
            ]
        surface = row["surface"]
        film = grid.areas[-1] * surface["liquid_water_kg_m2"]
        values[0] += film
        values[2] += (
            film
            * core.wa.state_Tp(surface["temperature_k"], config.dry.pressure_pa, "liquid").u_mass
        )
        return values

    if data.get("parent_checkpoint"):
        parent_path = Path(data["parent_checkpoint"])
        if hashlib.sha256(parent_path.read_bytes()).hexdigest() != data["parent_sha256"]:
            raise ValueError("parent bytes changed")
        initial = inventory(data["history"][0])
    clock = data["initial_time_s"]
    streams = [[], [], []]
    maxima = [0.0, 0.0, 0.0]
    checked = []
    for row in data["history"]:
        if row["kind"] == "resumed":
            continue
        dt = row["time_s"] - clock
        if dt <= 0:
            raise ValueError("accepted chronology refused")
        values = inventory(row)
        for s, key in zip(
            streams,
            (
                "external_water_flux_kg_m2_s",
                "external_hexane_flux_kg_m2_s",
                "external_energy_flux_w_m2",
            ),
        ):
            s.append(dt * grid.areas[-1] * row["surface"][key])
        residual = [
            abs(math.fsum((a, -b, *s))) / scale
            for a, b, s, scale in zip(values, initial, streams, scales)
        ]
        maxima = [max(a, b) for a, b in zip(maxima, residual)]
        checked.append({"time_s": row["time_s"], "kind": row["kind"], "normalized_W_H_U": residual})
        clock = row["time_s"]
    return {
        "passed": bool(checked) and max(maxima) <= 1e-10 and data["sources_unchanged_during_run"],
        "regime_b_completed": data["status"] == "REGIME_B_COMPLETED",
        "bath_exit_reached": bool(dry_path),
        "dry_result_sha256": (
            hashlib.sha256(dry_path.read_bytes()).hexdigest() if dry_path else None
        ),
        "auditor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "declared_tolerance": 1e-10,
        "maximum_normalized_W_H_U": maxima,
        "physical_error_bounds_per_kg_dry_W_H_U": [
            r * s / (config.dry.pore.rho_dm_p * math.fsum(grid.volumes))
            for r, s in zip(maxima, scales)
        ],
        "accepted_steps": len(checked),
        "events": [r["kind"] for r in checked if r["kind"] not in ("partial", "positive_target")],
        "result_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "checkpoint_sha256": hashlib.sha256(
            path.with_suffix(".accepted.pkl").read_bytes()
        ).hexdigest(),
        "linked_results": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (parts or [path])
        },
        "checks": checked,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--parts", type=Path, nargs="+")
    parser.add_argument("--dry", type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError("refusing to overwrite evidence")
    report = audit(args.result, args.parts, args.dry)
    args.out.write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "checks"}))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
