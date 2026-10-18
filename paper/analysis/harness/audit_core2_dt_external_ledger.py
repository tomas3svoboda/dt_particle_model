"""Independently integrate archived exterior rates and native state functions."""

from dataclasses import replace
from pathlib import Path
import argparse
import hashlib
import json
import math
import pickle

import experiment_core2_dt_feed as feed_driver
import experiment_core2_dt_birth as radial
import experiment_core2_dt_face as boundary_driver

core = radial.core
ROOT = radial.REPO_ROOT


def inventory(transport, config):
    """Evaluate storage from archived primitives, independently of ledger deltas."""
    values = [[], [], []]
    wet = iter(zip(transport["wet_temperatures_k"], transport["wet_retained_water_loadings"]))
    dry = iter(zip(transport["dry_temperatures_k"], transport["dry_y_hexane"]))
    dual = iter(
        transport.get("wet_retained_water_capacity_duals_over_rt", ())
        or [0.0] * len(transport["wet_temperatures_k"])
    )
    temperatures = []
    for i, cell in enumerate(transport["geometry"]["cells"]):
        volume = cell["wet_volume_m3"]
        if volume:
            t, w = next(wet)
            state = radial.ww.evaluate_cell(
                t,
                w,
                transport["historical_hexane_loadings"][i],
                transport["oil_fraction_labels"][i],
                config.wet,
                retained_water_capacity_dual_over_rt=next(dual),
            )
            densities = (
                state.retained_water_concentration_mol_m3,
                state.total_hexane_concentration_kg_m3 / core.hx.M,
                state.energy_density_j_m3,
            )
            for stream, density in zip(values, densities):
                stream.append(volume * density)
            temperatures.append(volume * t)
        volume = cell["dry_volume_m3"]
        if volume:
            t, y = next(dry)
            state = core.cp.evaluate_equilibrium(
                t,
                config.dry.pressure_pa,
                y,
                replace(config.dry.pore, w_o=transport["oil_fraction_labels"][i]),
            )
            densities = (
                state.total_water_concentration_mol_m3,
                state.total_hexane_concentration_mol_m3,
                state.energy_density_j_m3,
            )
            for stream, density in zip(values, densities):
                stream.append(volume * density)
            temperatures.append(volume * t)
    return [math.fsum(s) for s in values], math.fsum(temperatures) / math.fsum(
        c["wet_volume_m3"] + c["dry_volume_m3"] for c in transport["geometry"]["cells"]
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        raise ValueError("refusing to overwrite audit evidence")
    manifest = json.loads(args.manifest.read_text())
    base = ROOT / "paper/analysis/results_2026-10-06"
    sources = {}

    def load(name):
        path = base / "repairs" / (name + ".json")
        data = json.loads(path.read_text())
        if (
            data.get("source_changed_during_run")
            or data.get("sources_unchanged_during_run") is False
        ):
            raise ValueError("source provenance failed")
        sources[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        return data

    feed_path = base / "shared_model/dt_feed_massfilm_gated_R1p5.json"
    data = json.loads(feed_path.read_text())
    sources[str(feed_path.relative_to(ROOT))] = hashlib.sha256(feed_path.read_bytes()).hexdigest()
    birth = load(manifest["birth"])
    checkpoint_path = base / "repairs" / (manifest["dry"][-1] + ".accepted.pkl")
    checkpoint = pickle.loads(checkpoint_path.read_bytes())
    config = checkpoint["cut_config"]
    bath_path = core.DESTINY_PATH.with_name("destiny_v3_exit105.json")
    bath = json.loads(bath_path.read_text())
    march = feed_driver.FeedMarch(bath, data["radius_m"], 2)
    mass, area = march.dry, march.area
    initial_row = data["rows"][0]
    initial = [
        initial_row["hexane_kg_kg_dry"] * mass,
        initial_row["retained_water_kg_kg_dry"] * mass,
        initial_row["total_particle_plus_film_energy_j"],
    ]
    scales = [initial[0], initial[1], mass * 2e5]
    streams = [[], [], []]
    maximum = [0.0] * 3
    records = []
    clock = 0.0

    def accept(time, dt, actual, rates, source):
        nonlocal clock
        if abs(time - clock - dt) > 1e-9 or dt <= 0:
            raise ValueError("trajectory is not consecutive actual physical time")
        for stream, rate in zip(streams, rates):
            stream.append(dt * rate)
        residuals = [math.fsum([a, -i, *s]) for a, i, s in zip(actual, initial, streams)]
        normalized = [abs(r) / s for r, s in zip(residuals, scales)]
        maximum[:] = [max(a, b) for a, b in zip(maximum, normalized)]
        clock = time
        records.append({"time_s": time, "source": source, "normalized_residuals": normalized})

    for old, row in zip(data["rows"], data["rows"][1:]):
        march.time = old["t_s"]
        boundary = boundary_driver.bath_boundary(bath, march.radius, row["t_s"])
        if old["t_s"] < march.pd_end:
            native = core.reduced_film_boundary(
                temperature_k=345.15,
                y_hexane=0.6644429974143522,
                pair=march.pair,
                forcing_timescale_s=0.075,
                label="unchanged pre-desolventizer dome audit",
            )
            boundary = replace(
                boundary, temperature_k=345.15, y_hexane=0.6644429974143522, film_boundary=native
            )
        total, *_ = radial.film_at_beta(
            boundary, config, march.radius, march.ti, row["beta"], surface_y_hexane=march.yi
        )
        native, _, _, energy, _, _, active_area = march.interface(
            row["temperature_k"], total, row["t_s"]
        )
        actual = [
            row["hexane_kg_kg_dry"] * mass,
            (row["retained_water_kg_kg_dry"] + row["external_water_kg_kg_dry"]) * mass,
            row["total_particle_plus_film_energy_j"],
        ]
        rates = [
            active_area * native.hexane_molar_flux_mol_m2_s * core.hx.M,
            active_area * native.water_molar_flux_mol_m2_s * core.wa.M,
            energy,
        ]
        accept(row["t_s"], row["t_s"] - old["t_s"], actual, rates, "feed")

    def surface_actual(time, dt, inv, surface, source):
        film = area * surface["liquid_water_kg_m2"]
        actual = [
            inv[1] * core.hx.M,
            inv[0] * core.wa.M + film,
            inv[2]
            + film
            * core.wa.state_Tp(surface["temperature_k"], config.dry.pressure_pa, "liquid").u_mass,
        ]
        rates = [
            area * surface[f"external_{name}_flux_{unit}"]
            for name, unit in [("hexane", "kg_m2_s"), ("water", "kg_m2_s"), ("energy", "w_m2")]
        ]
        accept(time, dt, actual, rates, source)

    inv, _ = inventory(birth["after_transport"], config)
    surface_actual(
        birth["after_transport"]["time_s"],
        birth["dt_s"],
        inv,
        birth["surface_water_node"],
        manifest["birth"],
    )

    def radial_rows(data, name):
        for row in data["radial_continuation"]:
            if row["accepted"]:
                inv, _ = inventory(row["after_transport"], config)
                surface_actual(row["time_s"], row["ledger"]["dt_s"], inv, row["surface"], name)

    radial_rows(birth, manifest["birth"])
    for name in manifest["radial_before_face"]:
        radial_rows(load(name), name)
    for name in manifest["faces"]:
        data = load(name)
        after = data["after"].get("transport", data["after"])
        inv, _ = inventory(after, config)
        surface_actual(after["time_s"], after["time_s"] - clock, inv, data["surface"], name)
    radial_rows(load(manifest["inner"]), manifest["inner"])
    # Target files omit intermediate primitive states. Integrate their exterior
    # rates, then independently check the end checkpoint; do not reconstruct a
    # purported state by subtracting ledger changes.
    for name in manifest["targets"]:
        data = load(name)
        for row in data["attempts"]:
            if row["accepted"]:
                surface = row["surface"]
                for stream, component, unit in zip(
                    streams, ["hexane", "water", "energy"], ["kg_m2_s", "kg_m2_s", "w_m2"]
                ):
                    stream.append(
                        row["ledger"]["dt_s"] * area * surface[f"external_{component}_flux_{unit}"]
                    )
                clock = row["time_s"]
        end = pickle.loads((base / "repairs" / (name + ".accepted.pkl")).read_bytes())["current"]
        snap = core.ci.inventory_snapshot(end.transport)
        film = (
            area
            * data["attempts"][sum(r["accepted"] for r in data["attempts"]) - 1]["surface"][
                "liquid_water_kg_m2"
            ]
        )
        surface = data["attempts"][sum(r["accepted"] for r in data["attempts"]) - 1]["surface"]
        actual = [
            snap.total_hexane_mol * core.hx.M,
            snap.total_water_mol * core.wa.M + film,
            snap.total_energy_j
            + film
            * core.wa.state_Tp(surface["temperature_k"], config.dry.pressure_pa, "liquid").u_mass,
        ]
        residual = [
            abs(math.fsum([a, -i, *s])) / scale
            for a, i, s, scale in zip(actual, initial, streams, scales)
        ]
        maximum[:] = [max(a, b) for a, b in zip(maximum, residual)]
        records.append(
            {
                "time_s": clock,
                "source": name,
                "normalized_residuals": residual,
                "intermediate_target_state_functions_archived": False,
            }
        )
    data = load(manifest["projection"])
    after = data["after"]

    def fully_dry(temperatures, compositions):
        cells = [
            core.cp.evaluate_equilibrium(
                t, config.dry.pressure_pa, y, replace(config.dry.pore, w_o=o)
            )
            for t, y, o in zip(
                temperatures, compositions, [config.dry.pore.w_o] * checkpoint["current"].grid.n
            )
        ]
        return [
            math.fsum(
                v * getattr(c, field) for v, c in zip(checkpoint["current"].grid.volumes, cells)
            )
            for field in (
                "total_water_concentration_mol_m3",
                "total_hexane_concentration_mol_m3",
                "energy_density_j_m3",
            )
        ]

    surface_actual(
        after["time_s"],
        after["time_s"] - clock,
        fully_dry(after["temperatures_k"], after["y_hexane"]),
        data["surface"],
        manifest["projection"],
    )
    for name in manifest["dry"]:
        data = load(name)
        for row in data["attempts"]:
            if row["accepted"]:
                surface_actual(
                    row["time_s"],
                    row["dt_s"],
                    fully_dry(row["temperature_k"], row["y_hexane"]),
                    row["surface"],
                    name,
                )
    report = {
        "passed": max(maximum) <= 2e-10,
        "declared_tolerance": 2e-10,
        "maximum_normalized_external_ledger_residuals_H_W_U": maximum,
        "physical_error_bounds_per_kg_dry_H_W_U": [a * b / mass for a, b in zip(maximum, scales)],
        "last_time_s": clock,
        "checks": records,
        "source_result_sha256": sources,
        "state_functions_independent_of_archived_native_ledger_changes": True,
        "auditor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "bath_sha256": hashlib.sha256(bath_path.read_bytes()).hexdigest(),
        "end_checkpoint_sha256": hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
    }
    args.out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(
        json.dumps({k: v for k, v in report.items() if k not in ("checks", "source_result_sha256")})
    )


if __name__ == "__main__":
    main()
