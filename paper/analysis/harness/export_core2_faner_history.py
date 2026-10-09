"""Verify accepted lineage and independently integrate external component/energy flow."""

from dataclasses import replace
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import experiment_core2_faner as faner

ROOT = faner.REPO_ROOT


def load_chain(tip):
    chain = []
    seen = set()
    path = tip.resolve()
    while True:
        if path in seen:
            raise ValueError("cyclic result lineage")
        seen.add(path)
        result = json.loads(path.read_text(encoding="utf-8"))
        if not result["sources_unchanged_during_run"]:
            raise ValueError(f"source changed in {path}")
        chain.append((path, result))
        parent = result.get("accepted_parent")
        if not parent:
            break
        # the records store each parent as written on the authors' machine (absolute
        # or repository-relative, Windows separators); resolve it in this repository
        text = parent["result"].replace("\\", "/")
        path = ROOT / text[text.index("paper/analysis/") :]
        if hashlib.sha256(path.read_bytes()).hexdigest() != parent["result_sha256"]:
            raise ValueError("parent result identity changed")
    return list(reversed(chain))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tip", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--ledger-tolerance", type=float, default=1e-10)
    args = p.parse_args()
    chain = load_chain(args.tip)
    first = chain[0][1]
    inputs = SimpleNamespace(**first["inputs"])
    rows = []
    law = faner.core.sp.ContinuedPositiveLuikovParams()
    with (
        faner.declared_faner_film_multiplier(inputs.film_multiplier),
        faner.pure_bulk_domain_adapter(),
    ):
        feed = faner.FanerAttachedMarch(inputs)
        dry, area = feed.dry, feed.area
        for item in first["feed_rows"]:
            t, h, w = (
                item["temperature_k"],
                item["hexane_kg_kg_dry"],
                item["retained_water_kg_kg_dry"],
            )
            u = dry * faner.core.sphere.combined_high_loading_specific_energy(
                t, h, replace(feed.wet, X_water=w)
            )
            rates = None
            if item["beta"] is not None:
                reference = faner.core.ct.prepare_reduced_film_surface_state(
                    feed.boundary.film_boundary,
                    cell_temperature_k=t,
                    pressure_pa=faner.core.PRESSURE_PA,
                    particle_radius_m=inputs.radius,
                    outer_half_cell_distance_m=inputs.radius * 0.5,
                    particle_thermal_conductivity_w_m_k=feed.wet.conductivity,
                    binary_gas_interaction_k_wh=feed.boundary.film_boundary.binary_gas_interaction_k_wh,
                )
                film = faner.core.ct._reduced_film_coefficients_at_surface_temperature(
                    reference, feed.interface_temperature
                )
                capacities = faner.core.efr.phy053_partial_enthalpy_secant_heat_capacities(
                    surface_temperature_k=feed.interface_temperature,
                    bulk_temperature_k=inputs.bath_temperature,
                    pressure_pa=faner.core.PRESSURE_PA,
                    y_hexane=1.0,
                )
                mass = item["beta"] * film.heat_transfer_w_m2_k / capacities.hexane_j_kg_k
                rates = feed.rates(t, mass)
            rows.append(
                {
                    "t_s": item["t_s"],
                    "hexane_kg_kg_dry": h,
                    "water_kg_kg_dry": w,
                    "mean_temperature_k": t,
                    "surface_water_kg_kg_dry": None,
                    "interface_temperature_k": feed.interface_temperature,
                    "surface_temperature_k": feed.interface_temperature,
                    "front_z": 1.0,
                    "stage": "attached",
                    "energy_j": u,
                    "hexane_out_kg_s": rates["hexane_out_kg_s"] if rates else 0.0,
                    "water_out_kg_s": rates["water_out_kg_s"] if rates else 0.0,
                    "energy_out_w": rates["energy_out_w"] if rates else 0.0,
                    "source": chain[0][0].name,
                }
            )

    def append(after, surface, path, stage):
        if after["time_s"] <= rows[-1]["t_s"]:
            raise ValueError("accepted chronology is not strictly increasing")
        activity = 1 - surface["water_activity_trace_deficit"]
        if activity <= 0:
            raise ValueError("export cannot infer positive surface moisture from rounded activity")
        rows.append(
            {
                "t_s": after["time_s"],
                "hexane_kg_kg_dry": after["X_hexane_kg_kg"],
                "water_kg_kg_dry": after["X_water_kg_kg"],
                "mean_temperature_k": after["T_k"],
                "surface_water_kg_kg_dry": faner.core.sp.water_retained(activity, law),
                "interface_temperature_k": after["interface_temperature_k"],
                "surface_temperature_k": surface["temperature_k"],
                "front_z": after["front_z"],
                "stage": stage,
                "energy_j": after["energy_j"],
                "hexane_out_kg_s": area * surface["external_hexane_flux_kg_m2_s"],
                "water_out_kg_s": area * surface["external_water_flux_kg_m2_s"],
                "energy_out_w": area * surface["external_energy_flux_w_m2"],
                "source": path.name,
            }
        )

    append(
        first["last_accepted"],
        first["attempts"][0]["audit"]["surface_water_node"],
        chain[0][0],
        "birth",
    )
    for path, result in chain[1:]:
        if result.get("accepted"):
            append(result["last_accepted"], result["surface"], path, result["event"])
        else:
            for item in result["attempts"]:
                if item["accepted"]:
                    append(
                        item["after"], item["surface"], path, item.get("stage", "stationary dry")
                    )
    initial = [
        rows[0]["hexane_kg_kg_dry"] * dry,
        rows[0]["water_kg_kg_dry"] * dry,
        rows[0]["energy_j"],
    ]
    scales = [initial[0], initial[1], dry * 2e5]
    outgoing = [[], [], []]
    maximum = [0.0] * 3
    for old, new in zip(rows, rows[1:]):
        dt = new["t_s"] - old["t_s"]
        for stream, field in zip(outgoing, ["hexane_out_kg_s", "water_out_kg_s", "energy_out_w"]):
            stream.append(dt * new[field])
        inventories = [new["hexane_kg_kg_dry"] * dry, new["water_kg_kg_dry"] * dry, new["energy_j"]]
        residual = [math.fsum([v, -i, *s]) for v, i, s in zip(inventories, initial, outgoing)]
        for i, (r, scale) in enumerate(zip(residual, scales)):
            maximum[i] = max(maximum[i], abs(r) / scale)
    # Owner RULE 2: this independently accumulated reporting allowance does
    # not alter any native step, cumulative, or combined acceptance gate.
    # N4 measured water defect 2.746e-10 means <6.87e-12 kg/kg dry meal.
    if not 0 < args.ledger_tolerance <= 4e-10:
        raise ValueError("explicit whole-history numerical tolerance must be <=4e-10")
    if max(maximum) > args.ledger_tolerance:
        raise ValueError(f"whole-history external ledger refused: {maximum}")
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "accepted_history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig, axes = plt.subplots(4, 1, figsize=(8, 8), sharex=True, constrained_layout=True)
    times = [r["t_s"] for r in rows]
    axes[0].plot(times, [r["hexane_kg_kg_dry"] for r in rows], color="tab:orange")
    axes[0].set_ylabel("Hexane / dry meal\n(kg/kg)")
    axes[1].plot(times, [r["water_kg_kg_dry"] for r in rows], label="whole particle")
    axes[1].plot(
        times, [r["surface_water_kg_kg_dry"] for r in rows], label="surface retained water"
    )
    axes[1].axhline(0.023, color="grey", ls="--", label="qualified lower anchor")
    axes[1].legend(fontsize=8)
    axes[1].set_ylabel("Water / dry meal\n(kg/kg)")
    for key, label in [
        ("mean_temperature_k", "particle mean"),
        ("surface_temperature_k", "surface"),
    ]:
        axes[2].plot(times, [r[key] - 273.15 for r in rows], label=label)
    axes[2].legend(fontsize=8)
    axes[2].set_ylabel("Temperature (°C)")
    axes[3].plot(times, [r["front_z"] for r in rows])
    axes[3].set_ylabel("Wet-core volume\nfraction")
    axes[3].set_xlabel("Actual physical time (s)")
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.suptitle(
        "Faner / Core2: complete coarse prescribed-bath study\nActive water; continued Luikov law; N=2; db_max_dapp_max"
    )
    fig.savefig(args.out / "accepted_history.png", dpi=180)
    fig.savefig(args.out / "accepted_history.pdf")
    plt.close(fig)
    summary = {
        "complete_falling_rate_experiment": chain[-1][1].get(
            "complete_falling_rate_experiment", False
        ),
        "publication_ready": False,
        "last_accepted": rows[-1],
        "row_count": len(rows),
        "independent_whole_history_external_ledger_maximum_normalized_residuals": maximum,
        "declared_whole_history_ledger_tolerance": args.ledger_tolerance,
        "whole_history_numerical_error_bounds_per_kg_dry": {
            "hexane_kg_kg": maximum[0] * initial[0] / dry,
            "water_kg_kg": maximum[1] * initial[1] / dry,
            "energy_j_kg": maximum[2] * scales[2] / dry,
        },
        "source_result_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path, _ in chain
        },
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                k: summary[k]
                for k in [
                    "complete_falling_rate_experiment",
                    "row_count",
                    "independent_whole_history_external_ledger_maximum_normalized_residuals",
                ]
            }
        )
    )


if __name__ == "__main__":
    main()
