"""Paper figures from accepted Core2 histories; no trajectory continuation."""

from pathlib import Path
import argparse
import csv
import hashlib
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return [
            {k: float(v) if v and k != "stage" and k != "source" else v for k, v in row.items()}
            for row in csv.DictReader(stream)
        ]


def measured(quantity):
    base = ROOT / "paper/analysis/datasets/faner2019_v2"
    with (base / "adjudication_per_point.csv").open(newline="", encoding="utf-8") as stream:
        bands = {r["point_id"]: float(r["uncertainty_abs"]) for r in csv.DictReader(stream)}
    with (base / "curves.csv").open(newline="", encoding="utf-8") as stream:
        return sorted(
            (float(r["time_s"]), float(r["value"]), bands[r["point_id"]])
            for r in csv.DictReader(stream)
            if r["figure"] == "2"
            and r["material"] == "soybean"
            and r["temperature_c"] == "120"
            and r["quantity"] == quantity
            and r["series"] == "experimental"
            and not r["superseded_by"]
        )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dt-history", type=Path, required=True)
    p.add_argument("--faner-history", type=Path, required=True)
    p.add_argument("--dt-cells", type=int, default=2)
    p.add_argument("--faner-cells", type=int, default=4)
    args = p.parse_args()
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "savefig.dpi": 200,
        }
    )
    dt, faner = rows(args.dt_history), rows(args.faner_history)
    if dt[-1]["time_s"] < 1766 or faner[-1]["t_s"] != 240:
        raise ValueError("complete accepted histories required")
    out = PAPER / "figures"
    fig, axes = plt.subplots(3, 1, figsize=(7, 6.5), sharex=True, constrained_layout=True)
    time = [r["time_s"] for r in dt]
    axes[0].plot(time, [r["hexane_kg_kg_dry"] for r in dt], color="#b85b12")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Hexane (kg/kg dry meal)")
    axes[1].plot(time, [r["particle_water_kg_kg_dry"] for r in dt], label="particle water")
    axes[1].plot(time, [r["external_water_kg_kg_dry"] for r in dt], label="exterior liquid water")
    axes[1].set_ylabel("Water (kg/kg dry meal)")
    axes[1].legend(fontsize=8)
    real = [r for r in dt if r["temperature_k"] != ""]
    axes[2].plot(
        [r["time_s"] for r in real],
        [r["temperature_k"] - 273.15 for r in real],
        label="volume-mean particle temperature",
    )
    axes[2].set_ylabel("Temperature (deg C)")
    axes[2].set_xlabel("Actual physical time (s)")
    inset = axes[0].inset_axes([0.55, 0.48, 0.4, 0.43])
    transition = [r for r in dt if 490 <= r["time_s"] <= 506]
    inset.plot([r["time_s"] for r in transition], [r["wet_volume_fraction"] for r in transition])
    inset.set(xlabel="time (s)", ylabel="wet volume fraction", xlim=(490, 506), ylim=(0, 1.05))
    inset.tick_params(labelsize=7)
    inset.xaxis.label.set_size(7)
    inset.yaxis.label.set_size(7)
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.savefig(out / "fig_dt_particle_history_core2.pdf")
    fig.savefig(out / "fig_dt_particle_history_core2.png")
    plt.close(fig)

    fig, axes = plt.subplots(4, 1, figsize=(7, 8), sharex=True, constrained_layout=True)
    time = [r["t_s"] for r in faner]
    axes[0].plot(
        time, [r["hexane_kg_kg_dry"] for r in faner], label="declared single sphere, active water"
    )
    axes[0].set_ylabel("Hexane (kg/kg dry meal)")
    axes[1].plot(
        time, [r["mean_temperature_k"] - 273.15 for r in faner], label="particle volume mean"
    )
    axes[1].set_ylabel("Temperature (deg C)")
    for ax, quantity, label in [
        (axes[0], "solvent_loading", "measured sample loading"),
        (axes[1], "particle_temperature", "measured holder thermocouple"),
    ]:
        data = measured(quantity)
        ax.errorbar(
            [r[0] for r in data],
            [r[1] for r in data],
            yerr=[r[2] for r in data],
            fmt="o",
            ms=3,
            capsize=2,
            color="#923d32",
            label=label,
        )
        ax.legend(fontsize=7)
    axes[2].plot(time, [r["water_kg_kg_dry"] for r in faner], label="particle mean")
    axes[2].plot(
        time,
        [
            r["surface_water_kg_kg_dry"] if r["surface_water_kg_kg_dry"] != "" else float("nan")
            for r in faner
        ],
        label="surface retained water",
    )
    axes[2].axhline(0.023, color=".6", ls="--", label="qualified lower moisture anchor")
    axes[2].legend(fontsize=7)
    axes[2].set_ylabel("Water (kg/kg dry meal)")
    axes[3].plot(time, [r["front_z"] for r in faner])
    axes[3].set_ylabel("Wet-core volume fraction")
    axes[3].set_xlabel("Actual physical time (s)")
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.savefig(out / "fig_faner_particle_history_core2.pdf")
    fig.savefig(out / "fig_faner_particle_history_core2.png")
    plt.close(fig)
    record = {
        "trajectory_interpolation_or_tail_generation": False,
        "numerical_resolution": {
            "dt_radial_cells": args.dt_cells,
            "faner_radial_cells": args.faner_cells,
        },
        "full_history_continuum_convergence_claimed": False,
        "producer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {
            str(x.relative_to(ROOT)): hashlib.sha256(x.read_bytes()).hexdigest()
            for x in [
                args.dt_history.resolve(),
                args.faner_history.resolve(),
                ROOT / "paper/analysis/datasets/faner2019_v2/curves.csv",
                ROOT / "paper/analysis/datasets/faner2019_v2/adjudication_per_point.csv",
            ]
        },
    }
    (out / "fig_core2_particle_baths_sources.json").write_text(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
