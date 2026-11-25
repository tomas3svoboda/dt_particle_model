"""Review-only replacement of the original three-panel DT particle history.

Uses only accepted native history rows. No trajectory interpolation, smoothing,
extrapolation or manuscript/figure-source replacement. Bath interpolation is
the declared prescribed environment, not a manufactured particle trajectory.
"""

# ruff: noqa: E402
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bath", type=Path, required=True)
    args = parser.parse_args()
    root, out = args.case.resolve(), args.out.resolve()
    case = load(root / "case/case.json")
    production = load(root / "production.json")
    audit = load(root / "independent_whole_history_audit.json")
    if not all((case["complete_history"], production["independent_whole_history_passed"],
                production["extension_sources_unchanged"], audit["passed"], audit["complete_history"])):
        raise ValueError("Full accepted and independently audited history is required")
    out.mkdir(exist_ok=False)
    feed, bath = load(root / "case/feed.json"), load(args.bath.resolve())
    if hashlib.sha256(args.bath.read_bytes()).hexdigest() != case["bath_sha256"]:
        raise ValueError("Prescribed bath identity differs from the completed case")
    used = [
        root / p
        for p in (
            "case/case.json",
            "production.json",
            "independent_whole_history_audit.json",
            "case/feed.json",

        )
    ]
    used.append(args.bath.resolve())
    pd_end = feed["pd_source"]["reference"]["residence_s"]
    rows = []

    def add(t, temp, h, retained, external, z, regime, bath_t=None, bath_h=None):
        if t == 0:
            # The producer's initial-state metadata precedes the PD override.
            # Correct only the prescribed environment label, never the state.
            bath_t = feed["pd_source"]["temperature_k"]
            bath_h = feed["pd_source"]["y_hexane"]
        w = retained + external
        denominator = 1 + w + h
        if min(h, retained, external, z) < 0 or z > 1:
            raise ValueError("Inadmissible accepted plot state")
        rows.append(
            dict(
                time_s=t,
                regime=regime,
                temperature_c=temp - 273.15,
                hexane_kg_kg_dry=h,
                retained_water_kg_kg_dry=retained,
                external_water_kg_kg_dry=external,
                hexane_ppm_wet_basis=1e6 * h / denominator,
                moisture_pct_wet_basis=100 * w / denominator,
                wet_core_volume_fraction=z,
                bath_temperature_c=None if bath_t is None else bath_t - 273.15,
                bath_y_hexane=bath_h,
            )
        )

    for r in feed["rows"]:
        add(
            r["t_s"],
            r["temperature_k"],
            r["hexane_kg_kg_dry"],
            r["retained_water_kg_kg_dry"],
            r["external_water_kg_kg_dry"],
            1.0,
            "A",
            r["bath_temperature_k"],
            r["bath_y_hexane"],
        )
    for name, regime in (("front.history.jsonl", "B"), ("dry.history.jsonl", "C")):
        path = root / "case" / name
        used.append(path)
        with path.open(encoding="utf-8") as stream:
            for line in stream:
                r = json.loads(line)
                if r.get("kind") == "refusal":
                    continue
                if r.get("kind") == "initial" and any(
                    old["time_s"] == r["time_s"] and old["regime"] == regime for old in rows
                ):
                    continue  # Exact restart datum is not another physical advance.
                add(
                    r["time_s"],
                    r["temperature_k"],
                    r["hexane_kg_kg_dry"],
                    r["particle_water_kg_kg_dry"],
                    r["external_water_kg_kg_dry"],
                    r["front_z"] if regime == "B" else 0.0,
                    regime,
                )
    by_time = {}
    for r in rows:
        if r["time_s"] in by_time and by_time[r["time_s"]] != r:
            raise ValueError("different physical plot states share an accepted clock")
        by_time[r["time_s"]] = r
    rows = [by_time[t] for t in sorted(by_time)]
    if abs(rows[-1]["time_s"] - 1766.5355008060146) > 1e-8:
        raise ValueError("History does not reach prescribed DT/toaster outlet")
    last = rows[-1]
    result = {"outlet": dict(hexane_ppm_wet_basis=last["hexane_ppm_wet_basis"],
        moisture_pct_wet_basis=last["moisture_pct_wet_basis"], temperature_c=last["temperature_c"])}
    waypoints = bath["waypoints"]
    bt = np.array([r["t_s"] for r in waypoints])
    btemp = np.array([r["T_bath_K"] - 273.15 for r in waypoints])
    bh = np.array([r["y_hexane"] for r in waypoints])
    for r in rows:
        if r["bath_temperature_c"] is None:
            r["bath_temperature_c"] = float(np.interp(r["time_s"], bt, btemp))
            r["bath_y_hexane"] = float(np.interp(r["time_s"], bt, bh))
    with (out / "accepted_particle_history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    t = np.array([r["time_s"] for r in rows])
    minutes = t / 60

    def field(name):
        return np.array([r[name] for r in rows])

    birth = min(r["time_s"] for r in rows if r["regime"] == "B")
    core_end = min(r["time_s"] for r in rows if r["regime"] == "C")
    depletion = min(r["time_s"] for r in rows if r["regime"] == "C" and r["external_water_kg_kg_dry"] == 0.)
    fig, axes = plt.subplots(3, 1, figsize=(9.2, 10.7), sharex=True)
    purple, blue, orange, grey = "#89528d", "#148eae", "#bf6030", "0.4"
    for ax in axes:
        ax.axvspan(0, birth / 60, color="0.965", zorder=-2)
        ax.axvspan(birth / 60, core_end / 60, color="#e5eef3", zorder=-2)
        ax.axvspan(core_end / 60, depletion / 60, color="#f5f4e9", zorder=-2)
        for event in (pd_end, birth, core_end, depletion):
            ax.axvline(event / 60, color="0.6", linestyle="--" if event in (pd_end, depletion) else ":", lw=0.8)
        ax.set_xlim(0, 30.8)
    axes[0].plot(
        minutes, field("hexane_ppm_wet_basis"), color=purple, lw=1.8, label="Particle model"
    )
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Hexane (ppm wet meal)")
    axes[0].set_ylim(min(80, 0.8 * field("hexane_ppm_wet_basis").min()), 4e5)
    axes[1].plot(minutes, field("moisture_pct_wet_basis"), color=blue, lw=1.8)
    axes[1].set_ylabel("Particle moisture (% wet basis)")
    axes[2].plot(minutes, field("temperature_c"), color=orange, lw=1.8)
    axes[2].set_ylabel("Mean particle temperature (°C)")
    axes[2].set_xlabel("Time from feed (min)")
    axes[2].set_ylim(45, 115)
    bath_fields = (field("bath_y_hexane"), 1 - field("bath_y_hexane"), field("bath_temperature_c"))
    labels = ("Bath hexane mole fraction", "Bath water mole fraction", "Bath temperature (°C)")
    for i, (ax, values, label) in enumerate(zip(axes, bath_fields, labels)):
        twin = ax.twinx()
        twin.plot(minutes, values, color=grey, ls="--", lw=1.1, label="Prescribed bath")
        twin.set_ylabel(label, color=grey)
        twin.tick_params(axis="y", colors=grey)
        if i < 2:
            twin.set_ylim(0, 1)
        else:
            twin.set_ylim(45, 115)
        if i == 0:
            handles = ax.get_legend_handles_labels()[0] + twin.get_legend_handles_labels()[0]
            ax.legend(handles=handles, loc="upper right", frameon=False)
    for ax, bounds, label in zip(
        axes, ((100, 500), (16, 22), (105, 110)), ("100–500 ppm", "16–22%", "105–110°C")
    ):
        x = 30.1
        ax.errorbar(
            x, sum(bounds) / 2, yerr=(bounds[1] - bounds[0]) / 2, color="#2e7d32", capsize=6, lw=2
        )
        ax.text(29.7, bounds[1], label, ha="right", va="bottom", color="#2e7d32", fontsize=8)
    inset = axes[0].inset_axes([0.45, 0.34, 0.33, 0.26])
    selected = (t >= pd_end) & (t <= core_end + 2)
    inset.plot(t[selected] - pd_end, field("wet_core_volume_fraction")[selected], color=purple)
    inset.set_title("Liquid-hexane core", fontsize=8)
    inset.set_xlabel("After PD exit (s)", fontsize=7)
    inset.set_ylabel("V/V", fontsize=7)
    inset.tick_params(labelsize=7)
    inset.set_ylim(0, 1.05)
    for event in (birth, core_end):
        inset.axvline(event - pd_end, ls=":", color=grey, lw=0.7)
    inset_w = axes[1].inset_axes([0.50, 0.18, 0.32, 0.26])
    positive = t[field("external_water_kg_kg_dry") > 1e-6]
    right = min(t[-1], positive[-1] + 20) if len(positive) else pd_end + 60
    selected = (t >= pd_end - 15) & (t <= right)
    inset_w.plot(minutes[selected], field("external_water_kg_kg_dry")[selected], color=blue)
    inset_w.set_title("External water film", fontsize=8)
    inset_w.set_xlabel("Time from feed (min)", fontsize=7)
    inset_w.set_ylabel("kg/kg dry", fontsize=7)
    inset_w.tick_params(labelsize=7)
    inset_w.axvline(depletion / 60, color=grey, ls="--", lw=.8)
    for i, ax in enumerate(axes):
        ax.text(0.02, 0.87, f"({chr(97+i)})", transform=ax.transAxes, fontweight="bold")
    for a, b, label in ((0, birth, "Regime A"), (birth, core_end, "B"), (core_end, depletion, "C1"), (depletion, t[-1], "Regime C2")):
        axes[0].text(
            (a + b) / 120,
            1.02,
            label,
            transform=axes[0].get_xaxis_transform(),
            ha="center",
            fontsize=10,
        )
    fig.tight_layout(rect=(0, .01, 1, .98), h_pad=1.5)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(out / f"particle_history_review.{extension}", dpi=190)
    used.append(Path(__file__).resolve())
    pd_rows = [r for r in rows if r["time_s"] <= pd_end]
    manifest = dict(
        status="complete_review_candidate",
        outlet=result["outlet"],
        birth_time_s=birth,
        core_extinction_time_s=core_end,
        pd_external_water_max_kg_kg_dry=max(r["external_water_kg_kg_dry"] for r in pd_rows),
        pd_external_water_exit_kg_kg_dry=pd_rows[-1]["external_water_kg_kg_dry"],
        independent_audit=audit,
        accepted_plot_rows=len(rows),
        particle_trajectory_interpolated=False,
        fresh_original_feed=production["fresh_original_feed"],
        single_runtime_source_generation=production["single_runtime_source_generation"],
        manuscript_figure_replaced=False,
        inputs_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in used},
    )
    (out / "review_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {k: v for k, v in manifest.items() if k not in ("inputs_sha256", "independent_audit")}
        )
    )


if __name__ == "__main__":
    main()
