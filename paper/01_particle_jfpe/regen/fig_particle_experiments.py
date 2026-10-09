"""Publication figures from complete particle histories and source readings."""

import argparse
import bisect
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fig_core2_particle_baths import measured, rows

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
BASE = ROOT / "paper/analysis/results_2026-10-06/repairs"
FINITE = ROOT / "paper/analysis/results_2026-10-07/finite_wetting"
CONTACT_032 = ROOT / "paper/analysis/results_2026-10-08/dt_contact_0p32"
# Completed finite-contact DT marches that this figure can draw. Each case names
# the accepted regime-B segments in checkpoint order, the regime-C result, the
# audits that bind them, and the coefficients every file must declare.
CASES = {
    "k0p1": {
        "segments": [FINITE / f"b_f1000_k0p1_dt0p2_v{i}.json" for i in range(2, 6)],
        "dry": FINITE / "c_f1000_k0p1_dt5.json",
        "through_audit": FINITE / "b_through_exit_dt5.audit.json",
        "b_audit": FINITE / "b_complete.audit.json",
        "summary": FINITE / "completion_summary.json",
        "K": 0.1,
        "factor": 1000,
    },
    "k0p32": {
        "segments": [CONTACT_032 / "b_02.json", CONTACT_032 / "r3_b_01.json"],
        "dry": CONTACT_032 / "c_dt5.json",
        "through_audit": CONTACT_032 / "through_exit.audit.json",
        "b_audit": CONTACT_032 / "b.audit.json",
        "summary": None,
        "K": 0.32,
        "factor": 1000,
    },
}
COLORS = {"hexane": "#885184", "water": "#1689b0", "temperature": "#b65e32"}
OUTLET_COLOR = "#2e7d32"
PD_END_S = 490.20224190324825
BATH_PATH = ROOT / "paper/analysis/results_2026-10-01/dtdc_particle_destiny/destiny_v3_exit105.json"
FEED_PATH = ROOT / "paper/analysis/results_2026-10-06/shared_model/dt_feed_massfilm_gated_R1p5.json"


def prescribed_bath(exit_s):
    """Reconstruct the consumed boundary, including the feed adapter's PD override."""
    bath = json.loads(BATH_PATH.read_text(encoding="utf-8"))
    waypoints = bath["waypoints"]
    clocks = [r["t_s"] for r in waypoints]

    def interpolate(clock):
        # Same piecewise-linear interpolation as dtdc_particle_destiny_march;
        # all source fractions lie strictly inside its protective bounds.
        i = min(max(bisect.bisect_left(clocks, clock), 1), len(clocks) - 1)
        left, right = waypoints[i - 1], waypoints[i]
        w = (clock - left["t_s"]) / (right["t_s"] - left["t_s"])
        temperature = left["T_bath_K"] + w * (right["T_bath_K"] - left["T_bath_K"])
        hexane = left["y_hexane"] + w * (right["y_hexane"] - left["y_hexane"])
        assert 1e-6 < hexane < 1 - 1e-6
        return temperature, hexane

    assert exit_s == clocks[-1]
    pd = (345.15, 0.6644429974143522)
    feed = json.loads(FEED_PATH.read_text(encoding="utf-8"))
    checks = 0
    max_temperature_error = max_fraction_error = 0.0
    for row in feed["rows"]:
        if row["bath_temperature_k"] is None:
            continue
        temperature, hexane = pd if row["t_s"] <= PD_END_S else interpolate(row["t_s"])
        max_temperature_error = max(
            max_temperature_error, abs(temperature - row["bath_temperature_k"])
        )
        max_fraction_error = max(max_fraction_error, abs(hexane - row["bath_y_hexane"]))
        checks += 1
    assert max_temperature_error < 1e-10 and max_fraction_error < 1e-12
    # Duplicate PD-exit clock represents the imposed change of contacted bath.
    # Include every later knot, including the sub-second steam-contact peak.
    history = [(0.0, *pd), (PD_END_S, *pd), (PD_END_S, *interpolate(PD_END_S))]
    history.extend(
        (r["t_s"], r["T_bath_K"], r["y_hexane"]) for r in waypoints if PD_END_S < r["t_s"] <= exit_s
    )
    audit = {
        "archived_feed_boundaries_checked": checks,
        "maximum_temperature_difference_k": max_temperature_error,
        "maximum_hexane_mole_fraction_difference": max_fraction_error,
        "pd_temperature_k": pd[0],
        "pd_hexane_mole_fraction": pd[1],
        "post_pd_interpolation": "piecewise linear through original bath knots",
        "water_mole_fraction": "1 - hexane mole fraction (binary bath)",
        "history": [
            {"time_s": t, "temperature_k": temp, "y_hexane": h, "y_water": 1 - h}
            for t, temp, h in history
        ],
    }
    return history, audit


def read_finite(case):
    paths = list(case["segments"])
    data = []
    parent = None
    for path in paths:
        result = json.loads(path.read_text())
        assert result["sources_unchanged_during_run"]
        assert result["shell_mobility_factor"] == case["factor"]
        assert result["contact_conductance_mol_m2_s"] == case["K"]
        if parent:
            # The accepted checkpoint a segment resumes from must be the state the
            # previous accepted segment left; recovery attempts that refused in
            # between re-saved that same state byte for byte.
            assert (
                result["parent_sha256"]
                == hashlib.sha256(parent.with_suffix(".accepted.pkl").read_bytes()).hexdigest()
            )
        for r in result["history"]:
            if r["kind"] == "resumed":
                assert r["time_s"] == data[-1]["time_s"]
                continue
            data.append(r)
        parent = path
    assert result["status"] == "REGIME_B_COMPLETED"
    dry_path = case["dry"]
    dry = json.loads(dry_path.read_text())
    assert dry["bath_exit_reached"] and dry["sources_unchanged_during_run"]
    assert (
        dry["parent_sha256"]
        == hashlib.sha256(parent.with_suffix(".accepted.pkl").read_bytes()).hexdigest()
    )
    assert dry["finite_liquid_contact_conductance_mol_m2_s"] == case["K"]
    assert dry["additional_parallel_path_total_to_base_mobility_factor"] == case["factor"]
    data.extend({**r, "front_z": 0.0} for r in dry["history"][1:])
    audit_path = case["through_audit"]
    audit = json.loads(audit_path.read_text())
    assert audit["passed"]
    assert audit["dry_result_sha256"] == hashlib.sha256(dry_path.read_bytes()).hexdigest()
    b_audit = json.loads(case["b_audit"].read_text())
    assert b_audit["passed"]
    for path in paths:
        # the audit keys are Windows-style relative paths; look them up the same way
        key = str(path.relative_to(ROOT)).replace("/", "\\")
        expected = b_audit["linked_results"][key]
        assert expected == hashlib.sha256(path.read_bytes()).hexdigest()
    paths.extend([dry_path, audit_path, case["b_audit"]])
    if case["summary"]:
        paths.append(case["summary"])
    return data, paths


def exit_endpoint(case, dry_path):
    """Exit observables the drawn history must reproduce, from the case's own record."""
    if case["summary"]:
        return json.loads(case["summary"].read_text())["dt5"]
    final = json.loads(dry_path.read_text())["final"]
    return {
        "particle_water_kg_kg_dry": final["particle_water_kg_kg_dry"],
        "hexane_mg_kg_dry": final["hexane_kg_kg_dry"] * 1e6,
        "temperature_c": final["temperature_k"] - 273.15,
    }


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(PAPER / "figures" / f"{stem}.{ext}", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=sorted(CASES), default="k0p32")
    args = parser.parse_args()
    case = CASES[args.case]
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "axes.linewidth": 0.7,
        }
    )
    faner_path = (
        ROOT
        / "paper/analysis/results_2026-10-07/faner_assumptions"
        / "volume_mean_N4/history/accepted_history.csv"
    )
    baseline_path = BASE / "dt_shared_complete_trial/accepted_history.csv"
    faner, baseline = rows(faner_path), rows(baseline_path)
    assert faner[0]["t_s"] == 0 and faner[-1]["t_s"] == 240
    finite, inputs = read_finite(case)
    regime_b_start_s = json.loads(inputs[0].read_text())["initial_time_s"]
    regime_c_start_s = next(r["time_s"] for r in finite if r.get("kind") == "extinction_projection")
    endpoint = exit_endpoint(case, case["dry"])
    assert (
        abs(finite[-1]["particle_water_kg_kg_dry"] - endpoint["particle_water_kg_kg_dry"]) < 1e-14
    )
    assert abs(finite[-1]["hexane_kg_kg_dry"] * 1e6 - endpoint["hexane_mg_kg_dry"]) < 1e-9
    assert abs(finite[-1]["temperature_k"] - 273.15 - endpoint["temperature_c"]) < 1e-10
    prefix = [r for r in baseline if r["time_s"] < finite[0]["time_s"]]
    finite = prefix + finite
    assert all(b["time_s"] > a["time_s"] for a, b in zip(finite, finite[1:]))
    assert finite[-1]["time_s"] == baseline[-1]["time_s"]
    bath, bath_audit = prescribed_bath(finite[-1]["time_s"])

    fig, axes = plt.subplots(2, 1, figsize=(6.8, 5.4), sharex=True, layout="constrained")
    t = [r["t_s"] for r in faner]
    axes[0].plot(
        t,
        [r["hexane_kg_kg_dry"] for r in faner],
        color=COLORS["hexane"],
        lw=1.7,
        label="Particle model",
    )
    axes[1].plot(
        t,
        [r["mean_temperature_k"] - 273.15 for r in faner],
        color=COLORS["temperature"],
        lw=1.7,
        label="Particle model",
    )
    for ax, quantity, label in zip(
        axes,
        ("solvent_loading", "particle_temperature"),
        ("Faner experiment", "Faner experiment"),
    ):
        data = measured(quantity)
        ax.errorbar(
            *zip(*[(r[0], r[1]) for r in data]),
            yerr=[r[2] for r in data],
            fmt="o",
            ms=3.6,
            capsize=2,
            color="#303840",
            lw=0.8,
            label=label,
        )
        ax.grid(False)
        ax.legend(loc="upper right" if ax is axes[0] else "lower right", frameon=False)
    axes[0].set(ylabel="Hexane (kg/kg dry meal)", ylim=(0, 0.51), xlim=(0, 240))
    axes[1].set(ylabel="Temperature (°C)", xlabel="Time from reported start (s)", ylim=(60, 124))
    axes[0].text(0.025, 0.92, "(a)", transform=axes[0].transAxes, weight="bold")
    axes[1].text(0.025, 0.81, "(b)", transform=axes[1].transAxes, weight="bold")
    axes[1].axhline(120, color="0.65", ls=":", lw=0.9, label="Prescribed vapour: 120°C")
    gas_readings = measured("gas_temperature")
    axes[1].plot(
        [r[0] for r in gas_readings],
        [r[1] for r in gas_readings],
        color="0.35",
        ls="--",
        lw=0.9,
        label="Measured vapour",
    )
    axes[1].legend(loc="lower right", frameon=False)
    save(fig, "fig_faner_particle_history_core2")

    # Water and core volume support the interpretation; they are not observed by Faner.
    fig, axes = plt.subplots(2, 1, figsize=(6.8, 4.2), sharex=True, layout="constrained")
    axes[0].plot(
        t, [r["water_kg_kg_dry"] for r in faner], color=COLORS["water"], label="Particle water"
    )
    axes[0].plot(
        t,
        [
            r["surface_water_kg_kg_dry"] if r["surface_water_kg_kg_dry"] != "" else math.nan
            for r in faner
        ],
        color=COLORS["water"],
        ls="--",
        label="Surface retained water",
    )
    axes[0].set_ylabel("Water (kg/kg dry meal)")
    axes[0].legend(frameon=False)
    axes[1].plot(t, [r["front_z"] for r in faner], color=COLORS["hexane"])
    axes[1].set(ylabel="Liquid-core volume fraction", xlabel="Time from reported start (s)")
    for ax in axes:
        ax.grid(alpha=0.15)
    save(fig, "fig_faner_supporting_states")

    fig, axes = plt.subplots(3, 1, figsize=(6.8, 7.6), sharex=True, layout="constrained")
    for data, linestyle, label in ((finite, "-", "Particle model"),):
        times = [r["time_s"] / 60 for r in data]
        water = [r["particle_water_kg_kg_dry"] for r in data]
        hexane = [r["hexane_kg_kg_dry"] for r in data]
        denominator = [1 + w + h for w, h in zip(water, hexane)]
        axes[0].plot(
            times,
            [1e6 * h / d for h, d in zip(hexane, denominator)],
            linestyle,
            color=COLORS["hexane"],
            lw=1.6,
            label=label,
        )
        axes[1].plot(
            times,
            [100 * w / d for w, d in zip(water, denominator)],
            linestyle,
            color=COLORS["water"],
            lw=1.6,
        )
        axes[2].plot(
            times,
            [
                (
                    r["temperature_k"] - 273.15
                    if r.get("temperature_k") not in (None, "")
                    else math.nan
                )
                for r in data
            ],
            linestyle,
            color=COLORS["temperature"],
            lw=1.6,
        )
    bath_times = [r[0] / 60 for r in bath]
    bath_quantities = (
        [r[2] for r in bath],
        [1 - r[2] for r in bath],
        [r[1] - 273.15 for r in bath],
    )
    for ax, quantity, ylabel in zip(
        axes,
        bath_quantities,
        ("Bath hexane mole fraction", "Bath water mole fraction", "Bath temperature (°C)"),
    ):
        ax.axvspan(0, PD_END_S / 60, color="0.9", alpha=0.35, zorder=0)
        for event in (regime_b_start_s, regime_c_start_s):
            ax.axvline(event / 60, color="0.65", lw=0.7, ls=":", zorder=1)
        secondary = ax.twinx()
        secondary.plot(bath_times, quantity, color="0.35", ls="--", lw=1.1)
        secondary.set_ylabel(ylabel, color="0.35")
        secondary.tick_params(axis="y", colors="0.35", labelsize=8)
        secondary.spines["right"].set_visible(True)
        secondary.spines["right"].set_color("0.55")
        secondary.set_ylim((45, 117) if ax is axes[2] else (0, 1))
        if ax is not axes[2]:
            secondary.set_yticks([0, 0.25, 0.5, 0.75, 1])
        # Keep the primary histories and white inset backgrounds above the bath.
        secondary.set_zorder(1)
        ax.set_zorder(2)
        ax.patch.set_visible(False)
    # Published outlet intervals are contextual references placed at 30 min.
    # The simulated clock and its actual 29.44-min endpoint are unchanged.
    benchmark_x = 30.0
    for ax, bounds, text in zip(
        axes, ((100, 500), (16, 22), (105, 110)), ("100–500 ppm", "16–22%", "105–110°C")
    ):
        mid = sum(bounds) / 2
        ax.errorbar(
            benchmark_x,
            mid,
            yerr=[[mid - bounds[0]], [bounds[1] - mid]],
            fmt="none",
            ecolor=OUTLET_COLOR,
            elinewidth=2.6,
            capsize=6,
            capthick=1.2,
            zorder=4,
        )
        ax.text(
            benchmark_x - 0.55,
            (bounds[1] * 1.08 if ax is axes[0] else bounds[1] + (0.15 if ax is axes[1] else 0.65)),
            text,
            ha="right",
            va="bottom",
            fontsize=8,
            color=OUTLET_COLOR,
        )
        ax.grid(False)
        ax.tick_params(axis="x", labelbottom=True)
    axes[0].set(yscale="log", ylabel="Hexane (ppm wet meal)", ylim=(80, 4e5))
    axes[0].plot([], [], color="0.35", ls="--", lw=1.1, label="Prescribed bath")
    axes[0].legend(frameon=False, loc="upper right", bbox_to_anchor=(0.72, 0.99), fontsize=7)
    axes[0].text(30, 1e5, "Published\noutlet ranges", ha="right", fontsize=8, color=OUTLET_COLOR)
    axes[1].set(ylabel="Particle moisture (% wet basis)", ylim=(5, 24))
    axes[2].set(ylabel="Mean temperature (°C)", xlabel="Time from feed (min)", ylim=(45, 117))
    axes[2].set_xlim(0, 31.5)
    axes[2].set_xticks(range(0, 31, 5))
    regime_axis = axes[0].secondary_xaxis("top")
    regime_axis.set_xticks(
        [
            regime_b_start_s / 120,
            (regime_b_start_s + regime_c_start_s) / 120,
            (regime_c_start_s + finite[-1]["time_s"]) / 120,
        ],
        labels=["Regime A", "Regime B", "Regime C"],
    )
    regime_axis.tick_params(length=3, labelsize=8, pad=3)
    regime_axis.spines["top"].set_visible(True)
    equipment_axis = axes[2].secondary_xaxis(-0.31)
    equipment_axis.set_xticks(
        [PD_END_S / 120, (PD_END_S + finite[-1]["time_s"]) / 120],
        labels=["Predesolv.", "Main trays"],
    )
    equipment_axis.tick_params(length=0, labelsize=9, pad=3)
    equipment_axis.spines["bottom"].set_visible(False)
    for start, stop in ((0, PD_END_S / 60), (PD_END_S / 60, finite[-1]["time_s"] / 60)):
        axes[2].plot(
            [start, stop],
            [-0.31, -0.31],
            color="0.55",
            lw=0.7,
            transform=axes[2].get_xaxis_transform(),
            clip_on=False,
        )
        for edge in (start, stop):
            axes[2].plot(
                [edge, edge],
                [-0.31, -0.28],
                color="0.55",
                lw=0.7,
                transform=axes[2].get_xaxis_transform(),
                clip_on=False,
            )
    for ax, letter in zip(axes, "abc"):
        ax.text(0.025, 0.89, f"({letter})", transform=ax.transAxes, weight="bold")
    # Empty lower-middle region keeps the full moisture history visible.
    film_ax = axes[1].inset_axes([0.48, 0.19, 0.25, 0.27])
    for data, ls in ((finite, "-"),):
        film_ax.plot(
            [r["time_s"] / 60 for r in data],
            [r["external_water_kg_kg_dry"] for r in data],
            ls,
            color=COLORS["water"],
            lw=1,
        )
    film_ax.set(xlim=(7.8, 14), ylim=(0, 0.13))
    film_ax.set_title("External water film", fontsize=7, pad=2)
    film_ax.set_xlabel("Time from feed (min)", fontsize=7, labelpad=1)
    film_ax.set_ylabel("kg/kg dry", fontsize=7, labelpad=2)
    film_ax.set_xticks([8, 10, 12, 14])
    film_ax.tick_params(labelsize=7)
    # Keep its centre in the same region while reducing the inset footprint.
    core_ax = axes[0].inset_axes([0.42, 0.395, 0.28, 0.22])
    short = [r for r in finite if 490 <= r["time_s"] <= 505]
    core_ax.plot(
        [r["time_s"] - PD_END_S for r in short],
        [r.get("front_z", r.get("wet_volume_fraction", 0)) for r in short],
        color=COLORS["hexane"],
        lw=1,
    )
    core_ax.set(xlim=(0, 15), ylim=(0, 1.06))
    core_ax.set_title("Liquid-hexane core", fontsize=7, pad=2)
    core_ax.set_ylabel("V/V", fontsize=7, labelpad=2)
    core_ax.set_xticks([0, 5, 10, 15])
    core_ax.set_yticks([0, 1])
    transition_b = regime_b_start_s - PD_END_S
    transition_c = regime_c_start_s - PD_END_S
    core_ax.axvspan(transition_b, transition_c, color="#e5edf3", zorder=0)
    for event in (transition_b, transition_c):
        core_ax.axvline(event, color="0.55", lw=0.7, ls=":")
    for start, stop, letter in (
        (0, transition_b, "A"),
        (transition_b, transition_c, "B"),
        (transition_c, 15, "C"),
    ):
        core_ax.text((start + stop) / 2, 0.90, letter, ha="center", va="top", fontsize=7)
    core_ax.set_xlabel("After PD exit (s)", fontsize=7, labelpad=1)
    core_ax.tick_params(labelsize=7)
    for inset in (core_ax, film_ax):
        for spine in inset.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(0.6)
    save(fig, "fig_dt_particle_history_core2")
    inputs.extend(
        [
            faner_path,
            faner_path.parent / "summary.json",
            ROOT / "paper/analysis/results_2026-10-07/faner_assumptions/source_mean_n4_plan.json",
            baseline_path,
            BATH_PATH,
            FEED_PATH,
            ROOT / "paper/analysis/harness/experiment_core2_dt_feed.py",
            ROOT / "paper/analysis/harness/experiment_core2_dt_face.py",
            ROOT / "paper/analysis/harness/dtdc_particle_destiny_march.py",
            Path(__file__),
            ROOT / "paper/analysis/datasets/faner2019_v2/curves.csv",
            ROOT / "paper/analysis/datasets/faner2019_v2/adjudication_per_point.csv",
        ]
    )
    manifest = {
        "source_sha256": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs
        },
        "dt_case": args.case,
        "dt_extra_contact_parameters": {
            "shell_mobility_factor": case["factor"],
            "K_mol_m2_s": case["K"],
        },
        "outlet_context": {
            "hexane_ppm_wet_basis": [100, 500],
            "moisture_percent_wet_basis": [16, 22],
            "temperature_c": [105, 110],
        },
        "outlet_ranges_are_transient_targets": False,
        "outlet_reference_position_min": 30.0,
        "dt_curves_shown": "proposed particle model including liquid-water imbibition",
        "dt_prescribed_bath": bath_audit,
        "dt_events_s": {
            "pd_exit": PD_END_S,
            "regime_a_to_b": regime_b_start_s,
            "regime_b_to_c": regime_c_start_s,
            "computed_exit": finite[-1]["time_s"],
        },
        "outlet_context_sources": {
            "moisture": "Witte, printed p.97: 16-22% at DT exit",
            "residual_and_temperature": "Kemper, printed pp.112-113: 100-500 ppm and 105-110 C; ppm mass basis not explicit",
        },
        "faner_clock_shift_s": 0,
        "faner_size_ensemble_claimed": False,
        "faner_radius_m": 0.000885,
        "faner_radial_cells": 4,
        "faner_radius_basis": "reported mean volume-equivalent diameter divided by two",
        "faner_primary_bath_temperature_k": 393.15,
        "faner_measured_vapour_trace_is_prescribed_in_primary": False,
    }
    (PAPER / "figures/fig_particle_experiments.sources.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
