"""Recalculate supplementary table values from the supplied numerical curves."""

import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/supplementary"


def read(name):
    return json.loads((DATA / name).read_text())


def crossing(t, x, level):
    for i in range(len(t) - 1):
        if x[i] >= level >= x[i + 1] and x[i] != x[i + 1]:
            return float(t[i] + (t[i + 1] - t[i]) * (x[i] - level) / (x[i] - x[i + 1]))
    raise ValueError("curve does not cross requested loading")


def main():
    observations = list(csv.DictReader((DATA / "faner2008_observations.csv").open()))
    scores = []
    for case in read("scalar_drying_curves.json")["cases"]:
        series = "X_girasol" if "sunflower" in case["condition"] else "X_soja"
        points = sorted(
            (r for r in observations if r["figure"] == "4.21" and r["series"] == series),
            key=lambda r: float(r["x"]),
        )
        measured_t = [float(p["x"]) for p in points]
        measured_x = [float(p["y"]) for p in points]
        start = crossing(measured_t, measured_x, 0.20)
        window = [p for p in points if float(p["x"]) >= start]
        span = max(float(p["y"]) for p in window) - min(float(p["y"]) for p in window)
        result = dict(case=case["case"], cells=case["cells"], time_level=case["time_level"])
        for mode in ("dec", "full"):
            indices = [i for i, n in enumerate(case["step_index"]) if mode == "full" or n % 5 == 0]
            if mode == "dec":
                indices.append(len(case["step_index"]) - 1)
            times = [0.0, case["constant_rate_duration_s"]] + [
                case["constant_rate_duration_s"] + case["tail_time_s"][i] for i in indices
            ]
            loads = [measured_x[0], case["initial_tail_loading"]] + [
                case["loading_kg_kg"][i] for i in indices
            ]
            predictions = [float(np.interp(float(p["x"]), times, loads)) for p in window]
            rms = (
                math.sqrt(
                    sum((v - float(p["y"])) ** 2 for v, p in zip(predictions, window)) / len(window)
                )
                / span
            )
            inside = sum(
                float(p["y_reading_lower"]) <= v <= float(p["y_reading_upper"])
                for v, p in zip(predictions, window)
            )
            expected = case["expected_scores"][mode]
            if abs(rms - expected["rms_over_span"]) > 1e-12 or inside != expected["in_band"]:
                raise ValueError(
                    "scalar score differs from its original calculation: " + case["case"]
                )
            result[mode] = dict(rms_over_span=rms, in_band=inside)
        scores.append(result)
    ensemble = read("size_ensemble.json")
    ensemble_results = {}
    for label, curve in ensemble["curves"].items():
        mean = curve["mean_curve"]
        interval = crossing(mean["t"], mean["X"], 0.10) - crossing(mean["t"], mean["X"], 0.20)
        source = ensemble["source"]["solvent_loading"]
        rms = math.sqrt(
            sum((float(np.interp(t, mean["t"], mean["X"])) - x) ** 2 for t, x in source)
            / len(source)
        )
        if (
            abs(interval - curve["falling_interval_s"]) > 1e-10
            or abs(rms - curve["rms_loading_25"]) > 1e-12
        ):
            raise ValueError("ensemble table differs from original calculation")
        ensemble_results[label] = dict(falling_interval_s=interval, rms_loading_kg_kg=rms)
    ratios = []
    for row in csv.DictReader((DATA / "uptake_coefficients.csv").open()):
        if row.get("duplicate_of"):
            continue
        ratio = float(row["converted_storage_consistent_m2_s"]) / float(
            row["tabulated_apparent_m2_s"]
        )
        if abs(ratio - float(row["measured_retardation"])) > 1e-10:
            raise ValueError("uptake storage conversion differs")
        ratios.append(ratio)
    result = dict(
        passed=True,
        scalar_cases=len(scores),
        scalar_scores=scores,
        ensemble=ensemble_results,
        uptake_conversion_ratio_range=[min(ratios), max(ratios)],
        scope="Recalculation from archived numerical curves; no fresh particle-model histories run.",
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
