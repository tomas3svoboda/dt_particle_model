"""Serial, unfitted Faner geometry/moisture experiments using existing marchers.

This controller chooses numerical intervals and invokes the existing native
feed, radial, face, extinction and dry-state solvers. It defines no physical
law, supplies no drying continuation and never selects a case by its score.
"""

from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import subprocess
import sys

from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr

ROOT = Path(__file__).resolve().parents[3]
HARNESS = ROOT / "paper/analysis/harness"
BASE = ROOT / "paper/analysis/results_2026-10-07/faner_assumptions"
ARCHIVE = ROOT / "docs/evidence/paper1_faner2019_march_step4_2026-09-30/evidence/step4"


def plan():
    """Freeze cases from source summaries, before reading any new scores."""
    source = json.loads((ARCHIVE / "derive_ensemble.json").read_text())
    key = "lognormal|p=0.951|trunc=0.5-4|basis=number|R=declared|closed_form|tail=hold"
    dist = source["traces"]["soybean"]["cases"][key]["distribution"]
    mu, sigma = dist["mu_ln_mm"], dist["sigma_ln"]
    number_norm = ndtr((math.log(4) - mu) / sigma) - ndtr((math.log(0.5) - mu) / sigma)

    def number_density(d):
        x = (math.log(d) - mu) / sigma
        return math.exp(-0.5 * x * x) / (d * sigma * math.sqrt(2 * math.pi) * number_norm)

    total = quad(lambda d: d**3 * number_density(d), 0.5, 4, epsabs=1e-11)[0]

    def cdf(d):
        return quad(lambda v: v**3 * number_density(v), 0.5, d, epsabs=1e-11)[0] / total

    cuts = [0.5] + [brentq(lambda d: cdf(d) - q, 0.5, 4) for q in (1 / 3, 2 / 3)] + [4.0]
    classes = []
    for lo, hi in zip(cuts, cuts[1:]):
        mass = quad(lambda d: d**3 * number_density(d), lo, hi, epsabs=1e-11)[0]
        diameter = quad(lambda d: d**4 * number_density(d), lo, hi, epsabs=1e-11)[0] / mass
        classes.append(
            {
                "dry_mass_weight": mass / total,
                "mean_diameter_mm": diameter,
                "diameter_bounds_mm": [lo, hi],
            }
        )
    cases = {
        "baseline": {
            "radius_m": 0.0008078927451503978,
            "water": 0.025,
            "basis": "inherited representative sphere",
        },
        "volume_mean": {
            "radius_m": 0.000885,
            "water": 0.025,
            "basis": "reported number-mean volume-equivalent diameter / 2",
        },
        "moisture_05": {
            "radius_m": 0.000885,
            "water": 0.05,
            "basis": "mean volume-equivalent sphere; declared moisture sensitivity",
        },
        "moisture_10": {
            "radius_m": 0.000885,
            "water": 0.10,
            "basis": "mean volume-equivalent sphere; declared moisture sensitivity",
        },
    }
    for i, cls in enumerate(classes, 1):
        for label, factor in (("surface", 0.74), ("volume", 1.0)):
            cases[f"psd_{label}_{i}"] = {
                "radius_m": factor * cls["mean_diameter_mm"] * 0.0005,
                "water": 0.025,
                "basis": f"{label} equivalent sphere; independently constrained lognormal",
                **cls,
            }
    return {
        "selection": "source geometry, two declared sphere mappings and predeclared moisture range; no new curve enters case selection",
        "size_family_is_measured_histogram": False,
        "size_quadrature": "three equal dry-mass bins; within-bin mass-mean diameter; coarse sensitivity, not quadrature convergence",
        "distribution": dist,
        "fixed": {
            "cells": 2,
            "initial_hexane": 0.475,
            "oil": 0.0195,
            "bath_temperature_k": 393.15,
            "film_multiplier": 0.14578543636088787,
            "coefficient_case": "db_max_dapp_max",
            "continued_positive_water": True,
        },
        "boundary_is_independent_of_Faner": False,
        "cases": cases,
    }


def invoke(folder, kind, args):
    outputs = sorted(
        p for p in folder.glob(f"{kind}_*.json") if not p.name.endswith(".command.json")
    )
    index = max((int(p.stem.rsplit("_", 1)[-1]) for p in outputs), default=0) + 1
    out = folder / f"{kind}_{index:02d}.json"
    if any(p.exists() for p in (out, out.with_suffix(".command.json"), out.with_suffix(".log"))):
        raise ValueError("campaign refuses to overwrite any execution evidence")
    script = {
        "birth": "experiment_core2_faner.py",
        "radial": "experiment_core2_faner.py",
        "event": "experiment_particle_event.py",
        "dry": "experiment_particle_dry.py",
    }[kind]
    command = [sys.executable, str(HARNESS / script), "--out", str(out), *map(str, args)]
    print(json.dumps({"case": folder.name, "stage": kind, "output": out.name}), flush=True)
    with out.with_suffix(".command.json").open("x") as stream:
        json.dump({"argv": command}, stream, indent=2)
    with out.with_suffix(".log").open("x", encoding="utf-8") as stream:
        result = subprocess.run(
            command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, check=False
        )
    if not out.exists():
        raise RuntimeError(
            f"native CLI exited {result.returncode}; inspect {out.with_suffix('.log')}"
        )
    record = json.loads(out.read_text(encoding="utf-8"))
    if not record.get("sources_unchanged_during_run"):
        raise RuntimeError("source identity changed during native execution")
    return out, record


def checkpoint(out):
    return out.with_suffix(".accepted.pkl")


def advance_case(
    name,
    cfg,
    tip=None,
    *,
    half_decades=False,
    regular_event_tolerance=1e-10,
    departure_dt=0.01,
    departure_tolerance=1e-10,
    departure_native_budget=None,
    regular_root_tolerance=1e-9,
    extended_regional_root=False,
):
    folder = BASE / name
    cells = cfg.get("cells", 2)
    departure_budget_arguments = (
        ["--native-budget", departure_native_budget] if departure_native_budget else []
    )
    folder.mkdir(parents=True, exist_ok=True)
    if tip is None:
        out, record = invoke(
            folder,
            "birth",
            [
                "--radius",
                cfg["radius_m"],
                "--water",
                cfg["water"],
                "--cells",
                cells,
                "--birth-dt",
                0.01,
                "--coefficient-case",
                "db_max_dapp_max",
                "--continued-water",
                "--leading-seed",
                "--steps",
                0,
            ],
        )
        if not checkpoint(out).exists():
            return {
                "complete": False,
                "tip": str(out),
                "error": record.get("error"),
                "stage": "birth",
            }
    else:
        out = tip.resolve()
        record = json.loads(out.read_text())
    if record.get("event") == "arrival":
        new, result = invoke(
            folder,
            "event",
            [
                "--resume",
                checkpoint(out),
                "--event",
                "departure",
                "--dt",
                departure_dt,
                "--regional-tolerance",
                departure_tolerance,
                "--tangent-tolerance",
                1e-9,
                *departure_budget_arguments,
            ],
        )
        if not result.get("accepted"):
            return {
                "complete": False,
                "tip": str(out),
                "failure": str(new),
                "error": result.get("error"),
                "stage": "departure",
            }
        out, record = new, result
    for tranche in range(60):
        state = record["last_accepted"]
        z, t = state["front_z"], state["time_s"]
        if z == 0 or t >= 240:
            break
        nw = state.get("wet_pieces", 1)
        next_z = ((nw - 1) / cells) ** 3
        if nw > 1 and z - next_z <= 0.065:
            event = "arrival"
            rate = front_rate(record)
            guess = min(8.0, max(0.01, (z - next_z) / rate))
            args = ["--event", event, "--dt", guess, "--upper", 10, "--regional-tolerance", 1e-10]
            new, result = invoke(folder, "event", ["--resume", checkpoint(out), *args])
            if not result.get("accepted"):
                return {
                    "complete": False,
                    "tip": str(out),
                    "failure": str(new),
                    "error": result.get("error"),
                    "stage": event,
                }
            out, record = new, result
            new, result = invoke(
                folder,
                "event",
                [
                    "--resume",
                    checkpoint(out),
                    "--event",
                    "departure",
                    "--dt",
                    departure_dt,
                    "--regional-tolerance",
                    departure_tolerance,
                    "--tangent-tolerance",
                    1e-9,
                    *departure_budget_arguments,
                ],
            )
            if not result.get("accepted"):
                return {
                    "complete": False,
                    "tip": str(out),
                    "failure": str(new),
                    "error": result.get("error"),
                    "stage": "departure",
                }
            out, record = new, result
        elif nw == 1 and z <= 0.02:
            targets = (
                1e-3,
                1e-4,
                1e-5,
                1e-6,
                1e-7,
                1e-8,
                1e-9,
                1e-10,
                1e-11,
                1e-12,
                1e-13,
                1e-14,
            )
            if half_decades:
                targets = tuple(10 ** (-k / 2) for k in range(6, 29))
            for target in targets:
                if target >= record["last_accepted"]["front_z"]:
                    continue
                guess = min(
                    5.0,
                    max(1e-5, (record["last_accepted"]["front_z"] - target) / front_rate(record)),
                )
                tolerance = (
                    2e-8 if record["last_accepted"]["front_z"] <= 1e-7 else regular_event_tolerance
                )
                new, result = invoke(
                    folder,
                    "event",
                    [
                        "--resume",
                        checkpoint(out),
                        "--event",
                        "target",
                        "--target",
                        target,
                        "--dt",
                        guess,
                        "--budget",
                        4000,
                        "--regional-tolerance",
                        tolerance,
                        "--condition-ceiling",
                        1e18,
                    ],
                )
                if not result.get("accepted"):
                    # A smaller positive-volume interval improves conditioning;
                    # it does not book the rejected candidate or change a law.
                    midpoint = math.sqrt(record["last_accepted"]["front_z"] * target)
                    mid, mid_result = invoke(
                        folder,
                        "event",
                        [
                            "--resume",
                            checkpoint(out),
                            "--event",
                            "target",
                            "--target",
                            midpoint,
                            "--dt",
                            max(1e-6, guess * 0.5),
                            "--budget",
                            12000,
                            "--regional-tolerance",
                            1e-10,
                            "--condition-ceiling",
                            1e18,
                        ],
                    )
                    if mid_result.get("accepted"):
                        out, record = mid, mid_result
                        new, result = invoke(
                            folder,
                            "event",
                            [
                                "--resume",
                                checkpoint(out),
                                "--event",
                                "target",
                                "--target",
                                target,
                                "--dt",
                                max(1e-6, mid_result["ledger"]["dt_s"] * 0.5),
                                "--budget",
                                12000,
                                "--regional-tolerance",
                                1e-10,
                                "--condition-ceiling",
                                1e18,
                            ],
                        )
                if not result.get("accepted"):
                    return {
                        "complete": False,
                        "tip": str(out),
                        "failure": str(new),
                        "error": result.get("error"),
                        "stage": "positive target",
                    }
                out, record = new, result
                if target <= 1e-10:
                    new, result = invoke(
                        folder,
                        "event",
                        [
                            "--resume",
                            checkpoint(out),
                            "--event",
                            "projection",
                            "--regional-tolerance",
                            2e-8,
                            "--condition-ceiling",
                            1e18,
                            "--event-time-tolerance",
                            0.02,
                        ],
                    )
                    if result.get("accepted"):
                        out, record = new, result
                        break
            if record["last_accepted"]["front_z"] != 0:
                return {
                    "complete": False,
                    "tip": str(out),
                    "stage": "extinction",
                    "error": "native projection not accepted",
                }
        else:
            rate = front_rate(record)
            duration = min(10.0, max(0.05, 0.65 * (z - next_z) / rate))
            dt = 0.15 if z > 0.94 else (0.2 if z < 0.03 else 0.5)
            previous_dt = record.get("ledger", {}).get("dt_s")
            accepted_steps = [a for a in record.get("attempts", []) if a.get("accepted")]
            if accepted_steps:
                previous_dt = accepted_steps[-1].get("dt_s", previous_dt)
            if previous_dt is not None and previous_dt > 0:
                dt = min(dt, previous_dt * 1.5)
            horizon = min(240.0, t + duration)
            new, result = invoke(
                folder,
                "radial",
                [
                    "--resume",
                    checkpoint(out),
                    "--continued-water",
                    "--regional-tolerance",
                    regular_root_tolerance,
                    *(["--extended-regional-root"] if extended_regional_root else []),
                    "--difference-mode",
                    "central",
                    "--steps",
                    100,
                    "--continuation-dt",
                    dt,
                    "--horizon",
                    horizon,
                    "--wall-budget",
                    180,
                ],
            )
            if not checkpoint(new).exists() or result["last_accepted"]["time_s"] <= t:
                return {
                    "complete": False,
                    "tip": str(out),
                    "failure": str(new),
                    "error": result.get("error", result.get("stop")),
                    "stage": "radial",
                }
            out, record = new, result
    if record["last_accepted"]["front_z"] == 0:
        new, result = invoke(
            folder,
            "dry",
            ["--resume", checkpoint(out), "--dt", 2.5, "--steps", 150, "--wall-budget", 180],
        )
        out, record = new, result
    complete = record["last_accepted"]["time_s"] >= 240
    result = {
        "complete": complete,
        "tip": str(out),
        "last_accepted": record["last_accepted"],
        "error": record.get("error"),
    }
    if complete:
        command = [
            sys.executable,
            str(HARNESS / "export_core2_faner_history.py"),
            "--tip",
            str(out),
            "--out",
            str(folder / "history"),
            "--ledger-tolerance",
            "4e-10",
        ]
        with (folder / "export.log").open("w", encoding="utf-8") as stream:
            rc = subprocess.run(
                command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, check=False
            )
        result["export_passed"] = rc.returncode == 0
    return result


def front_rate(record):
    accepted = [a for a in record.get("attempts", []) if a.get("accepted") and a.get("after")]
    if len(accepted) >= 2:
        a, b = accepted[-2:]
        return max(1e-8, (a["after"]["front_z"] - b["after"]["front_z"]) / b["dt_s"])
    state = record["last_accepted"]
    z = state["front_z"]
    if "ledger" in record and record.get("accepted_parent"):
        parent = Path(record["accepted_parent"]["result"])
        if not parent.is_absolute():
            parent = ROOT / parent
        old = json.loads(parent.read_text())["last_accepted"]
        duration = record["ledger"]["dt_s"]
        if duration > 0 and old["front_z"] > z:
            return max(1e-8, (old["front_z"] - z) / duration)
    return max(1e-8, z ** (2 / 3) * 0.05)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", nargs="+")
    parser.add_argument("--resume-tip", type=Path)
    parser.add_argument("--target-half-decades", action="store_true")
    parser.add_argument("--regular-event-tolerance", type=float, default=1e-10)
    parser.add_argument("--extra-case-plan", type=Path)
    parser.add_argument("--departure-dt", type=float, default=0.01)
    parser.add_argument("--departure-tolerance", type=float, default=1e-10)
    parser.add_argument("--departure-native-budget", type=int)
    parser.add_argument("--regular-root-tolerance", type=float, default=1e-9)
    parser.add_argument("--extended-regional-root", action="store_true")
    args = parser.parse_args()
    if not 0 < args.regular_event_tolerance <= 1e-9:
        raise ValueError("local root tolerance must remain <=1e-9; native ledgers are unchanged")
    controller_entry_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for key in (
        "PYTHONHASHSEED",
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        if os.environ.get(key) != "1":
            raise ValueError(f"single-worker environment requires {key}=1")
    BASE.mkdir(parents=True, exist_ok=True)
    manifest = plan()
    path = BASE / "plan.json"
    if path.exists():
        if json.loads(path.read_text()) != manifest:
            raise ValueError("frozen campaign plan changed")
    else:
        path.write_text(json.dumps(manifest, indent=2))
    extra_plan_hash = None
    if args.extra_case_plan:
        extra = json.loads(args.extra_case_plan.read_text(encoding="utf-8"))
        if set(extra) & set(manifest["cases"]):
            raise ValueError("additional declared cases cannot replace frozen cases")
        manifest["cases"].update(extra)
        extra_plan_hash = hashlib.sha256(args.extra_case_plan.read_bytes()).hexdigest()
    outcomes = BASE / "outcomes.json"
    results = json.loads(outcomes.read_text()) if outcomes.exists() else {}
    for name in args.case or []:
        result = advance_case(
            name,
            manifest["cases"][name],
            args.resume_tip,
            half_decades=args.target_half_decades,
            regular_event_tolerance=args.regular_event_tolerance,
            departure_dt=args.departure_dt,
            departure_tolerance=args.departure_tolerance,
            departure_native_budget=args.departure_native_budget,
            regular_root_tolerance=args.regular_root_tolerance,
            extended_regional_root=args.extended_regional_root,
        )
        result["controller_sha256_at_entry"] = controller_entry_hash
        result["controller_sha256_at_exit"] = hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest()
        result["controller_arguments"] = {
            "target_half_decades": args.target_half_decades,
            "regular_event_tolerance": args.regular_event_tolerance,
            "resume_tip": str(args.resume_tip),
            "extra_case_plan": str(args.extra_case_plan),
            "departure_dt": args.departure_dt,
            "departure_tolerance": args.departure_tolerance,
            "departure_native_budget": args.departure_native_budget,
            "regular_root_tolerance": args.regular_root_tolerance,
            "extended_regional_root": args.extended_regional_root,
        }
        result["additional_plan_sha256"] = extra_plan_hash
        results[name] = result
        outcomes.write_text(json.dumps(results, indent=2))
        print(json.dumps({"case": name, **result}), flush=True)


if __name__ == "__main__":
    main()
