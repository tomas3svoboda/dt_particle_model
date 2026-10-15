"""Serial full-history property cases using existing native paper controllers."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
HARNESS = ROOT / "paper/analysis/harness"
RUNS = ROOT / "paper/analysis/temperature_property_runs"


def command(driver, args, reference, constant_D=False):
    return [
        sys.executable,
        str(HARNESS / "run_temperature_property_driver.py"),
        "--Db-ref",
        str(reference),
        *(["--constant-D"] if constant_D else []),
        "--driver",
        str(driver),
        "--",
        *map(str, args),
    ]


def invoke(folder, label, driver, arguments, reference, constant_D=False):
    index = len(list(folder.glob(f"{label}_*.command.json"))) + 1
    out = folder / f"{label}_{index:02d}.json"
    argv = command(driver, ["--out", out, *arguments], reference, constant_D)
    out.with_suffix(".command.json").write_text(json.dumps({"argv": argv}, indent=2))
    print(json.dumps({"stage": label, "out": str(out)}), flush=True)
    with out.with_suffix(".log").open("x", encoding="utf-8") as log:
        rc = subprocess.run(argv, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    if not out.exists():
        raise RuntimeError(
            f"driver returned {rc.returncode} without result: {out.with_suffix('.log')}"
        )
    return out, json.loads(out.read_text())


def faner(
    folder,
    reference,
    constant_D,
    resume,
    radius,
    water,
    cells,
    departure_dt,
    departure_tolerance,
    half_decades,
    terminal_entry_tolerance,
    terminal_duration_from_history,
    target_root_tolerance,
):
    import faner_assumption_campaign as existing

    existing.BASE = folder.parent
    original_run = subprocess.run

    def wrapped_run(argv, *args, **kwargs):
        if (
            isinstance(argv, list)
            and len(argv) > 1
            and argv[0] == sys.executable
            and Path(argv[1]).resolve().is_relative_to(HARNESS)
        ):
            argv = list(argv)
            if (
                Path(argv[1]).name == "experiment_particle_event.py"
                and "--event" in argv
                and argv[argv.index("--event") + 1] == "target"
            ):
                parent = Path(argv[argv.index("--resume") + 1])
                parent = parent.with_name(parent.name.replace(".accepted.pkl", ".json"))
                parent_record = json.loads(parent.read_text())
                fraction = parent_record["last_accepted"]["front_z"]
                if target_root_tolerance is not None and fraction > 1e-7:
                    argv[argv.index("--regional-tolerance") + 1] = str(target_root_tolerance)
                if terminal_entry_tolerance is not None and 1e-8 <= fraction <= 1e-7:
                    index = argv.index("--regional-tolerance") + 1
                    argv[index] = str(min(float(argv[index]), terminal_entry_tolerance))
                if terminal_duration_from_history and fraction <= 0.02:
                    prior = parent_record.get("accepted_parent")
                    if prior:
                        prior_path = Path(prior["result"])
                        if not prior_path.is_absolute():
                            prior_path = ROOT / prior_path
                        assert (
                            hashlib.sha256(prior_path.read_bytes()).hexdigest()
                            == prior["result_sha256"]
                        )
                        previous = json.loads(prior_path.read_text())
                        old_fraction = previous["last_accepted"]["front_z"]
                        duration = parent_record.get("ledger", {}).get("dt_s", 0)
                        old_duration = previous.get("ledger", {}).get("dt_s", 0)
                        target = float(argv[argv.index("--target") + 1])
                        if (
                            old_duration > 0
                            and duration > 0
                            and old_fraction > fraction
                            and math.isclose(
                                fraction / old_fraction, target / fraction, rel_tol=1e-10
                            )
                        ):
                            ratio = duration / old_duration
                            if 0 < ratio < 1:
                                argv[argv.index("--dt") + 1] = str(duration * ratio)
            argv = command(Path(argv[1]), argv[2:], reference, constant_D)
            if "--out" in argv:
                path = Path(argv[argv.index("--out") + 1])
                sidecar = path.parent / f"executed_{path.stem}.command.json"
                with sidecar.open("x", encoding="utf-8") as stream:
                    json.dump(
                        {"argv": argv, "duration_is_numerical_guess_only": True}, stream, indent=2
                    )
        return original_run(argv, *args, **kwargs)

    subprocess.run = wrapped_run
    try:
        return existing.advance_case(
            folder.name,
            {"radius_m": radius, "water": water, "cells": cells},
            tip=resume,
            departure_dt=departure_dt,
            departure_tolerance=departure_tolerance,
            half_decades=half_decades,
        )
    finally:
        subprocess.run = original_run


def dt(folder, reference, constant_D, resume, b_driver, birth_seed):
    tip = resume
    seed = birth_seed
    parts = []
    if resume:
        previous = json.loads(resume.read_text())
        # Retain the ordered B chain by following checkpoint parents.
        path = resume
        while True:
            parts.append(path)
            parent = previous.get("parent_checkpoint")
            if not parent:
                break
            path = Path(parent).with_suffix("")
            path = path.with_name(path.name.replace(".accepted", "") + ".json")
            previous = json.loads(path.read_text())
        parts.reverse()
    for _ in range(30):
        arguments = [
            "--factor",
            "1000",
            "--contact-conductance",
            ".32",
            "--dt",
            ".1",
            "--steps",
            "300",
            "--until",
            "510",
        ]
        if tip:
            arguments += ["--resume", tip.with_suffix(".accepted.pkl")]
        if seed:
            arguments += ["--birth-seed", seed]
        out, result = invoke(folder, "b", ROOT / b_driver, arguments, reference, constant_D)
        if out.with_suffix(".accepted.pkl").exists():
            tip = out
            parts.append(out)
        if result["status"] == "REGIME_B_COMPLETED":
            break
        if result.get("error"):
            debug = out.with_suffix(".birth_debug.json")
            if tip is None and seed is None and debug.exists():
                seed = debug
                continue
            return {
                "complete": False,
                "tip": str(tip),
                "failure": str(out),
                "error": result["error"],
                "parts": list(map(str, parts)),
            }
    else:
        return {
            "complete": False,
            "tip": str(tip),
            "error": "B tranche count exhausted",
            "parts": list(map(str, parts)),
        }
    out, result = invoke(
        folder,
        "c",
        HARNESS / "experiment_core2_radial_fast_wetting.py",
        [
            "--resume",
            tip.with_suffix(".accepted.pkl"),
            "--factor",
            "1000",
            "--contact-conductance",
            ".32",
            "--dt",
            "5",
            "--until",
            "1766.5355008060146",
            "--steps",
            "1500",
            "--max-halvings",
            "12",
        ],
        reference,
        constant_D,
    )
    if not result.get("bath_exit_reached"):
        return {
            "complete": False,
            "tip": str(out),
            "error": result.get("error"),
            "parts": list(map(str, parts)),
        }
    audit, report = invoke(
        folder,
        "audit",
        HARNESS / "audit_finite_wetting_b.py",
        ["--result", tip, "--parts", *parts, "--dry", out],
        reference,
        constant_D,
    )
    return {
        "complete": True,
        "audit_passed": report["passed"],
        "tip": str(out),
        "b_tip": str(tip),
        "parts": list(map(str, parts)),
        "full_audit": str(audit),
        "final": result["final"],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study", choices=["dt", "faner"], required=True)
    p.add_argument("--Db-ref", type=float, required=True)
    p.add_argument("--constant-D", action="store_true")
    p.add_argument("--resume-tip", type=Path)
    p.add_argument("--birth-seed", type=Path)
    p.add_argument("--radius", type=float, default=0.000885)
    p.add_argument("--water", type=float, default=0.025)
    p.add_argument("--cells", type=int, default=4)
    p.add_argument("--case-name")
    p.add_argument("--departure-dt", type=float, default=0.01)
    p.add_argument("--departure-tolerance", type=float, default=1e-10)
    p.add_argument("--target-half-decades", action="store_true")
    p.add_argument("--terminal-entry-tolerance", type=float)
    p.add_argument("--terminal-duration-from-history", action="store_true")
    p.add_argument("--target-root-tolerance", type=float)
    p.add_argument(
        "--b-driver", default="scratchpad/dt_contact_sensitivity_2026_10_08/run_b_032.py"
    )
    args = p.parse_args()
    if not 0 < args.departure_tolerance <= 1e-10:
        raise ValueError("departure refinement may tighten the root tolerance, never relax it")
    if args.terminal_entry_tolerance is not None and not (
        0 < args.terminal_entry_tolerance <= 2e-8
    ):
        raise ValueError("terminal entry may tighten the existing root policy, never relax it")
    if args.target_root_tolerance is not None and not (0 < args.target_root_tolerance <= 1e-8):
        raise ValueError("target roots must stay within the published regional allowance")
    name = f"{args.study}_Db{args.Db_ref/1e-10:g}_{'constant' if args.constant_D else 'Fuller'}"
    if args.case_name:
        name += "_" + args.case_name
    folder = RUNS / name
    folder.mkdir(parents=True, exist_ok=True)
    result = (
        faner(
            folder,
            args.Db_ref,
            args.constant_D,
            args.resume_tip,
            args.radius,
            args.water,
            args.cells,
            args.departure_dt,
            args.departure_tolerance,
            args.target_half_decades,
            args.terminal_entry_tolerance,
            args.terminal_duration_from_history,
            args.target_root_tolerance,
        )
        if args.study == "faner"
        else dt(
            folder, args.Db_ref, args.constant_D, args.resume_tip, args.b_driver, args.birth_seed
        )
    )
    result["controller_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result["scope"] = (
        "from original feed, oil cp(T), original bath and native event/acceptance policies"
    )
    path = folder / f"outcome_{len(list(folder.glob('outcome_*.json')))+1:02d}.json"
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result), flush=True)
    if not result.get("complete"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
