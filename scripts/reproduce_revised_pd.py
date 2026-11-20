"""Verify, inspect or rerun the explicitly unrefined revised PD reference."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import sys
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "revised_pd"
RUNTIME = BUNDLE / "runtime"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    checked = 0
    for line in (BUNDLE / "MANIFEST.sha256").read_text().splitlines():
        digest, relative = line.split("  ", 1)
        path = BUNDLE / relative
        if not path.is_file() or sha(path) != digest:
            raise ValueError("bundle identity mismatch: " + relative)
        checked += 1
    for relative, digest in json.loads((BUNDLE / "SOURCE_MANIFEST.json").read_text()).items():
        if sha(RUNTIME / relative) != digest:
            raise ValueError("producer source identity mismatch: " + relative)
    with (BUNDLE / "data/accepted_particle_history.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    stats = json.loads((BUNDLE / "data/reported_values.json").read_text())
    previous = -1.0
    pd_source = json.loads((BUNDLE / "data/pd_source.json").read_text())
    for row in rows:
        t = float(row["time_s"])
        h = float(row["hexane_kg_kg_dry"])
        w = float(row["retained_water_kg_kg_dry"])
        film = float(row["external_water_kg_kg_dry"])
        if (
            not all(math.isfinite(x) for x in (t, h, w, film))
            or min(h, w, film) < 0
            or t < previous
        ):
            raise ValueError("inadmissible CSV inventory/time")
        previous = t
        if t <= pd_source["reference"]["residence_s"]:
            if (
                abs(float(row["bath_temperature_c"]) - (pd_source["temperature_k"] - 273.15))
                > 1e-10
                or abs(float(row["bath_y_hexane"]) - pd_source["y_hexane"]) > 1e-12
            ):
                raise ValueError("PD CSV environment differs from operative prescription")
        denominator = 1 + h + w + film
        if abs(float(row["hexane_ppm_wet_basis"]) - 1e6 * h / denominator) > 1e-7:
            raise ValueError("CSV hexane basis mismatch")
        if abs(float(row["moisture_pct_wet_basis"]) - 100 * (w + film) / denominator) > 1e-9:
            raise ValueError("CSV moisture basis mismatch")
    last = rows[-1]
    for column, key in (
        ("temperature_c", "temperature_c"),
        ("hexane_ppm_wet_basis", "hexane_ppm_wet_basis"),
        ("moisture_pct_wet_basis", "moisture_pct_wet_basis"),
    ):
        if abs(float(last[column]) - stats["outlet"][key]) > 1e-9:
            raise ValueError("reported outlet disagrees with CSV")
    print(
        json.dumps(
            {
                "manifest_entries": checked,
                "csv_states": len(rows),
                "recorded_outlet": stats["outlet"],
                "all_csv_wet_basis_values_recomputed": True,
                "full_history_refinement_qualified": False,
                "fresh_full_history_rerun_performed_by_verify": False,
            },
            indent=2,
        )
    )


def native_archive():
    entry = json.loads((BUNDLE / "NATIVE_ARCHIVE.json").read_text())
    target = ROOT / "reruns/downloads" / entry["name"]
    if not target.exists():
        expected = (
            "https://github.com/tomas3svoboda/dt_particle_model/releases/download/v1.5.0-pd-reference/"
            + entry["name"]
        )
        if entry["download_url"] != expected:
            raise ValueError("unexpected native archive URL")
        target.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(expected, timeout=60) as response:
            content = response.read()
        if hashlib.sha256(content).hexdigest() != entry["sha256"]:
            raise ValueError("native download digest mismatch")
        with target.open("xb") as handle:
            handle.write(content)
    if sha(target) != entry["sha256"]:
        raise ValueError("native archive identity mismatch")
    with zipfile.ZipFile(target) as archive:
        expected_names = set(entry["members"]) | {"NATIVE_MANIFEST.json"}
        if set(archive.namelist()) != expected_names:
            raise ValueError("native member inventory differs")
        for name, digest in entry["members"].items():
            if hashlib.sha256(archive.read(name)).hexdigest() != digest:
                raise ValueError("native member digest mismatch: " + name)
    print(json.dumps({"native_archive_verified": str(target), "members": len(entry["members"])}))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        # The unchanged Windows guard sets PYTHONPATH to the public legacy src.
        # Select the isolated exact producer closure explicitly before imports.
        sys.path[:0] = [str(RUNTIME / "src"), str(RUNTIME / "paper/analysis/harness")]
        sys.argv = sys.argv[2:]
        runpy.run_path(sys.argv[0], run_name="__main__")
        return
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--verify-native", action="store_true")
    mode.add_argument("--probe", action="store_true", help="only a 0.1-s cold-feed smoke check")
    mode.add_argument("--run", action="store_true", help="fresh complete reference; may take hours")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    verify()
    if args.verify_native:
        native_archive()
        return
    if not (args.probe or args.run):
        return
    if sys.version_info < (3, 14):
        raise SystemExit("use Python 3.14, the producer environment")
    out = (
        args.out or ROOT / "reruns" / ("pd_probe" if args.probe else "revised_pd_reference")
    ).resolve()
    if out.exists():
        raise SystemExit("refuse to overwrite any existing output")
    stop = out.parent / (out.name + ".STOP")
    if stop.exists():
        raise SystemExit("cooperative stop file already exists")
    env = dict(os.environ)
    for key in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "PYTHONHASHSEED",
    ):
        env[key] = "1"
    env["PYTHONPATH"] = os.pathsep.join(str(RUNTIME / p) for p in ("src", "paper/analysis/harness"))
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        str(RUNTIME / "paper/analysis/harness/run_pd_prescribed_particle_case.py"),
        "--out",
        str(out),
        "--stop-file",
        str(stop),
        "--pd-heating",
        "constant",
        "--cells",
        "64",
        "--time-level",
        "0",
        "--water-film-energy",
        "body",
    ]
    if args.probe:
        command += ["--probe-seconds", "0.1"]
    if os.name == "nt":
        command = [
            sys.executable,
            str(ROOT / "scripts/run_with_resource_budget.py"),
            "--stop-file",
            str(stop),
            "--",
        ] + command
    elif hasattr(os, "sched_getaffinity"):
        cpus = sorted(os.sched_getaffinity(0))
        os.sched_setaffinity(0, cpus[: max(1, len(cpus) // 2)])
        os.nice(10)
    else:
        raise SystemExit("no supported CPU-reserve mechanism on this platform")
    subprocess.run(command, cwd=RUNTIME, env=env, check=True)
    if args.run:
        result = json.loads((out / "result.json").read_text())
        if result.get("status") != "complete" or not result.get("independent_audit_passed"):
            raise SystemExit("reference incomplete; inspect result.json, no completion claimed")
        expected = json.loads((BUNDLE / "data/reported_values.json").read_text())["outlet"]
        differences = {key: result["outlet"][key] - value for key, value in expected.items()}
        print(
            json.dumps(
                {"outlet_difference_from_record": differences, "independent_audit_passed": True}
            )
        )
    else:
        print(json.dumps({"probe_completed": True, "full_history_rerun": False, "out": str(out)}))


if __name__ == "__main__":
    main()
