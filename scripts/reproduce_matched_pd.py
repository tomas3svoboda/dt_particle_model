"""Verify article data, extract native records, or rerun one prescribed-DT case."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import stat
import subprocess
import sys
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT/"matched_pd"
RUNTIME = BUNDLE/"runtime"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    entries = (BUNDLE/"MANIFEST.sha256").read_text().splitlines()
    for entry in entries:
        digest, name = entry.split("  ", 1)
        path = (BUNDLE/name).resolve()
        if not path.is_relative_to(BUNDLE.resolve()) or not path.is_file() or sha(path) != digest:
            raise ValueError("bundle identity mismatch: "+name)
    for name, digest in read(BUNDLE/"SOURCE_MANIFEST.json").items():
        if sha(RUNTIME/name) != digest:
            raise ValueError("numerical source identity mismatch: "+name)
    with (BUNDLE/"data/accepted_particle_history.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    previous = -1.0
    for row in rows:
        t, h, w, film = (float(row[key]) for key in (
            "time_s", "hexane_kg_kg_dry", "retained_water_kg_kg_dry", "external_water_kg_kg_dry"))
        if not all(math.isfinite(x) for x in (t,h,w,film)) or min(h,w,film) < 0 or t <= previous:
            raise ValueError("inadmissible CSV state or clock")
        previous = t
        denominator = 1+h+w+film
        if abs(float(row["hexane_ppm_wet_basis"])-1e6*h/denominator) > 1e-7:
            raise ValueError("CSV hexane mass basis mismatch")
        if abs(float(row["moisture_pct_wet_basis"])-100*(w+film)/denominator) > 1e-9:
            raise ValueError("CSV moisture mass basis mismatch")
    reported = read(BUNDLE/"data/reported_values.json")["outlet"]
    for key, value in reported.items():
        if abs(float(rows[-1][key])-value) > 1e-9:
            raise ValueError("reported outlet differs from accepted CSV")
    audit = read(BUNDLE/"data/independent_whole_history_audit.json")
    if not (audit["passed"] and audit["complete_history"]):
        raise ValueError("archived reference audit has not passed")
    provenance = read(BUNDLE/"PROVENANCE.json")
    if sha(BUNDLE/provenance["operative_bath"]) != provenance["operative_bath_sha256"]:
        raise ValueError("operative bath identity mismatch")
    print(json.dumps(dict(manifest_files=len(entries), accepted_csv_states=len(rows),
        recorded_outlet=reported, csv_mass_bases_verified=True, archived_whole_history_audit_passed=True), indent=2))


def extract_native(out):
    record = read(BUNDLE/"NATIVE_ARCHIVE.json")
    tag = read(BUNDLE/"PROVENANCE.json")["release"]
    expected_url = "https://github.com/tomas3svoboda/dt_particle_model/releases/download/"+tag+"/"+record["name"]
    if record["download_url"] != expected_url:
        raise ValueError("unexpected native archive URL")
    archive = ROOT/"reruns/downloads"/record["name"]
    if not archive.exists():
        with urllib.request.urlopen(expected_url, timeout=60) as response:
            content = response.read()
        if hashlib.sha256(content).hexdigest() != record["sha256"]:
            raise ValueError("native archive download hash mismatch")
        archive.parent.mkdir(parents=True, exist_ok=True)
        with archive.open("xb") as stream:
            stream.write(content)
    if sha(archive) != record["sha256"]:
        raise ValueError("native archive identity mismatch")
    if out.exists():
        raise ValueError("use a new native-record destination")
    with zipfile.ZipFile(archive) as handle:
        if set(handle.namelist()) != set(record["members"])|{"NATIVE_MANIFEST.json"}:
            raise ValueError("native archive member list differs")
        for name, digest in record["members"].items():
            info = handle.getinfo(name)
            target = (out/name).resolve()
            if not target.is_relative_to(out) or stat.S_ISLNK(info.external_attr >> 16):
                raise ValueError("invalid native record path")
            if hashlib.sha256(handle.read(name)).hexdigest() != digest:
                raise ValueError("native record identity mismatch: "+name)
        out.mkdir(parents=True)
        for name in record["members"]:
            target = out/name
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("xb") as stream:
                stream.write(handle.read(name))
    print(json.dumps(dict(native_records=str(out), members=len(record["members"]))))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        sys.path[:0] = [str(RUNTIME/"src"), str(RUNTIME/"paper/analysis/harness")]
        sys.argv = sys.argv[2:]
        runpy.run_path(sys.argv[0], run_name="__main__")
        return
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--extract-native", action="store_true")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--case", default="fresh_n64_a16_l1")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--compiled-eos", action="store_true")
    args = parser.parse_args()
    verify()
    if args.extract_native:
        extract_native((args.out or RUNTIME/"scratchpad/archived_reference").resolve())
        return
    if not args.run:
        return
    if sys.version_info[:2] != (3,14):
        raise SystemExit("Use Python 3.14, the numerical producer environment")
    cases = {c["name"]:c for c in read(BUNDLE/"data/cases.json")["cases"]}
    if args.case not in cases:
        raise SystemExit("Unknown case; available: "+", ".join(cases))
    case = cases[args.case]
    out = (args.out or RUNTIME/"scratchpad/reproductions"/args.case).resolve()
    if not out.is_relative_to((RUNTIME/"scratchpad").resolve()) or out.exists():
        raise SystemExit("Use a new output folder inside matched_pd/runtime/scratchpad")
    stop = out.with_name(out.name+".STOP")
    if stop.exists():
        raise SystemExit("Cooperative stop marker already exists")
    env = dict(os.environ)
    for name in ("PYTHONHASHSEED", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[name] = "1"
    command = [sys.executable, str(Path(__file__).resolve()), "--worker",
        str(RUNTIME/"paper/analysis/harness/run_matched_particle_case.py"),
        "--out", str(out), "--stop-file", str(stop), "--cells", str(case["cells"]),
        "--main-shells", str(case["main_shells"]), "--time-level", str(case["time_level"]),
        "--geometry", case["geometry"], "--bath", str(RUNTIME/"inputs/publication_bath.json")]
    if args.compiled_eos:
        command.append("--compiled-eos")
    if os.name == "nt":
        command = [sys.executable, str(ROOT/"scripts/run_with_resource_budget.py"),
            "--stop-file", str(stop), "--"]+command
    elif hasattr(os, "sched_getaffinity"):
        cpus = sorted(os.sched_getaffinity(0))
        os.sched_setaffinity(0, cpus[:max(1,len(cpus)//2)])
        os.nice(10)
    else:
        raise SystemExit("No supported CPU-reserve mechanism on this platform")
    subprocess.run(command, cwd=RUNTIME, env=env, check=True)
    result = read(out/"production.json")
    if not (result.get("complete_history") and result.get("independent_whole_history_passed")
            and result.get("extension_sources_unchanged")):
        raise SystemExit("History incomplete or unaudited; inspect production.json")
    last = read(out/"independent_whole_history_audit.json")
    print(json.dumps(dict(case=args.case, complete_history=True, independent_audit=last,
        output=str(out)), indent=2))


if __name__ == "__main__":
    main()
