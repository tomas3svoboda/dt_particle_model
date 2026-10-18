"""Rerun the runtime-151 reference cases and compare them with the archived originals.

The reference cases are the convergence-ladder rows Faner_N8_L1 and DT_N8_L1 of the
frozen solver runtime 151, run here with exactly the ladder's argv (one thread, hash
seed 1).  The archived originals are extracted from data/runtime151_reference.tar.xz
(every member checked against its SHA-256) to their producer-relative, Git-ignored
paths, and each fresh case is compared value by value with compare_runtime151_case.py.

    python scripts/reproduce_runtime151_reference.py                 # both cases
    python scripts/reproduce_runtime151_reference.py --case Faner_N8_L1
    python scripts/reproduce_runtime151_reference.py --extract-only   # every runtime-151 archive, no run

Fresh outputs go to reruns/runtime151/<case>/ (Git-ignored); an existing output
directory is refused, never overwritten.  Expected wall time on one core: Faner about
10 minutes, DT about 20 minutes.
"""

import argparse
import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

REFERENCE = ROOT / "paper/analysis/results_2026-10-17/ladder_151"
COMMON = ["--per-area-qualification", "--dry-local-tolerance", "1e-10", "--conservative-components"]
CASES = {
    "Faner_N8_L1": ["--study", "Faner", "--cells", "8", "--time-level", "1", *COMMON, "--memoize-proof-points"],
    "DT_N8_L1": ["--study", "DT", "--cells", "8", "--time-level", "1", *COMMON,
                 "--event-first-film-depletion", "--eliminate-event-absent-components",
                 "--curvature-corrected-depletion-prediction", "--memoize-proof-points",
                 "--bath", "paper/analysis/inputs/qualification_DT_bath.json"],
}


def extract(everything=False):
    description = json.loads((ROOT / "data/runtime151.json").read_text(encoding="utf-8"))
    keys = [k for k in description if k.endswith("_archive")] if everything else ["reference_archive"]
    for key in keys:
        extract_one(ROOT / description[key]["archive"], description[key]["sha256"])


def extract_one(path, sha256):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == sha256, f"archive hash mismatch: {path.name}"
    with tarfile.open(path, "r:xz") as archive:
        manifest = archive.extractfile("EVIDENCE_MANIFEST.sha256").read().decode("utf-8")
        expected = {n: d for d, n in (line.split("  ", 1) for line in manifest.splitlines())}
        for member in archive.getmembers():
            if member.name == "EVIDENCE_MANIFEST.sha256":
                continue
            target = (ROOT / member.name).resolve()
            if not member.isfile() or not target.is_relative_to(ROOT.resolve()):
                raise ValueError("archive member is not a regular file inside the repository")
            content = archive.extractfile(member).read()
            assert hashlib.sha256(content).hexdigest() == expected[member.name], member.name
            if target.exists():
                if target.read_bytes() != content:
                    raise ValueError(f"modified extracted reference file; preserve it and retry: {member.name}")
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)


def environment():
    env = dict(os.environ)
    env.update({"PYTHONHASHSEED": "1", "PYTHONDONTWRITEBYTECODE": "1", "PYTHONUTF8": "1",
                "PYTHONPATH": os.pathsep.join(["src", "paper/analysis/harness"])})
    env.update({k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")})
    return env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=sorted(CASES), action="append")
    parser.add_argument("--extract-only", action="store_true")
    args = parser.parse_args()
    extract(everything=args.extract_only)
    if args.extract_only:
        print("All runtime-151 archives extracted and verified (Git-ignored paths)")
        return 0
    failures = 0
    for case in args.case or list(CASES):
        out = Path("reruns/runtime151") / case
        if (ROOT / out).exists():
            raise SystemExit(f"{out} exists; move it aside first")
        (ROOT / out).parent.mkdir(parents=True, exist_ok=True)
        argv = [sys.executable, "paper/analysis/harness/run_particle_qualification_case.py", *CASES[case], "--out", str(out)]
        print("running", " ".join(argv[1:]), flush=True)
        subprocess.run(argv, cwd=ROOT, env=environment(), check=False)
        result = subprocess.run([sys.executable, "scripts/compare_runtime151_case.py", str(out), str(REFERENCE / case),
                                 "--json", str(out) + ".comparison.json"], cwd=ROOT, env=environment())
        failures += result.returncode != 0
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
