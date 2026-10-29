"""Verify, extract and rerun the archived histories of this release, and compare them with the originals.

data/final.json pins the archives of this release (SHA-256 of each archive and an inner
EVIDENCE_MANIFEST.sha256 of every member) and records, for every archived case, the exact
runner arguments and the declared environment switches it was produced with.

    python scripts/reproduce_final_reference.py --verify-manifest   # check every file against MANIFEST.sha256
    python scripts/reproduce_final_reference.py --extract-only      # verify and unpack all archives, no run
    python scripts/reproduce_final_reference.py                     # rerun the reference case and compare
    python scripts/reproduce_final_reference.py --case DT_N20_L1    # rerun any archived case and compare

Archives are extracted to their producer-relative, Git-ignored paths (paper/analysis/results_*/);
a modified extracted file is refused, never overwritten.  Fresh outputs go to
reruns/<release>/<case>/ (Git-ignored); an existing output directory is refused.

A case with a declared source variant (data/variants/<case>.json: exact text replacements, e.g. the
Faner radius and tray exposure ratio, or a binary-diffusivity value) is rerun from a scratch copy of
this tree, reruns/<release>/<case>.tree/, with the replacements applied there; each replacement must
occur exactly its declared number of times.  The tree itself is never edited.  The comparison
(scripts/compare_final_case.py) is exact: only wall-clock fields, output paths and source-digest
maps are excused.  Run from the repository root with one numerical worker; the script sets
PYTHONHASHSEED=1, single-threaded BLAS and PYTHONPATH=src:paper/analysis/harness itself.
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import shutil
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "data" / "final.json"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def verify_manifest():
    bad, count = [], 0
    for line in (ROOT / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        count += 1
        path = ROOT / name
        if not path.is_file() or sha256(path.read_bytes()) != digest:
            bad.append(name)
    print(json.dumps({"manifest_entries": count, "mismatched_or_missing": bad}))
    return not bad



def fetch_archive_if_missing(path, digest):
    """Fetch an explicitly requested public release archive and verify it before writing."""
    if path.exists():
        return
    import urllib.request
    catalog = json.loads(FINAL.read_text(encoding="utf-8"))
    entry = next(a for a in catalog["archives"].values() if a["archive"] == path.relative_to(ROOT).as_posix())
    url = entry["download_url"]
    prefix = "https://github.com/tomas3svoboda/dt_particle_model/releases/download/" + catalog["release"] + "/"
    if not url.startswith(prefix) or url != prefix + path.name:
        raise ValueError("archive URL is outside this repository's declared release")
    with urllib.request.urlopen(url, timeout=120) as response:
        content = response.read()
    if sha256(content) != digest:
        raise ValueError("downloaded archive hash mismatch: " + path.name)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(content)

def extract_one(path, digest):
    fetch_archive_if_missing(path, digest)
    if sha256(path.read_bytes()) != digest:
        raise ValueError(f"archive hash mismatch: {path.name}")
    with tarfile.open(path, "r:xz") as archive:
        manifest = archive.extractfile("EVIDENCE_MANIFEST.sha256").read().decode("utf-8")
        expected = {n: d for d, n in (line.split("  ", 1) for line in manifest.splitlines())}
        members = 0
        for member in archive.getmembers():
            if member.name == "EVIDENCE_MANIFEST.sha256":
                continue
            target = (ROOT / member.name).resolve()
            if not member.isfile() or not target.is_relative_to(ROOT.resolve()):
                raise ValueError("archive member is not a regular file inside the repository")
            content = archive.extractfile(member).read()
            if sha256(content) != expected[member.name]:
                raise ValueError(f"member hash mismatch: {member.name}")
            members += 1
            if target.exists():
                if target.read_bytes() != content:
                    raise ValueError(f"modified extracted file; preserve it and retry: {member.name}")
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
    if members != len(expected):
        raise ValueError(f"member count differs from the inner manifest: {path.name}")
    return members


def variant_tree(case, out):
    """A scratch copy of the tree with the case's declared source variant applied; returns its root."""
    spec = json.loads((ROOT / "data" / "variants" / case["source_variant"]["spec"]).read_text(encoding="utf-8"))
    tree = ROOT / (str(out) + ".tree")
    if tree.exists():
        raise SystemExit(f"{tree.relative_to(ROOT)} exists; move it aside first")
    skip = {".git", "data", "reruns", "__pycache__"}
    for top in ROOT.iterdir():
        if top.name in skip:
            continue
        if top.is_dir():
            shutil.copytree(top, tree / top.name,
                            ignore=lambda d, names: [n for n in names if n == "__pycache__" or n.startswith("results_v")])
        else:
            tree.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(top, tree / top.name)
    by_file = {}
    for r in spec["replacements"]:
        by_file.setdefault(r["path"], []).append(r)
    for rel, reps in by_file.items():
        target = tree / rel
        data = target.read_bytes()
        if sha256(data) != reps[0]["file_sha256_shipped"]:
            raise SystemExit(f"{rel}: not the shipped bytes; the variant is declared for the shipped file only")
        crlf = b"\r\n" in data
        text = data.replace(b"\r\n", b"\n").decode("utf-8")
        for r in reps:
            if text.count(r["old"]) != r["count"]:
                raise SystemExit(f"{rel}: replacement occurs {text.count(r['old'])} times, declared {r['count']}")
            text = text.replace(r["old"], r["new"])
        new = text.encode("utf-8")
        if crlf:
            new = new.replace(b"\n", b"\r\n")
        if sha256(new) != reps[0]["file_sha256_with_variant"]:
            raise SystemExit(f"{rel}: the variant does not give the declared bytes")
        target.write_bytes(new)
    print(f"source variant of {case['name']} applied in {tree.relative_to(ROOT)}: {spec['summary']}")
    return tree


def environment(declared):
    env = {k: v for k, v in os.environ.items() if not k.upper().startswith("DTDC_")}
    env.update({"PYTHONHASHSEED": "1", "PYTHONDONTWRITEBYTECODE": "1", "PYTHONUTF8": "1",
                "PYTHONPATH": os.pathsep.join(["src", "paper/analysis/harness"])})
    env.update({k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")})
    env.update(declared)
    return env


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", action="append", help="archived case name (default: the reference case)")
    parser.add_argument("--extract-only", action="store_true")
    parser.add_argument("--verify-manifest", action="store_true")
    args = parser.parse_args()
    if args.verify_manifest:
        return 0 if verify_manifest() else 1
    final = json.loads(FINAL.read_text(encoding="utf-8"))
    archives = final["archives"]
    if args.extract_only:
        counts = {key: extract_one(ROOT / a["archive"], a["sha256"]) for key, a in archives.items()}
        print(json.dumps({"extracted_members": counts}))
        print(f"All {final['release']} archives extracted and verified (Git-ignored paths)")
        return 0
    cases = {c["name"]: c for c in final["cases"]}
    names = args.case or [final["reference_case"]]
    for name in names:
        if name not in cases:
            raise SystemExit(f"unknown case {name}; archived cases: {sorted(cases)}")
        if not cases[name].get("argv"):
            raise SystemExit(f"{name}: no runner arguments are recorded for this case; it cannot be rerun by this script"
                             + (f" ({cases[name]['not_rerunnable']})" if cases[name].get("not_rerunnable") else ""))
    extract_one(ROOT / archives["reference"]["archive"], archives["reference"]["sha256"])
    for key, a in archives.items():
        if any(n in a.get("cases", []) for n in names if n != final["reference_case"]):
            extract_one(ROOT / a["archive"], a["sha256"])
    failures = 0
    for name in names:
        case = cases[name]
        out = Path("reruns") / final["release"] / name
        if (ROOT / out).exists():
            raise SystemExit(f"{out} exists; move it aside first")
        (ROOT / out).parent.mkdir(parents=True, exist_ok=True)
        cwd = ROOT
        if case.get("source_variant"):
            cwd = variant_tree(case, out)
        argv = [sys.executable, "paper/analysis/harness/run_particle_qualification_case.py", *case["argv"],
                "--out", str((ROOT / out).resolve())]
        print("running", " ".join(argv[1:]), "in", cwd, flush=True)
        subprocess.run(argv, cwd=cwd, env=environment(case.get("environment", {})), check=False)
        result = subprocess.run([sys.executable, "scripts/compare_final_case.py", str(out), case["archived_path"],
                                 "--json", str(out) + ".comparison.json"], cwd=ROOT, env=environment({}))
        failures += result.returncode != 0
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
