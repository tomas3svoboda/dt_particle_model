"""Check the current released histories and regenerate the paper outputs.

Runs serially. Native constitutive inventories and exterior flow ledgers are
re-evaluated, not inferred from the solver's convergence flag. Existing raw
JSON and accepted checkpoints are never rewritten.
"""

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "paper/analysis/harness"
PAPER = ROOT / "paper/01_particle_jfpe"
STUDY = ROOT / "paper/analysis/submission_2026-10-15"
sys.path[:0] = [str(ROOT / "src"), str(HARNESS)]

from release_paths import resolve_record_path as resolve  # noqa: E402
import particle_temperature_properties as properties  # noqa: E402


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_manifest():
    manifest = ROOT / "MANIFEST.sha256"
    if not manifest.exists():
        return
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        assert sha(ROOT / name) == digest, name
    print("SHA-256 manifest: passed", flush=True)


def audit_dt(summary, initialize=False):
    import audit_finite_wetting_b as native

    results = {}
    for name, case in summary["dt_cases"].items():
        dry = resolve(case["result"])
        record = json.loads(dry.read_text(encoding="utf-8"))
        parent = resolve(record["parent_checkpoint"])
        part = parent.with_name(parent.name.replace(".accepted.pkl", ".json"))
        assert sha(parent) == record["parent_sha256"]
        law = case["property_variant"]
        with properties.installed(
            law["D_b_reference_m2_s"],
            temperature_law=law["binary_temperature_law"] != "constant",
        ):
            report = native.audit(part, [part], dry)
        assert report["passed"] and report["regime_b_completed"]
        assert report["dry_result_sha256"] == case["result_sha256"]
        assert report["maximum_normalized_W_H_U"] == case["maximum_normalized_W_H_U"]
        # Canonical public holder keys; the raw private audit remains untouched.
        report["linked_results"] = {
            part.relative_to(ROOT).as_posix(): sha(part)
        }
        if initialize:
            out = resolve(case["full_audit"])
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(report, indent=2), encoding="utf-8")
        results[name] = {"passed": True, "maximum_normalized_W_H_U": report["maximum_normalized_W_H_U"]}
        print(f"Native DT full-history audit: {name} passed", flush=True)
    return results


def audit_faner(summary):
    exporter = HARNESS / "export_core2_faner_history.py"
    results = {}
    with tempfile.TemporaryDirectory(prefix="particle_faner_check_") as temporary:
        for name, case in summary["faner_cases"].items():
            out = Path(temporary) / name
            argv = [sys.executable, str(HARNESS / "run_temperature_property_driver.py"),
                    "--Db-ref", str(int(name.split("Db")[1].split("_")[0]) * 1e-10),
                    "--driver", str(exporter), "--", "--tip", str(resolve(case["tip"])),
                    "--out", str(out), "--ledger-tolerance", "4e-10"]
            result = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
            assert result.returncode == 0, result.stdout + result.stderr
            original = json.loads(resolve(case["audit"]).read_text(encoding="utf-8"))
            actual = json.loads((out / "summary.json").read_text(encoding="utf-8"))
            for key in ("last_accepted", "row_count", "independent_whole_history_external_ledger_maximum_normalized_residuals"):
                assert actual[key] == original[key], (name, key)
            assert (out / "accepted_history.csv").read_bytes() == resolve(case["history"]).read_bytes()
            results[name] = {"passed": True, "maximum_normalized_H_W_U": actual["independent_whole_history_external_ledger_maximum_normalized_residuals"]}
            print(f"Native Faner lineage and exterior-ledger audit: {name} passed", flush=True)
    return results


def render(summary):
    import render_temperature_property_paper as adapter

    adapter.figures(summary, PAPER)
    adapter.numbers(summary, PAPER)
    # Restore canonical slashes so this holder is usable on Windows and POSIX.
    for name in ("core2_particle_numbers_audit.json", "figures/fig_particle_experiments.sources.json"):
        path = PAPER / name
        data = json.loads(path.read_text(encoding="utf-8"))
        for key in ("input_sha256", "source_sha256"):
            if key in data:
                data[key] = {str(k).replace("\\", "/"): v for k, v in data[key].items()}
        path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    graphical = PAPER / "regen/fig_graphical_abstract_core2.py"
    subprocess.run([sys.executable, str(graphical)], cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initialize", action="store_true", help="Create derived public audit holders at export time")
    parser.add_argument("--render", action="store_true", help="Regenerate figures, numeric macros and graphical abstract")
    parser.add_argument("--skip-faner", action="store_true", help="Check only DT (not the complete submission check)")
    args = parser.parse_args()
    os.chdir(ROOT)
    import sys
    sys.path.insert(0, str(ROOT / "paper/analysis/harness"))
    from submission_evidence import ensure_evidence
    ensure_evidence()
    summary = json.loads((STUDY / "study_summary.json").read_text(encoding="utf-8"))
    if not args.initialize:
        verify_manifest()
    results = {"dt": audit_dt(summary, args.initialize)}
    if not args.skip_faner:
        results["faner"] = audit_faner(summary)
    if args.render:
        render(summary)
    path = ROOT / "reruns/submission_check.json"
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("Submission checks passed; results in reruns/submission_check.json", flush=True)


if __name__ == "__main__":
    main()
