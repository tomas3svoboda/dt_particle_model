"""Reuse the existing manuscript figure and numeric-checking implementation.

The summary explicitly supplies new audited histories. Original records and
generators remain available. Plot layout, measured readings, bath reconstruction
and endpoint checks come from the repository's existing figure generator.
"""

import argparse
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import sys
from release_paths import resolve_record_path as record_path

ROOT = Path(__file__).resolve().parents[3]


def module(path, name):
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def figures(summary, draft):
    generator = ROOT / "paper/01_particle_jfpe/regen/fig_particle_experiments.py"
    existing = module(generator, "temperature_study_figures")
    existing._temperature_adapter_active = True
    existing.record_path = record_path
    existing.PAPER = draft
    existing.FEED_PATH = record_path(summary["dt_feed"])
    existing.CASES = {
        "temperature_study": {
            "segments": list(map(record_path, summary["dt_b_parts"])),
            "dry": record_path(summary["dt_c"]),
            "through_audit": record_path(summary["dt_full_audit"]),
            "b_audit": record_path(summary["dt_full_audit"]),
            "summary": None,
            "K": 0.32,
            "factor": summary.get("liquid_water_mobility_factor", 1000),
        }
    }
    original_reader = existing.rows
    old_faner = (
        ROOT
        / "paper/analysis/results_2026-10-07/faner_assumptions/volume_mean_N4/history/accepted_history.csv"
    )
    old_dt = existing.BASE / "dt_shared_complete_trial/accepted_history.csv"

    def rows(path):
        if path == old_faner:
            path = record_path(summary["faner_history"])
        elif path == old_dt:
            path = record_path(summary["dt_history"])
        return original_reader(path)

    existing.rows = rows
    source = inspect.getsource(existing.read_finite)
    old = 'b_audit["linked_results"][path.relative_to(ROOT).as_posix()]'
    new = 'b_audit["linked_results"].get(str(path), b_audit["linked_results"].get(path.relative_to(ROOT).as_posix()))'
    assert source.count(old) == 1
    exec(compile(source.replace(old, new), str(generator), "exec"), existing.__dict__)
    main_source = inspect.getsource(existing.main)
    # Redirect recorded input paths too, not just the plotted arrays.
    prefix = "    faner_path = (\n"
    start = main_source.index(prefix)
    end = main_source.index("    faner, baseline = rows(faner_path), rows(baseline_path)", start)
    main_source = (
        main_source[:start]
        + '    faner_path = record_path(_temperature_summary["faner_history"])\n'
        + '    baseline_path = record_path(_temperature_summary["dt_history"])\n'
        + main_source[end:]
    )
    existing._temperature_summary = summary
    exec(compile(main_source, str(generator), "exec"), existing.__dict__)
    sys.argv = [str(generator), "--case", "temperature_study"]
    existing.main()
    path = draft / "figures/fig_particle_experiments.sources.json"
    provenance = json.loads(path.read_text())
    provenance["temperature_property_study"] = {
        "summary": "paper/analysis/submission_2026-10-15/study_summary.json",
        "dt_c": summary["dt_c"],
        "faner_history": summary["faner_history"],
    }
    provenance["adapter_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")


def numbers(summary, draft):
    checker = ROOT / "paper/01_particle_jfpe/regen/check_core2_particle_results.py"
    existing = module(checker, "temperature_study_numbers")
    existing._temperature_adapter_active = True
    existing.record_path = record_path
    existing._expected_liquid_factor = summary.get("liquid_water_mobility_factor", 1000)
    existing.MANUSCRIPT = draft
    existing.FANER = record_path(summary["faner_checker_folder"])
    existing.DT032 = record_path(summary["dt_checker_folder"])
    original_path = existing.repo_path
    snapshot = record_path(summary["dt_feed"]).parents[4]

    def repo_path(value):
        path = Path(value.replace("\\", "/"))
        if not path.is_absolute() and str(path).replace("\\", "/").startswith(
            "paper/analysis/temperature_property_runs/"
        ):
            return snapshot / path
        return original_path(value)

    existing.repo_path = repo_path
    existing.main()
    path = draft / "core2_particle_numbers_audit.json"
    provenance = json.loads(path.read_text())
    provenance["temperature_property_study"] = {
        "summary": "paper/analysis/submission_2026-10-15/study_summary.json",
        "dt_c": summary["dt_c"],
        "faner_history": summary["faner_history"],
    }
    provenance["adapter_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--manuscript", type=Path, required=True)
    args = p.parse_args()
    summary = json.loads(args.summary.read_text())
    assert summary["nominal_complete_and_audited"]
    figures(summary, args.manuscript.resolve())
    numbers(summary, args.manuscript.resolve())


if __name__ == "__main__":
    main()
