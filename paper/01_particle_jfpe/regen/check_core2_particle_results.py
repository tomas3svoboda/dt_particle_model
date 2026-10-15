"""Bind current manuscript numbers to accepted shared-model histories.

Writes core2_particle_results.tex (the macros both documents input) and
core2_particle_numbers_audit.json (the macros, the sha256 of every record read
and of this checker) beside the manuscript.

2026-10-09 (REFEREE_AUDIT_8 part 9, section 8.T): the DT macros are read from
the adopted calculation, K_l = 0.32 mol/(m2 s) with total mobility factor 1000
(paper/analysis/results_2026-10-08/dt_contact_0p32/outlet_metrics.json), whose
regime-C result is checked against its recorded sha256 and its independent
through-exit audit; they were previously read from the earlier K_l = 0.10 case
(results_2026-10-07/finite_wetting/completion_summary.json, "dt5") and the
printed values had been entered by hand.  The Faner block is unchanged.
``DTDC_PAPER_RECORDS_ROOT``, when set, names the ``paper/`` folder whose
``analysis/`` records are read (so that a copy of the manuscript folder can be
regenerated against the repository's records); unset, it is this script's
``paper/`` folder, as before.  Record keys in the audit are written relative
to the repository root.

Run from the repository root:
  .venv\\Scripts\\python.exe paper\\01_particle_jfpe\\regen\\check_core2_particle_results.py
"""

import csv
import hashlib
import json
import os
from pathlib import Path

MANUSCRIPT = Path(__file__).resolve().parents[1]
PAPER = (
    Path(os.environ["DTDC_PAPER_RECORDS_ROOT"])
    if os.environ.get("DTDC_PAPER_RECORDS_ROOT")
    else MANUSCRIPT.parent
)
REPO = PAPER.parent
BASE = PAPER / "analysis/results_2026-10-06/repairs"
FANER = PAPER / "analysis/results_2026-10-07/faner_assumptions"
DT032 = PAPER / "analysis/results_2026-10-08/dt_contact_0p32"
BATH_EXIT_S = 1766.5355008060146


def repo_path(stored):
    """A path stored in a record: absolute, or relative to the repository root."""
    from release_paths import resolve_record_path
    return resolve_record_path(stored)


def shown(path):
    """The key under which a record is listed in the audit."""
    try:
        return Path(path).resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return str(path)


def main():
    adoption = MANUSCRIPT / "collaterals/particle_temperature_property_adoption_2026-10-15.json"
    if not adoption.exists():
        adoption = MANUSCRIPT / "collaterals/particle_temperature_property_adoption_2026-10-14.json"
    if adoption.exists() and not globals().get("_temperature_adapter_active", False):
        import sys

        sys.path.insert(0, str(REPO / "paper/analysis/harness"))
        import render_temperature_property_paper as adapter

        manifest = json.loads(adoption.read_text(encoding="utf-8"))
        summary = json.loads(repo_path(manifest["summary"]).read_text(encoding="utf-8"))
        adapter.numbers(summary, MANUSCRIPT)
        return
    faner_path = FANER / "volume_mean_N4/history/accepted_history.csv"
    with faner_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    times = [float(r["t_s"]) for r in rows]
    assert times[0] == 0 and times[-1] == 240
    assert all(b > a for a, b in zip(times, times[1:]))
    assert all(float(r["water_kg_kg_dry"]) > 0 for r in rows)
    fronts = [float(r["front_z"]) for r in rows]
    assert all(0 <= z <= 1 for z in fronts)
    assert all(b <= a for a, b in zip(fronts, fronts[1:]))
    assert float(rows[-1]["front_z"]) == 0

    def crossing(level):
        for left, right in zip(rows, rows[1:]):
            a, b = (float(r["hexane_kg_kg_dry"]) for r in (left, right))
            if a >= level >= b:
                t, u = float(left["t_s"]), float(right["t_s"])
                return t + (u - t) * (a - level) / (a - b)
        raise ValueError("loading crossing missing")

    last = rows[-1]
    dry_time = next(float(r["t_s"]) for r in rows if float(r["front_z"]) == 0)
    values = {
        "FanerDryTime": f"{dry_time:.1f}",
        "FanerCrossTwenty": f"{crossing(.20):.1f}",
        "FanerFallingInterval": f"{crossing(.10) - crossing(.20):.1f}",
        "FanerExitHexane": f"{float(last['hexane_kg_kg_dry']):.3g}",
        "FanerExitWater": f"{float(last['water_kg_kg_dry']):.3g}",
        "FanerExitTemperature": f"{float(last['mean_temperature_k']) - 273.15:.1f}",
        "FanerFreeEnd": f"{max(float(r['t_s']) for r in rows if r['stage'] == 'attached'):.1f}",
        "FanerFallingShortfall": f"{100 * (1 - (crossing(.10) - crossing(.20)) / 25.350441949123166):.0f}",
    }
    outcomes_path = FANER / "outcomes.json"
    outcome = json.loads(outcomes_path.read_text())["volume_mean_N4"]
    assert outcome["complete"] and outcome["export_passed"]
    source_path = repo_path(outcome["tip"])
    betas = []
    lineage = []
    while True:
        record = json.loads(source_path.read_text())
        assert record["sources_unchanged_during_run"]
        lineage.append(source_path)
        betas.extend(abs(r["beta"]) for r in record.get("feed_rows", []) if r["beta"] is not None)
        if record.get("surface"):
            betas.append(abs(record["surface"]["beta"]))
        for attempt in record.get("attempts", []):
            if not attempt.get("accepted"):
                continue
            surface = attempt.get("surface") or attempt.get("audit", {}).get("surface_water_node")
            if surface:
                betas.append(abs(surface["beta"]))
        parent = record.get("accepted_parent")
        if not parent:
            assert record["inputs"]["radius"] == 0.000885
            assert record["inputs"]["cells"] == 4
            assert record["inputs"]["water"] == 0.025
            assert record["inputs"]["coefficient_case"] == "db_max_dapp_max"
            break
        source_path = repo_path(parent["result"])
        assert hashlib.sha256(source_path.read_bytes()).hexdigest() == parent["result_sha256"]
    values["FanerMaxBeta"] = f"{max(betas):.3f}"
    for left, right in zip(rows, rows[1:]):
        a, b = float(left["t_s"]), float(right["t_s"])
        if a <= 235.058 <= b:
            weight = (235.058 - a) / (b - a)
            temperature = float(left["mean_temperature_k"]) + weight * (
                float(right["mean_temperature_k"]) - float(left["mean_temperature_k"])
            )
            values["FanerLastHolderClockTemperature"] = f"{temperature - 273.15:.1f}"
            break
    else:
        raise ValueError("last measured holder clock is outside the computed history")
    summary_path = FANER / "volume_mean_N4/history/summary.json"
    summary = json.loads(summary_path.read_text())
    assert summary["complete_falling_rate_experiment"]
    assert (
        max(summary["independent_whole_history_external_ledger_maximum_normalized_residuals"])
        <= 4e-10
    )

    # The adopted DT calculation: K_l = 0.32 mol/(m2 s), total mobility factor 1000.
    adopted_path = DT032 / "outlet_metrics.json"
    adopted = [
        e for e in json.loads(adopted_path.read_text(encoding="utf-8")) if e["K_l_mol_m2_s"] == 0.32
    ]
    assert len(adopted) == 1
    adopted = adopted[0]
    assert adopted["external_water_kg_kg_dry"] == 0
    assert abs(adopted["time_s"] - BATH_EXIT_S) < 1e-9
    result_path = repo_path(adopted["result"])
    assert hashlib.sha256(result_path.read_bytes()).hexdigest() == adopted["result_sha256"]
    regime_c = json.loads(result_path.read_text(encoding="utf-8"))
    assert regime_c["finite_liquid_contact_conductance_mol_m2_s"] == 0.32
    assert regime_c["additional_parallel_path_total_to_base_mobility_factor"] == globals().get(
        "_expected_liquid_factor", 1000.0
    )
    audit_path = repo_path(adopted["full_audit"])
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    assert audit["passed"] and audit["bath_exit_reached"] and audit["regime_b_completed"]
    assert audit["dry_result_sha256"] == adopted["result_sha256"]
    assert audit["declared_tolerance"] == 1e-10
    assert max(audit["maximum_normalized_W_H_U"]) <= 1e-10
    assert audit["maximum_normalized_W_H_U"] == adopted["maximum_normalized_W_H_U"]
    values.update(
        {
            "DTFiniteHexaneDry": f"{adopted['hexane_mg_kg_dry']:.0f}",
            "DTFiniteHexaneWet": f"{adopted['hexane_ppm_wet_basis']:.0f}",
            "DTFiniteMoisture": f"{adopted['moisture_percent_wet_basis']:.1f}",
            "DTFiniteTemperature": f"{adopted['mean_temperature_c']:.1f}",
        }
    )
    macros = "% Generated by regen/check_core2_particle_results.py.\n"
    macros += "".join(f"\\newcommand{{\\{key}}}{{{value}}}\n" for key, value in values.items())
    (MANUSCRIPT / "core2_particle_results.tex").write_text(macros, encoding="utf-8")
    report = {
        "passed": True,
        "scientific_submission_release": False,
        "macros": values,
        "input_sha256": {
            shown(p): hashlib.sha256(Path(p).read_bytes()).hexdigest()
            for p in (
                faner_path,
                summary_path,
                adopted_path,
                result_path,
                audit_path,
                *lineage,
            )
        },
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (MANUSCRIPT / "core2_particle_numbers_audit.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report["macros"]))


if __name__ == "__main__":
    main()
