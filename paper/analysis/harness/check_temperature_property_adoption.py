"""Check the current property law, source identities and manuscript adoption."""

import argparse
import hashlib
import json
from pathlib import Path

from particle_temperature_properties import declaration

ROOT = Path(__file__).resolve().parents[3]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(summary_path, manuscript):
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["nominal_complete_and_audited"]
    for path, expected in summary["source_sha256"].items():
        assert sha(path) == expected, path
    expected = json.loads(json.dumps(declaration(4e-10, True, True)))
    law_hash = sha(ROOT / "paper/analysis/harness/particle_temperature_properties.py")
    nominal = summary["dt_cases"]["dt_Db4_Fuller"]
    assert nominal["property_variant"] == expected
    factor = summary.get("liquid_water_mobility_factor", 1000)
    for path in summary["dt_b_parts"]:
        record = json.loads(Path(path).read_text(encoding="utf-8"))
        assert record["shell_mobility_factor"] == factor
        assert record["contact_conductance_mol_m2_s"] == 0.32
    dry = json.loads(Path(summary["dt_c"]).read_text(encoding="utf-8"))
    assert dry["additional_parallel_path_total_to_base_mobility_factor"] == factor
    assert dry["finite_liquid_contact_conductance_mol_m2_s"] == 0.32
    checked = []
    for path in [
        summary["dt_feed"],
        *summary["dt_b_parts"],
        summary["dt_c"],
        summary["faner_cases"]["faner_Db4_Fuller"]["tip"],
    ]:
        record = json.loads(Path(path).read_text(encoding="utf-8"))
        assert record["temperature_property_variant"] == expected
        assert record["property_variant_sources_unchanged"]
        hashes = record["property_variant_source_sha256"]
        assert any(
            Path(p.replace("\\", "/")).name == "particle_temperature_properties.py" and h == law_hash
            for p, h in hashes.items()
        )
        checked.append({"path": str(path), "sha256": sha(path)})
    assert nominal["maximum_normalized_W_H_U"] and max(nominal["maximum_normalized_W_H_U"]) <= 1e-10
    text = (manuscript / "sec_formulation.tex").read_text(encoding="utf-8")
    for token in [
        r"4.0\times10^{-10}",
        "1024+3T",
        "1.75",
        "garciarojas2013thermophysical",
        "fuller1966new",
        "unvalidated",
    ]:
        assert token in text, token
    if factor == 100:
        for token in ["99", r"8.39\times10^{-5}", "hemmingsen2008water"]:
            assert token in text, token
    for name in ("main.tex", "supplementary_submission.tex"):
        assert "Henke" not in (manuscript / name).read_text(encoding="utf-8")
    numbers = json.loads(
        (manuscript / "core2_particle_numbers_audit.json").read_text(encoding="utf-8")
    )
    assert numbers["passed"]
    for path, expected_hash in numbers["input_sha256"].items():
        candidate = Path(path)
        assert sha(candidate if candidate.is_absolute() else ROOT / candidate) == expected_hash
    macros = numbers["macros"]
    assert macros["DTFiniteHexaneWet"] == f"{nominal['hexane_ppm_wet_basis']:.0f}"
    assert macros["DTFiniteMoisture"] == f"{nominal['moisture_percent_wet_basis']:.1f}"
    figures = json.loads(
        (manuscript / "figures/fig_particle_experiments.sources.json").read_text(encoding="utf-8")
    )
    for path, expected_hash in figures["source_sha256"].items():
        assert sha(ROOT / path) == expected_hash, path
    assert figures["temperature_property_study"]["dt_c"] == summary["dt_c"]
    assert figures["temperature_property_study"]["faner_history"] == summary["faner_history"]
    figure_hashes = {
        name + "." + extension: sha(manuscript / "figures" / (name + "." + extension))
        for name in (
            "fig_dt_particle_history_core2",
            "fig_faner_particle_history_core2",
            "fig_faner_supporting_states",
        )
        for extension in ("pdf", "png")
    }
    return {
        "passed": True,
        "summary": summary_path.resolve().relative_to(ROOT).as_posix(),
        "summary_sha256": sha(summary_path),
        "property_law": expected,
        "property_law_sha256": law_hash,
        "accepted_nominal_records": checked,
        "checker_sha256": sha(__file__),
        "generated_figure_sha256": figure_hashes,
        "source_coefficient_convention_conversion_completed": False,
        "native_core_defaults_changed": False,
        "liquid_water_mobility_factor": factor,
        "liquid_contact_conductance_mol_m2_s": 0.32,
        "property_variant_flags_scope": "Flags describe the property adapter. Caller numerical recoveries are recorded separately in native result sidecars; source-pinned controllers and their refused attempts are retained with each history.",
        "native_constant_cp_fields": "Retained numerical reference fields; physical cp(T) and its integrated energy are supplied by the declared property adapter.",
        "physically_qualifying": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--manuscript", type=Path, required=True)
    parser.add_argument(
        "--record-name", default="particle_temperature_property_adoption_2026-10-14.json"
    )
    args = parser.parse_args()
    record = check(args.summary.resolve(), args.manuscript.resolve())
    assert Path(args.record_name).name == args.record_name
    path = args.manuscript / "collaterals" / args.record_name
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(f"PASS: {path}")


if __name__ == "__main__":
    main()
