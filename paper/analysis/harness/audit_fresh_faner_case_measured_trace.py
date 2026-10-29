"""Independent Faner original-feed W/H/U audit using native EOS and film rates, for a case run with
kernel-envelope FANER_MEASURED_VAPOUR_TRACE armed.

The campaign's Faner auditor (paper/analysis/results_2026-10-16/solver_qualification/
audit_fresh_Faner_case_51.py, sha256 e1a96083... as checked out, unchanged and still the auditor of
every constant-bath case) with one change, the measurement variant's (record
docs/GT_PS2_PAPER1_FANER_MEASURED_VAPOUR_TRACE_AND_DT_BATH_CONSTANTS_2026-10-28.md, its variant
auditor audit_fresh_Faner_case_51_measured_trace.py): each feed row's exterior rates are
reconstructed with the measured trace at the row's own clock (faner_measured_trace.boundary_at).
Built for landing by the three-knob build (docs/GT_PS2_PAPER1_THREE_HARNESS_KNOBS_BUILD_2026-10-29.md):
a case whose birth does not carry the trace stamp, or carries a trace other than the runtime's
declared one, is refused. The radial rows are audited exactly as the campaign's auditor does.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import pickle
import sys
from types import SimpleNamespace

parser = argparse.ArgumentParser()
parser.add_argument("--runtime", type=Path, required=True)
parser.add_argument("--case", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
runtime, folder, target = (p.resolve() for p in (args.runtime, args.case, args.out))
if target.exists():
    raise ValueError("refusing to overwrite independent evidence")
sys.path[:0] = [str(runtime / "src"), str(runtime / "paper/analysis/harness")]
import audit_particle_qualification_case as shared  # noqa: E402
import experiment_core2_faner as faner  # noqa: E402
import faner_measured_trace as measured_trace  # noqa: E402
from dtdc_simulator.core2.particle import grid as spherical_grid  # noqa: E402

core = shared.core
record = json.loads((folder / "case.json").read_text())
if (
    record["study"] != "Faner"
    or not record["complete_history"]
    or not record["sources_unchanged_during_run"]
):
    raise ValueError("completed unchanged-source Faner original-feed case required")
pins = {p: h for p, h in record["source_sha256_at_entry"].items() if p.startswith("src/")}
if not pins or any(
    hashlib.sha256((runtime / p).read_bytes()).hexdigest() != h for p, h in pins.items()
):
    raise ValueError("native audit source differs from producer")
birth = json.loads((folder / "birth.json").read_text())
if birth.get("faner_measured_vapour_trace") != measured_trace.declared_row():
    raise ValueError(
        "a measured-vapour-trace Faner case required (the birth's trace stamp equal to the declared trace)"
    )
settings = SimpleNamespace(**birth["inputs"])
checkpoint = pickle.loads((folder / "accepted.pkl").read_bytes())
config = checkpoint["config"]
with (
    faner.declared_faner_film_multiplier(settings.film_multiplier),
    faner.pure_bulk_domain_adapter(),
):
    march = faner.FanerAttachedMarch(settings)
    march.measured_trace = True  # the case's trace, whatever this process's knob
    area = march.area

    def feed_storage(row):
        critical = core.sphere.critical_loading(row["temperature_k"], march.params)
        total = (
            settings.initial_hexane
            if row["t_s"] == 0
            else critical + row["attached_hexane_kg_kg_dry"]
        )
        reported = row["hexane_kg_kg_dry"]
        if abs(total - reported) > 4 * max(math.ulp(total), math.ulp(reported)):
            raise ValueError(
                "primary attached loading and derived total exceed their serialization-roundoff budget"
            )
        state = core.sphere.initialize_qualified_feed(
            march.params,
            total_hexane_loading=total,
            temperature_k=row["temperature_k"],
            retained_water_loading=row["retained_water_kg_kg_dry"],
            oil_fraction=settings.oil,
        )
        return [
            state.retained_water_mass_kg,
            state.total_hexane_mass_kg,
            state.total_internal_energy_j,
        ]

    initial = feed_storage(birth["feed_rows"][0])
    balance = shared.OriginalFeedBalance(
        initial,
        [initial[0], initial[1], checkpoint["carry"]["import_capacity_scale_j"]],
        config.dry.ledger_tolerance,
    )
    for row in birth["feed_rows"][1:]:
        march.trace_clock = row["t_s"]
        ts = march.interface_temperature
        boundary = measured_trace.boundary_at(march.boundary, row["t_s"])
        reference = core.ct.prepare_reduced_film_surface_state(
            boundary.film_boundary,
            cell_temperature_k=row["temperature_k"],
            pressure_pa=core.PRESSURE_PA,
            particle_radius_m=march.params.particle_radius_m,
            outer_half_cell_distance_m=0.5 * march.params.particle_radius_m,
            particle_thermal_conductivity_w_m_k=march.wet.conductivity,
            binary_gas_interaction_k_wh=boundary.film_boundary.binary_gas_interaction_k_wh,
        )
        film = core.ct._reduced_film_coefficients_at_surface_temperature(reference, ts)
        capacities = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
            surface_temperature_k=ts,
            bulk_temperature_k=boundary.temperature_k,
            pressure_pa=core.PRESSURE_PA,
            y_hexane=1.0,
        )
        # Native Ackermann beta=(j_w cp_w+j_h cp_h)/h. Before attached-film
        # depletion j_w=0 and both endpoint gases are pure hexane.
        total = row["beta"] * film.heat_transfer_w_m2_k / capacities.hexane_j_kg_k
        rates = march.rates(row["temperature_k"], total)
        if abs(rates["beta"] - row["beta"]) > 1e-12 or rates["water_out_kg_s"] != 0:
            raise ValueError("native pure-hexane feed-rate inverse does not close")
        balance.accept(
            row["t_s"],
            feed_storage(row),
            [rates[k] for k in ("water_out_kg_s", "hexane_out_kg_s", "energy_out_w")],
            "fresh Faner feed; native EOS and reconstructed exterior energy/component rates",
        )


def accept_radial(transport, surface, kind):
    if "geometry" in transport:
        inv, _ = shared.independent.inventory(transport, config)
    else:
        grid = spherical_grid._from_faces(tuple(transport["grid"]["faces"]), transport["grid"]["R"])
        if tuple(transport["grid"]["volumes"]) != grid.volumes:
            raise ValueError("native dry geometry differs")
        cells = core.ct._equilibrium_states(
            transport["temperatures_k"], transport["y_hexane"], config.dry
        )
        native = core.ct._inventories(grid, cells)
        inv = [native.total_water_mol, native.total_hexane_mol, native.total_energy_j]
    film = area * surface["liquid_water_kg_m2"]
    actual = [
        inv[0] * core.wa.M + film,
        inv[1] * core.hx.M,
        inv[2]
        + film * core.wa.state_Tp(surface["temperature_k"], core.PRESSURE_PA, "liquid").u_mass,
    ]
    balance.accept(
        transport["time_s"],
        actual,
        [
            area * surface[k]
            for k in (
                "external_water_flux_kg_m2_s",
                "external_hexane_flux_kg_m2_s",
                "external_energy_flux_w_m2",
            )
        ],
        kind,
    )


native_birth = next(
    a["audit"] for a in birth["attempts"] if a.get("stage") == "birth" and a.get("accepted")
)
accept_radial(
    native_birth["after_transport"], native_birth["surface_water_node"], "native Faner birth"
)
for line in (folder / "front.history.jsonl").read_text().splitlines():
    row = json.loads(line)
    if row["kind"] not in ("initial", "refusal"):
        accept_radial(row["after_transport"], row["surface"], "front " + row["kind"])
endpoint = json.loads((folder / "endpoint.json").read_text())
accept_radial(
    endpoint["after_transport"], endpoint["surface"], "native endpoint and base-face restoration"
)
for line in (folder / "dry.history.jsonl").read_text().splitlines():
    row = json.loads(line)
    if row["kind"] not in ("initial", "refusal"):
        accept_radial(row["after_transport"], row["surface"], "dry " + row["kind"])
if balance.clock != 240.0:
    raise ValueError("prescribed Faner outlet time missing")
report = {
    "passed": True,
    "complete_history": True,
    "native_source_pins_checked": len(pins),
    "physical_rows_audited": len(balance.rows),
    "maximum_independent_normalized_W_H_U": balance.maximum,
    "original_feed_W_H_U": balance.initial,
    "original_scales": balance.scales,
    "outlet_time_s": balance.clock,
    "case_json_sha256": hashlib.sha256((folder / "case.json").read_bytes()).hexdigest(),
    "auditor_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "faner_measured_vapour_trace": birth["faner_measured_vapour_trace"],
    "feed_reconstruction": "declared initial total, then native critical+primary attached loading; validate derived total within four ulps; never clamp subcritical totals",
    "scope": "independent balance audit; native EOS and feed-film-rate reconstruction; radial archived exterior rates; not independent physical validation",
    "convergence_established": False,
    "submission_released": False,
    "rows": balance.rows,
}
with target.open("x") as stream:
    json.dump(report, stream, indent=2)
print(json.dumps({k: v for k, v in report.items() if k != "rows"}), flush=True)
