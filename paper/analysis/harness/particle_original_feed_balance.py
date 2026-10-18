"""Carry the original-feed conservation datum and integrals across C activation."""

import hashlib
import json
import math
import pickle
from types import SimpleNamespace

import audit_particle_qualification_case as shared
import experiment_core2_faner as faner

core = shared.core


def prepare_prefix(folder, study, bath_path, live_source_hashes):
    if study == "DT":
        return shared.audit(
            folder,
            bath_path,
            accepted_prefix=True,
            live_source_hashes=live_source_hashes,
            return_balance=True,
            until_core_endpoint=True,
        )
    if study != "Faner":
        raise ValueError("unsupported original-feed particle study")
    record = json.loads((folder / "case.json").read_text())
    shared.validate_producer_sources(record, live_source_hashes)
    if not record["fresh_original_feed"] or not record["single_runtime_source_generation"]:
        raise ValueError("original-feed single-runtime prefix required")
    native_pins = {
        p: h for p, h in record["source_sha256_at_entry"].items() if p.startswith("src/")
    }
    if not native_pins or any(
        hashlib.sha256((faner.REPO_ROOT / p).read_bytes()).hexdigest() != h
        for p, h in native_pins.items()
    ):
        raise ValueError("native prefix source differs from producer")
    birth = json.loads((folder / "birth.json").read_text())
    checkpoint = pickle.loads((folder / "birth.accepted.pkl").read_bytes())
    config = checkpoint["config"]
    settings = SimpleNamespace(**birth["inputs"])
    with (
        faner.declared_faner_film_multiplier(settings.film_multiplier),
        faner.pure_bulk_domain_adapter(),
    ):
        march = faner.FanerAttachedMarch(settings)
        area = march.area

        def storage(row):
            critical = core.sphere.critical_loading(row["temperature_k"], march.params)
            total = (
                settings.initial_hexane
                if row["t_s"] == 0
                else critical + row["attached_hexane_kg_kg_dry"]
            )
            reported = row["hexane_kg_kg_dry"]
            if abs(total - reported) > 4 * max(math.ulp(total), math.ulp(reported)):
                raise ValueError("primary/derived feed inventory codec discrepancy")
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

        initial = storage(birth["feed_rows"][0])
        balance = shared.OriginalFeedBalance(
            initial,
            [initial[0], initial[1], checkpoint["carry"]["import_capacity_scale_j"]],
            config.dry.ledger_tolerance,
        )
        for row in birth["feed_rows"][1:]:
            ts = march.interface_temperature
            boundary = march.boundary
            ref = core.ct.prepare_reduced_film_surface_state(
                boundary.film_boundary,
                cell_temperature_k=row["temperature_k"],
                pressure_pa=core.PRESSURE_PA,
                particle_radius_m=march.params.particle_radius_m,
                outer_half_cell_distance_m=0.5 * march.params.particle_radius_m,
                particle_thermal_conductivity_w_m_k=march.wet.conductivity,
                binary_gas_interaction_k_wh=boundary.film_boundary.binary_gas_interaction_k_wh,
            )
            film = core.ct._reduced_film_coefficients_at_surface_temperature(ref, ts)
            cp = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
                surface_temperature_k=ts,
                bulk_temperature_k=boundary.temperature_k,
                pressure_pa=core.PRESSURE_PA,
                y_hexane=1.0,
            )
            rates = march.rates(
                row["temperature_k"], row["beta"] * film.heat_transfer_w_m2_k / cp.hexane_j_kg_k
            )
            if abs(rates["beta"] - row["beta"]) > 1e-12 or rates["water_out_kg_s"] != 0:
                raise ValueError("native pure-hexane feed-rate reconstruction failed")
            balance.accept(
                row["t_s"],
                storage(row),
                [rates[k] for k in ("water_out_kg_s", "hexane_out_kg_s", "energy_out_w")],
                "native Faner feed",
            )

    def accept(transport, surface, kind):
        balance.accept(
            transport["time_s"],
            shared.radial_storage(transport, surface, config, area),
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
    accept(
        native_birth["after_transport"], native_birth["surface_water_node"], "native Faner birth"
    )
    for line in (folder / "front.history.jsonl").read_text().splitlines():
        row = json.loads(line)
        if row["kind"] not in ("initial", "refusal"):
            accept(row["after_transport"], row["surface"], "native B " + row["kind"])
    endpoint = json.loads((folder / "endpoint.json").read_text())
    accept(
        endpoint["after_transport"],
        endpoint["surface"],
        "native endpoint and base-face restoration",
    )
    return balance, {
        "passed": True,
        "complete_history": False,
        "scope": "original-feed prefix through C activation; not a complete history",
        "native_source_pins_checked": len(native_pins),
        "physical_rows_audited": len(balance.rows),
        "last_time_s": balance.clock,
        "maximum_independent_normalized_W_H_U": balance.maximum,
        "original_feed_W_H_U": balance.initial,
        "original_scales": balance.scales,
        "last_actual_W_H_U": balance.last_actual,
        "last_signed_defects_W_H_U": balance.last_signed_defects,
        "physically_qualifying": False,
        "convergence_established": False,
        "submission_released": False,
    }


def component_error_budgets_per_advance(balance, C_initial, maximum_advances):
    """Reserve half both existing margins; a spent original margin refuses."""
    if not isinstance(maximum_advances, int) or maximum_advances <= 0:
        raise ValueError("positive declared maximum physical advances required")
    margins = [
        balance.tolerance * balance.scales[i] - abs(balance.last_signed_defects[i]) for i in (0, 1)
    ]
    if min(margins) <= 0 or not all(math.isfinite(x) and x > 0 for x in C_initial[:2]):
        raise ValueError("no positive original-feed component margin available")
    return tuple(
        min(0.5 * margins[i], 0.5 * balance.tolerance * C_initial[i]) / maximum_advances
        for i in (0, 1)
    )
