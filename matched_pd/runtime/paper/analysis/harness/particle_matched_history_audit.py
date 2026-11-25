"""Rebuild cold-feed storage and signed exterior integrals across matched A/B/C.

Uses archived primitives and native state functions, not producer ledger deltas.
The mixed feed has uniform PD rows followed by resolved thermal main-A rows.
Native exact-face states and separate restart segments remain distinguishable.
"""

# Harness imports follow the explicit repository search path.
# ruff: noqa: E402
import argparse
import hashlib
import json
import math
from pathlib import Path
import pickle
import sys

from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "paper/analysis/harness")]
import audit_particle_qualification_case as shared
import pd_prescribed_bed_bath as model
import particle_matched_feed as pd_run
import particle_pd_diffusion_boundary as transport

core = model.core


def inverse_beta(pair, tb, yb, ti, yi, pressure, radius, beta, forcing):
    boundary = core.reduced_film_boundary(
        temperature_k=tb,
        y_hexane=yb,
        pair=pair,
        forcing_timescale_s=forcing,
        label="independent matched-history film reconstruction",
    )
    reference = core.ct.prepare_reduced_film_surface_state(
        boundary,
        cell_temperature_k=ti,
        pressure_pa=pressure,
        particle_radius_m=radius,
        outer_half_cell_distance_m=radius / 2,
        particle_thermal_conductivity_w_m_k=0.24,
        binary_gas_interaction_k_wh=boundary.binary_gas_interaction_k_wh,
    )
    film = core.ct._reduced_film_coefficients_at_surface_temperature(reference, ti)
    cp = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=ti,
        bulk_temperature_k=tb,
        pressure_pa=pressure,
        y_hexane=yi,
    )

    def evaluate(total):
        flux = core.efr.component_fluxes_for_supplied_total_mass_flux(
            film,
            total_mass_flux_kg_m2_s=total,
            interface_hexane_mass_fraction=core._hexane_mass_fraction(yi),
            bulk_hexane_mass_fraction=core._hexane_mass_fraction(yb),
        )
        ack = core.efr.ackermann_heat_transfer(
            film,
            fluxes=flux,
            water_film_heat_capacity_j_kg_k=cp.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=cp.hexane_j_kg_k,
            enforce_production_domain=False,
        )
        return flux, ack

    delta = (beta - evaluate(0.0)[1].beta) * film.heat_transfer_w_m2_k
    lo, hi = sorted((delta / cp.water_j_kg_k, delta / cp.hexane_j_kg_k))
    if delta == 0.0 or lo == hi:
        total = lo
    else:
        padding = 64 * math.ulp(max(abs(lo), abs(hi), 1.0))
        total = brentq(
            lambda value: evaluate(value)[1].beta - beta,
            lo - padding,
            hi + padding,
            xtol=1e-14,
            rtol=1e-14,
        )
    flux, ack = evaluate(total)
    if abs(ack.beta - beta) > 1e-10 or abs(ack.beta) > 3.0:
        raise ValueError("independent exact-film beta reconstruction failed")
    gas_heat_in = film.heat_transfer_w_m2_k * ack.factor * (tb - ti)
    hh = model.nc.hexane_partial_enthalpy(ti, pressure, yi).value_j_mol / core.hx.M
    hw = model.nc.water_partial_enthalpy(ti, pressure, yi).value_j_mol / core.wa.M
    return flux, flux.hexane_kg_m2_s * hh + flux.water_kg_m2_s * hw - gas_heat_in


def audit_prefix(case, histories, *, endpoint=None, dry_histories=(), stop_file=None):
    auditor_path = Path(__file__)
    auditor_entry_sha256 = hashlib.sha256(auditor_path.read_bytes()).hexdigest()

    def check_stop():
        if stop_file is not None and stop_file.exists():
            raise SystemExit("independent audit stopped cooperatively; no passed report written")

    feed = json.loads((case / "feed.json").read_text())
    if feed.get("pd_boundary_codec") != pd_run.FEED_CODEC:
        raise ValueError(
            "matched production feed codec required; diagnostic generations are distinct"
        )
    checkpoint = pickle.loads((case / "birth.accepted.pkl").read_bytes())
    record = json.loads((case / "case.json").read_text())
    pins = record["source_sha256_at_entry"]
    if not pins or any(
        hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != digest
        for path, digest in pins.items()
    ):
        raise ValueError("matched producer source differs from its recorded generation")
    birth = json.loads((case / "birth.json").read_text())
    if birth["feed_sha256"] != hashlib.sha256((case / "feed.json").read_bytes()).hexdigest():
        raise ValueError("native birth feed hash differs from archived cold-feed history")
    config = checkpoint["config"]
    bath_path = Path(checkpoint["bath_path"])
    bath = json.loads(bath_path.read_text())
    if record["bath_sha256"] != hashlib.sha256(bath_path.read_bytes()).hexdigest():
        raise ValueError("native case bath hash differs from the archived bath")
    source = feed["pd_source"]
    if source != model.source_gas(bath, model.PdReference(**source["reference"])):
        raise ValueError("PD gas/source construction differs from its archived inputs")
    pd_run.MatchedPD.geometry = feed["pd_diffusion_geometry"]
    with shared.declared_retained_water_law(record), transport.installed_pd_diffusion():
        pd = pd_run.MatchedPD(
            bath,
            source,
            forcing_timescale_s=feed["pd_film_forcing_timescale_s"],
            transport_mode="explicit_native_property_diffusion_cell",
            bed_boiling=True,
            water_film_energy="surface",
        )
        if pd_run.transfer_pair_record(pd.pd_pair) != feed["pd_transfer_pair"]:
            raise ValueError(
                "PD property/geometry pair differs from the independently rebuilt pair"
            )
        if pd.pd_contact_fraction != feed["pd_exposed_area_fraction"]:
            raise ValueError("PD top-area exposure differs")
        material = core.sphere._material_wet_params(pd.params, pd.wet.X_water, pd.wet.w_o)
        dry, area, pressure = pd.dry, pd.area, pd.pressure

        def storage(row):
            water = dry * row["external_water_kg_kg_dry"]
            if water and "interface_temperature_k" not in row:
                raise ValueError("nonzero external-water stock needs its archived temperature")
            # The cold t=0 datum has no solved surface and exactly zero film.
            # Its film energy is zero; no initial interface phase is invented.
            ti = row.get("interface_temperature_k", row["temperature_k"])
            if "shell_temperatures_k" not in row:
                state = core.sphere.initialize_qualified_feed(
                    pd.params,
                    total_hexane_loading=shared.feed_row_total_hexane(
                        row, lambda t: core.sphere.critical_loading(t, pd.params)
                    ),
                    temperature_k=row["temperature_k"],
                    retained_water_loading=row["retained_water_kg_kg_dry"],
                    oil_fraction=pd.wet.w_o,
                )
                h, energy = state.total_hexane_mass_kg, state.total_internal_energy_j
            else:
                temperatures = row["shell_temperatures_k"]
                grid = core.q.backend.uniform_grid(len(temperatures), pd.radius)
                masses = [v * material.rho_dm_p for v in grid.volumes]
                critical = [core.sphere.critical_loading(t, pd.params) for t in temperatures]
                attached = row["attached_hexane_kg"]
                h = math.fsum([*(m * x for m, x in zip(masses, critical)), attached])
                energy = math.fsum(
                    [
                        *(
                            m * core.sphere.combined_high_loading_specific_energy(t, x, material)
                            for m, t, x in zip(masses, temperatures, critical)
                        ),
                        attached * core.hx.state_Tp(temperatures[-1], pressure, "liquid").u_mass,
                    ]
                )
            total = [
                dry * row["retained_water_kg_kg_dry"] + water,
                h,
                energy + water * core.wa.state_Tp(ti, pressure, "liquid").u_mass,
            ]
            if abs(total[1] - dry * row["hexane_kg_kg_dry"]) / (dry * 0.4) > 1e-10:
                raise ValueError("feed primitive hexane storage differs from reported total")
            if abs(total[2] - row["total_particle_plus_film_energy_j"]) / (dry * 2e5) > 1e-10:
                raise ValueError("feed primitive energy storage differs from reported total")
            return total

        rows = feed["rows"]
        if rows[0]["t_s"] != 0.0:
            raise ValueError("a cold-feed datum at t=0 is required")
        initial = storage(rows[0])
        balance = shared.OriginalFeedBalance(
            initial,
            [initial[0], initial[1], checkpoint["carry"]["import_capacity_scale_j"]],
            1e-10,
        )
        counts = dict(PD=0, main_A=0, birth=0, B=0, endpoint=0, C=0)
        pd_uptake, pd_gas_release, pd_boiling_release = [], [], []

        def accept(clock, actual, rates, kind, **options):
            try:
                balance.accept(clock, actual, rates, kind, **options)
            except Exception as error:
                error.add_note(f"independent matched-history audit: {kind}, t={clock!r} s")
                raise

        for previous, row in zip(rows, rows[1:]):
            check_stop()
            in_pd = previous["t_s"] < pd.pd_end
            ti, yi = row["interface_temperature_k"], row["interface_y_hexane"]
            gas = core.bg.state(ti, pressure, 1 - yi, yi)
            ah = gas.fugacity_hexane_pa / core.hx.state_Tp(ti, pressure, "liquid").fugacity
            aw = gas.fugacity_water_pa / core.wa.state_Tp(ti, pressure, "liquid").fugacity
            if abs(math.log(ah)) > 1e-10:
                raise ValueError("feed surface H phase equality failed")
            if row["external_water_kg_kg_dry"] > 0.0:
                if abs(math.log(aw)) > 1e-10:
                    raise ValueError("present external-water phase equality failed")
            elif aw > 1.0 + 1e-10:
                raise ValueError("absent external-water phase is supersaturated")
            if in_pd:
                tb, yb = source["temperature_k"], source["y_hexane"]
                pair, exposure, forcing = pd.pd_pair, pd.pd_contact_fraction, pd.forcing_timescale_s
            else:
                tb, yb = core.interpolate_bath(bath["waypoints"], row["t_s"])
                pair, exposure, forcing = pd.pair, 1.0, 0.075
            flux, energy = inverse_beta(
                pair, tb, yb, ti, yi, pressure, pd.radius, row["beta"], forcing
            )
            energy_rate = area * exposure * energy - (
                pd.dry * pd.reference.mean_duty_w_kg_dry if in_pd else 0.0
            )
            boiled = row.get("bed_boiling_hexane_kg", 0.0) if in_pd else 0.0
            if in_pd:
                signed_gas_amount = (
                    (row["t_s"] - previous["t_s"]) * area * exposure * flux.hexane_kg_m2_s
                )
                pd_uptake.append(max(-signed_gas_amount, 0.0))
                pd_gas_release.append(max(signed_gas_amount, 0.0))
                pd_boiling_release.append(boiled)
            hg = (
                model.nc.hexane_partial_enthalpy(
                    pd.pure_hexane_temperature_k, pressure, 1.0
                ).value_j_mol
                / core.hx.M
            )
            accept(
                row["t_s"],
                storage(row),
                [
                    area * exposure * flux.water_kg_m2_s,
                    area * exposure * flux.hexane_kg_m2_s,
                    energy_rate,
                ],
                "independent PD" if in_pd else "independent resolved main A",
                event_terms=[[], [boiled], [boiled * hg]],
            )
            counts["PD" if in_pd else "main_A"] += 1

        def radial(state, surface, kind, *, events=()):
            accept(
                state["time_s"],
                shared.radial_storage(state, surface, config, area, pressure),
                [
                    area * surface[key]
                    for key in (
                        "external_water_flux_kg_m2_s",
                        "external_hexane_flux_kg_m2_s",
                        "external_energy_flux_w_m2",
                    )
                ],
                kind,
                **(
                    {
                        "event_terms": [
                            [event["stream_W_H_U"][k] for event in events] for k in range(3)
                        ]
                    }
                    if events
                    else {}
                ),
            )

        radial(
            birth["audit"]["after_transport"], birth["audit"]["surface_water_node"], "native birth"
        )
        counts["birth"] += 1
        for path in histories:
            with path.open() as stream:
                for line in stream:
                    check_stop()
                    row = json.loads(line)
                    if row["kind"] in ("initial", "refusal"):
                        continue
                    radial(row["after_transport"], row["surface"], "native B " + row["kind"])
                    counts["B"] += 1
        if endpoint is not None:
            radial(endpoint["after_transport"], endpoint["surface"], "native core endpoint")
            counts["endpoint"] += 1
        for path in dry_histories:
            with path.open() as stream:
                for line in stream:
                    check_stop()
                    row = json.loads(line)
                    if row["kind"] in ("initial", "refusal"):
                        continue
                    events = row.get("kernel_envelope_hexane_exhausted_regime", {}).get(
                        "events", ()
                    )
                    radial(
                        row["after_transport"],
                        row["surface"],
                        "native C " + row["kind"],
                        events=events,
                    )
                    counts["C"] += 1
        if auditor_entry_sha256 != hashlib.sha256(auditor_path.read_bytes()).hexdigest():
            raise ValueError("independent auditor source changed during execution")
        report = dict(
            passed=True,
            clock_s=balance.clock,
            counts=counts,
            maximum_normalized_W_H_U=balance.maximum,
            last_signed_defects_W_H_U=balance.last_signed_defects,
            auditor_sha256=auditor_entry_sha256,
            complete_history=balance.clock == shared.bath_outlet_time_s(bath),
            PD_signed_exchange=dict(
                gross_inward_gas_hexane_g_kg_dry=1000 * math.fsum(pd_uptake) / dry,
                gross_outward_gas_hexane_g_kg_dry=1000 * math.fsum(pd_gas_release) / dry,
                external_liquid_boiling_hexane_g_kg_dry=1000 * math.fsum(pd_boiling_release) / dry,
                definition="positive/negative parts of native signed film integrals, plus the separately booked boiling export; no state or flux is clamped",
            ),
            scope="independent primitive storage and cold-feed exterior/source integrals; no mesh/time qualification inferred",
        )
        return balance, report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--history", type=Path, action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--stop-file", type=Path)
    parser.add_argument("--endpoint", type=Path)
    parser.add_argument("--dry-history", type=Path, action="append", default=[])
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError("refuse overwrite of independent audit evidence")
    _, report = audit_prefix(
        args.case,
        args.history,
        stop_file=args.stop_file,
        endpoint=json.loads(args.endpoint.read_text()) if args.endpoint else None,
        dry_histories=args.dry_history,
    )
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
