"""Fresh 64-cell DT history with the declared PD bed/boiling prescription.

Reuse the original case producer, native birth/front/dry operators and independent
OriginalFeedBalance. Only the explicitly stamped feed recipe and its independent
exterior-rate reconstruction differ. No native solver or acceptance gate changes.
"""

import argparse
from dataclasses import asdict
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import pickle
import sys
import time
import types

from scipy.optimize import brentq

import pd_prescribed_bed_bath as model
import run_particle_qualification_case as producer
import audit_particle_qualification_case as independent
import particle_original_feed_balance as prefix

ROOT = Path(__file__).resolve().parents[3]
CODEC = "paper1-pd-equilibrium-source-k2-native-boiling-v1"
CONSTANT_HEAT_CODEC = "paper1-pd-equilibrium-source-constant-heat-native-boiling-v1"
SURFACE_FILM_CODEC = "paper1-pd-constant-heat-surface-water-energy-native-boiling-v2"
core = model.core


class CooperativeStop(BaseException):
    pass


def inverse_film_beta(march, row):
    """Existing independent nonlinear Stefan/Ackermann inverse, no linear partition."""
    if march.time < march.pd_end:
        boundary = march.pd_boundary
        ti, yi = row["interface_temperature_k"], row["interface_y_hexane"]
    else:
        tb, yb = core.interpolate_bath(march.bath["waypoints"], row["t_s"])
        boundary = core.reduced_film_boundary(
            temperature_k=tb,
            y_hexane=yb,
            pair=march.pair,
            forcing_timescale_s=0.075,
            label="original main-tray boundary audit",
        )
        ti, yi = march.two_liquid_trace
    ref = core.ct.prepare_reduced_film_surface_state(
        boundary,
        cell_temperature_k=row["temperature_k"],
        pressure_pa=march.pressure,
        particle_radius_m=march.radius,
        outer_half_cell_distance_m=march.radius / 2,
        particle_thermal_conductivity_w_m_k=core.CONDUCTIVITY,
        binary_gas_interaction_k_wh=boundary.binary_gas_interaction_k_wh,
    )
    film = core.ct._reduced_film_coefficients_at_surface_temperature(ref, ti)
    cp = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=ti,
        bulk_temperature_k=boundary.temperature_k,
        pressure_pa=march.pressure,
        y_hexane=yi,
    )
    wi, wb = core._hexane_mass_fraction(yi), core._hexane_mass_fraction(boundary.y_hexane)

    def beta(total):
        parts = core.efr.component_fluxes_for_supplied_total_mass_flux(
            film,
            total_mass_flux_kg_m2_s=total,
            interface_hexane_mass_fraction=wi,
            bulk_hexane_mass_fraction=wb,
        )
        return core.efr.ackermann_heat_transfer(
            film,
            fluxes=parts,
            water_film_heat_capacity_j_kg_k=cp.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=cp.hexane_j_kg_k,
            enforce_production_domain=False,
        ).beta

    increment = (row["beta"] - beta(0.0)) * film.heat_transfer_w_m2_k
    lo, hi = sorted(
        (
            increment / min(cp.water_j_kg_k, cp.hexane_j_kg_k),
            increment / max(cp.water_j_kg_k, cp.hexane_j_kg_k),
        )
    )
    rounding = 64 * math.ulp(max(abs(lo), abs(hi), 1.0))
    total = (
        0.0
        if increment == 0
        else brentq(
            lambda m: beta(m) - row["beta"],
            lo - rounding,
            hi + rounding,
            xtol=1e-14,
            rtol=1e-14,
            maxiter=80,
        )
    )
    if abs(beta(total) - row["beta"]) > 1e-12:
        raise ValueError("independent nonlinear beta reconstruction failed")
    return total, ti, yi


class CaseAdapter:
    def __init__(self, out, stop_file, *, pd_heating="circulating", water_film_energy="body"):
        if pd_heating not in ("constant", "circulating"):
            raise ValueError("declare constant or circulating PD heating")
        self.out = out
        self.stop_file = stop_file
        self.started = time.monotonic()
        self.last_check = 0.0
        self.source = None
        self.members = None
        self.pd_heating = pd_heating
        if water_film_energy not in ("body", "surface"):
            raise ValueError("declare body or surface external-water energy")
        if water_film_energy == "surface" and pd_heating != "constant":
            raise ValueError("surface-energy repair is declared for the constant-heating case")
        self.water_film_energy = water_film_energy
        self.codec = CONSTANT_HEAT_CODEC if pd_heating == "constant" else CODEC
        if water_film_energy == "surface":
            self.codec = SURFACE_FILM_CODEC

    def write(self, name, value):
        path = self.out / name
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(
            json.dumps(value, indent=2, default=producer.front.json_default) + "\n",
            encoding="utf-8",
        )
        temporary.replace(path)

    def check_stop(self, *unused):
        now = time.monotonic()
        if now - self.last_check < 0.1:
            return
        self.last_check = now
        if self.stop_file.exists() or now - self.started > 8 * 3600:
            raise CooperativeStop("cooperative resource or eight-hour wall stop")

    def run_feed(self, args):
        if args.radius_m != 0.0015 or args.cells < 2 or not args.paper_properties:
            raise ValueError("this prescription requires Paper-1 R=1.5mm and radial cells >=2")
        bath = json.loads(args.bath.read_text())
        reference = model.PdReference()
        if not all(w.get("pressure_pa") == reference.pressure_pa for w in bath["waypoints"]):
            raise ValueError("this recipe requires the declared constant pressure schedule")
        source = model.source_gas(bath, reference)
        self.source = source
        renewal = model.BedRenewal(2.0) if self.pd_heating == "circulating" else None
        forcing_time = renewal.period_s if renewal else reference.residence_s / reference.tray_count
        members = [
            model.PrescribedPdFeed(
                bath,
                source,
                forcing_timescale_s=forcing_time,
                transport_mode="canonical_reference",
                bed_boiling=True,
                water_film_energy=self.water_film_energy,
            )
            for _ in range(2 if renewal else 1)
        ]
        self.members = members
        errors = []
        activation = None
        last_print = -1
        pd_horizon = min(reference.residence_s, args.max_time_s)
        while members[0].time < pd_horizon:
            self.check_stop()
            clock = members[0].time
            tray_time = reference.residence_s / 3
            dt = min(
                args.dt_s,
                pd_horizon - clock,
                (math.floor(clock / tray_time) + 1) * tray_time - clock,
            )
            if renewal:
                half = renewal.period_s / 2
                dt = min(dt, (math.floor(clock / half) + 1) * half - clock)
                heat = renewal.heat_inputs(
                    clock, [m.state.temperature_k for m in members], reference
                )
                inputs = heat["net_w_kg_dry"]
            else:
                inputs = [reference.mean_duty_w_kg_dry]
            for m, q in zip(members, inputs):
                m.net_heat_w_kg_dry = q
            trials = [m.trial(dt) for m in members]
            actual_dt = min(s["dt"] for s in trials)
            trials = [
                s if s["dt"] == actual_dt else m.trial(actual_dt) for m, s in zip(members, trials)
            ]
            for m, s in zip(members, trials):
                m.commit_trial(s)
            if int(members[0].time // 25) > last_print:
                last_print = int(members[0].time // 25)
                progress = dict(
                    stage="fresh_PD",
                    time_s=members[0].time,
                    wall_s=time.monotonic() - self.started,
                    h_kg_kg_dry=[m.state.total_hexane_loading for m in members],
                )
                self.write("progress.json", progress)
                print(json.dumps(progress), flush=True)
        march = members[0]  # One actual representative, never an averaged state.
        if march.time >= reference.residence_s:
            pd_rows = march.rows
            initial_h = pd_rows[0]["hexane_kg_kg_dry"]
            peak = max(pd_rows, key=lambda row: row["hexane_kg_kg_dry"])
            boiling = [row for row in pd_rows if row.get("bed_boiling_active", False)]
            self.write(
                "pd_summary.json",
                dict(
                    status="pd_complete",
                    pd_heating=self.pd_heating,
                    time_s=march.time,
                    initial_hexane_kg_kg_dry=initial_h,
                    peak_hexane_kg_kg_dry=peak["hexane_kg_kg_dry"],
                    peak_age_s=peak["t_s"],
                    exit_hexane_kg_kg_dry=march.state.total_hexane_loading,
                    net_removed_fraction=1 - march.state.total_hexane_loading / initial_h,
                    first_bed_boiling_s=boiling[0]["t_s"] if boiling else None,
                    exit_temperature_c=march.state.temperature_k - 273.15,
                    specific_indirect_heat_w_kg_dry=reference.mean_duty_w_kg_dry,
                    radial_front_activated_in_pd=False,
                    scope="member 0 PD prefix; complete DT history and independent audit still pending",
                ),
            )
            self.write(
                "progress.json", dict(stage="PD_complete_entering_main_trays", time_s=march.time)
            )
        while march.time < args.max_time_s:
            self.check_stop()
            dt = min(args.steam_dt_s, args.max_time_s - march.time)
            try:
                march.solve(dt)
            except Exception as exc:
                errors.append(dict(time_s=march.time, error=str(exc)))
                try:
                    march.solve(dt, event=True)
                    radial, ledger = core.sphere.activate_radially(march.state, march.params)
                    activation = dict(
                        radial=asdict(radial),
                        ledger=asdict(ledger),
                        external_water_mass_kg_carried_separately=march.free_water,
                        external_water_energy_j_carried_separately=march.total_energy
                        - march.state.total_internal_energy_j,
                    )
                except Exception as event_error:
                    errors[-1]["event_error"] = str(event_error)
                break
        return dict(
            status=(
                "ACTUAL_FEED_TO_RADIAL_ACTIVATION"
                if activation
                else ("ACTUAL_FEED_MARCH_REFUSED" if errors else "ACTUAL_FEED_PARTIAL_HORIZON")
            ),
            complete_dt_destiny=False,
            physically_qualifying=False,
            bath_was_modified=True,
            boundary_history="explicit PD bed/boiling source, then original main-tray gas waypoints",
            pd_boundary_codec=self.codec,
            pd_heat_prescription=self.pd_heating,
            pd_film_forcing_timescale_s=forcing_time,
            pd_source=source,
            pd_renewal=asdict(renewal) if renewal else None,
            pd_transport_mode="canonical_reference",
            pd_bed_boiling=True,
            selected_member=0,
            companion_pd_rows=members[1].rows if renewal else [],
            variant_source_sha256={
                str(Path(p).relative_to(ROOT)): hashlib.sha256(Path(p).read_bytes()).hexdigest()
                for p in (__file__, model.__file__, model.feed.__file__)
            },
            dt_s=args.dt_s,
            radius_m=args.radius_m,
            high_loading_thermal_cells=1,
            property_parameters=asdict(model.ParticleModelParameters()),
            surface_to_bulk_conductance_w_m2_k=march.internal_conductance_w_m2_k,
            external_water_caloric=f"native liquid internal energy at {self.water_film_energy} temperature",
            external_water_energy_temperature=self.water_film_energy,
            external_water_temperature_k=(
                march.ti if self.water_film_energy == "surface" else march.state.temperature_k
            ),
            interface_temperature_k=march.ti,
            interface_y_hexane=march.yi,
            kernel_envelope_feed_pressure_pa=reference.pressure_pa,
            rows=march.rows,
            activation=activation,
            errors=errors,
        )

    def audit(self, folder, bath_path, live_source_hashes=None, until_core_endpoint=False):
        record = json.loads((folder / "case.json").read_text())
        independent.validate_producer_sources(record, live_source_hashes)
        native_pins = {
            p: h for p, h in record["source_sha256_at_entry"].items() if p.startswith("src/")
        }
        if not native_pins or any(
            hashlib.sha256((ROOT / p).read_bytes()).hexdigest() != h for p, h in native_pins.items()
        ):
            raise ValueError("native source identity differs from fresh producer")
        if not record["fresh_original_feed"] or not record["single_runtime_source_generation"]:
            raise ValueError("actual fresh original-feed single-runtime case required")
        feed = json.loads((folder / "feed.json").read_text())
        if feed.get("pd_boundary_codec") != self.codec or feed.get("selected_member") != 0:
            raise ValueError("PD feed codec/representative differs")
        for p, digest in feed["variant_source_sha256"].items():
            if hashlib.sha256((ROOT / p).read_bytes()).hexdigest() != digest:
                raise ValueError("PD variant source identity changed")
        bath = json.loads(bath_path.read_text())
        if record["bath_sha256"] != hashlib.sha256(bath_path.read_bytes()).hexdigest():
            raise ValueError("bath hash differs")
        source = model.source_gas(bath)
        if source != feed["pd_source"]:
            raise ValueError("offline PD source differs from exact declared source")
        birth = json.loads((folder / "birth.json").read_text())
        if birth["feed_sha256"] != hashlib.sha256((folder / "feed.json").read_bytes()).hexdigest():
            raise ValueError("birth parent feed differs")
        checkpoint = pickle.loads((folder / "birth.accepted.pkl").read_bytes())
        config = checkpoint["config"]
        march = model.PrescribedPdFeed(
            bath,
            source,
            forcing_timescale_s=feed["pd_film_forcing_timescale_s"],
            transport_mode="canonical_reference",
            bed_boiling=True,
            water_film_energy=feed.get("external_water_energy_temperature", "body"),
        )
        if (
            config.dry.pressure_pa != march.pressure
            or config.wet.wet.oil_heat_capacity_law is not march.wet.oil_heat_capacity_law
            or config.dry.pore.oil_heat_capacity_law is not march.wet.oil_heat_capacity_law
            or config.dry.binary_diffusivity.temperature_law
            is not model.ParticleModelParameters().binary_law
            or config.dry.binary_diffusivity.value_m2_s
            != model.ParticleModelParameters().binary_reference_m2_s
        ):
            raise ValueError("native pressure/property identities differ")
        mass, area = march.dry, march.area

        def body(row):
            return core.sphere.initialize_qualified_feed(
                march.params,
                total_hexane_loading=independent.feed_row_total_hexane(
                    row, lambda t: core.sphere.critical_loading(t, march.params)
                ),
                temperature_k=row["temperature_k"],
                retained_water_loading=row["retained_water_kg_kg_dry"],
                oil_fraction=march.wet.w_o,
            )

        def storage(row):
            state = body(row)
            water = mass * row["external_water_kg_kg_dry"]
            return [
                state.retained_water_mass_kg + water,
                state.total_hexane_mass_kg,
                state.total_internal_energy_j
                + water
                * core.wa.state_Tp(
                    (
                        row["external_water_temperature_k"]
                        if self.water_film_energy == "surface"
                        else row["temperature_k"]
                    ),
                    march.pressure,
                    "liquid",
                ).u_mass,
            ]

        initial = storage(feed["rows"][0])
        balance = independent.OriginalFeedBalance(
            initial,
            [initial[0], initial[1], checkpoint["carry"]["import_capacity_scale_j"]],
            config.dry.ledger_tolerance,
        )
        companion = {r["t_s"]: r for r in feed["companion_pd_rows"]}
        renewal = model.BedRenewal(**feed["pd_renewal"]) if feed["pd_renewal"] else None
        if feed["pd_heat_prescription"] != self.pd_heating:
            raise ValueError("PD heating prescription differs")
        for old, row in zip(feed["rows"], feed["rows"][1:]):
            self.check_stop()
            march.time = old["t_s"]
            dt = row["t_s"] - old["t_s"]
            in_pd = march.time < march.pd_end
            if in_pd:
                if renewal:
                    heat = renewal.heat_inputs(
                        march.time,
                        [old["temperature_k"], companion[march.time]["temperature_k"]],
                        march.reference,
                    )
                    march.net_heat_w_kg_dry = heat["net_w_kg_dry"][0]
                else:
                    march.net_heat_w_kg_dry = march.reference.mean_duty_w_kg_dry
                if row["net_indirect_and_bed_heat_w_kg_dry"] != march.net_heat_w_kg_dry:
                    raise ValueError("PD source duty differs from independent bed reconstruction")
            total, ti, yi = inverse_film_beta(march, row)
            if in_pd:
                if row["temperature_k"] > march.pure_hexane_temperature_k:
                    raise ValueError("PD body-liquid superheat in the phase-limited recipe")
                phase = core.bg.state(ti, march.pressure, 1 - yi, yi)
                h_ratio = (
                    phase.fugacity_hexane_pa
                    / core.hx.state_Tp(ti, march.pressure, "liquid").fugacity
                )
                w_ratio = (
                    phase.fugacity_water_pa
                    / core.wa.state_Tp(ti, march.pressure, "liquid").fugacity
                )
                if abs(math.log(h_ratio)) > 1e-10:
                    raise ValueError("PD surface H phase equality failed")
                if (row["external_water_kg_kg_dry"] > 0 and abs(math.log(w_ratio)) > 1e-10) or (
                    row["external_water_kg_kg_dry"] == 0 and w_ratio > 1
                ):
                    raise ValueError("PD water phase selection failed")
            march.ti, march.yi = ti, yi
            native, ack, _, energy, _, _, active = march.interface(
                row["temperature_k"], total, row["t_s"]
            )
            boiled = row.get("bed_boiling_hexane_kg", 0.0) if in_pd else 0.0
            if boiled < 0:
                raise ValueError("negative generated boiling stream")
            if row.get("bed_boiling_active", False):
                if (
                    row["temperature_k"] != march.pure_hexane_temperature_k
                    or row["external_water_kg_kg_dry"] != 0
                ):
                    raise ValueError("bed boiling phase equality/absent-water condition differs")
                old_u = storage(old)[2]
                film_h = dt * active * native.hexane_molar_flux_mol_m2_s * core.hx.M
                reference_state = core.sphere.initialize_qualified_feed(
                    march.params,
                    total_hexane_loading=old["hexane_kg_kg_dry"] - film_h / mass,
                    temperature_k=march.pure_hexane_temperature_k,
                    retained_water_loading=old["retained_water_kg_kg_dry"],
                    oil_fraction=march.wet.w_o,
                )
                hg = (
                    model.nc.hexane_partial_enthalpy(
                        march.pure_hexane_temperature_k, march.pressure, 1.0
                    ).value_j_mol
                    / core.hx.M
                )
                ul = core.hx.state_Tp(
                    march.pure_hexane_temperature_k, march.pressure, "liquid"
                ).u_mass
                expected = (old_u - reference_state.total_internal_energy_j - dt * energy) / (
                    hg - ul
                )
                if abs(expected - boiled) / (mass * feed["rows"][0]["hexane_kg_kg_dry"]) > 1e-10:
                    raise ValueError(
                        "boiling stream differs from independently reconstructed native heat balance"
                    )
                energy += boiled / dt * hg
            elif boiled != 0:
                raise ValueError("inactive boiling branch carries vapour")
            balance.accept(
                row["t_s"],
                storage(row),
                [
                    active * native.water_molar_flux_mol_m2_s * core.wa.M,
                    active * native.hexane_molar_flux_mol_m2_s * core.hx.M + boiled / dt,
                    energy,
                ],
                "independent PD/main phase, beta, bed heat and separate boiling rates",
            )

        def accept(transport, surface, kind):
            balance.accept(
                transport["time_s"],
                independent.radial_storage(
                    transport, surface, config, area, pressure_pa=march.pressure
                ),
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

        accept(birth["audit"]["after_transport"], birth["audit"]["surface_water_node"], "birth")
        for name in (
            "front.history.jsonl",
            *([] if until_core_endpoint else ["dry.history.jsonl"]),
        ):
            path = folder / name
            if path.exists():
                with path.open() as stream:
                    for line in stream:
                        self.check_stop()
                        row = json.loads(line)
                        if row["kind"] not in ("initial", "refusal"):
                            accept(row["after_transport"], row["surface"], name + " " + row["kind"])
                if name == "front.history.jsonl" and (folder / "endpoint.json").exists():
                    row = json.loads((folder / "endpoint.json").read_text())
                    accept(row["after_transport"], row["surface"], "core endpoint")
        report = dict(
            passed=True,
            maximum_independent_normalized_W_H_U=balance.maximum,
            last_time_s=balance.clock,
            physical_rows_audited=len(balance.rows),
            complete_history=record["complete_history"] and not until_core_endpoint,
            source_codec=self.codec,
            physically_qualifying=False,
        )
        return balance, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--stop-file", type=Path, required=True)
    parser.add_argument("--probe-seconds", type=float)
    parser.add_argument("--pd-heating", choices=("constant", "circulating"), default="constant")
    parser.add_argument("--cells", type=int, choices=(32, 64, 128), default=64)
    parser.add_argument("--time-level", type=int, choices=(0, 1, 2), default=1)
    parser.add_argument("--inlet-temperature-c", type=float, default=49.0)
    parser.add_argument("--water-film-energy", choices=("body", "surface"), default="body")
    args = parser.parse_args()
    args.out = args.out.resolve()
    args.stop_file = args.stop_file.resolve()
    if args.out.exists():
        raise SystemExit("refuse overwrite")
    args.out.mkdir(parents=True)
    adapter = CaseAdapter(
        args.out,
        args.stop_file,
        pd_heating=args.pd_heating,
        water_film_energy=args.water_film_energy,
    )
    original = (
        ROOT / "paper/analysis/results_2026-10-01/dtdc_particle_destiny/destiny_v3_exit105.json"
    )
    bath = json.loads(original.read_text())
    if not 49.0 <= args.inlet_temperature_c <= 58.0:
        raise ValueError("inlet sensitivity restricted to the declared 49–58 C range")
    bath["particle"]["feed_temperature_K"] = args.inlet_temperature_c + 273.15
    bath["particle"]["radius_m"] = 0.0015
    bath["radii_m"] = [0.0015]
    bath["radii_origin"] = (
        "native Paper 1 DT reference radius; explicit parent-bath metadata override"
    )
    for diagnostic in ("Xc_at_feed_kg_kg_dry", "attached_hexane_at_feed_kg_kg_dry"):
        bath["particle"].pop(diagnostic, None)
    bath["initial_geometry_convention"] = (
        "computed natively at actual pressure, radius and feed state"
    )
    for waypoint in bath["waypoints"]:
        waypoint["pressure_pa"] = model.PdReference().pressure_pa
    bath["pd_boundary_codec"] = adapter.codec
    bath_path = args.out / "bath.json"
    bath_path.write_text(json.dumps(bath, indent=2) + "\n", encoding="utf-8")
    if args.probe_seconds is not None:
        settings = types.SimpleNamespace(
            bath=bath_path,
            radius_m=0.0015,
            cells=args.cells,
            paper_properties=True,
            dt_s=0.1 / 2**args.time_level,
            steam_dt_s=0.02 / 2**args.time_level,
            max_time_s=args.probe_seconds,
        )
        adapter.write("feed.json", adapter.run_feed(settings))
        return
    unexpected = {
        k: v
        for k, v in os.environ.items()
        if k.startswith("DTDC_")
        and k not in ("DTDC_RESOURCE_AFFINITY_MASK", "DTDC_KERNEL_ENVELOPE_PRESSURE_SCHEDULE")
    }
    if unexpected:
        raise ValueError(f"undeclared DTDC environment for this case: {sorted(unexpected)}")
    os.environ["DTDC_KERNEL_ENVELOPE_PRESSURE_SCHEDULE"] = "1"
    command = [
        "run_particle_qualification_case.py",
        "--study",
        "DT",
        "--cells",
        str(args.cells),
        "--time-level",
        str(args.time_level),
        "--bath",
        str(bath_path),
        "--out",
        str(args.out / "case"),
        "--per-area-qualification",
        "--dry-local-tolerance",
        "1e-10",
        "--conservative-components",
        "--event-first-film-depletion",
        "--eliminate-event-absent-components",
        "--curvature-corrected-depletion-prediction",
        "--memoize-proof-points",
        "--dry-time-schedule",
        "uniform",
    ]
    adapter.write(
        "run_control.json",
        dict(
            codec=adapter.codec,
            pd_heating=args.pd_heating,
            inlet_temperature_c=args.inlet_temperature_c,
            water_film_energy=args.water_film_energy,
            command=command,
            stop_file=str(args.stop_file),
            source_hashes={
                str(Path(p).relative_to(ROOT)): hashlib.sha256(Path(p).read_bytes()).hexdigest()
                for p in (__file__, model.__file__, model.feed.__file__)
            },
            original_main_bath_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
            one_worker=True,
            maximum_wall_hours=8,
            physically_qualifying=False,
        ),
    )
    old_run, old_prefix = producer.feed_driver.run, prefix.prepare_prefix
    old_birth = producer.engine.driver.advance_birth

    def observed_birth(*values, **keywords):
        try:
            result = old_birth(*values, **keywords)
            adapter.write(
                "birth_diagnostic.json", dict(status="accepted", unchanged_native_call=True)
            )
            return result
        except Exception as error:
            debug = getattr(error, "study_debug", {}) or {}
            solver = debug.get("nonlinear_solver", {}) or {}
            adapter.write(
                "birth_diagnostic.json",
                dict(
                    status="refused",
                    error=str(error),
                    unchanged_native_call=True,
                    exact_rollback=getattr(error, "rollback_state", None) is values[0],
                    maximum_scaled_residual=getattr(error, "maximum_scaled_residual", None),
                    total_residual_assembly_attempts=getattr(
                        error, "total_residual_assembly_attempts", None
                    ),
                    external_beta_coordinate=debug.get("external_beta_coordinate"),
                    nonlinear_solver_success=solver.get("success"),
                    nonlinear_solver_message=solver.get("message"),
                    birth_condition_certificate=debug.get("birth_condition_certificate"),
                ),
            )
            raise

    tool_id = next(i for i in (3, 4) if sys.monitoring.get_tool(i) is None)
    sys.monitoring.use_tool_id(tool_id, "pd-prescribed-particle-resource-stop")
    sys.monitoring.register_callback(tool_id, sys.monitoring.events.PY_START, adapter.check_stop)
    for module in list(sys.modules.values()):
        if isinstance(module, types.ModuleType):
            for value in list(vars(module).values()):
                if inspect.isfunction(value) and "advance" in value.__name__:
                    filename = Path(value.__code__.co_filename)
                    if filename.is_relative_to(ROOT / "src") or filename.is_relative_to(
                        ROOT / "paper/analysis/harness"
                    ):
                        sys.monitoring.set_local_events(
                            tool_id, value.__code__, sys.monitoring.events.PY_START
                        )
    sys.monitoring.set_local_events(
        tool_id, producer.write_json.__code__, sys.monitoring.events.PY_START
    )
    sys.monitoring.set_local_events(
        tool_id, producer.HistoryTimeControl.dry_request.__code__, sys.monitoring.events.PY_START
    )

    def prepare(folder, study, bath_path, live_source_hashes):
        if study != "DT":
            raise ValueError("this explicit PD recipe is DT only")
        record = json.loads((folder / "case.json").read_text())
        with independent.declared_retained_water_law(record):
            return adapter.audit(folder, bath_path, live_source_hashes, until_core_endpoint=True)

    try:
        producer.feed_driver.run = adapter.run_feed
        prefix.prepare_prefix = prepare
        producer.engine.driver.advance_birth = observed_birth
        sys.argv = command
        producer.main()
        record = json.loads((args.out / "case/case.json").read_text())
        if not record["complete_history"]:
            adapter.write(
                "result.json",
                dict(status="incomplete", case_status=record["status"], error=record.get("error")),
            )
            return
        with independent.declared_retained_water_law(record):
            _, report = adapter.audit(args.out / "case", bath_path)
        adapter.write("independent_audit.json", report)
        last = record["last_accepted"]
        w, h = last["particle_water_kg_kg_dry"], last["hexane_kg_kg_dry"]
        denominator = 1 + w + h
        adapter.write(
            "result.json",
            dict(
                status="complete",
                independent_audit_passed=True,
                outlet=dict(
                    hexane_ppm_wet_basis=1e6 * h / denominator,
                    moisture_pct_wet_basis=100 * w / denominator,
                    temperature_c=last["temperature_k"] - 273.15,
                ),
                physically_qualifying=False,
            ),
        )
    except CooperativeStop as exc:
        adapter.write("result.json", dict(status="stopped_incomplete", reason=str(exc)))
    except Exception as exc:
        adapter.write("result.json", dict(status="failed_incomplete", reason=repr(exc)))
        raise
    finally:
        if adapter.members is not None and not (args.out / "case/feed.json").exists():
            adapter.write(
                "feed_prefix_snapshot.json",
                dict(
                    codec=adapter.codec,
                    scope="feed prefix only; later native checkpoints remain under case/",
                    members=[
                        dict(
                            state=asdict(m.state),
                            params=asdict(m.params),
                            time_s=m.time,
                            water_kg=m.free_water,
                            energy_j=m.total_energy,
                            cumulative=m.cumulative.tolist(),
                            rows=m.rows,
                        )
                        for m in adapter.members
                    ],
                ),
            )
        producer.feed_driver.run = old_run
        prefix.prepare_prefix = old_prefix
        producer.engine.driver.advance_birth = old_birth
        sys.monitoring.free_tool_id(tool_id)


if __name__ == "__main__":
    main()
