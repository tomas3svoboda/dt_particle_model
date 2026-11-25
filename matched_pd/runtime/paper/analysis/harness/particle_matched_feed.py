"""Prescribed contact-heated PD and conservative thermal matching into B.

PD is a uniform thermal reduction. Main A uses the existing radial Fourier
operator until exact external-hexane depletion, preserving its material labels.
Surface water storage and its interface caloric donor use the same temperature.
"""

from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import particle_feed_shells as shells
import particle_pd_diffusion_boundary as transport
import pd_prescribed_bed_bath as model
from run_pd_prescribed_particle_case import CaseAdapter

FEED_CODEC = "paper1-pd-actual-properties-matched-film-v4"
ROOT = Path(__file__).resolve().parents[3]


def transfer_pair_record(pair):
    """Encode the unbounded geometry explicitly in strict JSON, not Infinity."""
    record = asdict(pair)
    if math.isinf(pair.outer_radius_m):
        if pair.geometry_assumption != "unconfined_sphere":
            raise ValueError("only the unconfined geometry has an infinite outer radius")
        record.update(outer_radius_m=None, outer_radius_limit="unbounded")
    else:
        record["outer_radius_limit"] = "finite_cell"
    return record


class MatchedPD(model.PrescribedPdFeed):
    geometry = "unconfined_sphere"

    def __init__(self, *args, **kwargs):
        original = model.pd_film_pair
        model.pd_film_pair = lambda source, radius_m=0.0015, *, transport_mode: transport.pair(
            source, radius_m, geometry=self.geometry
        )
        try:
            super().__init__(*args, **kwargs)
        finally:
            model.pd_film_pair = original
        solid_fraction = 1 - self.pd_pair.state.bed_void_fraction
        self.pd_contact_fraction = self.radius / (3 * solid_fraction * self.reference.bed_depth_m)

    def interface(self, temperature, total_mass_flux, clock):
        native, ack, residual, out, tb, yb, area = super().interface(
            temperature, total_mass_flux, clock
        )
        surface_u = (
            model.core.wa.state_Tp(self.ti, self.pressure, "liquid").u_mass * model.core.wa.M
        )
        residual += (
            area
            * native.water_molar_flux_mol_m2_s
            * (surface_u - native.native_water_liquid_internal_energy_j_mol)
        )
        return (
            replace(
                native,
                native_water_liquid_internal_energy_j_mol=surface_u,
                interface_energy_residual_w=residual,
            ),
            ack,
            residual,
            out,
            tb,
            yb,
            area,
        )


class MatchedRadialFeed(shells.MultiShellFeedMarch):
    def row(self, step, ledger):
        row = super().row(step, ledger)
        row.update(
            interface_temperature_k=self.ti,
            interface_y_hexane=self.yi,
            external_water_temperature_k=self.ti,
            external_water_energy_temperature="surface",
            surface_topology="two_liquid" if self.free_water > 0 else "hexane_only",
        )
        if step is not None:
            row["separate_stock_elimination"] = step.get("eliminated_stock_diagnostic")
        return row

    def solve(self, requested_dt, event=False):
        from particle_a_component_elimination import solve

        if self.free_water > 0:
            return solve(self, requested_dt, event=event)
        errors = []
        for topology in ("hexane_only", "two_liquid"):
            try:
                return solve(self, requested_dt, event=event, topology=topology)
            except RuntimeError as exc:
                errors.append(
                    dict(topology=topology, error=str(exc), details=getattr(exc, "details", None))
                )
        failure = RuntimeError(f"radial A physical candidates refused: {errors}")
        failure.details = errors
        raise failure

    def liquid_water_u(self, temperature):
        return model.core.wa.state_Tp(self.ti, self.pressure, "liquid").u_mass

    def interface(self, temperature, total_mass_flux, clock):
        native, ack, residual, out, tb, yb, area = super().interface(
            temperature, total_mass_flux, clock
        )
        surface_u = self.liquid_water_u(self.ti) * model.core.wa.M
        residual += (
            area
            * native.water_molar_flux_mol_m2_s
            * (surface_u - native.native_water_liquid_internal_energy_j_mol)
        )
        return (
            replace(
                native,
                interface_energy_residual_w=residual,
                native_water_liquid_internal_energy_j_mol=surface_u,
            ),
            ack,
            residual,
            out,
            tb,
            yb,
            area,
        )


class FreshMatchedFeed:
    def __init__(self, out, stop_file, *, geometry, main_shells):
        self.out = out
        self.adapter = CaseAdapter(
            out, stop_file, pd_heating="constant", water_film_energy="surface"
        )
        self.adapter.codec = FEED_CODEC
        self.geometry = geometry
        self.main_shells = main_shells

    def run(self, settings):
        if self.main_shells < 2 or settings.cells % self.main_shells:
            raise ValueError("main-A shells must subdivide the B/C master grid exactly")
        declared_bath = json.loads(Path(settings.bath).read_text())
        declared_reference = declared_bath.get("pd_reference", asdict(model.PdReference()))
        if declared_reference != asdict(model.PdReference()):
            raise ValueError("this source snapshot requires its exact declared PD reference")
        pd_end = model.PdReference().residence_s
        pd_settings = SimpleNamespace(**(vars(settings) | {"max_time_s": pd_end}))
        original, previous_geometry = model.PrescribedPdFeed, MatchedPD.geometry
        MatchedPD.geometry = self.geometry
        try:
            model.PrescribedPdFeed = MatchedPD
            with transport.installed_pd_diffusion():
                feed = self.adapter.run_feed(pd_settings)
        finally:
            model.PrescribedPdFeed, MatchedPD.geometry = original, previous_geometry
        pd = self.adapter.members[0]
        if feed["errors"] or pd.time != pd_end or pd.free_water != 0.0:
            raise ValueError("fresh PD did not reach the declared admissible exit")
        feed.update(
            pd_transport_mode="actual_native_rho_cp_nist_wilke_cantera_fuller_original",
            pd_diffusion_geometry=self.geometry,
            pd_transfer_pair=transfer_pair_record(pd.pd_pair),
            pd_exposed_area_fraction=pd.pd_contact_fraction,
            pd_exposure_definition="top tray area divided by total sphere area: R/(3*(1-eps)*H)",
            pd_boundary_codec=self.adapter.codec,
        )
        feed["rows"][0].update(
            bath_temperature_k=feed["pd_source"]["temperature_k"],
            bath_y_hexane=feed["pd_source"]["y_hexane"],
        )
        feed["variant_source_sha256"].update(
            {
                path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (
                    Path(__file__),
                    Path(transport.__file__),
                    Path(__file__).with_name("particle_a_component_elimination.py"),
                )
            }
        )
        self.adapter.write("fresh_pd_feed.json", feed)
        last = feed["rows"][-1]
        bath = json.loads(Path(settings.bath).read_text())
        bath["particle"].update(
            feed_temperature_K=last["temperature_k"],
            feed_hexane_kg_kg_dry=last["hexane_kg_kg_dry"],
            feed_water_kg_kg_dry=last["retained_water_kg_kg_dry"],
        )
        march = MatchedRadialFeed(
            bath,
            settings.radius_m,
            settings.cells,
            shells=self.main_shells,
            property_parameters=model.ParticleModelParameters(),
            pressure_pa=feed["pd_source"]["reference"]["pressure_pa"],
        )
        march.time = pd_end
        march.two_liquid_trace = march.ti, march.yi
        # The native pure-hexane root was already independently solved for PD.
        march.pure_hexane_temperature_k = pd.pure_hexane_temperature_k
        march.rows = [march.row(None, None)]
        defect = abs(march.total_energy - last["total_particle_plus_film_energy_j"]) / (
            march.dry * 2e5
        )
        if defect > 1e-10:
            raise ValueError("uniform PD -> radial main-A storage import failed")
        failures = []
        activation = None
        while march.time < settings.max_time_s:
            self.adapter.check_stop()
            duration = min(settings.steam_dt_s, settings.max_time_s - march.time)
            try:
                march.solve(duration)
            except RuntimeError as error:
                failures.append(
                    dict(
                        time_s=march.time, error=str(error), details=getattr(error, "details", None)
                    )
                )
                # Exact absent-H stock is imposed IN the event solve, not clipped.
                march.solve(duration, event=True)
                activation = march.activate()
                break
        if activation is None:
            raise ValueError("fresh resolved main-A did not reach native radial activation")
        if activation["ledger"]["max_normalized_residual"] > 1e-10:
            raise ValueError("native fresh main-A activation ledger failed")
        if activation["radial"]["attached_hexane_mass_kg"] != 0.0:
            raise ValueError("native fresh main-A activation retains external hexane")
        feed.update(
            status="ACTUAL_FEED_TO_RADIAL_ACTIVATION",
            rows=[*feed["rows"], *march.rows[1:]],
            activation=activation,
            kernel_envelope_feed_cells=self.main_shells,
            high_loading_thermal_cells=dict(PD=1, main_A=self.main_shells),
            main_a_suffix_rows=march.rows,
            main_a_thermal_reduction="resolved main-tray A after a uniform PD body",
            pd_prefix_reused_checked=False,
            fresh_original_feed=True,
            interface_temperature_k=march.ti,
            interface_y_hexane=march.yi,
            external_water_temperature_k=march.ti,
            external_water_energy_temperature="surface",
            external_water_caloric="native surface water internal energy; consistent interface liquid donor",
            main_a_rejected_ordinary_trials=failures,
            uniform_to_radial_import_normalized_energy_defect=defect,
            surface_to_bulk_conductance_w_m2_k=march.internal_conductance_w_m2_k,
        )
        self.adapter.write(
            "fresh_main_a_summary.json",
            dict(
                time_s=march.time,
                shells=self.main_shells,
                mean_temperature_k=march.mass_mean_temperature_k,
                shell_temperature_bounds_k=[min(march.temperatures), max(march.temperatures)],
                external_water_kg_kg_dry=march.free_water / march.dry,
                maximum_main_A_normalized_ledger=max(
                    max(r["normalized_cumulative_ledger"]) for r in march.rows[1:]
                ),
            ),
        )
        return feed
