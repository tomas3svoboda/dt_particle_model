"""Core2 uniform high-loading DT experiment on the inherited PD/steam bath.

Reuses native sphere calorics, K2 phase/Stefan rows and signed PHY-053
heat capacities. Internal heat conduction is the native ONE-CELL Fourier
discretization, A*k/(R-r_center), not an empirical internal conductance.
Exterior water shares body temperature as in the native high-loading packet.
This is an explicitly coarse pre-activation study, not a full destiny result.
"""

from dataclasses import asdict, replace
import argparse
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import least_squares

import dtdc_particle_destiny_march as core
from dtdc_simulator.core2 import cell_native_caloric as nc
from dtdc_simulator.core2 import qsc_k2_native_source_kernel as kernel
from dtdc_simulator.core2.particle.surface_active_set import HexaneAccessOracle
from era_stamp import REPO_ROOT, era_stamp


class FeedMarch:
    def __init__(self, bath, radius, cells):
        self.bath = bath
        feed = bath["particle"]
        self.wet = replace(
            core.wet_core.WetCoreParams(),
            pressure_pa=core.PRESSURE_PA,
            X_water=feed["feed_water_kg_kg_dry"],
            w_o=feed["feed_oil_kg_kg_dry"],
            conductivity=core.CONDUCTIVITY,
        )
        self.params = core.sphere.SphereParams(
            particle_radius_m=radius, radial_cells=cells, wet_core=self.wet
        )
        self.state = core.sphere.initialize_qualified_feed(
            self.params,
            total_hexane_loading=feed["feed_hexane_kg_kg_dry"],
            temperature_k=feed["feed_temperature_K"],
            retained_water_loading=self.wet.X_water,
            oil_fraction=self.wet.w_o,
        )
        self.radius = radius
        self.area = 4 * math.pi * radius**2
        self.dry = self.state.dry_meal_mass_kg
        self.free_water = 0.0
        self.time = 0.0
        self.pair = core.film_pair_at_radius(radius)
        thermal_grid = core.q.backend.uniform_grid(1, radius)
        self.internal_conductance_w_m2_k = core.CONDUCTIVITY / (radius - thermal_grid.centers[-1])
        self.pd_end = 490.20224190324825
        self.pd_contact_fraction = 4.081e-3 * radius / 0.0015
        self.ti, self.yi = kernel._transparent_phase_root(
            pressure_pa=core.PRESSURE_PA, continuation=(334.5, 0.8)
        )
        self.codec_view = SimpleNamespace(
            component_datum_adapter=SimpleNamespace(
                hexane_native_to_common_offset_j_mol=0.0, water_native_to_common_offset_j_mol=0.0
            )
        )
        self.model_view = SimpleNamespace(
            properties=nc.NativeGasProperties(energy_datum_id=nc.NATIVE_ZERO_GAUGE_ENERGY_DATUM_ID),
            transfer=SimpleNamespace(
                wall_to_interface_ua_w_k=0.0,
                solid_to_interface_ua_w_k=self.area * self.internal_conductance_w_m2_k,
            ),
        )
        self.total_energy = self.state.total_internal_energy_j
        self.initial = (self.state.total_hexane_mass_kg, 0.0, self.total_energy)
        self.cumulative = np.zeros(3)
        self.rows = [self.row(None, None)]

    def interface(self, temperature, total_mass_flux, clock):
        in_pd = self.time < self.pd_end
        tb, yb = (
            (345.15, 0.6644429974143522)
            if in_pd
            else core.interpolate_bath(self.bath["waypoints"], clock)
        )
        contact = self.pd_contact_fraction if in_pd else 1.0
        active_area = self.area * contact
        pair = self.pair
        kappa = core.PRESSURE_PA / (core.bg.R * tb) * pair.binary_mass_transfer_m_s
        coefficients = kernel.NativeFilmCoefficients(
            3.9,
            pair.state.gas_density_kg_m3,
            pair.state.gas_superficial_velocity_m_s,
            pair.state.bulk_binary_diffusivity_m2_s,
            pair.reynolds_voidage,
            pair.prandtl,
            pair.schmidt,
            pair.nusselt_voidage,
            pair.heat_transfer_w_m2_k,
            pair.binary_mass_transfer_m_s,
            kappa,
        )
        # Pure-kernel field carriers: these are NOT host-issued accepted packets.
        basis = SimpleNamespace(
            request=SimpleNamespace(hexane_access_oracle=HexaneAccessOracle.TRANSPARENT_INTERFACE),
            complete=SimpleNamespace(state=SimpleNamespace(temperature_k=float(temperature))),
            packets=(SimpleNamespace(pressure_pa=core.PRESSURE_PA),),
            active_area_m2=active_area,
            dry_mass_fraction=contact,
        )
        native = kernel._group_trial(
            basis,
            model=self.model_view,
            coefficients=coefficients,
            bulk_temperature_k=tb,
            bulk_hexane_mole_fraction=yb,
            local_hydraulic_pressure_pa=core.PRESSURE_PA,
            shared_wall_temperature_k=self.ti,
            interface_temperature_k=self.ti,
            interface_hexane_mole_fraction=self.yi,
            total_molar_flux_mol_m2_s=0.0,
            codec=self.codec_view,
        )
        boundary = core.reduced_film_boundary(
            temperature_k=tb,
            y_hexane=yb,
            pair=pair,
            forcing_timescale_s=0.075,
            label="unchanged prescribed DT bath",
        )
        ref = core.ct.prepare_reduced_film_surface_state(
            boundary,
            cell_temperature_k=float(temperature),
            pressure_pa=core.PRESSURE_PA,
            particle_radius_m=self.radius,
            outer_half_cell_distance_m=0.5 * self.radius,
            particle_thermal_conductivity_w_m_k=core.CONDUCTIVITY,
            binary_gas_interaction_k_wh=boundary.binary_gas_interaction_k_wh,
        )
        film = core.ct._reduced_film_coefficients_at_surface_temperature(ref, self.ti)
        capacities = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
            surface_temperature_k=self.ti,
            bulk_temperature_k=tb,
            pressure_pa=core.PRESSURE_PA,
            y_hexane=self.yi,
        )
        fluxes = core.efr.component_fluxes_for_supplied_total_mass_flux(
            film,
            total_mass_flux_kg_m2_s=total_mass_flux,
            interface_hexane_mass_fraction=core._hexane_mass_fraction(self.yi),
            bulk_hexane_mass_fraction=core._hexane_mass_fraction(yb),
        )
        nh = fluxes.hexane_kg_m2_s / core.hx.M
        nw = fluxes.water_kg_m2_s / core.wa.M
        native = replace(
            native,
            hexane_molar_flux_mol_m2_s=nh,
            water_molar_flux_mol_m2_s=nw,
            total_molar_flux_mol_m2_s=nh + nw,
            raw_hexane_molar_flux_mol_m2_s=nh,
        )
        ack = core.efr.ackermann_heat_transfer(
            film,
            fluxes=fluxes,
            water_film_heat_capacity_j_kg_k=capacities.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=capacities.hexane_j_kg_k,
            enforce_production_domain=False,
        )
        qgas = active_area * film.heat_transfer_w_m2_k * ack.factor * (tb - self.ti)
        phase_energy = active_area * math.fsum(
            (
                nh
                * (
                    native.native_hexane_partial_enthalpy_j_mol
                    - native.native_hexane_liquid_internal_energy_j_mol
                ),
                nw
                * (
                    native.native_water_partial_enthalpy_j_mol
                    - native.native_water_liquid_internal_energy_j_mol
                ),
            )
        )
        interface_energy = qgas + native.particle_to_interface_w - phase_energy
        native = replace(
            native,
            gas_to_interface_w=qgas,
            interface_energy_residual_w=interface_energy,
            interface_energy_residual_scaled=interface_energy
            / max(abs(qgas), abs(native.particle_to_interface_w), abs(phase_energy), 1.0),
        )
        energy_out = (
            active_area
            * (
                native.hexane_molar_flux_mol_m2_s * native.native_hexane_partial_enthalpy_j_mol
                + native.water_molar_flux_mol_m2_s * native.native_water_partial_enthalpy_j_mol
            )
            - qgas
        )
        source = 85.0 * self.dry if in_pd else 0.0
        return native, ack, interface_energy, energy_out - source, tb, yb, active_area

    def solve(self, requested_dt, event=False):
        old = self.state
        old_h = old.total_hexane_loading
        old_w = self.free_water / self.dry
        old_u = self.total_energy
        history = {}

        def evaluate(x):
            temperature = float(x[0])
            attached = 0.0 if event else float(x[1])
            free_w = float(x[1] if event else x[2])
            total_flux = float(x[2] if event else x[3])
            dt = float(x[3]) if event else requested_dt
            loading = core.sphere.critical_loading(temperature, self.params) + attached
            state = core.sphere.initialize_qualified_feed(
                self.params,
                total_hexane_loading=loading,
                temperature_k=temperature,
                retained_water_loading=old.retained_water_loading,
                oil_fraction=old.oil_fraction,
            )
            energy = (
                state.total_internal_energy_j
                + self.dry
                * free_w
                * core.wa.state_Tp(temperature, core.PRESSURE_PA, "liquid").u_mass
            )
            native, ack, interface_e, energy_out, tb, yb, active_area = self.interface(
                temperature, total_flux, self.time + dt
            )
            residual = np.array(
                [
                    loading
                    - old_h
                    + dt * active_area / self.dry * native.hexane_molar_flux_mol_m2_s * core.hx.M,
                    free_w
                    - old_w
                    + dt * active_area / self.dry * native.water_molar_flux_mol_m2_s * core.wa.M,
                    (energy - old_u + dt * energy_out) / (self.dry * 2e5),
                    interface_e / (self.area * 1e4),
                ]
            )
            history.update(
                state=state,
                water=self.dry * free_w,
                energy=energy,
                dt=dt,
                native=native,
                ack=ack,
                energy_out=energy_out,
                residual=residual,
                tb=tb,
                yb=yb,
                active_area=active_area,
            )
            return residual

        attached0 = old.attached_hexane_mass_kg / self.dry
        x0 = (
            [old.temperature_k, old_w + 1e-4, 0.2, requested_dt * 0.5]
            if event
            else [old.temperature_k, attached0, old_w + 1e-4, 0.2]
        )
        lower = (
            [self.wet.T_min + 1e-4, 0.0, -100.0, 1e-9]
            if event
            else [self.wet.T_min + 1e-4, 0.0, 0.0, -100.0]
        )
        upper = (
            [self.wet.T_max - 1e-4, 1.0, 100.0, requested_dt]
            if event
            else [self.wet.T_max - 1e-4, 1.0, 1.0, 100.0]
        )
        solution = least_squares(
            evaluate,
            x0,
            bounds=(lower, upper),
            x_scale="jac",
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
            max_nfev=150,
            jac="3-point",
        )
        residual = evaluate(solution.x)
        if not solution.success or max(abs(residual)) > 1e-10:
            error = RuntimeError(
                f"native high-loading {'event' if event else 'step'} not closed: {max(abs(residual)):.6g}"
            )
            error.details = {
                "residual": residual.tolist(),
                "coordinates": solution.x.tolist(),
                "nfev": solution.nfev,
            }
            raise error
        if abs(history["ack"].beta) > 3.0:
            raise RuntimeError(
                f"closed trial outside declared signed Ackermann study band: {history['ack'].beta}"
            )
        if event and history["native"].hexane_molar_flux_mol_m2_s <= 0:
            raise RuntimeError("attached-film depletion event lacks outward hexane flux")
        # No positive film is imported: the appearance step must condense water.
        if (
            self.free_water == 0
            and history["water"] > 0
            and history["native"].water_molar_flux_mol_m2_s >= 0
        ):
            raise RuntimeError("water appearance lacks inward component flux")
        dt = history["dt"]
        increments = np.array(
            [
                dt
                * history["active_area"]
                * history["native"].hexane_molar_flux_mol_m2_s
                * core.hx.M,
                dt
                * history["active_area"]
                * history["native"].water_molar_flux_mol_m2_s
                * core.wa.M,
                dt * history["energy_out"],
            ]
        )
        trial_cumulative = self.cumulative + increments
        actual = np.array(
            [history["state"].total_hexane_mass_kg, history["water"], history["energy"]]
        )
        ledger = actual - np.array(self.initial) + trial_cumulative
        ledger_scales = np.array(
            [
                max(self.initial[0], actual[0], abs(trial_cumulative[0])),
                self.dry * old.retained_water_loading + max(actual[1], abs(trial_cumulative[1])),
                self.dry * 2e5,
            ]
        )
        normalized = np.abs(ledger) / ledger_scales
        if np.max(normalized) > 1e-10:
            raise RuntimeError(
                f"whole particle/film cumulative ledger refuses before commit: {normalized.tolist()}"
            )
        history["normalized_cumulative_ledger"] = normalized.tolist()
        self.cumulative = trial_cumulative
        self.state = history["state"]
        self.free_water = history["water"]
        self.total_energy = history["energy"]
        self.time += dt
        self.rows.append(self.row(history, ledger))
        return history

    def row(self, step, ledger):
        return {
            "t_s": self.time,
            "temperature_k": self.state.temperature_k,
            "hexane_kg_kg_dry": self.state.total_hexane_loading,
            "attached_hexane_kg_kg_dry": self.state.attached_hexane_mass_kg / self.dry,
            "retained_water_kg_kg_dry": self.state.retained_water_loading,
            "external_water_kg_kg_dry": self.free_water / self.dry,
            "total_particle_plus_film_energy_j": self.total_energy,
            "beta": step["ack"].beta if step else None,
            "cumulative_component_energy_residuals": (
                ledger.tolist() if ledger is not None else None
            ),
            "normalized_cumulative_ledger": (
                step.get("normalized_cumulative_ledger") if step else None
            ),
            "bath_temperature_k": step["tb"] if step else None,
            "bath_y_hexane": step["yb"] if step else None,
        }


def run(args):
    bath = json.loads(args.bath.read_text(encoding="utf-8"))
    entry_stamp = era_stamp([Path(__file__).resolve(), args.bath.resolve()])
    march = FeedMarch(bath, args.radius_m, args.cells)
    errors = []
    activation = None
    while march.time < args.max_time_s:
        try:
            step_dt = args.dt_s if march.time < march.pd_end else args.steam_dt_s
            dt = min(step_dt, args.max_time_s - march.time)
            if march.time < march.pd_end:
                dt = min(dt, march.pd_end - march.time)
            march.solve(dt)
        except Exception as step_error:
            errors.append(
                {
                    "at_s": march.time,
                    "step_error": str(step_error),
                    "details": getattr(step_error, "details", None),
                }
            )
            try:
                march.solve(dt, event=True)
                radial, ledger = core.sphere.activate_radially(march.state, march.params)
                activation = {
                    "radial": asdict(radial),
                    "ledger": asdict(ledger),
                    "external_water_mass_kg_carried_separately": march.free_water,
                    "external_water_energy_j_carried_separately": march.total_energy
                    - march.state.total_internal_energy_j,
                }
            except Exception as event_error:
                errors[-1].update(
                    event_error=str(event_error),
                    event_details=getattr(event_error, "details", None),
                )
            break
    return {
        "status": (
            "ACTUAL_FEED_TO_RADIAL_ACTIVATION"
            if activation
            else ("ACTUAL_FEED_MARCH_REFUSED" if errors else "ACTUAL_FEED_PARTIAL_HORIZON")
        ),
        "complete_dt_destiny": False,
        "physically_qualifying": False,
        "bath_was_modified": False,
        "boundary_history": "inherited PD dome/deck/exposed-area, then original flash/toaster waypoints",
        "phase_branch": "native condensing free-film appearance / transparent two-liquid interface",
        "external_film_law": "native EFr mass-coordinate Stefan and full signed-component PHY053 Ackermann, shared with radial boundary",
        "dt_s": args.dt_s,
        "radius_m": args.radius_m,
        "high_loading_thermal_cells": 1,
        "surface_to_bulk_conductance_w_m2_k": march.internal_conductance_w_m2_k,
        "external_water_caloric": "native liquid internal energy at body temperature",
        "interface_temperature_k": march.ti,
        "interface_y_hexane": march.yi,
        "rows": march.rows,
        "activation": activation,
        "errors": errors,
        "entry_source_stamp": entry_stamp,
        "stamp": era_stamp(
            [
                Path(__file__),
                args.bath,
                Path(kernel.__file__),
                Path(core.sphere.__file__),
                Path(core.efr.__file__),
                Path(core.bg.__file__),
                Path(core.hx.__file__),
                Path(core.wa.__file__),
            ]
        ),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--bath", type=Path, default=core.DESTINY_PATH.with_name("destiny_v3_exit105.json")
    )
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--dt-s", type=float, default=0.05)
    p.add_argument("--steam-dt-s", type=float, default=0.02)
    p.add_argument("--max-time-s", type=float, default=10)
    p.add_argument("--radius-m", type=float, default=0.000885)
    p.add_argument("--cells", type=int, default=2)
    args = p.parse_args()
    if args.out.exists() or not args.out.resolve().is_relative_to(REPO_ROOT / "paper"):
        p.error("use a new output under paper")
    result = run(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {"status": result["status"], "last": result["rows"][-1], "errors": result["errors"]},
            indent=2,
        )
    )
