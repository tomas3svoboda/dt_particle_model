"""Declared preheated Faner charge -> native activation -> binary radial birth.

The attached-hexane stage uses the existing pure-component Stefan limit.
Its interface energy and the native particle inventories close total flow.
This is a prescribed-bath numerical experiment, not a reconstructed measured
preheating history or a qualified prediction of the missing initial moisture.
"""

from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path
import argparse
import hashlib
import json
import math
import pickle
import sys
import time
import traceback

import numpy as np
from scipy.optimize import brentq, least_squares

import experiment_core2_dt_birth as shared
import finite_film_core2_study as study
import positive_water_continuation as continued
from era_stamp import REPO_ROOT

FANER_HELPERS = REPO_ROOT / "scratchpad/particle_destiny_2026-10-06/faner_recovery/active_faner"
sys.path.insert(0, str(FANER_HELPERS))
from attempt_active_faner import active_binary_controls  # noqa: E402
from bulk_endpoint import (  # noqa: E402
    INHERITED_FANER_COMMON_EXPOSURE_MULTIPLIER,
    declared_faner_film_multiplier,
    make_pure_bulk_boundary,
    pure_bulk_domain_adapter,
)

core = study.core


@dataclass(frozen=True)
class FanerStudyControls(core.ci.CutSolverControls):
    """Declared wider numerical chart; original acceptance fields validated.

    The native constructor additionally caps a solver coordinate at 35.
    The owner's RULE 2 authorizes this separately named study range through
    70. The validation probe advances nothing and is never used in a solve.
    """

    def __post_init__(self):
        if (
            not math.isfinite(self.maximum_logit_magnitude)
            or not 10 <= self.maximum_logit_magnitude <= 70
        ):
            raise ValueError("declared Faner coordinate magnitude must lie in [10,70]")
        if (
            not math.isfinite(self.nonlinear_residual_tolerance)
            or not 0 < self.nonlinear_residual_tolerance <= 1e-9
        ):
            raise ValueError("declared local numerical residual tolerance cannot exceed 1e-9")
        kwargs = {
            f.name: getattr(self, f.name) for f in fields(core.ci.CutSolverControls) if f.init
        }
        kwargs["maximum_logit_magnitude"] = min(self.maximum_logit_magnitude, 35.0)
        kwargs["nonlinear_residual_tolerance"] = min(self.nonlinear_residual_tolerance, 1e-10)
        core.ci.CutSolverControls(**kwargs)


def wider_numerical_controls(controls, magnitude, regional_tolerance=None, *, extended=False):
    kwargs = {
        f.name: getattr(controls, f.name) for f in fields(core.ci.CutSolverControls) if f.init
    }
    if regional_tolerance is not None:
        kwargs["nonlinear_residual_tolerance"] = regional_tolerance
    kind = continued.RegionalRootStudyControls if extended else FanerStudyControls
    return kind(**{**kwargs, "maximum_logit_magnitude": magnitude})


def leading_birth_seed(before, config, controls, boundary):
    """Native zero-thickness component/energy jumps as a seed, never a march.

    The original wet donor is retained, while the limiting dry state is the
    native saturated-H interface. A failed seed root is not a physical state.
    """
    old = shared.ww._evaluate_state_cells(before)[-1]
    distance = before.grid.R - before.grid.centers[-1]
    last = {}

    def residual(x):
        t, beta = map(float, x)
        interface = core.cut.evaluate_interface_state(
            t, config, before.historical_hexane_loadings[-1], before.oil_fraction_labels[-1]
        )
        wet_flux = core.cut._wet_flux(
            old.temperature_k,
            old.retained_water_loading,
            t,
            interface.wet.retained_water_loading,
            distance,
            config,
            "wet cell",
            "interface",
            left_capacity_dual_over_rt=old.retained_water_capacity_dual_over_rt,
            right_capacity_dual_over_rt=interface.retained_water_trace_capacity_dual_over_rt,
        )
        _, components, ack, film = shared.film_at_beta(
            boundary, config, before.grid.R, t, beta, surface_y_hexane=interface.y_hexane
        )
        jump_h = (
            old.total_hexane_concentration_kg_m3 / core.hx.M
            - interface.dry.total_hexane_concentration_mol_m3
        )
        speed = components.hexane_kg_m2_s / core.hx.M / jump_h
        qout = film.heat_transfer_w_m2_k * ack.factor * (t - boundary.temperature_k)
        energy_out = math.fsum(
            (
                qout,
                components.hexane_kg_m2_s
                / core.hx.M
                * interface.dry.hexane_gas_partial_enthalpy_j_mol,
                components.water_kg_m2_s
                / core.wa.M
                * interface.dry.water_gas_partial_enthalpy_j_mol,
            )
        )
        water_r = (
            components.water_kg_m2_s / core.wa.M
            - wet_flux.retained_water_flux_mol_m2_s
            - speed
            * (
                old.retained_water_concentration_mol_m3
                - interface.dry.total_water_concentration_mol_m3
            )
        )
        energy_r = (
            energy_out
            - wet_flux.total_energy_flux_w_m2
            - speed * (old.energy_density_j_m3 - interface.dry.energy_density_j_m3)
        )
        last.update(
            interface_temperature_k=t,
            beta=beta,
            speed_m_s=speed,
            water_residual_mol_m2_s=water_r,
            energy_residual_w_m2=energy_r,
            hexane_out_kg_m2_s=components.hexane_kg_m2_s,
            water_out_kg_m2_s=components.water_kg_m2_s,
            wet_water_out_mol_m2_s=wet_flux.retained_water_flux_mol_m2_s,
            interface_retained_water_loading=interface.wet.retained_water_loading,
        )
        return np.array([water_r / 0.02, energy_r / 5000.0])

    t0 = shared.interface_water_seed(before, config, controls)[0] + 0.1
    lo, hi = controls.interface_temperature_bounds_k
    solution = least_squares(
        residual,
        [t0, 1.0],
        bounds=([lo + 1e-6, 0.01], [hi - 1e-6, 3.0]),
        jac="3-point",
        x_scale="jac",
        max_nfev=150,
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    r = residual(solution.x)
    last.update(
        seed_only=True,
        actual_time_advanced_s=0.0,
        success=bool(solution.success and max(abs(r)) < 1e-9 and last["speed_m_s"] > 0),
        maximum_scaled_residual=float(max(abs(r))),
    )
    return last


def leading_pore_profile(before, config, boundary, leading):
    """Shoot the native pore face against the limiting external flux, seed only."""
    gamma = leading["interface_temperature_k"]
    interface = core.cut.evaluate_interface_state(
        gamma, config, before.historical_hexane_loadings[-1], before.oil_fraction_labels[-1]
    )
    pore = replace(config.dry.pore, w_o=before.oil_fraction_labels[-1])
    _, components, ack, film = shared.film_at_beta(
        boundary, config, before.grid.R, gamma, leading["beta"], surface_y_hexane=interface.y_hexane
    )
    qgas = film.heat_transfer_w_m2_k * ack.factor * (boundary.temperature_k - gamma)
    targets = np.array(
        [
            components.water_kg_m2_s / core.wa.M,
            components.hexane_kg_m2_s / core.hx.M,
            math.fsum(
                (
                    components.water_kg_m2_s
                    / core.wa.M
                    * interface.dry.water_gas_partial_enthalpy_j_mol,
                    components.hexane_kg_m2_s
                    / core.hx.M
                    * interface.dry.hexane_gas_partial_enthalpy_j_mol,
                    -qgas,
                )
            ),
        ]
    )
    distance = 1e-6
    last = {}

    def residual(x):
        g_t, g_y, stefan = map(float, x)
        face = core.cut._dry_flux(
            gamma,
            interface.y_hexane,
            gamma + distance * g_t,
            interface.y_hexane + distance * g_y,
            distance,
            stefan,
            pore,
            pore,
            config,
            thermodynamic_force_temperature_k=gamma,
        )
        actual = np.array(
            [
                face.component.conserved_water_flux_mol_m2_s,
                face.component.conserved_hexane_flux_mol_m2_s,
                face.energy.total_energy_flux_w_m2,
            ]
        )
        last.update(
            temperature_gradient_k_m=g_t,
            composition_gradient_m_inv=g_y,
            total_stefan_flux_mol_m2_s=stefan,
            native_face_fluxes=actual.tolist(),
        )
        return (actual - targets) / np.array([0.02, 0.1, 5000.0])

    linear_diagnostics = []
    for probe_distance in (1e-5, 1e-6, 1e-7):
        distance = probe_distance
        base = np.array([qgas / config.dry.thermal_conductivity.value_w_m_k, 0.0, sum(targets[:2])])
        base_r = residual(base)
        columns = []
        for column, h in enumerate((0.1, 0.1, 1e-5)):
            right, left = base.copy(), base.copy()
            right[column] += h
            left[column] -= h
            columns.append((residual(right) - residual(left)) / (2 * h))
        jacobian = np.column_stack(columns)
        required = base - np.linalg.solve(jacobian, base_r)
        right_t = gamma + distance * required[0]
        right_y = interface.y_hexane + distance * required[1]
        activity = None
        if 0 < right_y < 1:
            gas = core.bg.state(right_t, config.dry.pressure_pa, 1 - right_y, right_y)
            activity = (
                gas.fugacity_hexane_pa
                / core.hx.state_Tp(right_t, config.dry.pressure_pa, "liquid").fugacity
            )
        saturation_slope = (
            core.cp.gas_only_composition_interval(
                gamma + 1e-3, config.dry.pressure_pa, pore
            ).upper_y_hexane
            - core.cp.gas_only_composition_interval(
                gamma - 1e-3, config.dry.pressure_pa, pore
            ).upper_y_hexane
        ) / 2e-3
        linear_diagnostics.append(
            {
                "probe_distance_m": distance,
                "required_linearized_gradient": required.tolist(),
                "hexane_saturation_curve_dy_dt": saturation_slope,
                "composition_gradient_minus_saturation_tangent": float(
                    required[1] - saturation_slope * required[0]
                ),
                "predicted_endpoint_hexane_activity": activity,
                "predicted_endpoint_y_hexane": float(right_y),
                "all_jacobian_probes_native_admissible": True,
                "extrapolated_state_accepted": False,
            }
        )
    distance = 1e-6
    try:
        solution = study.numerics.domain_aware_newton(
            residual,
            [qgas / config.dry.thermal_conductivity.value_w_m_k, 0.0, sum(targets[:2])],
            bounds=([-1e7, -1e6, -10], [1e7, 1e6, 10]),
            max_nfev=600,
            gtol=1e-8,
        )
        values = residual(solution.x)
        last.update(
            success=bool(solution.success),
            message=str(solution.message),
            maximum_scaled_residual=float(max(abs(values))),
        )
    except Exception as error:
        last.update(
            success=False, error=str(error), optimizer_debug=getattr(error, "optimizer_debug", None)
        )
    last.update(
        seed_only=True,
        actual_time_advanced_s=0.0,
        probe_distance_m=distance,
        required_face_fluxes=targets.tolist(),
        linear_response_diagnostics=linear_diagnostics,
    )
    return last


class FanerAttachedMarch:
    """Existing HEXANE_ONLY surface selection with native high-loading storage.

    Retained water remains a conserved material component. The attached pure
    hexane film selects zero external water flux before its depletion; exposed
    pore water is solved dynamically by the shared binary node after birth.
    """

    def __init__(self, args):
        self.wet = replace(core.wet_core.WetCoreParams(), X_water=args.water, w_o=args.oil)
        self.params = core.sphere.SphereParams(
            particle_radius_m=args.radius, radial_cells=args.cells, wet_core=self.wet
        )
        self.area = 4 * math.pi * args.radius**2
        self.dry = self.params.dry_meal_mass_kg
        self.pair = core.film_pair_at_radius(args.radius)
        self.boundary = make_pure_bulk_boundary(
            self.pair,
            temperature_k=args.bath_temperature,
            forcing_timescale_s=args.forcing_timescale,
            beta_limit=3.0,
            stress_multiplier=args.film_multiplier,
        )
        self.interface_temperature = brentq(
            lambda t: math.log(
                core.bg.state(t, core.PRESSURE_PA, 0.0, 1.0).fugacity_hexane_pa
                / core.hx.state_Tp(t, core.PRESSURE_PA, "liquid").fugacity
            ),
            330.0,
            350.0,
            xtol=1e-11,
        )
        initial_t = (
            self.interface_temperature
            if args.initial_temperature is None
            else args.initial_temperature
        )
        self.state = core.sphere.initialize_qualified_feed(
            self.params,
            total_hexane_loading=args.initial_hexane,
            temperature_k=initial_t,
            retained_water_loading=args.water,
            oil_fraction=args.oil,
        )
        self.time = 0.0
        self.cumulative = np.zeros(3)
        self.initial = np.array(
            [
                self.state.total_hexane_mass_kg,
                self.state.retained_water_mass_kg,
                self.state.total_internal_energy_j,
            ]
        )
        # Same native one-cell Fourier approximation as the coarse DT feed.
        grid = core.backend.uniform_grid(1, args.radius)
        self.conductance = self.wet.conductivity / (args.radius - grid.centers[-1])
        self.rows = [self.row(None)]

    def rates(self, body_t, total_mass):
        ts, boundary = self.interface_temperature, self.boundary
        reference = core.ct.prepare_reduced_film_surface_state(
            boundary.film_boundary,
            cell_temperature_k=body_t,
            pressure_pa=core.PRESSURE_PA,
            particle_radius_m=self.params.particle_radius_m,
            outer_half_cell_distance_m=0.5 * self.params.particle_radius_m,
            particle_thermal_conductivity_w_m_k=self.wet.conductivity,
            binary_gas_interaction_k_wh=boundary.film_boundary.binary_gas_interaction_k_wh,
        )
        film = core.ct._reduced_film_coefficients_at_surface_temperature(reference, ts)
        components = core.efr.component_fluxes_for_supplied_total_mass_flux(
            film,
            total_mass_flux_kg_m2_s=total_mass,
            interface_hexane_mass_fraction=1.0,
            bulk_hexane_mass_fraction=1.0,
        )
        cp = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
            surface_temperature_k=ts,
            bulk_temperature_k=boundary.temperature_k,
            pressure_pa=core.PRESSURE_PA,
            y_hexane=1.0,
        )
        ack = core.efr.ackermann_heat_transfer(
            film,
            fluxes=components,
            water_film_heat_capacity_j_kg_k=cp.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=cp.hexane_j_kg_k,
            enforce_production_domain=False,
        )
        gas = core.bg.state(ts, core.PRESSURE_PA, 0.0, 1.0)
        hg = (core.hx.h_ideal(ts) + gas.partial_residual_enthalpy_hexane_molar) / core.hx.M
        ul = core.hx.state_Tp(body_t, core.PRESSURE_PA, "liquid").u_mass
        qgas = film.heat_transfer_w_m2_k * ack.factor * (boundary.temperature_k - ts)
        qbody = self.conductance * (body_t - ts)
        interface_residual = math.fsum((qgas, qbody, -components.hexane_kg_m2_s * (hg - ul)))
        return {
            "hexane_out_kg_s": self.area * components.hexane_kg_m2_s,
            "water_out_kg_s": self.area * components.water_kg_m2_s,
            "energy_out_w": self.area * (components.hexane_kg_m2_s * hg - qgas),
            "interface_energy_residual_w_m2": interface_residual,
            "beta": ack.beta,
            "interface_y_hexane": 1.0,
            "bulk_y_hexane": 1.0,
        }

    def solve(self, requested_dt, *, event=False):
        old = self.state
        latest = {}

        def residual(x):
            t = old.temperature_k + float(x[0])
            attached = 0.0 if event else float(x[1])
            dt = float(x[1]) if event else requested_dt
            rates = self.rates(t, float(x[2]))
            state = core.sphere.initialize_qualified_feed(
                self.params,
                total_hexane_loading=core.sphere.critical_loading(t, self.params) + attached,
                temperature_k=t,
                retained_water_loading=old.retained_water_loading,
                oil_fraction=old.oil_fraction,
            )
            r = np.array(
                [
                    (
                        state.total_hexane_mass_kg
                        - old.total_hexane_mass_kg
                        + dt * rates["hexane_out_kg_s"]
                    )
                    / self.dry,
                    (
                        state.total_internal_energy_j
                        - old.total_internal_energy_j
                        + dt * rates["energy_out_w"]
                    )
                    / (self.dry * 2e5),
                    rates["interface_energy_residual_w_m2"] / 1e4,
                ]
            )
            latest.update(state=state, dt=dt, rates=rates, residual=r)
            return r

        attached = old.attached_hexane_mass_kg / self.dry
        seed = [0.0, requested_dt * 0.5 if event else attached, 0.01]
        low = [self.wet.T_min - old.temperature_k + 1e-6, 1e-9 if event else 0.0, 0.0]
        high = [self.wet.T_max - old.temperature_k - 1e-6, requested_dt if event else 1.0, 0.1]
        result = least_squares(
            residual,
            seed,
            bounds=(low, high),
            jac="3-point",
            x_scale="jac",
            max_nfev=120,
            ftol=1e-13,
            xtol=1e-13,
            gtol=1e-13,
        )
        r = residual(result.x)
        if not result.success or max(abs(r)) > 1e-10:
            raise RuntimeError(
                f"attached {'event' if event else 'step'} residual {max(abs(r))}; x={result.x.tolist()}"
            )
        rates, dt, state = latest["rates"], latest["dt"], latest["state"]
        if (
            rates["water_out_kg_s"] != 0.0
            or rates["hexane_out_kg_s"] <= 0.0
            or abs(rates["beta"]) > 3
        ):
            raise RuntimeError("native pure-component surface branch not admissible")
        cumulative = self.cumulative + dt * np.array(
            [rates["hexane_out_kg_s"], rates["water_out_kg_s"], rates["energy_out_w"]]
        )
        actual = np.array(
            [
                state.total_hexane_mass_kg,
                state.retained_water_mass_kg,
                state.total_internal_energy_j,
            ]
        )
        ledger = actual - self.initial + cumulative
        normalized = abs(ledger) / np.array([self.initial[0], self.initial[1], self.dry * 2e5])
        if max(normalized) > 1e-10:
            raise RuntimeError(f"attached cumulative conservation refused: {normalized.tolist()}")
        self.state, self.cumulative, self.time = state, cumulative, self.time + dt
        latest["normalized_cumulative_ledger"] = normalized.tolist()
        self.rows.append(self.row(latest))
        return latest

    def row(self, step):
        return {
            "t_s": self.time,
            "temperature_k": self.state.temperature_k,
            "hexane_kg_kg_dry": self.state.total_hexane_loading,
            "attached_hexane_kg_kg_dry": self.state.attached_hexane_mass_kg / self.dry,
            "retained_water_kg_kg_dry": self.state.retained_water_loading,
            "external_water_kg_kg_dry": 0.0,
            "total_particle_energy_j": self.state.total_internal_energy_j,
            "dry_meal_mass_kg": self.dry,
            "outward_component_energy_rates": step["rates"] if step else None,
            "beta": step["rates"]["beta"] if step else None,
            "normalized_cumulative_ledger": step["normalized_cumulative_ledger"] if step else None,
        }


def continue_radial(checkpoint, surface, front_rate, record, args):
    """Advance only accepted native states; record and halve refused intervals."""
    before, config, carry, boundary = (
        checkpoint[name] for name in ("before", "config", "carry", "boundary")
    )
    dry_mass = (
        config.wet.wet.rho_dm_p * 4 / 3 * math.pi * before.transport.geometry.master_grid.R**3
    )
    dt = args.continuation_dt
    with (
        declared_faner_film_multiplier(boundary.film_boundary.stress_multiplier),
        pure_bulk_domain_adapter(),
        shared.bounded_path_proof(3000000, args.wall_budget),
        shared.path_bounds.monotone_hexane_path_enclosure([]),
        study.magnitude_storage_scale(),
    ):
        for index in range(args.steps):
            if before.transport.time_s >= args.horizon:
                record["requested_horizon_reached"] = True
                record["complete_falling_rate_experiment"] = args.horizon >= 240.0
                break
            trial_dt = min(dt, args.horizon - before.transport.time_s)
            committed = before
            while trial_dt >= args.minimum_dt:
                candidate = core.cut.candidate_from_state(
                    before.transport,
                    dry_total_stefan_fluxes_mol_m2_s=before.last_total_stefan_fluxes_mol_m2_s,
                    interface_temperature_k=before.last_interface_temperature_k,
                )
                candidate = replace(candidate, front_z=candidate.front_z - front_rate * trial_dt)
                seed = core.backend._seed(
                    candidate, "native Faner continuation from accepted birth"
                )
                pore = replace(config.dry.pore, w_o=before.transport.oil_fraction_labels[-1])
                lo, hi = core.cp.gas_only_composition_interval(
                    surface.temperature_k, config.dry.pressure_pa, pore
                ).y_hexane_bounds
                try:
                    step = shared.node.advance_water_node(
                        before,
                        trial_dt,
                        boundary,
                        seed,
                        old_loading=0.0,
                        old_temperature=surface.temperature_k,
                        surface_seed=[
                            surface.temperature_k,
                            surface.beta,
                            (surface.surface_y_hexane - lo) / (hi - lo),
                        ],
                        physical_solver=True,
                        residual_budget=args.budget,
                        gas_only_surface=True,
                        maximum_restart_trial_residual=0.05,
                        difference_mode=args.difference_mode,
                        beta_decoder=lambda t, b, y: shared.film_at_beta(
                            boundary,
                            config,
                            before.transport.geometry.master_grid.R,
                            t,
                            b,
                            surface_y_hexane=y,
                        )[0],
                        energy_context=lambda holder: shared.liquid_internal_energy_node(
                            holder, trial_dt, 0.0, surface.temperature_k, config.dry.pressure_pa
                        ),
                    )
                    audit = step.surface_film_audit
                    combined = shared.combined_ledger(before, config, carry, step, audit)
                    if combined["maximum_normalized_residual"] > before.controls.ledger_tolerance:
                        raise RuntimeError("Faner particle-plus-surface conservation refused")
                    front_rate = (
                        before.transport.geometry.front.z - step.after.transport.geometry.front.z
                    ) / trial_dt
                    before, surface = step.after, audit
                    carry = {**carry, **combined}
                    checkpoint.update(
                        before=before,
                        carry=carry,
                        surface=surface,
                        front_rate=front_rate,
                        source_record=str(args.out),
                    )
                    record["attempts"].append(
                        {
                            "accepted": True,
                            "stage": "radial",
                            "index": index,
                            "dt_s": trial_dt,
                            "after": core.receding_summary(before, dry_mass),
                            "ledger": asdict(step.ledger),
                            "surface": asdict(audit),
                            "combined_conservation": combined,
                        }
                    )
                    with args.out.with_suffix(".accepted.pkl").open("wb") as stream:
                        pickle.dump(checkpoint, stream)
                    record["last_accepted"] = core.receding_summary(before, dry_mass)
                    args.out.write_text(
                        json.dumps(record, indent=2, allow_nan=True), encoding="utf-8"
                    )
                    dt = min(args.continuation_dt, trial_dt * 1.5)
                    break
                except Exception as error:
                    record["attempts"].append(
                        {
                            "accepted": False,
                            "stage": "radial",
                            "index": index,
                            "at_s": before.transport.time_s,
                            "dt_s": trial_dt,
                            "before_identity_preserved": before is committed,
                            "error": str(error),
                            "debug": getattr(error, "study_debug", None),
                            "optimizer_debug": getattr(error, "optimizer_debug", None),
                        }
                    )
                    if "computational budget exhausted" in str(error):
                        record["stop"] = "computation budget reached; accepted checkpoint preserved"
                        return
                    trial_dt *= 0.5
            else:
                record["stop"] = "declared minimum interval reached after exact rollback"
                break
    if before.transport.time_s >= args.horizon:
        record["requested_horizon_reached"] = True
        record["complete_falling_rate_experiment"] = args.horizon >= 240.0


def run(args):
    sources = [
        Path(__file__),
        *sorted((REPO_ROOT / "src/dtdc_simulator/core2").rglob("*.py")),
        *sorted(FANER_HELPERS.glob("*.py")),
        *(
            REPO_ROOT / "paper/analysis/harness" / name
            for name in [
                "experiment_core2_dt_birth.py",
                "core2_water_film_node.py",
                "finite_film_core2_study.py",
                "particle_numerical_repairs.py",
                "native_water_path_bounds.py",
                "dtdc_particle_destiny_march.py",
                "positive_water_continuation.py",
            ]
        ),
    ]

    def hashes():
        return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}

    entry = hashes()
    record = {
        "experiment": "native Faner attached-film to active-water birth",
        "inputs": vars(args).copy(),
        "complete_falling_rate_experiment": False,
        "preheat_history_measured": False,
        "initial_hexane_basis": "declared kg/kg dry meal; source basis unspecified",
        "source_sha256_at_entry": entry,
        "attempts": [],
    }
    record["inputs"]["out"] = str(args.out)
    record["inputs"]["root_seed"] = str(args.root_seed) if args.root_seed else None
    record["inputs"]["resume"] = str(args.resume) if args.resume else None
    started = time.perf_counter()
    if args.resume:
        try:
            parent_path = args.resume.with_name(args.resume.name.replace(".accepted.pkl", ".json"))
            parent = json.loads(parent_path.read_text(encoding="utf-8"))
            if not parent["sources_unchanged_during_run"]:
                raise ValueError("parent sources changed during its run")
            approved_changes = []
            for path, expected in parent["source_sha256_at_exit"].items():
                if (
                    str(REPO_ROOT / "src/dtdc_simulator/core2") in path
                    and hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected
                ):
                    relative = (
                        Path(path).relative_to(REPO_ROOT / "src/dtdc_simulator/core2").as_posix()
                    )
                    if (
                        not args.continued_water
                        or relative not in continued.APPROVED_NATIVE_CHANGES
                    ):
                        raise ValueError(
                            "native Core2 source identity changed outside approved continuation"
                        )
                    approved_changes.append(
                        {"path": relative, "parent_sha256": expected, "current_sha256": entry[path]}
                    )
            with args.resume.open("rb") as stream:
                checkpoint = pickle.load(stream)
            before = checkpoint["before"]
            if (
                not isinstance(before, core.ci.CutIntegratorState)
                or before.transport.config is not checkpoint["config"]
            ):
                raise ValueError("checkpoint lost native radial state/config identity")
            controls = wider_numerical_controls(
                before.controls,
                args.coordinate_limit,
                args.regional_tolerance,
                extended=args.extended_regional_root,
            )
            numerical_after = replace(before, controls=controls)
            if numerical_after.transport is not before.transport:
                raise ValueError("numerical declaration changed the accepted physical state")
            checkpoint["before"] = before = numerical_after
            if args.continued_water:
                record["model_adoption"] = continued.adopt_checkpoint(checkpoint)
                record["approved_native_source_changes"] = approved_changes
                before = checkpoint["before"]
            record["numerical_contract"] = {
                "regional_residual_tolerance": controls.nonlinear_residual_tolerance,
                "step_cumulative_and_combined_ledger_tolerance": controls.ledger_tolerance,
                "physical_state_and_phase_laws_changed": False,
                "authorization": "owner RULE 2 explicit negligible-numerics allowance",
            }
            if "surface" in checkpoint:
                surface, front_rate = checkpoint["surface"], checkpoint["front_rate"]
            else:
                birth_audit = parent["attempts"][0]["audit"]
                surface = shared.node.WaterFilmNodeAudit(**birth_audit["surface_water_node"])
                front_rate = (1 - before.transport.geometry.front.z) / birth_audit[
                    "native_birth_ledger"
                ]["dt_s"]
            record["accepted_parent"] = {
                "result": str(parent_path),
                "checkpoint": str(args.resume),
                "checkpoint_sha256": hashlib.sha256(args.resume.read_bytes()).hexdigest(),
                "result_sha256": hashlib.sha256(parent_path.read_bytes()).hexdigest(),
            }
            record["effective_resumed_inputs"] = {
                "bulk_temperature_k": checkpoint["boundary"].temperature_k,
                "bulk_y_hexane": checkpoint["boundary"].y_hexane,
                "native_pore_binary_diffusivity_m2_s": checkpoint[
                    "config"
                ].dry.binary_diffusivity.value_m2_s,
                "retained_water_mobility_mol_m_s": checkpoint[
                    "config"
                ].dry.retained_water_mobility.value_mol_m_s,
                "parent_record": str(parent_path),
                "parsed_initialization_options_used_to_reset_state": False,
            }
            record["last_accepted"] = parent["last_accepted"]
            continue_radial(checkpoint, surface, front_rate, record, args)
        except Exception as error:
            record.update(
                error=str(error), error_type=type(error).__name__, traceback=traceback.format_exc()
            )
        record["wall_seconds"] = time.perf_counter() - started
        record["source_sha256_at_exit"] = hashes()
        record["sources_unchanged_during_run"] = entry == record["source_sha256_at_exit"]
        return record
    try:
        with declared_faner_film_multiplier(args.film_multiplier), pure_bulk_domain_adapter():
            feed = FanerAttachedMarch(args)
            record["initial_native_pure_hexane_phase_temperature_k"] = feed.interface_temperature
            while feed.time < args.horizon:
                dt = min(args.feed_dt, args.horizon - feed.time)
                try:
                    feed.solve(dt)
                except Exception as error:
                    record["attached_fixed_step_refusal"] = str(error)
                    feed.solve(dt, event=True)
                    break
            record["feed_rows"] = feed.rows
            if feed.state.attached_hexane_mass_kg != 0.0:
                raise RuntimeError("horizon reached before attached-hexane depletion")
            radial, activation_ledger = core.sphere.activate_radially(feed.state, feed.params)
            activation = {
                "status": "ACTUAL_FEED_TO_RADIAL_ACTIVATION",
                "rows": feed.rows,
                "interface_temperature_k": feed.interface_temperature,
                "activation": {
                    "radial": asdict(radial),
                    "ledger": asdict(activation_ledger),
                    "external_water_mass_kg_carried_separately": 0.0,
                    "external_water_energy_j_carried_separately": 0.0,
                },
            }
            record["activation"] = activation
            root_seed = None
            if args.root_seed:
                previous = json.loads(args.root_seed.read_text(encoding="utf-8"))
                if (
                    previous["activation"] != json.loads(json.dumps(activation))
                    or previous["inputs"]["birth_dt"] != args.birth_dt
                ):
                    raise ValueError(
                        "restart seed is not from this same activation and time interval"
                    )
                if not previous["sources_unchanged_during_run"]:
                    raise ValueError("restart source changed during its run")
                root_seed = {"study_debug": previous["debug"]}
                record["same_interval_seed_sha256"] = hashlib.sha256(
                    args.root_seed.read_bytes()
                ).hexdigest()
            before, config, controls, carry = shared.initialize_from_feed(
                activation, coefficient_case=args.coefficient_case
            )
            record["coefficient_case"] = {
                "identifier": args.coefficient_case,
                "pore_binary_diffusivity_m2_s": config.dry.binary_diffusivity.value_m2_s,
                "retained_water_mobility_mol_m_s": config.dry.retained_water_mobility.value_mol_m_s,
                "trajectory_parameter_fitted": False,
            }
            controls, config, controls_audit = active_binary_controls(controls, config)
            record["active_binary_chart"] = controls_audit
            if args.continued_water:
                config = continued.continued_configuration(config)
                before = replace(before, model=config.wet)
                controls = continued.continued_controls(controls, config)
                record["continued_water_law"] = type(config.wet.luikov).__name__
            original_coordinate_limit = controls.maximum_logit_magnitude
            controls = wider_numerical_controls(
                controls, args.coordinate_limit, args.regional_tolerance
            )
            record["native_coordinate_chart"] = {
                "original_magnitude": original_coordinate_limit,
                "declared_magnitude": args.coordinate_limit,
                "phase_and_conservation_gates_changed": False,
                "reason": "native gas-only barrier can exceed 30 in admissible low-moisture newborn states",
            }
            seed_options = {}
            if args.leading_seed:
                leading = leading_birth_seed(before, config, controls, feed.boundary)
                record["leading_birth_seed"] = leading
                if not leading["success"]:
                    raise RuntimeError("native leading-balance seed root did not close")
                seed_options.update(
                    interface_seed_temperature=leading["interface_temperature_k"],
                    film_beta_seed=leading["beta"],
                )
                profile = leading_pore_profile(before, config, feed.boundary, leading)
                record["leading_pore_profile"] = profile
                if not profile["success"]:
                    raise RuntimeError("native leading pore-profile shooting did not close")
                seed_options["leading_profile_seed"] = {**leading, **profile}
            checkpoint = {
                "before": before,
                "config": config,
                "controls": controls,
                "carry": carry,
                "boundary": feed.boundary,
                "feed": activation,
            }
            with args.out.with_suffix(".activation.pkl").open("wb") as stream:
                pickle.dump(checkpoint, stream)
            journal = []
            with shared.path_bounds.monotone_hexane_path_enclosure(journal):
                step, audit = shared.advance_birth(
                    before,
                    config,
                    controls,
                    carry,
                    feed.boundary,
                    args.birth_dt,
                    args.budget,
                    solver_policy="physical",
                    gas_only_surface=True,
                    root_seed=root_seed,
                    **seed_options,
                )
            record["attempts"].append(
                {"accepted": True, "stage": "birth", "audit": audit, "path_enclosures": journal}
            )
            checkpoint.update(
                before=step.after, carry={**carry, **audit["whole_particle_and_external_water"]}
            )
            with args.out.with_suffix(".accepted.pkl").open("wb") as stream:
                pickle.dump(checkpoint, stream)
            record["last_accepted"] = core.receding_summary(step.after, feed.dry)
            if args.steps:
                continue_radial(
                    checkpoint,
                    step.assembly.dry_face_fluxes[-1].surface_film_audit,
                    (1 - step.after.transport.geometry.front.z) / args.birth_dt,
                    record,
                    args,
                )
    except Exception as error:
        record["error"] = str(error)
        record["error_type"] = type(error).__name__
        record["debug"] = getattr(error, "study_debug", None)
        record["traceback"] = traceback.format_exc()
    record["wall_seconds"] = time.perf_counter() - started
    record["source_sha256_at_exit"] = hashes()
    record["sources_unchanged_during_run"] = entry == record["source_sha256_at_exit"]
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--water", type=float, default=0.05)
    parser.add_argument("--oil", type=float, default=0.0195)
    parser.add_argument("--radius", type=float, default=0.0008078927451503978)
    parser.add_argument("--cells", type=int, default=2)
    parser.add_argument("--initial-hexane", type=float, default=0.475)
    parser.add_argument("--initial-temperature", type=float)
    parser.add_argument("--bath-temperature", type=float, default=393.15)
    parser.add_argument(
        "--film-multiplier", type=float, default=INHERITED_FANER_COMMON_EXPOSURE_MULTIPLIER
    )
    parser.add_argument("--forcing-timescale", type=float, default=29.74922457229942)
    parser.add_argument("--feed-dt", type=float, default=0.25)
    parser.add_argument("--birth-dt", type=float, default=0.001)
    parser.add_argument("--budget", type=int, default=2400)
    parser.add_argument("--root-seed", type=Path)
    parser.add_argument("--leading-seed", action="store_true")
    parser.add_argument("--coordinate-limit", type=float, default=70.0)
    parser.add_argument(
        "--coefficient-case",
        choices=[case.identifier for case in core.p1ef.engineering_coefficient_cases()],
        default="nominal_log_center",
    )
    parser.add_argument("--horizon", type=float, default=240.0)
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--steps", type=int, default=0)
    parser.add_argument("--continuation-dt", type=float, default=0.01)
    parser.add_argument("--minimum-dt", type=float, default=1e-5)
    parser.add_argument("--wall-budget", type=float, default=120.0)
    parser.add_argument("--difference-mode", choices=["forward", "central"], default="forward")
    parser.add_argument("--regional-tolerance", type=float, default=1e-10)
    parser.add_argument("--continued-water", action="store_true")
    parser.add_argument("--extended-regional-root", action="store_true")
    args = parser.parse_args()
    if args.extended_regional_root and not args.resume:
        raise SystemExit("extended root study requires an actual accepted continuation checkpoint")
    if args.out.exists() or args.out.with_suffix(".activation.pkl").exists():
        raise SystemExit("refusing to overwrite existing evidence")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    record = run(args)
    args.out.write_text(json.dumps(record, indent=2, allow_nan=True), encoding="utf-8")
    print(
        json.dumps(
            {
                "out": str(args.out),
                "accepted_birth": any(
                    a.get("accepted") and a.get("stage") == "birth" for a in record["attempts"]
                ),
                "accepted_steps": sum(bool(a.get("accepted")) for a in record["attempts"]),
                "last_accepted": record.get("last_accepted"),
                "error": record.get("error"),
                "wall_seconds": record["wall_seconds"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
