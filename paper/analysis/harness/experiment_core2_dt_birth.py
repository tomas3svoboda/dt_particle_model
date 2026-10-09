"""Bounded native zero-old-volume DT birth from an actual feed activation.

Adds the existing transparent surface-water node's three conservative rows
to Core2's native birth residual. The incoming water energy remains the
native liquid INTERNAL energy at the feed event's body temperature; the
new node stores liquid internal energy at its solved surface temperature.
This is a candidate boundary-integration study, not universal qualification.
"""

from contextlib import contextmanager
from dataclasses import asdict, replace
import argparse
import hashlib
import json
import math
import pickle
from pathlib import Path
import time
import traceback

import numpy as np
from scipy.optimize import brentq, least_squares

import core2_water_film_node as node
import finite_film_core2_study as study
import native_water_path_bounds as path_bounds
from dtdc_simulator.core2.particle import cut_birth_event as birth
from dtdc_simulator.core2.particle import cut_birth_integrator as bi
from dtdc_simulator.core2.particle import wet_water as ww
from era_stamp import REPO_ROOT, era_stamp

core = study.core


class BirthAssemblyBudget(RuntimeError):
    """Includes seed, finite-difference, invalid trial and final assemblies."""


class UnresolvedTrialPath(ValueError):
    """Reject a trial whose enclosure work exceeds the per-path allowance."""


@contextmanager
def bounded_path_proof(point_budget, wall_s):
    """Stop unresolved native path proofs without admitting any unchecked path."""
    authority = core.cut.acfp._PathActivities
    original = authority._evaluate_coordinates
    started = time.perf_counter()
    audit = {"point_budget": point_budget, "wall_budget_s": wall_s, "points_started": 0}

    def evaluate(self, temperature_k, y_hexane):
        if audit["points_started"] >= point_budget or time.perf_counter() - started >= wall_s:
            audit["computational_budget_exhausted"] = True
            raise BirthAssemblyBudget(
                "native whole-path proof computational budget exhausted; no path admitted"
            )
        audit["points_started"] += 1
        path_work = getattr(self, "_study_trial_point_count", 0) + 1
        self._study_trial_point_count = path_work
        if path_work > 1200:
            audit["unresolved_trial_paths"] = audit.get("unresolved_trial_paths", 0) + 1
            raise UnresolvedTrialPath("trial path enclosure unresolved within 1200 evaluations")
        audit["last_path"] = {
            "left_temperature_k": self.left_temperature_k,
            "left_y_hexane": self.left_y_hexane,
            "delta_temperature_k": self.delta_temperature_k,
            "delta_y_hexane": self.delta_y_hexane,
        }
        audit["last_temperature_k"] = temperature_k
        audit["last_y_hexane"] = y_hexane
        return original(self, temperature_k, y_hexane)

    authority._evaluate_coordinates = evaluate
    try:
        yield audit
    finally:
        authority._evaluate_coordinates = original
        audit["wall_s"] = time.perf_counter() - started
        audit["native_callable_restored"] = authority._evaluate_coordinates is original


@contextmanager
def liquid_internal_energy_node(holder, dt, old_loading, old_temperature, pressure):
    """Reuse all node fluxes, replacing only its stored h by the chosen u.

    The caller already carries incoming u at old_temperature. This scoped
    substitution removes a caloric handoff change; it adds no latent source.
    """
    original = node.node_balances

    def balances(
        step_dt,
        liquid,
        old_liquid,
        unused_liquid_h,
        unused_old_h,
        pore_water,
        pore_hexane,
        pore_energy,
        film_water,
        film_hexane,
        film_energy,
    ):
        temperature = float(holder["x"][0])
        liquid_u = core.wa.state_Tp(temperature, pressure, "liquid").u_mass
        old_u = core.wa.state_Tp(old_temperature, pressure, "liquid").u_mass
        return original(
            step_dt,
            liquid,
            old_liquid,
            liquid_u,
            old_u,
            pore_water,
            pore_hexane,
            pore_energy,
            film_water,
            film_hexane,
            film_energy,
        )

    node.node_balances = balances
    try:
        with node.node_adapter(holder, dt, old_loading, old_temperature):
            yield
    finally:
        node.node_balances = original


def initialize_from_feed(
    feed_result, *, cells=None, coefficient_case="nominal_log_center", continued_water=False
):
    if feed_result["status"] != "ACTUAL_FEED_TO_RADIAL_ACTIVATION":
        raise ValueError("birth requires an actual accepted feed activation")
    event = feed_result["activation"]
    radial = event["radial"]
    radius = radial["grid"]["R"]
    count = len(radial["temperatures_k"]) if cells is None else cells
    temperature = radial["temperature_k"]
    if radial["coordinate"]["z"] != 1.0 or radial["attached_hexane_mass_kg"] != 0.0:
        raise ValueError("incoming radial endpoint must be exactly fully wet, zero attached H")
    for name in ("temperatures_k", "retained_water_loadings", "oil_fraction_labels"):
        if len(set(radial[name])) != 1:
            raise ValueError("this feed-event importer requires uniform native activation")
    if event["ledger"]["max_normalized_residual"] > 1e-10:
        raise ValueError("incoming activation did not conserve its material/energy")
    installed = core.q.install_engineering_case(
        coefficient_case,
        stress_multiplier=1.0,
        structural_conductivity_w_m_k=core.CONDUCTIVITY,
    )
    oil = radial["oil_fraction_labels"][0]
    water = radial["retained_water_loadings"][0]
    wet = replace(
        installed.configuration.wet.wet,
        X_water=water,
        w_o=oil,
        pressure_pa=core.PRESSURE_PA,
    )
    wet_model = replace(installed.configuration.wet, wet=wet)
    config = replace(
        installed.configuration.cut,
        wet=wet_model,
        dry=replace(
            installed.configuration.dry,
            pore=replace(installed.configuration.dry.pore, w_o=oil),
        ),
        moving_interface_composition_force_authority=(
            core.cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH
        ),
    )
    if continued_water:
        import positive_water_continuation as continuation

        config = continuation.continued_configuration(config)
        wet_model = config.wet
    grid = core.backend.uniform_grid(count, radius)
    # Uniform activation permits lossless grid refinement; these are the
    # actual event labels, never recalculated from a later temperature.
    historical = radial["wet_inventory"]["loadings"][0]
    before = ww.initialize(
        grid,
        (temperature,) * count,
        (water,) * count,
        (historical,) * count,
        wet_model,
        oil_fraction_labels=(oil,) * count,
    )
    before = replace(before, time_s=feed_result["rows"][-1]["t_s"])
    water_mass = event["external_water_mass_kg_carried_separately"]
    old_film_energy = event["external_water_energy_j_carried_separately"]
    predicted_film_energy = (
        water_mass * core.wa.state_Tp(temperature, core.PRESSURE_PA, "liquid").u_mass
    )
    input_energy_residual = math.fsum(
        (
            *before.reference_energy_cell_j,
            -radial["total_internal_energy_j"],
        )
    )
    scale = math.fsum(
        bi._wet_capacity_energy(v, c, config)
        for v, c in zip(grid.volumes, ww._evaluate_state_cells(before), strict=True)
    )
    if abs(input_energy_residual) / scale > 1e-10:
        raise ValueError("activation import changed native material energy")
    if abs(predicted_film_energy - old_film_energy) / scale > 1e-10:
        raise ValueError("incoming exterior-water energy does not match native u(Tbody)")
    carry = {
        "external_water_mass_kg": water_mass,
        "external_water_energy_j": old_film_energy,
        "external_water_temperature_k": temperature,
        "external_water_loading_kg_m2": water_mass / grid.areas[-1],
        "import_material_energy_residual_j": input_energy_residual,
        "import_film_energy_residual_j": predicted_film_energy - old_film_energy,
        "import_capacity_scale_j": scale,
        "incoming_film_beta": feed_result["rows"][-1]["beta"],
        "incoming_two_liquid_interface_temperature_k": feed_result["interface_temperature_k"],
    }
    original_controls = core.backend.CONTROLS
    # The historic 334.5 K solver bracket lies 0.0196 K above the native
    # binary two-liquid endpoint. Admit its already-defined cap branch in
    # a prospectively recorded numerical chart; phase checks remain native.
    lower = min(
        original_controls.interface_temperature_bounds_k[0],
        feed_result["interface_temperature_k"] + 1e-6,
    )
    controls = replace(
        original_controls,
        interface_temperature_bounds_k=(lower, original_controls.interface_temperature_bounds_k[1]),
    )
    carry["original_interface_solver_chart_k"] = original_controls.interface_temperature_bounds_k
    carry["birth_interface_solver_chart_k"] = controls.interface_temperature_bounds_k
    carry["interface_chart_amendment"] = (
        "numerical lower-bracket expansion to native two-liquid temperature +1e-6 K; "
        "phase, cap, balances and conservation unchanged"
    )
    if continued_water:
        controls = continuation.continued_controls(controls, config)
    return before, config, controls, carry


def interface_water_seed(before, config, controls):
    """Root the native saturated-H trace at the actual wet water loading."""
    target = before.retained_water_loadings[-1]
    historical = before.historical_hexane_loadings[-1]
    oil = before.oil_fraction_labels[-1]

    def difference(temperature):
        interface = core.cut.evaluate_interface_state(temperature, config, historical, oil)
        return interface.wet.retained_water_loading - target

    lo, hi = controls.interface_temperature_bounds_k
    root = brentq(difference, lo, hi, xtol=1e-11)
    return root, difference(root)


def film_at_beta(boundary, config, radius, temperature, beta, surface_y_hexane=None):
    """Invert the unchanged native film in its bounded heat-beta coordinate.

    The Stefan component derivative lies between bulk and surface mass
    fractions. Positive component heat capacities therefore bound d(beta)/dj
    between min(cp)/h and max(cp)/h: this scalar inverse is unique. No flux
    or physical residual is clipped, and every trial retains |beta|<=3.
    """
    if abs(beta) > boundary.accepted_beta_limit:
        raise ValueError("trial beta coordinate left its declared study band")
    pore = config.dry.pore
    y = (
        core.cp.decode_gas_only_y(temperature, config.dry.pressure_pa, -32.0, pore)
        if surface_y_hexane is None
        else surface_y_hexane
    )
    reference = core.ct.prepare_reduced_film_surface_state(
        boundary.film_boundary,
        cell_temperature_k=temperature,
        pressure_pa=config.dry.pressure_pa,
        particle_radius_m=radius,
        outer_half_cell_distance_m=0.5 * radius,
        particle_thermal_conductivity_w_m_k=config.dry.thermal_conductivity.value_w_m_k,
        binary_gas_interaction_k_wh=pore.k_wh,
    )
    film = core.ct._reduced_film_coefficients_at_surface_temperature(reference, temperature)
    capacities = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=temperature,
        bulk_temperature_k=boundary.temperature_k,
        pressure_pa=config.dry.pressure_pa,
        y_hexane=y,
        k_wh=pore.k_wh,
    )

    def evaluate(total_mass):
        components = core.efr.component_fluxes_for_supplied_total_mass_flux(
            film,
            total_mass_flux_kg_m2_s=total_mass,
            interface_hexane_mass_fraction=core._hexane_mass_fraction(y),
            bulk_hexane_mass_fraction=core._hexane_mass_fraction(boundary.y_hexane),
        )
        ack = core.efr.ackermann_heat_transfer(
            film,
            fluxes=components,
            water_film_heat_capacity_j_kg_k=capacities.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=capacities.hexane_j_kg_k,
            enforce_production_domain=False,
        )
        return components, ack

    beta_zero = evaluate(0.0)[1].beta
    increment = (beta - beta_zero) * film.heat_transfer_w_m2_k
    cp_min = min(capacities.water_j_kg_k, capacities.hexane_j_kg_k)
    cp_max = max(capacities.water_j_kg_k, capacities.hexane_j_kg_k)
    bounds = sorted((increment / cp_min, increment / cp_max))
    if increment == 0.0 or bounds[0] == bounds[1]:
        mass = bounds[0]
    else:
        # Expand the mathematically valid bracket only to retain endpoint
        # signs under floating-point evaluation of the analytic partition.
        roundoff = 64 * math.ulp(max(abs(bounds[0]), abs(bounds[1]), 1.0))
        mass = brentq(
            lambda value: evaluate(value)[1].beta - beta,
            bounds[0] - roundoff,
            bounds[1] + roundoff,
            xtol=1e-14,
            rtol=1e-14,
            maxiter=80,
        )
    components, ack = evaluate(mass)
    if abs(ack.beta - beta) > 1e-10:
        raise RuntimeError("native beta-to-mass inverse did not conserve its coordinate")
    return mass, components, ack, film


def advance_birth(
    before,
    config,
    controls,
    carry,
    boundary,
    dt,
    budget,
    *,
    front_z=None,
    solver_policy="scipy",
    root_seed=None,
    gas_only_surface=False,
    interface_seed_temperature=None,
    film_beta_seed=None,
    leading_profile_seed=None,
):
    layout = birth.layout_for_birth(before)
    bi._validate_birth_chart(before, config, controls)
    chart = bi._coordinate_chart(before, config, controls)
    gamma_water_root, seed_water_error = interface_water_seed(before, config, controls)
    # The thin newborn shell's liquid-active interface tends to the native
    # two-liquid endpoint. This is a seed only, never an imposed T_Gamma.
    # The old wet bulk temperatures remain the actual event temperatures.
    gamma = (
        gamma_water_root
        if gas_only_surface
        else carry["incoming_two_liquid_interface_temperature_k"] + 1e-4
    )
    if interface_seed_temperature is not None:
        if not gas_only_surface:
            raise ValueError("declared interface seed override is limited to the gas-only study")
        gamma = float(interface_seed_temperature)
    if not chart.interface_temperature[0] < gamma < chart.interface_temperature[1]:
        raise ValueError("actual event thermal trace is outside native interface chart")
    pore = replace(config.dry.pore, w_o=before.oil_fraction_labels[-1])
    seed_beta = float(carry["incoming_film_beta"])
    if film_beta_seed is not None:
        if not gas_only_surface:
            raise ValueError("declared beta seed override is limited to the gas-only study")
        seed_beta = float(film_beta_seed)
    gamma_interval = (
        core.cp.gas_only_composition_interval(gamma + 1e-4, config.dry.pressure_pa, pore)
        if gas_only_surface
        else None
    )
    seed_surface_y = sum(gamma_interval.y_hexane_bounds) / 2 if gas_only_surface else None
    _, components, ack, film = film_at_beta(
        boundary,
        config,
        before.grid.R,
        gamma + (1e-4 if gas_only_surface else 0.0),
        seed_beta,
        surface_y_hexane=seed_surface_y,
    )
    gas_heat = film.heat_transfer_w_m2_k * ack.factor * (boundary.temperature_k - gamma)
    interface = core.cut.evaluate_interface_state(
        gamma, config, before.historical_hexane_loadings[-1], before.oil_fraction_labels[-1]
    )
    phase_energy = (
        interface.dry.hexane_gas_partial_enthalpy_j_mol / core.hx.M
        - core.hx.state_Tp(gamma, config.dry.pressure_pa, "liquid").u_mass
    )
    storage = config.wet.wet.rho_dm_p * before.historical_hexane_loadings[-1]
    wet_release = (
        config.wet.wet.conductivity
        * (before.temperatures_k[-1] - gamma)
        / (before.grid.R - before.grid.centers[-1])
    )
    thickness = (gas_heat + wet_release) * dt / (storage * phase_energy)
    if leading_profile_seed is not None:
        if not gas_only_surface:
            raise ValueError("limiting pore-profile seed is limited to the gas-only study")
        thickness = leading_profile_seed["speed_m_s"] * dt
    if not 0 < thickness < before.grid.R - before.grid.faces[-2]:
        raise ValueError("energy-based birth seed did not fit the outer material cell")
    if front_z is None:
        front_z = ((before.grid.R - thickness) / before.grid.R) ** 3
    else:
        thickness = before.grid.R * (1 - front_z ** (1 / 3))
    dry_temperature = gamma + gas_heat * thickness / (2 * core.CONDUCTIVITY)
    if leading_profile_seed is not None:
        dry_temperature = gamma + leading_profile_seed["temperature_gradient_k_m"] * thickness / 2
    if not chart.dry_temperature[0] < dry_temperature < chart.dry_temperature[1]:
        raise ValueError("native water-trace seed does not admit a dry thermal state")
    interval = core.cp.gas_only_composition_interval(dry_temperature, config.dry.pressure_pa, pore)
    dry_y = 0.5 * (interval.lower_y_hexane + interval.upper_y_hexane)
    if leading_profile_seed is not None:
        dry_y = (
            interface.y_hexane + leading_profile_seed["composition_gradient_m_inv"] * thickness / 2
        )
    stefan = components.water_kg_m2_s / core.wa.M + components.hexane_kg_m2_s / core.hx.M
    if leading_profile_seed is not None:
        stefan = leading_profile_seed["total_stefan_flux_mol_m2_s"]
    seed = bi.BirthStepSeed(
        birth.BirthUnknowns(
            before.temperatures_k,
            before.retained_water_loadings,
            dry_temperature,
            dry_y,
            (stefan, stefan),
            front_z,
            gamma,
        ),
        (1.0, 1.0),
        "native two-liquid endpoint seed; actual old wet state; energy-scaled newborn thickness",
    )
    bi._validate_seed(layout, seed)
    native_x = bi._encode_candidate(seed.candidate, chart, seed.stefan_flux_scales_mol_m2_s)
    lower, upper = core.ci._coordinate_solver_bounds(layout, controls.maximum_logit_magnitude)
    old_loading = carry["external_water_loading_kg_m2"]
    if gas_only_surface and old_loading != 0.0:
        raise ValueError("gas-only birth requires exactly absent incoming external water")
    surface_temperature = gamma + gas_heat * thickness / core.CONDUCTIVITY
    if leading_profile_seed is not None:
        surface_temperature = gamma + leading_profile_seed["temperature_gradient_k_m"] * thickness
    surface_interval = (
        core.cp.gas_only_composition_interval(surface_temperature, config.dry.pressure_pa, pore)
        if gas_only_surface
        else None
    )
    seed_surface_y = sum(surface_interval.y_hexane_bounds) / 2 if gas_only_surface else None
    if leading_profile_seed is not None:
        seed_surface_y = (
            interface.y_hexane + leading_profile_seed["composition_gradient_m_inv"] * thickness
        )
    initial_mass, components, _, _ = film_at_beta(
        boundary,
        config,
        before.grid.R,
        surface_temperature,
        seed_beta,
        surface_y_hexane=seed_surface_y,
    )
    initial_liquid = 0.0 if gas_only_surface else old_loading - dt * components.water_kg_m2_s
    if not gas_only_surface and initial_liquid <= 0:
        raise ValueError("birth guess exhausted the carried exterior water")
    holder = {"x": np.array([surface_temperature, initial_mass, initial_liquid])}
    if gas_only_surface:
        holder["y_hexane"] = seed_surface_y
    solver_surface_seed = np.array(
        [
            surface_temperature,
            seed_beta,
            (
                (
                    (seed_surface_y - surface_interval.lower_y_hexane)
                    / (surface_interval.upper_y_hexane - surface_interval.lower_y_hexane)
                )
                if gas_only_surface
                else initial_liquid
            ),
        ]
    )
    if root_seed is not None:
        candidate_data = root_seed["study_debug"]["candidate"]
        previous = birth.BirthUnknowns(
            **{
                key: tuple(value) if isinstance(value, list) else value
                for key, value in candidate_data.items()
                if key != "physically_qualifying"
            }
        )
        seed = bi.BirthStepSeed(
            previous,
            seed.stefan_flux_scales_mol_m2_s,
            "unaccepted same-IVP nonlinear iterate; no trajectory state imported",
        )
        native_x = bi._encode_candidate(seed.candidate, chart, seed.stefan_flux_scales_mol_m2_s)
        surface_data = root_seed["study_debug"]["surface"]
        solver_surface_seed = np.array(
            [
                surface_data["temperature_k"],
                surface_data["beta"],
                (
                    (
                        surface_data["surface_y_hexane"]
                        - core.cp.gas_only_composition_interval(
                            surface_data["temperature_k"], config.dry.pressure_pa, pore
                        ).lower_y_hexane
                    )
                    / (
                        core.cp.gas_only_composition_interval(
                            surface_data["temperature_k"], config.dry.pressure_pa, pore
                        ).upper_y_hexane
                        - core.cp.gas_only_composition_interval(
                            surface_data["temperature_k"], config.dry.pressure_pa, pore
                        ).lower_y_hexane
                    )
                    if gas_only_surface
                    else surface_data["liquid_water_kg_m2"]
                ),
            ]
        )
        if gas_only_surface:
            holder["y_hexane"] = surface_data["surface_y_hexane"]
        seed_mass = film_at_beta(
            boundary,
            config,
            before.grid.R,
            solver_surface_seed[0],
            solver_surface_seed[1],
            surface_y_hexane=holder.get("y_hexane"),
        )[0]
        holder["x"] = np.array(
            [solver_surface_seed[0], seed_mass, 0.0 if gas_only_surface else solver_surface_seed[2]]
        )
        dry_temperature = previous.dry_temperature_k
    counter = {"attempts": 0, "rejected": 0}
    last = {
        "encoded_initial_candidate": native_x.tolist(),
        "native_lower_bounds": lower.tolist(),
        "native_upper_bounds": upper.tolist(),
        "initial_surface_coordinates": solver_surface_seed.tolist(),
    }

    def assemble(candidate):
        if counter["attempts"] >= budget:
            raise BirthAssemblyBudget(f"total birth residual-assembly budget {budget} exhausted")
        counter["attempts"] += 1
        try:
            return birth.assemble_birth_backward_euler(
                before,
                config,
                candidate,
                dt,
                boundary,
                enforce_reduced_film_thresholds=False,
            )
        except Exception:
            counter["rejected"] += 1
            raise

    try:
        with study.magnitude_storage_scale() as capacity_journal:
            with study.installed_boundary_adapter(numerical_cap_seam=True):
                with liquid_internal_energy_node(
                    holder,
                    dt,
                    old_loading,
                    carry["external_water_temperature_k"],
                    config.dry.pressure_pa,
                ):
                    initial = assemble(seed.candidate)
                    scales = bi._residual_scales(initial)
                    scale_vector = np.asarray(scales.vector)
                    reference = core.ct.prepare_reduced_film_surface_state(
                        boundary.film_boundary,
                        cell_temperature_k=dry_temperature,
                        pressure_pa=config.dry.pressure_pa,
                        particle_radius_m=before.grid.R,
                        outer_half_cell_distance_m=(
                            before.grid.R
                            - core.cut._dry_piece_center(
                                initial.candidate_geometry.cells[-1],
                                initial.candidate_geometry.front.radius_m,
                            )
                        ),
                        particle_thermal_conductivity_w_m_k=(
                            config.dry.thermal_conductivity.value_w_m_k
                        ),
                        binary_gas_interaction_k_wh=pore.k_wh,
                    )
                    mass_scale = max(
                        reference.coefficients.film_density_kg_m3
                        * reference.coefficients.binary_mass_transfer_m_s,
                        old_loading / dt,
                        1e-8,
                    )
                    heat_scale = max(
                        reference.coefficients.heat_transfer_w_m2_k
                        * abs(boundary.temperature_k - dry_temperature),
                        old_loading
                        * abs(
                            core.wa.state_Tp(
                                dry_temperature, config.dry.pressure_pa, "liquid"
                            ).u_mass
                        )
                        / dt,
                        1.0,
                    )
                    row_weights = np.r_[
                        np.full(
                            layout.unknown_count,
                            min(
                                1.0,
                                controls.nonlinear_residual_tolerance
                                / bi.BIRTH_REGIONAL_RH_ACCEPTANCE,
                            ),
                        ),
                        np.ones(3),
                    ]

                    def residual(x):
                        if gas_only_surface:
                            if not 0.0 < x[-1] < 1.0:
                                raise ValueError(
                                    "birth surface left its native open gas-only interval"
                                )
                            lo, hi = core.cp.gas_only_composition_interval(
                                float(x[-3]), config.dry.pressure_pa, pore
                            ).y_hexane_bounds
                            holder["y_hexane"] = lo + (hi - lo) * float(x[-1])
                        total_mass, _, _, _ = film_at_beta(
                            boundary,
                            config,
                            before.grid.R,
                            float(x[-3]),
                            float(x[-2]),
                            surface_y_hexane=holder.get("y_hexane"),
                        )
                        holder["x"] = np.array(
                            [x[-3], total_mass, 0.0 if gas_only_surface else x[-1]]
                        )
                        candidate = bi._decode_candidate(
                            x[:-3], layout, chart, seed.stefan_flux_scales_mol_m2_s
                        )
                        assembly = assemble(candidate)
                        audit = assembly.dry_face_fluxes[-1].surface_film_audit
                        native = (
                            np.asarray(birth.datum_covariant_residual_vector(assembly.residuals))
                            / scale_vector
                        )
                        values = np.r_[
                            native,
                            audit.water_node_residual_kg_m2_s / mass_scale,
                            audit.hexane_node_residual_kg_m2_s / mass_scale,
                            audit.energy_node_residual_w_m2 / heat_scale,
                        ]
                        last.update(
                            candidate=asdict(candidate),
                            surface=asdict(audit),
                            scaled_residuals=values.tolist(),
                            maximum_scaled_residual=float(np.max(np.abs(values))),
                            external_beta_coordinate=float(x[-2]),
                        )
                        return values * row_weights

                    temperature_lo, temperature_hi = chart.dry_temperature
                    film_loading_hi = 1.0 if gas_only_surface else max(0.05, 4.0 * old_loading)

                    def physical_solver(fun, x0, **kwargs):
                        def interval(temperature):
                            return core.cp.gas_only_composition_interval(
                                temperature, config.dry.pressure_pa, chart.newborn_surface_pore
                            ).y_hexane_bounds

                        def decode(values):
                            candidate = bi._decode_candidate(
                                values[:-3], layout, chart, seed.stefan_flux_scales_mol_m2_s
                            )
                            lo, hi = interval(candidate.dry_temperature_k)
                            candidate = replace(
                                candidate, dry_y_hexane=(candidate.dry_y_hexane - lo) / (hi - lo)
                            )
                            return np.r_[candidate.vector(), values[-3:]]

                        def encode(values):
                            candidate = birth.BirthUnknowns.from_vector(layout, values[:-3])
                            if not 0 < candidate.dry_y_hexane < 1:
                                raise ValueError("gas interval fraction outside its open domain")
                            lo, hi = interval(candidate.dry_temperature_k)
                            candidate = replace(
                                candidate, dry_y_hexane=lo + (hi - lo) * candidate.dry_y_hexane
                            )
                            encoded = bi._encode_candidate(
                                candidate, chart, seed.stefan_flux_scales_mol_m2_s
                            )
                            return np.r_[encoded, values[-3:]]

                        n = layout.wet_piece_count
                        pairs = (
                            [chart.wet_temperature] * n
                            + [chart.wet_water] * n
                            + [
                                chart.dry_temperature,
                                (0.0, 1.0),
                                (-math.inf, math.inf),
                                (-math.inf, math.inf),
                                chart.front_z,
                                chart.interface_temperature,
                                (temperature_lo, min(temperature_hi, boundary.temperature_k)),
                                (-boundary.accepted_beta_limit, boundary.accepted_beta_limit),
                                (0.0, film_loading_hi),
                            ]
                        )
                        return study.numerics.physical_chart_solve(
                            fun,
                            x0,
                            decode=decode,
                            encode=encode,
                            physical_bounds=tuple(zip(*pairs)),
                            temperature_indices=[*range(n), 2 * n, 2 * n + 5, 2 * n + 6],
                            **{**kwargs, "gtol": controls.nonlinear_residual_tolerance},
                        )

                    def solve(solver, fun, x0, **kwargs):
                        if solver_policy == "physical":
                            kwargs.pop("temperature_index")
                            return physical_solver(fun, x0, **kwargs)
                        return study.numerics.centered_temperature_solve(solver, fun, x0, **kwargs)

                    solution = solve(
                        (
                            physical_solver
                            if solver_policy == "physical"
                            else (
                                study.numerics.domain_aware_newton
                                if solver_policy == "newton"
                                else least_squares
                            )
                        ),
                        residual,
                        np.r_[native_x, solver_surface_seed],
                        temperature_index=-3,
                        bounds=(
                            np.r_[lower, temperature_lo, -boundary.accepted_beta_limit, 0.0],
                            np.r_[
                                upper,
                                min(temperature_hi, boundary.temperature_k),
                                boundary.accepted_beta_limit,
                                film_loading_hi,
                            ],
                        ),
                        method="trf",
                        jac="2-point",
                        x_scale="jac",
                        diff_step=None,
                        max_nfev=budget,
                        ftol=controls.nonlinear_step_tolerance,
                        xtol=controls.nonlinear_step_tolerance,
                        gtol=controls.nonlinear_step_tolerance,
                    )
                    values = residual(solution.x) / row_weights
                    last["nonlinear_solver"] = {
                        "success": bool(solution.success),
                        "message": str(solution.message),
                        "residual_calls": int(solution.nfev),
                        "active_mask": solution.active_mask.tolist(),
                        "physical_coordinates": getattr(
                            solution, "physical_coordinates", np.array([])
                        ).tolist(),
                    }
                    candidate = bi._decode_candidate(
                        solution.x[:-3], layout, chart, seed.stefan_flux_scales_mol_m2_s
                    )
                    assembly = assemble(candidate)
                    audit = assembly.dry_face_fluxes[-1].surface_film_audit
                    topology = bi._topology_diagnostics(assembly, chart)
                    bi._bounded_fractions(candidate, chart, topology)
                    condition = core.ci._jacobian_condition_proxy(
                        getattr(solution, "native_jac", solution.jac)
                        if gas_only_surface
                        else solution.jac
                    )
                    native_acceptance = max(
                        controls.nonlinear_residual_tolerance, bi.BIRTH_REGIONAL_RH_ACCEPTANCE
                    )
                    if not solution.success or np.any(solution.active_mask != 0):
                        raise RuntimeError(
                            "joint birth root did not certify success in open bounds"
                        )
                    if np.max(np.abs(values[:-3])) > native_acceptance:
                        raise RuntimeError(
                            "native birth regional/RH residual exceeds unchanged contract"
                        )
                    if np.max(np.abs(values[-3:])) > controls.nonlinear_residual_tolerance:
                        raise RuntimeError("external-water node residual exceeds 1e-10 contract")
                    if not math.isfinite(condition) or condition > controls.maximum_condition_proxy:
                        raise RuntimeError(f"joint Jacobian condition not certified: {condition}")
                    if abs(audit.beta) > boundary.accepted_beta_limit:
                        raise RuntimeError(
                            "surface root exceeds declared signed Ackermann study band"
                        )
                    if (
                        not gas_only_surface
                        and not 0.0 <= audit.water_activity_trace_deficit <= 2e-6
                    ):
                        raise RuntimeError(
                            "surface free-water activity approximation exceeds declared audit"
                        )
                    step = bi._accept_birth(
                        before,
                        controls,
                        boundary,
                        seed,
                        assembly,
                        scales,
                        tuple(values[:-3]),
                        counter["attempts"],
                        counter["rejected"],
                        str(solution.message),
                        condition,
                        chart,
                        topology,
                        None,
                    )
                    ledger = step.ledger
                    if (
                        max(
                            ledger.maximum_step_ledger_residual,
                            ledger.maximum_cumulative_ledger_residual,
                            ledger.birth_telescoping_residual,
                        )
                        > controls.ledger_tolerance
                    ):
                        raise RuntimeError(
                            "native birth conservation/recession ledger exceeded contract"
                        )
                    if ledger.relative_immediate_recession_error > controls.ledger_tolerance:
                        raise RuntimeError(
                            "birth kinematic recession error exceeds native tolerance"
                        )
                    if assembly.ledger.old_dry_volume_m3 != 0.0:
                        raise RuntimeError("birth imported nonzero old dry volume")
                    if assembly.ledger.minimum_entropy_production_w_m3_k < 0:
                        raise RuntimeError("native birth has negative entropy production")
                    for water, dual in (
                        (before.retained_water_loadings, (0.0,) * before.grid.n),
                        (
                            step.after.transport.wet_retained_water_loadings,
                            step.after.transport.effective_wet_retained_water_capacity_duals_over_rt,
                        ),
                    ):
                        core.ci.certify_explicit_wet_retained_water_applicability(
                            controls, water, dual, config.wet.luikov
                        )
                    combined = combined_ledger(before, config, carry, step, audit)
                    if combined["maximum_normalized_residual"] > controls.ledger_tolerance:
                        raise RuntimeError(
                            "combined body/surface film conservation exceeded contract"
                        )
                    return step, {
                        "accepted": True,
                        "gas_only_surface": gas_only_surface,
                        "total_residual_assembly_attempts": counter["attempts"],
                        "seed_interface_water_root_temperature_k": gamma_water_root,
                        "seed_interface_water_root_loading_residual": seed_water_error,
                        "seed_interface_temperature_k": gamma,
                        "seed_newborn_thickness_m": thickness,
                        "surface_flux_coordinate": "exact native heat-beta inverse bounded to declared study band",
                        "native_birth_ledger": asdict(ledger),
                        "surface_water_node": asdict(audit),
                        "whole_particle_and_external_water": combined,
                        "signed_numerical_capacity_samples": capacity_journal,
                        "after_transport": asdict(step.after.transport),
                        "after_inventory": asdict(core.ci.inventory_snapshot(step.after.transport)),
                    }
    except Exception as exc:
        failure = bi.BirthIntegratorStepError(
            f"actual-feed joint birth refused with exact fully-wet rollback: {exc}",
            before,
            nonlinear_evaluations=counter["attempts"],
            rejected_trial_evaluations=counter["rejected"],
            maximum_scaled_residual=last.get("maximum_scaled_residual", math.inf),
        )
        failure.study_debug = {
            **last,
            "initial_seed": asdict(seed),
            "seed_surface_temperature_k": surface_temperature,
            "seed_newborn_thickness_m": thickness,
            "seed_gas_heat_w_m2": gas_heat,
            "seed_wet_release_w_m2": wet_release,
            "seed_native_water_trace_root_temperature_k": gamma_water_root,
            "surface_flux_coordinate": "exact native heat-beta inverse bounded to declared study band",
            "optimizer_debug": getattr(exc, "optimizer_debug", None),
        }
        failure.total_residual_assembly_attempts = counter["attempts"]
        raise failure from exc


def combined_ledger(before, config, carry, step, audit):
    if isinstance(before, ww.WetWaterState):
        area = before.grid.areas[-1]
        old = bi._fully_wet_current_inventory(before, config)
    elif hasattr(before, "transport"):
        area = before.transport.geometry.master_grid.areas[-1]
        old = core.ci.inventory_snapshot(before.transport)
    else:
        area = before.geometry.master_grid.areas[-1]
        old = before.current_inventory
    if hasattr(step.after, "transport"):
        dt = step.ledger.dt_s
        new = core.ci.inventory_snapshot(step.after.transport)
        new_capacity = core.ci.capacity_energy_scale(step.after.transport)
    elif hasattr(step.after, "current_inventory"):
        from dtdc_simulator.core2.particle import cut_face_event_integrator as fi

        dt = step.ledger.event_time_s
        new = step.after.current_inventory
        new_capacity = fi._event_capacity_scale(step.assembly)
    else:
        dt = step.ledger.event_duration_s
        new = step.assembly.current_inventory
        new_capacity = math.fsum(
            a.capacity_j_m3_k * v * max(abs(t), 1.0)
            for a, v, t in zip(
                step.ledger.residual_scales.dry_energy_capacity_stencil_audits,
                step.after.grid.volumes,
                step.after.temperatures_k,
            )
        )
    water_mass = area * audit.liquid_water_kg_m2
    water_energy = (
        water_mass * core.wa.state_Tp(audit.temperature_k, config.dry.pressure_pa, "liquid").u_mass
    )
    water_out = dt * area * audit.external_water_flux_kg_m2_s
    hexane_out = dt * area * audit.external_hexane_flux_kg_m2_s
    energy_out = dt * area * audit.external_energy_flux_w_m2
    water_residual = math.fsum(
        (
            (new.total_water_mol - old.total_water_mol) * core.wa.M,
            water_mass,
            -carry["external_water_mass_kg"],
            water_out,
        )
    )
    hexane_residual = math.fsum(
        ((new.total_hexane_mol - old.total_hexane_mol) * core.hx.M, hexane_out)
    )
    energy_residual = math.fsum(
        (
            new.total_energy_j,
            -old.total_energy_j,
            water_energy,
            -carry["external_water_energy_j"],
            energy_out,
        )
    )
    normalized = (
        abs(water_residual)
        / max(old.total_water_mol * core.wa.M + carry["external_water_mass_kg"], abs(water_out)),
        abs(hexane_residual) / max(old.total_hexane_mol * core.hx.M, abs(hexane_out)),
        abs(energy_residual) / max(carry["import_capacity_scale_j"], new_capacity),
    )
    return {
        "water_residual_kg": water_residual,
        "hexane_residual_kg": hexane_residual,
        "energy_residual_j": energy_residual,
        "normalized_residuals": normalized,
        "maximum_normalized_residual": max(normalized),
        "external_water_mass_kg": water_mass,
        "external_water_energy_j": water_energy,
        "external_water_temperature_k": audit.temperature_k,
        "old_film_energy_j_carried_without_reset": carry["external_water_energy_j"],
        "exterior_water_caloric": "liquid internal energy; old Tbody, solved new Tsurface",
    }


def run(args):
    args.feed = args.feed.resolve()
    args.bath = args.bath.resolve()
    sources = [
        Path(__file__).resolve(),
        args.feed,
        args.bath,
        Path(__file__).with_name("positive_water_continuation.py"),
        *(
            Path(module.__file__).resolve()
            for module in (
                node,
                study.numerics,
                path_bounds,
                study,
                birth,
                bi,
                ww,
                core.cut,
                core.ci,
                core.ct,
                core.cp,
                core.efr,
                core.hx,
                core.wa,
                core.bg,
                core.sp,
            )
        ),
    ]
    result = {
        "status": "BIRTH_INPUT_SETUP_REFUSED",
        "accepted": False,
        "complete_dt_destiny": False,
        "physically_qualifying": False,
        "universal_model_alignment_complete": False,
        "remaining_alignment": "whole-trajectory and manuscript equation audit remain required",
        "feed_result": str(args.feed),
        "dt_s": args.dt_s,
        "solver_policy": args.solver,
        "birth_only_numerical_acceptance": bi.BIRTH_REGIONAL_RH_ACCEPTANCE,
        "residual_assembly_budget": args.budget,
        "entry_source_sha256": {},
    }
    before = None
    started = time.perf_counter()
    try:
        for path in sources:
            result["entry_source_sha256"][str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        feed = json.loads(args.feed.read_text(encoding="utf-8"))
        bath = json.loads(args.bath.read_text(encoding="utf-8"))
        before, config, controls, carry = initialize_from_feed(
            feed,
            cells=args.cells,
            coefficient_case=args.coefficient_case,
            continued_water=args.continued_water,
        )
        result["coefficient_case"] = args.coefficient_case
        result["water_law"] = type(config.wet.luikov).__name__
        root_seed = None
        if args.root_seed is not None:
            root_seed = json.loads(args.root_seed.read_text(encoding="utf-8"))
            if root_seed["dt_s"] != args.dt_s or root_seed["before_time_s"] != before.time_s:
                raise ValueError("nonlinear seed does not belong to the same physical-time IVP")
            previous_hashes = root_seed["entry_source_sha256"]
            for source in (args.feed, args.bath, *[p for p in sources if "src" in p.parts]):
                if (
                    previous_hashes.get(str(source))
                    != hashlib.sha256(source.read_bytes()).hexdigest()
                ):
                    raise ValueError(
                        "nonlinear seed physical input/source differs from current IVP"
                    )
            result["same_ivp_unaccepted_root_seed"] = {
                "path": str(args.root_seed),
                "sha256": hashlib.sha256(args.root_seed.read_bytes()).hexdigest(),
                "trajectory_state_imported": False,
            }
        clock = before.time_s + args.dt_s
        temperature, y_hexane = core.interpolate_bath(bath["waypoints"], clock)
        native_boundary = core.reduced_film_boundary(
            temperature_k=temperature,
            y_hexane=y_hexane,
            pair=core.film_pair_at_radius(before.grid.R),
            forcing_timescale_s=max(args.dt_s, 0.075),
            label=f"actual feed birth on unchanged prescribed DT bath at t={clock}",
        )
        boundary = study.FiniteFilmStudyBoundary(
            temperature,
            y_hexane,
            "candidate conservative transparent exterior-water boundary",
            film_boundary=native_boundary,
            accepted_beta_limit=3.0,
            actual_surface_path=True,
            numerical_cap_seam=True,
        )
        result.update(
            feed_external_film_law=feed.get("external_film_law"),
            shared_external_film_law=str(feed.get("external_film_law", "")).startswith(
                "native EFr mass-coordinate"
            ),
            before_time_s=before.time_s,
            before_grid_cells=before.grid.n,
            before_old_dry_volume_m3=0.0,
            incoming_external_water=carry,
            boundary_temperature_k=temperature,
            boundary_y_hexane=y_hexane,
        )
        path_journal = []
        result["native_water_monotonicity_certificates"] = path_journal
        with (
            bounded_path_proof(args.path_point_budget, args.wall_budget_s) as proof_audit,
            path_bounds.monotone_water_path_enclosure(path_journal),
            path_bounds.monotone_hexane_path_enclosure(path_journal),
        ):
            result["birth_numerical_allowance_audit"] = {
                "native_regional_rh_acceptance": bi.BIRTH_REGIONAL_RH_ACCEPTANCE,
                "strict_native_certification": True,
                "native_constants_mutated": False,
            }
            result["native_path_proof_budget"] = proof_audit
            if args.resume is None:
                birth_step, record = advance_birth(
                    before,
                    config,
                    controls,
                    carry,
                    boundary,
                    args.dt_s,
                    args.budget,
                    front_z=args.front_z,
                    solver_policy=args.solver,
                    root_seed=root_seed,
                )
                result.update(record, status="ACTUAL_FEED_ZERO_OLD_DRY_VOLUME_BIRTH_ACCEPTED")
                result["radial_continuation"] = []
                current = birth_step.after
                surface = birth_step.assembly.dry_face_fluxes[-1].surface_film_audit
                combined = record["whole_particle_and_external_water"]
                carry = {
                    **carry,
                    **combined,
                    "external_water_loading_kg_m2": surface.liquid_water_kg_m2,
                }
            else:
                checkpoint = pickle.loads(args.resume.read_bytes())
                for source in (args.feed, args.bath, *[p for p in sources if "src" in p.parts]):
                    if (
                        checkpoint["source_sha256"].get(str(source))
                        != hashlib.sha256(source.read_bytes()).hexdigest()
                    ):
                        raise ValueError("accepted checkpoint physical input/source mismatch")
                current = checkpoint["current"]
                config = current.transport.config
                surface = checkpoint["surface"]
                carry = checkpoint["carry"]
                result.update(
                    status="DT_ACCEPTED_STATE_CONTINUATION",
                    accepted=True,
                    resumed_accepted_checkpoint={
                        "path": str(args.resume),
                        "time_s": current.transport.time_s,
                        "sha256": hashlib.sha256(args.resume.read_bytes()).hexdigest(),
                    },
                    radial_continuation=[],
                )
            for index in range(args.continuation_steps):
                previous = current
                try:
                    dt = args.continuation_dt
                    clock = current.transport.time_s + dt
                    temperature, y_hexane = core.interpolate_bath(bath["waypoints"], clock)
                    native_boundary = core.reduced_film_boundary(
                        temperature_k=temperature,
                        y_hexane=y_hexane,
                        pair=core.film_pair_at_radius(before.grid.R),
                        forcing_timescale_s=max(dt, 0.075),
                        label=f"same prescribed DT bath, continued at t={clock}",
                    )
                    next_boundary = study.FiniteFilmStudyBoundary(
                        temperature,
                        y_hexane,
                        "conservative transparent external-water continuation",
                        film_boundary=native_boundary,
                        accepted_beta_limit=3.0,
                        actual_surface_path=True,
                        numerical_cap_seam=True,
                    )
                    candidate = core.cut.candidate_from_state(
                        current.transport,
                        dry_total_stefan_fluxes_mol_m2_s=current.last_total_stefan_fluxes_mol_m2_s,
                        interface_temperature_k=current.last_interface_temperature_k,
                    )
                    candidate = replace(
                        candidate, front_z=candidate.front_z * (1 - min(0.001, 0.001 * dt / 0.01))
                    )
                    seed = core.backend._seed(
                        candidate, "same accepted DT particle; numerical continuation seed"
                    )
                    old_loading = carry["external_water_loading_kg_m2"]
                    old_temperature = carry["external_water_temperature_k"]

                    def energy_context(holder):
                        return liquid_internal_energy_node(
                            holder, dt, old_loading, old_temperature, config.dry.pressure_pa
                        )

                    def decoder(ts, beta):
                        return film_at_beta(next_boundary, config, before.grid.R, ts, beta)[0]

                    with study.magnitude_storage_scale() as scale_journal:
                        next_step = node.advance_water_node(
                            current,
                            dt,
                            next_boundary,
                            seed,
                            old_loading=old_loading,
                            old_temperature=old_temperature,
                            surface_seed=[
                                surface.temperature_k,
                                surface.external_water_flux_kg_m2_s
                                + surface.external_hexane_flux_kg_m2_s,
                                surface.liquid_water_kg_m2,
                            ],
                            physical_solver=True,
                            residual_budget=args.budget,
                            beta_decoder=decoder,
                            energy_context=energy_context,
                        )
                    next_surface = next_step.surface_film_audit
                    next_combined = combined_ledger(current, config, carry, next_step, next_surface)
                    if (
                        next_combined["maximum_normalized_residual"]
                        > current.controls.ledger_tolerance
                    ):
                        raise RuntimeError(
                            "combined particle/water-film continuation ledger refused"
                        )
                    result["radial_continuation"].append(
                        {
                            "index": index,
                            "accepted": True,
                            "time_s": next_step.after.transport.time_s,
                            "ledger": asdict(next_step.ledger),
                            "surface": asdict(next_surface),
                            "combined": next_combined,
                            "after_transport": asdict(next_step.after.transport),
                            "scale_audit": scale_journal,
                        }
                    )
                    current = next_step.after
                    surface = next_surface
                    carry = {
                        **carry,
                        **next_combined,
                        "external_water_loading_kg_m2": surface.liquid_water_kg_m2,
                    }
                except Exception as exc:
                    result["radial_continuation"].append(
                        {
                            "index": index,
                            "accepted": False,
                            "time_s": previous.transport.time_s,
                            "error": str(exc),
                            "traceback": traceback.format_exc(),
                            "exact_accepted_state_preserved": current is previous,
                            "study_debug": getattr(exc, "study_debug", None),
                            "optimizer_debug": getattr(exc, "optimizer_debug", None),
                        }
                    )
                    result["status"] = "DT_BIRTH_ACCEPTED_RADIAL_CONTINUATION_REFUSED"
                    break
            checkpoint_path = args.out.with_suffix(".accepted.pkl")
            if checkpoint_path.exists():
                raise ValueError("refusing to overwrite an accepted checkpoint")
            checkpoint_path.write_bytes(
                pickle.dumps(
                    {
                        "current": current,
                        "surface": surface,
                        "carry": carry,
                        "source_sha256": result["entry_source_sha256"],
                        "result_path": str(args.out),
                    },
                    protocol=5,
                )
            )
            result["accepted_checkpoint"] = {
                "path": str(checkpoint_path),
                "time_s": current.transport.time_s,
                "sha256": hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
            }
    except Exception as exc:
        if before is not None and not isinstance(exc, bi.BirthIntegratorStepError):
            wrapped = bi.BirthIntegratorStepError(
                f"birth preparation refused with exact fully-wet rollback: {exc}",
                before,
                nonlinear_evaluations=0,
                rejected_trial_evaluations=0,
            )
            wrapped.__cause__ = exc
            exc = wrapped
        result.update(
            status=(
                "ACTUAL_FEED_ZERO_OLD_DRY_VOLUME_BIRTH_REFUSED"
                if before is not None
                else "BIRTH_INPUT_SETUP_REFUSED"
            ),
            error_type=type(exc).__name__,
            error=str(exc),
            exact_rollback_identity=(
                getattr(exc, "rollback_state", None) is before if before is not None else None
            ),
            total_residual_assembly_attempts=getattr(exc, "total_residual_assembly_attempts", 0),
            study_debug=getattr(exc, "study_debug", None),
            traceback=traceback.format_exc(),
        )
    result["wall_s"] = time.perf_counter() - started
    # A provenance failure must preserve the completed/refused numerical
    # evidence; it must not discard the run by raising after computation.
    try:
        result["stamp"] = era_stamp(sources)
        result["source_changed_during_run"] = [
            str(path)
            for path in sources
            if result["entry_source_sha256"].get(str(path))
            != hashlib.sha256(path.read_bytes()).hexdigest()
        ]
        if result["source_changed_during_run"]:
            result["numerical_status_before_provenance_refusal"] = result["status"]
            result.update(accepted=False, status="BIRTH_SOURCE_CHANGED_DURING_RUN")
    except Exception as exc:
        result["provenance_error"] = f"{type(exc).__name__}: {exc}"
        result["provenance_traceback"] = traceback.format_exc()
        result["numerical_status_before_provenance_refusal"] = result["status"]
        result.update(accepted=False, status="BIRTH_PROVENANCE_REFUSED")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feed", type=Path, required=True)
    parser.add_argument(
        "--bath", type=Path, default=core.DESTINY_PATH.with_name("destiny_v3_exit105.json")
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dt-s", type=float, default=0.001)
    parser.add_argument("--cells", type=int)
    parser.add_argument(
        "--coefficient-case",
        choices=[c.identifier for c in core.p1ef.engineering_coefficient_cases()],
        default="nominal_log_center",
    )
    parser.add_argument("--continued-water", action="store_true")
    parser.add_argument("--budget", type=int, default=400)
    parser.add_argument("--front-z", type=float)
    parser.add_argument("--solver", choices=("scipy", "newton", "physical"), default="scipy")
    parser.add_argument(
        "--resume", type=Path, help="trusted local accepted checkpoint from this driver"
    )
    parser.add_argument("--root-seed", type=Path)
    parser.add_argument("--continuation-steps", type=int, default=0)
    parser.add_argument("--continuation-dt", type=float, default=0.01)
    parser.add_argument("--path-point-budget", type=int, default=20000)
    parser.add_argument("--wall-budget-s", type=float, default=60.0)
    options = parser.parse_args()
    options.feed = options.feed.resolve()
    options.bath = options.bath.resolve()
    options.out = options.out.resolve()
    if options.out.exists() or not options.out.resolve().is_relative_to(REPO_ROOT / "paper"):
        parser.error("write to a new output under paper")
    if (
        options.dt_s <= 0
        or options.budget < 2
        or options.path_point_budget < 1
        or options.wall_budget_s <= 0
    ):
        parser.error("positive duration and residual-assembly budget >=2 required")
    output = run(options)
    options.out.parent.mkdir(parents=True, exist_ok=True)
    options.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                key: output.get(key)
                for key in (
                    "status",
                    "accepted",
                    "error",
                    "total_residual_assembly_attempts",
                    "wall_s",
                )
            },
            indent=2,
        )
    )
