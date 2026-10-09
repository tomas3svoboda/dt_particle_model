"""Conservative, zero-thickness external water-film onset study.

Independent film and pore Stefan flows meet at a surface node with liquid-water
mass and common-datum energy storage. Hexane storage/aqueous resistance is
neglected: this is the fully accessible thin-film limit, not a qualified law.
"""

from contextlib import contextmanager, nullcontext
from dataclasses import asdict, dataclass, replace
import math

import numpy as np
from scipy.optimize import least_squares

import finite_film_core2_study as study

core = study.core


def node_balances(
    dt,
    liquid,
    old_liquid,
    liquid_h,
    old_h,
    pore_water,
    pore_hexane,
    pore_energy,
    film_water,
    film_hexane,
    film_energy,
):
    """Mass and common-datum energy conservation of the explicit surface node."""
    return (
        (liquid - old_liquid) / dt - pore_water + film_water,
        film_hexane - pore_hexane,
        math.fsum(((liquid * liquid_h - old_liquid * old_h) / dt, -pore_energy, film_energy)),
    )


@dataclass(frozen=True)
class WaterFilmNodeAudit:
    temperature_k: float
    surface_y_hexane: float
    liquid_water_kg_m2: float
    external_water_flux_kg_m2_s: float
    external_hexane_flux_kg_m2_s: float
    external_energy_flux_w_m2: float
    water_node_residual_kg_m2_s: float
    hexane_node_residual_kg_m2_s: float
    energy_node_residual_w_m2: float
    water_activity_trace_deficit: float
    beta: float
    physically_qualifying: bool = False


def evaluate_water_surface(
    holder,
    dt,
    old_loading,
    old_temperature,
    cell,
    distance,
    stefan,
    pore,
    config,
    radius,
    boundary,
    surface_y_hexane=None,
    surface_liquid_transfer_adapter=None,
    surface_liquid_transfer_fraction=None,
):
    ts, total_mass_flux, liquid = holder["x"]
    ys = (
        core.cp.decode_gas_only_y(ts, config.dry.pressure_pa, -32.0, pore)
        if surface_y_hexane is None
        else surface_y_hexane
    )
    thermo = core.cp.evaluate_equilibrium(ts, config.dry.pressure_pa, ys, pore)
    exterior = core.cut._dry_flux(
        cell.temperature_k,
        cell.y_hexane,
        ts,
        ys,
        distance,
        stefan,
        pore,
        pore,
        replace(
            config,
            moving_interface_composition_force_authority=core.cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH,
        ),
        thermodynamic_force_temperature_k=cell.temperature_k,
    )
    reference = core.ct.prepare_reduced_film_surface_state(
        boundary.film_boundary,
        cell_temperature_k=cell.temperature_k,
        pressure_pa=config.dry.pressure_pa,
        particle_radius_m=radius,
        outer_half_cell_distance_m=distance,
        particle_thermal_conductivity_w_m_k=config.dry.thermal_conductivity.value_w_m_k,
        binary_gas_interaction_k_wh=pore.k_wh,
    )
    film = core.ct._reduced_film_coefficients_at_surface_temperature(reference, ts)
    components = core.efr.component_fluxes_for_supplied_total_mass_flux(
        film,
        total_mass_flux_kg_m2_s=total_mass_flux,
        interface_hexane_mass_fraction=core._hexane_mass_fraction(ys),
        bulk_hexane_mass_fraction=core._hexane_mass_fraction(boundary.y_hexane),
    )
    capacities = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=ts,
        bulk_temperature_k=boundary.temperature_k,
        pressure_pa=config.dry.pressure_pa,
        y_hexane=ys,
        k_wh=pore.k_wh,
    )
    ack = core.efr.ackermann_heat_transfer(
        film,
        fluxes=components,
        water_film_heat_capacity_j_kg_k=capacities.water_j_kg_k,
        hexane_film_heat_capacity_j_kg_k=capacities.hexane_j_kg_k,
        enforce_production_domain=False,
    )
    q_out = film.heat_transfer_w_m2_k * ack.factor * (ts - boundary.temperature_k)
    external_energy = math.fsum(
        (
            q_out,
            components.water_kg_m2_s * thermo.water_gas_partial_enthalpy_j_mol / core.wa.M,
            components.hexane_kg_m2_s * thermo.hexane_gas_partial_enthalpy_j_mol / core.hx.M,
        )
    )
    liquid_h = core.wa.state_Tp(ts, config.dry.pressure_pa, "liquid").h_mass
    old_h = core.wa.state_Tp(old_temperature, config.dry.pressure_pa, "liquid").h_mass
    if surface_liquid_transfer_adapter is not None:
        exterior = surface_liquid_transfer_adapter(
            exterior,
            old_liquid_loading_kg_m2=old_loading,
            dt_s=dt() if callable(dt) else dt,
            external_water_flux_kg_m2_s=components.water_kg_m2_s,
            requested_share=surface_liquid_transfer_fraction,
        )
    pore_water = exterior.component.conserved_water_flux_mol_m2_s * core.wa.M
    pore_hexane = exterior.component.conserved_hexane_flux_mol_m2_s * core.hx.M
    r_water, r_hexane, r_energy = node_balances(
        dt() if callable(dt) else dt,
        liquid,
        old_loading,
        liquid_h,
        old_h,
        pore_water,
        pore_hexane,
        exterior.energy.total_energy_flux_w_m2,
        components.water_kg_m2_s,
        components.hexane_kg_m2_s,
        external_energy,
    )
    audit = WaterFilmNodeAudit(
        ts,
        ys,
        liquid,
        components.water_kg_m2_s,
        components.hexane_kg_m2_s,
        external_energy,
        r_water,
        r_hexane,
        r_energy,
        1.0 - thermo.water_activity,
        ack.beta,
    )
    return replace(exterior, surface_film_audit=audit)


@contextmanager
def node_adapter(holder, dt, old_loading=0.0, old_temperature=343.0):
    original = core.cut._dry_face_fluxes

    def faces(geometry, layout, cells, interface, stefan, boundary, oil_labels, config, **kwargs):
        if not isinstance(boundary, study.FiniteFilmStudyBoundary):
            return original(
                geometry, layout, cells, interface, stefan, boundary, oil_labels, config, **kwargs
            )
        last = layout.dry_cell_indices[-1]
        pore = replace(config.dry.pore, w_o=oil_labels[last])
        distance = geometry.master_grid.R - core.cut._dry_piece_center(
            geometry.cells[last], geometry.front.radius_m
        )
        exterior = evaluate_water_surface(
            holder,
            dt,
            old_loading,
            old_temperature,
            cells[-1],
            distance,
            stefan[-1],
            pore,
            config,
            geometry.master_grid.R,
            boundary,
            surface_y_hexane=holder.get("y_hexane"),
        )
        resolved = core.ct.DirichletPoreBoundary(
            cells[-1].temperature_k,
            cells[-1].y_hexane,
            "native interior assembly; exterior replaced by water node",
        )
        native = original(
            geometry, layout, cells, interface, stefan, resolved, oil_labels, config, **kwargs
        )
        return (*native[:-1], exterior)

    core.cut._dry_face_fluxes = faces
    try:
        yield
    finally:
        core.cut._dry_face_fluxes = original


def advance_water_node(
    before,
    dt,
    boundary,
    seed,
    *,
    old_loading=0.0,
    old_temperature=343.0,
    surface_seed=None,
    physical_solver=False,
    residual_budget=1200,
    beta_decoder=None,
    energy_context=None,
    target_z=None,
    duration_bounds=(1e-8, 10.0),
    boundary_factory=None,
    gas_only_surface=False,
    maximum_restart_trial_residual=1e-5,
    difference_mode="forward",
):
    chart = core.ci._coordinate_chart(before)
    x_native = core.ci._encode_candidate_for_chart(
        seed.candidate, chart, seed.stefan_flux_scales_mol_m2_s
    )
    tc = seed.candidate.dry_temperatures_k[-1]
    holder = {
        "x": np.array([tc + 0.6, 0.03, 1e-5] if surface_seed is None else surface_seed),
        "dt": dt,
    }
    if gas_only_surface:
        if old_loading != 0.0:
            raise ValueError("absent-film surface requires exactly zero incoming liquid")
        if surface_seed is None:
            raise ValueError("absent-film surface requires a declared numerical seed")
        pore = replace(
            before.transport.config.dry.pore, w_o=before.transport.oil_fraction_labels[-1]
        )

        def configure_gas_surface(values):
            ts, flux_or_beta, q = values
            if not 0 < q < 1:
                raise ValueError("surface composition fraction outside open gas phase chart")
            lo, hi = core.cp.gas_only_composition_interval(
                ts, before.transport.config.dry.pressure_pa, pore
            ).y_hexane_bounds
            holder["y_hexane"] = lo + (hi - lo) * q
            holder["x"] = np.array([ts, flux_or_beta, 0.0])

        configure_gas_surface(surface_seed)
    fixed_front_coordinate = float(x_native[-2])
    if target_z is not None:
        if (
            seed.candidate.front_z != target_z
            or not 0 < target_z < before.transport.geometry.front.z
        ):
            raise ValueError("positive target must match seed and recede from accepted state")
        x_native[-2] = core.ci._encode_open(dt, duration_bounds)
    active_boundary = boundary
    if beta_decoder is not None and surface_seed is not None:
        holder["x"][1] = (
            beta_decoder(float(surface_seed[0]), float(surface_seed[1]), holder["y_hexane"])
            if gas_only_surface
            else (
                beta_decoder(float(surface_seed[0]), float(surface_seed[1]))
                if boundary_factory is None
                else beta_decoder(float(surface_seed[0]), float(surface_seed[1]), active_boundary)
            )
        )

    def decode_candidate(values):
        encoded = np.array(values, copy=True)
        if target_z is not None:
            holder["dt"] = core.ci._decode_open(encoded[-2], duration_bounds)
            encoded[-2] = fixed_front_coordinate
        c = core.ci._decode_candidate_for_chart(
            encoded, before.transport.layout, chart, seed.stefan_flux_scales_mol_m2_s
        )
        return c if target_z is None else replace(c, front_z=target_z)

    work = 0
    with (
        study.installed_boundary_adapter(numerical_cap_seam=True),
        nullcontext() if energy_context is None else energy_context(holder),
    ):
        with node_adapter(holder, lambda: holder["dt"], old_loading, old_temperature):
            initial = core.cut.assemble_backward_euler(
                before.transport,
                seed.candidate,
                dt,
                boundary,
                enforce_reduced_film_thresholds=False,
            )
            scales = core.ci._residual_scales(before.transport, initial, dt)
            scale_vector = np.asarray(scales.vector)
            mass_scale = max(old_loading / dt, 0.05)
            heat_scale = max(
                old_loading
                * abs(
                    core.wa.state_Tp(
                        old_temperature, before.transport.config.dry.pressure_pa, "liquid"
                    ).u_mass
                )
                / dt,
                5000.0,
            )

            def residual(x):
                nonlocal work, active_boundary
                candidate = decode_candidate(x[:-3])
                if boundary_factory is not None:
                    active_boundary = boundary_factory(holder["dt"])
                if gas_only_surface:
                    configure_gas_surface(x[-3:])
                else:
                    holder["x"] = np.array(x[-3:], copy=True)
                if beta_decoder is not None:
                    holder["x"][1] = (
                        beta_decoder(float(x[-3]), float(x[-2]), holder["y_hexane"])
                        if gas_only_surface
                        else (
                            beta_decoder(float(x[-3]), float(x[-2]))
                            if boundary_factory is None
                            else beta_decoder(float(x[-3]), float(x[-2]), active_boundary)
                        )
                    )
                assembly = core.cut.assemble_backward_euler(
                    before.transport,
                    candidate,
                    holder["dt"],
                    active_boundary,
                    enforce_reduced_film_thresholds=False,
                )
                work += 1
                a = assembly.surface_film_audit
                native = (
                    np.asarray(core.ci._datum_covariant_residual_vector(assembly.residuals))
                    / scale_vector
                )
                return np.r_[
                    native,
                    a.water_node_residual_kg_m2_s / mass_scale,
                    a.hexane_node_residual_kg_m2_s / mass_scale,
                    a.energy_node_residual_w_m2 / heat_scale,
                ]

            lb, ub = core.ci._coordinate_solver_bounds_for_chart(
                before.transport.layout, before.controls.maximum_logit_magnitude, chart
            )
            tl, tu = before.transport.config.dry.conditioned_temperature_domain.solver_bounds_k
            solver = least_squares
            surface_initial = holder["x"].copy()
            if gas_only_surface:
                surface_initial = np.asarray(surface_seed, dtype=float).copy()
            if beta_decoder is not None:
                surface_initial[1] = initial.surface_film_audit.beta
            liquid_hi = 1.0 if gas_only_surface else max(0.05, 4 * old_loading)
            coordinate_limit = boundary.accepted_beta_limit if beta_decoder is not None else 1.0
            if physical_solver:
                layout = before.transport.layout
                nw, nd, nf = layout.wet_piece_count, layout.dry_piece_count, layout.dry_face_count
                pressure = before.transport.config.dry.pressure_pa

                def interval(t, p):
                    return core.cp.gas_only_composition_interval(t, pressure, p).y_hexane_bounds

                def decode(v):
                    c = decode_candidate(v[:-3])
                    fractions = []
                    for t, y, p in zip(c.dry_temperatures_k, c.dry_y_hexane, chart.dry_pores):
                        lo, hi = interval(t, p)
                        fractions.append((y - lo) / (hi - lo))
                    c = replace(c, dry_y_hexane=tuple(fractions))
                    if target_z is not None:
                        c = replace(c, front_z=holder["dt"])
                    return np.r_[c.vector(), v[-3:]]

                def encode(v):
                    c = core.cut.CutTransportUnknowns.from_vector(layout, v[:-3])
                    duration = c.front_z
                    if target_z is not None:
                        c = replace(c, front_z=target_z)
                    ys = []
                    for t, q, p in zip(c.dry_temperatures_k, c.dry_y_hexane, chart.dry_pores):
                        if not 0 < q < 1:
                            raise ValueError("gas interval fraction outside its open domain")
                        lo, hi = interval(t, p)
                        ys.append(lo + (hi - lo) * q)
                    c = replace(c, dry_y_hexane=tuple(ys))
                    encoded = core.ci._encode_candidate_for_chart(
                        c, chart, seed.stefan_flux_scales_mol_m2_s
                    )
                    if target_z is not None:
                        encoded[-2] = core.ci._encode_open(duration, duration_bounds)
                    return np.r_[encoded, v[-3:]]

                pairs = (
                    [chart.wet_temperature] * nw
                    + [chart.wet_water] * nw
                    + [chart.dry_temperature] * nd
                    + [(0.0, 1.0)] * nd
                    + [(-math.inf, math.inf)] * nf
                    + [
                        chart.front_z if target_z is None else duration_bounds,
                        chart.interface_temperature,
                    ]
                    + [
                        (tl, min(tu, boundary.temperature_k)),
                        (-coordinate_limit, coordinate_limit),
                        (0.0, liquid_hi),
                    ]
                )

                def solver(fun, x0, **kwargs):
                    return study.numerics.physical_chart_solve(
                        fun,
                        x0,
                        decode=decode,
                        encode=encode,
                        physical_bounds=tuple(zip(*pairs)),
                        temperature_indices=[
                            *range(nw),
                            *range(2 * nw, 2 * nw + nd),
                            len(pairs) - 4,
                            len(pairs) - 3,
                        ],
                        **{
                            **kwargs,
                            "max_nfev": residual_budget,
                            "difference_mode": difference_mode,
                            "gtol": before.controls.nonlinear_residual_tolerance,
                        },
                    )

            solve_kwargs = dict(
                bounds=(
                    np.r_[lb, tl, -coordinate_limit, 0.0],
                    np.r_[ub, min(tu, boundary.temperature_k), coordinate_limit, liquid_hi],
                ),
                x_scale="jac",
                jac="3-point",
                diff_step=1e-4,
                max_nfev=400,
                xtol=1e-11,
                ftol=1e-11,
                gtol=1e-11,
            )
            restart_audit = []
            start = np.r_[x_native, surface_initial]
            total_calls = 0
            # A restart is another search for the SAME backward-Euler root.
            # It imports no accepted state or elapsed physical time. Native
            # normalization is rebuilt from this interval's trial assembly.
            for cycle in range(3 if physical_solver else 1):
                remaining = residual_budget - total_calls
                if physical_solver:
                    cycle_budget = (
                        min(remaining, max(200, residual_budget // 3)) if cycle < 2 else remaining
                    )

                    def cycle_solver(fun, x0, **kwargs):
                        return study.numerics.physical_chart_solve(
                            fun,
                            x0,
                            decode=decode,
                            encode=encode,
                            physical_bounds=tuple(zip(*pairs)),
                            temperature_indices=[
                                *range(nw),
                                *range(2 * nw, 2 * nw + nd),
                                len(pairs) - 4,
                                len(pairs) - 3,
                            ],
                            **{
                                **kwargs,
                                "max_nfev": cycle_budget,
                                "difference_mode": difference_mode,
                                "gtol": before.controls.nonlinear_residual_tolerance,
                            },
                        )

                else:
                    cycle_solver = solver
                try:
                    solution = cycle_solver(residual, start, **solve_kwargs)
                    total_calls += int(solution.nfev)
                    if solution.success:
                        break
                    trial = solution.x
                    maximum = float(np.max(np.abs(solution.fun)))
                except RuntimeError as exc:
                    if "computational budget exhausted" in str(exc):
                        raise
                    debug = getattr(exc, "optimizer_debug", None)
                    if not physical_solver or debug is None:
                        raise
                    total_calls += debug["residual_calls"]
                    trial = np.asarray(debug["best_native_coordinates"])
                    maximum = debug["best_maximum_scaled_residual"]
                    if cycle == 2 or total_calls >= residual_budget:
                        debug["same_interval_restarts"] = restart_audit
                        raise
                if not physical_solver or cycle == 2 or total_calls >= residual_budget:
                    break
                if maximum > maximum_restart_trial_residual:
                    raise RuntimeError(f"same-interval restart refused far from root: {maximum}")
                residual(trial)
                trial_candidate = decode_candidate(trial[:-3])
                trial_assembly = core.cut.assemble_backward_euler(
                    before.transport,
                    trial_candidate,
                    holder["dt"],
                    active_boundary,
                    enforce_reduced_film_thresholds=False,
                )
                scales = core.ci._residual_scales(before.transport, trial_assembly, holder["dt"])
                scale_vector = np.asarray(scales.vector)
                restart_audit.append(
                    {
                        "cycle": cycle,
                        "prior_maximum_scaled_residual": maximum,
                        "declared_maximum_trial_residual_for_seed_reuse": maximum_restart_trial_residual,
                        "physical_before_state_changed": False,
                        "residual_calls_so_far": total_calls,
                        "normalization": list(scales.vector),
                    }
                )
                start = trial.copy()
            solution.same_interval_restarts = restart_audit
            solution.nfev = total_calls
            values = residual(solution.x)
            maximum = float(np.max(np.abs(values)))
            candidate = decode_candidate(solution.x[:-3])
            assembly = core.cut.assemble_backward_euler(
                before.transport,
                candidate,
                holder["dt"],
                active_boundary,
                enforce_reduced_film_thresholds=False,
            )
            if (
                not solution.success
                or np.max(np.abs(values)) > before.controls.nonlinear_residual_tolerance
            ):
                error = RuntimeError(
                    f"water-node study not converged: residual={np.max(np.abs(values))}; work={work}"
                )
                error.study_debug = {
                    "candidate": asdict(candidate),
                    "surface": asdict(assembly.surface_film_audit),
                    "scaled_residuals": values.tolist(),
                    "work": work,
                    "same_interval_restarts": restart_audit,
                }
                raise error
            if abs(assembly.surface_film_audit.beta) > boundary.accepted_beta_limit:
                raise RuntimeError("resolved water-node root exceeds explicit beta study band")
            condition = float(np.linalg.cond(getattr(solution, "native_jac", solution.jac)))
            if condition > before.controls.maximum_condition_proxy:
                raise RuntimeError(
                    f"water-node Jacobian condition exceeds native limit: {condition}"
                )
            step = core.ci._accept_step(
                before,
                assembly,
                active_boundary,
                seed,
                scales,
                scales,
                tuple(values[:-3]),
                holder["dt"],
                work,
                0,
                str(solution.message),
                condition,
                chart,
                None,
                "joint native radial / conservative external water node study",
                "trf",
            )
            if (
                max(
                    step.ledger.maximum_step_ledger_residual,
                    step.ledger.maximum_cumulative_ledger_residual,
                )
                > before.controls.ledger_tolerance
            ):
                error = RuntimeError(
                    "water-node native particle ledger exceeds unchanged tolerances"
                )
                error.study_debug = {
                    "native_particle_ledger": asdict(step.ledger),
                    "maximum_scaled_root_residual": maximum,
                }
                raise error
            if assembly.ledger.minimum_entropy_production_w_m3_k < 0:
                raise RuntimeError("negative native entropy production")
            if physical_solver:
                step = replace(
                    step,
                    ledger=replace(
                        step.ledger,
                        nonlinear_algorithm="domain_aware_damped_newton_native_gas_fraction_chart",
                    ),
                )
            for state in (before.transport, step.after.transport):
                core.ci.certify_explicit_wet_retained_water_applicability(
                    before.controls,
                    state.wet_retained_water_loadings,
                    state.effective_wet_retained_water_capacity_duals_over_rt,
                    state.config.wet.luikov,
                )
            return step
