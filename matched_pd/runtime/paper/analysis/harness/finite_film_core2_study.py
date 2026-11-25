"""Unqualified finite-resistance boundary study around native Core2 cut balances.

The Dirichlet subclass is an integration carrier only: the context adapter
always replaces its exterior face by a solved binary film. It is deliberately
not a ReducedFilmPoreBoundary and does not claim its fast-mass audit. Internal
faces, inventories, RH equations, and acceptance/conservation checks are native.
"""

from contextlib import contextmanager
from dataclasses import dataclass, replace
import math
import inspect

import numpy as np
from scipy.optimize import least_squares

import dtdc_particle_destiny_march as core
import particle_numerical_repairs as numerics


@contextmanager
def magnitude_storage_scale():
    """Match fully-dry Core2's magnitude policy for a coordinate derivative.

    Only the GasOnlyPrimitiveDomain scale's final positivity check changes.
    Sampling, state topology checks, every balance and actual calorics remain
    unchanged. The compiled in-memory variant cannot update a production pin.
    """
    original = core.ci._dry_energy_capacity_with_audit
    source = inspect.getsource(original)
    needle = "\n    if not math.isfinite(capacity) or capacity <= 0.0:"
    if source.count(needle) != 1:
        raise RuntimeError("native storage scale changed; review the adapter")
    replacement = (
        "\n    _signed_storage_scale_audit.append((temperature, y_hexane, capacity))\n    capacity = abs(capacity)"
        + needle
    )
    namespace = dict(original.__globals__)
    journal = []
    namespace["_signed_storage_scale_audit"] = journal
    exec(
        compile(
            source.replace(needle, replacement), "<paper storage-scale magnitude adapter>", "exec"
        ),
        namespace,
    )
    core.ci._dry_energy_capacity_with_audit = namespace[original.__name__]
    try:
        yield journal
    finally:
        core.ci._dry_energy_capacity_with_audit = original


@contextmanager
def engineering_water_bridge():
    """Study reuse of C8 high-loading water activity, liquid-like calorics.

    Only explicitly tagged LuikovParams(W_cap=0.25) use this bridge. The
    qualified low-loading correlation and all production parameters remain
    unchanged. Both benchmarks must use this same declared material law.
    """
    activity_original = core.sp.water_activity
    retained_original = core.sp.water_retained
    original_cap = core.sp.LuikovParams().W_cap

    def activity(water, params):
        if params.W_cap != 0.25:
            return activity_original(water, params)
        base = replace(params, W_cap=original_cap)
        anchor = activity_original(original_cap, base)
        if water <= original_cap:
            return activity_original(water, base)
        if not original_cap < water <= params.W_cap:
            raise ValueError("water bridge loading outside its declared finite interval")
        return anchor + (1 - anchor) * (water - original_cap) / (params.W_cap - original_cap)

    def retained(aw, params):
        if params.W_cap != 0.25:
            return retained_original(aw, params)
        base = replace(params, W_cap=original_cap)
        anchor = activity_original(original_cap, base)
        if aw <= anchor:
            return retained_original(aw, base)
        if not anchor < aw <= 1.0:
            raise ValueError("water bridge activity outside (anchor,1]")
        return original_cap + (params.W_cap - original_cap) * (aw - anchor) / (1 - anchor)

    core.cp.clear_equilibrium_cache()
    core.sp.water_activity = activity
    core.sp.water_retained = retained
    try:
        yield
    finally:
        core.sp.water_activity = activity_original
        core.sp.water_retained = retained_original
        core.cp.clear_equilibrium_cache()


@dataclass(frozen=True)
class FiniteFilmStudyBoundary(core.ct.DirichletPoreBoundary):
    film_boundary: core.ct.ReducedFilmPoreBoundary | None = None
    accepted_beta_limit: float = 1.0
    actual_surface_path: bool = False
    numerical_cap_seam: bool = False

    def __post_init__(self):
        super().__post_init__()
        if self.film_boundary is None:
            raise ValueError("a declared bulk film boundary is required")
        if not 1.0 <= self.accepted_beta_limit <= 3.0:
            raise ValueError("explicit study beta limit must lie in [1,3]")
        if (self.temperature_k, self.y_hexane) != (
            self.film_boundary.temperature_k,
            self.film_boundary.y_hexane,
        ):
            raise ValueError("study carrier must retain exact bulk conditions")


@dataclass(frozen=True)
class FiniteFilmStudyAudit:
    boundary_model: str
    bulk_temperature_k: float
    bulk_y_hexane: float
    surface_temperature_k: float
    surface_y_hexane: float
    water_mass_flux_out_kg_m2_s: float
    hexane_mass_flux_out_kg_m2_s: float
    binary_mass_residual_kg_m2_s: float
    heat_residual_w_m2: float
    beta: float
    ackermann_factor: float
    nonlinear_evaluations: int
    normalized_residual: float
    fast_mass_approximation_used: bool = False
    physically_qualifying: bool = False
    endpoint_roundoff_audit: dict | None = None


def solve_surface(
    boundary, cell, distance, stefan, pore, config, radius, trial_coordinates=None, enforce=True
):
    """Solve Ts and an admissible gas-only surface coordinate, with fixed Ntot.

    Trial mass residual uses the algebraic partition helper, avoiding the
    logarithmic profile audit on unconverged trial fluxes. The accepted root
    separately passes that stricter profile audit and the full beta guard.
    """
    reference = core.ct.prepare_reduced_film_surface_state(
        boundary.film_boundary,
        cell_temperature_k=cell.temperature_k,
        pressure_pa=config.dry.pressure_pa,
        particle_radius_m=radius,
        outer_half_cell_distance_m=distance,
        particle_thermal_conductivity_w_m_k=config.dry.thermal_conductivity.value_w_m_k,
        binary_gas_interaction_k_wh=pore.k_wh,
    )
    pressure = config.dry.pressure_pa
    lower, upper = config.dry.conditioned_temperature_domain.solver_bounds_k
    # Finite gas-only intervals can disappear at low T. Start at the native
    # cell, whose equilibrium has already been validated by the radial solver.
    start_t = cell.temperature_k
    start_c = core.ci._encode_face_conditioned_dry_y(
        (start_t, start_t), pressure, cell.y_hexane, pore
    )
    mass_scale = max(
        reference.coefficients.film_density_kg_m3 * reference.coefficients.binary_mass_transfer_m_s,
        1e-8,
    )
    heat_scale = max(
        reference.coefficients.heat_transfer_w_m2_k
        * abs(boundary.temperature_k - cell.temperature_k),
        1.0,
    )

    def evaluate(x):
        ts, coordinate = map(float, x)
        ys = (
            core.cp.decode_gas_only_y(ts, pressure, coordinate, pore)
            if boundary.actual_surface_path
            else core.ci._decode_face_conditioned_dry_y(
                (ts, 0.5 * (ts + cell.temperature_k)), pressure, coordinate, pore
            )
        )
        face_config = (
            replace(
                config,
                moving_interface_composition_force_authority=core.cut.MovingInterfaceCompositionForceAuthority.ACTUAL_LINEAR_T_Y_PATH,
            )
            if boundary.actual_surface_path
            else config
        )
        face = core.cut._dry_flux(
            cell.temperature_k,
            cell.y_hexane,
            ts,
            ys,
            distance,
            stefan,
            pore,
            pore,
            face_config,
            thermodynamic_force_temperature_k=(
                cell.temperature_k if boundary.actual_surface_path else None
            ),
        )
        film = core.ct._reduced_film_coefficients_at_surface_temperature(reference, ts)
        flux = core.efr.ComponentMassFluxes(
            water_kg_m2_s=face.component.conserved_water_flux_mol_m2_s * core.wa.M,
            hexane_kg_m2_s=face.component.conserved_hexane_flux_mol_m2_s * core.hx.M,
        )
        predicted = core.efr.component_fluxes_for_supplied_total_mass_flux(
            film,
            total_mass_flux_kg_m2_s=flux.total_kg_m2_s,
            interface_hexane_mass_fraction=core._hexane_mass_fraction(ys),
            bulk_hexane_mass_fraction=core._hexane_mass_fraction(boundary.y_hexane),
        )
        capacities = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
            surface_temperature_k=ts,
            bulk_temperature_k=boundary.temperature_k,
            pressure_pa=pressure,
            y_hexane=ys,
            k_wh=pore.k_wh,
        )
        ack = core.efr.ackermann_heat_transfer(
            film,
            fluxes=flux,
            water_film_heat_capacity_j_kg_k=capacities.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=capacities.hexane_j_kg_k,
            # Study-only nonlinear continuation of the exact analytic factor.
            # The resolved root MUST pass the unchanged beta<=1 guard below.
            enforce_production_domain=False,
        )
        mass = flux.hexane_kg_m2_s - predicted.hexane_kg_m2_s
        heat = face.reduced_heat_flux_w_m2 + film.heat_transfer_w_m2_k * ack.factor * (
            boundary.temperature_k - ts
        )
        return np.array([mass / mass_scale, heat / heat_scale]), (
            ts,
            ys,
            face,
            film,
            flux,
            ack,
            mass,
            heat,
        )

    if trial_coordinates is None:
        solved = least_squares(
            lambda x: evaluate(x)[0],
            [start_t, start_c],
            bounds=(
                [max(lower, min(start_t, reference.surface_temperature_k)), -35.0],
                [min(upper, max(start_t, boundary.temperature_k)), 35.0],
            ),
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
            max_nfev=100,
            x_scale=[10.0, 1.0],
        )
        coordinates = solved.x
        evaluations = solved.nfev
    else:
        solved = None
        coordinates = trial_coordinates
        evaluations = 1
    residual, payload = evaluate(coordinates)
    maximum = float(np.max(np.abs(residual)))
    if enforce and ((solved is not None and not solved.success) or maximum > 1e-8):
        raise RuntimeError(
            f"finite-film surface failed: {solved.message}; scaled residual={maximum:.3g}; x={solved.x}"
        )
    ts, ys, face, film, flux, ack, mass, heat = payload
    capacities = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=ts,
        bulk_temperature_k=boundary.temperature_k,
        pressure_pa=pressure,
        y_hexane=ys,
        k_wh=pore.k_wh,
    )
    try:
        ack = core.efr.ackermann_heat_transfer(
            film,
            fluxes=flux,
            water_film_heat_capacity_j_kg_k=capacities.water_j_kg_k,
            hexane_film_heat_capacity_j_kg_k=capacities.hexane_j_kg_k,
            enforce_production_domain=True,
        )
    except core.efr.ReducedFilmValidityError as exc:
        if enforce and abs(ack.beta) > boundary.accepted_beta_limit:
            raise RuntimeError(
                f"resolved finite-film root fails explicit study band: {exc}; beta={ack.beta:.12g}, limit={boundary.accepted_beta_limit}, Ts={ts:.12g}, ys={ys:.12g}, jw={flux.water_kg_m2_s:.12g}, jh={flux.hexane_kg_m2_s:.12g}, scaled_residual={maximum:.3g}"
            ) from exc
    endpoint_roundoff = None
    if enforce:
        endpoint_roundoff = numerics.audit_film_profile(
            film,
            interface_hexane_mass_fraction=core._hexane_mass_fraction(ys),
            bulk_hexane_mass_fraction=core._hexane_mass_fraction(boundary.y_hexane),
            fluxes=flux,
            closure_residual_budget_kg_m2_s=1e-12*film.film_density_kg_m3*film.binary_mass_transfer_m_s,
        )
    thermo = core.cp.evaluate_equilibrium(ts, pressure, ys, pore)
    # Every conserved component exits the particle into the gas bath. Its
    # common-datum gas enthalpy is transported exactly once, as in native film.
    water_h = face.component.conserved_water_flux_mol_m2_s * thermo.water_gas_partial_enthalpy_j_mol
    hexane_h = (
        face.component.conserved_hexane_flux_mol_m2_s * thermo.hexane_gas_partial_enthalpy_j_mol
    )
    energy = core.cp.ComponentEnergyFlux(
        conductive_heat_flux_w_m2=face.reduced_heat_flux_w_m2,
        gas_water_enthalpy_flux_w_m2=water_h,
        gas_hexane_enthalpy_flux_w_m2=hexane_h,
        retained_water_enthalpy_flux_w_m2=0.0,
        total_energy_flux_w_m2=math.fsum((face.reduced_heat_flux_w_m2, water_h, hexane_h)),
    )
    audit = FiniteFilmStudyAudit(
        "finite binary Stefan / exact PHY053 Ackermann study",
        boundary.temperature_k,
        boundary.y_hexane,
        ts,
        ys,
        flux.water_kg_m2_s,
        flux.hexane_kg_m2_s,
        mass,
        heat,
        ack.beta,
        ack.factor,
        evaluations,
        maximum,
        endpoint_roundoff_audit=endpoint_roundoff,
    )
    return replace(face, energy=energy, surface_film_audit=audit)


@contextmanager
def installed_boundary_adapter(surface_coordinates=None, numerical_cap_seam=False):
    original = core.cut._dry_face_fluxes
    path_module = core.cut.acfp
    original_seam = path_module._bisect_cap_seam

    def representable_seam(activities, lower, upper, activity_at_cap, isolation_depth):
        try:
            return original_seam(activities, lower, upper, activity_at_cap, isolation_depth)
        except core.cp.CoupledPoreTopologyError as exc:
            if "representable path point" not in str(exc):
                raise
            lv = path_module._cap_residual(activities, lower, activity_at_cap)
            uv = path_module._cap_residual(activities, upper, activity_at_cap)
            if lv * uv >= 0:
                raise
            best = min(((lower, lv), (upper, uv)), key=lambda p: abs(p[1]))
            for iterations in range(1, 81):
                mid = lower + 0.5 * (upper - lower)
                if mid == lower or mid == upper:
                    break
                mv = path_module._cap_residual(activities, mid, activity_at_cap)
                if abs(mv) < abs(best[1]):
                    best = (mid, mv)
                if lv * mv < 0:
                    upper, uv = mid, mv
                else:
                    lower, lv = mid, mv
            if upper - lower > 1e-12 or abs(best[1]) > 1e-10:
                raise core.cp.CoupledPoreTopologyError(
                    f"study seam remains unresolved: activity residual={best[1]}, width={upper-lower}"
                ) from exc
            # Keep the kink and branch partition. Only the located seam's
            # representability criterion changes, with a recorded residual.
            return path_module.CapSeamAudit(
                best[0], best[1], upper - lower, isolation_depth, iterations
            )

    def adapted(
        geometry, layout, cells, interface, stefan_fluxes, boundary, oil_labels, config, **kwargs
    ):
        if not isinstance(boundary, FiniteFilmStudyBoundary):
            return original(
                geometry,
                layout,
                cells,
                interface,
                stefan_fluxes,
                boundary,
                oil_labels,
                config,
                **kwargs,
            )
        index = layout.dry_cell_indices[-1]
        pore = replace(config.dry.pore, w_o=oil_labels[index])
        distance = geometry.master_grid.R - core.cut._dry_piece_center(
            geometry.cells[index], geometry.front.radius_m
        )
        exterior = solve_surface(
            boundary,
            cells[-1],
            distance,
            stefan_fluxes[-1],
            pore,
            config,
            geometry.master_grid.R,
            trial_coordinates=(surface_coordinates["x"] if surface_coordinates else None),
            enforce=(
                surface_coordinates.get("enforce", False)
                if surface_coordinates
                else kwargs.get("enforce_reduced_film_thresholds", True)
            ),
        )
        audit = exterior.surface_film_audit
        resolved = core.ct.DirichletPoreBoundary(
            (
                cells[-1].temperature_k
                if boundary.actual_surface_path
                else audit.surface_temperature_k
            ),
            cells[-1].y_hexane if boundary.actual_surface_path else audit.surface_y_hexane,
            "resolved finite-film surface; native interior faces",
        )
        faces = original(
            geometry,
            layout,
            cells,
            interface,
            stefan_fluxes,
            resolved,
            oil_labels,
            config,
            **kwargs,
        )
        return (*faces[:-1], exterior)

    core.cut._dry_face_fluxes = adapted
    if numerical_cap_seam:
        path_module._bisect_cap_seam = representable_seam
    try:
        yield
    finally:
        core.cut._dry_face_fluxes = original
        path_module._bisect_cap_seam = original_seam


def advance_joint_study(before, dt, boundary, seed, *, surface_seed=None):
    """Solve native radial unknowns and two finite-film unknowns together.

    This avoids demanding an exact surface root for every arbitrary radial
    trial. No failed surface residual is replaced by a successful audit.
    """
    chart = core.ci._coordinate_chart(before)
    native_x = core.ci._encode_candidate_for_chart(
        seed.candidate, chart, seed.stefan_flux_scales_mol_m2_s
    )
    # A finite-resistance surface root is just an initial guess here.
    holder = {"x": None, "enforce": False}
    with installed_boundary_adapter(numerical_cap_seam=boundary.numerical_cap_seam):
        initial = core.cut.assemble_backward_euler(
            before.transport, seed.candidate, dt, boundary, enforce_reduced_film_thresholds=False
        )
    audit = initial.surface_film_audit
    pore = replace(before.transport.config.dry.pore, w_o=before.transport.oil_fraction_labels[-1])
    tc = seed.candidate.dry_temperatures_k[-1]
    cs = (
        core.cp.encode_gas_only_y(
            audit.surface_temperature_k,
            before.transport.config.dry.pressure_pa,
            audit.surface_y_hexane,
            pore,
        )
        if boundary.actual_surface_path
        else core.ci._encode_face_conditioned_dry_y(
            (audit.surface_temperature_k, 0.5 * (audit.surface_temperature_k + tc)),
            before.transport.config.dry.pressure_pa,
            audit.surface_y_hexane,
            pore,
        )
    )
    if surface_seed is not None:
        surface_temperature, surface_y = surface_seed
        cs = core.cp.encode_gas_only_y(surface_temperature,
              before.transport.config.dry.pressure_pa,surface_y,pore)
    else:
        surface_temperature = audit.surface_temperature_k
    x0 = np.r_[native_x, surface_temperature, cs]
    scales = core.ci._residual_scales(before.transport, initial, dt)
    base_boundary = boundary.film_boundary
    pair = base_boundary.correlated_pair
    mass_scale = max(pair.binary_mass_transfer_m_s, 1e-8)
    heat_scale = max(pair.heat_transfer_w_m2_k * abs(boundary.temperature_k - tc), 1.0)
    work = 0

    def residual(x):
        nonlocal work
        holder["x"] = x[-2:]
        candidate = core.ci._decode_candidate_for_chart(
            x[:-2], before.transport.layout, chart, seed.stefan_flux_scales_mol_m2_s
        )
        assembly = core.cut.assemble_backward_euler(
            before.transport, candidate, dt, boundary, enforce_reduced_film_thresholds=False
        )
        work += 1
        film = assembly.surface_film_audit
        native = (
            np.asarray(core.ci._datum_covariant_residual_vector(assembly.residuals)) / scales.vector
        )
        return np.r_[
            native,
            film.binary_mass_residual_kg_m2_s / mass_scale,
            film.heat_residual_w_m2 / heat_scale,
        ]

    lb, ub = core.ci._coordinate_solver_bounds_for_chart(
        before.transport.layout, before.controls.maximum_logit_magnitude, chart
    )
    tl, tu = before.transport.config.dry.conditioned_temperature_domain.solver_bounds_k
    with installed_boundary_adapter(holder, numerical_cap_seam=boundary.numerical_cap_seam):
        solution = least_squares(
            residual,
            x0,
            bounds=(np.r_[lb, tc, -30.0], np.r_[ub, min(tu, boundary.temperature_k), 30.0]),
            x_scale="jac",
            max_nfev=400,
            jac="3-point",
            diff_step=1e-4,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-11,
        )
        values = residual(solution.x)
        maximum = float(np.max(np.abs(values)))
        if not solution.success or maximum > before.controls.nonlinear_residual_tolerance:
            error = RuntimeError(
                f"joint radial/film study not converged; residual={maximum:.12g}; evaluations={work}; {solution.message}"
            )
            from dataclasses import asdict

            error.study_debug = {
                "scaled_residuals": values.tolist(),
                "surface_unknowns": solution.x[-2:].tolist(),
                "candidate": asdict(
                    core.ci._decode_candidate_for_chart(
                        solution.x[:-2],
                        before.transport.layout,
                        chart,
                        seed.stefan_flux_scales_mol_m2_s,
                    )
                ),
                "work": work,
                "jacobian_singular_values": np.linalg.svd(solution.jac, compute_uv=False).tolist(),
                "optimizer_optimality": float(solution.optimality),
                "active_bounds": solution.active_mask.tolist(),
                "coordinates": solution.x.tolist(),
                "residual_scales": list(scales.vector),
            }
            raise error
        holder["enforce"] = True
        candidate = core.ci._decode_candidate_for_chart(
            solution.x[:-2], before.transport.layout, chart, seed.stefan_flux_scales_mol_m2_s
        )
        assembly = core.cut.assemble_backward_euler(before.transport, candidate, dt, boundary)
        condition = float(np.linalg.cond(solution.jac))
        if condition > before.controls.maximum_condition_proxy:
            raise RuntimeError(f"joint study Jacobian condition exceeds native limit: {condition}")
        step = core.ci._accept_step(
            before,
            assembly,
            boundary,
            seed,
            scales,
            scales,
            tuple(values[:-2]),
            dt,
            work,
            0,
            str(solution.message),
            condition,
            chart,
            None,
            "explicit joint finite-film dense study",
            "trf",
        )
        if (
            max(
                step.ledger.maximum_step_ledger_residual,
                step.ledger.maximum_cumulative_ledger_residual,
            )
            > before.controls.ledger_tolerance
        ):
            raise RuntimeError("joint study fails unchanged native conservation tolerances")
        if assembly.ledger.minimum_entropy_production_w_m3_k < 0.0:
            raise RuntimeError("joint study has negative native entropy production")
        if hasattr(solution,"physical_coordinate_audit"):
            step = replace(step,ledger=replace(step.ledger,
                    nonlinear_algorithm="domain_aware_damped_newton_physical_chart"))
        for state in (before.transport, step.after.transport):
            core.ci.certify_explicit_wet_retained_water_applicability(
                before.controls,
                state.wet_retained_water_loadings,
                state.effective_wet_retained_water_capacity_duals_over_rt,
                state.config.wet.luikov,
            )
        return step
