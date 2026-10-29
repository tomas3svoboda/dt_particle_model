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
import particle_kernel_envelope as envelope
import particle_kernel_pd_boundary as pd_boundary

core = study.core

# Kernel-envelope FRONT_STAGE_NOISE_AWARE_ENERGY_GATE (armed only; Class B, off by default, NOT
# adopted: open decision 1 of GT_PS2_UNIVERSAL_KERNEL_RULING_PACKET_8_FRONT_STAGE_NUMERICS_AND_
# STALLED_TRANSIENT_2026-10-28; the three-knob build of 2026-10-28/29, record
# docs/GT_PS2_PAPER1_THREE_HARNESS_KNOBS_BUILD_2026-10-29.md, built from variant (A) of
# docs/GT_PS2_UNIVERSAL_KERNEL_THIN_BIRTH_FIRST_STEP_MECHANISM_2026-10-28.md section 12). The owner's
# ruling GT_PS2_OWNER_RULING_PAPER1_GATE_UNITS_2026-10-17 point 2 (the dry stage's noise-aware
# energy-row gate, experiment_core2_dt_dry.NoiseAwareEnergyRowAcceptance) ported to the front-stage
# film-node search: the declared multiple of the row's representation floor.
FRONT_STAGE_ENERGY_GATE_MULTIPLE = 2.0  # WAS absent until 2026-10-28 (every row at its target)


def front_energy_rows(nw, nd, n_rows, stall_row):
    """The front search's energy rows: the dry-piece (shell) energy rows, the Rankine-Hugoniot
    (jump) energy row and the surface node's energy row, in the full row layout of n_rows rows
    (the native rows plus the three node rows); with a stalled front, renumbered for the searched
    rows (the consistency row removed)."""
    rows = [2 * nw + 2 * nd + i for i in range(nd)] + [2 * nw + 3 * nd + 2, n_rows - 1]
    if stall_row is None:
        return rows
    return [r - (1 if r > stall_row else 0) for r in rows if r != stall_row]


class FrontNoiseAwareEnergyRowAcceptance:
    """FRONT_STAGE_NOISE_AWARE_ENERGY_GATE's per-row acceptance clause for the domain-aware Newton.

    Every searched row passes at the target, except the energy rows, which pass at
    max(target, m x floor_i), floor_i = sum_j |J_ij| numpy.spacing(|x_j|), J the search's physical
    Jacobian at the iterate, x = the iterate plus the temperature origins (the physical
    coordinates). The interface of experiment_core2_dt_dry.NoiseAwareEnergyRowAcceptance.
    """

    def __init__(self, target, energy_rows, origins, multiple=FRONT_STAGE_ENERGY_GATE_MULTIPLE):
        if not math.isfinite(multiple) or not multiple > 0:
            raise ValueError("positive finite declared gate multiple required")
        self.target = float(target)
        self.energy_rows = list(energy_rows)
        self.origins = np.asarray(origins, dtype=float)
        self.multiple = float(multiple)

    def gates(self, x, jac):
        jac = np.asarray(jac, dtype=float)
        physical = np.asarray(x, dtype=float) + self.origins
        floor = np.abs(jac[self.energy_rows]) @ np.spacing(np.abs(physical))
        gates = np.full(jac.shape[0], self.target)
        gates[self.energy_rows] = np.maximum(self.target, self.multiple * floor)
        return gates, floor

    def candidate(self, residual):
        fixed = np.ones(np.asarray(residual).size, dtype=bool)
        fixed[self.energy_rows] = False
        return bool(np.all(np.abs(np.asarray(residual)[fixed]) <= self.target))

    def accepts(self, x, residual, jac):
        gates, _ = self.gates(x, jac)
        return bool(np.all(np.abs(np.asarray(residual)) <= gates))


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
    deposition = 0.0
    if isinstance(boundary, pd_boundary.KernelEnvelopeBoundary):
        # Kernel-envelope pre-desolventizer items (armed only): the jacket's conduction
        # in place of the sensible gas term; the dew-line condensate onto the film.
        q_out = pd_boundary.surface_gas_heat_w_m2(boundary, q_out, ts)
        deposition = boundary.dew_deposition_kg_m2_s
        if deposition != 0.0 and liquid == 0.0 and old_loading == 0.0:
            raise pd_boundary.DewLineCondensateOnFilmFreeSurface(
                "the dew-line projection books condensate onto the external film, but the "
                f"surface carries none ({deposition!r} kg/m2/s); a film's re-appearance on a "
                "film-free surface is not carried"
            )
    external_energy = math.fsum(
        (
            q_out,
            components.water_kg_m2_s * thermo.water_gas_partial_enthalpy_j_mol / core.wa.M,
            components.hexane_kg_m2_s * thermo.hexane_gas_partial_enthalpy_j_mol / core.hx.M,
        )
    )
    external_water = components.water_kg_m2_s
    if deposition != 0.0:
        external_water = components.water_kg_m2_s - deposition
        external_energy = math.fsum(
            (external_energy, -deposition * boundary.dew_water_vapour_enthalpy_j_kg)
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
        external_water,
        components.hexane_kg_m2_s,
        external_energy,
    )
    audit = WaterFilmNodeAudit(
        ts,
        ys,
        liquid,
        external_water,
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
            **(
                {}
                if "surface_liquid_transfer_fraction" not in holder
                else {
                    "surface_liquid_transfer_fraction": holder[
                        "surface_liquid_transfer_fraction"
                    ]
                }
            ),
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
    search_tolerance=None,
    experimental_interface_area_chart_audit=None,
    film_depletion_event=False,
    scaled_condition_cap=None,
    surface_temperature_upper_k=None,
    target_config_factory=None,
    envelope_audit=None,
    stalled_front=False,
    residual_budget_per_unknown=False,
):
    # Counterfactual chart experiment only. The F9 frozen proxy definition
    # remains the default and requires owner approval before replacement.
    if experimental_interface_area_chart_audit is not None and (
        not isinstance(experimental_interface_area_chart_audit, dict) or target_z is None
    ):
        raise ValueError(
            "experimental area chart requires a positive-core target and audit dictionary"
        )
    # Kernel-envelope options (particle_kernel_envelope; all default off, the former
    # path executes unchanged): film_depletion_event replaces the film unknown by the
    # event duration with the film fixed at exactly zero and the front free;
    # scaled_condition_cap gates the searched (column-scaled physical) matrix;
    # surface_temperature_upper_k replaces the bath as the surface bound;
    # target_config_factory(duration) gives the native pressure-only target config.
    if film_depletion_event and (
        target_z is not None
        or gas_only_surface
        or not physical_solver
        or boundary_factory is None
        or not old_loading > 0.0
    ):
        raise ValueError(
            "a front-stage film depletion event needs a positive incoming film, a free "
            "front, the physical solver and a duration-dependent boundary"
        )
    if (
        isinstance(boundary, pd_boundary.KernelEnvelopeBoundary)
        and boundary.dew_deposition_kg_m2_s != 0.0
        and old_loading == 0.0
    ):
        raise pd_boundary.DewLineCondensateOnFilmFreeSurface(
            "the dew-line projection books condensate onto the external film, but the "
            f"surface carries none ({boundary.dew_deposition_kg_m2_s!r} kg/m2/s at "
            f"bath y {boundary.dew_bath_y_hexane!r}); a film's re-appearance on a film-free "
            "surface is not carried"
        )
    if surface_temperature_upper_k is not None and not math.isfinite(surface_temperature_upper_k):
        raise ValueError("declared surface temperature bound must be finite")

    def target_kwargs(duration):
        if target_config_factory is None:
            return {}
        target = target_config_factory(duration)
        return {} if target is None else {"target_config": target}

    surface_upper = (
        boundary.temperature_k
        if surface_temperature_upper_k is None
        else surface_temperature_upper_k
    )
    if surface_temperature_upper_k is None and envelope.armed(envelope.BIRTH_SURFACE_BOUND):
        surface_upper = envelope.surface_temperature_upper_k(
            boundary.temperature_k,
            (*before.transport.wet_temperatures_k, *before.transport.dry_temperatures_k),
        )
    root_target = before.controls.nonlinear_residual_tolerance
    # Kernel-envelope FRONT_STAGE_NOISE_AWARE_ENERGY_GATE (armed only; Class B, not adopted): the
    # physical search's energy rows accepted at max(target, 2 x floor); off, the former gates.
    front_gate = envelope.armed(envelope.FRONT_STAGE_NOISE_AWARE_ENERGY_GATE)
    if search_tolerance is not None:
        if not math.isfinite(search_tolerance) or not 0 < search_tolerance <= root_target:
            raise ValueError(
                "search tolerance must be positive and no looser than native acceptance"
            )
        root_target = search_tolerance
    chart = core.ci._coordinate_chart(before)
    x_native = core.ci._encode_candidate_for_chart(
        seed.candidate, chart, seed.stefan_flux_scales_mol_m2_s
    )
    tc = seed.candidate.dry_temperatures_k[-1]
    seed_target = target_kwargs(dt).get("target_config")
    chart_pressure = (
        before.transport.config.dry.pressure_pa
        if seed_target is None
        else seed_target.dry.pressure_pa
    )
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
                ts, chart_pressure, pore
            ).y_hexane_bounds
            holder["y_hexane"] = lo + (hi - lo) * q
            holder["x"] = np.array([ts, flux_or_beta, 0.0])

        configure_gas_surface(surface_seed)
    if film_depletion_event:
        if surface_seed is None:
            raise ValueError("a film depletion event requires a declared surface seed")
        lo_d, hi_d = duration_bounds
        if not lo_d < dt < hi_d:
            raise ValueError("depletion event duration seed outside its declared bounds")
        holder["x"] = np.array([float(surface_seed[0]), float(surface_seed[1]), 0.0])
        holder["dt"] = dt
        # The exact-zero film event is the instant the D9 film-to-shell drain
        # min(film, capacity) changes branch: the film supplies exactly the
        # opened capacity (transfer share 1) and is gone. Past it the drain
        # takes whatever film remains and the zero is no longer an event of
        # the duration (the degenerate full-drain branch). Without a surface
        # transfer law the share is not read.
        holder["surface_liquid_transfer_fraction"] = 1.0
    fixed_front_coordinate = float(x_native[-2])
    if target_z is not None:
        if (
            seed.candidate.front_z != target_z
            or not 0 < target_z < before.transport.geometry.front.z
        ):
            raise ValueError("positive target must match seed and recede from accepted state")
        x_native[-2] = core.ci._encode_open(dt, duration_bounds)
    # Kernel-envelope stalled front (default off): the native interior-stall
    # construction (cut_interior_stall.py, items 1-5, ratified 2026-08-21) at the harness
    # film node. The front is a literal at the accepted position (never encoded, zero
    # sweep); its search slot leaves the unknowns and the Rankine-Hugoniot water row
    # (at q = 0 the interface water-flux continuity) leaves the searched rows; the
    # hexane row (zero net hexane release: no driving force for the front) and the
    # energy row stay solved; afterwards every row, the water row included, is judged
    # at the unchanged acceptance gate.
    stall_slot = stall_row = stall_z = None
    # Kernel-envelope STALL_FILM_DEPLETION_EVENT (armed only, with STALLED_FRONT): the
    # film-depletion event of a held front, the stalled construction above with the film
    # unknown replaced by the event duration exactly as film_depletion_event does it (the film
    # fixed at exactly zero, the drain on its capacity branch, transfer share 1); the water
    # jump row stays the q = 0 consistency row judged at the unchanged gate.
    stalled_event = (
        stalled_front
        and film_depletion_event
        and envelope.armed(envelope.STALL_FILM_DEPLETION_EVENT)
    )
    if stalled_front:
        if (
            target_z is not None
            or (film_depletion_event and not stalled_event)
            or not physical_solver
        ):
            raise ValueError(
                "a stalled-front step needs a free duration, no depletion event and the "
                "physical solver"
            )
        stall_slot = x_native.size - 2
        stall_row = (
            2 * before.transport.layout.wet_piece_count
            + 3 * before.transport.layout.dry_piece_count
        )
        stall_z = before.transport.geometry.front.z

    def expand(values):
        values = np.asarray(values, dtype=float)
        return np.insert(values, stall_slot, fixed_front_coordinate)

    def stall_diagnostic_chart():
        # the native interior stall's own accept chart (cut_interior_stall.py
        # _stall_diagnostic_chart): the pinned front judged against its host cell's
        # geometric span, never a solve chart (the stall has no front coordinate)
        geometry = before.transport.geometry
        grid = geometry.master_grid
        cut_index = before.transport.layout.cut_cell_index
        inner_z = (grid.faces[cut_index] / grid.R) ** 3
        outer_z = (grid.faces[cut_index + 1] / grid.R) ** 3
        if not inner_z < stall_z < outer_z:
            raise RuntimeError("the pinned front is not strictly inside its own host cell")
        return replace(chart, front_z=(inner_z, outer_z))

    def stall_stamp(assembly, values):
        if stalled_event:
            # the held-front event: its located duration and the consistency row in its own
            # units (mol/s of water through the pinned interface, unbalanced at q = 0)
            return {
                **stall_stamp_common(assembly, values),
                "stalled_film_depletion_event": True,
                "event_duration_s": holder["dt"],
                "consistency_row_mol_s": assembly.residuals.rh_water_mol_s,
                "film_after_kg_m2": assembly.surface_film_audit.liquid_water_kg_m2,
            }
        return stall_stamp_common(assembly, values)

    def stall_stamp_common(assembly, values):
        return {
            "front_z": stall_z,
            "swept_volume_rate_m3_s": assembly.swept_geometry.interface_swept_volume_rate_m3_s,
            "interface_temperature_k": assembly.interface.temperature_k,
            "interface_hexane_flux_mol_m2_s": (
                assembly.dry_face_fluxes[0].component.conserved_hexane_flux_mol_m2_s
            ),
            "consistency_row": "rh_water at q = 0 (interface water-flux continuity)",
            "consistency_scaled_residual": float(values[stall_row]),
            "jump_rows_scaled": [float(v) for v in values[stall_row : stall_row + 3]],
            "gate": before.controls.nonlinear_residual_tolerance,
        }

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
        if stall_slot is not None and encoded.size == x_native.size - 1:
            encoded = expand(encoded)
        if target_z is not None:
            holder["dt"] = core.ci._decode_open(encoded[-2], duration_bounds)
            encoded[-2] = fixed_front_coordinate
        c = core.ci._decode_candidate_for_chart(
            encoded, before.transport.layout, chart, seed.stefan_flux_scales_mol_m2_s
        )
        if stall_z is not None:
            return replace(c, front_z=stall_z)
        return c if target_z is None else replace(c, front_z=target_z)

    work = 0
    with (
        study.installed_boundary_adapter(numerical_cap_seam=True),
        nullcontext() if energy_context is None else energy_context(holder),
    ):
        with node_adapter(holder, lambda: holder["dt"], old_loading, old_temperature):
            initial = core.cut.assemble_backward_euler(
                before.transport,
                seed.candidate if stall_z is None else replace(seed.candidate, front_z=stall_z),
                dt,
                boundary,
                enforce_reduced_film_thresholds=False,
                **target_kwargs(dt),
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

            row_multiplier = np.ones(before.transport.layout.unknown_count + 3)
            if experimental_interface_area_chart_audit is not None:
                if not 0.0 < target_z < 1.0:
                    raise ValueError(
                        "per-area jump chart requires a finite positive interface area"
                    )
                first_rh_row = (
                    2 * before.transport.layout.wet_piece_count
                    + 3 * before.transport.layout.dry_piece_count
                )
                row_multiplier[first_rh_row : first_rh_row + 3] = target_z ** (-2.0 / 3.0)

            def residual(x):
                nonlocal work, active_boundary
                if film_depletion_event:
                    holder["dt"] = float(x[-1])
                candidate = decode_candidate(x[:-3])
                if boundary_factory is not None:
                    active_boundary = boundary_factory(holder["dt"])
                if gas_only_surface:
                    configure_gas_surface(x[-3:])
                elif film_depletion_event:
                    holder["x"] = np.array([x[-3], x[-2], 0.0])
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
                    **target_kwargs(holder["dt"]),
                )
                work += 1
                a = assembly.surface_film_audit
                native = (
                    np.asarray(core.ci._datum_covariant_residual_vector(assembly.residuals))
                    / scale_vector
                )
                rows = row_multiplier * np.r_[
                    native,
                    a.water_node_residual_kg_m2_s / mass_scale,
                    a.hexane_node_residual_kg_m2_s / mass_scale,
                    a.energy_node_residual_w_m2 / heat_scale,
                ]
                if stall_row is None:
                    return rows
                holder["stall_full_rows"] = rows
                return np.delete(rows, stall_row)

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
            if film_depletion_event:
                surface_initial[2] = dt
            liquid_hi = 1.0 if gas_only_surface else max(0.05, 4 * old_loading)
            coordinate_limit = boundary.accepted_beta_limit if beta_decoder is not None else 1.0
            if physical_solver:
                layout = before.transport.layout
                nw, nd, nf = layout.wet_piece_count, layout.dry_piece_count, layout.dry_face_count
                pressure = chart_pressure

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
                        (tl, min(tu, surface_upper)),
                        (-coordinate_limit, coordinate_limit),
                        duration_bounds if film_depletion_event else (0.0, liquid_hi),
                    ]
                )
                if stall_slot is not None:
                    # The stalled front's slot leaves the physical chart as well; the
                    # placeholder front handed to the full encoder is the seed's own
                    # (encodable) prediction and leaves with its slot.
                    full_decode, full_encode = decode, encode
                    placeholder_z = seed.candidate.front_z

                    def decode(v):
                        return np.delete(full_decode(expand(v)), stall_slot)

                    def encode(v):
                        trial = np.insert(np.asarray(v, dtype=float), stall_slot, placeholder_z)
                        return np.delete(full_encode(trial), stall_slot)

                    pairs = pairs[:stall_slot] + pairs[stall_slot + 1 :]

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
                            "gtol": root_target,
                        },
                    )

            solve_kwargs = dict(
                bounds=(
                    np.r_[lb, tl, -coordinate_limit, 0.0],
                    np.r_[ub, min(tu, surface_upper), coordinate_limit, liquid_hi],
                ),
                x_scale="jac",
                jac="3-point",
                diff_step=1e-4,
                max_nfev=400,
                xtol=1e-11,
                ftol=1e-11,
                gtol=1e-11,
            )
            if stalled_event:
                # the event duration's native box is its declared duration bounds, as its
                # physical chart already states (the liquid-loading box does not apply to a
                # duration; a held front's film runs out over tenths of a second)
                lo_box, hi_box = (np.asarray(b, dtype=float).copy() for b in solve_kwargs["bounds"])
                lo_box[-1], hi_box[-1] = duration_bounds
                solve_kwargs["bounds"] = (lo_box, hi_box)
            restart_audit = []
            start = np.r_[x_native, surface_initial]
            if stall_slot is not None:
                start = np.delete(start, stall_slot)
                solve_kwargs["bounds"] = tuple(
                    np.delete(np.asarray(b, dtype=float), stall_slot)
                    for b in solve_kwargs["bounds"]
                )
            total_calls = 0
            # Ruling packet 6 item 1 (DTDC_KERNEL_ENVELOPE_BUDGETS_PER_UNKNOWN, armed by the caller;
            # GT_PS2_PAPER1_BUDGETS_PER_UNKNOWN_BUILD_2026-10-24.md): the physical search's total
            # max(residual_budget, 12 (n + 1)) and its first two cycles max(present, 4 (n + 1) + 1),
            # n the search unknowns. Off: the caller's residual_budget and the present split.
            search_budget = residual_budget
            if residual_budget_per_unknown and physical_solver:
                search_budget = envelope.front_search_total_budget(residual_budget, start.size)
            # A restart is another search for the SAME backward-Euler root.
            # It imports no accepted state or elapsed physical time. Native
            # normalization is rebuilt from this interval's trial assembly.
            for cycle in range(3 if physical_solver else 1):
                remaining = search_budget - total_calls
                if physical_solver:
                    split = max(200, residual_budget // 3)
                    if residual_budget_per_unknown:
                        split = envelope.front_search_cycle_budget(split, start.size)
                    cycle_budget = min(remaining, split) if cycle < 2 else remaining

                    def cycle_solver(fun, x0, **kwargs):
                        if front_gate:
                            # FRONT_STAGE_NOISE_AWARE_ENERGY_GATE (armed only): the per-row
                            # acceptance clause, its origins those physical_chart_solve takes
                            gate_initial = np.asarray(decode(x0), dtype=float)
                            gate_indices = [
                                *range(nw),
                                *range(2 * nw, 2 * nw + nd),
                                len(pairs) - 4,
                                len(pairs) - 3,
                            ]
                            gate_origin = np.zeros(gate_initial.size)
                            gate_origin[gate_indices] = gate_initial[gate_indices]
                            gate_rows = front_energy_rows(
                                nw, nd, before.transport.layout.unknown_count + 3, stall_row
                            )
                            holder["front_gate_rows"] = gate_rows
                            kwargs["acceptance"] = FrontNoiseAwareEnergyRowAcceptance(
                                root_target, gate_rows, gate_origin
                            )
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
                                "gtol": root_target,
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
                    if cycle == 2 or total_calls >= search_budget:
                        debug["same_interval_restarts"] = restart_audit
                        raise
                if not physical_solver or cycle == 2 or total_calls >= search_budget:
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
                    **target_kwargs(holder["dt"]),
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
            if residual_budget_per_unknown and physical_solver:
                solution.kernel_envelope_search_budget = {
                    "unknowns": int(start.size),
                    "total": search_budget,
                    "first_two_cycles": envelope.front_search_cycle_budget(
                        max(200, residual_budget // 3), start.size
                    ),
                }
            solution.nfev = total_calls
            search_values = residual(solution.x)
            if stall_row is not None:
                # every row, the water consistency row included, is judged at the gate
                search_values = holder["stall_full_rows"]
            values = search_values / row_multiplier
            maximum = float(np.max(np.abs(values)))
            candidate = decode_candidate(solution.x[:-3])
            assembly = core.cut.assemble_backward_euler(
                before.transport,
                candidate,
                holder["dt"],
                active_boundary,
                enforce_reduced_film_thresholds=False,
                **target_kwargs(holder["dt"]),
            )
            gate_stamp = None
            if front_gate:
                # FRONT_STAGE_NOISE_AWARE_ENERGY_GATE (armed only): the post-search check at the
                # per-row gates, the energy rows' lower bound the native tolerance, every other
                # row at the unchanged native tolerance.
                check_gates = np.full(values.size, before.controls.nonlinear_residual_tolerance)
                if physical_solver and "front_gate_rows" in holder and solution.success:
                    full_rows = front_energy_rows(
                        nw, nd, before.transport.layout.unknown_count + 3, None
                    )
                    search_rows = holder["front_gate_rows"]
                    physical = np.asarray(getattr(solution, "physical_coordinates", solution.x))
                    jac_rows = np.asarray(solution.jac, dtype=float)
                    floor = np.abs(jac_rows[search_rows]) @ np.spacing(np.abs(physical))
                    value_rows = (
                        full_rows if stall_row is None else [r for r in full_rows if r != stall_row]
                    )
                    gate_values = np.maximum(
                        before.controls.nonlinear_residual_tolerance,
                        FRONT_STAGE_ENERGY_GATE_MULTIPLE * floor,
                    )
                    check_gates[value_rows] = gate_values
                    gate_stamp = {
                        "knob": envelope.PREFIX + envelope.FRONT_STAGE_NOISE_AWARE_ENERGY_GATE,
                        "declared_multiple": FRONT_STAGE_ENERGY_GATE_MULTIPLE,
                        "dt_s": float(holder["dt"]),
                        "stalled_front": stall_row is not None,
                        "film_depletion_event": bool(film_depletion_event),
                        "energy_rows": list(map(int, value_rows)),
                        "energy_row_residuals": [float(values[r]) for r in value_rows],
                        "energy_row_floors": floor.tolist(),
                        "energy_row_gates": gate_values.tolist(),
                        "absorbed_by_floor": [
                            bool(abs(values[r]) > before.controls.nonlinear_residual_tolerance)
                            for r in value_rows
                        ],
                        "non_energy_max_abs": float(np.max(np.abs(np.delete(values, value_rows)))),
                        "all_residuals": values.tolist(),
                        "target": float(root_target),
                        "native_tolerance": float(before.controls.nonlinear_residual_tolerance),
                        "front_z": float(candidate.front_z),
                        "solver_message": str(solution.message),
                        "search_calls": int(total_calls),
                        "physically_qualifying": False,
                    }
                refused = not solution.success or bool(np.any(np.abs(values) > check_gates))
            else:
                refused = (
                    not solution.success
                    or np.max(np.abs(values)) > before.controls.nonlinear_residual_tolerance
                )
            if refused:
                error = RuntimeError(
                    f"water-node study not converged: residual={np.max(np.abs(values))}; work={work}"
                    if stall_row is None
                    else "stalled-front step refused: searched rows "
                    f"{'converged' if solution.success else 'not converged'}, water consistency "
                    f"row {float(values[stall_row])!r} against the gate "
                    f"{before.controls.nonlinear_residual_tolerance!r}; work={work}"
                )
                if stall_row is not None:
                    error.stall_stamp = stall_stamp(assembly, values)
                error.study_debug = {
                    "candidate": asdict(candidate),
                    "surface": asdict(assembly.surface_film_audit),
                    "scaled_residuals": values.tolist(),
                    "work": work,
                    "same_interval_restarts": restart_audit,
                }
                if front_gate:
                    error.study_debug["front_stage_noise_aware_energy_gate"] = gate_stamp
                raise error
            if abs(assembly.surface_film_audit.beta) > boundary.accepted_beta_limit:
                raise RuntimeError("resolved water-node root exceeds explicit beta study band")
            if gate_stamp is not None and envelope_audit is not None:
                envelope_audit["front_stage_noise_aware_energy_gate"] = gate_stamp
            condition = float(np.linalg.cond(getattr(solution, "native_jac", solution.jac)))
            if experimental_interface_area_chart_audit is not None:
                native_matrix = getattr(solution, "native_jac", solution.jac)
                legacy_condition = float(np.linalg.cond(native_matrix / row_multiplier[:, None]))
                experimental_interface_area_chart_audit.update(
                    experimental_counterfactual=True,
                    owner_adoption_required=True,
                    physically_qualifying=False,
                    chart="native dimensionless gas-barrier columns; RH rows multiplied by A_outer/A_Gamma",
                    row_multipliers=row_multiplier.tolist(),
                    original_integrated_row_condition=legacy_condition,
                    original_condition_limit=before.controls.maximum_condition_proxy,
                    original_condition_gate_would_refuse=legacy_condition
                    > before.controls.maximum_condition_proxy,
                    per_area_chart_condition=condition,
                    maximum_original_scaled_residual=float(max(abs(values))),
                    maximum_per_area_scaled_residual=float(max(abs(search_values))),
                    physical_equations_unchanged=True,
                )
            gated_condition, gated_cap = condition, before.controls.maximum_condition_proxy
            if scaled_condition_cap is not None:
                gated_condition = study.numerics.searched_matrix_condition(solution.jac)
                gated_cap = scaled_condition_cap
            if envelope_audit is not None:
                envelope_audit.update(
                    film_node_certificate={
                        "gated_matrix": (
                            "native pullback (former)"
                            if scaled_condition_cap is None
                            else "column-2-norm-scaled physical film-node Jacobian"
                        ),
                        "gated_condition": gated_condition,
                        "cap": gated_cap,
                        "native_pullback_condition": condition,
                        "native_pullback_cap_reported": before.controls.maximum_condition_proxy,
                        "scaled_physical_condition": study.numerics.searched_matrix_condition(
                            solution.jac
                        ),
                    },
                    surface_temperature_upper_k=surface_upper,
                    film_depletion_event=film_depletion_event,
                    event_duration_s=holder["dt"] if film_depletion_event else None,
                    target_pressure_pa=assembly.target_config.dry.pressure_pa,
                )
                if stall_z is not None:
                    envelope_audit["stalled_front"] = stall_stamp(assembly, values)
            if gated_condition > gated_cap or (
                scaled_condition_cap is not None and not math.isfinite(gated_condition)
            ):
                error = RuntimeError(
                    f"water-node Jacobian condition exceeds native limit: {condition}"
                    if scaled_condition_cap is None
                    else "water-node condition not certified on the searched "
                    f"(column-2-norm-scaled) matrix: {gated_condition} "
                    f"(native pullback {condition})"
                )
                error.optimizer_debug = {
                    "native_jacobian": getattr(solution, "native_jac", solution.jac).tolist(),
                    "physical_jacobian": solution.jac.tolist(),
                    "physical_coordinates": getattr(
                        solution, "physical_coordinates", solution.x
                    ).tolist(),
                    "native_coordinates": solution.x.tolist(),
                    "scaled_residual": values.tolist(),
                    "native_row_scales": list(scales.vector),
                    "surface_row_scales": [mass_scale, mass_scale, heat_scale],
                    "condition_proxy": condition,
                    "original_condition_limit": before.controls.maximum_condition_proxy,
                    "physical_acceptance_override": False,
                }
                raise error
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
                chart if stall_z is None else stall_diagnostic_chart(),
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
            if stall_z is not None and (
                step.after.transport.geometry.front.z != stall_z
                or assembly.swept_geometry.interface_swept_volume_rate_m3_s != 0.0
            ):
                raise RuntimeError("stalled-front step moved the front")
            return step
