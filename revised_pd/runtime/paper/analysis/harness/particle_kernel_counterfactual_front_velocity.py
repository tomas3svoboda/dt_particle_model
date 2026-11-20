"""COUNTERFACTUAL EVALUATION ONLY (evidence 170). Not a build path; physically_qualifying false.

Harness copy for the STALL_FILM_DEPLETION_EVENT knob's materiality row (evidence
stall_depletion_171): evidence pd_advance_demand_170's scripts/counterfactual_front_velocity.py
(sha256 fcac0126295da222510b5bd450e5959e8ad66c4b4a4c52c627596bed873b6d8b), formatted by Black
and this paragraph added; its code is the evidence-170 construction unchanged (the abstract
syntax trees equal but for this docstring, evidence 171 records/cf_module_identity.json). It is
called only by particle_kernel_stall_materiality, in a scratch overlay of the film node's
advance, after a stalled step is accepted; nothing it returns is adopted into any history.

The front velocity as a free, SIGNED unknown at a pinned front: the jump rows evaluated around
the frozen PHY-014 latch, never with it removed.

What this is. A verbatim copy of the harness's
``core2_water_film_node.advance_water_node`` (the file is not edited; this copy lives in the
evidence folder) in which the ``stalled_front=True`` branch is re-purposed:

* the front is a literal at the accepted position exactly as in the native interior-stall
  construction (zero swept volume, the geometry and every storage term frozen; PHY-014's
  vetoes at cut_geometry/front are never reached because nothing moves);
* the vacated search slot carries v, the front's normal velocity in micrometres per second,
  POSITIVE = outward = ADVANCE (the forbidden direction), unbounded in both signs;
* every Rankine-Hugoniot row is kept (none removed) and receives the jump term a moving front
  would add, with q = dV_wet/dt = A_Gamma * v (the native rows' own signed expansion,
  cut_transport._residual_blocks):
      rh_water  -= q (C_w,dry - C_w,donor)      [O9a donor, cut_integrator helper]
      rh_hexane -= q (C_h,dry - C_h,wet)        [C_h,wet = rho_dm_p X_hist / M_h]
      rh_energy -= q (e_dry  - e_donor)         [O10a donor, cut_integrator helper]
  the datum-covariant row basis then applied unchanged;
* no step is accepted: the solve's root (v, T_Gamma, the interior and the surface node) is
  returned as evidence with its residual by block. The regional rows keep their frozen-front
  form (the leading-order ALE limit: a moving front changes them only through the swept
  storage, which is excluded by construction - the same declared scope limit as the native
  cut_front_demand criterion). The material the jump terms move across the front is
  therefore NOT booked into a swept wet volume; its integral is reported, never hidden.
"""

from contextlib import nullcontext
from dataclasses import asdict, replace
import math

import numpy as np
from scipy.optimize import least_squares

import finite_film_core2_study as study
import particle_kernel_envelope as envelope
import particle_kernel_pd_boundary as pd_boundary
from core2_water_film_node import node_adapter

core = study.core
VELOCITY_UNIT_M_S = 1e-6
LAST = {}


def jump_coefficients(assembly):
    """(A_Gamma, dC_water, dC_hexane, de) of the native RH rows at this assembly."""
    before = assembly.before
    interface = assembly.interface
    area = 4.0 * math.pi * assembly.candidate_geometry.front.radius_m**2
    donor_w = core.ci._assembly_sweep_donor_concentration(before)
    donor_e = core.ci._assembly_sweep_donor_energy_density(before)
    wet_h = before.config.wet.wet.rho_dm_p * interface.historical_hexane_loading / core.hx.M
    return {
        "area_m2": area,
        "dry_water_mol_m3": interface.dry.total_water_concentration_mol_m3,
        "donor_water_mol_m3": donor_w,
        "dry_hexane_mol_m3": interface.dry.total_hexane_concentration_mol_m3,
        "wet_hexane_mol_m3": wet_h,
        "dry_energy_j_m3": interface.dry.energy_density_j_m3,
        "donor_energy_j_m3": donor_e,
        "jump_water_mol_m3": math.fsum((interface.dry.total_water_concentration_mol_m3, -donor_w)),
        "jump_hexane_mol_m3": math.fsum((interface.dry.total_hexane_concentration_mol_m3, -wet_h)),
        "jump_energy_j_m3": math.fsum((interface.dry.energy_density_j_m3, -donor_e)),
    }


def moved_residuals(assembly, velocity):
    """The native residual blocks with the moving-front jump terms at q = A v."""
    j = jump_coefficients(assembly)
    q = j["area_m2"] * velocity * VELOCITY_UNIT_M_S
    r = assembly.residuals
    return (
        replace(
            r,
            rh_water_mol_s=math.fsum((r.rh_water_mol_s, -q * j["jump_water_mol_m3"])),
            rh_hexane_mol_s=math.fsum((r.rh_hexane_mol_s, -q * j["jump_hexane_mol_m3"])),
            rh_energy_w=math.fsum((r.rh_energy_w, -q * j["jump_energy_j_m3"])),
        ),
        q,
        j,
    )


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
            lo, hi = core.cp.gas_only_composition_interval(ts, chart_pressure, pore).y_hexane_bounds
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
    if stalled_front:
        if target_z is not None or film_depletion_event or not physical_solver:
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

    # COUNTERFACTUAL (evidence 170): the slot is kept and carries the signed velocity.
    holder["v"] = float(LAST.get("velocity_seed_um_s", 0.0))

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
        if stall_slot is not None:
            holder["v"] = float(encoded[stall_slot])
            encoded[stall_slot] = fixed_front_coordinate
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
                blocks = assembly.residuals
                if stall_slot is not None:
                    blocks, _, _ = moved_residuals(assembly, holder["v"])
                native = np.asarray(core.ci._datum_covariant_residual_vector(blocks)) / scale_vector
                rows = (
                    row_multiplier
                    * np.r_[
                        native,
                        a.water_node_residual_kg_m2_s / mass_scale,
                        a.hexane_node_residual_kg_m2_s / mass_scale,
                        a.energy_node_residual_w_m2 / heat_scale,
                    ]
                )
                return rows

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
                    # COUNTERFACTUAL: the slot is the signed velocity (um/s) in both charts
                    # (identity); the placeholder front handed to the full encoder is the
                    # seed's own (encodable) prediction, the decoded front is the literal.
                    full_decode, full_encode = decode, encode
                    placeholder_z = seed.candidate.front_z

                    def decode(v):
                        v = np.asarray(v, dtype=float)
                        velocity = float(v[stall_slot])
                        trial = v.copy()
                        trial[stall_slot] = fixed_front_coordinate
                        out = np.asarray(full_decode(trial), dtype=float)
                        out[stall_slot] = velocity
                        return out

                    def encode(v):
                        v = np.asarray(v, dtype=float)
                        velocity = float(v[stall_slot])
                        trial = v.copy()
                        trial[stall_slot] = placeholder_z
                        out = np.asarray(full_encode(trial), dtype=float)
                        out[stall_slot] = velocity
                        return out

                    pairs = pairs[:stall_slot] + [(-math.inf, math.inf)] + pairs[stall_slot + 1 :]

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
            restart_audit = []
            start = np.r_[x_native, surface_initial]
            if stall_slot is not None:
                start[stall_slot] = holder["v"]
                lo_b, hi_b = (np.asarray(b, dtype=float).copy() for b in solve_kwargs["bounds"])
                lo_b[stall_slot], hi_b[stall_slot] = -1e12, 1e12
                solve_kwargs["bounds"] = (lo_b, hi_b)
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
            solution.nfev = total_calls
            # COUNTERFACTUAL RETURN (evidence 170): nothing is accepted, nothing is
            # committed; the root and its rows are evidence only.
            values = residual(solution.x) / row_multiplier
            candidate = decode_candidate(solution.x[:-3])
            assembly = core.cut.assemble_backward_euler(
                before.transport,
                candidate,
                holder["dt"],
                active_boundary,
                enforce_reduced_film_thresholds=False,
                **target_kwargs(holder["dt"]),
            )
            moved, q, jumps = moved_residuals(assembly, holder["v"])
            nw_, nd_ = (
                before.transport.layout.wet_piece_count,
                before.transport.layout.dry_piece_count,
            )
            first_rh = 2 * nw_ + 3 * nd_
            dry_face = assembly.dry_face_fluxes[0]
            wet_face = assembly.wet_face_fluxes[-1]
            evidence = {
                "counterfactual": True,
                "physically_qualifying": False,
                "latch_removed": False,
                "construction": "pinned front, signed free velocity in the RH jump terms only",
                "solver_success": bool(solution.success),
                "solver_message": str(solution.message),
                "residual_calls": int(solution.nfev),
                "max_scaled_residual": float(np.max(np.abs(values))),
                "gate": before.controls.nonlinear_residual_tolerance,
                "scaled_rh_rows": [float(x) for x in values[first_rh : first_rh + 3]],
                "scaled_node_rows": [float(x) for x in values[-3:]],
                "dt_s": holder["dt"],
                "start_time_s": before.transport.time_s,
                "velocity_um_s": holder["v"],
                "q_m3_s": q,
                "jumps": jumps,
                "interface_temperature_k": assembly.interface.temperature_k,
                "interface_y_hexane": assembly.interface.y_hexane,
                "interface_historical_hexane_loading": assembly.interface.historical_hexane_loading,
                "interface_gas_water_activity": assembly.interface.dry.water_activity,
                "interface_gas_hexane_activity": assembly.interface.dry.hexane_activity,
                "interface_retained_hexane_loading": assembly.interface.dry.retained_hexane_loading,
                "interface_retained_water_loading": assembly.interface.dry.retained_water_loading,
                "dry_interface_water_flux_mol_m2_s": dry_face.component.conserved_water_flux_mol_m2_s,
                "dry_interface_hexane_flux_mol_m2_s": dry_face.component.conserved_hexane_flux_mol_m2_s,
                "dry_interface_energy_flux_w_m2": dry_face.energy.total_energy_flux_w_m2,
                "wet_interface_water_flux_mol_m2_s": wet_face.retained_water_flux_mol_m2_s,
                "wet_interface_energy_flux_w_m2": wet_face.total_energy_flux_w_m2,
                "moved_rh_raw": [moved.rh_water_mol_s, moved.rh_hexane_mol_s, moved.rh_energy_w],
                "frozen_rh_raw": [
                    assembly.residuals.rh_water_mol_s,
                    assembly.residuals.rh_hexane_mol_s,
                    assembly.residuals.rh_energy_w,
                ],
                "candidate": asdict(candidate),
                "surface": asdict(assembly.surface_film_audit),
                "boundary": {
                    "temperature_k": active_boundary.temperature_k,
                    "y_hexane": active_boundary.y_hexane,
                    "jacket_wall_temperature_k": getattr(
                        active_boundary, "jacket_wall_temperature_k", None
                    ),
                    "jacket_contact_conductance_w_m2_k": getattr(
                        active_boundary, "jacket_contact_conductance_w_m2_k", None
                    ),
                    "dew_deposition_kg_m2_s": getattr(
                        active_boundary, "dew_deposition_kg_m2_s", 0.0
                    ),
                },
                "pressure_pa": assembly.target_config.dry.pressure_pa,
                "assembly": assembly,
            }
            LAST["evidence"] = evidence
            return evidence
