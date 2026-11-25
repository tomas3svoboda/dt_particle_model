"""Existing DT/Faner dry equations with conservative surface-phase dispatch.

The caller owns the single proof/assembly allowance for a physical advance,
including all rejected active sets and time-step halvings. This dispatcher
never changes the old physical state, a constitutive law or an acceptance gate.
"""

from dataclasses import dataclass, field, replace
import math

import experiment_core2_dt_dry as dry
import experiment_core2_radial_fast_wetting as liquid
import particle_kernel_envelope as envelope
from dtdc_simulator.core2.particle.numerical_work import StepBudgetExceeded


class DryPhaseAdvanceError(RuntimeError):
    def __init__(self, message, before, rejected_trials):
        super().__init__(message)
        self.rollback_state = before
        self.rejected_trials = tuple(rejected_trials)


@dataclass(frozen=True)
class DryPhaseAdvance:
    after: object
    surface: object
    carry: dict
    native_report: dict
    phase: str
    depletion_event: bool
    liquid_supply_audit: dict | None
    rejected_trials: tuple[dict, ...]
    # Kernel-envelope HEXANE_EXHAUSTED_REGIME (armed only): the stamped release or re-entry event of
    # this advance, each with its outward stream (W kg, H kg, U J) for the caller's stream integral.
    # Empty on every advance of the binary stage without an event (always, with the knob absent).
    regime_events: tuple[dict, ...] = field(default=())


def _absent_seed_lift(before, seed, boundary):
    """ABSENT_SEED_ABOVE_DEW's lifted seed and stamp, else None: the seed's zero-flux interface
    (the presented gas composition) lies at or past the water-saturated end of the gas-only chart
    at the seed temperature."""
    import particle_steam_condensing_film as condensing_film

    pressure = float(before.config.pressure_pa)
    gas_y = condensing_film._presented_y(boundary)
    try:
        lower = condensing_film.interface_y_lower_bound(
            float(seed.temperature_k), pressure, before.config.pore
        )
    except Exception:  # no chart at the seed temperature: the seed is kept
        return None
    if gas_y > lower * (1.0 + envelope.ABSENT_SEED_ABOVE_DEW_CHART_MARGIN):
        return None
    lifted = float(boundary.temperature_k) + envelope.ABSENT_SEED_ABOVE_DEW_OFFSET_K
    if not lifted > float(seed.temperature_k):
        return None
    return {
        "former_seed_temperature_k": float(seed.temperature_k),
        "former_seed_interface_y_lower_bound": lower,
        "chart_margin": envelope.ABSENT_SEED_ABOVE_DEW_CHART_MARGIN,
        "gas_y_hexane": gas_y,
        "seed_temperature_k": lifted,
        "seed_y_hexane": float(boundary.y_hexane),
        "offset_k": envelope.ABSENT_SEED_ABOVE_DEW_OFFSET_K,
        "order": "unchanged",
        "record": envelope.CONVERGENCE_ROOT_CAUSES_RECORD,
        "physically_qualifying": False,
    }


def _film_free_wet_set_skip(before, surface, duration, boundary_factory):
    """FILM_FREE_WET_SET_SKIP's condition and stamp, else None: the presented gas at both ends of
    the span is leaner in water than the kernel's water-saturated interface at the surface
    temperature minus the declared margin (no water condenses there)."""
    import particle_steam_condensing_film as condensing_film

    pressure = float(before.config.pressure_pa)
    colder = float(surface.temperature_k) - envelope.FILM_FREE_WET_SET_SKIP_MARGIN_K
    try:
        lower = condensing_film.interface_y_lower_bound(colder, pressure, before.config.pore)
    except Exception:  # no chart at that temperature: every set is tried
        return None
    ys = [
        condensing_film._presented_y(boundary_factory(t))
        for t in (before.time_s, before.time_s + duration)
    ]
    if not all(y >= lower for y in ys):
        return None
    return {
        "surface_temperature_k": float(surface.temperature_k),
        "margin_k": envelope.FILM_FREE_WET_SET_SKIP_MARGIN_K,
        "interface_y_lower_bound_at_margin": lower,
        "gas_y_hexane_start_end": ys,
        "sets_tried": ["absent"],
        "record": envelope.CONVERGENCE_ROOT_CAUSES_RECORD,
        "physically_qualifying": False,
    }


def advance_dry_phases(
    before,
    config,
    carry,
    surface,
    duration,
    boundary_factory,
    *,
    previous_phase,
    water_inventory_scale_kg,
    assembly_budget,
    liquid_mobility_factor=envelope.LIQUID_MOBILITY_TOTAL_FACTOR,
    contact_conductance=None,
    search_tolerance=None,
    global_component_mass_error_budgets_kg=None,
    use_native_conduction_predictor=False,
    conservative_liquid_elimination=False,
    film_depletion_event_first=False,
    predicted_depletion_duration_s=None,
    stop_positive_film_at_nonpositive_liquid=False,
    conservative_event_elimination=False,
    conservative_absent_elimination=False,
    previous_accepted_profile=None,
    target_pressure_pa=None,
    _hexane_regime_inner=False,
    _steam_quench_inner=False,
    _condensing_film_inner=False,
):
    """Try the repository's liquid/contact/absent active sets on one old state.

    Positive liquid may end only at a solved exact-zero-inventory event.
    Condensate-limited contact must close its native supply equation. A spent
    global allowance is fatal and never changes phase or retries a root.

    Declared search options (all off by default; off is the former code path):
    film_depletion_event_first tries the exact-zero event first, seeded with
    the caller's predicted duration, when that prediction lies inside the
    attempted interval; stop_positive_film_at_nonpositive_liquid ends the
    eliminated positive-film search at its first nonpositive-film trial;
    conservative_event/absent_elimination fix the component fluxes of those
    active sets algebraically from the cell inventories. No equation, gate,
    tolerance or allowance changes. previous_accepted_profile, a tuple
    (time_s, temperatures_k, y_hexane) of the accepted state before this one,
    seeds an event's cells by linear extrapolation to its seed duration.

    contact_conductance (default None): the liquid contact conductance K_l of the adopted
    film-covered surface law, read at call time from the one declared constant
    particle_kernel_envelope.K_L_MOL_M2_S when not given.
    """
    if not math.isfinite(duration) or duration <= 0:
        raise ValueError("positive finite dry advance required")
    if not _condensing_film_inner and envelope.armed(envelope.STEAM_CONDENSING_FILM):
        # Kernel-envelope STEAM_CONDENSING_FILM (owner ruling 2026-11-04, armed only): the true-gas
        # condensing-film branch of a binary dry state in hexane-lean superheated steam; None leaves
        # this advance to the dispatch below, unchanged.
        import particle_steam_condensing_film as condensing_film

        def former_dispatch():
            return advance_dry_phases(
                before,
                config,
                carry,
                surface,
                duration,
                boundary_factory,
                previous_phase=previous_phase,
                water_inventory_scale_kg=water_inventory_scale_kg,
                assembly_budget=assembly_budget,
                liquid_mobility_factor=liquid_mobility_factor,
                contact_conductance=contact_conductance,
                search_tolerance=search_tolerance,
                global_component_mass_error_budgets_kg=global_component_mass_error_budgets_kg,
                use_native_conduction_predictor=use_native_conduction_predictor,
                conservative_liquid_elimination=conservative_liquid_elimination,
                film_depletion_event_first=film_depletion_event_first,
                predicted_depletion_duration_s=predicted_depletion_duration_s,
                stop_positive_film_at_nonpositive_liquid=stop_positive_film_at_nonpositive_liquid,
                conservative_event_elimination=conservative_event_elimination,
                conservative_absent_elimination=conservative_absent_elimination,
                previous_accepted_profile=previous_accepted_profile,
                target_pressure_pa=target_pressure_pa,
                _hexane_regime_inner=_hexane_regime_inner,
                _steam_quench_inner=_steam_quench_inner,
                _condensing_film_inner=True,
            )

        condensed = condensing_film.route(
            before,
            config,
            carry,
            surface,
            duration,
            boundary_factory,
            previous_phase=previous_phase,
            dispatch=former_dispatch,
            advance_type=DryPhaseAdvance,
            error_type=DryPhaseAdvanceError,
            options={
                "assembly_budget": assembly_budget,
                "search_tolerance": search_tolerance,
                "global_component_mass_error_budgets_kg": global_component_mass_error_budgets_kg,
                "target_pressure_pa": target_pressure_pa,
                "liquid_mobility_factor": liquid_mobility_factor,
                "contact_conductance": contact_conductance,
                "predicted_depletion_duration_s": predicted_depletion_duration_s,
            },
        )
        if condensed is not None:
            return condensed
    if not _steam_quench_inner and envelope.armed(envelope.STEAM_QUENCH_ARRIVAL):
        # Kernel-envelope STEAM_QUENCH_ARRIVAL (2026-11-03, armed only): the stamped condensation
        # event of a sub-saturation particle in hexane-free superheated steam; None leaves this
        # advance to the dispatch below, unchanged.
        import particle_steam_quench_arrival as steam_quench

        def quenched_stage(state, next_carry, seed_surface, phase):
            return advance_dry_phases(
                state,
                config,
                next_carry,
                seed_surface,
                duration,
                boundary_factory,
                previous_phase=phase,
                water_inventory_scale_kg=water_inventory_scale_kg,
                assembly_budget=assembly_budget,
                liquid_mobility_factor=liquid_mobility_factor,
                contact_conductance=contact_conductance,
                search_tolerance=search_tolerance,
                global_component_mass_error_budgets_kg=global_component_mass_error_budgets_kg,
                use_native_conduction_predictor=use_native_conduction_predictor,
                conservative_liquid_elimination=conservative_liquid_elimination,
                film_depletion_event_first=film_depletion_event_first,
                predicted_depletion_duration_s=None,
                stop_positive_film_at_nonpositive_liquid=stop_positive_film_at_nonpositive_liquid,
                conservative_event_elimination=conservative_event_elimination,
                conservative_absent_elimination=conservative_absent_elimination,
                previous_accepted_profile=None,
                target_pressure_pa=target_pressure_pa,
                _hexane_regime_inner=_hexane_regime_inner,
                _steam_quench_inner=True,
                _condensing_film_inner=_condensing_film_inner,
            )

        quenched = steam_quench.route(
            before, carry, surface, duration, boundary_factory, binary=quenched_stage
        )
        if quenched is not None:
            return quenched
    if not _hexane_regime_inner and envelope.armed(envelope.HEXANE_EXHAUSTED_REGIME):
        # Kernel-envelope HEXANE_EXHAUSTED_REGIME (owner ruling 2026-11-01, armed only): the
        # n-hexane-exhausted regime's entry, vertex and re-entry; None leaves this advance to the
        # binary active sets below, unchanged.
        import particle_hexane_exhausted_regime as hexane_regime

        def binary_stage(state, seed_surface, phase):
            return advance_dry_phases(
                state,
                config,
                carry,
                seed_surface,
                duration,
                boundary_factory,
                previous_phase=phase,
                water_inventory_scale_kg=water_inventory_scale_kg,
                assembly_budget=assembly_budget,
                liquid_mobility_factor=liquid_mobility_factor,
                contact_conductance=contact_conductance,
                search_tolerance=search_tolerance,
                global_component_mass_error_budgets_kg=global_component_mass_error_budgets_kg,
                use_native_conduction_predictor=use_native_conduction_predictor,
                conservative_liquid_elimination=conservative_liquid_elimination,
                film_depletion_event_first=film_depletion_event_first,
                predicted_depletion_duration_s=predicted_depletion_duration_s,
                stop_positive_film_at_nonpositive_liquid=stop_positive_film_at_nonpositive_liquid,
                conservative_event_elimination=conservative_event_elimination,
                conservative_absent_elimination=conservative_absent_elimination,
                previous_accepted_profile=None,
                target_pressure_pa=target_pressure_pa,
                _hexane_regime_inner=True,
                _condensing_film_inner=_condensing_film_inner,
            )

        routed = hexane_regime.route(
            before,
            config,
            carry,
            surface,
            duration,
            boundary_factory,
            advance_type=DryPhaseAdvance,
            binary=binary_stage,
            assembly_budget=assembly_budget,
            target_pressure_pa=target_pressure_pa,
            search_tolerance=search_tolerance,
            global_component_mass_error_budgets_kg=global_component_mass_error_budgets_kg,
        )
        if routed is not None:
            return routed
    if contact_conductance is None:
        contact_conductance = envelope.surface_contact_conductance()
    # Kernel-envelope pressure schedule (default None: nothing passed): the
    # pressure of the attempted interval's end, for every active set.
    pressure_option = (
        {} if target_pressure_pa is None else {"target_pressure_pa": target_pressure_pa}
    )
    # LIQUID_MOBILITY_RULE (study knob, default unset: nothing passed): the interior liquid-path
    # rule, every active set.
    contact_option = {}
    mobility_rule = envelope.liquid_mobility_rule()
    if mobility_rule is not None:
        contact_option = {**contact_option, "liquid_mobility_rule": mobility_rule}
    if mobility_rule == envelope.LIQUID_MOBILITY_LOCAL_STATE_LAGGED:
        # local_state_lagged reads the rule on this advance's accepted start state.
        contact_option = {
            **contact_option,
            "liquid_mobility_start_state": (before.temperatures_k, before.y_hexane, before.config),
        }
    if not math.isfinite(water_inventory_scale_kg) or water_inventory_scale_kg <= 0:
        raise ValueError("positive finite water ledger scale required")
    if not isinstance(assembly_budget, dict):
        raise TypeError("one shared physical-advance assembly allowance is required")
    # Ruling GT_PS2_OWNER_RULING_PAPER1_GATE_UNITS_2026-10-17 point 1: a search
    # target other than the declared gate is refused here, before any phase
    # handler could absorb the typed refusal as a rejected active set.
    if search_tolerance is not None:
        dry.declared_dry_search_target(before.config.nonlinear_residual_tolerance, search_tolerance)
    boundary = boundary_factory(before.time_s + duration)
    present = carry["external_water_loading_kg_m2"] > 0
    absent_guess = surface
    following_wet_contact = previous_phase in ("liquid", "contact")
    film_free_seed = None
    gas_fallback_guess = None
    if not present and following_wet_contact:
        if previous_phase == "liquid" and envelope.armed(envelope.FILM_FREE_SEED_STORED_SURFACE):
            # Kernel-envelope film-free seed (the lead's record of 2026-10-20 on the 378 K
            # dry wall, section 3 item 1): the first film-free advance after the exterior
            # film depletes is seeded at the stored surface of the last accepted row, not at
            # the bath gas state. A seed only; the path proof and every gate are unchanged.
            film_free_seed = {
                "seed": "stored_surface",
                "seed_temperature_k": surface.temperature_k,
                "seed_y_hexane": surface.surface_y_hexane,
                "gas_state_temperature_k": boundary.temperature_k,
                "gas_state_y_hexane": boundary.y_hexane,
                "previous_phase": previous_phase,
                "record": envelope.SEED_RECORD,
                "physically_qualifying": False,
            }
            if envelope.armed(envelope.FILM_FREE_SEED_GAS_FALLBACK):
                # Kernel-envelope seed order (the lead's numerics queue of 2026-10-20): on the
                # conservative absent elimination's typed chart refusal of the stored seed, the
                # absent set is tried once more from the former bath gas state seed.
                gas_fallback_guess = replace(
                    surface,
                    temperature_k=boundary.temperature_k,
                    surface_y_hexane=boundary.y_hexane,
                )
                film_free_seed["order"] = "stored_surface_then_gas_state"
                film_free_seed["order_record"] = envelope.QUEUE_RECORD
        else:
            absent_guess = replace(
                surface, temperature_k=boundary.temperature_k, surface_y_hexane=boundary.y_hexane
            )
    phases = (
        ("liquid",)
        if present
        else (
            ("contact", "absent", "liquid")
            if surface.external_water_flux_kg_m2_s < 0 and following_wet_contact
            else ("absent", "contact", "liquid")
        )
    )
    # Kernel-envelope ABSENT_SEED_ABOVE_DEW (armed only): an absent seed below the gas's water dew
    # point is lifted above it (a seed only; the set order unchanged)
    absent_seed_lift = None
    if not present and envelope.armed(envelope.ABSENT_SEED_ABOVE_DEW):
        absent_seed_lift = _absent_seed_lift(before, absent_guess, boundary)
        if absent_seed_lift is not None:
            absent_guess = replace(
                absent_guess,
                temperature_k=absent_seed_lift["seed_temperature_k"],
                surface_y_hexane=absent_seed_lift["seed_y_hexane"],
            )
    # Kernel-envelope FILM_FREE_WET_SET_SKIP (armed only): a film-free surface well above the
    # presented gas's water dew point tries the absent set only
    wet_set_skip = None
    if not present and envelope.armed(envelope.FILM_FREE_WET_SET_SKIP):
        wet_set_skip = _film_free_wet_set_skip(before, surface, duration, boundary_factory)
        if wet_set_skip is not None:
            phases = ("absent",)
    rejected = []
    liquid_options = (
        {"stop_at_nonpositive_liquid": True}
        if stop_positive_film_at_nonpositive_liquid and conservative_liquid_elimination and present
        else {}
    )
    event_options = (
        {"conservative_event_elimination": True} if conservative_event_elimination else {}
    )
    absent_options = (
        {"conservative_absent_elimination": True} if conservative_absent_elimination else {}
    )

    def event_profile(seed):
        if previous_accepted_profile is None:
            return {}
        clock, temperatures, compositions = previous_accepted_profile
        gap = before.time_s - clock
        if not (math.isfinite(gap) and gap > 0):
            return {}
        share = seed / gap
        return {
            "event_initial_profile": (
                [a + share * (a - b) for a, b in zip(before.temperatures_k, temperatures)],
                [a + share * (a - b) for a, b in zip(before.y_hexane, compositions)],
            )
        }

    if (
        film_depletion_event_first
        and present
        and predicted_depletion_duration_s is not None
        and math.isfinite(predicted_depletion_duration_s)
        and 0 < predicted_depletion_duration_s <= duration
    ):
        with liquid.liquid_connected_fast_path(
            liquid_mobility_factor,
            exterior_liquid_present=True,
            contact_conductance=contact_conductance,
            **contact_option,
        ) as limiter:
            try:
                result = dry.advance(
                    before,
                    config,
                    surface,
                    carry,
                    duration,
                    boundary,
                    zero_film_event=True,
                    boundary_factory=lambda dt: boundary_factory(before.time_s + dt),
                    search_tolerance=search_tolerance,
                    physical_advance_assembly_budget=assembly_budget,
                    global_component_mass_error_budgets_kg=global_component_mass_error_budgets_kg,
                    use_native_conduction_predictor=use_native_conduction_predictor,
                    event_duration_seed=predicted_depletion_duration_s,
                    **event_options,
                    **event_profile(predicted_depletion_duration_s),
                    **pressure_option,
                )
            except StepBudgetExceeded:
                raise
            except Exception as event_error:
                result = None
                rejected.append(
                    {
                        "phase": "depletion event (predicted first)",
                        "error": str(event_error),
                        "optimizer_debug": getattr(event_error, "optimizer_debug", None),
                    }
                )
            if result is not None:
                after, audit, next_carry, report = result
                native_duration = report["dt_s"]
                if not (
                    math.isfinite(native_duration)
                    and 0 < native_duration <= duration
                    and after.time_s > before.time_s
                    and after.time_s == before.time_s + native_duration
                ):
                    raise DryPhaseAdvanceError(
                        "dry/event duration is outside its declared attempted interval",
                        before,
                        rejected,
                    )
                next_carry = {
                    **next_carry,
                    "external_water_mass_kg": after.grid.areas[-1] * audit.liquid_water_kg_m2,
                }
                return DryPhaseAdvance(
                    after,
                    audit,
                    next_carry,
                    report,
                    "liquid",
                    True,
                    limiter.last_audit,
                    tuple(rejected),
                )
    for phase in phases:
        wet = phase in ("liquid", "contact")
        event = False
        result = None
        with liquid.liquid_connected_fast_path(
            liquid_mobility_factor,
            exterior_liquid_present=wet,
            contact_conductance=contact_conductance,
            **contact_option,
        ) as limiter:
            try:
                result = dry.advance(
                    before,
                    config,
                    surface if wet else absent_guess,
                    carry,
                    duration,
                    boundary,
                    surface_phase_override=None if present else phase,
                    surface_liquid_transfer_adapter=limiter if phase == "contact" else None,
                    search_tolerance=search_tolerance,
                    physical_advance_assembly_budget=assembly_budget,
                    global_component_mass_error_budgets_kg=global_component_mass_error_budgets_kg,
                    use_native_conduction_predictor=use_native_conduction_predictor,
                    conservative_liquid_elimination=conservative_liquid_elimination and present,
                    **(liquid_options if present else {}),
                    **(absent_options if phase == "absent" and not present else {}),
                    **pressure_option,
                )
                if phase == "contact":
                    audit = limiter.last_audit
                    if result[1].liquid_water_kg_m2 != 0.0 or audit is None:
                        raise RuntimeError(
                            "wet contact requires exact zero storage and active resolved supply"
                        )
                    mismatch = abs(
                        audit["transferred_liquid_uptake_kg_m2_s"]
                        - audit["available_water_flux_kg_m2_s"]
                    )
                    if (
                        mismatch * duration * before.grid.areas[-1] / water_inventory_scale_kg
                        > config.dry.ledger_tolerance
                    ):
                        raise RuntimeError("native condensate/contact supply equation refused")
            except StepBudgetExceeded:
                raise
            except Exception as error:
                result = None
                rejected.append(
                    {
                        "phase": phase,
                        "error": str(error),
                        "optimizer_debug": getattr(error, "optimizer_debug", None),
                        **(
                            {"kernel_envelope_film_free_seed": film_free_seed}
                            if film_free_seed is not None and phase == "absent"
                            else {}
                        ),
                    }
                )
                if present:
                    try:
                        result = dry.advance(
                            before,
                            config,
                            surface,
                            carry,
                            duration,
                            boundary,
                            zero_film_event=True,
                            boundary_factory=lambda dt: boundary_factory(before.time_s + dt),
                            search_tolerance=search_tolerance,
                            physical_advance_assembly_budget=assembly_budget,
                            global_component_mass_error_budgets_kg=global_component_mass_error_budgets_kg,
                            use_native_conduction_predictor=use_native_conduction_predictor,
                            **event_options,
                            **event_profile(duration),
                            **pressure_option,
                        )
                        event = True
                    except StepBudgetExceeded:
                        raise
                    except Exception as event_error:
                        result = None
                        rejected.append(
                            {
                                "phase": "depletion event",
                                "error": str(event_error),
                                "optimizer_debug": getattr(event_error, "optimizer_debug", None),
                            }
                        )
                elif (
                    phase == "absent"
                    and envelope.armed(envelope.FILM_FREE_SEED_GAS_FALLBACK)
                    and str(error) == envelope.ABSENT_FILM_CHART_REFUSAL
                ):
                    if gas_fallback_guess is None:
                        gas_fallback_guess = replace(
                            absent_guess,
                            temperature_k=boundary.temperature_k,
                            surface_y_hexane=boundary.y_hexane,
                        )
                        film_free_seed = {
                            "seed": "dispatch_surface",
                            "seed_temperature_k": absent_guess.temperature_k,
                            "seed_y_hexane": absent_guess.surface_y_hexane,
                            "previous_phase": previous_phase,
                            "order": "dispatch_surface_then_gas_state",
                            "order_record": "GT_PS2_ARRIVAL_CHART_SEED_REPAIR_2026-11-04.md",
                            "physically_qualifying": False,
                        }
                    try:
                        result = dry.advance(
                            before,
                            config,
                            gas_fallback_guess,
                            carry,
                            duration,
                            boundary,
                            surface_phase_override=phase,
                            surface_liquid_transfer_adapter=None,
                            search_tolerance=search_tolerance,
                            physical_advance_assembly_budget=assembly_budget,
                            global_component_mass_error_budgets_kg=global_component_mass_error_budgets_kg,
                            use_native_conduction_predictor=use_native_conduction_predictor,
                            conservative_liquid_elimination=False,
                            **absent_options,
                            **pressure_option,
                        )
                        film_free_seed = {
                            **film_free_seed,
                            "accepted_seed": "gas_state_fallback",
                            "first_seed_refusal": str(error),
                        }
                    except StepBudgetExceeded:
                        raise
                    except Exception as fallback_error:
                        result = None
                        rejected.append(
                            {
                                "phase": "absent (gas-state seed after the stored seed)",
                                "error": str(fallback_error),
                                "optimizer_debug": getattr(fallback_error, "optimizer_debug", None),
                            }
                        )
            if result is None:
                continue
            after, audit, next_carry, report = result
            native_duration = report["dt_s"]
            if not (
                math.isfinite(native_duration)
                and 0 < native_duration <= duration
                and after.time_s > before.time_s
                and after.time_s == before.time_s + native_duration
            ):
                raise DryPhaseAdvanceError(
                    "dry/event duration is outside its declared attempted interval",
                    before,
                    rejected,
                )
            next_carry = {
                **next_carry,
                "external_water_mass_kg": after.grid.areas[-1] * audit.liquid_water_kg_m2,
            }
            if absent_seed_lift is not None:
                report = {
                    **report,
                    "kernel_envelope_absent_seed_above_dew": {
                        **absent_seed_lift,
                        "accepted_phase": phase,
                    },
                }
            if wet_set_skip is not None:
                report = {**report, "kernel_envelope_film_free_wet_set_skip": wet_set_skip}
            if film_free_seed is not None:
                seed_stamp = {**film_free_seed, "accepted_phase": phase}
                if "order" in seed_stamp and "accepted_seed" not in seed_stamp:
                    seed_stamp["accepted_seed"] = (
                        "stored_surface" if phase == "absent" else "not_seeded_wet_set"
                    )
                report = {**report, "kernel_envelope_film_free_seed": seed_stamp}
            return DryPhaseAdvance(
                after, audit, next_carry, report, phase, event, limiter.last_audit, tuple(rejected)
            )
    raise DryPhaseAdvanceError(
        "all native surface active sets refused on the exact old state", before, rejected
    )
