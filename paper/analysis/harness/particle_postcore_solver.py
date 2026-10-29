"""Existing DT/Faner dry equations with conservative surface-phase dispatch.

The caller owns the single proof/assembly allowance for a physical advance,
including all rejected active sets and time-step halvings. This dispatcher
never changes the old physical state, a constitutive law or an acceptance gate.
"""

from dataclasses import dataclass, replace
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
    liquid_mobility_factor=100.0,
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
    if contact_conductance is None:
        contact_conductance = envelope.surface_contact_conductance()
    # Kernel-envelope pressure schedule (default None: nothing passed): the
    # pressure of the attempted interval's end, for every active set.
    pressure_option = {} if target_pressure_pa is None else {"target_pressure_pa": target_pressure_pa}
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
                    and gas_fallback_guess is not None
                    and str(error) == envelope.ABSENT_FILM_CHART_REFUSAL
                ):
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
