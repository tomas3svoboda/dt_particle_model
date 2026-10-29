"""Conditional radial fast-wetting comparison, starting at accepted core extinction.

D9 supplies fast/no-drain limits, not a measured liquid uptake coefficient.
This experiment approaches fast chemical-potential relaxation in the opened
radial shell while exterior liquid is present. It is NOT the aggregate D9
instantaneous allocation, a Darcy/Washburn law, or an identified rate.
The base native retained-water flux remains; a parallel dissipative path
uses the actual cells' retained-water activities, with (factor - 1) times
the base mobility. Native component and energy assembly,
capacity active sets, film storage and accepted-step ledgers are reused.
No additional water storage or sorption capacity is introduced.
"""

import argparse
from contextlib import contextmanager
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import pickle
import sys
import time
import traceback
import weakref

import experiment_core2_dt_dry as dry
import particle_kernel_envelope as envelope
from dtdc_simulator.core2 import through_bed_sorbed_water_law as d9


def inventory_bounded_transfer_fraction(available_flux, demand_flux, dt_s):
    """Normalize the D9 drain in mass units, preserving its exact full branch.

    Converting a full drain back to a rate before normalizing can produce
    1 + one ulp. Dividing the drained mass by the requested mass returns
    exactly one when D9 supplies the complete request, without a clamp.
    """
    if available_flux <= 0:
        return 0.0
    requested_mass = demand_flux * dt_s
    drained_mass = d9.film_drain_kg(
        film_kg=float(available_flux * dt_s),
        capacity_kg=float(requested_mass),
        bracket_end=d9.DRAIN_BRACKET_END_NOMINAL,
    )
    return drained_mass / requested_mass


def add_retained_transfer(native, addition, enthalpy, dissipation):
    """Book one additional liquid flux in component, energy and entropy rows."""
    component = replace(
        native.component,
        retained_water_flux_mol_m2_s=native.component.retained_water_flux_mol_m2_s + addition,
        conserved_water_flux_mol_m2_s=native.component.conserved_water_flux_mol_m2_s + addition,
    )
    energy = replace(
        native.energy,
        retained_water_enthalpy_flux_w_m2=native.energy.retained_water_enthalpy_flux_w_m2
        + addition * enthalpy,
        total_energy_flux_w_m2=native.energy.total_energy_flux_w_m2 + addition * enthalpy,
    )
    return replace(
        native,
        component=component,
        energy=energy,
        retained_water_flux_mol_m2_s=native.retained_water_flux_mol_m2_s + addition,
        retained_water_entropy_production_w_m3_k=native.retained_water_entropy_production_w_m3_k
        + dissipation,
    )


def register_while_alive(registry, face, entry):
    """Kernel-envelope FRONT_RELEASE_SOLVE_MEMORY: registry[id(face)] while face is alive.

    The registry is read only through a live face object (supply_limited_surface looks up the
    face it is given), so an entry whose face has gone can never be read again. The entry is
    removed when the face is released, before its id can be reused, and the face itself is not
    held (the former tuple's first slot, which nothing reads, is None). Every lookup of a live
    face returns the entry the former registry returns.
    """
    key = id(face)
    registry[key] = (None, *entry)
    weakref.finalize(face, registry.pop, key, None)


class FiniteContactUndefined(ValueError):
    """The film-covered surface law is undefined at this trial point.

    surface_contact_face divides by K_l (psi_s - psi_p) minus the resisted outward flux; the law
    is defined only where that denominator has the sign of the force (the resisted paths' secant
    conductance above -K_l) and the force is nonzero. A trial point outside is refused as a
    ValueError, as the dry assembly refuses every other point outside its chart (the domain-aware
    Newton backtracks); it is never clamped.
    """


# The surface law's report-only holder: the path decomposition of the most recent film-covered
# surface-face evaluation (the last evaluation of an accepted step is its accepted root:
# experiment_core2_dt_dry.advance re-evaluates result.x), and counters (the evaluations on each
# side of the outward Stefan exemption and the switches between them within a step or an
# advance). Never read by an equation.
FINITE_CONTACT_LAST = {}
FINITE_CONTACT_COUNTS = {
    "surface_evaluations": 0,
    "domain_refusals": 0,
    "outward_stefan_evaluations": 0,
    "inward_or_zero_stefan_evaluations": 0,
    "stefan_branch_switches": 0,
}
_STEFAN_BRANCH = {"previous": None}


def _count_stefan_branch(branch):
    """Count one surface evaluation on ``branch`` and a switch from the previous evaluation's
    branch within the same step or advance (report only)."""
    key = (
        "outward_stefan_evaluations" if branch == "outward" else "inward_or_zero_stefan_evaluations"
    )
    FINITE_CONTACT_COUNTS[key] += 1
    if _STEFAN_BRANCH["previous"] is not None and _STEFAN_BRANCH["previous"] != branch:
        FINITE_CONTACT_COUNTS["stefan_branch_switches"] += 1
    _STEFAN_BRANCH["previous"] = branch


def surface_contact_face(
    native, left, right, aw_left, aw_right, distance, free_mobility, contact_conductance, retained
):
    """The film-covered surface face: the contact on every path, the outward Stefan water flux
    exempt (the model since the final adoption of 2026-10-25, record
    envelope.SURFACE_LAW_FINAL_RECORD).

    Signs are the native ones: a flux is positive outward (cell -> surface). psi = ln a_w,ret
    from the actual retained states, F = psi_s - psi_p = ln(aw_right / aw_left), d = ``distance``
    the outer half cell, K = K_l = ``contact_conductance`` (particle_kernel_envelope
    .K_L_MOL_M2_S, passed by the scope), L_ret = ``retained`` (the native retained-water
    mobility), L_l = ``free_mobility`` = (factor - 1) L_ret. The native face supplies the
    gas-frame water flux N_gas = y_w,adv N_t + J_w, split into its advective (Stefan) part
    N_stefan = N_gas - J_w (= y_w,adv N_t) and its diffusive part J_w.

    While an external water film covers the surface, the water taken up from the film across the
    surface face is in series with the contact conductance; only water vapour that LEAVES the
    pores with the outward gas flow is exempt (it enters the film and condenses there without
    passing the wetting contact):

        N_exempt = max(0, N_stefan),
        N_total = N_exempt + s R,  R = N_nat,free + N_liq,free + N_gas - N_exempt,
        s = K F / (K F - R),  i.e.  N_total = N_exempt - F / (d / L_eff + 1 / K),
        L_eff / d = -R / F,

    with N_nat,free = -(L_ret / d) F the native retained path on the SAME force F as the added
    path (the single retained-flux law of the cells read on the bridged isotherm of the actual
    retained states), N_liq,free = -(L_l / d) F. So where N_stefan <= 0 (the pore gas flowing
    inward, or at rest) every path, the whole gas-frame water flux included, is multiplied by s
    (the law adopted on 2026-10-25 as built and measured as V4, statement for statement); where
    N_stefan > 0 the Stefan part crosses unresisted and J_w, the native retained path and the added
    path are multiplied by s. N_exempt is continuous in N_t (zero at N_t = 0, where both sides give
    s (N_nat,free + N_liq,free + J_w)); its derivative jumps there and is not smoothed. As d -> 0
    the series part tends to -K F, never above; the exempt part is the outward Stefan flux only.

    Booking: the native retained row s N_nat,free (top-level and component), the conserved water
    (gas water) + s N_nat,free, the retained enthalpy flux s N_nat,free times the donor actual
    state's retained enthalpy, the retained dissipation R s L_ret (F / d)^2; the gas-frame water
    row s N_gas with its enthalpy flux s times the native (N_stefan <= 0), or N_stefan + s J_w
    with its enthalpy flux the native times (N_stefan + s J_w) / N_gas (the native upwind gas
    water partial enthalpy for both parts; where N_gas is exactly zero, the upwind cell state's
    partial enthalpy by the new flux's sign); the total energy re-composed by the native fsum;
    the added liquid s N_liq,free returned for add_retained_transfer (donor retained enthalpy,
    dissipation R s L_l (F / d)^2), so the film-stock limiter acts on it as before. The hexane
    row, J_w, N_t (the Newton's unknown, pinned by the hexane rows), the conduction and the gas
    entropy production are unchanged.

    The arithmetic is, statement for statement, the V4 build's (record envelope.SURFACE_LAW_RECORD)
    on the inward side and the outward-carried-exempt candidate C3's (record
    envelope.SURFACE_LAW_STEFAN_C3_RECORD) on the outward side. The holder's ``liquid_former`` and
    ``total_water_former`` report the former face on the same state (no equation reads them).
    """
    k = contact_conductance
    force = math.log(aw_right / aw_left)
    gradient = math.log(aw_right / aw_left) / distance
    former_mobility = distance / (distance / free_mobility + 1.0 / k)
    component = native.component
    n_gas = component.gas_water_flux_mol_m2_s
    n_diffusive = component.independent_water_flux_mol_m2_s
    n_stefan = n_gas - n_diffusive
    outward = n_stefan > 0.0
    n_nat_free = -retained * gradient
    n_liq_free = -free_mobility * gradient
    FINITE_CONTACT_COUNTS["surface_evaluations"] += 1
    _count_stefan_branch("outward" if outward else "inward_or_zero")
    if outward:
        # the outward Stefan water flux exempt: J_w alone of the gas route in the series
        exempt = n_stefan
        resisted = math.fsum((n_nat_free, n_liq_free, n_diffusive))
    else:
        exempt = 0.0
        resisted = math.fsum((n_nat_free, n_liq_free, n_gas))
    denominator = k * force - resisted
    if not (force != 0.0 and math.isfinite(denominator) and denominator * force > 0.0):
        FINITE_CONTACT_COUNTS["domain_refusals"] += 1
        raise FiniteContactUndefined(
            f"finite-contact {envelope.SURFACE_LAW}: the series law is undefined at "
            f"F={force!r} with resisted outward flux {resisted!r} (exempt outward Stefan flux "
            f"{exempt!r}) and K_l={k!r}"
        )
    factor = k * force / denominator
    native_retained = factor * n_nat_free
    mobility = factor * free_mobility
    addition = factor * n_liq_free
    donor = left if addition >= 0 else right
    energy = native.energy
    enthalpy_source = "native_upwind"
    if outward:
        diffusive_applied = factor * n_diffusive
        stefan_applied = n_stefan
        new_gas_water = math.fsum((n_stefan, diffusive_applied))
        if n_gas != 0.0:
            gas_water_energy = energy.gas_water_enthalpy_flux_w_m2 * (new_gas_water / n_gas)
        else:
            upwind = left if new_gas_water > 0 else right
            gas_water_energy = new_gas_water * upwind.water_gas_partial_enthalpy_j_mol
            enthalpy_source = "cell_upwind_native_flux_zero"
    else:
        diffusive_applied = factor * n_diffusive
        stefan_applied = factor * (n_gas - n_diffusive)
        new_gas_water = factor * n_gas
        gas_water_energy = factor * energy.gas_water_enthalpy_flux_w_m2
    retained_energy = native_retained * donor.retained_water_enthalpy_j_mol
    base = replace(
        native,
        component=replace(
            component,
            gas_water_flux_mol_m2_s=new_gas_water,
            retained_water_flux_mol_m2_s=native_retained,
            conserved_water_flux_mol_m2_s=new_gas_water + native_retained,
        ),
        energy=replace(
            energy,
            gas_water_enthalpy_flux_w_m2=gas_water_energy,
            retained_water_enthalpy_flux_w_m2=retained_energy,
            total_energy_flux_w_m2=math.fsum(
                (
                    energy.conductive_heat_flux_w_m2,
                    gas_water_energy,
                    energy.gas_hexane_enthalpy_flux_w_m2,
                    retained_energy,
                )
            ),
        ),
        retained_water_flux_mol_m2_s=native_retained,
        retained_water_entropy_production_w_m3_k=dry.core.bg.R
        * factor
        * retained
        * gradient
        * gradient,
    )
    dissipation = dry.core.bg.R * mobility * gradient * gradient
    path = native.composition_force_path
    FINITE_CONTACT_LAST.clear()
    FINITE_CONTACT_LAST.update(
        variant=envelope.SURFACE_LAW,
        sign="outward positive (cell to surface), mol/(m2 s)",
        K_l_mol_m2_s=k,
        d_m=distance,
        aw_outer=aw_left,
        aw_surface=aw_right,
        force_psi_s_minus_psi_p=force,
        native_retained=native_retained,
        native_retained_former=component.retained_water_flux_mol_m2_s,
        native_retained_free=n_nat_free,
        native_retained_force_gradient_m_inv=native.retained_water_chemical_force_gradient_m_inv,
        native_actual_path_force=(
            None if path is None else path.retained_force_difference_right_minus_left
        ),
        total_stefan_unknown=component.total_stefan_flux_mol_m2_s,
        surface_stefan_direction=(
            "outward"
            if component.total_stefan_flux_mol_m2_s > 0
            else "inward" if component.total_stefan_flux_mol_m2_s < 0 else "zero"
        ),
        surface_stefan_magnitude_mol_m2_s=abs(component.total_stefan_flux_mol_m2_s),
        # the physical total gas flux under the law (the ledger's recovered one): N_t - (1 - s) N_gas
        # where N_stefan <= 0, N_t - (1 - s) J_w where N_stefan > 0
        surface_physical_total_gas_flux_mol_m2_s=math.fsum(
            (new_gas_water, component.gas_hexane_flux_mol_m2_s)
        ),
        gas_hexane=component.gas_hexane_flux_mol_m2_s,
        gas_hexane_direction=(
            "outward"
            if component.gas_hexane_flux_mol_m2_s > 0
            else "inward" if component.gas_hexane_flux_mol_m2_s < 0 else "zero"
        ),
        gas_water_native=n_gas,
        gas_water_stefan_native=n_gas - n_diffusive,
        gas_water_diffusive_native=n_diffusive,
        liquid_free=n_liq_free,
        liquid_former=-former_mobility * gradient,
        resisted_outward_flux=resisted,
        gas_factor=factor,
        gas_water_applied=new_gas_water,
        gas_water_stefan_applied=stefan_applied,
        gas_water_diffusive_applied=diffusive_applied,
        liquid_applied=addition,
        total_water_applied=math.fsum((native_retained, new_gas_water, addition)),
        total_water_former=math.fsum(
            (component.retained_water_flux_mol_m2_s, n_gas, -former_mobility * gradient)
        ),
        contact_limit_minus_K_F=-k * force,
        stefan_branch="outward" if outward else "inward_or_zero",
        gas_water_exempt=exempt,
        series_total=math.fsum(
            (native_retained, diffusive_applied if outward else new_gas_water, addition)
        ),
        limit_exempt_minus_K_F=exempt - k * force,
        gas_water_enthalpy_source=enthalpy_source,
    )
    return base, addition, donor.retained_water_enthalpy_j_mol, dissipation


def surface_contact_advance_stamp(native_report, liquid_supply_audit):
    """The surface law's per-advance dry-history stamp (report only).

    The accepted root's surface paths (FINITE_CONTACT_LAST, None when the advance carried no
    exterior liquid), this advance's evaluation and domain-refusal counts, the film-stock supply
    limiter's audit (the contact active set only; the positive-film set closes the film stock by
    the conservative elimination and the depletion event), and the identity of the stamped paths
    with the accepted ledger's surface face (retained: exact; combined: storage-telescoped). The
    keys and their order are the V4 build's; ``path_forces`` names the path-face rules in force
    (particle_kernel_envelope.path_forces_in_force).
    """
    paths = dict(FINITE_CONTACT_LAST) if FINITE_CONTACT_LAST else None
    check = None
    if paths is not None:
        ledger = native_report["ledger"]
        recovery = ledger["stefan_recovery"] if isinstance(ledger, dict) else ledger.stefan_recovery
        get = recovery.get if isinstance(recovery, dict) else lambda key: getattr(recovery, key)
        retained = get("retained_water_fluxes_mol_m2_s")[-1]
        combined = get("combined_component_fluxes_mol_m2_s")[-1]
        stefan = get("total_stefan_fluxes_mol_m2_s")[-1]
        stamped_retained = math.fsum((paths["native_retained"], paths["liquid_applied"]))
        stamped_combined = math.fsum((paths["total_water_applied"], paths["gas_hexane"]))
        check = {
            "ledger_surface_retained": retained,
            "stamped_retained": stamped_retained,
            "retained_equal": stamped_retained == retained,
            "ledger_surface_combined": combined,
            "stamped_combined": stamped_combined,
            "combined_relative_difference": (stamped_combined - combined)
            / max(abs(combined), 1e-300),
            "ledger_recovered_total_stefan": stefan,
            "stamped_physical_total_gas": math.fsum(
                (paths["gas_water_applied"], paths["gas_hexane"])
            ),
        }
    return {
        "variant": envelope.SURFACE_LAW,
        "K_l_mol_m2_s": None if paths is None else paths["K_l_mol_m2_s"],
        "paths": paths,
        "accepted_face_check": check,
        "counts": dict(FINITE_CONTACT_COUNTS),
        "liquid_supply_limiter_audit": liquid_supply_audit,
        "path_forces": envelope.path_forces_in_force(),
        "record": envelope.SURFACE_LAW_FINAL_RECORD,
        "physically_qualifying": False,
    }


def surface_contact_front_stamp(old_liquid_kg_m2, elapsed_s, surface, limiter_audit):
    """The surface law's per-step front-stage stamp (report only; front.history.jsonl).

    The paths of the last surface-face evaluation (FINITE_CONTACT_LAST), this step's evaluation
    and domain-refusal counts, the limiter's last audit, and the identity of the stamped total
    water with the accepted surface node's water balance: the node books the pore water
    j_p = (L - L_old) / dt + j_ext (kg/m2 s, outward positive), which equals the stamped total
    times M_w when the last evaluation is the accepted root and the limiter transferred the
    whole added path. Keys and order as the V4 build's, plus ``path_forces``.
    """
    paths = dict(FINITE_CONTACT_LAST) if FINITE_CONTACT_LAST else None
    check = None
    if paths is not None and hasattr(surface, "liquid_water_kg_m2"):
        node = math.fsum(
            (
                (surface.liquid_water_kg_m2 - old_liquid_kg_m2) / elapsed_s,
                surface.external_water_flux_kg_m2_s,
            )
        )
        stamped = paths["total_water_applied"] * dry.core.wa.M
        check = {
            "node_pore_water_kg_m2_s": node,
            "stamped_total_water_kg_m2_s": stamped,
            "relative_difference": (stamped - node) / max(abs(node), 1e-300),
        }
    return {
        "variant": envelope.SURFACE_LAW,
        "paths": paths,
        "accepted_node_check": check,
        "counts": dict(FINITE_CONTACT_COUNTS),
        "liquid_supply_limiter_audit": limiter_audit,
        "path_forces": envelope.path_forces_in_force(),
        "record": envelope.SURFACE_LAW_FINAL_RECORD,
        "physically_qualifying": False,
    }


def surface_contact_reset():
    """Clear the report-only holder and the counters before a step or an advance."""
    FINITE_CONTACT_LAST.clear()
    for key in FINITE_CONTACT_COUNTS:
        FINITE_CONTACT_COUNTS[key] = 0
    _STEFAN_BRANCH["previous"] = None


# LIQUID_MOBILITY_RULE (study knob, armed only; record envelope.SURFACE_RULES_RECORD). Report only;
# never read by an equation.
SURFACE_RULES_LAST = {}
SURFACE_RULES_COUNTS = {
    "interior_evaluations": 0,
    "interior_evaluations_liquid_active": 0,
}
# local_state_lagged only: interior faces the lagged rule could not match to their two cells (they
# fall back to the trial-state rule); report only, stamped only under local_state_lagged.
LAGGED_COUNTS = {"interior_evaluations_unmatched": 0}


def surface_rules_reset():
    """Clear the report-only holder and the counters before an advance (study knob only)."""
    SURFACE_RULES_LAST.clear()
    for key in SURFACE_RULES_COUNTS:
        SURFACE_RULES_COUNTS[key] = 0
    LAGGED_COUNTS["interior_evaluations_unmatched"] = 0


def surface_rules_advance_stamp(mobility_rule, native_report, after, config):
    """The liquid-mobility study knob's per-advance history stamp (report only).

    local_state: the interior faces whose rule holds at the accepted state (cells' actual retained
    loadings against W_cap) and the evaluation counts; local_state_lagged adds the start state's
    active faces and the unmatched count. ``native_report`` is accepted for the former signature
    and not read.
    """
    stamp = {
        "liquid_mobility_rule": mobility_rule,
        "counts": dict(SURFACE_RULES_COUNTS),
        "record": envelope.SURFACE_RULES_RECORD,
        "physically_qualifying": False,
    }
    pore = config.dry.pore
    cells = [
        dry.cp.evaluate_equilibrium(t, config.dry.pressure_pa, y, pore)
        for t, y in zip(after.temperatures_k, after.y_hexane)
    ]
    loadings = [c.retained_water_loading for c in cells]
    above = [w > pore.luikov.W_cap for w in loadings]
    active = [i for i in range(len(cells) - 1) if above[i] or above[i + 1]]
    stamp["local_state"] = {
        "W_cap": pore.luikov.W_cap,
        "cells_above_cap": [i for i, flag in enumerate(above) if flag],
        "interior_faces": len(cells) - 1,
        "interior_faces_active_at_accepted_state": active,
        "outer_cell_loading": loadings[-1],
        "max_loading": max(loadings),
    }
    if mobility_rule == envelope.LIQUID_MOBILITY_LOCAL_STATE_LAGGED:
        start = SURFACE_RULES_LAST.get("start_state_loadings") or []
        start_above = [w > pore.luikov.W_cap for w in start]
        stamp["local_state"].update(
            state_read="the advance's accepted start state",
            start_state_cells_above_cap=[i for i, flag in enumerate(start_above) if flag],
            interior_faces_active_in_the_advance=[
                i for i in range(len(start) - 1) if start_above[i] or start_above[i + 1]
            ],
            interior_evaluations_unmatched=LAGGED_COUNTS["interior_evaluations_unmatched"],
        )
    return stamp


# ---- The path-face rules (the final adoption of 2026-10-25; record envelope.SURFACE_LAW_FINAL_RECORD,
# measured as K1 and K2 of the cap-and-liquid-force study, record envelope.CAP_FORCE_RECORD) ------
#
# K1, the cap rule: the native actual-path classifier (actual_composition_force_path
# ._active_set_on_segment) books a path segment whose water activity lies above the Luikov cap as
# RETAINED_CAP, a zero retained force (the pre-bridge capacity graph). Under the ruled imbibition
# bridge the retained water above the cap follows the gas activity, which every interior face's
# native retained flux carries. Inside every liquid-connected scope the classifier returns
# LUIKOV_SMOOTH for such a segment, so the native path integrates the smooth retained potential
# over it (the ideal part plus the second-virial quadrature) on every actual-path face: the moving
# wet-front interface face, the film-free surface face and the newborn face tangent. The native
# classifier is always called first, so its seam-resolution refusal and its probes are unchanged.
# The film-covered surface face's native evaluation keeps the native classification (the surface
# law reads only its gas rows; its own retained path is N_nat,free on the isotherm force F), so
# that face is the measured one bit for bit. The src classifier is not edited: the rule acts only
# inside the harness's scopes (the tower path is not moved).
#
# K2, the liquid force: on an interior face the added liquid path's force ln(a_R / a_L) is read
# from both cells' equilibrium retained loadings at the common face temperature (T_L + T_R) / 2,
# the temperature the native interior path uses, each at its own composition, in place of each
# cell's own temperature (whose thermal term (d ln a_w / dT)_y Delta T drove the same water by a
# second force). The donor's retained enthalpy, the mobility and the dissipation's form are
# unchanged; the surface face keeps its one force F.
#
# The suppression depth is the only state an equation reads: nonzero exactly while the wrapper
# evaluates the native face of a film-covered surface face.
_CAP_SUPPRESS = {"depth": 0}


def bridged_cap_segment_classifier(native_classifier):
    """K1's actual-path segment classifier (RETAINED_CAP read as LUIKOV_SMOOTH)."""
    smooth = dry.cp.RetainedWaterActiveSet.LUIKOV_SMOOTH
    capped = dry.cp.RetainedWaterActiveSet.RETAINED_CAP

    def classifier(activities, start, end, activity_at_cap):
        active_set = native_classifier(activities, start, end, activity_at_cap)
        if active_set is capped and not _CAP_SUPPRESS["depth"]:
            return smooth
        return active_set

    classifier.native = native_classifier
    classifier.cap_rule = envelope.PATH_FORCES
    return classifier


@contextmanager
def bridged_cap_rule_scope(active):
    """Install K1's classifier for the scope (in memory only; restored on exit). ``active`` False
    (the former path forces) installs nothing; an enclosing scope's classifier is kept as it is."""
    if not active:
        yield False
        return
    acfp = dry.core.cut.acfp
    current = acfp._active_set_on_segment
    if getattr(current, "cap_rule", None) == envelope.PATH_FORCES:
        yield True
        return
    acfp._active_set_on_segment = bridged_cap_segment_classifier(current)
    try:
        yield True
    finally:
        acfp._active_set_on_segment = current


@contextmanager
def native_cap_classification():
    """Inside a K1 scope: the native segment classification for this evaluation (the film-covered
    surface face's native evaluation)."""
    _CAP_SUPPRESS["depth"] += 1
    try:
        yield
    finally:
        _CAP_SUPPRESS["depth"] -= 1


def common_temperature_activities(args):
    """K2: the two cells' retained water activities for the added liquid path's force on an
    interior face, from their equilibrium retained loadings at the common face temperature."""
    config = args[8]
    face_temperature = 0.5 * (args[0] + args[2])
    force_left = dry.cp.evaluate_equilibrium(
        face_temperature, config.dry.pressure_pa, args[1], args[6]
    )
    force_right = dry.cp.evaluate_equilibrium(
        face_temperature, config.dry.pressure_pa, args[3], args[7]
    )
    aw_left = dry.core.sp.water_activity(force_left.retained_water_loading, args[6].luikov)
    aw_right = dry.core.sp.water_activity(force_right.retained_water_loading, args[7].luikov)
    return aw_left, aw_right


@contextmanager
def liquid_connected_fast_path(
    factor,
    *,
    exterior_liquid_present,
    contact_conductance=None,
    exclude_wet_front=False,
    liquid_mobility_rule=None,
    liquid_mobility_start_state=None,
):
    """Add a named parallel path only on an exterior-liquid-connected dry grid.

    contact_conductance (K_l; the Paper-1 scopes pass particle_kernel_envelope.K_L_MOL_M2_S): the
    surface face under the exterior liquid is then the film-covered surface law,
    surface_contact_face (the contact in series with every water path across the outer half
    cell, the outward Stefan water flux exempt). Without a contact conductance (the fast-wetting
    comparison of this module's main) the surface face carries the added path unresisted, as
    every interior face does.

    The path-face rules (K1 the cap rule, K2 the liquid force; see the block above) are the
    model's: read once at scope entry through particle_kernel_envelope.former_path_forces(),
    which is False unless the one-line revert MODEL_REVERT_PATH_FORCES or the labelled historical
    switch's value former_cap_rule_and_own_temperature_force restores the former statements.

    liquid_mobility_rule (default None, the former rule): a declared value of the study knob
    LIQUID_MOBILITY_RULE, passed only by the dry-stage dispatcher; armed, every face is assembled
    by rules_flux, which applies the added path on an interior face where at least one cell's
    actual retained loading exceeds W_cap, film or not; unset, flux executes its statements.

    liquid_mobility_start_state (local_state_lagged only, required there): the advance's accepted
    start state as (temperatures_k, y_hexane, its transport config); the rule reads its cells'
    loadings, so the active set is fixed for every trial of the advance.
    """
    if not math.isfinite(factor) or factor < 1:
        raise ValueError("fast-wetting multiplier must be finite and >= 1")
    if contact_conductance is not None and (
        not math.isfinite(contact_conductance) or contact_conductance <= 0
    ):
        raise ValueError("liquid contact conductance must be finite and positive")
    # LIQUID_MOBILITY_RULE (study knob; default None, the former rule): declared values only,
    # passed only by the dry-stage dispatcher.
    rules = liquid_mobility_rule is not None
    if rules and liquid_mobility_rule not in envelope.VALUED_ITEMS[envelope.LIQUID_MOBILITY_RULE]:
        raise envelope.EnvelopeKnobRefused(
            f"liquid mobility rule {liquid_mobility_rule!r} is not a declared value"
        )
    lagged = liquid_mobility_rule == envelope.LIQUID_MOBILITY_LOCAL_STATE_LAGGED
    if lagged != (liquid_mobility_start_state is not None):
        raise envelope.EnvelopeKnobRefused(
            "LIQUID_MOBILITY_RULE=local_state_lagged needs the advance's start state, and only it"
        )
    if lagged:
        # local_state_lagged: the start state's cell loadings, read once, before any trial.
        start_temperatures, start_compositions, start_config = liquid_mobility_start_state
        start_loadings = [
            cell.retained_water_loading
            for cell in dry.ct._equilibrium_states(
                start_temperatures, start_compositions, start_config
            )
        ]
        trial = {"temperatures": (), "compositions": ()}
        native_states = dry.ct._equilibrium_states

        def recording_states(temperatures, compositions, *rest, **options):
            # The trial's cells (report of the face identification only; values unchanged).
            trial["temperatures"] = tuple(temperatures)
            trial["compositions"] = tuple(compositions)
            return native_states(temperatures, compositions, *rest, **options)

        SURFACE_RULES_LAST["start_state_loadings"] = list(start_loadings)
    # The former surface face (the added liquid path alone in series with the contact): only
    # through the one-line revert MODEL_REVERT_SURFACE_LAW or the labelled historical switch's
    # value former_unresisted_routes (never a model option).
    former_surface_law = envelope.former_surface_face()
    # The path-face rules K1 and K2 (the model), unless reverted (MODEL_REVERT_PATH_FORCES or the
    # historical switch's value former_cap_rule_and_own_temperature_force).
    bridged_cap = common_temperature = not envelope.former_path_forces()
    original = dry.core.cut._dry_flux
    surface_pairs = {}
    # Kernel-envelope FRONT_RELEASE_SOLVE_MEMORY (armed only; result-neutral): an entry lives
    # only while its face does. Off, every entry lives as long as this scope.
    release = envelope.armed(envelope.FRONT_RELEASE_SOLVE_MEMORY)

    def register(augmented, entry):
        if release:
            register_while_alive(surface_pairs, augmented, entry)
        else:
            surface_pairs[id(augmented)] = (augmented, *entry)

    def contact_face(native, left, right, aw_left, aw_right, distance, mobility, config):
        # The adopted film-covered surface law on the surface face; registered for the limiter.
        base, addition, enthalpy, dissipation = surface_contact_face(
            native,
            left,
            right,
            aw_left,
            aw_right,
            distance,
            mobility,
            contact_conductance,
            config.dry.retained_water_mobility.value_mol_m_s,
        )
        augmented = add_retained_transfer(base, addition, enthalpy, dissipation)
        register(augmented, (base, addition, enthalpy, dissipation))
        return augmented

    def rules_flux(*args, **kwargs):
        # LIQUID_MOBILITY_RULE (study knob, armed only). Every branch the armed knob does not name
        # executes the statements of flux below.
        is_wet_front = (
            "temperature_gradient_k_m" in kwargs or "composition_gradient_m_inv" in kwargs
        )
        tagged = "thermodynamic_force_temperature_k" in kwargs
        if factor == 1 or (exclude_wet_front and is_wet_front):
            return original(*args, **kwargs)
        local = not tagged
        if not exterior_liquid_present and not local:
            return original(*args, **kwargs)
        config = args[8]
        if bridged_cap and contact_conductance is not None and tagged and not former_surface_law:
            # K1: the film-covered surface face's native evaluation keeps the native classification
            with native_cap_classification():
                native = original(*args, **kwargs)
        else:
            native = original(*args, **kwargs)
        left = dry.cp.evaluate_equilibrium(args[0], config.dry.pressure_pa, args[1], args[6])
        right = dry.cp.evaluate_equilibrium(args[2], config.dry.pressure_pa, args[3], args[7])
        if local:
            # local_state: the added path acts where at least one cell's actual retained loading
            # exceeds its W_cap, with or without an exterior film.
            SURFACE_RULES_COUNTS["interior_evaluations"] += 1
            left_loading = left.retained_water_loading
            right_loading = right.retained_water_loading
            if lagged:
                # local_state_lagged: the face's two cells on the advance's start state.
                temperatures, compositions = trial["temperatures"], trial["compositions"]
                index = next(
                    (
                        j
                        for j in range(len(temperatures) - 1)
                        if temperatures[j] == args[0]
                        and compositions[j] == args[1]
                        and temperatures[j + 1] == args[2]
                        and compositions[j + 1] == args[3]
                    ),
                    None,
                )
                if index is None or len(temperatures) != len(start_loadings):
                    LAGGED_COUNTS["interior_evaluations_unmatched"] += 1
                else:
                    left_loading = start_loadings[index]
                    right_loading = start_loadings[index + 1]
            if not (left_loading > args[6].luikov.W_cap or right_loading > args[7].luikov.W_cap):
                return native
            SURFACE_RULES_COUNTS["interior_evaluations_liquid_active"] += 1
        if common_temperature and local:
            # K2: the interior added path's force at the common face temperature
            aw_left, aw_right = common_temperature_activities(args)
        else:
            aw_left = dry.core.sp.water_activity(left.retained_water_loading, args[6].luikov)
            aw_right = dry.core.sp.water_activity(right.retained_water_loading, args[7].luikov)
        gradient = math.log(aw_right / aw_left) / args[4]
        mobility = (factor - 1) * config.dry.retained_water_mobility.value_mol_m_s
        if contact_conductance is not None and tagged:
            if not former_surface_law:
                return contact_face(
                    native, left, right, aw_left, aw_right, args[4], mobility, config
                )
            mobility = args[4] / (args[4] / mobility + 1.0 / contact_conductance)
        addition = -mobility * gradient
        donor = left if addition >= 0 else right
        dissipation = dry.core.bg.R * mobility * gradient * gradient
        augmented = add_retained_transfer(
            native, addition, donor.retained_water_enthalpy_j_mol, dissipation
        )
        if tagged:
            register(
                augmented, (native, addition, donor.retained_water_enthalpy_j_mol, dissipation)
            )
        return augmented

    def flux(*args, **kwargs):
        if rules:
            return rules_flux(*args, **kwargs)
        is_wet_front = (
            "temperature_gradient_k_m" in kwargs or "composition_gradient_m_inv" in kwargs
        )
        if factor == 1 or not exterior_liquid_present or (exclude_wet_front and is_wet_front):
            return original(*args, **kwargs)
        config = args[8]
        tagged = "thermodynamic_force_temperature_k" in kwargs
        if bridged_cap and contact_conductance is not None and tagged and not former_surface_law:
            # K1: the film-covered surface face's native evaluation keeps the native classification
            with native_cap_classification():
                native = original(*args, **kwargs)
        else:
            native = original(*args, **kwargs)
        left = dry.cp.evaluate_equilibrium(args[0], config.dry.pressure_pa, args[1], args[6])
        right = dry.cp.evaluate_equilibrium(args[2], config.dry.pressure_pa, args[3], args[7])
        if common_temperature and not tagged:
            # K2: the interior added path's force from both cells' equilibrium retained loadings at
            # the common face temperature (the native interior path's), each at its own composition.
            # The former reading below (actual W at each cell's own temperature) was chosen because
            # a common face T can reverse which cell is capped; at a manufactured thermal-cap corner
            # without the imbibition bridge K2 points the added path into the capped cell (a finding
            # of the final adoption's record; the corner's test runs on the former path forces).
            aw_left, aw_right = common_temperature_activities(args)
        else:
            # The surface face (and, reverted, every face): the actual cells' retained activities.
            # D9's available capacity is that of the ACTUAL retained inventory. Re-evaluating it at
            # a fictitious common face T can reverse which cell is capped; the added liquid path
            # therefore reads actual W, and removes pure-liquid T dependence by using water activity.
            aw_left = dry.core.sp.water_activity(left.retained_water_loading, args[6].luikov)
            aw_right = dry.core.sp.water_activity(right.retained_water_loading, args[7].luikov)
        gradient = math.log(aw_right / aw_left) / args[4]
        mobility = (factor - 1) * config.dry.retained_water_mobility.value_mol_m_s
        # The exterior liquid contact: the film-covered surface law on the surface face (the cut
        # assembly marks its wet-core trace with the gradient keywords).
        if contact_conductance is not None and "thermodynamic_force_temperature_k" in kwargs:
            if not former_surface_law:
                return contact_face(
                    native, left, right, aw_left, aw_right, args[4], mobility, config
                )
            # Reverted (MODEL_REVERT_SURFACE_LAW) or the historical switch's former face: the
            # added liquid path alone in series with the contact.
            mobility = args[4] / (args[4] / mobility + 1.0 / contact_conductance)
        addition = -mobility * gradient
        donor = left if addition >= 0 else right
        dissipation = dry.core.bg.R * mobility * gradient * gradient
        augmented = add_retained_transfer(
            native, addition, donor.retained_water_enthalpy_j_mol, dissipation
        )
        if "thermodynamic_force_temperature_k" in kwargs:
            register(
                augmented, (native, addition, donor.retained_water_enthalpy_j_mol, dissipation)
            )
        return augmented

    def supply_limited_surface(
        face,
        *,
        old_liquid_loading_kg_m2,
        dt_s,
        external_water_flux_kg_m2_s,
        requested_share=None,
    ):
        if factor == 1 or not exterior_liquid_present:
            return face
        _, native, addition, enthalpy, dissipation = surface_pairs[id(face)]
        available = (
            old_liquid_loading_kg_m2 / dt_s
            + native.component.conserved_water_flux_mol_m2_s * dry.core.wa.M
            - external_water_flux_kg_m2_s
        )
        demand = -addition * dry.core.wa.M
        if demand <= 0:
            return native  # D9's boundary transfer is film -> shell only.
        bounded_share = inventory_bounded_transfer_fraction(available, demand, dt_s)
        share = bounded_share if requested_share is None else requested_share
        if not 0 <= share <= 1:
            raise ValueError("wet-contact transfer fraction outside [0,1]")
        transferred = demand * share
        supply_limited_surface.last_audit = {
            "available_water_flux_kg_m2_s": float(available),
            "requested_liquid_uptake_kg_m2_s": float(demand),
            "transferred_liquid_uptake_kg_m2_s": float(transferred),
            "fraction_of_requested_path": float(share),
            "D9_supply_bounded_fraction": float(bounded_share),
            "old_liquid_loading_kg_m2": float(old_liquid_loading_kg_m2),
            "dt_s": float(dt_s),
            "D9_inventory_bounded_transfer_used": True,
        }
        return add_retained_transfer(native, addition * share, enthalpy, dissipation * share)

    supply_limited_surface.last_audit = None

    if lagged:
        # LIQUID_MOBILITY_RULE=local_state_lagged (armed only): record each _equilibrium_states
        # call's arrays so an interior face can be matched to its two cells.
        dry.ct._equilibrium_states = recording_states
    # K1 (the model, unless reverted): the bridged cap classifier for the whole scope.
    with bridged_cap_rule_scope(bridged_cap):
        dry.core.cut._dry_flux = flux
        try:
            yield supply_limited_surface
        finally:
            dry.core.cut._dry_flux = original
            if lagged:
                dry.ct._equilibrium_states = native_states


def inventory(state):
    values = dry.ct._inventories(
        state.grid,
        dry.ct._equilibrium_states(state.temperatures_k, state.y_hexane, state.config),
    )
    return (
        values.total_water_mol * dry.core.wa.M,
        values.total_hexane_mol * dry.core.hx.M,
        values.total_energy_j,
    )


def row(state, surface, carry):
    w, h, u = inventory(state)
    volume = math.fsum(state.grid.volumes)
    mass = state.config.pore.rho_dm_p * volume
    film = state.grid.areas[-1] * carry["external_water_loading_kg_m2"]
    cells = dry.ct._equilibrium_states(state.temperatures_k, state.y_hexane, state.config)
    return {
        "time_s": state.time_s,
        "temperature_k": math.fsum(t * v for t, v in zip(state.temperatures_k, state.grid.volumes))
        / volume,
        "particle_water_kg_kg_dry": w / mass,
        "external_water_kg_kg_dry": film / mass,
        "hexane_kg_kg_dry": h / mass,
        "moisture_wet_basis": (w + film) / (mass + w + film + h),
        "retained_water_cells_kg_kg_dry": [c.retained_water_loading for c in cells],
        "y_hexane_cells": list(state.y_hexane),
        "combined_W_H_U": [w + film, h, u + carry["external_water_energy_j"]],
        "surface": asdict(surface),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--factor", type=float, required=True)
    parser.add_argument("--contact-conductance", type=float)
    parser.add_argument("--dt", type=float, default=2.5)
    parser.add_argument("--until", type=float, default=1766.5355008060146)
    parser.add_argument("--steps", type=int, default=1500)
    parser.add_argument("--max-halvings", type=int, default=12)
    args = parser.parse_args()
    if args.dt <= 0 or args.steps <= 0 or args.max_halvings < 0:
        raise ValueError("positive timestep and step budget required")
    if args.out.exists() or args.out.with_suffix(".accepted.pkl").exists():
        raise ValueError("refusing to overwrite evidence")
    checkpoint = pickle.loads(args.resume.read_bytes())
    current = checkpoint["current"]
    if hasattr(current, "to_uniform_fully_dry_transport_state"):
        current = current.to_uniform_fully_dry_transport_state()
    if not isinstance(current, dry.ct.FullyDryTransportState):
        raise TypeError("this first comparison requires an accepted fully opened radial grid")
    current = replace(current, config=replace(current.config, nonlinear_residual_tolerance=1e-10))
    config, surface, carry = checkpoint["cut_config"], checkpoint["surface"], checkpoint["carry"]
    bath_path = dry.core.DESTINY_PATH.with_name("destiny_v3_exit105.json")
    bath = json.loads(bath_path.read_text())
    root = dry.driver.REPO_ROOT.resolve()
    sources = {Path(__file__).resolve(), bath_path.resolve()}
    sources.update(root.joinpath("src/dtdc_simulator/core2").rglob("*.py"))
    for module in list(sys.modules.values()):
        name = getattr(module, "__file__", None)
        if name and Path(name).suffix == ".py" and Path(name).resolve().is_relative_to(root):
            sources.add(Path(name).resolve())
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)}
    for name, expected in checkpoint["source_sha256"].items():
        path = Path(name)
        if path.is_relative_to(root / "src/dtdc_simulator/core2"):
            if hashes.get(str(path.resolve())) != expected:
                raise ValueError(f"parent native source changed: {path}")
    initial = row(current, surface, carry)
    mass = current.config.pore.rho_dm_p * math.fsum(current.grid.volumes)
    area = current.grid.areas[-1]
    integrated = [0.0, 0.0, 0.0]
    maximum = [0.0, 0.0, 0.0]
    scales = [max(abs(v), 1e-300) for v in initial["combined_W_H_U"]]
    # Energy uses an inventory/capacity scale rather than an arbitrary datum cancellation.
    scales[2] = max(scales[2], mass * current.config.pore.cp_dry_meal)
    record = {
        "experiment": "post-front radial liquid-connected fast-wetting limit approach",
        "physically_qualifying": False,
        "full_feed_to_exit_fast_imbibition_march": False,
        "aggregate_D9_instantaneous_drain_implemented": False,
        "base_transport_coefficient_changed": False,
        "additional_parallel_path_total_to_base_mobility_factor": args.factor,
        "finite_liquid_contact_conductance_mol_m2_s": args.contact_conductance,
        "parallel_path_force": "gradient of log retained-water activity from actual local retained inventory; no explicit thermal-diffusion term",
        "parallel_path_caloric": "native retained-water partial enthalpy at actual donor temperature; zero native excess binding heat",
        "activation": "positive exterior water or supply-limited saturated contact with resolved condensate; fully opened radial shell only",
        "native_capacity_kg_kg_dry": config.dry.pore.luikov.W_cap,
        "capacity_extension_to_D9_x_bf_used": False,
        "parent_checkpoint": str(args.resume),
        "parent_sha256": hashlib.sha256(args.resume.read_bytes()).hexdigest(),
        "source_sha256_at_entry": hashes,
        "initial": initial,
        "attempts": [],
        "history": [initial],
        "dt_requested_s": args.dt,
        "step_halving": f"on rejected step, exact rollback; at most {args.max_halvings} halvings",
        "absent_surface_initial_guess": "bulk bath T and y after depletion or wet contact; otherwise old surface; physical particle and old film storage unchanged",
        "radial_fast_wetting_force_version": "actual_retained_activity_v2",
        "zero_inventory_wet_contact": "water activity unity; extra liquid uptake bounded by resolved condensate/base-pore supply using D9 film_drain_kg; may form new positive exterior film",
        "native_and_combined_ledger_tolerance": 1e-10,
    }
    started = time.perf_counter()
    fatal = None
    for index in range(args.steps):
        if current.time_s >= args.until:
            break
        requested = min(args.dt, args.until - current.time_s)
        accepted = False
        for halving in range(args.max_halvings + 1):
            duration = requested / (2**halving)
            present = bool(carry["external_water_loading_kg_m2"] > 0)
            event = False
            try:
                boundary = dry.boundary_driver.bath_boundary(
                    bath, current.grid.R, current.time_s + duration
                )
                surface_guess = surface
                just_depleted = (
                    len(record["history"]) >= 2
                    and record["history"][-1]["external_water_kg_kg_dry"] == 0
                    and record["history"][-2]["external_water_kg_kg_dry"] > 0
                )
                previous_contact = any(
                    a.get("accepted") and a.get("zero_inventory_wet_contact")
                    for a in record["attempts"][-1:]
                )
                if not present and (just_depleted or previous_contact):
                    # At film depletion the old saturated trace can put every
                    # finite-difference trial on the wrong side of the phase
                    # boundary. The new gas-only node has independent T/y:
                    # initialize those unknowns from the admissible bath.
                    # This surface argument is ONLY a nonlinear initial guess;
                    # old film energy/T remain in carry, and current is exact.
                    surface_guess = replace(
                        surface,
                        temperature_k=boundary.temperature_k,
                        surface_y_hexane=boundary.y_hexane,
                    )
                phases = (
                    ("liquid",)
                    if present
                    else (
                        ("contact", "absent", "liquid")
                        if surface.external_water_flux_kg_m2_s < 0
                        and (just_depleted or previous_contact)
                        else ("absent", "contact", "liquid")
                    )
                )
                result = None
                last_error = None
                limiter_audit = None
                for phase in phases:
                    wet_contact = phase in ("liquid", "contact")
                    guess = surface if wet_contact else surface_guess
                    with (
                        liquid_connected_fast_path(
                            args.factor,
                            exterior_liquid_present=wet_contact,
                            contact_conductance=args.contact_conductance,
                        ) as limiter,
                        dry.driver.bounded_path_proof(3000000, 120.0) as path_audit,
                        dry.driver.path_bounds.monotone_water_path_enclosure([]),
                        dry.study.magnitude_storage_scale(),
                    ):
                        try:
                            result = dry.advance(
                                current,
                                config,
                                guess,
                                carry,
                                duration,
                                boundary,
                                surface_phase_override=None if present else phase,
                                surface_liquid_transfer_adapter=(
                                    limiter if phase == "contact" else None
                                ),
                            )
                            if phase == "contact":
                                contact_audit = limiter.last_audit
                                if result[1].liquid_water_kg_m2 != 0.0:
                                    raise RuntimeError(
                                        "wet contact must resolve exactly zero film storage"
                                    )
                                if contact_audit is None:
                                    raise RuntimeError(
                                        "wet contact has no active inward liquid transfer"
                                    )
                                deficit = abs(
                                    contact_audit["transferred_liquid_uptake_kg_m2_s"]
                                    - contact_audit["available_water_flux_kg_m2_s"]
                                )
                                if deficit * duration * area / scales[0] > 1e-10:
                                    raise RuntimeError("wet-contact resolved supply mismatch")
                        except Exception as root_error:
                            result = None
                            last_error = root_error
                            record["attempts"].append(
                                {
                                    "accepted": False,
                                    "time_s": current.time_s,
                                    "attempt_dt_s": duration,
                                    "surface_phase_candidate": phase,
                                    "error": str(root_error),
                                    "traceback": traceback.format_exc(),
                                    "optimizer_debug": getattr(root_error, "optimizer_debug", None),
                                    "path_audit": dict(path_audit),
                                }
                            )
                            if present:
                                result = dry.advance(
                                    current,
                                    config,
                                    surface,
                                    carry,
                                    duration,
                                    boundary,
                                    zero_film_event=True,
                                    boundary_factory=lambda d: dry.boundary_driver.bath_boundary(
                                        bath, current.grid.R, current.time_s + d
                                    ),
                                )
                                event = True
                                if not 0 < result[0].time_s - current.time_s <= duration:
                                    raise ValueError(
                                        "depletion root lies outside attempted interval"
                                    )
                        if result is not None:
                            limiter_audit = limiter.last_audit
                            break
                if result is None:
                    raise last_error
                after, audit, next_carry, report = result
                report["fast_path_enabled"] = wet_contact and args.factor > 1
                report["fast_path_factor"] = args.factor if wet_contact else 1.0
                report["zero_inventory_wet_contact"] = wet_contact and not present
                report["surface_liquid_supply_limiter"] = limiter_audit
                if wet_contact and not present:
                    report["external_water_phase"] = "wet contact with resolved supply"
                report["native_config_branch_is_coefficient_metadata_only"] = False
                candidate = row(after, audit, next_carry)
                dt_actual = after.time_s - current.time_s
                flows = [
                    audit.external_water_flux_kg_m2_s,
                    audit.external_hexane_flux_kg_m2_s,
                    audit.external_energy_flux_w_m2,
                ]
                next_integrated = [s + dt_actual * area * f for s, f in zip(integrated, flows)]
                defects = [
                    abs(v - i + s) / scale
                    for v, i, s, scale in zip(
                        candidate["combined_W_H_U"],
                        initial["combined_W_H_U"],
                        next_integrated,
                        scales,
                    )
                ]
                if max(defects) > 1e-10:
                    raise RuntimeError("independent cumulative exterior-flux ledger refused")
                integrated = next_integrated
                maximum = [max(a, b) for a, b in zip(maximum, defects)]
                report["independent_cumulative_W_H_U"] = defects
                record["attempts"].append(report)
                record["history"].append(candidate)
                current, surface, carry = after, audit, next_carry
                accepted = True
                args.out.parent.mkdir(parents=True, exist_ok=True)
                args.out.with_suffix(".progress.pkl").write_bytes(
                    pickle.dumps(
                        {
                            **checkpoint,
                            "current": current,
                            "surface": surface,
                            "carry": carry,
                            "source_sha256": hashes,
                            "result_path": str(args.out),
                            "fast_wetting_study_factor": args.factor,
                            "radial_fast_wetting_force_version": "actual_retained_activity_v2",
                        },
                        protocol=5,
                    )
                )
                if index % 5 == 0 or event:
                    args.out.with_suffix(".progress.json").write_text(
                        json.dumps(
                            {
                                **record,
                                "final": candidate,
                                "status": "RUNNING; not final evidence",
                            },
                            indent=2,
                        ),
                        encoding="utf-8",
                    )
                if index % 5 == 0 or event:
                    print(
                        json.dumps(
                            {
                                "step": index + 1,
                                "time_s": current.time_s,
                                "W": candidate["particle_water_kg_kg_dry"],
                                "L": candidate["external_water_kg_kg_dry"],
                                "event": event,
                                "phase": phase,
                            }
                        ),
                        flush=True,
                    )
                break
            except Exception as exc:
                record["attempts"].append(
                    {
                        "accepted": False,
                        "time_s": current.time_s,
                        "attempt_dt_s": duration,
                        "error": str(exc),
                        "traceback": traceback.format_exc(),
                    }
                )
                print(
                    json.dumps(
                        {
                            "rejected_at_s": current.time_s,
                            "attempt_dt_s": duration,
                            "error": str(exc),
                        }
                    ),
                    flush=True,
                )
        if not accepted:
            fatal = record["attempts"][-1]["error"]
            break
    record.update(
        final=row(current, surface, carry),
        requested_window_reached=bool(current.time_s >= args.until),
        bath_exit_reached=bool(current.time_s >= bath["waypoints"][-1]["t_s"]),
        error=fatal,
        wall_s=time.perf_counter() - started,
        independent_maximum_normalized_W_H_U=maximum,
        independent_physical_error_bounds_per_kg_dry_W_H_U=[
            v * s / mass for v, s in zip(maximum, scales)
        ],
        source_sha256_at_exit={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)
        },
    )
    record["sources_unchanged_during_run"] = record["source_sha256_at_exit"] == hashes
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.with_suffix(".accepted.pkl").write_bytes(
        pickle.dumps(
            {
                **checkpoint,
                "current": current,
                "surface": surface,
                "carry": carry,
                "result_path": str(args.out),
                "fast_wetting_study_factor": args.factor,
                "radial_fast_wetting_force_version": "actual_retained_activity_v2",
                "source_sha256": hashes,
            },
            protocol=5,
        )
    )
    args.out.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "bath_exit_reached": record["bath_exit_reached"],
                "final": record["final"],
                "error": fatal,
                "wall_s": record["wall_s"],
            }
        ),
        flush=True,
    )
    if fatal:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
