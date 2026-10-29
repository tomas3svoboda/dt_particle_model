"""The declined advance of a held front, stamped as a typed materiality row (2026-10-20).

Kernel-envelope STALL_FILM_DEPLETION_EVENT (particle_kernel_envelope; the lead's record
GT_PS2_UNIVERSAL_KERNEL_PD_ADVANCE_DEMAND_2026-10-20 section 2). A stalled step holds the front
(PHY-014: the front never advances) while the jump rows may demand an advance. The hexane that
the held front declines to condense is not booked anywhere in the accepted history; this module
measures it and stamps it, so that it is never silently dropped.

The measurement is evidence pd_advance_demand_170's construction
(particle_kernel_counterfactual_front_velocity): from the accepted step's start state, over the
step's duration, the front a literal at its accepted position and a SIGNED free velocity v in
the Rankine-Hugoniot jump terms only, q = A_Gamma v, every jump row kept; solved at the
unchanged search target, evaluated in a scratch overlay of the film node's advance and never
adopted (nothing it returns enters any state, ledger or history). The declined hexane rate is
q (C_h,wet - C_h,dry) (positive for an advance demand), C_h,wet = rho_dm_p X_hist / M_h.

Per stalled step: the rate at the step's start (the previous stalled step's end rate when that
step was stalled from the same state, else a declared probe of DEMAND_START_PROBE_S) and at its
end (the step's own duration), integrated by the trapezoid; per kg of dry meal (the particle's
dry mass rho_dm_p times its grid volume); its cumulative value per packet; its ratio to the
particle's hexane at the step's end. A refused held-front event's window (its located duration,
from the event's start state) is stamped on the event's refusal record by window_row, from a
declared ladder of rungs. Report only. physically_qualifying false.
"""

from __future__ import annotations

import math

import particle_kernel_counterfactual_front_velocity as cf
import particle_kernel_envelope as envelope

# The declared probe duration of the start-of-step rate when no earlier stalled rate from the
# same state exists: evidence 170's shortest ladder rung (at the stall it reads +0.0002 um/s).
# A measurement instrument's setting, not a model value.
DEMAND_START_PROBE_S = 1e-5
# The declared rungs of a refused held-front event's window, as fractions of its located
# duration (a measurement instrument's setting): evidence 170's held-front ladder converged to
# within about 0.05 s of the film's run-out and refused beyond it.
WINDOW_LADDER_FRACTIONS = (0.25, 0.5, 0.75, 0.9, 0.95, 1.0)


def particle_dry_mass_kg(transport):
    """The particle's dry mass (kg): rho_dm_p times the grid volume (primitive_row's mass)."""
    grid = transport.geometry.master_grid
    return transport.config.dry.pore.rho_dm_p * math.fsum(grid.volumes)


def counterfactual_demand(
    engine,
    current,
    carry,
    surface,
    bath,
    duration,
    *,
    initial_front_z,
    boundary_factory,
    gas_only_surface,
    options,
    search_tolerance,
):
    """One counterfactual evaluation from current over duration; a dict, never raises.

    The film node's advance is swapped for the counterfactual one inside this call only and
    restored in every case. options are the dispatcher's envelope options for the step (a fresh
    scratch audit dictionary is used).
    """
    node = engine.driver.node
    original = node.advance_water_node
    node.advance_water_node = cf.advance_water_node
    cf.LAST.clear()
    cf.LAST["velocity_seed_um_s"] = 0.0
    scratch = {k: v for k, v in options.items() if k != "envelope_audit"}
    scratch["envelope_audit"] = {}
    try:
        evidence = engine.advance_partial(
            current,
            carry,
            surface,
            bath,
            duration,
            initial_front_z=initial_front_z,
            search_tolerance=search_tolerance,
            boundary_factory=boundary_factory,
            gas_only_surface=gas_only_surface,
            stalled_front=True,
            **scratch,
        )
    except Exception as error:  # a measurement: its failure is reported, never raised
        return {"evaluated": False, "duration_s": duration, "error": str(error)[:600]}
    finally:
        node.advance_water_node = original
        cf.LAST.clear()
    jumps = evidence["jumps"]
    q = evidence["q_m3_s"]
    rate_mol_s = -q * jumps["jump_hexane_mol_m3"]
    return {
        "evaluated": True,
        "duration_s": evidence["dt_s"],
        "converged_at_gate": bool(
            evidence["solver_success"] and evidence["max_scaled_residual"] <= evidence["gate"]
        ),
        "max_scaled_residual": evidence["max_scaled_residual"],
        "gate": evidence["gate"],
        "velocity_um_s": evidence["velocity_um_s"],
        "q_m3_s": q,
        "declined_hexane_rate_kg_s": rate_mol_s * cf.core.hx.M,
        "water_jump_term_mol_s": q * jumps["jump_water_mol_m3"],
        "energy_jump_term_w": q * jumps["jump_energy_j_m3"],
        "film_at_end_kg_m2": evidence["surface"]["liquid_water_kg_m2"],
        "interface_temperature_k": evidence["interface_temperature_k"],
        "residual_calls": evidence["residual_calls"],
    }


def window_row(start, rungs, located_s, dry_mass_kg, particle_hexane_kg_kg_dry):
    """The declined advance over a refused held-front event's window (report only).

    start: the rate at the window's start; rungs: counterfactual evaluations at the declared
    WINDOW_LADDER_FRACTIONS of the located duration. Evidence 170's construction: the trapezoid
    over the converged rungs, then the last two converged rates extended linearly to the located
    instant (the counterfactual refuses where its own film would go below zero, just before the
    held front's film runs out); the extrapolated share is stamped.
    """
    good = [r for r in rungs if r.get("converged_at_gate")]
    row = {
        "kind": "DECLARED_ENGINEERING_MATERIALITY_ROW",
        "quantity": "hexane the held front would decline over the refused event's window "
        "(counterfactual advance demand; nothing accepted)",
        "record": envelope.STALL_RECORD,
        "located_duration_s": located_s,
        "start": start,
        "rungs": rungs,
        "quadrature": "trapezoid over the converged rungs, linear extension of the last two "
        "rates to the located instant (evidence pd_advance_demand_170 section 3.2)",
        "declined_hexane_kg_kg_dry": None,
        "extrapolated_share": None,
        "particle_hexane_kg_kg_dry": particle_hexane_kg_kg_dry,
        "ratio_to_particle_hexane": None,
        "dry_mass_kg": dry_mass_kg,
        "evaluated": False,
        "report_only": True,
        "physically_qualifying": False,
    }
    if not (start and start.get("converged_at_gate")) or len(good) < 2:
        return row
    points = [(0.0, start["declined_hexane_rate_kg_s"])] + [
        (r["duration_s"], r["declined_hexane_rate_kg_s"]) for r in good
    ]
    integral = sum((b[0] - a[0]) * (a[1] + b[1]) / 2 for a, b in zip(points, points[1:]))
    tail = 0.0
    (t1, m1), (t2, m2) = points[-2], points[-1]
    if located_s > t2:
        m_end = m2 + (m2 - m1) / (t2 - t1) * (located_s - t2)
        tail = (located_s - t2) * (m2 + m_end) / 2
    total = (integral + tail) / dry_mass_kg
    row.update(
        declined_hexane_kg_kg_dry=total,
        extrapolated_share=tail / (integral + tail) if integral + tail != 0.0 else None,
        ratio_to_particle_hexane=total / particle_hexane_kg_kg_dry,
        last_converged_duration_s=t2,
        evaluated=True,
    )
    return row


def materiality_row(start, end, elapsed_s, dry_mass_kg, particle_hexane_kg_kg_dry, cumulative):
    """The typed row of one held-front step; cumulative is the packet's running total (kg/kg)."""
    usable = (
        start is not None
        and end is not None
        and start.get("converged_at_gate")
        and end.get("converged_at_gate")
    )
    declined = None
    if usable:
        declined = (
            0.5
            * (start["declined_hexane_rate_kg_s"] + end["declined_hexane_rate_kg_s"])
            * elapsed_s
            / dry_mass_kg
        )
        cumulative = cumulative + declined
    return {
        "kind": "DECLARED_ENGINEERING_MATERIALITY_ROW",
        "quantity": "hexane the held front declined to condense (counterfactual advance demand)",
        "record": envelope.STALL_RECORD,
        "construction": "evidence pd_advance_demand_170: pinned front, signed free velocity in "
        "the jump terms only, evaluated in a scratch overlay, never adopted",
        "elapsed_s": elapsed_s,
        "start": start,
        "end": end,
        "quadrature": "trapezoid over the step",
        "declined_hexane_kg_kg_dry": declined,
        "cumulative_declined_hexane_kg_kg_dry": cumulative,
        "particle_hexane_kg_kg_dry": particle_hexane_kg_kg_dry,
        "ratio_to_particle_hexane": (
            None if declined is None else cumulative / particle_hexane_kg_kg_dry
        ),
        "dry_mass_kg": dry_mass_kg,
        "evaluated": bool(usable),
        "report_only": True,
        "physically_qualifying": False,
    }
