"""Kernel-envelope STEAM_QUENCH_ARRIVAL: the stamped condensation-to-saturation arrival event (armed only).

The lead's decision of 2026-11-03 under the owner's direction of the same day (a study-lane engineering
knob of the tower lane, Class B, built on the lead's authority as a stamped interim knob), on the
diagnosis docs/GT_PS2_DTDC_COHORT_TOWER_NEW_CELL_ARRIVAL_REFUSAL_DIAGNOSIS_2026-11-03.md; record
docs/GT_PS2_UNIVERSAL_KERNEL_STEAM_QUENCH_ARRIVAL_REGIME_BUILD_2026-11-03.md.  Imported only by
particle_postcore_solver.advance_dry_phases when the knob is armed; nothing on the knob-absent path
imports it.

The physics.  A particle colder than the saturation temperature of water at its pressure, T_sat(P),
that meets superheated steam carrying (almost) no n-hexane first condenses steam on itself until its
surface reaches T_sat(P), the first period of drying in superheated steam.  The kernel's surface laws
hold a water-saturated surface vapour whose n-hexane makes up the rest of the pressure, so below T_sat
that vapour must carry n-hexane (31 mol% at 363 K) and against a hexane-free gas neither the film law
(|beta| <= 3) nor the gas-only chart has a root: the diagnosis's refusal.

What it does, at the start of every dry advance while armed:

* **The trigger** (``trigger``): a binary dry state whose stored surface temperature lies below the
  kernel's own T_sat at the state's pressure (props.water.saturation_temperature), under a presented
  gas hotter than T_sat (superheated steam: the condensation can heat the surface to T_sat) whose n-hexane lies below the declared threshold at both ends of the span (one threshold for
  film-bearing and one for film-free particles, the diagnosis's measured grid) and whose water
  activity on the kernel's own saturation law, (1 - y) P / p_sat(T_gas), lies below 1 at both ends
  (superheated: it does not condense by its own law).
* **The event** (``quench_event``): steam from the presented gas condenses onto the particle until the
  declared heated depth, the outer cell and the exterior film, stands at T_sat(P): the outer cell's
  temperature is set to T_sat(P) (only ever raised) with its n-hexane inventory held (its pore-gas y
  solved by bisection on the kernel's own equilibrium so that the cell's n-hexane is unchanged; the
  bisection's last-bit residual booked exactly), the water the cell's sorption gives up or takes up at
  that state going to or from the film, and the film (re-created when the particle is film-free) at
  T_sat(P) on the kernel's own film caloric (props.water liquid internal energy, as the dry stage's
  carry).  The condensed mass m_c follows from the energy balance in closed form:
  m_c h_v(T_gas, P) = dU_cell + F1 u_l(T_sat) - E0 with F1 = F0 + m_c - dW_cell, h_v the props.water
  vapour enthalpy of the presented gas.  The combined (particle plus film) storage change is booked as
  one inward stream of the advance (W, H, U) in the stream integral and the original-feed ledger, as
  the n-hexane-exhausted regime's events are; the particle state's own cumulative ledger carries the
  cell's change as boundary transfer to the film.  No clamp, no presented-gas change, no gate, law,
  tolerance or allowance moved; the binary liquid set then advances the event's state under the cell's
  own gas from the stored surface at T_sat(P).

physically_qualifying: false.
"""

from __future__ import annotations

from dataclasses import replace
import math

import experiment_core2_dt_dry as dry
import particle_kernel_envelope as envelope

core = dry.core
ct = dry.ct
cp = dry.cp

RECORD = envelope.STEAM_QUENCH_ARRIVAL_RECORD
REGIME = "steam_quench_arrival"
EVENT = "steam_quench_arrival_condensation"
#: the bisection of the outer cell's pore-gas y (n-hexane held): iterations on a bracket of width 1
BISECTION_ITERATIONS = 200
#: labelled measurement variants of the build (module attributes; the declared rule is the default)
HEATED_DEPTH = "outer_cell_held_hexane"
SEED_RULE = "t_sat"


class SteamQuenchArrivalRefused(RuntimeError):
    """The event's state is not admissible (a negative film, a cell state the chart refuses)."""


def _presented_y(boundary):
    bath = getattr(boundary, "dew_bath_y_hexane", None)
    return float(boundary.y_hexane if bath is None else bath)


def gas_water_activity(boundary, pressure_pa):
    """The presented gas's water activity on the kernel's saturation law, (1 - y) P / p_sat(T)."""
    t = float(boundary.temperature_k)
    if t >= core.wa.TC:
        return 0.0
    return (1.0 - _presented_y(boundary)) * pressure_pa / core.wa.saturation_pressure(t)


def _film_kg(carry, area):
    return area * carry["external_water_loading_kg_m2"]


def trigger(before, carry, surface, boundary_start, boundary_end):
    """The readings when the trigger holds, else None."""
    if not isinstance(before, ct.FullyDryTransportState):
        return None
    pressure = before.config.pressure_pa
    t_sat = core.wa.saturation_temperature(pressure)
    t_s = float(surface.temperature_k)
    if not t_s < t_sat:
        return None
    film = float(carry["external_water_loading_kg_m2"])
    threshold = (
        envelope.STEAM_QUENCH_GAS_Y_HEXANE_FILM
        if film > 0.0
        else envelope.STEAM_QUENCH_GAS_Y_HEXANE_FILM_FREE
    )
    gas_t = (float(boundary_start.temperature_k), float(boundary_end.temperature_k))
    if not all(t > t_sat for t in gas_t):
        return None
    ys = (_presented_y(boundary_start), _presented_y(boundary_end))
    if not all(y < threshold for y in ys):
        return None
    activities = tuple(gas_water_activity(b, pressure) for b in (boundary_start, boundary_end))
    if not all(a < 1.0 for a in activities):
        return None
    return {
        "surface_temperature_k": t_s,
        "t_sat_k": t_sat,
        "pressure_pa": pressure,
        "film_kg_m2": film,
        "threshold_y_hexane": threshold,
        "gas_y_hexane_start_end": list(ys),
        "gas_water_activity_start_end": list(activities),
        "gas_temperature_k_start_end": list(gas_t),
    }


def _inventories(grid, temperatures, y, config):
    return ct._inventories(grid, ct._equilibrium_states(temperatures, y, config))


def _outer_cell_at(before, t_new):
    """The outer cell at t_new with its n-hexane inventory held: (temperatures, y) of the new state."""
    n = before.grid.n
    temps = list(before.temperatures_k)
    ys = list(before.y_hexane)
    if t_new <= temps[-1] or HEATED_DEPTH == "film":
        return tuple(temps), tuple(ys), 0
    if HEATED_DEPTH == "outer_cell_held_y":
        temps[-1] = t_new
        return tuple(temps), tuple(ys), 0
    target = _inventories(before.grid, before.temperatures_k, before.y_hexane, before.config)
    target_h = target.hexane_mol[-1]
    temps[-1] = t_new

    def hexane(y):
        trial = list(ys)
        trial[-1] = y
        return _inventories(before.grid, temps, trial, before.config).hexane_mol[n - 1]

    lo, hi = ys[-1], ys[-1]
    # bracket: heating at fixed y lowers the sorbed n-hexane, so the held inventory needs a larger y
    step = 1e-3
    while hexane(hi) < target_h:
        lo = hi
        hi = min(hi + step, 1.0 - 1e-12)
        step *= 2.0
        if hi >= 1.0 - 1e-12 and hexane(hi) < target_h:
            raise SteamQuenchArrivalRefused("no outer-cell y holds the n-hexane inventory")
    count = 0
    for count in range(1, BISECTION_ITERATIONS + 1):
        mid = 0.5 * (lo + hi)
        if mid in (lo, hi):
            break
        if hexane(mid) < target_h:
            lo = mid
        else:
            hi = mid
    ys[-1] = lo if abs(hexane(lo) - target_h) <= abs(hexane(hi) - target_h) else hi
    return tuple(temps), tuple(ys), count


def zero_flux_temperature(before, boundary, t_sat):
    """The temperature above t_sat where the kernel's water-saturated interface n-hexane mass
    fraction (cp.decode_gas_only_y at the film node's coordinate -32) equals the presented gas's."""
    pressure, pore = before.config.pressure_pa, before.config.pore
    bulk = core._hexane_mass_fraction(float(boundary.y_hexane))

    def excess(t):
        return core._hexane_mass_fraction(cp.decode_gas_only_y(t, pressure, -32.0, pore)) - bulk

    lo, hi = t_sat, t_sat + 1.0
    while excess(hi) > 0.0:
        lo, hi = hi, hi + 1.0
    for _ in range(BISECTION_ITERATIONS):
        mid = 0.5 * (lo + hi)
        if mid in (lo, hi):
            break
        if excess(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return hi


def _differences(old, new):
    return (
        math.fsum(b - a for a, b in zip(old.water_mol, new.water_mol)),
        math.fsum(b - a for a, b in zip(old.hexane_mol, new.hexane_mol)),
        math.fsum(b - a for a, b in zip(old.energy_j, new.energy_j)),
    )


def quench_event(before, carry, surface, boundary, readings):
    """(after state, carry, seed surface, event): the condensation-to-T_sat(P) event."""
    pressure = readings["pressure_pa"]
    t_sat = readings["t_sat_k"]
    grid = before.grid
    area = grid.areas[-1]
    old = _inventories(grid, before.temperatures_k, before.y_hexane, before.config)
    temps, ys, iterations = _outer_cell_at(before, t_sat)
    new = _inventories(grid, temps, ys, before.config)
    d_w_mol, d_h_mol, d_u = _differences(old, new)
    d_w = d_w_mol * core.wa.M
    f0 = _film_kg(carry, area)
    e0 = float(carry["external_water_energy_j"])
    u_l = core.wa.state_Tp(t_sat, pressure, "liquid").u_mass
    t_gas = float(boundary.temperature_k)
    h_v = core.wa.state_Tp(t_gas, pressure, "vapor").h_mass
    condensed = (d_u + (f0 - d_w) * u_l - e0) / (h_v - u_l)
    f1 = f0 + condensed - d_w
    if not (math.isfinite(condensed) and condensed >= 0.0 and f1 > 0.0):
        raise SteamQuenchArrivalRefused(
            f"the event's film is not admissible (condensed {condensed!r} kg, film {f1!r} kg)"
        )
    e1 = f1 * u_l
    # the outward stream of the advance: minus the combined (particle plus film) storage change
    stream = [
        -math.fsum((d_w, f1, -f0)),
        -d_h_mol * core.hx.M,
        -math.fsum((d_u, e1, -e0)),
    ]
    after = replace(
        before,
        temperatures_k=temps,
        y_hexane=ys,
        cumulative_boundary_water_out_mol=before.cumulative_boundary_water_out_mol - d_w_mol,
        cumulative_boundary_hexane_out_mol=before.cumulative_boundary_hexane_out_mol - d_h_mol,
        cumulative_boundary_energy_out_j=before.cumulative_boundary_energy_out_j - d_u,
        cumulative_water_inventory_change_mol=math.fsum(
            (before.cumulative_water_inventory_change_mol, d_w_mol)
        ),
        cumulative_hexane_inventory_change_mol=math.fsum(
            (before.cumulative_hexane_inventory_change_mol, d_h_mol)
        ),
        cumulative_energy_inventory_change_j=math.fsum(
            (before.cumulative_energy_inventory_change_j, d_u)
        ),
        cumulative_absolute_water_transfer_mol=before.cumulative_absolute_water_transfer_mol
        + abs(d_w_mol),
        cumulative_absolute_hexane_transfer_mol=before.cumulative_absolute_hexane_transfer_mol
        + abs(d_h_mol),
        cumulative_absolute_energy_transfer_j=before.cumulative_absolute_energy_transfer_j
        + abs(d_u),
    )
    next_carry = {
        **carry,
        "external_water_mass_kg": f1,
        "external_water_loading_kg_m2": f1 / area,
        "external_water_energy_j": e1,
        "external_water_temperature_k": t_sat,
    }
    seed_t = t_sat if SEED_RULE == "t_sat" else zero_flux_temperature(before, boundary, t_sat)
    seed = replace(surface, temperature_k=seed_t, liquid_water_kg_m2=f1 / area)
    event = {
        "kind": EVENT,
        "clock_s": before.time_s,
        "stream_W_H_U": stream,
        "condensed_water_kg": condensed,
        "steam_enthalpy_j_kg": h_v,
        "film_internal_energy_j_kg": u_l,
        "film_kg_before_after": [f0, f1],
        "film_energy_j_before_after": [e0, e1],
        "surface_temperature_k_before_after": [readings["surface_temperature_k"], t_sat],
        "outer_cell_temperature_k_before_after": [
            float(before.temperatures_k[-1]),
            float(temps[-1]),
        ],
        "outer_cell_y_hexane_before_after": [float(before.y_hexane[-1]), float(ys[-1])],
        "outer_cell_bisection_iterations": iterations,
        "cell_change_W_H_U": [d_w, d_h_mol * core.hx.M, d_u],
        "hexane_residual_booked_kg": d_h_mol * core.hx.M,
        "heated_depth": HEATED_DEPTH,
        "seed_rule": SEED_RULE,
        "seed_surface_temperature_k": seed_t,
        "readings": readings,
        "record": RECORD,
        "physically_qualifying": False,
    }
    return after, next_carry, seed, event


def regime_stamp(event):
    return {
        "regime": REGIME,
        "kind": EVENT,
        "events": [event["kind"]],
        "record": RECORD,
        "physically_qualifying": False,
    }


def route(before, carry, surface, duration, boundary_factory, *, binary):
    """The armed dispatch at the start of a dry advance: None (no event; the caller continues
    unchanged) or the binary liquid set's DryPhaseAdvance from the event's state, the event stamped.
    """
    start = boundary_factory(before.time_s)
    end = boundary_factory(before.time_s + duration)
    readings = trigger(before, carry, surface, start, end)
    if readings is None:
        return None
    after, next_carry, seed, event = quench_event(before, carry, surface, start, readings)
    outcome = binary(after, next_carry, seed, "liquid")
    return replace(
        outcome,
        native_report={
            **outcome.native_report,
            "kernel_envelope_steam_quench_arrival": regime_stamp(event),
        },
        regime_events=(event, *outcome.regime_events),
    )


def declared_row():
    return envelope.steam_quench_arrival_declared_row()
