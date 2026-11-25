"""Pre-desolventizer boundary items of the universal-kernel program (2026-10-19).

Owner ruling GT_PS2_OWNER_RULING_UNIVERSAL_KERNEL_ENVELOPE_2026-10-19, item 3 (the conduction
boundary and the stalled-front mode) and item 2's interim rule (the projection of a
water-supersaturated bath gas onto the kernel's dew line). Each behaviour is behind its own
knob of particle_kernel_envelope (CONDUCTION_BOUNDARY, STALLED_FRONT, GAS_DEW_LINE_PROJECTION),
off by default; nothing here runs unless a caller armed it. physically_qualifying is false.

* Conduction boundary: the bath carries, per waypoint, ``jacket_mode`` ("conduction" in the
  tower's jacket-heated pre-desolventizer trays, "gas_contact" elsewhere),
  ``jacket_wall_temperature_K`` and ``jacket_contact_conductance_w_m2_k`` (per unit outer
  particle area), taken per interval (the right waypoint's value on (t_{k-1}, t_k], as the
  tower books the jacket heat per interval). In a conduction interval the surface energy
  balance's sensible gas term h_film A (T_s - T_b) is replaced by h_c (T_s - T_wall) (heat
  leaving the surface positive), the tower's floor-to-bottom-solid conductance shared by dry
  mass with its gas contact fraction 0 (no sensible gas-to-meal exchange); the vapour leaving
  or condensing keeps the kernel's own film mass law against the bath gas.
* Dew-line projection: a bath gas with y_h below the lower end of the kernel's gas-only
  interval at (T_b, P) (water activity above one: water-supersaturated) is presented to the
  kernel at the dew line y_dew (as presented: the first representable gas state strictly inside
  the open gas-only chart, found by an ulp ladder from the bisected end and stamped with its
  offset and water activity; a representational choice, no physical offset); the water the projection removes per mole of gas,
  1 - y_h / y_dew, times the gas molar density P / (R T_b) and the water molar mass, carried
  to the surface at the film's own scaled mass-transfer coefficient, is booked as a water
  deposition flux onto the external film (the film node's water row) with the vapour's partial
  enthalpy at the bath state (the node's energy row). A film-free surface cannot take it
  (typed refusal).
* Stalled front (core2_water_film_node.advance_water_node(stalled_front=True), driven from
  qualify_sliding_front_history): the native interior-stall construction at the film node.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

import finite_film_core2_study as study

core = study.core
R_GAS = 8.31446261815324


class DewLineCondensateOnFilmFreeSurface(RuntimeError):
    """The dew-line projection books condensate, but the surface carries no film."""


class NoGasOnlyIntervalAtBath(RuntimeError):
    """The bath (T, P) has no gas-only interval at all (below the co-boiling root)."""


@dataclass(frozen=True)
class KernelEnvelopeBoundary(study.FiniteFilmStudyBoundary):
    """The study boundary plus the armed pre-desolventizer items (floats only, hashable)."""

    jacket_wall_temperature_k: float | None = None
    jacket_contact_conductance_w_m2_k: float | None = None
    dew_bath_y_hexane: float | None = None
    dew_removed_mol_per_mol: float = 0.0
    dew_deposition_kg_m2_s: float = 0.0
    dew_water_vapour_enthalpy_j_kg: float = 0.0
    dew_interior_offset_ulps: int | None = None
    dew_presented_water_activity: float | None = None
    dew_bisected_y_hexane: float | None = None


def bath_jacket(bath, t_s):
    """(T_wall, h_c) of the interval holding t_s, or None (no field / gas contact)."""
    waypoints = bath["waypoints"]
    carried = ["jacket_mode" in w for w in waypoints]
    if not any(carried):
        return None
    if not all(carried):
        raise ValueError("bath jacket field must be carried by every waypoint or by none")
    chosen = waypoints[-1]
    if t_s <= float(waypoints[0]["t_s"]):
        chosen = waypoints[0]
    else:
        for w in waypoints[1:]:
            if t_s <= float(w["t_s"]):
                chosen = w
                break
    mode = chosen["jacket_mode"]
    if mode == "gas_contact":
        return None
    if mode != "conduction":
        raise ValueError(f"unknown jacket mode {mode!r}")
    wall, h = float(chosen["jacket_wall_temperature_K"]), float(
        chosen["jacket_contact_conductance_w_m2_k"]
    )
    if not (math.isfinite(wall) and wall > 0 and math.isfinite(h) and h >= 0):
        raise ValueError("jacket wall temperature and conductance must be finite, T > 0, h >= 0")
    return wall, h


# The ulp ladder's search limit (a representational limit, not a physical offset): evidence
# 167 measured the newborn's dependence on the offset as immaterial up to 4,096 ulps of y_h
# (at most 4e-11 K of surface temperature); the first interior point sat at 16 to 64 ulps.
DEW_INTERIOR_MAXIMUM_ULPS = 4096


def _dew_line_interior(temperature_k, pressure_pa, k_wh):
    """The ulp ladder at the bath (evidence 167 section 3; the lead's record of 2026-10-20,
    section 3 item 3): from the bisected lower end of the kernel's gas-only interval
    (water saturation, a_w = 1) at (T, P), the smallest number of binary64 steps of y_h
    towards the interior at which the kernel's own water activity rounds strictly below 1,
    so that the open gas-only chart's strict check 1 - a_w > 0 holds. A representational
    choice: it declares no physical offset from the dew line."""
    params = replace(core.cp.CoupledPoreParams(), k_wh=k_wh)
    try:
        interval = core.cp.gas_only_composition_interval(temperature_k, pressure_pa, params)
    except core.cp.CoupledPoreTopologyError as error:
        raise NoGasOnlyIntervalAtBath(
            f"the bath gas at {temperature_k!r} K and {pressure_pa!r} Pa has no gas-only "
            f"interval on the kernel's equilibrium ({error}); no dew line to project onto"
        ) from error
    end = float(interval.y_hexane_bounds[0])
    y = end
    for offset in range(DEW_INTERIOR_MAXIMUM_ULPS + 1):
        try:
            state = core.cp.evaluate_equilibrium(temperature_k, pressure_pa, y, params)
        except core.cp.CoupledPoreTopologyError:
            state = None
        if state is not None and state.water_activity < 1.0:
            return y, interval, offset, float(state.water_activity), end
        y = math.nextafter(y, 1.0)
    raise ValueError(
        "no gas state strictly inside the open gas-only chart within "
        f"{DEW_INTERIOR_MAXIMUM_ULPS} ulps of the dew line at {temperature_k!r} K, "
        f"{pressure_pa!r} Pa"
    )


_DEW_LINE_MEMO = {}


def dew_line_point(temperature_k, pressure_pa, k_wh):
    """(y_presented, interval, offset in ulps from the bisected end, a_w at y_presented,
    the bisected end), memoized per exact (T, P, k_wh): the ladder is a pure function."""
    key = (float(temperature_k), float(pressure_pa), float(k_wh))
    if key not in _DEW_LINE_MEMO:
        _DEW_LINE_MEMO[key] = _dew_line_interior(temperature_k, pressure_pa, k_wh)
    return _DEW_LINE_MEMO[key]


def dew_line_y_hexane(temperature_k, pressure_pa, k_wh):
    """The dew line as presented to the kernel: the first representable gas state strictly
    inside the open gas-only chart (dew_line_point). Returns (y, interval)."""
    y, interval, _, _, _ = dew_line_point(temperature_k, pressure_pa, k_wh)
    return y, interval


def dew_projection(temperature_k, y_hexane, pressure_pa, film_boundary):
    """None inside the gas-only interval's water side; else the projection's quantities."""
    k_wh = film_boundary.binary_gas_interaction_k_wh
    y_dew, _, offset_ulps, presented_activity, bisected_end = dew_line_point(
        temperature_k, pressure_pa, k_wh
    )
    if not y_hexane < y_dew:
        return None
    removed = 1.0 - y_hexane / y_dew
    mass_transfer = (
        film_boundary.stress_multiplier * film_boundary.correlated_pair.binary_mass_transfer_m_s
    )
    density = pressure_pa / (R_GAS * temperature_k)
    deposition = mass_transfer * density * removed * core.wa.M
    params = replace(core.cp.CoupledPoreParams(), k_wh=k_wh)
    vapour = core.cp.evaluate_equilibrium(temperature_k, pressure_pa, y_dew, params)
    return {
        "y_dew": y_dew,
        "removed_mol_per_mol": removed,
        "deposition_kg_m2_s": deposition,
        "water_vapour_enthalpy_j_kg": vapour.water_gas_partial_enthalpy_j_mol / core.wa.M,
        "mass_transfer_m_s": mass_transfer,
        "gas_molar_density_mol_m3": density,
        "interior_offset_ulps": offset_ulps,
        "presented_water_activity": presented_activity,
        "bisected_dew_line_y_hexane": bisected_end,
    }


def envelope_boundary(base, *, jacket=None, projection=None):
    """The armed boundary: base's values with the jacket and the projected gas."""
    fields = dict(
        temperature_k=base.temperature_k,
        y_hexane=base.y_hexane,
        label=base.label,
        film_boundary=base.film_boundary,
        accepted_beta_limit=base.accepted_beta_limit,
        actual_surface_path=base.actual_surface_path,
        numerical_cap_seam=base.numerical_cap_seam,
    )
    extra = {}
    if jacket is not None:
        extra.update(
            jacket_wall_temperature_k=jacket[0], jacket_contact_conductance_w_m2_k=jacket[1]
        )
    if projection is not None:
        film = replace(base.film_boundary, y_hexane=projection["y_dew"])
        fields.update(y_hexane=projection["y_dew"], film_boundary=film)
        extra.update(
            dew_bath_y_hexane=base.y_hexane,
            dew_removed_mol_per_mol=projection["removed_mol_per_mol"],
            dew_deposition_kg_m2_s=projection["deposition_kg_m2_s"],
            dew_water_vapour_enthalpy_j_kg=projection["water_vapour_enthalpy_j_kg"],
            dew_interior_offset_ulps=projection.get("interior_offset_ulps"),
            dew_presented_water_activity=projection.get("presented_water_activity"),
            dew_bisected_y_hexane=projection.get("bisected_dew_line_y_hexane"),
        )
    return KernelEnvelopeBoundary(**fields, **extra)


def surface_gas_heat_w_m2(boundary, film_heat_w_m2, surface_temperature_k):
    """Heat leaving the surface to the outside: the film's term, or the jacket's."""
    wall = getattr(boundary, "jacket_wall_temperature_k", None)
    if wall is None:
        return film_heat_w_m2
    return boundary.jacket_contact_conductance_w_m2_k * (surface_temperature_k - wall)


# ---- per-step stamps ---------------------------------------------------------------------


def step_stamp(boundary, surface_temperature_k, elapsed_s):
    """The armed items at an accepted step's end-clock boundary; None when neither acts."""
    if not isinstance(boundary, KernelEnvelopeBoundary):
        return None
    stamp = {}
    if boundary.jacket_wall_temperature_k is not None:
        h = boundary.jacket_contact_conductance_w_m2_k
        stamp["jacket"] = {
            "wall_temperature_k": boundary.jacket_wall_temperature_k,
            "contact_conductance_w_m2_k": h,
            "heat_into_surface_w_m2": h
            * (boundary.jacket_wall_temperature_k - surface_temperature_k),
            "surface_temperature_k": surface_temperature_k,
        }
    if boundary.dew_bath_y_hexane is not None:
        stamp["dew_projection"] = {
            "bath_y_hexane": boundary.dew_bath_y_hexane,
            "projected_y_hexane": boundary.y_hexane,
            "removed_mol_per_mol": boundary.dew_removed_mol_per_mol,
            "deposition_kg_m2_s": boundary.dew_deposition_kg_m2_s,
            "deposited_kg_m2": boundary.dew_deposition_kg_m2_s * elapsed_s,
            "interior_offset_ulps": boundary.dew_interior_offset_ulps,
            "presented_water_activity": boundary.dew_presented_water_activity,
            "bisected_dew_line_y_hexane": boundary.dew_bisected_y_hexane,
        }
    return stamp or None


def empty_summary():
    return {
        "conduction_steps": 0,
        "projection_steps": 0,
        "projected_condensate_kg_m2": 0.0,
        "first_projection_time_s": None,
    }


def add_to_summary(summary, stamp, time_s):
    if stamp.get("jacket") is not None:
        summary["conduction_steps"] += 1
    dew = stamp.get("dew_projection")
    if dew is not None:
        summary["projection_steps"] += 1
        summary["projected_condensate_kg_m2"] += dew["deposited_kg_m2"]
        if summary["first_projection_time_s"] is None:
            summary["first_projection_time_s"] = time_s
    return summary
