"""Kernel-envelope HEXANE_EXHAUSTED_REGIME: the declared n-hexane-exhausted dry regime (armed only).

Owner ruling 2026-11-01 ("I think this 'ghost' regime for hexane free vapors is completely fine.
Let's build it."), on the diagnosis
docs/GT_PS2_UNIVERSAL_KERNEL_STRIPPED_PARTICLE_IN_SUPERHEATED_STEAM_DIAGNOSIS_2026-11-01.md; record
docs/GT_PS2_UNIVERSAL_KERNEL_HEXANE_EXHAUSTED_REGIME_BUILD_2026-11-01.md.  Imported only by
particle_postcore_solver.advance_dry_phases when the knob is armed (and by the row helpers when a
state of this module exists); nothing on the knob-absent path imports it.

What it does, at the start of every dry advance while armed:

* **Entry (the release event).**  A binary dry state holding less than
  HEXANE_EXHAUSTED_LOADING_KG_KG of n-hexane per kg dry matter, with no exterior film, under a
  presented gas whose n-hexane is at or below HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP at both ends of the
  span, is converted cell by cell to the pure-water vertex state (temperature and pressure carried;
  the vertex has no composition degree of freedom, so each cell's water follows from its temperature).
  The difference of the two storages, computed by the kernel's own caloric for every cell (water,
  n-hexane and common-datum energy), is released to the gas stream: booked as one stamped event in
  the stream integral and in the state's cumulative ledger, so every ledger closes exactly.
* **The vertex stage.**  Water and energy only (2N + 1 unknowns: N cell temperatures, N total Stefan
  fluxes, the surface temperature), y_hexane exactly 0 in every cell and every evaluation
  (dtdc_simulator.core2.particle.pure_water_pore_gas), the presented gas taken as hexane-free and
  stamped so; the same domain-aware Newton, allowance, gate units, noise-aware energy-row acceptance
  and 1e-10 ledgers as the binary film-free set.
* **Re-entry.**  A vertex state under a presented gas above the clamp (at either end of the span)
  returns to the binary stage: every cell's pore gas set to the presented gas's own n-hexane fraction,
  the sorbed n-hexane at the frozen isotherm's equilibrium with it, the n-hexane so dosed taken from
  the gas stream and booked as one stamped event; the binary dry stage then runs unchanged.  A gas
  state the binary chart refuses is the typed refusal KernelHexaneExhaustedRegimeGasHexaneReturned.

physically_qualifying: false.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import math

import numpy as np

import core2_water_film_node as node
import experiment_core2_dt_dry as dry
import particle_kernel_envelope as envelope
import particle_kernel_pd_boundary as pd_boundary
from dtdc_simulator.core2.particle import cut_extinction_event as extinction
from dtdc_simulator.core2.particle import pure_water_pore_gas as pw
from dtdc_simulator.core2.particle.numerical_work import StepBudgetExceeded

core = dry.core
ct = dry.ct
cp = dry.cp
study = dry.study

RECORD = envelope.HEXANE_EXHAUSTED_REGIME_RECORD
REGIME = "hexane_exhausted"
#: The DryPhaseAdvance phase of a vertex advance (a film-free set).
PHASE = "hexane_exhausted"
RELEASE = "hexane_exhausted_release"
REENTRY = "hexane_exhausted_reentry"


class KernelHexaneExhaustedRegimeGasHexaneReturned(RuntimeError):
    """The presented gas's n-hexane rose above the clamp at a vertex state and the binary chart
    refuses the re-entry state (a gas state neither branch admits)."""

    def __init__(self, message, row):
        super().__init__(message)
        self.row = row


class HexaneExhaustedVertexRefused(RuntimeError):
    """The vertex stage's root or ledger refused on this span (a rejected trial; the caller halves)."""


@dataclass(frozen=True)
class VertexDryTransportState:
    """The dry particle at the pure-water vertex: the binary FullyDryTransportState's fields with
    y_hexane exactly 0 in every cell (its own type, so no binary validator is asked to accept 0)."""

    grid: object
    config: object
    time_s: float
    temperatures_k: tuple
    reference_water_mol: float
    reference_hexane_mol: float
    reference_energy_j: float
    cumulative_boundary_water_out_mol: float = 0.0
    cumulative_boundary_hexane_out_mol: float = 0.0
    cumulative_boundary_energy_out_j: float = 0.0
    cumulative_water_inventory_change_mol: float = 0.0
    cumulative_hexane_inventory_change_mol: float = 0.0
    cumulative_energy_inventory_change_j: float = 0.0
    cumulative_absolute_water_transfer_mol: float = 0.0
    cumulative_absolute_hexane_transfer_mol: float = 0.0
    cumulative_absolute_energy_transfer_j: float = 0.0
    y_hexane: tuple = ()
    regime: str = "pure_water_vertex"
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self):
        n = self.grid.n
        if not math.isfinite(self.time_s) or self.time_s < 0.0:
            raise ValueError("state time must be finite and non-negative")
        if len(self.temperatures_k) != n or len(self.y_hexane) != n:
            raise ValueError("vertex state primitives must align with the stationary grid")
        if any(y != 0.0 for y in self.y_hexane):
            raise pw.PureWaterVertexError(
                "a vertex dry state holds y_hexane = 0 exactly in every cell"
            )
        low, high = ct._solver_temperature_bounds(self.config)
        if not all(math.isfinite(t) and low < t < high for t in self.temperatures_k):
            raise ValueError(
                "vertex state temperature must lie strictly inside its pressure-conditioned "
                "solver domain"
            )
        if self.regime != "pure_water_vertex":
            raise ValueError("unknown vertex regime label")


def is_vertex(state):
    return isinstance(state, VertexDryTransportState)


def vertex_states(temperatures_k, config):
    return tuple(
        pw.evaluate_pure_water_state(t, config.pressure_pa, config.pore) for t in temperatures_k
    )


def inventories(state):
    """The native _Inventories of a vertex state (ct._inventories on the vertex cell states)."""
    return ct._inventories(state.grid, vertex_states(state.temperatures_k, state.config))


def particle_inventory(state):
    """(water kg, n-hexane kg, energy J) of a vertex state, as radial_fast_wetting.inventory."""
    values = inventories(state)
    return (
        values.total_water_mol * core.wa.M,
        values.total_hexane_mol * core.hx.M,
        values.total_energy_j,
    )


def cell_retained_water(state):
    return [c.retained_water_loading for c in vertex_states(state.temperatures_k, state.config)]


def presented_gas_y(boundary):
    """The presented bath gas's n-hexane (the projection's bath value when the dew projection moved it)."""
    bath = getattr(boundary, "dew_bath_y_hexane", None)
    return float(boundary.y_hexane if bath is None else bath)


def _binary_inventories(state):
    return ct._inventories(
        state.grid, ct._equilibrium_states(state.temperatures_k, state.y_hexane, state.config)
    )


def dry_mass_kg(state):
    return state.config.pore.rho_dm_p * math.fsum(state.grid.volumes)


def binary_loading_kg_kg(state):
    return _binary_inventories(state).total_hexane_mol * core.hx.M / dry_mass_kg(state)


def _event_ledger(
    before, water_out, hexane_out, energy_out, water_change, hexane_change, energy_change
):
    """The state's cumulative ledger fields after an event booked as boundary outflow."""
    return dict(
        reference_water_mol=before.reference_water_mol,
        reference_hexane_mol=before.reference_hexane_mol,
        reference_energy_j=before.reference_energy_j,
        cumulative_boundary_water_out_mol=before.cumulative_boundary_water_out_mol + water_out,
        cumulative_boundary_hexane_out_mol=before.cumulative_boundary_hexane_out_mol + hexane_out,
        cumulative_boundary_energy_out_j=before.cumulative_boundary_energy_out_j + energy_out,
        cumulative_water_inventory_change_mol=math.fsum(
            (before.cumulative_water_inventory_change_mol, water_change)
        ),
        cumulative_hexane_inventory_change_mol=math.fsum(
            (before.cumulative_hexane_inventory_change_mol, hexane_change)
        ),
        cumulative_energy_inventory_change_j=math.fsum(
            (before.cumulative_energy_inventory_change_j, energy_change)
        ),
        cumulative_absolute_water_transfer_mol=before.cumulative_absolute_water_transfer_mol
        + abs(water_out),
        cumulative_absolute_hexane_transfer_mol=before.cumulative_absolute_hexane_transfer_mol
        + abs(hexane_out),
        cumulative_absolute_energy_transfer_j=before.cumulative_absolute_energy_transfer_j
        + abs(energy_out),
    )


def _differences(old, new):
    """Shell-differenced changes new - old (mol, mol, J), summed exactly (as _accept_step)."""
    return (
        math.fsum(b - a for a, b in zip(old.water_mol, new.water_mol)),
        math.fsum(b - a for a, b in zip(old.hexane_mol, new.hexane_mol)),
        math.fsum(b - a for a, b in zip(old.energy_j, new.energy_j)),
    )


def _gas_stamp(boundary):
    return {
        "temperature_k": float(boundary.temperature_k),
        "presented_y_hexane": presented_gas_y(boundary),
        "boundary_y_hexane": float(boundary.y_hexane),
        "clamp_y_hexane": envelope.HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP,
    }


def trigger(before, carry, boundary_start, boundary_end):
    """The entry condition (None when it does not hold, else the measured loading)."""
    if is_vertex(before) or carry["external_water_loading_kg_m2"] != 0.0:
        return None
    clamp = envelope.HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP
    if presented_gas_y(boundary_start) > clamp or presented_gas_y(boundary_end) > clamp:
        return None
    loading = binary_loading_kg_kg(before)
    if not loading < envelope.HEXANE_EXHAUSTED_LOADING_KG_KG:
        return None
    return loading


def release_event(before, surface, boundary, loading):
    """Binary state -> vertex state at the same temperatures; the storage difference released."""
    old = _binary_inventories(before)
    vertex_config = before.config
    new = ct._inventories(before.grid, vertex_states(before.temperatures_k, vertex_config))
    water_change, hexane_change, energy_change = _differences(old, new)
    after = VertexDryTransportState(
        grid=before.grid,
        config=vertex_config,
        time_s=before.time_s,
        temperatures_k=tuple(before.temperatures_k),
        y_hexane=tuple(0.0 for _ in before.temperatures_k),
        **_event_ledger(
            before,
            -water_change,
            -hexane_change,
            -energy_change,
            water_change,
            hexane_change,
            energy_change,
        ),
    )
    dry_mass = dry_mass_kg(before)
    stream = [-water_change * core.wa.M, -hexane_change * core.hx.M, -energy_change]
    cross_check = None
    try:
        surface_state = cp.evaluate_equilibrium(
            surface.temperature_k,
            vertex_config.pressure_pa,
            surface.surface_y_hexane,
            vertex_config.pore,
        )
        enthalpy = stream[1] / core.hx.M * surface_state.hexane_gas_partial_enthalpy_j_mol
        cross_check = {
            "released_hexane_times_surface_gas_partial_enthalpy_j": enthalpy,
            "booked_energy_minus_that_j": stream[2] - enthalpy,
            "surface_temperature_k": surface.temperature_k,
        }
    except (ValueError, RuntimeError) as error:
        cross_check = {"not_evaluated": str(error)[:200]}
    event = {
        "kind": RELEASE,
        "clock_s": before.time_s,
        "loading_before_kg_kg_dry": loading,
        "released_hexane_kg": stream[1],
        "released_hexane_kg_kg_dry": stream[1] / dry_mass,
        "cell_hexane_kg_kg_dry_before": [
            h * core.hx.M / (vertex_config.pore.rho_dm_p * v)
            for h, v in zip(old.hexane_mol, before.grid.volumes)
        ],
        "cell_y_hexane_before": [float(y) for y in before.y_hexane],
        "cell_temperatures_k": [float(t) for t in before.temperatures_k],
        "stream_W_H_U": stream,
        "water_booked_kg": stream[0],
        "energy_booked_j": stream[2],
        "caloric": "the storage difference of every cell at its own temperature and pressure by the "
        "kernel's own equilibrium caloric (coupled_pore's binary state at the cell's y_hexane minus "
        "the vertex state, the common datum of props.hexane H_REF and props.water); the water term is "
        "the vertex's O(y) pore-gas and sorbed-water difference at y_water = 1",
        "release_enthalpy_cross_check": cross_check,
        "gas": _gas_stamp(boundary),
        "exhaustion_loading_kg_kg": envelope.HEXANE_EXHAUSTED_LOADING_KG_KG,
        "record": RECORD,
        "physically_qualifying": False,
    }
    return after, event


def reentry_event(before, boundary):
    """Vertex state -> binary state at the presented gas's y in every cell; the dose from the stream."""
    y = presented_gas_y(boundary)
    clamp = envelope.HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP
    row = {
        "kind": REENTRY,
        "clock_s": before.time_s,
        "gas": _gas_stamp(boundary),
        "cell_temperatures_k": [float(t) for t in before.temperatures_k],
    }
    old = inventories(before)
    try:
        if not clamp < y < 1.0:
            raise ValueError("re-entry needs a presented n-hexane above the clamp")
        cells = tuple(y for _ in before.temperatures_k)
        new = ct._inventories(
            before.grid, ct._equilibrium_states(before.temperatures_k, cells, before.config)
        )
        water_change, hexane_change, energy_change = _differences(old, new)
        after = ct.FullyDryTransportState(
            grid=before.grid,
            config=before.config,
            time_s=before.time_s,
            temperatures_k=tuple(before.temperatures_k),
            y_hexane=cells,
            **_event_ledger(
                before,
                -water_change,
                -hexane_change,
                -energy_change,
                water_change,
                hexane_change,
                energy_change,
            ),
        )
    except (ValueError, RuntimeError) as error:
        row["refusal"] = str(error)[:300]
        raise KernelHexaneExhaustedRegimeGasHexaneReturned(
            "KernelHexaneExhaustedRegimeGasHexaneReturned: the presented gas's n-hexane "
            f"{y!r} rose above the clamp {clamp!r} at a vertex state and the binary chart refuses "
            f"the re-entry state: {error}",
            row,
        ) from error
    dry_mass = dry_mass_kg(before)
    stream = [-water_change * core.wa.M, -hexane_change * core.hx.M, -energy_change]
    row.update(
        dosed_hexane_kg=-stream[1],
        dosed_hexane_kg_kg_dry=-stream[1] / dry_mass,
        dosed_pore_gas_hexane_kg=math.fsum(
            v * c.pore_hexane_concentration_kg_m3
            for v, c in zip(
                before.grid.volumes,
                ct._equilibrium_states(before.temperatures_k, cells, before.config),
            )
        ),
        dosed_sorbed_hexane_kg=math.fsum(
            v * c.retained_hexane_concentration_kg_m3
            for v, c in zip(
                before.grid.volumes,
                ct._equilibrium_states(before.temperatures_k, cells, before.config),
            )
        ),
        cell_y_hexane_after=list(cells),
        stream_W_H_U=stream,
        water_booked_kg=stream[0],
        energy_booked_j=stream[2],
        caloric="the storage difference of every cell at its own temperature and pressure by the "
        "kernel's own equilibrium caloric (the vertex state minus coupled_pore's binary state at the "
        "presented gas's y_hexane, sorbed n-hexane at the frozen isotherm's equilibrium)",
        record=RECORD,
        physically_qualifying=False,
    )
    return after, row


def _energy_capacity(temperature_k, config):
    """ct._equilibrium_energy_capacity's stencil on the vertex states (scale only; the n-hexane datum
    term is zero at the vertex)."""
    t_lo, t_hi = ct._solver_temperature_bounds(config)
    difference = min(1.0e-3, 0.25 * (temperature_k - t_lo), 0.25 * (t_hi - temperature_k))
    if difference <= math.sqrt(np.finfo(float).eps) * temperature_k:
        raise ValueError("temperature is too close to the predictor stencil boundary")
    hotter = pw.evaluate_pure_water_state(
        temperature_k + difference, config.pressure_pa, config.pore
    )
    colder = pw.evaluate_pure_water_state(
        temperature_k - difference, config.pressure_pa, config.pore
    )
    derivative = (hotter.energy_density_j_m3 - colder.energy_density_j_m3) / (2.0 * difference)
    if not math.isfinite(derivative) or derivative == 0.0:
        raise ValueError("equilibrium energy storage derivative must be finite and nonzero")
    return abs(derivative)


def _vertex_surface(
    cell,
    cell_temperature,
    distance,
    stefan,
    surface_t,
    config,
    radius,
    boundary,
    target_water,
    dt,
    carry,
):
    """The film-free surface node at the vertex (core2_water_film_node.evaluate_water_surface with the
    conservative absent elimination, y_s = y_b = 0): the reduced film's pure-water branch."""
    pressure = config.dry.pressure_pa
    pore = config.dry.pore
    face, left, right = pw.surface_face_flux(
        cell_temperature, surface_t, distance, stefan, pore, config
    )
    reference = ct.prepare_reduced_film_surface_state(
        boundary.film_boundary,
        cell_temperature_k=cell_temperature,
        pressure_pa=pressure,
        particle_radius_m=radius,
        outer_half_cell_distance_m=distance,
        particle_thermal_conductivity_w_m_k=config.dry.thermal_conductivity.value_w_m_k,
        binary_gas_interaction_k_wh=pore.k_wh,
    )
    film = ct._reduced_film_coefficients_at_surface_temperature(reference, surface_t)
    components = core.efr.component_fluxes_for_supplied_total_mass_flux(
        film,
        total_mass_flux_kg_m2_s=target_water,
        interface_hexane_mass_fraction=0.0,
        bulk_hexane_mass_fraction=0.0,
    )
    capacities = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
        surface_temperature_k=surface_t,
        bulk_temperature_k=boundary.temperature_k,
        pressure_pa=pressure,
        y_hexane=0.0,
        k_wh=pore.k_wh,
    )
    ack = core.efr.ackermann_heat_transfer(
        film,
        fluxes=components,
        water_film_heat_capacity_j_kg_k=capacities.water_j_kg_k,
        hexane_film_heat_capacity_j_kg_k=capacities.hexane_j_kg_k,
        enforce_production_domain=False,
    )
    if abs(ack.beta) > boundary.accepted_beta_limit:
        raise ValueError("vertex surface requires a beta outside the unchanged band")
    q_out = film.heat_transfer_w_m2_k * ack.factor * (surface_t - boundary.temperature_k)
    if isinstance(boundary, pd_boundary.KernelEnvelopeBoundary):
        q_out = pd_boundary.surface_gas_heat_w_m2(boundary, q_out, surface_t)
        if boundary.dew_deposition_kg_m2_s != 0.0:
            raise pd_boundary.DewLineCondensateOnFilmFreeSurface(
                "the dew-line projection books condensate onto the external film, but the "
                "vertex surface carries none"
            )
    external_energy = math.fsum(
        (
            q_out,
            components.water_kg_m2_s * right.water_gas_partial_enthalpy_j_mol / core.wa.M,
            components.hexane_kg_m2_s * right.hexane_gas_partial_enthalpy_j_mol / core.hx.M,
        )
    )
    pore_water = face.component.conserved_water_flux_mol_m2_s * core.wa.M
    pore_hexane = face.component.conserved_hexane_flux_mol_m2_s * core.hx.M
    liquid_h = core.wa.state_Tp(surface_t, pressure, "liquid").h_mass
    old_h = core.wa.state_Tp(carry["external_water_temperature_k"], pressure, "liquid").h_mass
    r_water, r_hexane, r_energy = node.node_balances(
        dt,
        0.0,
        carry["external_water_loading_kg_m2"],
        liquid_h,
        old_h,
        pore_water,
        pore_hexane,
        face.energy.total_energy_flux_w_m2,
        components.water_kg_m2_s,
        components.hexane_kg_m2_s,
        external_energy,
    )
    audit = node.WaterFilmNodeAudit(
        surface_t,
        0.0,
        0.0,
        components.water_kg_m2_s,
        components.hexane_kg_m2_s,
        external_energy,
        r_water,
        r_hexane,
        r_energy,
        1.0 - right.water_activity,
        ack.beta,
    )
    return face, audit


def vertex_advance(
    before,
    config,
    carry,
    surface,
    dt,
    boundary,
    *,
    assembly_budget,
    target_pressure_pa=None,
    search_tolerance=None,
    global_component_mass_error_budgets_kg=None,
):
    """One vertex dry advance over dt (water and energy rows, the absent film eliminated)."""
    root_target = dry.declared_dry_search_target(
        before.config.nonlinear_residual_tolerance, search_tolerance
    )
    grid = before.grid
    n = grid.n
    target_dry = before.config
    if target_pressure_pa is not None and target_pressure_pa != before.config.pressure_pa:
        target_dry = ct.pressure_only_target_config(before.config, target_pressure_pa)
        config = core.cut.pressure_only_target_config(config, target_pressure_pa)
    pw._require_vertex_mode(config)
    if carry["external_water_loading_kg_m2"] != 0.0:
        raise HexaneExhaustedVertexRefused("the vertex stage carries no exterior film")
    old = inventories(before)
    area = grid.areas[-1]
    distance = grid.R - grid.centers[-1]
    holder = {}

    def evaluate(values):
        if assembly_budget is not None:
            counter = assembly_budget
            if counter["used"] >= counter["limit"]:
                raise StepBudgetExceeded("dry physical-advance residual assembly budget exhausted")
            abandon_at = counter.get(envelope.HALF_STEP_ABANDON_KEY)
            if abandon_at is not None and counter["used"] >= abandon_at:
                raise envelope.DryFullStepAbandoned(
                    "kernel-envelope half step first: the film-bearing full step reached the "
                    "declared assembly threshold",
                    used=counter["used"],
                    threshold=abandon_at,
                )
            counter["used"] += 1
        ts = values[:n]
        stefan = values[n : 2 * n]
        surf_t = float(values[-1])
        cells = vertex_states(ts, target_dry)
        new = ct._inventories(grid, cells)
        fluxes = [extinction._zero_center_flux()]
        for k in range(1, n):
            fluxes.append(
                pw.interior_face_flux(
                    ts[k - 1],
                    ts[k],
                    grid.centers[k] - grid.centers[k - 1],
                    stefan[k - 1],
                    config.dry.pore,
                    config,
                )
            )
        target_water = (
            -math.fsum(b - a for a, b in zip(old.water_mol, new.water_mol))
            * core.wa.M
            / (area * dt)
        )
        face, audit = _vertex_surface(
            cells[-1],
            ts[-1],
            distance,
            stefan[-1],
            surf_t,
            config,
            grid.R,
            boundary,
            target_water,
            dt,
            carry,
        )
        fluxes.append(face)
        changes = extinction.ExtinctionInventoryChanges(
            tuple(b - a for a, b in zip(old.water_mol, new.water_mol)),
            tuple(b - a for a, b in zip(old.hexane_mol, new.hexane_mol)),
            tuple(b - a for a, b in zip(old.energy_j, new.energy_j)),
        )
        blocks = extinction._residual_blocks(grid, changes, fluxes, dt)
        if (
            any(value != 0.0 for value in blocks.dry_hexane_mol_s)
            or audit.hexane_node_residual_kg_m2_s != 0.0
        ):
            raise pw.PureWaterVertexError("a vertex n-hexane row is not exactly zero")
        holder.update(
            cells=cells, new=new, fluxes=tuple(fluxes), audit=audit, blocks=blocks, changes=changes
        )
        original = np.r_[
            blocks.dry_water_mol_s,
            blocks.dry_datum_covariant_energy_w,
            audit.energy_node_residual_w_m2,
        ]
        holder["original"] = original
        return original

    initial = np.r_[
        np.asarray(before.temperatures_k, dtype=float), np.zeros(n), surface.temperature_k
    ]
    tl, tu = config.dry.conditioned_temperature_domain.solver_bounds_k
    origins = np.r_[
        np.asarray(before.temperatures_k, dtype=float), np.zeros(n), surface.temperature_k
    ]
    bounds = ([tl] * n + [-math.inf] * n + [tl], [tu] * n + [math.inf] * n + [tu])
    raw = evaluate(initial)
    capacities = [
        _energy_capacity(t, before.config) * v for t, v in zip(before.temperatures_k, grid.volumes)
    ]
    scales = np.r_[
        np.maximum(np.asarray(old.water_mol) / dt, np.abs(raw[:n])),
        np.maximum(np.asarray(capacities) / dt, np.abs(raw[n : 2 * n])),
        max(abs(carry["external_water_energy_j"]) / (area * dt), 5000.0),
    ]
    component_weight = root_target / min(root_target, dry.NATIVE_COMPONENT_RESIDUAL_LIMIT)
    weights = np.ones(2 * n + 1)
    weights[:n] = component_weight

    def residual(values):
        return evaluate(values + origins) / scales * weights

    energy_rows = slice(n, 2 * n)
    acceptance = dry.NoiseAwareEnergyRowAcceptance(root_target, energy_rows, origins)
    solve_allowance = envelope.dry_assembly_allowance(np.asarray(initial).size)
    result = study.numerics.domain_aware_newton(
        residual,
        initial - origins,
        bounds=(np.asarray(bounds[0]) - origins, np.asarray(bounds[1]) - origins),
        max_nfev=solve_allowance,
        gtol=root_target,
        acceptance=acceptance,
    )
    search_final = residual(result.x)
    final = holder["original"] / scales
    search_gates, energy_floor = acceptance.gates(result.x, result.jac)
    final_gates = np.full(final.size, 1e-10)
    final_gates[energy_rows] = np.maximum(1e-10, search_gates[energy_rows])
    energy_gate_stamp = {
        "ruling": dry.GATE_UNITS_RULING,
        "declared_multiple": acceptance.multiple,
        "declared_lower_bound": root_target,
        "floor_definition": "sum_j |J_ij| numpy.spacing(x_j), J the domain-aware Newton's dense "
        "forward-difference Jacobian at the returned iterate, x the physical unknowns",
        "rows": [f"energy[{k}]" for k in range(n)],
        "residual": final[energy_rows].tolist(),
        "floor": energy_floor.tolist(),
        "gate": final_gates[energy_rows].tolist(),
        "absorbed_by_floor": (final_gates[energy_rows] > 1e-10).tolist(),
        "absorbed_rows": int(np.count_nonzero(final_gates[energy_rows] > 1e-10)),
        "maximum_residual_over_gate": float(
            np.max(np.abs(final[energy_rows]) / final_gates[energy_rows])
        ),
        "other_rows_gate": 1e-10,
        "search_target": root_target,
        "maximum_other_row_scaled_residual": float(
            np.max(np.abs(np.r_[final[:n], final[2 * n :]]))
        ),
        "native_ledger_energy_scale_factors": (search_gates[energy_rows] / root_target).tolist(),
    }
    if (
        not result.success
        or np.any(np.abs(search_final) > search_gates)
        or np.any(np.abs(final) > final_gates)
    ):
        refusal = HexaneExhaustedVertexRefused("pure-water vertex dry root refused")
        refusal.optimizer_debug = {
            "message": str(result.message),
            "maximum_scaled_residual": float(max(abs(final))),
            "scaled_residual": final.tolist(),
            "physical_coordinates": (result.x + origins).tolist(),
            "residual_calls": int(result.nfev),
            "noise_aware_energy_row_gate": energy_gate_stamp,
        }
        raise refusal
    cells, new, fluxes, audit = holder["cells"], holder["new"], holder["fluxes"], holder["audit"]
    ts = tuple(float(c.temperature_k) for c in cells)
    # The vertex's own step and cumulative ledger (coupled_transport._accept_step's arithmetic on
    # the water and energy components; the n-hexane components exactly zero).
    boundary_flux = fluxes[-1]
    water_out = dt * area * boundary_flux.component.conserved_water_flux_mol_m2_s
    hexane_out = dt * area * boundary_flux.component.conserved_hexane_flux_mol_m2_s
    energy_out = dt * area * boundary_flux.energy.total_energy_flux_w_m2
    water_change, hexane_change, energy_change = _differences(old, new)
    if hexane_out != 0.0 or hexane_change != 0.0:
        raise pw.PureWaterVertexError("a vertex step moved n-hexane")
    water_step = math.fsum((water_change, water_out))
    energy_step = math.fsum((energy_change, energy_out))
    normalized_water_step = ct._normalized_residual(
        water_step, old.total_water_mol, new.total_water_mol, water_out
    )
    normalized_energy_step = ct._normalized_residual(
        energy_step, old.total_energy_j, new.total_energy_j, energy_out
    )
    ledger_fields = _event_ledger(
        before, water_out, 0.0, energy_out, water_change, 0.0, energy_change
    )
    water_cumulative = math.fsum(
        (
            ledger_fields["cumulative_water_inventory_change_mol"],
            ledger_fields["cumulative_boundary_water_out_mol"],
        )
    )
    energy_cumulative = math.fsum(
        (
            ledger_fields["cumulative_energy_inventory_change_j"],
            ledger_fields["cumulative_boundary_energy_out_j"],
        )
    )
    normalized_water_cumulative = ct._normalized_residual(
        water_cumulative,
        before.reference_water_mol,
        new.total_water_mol,
        ledger_fields["cumulative_water_inventory_change_mol"],
        ledger_fields["cumulative_boundary_water_out_mol"],
        ledger_fields["cumulative_absolute_water_transfer_mol"],
    )
    normalized_energy_cumulative = ct._normalized_residual(
        energy_cumulative,
        before.reference_energy_j,
        new.total_energy_j,
        ledger_fields["cumulative_energy_inventory_change_j"],
        ledger_fields["cumulative_boundary_energy_out_j"],
        ledger_fields["cumulative_absolute_energy_transfer_j"],
    )
    blocks = holder["blocks"]
    max_scaled_water = max(abs(v) / s for v, s in zip(blocks.dry_water_mol_s, scales[:n]))
    energy_scales = scales[n : 2 * n] * (search_gates[energy_rows] / root_target)
    max_scaled_energy = max(abs(v) / s for v, s in zip(blocks.dry_energy_w, energy_scales))
    recovered = cp.recover_total_stefan_flux(
        grid,
        old.water_mol,
        new.water_mol,
        old.hexane_mol,
        new.hexane_mol,
        tuple(f.retained_water_flux_mol_m2_s for f in fluxes),
        dt,
    )
    minimum_entropy = min(f.total_entropy_production_w_m3_k for f in fluxes[1:])
    tolerance = target_dry.ledger_tolerance
    vertex_ledger = {
        "water_step_residual_mol": water_step,
        "energy_step_residual_j": energy_step,
        "normalized_water_step_residual": normalized_water_step,
        "normalized_energy_step_residual": normalized_energy_step,
        "normalized_water_cumulative_residual": normalized_water_cumulative,
        "normalized_energy_cumulative_residual": normalized_energy_cumulative,
        "max_scaled_water_residual": max_scaled_water,
        "max_scaled_energy_residual": max_scaled_energy,
        "stefan_normalized_telescope_residual": recovered.normalized_telescope_residual,
        "minimum_total_entropy_production_w_m3_k": minimum_entropy,
        "boundary_water_out_mol": water_out,
        "boundary_hexane_out_mol": hexane_out,
        "boundary_energy_out_j": energy_out,
        "ledger_tolerance": tolerance,
        "component_limit": min(
            dry.NATIVE_COMPONENT_RESIDUAL_LIMIT, target_dry.nonlinear_residual_tolerance
        ),
    }
    if (
        max_scaled_water > vertex_ledger["component_limit"]
        or max_scaled_energy > target_dry.nonlinear_residual_tolerance
        or max(normalized_water_step, normalized_energy_step) > tolerance
        or max(normalized_water_cumulative, normalized_energy_cumulative) > tolerance
        or recovered.normalized_telescope_residual > tolerance
        or minimum_entropy < 0.0
    ):
        refusal = HexaneExhaustedVertexRefused(
            "pure-water vertex step failed its acceptance ledger"
        )
        refusal.optimizer_debug = {"vertex_ledger": vertex_ledger}
        raise refusal
    combined_raw = [
        (new.total_water_mol - old.total_water_mol) * core.wa.M
        + dt * area * audit.external_water_flux_kg_m2_s,
        (new.total_hexane_mol - old.total_hexane_mol) * core.hx.M
        + dt * area * audit.external_hexane_flux_kg_m2_s,
        new.total_energy_j
        - old.total_energy_j
        - carry["external_water_energy_j"]
        + dt * area * audit.external_energy_flux_w_m2,
    ]
    combined_scales = [
        max(old.total_water_mol * core.wa.M + area * carry["external_water_loading_kg_m2"], 1e-300),
        max(old.total_hexane_mol * core.hx.M, 1e-300),
        max(
            abs(old.total_energy_j) + abs(carry["external_water_energy_j"]), sum(capacities), 1e-300
        ),
    ]
    normalized = [abs(r) / s for r, s in zip(combined_raw, combined_scales)]
    if max(normalized) > 1e-10:
        refusal = HexaneExhaustedVertexRefused("vertex particle-plus-water-film dry ledger refused")
        refusal.optimizer_debug = {"combined_normalized_residuals": normalized}
        raise refusal
    if global_component_mass_error_budgets_kg is not None:
        defect = math.fsum((*holder["original"][:n],)) * core.wa.M * dt
        if (
            abs(combined_raw[0]) > global_component_mass_error_budgets_kg[0]
            or abs(defect) > global_component_mass_error_budgets_kg[0]
        ):
            raise HexaneExhaustedVertexRefused(
                "vertex water closure exceeds the declared per-advance budget"
            )
    after = VertexDryTransportState(
        grid=grid,
        config=target_dry,
        time_s=before.time_s + dt,
        temperatures_k=ts,
        y_hexane=tuple(0.0 for _ in ts),
        **ledger_fields,
    )
    report = {
        "time_s": after.time_s,
        "dt_s": dt,
        "accepted": True,
        "maximum_scaled_residual": float(max(abs(final))),
        "maximum_search_scaled_residual": float(max(abs(search_final))),
        "dry_search_target": root_target,
        "noise_aware_energy_row_gate": energy_gate_stamp,
        "line_search": getattr(result, "line_search", None),
        "newton_residual_calls": int(result.nfev),
        "newton_certificate_calls_outside_search_budget": int(
            getattr(result, "certificate_jacobian_calls_outside_search_budget", 0)
        ),
        "newton_message": str(result.message),
        "surface_search_basis": "algebraic pure-water film elimination (vertex, absent film)",
        "ledger": vertex_ledger,
        "surface": asdict(audit),
        "combined_normalized_residuals": normalized,
        "external_water_phase": "absent",
        "temperature_k": list(ts),
        "y_hexane": [0.0] * n,
        "water_mol": new.total_water_mol,
        "hexane_mol": new.total_hexane_mol,
        "energy_j": new.total_energy_j,
        "unknowns": int(np.asarray(initial).size),
        "solve_allowance": solve_allowance,
    }
    next_carry = {
        **carry,
        "external_water_loading_kg_m2": audit.liquid_water_kg_m2,
        "external_water_energy_j": 0.0,
        "external_water_temperature_k": audit.temperature_k,
        "external_water_mass_kg": area * audit.liquid_water_kg_m2,
    }
    return after, audit, next_carry, report


def regime_stamp(kind, boundary, extra=None):
    return {
        "regime": REGIME,
        "kind": kind,
        "presented_gas": _gas_stamp(boundary),
        "gas_taken_as_hexane_free": kind == "vertex_advance",
        "exhaustion_loading_kg_kg": envelope.HEXANE_EXHAUSTED_LOADING_KG_KG,
        "record": RECORD,
        "physically_qualifying": False,
        **(extra or {}),
    }


def route(
    before,
    config,
    carry,
    surface,
    duration,
    boundary_factory,
    *,
    advance_type,
    binary,
    assembly_budget,
    target_pressure_pa,
    search_tolerance,
    global_component_mass_error_budgets_kg,
):
    """The armed dispatch at the start of a dry advance: None (the binary stage runs unchanged) or
    the DryPhaseAdvance of a vertex advance (with the release event when entering) or of the binary
    stage after a re-entry event."""
    start = boundary_factory(before.time_s)
    end = boundary_factory(before.time_s + duration)
    clamp = envelope.HEXANE_EXHAUSTED_GAS_Y_HEXANE_CLAMP
    events = []
    ignored = None
    if is_vertex(before) and envelope.armed(envelope.HEXANE_EXHAUSTED_REENTRY_HYSTERESIS):
        # the re-entry hysteresis (armed only): re-entry above the declared fraction; between the
        # clamp and it the vertex continues, the gas's y ignored and stamped
        reentry = envelope.HEXANE_EXHAUSTED_REENTRY_GAS_Y_HEXANE
        ys = (presented_gas_y(start), presented_gas_y(end))
        if max(ys) > clamp and max(ys) <= reentry:
            ignored = {"ignored_gas_y_hexane": list(ys), "reentry_gas_y_hexane": reentry}
            clamp = reentry
    if is_vertex(before):
        if presented_gas_y(start) > clamp or presented_gas_y(end) > clamp:
            converted, event = reentry_event(before, end)
            seed = replace(surface, surface_y_hexane=presented_gas_y(end))
            outcome = binary(converted, seed, PHASE)
            stamp = regime_stamp(REENTRY, end, {"event": event})
            return replace(
                outcome,
                native_report={
                    **outcome.native_report,
                    "kernel_envelope_hexane_exhausted_regime": stamp,
                },
                regime_events=(event,),
            )
        state = before
    else:
        loading = trigger(before, carry, start, end)
        if loading is None:
            return None
        state, event = release_event(before, surface, end, loading)
        events.append(event)
    try:
        after, audit, next_carry, report = vertex_advance(
            state,
            config,
            carry,
            surface,
            duration,
            end,
            assembly_budget=assembly_budget,
            target_pressure_pa=target_pressure_pa,
            search_tolerance=search_tolerance,
            global_component_mass_error_budgets_kg=global_component_mass_error_budgets_kg,
        )
    except StepBudgetExceeded:
        raise
    except Exception as error:
        from particle_postcore_solver import DryPhaseAdvanceError

        raise DryPhaseAdvanceError(
            "pure-water vertex advance refused",
            before,
            [
                {
                    "phase": PHASE,
                    "error": f"{type(error).__name__}: {error}"[:400],
                    "optimizer_debug": getattr(error, "optimizer_debug", None),
                    "release_event": bool(events),
                }
            ],
        ) from error
    report["kernel_envelope_hexane_exhausted_regime"] = regime_stamp(
        "vertex_advance",
        end,
        {
            "events": [e["kind"] for e in events],
            **({} if ignored is None else {"reentry_hysteresis": ignored}),
        },
    )
    return advance_type(
        after, audit, next_carry, report, PHASE, False, None, (), regime_events=tuple(events)
    )


def declared_row():
    return envelope.hexane_exhausted_regime_declared_row()
