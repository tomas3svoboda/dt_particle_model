"""Existing conservative core endpoint with explicit conditional time control.

The measured positive-target sequence supplies a numerical truncation estimate,
not a rigorous time enclosure. This module neither changes a native gate nor
substitutes a zero-width state for an unresolved positive-core root.
"""

from contextlib import contextmanager
from dataclasses import dataclass, asdict, replace
import math

import numpy as np

import experiment_core2_dt_projection as projection
from face_volume_coalescence import restore_fully_dry_face


@contextmanager
def conservative_water_surface_at_dry_endpoint(before, carry, surface, dt, boundary):
    """Positive-film endpoint surface with algebraic component elimination.

    Declared option of the core-extinction projection (default off). For each
    projection trial the native exterior pore fluxes at a trial surface
    temperature fix the film hexane flux (unchanged native Stefan film inverted
    at the water-saturated surface composition over its beta band), and the
    node water balance fixes the liquid inventory explicitly. The surface
    temperature then solves the node energy row. The node water and hexane
    rows hold by construction; all three original node rows are still checked
    at the former 1e-11 surface tolerance, and the projection's native and
    combined ledgers are checked unchanged afterwards.
    """
    from particle_conservative_surface import total_flux_for_hexane_target

    core, driver, study = projection.core, projection.driver, projection.study
    extinction = projection.extinction
    config = before.transport.config
    grid = before.transport.geometry.master_grid
    old_liquid = carry["external_water_loading_kg_m2"]
    old_temperature = carry["external_water_temperature_k"]
    if not old_liquid > 0:
        raise ValueError("conservative endpoint surface requires a positive incoming film")
    pressure = config.dry.pressure_pa
    holder = {"x": np.array([surface.temperature_k, 0.0, old_liquid])}
    mass_scale = max(old_liquid / dt, 0.05)
    heat_scale = max(abs(carry["external_water_energy_j"]) / (grid.areas[-1] * dt), 5000.0)
    tl, tu = config.dry.conditioned_temperature_domain.solver_bounds_k
    original = extinction._fully_dry_face_fluxes
    work = {
        "surface_solves": 0,
        "surface_residual_calls": 0,
        "maximum_surface_residual": 0.0,
        "surface_basis": "algebraic native-film component elimination",
    }
    guess = np.array([surface.temperature_k])
    with driver.liquid_internal_energy_node(holder, dt, old_liquid, old_temperature, pressure):

        def faces(old, candidate, cells, actual_boundary, **kwargs):
            nonlocal guess
            if not isinstance(actual_boundary, study.FiniteFilmStudyBoundary):
                return original(old, candidate, cells, actual_boundary, **kwargs)
            cell = cells[-1]
            pore = replace(config.dry.pore, w_o=before.transport.oil_fraction_labels[-1])
            stefan = candidate.positive_area_total_stefan_fluxes_mol_m2_s[-1]
            distance = grid.R - grid.centers[-1]
            face = None

            def surface_at(ts):
                y = core.cp.decode_gas_only_y(ts, pressure, -32.0, pore)
                holder["x"] = np.array([ts, 0.0, old_liquid])
                probe = driver.node.evaluate_water_surface(
                    holder,
                    dt,
                    old_liquid,
                    old_temperature,
                    cell,
                    distance,
                    stefan,
                    pore,
                    config,
                    grid.R,
                    actual_boundary,
                    surface_y_hexane=y,
                )
                pore_hexane = probe.component.conserved_hexane_flux_mol_m2_s * core.hx.M
                pore_water = probe.component.conserved_water_flux_mol_m2_s * core.wa.M
                reference = core.ct.prepare_reduced_film_surface_state(
                    actual_boundary.film_boundary,
                    cell_temperature_k=cell.temperature_k,
                    pressure_pa=pressure,
                    particle_radius_m=grid.R,
                    outer_half_cell_distance_m=distance,
                    particle_thermal_conductivity_w_m_k=config.dry.thermal_conductivity.value_w_m_k,
                    binary_gas_interaction_k_wh=pore.k_wh,
                )
                film = core.ct._reduced_film_coefficients_at_surface_temperature(reference, ts)
                cp = core.efr.phy053_partial_enthalpy_secant_heat_capacities(
                    surface_temperature_k=ts,
                    bulk_temperature_k=actual_boundary.temperature_k,
                    pressure_pa=pressure,
                    y_hexane=y,
                    k_wh=pore.k_wh,
                )
                total, components, _, _ = total_flux_for_hexane_target(
                    film,
                    pore_hexane,
                    core._hexane_mass_fraction(y),
                    core._hexane_mass_fraction(actual_boundary.y_hexane),
                    cp.water_j_kg_k,
                    cp.hexane_j_kg_k,
                    actual_boundary.accepted_beta_limit,
                )
                liquid = math.fsum((old_liquid, dt * pore_water, -dt * components.water_kg_m2_s))
                if not 0 < liquid < max(0.05, 4 * old_liquid):
                    raise ValueError("endpoint liquid inventory left its original open bounds")
                holder["x"] = np.array([ts, total, liquid])
                return driver.node.evaluate_water_surface(
                    holder,
                    dt,
                    old_liquid,
                    old_temperature,
                    cell,
                    distance,
                    stefan,
                    pore,
                    config,
                    grid.R,
                    actual_boundary,
                    surface_y_hexane=y,
                )

            def residual(values):
                nonlocal face
                face = surface_at(float(values[0]))
                work["surface_residual_calls"] += 1
                return np.array([face.surface_film_audit.energy_node_residual_w_m2 / heat_scale])

            solution = study.numerics.centered_temperature_solve(
                study.numerics.domain_aware_newton,
                residual,
                guess,
                temperature_index=0,
                bounds=([tl], [tu]),
                max_nfev=160,
                gtol=1e-11,
            )
            residual(solution.x)
            s = face.surface_film_audit
            rows = (
                s.water_node_residual_kg_m2_s / mass_scale,
                s.hexane_node_residual_kg_m2_s / mass_scale,
                s.energy_node_residual_w_m2 / heat_scale,
            )
            maximum = float(max(abs(r) for r in rows))
            if not solution.success or maximum > 1e-11:
                raise RuntimeError("native endpoint water node refused")
            guess = np.array(solution.x, copy=True)
            work["surface_solves"] += 1
            work["maximum_surface_residual"] = max(work["maximum_surface_residual"], maximum)
            resolved = core.ct.DirichletPoreBoundary(
                cell.temperature_k, cell.y_hexane, "native internal faces only"
            )
            internal = original(old, candidate, cells, resolved, **kwargs)
            return (*internal[:-1], face)

        extinction._fully_dry_face_fluxes = faces
        try:
            yield work
        finally:
            extinction._fully_dry_face_fluxes = original


class CoreTimeNotResolved(ValueError):
    """A valid linked trajectory still needs more accepted positive targets."""


@dataclass(frozen=True)
class PositiveCoreTargetEvidence:
    before: object
    after: object
    dt_s: float

    def __post_init__(self):
        ci = projection.core.ci
        if not isinstance(self.before, ci.CutIntegratorState) or not isinstance(
            self.after, ci.CutIntegratorState
        ):
            raise TypeError("positive-core evidence requires native cut-integrator states")
        old, new = self.before.transport, self.after.transport
        if not math.isfinite(self.dt_s) or self.dt_s <= 0 or new.time_s != old.time_s + self.dt_s:
            raise ValueError("target evidence must preserve the exact native physical clock")
        if new.geometry.front.z != 0.5 * old.geometry.front.z:
            raise ValueError("terminal sequence requires its declared positive half-volume targets")
        if (
            old.config is not new.config
            or self.before.controls is not self.after.controls
            or old.geometry.master_grid is not new.geometry.master_grid
            or self.before.reference_inventory is not self.after.reference_inventory
        ):
            raise ValueError(
                "target evidence changed native config/grid/control/reference identities"
            )


@dataclass(frozen=True)
class ConditionalCoreTime:
    last_four_durations_s: tuple[float, ...]
    observed_ratios: tuple[float, ...]
    assumed_future_ratio_upper: float
    estimated_remaining_time_s: float
    conditional_remaining_time_bound_s: float
    requested_time_tolerance_s: float
    last_wet_volume_fraction: float
    rigorous_time_enclosure_claimed: bool = False


def estimate_core_time(before, targets, maximum_remaining_time_s=0.002):
    """Reuse the predeclared four-target, 0.82-contraction, 2-ms ceiling.

    Identity and clock checks prevent stitching unrelated states into a
    duration sequence. A tighter time ceiling may be used by a refinement
    ladder. A future-ratio assumption remains explicit in every result.
    """
    if not math.isfinite(maximum_remaining_time_s) or not 0 < maximum_remaining_time_s <= 0.002:
        raise ValueError("terminal time control may preserve or tighten the existing 2-ms ceiling")
    if len(targets) < 4:
        raise CoreTimeNotResolved("four consecutive accepted positive targets required")
    recent = tuple(targets[-4:])
    if not all(isinstance(t, PositiveCoreTargetEvidence) for t in recent):
        raise TypeError("native positive-target evidence records required")
    if recent[-1].after is not before or any(
        a.after is not b.before for a, b in zip(recent, recent[1:])
    ):
        raise ValueError("positive-target sequence does not preserve exact linked state identities")
    z = before.transport.geometry.front.z
    if not 0 < z <= 1e-12:
        raise CoreTimeNotResolved("positive core exceeds the existing terminal volume ceiling")
    durations = tuple(t.dt_s for t in recent)
    ratios = tuple(b / a for a, b in zip(durations, durations[1:]))
    if not all(0 < q < 0.82 for q in ratios):
        raise CoreTimeNotResolved("measured durations violate the predeclared 0.82 contraction")
    bound = durations[-1] * 0.82 / (1.0 - 0.82)
    if bound > maximum_remaining_time_s:
        raise CoreTimeNotResolved(
            "conditional remaining time exceeds the declared refinement tolerance"
        )
    estimate = durations[-1] * max(ratios) / (1.0 - max(ratios))
    return ConditionalCoreTime(
        durations, ratios, 0.82, estimate, bound, maximum_remaining_time_s, z
    )


@dataclass(frozen=True)
class CoreEndpointAdvance:
    after: object
    surface: object
    carry: dict
    native_step: object
    time_control: ConditionalCoreTime
    base_face_restoration: dict | None


def project_core_endpoint(
    before,
    carry,
    surface,
    targets,
    missing_face,
    boundary_factory,
    *,
    maximum_remaining_time_s=0.002,
    conservative_surface=False,
):
    """Apply the existing native fixed-time projection under caller work scope.

    conservative_surface (declared option, default off) replaces the endpoint
    surface sub-solve by conservative_water_surface_at_dry_endpoint when an
    exterior film is present; gas-only endpoints keep the former sub-solve.
    """
    estimate = estimate_core_time(before, targets, maximum_remaining_time_s)
    dt = estimate.estimated_remaining_time_s
    boundary = boundary_factory(before.transport.time_s + dt)
    candidate = projection.extinction.ExtinctionUnknowns(
        before.transport.dry_temperatures_k,
        before.transport.dry_y_hexane,
        before.last_total_stefan_fluxes_mol_m2_s[1:],
        dt,
    )
    seed = projection.projection.ExtinctionProjectionSeed(
        candidate,
        tuple(max(abs(v), 0.01) for v in candidate.positive_area_total_stefan_fluxes_mol_m2_s),
        "consecutive native positive targets; disclosed conditional event-time truncation",
    )
    surface_context = (
        conservative_water_surface_at_dry_endpoint
        if conservative_surface and carry["external_water_loading_kg_m2"] > 0
        else projection.water_surface_at_dry_endpoint
    )
    with surface_context(before, carry, surface, dt, boundary):
        step = projection.projection.solve_fixed_time_extinction_projection(
            before, dt, boundary, seed
        )
    new_surface = step.assembly.dry_face_fluxes[-1].surface_film_audit
    combined = projection.driver.combined_ledger(
        before, before.transport.config, carry, step, new_surface
    )
    if combined["maximum_normalized_residual"] > before.controls.ledger_tolerance:
        raise RuntimeError("native combined conservation refused the core endpoint")
    after = step.after.to_uniform_fully_dry_transport_state()
    restoration = None
    if missing_face is not None:
        after, restoration = restore_fully_dry_face(after, missing_face)
    next_carry = {
        **carry,
        **combined,
        "external_water_loading_kg_m2": new_surface.liquid_water_kg_m2,
        "core_endpoint_time_control": asdict(estimate),
    }
    return CoreEndpointAdvance(after, new_surface, next_carry, step, estimate, restoration)
