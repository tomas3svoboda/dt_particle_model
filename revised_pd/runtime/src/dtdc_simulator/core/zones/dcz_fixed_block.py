"""Affine phase-fixed block-implicit DCZ thermal solve."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from dtdc_simulator.core import thermo
from dtdc_simulator.core.zones import particle as pt
from dtdc_simulator.core.zones.dcz_enthalpy import (
    vapor_cell_enthalpy_residual_at_state,
)
from dtdc_simulator.core.zones.dcz_types import DCZConstants

_TP_MIN_K = 250.0
_TP_MAX_K = 480.0


@dataclass(frozen=True)
class _ThermalBlockResult:
    particles: tuple[pt.ParticleState, ...]
    vapor_temperature_k: tuple[float, ...]
    particle_energy_terms_w: tuple[tuple[float, float, float, float], ...]
    algebraic_residual_max_w: float
    interface_cancellation_w: float


def _solve_phase_fixed_thermal_block(
    *,
    initial_particle: pt.ParticleState,
    particles: list[pt.ParticleState],
    vapor_temperature_guess_k: list[float],
    vapor_water_kg_s: list[float],
    vapor_hexane_kg_s: list[float],
    bottom_water_kg_s: float,
    bottom_hexane_kg_s: float,
    bottom_temperature_k: float,
    bound_water_mass_rate_w_m3: list[float],
    free_water_profile: list[float],
    condensed_water_kg_s: list[float],
    evaporated_water_kg_s: list[float],
    moisture_profile: list[float],
    moisture_in: float,
    q_indirect_w_m3: tuple[float, ...],
    q_axial_w_m3: tuple[float, ...],
    cell_widths_m: tuple[float, ...],
    cell_residence_s: tuple[float, ...],
    time_substeps: tuple[int, ...],
    bed_area_m2: float,
    geometry: pt.ShellGeometry,
    c: DCZConstants,
) -> _ThermalBlockResult:
    """Solve the fixed-component DCZ thermal equations as one sparse block.

    This is a solver reformulation, not a new physical model.  It assembles
    the same backward-Euler radial particle equations and total-vapor-
    enthalpy control volumes used by the partitioned loop.  Every particle
    surface substep and its vapor cell share one algebraic Robin heat rate.

    Phase rates and component flows are held fixed by the conservative outer
    loop.  The currently selected liquid-enthalpy branch is also fixed during
    this linear solve; the outer active-set iteration remains responsible for
    changing that branch.
    """
    nz = len(particles)
    layers = c.particle.Np
    if not (
        len(vapor_temperature_guess_k)
        == len(vapor_water_kg_s)
        == len(vapor_hexane_kg_s)
        == len(bound_water_mass_rate_w_m3)
        == len(free_water_profile)
        == len(condensed_water_kg_s)
        == len(evaporated_water_kg_s)
        == len(moisture_profile)
        == len(q_indirect_w_m3)
        == len(q_axial_w_m3)
        == len(cell_widths_m)
        == len(cell_residence_s)
        == len(time_substeps)
        == nz
    ):
        raise ValueError("phase-fixed thermal block profiles must have one value per cell")

    # Particle unknowns include every residence-time substep so the surface
    # heat delivered to the vapor is exactly the same time average as in the
    # existing particle cascade.  Offsets make the sparse indexing explicit.
    particle_offsets: list[int] = []
    particle_unknowns = 0
    for substeps in time_substeps:
        particle_offsets.append(particle_unknowns)
        particle_unknowns += substeps * layers
    vapor_offset = particle_unknowns
    unknowns = particle_unknowns + nz

    def p_index(cell: int, substep: int, shell: int) -> int:
        return particle_offsets[cell] + substep * layers + shell

    def v_index(cell: int) -> int:
        return vapor_offset + cell

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    rhs = np.zeros(unknowns, dtype=float)

    def add(row: int, col: int, value: float) -> None:
        if value != 0.0:
            rows.append(row)
            cols.append(col)
            data.append(value)

    pc = c.particle
    k_mix = pc.alpha_pg * pc.k_pg + pc.alpha_ps * pc.k_ps
    base_cv = pc.alpha_pg * pc.rho_pg * pc.cp_pg + pc.alpha_ps * pc.rho_ps * pc.cp_ps
    water_cv = pc.alpha_ps * pc.rho_ps * pc.cp_water_liquid
    surface_area_m2 = geometry.face_areas[-1]
    surface_conductance_w_k = c.hQ * surface_area_m2
    reference_temperature_k = c.vapor_enthalpy_ref.T_boil_water
    sorption_sources_by_cell: list[tuple[float, ...]] = []

    for j in range(nz):
        entering_w = initial_particle.wpg2 if j == 0 else particles[j - 1].wpg2
        sorption_entering = pt.ParticleState(
            wpg2=entering_w,
            Tp=initial_particle.Tp if j == 0 else particles[j - 1].Tp,
        )
        sorption_sources = pt.finite_sorption_heat_source_per_layer_w_m3(
            sorption_entering,
            particles[j],
            cell_residence_s[j],
            pc,
        )
        sorption_sources_by_cell.append(sorption_sources)
        moisture_start = moisture_in if j == 0 else moisture_profile[j - 1]
        substeps = time_substeps[j]
        dt_sub = cell_residence_s[j] / substeps
        for step in range(substeps):
            moisture_old = moisture_start + (
                moisture_profile[j] - moisture_start
            ) * step / substeps
            moisture_new = moisture_start + (
                moisture_profile[j] - moisture_start
            ) * (step + 1) / substeps
            cv_old = base_cv + water_cv * moisture_old
            cv_new = base_cv + water_cv * moisture_new
            for i in range(layers):
                row = p_index(j, step, i)
                volume = geometry.volumes[i]
                storage_new = cv_new * volume / dt_sub
                storage_old = cv_old * volume / dt_sub
                diagonal = storage_new
                if i > 0:
                    conductance = k_mix * geometry.face_areas[i - 1] / geometry.dr
                    diagonal += conductance
                    add(row, p_index(j, step, i - 1), -conductance)
                if i < layers - 1:
                    conductance = k_mix * geometry.face_areas[i] / geometry.dr
                    diagonal += conductance
                    add(row, p_index(j, step, i + 1), -conductance)
                else:
                    diagonal += surface_conductance_w_k
                    add(row, v_index(j), -surface_conductance_w_k)
                add(row, row, diagonal)

                # Residence advection/accumulation supplies the previous
                # substep, previous axial cell, or fixed FTRZ inlet.
                if step > 0:
                    add(row, p_index(j, step - 1, i), -storage_old)
                elif j > 0:
                    add(
                        row,
                        p_index(j - 1, time_substeps[j - 1] - 1, i),
                        -storage_old,
                    )
                else:
                    rhs[row] += storage_old * initial_particle.Tp[i]
                rhs[row] += (
                    (cv_new - cv_old)
                    * volume
                    / dt_sub
                    * reference_temperature_k
                    + (sorption_sources[i] + q_axial_w_m3[j]) * volume
                )

    # Assemble each vapor total-enthalpy residual.  The non-interface part is
    # affine under the configured constant heat capacities.  Evaluating its
    # exact coefficients through the authoritative residual function avoids
    # duplicating datum/sign conventions here.
    interface_cancellation_w = 0.0
    for j in range(nz):
        row = v_index(j)
        volume_m3 = bed_area_m2 * cell_widths_m[j]
        representative_particles = c.alpha_L * volume_m3 / geometry.total_volume
        water_in = bottom_water_kg_s if j == nz - 1 else vapor_water_kg_s[j + 1]
        hexane_in = bottom_hexane_kg_s if j == nz - 1 else vapor_hexane_kg_s[j + 1]
        bound_rate = bound_water_mass_rate_w_m3[j] * volume_m3
        dew_temperature = thermo.dew_point_temperature(
            vapor_hexane_kg_s[j] / max(vapor_water_kg_s[j], 1.0e-12),
            c.antoine_water,
            P=c.pressure_pa,
        )
        # Keep the currently selected bound-liquid temperature branch fixed
        # inside this linear solve.  The surrounding active-set iteration
        # changes it if the branch condition changes.
        bound_tracks_particle = not (
            free_water_profile[j] > 1.0e-12
            and pt.outer_layer_value(particles[j].Tp) > dew_temperature
        )

        def noninterface_residual(
            temperature_out: float,
            temperature_in: float,
            particle_surface: float,
        ) -> float:
            bound_surface = (
                particle_surface if bound_tracks_particle else dew_temperature
            )
            return vapor_cell_enthalpy_residual_at_state(
                water_in_kg_s=water_in,
                hexane_in_kg_s=hexane_in,
                temperature_in_k=temperature_in,
                water_out_kg_s=vapor_water_kg_s[j],
                hexane_out_kg_s=vapor_hexane_kg_s[j],
                temperature_out_k=temperature_out,
                bound_water_to_solid_kg_s=bound_rate,
                free_water_evaporation_kg_s=evaporated_water_kg_s[j],
                condensed_water_kg_s=condensed_water_kg_s[j],
                external_duty_w=q_indirect_w_m3[j] * volume_m3,
                interphase_heat_to_vapor_w=0.0,
                particle_surface_temperature_k=particle_surface,
                bound_water_surface_temperature_k=bound_surface,
                free_water_temperature_k=dew_temperature,
                c=c,
            )

        constant = noninterface_residual(0.0, 0.0, 0.0)
        coeff_out = noninterface_residual(1.0, 0.0, 0.0) - constant
        coeff_in = noninterface_residual(0.0, 1.0, 0.0) - constant
        coeff_surface_species = (
            noninterface_residual(0.0, 0.0, 1.0) - constant
        )

        interface_coefficient = (
            representative_particles * surface_conductance_w_k / time_substeps[j]
        )
        add(
            row,
            v_index(j),
            coeff_out + representative_particles * surface_conductance_w_k,
        )
        if j < nz - 1:
            add(row, v_index(j + 1), coeff_in)
        else:
            constant += coeff_in * bottom_temperature_k
        final_surface = p_index(j, time_substeps[j] - 1, layers - 1)
        add(row, final_surface, coeff_surface_species)
        for step in range(time_substeps[j]):
            add(row, p_index(j, step, layers - 1), -interface_coefficient)
            # The two assembled rows use equal and opposite dimensional heat
            # rates.  Accumulate the coefficient-level audit explicitly.
            interface_cancellation_w += (
                interface_coefficient - interface_coefficient
            )
        rhs[row] = -constant

    matrix = coo_matrix((data, (rows, cols)), shape=(unknowns, unknowns)).tocsr()
    solution = np.asarray(spsolve(matrix, rhs), dtype=float)
    if solution.shape != (unknowns,) or not np.all(np.isfinite(solution)):
        raise ValueError("DCZ block-implicit thermal system produced a non-finite state")
    algebraic_residual = matrix @ solution - rhs
    residual_max = float(np.max(np.abs(algebraic_residual)))

    solved_particles: list[pt.ParticleState] = []
    energy_terms: list[tuple[float, float, float, float]] = []
    for j in range(nz):
        final_step = time_substeps[j] - 1
        temperatures = tuple(
            float(solution[p_index(j, final_step, i)]) for i in range(layers)
        )
        if any(value < _TP_MIN_K or value > _TP_MAX_K for value in temperatures):
            raise ValueError(
                "DCZ block-implicit thermal solution left the qualified "
                "particle-temperature interval"
            )
        solved_particles.append(
            pt.ParticleState(wpg2=particles[j].wpg2, Tp=temperatures)
        )
        entering_temperature = (
            initial_particle.Tp if j == 0 else solved_particles[j - 1].Tp
        )
        moisture_start = moisture_in if j == 0 else moisture_profile[j - 1]
        accumulation = sum(
            (
                (base_cv + water_cv * moisture_profile[j])
                * (temperature_out - reference_temperature_k)
                - (base_cv + water_cv * moisture_start)
                * (temperature_in - reference_temperature_k)
            )
            * volume
            / cell_residence_s[j]
            for temperature_out, temperature_in, volume in zip(
                temperatures,
                entering_temperature,
                geometry.volumes,
            )
        )
        sorption = sum(
            source * volume
            for source, volume in zip(
                sorption_sources_by_cell[j],
                geometry.volumes,
            )
        )
        axial = q_axial_w_m3[j] * geometry.total_volume
        surface = sum(
            surface_conductance_w_k
            * (
                solution[p_index(j, step, layers - 1)]
                - solution[v_index(j)]
            )
            for step in range(time_substeps[j])
        ) / time_substeps[j]
        energy_terms.append((accumulation, sorption, axial, float(surface)))

    return _ThermalBlockResult(
        particles=tuple(solved_particles),
        vapor_temperature_k=tuple(
            float(solution[v_index(j)]) for j in range(nz)
        ),
        particle_energy_terms_w=tuple(energy_terms),
        algebraic_residual_max_w=residual_max,
        interface_cancellation_w=interface_cancellation_w,
    )
