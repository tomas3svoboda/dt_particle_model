"""Optional Numba backend for complete DCZ particle cascades."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from dtdc_simulator.core.zones.particle import ParticleConstants, ParticleState, ShellGeometry

try:
    from numba import njit
except ImportError:  # pragma: no cover - exercised in environments without the optional extra
    njit = None

NUMBA_AVAILABLE = njit is not None
JIT_DISABLED = os.environ.get("DTDC_DISABLE_JIT", "").lower() in {"1", "true", "yes"}


@dataclass(frozen=True)
class CascadeInvariants:
    initial_w: np.ndarray
    initial_T: np.ndarray
    volumes: np.ndarray
    face_areas: np.ndarray


def build_invariants(
    initial_particle: "ParticleState",
    geometry: "ShellGeometry",
) -> CascadeInvariants:
    """Materialize immutable cascade arrays once for a complete DCZ solve."""
    arrays = (
        np.asarray(initial_particle.wpg2, dtype=np.float64),
        np.asarray(initial_particle.Tp, dtype=np.float64),
        np.asarray(geometry.volumes, dtype=np.float64),
        np.asarray(geometry.face_areas, dtype=np.float64),
    )
    for array in arrays:
        array.setflags(write=False)
    return CascadeInvariants(*arrays)


if njit is not None:

    @njit(cache=True)
    def _gab_loading(
        activity: float,
        temperature: float,
        gab_Xm: float,
        gab_C0: float,
        gab_dHC_R: float,
        gab_K0: float,
        gab_dHK_R: float,
    ) -> float:
        C = gab_C0 * np.exp(gab_dHC_R / temperature)
        K = gab_K0 * np.exp(gab_dHK_R / temperature)
        clamped_activity = activity
        if K > 0.0 and K * clamped_activity > 0.999:
            clamped_activity = 0.999 / K
        return (
            gab_Xm
            * C
            * K
            * clamped_activity
            / (
                (1.0 - K * clamped_activity)
                * (1.0 - K * clamped_activity + C * K * clamped_activity)
            )
        )

    @njit(cache=True)
    def _sorption_enthalpy_potential(
        loading: float,
        gab_Xm: float,
        dH_vap_hexane: float,
        sorption_C0: float,
        sorption_C1: float,
    ) -> float:
        W = max(loading, 0.0)
        floor = 0.02 * gab_Xm
        if abs(sorption_C1 + 1.0) < 1.0e-12:
            if W <= floor:
                excess_integral = W / floor
            else:
                excess_integral = 1.0 + np.log(W / floor)
        elif W <= floor:
            excess_integral = floor**sorption_C1 * W
        else:
            floor_term = floor ** (sorption_C1 + 1.0)
            excess_integral = floor_term + (
                W ** (sorption_C1 + 1.0) - floor_term
            ) / (sorption_C1 + 1.0)
        return dH_vap_hexane * W + sorption_C0 * excess_integral

    @njit(cache=True)
    def _energy_cascade_kernel(
        initial_w: np.ndarray,
        initial_T: np.ndarray,
        particle_w: np.ndarray,
        particle_T: np.ndarray,
        vapor_T: np.ndarray,
        axial_sources: np.ndarray,
        moisture: np.ndarray,
        volumes: np.ndarray,
        face_areas: np.ndarray,
        dr: float,
        cell_dt: np.ndarray,
        time_substeps: np.ndarray,
        hQ: float,
        alpha_pg: float,
        rho_pg: float,
        cp_pg: float,
        alpha_ps: float,
        rho_ps: float,
        cp_ps: float,
        X3: float,
        cp_water_liquid: float,
        k_pg: float,
        k_ps: float,
        gab_Xm: float,
        gab_C0: float,
        gab_dHC_R: float,
        gab_K0: float,
        gab_dHK_R: float,
        oil_A0: float,
        oil_B: float,
        dH_vap_hexane: float,
        sorption_C0: float,
        sorption_C1: float,
        minimum_T: float,
        maximum_T: float,
        moisture_in: float,
        water_enthalpy_reference_k: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        nz, layers = particle_w.shape
        output_T = np.empty((nz, layers), dtype=np.float64)
        # Per representative particle and cell:
        # accumulation, sorption source, axial source, surface loss [W].
        energy_terms = np.zeros((nz, 4), dtype=np.float64)
        k_mix = alpha_pg * k_pg + alpha_ps * k_ps

        for j in range(nz):
            dt = cell_dt[j]
            substeps = time_substeps[j]
            particle_dt = dt / substeps
            entering_T = np.empty(layers, dtype=np.float64)
            sorption_sources = np.empty(layers, dtype=np.float64)
            for i in range(layers):
                entering_T[i] = initial_T[i] if j == 0 else output_T[j - 1, i]
                entering_w = initial_w[i] if j == 0 else particle_w[j - 1, i]
                outgoing_w = particle_w[j, i]
                temperature = particle_T[j, i]
                W_in = _gab_loading(
                    entering_w,
                    temperature,
                    gab_Xm,
                    gab_C0,
                    gab_dHC_R,
                    gab_K0,
                    gab_dHK_R,
                )
                W_out = _gab_loading(
                    outgoing_w,
                    temperature,
                    gab_Xm,
                    gab_C0,
                    gab_dHC_R,
                    gab_K0,
                    gab_dHK_R,
                )
                q_in = oil_A0 * entering_w**oil_B
                q_out = oil_A0 * outgoing_w**oil_B
                delta_h = (
                    _sorption_enthalpy_potential(
                        W_out,
                        gab_Xm,
                        dH_vap_hexane,
                        sorption_C0,
                        sorption_C1,
                    )
                    - _sorption_enthalpy_potential(
                        W_in,
                        gab_Xm,
                        dH_vap_hexane,
                        sorption_C0,
                        sorption_C1,
                    )
                    + X3 * dH_vap_hexane * (q_out - q_in)
                )
                sorption_sources[i] = alpha_ps * rho_ps * delta_h / dt
                energy_terms[j, 1] += sorption_sources[i] * volumes[i]
            moisture_start = moisture_in if j == 0 else moisture[j - 1]
            for step in range(substeps):
                moisture_old = moisture_start + (
                    moisture[j] - moisture_start
                ) * step / substeps
                moisture_new = moisture_start + (
                    moisture[j] - moisture_start
                ) * (step + 1) / substeps
                Cv_old = (
                    alpha_pg * rho_pg * cp_pg
                    + alpha_ps * rho_ps * cp_ps
                    + alpha_ps * rho_ps * moisture_old * cp_water_liquid
                )
                Cv_new = (
                    alpha_pg * rho_pg * cp_pg
                    + alpha_ps * rho_ps * cp_ps
                    + alpha_ps * rho_ps * moisture_new * cp_water_liquid
                )
                sub = np.zeros(layers, dtype=np.float64)
                diagonal = np.empty(layers, dtype=np.float64)
                sup = np.zeros(layers, dtype=np.float64)
                rhs = np.empty(layers, dtype=np.float64)

                for i in range(layers):
                    if step == 0:
                        old_T = initial_T[i] if j == 0 else output_T[j - 1, i]
                    else:
                        old_T = output_T[j, i]
                    source = sorption_sources[i] + axial_sources[j]

                    volume = volumes[i]
                    diag_i = Cv_new * volume / particle_dt
                    rhs_i = (
                        Cv_old * volume / particle_dt * old_T
                        + (Cv_new - Cv_old)
                        * volume
                        / particle_dt
                        * water_enthalpy_reference_k
                        + source * volume
                    )
                    if i > 0:
                        coeff_in = k_mix * face_areas[i - 1] / dr
                        sub[i] = -coeff_in
                        diag_i += coeff_in
                    if i < layers - 1:
                        coeff_out = k_mix * face_areas[i] / dr
                        sup[i] = -coeff_out
                        diag_i += coeff_out
                    else:
                        coeff_surface = hQ * face_areas[layers - 1]
                        diag_i += coeff_surface
                        rhs_i += coeff_surface * vapor_T[j]
                    diagonal[i] = diag_i
                    rhs[i] = rhs_i

                cp = np.zeros(layers, dtype=np.float64)
                dp = np.empty(layers, dtype=np.float64)
                beta = diagonal[0]
                cp[0] = sup[0] / beta
                dp[0] = rhs[0] / beta
                for i in range(1, layers):
                    beta = diagonal[i] - sub[i] * cp[i - 1]
                    if i < layers - 1:
                        cp[i] = sup[i] / beta
                    dp[i] = (rhs[i] - sub[i] * dp[i - 1]) / beta
                output_T[j, layers - 1] = dp[layers - 1]
                for i in range(layers - 2, -1, -1):
                    output_T[j, i] = dp[i] - cp[i] * output_T[j, i + 1]
                for i in range(layers):
                    output_T[j, i] = min(max(output_T[j, i], minimum_T), maximum_T)
                energy_terms[j, 3] += (
                    hQ
                    * face_areas[layers - 1]
                    * (output_T[j, layers - 1] - vapor_T[j])
                    / substeps
                )
            accumulation = 0.0
            axial = 0.0
            base_cv = alpha_pg * rho_pg * cp_pg + alpha_ps * rho_ps * cp_ps
            for i in range(layers):
                accumulation += (
                    (
                        (
                            base_cv
                            + alpha_ps * rho_ps * moisture[j] * cp_water_liquid
                        )
                        * (output_T[j, i] - water_enthalpy_reference_k)
                        - (
                            base_cv
                            + alpha_ps * rho_ps * moisture_start * cp_water_liquid
                        )
                        * (entering_T[i] - water_enthalpy_reference_k)
                    )
                    * volumes[i]
                    / dt
                )
                axial += axial_sources[j] * volumes[i]
            energy_terms[j, 0] = accumulation
            energy_terms[j, 2] = axial
        return output_T, energy_terms

    @njit(cache=True)
    def _mass_cascade_kernel(
        initial_w: np.ndarray,
        particle_T: np.ndarray,
        vapor_w: np.ndarray,
        volumes: np.ndarray,
        face_areas: np.ndarray,
        dr: float,
        cell_dt: np.ndarray,
        time_substeps: np.ndarray,
        hM: float,
        rho_V: float,
        D_eff: float,
        alpha_pg: float,
        rho_pg: float,
        alpha_ps: float,
        rho_ps: float,
        X3: float,
        gab_Xm: float,
        gab_C0: float,
        gab_dHC_R: float,
        gab_K0: float,
        gab_dHK_R: float,
        oil_A0: float,
        oil_B: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        nz, layers = particle_T.shape
        output_w = np.empty((nz, layers), dtype=np.float64)
        rates = np.empty((nz, layers), dtype=np.float64)
        x2_bulk = np.empty(nz, dtype=np.float64)
        total_volume = 0.0
        for i in range(layers):
            total_volume += volumes[i]
        const_diff = alpha_pg * rho_pg * D_eff

        for j in range(nz):
            dt = cell_dt[j]
            substeps = time_substeps[j]
            particle_dt = dt / substeps
            for step in range(substeps):
                sub = np.zeros(layers, dtype=np.float64)
                diagonal = np.empty(layers, dtype=np.float64)
                sup = np.zeros(layers, dtype=np.float64)
                rhs = np.empty(layers, dtype=np.float64)

                for i in range(layers):
                    if step == 0:
                        old_w = initial_w[i] if j == 0 else output_w[j - 1, i]
                    else:
                        old_w = output_w[j, i]
                    temperature = particle_T[j, i]
                    C = gab_C0 * np.exp(gab_dHC_R / temperature)
                    K = gab_K0 * np.exp(gab_dHK_R / temperature)
                    activity = old_w
                    if K > 0.0 and K * activity > 0.999:
                        activity = 0.999 / K
                    x = K * activity
                    denominator = (1.0 - x) * (1.0 - x + C * x)
                    gab_slope = (
                        gab_Xm
                        * C
                        * K
                        * (1.0 + (C - 1.0) * x * x)
                        / (denominator * denominator)
                    )
                    oil_slope = (
                        oil_A0 * oil_B * old_w ** (oil_B - 1.0)
                        if old_w > 0.0
                        else 0.0
                    )
                    slope = gab_slope + X3 * oil_slope
                    accumulation = alpha_pg * rho_pg + alpha_ps * rho_ps * slope
                    accumulation_volume_rate = volumes[i] * accumulation / particle_dt
                    diag_i = accumulation_volume_rate
                    rhs[i] = accumulation_volume_rate * old_w
                    if i > 0:
                        coeff_in = const_diff * face_areas[i - 1] / dr
                        sub[i] = -coeff_in
                        diag_i += coeff_in
                    if i < layers - 1:
                        coeff_out = const_diff * face_areas[i] / dr
                        sup[i] = -coeff_out
                        diag_i += coeff_out
                    else:
                        coeff_surface = hM * rho_V * face_areas[layers - 1]
                        diag_i += coeff_surface
                        rhs[i] += coeff_surface * vapor_w[j]
                    diagonal[i] = diag_i

                # Thomas solve, identical operation order to the Python reference.
                cp = np.zeros(layers, dtype=np.float64)
                dp = np.empty(layers, dtype=np.float64)
                beta = diagonal[0]
                cp[0] = sup[0] / beta
                dp[0] = rhs[0] / beta
                for i in range(1, layers):
                    beta = diagonal[i] - sub[i] * cp[i - 1]
                    if i < layers - 1:
                        cp[i] = sup[i] / beta
                    dp[i] = (rhs[i] - sub[i] * dp[i - 1]) / beta
                output_w[j, layers - 1] = dp[layers - 1]
                for i in range(layers - 2, -1, -1):
                    output_w[j, i] = dp[i] - cp[i] * output_w[j, i + 1]
                for i in range(layers):
                    value = output_w[j, i]
                    if value < 0.0:
                        value = 0.0
                    elif value > 1.0:
                        value = 1.0
                    output_w[j, i] = value

            bulk_sum = 0.0
            for i in range(layers):
                entering_w = initial_w[i] if j == 0 else output_w[j - 1, i]
                value = output_w[j, i]
                rates[j, i] = (value - entering_w) / dt

                temperature = particle_T[j, i]
                C = gab_C0 * np.exp(gab_dHC_R / temperature)
                K = gab_K0 * np.exp(gab_dHK_R / temperature)
                activity = value
                if K > 0.0 and K * activity > 0.999:
                    activity = 0.999 / K
                gab_value = (
                    gab_Xm
                    * C
                    * K
                    * activity
                    / ((1.0 - K * activity) * (1.0 - K * activity + C * K * activity))
                )
                oil_value = oil_A0 * value**oil_B
                bulk_sum += (gab_value + X3 * oil_value) * volumes[i]
            x2_bulk[j] = bulk_sum / total_volume
        return output_w, rates, x2_bulk


def _cell_step_arrays(
    count: int,
    dt: float | tuple[float, ...],
    time_substeps: int | tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray]:
    cell_dt = (
        np.full(count, dt, dtype=np.float64)
        if isinstance(dt, (int, float))
        else np.asarray(dt, dtype=np.float64)
    )
    substeps = (
        np.full(count, time_substeps, dtype=np.int64)
        if isinstance(time_substeps, int)
        else np.asarray(time_substeps, dtype=np.int64)
    )
    if len(cell_dt) != count or len(substeps) != count:
        raise ValueError("particle time-step profiles must match the axial cell count")
    if np.any(~np.isfinite(cell_dt)) or np.any(cell_dt <= 0.0) or np.any(substeps < 1):
        raise ValueError("particle time steps and substep counts must be positive")
    return cell_dt, substeps


def energy_cascade(
    initial_particle: "ParticleState",
    particles: list["ParticleState"],
    previous_rates: list[tuple[float, ...]],
    vapor_T: list[float],
    axial_sources: tuple[float, ...],
    moisture: list[float],
    hQ: float,
    dt: float | tuple[float, ...],
    geometry: "ShellGeometry",
    constants: "ParticleConstants",
    minimum_T: float,
    maximum_T: float,
    invariants: CascadeInvariants | None = None,
    time_substeps: int | tuple[int, ...] = 1,
    moisture_in: float = 0.0,
    water_enthalpy_reference_k: float = 373.15,
) -> np.ndarray | None:
    """Return compiled particle temperatures, or ``None`` when unavailable.

    ``previous_rates`` is retained for API compatibility; the conservative
    finite sorption-enthalpy increment is evaluated from adjacent cell states.
    """
    if not NUMBA_AVAILABLE or JIT_DISABLED:
        return None
    cell_dt, substeps = _cell_step_arrays(len(particles), dt, time_substeps)
    cached = invariants or build_invariants(initial_particle, geometry)
    result = _energy_cascade_kernel(
        cached.initial_w,
        cached.initial_T,
        np.asarray([particle.wpg2 for particle in particles], dtype=np.float64),
        np.asarray([particle.Tp for particle in particles], dtype=np.float64),
        np.asarray(vapor_T, dtype=np.float64),
        np.asarray(axial_sources, dtype=np.float64),
        np.asarray(moisture, dtype=np.float64),
        cached.volumes,
        cached.face_areas,
        geometry.dr,
        cell_dt,
        substeps,
        hQ,
        constants.alpha_pg,
        constants.rho_pg,
        constants.cp_pg,
        constants.alpha_ps,
        constants.rho_ps,
        constants.cp_ps,
        constants.X3,
        constants.cp_water_liquid,
        constants.k_pg,
        constants.k_ps,
        constants.gab.Xm,
        constants.gab.C0,
        constants.gab.dHC_R,
        constants.gab.K0,
        constants.gab.dHK_R,
        constants.oil.A0,
        constants.oil.B,
        constants.dH_vap_hexane,
        constants.sorption_C0,
        constants.sorption_C1,
        minimum_T,
        maximum_T,
        moisture_in,
        water_enthalpy_reference_k,
    )
    return result[0]


def energy_cascade_with_diagnostics(
    initial_particle: "ParticleState",
    particles: list["ParticleState"],
    previous_rates: list[tuple[float, ...]],
    vapor_T: list[float],
    axial_sources: tuple[float, ...],
    moisture: list[float],
    hQ: float,
    dt: float | tuple[float, ...],
    geometry: "ShellGeometry",
    constants: "ParticleConstants",
    minimum_T: float,
    maximum_T: float,
    invariants: CascadeInvariants | None = None,
    time_substeps: int | tuple[int, ...] = 1,
    moisture_in: float = 0.0,
    water_enthalpy_reference_k: float = 373.15,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Return temperatures and exact discrete particle-energy terms.

    Diagnostic columns are accumulation, sorption source, axial source, and
    surface heat loss, all in W per representative particle and averaged
    over the full axial residence step.
    """
    if not NUMBA_AVAILABLE or JIT_DISABLED:
        return None
    cell_dt, substeps = _cell_step_arrays(len(particles), dt, time_substeps)
    cached = invariants or build_invariants(initial_particle, geometry)
    return _energy_cascade_kernel(
        cached.initial_w,
        cached.initial_T,
        np.asarray([particle.wpg2 for particle in particles], dtype=np.float64),
        np.asarray([particle.Tp for particle in particles], dtype=np.float64),
        np.asarray(vapor_T, dtype=np.float64),
        np.asarray(axial_sources, dtype=np.float64),
        np.asarray(moisture, dtype=np.float64),
        cached.volumes,
        cached.face_areas,
        geometry.dr,
        cell_dt,
        substeps,
        hQ,
        constants.alpha_pg,
        constants.rho_pg,
        constants.cp_pg,
        constants.alpha_ps,
        constants.rho_ps,
        constants.cp_ps,
        constants.X3,
        constants.cp_water_liquid,
        constants.k_pg,
        constants.k_ps,
        constants.gab.Xm,
        constants.gab.C0,
        constants.gab.dHC_R,
        constants.gab.K0,
        constants.gab.dHK_R,
        constants.oil.A0,
        constants.oil.B,
        constants.dH_vap_hexane,
        constants.sorption_C0,
        constants.sorption_C1,
        minimum_T,
        maximum_T,
        moisture_in,
        water_enthalpy_reference_k,
    )


def mass_cascade(
    initial_particle: "ParticleState",
    particles: list["ParticleState"],
    vapor_wV2: list[float],
    hM: float,
    rho_V: float,
    dt: float | tuple[float, ...],
    geometry: "ShellGeometry",
    constants: "ParticleConstants",
    invariants: CascadeInvariants | None = None,
    time_substeps: int | tuple[int, ...] = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """The frozen-storage compiled kernel is retired; use the exact Python map."""

    return None


def warm_up(
    initial_particle: "ParticleState",
    geometry: "ShellGeometry",
    constants: "ParticleConstants",
) -> bool:
    """Compile the backend once during initialization; return whether it is active."""
    result = mass_cascade(
        initial_particle,
        [initial_particle],
        [0.5],
        1.0e-3,
        1.0,
        1.0,
        geometry,
        constants,
    )
    return result is not None
