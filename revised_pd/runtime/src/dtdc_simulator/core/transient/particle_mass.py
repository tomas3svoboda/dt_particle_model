"""Exactly conservative transient radial hexane balance for a DCZ particle.

This is the first verification slice of the transient reference model.  It
uses Coletto's accumulated quantity from eq. A.22 directly,

``C_2 = alpha_pg*rho_pg*a_h + alpha_ps*rho_ps*X2_so(a_h, T)``,

as the nonlinear backward-Euler storage term.  Internal radial fluxes are
single shared face rates, and the outer Robin rate is the same quantity used
by the particle and its external ledger.

The production steady solver's historical particle march linearizes the
storage term.  That is appropriate as an efficient fixed-point map but leaves
a finite-step component mismatch.  A transient reference step cannot hide
that mismatch, so this module solves the exact nonlinear inventory residual
with an analytic Jacobian and a bounded Newton line search.

Temperature is fixed during this first slice.  The coupled enthalpy state is a
later milestone; naming the limitation explicitly prevents this verified mass
kernel from being mistaken for the complete transient DCZ.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from dtdc_simulator.core import thermo
from dtdc_simulator.core.transient.event_schedule import event_aligned_steps
from dtdc_simulator.core.zones import particle as pt


@dataclass(frozen=True)
class ParticleMassNumerics:
    """Numerical controls for one nonlinear backward-Euler step."""

    relative_tolerance: float = 1.0e-10
    absolute_tolerance_kg_s: float = 1.0e-22
    maximum_newton_iterations: int = 30
    minimum_line_search: float = 2.0**-16

    def __post_init__(self) -> None:
        if self.relative_tolerance <= 0.0:
            raise ValueError("relative tolerance must be positive")
        if self.absolute_tolerance_kg_s <= 0.0:
            raise ValueError("absolute tolerance must be positive")
        if self.maximum_newton_iterations <= 0:
            raise ValueError("maximum Newton iterations must be positive")
        if not 0.0 < self.minimum_line_search <= 1.0:
            raise ValueError("minimum line-search fraction must be in (0, 1]")


@dataclass(frozen=True)
class ParticleMassStepDiagnostics:
    """Executable component ledger for one accepted particle step."""

    newton_iterations: int
    nonlinear_residual_max_kg_s: float
    shell_ledger_residual_kg: tuple[float, ...]
    internal_face_flux_kg_s: tuple[float, ...]
    surface_flux_to_vapor_kg_s: float
    total_inventory_before_kg: float
    total_inventory_after_kg: float
    total_ledger_residual_kg: float


@dataclass(frozen=True)
class BoundaryEvent:
    """Piecewise-constant external vapor activity applied from ``time_s``."""

    time_s: float
    vapor_hexane_mass_fraction: float

    def __post_init__(self) -> None:
        if self.time_s < 0.0:
            raise ValueError("boundary-event time must be non-negative")
        if not 0.0 <= self.vapor_hexane_mass_fraction <= 1.0:
            raise ValueError("vapor hexane mass fraction must be in [0, 1]")


@dataclass(frozen=True)
class ParticleMassSnapshot:
    time_s: float
    state: pt.ParticleState
    step_boundary_wV2: float | None
    step_ledger_residual_kg: float


@dataclass(frozen=True)
class ParticleMassTrajectory:
    snapshots: tuple[ParticleMassSnapshot, ...]
    cumulative_ledger_residual_kg: float
    accepted_steps: int
    event_times_s: tuple[float, ...]

    @property
    def final_state(self) -> pt.ParticleState:
        return self.snapshots[-1].state


def _validate_geometry(
    state: pt.ParticleState,
    geometry: pt.ShellGeometry,
    c: pt.ParticleConstants,
) -> None:
    if c.Np <= 0:
        raise ValueError("particle shell count must be positive")
    if not (
        len(state.wpg2)
        == len(state.Tp)
        == len(geometry.volumes)
        == len(geometry.face_areas)
        == c.Np
    ):
        raise ValueError("particle state, geometry, and constants disagree")
    if any(not 0.0 <= value <= 1.0 for value in state.wpg2):
        raise ValueError("particle activity must remain in [0, 1]")
    if any(value <= 0.0 or not np.isfinite(value) for value in state.Tp):
        raise ValueError("particle temperature must be finite and positive")


def _inventory_density_and_slope(
    activity: float,
    temperature_k: float,
    c: pt.ParticleConstants,
) -> tuple[float, float]:
    x2_so, slope = thermo.x2_so_and_slope(
        activity,
        temperature_k,
        c.X3,
        c.gab,
        c.oil,
    )
    density = (
        c.alpha_pg * c.rho_pg * activity
        + c.alpha_ps * c.rho_ps * x2_so
    )
    derivative = c.alpha_pg * c.rho_pg + c.alpha_ps * c.rho_ps * slope
    return density, derivative


def particle_hexane_inventory_kg(
    state: pt.ParticleState,
    geometry: pt.ShellGeometry,
    c: pt.ParticleConstants,
) -> float:
    """Total free plus sorbed/absorbed hexane in one representative particle."""

    _validate_geometry(state, geometry, c)
    return sum(
        _inventory_density_and_slope(activity, temperature, c)[0] * volume
        for activity, temperature, volume in zip(
            state.wpg2,
            state.Tp,
            geometry.volumes,
        )
    )


def _residual_and_tridiagonal_jacobian(
    activity: np.ndarray,
    old_density: np.ndarray,
    temperature: np.ndarray,
    *,
    vapor_wV2: float,
    hM: float,
    rho_vapor: float,
    dt_s: float,
    geometry: pt.ShellGeometry,
    c: pt.ParticleConstants,
) -> tuple[
    np.ndarray,
    tuple[np.ndarray, np.ndarray, np.ndarray],
    tuple[float, ...],
    float,
]:
    shells = c.Np
    residual = np.zeros(shells, dtype=float)
    lower = np.zeros(shells, dtype=float)
    diagonal = np.zeros(shells, dtype=float)
    upper = np.zeros(shells, dtype=float)
    internal_flux = [0.0] * max(shells - 1, 0)
    diffusion = c.alpha_pg * c.rho_pg * c.D_eff

    for face in range(shells - 1):
        conductance = diffusion * geometry.face_areas[face] / geometry.dr
        internal_flux[face] = conductance * (
            activity[face] - activity[face + 1]
        )

    surface_conductance = hM * rho_vapor * geometry.face_areas[-1]
    surface_flux = surface_conductance * (activity[-1] - vapor_wV2)

    for shell in range(shells):
        density, density_slope = _inventory_density_and_slope(
            float(activity[shell]),
            float(temperature[shell]),
            c,
        )
        volume = geometry.volumes[shell]
        residual[shell] = volume * (density - old_density[shell]) / dt_s
        diagonal[shell] = volume * density_slope / dt_s

        if shell > 0:
            conductance = (
                diffusion * geometry.face_areas[shell - 1] / geometry.dr
            )
            residual[shell] -= internal_flux[shell - 1]
            lower[shell] -= conductance
            diagonal[shell] += conductance
        if shell < shells - 1:
            conductance = diffusion * geometry.face_areas[shell] / geometry.dr
            residual[shell] += internal_flux[shell]
            diagonal[shell] += conductance
            upper[shell] -= conductance
        else:
            residual[shell] += surface_flux
            diagonal[shell] += surface_conductance

    return (
        residual,
        (lower, diagonal, upper),
        tuple(internal_flux),
        float(surface_flux),
    )


def _solve_tridiagonal(
    lower: np.ndarray,
    diagonal: np.ndarray,
    upper: np.ndarray,
    rhs: np.ndarray,
) -> np.ndarray:
    """Thomas solve for the radial finite-volume Newton matrix."""

    count = len(rhs)
    modified_upper = np.empty(count, dtype=float)
    modified_rhs = np.empty(count, dtype=float)
    pivot = diagonal[0]
    if pivot == 0.0:
        raise ValueError("particle transient Jacobian is singular")
    modified_upper[0] = upper[0] / pivot
    modified_rhs[0] = rhs[0] / pivot
    for index in range(1, count):
        pivot = diagonal[index] - lower[index] * modified_upper[index - 1]
        if pivot == 0.0:
            raise ValueError("particle transient Jacobian is singular")
        modified_upper[index] = (
            upper[index] / pivot if index < count - 1 else 0.0
        )
        modified_rhs[index] = (
            rhs[index] - lower[index] * modified_rhs[index - 1]
        ) / pivot
    solution = np.empty(count, dtype=float)
    solution[-1] = modified_rhs[-1]
    for index in range(count - 2, -1, -1):
        solution[index] = (
            modified_rhs[index]
            - modified_upper[index] * solution[index + 1]
        )
    return solution


def step_isothermal_particle_hexane(
    state: pt.ParticleState,
    *,
    vapor_wV2: float,
    hM: float,
    rho_vapor: float,
    dt_s: float,
    geometry: pt.ShellGeometry,
    c: pt.ParticleConstants,
    numerics: ParticleMassNumerics = ParticleMassNumerics(),
) -> tuple[pt.ParticleState, ParticleMassStepDiagnostics]:
    """Advance the exact isothermal particle inventory by one BE step."""

    _validate_geometry(state, geometry, c)
    if not 0.0 <= vapor_wV2 <= 1.0:
        raise ValueError("external vapor hexane fraction must be in [0, 1]")
    if hM < 0.0 or rho_vapor <= 0.0 or dt_s <= 0.0:
        raise ValueError("hM must be non-negative; density and step positive")

    old_activity = np.asarray(state.wpg2, dtype=float)
    temperature = np.asarray(state.Tp, dtype=float)
    old_density = np.asarray(
        [
            _inventory_density_and_slope(activity, temp, c)[0]
            for activity, temp in zip(old_activity, temperature)
        ],
        dtype=float,
    )
    activity = old_activity.copy()

    residual = np.empty(c.Np, dtype=float)
    internal_flux: tuple[float, ...] = ()
    surface_flux = 0.0
    converged = False
    iterations = 0
    for iterations in range(1, numerics.maximum_newton_iterations + 1):
        residual, jacobian, internal_flux, surface_flux = (
            _residual_and_tridiagonal_jacobian(
                activity,
                old_density,
                temperature,
                vapor_wV2=vapor_wV2,
                hM=hM,
                rho_vapor=rho_vapor,
                dt_s=dt_s,
                geometry=geometry,
                c=c,
            )
        )
        residual_norm = float(np.max(np.abs(residual)))
        scale = max(
            float(np.max(np.abs(old_density * np.asarray(geometry.volumes))))
            / dt_s,
            abs(surface_flux),
            numerics.absolute_tolerance_kg_s,
        )
        tolerance = (
            numerics.absolute_tolerance_kg_s
            + numerics.relative_tolerance * scale
        )
        if residual_norm <= tolerance:
            converged = True
            break

        direction = _solve_tridiagonal(*jacobian, -residual)

        fraction = 1.0
        accepted = False
        while fraction >= numerics.minimum_line_search:
            trial = activity + fraction * direction
            if np.all((trial >= 0.0) & (trial <= 1.0)):
                trial_residual, _, _, _ = (
                    _residual_and_tridiagonal_jacobian(
                    trial,
                    old_density,
                    temperature,
                    vapor_wV2=vapor_wV2,
                    hM=hM,
                    rho_vapor=rho_vapor,
                    dt_s=dt_s,
                    geometry=geometry,
                    c=c,
                    )
                )
                if float(np.max(np.abs(trial_residual))) < residual_norm:
                    activity = trial
                    accepted = True
                    break
            fraction *= 0.5
        if not accepted:
            raise ValueError(
                "particle transient Newton line search failed within [0, 1]"
            )

    if not converged:
        raise ValueError(
            "particle transient Newton solve exceeded its iteration limit"
        )

    new_state = pt.ParticleState(
        wpg2=tuple(float(value) for value in activity),
        Tp=state.Tp,
    )
    new_density = np.asarray(
        [
            _inventory_density_and_slope(activity_i, temp, c)[0]
            for activity_i, temp in zip(activity, temperature)
        ],
        dtype=float,
    )
    inventory_change = (
        np.asarray(geometry.volumes) * (new_density - old_density)
    )
    shell_ledger = np.zeros(c.Np, dtype=float)
    for shell in range(c.Np):
        flux_in = internal_flux[shell - 1] if shell > 0 else 0.0
        flux_out = (
            internal_flux[shell]
            if shell < c.Np - 1
            else surface_flux
        )
        shell_ledger[shell] = (
            inventory_change[shell] + dt_s * (flux_out - flux_in)
        )

    total_before = float(
        np.sum(np.asarray(geometry.volumes) * old_density)
    )
    total_after = float(
        np.sum(np.asarray(geometry.volumes) * new_density)
    )
    total_ledger = total_after - total_before + dt_s * surface_flux
    diagnostics = ParticleMassStepDiagnostics(
        newton_iterations=iterations,
        nonlinear_residual_max_kg_s=float(np.max(np.abs(residual))),
        shell_ledger_residual_kg=tuple(float(value) for value in shell_ledger),
        internal_face_flux_kg_s=internal_flux,
        surface_flux_to_vapor_kg_s=surface_flux,
        total_inventory_before_kg=total_before,
        total_inventory_after_kg=total_after,
        total_ledger_residual_kg=float(total_ledger),
    )
    return new_state, diagnostics


def integrate_isothermal_particle_hexane(
    initial_state: pt.ParticleState,
    *,
    end_time_s: float,
    maximum_step_s: float,
    boundary_events: tuple[BoundaryEvent, ...],
    hM: float,
    rho_vapor: float,
    geometry: pt.ShellGeometry,
    c: pt.ParticleConstants,
    numerics: ParticleMassNumerics = ParticleMassNumerics(),
) -> ParticleMassTrajectory:
    """Integrate without allowing a numerical step to straddle an event."""

    _validate_geometry(initial_state, geometry, c)
    snapshots = [
        ParticleMassSnapshot(
            time_s=0.0,
            state=initial_state,
            step_boundary_wV2=None,
            step_ledger_residual_kg=0.0,
        )
    ]
    state = initial_state
    cumulative_ledger = 0.0
    accepted_steps = 0

    for step_end, dt_s, event in event_aligned_steps(
        end_time_s=end_time_s,
        maximum_step_s=maximum_step_s,
        events=boundary_events,
        event_name="boundary",
        horizon_name="the integration horizon",
    ):
        boundary = event.vapor_hexane_mass_fraction
        state, diagnostics = step_isothermal_particle_hexane(
            state,
            vapor_wV2=boundary,
            hM=hM,
            rho_vapor=rho_vapor,
            dt_s=dt_s,
            geometry=geometry,
            c=c,
            numerics=numerics,
        )
        accepted_steps += 1
        cumulative_ledger += diagnostics.total_ledger_residual_kg
        snapshots.append(
            ParticleMassSnapshot(
                time_s=step_end,
                state=state,
                step_boundary_wV2=boundary,
                step_ledger_residual_kg=diagnostics.total_ledger_residual_kg,
            )
        )

    return ParticleMassTrajectory(
        snapshots=tuple(snapshots),
        cumulative_ledger_residual_kg=float(cumulative_ledger),
        accepted_steps=accepted_steps,
        event_times_s=tuple(event.time_s for event in boundary_events),
    )
