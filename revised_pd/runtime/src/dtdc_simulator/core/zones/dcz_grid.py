"""Conservative axial-grid utilities for the DCZ."""

from __future__ import annotations

import math

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

from dtdc_simulator.core.zones import particle as pt
from dtdc_simulator.core.zones.dcz_types import DCZWarmStart


def _cell_centres(
    widths: tuple[float, ...],
    count: int,
    height_m: float,
) -> np.ndarray:
    if not widths:
        return (np.arange(count, dtype=float) + 0.5) * height_m / count
    if len(widths) != count:
        raise ValueError("DCZ warm-start cell widths must match profile length")
    faces = np.concatenate(([0.0], np.cumsum(np.asarray(widths, dtype=float))))
    return 0.5 * (faces[:-1] + faces[1:])


def interpolate_profile(
    values: tuple[float, ...],
    source_height_m: float,
    target_n: int,
    target_height_m: float,
    source_widths_m: tuple[float, ...] = (),
    target_widths_m: tuple[float, ...] = (),
) -> list[float]:
    """Interpolate cell-centred warm data after a DCZ boundary remesh."""
    source_n = len(values)
    if source_n == 0:
        raise ValueError("DCZ warm-start profiles must be nonempty")
    source_z = _cell_centres(source_widths_m, source_n, source_height_m)
    target_z = _cell_centres(target_widths_m, target_n, target_height_m)
    return [
        float(value)
        for value in np.interp(
            target_z,
            source_z,
            np.asarray(values, dtype=float),
            left=values[0],
            right=values[-1],
        )
    ]


def resample_warm_start(
    warm: DCZWarmStart,
    target_n: int,
    target_height_m: float,
    particle_layers: int,
    target_widths_m: tuple[float, ...] = (),
) -> DCZWarmStart:
    """Preserve a lagged DCZ iterate across axial and radial remeshing."""
    if len(warm.particles) == 0:
        raise ValueError("DCZ warm-start particle profile must be nonempty")
    source_layers = len(warm.particles[0].wpg2)
    if source_layers == 0 or any(
        len(particle.wpg2) != source_layers
        or len(particle.Tp) != source_layers
        for particle in warm.particles
    ):
        raise ValueError("DCZ warm-start particle profiles must have a consistent shape")
    if any(len(rates) != source_layers for rates in warm.dwpg2_dt):
        raise ValueError("DCZ warm-start particle rates must match the particle shape")

    def particle_component(name: str, layer: int) -> tuple[float, ...]:
        return tuple(getattr(particle, name)[layer] for particle in warm.particles)

    def interpolate(values: tuple[float, ...]) -> list[float]:
        return interpolate_profile(
            values,
            warm.bed_height_m,
            target_n,
            target_height_m,
            warm.cell_widths_m,
            target_widths_m,
        )

    wpg_axial = [
        interpolate(particle_component("wpg2", layer))
        for layer in range(source_layers)
    ]
    temperature_axial = [
        interpolate(particle_component("Tp", layer))
        for layer in range(source_layers)
    ]
    rate_axial = [
        interpolate(tuple(rates[layer] for rates in warm.dwpg2_dt))
        for layer in range(source_layers)
    ]

    source_r = (np.arange(source_layers, dtype=float) + 0.5) / source_layers
    target_r = (np.arange(particle_layers, dtype=float) + 0.5) / particle_layers

    def radial_profile(
        axial_layers: list[list[float]],
        cell: int,
    ) -> tuple[float, ...]:
        values = np.asarray(
            [axial_layers[layer][cell] for layer in range(source_layers)],
            dtype=float,
        )
        return tuple(
            float(value)
            for value in np.interp(
                target_r,
                source_r,
                values,
                left=values[0],
                right=values[-1],
            )
        )

    particles = tuple(
        pt.ParticleState(
            wpg2=radial_profile(wpg_axial, cell),
            Tp=radial_profile(temperature_axial, cell),
        )
        for cell in range(target_n)
    )
    return DCZWarmStart(
        bed_height_m=target_height_m,
        particles=particles,
        vapor_T=tuple(interpolate(warm.vapor_T)),
        vapor_water_kg_s=tuple(interpolate(warm.vapor_water_kg_s)),
        vapor_hexane_kg_s=tuple(interpolate(warm.vapor_hexane_kg_s)),
        X1=tuple(interpolate(warm.X1)),
        free_water=tuple(interpolate(warm.free_water)),
        condensed_water_kg_s=tuple(interpolate(warm.condensed_water_kg_s)),
        evaporated_water_kg_s=tuple(interpolate(warm.evaporated_water_kg_s)),
        bound_water_mass_rate_w_m3=tuple(
            interpolate(warm.bound_water_mass_rate_w_m3)
        ),
        dwpg2_dt=tuple(
            radial_profile(rate_axial, cell) for cell in range(target_n)
        ),
        temperature_relaxation=warm.temperature_relaxation,
        hexane_relaxation=warm.hexane_relaxation,
        water_relaxation=warm.water_relaxation,
        water_active_set=(
            warm.water_active_set
            if len(warm.particles) == target_n
            and source_layers == particle_layers
            else None
        ),
        cell_widths_m=target_widths_m,
    )


def axial_laplacian(
    profile: tuple[float, ...],
    dz: float | tuple[float, ...],
) -> tuple[float, ...]:
    """Conservative cell-centred Laplacian with zero-flux end faces."""
    n = len(profile)
    widths = tuple(dz for _ in range(n)) if isinstance(dz, (int, float)) else dz
    if len(widths) != n or any(
        not math.isfinite(value) or value <= 0.0 for value in widths
    ):
        raise ValueError(
            "axial cell widths must be finite, positive, and match the profile"
        )
    face_gradient = [0.0] * (n + 1)
    for face in range(1, n):
        centre_distance = 0.5 * (widths[face - 1] + widths[face])
        face_gradient[face] = (
            profile[face] - profile[face - 1]
        ) / centre_distance
    return tuple(
        (face_gradient[cell + 1] - face_gradient[cell]) / widths[cell]
        for cell in range(n)
    )


def solve_vapor_hexane_dispersion(
    *,
    vapor_water_kg_s: list[float] | tuple[float, ...],
    hexane_release_kg_s: list[float] | tuple[float, ...],
    bottom_hexane_kg_s: float,
    cell_widths_m: tuple[float, ...],
    bed_area_m2: float,
    alpha_v: float,
    axial_diffusivity_m2_s: float,
    vapor_density_kg_m3: float,
    initial_hexane_kg_s: list[float] | tuple[float, ...] | None = None,
    tolerance_kg_s: float = 1.0e-11,
    maximum_iterations: int = 20,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Solve the conservative counter-current hexane balance implicitly."""
    n = len(cell_widths_m)
    if (
        n == 0
        or len(vapor_water_kg_s) != n
        or len(hexane_release_kg_s) != n
    ):
        raise ValueError("hexane-dispersion profiles must have one value per cell")
    if any(
        not math.isfinite(value) or value <= 0.0 for value in cell_widths_m
    ):
        raise ValueError("hexane-dispersion cell widths must be finite and positive")
    if any(
        not math.isfinite(value) or value < 0.0 for value in vapor_water_kg_s
    ):
        raise ValueError("vapor water flows must be finite and non-negative")
    if (
        not math.isfinite(bottom_hexane_kg_s)
        or bottom_hexane_kg_s < 0.0
        or not math.isfinite(bed_area_m2)
        or bed_area_m2 <= 0.0
        or not math.isfinite(alpha_v)
        or alpha_v < 0.0
        or not math.isfinite(axial_diffusivity_m2_s)
        or axial_diffusivity_m2_s < 0.0
        or not math.isfinite(vapor_density_kg_m3)
        or vapor_density_kg_m3 < 0.0
        or not math.isfinite(tolerance_kg_s)
        or tolerance_kg_s <= 0.0
        or maximum_iterations <= 0
    ):
        raise ValueError("invalid hexane-dispersion solve parameters")

    lap_rows: list[int] = []
    lap_cols: list[int] = []
    lap_data: list[float] = []
    for face in range(1, n):
        left = face - 1
        right = face
        centre_distance = 0.5 * (
            cell_widths_m[left] + cell_widths_m[right]
        )
        left_coefficient = 1.0 / (
            cell_widths_m[left] * centre_distance
        )
        right_coefficient = 1.0 / (
            cell_widths_m[right] * centre_distance
        )
        lap_rows.extend((left, left, right, right))
        lap_cols.extend((left, right, left, right))
        lap_data.extend(
            (
                -left_coefficient,
                left_coefficient,
                right_coefficient,
                -right_coefficient,
            )
        )
    laplacian_matrix = coo_matrix(
        (lap_data, (lap_rows, lap_cols)),
        shape=(n, n),
    ).tocsc()

    base_matrix = diags(
        (np.ones(n), -np.ones(max(n - 1, 0))),
        offsets=(0, 1),
        shape=(n, n),
        format="csc",
    )
    coefficient_kg_s = (
        alpha_v
        * axial_diffusivity_m2_s
        * vapor_density_kg_m3
        * bed_area_m2
        * np.asarray(cell_widths_m, dtype=float)
    )
    water = np.asarray(vapor_water_kg_s, dtype=float)
    release = np.asarray(hexane_release_kg_s, dtype=float)

    if initial_hexane_kg_s is None:
        initial = np.empty(n, dtype=float)
        running = bottom_hexane_kg_s
        for cell in range(n - 1, -1, -1):
            running += release[cell]
            if running < 0.0:
                raise ValueError(
                    "particle hexane uptake exceeds the available vapor inventory"
                )
            initial[cell] = running
        hexane = initial
    else:
        if len(initial_hexane_kg_s) != n:
            raise ValueError("initial hexane profile must match the axial grid")
        hexane = np.maximum(
            np.asarray(initial_hexane_kg_s, dtype=float),
            0.0,
        )
        if not np.all(np.isfinite(hexane)):
            raise ValueError("initial hexane profile must be finite")

    bottom_vector = np.zeros(n, dtype=float)
    bottom_vector[-1] = bottom_hexane_kg_s

    def residual_at(candidate: np.ndarray) -> np.ndarray:
        total = water + candidate
        if np.any(total <= 0.0):
            return np.full(n, math.inf)
        mass_fraction = candidate / total
        return (
            np.asarray(base_matrix @ candidate, dtype=float)
            - bottom_vector
            - release
            - coefficient_kg_s
            * np.asarray(laplacian_matrix @ mass_fraction, dtype=float)
        )

    residual = residual_at(hexane)
    for _ in range(maximum_iterations):
        residual_norm = float(np.max(np.abs(residual)))
        if residual_norm <= tolerance_kg_s:
            break
        total = water + hexane
        mass_fraction_derivative = water / (total * total)
        jacobian = (
            base_matrix
            - diags(coefficient_kg_s, format="csc")
            @ laplacian_matrix
            @ diags(mass_fraction_derivative, format="csc")
        )
        step = np.asarray(spsolve(jacobian, -residual), dtype=float)
        if not np.all(np.isfinite(step)):
            raise ValueError("hexane-dispersion Newton step is non-finite")

        merit = residual_norm
        accepted = False
        step_fraction = 1.0
        for _ in range(40):
            candidate = hexane + step_fraction * step
            if np.min(candidate) >= 0.0:
                candidate_residual = residual_at(candidate)
                candidate_merit = float(
                    np.max(np.abs(candidate_residual))
                )
                if candidate_merit < merit:
                    hexane = candidate
                    residual = candidate_residual
                    accepted = True
                    break
            step_fraction *= 0.5
        if not accepted:
            raise ValueError(
                "hexane-dispersion Newton line search could not reduce the residual"
            )
    else:
        raise ValueError(
            "hexane-dispersion solve did not converge "
            f"(residual {float(np.max(np.abs(residual))):.6g} kg/s)"
        )

    mass_fraction = hexane / (water + hexane)
    axial_source_kg_m3_s = (
        alpha_v
        * axial_diffusivity_m2_s
        * vapor_density_kg_m3
        * np.asarray(
            axial_laplacian(tuple(mass_fraction), cell_widths_m),
            dtype=float,
        )
    )
    return (
        tuple(float(value) for value in hexane),
        tuple(float(value) for value in axial_source_kg_m3_s),
    )
