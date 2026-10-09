"""Conservative spherical finite-volume geometry for the resolved particle.

Cell-centered finite volumes on ``r_inner <= r <= R``. Face radii bound each
cell; cell "centers" are the geometric face midpoints (used for gradient
stencils). Cell volumes and face areas are the exact spherical measures, so a
flux-form update conserves the integrated quantity by construction
(telescoping face fluxes).

    d/dt( V_i C_i ) = A_{i-1/2} F_{i-1/2} - A_{i+1/2} F_{i+1/2} + V_i S_i

with F the radial flux (outward positive), A = 4*pi*r_face^2, and the innermost
left face at r=0 (A=0) giving automatic zero-flux symmetry.
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class SphericalGrid:
    R: float
    faces: tuple[float, ...]     # N+1 face radii, faces[0]>=0, faces[-1]=R
    centers: tuple[float, ...]   # N cell centers (face midpoints)
    volumes: tuple[float, ...]   # N cell volumes (4/3 pi (rR^3 - rL^3))
    areas: tuple[float, ...]     # N+1 face areas (4 pi r^2)

    @property
    def n(self) -> int:
        return len(self.centers)

    @property
    def inner_radius(self) -> float:
        return self.faces[0]

    @property
    def outer_radius(self) -> float:
        return self.R

    @property
    def total_volume(self) -> float:
        return 4.0 / 3.0 * math.pi * (self.R**3 - self.inner_radius**3)


def uniform_grid(n: int, R: float) -> SphericalGrid:
    """Uniform radial spacing in r (faces equally spaced on [0, R])."""
    return annular_grid(n, 0.0, R)


def annular_grid(n: int, inner_radius: float, outer_radius: float) -> SphericalGrid:
    """Uniform radial spacing on an exact, non-empty spherical annulus.

    ``inner_radius=0`` is the exact dry-sphere endpoint.  A zero-width annulus
    is rejected: the fully wet endpoint has no dry-region cells and must be
    represented by the regime state, never by an epsilon-thickness domain.
    """
    try:
        cell_count = operator.index(n)
    except TypeError as exc:
        raise ValueError("need an integer number of cells") from exc
    if isinstance(n, bool) or cell_count < 1:
        raise ValueError("need at least one cell")
    radii = (inner_radius, outer_radius)
    if not all(isinstance(radius, (int, float)) and math.isfinite(radius) for radius in radii):
        raise ValueError("grid radii must be finite real numbers")
    if inner_radius < 0.0:
        raise ValueError("inner radius must be non-negative")
    if outer_radius <= inner_radius:
        raise ValueError("outer radius must exceed inner radius; no epsilon domains")
    width = outer_radius - inner_radius
    face_list = [inner_radius + width * k / cell_count for k in range(cell_count + 1)]
    # Preserve endpoint identity for regime/event checks despite floating-point
    # evaluation of ``inner + width`` at k=n.
    face_list[0] = inner_radius
    face_list[-1] = outer_radius
    return _from_faces(tuple(face_list), outer_radius)


def _from_faces(faces: tuple[float, ...], R: float) -> SphericalGrid:
    if len(faces) < 2:
        raise ValueError("a spherical grid needs at least two faces")
    if not all(math.isfinite(face) for face in faces):
        raise ValueError("grid faces must be finite")
    if faces[0] < 0.0 or any(right <= left for left, right in zip(faces, faces[1:])):
        raise ValueError("grid faces must be non-negative and strictly increasing")
    if faces[-1] != R:
        raise ValueError("outer face must equal the declared particle radius")
    centers = tuple(0.5 * (faces[i] + faces[i + 1]) for i in range(len(faces) - 1))
    volumes = tuple(
        4.0 / 3.0 * math.pi * (faces[i + 1] ** 3 - faces[i] ** 3)
        for i in range(len(faces) - 1)
    )
    areas = tuple(4.0 * math.pi * rf * rf for rf in faces)
    return SphericalGrid(R=R, faces=faces, centers=centers, volumes=volumes, areas=areas)


def integrate(grid: SphericalGrid, cell_values: Sequence[float]) -> float:
    """Volume integral sum_i V_i * value_i (extensive total)."""
    if len(cell_values) != grid.n:
        raise ValueError("cell values must align with the grid")
    return sum(v * c for v, c in zip(grid.volumes, cell_values))
