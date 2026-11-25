"""Surface-graded radial face set for the Paper-1 particle kernel (debt D2, 2026-10-22).

The surface-layer mesh packet (docs/GT_PS2_PAPER1_SURFACE_LAYER_MESH_PACKET_2026-10-18.md,
evidence surface_layer_mesh_scoping_151) found the DT outlet moisture's mesh error in a thin
diffusion layer under the ruled film-free surface and scoped a radial mesh with a uniform
interior and a geometrically refined band under the surface. Owner ruling
GT_PS2_OWNER_RULING_UNIVERSAL_KERNEL_ENVELOPE_2026-10-19 item 7 cleared it to be run. This
module builds that face set; it is armed by the knob
``DTDC_KERNEL_ENVELOPE_SURFACE_GRADED_MESH`` (particle_kernel_envelope) together with the
runner's three declared parameters, and with the knob absent nothing here is called.

The declared parameters are the packet's: the interior cell count m, the band's cell count n
and the band's ratio q > 1. The faces are the image of the uniform computational grid
xi_i = i / (m + n) under one fixed C1 map of [0, 1] onto [0, R]:

* interior, xi <= xi_b = m / (m + n): r = s xi (uniform cells of width w = s / (m + n));
* band, xi >= xi_b: the depth below the surface d = R - r = A (q^((m + n)(1 - xi)) - 1),
  so the band's cells are exactly geometric, A (q - 1) q^j for j = 0 .. n - 1 counted from
  the surface;
* value and slope continuous at xi_b, which fixes A = R / (m q^n ln q + q^n - 1) and
  w = A q^n ln q (the interior width is the band's progression continued smoothly; the
  width ratio across the junction is q ln q / (q - 1), within (q - 1) / 2 of one).

The family at refinement ratio 2 (the order and GCI formulas of the convergence ladder
assume it): (m, n, q) -> (2 m, 2 n, q^(1/2)) leaves q^n, m ln q, A, w (m + n) and the map
itself unchanged, so every face of a rung is a face of the next finer rung (up to rounding)
and every cell is split in two, the interior cells at their midpoints and the band cells in
the band's own geometric proportion. The packet's own plan wrote "the first cell halved, the
ratio kept, the interior halved"; that family thins the band at every rung (its cell count
does not grow), so the local width ratio between rungs is not 2 everywhere; the nested family
here is the refinement that the packet's stated purpose ("so the order formulas hold")
requires, and its finest rung is the packet's mesh.

Only the master grid's faces change: the native ``SphericalGrid`` is built by the native
``_from_faces`` (exact spherical volumes and areas), and every invariant the kernel checks
(the dry region's exact metrics, the cut geometry's master grid, the ledgers' area terms) is
checked on it by ``validate_master_grid`` before it is used. ``physically_qualifying`` is
false.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import math

RECORD = "docs/GT_PS2_PAPER1_SURFACE_GRADED_LADDER_2026-10-22.md"
PACKET = "docs/GT_PS2_PAPER1_SURFACE_LAYER_MESH_PACKET_2026-10-18.md"
RULING = "GT_PS2_OWNER_RULING_UNIVERSAL_KERNEL_ENVELOPE_2026-10-19 item 7"


class SurfaceGradedMeshRefused(ValueError):
    """A declared face set that is not admissible (refused, never repaired)."""


@dataclass(frozen=True)
class SurfaceGradedFaceSet:
    """The declared surface-graded face set: m interior cells, n band cells at ratio q."""

    interior_cells: int
    band_cells: int
    band_ratio: float

    def __post_init__(self):
        for name in ("interior_cells", "band_cells"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise SurfaceGradedMeshRefused(f"{name} must be a positive integer, not {value!r}")
        q = self.band_ratio
        if isinstance(q, bool) or not isinstance(q, (int, float)):
            raise SurfaceGradedMeshRefused(f"band_ratio must be a real number, not {q!r}")
        if not (math.isfinite(q) and q > 1.0):
            raise SurfaceGradedMeshRefused(
                f"band_ratio must be finite and above 1 (a refined band), not {q!r}"
            )

    @property
    def cells(self):
        return self.interior_cells + self.band_cells

    def amplitude_m(self, radius_m):
        """A of the band map: the depth below the surface is A (q^j - 1) at band face j."""
        q, n = float(self.band_ratio), self.band_cells
        growth = q**n
        return radius_m / (self.interior_cells * growth * math.log(q) + growth - 1.0)

    def band_inner_radius_m(self, radius_m):
        q, n = float(self.band_ratio), self.band_cells
        return radius_m - self.amplitude_m(radius_m) * (q**n - 1.0)

    def faces(self, radius_m):
        """The N + 1 face radii, from the centre (exactly 0) to the surface (exactly R)."""
        if not (isinstance(radius_m, float) and math.isfinite(radius_m) and radius_m > 0.0):
            raise SurfaceGradedMeshRefused(f"radius must be a positive finite float: {radius_m!r}")
        q, m, n = float(self.band_ratio), self.interior_cells, self.band_cells
        amplitude = self.amplitude_m(radius_m)
        inner = self.band_inner_radius_m(radius_m)
        if not (0.0 < inner < radius_m):
            raise SurfaceGradedMeshRefused("the band does not leave a positive uniform interior")
        interior = [inner * i / m for i in range(m)]
        band = [radius_m - amplitude * (q**j - 1.0) for j in range(n, 0, -1)]
        faces = (0.0, *interior[1:], *band, radius_m)
        if len(faces) != m + n + 1:
            raise SurfaceGradedMeshRefused("face count does not match the declared cells")
        if any(right <= left for left, right in zip(faces, faces[1:])):
            raise SurfaceGradedMeshRefused("faces are not strictly increasing")
        return faces

    def refined(self):
        """The next finer rung of the ratio-2 family: (2 m, 2 n, q^(1/2))."""
        return SurfaceGradedFaceSet(
            2 * self.interior_cells, 2 * self.band_cells, math.sqrt(float(self.band_ratio))
        )

    def coarsened(self):
        """The next coarser rung of the ratio-2 family: (m / 2, n / 2, q^2); refused when odd."""
        if self.interior_cells % 2 or self.band_cells % 2:
            raise SurfaceGradedMeshRefused("an odd cell count has no ratio-2 coarser rung")
        return SurfaceGradedFaceSet(
            self.interior_cells // 2, self.band_cells // 2, float(self.band_ratio) ** 2
        )

    def describe(self, radius_m):
        faces = self.faces(radius_m)
        widths = [right - left for left, right in zip(faces, faces[1:])]
        m = self.interior_cells
        return {
            "kind": "DECLARED_ENGINEERING_ASSUMPTION",
            "name": "surface_graded_radial_master_grid",
            "record": RECORD,
            "packet": PACKET,
            "authority": RULING + " (the graded ladder) and the packet's section 4 plan "
            "(Class A numerics); the knob off reproduces the uniform grid",
            "interior_cells": m,
            "band_cells": self.band_cells,
            "band_ratio": float(self.band_ratio),
            "cells": self.cells,
            "radius_m": radius_m,
            "map": "uniform interior r = s xi; band depth A (q^((m+n)(1-xi)) - 1); C1 at "
            "xi_b = m/(m+n)",
            "band_amplitude_m": self.amplitude_m(radius_m),
            "band_inner_radius_m": faces[m],
            "band_thickness_m": radius_m - faces[m],
            "first_cell_width_m": widths[-1],
            "last_band_cell_width_m": widths[m],
            "interior_cell_width_m": widths[0],
            "junction_width_ratio": widths[m - 1] / widths[m],
            "largest_neighbour_width_ratio": max(
                max(a, b) / min(a, b) for a, b in zip(widths, widths[1:])
            ),
            "faces_m": list(faces),
            "refinement_family": "(m, n, q) -> (2m, 2n, sqrt(q)): nested, every cell split in two",
            "physically_qualifying": False,
        }


def parse(interior_cells, band_cells, band_ratio):
    return SurfaceGradedFaceSet(int(interior_cells), int(band_cells), float(band_ratio))


def native_grid(face_set, radius_m):
    """The native SphericalGrid on the declared faces, checked by the kernel's invariants."""
    from dtdc_simulator.core2.particle import grid as spherical_grid

    grid = spherical_grid._from_faces(face_set.faces(radius_m), radius_m)
    validate_master_grid(grid)
    return grid


def validate_master_grid(grid):
    """Every mesh invariant the kernel checks, applied to a master grid (refuses typed).

    The native validators (the dry region's exact spherical metrics at 2e-15, the cut
    geometry's master grid at 128 ulp, its topology from r = 0 to R) and the ledgers' area
    terms: the outer face area is exactly 4 pi R^2 as the native constructor computes it, the
    face areas telescope the volume (sum V = 4/3 pi R^3 to rounding), and no cell is empty.
    """
    from dtdc_simulator.core2.particle import cut_geometry, dry_region

    try:
        dry_region._validate_grid_metrics(grid)
        cut_geometry._validate_master_grid(grid)
    except Exception as error:  # the native types differ; the refusal is re-typed here
        raise SurfaceGradedMeshRefused(f"native grid invariant refused: {error}") from error
    if grid.areas[-1] != 4.0 * math.pi * grid.R * grid.R:
        raise SurfaceGradedMeshRefused("the outer face area is not the native 4 pi R^2")
    total = math.fsum(grid.volumes)
    sphere = 4.0 / 3.0 * math.pi * grid.R**3
    if abs(total - sphere) > 64.0 * math.ulp(sphere) * grid.n:
        raise SurfaceGradedMeshRefused("cell volumes do not sum to the sphere")
    if any(not (v > 0.0) for v in grid.volumes):
        raise SurfaceGradedMeshRefused("a cell has no positive volume")
    if any(
        not (left < centre < right)
        for left, centre, right in zip(grid.faces, grid.centers, grid.faces[1:])
    ):
        raise SurfaceGradedMeshRefused("a cell centre lies outside its faces")
    return True


_INSTALLED = None


@contextmanager
def installed(face_set):
    """Make the birth's grid construction use the declared face set (None: no change)."""
    global _INSTALLED
    if face_set is not None and not isinstance(face_set, SurfaceGradedFaceSet):
        raise SurfaceGradedMeshRefused("installed() takes a SurfaceGradedFaceSet or None")
    previous = _INSTALLED
    _INSTALLED = face_set
    try:
        yield face_set
    finally:
        _INSTALLED = previous


def installed_face_set():
    return _INSTALLED


def master_grid_or_none(cells, radius_m):
    """The graded master grid when a face set is installed, else None (the uniform path).

    The installed face set must declare exactly the requested cell count: the runner's
    ``--cells`` is the graded total, so the feed, the activation and every per-cell budget
    see the same N as the grid.
    """
    face_set = _INSTALLED
    if face_set is None:
        return None
    if face_set.cells != cells:
        raise SurfaceGradedMeshRefused(
            f"the installed face set has {face_set.cells} cells, the birth asked for {cells}"
        )
    return native_grid(face_set, float(radius_m))
