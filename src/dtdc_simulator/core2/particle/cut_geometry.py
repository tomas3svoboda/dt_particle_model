"""Fixed-master-grid cut geometry for the Gate-1g numerical foundation.

This module partitions stationary spherical material cells with one sharp
front.  It is intentionally geometry only: it contains no transport
coefficient, storage closure, Rankine--Hugoniot solve, or constitutive phase
fraction.  In particular, ``coupled_pore.gas_accessible_fraction`` is a
constitutive pore-accessibility input and **must never** be used as the
geometric wet/dry fraction of a material cell.  The only geometric measures
exposed here are exact spherical intersection volumes.

The front uses ``z=(s/R)^3`` together with its canonical radius ``s``.  A
front within a few floating-point ULPs of a master face is snapped to that
face.  This is an identity canonicalisation, not a physical epsilon layer:
at a canonical face or endpoint there is no cut cell and all zero-volume
pieces are exactly zero.

Every result is numerical and reports ``physically_qualifying=False``.  A
transport/RH integrator must qualify the physics separately.
"""

from __future__ import annotations

import bisect
import math
import operator
from dataclasses import dataclass, field
from enum import Enum

from dtdc_simulator.core2.particle.grid import SphericalGrid


DEFAULT_FACE_SNAP_ULPS = 8


class CutGeometryError(ValueError):
    """The requested fixed-grid/front geometry is numerically inadmissible."""


class PrimaryDrainageError(CutGeometryError):
    """A proposed front step violates monotone primary drainage."""


class CellRegion(str, Enum):
    """Geometric phase occupancy of one stationary material cell."""

    WET = "wet"
    DRY = "dry"
    CUT = "cut"


class FrontStepClassification(str, Enum):
    """Numerical classification of one proposed primary-drainage segment."""

    NO_MOTION = "no_motion"
    ACTIVATION_BIRTH = "activation_birth"
    WITHIN_CELL = "within_cell"
    FACE_ARRIVAL = "face_arrival"
    FACE_CROSSING = "face_crossing"
    MULTI_FACE_REFINEMENT = "multi_face_refinement"
    EXTINCTION = "extinction"


@dataclass(frozen=True)
class CanonicalFront:
    """One front represented simultaneously by exact endpoint-aware ``z`` and ``s``."""

    particle_radius_m: float
    z: float
    radius_m: float
    master_face_index: int | None
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (self.particle_radius_m, self.z, self.radius_m)
        if not all(math.isfinite(value) for value in values):
            raise CutGeometryError("front coordinates must be finite")
        if self.particle_radius_m <= 0.0:
            raise CutGeometryError("particle radius must be positive")
        if not 0.0 <= self.z <= 1.0:
            raise CutGeometryError("front volume coordinate z must lie in [0, 1]")
        if not 0.0 <= self.radius_m <= self.particle_radius_m:
            raise CutGeometryError("front radius s must lie in [0, R]")
        if self.z == 0.0 and self.radius_m != 0.0:
            raise CutGeometryError("z=0 must have the exact extinction radius s=0")
        if self.z == 1.0 and self.radius_m != self.particle_radius_m:
            raise CutGeometryError("z=1 must have the exact activation radius s=R")
        if self.radius_m == 0.0 and self.z != 0.0:
            raise CutGeometryError("s=0 must have the exact extinction coordinate z=0")
        if self.radius_m == self.particle_radius_m and self.z != 1.0:
            raise CutGeometryError("s=R must have the exact activation coordinate z=1")
        if self.master_face_index is not None and self.master_face_index < 0:
            raise CutGeometryError("master face index must be non-negative")

    @property
    def regime(self) -> str:
        if self.z == 1.0:
            return "fully_wet"
        if self.z == 0.0:
            return "fully_dry"
        return "partial"

    @property
    def at_master_face(self) -> bool:
        return self.master_face_index is not None


@dataclass(frozen=True)
class MaterialCellPartition:
    """Exact spherical wet/dry intersection volumes in one master cell."""

    index: int
    inner_radius_m: float
    outer_radius_m: float
    wet_volume_m3: float
    dry_volume_m3: float
    region: CellRegion
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.inner_radius_m,
            self.outer_radius_m,
            self.wet_volume_m3,
            self.dry_volume_m3,
        )
        if self.index < 0:
            raise CutGeometryError("material-cell index must be non-negative")
        if not all(math.isfinite(value) for value in values):
            raise CutGeometryError("material-cell geometry must be finite")
        if self.inner_radius_m < 0.0 or self.outer_radius_m <= self.inner_radius_m:
            raise CutGeometryError("material-cell radii must be positive-width and ordered")
        if self.wet_volume_m3 < 0.0 or self.dry_volume_m3 < 0.0:
            raise CutGeometryError("material-cell subvolumes must be non-negative")
        if self.region is CellRegion.WET and self.dry_volume_m3 != 0.0:
            raise CutGeometryError("a wet cell must have exactly zero dry volume")
        if self.region is CellRegion.DRY and self.wet_volume_m3 != 0.0:
            raise CutGeometryError("a dry cell must have exactly zero wet volume")
        if self.region is CellRegion.CUT and not (
            self.wet_volume_m3 > 0.0 and self.dry_volume_m3 > 0.0
        ):
            raise CutGeometryError("a cut cell must have two strictly positive pieces")

    @property
    def volume_m3(self) -> float:
        return math.fsum((self.wet_volume_m3, self.dry_volume_m3))


@dataclass(frozen=True)
class CutGeometry:
    """A sharp front intersecting one fixed, full-sphere master grid."""

    master_grid: SphericalGrid
    front: CanonicalFront
    cells: tuple[MaterialCellPartition, ...]
    cut_cell_index: int | None
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if len(self.cells) != self.master_grid.n:
            raise CutGeometryError("one partition is required for every material cell")
        if self.front.particle_radius_m != self.master_grid.R:
            raise CutGeometryError("front and master-grid particle radii must be identical")
        cut_indices = tuple(cell.index for cell in self.cells if cell.region is CellRegion.CUT)
        expected = () if self.cut_cell_index is None else (self.cut_cell_index,)
        if cut_indices != expected:
            raise CutGeometryError("there must be exactly one declared cut cell or none")
        if self.front.at_master_face and self.cut_cell_index is not None:
            raise CutGeometryError("a front on a master face cannot leave a cut cell")

    @property
    def wet_volume_m3(self) -> float:
        return math.fsum(cell.wet_volume_m3 for cell in self.cells)

    @property
    def dry_volume_m3(self) -> float:
        return math.fsum(cell.dry_volume_m3 for cell in self.cells)

    @property
    def particle_volume_m3(self) -> float:
        return self.master_grid.total_volume

    @property
    def spherical_wet_volume_m3(self) -> float:
        return _sphere_volume(self.front.radius_m)

    @property
    def spherical_dry_volume_m3(self) -> float:
        return math.fsum((self.particle_volume_m3, -self.spherical_wet_volume_m3))

    @property
    def volume_closure_residual_m3(self) -> float:
        return math.fsum((self.wet_volume_m3, self.dry_volume_m3, -self.particle_volume_m3))

    @property
    def max_cell_closure_residual_m3(self) -> float:
        return max(
            abs(
                math.fsum(
                    (
                        cell.wet_volume_m3,
                        cell.dry_volume_m3,
                        -self.master_grid.volumes[cell.index],
                    )
                )
            )
            for cell in self.cells
        )


@dataclass(frozen=True)
class SweptCutGeometry:
    """Exact finite-step spherical sweep and cellwise geometric conservation data."""

    previous: CutGeometry
    current: CutGeometry
    dt_s: float
    wet_to_dry_volume_m3: float
    cell_wet_to_dry_volumes_m3: tuple[float, ...]
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.dt_s) or self.dt_s <= 0.0:
            raise CutGeometryError("geometry time step must be positive and finite")
        if not math.isfinite(self.wet_to_dry_volume_m3) or self.wet_to_dry_volume_m3 < 0.0:
            raise CutGeometryError("swept wet-to-dry volume must be finite and non-negative")
        if len(self.cell_wet_to_dry_volumes_m3) != self.previous.master_grid.n:
            raise CutGeometryError("one swept volume is required for every material cell")
        if not all(
            math.isfinite(value) and value >= 0.0
            for value in self.cell_wet_to_dry_volumes_m3
        ):
            raise CutGeometryError("cell swept volumes must be finite and non-negative")

    @property
    def signed_wet_volume_change_m3(self) -> float:
        return -self.wet_to_dry_volume_m3

    @property
    def signed_dry_volume_change_m3(self) -> float:
        return self.wet_to_dry_volume_m3

    @property
    def interface_swept_volume_rate_m3_s(self) -> float:
        """Signed ``dV_w/dt``; negative for a receding primary-drainage front."""

        return -self.wet_to_dry_volume_m3 / self.dt_s

    @property
    def cell_wet_volume_changes_m3(self) -> tuple[float, ...]:
        return tuple(-value for value in self.cell_wet_to_dry_volumes_m3)

    @property
    def cell_dry_volume_changes_m3(self) -> tuple[float, ...]:
        return self.cell_wet_to_dry_volumes_m3

    @property
    def cell_sweep_closure_residual_m3(self) -> float:
        return math.fsum((*self.cell_wet_to_dry_volumes_m3, -self.wet_to_dry_volume_m3))

    @property
    def max_cell_gcl_residual_m3(self) -> float:
        residuals: list[float] = []
        for old, new, swept in zip(
            self.previous.cells,
            self.current.cells,
            self.cell_wet_to_dry_volumes_m3,
        ):
            residuals.append(
                abs(math.fsum((new.wet_volume_m3, -old.wet_volume_m3, swept)))
            )
            residuals.append(
                abs(math.fsum((new.dry_volume_m3, -old.dry_volume_m3, -swept)))
            )
        return max(residuals, default=0.0)


@dataclass(frozen=True)
class RemainingStepMetadata:
    """Unresolved part of a proposal after landing on its first master-face event.

    The time fraction is linear interpolation in the proposed volume coordinate
    ``z``.  It is numerical proposal metadata only; it is not a physical claim
    that the true front velocity remains constant.  The coupled nonlinear
    problem must be solved again from ``restart``.
    """

    restart: CutGeometry
    requested_target: CutGeometry
    duration_s: float
    proposal_fraction: float
    unresolved_face_indices: tuple[int, ...]
    must_resolve_again: bool = field(default=True, init=False)
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not math.isfinite(self.duration_s) or self.duration_s <= 0.0:
            raise CutGeometryError("remaining event-step duration must be positive and finite")
        if not math.isfinite(self.proposal_fraction) or not 0.0 < self.proposal_fraction < 1.0:
            raise CutGeometryError("remaining proposal fraction must lie strictly in (0, 1)")


@dataclass(frozen=True)
class FrontStepPlan:
    """Canonical first-event plan for one monotone proposed front segment.

    For a proposal that crosses a master face, ``accepted_dt_s`` is only a
    volume-coordinate interpolation used to seed an event solve.  It is not a
    converged physical event time and the state at ``accepted`` must not be
    committed until the coupled ALE/RH residual has been solved again with
    that face as its endpoint.
    """

    previous: CutGeometry
    requested: CutGeometry
    accepted: CutGeometry
    classification: FrontStepClassification
    requested_dt_s: float
    accepted_dt_s: float
    accepted_proposal_fraction: float
    encountered_face_indices: tuple[int, ...]
    first_event_face_index: int | None
    remaining: RemainingStepMetadata | None
    starts_at_activation_endpoint: bool
    requests_extinction_endpoint: bool
    requires_refinement: bool
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def remaining_dt_s(self) -> float:
        return 0.0 if self.remaining is None else self.remaining.duration_s

    @property
    def has_remaining_step(self) -> bool:
        return self.remaining is not None

    @property
    def must_resolve_accepted_event(self) -> bool:
        """Whether the interpolated face arrival still needs a coupled solve."""

        return self.classification in (
            FrontStepClassification.FACE_CROSSING,
            FrontStepClassification.MULTI_FACE_REFINEMENT,
        )


def partition_master_grid(
    master_grid: SphericalGrid,
    *,
    z: float | None = None,
    radius_m: float | None = None,
    face_snap_ulps: int = DEFAULT_FACE_SNAP_ULPS,
) -> CutGeometry:
    """Intersect a fixed full-sphere master grid with one canonical sharp front.

    Exactly one of ``z`` or ``radius_m`` must be supplied.  No porosity,
    accessibility, saturation, or other constitutive fraction is accepted by
    this geometry API.
    """

    _validate_master_grid(master_grid)
    snap_ulps = _validate_face_snap_ulps(face_snap_ulps)
    front = _canonical_front(master_grid, z=z, radius_m=radius_m, snap_ulps=snap_ulps)

    if front.master_face_index is None:
        cut_cell_index = bisect.bisect_right(master_grid.faces, front.radius_m) - 1
        if not 0 <= cut_cell_index < master_grid.n:
            raise CutGeometryError("interior front did not resolve to one master cell")
    else:
        cut_cell_index = None

    cells: list[MaterialCellPartition] = []
    for index, (inner, outer, volume) in enumerate(
        zip(master_grid.faces, master_grid.faces[1:], master_grid.volumes)
    ):
        if front.master_face_index is not None:
            if index < front.master_face_index:
                wet, dry, region = volume, 0.0, CellRegion.WET
            else:
                wet, dry, region = 0.0, volume, CellRegion.DRY
        elif index < cut_cell_index:
            wet, dry, region = volume, 0.0, CellRegion.WET
        elif index > cut_cell_index:
            wet, dry, region = 0.0, volume, CellRegion.DRY
        else:
            wet, dry = _cut_cell_volumes(inner, front.radius_m, outer, volume)
            region = CellRegion.CUT
        cells.append(
            MaterialCellPartition(
                index=index,
                inner_radius_m=inner,
                outer_radius_m=outer,
                wet_volume_m3=wet,
                dry_volume_m3=dry,
                region=region,
            )
        )

    result = CutGeometry(
        master_grid=master_grid,
        front=front,
        cells=tuple(cells),
        cut_cell_index=cut_cell_index,
    )
    _validate_partition_closure(result)
    return result


def swept_cut_geometry(
    previous: CutGeometry,
    current: CutGeometry,
    dt_s: float,
) -> SweptCutGeometry:
    """Return the exact spherical wet-to-dry sweep for one primary-drainage step."""

    _validate_geometry_pair(previous, current)
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise CutGeometryError("geometry time step must be positive and finite")
    if current.front.radius_m > previous.front.radius_m:
        raise PrimaryDrainageError("primary-drainage front cannot advance outward")

    swept = _shell_volume(current.front.radius_m, previous.front.radius_m)
    cell_sweeps = tuple(
        _interval_shell_volume(
            max(cell.inner_radius_m, current.front.radius_m),
            min(cell.outer_radius_m, previous.front.radius_m),
        )
        for cell in previous.cells
    )
    result = SweptCutGeometry(
        previous=previous,
        current=current,
        dt_s=dt_s,
        wet_to_dry_volume_m3=swept,
        cell_wet_to_dry_volumes_m3=cell_sweeps,
    )
    scale = max(previous.master_grid.total_volume, math.ulp(previous.master_grid.total_volume))
    if abs(result.cell_sweep_closure_residual_m3) > 64.0 * math.ulp(scale):
        raise CutGeometryError("cell sweeps do not close to the exact spherical front sweep")
    if result.max_cell_gcl_residual_m3 > 64.0 * math.ulp(scale):
        raise CutGeometryError("cell partitions violate the geometric conservation law")
    return result


def plan_primary_drainage_step(
    previous: CutGeometry,
    requested: CutGeometry,
    dt_s: float,
) -> FrontStepPlan:
    """Canonicalise a proposal at its first encountered master-face event.

    A proposal that crosses a face is accepted only to the first face.  The
    rest is returned as :class:`RemainingStepMetadata` and must be solved
    again.  Crossing more than one face in the original proposal is explicitly
    classified as requiring refinement.
    """

    _validate_geometry_pair(previous, requested)
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise CutGeometryError("proposed step duration must be positive and finite")
    s_previous = previous.front.radius_m
    s_requested = requested.front.radius_m
    if s_requested > s_previous:
        raise PrimaryDrainageError("primary-drainage front cannot advance outward")

    encountered = tuple(
        index
        for index in range(previous.master_grid.n, -1, -1)
        if s_requested <= previous.master_grid.faces[index] < s_previous
    )
    starts_at_activation = previous.front.z == 1.0
    requests_extinction = requested.front.z == 0.0

    if s_requested == s_previous:
        return FrontStepPlan(
            previous=previous,
            requested=requested,
            accepted=requested,
            classification=FrontStepClassification.NO_MOTION,
            requested_dt_s=dt_s,
            accepted_dt_s=dt_s,
            accepted_proposal_fraction=1.0,
            encountered_face_indices=(),
            first_event_face_index=None,
            remaining=None,
            starts_at_activation_endpoint=starts_at_activation,
            requests_extinction_endpoint=requests_extinction,
            requires_refinement=False,
        )

    if not encountered:
        classification = (
            FrontStepClassification.ACTIVATION_BIRTH
            if starts_at_activation
            else FrontStepClassification.WITHIN_CELL
        )
        return FrontStepPlan(
            previous=previous,
            requested=requested,
            accepted=requested,
            classification=classification,
            requested_dt_s=dt_s,
            accepted_dt_s=dt_s,
            accepted_proposal_fraction=1.0,
            encountered_face_indices=(),
            first_event_face_index=None,
            remaining=None,
            starts_at_activation_endpoint=starts_at_activation,
            requests_extinction_endpoint=requests_extinction,
            requires_refinement=False,
        )

    first_face = encountered[0]
    event = partition_master_grid(
        previous.master_grid,
        radius_m=previous.master_grid.faces[first_face],
    )
    if event.front.radius_m == s_requested:
        classification = (
            FrontStepClassification.EXTINCTION
            if first_face == 0
            else FrontStepClassification.FACE_ARRIVAL
        )
        return FrontStepPlan(
            previous=previous,
            requested=requested,
            accepted=requested,
            classification=classification,
            requested_dt_s=dt_s,
            accepted_dt_s=dt_s,
            accepted_proposal_fraction=1.0,
            encountered_face_indices=encountered,
            first_event_face_index=first_face,
            remaining=None,
            starts_at_activation_endpoint=starts_at_activation,
            requests_extinction_endpoint=requests_extinction,
            requires_refinement=False,
        )

    total_delta_z = previous.front.z - requested.front.z
    accepted_delta_z = previous.front.z - event.front.z
    accepted_fraction = accepted_delta_z / total_delta_z
    if not 0.0 < accepted_fraction < 1.0:
        raise CutGeometryError("first-face proposal fraction must lie strictly in (0, 1)")
    accepted_dt = dt_s * accepted_fraction
    remaining_dt = dt_s - accepted_dt
    remaining_fraction = 1.0 - accepted_fraction
    if accepted_dt <= 0.0 or remaining_dt <= 0.0:
        raise CutGeometryError("face event has no representable positive substep")
    multi_face = len(encountered) > 1
    classification = (
        FrontStepClassification.MULTI_FACE_REFINEMENT
        if multi_face
        else FrontStepClassification.FACE_CROSSING
    )
    remaining = RemainingStepMetadata(
        restart=event,
        requested_target=requested,
        duration_s=remaining_dt,
        proposal_fraction=remaining_fraction,
        unresolved_face_indices=encountered[1:],
    )
    return FrontStepPlan(
        previous=previous,
        requested=requested,
        accepted=event,
        classification=classification,
        requested_dt_s=dt_s,
        accepted_dt_s=accepted_dt,
        accepted_proposal_fraction=accepted_fraction,
        encountered_face_indices=encountered,
        first_event_face_index=first_face,
        remaining=remaining,
        starts_at_activation_endpoint=starts_at_activation,
        requests_extinction_endpoint=requests_extinction,
        requires_refinement=multi_face,
    )


def _canonical_front(
    master_grid: SphericalGrid,
    *,
    z: float | None,
    radius_m: float | None,
    snap_ulps: int,
) -> CanonicalFront:
    if (z is None) == (radius_m is None):
        raise CutGeometryError("supply exactly one of z or radius_m")
    if z is not None:
        if isinstance(z, bool) or not isinstance(z, (int, float)) or not math.isfinite(z):
            raise CutGeometryError("front volume coordinate z must be a finite real number")
        if not 0.0 <= z <= 1.0:
            raise CutGeometryError("front volume coordinate z must lie in [0, 1]")
        z_value = float(z)
        if z_value == 0.0:
            radius = 0.0
        elif z_value == 1.0:
            radius = master_grid.R
        else:
            radius = master_grid.R * z_value ** (1.0 / 3.0)
    else:
        if (
            isinstance(radius_m, bool)
            or not isinstance(radius_m, (int, float))
            or not math.isfinite(radius_m)
        ):
            raise CutGeometryError("front radius s must be a finite real number")
        if not 0.0 <= radius_m <= master_grid.R:
            raise CutGeometryError("front radius s must lie in [0, R]")
        radius = float(radius_m)
        z_value = _z_from_radius(radius, master_grid.R)

    face_index = _near_face_index(master_grid.faces, radius, snap_ulps)
    if face_index is not None:
        radius = master_grid.faces[face_index]
        z_value = _z_from_radius(radius, master_grid.R)
    return CanonicalFront(
        particle_radius_m=master_grid.R,
        z=z_value,
        radius_m=radius,
        master_face_index=face_index,
    )


def _near_face_index(faces: tuple[float, ...], radius_m: float, snap_ulps: int) -> int | None:
    insertion = bisect.bisect_left(faces, radius_m)
    candidates = {
        index for index in (insertion - 1, insertion) if 0 <= index < len(faces)
    }
    near: list[tuple[float, int]] = []
    for index in candidates:
        face = faces[index]
        distance = abs(radius_m - face)
        tolerance = snap_ulps * max(math.ulp(radius_m), math.ulp(face))
        if distance <= tolerance:
            near.append((distance, index))
    if not near:
        return None
    return min(near)[1]


def _cut_cell_volumes(
    inner: float,
    front: float,
    outer: float,
    master_volume: float,
) -> tuple[float, float]:
    if not inner < front < outer:
        raise CutGeometryError("a cut-cell front must lie strictly between its master faces")
    # Evaluate the smaller radial piece directly with the cancellation-free
    # shell formula, and obtain the other as its exact material-cell complement.
    if front - inner <= outer - front:
        wet = _shell_volume(inner, front)
        dry = master_volume - wet
    else:
        dry = _shell_volume(front, outer)
        wet = master_volume - dry
    if not 0.0 < wet < master_volume or not 0.0 < dry < master_volume:
        raise CutGeometryError("cut-cell subvolumes lost representable positive measure")
    return wet, dry


def _validate_master_grid(master_grid: SphericalGrid) -> None:
    if not isinstance(master_grid, SphericalGrid):
        raise CutGeometryError("master_grid must be a SphericalGrid")
    if not math.isfinite(master_grid.R) or master_grid.R <= 0.0:
        raise CutGeometryError("master-grid radius must be positive and finite")
    if master_grid.n < 1:
        raise CutGeometryError("fixed material master grid needs at least one cell")
    if (
        len(master_grid.faces) != master_grid.n + 1
        or len(master_grid.volumes) != master_grid.n
        or len(master_grid.areas) != master_grid.n + 1
    ):
        raise CutGeometryError("master-grid metric arrays are not aligned")
    if master_grid.inner_radius != 0.0:
        raise CutGeometryError("fixed material master grid must cover the full sphere from r=0")
    if master_grid.faces[-1] != master_grid.R:
        raise CutGeometryError("master-grid outer face must equal R exactly")
    if not all(math.isfinite(face) for face in master_grid.faces):
        raise CutGeometryError("master-grid faces must be finite")
    if master_grid.faces[0] < 0.0 or any(
        right <= left for left, right in zip(master_grid.faces, master_grid.faces[1:])
    ):
        raise CutGeometryError("master-grid faces must be strictly increasing")
    if any(
        not math.isfinite(volume) or volume <= 0.0 for volume in master_grid.volumes
    ):
        raise CutGeometryError("master material cells must have positive finite volume")
    if not all(math.isfinite(center) for center in master_grid.centers):
        raise CutGeometryError("master-grid centers must be finite")
    if not all(math.isfinite(area) and area >= 0.0 for area in master_grid.areas):
        raise CutGeometryError("master-grid face areas must be finite and non-negative")

    def metric_matches(actual: float, expected: float) -> bool:
        scale = max(abs(actual), abs(expected))
        if scale == 0.0:
            return actual == expected
        return abs(actual - expected) <= 128.0 * math.ulp(scale)

    for index, (inner, outer, center, volume) in enumerate(
        zip(
            master_grid.faces,
            master_grid.faces[1:],
            master_grid.centers,
            master_grid.volumes,
        )
    ):
        if not metric_matches(center, 0.5 * (inner + outer)):
            raise CutGeometryError(
                f"master-grid center {index} is inconsistent with its faces"
            )
        if not metric_matches(volume, _shell_volume(inner, outer)):
            raise CutGeometryError(
                f"master-grid volume {index} is inconsistent with its spherical faces"
            )
    for index, (face, area) in enumerate(zip(master_grid.faces, master_grid.areas)):
        if not metric_matches(area, 4.0 * math.pi * face * face):
            raise CutGeometryError(
                f"master-grid area {index} is inconsistent with its spherical face"
            )


def _validate_face_snap_ulps(face_snap_ulps: int) -> int:
    try:
        value = operator.index(face_snap_ulps)
    except TypeError as exc:
        raise CutGeometryError("face_snap_ulps must be an integer") from exc
    if isinstance(face_snap_ulps, bool) or not 0 <= value <= 1024:
        raise CutGeometryError("face_snap_ulps must be an integer in [0, 1024]")
    return value


def _validate_partition_closure(geometry: CutGeometry) -> None:
    scale = max(geometry.particle_volume_m3, math.ulp(geometry.particle_volume_m3))
    tolerance = 64.0 * math.ulp(scale)
    if abs(geometry.volume_closure_residual_m3) > tolerance:
        raise CutGeometryError("wet/dry material-cell volumes do not close to the sphere")
    if geometry.max_cell_closure_residual_m3 > tolerance:
        raise CutGeometryError("wet/dry pieces do not close to their material cells")
    if abs(geometry.wet_volume_m3 - geometry.spherical_wet_volume_m3) > tolerance:
        raise CutGeometryError("cellwise wet volume does not match the spherical front volume")


def _validate_geometry_pair(previous: CutGeometry, current: CutGeometry) -> None:
    if not isinstance(previous, CutGeometry) or not isinstance(current, CutGeometry):
        raise CutGeometryError("front steps require two CutGeometry states")
    if previous.master_grid != current.master_grid:
        raise CutGeometryError("front step must retain one identical fixed material master grid")


def _sphere_volume(radius_m: float) -> float:
    return 4.0 / 3.0 * math.pi * radius_m**3


def _shell_volume(inner_radius_m: float, outer_radius_m: float) -> float:
    if outer_radius_m < inner_radius_m:
        raise CutGeometryError("spherical shell radii must be ordered")
    if outer_radius_m == inner_radius_m:
        return 0.0
    return (
        4.0
        / 3.0
        * math.pi
        * (outer_radius_m - inner_radius_m)
        * math.fsum(
            (
                outer_radius_m * outer_radius_m,
                outer_radius_m * inner_radius_m,
                inner_radius_m * inner_radius_m,
            )
        )
    )


def _interval_shell_volume(inner_radius_m: float, outer_radius_m: float) -> float:
    if outer_radius_m <= inner_radius_m:
        return 0.0
    return _shell_volume(inner_radius_m, outer_radius_m)


def _z_from_radius(radius_m: float, particle_radius_m: float) -> float:
    if radius_m == 0.0:
        return 0.0
    if radius_m == particle_radius_m:
        return 1.0
    return (radius_m / particle_radius_m) ** 3
