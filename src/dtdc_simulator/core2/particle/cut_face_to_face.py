r"""Two-ended exact master-face event residual for primary drainage.

The immutable before state lies on interior master face ``K > 1``.  The event
endpoint lies on the immediately adjacent inner face ``k = K - 1``.  Material
cell ``k`` is therefore fully wet at the start and fully dry at the endpoint.
Neither its zero-volume start-dry piece nor its zero-volume endpoint-wet piece
is represented by a primitive or balance row.

The rank reduction is the exact full-cell analogue of
:mod:`cut_face_event`: the first endpoint dry-cell balance compares its full
dry inventory with the complete old wet inventory, and uses the fixed wet
inner-face rate plus the endpoint dry outer-face rate.  Three independent
ALE Rankine--Hugoniot rows retain the endpoint interface jump.  The swept
volume rate is the exact spherical volume of material cell ``k`` divided by
the unknown event duration.

For ``n_w=k`` endpoint wet cells, ``n_d=N-k`` endpoint dry cells, and
``n_f=n_d+1`` dry faces, the unknowns are

``2*n_w wet + 2*n_d dry + n_f Stefan + duration + interface T``

and the residuals are

``2*n_w wet + 3*n_d dry + 3 RH``.

Both counts equal ``2*n_w + 3*n_d + 3``.  No overlap duration, epsilon shell,
new physics equation, or endpoint pinning condition occurs in this residual.
All coefficients remain nonqualifying numerical oracles.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from dtdc_simulator.core2.particle import conditioning
from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_face_event as arrival
from dtdc_simulator.core2.particle import cut_face_event_integrator as fi
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.props import hexane as hx


class FaceToFaceTopologyError(ValueError):
    """Proposal is not one adjacent, inward, exact-face event."""


class FaceToFaceStepError(RuntimeError):
    """Rejected residual carrying the exact immutable face rollback state."""

    def __init__(
        self,
        message: str,
        rollback_state: fi.FaceArrivalCommittedState,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state


@dataclass(frozen=True)
class FaceToFaceLayout:
    """Square rank after eliminating both endpoint zero-volume pieces."""

    wet_cell_indices: tuple[int, ...]
    dry_cell_indices: tuple[int, ...]
    departure_face_index: int
    arrival_face_index: int
    transformed_cell_index: int
    master_cell_count: int
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        start = self.departure_face_index
        end = self.arrival_face_index
        if not 1 < start < self.master_cell_count:
            raise ValueError("face-to-face start requires 1 < K < N")
        if end != start - 1 or self.transformed_cell_index != end:
            raise ValueError("face-to-face endpoint must be the adjacent inner face")
        if self.wet_cell_indices != tuple(range(end)):
            raise ValueError("endpoint wet cells must be exactly 0..k-1")
        if self.dry_cell_indices != tuple(range(end, self.master_cell_count)):
            raise ValueError("endpoint dry cells must be exactly k..N-1")
        if not self.is_square:
            raise ValueError("face-to-face rank formula is not square")

    @property
    def wet_piece_count(self) -> int:
        return len(self.wet_cell_indices)

    @property
    def dry_piece_count(self) -> int:
        return len(self.dry_cell_indices)

    @property
    def dry_face_count(self) -> int:
        return self.dry_piece_count + 1

    @property
    def unknown_blocks(self) -> tuple[tuple[str, int], ...]:
        return (
            ("wet_temperature", self.wet_piece_count),
            ("wet_retained_water", self.wet_piece_count),
            ("dry_temperature", self.dry_piece_count),
            ("dry_y_hexane", self.dry_piece_count),
            ("dry_total_stefan_flux", self.dry_face_count),
            ("event_time", 1),
            ("interface_temperature", 1),
        )

    @property
    def residual_blocks(self) -> tuple[tuple[str, int], ...]:
        return (
            ("wet_water", self.wet_piece_count),
            ("wet_energy", self.wet_piece_count),
            ("dry_water", self.dry_piece_count),
            ("dry_hexane", self.dry_piece_count),
            ("dry_energy", self.dry_piece_count),
            ("rh_water", 1),
            ("rh_hexane", 1),
            ("rh_energy", 1),
        )

    @property
    def unknown_count(self) -> int:
        return sum(size for _, size in self.unknown_blocks)

    @property
    def residual_count(self) -> int:
        return sum(size for _, size in self.residual_blocks)

    @property
    def rank_formula(self) -> int:
        return 2 * self.wet_piece_count + 3 * self.dry_piece_count + 3

    @property
    def eliminated_zero_volume_primitive_count(self) -> int:
        # Start-dry (T,y) and endpoint-wet (T,W) for transformed cell k.
        return 4

    @property
    def is_square(self) -> bool:
        return self.unknown_count == self.residual_count == self.rank_formula


@dataclass(frozen=True)
class FaceToFaceUnknowns:
    """Exact inner-face endpoint variables in declared block order."""

    wet_temperatures_k: tuple[float, ...]
    wet_retained_water_loadings: tuple[float, ...]
    dry_temperatures_k: tuple[float, ...]
    dry_y_hexane: tuple[float, ...]
    dry_total_stefan_fluxes_mol_m2_s: tuple[float, ...]
    event_time_s: float
    interface_temperature_k: float
    physically_qualifying: bool = field(default=False, init=False)

    def vector(self) -> tuple[float, ...]:
        return (
            *self.wet_temperatures_k,
            *self.wet_retained_water_loadings,
            *self.dry_temperatures_k,
            *self.dry_y_hexane,
            *self.dry_total_stefan_fluxes_mol_m2_s,
            self.event_time_s,
            self.interface_temperature_k,
        )

    @classmethod
    def from_vector(
        cls,
        layout: FaceToFaceLayout,
        values: Sequence[float],
    ) -> "FaceToFaceUnknowns":
        vector = tuple(float(value) for value in values)
        if len(vector) != layout.unknown_count:
            raise ValueError("face-to-face vector has the wrong rank")
        if not all(math.isfinite(value) for value in vector):
            raise ValueError("face-to-face vector must remain finite")
        nw = layout.wet_piece_count
        nd = layout.dry_piece_count
        nf = layout.dry_face_count
        cursor = 0

        def take(count: int) -> tuple[float, ...]:
            nonlocal cursor
            result = vector[cursor : cursor + count]
            cursor += count
            return result

        return cls(
            wet_temperatures_k=take(nw),
            wet_retained_water_loadings=take(nw),
            dry_temperatures_k=take(nd),
            dry_y_hexane=take(nd),
            dry_total_stefan_fluxes_mol_m2_s=take(nf),
            event_time_s=take(1)[0],
            interface_temperature_k=take(1)[0],
        )


@dataclass(frozen=True)
class FaceToFaceInventoryChanges:
    """Per-cell changes with no endpoint zero-volume subtraction."""

    wet_water_mol: tuple[float, ...]
    wet_energy_j: tuple[float, ...]
    dry_water_mol: tuple[float, ...]
    dry_hexane_mol: tuple[float, ...]
    dry_energy_j: tuple[float, ...]
    transformed_old_dry_volume_m3: float = field(default=0.0, init=False)
    transformed_new_wet_volume_m3: float = field(default=0.0, init=False)
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceToFaceRankAudit:
    """Structural count and unit-invariant local physical-Jacobian audit."""

    unknown_count: int
    residual_count: int
    rank_formula: int
    eliminated_zero_volume_primitive_count: int
    physical_variable_labels: tuple[str, ...]
    physical_variable_scales: tuple[float, ...]
    physical_jacobian_perturbations: tuple[float, ...]
    jacobian_certificate: conditioning.UnitInvariantJacobianCertificate
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def structurally_square(self) -> bool:
        return self.unknown_count == self.residual_count == self.rank_formula

    @property
    def numerically_full_rank(self) -> bool:
        return self.jacobian_certificate.full_rank


@dataclass(frozen=True)
class FaceToFaceAssembly:
    """Complete adjacent-face residual evaluation without commitment."""

    before: fi.FaceArrivalCommittedState
    candidate: FaceToFaceUnknowns
    layout: FaceToFaceLayout
    event_geometry: cg.CutGeometry
    swept_geometry: cg.SweptCutGeometry
    interface: arrival.FaceArrivalInterfaceState
    old_wet_cells: tuple[ww.WetWaterCellState, ...]
    old_dry_cells: tuple[cp.EquilibriumPoreState, ...]
    wet_cells: tuple[ww.WetWaterCellState, ...]
    dry_cells: tuple[cp.EquilibriumPoreState, ...]
    wet_face_fluxes: tuple[ww.WetWaterFaceFlux, ...]
    dry_face_fluxes: tuple[cut.DryFaceFlux, ...]
    inventory_changes: FaceToFaceInventoryChanges
    residuals: arrival.FaceArrivalResidualBlocks
    ledger: arrival.FaceArrivalLedger
    surface_boundary: ct.PoreBoundary
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.layout.is_square:
            raise RuntimeError("face-to-face assembly lost square structural rank")
        if len(self.candidate.vector()) != self.layout.unknown_count:
            raise RuntimeError("face-to-face candidate lost declared rank")
        if len(self.residuals.vector) != self.layout.residual_count:
            raise RuntimeError("face-to-face residual lost declared rank")


def layout_for_face_to_face(
    before: fi.FaceArrivalCommittedState,
) -> FaceToFaceLayout:
    """Return the adjacent inner-face endpoint layout from exact face ``K``."""

    if not isinstance(before, fi.FaceArrivalCommittedState):
        raise TypeError("face-to-face event requires FaceArrivalCommittedState")
    start = before.arrival_face_index
    n = before.geometry.master_grid.n
    if (
        before.geometry.front.regime != "partial"
        or not before.geometry.front.at_master_face
        or before.geometry.front.master_face_index != start
        or before.geometry.cut_cell_index is not None
    ):
        raise FaceToFaceTopologyError("before state is not an exact interior face")
    if start <= 1:
        raise FaceToFaceTopologyError(
            "face K=1 leads to extinction, not another interior face"
        )
    end = start - 1
    return FaceToFaceLayout(
        wet_cell_indices=tuple(range(end)),
        dry_cell_indices=tuple(range(end, n)),
        departure_face_index=start,
        arrival_face_index=end,
        transformed_cell_index=end,
        master_cell_count=n,
    )


def assemble_face_to_face(
    before: fi.FaceArrivalCommittedState,
    candidate: FaceToFaceUnknowns,
    surface_boundary: ct.PoreBoundary,
    *,
    enforce_reduced_film_thresholds: bool = True,
) -> FaceToFaceAssembly:
    """Assemble one exact adjacent-face event or roll back exactly."""

    try:
        return _assemble_face_to_face(
            before,
            candidate,
            surface_boundary,
            enforce_reduced_film_thresholds,
        )
    except FaceToFaceStepError:
        raise
    except Exception as exc:
        if not isinstance(before, fi.FaceArrivalCommittedState):
            raise
        raise FaceToFaceStepError(
            f"face-to-face residual rejected with exact rollback: {exc}",
            before,
        ) from exc


def _assemble_face_to_face(
    before: fi.FaceArrivalCommittedState,
    candidate: FaceToFaceUnknowns,
    surface_boundary: ct.PoreBoundary,
    enforce_reduced_film_thresholds: bool,
) -> FaceToFaceAssembly:
    if not isinstance(candidate, FaceToFaceUnknowns):
        raise TypeError("candidate must be FaceToFaceUnknowns")
    if not isinstance(
        surface_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("face-to-face event requires a supported pore boundary")
    if not isinstance(enforce_reduced_film_thresholds, bool):
        raise TypeError("reduced-film threshold enforcement flag must be boolean")
    layout = layout_for_face_to_face(before)
    if len(candidate.vector()) != layout.unknown_count:
        raise ValueError("face-to-face candidate does not match endpoint rank")
    if not all(math.isfinite(value) for value in candidate.vector()):
        raise ValueError("face-to-face candidate must remain finite")
    if candidate.event_time_s <= 0.0:
        raise FaceToFaceTopologyError("face-to-face duration must be positive")
    _require_length(
        "endpoint wet temperatures",
        candidate.wet_temperatures_k,
        layout.wet_piece_count,
    )
    _require_length(
        "endpoint wet water",
        candidate.wet_retained_water_loadings,
        layout.wet_piece_count,
    )
    _require_length(
        "endpoint dry temperatures",
        candidate.dry_temperatures_k,
        layout.dry_piece_count,
    )
    _require_length(
        "endpoint dry compositions",
        candidate.dry_y_hexane,
        layout.dry_piece_count,
    )
    _require_length(
        "endpoint Stefan fluxes",
        candidate.dry_total_stefan_fluxes_mol_m2_s,
        layout.dry_face_count,
    )
    cut._validate_dry_boundary(surface_boundary, before.config)  # noqa: SLF001

    grid = before.geometry.master_grid
    start = layout.departure_face_index
    end = layout.arrival_face_index
    event_geometry = cg.partition_master_grid(grid, radius_m=grid.faces[end])
    if (
        event_geometry.front.master_face_index != end
        or event_geometry.cut_cell_index is not None
    ):
        raise RuntimeError("face-to-face endpoint is not the exact inner face")
    swept = cg.swept_cut_geometry(
        before.geometry,
        event_geometry,
        candidate.event_time_s,
    )

    start_layout = arrival.FaceArrivalLayout(
        wet_cell_indices=tuple(range(start)),
        dry_cell_indices=tuple(range(start, grid.n)),
        arrival_face_index=start,
        master_cell_count=grid.n,
    )
    old_wet = cut._evaluate_wet_piece_states(  # noqa: SLF001
        start_layout,
        before.wet_temperatures_k,
        before.wet_retained_water_loadings,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
        before.config,
    )
    old_dry = cut._evaluate_dry_piece_states(  # noqa: SLF001
        start_layout,
        before.dry_temperatures_k,
        before.dry_y_hexane,
        before.oil_fraction_labels,
        before.config,
    )
    wet_cells = cut._evaluate_wet_piece_states(  # noqa: SLF001
        layout,
        candidate.wet_temperatures_k,
        candidate.wet_retained_water_loadings,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
        before.config,
    )
    dry_cells = cut._evaluate_dry_piece_states(  # noqa: SLF001
        layout,
        candidate.dry_temperatures_k,
        candidate.dry_y_hexane,
        before.oil_fraction_labels,
        before.config,
    )
    interface = arrival.evaluate_face_interface_state(
        candidate.interface_temperature_k,
        before.config,
        before.historical_hexane_loadings[end],
        before.oil_fraction_labels[end],
        before.oil_fraction_labels[end],
    )
    wet_fluxes = cut._wet_face_fluxes(  # noqa: SLF001
        event_geometry,
        layout,
        wet_cells,
        interface,
        before.config,
    )
    dry_fluxes = cut._dry_face_fluxes(  # noqa: SLF001
        event_geometry,
        layout,
        dry_cells,
        interface,
        candidate.dry_total_stefan_fluxes_mol_m2_s,
        surface_boundary,
        before.oil_fraction_labels,
        before.config,
        enforce_reduced_film_thresholds=enforce_reduced_film_thresholds,
    )
    changes = _inventory_changes(
        before,
        layout,
        old_wet,
        old_dry,
        wet_cells,
        dry_cells,
    )
    arrival_changes = arrival.FaceArrivalInventoryChanges(
        wet_water_mol=changes.wet_water_mol,
        wet_energy_j=changes.wet_energy_j,
        dry_water_mol=changes.dry_water_mol,
        dry_hexane_mol=changes.dry_hexane_mol,
        dry_energy_j=changes.dry_energy_j,
    )
    residuals = arrival._residual_blocks(  # noqa: SLF001
        grid,
        layout,
        arrival_changes,
        interface,
        wet_fluxes,
        dry_fluxes,
        swept,
        candidate.event_time_s,
        # O9a/O10a, applied together as the O10a family treatment (rulings
        # docs/GT_PS2_ALE_DONOR_CORRECTION_RULING_PACKET_2026-08-21.md and
        # docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md).  The
        # Rankine-Hugoniot jump conditions must convect mass and enthalpy OF
        # THE SAME MATERIAL STATE - an internal-consistency requirement of the
        # balance laws - and that state is the bulk, undrained one the O9
        # traverse measured (Peclet 3.4e3).  This chart consumes the WHOLE wet
        # material cell ``k`` within the event duration, so both donors are
        # that cell's own old bulk state, ``old_wet[-1]`` - the identical state
        # ``_inventory_changes`` compares the transformed full dry inventory
        # against.
        old_wet[-1].retained_water_concentration_mol_m3,
        old_wet[-1].energy_density_j_m3,
    )
    ledger = arrival._ledger(  # noqa: SLF001
        arrival_changes,
        residuals,
        dry_fluxes[-1],
        grid.areas[-1],
        swept,
        candidate.event_time_s,
        (*wet_fluxes, *dry_fluxes),
    )
    return FaceToFaceAssembly(
        before=before,
        candidate=candidate,
        layout=layout,
        event_geometry=event_geometry,
        swept_geometry=swept,
        interface=interface,
        old_wet_cells=old_wet,
        old_dry_cells=old_dry,
        wet_cells=wet_cells,
        dry_cells=dry_cells,
        wet_face_fluxes=wet_fluxes,
        dry_face_fluxes=dry_fluxes,
        inventory_changes=changes,
        residuals=residuals,
        ledger=ledger,
        surface_boundary=surface_boundary,
    )


def _inventory_changes(
    before: fi.FaceArrivalCommittedState,
    layout: FaceToFaceLayout,
    old_wet: Sequence[ww.WetWaterCellState],
    old_dry: Sequence[cp.EquilibriumPoreState],
    new_wet: Sequence[ww.WetWaterCellState],
    new_dry: Sequence[cp.EquilibriumPoreState],
) -> FaceToFaceInventoryChanges:
    grid = before.geometry.master_grid
    transformed = layout.transformed_cell_index
    wet_water = tuple(
        grid.volumes[index]
        * (
            new_cell.retained_water_concentration_mol_m3
            - old_cell.retained_water_concentration_mol_m3
        )
        for index, old_cell, new_cell in zip(
            layout.wet_cell_indices,
            old_wet[: layout.wet_piece_count],
            new_wet,
        )
    )
    wet_energy = tuple(
        grid.volumes[index]
        * (new_cell.energy_density_j_m3 - old_cell.energy_density_j_m3)
        for index, old_cell, new_cell in zip(
            layout.wet_cell_indices,
            old_wet[: layout.wet_piece_count],
            new_wet,
        )
    )

    old_transformed = old_wet[-1]
    new_transformed = new_dry[0]
    volume = grid.volumes[transformed]
    old_hexane = (
        before.config.wet.wet.rho_dm_p
        * before.historical_hexane_loadings[transformed]
        / hx.M
    )
    dry_water = [
        volume
        * (
            new_transformed.total_water_concentration_mol_m3
            - old_transformed.retained_water_concentration_mol_m3
        )
    ]
    dry_hexane = [
        volume
        * (new_transformed.total_hexane_concentration_mol_m3 - old_hexane)
    ]
    dry_energy = [
        volume
        * (new_transformed.energy_density_j_m3 - old_transformed.energy_density_j_m3)
    ]
    for cell_index, old_cell, new_cell in zip(
        layout.dry_cell_indices[1:],
        old_dry,
        new_dry[1:],
    ):
        cell_volume = grid.volumes[cell_index]
        dry_water.append(
            cell_volume
            * (
                new_cell.total_water_concentration_mol_m3
                - old_cell.total_water_concentration_mol_m3
            )
        )
        dry_hexane.append(
            cell_volume
            * (
                new_cell.total_hexane_concentration_mol_m3
                - old_cell.total_hexane_concentration_mol_m3
            )
        )
        dry_energy.append(
            cell_volume
            * (new_cell.energy_density_j_m3 - old_cell.energy_density_j_m3)
        )
    return FaceToFaceInventoryChanges(
        wet_water_mol=wet_water,
        wet_energy_j=wet_energy,
        dry_water_mol=tuple(dry_water),
        dry_hexane_mol=tuple(dry_hexane),
        dry_energy_j=tuple(dry_energy),
    )


def datum_covariant_residual_vector(
    residuals: arrival.FaceArrivalResidualBlocks,
) -> tuple[float, ...]:
    """Return the unchanged invertible common-datum energy-row basis."""

    return fi._datum_covariant_residual_vector(residuals)  # noqa: SLF001


def audit_face_to_face_rank(
    assembly: FaceToFaceAssembly,
    residual_scales: Sequence[float],
    physical_variable_labels: Sequence[str],
    physical_variable_scales: Sequence[float],
    *,
    relative_probe: float = 1.0e-5,
    rank_relative_tolerance: float = 1.0e-13,
) -> FaceToFaceRankAudit:
    """Audit the physical point Jacobian without claiming a time enclosure."""

    try:
        if not isinstance(assembly, FaceToFaceAssembly):
            raise TypeError("rank audit requires FaceToFaceAssembly")
        layout = assembly.layout
        row_scales = tuple(float(value) for value in residual_scales)
        labels = tuple(str(value) for value in physical_variable_labels)
        variable_scales = tuple(float(value) for value in physical_variable_scales)
        if len(row_scales) != layout.residual_count or not all(
            math.isfinite(value) and value > 0.0 for value in row_scales
        ):
            raise ValueError("rank audit residual scales lost row rank")
        if (
            len(labels) != layout.unknown_count
            or len(variable_scales) != layout.unknown_count
            or not all(
                math.isfinite(value) and value > 0.0
                for value in variable_scales
            )
        ):
            raise ValueError("rank audit variable contract lost unknown rank")
        if not math.isfinite(relative_probe) or relative_probe <= 0.0:
            raise ValueError("rank probe must be positive")

        center = np.asarray(assembly.candidate.vector(), dtype=float)
        matrix = np.empty((layout.residual_count, layout.unknown_count), dtype=float)
        perturbations = np.empty(layout.unknown_count, dtype=float)
        for column, characteristic in enumerate(variable_scales):
            delta = relative_probe * characteristic
            for _ in range(20):
                upper = center.copy()
                lower = center.copy()
                upper[column] += delta
                lower[column] -= delta
                try:
                    upper_candidate = FaceToFaceUnknowns.from_vector(layout, upper)
                    lower_candidate = FaceToFaceUnknowns.from_vector(layout, lower)
                    upper_residual = np.asarray(
                        datum_covariant_residual_vector(
                            assemble_face_to_face(
                                assembly.before,
                                upper_candidate,
                                assembly.surface_boundary,
                                enforce_reduced_film_thresholds=False,
                            ).residuals
                        ),
                        dtype=float,
                    )
                    lower_residual = np.asarray(
                        datum_covariant_residual_vector(
                            assemble_face_to_face(
                                assembly.before,
                                lower_candidate,
                                assembly.surface_boundary,
                                enforce_reduced_film_thresholds=False,
                            ).residuals
                        ),
                        dtype=float,
                    )
                except (ValueError, RuntimeError, FaceToFaceStepError):
                    delta *= 0.5
                    continue
                matrix[:, column] = (
                    upper_residual - lower_residual
                ) / (2.0 * delta)
                perturbations[column] = delta
                break
            else:
                raise ValueError(
                    "face-to-face physical Jacobian has no phase-safe central probe"
                )
        certificate = conditioning.certify_unit_invariant_jacobian(
            matrix,
            row_scales,
            variable_scales,
            rank_relative_tolerance=rank_relative_tolerance,
        )
        return FaceToFaceRankAudit(
            unknown_count=layout.unknown_count,
            residual_count=layout.residual_count,
            rank_formula=layout.rank_formula,
            eliminated_zero_volume_primitive_count=(
                layout.eliminated_zero_volume_primitive_count
            ),
            physical_variable_labels=labels,
            physical_variable_scales=variable_scales,
            physical_jacobian_perturbations=tuple(
                float(value) for value in perturbations
            ),
            jacobian_certificate=certificate,
        )
    except Exception as exc:
        before = getattr(assembly, "before", None)
        if not isinstance(before, fi.FaceArrivalCommittedState):
            raise
        raise FaceToFaceStepError(
            f"face-to-face rank audit rejected with exact rollback: {exc}",
            before,
        ) from exc


def _require_length(name: str, values: Sequence[float], expected: int) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} has length {len(values)}, expected {expected}")


__all__ = [
    "FaceToFaceAssembly",
    "FaceToFaceInventoryChanges",
    "FaceToFaceLayout",
    "FaceToFaceRankAudit",
    "FaceToFaceStepError",
    "FaceToFaceTopologyError",
    "FaceToFaceUnknowns",
    "assemble_face_to_face",
    "audit_face_to_face_rank",
    "datum_covariant_residual_vector",
    "layout_for_face_to_face",
]
