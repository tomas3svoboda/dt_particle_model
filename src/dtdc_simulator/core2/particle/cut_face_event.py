r"""Exact inner-master-face arrival residual for Gate 1g.

The incoming state has a strict partial front in material cell ``k > 0``.
The event state lies exactly on inner master face ``k``: cells ``0..k-1`` are
full wet cells and cells ``k..N-1`` are full dry cells.  No epsilon phase
piece is retained.

The wet cut piece that disappears at arrival has no event-state primitive and
therefore no event-state row.  Its old inventory is not discarded.  The first
new full-dry-cell balance is the conservative algebraic sum of the vanished
wet-piece balance, the old dry-piece balance, and the corresponding jump: its
new full-cell inventory is compared with *both* old subpiece inventories, and
its boundary rates are the fixed wet material-face rate and the outer dry
rate.  The three ALE Rankine--Hugoniot equations remain as independent event
conditions, using the exact spherical sweep divided by the unknown event
time.

This is a residual foundation, not an event solver or performance
qualification.  All coefficients and boundary traces retain their numerical,
non-qualifying status.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.particle.grid import SphericalGrid
from dtdc_simulator.core2.props import hexane as hx


class FaceArrivalTopologyError(ValueError):
    """The proposal is not the exact inner-face event chart."""


class FaceArrivalStepError(RuntimeError):
    """Rejected event residual carrying the exact immutable rollback state."""

    def __init__(
        self,
        message: str,
        rollback_state: cut.CutTransportState,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state


@dataclass(frozen=True)
class FaceArrivalLayout:
    """Rank-reduced positive-region layout at inner face ``k``."""

    wet_cell_indices: tuple[int, ...]
    dry_cell_indices: tuple[int, ...]
    arrival_face_index: int
    master_cell_count: int
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        k = self.arrival_face_index
        if not 0 < k < self.master_cell_count:
            raise ValueError("inner-face arrival requires 0 < k < N")
        if self.wet_cell_indices != tuple(range(k)):
            raise ValueError("event wet cells must be exactly 0..k-1")
        if self.dry_cell_indices != tuple(range(k, self.master_cell_count)):
            raise ValueError("event dry cells must be exactly k..N-1")

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
    def is_square(self) -> bool:
        return self.unknown_count == self.residual_count == self.rank_formula


@dataclass(frozen=True)
class FaceArrivalUnknowns:
    """Event-state primitives in the documented rank ordering."""

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
        layout: FaceArrivalLayout,
        values: Sequence[float],
    ) -> "FaceArrivalUnknowns":
        vector = tuple(values)
        if len(vector) != layout.unknown_count:
            raise ValueError(
                f"event vector has length {len(vector)}, expected {layout.unknown_count}"
            )
        if not all(math.isfinite(value) for value in vector):
            raise ValueError("event candidate vector must be finite")
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
class FaceArrivalInterfaceState:
    """Split material-label traces at the exact master face."""

    temperature_k: float
    y_hexane: float
    log_hexane_fugacity_residual: float
    log_retained_water_fugacity_residual: float
    retained_water_trace_capacity_dual_over_rt: float
    dry_actual_water_potential_isothermal: float
    wet: ww.WetWaterCellState
    dry: cp.EquilibriumPoreState
    wet_historical_hexane_loading: float
    wet_oil_fraction_label: float
    dry_oil_fraction_label: float
    interface_composition_root_audit: cut.InterfaceCompositionRootAudit
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceArrivalInventoryChanges:
    """Cancellation-safe material inventory changes over the event interval."""

    wet_water_mol: tuple[float, ...]
    wet_energy_j: tuple[float, ...]
    dry_water_mol: tuple[float, ...]
    dry_hexane_mol: tuple[float, ...]
    dry_energy_j: tuple[float, ...]
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceArrivalResidualBlocks:
    """Square rank-reduced regional balances and three ALE jumps."""

    wet_water_mol_s: tuple[float, ...]
    wet_energy_w: tuple[float, ...]
    dry_water_mol_s: tuple[float, ...]
    dry_hexane_mol_s: tuple[float, ...]
    dry_energy_w: tuple[float, ...]
    rh_water_mol_s: float
    rh_hexane_mol_s: float
    rh_energy_w: float
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def vector(self) -> tuple[float, ...]:
        return (
            *self.wet_water_mol_s,
            *self.wet_energy_w,
            *self.dry_water_mol_s,
            *self.dry_hexane_mol_s,
            *self.dry_energy_w,
            self.rh_water_mol_s,
            self.rh_hexane_mol_s,
            self.rh_energy_w,
        )


@dataclass(frozen=True)
class FaceArrivalLedger:
    """Independent global telescoping and geometry checks."""

    water_global_from_regions_mol_s: float
    water_global_direct_mol_s: float
    water_telescoping_error_mol_s: float
    hexane_global_from_regions_mol_s: float
    hexane_global_direct_mol_s: float
    hexane_telescoping_error_mol_s: float
    energy_global_from_regions_w: float
    energy_global_direct_w: float
    energy_telescoping_error_w: float
    fixed_wet_hexane_flux_mol_s: float
    geometry_cell_gcl_residual_m3: float
    minimum_entropy_production_w_m3_k: float
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceArrivalAssembly:
    """Complete exact-face residual evaluation without state commitment."""

    before: cut.CutTransportState
    candidate: FaceArrivalUnknowns
    layout: FaceArrivalLayout
    event_geometry: cg.CutGeometry
    swept_geometry: cg.SweptCutGeometry
    interface: FaceArrivalInterfaceState
    old_wet_cells: tuple[ww.WetWaterCellState, ...]
    old_dry_cells: tuple[cp.EquilibriumPoreState, ...]
    wet_cells: tuple[ww.WetWaterCellState, ...]
    dry_cells: tuple[cp.EquilibriumPoreState, ...]
    wet_face_fluxes: tuple[ww.WetWaterFaceFlux, ...]
    dry_face_fluxes: tuple[cut.DryFaceFlux, ...]
    inventory_changes: FaceArrivalInventoryChanges
    residuals: FaceArrivalResidualBlocks
    ledger: FaceArrivalLedger
    surface_boundary: ct.PoreBoundary
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.layout.is_square:
            raise RuntimeError("face-arrival residual lost its square rank")
        if len(self.candidate.vector()) != self.layout.unknown_count:
            raise RuntimeError("face-arrival candidate does not match declared rank")
        if len(self.residuals.vector) != self.layout.residual_count:
            raise RuntimeError("face-arrival residual does not match declared rank")


def layout_for_arrival(before: cut.CutTransportState) -> FaceArrivalLayout:
    """Return the exact-face rank after a strict cut in cell ``k > 0``."""

    if not isinstance(before, cut.CutTransportState):
        raise TypeError("face arrival requires a CutTransportState")
    geometry = before.geometry
    k = geometry.cut_cell_index
    if geometry.front.regime != "partial" or k is None:
        raise FaceArrivalTopologyError(
            "face arrival requires a strict positive-volume cut before state"
        )
    if k <= 0:
        raise FaceArrivalTopologyError(
            "cell-zero arrival is the wet-domain extinction event, not this chart"
        )
    layout = FaceArrivalLayout(
        wet_cell_indices=tuple(range(k)),
        dry_cell_indices=tuple(range(k, geometry.master_grid.n)),
        arrival_face_index=k,
        master_cell_count=geometry.master_grid.n,
    )
    if not layout.is_square:
        raise RuntimeError("face-arrival rank formula is internally inconsistent")
    return layout


def evaluate_face_interface_state(
    temperature_k: float,
    config: cut.CutTransportConfig,
    wet_historical_hexane_loading: float,
    wet_oil_fraction_label: float,
    dry_oil_fraction_label: float,
) -> FaceArrivalInterfaceState:
    """Evaluate adjacent wet/dry traces without conflating material labels."""

    values = (
        temperature_k,
        wet_historical_hexane_loading,
        wet_oil_fraction_label,
        dry_oil_fraction_label,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("face-interface temperature and labels must be finite")
    if wet_historical_hexane_loading < 0.0:
        raise ValueError("wet historical n-hexane label must be non-negative")
    if wet_oil_fraction_label != dry_oil_fraction_label:
        raise FaceArrivalTopologyError(
            "exact arrival wet and dry traces must use one consumed-shell oil label"
        )
    if not config.wet.wet.T_min <= temperature_k <= config.wet.wet.T_max:
        raise FaceArrivalTopologyError(
            "face-interface temperature left the wet caloric bracket"
        )
    pore_lo, pore_hi = config.dry.pore.temperature_bounds_k
    if not pore_lo <= temperature_k <= pore_hi:
        raise FaceArrivalTopologyError(
            "face-interface temperature left the pore-state bracket"
        )

    # Use the same frozen root authority as the within-cell cut chart.
    root = cut._solve_interface_composition_with_audit(  # noqa: SLF001
        temperature_k,
        config,
        dry_oil_fraction_label,
    )
    y_hexane = root.y_hexane
    hexane_residual = root.log_fugacity_residual
    dry = root.dry
    trace_capacity_dual, dry_actual_water_potential = (
        cut._retained_water_interface_trace_kkt(dry, config)  # noqa: SLF001
    )
    wet = ww.evaluate_cell(
        temperature_k,
        dry.retained_water_loading,
        wet_historical_hexane_loading,
        wet_oil_fraction_label,
        config.wet,
        retained_water_capacity_dual_over_rt=trace_capacity_dual,
    )
    water_residual = dry_actual_water_potential - wet.retained_water_potential
    if abs(water_residual) > 1.0e-10:
        raise FaceArrivalTopologyError(
            "split-label wet/dry retained-water potentials lost their shared authority"
        )
    return FaceArrivalInterfaceState(
        temperature_k=temperature_k,
        y_hexane=y_hexane,
        log_hexane_fugacity_residual=hexane_residual,
        log_retained_water_fugacity_residual=water_residual,
        retained_water_trace_capacity_dual_over_rt=trace_capacity_dual,
        dry_actual_water_potential_isothermal=dry_actual_water_potential,
        wet=wet,
        dry=dry,
        wet_historical_hexane_loading=wet_historical_hexane_loading,
        wet_oil_fraction_label=wet_oil_fraction_label,
        dry_oil_fraction_label=dry_oil_fraction_label,
        interface_composition_root_audit=root.audit,
    )


def assemble_face_arrival(
    before: cut.CutTransportState,
    candidate: FaceArrivalUnknowns,
    surface_boundary: ct.PoreBoundary,
    *,
    enforce_reduced_film_thresholds: bool = True,
) -> FaceArrivalAssembly:
    """Assemble the exact-face event residual or reject with exact rollback."""

    try:
        return _assemble_face_arrival(
            before,
            candidate,
            surface_boundary,
            enforce_reduced_film_thresholds,
        )
    except FaceArrivalStepError:
        raise
    except Exception as exc:
        if not isinstance(before, cut.CutTransportState):
            raise
        raise FaceArrivalStepError(
            f"face-arrival residual rejected with exact rollback: {exc}",
            before,
        ) from exc


def _assemble_face_arrival(
    before: cut.CutTransportState,
    candidate: FaceArrivalUnknowns,
    surface_boundary: ct.PoreBoundary,
    enforce_reduced_film_thresholds: bool,
) -> FaceArrivalAssembly:
    if not isinstance(before, cut.CutTransportState):
        raise TypeError("before must be a CutTransportState")
    if not isinstance(candidate, FaceArrivalUnknowns):
        raise TypeError("candidate must be FaceArrivalUnknowns")
    if not isinstance(
        surface_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("surface boundary must be a supported pore boundary")
    if not isinstance(enforce_reduced_film_thresholds, bool):
        raise TypeError("reduced-film threshold enforcement flag must be boolean")
    layout = layout_for_arrival(before)
    if len(candidate.vector()) != layout.unknown_count:
        raise ValueError("event candidate does not match the exact-face rank")
    if not all(math.isfinite(value) for value in candidate.vector()):
        raise ValueError("event candidate values must be finite")
    if candidate.event_time_s <= 0.0:
        raise FaceArrivalTopologyError("event time must be strictly positive")
    _require_length(
        "event wet temperatures",
        candidate.wet_temperatures_k,
        layout.wet_piece_count,
    )
    _require_length(
        "event wet retained water",
        candidate.wet_retained_water_loadings,
        layout.wet_piece_count,
    )
    _require_length(
        "event dry temperatures",
        candidate.dry_temperatures_k,
        layout.dry_piece_count,
    )
    _require_length(
        "event dry compositions",
        candidate.dry_y_hexane,
        layout.dry_piece_count,
    )
    _require_length(
        "event dry Stefan fluxes",
        candidate.dry_total_stefan_fluxes_mol_m2_s,
        layout.dry_face_count,
    )
    cut._validate_dry_boundary(surface_boundary, before.config)  # noqa: SLF001
    if any(
        dual != 0.0
        for dual in before.effective_wet_retained_water_capacity_duals_over_rt
    ):
        raise FaceArrivalTopologyError(
            "exact face arrival cannot yet continue a positive-volume "
            "capacity-active wet state; semismooth face-event coordinates are required"
        )

    grid = before.geometry.master_grid
    k = layout.arrival_face_index
    event_geometry = cg.partition_master_grid(grid, radius_m=grid.faces[k])
    if (
        not event_geometry.front.at_master_face
        or event_geometry.front.master_face_index != k
        or event_geometry.cut_cell_index is not None
    ):
        raise RuntimeError("event geometry is not the exact requested master face")
    swept = cg.swept_cut_geometry(
        before.geometry,
        event_geometry,
        candidate.event_time_s,
    )

    before_layout = before.layout
    old_wet = cut._evaluate_wet_piece_states(  # noqa: SLF001
        before_layout,
        before.wet_temperatures_k,
        before.wet_retained_water_loadings,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
        before.config,
    )
    old_dry = cut._evaluate_dry_piece_states(  # noqa: SLF001
        before_layout,
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
    # Arrival is the one-sided limit of the front while it consumes material
    # shell k.  Preserve that disappearing shell's historical wet trace in the
    # RH state; shell k-1 becomes authoritative only after departure into it.
    interface = evaluate_face_interface_state(
        candidate.interface_temperature_k,
        before.config,
        before.historical_hexane_loadings[k],
        before.oil_fraction_labels[k],
        before.oil_fraction_labels[k],
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
    residuals = _residual_blocks(
        grid,
        layout,
        changes,
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
        # traverse measured (Peclet 3.4e3).  The front consumes the whole
        # remaining wet subpiece of shell ``k`` within the event time, so both
        # donors are that subpiece's OWN old bulk state, which is exactly the
        # state ``_inventory_changes`` retires into the first full dry cell.
        old_wet[-1].retained_water_concentration_mol_m3,
        old_wet[-1].energy_density_j_m3,
    )
    ledger = _ledger(
        changes,
        residuals,
        dry_fluxes[-1],
        grid.areas[-1],
        swept,
        candidate.event_time_s,
        (*wet_fluxes, *dry_fluxes),
    )
    return FaceArrivalAssembly(
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
    before: cut.CutTransportState,
    layout: FaceArrivalLayout,
    old_wet: Sequence[ww.WetWaterCellState],
    old_dry: Sequence[cp.EquilibriumPoreState],
    new_wet: Sequence[ww.WetWaterCellState],
    new_dry: Sequence[cp.EquilibriumPoreState],
) -> FaceArrivalInventoryChanges:
    grid = before.geometry.master_grid
    k = layout.arrival_face_index

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

    old_cut_wet_volume = before.geometry.cells[k].wet_volume_m3
    old_cut_dry_volume = before.geometry.cells[k].dry_volume_m3
    old_cut_wet = old_wet[-1]
    old_cut_dry = old_dry[0]
    new_first_dry = new_dry[0]
    old_cut_wet_hexane = (
        before.config.wet.wet.rho_dm_p
        * before.historical_hexane_loadings[k]
        / hx.M
    )

    def combined_change(
        new_concentration: float,
        old_wet_concentration: float,
        old_dry_concentration: float,
    ) -> float:
        # V_full*C_new - V_w*C_w_old - V_d*C_d_old, evaluated without
        # subtracting nearly equal full-cell inventories as V_w -> 0.
        return math.fsum(
            (
                old_cut_wet_volume
                * (new_concentration - old_wet_concentration),
                old_cut_dry_volume
                * (new_concentration - old_dry_concentration),
            )
        )

    first_water = combined_change(
        new_first_dry.total_water_concentration_mol_m3,
        old_cut_wet.retained_water_concentration_mol_m3,
        old_cut_dry.total_water_concentration_mol_m3,
    )
    first_hexane = combined_change(
        new_first_dry.total_hexane_concentration_mol_m3,
        old_cut_wet_hexane,
        old_cut_dry.total_hexane_concentration_mol_m3,
    )
    first_energy = combined_change(
        new_first_dry.energy_density_j_m3,
        old_cut_wet.energy_density_j_m3,
        old_cut_dry.energy_density_j_m3,
    )

    dry_water = [first_water]
    dry_hexane = [first_hexane]
    dry_energy = [first_energy]
    for cell_index, old_cell, new_cell in zip(
        layout.dry_cell_indices[1:],
        old_dry[1:],
        new_dry[1:],
    ):
        volume = grid.volumes[cell_index]
        dry_water.append(
            volume
            * (
                new_cell.total_water_concentration_mol_m3
                - old_cell.total_water_concentration_mol_m3
            )
        )
        dry_hexane.append(
            volume
            * (
                new_cell.total_hexane_concentration_mol_m3
                - old_cell.total_hexane_concentration_mol_m3
            )
        )
        dry_energy.append(
            volume * (new_cell.energy_density_j_m3 - old_cell.energy_density_j_m3)
        )
    return FaceArrivalInventoryChanges(
        wet_water_mol=wet_water,
        wet_energy_j=wet_energy,
        dry_water_mol=tuple(dry_water),
        dry_hexane_mol=tuple(dry_hexane),
        dry_energy_j=tuple(dry_energy),
    )


def _residual_blocks(
    grid: SphericalGrid,
    layout: FaceArrivalLayout,
    changes: FaceArrivalInventoryChanges,
    interface: FaceArrivalInterfaceState,
    wet_fluxes: Sequence[ww.WetWaterFaceFlux],
    dry_fluxes: Sequence[cut.DryFaceFlux],
    swept: cg.SweptCutGeometry,
    event_time_s: float,
    swept_wet_water_donor_concentration_mol_m3: float,
    swept_wet_energy_donor_density_j_m3: float,
) -> FaceArrivalResidualBlocks:
    areas = grid.areas
    k = layout.arrival_face_index
    area_gamma = areas[k]
    q_gamma = swept.interface_swept_volume_rate_m3_s

    # O9a/O10a ALE FRONT-DONOR CORRECTION AND ITS ENERGY COMPLETION, ruled
    # 2026-08-21 in docs/GT_PS2_ALE_DONOR_CORRECTION_RULING_PACKET_2026-08-21.md
    # and docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md.  O10a
    # rules the two together as ONE COHERENT FAMILY TREATMENT, which is what
    # reaches this chart: Stage 1b measured that moving the mass donor alone
    # here destroys the water-datum gauge covariance of ``rh_energy``, because
    # the two RH rows would then convect different material states.
    #
    # PHYSICS GROUND, as recorded in the O10a ruling: two-phase moving-boundary
    # (Rankine-Hugoniot) jump conditions must convect mass and enthalpy OF THE
    # SAME MATERIAL STATE - an internal-consistency requirement of the balance
    # laws, not an empirical claim - and WHICH state it is, the bulk undrained
    # one, is the O9-measured fact (Peclet 3.4e3 forbids pre-drainage).  These
    # three rows ARE those jump conditions, so both donors move together.
    #
    # THE CONSUMED MATERIAL IS THIS CHART'S OWN, NOT A NEW CHOICE.  At arrival
    # the front consumes the whole remaining wet subpiece of shell ``k`` inside
    # the event time, so ``q_gamma`` is exactly ``-V_w_old / t`` and the swept
    # material is that subpiece at its OWN old bulk state - the identical state
    # ``_inventory_changes`` already books through ``combined_change`` when it
    # retires ``V_w_old * C_w_old`` into the first full dry cell.  Both donors
    # are therefore read from ``old_wet[-1]``, at the same time level and from
    # the same material state as each other and as the storage side.
    #
    # ROLE SEPARATION IS UNCHANGED.  The zero-volume equilibrium trace keeps
    # its potential/flux role in full: ``interface.wet`` is still the right
    # state of the front-face diffusive water and energy fluxes, and the
    # ``rh_hexane`` row keeps the disappearing shell's historical hexane trace
    # (which is that shell's own bulk loading anyway).
    wet_water_donor_concentration = swept_wet_water_donor_concentration_mol_m3
    wet_energy_donor_density = swept_wet_energy_donor_density_j_m3

    wet_interface_water_lab = (
        area_gamma * wet_fluxes[-1].retained_water_flux_mol_m2_s
    )
    wet_interface_energy_lab = area_gamma * wet_fluxes[-1].total_energy_flux_w_m2
    wet_hexane_concentration = (
        interface.wet.total_hexane_concentration_kg_m3 / hx.M
    )

    dry_interface = dry_fluxes[0]
    rh_water = math.fsum(
        (
            area_gamma
            * dry_interface.component.conserved_water_flux_mol_m2_s,
            -wet_interface_water_lab,
            -q_gamma * interface.dry.total_water_concentration_mol_m3,
            # O9a, carried here by the O10a family treatment.
            q_gamma * wet_water_donor_concentration,
        )
    )
    rh_hexane = math.fsum(
        (
            area_gamma
            * dry_interface.component.conserved_hexane_flux_mol_m2_s,
            -q_gamma * interface.dry.total_hexane_concentration_mol_m3,
            q_gamma * wet_hexane_concentration,
        )
    )
    rh_energy = math.fsum(
        (
            area_gamma * dry_interface.energy.total_energy_flux_w_m2,
            -wet_interface_energy_lab,
            -q_gamma * interface.dry.energy_density_j_m3,
            # O10a: the SAME material state as ``rh_water``'s donor above, so
            # the water caloric datum stays an exact gauge freedom of this row.
            q_gamma * wet_energy_donor_density,
        )
    )

    wet_water_rates = [0.0]
    wet_energy_rates = [0.0]
    for face_position in range(1, layout.wet_piece_count):
        wet_water_rates.append(
            areas[face_position]
            * wet_fluxes[face_position].retained_water_flux_mol_m2_s
        )
        wet_energy_rates.append(
            areas[face_position] * wet_fluxes[face_position].total_energy_flux_w_m2
        )
    wet_water_rates.append(wet_interface_water_lab)
    wet_energy_rates.append(wet_interface_energy_lab)
    wet_water_residuals = tuple(
        math.fsum(
            (
                change / event_time_s,
                wet_water_rates[index + 1],
                -wet_water_rates[index],
            )
        )
        for index, change in enumerate(changes.wet_water_mol)
    )
    wet_energy_residuals = tuple(
        math.fsum(
            (
                change / event_time_s,
                wet_energy_rates[index + 1],
                -wet_energy_rates[index],
            )
        )
        for index, change in enumerate(changes.wet_energy_j)
    )

    dry_water_rates = tuple(
        areas[k + face]
        * flux.component.conserved_water_flux_mol_m2_s
        for face, flux in enumerate(dry_fluxes)
    )
    dry_hexane_rates = tuple(
        areas[k + face]
        * flux.component.conserved_hexane_flux_mol_m2_s
        for face, flux in enumerate(dry_fluxes)
    )
    dry_energy_rates = tuple(
        areas[k + face] * flux.energy.total_energy_flux_w_m2
        for face, flux in enumerate(dry_fluxes)
    )

    def rank_reduced_region(
        inventory_changes: Sequence[float],
        dry_rates: Sequence[float],
        fixed_wet_inner_rate: float,
    ) -> tuple[float, ...]:
        residuals = [
            math.fsum(
                (
                    inventory_changes[0] / event_time_s,
                    dry_rates[1],
                    -fixed_wet_inner_rate,
                )
            )
        ]
        residuals.extend(
            math.fsum(
                (
                    inventory_changes[index] / event_time_s,
                    dry_rates[index + 1],
                    -dry_rates[index],
                )
            )
            for index in range(1, layout.dry_piece_count)
        )
        return tuple(residuals)

    dry_water_residuals = rank_reduced_region(
        changes.dry_water_mol,
        dry_water_rates,
        wet_interface_water_lab,
    )
    dry_hexane_residuals = rank_reduced_region(
        changes.dry_hexane_mol,
        dry_hexane_rates,
        0.0,
    )
    dry_energy_residuals = rank_reduced_region(
        changes.dry_energy_j,
        dry_energy_rates,
        wet_interface_energy_lab,
    )
    return FaceArrivalResidualBlocks(
        wet_water_mol_s=wet_water_residuals,
        wet_energy_w=wet_energy_residuals,
        dry_water_mol_s=dry_water_residuals,
        dry_hexane_mol_s=dry_hexane_residuals,
        dry_energy_w=dry_energy_residuals,
        rh_water_mol_s=rh_water,
        rh_hexane_mol_s=rh_hexane,
        rh_energy_w=rh_energy,
    )


def _ledger(
    changes: FaceArrivalInventoryChanges,
    residuals: FaceArrivalResidualBlocks,
    surface: cut.DryFaceFlux,
    surface_area_m2: float,
    swept: cg.SweptCutGeometry,
    event_time_s: float,
    all_fluxes: Sequence[object],
) -> FaceArrivalLedger:
    water_regions = math.fsum(
        (*residuals.wet_water_mol_s, *residuals.dry_water_mol_s)
    )
    hexane_regions = math.fsum(residuals.dry_hexane_mol_s)
    energy_regions = math.fsum(
        (*residuals.wet_energy_w, *residuals.dry_energy_w)
    )
    water_direct = math.fsum(
        (
            math.fsum((*changes.wet_water_mol, *changes.dry_water_mol))
            / event_time_s,
            surface_area_m2 * surface.component.conserved_water_flux_mol_m2_s,
        )
    )
    hexane_direct = math.fsum(
        (
            math.fsum(changes.dry_hexane_mol) / event_time_s,
            surface_area_m2 * surface.component.conserved_hexane_flux_mol_m2_s,
        )
    )
    energy_direct = math.fsum(
        (
            math.fsum((*changes.wet_energy_j, *changes.dry_energy_j))
            / event_time_s,
            surface_area_m2 * surface.energy.total_energy_flux_w_m2,
        )
    )
    entropy_values: list[float] = []
    for flux in all_fluxes:
        if isinstance(flux, ww.WetWaterFaceFlux):
            entropy_values.append(flux.total_entropy_production_w_m3_k)
        elif isinstance(flux, cut.DryFaceFlux):
            entropy_values.append(flux.total_entropy_production_w_m3_k)
    return FaceArrivalLedger(
        water_global_from_regions_mol_s=water_regions,
        water_global_direct_mol_s=water_direct,
        water_telescoping_error_mol_s=water_regions - water_direct,
        hexane_global_from_regions_mol_s=hexane_regions,
        hexane_global_direct_mol_s=hexane_direct,
        hexane_telescoping_error_mol_s=hexane_regions - hexane_direct,
        energy_global_from_regions_w=energy_regions,
        energy_global_direct_w=energy_direct,
        energy_telescoping_error_w=energy_regions - energy_direct,
        fixed_wet_hexane_flux_mol_s=0.0,
        geometry_cell_gcl_residual_m3=swept.max_cell_gcl_residual_m3,
        minimum_entropy_production_w_m3_k=min(entropy_values, default=0.0),
    )


def _require_length(name: str, values: Sequence[float], expected: int) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} has length {len(values)}, expected {expected}")


__all__ = [
    "FaceArrivalAssembly",
    "FaceArrivalInterfaceState",
    "FaceArrivalInventoryChanges",
    "FaceArrivalLayout",
    "FaceArrivalLedger",
    "FaceArrivalResidualBlocks",
    "FaceArrivalStepError",
    "FaceArrivalTopologyError",
    "FaceArrivalUnknowns",
    "assemble_face_arrival",
    "evaluate_face_interface_state",
    "layout_for_arrival",
]
