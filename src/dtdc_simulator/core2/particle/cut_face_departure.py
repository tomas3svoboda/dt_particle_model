r"""Exact inner-master-face departure residual for Gate 1g.

The before state lies exactly on master face ``k`` with full wet material cells
``0..k-1`` and full dry cells ``k..N-1``.  A fixed positive substep moves the
front strictly into cell ``k-1``.  Cell ``k-1`` then contains its retained wet
subpiece and a first dry subpiece whose old dry volume and all old dry
inventories are exactly zero.  Existing outer dry cells keep their prior
inventories and material labels.

The first dry-subpiece storage reuses
:func:`cut_birth_event.zero_old_dry_inventories`; no epsilon state is created.
The candidate ``z`` is an ordinary nonlinear unknown.  Water, n-hexane, and
common-datum energy use the same laboratory/common-frame ALE balances and
three Rankine--Hugoniot jumps as the strict cut chart.

The narrow immutable :class:`FaceDepartureState` has two explicit adapters.
``adapt_face_arrival_assembly`` is retained for residual/manufactured tests and
does not claim acceptance.  ``adapt_face_arrival_committed_state`` accepts only
the safeguarded event integrator's certified exact-face handoff.  Neither
adapter manufactures a strict-cut state or an epsilon phase piece.

Since the O9a front-donor ruling
(docs/GT_PS2_ALE_DONOR_CORRECTION_RULING_PACKET_2026-08-21.md, ruled
2026-08-21, carried to this sibling chart in Stage 1b) the zero-volume
equilibrium interface trace is NOT the ALE mass donor here either: the receding
front consumes wet material at that material's own bulk retained-water loading,
and the excess over the dry-side equilibrium leaves through the unchanged
Rankine--Hugoniot water jump.  The trace keeps its potential/flux role in full.
See :func:`_residual_blocks`.

Nor is it the ALE ENERGY donor.  Since the O10a completion
(docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md, ruled
2026-08-21) the swept material's enthalpy is convected at the SAME (bulk,
old-time) material state as its mass, because the two-phase moving-boundary
Rankine--Hugoniot jump conditions must convect mass and enthalpy of the same
material state - an internal-consistency requirement of the balance laws - and
that state is the bulk, undrained one the O9 traverse measured (Peclet 3.4e3).

This module is a residual foundation, not an event integrator, contact/stall
law, or performance qualification.  All coefficients remain numerical and
``physically_qualifying=False``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Sequence

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_birth_event as birth
from dtdc_simulator.core2.particle import cut_face_event as arrival
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.props import hexane as hx


class FaceDepartureTopologyError(ValueError):
    """Proposal is not a strict departure into cell ``k-1``."""


class FaceDepartureStepError(RuntimeError):
    """Rejected departure residual carrying the exact immutable face state."""

    def __init__(self, message: str, rollback_state: "FaceDepartureState") -> None:
        super().__init__(message)
        self.rollback_state = rollback_state


@dataclass(frozen=True)
class FaceDepartureState:
    """Narrow exact-face state protocol pending a face-integrator commit API."""

    geometry: cg.CutGeometry
    config: cut.CutTransportConfig
    time_s: float
    wet_temperatures_k: tuple[float, ...]
    wet_retained_water_loadings: tuple[float, ...]
    dry_temperatures_k: tuple[float, ...]
    dry_y_hexane: tuple[float, ...]
    historical_hexane_loadings: tuple[float, ...]
    oil_fraction_labels: tuple[float, ...]
    source: str = "face-arrival immutable adapter"
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.geometry, cg.CutGeometry):
            raise TypeError("face state requires an exact CutGeometry")
        if not isinstance(self.config, cut.CutTransportConfig):
            raise TypeError("face state requires a CutTransportConfig")
        if not math.isfinite(self.time_s) or self.time_s < 0.0:
            raise ValueError("face-state time must be finite and non-negative")
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("face-state adapter needs provenance")
        front = self.geometry.front
        k = front.master_face_index
        if (
            front.regime != "partial"
            or not front.at_master_face
            or k is None
            or not 0 < k < self.geometry.master_grid.n
            or self.geometry.cut_cell_index is not None
        ):
            raise FaceDepartureTopologyError(
                "departure state must be an exact interior master face"
            )
        _require_length("face wet temperatures", self.wet_temperatures_k, k)
        _require_length(
            "face wet retained water",
            self.wet_retained_water_loadings,
            k,
        )
        _require_length(
            "face dry temperatures",
            self.dry_temperatures_k,
            self.geometry.master_grid.n - k,
        )
        _require_length(
            "face dry compositions",
            self.dry_y_hexane,
            self.geometry.master_grid.n - k,
        )
        _require_length(
            "face historical n-hexane labels",
            self.historical_hexane_loadings,
            self.geometry.master_grid.n,
        )
        _require_length(
            "face oil-fraction labels",
            self.oil_fraction_labels,
            self.geometry.master_grid.n,
        )
        if not all(
            math.isfinite(value) and value >= 0.0
            for value in self.historical_hexane_loadings
        ):
            raise ValueError("historical n-hexane labels must be finite and non-negative")
        if not all(math.isfinite(value) for value in self.oil_fraction_labels):
            raise ValueError("oil-fraction labels must remain finite")
        _evaluate_old_wet_cells(self)
        _evaluate_old_dry_cells(self)

    @property
    def departure_face_index(self) -> int:
        index = self.geometry.front.master_face_index
        if index is None:
            raise RuntimeError("validated face state lost its master-face index")
        return index


@dataclass(frozen=True)
class FaceDepartureLayout:
    """Strict cut layout immediately inside exact face ``k``."""

    wet_cell_indices: tuple[int, ...]
    dry_cell_indices: tuple[int, ...]
    cut_cell_index: int
    departure_face_index: int
    master_cell_count: int
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        k = self.departure_face_index
        n = self.master_cell_count
        if not 0 < k < n:
            raise ValueError("inner-face departure requires 0 < k < N")
        if self.cut_cell_index != k - 1:
            raise ValueError("departure cut cell must be k-1")
        if self.wet_cell_indices != tuple(range(k)):
            raise ValueError("departure wet cells must be exactly 0..k-1")
        if self.dry_cell_indices != tuple(range(k - 1, n)):
            raise ValueError("departure dry cells must be exactly k-1..N-1")

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
            ("front_z", 1),
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
class NewbornDryTraceGradients:
    """Cancellation-free one-sided gradients for the coalescing dry trace."""

    temperature_gradient_k_m: float
    composition_gradient_m_inv: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not all(
            math.isfinite(value)
            for value in (
                self.temperature_gradient_k_m,
                self.composition_gradient_m_inv,
            )
        ):
            raise ValueError("newborn dry trace gradients must be finite")


@dataclass(frozen=True)
class FaceDepartureUnknowns:
    """One strict post-face candidate in documented block order."""

    wet_temperatures_k: tuple[float, ...]
    wet_retained_water_loadings: tuple[float, ...]
    dry_temperatures_k: tuple[float, ...]
    dry_y_hexane: tuple[float, ...]
    dry_total_stefan_fluxes_mol_m2_s: tuple[float, ...]
    front_z: float
    interface_temperature_k: float
    physically_qualifying: bool = field(default=False, init=False)

    def vector(self) -> tuple[float, ...]:
        return (
            *self.wet_temperatures_k,
            *self.wet_retained_water_loadings,
            *self.dry_temperatures_k,
            *self.dry_y_hexane,
            *self.dry_total_stefan_fluxes_mol_m2_s,
            self.front_z,
            self.interface_temperature_k,
        )

    @classmethod
    def from_vector(
        cls,
        layout: FaceDepartureLayout,
        values: Sequence[float],
    ) -> "FaceDepartureUnknowns":
        vector = tuple(values)
        if len(vector) != layout.unknown_count:
            raise ValueError(
                f"departure vector has length {len(vector)}, expected {layout.unknown_count}"
            )
        if not all(math.isfinite(value) for value in vector):
            raise ValueError("departure candidate vector must remain finite")
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
            front_z=take(1)[0],
            interface_temperature_k=take(1)[0],
        )


@dataclass(frozen=True)
class FaceDepartureInventoryChanges:
    """Stable changes with one exact-zero-old dry subpiece."""

    wet_water_mol: tuple[float, ...]
    wet_hexane_mol: tuple[float, ...]
    wet_energy_j: tuple[float, ...]
    dry_water_mol: tuple[float, ...]
    dry_hexane_mol: tuple[float, ...]
    dry_energy_j: tuple[float, ...]
    newborn_old_dry_water_mol: float = field(default=0.0, init=False)
    newborn_old_dry_hexane_mol: float = field(default=0.0, init=False)
    newborn_old_dry_energy_j: float = field(default=0.0, init=False)
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class FaceDepartureLedger:
    """Global telescope, GCL, newborn, and material-label checks."""

    wet_hexane_sweep_identity_mol_s: float
    water_global_from_blocks_mol_s: float
    water_global_direct_mol_s: float
    water_telescoping_error_mol_s: float
    hexane_global_from_blocks_mol_s: float
    hexane_global_direct_mol_s: float
    hexane_telescoping_error_mol_s: float
    energy_global_from_blocks_w: float
    energy_global_direct_w: float
    energy_telescoping_error_w: float
    water_telescoping_relative_error: float
    hexane_telescoping_relative_error: float
    energy_telescoping_relative_error: float
    newborn_old_dry_volume_m3: float
    newborn_dry_volume_m3: float
    newborn_swept_volume_m3: float
    newborn_water_inventory_identity_mol: float
    newborn_hexane_inventory_identity_mol: float
    newborn_energy_inventory_identity_j: float
    maximum_existing_dry_volume_change_m3: float
    historical_wet_hexane_label_change: float
    geometry_cell_gcl_residual_m3: float
    minimum_entropy_production_w_m3_k: float
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def maximum_normalized_telescoping_error(self) -> float:
        return max(
            self.water_telescoping_relative_error,
            self.hexane_telescoping_relative_error,
            self.energy_telescoping_relative_error,
        )


@dataclass(frozen=True)
class FaceDepartureAssembly:
    """Complete fixed-dt inner-face departure residual evaluation."""

    before: FaceDepartureState
    candidate: FaceDepartureUnknowns
    layout: FaceDepartureLayout
    candidate_geometry: cg.CutGeometry
    swept_geometry: cg.SweptCutGeometry
    old_wet_cells: tuple[ww.WetWaterCellState, ...]
    old_dry_cells: tuple[cp.EquilibriumPoreState, ...]
    wet_cells: tuple[ww.WetWaterCellState, ...]
    dry_cells: tuple[cp.EquilibriumPoreState, ...]
    interface: cut.CutInterfaceState
    wet_face_fluxes: tuple[ww.WetWaterFaceFlux, ...]
    dry_face_fluxes: tuple[cut.DryFaceFlux, ...]
    inventory_changes: FaceDepartureInventoryChanges
    vanishing_wet_cut_treatment_audit: cut.VanishingWetCutTreatmentAudit
    residuals: cut.CutResidualBlocks
    ledger: FaceDepartureLedger
    surface_boundary: ct.PoreBoundary
    dt_s: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.layout.is_square:
            raise RuntimeError("face-departure residual lost its square rank")
        if len(self.candidate.vector()) != self.layout.unknown_count:
            raise RuntimeError("departure candidate does not match declared rank")
        if len(self.residuals.vector) != self.layout.residual_count:
            raise RuntimeError("departure residual does not match declared rank")
        if not isinstance(
            self.vanishing_wet_cut_treatment_audit,
            cut.VanishingWetCutTreatmentAudit,
        ):
            raise TypeError("face departure needs a wet-cut treatment audit")
        audit = self.vanishing_wet_cut_treatment_audit
        if (
            audit.threshold_authority
            is not self.before.config.vanishing_wet_cut_threshold_authority
        ):
            raise RuntimeError("departure wet-cut audit lost its config authority")
        if audit.cut_master_cell_index != self.layout.cut_cell_index:
            raise RuntimeError("departure wet-cut audit names the wrong cut cell")
        partition = self.candidate_geometry.cells[self.layout.cut_cell_index]
        expected_fraction = (
            partition.wet_volume_m3
            / self.candidate_geometry.master_grid.volumes[self.layout.cut_cell_index]
        )
        if audit.wet_volume_fraction != expected_fraction:
            raise RuntimeError("departure wet-cut audit lost its candidate geometry")
        if (
            audit.transformed_cut_water_residual_mol_s
            != self.residuals.wet_water_mol_s[-1]
            or audit.transformed_cut_common_energy_residual_w
            != self.residuals.wet_energy_w[-1]
        ):
            raise RuntimeError("departure wet-cut audit lost its transformed rows")


def adapt_face_arrival_assembly(
    source: arrival.FaceArrivalAssembly,
) -> FaceDepartureState:
    """Adapt exact arrival event primitives to the narrow departure protocol.

    The arrival assembly remains the provenance object, but acceptance and
    nonlinear certification are deliberately owned by its future integrator.
    """

    if not isinstance(source, arrival.FaceArrivalAssembly):
        raise TypeError("departure adapter requires a FaceArrivalAssembly")
    return FaceDepartureState(
        geometry=source.event_geometry,
        config=source.before.config,
        time_s=source.before.time_s + source.candidate.event_time_s,
        wet_temperatures_k=source.candidate.wet_temperatures_k,
        wet_retained_water_loadings=(
            source.candidate.wet_retained_water_loadings
        ),
        dry_temperatures_k=source.candidate.dry_temperatures_k,
        dry_y_hexane=source.candidate.dry_y_hexane,
        historical_hexane_loadings=source.before.historical_hexane_loadings,
        oil_fraction_labels=source.before.oil_fraction_labels,
        source="adapted exact face-arrival residual state; acceptance external",
    )


def adapt_face_arrival_committed_state(source: object) -> FaceDepartureState:
    """Adapt only the safeguarded event integrator's certified face commit.

    The import is local so the residual foundation remains independent of the
    event solver at module-import time.  Full cumulative ledger fields stay on
    the immutable committed source and are consumed by a future departure
    integrator; this narrow state supplies only the constitutive before data
    needed by :func:`assemble_face_departure`.
    """

    from dtdc_simulator.core2.particle import cut_face_event_integrator as event_solver

    if not isinstance(source, event_solver.FaceArrivalCommittedState):
        raise TypeError(
            "certified departure adapter requires FaceArrivalCommittedState"
        )
    return FaceDepartureState(
        geometry=source.geometry,
        config=source.config,
        time_s=source.time_s,
        wet_temperatures_k=source.wet_temperatures_k,
        wet_retained_water_loadings=source.wet_retained_water_loadings,
        dry_temperatures_k=source.dry_temperatures_k,
        dry_y_hexane=source.dry_y_hexane,
        historical_hexane_loadings=source.historical_hexane_loadings,
        oil_fraction_labels=source.oil_fraction_labels,
        source=(
            "adapted safeguarded FaceArrivalCommittedState; cumulative ledgers "
            "remain on the immutable source commit"
        ),
    )


def layout_for_departure(before: FaceDepartureState) -> FaceDepartureLayout:
    """Return the unique strict cut chart immediately inside face ``k``."""

    if not isinstance(before, FaceDepartureState):
        raise TypeError("departure requires a FaceDepartureState")
    k = before.departure_face_index
    n = before.geometry.master_grid.n
    layout = FaceDepartureLayout(
        wet_cell_indices=tuple(range(k)),
        dry_cell_indices=tuple(range(k - 1, n)),
        cut_cell_index=k - 1,
        departure_face_index=k,
        master_cell_count=n,
    )
    if not layout.is_square:
        raise RuntimeError("face-departure rank formula is internally inconsistent")
    return layout


def assemble_face_departure(
    before: FaceDepartureState,
    candidate: FaceDepartureUnknowns,
    dt_s: float,
    surface_boundary: ct.PoreBoundary,
    newborn_trace_gradients: NewbornDryTraceGradients | None = None,
    *,
    enforce_reduced_film_thresholds: bool = True,
) -> FaceDepartureAssembly:
    """Assemble one strict fixed-dt departure residual or roll back exactly."""

    try:
        return _assemble_face_departure(
            before,
            candidate,
            dt_s,
            surface_boundary,
            newborn_trace_gradients,
            enforce_reduced_film_thresholds,
        )
    except FaceDepartureStepError:
        raise
    except Exception as exc:
        if not isinstance(before, FaceDepartureState):
            raise
        raise FaceDepartureStepError(
            f"inner-face departure residual rejected with exact rollback: {exc}",
            before,
        ) from exc


def _assemble_face_departure(
    before: FaceDepartureState,
    candidate: FaceDepartureUnknowns,
    dt_s: float,
    surface_boundary: ct.PoreBoundary,
    newborn_trace_gradients: NewbornDryTraceGradients | None,
    enforce_reduced_film_thresholds: bool,
) -> FaceDepartureAssembly:
    if not isinstance(before, FaceDepartureState):
        raise TypeError("before must be a FaceDepartureState")
    if not isinstance(candidate, FaceDepartureUnknowns):
        raise TypeError("candidate must be FaceDepartureUnknowns")
    if not isinstance(
        surface_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("surface boundary must be a supported pore boundary")
    if not isinstance(enforce_reduced_film_thresholds, bool):
        raise TypeError("reduced-film threshold enforcement flag must be boolean")
    if newborn_trace_gradients is not None and not isinstance(
        newborn_trace_gradients,
        NewbornDryTraceGradients,
    ):
        raise TypeError("newborn trace gradients have the wrong type")
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError("departure substep duration must be positive and finite")
    layout = layout_for_departure(before)
    if len(candidate.vector()) != layout.unknown_count:
        raise ValueError("departure candidate does not match its strict-cut rank")
    if not all(math.isfinite(value) for value in candidate.vector()):
        raise ValueError("departure candidate values must remain finite")
    _require_length(
        "departure wet temperatures",
        candidate.wet_temperatures_k,
        layout.wet_piece_count,
    )
    _require_length(
        "departure wet retained water",
        candidate.wet_retained_water_loadings,
        layout.wet_piece_count,
    )
    _require_length(
        "departure dry temperatures",
        candidate.dry_temperatures_k,
        layout.dry_piece_count,
    )
    _require_length(
        "departure dry compositions",
        candidate.dry_y_hexane,
        layout.dry_piece_count,
    )
    _require_length(
        "departure dry Stefan fluxes",
        candidate.dry_total_stefan_fluxes_mol_m2_s,
        layout.dry_face_count,
    )
    cut._validate_dry_boundary(surface_boundary, before.config)  # noqa: SLF001

    geometry = cg.partition_master_grid(
        before.geometry.master_grid,
        z=candidate.front_z,
    )
    cut_layout = cut.layout_for_geometry(geometry)
    if (
        cut_layout.wet_cell_indices != layout.wet_cell_indices
        or cut_layout.dry_cell_indices != layout.dry_cell_indices
        or cut_layout.cut_cell_index != layout.cut_cell_index
    ):
        raise FaceDepartureTopologyError(
            "departure must remain a strict cut in cell k-1; face crossing rejected"
        )
    if geometry.front.radius_m >= before.geometry.front.radius_m:
        raise FaceDepartureTopologyError(
            "departure requires strict inward recession, not an interior stall"
        )
    swept = cg.swept_cut_geometry(before.geometry, geometry, dt_s)
    old_wet = _evaluate_old_wet_cells(before)
    old_dry = _evaluate_old_dry_cells(before)
    wet_cells = cut._evaluate_wet_piece_states(  # noqa: SLF001
        cut_layout,
        candidate.wet_temperatures_k,
        candidate.wet_retained_water_loadings,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
        before.config,
    )
    dry_cells = cut._evaluate_dry_piece_states(  # noqa: SLF001
        cut_layout,
        candidate.dry_temperatures_k,
        candidate.dry_y_hexane,
        before.oil_fraction_labels,
        before.config,
    )
    interface = cut.evaluate_interface_state(
        candidate.interface_temperature_k,
        before.config,
        before.historical_hexane_loadings[layout.cut_cell_index],
        before.oil_fraction_labels[layout.cut_cell_index],
    )
    wet_fluxes = cut._wet_face_fluxes(  # noqa: SLF001
        geometry,
        cut_layout,
        wet_cells,
        interface,
        before.config,
    )
    dry_fluxes = cut._dry_face_fluxes(  # noqa: SLF001
        geometry,
        cut_layout,
        dry_cells,
        interface,
        candidate.dry_total_stefan_fluxes_mol_m2_s,
        surface_boundary,
        before.oil_fraction_labels,
        before.config,
        newborn_temperature_gradient_k_m=(
            None
            if newborn_trace_gradients is None
            else newborn_trace_gradients.temperature_gradient_k_m
        ),
        newborn_composition_gradient_m_inv=(
            None
            if newborn_trace_gradients is None
            else newborn_trace_gradients.composition_gradient_m_inv
        ),
        enforce_reduced_film_thresholds=enforce_reduced_film_thresholds,
    )
    changes = _inventory_changes(
        before,
        geometry,
        cut_layout,
        old_wet,
        old_dry,
        wet_cells,
        dry_cells,
        swept,
    )
    residuals, wet_hexane_identity = _residual_blocks(
        before.config,
        geometry,
        cut_layout,
        changes,
        interface,
        wet_fluxes,
        dry_fluxes,
        swept,
        dt_s,
        # O9a front-donor correction: the piece the departing front is
        # consuming is the OUTERMOST wet piece (the new cut cell's wet part),
        # and the donor is that piece's own bulk retained-water concentration.
        #
        # THE TIME LEVEL IS THE MODULE'S OWN, NOT A NEW CHOICE.
        # ``_inventory_changes`` books the matching geometric storage term as
        # ``Delta(V)*C_old`` in ``stable_wet_change``; the ALE sweep is that
        # identical swept volume seen from the flux side, so the donor is read
        # at the same time level - ``C_old`` of the consumed piece.  The two
        # then cancel term by term and the row keeps exactly the ``V_new``
        # storage form the departure decomposition is written in.
        old_wet[-1].retained_water_concentration_mol_m3,
        # O10a energy-donor completion (ruled 2026-08-21,
        # docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md).  The
        # Rankine-Hugoniot jump conditions must convect mass and enthalpy OF
        # THE SAME MATERIAL STATE - an internal-consistency requirement of the
        # balance laws - and that state is the bulk, undrained one the O9
        # traverse measured (Peclet 3.4e3).  So the enthalpy donor is the SAME
        # consumed piece read at the SAME time level as its mass donor one line
        # above: ``_inventory_changes`` books the paired geometric energy
        # storage term as ``Delta(V)*rho_e_old`` from exactly this state.
        old_wet[-1].energy_density_j_m3,
    )
    # A departure candidate is always a strict positive cut in cell k-1.  Apply
    # the same prospectively selected Amendment-12 finite-volume operator used
    # by ordinary same-cell transport; its geometric weight is identically zero
    # outside the declared vanishing-wet-cut band.  The exact zero-volume face
    # state itself remains exclusively owned by the face-arrival transaction.
    residuals, wet_cut_treatment_audit = (
        cut._apply_vanishing_wet_cut_treatment(  # noqa: SLF001
            geometry,
            cut_layout,
            wet_cells,
            interface,
            wet_fluxes,
            residuals,
            before.config,
        )
    )
    ledger = _ledger(
        before,
        geometry,
        changes,
        residuals,
        wet_hexane_identity,
        dry_fluxes[-1],
        swept,
        dt_s,
        (*wet_fluxes, *dry_fluxes),
        dry_cells[0],
    )
    return FaceDepartureAssembly(
        before=before,
        candidate=candidate,
        layout=layout,
        candidate_geometry=geometry,
        swept_geometry=swept,
        old_wet_cells=old_wet,
        old_dry_cells=old_dry,
        wet_cells=wet_cells,
        dry_cells=dry_cells,
        interface=interface,
        wet_face_fluxes=wet_fluxes,
        dry_face_fluxes=dry_fluxes,
        inventory_changes=changes,
        vanishing_wet_cut_treatment_audit=wet_cut_treatment_audit,
        residuals=residuals,
        ledger=ledger,
        surface_boundary=surface_boundary,
        dt_s=dt_s,
    )


def _evaluate_old_wet_cells(
    state: FaceDepartureState,
) -> tuple[ww.WetWaterCellState, ...]:
    k = state.departure_face_index
    return tuple(
        ww.evaluate_cell(T, W, Xh, wo, state.config.wet)
        for T, W, Xh, wo in zip(
            state.wet_temperatures_k,
            state.wet_retained_water_loadings,
            state.historical_hexane_loadings[:k],
            state.oil_fraction_labels[:k],
        )
    )


def _evaluate_old_dry_cells(
    state: FaceDepartureState,
) -> tuple[cp.EquilibriumPoreState, ...]:
    k = state.departure_face_index
    cells = []
    for offset, (temperature, y_hexane) in enumerate(
        zip(state.dry_temperatures_k, state.dry_y_hexane)
    ):
        pore = replace(
            state.config.dry.pore,
            w_o=state.oil_fraction_labels[k + offset],
        )
        cut._validate_open_dry_band(  # noqa: SLF001
            temperature,
            y_hexane,
            state.config,
            pore=pore,
        )
        cells.append(
            cp.evaluate_equilibrium(
                temperature,
                state.config.dry.pressure_pa,
                y_hexane,
                pore,
            )
        )
    return tuple(cells)


def _inventory_changes(
    before: FaceDepartureState,
    geometry: cg.CutGeometry,
    layout: cut.CutTransportLayout,
    old_wet: Sequence[ww.WetWaterCellState],
    old_dry: Sequence[cp.EquilibriumPoreState],
    new_wet: Sequence[ww.WetWaterCellState],
    new_dry: Sequence[cp.EquilibriumPoreState],
    swept: cg.SweptCutGeometry,
) -> FaceDepartureInventoryChanges:
    wet_volumes = tuple(
        geometry.cells[index].wet_volume_m3 for index in layout.wet_cell_indices
    )
    wet_volume_changes = tuple(
        swept.cell_wet_volume_changes_m3[index]
        for index in layout.wet_cell_indices
    )

    def stable_wet_change(
        old_values: Sequence[float],
        new_values: Sequence[float],
    ) -> tuple[float, ...]:
        return tuple(
            math.fsum((volume * (new - old), volume_change * old))
            for volume, volume_change, old, new in zip(
                wet_volumes,
                wet_volume_changes,
                old_values,
                new_values,
            )
        )

    wet_hexane_concentrations = tuple(
        before.config.wet.wet.rho_dm_p
        * before.historical_hexane_loadings[index]
        / hx.M
        for index in layout.wet_cell_indices
    )
    wet_water = stable_wet_change(
        tuple(cell.retained_water_concentration_mol_m3 for cell in old_wet),
        tuple(cell.retained_water_concentration_mol_m3 for cell in new_wet),
    )
    wet_hexane = stable_wet_change(
        wet_hexane_concentrations,
        wet_hexane_concentrations,
    )
    wet_energy = stable_wet_change(
        tuple(cell.energy_density_j_m3 for cell in old_wet),
        tuple(cell.energy_density_j_m3 for cell in new_wet),
    )

    newborn_index = layout.cut_cell_index
    newborn_volume = geometry.cells[newborn_index].dry_volume_m3
    newborn_swept = swept.cell_dry_volume_changes_m3[newborn_index]
    tolerance = 64.0 * math.ulp(geometry.master_grid.total_volume)
    if abs(newborn_volume - newborn_swept) > tolerance:
        raise RuntimeError("newborn departure subpiece disagrees with exact ALE sweep")
    first_water, first_hexane, first_energy = birth.zero_old_dry_inventories(
        newborn_volume,
        new_dry[0],
    )
    dry_water = [first_water]
    dry_hexane = [first_hexane]
    dry_energy = [first_energy]
    for cell_index, old_cell, new_cell in zip(
        layout.dry_cell_indices[1:],
        old_dry,
        new_dry[1:],
    ):
        volume = geometry.master_grid.volumes[cell_index]
        volume_change = swept.cell_dry_volume_changes_m3[cell_index]
        if abs(volume_change) > tolerance:
            raise RuntimeError("an existing outer dry cell changed geometric volume")
        dry_water.append(
            math.fsum(
                (
                    volume
                    * (
                        new_cell.total_water_concentration_mol_m3
                        - old_cell.total_water_concentration_mol_m3
                    ),
                    volume_change * old_cell.total_water_concentration_mol_m3,
                )
            )
        )
        dry_hexane.append(
            math.fsum(
                (
                    volume
                    * (
                        new_cell.total_hexane_concentration_mol_m3
                        - old_cell.total_hexane_concentration_mol_m3
                    ),
                    volume_change * old_cell.total_hexane_concentration_mol_m3,
                )
            )
        )
        dry_energy.append(
            math.fsum(
                (
                    volume * (new_cell.energy_density_j_m3 - old_cell.energy_density_j_m3),
                    volume_change * old_cell.energy_density_j_m3,
                )
            )
        )
    return FaceDepartureInventoryChanges(
        wet_water_mol=wet_water,
        wet_hexane_mol=wet_hexane,
        wet_energy_j=wet_energy,
        dry_water_mol=tuple(dry_water),
        dry_hexane_mol=tuple(dry_hexane),
        dry_energy_j=tuple(dry_energy),
    )


def _residual_blocks(
    config: cut.CutTransportConfig,
    geometry: cg.CutGeometry,
    layout: cut.CutTransportLayout,
    changes: FaceDepartureInventoryChanges,
    interface: cut.CutInterfaceState,
    wet_fluxes: Sequence[ww.WetWaterFaceFlux],
    dry_fluxes: Sequence[cut.DryFaceFlux],
    swept: cg.SweptCutGeometry,
    dt_s: float,
    swept_wet_water_donor_concentration_mol_m3: float,
    swept_wet_energy_donor_density_j_m3: float,
) -> tuple[cut.CutResidualBlocks, float]:
    grid = geometry.master_grid
    area_gamma = 4.0 * math.pi * geometry.front.radius_m**2
    q_gamma = swept.interface_swept_volume_rate_m3_s

    # O9a ALE FRONT-DONOR CORRECTION, ruled 2026-08-21 in
    # docs/GT_PS2_ALE_DONOR_CORRECTION_RULING_PACKET_2026-08-21.md.  Carried to
    # this sibling chart in Stage 1b for formulation consistency: the departure
    # chart books the same ALE sweep as ``cut_transport._residual_blocks``, so
    # it must book the same donor.  The receding front consumes wet material at
    # that material's OWN BULK LOADING, and the sweep is therefore debited at
    # the consumed piece's bulk retained-water concentration rather than at the
    # zero-volume equilibrium interface trace.
    #
    # REMOVED ASSUMPTION (a removal, not an addition): the trace donor
    # implicitly asserted that retained water arriving at the front has
    # PRE-DRAINED to the dry-pore equilibrium loading before conversion, which
    # the traverse's measured Peclet 3.4e3 forbids.
    #
    # ROLE SEPARATION - only the MASS-BOOKKEEPING donor moves.  The interface
    # trace keeps its potential/flux role in full: ``cut._wet_face_fluxes``
    # still evaluates the front-face diffusive flux against
    # ``interface.wet.retained_water_loading`` and the trace KKT dual, and the
    # hexane and energy rows are untouched.
    #
    # The released excess leaves through the UNCHANGED Rankine-Hugoniot balance
    # below: ``rh_water`` is the signed expansion of ``dry_h - wet_h``, so its
    # ``C_wet`` is this same donor by construction and the jump row keeps its
    # existing structure.
    wet_water_donor_concentration = swept_wet_water_donor_concentration_mol_m3
    # O10a ENERGY-DONOR COMPLETION, ruled 2026-08-21 in
    # docs/GT_PS2_O10_DONOR_COMPLETION_RULING_PACKET_2026-08-21.md and carried
    # to this sibling chart in the same family pass.
    #
    # PHYSICS GROUND, as recorded in that ruling: two-phase moving-boundary
    # (Rankine-Hugoniot) jump conditions must convect mass and enthalpy OF THE
    # SAME MATERIAL STATE.  That is an internal-consistency requirement of the
    # balance laws themselves, not an empirical claim; WHICH material state it
    # is - the bulk, undrained one - is the O9-measured fact (Peclet 3.4e3
    # forbids pre-drainage).  With the mass donor on the consumed piece's bulk
    # loading, an energy donor left on the zero-volume equilibrium trace would
    # have the same swept material carry water mass at the bulk loading and
    # water enthalpy at the trace loading.
    #
    # The donor is the consumed piece's OWN bulk energy density at the SAME
    # time level as its mass donor, which is the level ``_inventory_changes``
    # books the paired geometric energy storage term ``Delta(V)*rho_e_old`` at
    # in ``stable_wet_change``, so the energy mis-booking term cancels term by
    # term exactly as the mass one does.  The trace keeps its potential/flux
    # role in full on the energy side too.
    wet_energy_donor_density = swept_wet_energy_donor_density_j_m3
    wet_water_h = math.fsum(
        (
            area_gamma * wet_fluxes[-1].retained_water_flux_mol_m2_s,
            -q_gamma * wet_water_donor_concentration,
        )
    )
    wet_hexane_concentration = (
        config.wet.wet.rho_dm_p * interface.historical_hexane_loading / hx.M
    )
    wet_hexane_h = -q_gamma * wet_hexane_concentration
    wet_energy_h = math.fsum(
        (
            area_gamma * wet_fluxes[-1].total_energy_flux_w_m2,
            # O10a: the swept material's enthalpy travels at the same (bulk,
            # old-time) state as its mass; see the block above.
            -q_gamma * wet_energy_donor_density,
        )
    )
    dry_component = dry_fluxes[0].component
    dry_water_h = math.fsum(
        (
            area_gamma * dry_component.conserved_water_flux_mol_m2_s,
            -q_gamma * interface.dry.total_water_concentration_mol_m3,
        )
    )
    dry_hexane_h = math.fsum(
        (
            area_gamma * dry_component.conserved_hexane_flux_mol_m2_s,
            -q_gamma * interface.dry.total_hexane_concentration_mol_m3,
        )
    )
    dry_energy_h = math.fsum(
        (
            area_gamma * dry_fluxes[0].energy.total_energy_flux_w_m2,
            -q_gamma * interface.dry.energy_density_j_m3,
        )
    )
    rh_water = math.fsum(
        (
            area_gamma * dry_component.conserved_water_flux_mol_m2_s,
            -area_gamma * wet_fluxes[-1].retained_water_flux_mol_m2_s,
            -q_gamma * interface.dry.total_water_concentration_mol_m3,
            # O9a: the same corrected wet-side donor as ``wet_water_h`` above,
            # so this row stays exactly ``dry_h - wet_h`` and the global water
            # telescoping identity is untouched.
            q_gamma * wet_water_donor_concentration,
        )
    )
    rh_hexane = math.fsum(
        (
            area_gamma * dry_component.conserved_hexane_flux_mol_m2_s,
            -q_gamma * interface.dry.total_hexane_concentration_mol_m3,
            q_gamma * wet_hexane_concentration,
        )
    )
    rh_energy = math.fsum(
        (
            area_gamma * dry_fluxes[0].energy.total_energy_flux_w_m2,
            -area_gamma * wet_fluxes[-1].total_energy_flux_w_m2,
            -q_gamma * interface.dry.energy_density_j_m3,
            # O10a: the same corrected wet-side energy donor as
            # ``wet_energy_h`` above, so this row stays exactly
            # ``dry_h - wet_h`` and the global energy telescoping identity is
            # untouched.  Because this donor's water content is now the SAME
            # bulk concentration ``rh_water`` convects, the water caloric datum
            # is again an exact gauge freedom of this row.
            q_gamma * wet_energy_donor_density,
        )
    )

    wet_water_rates = [0.0]
    wet_energy_rates = [0.0]
    for face_position, cell_index in enumerate(
        layout.wet_cell_indices[1:],
        start=1,
    ):
        area = grid.areas[cell_index]
        wet_water_rates.append(
            area * wet_fluxes[face_position].retained_water_flux_mol_m2_s
        )
        wet_energy_rates.append(
            area * wet_fluxes[face_position].total_energy_flux_w_m2
        )
    wet_water_rates.append(wet_water_h)
    wet_energy_rates.append(wet_energy_h)
    wet_water_residuals = tuple(
        math.fsum(
            (
                change / dt_s,
                wet_water_rates[index + 1],
                -wet_water_rates[index],
            )
        )
        for index, change in enumerate(changes.wet_water_mol)
    )
    wet_energy_residuals = tuple(
        math.fsum(
            (
                change / dt_s,
                wet_energy_rates[index + 1],
                -wet_energy_rates[index],
            )
        )
        for index, change in enumerate(changes.wet_energy_j)
    )

    dry_water_rates = [dry_water_h]
    dry_hexane_rates = [dry_hexane_h]
    dry_energy_rates = [dry_energy_h]
    for face_position in range(1, layout.dry_face_count):
        area = (
            grid.areas[-1]
            if face_position == layout.dry_piece_count
            else grid.areas[layout.dry_cell_indices[face_position]]
        )
        flux = dry_fluxes[face_position]
        dry_water_rates.append(
            area * flux.component.conserved_water_flux_mol_m2_s
        )
        dry_hexane_rates.append(
            area * flux.component.conserved_hexane_flux_mol_m2_s
        )
        dry_energy_rates.append(area * flux.energy.total_energy_flux_w_m2)

    def regional(
        inventory_changes: Sequence[float],
        rates: Sequence[float],
    ) -> tuple[float, ...]:
        return tuple(
            math.fsum(
                (
                    change / dt_s,
                    rates[index + 1],
                    -rates[index],
                )
            )
            for index, change in enumerate(inventory_changes)
        )

    wet_hexane_identity = math.fsum(
        (math.fsum(changes.wet_hexane_mol) / dt_s, wet_hexane_h)
    )
    return (
        cut.CutResidualBlocks(
            wet_water_mol_s=wet_water_residuals,
            wet_energy_w=wet_energy_residuals,
            dry_water_mol_s=regional(changes.dry_water_mol, dry_water_rates),
            dry_hexane_mol_s=regional(changes.dry_hexane_mol, dry_hexane_rates),
            dry_energy_w=regional(changes.dry_energy_j, dry_energy_rates),
            rh_water_mol_s=rh_water,
            rh_hexane_mol_s=rh_hexane,
            rh_energy_w=rh_energy,
        ),
        wet_hexane_identity,
    )


def _ledger(
    before: FaceDepartureState,
    geometry: cg.CutGeometry,
    changes: FaceDepartureInventoryChanges,
    residuals: cut.CutResidualBlocks,
    wet_hexane_identity: float,
    surface: cut.DryFaceFlux,
    swept: cg.SweptCutGeometry,
    dt_s: float,
    all_fluxes: Sequence[object],
    newborn_dry: cp.EquilibriumPoreState,
) -> FaceDepartureLedger:
    water_blocks = math.fsum(
        (
            *residuals.wet_water_mol_s,
            *residuals.dry_water_mol_s,
            residuals.rh_water_mol_s,
        )
    )
    hexane_blocks = math.fsum(
        (
            wet_hexane_identity,
            *residuals.dry_hexane_mol_s,
            residuals.rh_hexane_mol_s,
        )
    )
    energy_blocks = math.fsum(
        (
            *residuals.wet_energy_w,
            *residuals.dry_energy_w,
            residuals.rh_energy_w,
        )
    )
    surface_area = geometry.master_grid.areas[-1]
    water_surface = surface_area * surface.component.conserved_water_flux_mol_m2_s
    hexane_surface = surface_area * surface.component.conserved_hexane_flux_mol_m2_s
    energy_surface = surface_area * surface.energy.total_energy_flux_w_m2
    water_direct = math.fsum(
        (
            math.fsum((*changes.wet_water_mol, *changes.dry_water_mol)) / dt_s,
            water_surface,
        )
    )
    hexane_direct = math.fsum(
        (
            math.fsum((*changes.wet_hexane_mol, *changes.dry_hexane_mol)) / dt_s,
            hexane_surface,
        )
    )
    energy_direct = math.fsum(
        (
            math.fsum((*changes.wet_energy_j, *changes.dry_energy_j)) / dt_s,
            energy_surface,
        )
    )
    water_error = water_blocks - water_direct
    hexane_error = hexane_blocks - hexane_direct
    energy_error = energy_blocks - energy_direct
    # Use transfer/inventory-change scales, not the nearly-zero residual itself.
    water_scale = max(
        math.fsum(abs(value) for value in (*changes.wet_water_mol, *changes.dry_water_mol))
        / dt_s,
        abs(water_surface),
        *(abs(value) for value in residuals.wet_water_mol_s),
        *(abs(value) for value in residuals.dry_water_mol_s),
        abs(residuals.rh_water_mol_s),
        1.0e-300,
    )
    hexane_scale = max(
        math.fsum(abs(value) for value in (*changes.wet_hexane_mol, *changes.dry_hexane_mol))
        / dt_s,
        abs(hexane_surface),
        abs(wet_hexane_identity),
        *(abs(value) for value in residuals.dry_hexane_mol_s),
        abs(residuals.rh_hexane_mol_s),
        1.0e-300,
    )
    energy_scale = max(
        math.fsum(abs(value) for value in (*changes.wet_energy_j, *changes.dry_energy_j))
        / dt_s,
        abs(energy_surface),
        *(abs(value) for value in residuals.wet_energy_w),
        *(abs(value) for value in residuals.dry_energy_w),
        abs(residuals.rh_energy_w),
        1.0e-300,
    )
    entropy_values: list[float] = []
    for flux in all_fluxes:
        if isinstance(flux, ww.WetWaterFaceFlux):
            entropy_values.append(flux.total_entropy_production_w_m3_k)
        elif isinstance(flux, cut.DryFaceFlux):
            entropy_values.append(flux.total_entropy_production_w_m3_k)
    newborn_index = before.departure_face_index - 1
    newborn_volume = geometry.cells[newborn_index].dry_volume_m3
    newborn_swept = swept.cell_dry_volume_changes_m3[newborn_index]
    existing_volume_changes = tuple(
        swept.cell_dry_volume_changes_m3[index]
        for index in range(before.departure_face_index, geometry.master_grid.n)
    )
    return FaceDepartureLedger(
        wet_hexane_sweep_identity_mol_s=wet_hexane_identity,
        water_global_from_blocks_mol_s=water_blocks,
        water_global_direct_mol_s=water_direct,
        water_telescoping_error_mol_s=water_error,
        hexane_global_from_blocks_mol_s=hexane_blocks,
        hexane_global_direct_mol_s=hexane_direct,
        hexane_telescoping_error_mol_s=hexane_error,
        energy_global_from_blocks_w=energy_blocks,
        energy_global_direct_w=energy_direct,
        energy_telescoping_error_w=energy_error,
        water_telescoping_relative_error=abs(water_error) / water_scale,
        hexane_telescoping_relative_error=abs(hexane_error) / hexane_scale,
        energy_telescoping_relative_error=abs(energy_error) / energy_scale,
        newborn_old_dry_volume_m3=0.0,
        newborn_dry_volume_m3=newborn_volume,
        newborn_swept_volume_m3=newborn_swept,
        newborn_water_inventory_identity_mol=math.fsum(
            (
                changes.dry_water_mol[0],
                -newborn_volume * newborn_dry.total_water_concentration_mol_m3,
            )
        ),
        newborn_hexane_inventory_identity_mol=math.fsum(
            (
                changes.dry_hexane_mol[0],
                -newborn_volume * newborn_dry.total_hexane_concentration_mol_m3,
            )
        ),
        newborn_energy_inventory_identity_j=math.fsum(
            (
                changes.dry_energy_j[0],
                -newborn_volume * newborn_dry.energy_density_j_m3,
            )
        ),
        maximum_existing_dry_volume_change_m3=max(
            (abs(value) for value in existing_volume_changes),
            default=0.0,
        ),
        historical_wet_hexane_label_change=0.0,
        geometry_cell_gcl_residual_m3=swept.max_cell_gcl_residual_m3,
        minimum_entropy_production_w_m3_k=min(entropy_values, default=0.0),
    )


def datum_covariant_residual_vector(
    residuals: cut.CutResidualBlocks,
) -> tuple[float, ...]:
    """Reuse the approved equation-equivalent common-datum row basis."""

    return birth.datum_covariant_residual_vector(residuals)


def _require_length(name: str, values: Sequence[float], expected: int) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} has length {len(values)}, expected {expected}")


__all__ = [
    "FaceDepartureAssembly",
    "FaceDepartureInventoryChanges",
    "FaceDepartureLayout",
    "FaceDepartureLedger",
    "FaceDepartureState",
    "FaceDepartureStepError",
    "FaceDepartureTopologyError",
    "FaceDepartureUnknowns",
    "NewbornDryTraceGradients",
    "adapt_face_arrival_assembly",
    "adapt_face_arrival_committed_state",
    "assemble_face_departure",
    "datum_covariant_residual_vector",
    "layout_for_departure",
]
