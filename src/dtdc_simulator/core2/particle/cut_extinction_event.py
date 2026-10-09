r"""Exact ``s=0`` endpoint residual and rank finding for Gate 1g.

The incoming state is a strict cut in material cell zero.  The endpoint is the
exact fully dry sphere: every stationary material cell is dry, the center face
has zero area and exact symmetry flux, and neither a wet epsilon-piece nor an
interface temperature is retained.

This endpoint is deliberately *not* presented as a solvable event integrator.
The frozen endpoint authority says that the wet-domain and interface equations
are absent at ``s=0``.  After eliminating the zero-area center Stefan flux, an
unknown-time endpoint candidate contains

``2*N dry primitives + N positive-area Stefan fluxes + event time = 3*N + 1``

unknowns, while water, n-hexane, and common-datum energy provide ``3*N``
regional balances.  The fixed-time projection is square, but replacing the
front coordinate by an unknown event time leaves one genuine right-null
degree of freedom.  Retaining the strict-cut ``T_Gamma`` and center-face
Stefan flux does not repair this: the former belongs to endpoint-absent
interface equations and the latter has an identically zero area column.

Consequently this module assembles and audits the conservative endpoint but
refuses to select or commit one member of the one-parameter family.  The event
time must instead be supplied by a separately justified integration and root
localisation of the positive-core chart.  Treating the finite-step swept rate
as a local ``s=0`` Rankine--Hugoniot rate would be an invented extinction
condition, because the physical endpoint interface is absent.

All numerical results remain non-qualifying coefficient oracles.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Sequence

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.props import hexane as hx


class ExtinctionTopologyError(ValueError):
    """The proposal is not the exact cell-zero wet-core extinction seam."""


class ExtinctionAuthorityRequired(RuntimeError):
    """A commit was requested despite the documented unknown-time rank loss."""


class ExtinctionStepError(RuntimeError):
    """Rejected endpoint evaluation carrying the exact immutable rollback."""

    def __init__(
        self,
        message: str,
        rollback_state: cut.CutTransportState,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state


@dataclass(frozen=True)
class ExtinctionLayout:
    """Endpoint rank after exact center symmetry and interface elimination."""

    dry_cell_indices: tuple[int, ...]
    positive_area_face_indices: tuple[int, ...]
    master_cell_count: int
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        n = self.master_cell_count
        if n < 1:
            raise ValueError("extinction needs at least one material cell")
        if self.dry_cell_indices != tuple(range(n)):
            raise ValueError("the extinction endpoint must make every cell dry")
        if self.positive_area_face_indices != tuple(range(1, n + 1)):
            raise ValueError("only positive-area dry faces may carry unknown fluxes")

    @property
    def dry_cell_count(self) -> int:
        return self.master_cell_count

    @property
    def positive_area_face_count(self) -> int:
        return self.master_cell_count

    @property
    def unknown_blocks(self) -> tuple[tuple[str, int], ...]:
        return (
            ("dry_temperature", self.dry_cell_count),
            ("dry_y_hexane", self.dry_cell_count),
            ("positive_area_total_stefan_flux", self.positive_area_face_count),
            ("event_time", 1),
        )

    @property
    def residual_blocks(self) -> tuple[tuple[str, int], ...]:
        return (
            ("dry_water", self.dry_cell_count),
            ("dry_hexane", self.dry_cell_count),
            ("dry_energy", self.dry_cell_count),
        )

    @property
    def unknown_count(self) -> int:
        return sum(size for _, size in self.unknown_blocks)

    @property
    def residual_count(self) -> int:
        return sum(size for _, size in self.residual_blocks)

    @property
    def fixed_time_unknown_count(self) -> int:
        return self.unknown_count - 1

    @property
    def right_nullity_lower_bound(self) -> int:
        return self.unknown_count - self.residual_count

    @property
    def is_unknown_time_square(self) -> bool:
        return self.unknown_count == self.residual_count

    @property
    def is_fixed_time_square(self) -> bool:
        return self.fixed_time_unknown_count == self.residual_count


@dataclass(frozen=True)
class ExtinctionRankAudit:
    """Machine-readable authority and structural-rank certificate."""

    layout: ExtinctionLayout
    center_face_area_m2: float
    endpoint_interface_equation_count: int = field(default=0, init=False)
    endpoint_interface_temperature_authorized: bool = field(
        default=False,
        init=False,
    )
    center_total_stefan_flux_is_unknown: bool = field(default=False, init=False)
    naive_continuation_center_flux_zero_columns: int = field(default=1, init=False)
    safeguarded_unknown_time_integrator_authorized: bool = field(
        default=False,
        init=False,
    )
    finding_code: str = field(default="GT-PS-2-G1G-EXT-RANK-01", init=False)
    finding: str = field(
        default=(
            "At s=0 the wet domain and interface equations are absent and the "
            "center flux is fixed by exact symmetry. Event time therefore leaves "
            "one unresolved scalar degree of freedom. Integrate/localize the "
            "positive-core chart; do not invent T_Gamma, an epsilon core, or an "
            "extinction condition."
        ),
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.center_face_area_m2 != 0.0:
            raise ValueError("the spherical center face must have exactly zero area")
        if self.layout.right_nullity_lower_bound != 1:
            raise RuntimeError("extinction endpoint rank deficit changed unexpectedly")

    @property
    def unknown_time_rank_deficit(self) -> int:
        return self.layout.right_nullity_lower_bound


@dataclass(frozen=True)
class ExtinctionUnknowns:
    """Fully dry primitives, positive-area fluxes, and unknown event time."""

    dry_temperatures_k: tuple[float, ...]
    dry_y_hexane: tuple[float, ...]
    positive_area_total_stefan_fluxes_mol_m2_s: tuple[float, ...]
    event_time_s: float
    physically_qualifying: bool = field(default=False, init=False)

    def vector(self) -> tuple[float, ...]:
        return (
            *self.dry_temperatures_k,
            *self.dry_y_hexane,
            *self.positive_area_total_stefan_fluxes_mol_m2_s,
            self.event_time_s,
        )

    @classmethod
    def from_vector(
        cls,
        layout: ExtinctionLayout,
        values: Sequence[float],
    ) -> "ExtinctionUnknowns":
        vector = tuple(values)
        if len(vector) != layout.unknown_count:
            raise ValueError(
                f"extinction vector has length {len(vector)}, expected {layout.unknown_count}"
            )
        if not all(math.isfinite(value) for value in vector):
            raise ValueError("extinction candidate vector must remain finite")
        n = layout.dry_cell_count
        return cls(
            dry_temperatures_k=vector[:n],
            dry_y_hexane=vector[n : 2 * n],
            positive_area_total_stefan_fluxes_mol_m2_s=vector[2 * n : 3 * n],
            event_time_s=vector[-1],
        )


@dataclass(frozen=True)
class ExtinctionMaterialInventory:
    """Per-material-cell component and common-datum energy inventory."""

    water_cell_mol: tuple[float, ...]
    hexane_cell_mol: tuple[float, ...]
    energy_cell_j: tuple[float, ...]
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def total_water_mol(self) -> float:
        return math.fsum(self.water_cell_mol)

    @property
    def total_hexane_mol(self) -> float:
        return math.fsum(self.hexane_cell_mol)

    @property
    def total_energy_j(self) -> float:
        return math.fsum(self.energy_cell_j)


@dataclass(frozen=True)
class ExtinctionInventoryChanges:
    """Cancellation-safe changes preserving both old cell-zero subpieces."""

    water_cell_mol: tuple[float, ...]
    hexane_cell_mol: tuple[float, ...]
    energy_cell_j: tuple[float, ...]
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class ExtinctionResidualBlocks:
    """Three full-dry regional balances; no endpoint interface rows."""

    dry_water_mol_s: tuple[float, ...]
    dry_hexane_mol_s: tuple[float, ...]
    dry_energy_w: tuple[float, ...]
    dry_datum_covariant_energy_w: tuple[float, ...]
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def vector(self) -> tuple[float, ...]:
        return (
            *self.dry_water_mol_s,
            *self.dry_hexane_mol_s,
            *self.dry_energy_w,
        )

    @property
    def datum_covariant_vector(self) -> tuple[float, ...]:
        return (
            *self.dry_water_mol_s,
            *self.dry_hexane_mol_s,
            *self.dry_datum_covariant_energy_w,
        )


@dataclass(frozen=True)
class ExtinctionLedger:
    """Endpoint telescope, GCL, symmetry, and rank diagnostics."""

    water_global_from_regions_mol_s: float
    water_global_direct_mol_s: float
    water_telescoping_error_mol_s: float
    hexane_global_from_regions_mol_s: float
    hexane_global_direct_mol_s: float
    hexane_telescoping_error_mol_s: float
    energy_global_from_regions_w: float
    energy_global_direct_w: float
    energy_telescoping_error_w: float
    datum_covariant_energy_global_from_regions_w: float
    datum_covariant_energy_global_direct_w: float
    datum_covariant_energy_telescoping_error_w: float
    water_telescoping_relative_error: float
    hexane_telescoping_relative_error: float
    energy_telescoping_relative_error: float
    datum_covariant_energy_telescoping_relative_error: float
    geometry_cell_gcl_residual_m3: float
    center_face_area_m2: float
    center_total_stefan_flux_mol_m2_s: float
    minimum_entropy_production_w_m3_k: float
    event_time_rank_deficit: int
    endpoint_interface_equation_count: int
    historical_label_change: float
    oil_label_change: float
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def maximum_normalized_telescoping_error(self) -> float:
        return max(
            self.water_telescoping_relative_error,
            self.hexane_telescoping_relative_error,
            self.energy_telescoping_relative_error,
            self.datum_covariant_energy_telescoping_relative_error,
        )


@dataclass(frozen=True)
class ExtinctionAssembly:
    """Conservative endpoint evaluation with an explicit no-commit finding."""

    before: cut.CutTransportState
    candidate: ExtinctionUnknowns
    layout: ExtinctionLayout
    rank_audit: ExtinctionRankAudit
    endpoint_geometry: cg.CutGeometry
    swept_geometry: cg.SweptCutGeometry
    old_wet_cells: tuple[ww.WetWaterCellState, ...]
    old_dry_cells: tuple[cp.EquilibriumPoreState, ...]
    dry_cells: tuple[cp.EquilibriumPoreState, ...]
    dry_face_fluxes: tuple[cut.DryFaceFlux, ...]
    old_inventory: ExtinctionMaterialInventory
    current_inventory: ExtinctionMaterialInventory
    inventory_changes: ExtinctionInventoryChanges
    residuals: ExtinctionResidualBlocks
    ledger: ExtinctionLedger
    surface_boundary: ct.PoreBoundary
    historical_hexane_loadings: tuple[float, ...]
    oil_fraction_labels: tuple[float, ...]
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        n = self.layout.dry_cell_count
        if self.endpoint_geometry.front.regime != "fully_dry":
            raise RuntimeError("extinction assembly lost the exact dry endpoint")
        if self.endpoint_geometry.cut_cell_index is not None:
            raise RuntimeError("extinction endpoint retained an epsilon cut cell")
        if len(self.candidate.vector()) != self.layout.unknown_count:
            raise RuntimeError("extinction candidate does not match declared rank")
        if len(self.residuals.vector) != self.layout.residual_count:
            raise RuntimeError("extinction residual does not match declared rank")
        if len(self.dry_face_fluxes) != n + 1:
            raise RuntimeError("fully dry endpoint needs center through surface faces")
        for values in (
            self.historical_hexane_loadings,
            self.oil_fraction_labels,
            self.old_inventory.water_cell_mol,
            self.current_inventory.water_cell_mol,
            self.inventory_changes.water_cell_mol,
        ):
            if len(values) != n:
                raise RuntimeError("extinction material arrays lost grid alignment")

    @property
    def can_commit_to_fully_dry_transport(self) -> bool:
        return False

    def to_fully_dry_transport_state(self) -> ct.FullyDryTransportState:
        """Refuse an underdetermined commit without a localized event time."""

        raise ExtinctionAuthorityRequired(
            f"{self.rank_audit.finding_code}: {self.rank_audit.finding} "
            "The conservative endpoint payload is not an accepted state."
        )


def layout_for_extinction(before: cut.CutTransportState) -> ExtinctionLayout:
    """Return the exact endpoint rank for a strict cut in material cell zero."""

    if not isinstance(before, cut.CutTransportState):
        raise TypeError("extinction requires a CutTransportState")
    geometry = before.geometry
    if (
        geometry.front.regime != "partial"
        or geometry.cut_cell_index != 0
        or before.layout.wet_cell_indices != (0,)
        or before.layout.dry_cell_indices != tuple(range(geometry.master_grid.n))
    ):
        raise ExtinctionTopologyError(
            "extinction requires a strict positive-volume wet/dry cut in cell zero"
        )
    n = geometry.master_grid.n
    return ExtinctionLayout(tuple(range(n)), tuple(range(1, n + 1)), n)


def audit_extinction_rank(before: cut.CutTransportState) -> ExtinctionRankAudit:
    """Certify the exact one-scalar rank deficit before any nonlinear solve."""

    layout = layout_for_extinction(before)
    return ExtinctionRankAudit(layout, before.geometry.master_grid.areas[0])


def assemble_extinction_endpoint(
    before: cut.CutTransportState,
    candidate: ExtinctionUnknowns,
    surface_boundary: ct.PoreBoundary,
    *,
    enforce_reduced_film_thresholds: bool = True,
) -> ExtinctionAssembly:
    """Assemble the exact fully dry endpoint or reject with exact rollback."""

    try:
        return _assemble_extinction_endpoint(
            before,
            candidate,
            surface_boundary,
            enforce_reduced_film_thresholds,
        )
    except ExtinctionStepError:
        raise
    except Exception as exc:
        if not isinstance(before, cut.CutTransportState):
            raise
        raise ExtinctionStepError(
            f"extinction endpoint residual rejected with exact rollback: {exc}",
            before,
        ) from exc


def _assemble_extinction_endpoint(
    before: cut.CutTransportState,
    candidate: ExtinctionUnknowns,
    surface_boundary: ct.PoreBoundary,
    enforce_reduced_film_thresholds: bool,
) -> ExtinctionAssembly:
    if not isinstance(before, cut.CutTransportState):
        raise TypeError("before must be a CutTransportState")
    if not isinstance(candidate, ExtinctionUnknowns):
        raise TypeError("candidate must be ExtinctionUnknowns")
    if not isinstance(
        surface_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("surface boundary must be a supported pore boundary")
    if not isinstance(enforce_reduced_film_thresholds, bool):
        raise TypeError("reduced-film threshold enforcement flag must be boolean")
    layout = layout_for_extinction(before)
    if len(candidate.vector()) != layout.unknown_count:
        raise ValueError("extinction candidate does not match its endpoint rank")
    if not all(math.isfinite(value) for value in candidate.vector()):
        raise ValueError("extinction candidate values must remain finite")
    if candidate.event_time_s <= 0.0:
        raise ExtinctionTopologyError("extinction event time must be strictly positive")
    _require_length(
        "extinction dry temperatures",
        candidate.dry_temperatures_k,
        layout.dry_cell_count,
    )
    _require_length(
        "extinction dry compositions",
        candidate.dry_y_hexane,
        layout.dry_cell_count,
    )
    _require_length(
        "extinction positive-area Stefan fluxes",
        candidate.positive_area_total_stefan_fluxes_mol_m2_s,
        layout.positive_area_face_count,
    )
    cut._validate_dry_boundary(surface_boundary, before.config)  # noqa: SLF001

    endpoint = cg.partition_master_grid(before.geometry.master_grid, z=0.0)
    if (
        endpoint.front.regime != "fully_dry"
        or endpoint.front.radius_m != 0.0
        or endpoint.front.master_face_index != 0
        or endpoint.cut_cell_index is not None
        or endpoint.wet_volume_m3 != 0.0
    ):
        raise RuntimeError("geometry authority did not produce exact extinction")
    swept = cg.swept_cut_geometry(
        before.geometry,
        endpoint,
        candidate.event_time_s,
    )
    old_wet = cut._evaluate_wet_piece_states(  # noqa: SLF001
        before.layout,
        before.wet_temperatures_k,
        before.wet_retained_water_loadings,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
        before.config,
    )
    old_dry = cut._evaluate_dry_piece_states(  # noqa: SLF001
        before.layout,
        before.dry_temperatures_k,
        before.dry_y_hexane,
        before.oil_fraction_labels,
        before.config,
    )
    dry_cells = _evaluate_endpoint_dry_cells(before, candidate)
    face_fluxes = _fully_dry_face_fluxes(
        before,
        candidate,
        dry_cells,
        surface_boundary,
        enforce_reduced_film_thresholds=enforce_reduced_film_thresholds,
    )
    old_inventory, current_inventory, changes = _material_inventories(
        before,
        old_wet,
        old_dry,
        dry_cells,
    )
    residuals = _residual_blocks(
        before.geometry.master_grid,
        changes,
        face_fluxes,
        candidate.event_time_s,
    )
    rank_audit = audit_extinction_rank(before)
    ledger = _ledger(
        before,
        residuals,
        changes,
        face_fluxes,
        swept,
        candidate.event_time_s,
        rank_audit,
    )
    return ExtinctionAssembly(
        before=before,
        candidate=candidate,
        layout=layout,
        rank_audit=rank_audit,
        endpoint_geometry=endpoint,
        swept_geometry=swept,
        old_wet_cells=old_wet,
        old_dry_cells=old_dry,
        dry_cells=dry_cells,
        dry_face_fluxes=face_fluxes,
        old_inventory=old_inventory,
        current_inventory=current_inventory,
        inventory_changes=changes,
        residuals=residuals,
        ledger=ledger,
        surface_boundary=surface_boundary,
        historical_hexane_loadings=before.historical_hexane_loadings,
        oil_fraction_labels=before.oil_fraction_labels,
    )


def _evaluate_endpoint_dry_cells(
    before: cut.CutTransportState,
    candidate: ExtinctionUnknowns,
) -> tuple[cp.EquilibriumPoreState, ...]:
    cells = []
    for index, (temperature, y_hexane) in enumerate(
        zip(candidate.dry_temperatures_k, candidate.dry_y_hexane)
    ):
        pore = replace(
            before.config.dry.pore,
            w_o=before.oil_fraction_labels[index],
        )
        cut._validate_open_dry_band(  # noqa: SLF001
            temperature,
            y_hexane,
            before.config,
            pore=pore,
        )
        cells.append(
            cp.evaluate_equilibrium(
                temperature,
                before.config.dry.pressure_pa,
                y_hexane,
                pore,
            )
        )
    return tuple(cells)


def _zero_center_flux() -> cut.DryFaceFlux:
    component = cp.ComponentMolarFluxes(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    energy = cp.ComponentEnergyFlux(0.0, 0.0, 0.0, 0.0, 0.0)
    return cut.DryFaceFlux(
        component,
        energy,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        None,
        cut.DryThermodynamicForceReconstruction.CENTER_SYMMETRY,
    )


def _fully_dry_face_fluxes(
    before: cut.CutTransportState,
    candidate: ExtinctionUnknowns,
    cells: Sequence[cp.EquilibriumPoreState],
    boundary: ct.PoreBoundary,
    *,
    enforce_reduced_film_thresholds: bool = True,
) -> tuple[cut.DryFaceFlux, ...]:
    grid = before.geometry.master_grid
    fluxes = [_zero_center_flux()]
    stefan = candidate.positive_area_total_stefan_fluxes_mol_m2_s
    for face_index in range(1, grid.n):
        left_index = face_index - 1
        right_index = face_index
        left = cells[left_index]
        right = cells[right_index]
        fluxes.append(
            cut._dry_flux(  # noqa: SLF001
                left.temperature_k,
                left.y_hexane,
                right.temperature_k,
                right.y_hexane,
                grid.centers[right_index] - grid.centers[left_index],
                stefan[face_index - 1],
                replace(
                    before.config.dry.pore,
                    w_o=before.oil_fraction_labels[left_index],
                ),
                replace(
                    before.config.dry.pore,
                    w_o=before.oil_fraction_labels[right_index],
                ),
                before.config,
            )
        )
    last = cells[-1]
    outer_distance = grid.R - grid.centers[-1]
    surface_pore = replace(
        before.config.dry.pore,
        w_o=before.oil_fraction_labels[-1],
    )
    if isinstance(boundary, ct.DirichletPoreBoundary):
        # Preserve the original numerical-oracle endpoint path exactly.
        fluxes.append(
            cut._dry_flux(  # noqa: SLF001
                last.temperature_k,
                last.y_hexane,
                boundary.temperature_k,
                boundary.y_hexane,
                outer_distance,
                stefan[-1],
                surface_pore,
                surface_pore,
                before.config,
            )
        )
    else:
        reference_surface = ct.prepare_reduced_film_surface_state(
            boundary,
            cell_temperature_k=last.temperature_k,
            pressure_pa=before.config.dry.pressure_pa,
            particle_radius_m=grid.R,
            outer_half_cell_distance_m=outer_distance,
            particle_thermal_conductivity_w_m_k=(
                before.config.dry.thermal_conductivity.value_w_m_k
            ),
            binary_gas_interaction_k_wh=before.config.dry.pore.k_wh,
        )

        def evaluate_surface_flux(
            surface_temperature_k: float,
            conductive_outward_w_m2: float,
        ) -> ct.ReducedFilmSurfaceFluxEvaluation:
            cut._validate_open_dry_band(  # noqa: SLF001
                surface_temperature_k,
                boundary.y_hexane,
                before.config,
                pore=surface_pore,
            )
            trial_flux = cut._dry_flux(  # noqa: SLF001
                last.temperature_k,
                last.y_hexane,
                surface_temperature_k,
                boundary.y_hexane,
                outer_distance,
                stefan[-1],
                surface_pore,
                surface_pore,
                before.config,
                conductive_heat_flux_override_w_m2=conductive_outward_w_m2,
            )
            return ct.ReducedFilmSurfaceFluxEvaluation(
                component=trial_flux.component,
                payload=trial_flux,
                non_fourier_reduced_heat_flux_w_m2=(
                    trial_flux.reciprocal_heat_of_transport_flux_w_m2
                ),
            )

        surface = ct.solve_reduced_film_surface_closure(
            reference_surface,
            evaluate_surface_flux,
            admissible_surface_temperature_bounds_k=(
                cut._non_fourier_surface_temperature_bounds(  # noqa: SLF001
                    before.config,
                    reference_surface.surface_y_hexane,
                    surface_pore,
                )
                if cut._uses_nonzero_thermal_factors(before.config)  # noqa: SLF001
                else None
            ),
        )
        surface_flux = surface.surface_flux_evaluation.payload
        if not isinstance(surface_flux, cut.DryFaceFlux):
            raise TypeError("extinction surface closure lost its dry-flux payload")
        surface_thermo = cp.evaluate_equilibrium(
            surface.surface_temperature_k,
            before.config.dry.pressure_pa,
            surface.surface_y_hexane,
            surface_pore,
        )
        surface_energy, audit = ct.reduced_film_surface_energy_and_audit(
            surface,
            water_gas_partial_enthalpy_j_mol=(surface_thermo.water_gas_partial_enthalpy_j_mol),
            hexane_gas_partial_enthalpy_j_mol=(surface_thermo.hexane_gas_partial_enthalpy_j_mol),
            enforce_fast_mass_threshold=enforce_reduced_film_thresholds,
        )
        if surface.fourier_heat_flux_outward_w_m2 != (surface_flux.fourier_heat_flux_w_m2):
            raise RuntimeError("extinction surface closure lost its Fourier face identity")
        if surface.reciprocal_heat_of_transport_flux_outward_w_m2 != (
            surface_flux.reciprocal_heat_of_transport_flux_w_m2
        ):
            raise RuntimeError("extinction surface closure lost its reciprocal-heat identity")
        if surface_energy.conductive_heat_flux_w_m2 != (surface_flux.reduced_heat_flux_w_m2):
            raise RuntimeError(
                "extinction surface energy did not use the locally closed reduced heat"
            )
        surface_flux = replace(surface_flux, energy=surface_energy)
        fluxes.append(replace(surface_flux, surface_film_audit=audit))
    return tuple(fluxes)


def _material_inventories(
    before: cut.CutTransportState,
    old_wet: Sequence[ww.WetWaterCellState],
    old_dry: Sequence[cp.EquilibriumPoreState],
    new_dry: Sequence[cp.EquilibriumPoreState],
) -> tuple[
    ExtinctionMaterialInventory,
    ExtinctionMaterialInventory,
    ExtinctionInventoryChanges,
]:
    grid = before.geometry.master_grid
    old_partition = before.geometry.cells[0]
    old_wet_cell = old_wet[0]
    old_dry_cell = old_dry[0]
    wet_hexane_concentration = (
        before.config.wet.wet.rho_dm_p * before.historical_hexane_loadings[0] / hx.M
    )

    old_water = [
        math.fsum(
            (
                old_partition.wet_volume_m3 * old_wet_cell.retained_water_concentration_mol_m3,
                old_partition.dry_volume_m3 * old_dry_cell.total_water_concentration_mol_m3,
            )
        )
    ]
    old_hexane = [
        math.fsum(
            (
                old_partition.wet_volume_m3 * wet_hexane_concentration,
                old_partition.dry_volume_m3 * old_dry_cell.total_hexane_concentration_mol_m3,
            )
        )
    ]
    old_energy = [
        math.fsum(
            (
                old_partition.wet_volume_m3 * old_wet_cell.energy_density_j_m3,
                old_partition.dry_volume_m3 * old_dry_cell.energy_density_j_m3,
            )
        )
    ]
    for index, cell in enumerate(old_dry[1:], start=1):
        volume = grid.volumes[index]
        old_water.append(volume * cell.total_water_concentration_mol_m3)
        old_hexane.append(volume * cell.total_hexane_concentration_mol_m3)
        old_energy.append(volume * cell.energy_density_j_m3)

    new_water = tuple(
        volume * cell.total_water_concentration_mol_m3
        for volume, cell in zip(grid.volumes, new_dry)
    )
    new_hexane = tuple(
        volume * cell.total_hexane_concentration_mol_m3
        for volume, cell in zip(grid.volumes, new_dry)
    )
    new_energy = tuple(
        volume * cell.energy_density_j_m3 for volume, cell in zip(grid.volumes, new_dry)
    )

    def combined_change(
        new_concentration: float,
        old_wet_concentration: float,
        old_dry_concentration: float,
    ) -> float:
        return math.fsum(
            (
                old_partition.wet_volume_m3 * (new_concentration - old_wet_concentration),
                old_partition.dry_volume_m3 * (new_concentration - old_dry_concentration),
            )
        )

    water_changes = [
        combined_change(
            new_dry[0].total_water_concentration_mol_m3,
            old_wet_cell.retained_water_concentration_mol_m3,
            old_dry_cell.total_water_concentration_mol_m3,
        )
    ]
    hexane_changes = [
        combined_change(
            new_dry[0].total_hexane_concentration_mol_m3,
            wet_hexane_concentration,
            old_dry_cell.total_hexane_concentration_mol_m3,
        )
    ]
    energy_changes = [
        combined_change(
            new_dry[0].energy_density_j_m3,
            old_wet_cell.energy_density_j_m3,
            old_dry_cell.energy_density_j_m3,
        )
    ]
    for index, (old_cell, new_cell) in enumerate(
        zip(old_dry[1:], new_dry[1:]),
        start=1,
    ):
        volume = grid.volumes[index]
        water_changes.append(
            volume
            * (
                new_cell.total_water_concentration_mol_m3
                - old_cell.total_water_concentration_mol_m3
            )
        )
        hexane_changes.append(
            volume
            * (
                new_cell.total_hexane_concentration_mol_m3
                - old_cell.total_hexane_concentration_mol_m3
            )
        )
        energy_changes.append(
            volume * (new_cell.energy_density_j_m3 - old_cell.energy_density_j_m3)
        )
    return (
        ExtinctionMaterialInventory(tuple(old_water), tuple(old_hexane), tuple(old_energy)),
        ExtinctionMaterialInventory(new_water, new_hexane, new_energy),
        ExtinctionInventoryChanges(
            tuple(water_changes),
            tuple(hexane_changes),
            tuple(energy_changes),
        ),
    )


def _residual_blocks(
    grid,
    changes: ExtinctionInventoryChanges,
    fluxes: Sequence[cut.DryFaceFlux],
    event_time_s: float,
) -> ExtinctionResidualBlocks:
    water_rates = tuple(
        area * flux.component.conserved_water_flux_mol_m2_s
        for area, flux in zip(grid.areas, fluxes)
    )
    hexane_rates = tuple(
        area * flux.component.conserved_hexane_flux_mol_m2_s
        for area, flux in zip(grid.areas, fluxes)
    )
    energy_rates = tuple(
        area * flux.energy.total_energy_flux_w_m2 for area, flux in zip(grid.areas, fluxes)
    )

    def regional(
        inventory_changes: Sequence[float],
        rates: Sequence[float],
    ) -> tuple[float, ...]:
        return tuple(
            math.fsum(
                (
                    change / event_time_s,
                    rates[index + 1],
                    -rates[index],
                )
            )
            for index, change in enumerate(inventory_changes)
        )

    water = regional(changes.water_cell_mol, water_rates)
    hexane = regional(changes.hexane_cell_mol, hexane_rates)
    energy = regional(changes.energy_cell_j, energy_rates)
    covariant_energy = tuple(
        math.fsum((energy_value, -hx.H_REF * hexane_value))
        for energy_value, hexane_value in zip(energy, hexane)
    )
    return ExtinctionResidualBlocks(water, hexane, energy, covariant_energy)


def _ledger(
    before: cut.CutTransportState,
    residuals: ExtinctionResidualBlocks,
    changes: ExtinctionInventoryChanges,
    fluxes: Sequence[cut.DryFaceFlux],
    swept: cg.SweptCutGeometry,
    event_time_s: float,
    rank_audit: ExtinctionRankAudit,
) -> ExtinctionLedger:
    grid = before.geometry.master_grid
    surface_area = grid.areas[-1]
    surface = fluxes[-1]
    water_surface = surface_area * surface.component.conserved_water_flux_mol_m2_s
    hexane_surface = surface_area * surface.component.conserved_hexane_flux_mol_m2_s
    energy_surface = surface_area * surface.energy.total_energy_flux_w_m2
    covariant_energy_surface = math.fsum((energy_surface, -hx.H_REF * hexane_surface))
    water_regions = math.fsum(residuals.dry_water_mol_s)
    hexane_regions = math.fsum(residuals.dry_hexane_mol_s)
    energy_regions = math.fsum(residuals.dry_energy_w)
    covariant_energy_regions = math.fsum(residuals.dry_datum_covariant_energy_w)
    water_direct = math.fsum((math.fsum(changes.water_cell_mol) / event_time_s, water_surface))
    hexane_direct = math.fsum((math.fsum(changes.hexane_cell_mol) / event_time_s, hexane_surface))
    energy_direct = math.fsum((math.fsum(changes.energy_cell_j) / event_time_s, energy_surface))
    covariant_energy_changes = tuple(
        math.fsum((energy, -hx.H_REF * hexane))
        for energy, hexane in zip(
            changes.energy_cell_j,
            changes.hexane_cell_mol,
        )
    )
    covariant_energy_direct = math.fsum(
        (
            math.fsum(covariant_energy_changes) / event_time_s,
            covariant_energy_surface,
        )
    )
    water_error = water_regions - water_direct
    hexane_error = hexane_regions - hexane_direct
    energy_error = energy_regions - energy_direct
    covariant_energy_error = covariant_energy_regions - covariant_energy_direct

    def relative(error: float, *values: float) -> float:
        return abs(error) / max(*(abs(value) for value in values), 1.0e-300)

    return ExtinctionLedger(
        water_global_from_regions_mol_s=water_regions,
        water_global_direct_mol_s=water_direct,
        water_telescoping_error_mol_s=water_error,
        hexane_global_from_regions_mol_s=hexane_regions,
        hexane_global_direct_mol_s=hexane_direct,
        hexane_telescoping_error_mol_s=hexane_error,
        energy_global_from_regions_w=energy_regions,
        energy_global_direct_w=energy_direct,
        energy_telescoping_error_w=energy_error,
        datum_covariant_energy_global_from_regions_w=(covariant_energy_regions),
        datum_covariant_energy_global_direct_w=covariant_energy_direct,
        datum_covariant_energy_telescoping_error_w=covariant_energy_error,
        water_telescoping_relative_error=relative(
            water_error,
            water_regions,
            water_direct,
        ),
        hexane_telescoping_relative_error=relative(
            hexane_error,
            hexane_regions,
            hexane_direct,
        ),
        energy_telescoping_relative_error=relative(
            energy_error,
            energy_regions,
            energy_direct,
        ),
        datum_covariant_energy_telescoping_relative_error=relative(
            covariant_energy_error,
            covariant_energy_regions,
            covariant_energy_direct,
        ),
        geometry_cell_gcl_residual_m3=swept.max_cell_gcl_residual_m3,
        center_face_area_m2=grid.areas[0],
        center_total_stefan_flux_mol_m2_s=(fluxes[0].component.total_stefan_flux_mol_m2_s),
        minimum_entropy_production_w_m3_k=min(
            flux.total_entropy_production_w_m3_k for flux in fluxes
        ),
        event_time_rank_deficit=rank_audit.unknown_time_rank_deficit,
        endpoint_interface_equation_count=(rank_audit.endpoint_interface_equation_count),
        historical_label_change=0.0,
        oil_label_change=0.0,
    )


def _require_length(name: str, values: Sequence[float], expected: int) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} has length {len(values)}, expected {expected}")


__all__ = [
    "ExtinctionAssembly",
    "ExtinctionAuthorityRequired",
    "ExtinctionInventoryChanges",
    "ExtinctionLayout",
    "ExtinctionLedger",
    "ExtinctionMaterialInventory",
    "ExtinctionRankAudit",
    "ExtinctionResidualBlocks",
    "ExtinctionStepError",
    "ExtinctionTopologyError",
    "ExtinctionUnknowns",
    "assemble_extinction_endpoint",
    "audit_extinction_rank",
    "layout_for_extinction",
]
