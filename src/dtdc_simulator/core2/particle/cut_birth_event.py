r"""Exact fully-wet to first strict-partial dry-shell birth foundation.

This module implements the owner-approved ``GT-PS-2-G1G-01-A1`` revision-2
zero-old-dry-volume limit.  The before state is an exact fully-wet stationary
material grid.  The candidate is a strict partial front in the outermost
material cell.  Its first nonzero dry water, n-hexane, and common-datum energy
inventories are accumulated in the same backward-Euler/ALE/Rankine--Hugoniot
residual as the wet-region update and surface fluxes.  There is no isolated
shell inversion, dry-state impulse, or epsilon dry seed.

The receding chart has ``2*N + 6`` unknowns and equations: two primitives for
each of ``N`` wet pieces; one dry temperature and composition; two algebraic
dry-face total Stefan fluxes; front ``z``; interface temperature; corresponding
wet water/energy, dry water/hexane/energy, and three RH rows.  Water remains on
the stationary material cell and is repartitioned locally only through these
coupled equations.  Every still-wet piece retains its immutable historical
n-hexane label.

The complete finite-thickness external-film law needed to *solve* the fully-wet
contact branch is not frozen.  Consequently this foundation does not fabricate
one.  :func:`classify_endpoint_branch` admits contact only from independently
converged evidence with exactly zero outward n-hexane transfer and closed wet
water/common-energy equations.  A positive outward flux requires immediate
recession; an inward tendency and every numerical failure are explicit
non-contact outcomes.

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

All boundary traces and coefficients used here remain numerical and
``physically_qualifying=False``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Sequence

from dtdc_simulator.core2.particle import coupled_pore as cp
from dtdc_simulator.core2.particle import coupled_transport as ct
from dtdc_simulator.core2.particle import cut_geometry as cg
from dtdc_simulator.core2.particle import cut_transport as cut
from dtdc_simulator.core2.particle import wet_water as ww
from dtdc_simulator.core2.props import hexane as hx


class BirthTopologyError(ValueError):
    """Candidate or endpoint evidence lies outside the frozen birth chart."""


class BirthStepError(RuntimeError):
    """Rejected birth residual carrying the exact immutable fully-wet state."""

    def __init__(self, message: str, rollback_state: ww.WetWaterState) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state


class BirthBranch(str, Enum):
    """Exact-sign endpoint active-set outcomes."""

    RECEDING = "receding_birth"
    CONTACT = "fully_wet_zero_flux_contact"
    INWARD_OUTSIDE_ENVELOPE = "inward_outside_primary_drainage_envelope"
    REJECTED_NUMERICAL_FAILURE = "rejected_numerical_failure"
    REJECTED_MISSING_BRANCH_SOLUTION = "rejected_missing_branch_solution"


@dataclass(frozen=True)
class BirthLayout:
    """Square first-birth chart: all wet cells plus one newborn dry piece."""

    wet_cell_indices: tuple[int, ...]
    dry_cell_indices: tuple[int, ...]
    cut_cell_index: int
    master_cell_count: int
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if self.master_cell_count < 1:
            raise ValueError("birth layout needs at least one material cell")
        expected_wet = tuple(range(self.master_cell_count))
        expected_dry = (self.master_cell_count - 1,)
        if self.wet_cell_indices != expected_wet:
            raise ValueError("first birth must retain every material cell on the wet side")
        if self.dry_cell_indices != expected_dry:
            raise ValueError("first birth must create only the outer-cell dry piece")
        if self.cut_cell_index != self.master_cell_count - 1:
            raise ValueError("first birth cut cell must be the outermost material cell")

    @property
    def wet_piece_count(self) -> int:
        return len(self.wet_cell_indices)

    @property
    def dry_piece_count(self) -> int:
        return 1

    @property
    def dry_face_count(self) -> int:
        return 2

    @property
    def unknown_blocks(self) -> tuple[tuple[str, int], ...]:
        return (
            ("wet_temperature", self.wet_piece_count),
            ("wet_retained_water", self.wet_piece_count),
            ("dry_temperature", 1),
            ("dry_y_hexane", 1),
            ("dry_total_stefan_flux", 2),
            ("front_z", 1),
            ("interface_temperature", 1),
        )

    @property
    def residual_blocks(self) -> tuple[tuple[str, int], ...]:
        return (
            ("wet_water", self.wet_piece_count),
            ("wet_energy", self.wet_piece_count),
            ("dry_water", 1),
            ("dry_hexane", 1),
            ("dry_energy", 1),
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
        return 2 * self.master_cell_count + 6

    @property
    def is_square(self) -> bool:
        return self.unknown_count == self.residual_count == self.rank_formula


@dataclass(frozen=True)
class BirthUnknowns:
    """First-partial state in the documented square-rank ordering."""

    wet_temperatures_k: tuple[float, ...]
    wet_retained_water_loadings: tuple[float, ...]
    dry_temperature_k: float
    dry_y_hexane: float
    dry_total_stefan_fluxes_mol_m2_s: tuple[float, float]
    front_z: float
    interface_temperature_k: float
    physically_qualifying: bool = field(default=False, init=False)

    def vector(self) -> tuple[float, ...]:
        return (
            *self.wet_temperatures_k,
            *self.wet_retained_water_loadings,
            self.dry_temperature_k,
            self.dry_y_hexane,
            *self.dry_total_stefan_fluxes_mol_m2_s,
            self.front_z,
            self.interface_temperature_k,
        )

    @classmethod
    def from_vector(
        cls,
        layout: BirthLayout,
        values: Sequence[float],
    ) -> "BirthUnknowns":
        vector = tuple(values)
        if len(vector) != layout.unknown_count:
            raise ValueError(
                f"birth vector has length {len(vector)}, expected {layout.unknown_count}"
            )
        if not all(math.isfinite(value) for value in vector):
            raise ValueError("birth candidate vector must remain finite")
        n = layout.wet_piece_count
        return cls(
            wet_temperatures_k=vector[:n],
            wet_retained_water_loadings=vector[n : 2 * n],
            dry_temperature_k=vector[2 * n],
            dry_y_hexane=vector[2 * n + 1],
            dry_total_stefan_fluxes_mol_m2_s=(
                vector[2 * n + 2],
                vector[2 * n + 3],
            ),
            front_z=vector[2 * n + 4],
            interface_temperature_k=vector[2 * n + 5],
        )

    def as_cut_unknowns(self) -> cut.CutTransportUnknowns:
        """Adapt to the shared strict-partial constitutive kernels."""

        return cut.CutTransportUnknowns(
            wet_temperatures_k=self.wet_temperatures_k,
            wet_retained_water_loadings=self.wet_retained_water_loadings,
            dry_temperatures_k=(self.dry_temperature_k,),
            dry_y_hexane=(self.dry_y_hexane,),
            dry_total_stefan_fluxes_mol_m2_s=(
                self.dry_total_stefan_fluxes_mol_m2_s
            ),
            front_z=self.front_z,
            interface_temperature_k=self.interface_temperature_k,
        )


@dataclass(frozen=True)
class BirthInventoryChanges:
    """Stable material-piece changes from an exactly zero old dry inventory."""

    wet_water_mol: tuple[float, ...]
    wet_hexane_mol: tuple[float, ...]
    wet_energy_j: tuple[float, ...]
    dry_water_mol: tuple[float]
    dry_hexane_mol: tuple[float]
    dry_energy_j: tuple[float]
    old_dry_water_mol: float = field(default=0.0, init=False)
    old_dry_hexane_mol: float = field(default=0.0, init=False)
    old_dry_energy_j: float = field(default=0.0, init=False)
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class BirthLedger:
    """Independent common-frame global and zero-old-volume checks."""

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
    old_dry_volume_m3: float
    newborn_dry_volume_m3: float
    swept_dry_volume_m3: float
    newborn_water_inventory_identity_mol: float
    newborn_hexane_inventory_identity_mol: float
    newborn_energy_inventory_identity_j: float
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
class BirthAssembly:
    """Complete first-birth residual evaluation without state commitment."""

    before: ww.WetWaterState
    config: cut.CutTransportConfig
    candidate: BirthUnknowns
    layout: BirthLayout
    before_geometry: cg.CutGeometry
    candidate_geometry: cg.CutGeometry
    swept_geometry: cg.SweptCutGeometry
    old_wet_cells: tuple[ww.WetWaterCellState, ...]
    wet_cells: tuple[ww.WetWaterCellState, ...]
    dry_cells: tuple[cp.EquilibriumPoreState]
    interface: cut.CutInterfaceState
    wet_face_fluxes: tuple[ww.WetWaterFaceFlux, ...]
    dry_face_fluxes: tuple[cut.DryFaceFlux, cut.DryFaceFlux]
    inventory_changes: BirthInventoryChanges
    residuals: cut.CutResidualBlocks
    ledger: BirthLedger
    surface_boundary: ct.PoreBoundary
    dt_s: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not self.layout.is_square:
            raise RuntimeError("birth residual lost its square rank")
        if len(self.candidate.vector()) != self.layout.unknown_count:
            raise RuntimeError("birth candidate does not match its declared rank")
        if len(self.residuals.vector) != self.layout.residual_count:
            raise RuntimeError("birth residual does not match its declared rank")
        if self.before_geometry.dry_volume_m3 != 0.0:
            raise RuntimeError("birth before geometry must have exactly zero dry volume")

    @property
    def average_recession_speed_m_s(self) -> float:
        return (
            self.before.grid.R - self.candidate_geometry.front.radius_m
        ) / self.dt_s

    @property
    def surface_hexane_flux_mol_m2_s(self) -> float:
        return self.dry_face_fluxes[-1].component.conserved_hexane_flux_mol_m2_s


# GT-PS-2-P1E-R15 (2026-08-05): birth-subsystem regional/RH acceptance.
# The newborn cell's near-zero volume amplifies the birth residual scale;
# measured convergence points at the pinned joint forcing straddle the
# shared 1e-10 gate (N1: 7.744e-11, N2: 1.104e-10).  This bound is >2x the
# worst measured point, applies ONLY to the birth solve's regional/RH
# acceptance and branch evidence (conservation ledgers and the shared
# trajectory tolerance are unchanged); the A21-grade scale-authority
# derivation is queued with A22.
BIRTH_REGIONAL_RH_ACCEPTANCE = 2.5e-10


@dataclass(frozen=True)
class RecedingBranchEvidence:
    """Evidence from a separately safeguarded nonlinear birth solve."""

    independently_converged: bool
    front_z: float
    maximum_scaled_residual: float
    nonlinear_tolerance: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.front_z,
            self.maximum_scaled_residual,
            self.nonlinear_tolerance,
        )
        if not all(math.isfinite(value) for value in values):
            raise ValueError("receding evidence must remain finite")
        if not 0.0 < self.front_z < 1.0:
            raise ValueError("receding evidence requires a strict partial front")
        if self.maximum_scaled_residual < 0.0:
            raise ValueError("scaled residual cannot be negative")
        if not 0.0 < self.nonlinear_tolerance <= BIRTH_REGIONAL_RH_ACCEPTANCE:
            raise ValueError(
                "birth nonlinear tolerance must lie in "
                f"(0,{BIRTH_REGIONAL_RH_ACCEPTANCE:g}]"
            )

    @property
    def certified(self) -> bool:
        return self.independently_converged and (
            self.maximum_scaled_residual <= self.nonlinear_tolerance
        )


@dataclass(frozen=True)
class FullyWetContactEvidence:
    """Complete independently solved contact evidence; never a fallback flag."""

    independently_converged: bool
    outward_hexane_flux_mol_m2_s: float
    outward_water_flux_mol_m2_s: float
    outward_energy_flux_w_m2: float
    hexane_thermodynamic_driving_force: float
    interfacial_hexane_access_multiplier: float
    scaled_water_residuals: tuple[float, ...]
    scaled_energy_residuals: tuple[float, ...]
    nonlinear_tolerance: float
    external_active_set_admissible: bool
    property_bounds_admissible: bool
    historical_hexane_labels_preserved: bool
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        residuals = (*self.scaled_water_residuals, *self.scaled_energy_residuals)
        if not residuals:
            raise ValueError("contact evidence needs water and energy balance rows")
        values = (
            self.outward_hexane_flux_mol_m2_s,
            self.outward_water_flux_mol_m2_s,
            self.outward_energy_flux_w_m2,
            self.hexane_thermodynamic_driving_force,
            self.interfacial_hexane_access_multiplier,
            self.nonlinear_tolerance,
        )
        if not all(math.isfinite(value) for value in (*values, *residuals)):
            raise ValueError("contact evidence must remain finite")
        if not 0.0 <= self.interfacial_hexane_access_multiplier <= 1.0:
            raise ValueError("contact n-hexane access multiplier must lie in [0,1]")
        if not 0.0 < self.nonlinear_tolerance <= 1.0e-10:
            raise ValueError("contact nonlinear tolerance must lie in (0,1e-10]")

    @property
    def maximum_scaled_residual(self) -> float:
        return max(
            abs(value)
            for value in (*self.scaled_water_residuals, *self.scaled_energy_residuals)
        )

    @property
    def certified(self) -> bool:
        return (
            self.independently_converged
            and self.outward_hexane_flux_mol_m2_s == 0.0
            and self.maximum_scaled_residual <= self.nonlinear_tolerance
            and self.external_active_set_admissible
            and self.property_bounds_admissible
            and self.historical_hexane_labels_preserved
        )


@dataclass(frozen=True)
class BirthBranchDecision:
    """Auditable exact-sign endpoint decision."""

    branch: BirthBranch
    outward_hexane_flux_mol_m2_s: float
    hexane_storage_jump_mol_m3: float
    recession_speed_m_s: float | None
    reason: str
    contact_evidence_used: bool
    physically_qualifying: bool = field(default=False, init=False)


def layout_for_birth(before: ww.WetWaterState) -> BirthLayout:
    """Return the unique first-partial chart for an exact fully-wet state."""

    if not isinstance(before, ww.WetWaterState):
        raise TypeError("birth requires a fully-wet WetWaterState")
    n = before.grid.n
    layout = BirthLayout(tuple(range(n)), (n - 1,), n - 1, n)
    if not layout.is_square:
        raise RuntimeError("birth rank formula is internally inconsistent")
    return layout


def assemble_birth_backward_euler(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
    candidate: BirthUnknowns,
    dt_s: float,
    surface_boundary: ct.PoreBoundary,
    *,
    enforce_reduced_film_thresholds: bool = True,
) -> BirthAssembly:
    """Assemble one exact zero-old-dry-volume birth residual or roll back."""

    try:
        return _assemble_birth_backward_euler(
            before,
            config,
            candidate,
            dt_s,
            surface_boundary,
            enforce_reduced_film_thresholds,
        )
    except BirthStepError:
        raise
    except Exception as exc:
        if not isinstance(before, ww.WetWaterState):
            raise
        raise BirthStepError(
            f"dry-shell birth residual rejected with exact rollback: {exc}",
            before,
        ) from exc


def _assemble_birth_backward_euler(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
    candidate: BirthUnknowns,
    dt_s: float,
    surface_boundary: ct.PoreBoundary,
    enforce_reduced_film_thresholds: bool,
) -> BirthAssembly:
    if not isinstance(before, ww.WetWaterState):
        raise TypeError("before must be a WetWaterState")
    if not isinstance(config, cut.CutTransportConfig):
        raise TypeError("birth requires a CutTransportConfig")
    if before.model != config.wet:
        raise ValueError("fully-wet state and birth chart must share one wet authority")
    if not isinstance(candidate, BirthUnknowns):
        raise TypeError("candidate must be BirthUnknowns")
    if not isinstance(
        surface_boundary,
        (ct.DirichletPoreBoundary, ct.ReducedFilmPoreBoundary),
    ):
        raise TypeError("surface boundary must be a supported pore boundary")
    if not isinstance(enforce_reduced_film_thresholds, bool):
        raise TypeError("reduced-film threshold enforcement flag must be boolean")
    if not math.isfinite(dt_s) or dt_s <= 0.0:
        raise ValueError("birth step duration must be positive and finite")
    if not math.isfinite(before.time_s) or before.time_s < 0.0:
        raise ValueError("fully-wet state time must be finite and non-negative")
    layout = layout_for_birth(before)
    for name, values in (
        ("fully-wet temperatures", before.temperatures_k),
        ("fully-wet retained water", before.retained_water_loadings),
        ("historical n-hexane labels", before.historical_hexane_loadings),
        ("oil-fraction labels", before.oil_fraction_labels),
    ):
        _require_length(name, values, layout.master_cell_count)
    if len(candidate.vector()) != layout.unknown_count:
        raise ValueError("birth candidate does not match the first-partial rank")
    if not all(math.isfinite(value) for value in candidate.vector()):
        raise ValueError("birth candidate values must remain finite")
    _require_length(
        "birth wet temperatures",
        candidate.wet_temperatures_k,
        layout.wet_piece_count,
    )
    _require_length(
        "birth wet retained water",
        candidate.wet_retained_water_loadings,
        layout.wet_piece_count,
    )
    cut._validate_dry_boundary(surface_boundary, config)  # noqa: SLF001

    before_geometry = cg.partition_master_grid(before.grid, z=1.0)
    if before_geometry.front.regime != "fully_wet":
        raise RuntimeError("birth before geometry lost its exact endpoint")
    candidate_geometry = cg.partition_master_grid(before.grid, z=candidate.front_z)
    cut_layout = cut.layout_for_geometry(candidate_geometry)
    if (
        cut_layout.wet_cell_indices != layout.wet_cell_indices
        or cut_layout.dry_cell_indices != layout.dry_cell_indices
        or cut_layout.cut_cell_index != layout.cut_cell_index
    ):
        raise BirthTopologyError(
            "first birth must be a strict partial cut in the outermost material cell"
        )
    swept = cg.swept_cut_geometry(before_geometry, candidate_geometry, dt_s)
    old_wet = tuple(
        ww.evaluate_cell(T, W, Xh, wo, config.wet)
        for T, W, Xh, wo in zip(
            before.temperatures_k,
            before.retained_water_loadings,
            before.historical_hexane_loadings,
            before.oil_fraction_labels,
        )
    )
    wet_cells = cut._evaluate_wet_piece_states(  # noqa: SLF001
        cut_layout,
        candidate.wet_temperatures_k,
        candidate.wet_retained_water_loadings,
        before.historical_hexane_loadings,
        before.oil_fraction_labels,
        config,
    )
    dry_cells = cut._evaluate_dry_piece_states(  # noqa: SLF001
        cut_layout,
        (candidate.dry_temperature_k,),
        (candidate.dry_y_hexane,),
        before.oil_fraction_labels,
        config,
    )
    interface = cut.evaluate_interface_state(
        candidate.interface_temperature_k,
        config,
        before.historical_hexane_loadings[-1],
        before.oil_fraction_labels[-1],
    )
    wet_fluxes = cut._wet_face_fluxes(  # noqa: SLF001
        candidate_geometry,
        cut_layout,
        wet_cells,
        interface,
        config,
    )
    dry_fluxes = cut._dry_face_fluxes(  # noqa: SLF001
        candidate_geometry,
        cut_layout,
        dry_cells,
        interface,
        candidate.dry_total_stefan_fluxes_mol_m2_s,
        surface_boundary,
        before.oil_fraction_labels,
        config,
        enforce_reduced_film_thresholds=enforce_reduced_film_thresholds,
    )
    changes = _inventory_changes(
        before,
        config,
        candidate_geometry,
        cut_layout,
        old_wet,
        wet_cells,
        dry_cells,
        swept,
    )
    residuals, wet_hexane_identity = _residual_blocks(
        config,
        candidate_geometry,
        cut_layout,
        changes,
        interface,
        wet_fluxes,
        dry_fluxes,
        swept,
        dt_s,
        # O9a front-donor correction: the piece the newborn front is consuming
        # is the OUTERMOST wet piece (the cut cell's wet part), and the donor is
        # that piece's own bulk retained-water concentration.
        #
        # THE TIME LEVEL IS THE MODULE'S OWN, NOT A NEW CHOICE.
        # ``_inventory_changes`` books the matching geometric storage term as
        # ``Delta(V)*C_old`` in ``stable_wet_change``; the ALE sweep is that
        # identical swept volume seen from the flux side, so the donor is read
        # at the same time level - ``C_old`` of the consumed piece.  The two
        # then cancel term by term and the row keeps exactly the ``V_new``
        # storage form the birth decomposition is written in.
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
    ledger = _ledger(
        before,
        changes,
        residuals,
        wet_hexane_identity,
        dry_fluxes[-1],
        candidate_geometry.master_grid.areas[-1],
        swept,
        dt_s,
        (*wet_fluxes, *dry_fluxes),
        candidate_geometry,
        dry_cells[0],
    )
    return BirthAssembly(
        before=before,
        config=config,
        candidate=candidate,
        layout=layout,
        before_geometry=before_geometry,
        candidate_geometry=candidate_geometry,
        swept_geometry=swept,
        old_wet_cells=old_wet,
        wet_cells=wet_cells,
        dry_cells=(dry_cells[0],),
        interface=interface,
        wet_face_fluxes=wet_fluxes,
        dry_face_fluxes=(dry_fluxes[0], dry_fluxes[1]),
        inventory_changes=changes,
        residuals=residuals,
        ledger=ledger,
        surface_boundary=surface_boundary,
        dt_s=dt_s,
    )


def _inventory_changes(
    before: ww.WetWaterState,
    config: cut.CutTransportConfig,
    geometry: cg.CutGeometry,
    layout: cut.CutTransportLayout,
    old_wet: Sequence[ww.WetWaterCellState],
    new_wet: Sequence[ww.WetWaterCellState],
    new_dry: Sequence[cp.EquilibriumPoreState],
    swept: cg.SweptCutGeometry,
) -> BirthInventoryChanges:
    wet_volumes = tuple(
        geometry.cells[index].wet_volume_m3 for index in layout.wet_cell_indices
    )
    wet_volume_changes = tuple(
        swept.cell_wet_volume_changes_m3[index]
        for index in layout.wet_cell_indices
    )

    def stable_wet_change(
        old_concentrations: Sequence[float],
        new_concentrations: Sequence[float],
    ) -> tuple[float, ...]:
        return tuple(
            math.fsum(
                (
                    volume * (new - old),
                    volume_change * old,
                )
            )
            for volume, volume_change, old, new in zip(
                wet_volumes,
                wet_volume_changes,
                old_concentrations,
                new_concentrations,
            )
        )

    old_hexane = tuple(
        config.wet.wet.rho_dm_p * loading / hx.M
        for loading in before.historical_hexane_loadings
    )
    wet_water = stable_wet_change(
        tuple(cell.retained_water_concentration_mol_m3 for cell in old_wet),
        tuple(cell.retained_water_concentration_mol_m3 for cell in new_wet),
    )
    wet_hexane = stable_wet_change(old_hexane, old_hexane)
    wet_energy = stable_wet_change(
        tuple(cell.energy_density_j_m3 for cell in old_wet),
        tuple(cell.energy_density_j_m3 for cell in new_wet),
    )
    dry_volume = geometry.cells[layout.cut_cell_index].dry_volume_m3
    swept_volume = swept.cell_dry_volume_changes_m3[layout.cut_cell_index]
    if dry_volume != swept_volume:
        tolerance = 64.0 * math.ulp(geometry.master_grid.total_volume)
        if abs(dry_volume - swept_volume) > tolerance:
            raise RuntimeError("newborn dry volume differs from the exact ALE sweep")
    dry = new_dry[0]
    dry_water, dry_hexane, dry_energy = zero_old_dry_inventories(
        dry_volume,
        dry,
    )
    return BirthInventoryChanges(
        wet_water_mol=wet_water,
        wet_hexane_mol=wet_hexane,
        wet_energy_j=wet_energy,
        dry_water_mol=(dry_water,),
        dry_hexane_mol=(dry_hexane,),
        dry_energy_j=(dry_energy,),
    )


def zero_old_dry_inventories(
    new_dry_volume_m3: float,
    new_dry_state: cp.EquilibriumPoreState,
) -> tuple[float, float, float]:
    """Return first dry inventories from the exact ``V_old,d=0`` limit.

    This local operation is independent of the outer-endpoint geometry and is
    intentionally reusable when a later inner master-face departure creates a
    zero-old-dry subpiece next to already existing outer dry cells.  It never
    evaluates a fictitious dry state at zero volume.
    """

    if not math.isfinite(new_dry_volume_m3) or new_dry_volume_m3 <= 0.0:
        raise BirthTopologyError("newborn dry volume must be strictly positive")
    if not isinstance(new_dry_state, cp.EquilibriumPoreState):
        raise TypeError("newborn inventory requires an equilibrium dry state")
    return (
        new_dry_volume_m3 * new_dry_state.total_water_concentration_mol_m3,
        new_dry_volume_m3 * new_dry_state.total_hexane_concentration_mol_m3,
        new_dry_volume_m3 * new_dry_state.energy_density_j_m3,
    )


def _residual_blocks(
    config: cut.CutTransportConfig,
    geometry: cg.CutGeometry,
    layout: cut.CutTransportLayout,
    changes: BirthInventoryChanges,
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
    # this sibling chart in Stage 1b for formulation consistency: the birth
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
    surface_area = grid.areas[-1]
    dry_water_rates = (
        dry_water_h,
        surface_area * dry_fluxes[-1].component.conserved_water_flux_mol_m2_s,
    )
    dry_hexane_rates = (
        dry_hexane_h,
        surface_area * dry_fluxes[-1].component.conserved_hexane_flux_mol_m2_s,
    )
    dry_energy_rates = (
        dry_energy_h,
        surface_area * dry_fluxes[-1].energy.total_energy_flux_w_m2,
    )

    def single_dry_residual(change: float, rates: Sequence[float]) -> float:
        return math.fsum((change / dt_s, rates[1], -rates[0]))

    wet_hexane_identity = math.fsum(
        (math.fsum(changes.wet_hexane_mol) / dt_s, wet_hexane_h)
    )
    return (
        cut.CutResidualBlocks(
            wet_water_mol_s=wet_water_residuals,
            wet_energy_w=wet_energy_residuals,
            dry_water_mol_s=(
                single_dry_residual(changes.dry_water_mol[0], dry_water_rates),
            ),
            dry_hexane_mol_s=(
                single_dry_residual(changes.dry_hexane_mol[0], dry_hexane_rates),
            ),
            dry_energy_w=(
                single_dry_residual(changes.dry_energy_j[0], dry_energy_rates),
            ),
            rh_water_mol_s=rh_water,
            rh_hexane_mol_s=rh_hexane,
            rh_energy_w=rh_energy,
        ),
        wet_hexane_identity,
    )


def _ledger(
    before: ww.WetWaterState,
    changes: BirthInventoryChanges,
    residuals: cut.CutResidualBlocks,
    wet_hexane_identity: float,
    surface: cut.DryFaceFlux,
    surface_area_m2: float,
    swept: cg.SweptCutGeometry,
    dt_s: float,
    all_fluxes: Sequence[object],
    geometry: cg.CutGeometry,
    dry_cell: cp.EquilibriumPoreState,
) -> BirthLedger:
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
    water_direct = math.fsum(
        (
            math.fsum((*changes.wet_water_mol, *changes.dry_water_mol)) / dt_s,
            surface_area_m2 * surface.component.conserved_water_flux_mol_m2_s,
        )
    )
    hexane_direct = math.fsum(
        (
            math.fsum((*changes.wet_hexane_mol, *changes.dry_hexane_mol)) / dt_s,
            surface_area_m2 * surface.component.conserved_hexane_flux_mol_m2_s,
        )
    )
    energy_direct = math.fsum(
        (
            math.fsum((*changes.wet_energy_j, *changes.dry_energy_j)) / dt_s,
            surface_area_m2 * surface.energy.total_energy_flux_w_m2,
        )
    )
    water_error = water_blocks - water_direct
    hexane_error = hexane_blocks - hexane_direct
    energy_error = energy_blocks - energy_direct
    water_scale = max(abs(water_blocks), abs(water_direct), 1.0e-300)
    hexane_scale = max(abs(hexane_blocks), abs(hexane_direct), 1.0e-300)
    energy_scale = max(abs(energy_blocks), abs(energy_direct), 1.0e-300)
    entropy_values: list[float] = []
    for flux in all_fluxes:
        if isinstance(flux, ww.WetWaterFaceFlux):
            entropy_values.append(flux.total_entropy_production_w_m3_k)
        elif isinstance(flux, cut.DryFaceFlux):
            entropy_values.append(flux.total_entropy_production_w_m3_k)
    dry_volume = geometry.cells[-1].dry_volume_m3
    swept_volume = swept.cell_dry_volume_changes_m3[-1]
    return BirthLedger(
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
        old_dry_volume_m3=0.0,
        newborn_dry_volume_m3=dry_volume,
        swept_dry_volume_m3=swept_volume,
        newborn_water_inventory_identity_mol=math.fsum(
            (
                changes.dry_water_mol[0],
                -dry_volume * dry_cell.total_water_concentration_mol_m3,
            )
        ),
        newborn_hexane_inventory_identity_mol=math.fsum(
            (
                changes.dry_hexane_mol[0],
                -dry_volume * dry_cell.total_hexane_concentration_mol_m3,
            )
        ),
        newborn_energy_inventory_identity_j=math.fsum(
            (
                changes.dry_energy_j[0],
                -dry_volume * dry_cell.energy_density_j_m3,
            )
        ),
        historical_wet_hexane_label_change=max(
            (
                abs(new - old)
                for new, old in zip(
                    before.historical_hexane_loadings,
                    before.historical_hexane_loadings,
                )
            ),
            default=0.0,
        ),
        geometry_cell_gcl_residual_m3=swept.max_cell_gcl_residual_m3,
        minimum_entropy_production_w_m3_k=min(entropy_values, default=0.0),
    )


def datum_covariant_residual_vector(
    residuals: cut.CutResidualBlocks,
) -> tuple[float, ...]:
    """Return the equation-equivalent ``R_E-H_REF*R_h`` energy basis."""

    dry_energy = tuple(
        math.fsum((energy, -hx.H_REF * hexane))
        for energy, hexane in zip(
            residuals.dry_energy_w,
            residuals.dry_hexane_mol_s,
        )
    )
    rh_energy = math.fsum(
        (residuals.rh_energy_w, -hx.H_REF * residuals.rh_hexane_mol_s)
    )
    return (
        *residuals.wet_water_mol_s,
        *residuals.wet_energy_w,
        *residuals.dry_water_mol_s,
        *residuals.dry_hexane_mol_s,
        *dry_energy,
        residuals.rh_water_mol_s,
        residuals.rh_hexane_mol_s,
        rh_energy,
    )


def classify_endpoint_branch(
    *,
    outward_hexane_flux_mol_m2_s: float,
    hexane_storage_jump_mol_m3: float,
    receding: RecedingBranchEvidence | None = None,
    contact: FullyWetContactEvidence | None = None,
    numerical_failure: str | None = None,
) -> BirthBranchDecision:
    """Apply the revision-2 exact-sign endpoint active set without pinning.

    ``numerical_failure`` takes precedence and can never select contact.  Flux
    sign is exact; no physical deadband is introduced.  Positive flux uses the
    frozen zero-volume relation ``v=N_h/Delta C_h`` and requires a separately
    certified strict-partial solve.  Exact zero requires independently
    certified fully-wet contact evidence.  Negative flux is outside the
    qualified monotonic primary-drainage envelope.
    """

    values = (outward_hexane_flux_mol_m2_s, hexane_storage_jump_mol_m3)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("endpoint flux and storage jump must remain finite")
    if hexane_storage_jump_mol_m3 <= 0.0:
        raise BirthTopologyError("birth requires a positive n-hexane storage jump")
    if numerical_failure is not None:
        if not isinstance(numerical_failure, str) or not numerical_failure.strip():
            raise ValueError("numerical failure needs a non-empty classification")
        return BirthBranchDecision(
            BirthBranch.REJECTED_NUMERICAL_FAILURE,
            outward_hexane_flux_mol_m2_s,
            hexane_storage_jump_mol_m3,
            None,
            f"rollback/reject; contact forbidden after {numerical_failure}",
            False,
        )
    if outward_hexane_flux_mol_m2_s > 0.0:
        speed = outward_hexane_flux_mol_m2_s / hexane_storage_jump_mol_m3
        if receding is None or not receding.certified:
            return BirthBranchDecision(
                BirthBranch.REJECTED_MISSING_BRANCH_SOLUTION,
                outward_hexane_flux_mol_m2_s,
                hexane_storage_jump_mol_m3,
                speed,
                "positive drive requires a converged simultaneous receding birth solve",
                False,
            )
        return BirthBranchDecision(
            BirthBranch.RECEDING,
            outward_hexane_flux_mol_m2_s,
            hexane_storage_jump_mol_m3,
            speed,
            "positive outward n-hexane flux requires immediate recession",
            False,
        )
    if outward_hexane_flux_mol_m2_s < 0.0:
        return BirthBranchDecision(
            BirthBranch.INWARD_OUTSIDE_ENVELOPE,
            outward_hexane_flux_mol_m2_s,
            hexane_storage_jump_mol_m3,
            None,
            "inward coincident-endpoint tendency is outside primary drainage",
            False,
        )
    if contact is None or not contact.certified:
        return BirthBranchDecision(
            BirthBranch.REJECTED_MISSING_BRANCH_SOLUTION,
            0.0,
            hexane_storage_jump_mol_m3,
            None,
            "exact-zero transfer still requires an independently converged contact solve",
            False,
        )
    return BirthBranchDecision(
        BirthBranch.CONTACT,
        0.0,
        hexane_storage_jump_mol_m3,
        0.0,
        "independent fully-wet exact-zero-transfer contact solve certified",
        True,
    )


def _require_length(name: str, values: Sequence[float], expected: int) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} has length {len(values)}, expected {expected}")


__all__ = [
    "BirthAssembly",
    "BirthBranch",
    "BirthBranchDecision",
    "BirthInventoryChanges",
    "BirthLayout",
    "BirthLedger",
    "BirthStepError",
    "BirthTopologyError",
    "BirthUnknowns",
    "FullyWetContactEvidence",
    "RecedingBranchEvidence",
    "assemble_birth_backward_euler",
    "classify_endpoint_branch",
    "datum_covariant_residual_vector",
    "layout_for_birth",
    "zero_old_dry_inventories",
]
