"""Bounded pre-acceptance K=2 native particle/cell source replacement.

This module is the first executable joint-source slice.  It starts from one
accepted K=2, single-wall tray state and solves the gas, pressure, shared wall,
binary film, and tag-1 packet endpoints in one trial.  The resulting native
packet endpoints are then presented once to the existing weighted host.

The accepted source-step path admits one homogeneous,
positive-attached-hexane/positive-external-water tag-1 group per K layer.  A
separate exact-``+0.0`` retained-upper path is assessment-only: it returns a
typed pre-event NO-GO with counterfactual boundary evidence and never enters
the weighted host.  One boundary assessment carries one explicit
``HexaneAccessOracle`` model-form case.  Transparent rows impose both
pure-liquid fugacity equalities.  Impermeable-water-film rows impose the water
equality and exact ``N_h = +0.0``; no inaccessible-hexane fugacity row is
added.  The access cases are evaluated independently from the same prior and
never share or duplicate physical layer area.

The accepted Coletto/FSG coefficient chain uses the local hydraulic pressure.
Native particle liquid calorics retain the packet caloric pressure.  Native
pressure-virial partial enthalpies close phase and gas energy.  Component datum
offsets are added to both gas enthalpy and packet liquid internal energy and
therefore cancel from the interface first law.

The legacy aggregate packet transfer, allocation, and accepted-flux callback
are not evaluated.  A private endpoint feeder exists solely because the
current weighted-host API consumes precomputed packet entries through that
protocol.  The host is called exactly once only on the accepted positive-water
source path and zero times on the absent-water assessment.  This bounded
result is engineering pre-acceptance evidence only; every qualification,
production, and plant-prediction claim remains false.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
from typing import ClassVar

import numpy as np
from scipy import optimize

from . import cell_closure as closure
from . import cell_engineering_feasibility as engineering
from . import cell_native_caloric as native_caloric
from . import high_loading_packet_adapter as high_loading
from . import k_cell_tray_host as tray_host

# W3-1 (owner ruling 2026-09-03, executed at HEAD f5a70bb): the shared
# telescoped-pressure-row resolution-bound construction lives beside the
# module that owns `_binary64_resolution` / `_residual_beyond_resolution`.
from . import sp1_k2_law2_tray_integration as legacy
from . import through_bed_k2_single_wall_tray_integration as single_wall
from . import through_bed_k_single_wall_thermohydraulic as arbitrary_k
from .particle import surface_active_set
from .particle.surface_active_set import HexaneAccessOracle
from .props import binary_gas
from .props import hexane
from .props import water
from .tray_particle_stateful import CellExternalTransfer, VerticalLayerId


QSC_K2_NATIVE_SOURCE_SCHEMA_ID = "dtdc-core2-qsc-k2-native-source-kernel-v1"
QSC_K2_NATIVE_SOURCE_SOLVER_ID = "scipy-phase-aware-nested-seed-scaled-joint-k2-native-source-v1"
QSC_K2_NATIVE_ENDPOINT_FEEDER_ID = "qsc-k2-native-preaccepted-endpoints-v1"
QSC_K2_NATIVE_SOURCE_ID = "native-virial-film-high-loading-k2-single-wall-v1"


class QSCK2NativeSourceError(ValueError):
    """The bounded native-source request or nonlinear trial is inadmissible."""


class K2NativeSourceEndpointEventKind(str, Enum):
    """First endpoint event that needs a different particle/film branch."""

    ATTACHED_HEXANE_DISAPPEARANCE = "attached_hexane_disappearance"
    EXTERNAL_WATER_DISAPPEARANCE = "external_water_disappearance"
    SHARED_FILM_AREA_CLOSURE = "shared_film_area_closure"


class _K2NativeWaterConditionAssessmentKind(str, Enum):
    """Exact-sign counterfactual assessment of absent external water."""

    ABSENT_CONDITIONS_RECORDED = "absent_conditions_recorded"


class K2NativeWaterBoundaryFluxSign(str, Enum):
    """Exact sign of outward-positive water flux at the W=0 boundary."""

    CONDENSATION = "condensation"
    EXACT_ZERO = "exact_zero"
    EVAPORATION = "evaporation"


class K2NativeWaterAppearanceRefusalKind(str, Enum):
    """Why the classifier packet cannot admit a positive-duration successor."""

    MIXED_LAYER_ACTIVE_SET_ORDERING = "mixed_layer_active_set_ordering"
    NONCONDENSING_LAYER = "noncondensing_layer"
    BOUNDARY_NONLINEAR_GATE = "boundary_nonlinear_gate"
    POSITIVE_SUCCESSOR_NOT_ADMITTED = "positive_successor_not_admitted"


class K2NativeSourceEndpointFirstEventRefusal(QSCK2NativeSourceError):
    """The requested endpoint reaches or crosses the current FilmAreaLaw domain."""

    def __init__(
        self,
        message: str,
        *,
        event_kind: K2NativeSourceEndpointEventKind,
        group_id: str,
        accepted_time_s: float,
        endpoint_time_s: float,
        endpoint_attached_hexane_saturation: float,
        endpoint_external_water_saturation: float,
    ) -> None:
        super().__init__(message)
        self.event_kind = event_kind
        self.group_id = group_id
        self.accepted_time_s = accepted_time_s
        self.endpoint_time_s = endpoint_time_s
        self.endpoint_attached_hexane_saturation = endpoint_attached_hexane_saturation
        self.endpoint_external_water_saturation = endpoint_external_water_saturation
        self.clamped = False
        self.post_event_branch_followed = False


@dataclass(frozen=True, slots=True)
class K2NativeWaterAppearanceCondition:
    """Native accepted-``t_n`` fugacity-sign evidence for one K layer."""

    group_id: str
    vertical_layer_id: VerticalLayerId
    packet_ids: tuple[str, ...]
    accepted_external_water_kg: float
    interface_temperature_k: float
    interface_hexane_mole_fraction: float
    local_hydraulic_pressure_pa: float
    phase_condition: surface_active_set.ComponentPhaseCondition

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class _K2NativeWaterConditionAssessment:
    """Access-neutral condition evidence; this object is not a transition event."""

    assessment_kind: _K2NativeWaterConditionAssessmentKind
    attempt_id: str
    accepted_state: single_wall.AcceptedThroughBedK2SingleWallTrayState = field(
        repr=False,
        compare=False,
    )
    accepted_state_digest: str
    accepted_tray_state_digest: str
    tray_step_definition_digest: str
    access_neutral_request_digest: str
    accepted_time_s: float
    event_time_s: float
    requested_end_time_s: float
    conditions: tuple[K2NativeWaterAppearanceCondition, ...]

    host_evaluation_count: ClassVar[int] = 0
    elapsed_time_s: ClassVar[float] = 0.0
    mass_change_kg: ClassVar[float] = 0.0
    energy_change_j: ClassVar[float] = 0.0
    k_rekey_count: ClassVar[int] = 0
    rtd_transition_count: ClassVar[int] = 0
    access_oracle_selected: ClassVar[bool] = False
    epsilon_inventory_used: ClassVar[bool] = False
    nucleation_law_used: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False

    @property
    def assessment_identity(self) -> str:
        payload = repr(
            (
                self.assessment_kind.value,
                self.attempt_id,
                self.accepted_state_digest,
                self.accepted_tray_state_digest,
                self.tray_step_definition_digest,
                self.access_neutral_request_digest,
                self.accepted_time_s,
                self.event_time_s,
                self.requested_end_time_s,
                self.conditions,
            )
        ).encode("utf-8")
        return "sha256:" + hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True, slots=True)
class K2NativeWaterBoundaryLayerEvidence:
    """One layer's zero-water boundary flux and backward-Euler direction."""

    group_id: str
    vertical_layer_id: VerticalLayerId
    hexane_access_oracle: HexaneAccessOracle
    water_molar_flux_mol_m2_s: float
    active_area_m2: float
    step_duration_s: float
    water_molar_mass_kg_mol: float
    physics_predictor_weighted_water_kg: float
    flux_sign: K2NativeWaterBoundaryFluxSign
    sign_is_stable_under_adjacent_binary64: bool

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


@dataclass(frozen=True, slots=True)
class K2NativeWaterAppearanceBoundaryRefusal:
    """Typed, host-free boundary result for an inadmissible successor ordering."""

    kind: K2NativeWaterAppearanceRefusalKind
    access_neutral_assessment_identity: str
    hexane_access_oracle: HexaneAccessOracle
    layer_evidence: tuple[
        K2NativeWaterBoundaryLayerEvidence,
        K2NativeWaterBoundaryLayerEvidence,
    ]
    maximum_scaled_residual: float
    nonlinear_limit: float
    boundary_nonlinear_gate_passed: bool

    host_evaluation_count: ClassVar[int] = 0
    positive_duration_accepted: ClassVar[bool] = False
    epsilon_inventory_used: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False


class K2NativeWaterPreEventNoGoKind(str, Enum):
    """First missing authority before a retained-upper transition may be issued."""

    MISSING_ACCEPTED_NATIVE_HEXANE_ONLY_LAYER_STATE = (
        "missing_accepted_native_hexane_only_layer_state"
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class K2NativeWaterPreEventNoGo:
    """Host-free refusal: no canonical GOV02 gas/wall/cell state exists to transition."""

    kind: K2NativeWaterPreEventNoGoKind
    attempt_id: str
    rollback_state: single_wall.AcceptedThroughBedK2SingleWallTrayState = field(repr=False)
    accepted_state_digest: str
    accepted_tray_state_digest: str
    accepted_time_s: float
    access_neutral_assessment_identity: str
    counterfactual_phase_conditions: tuple[K2NativeWaterAppearanceCondition, ...]
    counterfactual_boundary_refusal: K2NativeWaterAppearanceBoundaryRefusal
    missing_authority: str

    host_evaluation_count: ClassVar[int] = 0
    event_issued: ClassVar[bool] = False
    positive_duration_accepted: ClassVar[bool] = False
    elapsed_time_s: ClassVar[float] = 0.0
    mass_change_kg: ClassVar[float] = 0.0
    energy_change_j: ClassVar[float] = 0.0
    k_rekey_count: ClassVar[int] = 0
    rtd_transition_count: ClassVar[int] = 0
    epsilon_inventory_used: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False


def _finite(label: str, value: float, *, positive: bool = False) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise QSCK2NativeSourceError(f"{label} must be a finite exact binary64")
    if positive and value <= 0.0:
        raise QSCK2NativeSourceError(f"{label} must be strictly positive")


def _positive_zero(value: float) -> bool:
    return type(value) is float and value == 0.0 and math.copysign(1.0, value) == 1.0


def _relative(residual: float, *terms: float) -> float:
    return abs(residual) / max(*(abs(value) for value in terms), 1.0e-300)


def _resolution(*values: float) -> float:
    return math.fsum(math.ulp(abs(value)) for value in values)


def _beyond_resolution(residual: float, bound: float) -> float:
    return max(abs(residual) - bound, 0.0)


@dataclass(frozen=True, slots=True, kw_only=True)
class K2NativeSourceGroup:
    """One explicit homogeneous tag-1 source group in one K layer."""

    group_id: str
    vertical_layer_id: VerticalLayerId
    packet_ids: tuple[str, ...]
    hexane_access_oracle: HexaneAccessOracle

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.group_id) is not str or not self.group_id.strip():
            raise QSCK2NativeSourceError("native source group ID must be nonblank")
        if type(self.vertical_layer_id) is not VerticalLayerId or (
            self.vertical_layer_id.value not in (1, 2)
        ):
            raise QSCK2NativeSourceError("native source group must name K1 or K2")
        if (
            type(self.packet_ids) is not tuple
            or not self.packet_ids
            or any(type(value) is not str or not value.strip() for value in self.packet_ids)
        ):
            raise QSCK2NativeSourceError("native source group requires packet IDs")
        if self.packet_ids != tuple(sorted(self.packet_ids)) or len(set(self.packet_ids)) != len(
            self.packet_ids
        ):
            raise QSCK2NativeSourceError(
                "native source group packet IDs must be unique and canonically ordered"
            )
        if self.hexane_access_oracle not in (
            HexaneAccessOracle.TRANSPARENT_INTERFACE,
            HexaneAccessOracle.IMPERMEABLE_WATER_FILM,
        ):
            raise QSCK2NativeSourceError(
                "native source group requires an explicit transparent or impermeable endpoint"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCK2NativeSourceRequest:
    """One immutable K=2 trial definition."""

    tray_step: single_wall.ThroughBedK2SingleWallStepRequest
    groups: tuple[K2NativeSourceGroup, ...]
    maximum_function_evaluations: int = 4_000

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.tray_step) is not single_wall.ThroughBedK2SingleWallStepRequest:
            raise QSCK2NativeSourceError("native source requires an exact single-wall K2 step")
        if (
            type(self.groups) is not tuple
            or not self.groups
            or any(type(group) is not K2NativeSourceGroup for group in self.groups)
        ):
            raise QSCK2NativeSourceError("native source requires exact homogeneous groups")
        expected_order = tuple(
            sorted(self.groups, key=lambda item: (item.vertical_layer_id.value, item.group_id))
        )
        if self.groups != expected_order:
            raise QSCK2NativeSourceError("native source groups must be in canonical K/group order")
        group_ids = tuple(group.group_id for group in self.groups)
        packet_ids = tuple(packet_id for group in self.groups for packet_id in group.packet_ids)
        if len(set(group_ids)) != len(group_ids) or len(set(packet_ids)) != len(packet_ids):
            raise QSCK2NativeSourceError("native source group or packet membership is duplicated")
        if len(self.groups) != 2 or tuple(
            group.vertical_layer_id.value for group in self.groups
        ) != (1, 2):
            raise QSCK2NativeSourceError(
                "one native source solve requires exactly one full-area group per K layer"
            )
        if len({group.hexane_access_oracle for group in self.groups}) != 1:
            raise QSCK2NativeSourceError(
                "transparent and impermeable access are separate model-form solves"
            )
        if (
            type(self.maximum_function_evaluations) is not int
            or self.maximum_function_evaluations <= 0
        ):
            raise QSCK2NativeSourceError("native source solve budget must be a positive integer")
        if self.tray_step.meal_mechanical_work_deposition is not None:
            raise QSCK2NativeSourceError(
                "the first native source slice does not yet admit packet mechanical deposition"
            )


@dataclass(frozen=True, slots=True)
class NativeFilmCoefficients:
    accepted_total_mass_flow_kg_s: float
    gas_density_kg_m3: float
    superficial_velocity_m_s: float
    binary_diffusivity_m2_s: float
    reynolds_voidage: float
    prandtl_number: float
    schmidt_number: float
    nusselt_voidage: float
    heat_transfer_coefficient_w_m2_k: float
    mass_transfer_coefficient_m_s: float
    molar_transport_coefficient_mol_m2_s: float


@dataclass(frozen=True, slots=True)
class NativeSourceGroupEndpoint:
    group_id: str
    vertical_layer_id: VerticalLayerId
    packet_ids: tuple[str, ...]
    hexane_access_oracle: HexaneAccessOracle
    dry_mass_fraction: float
    representative_particles: float
    accepted_t_n_attached_hexane_saturation: float
    accepted_t_n_external_water_saturation: float
    endpoint_attached_hexane_saturation: float
    endpoint_external_water_saturation: float
    accepted_group_bed_subvolume_m3: float
    hexane_inventory_capacity_kg: float
    water_inventory_capacity_kg: float
    active_area_m2: float
    interface_temperature_k: float
    interface_hexane_mole_fraction: float
    total_molar_flux_mol_m2_s: float
    hexane_molar_flux_mol_m2_s: float
    water_molar_flux_mol_m2_s: float
    raw_impermeable_hexane_flux_mol_m2_s: float
    gas_to_interface_w: float
    wall_to_interface_w: float
    particle_to_interface_w: float
    native_hexane_partial_enthalpy_j_mol: float
    native_water_partial_enthalpy_j_mol: float
    native_hexane_liquid_internal_energy_j_mol: float
    native_water_liquid_internal_energy_j_mol: float
    hexane_fugacity_residual_pa: float | None
    water_fugacity_residual_pa: float
    interface_energy_residual_w: float
    interface_energy_resolution_bound_w: float
    endpoint_payload_identity: str
    endpoint_attached_hexane_kg: float
    endpoint_external_water_kg: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    aggregate_transfer_allocation_used: ClassVar[bool] = False
    endpoint_inventory_basis: ClassVar[str] = "weighted-accepted-group"
    active_area_is_lagged_from_accepted_t_n: ClassVar[bool] = True
    endpoint_area_feedback_used_within_step: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        access_passed = (
            self.hexane_fugacity_residual_pa is not None
            if self.hexane_access_oracle is HexaneAccessOracle.TRANSPARENT_INTERFACE
            else _positive_zero(self.hexane_molar_flux_mol_m2_s)
        )
        return (
            access_passed
            and self.endpoint_attached_hexane_kg > 0.0
            and self.endpoint_external_water_kg > 0.0
            and self.endpoint_attached_hexane_saturation > 0.0
            and self.endpoint_external_water_saturation > 0.0
            and math.fsum(
                (
                    self.endpoint_attached_hexane_saturation,
                    self.endpoint_external_water_saturation,
                )
            )
            < 1.0
        )


@dataclass(frozen=True, slots=True)
class NativeSourceLayerResult:
    vertical_layer_id: VerticalLayerId
    inlet_molar_flow_mol_s: float
    inlet_hexane_mole_fraction: float
    inlet_temperature_k: float
    inlet_face_pressure_pa: float
    outlet_molar_flow_mol_s: float
    outlet_hexane_mole_fraction: float
    outlet_temperature_k: float
    outlet_face_pressure_pa: float
    local_hydraulic_pressure_pa: float
    pressure_drop_pa: float
    component_hexane_residual_mol_s: float
    independent_water_residual_mol_s: float
    gas_energy_residual_w: float
    pressure_residual_pa: float
    wall_to_gas_w: float
    coefficients: NativeFilmCoefficients
    external_cell_transfer: CellExternalTransfer

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


@dataclass(frozen=True, slots=True)
class NativeSourceConservationLedger:
    layer_hexane_flux_kg: tuple[float, float]
    layer_weighted_packet_hexane_decrement_kg: tuple[float, float]
    layer_water_flux_kg: tuple[float, float]
    layer_weighted_packet_water_decrement_kg: tuple[float, float]
    layer_particle_energy_flux_j: tuple[float, float]
    layer_weighted_packet_energy_decrement_j: tuple[float, float]
    cell_water_storage_residual_kg: tuple[float, float]
    cell_hexane_storage_residual_kg: tuple[float, float]
    cell_energy_storage_residual_j: tuple[float, float]
    expected_cell_energy_storage_j: tuple[float, float]
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    aggregate_transfer_allocation_used: ClassVar[bool] = False
    particle_delta_counted_from_before_after_only: ClassVar[bool] = True

    @property
    def passed(self) -> bool:
        return self.maximum_relative_residual <= self.relative_limit


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedQSCK2NativeSourceStep:
    request: QSCK2NativeSourceRequest
    host_transaction: tray_host.ValidatedKCellTrayStep
    layer_results: tuple[NativeSourceLayerResult, NativeSourceLayerResult]
    group_endpoints: tuple[NativeSourceGroupEndpoint, ...]
    shared_wall_temperature_k: float
    shared_wall_residual_w: float
    shared_wall_resolution_bound_w: float
    maximum_scaled_residual: float
    scipy_function_evaluations: int
    aggregate_pressure: arbitrary_k.AggregatePhysicalTrayPressureLedger
    conservation_ledger: NativeSourceConservationLedger
    endpoint_entries: tuple[tray_host.PacketPlanEntry, ...] = field(repr=False)

    schema_id: ClassVar[str] = QSC_K2_NATIVE_SOURCE_SCHEMA_ID
    host_evaluation_count: ClassVar[int] = 1
    aggregate_packet_exchange_used: ClassVar[bool] = False
    aggregate_packet_allocation_used: ClassVar[bool] = False
    legacy_accepted_flux_callback_used: ClassVar[bool] = False
    precomputed_endpoint_feeder_used: ClassVar[bool] = True
    fixed_packet_caloric_pressure_retained: ClassVar[bool] = True
    local_hydraulic_pressure_used_for_film_and_fugacity: ClassVar[bool] = True
    native_virial_partial_enthalpy_used: ClassVar[bool] = True
    physical_wall_count: ClassVar[int] = 1
    coupled_positive_outlet_flux_chart_used: ClassVar[bool] = True
    static_rectangular_total_flux_bounds_used: ClassVar[bool] = False
    physical_total_molar_flux_preserved_in_evidence: ClassVar[bool] = True
    local_per_cell_pressure_gate_used: ClassVar[bool] = False
    aggregate_physical_tray_pressure_gate_used: ClassVar[bool] = True
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    industrially_feasible: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCK2NativeSourceRejection:
    attempt_id: str
    reason: str
    rollback_state: single_wall.AcceptedThroughBedK2SingleWallTrayState
    endpoint_first_event: K2NativeSourceEndpointFirstEventRefusal | None = field(
        default=None,
        repr=False,
    )
    water_appearance_boundary_refusal: K2NativeWaterAppearanceBoundaryRefusal | None = field(
        default=None,
        repr=False,
    )

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


@dataclass(frozen=True, slots=True)
class _GroupBasis:
    request: K2NativeSourceGroup
    packets: tuple[tray_host.AcceptedTrayPacket, ...]
    complete: high_loading.CompleteHighLoadingPacket
    dry_mass_fraction: float
    representative_particles: float
    attached_hexane_saturation: float
    external_water_saturation: float
    accepted_group_bed_subvolume_m3: float
    hexane_inventory_capacity_kg: float
    water_inventory_capacity_kg: float
    active_area_m2: float


@dataclass(frozen=True, slots=True)
class _GroupTrial:
    basis: _GroupBasis
    interface_temperature_k: float
    interface_hexane_mole_fraction: float
    total_molar_flux_mol_m2_s: float
    raw_hexane_molar_flux_mol_m2_s: float
    hexane_molar_flux_mol_m2_s: float
    water_molar_flux_mol_m2_s: float
    gas_to_interface_w: float
    wall_to_interface_w: float
    particle_to_interface_w: float
    native_hexane_partial_enthalpy_j_mol: float
    native_water_partial_enthalpy_j_mol: float
    common_hexane_partial_enthalpy_j_mol: float
    common_water_partial_enthalpy_j_mol: float
    native_hexane_liquid_internal_energy_j_mol: float
    native_water_liquid_internal_energy_j_mol: float
    common_hexane_liquid_internal_energy_j_mol: float
    common_water_liquid_internal_energy_j_mol: float
    hexane_fugacity_residual_pa: float | None
    water_fugacity_residual_pa: float
    access_residual_scaled: float
    water_fugacity_residual_scaled: float
    interface_energy_residual_w: float
    interface_energy_residual_scaled: float
    interface_energy_resolution_bound_w: float


@dataclass(frozen=True, slots=True)
class _LayerTrial:
    result: NativeSourceLayerResult
    groups: tuple[_GroupTrial, ...]
    component_hexane_residual_scaled: float
    gas_energy_residual_scaled: float
    pressure_residual_scaled: float
    phase_energy_common_w: float
    particle_energy_loss_w: float
    wall_to_interface_w: float


@dataclass(frozen=True, slots=True)
class _JointTrial:
    layers: tuple[_LayerTrial, _LayerTrial]
    shared_wall_temperature_k: float
    shared_wall_residual_w: float
    shared_wall_residual_scaled: float
    shared_wall_resolution_bound_w: float

    @property
    def scaled_residuals(self) -> tuple[float, ...]:
        rows: list[float] = []
        for layer in self.layers:
            for group in layer.groups:
                rows.extend(
                    (
                        group.access_residual_scaled,
                        group.water_fugacity_residual_scaled,
                        group.interface_energy_residual_scaled,
                    )
                )
            rows.extend(
                (
                    layer.component_hexane_residual_scaled,
                    layer.gas_energy_residual_scaled,
                    layer.pressure_residual_scaled,
                )
            )
        rows.append(self.shared_wall_residual_scaled)
        return tuple(rows)


@dataclass(frozen=True, slots=True)
class _EndpointFeeder:
    entries: dict[str, tray_host.PacketPlanEntry]
    maximum_scaled_residual: float
    callback_id: str = QSC_K2_NATIVE_ENDPOINT_FEEDER_ID
    source_identity: str = QSC_K2_NATIVE_SOURCE_ID

    physically_qualifying: ClassVar[bool] = False
    aggregate_packet_exchange_used: ClassVar[bool] = False

    def propose(
        self, request: tray_host.ManufacturedPacketAdvanceRequest
    ) -> tray_host.ManufacturedPacketAdvanceProposal:
        return tray_host.ManufacturedPacketAdvanceProposal(
            packet_id=request.packet.packet_id,
            entry_after=self.entries[request.packet.packet_id],
            callback_id=self.callback_id,
            source_identity=self.source_identity,
            residual_contract_passed=True,
            maximum_scaled_residual=self.maximum_scaled_residual,
        )


@lru_cache(maxsize=4_096)
def _pure_liquid_fugacity_pa(module: object, temperature_k: float, pressure_pa: float) -> float:
    # Evaluate the pure compressed liquid at the actual local hydraulic
    # pressure.  The explicit critical-density check prevents either property
    # package from silently returning its vapor root.  This is the existing
    # heteroazeotrope/PD convention and avoids a numerically fragile auxiliary
    # saturation solve at every finite-difference point.
    density = module.density(temperature_k, pressure_pa, "liquid")
    if not density > module.RHOC:
        raise QSCK2NativeSourceError(
            f"{module.__name__.rsplit('.', 1)[-1]} pure-liquid root collapsed"
        )
    value = math.exp(module.ln_fugacity(temperature_k, density))
    _finite("pure-liquid fugacity", value, positive=True)
    return value


def _partial_enthalpies(
    temperature_k: float,
    pressure_pa: float,
    hexane_mole_fraction: float,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> tuple[float, float, float, float, binary_gas.BinaryVirialState]:
    gas = binary_gas.state(
        temperature_k,
        pressure_pa,
        1.0 - hexane_mole_fraction,
        hexane_mole_fraction,
    )
    native_hexane = math.fsum(
        (hexane.h_ideal(temperature_k), gas.partial_residual_enthalpy_hexane_molar)
    )
    native_water = math.fsum(
        (
            binary_gas.water_ideal_enthalpy_molar(temperature_k),
            gas.partial_residual_enthalpy_water_molar,
        )
    )
    datum = codec.component_datum_adapter
    return (
        native_hexane,
        native_water,
        native_hexane + datum.hexane_native_to_common_offset_j_mol,
        native_water + datum.water_native_to_common_offset_j_mol,
        gas,
    )


def _bulk_common_enthalpy_j_mol(
    temperature_k: float,
    pressure_pa: float,
    hexane_mole_fraction: float,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    _, _, hexane_common, water_common, _ = _partial_enthalpies(
        temperature_k,
        pressure_pa,
        hexane_mole_fraction,
        codec=codec,
    )
    return math.fsum(
        (
            hexane_mole_fraction * hexane_common,
            (1.0 - hexane_mole_fraction) * water_common,
        )
    )


def _film_coefficients(
    model: single_wall.ThroughBedKLayerEngineeringModel,
    *,
    outlet_molar_flow_mol_s: float,
    bulk_temperature_k: float,
    bulk_hexane_mole_fraction: float,
    local_hydraulic_pressure_pa: float,
) -> NativeFilmCoefficients:
    transport = model.transport
    bed_void = model.film_area_law.bed_void_fraction
    mixture_molar_mass = math.fsum(
        (
            bulk_hexane_mole_fraction * hexane.M,
            (1.0 - bulk_hexane_mole_fraction) * water.M,
        )
    )
    density = (
        local_hydraulic_pressure_pa
        * mixture_molar_mass
        / (engineering.UNIVERSAL_GAS_CONSTANT_J_MOL_K * bulk_temperature_k)
    )
    total_mass_flow = outlet_molar_flow_mol_s * mixture_molar_mass
    velocity = total_mass_flow / (density * model.geometry.cross_section_m2)
    diffusivity = (
        transport.fsg_diffusivity_multiplier
        * engineering.FSG_DIFFUSIVITY_PREFACTOR
        * bulk_temperature_k**1.75
        / local_hydraulic_pressure_pa
    )
    reynolds = (
        density
        * velocity
        * transport.hydraulic_particle_diameter_m
        / (transport.gas_dynamic_viscosity_pa_s * bed_void)
    )
    prandtl = (
        transport.gas_mass_heat_capacity_j_kg_k
        * transport.gas_dynamic_viscosity_pa_s
        / transport.gas_thermal_conductivity_w_m_k
    )
    schmidt = transport.gas_dynamic_viscosity_pa_s / (density * diffusivity)
    nusselt = 0.6949 * reynolds**0.579 * prandtl ** (1.0 / 3.0)
    heat = (
        nusselt
        * transport.gas_thermal_conductivity_w_m_k
        * (1.0 - bed_void)
        / (transport.hydraulic_particle_diameter_m * bed_void)
    )
    mass = (
        heat
        / (density * transport.gas_mass_heat_capacity_j_kg_k)
        * (prandtl / schmidt) ** (2.0 / 3.0)
    )
    kappa = (
        local_hydraulic_pressure_pa
        / (engineering.UNIVERSAL_GAS_CONSTANT_J_MOL_K * bulk_temperature_k)
        * mass
    )
    values = (
        total_mass_flow,
        density,
        velocity,
        diffusivity,
        reynolds,
        prandtl,
        schmidt,
        nusselt,
        heat,
        mass,
        kappa,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise QSCK2NativeSourceError("native Coletto/FSG coefficient chain left its domain")
    return NativeFilmCoefficients(*values)


def _group_bases(
    request: QSCK2NativeSourceRequest,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> tuple[_GroupBasis, ...]:
    tray = request.tray_step.prior.tray
    by_id = {packet.packet_id: packet for packet in tray.packets}
    declared_ids = tuple(packet_id for group in request.groups for packet_id in group.packet_ids)
    accepted_ids = tuple(packet.packet_id for packet in tray.packets)
    if set(declared_ids) != set(accepted_ids) or len(declared_ids) != len(accepted_ids):
        raise QSCK2NativeSourceError(
            "native source groups must cover every accepted packet exactly once"
        )
    total_dry_by_layer = {
        layer: math.fsum(
            packet.weighted_totals.totals[0]
            for packet in tray.packets
            if packet.owner_key.vertical_layer_id.value == layer
        )
        for layer in (1, 2)
    }
    bases: list[_GroupBasis] = []
    for group in request.groups:
        layer = group.vertical_layer_id.value
        packets = tuple(by_id[packet_id] for packet_id in group.packet_ids)
        if any(packet.owner_key.vertical_layer_id.value != layer for packet in packets):
            raise QSCK2NativeSourceError(f"group {group.group_id} crosses a K-layer boundary")
        decoded = tuple(codec.decode(packet.entry.payload_bytes) for packet in packets)
        first = decoded[0]
        if any(item.payload_bytes != first.payload_bytes for item in decoded[1:]):
            raise QSCK2NativeSourceError(
                f"group {group.group_id} is not an exact homogeneous tag-1 group"
            )
        if first.state.attached_hexane_mass_kg <= 0.0:
            raise QSCK2NativeSourceError(f"group {group.group_id} lacks positive attached hexane")
        if first.surface.free_water_kg == 0.0 and not _positive_zero(first.surface.free_water_kg):
            raise QSCK2NativeSourceError(
                f"group {group.group_id} external water must use canonical +0.0"
            )
        representatives = math.fsum(
            packet.entry.weight.representative_particles for packet in packets
        )
        weighted_dry = math.fsum(packet.weighted_totals.totals[0] for packet in packets)
        dry_fraction = weighted_dry / total_dry_by_layer[layer]
        model = request.tray_step.layer_models[layer - 1]
        law = model.film_area_law
        subvolume = dry_fraction * model.geometry.bed_volume_m3
        hexane_capacity = law.hexane_inventory_capacity_kg(subvolume)
        water_capacity = law.water_inventory_capacity_kg(subvolume)
        attached = math.fsum(packet.weighted_totals.totals[2] for packet in packets)
        external_water = math.fsum(packet.weighted_totals.totals[4] for packet in packets)
        saturation_h = attached / hexane_capacity
        saturation_w = external_water / water_capacity
        water_admissible = saturation_w > 0.0 or _positive_zero(saturation_w)
        if saturation_h <= 0.0 or not water_admissible or saturation_h + saturation_w >= 1.0:
            raise QSCK2NativeSourceError(
                f"group {group.group_id} has an inadmissible positive-H/nonnegative-W saturation"
            )
        active_area = (
            dry_fraction
            * model.geometry.gas_side_reference_area_m2
            * (1.0 - saturation_h - saturation_w) ** model.film_area_law.exponent.value
        )
        _finite("native group active area", active_area, positive=True)
        bases.append(
            _GroupBasis(
                request=group,
                packets=packets,
                complete=first,
                dry_mass_fraction=dry_fraction,
                representative_particles=representatives,
                attached_hexane_saturation=saturation_h,
                external_water_saturation=saturation_w,
                accepted_group_bed_subvolume_m3=subvolume,
                hexane_inventory_capacity_kg=hexane_capacity,
                water_inventory_capacity_kg=water_capacity,
                active_area_m2=active_area,
            )
        )
    for layer in (1, 2):
        fraction_sum = math.fsum(
            basis.dry_mass_fraction
            for basis in bases
            if basis.request.vertical_layer_id.value == layer
        )
        if fraction_sum != 1.0:
            raise QSCK2NativeSourceError(
                f"K{layer} homogeneous group dry fractions do not close exactly"
            )
    return tuple(bases)


def _access_neutral_request_digest(request: QSCK2NativeSourceRequest) -> str:
    payload = repr(
        (
            request.tray_step.definition_digest,
            tuple(
                (
                    group.group_id,
                    group.vertical_layer_id.value,
                    group.packet_ids,
                )
                for group in request.groups
            ),
            request.maximum_function_evaluations,
        )
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def _water_condition_assessment(
    request: QSCK2NativeSourceRequest,
    bases: tuple[_GroupBasis, ...],
) -> _K2NativeWaterConditionAssessment:
    if not all(_positive_zero(basis.external_water_saturation) for basis in bases):
        raise QSCK2NativeSourceError(
            "water appearance classification requires exact accepted +0.0 external water"
        )
    fields = request.tray_step.prior.cell_field_states
    conditions: list[K2NativeWaterAppearanceCondition] = []
    for basis in bases:
        field_state = fields[basis.request.vertical_layer_id.value - 1]
        gas = binary_gas.state(
            field_state.interface_temperature_k,
            field_state.layer_pressure_pa,
            1.0 - field_state.interface_hexane_mole_fraction,
            field_state.interface_hexane_mole_fraction,
        )
        liquid_fugacity = _pure_liquid_fugacity_pa(
            water,
            field_state.interface_temperature_k,
            field_state.layer_pressure_pa,
        )
        phase_condition = surface_active_set.classify_phase_condition(
            surface_active_set.Component.WATER,
            0.0,
            gas.fugacity_water_pa,
            liquid_fugacity,
        )
        conditions.append(
            K2NativeWaterAppearanceCondition(
                group_id=basis.request.group_id,
                vertical_layer_id=basis.request.vertical_layer_id,
                packet_ids=basis.request.packet_ids,
                accepted_external_water_kg=0.0,
                interface_temperature_k=field_state.interface_temperature_k,
                interface_hexane_mole_fraction=field_state.interface_hexane_mole_fraction,
                local_hydraulic_pressure_pa=field_state.layer_pressure_pa,
                phase_condition=phase_condition,
            )
        )
    prior = request.tray_step.prior
    accepted_time = prior.tray.time_s
    return _K2NativeWaterConditionAssessment(
        assessment_kind=(_K2NativeWaterConditionAssessmentKind.ABSENT_CONDITIONS_RECORDED),
        attempt_id=request.tray_step.attempt_id,
        accepted_state=prior,
        accepted_state_digest=prior.state_digest,
        accepted_tray_state_digest=prior.tray.state_digest,
        tray_step_definition_digest=request.tray_step.definition_digest,
        access_neutral_request_digest=_access_neutral_request_digest(request),
        accepted_time_s=accepted_time,
        event_time_s=accepted_time,
        requested_end_time_s=request.tray_step.end_time_s,
        conditions=tuple(conditions),
    )


def _group_trial(
    basis: _GroupBasis,
    *,
    model: single_wall.ThroughBedKLayerEngineeringModel,
    coefficients: NativeFilmCoefficients,
    bulk_temperature_k: float,
    bulk_hexane_mole_fraction: float,
    local_hydraulic_pressure_pa: float,
    shared_wall_temperature_k: float,
    interface_temperature_k: float,
    interface_hexane_mole_fraction: float,
    total_molar_flux_mol_m2_s: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> _GroupTrial:
    kappa = coefficients.molar_transport_coefficient_mol_m2_s
    mass_bernoulli, _ = engineering.bernoulli_pair(float(total_molar_flux_mol_m2_s / kappa))
    raw_hexane_flux = math.fsum(
        (
            interface_hexane_mole_fraction * total_molar_flux_mol_m2_s,
            kappa * mass_bernoulli * (interface_hexane_mole_fraction - bulk_hexane_mole_fraction),
        )
    )
    if basis.request.hexane_access_oracle is HexaneAccessOracle.IMPERMEABLE_WATER_FILM:
        hexane_flux = 0.0
    else:
        hexane_flux = raw_hexane_flux
    water_flux = total_molar_flux_mol_m2_s - hexane_flux
    (
        native_h_h,
        native_h_w,
        common_h_h,
        common_h_w,
        interface_gas,
    ) = _partial_enthalpies(
        interface_temperature_k,
        local_hydraulic_pressure_pa,
        interface_hexane_mole_fraction,
        codec=codec,
    )
    particle_temperature = basis.complete.state.temperature_k
    caloric_pressure = basis.packets[0].pressure_pa
    liquid_h_native = (
        hexane.state_Tp(particle_temperature, caloric_pressure, "liquid").u_mass * hexane.M
    )
    liquid_w_native = (
        water.state_Tp(particle_temperature, caloric_pressure, "liquid").u_mass * water.M
    )
    datum = codec.component_datum_adapter
    liquid_h_common = liquid_h_native + datum.hexane_native_to_common_offset_j_mol
    liquid_w_common = liquid_w_native + datum.water_native_to_common_offset_j_mol
    area = basis.active_area_m2
    # QB-16 (ruled 2026-08-24): the Ackermann/Colburn-Drew blowing correction
    # is the ONLY cp consumer on this path, and its argument N*cp/h is
    # dimensionless only if cp is the MOLAR heat capacity of the film gas
    # mixture.  The declared arm reads the declared row exactly as before -
    # the 40.0 row stays issued, labelled and canonical, so this is not the
    # E-7 canonical-pack flip.  The native arm has no such row (CELL-02d
    # NativeGasProperties carries no numbers), so the same frozen law the
    # certified Law-2 native arm uses supplies it, at the evaluation state
    # three landed precedents fix rather than leave free: the interface
    # triple this function already hands to `_partial_enthalpies` above.
    # A plain float evaluation is derivative-consistent here because this
    # kernel passes scipy no `jac=`; the Jacobian is finite-differenced.
    if type(model.properties) is closure.DeclaredGasProperties:
        ackermann_heat_capacity_j_mol_k = model.properties.molar_heat_capacity_j_mol_k
    else:
        ackermann_heat_capacity_j_mol_k = native_caloric.mixture_heat_capacity(
            interface_temperature_k,
            local_hydraulic_pressure_pa,
            interface_hexane_mole_fraction,
        ).value_j_mol_k
    heat_bernoulli, _ = engineering.bernoulli_pair(
        float(
            total_molar_flux_mol_m2_s
            * ackermann_heat_capacity_j_mol_k
            / coefficients.heat_transfer_coefficient_w_m2_k
        )
    )
    gas_heat = (
        coefficients.heat_transfer_coefficient_w_m2_k
        * area
        * (bulk_temperature_k - interface_temperature_k)
        * heat_bernoulli
    )
    wall_heat = (
        basis.dry_mass_fraction
        * model.transfer.wall_to_interface_ua_w_k
        * (shared_wall_temperature_k - interface_temperature_k)
    )
    particle_heat = (
        basis.dry_mass_fraction
        * model.transfer.solid_to_interface_ua_w_k
        * (particle_temperature - interface_temperature_k)
    )
    phase_energy = area * math.fsum(
        (
            hexane_flux * (native_h_h - liquid_h_native),
            water_flux * (native_h_w - liquid_w_native),
        )
    )
    interface_energy = math.fsum((gas_heat, wall_heat, particle_heat, -phase_energy))
    interface_scale = max(abs(gas_heat), abs(wall_heat), abs(particle_heat), abs(phase_energy), 1.0)
    water_fugacity = interface_gas.fugacity_water_pa - _pure_liquid_fugacity_pa(
        water,
        interface_temperature_k,
        local_hydraulic_pressure_pa,
    )
    if basis.request.hexane_access_oracle is HexaneAccessOracle.TRANSPARENT_INTERFACE:
        hexane_fugacity: float | None = interface_gas.fugacity_hexane_pa - _pure_liquid_fugacity_pa(
            hexane,
            interface_temperature_k,
            local_hydraulic_pressure_pa,
        )
        access_scaled = hexane_fugacity / local_hydraulic_pressure_pa
    else:
        hexane_fugacity = None
        access_scaled = raw_hexane_flux / max(kappa, 1.0e-300)
    energy_resolution = _resolution(
        gas_heat,
        wall_heat,
        particle_heat,
        phase_energy,
        interface_energy,
    )
    return _GroupTrial(
        basis=basis,
        interface_temperature_k=interface_temperature_k,
        interface_hexane_mole_fraction=interface_hexane_mole_fraction,
        total_molar_flux_mol_m2_s=total_molar_flux_mol_m2_s,
        raw_hexane_molar_flux_mol_m2_s=raw_hexane_flux,
        hexane_molar_flux_mol_m2_s=hexane_flux,
        water_molar_flux_mol_m2_s=water_flux,
        gas_to_interface_w=gas_heat,
        wall_to_interface_w=wall_heat,
        particle_to_interface_w=particle_heat,
        native_hexane_partial_enthalpy_j_mol=native_h_h,
        native_water_partial_enthalpy_j_mol=native_h_w,
        common_hexane_partial_enthalpy_j_mol=common_h_h,
        common_water_partial_enthalpy_j_mol=common_h_w,
        native_hexane_liquid_internal_energy_j_mol=liquid_h_native,
        native_water_liquid_internal_energy_j_mol=liquid_w_native,
        common_hexane_liquid_internal_energy_j_mol=liquid_h_common,
        common_water_liquid_internal_energy_j_mol=liquid_w_common,
        hexane_fugacity_residual_pa=hexane_fugacity,
        water_fugacity_residual_pa=water_fugacity,
        access_residual_scaled=access_scaled,
        water_fugacity_residual_scaled=(water_fugacity / local_hydraulic_pressure_pa),
        interface_energy_residual_w=interface_energy,
        interface_energy_residual_scaled=(interface_energy / interface_scale),
        interface_energy_resolution_bound_w=energy_resolution,
    )


def _positive_outlet_interval_mol_s(
    *,
    inlet_flow_mol_s: float,
    active_area_m2: float,
) -> tuple[float, float, float]:
    law_lower, law_upper = engineering.LAW2_TOTAL_FLUX_DOMAIN_MOL_M2_S
    law_lower_outlet = math.fsum((inlet_flow_mol_s, active_area_m2 * law_lower))
    lower_flux = law_lower
    if law_lower_outlet <= 0.0:
        lower_flux = max(law_lower, -inlet_flow_mol_s / active_area_m2)
        lower_outlet = math.fsum((inlet_flow_mol_s, active_area_m2 * lower_flux))
        while lower_outlet <= 0.0:
            next_flux = math.nextafter(lower_flux, math.inf)
            if next_flux == lower_flux or next_flux > law_upper:
                raise QSCK2NativeSourceError(
                    "Law-2 signed-flux domain has no binary64 point with positive outlet"
                )
            lower_flux = next_flux
            lower_outlet = math.fsum((inlet_flow_mol_s, active_area_m2 * lower_flux))
    else:
        lower_outlet = law_lower_outlet
    law_upper_outlet = math.fsum((inlet_flow_mol_s, active_area_m2 * law_upper))
    if not math.isfinite(law_upper_outlet) or law_upper_outlet <= lower_outlet:
        raise QSCK2NativeSourceError(
            "Law-2 signed-flux domain has no strictly positive outlet intersection"
        )
    return lower_flux, lower_outlet, law_upper_outlet


def _decode_positive_outlet_flux_chart(
    chart_vector: tuple[float, ...],
    *,
    request: QSCK2NativeSourceRequest,
    bases: tuple[_GroupBasis, ...],
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    physical = list(chart_vector)
    inlet_flow = request.tray_step.bottom_gas_boundary.inlet_molar_flow_mol_s
    outlet_flows: list[float] = []
    law_lower, law_upper = engineering.LAW2_TOTAL_FLUX_DOMAIN_MOL_M2_S
    for index, basis in enumerate(bases):
        coordinate_index = 9 + 3 * index
        coordinate = chart_vector[coordinate_index]
        if not 0.0 <= coordinate <= 1.0:
            raise QSCK2NativeSourceError(
                f"group {basis.request.group_id} positive-outlet chart coordinate left [0, 1]"
            )
        lower_flux, lower_outlet, upper_outlet = _positive_outlet_interval_mol_s(
            inlet_flow_mol_s=inlet_flow,
            active_area_m2=basis.active_area_m2,
        )
        if coordinate == 0.0:
            total_flux = lower_flux
        elif coordinate == 1.0:
            total_flux = law_upper
        else:
            total_flux = math.fsum((lower_flux, coordinate * (law_upper - lower_flux)))
        outlet_flow = math.fsum((inlet_flow, basis.active_area_m2 * total_flux))
        if (
            outlet_flow <= 0.0
            or not lower_flux <= total_flux <= law_upper
            or not lower_outlet <= outlet_flow <= upper_outlet
        ):
            raise QSCK2NativeSourceError(
                f"group {basis.request.group_id} chart decoded outside its signed-flux domain"
            )
        physical[coordinate_index] = total_flux
        outlet_flows.append(outlet_flow)
        inlet_flow = outlet_flow
    return tuple(physical), tuple(outlet_flows)


def _encode_positive_outlet_flux_chart(
    physical_vector: np.ndarray,
    *,
    request: QSCK2NativeSourceRequest,
    bases: tuple[_GroupBasis, ...],
) -> np.ndarray:
    chart = physical_vector.copy()
    inlet_flow = request.tray_step.bottom_gas_boundary.inlet_molar_flow_mol_s
    law_lower, law_upper = engineering.LAW2_TOTAL_FLUX_DOMAIN_MOL_M2_S
    for index, basis in enumerate(bases):
        coordinate_index = 9 + 3 * index
        total_flux = float(physical_vector[coordinate_index])
        if not law_lower <= total_flux <= law_upper:
            raise QSCK2NativeSourceError(
                f"group {basis.request.group_id} seed left the full Law-2 flux domain"
            )
        outlet_flow = math.fsum((inlet_flow, basis.active_area_m2 * total_flux))
        if outlet_flow <= 0.0:
            raise QSCK2NativeSourceError(
                f"group {basis.request.group_id} seed has nonpositive sequential outlet flow"
            )
        lower_flux, _, _ = _positive_outlet_interval_mol_s(
            inlet_flow_mol_s=inlet_flow,
            active_area_m2=basis.active_area_m2,
        )
        if total_flux == lower_flux:
            coordinate = 0.0
        elif total_flux == law_upper:
            coordinate = 1.0
        else:
            coordinate = (total_flux - lower_flux) / (law_upper - lower_flux)
        if not 0.0 <= coordinate <= 1.0:
            raise QSCK2NativeSourceError(
                f"group {basis.request.group_id} seed cannot enter the positive-outlet chart"
            )
        chart[coordinate_index] = coordinate
        inlet_flow = outlet_flow
    return chart


def _positive_outlet_chart_bounds(
    physical_lower: np.ndarray,
    physical_upper: np.ndarray,
    *,
    bases: tuple[_GroupBasis, ...],
) -> tuple[np.ndarray, np.ndarray]:
    lower = physical_lower.copy()
    upper = physical_upper.copy()
    for index, _ in enumerate(bases):
        coordinate_index = 9 + 3 * index
        lower[coordinate_index] = 0.0
        upper[coordinate_index] = 1.0
    return lower, upper


def _joint_trial(
    vector: tuple[float, ...],
    *,
    request: QSCK2NativeSourceRequest,
    bases: tuple[_GroupBasis, ...],
    codec: high_loading.HighLoadingV1PacketAdapter,
    positive_outlet_flux_chart: bool = False,
) -> _JointTrial:
    values = tuple(float(value) for value in vector)
    expected_length = 7 + 3 * len(bases)
    if len(values) != expected_length or any(not math.isfinite(value) for value in values):
        raise QSCK2NativeSourceError("joint native-source vector is malformed")
    charted_outlet_flows: tuple[float, ...] | None = None
    if positive_outlet_flux_chart:
        values, charted_outlet_flows = _decode_positive_outlet_flux_chart(
            values,
            request=request,
            bases=bases,
        )
    bulk_values = ((values[0], values[1], values[2]), (values[3], values[4], values[5]))
    wall_temperature = values[6]
    group_vectors = tuple(
        (values[index], values[index + 1], values[index + 2]) for index in range(7, len(values), 3)
    )
    step = request.tray_step
    inlet_flow = step.bottom_gas_boundary.inlet_molar_flow_mol_s
    inlet_hexane = step.bottom_gas_boundary.inlet_hexane_mole_fraction
    inlet_temperature = step.bottom_gas_boundary.inlet_temperature_k
    inlet_face_pressure = bulk_values[0][2]
    layers: list[_LayerTrial] = []
    for layer in (1, 2):
        bulk_temperature, bulk_hexane, hydraulic_pressure = bulk_values[layer - 1]
        model = step.layer_models[layer - 1]
        selected = tuple(
            (basis, group_vector)
            for basis, group_vector in zip(bases, group_vectors, strict=True)
            if basis.request.vertical_layer_id.value == layer
        )
        outlet_flow = (
            charted_outlet_flows[layer - 1]
            if charted_outlet_flows is not None
            else math.fsum(
                (
                    inlet_flow,
                    math.fsum(
                        basis.active_area_m2 * group_vector[2] for basis, group_vector in selected
                    ),
                )
            )
        )
        if outlet_flow <= 0.0:
            raise QSCK2NativeSourceError("joint native-source trial has nonpositive outlet flow")
        coefficients = _film_coefficients(
            model,
            outlet_molar_flow_mol_s=outlet_flow,
            bulk_temperature_k=bulk_temperature,
            bulk_hexane_mole_fraction=bulk_hexane,
            local_hydraulic_pressure_pa=hydraulic_pressure,
        )
        group_trials = tuple(
            _group_trial(
                basis,
                model=model,
                coefficients=coefficients,
                bulk_temperature_k=bulk_temperature,
                bulk_hexane_mole_fraction=bulk_hexane,
                local_hydraulic_pressure_pa=hydraulic_pressure,
                shared_wall_temperature_k=wall_temperature,
                interface_temperature_k=group_vector[0],
                interface_hexane_mole_fraction=group_vector[1],
                total_molar_flux_mol_m2_s=group_vector[2],
                codec=codec,
            )
            for basis, group_vector in selected
        )
        hexane_source = math.fsum(
            group.basis.active_area_m2 * group.hexane_molar_flux_mol_m2_s for group in group_trials
        )
        water_source = math.fsum(
            group.basis.active_area_m2 * group.water_molar_flux_mol_m2_s for group in group_trials
        )
        hexane_residual = math.fsum(
            (
                outlet_flow * bulk_hexane,
                -inlet_flow * inlet_hexane,
                -hexane_source,
            )
        )
        water_residual = math.fsum(
            (
                outlet_flow * (1.0 - bulk_hexane),
                -inlet_flow * (1.0 - inlet_hexane),
                -water_source,
            )
        )
        outlet_face_pressure = bulk_values[1][2] if layer == 1 else step.top_boundary_pressure_pa
        inlet_enthalpy = _bulk_common_enthalpy_j_mol(
            inlet_temperature,
            inlet_face_pressure,
            inlet_hexane,
            codec=codec,
        )
        outlet_enthalpy = _bulk_common_enthalpy_j_mol(
            bulk_temperature,
            outlet_face_pressure,
            bulk_hexane,
            codec=codec,
        )
        phase_energy_common = math.fsum(
            group.basis.active_area_m2
            * math.fsum(
                (
                    group.hexane_molar_flux_mol_m2_s * group.common_hexane_partial_enthalpy_j_mol,
                    group.water_molar_flux_mol_m2_s * group.common_water_partial_enthalpy_j_mol,
                )
            )
            for group in group_trials
        )
        particle_energy_loss = math.fsum(
            (
                group.particle_to_interface_w
                + group.basis.active_area_m2
                * math.fsum(
                    (
                        group.hexane_molar_flux_mol_m2_s
                        * group.common_hexane_liquid_internal_energy_j_mol,
                        group.water_molar_flux_mol_m2_s
                        * group.common_water_liquid_internal_energy_j_mol,
                    )
                )
            )
            for group in group_trials
        )
        wall_to_gas = model.transfer.wall_to_gas_ua_w_k * (wall_temperature - bulk_temperature)
        gas_heat_to_interface = math.fsum(group.gas_to_interface_w for group in group_trials)
        gas_energy_residual = math.fsum(
            (
                outlet_flow * outlet_enthalpy,
                -inlet_flow * inlet_enthalpy,
                -phase_energy_common,
                -wall_to_gas,
                gas_heat_to_interface,
            )
        )
        pressure_drop = closure.layer_pressure_drop(
            series=model.hydraulics.series,
            superficial_velocity_m_s=coefficients.superficial_velocity_m_s,
            gas_density_kg_m3=coefficients.gas_density_kg_m3,
        ).drop_pa
        pressure_residual = hydraulic_pressure - outlet_face_pressure - pressure_drop
        hexane_scale = max(abs(inlet_flow), abs(outlet_flow), abs(hexane_source), 1.0e-12)
        energy_scale = max(
            abs(outlet_flow * outlet_enthalpy),
            abs(inlet_flow * inlet_enthalpy),
            abs(phase_energy_common),
            abs(wall_to_gas),
            abs(gas_heat_to_interface),
            1.0,
        )
        q_steam = model.wall.steam_side_ua_w_k * (model.wall.steam_temperature_k - wall_temperature)
        q_ambient = model.wall.ambient_ua_w_k * (
            wall_temperature - model.wall.ambient_temperature_k
        )
        dt = step.end_time_s - step.prior.tray.time_s
        external = CellExternalTransfer(
            water_to_cell_kg=(
                dt
                * math.fsum(
                    (
                        inlet_flow * (1.0 - inlet_hexane) * water.M,
                        -outlet_flow * (1.0 - bulk_hexane) * water.M,
                    )
                )
            ),
            hexane_to_cell_kg=(
                dt
                * math.fsum(
                    (
                        inlet_flow * inlet_hexane * hexane.M,
                        -outlet_flow * bulk_hexane * hexane.M,
                    )
                )
            ),
            common_datum_energy_to_cell_j=(
                dt
                * math.fsum(
                    (
                        inlet_flow * inlet_enthalpy,
                        -outlet_flow * outlet_enthalpy,
                        q_steam,
                        -q_ambient,
                    )
                )
            ),
            energy_datum_id=step.prior.tray.energy_datum_id,
        )
        result = NativeSourceLayerResult(
            vertical_layer_id=VerticalLayerId(layer),
            inlet_molar_flow_mol_s=inlet_flow,
            inlet_hexane_mole_fraction=inlet_hexane,
            inlet_temperature_k=inlet_temperature,
            inlet_face_pressure_pa=inlet_face_pressure,
            outlet_molar_flow_mol_s=outlet_flow,
            outlet_hexane_mole_fraction=bulk_hexane,
            outlet_temperature_k=bulk_temperature,
            outlet_face_pressure_pa=outlet_face_pressure,
            local_hydraulic_pressure_pa=hydraulic_pressure,
            pressure_drop_pa=pressure_drop,
            component_hexane_residual_mol_s=hexane_residual,
            independent_water_residual_mol_s=water_residual,
            gas_energy_residual_w=gas_energy_residual,
            pressure_residual_pa=pressure_residual,
            wall_to_gas_w=wall_to_gas,
            coefficients=coefficients,
            external_cell_transfer=external,
        )
        layers.append(
            _LayerTrial(
                result=result,
                groups=group_trials,
                component_hexane_residual_scaled=(hexane_residual / hexane_scale),
                gas_energy_residual_scaled=(gas_energy_residual / energy_scale),
                pressure_residual_scaled=(pressure_residual / 10_000.0),
                phase_energy_common_w=phase_energy_common,
                particle_energy_loss_w=particle_energy_loss,
                wall_to_interface_w=math.fsum(group.wall_to_interface_w for group in group_trials),
            )
        )
        inlet_flow = outlet_flow
        inlet_hexane = bulk_hexane
        inlet_temperature = bulk_temperature
        inlet_face_pressure = outlet_face_pressure
    lower_wall = step.layer_models[0].wall
    upper_wall = step.layer_models[1].wall
    if any(
        value != 0.0
        for value in (
            upper_wall.steam_side_ua_w_k,
            upper_wall.ambient_ua_w_k,
            step.layer_models[1].transfer.wall_to_gas_ua_w_k,
            step.layer_models[1].transfer.wall_to_interface_ua_w_k,
        )
    ):
        raise QSCK2NativeSourceError("native K2 trial requires the existing one-wall ownership")
    dt = step.end_time_s - step.prior.tray.time_s
    storage_w = (
        lower_wall.thermal_capacity_j_k * (wall_temperature - step.prior.wall_temperature_k) / dt
    )
    steam_w = lower_wall.steam_side_ua_w_k * (lower_wall.steam_temperature_k - wall_temperature)
    ambient_w = lower_wall.ambient_ua_w_k * (wall_temperature - lower_wall.ambient_temperature_k)
    wall_to_gas_w = math.fsum(layer.result.wall_to_gas_w for layer in layers)
    wall_to_interface_w = math.fsum(layer.wall_to_interface_w for layer in layers)
    wall_residual = math.fsum((storage_w, -steam_w, ambient_w, wall_to_gas_w, wall_to_interface_w))
    wall_scale = max(
        abs(storage_w),
        abs(steam_w),
        abs(ambient_w),
        abs(wall_to_gas_w),
        abs(wall_to_interface_w),
        1.0,
    )
    wall_resolution = _resolution(
        storage_w,
        steam_w,
        ambient_w,
        wall_to_gas_w,
        wall_to_interface_w,
        wall_residual,
    )
    return _JointTrial(
        layers=(layers[0], layers[1]),
        shared_wall_temperature_k=wall_temperature,
        shared_wall_residual_w=wall_residual,
        shared_wall_residual_scaled=(wall_residual / wall_scale),
        shared_wall_resolution_bound_w=wall_resolution,
    )


def _aggregate_physical_tray_pressure_ledger(
    request: QSCK2NativeSourceRequest,
    trial: _JointTrial,
) -> arbitrary_k.AggregatePhysicalTrayPressureLedger:
    layer_results = tuple(layer.result for layer in trial.layers)
    drops = tuple(
        closure.layer_pressure_drop(
            series=model.hydraulics.series,
            superficial_velocity_m_s=result.coefficients.superficial_velocity_m_s,
            gas_density_kg_m3=result.coefficients.gas_density_kg_m3,
        )
        for model, result in zip(
            request.tray_step.layer_models,
            layer_results,
            strict=True,
        )
    )
    bed = tuple(drop.contribution_for(closure.BedErgunElement.label).drop_pa for drop in drops)
    floor = tuple(
        drop.contribution_for(closure.TaperedBoreFloorPassageElement.label).drop_pa
        for drop in drops
    )
    ordered = math.fsum(
        value
        for bed_value, floor_value in zip(bed, floor, strict=True)
        for value in (bed_value, floor_value)
    )
    cell_sum = math.fsum(result.pressure_drop_pa for result in layer_results)
    endpoint = (
        layer_results[0].local_hydraulic_pressure_pa - request.tray_step.top_boundary_pressure_pa
    )
    cell_residuals = tuple(
        drop.drop_pa - result.pressure_drop_pa
        for drop, result in zip(drops, layer_results, strict=True)
    )
    aggregate_residual = endpoint - ordered
    # W3-1 (owner ruling 2026-09-03, executed at HEAD f5a70bb): gate the folded
    # pressure residuals on the Newton pressure rows' own guaranteed resolution
    # (K rows x newton_tolerance x 1.0e4 Pa row scale) instead of on a fixed
    # relative literal divided by the tray drop.  The joint K=2 native trial
    # carries the internal face as a joint unknown, so no face-Picard slip
    # enters the telescoped sum.  `relative_limit` is untouched.
    resolution_bound = legacy.aggregate_pressure_binary64_resolution_bound(
        newton_tolerance=request.tray_step.newton_tolerance,
        pressure_row_count=len(layer_results),
        summed_terms_pa=(*bed, *floor, ordered, cell_sum, endpoint),
    )
    maximum = max(
        *(
            _relative(
                legacy._residual_beyond_resolution(residual, resolution_bound),
                drop.drop_pa,
                result.pressure_drop_pa,
            )
            for residual, drop, result in zip(
                cell_residuals,
                drops,
                layer_results,
                strict=True,
            )
        ),
        _relative(
            legacy._residual_beyond_resolution(ordered - cell_sum, resolution_bound),
            ordered,
            cell_sum,
        ),
        _relative(
            legacy._residual_beyond_resolution(aggregate_residual, resolution_bound),
            endpoint,
            ordered,
        ),
    )
    low, high = arbitrary_k.AGGREGATE_PHYSICAL_TRAY_DROP_HULL_PA
    ledger = arbitrary_k.AggregatePhysicalTrayPressureLedger(
        bed_drop_by_layer_pa=bed,
        floor_drop_by_layer_pa=floor,
        active_physical_floor_count=sum(value > 0.0 for value in floor),
        ordered_element_sum_pa=ordered,
        cell_drop_sum_pa=cell_sum,
        bottom_to_top_pressure_drop_pa=endpoint,
        cell_series_residuals_pa=cell_residuals,
        aggregate_pressure_residual_pa=aggregate_residual,
        aggregate_gate_band_pa=(low, high),
        inside_aggregate_gate=low <= ordered <= high,
        maximum_relative_residual=maximum,
        relative_limit=request.tray_step.relative_limit,
        pressure_authority_id=(arbitrary_k.AGGREGATE_PHYSICAL_TRAY_PRESSURE_AUTHORITY_ID),
    )
    if not ledger.passed:
        raise QSCK2NativeSourceError(
            "aggregate physical-tray bed-plus-one-floor pressure gate did not pass"
        )
    return ledger


def _seed_and_bounds(
    request: QSCK2NativeSourceRequest,
    bases: tuple[_GroupBasis, ...],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    prior = request.tray_step.prior
    lower, upper = prior.cell_field_states
    seed: list[float] = [
        lower.gas_temperature_k,
        lower.gas_hexane_mole_fraction,
        lower.layer_pressure_pa,
        upper.gas_temperature_k,
        upper.gas_hexane_mole_fraction,
        upper.layer_pressure_pa,
        prior.wall_temperature_k,
    ]
    lower_bounds: list[float] = [
        engineering.LAW2_TEMPERATURE_DOMAIN_K[0],
        0.0,
        engineering.LAW2_PRESSURE_DOMAIN_PA[0],
        engineering.LAW2_TEMPERATURE_DOMAIN_K[0],
        0.0,
        engineering.LAW2_PRESSURE_DOMAIN_PA[0],
        250.0,
    ]
    upper_bounds: list[float] = [
        engineering.LAW2_TEMPERATURE_DOMAIN_K[1],
        1.0,
        engineering.LAW2_PRESSURE_DOMAIN_PA[1],
        engineering.LAW2_TEMPERATURE_DOMAIN_K[1],
        1.0,
        engineering.LAW2_PRESSURE_DOMAIN_PA[1],
        500.0,
    ]
    by_layer = {1: lower, 2: upper}
    # These are physical Law-2 limits for the phase-aware scalar seeds, not
    # rectangular solver bounds.  The nonlinear solve uses the coupled
    # sequential positive-outlet chart defined above.
    for basis in bases:
        prior_field = by_layer[basis.request.vertical_layer_id.value]
        seed.extend(
            (
                prior_field.interface_temperature_k,
                prior_field.interface_hexane_mole_fraction,
                prior_field.total_molar_flux_mol_m2_s,
            )
        )
        lower_bounds.extend(
            (
                engineering.LAW2_TEMPERATURE_DOMAIN_K[0],
                0.0,
                engineering.LAW2_TOTAL_FLUX_DOMAIN_MOL_M2_S[0],
            )
        )
        upper_bounds.extend(
            (
                engineering.LAW2_TEMPERATURE_DOMAIN_K[1],
                1.0,
                engineering.LAW2_TOTAL_FLUX_DOMAIN_MOL_M2_S[1],
            )
        )
    seed_array = np.asarray(seed, dtype=float)
    lower_array = np.asarray(lower_bounds, dtype=float)
    upper_array = np.asarray(upper_bounds, dtype=float)
    seed_array = np.minimum(np.maximum(seed_array, lower_array), upper_array)
    return seed_array, lower_array, upper_array


def _signed_flux_grid(
    *,
    inlet_flow_mol_s: float,
    area_m2: float,
    lower_bound: float,
    upper_bound: float,
    continuation_flux: float,
) -> tuple[float, ...]:
    dynamic_lower = max(lower_bound, -inlet_flow_mol_s / area_m2)
    while math.fsum((inlet_flow_mol_s, area_m2 * dynamic_lower)) <= 0.0:
        next_lower = math.nextafter(dynamic_lower, math.inf)
        if next_lower == dynamic_lower or next_lower > upper_bound:
            raise QSCK2NativeSourceError(
                "signed-flux seed has no binary64 point with positive outlet flow"
            )
        dynamic_lower = next_lower
    candidates = [
        dynamic_lower,
        0.999 * dynamic_lower,
        0.9 * dynamic_lower,
        0.75 * dynamic_lower,
        0.5 * dynamic_lower,
        0.25 * dynamic_lower,
        0.1 * dynamic_lower,
        continuation_flux,
        0.0,
        1.0e-5,
        1.0e-4,
        1.0e-3,
        3.0e-3,
        1.0e-2,
        3.0e-2,
        1.0e-1,
        3.0e-1,
        1.0,
        3.0,
        upper_bound,
    ]
    return tuple(
        sorted({min(max(float(value), dynamic_lower), upper_bound) for value in candidates})
    )


def _bracketed_scalar_roots(
    function: object,
    points: tuple[float, ...],
) -> tuple[float, ...]:
    values = tuple(float(function(point)) for point in points)  # type: ignore[operator]
    if any(not math.isfinite(value) for value in values):
        raise QSCK2NativeSourceError("phase-aware scalar seed left its finite domain")
    roots: list[float] = []
    for left, right, value_left, value_right in zip(
        points[:-1],
        points[1:],
        values[:-1],
        values[1:],
        strict=True,
    ):
        if value_left == 0.0:
            candidate = left
        elif value_left * value_right < 0.0:
            candidate = float(
                optimize.brentq(
                    function,  # type: ignore[arg-type]
                    left,
                    right,
                    xtol=1.0e-13,
                    rtol=1.0e-13,
                )
            )
        else:
            continue
        if not roots or abs(candidate - roots[-1]) > 1.0e-12:
            roots.append(candidate)
    if values[-1] == 0.0 and (not roots or abs(points[-1] - roots[-1]) > 1.0e-12):
        roots.append(points[-1])
    return tuple(roots)


def _transparent_phase_root(
    *,
    pressure_pa: float,
    continuation: tuple[float, float],
) -> tuple[float, float]:
    def residual(vector: np.ndarray) -> np.ndarray:
        temperature, hexane_fraction = (float(vector[0]), float(vector[1]))
        gas = binary_gas.state(
            temperature,
            pressure_pa,
            1.0 - hexane_fraction,
            hexane_fraction,
        )
        return np.asarray(
            (
                (
                    gas.fugacity_hexane_pa
                    - _pure_liquid_fugacity_pa(hexane, temperature, pressure_pa)
                )
                / pressure_pa,
                (gas.fugacity_water_pa - _pure_liquid_fugacity_pa(water, temperature, pressure_pa))
                / pressure_pa,
            ),
            dtype=float,
        )

    solved = optimize.root(residual, np.asarray(continuation, dtype=float), tol=1.0e-11)
    values = residual(solved.x)
    temperature, hexane_fraction = (float(solved.x[0]), float(solved.x[1]))
    # Acceptance is the residual-plus-domain criterion, written NaN-safe. The
    # solver's progress flag is advisory only: when the continuation lands
    # within ULPs of the root, MINPACK reports "no progress" and would falsely
    # condemn a machine-precision root (2026-08-18 cleanup evidence: residuals
    # at 1e-15 with success=False). A diverged or NaN solve still refuses
    # through the residual and domain checks below.
    residual_magnitude = max(abs(value) for value in values)
    if (
        not residual_magnitude <= 1.0e-9
        or not engineering.LAW2_TEMPERATURE_DOMAIN_K[0]
        <= temperature
        <= engineering.LAW2_TEMPERATURE_DOMAIN_K[1]
        or not 0.0 <= hexane_fraction <= 1.0
    ):
        raise QSCK2NativeSourceError(
            "transparent two-fugacity phase seed did not produce an in-domain root"
        )
    return temperature, hexane_fraction


def _transparent_interface_seed(
    *,
    basis: _GroupBasis,
    model: single_wall.ThroughBedKLayerEngineeringModel,
    inlet_flow_mol_s: float,
    bulk_temperature_k: float,
    bulk_hexane_fraction: float,
    hydraulic_pressure_pa: float,
    wall_temperature_k: float,
    continuation: tuple[float, float, float],
    flux_lower_bound: float,
    flux_upper_bound: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> tuple[float, float, float]:
    interface_temperature, interface_hexane = _transparent_phase_root(
        pressure_pa=hydraulic_pressure_pa,
        continuation=(continuation[0], continuation[1]),
    )

    def energy_residual(total_flux: float) -> float:
        outlet_flow = inlet_flow_mol_s + basis.active_area_m2 * total_flux
        coefficients = _film_coefficients(
            model,
            outlet_molar_flow_mol_s=outlet_flow,
            bulk_temperature_k=bulk_temperature_k,
            bulk_hexane_mole_fraction=bulk_hexane_fraction,
            local_hydraulic_pressure_pa=hydraulic_pressure_pa,
        )
        return _group_trial(
            basis,
            model=model,
            coefficients=coefficients,
            bulk_temperature_k=bulk_temperature_k,
            bulk_hexane_mole_fraction=bulk_hexane_fraction,
            local_hydraulic_pressure_pa=hydraulic_pressure_pa,
            shared_wall_temperature_k=wall_temperature_k,
            interface_temperature_k=interface_temperature,
            interface_hexane_mole_fraction=interface_hexane,
            total_molar_flux_mol_m2_s=total_flux,
            codec=codec,
        ).interface_energy_residual_scaled

    roots = _bracketed_scalar_roots(
        energy_residual,
        _signed_flux_grid(
            inlet_flow_mol_s=inlet_flow_mol_s,
            area_m2=basis.active_area_m2,
            lower_bound=flux_lower_bound,
            upper_bound=flux_upper_bound,
            continuation_flux=continuation[2],
        ),
    )
    if not roots:
        raise QSCK2NativeSourceError(
            f"group {basis.request.group_id} transparent interface energy has no signed root"
        )
    total_flux = min(roots, key=lambda value: abs(value - continuation[2]))
    return interface_temperature, interface_hexane, total_flux


def _transparent_nested_seed(
    request: QSCK2NativeSourceRequest,
    bases: tuple[_GroupBasis, ...],
    seed: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> tuple[np.ndarray, int]:
    continuation = (
        (float(seed[7]), float(seed[8]), float(seed[9])),
        (float(seed[10]), float(seed[11]), float(seed[12])),
    )

    def interfaces(outer: np.ndarray) -> np.ndarray:
        vector = seed.copy()
        vector[:7] = outer
        inlet_flow = request.tray_step.bottom_gas_boundary.inlet_molar_flow_mol_s
        for index, basis in enumerate(bases):
            offset = 7 + 3 * index
            bulk_offset = 3 * index
            interface = _transparent_interface_seed(
                basis=basis,
                model=request.tray_step.layer_models[index],
                inlet_flow_mol_s=inlet_flow,
                bulk_temperature_k=float(outer[bulk_offset]),
                bulk_hexane_fraction=float(outer[bulk_offset + 1]),
                hydraulic_pressure_pa=float(outer[bulk_offset + 2]),
                wall_temperature_k=float(outer[6]),
                continuation=continuation[index],
                flux_lower_bound=float(lower[offset + 2]),
                flux_upper_bound=float(upper[offset + 2]),
                codec=codec,
            )
            vector[offset : offset + 3] = interface
            inlet_flow += basis.active_area_m2 * interface[2]
        return vector

    outer_rows = (3, 4, 5, 9, 10, 11, 12)

    def residual(outer: np.ndarray) -> np.ndarray:
        trial = _joint_trial(
            tuple(float(value) for value in interfaces(outer)),
            request=request,
            bases=bases,
            codec=codec,
        )
        rows = trial.scaled_residuals
        return np.asarray(tuple(rows[index] for index in outer_rows), dtype=float)

    solved = optimize.least_squares(
        residual,
        seed[:7],
        bounds=(lower[:7], upper[:7]),
        x_scale="jac",
        ftol=1.0e-13,
        xtol=1.0e-13,
        gtol=1.0e-13,
        max_nfev=request.maximum_function_evaluations,
    )
    nested_seed = interfaces(solved.x)
    maximum = max(abs(value) for value in residual(solved.x))
    seed_tolerance = max(10.0 * request.tray_step.newton_tolerance, 1.0e-9)
    if not solved.success or maximum > seed_tolerance:
        raise QSCK2NativeSourceError(
            "transparent nested bulk/wall seed did not close its seven rows: "
            f"success={solved.success}, nfev={solved.nfev}, maximum={maximum:.17g}, "
            f"seed_limit={seed_tolerance:.17g}"
        )
    return nested_seed, int(solved.nfev)


def _impermeable_interface_seed(
    *,
    basis: _GroupBasis,
    model: single_wall.ThroughBedKLayerEngineeringModel,
    inlet_flow_mol_s: float,
    bulk_temperature_k: float,
    bulk_hexane_fraction: float,
    hydraulic_pressure_pa: float,
    wall_temperature_k: float,
    continuation: tuple[float, float, float],
    flux_lower_bound: float,
    flux_upper_bound: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> tuple[float, float, float]:
    def state_at_temperature(
        interface_temperature: float,
    ) -> tuple[float, float, _GroupTrial]:
        liquid_water_fugacity = _pure_liquid_fugacity_pa(
            water,
            interface_temperature,
            hydraulic_pressure_pa,
        )

        def water_residual(interface_hexane: float) -> float:
            gas = binary_gas.state(
                interface_temperature,
                hydraulic_pressure_pa,
                1.0 - interface_hexane,
                interface_hexane,
            )
            return gas.fugacity_water_pa - liquid_water_fugacity

        at_zero = water_residual(0.0)
        at_one = water_residual(1.0)
        if at_zero == 0.0:
            interface_hexane = 0.0
        elif at_one == 0.0:
            interface_hexane = 1.0
        elif at_zero * at_one < 0.0:
            interface_hexane = float(
                optimize.brentq(water_residual, 0.0, 1.0, xtol=1.0e-13, rtol=1.0e-13)
            )
        else:
            raise QSCK2NativeSourceError("water-film phase seed has no composition root")

        def group_at_flux(total_flux: float) -> _GroupTrial:
            outlet_flow = inlet_flow_mol_s + basis.active_area_m2 * total_flux
            coefficients = _film_coefficients(
                model,
                outlet_molar_flow_mol_s=outlet_flow,
                bulk_temperature_k=bulk_temperature_k,
                bulk_hexane_mole_fraction=bulk_hexane_fraction,
                local_hydraulic_pressure_pa=hydraulic_pressure_pa,
            )
            return _group_trial(
                basis,
                model=model,
                coefficients=coefficients,
                bulk_temperature_k=bulk_temperature_k,
                bulk_hexane_mole_fraction=bulk_hexane_fraction,
                local_hydraulic_pressure_pa=hydraulic_pressure_pa,
                shared_wall_temperature_k=wall_temperature_k,
                interface_temperature_k=interface_temperature,
                interface_hexane_mole_fraction=interface_hexane,
                total_molar_flux_mol_m2_s=total_flux,
                codec=codec,
            )

        roots = _bracketed_scalar_roots(
            lambda total_flux: group_at_flux(total_flux).raw_hexane_molar_flux_mol_m2_s,
            _signed_flux_grid(
                inlet_flow_mol_s=inlet_flow_mol_s,
                area_m2=basis.active_area_m2,
                lower_bound=flux_lower_bound,
                upper_bound=flux_upper_bound,
                continuation_flux=continuation[2],
            ),
        )
        if not roots:
            raise QSCK2NativeSourceError("water-film phase seed has no zero-hexane-flux root")
        total_flux = min(roots, key=lambda value: abs(value - continuation[2]))
        return interface_hexane, total_flux, group_at_flux(total_flux)

    temperatures = tuple(
        sorted(
            {
                float(continuation[0]),
                *(
                    float(value)
                    for value in np.linspace(
                        engineering.LAW2_TEMPERATURE_DOMAIN_K[0],
                        engineering.LAW2_TEMPERATURE_DOMAIN_K[1],
                        65,
                    )
                ),
            }
        )
    )
    admissible: list[tuple[float, float]] = []
    for temperature in temperatures:
        try:
            group = state_at_temperature(temperature)[2]
        except QSCK2NativeSourceError:
            continue
        admissible.append((temperature, group.interface_energy_residual_scaled))
    brackets = tuple(
        (left[0], right[0])
        for left, right in zip(admissible[:-1], admissible[1:], strict=True)
        if left[1] == 0.0 or left[1] * right[1] < 0.0
    )
    if not brackets:
        raise QSCK2NativeSourceError("water-film interface energy has no signed phase root")
    bracket = min(
        brackets,
        key=lambda values: min(abs(values[0] - continuation[0]), abs(values[1] - continuation[0])),
    )
    interface_temperature = float(
        optimize.brentq(
            lambda temperature: (
                state_at_temperature(temperature)[2].interface_energy_residual_scaled
            ),
            bracket[0],
            bracket[1],
            xtol=1.0e-12,
            rtol=1.0e-12,
        )
    )
    interface_hexane, total_flux, _ = state_at_temperature(interface_temperature)
    return interface_temperature, interface_hexane, total_flux


def _phase_aware_seed(
    request: QSCK2NativeSourceRequest,
    bases: tuple[_GroupBasis, ...],
    seed: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> tuple[np.ndarray, int]:
    access = bases[0].request.hexane_access_oracle
    if access is HexaneAccessOracle.TRANSPARENT_INTERFACE:
        return _transparent_nested_seed(
            request,
            bases,
            seed,
            lower,
            upper,
            codec=codec,
        )
    vector = seed.copy()
    inlet_flow = request.tray_step.bottom_gas_boundary.inlet_molar_flow_mol_s
    for index, basis in enumerate(bases):
        offset = 7 + 3 * index
        bulk_offset = 3 * index
        interface = _impermeable_interface_seed(
            basis=basis,
            model=request.tray_step.layer_models[index],
            inlet_flow_mol_s=inlet_flow,
            bulk_temperature_k=float(seed[bulk_offset]),
            bulk_hexane_fraction=float(seed[bulk_offset + 1]),
            hydraulic_pressure_pa=float(seed[bulk_offset + 2]),
            wall_temperature_k=float(seed[6]),
            continuation=(
                float(seed[offset]),
                float(seed[offset + 1]),
                float(seed[offset + 2]),
            ),
            flux_lower_bound=float(lower[offset + 2]),
            flux_upper_bound=float(upper[offset + 2]),
            codec=codec,
        )
        vector[offset : offset + 3] = interface
        inlet_flow += basis.active_area_m2 * interface[2]
    return vector, 0


def _solve_joint_root(
    request: QSCK2NativeSourceRequest,
    bases: tuple[_GroupBasis, ...],
    seed: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    enforce_scaled_gate: bool = True,
) -> tuple[np.ndarray, _JointTrial, int]:
    chart_seed = _encode_positive_outlet_flux_chart(
        seed,
        request=request,
        bases=bases,
    )
    chart_lower, chart_upper = _positive_outlet_chart_bounds(
        lower,
        upper,
        bases=bases,
    )

    def residual(vector: np.ndarray) -> np.ndarray:
        trial = _joint_trial(
            tuple(float(value) for value in vector),
            request=request,
            bases=bases,
            codec=codec,
            positive_outlet_flux_chart=True,
        )
        return np.asarray(trial.scaled_residuals, dtype=float)

    solved = optimize.least_squares(
        residual,
        chart_seed,
        bounds=(chart_lower, chart_upper),
        x_scale="jac",
        ftol=1.0e-13,
        xtol=1.0e-13,
        gtol=1.0e-13,
        max_nfev=request.maximum_function_evaluations,
    )
    physical_solution, _ = _decode_positive_outlet_flux_chart(
        tuple(float(value) for value in solved.x),
        request=request,
        bases=bases,
    )
    trial = _joint_trial(
        physical_solution,
        request=request,
        bases=bases,
        codec=codec,
    )
    maximum = max(abs(value) for value in trial.scaled_residuals)
    total_nfev = int(solved.nfev)
    if solved.success and maximum > request.tray_step.newton_tolerance:
        # Escalated retry (RW-1): the owner-ruled solve-deepening pattern
        # (90b457c) applied to this joint solve.  The 1e-13 step criteria
        # can stop the solver at row floors above the acceptance limit
        # even though the same bounded system admits a deeper root
        # (proof record GT_PS2_S7C_MULTIGROUP_TERMINATION_FLOOR_PROOF_
        # 2026-08-19; the retained-upper boundary assessment terminates
        # at 1.179e-10 against the 1e-10 gate).  One deeper round each
        # from the terminated candidate and the seed, step-limited only
        # by the binary64 chart and the request's declared budget; a
        # deep trial is adopted only when it strictly improves and no
        # acceptance limit moves.  Configurations that meet the limit on
        # attempt one never reach this block.
        for start in (solved.x, chart_seed):
            deep = optimize.least_squares(
                residual,
                start,
                bounds=(chart_lower, chart_upper),
                x_scale="jac",
                ftol=None,
                xtol=3.0e-16,
                gtol=None,
                max_nfev=request.maximum_function_evaluations,
            )
            total_nfev += int(deep.nfev)
            if not deep.success:
                continue
            deep_solution, _ = _decode_positive_outlet_flux_chart(
                tuple(float(value) for value in deep.x),
                request=request,
                bases=bases,
            )
            deep_trial = _joint_trial(
                deep_solution,
                request=request,
                bases=bases,
                codec=codec,
            )
            deep_maximum = max(abs(value) for value in deep_trial.scaled_residuals)
            if deep_maximum < maximum:
                physical_solution = deep_solution
                trial = deep_trial
                maximum = deep_maximum
            if maximum <= request.tray_step.newton_tolerance:
                break
    if not solved.success or (enforce_scaled_gate and maximum > request.tray_step.newton_tolerance):
        raise QSCK2NativeSourceError(
            "joint native-source root did not satisfy the declared scaled residual limit: "
            f"success={solved.success}, nfev={total_nfev}, maximum={maximum:.17g}, "
            f"limit={request.tray_step.newton_tolerance:.17g}"
        )
    return np.asarray(physical_solution), trial, total_nfev


def _water_boundary_flux_sign(value: float) -> K2NativeWaterBoundaryFluxSign:
    if value < 0.0:
        return K2NativeWaterBoundaryFluxSign.CONDENSATION
    if value > 0.0:
        return K2NativeWaterBoundaryFluxSign.EVAPORATION
    return K2NativeWaterBoundaryFluxSign.EXACT_ZERO


def _water_appearance_boundary_refusal(
    request: QSCK2NativeSourceRequest,
    trial: _JointTrial,
    *,
    assessment_identity: str,
) -> K2NativeWaterAppearanceBoundaryRefusal:
    dt = request.tray_step.end_time_s - request.tray_step.prior.tray.time_s
    evidence: list[K2NativeWaterBoundaryLayerEvidence] = []
    for layer in trial.layers:
        group = layer.groups[0]
        flux = group.water_molar_flux_mol_m2_s
        sign = _water_boundary_flux_sign(flux)
        adjacent_signs = (
            _water_boundary_flux_sign(math.nextafter(flux, -math.inf)),
            _water_boundary_flux_sign(math.nextafter(flux, math.inf)),
        )
        evidence.append(
            K2NativeWaterBoundaryLayerEvidence(
                group_id=group.basis.request.group_id,
                vertical_layer_id=group.basis.request.vertical_layer_id,
                hexane_access_oracle=group.basis.request.hexane_access_oracle,
                water_molar_flux_mol_m2_s=flux,
                active_area_m2=group.basis.active_area_m2,
                step_duration_s=dt,
                water_molar_mass_kg_mol=water.M,
                physics_predictor_weighted_water_kg=(
                    -dt * water.M * group.basis.active_area_m2 * flux
                ),
                flux_sign=sign,
                sign_is_stable_under_adjacent_binary64=(adjacent_signs == (sign, sign)),
            )
        )
    signs = tuple(item.flux_sign for item in evidence)
    maximum = max(abs(value) for value in trial.scaled_residuals)
    if any(sign is not K2NativeWaterBoundaryFluxSign.CONDENSATION for sign in signs):
        kind = (
            K2NativeWaterAppearanceRefusalKind.MIXED_LAYER_ACTIVE_SET_ORDERING
            if len(set(signs)) > 1
            else K2NativeWaterAppearanceRefusalKind.NONCONDENSING_LAYER
        )
    elif maximum > request.tray_step.newton_tolerance:
        kind = K2NativeWaterAppearanceRefusalKind.BOUNDARY_NONLINEAR_GATE
    else:
        kind = K2NativeWaterAppearanceRefusalKind.POSITIVE_SUCCESSOR_NOT_ADMITTED
    return K2NativeWaterAppearanceBoundaryRefusal(
        kind=kind,
        access_neutral_assessment_identity=assessment_identity,
        hexane_access_oracle=request.groups[0].hexane_access_oracle,
        layer_evidence=(evidence[0], evidence[1]),
        maximum_scaled_residual=maximum,
        nonlinear_limit=request.tray_step.newton_tolerance,
        boundary_nonlinear_gate_passed=(maximum <= request.tray_step.newton_tolerance),
    )


def _endpoint_entries(
    trial: _JointTrial,
    *,
    request: QSCK2NativeSourceRequest,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> tuple[dict[str, tray_host.PacketPlanEntry], tuple[NativeSourceGroupEndpoint, ...]]:
    dt = request.tray_step.end_time_s - request.tray_step.prior.tray.time_s
    entries: dict[str, tray_host.PacketPlanEntry] = {}
    endpoints: list[NativeSourceGroupEndpoint] = []
    for layer in trial.layers:
        for group in layer.groups:
            basis = group.basis
            count = basis.representative_particles
            area = basis.active_area_m2
            hexane_decrement_per_particle = (
                dt * area * group.hexane_molar_flux_mol_m2_s * hexane.M / count
            )
            water_decrement_per_particle = (
                dt * area * group.water_molar_flux_mol_m2_s * water.M / count
            )
            energy_decrement_per_particle = (
                dt
                * math.fsum(
                    (
                        group.particle_to_interface_w,
                        area
                        * math.fsum(
                            (
                                group.hexane_molar_flux_mol_m2_s
                                * group.common_hexane_liquid_internal_energy_j_mol,
                                group.water_molar_flux_mol_m2_s
                                * group.common_water_liquid_internal_energy_j_mol,
                            )
                        ),
                    )
                )
                / count
            )
            before = basis.complete
            total_hexane_after = math.fsum(
                (before.state.total_hexane_mass_kg, -hexane_decrement_per_particle)
            )
            external_water_after = math.fsum(
                (before.surface.free_water_kg, -water_decrement_per_particle)
            )
            energy_after = math.fsum(
                (before.inventory.common_datum_energy_j, -energy_decrement_per_particle)
            )
            complete_after = codec.advance(
                before.payload_bytes,
                total_hexane_mass_after_kg=total_hexane_after,
                external_free_water_mass_after_kg=external_water_after,
                envelope_energy_after_j=energy_after,
            )
            endpoint_attached_hexane = math.fsum(
                packet.entry.weight.representative_particles
                * complete_after.state.attached_hexane_mass_kg
                for packet in basis.packets
            )
            endpoint_external_water = math.fsum(
                packet.entry.weight.representative_particles * complete_after.surface.free_water_kg
                for packet in basis.packets
            )
            endpoint_saturation_h = endpoint_attached_hexane / basis.hexane_inventory_capacity_kg
            endpoint_saturation_w = endpoint_external_water / basis.water_inventory_capacity_kg
            endpoint_saturation_sum = math.fsum((endpoint_saturation_h, endpoint_saturation_w))
            event_kind: K2NativeSourceEndpointEventKind | None = None
            if endpoint_saturation_h <= 0.0:
                event_kind = K2NativeSourceEndpointEventKind.ATTACHED_HEXANE_DISAPPEARANCE
            elif endpoint_saturation_w <= 0.0:
                event_kind = K2NativeSourceEndpointEventKind.EXTERNAL_WATER_DISAPPEARANCE
            elif endpoint_saturation_sum >= 1.0:
                event_kind = K2NativeSourceEndpointEventKind.SHARED_FILM_AREA_CLOSURE
            if event_kind is not None:
                raise K2NativeSourceEndpointFirstEventRefusal(
                    f"group {basis.request.group_id} reaches or crosses endpoint FilmAreaLaw "
                    f"event {event_kind.value}; reject, do not clamp",
                    event_kind=event_kind,
                    group_id=basis.request.group_id,
                    accepted_time_s=request.tray_step.prior.tray.time_s,
                    endpoint_time_s=request.tray_step.end_time_s,
                    endpoint_attached_hexane_saturation=endpoint_saturation_h,
                    endpoint_external_water_saturation=endpoint_saturation_w,
                )
            for packet in basis.packets:
                entries[packet.packet_id] = tray_host.PacketPlanEntry(
                    packet_id=packet.packet_id,
                    payload_bytes=complete_after.payload_bytes,
                    inventory=complete_after.inventory,
                    weight=packet.entry.weight,
                    surface=complete_after.surface,
                )
            endpoints.append(
                NativeSourceGroupEndpoint(
                    group_id=basis.request.group_id,
                    vertical_layer_id=basis.request.vertical_layer_id,
                    packet_ids=basis.request.packet_ids,
                    hexane_access_oracle=basis.request.hexane_access_oracle,
                    dry_mass_fraction=basis.dry_mass_fraction,
                    representative_particles=basis.representative_particles,
                    accepted_t_n_attached_hexane_saturation=(basis.attached_hexane_saturation),
                    accepted_t_n_external_water_saturation=(basis.external_water_saturation),
                    endpoint_attached_hexane_saturation=endpoint_saturation_h,
                    endpoint_external_water_saturation=endpoint_saturation_w,
                    accepted_group_bed_subvolume_m3=(basis.accepted_group_bed_subvolume_m3),
                    hexane_inventory_capacity_kg=basis.hexane_inventory_capacity_kg,
                    water_inventory_capacity_kg=basis.water_inventory_capacity_kg,
                    active_area_m2=basis.active_area_m2,
                    interface_temperature_k=group.interface_temperature_k,
                    interface_hexane_mole_fraction=group.interface_hexane_mole_fraction,
                    total_molar_flux_mol_m2_s=group.total_molar_flux_mol_m2_s,
                    hexane_molar_flux_mol_m2_s=group.hexane_molar_flux_mol_m2_s,
                    water_molar_flux_mol_m2_s=group.water_molar_flux_mol_m2_s,
                    raw_impermeable_hexane_flux_mol_m2_s=(group.raw_hexane_molar_flux_mol_m2_s),
                    gas_to_interface_w=group.gas_to_interface_w,
                    wall_to_interface_w=group.wall_to_interface_w,
                    particle_to_interface_w=group.particle_to_interface_w,
                    native_hexane_partial_enthalpy_j_mol=(
                        group.native_hexane_partial_enthalpy_j_mol
                    ),
                    native_water_partial_enthalpy_j_mol=(group.native_water_partial_enthalpy_j_mol),
                    native_hexane_liquid_internal_energy_j_mol=(
                        group.native_hexane_liquid_internal_energy_j_mol
                    ),
                    native_water_liquid_internal_energy_j_mol=(
                        group.native_water_liquid_internal_energy_j_mol
                    ),
                    hexane_fugacity_residual_pa=group.hexane_fugacity_residual_pa,
                    water_fugacity_residual_pa=group.water_fugacity_residual_pa,
                    interface_energy_residual_w=group.interface_energy_residual_w,
                    interface_energy_resolution_bound_w=(group.interface_energy_resolution_bound_w),
                    endpoint_payload_identity=complete_after.canonical_identity,
                    endpoint_attached_hexane_kg=endpoint_attached_hexane,
                    endpoint_external_water_kg=endpoint_external_water,
                )
            )
    return entries, tuple(endpoints)


def _conservation_ledger(
    trial: _JointTrial,
    *,
    request: QSCK2NativeSourceRequest,
    host_transaction: tray_host.ValidatedKCellTrayStep,
) -> NativeSourceConservationLedger:
    step = request.tray_step
    dt = step.end_time_s - step.prior.tray.time_s
    before_by_id = {packet.packet_id: packet for packet in step.prior.tray.packets}
    after_by_id = {packet.packet_id: packet for packet in host_transaction.candidate.packets}
    h_flux: list[float] = []
    h_packets: list[float] = []
    w_flux: list[float] = []
    w_packets: list[float] = []
    e_flux: list[float] = []
    e_packets: list[float] = []
    cell_w: list[float] = []
    cell_h: list[float] = []
    cell_e: list[float] = []
    expected_e: list[float] = []
    residual_rows: list[tuple[float, tuple[float, ...]]] = []
    for layer_index, layer in enumerate(trial.layers, start=1):
        groups = layer.groups
        flux_h = dt * math.fsum(
            group.basis.active_area_m2 * group.hexane_molar_flux_mol_m2_s * hexane.M
            for group in groups
        )
        flux_w = dt * math.fsum(
            group.basis.active_area_m2 * group.water_molar_flux_mol_m2_s * water.M
            for group in groups
        )
        flux_e = dt * layer.particle_energy_loss_w
        ids = tuple(
            packet.packet_id
            for packet in step.prior.tray.packets
            if packet.owner_key.vertical_layer_id.value == layer_index
        )
        before_h = math.fsum(
            before_by_id[packet_id].weighted_totals.totals[2]
            + before_by_id[packet_id].weighted_totals.totals[3]
            for packet_id in ids
        )
        after_h = math.fsum(
            after_by_id[packet_id].weighted_totals.totals[2]
            + after_by_id[packet_id].weighted_totals.totals[3]
            for packet_id in ids
        )
        before_w = math.fsum(before_by_id[packet_id].weighted_totals.totals[4] for packet_id in ids)
        after_w = math.fsum(after_by_id[packet_id].weighted_totals.totals[4] for packet_id in ids)
        before_e = math.fsum(before_by_id[packet_id].weighted_totals.totals[6] for packet_id in ids)
        after_e = math.fsum(after_by_id[packet_id].weighted_totals.totals[6] for packet_id in ids)
        packet_h = math.fsum((before_h, -after_h))
        packet_w = math.fsum((before_w, -after_w))
        packet_e = math.fsum((before_e, -after_e))
        h_flux.append(flux_h)
        h_packets.append(packet_h)
        w_flux.append(flux_w)
        w_packets.append(packet_w)
        e_flux.append(flux_e)
        e_packets.append(packet_e)
        residual_rows.extend(
            (
                (packet_h - flux_h, (before_h, after_h, packet_h, flux_h)),
                (packet_w - flux_w, (before_w, after_w, packet_w, flux_w)),
                (packet_e - flux_e, (before_e, after_e, packet_e, flux_e)),
            )
        )
        before_cell = step.prior.tray.cells[layer_index - 1].conserved_state.totals
        after_cell = host_transaction.candidate.cells[layer_index - 1].conserved_state.totals
        water_storage = after_cell.water_kg - before_cell.water_kg
        hexane_storage = after_cell.hexane_kg - before_cell.hexane_kg
        energy_storage = after_cell.common_datum_energy_j - before_cell.common_datum_energy_j
        expected_energy = (
            step.layer_models[0].wall.thermal_capacity_j_k
            * (trial.shared_wall_temperature_k - step.prior.wall_temperature_k)
            if layer_index == 1
            else 0.0
        )
        cell_w.append(water_storage)
        cell_h.append(hexane_storage)
        cell_e.append(energy_storage - expected_energy)
        expected_e.append(expected_energy)
        residual_rows.extend(
            (
                (water_storage, (after_cell.water_kg, before_cell.water_kg)),
                (hexane_storage, (after_cell.hexane_kg, before_cell.hexane_kg)),
                (
                    energy_storage - expected_energy,
                    (energy_storage, expected_energy, after_cell.common_datum_energy_j),
                ),
            )
        )
    maximum = max(_relative(residual, *terms) for residual, terms in residual_rows)
    return NativeSourceConservationLedger(
        layer_hexane_flux_kg=(h_flux[0], h_flux[1]),
        layer_weighted_packet_hexane_decrement_kg=(h_packets[0], h_packets[1]),
        layer_water_flux_kg=(w_flux[0], w_flux[1]),
        layer_weighted_packet_water_decrement_kg=(w_packets[0], w_packets[1]),
        layer_particle_energy_flux_j=(e_flux[0], e_flux[1]),
        layer_weighted_packet_energy_decrement_j=(e_packets[0], e_packets[1]),
        cell_water_storage_residual_kg=(cell_w[0], cell_w[1]),
        cell_hexane_storage_residual_kg=(cell_h[0], cell_h[1]),
        cell_energy_storage_residual_j=(cell_e[0], cell_e[1]),
        expected_cell_energy_storage_j=(expected_e[0], expected_e[1]),
        maximum_relative_residual=maximum,
        relative_limit=step.relative_limit,
    )


def evaluate_qsc_k2_native_source_step(
    request: QSCK2NativeSourceRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner,
    integration_law_authority: single_wall.ThroughBedK2SingleWallIntegrationLawAuthority,
    expected_integration_law_definition_digest: str,
    expected_layer_model_configuration_digest: str,
) -> K2NativeWaterPreEventNoGo | ValidatedQSCK2NativeSourceStep | QSCK2NativeSourceRejection:
    """Solve one immutable K=2 native-source candidate and invoke the host once."""

    if type(request) is not QSCK2NativeSourceRequest:
        raise TypeError("native source evaluator requires an exact request")
    step = request.tray_step
    try:
        single_wall._validate_request_law_and_model_pins(
            step,
            integration_law_authority=integration_law_authority,
            expected_integration_law_definition_digest=(expected_integration_law_definition_digest),
            expected_layer_model_configuration_digest=(expected_layer_model_configuration_digest),
        )
        single_wall._validate_codec_and_k_planner(step, codec, k_planner)
        bases = _group_bases(request, codec)
        absent_water = tuple(_positive_zero(basis.external_water_saturation) for basis in bases)
        if any(absent_water) and not all(absent_water):
            raise QSCK2NativeSourceError(
                "one solve cannot mix absent- and present-external-water active sets"
            )
        seed, lower, upper = _seed_and_bounds(request, bases)
        seed_function_evaluations = 0
        if all(absent_water):
            condition_assessment = _water_condition_assessment(request, bases)
            seed, _phase_nfev = _phase_aware_seed(
                request,
                bases,
                seed,
                lower,
                upper,
                codec=codec,
            )
            _boundary_solution, boundary_trial, _boundary_nfev = _solve_joint_root(
                request,
                bases,
                seed,
                lower,
                upper,
                codec=codec,
                enforce_scaled_gate=False,
            )
            refusal = _water_appearance_boundary_refusal(
                request,
                boundary_trial,
                assessment_identity=condition_assessment.assessment_identity,
            )
            return K2NativeWaterPreEventNoGo(
                kind=(
                    K2NativeWaterPreEventNoGoKind.MISSING_ACCEPTED_NATIVE_HEXANE_ONLY_LAYER_STATE
                ),
                attempt_id=step.attempt_id,
                rollback_state=step.prior,
                accepted_state_digest=step.prior.state_digest,
                accepted_tray_state_digest=step.prior.tray.state_digest,
                accepted_time_s=step.prior.tray.time_s,
                access_neutral_assessment_identity=(condition_assessment.assessment_identity),
                counterfactual_phase_conditions=condition_assessment.conditions,
                counterfactual_boundary_refusal=refusal,
                missing_authority=(
                    "canonical GOV02 initialization has no accepted gas/wall/cell native "
                    "HEXANE_ONLY state; manufactured SP1 K1 separately refuses its zero-water "
                    "root outside the Law-2 Stefan-flux domain"
                ),
            )
        else:
            seed, phase_nfev = _phase_aware_seed(
                request,
                bases,
                seed,
                lower,
                upper,
                codec=codec,
            )
            seed_function_evaluations += phase_nfev

        _, trial, solve_function_evaluations = _solve_joint_root(
            request,
            bases,
            seed,
            lower,
            upper,
            codec=codec,
        )
        maximum = max(abs(value) for value in trial.scaled_residuals)
        aggregate_pressure = _aggregate_physical_tray_pressure_ledger(
            request,
            trial,
        )
        entries, endpoints = _endpoint_entries(trial, request=request, codec=codec)
        if not all(endpoint.passed for endpoint in endpoints):
            raise QSCK2NativeSourceError("a native group endpoint failed its phase/access gate")
        feeder = _EndpointFeeder(entries=entries, maximum_scaled_residual=step.newton_tolerance)
        host_request = tray_host.KCellTrayStepRequest(
            attempt_id=step.attempt_id,
            prior=step.prior.tray,
            end_time_s=step.end_time_s,
            flow_segments=step.flow_segments,
            external_cell_transfers=tuple(
                layer.result.external_cell_transfer for layer in trial.layers
            ),
            relative_limit=step.relative_limit,
        )
        host_result = tray_host.evaluate_k_cell_tray_step(
            host_request,
            packet_callback=feeder,
            k_planner=k_planner,
            auditor=codec,
        )
        if type(host_result) is tray_host.KCellTrayStepRejection:
            raise QSCK2NativeSourceError(
                f"weighted host rejected native endpoints: {host_result.reason}"
            )
        if type(host_result) is not tray_host.ValidatedKCellTrayStep:
            raise QSCK2NativeSourceError("weighted host returned a foreign transaction")
        ledger = _conservation_ledger(trial, request=request, host_transaction=host_result)
        if not ledger.passed:
            raise QSCK2NativeSourceError(
                "native before/after, cell-storage, or no-double-count ledger did not pass: "
                f"maximum={ledger.maximum_relative_residual:.17g}, "
                f"limit={ledger.relative_limit:.17g}"
            )
        endpoint_entries = tuple(entries[packet.packet_id] for packet in step.prior.tray.packets)
        return ValidatedQSCK2NativeSourceStep(
            request=request,
            host_transaction=host_result,
            layer_results=(trial.layers[0].result, trial.layers[1].result),
            group_endpoints=endpoints,
            shared_wall_temperature_k=trial.shared_wall_temperature_k,
            shared_wall_residual_w=trial.shared_wall_residual_w,
            shared_wall_resolution_bound_w=trial.shared_wall_resolution_bound_w,
            maximum_scaled_residual=maximum,
            scipy_function_evaluations=(seed_function_evaluations + solve_function_evaluations),
            aggregate_pressure=aggregate_pressure,
            conservation_ledger=ledger,
            endpoint_entries=endpoint_entries,
        )
    except (
        QSCK2NativeSourceError,
        single_wall.ThroughBedK2SingleWallIntegrationError,
        high_loading.HighLoadingPacketAdapterError,
        tray_host.KCellTrayHostError,
        closure.CellClosureError,
        engineering.EngineeringFeasibilityError,
        ValueError,
        OverflowError,
        ZeroDivisionError,
    ) as error:
        return QSCK2NativeSourceRejection(
            attempt_id=step.attempt_id,
            reason=f"native K2 source step failed: {type(error).__name__}: {error}",
            rollback_state=step.prior,
            endpoint_first_event=(
                error if isinstance(error, K2NativeSourceEndpointFirstEventRefusal) else None
            ),
        )


__all__ = (
    "K2NativeSourceEndpointEventKind",
    "K2NativeSourceEndpointFirstEventRefusal",
    "K2NativeSourceGroup",
    "K2NativeWaterAppearanceBoundaryRefusal",
    "K2NativeWaterAppearanceCondition",
    "K2NativeWaterAppearanceRefusalKind",
    "K2NativeWaterBoundaryFluxSign",
    "K2NativeWaterBoundaryLayerEvidence",
    "K2NativeWaterPreEventNoGo",
    "K2NativeWaterPreEventNoGoKind",
    "NativeFilmCoefficients",
    "NativeSourceConservationLedger",
    "NativeSourceGroupEndpoint",
    "NativeSourceLayerResult",
    "QSC_K2_NATIVE_ENDPOINT_FEEDER_ID",
    "QSC_K2_NATIVE_SOURCE_SCHEMA_ID",
    "QSC_K2_NATIVE_SOURCE_SOLVER_ID",
    "QSCK2NativeSourceError",
    "QSCK2NativeSourceRejection",
    "QSCK2NativeSourceRequest",
    "ValidatedQSCK2NativeSourceStep",
    "evaluate_qsc_k2_native_source_step",
)
