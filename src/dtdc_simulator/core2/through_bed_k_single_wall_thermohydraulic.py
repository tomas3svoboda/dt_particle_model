"""Arbitrary-K through-bed thermohydraulic successor with one physical wall.

This bounded successor leaves the accepted packet population, K owners, RTD
clocks, and host cell conserved states untouched.  It reads packet-owned layer
inventories through the existing tag-1 adapter, solves the quasi-steady Law-2
gas cells in bottom-to-top series, and returns only candidate fast states and
physical ledgers for a later packet/host transaction.

The resolution rules are physical rather than per-numerical-cell:

* every K cell carries its own full two-term Ergun bed segment;
* exactly K1 carries the one floor, indirect-duty record, and dynamic wall;
* the printed 1000--5600 Pa hull is applied once to the aggregate physical-tray
  bed-plus-one-floor drop; local cell pressure corroboration is disabled;
* adjacent numerical layers exchange simultaneous-t_n bed mixture energy through
  ``q = k_mixL A (T_lower - T_upper) / dz`` with the D-T2-1 locked
  effective axial ``k_mixL = 0.24 W m-1 K-1`` authority.  This is neither a
  particle-material conductivity nor an added gas dispersed-enthalpy flux.
  Face energy is exactly antisymmetric, and a step that would reverse any old
  adjacent temperature ordering is refused for substepping.

``ZERO_FLUX_PARITY`` is a diagnostic nested limit used only to compare this
successor with the committed K=2 compatibility oracle.  It does not replace
the locked active axial mixture-energy model.

K=4 resolves a physical SP1 top-cell root below the historical 102325 Pa
six-layer admission.  The distinct noncertifying arbitrary-K Law-2 runtime
admission reaches the 101325 Pa downstream boundary and now admits that root.
The frozen single-tray certificate remains unchanged; this module never clips,
widens, or seeds around the runtime domain.

This is an engineering prequalification seam, not QSC-10 completion, physical
qualification, calibration, plant prediction, or production wiring.
"""

from __future__ import annotations

import enum
import math
import struct
from dataclasses import dataclass, replace
from typing import ClassVar

from . import cell_closure as cc
from . import cell_engineering_feasibility as ef
from . import engineering_dual_variant_packet_adapter as dual_variant
from . import high_loading_packet_adapter as high_loading
from . import k_cell_tray_host as tray_host
from . import sp1_k2_law2_tray_integration as legacy
from . import through_bed_k2_law2_kernel as k2_kernel
from .dtdc_stack import REFERENCE_TRAYS
from .tray_type import reference_tray_type

SUPPORTED_REFERENCE_TRAYS = ("MN1", "MN2", "SP1")
AGGREGATE_PHYSICAL_TRAY_DROP_HULL_PA = cc.CORROBORATED_LAYER_DROP_HULL_PA
AGGREGATE_PHYSICAL_TRAY_PRESSURE_AUTHORITY_ID = (
    "frozen-phy017-physical-tray-printed-drop-hull-1000-5600-pa"
)
THROUGH_BED_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K = 0.24
THROUGH_BED_INTER_K_CONDUCTIVITY_AUTHORITY_ID = (
    "frozen-t2-k-mixl-effective-bed-conductivity-0p24-w-m-k"
)
THROUGH_BED_INTER_K_DECISION_ID = "D-T2-1-LOCKED-2026-08-06"
THROUGH_BED_INTER_K_ALLOCATION_LAW_ID = "static-layer-representative-caloric-guard-v1"
K4_EXECUTABLE_DEMONSTRATED = True
LAW2_LOW_PRESSURE_EXTENSION_REQUIRED_FOR_K4 = False


class ThroughBedKThermohydraulicError(ValueError):
    """The variable-K state, model, or thermohydraulic step is inadmissible."""


class ThroughBedKConfigurationError(ThroughBedKThermohydraulicError):
    """Static topology, geometry, or boundary declarations are inconsistent."""


class ThroughBedKStepError(ThroughBedKThermohydraulicError):
    """A reversible thermohydraulic trial needs rejection or substepping."""


class AdjacentMixtureEnergyMode(enum.Enum):
    ZERO_FLUX_PARITY = "zero_flux_k2_compatibility_parity"
    LOCKED_T2_EFFECTIVE_CONDUCTION = "locked_t2_effective_conduction"


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedInterKEnergyConductivityAuthority:
    decision_id: str = THROUGH_BED_INTER_K_DECISION_ID
    authority_id: str = THROUGH_BED_INTER_K_CONDUCTIVITY_AUTHORITY_ID
    effective_conductivity_w_m_k: float = THROUGH_BED_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K
    equation: str = "q_lower_to_upper=k_mixL*A*(T_lower-T_upper)/delta_z"
    allocation_law_id: str = THROUGH_BED_INTER_K_ALLOCATION_LAW_ID

    source_decision_locked: ClassVar[bool] = True
    lambda_mix_used: ClassVar[bool] = False
    particle_solid_conductivity_used: ClassVar[bool] = False
    gas_dispersed_enthalpy_term_added: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if (
            self.decision_id != THROUGH_BED_INTER_K_DECISION_ID
            or self.authority_id != THROUGH_BED_INTER_K_CONDUCTIVITY_AUTHORITY_ID
            or self.effective_conductivity_w_m_k != THROUGH_BED_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K
            or self.equation != "q_lower_to_upper=k_mixL*A*(T_lower-T_upper)/delta_z"
            or self.allocation_law_id != THROUGH_BED_INTER_K_ALLOCATION_LAW_ID
        ):
            raise ThroughBedKConfigurationError("inter-K conduction authority drifted")


THROUGH_BED_INTER_K_ENERGY_CONDUCTIVITY_AUTHORITY = ThroughBedInterKEnergyConductivityAuthority()


def _require_finite(name: str, value: float, *, positive: bool = False) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise ThroughBedKConfigurationError(f"{name} must be a finite exact binary64")
    if positive and value <= 0.0:
        raise ThroughBedKConfigurationError(f"{name} must be strictly positive")


def _same_binary64(left: float, right: float) -> bool:
    return struct.pack(">d", left) == struct.pack(">d", right)


def _relative_residual(residual: float, *terms: float) -> float:
    return abs(residual) / max(1.0, *(abs(value) for value in terms))


def _physical_topology(state: tray_host.AcceptedKCellTrayState):
    return next(
        item for item in state.topology.trays if item.physical_tray_id == state.physical_tray_id
    )


def inactive_k_wall_parameters(physical_wall: cc.WallNodeParameters) -> cc.WallNodeParameters:
    """Return the exact zero-storage/zero-boundary placeholder for K>1."""

    if type(physical_wall) is not cc.WallNodeParameters:
        raise TypeError("physical wall must be exact WallNodeParameters")
    return replace(
        physical_wall,
        steam_side_ua_w_k=0.0,
        ambient_ua_w_k=0.0,
        closure=cc.WallClosure.QUASI_STEADY_LIMIT,
    )


def _validate_layer_models(
    physical_tray_id: str,
    models: tuple[legacy.SP1KLayerEngineeringModel, ...],
) -> None:
    if physical_tray_id not in SUPPORTED_REFERENCE_TRAYS:
        raise ThroughBedKConfigurationError("variable-K through-bed tray must be MN1, MN2, or SP1")
    if (
        type(models) is not tuple
        or len(models) < 2
        or any(type(model) is not legacy.SP1KLayerEngineeringModel for model in models)
    ):
        raise ThroughBedKConfigurationError(
            "variable-K integration requires at least two exact layer models"
        )

    reference = next(item for item in REFERENCE_TRAYS if item.tray_id == physical_tray_id)
    tray_type = reference_tray_type(reference)
    physical_wall = models[0].wall
    inactive_wall = inactive_k_wall_parameters(physical_wall)
    floors: list[cc.TaperedBoreFloorPassageElement] = []
    properties = models[0].properties
    for layer, model in enumerate(models, start=1):
        tray = model.geometry.tray
        if tray.tray_id != f"{physical_tray_id}:K{layer}":
            raise ThroughBedKConfigurationError(
                "K models are not canonically identified and ordered"
            )
        if model.geometry.tray_type != tray_type:
            raise ThroughBedKConfigurationError("K layer changed the inherited tray type")
        if (
            tray.role != reference.role
            or tray.contact is not reference.contact
            or tray.diameter_m != reference.diameter_m
        ):
            raise ThroughBedKConfigurationError(
                "K layer changed its physical tray role or diameter"
            )
        bed = model.hydraulics.series.bed_element
        if bed is None or bed.flow_length_m != tray.loaded_depth_m or bed.resistance_scale != 1.0:
            raise ThroughBedKConfigurationError(
                "each K cell requires one active Ergun segment over its own depth"
            )
        floor_rows = tuple(
            item
            for item in model.hydraulics.series.elements
            if type(item) is cc.TaperedBoreFloorPassageElement
        )
        if len(floor_rows) != 1:
            raise ThroughBedKConfigurationError("each K series requires one typed floor coordinate")
        floors.append(floor_rows[0])
        if model.properties != properties:
            raise ThroughBedKConfigurationError("all K cells require one gas caloric convention")
        if layer == 1:
            if physical_wall.closure is not cc.WallClosure.DYNAMIC_NODE:
                raise ThroughBedKConfigurationError("K1 must own the sole dynamic physical wall")
            if physical_wall.steam_side_ua_w_k <= 0.0:
                raise ThroughBedKConfigurationError("K1 physical wall requires steam conductance")
            if (
                model.transfer.wall_to_gas_ua_w_k == 0.0
                and model.transfer.wall_to_interface_ua_w_k == 0.0
            ):
                raise ThroughBedKConfigurationError("K1 requires an active wall heat path")
        else:
            if model.wall != inactive_wall:
                raise ThroughBedKConfigurationError(
                    "every K>1 wall must be the exact inactive K1-derived placeholder"
                )
            if (
                model.transfer.wall_to_gas_ua_w_k,
                model.transfer.wall_to_interface_ua_w_k,
            ) != (0.0, 0.0):
                raise ThroughBedKConfigurationError("only K1 may own wall transfer conductance")

    floor_scales = tuple(floor.resistance_scale for floor in floors)
    if floor_scales != (1.0, *(0.0 for _ in models[1:])):
        raise ThroughBedKConfigurationError("exactly K1 must own the one physical floor")
    inactive_floor = replace(floors[0], resistance_scale=0.0)
    if any(floor != inactive_floor for floor in floors[1:]):
        raise ThroughBedKConfigurationError(
            "K>1 floor coordinates must be exact zero-resistance physical-floor placeholders"
        )
    duties = tuple(model.geometry.tray.indirect_duty_w for model in models)
    if duties != (reference.indirect_duty_w, *(0.0 for _ in models[1:])):
        raise ThroughBedKConfigurationError("exactly K1 must carry the one indirect-duty record")
    if math.fsum(model.geometry.tray.loaded_depth_m for model in models) != (
        reference.loaded_depth_m
    ):
        raise ThroughBedKConfigurationError("K depths do not partition the physical bed")
    reference_geometry = cc.CellGeometry(tray=reference, tray_type=tray_type)
    if math.fsum(model.geometry.gas_side_reference_area_m2 for model in models) != (
        reference_geometry.gas_side_reference_area_m2
    ):
        raise ThroughBedKConfigurationError("K gas-side areas do not partition the physical bed")
    cross_section = models[0].geometry.cross_section_m2
    if any(model.geometry.cross_section_m2 != cross_section for model in models[1:]):
        raise ThroughBedKConfigurationError("adjacent K layers must share one physical face area")


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedThroughBedKSingleWallState:
    tray: tray_host.AcceptedKCellTrayState
    cell_field_states: tuple[ef.BinaryNoInertCellFieldState, ...]
    wall_temperature_k: float

    packet_population_and_rtd_state_are_host_owned: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False
    executable_k4_demonstrated: ClassVar[bool] = True
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.tray) is not tray_host.AcceptedKCellTrayState:
            raise ThroughBedKConfigurationError("tray must be an exact accepted K-cell state")
        physical = _physical_topology(self.tray)
        # D9-e: the accepted fast states may now be the TEN-unknown dry-shell
        # state as well as the shipped binary one, so the gate tests exact
        # membership in the cell's own exported tuple rather than one name.
        # Exactly as before it is an EXACT-type test: an ``isinstance`` here
        # would admit any future subclass silently, and the dry-shell state is
        # not a subclass of the binary one at all, so only naming it admits it.
        # A genuinely foreign type still refuses typed, on this same line.
        if (
            type(self.cell_field_states) is not tuple
            or len(self.cell_field_states) != physical.vertical_layer_count
            or any(
                type(state) not in ef.ACCEPTED_FAST_BLOCK_STATE_TYPES
                for state in self.cell_field_states
            )
        ):
            raise ThroughBedKConfigurationError("accepted fast states must match arbitrary K")
        _require_finite(
            "accepted physical-wall temperature", self.wall_temperature_k, positive=True
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedKSingleWallStepRequest:
    prior: AcceptedThroughBedKSingleWallState
    end_time_s: float
    layer_models: tuple[legacy.SP1KLayerEngineeringModel, ...]
    bottom_gas_boundary: ef.EngineeringGasBoundary
    bottom_gas_binding: k2_kernel.ThroughBedGasBoundaryBinding
    top_boundary_pressure_pa: float
    initial_face_pressures_pa: tuple[float, ...]
    adjacent_mixture_energy_mode: AdjacentMixtureEnergyMode
    newton_tolerance: float = ef.DEFAULT_NEWTON_TOLERANCE
    face_tolerance: float = legacy.DEFAULT_FACE_TOLERANCE
    face_max_iterations: int = legacy.DEFAULT_FACE_MAX_ITERATIONS
    relative_limit: float = legacy.DEFAULT_RELATIVE_LIMIT

    physically_qualifying: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False
    executable_k4_demonstrated: ClassVar[bool] = True
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.prior) is not AcceptedThroughBedKSingleWallState:
            raise TypeError("variable-K request prior has a foreign exact type")
        AcceptedThroughBedKSingleWallState.__post_init__(self.prior)
        _require_finite("end time", self.end_time_s)
        if self.end_time_s <= self.prior.tray.time_s:
            raise ThroughBedKConfigurationError("thermohydraulic macro-step must advance time")
        _validate_layer_models(self.prior.tray.physical_tray_id, self.layer_models)
        k_count = len(self.layer_models)
        physical = _physical_topology(self.prior.tray)
        if physical.vertical_layer_count != k_count:
            raise ThroughBedKConfigurationError("model K differs from accepted host K")
        if len(self.prior.cell_field_states) != k_count:
            raise ThroughBedKConfigurationError("accepted seed K differs from model K")
        if type(self.bottom_gas_boundary) is not ef.EngineeringGasBoundary:
            raise TypeError("bottom gas boundary has a foreign exact type")
        if self.bottom_gas_boundary.carrier_topology is not ef.CarrierTopology.BINARY_NO_INERT:
            raise ThroughBedKConfigurationError("through-bed successor requires Law 2")
        if type(self.bottom_gas_binding) is not k2_kernel.ThroughBedGasBoundaryBinding:
            raise TypeError("bottom gas binding has a foreign exact type")
        if self.bottom_gas_binding.physical_tray_id != self.prior.tray.physical_tray_id:
            raise ThroughBedKConfigurationError("bottom gas binding belongs to another tray")
        _require_finite("top boundary pressure", self.top_boundary_pressure_pa, positive=True)
        if (
            type(self.initial_face_pressures_pa) is not tuple
            or len(self.initial_face_pressures_pa) != k_count - 1
        ):
            raise ThroughBedKConfigurationError("arbitrary K requires exactly K-1 pressure seeds")
        for pressure in self.initial_face_pressures_pa:
            _require_finite("initial inter-K face pressure", pressure, positive=True)
        if type(self.adjacent_mixture_energy_mode) is not AdjacentMixtureEnergyMode:
            raise TypeError("adjacent mixture-energy mode has a foreign exact type")
        _require_finite("Newton tolerance", self.newton_tolerance, positive=True)
        _require_finite("face tolerance", self.face_tolerance, positive=True)
        _require_finite("relative limit", self.relative_limit, positive=True)
        if type(self.face_max_iterations) is not int or self.face_max_iterations <= 0:
            raise ThroughBedKConfigurationError("face iteration limit must be a positive int")
        for model in self.layer_models:
            if model.properties.energy_datum_id != self.prior.tray.energy_datum_id:
                raise ThroughBedKConfigurationError(
                    "cell and accepted packets use different datums"
                )

    @property
    def duration_s(self) -> float:
        return self.end_time_s - self.prior.tray.time_s


@dataclass(frozen=True, slots=True, kw_only=True)
class InterKGasFaceLedger:
    lower_vertical_layer: int
    upper_vertical_layer: int
    lower_outlet_molar_flow_mol_s: float
    upper_inlet_molar_flow_mol_s: float
    lower_outlet_hexane_mole_fraction: float
    upper_inlet_hexane_mole_fraction: float
    lower_outlet_temperature_k: float
    upper_inlet_temperature_k: float
    lower_downstream_pressure_pa: float
    upper_layer_pressure_pa: float
    pressure_residual_pa: float
    pressure_relative_residual: float
    pressure_relative_limit: float
    iterations: int

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.upper_vertical_layer == self.lower_vertical_layer + 1
            and self.lower_outlet_molar_flow_mol_s == self.upper_inlet_molar_flow_mol_s
            and self.lower_outlet_hexane_mole_fraction == self.upper_inlet_hexane_mole_fraction
            and self.lower_outlet_temperature_k == self.upper_inlet_temperature_k
            and self.pressure_relative_residual <= self.pressure_relative_limit
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AggregatePhysicalTrayPressureLedger:
    bed_drop_by_layer_pa: tuple[float, ...]
    floor_drop_by_layer_pa: tuple[float, ...]
    active_physical_floor_count: int
    ordered_element_sum_pa: float
    cell_drop_sum_pa: float
    bottom_to_top_pressure_drop_pa: float
    cell_series_residuals_pa: tuple[float, ...]
    aggregate_pressure_residual_pa: float
    aggregate_gate_band_pa: tuple[float, float]
    inside_aggregate_gate: bool
    maximum_relative_residual: float
    relative_limit: float
    pressure_authority_id: str

    local_per_cell_pressure_gate_used: ClassVar[bool] = False
    aggregate_physical_tray_gate_used: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.active_physical_floor_count == 1
            and self.floor_drop_by_layer_pa[0] > 0.0
            and all(value == 0.0 for value in self.floor_drop_by_layer_pa[1:])
            and all(value > 0.0 for value in self.bed_drop_by_layer_pa)
            and self.inside_aggregate_gate
            and self.maximum_relative_residual <= self.relative_limit
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AdjacentMixtureEnergyFaceLedger:
    lower_vertical_layer: int
    upper_vertical_layer: int
    mode: AdjacentMixtureEnergyMode
    authority_id: str
    governed_effective_conductivity_w_m_k: float
    active_effective_conductivity_w_m_k: float
    face_area_m2: float
    gradient_distance_m: float
    conductance_w_k: float
    lower_temperature_before_k: float
    upper_temperature_before_k: float
    lower_temperature_after_k: float
    upper_temperature_after_k: float
    heat_lower_to_upper_j: float
    energy_to_lower_layer_j: float
    energy_to_upper_layer_j: float
    antisymmetry_residual_j: float
    entropy_production_j_k: float
    no_explicit_equilibrium_crossing: bool

    internal_energy_transfer_only: ClassVar[bool] = True
    particle_solid_conductivity_used: ClassVar[bool] = False
    gas_dispersed_enthalpy_term_added: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.upper_vertical_layer == self.lower_vertical_layer + 1
            and self.energy_to_lower_layer_j == -self.energy_to_upper_layer_j
            and self.antisymmetry_residual_j == 0.0
            and self.entropy_production_j_k >= 0.0
            and self.no_explicit_equilibrium_crossing
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AdjacentMixtureEnergyLedger:
    mode: AdjacentMixtureEnergyMode
    authority: ThroughBedInterKEnergyConductivityAuthority
    faces: tuple[AdjacentMixtureEnergyFaceLedger, ...]
    layer_temperatures_before_k: tuple[float, ...]
    layer_temperatures_after_k: tuple[float, ...]
    exact_antisymmetry_residual_j: float

    packet_payloads_advanced: ClassVar[bool] = False
    gas_energy_storage_advanced: ClassVar[bool] = False
    particle_solid_conductivity_used: ClassVar[bool] = False
    gas_dispersed_enthalpy_term_added: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            len(self.faces) == len(self.layer_temperatures_before_k) - 1
            and self.exact_antisymmetry_residual_j == 0.0
            and all(face.passed for face in self.faces)
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedKThermohydraulicStep:
    request: ThroughBedKSingleWallStepRequest
    layer_inventories: tuple[legacy.PacketDerivedLayerInventory, ...]
    cell_inputs: tuple[ef.EngineeringCellInputs, ...]
    cell_fast_solves: tuple[ef.BinaryNoInertCellSolve, ...]
    accepted_fast_steps: tuple[ef.AcceptedEngineeringFastStep, ...]
    inter_k_gas_faces: tuple[InterKGasFaceLedger, ...]
    aggregate_pressure: AggregatePhysicalTrayPressureLedger
    physical_wall: ef.EngineeringWallNodeSolve
    adjacent_mixture_energy: AdjacentMixtureEnergyLedger

    accepted_packet_population_advanced: ClassVar[bool] = False
    accepted_k_rtd_clocks_advanced: ClassVar[bool] = False
    host_commit_implemented: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False
    executable_k4_demonstrated: ClassVar[bool] = True
    law2_low_pressure_extension_required_for_k4: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def candidate_cell_field_states(self) -> tuple[ef.BinaryNoInertCellFieldState, ...]:
        return tuple(solve.state for solve in self.cell_fast_solves)

    @property
    def candidate_wall_temperature_k(self) -> float:
        return self.physical_wall.wall_temperature_k


def _layer_packets(
    tray: tray_host.AcceptedKCellTrayState,
    layer: int,
) -> tuple[tray_host.AcceptedTrayPacket, ...]:
    packets = tuple(
        packet for packet in tray.packets if packet.owner_key.vertical_layer_id.value == layer
    )
    if not packets:
        raise ThroughBedKStepError(f"K{layer} has no accepted packet population")
    return packets


def _packet_energy_at_temperature_j(
    packet: tray_host.AcceptedTrayPacket,
    temperature_k: float,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    complete = codec.decode(packet.entry.payload_bytes)
    inventory = complete.inventory
    dry = inventory.dry_matter_kg
    total_hexane = math.fsum((inventory.attached_hexane_kg, inventory.internal_hexane_kg))
    rebuilt = codec.initialize(
        total_hexane_loading=total_hexane / dry,
        temperature_k=temperature_k,
        external_free_water_loading=inventory.external_water_kg / dry,
        retained_water_loading=inventory.retained_water_kg / dry,
        oil_fraction=inventory.residual_oil_label_kg / dry,
    )
    return rebuilt.inventory.common_datum_energy_j * packet.entry.weight.representative_particles


def _layer_uniform_energy_j(
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
    temperature_k: float,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    return math.fsum(
        _packet_energy_at_temperature_j(packet, temperature_k, codec=codec) for packet in packets
    )


def _temperature_after_internal_heat(
    packets: tuple[tray_host.AcceptedTrayPacket, ...],
    *,
    temperature_before_k: float,
    internal_heat_j: float,
    global_minimum_temperature_k: float,
    global_maximum_temperature_k: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> float:
    if internal_heat_j == 0.0:
        return temperature_before_k
    energy_before = _layer_uniform_energy_j(
        packets,
        temperature_before_k,
        codec=codec,
    )
    target = energy_before + internal_heat_j
    if internal_heat_j > 0.0:
        lower = temperature_before_k
        upper = global_maximum_temperature_k
    else:
        lower = global_minimum_temperature_k
        upper = temperature_before_k
    if lower == upper:
        raise ThroughBedKStepError(
            "explicit inter-K conduction would leave the old temperature envelope; "
            "retry at a shorter endpoint"
        )
    energy_lower = _layer_uniform_energy_j(packets, lower, codec=codec)
    energy_upper = _layer_uniform_energy_j(packets, upper, codec=codec)
    if not energy_lower <= target <= energy_upper:
        raise ThroughBedKStepError(
            "explicit inter-K conduction would cross the old layer equilibrium; "
            "retry at a shorter endpoint"
        )
    for _ in range(120):
        middle = 0.5 * (lower + upper)
        if middle == lower or middle == upper:
            break
        if _layer_uniform_energy_j(packets, middle, codec=codec) < target:
            lower = middle
        else:
            upper = middle
    return 0.5 * (lower + upper)


def _adjacent_mixture_energy_ledger(
    request: ThroughBedKSingleWallStepRequest,
    inventories: tuple[legacy.PacketDerivedLayerInventory, ...],
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> AdjacentMixtureEnergyLedger:
    """Book the ruled axial bed mixture-energy term without changing packets.

    The quasi-steady gas block owns no transient gas-energy inventory.  The
    candidate layer temperatures therefore use the packet common-datum energy
    map solely as the available transient layer-energy coordinate.  That
    allocation does not reinterpret ``k_mixL`` as particle-solid conductivity
    and does not add a dispersed-enthalpy term to the serial gas boundary.
    """

    authority = THROUGH_BED_INTER_K_ENERGY_CONDUCTIVITY_AUTHORITY
    temperatures = tuple(item.solid_temperature_k for item in inventories)
    active_conductivity = (
        0.0
        if request.adjacent_mixture_energy_mode is AdjacentMixtureEnergyMode.ZERO_FLUX_PARITY
        else authority.effective_conductivity_w_m_k
    )
    duration = request.duration_s
    raw: list[tuple[float, float, float, float, float]] = []
    layer_internal_heat = [0.0 for _ in temperatures]
    for index, (lower_model, upper_model) in enumerate(
        zip(request.layer_models[:-1], request.layer_models[1:], strict=True)
    ):
        area = lower_model.geometry.cross_section_m2
        distance = 0.5 * math.fsum(
            (
                lower_model.geometry.tray.loaded_depth_m,
                upper_model.geometry.tray.loaded_depth_m,
            )
        )
        conductance = active_conductivity * area / distance
        delta = temperatures[index] - temperatures[index + 1]
        heat = conductance * delta * duration
        if heat == 0.0:
            heat = 0.0
        layer_internal_heat[index] = math.fsum((layer_internal_heat[index], -heat))
        layer_internal_heat[index + 1] = math.fsum((layer_internal_heat[index + 1], heat))
        entropy = (
            conductance * duration * delta * delta / (temperatures[index] * temperatures[index + 1])
        )
        raw.append((area, distance, conductance, heat, entropy))

    if request.adjacent_mixture_energy_mode is AdjacentMixtureEnergyMode.ZERO_FLUX_PARITY:
        temperatures_after = temperatures
    else:
        minimum = min(temperatures)
        maximum = max(temperatures)
        try:
            temperatures_after = tuple(
                _temperature_after_internal_heat(
                    _layer_packets(request.prior.tray, layer),
                    temperature_before_k=temperatures[layer - 1],
                    internal_heat_j=layer_internal_heat[layer - 1],
                    global_minimum_temperature_k=minimum,
                    global_maximum_temperature_k=maximum,
                    codec=codec,
                )
                for layer in range(1, len(temperatures) + 1)
            )
        except high_loading.HighLoadingPacketAdapterError as error:
            raise ThroughBedKStepError(
                "inter-K caloric guard left the accepted tag-1 branch; retry at a shorter endpoint"
            ) from error

    faces: list[AdjacentMixtureEnergyFaceLedger] = []
    for index, (area, distance, conductance, heat, entropy) in enumerate(raw):
        before_difference = temperatures[index] - temperatures[index + 1]
        after_difference = temperatures_after[index] - temperatures_after[index + 1]
        crossed = (before_difference > 0.0 > after_difference) or (
            before_difference < 0.0 < after_difference
        )
        if crossed:
            raise ThroughBedKStepError(
                "explicit inter-K conduction reversed an old adjacent temperature ordering; "
                "retry at a shorter endpoint"
            )
        to_upper = heat
        to_lower = 0.0 if heat == 0.0 else -heat
        residual = math.fsum((to_lower, to_upper))
        face = AdjacentMixtureEnergyFaceLedger(
            lower_vertical_layer=index + 1,
            upper_vertical_layer=index + 2,
            mode=request.adjacent_mixture_energy_mode,
            authority_id=authority.authority_id,
            governed_effective_conductivity_w_m_k=authority.effective_conductivity_w_m_k,
            active_effective_conductivity_w_m_k=active_conductivity,
            face_area_m2=area,
            gradient_distance_m=distance,
            conductance_w_k=conductance,
            lower_temperature_before_k=temperatures[index],
            upper_temperature_before_k=temperatures[index + 1],
            lower_temperature_after_k=temperatures_after[index],
            upper_temperature_after_k=temperatures_after[index + 1],
            heat_lower_to_upper_j=heat,
            energy_to_lower_layer_j=to_lower,
            energy_to_upper_layer_j=to_upper,
            antisymmetry_residual_j=residual,
            entropy_production_j_k=entropy,
            no_explicit_equilibrium_crossing=not crossed,
        )
        if not face.passed:
            raise ThroughBedKStepError("adjacent mixture-energy face ledger did not pass")
        faces.append(face)
    exact_residual = math.fsum(
        value
        for face in faces
        for value in (face.energy_to_lower_layer_j, face.energy_to_upper_layer_j)
    )
    ledger = AdjacentMixtureEnergyLedger(
        mode=request.adjacent_mixture_energy_mode,
        authority=authority,
        faces=tuple(faces),
        layer_temperatures_before_k=temperatures,
        layer_temperatures_after_k=temperatures_after,
        exact_antisymmetry_residual_j=exact_residual,
    )
    if not ledger.passed:
        raise ThroughBedKStepError("whole-tray adjacent mixture-energy ledger did not pass")
    return ledger


def _solve_serial_k_cells(
    request: ThroughBedKSingleWallStepRequest,
    inventories: tuple[legacy.PacketDerivedLayerInventory, ...],
) -> tuple[
    tuple[ef.EngineeringCellInputs, ...],
    tuple[ef.BinaryNoInertCellSolve, ...],
    tuple[ef.AcceptedEngineeringFastStep, ...],
    tuple[InterKGasFaceLedger, ...],
]:
    k_count = len(request.layer_models)
    face_pressures = request.initial_face_pressures_pa
    seeds = list(request.prior.cell_field_states)
    for iteration in range(1, request.face_max_iterations + 1):
        inputs: list[ef.EngineeringCellInputs] = []
        solves: list[ef.BinaryNoInertCellSolve] = []
        gas_boundaries: list[ef.EngineeringGasBoundary] = []
        for index in range(k_count):
            downstream_pressure = (
                face_pressures[index] if index < k_count - 1 else request.top_boundary_pressure_pa
            )
            if index == 0:
                boundary = replace(
                    request.bottom_gas_boundary,
                    downstream_boundary_pressure_pa=downstream_pressure,
                )
            else:
                lower = solves[index - 1]
                boundary = ef.EngineeringGasBoundary(
                    inlet_molar_flow_mol_s=lower.evaluation.outlet_molar_flow_mol_s,
                    inlet_hexane_mole_fraction=lower.state.gas_hexane_mole_fraction,
                    inlet_water_mole_fraction=lower.state.gas_water_mole_fraction,
                    inlet_temperature_k=lower.state.gas_temperature_k,
                    downstream_boundary_pressure_pa=downstream_pressure,
                    carrier_topology=ef.CarrierTopology.BINARY_NO_INERT,
                )
            cell_inputs = replace(
                legacy._cell_inputs(
                    model=request.layer_models[index],
                    inventory=inventories[index],
                    gas_boundary=boundary,
                    wall_temperature_k=request.prior.wall_temperature_k,
                    macro_step_s=request.duration_s,
                ),
                require_corroborated_pressure_drop=False,
            )
            solve = ef.solve_binary_no_inert_fast_block(
                cell_inputs,
                seed=seeds[index],
                tolerance=request.newton_tolerance,
            )
            seeds[index] = solve.state
            inputs.append(cell_inputs)
            solves.append(solve)
            gas_boundaries.append(boundary)

        computed_faces = tuple(
            solves[index + 1].state.layer_pressure_pa for index in range(k_count - 1)
        )
        relative_residuals = tuple(
            _relative_residual(computed - guessed, computed, guessed)
            for computed, guessed in zip(computed_faces, face_pressures, strict=True)
        )
        if max(relative_residuals) <= request.face_tolerance:
            accepted = tuple(
                ef.accept_engineering_fast_step(cell_inputs, solve)
                for cell_inputs, solve in zip(inputs, solves, strict=True)
            )
            if any(
                not _same_binary64(step.old_wall_temperature_k, request.prior.wall_temperature_k)
                for step in accepted
            ):
                raise ThroughBedKStepError("a K solve did not retain the one accepted old wall")
            faces = tuple(
                InterKGasFaceLedger(
                    lower_vertical_layer=index + 1,
                    upper_vertical_layer=index + 2,
                    lower_outlet_molar_flow_mol_s=solves[index].evaluation.outlet_molar_flow_mol_s,
                    upper_inlet_molar_flow_mol_s=gas_boundaries[index + 1].inlet_molar_flow_mol_s,
                    lower_outlet_hexane_mole_fraction=solves[index].state.gas_hexane_mole_fraction,
                    upper_inlet_hexane_mole_fraction=(
                        gas_boundaries[index + 1].inlet_hexane_mole_fraction
                    ),
                    lower_outlet_temperature_k=solves[index].state.gas_temperature_k,
                    upper_inlet_temperature_k=gas_boundaries[index + 1].inlet_temperature_k,
                    lower_downstream_pressure_pa=face_pressures[index],
                    upper_layer_pressure_pa=computed_faces[index],
                    pressure_residual_pa=computed_faces[index] - face_pressures[index],
                    pressure_relative_residual=relative_residuals[index],
                    pressure_relative_limit=request.face_tolerance,
                    iterations=iteration,
                )
                for index in range(k_count - 1)
            )
            if not all(face.passed for face in faces):
                raise ThroughBedKStepError("an inter-K gas face did not close")
            return tuple(inputs), tuple(solves), accepted, faces
        face_pressures = computed_faces
    raise ThroughBedKStepError(
        "the arbitrary-K internal pressure faces did not converge without clipping"
    )


def _aggregate_pressure_ledger(
    request: ThroughBedKSingleWallStepRequest,
    inputs: tuple[ef.EngineeringCellInputs, ...],
    solves: tuple[ef.BinaryNoInertCellSolve, ...],
) -> AggregatePhysicalTrayPressureLedger:
    drops = tuple(
        cc.layer_pressure_drop(
            series=cell_inputs.hydraulics.series,
            superficial_velocity_m_s=solve.evaluation.superficial_velocity_m_s,
            gas_density_kg_m3=solve.evaluation.gas_density_kg_m3,
        )
        for cell_inputs, solve in zip(inputs, solves, strict=True)
    )
    bed = tuple(drop.contribution_for(cc.BedErgunElement.label).drop_pa for drop in drops)
    floor = tuple(
        drop.contribution_for(cc.TaperedBoreFloorPassageElement.label).drop_pa for drop in drops
    )
    ordered = math.fsum(
        value
        for bed_value, floor_value in zip(bed, floor, strict=True)
        for value in (bed_value, floor_value)
    )
    cell_sum = math.fsum(solve.evaluation.layer_pressure_drop_pa for solve in solves)
    endpoint = solves[0].state.layer_pressure_pa - request.top_boundary_pressure_pa
    cell_residuals = tuple(
        drop.drop_pa - solve.evaluation.layer_pressure_drop_pa
        for drop, solve in zip(drops, solves, strict=True)
    )
    aggregate_residual = endpoint - ordered
    # W3-1 (owner ruling 2026-09-03, executed at HEAD f5a70bb): the aggregate
    # residual telescopes K Newton pressure rows plus the K-1 internal face
    # Picard slips, so it is gated on that derived absolute resolution instead
    # of on a fixed relative literal divided by the tray drop.  Same intent
    # (pressure closure across the tray), same `relative_limit`.
    resolution_bound = legacy.aggregate_pressure_binary64_resolution_bound(
        newton_tolerance=request.newton_tolerance,
        pressure_row_count=len(solves),
        face_tolerance=request.face_tolerance,
        face_pressures_pa=tuple(solve.state.layer_pressure_pa for solve in solves[1:]),
        summed_terms_pa=(*bed, *floor, ordered, cell_sum, endpoint),
    )
    relative_rows = tuple(
        _relative_residual(
            legacy._residual_beyond_resolution(residual, resolution_bound),
            drop.drop_pa,
            solve.evaluation.layer_pressure_drop_pa,
        )
        for residual, drop, solve in zip(cell_residuals, drops, solves, strict=True)
    )
    maximum = max(
        *relative_rows,
        _relative_residual(
            legacy._residual_beyond_resolution(ordered - cell_sum, resolution_bound),
            ordered,
            cell_sum,
        ),
        _relative_residual(
            legacy._residual_beyond_resolution(aggregate_residual, resolution_bound),
            endpoint,
            ordered,
        ),
    )
    low, high = AGGREGATE_PHYSICAL_TRAY_DROP_HULL_PA
    ledger = AggregatePhysicalTrayPressureLedger(
        bed_drop_by_layer_pa=bed,
        floor_drop_by_layer_pa=floor,
        active_physical_floor_count=sum(value > 0.0 for value in floor),
        ordered_element_sum_pa=ordered,
        cell_drop_sum_pa=cell_sum,
        bottom_to_top_pressure_drop_pa=endpoint,
        cell_series_residuals_pa=cell_residuals,
        aggregate_pressure_residual_pa=aggregate_residual,
        aggregate_gate_band_pa=AGGREGATE_PHYSICAL_TRAY_DROP_HULL_PA,
        inside_aggregate_gate=low <= ordered <= high,
        maximum_relative_residual=maximum,
        relative_limit=request.relative_limit,
        pressure_authority_id=AGGREGATE_PHYSICAL_TRAY_PRESSURE_AUTHORITY_ID,
    )
    if not ledger.passed:
        raise ThroughBedKStepError(
            "aggregate physical-tray bed-plus-one-floor pressure gate did not pass"
        )
    return ledger


def _advance_one_physical_wall(
    request: ThroughBedKSingleWallStepRequest,
    accepted: tuple[ef.AcceptedEngineeringFastStep, ...],
) -> ef.EngineeringWallNodeSolve:
    if len(accepted) != len(request.layer_models):
        raise ThroughBedKStepError("physical wall and K fast-step counts differ")
    if tuple(step.wall_storage_scale for step in accepted) != (
        1.0,
        *(0.0 for _ in accepted[1:]),
    ):
        raise ThroughBedKStepError("exactly K1 may own the physical wall capacity")
    if any((step.wall_to_gas_w, step.wall_to_interface_w) != (0.0, 0.0) for step in accepted[1:]):
        raise ThroughBedKStepError("a K>1 placeholder acquired physical wall heat")
    wall = request.layer_models[0].wall
    dt = request.duration_s
    storage_conductance = wall.thermal_capacity_j_k / dt
    denominator = storage_conductance + wall.steam_side_ua_w_k + wall.ambient_ua_w_k
    to_gas = math.fsum(tuple(step.wall_to_gas_w for step in accepted))
    to_interface = math.fsum(tuple(step.wall_to_interface_w for step in accepted))
    numerator = math.fsum(
        (
            storage_conductance * request.prior.wall_temperature_k,
            wall.steam_side_ua_w_k * wall.steam_temperature_k,
            wall.ambient_ua_w_k * wall.ambient_temperature_k,
            -to_gas,
            -to_interface,
        )
    )
    temperature = numerator / denominator
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ThroughBedKStepError("physical wall update produced an invalid temperature")
    for layer, step in enumerate(accepted, start=1):
        for receiver, conductance, receiver_temperature in (
            ("gas", step.transfer.wall_to_gas_ua_w_k, step.state.gas_temperature_k),
            (
                "interface",
                step.transfer.wall_to_interface_ua_w_k,
                step.state.interface_temperature_k,
            ),
        ):
            if conductance == 0.0:
                continue
            old_delta = request.prior.wall_temperature_k - receiver_temperature
            new_delta = temperature - receiver_temperature
            if (old_delta < 0.0 < new_delta) or (new_delta < 0.0 < old_delta):
                raise ThroughBedKStepError(
                    f"physical wall would reverse K{layer} wall-to-{receiver} ordering; "
                    "retry at a shorter endpoint"
                )
    steam = wall.steam_side_ua_w_k * (wall.steam_temperature_k - temperature)
    ambient = wall.ambient_ua_w_k * (temperature - wall.ambient_temperature_k)
    storage = storage_conductance * (temperature - request.prior.wall_temperature_k)
    residual = math.fsum((storage, -steam, ambient, to_gas, to_interface))
    return ef.EngineeringWallNodeSolve(
        wall_temperature_k=temperature,
        old_wall_temperature_k=request.prior.wall_temperature_k,
        steam_side_w=steam,
        to_ambient_w=ambient,
        shared_to_gas_w=to_gas,
        shared_to_interface_w=to_interface,
        storage_w=storage,
        residual_w=residual,
        row_derivative_w_k=denominator,
    )


def evaluate_through_bed_k_single_wall_thermohydraulic_step(
    request: ThroughBedKSingleWallStepRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> ThroughBedKThermohydraulicStep:
    """Evaluate a reversible arbitrary-K gas/wall/axial-energy candidate.

    No accepted packet, K owner, RTD clock, or host cell is advanced here.
    K=4 remains subject to the native Law-2 state-domain guard; this function
    deliberately provides no pressure clipping or resolution-specific bypass.
    """

    if type(request) is not ThroughBedKSingleWallStepRequest:
        raise TypeError("variable-K evaluation requires an exact request")
    ThroughBedKSingleWallStepRequest.__post_init__(request)
    if type(codec) is not high_loading.HighLoadingV1PacketAdapter and (
        type(codec) is not dual_variant.EngineeringDualVariantPacketAdapter
    ):
        raise TypeError(
            "variable-K evaluation requires the exact tag-1 adapter "
            "or the owner-ruled dual-variant adapter"
        )
    codec.require_host_binding(
        configuration_digest=request.prior.tray.configuration_digest,
        energy_datum_id=request.prior.tray.energy_datum_id,
        pressure_pa=request.prior.tray.pressure_pa,
    )
    if (
        codec.codec_registry_digest != request.prior.tray.codec_registry_digest
        or codec.auditor_identity_digest != request.prior.tray.auditor_identity_digest
    ):
        raise ThroughBedKConfigurationError("codec differs from accepted packet host")
    inventories = tuple(
        legacy._layer_inventory(
            request.prior.tray,
            layer,
            model=request.layer_models[layer - 1],
            codec=codec,
        )
        for layer in range(1, len(request.layer_models) + 1)
    )
    adjacent_energy = _adjacent_mixture_energy_ledger(request, inventories, codec=codec)
    inputs, solves, accepted, gas_faces = _solve_serial_k_cells(request, inventories)
    pressure = _aggregate_pressure_ledger(request, inputs, solves)
    wall = _advance_one_physical_wall(request, accepted)
    return ThroughBedKThermohydraulicStep(
        request=request,
        layer_inventories=inventories,
        cell_inputs=inputs,
        cell_fast_solves=solves,
        accepted_fast_steps=accepted,
        inter_k_gas_faces=gas_faces,
        aggregate_pressure=pressure,
        physical_wall=wall,
        adjacent_mixture_energy=adjacent_energy,
    )


__all__ = (
    "AGGREGATE_PHYSICAL_TRAY_DROP_HULL_PA",
    "AGGREGATE_PHYSICAL_TRAY_PRESSURE_AUTHORITY_ID",
    "AcceptedThroughBedKSingleWallState",
    "AdjacentMixtureEnergyFaceLedger",
    "AdjacentMixtureEnergyLedger",
    "AdjacentMixtureEnergyMode",
    "AggregatePhysicalTrayPressureLedger",
    "InterKGasFaceLedger",
    "K4_EXECUTABLE_DEMONSTRATED",
    "LAW2_LOW_PRESSURE_EXTENSION_REQUIRED_FOR_K4",
    "SUPPORTED_REFERENCE_TRAYS",
    "THROUGH_BED_INTER_K_ALLOCATION_LAW_ID",
    "THROUGH_BED_INTER_K_CONDUCTIVITY_AUTHORITY_ID",
    "THROUGH_BED_INTER_K_DECISION_ID",
    "THROUGH_BED_INTER_K_ENERGY_CONDUCTIVITY_AUTHORITY",
    "THROUGH_BED_INTER_K_EFFECTIVE_CONDUCTIVITY_W_M_K",
    "ThroughBedInterKEnergyConductivityAuthority",
    "ThroughBedKConfigurationError",
    "ThroughBedKSingleWallStepRequest",
    "ThroughBedKStepError",
    "ThroughBedKThermohydraulicError",
    "ThroughBedKThermohydraulicStep",
    "evaluate_through_bed_k_single_wall_thermohydraulic_step",
    "inactive_k_wall_parameters",
)
