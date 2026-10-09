"""Shared nonqualifying K=2 / RTD=8 Law-2 kernel for reference through-bed trays.

This module is deliberately an internal normalization seam.  It lets the
legacy SP1 public transaction and the generic MN1/MN2/SP1 public transaction
execute one implementation of the cell, wall, packet, pressure, and host
algorithms.  The legacy public types and their digest domain remain in
``sp1_k2_law2_tray_integration``.

The reference ``indirect_duty_w`` values are geometry/configuration records
only.  Actual wall heat continues to come from the declared steam-side wall
closure; this kernel neither prescribes nor claims realization of 140/280 kW.
All qualification, calibration, production, and plant-prediction claims remain
false.
"""

from __future__ import annotations

import enum
import hashlib
import math
from dataclasses import dataclass, replace
from typing import ClassVar

from . import cell_closure as cc
from . import cell_engineering_feasibility as ef
from . import component_energy_datum_adapter as component_datum
from . import engineering_physical_k_event_scheduler as k_scheduler
from . import engineering_dual_variant_packet_adapter as dual_variant
from . import high_loading_packet_adapter as high_loading
from . import k_cell_tray_host as tray_host
from . import sp1_k2_law2_tray_integration as legacy
from .dtdc_stack import REFERENCE_TRAYS, VaporContact
from .tray_type import (
    GasSource,
    PressureLossRelation,
    TrayFloorFamily,
    reference_tray_type,
)

THROUGH_BED_K2_KERNEL_SCHEMA_ID = "dtdc-core2-through-bed-k2-rtd8-law2-kernel-v1"
THROUGH_BED_LAYER_COUNT = 2
THROUGH_BED_RTD_STAGE_COUNT = 8
SUPPORTED_REFERENCE_TRAYS = ("MN1", "MN2", "SP1")


class ThroughBedK2KernelError(ValueError):
    """The normalized through-bed request or calculation is inadmissible."""


class ThroughBedGasBoundaryKind(enum.Enum):
    """Declared origin of one physical tray's lower gas boundary."""

    DIRECT_STEAM_SOURCE = "direct_steam_source"
    LOWER_TRAY_OUTLET = "lower_tray_outlet"


_EXPECTED_BOUNDARIES = {
    "SP1": (ThroughBedGasBoundaryKind.DIRECT_STEAM_SOURCE, None),
    "MN2": (ThroughBedGasBoundaryKind.LOWER_TRAY_OUTLET, "SP1"),
    "MN1": (ThroughBedGasBoundaryKind.LOWER_TRAY_OUTLET, "MN2"),
}


def _require_sha256(name: str, value: str) -> None:
    if (
        type(value) is not str
        or len(value) != 71
        or not value.startswith("sha256:")
        or any(character not in "0123456789abcdef" for character in value[7:])
    ):
        raise ThroughBedK2KernelError(f"{name} must be a canonical sha256 identity")


def _reference_tray(physical_tray_id: str):
    try:
        return next(item for item in REFERENCE_TRAYS if item.tray_id == physical_tray_id)
    except StopIteration as error:  # pragma: no cover - guarded by the supported set
        raise ThroughBedK2KernelError(
            f"{physical_tray_id!r} is not a supported reference through-bed tray"
        ) from error


@dataclass(frozen=True, slots=True, kw_only=True)
class ReferenceThroughBedK2Contract:
    """Pinned reference identity and lower-boundary role for one K2 tray."""

    physical_tray_id: str
    lower_gas_boundary_kind: ThroughBedGasBoundaryKind
    expected_lower_gas_source_tray_id: str | None

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    reference_indirect_duty_is_record_only: ClassVar[bool] = True
    reference_indirect_duty_realized: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self.physical_tray_id not in SUPPORTED_REFERENCE_TRAYS:
            raise ThroughBedK2KernelError(
                "the K2 kernel admits only the MN1, MN2, and SP1 reference trays"
            )
        if type(self.lower_gas_boundary_kind) is not ThroughBedGasBoundaryKind:
            raise ThroughBedK2KernelError("lower gas-boundary kind has a foreign type")
        expected = _EXPECTED_BOUNDARIES[self.physical_tray_id]
        if (self.lower_gas_boundary_kind, self.expected_lower_gas_source_tray_id) != expected:
            raise ThroughBedK2KernelError(
                "reference tray and lower gas-boundary provenance disagree"
            )
        reference = _reference_tray(self.physical_tray_id)
        tray_type = reference_tray_type(reference)
        if (
            reference.contact is not VaporContact.THROUGH_BED
            or tray_type.family is not TrayFloorFamily.STRIPPING_DECK
            or tray_type.pressure_loss_relation is not PressureLossRelation.PACKED_BED_ERGUN
        ):
            raise ThroughBedK2KernelError(
                "reference contract is not a family-2 through-bed Ergun tray"
            )
        expected_source = (
            GasSource.DIRECT_STEAM
            if self.lower_gas_boundary_kind is ThroughBedGasBoundaryKind.DIRECT_STEAM_SOURCE
            else GasSource.NONE
        )
        if tray_type.gas_source is not expected_source:
            raise ThroughBedK2KernelError("reference tray gas-source type disagrees with routing")

    @property
    def reference_indirect_duty_w(self) -> float:
        return _reference_tray(self.physical_tray_id).indirect_duty_w

    @property
    def definition_digest(self) -> str:
        reference = _reference_tray(self.physical_tray_id)
        tray_type = reference_tray_type(reference)
        digest = hashlib.sha256()
        for part in (
            THROUGH_BED_K2_KERNEL_SCHEMA_ID.encode("ascii"),
            self.physical_tray_id.encode("ascii"),
            self.lower_gas_boundary_kind.value.encode("ascii"),
            (self.expected_lower_gas_source_tray_id or "<external>").encode("ascii"),
            repr(reference).encode("utf-8"),
            repr(tray_type).encode("utf-8"),
            b"reference-indirect-duty-record-only",
        ):
            digest.update(len(part).to_bytes(8, "big"))
            digest.update(part)
        return "sha256:" + digest.hexdigest()


def reference_through_bed_k2_contract(
    physical_tray_id: str,
) -> ReferenceThroughBedK2Contract:
    """Return the exact reference contract for MN1, MN2, or SP1."""

    if type(physical_tray_id) is not str or physical_tray_id not in _EXPECTED_BOUNDARIES:
        raise ThroughBedK2KernelError(
            "reference through-bed contract requires exact tray ID MN1, MN2, or SP1"
        )
    kind, source = _EXPECTED_BOUNDARIES[physical_tray_id]
    return ReferenceThroughBedK2Contract(
        physical_tray_id=physical_tray_id,
        lower_gas_boundary_kind=kind,
        expected_lower_gas_source_tray_id=source,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedGasBoundaryBinding:
    """Provenance declaration for one supplied lower gas boundary.

    This record is not an inter-tray equality proof.  For MN trays the source
    state digest is mandatory so a later tower orchestrator can bind and replay
    the exact lower-tray outlet transaction.
    """

    physical_tray_id: str
    boundary_kind: ThroughBedGasBoundaryKind
    source_physical_tray_id: str | None
    authority_id: str
    source_state_digest: str | None = None

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    intertray_equality_proved: ClassVar[bool] = False
    authority_authenticated: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.authority_id) is not str or not self.authority_id.strip():
            raise ThroughBedK2KernelError("gas-boundary authority ID must be nonblank")
        contract = reference_through_bed_k2_contract(self.physical_tray_id)
        if self.boundary_kind is not contract.lower_gas_boundary_kind:
            raise ThroughBedK2KernelError("gas-boundary kind disagrees with the reference tray")
        if self.source_physical_tray_id != contract.expected_lower_gas_source_tray_id:
            raise ThroughBedK2KernelError(
                "gas-boundary source tray disagrees with the reference countercurrent order"
            )
        if self.boundary_kind is ThroughBedGasBoundaryKind.LOWER_TRAY_OUTLET:
            if self.source_state_digest is None:
                raise ThroughBedK2KernelError(
                    "an MN lower-tray gas boundary requires its source state digest"
                )
            _require_sha256("gas-boundary source state", self.source_state_digest)
        elif self.source_state_digest is not None:
            raise ThroughBedK2KernelError(
                "the external SP1 steam-source binding must not name a lower-tray state"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class KernelAcceptedThroughBedK2State:
    contract: ReferenceThroughBedK2Contract
    tray: tray_host.AcceptedKCellTrayState
    layer_model_configuration_digest: str
    cell_field_states: tuple[ef.BinaryNoInertCellFieldState, ef.BinaryNoInertCellFieldState]
    wall_temperatures_k: tuple[float, float]


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedK2KernelRequest:
    attempt_id: str
    prior: KernelAcceptedThroughBedK2State
    end_time_s: float
    layer_models: tuple[legacy.SP1KLayerEngineeringModel, legacy.SP1KLayerEngineeringModel]
    layer_model_configuration_digest: str
    bottom_gas_boundary: ef.EngineeringGasBoundary
    bottom_gas_binding: ThroughBedGasBoundaryBinding
    top_boundary_pressure_pa: float
    flow_segments: tuple[tray_host.AcceptedDryMatterFlowSegment, ...]
    newton_tolerance: float
    face_tolerance: float
    face_max_iterations: int
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.attempt_id) is not str or not self.attempt_id.strip():
            raise ThroughBedK2KernelError("kernel attempt ID must be nonblank")
        if type(self.prior) is not KernelAcceptedThroughBedK2State:
            raise TypeError("kernel prior has a foreign exact type")
        if self.prior.contract != reference_through_bed_k2_contract(
            self.prior.tray.physical_tray_id
        ):
            raise ThroughBedK2KernelError("kernel prior tray and reference contract disagree")
        if self.bottom_gas_binding.physical_tray_id != self.prior.tray.physical_tray_id:
            raise ThroughBedK2KernelError("gas-boundary binding belongs to another tray")
        if (
            type(self.layer_models) is not tuple
            or len(self.layer_models) != THROUGH_BED_LAYER_COUNT
            or any(type(item) is not legacy.SP1KLayerEngineeringModel for item in self.layer_models)
        ):
            raise ThroughBedK2KernelError("kernel requires two exact K-layer models")
        if self.layer_model_configuration_digest != self.prior.layer_model_configuration_digest:
            raise ThroughBedK2KernelError("accepted K-layer model digest drifted")
        if type(self.bottom_gas_boundary) is not ef.EngineeringGasBoundary:
            raise TypeError("kernel lower gas boundary has a foreign exact type")
        if self.bottom_gas_boundary.carrier_topology is not ef.CarrierTopology.BINARY_NO_INERT:
            raise ThroughBedK2KernelError("through-bed K2 kernel requires Law 2")
        if type(self.flow_segments) is not tuple or not self.flow_segments:
            raise ThroughBedK2KernelError("kernel requires nonempty accepted-flow history")
        if self.flow_segments[0].start_time_s != self.prior.tray.time_s:
            raise ThroughBedK2KernelError("kernel flow history starts at another time")
        if self.flow_segments[-1].end_time_s != self.end_time_s:
            raise ThroughBedK2KernelError("kernel flow history ends at another time")
        _validate_through_bed_layer_geometry(self.prior.contract, self.layer_models)
        _validate_one_gas_enthalpy_convention(self.layer_models)


def _validate_one_gas_enthalpy_convention(
    models: tuple[legacy.SP1KLayerEngineeringModel, legacy.SP1KLayerEngineeringModel],
) -> None:
    if models[0].properties != models[1].properties:
        raise ThroughBedK2KernelError(
            "both K cells must bind one exact gas-property and energy-datum convention"
        )


def _validate_through_bed_layer_geometry(
    contract: ReferenceThroughBedK2Contract,
    models: tuple[legacy.SP1KLayerEngineeringModel, legacy.SP1KLayerEngineeringModel],
) -> None:
    """Require an exact K2 partition with one physical lower floor."""

    reference = _reference_tray(contract.physical_tray_id)
    tray_type = reference_tray_type(reference)
    expected_ids = tuple(
        f"{contract.physical_tray_id}:K{layer}" for layer in range(1, THROUGH_BED_LAYER_COUNT + 1)
    )
    floors: list[cc.TaperedBoreFloorPassageElement] = []
    for expected_id, model in zip(expected_ids, models, strict=True):
        tray = model.geometry.tray
        if tray.tray_id != expected_id:
            raise ThroughBedK2KernelError(
                f"through-bed layer geometry must be explicitly identified as {expected_id}"
            )
        if model.geometry.tray_type != tray_type:
            raise ThroughBedK2KernelError(
                f"{contract.physical_tray_id} K layers require their inherited tray type"
            )
        if tray.role != reference.role or tray.contact is not reference.contact:
            raise ThroughBedK2KernelError("K layer changed its reference role or contact mode")
        if tray.diameter_m != reference.diameter_m:
            raise ThroughBedK2KernelError("K layer changed its physical tray diameter")
        bed = model.hydraulics.series.bed_element
        if bed is None or bed.flow_length_m != tray.loaded_depth_m:
            raise ThroughBedK2KernelError(
                "each K-layer Ergun length must equal that layer's loaded depth"
            )
        if bed.resistance_scale != 1.0:
            raise ThroughBedK2KernelError(
                "each K layer must retain its active physical bed resistance"
            )
        floor_rows = tuple(
            item
            for item in model.hydraulics.series.elements
            if type(item) is cc.TaperedBoreFloorPassageElement
        )
        if len(floor_rows) != 1:
            raise ThroughBedK2KernelError(
                "each typed K series must expose exactly one floor placeholder"
            )
        floors.append(floor_rows[0])
    if tuple(item.resistance_scale for item in floors) != (1.0, 0.0):
        raise ThroughBedK2KernelError(
            "K=2 requires one active physical floor at the lower boundary and an exact "
            "zero-resistance upper placeholder"
        )
    if replace(floors[0], resistance_scale=0.0) != floors[1]:
        raise ThroughBedK2KernelError(
            "inactive upper floor must be the zero-resistance copy of the physical floor"
        )
    if math.fsum(model.geometry.tray.loaded_depth_m for model in models) != (
        reference.loaded_depth_m
    ):
        raise ThroughBedK2KernelError("K layers must partition reference bed depth exactly")
    if math.fsum(model.geometry.tray.indirect_duty_w for model in models) != (
        reference.indirect_duty_w
    ):
        raise ThroughBedK2KernelError(
            "K layers must partition record-only reference indirect duty exactly"
        )
    reference_geometry = cc.CellGeometry(tray=reference, tray_type=tray_type)
    if math.fsum(model.geometry.gas_side_reference_area_m2 for model in models) != (
        reference_geometry.gas_side_reference_area_m2
    ):
        raise ThroughBedK2KernelError("K layers must partition gas-side area exactly")


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedK2KernelRejection:
    attempt_id: str
    reason: str
    rollback_state: KernelAcceptedThroughBedK2State

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedThroughBedK2KernelStep:
    request: ThroughBedK2KernelRequest
    candidate: KernelAcceptedThroughBedK2State
    host_transaction: tray_host.ValidatedKCellTrayStep
    layer_inventories: tuple[legacy.PacketDerivedLayerInventory, ...]
    cell_inputs: tuple[ef.EngineeringCellInputs, ...]
    cell_macro_steps: tuple[ef.EngineeringCellMacroStep, ...]
    layer_transfers: tuple[legacy.LayerAcceptedTransfers, ...]
    inter_k_face: legacy.InterKGasFaceLedger
    tray_pressure: legacy.AggregateSP1PressureLedger
    whole_tray_energy: (
        legacy.WholeTrayExternalEnergyLedger | legacy.NativeWholeTrayExternalEnergyLedger
    )
    packet_exchange_ledgers: tuple[legacy.LayerPacketExchangeLedger, ...]
    coupled_ledger: legacy.CoupledTrayLedger
    codec_registry_digest: str
    auditor_identity_digest: str
    auditor_source_identity: str

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False


def evaluate_through_bed_k2_kernel(
    request: ThroughBedK2KernelRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner,
) -> ValidatedThroughBedK2KernelStep | ThroughBedK2KernelRejection:
    """Evaluate one normalized K2 tray step without mutating accepted state."""

    if type(request) is not ThroughBedK2KernelRequest:
        raise TypeError("through-bed kernel requires an exact normalized request")
    try:
        if type(codec) is not high_loading.HighLoadingV1PacketAdapter and (
            type(codec) is not dual_variant.EngineeringDualVariantPacketAdapter
        ):
            # TAG-2 (owner-ruled GO 2026-09-01): the dual-variant adapter is
            # the second admitted exact type - same closed-union pattern the
            # native-datum V2 gates use; nothing else is admitted.
            raise ThroughBedK2KernelError(
                "through-bed kernel requires the exact detached tag-1 adapter "
                "or the owner-ruled dual-variant adapter"
            )
        planner_type = type(k_planner)
        if planner_type not in (
            tray_host.NoKRekeyPlanner,
            k_scheduler.EngineeringPhysicalKRekeyPlanner,
        ):
            raise ThroughBedK2KernelError(
                "the fixed-pressure tag-1 slice admits only the exact zero-motion "
                "or replay-validated engineering K-event planner"
            )
        if planner_type is k_scheduler.EngineeringPhysicalKRekeyPlanner:
            schedule_step = k_planner.step
            if schedule_step.prior.time_s != request.prior.tray.time_s:
                raise ThroughBedK2KernelError(
                    "engineering K planner prior time differs from accepted tray"
                )
            if schedule_step.candidate.time_s != request.end_time_s:
                raise ThroughBedK2KernelError(
                    "engineering K planner endpoint differs from tray endpoint"
                )
            if schedule_step.prior.physical_tray_id != request.prior.tray.physical_tray_id:
                raise ThroughBedK2KernelError("engineering K planner belongs to another tray")
        for name, value in (
            ("codec registry", getattr(codec, "codec_registry_digest", None)),
            ("auditor identity", getattr(codec, "auditor_identity_digest", None)),
            ("auditor source", getattr(codec, "source_identity", None)),
        ):
            if type(value) is not str or not value.strip():
                raise ThroughBedK2KernelError(f"{name} is missing or blank")
        if codec.codec_registry_digest != request.prior.tray.codec_registry_digest:
            raise ThroughBedK2KernelError("codec violates host registry pin")
        if codec.auditor_identity_digest != request.prior.tray.auditor_identity_digest:
            raise ThroughBedK2KernelError("codec violates host auditor pin")
        codec.require_host_binding(
            configuration_digest=request.prior.tray.configuration_digest,
            energy_datum_id=request.prior.tray.energy_datum_id,
            pressure_pa=request.prior.tray.pressure_pa,
        )
        # CELL-02d W2: the declared arm keeps the pinned-adapter drift gate;
        # the native arm's coherence is the datum-identity triple (its own
        # post-init pin, the host gates above, and the codec datum here).
        cell_properties = request.layer_models[0].properties
        if type(cell_properties) is cc.DeclaredGasProperties:
            codec.component_datum_adapter.require_gas_properties(cell_properties)
        elif cell_properties.energy_datum_id != codec.energy_datum_id:
            raise ThroughBedK2KernelError(
                "native cell properties and codec use different energy datums"
            )
        inventories = tuple(
            legacy._layer_inventory(
                request.prior.tray,
                layer,
                model=request.layer_models[layer - 1],
                codec=codec,
            )
            for layer in range(1, THROUGH_BED_LAYER_COUNT + 1)
        )
        inputs, macro_steps, face = legacy._solve_serial_k_cells(request, inventories)
        if not face.passed:
            raise ThroughBedK2KernelError("internal K gas-face ledger did not pass")
        tray_pressure = legacy._aggregate_pressure_ledger(request, inputs, macro_steps)
        if not tray_pressure.passed:
            raise ThroughBedK2KernelError("aggregate one-floor pressure ledger did not pass")
        transfers = tuple(
            legacy._accepted_layer_transfers(
                inputs[layer - 1],
                macro_steps[layer - 1],
                layer=layer,
                datum_adapter=codec.component_datum_adapter,
            )
            for layer in range(1, THROUGH_BED_LAYER_COUNT + 1)
        )
        entries, packet_exchange_ledgers = legacy._packet_entries_after(
            request,
            transfers,
            inventories,
            macro_steps,
            codec=codec,
        )
        if not all(ledger.passed for ledger in packet_exchange_ledgers):
            raise ThroughBedK2KernelError("decoded packet allocation ledger did not pass")
        callback = legacy._AcceptedLaw2PacketAdvance(
            entries=entries,
            residual_by_layer=tuple(
                step.fast.evaluation.scaled_residual_norm for step in macro_steps
            ),
            maximum_scaled_residual=request.newton_tolerance,
        )
        host_request = tray_host.KCellTrayStepRequest(
            attempt_id=request.attempt_id,
            prior=request.prior.tray,
            end_time_s=request.end_time_s,
            flow_segments=request.flow_segments,
            external_cell_transfers=tuple(
                tray_host.CellExternalTransfer(
                    water_to_cell_kg=transfer.gas_boundary_water_to_cell_kg,
                    hexane_to_cell_kg=transfer.gas_boundary_hexane_to_cell_kg,
                    common_datum_energy_to_cell_j=transfer.gas_and_wall_boundary_energy_to_cell_j,
                    energy_datum_id=request.prior.tray.energy_datum_id,
                )
                for transfer in transfers
            ),
            relative_limit=request.relative_limit,
        )
        host_result = tray_host.evaluate_k_cell_tray_step(
            host_request,
            packet_callback=callback,
            k_planner=k_planner,
            auditor=codec,
        )
        if type(host_result) is tray_host.KCellTrayStepRejection:
            raise ThroughBedK2KernelError(f"K-cell host rejected step: {host_result.reason}")
        if type(host_result) is not tray_host.ValidatedKCellTrayStep:
            raise ThroughBedK2KernelError("K-cell host returned a foreign transaction")
        whole_tray_energy = legacy._whole_tray_external_energy_ledger(
            request,
            inputs,
            macro_steps,
            transfers,
            host_result,
            request_codec=codec,
        )
        if not whole_tray_energy.passed:
            raise ThroughBedK2KernelError("whole-tray external-only energy ledger did not pass")
        candidate = KernelAcceptedThroughBedK2State(
            contract=request.prior.contract,
            tray=host_result.candidate,
            layer_model_configuration_digest=request.prior.layer_model_configuration_digest,
            cell_field_states=tuple(step.fast.state for step in macro_steps),
            wall_temperatures_k=tuple(step.wall.wall_temperature_k for step in macro_steps),
        )
        coupled = legacy._coupled_ledger(
            request,
            host_result,
            macro_steps,
            transfers,
            whole_tray_energy,
            tray_pressure,
            packet_exchange_ledgers,
        )
        if not coupled.passed:
            raise ThroughBedK2KernelError("coupled packet/cell/wall ledger did not pass")
        return ValidatedThroughBedK2KernelStep(
            request=request,
            candidate=candidate,
            host_transaction=host_result,
            layer_inventories=inventories,
            cell_inputs=inputs,
            cell_macro_steps=macro_steps,
            layer_transfers=transfers,
            inter_k_face=face,
            tray_pressure=tray_pressure,
            whole_tray_energy=whole_tray_energy,
            packet_exchange_ledgers=packet_exchange_ledgers,
            coupled_ledger=coupled,
            codec_registry_digest=codec.codec_registry_digest,
            auditor_identity_digest=codec.auditor_identity_digest,
            auditor_source_identity=codec.source_identity,
        )
    except (
        ThroughBedK2KernelError,
        legacy.SP1K2Law2IntegrationError,
        ef.EngineeringFeasibilityError,
        high_loading.HighLoadingPacketAdapterError,
        component_datum.ComponentEnergyDatumError,
    ) as error:
        return ThroughBedK2KernelRejection(
            attempt_id=request.attempt_id,
            reason=str(error),
            rollback_state=request.prior,
        )
    except Exception as error:
        return ThroughBedK2KernelRejection(
            attempt_id=request.attempt_id,
            reason=f"through-bed K2 kernel failed: {type(error).__name__}: {error}",
            rollback_state=request.prior,
        )


__all__ = (
    "KernelAcceptedThroughBedK2State",
    "ReferenceThroughBedK2Contract",
    "SUPPORTED_REFERENCE_TRAYS",
    "THROUGH_BED_K2_KERNEL_SCHEMA_ID",
    "THROUGH_BED_LAYER_COUNT",
    "THROUGH_BED_RTD_STAGE_COUNT",
    "ThroughBedGasBoundaryBinding",
    "ThroughBedGasBoundaryKind",
    "ThroughBedK2KernelError",
    "ThroughBedK2KernelRejection",
    "ThroughBedK2KernelRequest",
    "ValidatedThroughBedK2KernelStep",
    "evaluate_through_bed_k2_kernel",
    "reference_through_bed_k2_contract",
)
