"""Generic nonqualifying K=2 / RTD=8 Law-2 reference through-bed tray seam.

MN1, MN2, and SP1 share the same family-2 through-bed floor geometry but retain
their exact reference depths, recorded indirect duties, roles, and gas-source
provenance.  This public seam normalizes all three into the shared kernel while
the legacy SP1 module retains its existing public API and digest domain.

The 140/280 kW reference duties are record-only configuration values.  Heat is
still supplied by the declared wall steam-side closure.  Inter-tray gas equality,
authenticated boundary histories, scheduler reconciliation after population
routing, durable atomic persistence, physical qualification, calibration, and
plant prediction are not implemented.
"""

from __future__ import annotations

import hashlib
import math
import struct
from dataclasses import dataclass, field
from typing import ClassVar

from . import cell_engineering_feasibility as ef
from . import high_loading_packet_adapter as high_loading
from . import k_cell_tray_host as tray_host
from . import sp1_k2_law2_tray_integration as legacy
from . import through_bed_k2_law2_kernel as kernel


THROUGH_BED_K2_INTEGRATION_DIGEST_DOMAIN = "GT-PS-2/through-bed-k2-rtd8-law2-tray-integration/v1"

ThroughBedKLayerEngineeringModel = legacy.SP1KLayerEngineeringModel
InterKGasFaceLedger = legacy.InterKGasFaceLedger
AggregateThroughBedPressureLedger = legacy.AggregateSP1PressureLedger
PacketDerivedLayerInventory = legacy.PacketDerivedLayerInventory
LayerAcceptedTransfers = legacy.LayerAcceptedTransfers
LayerPacketExchangeLedger = legacy.LayerPacketExchangeLedger
WholeTrayExternalEnergyLedger = legacy.WholeTrayExternalEnergyLedger
CoupledTrayLedger = legacy.CoupledTrayLedger

ReferenceThroughBedK2Contract = kernel.ReferenceThroughBedK2Contract
ThroughBedGasBoundaryBinding = kernel.ThroughBedGasBoundaryBinding
ThroughBedGasBoundaryKind = kernel.ThroughBedGasBoundaryKind
reference_through_bed_k2_contract = kernel.reference_through_bed_k2_contract


class ThroughBedK2Law2IntegrationError(ValueError):
    """Base refusal for a malformed generic through-bed transaction."""


class ThroughBedK2Law2ConfigurationError(ThroughBedK2Law2IntegrationError):
    """Reference geometry, state, boundary, or closure declarations disagree."""


class ThroughBedK2Law2StepError(ThroughBedK2Law2IntegrationError):
    """A generic through-bed macro-step is inadmissible."""


class StaleThroughBedK2Law2StateError(ThroughBedK2Law2IntegrationError):
    """Commit was attempted against another accepted state."""


def _require_finite(name: str, value: float, *, positive: bool = False) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise ThroughBedK2Law2ConfigurationError(f"{name} must be a finite exact binary64 value")
    if positive and value <= 0.0:
        raise ThroughBedK2Law2ConfigurationError(f"{name} must be strictly positive")


def _require_sha256(name: str, value: str) -> None:
    if (
        type(value) is not str
        or len(value) != 71
        or not value.startswith("sha256:")
        or any(character not in "0123456789abcdef" for character in value[7:])
    ):
        raise ThroughBedK2Law2ConfigurationError(f"{name} must be a canonical sha256 identity")


def canonical_through_bed_k2_layer_model_digest(
    contract: ReferenceThroughBedK2Contract,
    models: tuple[ThroughBedKLayerEngineeringModel, ThroughBedKLayerEngineeringModel],
) -> str:
    """Bind one tray contract and both layer models in canonical K order."""

    if type(contract) is not ReferenceThroughBedK2Contract:
        raise ThroughBedK2Law2ConfigurationError("model digest requires an exact tray contract")
    if (
        type(models) is not tuple
        or len(models) != kernel.THROUGH_BED_LAYER_COUNT
        or any(type(item) is not ThroughBedKLayerEngineeringModel for item in models)
    ):
        raise ThroughBedK2Law2ConfigurationError(
            "model digest requires two exact layer models in canonical K order"
        )
    try:
        kernel._validate_through_bed_layer_geometry(contract, models)
        kernel._validate_one_gas_enthalpy_convention(models)
    except kernel.ThroughBedK2KernelError as error:
        raise ThroughBedK2Law2ConfigurationError(str(error)) from error
    digest = hashlib.sha256()
    for part in (
        (THROUGH_BED_K2_INTEGRATION_DIGEST_DOMAIN + "/layer-models").encode("ascii"),
        contract.definition_digest.encode("ascii"),
        *(repr(model).encode("utf-8") for model in models),
    ):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedThroughBedK2Law2TrayState:
    """One accepted single-tray host, two fast seeds, and two wall nodes."""

    contract: ReferenceThroughBedK2Contract
    tray: tray_host.AcceptedKCellTrayState
    layer_model_configuration_digest: str
    cell_field_states: tuple[ef.BinaryNoInertCellFieldState, ef.BinaryNoInertCellFieldState]
    wall_temperatures_k: tuple[float, float]

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    complete_intertray_model: ClassVar[bool] = False
    complete_particle_thermodynamic_state_coupling: ClassVar[bool] = False
    reference_indirect_duty_is_record_only: ClassVar[bool] = True
    reference_indirect_duty_realized: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.contract) is not ReferenceThroughBedK2Contract:
            raise ThroughBedK2Law2ConfigurationError("state contract has a foreign exact type")
        if type(self.tray) is not tray_host.AcceptedKCellTrayState:
            raise ThroughBedK2Law2ConfigurationError("tray must be an accepted K-cell host state")
        if self.tray.physical_tray_id != self.contract.physical_tray_id:
            raise ThroughBedK2Law2ConfigurationError("accepted tray and contract IDs disagree")
        expected_contract = reference_through_bed_k2_contract(self.tray.physical_tray_id)
        if self.contract != expected_contract:
            raise ThroughBedK2Law2ConfigurationError("accepted tray contract is not canonical")
        _require_sha256(
            "accepted layer-model configuration",
            self.layer_model_configuration_digest,
        )
        physical = next(
            item
            for item in self.tray.topology.trays
            if item.physical_tray_id == self.tray.physical_tray_id
        )
        if (
            physical.vertical_layer_count != kernel.THROUGH_BED_LAYER_COUNT
            or physical.rtd_stage_count != kernel.THROUGH_BED_RTD_STAGE_COUNT
        ):
            raise ThroughBedK2Law2ConfigurationError(
                "generic through-bed integration requires exactly K=2 and RTD=8"
            )
        if type(self.cell_field_states) is not tuple or len(self.cell_field_states) != 2:
            raise ThroughBedK2Law2ConfigurationError("accepted state requires two fast states")
        # D9-e: the accepted fast states may now be the TEN-unknown dry-shell
        # state as well as the shipped binary one, so the gate tests exact
        # membership in the cell's own exported tuple rather than one name.
        # Exactly as before it is an EXACT-type test: an ``isinstance`` here
        # would admit any future subclass silently, and the dry-shell state is
        # not a subclass of the binary one at all, so only naming it admits it.
        # A genuinely foreign type still refuses typed, on this same line.
        if any(
            type(item) not in ef.ACCEPTED_FAST_BLOCK_STATE_TYPES for item in self.cell_field_states
        ):
            raise ThroughBedK2Law2ConfigurationError("accepted fast state has a foreign type")
        if type(self.wall_temperatures_k) is not tuple or len(self.wall_temperatures_k) != 2:
            raise ThroughBedK2Law2ConfigurationError(
                "accepted state requires two wall temperatures"
            )
        for value in self.wall_temperatures_k:
            _require_finite("accepted wall temperature", value, positive=True)

    @property
    def state_digest(self) -> str:
        digest = hashlib.sha256()
        for part in (
            THROUGH_BED_K2_INTEGRATION_DIGEST_DOMAIN.encode("ascii"),
            self.contract.definition_digest.encode("ascii"),
            self.tray.state_digest.encode("ascii"),
            self.layer_model_configuration_digest.encode("ascii"),
            repr(self.cell_field_states).encode("utf-8"),
            *(struct.pack(">d", value) for value in self.wall_temperatures_k),
        ):
            digest.update(len(part).to_bytes(8, "big"))
            digest.update(part)
        return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedK2Law2StepRequest:
    attempt_id: str
    prior: AcceptedThroughBedK2Law2TrayState
    end_time_s: float
    layer_models: tuple[ThroughBedKLayerEngineeringModel, ThroughBedKLayerEngineeringModel]
    bottom_gas_boundary: ef.EngineeringGasBoundary
    bottom_gas_binding: ThroughBedGasBoundaryBinding
    top_boundary_pressure_pa: float
    flow_segments: tuple[tray_host.AcceptedDryMatterFlowSegment, ...]
    newton_tolerance: float = ef.DEFAULT_NEWTON_TOLERANCE
    face_tolerance: float = legacy.DEFAULT_FACE_TOLERANCE
    face_max_iterations: int = legacy.DEFAULT_FACE_MAX_ITERATIONS
    relative_limit: float = legacy.DEFAULT_RELATIVE_LIMIT

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    intertray_gas_equality_proved: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.attempt_id) is not str or not self.attempt_id.strip():
            raise ThroughBedK2Law2ConfigurationError("attempt ID must be nonblank")
        if type(self.prior) is not AcceptedThroughBedK2Law2TrayState:
            raise ThroughBedK2Law2ConfigurationError("prior has a foreign exact type")
        _require_finite("end time", self.end_time_s)
        if self.end_time_s <= self.prior.tray.time_s:
            raise ThroughBedK2Law2ConfigurationError("macro-step must advance accepted time")
        if (
            type(self.layer_models) is not tuple
            or len(self.layer_models) != 2
            or any(type(item) is not ThroughBedKLayerEngineeringModel for item in self.layer_models)
        ):
            raise ThroughBedK2Law2ConfigurationError("request requires two exact layer models")
        if type(self.bottom_gas_boundary) is not ef.EngineeringGasBoundary:
            raise ThroughBedK2Law2ConfigurationError("lower gas boundary has a foreign type")
        if self.bottom_gas_boundary.carrier_topology is not ef.CarrierTopology.BINARY_NO_INERT:
            raise ThroughBedK2Law2ConfigurationError("generic through-bed seam requires Law 2")
        if type(self.bottom_gas_binding) is not ThroughBedGasBoundaryBinding:
            raise ThroughBedK2Law2ConfigurationError("gas-boundary binding has a foreign type")
        if self.bottom_gas_binding.physical_tray_id != self.prior.tray.physical_tray_id:
            raise ThroughBedK2Law2ConfigurationError("gas-boundary binding belongs to another tray")
        _require_finite("top boundary pressure", self.top_boundary_pressure_pa, positive=True)
        _require_finite("Newton tolerance", self.newton_tolerance, positive=True)
        _require_finite("face tolerance", self.face_tolerance, positive=True)
        _require_finite("relative limit", self.relative_limit, positive=True)
        if type(self.face_max_iterations) is not int or self.face_max_iterations <= 0:
            raise ThroughBedK2Law2ConfigurationError(
                "face iteration limit must be a positive exact integer"
            )
        if type(self.flow_segments) is not tuple or not self.flow_segments:
            raise ThroughBedK2Law2ConfigurationError("accepted-flow history must be nonempty")
        if any(
            type(item) is not tray_host.AcceptedDryMatterFlowSegment for item in self.flow_segments
        ):
            raise ThroughBedK2Law2ConfigurationError("flow history contains a foreign segment")
        if self.flow_segments[0].start_time_s != self.prior.tray.time_s:
            raise ThroughBedK2Law2ConfigurationError("flow history starts at another time")
        if self.flow_segments[-1].end_time_s != self.end_time_s:
            raise ThroughBedK2Law2ConfigurationError("flow history ends at another time")
        observed_digest = canonical_through_bed_k2_layer_model_digest(
            self.prior.contract,
            self.layer_models,
        )
        if observed_digest != self.prior.layer_model_configuration_digest:
            raise ThroughBedK2Law2ConfigurationError(
                "layer-model configuration drifted from accepted state"
            )
        for model in self.layer_models:
            if model.properties.energy_datum_id != self.prior.tray.energy_datum_id:
                raise ThroughBedK2Law2ConfigurationError(
                    "cell model and packet host use different energy datums"
                )


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedK2Law2StepRejection:
    attempt_id: str
    reason: str
    rollback_state: AcceptedThroughBedK2Law2TrayState

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


_VALIDATED_GENERIC_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedThroughBedK2Law2TrayStep:
    request: ThroughBedK2Law2StepRequest
    candidate: AcceptedThroughBedK2Law2TrayState
    host_transaction: tray_host.ValidatedKCellTrayStep
    layer_inventories: tuple[PacketDerivedLayerInventory, ...]
    cell_inputs: tuple[ef.EngineeringCellInputs, ...]
    cell_macro_steps: tuple[ef.EngineeringCellMacroStep, ...]
    layer_transfers: tuple[LayerAcceptedTransfers, ...]
    inter_k_face: InterKGasFaceLedger
    tray_pressure: AggregateThroughBedPressureLedger
    whole_tray_energy: WholeTrayExternalEnergyLedger
    packet_exchange_ledgers: tuple[LayerPacketExchangeLedger, ...]
    coupled_ledger: CoupledTrayLedger
    codec_registry_digest: str
    auditor_identity_digest: str
    auditor_source_identity: str
    _seal: object = field(repr=False, compare=False)

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    complete_intertray_model: ClassVar[bool] = False
    reference_indirect_duty_realized: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _VALIDATED_GENERIC_SEAL:
            raise TypeError("validated generic steps are issued only by the evaluator")


def _kernel_state(
    state: AcceptedThroughBedK2Law2TrayState,
) -> kernel.KernelAcceptedThroughBedK2State:
    return kernel.KernelAcceptedThroughBedK2State(
        contract=state.contract,
        tray=state.tray,
        layer_model_configuration_digest=state.layer_model_configuration_digest,
        cell_field_states=state.cell_field_states,
        wall_temperatures_k=state.wall_temperatures_k,
    )


def evaluate_through_bed_k2_law2_tray_step(
    request: ThroughBedK2Law2StepRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner,
) -> ValidatedThroughBedK2Law2TrayStep | ThroughBedK2Law2StepRejection:
    """Evaluate one MN1, MN2, or SP1 K2 step without mutating accepted state."""

    if type(request) is not ThroughBedK2Law2StepRequest:
        raise TypeError("generic integration requires an exact step request")
    normalized = kernel.ThroughBedK2KernelRequest(
        attempt_id=request.attempt_id,
        prior=_kernel_state(request.prior),
        end_time_s=request.end_time_s,
        layer_models=request.layer_models,
        layer_model_configuration_digest=canonical_through_bed_k2_layer_model_digest(
            request.prior.contract,
            request.layer_models,
        ),
        bottom_gas_boundary=request.bottom_gas_boundary,
        bottom_gas_binding=request.bottom_gas_binding,
        top_boundary_pressure_pa=request.top_boundary_pressure_pa,
        flow_segments=request.flow_segments,
        newton_tolerance=request.newton_tolerance,
        face_tolerance=request.face_tolerance,
        face_max_iterations=request.face_max_iterations,
        relative_limit=request.relative_limit,
    )
    result = kernel.evaluate_through_bed_k2_kernel(
        normalized,
        codec=codec,
        k_planner=k_planner,
    )
    if type(result) is kernel.ThroughBedK2KernelRejection:
        return ThroughBedK2Law2StepRejection(
            attempt_id=request.attempt_id,
            reason=result.reason,
            rollback_state=request.prior,
        )
    if type(result) is not kernel.ValidatedThroughBedK2KernelStep:
        return ThroughBedK2Law2StepRejection(
            attempt_id=request.attempt_id,
            reason="shared through-bed kernel returned a foreign result",
            rollback_state=request.prior,
        )
    candidate = AcceptedThroughBedK2Law2TrayState(
        contract=request.prior.contract,
        tray=result.candidate.tray,
        layer_model_configuration_digest=request.prior.layer_model_configuration_digest,
        cell_field_states=result.candidate.cell_field_states,
        wall_temperatures_k=result.candidate.wall_temperatures_k,
    )
    return ValidatedThroughBedK2Law2TrayStep(
        request=request,
        candidate=candidate,
        host_transaction=result.host_transaction,
        layer_inventories=result.layer_inventories,
        cell_inputs=result.cell_inputs,
        cell_macro_steps=result.cell_macro_steps,
        layer_transfers=result.layer_transfers,
        inter_k_face=result.inter_k_face,
        tray_pressure=result.tray_pressure,
        whole_tray_energy=result.whole_tray_energy,
        packet_exchange_ledgers=result.packet_exchange_ledgers,
        coupled_ledger=result.coupled_ledger,
        codec_registry_digest=result.codec_registry_digest,
        auditor_identity_digest=result.auditor_identity_digest,
        auditor_source_identity=result.auditor_source_identity,
        _seal=_VALIDATED_GENERIC_SEAL,
    )


def commit_through_bed_k2_law2_tray_step(
    current: AcceptedThroughBedK2Law2TrayState,
    transaction: ValidatedThroughBedK2Law2TrayStep,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner | None = None,
) -> AcceptedThroughBedK2Law2TrayState:
    """Return the generic candidate after deterministic commit-time replay."""

    if type(current) is not AcceptedThroughBedK2Law2TrayState:
        raise TypeError("generic commit requires an exact accepted state")
    if (
        type(transaction) is not ValidatedThroughBedK2Law2TrayStep
        or transaction._seal is not _VALIDATED_GENERIC_SEAL
    ):
        raise ThroughBedK2Law2IntegrationError("generic commit requires a validated transaction")
    if current is not transaction.request.prior:
        raise StaleThroughBedK2Law2StateError("generic transaction is stale or foreign")
    if current.state_digest != transaction.request.prior.state_digest:
        raise StaleThroughBedK2Law2StateError("generic prior digest changed before commit")
    if (
        codec.codec_registry_digest != transaction.codec_registry_digest
        or codec.auditor_identity_digest != transaction.auditor_identity_digest
        or codec.source_identity != transaction.auditor_source_identity
    ):
        raise ThroughBedK2Law2IntegrationError("generic commit codec violates transaction pins")
    if k_planner is None:
        k_planner = tray_host.NoKRekeyPlanner()
    reevaluated = evaluate_through_bed_k2_law2_tray_step(
        transaction.request,
        codec=codec,
        k_planner=k_planner,
    )
    if type(reevaluated) is not ValidatedThroughBedK2Law2TrayStep or reevaluated != transaction:
        raise ThroughBedK2Law2IntegrationError(
            "generic transaction differs from a complete deterministic commit-time reevaluation"
        )
    committed_tray = tray_host.commit_k_cell_tray_step(
        current.tray,
        transaction.host_transaction,
        auditor=codec,
    )
    if committed_tray is not transaction.candidate.tray:
        raise ThroughBedK2Law2IntegrationError("host and generic candidates diverged")
    if (
        not transaction.inter_k_face.passed
        or not transaction.tray_pressure.passed
        or not transaction.whole_tray_energy.passed
        or not all(ledger.passed for ledger in transaction.packet_exchange_ledgers)
        or not transaction.coupled_ledger.passed
    ):
        raise ThroughBedK2Law2IntegrationError("generic acceptance ledgers no longer pass")
    return transaction.candidate


__all__ = (
    "AcceptedThroughBedK2Law2TrayState",
    "AggregateThroughBedPressureLedger",
    "CoupledTrayLedger",
    "InterKGasFaceLedger",
    "LayerAcceptedTransfers",
    "LayerPacketExchangeLedger",
    "PacketDerivedLayerInventory",
    "ReferenceThroughBedK2Contract",
    "StaleThroughBedK2Law2StateError",
    "THROUGH_BED_K2_INTEGRATION_DIGEST_DOMAIN",
    "ThroughBedGasBoundaryBinding",
    "ThroughBedGasBoundaryKind",
    "ThroughBedK2Law2ConfigurationError",
    "ThroughBedK2Law2IntegrationError",
    "ThroughBedK2Law2StepError",
    "ThroughBedK2Law2StepRejection",
    "ThroughBedK2Law2StepRequest",
    "ThroughBedKLayerEngineeringModel",
    "ValidatedThroughBedK2Law2TrayStep",
    "WholeTrayExternalEnergyLedger",
    "canonical_through_bed_k2_layer_model_digest",
    "commit_through_bed_k2_law2_tray_step",
    "evaluate_through_bed_k2_law2_tray_step",
    "reference_through_bed_k2_contract",
)
