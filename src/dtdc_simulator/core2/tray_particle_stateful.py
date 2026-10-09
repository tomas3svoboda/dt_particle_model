"""Stateful particle-to-cell transaction semantics for the T3C host.

This module deliberately imports no resolved-particle implementation.  It
defines a cell/cohort-total, proposal/validation/commit seam around an opaque
particle payload.  The cell owns accepted state; drivers only propose.

All component and energy transfers are signed positive from the represented
particle cohort to the enclosing cell.  Every result is architecture-only and
``physically_qualifying`` remains false until the real F4 physics gates pass.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from typing import ClassVar, Generic, Protocol, TypeAlias, TypeVar


CONSERVATION_RELATIVE_LIMIT = 1.0e-10
_TIME_CLOSURE_ULPS = 8.0
_COMPOSITION_CLOSURE_ULPS = 8.0
_DRY_AIR_MOLAR_MASS_KG_MOL = 0.02896546
_WATER_MOLAR_MASS_KG_MOL = 0.018015268
_HEXANE_MOLAR_MASS_KG_MOL = 0.08617536

PayloadT = TypeVar("PayloadT")
CanonicalPart: TypeAlias = str | int | float | bool


def _require_finite(name: str, *values: float) -> None:
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"{name} must contain only finite values")


def _relative_residual(residual: float, *terms: float) -> float:
    scale = max((abs(term) for term in terms), default=0.0)
    if scale == 0.0:
        return 0.0 if residual == 0.0 else math.inf
    return abs(residual) / scale


def _time_roundoff_s(time_before_s: float, time_after_s: float, duration_s: float) -> float:
    return _TIME_CLOSURE_ULPS * max(
        math.ulp(max(abs(time_before_s), 1.0)),
        math.ulp(max(abs(time_after_s), 1.0)),
        math.ulp(max(abs(duration_s), 1.0)),
    )


def _canonical_part(part: CanonicalPart) -> tuple[bytes, bytes]:
    """Return an unambiguous type tag and canonical byte representation."""

    if isinstance(part, bool):
        return b"b", b"1" if part else b"0"
    if isinstance(part, int):
        return b"i", str(part).encode("ascii")
    if isinstance(part, float):
        if not math.isfinite(part):
            raise ValueError("canonical digest floats must be finite")
        return b"f", part.hex().encode("ascii")
    if isinstance(part, str):
        return b"s", part.encode("utf-8")
    raise TypeError(f"unsupported canonical digest part: {type(part).__name__}")


def _digest(parts: list[CanonicalPart]) -> str:
    """Hash a typed, length-prefixed sequence; delimiters inside strings are safe."""

    digest = hashlib.sha256()
    digest.update(b"DTDC-CANONICAL-DIGEST-V3\0")
    digest.update(len(parts).to_bytes(8, "big"))
    for part in parts:
        tag, payload = _canonical_part(part)
        digest.update(tag)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, order=True)
class VerticalLayerId:
    value: int

    def __post_init__(self) -> None:
        if self.value <= 0:
            raise ValueError("vertical layer ID must be positive")


@dataclass(frozen=True, slots=True, order=True)
class RTDStageId:
    value: int

    def __post_init__(self) -> None:
        if self.value <= 0:
            raise ValueError("RTD stage ID must be positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class PhysicalTrayTopology:
    physical_tray_id: str
    vertical_layer_count: int
    rtd_stage_count: int

    def __post_init__(self) -> None:
        if not self.physical_tray_id:
            raise ValueError("physical_tray_id must not be empty")
        if self.vertical_layer_count <= 0 or self.rtd_stage_count <= 0:
            raise ValueError("tray K and RTD counts must be positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleTopologyAuthority:
    """Host authority that keeps physical tray, K layer, and RTD stage distinct."""

    authority_id: str
    tower_id: str
    trays: tuple[PhysicalTrayTopology, ...]

    def __post_init__(self) -> None:
        if not self.authority_id or not self.tower_id:
            raise ValueError("topology authority and tower IDs must not be empty")
        if not self.trays:
            raise ValueError("topology authority must declare at least one tray")
        identifiers = tuple(tray.physical_tray_id for tray in self.trays)
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("physical tray IDs must be unique")

    @property
    def digest(self) -> str:
        parts: list[CanonicalPart] = [
            "particle-topology-v3",
            self.authority_id,
            self.tower_id,
        ]
        for tray in self.trays:
            parts.extend(
                (
                    tray.physical_tray_id,
                    tray.vertical_layer_count,
                    tray.rtd_stage_count,
                )
            )
        return _digest(parts)

    def validate_owner(self, owner: ParticlePopulationKey) -> None:
        if owner.tower_id != self.tower_id:
            raise ValueError("particle owner belongs to a different tower")
        tray = next(
            (item for item in self.trays if item.physical_tray_id == owner.physical_tray_id),
            None,
        )
        if tray is None:
            raise ValueError("particle owner references an undeclared physical tray")
        if owner.vertical_layer_id.value > tray.vertical_layer_count:
            raise ValueError("particle owner vertical layer exceeds declared K")
        if owner.rtd_stage_id.value > tray.rtd_stage_count:
            raise ValueError("particle owner RTD stage exceeds declared order")


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticlePopulationKey:
    tower_id: str
    physical_tray_id: str
    vertical_layer_id: VerticalLayerId
    rtd_stage_id: RTDStageId
    cohort_id: str

    def __post_init__(self) -> None:
        if not self.tower_id or not self.physical_tray_id or not self.cohort_id:
            raise ValueError("tower, tray, and cohort IDs must not be empty")
        if not isinstance(self.vertical_layer_id, VerticalLayerId):
            raise TypeError("vertical_layer_id must be a VerticalLayerId")
        if not isinstance(self.rtd_stage_id, RTDStageId):
            raise TypeError("rtd_stage_id must be an RTDStageId")


@dataclass(frozen=True, slots=True, kw_only=True)
class CellCohortTotalBasis:
    """Explicit extensive basis; per-particle quantities are forbidden at this seam."""

    basis_id: str
    represented_dry_matter_kg: float
    represented_residual_oil_label_kg: float

    kind: ClassVar[str] = "cell_cohort_total"

    def __post_init__(self) -> None:
        _require_finite(
            "cell-cohort-total basis",
            self.represented_dry_matter_kg,
            self.represented_residual_oil_label_kg,
        )
        if not self.basis_id:
            raise ValueError("basis_id must not be empty")
        if self.represented_dry_matter_kg <= 0.0:
            raise ValueError("represented dry matter must be strictly positive")
        if not 0.0 <= self.represented_residual_oil_label_kg <= (
            self.represented_dry_matter_kg
        ):
            raise ValueError("represented residual-oil label must lie within dry meal")


@dataclass(frozen=True, slots=True)
class ParticleInventories:
    dry_matter_kg: float
    residual_oil_label_kg: float
    water_kg: float
    hexane_kg: float
    common_datum_energy_j: float

    def __post_init__(self) -> None:
        _require_finite(
            "particle inventories",
            self.dry_matter_kg,
            self.residual_oil_label_kg,
            self.water_kg,
            self.hexane_kg,
            self.common_datum_energy_j,
        )
        if self.dry_matter_kg <= 0.0:
            raise ValueError("dry_matter_kg must be strictly positive")
        if not 0.0 <= self.residual_oil_label_kg <= self.dry_matter_kg:
            raise ValueError("residual-oil label must lie within dry meal")
        if self.water_kg < 0.0 or self.hexane_kg < 0.0:
            raise ValueError("water and hexane inventories must be non-negative")


@dataclass(frozen=True, slots=True)
class ParticlePayloadAudit:
    """Host observation of opaque state, including semantic coordinates."""

    payload_digest: str
    time_s: float
    pressure_pa: float
    configuration_digest: str
    energy_datum_id: str
    regime: str
    inventories: ParticleInventories

    def __post_init__(self) -> None:
        _require_finite("payload coordinates", self.time_s, self.pressure_pa)
        if not self.payload_digest:
            raise ValueError("payload_digest must not be empty")
        if self.time_s < 0.0 or self.pressure_pa <= 0.0:
            raise ValueError("payload time/pressure is outside its domain")
        if not self.configuration_digest or not self.energy_datum_id or not self.regime:
            raise ValueError("payload configuration, datum, and regime must not be empty")


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleAuditorPin:
    """Host authority for the opaque-payload observation schema and implementation."""

    auditor_id: str
    schema_id: str
    source_identity: str

    def __post_init__(self) -> None:
        if not self.auditor_id or not self.schema_id or not self.source_identity:
            raise ValueError("auditor, schema, and source identities must not be empty")


def _validate_basis(audit: ParticlePayloadAudit, basis: CellCohortTotalBasis) -> None:
    inventories = audit.inventories
    if inventories.dry_matter_kg != basis.represented_dry_matter_kg:
        raise ValueError("payload dry matter does not match cell-cohort-total basis")
    if inventories.residual_oil_label_kg != basis.represented_residual_oil_label_kg:
        raise ValueError("payload oil label does not match cell-cohort-total basis")


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedParticleState(Generic[PayloadT]):
    owner_key: ParticlePopulationKey
    topology_digest: str
    basis: CellCohortTotalBasis
    revision: int
    time_s: float
    pressure_pa: float
    configuration_digest: str
    energy_datum_id: str
    payload: PayloadT
    payload_audit: ParticlePayloadAudit

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if not isinstance(self.basis, CellCohortTotalBasis):
            raise TypeError("particle extensive basis must be CellCohortTotalBasis")
        _require_finite("accepted particle coordinates", self.time_s, self.pressure_pa)
        if not self.topology_digest:
            raise ValueError("topology_digest must not be empty")
        if self.revision < 0 or self.time_s < 0.0 or self.pressure_pa <= 0.0:
            raise ValueError("accepted particle revision/time/pressure is outside its domain")
        if not self.configuration_digest or not self.energy_datum_id:
            raise ValueError("accepted configuration and datum must not be empty")
        audit = self.payload_audit
        if (
            audit.time_s != self.time_s
            or audit.pressure_pa != self.pressure_pa
            or audit.configuration_digest != self.configuration_digest
            or audit.energy_datum_id != self.energy_datum_id
        ):
            raise ValueError("accepted wrapper does not match observed payload coordinates")
        _validate_basis(audit, self.basis)


class ParticleCellContractError(RuntimeError, Generic[PayloadT]):
    def __init__(self, message: str, rollback_state: AcceptedParticleState[PayloadT]) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedGasFlowAudit:
    """Accepted local gas flow from which superficial velocity is reconstructed."""

    flow_authority_id: str
    water_mass_flow_kg_s: float
    hexane_mass_flow_kg_s: float
    dry_air_mass_flow_kg_s: float
    mixture_density_kg_m3: float
    open_flow_area_m2: float
    direction: int = 1

    def __post_init__(self) -> None:
        _require_finite(
            "accepted gas flow",
            self.water_mass_flow_kg_s,
            self.hexane_mass_flow_kg_s,
            self.dry_air_mass_flow_kg_s,
            self.mixture_density_kg_m3,
            self.open_flow_area_m2,
        )
        if not self.flow_authority_id:
            raise ValueError("flow_authority_id must not be empty")
        if min(
            self.water_mass_flow_kg_s,
            self.hexane_mass_flow_kg_s,
            self.dry_air_mass_flow_kg_s,
        ) < 0.0:
            raise ValueError("accepted component gas flows must be non-negative")
        if self.mixture_density_kg_m3 <= 0.0 or self.open_flow_area_m2 <= 0.0:
            raise ValueError("gas density and open flow area must be strictly positive")
        if self.direction not in (-1, 1):
            raise ValueError("gas-flow direction must be -1 or +1")

    @property
    def total_mass_flow_kg_s(self) -> float:
        return math.fsum(
            (
                self.water_mass_flow_kg_s,
                self.hexane_mass_flow_kg_s,
                self.dry_air_mass_flow_kg_s,
            )
        )

    @property
    def superficial_velocity_m_s(self) -> float:
        magnitude = self.total_mass_flow_kg_s / (
            self.mixture_density_kg_m3 * self.open_flow_area_m2
        )
        return self.direction * magnitude

    @property
    def superficial_speed_m_s(self) -> float:
        return abs(self.superficial_velocity_m_s)

    def mole_fractions(self) -> tuple[float, float, float] | None:
        moles = (
            self.water_mass_flow_kg_s / _WATER_MOLAR_MASS_KG_MOL,
            self.hexane_mass_flow_kg_s / _HEXANE_MOLAR_MASS_KG_MOL,
            self.dry_air_mass_flow_kg_s / _DRY_AIR_MOLAR_MASS_KG_MOL,
        )
        total = math.fsum(moles)
        if total == 0.0:
            return None
        water_moles, hexane_moles, dry_air_moles = moles
        return (
            water_moles / total,
            hexane_moles / total,
            dry_air_moles / total,
        )


def _validate_segment(
    *,
    time_before_s: float,
    time_after_s: float,
    gas_temperature_k: float,
    y_water: float,
    y_hexane: float,
    y_dry_air: float,
    gas_flow: AcceptedGasFlowAudit,
) -> None:
    _require_finite(
        "gas boundary segment",
        time_before_s,
        time_after_s,
        gas_temperature_k,
        y_water,
        y_hexane,
        y_dry_air,
    )
    if time_before_s < 0.0 or time_after_s <= time_before_s:
        raise ValueError("boundary segment must be strictly forward in time")
    if gas_temperature_k <= 0.0:
        raise ValueError("gas_temperature_k must be strictly positive")
    fractions = (y_water, y_hexane, y_dry_air)
    if any(not 0.0 <= value <= 1.0 for value in fractions):
        raise ValueError("gas mole fractions must lie in [0, 1]")
    if abs(math.fsum((*fractions, -1.0))) > (
        _COMPOSITION_CLOSURE_ULPS * math.ulp(1.0)
    ):
        raise ValueError("gas mole fractions must sum to one")
    flow_fractions = gas_flow.mole_fractions()
    if flow_fractions is not None:
        mismatch = max(abs(left - right) for left, right in zip(fractions, flow_fractions))
        if mismatch > 1.0e-12:
            raise ValueError("boundary composition does not match accepted component gas flow")


@dataclass(frozen=True, slots=True, kw_only=True)
class GasBoundarySegment:
    time_before_s: float
    time_after_s: float
    gas_temperature_k: float
    y_water: float
    y_hexane: float
    y_dry_air: float
    gas_flow: AcceptedGasFlowAudit

    def __post_init__(self) -> None:
        _validate_segment(
            time_before_s=self.time_before_s,
            time_after_s=self.time_after_s,
            gas_temperature_k=self.gas_temperature_k,
            y_water=self.y_water,
            y_hexane=self.y_hexane,
            y_dry_air=self.y_dry_air,
            gas_flow=self.gas_flow,
        )

    @property
    def duration_s(self) -> float:
        return self.time_after_s - self.time_before_s


@dataclass(frozen=True, slots=True, kw_only=True)
class BinaryGasBoundarySegment:
    """Real sealed-DT binary boundary; dry-air input cannot be smuggled through."""

    time_before_s: float
    time_after_s: float
    gas_temperature_k: float
    y_water: float
    y_hexane: float
    gas_flow: AcceptedGasFlowAudit

    y_dry_air: ClassVar[float] = 0.0

    def __post_init__(self) -> None:
        if self.gas_flow.dry_air_mass_flow_kg_s != 0.0:
            raise ValueError("binary DT boundary requires exactly zero accepted dry-air flow")
        _validate_segment(
            time_before_s=self.time_before_s,
            time_after_s=self.time_after_s,
            gas_temperature_k=self.gas_temperature_k,
            y_water=self.y_water,
            y_hexane=self.y_hexane,
            y_dry_air=0.0,
            gas_flow=self.gas_flow,
        )

    @property
    def duration_s(self) -> float:
        return self.time_after_s - self.time_before_s


BoundarySegment = GasBoundarySegment | BinaryGasBoundarySegment


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleBoundaryHistory:
    pressure_pa: float
    segments: tuple[BoundarySegment, ...]

    def __post_init__(self) -> None:
        _require_finite("boundary-history pressure", self.pressure_pa)
        if self.pressure_pa <= 0.0:
            raise ValueError("boundary-history pressure must be strictly positive")
        if not self.segments:
            raise ValueError("boundary history must contain at least one segment")
        for previous, current in zip(self.segments, self.segments[1:]):
            if current.time_before_s != previous.time_after_s:
                raise ValueError("boundary-history segments must be exactly contiguous")
        duration_sum = math.fsum(segment.duration_s for segment in self.segments)
        endpoint_duration = math.fsum((self.time_after_s, -self.time_before_s))
        if abs(math.fsum((duration_sum, -endpoint_duration))) > _time_roundoff_s(
            self.time_before_s,
            self.time_after_s,
            endpoint_duration,
        ):
            raise ValueError("boundary-history durations do not telescope")

    @property
    def time_before_s(self) -> float:
        return self.segments[0].time_before_s

    @property
    def time_after_s(self) -> float:
        return self.segments[-1].time_after_s

    @property
    def duration_s(self) -> float:
        return math.fsum(segment.duration_s for segment in self.segments)


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleDriverPin:
    driver_id: str
    source_identity: str
    maximum_scaled_residual: float

    def __post_init__(self) -> None:
        _require_finite("driver residual pin", self.maximum_scaled_residual)
        if not self.driver_id or not self.source_identity:
            raise ValueError("driver and source identities must not be empty")
        if self.maximum_scaled_residual <= 0.0:
            raise ValueError("driver residual pin must be strictly positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleAttemptDiagnostics:
    driver_id: str
    source_identity: str
    residual_contract_passed: bool
    maximum_scaled_residual: float
    accepted_nonlinear_evaluations: int
    rejected_nonlinear_evaluations: int
    message: str = ""

    def __post_init__(self) -> None:
        _require_finite("particle attempt diagnostics", self.maximum_scaled_residual)
        if not self.driver_id or not self.source_identity:
            raise ValueError("driver and source identities must not be empty")
        if self.maximum_scaled_residual < 0.0:
            raise ValueError("maximum_scaled_residual must be non-negative")
        if self.accepted_nonlinear_evaluations < 0 or self.rejected_nonlinear_evaluations < 0:
            raise ValueError("nonlinear evaluation counts must be non-negative")


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleAdvanceRequest(Generic[PayloadT]):
    attempt_id: str
    prior: AcceptedParticleState[PayloadT]
    boundary_history: ParticleBoundaryHistory
    topology: ParticleTopologyAuthority
    driver_pin: ParticleDriverPin
    auditor_pin: ParticleAuditorPin

    def __post_init__(self) -> None:
        if not self.attempt_id:
            raise ValueError("attempt_id must not be empty")
        self.topology.validate_owner(self.prior.owner_key)
        if self.prior.topology_digest != self.topology.digest:
            raise ValueError("accepted state was created under a different topology authority")
        if self.boundary_history.time_before_s != self.prior.time_s:
            raise ValueError("boundary history must start at the exact accepted state time")
        if self.boundary_history.pressure_pa != self.prior.pressure_pa:
            raise ValueError(
                "positive-duration pressure changes are forbidden; use a pressure remap"
            )

    @property
    def request_digest(self) -> str:
        return canonical_advance_request_digest(self)


def _append_audit_identity(
    parts: list[CanonicalPart],
    audit: ParticlePayloadAudit,
) -> None:
    inventories = audit.inventories
    parts.extend(
        (
            "payload-audit-v3",
            audit.payload_digest,
            audit.time_s,
            audit.pressure_pa,
            audit.configuration_digest,
            audit.energy_datum_id,
            audit.regime,
            inventories.dry_matter_kg,
            inventories.residual_oil_label_kg,
            inventories.water_kg,
            inventories.hexane_kg,
            inventories.common_datum_energy_j,
        )
    )


def _append_auditor_pin(
    parts: list[CanonicalPart],
    pin: ParticleAuditorPin,
) -> None:
    parts.extend(("auditor-pin-v1", pin.auditor_id, pin.schema_id, pin.source_identity))


def _append_state_identity(
    parts: list[CanonicalPart],
    state: AcceptedParticleState[object],
) -> None:
    owner = state.owner_key
    parts.extend(
        (
            owner.tower_id,
            owner.physical_tray_id,
            owner.vertical_layer_id.value,
            owner.rtd_stage_id.value,
            owner.cohort_id,
            state.topology_digest,
            state.basis.kind,
            state.basis.basis_id,
            state.basis.represented_dry_matter_kg,
            state.basis.represented_residual_oil_label_kg,
            state.revision,
            state.time_s,
            state.pressure_pa,
            state.configuration_digest,
            state.energy_datum_id,
        )
    )
    _append_audit_identity(parts, state.payload_audit)


def canonical_accepted_particle_state_digest(
    state: AcceptedParticleState[object],
) -> str:
    """Canonical identity of the accepted wrapper, basis, and complete audit."""

    parts: list[CanonicalPart] = ["accepted-particle-state-identity-v1"]
    _append_state_identity(parts, state)
    return _digest(parts)


def canonical_advance_request_digest(request: ParticleAdvanceRequest[object]) -> str:
    parts: list[CanonicalPart] = ["particle-advance-request-v3", request.attempt_id]
    _append_state_identity(parts, request.prior)
    _append_auditor_pin(parts, request.auditor_pin)
    parts.extend(
        (
            request.driver_pin.driver_id,
            request.driver_pin.source_identity,
            request.driver_pin.maximum_scaled_residual,
            request.boundary_history.pressure_pa,
        )
    )
    for index, segment in enumerate(request.boundary_history.segments):
        flow = segment.gas_flow
        parts.extend(
            (
                index,
                type(segment).__name__,
                segment.time_before_s,
                segment.time_after_s,
                segment.gas_temperature_k,
                segment.y_water,
                segment.y_hexane,
                segment.y_dry_air,
                flow.flow_authority_id,
                flow.water_mass_flow_kg_s,
                flow.hexane_mass_flow_kg_s,
                flow.dry_air_mass_flow_kg_s,
                flow.mixture_density_kg_m3,
                flow.open_flow_area_m2,
                flow.direction,
            )
        )
    return _digest(parts)


@dataclass(frozen=True, slots=True)
class ParticleBoundaryTransfer:
    """Signed cell/cohort-total transfer, positive particle -> cell."""

    water_to_cell_kg: float
    hexane_to_cell_kg: float
    common_datum_energy_to_cell_j: float
    energy_datum_id: str

    def __post_init__(self) -> None:
        _require_finite(
            "particle boundary transfer",
            self.water_to_cell_kg,
            self.hexane_to_cell_kg,
            self.common_datum_energy_to_cell_j,
        )
        if not self.energy_datum_id:
            raise ValueError("transfer energy datum must not be empty")


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleBalanceInterval(Generic[PayloadT]):
    time_before_s: float
    time_after_s: float
    boundary_segment_index: int
    before_payload: PayloadT
    before_audit: ParticlePayloadAudit
    after_payload: PayloadT
    after_audit: ParticlePayloadAudit
    transfer: ParticleBoundaryTransfer

    def __post_init__(self) -> None:
        _require_finite("particle balance interval", self.time_before_s, self.time_after_s)
        if self.time_after_s <= self.time_before_s:
            raise ValueError("balance interval must be strictly forward in time")
        if self.boundary_segment_index < 0:
            raise ValueError("boundary_segment_index must be non-negative")

    @property
    def duration_s(self) -> float:
        return self.time_after_s - self.time_before_s


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleEventRecord(Generic[PayloadT]):
    event_time_s: float
    event_kind: str
    regime_before: str
    regime_after: str
    regime_changed: bool
    before_payload: PayloadT
    before_audit: ParticlePayloadAudit
    after_payload: PayloadT
    after_audit: ParticlePayloadAudit

    def __post_init__(self) -> None:
        _require_finite("particle event time", self.event_time_s)
        if not self.event_kind or not self.regime_before or not self.regime_after:
            raise ValueError("event kind and regimes must not be empty")
        if self.regime_changed != (self.regime_before != self.regime_after):
            raise ValueError("event regime flag and regime transition disagree")


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleAdvanceProposal(Generic[PayloadT]):
    request_digest: str
    base_revision: int
    base_payload_digest: str
    candidate_payload: PayloadT
    candidate_audit: ParticlePayloadAudit
    transfer: ParticleBoundaryTransfer
    balance_intervals: tuple[ParticleBalanceInterval[PayloadT], ...]
    accepted_events: tuple[ParticleEventRecord[PayloadT], ...]
    diagnostics: ParticleAttemptDiagnostics


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticleAdvanceRejection:
    request_digest: str
    base_revision: int
    base_payload_digest: str
    failure_class: str
    retryable: bool
    suggested_dt_s: float | None
    diagnostics: ParticleAttemptDiagnostics

    def __post_init__(self) -> None:
        if not self.request_digest or not self.base_payload_digest or not self.failure_class:
            raise ValueError("rejection request/base/failure identity must not be empty")
        if self.suggested_dt_s is not None:
            _require_finite("suggested_dt_s", self.suggested_dt_s)
            if self.suggested_dt_s <= 0.0:
                raise ValueError("suggested_dt_s must be positive when supplied")


ParticleAdvanceOutcome = ParticleAdvanceProposal[PayloadT] | ParticleAdvanceRejection


class StatefulParticleDriver(Protocol[PayloadT]):
    def propose(self, request: ParticleAdvanceRequest[PayloadT]) -> ParticleAdvanceOutcome[PayloadT]:
        ...


class ParticlePayloadAuditor(Protocol[PayloadT]):
    auditor_pin: ParticleAuditorPin

    def audit(self, payload: PayloadT) -> ParticlePayloadAudit:
        ...


@dataclass(frozen=True, slots=True)
class ParticleLedger:
    dry_matter_residual_kg: float
    residual_oil_label_residual_kg: float
    water_residual_kg: float
    hexane_residual_kg: float
    common_datum_energy_residual_j: float
    dry_matter_residual_rel: float
    residual_oil_label_residual_rel: float
    water_residual_rel: float
    hexane_residual_rel: float
    common_datum_energy_residual_rel: float
    relative_limit: float

    @property
    def passed(self) -> bool:
        return max(
            self.dry_matter_residual_rel,
            self.residual_oil_label_residual_rel,
            self.water_residual_rel,
            self.hexane_residual_rel,
            self.common_datum_energy_residual_rel,
        ) <= self.relative_limit


@dataclass(frozen=True, slots=True)
class ValidatedParticleInterval:
    interval_index: int
    boundary_segment_index: int
    time_before_s: float
    time_after_s: float
    regime: str
    before_payload_digest: str
    after_payload_digest: str
    transfer: ParticleBoundaryTransfer
    ledger: ParticleLedger


@dataclass(frozen=True, slots=True)
class ValidatedParticleEvent:
    event_index: int
    event_time_s: float
    event_kind: str
    regime_before: str
    regime_after: str
    regime_changed: bool
    before_payload_digest: str
    after_payload_digest: str
    ledger: ParticleLedger


_PARTICLE_VALIDATION_SEAL = object()


@dataclass(frozen=True, slots=True)
class ValidatedParticleAdvance(Generic[PayloadT]):
    prior: AcceptedParticleState[PayloadT]
    candidate: AcceptedParticleState[PayloadT]
    request_digest: str
    auditor_pin: ParticleAuditorPin
    transfer: ParticleBoundaryTransfer
    ledger: ParticleLedger
    intervals: tuple[ValidatedParticleInterval, ...]
    events: tuple[ValidatedParticleEvent, ...]
    diagnostics: ParticleAttemptDiagnostics
    _seal: object = field(repr=False)

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _PARTICLE_VALIDATION_SEAL:
            raise ValueError("validated particle advance must be issued by the evaluator")

    @property
    def transaction_digest(self) -> str:
        return canonical_validated_particle_transaction_digest(self)


def canonical_validated_particle_transaction_digest(
    particle: ValidatedParticleAdvance[object],
) -> str:
    """Bind request, full accepted candidate identity, and exact boundary transfer."""

    transfer = particle.transfer
    return _digest(
        [
            "validated-particle-transaction-v1",
            particle.request_digest,
            canonical_accepted_particle_state_digest(particle.candidate),
            transfer.water_to_cell_kg,
            transfer.hexane_to_cell_kg,
            transfer.common_datum_energy_to_cell_j,
            transfer.energy_datum_id,
        ]
    )


def _build_ledger(
    before: ParticleInventories,
    after: ParticleInventories,
    transfer: ParticleBoundaryTransfer,
    relative_limit: float,
) -> ParticleLedger:
    dry = math.fsum((after.dry_matter_kg, -before.dry_matter_kg))
    oil = math.fsum((after.residual_oil_label_kg, -before.residual_oil_label_kg))
    water = math.fsum((after.water_kg, -before.water_kg, transfer.water_to_cell_kg))
    hexane = math.fsum((after.hexane_kg, -before.hexane_kg, transfer.hexane_to_cell_kg))
    energy = math.fsum(
        (
            after.common_datum_energy_j,
            -before.common_datum_energy_j,
            transfer.common_datum_energy_to_cell_j,
        )
    )
    return ParticleLedger(
        dry_matter_residual_kg=dry,
        residual_oil_label_residual_kg=oil,
        water_residual_kg=water,
        hexane_residual_kg=hexane,
        common_datum_energy_residual_j=energy,
        dry_matter_residual_rel=_relative_residual(
            dry, before.dry_matter_kg, after.dry_matter_kg
        ),
        residual_oil_label_residual_rel=_relative_residual(
            oil, before.residual_oil_label_kg, after.residual_oil_label_kg
        ),
        water_residual_rel=_relative_residual(
            water, before.water_kg, after.water_kg, transfer.water_to_cell_kg
        ),
        hexane_residual_rel=_relative_residual(
            hexane, before.hexane_kg, after.hexane_kg, transfer.hexane_to_cell_kg
        ),
        common_datum_energy_residual_rel=_relative_residual(
            energy,
            before.common_datum_energy_j,
            after.common_datum_energy_j,
            transfer.common_datum_energy_to_cell_j,
        ),
        relative_limit=relative_limit,
    )


def _zero_transfer(datum: str) -> ParticleBoundaryTransfer:
    return ParticleBoundaryTransfer(0.0, 0.0, 0.0, datum)


def _validate_relative_limit(relative_limit: float) -> None:
    if (
        not math.isfinite(relative_limit)
        or relative_limit <= 0.0
        or relative_limit > CONSERVATION_RELATIVE_LIMIT
    ):
        raise ValueError(
            "relative_limit must be finite, positive, and no looser than "
            f"{CONSERVATION_RELATIVE_LIMIT:.1e}"
        )


def _validate_auditor_pin(
    auditor: ParticlePayloadAuditor[PayloadT],
    expected: ParticleAuditorPin,
    prior: AcceptedParticleState[PayloadT],
) -> None:
    if getattr(auditor, "auditor_pin", None) != expected:
        raise ParticleCellContractError("payload auditor violates host authority pin", prior)


def _audit_payload(
    payload: PayloadT,
    declared: ParticlePayloadAudit,
    auditor: ParticlePayloadAuditor[PayloadT],
    prior: AcceptedParticleState[PayloadT],
    label: str,
) -> ParticlePayloadAudit:
    observed = auditor.audit(payload)
    if observed != declared:
        raise ParticleCellContractError(f"{label} payload audit mismatch", prior)
    return observed


def _validate_observation(
    audit: ParticlePayloadAudit,
    *,
    time_s: float,
    pressure_pa: float,
    configuration_digest: str,
    energy_datum_id: str,
    basis: CellCohortTotalBasis,
    prior: AcceptedParticleState[object],
    label: str,
) -> None:
    if (
        audit.time_s != time_s
        or audit.pressure_pa != pressure_pa
        or audit.configuration_digest != configuration_digest
        or audit.energy_datum_id != energy_datum_id
    ):
        raise ParticleCellContractError(f"{label} semantic coordinates mismatch", prior)
    try:
        _validate_basis(audit, basis)
    except ValueError as exc:
        raise ParticleCellContractError(f"{label} violates cell-cohort-total basis", prior) from exc


def _validate_diagnostics(
    diagnostics: ParticleAttemptDiagnostics,
    pin: ParticleDriverPin,
    prior: AcceptedParticleState[object],
    *,
    require_pass: bool,
) -> None:
    if diagnostics.driver_id != pin.driver_id or diagnostics.source_identity != pin.source_identity:
        raise ParticleCellContractError("particle diagnostics violate driver/source pin", prior)
    if require_pass:
        if not diagnostics.residual_contract_passed:
            raise ParticleCellContractError("particle residual contract did not pass", prior)
        if diagnostics.maximum_scaled_residual > pin.maximum_scaled_residual:
            raise ParticleCellContractError("particle scaled residual exceeds host pin", prior)


def _assert_prior_unchanged(
    prior: AcceptedParticleState[PayloadT],
    auditor: ParticlePayloadAuditor[PayloadT],
    original: ParticlePayloadAudit,
) -> None:
    if auditor.audit(prior.payload) != original:
        raise ParticleCellContractError("particle driver mutated accepted payload", prior)


def _transfer_sum(
    intervals: tuple[ParticleBalanceInterval[object], ...], datum: str
) -> ParticleBoundaryTransfer:
    return ParticleBoundaryTransfer(
        math.fsum(item.transfer.water_to_cell_kg for item in intervals),
        math.fsum(item.transfer.hexane_to_cell_kg for item in intervals),
        math.fsum(item.transfer.common_datum_energy_to_cell_j for item in intervals),
        datum,
    )


def _transfers_match(
    left: ParticleBoundaryTransfer,
    right: ParticleBoundaryTransfer,
    relative_limit: float,
) -> bool:
    if left.energy_datum_id != right.energy_datum_id:
        return False
    for a, b in (
        (left.water_to_cell_kg, right.water_to_cell_kg),
        (left.hexane_to_cell_kg, right.hexane_to_cell_kg),
        (left.common_datum_energy_to_cell_j, right.common_datum_energy_to_cell_j),
    ):
        residual = math.fsum((a, -b))
        if _relative_residual(residual, a, b) > relative_limit:
            return False
    return True


def _validate_interval_and_event_chain(
    request: ParticleAdvanceRequest[PayloadT],
    proposal: ParticleAdvanceProposal[PayloadT],
    auditor: ParticlePayloadAuditor[PayloadT],
    observed_prior: ParticlePayloadAudit,
    observed_candidate: ParticlePayloadAudit,
    relative_limit: float,
) -> tuple[tuple[ValidatedParticleInterval, ...], tuple[ValidatedParticleEvent, ...]]:
    intervals = proposal.balance_intervals
    if not intervals:
        raise ParticleCellContractError("proposal has no accepted balance intervals", request.prior)
    if intervals[0].time_before_s != request.boundary_history.time_before_s:
        raise ParticleCellContractError("balance partition does not start at exact t0", request.prior)
    if intervals[-1].time_after_s != request.boundary_history.time_after_s:
        raise ParticleCellContractError("balance partition does not end at exact t1", request.prior)

    events = proposal.accepted_events
    event_keys: set[tuple[float, str]] = set()
    previous_event_time = request.boundary_history.time_before_s
    for event in events:
        if not request.boundary_history.time_before_s < event.event_time_s <= (
            request.boundary_history.time_after_s
        ):
            raise ParticleCellContractError("accepted event lies outside (t0, t1]", request.prior)
        if event.event_time_s < previous_event_time:
            raise ParticleCellContractError("accepted events are out of order", request.prior)
        key = (event.event_time_s, event.event_kind)
        if key in event_keys:
            raise ParticleCellContractError("duplicate accepted event record", request.prior)
        event_keys.add(key)
        previous_event_time = event.event_time_s

    segment_intervals: dict[int, list[ParticleBalanceInterval[PayloadT]]] = {
        index: [] for index in range(len(request.boundary_history.segments))
    }
    for interval in intervals:
        if interval.boundary_segment_index not in segment_intervals:
            raise ParticleCellContractError("balance interval references unknown forcing segment", request.prior)
        segment = request.boundary_history.segments[interval.boundary_segment_index]
        if (
            interval.time_before_s < segment.time_before_s
            or interval.time_after_s > segment.time_after_s
        ):
            raise ParticleCellContractError("balance interval crosses a forcing discontinuity", request.prior)
        segment_intervals[interval.boundary_segment_index].append(interval)
    for index, segment in enumerate(request.boundary_history.segments):
        owned = segment_intervals[index]
        if not owned:
            raise ParticleCellContractError("forcing segment has no accepted balance interval", request.prior)
        if owned[0].time_before_s != segment.time_before_s:
            raise ParticleCellContractError("forcing-segment balance coverage starts late", request.prior)
        if owned[-1].time_after_s != segment.time_after_s:
            raise ParticleCellContractError("forcing-segment balance coverage ends early", request.prior)
        for previous, current in zip(owned, owned[1:]):
            if current.time_before_s != previous.time_after_s:
                raise ParticleCellContractError("forcing-segment intervals contain a gap", request.prior)

    cursor_time = request.boundary_history.time_before_s
    cursor_audit = observed_prior
    event_index = 0
    validated_intervals: list[ValidatedParticleInterval] = []
    validated_events: list[ValidatedParticleEvent] = []

    for interval_index, interval in enumerate(intervals):
        if interval.time_before_s != cursor_time:
            raise ParticleCellContractError("accepted balance intervals contain a gap or overlap", request.prior)
        before = _audit_payload(
            interval.before_payload,
            interval.before_audit,
            auditor,
            request.prior,
            "balance-before",
        )
        after = _audit_payload(
            interval.after_payload,
            interval.after_audit,
            auditor,
            request.prior,
            "balance-after",
        )
        for audit, time_s, label in (
            (before, interval.time_before_s, "balance-before"),
            (after, interval.time_after_s, "balance-after"),
        ):
            _validate_observation(
                audit,
                time_s=time_s,
                pressure_pa=request.prior.pressure_pa,
                configuration_digest=request.prior.configuration_digest,
                energy_datum_id=request.prior.energy_datum_id,
                basis=request.prior.basis,
                prior=request.prior,
                label=label,
            )
        if before != cursor_audit:
            raise ParticleCellContractError("balance state chain is discontinuous", request.prior)
        if before.regime != after.regime:
            raise ParticleCellContractError("regime changed without an explicit zero-time event", request.prior)
        if interval.transfer.energy_datum_id != request.prior.energy_datum_id:
            raise ParticleCellContractError("balance interval changed energy datum", request.prior)
        ledger = _build_ledger(before.inventories, after.inventories, interval.transfer, relative_limit)
        if not ledger.passed:
            raise ParticleCellContractError("accepted balance interval ledger failed", request.prior)
        validated_intervals.append(
            ValidatedParticleInterval(
                interval_index=interval_index,
                boundary_segment_index=interval.boundary_segment_index,
                time_before_s=interval.time_before_s,
                time_after_s=interval.time_after_s,
                regime=before.regime,
                before_payload_digest=before.payload_digest,
                after_payload_digest=after.payload_digest,
                transfer=interval.transfer,
                ledger=ledger,
            )
        )
        cursor_time = interval.time_after_s
        cursor_audit = after

        while event_index < len(events) and events[event_index].event_time_s == cursor_time:
            event = events[event_index]
            event_before = _audit_payload(
                event.before_payload,
                event.before_audit,
                auditor,
                request.prior,
                "event-before",
            )
            event_after = _audit_payload(
                event.after_payload,
                event.after_audit,
                auditor,
                request.prior,
                "event-after",
            )
            for audit, label in ((event_before, "event-before"), (event_after, "event-after")):
                _validate_observation(
                    audit,
                    time_s=event.event_time_s,
                    pressure_pa=request.prior.pressure_pa,
                    configuration_digest=request.prior.configuration_digest,
                    energy_datum_id=request.prior.energy_datum_id,
                    basis=request.prior.basis,
                    prior=request.prior,
                    label=label,
                )
            if event_before != cursor_audit:
                raise ParticleCellContractError("event state chain is discontinuous", request.prior)
            if event_before.regime != event.regime_before or event_after.regime != event.regime_after:
                raise ParticleCellContractError("event regime chain contradicts payload observations", request.prior)
            event_ledger = _build_ledger(
                event_before.inventories,
                event_after.inventories,
                _zero_transfer(request.prior.energy_datum_id),
                relative_limit,
            )
            if not event_ledger.passed:
                raise ParticleCellContractError("zero-transfer event ledger failed", request.prior)
            validated_events.append(
                ValidatedParticleEvent(
                    event_index=event_index,
                    event_time_s=event.event_time_s,
                    event_kind=event.event_kind,
                    regime_before=event.regime_before,
                    regime_after=event.regime_after,
                    regime_changed=event.regime_changed,
                    before_payload_digest=event_before.payload_digest,
                    after_payload_digest=event_after.payload_digest,
                    ledger=event_ledger,
                )
            )
            cursor_audit = event_after
            event_index += 1

        if event_index < len(events) and events[event_index].event_time_s < cursor_time:
            raise ParticleCellContractError("accepted event lacks an exact partition boundary", request.prior)

    if event_index != len(events):
        raise ParticleCellContractError("accepted event lacks a complete post-event partition", request.prior)
    if cursor_audit != observed_candidate:
        raise ParticleCellContractError("candidate is not the endpoint of the audited state chain", request.prior)

    duration = math.fsum(interval.duration_s for interval in intervals)
    requested = request.boundary_history.duration_s
    if abs(math.fsum((duration, -requested))) > _time_roundoff_s(
        request.boundary_history.time_before_s,
        request.boundary_history.time_after_s,
        requested,
    ):
        raise ParticleCellContractError("accepted intervals fail eight-ULP duration closure", request.prior)
    return tuple(validated_intervals), tuple(validated_events)


def evaluate_particle_advance(
    driver: StatefulParticleDriver[PayloadT],
    auditor: ParticlePayloadAuditor[PayloadT],
    request: ParticleAdvanceRequest[PayloadT],
    *,
    relative_limit: float = CONSERVATION_RELATIVE_LIMIT,
) -> ValidatedParticleAdvance[PayloadT] | ParticleAdvanceRejection:
    _validate_relative_limit(relative_limit)
    prior = request.prior
    _validate_auditor_pin(auditor, request.auditor_pin, prior)
    observed_prior = _audit_payload(
        prior.payload, prior.payload_audit, auditor, prior, "accepted-prior"
    )
    try:
        outcome = driver.propose(request)
    except Exception as exc:
        try:
            _assert_prior_unchanged(prior, auditor, observed_prior)
        except ParticleCellContractError as mutation:
            raise mutation from exc
        raise
    _validate_auditor_pin(auditor, request.auditor_pin, prior)
    _assert_prior_unchanged(prior, auditor, observed_prior)

    if outcome.request_digest != request.request_digest:
        raise ParticleCellContractError("particle outcome is bound to another host request", prior)
    if outcome.base_revision != prior.revision:
        raise ParticleCellContractError("particle outcome used a stale base revision", prior)
    if outcome.base_payload_digest != prior.payload_audit.payload_digest:
        raise ParticleCellContractError("particle outcome used the wrong base payload", prior)
    if isinstance(outcome, ParticleAdvanceRejection):
        _validate_diagnostics(outcome.diagnostics, request.driver_pin, prior, require_pass=False)
        if outcome.suggested_dt_s is not None and outcome.suggested_dt_s >= (
            request.boundary_history.duration_s
        ):
            raise ParticleCellContractError("retry timestep is not shorter than rejected interval", prior)
        return outcome

    proposal = outcome
    _validate_diagnostics(proposal.diagnostics, request.driver_pin, prior, require_pass=True)
    observed_candidate = _audit_payload(
        proposal.candidate_payload,
        proposal.candidate_audit,
        auditor,
        prior,
        "candidate",
    )
    _validate_observation(
        observed_candidate,
        time_s=request.boundary_history.time_after_s,
        pressure_pa=prior.pressure_pa,
        configuration_digest=prior.configuration_digest,
        energy_datum_id=prior.energy_datum_id,
        basis=prior.basis,
        prior=prior,
        label="candidate",
    )
    if proposal.transfer.energy_datum_id != prior.energy_datum_id:
        raise ParticleCellContractError("proposal transfer changed energy datum", prior)
    interval_evidence, event_evidence = _validate_interval_and_event_chain(
        request,
        proposal,
        auditor,
        observed_prior,
        observed_candidate,
        relative_limit,
    )
    interval_total = _transfer_sum(proposal.balance_intervals, prior.energy_datum_id)
    if not _transfers_match(interval_total, proposal.transfer, relative_limit):
        raise ParticleCellContractError("balance intervals do not telescope to transfer", prior)
    ledger = _build_ledger(
        observed_prior.inventories,
        observed_candidate.inventories,
        proposal.transfer,
        relative_limit,
    )
    if not ledger.passed:
        raise ParticleCellContractError("particle macro endpoint ledger failed", prior)

    candidate = AcceptedParticleState(
        owner_key=prior.owner_key,
        topology_digest=prior.topology_digest,
        basis=prior.basis,
        revision=prior.revision + 1,
        time_s=observed_candidate.time_s,
        pressure_pa=observed_candidate.pressure_pa,
        configuration_digest=observed_candidate.configuration_digest,
        energy_datum_id=observed_candidate.energy_datum_id,
        payload=proposal.candidate_payload,
        payload_audit=observed_candidate,
    )
    return ValidatedParticleAdvance(
        prior=prior,
        candidate=candidate,
        request_digest=request.request_digest,
        auditor_pin=request.auditor_pin,
        transfer=proposal.transfer,
        ledger=ledger,
        intervals=interval_evidence,
        events=event_evidence,
        diagnostics=proposal.diagnostics,
        _seal=_PARTICLE_VALIDATION_SEAL,
    )


@dataclass(frozen=True, slots=True)
class CellConservedTotals:
    """Receiving-cell totals excluding the represented particle cohort.

    The particle cohort has its own independently audited inventories.  These
    totals cover the complementary enclosing-cell control volume so that the
    signed particle transfer appears exactly once in the cancellation ledger.
    """

    water_kg: float
    hexane_kg: float
    common_datum_energy_j: float

    def __post_init__(self) -> None:
        _require_finite(
            "cell conserved totals", self.water_kg, self.hexane_kg, self.common_datum_energy_j
        )
        if self.water_kg < 0.0 or self.hexane_kg < 0.0:
            raise ValueError("cell component inventories must be non-negative")


@dataclass(frozen=True, slots=True)
class CellExternalTransfer:
    """Signed transfer into the cell from every boundary except this particle cohort."""

    water_to_cell_kg: float
    hexane_to_cell_kg: float
    common_datum_energy_to_cell_j: float
    energy_datum_id: str

    def __post_init__(self) -> None:
        _require_finite(
            "cell external transfer",
            self.water_to_cell_kg,
            self.hexane_to_cell_kg,
            self.common_datum_energy_to_cell_j,
        )
        if not self.energy_datum_id:
            raise ValueError("cell external transfer datum must not be empty")


@dataclass(frozen=True, slots=True, kw_only=True)
class HostCellConservedState:
    """Host-owned identity and conserved totals for one accepted/proposed cell revision."""

    cell_id: str
    revision: int
    state_digest: str
    energy_datum_id: str
    totals: CellConservedTotals

    def __post_init__(self) -> None:
        if not self.cell_id or not self.state_digest or not self.energy_datum_id:
            raise ValueError("cell, state-digest, and datum identities must not be empty")
        if self.revision < 0:
            raise ValueError("cell revision must be non-negative")


@dataclass(frozen=True, slots=True, kw_only=True)
class CellResidualPin:
    cell_solver_id: str
    source_identity: str
    maximum_scaled_residual: float

    def __post_init__(self) -> None:
        _require_finite("cell residual pin", self.maximum_scaled_residual)
        if (
            not self.cell_solver_id
            or not self.source_identity
            or self.maximum_scaled_residual <= 0.0
        ):
            raise ValueError(
                "cell solver/source identities and positive residual pin are required"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class HostCellTransactionRequest:
    """Immutable declared authority for one cell/particle cancellation attempt.

    Its digest binds exact cell identities, revisions, state digests, totals,
    external transfer, the validated particle transaction, and solver pin.
    The tower must still establish that these declared cell states are its
    authoritative states.
    """

    attempt_id: str
    cell_before: HostCellConservedState
    cell_after: HostCellConservedState
    other_external_transfer: CellExternalTransfer
    particle_transaction_digest: str
    residual_pin: CellResidualPin

    def __post_init__(self) -> None:
        if not self.attempt_id or not self.particle_transaction_digest:
            raise ValueError("cell attempt and particle transaction identities must not be empty")
        if self.cell_before.cell_id != self.cell_after.cell_id:
            raise ValueError("cell transaction cannot change cell identity")
        if self.cell_after.revision != self.cell_before.revision + 1:
            raise ValueError("cell transaction must advance revision exactly once")
        if self.cell_before.energy_datum_id != self.cell_after.energy_datum_id:
            raise ValueError("cell transaction cannot change energy datum")
        if self.other_external_transfer.energy_datum_id != (
            self.cell_before.energy_datum_id
        ):
            raise ValueError("cell external transfer uses a different energy datum")

    @property
    def request_digest(self) -> str:
        before = self.cell_before
        after = self.cell_after
        external = self.other_external_transfer
        pin = self.residual_pin
        parts: list[CanonicalPart] = [
            "host-cell-particle-transaction-v1",
            self.attempt_id,
            before.cell_id,
            before.revision,
            before.state_digest,
            before.energy_datum_id,
            before.totals.water_kg,
            before.totals.hexane_kg,
            before.totals.common_datum_energy_j,
            after.cell_id,
            after.revision,
            after.state_digest,
            after.energy_datum_id,
            after.totals.water_kg,
            after.totals.hexane_kg,
            after.totals.common_datum_energy_j,
            external.water_to_cell_kg,
            external.hexane_to_cell_kg,
            external.common_datum_energy_to_cell_j,
            external.energy_datum_id,
            self.particle_transaction_digest,
            pin.cell_solver_id,
            pin.source_identity,
            pin.maximum_scaled_residual,
        ]
        return _digest(parts)


@dataclass(frozen=True, slots=True, kw_only=True)
class CellAttemptDiagnostics:
    attempt_id: str
    cell_solver_id: str
    source_identity: str
    residual_contract_passed: bool
    maximum_scaled_residual: float

    def __post_init__(self) -> None:
        _require_finite("cell attempt diagnostics", self.maximum_scaled_residual)
        if (
            not self.attempt_id
            or not self.cell_solver_id
            or not self.source_identity
            or self.maximum_scaled_residual < 0.0
        ):
            raise ValueError(
                "cell attempt/solver/source identities and nonnegative residual are required"
            )


@dataclass(frozen=True, slots=True)
class CellLedger:
    water_residual_kg: float
    hexane_residual_kg: float
    common_datum_energy_residual_j: float
    water_residual_rel: float
    hexane_residual_rel: float
    common_datum_energy_residual_rel: float
    relative_limit: float

    @property
    def passed(self) -> bool:
        return max(
            self.water_residual_rel,
            self.hexane_residual_rel,
            self.common_datum_energy_residual_rel,
        ) <= self.relative_limit


_CELL_VALIDATION_SEAL = object()


@dataclass(frozen=True, slots=True)
class CellValidatedParticleTransaction(Generic[PayloadT]):
    """Sealed arithmetic certificate; not a lock or persistence transaction."""

    particle: ValidatedParticleAdvance[PayloadT]
    cell_request: HostCellTransactionRequest
    cell_request_digest: str
    cell_ledger: CellLedger
    diagnostics: CellAttemptDiagnostics
    _seal: object = field(repr=False)

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _CELL_VALIDATION_SEAL:
            raise ValueError("cell transaction token must be issued by the cell validator")


def validate_cell_particle_transaction(
    particle: ValidatedParticleAdvance[PayloadT],
    *,
    cell_request: HostCellTransactionRequest,
    diagnostics: CellAttemptDiagnostics,
    relative_limit: float = CONSERVATION_RELATIVE_LIMIT,
) -> CellValidatedParticleTransaction[PayloadT]:
    _validate_relative_limit(relative_limit)
    if particle._seal is not _PARTICLE_VALIDATION_SEAL:
        raise ValueError("particle proposal lacks evaluator validation")
    if cell_request.particle_transaction_digest != particle.transaction_digest:
        raise ParticleCellContractError(
            "cell transaction is bound to another particle transaction",
            particle.prior,
        )
    datum = particle.candidate.energy_datum_id
    if cell_request.cell_before.energy_datum_id != datum:
        raise ParticleCellContractError("cell and particle energy datums differ", particle.prior)
    pin = cell_request.residual_pin
    if diagnostics.attempt_id != cell_request.attempt_id:
        raise ParticleCellContractError("cell diagnostics violate attempt pin", particle.prior)
    if (
        diagnostics.cell_solver_id != pin.cell_solver_id
        or diagnostics.source_identity != pin.source_identity
    ):
        raise ParticleCellContractError(
            "cell diagnostics violate solver/source pin",
            particle.prior,
        )
    if not diagnostics.residual_contract_passed:
        raise ParticleCellContractError("cell residual contract did not pass", particle.prior)
    if diagnostics.maximum_scaled_residual > pin.maximum_scaled_residual:
        raise ParticleCellContractError("cell scaled residual exceeds host pin", particle.prior)

    transfer = particle.transfer
    cell_before = cell_request.cell_before.totals
    cell_after = cell_request.cell_after.totals
    other_external_transfer = cell_request.other_external_transfer
    water = math.fsum(
        (
            cell_after.water_kg,
            -cell_before.water_kg,
            -transfer.water_to_cell_kg,
            -other_external_transfer.water_to_cell_kg,
        )
    )
    hexane = math.fsum(
        (
            cell_after.hexane_kg,
            -cell_before.hexane_kg,
            -transfer.hexane_to_cell_kg,
            -other_external_transfer.hexane_to_cell_kg,
        )
    )
    energy = math.fsum(
        (
            cell_after.common_datum_energy_j,
            -cell_before.common_datum_energy_j,
            -transfer.common_datum_energy_to_cell_j,
            -other_external_transfer.common_datum_energy_to_cell_j,
        )
    )
    ledger = CellLedger(
        water_residual_kg=water,
        hexane_residual_kg=hexane,
        common_datum_energy_residual_j=energy,
        water_residual_rel=_relative_residual(
            water,
            cell_before.water_kg,
            cell_after.water_kg,
            transfer.water_to_cell_kg,
            other_external_transfer.water_to_cell_kg,
        ),
        hexane_residual_rel=_relative_residual(
            hexane,
            cell_before.hexane_kg,
            cell_after.hexane_kg,
            transfer.hexane_to_cell_kg,
            other_external_transfer.hexane_to_cell_kg,
        ),
        common_datum_energy_residual_rel=_relative_residual(
            energy,
            cell_before.common_datum_energy_j,
            cell_after.common_datum_energy_j,
            transfer.common_datum_energy_to_cell_j,
            other_external_transfer.common_datum_energy_to_cell_j,
        ),
        relative_limit=relative_limit,
    )
    if not ledger.passed:
        raise ParticleCellContractError("receiving-cell cancellation ledger failed", particle.prior)
    return CellValidatedParticleTransaction(
        particle=particle,
        cell_request=cell_request,
        cell_request_digest=cell_request.request_digest,
        cell_ledger=ledger,
        diagnostics=diagnostics,
        _seal=_CELL_VALIDATION_SEAL,
    )


def commit_particle_advance(
    current: AcceptedParticleState[PayloadT],
    transaction: CellValidatedParticleTransaction[PayloadT],
    auditor: ParticlePayloadAuditor[PayloadT],
    *,
    current_cell: HostCellConservedState,
) -> AcceptedParticleState[PayloadT]:
    """Return a commit-eligible particle state after exact pre-commit checks.

    The caller remains responsible for holding the host lock and atomically
    persisting both ``cell_request.cell_after`` and the returned particle state.
    """

    if (
        not isinstance(transaction, CellValidatedParticleTransaction)
        or transaction._seal is not _CELL_VALIDATION_SEAL
    ):
        raise ParticleCellContractError("particle commit lacks cell validation", current)
    proposal = transaction.particle
    if transaction.cell_request.request_digest != transaction.cell_request_digest:
        raise ParticleCellContractError("cell transaction identity changed", current)
    if transaction.cell_request.particle_transaction_digest != proposal.transaction_digest:
        raise ParticleCellContractError(
            "cell certificate has foreign particle transaction",
            current,
        )
    if current_cell != transaction.cell_request.cell_before:
        raise ParticleCellContractError("stale or foreign receiving-cell transaction", current)
    if current is not proposal.prior:
        raise ParticleCellContractError("stale or foreign particle transaction", current)
    _validate_auditor_pin(auditor, proposal.auditor_pin, current)
    if auditor.audit(current.payload) != current.payload_audit:
        raise ParticleCellContractError("accepted payload changed before commit", current)
    if auditor.audit(proposal.candidate.payload) != proposal.candidate.payload_audit:
        raise ParticleCellContractError("candidate payload changed before commit", current)
    if current.owner_key != proposal.candidate.owner_key:
        raise ParticleCellContractError("particle owner changed", current)
    if current.basis != proposal.candidate.basis:
        raise ParticleCellContractError("cell-cohort-total basis changed", current)
    if proposal.candidate.revision != current.revision + 1:
        raise ParticleCellContractError("particle revision did not advance exactly once", current)
    return proposal.candidate


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticlePressureRemapRequest(Generic[PayloadT]):
    attempt_id: str
    prior: AcceptedParticleState[PayloadT]
    target_pressure_pa: float
    target_configuration_digest: str
    topology: ParticleTopologyAuthority
    driver_pin: ParticleDriverPin
    auditor_pin: ParticleAuditorPin

    def __post_init__(self) -> None:
        _require_finite("pressure-remap target", self.target_pressure_pa)
        if not self.attempt_id or not self.target_configuration_digest:
            raise ValueError("pressure-remap attempt/configuration identity is required")
        if self.target_pressure_pa <= 0.0 or self.target_pressure_pa == self.prior.pressure_pa:
            raise ValueError("pressure remap requires a distinct positive target pressure")
        self.topology.validate_owner(self.prior.owner_key)
        if self.topology.digest != self.prior.topology_digest:
            raise ValueError("pressure remap uses a different topology authority")

    @property
    def request_digest(self) -> str:
        parts: list[CanonicalPart] = [
            "particle-pressure-remap-request-v3",
            self.attempt_id,
        ]
        _append_state_identity(parts, self.prior)
        _append_auditor_pin(parts, self.auditor_pin)
        parts.extend(
            (
                self.target_pressure_pa,
                self.target_configuration_digest,
                self.driver_pin.driver_id,
                self.driver_pin.source_identity,
                self.driver_pin.maximum_scaled_residual,
            )
        )
        return _digest(parts)


@dataclass(frozen=True, slots=True, kw_only=True)
class ParticlePressureRemapProposal(Generic[PayloadT]):
    request_digest: str
    base_revision: int
    base_payload_digest: str
    target_pressure_pa: float
    target_configuration_digest: str
    candidate_payload: PayloadT
    candidate_audit: ParticlePayloadAudit
    diagnostics: ParticleAttemptDiagnostics


_PRESSURE_VALIDATION_SEAL = object()


@dataclass(frozen=True, slots=True)
class ValidatedPressureRemap(Generic[PayloadT]):
    prior: AcceptedParticleState[PayloadT]
    candidate: AcceptedParticleState[PayloadT]
    request_digest: str
    auditor_pin: ParticleAuditorPin
    diagnostics: ParticleAttemptDiagnostics
    _seal: object = field(repr=False)

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _PRESSURE_VALIDATION_SEAL:
            raise ValueError("validated pressure remap must be issued by the remap validator")


def validate_pressure_remap(
    request: ParticlePressureRemapRequest[PayloadT],
    proposal: ParticlePressureRemapProposal[PayloadT],
    auditor: ParticlePayloadAuditor[PayloadT],
) -> ValidatedPressureRemap[PayloadT]:
    prior = request.prior
    _validate_auditor_pin(auditor, request.auditor_pin, prior)
    observed_prior = _audit_payload(
        prior.payload, prior.payload_audit, auditor, prior, "pressure-remap prior"
    )
    if proposal.request_digest != request.request_digest:
        raise ParticleCellContractError("pressure remap is bound to another request", prior)
    if proposal.base_revision != prior.revision:
        raise ParticleCellContractError("pressure remap used a stale revision", prior)
    if proposal.base_payload_digest != prior.payload_audit.payload_digest:
        raise ParticleCellContractError("pressure remap used the wrong base payload", prior)
    if proposal.target_pressure_pa != request.target_pressure_pa:
        raise ParticleCellContractError("pressure remap changed the host target pressure", prior)
    if proposal.target_configuration_digest != request.target_configuration_digest:
        raise ParticleCellContractError("pressure remap changed the host target configuration", prior)
    _validate_diagnostics(proposal.diagnostics, request.driver_pin, prior, require_pass=True)
    observed_candidate = _audit_payload(
        proposal.candidate_payload,
        proposal.candidate_audit,
        auditor,
        prior,
        "pressure-remap candidate",
    )
    _validate_observation(
        observed_candidate,
        time_s=prior.time_s,
        pressure_pa=request.target_pressure_pa,
        configuration_digest=request.target_configuration_digest,
        energy_datum_id=prior.energy_datum_id,
        basis=prior.basis,
        prior=prior,
        label="pressure-remap candidate",
    )
    if observed_candidate.regime != observed_prior.regime:
        raise ParticleCellContractError(
            "pressure remap changed regime without an explicit ruling",
            prior,
        )
    if observed_candidate.inventories != observed_prior.inventories:
        raise ParticleCellContractError(
            "pressure remap did not preserve all extensive inventories exactly", prior
        )
    candidate = AcceptedParticleState(
        owner_key=prior.owner_key,
        topology_digest=prior.topology_digest,
        basis=prior.basis,
        revision=prior.revision + 1,
        time_s=prior.time_s,
        pressure_pa=request.target_pressure_pa,
        configuration_digest=request.target_configuration_digest,
        energy_datum_id=prior.energy_datum_id,
        payload=proposal.candidate_payload,
        payload_audit=observed_candidate,
    )
    return ValidatedPressureRemap(
        prior=prior,
        candidate=candidate,
        request_digest=request.request_digest,
        auditor_pin=request.auditor_pin,
        diagnostics=proposal.diagnostics,
        _seal=_PRESSURE_VALIDATION_SEAL,
    )


def commit_pressure_remap(
    current: AcceptedParticleState[PayloadT],
    remap: ValidatedPressureRemap[PayloadT],
    auditor: ParticlePayloadAuditor[PayloadT],
) -> AcceptedParticleState[PayloadT]:
    if not isinstance(remap, ValidatedPressureRemap) or remap._seal is not (
        _PRESSURE_VALIDATION_SEAL
    ):
        raise ParticleCellContractError("pressure remap lacks host validation", current)
    if current is not remap.prior:
        raise ParticleCellContractError("stale or foreign pressure remap", current)
    _validate_auditor_pin(auditor, remap.auditor_pin, current)
    if auditor.audit(current.payload) != current.payload_audit:
        raise ParticleCellContractError("pressure-remap prior changed before commit", current)
    if auditor.audit(remap.candidate.payload) != remap.candidate.payload_audit:
        raise ParticleCellContractError("pressure-remap candidate changed before commit", current)
    if remap.candidate.revision != current.revision + 1:
        raise ParticleCellContractError("pressure-remap revision did not advance once", current)
    return remap.candidate


__all__ = [
    "AcceptedGasFlowAudit",
    "AcceptedParticleState",
    "BinaryGasBoundarySegment",
    "CellAttemptDiagnostics",
    "CellCohortTotalBasis",
    "CellConservedTotals",
    "CellExternalTransfer",
    "CellResidualPin",
    "CellValidatedParticleTransaction",
    "CONSERVATION_RELATIVE_LIMIT",
    "GasBoundarySegment",
    "HostCellConservedState",
    "HostCellTransactionRequest",
    "ParticleAdvanceProposal",
    "ParticleAdvanceRejection",
    "ParticleAdvanceRequest",
    "ParticleAttemptDiagnostics",
    "ParticleAuditorPin",
    "ParticleBalanceInterval",
    "ParticleBoundaryHistory",
    "ParticleBoundaryTransfer",
    "ParticleCellContractError",
    "ParticleDriverPin",
    "ParticleEventRecord",
    "ParticleInventories",
    "ParticlePayloadAudit",
    "ParticlePopulationKey",
    "ParticlePressureRemapProposal",
    "ParticlePressureRemapRequest",
    "ParticleTopologyAuthority",
    "PhysicalTrayTopology",
    "RTDStageId",
    "StatefulParticleDriver",
    "ValidatedParticleAdvance",
    "ValidatedPressureRemap",
    "VerticalLayerId",
    "canonical_accepted_particle_state_digest",
    "canonical_advance_request_digest",
    "canonical_validated_particle_transaction_digest",
    "commit_particle_advance",
    "commit_pressure_remap",
    "evaluate_particle_advance",
    "validate_cell_particle_transaction",
    "validate_pressure_remap",
]
