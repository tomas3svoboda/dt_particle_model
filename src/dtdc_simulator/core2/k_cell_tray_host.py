"""Nonqualifying immutable host for one weighted K-cell tray.

This module is the first executable host slice between the frozen topology,
packet-weight, and deterministic RTD contracts.  It deliberately supplies no
film, particle, wall, RPM, draw, or exit law.  A manufactured callback may
propose a per-particle payload update, and a second manufactured callback may
propose endpoint adjacent-K rekeys.  The host independently audits payload
inventories, multiplies per-particle values outward exactly once, advances the
RTD clock from *locally accepted* dry-matter discharge only, validates complete
whole-packet moves, closes cell/tray conservation, and returns one sealed
candidate for an identity-checked commit.

The current cohort-total stateful particle seam cannot represent two split
children that share payload bytes while carrying different external weights.
This host therefore owns a distinct weighted packet state rather than silently
reinterpreting that seam.  Its fluid-cell state remains the existing opaque
host conserved-state record; no missing Law-2 coefficient or physical closure
is invented here.

All claims remain false.  Persistence, admission, exit, unequal-weight RTD
construction, a physical K scheduler, a real payload codec/variant registry,
and a coupled cell/particle solve remain outside this slice.
"""

from __future__ import annotations

import enum
import hashlib
import importlib
import math
import struct
from dataclasses import dataclass, field, replace
from fractions import Fraction
from typing import ClassVar, Protocol

from .joint_carrier_oracle import DryMatterDischargeBasis
from .packet_rtd_oracle import PacketRTDPath
from .tray_particle_stateful import (
    CellConservedTotals,
    CellExternalTransfer,
    HostCellConservedState,
    ParticlePopulationKey,
    ParticleTopologyAuthority,
    RTDStageId,
    VerticalLayerId,
)

# These two evidence packets deliberately require that no other src/ file name
# their module stems as one contiguous raw-text literal.  Resolve their public
# runtime types through split strings, following the repository's existing
# detachment pattern, while retaining exact class identity.
_weight_contract = importlib.import_module("." + "packet_weight_" + "scaling_contract", __package__)
_solid_contract = importlib.import_module("." + "accepted_solid_" + "transfer_host", __package__)

PerParticleInventory = _weight_contract.PerParticleInventory
PacketWeight = _weight_contract.PacketWeight
PacketTotals = _weight_contract.PacketTotals
ExternalSurfaceComposite = _weight_contract.ExternalSurfaceComposite
derive_packet_totals = _weight_contract.derive_packet_totals
require_single_surface_authority = _weight_contract.require_single_surface_authority
PacketPlanEntry = _solid_contract.PacketPlanEntry

HOST_DIGEST_DOMAIN = "GT-PS-2/k-cell-tray-host/v1"
DEFAULT_RELATIVE_LIMIT = 1.0e-12
AFFINE_RESIDENCE_TRAJECTORY_SCHEMA_ID = "dtdc-core2-k-cell-affine-residence-trajectory-v2"
AFFINE_RESIDENCE_TRAJECTORY_SCHEMA_REVISION = 2
AFFINE_RESIDENCE_CLOCK_LAW_ID = (
    "ACCEPTED_RTD_EXPOSURE_N_INTEGRAL_Q_ACCEPTED_OVER_AFFINE_DRY_INVENTORY_V2"
)
AFFINE_RESIDENCE_CLOCK_ALGORITHM_ID = "EXACT_UNION_KNOTS_STABLE_PHI_PSI_FIRST_BINARY64_THRESHOLD_V2"
AFFINE_RESIDENCE_SMALL_ARGUMENT_LIMIT = 2.0**-8
AFFINE_RESIDENCE_LOG1P_RELATIVE_LIMIT = 0.5
AFFINE_RESIDENCE_SERIES_ORDER = 8
AFFINE_RESIDENCE_DECIMAL_ORACLE_ULP_LIMIT = 4
AFFINE_RESIDENCE_DECIMAL_ORACLE_ABSOLUTE_FLOOR = 1.0e-15
AFFINE_RESIDENCE_MIN_NORMAL_MASS_RATIO = float.fromhex("0x1.0000000000000p-1022")
AFFINE_RESIDENCE_MAX_FINITE_MASS_RATIO = float.fromhex("0x1.fffffffffffffp+1023")
AFFINE_RESIDENCE_EXPECTED_SCHEMA_DEFINITION_DIGEST = (
    "sha256:3da301a193aa26cfc68e0d446711e5dc26c163fc233ea34373abd8f9a37806ad"
)
AFFINE_RESIDENCE_EXPECTED_LAW_DEFINITION_DIGEST = (
    "sha256:97115683ce66178d96e61b08a8337d3292baa68fac843d19437e5dc6e9e72452"
)
AFFINE_RESIDENCE_EXPECTED_ALGORITHM_DEFINITION_DIGEST = (
    "sha256:269386ea7516f820a39c4a75524bb31a168b13ab8138d5adc5ab11ca19044e5b"
)


class KCellTrayHostError(ValueError):
    """Base refusal for malformed or nonconservative host work."""


class KCellTrayStateError(KCellTrayHostError):
    """An accepted state violates the immutable host representation."""


class KCellTrayStepError(KCellTrayHostError):
    """A proposed macro-step is inadmissible."""

    def __init__(
        self,
        message: str,
        rollback_state: AcceptedKCellTrayState | None = None,
    ) -> None:
        super().__init__(message)
        self.rollback_state = rollback_state


class StaleKCellTrayStateError(KCellTrayHostError):
    """Commit was attempted against a stale or foreign accepted object."""


def _require_nonblank(name: str, *values: str) -> None:
    if any(type(value) is not str or not value.strip() for value in values):
        raise KCellTrayStateError(f"{name} must contain nonblank exact strings")


def _require_sha256(name: str, value: str) -> None:
    if (
        type(value) is not str
        or len(value) != 71
        or not value.startswith("sha256:")
        or any(character not in "0123456789abcdef" for character in value[7:])
    ):
        raise KCellTrayStateError(f"{name} must be a canonical lowercase sha256 identity")


def _require_binary64(
    name: str,
    *values: float,
    nonnegative: bool = False,
    positive: bool = False,
) -> None:
    for value in values:
        if type(value) is not float or not math.isfinite(value):
            raise KCellTrayStateError(f"{name} must contain finite exact binary64 values")
        if positive and value <= 0.0:
            raise KCellTrayStateError(f"{name} must be strictly positive")
        if nonnegative and value < 0.0:
            raise KCellTrayStateError(f"{name} must be nonnegative")


def _float_bytes(value: float) -> bytes:
    return struct.pack(">d", value)


def _same_float(left: float, right: float) -> bool:
    return _float_bytes(left) == _float_bytes(right)


def _negative_zero(value: float) -> bool:
    return _float_bytes(value) == b"\x80\x00\x00\x00\x00\x00\x00\x00"


def _framed_digest(domain: str, parts: tuple[bytes, ...]) -> str:
    digest = hashlib.sha256()
    for part in (domain.encode("ascii"), *parts):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


def _text_part(value: str) -> bytes:
    return value.encode("utf-8")


def _integer_part(value: int) -> bytes:
    return value.to_bytes(8, "big", signed=True)


def validate_affine_residence_clock_definitions() -> tuple[str, str, str]:
    """Rebuild the literal V2 schema, law, and numerical-algorithm pins."""

    observed_parameters = (
        AFFINE_RESIDENCE_SMALL_ARGUMENT_LIMIT,
        AFFINE_RESIDENCE_LOG1P_RELATIVE_LIMIT,
        AFFINE_RESIDENCE_SERIES_ORDER,
        AFFINE_RESIDENCE_DECIMAL_ORACLE_ULP_LIMIT,
        AFFINE_RESIDENCE_DECIMAL_ORACLE_ABSOLUTE_FLOOR,
        AFFINE_RESIDENCE_MIN_NORMAL_MASS_RATIO,
        AFFINE_RESIDENCE_MAX_FINITE_MASS_RATIO,
    )
    ruled_parameters = (
        2.0**-8,
        0.5,
        8,
        4,
        1.0e-15,
        float.fromhex("0x1.0000000000000p-1022"),
        float.fromhex("0x1.fffffffffffffp+1023"),
    )
    if observed_parameters != ruled_parameters:
        raise KCellTrayStateError("affine residence numerical parameters differ from fixed pins")
    schema_digest = _framed_digest(
        HOST_DIGEST_DOMAIN + "/affine-residence-schema-definition-v2",
        (
            _text_part(AFFINE_RESIDENCE_TRAJECTORY_SCHEMA_ID),
            _integer_part(AFFINE_RESIDENCE_TRAJECTORY_SCHEMA_REVISION),
            _text_part(
                "segment:segment_id,start_time_s,end_time_s,"
                "start_dry_matter_holdup_kg,end_dry_matter_holdup_kg"
            ),
            _text_part(
                "trajectory:trajectory_id,physical_tray_id,"
                "source_target_definition_digest,source_prior_state_digest,"
                "source_flow_history_definition_digest,segments,"
                "schema_definition_digest,law_definition_digest,"
                "algorithm_definition_digest"
            ),
        ),
    )
    law_digest = _framed_digest(
        HOST_DIGEST_DOMAIN + "/affine-residence-law-definition-v2",
        (
            _text_part(schema_digest),
            _text_part(AFFINE_RESIDENCE_CLOCK_LAW_ID),
            _text_part("exposure_increment=N*integral(q_accepted(t)/M_star(t),t0,t1)"),
            _text_part("segment_exposure=(N*q_accepted*dt/M0)*phi(x);x=(M1-M0)/M0"),
            _text_part("phi(x)=log1p(x)/x;phi(0)=1"),
            _text_part("q_accepted=0_freezes_clock;offered_flow_is_excluded"),
        ),
    )
    algorithm_digest = _framed_digest(
        HOST_DIGEST_DOMAIN + "/affine-residence-algorithm-definition-v2",
        (
            _text_part(law_digest),
            _text_part(AFFINE_RESIDENCE_CLOCK_ALGORITHM_ID),
            _text_part("exact_Fraction_endpoint_reconstruction_over_flow_and_target_union_knots"),
            _text_part("phi_small_x_alternating_series_through_x8_for_abs_x_le_2^-8"),
            _text_part(
                "medium_phi=log1p(binary64(exact_Fraction_x))/"
                "binary64(exact_Fraction_x)_for_abs_exact_x_le_0.5"
            ),
            _text_part(
                "outer_segment_exposure=(N*q*dt/(M1-M0))*"
                "log(normal_binary64(exact_Fraction_M1_over_M0))"
            ),
            _text_part(
                "ratio_projection_domain=[0x1.0000000000000p-1022,"
                "0x1.fffffffffffffp+1023];otherwise_refuse"
            ),
            _text_part("inverse_tau=(needed/base_rate)*psi(z);psi(z)=expm1(z)/z"),
            _text_part("psi_small_z_series_through_z8_for_abs_z_le_2^-8"),
            _text_part(
                "analytic_inverse_seed_is_optional;"
                "conditioning_failure_uses_complete_binary64_bisection"
            ),
            _text_part("first_binary64_time_reaches_threshold_and_predecessor_does_not"),
            _text_part("decimal_oracle_acceptance=max(4_ULP,1e-15_absolute)"),
            _text_part("transcendental_binary64_arithmetic_not_exact"),
            _text_part(
                "outer_public_evaluator_replay_required;host_transaction_not_self_authenticating"
            ),
        ),
    )
    observed = (schema_digest, law_digest, algorithm_digest)
    expected = (
        AFFINE_RESIDENCE_EXPECTED_SCHEMA_DEFINITION_DIGEST,
        AFFINE_RESIDENCE_EXPECTED_LAW_DEFINITION_DIGEST,
        AFFINE_RESIDENCE_EXPECTED_ALGORITHM_DEFINITION_DIGEST,
    )
    if observed != expected:
        raise KCellTrayStateError("affine residence clock definitions differ from fixed pins")
    return observed


def _payload_identity(payload_bytes: bytes) -> str:
    return _framed_digest(HOST_DIGEST_DOMAIN + "/payload", (payload_bytes,))


def _inventory_parts(inventory: PerParticleInventory) -> tuple[bytes, ...]:
    return (
        _text_part(inventory.energy_datum_id),
        *(_float_bytes(value) for value in inventory.as_tuple()),
    )


def _same_inventory(left: PerParticleInventory, right: PerParticleInventory) -> bool:
    return (
        type(left) is PerParticleInventory
        and type(right) is PerParticleInventory
        and left.energy_datum_id == right.energy_datum_id
        and all(_same_float(a, b) for a, b in zip(left.as_tuple(), right.as_tuple(), strict=True))
    )


def _same_weight(left: PacketWeight, right: PacketWeight) -> bool:
    return (
        type(left) is PacketWeight
        and type(right) is PacketWeight
        and left.energy_datum_id == right.energy_datum_id
        and _same_float(left.representative_particles, right.representative_particles)
    )


def _packet_sort_key(packet: AcceptedTrayPacket) -> tuple[str, str, int, int, str]:
    owner = packet.owner_key
    return (
        owner.tower_id,
        owner.physical_tray_id,
        owner.vertical_layer_id.value,
        owner.rtd_stage_id.value,
        owner.cohort_id,
    )


def _stage_for_exposure(path: PacketRTDPath, exposure: float) -> int:
    # The final threshold is a ranking/exit reference, not an authorization to
    # remove material.  Crossing it therefore leaves the packet in the final
    # represented RTD stage.
    crossed_internal = sum(exposure >= value for value in path.transition_exposures[:-1])
    return crossed_internal + 1


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedDryMatterFlowSegment:
    """Piecewise-constant solid-flow declaration used by the host clock.

    ``offered_dry_matter_discharge_kg_s`` is recorded for the negative gate and
    never enters the exposure arithmetic.  This type authenticates no sensor or
    controller; its accepted authority remains a manufactured callback input.
    """

    start_time_s: float
    end_time_s: float
    accepted_flow_authority_id: str
    accepted_dry_matter_discharge_kg_s: float
    offered_dry_matter_discharge_kg_s: float
    discharge_basis: DryMatterDischargeBasis = DryMatterDischargeBasis.LOCALLY_ACCEPTED_DISCHARGE

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    accepted_flow_instrumentation_authenticated: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_binary64(
            "solid-flow interval",
            self.start_time_s,
            self.end_time_s,
            self.accepted_dry_matter_discharge_kg_s,
            self.offered_dry_matter_discharge_kg_s,
        )
        _require_nonblank("accepted flow authority", self.accepted_flow_authority_id)
        if self.start_time_s < 0.0 or self.end_time_s <= self.start_time_s:
            raise KCellTrayStateError("solid-flow interval must be strictly forward")
        if (
            self.accepted_dry_matter_discharge_kg_s < 0.0
            or self.offered_dry_matter_discharge_kg_s < 0.0
        ):
            raise KCellTrayStateError("solid-flow rates must be nonnegative")
        if type(self.discharge_basis) is not DryMatterDischargeBasis:
            raise KCellTrayStateError("solid-flow basis has the wrong type")
        if self.discharge_basis is not DryMatterDischargeBasis.LOCALLY_ACCEPTED_DISCHARGE:
            raise KCellTrayStateError("only locally accepted discharge may advance RTD exposure")

    @property
    def duration_s(self) -> float:
        return self.end_time_s - self.start_time_s


def accepted_dry_matter_flow_history_digest(
    *,
    physical_tray_id: str,
    segments: tuple[AcceptedDryMatterFlowSegment, ...],
) -> str:
    """Bind one accepted flow history independently of a tray-step request."""

    _require_nonblank("accepted flow-history tray", physical_tray_id)
    if type(segments) is not tuple or not segments:
        raise KCellTrayStateError("accepted flow-history digest requires exact segments")
    if any(type(item) is not AcceptedDryMatterFlowSegment for item in segments):
        raise KCellTrayStateError("accepted flow-history digest contains a foreign segment")
    parts: list[bytes] = [_text_part(physical_tray_id), _integer_part(len(segments))]
    for item in segments:
        AcceptedDryMatterFlowSegment.__post_init__(item)
        parts.extend(
            (
                _float_bytes(item.start_time_s),
                _float_bytes(item.end_time_s),
                _text_part(item.accepted_flow_authority_id),
                _float_bytes(item.accepted_dry_matter_discharge_kg_s),
                _float_bytes(item.offered_dry_matter_discharge_kg_s),
                _text_part(item.discharge_basis.value),
            )
        )
    return _framed_digest(HOST_DIGEST_DOMAIN + "/accepted-flow-history-v2", tuple(parts))


@dataclass(frozen=True, slots=True, kw_only=True)
class AffineDryInventoryResidenceSegment:
    """One exact piecewise-affine dry-inventory coordinate for RTD clocks."""

    segment_id: str
    start_time_s: float
    end_time_s: float
    start_dry_matter_holdup_kg: float
    end_dry_matter_holdup_kg: float

    architecture_only: ClassVar[bool] = True
    accepted_history_authenticated: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank("affine residence segment", self.segment_id)
        _require_binary64(
            "affine residence segment",
            self.start_time_s,
            self.end_time_s,
            self.start_dry_matter_holdup_kg,
            self.end_dry_matter_holdup_kg,
        )
        if self.start_time_s < 0.0 or self.end_time_s <= self.start_time_s:
            raise KCellTrayStateError("affine residence segment must advance nonnegative time")
        if _negative_zero(self.start_time_s) or _negative_zero(self.end_time_s):
            raise KCellTrayStateError("affine residence segment time rejects negative zero")
        if self.start_dry_matter_holdup_kg <= 0.0 or self.end_dry_matter_holdup_kg <= 0.0:
            raise KCellTrayStateError("affine residence inventory must remain strictly positive")

    @property
    def exact_duration_s(self) -> Fraction:
        return Fraction.from_float(self.end_time_s) - Fraction.from_float(self.start_time_s)

    @property
    def exact_slope_kg_s(self) -> Fraction:
        return (
            Fraction.from_float(self.end_dry_matter_holdup_kg)
            - Fraction.from_float(self.start_dry_matter_holdup_kg)
        ) / self.exact_duration_s

    def exact_mass_at(self, time_s: float) -> Fraction:
        _require_binary64("affine residence query time", time_s)
        if not self.start_time_s <= time_s <= self.end_time_s:
            raise KCellTrayStepError("affine residence query lies outside its segment")
        elapsed = Fraction.from_float(time_s) - Fraction.from_float(self.start_time_s)
        return Fraction.from_float(self.start_dry_matter_holdup_kg) + (
            self.exact_slope_kg_s * elapsed
        )

    @property
    def definition_digest(self) -> str:
        return _framed_digest(
            HOST_DIGEST_DOMAIN + "/affine-residence-segment-v2",
            (
                _text_part(AFFINE_RESIDENCE_EXPECTED_SCHEMA_DEFINITION_DIGEST),
                _text_part(self.segment_id),
                _float_bytes(self.start_time_s),
                _float_bytes(self.end_time_s),
                _float_bytes(self.start_dry_matter_holdup_kg),
                _float_bytes(self.end_dry_matter_holdup_kg),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AffineDryInventoryResidenceTrajectory:
    """Fully reconstructive optional V2 coordinate for one tray RTD step."""

    trajectory_id: str
    physical_tray_id: str
    source_target_definition_digest: str
    source_prior_state_digest: str
    source_flow_history_definition_digest: str
    segments: tuple[AffineDryInventoryResidenceSegment, ...]
    schema_definition_digest: str = AFFINE_RESIDENCE_EXPECTED_SCHEMA_DEFINITION_DIGEST
    law_definition_digest: str = AFFINE_RESIDENCE_EXPECTED_LAW_DEFINITION_DIGEST
    algorithm_definition_digest: str = AFFINE_RESIDENCE_EXPECTED_ALGORITHM_DEFINITION_DIGEST

    exact_prior_and_flow_pins: ClassVar[bool] = True
    fixed_prior_inventory_approximation_used: ClassVar[bool] = False
    architecture_only: ClassVar[bool] = True
    durable_persistence_implemented: ClassVar[bool] = False
    host_transaction_self_authenticating: ClassVar[bool] = False
    outer_public_evaluator_replay_required: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        validate_affine_residence_clock_definitions()
        _require_nonblank("affine residence trajectory", self.trajectory_id, self.physical_tray_id)
        for label, value in (
            ("affine residence target", self.source_target_definition_digest),
            ("affine residence prior", self.source_prior_state_digest),
            ("affine residence flow", self.source_flow_history_definition_digest),
            ("affine residence schema", self.schema_definition_digest),
            ("affine residence law", self.law_definition_digest),
            ("affine residence algorithm", self.algorithm_definition_digest),
        ):
            _require_sha256(label, value)
        if (
            self.schema_definition_digest != AFFINE_RESIDENCE_EXPECTED_SCHEMA_DEFINITION_DIGEST
            or self.law_definition_digest != AFFINE_RESIDENCE_EXPECTED_LAW_DEFINITION_DIGEST
            or self.algorithm_definition_digest
            != AFFINE_RESIDENCE_EXPECTED_ALGORITHM_DEFINITION_DIGEST
        ):
            raise KCellTrayStateError("affine residence trajectory has foreign definition pins")
        if type(self.segments) is not tuple or not self.segments:
            raise KCellTrayStateError("affine residence trajectory requires exact segments")
        if any(type(item) is not AffineDryInventoryResidenceSegment for item in self.segments):
            raise KCellTrayStateError("affine residence trajectory contains a foreign segment")
        for item in self.segments:
            AffineDryInventoryResidenceSegment.__post_init__(item)
        segment_ids = tuple(item.segment_id for item in self.segments)
        if len(set(segment_ids)) != len(segment_ids):
            raise KCellTrayStateError("affine residence segment identities must be unique")
        if tuple(self.segments) != tuple(sorted(self.segments, key=lambda item: item.start_time_s)):
            raise KCellTrayStateError("affine residence segments must be chronological")
        for left, right in zip(self.segments, self.segments[1:]):
            if not _same_float(left.end_time_s, right.start_time_s):
                raise KCellTrayStateError("affine residence segments must be exactly contiguous")
            if not _same_float(
                left.end_dry_matter_holdup_kg,
                right.start_dry_matter_holdup_kg,
            ):
                raise KCellTrayStateError(
                    "affine residence endpoint inventories must join bit for bit"
                )

    @property
    def start_time_s(self) -> float:
        return self.segments[0].start_time_s

    @property
    def end_time_s(self) -> float:
        return self.segments[-1].end_time_s

    @property
    def start_dry_matter_holdup_kg(self) -> float:
        return self.segments[0].start_dry_matter_holdup_kg

    def segment_covering(
        self,
        start_time_s: float,
        end_time_s: float,
    ) -> AffineDryInventoryResidenceSegment:
        matches = tuple(
            item
            for item in self.segments
            if item.start_time_s <= start_time_s and end_time_s <= item.end_time_s
        )
        if len(matches) != 1:
            raise KCellTrayStepError(
                "affine residence integration interval lacks one covering segment"
            )
        return matches[0]

    @property
    def definition_digest(self) -> str:
        return _framed_digest(
            HOST_DIGEST_DOMAIN + "/affine-residence-trajectory-v2",
            (
                _text_part(self.trajectory_id),
                _text_part(self.physical_tray_id),
                _text_part(self.source_target_definition_digest),
                _text_part(self.source_prior_state_digest),
                _text_part(self.source_flow_history_definition_digest),
                _text_part(self.schema_definition_digest),
                _text_part(self.law_definition_digest),
                _text_part(self.algorithm_definition_digest),
                _integer_part(len(self.segments)),
                *(_text_part(item.definition_digest) for item in self.segments),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedKCell:
    """One K cell's complementary fluid/wall conserved control volume."""

    vertical_layer_id: VerticalLayerId
    conserved_state: HostCellConservedState

    physically_qualifying: ClassVar[bool] = False
    complete_cell_physics_present: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.vertical_layer_id) is not VerticalLayerId:
            raise KCellTrayStateError("K-cell layer ID has the wrong type")
        if type(self.conserved_state) is not HostCellConservedState:
            raise KCellTrayStateError("K-cell conserved state has the wrong type")


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedTrayPacket:
    """Host-owned weighted packet with immutable payload bytes and RTD clock."""

    owner_key: ParticlePopulationKey
    entry: PacketPlanEntry
    revision: int
    time_s: float
    pressure_pa: float
    configuration_digest: str
    codec_registry_digest: str
    auditor_identity_digest: str
    rtd_path: PacketRTDPath
    rtd_schedule_digest: str
    cumulative_rtd_exposure: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.owner_key) is not ParticlePopulationKey:
            raise KCellTrayStateError("packet owner key has the wrong type")
        if type(self.entry) is not PacketPlanEntry:
            raise KCellTrayStateError("packet entry has the wrong type")
        if self.entry.packet_id != self.owner_key.cohort_id:
            raise KCellTrayStateError("packet ID and cohort identity must be one authority")
        if type(self.revision) is not int or self.revision < 0:
            raise KCellTrayStateError("packet revision must be a nonnegative exact integer")
        _require_binary64(
            "packet time, pressure, and exposure",
            self.time_s,
            self.pressure_pa,
            self.cumulative_rtd_exposure,
        )
        if self.time_s < 0.0 or self.pressure_pa <= 0.0 or self.cumulative_rtd_exposure < 0.0:
            raise KCellTrayStateError(
                "packet time/exposure must be nonnegative and pressure positive"
            )
        _require_nonblank(
            "packet configuration identities",
            self.configuration_digest,
            self.codec_registry_digest,
            self.auditor_identity_digest,
            self.rtd_schedule_digest,
        )
        if type(self.rtd_path) is not PacketRTDPath:
            raise KCellTrayStateError("packet RTD path has the wrong type")
        if self.entry.inventory.energy_datum_id != self.entry.weight.energy_datum_id:
            raise KCellTrayStateError("packet inventory and weight use different energy datums")
        if self.entry.surface is not None:
            require_single_surface_authority(self.entry.inventory, self.entry.surface)
        derive_packet_totals(self.entry.inventory, self.entry.weight)

    @property
    def packet_id(self) -> str:
        return self.entry.packet_id

    @property
    def payload_identity(self) -> str:
        return _payload_identity(self.entry.payload_bytes)

    @property
    def weighted_totals(self) -> PacketTotals:
        return derive_packet_totals(self.entry.inventory, self.entry.weight)


def _cell_digest(
    *,
    attempt_id: str,
    prior: HostCellConservedState,
    layer: VerticalLayerId,
    time_after_s: float,
    totals: CellConservedTotals,
) -> str:
    return _framed_digest(
        HOST_DIGEST_DOMAIN + "/cell",
        (
            _text_part(attempt_id),
            _text_part(prior.cell_id),
            _integer_part(prior.revision),
            _text_part(prior.state_digest),
            _integer_part(layer.value),
            _float_bytes(time_after_s),
            _text_part(prior.energy_datum_id),
            _float_bytes(totals.water_kg),
            _float_bytes(totals.hexane_kg),
            _float_bytes(totals.common_datum_energy_j),
        ),
    )


def canonical_k_cell_tray_state_digest(state: AcceptedKCellTrayState) -> str:
    parts: list[bytes] = [
        _text_part(state.topology.digest),
        _text_part(state.physical_tray_id),
        _integer_part(state.revision),
        _float_bytes(state.time_s),
        _float_bytes(state.pressure_pa),
        _text_part(state.configuration_digest),
        _text_part(state.codec_registry_digest),
        _text_part(state.auditor_identity_digest),
        _text_part(state.energy_datum_id),
        _text_part(state.rtd_schedule_digest),
    ]
    for cell in state.cells:
        conserved = cell.conserved_state
        parts.extend(
            (
                _integer_part(cell.vertical_layer_id.value),
                _text_part(conserved.cell_id),
                _integer_part(conserved.revision),
                _text_part(conserved.state_digest),
                _text_part(conserved.energy_datum_id),
                _float_bytes(conserved.totals.water_kg),
                _float_bytes(conserved.totals.hexane_kg),
                _float_bytes(conserved.totals.common_datum_energy_j),
            )
        )
    for packet in state.packets:
        owner = packet.owner_key
        parts.extend(
            (
                _text_part(owner.tower_id),
                _text_part(owner.physical_tray_id),
                _integer_part(owner.vertical_layer_id.value),
                _integer_part(owner.rtd_stage_id.value),
                _text_part(owner.cohort_id),
                _integer_part(packet.revision),
                _float_bytes(packet.time_s),
                _float_bytes(packet.pressure_pa),
                _text_part(packet.configuration_digest),
                _text_part(packet.codec_registry_digest),
                _text_part(packet.auditor_identity_digest),
                _text_part(packet.payload_identity),
                *_inventory_parts(packet.entry.inventory),
                _float_bytes(packet.entry.weight.representative_particles),
                _text_part(packet.entry.weight.energy_datum_id),
                _integer_part(packet.rtd_path.packet_index),
                *(_float_bytes(value) for value in packet.rtd_path.transition_exposures),
                _text_part(packet.rtd_schedule_digest),
                _float_bytes(packet.cumulative_rtd_exposure),
            )
        )
    return _framed_digest(HOST_DIGEST_DOMAIN + "/state", tuple(parts))


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedKCellTrayState:
    """One immutable, accepted, single-physical-tray host state."""

    topology: ParticleTopologyAuthority
    physical_tray_id: str
    revision: int
    time_s: float
    pressure_pa: float
    configuration_digest: str
    codec_registry_digest: str
    auditor_identity_digest: str
    energy_datum_id: str
    rtd_schedule_digest: str
    cells: tuple[AcceptedKCell, ...]
    packets: tuple[AcceptedTrayPacket, ...]

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    production_codec_implemented: ClassVar[bool] = False
    host_atomic_persistence_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.topology) is not ParticleTopologyAuthority:
            raise KCellTrayStateError("tray topology has the wrong type")
        _require_nonblank(
            "tray identities",
            self.physical_tray_id,
            self.configuration_digest,
            self.codec_registry_digest,
            self.auditor_identity_digest,
            self.energy_datum_id,
            self.rtd_schedule_digest,
        )
        if type(self.revision) is not int or self.revision < 0:
            raise KCellTrayStateError("tray revision must be a nonnegative exact integer")
        _require_binary64("tray time and pressure", self.time_s, self.pressure_pa)
        if self.time_s < 0.0 or self.pressure_pa <= 0.0:
            raise KCellTrayStateError("tray time must be nonnegative and pressure positive")
        tray = next(
            (
                item
                for item in self.topology.trays
                if item.physical_tray_id == self.physical_tray_id
            ),
            None,
        )
        if tray is None:
            raise KCellTrayStateError("accepted tray is absent from its topology authority")
        if type(self.cells) is not tuple or len(self.cells) != tray.vertical_layer_count:
            raise KCellTrayStateError("accepted tray must carry exactly one cell per K layer")
        if type(self.packets) is not tuple or not self.packets:
            raise KCellTrayStateError("accepted tray must carry a nonempty packet population")
        if any(type(cell) is not AcceptedKCell for cell in self.cells):
            raise KCellTrayStateError("accepted cells must have exact AcceptedKCell type")
        if any(type(packet) is not AcceptedTrayPacket for packet in self.packets):
            raise KCellTrayStateError("accepted packets must have exact AcceptedTrayPacket type")
        expected_layers = tuple(range(1, tray.vertical_layer_count + 1))
        if tuple(cell.vertical_layer_id.value for cell in self.cells) != expected_layers:
            raise KCellTrayStateError("K cells must be complete and canonically layer-ordered")
        if tuple(sorted(self.packets, key=_packet_sort_key)) != self.packets:
            raise KCellTrayStateError("packets must be in canonical owner order")
        packet_ids = tuple(packet.packet_id for packet in self.packets)
        if len(set(packet_ids)) != len(packet_ids):
            raise KCellTrayStateError("packet/cohort identities must be unique")
        cell_ids = tuple(cell.conserved_state.cell_id for cell in self.cells)
        if len(set(cell_ids)) != len(cell_ids):
            raise KCellTrayStateError("K-cell identities must be unique")
        for cell in self.cells:
            conserved = cell.conserved_state
            if conserved.revision != self.revision:
                raise KCellTrayStateError("every K cell must share the accepted host revision")
            if conserved.energy_datum_id != self.energy_datum_id:
                raise KCellTrayStateError("K cell and tray use different energy datums")
        for packet in self.packets:
            self.topology.validate_owner(packet.owner_key)
            if packet.owner_key.physical_tray_id != self.physical_tray_id:
                raise KCellTrayStateError("packet belongs to another physical tray")
            if packet.revision != self.revision or not _same_float(packet.time_s, self.time_s):
                raise KCellTrayStateError("packet revision/time differs from the accepted tray")
            if not _same_float(packet.pressure_pa, self.pressure_pa):
                raise KCellTrayStateError("packet pressure differs from the accepted tray")
            if packet.configuration_digest != self.configuration_digest:
                raise KCellTrayStateError("packet configuration differs from the accepted tray")
            if packet.codec_registry_digest != self.codec_registry_digest:
                raise KCellTrayStateError("packet codec identity differs from the accepted tray")
            if packet.auditor_identity_digest != self.auditor_identity_digest:
                raise KCellTrayStateError("packet auditor identity differs from the accepted tray")
            if packet.entry.inventory.energy_datum_id != self.energy_datum_id:
                raise KCellTrayStateError("packet and tray use different energy datums")
            if packet.rtd_schedule_digest != self.rtd_schedule_digest:
                raise KCellTrayStateError("packet and tray bind different RTD schedules")
            if len(packet.rtd_path.transition_exposures) != tray.rtd_stage_count:
                raise KCellTrayStateError("packet RTD path order differs from tray topology")
            expected_stage = _stage_for_exposure(packet.rtd_path, packet.cumulative_rtd_exposure)
            if packet.owner_key.rtd_stage_id.value != expected_stage:
                raise KCellTrayStateError("packet RTD owner disagrees with its accepted exposure")

    @property
    def state_digest(self) -> str:
        return canonical_k_cell_tray_state_digest(self)

    @property
    def dry_matter_holdup_kg(self) -> float:
        return math.fsum(packet.weighted_totals.totals[0] for packet in self.packets)


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedPacketAdvanceRequest:
    attempt_id: str
    packet: AcceptedTrayPacket
    cell_before: AcceptedKCell
    time_after_s: float

    physically_qualifying: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedPacketAdvanceProposal:
    packet_id: str
    entry_after: PacketPlanEntry
    callback_id: str
    source_identity: str
    residual_contract_passed: bool
    maximum_scaled_residual: float

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank(
            "packet advance proposal identities",
            self.packet_id,
            self.callback_id,
            self.source_identity,
        )
        if type(self.entry_after) is not PacketPlanEntry:
            raise KCellTrayStateError("packet advance candidate entry has the wrong type")
        if type(self.residual_contract_passed) is not bool:
            raise KCellTrayStateError("packet residual disposition must be an exact bool")
        _require_binary64(
            "packet maximum scaled residual",
            self.maximum_scaled_residual,
            nonnegative=True,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedPacketAdvanceRejection:
    packet_id: str
    callback_id: str
    source_identity: str
    reason: str

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank(
            "packet advance rejection",
            self.packet_id,
            self.callback_id,
            self.source_identity,
            self.reason,
        )


class ManufacturedPacketAdvanceCallback(Protocol):
    callback_id: str
    source_identity: str
    maximum_scaled_residual: float

    def propose(
        self, request: ManufacturedPacketAdvanceRequest
    ) -> ManufacturedPacketAdvanceProposal | ManufacturedPacketAdvanceRejection: ...


class WeightedPacketPayloadAuditor(Protocol):
    auditor_identity_digest: str
    source_identity: str

    def audit(self, payload_bytes: bytes) -> PerParticleInventory: ...


@dataclass(frozen=True, slots=True, kw_only=True)
class KLayerRekeyProposal:
    packet_id: str
    target_layer_id: VerticalLayerId
    authority_id: str

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank("K-rekey proposal", self.packet_id, self.authority_id)
        if type(self.target_layer_id) is not VerticalLayerId:
            raise KCellTrayStateError("K-rekey target has the wrong type")


@dataclass(frozen=True, slots=True, kw_only=True)
class KLayerRekeyPlanningRequest:
    attempt_id: str
    event_time_s: float
    topology: ParticleTopologyAuthority
    physical_tray_id: str
    packets: tuple[AcceptedTrayPacket, ...]

    physically_qualifying: ClassVar[bool] = False


class ManufacturedKRekeyPlanner(Protocol):
    planner_id: str
    source_identity: str

    def propose(self, request: KLayerRekeyPlanningRequest) -> tuple[KLayerRekeyProposal, ...]: ...


@dataclass(frozen=True, slots=True)
class NoKRekeyPlanner:
    """Explicit zero-motion manufactured planner for the first host slice."""

    planner_id: str = "manufactured-zero-k-rekey-v1"
    source_identity: str = "declared-zero-motion"

    physically_qualifying: ClassVar[bool] = False

    def propose(self, request: KLayerRekeyPlanningRequest) -> tuple[KLayerRekeyProposal, ...]:
        del request
        return ()


@dataclass(frozen=True, slots=True, kw_only=True)
class RoundedPacketCellEnergyTerms:
    """Opt-in exact ownership decomposition for one K-cell energy update.

    The physical boundary transfer remains in ``CellExternalTransfer`` and is
    independently reconciled by the caller.  These terms only tell the host to
    cancel the one already-rounded packet-to-cell transaction by its exact
    unary inverse, then book the physical wall storage owned by this K cell.
    """

    packet_energy_to_cell_j: float
    packet_energy_exact_inverse_j: float
    wall_storage_energy_j: float
    energy_datum_id: str

    physically_qualifying: ClassVar[bool] = False
    changes_physical_boundary_transfer: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_binary64(
            "rounded packet-cell energy terms",
            self.packet_energy_to_cell_j,
            self.packet_energy_exact_inverse_j,
            self.wall_storage_energy_j,
        )
        _require_nonblank("rounded packet-cell energy datum", self.energy_datum_id)
        if not _same_float(
            self.packet_energy_exact_inverse_j,
            -self.packet_energy_to_cell_j,
        ):
            raise KCellTrayStateError(
                "rounded packet-cell energy inverse is not the exact unary inverse"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class KCellTrayStepRequest:
    attempt_id: str
    prior: AcceptedKCellTrayState
    end_time_s: float
    flow_segments: tuple[AcceptedDryMatterFlowSegment, ...]
    external_cell_transfers: tuple[CellExternalTransfer, ...]
    relative_limit: float = DEFAULT_RELATIVE_LIMIT
    affine_residence_trajectory: AffineDryInventoryResidenceTrajectory | None = None
    rounded_packet_cell_energy_terms: tuple[RoundedPacketCellEnergyTerms, ...] | None = None

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank("K-cell step attempt", self.attempt_id)
        if type(self.prior) is not AcceptedKCellTrayState:
            raise KCellTrayStateError("K-cell step prior has the wrong type")
        _require_binary64("K-cell step end time", self.end_time_s)
        if self.end_time_s <= self.prior.time_s:
            raise KCellTrayStateError("K-cell step must advance time")
        if type(self.flow_segments) is not tuple or not self.flow_segments:
            raise KCellTrayStateError("K-cell step requires a nonempty flow history")
        if any(type(item) is not AcceptedDryMatterFlowSegment for item in self.flow_segments):
            raise KCellTrayStateError("flow history contains a foreign segment")
        if not _same_float(self.flow_segments[0].start_time_s, self.prior.time_s):
            raise KCellTrayStateError("flow history must start at accepted tray time")
        if not _same_float(self.flow_segments[-1].end_time_s, self.end_time_s):
            raise KCellTrayStateError("flow history must end at requested time")
        for left, right in zip(self.flow_segments, self.flow_segments[1:]):
            if not _same_float(left.end_time_s, right.start_time_s):
                raise KCellTrayStateError("flow history must be exactly contiguous")
        if type(self.external_cell_transfers) is not tuple or len(
            self.external_cell_transfers
        ) != len(self.prior.cells):
            raise KCellTrayStateError("step requires one external transfer per K cell")
        if any(type(item) is not CellExternalTransfer for item in self.external_cell_transfers):
            raise KCellTrayStateError("cell external transfer has the wrong type")
        if any(
            item.energy_datum_id != self.prior.energy_datum_id
            for item in self.external_cell_transfers
        ):
            raise KCellTrayStateError("cell external transfer uses a foreign energy datum")
        energy_terms = self.rounded_packet_cell_energy_terms
        if energy_terms is not None:
            if type(energy_terms) is not tuple or len(energy_terms) != len(self.prior.cells):
                raise KCellTrayStateError(
                    "rounded packet-cell energy decomposition requires one row per K cell"
                )
            if any(type(item) is not RoundedPacketCellEnergyTerms for item in energy_terms):
                raise KCellTrayStateError(
                    "rounded packet-cell energy decomposition has a foreign row"
                )
            for item in energy_terms:
                RoundedPacketCellEnergyTerms.__post_init__(item)
                if item.energy_datum_id != self.prior.energy_datum_id:
                    raise KCellTrayStateError(
                        "rounded packet-cell energy decomposition uses a foreign datum"
                    )
            if any(not _same_float(item.wall_storage_energy_j, 0.0) for item in energy_terms[1:]):
                raise KCellTrayStateError("rounded packet-cell wall storage is owned by K1 only")
        _require_binary64("K-cell relative limit", self.relative_limit, positive=True)
        trajectory = self.affine_residence_trajectory
        if trajectory is None:
            return
        if type(trajectory) is not AffineDryInventoryResidenceTrajectory:
            raise KCellTrayStateError("affine residence request has a foreign trajectory")
        AcceptedKCellTrayState.__post_init__(self.prior)
        for segment in self.flow_segments:
            AcceptedDryMatterFlowSegment.__post_init__(segment)
            if _negative_zero(segment.start_time_s) or _negative_zero(segment.end_time_s):
                raise KCellTrayStateError("affine residence flow time rejects negative zero")
        AffineDryInventoryResidenceTrajectory.__post_init__(trajectory)
        if trajectory.physical_tray_id != self.prior.physical_tray_id:
            raise KCellTrayStateError("affine residence trajectory names a foreign tray")
        if trajectory.source_prior_state_digest != self.prior.state_digest:
            raise KCellTrayStateError("affine residence trajectory has a foreign prior pin")
        expected_flow_digest = accepted_dry_matter_flow_history_digest(
            physical_tray_id=self.prior.physical_tray_id,
            segments=self.flow_segments,
        )
        if trajectory.source_flow_history_definition_digest != expected_flow_digest:
            raise KCellTrayStateError("affine residence trajectory has a foreign flow pin")
        if not _same_float(trajectory.start_time_s, self.prior.time_s):
            raise KCellTrayStateError("affine residence trajectory must start at prior time")
        if not _same_float(trajectory.end_time_s, self.end_time_s):
            raise KCellTrayStateError("affine residence trajectory must end at requested time")
        if not _same_float(
            trajectory.start_dry_matter_holdup_kg,
            self.prior.dry_matter_holdup_kg,
        ):
            raise KCellTrayStateError(
                "affine residence trajectory must reconstruct exact prior dry inventory"
            )


class WholePacketRekeyMechanism(enum.Enum):
    RTD_STAGE = "rtd_stage_rekey"
    K_LAYER = "k_layer_rekey"


@dataclass(frozen=True, slots=True, kw_only=True)
class WholePacketRekeyRecord:
    packet_id: str
    mechanism: WholePacketRekeyMechanism
    event_time_s: float
    source_owner: ParticlePopulationKey
    target_owner: ParticlePopulationKey
    payload_identity: str
    weighted_totals: tuple[float, ...]
    authority_id: str

    physically_qualifying: ClassVar[bool] = False
    payload_bytes_unchanged_by_motion: ClassVar[bool] = True


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedResidenceClockLedger:
    dry_matter_holdup_kg: float
    rtd_stage_count: int
    accepted_dry_matter_throughput_kg: float
    offered_dry_matter_throughput_kg: float
    exposure_increment: float

    physically_qualifying: ClassVar[bool] = False
    offered_flow_changes_clock: ClassVar[bool] = False
    affine_trajectory_reconstructed_by_host: ClassVar[bool] = True
    supplied_exposure_or_crossing_trusted: ClassVar[bool] = False
    transcendental_binary64_arithmetic_exact: ClassVar[bool] = False
    normal_binary64_mass_ratio_projection_required: ClassVar[bool] = True
    subnormal_or_nonfinite_mass_ratio_supported: ClassVar[bool] = False

    @property
    def zero_accepted_flow_froze_clock(self) -> bool:
        return self.accepted_dry_matter_throughput_kg == 0.0 and self.exposure_increment == 0.0


def _relative_residual(residual: float, *terms: float) -> float:
    scale = max(1.0, *(abs(value) for value in terms))
    return abs(residual) / scale


@dataclass(frozen=True, slots=True, kw_only=True)
class WeightedTrayLedger:
    packet_before: tuple[float, ...]
    packet_after_advance: tuple[float, ...]
    packet_after_motion: tuple[float, ...]
    cell_before: tuple[float, float, float]
    cell_after: tuple[float, float, float]
    external_to_cells: tuple[float, float, float]
    dry_matter_residual_kg: float
    residual_oil_label_residual_kg: float
    water_residual_kg: float
    hexane_residual_kg: float
    common_datum_energy_residual_j: float
    motion_residuals: tuple[float, ...]
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return self.maximum_relative_residual <= self.relative_limit and all(
            value == 0.0 for value in self.motion_residuals
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class KCellTrayStepRejection:
    attempt_id: str
    reason: str
    rollback_state: AcceptedKCellTrayState

    physically_qualifying: ClassVar[bool] = False


_VALIDATED_STEP_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedKCellTrayStep:
    request: KCellTrayStepRequest
    request_digest: str
    candidate: AcceptedKCellTrayState
    candidate_digest: str
    clock_ledger: AcceptedResidenceClockLedger
    tray_ledger: WeightedTrayLedger
    rekeys: tuple[WholePacketRekeyRecord, ...]
    packet_callback_id: str
    packet_callback_source_identity: str
    k_planner_id: str
    k_planner_source_identity: str
    auditor_identity_digest: str
    auditor_source_identity: str
    _seal: object = field(repr=False, compare=False)

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    self_authenticating: ClassVar[bool] = False
    outer_public_evaluator_replay_required: ClassVar[bool] = True

    def __post_init__(self) -> None:
        if self._seal is not _VALIDATED_STEP_SEAL:
            raise TypeError("validated K-cell tray steps are issued only by the host evaluator")


def _request_digest(
    request: KCellTrayStepRequest,
    *,
    packet_callback_id: str,
    packet_callback_source_identity: str,
    planner_id: str,
    planner_source_identity: str,
    auditor_identity_digest: str,
    auditor_source_identity: str,
) -> str:
    parts: list[bytes] = [
        _text_part(request.attempt_id),
        _text_part(request.prior.state_digest),
        _float_bytes(request.end_time_s),
        _float_bytes(request.relative_limit),
        _text_part(packet_callback_id),
        _text_part(packet_callback_source_identity),
        _text_part(planner_id),
        _text_part(planner_source_identity),
        _text_part(auditor_identity_digest),
        _text_part(auditor_source_identity),
    ]
    for segment in request.flow_segments:
        parts.extend(
            (
                _float_bytes(segment.start_time_s),
                _float_bytes(segment.end_time_s),
                _text_part(segment.accepted_flow_authority_id),
                _float_bytes(segment.accepted_dry_matter_discharge_kg_s),
                _float_bytes(segment.offered_dry_matter_discharge_kg_s),
                _text_part(segment.discharge_basis.value),
            )
        )
    for transfer in request.external_cell_transfers:
        parts.extend(
            (
                _float_bytes(transfer.water_to_cell_kg),
                _float_bytes(transfer.hexane_to_cell_kg),
                _float_bytes(transfer.common_datum_energy_to_cell_j),
                _text_part(transfer.energy_datum_id),
            )
        )
    if request.rounded_packet_cell_energy_terms is not None:
        parts.append(_text_part("rounded-packet-cell-energy-terms-v1"))
        for item in request.rounded_packet_cell_energy_terms:
            parts.extend(
                (
                    _float_bytes(item.packet_energy_to_cell_j),
                    _float_bytes(item.packet_energy_exact_inverse_j),
                    _float_bytes(item.wall_storage_energy_j),
                    _text_part(item.energy_datum_id),
                )
            )
    if request.affine_residence_trajectory is not None:
        parts.append(_text_part(request.affine_residence_trajectory.definition_digest))
    return _framed_digest(HOST_DIGEST_DOMAIN + "/request", tuple(parts))


def _audit_entry(
    packet: AcceptedTrayPacket,
    auditor: WeightedPacketPayloadAuditor,
    *,
    label: str,
) -> None:
    observed = auditor.audit(packet.entry.payload_bytes)
    if type(observed) is not PerParticleInventory:
        raise KCellTrayStepError(f"{label} payload auditor returned a foreign inventory type")
    if not _same_inventory(observed, packet.entry.inventory):
        raise KCellTrayStepError(f"{label} payload bytes disagree with declared inventory")
    if packet.entry.surface is not None:
        require_single_surface_authority(packet.entry.inventory, packet.entry.surface)


def _sum_packet_totals(packets: tuple[AcceptedTrayPacket, ...]) -> tuple[float, ...]:
    rows = tuple(packet.weighted_totals.totals for packet in packets)
    return tuple(math.fsum(row[index] for row in rows) for index in range(7))


def _sum_cell_totals(cells: tuple[AcceptedKCell, ...]) -> tuple[float, float, float]:
    return (
        math.fsum(cell.conserved_state.totals.water_kg for cell in cells),
        math.fsum(cell.conserved_state.totals.hexane_kg for cell in cells),
        math.fsum(cell.conserved_state.totals.common_datum_energy_j for cell in cells),
    )


def _fraction_to_binary64(name: str, value: Fraction) -> float:
    try:
        result = float(value)
    except OverflowError as error:
        raise KCellTrayStepError(f"{name} overflows binary64") from error
    if not math.isfinite(result):
        raise KCellTrayStepError(f"{name} is not finite binary64")
    return result


def _phi_small_relative_change(value: float) -> float:
    """Return the certified local series for log1p(x)/x."""

    if value == 0.0:
        return 1.0
    if not math.isfinite(value) or abs(value) > AFFINE_RESIDENCE_SMALL_ARGUMENT_LIMIT:
        raise KCellTrayStepError("affine residence small-change series is out of domain")
    # Sum through x**8.  At the selected limit the absolute omitted tail is
    # below |x|**9 / (10 * (1 - |x|)), well beneath binary64 roundoff.
    return 1.0 + value * (
        -1.0 / 2.0
        + value
        * (
            1.0 / 3.0
            + value
            * (
                -1.0 / 4.0
                + value
                * (
                    1.0 / 5.0
                    + value
                    * (
                        -1.0 / 6.0
                        + value * (1.0 / 7.0 + value * (-1.0 / 8.0 + value * (1.0 / 9.0)))
                    )
                )
            )
        )
    )


def _normal_binary64_mass_ratio(value: Fraction) -> float:
    """Project an exact mass ratio only inside the pinned normal-float domain."""

    if value <= 0:
        raise KCellTrayStepError("affine residence mass ratio must be strictly positive")
    try:
        rounded_ratio = float(value)
    except OverflowError as error:
        raise KCellTrayStepError(
            "affine residence mass ratio exceeds the pinned normal binary64 domain"
        ) from error
    if (
        not math.isfinite(rounded_ratio)
        or rounded_ratio < AFFINE_RESIDENCE_MIN_NORMAL_MASS_RATIO
        or rounded_ratio > AFFINE_RESIDENCE_MAX_FINITE_MASS_RATIO
    ):
        raise KCellTrayStepError(
            "affine residence mass ratio is outside the pinned normal binary64 domain"
        )
    return rounded_ratio


def _psi_expm1_over_z(value: float) -> float:
    """Return expm1(z)/z with its exact zero-slope limit."""

    if value == 0.0:
        return 1.0
    if not math.isfinite(value):
        raise KCellTrayStepError("affine residence inverse argument is non-finite")
    if abs(value) <= AFFINE_RESIDENCE_SMALL_ARGUMENT_LIMIT:
        return 1.0 + value * (
            1.0 / 2.0
            + value
            * (
                1.0 / 6.0
                + value
                * (
                    1.0 / 24.0
                    + value
                    * (
                        1.0 / 120.0
                        + value
                        * (
                            1.0 / 720.0
                            + value
                            * (1.0 / 5040.0 + value * (1.0 / 40320.0 + value * (1.0 / 362880.0)))
                        )
                    )
                )
            )
        )
    return math.expm1(value) / value


def _affine_residence_exact_interval_exposure(
    *,
    stage_count: int,
    accepted_flow: Fraction,
    duration: Fraction,
    mass_before: Fraction,
    mass_after: Fraction,
) -> float:
    """Evaluate one exact-data affine interval inside the pinned float domain."""

    if type(stage_count) is not int or stage_count <= 0:
        raise KCellTrayStepError("affine residence stage count must be a positive integer")
    if accepted_flow < 0 or duration <= 0 or mass_before <= 0 or mass_after <= 0:
        raise KCellTrayStepError("affine residence exact interval is inadmissible")
    if accepted_flow == 0:
        return 0.0
    exact_relative_change = (mass_after - mass_before) / mass_before
    absolute_change = abs(exact_relative_change)
    if absolute_change <= Fraction.from_float(AFFINE_RESIDENCE_LOG1P_RELATIVE_LIMIT):
        relative_change = _fraction_to_binary64(
            "affine residence relative change",
            exact_relative_change,
        )
        if exact_relative_change != 0 and relative_change == 0.0:
            raise KCellTrayStepError("affine residence relative change underflowed binary64")
        base = _fraction_to_binary64(
            "affine residence constant-limit exposure",
            Fraction(stage_count) * accepted_flow * duration / mass_before,
        )
        if base == 0.0:
            raise KCellTrayStepError("positive affine residence exposure underflowed binary64")
        if absolute_change <= Fraction.from_float(AFFINE_RESIDENCE_SMALL_ARGUMENT_LIMIT):
            phi = _phi_small_relative_change(relative_change)
        else:
            phi = math.log1p(relative_change) / relative_change
        exposure = base * phi
    else:
        mass_ratio = _normal_binary64_mass_ratio(mass_after / mass_before)
        coefficient = _fraction_to_binary64(
            "affine residence logarithmic coefficient",
            Fraction(stage_count) * accepted_flow * duration / (mass_after - mass_before),
        )
        exposure = coefficient * math.log(mass_ratio)
    if not math.isfinite(exposure) or exposure <= 0.0:
        raise KCellTrayStepError("affine residence exposure is not positive finite binary64")
    return exposure


def _affine_residence_partition(
    request: KCellTrayStepRequest,
    *,
    start_time_s: float,
    end_time_s: float,
) -> tuple[float, ...]:
    trajectory = request.affine_residence_trajectory
    if type(trajectory) is not AffineDryInventoryResidenceTrajectory:
        raise KCellTrayStepError("affine residence partition requires an exact trajectory")
    _require_binary64("affine residence integration time", start_time_s, end_time_s)
    if (
        start_time_s < request.prior.time_s
        or end_time_s > request.end_time_s
        or end_time_s < start_time_s
    ):
        raise KCellTrayStepError("affine residence integration lies outside the request")
    knots = {start_time_s, end_time_s}
    for segment in request.flow_segments:
        if start_time_s < segment.start_time_s < end_time_s:
            knots.add(segment.start_time_s)
        if start_time_s < segment.end_time_s < end_time_s:
            knots.add(segment.end_time_s)
    for segment in trajectory.segments:
        if start_time_s < segment.start_time_s < end_time_s:
            knots.add(segment.start_time_s)
        if start_time_s < segment.end_time_s < end_time_s:
            knots.add(segment.end_time_s)
    return tuple(sorted(knots))


def _flow_segment_covering(
    request: KCellTrayStepRequest,
    *,
    start_time_s: float,
    end_time_s: float,
) -> AcceptedDryMatterFlowSegment:
    matches = tuple(
        segment
        for segment in request.flow_segments
        if segment.start_time_s <= start_time_s and end_time_s <= segment.end_time_s
    )
    if len(matches) != 1:
        raise KCellTrayStepError(
            "affine residence integration interval lacks one covering flow segment"
        )
    return matches[0]


def _affine_residence_interval_exposure(
    request: KCellTrayStepRequest,
    *,
    start_time_s: float,
    end_time_s: float,
    stage_count: int,
) -> float:
    trajectory = request.affine_residence_trajectory
    if type(trajectory) is not AffineDryInventoryResidenceTrajectory:
        raise KCellTrayStepError("affine residence exposure requires an exact trajectory")
    target_segment = trajectory.segment_covering(start_time_s, end_time_s)
    flow_segment = _flow_segment_covering(
        request,
        start_time_s=start_time_s,
        end_time_s=end_time_s,
    )
    accepted_flow = Fraction.from_float(flow_segment.accepted_dry_matter_discharge_kg_s)
    duration = Fraction.from_float(end_time_s) - Fraction.from_float(start_time_s)
    mass_before = target_segment.exact_mass_at(start_time_s)
    mass_after = target_segment.exact_mass_at(end_time_s)
    return _affine_residence_exact_interval_exposure(
        stage_count=stage_count,
        accepted_flow=accepted_flow,
        duration=duration,
        mass_before=mass_before,
        mass_after=mass_after,
    )


def _affine_residence_exposure_between(
    request: KCellTrayStepRequest,
    *,
    start_time_s: float,
    end_time_s: float,
    stage_count: int,
) -> float:
    if type(stage_count) is not int or stage_count <= 0:
        raise KCellTrayStepError("affine residence stage count must be a positive integer")
    knots = _affine_residence_partition(
        request,
        start_time_s=start_time_s,
        end_time_s=end_time_s,
    )
    if len(knots) == 1:
        return 0.0
    contributions = tuple(
        _affine_residence_interval_exposure(
            request,
            start_time_s=left,
            end_time_s=right,
            stage_count=stage_count,
        )
        for left, right in zip(knots, knots[1:])
    )
    exposure = math.fsum(contributions)
    if not math.isfinite(exposure) or exposure < 0.0:
        raise KCellTrayStepError("affine residence exposure sum is inadmissible")
    return exposure


def _optional_affine_inverse_elapsed(
    *,
    needed_exposure: float,
    stage_count: int,
    accepted_flow: Fraction,
    mass_before: Fraction,
    exact_slope: Fraction,
) -> float | None:
    """Return a conditioned analytic seed only; never decide admissibility."""

    if (
        type(needed_exposure) is not float
        or not math.isfinite(needed_exposure)
        or needed_exposure <= 0.0
        or type(stage_count) is not int
        or stage_count <= 0
        or accepted_flow <= 0
        or mass_before <= 0
    ):
        return None
    try:
        throughput_factor = Fraction(stage_count) * accepted_flow
        exact_base_rate = throughput_factor / mass_before
        base_rate = _fraction_to_binary64(
            "affine residence inverse base rate",
            exact_base_rate,
        )
        if base_rate <= 0.0 or (exact_base_rate != 0 and base_rate == 0.0):
            return None
        exact_slope_ratio = exact_slope / throughput_factor
        slope_ratio = _fraction_to_binary64(
            "affine residence inverse slope ratio",
            exact_slope_ratio,
        )
        if exact_slope_ratio != 0 and slope_ratio == 0.0:
            return None
        inverse_argument = needed_exposure * slope_ratio
        if not math.isfinite(inverse_argument):
            return None
        zero_slope_elapsed = needed_exposure / base_rate
        if not math.isfinite(zero_slope_elapsed) or zero_slope_elapsed <= 0.0:
            return None
        elapsed = zero_slope_elapsed * _psi_expm1_over_z(inverse_argument)
    except ArithmeticError, KCellTrayStepError, ValueError:
        return None
    if not math.isfinite(elapsed) or elapsed <= 0.0:
        return None
    return elapsed


def _stable_affine_crossing_seed(
    request: KCellTrayStepRequest,
    *,
    stage_count: int,
    exposure_before: float,
    threshold: float,
) -> float | None:
    trajectory = request.affine_residence_trajectory
    if type(trajectory) is not AffineDryInventoryResidenceTrajectory:
        raise KCellTrayStepError("affine residence inverse requires an exact trajectory")
    knots = _affine_residence_partition(
        request,
        start_time_s=request.prior.time_s,
        end_time_s=request.end_time_s,
    )
    intervals = tuple(zip(knots, knots[1:]))
    contributions = tuple(
        _affine_residence_interval_exposure(
            request,
            start_time_s=left,
            end_time_s=right,
            stage_count=stage_count,
        )
        for left, right in intervals
    )
    for index, (left, right) in enumerate(intervals):
        exposure = math.fsum((exposure_before, math.fsum(contributions[:index])))
        exposure_after = math.fsum((exposure_before, math.fsum(contributions[: index + 1])))
        if not exposure < threshold <= exposure_after:
            continue
        flow_segment = _flow_segment_covering(
            request,
            start_time_s=left,
            end_time_s=right,
        )
        accepted_flow = Fraction.from_float(flow_segment.accepted_dry_matter_discharge_kg_s)
        if accepted_flow == 0:  # pragma: no cover - excluded by the crossing inequality
            raise KCellTrayStepError("positive affine threshold crossing has zero flow")
        target_segment = trajectory.segment_covering(left, right)
        mass_before = target_segment.exact_mass_at(left)
        exact_slope = target_segment.exact_slope_kg_s
        needed = math.fsum((threshold, -exposure))
        elapsed = _optional_affine_inverse_elapsed(
            needed_exposure=needed,
            stage_count=stage_count,
            accepted_flow=accepted_flow,
            mass_before=mass_before,
            exact_slope=exact_slope,
        )
        if elapsed is None:
            return None
        seed = left + elapsed
        if not math.isfinite(seed) or seed <= left or seed > right:
            return None
        # The seed only accelerates the certified binary64 search below.  It is
        # never accepted as the event time without the forward/predecessor test.
        return seed
    # Reachability is certified independently from the same global fsum path.
    # A missing numerical seed therefore falls back to the complete binary64
    # key search instead of creating a false rejection.
    return None


def _nonnegative_binary64_key(value: float) -> int:
    _require_binary64("affine residence binary64 time", value, nonnegative=True)
    if value == 0.0:
        return 0
    return struct.unpack(">Q", _float_bytes(value))[0]


def _binary64_from_nonnegative_key(value: int) -> float:
    return struct.unpack(">d", struct.pack(">Q", value))[0]


def _first_binary64_affine_crossing(
    request: KCellTrayStepRequest,
    *,
    stage_count: int,
    exposure_before: float,
    threshold: float,
) -> float:
    start = request.prior.time_s
    end = request.end_time_s

    def reached(time_s: float) -> bool:
        increment = _affine_residence_exposure_between(
            request,
            start_time_s=start,
            end_time_s=time_s,
            stage_count=stage_count,
        )
        return math.fsum((exposure_before, increment)) >= threshold

    if reached(start):
        raise KCellTrayStepError("affine residence threshold was already reached at step start")
    if not reached(end):
        raise KCellTrayStepError("affine residence threshold is not reached by step end")
    seed = _stable_affine_crossing_seed(
        request,
        stage_count=stage_count,
        exposure_before=exposure_before,
        threshold=threshold,
    )
    lower = _nonnegative_binary64_key(start)
    upper = _nonnegative_binary64_key(end)
    if seed is not None:
        seed_key = _nonnegative_binary64_key(seed)
        if lower < seed_key < upper:
            if reached(seed):
                upper = seed_key
            else:
                lower = seed_key
    while upper - lower > 1:
        midpoint = (lower + upper) // 2
        if reached(_binary64_from_nonnegative_key(midpoint)):
            upper = midpoint
        else:
            lower = midpoint
    candidate = _binary64_from_nonnegative_key(upper)
    predecessor = _binary64_from_nonnegative_key(upper - 1)
    if not reached(candidate) or (
        predecessor >= start and reached(predecessor)
    ):  # pragma: no cover - final fail-closed certificate
        raise KCellTrayStepError("affine residence first-binary64 crossing certificate failed")
    return candidate


def _clock_ledger(request: KCellTrayStepRequest) -> AcceptedResidenceClockLedger:
    if request.affine_residence_trajectory is None:
        holdup = request.prior.dry_matter_holdup_kg
        if holdup <= 0.0:
            raise KCellTrayStepError("RTD exposure requires positive weighted dry-matter holdup")
        tray = next(
            item
            for item in request.prior.topology.trays
            if item.physical_tray_id == request.prior.physical_tray_id
        )
        accepted = math.fsum(
            segment.accepted_dry_matter_discharge_kg_s * segment.duration_s
            for segment in request.flow_segments
        )
        offered = math.fsum(
            segment.offered_dry_matter_discharge_kg_s * segment.duration_s
            for segment in request.flow_segments
        )
        exposure = math.fsum(
            tray.rtd_stage_count
            * segment.accepted_dry_matter_discharge_kg_s
            / holdup
            * segment.duration_s
            for segment in request.flow_segments
        )
        if accepted == 0.0 and exposure != 0.0:
            raise KCellTrayStepError("zero accepted discharge did not freeze RTD exposure")
        _require_binary64(
            "accepted residence clock",
            holdup,
            accepted,
            offered,
            exposure,
            nonnegative=True,
        )
        return AcceptedResidenceClockLedger(
            dry_matter_holdup_kg=holdup,
            rtd_stage_count=tray.rtd_stage_count,
            accepted_dry_matter_throughput_kg=accepted,
            offered_dry_matter_throughput_kg=offered,
            exposure_increment=exposure,
        )

    holdup = request.prior.dry_matter_holdup_kg
    if holdup <= 0.0:
        raise KCellTrayStepError("RTD exposure requires positive weighted dry-matter holdup")
    tray = next(
        item
        for item in request.prior.topology.trays
        if item.physical_tray_id == request.prior.physical_tray_id
    )
    accepted = math.fsum(
        segment.accepted_dry_matter_discharge_kg_s * segment.duration_s
        for segment in request.flow_segments
    )
    offered = math.fsum(
        segment.offered_dry_matter_discharge_kg_s * segment.duration_s
        for segment in request.flow_segments
    )
    exposure = _affine_residence_exposure_between(
        request,
        start_time_s=request.prior.time_s,
        end_time_s=request.end_time_s,
        stage_count=tray.rtd_stage_count,
    )
    if accepted == 0.0 and exposure != 0.0:
        raise KCellTrayStepError("zero accepted discharge did not freeze RTD exposure")
    _require_binary64(
        "accepted affine residence clock",
        holdup,
        accepted,
        offered,
        exposure,
        nonnegative=True,
    )
    return AcceptedResidenceClockLedger(
        dry_matter_holdup_kg=holdup,
        rtd_stage_count=tray.rtd_stage_count,
        accepted_dry_matter_throughput_kg=accepted,
        offered_dry_matter_throughput_kg=offered,
        exposure_increment=exposure,
    )


def _crossing_time(
    *,
    request: KCellTrayStepRequest,
    holdup_kg: float,
    stage_count: int,
    exposure_before: float,
    threshold: float,
) -> float:
    if request.affine_residence_trajectory is None:
        exposure = exposure_before
        for segment in request.flow_segments:
            rate = stage_count * segment.accepted_dry_matter_discharge_kg_s / holdup_kg
            contribution = rate * segment.duration_s
            exposure_after = math.fsum((exposure, contribution))
            if exposure < threshold <= exposure_after:
                if rate <= 0.0:  # pragma: no cover - guarded by the inequality
                    raise KCellTrayStepError("positive RTD threshold crossing has zero rate")
                return segment.start_time_s + (threshold - exposure) / rate
            exposure = exposure_after
        raise KCellTrayStepError("RTD threshold crossing could not be located in accepted history")

    return _first_binary64_affine_crossing(
        request,
        stage_count=stage_count,
        exposure_before=exposure_before,
        threshold=threshold,
    )


def _advance_packets(
    request: KCellTrayStepRequest,
    *,
    callback: ManufacturedPacketAdvanceCallback,
    auditor: WeightedPacketPayloadAuditor,
) -> tuple[AcceptedTrayPacket, ...] | KCellTrayStepRejection:
    by_layer = {cell.vertical_layer_id.value: cell for cell in request.prior.cells}
    candidates: list[AcceptedTrayPacket] = []
    for packet in request.prior.packets:
        _audit_entry(packet, auditor, label="accepted-prior")
        outcome = callback.propose(
            ManufacturedPacketAdvanceRequest(
                attempt_id=request.attempt_id,
                packet=packet,
                cell_before=by_layer[packet.owner_key.vertical_layer_id.value],
                time_after_s=request.end_time_s,
            )
        )
        _audit_entry(packet, auditor, label="accepted-prior-after-callback")
        if type(outcome) is ManufacturedPacketAdvanceRejection:
            if (
                outcome.packet_id != packet.packet_id
                or outcome.callback_id != callback.callback_id
                or outcome.source_identity != callback.source_identity
            ):
                raise KCellTrayStepError("packet rejection carries foreign identities")
            return KCellTrayStepRejection(
                attempt_id=request.attempt_id,
                reason=outcome.reason,
                rollback_state=request.prior,
            )
        if type(outcome) is not ManufacturedPacketAdvanceProposal:
            raise KCellTrayStepError("packet callback returned a foreign outcome type")
        if (
            outcome.packet_id != packet.packet_id
            or outcome.entry_after.packet_id != packet.packet_id
        ):
            raise KCellTrayStepError("packet callback changed packet identity")
        if (
            outcome.callback_id != callback.callback_id
            or outcome.source_identity != callback.source_identity
        ):
            raise KCellTrayStepError("packet callback result violates its source pin")
        if not outcome.residual_contract_passed:
            raise KCellTrayStepError("packet callback residual contract did not pass")
        if outcome.maximum_scaled_residual > callback.maximum_scaled_residual:
            raise KCellTrayStepError("packet callback exceeded its residual pin")
        if not _same_weight(outcome.entry_after.weight, packet.entry.weight):
            raise KCellTrayStepError("local packet advance changed host-carried packet weight")
        before_inventory = packet.entry.inventory
        after_inventory = outcome.entry_after.inventory
        if not _same_float(before_inventory.dry_matter_kg, after_inventory.dry_matter_kg):
            raise KCellTrayStepError("local packet advance changed dry-matter identity")
        if not _same_float(
            before_inventory.residual_oil_label_kg,
            after_inventory.residual_oil_label_kg,
        ):
            raise KCellTrayStepError("local packet advance changed residual-oil identity")
        provisional = replace(
            packet,
            entry=outcome.entry_after,
            revision=packet.revision + 1,
            time_s=request.end_time_s,
        )
        _audit_entry(provisional, auditor, label="packet-candidate")
        candidates.append(provisional)
    return tuple(candidates)


def _apply_rtd_rekeys(
    request: KCellTrayStepRequest,
    packets: tuple[AcceptedTrayPacket, ...],
    clock: AcceptedResidenceClockLedger,
) -> tuple[tuple[AcceptedTrayPacket, ...], tuple[WholePacketRekeyRecord, ...]]:
    moved: list[AcceptedTrayPacket] = []
    records: list[WholePacketRekeyRecord] = []
    for packet in packets:
        exposure_before = next(
            item.cumulative_rtd_exposure
            for item in request.prior.packets
            if item.packet_id == packet.packet_id
        )
        exposure_after = math.fsum((exposure_before, clock.exposure_increment))
        _require_binary64("packet accepted RTD exposure", exposure_after, nonnegative=True)
        owner = packet.owner_key
        current_owner = owner
        for stage_index, threshold in enumerate(packet.rtd_path.transition_exposures[:-1], start=1):
            if not exposure_before < threshold <= exposure_after:
                continue
            target_owner = replace(current_owner, rtd_stage_id=RTDStageId(stage_index + 1))
            records.append(
                WholePacketRekeyRecord(
                    packet_id=packet.packet_id,
                    mechanism=WholePacketRekeyMechanism.RTD_STAGE,
                    event_time_s=_crossing_time(
                        request=request,
                        holdup_kg=clock.dry_matter_holdup_kg,
                        stage_count=clock.rtd_stage_count,
                        exposure_before=exposure_before,
                        threshold=threshold,
                    ),
                    source_owner=current_owner,
                    target_owner=target_owner,
                    payload_identity=packet.payload_identity,
                    weighted_totals=packet.weighted_totals.totals,
                    authority_id=packet.rtd_schedule_digest,
                )
            )
            current_owner = target_owner
        expected_stage = _stage_for_exposure(packet.rtd_path, exposure_after)
        if current_owner.rtd_stage_id.value != expected_stage:
            current_owner = replace(current_owner, rtd_stage_id=RTDStageId(expected_stage))
        moved.append(
            replace(
                packet,
                owner_key=current_owner,
                cumulative_rtd_exposure=exposure_after,
            )
        )
    return tuple(moved), tuple(records)


def _apply_k_rekeys(
    request: KCellTrayStepRequest,
    packets: tuple[AcceptedTrayPacket, ...],
    *,
    planner: ManufacturedKRekeyPlanner,
) -> tuple[tuple[AcceptedTrayPacket, ...], tuple[WholePacketRekeyRecord, ...]]:
    planning_request = KLayerRekeyPlanningRequest(
        attempt_id=request.attempt_id,
        event_time_s=request.end_time_s,
        topology=request.prior.topology,
        physical_tray_id=request.prior.physical_tray_id,
        packets=packets,
    )
    proposals = planner.propose(planning_request)
    if type(proposals) is not tuple or any(
        type(item) is not KLayerRekeyProposal for item in proposals
    ):
        raise KCellTrayStepError("K planner must return an exact tuple of exact proposals")
    proposal_ids = tuple(item.packet_id for item in proposals)
    if proposal_ids != tuple(sorted(proposal_ids)):
        raise KCellTrayStepError("K planner proposals must be canonically packet-ID ordered")
    if len(set(proposal_ids)) != len(proposal_ids):
        raise KCellTrayStepError("K planner proposed more than one move for one packet")
    by_id = {packet.packet_id: packet for packet in packets}
    for proposal in proposals:
        if proposal.authority_id != planner.planner_id:
            raise KCellTrayStepError("K-rekey proposal violates planner identity")
        if proposal.packet_id not in by_id:
            raise KCellTrayStepError("K planner named a foreign packet")
        packet = by_id[proposal.packet_id]
        source_layer = packet.owner_key.vertical_layer_id.value
        target_layer = proposal.target_layer_id.value
        if abs(target_layer - source_layer) != 1:
            raise KCellTrayStepError("K rekeys must cross exactly one adjacent layer face")
        target_owner = replace(packet.owner_key, vertical_layer_id=proposal.target_layer_id)
        request.prior.topology.validate_owner(target_owner)
        candidate = replace(packet, owner_key=target_owner)
        if candidate.payload_identity != packet.payload_identity:
            raise KCellTrayStepError("K rekey changed complete payload identity")
        if candidate.weighted_totals.totals != packet.weighted_totals.totals:
            raise KCellTrayStepError("K rekey changed weighted packet inventory")
        by_id[proposal.packet_id] = candidate
    ordered = tuple(sorted(by_id.values(), key=_packet_sort_key))
    records = tuple(
        WholePacketRekeyRecord(
            packet_id=proposal.packet_id,
            mechanism=WholePacketRekeyMechanism.K_LAYER,
            event_time_s=request.end_time_s,
            source_owner=next(
                packet.owner_key for packet in packets if packet.packet_id == proposal.packet_id
            ),
            target_owner=next(
                packet.owner_key for packet in ordered if packet.packet_id == proposal.packet_id
            ),
            payload_identity=next(
                packet.payload_identity
                for packet in packets
                if packet.packet_id == proposal.packet_id
            ),
            weighted_totals=next(
                packet.weighted_totals.totals
                for packet in packets
                if packet.packet_id == proposal.packet_id
            ),
            authority_id=proposal.authority_id,
        )
        for proposal in proposals
    )
    return ordered, records


def _candidate_cells(
    request: KCellTrayStepRequest,
    advanced_packets: tuple[AcceptedTrayPacket, ...],
) -> tuple[AcceptedKCell, ...]:
    prior_by_id = {packet.packet_id: packet for packet in request.prior.packets}
    candidates: list[AcceptedKCell] = []
    for index, (cell, external) in enumerate(
        zip(request.prior.cells, request.external_cell_transfers, strict=True)
    ):
        layer = cell.vertical_layer_id.value
        prior_rows = tuple(
            prior_by_id[packet.packet_id].weighted_totals.totals
            for packet in advanced_packets
            if prior_by_id[packet.packet_id].owner_key.vertical_layer_id.value == layer
        )
        after_rows = tuple(
            packet.weighted_totals.totals
            for packet in advanced_packets
            if prior_by_id[packet.packet_id].owner_key.vertical_layer_id.value == layer
        )
        before_water = math.fsum(row[4] + row[5] for row in prior_rows)
        after_water = math.fsum(row[4] + row[5] for row in after_rows)
        before_hexane = math.fsum(row[2] + row[3] for row in prior_rows)
        after_hexane = math.fsum(row[2] + row[3] for row in after_rows)
        before_energy = math.fsum(row[6] for row in prior_rows)
        after_energy = math.fsum(row[6] for row in after_rows)
        packet_energy_to_cell = math.fsum((before_energy, -after_energy))
        old = cell.conserved_state
        packet_water_to_cell = math.fsum((before_water, -after_water))
        packet_hexane_to_cell = math.fsum((before_hexane, -after_hexane))
        # Reuse the rounded packet source so its exact unary inverse cancels.
        proposed_water = math.fsum(
            (old.totals.water_kg, packet_water_to_cell, external.water_to_cell_kg)
        )
        proposed_hexane = math.fsum(
            (
                old.totals.hexane_kg,
                packet_hexane_to_cell,
                external.hexane_to_cell_kg,
            )
        )
        details: list[str] = []
        for (
            component,
            prior_cell,
            packet_before,
            packet_after,
            packet_source,
            external_source,
            proposed,
        ) in (
            (
                "water",
                old.totals.water_kg,
                before_water,
                after_water,
                packet_water_to_cell,
                external.water_to_cell_kg,
                proposed_water,
            ),
            (
                "hexane",
                old.totals.hexane_kg,
                before_hexane,
                after_hexane,
                packet_hexane_to_cell,
                external.hexane_to_cell_kg,
                proposed_hexane,
            ),
        ):
            if proposed < 0.0:
                details.append(
                    f"tray={request.prior.physical_tray_id}, K_layer=K{layer}, "
                    f"component={component}, prior_cell_kg={prior_cell:.17g}, "
                    f"aggregate_packet_before_kg={packet_before:.17g}, "
                    f"aggregate_packet_after_kg={packet_after:.17g}, "
                    f"aggregate_packet_to_cell_source_kg={packet_source:.17g}, "
                    f"external_to_cell_kg={external_source:.17g}, "
                    f"proposed_cell_kg={proposed:.17g}, deficit_kg={-proposed:.17g}"
                )
        if details:
            raise KCellTrayStepError(
                "negative complementary cell component inventory: " + " | ".join(details),
                request.prior,
            )
        energy_terms = request.rounded_packet_cell_energy_terms
        if energy_terms is None:
            proposed_energy = math.fsum(
                (
                    old.totals.common_datum_energy_j,
                    before_energy,
                    -after_energy,
                    external.common_datum_energy_to_cell_j,
                )
            )
        else:
            terms = energy_terms[index]
            if not _same_float(terms.packet_energy_to_cell_j, packet_energy_to_cell):
                raise KCellTrayStepError(
                    f"rounded packet-cell energy source changed at K{layer}",
                    request.prior,
                )
            if not _same_float(
                terms.packet_energy_exact_inverse_j,
                -packet_energy_to_cell,
            ):
                raise KCellTrayStepError(
                    f"rounded packet-cell energy inverse changed at K{layer}",
                    request.prior,
                )
            proposed_energy = math.fsum(
                (
                    old.totals.common_datum_energy_j,
                    packet_energy_to_cell,
                    terms.packet_energy_exact_inverse_j,
                    terms.wall_storage_energy_j,
                )
            )
            expected_energy = math.fsum(
                (old.totals.common_datum_energy_j, terms.wall_storage_energy_j)
            )
            if not _same_float(proposed_energy, expected_energy):
                raise KCellTrayStepError(
                    f"rounded packet-cell energy terms did not cancel exactly at K{layer}",
                    request.prior,
                )
        totals = CellConservedTotals(
            water_kg=proposed_water,
            hexane_kg=proposed_hexane,
            common_datum_energy_j=proposed_energy,
        )
        conserved = HostCellConservedState(
            cell_id=old.cell_id,
            revision=old.revision + 1,
            state_digest=_cell_digest(
                attempt_id=request.attempt_id,
                prior=old,
                layer=cell.vertical_layer_id,
                time_after_s=request.end_time_s,
                totals=totals,
            ),
            energy_datum_id=old.energy_datum_id,
            totals=totals,
        )
        candidates.append(
            AcceptedKCell(vertical_layer_id=cell.vertical_layer_id, conserved_state=conserved)
        )
    return tuple(candidates)


def _tray_ledger(
    request: KCellTrayStepRequest,
    *,
    advanced_packets: tuple[AcceptedTrayPacket, ...],
    moved_packets: tuple[AcceptedTrayPacket, ...],
    cells_after: tuple[AcceptedKCell, ...],
) -> WeightedTrayLedger:
    packet_before = _sum_packet_totals(request.prior.packets)
    packet_after_advance = _sum_packet_totals(advanced_packets)
    packet_after_motion = _sum_packet_totals(moved_packets)
    cell_before = _sum_cell_totals(request.prior.cells)
    cell_after = _sum_cell_totals(cells_after)
    external = (
        math.fsum(item.water_to_cell_kg for item in request.external_cell_transfers),
        math.fsum(item.hexane_to_cell_kg for item in request.external_cell_transfers),
        math.fsum(item.common_datum_energy_to_cell_j for item in request.external_cell_transfers),
    )
    dry = math.fsum((packet_after_advance[0], -packet_before[0]))
    oil = math.fsum((packet_after_advance[1], -packet_before[1]))
    water = math.fsum(
        (
            packet_after_advance[4],
            packet_after_advance[5],
            cell_after[0],
            -packet_before[4],
            -packet_before[5],
            -cell_before[0],
            -external[0],
        )
    )
    hexane = math.fsum(
        (
            packet_after_advance[2],
            packet_after_advance[3],
            cell_after[1],
            -packet_before[2],
            -packet_before[3],
            -cell_before[1],
            -external[1],
        )
    )
    energy = math.fsum(
        (
            packet_after_advance[6],
            cell_after[2],
            -packet_before[6],
            -cell_before[2],
            -external[2],
        )
    )
    motion = tuple(
        math.fsum((after, -before))
        for before, after in zip(packet_after_advance, packet_after_motion, strict=True)
    )
    relative = max(
        _relative_residual(dry, packet_before[0], packet_after_advance[0]),
        _relative_residual(oil, packet_before[1], packet_after_advance[1]),
        _relative_residual(
            water,
            packet_before[4],
            packet_before[5],
            packet_after_advance[4],
            packet_after_advance[5],
            cell_before[0],
            cell_after[0],
            external[0],
        ),
        _relative_residual(
            hexane,
            packet_before[2],
            packet_before[3],
            packet_after_advance[2],
            packet_after_advance[3],
            cell_before[1],
            cell_after[1],
            external[1],
        ),
        _relative_residual(
            energy,
            packet_before[6],
            packet_after_advance[6],
            cell_before[2],
            cell_after[2],
            external[2],
        ),
    )
    return WeightedTrayLedger(
        packet_before=packet_before,
        packet_after_advance=packet_after_advance,
        packet_after_motion=packet_after_motion,
        cell_before=cell_before,
        cell_after=cell_after,
        external_to_cells=external,
        dry_matter_residual_kg=dry,
        residual_oil_label_residual_kg=oil,
        water_residual_kg=water,
        hexane_residual_kg=hexane,
        common_datum_energy_residual_j=energy,
        motion_residuals=motion,
        maximum_relative_residual=relative,
        relative_limit=request.relative_limit,
    )


def evaluate_k_cell_tray_step(
    request: KCellTrayStepRequest,
    *,
    packet_callback: ManufacturedPacketAdvanceCallback,
    k_planner: ManufacturedKRekeyPlanner,
    auditor: WeightedPacketPayloadAuditor,
) -> ValidatedKCellTrayStep | KCellTrayStepRejection:
    """Evaluate one immutable weighted tray step without mutating accepted state."""

    if type(request) is not KCellTrayStepRequest:
        raise TypeError("K-cell host requires an exact KCellTrayStepRequest")
    for name, value in (
        ("packet callback ID", getattr(packet_callback, "callback_id", None)),
        ("packet callback source", getattr(packet_callback, "source_identity", None)),
        ("K planner ID", getattr(k_planner, "planner_id", None)),
        ("K planner source", getattr(k_planner, "source_identity", None)),
        ("auditor identity", getattr(auditor, "auditor_identity_digest", None)),
        ("auditor source", getattr(auditor, "source_identity", None)),
    ):
        if type(value) is not str or not value.strip():
            return KCellTrayStepRejection(
                attempt_id=request.attempt_id,
                reason=f"{name} is missing or blank",
                rollback_state=request.prior,
            )
    maximum = getattr(packet_callback, "maximum_scaled_residual", None)
    if type(maximum) is not float or not math.isfinite(maximum) or maximum <= 0.0:
        return KCellTrayStepRejection(
            attempt_id=request.attempt_id,
            reason="packet callback residual pin is not a positive binary64",
            rollback_state=request.prior,
        )
    if auditor.auditor_identity_digest != request.prior.auditor_identity_digest:
        return KCellTrayStepRejection(
            attempt_id=request.attempt_id,
            reason="payload auditor differs from accepted tray authority",
            rollback_state=request.prior,
        )
    try:
        clock = _clock_ledger(request)
        advanced = _advance_packets(request, callback=packet_callback, auditor=auditor)
        if type(advanced) is KCellTrayStepRejection:
            return advanced
        after_rtd, rtd_records = _apply_rtd_rekeys(request, advanced, clock)
        moved, k_records = _apply_k_rekeys(request, after_rtd, planner=k_planner)
        cells_after = _candidate_cells(request, advanced)
        ledger = _tray_ledger(
            request,
            advanced_packets=advanced,
            moved_packets=moved,
            cells_after=cells_after,
        )
        if not ledger.passed:
            raise KCellTrayStepError("weighted tray conservation ledger failed", request.prior)
        candidate = AcceptedKCellTrayState(
            topology=request.prior.topology,
            physical_tray_id=request.prior.physical_tray_id,
            revision=request.prior.revision + 1,
            time_s=request.end_time_s,
            pressure_pa=request.prior.pressure_pa,
            configuration_digest=request.prior.configuration_digest,
            codec_registry_digest=request.prior.codec_registry_digest,
            auditor_identity_digest=request.prior.auditor_identity_digest,
            energy_datum_id=request.prior.energy_datum_id,
            rtd_schedule_digest=request.prior.rtd_schedule_digest,
            cells=cells_after,
            packets=moved,
        )
        digest = _request_digest(
            request,
            packet_callback_id=packet_callback.callback_id,
            packet_callback_source_identity=packet_callback.source_identity,
            planner_id=k_planner.planner_id,
            planner_source_identity=k_planner.source_identity,
            auditor_identity_digest=auditor.auditor_identity_digest,
            auditor_source_identity=auditor.source_identity,
        )
        return ValidatedKCellTrayStep(
            request=request,
            request_digest=digest,
            candidate=candidate,
            candidate_digest=candidate.state_digest,
            clock_ledger=clock,
            tray_ledger=ledger,
            rekeys=(*rtd_records, *k_records),
            packet_callback_id=packet_callback.callback_id,
            packet_callback_source_identity=packet_callback.source_identity,
            k_planner_id=k_planner.planner_id,
            k_planner_source_identity=k_planner.source_identity,
            auditor_identity_digest=auditor.auditor_identity_digest,
            auditor_source_identity=auditor.source_identity,
            _seal=_VALIDATED_STEP_SEAL,
        )
    except KCellTrayHostError as error:
        return KCellTrayStepRejection(
            attempt_id=request.attempt_id,
            reason=str(error),
            rollback_state=request.prior,
        )
    except Exception as error:  # manufactured collaborators are untrusted host boundaries
        return KCellTrayStepRejection(
            attempt_id=request.attempt_id,
            reason=f"host evaluation failed: {type(error).__name__}: {error}",
            rollback_state=request.prior,
        )


def commit_k_cell_tray_step(
    current: AcceptedKCellTrayState,
    transaction: ValidatedKCellTrayStep,
    *,
    auditor: WeightedPacketPayloadAuditor,
) -> AcceptedKCellTrayState:
    """Return the sealed candidate after exact stale-state and payload re-audits.

    The caller still owns the lock and the eventual durable all-or-nothing write.
    This stateless function performs every reversible check before returning the
    candidate and never mutates ``current`` or any callback-owned value.
    """

    if type(current) is not AcceptedKCellTrayState:
        raise TypeError("K-cell commit requires an exact accepted tray state")
    if (
        type(transaction) is not ValidatedKCellTrayStep
        or transaction._seal is not _VALIDATED_STEP_SEAL
    ):
        raise KCellTrayHostError("K-cell commit requires a host-validated transaction")
    if current is not transaction.request.prior:
        raise StaleKCellTrayStateError("K-cell transaction is stale or foreign")
    if current.state_digest != transaction.request.prior.state_digest:
        raise StaleKCellTrayStateError("accepted prior digest changed before commit")
    if auditor.auditor_identity_digest != transaction.auditor_identity_digest:
        raise KCellTrayHostError("commit payload auditor violates transaction pin")
    if auditor.source_identity != transaction.auditor_source_identity:
        raise KCellTrayHostError("commit payload auditor source violates transaction pin")
    for packet in current.packets:
        _audit_entry(packet, auditor, label="commit-prior")
    for packet in transaction.candidate.packets:
        _audit_entry(packet, auditor, label="commit-candidate")
    expected_request_digest = _request_digest(
        transaction.request,
        packet_callback_id=transaction.packet_callback_id,
        packet_callback_source_identity=transaction.packet_callback_source_identity,
        planner_id=transaction.k_planner_id,
        planner_source_identity=transaction.k_planner_source_identity,
        auditor_identity_digest=transaction.auditor_identity_digest,
        auditor_source_identity=transaction.auditor_source_identity,
    )
    if expected_request_digest != transaction.request_digest:
        raise KCellTrayHostError("K-cell request identity changed after evaluation")
    if transaction.candidate.state_digest != transaction.candidate_digest:
        raise KCellTrayHostError("K-cell candidate identity changed after evaluation")
    if transaction.candidate.revision != current.revision + 1:
        raise KCellTrayHostError("K-cell candidate revision did not advance exactly once")
    if not transaction.tray_ledger.passed:
        raise KCellTrayHostError("K-cell candidate conservation certificate failed")
    return transaction.candidate


__all__ = (
    "AFFINE_RESIDENCE_CLOCK_ALGORITHM_ID",
    "AFFINE_RESIDENCE_CLOCK_LAW_ID",
    "AFFINE_RESIDENCE_DECIMAL_ORACLE_ABSOLUTE_FLOOR",
    "AFFINE_RESIDENCE_DECIMAL_ORACLE_ULP_LIMIT",
    "AFFINE_RESIDENCE_EXPECTED_ALGORITHM_DEFINITION_DIGEST",
    "AFFINE_RESIDENCE_EXPECTED_LAW_DEFINITION_DIGEST",
    "AFFINE_RESIDENCE_EXPECTED_SCHEMA_DEFINITION_DIGEST",
    "AFFINE_RESIDENCE_LOG1P_RELATIVE_LIMIT",
    "AFFINE_RESIDENCE_MAX_FINITE_MASS_RATIO",
    "AFFINE_RESIDENCE_MIN_NORMAL_MASS_RATIO",
    "AFFINE_RESIDENCE_SERIES_ORDER",
    "AFFINE_RESIDENCE_SMALL_ARGUMENT_LIMIT",
    "AFFINE_RESIDENCE_TRAJECTORY_SCHEMA_ID",
    "AFFINE_RESIDENCE_TRAJECTORY_SCHEMA_REVISION",
    "AcceptedDryMatterFlowSegment",
    "AcceptedKCell",
    "AcceptedKCellTrayState",
    "AcceptedResidenceClockLedger",
    "AcceptedTrayPacket",
    "AffineDryInventoryResidenceSegment",
    "AffineDryInventoryResidenceTrajectory",
    "DEFAULT_RELATIVE_LIMIT",
    "HOST_DIGEST_DOMAIN",
    "KCellTrayHostError",
    "KCellTrayStateError",
    "KCellTrayStepError",
    "KCellTrayStepRejection",
    "KCellTrayStepRequest",
    "KLayerRekeyPlanningRequest",
    "KLayerRekeyProposal",
    "ManufacturedKRekeyPlanner",
    "ManufacturedPacketAdvanceCallback",
    "ManufacturedPacketAdvanceProposal",
    "ManufacturedPacketAdvanceRejection",
    "ManufacturedPacketAdvanceRequest",
    "NoKRekeyPlanner",
    "RoundedPacketCellEnergyTerms",
    "StaleKCellTrayStateError",
    "ValidatedKCellTrayStep",
    "WeightedPacketPayloadAuditor",
    "WeightedTrayLedger",
    "WholePacketRekeyMechanism",
    "WholePacketRekeyRecord",
    "canonical_k_cell_tray_state_digest",
    "commit_k_cell_tray_step",
    "evaluate_k_cell_tray_step",
    "accepted_dry_matter_flow_history_digest",
    "validate_affine_residence_clock_definitions",
)
