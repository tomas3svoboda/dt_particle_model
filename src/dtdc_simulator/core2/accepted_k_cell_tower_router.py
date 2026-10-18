"""In-memory conservative population routing for accepted weighted K-cell trays.

This companion host closes only the population-ownership seam left open by
``k_cell_tray_host``.  It removes final-RTD-stage packet weight from each source,
admits the identical payload/inventory/weight at the next serial tray with an
explicit manufactured RTD assignment, and exposes the last-tray debit as a sink
packet plan.  Boundary splits retain the source packet as the remainder and issue
one deterministically named transfer child with identical payload bytes.

The implementation is deliberately narrow.  It is in-memory, single-host, and
same-context/same-pressure only.  It supplies no solid-motion law, accepted-flow
instrumentation, physical K scheduler, mixed-temperature closure, pressure remap,
PD dome closure, variant handoff, gas composition, or durable persistence.  Every
physical and production claim therefore remains false.
"""

from __future__ import annotations

import hashlib
import math
import struct
from dataclasses import dataclass, replace
from enum import Enum
from fractions import Fraction
from typing import ClassVar, Protocol

from .accepted_solid_transfer_host import (
    PacketPlanEntry,
    TransferAuthorityClass,
    TransferTopology,
)
from .k_cell_tray_host import (
    AcceptedKCell,
    AcceptedKCellTrayState,
    AcceptedTrayPacket,
    WeightedPacketPayloadAuditor,
)
from .packet_rtd_oracle import PacketRTDPath
from .packet_weight_scaling_contract import (
    PacketWeight,
    PerParticleInventory,
    derive_packet_totals,
    require_single_surface_authority,
)
from .tray_particle_stateful import ParticlePopulationKey, RTDStageId, VerticalLayerId

ROUTER_ID = "dtdc-core2-accepted-k-cell-tower-router-v1"
ROUTER_CONSERVATION_AUDIT_LAW_ID = (
    "fraction-from-float-signed-entry-totals-plus-certified-split-correction-v1"
)
#: Version bump ruled by GT_PS2_WALLED_FLOOR_AND_V3_APPROVAL_RULING_2026-08-29:
#: the bump deliberately changes every transaction digest, including
#: pure-legacy cases, so regression evidence compares numerical route and
#: split selections, never digests.
QSC_K_SHARE_APPORTIONMENT_POLICY_ID = "qsc-k-share-apportionment-representability-v3"
#: Sealed prior-lineage identifier: v2-era persisted history rows remain
#: readable under this identifier per the versioned history-lineage
#: migration rule of the same ruling; it never labels a new transaction.
QSC_K_SHARE_APPORTIONMENT_POLICY_ID_V2 = "qsc-k-share-apportionment-representability-v2"
#: Owner-selected closing-candidate order law (ruling question b):
#: nearest-to-ideal exact-rational distance, ties broken by the recorded
#: encounter rank of the typed exhaustion's one-ULP fallback candidates.
QSC_K_SHARE_CLOSING_EXCHANGE_ORDER_LAW = (
    "nearest-to-ideal-exact-rational-distance-ties-recorded-encounter-rank"
)
#: Owner-declared FINAL per-lineage drift-tripwire budget
#: (GT_PS2_V3_TRIPWIRE_BUDGET_RULING_2026-08-29: "Ok, hereby I declare 3
#: final."); any re-declaration is a recorded owner amendment that re-runs
#: all evidence delivered under the superseded number.
QSC_K_SHARE_RULED_DRIFT_TRIPWIRE_ULP = 3  # GT_PS2_V3_TRIPWIRE_BUDGET_RULING_2026-08-29
_DIGEST_DOMAIN = ROUTER_ID.encode("ascii")
_VALIDATED_SEAL = object()


class KShareApportionmentDisposition(str, Enum):
    """Owner-ruled disposition for one tray's represented K-share vector.

    ``CLOSING_ONE_ULP_EXCHANGE`` is the v3 closing exchange approved by
    GT_PS2_WALLED_FLOOR_AND_V3_APPROVAL_RULING_2026-08-29: the declared
    closing share moved by exactly one ULP with exactly one nonclosing
    layer carrying the exact rational complement.
    """

    LEGACY_ZERO = "legacy_zero"
    ONE_NONCLOSING_ONE_ULP = "one_nonclosing_one_ulp"
    CLOSING_ONE_ULP_EXCHANGE = "closing_one_ulp_exchange"


class KShareApportionmentRole(str, Enum):
    """Role of one K row in the owner-ruled share construction.

    Under the v3 closing exchange the ``DECLARED_CLOSING`` role may carry a
    nonzero signed displacement from its legacy represented target; in every
    other disposition it remains at zero displacement.
    """

    UNCHANGED = "unchanged"
    ADJUSTED_NONCLOSING = "adjusted_nonclosing"
    DECLARED_CLOSING = "declared_closing"


@dataclass(frozen=True, slots=True, kw_only=True)
class KShareClosingExchangePolicy:
    """Owner-issued immutable v3 closing-exchange policy: no policy, no v3.

    Approved by GT_PS2_WALLED_FLOOR_AND_V3_APPROVAL_RULING_2026-08-29.  The
    policy identity is bound into the routing request digest, so a v3
    construction without an issued policy is replay-provably impossible.
    The tripwire budget is an UNDERIVED owner-declared per-campaign cap on
    any single K's cumulative signed discrepancy magnitude in closing-share
    ULP counts; exceeding it is a typed escalation on the history surface.
    """

    policy_id: str
    candidate_order_law: str
    ruling_provenance: str
    cumulative_signed_drift_tripwire_ulp: int

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.policy_id) is not str or self.policy_id != QSC_K_SHARE_APPORTIONMENT_POLICY_ID:
            raise KCellTowerRouterConfigurationError(
                "closing-exchange policy identity differs from the ruled v3 identifier"
            )
        if (
            type(self.candidate_order_law) is not str
            or self.candidate_order_law != QSC_K_SHARE_CLOSING_EXCHANGE_ORDER_LAW
        ):
            raise KCellTowerRouterConfigurationError(
                "closing-exchange candidate order differs from the owner-selected law"
            )
        if type(self.ruling_provenance) is not str or not self.ruling_provenance.strip():
            raise KCellTowerRouterConfigurationError(
                "closing-exchange policy requires its ruling provenance"
            )
        if (
            type(self.cumulative_signed_drift_tripwire_ulp) is not int
            or self.cumulative_signed_drift_tripwire_ulp < 1
        ):
            raise KCellTowerRouterConfigurationError(
                "closing-exchange drift tripwire must be a positive exact integer"
            )


def make_ruled_closing_exchange_policy(
    *,
    cumulative_signed_drift_tripwire_ulp: int = QSC_K_SHARE_RULED_DRIFT_TRIPWIRE_ULP,
) -> KShareClosingExchangePolicy:
    return KShareClosingExchangePolicy(
        policy_id=QSC_K_SHARE_APPORTIONMENT_POLICY_ID,
        candidate_order_law=QSC_K_SHARE_CLOSING_EXCHANGE_ORDER_LAW,
        ruling_provenance="GT_PS2_WALLED_FLOOR_AND_V3_APPROVAL_RULING_2026-08-29",
        cumulative_signed_drift_tripwire_ulp=cumulative_signed_drift_tripwire_ulp,
    )


class KCellTowerRouterError(RuntimeError):
    """Base error for the in-memory serial population router."""


class KCellTowerRouterConfigurationError(KCellTowerRouterError):
    """The requested serial route is outside the bounded router contract."""


class KCellTowerRouterConservationError(KCellTowerRouterError):
    """A split, exit, edge, or whole-router extensive ledger did not close."""


class _KShareConstructionExhaustion(KCellTowerRouterConservationError):
    """One K target exhausted the unchanged exact one-child construction."""

    def __init__(
        self,
        *,
        layer_id: VerticalLayerId,
        represented_target_dry_kg: float,
    ) -> None:
        super().__init__(
            "boundary split cannot represent the exact K-share construction with "
            "exact parent-weight closure"
        )
        self.layer_id = layer_id
        self.represented_target_dry_kg = represented_target_dry_kg


class _KShareRepresentabilityExhaustion(_KShareConstructionExhaustion):
    """One nonclosing legacy target exposes ordered adjacent represented targets."""

    def __init__(
        self,
        *,
        layer_id: VerticalLayerId,
        represented_target_dry_kg: float,
        one_ulp_fallback_candidates: tuple[tuple[float, int], ...],
    ) -> None:
        if not one_ulp_fallback_candidates:
            raise ValueError("representability fallback exhaustion requires adjacent targets")
        super().__init__(
            layer_id=layer_id,
            represented_target_dry_kg=represented_target_dry_kg,
        )
        self.one_ulp_fallback_candidates = one_ulp_fallback_candidates


class StaleKCellTowerRoutingError(KCellTowerRouterError):
    """The current tray tuple is not the tuple against which routing was evaluated."""


class ReplayedKCellTowerRoutingError(KCellTowerRouterError):
    """An already consumed in-memory routing transaction was presented again."""


def _require_nonblank(label: str, *values: str) -> None:
    if any(type(value) is not str or not value.strip() for value in values):
        raise KCellTowerRouterConfigurationError(f"{label} must be nonblank exact strings")


def _require_binary64(label: str, *values: float, nonnegative: bool = False) -> None:
    if any(type(value) is not float or not math.isfinite(value) for value in values):
        raise KCellTowerRouterConfigurationError(f"{label} must be finite binary64 values")
    if nonnegative and any(value < 0.0 or _is_negative_zero(value) for value in values):
        raise KCellTowerRouterConfigurationError(
            f"{label} must be nonnegative canonical binary64 values"
        )


def _is_negative_zero(value: float) -> bool:
    return value == 0.0 and math.copysign(1.0, value) < 0.0


def _same_float(left: float, right: float) -> bool:
    return struct.pack(">d", left) == struct.pack(">d", right)


def _signed_nonnegative_ulp_distance(value: float, reference: float) -> int:
    """Return the exact binary64 rank displacement on the nonnegative domain."""

    _require_binary64("ULP-distance operand", value, reference, nonnegative=True)
    value_bits = struct.unpack(">Q", struct.pack(">d", value))[0]
    reference_bits = struct.unpack(">Q", struct.pack(">d", reference))[0]
    return value_bits - reference_bits


def _float_part(value: float) -> bytes:
    return struct.pack(">d", value)


def _integer_part(value: int) -> bytes:
    return value.to_bytes(8, "big", signed=False)


def _signed_integer_part(value: int) -> bytes:
    return value.to_bytes(8, "big", signed=True)


def _text_part(value: str) -> bytes:
    return value.encode("utf-8")


def _fraction_part(value: Fraction) -> bytes:
    if type(value) is not Fraction:
        raise TypeError("exact conservation digest values must be exact Fractions")
    return f"{value.numerator}/{value.denominator}".encode("ascii")


def _digest(*parts: bytes) -> str:
    digest = hashlib.sha256()
    digest.update(_DIGEST_DOMAIN)
    digest.update(len(parts).to_bytes(8, "big"))
    for part in parts:
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


def _payload_identity(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def _entry_parts(entry: PacketPlanEntry) -> tuple[bytes, ...]:
    surface_parts = (
        (
            b"surface",
            _float_part(entry.surface.attached_hexane_kg),
            _float_part(entry.surface.free_water_kg),
        )
        if entry.surface is not None
        else (b"no-surface",)
    )
    return (
        _text_part(entry.packet_id),
        _text_part(_payload_identity(entry.payload_bytes)),
        *(_float_part(value) for value in entry.inventory.as_tuple()),
        _text_part(entry.inventory.energy_datum_id),
        _float_part(entry.weight.representative_particles),
        _text_part(entry.weight.energy_datum_id),
        *surface_parts,
    )


def _packet_sort_key(packet: AcceptedTrayPacket) -> tuple[int, int, str]:
    return (
        packet.owner_key.vertical_layer_id.value,
        packet.owner_key.rtd_stage_id.value,
        packet.packet_id,
    )


def _sum_entries(entries: tuple[PacketPlanEntry, ...]) -> tuple[float, ...]:
    columns: tuple[list[float], ...] = tuple([] for _ in range(7))
    for entry in entries:
        totals = derive_packet_totals(entry.inventory, entry.weight).totals
        for column, value in zip(columns, totals, strict=True):
            column.append(value)
    return tuple(math.fsum(column) for column in columns)


def _sum_trays(trays: tuple[AcceptedKCellTrayState, ...]) -> tuple[float, ...]:
    return _sum_entries(tuple(packet.entry for tray in trays for packet in tray.packets))


def _zero_extensives() -> tuple[float, ...]:
    return (0.0,) * 7


def _binary64_regrouping_residual(
    *, credits: tuple[tuple[float, ...], ...], debits: tuple[tuple[float, ...], ...]
) -> tuple[float, ...]:
    return tuple(
        math.fsum(
            (
                *(row[index] for row in credits),
                *(-row[index] for row in debits),
            )
        )
        for index in range(7)
    )


def _exact_signed_entry_residual(
    *,
    credits: tuple[tuple[PacketPlanEntry, ...], ...],
    debits: tuple[tuple[PacketPlanEntry, ...], ...],
) -> tuple[Fraction, ...]:
    """Sum the signed packet-total binary64 operands without regrouping.

    ``derive_packet_totals`` remains the sole multiplicity authority.  Each of
    its already-rounded binary64 outputs is lifted losslessly with
    ``Fraction.from_float`` before any addition.  Consequently, ownership
    regrouping cannot alter the authoritative residual while a one-ULP change
    to any represented packet total remains visible.
    """

    residual = [Fraction(0, 1) for _ in range(7)]
    for sign, groups in ((1.0, credits), (-1.0, debits)):
        for entries in groups:
            for entry in entries:
                for index, value in enumerate(
                    derive_packet_totals(entry.inventory, entry.weight).totals
                ):
                    residual[index] += Fraction.from_float(sign * value)
    return tuple(residual)


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedFeedPacket:
    """One externally prepared packet admitted at the first serial tray."""

    entry: PacketPlanEntry
    pressure_pa: float
    configuration_digest: str
    codec_registry_digest: str
    auditor_identity_digest: str

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    feed_instrumentation_authenticated: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.entry) is not PacketPlanEntry:
            raise KCellTowerRouterConfigurationError("feed entry must be an exact PacketPlanEntry")
        _require_binary64("feed pressure", self.pressure_pa)
        if self.pressure_pa <= 0.0:
            raise KCellTowerRouterConfigurationError("feed pressure must be positive")
        _require_nonblank(
            "feed context identities",
            self.configuration_digest,
            self.codec_registry_digest,
            self.auditor_identity_digest,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedTrayExitDebit:
    """Integrated accepted dry-matter debit for one tray routing boundary."""

    physical_tray_id: str
    accepted_dry_matter_kg: float
    closing_layer_cursor: int = 1

    physically_qualifying: ClassVar[bool] = False
    accepted_flow_instrumentation_authenticated: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank("exit tray ID", self.physical_tray_id)
        _require_binary64(
            "accepted dry-matter debit", self.accepted_dry_matter_kg, nonnegative=True
        )
        if type(self.closing_layer_cursor) is not int or self.closing_layer_cursor <= 0:
            raise KCellTowerRouterConfigurationError(
                "closing_layer_cursor must be a positive exact integer"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class AdmissionPlanningRequest:
    """Complete incoming plan presented to one manufactured admission planner."""

    attempt_id: str
    destination_prior: AcceptedKCellTrayState
    incoming_entries: tuple[PacketPlanEntry, ...]

    physically_qualifying: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class AdmissionAssignment:
    """Explicit destination K layer and new tray-local RTD path for one packet."""

    packet_id: str
    destination_layer_id: VerticalLayerId
    destination_rtd_path: PacketRTDPath
    authority_id: str

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank("admission packet/authority", self.packet_id, self.authority_id)
        if type(self.destination_layer_id) is not VerticalLayerId:
            raise KCellTowerRouterConfigurationError(
                "admission layer must be an exact VerticalLayerId"
            )
        if type(self.destination_rtd_path) is not PacketRTDPath:
            raise KCellTowerRouterConfigurationError(
                "admission RTD path must be an exact PacketRTDPath"
            )


class ManufacturedAdmissionPlanner(Protocol):
    planner_id: str
    source_identity: str

    def assign(self, request: AdmissionPlanningRequest) -> tuple[AdmissionAssignment, ...]: ...


@dataclass(frozen=True, slots=True, kw_only=True)
class SerialTowerRoutingRequest:
    """One simultaneous feed/inter-tray/final-sink population routing attempt."""

    attempt_id: str
    event_time_s: float
    prior_trays: tuple[AcceptedKCellTrayState, ...]
    exit_debits: tuple[AcceptedTrayExitDebit, ...]
    feed_packets: tuple[AcceptedFeedPacket, ...] = ()
    #: Owner-issued v3 closing-exchange authority; absent means the closing
    #: exhaustion keeps its typed refusal (no policy, no v3).
    k_share_closing_exchange_policy: KShareClosingExchangePolicy | None = None

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank("routing attempt ID", self.attempt_id)
        _require_binary64("routing event time", self.event_time_s, nonnegative=True)
        if self.k_share_closing_exchange_policy is not None:
            if type(self.k_share_closing_exchange_policy) is not KShareClosingExchangePolicy:
                raise KCellTowerRouterConfigurationError(
                    "closing-exchange policy carrier has a foreign exact type"
                )
            self.k_share_closing_exchange_policy.__post_init__()
        if type(self.prior_trays) is not tuple or not self.prior_trays:
            raise KCellTowerRouterConfigurationError("routing requires a nonempty exact tray tuple")
        if any(type(tray) is not AcceptedKCellTrayState for tray in self.prior_trays):
            raise KCellTowerRouterConfigurationError(
                "routing prior must contain exact AcceptedKCellTrayState values"
            )
        if type(self.exit_debits) is not tuple or len(self.exit_debits) != len(self.prior_trays):
            raise KCellTowerRouterConfigurationError("routing needs exactly one debit per tray")
        if any(type(debit) is not AcceptedTrayExitDebit for debit in self.exit_debits):
            raise KCellTowerRouterConfigurationError(
                "exit_debits must contain exact AcceptedTrayExitDebit values"
            )
        if type(self.feed_packets) is not tuple or any(
            type(packet) is not AcceptedFeedPacket for packet in self.feed_packets
        ):
            raise KCellTowerRouterConfigurationError(
                "feed_packets must contain exact AcceptedFeedPacket values"
            )
        tray_ids = tuple(tray.physical_tray_id for tray in self.prior_trays)
        if len(set(tray_ids)) != len(tray_ids):
            raise KCellTowerRouterConfigurationError("serial tray IDs must be unique")
        if tuple(debit.physical_tray_id for debit in self.exit_debits) != tray_ids:
            raise KCellTowerRouterConfigurationError(
                "exit debits must follow the exact declared serial tray order"
            )
        first = self.prior_trays[0]
        for tray, debit in zip(self.prior_trays, self.exit_debits, strict=True):
            if not _same_float(tray.time_s, self.event_time_s):
                raise KCellTowerRouterConfigurationError(
                    "routing occurs at an already accepted common tray time"
                )
            physical = next(
                item
                for item in tray.topology.trays
                if item.physical_tray_id == tray.physical_tray_id
            )
            if debit.closing_layer_cursor > physical.vertical_layer_count:
                raise KCellTowerRouterConfigurationError(
                    "closing-layer cursor exceeds the source tray K"
                )
            if tray.topology.digest != first.topology.digest:
                raise KCellTowerRouterConfigurationError(
                    "all serial trays must bind one exact topology authority"
                )
            for label, observed, expected in (
                ("pressure", tray.pressure_pa, first.pressure_pa),
                (
                    "configuration",
                    tray.configuration_digest,
                    first.configuration_digest,
                ),
                (
                    "codec registry",
                    tray.codec_registry_digest,
                    first.codec_registry_digest,
                ),
                (
                    "auditor",
                    tray.auditor_identity_digest,
                    first.auditor_identity_digest,
                ),
                ("energy datum", tray.energy_datum_id, first.energy_datum_id),
            ):
                equal = (
                    _same_float(observed, expected) if label == "pressure" else observed == expected
                )
                if not equal:
                    raise KCellTowerRouterConfigurationError(
                        f"bounded router refuses cross-tray {label} remapping"
                    )
        for feed in self.feed_packets:
            if not _same_float(feed.pressure_pa, first.pressure_pa):
                raise KCellTowerRouterConfigurationError(
                    "bounded router refuses feed pressure remapping"
                )
            if (
                feed.configuration_digest != first.configuration_digest
                or feed.codec_registry_digest != first.codec_registry_digest
                or feed.auditor_identity_digest != first.auditor_identity_digest
                or feed.entry.inventory.energy_datum_id != first.energy_datum_id
            ):
                raise KCellTowerRouterConfigurationError(
                    "feed packet context differs from the first tray context"
                )
        identifiers = tuple(
            packet.packet_id for tray in self.prior_trays for packet in tray.packets
        ) + tuple(packet.entry.packet_id for packet in self.feed_packets)
        if len(identifiers) != len(set(identifiers)):
            raise KCellTowerRouterConfigurationError(
                "packet identities must be globally unique across trays and feed"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class PacketSplitLineage:
    parent_packet_id: str
    source_remainder_packet_id: str
    transfer_child_packet_id: str
    source_tray_id: str
    source_layer_id: VerticalLayerId
    payload_identity: str
    parent_weight: float
    remainder_weight: float
    transfer_weight: float

    physically_qualifying: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class KShareExitLedger:
    layer_id: VerticalLayerId
    eligible_dry_matter_kg: float
    ideal_dry_matter_debit_kg: float
    legacy_represented_dry_matter_debit_kg: float
    represented_dry_matter_debit_kg: float
    represented_minus_ideal_kg: float
    signed_binary64_steps_from_legacy: int
    apportionment_role: KShareApportionmentRole
    selected_boundary_weight_candidate_rank: int | None
    whole_packet_count: int
    boundary_split_count: int

    physically_qualifying: ClassVar[bool] = False

    @property
    def exact_represented_minus_ideal_kg(self) -> Fraction:
        """Return the exact signed discrepancy between two binary64 operands."""

        return Fraction.from_float(self.represented_dry_matter_debit_kg) - Fraction.from_float(
            self.ideal_dry_matter_debit_kg
        )

    @property
    def exact_absolute_represented_minus_ideal_kg(self) -> Fraction:
        """Return the exact magnitude of the represented-share discrepancy."""

        return abs(self.exact_represented_minus_ideal_kg)


@dataclass(frozen=True, slots=True, kw_only=True)
class TrayExitLedger:
    physical_tray_id: str
    apportionment_policy_id: str
    apportionment_disposition: KShareApportionmentDisposition
    adjusted_nonclosing_layer_id: VerticalLayerId | None
    selected_nonclosing_signed_ulp_displacement: int
    inherited_direction_candidate_rank: int | None
    accepted_dry_matter_debit_kg: float
    represented_dry_matter_debit_kg: float
    eligible_dry_matter_kg: float
    closing_layer_id: VerticalLayerId | None
    requested_closing_layer_cursor: int
    next_closing_layer_cursor: int
    k_shares: tuple[KShareExitLedger, ...]
    passed: bool
    #: v3 closing-exchange selection evidence
    #: (GT_PS2_WALLED_FLOOR_AND_V3_APPROVAL_RULING_2026-08-29); both stay at
    #: their zero/absent values for every legacy and v2 disposition.
    selected_closing_signed_ulp_displacement: int = 0
    inherited_closing_candidate_rank: int | None = None

    accepted_debit_is_booked_boundary_flow: ClassVar[bool] = True
    flat_packet_debit_is_diagnostic_only: ClassVar[bool] = True
    whole_router_exact_zero_is_commit_authority: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False

    @property
    def exact_flattened_minus_accepted_dry_matter_kg(self) -> Fraction:
        """Return the exact binary64 association diagnostic ``F - A``."""

        return Fraction.from_float(self.represented_dry_matter_debit_kg) - Fraction.from_float(
            self.accepted_dry_matter_debit_kg
        )

    @property
    def representation_fallback_layer_id(self) -> VerticalLayerId | None:
        """Identify the sole nonclosing K share moved by the v2 policy."""

        displaced = tuple(
            share.layer_id
            for share in self.k_shares
            if share.layer_id != self.closing_layer_id
            and share.exact_represented_minus_ideal_kg != 0
        )
        if len(displaced) > 1:
            raise KCellTowerRouterConservationError(
                "K-share ledger declares more than one displaced nonclosing layer"
            )
        observed = displaced[0] if displaced else None
        if observed != self.adjusted_nonclosing_layer_id:
            raise KCellTowerRouterConservationError(
                "K-share fallback layer differs from the recorded apportionment selection"
            )
        return observed

    @property
    def representation_fallback_applied(self) -> bool:
        """Report whether the bounded v2 nonclosing fallback was used."""

        return self.representation_fallback_layer_id is not None


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedPopulationRoute:
    authority_class: TransferAuthorityClass
    topology: TransferTopology
    packet_plan: tuple[PacketPlanEntry, ...]
    extensive_totals: tuple[float, ...]

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    durable_commit_selected: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class WholeRouterLedger:
    audit_law_id: str
    tray_before_totals: tuple[float, ...]
    feed_totals: tuple[float, ...]
    tray_after_totals: tuple[float, ...]
    final_sink_totals: tuple[float, ...]
    raw_exact_signed_entry_residual: tuple[Fraction, ...]
    certified_split_rounding_correction: tuple[Fraction, ...]
    conservation_residual: tuple[Fraction, ...]
    binary64_regrouping_residual: tuple[float, ...]
    route_count: int
    passed: bool

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self.audit_law_id != ROUTER_CONSERVATION_AUDIT_LAW_ID:
            raise KCellTowerRouterConfigurationError(
                "whole-router ledger carries a foreign conservation audit law"
            )
        for label, values in (
            ("tray-before totals", self.tray_before_totals),
            ("feed totals", self.feed_totals),
            ("tray-after totals", self.tray_after_totals),
            ("final-sink totals", self.final_sink_totals),
            ("binary64 regrouping residual", self.binary64_regrouping_residual),
        ):
            if type(values) is not tuple or len(values) != 7:
                raise KCellTowerRouterConfigurationError(
                    f"{label} must be an exact seven-value tuple"
                )
            _require_binary64(label, *values)
        for label, values in (
            ("raw exact signed-entry residual", self.raw_exact_signed_entry_residual),
            (
                "certified split-rounding correction",
                self.certified_split_rounding_correction,
            ),
            ("authoritative conservation residual", self.conservation_residual),
        ):
            if (
                type(values) is not tuple
                or len(values) != 7
                or any(type(value) is not Fraction for value in values)
            ):
                raise KCellTowerRouterConfigurationError(f"{label} must be seven exact Fractions")
        recomposed = tuple(
            raw + correction
            for raw, correction in zip(
                self.raw_exact_signed_entry_residual,
                self.certified_split_rounding_correction,
                strict=True,
            )
        )
        if self.conservation_residual != recomposed:
            raise KCellTowerRouterConfigurationError(
                "authoritative residual must equal raw exact residual plus certified correction"
            )
        if type(self.route_count) is not int or self.route_count < 0:
            raise KCellTowerRouterConfigurationError(
                "whole-router route count must be a nonnegative exact integer"
            )
        exact_zero = all(value == 0 for value in self.conservation_residual)
        if type(self.passed) is not bool or self.passed is not exact_zero:
            raise KCellTowerRouterConfigurationError(
                "whole-router pass flag must equal the authoritative exact-zero gate"
            )


def _same_split_envelope(left: PacketPlanEntry, right: PacketPlanEntry) -> bool:
    if left.payload_bytes != right.payload_bytes:
        return False
    if left.inventory.energy_datum_id != right.inventory.energy_datum_id or any(
        not _same_float(a, b)
        for a, b in zip(left.inventory.as_tuple(), right.inventory.as_tuple(), strict=True)
    ):
        return False
    if (left.surface is None) != (right.surface is None):
        return False
    if left.surface is not None and right.surface is not None:
        if not _same_float(left.surface.attached_hexane_kg, right.surface.attached_hexane_kg):
            return False
        if not _same_float(left.surface.free_water_kg, right.surface.free_water_kg):
            return False
    return left.weight.energy_datum_id == right.weight.energy_datum_id


def _deterministic_split_totals(entry: PacketPlanEntry) -> tuple[float, ...]:
    """Return, and independently verify, the seven ruled weight products."""

    totals = derive_packet_totals(entry.inventory, entry.weight).totals
    expected = tuple(
        value * entry.weight.representative_particles for value in entry.inventory.as_tuple()
    )
    if totals != expected:
        raise KCellTowerRouterConservationError(
            "split packet totals are not the deterministic inventory-times-weight products"
        )
    return totals


def _certified_split_rounding_correction(
    *,
    attempt_id: str,
    prior_trays: tuple[AcceptedKCellTrayState, ...],
    candidate_trays: tuple[AcceptedKCellTrayState, ...],
    final_sink_entries: tuple[PacketPlanEntry, ...],
    routes: tuple[AcceptedPopulationRoute, ...],
    exit_ledgers: tuple[TrayExitLedger, ...],
    split_lineage: tuple[PacketSplitLineage, ...],
) -> tuple[Fraction, ...]:
    """Derive, never declare, exact corrections for certified float splits."""

    declared_counts = {
        (ledger.physical_tray_id, share.layer_id.value): share.boundary_split_count
        for ledger in exit_ledgers
        for share in ledger.k_shares
    }
    observed_counts = {key: 0 for key in declared_counts}
    correction = [Fraction(0, 1) for _ in range(7)]
    seen_parents: set[tuple[str, str]] = set()
    seen_children: set[str] = set()
    for lineage in split_lineage:
        if type(lineage) is not PacketSplitLineage:
            raise KCellTowerRouterConservationError(
                "split correction requires exact PacketSplitLineage records"
            )
        parent_key = (lineage.source_tray_id, lineage.parent_packet_id)
        if parent_key in seen_parents or lineage.transfer_child_packet_id in seen_children:
            raise KCellTowerRouterConservationError(
                "split correction refuses duplicate parent or transfer-child lineage"
            )
        seen_parents.add(parent_key)
        seen_children.add(lineage.transfer_child_packet_id)
        count_key = (lineage.source_tray_id, lineage.source_layer_id.value)
        if count_key not in observed_counts:
            raise KCellTowerRouterConservationError(
                "split lineage names a source tray/layer absent from its exit ledger"
            )
        observed_counts[count_key] += 1
        if lineage.source_remainder_packet_id != lineage.parent_packet_id:
            raise KCellTowerRouterConservationError(
                "certified split remainder must retain the parent packet identity"
            )
        expected_child_id = _split_child_id(
            attempt_id=attempt_id,
            tray_id=lineage.source_tray_id,
            packet_id=lineage.parent_packet_id,
        )
        if lineage.transfer_child_packet_id != expected_child_id:
            raise KCellTowerRouterConservationError(
                "split transfer-child identity is not the canonical router identity"
            )
        parents = tuple(
            packet
            for tray in prior_trays
            if tray.physical_tray_id == lineage.source_tray_id
            for packet in tray.packets
            if packet.packet_id == lineage.parent_packet_id
        )
        remainders = tuple(
            packet
            for tray in candidate_trays
            if tray.physical_tray_id == lineage.source_tray_id
            for packet in tray.packets
            if packet.packet_id == lineage.source_remainder_packet_id
        )
        candidate_children = tuple(
            (tray.physical_tray_id, packet.entry)
            for tray in candidate_trays
            for packet in tray.packets
            if packet.packet_id == lineage.transfer_child_packet_id
        )
        sink_children = tuple(
            ("external-final-sink", entry)
            for entry in final_sink_entries
            if entry.packet_id == lineage.transfer_child_packet_id
        )
        owned_children = candidate_children + sink_children
        route_children = tuple(
            (route, entry)
            for route in routes
            for entry in route.packet_plan
            if entry.packet_id == lineage.transfer_child_packet_id
        )
        if not len(parents) == len(remainders) == len(owned_children) == len(route_children) == 1:
            raise KCellTowerRouterConservationError(
                "split correction requires one exact parent, remainder, owned child, and route"
            )
        parent_packet = parents[0]
        remainder_packet = remainders[0]
        child_destination, child_entry = owned_children[0]
        child_route, routed_child_entry = route_children[0]
        if (
            parent_packet.owner_key.vertical_layer_id != lineage.source_layer_id
            or remainder_packet.owner_key.vertical_layer_id != lineage.source_layer_id
        ):
            raise KCellTowerRouterConservationError(
                "split parent or remainder differs from the certified source layer"
            )
        parent_entry = parent_packet.entry
        remainder_entry = remainder_packet.entry
        if not (
            _same_split_envelope(parent_entry, remainder_entry)
            and _same_split_envelope(parent_entry, child_entry)
            and _same_split_envelope(child_entry, routed_child_entry)
        ):
            raise KCellTowerRouterConservationError(
                "split parent, remainder, owned child, and routed child envelopes differ"
            )
        if lineage.payload_identity != _payload_identity(parent_entry.payload_bytes):
            raise KCellTowerRouterConservationError(
                "split lineage payload identity differs from its parent"
            )
        if not (
            _same_float(lineage.parent_weight, parent_entry.weight.representative_particles)
            and _same_float(
                lineage.remainder_weight,
                remainder_entry.weight.representative_particles,
            )
            and _same_float(
                lineage.transfer_weight,
                child_entry.weight.representative_particles,
            )
            and _same_float(
                child_entry.weight.representative_particles,
                routed_child_entry.weight.representative_particles,
            )
        ):
            raise KCellTowerRouterConservationError(
                "split lineage weights differ from the certified packet entries"
            )
        if math.fsum((lineage.remainder_weight, lineage.transfer_weight)) != (
            lineage.parent_weight
        ):
            raise KCellTowerRouterConservationError(
                "split weights do not satisfy the exact certified binary64 closure"
            )
        expected_authority = (
            TransferAuthorityClass.DT_TO_DC_HANDOFF
            if child_destination == "external-final-sink"
            else TransferAuthorityClass.INTER_TRAY
        )
        if (
            child_route.authority_class is not expected_authority
            or child_route.topology.source_cell_id != lineage.source_tray_id
            or child_route.topology.destination_cell_id != child_destination
        ):
            raise KCellTowerRouterConservationError(
                "split child route topology differs from certified ownership"
            )
        parent_totals = _deterministic_split_totals(parent_entry)
        remainder_totals = _deterministic_split_totals(remainder_entry)
        child_totals = _deterministic_split_totals(child_entry)
        routed_child_totals = _deterministic_split_totals(routed_child_entry)
        if child_totals != routed_child_totals:
            raise KCellTowerRouterConservationError(
                "owned and routed split-child deterministic totals differ"
            )
        for index, (parent, remainder, child) in enumerate(
            zip(parent_totals, remainder_totals, child_totals, strict=True)
        ):
            correction[index] += (
                Fraction.from_float(remainder)
                + Fraction.from_float(child)
                - Fraction.from_float(parent)
            )
    if observed_counts != declared_counts:
        raise KCellTowerRouterConservationError(
            "split lineage count differs from per-layer exit-ledger certification"
        )
    return tuple(correction)


def _compute_whole_router_ledger(
    *,
    attempt_id: str,
    prior_trays: tuple[AcceptedKCellTrayState, ...],
    feed_entries: tuple[PacketPlanEntry, ...],
    candidate_trays: tuple[AcceptedKCellTrayState, ...],
    final_sink_entries: tuple[PacketPlanEntry, ...],
    routes: tuple[AcceptedPopulationRoute, ...],
    exit_ledgers: tuple[TrayExitLedger, ...],
    split_lineage: tuple[PacketSplitLineage, ...],
) -> WholeRouterLedger:
    before_entries = tuple(packet.entry for tray in prior_trays for packet in tray.packets)
    after_entries = tuple(packet.entry for tray in candidate_trays for packet in tray.packets)
    before_totals = _sum_entries(before_entries)
    feed_totals = _sum_entries(feed_entries) if feed_entries else _zero_extensives()
    after_totals = _sum_entries(after_entries)
    sink_totals = _sum_entries(final_sink_entries) if final_sink_entries else _zero_extensives()
    raw_exact_residual = _exact_signed_entry_residual(
        credits=(before_entries, feed_entries),
        debits=(after_entries, final_sink_entries),
    )
    split_correction = _certified_split_rounding_correction(
        attempt_id=attempt_id,
        prior_trays=prior_trays,
        candidate_trays=candidate_trays,
        final_sink_entries=final_sink_entries,
        routes=routes,
        exit_ledgers=exit_ledgers,
        split_lineage=split_lineage,
    )
    authoritative_residual = tuple(
        raw + correction
        for raw, correction in zip(raw_exact_residual, split_correction, strict=True)
    )
    binary64_residual = _binary64_regrouping_residual(
        credits=(before_totals, feed_totals),
        debits=(after_totals, sink_totals),
    )
    return WholeRouterLedger(
        audit_law_id=ROUTER_CONSERVATION_AUDIT_LAW_ID,
        tray_before_totals=before_totals,
        feed_totals=feed_totals,
        tray_after_totals=after_totals,
        final_sink_totals=sink_totals,
        raw_exact_signed_entry_residual=raw_exact_residual,
        certified_split_rounding_correction=split_correction,
        conservation_residual=authoritative_residual,
        binary64_regrouping_residual=binary64_residual,
        route_count=len(routes),
        passed=all(value == 0 for value in authoritative_residual),
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class TowerRoutingRejection:
    attempt_id: str
    reason: str
    rollback_trays: tuple[AcceptedKCellTrayState, ...]

    physically_qualifying: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedSerialTowerRouting:
    request: SerialTowerRoutingRequest
    candidate_trays: tuple[AcceptedKCellTrayState, ...]
    exit_ledgers: tuple[TrayExitLedger, ...]
    routes: tuple[AcceptedPopulationRoute, ...]
    split_lineage: tuple[PacketSplitLineage, ...]
    final_sink_packet_plan: tuple[PacketPlanEntry, ...]
    whole_router_ledger: WholeRouterLedger
    auditor_pin: tuple[str, str]
    planner_pins: tuple[tuple[str, str], ...]
    request_digest: str
    transaction_digest: str
    replay_token: str
    _router_token: object
    _seal: object

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _VALIDATED_SEAL:
            raise TypeError("validated tower routes are issued only by the router")


@dataclass(frozen=True, slots=True)
class _ExitResult:
    remaining_packets: tuple[AcceptedTrayPacket, ...]
    outgoing_entries: tuple[PacketPlanEntry, ...]
    ledger: TrayExitLedger
    split_lineage: tuple[PacketSplitLineage, ...]


def _neighboring_floats(center: float, *, limit: int = 128) -> tuple[float, ...]:
    values = [center]
    lower = center
    upper = center
    for _ in range(limit):
        lower = math.nextafter(lower, -math.inf)
        upper = math.nextafter(upper, math.inf)
        values.extend((lower, upper))
    return tuple(values)


def _closing_share(
    *,
    target_total: float,
    other_shares: tuple[float, ...],
    upper_bound: float,
    exact_offset_from: float | None = None,
    required_exact_offset: Fraction | None = None,
) -> float:
    if (exact_offset_from is None) != (required_exact_offset is None):
        raise KCellTowerRouterConfigurationError(
            "closing-share exact-exchange operands must be supplied together"
        )
    if exact_offset_from is not None:
        _require_binary64("closing-share exact-exchange baseline", exact_offset_from)
        if type(required_exact_offset) is not Fraction:
            raise KCellTowerRouterConfigurationError(
                "closing-share required exact exchange must be a Fraction"
            )
    center = math.fsum((target_total, -math.fsum(other_shares)))
    for candidate in _neighboring_floats(center):
        if candidate < 0.0 or candidate > upper_bound or not math.isfinite(candidate):
            continue
        if math.fsum((*other_shares, candidate)) != target_total:
            continue
        if exact_offset_from is not None and (
            Fraction.from_float(candidate) - Fraction.from_float(exact_offset_from)
            != required_exact_offset
        ):
            continue
        return 0.0 if candidate == 0.0 else candidate
    message = (
        "accepted K-share closing residual is not exactly representable"
        if exact_offset_from is None
        else (
            "accepted K-share closing residual and required exact exchange are not "
            "jointly representable"
        )
    )
    raise KCellTowerRouterConservationError(message)


def _ruled_initial_closing_baseline(
    *,
    accepted: float,
    ideals: dict[int, float],
    eligible_mass: dict[int, float],
    positive_layers: tuple[int, ...],
    closing_layer: int,
    k_count: int,
    policy: "KShareClosingExchangePolicy | None",
) -> tuple[float, int | None, float | None, int, int | None]:
    """Legacy closing baseline with the ruled initial-closing extension.

    GT_PS2_V3_EXTENSION_AND_RPM_AUTHORITY_RULING_2026-08-30 Ruling A: when
    binary64 cannot close the ideal split at all and an owner-issued
    closing-exchange policy is present, displace exactly ONE nonclosing
    ideal by one ULP and re-derive the closing share.  Candidates are
    enumerated absorber-cyclically (layers after the closing layer first,
    wrapping, never the closing layer), direction-minor with the downward
    ULP before the upward; the encounter rank is the enumeration index.
    Selection follows the ruled order law: nearest-to-ideal exact-rational
    distance of the resulting closing value, ties by encounter rank.  The
    selection is pure arithmetic so the coherence validator re-derives it
    independently.  Without a policy the sealed refusal re-raises
    byte-identically.

    Returns ``(closing_value, absorber_layer, displaced_value,
    signed_ulp_displacement, encounter_rank)``; the sealed path returns
    ``(value, None, None, 0, None)``.
    """

    other_shares = tuple(ideals[layer] for layer in range(1, k_count + 1) if layer != closing_layer)
    try:
        value = _closing_share(
            target_total=accepted,
            other_shares=other_shares,
            upper_bound=eligible_mass[closing_layer],
        )
    except KCellTowerRouterConservationError:
        if policy is None:
            raise
        absorber_layers = tuple(
            layer for layer in positive_layers if layer > closing_layer
        ) + tuple(layer for layer in positive_layers if layer < closing_layer)
        candidates: list[tuple[float, int, float, int, int]] = []
        rank = -1
        for absorber in absorber_layers:
            for direction, signed_ulp in ((-math.inf, -1), (math.inf, 1)):
                rank += 1
                displaced = math.nextafter(ideals[absorber], direction)
                if (
                    displaced < 0.0
                    or displaced > eligible_mass[absorber]
                    or not math.isfinite(displaced)
                ):
                    continue
                displaced_shares = tuple(
                    (displaced if layer == absorber else ideals[layer])
                    for layer in range(1, k_count + 1)
                    if layer != closing_layer
                )
                try:
                    closing_value = _closing_share(
                        target_total=accepted,
                        other_shares=displaced_shares,
                        upper_bound=eligible_mass[closing_layer],
                    )
                except KCellTowerRouterConservationError:
                    continue
                candidates.append((closing_value, absorber, displaced, signed_ulp, rank))
        if not candidates:
            raise
        ideal_closing_exact = Fraction.from_float(ideals[closing_layer])
        selected = min(
            candidates,
            key=lambda item: (
                abs(Fraction.from_float(item[0]) - ideal_closing_exact),
                item[4],
            ),
        )
        return (selected[0], selected[1], selected[2], selected[3], selected[4])
    return (value, None, None, 0, None)


def _split_child_id(*, attempt_id: str, tray_id: str, packet_id: str) -> str:
    identity = _digest(
        b"boundary-split-child",
        _text_part(attempt_id),
        _text_part(tray_id),
        _text_part(packet_id),
    )
    return f"{packet_id}::transfer::{identity[7:31]}"


def _split_candidate(
    *,
    entry: PacketPlanEntry,
    child_id: str,
    child_weight_value: float,
    parent_weight: float,
    tray_id: str,
    layer_id: VerticalLayerId,
) -> tuple[PacketPlanEntry, PacketPlanEntry, PacketSplitLineage, tuple[float, ...]] | None:
    if not 0.0 < child_weight_value < parent_weight:
        return None
    remainder_weight_value = parent_weight - child_weight_value
    if not 0.0 < remainder_weight_value < parent_weight:
        return None
    if math.fsum((child_weight_value, remainder_weight_value)) != parent_weight:
        return None
    child_weight = PacketWeight(
        representative_particles=child_weight_value,
        energy_datum_id=entry.weight.energy_datum_id,
    )
    remainder_weight = PacketWeight(
        representative_particles=remainder_weight_value,
        energy_datum_id=entry.weight.energy_datum_id,
    )
    child = PacketPlanEntry(
        packet_id=child_id,
        payload_bytes=entry.payload_bytes,
        inventory=entry.inventory,
        weight=child_weight,
        surface=entry.surface,
    )
    remainder = PacketPlanEntry(
        packet_id=entry.packet_id,
        payload_bytes=entry.payload_bytes,
        inventory=entry.inventory,
        weight=remainder_weight,
        surface=entry.surface,
    )
    child_totals = _deterministic_split_totals(child)
    _deterministic_split_totals(remainder)
    lineage = PacketSplitLineage(
        parent_packet_id=entry.packet_id,
        source_remainder_packet_id=entry.packet_id,
        transfer_child_packet_id=child_id,
        source_tray_id=tray_id,
        source_layer_id=layer_id,
        payload_identity=_payload_identity(entry.payload_bytes),
        parent_weight=parent_weight,
        remainder_weight=remainder_weight_value,
        transfer_weight=child_weight_value,
    )
    return remainder, child, lineage, child_totals


def _split_entry_exact(
    *,
    entry: PacketPlanEntry,
    requested_child_dry_kg: float,
    prior_layer_outgoing_entries: tuple[PacketPlanEntry, ...],
    represented_layer_target_dry_kg: float,
    attempt_id: str,
    tray_id: str,
    layer_id: VerticalLayerId,
) -> tuple[PacketPlanEntry, PacketPlanEntry, PacketSplitLineage, int]:
    _deterministic_split_totals(entry)
    per_particle_dry = entry.inventory.dry_matter_kg
    if per_particle_dry <= 0.0:
        raise KCellTowerRouterConservationError(
            "a positive boundary split requires positive per-particle dry matter"
        )
    replayed_remaining = represented_layer_target_dry_kg
    for prior_entry in prior_layer_outgoing_entries:
        replayed_remaining = math.fsum(
            (replayed_remaining, -_deterministic_split_totals(prior_entry)[0])
        )
    if replayed_remaining != requested_child_dry_kg:
        raise KCellTowerRouterConservationError(
            "boundary split search center differs from the canonical K-share subtraction tree"
        )
    parent_weight = entry.weight.representative_particles
    center = requested_child_dry_kg / per_particle_dry
    child_id = _split_child_id(
        attempt_id=attempt_id,
        tray_id=tray_id,
        packet_id=entry.packet_id,
    )
    candidates = _neighboring_floats(center)
    lower_target = math.nextafter(represented_layer_target_dry_kg, -math.inf)
    upper_target = math.nextafter(represented_layer_target_dry_kg, math.inf)
    fallback_candidates: list[tuple[float, int]] = []
    # Preserve the bit identity of every existing exact-child construction.
    # When the sequential subtraction remainder is outside the image of the
    # binary64 inventory-times-weight product, the ruled construction authority
    # is instead the actual encoded K-leaf group closing its represented share.
    for require_exact_child_product in (True, False):
        for candidate_rank, child_weight_value in enumerate(candidates):
            candidate = _split_candidate(
                entry=entry,
                child_id=child_id,
                child_weight_value=child_weight_value,
                parent_weight=parent_weight,
                tray_id=tray_id,
                layer_id=layer_id,
            )
            if candidate is None:
                continue
            remainder, child, lineage, child_totals = candidate
            represented_group = _sum_entries((*prior_layer_outgoing_entries, child))[0]
            if (
                not require_exact_child_product
                and represented_group in (lower_target, upper_target)
                and all(target != represented_group for target, _ in fallback_candidates)
            ):
                fallback_candidates.append((represented_group, candidate_rank))
            if require_exact_child_product and child_totals[0] != requested_child_dry_kg:
                continue
            if represented_group != represented_layer_target_dry_kg:
                continue
            # The complete K-group dry debit and external weight partition are
            # the exact routing constraints.  The sequential subtraction
            # remainder is only the center of the canonical weight search.
            # Seven-component product differences remain derived by the
            # lineage-certified whole-router split correction.
            return remainder, child, lineage, candidate_rank

    if fallback_candidates:
        raise _KShareRepresentabilityExhaustion(
            layer_id=layer_id,
            represented_target_dry_kg=represented_layer_target_dry_kg,
            one_ulp_fallback_candidates=tuple(fallback_candidates),
        )
    raise _KShareConstructionExhaustion(
        layer_id=layer_id,
        represented_target_dry_kg=represented_layer_target_dry_kg,
    )


def _final_stage_count(tray: AcceptedKCellTrayState) -> int:
    return next(
        item.rtd_stage_count
        for item in tray.topology.trays
        if item.physical_tray_id == tray.physical_tray_id
    )


def _layer_count(tray: AcceptedKCellTrayState) -> int:
    return next(
        item.vertical_layer_count
        for item in tray.topology.trays
        if item.physical_tray_id == tray.physical_tray_id
    )


def _validate_k_share_apportionment_policy(
    ledgers: tuple[TrayExitLedger, ...],
    *,
    closing_exchange_policy: "KShareClosingExchangePolicy | None" = None,
) -> None:
    """Recompute the v2 share vector and its explicit transaction record.

    ``closing_exchange_policy`` gates the ruled initial-closing extension
    (2026-08-30 Ruling A): without it, an activated initial-closing record
    fails closed exactly like every other unexplained construction.
    """

    for ledger in ledgers:
        if type(ledger) is not TrayExitLedger:
            raise KCellTowerRouterConservationError(
                "K-share apportionment validation requires exact tray exit ledgers"
            )
        if ledger.apportionment_policy_id != QSC_K_SHARE_APPORTIONMENT_POLICY_ID:
            raise KCellTowerRouterConservationError(
                "tray exit ledger carries a foreign K-share apportionment policy"
            )
        if type(ledger.apportionment_disposition) is not KShareApportionmentDisposition:
            raise KCellTowerRouterConservationError(
                "tray exit ledger carries a foreign K-share apportionment disposition"
            )
        if (
            ledger.adjusted_nonclosing_layer_id is not None
            and type(ledger.adjusted_nonclosing_layer_id) is not VerticalLayerId
        ):
            raise KCellTowerRouterConservationError(
                "tray exit ledger carries a foreign adjusted K-layer identity"
            )
        if (
            ledger.closing_layer_id is not None
            and type(ledger.closing_layer_id) is not VerticalLayerId
        ):
            raise KCellTowerRouterConservationError(
                "tray exit ledger carries a foreign closing K-layer identity"
            )
        if type(ledger.selected_nonclosing_signed_ulp_displacement) is not int:
            raise KCellTowerRouterConservationError(
                "tray exit ledger carries a non-integer selected K-share displacement"
            )
        if ledger.inherited_direction_candidate_rank is not None and (
            type(ledger.inherited_direction_candidate_rank) is not int
            or ledger.inherited_direction_candidate_rank < 0
        ):
            raise KCellTowerRouterConservationError(
                "tray exit ledger carries a foreign inherited candidate rank"
            )
        _require_nonblank("tray exit ledger physical ID", ledger.physical_tray_id)
        _require_binary64(
            "tray exit ledger nonnegative operand",
            ledger.accepted_dry_matter_debit_kg,
            ledger.represented_dry_matter_debit_kg,
            ledger.eligible_dry_matter_kg,
            nonnegative=True,
        )
        if ledger.passed is not True:
            raise KCellTowerRouterConservationError(
                "tray exit ledger does not carry an exact passed selection record"
            )
        if (
            type(ledger.requested_closing_layer_cursor) is not int
            or ledger.requested_closing_layer_cursor < 1
            or type(ledger.next_closing_layer_cursor) is not int
            or ledger.next_closing_layer_cursor < 1
        ):
            raise KCellTowerRouterConservationError(
                "K-share apportionment cursors must be positive exact integers"
            )
        if type(ledger.k_shares) is not tuple or not ledger.k_shares:
            raise KCellTowerRouterConservationError(
                "K-share apportionment requires an exact nonempty K-row tuple"
            )
        if any(type(share) is not KShareExitLedger for share in ledger.k_shares):
            raise KCellTowerRouterConservationError(
                "K-share apportionment validation requires exact K-share ledgers"
            )
        if any(type(share.layer_id) is not VerticalLayerId for share in ledger.k_shares):
            raise KCellTowerRouterConservationError(
                "K-share apportionment rows carry foreign K-layer identities"
            )
        expected_layers = tuple(range(1, len(ledger.k_shares) + 1))
        if tuple(share.layer_id.value for share in ledger.k_shares) != expected_layers:
            raise KCellTowerRouterConservationError(
                "K-share apportionment rows differ from canonical unique K order"
            )
        for share in ledger.k_shares:
            _require_binary64(
                "K-share nonnegative operand",
                share.eligible_dry_matter_kg,
                share.ideal_dry_matter_debit_kg,
                share.legacy_represented_dry_matter_debit_kg,
                share.represented_dry_matter_debit_kg,
                nonnegative=True,
            )
            _require_binary64(
                "K-share represented-minus-ideal diagnostic",
                share.represented_minus_ideal_kg,
            )
            expected_discrepancy = math.fsum(
                (
                    share.represented_dry_matter_debit_kg,
                    -share.ideal_dry_matter_debit_kg,
                )
            )
            if not _same_float(share.represented_minus_ideal_kg, expected_discrepancy):
                raise KCellTowerRouterConservationError(
                    "K-share represented-minus-ideal diagnostic is not derived"
                )
            if type(
                share.signed_binary64_steps_from_legacy
            ) is not int or share.signed_binary64_steps_from_legacy != _signed_nonnegative_ulp_distance(
                share.represented_dry_matter_debit_kg,
                share.legacy_represented_dry_matter_debit_kg,
            ):
                raise KCellTowerRouterConservationError(
                    "K-share signed binary64 displacement is not derived"
                )
            if type(share.apportionment_role) is not KShareApportionmentRole:
                raise KCellTowerRouterConservationError(
                    "K-share row carries a foreign apportionment role"
                )
            if (
                type(share.whole_packet_count) is not int
                or share.whole_packet_count < 0
                or type(share.boundary_split_count) is not int
                or share.boundary_split_count not in (0, 1)
            ):
                raise KCellTowerRouterConservationError(
                    "K-share packet counts are outside the one-split construction"
                )
            if share.boundary_split_count == 0:
                if share.selected_boundary_weight_candidate_rank is not None:
                    raise KCellTowerRouterConservationError(
                        "unsplit K share declares a boundary-weight candidate rank"
                    )
            elif (
                type(share.selected_boundary_weight_candidate_rank) is not int
                or share.selected_boundary_weight_candidate_rank < 0
            ):
                raise KCellTowerRouterConservationError(
                    "split K share omits its selected boundary-weight candidate rank"
                )
        represented_total = math.fsum(
            share.represented_dry_matter_debit_kg for share in ledger.k_shares
        )
        if represented_total != ledger.accepted_dry_matter_debit_kg:
            raise KCellTowerRouterConservationError(
                "K-share apportionment does not close the booked accepted debit"
            )
        if ledger.accepted_dry_matter_debit_kg == 0.0:
            if (
                ledger.apportionment_disposition is not KShareApportionmentDisposition.LEGACY_ZERO
                or ledger.closing_layer_id is not None
                or ledger.adjusted_nonclosing_layer_id is not None
                or ledger.selected_nonclosing_signed_ulp_displacement != 0
                or ledger.inherited_direction_candidate_rank is not None
                or ledger.selected_closing_signed_ulp_displacement != 0
                or ledger.inherited_closing_candidate_rank is not None
                or ledger.representation_fallback_applied
                or ledger.next_closing_layer_cursor != ledger.requested_closing_layer_cursor
            ):
                raise KCellTowerRouterConservationError(
                    "zero-debit K-share apportionment carries nonlegacy selection state"
                )
            if any(
                share.ideal_dry_matter_debit_kg != 0.0
                or share.legacy_represented_dry_matter_debit_kg != 0.0
                or share.represented_dry_matter_debit_kg != 0.0
                or share.exact_represented_minus_ideal_kg != 0
                or share.signed_binary64_steps_from_legacy != 0
                or share.apportionment_role is not KShareApportionmentRole.UNCHANGED
                for share in ledger.k_shares
            ):
                raise KCellTowerRouterConservationError(
                    "zero-debit K-share apportionment changes a zero target"
                )
            continue
        if ledger.closing_layer_id is None:
            raise KCellTowerRouterConservationError(
                "positive K-share apportionment omits its declared closing layer"
            )
        closing_rows = tuple(
            share for share in ledger.k_shares if share.layer_id == ledger.closing_layer_id
        )
        if len(closing_rows) != 1:
            raise KCellTowerRouterConservationError(
                "K-share apportionment does not contain exactly one declared closing share"
            )
        closing = closing_rows[0]
        nonclosing = tuple(
            share for share in ledger.k_shares if share.layer_id != ledger.closing_layer_id
        )
        positive_layers = tuple(
            share.layer_id.value for share in ledger.k_shares if share.eligible_dry_matter_kg > 0.0
        )
        if not positive_layers:
            raise KCellTowerRouterConservationError(
                "positive K-share apportionment has no eligible K layer"
            )
        expected_closing_layer = next(
            (layer for layer in positive_layers if layer >= ledger.requested_closing_layer_cursor),
            positive_layers[0],
        )
        if ledger.closing_layer_id.value != expected_closing_layer:
            raise KCellTowerRouterConservationError(
                "K-share closing layer differs from the requested balanced cursor"
            )
        if ledger.next_closing_layer_cursor != expected_closing_layer % len(ledger.k_shares) + 1:
            raise KCellTowerRouterConservationError(
                "K-share next closing cursor is not the canonical successor"
            )
        (
            legacy_closing,
            derived_initial_absorber,
            derived_initial_absorber_value,
            derived_initial_signed_ulp,
            derived_initial_rank,
        ) = _ruled_initial_closing_baseline(
            accepted=ledger.accepted_dry_matter_debit_kg,
            ideals={
                share.layer_id.value: share.ideal_dry_matter_debit_kg for share in ledger.k_shares
            },
            eligible_mass={
                share.layer_id.value: share.eligible_dry_matter_kg for share in ledger.k_shares
            },
            positive_layers=positive_layers,
            closing_layer=ledger.closing_layer_id.value,
            k_count=len(ledger.k_shares),
            policy=closing_exchange_policy,
        )
        initial_extension_active = derived_initial_absorber is not None
        for share in ledger.k_shares:
            expected_legacy = (
                legacy_closing
                if share.layer_id == ledger.closing_layer_id
                else share.ideal_dry_matter_debit_kg
            )
            if not _same_float(share.legacy_represented_dry_matter_debit_kg, expected_legacy):
                raise KCellTowerRouterConservationError(
                    "K-share legacy represented target is not the canonical construction"
                )
        fallback_layer_id = ledger.representation_fallback_layer_id
        if fallback_layer_id is None:
            if (
                ledger.apportionment_disposition is not KShareApportionmentDisposition.LEGACY_ZERO
                or ledger.adjusted_nonclosing_layer_id is not None
                or ledger.selected_nonclosing_signed_ulp_displacement != 0
                or ledger.inherited_direction_candidate_rank is not None
                or ledger.selected_closing_signed_ulp_displacement != 0
                or ledger.inherited_closing_candidate_rank is not None
                or any(
                    not _same_float(
                        share.represented_dry_matter_debit_kg,
                        share.legacy_represented_dry_matter_debit_kg,
                    )
                    or share.signed_binary64_steps_from_legacy != 0
                    or share.apportionment_role
                    is not (
                        KShareApportionmentRole.DECLARED_CLOSING
                        if share.layer_id == ledger.closing_layer_id
                        else KShareApportionmentRole.UNCHANGED
                    )
                    for share in ledger.k_shares
                )
            ):
                raise KCellTowerRouterConservationError(
                    "legacy K-share record differs from its canonical zero-displacement result"
                )
            continue
        fallback = next(share for share in nonclosing if share.layer_id == fallback_layer_id)
        if any(
            share.layer_id != fallback_layer_id
            and not _same_float(
                share.represented_dry_matter_debit_kg,
                share.legacy_represented_dry_matter_debit_kg,
            )
            for share in nonclosing
        ):
            raise KCellTowerRouterConservationError(
                "more than one nonclosing K share changed under the fallback policy"
            )
        if any(
            share.layer_id != fallback_layer_id
            and (
                share.signed_binary64_steps_from_legacy != 0
                or share.apportionment_role is not KShareApportionmentRole.UNCHANGED
            )
            for share in nonclosing
        ):
            raise KCellTowerRouterConservationError(
                "unchanged nonclosing K-share record carries fallback evidence"
            )
        if (
            fallback.apportionment_role is not KShareApportionmentRole.ADJUSTED_NONCLOSING
            or closing.apportionment_role is not KShareApportionmentRole.DECLARED_CLOSING
        ):
            raise KCellTowerRouterConservationError(
                "fallback K-share roles or signed displacement contradict the tray record"
            )
        exact_nonclosing_offset = Fraction.from_float(
            fallback.represented_dry_matter_debit_kg
        ) - Fraction.from_float(fallback.legacy_represented_dry_matter_debit_kg)
        exact_closing_offset = Fraction.from_float(
            closing.represented_dry_matter_debit_kg
        ) - Fraction.from_float(legacy_closing)
        if initial_extension_active:
            # Ruled initial-closing extension (2026-08-30 Ruling A): the
            # baseline itself carries the one-ULP displacement, so the
            # zero-sum exchange law is replaced by exact equality with the
            # independent re-derivation — a forged absorber, direction,
            # rank, or share cannot validate.
            if (
                ledger.apportionment_disposition
                is not KShareApportionmentDisposition.ONE_NONCLOSING_ONE_ULP
                or ledger.adjusted_nonclosing_layer_id is None
                or ledger.adjusted_nonclosing_layer_id.value != derived_initial_absorber
                or ledger.selected_nonclosing_signed_ulp_displacement != derived_initial_signed_ulp
                or ledger.inherited_direction_candidate_rank != derived_initial_rank
                or ledger.selected_closing_signed_ulp_displacement != 0
                or ledger.inherited_closing_candidate_rank != derived_initial_rank
                or not _same_float(
                    fallback.represented_dry_matter_debit_kg,
                    derived_initial_absorber_value,
                )
                or not _same_float(closing.represented_dry_matter_debit_kg, legacy_closing)
                or exact_closing_offset != 0
            ):
                raise KCellTowerRouterConservationError(
                    "initial-closing extension record differs from its ruled re-derivation"
                )
            continue
        if exact_nonclosing_offset + exact_closing_offset != 0:
            raise KCellTowerRouterConservationError(
                "fallback nonclosing/closing K-share exchange is not exact zero-sum"
            )
        if (
            ledger.apportionment_disposition
            is KShareApportionmentDisposition.ONE_NONCLOSING_ONE_ULP
        ):
            if (
                ledger.adjusted_nonclosing_layer_id != fallback_layer_id
                or ledger.selected_nonclosing_signed_ulp_displacement not in (-1, 1)
                or type(ledger.inherited_direction_candidate_rank) is not int
                or ledger.inherited_direction_candidate_rank < 0
                or ledger.selected_closing_signed_ulp_displacement != 0
                or ledger.inherited_closing_candidate_rank is not None
            ):
                raise KCellTowerRouterConservationError(
                    "fallback tray record omits its bounded v2 selection evidence"
                )
            lower = math.nextafter(fallback.ideal_dry_matter_debit_kg, -math.inf)
            upper = math.nextafter(fallback.ideal_dry_matter_debit_kg, math.inf)
            if fallback.represented_dry_matter_debit_kg not in (lower, upper):
                raise KCellTowerRouterConservationError(
                    "nonclosing K-share fallback exceeds the approved one-ULP domain"
                )
            if (
                fallback.signed_binary64_steps_from_legacy
                != ledger.selected_nonclosing_signed_ulp_displacement
            ):
                raise KCellTowerRouterConservationError(
                    "fallback K-share roles or signed displacement contradict the tray record"
                )
            expected_closing = _closing_share(
                target_total=ledger.accepted_dry_matter_debit_kg,
                other_shares=tuple(share.represented_dry_matter_debit_kg for share in nonclosing),
                upper_bound=closing.eligible_dry_matter_kg,
                exact_offset_from=legacy_closing,
                required_exact_offset=-exact_nonclosing_offset,
            )
            if not _same_float(closing.represented_dry_matter_debit_kg, expected_closing):
                raise KCellTowerRouterConservationError(
                    "fallback did not recompute only the declared closing K share"
                )
        elif (
            ledger.apportionment_disposition
            is KShareApportionmentDisposition.CLOSING_ONE_ULP_EXCHANGE
        ):
            # v3 closing exchange
            # (GT_PS2_WALLED_FLOOR_AND_V3_APPROVAL_RULING_2026-08-29).
            if (
                ledger.adjusted_nonclosing_layer_id != fallback_layer_id
                or type(ledger.selected_nonclosing_signed_ulp_displacement) is not int
                or ledger.selected_nonclosing_signed_ulp_displacement == 0
                or ledger.inherited_direction_candidate_rank is not None
                or ledger.selected_closing_signed_ulp_displacement not in (-1, 1)
                or type(ledger.inherited_closing_candidate_rank) is not int
                or ledger.inherited_closing_candidate_rank < 0
            ):
                raise KCellTowerRouterConservationError(
                    "closing-exchange tray record omits its bounded v3 selection evidence"
                )
            closing_lower = math.nextafter(legacy_closing, -math.inf)
            closing_upper = math.nextafter(legacy_closing, math.inf)
            if closing.represented_dry_matter_debit_kg not in (closing_lower, closing_upper):
                raise KCellTowerRouterConservationError(
                    "closing K-share exchange exceeds the approved one-ULP domain"
                )
            if (
                closing.signed_binary64_steps_from_legacy
                != ledger.selected_closing_signed_ulp_displacement
                or fallback.signed_binary64_steps_from_legacy
                != ledger.selected_nonclosing_signed_ulp_displacement
            ):
                raise KCellTowerRouterConservationError(
                    "closing-exchange K-share displacements contradict the tray record"
                )
            expected_absorber = _closing_share(
                target_total=ledger.accepted_dry_matter_debit_kg,
                other_shares=tuple(
                    share.represented_dry_matter_debit_kg
                    for share in ledger.k_shares
                    if share.layer_id != fallback_layer_id
                ),
                upper_bound=fallback.eligible_dry_matter_kg,
                exact_offset_from=(fallback.legacy_represented_dry_matter_debit_kg),
                required_exact_offset=-exact_closing_offset,
            )
            if not _same_float(fallback.represented_dry_matter_debit_kg, expected_absorber):
                raise KCellTowerRouterConservationError(
                    "closing exchange did not recompute only the single absorber K share"
                )
        else:
            raise KCellTowerRouterConservationError(
                "fallback tray record names an unknown apportionment disposition"
            )


def _construct_tray_exit(
    *,
    request: SerialTowerRoutingRequest,
    tray: AcceptedKCellTrayState,
    k_count: int,
    eligible: dict[int, tuple[AcceptedTrayPacket, ...]],
    eligible_mass: dict[int, float],
    total_eligible: float,
    accepted: float,
    closing_layer: int,
    requested_closing_layer_cursor: int,
    ideals: dict[int, float],
    legacy_represented: dict[int, float],
    represented: dict[int, float],
    apportionment_disposition: KShareApportionmentDisposition,
    adjusted_nonclosing_layer: int | None,
    selected_nonclosing_signed_ulp_displacement: int,
    inherited_direction_candidate_rank: int | None,
    selected_closing_signed_ulp_displacement: int = 0,
    inherited_closing_candidate_rank: int | None = None,
) -> _ExitResult:
    """Build one immutable exit candidate from the untouched accepted tray."""

    by_id = {packet.packet_id: packet for packet in tray.packets}
    outgoing: list[PacketPlanEntry] = []
    lineage: list[PacketSplitLineage] = []
    share_ledgers: list[KShareExitLedger] = []
    for layer in range(1, k_count + 1):
        target = represented[layer]
        remaining = target
        layer_outgoing: list[PacketPlanEntry] = []
        whole_count = 0
        split_count = 0
        selected_boundary_weight_candidate_rank: int | None = None
        for packet in eligible[layer]:
            if remaining == 0.0:
                break
            packet_dry = packet.weighted_totals.totals[0]
            if packet_dry <= remaining:
                layer_outgoing.append(packet.entry)
                del by_id[packet.packet_id]
                whole_count += 1
                remaining = math.fsum((remaining, -packet_dry))
                if remaining == 0.0:
                    remaining = 0.0
                continue
            (
                remainder_entry,
                child_entry,
                split,
                selected_boundary_weight_candidate_rank,
            ) = _split_entry_exact(
                entry=packet.entry,
                requested_child_dry_kg=remaining,
                prior_layer_outgoing_entries=tuple(layer_outgoing),
                represented_layer_target_dry_kg=target,
                attempt_id=request.attempt_id,
                tray_id=tray.physical_tray_id,
                layer_id=VerticalLayerId(layer),
            )
            by_id[packet.packet_id] = replace(packet, entry=remainder_entry)
            layer_outgoing.append(child_entry)
            lineage.append(split)
            split_count = 1
            remaining = 0.0
            break
        if remaining != 0.0:
            raise _KShareConstructionExhaustion(
                layer_id=VerticalLayerId(layer),
                represented_target_dry_kg=target,
            )
        represented_total = _sum_entries(tuple(layer_outgoing))[0]
        if represented_total != target:
            raise KCellTowerRouterConservationError(
                f"{tray.physical_tray_id}:K{layer} exit differs from represented share"
            )
        outgoing.extend(layer_outgoing)
        share_ledgers.append(
            KShareExitLedger(
                layer_id=VerticalLayerId(layer),
                eligible_dry_matter_kg=eligible_mass[layer],
                ideal_dry_matter_debit_kg=ideals[layer],
                legacy_represented_dry_matter_debit_kg=legacy_represented[layer],
                represented_dry_matter_debit_kg=represented_total,
                represented_minus_ideal_kg=math.fsum((represented_total, -ideals[layer])),
                signed_binary64_steps_from_legacy=(
                    _signed_nonnegative_ulp_distance(
                        represented_total,
                        legacy_represented[layer],
                    )
                ),
                apportionment_role=(
                    KShareApportionmentRole.DECLARED_CLOSING
                    if layer == closing_layer
                    else (
                        KShareApportionmentRole.ADJUSTED_NONCLOSING
                        if layer == adjusted_nonclosing_layer
                        else KShareApportionmentRole.UNCHANGED
                    )
                ),
                selected_boundary_weight_candidate_rank=(selected_boundary_weight_candidate_rank),
                whole_packet_count=whole_count,
                boundary_split_count=split_count,
            )
        )
    if math.fsum(share.represented_dry_matter_debit_kg for share in share_ledgers) != accepted:
        raise KCellTowerRouterConservationError(
            f"{tray.physical_tray_id} grouped K shares differ from its booked accepted debit"
        )
    outgoing_tuple = tuple(sorted(outgoing, key=lambda entry: entry.packet_id))
    represented_total = _sum_entries(outgoing_tuple)[0]
    remaining_packets = tuple(sorted(by_id.values(), key=_packet_sort_key))
    if not remaining_packets:
        raise KCellTowerRouterConfigurationError(
            f"{tray.physical_tray_id} route would leave an empty accepted population"
        )
    next_cursor = closing_layer % k_count + 1
    ledger = TrayExitLedger(
        physical_tray_id=tray.physical_tray_id,
        apportionment_policy_id=QSC_K_SHARE_APPORTIONMENT_POLICY_ID,
        apportionment_disposition=apportionment_disposition,
        adjusted_nonclosing_layer_id=(
            VerticalLayerId(adjusted_nonclosing_layer)
            if adjusted_nonclosing_layer is not None
            else None
        ),
        selected_nonclosing_signed_ulp_displacement=(selected_nonclosing_signed_ulp_displacement),
        inherited_direction_candidate_rank=inherited_direction_candidate_rank,
        selected_closing_signed_ulp_displacement=(selected_closing_signed_ulp_displacement),
        inherited_closing_candidate_rank=inherited_closing_candidate_rank,
        accepted_dry_matter_debit_kg=accepted,
        represented_dry_matter_debit_kg=represented_total,
        eligible_dry_matter_kg=total_eligible,
        closing_layer_id=VerticalLayerId(closing_layer),
        requested_closing_layer_cursor=requested_closing_layer_cursor,
        next_closing_layer_cursor=next_cursor,
        k_shares=tuple(share_ledgers),
        passed=True,
    )
    return _ExitResult(
        remaining_packets=remaining_packets,
        outgoing_entries=outgoing_tuple,
        ledger=ledger,
        split_lineage=tuple(lineage),
    )


def _legacy_representability_refusal(message: str) -> KCellTowerRouterConservationError:
    return KCellTowerRouterConservationError(
        "boundary split cannot represent the exact K-share construction with exact "
        f"parent-weight closure: {message}"
    )


def _exit_from_tray(
    *,
    request: SerialTowerRoutingRequest,
    tray: AcceptedKCellTrayState,
    debit: AcceptedTrayExitDebit,
) -> _ExitResult:
    k_count = _layer_count(tray)
    final_stage = _final_stage_count(tray)
    eligible: dict[int, tuple[AcceptedTrayPacket, ...]] = {}
    eligible_mass: dict[int, float] = {}
    for layer in range(1, k_count + 1):
        packets = tuple(
            sorted(
                (
                    packet
                    for packet in tray.packets
                    if packet.owner_key.vertical_layer_id.value == layer
                    and packet.owner_key.rtd_stage_id.value == final_stage
                    and packet.weighted_totals.totals[0] > 0.0
                ),
                key=lambda packet: (
                    packet.rtd_path.transition_exposures[-1],
                    packet.packet_id,
                ),
            )
        )
        eligible[layer] = packets
        eligible_mass[layer] = math.fsum(packet.weighted_totals.totals[0] for packet in packets)
    total_eligible = math.fsum(eligible_mass.values())
    accepted = debit.accepted_dry_matter_kg
    if accepted > total_eligible:
        raise KCellTowerRouterConfigurationError(
            f"{tray.physical_tray_id} accepted debit exceeds final-stage eligible inventory"
        )
    if accepted == 0.0:
        shares = tuple(
            KShareExitLedger(
                layer_id=VerticalLayerId(layer),
                eligible_dry_matter_kg=eligible_mass[layer],
                ideal_dry_matter_debit_kg=0.0,
                legacy_represented_dry_matter_debit_kg=0.0,
                represented_dry_matter_debit_kg=0.0,
                represented_minus_ideal_kg=0.0,
                signed_binary64_steps_from_legacy=0,
                apportionment_role=KShareApportionmentRole.UNCHANGED,
                selected_boundary_weight_candidate_rank=None,
                whole_packet_count=0,
                boundary_split_count=0,
            )
            for layer in range(1, k_count + 1)
        )
        return _ExitResult(
            remaining_packets=tray.packets,
            outgoing_entries=(),
            ledger=TrayExitLedger(
                physical_tray_id=tray.physical_tray_id,
                apportionment_policy_id=QSC_K_SHARE_APPORTIONMENT_POLICY_ID,
                apportionment_disposition=(KShareApportionmentDisposition.LEGACY_ZERO),
                adjusted_nonclosing_layer_id=None,
                selected_nonclosing_signed_ulp_displacement=0,
                inherited_direction_candidate_rank=None,
                accepted_dry_matter_debit_kg=0.0,
                represented_dry_matter_debit_kg=0.0,
                eligible_dry_matter_kg=total_eligible,
                closing_layer_id=None,
                requested_closing_layer_cursor=debit.closing_layer_cursor,
                next_closing_layer_cursor=debit.closing_layer_cursor,
                k_shares=shares,
                passed=True,
            ),
            split_lineage=(),
        )
    positive_layers = tuple(layer for layer in range(1, k_count + 1) if eligible_mass[layer] > 0.0)
    if not positive_layers:
        raise KCellTowerRouterConfigurationError(
            f"{tray.physical_tray_id} has no positive final-stage eligible inventory"
        )
    closing_layer = next(
        (layer for layer in positive_layers if layer >= debit.closing_layer_cursor),
        positive_layers[0],
    )
    ideals = {
        layer: accepted * (eligible_mass[layer] / total_eligible) for layer in range(1, k_count + 1)
    }
    (
        initial_closing_value,
        initial_absorber_layer,
        initial_absorber_value,
        initial_absorber_signed_ulp,
        initial_closing_candidate_rank,
    ) = _ruled_initial_closing_baseline(
        accepted=accepted,
        ideals=ideals,
        eligible_mass=eligible_mass,
        positive_layers=positive_layers,
        closing_layer=closing_layer,
        k_count=k_count,
        policy=request.k_share_closing_exchange_policy,
    )
    represented = dict(ideals)
    represented[closing_layer] = initial_closing_value
    legacy_represented = dict(represented)
    if initial_absorber_layer is not None:
        represented[initial_absorber_layer] = initial_absorber_value

    def build(
        targets: dict[int, float],
        *,
        disposition: KShareApportionmentDisposition,
        adjusted_nonclosing_layer: int | None,
        selected_nonclosing_signed_ulp_displacement: int,
        inherited_direction_candidate_rank: int | None,
        selected_closing_signed_ulp_displacement: int = 0,
        inherited_closing_candidate_rank: int | None = None,
    ) -> _ExitResult:
        return _construct_tray_exit(
            request=request,
            tray=tray,
            k_count=k_count,
            eligible=eligible,
            eligible_mass=eligible_mass,
            total_eligible=total_eligible,
            accepted=accepted,
            closing_layer=closing_layer,
            requested_closing_layer_cursor=debit.closing_layer_cursor,
            ideals=ideals,
            legacy_represented=legacy_represented,
            represented=targets,
            apportionment_disposition=disposition,
            adjusted_nonclosing_layer=adjusted_nonclosing_layer,
            selected_nonclosing_signed_ulp_displacement=(
                selected_nonclosing_signed_ulp_displacement
            ),
            inherited_direction_candidate_rank=inherited_direction_candidate_rank,
            selected_closing_signed_ulp_displacement=(selected_closing_signed_ulp_displacement),
            inherited_closing_candidate_rank=inherited_closing_candidate_rank,
        )

    if initial_absorber_layer is not None:
        # Ruled initial-closing extension (2026-08-30 Ruling A): the legacy
        # split was not representable at all, so the policy-derived baseline
        # already carries one nonclosing one-ULP displacement.  Never nest a
        # second fallback: a build exhaustion here is a typed refusal.
        try:
            return build(
                represented,
                disposition=KShareApportionmentDisposition.ONE_NONCLOSING_ONE_ULP,
                adjusted_nonclosing_layer=initial_absorber_layer,
                selected_nonclosing_signed_ulp_displacement=initial_absorber_signed_ulp,
                inherited_direction_candidate_rank=initial_closing_candidate_rank,
                selected_closing_signed_ulp_displacement=0,
                inherited_closing_candidate_rank=initial_closing_candidate_rank,
            )
        except _KShareConstructionExhaustion as gated_failure:
            raise _legacy_representability_refusal(
                "the ruled initial-closing construction exhausted its selected candidate"
            ) from gated_failure

    # Candidate zero is the complete legacy plan.  It is built and returned
    # before the v2 policy can inspect any representability alternative.
    try:
        return build(
            represented,
            disposition=KShareApportionmentDisposition.LEGACY_ZERO,
            adjusted_nonclosing_layer=None,
            selected_nonclosing_signed_ulp_displacement=0,
            inherited_direction_candidate_rank=None,
        )
    except _KShareRepresentabilityExhaustion as caught_failure:
        legacy_failure = caught_failure
        failed_layer = caught_failure.layer_id.value
        if failed_layer == closing_layer:
            # v3 closing exchange, approved by
            # GT_PS2_WALLED_FLOOR_AND_V3_APPROVAL_RULING_2026-08-29: only the
            # declared closing layer's adjacent-candidate exhaustion, only
            # under an owner-issued policy; every other configuration keeps
            # the typed refusal below unchanged.
            policy = request.k_share_closing_exchange_policy
            if policy is None:
                raise _legacy_representability_refusal(
                    "the exhausted layer is the declared closing share"
                ) from caught_failure
            if not _same_float(
                caught_failure.represented_target_dry_kg,
                represented[closing_layer],
            ):
                raise _legacy_representability_refusal(
                    "the exhausted closing target differs from the legacy represented share"
                ) from caught_failure
            legacy_closing_target = represented[closing_layer]
            closing_lower = math.nextafter(legacy_closing_target, -math.inf)
            closing_upper = math.nextafter(legacy_closing_target, math.inf)
            ideal_closing_exact = Fraction.from_float(ideals[closing_layer])
            # Owner-selected order law (ruling question b): nearest-to-ideal
            # exact-rational distance, ties by recorded encounter rank.
            ordered_candidates = tuple(
                sorted(
                    legacy_failure.one_ulp_fallback_candidates,
                    key=lambda item: (
                        abs(Fraction.from_float(item[0]) - ideal_closing_exact),
                        item[1],
                    ),
                )
            )
            # First positive-eligible nonclosing layer strictly after the
            # resolved closing layer in cyclic order, wrapping once, never
            # revisiting the closing layer.
            absorber_layers = tuple(
                layer for layer in positive_layers if layer > closing_layer
            ) + tuple(layer for layer in positive_layers if layer < closing_layer)
            transcript: list[tuple[float, int, int | None, str]] = []
            for closing_target, encounter_rank in ordered_candidates:
                if closing_target not in (closing_lower, closing_upper):
                    transcript.append(
                        (closing_target, encounter_rank, None, "not one ULP adjacent")
                    )
                    continue
                if closing_target < 0.0 or closing_target > eligible_mass[closing_layer]:
                    transcript.append(
                        (closing_target, encounter_rank, None, "outside eligible inventory")
                    )
                    continue
                closing_offset = Fraction.from_float(closing_target) - Fraction.from_float(
                    legacy_closing_target
                )
                for absorber in absorber_layers:
                    adjusted = dict(represented)
                    adjusted[closing_layer] = closing_target
                    adjusted_other_shares = tuple(
                        adjusted[layer] for layer in range(1, k_count + 1) if layer != absorber
                    )
                    try:
                        absorber_target = _closing_share(
                            target_total=accepted,
                            other_shares=adjusted_other_shares,
                            upper_bound=eligible_mass[absorber],
                            exact_offset_from=represented[absorber],
                            required_exact_offset=-closing_offset,
                        )
                    except KCellTowerRouterConservationError:
                        transcript.append(
                            (
                                closing_target,
                                encounter_rank,
                                absorber,
                                "no exact representable complement",
                            )
                        )
                        continue
                    absorber_offset = Fraction.from_float(absorber_target) - Fraction.from_float(
                        represented[absorber]
                    )
                    if closing_offset + absorber_offset != 0:
                        transcript.append(
                            (closing_target, encounter_rank, absorber, "exchange not zero-sum")
                        )
                        continue
                    adjusted[absorber] = absorber_target
                    if math.fsum(adjusted[layer] for layer in range(1, k_count + 1)) != accepted:
                        transcript.append(
                            (
                                closing_target,
                                encounter_rank,
                                absorber,
                                "grouped shares do not close",
                            )
                        )
                        continue
                    try:
                        return build(
                            adjusted,
                            disposition=(KShareApportionmentDisposition.CLOSING_ONE_ULP_EXCHANGE),
                            adjusted_nonclosing_layer=absorber,
                            selected_nonclosing_signed_ulp_displacement=(
                                _signed_nonnegative_ulp_distance(
                                    absorber_target,
                                    represented[absorber],
                                )
                            ),
                            inherited_direction_candidate_rank=None,
                            selected_closing_signed_ulp_displacement=(
                                _signed_nonnegative_ulp_distance(
                                    closing_target,
                                    legacy_closing_target,
                                )
                            ),
                            inherited_closing_candidate_rank=encounter_rank,
                        )
                    except _KShareConstructionExhaustion as joint_failure:
                        joint_layer = joint_failure.layer_id.value
                        if joint_layer not in (closing_layer, absorber):
                            raise _legacy_representability_refusal(
                                "a v3 rebuild exhausted an unadjusted layer"
                            ) from joint_failure
                        # Never nest another fallback: discard this candidate.
                        transcript.append(
                            (closing_target, encounter_rank, absorber, "joint rebuild exhausted")
                        )
                        continue
            no_exchange = _legacy_representability_refusal(
                "no exact joint one-ULP closing/nonclosing exchange exists"
            )
            no_exchange.v3_search_transcript = tuple(transcript)
            raise no_exchange from legacy_failure
        if represented[failed_layer] != ideals[failed_layer]:
            raise _legacy_representability_refusal(
                "the exhausted nonclosing layer was not at its legacy ideal target"
            ) from caught_failure
        if not _same_float(
            caught_failure.represented_target_dry_kg,
            represented[failed_layer],
        ):
            raise _legacy_representability_refusal(
                "the exhausted target differs from the legacy represented share"
            ) from caught_failure

    legacy_closing_target = represented[closing_layer]
    lower_target = math.nextafter(represented[failed_layer], -math.inf)
    upper_target = math.nextafter(represented[failed_layer], math.inf)
    for (
        fallback_target,
        inherited_direction_candidate_rank,
    ) in legacy_failure.one_ulp_fallback_candidates:
        if fallback_target not in (lower_target, upper_target):
            continue
        if fallback_target < 0.0 or fallback_target > eligible_mass[failed_layer]:
            continue
        adjusted = dict(represented)
        adjusted[failed_layer] = fallback_target
        adjusted_other_shares = tuple(
            adjusted[layer] for layer in range(1, k_count + 1) if layer != closing_layer
        )
        nonclosing_offset = Fraction.from_float(fallback_target) - Fraction.from_float(
            represented[failed_layer]
        )
        try:
            adjusted_closing_target = _closing_share(
                target_total=accepted,
                other_shares=adjusted_other_shares,
                upper_bound=eligible_mass[closing_layer],
                exact_offset_from=legacy_closing_target,
                required_exact_offset=-nonclosing_offset,
            )
        except KCellTowerRouterConservationError:
            continue
        closing_offset = Fraction.from_float(adjusted_closing_target) - Fraction.from_float(
            legacy_closing_target
        )
        if nonclosing_offset + closing_offset != 0:
            continue
        adjusted[closing_layer] = adjusted_closing_target
        if math.fsum(adjusted[layer] for layer in range(1, k_count + 1)) != accepted:
            continue
        try:
            return build(
                adjusted,
                disposition=KShareApportionmentDisposition.ONE_NONCLOSING_ONE_ULP,
                adjusted_nonclosing_layer=failed_layer,
                selected_nonclosing_signed_ulp_displacement=(
                    _signed_nonnegative_ulp_distance(
                        fallback_target,
                        represented[failed_layer],
                    )
                ),
                inherited_direction_candidate_rank=inherited_direction_candidate_rank,
            )
        except _KShareConstructionExhaustion as joint_failure:
            joint_layer = joint_failure.layer_id.value
            if joint_layer not in (failed_layer, closing_layer):
                raise _legacy_representability_refusal(
                    "more than one nonclosing layer exhausted representability"
                ) from joint_failure
            continue
    raise _legacy_representability_refusal(
        "no exact joint one-ULP nonclosing/closing candidate exists"
    ) from legacy_failure


def _route_totals(entries: tuple[PacketPlanEntry, ...]) -> tuple[float, ...]:
    totals = _sum_entries(entries)
    if any(not math.isfinite(value) for value in totals):
        raise KCellTowerRouterConservationError("route totals are not finite")
    return totals


def _validate_k_share_selection_replay(
    *,
    request: SerialTowerRoutingRequest,
    exit_ledgers: tuple[TrayExitLedger, ...],
    split_lineage: tuple[PacketSplitLineage, ...],
) -> None:
    """Replay the exact v2 tray construction and compare its complete selection record."""

    if (
        type(exit_ledgers) is not tuple
        or type(split_lineage) is not tuple
        or len(exit_ledgers) != len(request.prior_trays)
        or len(request.exit_debits) != len(request.prior_trays)
    ):
        raise KCellTowerRouterConservationError(
            "K-share selection replay differs from the request tray count"
        )
    replayed = tuple(
        _exit_from_tray(request=request, tray=tray, debit=debit)
        for tray, debit in zip(request.prior_trays, request.exit_debits, strict=True)
    )
    expected_ledgers = tuple(result.ledger for result in replayed)
    if exit_ledgers != expected_ledgers:
        raise KCellTowerRouterConservationError(
            "K-share selection record differs from canonical request replay"
        )
    expected_lineage = tuple(split for result in replayed for split in result.split_lineage)
    if split_lineage != expected_lineage:
        raise KCellTowerRouterConservationError(
            "K-share split lineage differs from canonical request replay"
        )


def _request_digest(
    request: SerialTowerRoutingRequest,
    auditor_pin: tuple[str, str],
    planner_pins: tuple[tuple[str, str], ...],
) -> str:
    parts: list[bytes] = [
        b"request",
        b"conservation-audit-law",
        _text_part(ROUTER_CONSERVATION_AUDIT_LAW_ID),
        _text_part(request.attempt_id),
        _float_part(request.event_time_s),
    ]
    for tray in request.prior_trays:
        parts.extend((_text_part(tray.physical_tray_id), _text_part(tray.state_digest)))
    for debit in request.exit_debits:
        parts.extend(
            (
                _text_part(debit.physical_tray_id),
                _float_part(debit.accepted_dry_matter_kg),
                _integer_part(debit.closing_layer_cursor),
            )
        )
    for feed in request.feed_packets:
        parts.extend(_entry_parts(feed.entry))
        parts.extend(
            (
                _float_part(feed.pressure_pa),
                _text_part(feed.configuration_digest),
                _text_part(feed.codec_registry_digest),
                _text_part(feed.auditor_identity_digest),
            )
        )
    policy = request.k_share_closing_exchange_policy
    if policy is None:
        parts.extend((b"k-share-closing-exchange-policy", b"absent"))
    else:
        parts.extend(
            (
                b"k-share-closing-exchange-policy",
                _text_part(policy.policy_id),
                _text_part(policy.candidate_order_law),
                _text_part(policy.ruling_provenance),
                _integer_part(policy.cumulative_signed_drift_tripwire_ulp),
            )
        )
    for planner_id, source_identity in (auditor_pin, *planner_pins):
        parts.extend((_text_part(planner_id), _text_part(source_identity)))
    return _digest(*parts)


def _transaction_digest(
    *,
    request_digest: str,
    candidates: tuple[AcceptedKCellTrayState, ...],
    exits: tuple[TrayExitLedger, ...],
    routes: tuple[AcceptedPopulationRoute, ...],
    splits: tuple[PacketSplitLineage, ...],
    sink: tuple[PacketPlanEntry, ...],
    ledger: WholeRouterLedger,
) -> str:
    if not exits:
        raise KCellTowerRouterConservationError(
            "routing transaction digest requires a nonempty exit-ledger tuple"
        )
    parts: list[bytes] = [
        b"transaction",
        _text_part(request_digest),
        b"K-share-apportionment-policy",
        _text_part(exits[0].apportionment_policy_id),
    ]
    parts.extend(_text_part(candidate.state_digest) for candidate in candidates)
    parts.append(_integer_part(len(exits)))
    for exit_ledger in exits:
        parts.extend(
            (
                _text_part(exit_ledger.physical_tray_id),
                _text_part(exit_ledger.apportionment_policy_id),
                _text_part(exit_ledger.apportionment_disposition.value),
                _integer_part(
                    exit_ledger.adjusted_nonclosing_layer_id.value
                    if exit_ledger.adjusted_nonclosing_layer_id is not None
                    else 0
                ),
                _signed_integer_part(exit_ledger.selected_nonclosing_signed_ulp_displacement),
                _signed_integer_part(
                    exit_ledger.inherited_direction_candidate_rank
                    if exit_ledger.inherited_direction_candidate_rank is not None
                    else -1
                ),
                _signed_integer_part(exit_ledger.selected_closing_signed_ulp_displacement),
                _signed_integer_part(
                    exit_ledger.inherited_closing_candidate_rank
                    if exit_ledger.inherited_closing_candidate_rank is not None
                    else -1
                ),
                _float_part(exit_ledger.accepted_dry_matter_debit_kg),
                _float_part(exit_ledger.represented_dry_matter_debit_kg),
                _float_part(exit_ledger.eligible_dry_matter_kg),
                _integer_part(
                    exit_ledger.closing_layer_id.value
                    if exit_ledger.closing_layer_id is not None
                    else 0
                ),
                _integer_part(exit_ledger.requested_closing_layer_cursor),
                _integer_part(exit_ledger.next_closing_layer_cursor),
                _integer_part(len(exit_ledger.k_shares)),
            )
        )
        for share in exit_ledger.k_shares:
            parts.extend(
                (
                    _integer_part(share.layer_id.value),
                    _float_part(share.eligible_dry_matter_kg),
                    _float_part(share.ideal_dry_matter_debit_kg),
                    _float_part(share.legacy_represented_dry_matter_debit_kg),
                    _float_part(share.represented_dry_matter_debit_kg),
                    _float_part(share.represented_minus_ideal_kg),
                    _signed_integer_part(share.signed_binary64_steps_from_legacy),
                    _text_part(share.apportionment_role.value),
                    _signed_integer_part(
                        share.selected_boundary_weight_candidate_rank
                        if share.selected_boundary_weight_candidate_rank is not None
                        else -1
                    ),
                    _integer_part(share.whole_packet_count),
                    _integer_part(share.boundary_split_count),
                )
            )
        parts.append(b"1" if exit_ledger.passed else b"0")
    for route in routes:
        parts.extend(
            (
                _text_part(route.authority_class.value),
                _text_part(route.topology.source_cell_id),
                _text_part(route.topology.destination_cell_id),
                *(_float_part(value) for value in route.extensive_totals),
            )
        )
        for entry in route.packet_plan:
            parts.extend(_entry_parts(entry))
    parts.append(_text_part(repr(splits)))
    for entry in sink:
        parts.extend(_entry_parts(entry))
    parts.extend(
        (
            _text_part(ledger.audit_law_id),
            *(_float_part(value) for value in ledger.tray_before_totals),
            *(_float_part(value) for value in ledger.feed_totals),
            *(_float_part(value) for value in ledger.tray_after_totals),
            *(_float_part(value) for value in ledger.final_sink_totals),
            *(_fraction_part(value) for value in ledger.raw_exact_signed_entry_residual),
            *(_fraction_part(value) for value in ledger.certified_split_rounding_correction),
            *(_fraction_part(value) for value in ledger.conservation_residual),
            *(_float_part(value) for value in ledger.binary64_regrouping_residual),
            _integer_part(ledger.route_count),
            b"1" if ledger.passed else b"0",
        )
    )
    return _digest(*parts)


def _same_binary64_tuple(left: tuple[float, ...], right: tuple[float, ...]) -> bool:
    return (
        type(left) is tuple
        and type(right) is tuple
        and len(left) == len(right)
        and all(_same_float(a, b) for a, b in zip(left, right, strict=True))
    )


def _validate_transaction_conservation(
    transaction: ValidatedSerialTowerRouting,
) -> None:
    """Recompute the packet-level audit and every declared route total.

    The helper intentionally accepts a structurally valid transaction even when
    it is not an evaluator-issued object so downstream sealed authorities can
    independently fail closed on digest-consistent forgeries.
    """

    if type(transaction) is not ValidatedSerialTowerRouting:
        raise TypeError("conservation validation requires an exact routing transaction")
    ledger = transaction.whole_router_ledger
    if type(ledger) is not WholeRouterLedger:
        raise KCellTowerRouterConservationError(
            "routing transaction carries a foreign whole-router ledger"
        )
    _validate_k_share_apportionment_policy(
        transaction.exit_ledgers,
        closing_exchange_policy=transaction.request.k_share_closing_exchange_policy,
    )
    _validate_k_share_selection_replay(
        request=transaction.request,
        exit_ledgers=transaction.exit_ledgers,
        split_lineage=transaction.split_lineage,
    )
    feed_entries = tuple(feed.entry for feed in transaction.request.feed_packets)
    expected = _compute_whole_router_ledger(
        attempt_id=transaction.request.attempt_id,
        prior_trays=transaction.request.prior_trays,
        feed_entries=feed_entries,
        candidate_trays=transaction.candidate_trays,
        final_sink_entries=transaction.final_sink_packet_plan,
        routes=transaction.routes,
        exit_ledgers=transaction.exit_ledgers,
        split_lineage=transaction.split_lineage,
    )
    if ledger.audit_law_id != expected.audit_law_id:
        raise KCellTowerRouterConservationError(
            "whole-router conservation audit law changed after evaluation"
        )
    for label, observed, recomputed in (
        ("tray-before totals", ledger.tray_before_totals, expected.tray_before_totals),
        ("feed totals", ledger.feed_totals, expected.feed_totals),
        ("tray-after totals", ledger.tray_after_totals, expected.tray_after_totals),
        ("final-sink totals", ledger.final_sink_totals, expected.final_sink_totals),
        (
            "binary64 regrouping residual",
            ledger.binary64_regrouping_residual,
            expected.binary64_regrouping_residual,
        ),
    ):
        if not _same_binary64_tuple(observed, recomputed):
            raise KCellTowerRouterConservationError(
                f"whole-router {label} changed after evaluation"
            )
    for label, observed, recomputed in (
        (
            "raw exact signed-entry residual",
            ledger.raw_exact_signed_entry_residual,
            expected.raw_exact_signed_entry_residual,
        ),
        (
            "certified split-rounding correction",
            ledger.certified_split_rounding_correction,
            expected.certified_split_rounding_correction,
        ),
        (
            "authoritative exact conservation residual",
            ledger.conservation_residual,
            expected.conservation_residual,
        ),
    ):
        if observed != recomputed:
            raise KCellTowerRouterConservationError(f"{label} changed after evaluation")
    if ledger.route_count != expected.route_count or ledger.passed is not expected.passed:
        raise KCellTowerRouterConservationError(
            "whole-router route count or exact-zero verdict changed after evaluation"
        )
    if not expected.passed:
        raise KCellTowerRouterConservationError(
            "whole-router authoritative exact conservation residual is nonzero"
        )
    for route in transaction.routes:
        if type(route) is not AcceptedPopulationRoute:
            raise KCellTowerRouterConservationError(
                "routing transaction contains a foreign route record"
            )
        recomputed = _route_totals(route.packet_plan)
        if not _same_binary64_tuple(route.extensive_totals, recomputed):
            raise KCellTowerRouterConservationError(
                "declared route totals differ from packet-level weighted totals"
            )


class AcceptedKCellTowerRouter:
    """Single-host, sealed, in-memory authority for one serial population route."""

    k_share_apportionment_policy_id: ClassVar[str] = QSC_K_SHARE_APPORTIONMENT_POLICY_ID
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    durable_replay_protection_selected: ClassVar[bool] = False
    in_memory_replay_protection_selected: ClassVar[bool] = True
    physical_k_scheduler_implemented: ClassVar[bool] = False
    mixed_temperature_closure_implemented: ClassVar[bool] = False
    cross_pressure_payload_remap_implemented: ClassVar[bool] = False
    pd_dome_bypass_closure_implemented: ClassVar[bool] = False
    arbitrary_variant_handoff_implemented: ClassVar[bool] = False

    def __init__(self, *, auditor: WeightedPacketPayloadAuditor) -> None:
        identity = getattr(auditor, "auditor_identity_digest", None)
        source = getattr(auditor, "source_identity", None)
        _require_nonblank("router auditor pins", identity, source)
        if not callable(getattr(auditor, "audit", None)):
            raise KCellTowerRouterConfigurationError("router auditor must expose audit(payload)")
        self._auditor = auditor
        self._auditor_identity = identity
        self._auditor_source = source
        self._router_token = object()
        self._consumed_replay_tokens: set[str] = set()
        self._consumed_prior_identities: set[str] = set()
        self._issued_transactions: dict[int, ValidatedSerialTowerRouting] = {}

    def _audit_entry(self, entry: PacketPlanEntry) -> None:
        try:
            observed = self._auditor.audit(entry.payload_bytes)
        except Exception as error:
            raise KCellTowerRouterConfigurationError(
                f"payload auditor refused packet {entry.packet_id!r}: {type(error).__name__}"
            ) from error
        if type(observed) is not PerParticleInventory:
            raise KCellTowerRouterConfigurationError(
                "payload auditor must return an exact PerParticleInventory"
            )
        if observed != entry.inventory:
            raise KCellTowerRouterConfigurationError(
                f"payload audit disagrees with packet {entry.packet_id!r} inventory"
            )
        if entry.surface is not None:
            require_single_surface_authority(entry.inventory, entry.surface)
        derive_packet_totals(entry.inventory, entry.weight)

    def evaluate(
        self,
        request: SerialTowerRoutingRequest,
        *,
        admission_planners: tuple[ManufacturedAdmissionPlanner, ...],
    ) -> ValidatedSerialTowerRouting | TowerRoutingRejection:
        """Evaluate a complete serial route without mutating accepted tray state."""

        if type(request) is not SerialTowerRoutingRequest:
            raise TypeError("tower router requires an exact SerialTowerRoutingRequest")
        try:
            if type(admission_planners) is not tuple or len(admission_planners) != len(
                request.prior_trays
            ):
                raise KCellTowerRouterConfigurationError(
                    "one manufactured admission planner is required per destination tray"
                )
            planner_pins: list[tuple[str, str]] = []
            for planner in admission_planners:
                planner_id = getattr(planner, "planner_id", None)
                source_identity = getattr(planner, "source_identity", None)
                _require_nonblank("admission planner pins", planner_id, source_identity)
                if not callable(getattr(planner, "assign", None)):
                    raise KCellTowerRouterConfigurationError(
                        "admission planner must expose assign(request)"
                    )
                planner_pins.append((planner_id, source_identity))
            pins = tuple(planner_pins)
            auditor_pin = (self._auditor_identity, self._auditor_source)
            if any(
                tray.auditor_identity_digest != self._auditor_identity
                for tray in request.prior_trays
            ):
                raise KCellTowerRouterConfigurationError(
                    "router auditor differs from the accepted tray authority"
                )
            for tray in request.prior_trays:
                for packet in tray.packets:
                    self._audit_entry(packet.entry)
            for feed in request.feed_packets:
                self._audit_entry(feed.entry)
            exits = tuple(
                _exit_from_tray(request=request, tray=tray, debit=debit)
                for tray, debit in zip(request.prior_trays, request.exit_debits, strict=True)
            )
            exit_ledger_tuple = tuple(exit_result.ledger for exit_result in exits)
            split_lineage = tuple(
                split for exit_result in exits for split in exit_result.split_lineage
            )
            _validate_k_share_apportionment_policy(
                exit_ledger_tuple,
                closing_exchange_policy=request.k_share_closing_exchange_policy,
            )
            _validate_k_share_selection_replay(
                request=request,
                exit_ledgers=exit_ledger_tuple,
                split_lineage=split_lineage,
            )
            incoming_by_destination: list[tuple[PacketPlanEntry, ...]] = [
                tuple(
                    sorted(
                        (feed.entry for feed in request.feed_packets),
                        key=lambda entry: entry.packet_id,
                    )
                )
            ]
            incoming_by_destination.extend(exit.outgoing_entries for exit in exits[:-1])
            assignments_by_destination: list[dict[str, AdmissionAssignment]] = []
            for tray, incoming, planner in zip(
                request.prior_trays,
                incoming_by_destination,
                admission_planners,
                strict=True,
            ):
                planning_request = AdmissionPlanningRequest(
                    attempt_id=request.attempt_id,
                    destination_prior=tray,
                    incoming_entries=incoming,
                )
                try:
                    assignments = planner.assign(planning_request)
                except Exception as error:
                    raise KCellTowerRouterConfigurationError(
                        f"manufactured admission planner failed: {type(error).__name__}: {error}"
                    ) from error
                if type(assignments) is not tuple or any(
                    type(item) is not AdmissionAssignment for item in assignments
                ):
                    raise KCellTowerRouterConfigurationError(
                        "admission planner must return an exact tuple of exact assignments"
                    )
                incoming_ids = tuple(entry.packet_id for entry in incoming)
                assignment_ids = tuple(item.packet_id for item in assignments)
                if assignment_ids != incoming_ids:
                    raise KCellTowerRouterConfigurationError(
                        "admission assignments must exactly follow canonical incoming packet order"
                    )
                physical = next(
                    item
                    for item in tray.topology.trays
                    if item.physical_tray_id == tray.physical_tray_id
                )
                for assignment in assignments:
                    if assignment.authority_id != planner.planner_id:
                        raise KCellTowerRouterConfigurationError(
                            "admission assignment violates its planner identity"
                        )
                    if assignment.destination_layer_id.value > physical.vertical_layer_count:
                        raise KCellTowerRouterConfigurationError(
                            "admission assignment exceeds destination K"
                        )
                    if (
                        len(assignment.destination_rtd_path.transition_exposures)
                        != physical.rtd_stage_count
                    ):
                        raise KCellTowerRouterConfigurationError(
                            "admission RTD path order differs from destination topology"
                        )
                assignments_by_destination.append(
                    {assignment.packet_id: assignment for assignment in assignments}
                )
            candidates: list[AcceptedKCellTrayState] = []
            for tray, exit_result, incoming, assignments in zip(
                request.prior_trays,
                exits,
                incoming_by_destination,
                assignments_by_destination,
                strict=True,
            ):
                revision = tray.revision + 1
                packets = [
                    replace(packet, revision=revision, time_s=request.event_time_s)
                    for packet in exit_result.remaining_packets
                ]
                for entry in incoming:
                    assignment = assignments[entry.packet_id]
                    packets.append(
                        AcceptedTrayPacket(
                            owner_key=ParticlePopulationKey(
                                tower_id=tray.topology.tower_id,
                                physical_tray_id=tray.physical_tray_id,
                                vertical_layer_id=assignment.destination_layer_id,
                                rtd_stage_id=RTDStageId(1),
                                cohort_id=entry.packet_id,
                            ),
                            entry=entry,
                            revision=revision,
                            time_s=request.event_time_s,
                            pressure_pa=tray.pressure_pa,
                            configuration_digest=tray.configuration_digest,
                            codec_registry_digest=tray.codec_registry_digest,
                            auditor_identity_digest=tray.auditor_identity_digest,
                            rtd_path=assignment.destination_rtd_path,
                            rtd_schedule_digest=tray.rtd_schedule_digest,
                            cumulative_rtd_exposure=0.0,
                        )
                    )
                ids = tuple(packet.packet_id for packet in packets)
                if len(ids) != len(set(ids)):
                    raise KCellTowerRouterConfigurationError(
                        f"routing creates duplicate packet identity in {tray.physical_tray_id}"
                    )
                cells = tuple(
                    AcceptedKCell(
                        vertical_layer_id=cell.vertical_layer_id,
                        conserved_state=replace(
                            cell.conserved_state,
                            revision=revision,
                            state_digest=_digest(
                                b"routing-cell-revision",
                                _text_part(request.attempt_id),
                                _text_part(cell.conserved_state.state_digest),
                                _text_part(tray.physical_tray_id),
                                _integer_part(cell.vertical_layer_id.value),
                                _integer_part(revision),
                            ),
                        ),
                    )
                    for cell in tray.cells
                )
                candidate = AcceptedKCellTrayState(
                    topology=tray.topology,
                    physical_tray_id=tray.physical_tray_id,
                    revision=revision,
                    time_s=request.event_time_s,
                    pressure_pa=tray.pressure_pa,
                    configuration_digest=tray.configuration_digest,
                    codec_registry_digest=tray.codec_registry_digest,
                    auditor_identity_digest=tray.auditor_identity_digest,
                    energy_datum_id=tray.energy_datum_id,
                    rtd_schedule_digest=tray.rtd_schedule_digest,
                    cells=cells,
                    packets=tuple(sorted(packets, key=_packet_sort_key)),
                )
                for packet in candidate.packets:
                    self._audit_entry(packet.entry)
                candidates.append(candidate)
            candidate_tuple = tuple(candidates)
            sink = exits[-1].outgoing_entries
            routed_ids = tuple(
                packet.packet_id for candidate in candidate_tuple for packet in candidate.packets
            ) + tuple(entry.packet_id for entry in sink)
            if len(routed_ids) != len(set(routed_ids)):
                raise KCellTowerRouterConfigurationError(
                    "routing creates a duplicate packet identity across the tower or final sink"
                )
            routes: list[AcceptedPopulationRoute] = []
            feed_entries = incoming_by_destination[0]
            if feed_entries:
                routes.append(
                    AcceptedPopulationRoute(
                        authority_class=TransferAuthorityClass.FEED_ADMISSION,
                        topology=TransferTopology(
                            source_cell_id="external-feed-buffer",
                            destination_cell_id=request.prior_trays[0].physical_tray_id,
                        ),
                        packet_plan=feed_entries,
                        extensive_totals=_route_totals(feed_entries),
                    )
                )
            for index, exit_result in enumerate(exits[:-1]):
                if not exit_result.outgoing_entries:
                    continue
                routes.append(
                    AcceptedPopulationRoute(
                        authority_class=TransferAuthorityClass.INTER_TRAY,
                        topology=TransferTopology(
                            source_cell_id=request.prior_trays[index].physical_tray_id,
                            destination_cell_id=request.prior_trays[index + 1].physical_tray_id,
                        ),
                        packet_plan=exit_result.outgoing_entries,
                        extensive_totals=_route_totals(exit_result.outgoing_entries),
                    )
                )
            if sink:
                routes.append(
                    AcceptedPopulationRoute(
                        authority_class=TransferAuthorityClass.DT_TO_DC_HANDOFF,
                        topology=TransferTopology(
                            source_cell_id=request.prior_trays[-1].physical_tray_id,
                            destination_cell_id="external-final-sink",
                        ),
                        packet_plan=sink,
                        extensive_totals=_route_totals(sink),
                    )
                )
            route_tuple = tuple(routes)
            ledger = _compute_whole_router_ledger(
                attempt_id=request.attempt_id,
                prior_trays=request.prior_trays,
                feed_entries=feed_entries,
                candidate_trays=candidate_tuple,
                final_sink_entries=sink,
                routes=route_tuple,
                exit_ledgers=exit_ledger_tuple,
                split_lineage=split_lineage,
            )
            if not ledger.passed:
                raise KCellTowerRouterConservationError(
                    "whole-router exact seven-extensive residual is "
                    f"{ledger.conservation_residual!r}"
                )
            request_identity = _request_digest(request, auditor_pin, pins)
            transaction_identity = _transaction_digest(
                request_digest=request_identity,
                candidates=candidate_tuple,
                exits=exit_ledger_tuple,
                routes=route_tuple,
                splits=split_lineage,
                sink=sink,
                ledger=ledger,
            )
            validated = ValidatedSerialTowerRouting(
                request=request,
                candidate_trays=candidate_tuple,
                exit_ledgers=exit_ledger_tuple,
                routes=route_tuple,
                split_lineage=split_lineage,
                final_sink_packet_plan=sink,
                whole_router_ledger=ledger,
                auditor_pin=auditor_pin,
                planner_pins=pins,
                request_digest=request_identity,
                transaction_digest=transaction_identity,
                replay_token=_digest(b"replay", _text_part(transaction_identity)),
                _router_token=self._router_token,
                _seal=_VALIDATED_SEAL,
            )
            self._issued_transactions[id(validated)] = validated
            return validated
        except (KCellTowerRouterError, TypeError, ValueError, OverflowError) as error:
            return TowerRoutingRejection(
                attempt_id=request.attempt_id,
                reason=f"{type(error).__name__}: {error}",
                rollback_trays=request.prior_trays,
            )

    def verify_issued_transaction(
        self,
        transaction: ValidatedSerialTowerRouting,
    ) -> ValidatedSerialTowerRouting:
        """Fail closed unless ``transaction`` is the exact evaluator-issued object."""

        if (
            type(transaction) is not ValidatedSerialTowerRouting
            or transaction._seal is not _VALIDATED_SEAL
        ):
            raise KCellTowerRouterError("verification requires a router-validated transaction")
        if transaction._router_token is not self._router_token:
            raise KCellTowerRouterError("routing transaction was issued by another router")
        if self._issued_transactions.get(id(transaction)) is not transaction:
            raise KCellTowerRouterError(
                "routing transaction is not the exact evaluator-issued in-memory object"
            )
        current_auditor_pin = (
            getattr(self._auditor, "auditor_identity_digest", None),
            getattr(self._auditor, "source_identity", None),
        )
        if current_auditor_pin != transaction.auditor_pin:
            raise KCellTowerRouterError("router payload auditor pin changed after evaluation")
        expected_request = _request_digest(
            transaction.request,
            transaction.auditor_pin,
            transaction.planner_pins,
        )
        if expected_request != transaction.request_digest:
            raise KCellTowerRouterError("routing request identity changed after evaluation")
        _validate_transaction_conservation(transaction)
        expected_transaction = _transaction_digest(
            request_digest=transaction.request_digest,
            candidates=transaction.candidate_trays,
            exits=transaction.exit_ledgers,
            routes=transaction.routes,
            splits=transaction.split_lineage,
            sink=transaction.final_sink_packet_plan,
            ledger=transaction.whole_router_ledger,
        )
        if expected_transaction != transaction.transaction_digest:
            raise KCellTowerRouterError("routing transaction identity changed after evaluation")
        expected_replay = _digest(b"replay", _text_part(expected_transaction))
        if transaction.replay_token != expected_replay:
            raise KCellTowerRouterError("routing replay identity changed after evaluation")
        return transaction

    def verify_commit_available(
        self,
        current_trays: tuple[AcceptedKCellTrayState, ...],
        transaction: ValidatedSerialTowerRouting,
    ) -> ValidatedSerialTowerRouting:
        """Verify one route can still commit without consuming any authority state."""

        self.verify_issued_transaction(transaction)
        if transaction.replay_token in self._consumed_replay_tokens:
            raise ReplayedKCellTowerRoutingError("routing replay token was already consumed")
        prior_identity = _digest(
            b"accepted-prior-tuple",
            *(_text_part(prior.state_digest) for prior in transaction.request.prior_trays),
        )
        if prior_identity in self._consumed_prior_identities:
            raise StaleKCellTowerRoutingError(
                "another routing transaction already consumed this accepted prior tuple"
            )
        if type(current_trays) is not tuple or len(current_trays) != len(
            transaction.request.prior_trays
        ):
            raise StaleKCellTowerRoutingError("current tray tuple shape differs from routing prior")
        if any(
            current is not prior
            for current, prior in zip(current_trays, transaction.request.prior_trays, strict=True)
        ):
            raise StaleKCellTowerRoutingError("routing transaction is stale or foreign")
        for current in current_trays:
            for packet in current.packets:
                self._audit_entry(packet.entry)
        for candidate, current in zip(transaction.candidate_trays, current_trays, strict=True):
            if candidate.revision != current.revision + 1:
                raise KCellTowerRouterError("routing candidate revision did not advance once")
            for packet in candidate.packets:
                self._audit_entry(packet.entry)
        _validate_transaction_conservation(transaction)
        if not transaction.whole_router_ledger.passed or any(
            value != 0 for value in transaction.whole_router_ledger.conservation_residual
        ):
            raise KCellTowerRouterConservationError("routing conservation certificate failed")
        return transaction

    def commit(
        self,
        current_trays: tuple[AcceptedKCellTrayState, ...],
        transaction: ValidatedSerialTowerRouting,
    ) -> tuple[AcceptedKCellTrayState, ...]:
        """Consume one sealed route and return its immutable candidate tray tuple."""

        self.verify_commit_available(current_trays, transaction)
        prior_identity = _digest(
            b"accepted-prior-tuple",
            *(_text_part(prior.state_digest) for prior in transaction.request.prior_trays),
        )
        self._consumed_replay_tokens.add(transaction.replay_token)
        self._consumed_prior_identities.add(prior_identity)
        return transaction.candidate_trays


__all__ = (
    "AcceptedFeedPacket",
    "AcceptedKCellTowerRouter",
    "AcceptedPopulationRoute",
    "AcceptedTrayExitDebit",
    "AdmissionAssignment",
    "AdmissionPlanningRequest",
    "KCellTowerRouterConfigurationError",
    "KCellTowerRouterConservationError",
    "KCellTowerRouterError",
    "KShareApportionmentDisposition",
    "KShareApportionmentRole",
    "KShareExitLedger",
    "ManufacturedAdmissionPlanner",
    "PacketSplitLineage",
    "QSC_K_SHARE_APPORTIONMENT_POLICY_ID",
    "ROUTER_ID",
    "ROUTER_CONSERVATION_AUDIT_LAW_ID",
    "ReplayedKCellTowerRoutingError",
    "SerialTowerRoutingRequest",
    "StaleKCellTowerRoutingError",
    "TowerRoutingRejection",
    "TrayExitLedger",
    "ValidatedSerialTowerRouting",
    "WholeRouterLedger",
)
