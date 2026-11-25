"""Sealed in-memory reconciliation of tower routing and engineering K clocks.

``accepted_k_cell_tower_router`` changes packet membership instantaneously at an
already accepted common time.  ``engineering_physical_k_event_scheduler`` owns
one immutable stochastic clock per packet, so a routed tower cannot retain its
pre-route scheduler tuple unchanged.  This module supplies that narrow seam:

* a packet that remains on the same tray retains its exact clock;
* a boundary-split remainder retains the parent clock explicitly;
* a whole packet or split child admitted to another tray receives a fresh,
  deterministic clock under the destination authority;
* feed packets receive the same destination-authority initialization; and
* whole-packet exits and final-sink split children carry no tower clock.

The reconciler accepts only a self-consistent, digest-valid router transaction.
It does not infer missing lineage, alter payloads, route mass, advance physical
time, authenticate RPM histories, or provide durable atomic persistence.  All
physical, predictive, production, and persistence claims remain false.
"""

from __future__ import annotations

import enum
import hashlib
import math
from dataclasses import dataclass, field, replace
from typing import ClassVar, Iterable

from . import accepted_k_cell_tower_router as tower_router
from . import engineering_physical_k_event_scheduler as k_scheduler
from . import k_cell_tray_host as tray_host
from .accepted_solid_transfer_host import TransferAuthorityClass


RECONCILER_SCHEMA_ID = "dtdc-core2-engineering-k-tower-route-reconciler-v1"
RECONCILER_SOURCE_IDENTITY = "sealed-router-lineage-to-destination-k-clock-v1"
_VALIDATED_RECONCILIATION_SEAL = object()


class EngineeringKTowerRouteReconciliationError(ValueError):
    """Base refusal for the bounded router/K-clock reconciliation seam."""


class EngineeringKTowerRouteConfigurationError(EngineeringKTowerRouteReconciliationError):
    """The route, authorities, or schedule tuple are mutually inconsistent."""


class StaleEngineeringKTowerRouteStateError(EngineeringKTowerRouteReconciliationError):
    """Commit was attempted against a foreign or stale prior tuple."""


class ReplayedEngineeringKTowerRouteError(EngineeringKTowerRouteReconciliationError):
    """A reconciliation transaction or its accepted prior was consumed already."""


def _require_identity(label: str, *values: str) -> None:
    if any(type(value) is not str or not value or value != value.strip() for value in values):
        raise EngineeringKTowerRouteConfigurationError(
            f"{label} must contain nonblank exact strings"
        )


def _digest(*, domain: str, parts: Iterable[str | int | bool]) -> str:
    _require_identity("digest domain", domain)
    digest = hashlib.sha256()
    digest.update(b"DTDC-ENGINEERING-K-TOWER-ROUTE-RECONCILER-V1\0")
    material = tuple(parts)
    for value in (domain, str(len(material))):
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    for part in material:
        if isinstance(part, bool):
            tag, payload = b"b", b"1" if part else b"0"
        elif type(part) is int:
            tag, payload = b"i", str(part).encode("ascii")
        elif type(part) is str:
            _require_identity("digest string", part)
            tag, payload = b"s", part.encode("utf-8")
        else:  # pragma: no cover - closed union guarded by callers
            raise TypeError(f"unsupported reconciler digest part {type(part).__name__}")
        digest.update(tag)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return "sha256:" + digest.hexdigest()


def _payload_identity(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def _packet_by_id(
    tray: tray_host.AcceptedKCellTrayState,
) -> dict[str, tray_host.AcceptedTrayPacket]:
    return {packet.packet_id: packet for packet in tray.packets}


def _clock_by_id(
    state: k_scheduler.EngineeringKScheduleState,
) -> dict[str, k_scheduler.EngineeringKPacketClock]:
    return {clock.packet_id: clock for clock in state.packet_clocks}


class KClockRouteAction(enum.Enum):
    """Exhaustive clock disposition at one instantaneous routing event."""

    PRESERVE_SURVIVOR = "preserve_survivor"
    PRESERVE_SPLIT_REMAINDER = "preserve_split_remainder"
    INITIALIZE_FEED = "initialize_feed"
    INITIALIZE_INTERTRAY_WHOLE = "initialize_intertray_whole"
    INITIALIZE_INTERTRAY_SPLIT_CHILD = "initialize_intertray_split_child"
    DROP_FINAL_SINK_WHOLE = "drop_final_sink_whole"
    DROP_FINAL_SINK_SPLIT_CHILD = "drop_final_sink_split_child"


@dataclass(frozen=True, slots=True, kw_only=True)
class KClockRouteRecord:
    """Auditable disposition of one surviving, admitted, or final-sink clock."""

    packet_id: str
    action: KClockRouteAction
    source_packet_id: str | None
    source_tray_id: str
    destination_tray_id: str
    source_clock_digest: str | None
    candidate_clock_digest: str | None

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity(
            "clock-route record",
            self.packet_id,
            self.source_tray_id,
            self.destination_tray_id,
        )
        if type(self.action) is not KClockRouteAction:
            raise TypeError("clock-route action must be an exact KClockRouteAction")
        if self.source_packet_id is not None:
            _require_identity("clock-route source packet", self.source_packet_id)
        if self.source_clock_digest is not None:
            _require_identity("clock-route source clock", self.source_clock_digest)
        if self.candidate_clock_digest is not None:
            _require_identity("clock-route candidate clock", self.candidate_clock_digest)

    @property
    def definition_digest(self) -> str:
        return _digest(
            domain="k-clock-route-record",
            parts=(
                self.packet_id,
                self.action.value,
                self.source_packet_id or "none:source-packet",
                self.source_tray_id,
                self.destination_tray_id,
                self.source_clock_digest or "none:source-clock",
                self.candidate_clock_digest or "none:candidate-clock",
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringKTowerClockLedger:
    prior_clock_count: int
    candidate_clock_count: int
    preserved_survivor_count: int
    preserved_split_remainder_count: int
    initialized_feed_count: int
    initialized_intertray_whole_count: int
    initialized_intertray_split_child_count: int
    dropped_final_sink_whole_count: int
    dropped_final_sink_split_child_count: int
    route_entry_count: int
    split_lineage_count: int
    candidate_packet_sets_match: bool
    packet_k_owners_match: bool
    collision_free: bool
    every_route_entry_accounted: bool
    every_prior_clock_accounted: bool

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.candidate_packet_sets_match
            and self.packet_k_owners_match
            and self.collision_free
            and self.every_route_entry_accounted
            and self.every_prior_clock_accounted
            and self.prior_clock_count
            + self.initialized_feed_count
            + self.initialized_intertray_split_child_count
            - self.dropped_final_sink_whole_count
            == self.candidate_clock_count
        )

    @property
    def definition_digest(self) -> str:
        return _digest(
            domain="engineering-k-tower-clock-ledger",
            parts=(
                self.prior_clock_count,
                self.candidate_clock_count,
                self.preserved_survivor_count,
                self.preserved_split_remainder_count,
                self.initialized_feed_count,
                self.initialized_intertray_whole_count,
                self.initialized_intertray_split_child_count,
                self.dropped_final_sink_whole_count,
                self.dropped_final_sink_split_child_count,
                self.route_entry_count,
                self.split_lineage_count,
                self.candidate_packet_sets_match,
                self.packet_k_owners_match,
                self.collision_free,
                self.every_route_entry_accounted,
                self.every_prior_clock_accounted,
                self.passed,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringKTowerRouteReconciliationRequest:
    """One evaluator-issued population route and its exact pre-route K states."""

    attempt_id: str
    routing_transaction: tower_router.ValidatedSerialTowerRouting
    prior_schedule_states: tuple[k_scheduler.EngineeringKScheduleState, ...]

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("reconciliation attempt", self.attempt_id)
        if type(self.routing_transaction) is not tower_router.ValidatedSerialTowerRouting:
            raise TypeError("routing_transaction must be an exact validated router transaction")
        if type(self.prior_schedule_states) is not tuple or not self.prior_schedule_states:
            raise TypeError("prior_schedule_states must be a nonempty exact tuple")
        if any(
            type(state) is not k_scheduler.EngineeringKScheduleState
            for state in self.prior_schedule_states
        ):
            raise TypeError("prior_schedule_states contain a foreign exact type")


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringKTowerRouteReconciliationRejection:
    attempt_id: str
    reason: str
    rollback_schedule_states: tuple[k_scheduler.EngineeringKScheduleState, ...]

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedEngineeringKTowerRouteReconciliation:
    request: EngineeringKTowerRouteReconciliationRequest
    candidate_schedule_states: tuple[k_scheduler.EngineeringKScheduleState, ...]
    records: tuple[KClockRouteRecord, ...]
    ledger: EngineeringKTowerClockLedger
    authority_definition_digests: tuple[str, ...]
    request_digest: str
    transaction_digest: str
    replay_token: str
    _reconciler_token: object = field(repr=False, compare=False)
    _seal: object = field(repr=False, compare=False)

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    durable_replay_protection_selected: ClassVar[bool] = False
    in_memory_replay_protection_selected: ClassVar[bool] = True

    def __post_init__(self) -> None:
        if self._seal is not _VALIDATED_RECONCILIATION_SEAL:
            raise TypeError("validated reconciliations are issued only by the reconciler")
        if not self.ledger.passed:
            raise EngineeringKTowerRouteConfigurationError(
                "validated reconciliation carries a failed clock ledger"
            )


def _validate_router_transaction(
    transaction: tower_router.ValidatedSerialTowerRouting,
) -> None:
    """Recompute the router's own canonical identities before trusting lineage."""

    try:
        tower_router._validate_transaction_conservation(transaction)  # noqa: SLF001
    except tower_router.KCellTowerRouterError as error:
        raise EngineeringKTowerRouteConfigurationError(
            f"router conservation audit failed: {error}"
        ) from error
    request_digest = tower_router._request_digest(  # noqa: SLF001
        transaction.request,
        transaction.auditor_pin,
        transaction.planner_pins,
    )
    if request_digest != transaction.request_digest:
        raise EngineeringKTowerRouteConfigurationError(
            "router request identity changed after evaluation"
        )
    transaction_digest = tower_router._transaction_digest(  # noqa: SLF001
        request_digest=request_digest,
        candidates=transaction.candidate_trays,
        exits=transaction.exit_ledgers,
        routes=transaction.routes,
        splits=transaction.split_lineage,
        sink=transaction.final_sink_packet_plan,
        ledger=transaction.whole_router_ledger,
    )
    if transaction_digest != transaction.transaction_digest:
        raise EngineeringKTowerRouteConfigurationError(
            "router transaction identity changed after evaluation"
        )
    replay_token = tower_router._digest(  # noqa: SLF001
        b"replay",
        tower_router._text_part(transaction_digest),  # noqa: SLF001
    )
    if replay_token != transaction.replay_token:
        raise EngineeringKTowerRouteConfigurationError(
            "router replay identity changed after evaluation"
        )
    if (
        transaction.whole_router_ledger.audit_law_id
        != tower_router.ROUTER_CONSERVATION_AUDIT_LAW_ID
        or not transaction.whole_router_ledger.passed
        or any(value != 0 for value in transaction.whole_router_ledger.conservation_residual)
    ):
        raise EngineeringKTowerRouteConfigurationError(
            "router transaction does not carry an exact passing seven-extensive ledger"
        )


def _request_digest(
    request: EngineeringKTowerRouteReconciliationRequest,
    authority_digests: tuple[str, ...],
) -> str:
    return _digest(
        domain="engineering-k-tower-route-reconciliation-request",
        parts=(
            RECONCILER_SCHEMA_ID,
            request.attempt_id,
            request.routing_transaction.transaction_digest,
            len(request.prior_schedule_states),
            *(state.definition_digest for state in request.prior_schedule_states),
            *authority_digests,
        ),
    )


def _transaction_digest(
    *,
    request_digest: str,
    candidates: tuple[k_scheduler.EngineeringKScheduleState, ...],
    records: tuple[KClockRouteRecord, ...],
    ledger: EngineeringKTowerClockLedger,
) -> str:
    return _digest(
        domain="engineering-k-tower-route-reconciliation-transaction",
        parts=(
            request_digest,
            len(candidates),
            *(state.definition_digest for state in candidates),
            len(records),
            *(record.definition_digest for record in records),
            ledger.definition_digest,
        ),
    )


class EngineeringKTowerRouteReconciler:
    """Single-host authority for one exact ordered tuple of tray K authorities."""

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    durable_replay_protection_selected: ClassVar[bool] = False
    in_memory_replay_protection_selected: ClassVar[bool] = True
    accepted_history_authenticated: ClassVar[bool] = False
    route_lineage_authenticated_beyond_router_digest: ClassVar[bool] = True

    def __init__(
        self,
        *,
        router: tower_router.AcceptedKCellTowerRouter,
        authorities: tuple[k_scheduler.EngineeringKEventAuthority, ...],
    ) -> None:
        if type(router) is not tower_router.AcceptedKCellTowerRouter:
            raise TypeError("router must be the exact issuing AcceptedKCellTowerRouter")
        if type(authorities) is not tuple or not authorities:
            raise TypeError("authorities must be a nonempty exact tuple")
        if any(
            type(authority) is not k_scheduler.EngineeringKEventAuthority
            for authority in authorities
        ):
            raise TypeError("authorities contain a foreign exact type")
        tray_ids = tuple(authority.mesh.physical_tray_id for authority in authorities)
        if len(set(tray_ids)) != len(tray_ids):
            raise EngineeringKTowerRouteConfigurationError(
                "destination K authorities must have unique physical tray IDs"
            )
        self._router = router
        self._authorities = authorities
        self._authority_digests = tuple(authority.definition_digest for authority in authorities)
        self._reconciler_token = object()
        self._consumed_replay_tokens: set[str] = set()
        self._consumed_prior_identities: set[str] = set()
        self._issued_transactions: dict[int, ValidatedEngineeringKTowerRouteReconciliation] = {}

    def _validate_shape_and_priors(
        self,
        request: EngineeringKTowerRouteReconciliationRequest,
    ) -> None:
        routing = request.routing_transaction
        priors = routing.request.prior_trays
        candidates = routing.candidate_trays
        schedules = request.prior_schedule_states
        count = len(self._authorities)
        if not len(priors) == len(candidates) == len(schedules) == count:
            raise EngineeringKTowerRouteConfigurationError(
                "router trays, K states, and destination authorities have different shapes"
            )
        expected_ids = tuple(tray.physical_tray_id for tray in priors)
        if tuple(tray.physical_tray_id for tray in candidates) != expected_ids:
            raise EngineeringKTowerRouteConfigurationError(
                "router candidate tray order differs from its accepted prior order"
            )
        if tuple(state.physical_tray_id for state in schedules) != expected_ids:
            raise EngineeringKTowerRouteConfigurationError(
                "pre-route K state order differs from the router tray order"
            )
        if (
            tuple(authority.mesh.physical_tray_id for authority in self._authorities)
            != expected_ids
        ):
            raise EngineeringKTowerRouteConfigurationError(
                "destination K authority order differs from the router tray order"
            )
        event_time = routing.request.event_time_s
        for tray, state, authority in zip(priors, schedules, self._authorities, strict=True):
            if state.time_s != event_time or tray.time_s != event_time:
                raise EngineeringKTowerRouteConfigurationError(
                    "routing and every pre-route K state must share one exact accepted time"
                )
            k_scheduler._validate_state_authority(state, authority)  # noqa: SLF001
            packet_map = _packet_by_id(tray)
            clock_map = _clock_by_id(state)
            if set(packet_map) != set(clock_map):
                raise EngineeringKTowerRouteConfigurationError(
                    f"{tray.physical_tray_id} router prior and K state have different packet sets"
                )
            if any(
                packet_map[packet_id].owner_key.vertical_layer_id != clock.vertical_layer_id
                for packet_id, clock in clock_map.items()
            ):
                raise EngineeringKTowerRouteConfigurationError(
                    f"{tray.physical_tray_id} router prior and K state disagree on packet K"
                )

    def _reconcile(
        self,
        request: EngineeringKTowerRouteReconciliationRequest,
    ) -> tuple[
        tuple[k_scheduler.EngineeringKScheduleState, ...],
        tuple[KClockRouteRecord, ...],
        EngineeringKTowerClockLedger,
    ]:
        routing = request.routing_transaction
        priors = routing.request.prior_trays
        candidates = routing.candidate_trays
        schedules = request.prior_schedule_states
        event_time_s = routing.request.event_time_s
        prior_by_tray = {tray.physical_tray_id: tray for tray in priors}
        candidate_by_tray = {tray.physical_tray_id: tray for tray in candidates}
        schedule_by_tray = {state.physical_tray_id: state for state in schedules}

        route_by_packet: dict[str, tower_router.AcceptedPopulationRoute] = {}
        entry_by_packet = {}
        for route in routing.routes:
            for entry in route.packet_plan:
                if entry.packet_id in route_by_packet:
                    raise EngineeringKTowerRouteConfigurationError(
                        "one packet appears on more than one simultaneous population route"
                    )
                route_by_packet[entry.packet_id] = route
                entry_by_packet[entry.packet_id] = entry

        lineage_by_child = {
            lineage.transfer_child_packet_id: lineage for lineage in routing.split_lineage
        }
        if len(lineage_by_child) != len(routing.split_lineage):
            raise EngineeringKTowerRouteConfigurationError(
                "router split lineage repeats a transfer-child identity"
            )
        lineage_by_parent = {
            (lineage.source_tray_id, lineage.parent_packet_id): lineage
            for lineage in routing.split_lineage
        }
        if len(lineage_by_parent) != len(routing.split_lineage):
            raise EngineeringKTowerRouteConfigurationError(
                "router split lineage repeats one source parent"
            )

        # Validate every split against the source prior, source remainder, and
        # one authenticated route entry before any clock is issued.
        for child_id, lineage in lineage_by_child.items():
            if lineage.parent_packet_id != lineage.source_remainder_packet_id:
                raise EngineeringKTowerRouteConfigurationError(
                    "the bounded router must preserve the parent ID for a split remainder"
                )
            source = prior_by_tray.get(lineage.source_tray_id)
            source_candidate = candidate_by_tray.get(lineage.source_tray_id)
            if source is None or source_candidate is None:
                raise EngineeringKTowerRouteConfigurationError(
                    "split lineage names a tray outside the routed tower"
                )
            parents = _packet_by_id(source)
            remainders = _packet_by_id(source_candidate)
            if (
                lineage.parent_packet_id not in parents
                or lineage.source_remainder_packet_id not in (remainders)
            ):
                raise EngineeringKTowerRouteConfigurationError(
                    "split lineage lacks its authenticated parent or source remainder"
                )
            parent = parents[lineage.parent_packet_id]
            remainder = remainders[lineage.source_remainder_packet_id]
            route = route_by_packet.get(child_id)
            child_entry = entry_by_packet.get(child_id)
            if route is None or child_entry is None:
                raise EngineeringKTowerRouteConfigurationError(
                    "split transfer child is absent from the authenticated route plan"
                )
            if route.topology.source_cell_id != lineage.source_tray_id:
                raise EngineeringKTowerRouteConfigurationError(
                    "split lineage and population route name different source trays"
                )
            if parent.owner_key.vertical_layer_id != lineage.source_layer_id:
                raise EngineeringKTowerRouteConfigurationError(
                    "split lineage and source parent name different K layers"
                )
            if (
                parent.packet_id != remainder.packet_id
                or parent.entry.payload_bytes != remainder.entry.payload_bytes
                or parent.entry.inventory != remainder.entry.inventory
                or parent.entry.payload_bytes != child_entry.payload_bytes
                or parent.entry.inventory != child_entry.inventory
                or lineage.payload_identity != _payload_identity(parent.entry.payload_bytes)
                or lineage.parent_weight != parent.entry.weight.representative_particles
                or lineage.remainder_weight != remainder.entry.weight.representative_particles
                or lineage.transfer_weight != child_entry.weight.representative_particles
                or math.fsum((lineage.remainder_weight, lineage.transfer_weight))
                != lineage.parent_weight
            ):
                raise EngineeringKTowerRouteConfigurationError(
                    "split lineage does not exactly reproduce parent, remainder, and child"
                )

        records: list[KClockRouteRecord] = []
        candidate_states: list[k_scheduler.EngineeringKScheduleState] = []
        accounted_route_ids: set[str] = set()
        accounted_prior_keys: set[tuple[str, str]] = set()
        candidate_packet_sets_match = True
        packet_k_owners_match = True

        for prior_tray, candidate_tray, prior_state, authority in zip(
            priors,
            candidates,
            schedules,
            self._authorities,
            strict=True,
        ):
            prior_packets = _packet_by_id(prior_tray)
            candidate_packets = _packet_by_id(candidate_tray)
            prior_clocks = _clock_by_id(prior_state)
            incoming = {
                packet_id: route
                for packet_id, route in route_by_packet.items()
                if route.topology.destination_cell_id == candidate_tray.physical_tray_id
            }
            survivor_ids = set(candidate_packets) & set(prior_packets)
            if set(candidate_packets) != survivor_ids | set(incoming):
                candidate_packet_sets_match = False
                raise EngineeringKTowerRouteConfigurationError(
                    f"{candidate_tray.physical_tray_id} candidate packets are neither "
                    "same-tray survivors nor authenticated admissions"
                )

            clocks_after: list[k_scheduler.EngineeringKPacketClock] = []
            for packet_id in sorted(candidate_packets):
                candidate_packet = candidate_packets[packet_id]
                if packet_id in survivor_ids:
                    prior_packet = prior_packets[packet_id]
                    clock = prior_clocks[packet_id]
                    lineage = lineage_by_parent.get((prior_tray.physical_tray_id, packet_id))
                    if lineage is None:
                        if candidate_packet.entry != prior_packet.entry:
                            raise EngineeringKTowerRouteConfigurationError(
                                "an untouched same-tray survivor changed its packet entry"
                            )
                        action = KClockRouteAction.PRESERVE_SURVIVOR
                    else:
                        if candidate_packet.entry.weight.representative_particles != (
                            lineage.remainder_weight
                        ):
                            raise EngineeringKTowerRouteConfigurationError(
                                "split remainder weight differs from authenticated lineage"
                            )
                        action = KClockRouteAction.PRESERVE_SPLIT_REMAINDER
                    if candidate_packet.owner_key.vertical_layer_id != clock.vertical_layer_id:
                        packet_k_owners_match = False
                        raise EngineeringKTowerRouteConfigurationError(
                            "same-tray survivor K differs from its preserved K clock"
                        )
                    clocks_after.append(clock)
                    accounted_prior_keys.add((prior_tray.physical_tray_id, packet_id))
                    records.append(
                        KClockRouteRecord(
                            packet_id=packet_id,
                            action=action,
                            source_packet_id=packet_id,
                            source_tray_id=prior_tray.physical_tray_id,
                            destination_tray_id=candidate_tray.physical_tray_id,
                            source_clock_digest=clock.definition_digest,
                            candidate_clock_digest=clock.definition_digest,
                        )
                    )
                    continue

                route = incoming[packet_id]
                # The full schedule initializer binds the initial packet
                # ensemble to the declared nominal K-mesh holdups.  At an
                # evolved route time, real whole-packet K motion and admitted
                # packets make those per-layer holdups dynamic.  Initialize
                # only this authenticated admission; do not reassert the t0
                # layer-holdup equality over the candidate population.
                initialized_clock = (
                    k_scheduler.initialize_engineering_k_admission_clock(
                        packet=candidate_packet,
                        authority=authority,
                        accepted_time_s=event_time_s,
                    )
                )
                if candidate_packet.entry != entry_by_packet[packet_id]:
                    raise EngineeringKTowerRouteConfigurationError(
                        "destination packet entry differs from the authenticated route plan"
                    )
                if candidate_packet.owner_key.vertical_layer_id != (
                    initialized_clock.vertical_layer_id
                ):
                    packet_k_owners_match = False
                    raise EngineeringKTowerRouteConfigurationError(
                        "destination packet K differs from its initialized K clock"
                    )
                if route.authority_class is TransferAuthorityClass.FEED_ADMISSION:
                    action = KClockRouteAction.INITIALIZE_FEED
                    source_packet_id = None
                    source_clock = None
                elif route.authority_class is TransferAuthorityClass.INTER_TRAY:
                    lineage = lineage_by_child.get(packet_id)
                    if lineage is None:
                        source_packet_id = packet_id
                        source_state = schedule_by_tray.get(route.topology.source_cell_id)
                        source_tray = prior_by_tray.get(route.topology.source_cell_id)
                        if source_state is None or source_tray is None:
                            raise EngineeringKTowerRouteConfigurationError(
                                "intertray admission names a source outside the routed tower"
                            )
                        source_packet = _packet_by_id(source_tray).get(packet_id)
                        if (
                            source_packet is None
                            or source_packet.entry != entry_by_packet[packet_id]
                        ):
                            raise EngineeringKTowerRouteConfigurationError(
                                "whole-packet route differs from its authenticated source entry"
                            )
                        source_clock = _clock_by_id(source_state).get(packet_id)
                        action = KClockRouteAction.INITIALIZE_INTERTRAY_WHOLE
                    else:
                        source_packet_id = lineage.parent_packet_id
                        source_state = schedule_by_tray.get(lineage.source_tray_id)
                        if source_state is None:
                            raise EngineeringKTowerRouteConfigurationError(
                                "split-child admission names a source outside the K-state tuple"
                            )
                        source_clock = _clock_by_id(source_state).get(lineage.parent_packet_id)
                        action = KClockRouteAction.INITIALIZE_INTERTRAY_SPLIT_CHILD
                    if source_clock is None:
                        raise EngineeringKTowerRouteConfigurationError(
                            "intertray admission has no authenticated source clock"
                        )
                    accounted_prior_keys.add((route.topology.source_cell_id, source_packet_id))
                else:
                    raise EngineeringKTowerRouteConfigurationError(
                        "a routed tower destination uses a non-admission authority class"
                    )
                clocks_after.append(initialized_clock)
                accounted_route_ids.add(packet_id)
                records.append(
                    KClockRouteRecord(
                        packet_id=packet_id,
                        action=action,
                        source_packet_id=source_packet_id,
                        source_tray_id=route.topology.source_cell_id,
                        destination_tray_id=candidate_tray.physical_tray_id,
                        source_clock_digest=(
                            None if source_clock is None else source_clock.definition_digest
                        ),
                        candidate_clock_digest=initialized_clock.definition_digest,
                    )
                )

            candidate_state = replace(
                prior_state,
                time_s=event_time_s,
                packet_clocks=tuple(sorted(clocks_after, key=lambda clock: clock.packet_id)),
            )
            k_scheduler._validate_state_authority(candidate_state, authority)  # noqa: SLF001
            if set(_clock_by_id(candidate_state)) != set(candidate_packets):
                candidate_packet_sets_match = False
                raise EngineeringKTowerRouteConfigurationError(
                    "reconciled K state and router candidate have different packet sets"
                )
            candidate_states.append(candidate_state)

        # Final sink routes remove source clocks.  A split child never owned a
        # clock; its record explicitly names the preserved source parent.
        for packet_id, route in sorted(route_by_packet.items()):
            if route.topology.destination_cell_id != "external-final-sink":
                continue
            if route.authority_class is not TransferAuthorityClass.DT_TO_DC_HANDOFF:
                raise EngineeringKTowerRouteConfigurationError(
                    "the final sink route has a foreign authority class"
                )
            lineage = lineage_by_child.get(packet_id)
            source_state = schedule_by_tray.get(route.topology.source_cell_id)
            if source_state is None:
                raise EngineeringKTowerRouteConfigurationError(
                    "final sink route names a source outside the scheduler tuple"
                )
            source_clocks = _clock_by_id(source_state)
            if lineage is None:
                source_packet_id = packet_id
                source_clock = source_clocks.get(packet_id)
                action = KClockRouteAction.DROP_FINAL_SINK_WHOLE
                source_tray = prior_by_tray[route.topology.source_cell_id]
                source_packet = _packet_by_id(source_tray).get(packet_id)
                if source_packet is None or source_packet.entry != entry_by_packet[packet_id]:
                    raise EngineeringKTowerRouteConfigurationError(
                        "whole-packet final sink route differs from its authenticated source"
                    )
            else:
                source_packet_id = lineage.parent_packet_id
                source_clock = source_clocks.get(source_packet_id)
                action = KClockRouteAction.DROP_FINAL_SINK_SPLIT_CHILD
            if source_clock is None:
                raise EngineeringKTowerRouteConfigurationError(
                    "final sink route has no authenticated source clock"
                )
            if packet_id in {packet.packet_id for tray in candidates for packet in tray.packets}:
                raise EngineeringKTowerRouteConfigurationError(
                    "a final-sink packet remains in a router candidate tray"
                )
            accounted_route_ids.add(packet_id)
            accounted_prior_keys.add((route.topology.source_cell_id, source_packet_id))
            records.append(
                KClockRouteRecord(
                    packet_id=packet_id,
                    action=action,
                    source_packet_id=source_packet_id,
                    source_tray_id=route.topology.source_cell_id,
                    destination_tray_id="external-final-sink",
                    source_clock_digest=source_clock.definition_digest,
                    candidate_clock_digest=None,
                )
            )

        # Every whole source exit must disappear from its source candidate and
        # every prior clock must either survive or source exactly one route.
        for source_tray in priors:
            source_candidate_ids = set(
                _packet_by_id(candidate_by_tray[source_tray.physical_tray_id])
            )
            for packet in source_tray.packets:
                key = (source_tray.physical_tray_id, packet.packet_id)
                if packet.packet_id in source_candidate_ids:
                    accounted_prior_keys.add(key)
                    continue
                route = route_by_packet.get(packet.packet_id)
                if route is None or route.topology.source_cell_id != source_tray.physical_tray_id:
                    raise EngineeringKTowerRouteConfigurationError(
                        "a removed source packet has no authenticated population route"
                    )
                accounted_prior_keys.add(key)

        all_route_ids = set(route_by_packet)
        every_route_entry_accounted = accounted_route_ids == all_route_ids
        if not every_route_entry_accounted:
            missing = tuple(sorted(all_route_ids - accounted_route_ids))
            raise EngineeringKTowerRouteConfigurationError(
                f"population route entries lack K-clock dispositions: {missing!r}"
            )
        all_prior_keys = {
            (tray.physical_tray_id, packet.packet_id) for tray in priors for packet in tray.packets
        }
        every_prior_clock_accounted = accounted_prior_keys == all_prior_keys
        if not every_prior_clock_accounted:
            raise EngineeringKTowerRouteConfigurationError(
                "one or more pre-route clocks lack an exact routing disposition"
            )

        candidate_tuple = tuple(candidate_states)
        all_candidate_ids = tuple(
            clock.packet_id for state in candidate_tuple for clock in state.packet_clocks
        )
        collision_free = len(all_candidate_ids) == len(set(all_candidate_ids))
        if not collision_free:
            raise EngineeringKTowerRouteConfigurationError(
                "reconciled K states contain a cross-tray packet collision"
            )
        action_counts = {
            action: sum(record.action is action for record in records)
            for action in KClockRouteAction
        }
        ledger = EngineeringKTowerClockLedger(
            prior_clock_count=sum(len(state.packet_clocks) for state in schedules),
            candidate_clock_count=sum(len(state.packet_clocks) for state in candidate_tuple),
            preserved_survivor_count=action_counts[KClockRouteAction.PRESERVE_SURVIVOR],
            preserved_split_remainder_count=action_counts[
                KClockRouteAction.PRESERVE_SPLIT_REMAINDER
            ],
            initialized_feed_count=action_counts[KClockRouteAction.INITIALIZE_FEED],
            initialized_intertray_whole_count=action_counts[
                KClockRouteAction.INITIALIZE_INTERTRAY_WHOLE
            ],
            initialized_intertray_split_child_count=action_counts[
                KClockRouteAction.INITIALIZE_INTERTRAY_SPLIT_CHILD
            ],
            dropped_final_sink_whole_count=action_counts[KClockRouteAction.DROP_FINAL_SINK_WHOLE],
            dropped_final_sink_split_child_count=action_counts[
                KClockRouteAction.DROP_FINAL_SINK_SPLIT_CHILD
            ],
            route_entry_count=len(route_by_packet),
            split_lineage_count=len(routing.split_lineage),
            candidate_packet_sets_match=candidate_packet_sets_match,
            packet_k_owners_match=packet_k_owners_match,
            collision_free=collision_free,
            every_route_entry_accounted=every_route_entry_accounted,
            every_prior_clock_accounted=every_prior_clock_accounted,
        )
        if not ledger.passed:
            raise EngineeringKTowerRouteConfigurationError(
                "whole-tower K-clock routing ledger failed"
            )
        return (
            candidate_tuple,
            tuple(sorted(records, key=lambda record: (record.packet_id, record.action.value))),
            ledger,
        )

    def evaluate(
        self,
        request: EngineeringKTowerRouteReconciliationRequest,
    ) -> (
        ValidatedEngineeringKTowerRouteReconciliation
        | EngineeringKTowerRouteReconciliationRejection
    ):
        """Reconcile a router transaction without mutating accepted state."""

        if type(request) is not EngineeringKTowerRouteReconciliationRequest:
            raise TypeError("reconciler requires an exact reconciliation request")
        try:
            self._router.verify_issued_transaction(request.routing_transaction)
            _validate_router_transaction(request.routing_transaction)
            self._validate_shape_and_priors(request)
            candidates, records, ledger = self._reconcile(request)
            request_digest = _request_digest(request, self._authority_digests)
            transaction_digest = _transaction_digest(
                request_digest=request_digest,
                candidates=candidates,
                records=records,
                ledger=ledger,
            )
            validated = ValidatedEngineeringKTowerRouteReconciliation(
                request=request,
                candidate_schedule_states=candidates,
                records=records,
                ledger=ledger,
                authority_definition_digests=self._authority_digests,
                request_digest=request_digest,
                transaction_digest=transaction_digest,
                replay_token=_digest(
                    domain="engineering-k-tower-route-replay",
                    parts=(transaction_digest,),
                ),
                _reconciler_token=self._reconciler_token,
                _seal=_VALIDATED_RECONCILIATION_SEAL,
            )
            self._issued_transactions[id(validated)] = validated
            return validated
        except (
            EngineeringKTowerRouteReconciliationError,
            tower_router.KCellTowerRouterError,
            k_scheduler.EngineeringKEventSchedulerError,
            TypeError,
            ValueError,
            OverflowError,
        ) as error:
            return EngineeringKTowerRouteReconciliationRejection(
                attempt_id=request.attempt_id,
                reason=f"{type(error).__name__}: {error}",
                rollback_schedule_states=request.prior_schedule_states,
            )

    @staticmethod
    def _prior_identity(
        transaction: ValidatedEngineeringKTowerRouteReconciliation,
    ) -> str:
        routing_priors = transaction.request.routing_transaction.request.prior_trays
        return _digest(
            domain="engineering-k-tower-route-prior",
            parts=(
                *(tray.state_digest for tray in routing_priors),
                *(state.definition_digest for state in transaction.request.prior_schedule_states),
            ),
        )

    def verify_issued_transaction(
        self,
        transaction: ValidatedEngineeringKTowerRouteReconciliation,
    ) -> ValidatedEngineeringKTowerRouteReconciliation:
        """Verify exact issuer, deterministic contents, and unused commit authority."""

        if (
            type(transaction) is not ValidatedEngineeringKTowerRouteReconciliation
            or transaction._seal is not _VALIDATED_RECONCILIATION_SEAL
        ):
            raise EngineeringKTowerRouteReconciliationError(
                "verification requires a reconciler-validated transaction"
            )
        if transaction._reconciler_token is not self._reconciler_token:
            raise EngineeringKTowerRouteReconciliationError(
                "reconciliation transaction was issued by another reconciler"
            )
        if self._issued_transactions.get(id(transaction)) is not transaction:
            raise EngineeringKTowerRouteReconciliationError(
                "reconciliation transaction is not the exact evaluator-issued in-memory object"
            )
        if transaction.authority_definition_digests != self._authority_digests:
            raise EngineeringKTowerRouteReconciliationError(
                "destination K authority pins changed after evaluation"
            )
        request = transaction.request
        self._router.verify_issued_transaction(request.routing_transaction)
        _validate_router_transaction(request.routing_transaction)
        self._validate_shape_and_priors(request)
        candidates, records, ledger = self._reconcile(request)
        expected_request = _request_digest(request, self._authority_digests)
        if expected_request != transaction.request_digest:
            raise EngineeringKTowerRouteReconciliationError(
                "reconciliation request identity changed after evaluation"
            )
        expected_transaction = _transaction_digest(
            request_digest=expected_request,
            candidates=candidates,
            records=records,
            ledger=ledger,
        )
        if (
            expected_transaction != transaction.transaction_digest
            or candidates != transaction.candidate_schedule_states
            or records != transaction.records
            or ledger != transaction.ledger
        ):
            raise EngineeringKTowerRouteReconciliationError(
                "reconciliation transaction changed after evaluation"
            )
        expected_replay = _digest(
            domain="engineering-k-tower-route-replay",
            parts=(expected_transaction,),
        )
        if expected_replay != transaction.replay_token:
            raise EngineeringKTowerRouteReconciliationError(
                "reconciliation replay identity changed after evaluation"
            )
        if transaction.replay_token in self._consumed_replay_tokens:
            raise ReplayedEngineeringKTowerRouteError(
                "reconciliation transaction was already consumed"
            )
        if self._prior_identity(transaction) in self._consumed_prior_identities:
            raise ReplayedEngineeringKTowerRouteError(
                "another reconciliation already consumed this accepted prior tuple"
            )
        return transaction

    def verify_commit_available(
        self,
        current_prior_trays: tuple[tray_host.AcceptedKCellTrayState, ...],
        current_schedule_states: tuple[k_scheduler.EngineeringKScheduleState, ...],
        transaction: ValidatedEngineeringKTowerRouteReconciliation,
    ) -> ValidatedEngineeringKTowerRouteReconciliation:
        """Verify a reconciliation can commit without consuming authority state."""

        self.verify_issued_transaction(transaction)
        routing_priors = transaction.request.routing_transaction.request.prior_trays
        if (
            type(current_prior_trays) is not tuple
            or len(current_prior_trays) != len(routing_priors)
            or any(
                current is not prior
                for current, prior in zip(current_prior_trays, routing_priors, strict=True)
            )
        ):
            raise StaleEngineeringKTowerRouteStateError(
                "current router prior tuple is stale or foreign"
            )
        schedule_priors = transaction.request.prior_schedule_states
        if (
            type(current_schedule_states) is not tuple
            or len(current_schedule_states) != len(schedule_priors)
            or any(
                current is not prior
                for current, prior in zip(current_schedule_states, schedule_priors, strict=True)
            )
        ):
            raise StaleEngineeringKTowerRouteStateError(
                "current K schedule tuple is stale or foreign"
            )
        return transaction

    def commit(
        self,
        current_prior_trays: tuple[tray_host.AcceptedKCellTrayState, ...],
        current_schedule_states: tuple[k_scheduler.EngineeringKScheduleState, ...],
        transaction: ValidatedEngineeringKTowerRouteReconciliation,
    ) -> tuple[k_scheduler.EngineeringKScheduleState, ...]:
        """Consume one reconciliation after deterministic stale/tamper replay."""

        self.verify_commit_available(current_prior_trays, current_schedule_states, transaction)
        request = transaction.request
        prior_identity = self._prior_identity(transaction)
        replayed = self.evaluate(request)
        if type(replayed) is EngineeringKTowerRouteReconciliationRejection:
            raise EngineeringKTowerRouteReconciliationError(
                f"reconciliation failed deterministic commit replay: {replayed.reason}"
            )
        if type(replayed) is not ValidatedEngineeringKTowerRouteReconciliation:
            raise EngineeringKTowerRouteReconciliationError(
                "reconciliation replay returned a foreign result type"
            )
        if (
            replayed.request_digest != transaction.request_digest
            or replayed.transaction_digest != transaction.transaction_digest
            or replayed.replay_token != transaction.replay_token
            or replayed.candidate_schedule_states != transaction.candidate_schedule_states
            or replayed.records != transaction.records
            or replayed.ledger != transaction.ledger
        ):
            raise EngineeringKTowerRouteReconciliationError(
                "reconciliation transaction differs from deterministic commit replay"
            )
        if not transaction.ledger.passed:
            raise EngineeringKTowerRouteReconciliationError(
                "reconciliation ledger no longer passes"
            )
        self._consumed_replay_tokens.add(transaction.replay_token)
        self._consumed_prior_identities.add(prior_identity)
        return transaction.candidate_schedule_states


__all__ = (
    "EngineeringKTowerClockLedger",
    "EngineeringKTowerRouteConfigurationError",
    "EngineeringKTowerRouteReconciliationError",
    "EngineeringKTowerRouteReconciliationRejection",
    "EngineeringKTowerRouteReconciliationRequest",
    "EngineeringKTowerRouteReconciler",
    "KClockRouteAction",
    "KClockRouteRecord",
    "RECONCILER_SCHEMA_ID",
    "RECONCILER_SOURCE_IDENTITY",
    "ReplayedEngineeringKTowerRouteError",
    "StaleEngineeringKTowerRouteStateError",
    "ValidatedEngineeringKTowerRouteReconciliation",
)
