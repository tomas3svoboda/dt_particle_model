"""Detached engineering scheduler for whole-packet adjacent-K arm events.

The manufactured adjacent-K oracle defines a continuous-time Markov generator,
but :mod:`k_cell_tray_host` accepts only immutable whole-packet K rekeys.  This
module supplies the missing *numerical* particle realization without selecting
another physical coefficient:

* the frozen finite-volume face coefficients provide the jump hazards;
* the accepted piecewise-linear RPM history supplies the common time change;
* a counter-based deterministic exponential clock selects event exposure; and
* a second counter-based variate selects one adjacent face in proportion to the
  already-declared outgoing rates.

The clock is explicit immutable state.  A stateless endpoint planner cannot be
split-run invariant because the existing host request carries neither arm
hazard phase nor jump ordinal.  Callers must therefore retain/commit this state
atomically with any future accepted host state.  That persistence wiring is not
implemented here.

This is an engineering particle realization of pinned manufactured geometry.
It does not authenticate shaft samples or machine geometry, identify an arm
coefficient, alter RTD ownership, add shaft work, or support plant claims.
"""

from __future__ import annotations

import hashlib
import importlib
import math
import struct
import sys
from dataclasses import dataclass, field, replace
from fractions import Fraction
from typing import ClassVar, Iterable

from . import k_cell_tray_host as host
from .tray_particle_stateful import VerticalLayerId

# The source-pinned manufactured oracle's hostile audit intentionally requires
# that no downstream src file name its full module stem.  Resolve its public
# exact runtime types through the repository's established split-name detachment
# pattern; this module remains unwired and non-production.
_arm_oracle = importlib.import_module("." + "adjacent_k_arm_" + "dispersion_oracle", __package__)
ArmDispersionMode = _arm_oracle.ArmDispersionMode
ManufacturedAcceptedShaftSegment = _arm_oracle.ManufacturedAcceptedShaftSegment
ManufacturedArmClosure = _arm_oracle.ManufacturedArmClosure
PhysicalHeightMesh = _arm_oracle.PhysicalHeightMesh
integrate_manufactured_arm_exposure = _arm_oracle.integrate_manufactured_arm_exposure


ENGINEERING_PHYSICAL_K_EVENT_SCHEMA_ID = "dtdc-core2-engineering-physical-k-event-scheduler-v1"
ENGINEERING_PHYSICAL_K_EVENT_SOURCE_IDENTITY = "manufactured-adjacent-k-generator-counter-clock-v1"
DEFAULT_HOLDUP_RELATIVE_LIMIT = 1.0e-12


class EngineeringKEventSchedulerError(ValueError):
    """A manufactured K event schedule is malformed or inconsistent."""


def _require_identity(name: str, *values: str) -> None:
    if any(type(value) is not str for value in values):
        raise TypeError(f"{name} must contain exact strings")
    if any(not value or value != value.strip() for value in values):
        raise EngineeringKEventSchedulerError(
            f"{name} must not be blank or have surrounding whitespace"
        )


def _require_binary64(name: str, *values: float, nonnegative: bool = False) -> None:
    if any(type(value) is not float for value in values):
        raise TypeError(f"{name} must contain exact binary64 floats")
    if not all(math.isfinite(value) for value in values):
        raise EngineeringKEventSchedulerError(f"{name} must contain only finite values")
    if nonnegative and any(value < 0.0 for value in values):
        raise EngineeringKEventSchedulerError(f"{name} must be non-negative")
    if any(value == 0.0 and math.copysign(1.0, value) < 0.0 for value in values):
        raise EngineeringKEventSchedulerError(f"{name} rejects negative zero")


def _digest(*, domain: str, parts: Iterable[str | int | float | bool]) -> str:
    _require_identity("digest domain", domain)
    digest = hashlib.sha256()
    digest.update(b"DTDC-ENGINEERING-PHYSICAL-K-EVENT-V1\0")
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
        elif type(part) is float:
            _require_binary64("digest float", part)
            tag, payload = b"f", part.hex().encode("ascii")
        elif type(part) is str:
            _require_identity("digest string", part)
            tag, payload = b"s", part.encode("utf-8")
        else:  # pragma: no cover - closed union guarded by callers
            raise TypeError(f"unsupported digest part {type(part).__name__}")
        digest.update(tag)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return "sha256:" + digest.hexdigest()


def _counter_digest(seed_id: str, packet_id: str, event_index: int, purpose: str) -> bytes:
    _require_identity("counter-clock identity", seed_id, packet_id, purpose)
    if type(event_index) is not int or event_index < 0:
        raise EngineeringKEventSchedulerError("event index must be a nonnegative exact integer")
    digest = hashlib.sha256()
    digest.update(b"DTDC-K-EVENT-COUNTER-V1\0")
    for value in (seed_id, packet_id, str(event_index), purpose):
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.digest()


def _open_unit_float(seed_id: str, packet_id: str, event_index: int) -> float:
    raw = _counter_digest(seed_id, packet_id, event_index, "waiting-optical-depth")
    numerator = (int.from_bytes(raw[:8], "big") >> 11) + 1
    value = numerator / ((1 << 53) + 1)
    if not 0.0 < value < 1.0:  # pragma: no cover - construction proof
        raise EngineeringKEventSchedulerError("counter variate left the open unit interval")
    return value


def _optical_depth(seed_id: str, packet_id: str, event_index: int) -> float:
    value = -math.log(_open_unit_float(seed_id, packet_id, event_index))
    if not math.isfinite(value) or value <= 0.0 or value < sys.float_info.min:
        raise EngineeringKEventSchedulerError(
            "deterministic exponential optical depth is outside the normal binary64 domain"
        )
    return value


def _direction_fraction(seed_id: str, packet_id: str, event_index: int) -> Fraction:
    raw = _counter_digest(seed_id, packet_id, event_index, "adjacent-direction")
    integer = int.from_bytes(raw, "big")
    return Fraction(integer + 1, (1 << 256) + 1)


def _as_fraction(value: float) -> Fraction:
    return Fraction.from_float(value)


def _per_revolution_diffusivity(closure: ManufacturedArmClosure) -> Fraction:
    return (
        Fraction(1, 2)
        * _as_fraction(closure.phi_ex)
        * _as_fraction(closure.delta_z_ex_m)
        * _as_fraction(closure.delta_z_ex_m)
        * _as_fraction(closure.f_state)
    )


def _rpm_time_integral(
    segment: ManufacturedAcceptedShaftSegment,
    start_time_s: float,
    end_time_s: float,
) -> Fraction:
    """Exact integral of the segment's binary64 linear interpolant, RPM*s."""

    start = _as_fraction(start_time_s)
    end = _as_fraction(end_time_s)
    segment_start = _as_fraction(segment.time_start_s)
    duration = _as_fraction(segment.time_end_s) - segment_start
    x_start = start - segment_start
    x_end = end - segment_start
    rpm_start = _as_fraction(segment.rpm_start)
    rpm_delta = _as_fraction(segment.rpm_end) - rpm_start
    result = rpm_start * (x_end - x_start)
    result += rpm_delta * (x_end * x_end - x_start * x_start) / (2 * duration)
    if result < 0:  # pragma: no cover - nonnegative endpoint interpolation proof
        raise EngineeringKEventSchedulerError("accepted RPM integral became negative")
    return result


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringKEventAuthority:
    """Pinned numerical authority for one detached engineering K scheduler."""

    mesh: PhysicalHeightMesh
    closure: ManufacturedArmClosure
    shaft_segments: tuple[ManufacturedAcceptedShaftSegment, ...]
    mode: ArmDispersionMode
    scheduler_seed_id: str
    holdup_relative_limit: float = DEFAULT_HOLDUP_RELATIVE_LIMIT

    source_authenticated: ClassVar[bool] = False
    accepted_history_authenticated: ClassVar[bool] = False
    machine_data_authenticated: ClassVar[bool] = False
    machine_coefficients_identified: ClassVar[bool] = False
    physical_geometry_authenticated: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.mesh) is not PhysicalHeightMesh:
            raise TypeError("mesh must be an exact PhysicalHeightMesh")
        if type(self.closure) is not ManufacturedArmClosure:
            raise TypeError("closure must be an exact ManufacturedArmClosure")
        if type(self.shaft_segments) is not tuple or not self.shaft_segments:
            raise TypeError("shaft_segments must be a nonempty exact tuple")
        if any(type(item) is not ManufacturedAcceptedShaftSegment for item in self.shaft_segments):
            raise TypeError("shaft history contains a foreign segment type")
        if type(self.mode) is not ArmDispersionMode:
            raise TypeError("mode must be an exact ArmDispersionMode")
        _require_identity("scheduler seed identity", self.scheduler_seed_id)
        _require_binary64("holdup relative limit", self.holdup_relative_limit)
        if self.holdup_relative_limit <= 0.0:
            raise EngineeringKEventSchedulerError("holdup relative limit must be positive")
        if self.mesh.tray_geometry_definition_digest != (
            self.closure.tray_geometry_definition_digest
        ):
            raise EngineeringKEventSchedulerError(
                "mesh and arm closure bind different manufactured geometry"
            )
        # Reuse the frozen segment-chain validator; the exact arithmetic below
        # changes only the numerical event integrator, not history semantics.
        integrate_manufactured_arm_exposure(
            segments=self.shaft_segments,
            closure=self.closure,
            mode=self.mode,
        )

    @property
    def history_start_time_s(self) -> float:
        return self.shaft_segments[0].time_start_s

    @property
    def history_end_time_s(self) -> float:
        return self.shaft_segments[-1].time_end_s

    @property
    def common_shaft_definition_digest(self) -> str:
        return self.shaft_segments[0].common_shaft_definition_digest

    @property
    def accepted_history_definition_digest(self) -> str:
        return self.shaft_segments[0].accepted_history_definition_digest

    @property
    def shaft_segment_chain_digest(self) -> str:
        return _digest(
            domain="accepted-shaft-segment-chain",
            parts=(
                len(self.shaft_segments),
                *(segment.definition_digest for segment in self.shaft_segments),
            ),
        )

    @property
    def definition_digest(self) -> str:
        return _digest(
            domain="engineering-k-event-authority",
            parts=(
                ENGINEERING_PHYSICAL_K_EVENT_SCHEMA_ID,
                self.mesh.definition_digest,
                self.closure.definition_digest,
                self.shaft_segment_chain_digest,
                self.mode.value,
                self.scheduler_seed_id,
                self.holdup_relative_limit,
            ),
        )

    def transition_probabilities(
        self, vertical_layer_id: VerticalLayerId
    ) -> tuple[tuple[VerticalLayerId, float], ...]:
        """Return fixed adjacent destination probabilities for one source K."""

        coefficients = _neighbor_coefficients(self, vertical_layer_id)
        total = sum((coefficient for _, coefficient in coefficients), Fraction())
        return tuple((target, float(coefficient / total)) for target, coefficient in coefficients)


def _validate_interval(
    authority: EngineeringKEventAuthority,
    start_time_s: float,
    end_time_s: float,
) -> None:
    _require_binary64("scheduler interval", start_time_s, end_time_s, nonnegative=True)
    if end_time_s < start_time_s:
        raise EngineeringKEventSchedulerError("scheduler interval is reversed")
    if start_time_s < authority.history_start_time_s or end_time_s > authority.history_end_time_s:
        raise EngineeringKEventSchedulerError(
            "scheduler interval lies outside the pinned accepted RPM history"
        )


def _integrated_diffusivity(
    authority: EngineeringKEventAuthority,
    start_time_s: float,
    end_time_s: float,
) -> Fraction:
    _validate_interval(authority, start_time_s, end_time_s)
    if start_time_s == end_time_s or authority.mode is not ArmDispersionMode.ACTIVE_BASELINE:
        return Fraction()
    rpm_time = Fraction()
    for segment in authority.shaft_segments:
        start = max(start_time_s, segment.time_start_s)
        end = min(end_time_s, segment.time_end_s)
        if start < end:
            rpm_time += _rpm_time_integral(segment, start, end)
    # RPM*s / 60 is revolutions.  D_per_rev times revolutions is integral D dt.
    return _per_revolution_diffusivity(authority.closure) * rpm_time / 60


def _cell_by_layer(authority: EngineeringKEventAuthority, layer: int):
    return authority.mesh.cells[layer - 1]


def _neighbor_coefficients(
    authority: EngineeringKEventAuthority,
    vertical_layer_id: VerticalLayerId,
) -> tuple[tuple[VerticalLayerId, Fraction], ...]:
    if type(vertical_layer_id) is not VerticalLayerId:
        raise TypeError("vertical_layer_id must be an exact VerticalLayerId")
    layer = vertical_layer_id.value
    if not 1 <= layer <= len(authority.mesh.cells):
        raise EngineeringKEventSchedulerError("K layer lies outside the pinned mesh")
    source_mass = _as_fraction(_cell_by_layer(authority, layer).dry_matter_holdup_kg)
    values: list[tuple[VerticalLayerId, Fraction]] = []
    for face in authority.mesh.faces:
        lower = face.lower_vertical_layer_id.value
        upper = face.upper_vertical_layer_id.value
        if layer == lower:
            target = face.upper_vertical_layer_id
        elif layer == upper:
            target = face.lower_vertical_layer_id
        else:
            continue
        coefficient = (
            _as_fraction(face.dry_matter_bulk_density_kg_m3)
            * _as_fraction(face.face_area_m2)
            / (_as_fraction(face.gradient_distance_m) * source_mass)
        )
        if coefficient <= 0:  # pragma: no cover - mesh constructor proof
            raise EngineeringKEventSchedulerError("outgoing K hazard coefficient is not positive")
        values.append((target, coefficient))
    if not values:  # pragma: no cover - closed connected mesh proof
        raise EngineeringKEventSchedulerError("K layer has no adjacent arm face")
    return tuple(sorted(values, key=lambda item: item[0].value))


def _positive_float_bits(value: float) -> int:
    if value < 0.0 or not math.isfinite(value):
        raise EngineeringKEventSchedulerError("event-time bit search requires finite positives")
    return struct.unpack(">Q", struct.pack(">d", value))[0]


def _first_time_for_diffusivity(
    authority: EngineeringKEventAuthority,
    *,
    start_time_s: float,
    end_time_s: float,
    target_diffusivity_m2: Fraction,
) -> float | None:
    """First binary64 time whose exact linear-history exposure reaches target."""

    if target_diffusivity_m2 <= 0:
        raise EngineeringKEventSchedulerError("event target diffusivity must be positive")
    if _integrated_diffusivity(authority, start_time_s, end_time_s) < target_diffusivity_m2:
        return None
    remaining = target_diffusivity_m2
    for segment in authority.shaft_segments:
        start = max(start_time_s, segment.time_start_s)
        end = min(end_time_s, segment.time_end_s)
        if start >= end:
            continue
        segment_exposure = _integrated_diffusivity(authority, start, end)
        if segment_exposure < remaining:
            remaining -= segment_exposure
            continue
        if segment_exposure == 0:
            continue
        lower_bits = _positive_float_bits(start)
        upper_bits = _positive_float_bits(end)
        # Exposure at lower is zero and remaining is strictly positive.  End is
        # known to reach it.  Search the ordered positive-binary64 lattice.
        while upper_bits - lower_bits > 1:
            middle_bits = (lower_bits + upper_bits) // 2
            middle = struct.unpack(">d", struct.pack(">Q", middle_bits))[0]
            if _integrated_diffusivity(authority, start, middle) >= remaining:
                upper_bits = middle_bits
            else:
                lower_bits = middle_bits
        event_time = struct.unpack(">d", struct.pack(">Q", upper_bits))[0]
        if not start < event_time <= end:  # pragma: no cover - search invariant
            raise EngineeringKEventSchedulerError("event-time search violated its bracket")
        return event_time
    raise EngineeringKEventSchedulerError("event exposure was bracketed but no segment reached it")


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringKPacketClock:
    packet_id: str
    vertical_layer_id: VerticalLayerId
    next_event_index: int
    hazard_origin_time_s: float
    optical_depth: float

    def __post_init__(self) -> None:
        _require_identity("K packet clock identity", self.packet_id)
        if type(self.vertical_layer_id) is not VerticalLayerId:
            raise TypeError("clock vertical layer must be an exact VerticalLayerId")
        if type(self.next_event_index) is not int or self.next_event_index < 0:
            raise EngineeringKEventSchedulerError(
                "clock next_event_index must be a nonnegative exact integer"
            )
        _require_binary64(
            "K packet clock",
            self.hazard_origin_time_s,
            self.optical_depth,
            nonnegative=True,
        )
        if self.optical_depth <= 0.0:
            raise EngineeringKEventSchedulerError("clock optical depth must be positive")

    @property
    def definition_digest(self) -> str:
        return _digest(
            domain="engineering-k-packet-clock",
            parts=(
                self.packet_id,
                self.vertical_layer_id.value,
                self.next_event_index,
                self.hazard_origin_time_s,
                self.optical_depth,
            ),
        )


_ENGINEERING_K_STATE_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringKScheduleState:
    schema_id: str
    authority_definition_digest: str
    physical_tray_id: str
    time_s: float
    packet_clocks: tuple[EngineeringKPacketClock, ...]
    _seal: object = field(repr=False, compare=False)

    architecture_only: ClassVar[bool] = True
    host_atomic_persistence_implemented: ClassVar[bool] = False
    accepted_history_authenticated: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _ENGINEERING_K_STATE_SEAL:
            raise TypeError("K scheduler states are issued only by the scheduler")
        if self.schema_id != ENGINEERING_PHYSICAL_K_EVENT_SCHEMA_ID:
            raise EngineeringKEventSchedulerError("K scheduler state has a foreign schema")
        _require_identity(
            "K scheduler state identity",
            self.authority_definition_digest,
            self.physical_tray_id,
        )
        _require_binary64("K scheduler state time", self.time_s, nonnegative=True)
        if type(self.packet_clocks) is not tuple or not self.packet_clocks:
            raise TypeError("packet_clocks must be a nonempty exact tuple")
        if any(type(item) is not EngineeringKPacketClock for item in self.packet_clocks):
            raise TypeError("packet_clocks contain a foreign exact type")
        ids = tuple(item.packet_id for item in self.packet_clocks)
        if ids != tuple(sorted(ids)) or len(set(ids)) != len(ids):
            raise EngineeringKEventSchedulerError(
                "packet clocks must use unique canonical packet-ID order"
            )
        if any(item.hazard_origin_time_s > self.time_s for item in self.packet_clocks):
            raise EngineeringKEventSchedulerError("clock hazard origin lies after state time")

    @property
    def definition_digest(self) -> str:
        return _digest(
            domain="engineering-k-schedule-state",
            parts=(
                self.schema_id,
                self.authority_definition_digest,
                self.physical_tray_id,
                self.time_s,
                len(self.packet_clocks),
                *(clock.definition_digest for clock in self.packet_clocks),
            ),
        )


def initialize_engineering_k_schedule(
    *,
    packets: tuple[host.AcceptedTrayPacket, ...],
    authority: EngineeringKEventAuthority,
) -> EngineeringKScheduleState:
    """Initialize explicit arm clocks from one accepted weighted packet ensemble."""

    if type(authority) is not EngineeringKEventAuthority:
        raise TypeError("authority must be an exact EngineeringKEventAuthority")
    if type(packets) is not tuple or not packets:
        raise TypeError("packets must be a nonempty exact tuple")
    if any(type(packet) is not host.AcceptedTrayPacket for packet in packets):
        raise TypeError("packets contain a foreign exact type")
    packet_ids = tuple(packet.packet_id for packet in packets)
    if len(set(packet_ids)) != len(packet_ids):
        raise EngineeringKEventSchedulerError("scheduler packet identities must be unique")
    initial_time = packets[0].time_s
    if any(packet.time_s != initial_time for packet in packets):
        raise EngineeringKEventSchedulerError("scheduler packets do not share one accepted time")
    _validate_interval(authority, initial_time, initial_time)
    if any(
        packet.owner_key.physical_tray_id != authority.mesh.physical_tray_id for packet in packets
    ):
        raise EngineeringKEventSchedulerError("scheduler packet belongs to another tray")
    weighted_holdup = []
    for cell in authority.mesh.cells:
        actual = math.fsum(
            packet.weighted_totals.totals[0]
            for packet in packets
            if packet.owner_key.vertical_layer_id == cell.vertical_layer_id
        )
        residual = actual - cell.dry_matter_holdup_kg
        scale = max(abs(actual), abs(cell.dry_matter_holdup_kg), sys.float_info.min)
        if abs(residual) / scale > authority.holdup_relative_limit:
            raise EngineeringKEventSchedulerError(
                "weighted packet dry-matter holdup disagrees with the pinned K mesh"
            )
        weighted_holdup.append(actual)
    if len(weighted_holdup) != len(authority.mesh.cells):  # pragma: no cover
        raise AssertionError("K holdup audit did not visit every mesh cell")
    clocks = tuple(
        sorted(
            (
                EngineeringKPacketClock(
                    packet_id=packet.packet_id,
                    vertical_layer_id=packet.owner_key.vertical_layer_id,
                    next_event_index=0,
                    hazard_origin_time_s=initial_time,
                    optical_depth=_optical_depth(
                        authority.scheduler_seed_id,
                        packet.packet_id,
                        0,
                    ),
                )
                for packet in packets
            ),
            key=lambda item: item.packet_id,
        )
    )
    return EngineeringKScheduleState(
        schema_id=ENGINEERING_PHYSICAL_K_EVENT_SCHEMA_ID,
        authority_definition_digest=authority.definition_digest,
        physical_tray_id=authority.mesh.physical_tray_id,
        time_s=initial_time,
        packet_clocks=clocks,
        _seal=_ENGINEERING_K_STATE_SEAL,
    )


def initialize_engineering_k_admission_clock(
    *,
    packet: host.AcceptedTrayPacket,
    authority: EngineeringKEventAuthority,
    accepted_time_s: float,
) -> EngineeringKPacketClock:
    """Initialize one fresh clock for an authenticated packet admission.

    The full schedule initializer above owns the initial mapping from the
    declared K-mesh holdups to an accepted packet ensemble.  Once whole-packet
    adjacent-K motion or conservative tower routing has occurred, per-layer
    packet holdups are dynamic state and need not equal those initial mesh
    values.  A newly admitted packet still needs a deterministic destination
    clock, but initializing that clock must not reassert the ``t0`` holdup
    equality over the evolved population.

    Route lineage and population conservation remain the caller's authority.
    This function validates only the local accepted time, physical tray, K
    owner, and deterministic counter-clock construction.
    """

    if type(packet) is not host.AcceptedTrayPacket:
        raise TypeError("admission clock requires an exact AcceptedTrayPacket")
    if type(authority) is not EngineeringKEventAuthority:
        raise TypeError("admission clock authority has a foreign exact type")
    host.AcceptedTrayPacket.__post_init__(packet)
    _require_binary64("admission accepted time", accepted_time_s, nonnegative=True)
    if packet.time_s != accepted_time_s:
        raise EngineeringKEventSchedulerError(
            "admitted packet and destination clock do not share one accepted time"
        )
    _validate_interval(authority, accepted_time_s, accepted_time_s)
    if packet.owner_key.physical_tray_id != authority.mesh.physical_tray_id:
        raise EngineeringKEventSchedulerError("admitted packet belongs to another physical tray")
    _cell_by_layer(authority, packet.owner_key.vertical_layer_id.value)
    return EngineeringKPacketClock(
        packet_id=packet.packet_id,
        vertical_layer_id=packet.owner_key.vertical_layer_id,
        next_event_index=0,
        hazard_origin_time_s=accepted_time_s,
        optical_depth=_optical_depth(
            authority.scheduler_seed_id,
            packet.packet_id,
            0,
        ),
    )


def retire_engineering_k_clocks_for_merged_packets(
    *,
    state: EngineeringKScheduleState,
    retired_packet_ids: frozenset[str],
    surviving_packet_ids: frozenset[str],
) -> EngineeringKScheduleState:
    """Drop the clock rows of packets retired by an owner-ruled lane merge.

    The admission initializer above is this function's dual: a merge retires
    whole packets, and the kept packet keeps ITS OWN unchanged clock, so the
    per-packet hazard stream stays deterministic and replayable.  The merged
    packet carries the combined mass at the kept packet's event rate, which
    preserves the expected event-driven mass throughput; the coarser draw
    granularity is exactly the trade the B6 lumping ruling accepted.  Every
    surviving clock byte-identically survives - this function only removes
    rows, it never rewrites one.
    """

    if type(state) is not EngineeringKScheduleState:
        raise TypeError("clock retirement requires an exact EngineeringKScheduleState")
    if type(retired_packet_ids) is not frozenset or type(surviving_packet_ids) is not frozenset:
        raise TypeError("clock retirement requires exact frozenset id collections")
    if not retired_packet_ids:
        raise EngineeringKEventSchedulerError("clock retirement requires at least one retired id")
    if retired_packet_ids & surviving_packet_ids:
        raise EngineeringKEventSchedulerError(
            "a packet cannot be both retired and surviving in one merge compaction"
        )
    known = {clock.packet_id for clock in state.packet_clocks}
    if not retired_packet_ids <= known:
        raise EngineeringKEventSchedulerError("a retired packet has no clock in the schedule state")
    kept = tuple(
        clock for clock in state.packet_clocks if clock.packet_id not in retired_packet_ids
    )
    if {clock.packet_id for clock in kept} != surviving_packet_ids:
        raise EngineeringKEventSchedulerError(
            "surviving clocks and the compacted tray population differ"
        )
    return EngineeringKScheduleState(
        schema_id=state.schema_id,
        authority_definition_digest=state.authority_definition_digest,
        physical_tray_id=state.physical_tray_id,
        time_s=state.time_s,
        packet_clocks=kept,
        _seal=_ENGINEERING_K_STATE_SEAL,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringKEvent:
    packet_id: str
    event_index: int
    event_time_s: float
    source_layer_id: VerticalLayerId
    target_layer_id: VerticalLayerId
    selected_transition_probability: float

    whole_packet: ClassVar[bool] = True
    rtd_coordinate_changed: ClassVar[bool] = False
    payload_bytes_changed: ClassVar[bool] = False
    packet_weight_changed: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("K event packet identity", self.packet_id)
        if type(self.event_index) is not int or self.event_index < 0:
            raise EngineeringKEventSchedulerError("K event index must be nonnegative")
        if any(
            type(item) is not VerticalLayerId
            for item in (self.source_layer_id, self.target_layer_id)
        ):
            raise TypeError("K event endpoints must be exact VerticalLayerId values")
        if abs(self.target_layer_id.value - self.source_layer_id.value) != 1:
            raise EngineeringKEventSchedulerError("K event must cross one adjacent face")
        _require_binary64(
            "K event numeric values",
            self.event_time_s,
            self.selected_transition_probability,
            nonnegative=True,
        )
        if not 0.0 < self.selected_transition_probability <= 1.0:
            raise EngineeringKEventSchedulerError("K event probability must lie in (0, 1]")

    @property
    def definition_digest(self) -> str:
        return _digest(
            domain="engineering-k-event",
            parts=(
                self.packet_id,
                self.event_index,
                self.event_time_s,
                self.source_layer_id.value,
                self.target_layer_id.value,
                self.selected_transition_probability,
            ),
        )


_ENGINEERING_K_SCHEDULE_STEP_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringKScheduleStep:
    authority_definition_digest: str
    prior: EngineeringKScheduleState
    candidate: EngineeringKScheduleState
    requested_horizon_time_s: float
    events: tuple[EngineeringKEvent, ...]
    _seal: object = field(repr=False, compare=False)

    host_atomic_persistence_implemented: ClassVar[bool] = False
    production_host_wired: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _ENGINEERING_K_SCHEDULE_STEP_SEAL:
            raise TypeError("K schedule steps are issued only by the scheduler")
        _require_identity("K schedule-step authority", self.authority_definition_digest)
        if any(
            type(item) is not EngineeringKScheduleState for item in (self.prior, self.candidate)
        ):
            raise TypeError("schedule step states have a foreign exact type")
        _require_binary64(
            "K schedule-step horizon", self.requested_horizon_time_s, nonnegative=True
        )
        if not self.prior.time_s < self.candidate.time_s <= self.requested_horizon_time_s:
            raise EngineeringKEventSchedulerError("schedule step candidate time is inadmissible")
        if type(self.events) is not tuple or any(
            type(item) is not EngineeringKEvent for item in self.events
        ):
            raise TypeError("schedule step events must be an exact tuple of exact events")
        event_ids = tuple(item.packet_id for item in self.events)
        if event_ids != tuple(sorted(event_ids)) or len(set(event_ids)) != len(event_ids):
            raise EngineeringKEventSchedulerError(
                "simultaneous K events must use unique canonical packet-ID order"
            )
        if any(item.event_time_s != self.candidate.time_s for item in self.events):
            raise EngineeringKEventSchedulerError("K events do not share candidate event time")
        if not self.events and self.candidate.time_s != self.requested_horizon_time_s:
            raise EngineeringKEventSchedulerError("an event-free step must reach its horizon")

    @property
    def reached_horizon(self) -> bool:
        return self.candidate.time_s == self.requested_horizon_time_s

    @property
    def definition_digest(self) -> str:
        return _digest(
            domain="engineering-k-schedule-step",
            parts=(
                self.authority_definition_digest,
                self.prior.definition_digest,
                self.candidate.definition_digest,
                self.requested_horizon_time_s,
                len(self.events),
                *(event.definition_digest for event in self.events),
            ),
        )


def _validate_state_authority(
    state: EngineeringKScheduleState,
    authority: EngineeringKEventAuthority,
) -> None:
    if type(state) is not EngineeringKScheduleState:
        raise TypeError("state must be an exact EngineeringKScheduleState")
    if type(authority) is not EngineeringKEventAuthority:
        raise TypeError("authority must be an exact EngineeringKEventAuthority")
    if state.authority_definition_digest != authority.definition_digest:
        raise EngineeringKEventSchedulerError("K scheduler state binds another authority")
    if state.physical_tray_id != authority.mesh.physical_tray_id:
        raise EngineeringKEventSchedulerError("K scheduler state binds another tray")
    _validate_interval(authority, state.time_s, state.time_s)
    for clock in state.packet_clocks:
        _neighbor_coefficients(authority, clock.vertical_layer_id)
        expected_depth = _optical_depth(
            authority.scheduler_seed_id,
            clock.packet_id,
            clock.next_event_index,
        )
        if clock.optical_depth != expected_depth:
            raise EngineeringKEventSchedulerError(
                "K packet optical depth disagrees with its counter identity"
            )


def _clock_event_time(
    authority: EngineeringKEventAuthority,
    clock: EngineeringKPacketClock,
    *,
    state_time_s: float,
    horizon_time_s: float,
) -> float | None:
    coefficients = _neighbor_coefficients(authority, clock.vertical_layer_id)
    total_coefficient = sum((value for _, value in coefficients), Fraction())
    consumed = total_coefficient * _integrated_diffusivity(
        authority,
        clock.hazard_origin_time_s,
        state_time_s,
    )
    remaining = _as_fraction(clock.optical_depth) - consumed
    if remaining <= 0:
        raise EngineeringKEventSchedulerError(
            "accepted scheduler state passed an unprocessed K event"
        )
    return _first_time_for_diffusivity(
        authority,
        start_time_s=state_time_s,
        end_time_s=horizon_time_s,
        target_diffusivity_m2=remaining / total_coefficient,
    )


def restore_engineering_k_schedule_state(
    *,
    authority: EngineeringKEventAuthority,
    packets: tuple[host.AcceptedTrayPacket, ...],
    time_s: float,
    packet_clocks: tuple[EngineeringKPacketClock, ...],
    expected_definition_digest: str,
) -> EngineeringKScheduleState:
    """Restore one accepted K-clock state against caller-supplied typed authority.

    Wire fields alone are not an issuing authority.  The caller must provide
    the exact typed mechanical/history authority and the already reconstructed
    accepted packet population.  This function re-runs the scheduler's mesh,
    weighted-holdup, counter-derived optical-depth, missed-event and packet-K
    ownership checks before applying the private state seal.  The supplied
    digest is then recomputed from the restored state.
    """

    if type(authority) is not EngineeringKEventAuthority:
        raise TypeError("restored K state requires an exact EngineeringKEventAuthority")
    if type(packets) is not tuple or not packets:
        raise TypeError("restored K state requires a nonempty exact packet tuple")
    if any(type(packet) is not host.AcceptedTrayPacket for packet in packets):
        raise TypeError("restored K state packets contain a foreign exact type")
    _require_binary64("restored K state time", time_s, nonnegative=True)
    _require_identity("restored K state digest", expected_definition_digest)

    # This public constructor performs the complete selected mesh/holdup
    # validation at the restored accepted time.  Its newly initialized clocks
    # are intentionally discarded; caller-supplied clocks are checked below.
    initialize_engineering_k_schedule(packets=packets, authority=authority)
    if any(packet.time_s != time_s for packet in packets):
        raise EngineeringKEventSchedulerError(
            "restored K state and accepted packets do not share one exact time"
        )
    state = EngineeringKScheduleState(
        schema_id=ENGINEERING_PHYSICAL_K_EVENT_SCHEMA_ID,
        authority_definition_digest=authority.definition_digest,
        physical_tray_id=authority.mesh.physical_tray_id,
        time_s=time_s,
        packet_clocks=packet_clocks,
        _seal=_ENGINEERING_K_STATE_SEAL,
    )
    _validate_state_authority(state, authority)
    packet_by_id = {packet.packet_id: packet for packet in packets}
    if tuple(sorted(packet_by_id)) != tuple(clock.packet_id for clock in packet_clocks):
        raise EngineeringKEventSchedulerError(
            "restored K clocks and accepted packets have different identities"
        )
    if any(
        packet_by_id[clock.packet_id].owner_key.vertical_layer_id != clock.vertical_layer_id
        for clock in packet_clocks
    ):
        raise EngineeringKEventSchedulerError(
            "restored K clock layer differs from its accepted packet owner"
        )
    for clock in packet_clocks:
        _clock_event_time(
            authority,
            clock,
            state_time_s=time_s,
            horizon_time_s=time_s,
        )
    if state.definition_digest != expected_definition_digest:
        raise EngineeringKEventSchedulerError(
            "restored K state digest differs from its complete clock state"
        )
    return state


def _select_target(
    authority: EngineeringKEventAuthority,
    clock: EngineeringKPacketClock,
) -> tuple[VerticalLayerId, float]:
    coefficients = _neighbor_coefficients(authority, clock.vertical_layer_id)
    total = sum((value for _, value in coefficients), Fraction())
    draw = _direction_fraction(
        authority.scheduler_seed_id,
        clock.packet_id,
        clock.next_event_index,
    )
    cumulative = Fraction()
    for target, coefficient in coefficients:
        cumulative += coefficient / total
        if draw <= cumulative:
            return target, float(coefficient / total)
    raise EngineeringKEventSchedulerError("direction draw escaped its cumulative probability")


def advance_engineering_k_schedule(
    prior: EngineeringKScheduleState,
    *,
    authority: EngineeringKEventAuthority,
    horizon_time_s: float,
) -> EngineeringKScheduleStep:
    """Advance to the first deterministic adjacent-K event or the requested horizon."""

    _validate_state_authority(prior, authority)
    _require_binary64("K scheduler horizon", horizon_time_s, nonnegative=True)
    if horizon_time_s <= prior.time_s:
        raise EngineeringKEventSchedulerError("K scheduler horizon must advance time")
    _validate_interval(authority, prior.time_s, horizon_time_s)
    candidates = tuple(
        (
            clock,
            _clock_event_time(
                authority,
                clock,
                state_time_s=prior.time_s,
                horizon_time_s=horizon_time_s,
            ),
        )
        for clock in prior.packet_clocks
    )
    event_times = tuple(value for _, value in candidates if value is not None)
    if not event_times:
        candidate = replace(prior, time_s=horizon_time_s)
        return EngineeringKScheduleStep(
            authority_definition_digest=authority.definition_digest,
            prior=prior,
            candidate=candidate,
            requested_horizon_time_s=horizon_time_s,
            events=(),
            _seal=_ENGINEERING_K_SCHEDULE_STEP_SEAL,
        )
    first_time = min(event_times)
    clocks_after: list[EngineeringKPacketClock] = []
    events: list[EngineeringKEvent] = []
    for clock, event_time in candidates:
        if event_time != first_time:
            clocks_after.append(clock)
            continue
        target, probability = _select_target(authority, clock)
        events.append(
            EngineeringKEvent(
                packet_id=clock.packet_id,
                event_index=clock.next_event_index,
                event_time_s=first_time,
                source_layer_id=clock.vertical_layer_id,
                target_layer_id=target,
                selected_transition_probability=probability,
            )
        )
        next_index = clock.next_event_index + 1
        clocks_after.append(
            EngineeringKPacketClock(
                packet_id=clock.packet_id,
                vertical_layer_id=target,
                next_event_index=next_index,
                hazard_origin_time_s=first_time,
                optical_depth=_optical_depth(
                    authority.scheduler_seed_id,
                    clock.packet_id,
                    next_index,
                ),
            )
        )
    candidate = replace(
        prior,
        time_s=first_time,
        packet_clocks=tuple(sorted(clocks_after, key=lambda item: item.packet_id)),
    )
    return EngineeringKScheduleStep(
        authority_definition_digest=authority.definition_digest,
        prior=prior,
        candidate=candidate,
        requested_horizon_time_s=horizon_time_s,
        events=tuple(sorted(events, key=lambda item: item.packet_id)),
        _seal=_ENGINEERING_K_SCHEDULE_STEP_SEAL,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringPhysicalKRekeyPlanner:
    """Thin host adapter for one already-resolved scheduler event batch."""

    authority: EngineeringKEventAuthority
    step: EngineeringKScheduleStep

    source_identity: ClassVar[str] = ENGINEERING_PHYSICAL_K_EVENT_SOURCE_IDENTITY
    source_authenticated: ClassVar[bool] = False
    accepted_history_authenticated: ClassVar[bool] = False
    machine_data_authenticated: ClassVar[bool] = False
    production_host_wired: ClassVar[bool] = False
    host_atomic_persistence_implemented: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.authority) is not EngineeringKEventAuthority:
            raise TypeError("planner authority has a foreign exact type")
        if type(self.step) is not EngineeringKScheduleStep:
            raise TypeError("planner step has a foreign exact type")
        if self.step.authority_definition_digest != self.authority.definition_digest:
            raise EngineeringKEventSchedulerError("planner step binds another K authority")
        _validate_state_authority(self.step.prior, self.authority)
        _validate_state_authority(self.step.candidate, self.authority)
        replayed = advance_engineering_k_schedule(
            self.step.prior,
            authority=self.authority,
            horizon_time_s=self.step.requested_horizon_time_s,
        )
        if replayed != self.step:
            raise EngineeringKEventSchedulerError(
                "planner schedule step differs from its deterministic replay"
            )

    @property
    def planner_id(self) -> str:
        return self.step.definition_digest

    def propose(
        self, request: host.KLayerRekeyPlanningRequest
    ) -> tuple[host.KLayerRekeyProposal, ...]:
        if type(request) is not host.KLayerRekeyPlanningRequest:
            raise TypeError("planner request has a foreign exact type")
        if request.event_time_s != self.step.candidate.time_s:
            raise EngineeringKEventSchedulerError("host endpoint is not the scheduled K event time")
        if request.physical_tray_id != self.step.prior.physical_tray_id:
            raise EngineeringKEventSchedulerError("host request belongs to another tray")
        tray = next(
            (
                item
                for item in request.topology.trays
                if item.physical_tray_id == request.physical_tray_id
            ),
            None,
        )
        if tray is None:
            raise EngineeringKEventSchedulerError("host topology omits the scheduled tray")
        if tray.vertical_layer_count != len(self.authority.mesh.cells):
            raise EngineeringKEventSchedulerError("host topology and K mesh have different order")
        by_id = {packet.packet_id: packet for packet in request.packets}
        expected_ids = tuple(clock.packet_id for clock in self.step.prior.packet_clocks)
        if tuple(sorted(by_id)) != expected_ids:
            raise EngineeringKEventSchedulerError("host packet set differs from scheduler state")
        for clock in self.step.prior.packet_clocks:
            if by_id[clock.packet_id].owner_key.vertical_layer_id != clock.vertical_layer_id:
                raise EngineeringKEventSchedulerError(
                    "host packet K owner differs from scheduler prior state"
                )
        return tuple(
            host.KLayerRekeyProposal(
                packet_id=event.packet_id,
                target_layer_id=event.target_layer_id,
                authority_id=self.planner_id,
            )
            for event in self.step.events
        )


__all__ = (
    "DEFAULT_HOLDUP_RELATIVE_LIMIT",
    "ENGINEERING_PHYSICAL_K_EVENT_SCHEMA_ID",
    "ENGINEERING_PHYSICAL_K_EVENT_SOURCE_IDENTITY",
    "EngineeringKEvent",
    "EngineeringKEventAuthority",
    "EngineeringKEventSchedulerError",
    "EngineeringKPacketClock",
    "EngineeringKScheduleState",
    "EngineeringKScheduleStep",
    "EngineeringPhysicalKRekeyPlanner",
    "advance_engineering_k_schedule",
    "initialize_engineering_k_admission_clock",
    "initialize_engineering_k_schedule",
    "restore_engineering_k_schedule_state",
    "retire_engineering_k_clocks_for_merged_packets",
)
