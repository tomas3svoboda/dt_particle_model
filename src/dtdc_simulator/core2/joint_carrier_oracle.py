"""Manufactured K--RTD carrier oracle for the T3C tower host.

This module fixes *mathematical plumbing*, not an unidentified process law.
It builds the separable reference operator

    G = G_RTD + G_arm

on the Cartesian product of vertical tray layer and RTD stage.  The frozen
tanks-in-series rate is derived from dry-matter holdup, locally accepted
dry-matter discharge, and RTD order.  Every arm-redistribution rate is supplied
explicitly by an external authority: this module never infers one from RPM,
geometry, or a numerical mesh.

The matrix convention is ``d m / dt = G m``.  Off-diagonal entries are
non-negative, and the incidence ledger is algebraically conservative when its
declared exits are included.  Binary64 audits expose (rather than hide) the
roundoff residual.  The same transfer objects are applied to every scalar
extensive inventory, including datum-bound common-datum energy.  Opaque
particle/cohort state transport, a timestepper, limiter, clamp, particle-state
averaging, and a production RPM closure are deliberately absent.

This is a manufactured/reference oracle only.  It is deliberately not wired
to :mod:`layered_tray`, :mod:`dtdc_stack`, or the certified particle solver,
and ``physically_qualifying`` is always false.
"""

from __future__ import annotations

import enum
import hashlib
import math
from dataclasses import dataclass
from typing import ClassVar, Iterable, Sequence

from .tray_particle_stateful import (
    PhysicalTrayTopology,
    RTDStageId,
    VerticalLayerId,
)


def _require_finite(name: str, *values: float) -> None:
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"{name} must contain only finite values")


def _require_binary64(name: str, *values: float) -> None:
    if any(type(value) is not float for value in values):
        raise TypeError(f"{name} must contain binary64 float values")
    _require_finite(name, *values)


def _require_nonblank_strings(name: str, *values: str) -> None:
    if any(type(value) is not str for value in values):
        raise TypeError(f"{name} must contain strings")
    if any(not value.strip() for value in values):
        raise ValueError(f"{name} must not contain blank values")


def _canonical_digest(*, domain: str, parts: Iterable[str | int | float | bool]) -> str:
    """Hash a typed, length-prefixed manufactured-oracle identity."""

    digest = hashlib.sha256()
    digest.update(b"DTDC-JOINT-CARRIER-DIGEST-V1\0")
    domain_bytes = domain.encode("utf-8")
    digest.update(len(domain_bytes).to_bytes(8, "big"))
    digest.update(domain_bytes)
    material = tuple(parts)
    digest.update(len(material).to_bytes(8, "big"))
    for part in material:
        if isinstance(part, bool):
            tag, payload = b"b", b"1" if part else b"0"
        elif isinstance(part, int):
            tag, payload = b"i", str(part).encode("ascii")
        elif isinstance(part, float):
            if not math.isfinite(part):
                raise ValueError("canonical digest floats must be finite")
            canonical_float = 0.0 if part == 0.0 else part
            tag, payload = b"f", canonical_float.hex().encode("ascii")
        elif isinstance(part, str):
            tag, payload = b"s", part.encode("utf-8")
        else:  # pragma: no cover - closed union guarded at call sites
            raise TypeError(f"unsupported canonical digest part: {type(part).__name__}")
        digest.update(tag)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, order=True, kw_only=True)
class JointCarrierNode:
    """One node in the independent ``vertical layer x RTD stage`` grid."""

    physical_tray_id: str
    vertical_layer_id: VerticalLayerId
    rtd_stage_id: RTDStageId

    def __post_init__(self) -> None:
        _require_nonblank_strings("carrier node physical_tray_id", self.physical_tray_id)
        if not isinstance(self.vertical_layer_id, VerticalLayerId):
            raise TypeError("vertical_layer_id must be a VerticalLayerId")
        if not isinstance(self.rtd_stage_id, RTDStageId):
            raise TypeError("rtd_stage_id must be an RTDStageId")
        if type(self.vertical_layer_id.value) is not int or isinstance(
            self.vertical_layer_id.value, bool
        ):
            raise TypeError("vertical_layer_id.value must be an integer")
        if type(self.rtd_stage_id.value) is not int or isinstance(self.rtd_stage_id.value, bool):
            raise TypeError("rtd_stage_id.value must be an integer")
        if self.vertical_layer_id.value <= 0 or self.rtd_stage_id.value <= 0:
            raise ValueError("carrier node coordinates must be positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class ArmLayerTransition:
    """One supplied directed arm-redistribution rate between K layers.

    A missing edge means exactly zero rate.  A present edge must therefore be
    strictly positive; silent numerical floors and signed reverse-flux tricks
    are not admissible.
    """

    source_vertical_layer_id: VerticalLayerId
    destination_vertical_layer_id: VerticalLayerId
    rate_s_inv: float

    def __post_init__(self) -> None:
        _require_binary64("arm redistribution rate", self.rate_s_inv)
        if not isinstance(self.source_vertical_layer_id, VerticalLayerId) or not isinstance(
            self.destination_vertical_layer_id, VerticalLayerId
        ):
            raise TypeError("arm transition endpoints must be VerticalLayerId values")
        if (
            type(self.source_vertical_layer_id.value) is not int
            or type(self.destination_vertical_layer_id.value) is not int
        ):
            raise TypeError("arm transition endpoint values must be integers")
        if isinstance(self.source_vertical_layer_id.value, bool) or isinstance(
            self.destination_vertical_layer_id.value, bool
        ):
            raise TypeError("arm transition endpoint values must be integers")
        if (
            self.source_vertical_layer_id.value <= 0
            or self.destination_vertical_layer_id.value <= 0
        ):
            raise ValueError("arm transition endpoint values must be positive")
        if self.source_vertical_layer_id == self.destination_vertical_layer_id:
            raise ValueError("arm redistribution needs distinct K layers")
        if self.rate_s_inv <= 0.0:
            raise ValueError("a declared arm transition rate must be strictly positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class ArmRedistributionAuthority:
    """Explicit input authority for an otherwise-unselected arm closure.

    The zero/positive-RPM structural rules are owner-approved.  The numerical
    relationship between RPM and the supplied rates remains outside this
    class and must identify its own ``closure_id`` and ``source_identity``.
    """

    physical_tray_id: str
    shaft_history_id: str
    shaft_speed_rpm: float
    tray_geometry_id: str
    closure_id: str
    source_identity: str
    active_bed: bool
    transitions: tuple[ArmLayerTransition, ...]

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_binary64("shaft speed", self.shaft_speed_rpm)
        if type(self.active_bed) is not bool:
            raise TypeError("active_bed must be a bool")
        if not isinstance(self.transitions, tuple):
            raise TypeError("arm authority transitions must be an immutable tuple")
        if any(not isinstance(item, ArmLayerTransition) for item in self.transitions):
            raise TypeError("arm authority transitions must be ArmLayerTransition values")
        _require_nonblank_strings(
            "arm authority identifiers",
            self.physical_tray_id,
            self.shaft_history_id,
            self.tray_geometry_id,
            self.closure_id,
            self.source_identity,
        )
        if self.shaft_speed_rpm < 0.0:
            raise ValueError("shaft_speed_rpm must be non-negative")
        directed_edges = tuple(
            (
                item.source_vertical_layer_id,
                item.destination_vertical_layer_id,
            )
            for item in self.transitions
        )
        if len(set(directed_edges)) != len(directed_edges):
            raise ValueError("arm authority contains a duplicate directed K edge")
        canonical_transitions = tuple(
            sorted(
                self.transitions,
                key=lambda item: (
                    item.source_vertical_layer_id.value,
                    item.destination_vertical_layer_id.value,
                ),
            )
        )
        if self.transitions != canonical_transitions:
            raise ValueError("arm authority transitions must use canonical K-edge order")
        if (self.shaft_speed_rpm == 0.0 or not self.active_bed) and self.transitions:
            raise ValueError("zero-RPM or inactive-bed arm redistribution must be exactly zero")
        if self.shaft_speed_rpm > 0.0 and self.active_bed and not self.transitions:
            raise ValueError("a positive-RPM active bed requires an explicit nonzero arm operator")


def validate_common_shaft_snapshot(
    authorities: Sequence[ArmRedistributionAuthority],
) -> None:
    """Require one physical shaft history and speed across all DT trays."""

    if not authorities:
        raise ValueError("at least one arm authority is required")
    if any(not isinstance(item, ArmRedistributionAuthority) for item in authorities):
        raise TypeError("shaft snapshot entries must be ArmRedistributionAuthority values")
    tray_ids = tuple(item.physical_tray_id for item in authorities)
    if len(set(tray_ids)) != len(tray_ids):
        raise ValueError("arm authorities must identify distinct physical trays")
    history_id = authorities[0].shaft_history_id
    speed_rpm = authorities[0].shaft_speed_rpm
    if any(item.shaft_history_id != history_id for item in authorities[1:]):
        raise ValueError("all DT trays must use one shaft-history authority")
    if any(item.shaft_speed_rpm != speed_rpm for item in authorities[1:]):
        raise ValueError("all DT trays must use the same instantaneous shaft speed")


def _validate_topology_runtime(topology: PhysicalTrayTopology) -> None:
    if not isinstance(topology, PhysicalTrayTopology):
        raise TypeError("topology must be a PhysicalTrayTopology")
    _require_nonblank_strings("carrier topology physical_tray_id", topology.physical_tray_id)
    if type(topology.vertical_layer_count) is not int or type(topology.rtd_stage_count) is not int:
        raise TypeError("carrier topology counts must be integers")
    if topology.vertical_layer_count <= 0 or topology.rtd_stage_count <= 0:
        raise ValueError("carrier topology counts must be positive")


def _validate_active_arm_graph(
    *,
    topology: PhysicalTrayTopology,
    authority: ArmRedistributionAuthority,
) -> None:
    """Require the supplied positive-RPM graph to connect every declared K layer.

    Connectivity is evaluated without assigning a direction to the physical
    redistribution law.  Requiring strong connectivity would select more
    physics than the owner has approved; allowing an isolated layer would fail
    the approved nonzero connected-graph acceptance condition.
    """

    if authority.shaft_speed_rpm == 0.0 or not authority.active_bed:
        return
    neighbours = {layer: set() for layer in range(1, topology.vertical_layer_count + 1)}
    for transition in authority.transitions:
        source = transition.source_vertical_layer_id.value
        destination = transition.destination_vertical_layer_id.value
        neighbours[source].add(destination)
        neighbours[destination].add(source)
    reached = {1}
    frontier = [1]
    while frontier:
        source = frontier.pop()
        for destination in neighbours[source] - reached:
            reached.add(destination)
            frontier.append(destination)
    if reached != set(neighbours):
        raise ValueError("a positive-RPM active-bed arm graph must connect every declared K layer")


class CarrierMechanism(enum.Enum):
    """Mechanisms present in the separable manufactured oracle."""

    RTD_ADVANCE = "RTD_ADVANCE"
    ARM_REDISTRIBUTION = "ARM_REDISTRIBUTION"


class DryMatterDischargeBasis(enum.Enum):
    """Admissible flow basis for the frozen PHY-009 residence identity."""

    LOCALLY_ACCEPTED_DISCHARGE = "LOCALLY_ACCEPTED_DRY_MATTER_DISCHARGE"


@dataclass(frozen=True, slots=True, kw_only=True)
class InternalCarrierTransition:
    source: JointCarrierNode
    destination: JointCarrierNode
    rate_s_inv: float
    mechanism: CarrierMechanism
    authority_id: str

    def __post_init__(self) -> None:
        _require_binary64("internal carrier rate", self.rate_s_inv)
        if not isinstance(self.source, JointCarrierNode) or not isinstance(
            self.destination, JointCarrierNode
        ):
            raise TypeError("carrier transition endpoints must be JointCarrierNode values")
        if not isinstance(self.mechanism, CarrierMechanism):
            raise TypeError("carrier transition mechanism must be a CarrierMechanism")
        if self.source == self.destination:
            raise ValueError("an internal carrier transition needs distinct nodes")
        if self.rate_s_inv <= 0.0:
            raise ValueError("an internal carrier transition rate must be positive")
        _require_nonblank_strings("carrier transition authority_id", self.authority_id)


@dataclass(frozen=True, slots=True, kw_only=True)
class CarrierExit:
    source: JointCarrierNode
    rate_s_inv: float
    authority_id: str
    physical_port_id: str

    def __post_init__(self) -> None:
        _require_binary64("carrier exit rate", self.rate_s_inv)
        if not isinstance(self.source, JointCarrierNode):
            raise TypeError("carrier exit source must be a JointCarrierNode")
        if self.rate_s_inv <= 0.0:
            raise ValueError("a carrier exit rate must be positive")
        _require_nonblank_strings(
            "carrier exit authority and physical port IDs",
            self.authority_id,
            self.physical_port_id,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class RTDMoments:
    mean_s: float
    variance_s2: float

    def __post_init__(self) -> None:
        _require_binary64("RTD moments", self.mean_s, self.variance_s2)
        if self.mean_s <= 0.0 or self.variance_s2 <= 0.0:
            raise ValueError("RTD mean and variance must be strictly positive")


@dataclass(frozen=True, slots=True, kw_only=True)
class RTDResidenceAuthority:
    """Declared PHY-009 residence snapshot for one physical tray.

    ``dry_matter_holdup_kg`` is total physical-tray dry-matter holdup, never a K
    layer or RTD-cell inventory.  The caller declares a locally accepted
    dry-meal discharge, not offered feed.  This oracle records and type-checks
    that basis but does not verify the still-open binding to a production
    cell/airlock transaction.
    """

    physical_tray_id: str
    snapshot_id: str
    authority_id: str
    source_identity: str
    holdup_state_id: str
    accepted_discharge_flow_id: str
    discharge_basis: DryMatterDischargeBasis
    dry_matter_holdup_kg: float
    accepted_dry_matter_discharge_kg_s: float

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_flow_binding_selected: ClassVar[bool] = False

    def __post_init__(self) -> None:
        identifiers = (
            self.physical_tray_id,
            self.snapshot_id,
            self.authority_id,
            self.source_identity,
            self.holdup_state_id,
            self.accepted_discharge_flow_id,
        )
        _require_nonblank_strings("RTD residence authority identifiers", *identifiers)
        _require_binary64(
            "RTD residence authority",
            self.dry_matter_holdup_kg,
            self.accepted_dry_matter_discharge_kg_s,
        )
        if not isinstance(self.discharge_basis, DryMatterDischargeBasis):
            raise TypeError("RTD residence discharge_basis has the wrong type")
        if self.discharge_basis is not DryMatterDischargeBasis.LOCALLY_ACCEPTED_DISCHARGE:
            raise ValueError("RTD residence requires locally accepted dry-matter discharge")
        if self.dry_matter_holdup_kg <= 0.0:
            raise ValueError("RTD dry-matter holdup must be strictly positive")
        if self.accepted_dry_matter_discharge_kg_s <= 0.0:
            raise ValueError("accepted dry-matter discharge must be strictly positive")
        try:
            mean_residence_time_s = (
                self.dry_matter_holdup_kg / self.accepted_dry_matter_discharge_kg_s
            )
        except OverflowError as error:
            raise ValueError("derived RTD residence coordinates must be representable") from error
        _require_finite(
            "derived RTD residence coordinates",
            mean_residence_time_s,
        )
        if mean_residence_time_s <= 0.0:
            raise ValueError("derived RTD residence coordinates must be strictly positive")

    @property
    def mean_residence_time_s(self) -> float:
        return self.dry_matter_holdup_kg / self.accepted_dry_matter_discharge_kg_s

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="rtd-residence-authority-v1",
            parts=(
                self.physical_tray_id,
                self.snapshot_id,
                self.authority_id,
                self.source_identity,
                self.holdup_state_id,
                self.accepted_discharge_flow_id,
                self.discharge_basis.value,
                self.dry_matter_holdup_kg,
                self.accepted_dry_matter_discharge_kg_s,
            ),
        )


def _derived_stage_rate_s_inv(*, stage_count: int, mean_residence_time_s: float) -> float:
    try:
        rate_s_inv = stage_count / mean_residence_time_s
    except (OverflowError, ZeroDivisionError) as error:
        raise ValueError("derived RTD stage rate must be representable") from error
    _require_finite("derived RTD stage rate", rate_s_inv)
    if rate_s_inv <= 0.0:
        raise ValueError("derived RTD stage rate must be strictly positive")
    return rate_s_inv


def tanks_in_series_moments(*, stage_count: int, mean_residence_time_s: float) -> RTDMoments:
    """Return the analytic Erlang moments for equal tanks in series."""

    _require_binary64("mean residence time", mean_residence_time_s)
    if not isinstance(stage_count, int) or isinstance(stage_count, bool):
        raise TypeError("stage_count must be an integer")
    if stage_count <= 0:
        raise ValueError("stage_count must be positive")
    if mean_residence_time_s <= 0.0:
        raise ValueError("mean_residence_time_s must be strictly positive")
    try:
        variance_s2 = mean_residence_time_s**2 / stage_count
    except OverflowError as error:
        raise ValueError("RTD moment coordinates must be representable") from error
    return RTDMoments(mean_s=mean_residence_time_s, variance_s2=variance_s2)


@dataclass(frozen=True, slots=True, kw_only=True)
class CarrierInventory:
    """Frozen PHY-009 extensive carrier basis at one oracle node."""

    dry_matter_kg: float
    residual_oil_label_kg: float
    attached_hexane_kg: float
    internal_hexane_kg: float
    external_water_kg: float
    retained_water_kg: float
    common_datum_energy_j: float
    energy_datum_id: str

    def __post_init__(self) -> None:
        values = self.as_tuple()
        _require_binary64("carrier inventory", *values)
        _require_nonblank_strings("carrier inventory energy_datum_id", self.energy_datum_id)
        masses = values[:-1]
        if any(value < 0.0 for value in masses):
            raise ValueError("carrier mass inventories must be non-negative")
        if self.residual_oil_label_kg > self.dry_matter_kg:
            raise ValueError("residual-oil label cannot exceed dry matter")

    def as_tuple(self) -> tuple[float, ...]:
        return (
            self.dry_matter_kg,
            self.residual_oil_label_kg,
            self.attached_hexane_kg,
            self.internal_hexane_kg,
            self.external_water_kg,
            self.retained_water_kg,
            self.common_datum_energy_j,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class CarrierInventoryRate:
    dry_matter_kg_s: float
    residual_oil_label_kg_s: float
    attached_hexane_kg_s: float
    internal_hexane_kg_s: float
    external_water_kg_s: float
    retained_water_kg_s: float
    common_datum_energy_w: float

    def __post_init__(self) -> None:
        _require_binary64("carrier inventory rate", *self.as_tuple())

    def as_tuple(self) -> tuple[float, ...]:
        return (
            self.dry_matter_kg_s,
            self.residual_oil_label_kg_s,
            self.attached_hexane_kg_s,
            self.internal_hexane_kg_s,
            self.external_water_kg_s,
            self.retained_water_kg_s,
            self.common_datum_energy_w,
        )

    @classmethod
    def from_values(cls, values: Iterable[float]) -> CarrierInventoryRate:
        items = tuple(values)
        if len(items) != 7:
            raise ValueError("a carrier rate requires exactly seven extensives")
        return cls(
            dry_matter_kg_s=items[0],
            residual_oil_label_kg_s=items[1],
            attached_hexane_kg_s=items[2],
            internal_hexane_kg_s=items[3],
            external_water_kg_s=items[4],
            retained_water_kg_s=items[5],
            common_datum_energy_w=items[6],
        )


def _input_inventory_digest(
    *,
    nodes: Sequence[JointCarrierNode],
    inventories: Sequence[CarrierInventory],
) -> str:
    parts: list[str | int | float | bool] = []
    for node, inventory in zip(nodes, inventories, strict=True):
        parts.extend(
            (
                node.physical_tray_id,
                node.vertical_layer_id.value,
                node.rtd_stage_id.value,
                inventory.energy_datum_id,
                *inventory.as_tuple(),
            )
        )
    return _canonical_digest(domain="carrier-input-inventory-v1", parts=parts)


@dataclass(frozen=True, slots=True, kw_only=True)
class InternalTransferRate:
    source: JointCarrierNode
    destination: JointCarrierNode
    mechanism: CarrierMechanism
    rate_s_inv: float
    authority_id: str
    transferred: CarrierInventoryRate

    def __post_init__(self) -> None:
        _require_binary64("internal transfer rate", self.rate_s_inv)
        if not isinstance(self.source, JointCarrierNode) or not isinstance(
            self.destination, JointCarrierNode
        ):
            raise TypeError("internal transfer endpoints must be JointCarrierNode values")
        if not isinstance(self.mechanism, CarrierMechanism):
            raise TypeError("internal transfer mechanism must be a CarrierMechanism")
        if self.source == self.destination:
            raise ValueError("an internal transfer needs distinct nodes")
        if self.source.physical_tray_id != self.destination.physical_tray_id:
            raise ValueError("an internal transfer cannot cross physical trays")
        if self.mechanism is CarrierMechanism.RTD_ADVANCE:
            if (
                self.source.vertical_layer_id != self.destination.vertical_layer_id
                or self.destination.rtd_stage_id.value != self.source.rtd_stage_id.value + 1
            ):
                raise ValueError("an RTD transfer must preserve K and advance one stage")
        elif self.mechanism is CarrierMechanism.ARM_REDISTRIBUTION:
            if (
                self.source.rtd_stage_id != self.destination.rtd_stage_id
                or self.source.vertical_layer_id == self.destination.vertical_layer_id
            ):
                raise ValueError("an arm transfer must preserve RTD stage and change K")
        else:  # pragma: no cover - guarded by the enum type check
            raise TypeError("internal transfer mechanism is unsupported")
        if self.rate_s_inv <= 0.0:
            raise ValueError("internal transfer rate must be positive")
        _require_nonblank_strings("internal transfer authority_id", self.authority_id)
        if not isinstance(self.transferred, CarrierInventoryRate):
            raise TypeError("internal transfer must carry a CarrierInventoryRate")


@dataclass(frozen=True, slots=True, kw_only=True)
class ExitTransferRate:
    source: JointCarrierNode
    rate_s_inv: float
    authority_id: str
    physical_port_id: str
    transferred: CarrierInventoryRate

    def __post_init__(self) -> None:
        _require_binary64("exit transfer rate", self.rate_s_inv)
        if not isinstance(self.source, JointCarrierNode):
            raise TypeError("exit transfer source must be a JointCarrierNode")
        if self.rate_s_inv <= 0.0:
            raise ValueError("exit transfer rate must be positive")
        _require_nonblank_strings(
            "exit transfer authority and physical port IDs",
            self.authority_id,
            self.physical_port_id,
        )
        if not isinstance(self.transferred, CarrierInventoryRate):
            raise TypeError("exit transfer must carry a CarrierInventoryRate")


@dataclass(frozen=True, slots=True, kw_only=True)
class CarrierRateEvaluation:
    energy_datum_id: str
    operator_definition_digest: str
    residence_authority_digest: str
    input_inventory_digest: str
    nodes: tuple[JointCarrierNode, ...]
    node_rates: tuple[CarrierInventoryRate, ...]
    internal_transfers: tuple[InternalTransferRate, ...]
    exit_transfers: tuple[ExitTransferRate, ...]

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank_strings(
            "carrier evaluation identities",
            self.energy_datum_id,
            self.operator_definition_digest,
            self.residence_authority_digest,
            self.input_inventory_digest,
        )
        if not all(
            isinstance(items, tuple)
            for items in (
                self.nodes,
                self.node_rates,
                self.internal_transfers,
                self.exit_transfers,
            )
        ):
            raise TypeError("carrier evaluation collections must be immutable tuples")
        if len(self.nodes) != len(self.node_rates):
            raise ValueError("carrier evaluation needs one rate for every canonical node")
        if len(set(self.nodes)) != len(self.nodes):
            raise ValueError("carrier evaluation contains duplicate nodes")
        if any(not isinstance(node, JointCarrierNode) for node in self.nodes):
            raise TypeError("carrier evaluation nodes must be JointCarrierNode values")
        if any(not isinstance(rate, CarrierInventoryRate) for rate in self.node_rates):
            raise TypeError("carrier evaluation node rates must be CarrierInventoryRate values")
        if any(
            not isinstance(transfer, InternalTransferRate) for transfer in self.internal_transfers
        ):
            raise TypeError("carrier evaluation internal transfers have the wrong type")
        if any(not isinstance(transfer, ExitTransferRate) for transfer in self.exit_transfers):
            raise TypeError("carrier evaluation exit transfers have the wrong type")
        node_set = set(self.nodes)
        if len({node.physical_tray_id for node in self.nodes}) != 1:
            raise ValueError("carrier evaluation nodes must belong to one physical tray")
        canonical_nodes = tuple(
            sorted(
                self.nodes,
                key=lambda node: (
                    node.physical_tray_id,
                    node.rtd_stage_id.value,
                    node.vertical_layer_id.value,
                ),
            )
        )
        if self.nodes != canonical_nodes:
            raise ValueError("carrier evaluation nodes must use canonical coordinate order")
        if any(
            transfer.source not in node_set or transfer.destination not in node_set
            for transfer in self.internal_transfers
        ):
            raise ValueError("carrier evaluation internal transfer endpoint is outside its nodes")
        if any(transfer.source not in node_set for transfer in self.exit_transfers):
            raise ValueError("carrier evaluation exit source is outside its nodes")
        canonical_internal = tuple(
            sorted(
                self.internal_transfers,
                key=lambda item: (
                    (
                        0,
                        item.source.vertical_layer_id.value,
                        item.source.rtd_stage_id.value,
                        item.destination.vertical_layer_id.value,
                        item.destination.rtd_stage_id.value,
                    )
                    if item.mechanism is CarrierMechanism.RTD_ADVANCE
                    else (
                        1,
                        item.source.rtd_stage_id.value,
                        item.source.vertical_layer_id.value,
                        item.destination.vertical_layer_id.value,
                        item.destination.rtd_stage_id.value,
                    )
                ),
            )
        )
        if self.internal_transfers != canonical_internal:
            raise ValueError("carrier evaluation internal transfers must use canonical order")
        canonical_exits = tuple(
            sorted(
                self.exit_transfers,
                key=lambda item: (
                    item.source.vertical_layer_id.value,
                    item.source.rtd_stage_id.value,
                    item.physical_port_id,
                ),
            )
        )
        if self.exit_transfers != canonical_exits:
            raise ValueError("carrier evaluation exits must use canonical order")

    @property
    def definition_digest(self) -> str:
        parts: list[str | int | float | bool] = [
            self.energy_datum_id,
            self.operator_definition_digest,
            self.residence_authority_digest,
            self.input_inventory_digest,
        ]
        for node, rate in zip(self.nodes, self.node_rates, strict=True):
            parts.extend(
                (
                    node.physical_tray_id,
                    node.vertical_layer_id.value,
                    node.rtd_stage_id.value,
                    *rate.as_tuple(),
                )
            )
        for transfer in self.internal_transfers:
            parts.extend(
                (
                    transfer.source.physical_tray_id,
                    transfer.source.vertical_layer_id.value,
                    transfer.source.rtd_stage_id.value,
                    transfer.destination.physical_tray_id,
                    transfer.destination.vertical_layer_id.value,
                    transfer.destination.rtd_stage_id.value,
                    transfer.mechanism.value,
                    transfer.rate_s_inv,
                    transfer.authority_id,
                    *transfer.transferred.as_tuple(),
                )
            )
        for transfer in self.exit_transfers:
            parts.extend(
                (
                    transfer.source.physical_tray_id,
                    transfer.source.vertical_layer_id.value,
                    transfer.source.rtd_stage_id.value,
                    transfer.rate_s_inv,
                    transfer.authority_id,
                    transfer.physical_port_id,
                    *transfer.transferred.as_tuple(),
                )
            )
        return _canonical_digest(domain="carrier-rate-evaluation-v1", parts=parts)

    @property
    def exit_rate(self) -> CarrierInventoryRate:
        return CarrierInventoryRate.from_values(
            math.fsum(transfer.transferred.as_tuple()[index] for transfer in self.exit_transfers)
            for index in range(7)
        )

    @property
    def conservation_residual(self) -> CarrierInventoryRate:
        exit_values = self.exit_rate.as_tuple()
        return CarrierInventoryRate.from_values(
            math.fsum(
                (
                    *(rate.as_tuple()[index] for rate in self.node_rates),
                    exit_values[index],
                )
            )
            for index in range(7)
        )

    @property
    def conservation_relative_residuals(self) -> tuple[float, ...]:
        """Return componentwise roundoff residuals scaled by audited rates.

        No acceptance tolerance is embedded here.  Qualification compares the
        exposed values with a declared operation-count/ulp bound and the
        applicable conservation requirement.
        """

        residuals = self.conservation_residual.as_tuple()
        exit_values = self.exit_rate.as_tuple()
        scaled: list[float] = []
        for index, residual in enumerate(residuals):
            scale = math.fsum(
                (
                    *(abs(rate.as_tuple()[index]) for rate in self.node_rates),
                    abs(exit_values[index]),
                )
            )
            if scale == 0.0:
                scaled.append(0.0 if residual == 0.0 else math.inf)
            else:
                scaled.append(abs(residual) / scale)
        return tuple(scaled)


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedJointCarrierOperator:
    """Separable K--RTD reference operator, not a production closure."""

    topology: PhysicalTrayTopology
    nodes: tuple[JointCarrierNode, ...]
    transitions: tuple[InternalCarrierTransition, ...]
    exits: tuple[CarrierExit, ...]
    residence_authority: RTDResidenceAuthority
    energy_datum_id: str
    exit_port_id: str
    arm_authority: ArmRedistributionAuthority

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_closure_selected: ClassVar[bool] = False
    production_flow_binding_selected: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False
    shaft_work_energy_model_selected: ClassVar[bool] = False
    time_integrator_selected: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _validate_topology_runtime(self.topology)
        if not isinstance(self.arm_authority, ArmRedistributionAuthority):
            raise TypeError("carrier arm authority has the wrong type")
        if not isinstance(self.residence_authority, RTDResidenceAuthority):
            raise TypeError("carrier residence authority has the wrong type")
        if not all(
            isinstance(items, tuple) for items in (self.nodes, self.transitions, self.exits)
        ):
            raise TypeError("carrier nodes, transitions, and exits must be immutable tuples")
        if any(not isinstance(item, JointCarrierNode) for item in self.nodes):
            raise TypeError("carrier nodes must be JointCarrierNode values")
        if any(not isinstance(item, InternalCarrierTransition) for item in self.transitions):
            raise TypeError("carrier transitions must be InternalCarrierTransition values")
        if any(not isinstance(item, CarrierExit) for item in self.exits):
            raise TypeError("carrier exits must be CarrierExit values")
        _require_nonblank_strings(
            "carrier operator energy datum and exit port IDs",
            self.energy_datum_id,
            self.exit_port_id,
        )
        if self.residence_authority.physical_tray_id != self.topology.physical_tray_id:
            raise ValueError("RTD residence authority belongs to another physical tray")
        _derived_stage_rate_s_inv(
            stage_count=self.topology.rtd_stage_count,
            mean_residence_time_s=self.residence_authority.mean_residence_time_s,
        )
        expected_nodes = tuple(
            JointCarrierNode(
                physical_tray_id=self.topology.physical_tray_id,
                vertical_layer_id=VerticalLayerId(k),
                rtd_stage_id=RTDStageId(n),
            )
            for n in range(1, self.topology.rtd_stage_count + 1)
            for k in range(1, self.topology.vertical_layer_count + 1)
        )
        if self.nodes != expected_nodes:
            raise ValueError("carrier nodes must be the complete canonical K x RTD grid")
        if self.arm_authority.physical_tray_id != self.topology.physical_tray_id:
            raise ValueError("arm authority belongs to another physical tray")
        node_set = set(self.nodes)
        identities: list[tuple[CarrierMechanism, JointCarrierNode, JointCarrierNode]] = []
        for transition in self.transitions:
            if transition.source not in node_set or transition.destination not in node_set:
                raise ValueError("carrier transition references a node outside the topology")
            identities.append((transition.mechanism, transition.source, transition.destination))
            source_k = transition.source.vertical_layer_id
            destination_k = transition.destination.vertical_layer_id
            source_n = transition.source.rtd_stage_id.value
            destination_n = transition.destination.rtd_stage_id.value
            if transition.mechanism is CarrierMechanism.RTD_ADVANCE:
                if source_k != destination_k or destination_n != source_n + 1:
                    raise ValueError("RTD_ADVANCE must preserve K and advance one RTD stage")
            elif transition.mechanism is CarrierMechanism.ARM_REDISTRIBUTION:
                if source_n != destination_n or source_k == destination_k:
                    raise ValueError("ARM_REDISTRIBUTION must preserve RTD stage and change K")
            else:  # pragma: no cover - guarded by InternalCarrierTransition
                raise TypeError("carrier operator contains an unsupported mechanism")
        if len(set(identities)) != len(identities):
            raise ValueError("carrier operator contains a duplicate transition")
        exit_sources = tuple(item.source for item in self.exits)
        if len(set(exit_sources)) != len(exit_sources):
            raise ValueError("carrier operator contains duplicate exits")
        if any(source not in node_set for source in exit_sources):
            raise ValueError("carrier exit references a node outside the topology")
        for source in self.nodes:
            outgoing_rates = [
                transition.rate_s_inv
                for transition in self.transitions
                if transition.source == source
            ]
            outgoing_rates.extend(
                exit_rate.rate_s_inv for exit_rate in self.exits if exit_rate.source == source
            )
            try:
                total_outgoing_rate = math.fsum(outgoing_rates)
            except OverflowError as error:
                raise ValueError("outgoing carrier rate sum must be representable") from error
            _require_finite("outgoing carrier rate sum", total_outgoing_rate)

        expected_rtd = {
            (
                JointCarrierNode(
                    physical_tray_id=self.topology.physical_tray_id,
                    vertical_layer_id=VerticalLayerId(k),
                    rtd_stage_id=RTDStageId(n),
                ),
                JointCarrierNode(
                    physical_tray_id=self.topology.physical_tray_id,
                    vertical_layer_id=VerticalLayerId(k),
                    rtd_stage_id=RTDStageId(n + 1),
                ),
                self.rtd_rate_s_inv,
                self.residence_authority.authority_id,
            )
            for k in range(1, self.topology.vertical_layer_count + 1)
            for n in range(1, self.topology.rtd_stage_count)
        }
        observed_rtd = {
            (
                item.source,
                item.destination,
                item.rate_s_inv,
                item.authority_id,
            )
            for item in self.transitions
            if item.mechanism is CarrierMechanism.RTD_ADVANCE
        }
        if observed_rtd != expected_rtd:
            raise ValueError("RTD transitions must be the complete K-independent Erlang chain")
        expected_exits = {
            (
                JointCarrierNode(
                    physical_tray_id=self.topology.physical_tray_id,
                    vertical_layer_id=VerticalLayerId(k),
                    rtd_stage_id=RTDStageId(self.topology.rtd_stage_count),
                ),
                self.rtd_rate_s_inv,
                self.residence_authority.authority_id,
                self.exit_port_id,
            )
            for k in range(1, self.topology.vertical_layer_count + 1)
        }
        observed_exits = {
            (
                item.source,
                item.rate_s_inv,
                item.authority_id,
                item.physical_port_id,
            )
            for item in self.exits
        }
        if observed_exits != expected_exits:
            raise ValueError("RTD exits must be complete, K-independent, and final-stage only")
        expected_arm = {
            (
                JointCarrierNode(
                    physical_tray_id=self.topology.physical_tray_id,
                    vertical_layer_id=item.source_vertical_layer_id,
                    rtd_stage_id=RTDStageId(n),
                ),
                JointCarrierNode(
                    physical_tray_id=self.topology.physical_tray_id,
                    vertical_layer_id=item.destination_vertical_layer_id,
                    rtd_stage_id=RTDStageId(n),
                ),
                item.rate_s_inv,
                self.arm_authority.source_identity,
            )
            for item in self.arm_authority.transitions
            for n in range(1, self.topology.rtd_stage_count + 1)
        }
        observed_arm = {
            (
                item.source,
                item.destination,
                item.rate_s_inv,
                item.authority_id,
            )
            for item in self.transitions
            if item.mechanism is CarrierMechanism.ARM_REDISTRIBUTION
        }
        if observed_arm != expected_arm:
            raise ValueError("arm rates must be age-preserving copies of the supplied K operator")
        canonical_transitions = tuple(
            sorted(
                self.transitions,
                key=lambda item: (
                    (
                        0,
                        item.source.vertical_layer_id.value,
                        item.source.rtd_stage_id.value,
                        item.destination.vertical_layer_id.value,
                        item.destination.rtd_stage_id.value,
                    )
                    if item.mechanism is CarrierMechanism.RTD_ADVANCE
                    else (
                        1,
                        item.source.rtd_stage_id.value,
                        item.source.vertical_layer_id.value,
                        item.destination.vertical_layer_id.value,
                        item.destination.rtd_stage_id.value,
                    )
                ),
            )
        )
        if self.transitions != canonical_transitions:
            raise ValueError("carrier transitions must use canonical mechanism/coordinate order")
        canonical_exits = tuple(
            sorted(
                self.exits,
                key=lambda item: (
                    item.source.vertical_layer_id.value,
                    item.source.rtd_stage_id.value,
                    item.physical_port_id,
                ),
            )
        )
        if self.exits != canonical_exits:
            raise ValueError("carrier exits must use canonical coordinate order")

    @property
    def rtd_mean_residence_time_s(self) -> float:
        return self.residence_authority.mean_residence_time_s

    @property
    def rtd_rate_s_inv(self) -> float:
        return _derived_stage_rate_s_inv(
            stage_count=self.topology.rtd_stage_count,
            mean_residence_time_s=self.rtd_mean_residence_time_s,
        )

    @property
    def rtd_moments(self) -> RTDMoments:
        return tanks_in_series_moments(
            stage_count=self.topology.rtd_stage_count,
            mean_residence_time_s=self.rtd_mean_residence_time_s,
        )

    @property
    def definition_digest(self) -> str:
        """Bind the complete manufactured operator and its external authorities."""

        parts: list[str | int | float | bool] = [
            self.topology.physical_tray_id,
            self.topology.vertical_layer_count,
            self.topology.rtd_stage_count,
            self.rtd_mean_residence_time_s,
            self.rtd_rate_s_inv,
            self.residence_authority.definition_digest,
            self.energy_datum_id,
            self.exit_port_id,
            self.arm_authority.physical_tray_id,
            self.arm_authority.shaft_history_id,
            self.arm_authority.shaft_speed_rpm,
            self.arm_authority.tray_geometry_id,
            self.arm_authority.closure_id,
            self.arm_authority.source_identity,
            self.arm_authority.active_bed,
        ]
        for transition in self.transitions:
            parts.extend(
                (
                    transition.source.physical_tray_id,
                    transition.source.vertical_layer_id.value,
                    transition.source.rtd_stage_id.value,
                    transition.destination.physical_tray_id,
                    transition.destination.vertical_layer_id.value,
                    transition.destination.rtd_stage_id.value,
                    transition.rate_s_inv,
                    transition.mechanism.value,
                    transition.authority_id,
                )
            )
        for exit_rate in self.exits:
            parts.extend(
                (
                    exit_rate.source.physical_tray_id,
                    exit_rate.source.vertical_layer_id.value,
                    exit_rate.source.rtd_stage_id.value,
                    exit_rate.rate_s_inv,
                    exit_rate.authority_id,
                    exit_rate.physical_port_id,
                )
            )
        return _canonical_digest(domain="manufactured-joint-carrier-operator-v1", parts=parts)

    def _generator_matrix_for(
        self,
        *,
        mechanisms: frozenset[CarrierMechanism],
        include_rtd_exits: bool,
    ) -> tuple[tuple[float, ...], ...]:
        indices = {node: index for index, node in enumerate(self.nodes)}
        columns: list[list[list[float]]] = [[[] for _ in self.nodes] for _ in self.nodes]
        outgoing: list[list[float]] = [[] for _ in self.nodes]
        for transition in self.transitions:
            if transition.mechanism not in mechanisms:
                continue
            source = indices[transition.source]
            destination = indices[transition.destination]
            columns[destination][source].append(transition.rate_s_inv)
            outgoing[source].append(transition.rate_s_inv)
        if include_rtd_exits:
            for exit_rate in self.exits:
                outgoing[indices[exit_rate.source]].append(exit_rate.rate_s_inv)
        for source, rates in enumerate(outgoing):
            columns[source][source].append(-math.fsum(rates))
        return tuple(tuple(math.fsum(contributions) for contributions in row) for row in columns)

    def generator_matrix(self) -> tuple[tuple[float, ...], ...]:
        """Return ``G`` in the declared column/source convention."""

        return self._generator_matrix_for(
            mechanisms=frozenset(CarrierMechanism),
            include_rtd_exits=True,
        )

    def rtd_joint_generator_matrix(self) -> tuple[tuple[float, ...], ...]:
        """Return the K-replicated RTD subgenerator, including physical exits."""

        return self._generator_matrix_for(
            mechanisms=frozenset((CarrierMechanism.RTD_ADVANCE,)),
            include_rtd_exits=True,
        )

    def arm_generator_matrix(self) -> tuple[tuple[float, ...], ...]:
        """Return the conservative age-preserving arm-redistribution generator."""

        return self._generator_matrix_for(
            mechanisms=frozenset((CarrierMechanism.ARM_REDISTRIBUTION,)),
            include_rtd_exits=False,
        )

    def splitting_commutator_residual(self) -> tuple[tuple[float, ...], ...]:
        """Return ``G_RTD G_arm - G_arm G_RTD`` without a pass tolerance.

        The manufactured separable oracle should commute algebraically.  Raw
        residuals are exposed because non-dyadic supplied rates may leave a few
        floating-point ulps; qualification chooses a scale-aware bound rather
        than hiding one in this architecture module.
        """

        rtd = self.rtd_joint_generator_matrix()
        arm = self.arm_generator_matrix()
        dimension = len(self.nodes)
        return tuple(
            tuple(
                math.fsum(
                    (
                        math.fsum(
                            rtd[row][inner] * arm[inner][column] for inner in range(dimension)
                        ),
                        -math.fsum(
                            arm[row][inner] * rtd[inner][column] for inner in range(dimension)
                        ),
                    )
                )
                for column in range(dimension)
            )
            for row in range(dimension)
        )

    def rtd_marginal_generator(self) -> tuple[tuple[float, ...], ...]:
        """Return the frozen Erlang subgenerator in the same convention."""

        stage_count = self.topology.rtd_stage_count
        return tuple(
            tuple(
                (
                    -self.rtd_rate_s_inv
                    if destination == source
                    else self.rtd_rate_s_inv
                    if destination == source + 1
                    else 0.0
                )
                for source in range(stage_count)
            )
            for destination in range(stage_count)
        )

    def strong_lumpability_residual(self) -> tuple[tuple[float, ...], ...]:
        """Evaluate ``L G - R L`` for independent numerical audit.

        Construction already enforces the stronger typed-edge proof.  This
        diagnostic exposes the corresponding floating-point matrix identity
        without introducing a pass tolerance.
        """

        generator = self.generator_matrix()
        rtd_generator = self.rtd_marginal_generator()
        stage_count = self.topology.rtd_stage_count
        layer_count = self.topology.vertical_layer_count
        residual: list[tuple[float, ...]] = []
        for destination_stage in range(stage_count):
            row: list[float] = []
            destination_indices = range(
                destination_stage * layer_count,
                (destination_stage + 1) * layer_count,
            )
            for source_index, source_node in enumerate(self.nodes):
                left = math.fsum(
                    generator[destination_index][source_index]
                    for destination_index in destination_indices
                )
                source_stage = source_node.rtd_stage_id.value - 1
                right = rtd_generator[destination_stage][source_stage]
                row.append(math.fsum((left, -right)))
            residual.append(tuple(row))
        return tuple(residual)

    def evaluate_rates(
        self,
        inventories: Sequence[CarrierInventory],
    ) -> CarrierRateEvaluation:
        """Apply one identical operator to every declared extensive.

        This method evaluates rates only.  It performs no time integration and
        therefore cannot conceal a positivity limiter, timestep, or clamp.
        """

        if len(inventories) != len(self.nodes):
            raise ValueError("one carrier inventory is required for every operator node")
        if any(not isinstance(inventory, CarrierInventory) for inventory in inventories):
            raise TypeError("operator inventories must be CarrierInventory values")
        energy_datum_ids = {inventory.energy_datum_id for inventory in inventories}
        if energy_datum_ids != {self.energy_datum_id}:
            raise ValueError("carrier inventories do not match the operator energy datum authority")
        indices = {node: index for index, node in enumerate(self.nodes)}
        contributions: list[list[list[float]]] = [[[] for _ in range(7)] for _ in self.nodes]
        internal: list[InternalTransferRate] = []
        exits: list[ExitTransferRate] = []
        for transition in self.transitions:
            source_index = indices[transition.source]
            destination_index = indices[transition.destination]
            amount = tuple(
                transition.rate_s_inv * value for value in inventories[source_index].as_tuple()
            )
            transfer = CarrierInventoryRate.from_values(amount)
            internal.append(
                InternalTransferRate(
                    source=transition.source,
                    destination=transition.destination,
                    mechanism=transition.mechanism,
                    rate_s_inv=transition.rate_s_inv,
                    authority_id=transition.authority_id,
                    transferred=transfer,
                )
            )
            for component, value in enumerate(amount):
                contributions[source_index][component].append(-value)
                contributions[destination_index][component].append(value)
        for exit_rate in self.exits:
            source_index = indices[exit_rate.source]
            amount = tuple(
                exit_rate.rate_s_inv * value for value in inventories[source_index].as_tuple()
            )
            exits.append(
                ExitTransferRate(
                    source=exit_rate.source,
                    rate_s_inv=exit_rate.rate_s_inv,
                    authority_id=exit_rate.authority_id,
                    physical_port_id=exit_rate.physical_port_id,
                    transferred=CarrierInventoryRate.from_values(amount),
                )
            )
            for component, value in enumerate(amount):
                contributions[source_index][component].append(-value)
        node_rates = tuple(
            CarrierInventoryRate.from_values(
                math.fsum(component_terms) for component_terms in node_terms
            )
            for node_terms in contributions
        )
        return CarrierRateEvaluation(
            energy_datum_id=self.energy_datum_id,
            operator_definition_digest=self.definition_digest,
            residence_authority_digest=self.residence_authority.definition_digest,
            input_inventory_digest=_input_inventory_digest(
                nodes=self.nodes,
                inventories=inventories,
            ),
            nodes=self.nodes,
            node_rates=node_rates,
            internal_transfers=tuple(internal),
            exit_transfers=tuple(exits),
        )


def build_separable_joint_carrier_oracle(
    *,
    topology: PhysicalTrayTopology,
    residence_authority: RTDResidenceAuthority,
    energy_datum_id: str,
    exit_port_id: str,
    arm_authority: ArmRedistributionAuthority,
) -> ManufacturedJointCarrierOperator:
    """Build the owner-bounded separable oracle without selecting arm physics."""

    _validate_topology_runtime(topology)
    if not isinstance(residence_authority, RTDResidenceAuthority):
        raise TypeError("residence_authority must be an RTDResidenceAuthority")
    if not isinstance(arm_authority, ArmRedistributionAuthority):
        raise TypeError("arm_authority must be an ArmRedistributionAuthority")
    _require_nonblank_strings(
        "carrier builder energy datum and exit port IDs",
        energy_datum_id,
        exit_port_id,
    )
    if residence_authority.physical_tray_id != topology.physical_tray_id:
        raise ValueError("RTD residence authority belongs to another physical tray")
    if arm_authority.physical_tray_id != topology.physical_tray_id:
        raise ValueError("arm authority belongs to another physical tray")
    for transition in arm_authority.transitions:
        if transition.source_vertical_layer_id.value > topology.vertical_layer_count or (
            transition.destination_vertical_layer_id.value > topology.vertical_layer_count
        ):
            raise ValueError("arm transition references a layer outside the declared K")
    _validate_active_arm_graph(topology=topology, authority=arm_authority)

    nodes = tuple(
        JointCarrierNode(
            physical_tray_id=topology.physical_tray_id,
            vertical_layer_id=VerticalLayerId(k),
            rtd_stage_id=RTDStageId(n),
        )
        for n in range(1, topology.rtd_stage_count + 1)
        for k in range(1, topology.vertical_layer_count + 1)
    )
    node_by_coordinate = {
        (node.vertical_layer_id.value, node.rtd_stage_id.value): node for node in nodes
    }
    mean_residence_time_s = residence_authority.mean_residence_time_s
    rtd_rate_s_inv = _derived_stage_rate_s_inv(
        stage_count=topology.rtd_stage_count,
        mean_residence_time_s=mean_residence_time_s,
    )
    transitions: list[InternalCarrierTransition] = []
    exits: list[CarrierExit] = []
    for k in range(1, topology.vertical_layer_count + 1):
        for n in range(1, topology.rtd_stage_count):
            transitions.append(
                InternalCarrierTransition(
                    source=node_by_coordinate[k, n],
                    destination=node_by_coordinate[k, n + 1],
                    rate_s_inv=rtd_rate_s_inv,
                    mechanism=CarrierMechanism.RTD_ADVANCE,
                    authority_id=residence_authority.authority_id,
                )
            )
        exits.append(
            CarrierExit(
                source=node_by_coordinate[k, topology.rtd_stage_count],
                rate_s_inv=rtd_rate_s_inv,
                authority_id=residence_authority.authority_id,
                physical_port_id=exit_port_id,
            )
        )
    for n in range(1, topology.rtd_stage_count + 1):
        for arm_transition in arm_authority.transitions:
            transitions.append(
                InternalCarrierTransition(
                    source=node_by_coordinate[arm_transition.source_vertical_layer_id.value, n],
                    destination=node_by_coordinate[
                        arm_transition.destination_vertical_layer_id.value, n
                    ],
                    rate_s_inv=arm_transition.rate_s_inv,
                    mechanism=CarrierMechanism.ARM_REDISTRIBUTION,
                    authority_id=arm_authority.source_identity,
                )
            )
    return ManufacturedJointCarrierOperator(
        topology=topology,
        nodes=nodes,
        transitions=tuple(transitions),
        exits=tuple(exits),
        residence_authority=residence_authority,
        energy_datum_id=energy_datum_id,
        exit_port_id=exit_port_id,
        arm_authority=arm_authority,
    )


__all__ = [
    "ArmLayerTransition",
    "ArmRedistributionAuthority",
    "CarrierExit",
    "CarrierInventory",
    "CarrierInventoryRate",
    "CarrierMechanism",
    "CarrierRateEvaluation",
    "DryMatterDischargeBasis",
    "ExitTransferRate",
    "InternalCarrierTransition",
    "InternalTransferRate",
    "JointCarrierNode",
    "ManufacturedJointCarrierOperator",
    "RTDMoments",
    "RTDResidenceAuthority",
    "build_separable_joint_carrier_oracle",
    "tanks_in_series_moments",
    "validate_common_shaft_snapshot",
]
