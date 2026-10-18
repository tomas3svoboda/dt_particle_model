"""Bounded MECH-01 accepted draw/exit coupling for the six-tray QSC.

The accepted dry-matter flow is an input history, not a consequence of shaft
speed, valve position, or offered feed.  This seam integrates those accepted
histories over one already synchronized tower step and presents the resulting
debits to :mod:`accepted_k_cell_tower_router`.  The router remains the owner of
final-RTD eligibility, proportional all-K/cohort allocation, exact packet
splits, lineage, and seven-extensive conservation.

The active exit measure is the owner-selected PR-07 reference measure.  The
K-selective floor-band alternative is retained only as an inactive physical
bracket because D6 supplies no authenticated footprint, controller, or
speed-to-flow law.  PR-10 is applied only as an eligibility veto: positive
accepted draw through the five common-shaft sweep ports is impossible during
a positive-duration exact-zero-RPM interval.  The independently driven SP1
rotary airlock is outside that common-shaft veto.

This is engineering prequalification infrastructure.  It authenticates no
instrument, selects no discharge constitutive law, and advances no QSC-10,
physical, plant, production, or industrial-feasibility claim.
"""

from __future__ import annotations

import enum
import hashlib
import math
import struct
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import ClassVar, Iterable

from . import accepted_k_cell_tower_router as tower_router
from . import adjacent_k_arm_dispersion_oracle as arm_oracle
from . import k_cell_tray_host as tray_host
from . import mechanical_authority as mechanical
from . import qsc_mechanical_runtime_authority as mech_runtime

QSC_ACCEPTED_DRAW_SCHEMA_ID = "dtdc-core2-qsc-accepted-draw-exit-v1"
QSC_ACCEPTED_DRAW_SCHEMA_REVISION = 1
QSC_ACCEPTED_DRAW_LAW_SCHEMA_ID = "dtdc-core2-qsc-accepted-draw-exit-law-v1"
QSC_ACCEPTED_DRAW_LAW_SCHEMA_REVISION = 1
QSC_ACCEPTED_DRAW_OWNER_SOURCE_PATH = (
    "docs/GT_PS2_PRODUCTION_RULINGS_OWNER_DISPOSITION_2026-08-09.md"
)
QSC_ACCEPTED_DRAW_OWNER_SOURCE_CANONICAL_LF_SHA256 = (
    "1d8dd6e1e39f98a72c70f89b5510bded3037d7037bdae07ba733c6d12aecab81"
)
QSC_ACCEPTED_DRAW_EXPECTED_LAW_DEFINITION_DIGEST = (
    "sha256:24c630670a46c738bdc4d47d509e54373d80fe06f4e0f46f3a35d9cfe954bce5"
)
QSC_ACCEPTED_DRAW_CAUSAL_ORDER_LAW_ID = (
    "ACCEPTED_TN_PRIOR_THEN_SYNCHRONIZED_LOCAL_ENDPOINT_THEN_DRAW_ROUTE_THEN_ATOMIC_COMMIT_V1"
)
QSC_ACCEPTED_DRAW_HISTORY_LAW_ID = "PIECEWISE_CONSTANT_LOCALLY_ACCEPTED_DRY_MATTER_KG_S_INTEGRAL_V1"
QSC_ACCEPTED_DRAW_PR10_LAW_ID = "ZERO_RPM_VETOES_POSITIVE_DURATION_COMMON_SHAFT_SWEEP_DRAW_ONLY_V1"

_ROOT = Path(__file__).resolve().parents[3]


class QSCAcceptedDrawExitError(ValueError):
    """A draw law, history, route, or commit refused closed."""


class QSCAcceptedDrawExitReplayError(QSCAcceptedDrawExitError):
    """An already consumed draw/exit transaction was presented again."""


class QSCAcceptedDrawFootprint(enum.Enum):
    """Selected reference measure and still-open physical bracket."""

    PR07_ALL_FINAL_RTD_ALL_K_COHORT_PROPORTIONAL = "PR07_ALL_FINAL_RTD_ALL_K_COHORT_PROPORTIONAL"
    FLOOR_BAND_K_SELECTIVE_REQUIRES_GEOMETRY_OR_TRACER = (
        "FLOOR_BAND_K_SELECTIVE_REQUIRES_GEOMETRY_OR_TRACER"
    )


class QSCCommonShaftDrawDisposition(enum.Enum):
    SWEEP_POSITIVE_RPM_ELIGIBLE = "SWEEP_POSITIVE_RPM_ELIGIBLE"
    SWEEP_ZERO_RPM_NO_ACCEPTED_DRAW = "SWEEP_ZERO_RPM_NO_ACCEPTED_DRAW"
    INDEPENDENT_SP1_OUTSIDE_COMMON_SHAFT_VETO = "INDEPENDENT_SP1_OUTSIDE_COMMON_SHAFT_VETO"


def _require_identity(name: str, *values: str) -> None:
    if any(type(value) is not str or not value or value != value.strip() for value in values):
        raise QSCAcceptedDrawExitError(
            f"{name} must contain nonblank exact strings without surrounding whitespace"
        )


def _require_sha256(name: str, value: str) -> None:
    if (
        type(value) is not str
        or len(value) != 71
        or not value.startswith("sha256:")
        or any(character not in "0123456789abcdef" for character in value[7:])
    ):
        raise QSCAcceptedDrawExitError(f"{name} must be a canonical lowercase sha256 identity")


def _require_float(
    name: str,
    *values: float,
    nonnegative: bool = False,
) -> None:
    if any(type(value) is not float or not math.isfinite(value) for value in values):
        raise QSCAcceptedDrawExitError(f"{name} must contain finite exact binary64 values")
    if any(value == 0.0 and math.copysign(1.0, value) < 0.0 for value in values):
        raise QSCAcceptedDrawExitError(f"{name} rejects negative zero")
    if nonnegative and any(value < 0.0 for value in values):
        raise QSCAcceptedDrawExitError(f"{name} must be nonnegative")


def _canonical_lf_sha256(path: Path) -> str:
    raw = path.read_bytes()
    return hashlib.sha256(raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")).hexdigest()


def _digest_part(value: str | int | float | bool) -> tuple[bytes, bytes]:
    if type(value) is str:
        return b"s", value.encode("utf-8")
    if type(value) is bool:
        return b"b", b"1" if value else b"0"
    if type(value) is int:
        return b"i", str(value).encode("ascii")
    if type(value) is float:
        return b"f", struct.pack(">d", value)
    raise TypeError(f"unsupported accepted-draw digest part: {type(value).__name__}")


def _digest(domain: str, parts: Iterable[str | int | float | bool]) -> str:
    digest = hashlib.sha256()
    digest.update(b"DTDC-QSC-ACCEPTED-DRAW-EXIT-V1\0")
    domain_bytes = domain.encode("ascii")
    digest.update(len(domain_bytes).to_bytes(8, "big"))
    digest.update(domain_bytes)
    material = tuple(parts)
    digest.update(len(material).to_bytes(8, "big"))
    for value in material:
        tag, payload = _digest_part(value)
        digest.update(tag)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCAcceptedDrawExitLawAuthority:
    schema_id: str
    schema_revision: int
    owner_source_path: str
    owner_source_canonical_lf_sha256: str
    source_runtime_definition_digest: str
    source_configuration_canonical_lf_sha256: str
    source_router_id: str
    source_router_conservation_law_id: str
    solid_tray_order: tuple[str, ...]
    active_reference_footprint: QSCAcceptedDrawFootprint
    inactive_physical_bracket: QSCAcceptedDrawFootprint
    causal_order_law_id: str
    history_integration_law_id: str
    pr10_eligibility_law_id: str
    speed_to_flow_law_selected: bool
    offered_flow_authorizes_motion: bool

    declared_engineering_law: ClassVar[bool] = True
    accepted_flow_instrumentation_authenticated: ClassVar[bool] = False
    draw_footprint_physically_selected: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False
    mech_01: ClassVar[bool] = False
    f2: ClassVar[bool] = False
    f3: ClassVar[bool] = False
    f4: ClassVar[bool] = False
    f5: ClassVar[bool] = False
    f6: ClassVar[bool] = False
    f7: ClassVar[bool] = False
    f8a: ClassVar[bool] = False
    f8b: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    industry_feasible_operating_envelope: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if (
            self.schema_id != QSC_ACCEPTED_DRAW_LAW_SCHEMA_ID
            or type(self.schema_revision) is not int
            or self.schema_revision != QSC_ACCEPTED_DRAW_LAW_SCHEMA_REVISION
        ):
            raise QSCAcceptedDrawExitError("accepted-draw law has a foreign schema")
        _require_identity("accepted-draw source", self.owner_source_path)
        if self.owner_source_path != QSC_ACCEPTED_DRAW_OWNER_SOURCE_PATH:
            raise QSCAcceptedDrawExitError("accepted-draw law binds another owner source")
        if (
            self.owner_source_canonical_lf_sha256
            != QSC_ACCEPTED_DRAW_OWNER_SOURCE_CANONICAL_LF_SHA256
        ):
            raise QSCAcceptedDrawExitError("accepted-draw owner-source pin drifted")
        _require_sha256("source runtime", self.source_runtime_definition_digest)
        if (
            self.source_runtime_definition_digest
            != mech_runtime.QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST
            or self.source_configuration_canonical_lf_sha256
            != mech_runtime.QSC_MECH01_CONFIG_CANONICAL_LF_SHA256
        ):
            raise QSCAcceptedDrawExitError("accepted-draw law binds a foreign D6 runtime")
        if (
            self.source_router_id != tower_router.ROUTER_ID
            or self.source_router_conservation_law_id
            != tower_router.ROUTER_CONSERVATION_AUDIT_LAW_ID
        ):
            raise QSCAcceptedDrawExitError("accepted-draw law binds another router law")
        if self.solid_tray_order != mech_runtime.QSC_MECH01_SOLID_TRAY_ORDER:
            raise QSCAcceptedDrawExitError("accepted-draw law has a foreign solid order")
        if (
            self.active_reference_footprint
            is not QSCAcceptedDrawFootprint.PR07_ALL_FINAL_RTD_ALL_K_COHORT_PROPORTIONAL
            or self.inactive_physical_bracket
            is not QSCAcceptedDrawFootprint.FLOOR_BAND_K_SELECTIVE_REQUIRES_GEOMETRY_OR_TRACER
        ):
            raise QSCAcceptedDrawExitError("accepted-draw footprint disposition changed")
        if (
            self.causal_order_law_id != QSC_ACCEPTED_DRAW_CAUSAL_ORDER_LAW_ID
            or self.history_integration_law_id != QSC_ACCEPTED_DRAW_HISTORY_LAW_ID
            or self.pr10_eligibility_law_id != QSC_ACCEPTED_DRAW_PR10_LAW_ID
        ):
            raise QSCAcceptedDrawExitError("accepted-draw causal or eligibility law changed")
        if self.speed_to_flow_law_selected is not False:
            raise QSCAcceptedDrawExitError("D6 supplies no speed-to-flow law")
        if self.offered_flow_authorizes_motion is not False:
            raise QSCAcceptedDrawExitError("offered flow cannot authorize solids motion")

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-accepted-draw-exit-law-authority",
            (
                self.schema_id,
                self.schema_revision,
                self.owner_source_path,
                self.owner_source_canonical_lf_sha256,
                self.source_runtime_definition_digest,
                self.source_configuration_canonical_lf_sha256,
                self.source_router_id,
                self.source_router_conservation_law_id,
                *self.solid_tray_order,
                self.active_reference_footprint.value,
                self.inactive_physical_bracket.value,
                self.causal_order_law_id,
                self.history_integration_law_id,
                self.pr10_eligibility_law_id,
                self.speed_to_flow_law_selected,
                self.offered_flow_authorizes_motion,
            ),
        )


QSC_ACCEPTED_DRAW_EXIT_LAW_AUTHORITY = QSCAcceptedDrawExitLawAuthority(
    schema_id=QSC_ACCEPTED_DRAW_LAW_SCHEMA_ID,
    schema_revision=QSC_ACCEPTED_DRAW_LAW_SCHEMA_REVISION,
    owner_source_path=QSC_ACCEPTED_DRAW_OWNER_SOURCE_PATH,
    owner_source_canonical_lf_sha256=QSC_ACCEPTED_DRAW_OWNER_SOURCE_CANONICAL_LF_SHA256,
    source_runtime_definition_digest=mech_runtime.QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST,
    source_configuration_canonical_lf_sha256=(mech_runtime.QSC_MECH01_CONFIG_CANONICAL_LF_SHA256),
    source_router_id=tower_router.ROUTER_ID,
    source_router_conservation_law_id=tower_router.ROUTER_CONSERVATION_AUDIT_LAW_ID,
    solid_tray_order=mech_runtime.QSC_MECH01_SOLID_TRAY_ORDER,
    active_reference_footprint=(
        QSCAcceptedDrawFootprint.PR07_ALL_FINAL_RTD_ALL_K_COHORT_PROPORTIONAL
    ),
    inactive_physical_bracket=(
        QSCAcceptedDrawFootprint.FLOOR_BAND_K_SELECTIVE_REQUIRES_GEOMETRY_OR_TRACER
    ),
    causal_order_law_id=QSC_ACCEPTED_DRAW_CAUSAL_ORDER_LAW_ID,
    history_integration_law_id=QSC_ACCEPTED_DRAW_HISTORY_LAW_ID,
    pr10_eligibility_law_id=QSC_ACCEPTED_DRAW_PR10_LAW_ID,
    speed_to_flow_law_selected=False,
    offered_flow_authorizes_motion=False,
)


def validate_qsc_accepted_draw_exit_law_authority(
    authority: QSCAcceptedDrawExitLawAuthority,
) -> str:
    if type(authority) is not QSCAcceptedDrawExitLawAuthority:
        raise TypeError("accepted-draw law must be an exact typed authority")
    QSCAcceptedDrawExitLawAuthority.__post_init__(authority)
    source = (_ROOT / authority.owner_source_path).resolve()
    if not source.is_relative_to(_ROOT.resolve()) or not source.is_file():
        raise QSCAcceptedDrawExitError("accepted-draw owner source is missing")
    if _canonical_lf_sha256(source) != authority.owner_source_canonical_lf_sha256:
        raise QSCAcceptedDrawExitError("accepted-draw owner source bytes changed")
    digest = authority.definition_digest
    if digest != QSC_ACCEPTED_DRAW_EXPECTED_LAW_DEFINITION_DIGEST:
        raise QSCAcceptedDrawExitError("accepted-draw law differs from its expected identity")
    return digest


def qsc_accepted_draw_flow_history_digest(
    *,
    physical_tray_id: str,
    segments: tuple[tray_host.AcceptedDryMatterFlowSegment, ...],
) -> str:
    _require_identity("flow-history tray", physical_tray_id)
    if type(segments) is not tuple or not segments:
        raise QSCAcceptedDrawExitError("accepted draw requires a nonempty exact history")
    if any(type(item) is not tray_host.AcceptedDryMatterFlowSegment for item in segments):
        raise TypeError("accepted draw history contains a foreign segment")
    authority_ids = tuple(item.accepted_flow_authority_id for item in segments)
    if len(set(authority_ids)) != 1:
        raise QSCAcceptedDrawExitError(
            "one tray step must bind one accepted flow/controller authority"
        )
    return _digest(
        "qsc-accepted-draw-flow-history",
        (
            physical_tray_id,
            authority_ids[0],
            len(segments),
            *(
                value
                for item in segments
                for value in (
                    item.start_time_s,
                    item.end_time_s,
                    item.accepted_flow_authority_id,
                    item.accepted_dry_matter_discharge_kg_s,
                    item.offered_dry_matter_discharge_kg_s,
                    item.discharge_basis.value,
                )
            ),
        ),
    )


def qsc_accepted_draw_shaft_history_digest(
    segments: tuple[arm_oracle.ManufacturedAcceptedShaftSegment, ...],
) -> str:
    if type(segments) is not tuple or not segments:
        raise QSCAcceptedDrawExitError("accepted draw requires a nonempty shaft history")
    if any(type(item) is not arm_oracle.ManufacturedAcceptedShaftSegment for item in segments):
        raise TypeError("accepted shaft history contains a foreign segment")
    return _digest(
        "qsc-accepted-draw-shaft-history",
        (len(segments), *(item.definition_digest for item in segments)),
    )


def _validate_contiguous_flow_history(
    *,
    physical_tray_id: str,
    segments: tuple[tray_host.AcceptedDryMatterFlowSegment, ...],
    start_time_s: float,
    end_time_s: float,
) -> None:
    qsc_accepted_draw_flow_history_digest(
        physical_tray_id=physical_tray_id,
        segments=segments,
    )
    if segments[0].start_time_s != start_time_s or segments[-1].end_time_s != end_time_s:
        raise QSCAcceptedDrawExitError("accepted flow history does not span the tower step")
    if any(
        left.end_time_s != right.start_time_s
        for left, right in zip(segments, segments[1:], strict=False)
    ):
        raise QSCAcceptedDrawExitError("accepted flow history must be exactly contiguous")


def _validate_shaft_history(
    *,
    segments: tuple[arm_oracle.ManufacturedAcceptedShaftSegment, ...],
    runtime_authority: mech_runtime.QSCMechanicalRuntimeAuthority,
    start_time_s: float,
    end_time_s: float,
) -> None:
    qsc_accepted_draw_shaft_history_digest(segments)
    common_digest = runtime_authority.common_shaft.definition_digest
    accepted_history_digest = (
        runtime_authority.common_shaft.accepted_shaft_history.definition_digest
    )
    for item in segments:
        arm_oracle.ManufacturedAcceptedShaftSegment.__post_init__(item)
        if (
            item.common_shaft_definition_digest != common_digest
            or item.accepted_history_definition_digest != accepted_history_digest
        ):
            raise QSCAcceptedDrawExitError("shaft history binds another common shaft or history")
    if tuple(segments) != tuple(sorted(segments, key=lambda item: item.time_start_s)):
        raise QSCAcceptedDrawExitError("shaft history must use chronological order")
    if segments[0].time_start_s > start_time_s or segments[-1].time_end_s < end_time_s:
        raise QSCAcceptedDrawExitError("shaft history does not cover the tower step")
    for left, right in zip(segments, segments[1:], strict=False):
        if left.time_end_s != right.time_start_s:
            raise QSCAcceptedDrawExitError("shaft history must be contiguous without overlap")


def _integrated_mass(
    segments: tuple[tray_host.AcceptedDryMatterFlowSegment, ...],
    *,
    accepted: bool,
) -> float:
    terms = tuple(
        (
            item.accepted_dry_matter_discharge_kg_s
            if accepted
            else item.offered_dry_matter_discharge_kg_s
        )
        * item.duration_s
        for item in segments
    )
    if any(not math.isfinite(value) for value in terms):
        raise QSCAcceptedDrawExitError("integrated draw mass is not representable")
    result = math.fsum(terms)
    _require_float("integrated draw mass", result, nonnegative=True)
    return result


def _zero_rpm_overlap_s(
    flow_segments: tuple[tray_host.AcceptedDryMatterFlowSegment, ...],
    shaft_segments: tuple[arm_oracle.ManufacturedAcceptedShaftSegment, ...],
) -> tuple[float, bool]:
    zero_overlap: list[float] = []
    positive_accepted_on_zero = False
    for flow in flow_segments:
        for shaft in shaft_segments:
            overlap_start = max(flow.start_time_s, shaft.time_start_s)
            overlap_end = min(flow.end_time_s, shaft.time_end_s)
            if overlap_end <= overlap_start:
                continue
            # With the selected linear segment semantics, exact zero occupies a
            # positive-duration interval only when both endpoints are zero.
            if shaft.rpm_start == 0.0 and shaft.rpm_end == 0.0:
                duration = overlap_end - overlap_start
                zero_overlap.append(duration)
                if flow.accepted_dry_matter_discharge_kg_s > 0.0:
                    positive_accepted_on_zero = True
    duration = math.fsum(zero_overlap)
    _require_float("zero-RPM overlap duration", duration, nonnegative=True)
    return duration, positive_accepted_on_zero


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCAcceptedDrawExitRequest:
    attempt_id: str
    start_time_s: float
    end_time_s: float
    accepted_prior_trays: tuple[tray_host.AcceptedKCellTrayState, ...]
    routing_endpoint_trays: tuple[tray_host.AcceptedKCellTrayState, ...]
    flow_histories: tuple[tuple[tray_host.AcceptedDryMatterFlowSegment, ...], ...]
    expected_flow_history_digests: tuple[str, ...]
    accepted_shaft_segments: tuple[arm_oracle.ManufacturedAcceptedShaftSegment, ...]
    expected_shaft_history_digest: str
    closing_layer_cursors: tuple[int, ...]
    feed_packets: tuple[tower_router.AcceptedFeedPacket, ...]
    #: Owner-issued v3 closing-exchange authority threaded into the router
    #: request (GT_PS2_WALLED_FLOOR_AND_V3_APPROVAL_RULING_2026-08-29);
    #: absent means a closing-layer representability exhaustion keeps its
    #: typed refusal (no policy, no v3).
    k_share_closing_exchange_policy: "tower_router.KShareClosingExchangePolicy | None" = None

    physically_qualifying: ClassVar[bool] = False
    accepted_histories_authenticated: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("accepted-draw attempt", self.attempt_id)
        _require_float("accepted-draw step time", self.start_time_s, self.end_time_s)
        if self.start_time_s < 0.0 or self.end_time_s <= self.start_time_s:
            raise QSCAcceptedDrawExitError("accepted-draw step must advance time")
        order = mech_runtime.QSC_MECH01_SOLID_TRAY_ORDER
        for label, trays, expected_time in (
            ("accepted prior", self.accepted_prior_trays, self.start_time_s),
            ("routing endpoint", self.routing_endpoint_trays, self.end_time_s),
        ):
            if type(trays) is not tuple or len(trays) != len(order):
                raise TypeError(f"{label} must contain six exact tray states")
            if any(type(item) is not tray_host.AcceptedKCellTrayState for item in trays):
                raise TypeError(f"{label} contains a foreign tray state")
            if tuple(item.physical_tray_id for item in trays) != order:
                raise QSCAcceptedDrawExitError(f"{label} differs from solid tray order")
            if any(item.time_s != expected_time for item in trays):
                raise QSCAcceptedDrawExitError(f"{label} does not share its declared time")
        for prior, endpoint in zip(
            self.accepted_prior_trays,
            self.routing_endpoint_trays,
            strict=True,
        ):
            if (
                endpoint.topology.digest != prior.topology.digest
                or endpoint.configuration_digest != prior.configuration_digest
                or endpoint.codec_registry_digest != prior.codec_registry_digest
                or endpoint.auditor_identity_digest != prior.auditor_identity_digest
                or endpoint.energy_datum_id != prior.energy_datum_id
                or endpoint.rtd_schedule_digest != prior.rtd_schedule_digest
                or endpoint.revision <= prior.revision
            ):
                raise QSCAcceptedDrawExitError(
                    "routing endpoint is not a forward local candidate of its accepted prior"
                )
        if type(self.flow_histories) is not tuple or len(self.flow_histories) != len(order):
            raise TypeError("accepted draw requires six exact flow histories")
        if type(self.expected_flow_history_digests) is not tuple or len(
            self.expected_flow_history_digests
        ) != len(order):
            raise TypeError("accepted draw requires six caller-held flow-history pins")
        observed = []
        for tray_id, history in zip(order, self.flow_histories, strict=True):
            _validate_contiguous_flow_history(
                physical_tray_id=tray_id,
                segments=history,
                start_time_s=self.start_time_s,
                end_time_s=self.end_time_s,
            )
            observed.append(
                qsc_accepted_draw_flow_history_digest(
                    physical_tray_id=tray_id,
                    segments=history,
                )
            )
        if tuple(observed) != self.expected_flow_history_digests:
            raise QSCAcceptedDrawExitError("accepted flow histories differ from caller-held pins")
        _require_sha256("expected shaft history", self.expected_shaft_history_digest)
        if qsc_accepted_draw_shaft_history_digest(self.accepted_shaft_segments) != (
            self.expected_shaft_history_digest
        ):
            raise QSCAcceptedDrawExitError("accepted shaft history differs from caller-held pin")
        if (
            type(self.closing_layer_cursors) is not tuple
            or len(self.closing_layer_cursors) != len(order)
            or any(type(item) is not int or item <= 0 for item in self.closing_layer_cursors)
        ):
            raise QSCAcceptedDrawExitError("accepted draw requires six positive K cursors")
        if type(self.feed_packets) is not tuple or any(
            type(item) is not tower_router.AcceptedFeedPacket for item in self.feed_packets
        ):
            raise TypeError("accepted draw feed contains a foreign packet")
        if self.k_share_closing_exchange_policy is not None:
            if (
                type(self.k_share_closing_exchange_policy)
                is not tower_router.KShareClosingExchangePolicy
            ):
                raise TypeError("accepted draw closing-exchange policy has a foreign exact type")
            self.k_share_closing_exchange_policy.__post_init__()

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-accepted-draw-exit-request",
            (
                self.attempt_id,
                self.start_time_s,
                self.end_time_s,
                *(item.state_digest for item in self.accepted_prior_trays),
                *(item.state_digest for item in self.routing_endpoint_trays),
                *self.expected_flow_history_digests,
                self.expected_shaft_history_digest,
                *self.closing_layer_cursors,
                (
                    "k-share-closing-exchange-policy:absent"
                    if self.k_share_closing_exchange_policy is None
                    else (
                        "k-share-closing-exchange-policy:"
                        f"{self.k_share_closing_exchange_policy.policy_id};"
                        f"{self.k_share_closing_exchange_policy.candidate_order_law};"
                        f"{self.k_share_closing_exchange_policy.ruling_provenance};"
                        f"{self.k_share_closing_exchange_policy.cumulative_signed_drift_tripwire_ulp}"
                    )
                ),
                len(self.feed_packets),
                *(
                    value
                    for item in self.feed_packets
                    for value in (
                        item.entry.packet_id,
                        hashlib.sha256(item.entry.payload_bytes).hexdigest(),
                        item.entry.weight.representative_particles,
                    )
                ),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCAcceptedTrayDrawLedger:
    physical_tray_id: str
    draw_device_id: str
    draw_geometry_definition_digest: str
    accepted_flow_controller_authority_id: str
    flow_history_definition_digest: str
    accepted_dry_matter_kg: float
    offered_dry_matter_kg_record_only: float
    eligible_final_rtd_dry_matter_kg: float
    represented_dry_matter_debit_kg: float
    zero_rpm_overlap_s: float
    disposition: QSCCommonShaftDrawDisposition
    active_reference_footprint: QSCAcceptedDrawFootprint
    inactive_physical_bracket: QSCAcceptedDrawFootprint
    boundary_split_count: int

    accepted_not_offered_is_authoritative: ClassVar[bool] = True
    accepted_debit_is_booked_boundary_flow: ClassVar[bool] = True
    flat_packet_debit_is_diagnostic_only: ClassVar[bool] = True
    whole_router_exact_zero_is_commit_authority: ClassVar[bool] = True
    speed_to_flow_law_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity(
            "accepted-draw ledger identity",
            self.physical_tray_id,
            self.draw_device_id,
            self.accepted_flow_controller_authority_id,
        )
        _require_sha256("draw geometry", self.draw_geometry_definition_digest)
        _require_sha256("flow history", self.flow_history_definition_digest)
        _require_float(
            "accepted-draw ledger mass/time",
            self.accepted_dry_matter_kg,
            self.offered_dry_matter_kg_record_only,
            self.eligible_final_rtd_dry_matter_kg,
            self.represented_dry_matter_debit_kg,
            self.zero_rpm_overlap_s,
            nonnegative=True,
        )
        if self.accepted_dry_matter_kg > self.eligible_final_rtd_dry_matter_kg:
            raise QSCAcceptedDrawExitError("accepted draw exceeds final-RTD eligibility")
        if type(self.disposition) is not QSCCommonShaftDrawDisposition:
            raise TypeError("accepted-draw disposition has a foreign type")
        if (
            self.active_reference_footprint
            is not QSCAcceptedDrawFootprint.PR07_ALL_FINAL_RTD_ALL_K_COHORT_PROPORTIONAL
            or self.inactive_physical_bracket
            is not QSCAcceptedDrawFootprint.FLOOR_BAND_K_SELECTIVE_REQUIRES_GEOMETRY_OR_TRACER
        ):
            raise QSCAcceptedDrawExitError("accepted-draw ledger changed footprint meaning")
        if type(self.boundary_split_count) is not int or self.boundary_split_count < 0:
            raise QSCAcceptedDrawExitError("boundary split count must be nonnegative")

    @property
    def exact_flattened_minus_accepted_dry_matter_kg(self) -> Fraction:
        """Return the exact binary64 association diagnostic ``F - A``."""

        return Fraction.from_float(self.represented_dry_matter_debit_kg) - Fraction.from_float(
            self.accepted_dry_matter_kg
        )

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-accepted-tray-draw-ledger",
            (
                self.physical_tray_id,
                self.draw_device_id,
                self.draw_geometry_definition_digest,
                self.accepted_flow_controller_authority_id,
                self.flow_history_definition_digest,
                self.accepted_dry_matter_kg,
                self.offered_dry_matter_kg_record_only,
                self.eligible_final_rtd_dry_matter_kg,
                self.represented_dry_matter_debit_kg,
                self.zero_rpm_overlap_s,
                self.disposition.value,
                self.active_reference_footprint.value,
                self.inactive_physical_bracket.value,
                self.boundary_split_count,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCAcceptedDrawExitRejection:
    attempt_id: str
    reason: str
    rollback_trays: tuple[tray_host.AcceptedKCellTrayState, ...]

    exact_accepted_prior_rollback: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False


_TRANSACTION_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedQSCAcceptedDrawExit:
    request: QSCAcceptedDrawExitRequest
    law_definition_digest: str
    runtime_definition_digest: str
    request_definition_digest: str
    tray_ledgers: tuple[QSCAcceptedTrayDrawLedger, ...]
    routing_transaction: tower_router.ValidatedSerialTowerRouting
    _evaluator_token: object = field(repr=False, compare=False)
    _seal: object = field(repr=False, compare=False)

    bounded_engineering_draw: ClassVar[bool] = True
    accepted_flow_instrumentation_authenticated: ClassVar[bool] = False
    draw_footprint_physically_selected: ClassVar[bool] = False
    durable_atomic_persistence_implemented: ClassVar[bool] = False
    mech_01: ClassVar[bool] = False
    f2: ClassVar[bool] = False
    f3: ClassVar[bool] = False
    f4: ClassVar[bool] = False
    f5: ClassVar[bool] = False
    f6: ClassVar[bool] = False
    f7: ClassVar[bool] = False
    f8a: ClassVar[bool] = False
    f8b: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _TRANSACTION_SEAL:
            raise TypeError("accepted draw transactions are evaluator-issued only")
        if type(self.request) is not QSCAcceptedDrawExitRequest:
            raise TypeError("accepted draw transaction has a foreign request")
        if type(self.routing_transaction) is not tower_router.ValidatedSerialTowerRouting:
            raise TypeError("accepted draw transaction has a foreign router transaction")
        for name in (
            "law_definition_digest",
            "runtime_definition_digest",
            "request_definition_digest",
        ):
            _require_sha256(name, getattr(self, name))
        if self.request_definition_digest != self.request.definition_digest:
            raise QSCAcceptedDrawExitError("accepted draw request identity changed")
        if (
            type(self.tray_ledgers) is not tuple
            or len(self.tray_ledgers) != len(mech_runtime.QSC_MECH01_SOLID_TRAY_ORDER)
            or any(type(item) is not QSCAcceptedTrayDrawLedger for item in self.tray_ledgers)
        ):
            raise TypeError("accepted draw transaction requires six typed ledgers")
        if tuple(item.physical_tray_id for item in self.tray_ledgers) != (
            mech_runtime.QSC_MECH01_SOLID_TRAY_ORDER
        ):
            raise QSCAcceptedDrawExitError("accepted draw ledgers have a foreign order")
        whole = self.routing_transaction.whole_router_ledger
        if not whole.passed or any(value != 0 for value in whole.conservation_residual):
            raise QSCAcceptedDrawExitError("accepted draw lacks exact seven-extensive closure")
        for draw_ledger, exit_ledger in zip(
            self.tray_ledgers,
            self.routing_transaction.exit_ledgers,
            strict=True,
        ):
            if (
                draw_ledger.physical_tray_id != exit_ledger.physical_tray_id
                or draw_ledger.accepted_dry_matter_kg != exit_ledger.accepted_dry_matter_debit_kg
                or draw_ledger.represented_dry_matter_debit_kg
                != exit_ledger.represented_dry_matter_debit_kg
                or draw_ledger.boundary_split_count
                != sum(item.boundary_split_count for item in exit_ledger.k_shares)
            ):
                raise QSCAcceptedDrawExitError(
                    "accepted draw ledger differs from its exact router exit booking"
                )

    @property
    def definition_digest(self) -> str:
        return _digest(
            "validated-qsc-accepted-draw-exit",
            (
                self.law_definition_digest,
                self.runtime_definition_digest,
                self.request_definition_digest,
                *(item.definition_digest for item in self.tray_ledgers),
                self.routing_transaction.transaction_digest,
            ),
        )


class QSCAcceptedDrawExitEvaluator:
    """One in-memory evaluate/commit owner for a coupled accepted draw step."""

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False
    durable_atomic_persistence_implemented: ClassVar[bool] = False
    nonmutating_commit_availability_preflight: ClassVar[bool] = True
    commit_availability_preflight_reserves_authority: ClassVar[bool] = False

    def __init__(
        self,
        *,
        population_router: tower_router.AcceptedKCellTowerRouter,
        runtime_authority: mech_runtime.QSCMechanicalRuntimeAuthority,
        law_authority: QSCAcceptedDrawExitLawAuthority,
        expected_runtime_definition_digest: str,
        expected_law_definition_digest: str,
    ) -> None:
        if type(population_router) is not tower_router.AcceptedKCellTowerRouter:
            raise TypeError("accepted draw requires an exact population router")
        if type(runtime_authority) is not mech_runtime.QSCMechanicalRuntimeAuthority:
            raise TypeError("accepted draw requires an exact MECH-01A runtime")
        if type(law_authority) is not QSCAcceptedDrawExitLawAuthority:
            raise TypeError("accepted draw requires an exact law authority")
        _require_sha256("expected runtime", expected_runtime_definition_digest)
        _require_sha256("expected accepted-draw law", expected_law_definition_digest)
        observed_runtime = mech_runtime.validate_qsc_mechanical_runtime_authority(runtime_authority)
        observed_law = validate_qsc_accepted_draw_exit_law_authority(law_authority)
        if (
            observed_runtime != expected_runtime_definition_digest
            or observed_law != expected_law_definition_digest
        ):
            raise QSCAcceptedDrawExitError("caller-held runtime or draw-law pin differs")
        self._router = population_router
        self._runtime = runtime_authority
        self._law = law_authority
        self._runtime_digest = observed_runtime
        self._law_digest = observed_law
        self._token = object()
        self._issued: dict[int, ValidatedQSCAcceptedDrawExit] = {}
        self._consumed: set[str] = set()

    def _validate_live_authorities(self) -> None:
        if mech_runtime.validate_qsc_mechanical_runtime_authority(self._runtime) != (
            self._runtime_digest
        ):
            raise QSCAcceptedDrawExitError("operative MECH-01A runtime changed")
        if validate_qsc_accepted_draw_exit_law_authority(self._law) != self._law_digest:
            raise QSCAcceptedDrawExitError("operative accepted-draw law changed")

    def evaluate(
        self,
        request: QSCAcceptedDrawExitRequest,
        *,
        admission_planners: tuple[tower_router.ManufacturedAdmissionPlanner, ...],
    ) -> ValidatedQSCAcceptedDrawExit | QSCAcceptedDrawExitRejection:
        if type(request) is not QSCAcceptedDrawExitRequest:
            raise TypeError("accepted draw evaluator requires an exact request")
        try:
            self._validate_live_authorities()
            QSCAcceptedDrawExitRequest.__post_init__(request)
            _validate_shaft_history(
                segments=request.accepted_shaft_segments,
                runtime_authority=self._runtime,
                start_time_s=request.start_time_s,
                end_time_s=request.end_time_s,
            )
            draw_map = {
                item.physical_tray_id: item.draw_geometry
                for item in self._runtime.common_shaft.trays
            }
            integrated: list[tuple[float, float, float, QSCCommonShaftDrawDisposition]] = []
            debits: list[tower_router.AcceptedTrayExitDebit] = []
            for tray_id, history, cursor in zip(
                self._law.solid_tray_order,
                request.flow_histories,
                request.closing_layer_cursors,
                strict=True,
            ):
                draw = draw_map[tray_id]
                accepted_kg = _integrated_mass(history, accepted=True)
                offered_kg = _integrated_mass(history, accepted=False)
                zero_overlap_s, positive_on_zero = _zero_rpm_overlap_s(
                    history,
                    request.accepted_shaft_segments,
                )
                if draw.conveyance_domain is (
                    mechanical.SolidsConveyanceDomain.COMMON_SHAFT_SWEEP_CONVEYED
                ):
                    if positive_on_zero:
                        raise QSCAcceptedDrawExitError(
                            f"PR-10 veto: {tray_id} reports positive accepted sweep draw "
                            "during a positive-duration exact-zero-RPM interval"
                        )
                    disposition = (
                        QSCCommonShaftDrawDisposition.SWEEP_ZERO_RPM_NO_ACCEPTED_DRAW
                        if zero_overlap_s > 0.0
                        else QSCCommonShaftDrawDisposition.SWEEP_POSITIVE_RPM_ELIGIBLE
                    )
                elif (
                    tray_id == "SP1"
                    and draw.conveyance_domain
                    is mechanical.SolidsConveyanceDomain.INDEPENDENT_DRIVE
                    and draw.actuation_domain is mechanical.DrawActuationDomain.INDEPENDENT_ACTUATOR
                ):
                    disposition = (
                        QSCCommonShaftDrawDisposition.INDEPENDENT_SP1_OUTSIDE_COMMON_SHAFT_VETO
                    )
                else:
                    raise QSCAcceptedDrawExitError(
                        f"{tray_id} draw domain differs from the ruled D6 configuration"
                    )
                integrated.append((accepted_kg, offered_kg, zero_overlap_s, disposition))
                debits.append(
                    tower_router.AcceptedTrayExitDebit(
                        physical_tray_id=tray_id,
                        accepted_dry_matter_kg=accepted_kg,
                        closing_layer_cursor=cursor,
                    )
                )
            route_request = tower_router.SerialTowerRoutingRequest(
                attempt_id=f"{request.attempt_id}:accepted-draw-route",
                event_time_s=request.end_time_s,
                prior_trays=request.routing_endpoint_trays,
                exit_debits=tuple(debits),
                feed_packets=request.feed_packets,
                k_share_closing_exchange_policy=(request.k_share_closing_exchange_policy),
            )
            routed = self._router.evaluate(
                route_request,
                admission_planners=admission_planners,
            )
            if type(routed) is tower_router.TowerRoutingRejection:
                raise QSCAcceptedDrawExitError(f"population router refused: {routed.reason}")
            if type(routed) is not tower_router.ValidatedSerialTowerRouting:
                raise TypeError("population router returned a foreign result")
            if not routed.whole_router_ledger.passed or any(
                value != 0 for value in routed.whole_router_ledger.conservation_residual
            ):
                raise QSCAcceptedDrawExitError(
                    "population draw lacks exact seven-extensive closure"
                )
            ledgers = []
            for tray_id, history, values, exit_ledger in zip(
                self._law.solid_tray_order,
                request.flow_histories,
                integrated,
                routed.exit_ledgers,
                strict=True,
            ):
                accepted_kg, offered_kg, zero_overlap_s, disposition = values
                draw = draw_map[tray_id]
                ledgers.append(
                    QSCAcceptedTrayDrawLedger(
                        physical_tray_id=tray_id,
                        draw_device_id=draw.draw_device_id,
                        draw_geometry_definition_digest=draw.definition_digest,
                        accepted_flow_controller_authority_id=(
                            history[0].accepted_flow_authority_id
                        ),
                        flow_history_definition_digest=(
                            qsc_accepted_draw_flow_history_digest(
                                physical_tray_id=tray_id,
                                segments=history,
                            )
                        ),
                        accepted_dry_matter_kg=accepted_kg,
                        offered_dry_matter_kg_record_only=offered_kg,
                        eligible_final_rtd_dry_matter_kg=(exit_ledger.eligible_dry_matter_kg),
                        represented_dry_matter_debit_kg=(
                            exit_ledger.represented_dry_matter_debit_kg
                        ),
                        zero_rpm_overlap_s=zero_overlap_s,
                        disposition=disposition,
                        active_reference_footprint=(self._law.active_reference_footprint),
                        inactive_physical_bracket=(self._law.inactive_physical_bracket),
                        boundary_split_count=sum(
                            item.boundary_split_count for item in exit_ledger.k_shares
                        ),
                    )
                )
            transaction = ValidatedQSCAcceptedDrawExit(
                request=request,
                law_definition_digest=self._law_digest,
                runtime_definition_digest=self._runtime_digest,
                request_definition_digest=request.definition_digest,
                tray_ledgers=tuple(ledgers),
                routing_transaction=routed,
                _evaluator_token=self._token,
                _seal=_TRANSACTION_SEAL,
            )
            self._issued[id(transaction)] = transaction
            return transaction
        except (
            QSCAcceptedDrawExitError,
            tower_router.KCellTowerRouterError,
            TypeError,
            ValueError,
            OverflowError,
        ) as error:
            return QSCAcceptedDrawExitRejection(
                attempt_id=request.attempt_id,
                reason=f"{type(error).__name__}: {error}",
                rollback_trays=request.accepted_prior_trays,
            )

    def verify_issued_transaction(
        self,
        transaction: ValidatedQSCAcceptedDrawExit,
    ) -> ValidatedQSCAcceptedDrawExit:
        if (
            type(transaction) is not ValidatedQSCAcceptedDrawExit
            or transaction._seal is not _TRANSACTION_SEAL
            or transaction._evaluator_token is not self._token
            or self._issued.get(id(transaction)) is not transaction
        ):
            raise QSCAcceptedDrawExitError(
                "verification requires the exact evaluator-issued draw transaction"
            )
        self._validate_live_authorities()
        ValidatedQSCAcceptedDrawExit.__post_init__(transaction)
        QSCAcceptedDrawExitRequest.__post_init__(transaction.request)
        self._router.verify_issued_transaction(transaction.routing_transaction)
        expected_debits = tuple(
            _integrated_mass(history, accepted=True)
            for history in transaction.request.flow_histories
        )
        if tuple(item.accepted_dry_matter_kg for item in transaction.tray_ledgers) != (
            expected_debits
        ):
            raise QSCAcceptedDrawExitError("accepted draw debit changed after evaluation")
        return transaction

    def commit(
        self,
        current_accepted_trays: tuple[tray_host.AcceptedKCellTrayState, ...],
        transaction: ValidatedQSCAcceptedDrawExit,
    ) -> tuple[tray_host.AcceptedKCellTrayState, ...]:
        self.verify_commit_available(current_accepted_trays, transaction)
        token = transaction.definition_digest
        candidate = self._router.commit(
            transaction.request.routing_endpoint_trays,
            transaction.routing_transaction,
        )
        self._consumed.add(token)
        return candidate

    def verify_commit_available(
        self,
        current_accepted_trays: tuple[tray_host.AcceptedKCellTrayState, ...],
        transaction: ValidatedQSCAcceptedDrawExit,
    ) -> ValidatedQSCAcceptedDrawExit:
        """Re-audit one draw commit without consuming draw or router authority.

        This is the read-only child preflight needed by a future composite
        publisher.  It authenticates the exact evaluator-issued draw, checks
        the outer accepted prior, and delegates the synchronized routing
        endpoint check to the router's own nonmutating availability audit.
        It deliberately does not reserve either authority and does not make a
        durable-persistence claim.
        """

        self.verify_issued_transaction(transaction)
        token = transaction.definition_digest
        if token in self._consumed:
            raise QSCAcceptedDrawExitReplayError("accepted draw transaction was already consumed")
        if (
            type(current_accepted_trays) is not tuple
            or len(current_accepted_trays) != len(transaction.request.accepted_prior_trays)
            or any(
                current is not prior
                for current, prior in zip(
                    current_accepted_trays,
                    transaction.request.accepted_prior_trays,
                    strict=True,
                )
            )
        ):
            raise QSCAcceptedDrawExitError("accepted draw transaction is stale or foreign")
        self._router.verify_commit_available(
            transaction.request.routing_endpoint_trays,
            transaction.routing_transaction,
        )
        return transaction


__all__ = (
    "QSC_ACCEPTED_DRAW_CAUSAL_ORDER_LAW_ID",
    "QSC_ACCEPTED_DRAW_EXPECTED_LAW_DEFINITION_DIGEST",
    "QSC_ACCEPTED_DRAW_EXIT_LAW_AUTHORITY",
    "QSC_ACCEPTED_DRAW_HISTORY_LAW_ID",
    "QSC_ACCEPTED_DRAW_OWNER_SOURCE_CANONICAL_LF_SHA256",
    "QSC_ACCEPTED_DRAW_OWNER_SOURCE_PATH",
    "QSC_ACCEPTED_DRAW_PR10_LAW_ID",
    "QSC_ACCEPTED_DRAW_SCHEMA_ID",
    "QSC_ACCEPTED_DRAW_SCHEMA_REVISION",
    "QSCAcceptedDrawExitError",
    "QSCAcceptedDrawExitEvaluator",
    "QSCAcceptedDrawExitLawAuthority",
    "QSCAcceptedDrawExitReplayError",
    "QSCAcceptedDrawExitRejection",
    "QSCAcceptedDrawExitRequest",
    "QSCAcceptedDrawFootprint",
    "QSCAcceptedTrayDrawLedger",
    "QSCCommonShaftDrawDisposition",
    "ValidatedQSCAcceptedDrawExit",
    "qsc_accepted_draw_flow_history_digest",
    "qsc_accepted_draw_shaft_history_digest",
    "validate_qsc_accepted_draw_exit_law_authority",
)
