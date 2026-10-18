"""Synchronize the six finite-speed K schedulers at one physical endpoint.

The common shaft does not make the six stochastic packet schedules share the
same first event.  A coupled tray step must therefore stop at the earliest of
the six candidate endpoints and advance every schedule to that exact time.

This module is deliberately a small engineering coordinator.  The underlying
K scheduler remains the authority for event times, target layers, optical
clocks, and planner replay.  This coordinator neither clips accepted flow nor
prorates feed: if an event shortens the requested interval, their owners must
reissue interval data for the shorter endpoint.
"""

from __future__ import annotations

import enum
import hashlib
import math
import struct
from dataclasses import dataclass
from typing import ClassVar

from . import engineering_physical_k_event_scheduler as k_scheduler
from . import qsc_mechanical_runtime_authority as mech_runtime


QSC_SIX_TRAY_K2_ENDPOINT_SCHEMA_ID = "dtdc-core2-qsc-six-tray-k2-event-endpoint-v2"
QSC_SIX_TRAY_K2_ENDPOINT_SCHEMA_REVISION = 2


class QSCSixTrayK2EventEndpointError(ValueError):
    """The six schedules could not be advanced on one common interval."""


class QSCSixTrayK2EventEndpointPlanningError(QSCSixTrayK2EventEndpointError):
    """A pure planning attempt failed and retained the exact input states."""

    exact_input_state_rollback: ClassVar[bool] = True

    def __init__(
        self,
        message: str,
        rollback_states: tuple[k_scheduler.EngineeringKScheduleState, ...],
    ) -> None:
        super().__init__(message)
        self.rollback_states = rollback_states


class QSCFlowFeedEndpointDisposition(enum.Enum):
    REQUESTED_HORIZON_REACHED = "REQUESTED_HORIZON_REACHED"
    SHORTENED_REISSUE_FLOW_AND_FEED = "SHORTENED_REISSUE_FLOW_AND_FEED"


def _require_finite_time(name: str, value: float) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise QSCSixTrayK2EventEndpointError(f"{name} must be a finite binary64")
    if struct.pack(">d", value) == b"\x80\x00\x00\x00\x00\x00\x00\x00":
        raise QSCSixTrayK2EventEndpointError(f"{name} rejects negative zero")


def _runtime_authorities(
    runtime_authority: mech_runtime.QSCMechanicalRuntimeAuthority,
) -> tuple[k_scheduler.EngineeringKEventAuthority, ...]:
    if type(runtime_authority) is not mech_runtime.QSCMechanicalRuntimeAuthority:
        raise TypeError("endpoint planning requires an exact QSC mechanical runtime")
    mech_runtime.validate_qsc_mechanical_runtime_authority(runtime_authority)
    authorities = runtime_authority.k_authorities
    if type(authorities) is not tuple or len(authorities) != 6:
        raise QSCSixTrayK2EventEndpointError("mechanical runtime must provide six K authorities")
    return authorities


def _validate_inputs(
    states: tuple[k_scheduler.EngineeringKScheduleState, ...],
    *,
    authorities: tuple[k_scheduler.EngineeringKEventAuthority, ...],
    requested_horizon_time_s: float,
) -> float:
    if type(states) is not tuple or len(states) != 6:
        raise TypeError("endpoint planning requires one exact six-state tuple")
    if any(type(state) is not k_scheduler.EngineeringKScheduleState for state in states):
        raise TypeError("endpoint planning received a foreign schedule-state type")
    expected_order = mech_runtime.QSC_MECH01_SOLID_TRAY_ORDER
    if tuple(state.physical_tray_id for state in states) != expected_order:
        raise QSCSixTrayK2EventEndpointError("schedule states are not in physical solid order")
    if tuple(authority.mesh.physical_tray_id for authority in authorities) != expected_order:
        raise QSCSixTrayK2EventEndpointError("K authorities are not in physical solid order")
    start_time_s = states[0].time_s
    if any(state.time_s != start_time_s for state in states):
        raise QSCSixTrayK2EventEndpointError("six schedules do not share one start time")
    _require_finite_time("requested endpoint horizon", requested_horizon_time_s)
    if requested_horizon_time_s <= start_time_s:
        raise QSCSixTrayK2EventEndpointError("endpoint horizon must advance time")
    for state, authority in zip(states, authorities, strict=True):
        if state.authority_definition_digest != authority.definition_digest:
            raise QSCSixTrayK2EventEndpointError(
                f"{state.physical_tray_id} schedule binds another K authority"
            )
    return start_time_s


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCTrayK2Endpoint:
    physical_tray_id: str
    prior_state: k_scheduler.EngineeringKScheduleState
    preview_step: k_scheduler.EngineeringKScheduleStep
    synchronized_step: k_scheduler.EngineeringKScheduleStep
    planner: k_scheduler.EngineeringPhysicalKRekeyPlanner

    def __post_init__(self) -> None:
        if self.physical_tray_id != self.prior_state.physical_tray_id:
            raise QSCSixTrayK2EventEndpointError("tray endpoint names another prior tray")
        if self.preview_step.prior is not self.prior_state:
            raise QSCSixTrayK2EventEndpointError("tray preview does not start from its prior")
        if self.synchronized_step.prior is not self.prior_state:
            raise QSCSixTrayK2EventEndpointError("tray synchronized step has another prior")
        if self.planner.step is not self.synchronized_step:
            raise QSCSixTrayK2EventEndpointError(
                "tray planner does not carry its synchronized step"
            )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCSixTrayK2EventEndpointPlan:
    common_start_time_s: float
    requested_horizon_time_s: float
    common_endpoint_time_s: float
    endpoint_disposition: QSCFlowFeedEndpointDisposition
    tray_endpoints: tuple[QSCTrayK2Endpoint, ...]

    host_atomic_persistence_implemented: ClassVar[bool] = False
    durable_evidence_artifact: ClassVar[bool] = False
    mech01_complete: ClassVar[bool] = False
    qsc10_complete: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    industrially_qualified: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_finite_time("common start time", self.common_start_time_s)
        _require_finite_time("requested horizon", self.requested_horizon_time_s)
        _require_finite_time("common endpoint", self.common_endpoint_time_s)
        if (
            not self.common_start_time_s
            < self.common_endpoint_time_s
            <= self.requested_horizon_time_s
        ):
            raise QSCSixTrayK2EventEndpointError("common endpoint lies outside its interval")
        if type(self.tray_endpoints) is not tuple or len(self.tray_endpoints) != 6:
            raise QSCSixTrayK2EventEndpointError("endpoint plan must contain six tray rows")
        if tuple(row.physical_tray_id for row in self.tray_endpoints) != (
            mech_runtime.QSC_MECH01_SOLID_TRAY_ORDER
        ):
            raise QSCSixTrayK2EventEndpointError("endpoint rows are not in physical solid order")
        if any(row.prior_state.time_s != self.common_start_time_s for row in self.tray_endpoints):
            raise QSCSixTrayK2EventEndpointError("endpoint rows do not share the common start")
        if any(
            row.synchronized_step.candidate.time_s != self.common_endpoint_time_s
            for row in self.tray_endpoints
        ):
            raise QSCSixTrayK2EventEndpointError("a synchronized schedule missed the endpoint")
        if self.common_endpoint_time_s != min(
            row.preview_step.candidate.time_s for row in self.tray_endpoints
        ):
            raise QSCSixTrayK2EventEndpointError("common endpoint is not the first tray endpoint")
        shortened = self.common_endpoint_time_s < self.requested_horizon_time_s
        expected = (
            QSCFlowFeedEndpointDisposition.SHORTENED_REISSUE_FLOW_AND_FEED
            if shortened
            else QSCFlowFeedEndpointDisposition.REQUESTED_HORIZON_REACHED
        )
        if self.endpoint_disposition is not expected:
            raise QSCSixTrayK2EventEndpointError("flow/feed endpoint disposition is inconsistent")

    @property
    def flow_inputs_reissue_required(self) -> bool:
        return self.endpoint_disposition is (
            QSCFlowFeedEndpointDisposition.SHORTENED_REISSUE_FLOW_AND_FEED
        )

    @property
    def feed_inputs_reissue_required(self) -> bool:
        return self.flow_inputs_reissue_required

    @property
    def endpoint_has_events(self) -> bool:
        return any(row.synchronized_step.events for row in self.tray_endpoints)

    @property
    def definition_digest(self) -> str:
        digest = hashlib.sha256()
        digest.update(QSC_SIX_TRAY_K2_ENDPOINT_SCHEMA_ID.encode("ascii"))
        for value in (
            self.common_start_time_s,
            self.requested_horizon_time_s,
            self.common_endpoint_time_s,
        ):
            digest.update(struct.pack(">d", value))
        digest.update(self.endpoint_disposition.value.encode("ascii"))
        for row in self.tray_endpoints:
            digest.update(row.physical_tray_id.encode("ascii"))
            digest.update(row.preview_step.definition_digest.encode("ascii"))
            digest.update(row.synchronized_step.definition_digest.encode("ascii"))
        return "sha256:" + digest.hexdigest()


def _build_plan(
    states: tuple[k_scheduler.EngineeringKScheduleState, ...],
    *,
    runtime_authority: mech_runtime.QSCMechanicalRuntimeAuthority,
    requested_horizon_time_s: float,
) -> QSCSixTrayK2EventEndpointPlan:
    authorities = _runtime_authorities(runtime_authority)
    start_time_s = _validate_inputs(
        states,
        authorities=authorities,
        requested_horizon_time_s=requested_horizon_time_s,
    )
    previews = tuple(
        k_scheduler.advance_engineering_k_schedule(
            state,
            authority=authority,
            horizon_time_s=requested_horizon_time_s,
        )
        for state, authority in zip(states, authorities, strict=True)
    )
    endpoint_time_s = min(step.candidate.time_s for step in previews)
    synchronized = tuple(
        k_scheduler.advance_engineering_k_schedule(
            state,
            authority=authority,
            horizon_time_s=endpoint_time_s,
        )
        for state, authority in zip(states, authorities, strict=True)
    )
    rows = tuple(
        QSCTrayK2Endpoint(
            physical_tray_id=state.physical_tray_id,
            prior_state=state,
            preview_step=preview,
            synchronized_step=step,
            planner=k_scheduler.EngineeringPhysicalKRekeyPlanner(
                authority=authority,
                step=step,
            ),
        )
        for state, authority, preview, step in zip(
            states,
            authorities,
            previews,
            synchronized,
            strict=True,
        )
    )
    return QSCSixTrayK2EventEndpointPlan(
        common_start_time_s=start_time_s,
        requested_horizon_time_s=requested_horizon_time_s,
        common_endpoint_time_s=endpoint_time_s,
        endpoint_disposition=(
            QSCFlowFeedEndpointDisposition.SHORTENED_REISSUE_FLOW_AND_FEED
            if endpoint_time_s < requested_horizon_time_s
            else QSCFlowFeedEndpointDisposition.REQUESTED_HORIZON_REACHED
        ),
        tray_endpoints=rows,
    )


def plan_qsc_six_tray_k2_event_endpoint(
    states: tuple[k_scheduler.EngineeringKScheduleState, ...],
    *,
    runtime_authority: mech_runtime.QSCMechanicalRuntimeAuthority,
    requested_horizon_time_s: float,
) -> QSCSixTrayK2EventEndpointPlan:
    """Return a pure common-endpoint plan or an exact-input rollback refusal."""

    rollback = states if type(states) is tuple else ()
    try:
        return _build_plan(
            states,
            runtime_authority=runtime_authority,
            requested_horizon_time_s=requested_horizon_time_s,
        )
    except Exception as error:
        if isinstance(error, QSCSixTrayK2EventEndpointPlanningError):
            raise
        raise QSCSixTrayK2EventEndpointPlanningError(
            f"six-tray K endpoint planning refused: {error}",
            rollback,
        ) from error


def validate_qsc_six_tray_k2_event_endpoint_plan(
    plan: QSCSixTrayK2EventEndpointPlan,
    *,
    states: tuple[k_scheduler.EngineeringKScheduleState, ...],
    runtime_authority: mech_runtime.QSCMechanicalRuntimeAuthority,
    requested_horizon_time_s: float,
) -> str:
    """Deterministically replay a plan and return its reproducibility digest."""

    if type(plan) is not QSCSixTrayK2EventEndpointPlan:
        raise TypeError("endpoint replay requires an exact endpoint plan")
    replayed = _build_plan(
        states,
        runtime_authority=runtime_authority,
        requested_horizon_time_s=requested_horizon_time_s,
    )
    if replayed != plan:
        raise QSCSixTrayK2EventEndpointError("endpoint plan differs from deterministic replay")
    return plan.definition_digest


__all__ = [
    "QSCFlowFeedEndpointDisposition",
    "QSCSixTrayK2EventEndpointError",
    "QSCSixTrayK2EventEndpointPlan",
    "QSCSixTrayK2EventEndpointPlanningError",
    "QSCTrayK2Endpoint",
    "plan_qsc_six_tray_k2_event_endpoint",
    "validate_qsc_six_tray_k2_event_endpoint_plan",
]
