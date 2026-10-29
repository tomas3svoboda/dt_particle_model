"""Evaluate one prepared six-tray K=2 mechanism-integration interval.

This coordinator composes already bounded engineering seams at one common
finite-speed endpoint: shaft-work allocation and meal deposition, the
single-physical-wall SP1 -> MN2 -> MN1 through-bed chain, the zero-contact
PD3 -> PD2 -> PD1 bypass chain, accepted draw routing, and reconciliation of
the routed populations with their engineering K clocks.

The accepted dry-matter histories remain caller inputs.  They may differ by
tray, so the result records the resulting routing accumulation explicitly.
This is neither a nominal nor a steady-flow tower model.  It adds no commit or
durable host transaction and advances no QSC-10, physical-qualification,
plant-prediction, calibration, production, or industrial-feasibility claim.
"""

from __future__ import annotations

import hashlib
import math
import os
from dataclasses import dataclass, field, replace
from typing import ClassVar

from . import accepted_k_cell_tower_router as tower_router
from . import cell_engineering_feasibility as ef
from . import cell_native_caloric as native_caloric
from . import engineering_k_tower_route_reconciler as k_reconciler
from . import engineering_physical_k_event_scheduler as k_scheduler
from . import engineering_dual_variant_packet_adapter as dual_variant
from . import high_loading_packet_adapter as high_loading
from . import k_cell_tray_host as tray_host
from . import pd_k_dome_bypass_tray_integration as pd
from . import qsc_accepted_draw_exit_coupling as accepted_draw
from . import qsc_meal_mechanical_work_deposition as meal_work
from . import qsc_mechanical_runtime_authority as mech_runtime
from . import qsc_pd_bypass_pressure_network as pd_network
from . import qsc_pd_zero_contact_network_coupling as pd_zero_contact
from . import qsc_shaft_work_energy_coupling as shaft_work
from . import qsc_six_tray_k2_event_endpoint as k_endpoint
from . import through_bed_k2_single_wall_tray_integration as through_bed
from .packet_weight_scaling_contract import derive_packet_totals

QSC_K2_MECHANISM_INTERVAL_SCHEMA_ID = "dtdc-core2-qsc-k2-mechanism-interval-v1"
SOLID_TRAY_ORDER = mech_runtime.QSC_MECH01_SOLID_TRAY_ORDER
THROUGH_BED_SOLID_ORDER = ("MN1", "MN2", "SP1")
PD_CAUSAL_ORDER = ("PD3", "PD2", "PD1")
MN1_TOP_PRESSURE_PA = 101_325.0
GAS_HANDOFF_AUTHORITY_ID = "qsc-k2-mechanism-same-tn-countercurrent-gas-v1"

#: EXPLORATION MODE (owner ruling of 2026-10-13, record
#: docs/GT_PS2_OWNER_RULING_EXPLORATION_MODE_2026-10-13.md): the breach-aware gate of the
#: interval's conservative-ledger check.  ``None`` on disk; ``study_lane/exploration_rams.py``
#: binds it when the inter-K or packet-ledger family is armed.  While it is None nothing changes.
EXPLORATION_LEDGER_GATE = None

#: THE TOP-FACE ENERGY BOOKING (knob, ON by default: owner ruling six of
#: 2026-09-29, recorded in ``docs/GT_PS2_OWNER_RULINGS_2026-09-29_CAMPAIGN_DAY3.md``
#: and its addendum
#: ``docs/GT_PS2_OWNER_RULINGS_2026-09-29_ENERGY_BOOKING_AND_COARSE_REGIME.md``,
#: on ``docs/GT_PS2_ENERGY_HANDOFF_DIAGNOSIS_2026-09-28.md`` section 3, item 3,
#: "an open booking item, the owner's call (Q-F2a class)"; built behind this knob
#: so that the ruling flips a default and nothing else).
#:
#: Criterion.  A through-bed tray's energy ledger books its top outlet at the
#: accepted upper layer's pressure p_L2 (``sp1_k2_law2_tray_integration.py``, the
#: W4 per-cell convention), while the tower hands the SAME stream (F, T, y bit
#: for bit) on at the top boundary pressure p_top.  On the native caloric arm the
#: enthalpy depends on pressure, so the step from p_L2 to p_top is an energy term
#: that no ledger books: -1,048 W on the MN1 row the diagnosis reproduced to
#: 3.9e-10 W, 1,731 W over the three through-bed trays on its S1 step 1.
#:
#: Value.  At every through-bed seam (SP1 -> MN2, MN2 -> MN1, MN1 -> PD3) the
#: handed point carries F [h(T, p_top, y) - h(T, p_L2, y)] (the exact frozen-law
#: difference, as Q-F2a books the tower-seam term), with p_L2 the accepted upper
#: step's layer pressure, which is also the pressure at which the tray ledger
#: quotes a composed layer's presented stream (the level-0b seam presents the
#: upper face mix on the aggregate root's state).  The term is checked by
#: RECONSTRUCTION: F h(T, p_L2, y) must equal the source tray ledger's own
#: ``top_outlet_energy_w`` within the ULP sum of the two operands and their
#: difference, else the booking refuses typed (the stream was not the one the
#: ledger booked).  The hand-off audit carries it in
#: ``booked_top_face_throttling_w`` (exact zero by default and off the knob) and
#: its energy residual becomes tray-ledger top outlet + top-face term - handed
#: destination + Q-F2a defect, so the chain from the tray ledger to the next
#: tray closes; ``tower_gas_enthalpy_telescope`` sums it over the tower.  On the
#: DECLARED caloric arm (pressure-free by model form) the term is exactly zero.
#: On a tray-local interval composed from N sub-steps the booking is composed
#: with the same duration weights as the outlet and the composite outlet
#: ledger's own energy bound enters the hand-off bound (the tray ledger belongs
#: to one sub-step, so each sub-step is reconstructed against its own ledger).
#:
#: Bracket.  Off (``DTDC_HANDOFF_TOPFACE_BOOKING=0``, the identity twin of every
#: record): the pre-ruling audit byte for byte (the term unbooked, the residual
#: quoted at p_top).  On (the default): the term booked at all three seams.
#:
#: Rejecting outcomes, written first.  IDENTITY: a knob-off lane changing a
#: committed byte.  DECLARED ARM: a nonzero term on the declared caloric arm.
#: RECONSTRUCTION: F h(T, p_L2, y) differing from the tray ledger's top outlet
#: beyond the ULP sum (refused typed).  CLOSURE: the booked hand-off residual
#: outside its ULP bound on an otherwise passing interval, or the telescope
#: residual with the booking outside its bound on an unsplit interval.
HANDOFF_TOPFACE_BOOKING_ENVIRONMENT_VARIABLE = "DTDC_HANDOFF_TOPFACE_BOOKING"
HANDOFF_TOPFACE_BOOKING_ARMED = (
    os.environ.get(HANDOFF_TOPFACE_BOOKING_ENVIRONMENT_VARIABLE, "1") != "0"
)
HANDOFF_TOPFACE_BOOKING_LAW_ID = "handoff-top-face-throttling-booking-v1"


class QSCK2MechanismIntervalError(ValueError):
    """A prepared K=2 mechanism interval could not be evaluated."""


class QSCK2MechanismPhaseExhaustionError(QSCK2MechanismIntervalError):
    """A tray refused on Law-2 layer phase exhaustion; the typed signal rides.

    B1 (owner-signed design): message text is exactly what the plain
    error would have carried, so every string matcher is unchanged;
    ``signal`` is the sp1 ``LayerPhaseExhaustionSignal`` for callers
    that re-request the shortened interval.
    """

    def __init__(self, message: str, *, signal: object) -> None:
        super().__init__(message)
        self.signal = signal


def _require_identity(label: str, value: str) -> None:
    if type(value) is not str or not value or value != value.strip():
        raise QSCK2MechanismIntervalError(f"{label} must be a nonblank exact string")


def _require_finite(label: str, value: float, *, positive: bool = False) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise QSCK2MechanismIntervalError(f"{label} must be a finite exact binary64")
    if positive and value <= 0.0:
        raise QSCK2MechanismIntervalError(f"{label} must be strictly positive")


def _host_trays(
    pd_states: tuple[pd.AcceptedPDKTrayState, ...],
    through_states: tuple[through_bed.AcceptedThroughBedK2SingleWallTrayState, ...],
) -> tuple[tray_host.AcceptedKCellTrayState, ...]:
    return tuple(item.tray for item in (*pd_states, *through_states))


def _validate_solid_order(trays: tuple[tray_host.AcceptedKCellTrayState, ...]) -> None:
    if tuple(item.physical_tray_id for item in trays) != SOLID_TRAY_ORDER:
        raise QSCK2MechanismIntervalError(
            "solid trays must be ordered PD1, PD2, PD3, MN1, MN2, SP1"
        )
    t_n_s = trays[0].time_s
    if any(item.time_s != t_n_s for item in trays):
        raise QSCK2MechanismIntervalError("all six local priors must share one accepted time")
    for tray in trays:
        physical = next(
            item for item in tray.topology.trays if item.physical_tray_id == tray.physical_tray_id
        )
        if physical.vertical_layer_count != 2 or physical.rtd_stage_count != 8:
            raise QSCK2MechanismIntervalError("mechanism interval requires exact K=2 / RTD=8")
        if (
            tray.topology.digest != trays[0].topology.digest
            or tray.configuration_digest != trays[0].configuration_digest
            or tray.codec_registry_digest != trays[0].codec_registry_digest
            or tray.auditor_identity_digest != trays[0].auditor_identity_digest
            or tray.energy_datum_id != trays[0].energy_datum_id
            or tray.rtd_schedule_digest != trays[0].rtd_schedule_digest
            or tray.pressure_pa != trays[0].pressure_pa
        ):
            raise QSCK2MechanismIntervalError(
                "six-tray host topology, codec, datum, RTD, or pressure context differs"
            )
    if trays[0].pressure_pa != MN1_TOP_PRESSURE_PA:
        raise QSCK2MechanismIntervalError(
            "the zero-contact PD mechanism requires fixed 101325 Pa packet pressure"
        )


def _validate_schedule_pairing(
    trays: tuple[tray_host.AcceptedKCellTrayState, ...],
    schedules: tuple[k_scheduler.EngineeringKScheduleState, ...],
) -> None:
    if type(schedules) is not tuple or len(schedules) != 6:
        raise QSCK2MechanismIntervalError("six exact engineering K states are required")
    if tuple(item.physical_tray_id for item in schedules) != SOLID_TRAY_ORDER:
        raise QSCK2MechanismIntervalError("K states differ from physical solid order")
    for tray, schedule in zip(trays, schedules, strict=True):
        packets = {item.packet_id: item for item in tray.packets}
        clocks = {item.packet_id: item for item in schedule.packet_clocks}
        if schedule.time_s != tray.time_s or set(packets) != set(clocks):
            raise QSCK2MechanismIntervalError("tray population and K-clock state differ")
        if any(
            packets[packet_id].owner_key.vertical_layer_id != clock.vertical_layer_id
            for packet_id, clock in clocks.items()
        ):
            raise QSCK2MechanismIntervalError("packet K owner and clock K owner differ")


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCK2MechanismIntervalRequest:
    """Prepared inputs for one common finite-speed K=2 interval."""

    attempt_id: str
    pd_states: tuple[pd.AcceptedPDKTrayState, ...]
    pd_models: tuple[pd.PDKTrayEngineeringModel, ...]
    through_bed_states: tuple[through_bed.AcceptedThroughBedK2SingleWallTrayState, ...]
    through_bed_layer_models: tuple[
        tuple[
            through_bed.ThroughBedKLayerEngineeringModel,
            through_bed.ThroughBedKLayerEngineeringModel,
        ],
        ...,
    ]
    endpoint_plan: k_endpoint.QSCSixTrayK2EventEndpointPlan
    sp1_bottom_gas_boundary: ef.EngineeringGasBoundary
    flow_histories: tuple[tuple[tray_host.AcceptedDryMatterFlowSegment, ...], ...]
    closing_layer_cursors: tuple[int, ...]
    feed_packets: tuple[tower_router.AcceptedFeedPacket, ...]
    shaft_power_w: float
    pressure_face_relative_limit: float = 1.0e-10
    pressure_root_max_evaluations: int = 20
    #: Optional ruled K-share closing-exchange policy, forwarded verbatim to
    #: the accepted-draw request (same optional-carrier pattern the
    #: 2026-08-29 v3 approval ruling added to the router and draw coupling).
    #: ``None`` preserves the sealed refuse-on-unrepresentable behavior.
    k_share_closing_exchange_policy: tower_router.KShareClosingExchangePolicy | None = None
    #: WS-0 (owner-ruled 2026-09-04, design packet
    #: ``GT_PS2_PD_PHASE_CHANGE_PROGRAM_DESIGN_2026-09-03.md`` S:2.1/S:2.0):
    #: a DECLARED, ADDITIVE knob overriding the vapour mole-fraction
    #: composition (hexane, water) the PD3 bypass inlet is evaluated
    #: against.  ``None`` (the default) preserves today's behavior EXACTLY
    #: -- the composition is read from the physics-chain-derived MN1 outlet
    #: (pure steam at the sparge boundary, undiluted by construction while
    #: ``flash_or_condensation_implemented`` stays ``False``).  When set, the
    #: pair overrides ONLY the mole fractions PD3's bypass inlet is built
    #: from; molar flow, pressure and temperature stay physics-derived from
    #: the MN1 outlet.  Applied ONCE, to the single ``mn1_outlet`` value
    #: every downstream reader of "what MN1 hands off to PD3" consumes --
    #: the actual PD3 bypass physics, the MN1->PD3 gas-handoff conservation
    #: ledger, and the zero-contact network's own independent cross-check --
    #: so all three see the identical (possibly overridden) boundary and
    #: none of them can disagree with another about it.
    pd3_inlet_vapour_composition_override: tuple[float, float] | None = None

    mechanism_integration_fixture: ClassVar[bool] = True
    nominal_flow_model: ClassVar[bool] = False
    steady_throughput_demonstrated: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("mechanism interval attempt", self.attempt_id)
        if (
            type(self.pd_states) is not tuple
            or len(self.pd_states) != 3
            or any(type(item) is not pd.AcceptedPDKTrayState for item in self.pd_states)
            or type(self.pd_models) is not tuple
            or len(self.pd_models) != 3
            or any(type(item) is not pd.PDKTrayEngineeringModel for item in self.pd_models)
        ):
            raise TypeError("mechanism interval requires three exact PD states and models")
        if (
            type(self.through_bed_states) is not tuple
            or len(self.through_bed_states) != 3
            or any(
                type(item) is not through_bed.AcceptedThroughBedK2SingleWallTrayState
                for item in self.through_bed_states
            )
            or type(self.through_bed_layer_models) is not tuple
            or len(self.through_bed_layer_models) != 3
            or any(
                type(item) is not tuple or len(item) != 2 for item in self.through_bed_layer_models
            )
        ):
            raise TypeError("mechanism interval requires three exact one-wall state/model pairs")
        trays = _host_trays(self.pd_states, self.through_bed_states)
        _validate_solid_order(trays)
        if type(self.endpoint_plan) is not k_endpoint.QSCSixTrayK2EventEndpointPlan:
            raise TypeError("mechanism interval requires an exact synchronized endpoint plan")
        if (
            self.endpoint_plan.common_start_time_s != trays[0].time_s
            or tuple(row.physical_tray_id for row in self.endpoint_plan.tray_endpoints)
            != SOLID_TRAY_ORDER
        ):
            raise QSCK2MechanismIntervalError("endpoint plan differs from the six local priors")
        _validate_schedule_pairing(
            trays,
            tuple(row.prior_state for row in self.endpoint_plan.tray_endpoints),
        )
        if type(self.sp1_bottom_gas_boundary) is not ef.EngineeringGasBoundary:
            raise TypeError("SP1 lower gas boundary has a foreign exact type")
        if type(self.flow_histories) is not tuple or len(self.flow_histories) != 6:
            raise TypeError("mechanism interval requires six accepted flow histories")
        endpoint = self.endpoint_plan.common_endpoint_time_s
        for history in self.flow_histories:
            if (
                type(history) is not tuple
                or not history
                or any(type(item) is not tray_host.AcceptedDryMatterFlowSegment for item in history)
                or history[0].start_time_s != trays[0].time_s
                or history[-1].end_time_s != endpoint
            ):
                raise QSCK2MechanismIntervalError(
                    "every accepted flow history must span the prepared endpoint"
                )
        if (
            type(self.closing_layer_cursors) is not tuple
            or len(self.closing_layer_cursors) != 6
            or any(type(item) is not int or item <= 0 for item in self.closing_layer_cursors)
        ):
            raise QSCK2MechanismIntervalError("six positive closing-layer cursors are required")
        if type(self.feed_packets) is not tuple or any(
            type(item) is not tower_router.AcceptedFeedPacket for item in self.feed_packets
        ):
            raise TypeError("mechanism interval feed contains a foreign packet")
        _require_finite("shaft power", self.shaft_power_w, positive=True)
        _require_finite(
            "pressure-face relative limit",
            self.pressure_face_relative_limit,
            positive=True,
        )
        if (
            type(self.pressure_root_max_evaluations) is not int
            or self.pressure_root_max_evaluations <= 0
        ):
            raise QSCK2MechanismIntervalError(
                "pressure-root evaluation limit must be a positive exact integer"
            )
        if self.k_share_closing_exchange_policy is not None:
            if (
                type(self.k_share_closing_exchange_policy)
                is not tower_router.KShareClosingExchangePolicy
            ):
                raise TypeError("closing-exchange policy has a foreign exact type")
            self.k_share_closing_exchange_policy.__post_init__()
        if self.pd3_inlet_vapour_composition_override is not None:
            override = self.pd3_inlet_vapour_composition_override
            if (
                type(override) is not tuple
                or len(override) != 2
                or any(type(item) is not float for item in override)
            ):
                raise TypeError(
                    "PD3 inlet vapour composition override must be an exact "
                    "(hexane_mole_fraction, water_mole_fraction) float pair"
                )
            hexane_fraction, water_fraction = override
            _require_finite("PD3 inlet override hexane mole fraction", hexane_fraction)
            _require_finite("PD3 inlet override water mole fraction", water_fraction)
            if not (0.0 <= hexane_fraction <= 1.0) or not (0.0 <= water_fraction <= 1.0):
                raise QSCK2MechanismIntervalError(
                    "PD3 inlet vapour composition override mole fractions must lie in [0, 1]"
                )
            if math.fsum((hexane_fraction, water_fraction)) != 1.0:
                raise QSCK2MechanismIntervalError(
                    "PD3 inlet vapour composition override must be exactly binary "
                    "(hexane + water mole fractions sum to 1)"
                )

    @property
    def prior_host_trays(self) -> tuple[tray_host.AcceptedKCellTrayState, ...]:
        return _host_trays(self.pd_states, self.through_bed_states)


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMechanismPressureSolveEvidence:
    evaluations: int
    sp1_mn2_face_pressure_pa: float
    mn2_mn1_face_pressure_pa: float
    sp1_mn2_residual_pa: float
    mn2_mn1_residual_pa: float
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return self.maximum_relative_residual <= self.relative_limit


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMechanismGasHandoffAudit:
    source_tray_id: str
    destination_tray_id: str
    pressure_residual_pa: float
    pressure_relative_residual: float
    hexane_molar_residual_mol_s: float
    water_molar_residual_mol_s: float
    common_datum_energy_residual_w: float
    hexane_molar_roundoff_bound_mol_s: float
    water_molar_roundoff_bound_mol_s: float
    common_datum_energy_roundoff_bound_w: float
    pressure_relative_limit: float
    #: Q-F2a (owner ruling 2026-09-02, batch commit 20f5028).  Exact
    #: positive zero on the DECLARED arm, whose stream enthalpy is
    #: pressure-free by model form, so every pre-flip audit reconstructs
    #: bit-identically.  On the NATIVE arm the tower's pressure-consistency
    #: solve leaves a RELATIVE pressure residual (this ledger's own
    #: ``pressure_relative_limit``), and a pressure-dependent enthalpy
    #: turns that residual into an energy-rate difference.  It is BOOKED
    #: here - the same treatment CELL-02d W5 gives the internal tray face -
    #: rather than absorbed into a widened bound.  Booked ONLY when the
    #: hand-off is otherwise verbatim (bit-identical temperature and
    #: composition); any other difference still bites the energy row.
    booked_pressure_throttling_defect_w: float = 0.0
    #: TOP-FACE BOOKING (``DTDC_HANDOFF_TOPFACE_BOOKING``, on by default): the
    #: source tray's step from its upper layer pressure to the handed pressure,
    #: F [h(T, p_top, y) - h(T, p_L2, y)].  Exact positive zero off the knob, on
    #: the declared arm and on a source that is not a through-bed tray.
    booked_top_face_throttling_w: float = 0.0
    #: TOP-FACE BOOKING: the audit's source-side energy (the tray ledger's top
    #: outlet, or its duration-weighted composite), the destination's handed
    #: energy and the composition count; None off the knob.
    source_tray_top_outlet_energy_w: float | None = None
    destination_energy_rate_w: float | None = None
    source_sub_steps: int | None = None

    physically_qualifying: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.pressure_relative_residual <= self.pressure_relative_limit
            and abs(self.hexane_molar_residual_mol_s) <= self.hexane_molar_roundoff_bound_mol_s
            and abs(self.water_molar_residual_mol_s) <= self.water_molar_roundoff_bound_mol_s
            and abs(self.common_datum_energy_residual_w)
            <= self.common_datum_energy_roundoff_bound_w
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMechanismSolidsAccumulationEvidence:
    feed_dry_matter_kg: float
    accepted_draw_dry_matter_kg: tuple[float, ...]
    tray_routing_accumulation_dry_matter_kg: tuple[float, ...]

    nominal_flow_model: ClassVar[bool] = False
    steady_throughput_demonstrated: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    @property
    def has_unequal_per_tray_draw(self) -> bool:
        return len(set(self.accepted_draw_dry_matter_kg)) > 1

    @property
    def has_routing_accumulation(self) -> bool:
        return any(value != 0.0 for value in self.tray_routing_accumulation_dry_matter_kg)


_PREPARED_CANDIDATE_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class PreparedQSCK2MechanismCandidate:
    pd_states: tuple[pd.AcceptedPDKTrayState, ...]
    through_bed_states: tuple[through_bed.AcceptedThroughBedK2SingleWallTrayState, ...]
    k_schedule_states: tuple[k_scheduler.EngineeringKScheduleState, ...]
    _seal: object = field(repr=False, compare=False)

    prepared_not_committed: ClassVar[bool] = True
    durable_atomic_persistence_implemented: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _PREPARED_CANDIDATE_SEAL:
            raise TypeError("prepared mechanism candidates are evaluator-issued only")
        trays = _host_trays(self.pd_states, self.through_bed_states)
        _validate_solid_order(trays)
        _validate_schedule_pairing(trays, self.k_schedule_states)

    @property
    def time_s(self) -> float:
        return self.pd_states[0].tray.time_s

    @property
    def host_trays(self) -> tuple[tray_host.AcceptedKCellTrayState, ...]:
        return _host_trays(self.pd_states, self.through_bed_states)


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCK2MechanismIntervalRejection:
    attempt_id: str
    reason: str
    rollback_pd_states: tuple[pd.AcceptedPDKTrayState, ...]
    rollback_through_bed_states: tuple[
        through_bed.AcceptedThroughBedK2SingleWallTrayState,
        ...,
    ]
    #: B1 (owner-signed): typed Law-2 phase-exhaustion transition data
    #: when that is the refusal; None otherwise (byte-compat default).
    phase_exhaustion: through_bed.legacy.LayerPhaseExhaustionSignal | None = None

    exact_input_rollback: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False


_VALIDATED_INTERVAL_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedQSCK2MechanismInterval:
    request: QSCK2MechanismIntervalRequest
    endpoint_plan: k_endpoint.QSCSixTrayK2EventEndpointPlan
    allocation_request: shaft_work.QSCMechanicalAllocationRequest
    mechanical_work_allocation: shaft_work.QSCMechanicalWorkAllocation
    deposition_bundles: tuple[meal_work.QSCMealMechanicalWorkDepositionBundle, ...]
    pressure_solve: QSCMechanismPressureSolveEvidence
    through_bed_transactions: tuple[
        through_bed.ValidatedThroughBedK2SingleWallTrayStep,
        ...,
    ]
    pd_transactions: tuple[pd.ValidatedPDKTrayStep, ...]
    gas_handoffs: tuple[QSCMechanismGasHandoffAudit, ...]
    pd_zero_contact_network: pd_zero_contact.QSCPDZeroContactNetworkCoupling
    accepted_draw_transaction: accepted_draw.ValidatedQSCAcceptedDrawExit
    reconciliation_transaction: k_reconciler.ValidatedEngineeringKTowerRouteReconciliation
    solids_accumulation: QSCMechanismSolidsAccumulationEvidence
    prepared_candidate: PreparedQSCK2MechanismCandidate
    _seal: object = field(repr=False, compare=False)

    mechanism_integration_fixture: ClassVar[bool] = True
    prepared_candidate_only: ClassVar[bool] = True
    nominal_flow_model: ClassVar[bool] = False
    steady_throughput_demonstrated: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    industrially_qualified: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _VALIDATED_INTERVAL_SEAL:
            raise TypeError("validated mechanism intervals are evaluator-issued only")

    @property
    def passed(self) -> bool:
        return (
            self.mechanical_work_allocation.ledger.authoritative_conservation
            and self.pressure_solve.passed
            and all(
                item.inter_k_face.passed
                and item.tray_pressure.passed
                and item.whole_tray_energy.passed
                and item.coupled_ledger.passed
                for item in self.through_bed_transactions
            )
            and all(
                item.inter_k_face.passed
                and item.pressure_ledger.passed
                and item.plenum_ledger.passed
                and item.wall_ledger.passed
                and item.system_ledger.passed
                for item in self.pd_transactions
            )
            and all(item.passed for item in self.gas_handoffs)
            and self.pd_zero_contact_network.passed
            and self.accepted_draw_transaction.routing_transaction.whole_router_ledger.passed
            and self.reconciliation_transaction.ledger.passed
        )


@dataclass(frozen=True, slots=True)
class _GasPoint:
    pressure_pa: float
    molar_flow_mol_s: float
    hexane_mole_fraction: float
    water_mole_fraction: float
    temperature_k: float
    common_datum_energy_rate_w: float
    #: Q-F2a (2026-09-02): whether this point's enthalpy came from the
    #: native frozen-law arm, whose enthalpy is a function of PRESSURE.
    #: The declared arm's enthalpy is pressure-free by model form.
    native_arm: bool = False
    #: TOP-FACE BOOKING (``DTDC_HANDOFF_TOPFACE_BOOKING``): None unless booked.
    top_face: _TopFaceBooking | None = None


@dataclass(frozen=True, slots=True)
class _TopFaceBooking:
    """The top-face term of one handed through-bed outlet (see the knob's row)."""

    #: F [h(T, p_top, y) - h(T, p_L2, y)] (W); exact 0.0 on the declared arm.
    term_w: float
    #: The source tray ledger's own top outlet energy (W), or the
    #: duration-weighted composite of its sub-steps' ledgers.
    tray_ledger_top_outlet_energy_w: float
    #: The accepted upper layer pressure p_L2 (Pa) of the (last) sub-step.
    upper_layer_pressure_pa: float
    #: Largest |F h(T, p_L2, y) - top_outlet| over the reconstructed ledgers,
    #: and the ULP bound it was checked against (W).
    reconstruction_residual_w: float
    reconstruction_bound_w: float
    #: The ULP resolution of the booked terms, and, on a composite, the
    #: composite outlet ledger's energy bound per second (W).
    resolution_bound_w: float
    sub_steps: int = 1


def _common_enthalpy_j_mol(
    properties: object,
    codec: high_loading.HighLoadingV1PacketAdapter,
    *,
    temperature_k: float,
    pressure_pa: float,
    hexane_mole_fraction: float,
    water_mole_fraction: float,
) -> float:
    """One tower-seam stream enthalpy, on whichever caloric arm the lane runs.

    Q-F2a (owner ruling 2026-09-02, batch commit ``20f5028``): the declared
    arm below is the pre-flip arithmetic VERBATIM and still ignores the
    pressure (its enthalpy is pressure-free by model form); the native arm
    evaluates the frozen binary mixture law at the stream's own pressure,
    exactly as ``sp1_k2_law2_tray_integration._stream_enthalpy_j_mol`` does
    at the tray seam, so the two seams cannot disagree about a stream.
    """

    codec.component_datum_adapter.require_gas_properties(properties)
    if type(properties) is native_caloric.NativeGasProperties:
        if water_mole_fraction != 1.0 - hexane_mole_fraction:
            raise QSCK2MechanismIntervalError(
                "the native binary stream requires the exact complementary water fraction"
            )
        return native_caloric.mixture_enthalpy(
            temperature_k, pressure_pa, hexane_mole_fraction
        ).value_j_mol
    hexane_shift, water_shift = codec.component_datum_adapter.common_datum_species_shifts_j_mol
    return math.fsum(
        (
            properties.molar_heat_capacity_j_mol_k
            * (temperature_k - properties.energy_datum_temperature_k),
            hexane_mole_fraction * (properties.hexane_latent_heat_j_mol + hexane_shift),
            water_mole_fraction * (properties.water_latent_heat_j_mol + water_shift),
        )
    )


def _gas_point(
    *,
    pressure_pa: float,
    molar_flow_mol_s: float,
    hexane_mole_fraction: float,
    water_mole_fraction: float,
    temperature_k: float,
    properties: object,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> _GasPoint:
    enthalpy = _common_enthalpy_j_mol(
        properties,
        codec,
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        hexane_mole_fraction=hexane_mole_fraction,
        water_mole_fraction=water_mole_fraction,
    )
    return _GasPoint(
        pressure_pa=pressure_pa,
        molar_flow_mol_s=molar_flow_mol_s,
        hexane_mole_fraction=hexane_mole_fraction,
        water_mole_fraction=water_mole_fraction,
        temperature_k=temperature_k,
        common_datum_energy_rate_w=molar_flow_mol_s * enthalpy,
        native_arm=type(properties) is native_caloric.NativeGasProperties,
    )


def _through_outlet(
    transaction: through_bed.ValidatedThroughBedK2SingleWallTrayStep,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> _GasPoint:
    step = transaction.accepted_fast_steps[1]
    state = step.state
    return _gas_point(
        pressure_pa=transaction.request.top_boundary_pressure_pa,
        molar_flow_mol_s=step.evaluation.outlet_molar_flow_mol_s,
        hexane_mole_fraction=state.gas_hexane_mole_fraction,
        water_mole_fraction=state.gas_water_mole_fraction,
        temperature_k=state.gas_temperature_k,
        properties=transaction.request.layer_models[1].properties,
        codec=codec,
    )


def _tray_ledger_top_outlet_energy_w(
    transaction: through_bed.ValidatedThroughBedK2SingleWallTrayStep,
) -> float:
    """The tray's THERMAL ledger top outlet, whichever wrapper carries it."""

    ledger = transaction.whole_tray_energy
    inner = getattr(ledger, "thermal_ledger", None)
    return (ledger if inner is None else inner).top_outlet_energy_w


def _book_top_face(
    point: _GasPoint,
    transaction: through_bed.ValidatedThroughBedK2SingleWallTrayStep,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> _GasPoint:
    """Attach the top-face term to a handed through-bed outlet (knob armed)."""

    upper_pressure_pa = transaction.accepted_fast_steps[1].state.layer_pressure_pa
    ledger_top_w = _tray_ledger_top_outlet_energy_w(transaction)
    properties = transaction.request.layer_models[1].properties
    at_upper_w = point.molar_flow_mol_s * _common_enthalpy_j_mol(
        properties,
        codec,
        temperature_k=point.temperature_k,
        pressure_pa=upper_pressure_pa,
        hexane_mole_fraction=point.hexane_mole_fraction,
        water_mole_fraction=point.water_mole_fraction,
    )
    reconstruction = math.fsum((at_upper_w, -ledger_top_w))
    reconstruction_bound = math.fsum(
        math.ulp(value) for value in (at_upper_w, ledger_top_w, reconstruction)
    )
    if not abs(reconstruction) <= reconstruction_bound:
        raise QSCK2MechanismIntervalError(
            f"{transaction.request.prior.tray.physical_tray_id} top-face booking does not "
            "reconstruct the tray ledger's top outlet (F h(T, p_L2, y) "
            f"{at_upper_w!r} W against "
            f"{ledger_top_w!r} W, residual {reconstruction!r} beyond {reconstruction_bound!r})"
        )
    if not point.native_arm and point.common_datum_energy_rate_w != at_upper_w:
        raise QSCK2MechanismIntervalError(
            "top-face booking: the declared caloric arm's stream enthalpy moved with pressure"
        )
    term = math.fsum((point.common_datum_energy_rate_w, -at_upper_w)) if point.native_arm else 0.0
    return replace(
        point,
        top_face=_TopFaceBooking(
            term_w=term,
            tray_ledger_top_outlet_energy_w=ledger_top_w,
            upper_layer_pressure_pa=upper_pressure_pa,
            reconstruction_residual_w=abs(reconstruction),
            reconstruction_bound_w=reconstruction_bound,
            resolution_bound_w=math.fsum(
                math.ulp(value) for value in (term, ledger_top_w, point.common_datum_energy_rate_w)
            ),
        ),
    )


def _handoff_source_outlet(
    transaction: through_bed.ValidatedThroughBedK2SingleWallTrayStep,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> _GasPoint:
    """The through-bed outlet a hand-off starts from.

    ``_through_outlet`` is resolved at call time, so a study-lane presentation
    wrapper is honoured; the top-face term is attached when
    ``DTDC_HANDOFF_TOPFACE_BOOKING`` is armed, and the very same point is
    returned otherwise.
    """

    point = _through_outlet(transaction, codec=codec)
    if not HANDOFF_TOPFACE_BOOKING_ARMED:
        return point
    return _book_top_face(point, transaction, codec=codec)


def compose_top_face_booking(
    points: tuple[_GasPoint, ...],
    durations: tuple[float, ...],
    *,
    composite: _GasPoint,
    composite_energy_bound_j: float,
) -> _GasPoint:
    """The duration-weighted top-face booking of a composed interval outlet.

    Used by the tray-local lane when an interval's outlet is composed from N
    sub-steps (``tray_local_events.compose_interval_outlet``).  A composite of
    points that carry no booking is returned unchanged.  The composite outlet
    ledger's energy bound, per second, enters the resolution bound: the
    composite's energy differs from the weighted sub-step energies by the
    composition's own bounded residual.
    """

    bookings = [point.top_face for point in points]
    if all(item is None for item in bookings):
        return composite
    if any(item is None for item in bookings):
        raise QSCK2MechanismIntervalError(
            "top-face booking: an interval's sub-step outlets are booked unevenly"
        )
    total = math.fsum(durations)
    term = math.fsum(d * item.term_w for d, item in zip(durations, bookings, strict=True)) / total
    ledger_top = (
        math.fsum(
            d * item.tray_ledger_top_outlet_energy_w
            for d, item in zip(durations, bookings, strict=True)
        )
        / total
    )
    return replace(
        composite,
        top_face=_TopFaceBooking(
            term_w=term,
            tray_ledger_top_outlet_energy_w=ledger_top,
            upper_layer_pressure_pa=bookings[-1].upper_layer_pressure_pa,
            reconstruction_residual_w=max(item.reconstruction_residual_w for item in bookings),
            reconstruction_bound_w=max(item.reconstruction_bound_w for item in bookings),
            resolution_bound_w=math.fsum(
                (
                    *(item.resolution_bound_w for item in bookings),
                    math.ulp(term),
                    math.ulp(ledger_top),
                    abs(composite_energy_bound_j) / total,
                )
            ),
            sub_steps=len(points),
        ),
    )


def _through_lower_inlet(
    transaction: through_bed.ValidatedThroughBedK2SingleWallTrayStep,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> _GasPoint:
    boundary = transaction.cell_inputs[0].gas_boundary
    return _gas_point(
        pressure_pa=transaction.accepted_fast_steps[0].state.layer_pressure_pa,
        molar_flow_mol_s=boundary.inlet_molar_flow_mol_s,
        hexane_mole_fraction=boundary.inlet_hexane_mole_fraction,
        water_mole_fraction=boundary.inlet_water_mole_fraction,
        temperature_k=boundary.inlet_temperature_k,
        properties=transaction.request.layer_models[0].properties,
        codec=codec,
    )


def _pd_inlet(
    transaction: pd.ValidatedPDKTrayStep,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> _GasPoint:
    stream = transaction.plenum_ledger.inlet
    return _gas_point(
        pressure_pa=stream.pressure_pa,
        molar_flow_mol_s=stream.molar_flow_mol_s,
        hexane_mole_fraction=stream.hexane_mole_fraction,
        water_mole_fraction=stream.water_mole_fraction,
        temperature_k=stream.temperature_k,
        properties=transaction.request.model.properties,
        codec=codec,
    )


def _handoff(
    source_id: str,
    destination_id: str,
    source: _GasPoint,
    destination: _GasPoint,
    *,
    pressure_limit: float,
) -> QSCMechanismGasHandoffAudit:
    source_hexane = source.molar_flow_mol_s * source.hexane_mole_fraction
    destination_hexane = destination.molar_flow_mol_s * destination.hexane_mole_fraction
    source_water = source.molar_flow_mol_s * source.water_mole_fraction
    destination_water = destination.molar_flow_mol_s * destination.water_mole_fraction
    pressure_residual = math.fsum((source.pressure_pa, -destination.pressure_pa))
    # Q-F2a: book the native arm's pressure-throttling term at the tower
    # seam, exactly as the tray books its internal-face term (CELL-02d W5).
    # The hand-off must be verbatim in every datum-invariant fact first:
    # temperature, both mole fractions and the molar flow bit-identical.
    # Then the ONLY thing that can move the stream's enthalpy is the
    # pressure, and the booked term is the exact frozen-law difference
    # ``F * [h(T, P_destination, y) - h(T, P_source, y)]``.
    throttling_defect = 0.0
    if (
        source.native_arm
        and destination.native_arm
        and destination.pressure_pa != source.pressure_pa
        and destination.temperature_k == source.temperature_k
        and destination.hexane_mole_fraction == source.hexane_mole_fraction
        and destination.water_mole_fraction == source.water_mole_fraction
        and destination.molar_flow_mol_s == source.molar_flow_mol_s
    ):
        throttling_defect = math.fsum(
            (
                destination.common_datum_energy_rate_w,
                -destination.molar_flow_mol_s
                * native_caloric.mixture_enthalpy(
                    destination.temperature_k,
                    source.pressure_pa,
                    destination.hexane_mole_fraction,
                ).value_j_mol,
            )
        )
    # TOP-FACE BOOKING (``DTDC_HANDOFF_TOPFACE_BOOKING``): when the source point
    # carries the term, the energy row starts from the source TRAY LEDGER's top
    # outlet and books the step to the handed pressure beside the Q-F2a defect.
    top_face = source.top_face
    if top_face is None:
        source_energy_w = source.common_datum_energy_rate_w
        top_face_term = 0.0
        energy_terms: tuple[float, ...] = (
            source_energy_w,
            -destination.common_datum_energy_rate_w,
            throttling_defect,
        )
    else:
        source_energy_w = top_face.tray_ledger_top_outlet_energy_w
        top_face_term = top_face.term_w
        energy_terms = (
            source_energy_w,
            top_face_term,
            -destination.common_datum_energy_rate_w,
            throttling_defect,
        )
    residuals = (
        math.fsum((source_hexane, -destination_hexane)),
        math.fsum((source_water, -destination_water)),
        math.fsum(energy_terms),
    )
    hexane_bound = math.fsum(
        math.ulp(value) for value in (source_hexane, destination_hexane, residuals[0])
    )
    water_bound = math.fsum(
        math.ulp(value) for value in (source_water, destination_water, residuals[1])
    )
    energy_bound = math.fsum(
        math.ulp(value)
        for value in (
            source.common_datum_energy_rate_w,
            destination.common_datum_energy_rate_w,
            residuals[2],
            # Q-F2a: the booked throttling term is itself a rounded
            # binary64 product; exact positive zero on the declared arm,
            # where ``math.ulp(0.0)`` is the subnormal minimum and the
            # bound is unmoved at every scale this ledger operates on.
            throttling_defect,
        )
    )
    if top_face is not None:
        # The booked terms' own resolution (and, on a composite, the
        # composition's bound); added only when booked, so a knob-off bound
        # is the same expression as before.
        energy_bound = math.fsum(
            (
                energy_bound,
                math.ulp(source.common_datum_energy_rate_w),
                math.ulp(top_face_term),
                top_face.resolution_bound_w,
            )
        )
    audit = QSCMechanismGasHandoffAudit(
        source_tray_id=source_id,
        destination_tray_id=destination_id,
        pressure_residual_pa=pressure_residual,
        pressure_relative_residual=abs(pressure_residual)
        / max(abs(source.pressure_pa), abs(destination.pressure_pa), 1.0),
        hexane_molar_residual_mol_s=residuals[0],
        water_molar_residual_mol_s=residuals[1],
        common_datum_energy_residual_w=residuals[2],
        hexane_molar_roundoff_bound_mol_s=hexane_bound,
        water_molar_roundoff_bound_mol_s=water_bound,
        common_datum_energy_roundoff_bound_w=energy_bound,
        pressure_relative_limit=pressure_limit,
        booked_pressure_throttling_defect_w=throttling_defect,
        **(
            {}
            if top_face is None
            else {
                "booked_top_face_throttling_w": top_face_term,
                "source_tray_top_outlet_energy_w": source_energy_w,
                "destination_energy_rate_w": destination.common_datum_energy_rate_w,
                "source_sub_steps": top_face.sub_steps,
            }
        ),
    )
    if not audit.passed:
        # A refusal that names its own magnitudes: which of the four rows
        # bit, by how much, against which bound.
        raise QSCK2MechanismIntervalError(
            f"gas handoff {source_id}->{destination_id} failed its causal ledger "
            f"(pressure {audit.pressure_relative_residual:.6e} of "
            f"{audit.pressure_relative_limit:.6e}; "
            f"hexane {audit.hexane_molar_residual_mol_s:.6e} of "
            f"{audit.hexane_molar_roundoff_bound_mol_s:.6e} mol/s; "
            f"water {audit.water_molar_residual_mol_s:.6e} of "
            f"{audit.water_molar_roundoff_bound_mol_s:.6e} mol/s; "
            f"energy {audit.common_datum_energy_residual_w:.6e} of "
            f"{audit.common_datum_energy_roundoff_bound_w:.6e} W)"
        )
    return audit


def _allocation_and_deposition(
    request: QSCK2MechanismIntervalRequest,
    *,
    runtime_authority: mech_runtime.QSCMechanicalRuntimeAuthority,
) -> tuple[
    shaft_work.QSCMechanicalAllocationRequest,
    shaft_work.QSCMechanicalWorkAllocation,
    meal_work.QSCMealMechanicalWorkValidationContext,
    tuple[meal_work.QSCMealMechanicalWorkDepositionBundle, ...],
]:
    trays = request.prior_host_trays
    rank = {tray_id: index for index, tray_id in enumerate(SOLID_TRAY_ORDER)}
    packets = tuple(
        sorted(
            (
                shaft_work.AcceptedQSCMechanicalPacket(
                    physical_tray_id=tray.physical_tray_id,
                    vertical_layer_id=packet.owner_key.vertical_layer_id,
                    packet_id=packet.packet_id,
                    particle_state_identity_digest=packet.payload_identity,
                    inventory=packet.entry.inventory,
                    weight=packet.entry.weight,
                )
                for tray in trays
                for packet in tray.packets
            ),
            key=lambda packet: (
                rank[packet.physical_tray_id],
                packet.vertical_layer_id.value,
                packet.packet_id,
            ),
        )
    )
    endpoint = request.endpoint_plan.common_endpoint_time_s
    accepted_step_digest = (
        "sha256:"
        + hashlib.sha256(
            "\0".join(
                (
                    QSC_K2_MECHANISM_INTERVAL_SCHEMA_ID,
                    request.attempt_id,
                    endpoint.hex(),
                    *(tray.state_digest for tray in trays),
                    request.endpoint_plan.definition_digest,
                )
            ).encode("utf-8")
        ).hexdigest()
    )
    # Ruled 2026-08-30: the allocator law binds the OPERATIVE runtime source,
    # so a constant-speed variant authority pairs with its own variant law;
    # the canonical 11 r/min path is byte-identical to the sealed behavior.
    allocator_law = shaft_work.make_engineering_allocator_law_authority(
        runtime_authority.declared_constant_shaft_rpm
    )
    allocation_request = shaft_work.QSCMechanicalAllocationRequest(
        accepted_step_identity_digest=accepted_step_digest,
        allocator_law_definition_digest=allocator_law.definition_digest,
        energy_datum_id=trays[0].energy_datum_id,
        dt_s=endpoint - request.endpoint_plan.common_start_time_s,
        shaft_power_w=request.shaft_power_w,
        boundary_partition=runtime_authority.nominal_boundary_partition,
        tray_distribution=runtime_authority.nominal_within_meal_distribution,
        packets=packets,
    )
    # Owner-approved B9 (2026-08-31): ONE full runtime validation per
    # interval; every downstream site pins this digest (a value-keyed
    # re-hash catches any mutation) instead of repeating the structural
    # re-walk 50+ times on the same constant object.
    preverified_runtime_digest = mech_runtime.validate_qsc_mechanical_runtime_authority(
        runtime_authority
    )
    allocation = shaft_work.build_qsc_shaft_work_allocation(
        allocation_request,
        runtime_authority=runtime_authority,
        allocator_law_authority=allocator_law,
        preverified_runtime_definition_digest=preverified_runtime_digest,
    )
    context = meal_work.QSCMealMechanicalWorkValidationContext(
        runtime_authority=runtime_authority,
        allocator_law_authority=allocator_law,
        deposition_law_authority=(meal_work.QSC_MEAL_MECHANICAL_WORK_DEPOSITION_LAW_AUTHORITY),
        expected_runtime_definition_digest=preverified_runtime_digest,
        expected_allocator_law_definition_digest=allocator_law.definition_digest,
        preverified_runtime_definition_digest=preverified_runtime_digest,
        expected_deposition_law_definition_digest=(
            meal_work.QSC_MEAL_DEPOSITION_EXPECTED_LAW_DEFINITION_DIGEST
        ),
        expected_allocation_request_definition_digest=allocation_request.definition_digest,
        expected_allocation_definition_digest=allocation.definition_digest,
    )
    bundles = tuple(
        meal_work.build_qsc_meal_mechanical_work_deposition_bundle(
            accepted_tray=tray,
            end_time_s=endpoint,
            allocation_request=allocation_request,
            mechanical_work_allocation=allocation,
            validation_context=context,
        )
        for tray in trays
    )
    return allocation_request, allocation, context, bundles


def _through_request(
    request: QSCK2MechanismIntervalRequest,
    *,
    tray_id: str,
    bottom: ef.EngineeringGasBoundary,
    top_pressure_pa: float,
    source_id: str | None,
    source_digest: str | None,
    bundle: meal_work.QSCMealMechanicalWorkDepositionBundle,
    bundle_digest: str | None = None,
) -> through_bed.ThroughBedK2SingleWallStepRequest:
    index = THROUGH_BED_SOLID_ORDER.index(tray_id)
    state = request.through_bed_states[index]
    return through_bed.ThroughBedK2SingleWallStepRequest(
        attempt_id=f"{request.attempt_id}:physics:{tray_id}",
        prior=state,
        integration_law_definition_digest=(
            through_bed.SINGLE_WALL_EXPECTED_INTEGRATION_LAW_DEFINITION_DIGEST
        ),
        end_time_s=request.endpoint_plan.common_endpoint_time_s,
        layer_models=request.through_bed_layer_models[index],  # type: ignore[arg-type]
        bottom_gas_boundary=bottom,
        bottom_gas_binding=through_bed.ThroughBedGasBoundaryBinding(
            physical_tray_id=tray_id,
            boundary_kind=state.contract.lower_gas_boundary_kind,
            source_physical_tray_id=source_id,
            authority_id=GAS_HANDOFF_AUTHORITY_ID,
            source_state_digest=source_digest,
        ),
        top_boundary_pressure_pa=top_pressure_pa,
        flow_segments=request.flow_histories[SOLID_TRAY_ORDER.index(tray_id)],
        meal_mechanical_work_deposition=bundle,
        preverified_meal_bundle_digest=bundle_digest,
    )


def _evaluate_through_chain(
    request: QSCK2MechanismIntervalRequest,
    *,
    sp1_mn2_pressure_pa: float,
    mn2_mn1_pressure_pa: float,
    codec: high_loading.HighLoadingV1PacketAdapter,
    planners: dict[str, k_scheduler.EngineeringPhysicalKRekeyPlanner],
    bundles: dict[str, meal_work.QSCMealMechanicalWorkDepositionBundle],
    bundle_digests: dict[str, str] | None = None,
    work_context: meal_work.QSCMealMechanicalWorkValidationContext,
    prepared_cache: dict[str, through_bed.PreparedThroughBedTrayContext] | None = None,
) -> tuple[through_bed.ValidatedThroughBedK2SingleWallTrayStep, ...]:
    def evaluate(
        tray_id: str,
        bottom: ef.EngineeringGasBoundary,
        top_pressure_pa: float,
        source_id: str | None,
        source_digest: str | None,
    ) -> through_bed.ValidatedThroughBedK2SingleWallTrayStep:
        index = THROUGH_BED_SOLID_ORDER.index(tray_id)
        state = request.through_bed_states[index]
        result = through_bed.evaluate_through_bed_k2_single_wall_tray_step(
            _through_request(
                request,
                tray_id=tray_id,
                bottom=bottom,
                top_pressure_pa=top_pressure_pa,
                source_id=source_id,
                source_digest=source_digest,
                bundle=bundles[tray_id],
                bundle_digest=(None if bundle_digests is None else bundle_digests.get(tray_id)),
            ),
            codec=codec,
            k_planner=planners[tray_id],
            integration_law_authority=through_bed.SINGLE_WALL_INTEGRATION_LAW_AUTHORITY,
            expected_integration_law_definition_digest=(
                through_bed.SINGLE_WALL_EXPECTED_INTEGRATION_LAW_DEFINITION_DIGEST
            ),
            expected_layer_model_configuration_digest=(state.layer_model_configuration_digest),
            meal_mechanical_work_validation_context=work_context,
            prepared_cache=prepared_cache,
        )
        if type(result) is not through_bed.ValidatedThroughBedK2SingleWallTrayStep:
            if getattr(result, "phase_exhaustion", None) is not None:
                raise QSCK2MechanismPhaseExhaustionError(
                    f"{tray_id} rejected: {result.reason}",
                    signal=result.phase_exhaustion,
                )
            raise QSCK2MechanismIntervalError(f"{tray_id} rejected: {result.reason}")
        return result

    sp1 = evaluate(
        "SP1",
        request.sp1_bottom_gas_boundary,
        sp1_mn2_pressure_pa,
        None,
        None,
    )
    sp1_out = _through_outlet(sp1, codec=codec)
    mn2 = evaluate(
        "MN2",
        ef.EngineeringGasBoundary(
            inlet_molar_flow_mol_s=sp1_out.molar_flow_mol_s,
            inlet_hexane_mole_fraction=sp1_out.hexane_mole_fraction,
            inlet_water_mole_fraction=sp1_out.water_mole_fraction,
            inlet_temperature_k=sp1_out.temperature_k,
            downstream_boundary_pressure_pa=sp1_mn2_pressure_pa,
            carrier_topology=ef.CarrierTopology.BINARY_NO_INERT,
        ),
        mn2_mn1_pressure_pa,
        "SP1",
        sp1.candidate.state_digest,
    )
    mn2_out = _through_outlet(mn2, codec=codec)
    mn1 = evaluate(
        "MN1",
        ef.EngineeringGasBoundary(
            inlet_molar_flow_mol_s=mn2_out.molar_flow_mol_s,
            inlet_hexane_mole_fraction=mn2_out.hexane_mole_fraction,
            inlet_water_mole_fraction=mn2_out.water_mole_fraction,
            inlet_temperature_k=mn2_out.temperature_k,
            downstream_boundary_pressure_pa=mn2_mn1_pressure_pa,
            carrier_topology=ef.CarrierTopology.BINARY_NO_INERT,
        ),
        MN1_TOP_PRESSURE_PA,
        "MN2",
        mn2.candidate.state_digest,
    )
    return mn1, mn2, sp1


#: The four clauses of the through-chain pressure admission test, named in the
#: order the test evaluates them.  REPORT-ONLY: the names exist so a refusal can
#: say WHICH clause it failed, and nothing reads them to decide anything.
ORDERED_STRUCTURAL_DOMAIN_CLAUSES = ("low_bound", "ordering", "high_bound", "mn1_top")


def _ordered_structural_domain_clauses(
    *,
    sp1_mn2_pressure_pa: float,
    mn2_mn1_pressure_pa: float,
    low_pa: float,
    high_pa: float,
) -> tuple[tuple[str, bool, float, str], ...]:
    """Split the admission test into its four clauses, each with its margin.

    WHY THIS EXISTS (2026-09-25, report-only).  Five production-point walls on
    five lanes since 2026-09-19 have been attributed to
    "pressure Picard iterate left the ordered structural domain" without anyone
    knowing which of the four comparisons was crossed or in which direction,
    because the refusal named none of the quantities it tested and no Picard
    PROPOSAL is recorded anywhere in a march record.  This function names them.

    It changes NOTHING.  The clauses below are the same comparisons
    :func:`_close_through_pressures` writes as one expression, taken one at a
    time; the caller still raises on that ONE unchanged expression, this
    function is never consulted to decide anything, and the exception type, the
    condition, the evaluation limit and the control flow are untouched.  Only
    the message text is longer.

    Each row is ``(clause name, whether it HOLDS, its admissibility margin in
    pascal, how it reads when it fails)``.  The margin is signed towards
    admissibility, so a failing clause's ``abs(margin)`` is its shortfall, and
    an exact tie on one of the two STRICT clauses reports a shortfall of
    exactly zero rather than being rounded into the band.  A non-finite iterate
    fails every comparison, so it is reported as failing every clause with a
    non-finite margin instead of being described as admissible.
    """

    return (
        (
            "low_bound",
            low_pa <= mn2_mn1_pressure_pa,
            mn2_mn1_pressure_pa - low_pa,
            "p(MN2|MN1) below the admitted low bound",
        ),
        (
            "ordering",
            mn2_mn1_pressure_pa < sp1_mn2_pressure_pa,
            sp1_mn2_pressure_pa - mn2_mn1_pressure_pa,
            "p(MN2|MN1) not below p(SP1|MN2)",
        ),
        (
            "high_bound",
            sp1_mn2_pressure_pa <= high_pa,
            high_pa - sp1_mn2_pressure_pa,
            "p(SP1|MN2) above the admitted high bound",
        ),
        (
            "mn1_top",
            MN1_TOP_PRESSURE_PA < mn2_mn1_pressure_pa,
            mn2_mn1_pressure_pa - MN1_TOP_PRESSURE_PA,
            "p(MN2|MN1) not above the MN1 top pressure",
        ),
    )


def _ordered_structural_domain_detail(
    *,
    sp1_mn2_pressure_pa: float,
    mn2_mn1_pressure_pa: float,
    low_pa: float,
    high_pa: float,
    evaluation: int,
    evaluation_limit: int,
) -> str:
    """The report-only suffix the ordered-domain refusal carries.

    Returns the evaluation index, both iterate values, the declared band and
    the MN1 top pressure, and every clause that failed with its direction and
    its shortfall.  Report-only and side-effect free; see
    :func:`_ordered_structural_domain_clauses` for what is and is not changed.
    """

    clauses = _ordered_structural_domain_clauses(
        sp1_mn2_pressure_pa=sp1_mn2_pressure_pa,
        mn2_mn1_pressure_pa=mn2_mn1_pressure_pa,
        low_pa=low_pa,
        high_pa=high_pa,
    )
    failed = [row for row in clauses if not row[1]]
    if failed:
        printed = "; ".join(
            f"{name} ({reading} by {abs(margin):.6e} Pa)" for name, _, margin, reading in failed
        )
    else:
        # Unreachable from the caller's own condition; printed rather than
        # asserted so a future divergence between the two is visible instead
        # of silent.
        printed = "no clause reported as failing"
    return (
        f" (evaluation {evaluation} of {evaluation_limit};"
        f" p(SP1|MN2)={sp1_mn2_pressure_pa:.6e} Pa,"
        f" p(MN2|MN1)={mn2_mn1_pressure_pa:.6e} Pa;"
        f" admitted band [{low_pa:.6e}, {high_pa:.6e}] Pa,"
        f" MN1 top {MN1_TOP_PRESSURE_PA:.6e} Pa;"
        f" failed: {printed})"
    )


def _close_through_pressures(
    request: QSCK2MechanismIntervalRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    planners: dict[str, k_scheduler.EngineeringPhysicalKRekeyPlanner],
    bundles: dict[str, meal_work.QSCMealMechanicalWorkDepositionBundle],
    bundle_digests: dict[str, str] | None = None,
    work_context: meal_work.QSCMealMechanicalWorkValidationContext,
) -> tuple[
    tuple[through_bed.ValidatedThroughBedK2SingleWallTrayStep, ...],
    QSCMechanismPressureSolveEvidence,
]:
    face_pressures = (
        request.through_bed_states[1].cell_field_states[0].layer_pressure_pa,
        request.through_bed_states[0].cell_field_states[0].layer_pressure_pa,
    )
    # B0 hoist: one per-attempt, per-tray prepared record, filled by the first
    # chain evaluation and value-re-verified by the rest.  Nothing in it
    # depends on the face pressures, so every iterate reads the identical
    # objects the un-hoisted path would have rebuilt.
    prepared_cache: dict[str, through_bed.PreparedThroughBedTrayContext] | None = (
        {} if through_bed.HOIST_ENABLED else None
    )
    low, high = ef.ENGINEERING_SIX_LAYER_SERIAL_LAW2_PRESSURE_ADMISSION_PA
    for evaluation in range(1, request.pressure_root_max_evaluations + 1):
        sp1_mn2, mn2_mn1 = face_pressures
        if not low <= mn2_mn1 < sp1_mn2 <= high or not MN1_TOP_PRESSURE_PA < mn2_mn1:
            raise QSCK2MechanismIntervalError(
                "pressure Picard iterate left the ordered structural domain"
                + _ordered_structural_domain_detail(
                    sp1_mn2_pressure_pa=sp1_mn2,
                    mn2_mn1_pressure_pa=mn2_mn1,
                    low_pa=low,
                    high_pa=high,
                    evaluation=evaluation,
                    evaluation_limit=request.pressure_root_max_evaluations,
                )
            )
        transactions = _evaluate_through_chain(
            request,
            sp1_mn2_pressure_pa=sp1_mn2,
            mn2_mn1_pressure_pa=mn2_mn1,
            codec=codec,
            planners=planners,
            bundles=bundles,
            bundle_digests=bundle_digests,
            work_context=work_context,
            prepared_cache=prepared_cache,
        )
        mn1, mn2, _sp1 = transactions
        proposed = (
            mn2.accepted_fast_steps[0].state.layer_pressure_pa,
            mn1.accepted_fast_steps[0].state.layer_pressure_pa,
        )
        residuals = (proposed[0] - sp1_mn2, proposed[1] - mn2_mn1)
        relative = max(
            abs(residuals[0]) / max(abs(proposed[0]), abs(sp1_mn2), 1.0),
            abs(residuals[1]) / max(abs(proposed[1]), abs(mn2_mn1), 1.0),
        )
        if relative <= request.pressure_face_relative_limit:
            evidence = QSCMechanismPressureSolveEvidence(
                evaluations=evaluation,
                sp1_mn2_face_pressure_pa=sp1_mn2,
                mn2_mn1_face_pressure_pa=mn2_mn1,
                sp1_mn2_residual_pa=residuals[0],
                mn2_mn1_residual_pa=residuals[1],
                maximum_relative_residual=relative,
                relative_limit=request.pressure_face_relative_limit,
            )
            return transactions, evidence
        face_pressures = proposed
    raise QSCK2MechanismIntervalError("pressure Picard iteration exceeded its evaluation limit")


def _apply_pd3_inlet_vapour_composition_override(
    request: QSCK2MechanismIntervalRequest,
    mn1_outlet: _GasPoint,
    *,
    properties: object,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> _GasPoint:
    """WS-0: apply the declared PD3 inlet vapour composition override, once.

    ``None`` (the default) returns ``mn1_outlet`` UNCHANGED -- today's
    behavior exactly.  A declared override replaces ONLY the mole
    fractions; molar flow, pressure and temperature stay physics-derived,
    and the common-datum enthalpy rate is RECOMPUTED (via the same
    ``_gas_point`` seam every other tray-boundary point uses) so the point
    stays internally self-consistent.  This is applied to the ONE
    ``mn1_outlet`` value the caller then threads through every consumer of
    "what MN1 hands off to PD3" -- the actual PD3 bypass physics
    (``_evaluate_pd_chain``), the MN1->PD3 gas-handoff conservation ledger,
    and the zero-contact network's independent cross-check -- so all three
    read the identical (possibly overridden) boundary and never disagree.
    """

    if request.pd3_inlet_vapour_composition_override is None:
        return mn1_outlet
    hexane_mole_fraction, water_mole_fraction = request.pd3_inlet_vapour_composition_override
    return _gas_point(
        pressure_pa=mn1_outlet.pressure_pa,
        molar_flow_mol_s=mn1_outlet.molar_flow_mol_s,
        hexane_mole_fraction=hexane_mole_fraction,
        water_mole_fraction=water_mole_fraction,
        temperature_k=mn1_outlet.temperature_k,
        properties=properties,
        codec=codec,
    )


def _evaluate_pd_chain(
    request: QSCK2MechanismIntervalRequest,
    *,
    mn1_outlet: _GasPoint,
    codec: high_loading.HighLoadingV1PacketAdapter,
    planners: dict[str, k_scheduler.EngineeringPhysicalKRekeyPlanner],
    bundles: dict[str, meal_work.QSCMealMechanicalWorkDepositionBundle],
    bundle_digests: dict[str, str] | None = None,
    work_context: meal_work.QSCMealMechanicalWorkValidationContext,
) -> tuple[pd.ValidatedPDKTrayStep, ...]:
    stream = pd.PDBinaryGasStream(
        pressure_pa=mn1_outlet.pressure_pa,
        molar_flow_mol_s=mn1_outlet.molar_flow_mol_s,
        hexane_mole_fraction=mn1_outlet.hexane_mole_fraction,
        water_mole_fraction=mn1_outlet.water_mole_fraction,
        temperature_k=mn1_outlet.temperature_k,
    )
    results: dict[str, pd.ValidatedPDKTrayStep] = {}
    for tray_id in PD_CAUSAL_ORDER:
        index = SOLID_TRAY_ORDER.index(tray_id)
        model = request.pd_models[index]
        pd_request = pd.PDKTrayStepRequest(
            attempt_id=f"{request.attempt_id}:physics:{tray_id}",
            prior=request.pd_states[index],
            end_time_s=request.endpoint_plan.common_endpoint_time_s,
            model=model,
            bypass_inlet=stream,
            prescribed_bypass_pressure_drop_pa=0.0,
            accepted_indirect_duty_w=model.geometry.tray.indirect_duty_w,
            flow_segments=request.flow_histories[index],
            meal_mechanical_work_deposition=bundles[tray_id],
            preverified_meal_bundle_digest=(
                None if bundle_digests is None else bundle_digests.get(tray_id)
            ),
        )
        result = pd.evaluate_pd_k_dome_bypass_tray_step(
            pd_request,
            codec=codec,
            k_planner=planners[tray_id],
            meal_mechanical_work_validation_context=work_context,
        )
        if type(result) is not pd.ValidatedPDKTrayStep:
            raise QSCK2MechanismIntervalError(f"{tray_id} rejected: {result.reason}")
        results[tray_id] = result
        stream = result.plenum_ledger.outlet
    return results["PD1"], results["PD2"], results["PD3"]


def _solids_accumulation(
    request: QSCK2MechanismIntervalRequest,
    draw_transaction: accepted_draw.ValidatedQSCAcceptedDrawExit,
) -> QSCMechanismSolidsAccumulationEvidence:
    feed_dry = math.fsum(
        derive_packet_totals(item.entry.inventory, item.entry.weight).totals[0]
        for item in request.feed_packets
    )
    draws = tuple(item.accepted_dry_matter_kg for item in draw_transaction.tray_ledgers)
    incoming = (feed_dry, *draws[:-1])
    accumulation = tuple(
        math.fsum((entry, -exit_mass)) for entry, exit_mass in zip(incoming, draws, strict=True)
    )
    return QSCMechanismSolidsAccumulationEvidence(
        feed_dry_matter_kg=feed_dry,
        accepted_draw_dry_matter_kg=draws,
        tray_routing_accumulation_dry_matter_kg=accumulation,
    )


def evaluate_qsc_six_tray_k2_mechanism_interval(
    request: QSCK2MechanismIntervalRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    runtime_authority: mech_runtime.QSCMechanicalRuntimeAuthority,
    draw_evaluator: accepted_draw.QSCAcceptedDrawExitEvaluator,
    reconciler: k_reconciler.EngineeringKTowerRouteReconciler,
    admission_planners: tuple[tower_router.ManufacturedAdmissionPlanner, ...],
) -> ValidatedQSCK2MechanismInterval | QSCK2MechanismIntervalRejection:
    """Evaluate one prepared interval and return a routed, uncommitted candidate."""

    if type(request) is not QSCK2MechanismIntervalRequest:
        raise TypeError("mechanism evaluator requires an exact request")
    try:
        QSCK2MechanismIntervalRequest.__post_init__(request)
        if type(codec) is not high_loading.HighLoadingV1PacketAdapter and (
            type(codec) is not dual_variant.EngineeringDualVariantPacketAdapter
        ):
            raise TypeError(
                "mechanism evaluator requires the exact tag-1 packet adapter "
                "or the owner-ruled dual-variant adapter"
            )
        if type(runtime_authority) is not mech_runtime.QSCMechanicalRuntimeAuthority:
            raise TypeError("mechanism evaluator requires an exact mechanical runtime")
        if type(draw_evaluator) is not accepted_draw.QSCAcceptedDrawExitEvaluator:
            raise TypeError("mechanism evaluator requires an exact accepted-draw evaluator")
        if type(reconciler) is not k_reconciler.EngineeringKTowerRouteReconciler:
            raise TypeError("mechanism evaluator requires an exact K reconciler")
        if runtime_authority.nominal_boundary_partition is not (
            mech_runtime.ShaftWorkBoundaryPartition.ALL_TO_MEAL
        ):
            raise QSCK2MechanismIntervalError(
                "current mechanism interval requires the runtime all-to-meal endpoint"
            )
        states = tuple(row.prior_state for row in request.endpoint_plan.tray_endpoints)
        k_endpoint.validate_qsc_six_tray_k2_event_endpoint_plan(
            request.endpoint_plan,
            states=states,
            runtime_authority=runtime_authority,
            requested_horizon_time_s=request.endpoint_plan.requested_horizon_time_s,
        )
        allocation_request, allocation, work_context, deposition = _allocation_and_deposition(
            request, runtime_authority=runtime_authority
        )
        planners = {
            row.physical_tray_id: row.planner for row in request.endpoint_plan.tray_endpoints
        }
        bundle_map = {item.physical_tray_id: item for item in deposition}
        # B7-family memo (2026-08-31 rulings): each bundle's digest is
        # derived ONCE per attempt here; every per-Picard-iteration
        # validation below pins this digest instead of re-running the
        # full reconstruction (a mutated bundle re-hashes and refuses).
        bundle_digests = {tray_id: item.definition_digest for tray_id, item in bundle_map.items()}
        through_transactions, pressure_solve = _close_through_pressures(
            request,
            codec=codec,
            planners=planners,
            bundles=bundle_map,
            bundle_digests=bundle_digests,
            work_context=work_context,
        )
        mn1, mn2, sp1 = through_transactions
        through_handoffs = (
            _handoff(
                "SP1",
                "MN2",
                _handoff_source_outlet(sp1, codec=codec),
                _through_lower_inlet(mn2, codec=codec),
                pressure_limit=request.pressure_face_relative_limit,
            ),
            _handoff(
                "MN2",
                "MN1",
                _handoff_source_outlet(mn2, codec=codec),
                _through_lower_inlet(mn1, codec=codec),
                pressure_limit=request.pressure_face_relative_limit,
            ),
        )
        mn1_outlet = _handoff_source_outlet(mn1, codec=codec)
        # WS-0 (owner-ruled AUTHORIZED 2026-09-04): apply the declared PD3
        # inlet vapour composition override ONCE, here, to the ONE
        # ``mn1_outlet`` value every downstream consumer of "what MN1 hands
        # off to PD3" reads -- the actual PD3 bypass physics below, the
        # MN1->PD3 gas-handoff conservation ledger, and the zero-contact
        # network's own independent cross-check.  ``None`` (the default)
        # returns ``mn1_outlet`` unchanged, so every one of those three
        # readers reproduces today's behavior EXACTLY; a declared override
        # moves all three identically, so none of them can disagree with
        # another about what PD3's inlet is.
        mn1_outlet = _apply_pd3_inlet_vapour_composition_override(
            request,
            mn1_outlet,
            properties=mn1.request.layer_models[1].properties,
            codec=codec,
        )
        pd_transactions = _evaluate_pd_chain(
            request,
            bundle_digests=bundle_digests,
            mn1_outlet=mn1_outlet,
            codec=codec,
            planners=planners,
            bundles=bundle_map,
            work_context=work_context,
        )
        pd1, pd2, pd3 = pd_transactions
        mn1_pd3 = _handoff(
            "MN1",
            "PD3",
            mn1_outlet,
            _pd_inlet(pd3, codec=codec),
            pressure_limit=request.pressure_face_relative_limit,
        )
        pd3_inlet = pd_network.QSCPDBinaryGasRate(
            pressure_pa=mn1_outlet.pressure_pa,
            water_molar_flow_mol_s=(mn1_outlet.molar_flow_mol_s * mn1_outlet.water_mole_fraction),
            hexane_molar_flow_mol_s=(
                mn1_outlet.molar_flow_mol_s * mn1_outlet.hexane_mole_fraction
            ),
            common_datum_enthalpy_rate_w=mn1_outlet.common_datum_energy_rate_w,
            energy_datum_id=request.prior_host_trays[0].energy_datum_id,
        )
        zero_contact_network = pd_zero_contact.couple_qsc_pd_zero_contact_network(
            (pd3, pd2, pd1),
            pd3_inlet=pd3_inlet,
        )
        local_endpoint_trays = (
            pd1.candidate.tray,
            pd2.candidate.tray,
            pd3.candidate.tray,
            mn1.candidate.tray,
            mn2.candidate.tray,
            sp1.candidate.tray,
        )
        shaft_segments = runtime_authority.k_authorities[0].shaft_segments
        draw_request = accepted_draw.QSCAcceptedDrawExitRequest(
            attempt_id=f"{request.attempt_id}:accepted-draw",
            start_time_s=request.endpoint_plan.common_start_time_s,
            end_time_s=request.endpoint_plan.common_endpoint_time_s,
            accepted_prior_trays=request.prior_host_trays,
            routing_endpoint_trays=local_endpoint_trays,
            flow_histories=request.flow_histories,
            expected_flow_history_digests=tuple(
                accepted_draw.qsc_accepted_draw_flow_history_digest(
                    physical_tray_id=tray_id,
                    segments=history,
                )
                for tray_id, history in zip(SOLID_TRAY_ORDER, request.flow_histories, strict=True)
            ),
            accepted_shaft_segments=shaft_segments,
            expected_shaft_history_digest=(
                accepted_draw.qsc_accepted_draw_shaft_history_digest(shaft_segments)
            ),
            closing_layer_cursors=request.closing_layer_cursors,
            feed_packets=request.feed_packets,
            k_share_closing_exchange_policy=request.k_share_closing_exchange_policy,
        )
        draw_transaction = draw_evaluator.evaluate(
            draw_request,
            admission_planners=admission_planners,
        )
        if type(draw_transaction) is not accepted_draw.ValidatedQSCAcceptedDrawExit:
            raise QSCK2MechanismIntervalError(
                f"accepted draw/router rejected: {draw_transaction.reason}"
            )
        reconciliation_request = k_reconciler.EngineeringKTowerRouteReconciliationRequest(
            attempt_id=f"{request.attempt_id}:k-route-reconciliation",
            routing_transaction=draw_transaction.routing_transaction,
            prior_schedule_states=tuple(
                row.synchronized_step.candidate for row in request.endpoint_plan.tray_endpoints
            ),
        )
        reconciliation = reconciler.evaluate(reconciliation_request)
        if type(reconciliation) is not (k_reconciler.ValidatedEngineeringKTowerRouteReconciliation):
            raise QSCK2MechanismIntervalError(
                f"K route reconciliation rejected: {reconciliation.reason}"
            )
        routed = draw_transaction.routing_transaction.candidate_trays
        candidate_pd = tuple(
            replace(transaction.candidate, tray=routed[index])
            for index, transaction in enumerate(pd_transactions)
        )
        candidate_through = tuple(
            replace(transaction.candidate, tray=routed[index + 3])
            for index, transaction in enumerate(through_transactions)
        )
        prepared_candidate = PreparedQSCK2MechanismCandidate(
            pd_states=candidate_pd,
            through_bed_states=candidate_through,
            k_schedule_states=reconciliation.candidate_schedule_states,
            _seal=_PREPARED_CANDIDATE_SEAL,
        )
        result = ValidatedQSCK2MechanismInterval(
            request=request,
            endpoint_plan=request.endpoint_plan,
            allocation_request=allocation_request,
            mechanical_work_allocation=allocation,
            deposition_bundles=deposition,
            pressure_solve=pressure_solve,
            through_bed_transactions=through_transactions,
            pd_transactions=pd_transactions,
            gas_handoffs=(*through_handoffs, mn1_pd3),
            pd_zero_contact_network=zero_contact_network,
            accepted_draw_transaction=draw_transaction,
            reconciliation_transaction=reconciliation,
            solids_accumulation=_solids_accumulation(request, draw_transaction),
            prepared_candidate=prepared_candidate,
            _seal=_VALIDATED_INTERVAL_SEAL,
        )
        if not result.passed and not (
            EXPLORATION_LEDGER_GATE is not None and EXPLORATION_LEDGER_GATE.interval_ok(result)
        ):
            raise QSCK2MechanismIntervalError(
                "prepared mechanism interval contains a failed conservative ledger"
            )
        return result
    except Exception as error:
        return QSCK2MechanismIntervalRejection(
            attempt_id=request.attempt_id,
            reason=f"{type(error).__name__}: {error}",
            rollback_pd_states=request.pd_states,
            rollback_through_bed_states=request.through_bed_states,
            phase_exhaustion=(
                error.signal if isinstance(error, QSCK2MechanismPhaseExhaustionError) else None
            ),
        )


def tower_gas_enthalpy_telescope(
    result: ValidatedQSCK2MechanismInterval,
) -> dict[str, object] | None:
    """REPORT-ONLY: the tower's gas enthalpy telescope, sparge to PD3 inlet.

    Present only when the hand-offs carry the top-face booking
    (``DTDC_HANDOFF_TOPFACE_BOOKING``); ``None`` otherwise, so a knob-off
    record gains nothing.  The chain, every term a number an audit or a tray
    ledger already holds:

        PD3 inlet = SP1 ledger bottom inlet
                    + sum over SP1, MN2, MN1 of (ledger top outlet - bottom inlet)
                    + sum over the three seams of (reconstruction + top-face term
                      + Q-F2a defect - hand-off residual)
                    + sum over the two internal seams of (next ledger's bottom
                      inlet - handed destination)

    The telescope residual is that identity's left minus right side, with its
    ULP bound; ``unbooked_residual_w`` is what the same telescope reads with the
    top-face terms left out (their sum, the diagnosis's term).  A seam whose
    source outlet was composed from several sub-steps is not comparable (its
    tray ledger belongs to one sub-step), and ``comparable`` says so.
    """

    handoffs = {
        (item.source_tray_id, item.destination_tray_id): item for item in result.gas_handoffs
    }
    seams = (("SP1", "MN2"), ("MN2", "MN1"), ("MN1", "PD3"))
    if any(seam not in handoffs for seam in seams) or any(
        handoffs[seam].source_tray_top_outlet_energy_w is None for seam in seams
    ):
        return None
    ledgers = {}
    for transaction in result.through_bed_transactions:
        whole = transaction.whole_tray_energy
        inner = getattr(whole, "thermal_ledger", None)
        ledgers[transaction.request.prior.tray.physical_tray_id] = whole if inner is None else inner
    sparge_in = ledgers["SP1"].bottom_inlet_energy_w
    pd3_in = handoffs[("MN1", "PD3")].destination_energy_rate_w
    gains = {
        tray_id: math.fsum(
            (ledgers[tray_id].top_outlet_energy_w, -ledgers[tray_id].bottom_inlet_energy_w)
        )
        for tray_id in ("SP1", "MN2", "MN1")
    }
    seam_rows: dict[str, dict[str, object]] = {}
    terms: list[float] = [pd3_in, -sparge_in, *(-value for value in gains.values())]
    top_face_terms: list[float] = []
    for source_id, destination_id in seams:
        audit = handoffs[(source_id, destination_id)]
        reconstruction = math.fsum(
            (audit.source_tray_top_outlet_energy_w, -ledgers[source_id].top_outlet_energy_w)
        )
        link = (
            0.0
            if destination_id == "PD3"
            else math.fsum(
                (ledgers[destination_id].bottom_inlet_energy_w, -audit.destination_energy_rate_w)
            )
        )
        terms.extend(
            (
                -reconstruction,
                -audit.booked_top_face_throttling_w,
                -audit.booked_pressure_throttling_defect_w,
                audit.common_datum_energy_residual_w,
                -link,
            )
        )
        top_face_terms.append(audit.booked_top_face_throttling_w)
        seam_rows[f"{source_id}->{destination_id}"] = {
            "top_face_w": audit.booked_top_face_throttling_w,
            "q_f2a_defect_w": audit.booked_pressure_throttling_defect_w,
            "handoff_residual_w": audit.common_datum_energy_residual_w,
            "handoff_bound_w": audit.common_datum_energy_roundoff_bound_w,
            "tray_ledger_reconstruction_w": reconstruction,
            "link_to_next_ledger_w": link,
            "source_sub_steps": audit.source_sub_steps,
        }
    residual = math.fsum(terms)
    top_face_sum = math.fsum(top_face_terms)
    bound = math.fsum(
        (
            *(math.ulp(value) for value in terms),
            *(math.ulp(ledger.top_outlet_energy_w) for ledger in ledgers.values()),
            *(math.ulp(ledger.bottom_inlet_energy_w) for ledger in ledgers.values()),
            math.ulp(residual),
        )
    )
    return {
        "law_id": HANDOFF_TOPFACE_BOOKING_LAW_ID,
        "report_only": True,
        "comparable": all(handoffs[seam].source_sub_steps == 1 for seam in seams),
        "sparge_inlet_w": sparge_in,
        "pd3_inlet_w": pd3_in,
        "tray_gain_w": gains,
        "seams": seam_rows,
        "top_face_sum_w": top_face_sum,
        "telescope_residual_w": residual,
        "telescope_bound_w": bound,
        "closed": abs(residual) <= bound,
        "unbooked_residual_w": math.fsum((residual, top_face_sum)),
    }


__all__ = (
    "GAS_HANDOFF_AUTHORITY_ID",
    "MN1_TOP_PRESSURE_PA",
    "PD_CAUSAL_ORDER",
    "PreparedQSCK2MechanismCandidate",
    "QSC_K2_MECHANISM_INTERVAL_SCHEMA_ID",
    "QSCK2MechanismIntervalError",
    "QSCK2MechanismIntervalRejection",
    "QSCK2MechanismIntervalRequest",
    "QSCMechanismGasHandoffAudit",
    "QSCMechanismPressureSolveEvidence",
    "QSCMechanismSolidsAccumulationEvidence",
    "SOLID_TRAY_ORDER",
    "THROUGH_BED_SOLID_ORDER",
    "ValidatedQSCK2MechanismInterval",
    "evaluate_qsc_six_tray_k2_mechanism_interval",
    "tower_gas_enthalpy_telescope",
)
