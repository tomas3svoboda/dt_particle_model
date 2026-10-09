"""Nonqualifying K=2 through-bed tray seam with one physical wall node.

MN1, MN2, and SP1 each have one physical floor-facing tray wall.  The lower
K cell therefore owns the full record-only reference duty, the active wall
transfer coefficients, and the sole dynamic wall capacity.  The upper K cell
uses an exact quasi-steady placeholder with zero wall, steam, and ambient
conductance.  Both fast Law-2 solves see the same accepted old wall
temperature, but only one aggregate wall row is advanced after both fast
solves have been accepted.

This is an engineering prequalification seam.  It does not identify wall
parameters, realize the record-only reference duty, add shaft work, prove an
inter-tray gas equality, authenticate machine data, or advance any physical,
industrial, production, calibration, release, or plant-predictive claim.
"""

from __future__ import annotations

import hashlib
import math
import os
import struct
from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import ClassVar, Iterable

from . import cell_closure as cc
from . import cell_engineering_feasibility as ef
from . import component_energy_datum_adapter as component_datum
from . import engineering_physical_k_event_scheduler as k_scheduler
from . import engineering_dual_variant_packet_adapter as dual_variant
from . import high_loading_packet_adapter as high_loading
from . import k_cell_tray_host as tray_host
from . import layer_phase_exhaustion as exhaustion
from . import qsc_meal_mechanical_work_deposition as meal_work
from . import sp1_k2_law2_tray_integration as legacy
from . import through_bed_k2_law2_kernel as kernel
from . import through_bed_k2_law2_tray_integration as two_wall

SINGLE_WALL_INTEGRATION_DIGEST_DOMAIN = (
    "GT-PS-2/through-bed-k2-rtd8-law2-single-physical-wall-integration/v1"
)
SINGLE_WALL_COUPLED_ENERGY_RESOLUTION_LAW_ID = (
    "through-bed-k2-single-wall-coupled-binary64-resolution-v1"
)
SINGLE_WALL_Q_OTHER_S_EXTERNAL_ENERGY_RESOLUTION_LAW_ID = (
    "through-bed-k2-single-wall-q-other-s-external-energy-binary64-resolution-v1"
)
SINGLE_WALL_INTEGRATION_LAW_SCHEMA_ID = "dtdc-core2-through-bed-k2-single-wall-integration-law-v1"
SINGLE_WALL_INTEGRATION_LAW_SCHEMA_REVISION = 1
SINGLE_WALL_EXPECTED_INTEGRATION_LAW_DEFINITION_DIGEST = (
    "sha256:a5626a1a44a309e17ea4561c5eb2454ff92496d0c6f7b076fe72c52108d51824"
)
#: CELL-02d W7: the native lane's frozen expected digest (additive; the
#: declared digest above is immutable).
SINGLE_WALL_EXPECTED_NATIVE_INTEGRATION_LAW_DEFINITION_DIGEST = (
    "sha256:8813de165bd120e2d0c538bb50d350e5c81c0de3a2b7902bef1570fc4dda367e"
)

#: B0 hoist (speed program 2026-09-02).  Every validation and inventory the
#: through-bed tray step derives from ``(request.prior, layer_models, codec,
#: k_planner, bundle, end_time_s)`` is IDENTICAL across every iteration of the
#: pressure Picard closure, because ``request.through_bed_states[index]`` and
#: the layer models are fixed for the whole ``_close_through_pressures`` call.
#: They are prepared ONCE per attempt and threaded in, with a value-keyed
#: re-hash on every use so a mutated prepared context still refuses.  Exactly
#: "1" disables the hoist and restores the un-hoisted path for the A/B.
HOIST_DISABLE_ENVIRONMENT_VARIABLE = "DTDC_DISABLE_HOIST"
HOIST_ENABLED = os.environ.get(HOIST_DISABLE_ENVIRONMENT_VARIABLE) != "1"
PREPARED_THROUGH_BED_CONTEXT_LAW_ID = (
    "through-bed-k2-single-wall-prepared-interval-context-v1"
)

ThroughBedKLayerEngineeringModel = legacy.SP1KLayerEngineeringModel
ReferenceThroughBedK2Contract = kernel.ReferenceThroughBedK2Contract
ThroughBedGasBoundaryBinding = kernel.ThroughBedGasBoundaryBinding
ThroughBedGasBoundaryKind = kernel.ThroughBedGasBoundaryKind
InterKGasFaceLedger = legacy.InterKGasFaceLedger
AggregateThroughBedPressureLedger = legacy.AggregateSP1PressureLedger
PacketDerivedLayerInventory = legacy.PacketDerivedLayerInventory
LayerAcceptedTransfers = legacy.LayerAcceptedTransfers
LayerPacketExchangeLedger = legacy.LayerPacketExchangeLedger
WholeTrayExternalEnergyLedger = legacy.WholeTrayExternalEnergyLedger
NativeWholeTrayExternalEnergyLedger = legacy.NativeWholeTrayExternalEnergyLedger
reference_through_bed_k2_contract = kernel.reference_through_bed_k2_contract

#: F-HULL-1, owner ruling 2026-09-06 ~13:55 ("I authorize the ruling, proceed
#: with implementation", verbatim in
#: ``docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md``): the printed pressure-drop
#: hull is a PER-TRAY figure.  The published 1-6 kPa is printed per
#: COUNTERCURRENT TRAY (envelope B17), and a through-bed tray of K layers sums
#: to it; enforcing the same numbers on each layer halves the admitted figure
#: at K = 2 and refused the ruled operating point itself (SP1 652 + 653 Pa,
#: tray 1305 Pa, inside the band).  The band is IMPORTED from the certified
#: holder and never retyped: this module introduces no numeric literal, no new
#: band, and no change to any law, rate, or output -- only the unit of account
#: of a gate.
PRINTED_TRAY_DROP_HULL_PA = cc.CORROBORATED_LAYER_DROP_HULL_PA
PRINTED_TRAY_DROP_HULL_UNIT_OF_ACCOUNT_LAW_ID = (
    "through-bed-printed-drop-hull-per-tray-sum-of-k-layers-v1"
)
PRINTED_TRAY_DROP_HULL_UNIT_OF_ACCOUNT_SOURCE = (
    "owner ruling 2026-09-06 ~13:55 (F-HULL-1), verbatim in "
    "docs/GT_PS2_OWNER_BATCH_QUEUE_2026-09-01.md; measured basis F-HULL-1 and "
    "F-SCALE-1 MEASURED in the same ledger (at the ruled operating point the "
    "per-tray drops are 1.30 / 2.01 / 3.21 kPa, inside the printed band, while "
    "the sparge tray's two layers at 0.65 kPa each fail the per-layer floor).  "
    "The band values themselves stay where they were certified: "
    "src/dtdc_simulator/core2/cell_closure.py:386 "
    "(CORROBORATED_LAYER_DROP_HULL_PA), the hull of the frozen Ergun band and "
    "the independently corroborated printed band."
)
#: D15 (owner ruling 2026-09-19, verbatim "D15 I rule as you suggest"; packet
#: ``docs/GT_PS2_RULING_PACKET_D15_2026-09-19.md`` section 2, built as its
#: adversarial audit's construction B,
#: ``docs/evidence/ruling_audit_2026-09-19/D15v2/AUDIT_D15v2.md``, on the
#: census ``docs/GT_PS2_T11_WALL_LEDGER_CENSUS_RECORD_2026-09-19.md``): the
#: multiple of the aggregate wall row's propagated binary64 resolution bound
#: that the aggregate physical-wall ledger admits.
#:
#: What the row is.  The wall temperature is ONE DIVISION, numerator over the
#: denominator ``C/dt + UA_steam + UA_ambient``, and that denominator is
#: exactly the row's derivative with respect to the wall temperature (the
#: audit checks it by finite difference to the last bit on the bit-for-bit
#: reproduced farm refusal), so the row residual divided by the derivative
#: times the wall temperature's ulp is the distance, in ulps of the wall
#: temperature, between the rounded quotient and the row's exact root.  The
#: UNMULTIPLIED bound built by ``_aggregate_wall_row_resolution_bound`` is
#: 2.000 of those derivative-ulps on every one of the census's 15,577 rows
#: (2.0000000032 to 2.0006612, at every stride decade from 1e-6 s to 1 s), so
#: this multiple admits a landing of four ulps.
#:
#: The row residual is a CHECK, not a correction: no state value at a refused
#: step changes when the admission widens, only whether the step is accepted.
#: The relative clause and ``request.relative_limit`` are untouched, as are
#: the two cancellation contributions, which refuse on their own exact-zero
#: requirement whatever this multiple is.
WALL_ROW_RESOLUTION_BOUND_MULTIPLE = 2.0
WALL_ROW_RESOLUTION_BOUND_MULTIPLE_LAW_ID = (
    "through-bed-k2-single-wall-aggregate-wall-row-resolution-bound-multiple-v1"
)
#: The bracket this value was picked from, as a multiple of the unmultiplied
#: bound.  Both ends are MEASUREMENTS, as D10-b's were
#: (``docs/GT_PS2_D10B_TEMPERATURE_ALLOWANCE_96_BUILD_RECORD_2026-09-14.md``):
#: 1.125 is the maximum landing the exactly reproduced row produced over
#: 6,260,000 evaluations (2.250 ulps; the smallest multiple that admits every
#: landing ON RECORD is 1.0117, the farm refusal), and 1.94 is the row's
#: ANALYTIC ceiling, two numerator roundings of at most 1.686 ulps each seen
#: through the derivative near 430 K plus the division's half ulp, 3.87 ulps
#: in all.  The declared value is the smallest WHOLE multiple above the
#: bracket, as 96 sat above D10-b's 90.17 to 97.54 by the same construction.
WALL_ROW_RESOLUTION_BOUND_MULTIPLE_BRACKET = (1.125, 1.94)
WALL_ROW_RESOLUTION_BOUND_MULTIPLE_SOURCE = (
    "owner ruling 2026-09-19, verbatim \"D15 I rule as you suggest\", on "
    "docs/GT_PS2_RULING_PACKET_D15_2026-09-19.md section 2 (D15-a, Class B "
    "under D3 of docs/GT_PS2_FALLING_RATE_CONTINUATION_PLAN_2026-09-12.md "
    "section 5, under P-1 of 2026-09-03 and the numerical-materiality "
    "philosophy of 2026-08-30).  Measured basis: the T11 census "
    "docs/GT_PS2_T11_WALL_LEDGER_CENSUS_RECORD_2026-09-19.md -- 15,577 "
    "aggregate ledgers stamped on MN1, MN2 and SP1 across two workstation "
    "lanes, none refused, ratio mean 0.23 to 0.25, p99 0.71 to 0.76, maximum "
    "0.9591, correlation with the macro step, the row derivative and the "
    "wall temperature all below 0.09 -- against the single refusal on record "
    "(farm 1, MN1, step 468 at 31.425 s: residual 2.9729e-4 W against a bound "
    "of 2.9386e-4 W, ratio 1.0117, a landing at 2.023 ulps), and the audit's "
    "6,260,000 reproduced evaluations of the same row whose maximum landing "
    "is 2.250 ulps with nothing above 2.5.  Materiality of the one refusal, "
    "in the philosophy's units: the residual beyond the bound is 3.4248e-6 W, "
    "over the 2.9015e-4 s step 9.937e-10 J and 1.325e-15 K on the wall's "
    "750,000 J/K capacity; the whole residual is 8.626e-8 J and 1.150e-13 K.  "
    "The pinned literals are untouched: DEFAULT_RELATIVE_LIMIT in "
    "src/dtdc_simulator/core2/sp1_k2_law2_tray_integration.py and everything "
    "in src/dtdc_simulator/core2/through_bed_k_single_wall_thermohydraulic.py "
    "are only read, and request.relative_limit is passed through unchanged."
)
DECLARED_ENGINEERING_ASSUMPTION = "DECLARED_ENGINEERING_ASSUMPTION"
DECLARED_ASSUMPTIONS = MappingProxyType(
    {
        "printed_drop_hull_unit_of_account": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "criterion": (
                    "the printed hull is a per-TRAY figure; the K layers of one "
                    "through-bed tray sum to it.  The gate is applied to the sum "
                    "of the K CONVERGED layer drops, never to an iterate and no "
                    "longer to a single layer"
                ),
                "unit_of_account": "one physical countercurrent tray",
                "value_pa": PRINTED_TRAY_DROP_HULL_PA,
                "law_id": PRINTED_TRAY_DROP_HULL_UNIT_OF_ACCOUNT_LAW_ID,
                "sources": PRINTED_TRAY_DROP_HULL_UNIT_OF_ACCOUNT_SOURCE,
                "rejecting_outcome": (
                    "a converged tray drop outside [1000.0, 5600.0] Pa is refused "
                    "typed as cell_closure.PressureBandError, naming the tray sum, "
                    "its K layer drops and the band; the step is rejected, never "
                    "clamped"
                ),
                "new_numeric_literal": False,
            }
        ),
        #: D15-a (owner ruling 2026-09-19).  The construction is the audit's
        #: construction B, the only one that passes the ledger's own
        #: self-consistency clause: ``_aggregate_physical_wall_ledger`` stores
        #: the MULTIPLIED bound in the existing field
        #: ``wall_row_binary64_resolution_bound_w``, from which both ``passed``
        #: properties re-derive the residual beyond resolution, and records the
        #: unmultiplied bound and the multiple beside it so the ledger says
        #: what was admitted against what.  The coupled ledger copies all
        #: three.  The stored bound therefore changes on every evaluation,
        #: relieved or not; the residual, the derivative, the temperature and
        #: every march-record leaf do not.
        "wall_row_resolution_bound_multiple": MappingProxyType(
            {
                "class": DECLARED_ENGINEERING_ASSUMPTION,
                "name": "aggregate_physical_wall_row_resolution_bound_multiple",
                "criterion": (
                    "the aggregate physical-wall ledger admits the wall row "
                    "when its residual lies within the multiple times the "
                    "row's two-ulp binary64 resolution bound, the multiple "
                    "being the smallest WHOLE one that admits every landing "
                    "the row's own arithmetic can produce.  The row is one "
                    "division, so its landing error has an arithmetic ceiling "
                    "of 3.87 ulps of the wall temperature (1.94 of the bound) "
                    "and a measured ceiling of 2.250 ulps (1.125 of the "
                    "bound) over 6,260,000 reproduced evaluations"
                ),
                "value": WALL_ROW_RESOLUTION_BOUND_MULTIPLE,
                "bracket_multiple_of_unmultiplied_bound": (
                    WALL_ROW_RESOLUTION_BOUND_MULTIPLE_BRACKET
                ),
                "law_id": WALL_ROW_RESOLUTION_BOUND_MULTIPLE_LAW_ID,
                "sources": WALL_ROW_RESOLUTION_BOUND_MULTIPLE_SOURCE,
                "rejecting_outcome": (
                    "written first, four of them.  (1) Any committed "
                    "march-record leaf differing on a lane that never met the "
                    "gate, measured on identity twins at the ceremony.  "
                    "(2) Any stamped landing above the ANALYTIC ceiling, 3.87 "
                    "ulps of the wall temperature, that is a ratio above 1.94 "
                    "of the unmultiplied bound, on any host: it says the row "
                    "is not the closed-form division this declaration "
                    "describes, or the host's arithmetic is not binary64 "
                    "round-to-nearest, and the handling is withdrawn with "
                    "that evaluation as the record.  At a multiple of 2 such "
                    "a landing is ADMITTED and the T11 capture stamps it, so "
                    "this outcome can fire and be seen.  (3) The census on a "
                    "relieved lane, run through the T11 knob on the same "
                    "host, with the ratio's mean outside 0.20 to 0.30 or its "
                    "99th percentile outside 0.65 to 0.80 per tray, which "
                    "would say the arithmetic changed, which this relief "
                    "cannot cause.  (4) A relieved lane refusing within a few "
                    "intervals of its former refusal on the coupled ledger's "
                    "OTHER sub-ledgers, whose relative clause passes at "
                    "exactly the limit on accepted ledgers, which would say "
                    "the item is not this row.  Nothing here is clamped: a "
                    "row outside the multiplied bound refuses the step"
                ),
                "new_numeric_literal": True,
            }
        ),
    }
)


class ThroughBedK2SingleWallIntegrationError(ValueError):
    """Base refusal for the one-wall through-bed transaction."""


class ThroughBedK2SingleWallConfigurationError(ThroughBedK2SingleWallIntegrationError):
    """Static topology, model, boundary, or accepted state is malformed."""


class ThroughBedK2SingleWallStepError(ThroughBedK2SingleWallIntegrationError):
    """A one-wall macro-step is inadmissible."""


class StaleThroughBedK2SingleWallStateError(ThroughBedK2SingleWallIntegrationError):
    """Commit was attempted against a different accepted state."""


def _require_finite(name: str, value: float, *, positive: bool = False) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise ThroughBedK2SingleWallConfigurationError(
            f"{name} must be a finite exact binary64 value"
        )
    if positive and value <= 0.0:
        raise ThroughBedK2SingleWallConfigurationError(f"{name} must be strictly positive")


def _require_sha256(name: str, value: str) -> None:
    if (
        type(value) is not str
        or len(value) != 71
        or not value.startswith("sha256:")
        or any(character not in "0123456789abcdef" for character in value[7:])
    ):
        raise ThroughBedK2SingleWallConfigurationError(
            f"{name} must be a canonical sha256 identity"
        )


def _digest_part(value: str | int | float | bool) -> tuple[bytes, bytes]:
    if type(value) is str:
        return b"s", value.encode("utf-8")
    if type(value) is bool:
        return b"b", b"1" if value else b"0"
    if type(value) is int:
        return b"i", str(value).encode("ascii")
    if type(value) is float:
        return b"f", struct.pack(">d", value)
    raise TypeError(f"unsupported single-wall digest part: {type(value).__name__}")


def _digest(domain: str, parts: Iterable[str | int | float | bool]) -> str:
    digest = hashlib.sha256()
    digest.update(b"DTDC-THROUGH-BED-K2-SINGLE-WALL-INTEGRATION-V1\0")
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


_CANONICAL_REFERENCE_CONTRACT_DIGESTS = tuple(
    kernel.reference_through_bed_k2_contract(tray_id).definition_digest
    for tray_id in ("MN1", "MN2", "SP1")
)


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedK2SingleWallIntegrationLawAuthority:
    """Typed source/law pins for the bounded one-wall integration seam.

    The numerical wall declarations remain caller-supplied engineering
    assumptions.  This authority binds only the operative topology, routing,
    transfer/allocation conventions, and binary64 resolution laws used to
    interpret those declarations.
    """

    schema_id: str
    schema_revision: int
    integration_digest_domain: str
    coupled_energy_resolution_law_id: str
    source_kernel_schema_id: str
    source_legacy_integration_digest_domain: str
    source_gas_boundary_enthalpy_convention: str
    source_packet_thermal_aggregation_law_id: str
    source_packet_energy_allocation_law_id: str
    source_packet_common_energy_resolution_law_id: str
    source_whole_tray_energy_resolution_law_id: str
    supported_reference_trays: tuple[str, ...]
    reference_contract_definition_digests: tuple[str, ...]
    through_bed_layer_count: int
    through_bed_rtd_stage_count: int
    physical_wall_host_layer: int
    physical_wall_count: int
    inactive_upper_wall_closure: str
    reference_indirect_duty_record_only: bool
    wall_parameters_declared_engineering_assumptions: bool
    wall_parameters_identified: bool

    declared_engineering_law: ClassVar[bool] = True
    authenticated_machine_data: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    industrially_feasible: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self.schema_id != ("dtdc-core2-through-bed-k2-single-wall-integration-law-v1"):
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall law authority has a foreign schema identity"
            )
        if type(self.schema_revision) is not int or self.schema_revision != 1:
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall law authority has a foreign schema revision"
            )
        expected_strings = (
            (
                self.integration_digest_domain,
                "GT-PS-2/through-bed-k2-rtd8-law2-single-physical-wall-integration/v1",
            ),
            (
                self.coupled_energy_resolution_law_id,
                "through-bed-k2-single-wall-coupled-binary64-resolution-v1",
            ),
            (
                self.source_kernel_schema_id,
                "dtdc-core2-through-bed-k2-rtd8-law2-kernel-v1",
            ),
            (
                self.source_legacy_integration_digest_domain,
                "GT-PS-2/sp1-k2-law2-tray-integration/v1",
            ),
            (
                self.source_packet_thermal_aggregation_law_id,
                "static_uniform_representative_particle_ua_v1",
            ),
            (
                self.source_packet_energy_allocation_law_id,
                "packet_local_conduction_plus_representative_particle_weight_share_remainder_v1",
            ),
            (
                self.source_packet_common_energy_resolution_law_id,
                "weighted_before_after_common_energy_nextafter_spacing_v1",
            ),
            (
                self.source_whole_tray_energy_resolution_law_id,
                "weighted_packet_cell_state_and_external_operation_nextafter_budget_v1",
            ),
            (self.inactive_upper_wall_closure, "quasi_steady_limit"),
        )
        # CELL-02d W7: the enthalpy convention is a closed two-member union -
        # the declared constant-property convention or the native property-law
        # convention - never anything else. The declared lane's digest and
        # checks are unchanged; the native lane is additive.
        if (
            self.source_gas_boundary_enthalpy_convention
            not in (
                legacy.GAS_BOUNDARY_ENTHALPY_CONVENTION,
                legacy.NATIVE_GAS_BOUNDARY_ENTHALPY_CONVENTION,
            )
            or type(self.source_gas_boundary_enthalpy_convention) is not str
        ):
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall law authority has a foreign gas-boundary enthalpy convention"
            )
        if any(
            type(observed) is not str or observed != expected
            for observed, expected in expected_strings
        ):
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall law authority has a foreign operative law identity"
            )
        if (
            type(self.supported_reference_trays) is not tuple
            or self.supported_reference_trays != ("MN1", "MN2", "SP1")
            or any(type(item) is not str for item in self.supported_reference_trays)
        ):
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall law authority has a foreign reference-tray set or order"
            )
        if (
            type(self.reference_contract_definition_digests) is not tuple
            or len(self.reference_contract_definition_digests) != 3
        ):
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall law authority requires three canonical contract pins"
            )
        for digest in self.reference_contract_definition_digests:
            _require_sha256("single-wall reference contract definition", digest)
        if self.reference_contract_definition_digests != _CANONICAL_REFERENCE_CONTRACT_DIGESTS:
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall law authority binds foreign reference contracts"
            )
        if (
            type(self.through_bed_layer_count) is not int
            or type(self.through_bed_rtd_stage_count) is not int
            or type(self.physical_wall_host_layer) is not int
            or type(self.physical_wall_count) is not int
            or (
                self.through_bed_layer_count,
                self.through_bed_rtd_stage_count,
                self.physical_wall_host_layer,
                self.physical_wall_count,
            )
            != (2, 8, 1, 1)
        ):
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall law authority topology or wall ownership drifted"
            )
        if (
            type(self.reference_indirect_duty_record_only) is not bool
            or type(self.wall_parameters_declared_engineering_assumptions) is not bool
            or type(self.wall_parameters_identified) is not bool
            or self.reference_indirect_duty_record_only is not True
            or self.wall_parameters_declared_engineering_assumptions is not True
            or self.wall_parameters_identified is not False
        ):
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall law authority overstates its duty or wall-parameter status"
            )

    @property
    def definition_digest(self) -> str:
        return _digest(
            "through-bed-k2-single-wall-integration-law-authority",
            (
                self.schema_id,
                self.schema_revision,
                self.integration_digest_domain,
                self.coupled_energy_resolution_law_id,
                self.source_kernel_schema_id,
                self.source_legacy_integration_digest_domain,
                self.source_gas_boundary_enthalpy_convention,
                self.source_packet_thermal_aggregation_law_id,
                self.source_packet_energy_allocation_law_id,
                self.source_packet_common_energy_resolution_law_id,
                self.source_whole_tray_energy_resolution_law_id,
                *self.supported_reference_trays,
                *self.reference_contract_definition_digests,
                self.through_bed_layer_count,
                self.through_bed_rtd_stage_count,
                self.physical_wall_host_layer,
                self.physical_wall_count,
                self.inactive_upper_wall_closure,
                self.reference_indirect_duty_record_only,
                self.wall_parameters_declared_engineering_assumptions,
                self.wall_parameters_identified,
            ),
        )


SINGLE_WALL_INTEGRATION_LAW_AUTHORITY = ThroughBedK2SingleWallIntegrationLawAuthority(
    schema_id=SINGLE_WALL_INTEGRATION_LAW_SCHEMA_ID,
    schema_revision=SINGLE_WALL_INTEGRATION_LAW_SCHEMA_REVISION,
    integration_digest_domain=SINGLE_WALL_INTEGRATION_DIGEST_DOMAIN,
    coupled_energy_resolution_law_id=SINGLE_WALL_COUPLED_ENERGY_RESOLUTION_LAW_ID,
    source_kernel_schema_id=kernel.THROUGH_BED_K2_KERNEL_SCHEMA_ID,
    source_legacy_integration_digest_domain=legacy.INTEGRATION_DIGEST_DOMAIN,
    source_gas_boundary_enthalpy_convention=legacy.GAS_BOUNDARY_ENTHALPY_CONVENTION,
    source_packet_thermal_aggregation_law_id=(
        legacy.UNIFORM_REPRESENTATIVE_PARTICLE_UA_AGGREGATION_LAW_ID
    ),
    source_packet_energy_allocation_law_id=(
        legacy.PACKET_LOCAL_CONDUCTION_WEIGHT_SHARE_REMAINDER_ENERGY_ALLOCATION_LAW_ID
    ),
    source_packet_common_energy_resolution_law_id=(
        legacy.PACKET_COMMON_ENERGY_BINARY64_RESOLUTION_LAW_ID
    ),
    source_whole_tray_energy_resolution_law_id=(
        legacy.WHOLE_TRAY_ENERGY_BINARY64_RESOLUTION_LAW_ID
    ),
    supported_reference_trays=kernel.SUPPORTED_REFERENCE_TRAYS,
    reference_contract_definition_digests=_CANONICAL_REFERENCE_CONTRACT_DIGESTS,
    through_bed_layer_count=kernel.THROUGH_BED_LAYER_COUNT,
    through_bed_rtd_stage_count=kernel.THROUGH_BED_RTD_STAGE_COUNT,
    physical_wall_host_layer=1,
    physical_wall_count=1,
    inactive_upper_wall_closure=cc.WallClosure.QUASI_STEADY_LIMIT.value,
    reference_indirect_duty_record_only=True,
    wall_parameters_declared_engineering_assumptions=True,
    wall_parameters_identified=False,
)

#: CELL-02d W7: the native-mode law authority - identical operative laws with
#: the native property-law enthalpy convention in the one convention slot. The
#: declared authority, its digest, and its checks are byte-untouched.
SINGLE_WALL_NATIVE_INTEGRATION_LAW_AUTHORITY = ThroughBedK2SingleWallIntegrationLawAuthority(
    schema_id=SINGLE_WALL_INTEGRATION_LAW_SCHEMA_ID,
    schema_revision=SINGLE_WALL_INTEGRATION_LAW_SCHEMA_REVISION,
    integration_digest_domain=SINGLE_WALL_INTEGRATION_DIGEST_DOMAIN,
    coupled_energy_resolution_law_id=SINGLE_WALL_COUPLED_ENERGY_RESOLUTION_LAW_ID,
    source_kernel_schema_id=kernel.THROUGH_BED_K2_KERNEL_SCHEMA_ID,
    source_legacy_integration_digest_domain=legacy.INTEGRATION_DIGEST_DOMAIN,
    source_gas_boundary_enthalpy_convention=legacy.NATIVE_GAS_BOUNDARY_ENTHALPY_CONVENTION,
    source_packet_thermal_aggregation_law_id=(
        legacy.UNIFORM_REPRESENTATIVE_PARTICLE_UA_AGGREGATION_LAW_ID
    ),
    source_packet_energy_allocation_law_id=(
        legacy.PACKET_LOCAL_CONDUCTION_WEIGHT_SHARE_REMAINDER_ENERGY_ALLOCATION_LAW_ID
    ),
    source_packet_common_energy_resolution_law_id=(
        legacy.PACKET_COMMON_ENERGY_BINARY64_RESOLUTION_LAW_ID
    ),
    source_whole_tray_energy_resolution_law_id=(
        legacy.WHOLE_TRAY_ENERGY_BINARY64_RESOLUTION_LAW_ID
    ),
    supported_reference_trays=kernel.SUPPORTED_REFERENCE_TRAYS,
    reference_contract_definition_digests=_CANONICAL_REFERENCE_CONTRACT_DIGESTS,
    through_bed_layer_count=kernel.THROUGH_BED_LAYER_COUNT,
    through_bed_rtd_stage_count=kernel.THROUGH_BED_RTD_STAGE_COUNT,
    physical_wall_host_layer=1,
    physical_wall_count=1,
    inactive_upper_wall_closure=cc.WallClosure.QUASI_STEADY_LIMIT.value,
    reference_indirect_duty_record_only=True,
    wall_parameters_declared_engineering_assumptions=True,
    wall_parameters_identified=False,
)


def validate_single_wall_integration_law_authority(
    authority: ThroughBedK2SingleWallIntegrationLawAuthority,
) -> str:
    """Refuse any typed-law reissue, mutation, or operative global substitution."""

    if type(authority) is not ThroughBedK2SingleWallIntegrationLawAuthority:
        raise TypeError("single-wall law must be an exact authority record")
    ThroughBedK2SingleWallIntegrationLawAuthority.__post_init__(authority)
    digest = authority.definition_digest
    if authority.source_gas_boundary_enthalpy_convention == (
        legacy.NATIVE_GAS_BOUNDARY_ENTHALPY_CONVENTION
    ):
        expected_digest = SINGLE_WALL_EXPECTED_NATIVE_INTEGRATION_LAW_DEFINITION_DIGEST
        expected_convention = legacy.NATIVE_GAS_BOUNDARY_ENTHALPY_CONVENTION
    else:
        expected_digest = SINGLE_WALL_EXPECTED_INTEGRATION_LAW_DEFINITION_DIGEST
        expected_convention = legacy.GAS_BOUNDARY_ENTHALPY_CONVENTION
    if digest != expected_digest:
        raise ThroughBedK2SingleWallConfigurationError(
            "single-wall law differs from its fixed expected definition digest"
        )
    live_law_values = (
        SINGLE_WALL_INTEGRATION_DIGEST_DOMAIN,
        SINGLE_WALL_COUPLED_ENERGY_RESOLUTION_LAW_ID,
        kernel.THROUGH_BED_K2_KERNEL_SCHEMA_ID,
        legacy.INTEGRATION_DIGEST_DOMAIN,
        expected_convention,
        legacy.UNIFORM_REPRESENTATIVE_PARTICLE_UA_AGGREGATION_LAW_ID,
        legacy.PACKET_LOCAL_CONDUCTION_WEIGHT_SHARE_REMAINDER_ENERGY_ALLOCATION_LAW_ID,
        legacy.PACKET_COMMON_ENERGY_BINARY64_RESOLUTION_LAW_ID,
        legacy.WHOLE_TRAY_ENERGY_BINARY64_RESOLUTION_LAW_ID,
    )
    pinned_law_values = (
        authority.integration_digest_domain,
        authority.coupled_energy_resolution_law_id,
        authority.source_kernel_schema_id,
        authority.source_legacy_integration_digest_domain,
        authority.source_gas_boundary_enthalpy_convention,
        authority.source_packet_thermal_aggregation_law_id,
        authority.source_packet_energy_allocation_law_id,
        authority.source_packet_common_energy_resolution_law_id,
        authority.source_whole_tray_energy_resolution_law_id,
    )
    if live_law_values != pinned_law_values:
        raise ThroughBedK2SingleWallConfigurationError(
            "operative single-wall integration or source-law identity was substituted"
        )
    if (
        kernel.SUPPORTED_REFERENCE_TRAYS != authority.supported_reference_trays
        or kernel.THROUGH_BED_LAYER_COUNT != authority.through_bed_layer_count
        or kernel.THROUGH_BED_RTD_STAGE_COUNT != authority.through_bed_rtd_stage_count
        or cc.WallClosure.QUASI_STEADY_LIMIT.value != authority.inactive_upper_wall_closure
    ):
        raise ThroughBedK2SingleWallConfigurationError(
            "operative single-wall topology or inactive-wall closure was substituted"
        )
    live_contract_digests = tuple(
        kernel.reference_through_bed_k2_contract(tray_id).definition_digest
        for tray_id in authority.supported_reference_trays
    )
    if live_contract_digests != authority.reference_contract_definition_digests:
        raise ThroughBedK2SingleWallConfigurationError(
            "operative single-wall reference contracts were substituted"
        )
    if (
        reference_through_bed_k2_contract is not kernel.reference_through_bed_k2_contract
        or ReferenceThroughBedK2Contract is not kernel.ReferenceThroughBedK2Contract
        or ThroughBedKLayerEngineeringModel is not legacy.SP1KLayerEngineeringModel
    ):
        raise ThroughBedK2SingleWallConfigurationError(
            "operative single-wall contract or model type was substituted"
        )
    return digest


def _validate_expected_single_wall_integration_law(
    authority: ThroughBedK2SingleWallIntegrationLawAuthority,
    expected_definition_digest: str,
) -> str:
    _require_sha256("expected single-wall integration-law definition", expected_definition_digest)
    observed = validate_single_wall_integration_law_authority(authority)
    if observed != expected_definition_digest:
        raise ThroughBedK2SingleWallConfigurationError(
            "single-wall law differs from the caller-held expected definition digest"
        )
    return observed


def _same_binary64(left: float, right: float) -> bool:
    return struct.pack(">d", left) == struct.pack(">d", right)


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedTrayPrintedDropCorroboration:
    """Where the TRAY drop sits against the printed hull, as the sources print it.

    The tray drop is the sum of the CONVERGED drops of the tray's K layers.
    Reported in full - the layer drops, their sum, the band, the verdict - so a
    reader never has to recover the unit of account from the verdict alone.
    """

    law_id: str
    layer_drop_pa: tuple[float, ...]
    tray_drop_sum_pa: float
    printed_hull_pa: tuple[float, float]
    inside_printed_hull: bool

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def classify_tray_pressure_drop(
    layer_drop_pa: tuple[float, ...],
) -> ThroughBedTrayPrintedDropCorroboration:
    """Report, never refuse: the tray sum against the printed hull."""

    if type(layer_drop_pa) is not tuple or not layer_drop_pa:
        raise ThroughBedK2SingleWallIntegrationError(
            "the per-tray printed-hull check reads a non-empty tuple of converged layer drops"
        )
    for drop in layer_drop_pa:
        if type(drop) is not float or not math.isfinite(drop):
            raise ThroughBedK2SingleWallIntegrationError(
                "every converged layer drop must be an exact finite float"
            )
    low, high = PRINTED_TRAY_DROP_HULL_PA
    # The same fsum, in the same canonical layer order, that the aggregate
    # pressure ledger already computes as ``layer_drop_sum_pa``; the step
    # asserts the two agree bit-for-bit rather than trusting the coincidence.
    tray_drop_sum_pa = math.fsum(layer_drop_pa)
    return ThroughBedTrayPrintedDropCorroboration(
        law_id=PRINTED_TRAY_DROP_HULL_UNIT_OF_ACCOUNT_LAW_ID,
        layer_drop_pa=layer_drop_pa,
        tray_drop_sum_pa=tray_drop_sum_pa,
        printed_hull_pa=PRINTED_TRAY_DROP_HULL_PA,
        inside_printed_hull=low <= tray_drop_sum_pa <= high,
    )


def require_corroborated_tray_drop(
    layer_drop_pa: tuple[float, ...],
) -> ThroughBedTrayPrintedDropCorroboration:
    """Fail closed outside the printed hull, on the TRAY (F-HULL-1).

    The successor of ``cell_closure.require_corroborated_layer_drop`` for the
    through-bed path, with the ruled unit of account: the same imported band,
    the same message family, the sum of the K converged layer drops in place of
    one layer's drop.
    """

    corroboration = classify_tray_pressure_drop(layer_drop_pa)
    if not corroboration.inside_printed_hull:
        low, high = PRINTED_TRAY_DROP_HULL_PA
        raise cc.PressureBandError(
            f"the converged tray drop {corroboration.tray_drop_sum_pa!r} Pa - the "
            f"sum of the {len(layer_drop_pa)} converged layer drops "
            f"{corroboration.layer_drop_pa!r} - falls outside the printed hull "
            f"[{low}, {high}] Pa; this tray introduces no new pressure band "
            f"({PRINTED_TRAY_DROP_HULL_UNIT_OF_ACCOUNT_LAW_ID})"
        )
    return corroboration


def _tray_unit_of_account_cell_inputs(
    inputs: ef.EngineeringCellInputs,
) -> ef.EngineeringCellInputs:
    """Retire the PER-LAYER printed-hull raise on the through-bed path.

    F-HULL-1 execution, composed around two pinned modules and editing neither.
    ``sp1_k2_law2_tray_integration._cell_inputs`` (PINNED, :1347) sets
    ``require_corroborated_pressure_drop`` to ``model.declared_lane_drop_band_pa
    is None``, and ``cell_engineering_feasibility`` (PINNED, :2191-2192) raises
    ``cell_closure.require_corroborated_layer_drop`` on that flag.  The flag is
    the kernel's OWN public opt-out (``cell_closure.py:79``), so the caller
    clears it on the exact dataclass the kernel reads - the identical idiom the
    frozen K-successor seam already uses at
    ``through_bed_k_single_wall_thermohydraulic.py:755-764``, where the printed
    band is likewise enforced on the tray aggregate.  No declared lane band is
    invented: the B8 band keeps its own meaning and still gates per layer
    wherever a model declares one.
    """

    if inputs.require_corroborated_pressure_drop is False:
        return inputs
    return replace(inputs, require_corroborated_pressure_drop=False)


def inactive_upper_wall_parameters(
    physical_wall: cc.WallNodeParameters,
) -> cc.WallNodeParameters:
    """Return the exact inert K2 placeholder derived from the physical wall.

    ``WallNodeParameters`` requires a positive capacity even for its
    quasi-steady limit.  The retained scalar is consequently only a schema
    coordinate: ``QUASI_STEADY_LIMIT`` makes its effective storage exactly
    zero, and this module never advances it.
    """

    if type(physical_wall) is not cc.WallNodeParameters:
        raise ThroughBedK2SingleWallConfigurationError(
            "physical wall must be an exact WallNodeParameters"
        )
    return replace(
        physical_wall,
        steam_side_ua_w_k=0.0,
        ambient_ua_w_k=0.0,
        closure=cc.WallClosure.QUASI_STEADY_LIMIT,
    )


def inactive_upper_transfer_coefficients(
    physical_transfer: cc.TransferCoefficients,
) -> cc.TransferCoefficients:
    """Return the exact K2 transfer declaration with both wall paths absent."""

    if type(physical_transfer) is not cc.TransferCoefficients:
        raise ThroughBedK2SingleWallConfigurationError(
            "physical transfer must be exact TransferCoefficients"
        )
    return replace(
        physical_transfer,
        wall_to_interface_ua_w_k=0.0,
        wall_to_gas_ua_w_k=0.0,
    )


def _validate_single_wall_layer_models(
    contract: ReferenceThroughBedK2Contract,
    models: tuple[ThroughBedKLayerEngineeringModel, ThroughBedKLayerEngineeringModel],
) -> None:
    if type(contract) is not ReferenceThroughBedK2Contract:
        raise ThroughBedK2SingleWallConfigurationError(
            "single-wall validation requires an exact reference contract"
        )
    if (
        type(models) is not tuple
        or len(models) != kernel.THROUGH_BED_LAYER_COUNT
        or any(type(item) is not ThroughBedKLayerEngineeringModel for item in models)
    ):
        raise ThroughBedK2SingleWallConfigurationError(
            "single-wall validation requires two exact layer models in canonical K order"
        )
    try:
        kernel._validate_through_bed_layer_geometry(contract, models)
        kernel._validate_one_gas_enthalpy_convention(models)
    except kernel.ThroughBedK2KernelError as error:
        raise ThroughBedK2SingleWallConfigurationError(str(error)) from error

    physical, upper = models
    reference = kernel._reference_tray(contract.physical_tray_id)
    duties = (
        physical.geometry.tray.indirect_duty_w,
        upper.geometry.tray.indirect_duty_w,
    )
    if duties != (reference.indirect_duty_w, 0.0):
        raise ThroughBedK2SingleWallConfigurationError(
            "the floor-facing K1 layer must carry the full record-only reference duty "
            "and K2 must carry exact zero"
        )
    if physical.wall.closure is not cc.WallClosure.DYNAMIC_NODE:
        raise ThroughBedK2SingleWallConfigurationError(
            "the sole physical K1 wall must use DYNAMIC_NODE"
        )
    if physical.wall.steam_side_ua_w_k <= 0.0:
        raise ThroughBedK2SingleWallConfigurationError(
            "the sole physical K1 wall requires positive steam-side conductance"
        )
    if (
        physical.transfer.wall_to_interface_ua_w_k == 0.0
        and physical.transfer.wall_to_gas_ua_w_k == 0.0
    ):
        raise ThroughBedK2SingleWallConfigurationError(
            "the sole physical K1 wall requires at least one active wall-transfer path"
        )
    expected_upper_wall = inactive_upper_wall_parameters(physical.wall)
    if upper.wall != expected_upper_wall:
        raise ThroughBedK2SingleWallConfigurationError(
            "K2 wall must be the exact inactive QUASI_STEADY_LIMIT placeholder "
            "derived from K1; duplicate capacity, steam, or ambient ownership is refused"
        )
    expected_upper_transfer = inactive_upper_transfer_coefficients(physical.transfer)
    if upper.transfer != expected_upper_transfer:
        raise ThroughBedK2SingleWallConfigurationError(
            "K2 transfer must be the exact K1-derived declaration with zero wall UAs"
        )


def canonical_single_wall_layer_model_digest(
    contract: ReferenceThroughBedK2Contract,
    models: tuple[ThroughBedKLayerEngineeringModel, ThroughBedKLayerEngineeringModel],
    *,
    integration_law_definition_digest: str,
) -> str:
    """Bind one reference contract and its ruled one-wall K1/K2 model pair."""

    _require_sha256("single-wall integration-law definition", integration_law_definition_digest)
    _validate_single_wall_layer_models(contract, models)
    digest = hashlib.sha256()
    for part in (
        (SINGLE_WALL_INTEGRATION_DIGEST_DOMAIN + "/layer-models").encode("ascii"),
        integration_law_definition_digest.encode("ascii"),
        contract.definition_digest.encode("ascii"),
        *(repr(model).encode("utf-8") for model in models),
    ):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    # B8: the lane band variant is repr-excluded on the shared model class;
    # bind it explicitly and only when declared (sealed digests unchanged).
    for index, model in enumerate(models):
        if model.declared_lane_drop_band_pa is not None:
            part = (f"lane-drop-band/{index}/{model.declared_lane_drop_band_pa!r}").encode("ascii")
            digest.update(len(part).to_bytes(8, "big"))
            digest.update(part)
    # C8 Tier 1: the sorbed-water arm is repr-excluded on the same shared model
    # class; bind it explicitly and only when declared (sealed digests
    # unchanged, and an armed model can never share a digest with a sealed one).
    for index, model in enumerate(models):
        if model.declared_sorbed_water_arm is not None:
            part = (f"sorbed-water-arm/{index}/{model.declared_sorbed_water_arm!r}").encode("ascii")
            digest.update(len(part).to_bytes(8, "big"))
            digest.update(part)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedThroughBedK2SingleWallTrayState:
    """Accepted K2 host and fast seeds with one scalar physical-wall state."""

    contract: ReferenceThroughBedK2Contract
    tray: tray_host.AcceptedKCellTrayState
    integration_law_definition_digest: str
    layer_model_configuration_digest: str
    cell_field_states: tuple[ef.BinaryNoInertCellFieldState, ef.BinaryNoInertCellFieldState]
    wall_temperature_k: float

    physically_qualifying: ClassVar[bool] = False
    industrially_feasible: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    complete_intertray_model: ClassVar[bool] = False
    complete_particle_thermodynamic_state_coupling: ClassVar[bool] = False
    reference_indirect_duty_is_record_only: ClassVar[bool] = True
    reference_indirect_duty_realized: ClassVar[bool] = False
    wall_parameters_declared_engineering_assumptions: ClassVar[bool] = True
    wall_parameters_identified: ClassVar[bool] = False
    physical_wall_count: ClassVar[int] = 1

    def __post_init__(self) -> None:
        if type(self.contract) is not ReferenceThroughBedK2Contract:
            raise ThroughBedK2SingleWallConfigurationError(
                "state contract has a foreign exact type"
            )
        if type(self.tray) is not tray_host.AcceptedKCellTrayState:
            raise ThroughBedK2SingleWallConfigurationError(
                "tray must be an accepted K-cell host state"
            )
        if self.tray.physical_tray_id != self.contract.physical_tray_id:
            raise ThroughBedK2SingleWallConfigurationError(
                "accepted tray and reference contract IDs disagree"
            )
        if self.contract != reference_through_bed_k2_contract(self.tray.physical_tray_id):
            raise ThroughBedK2SingleWallConfigurationError(
                "accepted tray contract is not canonical"
            )
        _require_sha256(
            "accepted single-wall integration-law definition",
            self.integration_law_definition_digest,
        )
        _require_sha256(
            "accepted single-wall layer-model configuration",
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
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall through-bed integration requires exactly K=2 and RTD=8"
            )
        if type(self.cell_field_states) is not tuple or len(self.cell_field_states) != 2:
            raise ThroughBedK2SingleWallConfigurationError(
                "accepted state requires two fast states"
            )
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
            raise ThroughBedK2SingleWallConfigurationError("accepted fast state has a foreign type")
        _require_finite(
            "accepted physical-wall temperature", self.wall_temperature_k, positive=True
        )

    @property
    def state_digest(self) -> str:
        digest = hashlib.sha256()
        for part in (
            SINGLE_WALL_INTEGRATION_DIGEST_DOMAIN.encode("ascii"),
            self.contract.definition_digest.encode("ascii"),
            self.tray.state_digest.encode("ascii"),
            self.integration_law_definition_digest.encode("ascii"),
            self.layer_model_configuration_digest.encode("ascii"),
            repr(self.cell_field_states).encode("utf-8"),
            struct.pack(">d", self.wall_temperature_k),
        ):
            digest.update(len(part).to_bytes(8, "big"))
            digest.update(part)
        return "sha256:" + digest.hexdigest()


def migrate_equal_legacy_wall_state(
    legacy_state: two_wall.AcceptedThroughBedK2Law2TrayState,
    *,
    layer_models: tuple[ThroughBedKLayerEngineeringModel, ThroughBedKLayerEngineeringModel],
    integration_law_authority: ThroughBedK2SingleWallIntegrationLawAuthority,
    expected_integration_law_definition_digest: str,
    expected_layer_model_configuration_digest: str,
) -> AcceptedThroughBedK2SingleWallTrayState:
    """Explicitly migrate only a legacy state whose two old walls are bit-identical.

    Evolved unequal wall temperatures cannot be assigned one physical history
    without inventing a mixing or averaging law, so they fail closed.
    """

    law_digest = _validate_expected_single_wall_integration_law(
        integration_law_authority,
        expected_integration_law_definition_digest,
    )
    if type(legacy_state) is not two_wall.AcceptedThroughBedK2Law2TrayState:
        raise ThroughBedK2SingleWallConfigurationError(
            "migration requires an exact legacy two-wall accepted state"
        )
    lower, upper = legacy_state.wall_temperatures_k
    if not _same_binary64(lower, upper):
        raise ThroughBedK2SingleWallConfigurationError(
            "unequal legacy wall temperatures cannot migrate to one physical wall"
        )
    model_digest = canonical_single_wall_layer_model_digest(
        legacy_state.contract,
        layer_models,
        integration_law_definition_digest=law_digest,
    )
    _require_sha256(
        "expected single-wall layer-model configuration",
        expected_layer_model_configuration_digest,
    )
    if model_digest != expected_layer_model_configuration_digest:
        raise ThroughBedK2SingleWallConfigurationError(
            "single-wall layer models differ from the caller-held expected configuration"
        )
    return AcceptedThroughBedK2SingleWallTrayState(
        contract=legacy_state.contract,
        tray=legacy_state.tray,
        integration_law_definition_digest=law_digest,
        layer_model_configuration_digest=model_digest,
        cell_field_states=legacy_state.cell_field_states,
        wall_temperature_k=lower,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedK2SingleWallStepRequest:
    attempt_id: str
    prior: AcceptedThroughBedK2SingleWallTrayState
    integration_law_definition_digest: str
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
    meal_mechanical_work_deposition: meal_work.QSCMealMechanicalWorkDepositionBundle | None = None
    #: B7-family memo (2026-08-31 rulings): the interval-preverified bundle
    #: digest; None keeps the sealed full-reconstruction behavior.
    preverified_meal_bundle_digest: str | None = None

    physically_qualifying: ClassVar[bool] = False
    industrially_feasible: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    intertray_gas_equality_proved: ClassVar[bool] = False
    shaft_work_energy_model_selected: ClassVar[bool] = False
    wall_parameters_declared_engineering_assumptions: ClassVar[bool] = True
    wall_parameters_identified: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.attempt_id) is not str or not self.attempt_id.strip():
            raise ThroughBedK2SingleWallConfigurationError("attempt ID must be nonblank")
        if type(self.prior) is not AcceptedThroughBedK2SingleWallTrayState:
            raise ThroughBedK2SingleWallConfigurationError("prior has a foreign exact type")
        AcceptedThroughBedK2SingleWallTrayState.__post_init__(self.prior)
        _require_sha256(
            "request single-wall integration-law definition",
            self.integration_law_definition_digest,
        )
        if self.integration_law_definition_digest != (self.prior.integration_law_definition_digest):
            raise ThroughBedK2SingleWallConfigurationError(
                "request and accepted state bind different single-wall laws"
            )
        _require_finite("end time", self.end_time_s)
        if self.end_time_s <= self.prior.tray.time_s:
            raise ThroughBedK2SingleWallConfigurationError("macro-step must advance accepted time")
        _validate_single_wall_layer_models(self.prior.contract, self.layer_models)
        if type(self.bottom_gas_boundary) is not ef.EngineeringGasBoundary:
            raise ThroughBedK2SingleWallConfigurationError("lower gas boundary has a foreign type")
        if self.bottom_gas_boundary.carrier_topology is not ef.CarrierTopology.BINARY_NO_INERT:
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall through-bed seam requires Law 2"
            )
        if type(self.bottom_gas_binding) is not ThroughBedGasBoundaryBinding:
            raise ThroughBedK2SingleWallConfigurationError(
                "gas-boundary binding has a foreign type"
            )
        if self.bottom_gas_binding.physical_tray_id != self.prior.tray.physical_tray_id:
            raise ThroughBedK2SingleWallConfigurationError(
                "gas-boundary binding belongs to another tray"
            )
        _require_finite("top boundary pressure", self.top_boundary_pressure_pa, positive=True)
        _require_finite("Newton tolerance", self.newton_tolerance, positive=True)
        _require_finite("face tolerance", self.face_tolerance, positive=True)
        _require_finite("relative limit", self.relative_limit, positive=True)
        if type(self.face_max_iterations) is not int or self.face_max_iterations <= 0:
            raise ThroughBedK2SingleWallConfigurationError(
                "face iteration limit must be a positive exact integer"
            )
        if type(self.flow_segments) is not tuple or not self.flow_segments:
            raise ThroughBedK2SingleWallConfigurationError("accepted-flow history must be nonempty")
        if any(
            type(item) is not tray_host.AcceptedDryMatterFlowSegment for item in self.flow_segments
        ):
            raise ThroughBedK2SingleWallConfigurationError(
                "flow history contains a foreign segment"
            )
        if self.flow_segments[0].start_time_s != self.prior.tray.time_s:
            raise ThroughBedK2SingleWallConfigurationError("flow history starts at another time")
        if self.flow_segments[-1].end_time_s != self.end_time_s:
            raise ThroughBedK2SingleWallConfigurationError("flow history ends at another time")
        observed_digest = canonical_single_wall_layer_model_digest(
            self.prior.contract,
            self.layer_models,
            integration_law_definition_digest=self.integration_law_definition_digest,
        )
        if observed_digest != self.prior.layer_model_configuration_digest:
            raise ThroughBedK2SingleWallConfigurationError(
                "single-wall layer-model configuration drifted from accepted state"
            )
        for model in self.layer_models:
            if model.properties.energy_datum_id != self.prior.tray.energy_datum_id:
                raise ThroughBedK2SingleWallConfigurationError(
                    "cell model and packet host use different energy datums"
                )
        if self.meal_mechanical_work_deposition is not None:
            bundle = self.meal_mechanical_work_deposition
            if type(bundle) is not meal_work.QSCMealMechanicalWorkDepositionBundle:
                raise ThroughBedK2SingleWallConfigurationError(
                    "meal mechanical-work deposition has a foreign exact type"
                )
            # B0 hoist: this request is rebuilt on EVERY Picard iteration and
            # its ``__post_init__`` runs three times per tray step, so the
            # bundle's full structural walk ran four times per step against a
            # bit-identical object (measured 5.6 ms each on a 96-packet tray,
            # 31 % of the closure).  When the caller threads the digest it
            # already derived from the SAME interval's full validation, the
            # walk collapses to the value-keyed re-derivation (0.22 ms) that
            # the B7 memo already uses one layer down: a mutated bundle
            # re-hashes and refuses.  This is the B9 precedent (one full
            # runtime validation per interval, digest pinned downstream).
            if HOIST_ENABLED and self.preverified_meal_bundle_digest is not None:
                _require_sha256(
                    "preverified meal bundle definition",
                    self.preverified_meal_bundle_digest,
                )
                if bundle._seal is not meal_work._QSC_MEAL_DEPOSITION_BUNDLE_SEAL:
                    raise ThroughBedK2SingleWallConfigurationError(
                        "meal mechanical-work deposition is not builder-issued"
                    )
                if bundle.definition_digest != self.preverified_meal_bundle_digest:
                    raise ThroughBedK2SingleWallConfigurationError(
                        "meal mechanical-work deposition differs from its "
                        "interval-preverified digest"
                    )
            else:
                meal_work.QSCMealMechanicalWorkDepositionBundle.__post_init__(bundle)
            if (
                bundle.prior_tray_state_digest != self.prior.tray.state_digest
                or bundle.physical_tray_id != self.prior.tray.physical_tray_id
                or bundle.t_n_s != self.prior.tray.time_s
                or bundle.end_time_s != self.end_time_s
                or bundle.energy_datum_id != self.prior.tray.energy_datum_id
            ):
                raise ThroughBedK2SingleWallConfigurationError(
                    "meal mechanical-work deposition differs from the local t_n tray step"
                )

    @property
    def definition_digest(self) -> str:
        return _digest(
            "through-bed-k2-single-wall-step-request",
            (
                self.attempt_id,
                self.prior.state_digest,
                self.integration_law_definition_digest,
                self.end_time_s,
                *(repr(model) for model in self.layer_models),
                repr(self.bottom_gas_boundary),
                repr(self.bottom_gas_binding),
                self.top_boundary_pressure_pa,
                *(repr(segment) for segment in self.flow_segments),
                self.newton_tolerance,
                self.face_tolerance,
                self.face_max_iterations,
                self.relative_limit,
                (
                    "NO_Q_OTHER_S"
                    if self.meal_mechanical_work_deposition is None
                    else self.meal_mechanical_work_deposition.definition_digest
                ),
            ),
        )


def _validate_request_law_and_model_pins(
    request: ThroughBedK2SingleWallStepRequest,
    *,
    integration_law_authority: ThroughBedK2SingleWallIntegrationLawAuthority,
    expected_integration_law_definition_digest: str,
    expected_layer_model_configuration_digest: str,
) -> tuple[str, str]:
    ThroughBedK2SingleWallStepRequest.__post_init__(request)
    law_digest = _validate_expected_single_wall_integration_law(
        integration_law_authority,
        expected_integration_law_definition_digest,
    )
    _require_sha256(
        "expected single-wall layer-model configuration",
        expected_layer_model_configuration_digest,
    )
    if request.integration_law_definition_digest != law_digest:
        raise ThroughBedK2SingleWallConfigurationError(
            "request binds a foreign caller-held single-wall law"
        )
    if request.prior.layer_model_configuration_digest != (
        expected_layer_model_configuration_digest
    ):
        raise ThroughBedK2SingleWallConfigurationError(
            "accepted state differs from the caller-held layer-model configuration"
        )
    observed_model_digest = canonical_single_wall_layer_model_digest(
        request.prior.contract,
        request.layer_models,
        integration_law_definition_digest=law_digest,
    )
    if observed_model_digest != expected_layer_model_configuration_digest:
        raise ThroughBedK2SingleWallConfigurationError(
            "request models differ from the caller-held layer-model configuration"
        )
    return law_digest, observed_model_digest


def _validate_meal_mechanical_work_deposition(
    request: ThroughBedK2SingleWallStepRequest,
    validation_context: meal_work.QSCMealMechanicalWorkValidationContext | None,
) -> meal_work.QSCMealMechanicalWorkDepositionBundle | None:
    bundle = request.meal_mechanical_work_deposition
    if bundle is None:
        if validation_context is not None:
            raise ThroughBedK2SingleWallConfigurationError(
                "a Q_other_s validation context is inadmissible without a deposition bundle"
            )
        return None
    if validation_context is None:
        raise ThroughBedK2SingleWallConfigurationError(
            "a meal mechanical-work deposition bundle requires caller-held validation pins"
        )
    try:
        meal_work.validate_qsc_meal_mechanical_work_deposition_bundle(
            bundle,
            accepted_tray=request.prior.tray,
            end_time_s=request.end_time_s,
            validation_context=validation_context,
            preverified_definition_digest=request.preverified_meal_bundle_digest,
        )
    except (TypeError, meal_work.QSCMealMechanicalWorkDepositionError) as error:
        raise ThroughBedK2SingleWallConfigurationError(
            f"meal mechanical-work deposition failed exact reconstruction: {error}"
        ) from error
    return bundle


@dataclass(frozen=True, slots=True, kw_only=True)
class PreparedThroughBedTrayContext:
    """B0: the Picard-invariant half of one through-bed tray step.

    Built once per attempt per tray by
    :func:`prepare_through_bed_k2_single_wall_tray_context` and threaded into
    every Picard iteration.  It carries no derived physics — only the
    verdicts and aggregates that are pure functions of objects that do not
    change while the face pressures iterate.  Every consumer re-verifies it
    against the live request by VALUE (tray state digest, bundle definition
    digest, inventory digest plus the aggregate identities the inventory's own
    construction guarantees), so a mutated context refuses instead of serving.
    """

    law_id: str
    physical_tray_id: str
    prior_tray_state_digest: str
    end_time_s: float
    law_digest: str
    model_digest: str
    codec_registry_digest: str
    auditor_identity_digest: str
    deposition_bundle: meal_work.QSCMealMechanicalWorkDepositionBundle | None
    deposition_bundle_digest: str | None
    inventories: tuple[PacketDerivedLayerInventory, PacketDerivedLayerInventory]
    inventory_digests: tuple[str, str]

    physically_qualifying: ClassVar[bool] = False


def _prepared_inventory_digest(inventory: PacketDerivedLayerInventory) -> str:
    parts: list[str | int | float | bool] = [
        inventory.layer,
        inventory.dry_matter_kg,
        inventory.attached_hexane_kg,
        inventory.internal_hexane_kg,
        inventory.residual_oil_label_kg,
        inventory.external_water_kg,
        inventory.common_datum_energy_j,
        inventory.solid_temperature_k,
        inventory.nil_film_reclassified_kg,
        inventory.solid_to_interface_ua_w_k,
        inventory.ua_share_sum_w_k,
        inventory.ua_partition_residual_w_k,
        inventory.ua_temperature_moment_residual_w,
        inventory.aggregation_law_id,
        inventory.energy_allocation_law_id,
        inventory.packet_count,
    ]
    # B1 stage 4 (2026-09-20): the external-water nil band's reclassified mass
    # binds this digest EXPLICITLY and CONDITIONALLY, exactly as the canonical
    # model digest binds the declared sorbed-water arm and the declared lane
    # drop band.  ``0.0`` - every certified run, and every arm-declared run
    # with no in-band layer - appends nothing and keeps the digest
    # byte-identical to the digest this function produced before the band
    # existed; a reclassifying layer carries the moved mass into the digest,
    # which is what makes the object pinned bit-for-bit on that step too.
    if inventory.nil_water_film_reclassified_kg != 0.0:
        parts.append(inventory.nil_water_film_reclassified_kg)
    for item in inventory.contributions:
        parts.append(item.packet_id)
        parts.append(item.representative_particles)
        parts.append(item.temperature_k)
    return _digest("through-bed-prepared-layer-inventory", parts)


def _prepared_inventory_self_consistent(inventory: PacketDerivedLayerInventory) -> bool:
    """The aggregate identities ``_layer_inventory`` established at build time.

    Every per-contribution field the digest above does not carry is an exact
    ``math.fsum`` summand of an aggregate field that it does carry, so the
    two together pin the whole object bit-for-bit.
    """

    contributions = inventory.contributions
    if len(contributions) != inventory.packet_count:
        return False
    attached = math.fsum(item.weighted_attached_hexane_kg for item in contributions)
    if inventory.nil_film_reclassified_kg != 0.0:
        if inventory.attached_hexane_kg != 0.0 or attached != inventory.nil_film_reclassified_kg:
            return False
    elif attached != inventory.attached_hexane_kg:
        return False
    # B1 stage 4 (2026-09-20): the external-water aggregate obeys the SAME
    # two-branch identity the attached film above obeys, for the same reason.
    # When the water nil band has fired, the view's external water is exactly
    # zero and the contributions sum to exactly the reclassified mass; when it
    # has not - every certified run - the pre-existing exact identity stands
    # unchanged.
    external = math.fsum(item.weighted_external_water_kg for item in contributions)
    if inventory.nil_water_film_reclassified_kg != 0.0:
        if (
            inventory.external_water_kg != 0.0
            or external != inventory.nil_water_film_reclassified_kg
        ):
            return False
    elif external != inventory.external_water_kg:
        return False
    moment = math.fsum(item.ua_share_w_k * item.temperature_k for item in contributions)
    return (
        math.fsum(item.weighted_dry_matter_kg for item in contributions)
        == inventory.dry_matter_kg
        and math.fsum(item.weighted_common_datum_energy_j for item in contributions)
        == inventory.common_datum_energy_j
        and math.fsum(item.ua_share_w_k for item in contributions) == inventory.ua_share_sum_w_k
        and (moment - inventory.solid_to_interface_ua_w_k * inventory.solid_temperature_k)
        == inventory.ua_temperature_moment_residual_w
    )


def _through_bed_layer_inventories(
    request: ThroughBedK2SingleWallStepRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
) -> tuple[PacketDerivedLayerInventory, PacketDerivedLayerInventory]:
    inventories = tuple(
        legacy._layer_inventory(
            request.prior.tray,
            layer,
            model=request.layer_models[layer - 1],
            codec=codec,
        )
        for layer in range(1, kernel.THROUGH_BED_LAYER_COUNT + 1)
    )
    if len(inventories) != 2:  # pragma: no cover - exact range guard
        raise AssertionError("unexpected through-bed inventory count")
    return (inventories[0], inventories[1])


def prepare_through_bed_k2_single_wall_tray_context(
    request: ThroughBedK2SingleWallStepRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner,
    integration_law_authority: ThroughBedK2SingleWallIntegrationLawAuthority,
    expected_integration_law_definition_digest: str,
    expected_layer_model_configuration_digest: str,
    meal_mechanical_work_validation_context: (
        meal_work.QSCMealMechanicalWorkValidationContext | None
    ) = None,
) -> PreparedThroughBedTrayContext:
    """Run the Picard-invariant half of one tray step exactly once."""

    if type(request) is not ThroughBedK2SingleWallStepRequest:
        raise TypeError("single-wall preparation requires an exact step request")
    law_digest, model_digest = _validate_request_law_and_model_pins(
        request,
        integration_law_authority=integration_law_authority,
        expected_integration_law_definition_digest=(expected_integration_law_definition_digest),
        expected_layer_model_configuration_digest=(expected_layer_model_configuration_digest),
    )
    _validate_codec_and_k_planner(request, codec, k_planner)
    deposition_bundle = _validate_meal_mechanical_work_deposition(
        request,
        meal_mechanical_work_validation_context,
    )
    inventories = _through_bed_layer_inventories(request, codec=codec)
    return PreparedThroughBedTrayContext(
        law_id=PREPARED_THROUGH_BED_CONTEXT_LAW_ID,
        physical_tray_id=request.prior.tray.physical_tray_id,
        prior_tray_state_digest=request.prior.tray.state_digest,
        end_time_s=request.end_time_s,
        law_digest=law_digest,
        model_digest=model_digest,
        codec_registry_digest=codec.codec_registry_digest,
        auditor_identity_digest=codec.auditor_identity_digest,
        deposition_bundle=deposition_bundle,
        deposition_bundle_digest=(
            None if deposition_bundle is None else deposition_bundle.definition_digest
        ),
        inventories=inventories,
        inventory_digests=(
            _prepared_inventory_digest(inventories[0]),
            _prepared_inventory_digest(inventories[1]),
        ),
    )


def _verify_prepared_through_bed_context(
    prepared: PreparedThroughBedTrayContext,
    request: ThroughBedK2SingleWallStepRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    expected_integration_law_definition_digest: str,
    expected_layer_model_configuration_digest: str,
) -> None:
    """Value-keyed re-hash: a mutated prepared context refuses, never serves."""

    if type(prepared) is not PreparedThroughBedTrayContext:
        raise TypeError("prepared through-bed context has a foreign exact type")
    if prepared.law_id != PREPARED_THROUGH_BED_CONTEXT_LAW_ID:
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context declares a foreign law"
        )
    if prepared.physical_tray_id != request.prior.tray.physical_tray_id:
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context belongs to another tray"
        )
    if prepared.prior_tray_state_digest != request.prior.tray.state_digest:
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context was prepared against another accepted tray"
        )
    if not _same_binary64(prepared.end_time_s, request.end_time_s):
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context was prepared for another endpoint"
        )
    if prepared.law_digest != expected_integration_law_definition_digest or (
        prepared.law_digest != request.integration_law_definition_digest
    ):
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context binds a foreign single-wall law"
        )
    if prepared.model_digest != expected_layer_model_configuration_digest or (
        prepared.model_digest != request.prior.layer_model_configuration_digest
    ):
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context binds a foreign layer-model configuration"
        )
    if prepared.codec_registry_digest != codec.codec_registry_digest or (
        prepared.auditor_identity_digest != codec.auditor_identity_digest
    ):
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context was prepared against another codec"
        )
    bundle = request.meal_mechanical_work_deposition
    if bundle is not prepared.deposition_bundle:
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context carries another deposition bundle"
        )
    if bundle is not None:
        if bundle.definition_digest != prepared.deposition_bundle_digest:
            raise ThroughBedK2SingleWallConfigurationError(
                "prepared through-bed deposition bundle re-hashes to another digest"
            )
    elif prepared.deposition_bundle_digest is not None:
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context declares a digest without a bundle"
        )
    if type(prepared.inventories) is not tuple or len(prepared.inventories) != 2:
        raise ThroughBedK2SingleWallConfigurationError(
            "prepared through-bed context has a foreign inventory pair"
        )
    for index, inventory in enumerate(prepared.inventories):
        if type(inventory) is not PacketDerivedLayerInventory:
            raise ThroughBedK2SingleWallConfigurationError(
                "prepared through-bed inventory has a foreign exact type"
            )
        if inventory.layer != index + 1:
            raise ThroughBedK2SingleWallConfigurationError(
                "prepared through-bed inventory is out of layer order"
            )
        if _prepared_inventory_digest(inventory) != prepared.inventory_digests[index]:
            raise ThroughBedK2SingleWallConfigurationError(
                "prepared through-bed inventory re-hashes to another digest"
            )
        if not _prepared_inventory_self_consistent(inventory):
            raise ThroughBedK2SingleWallConfigurationError(
                "prepared through-bed inventory lost its own aggregate identities"
            )
        if inventory.solid_to_interface_ua_w_k != (
            request.layer_models[index].transfer.solid_to_interface_ua_w_k
        ):
            raise ThroughBedK2SingleWallConfigurationError(
                "prepared through-bed inventory disagrees with the live layer model"
            )


def _validate_accepted_wall_pair(
    contract: ReferenceThroughBedK2Contract,
    models: tuple[ThroughBedKLayerEngineeringModel, ThroughBedKLayerEngineeringModel],
    accepted: tuple[ef.AcceptedEngineeringFastStep, ef.AcceptedEngineeringFastStep],
) -> None:
    _validate_single_wall_layer_models(contract, models)
    if (
        type(accepted) is not tuple
        or len(accepted) != 2
        or any(type(item) is not ef.AcceptedEngineeringFastStep for item in accepted)
    ):
        raise ThroughBedK2SingleWallStepError(
            "aggregate wall requires two exact accepted Law-2 fast steps"
        )
    lower, upper = accepted
    if lower.wall_parameters != models[0].wall or upper.wall_parameters != models[1].wall:
        raise ThroughBedK2SingleWallStepError(
            "accepted fast-step wall declarations differ from the ruled model pair"
        )
    if lower.transfer != models[0].transfer or upper.transfer != models[1].transfer:
        raise ThroughBedK2SingleWallStepError(
            "accepted fast-step transfer declarations differ from the ruled model pair"
        )
    if not _same_binary64(lower.old_wall_temperature_k, upper.old_wall_temperature_k):
        raise ThroughBedK2SingleWallStepError(
            "both K fast solves must see the bit-identical old physical-wall temperature"
        )
    if lower.macro_step_s != upper.macro_step_s:
        raise ThroughBedK2SingleWallStepError("both K fast solves must use one macro-step")
    if (lower.wall_storage_scale, upper.wall_storage_scale) != (1.0, 0.0):
        raise ThroughBedK2SingleWallStepError("exactly K1 may own the sole active wall capacity")
    if (upper.wall_to_gas_w, upper.wall_to_interface_w) != (0.0, 0.0):
        raise ThroughBedK2SingleWallStepError(
            "the inactive K2 wall placeholder cannot own an accepted wall transfer"
        )


def advance_aggregate_through_bed_wall_node(
    contract: ReferenceThroughBedK2Contract,
    models: tuple[ThroughBedKLayerEngineeringModel, ThroughBedKLayerEngineeringModel],
    accepted: tuple[ef.AcceptedEngineeringFastStep, ef.AcceptedEngineeringFastStep],
) -> ef.EngineeringWallNodeSolve:
    """Advance the sole wall once from both sealed shared-old-wall fast solves."""

    _validate_accepted_wall_pair(contract, models, accepted)
    physical_wall = models[0].wall
    lower, upper = accepted
    dt = lower.macro_step_s
    storage_conductance = physical_wall.thermal_capacity_j_k / dt
    denominator = (
        storage_conductance + physical_wall.steam_side_ua_w_k + physical_wall.ambient_ua_w_k
    )
    to_gas = math.fsum((lower.wall_to_gas_w, upper.wall_to_gas_w))
    to_interface = math.fsum((lower.wall_to_interface_w, upper.wall_to_interface_w))
    numerator = math.fsum(
        (
            storage_conductance * lower.old_wall_temperature_k,
            physical_wall.steam_side_ua_w_k * physical_wall.steam_temperature_k,
            physical_wall.ambient_ua_w_k * physical_wall.ambient_temperature_k,
            -to_gas,
            -to_interface,
        )
    )
    temperature = numerator / denominator
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ThroughBedK2SingleWallStepError(
            "the aggregate physical-wall update produced an invalid temperature"
        )
    for layer, step in enumerate(accepted, start=1):
        for receiver_name, conductance, receiver_temperature in (
            (
                "gas",
                step.transfer.wall_to_gas_ua_w_k,
                step.state.gas_temperature_k,
            ),
            (
                "interface",
                step.transfer.wall_to_interface_ua_w_k,
                step.state.interface_temperature_k,
            ),
        ):
            if conductance == 0.0:
                continue
            old_delta = step.old_wall_temperature_k - receiver_temperature
            new_delta = temperature - receiver_temperature
            reversed_order = (old_delta < 0.0 < new_delta) or (new_delta < 0.0 < old_delta)
            if reversed_order:
                raise ThroughBedK2SingleWallStepError(
                    "aggregate physical-wall update would reverse the shared-old-wall "
                    f"K{layer} wall-to-{receiver_name} temperature ordering; substep "
                    "required, never clamped"
                )
    steam_side = physical_wall.steam_side_ua_w_k * (physical_wall.steam_temperature_k - temperature)
    ambient = physical_wall.ambient_ua_w_k * (temperature - physical_wall.ambient_temperature_k)
    storage = storage_conductance * (temperature - lower.old_wall_temperature_k)
    residual = math.fsum((storage, -steam_side, ambient, to_gas, to_interface))
    return ef.EngineeringWallNodeSolve(
        wall_temperature_k=temperature,
        old_wall_temperature_k=lower.old_wall_temperature_k,
        steam_side_w=steam_side,
        to_ambient_w=ambient,
        shared_to_gas_w=to_gas,
        shared_to_interface_w=to_interface,
        storage_w=storage,
        residual_w=residual,
        row_derivative_w_k=denominator,
    )


def _aggregate_wall_row_resolution_bound(
    accepted: tuple[ef.AcceptedEngineeringFastStep, ef.AcceptedEngineeringFastStep],
    wall: ef.EngineeringWallNodeSolve,
) -> float:
    lower, upper = accepted
    parameters = lower.wall_parameters
    dt = lower.macro_step_s
    capacity = parameters.thermal_capacity_j_k
    storage_conductance = capacity / dt
    conductance_resolution = math.fsum(
        (
            abs(1.0 / dt) * legacy._binary64_resolution(capacity),
            abs(storage_conductance / dt) * legacy._binary64_resolution(dt),
            legacy._binary64_resolution(storage_conductance),
        )
    )
    temperature_delta = wall.wall_temperature_k - wall.old_wall_temperature_k
    temperature_delta_resolution = math.fsum(
        (
            legacy._binary64_resolution(wall.wall_temperature_k),
            legacy._binary64_resolution(wall.old_wall_temperature_k),
            legacy._binary64_resolution(temperature_delta),
        )
    )
    storage_resolution = math.fsum(
        (
            abs(storage_conductance) * temperature_delta_resolution,
            abs(temperature_delta) * conductance_resolution,
            legacy._binary64_resolution(wall.storage_w),
        )
    )
    steam_delta = parameters.steam_temperature_k - wall.wall_temperature_k
    steam_resolution = math.fsum(
        (
            abs(steam_delta) * legacy._binary64_resolution(parameters.steam_side_ua_w_k),
            abs(parameters.steam_side_ua_w_k)
            * math.fsum(
                (
                    legacy._binary64_resolution(parameters.steam_temperature_k),
                    legacy._binary64_resolution(wall.wall_temperature_k),
                )
            ),
            legacy._binary64_resolution(wall.steam_side_w),
        )
    )
    ambient_delta = wall.wall_temperature_k - parameters.ambient_temperature_k
    ambient_resolution = math.fsum(
        (
            abs(ambient_delta) * legacy._binary64_resolution(parameters.ambient_ua_w_k),
            abs(parameters.ambient_ua_w_k)
            * math.fsum(
                (
                    legacy._binary64_resolution(wall.wall_temperature_k),
                    legacy._binary64_resolution(parameters.ambient_temperature_k),
                )
            ),
            legacy._binary64_resolution(wall.to_ambient_w),
        )
    )
    transfer_resolution = math.fsum(
        legacy._binary64_resolution(value)
        for value in (
            lower.wall_to_gas_w,
            upper.wall_to_gas_w,
            wall.shared_to_gas_w,
            lower.wall_to_interface_w,
            upper.wall_to_interface_w,
            wall.shared_to_interface_w,
        )
    )
    row_operation_resolution = math.fsum(
        legacy._binary64_resolution(value)
        for value in (
            wall.storage_w,
            wall.steam_side_w,
            wall.to_ambient_w,
            wall.shared_to_gas_w,
            wall.shared_to_interface_w,
            wall.residual_w,
        )
    )
    return math.fsum(
        (
            storage_resolution,
            steam_resolution,
            ambient_resolution,
            transfer_resolution,
            row_operation_resolution,
        )
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class AggregatePhysicalWallLedger:
    physical_wall_count: int
    active_wall_capacity_count: int
    old_wall_temperature_bits_identical: bool
    inactive_upper_placeholder_advanced: bool
    lower_reference_indirect_duty_w: float
    upper_reference_indirect_duty_w: float
    lower_wall_storage_scale: float
    upper_wall_storage_scale: float
    lower_wall_to_gas_w: float
    upper_wall_to_gas_w: float
    accepted_wall_to_gas_sum_w: float
    wall_shared_to_gas_w: float
    wall_to_gas_cancellation_w: float
    lower_wall_to_interface_w: float
    upper_wall_to_interface_w: float
    accepted_wall_to_interface_sum_w: float
    wall_shared_to_interface_w: float
    wall_to_interface_cancellation_w: float
    wall_row_residual_w: float
    #: D15-a: the MULTIPLIED bound the ledger admits against; both
    #: ``passed`` properties re-derive the residual beyond resolution
    #: from THIS field.
    wall_row_binary64_resolution_bound_w: float
    #: D15-a: the bound as ``_aggregate_wall_row_resolution_bound``
    #: builds it, and the declared multiple applied to it, recorded so
    #: the ledger says what was admitted against what.
    wall_row_unmultiplied_resolution_bound_w: float
    wall_row_resolution_bound_multiple: float
    wall_row_residual_beyond_resolution_w: float
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    industrially_feasible: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    reference_indirect_duty_realized: ClassVar[bool] = False
    shaft_work_included: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        return (
            self.physical_wall_count == 1
            and self.active_wall_capacity_count == 1
            and self.old_wall_temperature_bits_identical
            and not self.inactive_upper_placeholder_advanced
            and self.upper_reference_indirect_duty_w == 0.0
            and (self.lower_wall_storage_scale, self.upper_wall_storage_scale) == (1.0, 0.0)
            and self.upper_wall_to_gas_w == 0.0
            and self.upper_wall_to_interface_w == 0.0
            and self.accepted_wall_to_gas_sum_w
            == math.fsum((self.lower_wall_to_gas_w, self.upper_wall_to_gas_w))
            and self.accepted_wall_to_interface_sum_w
            == math.fsum((self.lower_wall_to_interface_w, self.upper_wall_to_interface_w))
            and self.wall_to_gas_cancellation_w
            == self.accepted_wall_to_gas_sum_w - self.wall_shared_to_gas_w
            and self.wall_to_interface_cancellation_w
            == self.accepted_wall_to_interface_sum_w - self.wall_shared_to_interface_w
            and self.wall_to_gas_cancellation_w == 0.0
            and self.wall_to_interface_cancellation_w == 0.0
            and self.wall_row_residual_beyond_resolution_w
            == legacy._residual_beyond_resolution(
                self.wall_row_residual_w,
                self.wall_row_binary64_resolution_bound_w,
            )
            and self.maximum_relative_residual <= self.relative_limit
        )


def _aggregate_physical_wall_ledger(
    request: ThroughBedK2SingleWallStepRequest,
    accepted: tuple[ef.AcceptedEngineeringFastStep, ef.AcceptedEngineeringFastStep],
    wall: ef.EngineeringWallNodeSolve,
) -> AggregatePhysicalWallLedger:
    lower, upper = accepted
    to_gas_sum = math.fsum((lower.wall_to_gas_w, upper.wall_to_gas_w))
    to_interface_sum = math.fsum((lower.wall_to_interface_w, upper.wall_to_interface_w))
    gas_cancellation = to_gas_sum - wall.shared_to_gas_w
    interface_cancellation = to_interface_sum - wall.shared_to_interface_w
    unmultiplied_row_bound = _aggregate_wall_row_resolution_bound(accepted, wall)
    row_bound = unmultiplied_row_bound * WALL_ROW_RESOLUTION_BOUND_MULTIPLE
    row_beyond = legacy._residual_beyond_resolution(wall.residual_w, row_bound)
    maximum = max(
        legacy._relative_residual(gas_cancellation, to_gas_sum, wall.shared_to_gas_w),
        legacy._relative_residual(
            interface_cancellation,
            to_interface_sum,
            wall.shared_to_interface_w,
        ),
        legacy._relative_residual(
            row_beyond,
            wall.storage_w,
            wall.steam_side_w,
            wall.to_ambient_w,
            wall.shared_to_gas_w,
            wall.shared_to_interface_w,
        ),
    )
    return AggregatePhysicalWallLedger(
        physical_wall_count=1,
        active_wall_capacity_count=1,
        old_wall_temperature_bits_identical=_same_binary64(
            lower.old_wall_temperature_k,
            upper.old_wall_temperature_k,
        ),
        inactive_upper_placeholder_advanced=False,
        lower_reference_indirect_duty_w=(request.layer_models[0].geometry.tray.indirect_duty_w),
        upper_reference_indirect_duty_w=(request.layer_models[1].geometry.tray.indirect_duty_w),
        lower_wall_storage_scale=lower.wall_storage_scale,
        upper_wall_storage_scale=upper.wall_storage_scale,
        lower_wall_to_gas_w=lower.wall_to_gas_w,
        upper_wall_to_gas_w=upper.wall_to_gas_w,
        accepted_wall_to_gas_sum_w=to_gas_sum,
        wall_shared_to_gas_w=wall.shared_to_gas_w,
        wall_to_gas_cancellation_w=gas_cancellation,
        lower_wall_to_interface_w=lower.wall_to_interface_w,
        upper_wall_to_interface_w=upper.wall_to_interface_w,
        accepted_wall_to_interface_sum_w=to_interface_sum,
        wall_shared_to_interface_w=wall.shared_to_interface_w,
        wall_to_interface_cancellation_w=interface_cancellation,
        wall_row_residual_w=wall.residual_w,
        wall_row_binary64_resolution_bound_w=row_bound,
        wall_row_unmultiplied_resolution_bound_w=unmultiplied_row_bound,
        wall_row_resolution_bound_multiple=WALL_ROW_RESOLUTION_BOUND_MULTIPLE,
        wall_row_residual_beyond_resolution_w=row_beyond,
        maximum_relative_residual=maximum,
        relative_limit=request.relative_limit,
    )


@dataclass(frozen=True, slots=True)
class _LayerWallEnergyAllocation:
    steam_side_w: float
    to_ambient_w: float
    storage_w: float


@dataclass(frozen=True, slots=True)
class _LayerFastStepView:
    fast: ef.BinaryNoInertCellSolve
    accepted: ef.AcceptedEngineeringFastStep
    wall: ef.EngineeringWallNodeSolve | _LayerWallEnergyAllocation


def _layer_views(
    solves: tuple[ef.BinaryNoInertCellSolve, ef.BinaryNoInertCellSolve],
    accepted: tuple[ef.AcceptedEngineeringFastStep, ef.AcceptedEngineeringFastStep],
    wall: ef.EngineeringWallNodeSolve,
) -> tuple[_LayerFastStepView, _LayerFastStepView]:
    return (
        _LayerFastStepView(fast=solves[0], accepted=accepted[0], wall=wall),
        _LayerFastStepView(
            fast=solves[1],
            accepted=accepted[1],
            wall=_LayerWallEnergyAllocation(
                steam_side_w=0.0,
                to_ambient_w=0.0,
                storage_w=0.0,
            ),
        ),
    )


def _solve_serial_k_cells(
    request: ThroughBedK2SingleWallStepRequest,
    inventories: tuple[PacketDerivedLayerInventory, PacketDerivedLayerInventory],
) -> tuple[
    tuple[ef.EngineeringCellInputs, ef.EngineeringCellInputs],
    tuple[ef.BinaryNoInertCellSolve, ef.BinaryNoInertCellSolve],
    tuple[ef.AcceptedEngineeringFastStep, ef.AcceptedEngineeringFastStep],
    InterKGasFaceLedger,
]:
    macro_step_s = request.end_time_s - request.prior.tray.time_s
    face_pressure = request.bottom_gas_boundary.downstream_boundary_pressure_pa
    seeds = list(request.prior.cell_field_states)
    for iteration in range(1, request.face_max_iterations + 1):
        lower_boundary = replace(
            request.bottom_gas_boundary,
            downstream_boundary_pressure_pa=face_pressure,
        )
        lower_inputs = _tray_unit_of_account_cell_inputs(
            legacy._cell_inputs(
                model=request.layer_models[0],
                inventory=inventories[0],
                gas_boundary=lower_boundary,
                wall_temperature_k=request.prior.wall_temperature_k,
                macro_step_s=macro_step_s,
            )
        )
        lower_solve = legacy._lane_band_checked(
            request.layer_models[0],
            ef.solve_binary_no_inert_fast_block(
                lower_inputs,
                seed=seeds[0],
                tolerance=request.newton_tolerance,
            ),
        )
        seeds[0] = lower_solve.state
        upper_boundary = ef.EngineeringGasBoundary(
            inlet_molar_flow_mol_s=lower_solve.evaluation.outlet_molar_flow_mol_s,
            inlet_hexane_mole_fraction=lower_solve.state.gas_hexane_mole_fraction,
            inlet_water_mole_fraction=lower_solve.state.gas_water_mole_fraction,
            inlet_temperature_k=lower_solve.state.gas_temperature_k,
            downstream_boundary_pressure_pa=request.top_boundary_pressure_pa,
            carrier_topology=ef.CarrierTopology.BINARY_NO_INERT,
        )
        upper_inputs = _tray_unit_of_account_cell_inputs(
            legacy._cell_inputs(
                model=request.layer_models[1],
                inventory=inventories[1],
                gas_boundary=upper_boundary,
                wall_temperature_k=request.prior.wall_temperature_k,
                macro_step_s=macro_step_s,
            )
        )
        upper_solve = legacy._lane_band_checked(
            request.layer_models[1],
            ef.solve_binary_no_inert_fast_block(
                upper_inputs,
                seed=seeds[1],
                tolerance=request.newton_tolerance,
            ),
        )
        seeds[1] = upper_solve.state
        next_face_pressure = upper_solve.state.layer_pressure_pa
        pressure_residual = next_face_pressure - face_pressure
        pressure_relative = legacy._relative_residual(
            pressure_residual,
            next_face_pressure,
            face_pressure,
        )
        if pressure_relative <= request.face_tolerance:
            solves = (lower_solve, upper_solve)
            accepted = (
                ef.accept_engineering_fast_step(lower_inputs, lower_solve),
                ef.accept_engineering_fast_step(upper_inputs, upper_solve),
            )
            if not _same_binary64(
                accepted[0].old_wall_temperature_k,
                request.prior.wall_temperature_k,
            ) or not _same_binary64(
                accepted[1].old_wall_temperature_k,
                request.prior.wall_temperature_k,
            ):
                raise ThroughBedK2SingleWallStepError(
                    "a K fast solve did not retain the accepted scalar old wall"
                )
            face = InterKGasFaceLedger(
                lower_outlet_molar_flow_mol_s=(lower_solve.evaluation.outlet_molar_flow_mol_s),
                upper_inlet_molar_flow_mol_s=upper_boundary.inlet_molar_flow_mol_s,
                lower_outlet_hexane_mole_fraction=(lower_solve.state.gas_hexane_mole_fraction),
                upper_inlet_hexane_mole_fraction=(upper_boundary.inlet_hexane_mole_fraction),
                lower_outlet_temperature_k=lower_solve.state.gas_temperature_k,
                upper_inlet_temperature_k=upper_boundary.inlet_temperature_k,
                lower_downstream_pressure_pa=face_pressure,
                upper_layer_pressure_pa=next_face_pressure,
                pressure_residual_pa=pressure_residual,
                pressure_relative_residual=pressure_relative,
                pressure_relative_limit=request.face_tolerance,
                iterations=iteration,
            )
            return (lower_inputs, upper_inputs), solves, accepted, face
        face_pressure = next_face_pressure
    raise ThroughBedK2SingleWallStepError(
        "the serial K-cell internal pressure face did not converge without clipping"
    )


def _layer_cell_energy_resolution_bound(
    *,
    before: tray_host.AcceptedKCell,
    after: tray_host.AcceptedKCell,
    expected_storage_j: float,
    wall: ef.EngineeringWallNodeSolve,
    dt: float,
    transfer: LayerAcceptedTransfers,
    packet_ledger: LayerPacketExchangeLedger,
) -> float:
    before_energy = before.conserved_state.totals.common_datum_energy_j
    after_energy = after.conserved_state.totals.common_datum_energy_j
    delta_energy = after_energy - before_energy
    storage_resolution = legacy._binary64_resolution(expected_storage_j)
    if expected_storage_j != 0.0:
        storage_resolution = legacy._product_binary64_resolution(
            wall.storage_w,
            dt,
            expected_storage_j,
        )
    return math.fsum(
        (
            packet_ledger.energy_binary64_resolution_bound_j,
            legacy._binary64_resolution(before_energy),
            legacy._binary64_resolution(after_energy),
            legacy._binary64_resolution(delta_energy),
            storage_resolution,
            legacy._binary64_resolution(transfer.energy_from_packets_j),
            legacy._binary64_resolution(transfer.gas_and_wall_boundary_energy_to_cell_j),
        )
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class SingleWallCoupledTrayLedger:
    layer_hexane_mass_residuals_kg: tuple[float, float]
    layer_water_mass_residuals_kg: tuple[float, float]
    layer_expected_wall_storage_j: tuple[float, float]
    layer_cell_energy_residuals_j: tuple[float, float]
    layer_cell_energy_binary64_resolution_bounds_j: tuple[float, float]
    layer_cell_energy_residuals_beyond_resolution_j: tuple[float, float]
    physical_wall_host_layer: int
    physical_wall_storage_j: float
    physical_wall_row_residual_w: float
    physical_wall_row_binary64_resolution_bound_w: float
    #: D15-a: copied from the aggregate ledger beside the multiplied
    #: bound this ledger's ``passed`` re-derives from.
    physical_wall_row_unmultiplied_resolution_bound_w: float
    physical_wall_row_resolution_bound_multiple: float
    physical_wall_row_residual_beyond_resolution_w: float
    energy_resolution_law_id: str
    maximum_relative_residual: float
    relative_limit: float

    physically_qualifying: ClassVar[bool] = False
    industrially_feasible: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    physical_wall_count: ClassVar[int] = 1
    shaft_work_included: ClassVar[bool] = False

    @property
    def passed(self) -> bool:
        expected_energy = tuple(
            legacy._residual_beyond_resolution(residual, bound)
            for residual, bound in zip(
                self.layer_cell_energy_residuals_j,
                self.layer_cell_energy_binary64_resolution_bounds_j,
                strict=True,
            )
        )
        return (
            self.physical_wall_host_layer == 1
            and self.layer_expected_wall_storage_j == (self.physical_wall_storage_j, 0.0)
            and self.layer_cell_energy_residuals_beyond_resolution_j == expected_energy
            and self.physical_wall_row_residual_beyond_resolution_w
            == legacy._residual_beyond_resolution(
                self.physical_wall_row_residual_w,
                self.physical_wall_row_binary64_resolution_bound_w,
            )
            and self.energy_resolution_law_id == SINGLE_WALL_COUPLED_ENERGY_RESOLUTION_LAW_ID
            and self.maximum_relative_residual <= self.relative_limit
        )


def _single_wall_coupled_ledger(
    request: ThroughBedK2SingleWallStepRequest,
    host_step: tray_host.ValidatedKCellTrayStep,
    solves: tuple[ef.BinaryNoInertCellSolve, ef.BinaryNoInertCellSolve],
    accepted: tuple[ef.AcceptedEngineeringFastStep, ef.AcceptedEngineeringFastStep],
    wall: ef.EngineeringWallNodeSolve,
    wall_ledger: AggregatePhysicalWallLedger,
    transfers: tuple[LayerAcceptedTransfers, LayerAcceptedTransfers],
    whole_tray_energy: (
        WholeTrayExternalEnergyLedger
        | NativeWholeTrayExternalEnergyLedger
        | SingleWallQOtherSExternalEnergyLedger
    ),
    tray_pressure: AggregateThroughBedPressureLedger,
    packet_ledgers: tuple[LayerPacketExchangeLedger, LayerPacketExchangeLedger],
) -> SingleWallCoupledTrayLedger:
    dt = accepted[0].macro_step_s
    physical_storage_j = wall.storage_w * dt
    expected_storage = (physical_storage_j, 0.0)
    hexane_residuals: list[float] = []
    water_residuals: list[float] = []
    energy_residuals: list[float] = []
    energy_bounds: list[float] = []
    energy_beyond: list[float] = []
    relative: list[float] = [
        host_step.tray_ledger.maximum_relative_residual,
        whole_tray_energy.maximum_relative_residual,
        tray_pressure.maximum_relative_residual,
        wall_ledger.maximum_relative_residual,
        *(ledger.maximum_relative_residual for ledger in packet_ledgers),
    ]
    for index, (before, after, solve, transfer, packet_ledger, expected) in enumerate(
        zip(
            request.prior.tray.cells,
            host_step.candidate.cells,
            solves,
            transfers,
            packet_ledgers,
            expected_storage,
            strict=True,
        )
    ):
        delta_hexane = (
            after.conserved_state.totals.hexane_kg - before.conserved_state.totals.hexane_kg
        )
        delta_water = after.conserved_state.totals.water_kg - before.conserved_state.totals.water_kg
        delta_energy = (
            after.conserved_state.totals.common_datum_energy_j
            - before.conserved_state.totals.common_datum_energy_j
        )
        energy_residual = delta_energy - expected
        energy_bound = _layer_cell_energy_resolution_bound(
            before=before,
            after=after,
            expected_storage_j=expected,
            wall=wall,
            dt=dt,
            transfer=transfer,
            packet_ledger=packet_ledger,
        )
        beyond = legacy._residual_beyond_resolution(energy_residual, energy_bound)
        hexane_residuals.append(delta_hexane)
        water_residuals.append(delta_water)
        energy_residuals.append(energy_residual)
        energy_bounds.append(energy_bound)
        energy_beyond.append(beyond)
        relative.extend(
            (
                legacy._relative_residual(
                    delta_hexane,
                    transfer.hexane_from_packets_kg,
                    transfer.gas_boundary_hexane_to_cell_kg,
                ),
                legacy._relative_residual(
                    delta_water,
                    transfer.water_from_packets_kg,
                    transfer.gas_boundary_water_to_cell_kg,
                ),
                legacy._relative_residual(beyond, delta_energy, expected),
                solve.evaluation.scaled_residual_norm,
            )
        )
        if index not in (0, 1):  # pragma: no cover - exact K=2 zip guard
            raise AssertionError("unexpected through-bed K-layer index")
    return SingleWallCoupledTrayLedger(
        layer_hexane_mass_residuals_kg=tuple(hexane_residuals),
        layer_water_mass_residuals_kg=tuple(water_residuals),
        layer_expected_wall_storage_j=expected_storage,
        layer_cell_energy_residuals_j=tuple(energy_residuals),
        layer_cell_energy_binary64_resolution_bounds_j=tuple(energy_bounds),
        layer_cell_energy_residuals_beyond_resolution_j=tuple(energy_beyond),
        physical_wall_host_layer=1,
        physical_wall_storage_j=physical_storage_j,
        physical_wall_row_residual_w=wall_ledger.wall_row_residual_w,
        physical_wall_row_binary64_resolution_bound_w=(
            wall_ledger.wall_row_binary64_resolution_bound_w
        ),
        physical_wall_row_unmultiplied_resolution_bound_w=(
            wall_ledger.wall_row_unmultiplied_resolution_bound_w
        ),
        physical_wall_row_resolution_bound_multiple=(
            wall_ledger.wall_row_resolution_bound_multiple
        ),
        physical_wall_row_residual_beyond_resolution_w=(
            wall_ledger.wall_row_residual_beyond_resolution_w
        ),
        energy_resolution_law_id=SINGLE_WALL_COUPLED_ENERGY_RESOLUTION_LAW_ID,
        maximum_relative_residual=max(relative),
        relative_limit=request.relative_limit,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class ThroughBedK2SingleWallStepRejection:
    attempt_id: str
    reason: str
    rollback_state: AcceptedThroughBedK2SingleWallTrayState
    #: B1 (owner-signed): typed transition data when the refusal is a
    #: Law-2 layer phase exhaustion; None for every other refusal, so
    #: existing constructions and matchers are unchanged.
    phase_exhaustion: legacy.LayerPhaseExhaustionSignal | None = None

    physically_qualifying: ClassVar[bool] = False
    industrially_feasible: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


_VALIDATED_SINGLE_WALL_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class ValidatedThroughBedK2SingleWallTrayStep:
    request: ThroughBedK2SingleWallStepRequest
    candidate: AcceptedThroughBedK2SingleWallTrayState
    integration_law_definition_digest: str
    layer_model_configuration_digest: str
    request_definition_digest: str
    host_transaction: tray_host.ValidatedKCellTrayStep
    layer_inventories: tuple[PacketDerivedLayerInventory, PacketDerivedLayerInventory]
    cell_inputs: tuple[ef.EngineeringCellInputs, ef.EngineeringCellInputs]
    cell_fast_solves: tuple[ef.BinaryNoInertCellSolve, ef.BinaryNoInertCellSolve]
    accepted_fast_steps: tuple[
        ef.AcceptedEngineeringFastStep,
        ef.AcceptedEngineeringFastStep,
    ]
    physical_wall: ef.EngineeringWallNodeSolve
    physical_wall_ledger: AggregatePhysicalWallLedger
    layer_transfers: tuple[LayerAcceptedTransfers, LayerAcceptedTransfers]
    inter_k_face: InterKGasFaceLedger
    tray_pressure: AggregateThroughBedPressureLedger
    #: F-HULL-1: the printed hull's ruled unit of account, carried on the
    #: transaction beside the pressure ledger it is derived from.
    tray_printed_hull: ThroughBedTrayPrintedDropCorroboration
    whole_tray_energy: (
        WholeTrayExternalEnergyLedger
        | NativeWholeTrayExternalEnergyLedger
        | SingleWallQOtherSExternalEnergyLedger
    )
    packet_exchange_ledgers: tuple[LayerPacketExchangeLedger, LayerPacketExchangeLedger]
    coupled_ledger: SingleWallCoupledTrayLedger
    meal_mechanical_work_deposition_ledger: meal_work.QSCMealMechanicalWorkDepositionLedger | None
    meal_mechanical_work_validation_context_digest: str | None
    codec_registry_digest: str
    auditor_identity_digest: str
    auditor_source_identity: str
    _seal: object = field(repr=False, compare=False)

    physically_qualifying: ClassVar[bool] = False
    industrially_feasible: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False
    atomic_persistence_implemented: ClassVar[bool] = False
    complete_intertray_model: ClassVar[bool] = False
    reference_indirect_duty_realized: ClassVar[bool] = False
    shaft_work_energy_model_selected: ClassVar[bool] = False
    wall_parameters_declared_engineering_assumptions: ClassVar[bool] = True
    wall_parameters_identified: ClassVar[bool] = False
    physical_wall_count: ClassVar[int] = 1

    def __post_init__(self) -> None:
        if self._seal is not _VALIDATED_SINGLE_WALL_SEAL:
            raise TypeError("validated single-wall steps are issued only by the evaluator")
        for name, value in (
            ("transaction integration-law definition", self.integration_law_definition_digest),
            ("transaction layer-model configuration", self.layer_model_configuration_digest),
            ("transaction request definition", self.request_definition_digest),
        ):
            _require_sha256(name, value)
        if type(self.request) is not ThroughBedK2SingleWallStepRequest:
            raise TypeError("validated single-wall request has a foreign exact type")
        if type(self.candidate) is not AcceptedThroughBedK2SingleWallTrayState:
            raise TypeError("validated single-wall candidate has a foreign exact type")
        ThroughBedK2SingleWallStepRequest.__post_init__(self.request)
        AcceptedThroughBedK2SingleWallTrayState.__post_init__(self.candidate)
        if (
            self.integration_law_definition_digest != self.request.integration_law_definition_digest
            or self.integration_law_definition_digest
            != self.candidate.integration_law_definition_digest
        ):
            raise ThroughBedK2SingleWallStepError(
                "validated transaction binds inconsistent single-wall laws"
            )
        if (
            self.layer_model_configuration_digest
            != self.request.prior.layer_model_configuration_digest
            or self.layer_model_configuration_digest
            != self.candidate.layer_model_configuration_digest
        ):
            raise ThroughBedK2SingleWallStepError(
                "validated transaction binds inconsistent layer-model configurations"
            )
        if self.request_definition_digest != self.request.definition_digest:
            raise ThroughBedK2SingleWallStepError(
                "validated transaction binds a foreign request definition"
            )
        if type(self.tray_printed_hull) is not ThroughBedTrayPrintedDropCorroboration:
            raise TypeError("the tray printed-hull row has a foreign exact type")
        if (
            self.tray_printed_hull.printed_hull_pa is not PRINTED_TRAY_DROP_HULL_PA
            or self.tray_printed_hull.law_id != PRINTED_TRAY_DROP_HULL_UNIT_OF_ACCOUNT_LAW_ID
        ):
            raise ThroughBedK2SingleWallStepError(
                "the tray printed-hull row carries a foreign band or law id"
            )
        if not self.tray_printed_hull.inside_printed_hull or not _same_binary64(
            self.tray_printed_hull.tray_drop_sum_pa,
            self.tray_pressure.layer_drop_sum_pa,
        ):
            raise ThroughBedK2SingleWallStepError(
                "validated transaction carries a tray drop that is outside the "
                "printed hull or disagrees with its own pressure ledger"
            )
        bundle = self.request.meal_mechanical_work_deposition
        if bundle is None:
            if (
                self.meal_mechanical_work_deposition_ledger is not None
                or self.meal_mechanical_work_validation_context_digest is not None
                or (
                    type(self.whole_tray_energy) is not WholeTrayExternalEnergyLedger
                    and type(self.whole_tray_energy) is not NativeWholeTrayExternalEnergyLedger
                )
            ):
                raise ThroughBedK2SingleWallStepError(
                    "zero-work transaction contains a mechanical-work artifact"
                )
        else:
            if (
                self.meal_mechanical_work_deposition_ledger != bundle.ledger
                or type(self.meal_mechanical_work_validation_context_digest) is not str
                or type(self.whole_tray_energy) is not SingleWallQOtherSExternalEnergyLedger
                or self.whole_tray_energy.deposition_ledger != bundle.ledger
            ):
                raise ThroughBedK2SingleWallStepError(
                    "work-bearing transaction and Q_other_s authorities disagree"
                )
            _require_sha256(
                "transaction meal-work validation context",
                self.meal_mechanical_work_validation_context_digest,
            )


def _validate_codec_and_k_planner(
    request: ThroughBedK2SingleWallStepRequest,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner,
) -> None:
    if type(codec) is not high_loading.HighLoadingV1PacketAdapter and (
        type(codec) is not dual_variant.EngineeringDualVariantPacketAdapter
    ):
        raise ThroughBedK2SingleWallStepError(
            "single-wall integration requires the exact detached tag-1 adapter "
            "or the owner-ruled dual-variant adapter"
        )
    planner_type = type(k_planner)
    if planner_type not in (
        tray_host.NoKRekeyPlanner,
        k_scheduler.EngineeringPhysicalKRekeyPlanner,
    ):
        raise ThroughBedK2SingleWallStepError(
            "the fixed-pressure tag-1 slice admits only the exact zero-motion or "
            "replay-validated engineering K-event planner"
        )
    if planner_type is tray_host.NoKRekeyPlanner:
        if (
            k_planner.planner_id != "manufactured-zero-k-rekey-v1"
            or k_planner.source_identity != "declared-zero-motion"
        ):
            raise ThroughBedK2SingleWallStepError(
                "zero-motion planner ID or source differs from the canonical NoK declaration"
            )
    else:
        try:
            k_scheduler.EngineeringPhysicalKRekeyPlanner.__post_init__(k_planner)
        except (TypeError, k_scheduler.EngineeringKEventSchedulerError) as error:
            raise ThroughBedK2SingleWallStepError(
                f"engineering K planner failed deterministic replay validation: {error}"
            ) from error
        schedule_step = k_planner.step
        if schedule_step.prior.time_s != request.prior.tray.time_s:
            raise ThroughBedK2SingleWallStepError(
                "engineering K planner prior time differs from accepted tray"
            )
        if schedule_step.candidate.time_s != request.end_time_s:
            raise ThroughBedK2SingleWallStepError(
                "engineering K planner endpoint differs from tray endpoint"
            )
        if schedule_step.prior.physical_tray_id != request.prior.tray.physical_tray_id:
            raise ThroughBedK2SingleWallStepError("engineering K planner belongs to another tray")
    for name, value in (
        ("codec registry", getattr(codec, "codec_registry_digest", None)),
        ("auditor identity", getattr(codec, "auditor_identity_digest", None)),
        ("auditor source", getattr(codec, "source_identity", None)),
    ):
        if type(value) is not str or not value.strip():
            raise ThroughBedK2SingleWallStepError(f"{name} is missing or blank")
    if codec.codec_registry_digest != request.prior.tray.codec_registry_digest:
        raise ThroughBedK2SingleWallStepError("codec violates host registry pin")
    if codec.auditor_identity_digest != request.prior.tray.auditor_identity_digest:
        raise ThroughBedK2SingleWallStepError("codec violates host auditor pin")
    codec.require_host_binding(
        configuration_digest=request.prior.tray.configuration_digest,
        energy_datum_id=request.prior.tray.energy_datum_id,
        pressure_pa=request.prior.tray.pressure_pa,
    )
    # CELL-02d W2: the declared arm keeps the pinned-adapter drift gate; the
    # native arm's coherence is the datum-identity triple (its own post-init
    # pin, the host gates above, and the codec datum here).
    cell_properties = request.layer_models[0].properties
    if type(cell_properties) is cc.DeclaredGasProperties:
        codec.component_datum_adapter.require_gas_properties(cell_properties)
    elif cell_properties.energy_datum_id != codec.energy_datum_id:
        raise ThroughBedK2SingleWallStepError(
            "native cell properties and codec use different energy datums"
        )


def _packet_entries_after_with_q_other_s(
    request: ThroughBedK2SingleWallStepRequest,
    transfers: tuple[LayerAcceptedTransfers, LayerAcceptedTransfers],
    inventories: tuple[PacketDerivedLayerInventory, PacketDerivedLayerInventory],
    views: tuple[ef.EngineeringCellMacroStep, ef.EngineeringCellMacroStep],
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    bundle: meal_work.QSCMealMechanicalWorkDepositionBundle | None,
) -> tuple[
    dict[str, tray_host.PacketPlanEntry],
    tuple[LayerPacketExchangeLedger, LayerPacketExchangeLedger],
]:
    """Solve every candidate from t_n with thermal transfer plus Q_other_s.

    The legacy solve supplies only the already-accepted mass-transfer and
    thermal-energy targets.  A work-bearing candidate is then solved again
    from the original t_n payload with the combined energy target; it is never
    formed by adding energy to, or otherwise patching, the thermal candidate.
    """

    thermal_entries, thermal_ledgers = legacy._packet_entries_after(
        request,
        transfers,
        inventories,
        views,
        codec=codec,
    )
    if bundle is None:
        return thermal_entries, thermal_ledgers

    entries: dict[str, tray_host.PacketPlanEntry] = {}
    rows_by_layer: dict[int, list[tuple[float, float, float, float]]] = {1: [], 2: []}
    for packet in request.prior.tray.packets:
        layer = packet.owner_key.vertical_layer_id.value
        deposition = bundle.packet_row(packet.packet_id)
        thermal = thermal_entries[packet.packet_id]
        total_hexane_after = math.fsum(
            (
                thermal.inventory.attached_hexane_kg,
                thermal.inventory.internal_hexane_kg,
            )
        )
        combined_energy = math.fsum(
            (
                thermal.inventory.common_datum_energy_j,
                deposition.mechanical_work_per_particle_j,
            )
        )
        # C8 Tier 1: the work-bearing candidate is re-solved from the ORIGINAL
        # payload, so it must carry the thermal candidate's SORBED water too or
        # the sorbed-arm sink would be silently dropped here.  Passed only when
        # the sorbed water actually moved, so every film-arm call reaches the
        # codec with its pre-2026-09-05 argument list exactly.
        retained_water_keyword = (
            {}
            if thermal.inventory.retained_water_kg == packet.entry.inventory.retained_water_kg
            else {"retained_water_mass_after_kg": thermal.inventory.retained_water_kg}
        )
        complete = codec.advance(
            packet.entry.payload_bytes,
            total_hexane_mass_after_kg=total_hexane_after,
            external_free_water_mass_after_kg=thermal.inventory.external_water_kg,
            envelope_energy_after_j=combined_energy,
            **retained_water_keyword,
        )
        representatives = packet.entry.weight.representative_particles
        before = packet.entry.inventory
        after = complete.inventory
        realized_hexane = representatives * math.fsum(
            (
                before.attached_hexane_kg,
                before.internal_hexane_kg,
                -after.attached_hexane_kg,
                -after.internal_hexane_kg,
            )
        )
        # C8 Tier 1: the realized water draw is the packet's TOTAL water change
        # (external plus sorbed), matching the legacy thermal ledger's own row.
        realized_water = representatives * math.fsum(
            (
                before.external_water_kg,
                before.retained_water_kg,
                -after.external_water_kg,
                -after.retained_water_kg,
            )
        )
        # Remove the external work before interpreting the remaining packet
        # energy decrement as packet-to-cell thermal transfer.
        realized_thermal_energy = representatives * math.fsum(
            (
                before.common_datum_energy_j,
                deposition.mechanical_work_per_particle_j,
                -after.common_datum_energy_j,
            )
        )
        resolution = meal_work.packet_q_other_s_roundoff_bound_j(
            deposition,
            per_particle_energy_before_j=before.common_datum_energy_j,
            per_particle_energy_after_j=after.common_datum_energy_j,
        )
        rows_by_layer[layer].append(
            (realized_hexane, realized_water, realized_thermal_energy, resolution)
        )
        entries[packet.packet_id] = tray_host.PacketPlanEntry(
            packet_id=packet.packet_id,
            payload_bytes=complete.payload_bytes,
            inventory=complete.inventory,
            weight=packet.entry.weight,
            surface=complete.surface,
        )

    ledgers: list[LayerPacketExchangeLedger] = []
    for layer, thermal_ledger in zip((1, 2), thermal_ledgers, strict=True):
        rows = rows_by_layer[layer]
        realized_hexane = math.fsum(row[0] for row in rows)
        realized_water = math.fsum(row[1] for row in rows)
        realized_energy = math.fsum(row[2] for row in rows)
        energy_residual = realized_energy - thermal_ledger.requested_energy_transfer_j
        energy_bound = math.fsum(
            (
                thermal_ledger.energy_binary64_resolution_bound_j,
                *(row[3] for row in rows),
                legacy._binary64_resolution(realized_energy),
                legacy._binary64_resolution(energy_residual),
            )
        )
        energy_beyond = legacy._residual_beyond_resolution(energy_residual, energy_bound)
        maximum = max(
            thermal_ledger.maximum_relative_residual,
            legacy._relative_residual(
                realized_hexane - thermal_ledger.requested_hexane_transfer_kg,
                realized_hexane,
                thermal_ledger.requested_hexane_transfer_kg,
            ),
            legacy._relative_residual(
                realized_water - thermal_ledger.requested_water_transfer_kg,
                realized_water,
                thermal_ledger.requested_water_transfer_kg,
            ),
            legacy._relative_residual(
                energy_beyond,
                realized_energy,
                thermal_ledger.requested_energy_transfer_j,
            ),
        )
        ledgers.append(
            replace(
                thermal_ledger,
                realized_hexane_transfer_kg=realized_hexane,
                realized_water_transfer_kg=realized_water,
                realized_energy_transfer_j=realized_energy,
                energy_transfer_residual_j=energy_residual,
                energy_binary64_resolution_bound_j=energy_bound,
                energy_transfer_residual_beyond_resolution_j=energy_beyond,
                maximum_relative_residual=maximum,
            )
        )
    return entries, (ledgers[0], ledgers[1])


@dataclass(frozen=True, slots=True, kw_only=True)
class SingleWallQOtherSExternalEnergyLedger:
    """Whole-tray external-energy closure with meal mechanical work explicit."""

    thermal_ledger: WholeTrayExternalEnergyLedger | NativeWholeTrayExternalEnergyLedger
    deposition_ledger: meal_work.QSCMealMechanicalWorkDepositionLedger
    thermal_external_energy_j: float
    q_other_s_external_energy_j: float
    total_external_energy_j: float
    whole_tray_system_energy_change_j: float
    common_energy_residual_j: float
    binary64_resolution_bound_j: float
    residual_beyond_resolution_j: float
    energy_resolution_law_id: str
    maximum_relative_residual: float
    relative_limit: float
    #: Q-F2a (owner ruling 2026-09-02, commit 20f5028).  On the DECLARED
    #: arm the internal gas face cancels exactly and this is positive zero,
    #: so every pre-flip ledger reconstructs bit-identically.  On the
    #: NATIVE arm the verbatim face hand-off books a throttling defect
    #: (CELL-02d W5) which is a real energy input to the tray, so it must
    #: appear on the EXTERNAL side of the whole-tray closure instead of
    #: being smoothed into the residual.
    booked_face_throttling_defect_energy_j: float = 0.0

    source_temperature_defined: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_wiring_added: ClassVar[bool] = False

    @property
    def face_cancellation_is_booked(self) -> bool:
        """The internal gas face is either exactly closed or exactly booked."""

        if type(self.thermal_ledger) is NativeWholeTrayExternalEnergyLedger:
            return (
                self.thermal_ledger.internal_face_cancellation_w
                == self.thermal_ledger.declared_face_throttling_defect_w
            )
        return self.thermal_ledger.internal_face_cancellation_w == 0.0

    @property
    def passed(self) -> bool:
        return (
            self.deposition_ledger.passed
            and self.face_cancellation_is_booked
            and self.thermal_external_energy_j == self.thermal_ledger.external_only_tray_energy_j
            and self.q_other_s_external_energy_j
            == self.deposition_ledger.emitted_packet_mechanical_work_j
            and self.total_external_energy_j
            == math.fsum(
                (
                    self.thermal_external_energy_j,
                    self.q_other_s_external_energy_j,
                    self.booked_face_throttling_defect_energy_j,
                )
            )
            and self.whole_tray_system_energy_change_j
            == self.thermal_ledger.whole_tray_system_energy_change_j
            and self.common_energy_residual_j
            == self.whole_tray_system_energy_change_j - self.total_external_energy_j
            and self.residual_beyond_resolution_j
            == legacy._residual_beyond_resolution(
                self.common_energy_residual_j,
                self.binary64_resolution_bound_j,
            )
            and self.energy_resolution_law_id
            == SINGLE_WALL_Q_OTHER_S_EXTERNAL_ENERGY_RESOLUTION_LAW_ID
            and self.maximum_relative_residual <= self.relative_limit
        )


def _q_other_s_whole_tray_external_energy_ledger(
    thermal_ledger: WholeTrayExternalEnergyLedger | NativeWholeTrayExternalEnergyLedger,
    host_step: tray_host.ValidatedKCellTrayStep,
    bundle: meal_work.QSCMealMechanicalWorkDepositionBundle,
    *,
    relative_limit: float,
    macro_step_s: float,
) -> SingleWallQOtherSExternalEnergyLedger:
    q_other_s = bundle.ledger.emitted_packet_mechanical_work_j
    thermal_external = thermal_ledger.external_only_tray_energy_j
    # Q-F2a: the native arm's booked internal-face throttling defect enters
    # the whole-tray closure as external energy, reconstructed from the
    # thermal ledger's own recorded defect POWER by the same ``dt * W``
    # operation ``sp1_k2_law2_tray_integration`` used to book it, so the
    # two are bit-identical.  Exact positive zero on the declared arm.
    native_arm = type(thermal_ledger) is NativeWholeTrayExternalEnergyLedger
    face_defect_energy = (
        macro_step_s * thermal_ledger.declared_face_throttling_defect_w if native_arm else 0.0
    )
    total_external = math.fsum((thermal_external, q_other_s, face_defect_energy))
    system_change = thermal_ledger.whole_tray_system_energy_change_j
    residual = system_change - total_external
    resolution = math.fsum(
        (
            thermal_ledger.whole_tray_energy_binary64_resolution_bound_j,
            bundle.ledger.binary64_roundoff_bound_j,
            legacy._binary64_resolution(q_other_s),
            legacy._binary64_resolution(total_external),
            legacy._binary64_resolution(residual),
        )
    )
    beyond = legacy._residual_beyond_resolution(residual, resolution)
    internal_face_residual = (
        thermal_ledger.defect_corrected_internal_face_energy_residual_j
        if native_arm
        else thermal_ledger.internal_face_energy_residual_j
    )
    face_cancellation_residual = (
        thermal_ledger.internal_face_cancellation_w
        - thermal_ledger.declared_face_throttling_defect_w
        if native_arm
        else thermal_ledger.internal_face_cancellation_w
    )
    maximum = max(
        host_step.tray_ledger.maximum_relative_residual,
        legacy._relative_residual(
            face_cancellation_residual,
            thermal_ledger.internal_face_inlet_energy_w,
            thermal_ledger.internal_face_outlet_energy_w,
        ),
        legacy._relative_residual(
            internal_face_residual,
            thermal_ledger.per_cell_external_energy_j,
            thermal_external,
        ),
        legacy._relative_residual(beyond, system_change, total_external),
        legacy._relative_residual(
            bundle.ledger.binary64_closure_residual_j,
            bundle.ledger.tray_declared_mechanical_work_j,
            q_other_s,
        ),
    )
    return SingleWallQOtherSExternalEnergyLedger(
        thermal_ledger=thermal_ledger,
        deposition_ledger=bundle.ledger,
        thermal_external_energy_j=thermal_external,
        q_other_s_external_energy_j=q_other_s,
        total_external_energy_j=total_external,
        whole_tray_system_energy_change_j=system_change,
        common_energy_residual_j=residual,
        binary64_resolution_bound_j=resolution,
        residual_beyond_resolution_j=beyond,
        energy_resolution_law_id=(SINGLE_WALL_Q_OTHER_S_EXTERNAL_ENERGY_RESOLUTION_LAW_ID),
        maximum_relative_residual=maximum,
        relative_limit=relative_limit,
        booked_face_throttling_defect_energy_j=face_defect_energy,
    )


def evaluate_through_bed_k2_single_wall_tray_step(
    request: ThroughBedK2SingleWallStepRequest,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    k_planner: tray_host.ManufacturedKRekeyPlanner,
    integration_law_authority: ThroughBedK2SingleWallIntegrationLawAuthority,
    expected_integration_law_definition_digest: str,
    expected_layer_model_configuration_digest: str,
    meal_mechanical_work_validation_context: (
        meal_work.QSCMealMechanicalWorkValidationContext | None
    ) = None,
    prepared_cache: dict[str, PreparedThroughBedTrayContext] | None = None,
) -> ValidatedThroughBedK2SingleWallTrayStep | ThroughBedK2SingleWallStepRejection:
    """Evaluate one MN1, MN2, or SP1 step without mutating accepted state.

    ``prepared_cache`` is an OPTIONAL caller-owned, per-attempt dictionary
    keyed by physical tray ID (B0 hoist).  The first call for a tray fills it
    from the ordinary un-hoisted path; later calls in the same attempt — the
    remaining pressure Picard iterations — re-verify the record by value and
    reuse it.  Passing ``None`` (the default) keeps every prior caller
    byte-identical, and so does ``DTDC_DISABLE_HOIST=1``.
    """

    if type(request) is not ThroughBedK2SingleWallStepRequest:
        raise TypeError("single-wall integration requires an exact step request")
    hoisting = prepared_cache is not None and HOIST_ENABLED
    if hoisting and type(prepared_cache) is not dict:
        raise TypeError("prepared through-bed cache must be an exact dict")
    try:
        prepared = (
            prepared_cache.get(request.prior.tray.physical_tray_id) if hoisting else None
        )
        if prepared is not None:
            # B0 hoist: the law/model pins, the codec/K-planner validation, the
            # deposition-bundle reconstruction, and both layer inventories are
            # pure functions of objects that do not move while the face
            # pressures iterate.  Re-verify the prepared record BY VALUE and
            # reuse it; the un-hoisted branch below stays the reference.
            _verify_prepared_through_bed_context(
                prepared,
                request,
                codec=codec,
                expected_integration_law_definition_digest=(
                    _validate_expected_single_wall_integration_law(
                        integration_law_authority,
                        expected_integration_law_definition_digest,
                    )
                ),
                expected_layer_model_configuration_digest=(
                    expected_layer_model_configuration_digest
                ),
            )
            law_digest = prepared.law_digest
            model_digest = prepared.model_digest
            deposition_bundle = prepared.deposition_bundle
            typed_inventories = prepared.inventories
        else:
            law_digest, model_digest = _validate_request_law_and_model_pins(
                request,
                integration_law_authority=integration_law_authority,
                expected_integration_law_definition_digest=(
                    expected_integration_law_definition_digest
                ),
                expected_layer_model_configuration_digest=(
                    expected_layer_model_configuration_digest
                ),
            )
            _validate_codec_and_k_planner(request, codec, k_planner)
            deposition_bundle = _validate_meal_mechanical_work_deposition(
                request,
                meal_mechanical_work_validation_context,
            )
            typed_inventories = _through_bed_layer_inventories(request, codec=codec)
            if hoisting:
                prepared_cache[request.prior.tray.physical_tray_id] = (
                    PreparedThroughBedTrayContext(
                        law_id=PREPARED_THROUGH_BED_CONTEXT_LAW_ID,
                        physical_tray_id=request.prior.tray.physical_tray_id,
                        prior_tray_state_digest=request.prior.tray.state_digest,
                        end_time_s=request.end_time_s,
                        law_digest=law_digest,
                        model_digest=model_digest,
                        codec_registry_digest=codec.codec_registry_digest,
                        auditor_identity_digest=codec.auditor_identity_digest,
                        deposition_bundle=deposition_bundle,
                        deposition_bundle_digest=(
                            None
                            if deposition_bundle is None
                            else deposition_bundle.definition_digest
                        ),
                        inventories=typed_inventories,
                        inventory_digests=(
                            _prepared_inventory_digest(typed_inventories[0]),
                            _prepared_inventory_digest(typed_inventories[1]),
                        ),
                    )
                )
        inputs, solves, accepted, face = _solve_serial_k_cells(
            request,
            typed_inventories,
        )
        if not face.passed:
            raise ThroughBedK2SingleWallStepError("internal K gas-face ledger did not pass")
        physical_wall = advance_aggregate_through_bed_wall_node(
            request.prior.contract,
            request.layer_models,
            accepted,
        )
        wall_ledger = _aggregate_physical_wall_ledger(
            request,
            accepted,
            physical_wall,
        )
        if not wall_ledger.passed:
            raise ThroughBedK2SingleWallStepError("aggregate physical-wall ledger did not pass")
        views = _layer_views(solves, accepted, physical_wall)
        tray_pressure = legacy._aggregate_pressure_ledger(request, inputs, views)
        if not tray_pressure.passed:
            raise ThroughBedK2SingleWallStepError(
                "aggregate one-floor pressure ledger did not pass"
            )
        # F-HULL-1: the printed hull, on the TRAY.  Both converged layer drops
        # are known here and nowhere earlier, which is why the ruled gate lives
        # in this seam and not inside the pinned per-cell solve.
        tray_printed_hull = require_corroborated_tray_drop(
            tuple(solve.evaluation.layer_pressure_drop_pa for solve in solves)
        )
        if not _same_binary64(
            tray_printed_hull.tray_drop_sum_pa,
            tray_pressure.layer_drop_sum_pa,
        ):
            raise ThroughBedK2SingleWallStepError(
                "the per-tray printed-hull sum and the aggregate pressure ledger's "
                "own layer-drop sum are not the same binary64"
            )
        transfers = tuple(
            legacy._accepted_layer_transfers(
                inputs[layer - 1],
                views[layer - 1],
                layer=layer,
                datum_adapter=codec.component_datum_adapter,
            )
            for layer in range(1, kernel.THROUGH_BED_LAYER_COUNT + 1)
        )
        typed_transfers = (transfers[0], transfers[1])
        # B2 (speed program 2026-09-02): the closed-form exhaustion constraint
        # set, evaluated ABOVE the packet loop in the kernel's canonical order
        # (populated layers ascending, hexane before external water).  Nothing
        # between the transfers above and the kernel gate inside
        # ``legacy._packet_entries_after`` can refuse, so this raises exactly
        # the refusal that gate would raise, with the identical payload.  The
        # kernel gate is RETAINED, byte-frozen, as defense in depth (its module
        # is source-pinned by the RS-2 control authority and is not edited);
        # the payload agreement is held by test, not by assertion.
        #
        # F-B2-1 (2026-09-06): the kernel's water row is C8 Tier-1 branched -
        # a film-free layer whose model DECLARES a sorbed-water arm draws
        # retained water above the qualified floor, not the (zero) external
        # film - so the DECLARED ARM of each layer is handed to the closed
        # form through the kernel's own reader.  Without it the hoist refused
        # armed evaporating layers the kernel would have marched.
        if HOIST_ENABLED:
            exhaustion.raise_first_layer_phase_exhaustion(
                exhaustion.layer_phase_exhaustion_predicates(
                    populated_layers=exhaustion.populated_layer_ids(request.prior.tray),
                    inventories={item.layer: item for item in typed_inventories},
                    transfers={item.layer: item for item in typed_transfers},
                    macro_steps={layer: views[layer - 1] for layer in (1, 2)},
                    sorbed_arms={
                        layer: legacy._layer_sorbed_arm(request, layer) for layer in (1, 2)
                    },
                ),
                interval_end_time_s=request.end_time_s,
            )
        entries, packet_ledgers = _packet_entries_after_with_q_other_s(
            request,
            typed_transfers,
            typed_inventories,
            views,
            codec=codec,
            bundle=deposition_bundle,
        )
        if not all(ledger.passed for ledger in packet_ledgers):
            raise ThroughBedK2SingleWallStepError("decoded packet allocation ledger did not pass")
        callback = legacy._AcceptedLaw2PacketAdvance(
            entries=entries,
            residual_by_layer=tuple(solve.evaluation.scaled_residual_norm for solve in solves),
            maximum_scaled_residual=request.newton_tolerance,
        )
        host_request = tray_host.KCellTrayStepRequest(
            attempt_id=request.attempt_id,
            prior=request.prior.tray,
            end_time_s=request.end_time_s,
            flow_segments=request.flow_segments,
            external_cell_transfers=tuple(
                tray_host.CellExternalTransfer(
                    water_to_cell_kg=transfer.gas_boundary_water_to_cell_kg,
                    hexane_to_cell_kg=transfer.gas_boundary_hexane_to_cell_kg,
                    common_datum_energy_to_cell_j=math.fsum(
                        (
                            transfer.gas_and_wall_boundary_energy_to_cell_j,
                            (
                                0.0
                                if deposition_bundle is None
                                else deposition_bundle.ledger.layer_mechanical_work_j[layer - 1]
                            ),
                        )
                    ),
                    energy_datum_id=request.prior.tray.energy_datum_id,
                )
                for layer, transfer in enumerate(typed_transfers, start=1)
            ),
            relative_limit=request.relative_limit,
        )
        host_result = tray_host.evaluate_k_cell_tray_step(
            host_request,
            packet_callback=callback,
            k_planner=k_planner,
            auditor=codec,
        )
        if type(host_result) is tray_host.KCellTrayStepRejection:
            raise ThroughBedK2SingleWallStepError(
                f"K-cell host rejected step: {host_result.reason}"
            )
        if type(host_result) is not tray_host.ValidatedKCellTrayStep:
            raise ThroughBedK2SingleWallStepError("K-cell host returned a foreign transaction")
        thermal_whole_tray_energy = legacy._whole_tray_external_energy_ledger(
            request,
            inputs,
            views,
            typed_transfers,
            host_result,
            request_codec=codec,
        )
        whole_tray_energy: (
            WholeTrayExternalEnergyLedger
            | NativeWholeTrayExternalEnergyLedger
            | SingleWallQOtherSExternalEnergyLedger
        )
        if deposition_bundle is None:
            whole_tray_energy = thermal_whole_tray_energy
        else:
            whole_tray_energy = _q_other_s_whole_tray_external_energy_ledger(
                thermal_whole_tray_energy,
                host_result,
                deposition_bundle,
                relative_limit=request.relative_limit,
                # Q-F2a: needed only to reconstruct the native arm's booked
                # face defect energy; inert on the declared arm.
                macro_step_s=inputs[0].macro_step_s,
            )
        if not whole_tray_energy.passed:
            raise ThroughBedK2SingleWallStepError(
                "whole-tray external-only energy ledger did not pass"
            )
        candidate = AcceptedThroughBedK2SingleWallTrayState(
            contract=request.prior.contract,
            tray=host_result.candidate,
            integration_law_definition_digest=law_digest,
            layer_model_configuration_digest=model_digest,
            cell_field_states=tuple(solve.state for solve in solves),
            wall_temperature_k=physical_wall.wall_temperature_k,
        )
        coupled = _single_wall_coupled_ledger(
            request,
            host_result,
            solves,
            accepted,
            physical_wall,
            wall_ledger,
            typed_transfers,
            whole_tray_energy,
            tray_pressure,
            packet_ledgers,
        )
        if not coupled.passed:
            raise ThroughBedK2SingleWallStepError(
                "single-wall coupled packet/cell/wall ledger did not pass"
            )
        return ValidatedThroughBedK2SingleWallTrayStep(
            request=request,
            candidate=candidate,
            integration_law_definition_digest=law_digest,
            layer_model_configuration_digest=model_digest,
            request_definition_digest=request.definition_digest,
            host_transaction=host_result,
            layer_inventories=typed_inventories,
            cell_inputs=inputs,
            cell_fast_solves=solves,
            accepted_fast_steps=accepted,
            physical_wall=physical_wall,
            physical_wall_ledger=wall_ledger,
            layer_transfers=typed_transfers,
            inter_k_face=face,
            tray_pressure=tray_pressure,
            tray_printed_hull=tray_printed_hull,
            whole_tray_energy=whole_tray_energy,
            packet_exchange_ledgers=packet_ledgers,
            coupled_ledger=coupled,
            meal_mechanical_work_deposition_ledger=(
                None if deposition_bundle is None else deposition_bundle.ledger
            ),
            meal_mechanical_work_validation_context_digest=(
                None
                if meal_mechanical_work_validation_context is None
                else meal_mechanical_work_validation_context.definition_digest
            ),
            codec_registry_digest=codec.codec_registry_digest,
            auditor_identity_digest=codec.auditor_identity_digest,
            auditor_source_identity=codec.source_identity,
            _seal=_VALIDATED_SINGLE_WALL_SEAL,
        )
    except (
        ThroughBedK2SingleWallIntegrationError,
        kernel.ThroughBedK2KernelError,
        legacy.SP1K2Law2IntegrationError,
        ef.EngineeringFeasibilityError,
        high_loading.HighLoadingPacketAdapterError,
        component_datum.ComponentEnergyDatumError,
        meal_work.QSCMealMechanicalWorkDepositionError,
    ) as error:
        return ThroughBedK2SingleWallStepRejection(
            attempt_id=request.attempt_id,
            reason=str(error),
            rollback_state=request.prior,
            phase_exhaustion=(
                error.signal if isinstance(error, legacy.SP1K2Law2PhaseExhaustionError) else None
            ),
        )
    except Exception as error:
        return ThroughBedK2SingleWallStepRejection(
            attempt_id=request.attempt_id,
            reason=f"single-wall through-bed step failed: {type(error).__name__}: {error}",
            rollback_state=request.prior,
        )


def commit_through_bed_k2_single_wall_tray_step(
    current: AcceptedThroughBedK2SingleWallTrayState,
    transaction: ValidatedThroughBedK2SingleWallTrayStep,
    *,
    codec: high_loading.HighLoadingV1PacketAdapter,
    integration_law_authority: ThroughBedK2SingleWallIntegrationLawAuthority,
    expected_integration_law_definition_digest: str,
    expected_layer_model_configuration_digest: str,
    k_planner: tray_host.ManufacturedKRekeyPlanner | None = None,
    meal_mechanical_work_validation_context: (
        meal_work.QSCMealMechanicalWorkValidationContext | None
    ) = None,
) -> AcceptedThroughBedK2SingleWallTrayState:
    """Return the candidate only after complete deterministic replay and host commit."""

    if type(current) is not AcceptedThroughBedK2SingleWallTrayState:
        raise TypeError("single-wall commit requires an exact accepted state")
    if (
        type(transaction) is not ValidatedThroughBedK2SingleWallTrayStep
        or transaction._seal is not _VALIDATED_SINGLE_WALL_SEAL
    ):
        raise ThroughBedK2SingleWallIntegrationError(
            "single-wall commit requires a validated transaction"
        )
    ValidatedThroughBedK2SingleWallTrayStep.__post_init__(transaction)
    AcceptedThroughBedK2SingleWallTrayState.__post_init__(current)
    law_digest, model_digest = _validate_request_law_and_model_pins(
        transaction.request,
        integration_law_authority=integration_law_authority,
        expected_integration_law_definition_digest=(expected_integration_law_definition_digest),
        expected_layer_model_configuration_digest=expected_layer_model_configuration_digest,
    )
    if (
        current.integration_law_definition_digest != law_digest
        or transaction.integration_law_definition_digest != law_digest
        or current.layer_model_configuration_digest != model_digest
        or transaction.layer_model_configuration_digest != model_digest
    ):
        raise ThroughBedK2SingleWallIntegrationError(
            "single-wall current state or transaction violates caller-held law/model pins"
        )
    deposition_bundle = _validate_meal_mechanical_work_deposition(
        transaction.request,
        meal_mechanical_work_validation_context,
    )
    context_digest = (
        None
        if meal_mechanical_work_validation_context is None
        else meal_mechanical_work_validation_context.definition_digest
    )
    if (
        transaction.meal_mechanical_work_validation_context_digest != context_digest
        or transaction.meal_mechanical_work_deposition_ledger
        != (None if deposition_bundle is None else deposition_bundle.ledger)
    ):
        raise ThroughBedK2SingleWallIntegrationError(
            "single-wall commit meal-work bundle or caller-held pins differ"
        )
    if current is not transaction.request.prior:
        raise StaleThroughBedK2SingleWallStateError("single-wall transaction is stale or foreign")
    if current.state_digest != transaction.request.prior.state_digest:
        raise StaleThroughBedK2SingleWallStateError(
            "single-wall prior digest changed before commit"
        )
    if (
        codec.codec_registry_digest != transaction.codec_registry_digest
        or codec.auditor_identity_digest != transaction.auditor_identity_digest
        or codec.source_identity != transaction.auditor_source_identity
    ):
        raise ThroughBedK2SingleWallIntegrationError(
            "single-wall commit codec violates transaction pins"
        )
    if k_planner is None:
        k_planner = tray_host.NoKRekeyPlanner()
    _validate_codec_and_k_planner(transaction.request, codec, k_planner)
    reevaluated = evaluate_through_bed_k2_single_wall_tray_step(
        transaction.request,
        codec=codec,
        k_planner=k_planner,
        integration_law_authority=integration_law_authority,
        expected_integration_law_definition_digest=(expected_integration_law_definition_digest),
        expected_layer_model_configuration_digest=(expected_layer_model_configuration_digest),
        meal_mechanical_work_validation_context=(meal_mechanical_work_validation_context),
    )
    if (
        type(reevaluated) is not ValidatedThroughBedK2SingleWallTrayStep
        or reevaluated != transaction
    ):
        raise ThroughBedK2SingleWallIntegrationError(
            "single-wall transaction differs from complete commit-time reevaluation"
        )
    committed_tray = tray_host.commit_k_cell_tray_step(
        current.tray,
        transaction.host_transaction,
        auditor=codec,
    )
    if committed_tray is not transaction.candidate.tray:
        raise ThroughBedK2SingleWallIntegrationError("host and single-wall candidates diverged")
    if (
        not transaction.inter_k_face.passed
        or not transaction.tray_pressure.passed
        or not transaction.physical_wall_ledger.passed
        or not transaction.whole_tray_energy.passed
        or not all(ledger.passed for ledger in transaction.packet_exchange_ledgers)
        or not transaction.coupled_ledger.passed
    ):
        raise ThroughBedK2SingleWallIntegrationError(
            "single-wall acceptance ledgers no longer pass"
        )
    return transaction.candidate


__all__ = (
    "AcceptedThroughBedK2SingleWallTrayState",
    "AggregatePhysicalWallLedger",
    "AggregateThroughBedPressureLedger",
    "DECLARED_ASSUMPTIONS",
    "InterKGasFaceLedger",
    "LayerAcceptedTransfers",
    "LayerPacketExchangeLedger",
    "PRINTED_TRAY_DROP_HULL_PA",
    "PRINTED_TRAY_DROP_HULL_UNIT_OF_ACCOUNT_LAW_ID",
    "PacketDerivedLayerInventory",
    "ReferenceThroughBedK2Contract",
    "SINGLE_WALL_COUPLED_ENERGY_RESOLUTION_LAW_ID",
    "SINGLE_WALL_Q_OTHER_S_EXTERNAL_ENERGY_RESOLUTION_LAW_ID",
    "SINGLE_WALL_EXPECTED_INTEGRATION_LAW_DEFINITION_DIGEST",
    "SINGLE_WALL_INTEGRATION_DIGEST_DOMAIN",
    "SINGLE_WALL_INTEGRATION_LAW_AUTHORITY",
    "SINGLE_WALL_INTEGRATION_LAW_SCHEMA_ID",
    "SINGLE_WALL_INTEGRATION_LAW_SCHEMA_REVISION",
    "SingleWallCoupledTrayLedger",
    "SingleWallQOtherSExternalEnergyLedger",
    "StaleThroughBedK2SingleWallStateError",
    "ThroughBedGasBoundaryBinding",
    "ThroughBedGasBoundaryKind",
    "ThroughBedK2SingleWallConfigurationError",
    "ThroughBedK2SingleWallIntegrationError",
    "ThroughBedK2SingleWallIntegrationLawAuthority",
    "ThroughBedK2SingleWallStepError",
    "ThroughBedK2SingleWallStepRejection",
    "ThroughBedK2SingleWallStepRequest",
    "ThroughBedKLayerEngineeringModel",
    "ThroughBedTrayPrintedDropCorroboration",
    "ValidatedThroughBedK2SingleWallTrayStep",
    "WALL_ROW_RESOLUTION_BOUND_MULTIPLE",
    "WALL_ROW_RESOLUTION_BOUND_MULTIPLE_BRACKET",
    "WALL_ROW_RESOLUTION_BOUND_MULTIPLE_LAW_ID",
    "WholeTrayExternalEnergyLedger",
    "advance_aggregate_through_bed_wall_node",
    "canonical_single_wall_layer_model_digest",
    "classify_tray_pressure_drop",
    "commit_through_bed_k2_single_wall_tray_step",
    "evaluate_through_bed_k2_single_wall_tray_step",
    "inactive_upper_transfer_coefficients",
    "inactive_upper_wall_parameters",
    "migrate_equal_legacy_wall_state",
    "reference_through_bed_k2_contract",
    "require_corroborated_tray_drop",
    "validate_single_wall_integration_law_authority",
)
