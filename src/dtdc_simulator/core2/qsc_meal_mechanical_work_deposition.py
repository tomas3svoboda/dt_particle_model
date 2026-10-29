"""Source-pinned MECH-01C3a meal mechanical-work deposition contract.

This bounded seam turns an already validated MECH-01C2 nominal
``ALL_TO_MEAL`` allocation into per-packet ``BAL-ENERGY-S:Q_other_s`` inputs
for one physical tray.  It does not choose shaft power, repartition work, or
invent a source temperature.  Every row is bound to the accepted packet at
``t_n`` (payload identity, per-particle inventory, external weight, K layer,
and packet ID), and the complete local packet set must be covered exactly.

The bundle is engineering prequalification evidence.  It is not production
wiring, physical qualification, calibration, plant prediction, or an
industrial feasibility result.  The local tray integrations remain
responsible for applying each input inside their packet energy/phase solve and
for exposing it in their host and whole-system external-energy ledgers.
"""

from __future__ import annotations

import hashlib
import math
import struct
from dataclasses import dataclass, field
from fractions import Fraction
from typing import ClassVar, Iterable

from . import k_cell_tray_host as tray_host
from . import qsc_mechanical_runtime_authority as runtime
from . import qsc_shaft_work_energy_coupling as allocation
from .tray_particle_stateful import VerticalLayerId

QSC_MEAL_DEPOSITION_SCHEMA_ID = "dtdc-core2-qsc-meal-mechanical-work-deposition-v1"
QSC_MEAL_DEPOSITION_SCHEMA_REVISION = 1
QSC_MEAL_DEPOSITION_LAW_SCHEMA_ID = "dtdc-core2-qsc-meal-mechanical-work-deposition-law-v1"
QSC_MEAL_DEPOSITION_LAW_SCHEMA_REVISION = 1
QSC_MEAL_DEPOSITION_SLOT_DECLARATION = "BAL-ENERGY-S:Q_other_s"
QSC_MEAL_DEPOSITION_PACKET_BASIS_LAW_ID = "T_N_ACCEPTED_PACKET_IDENTITY_AND_EXTERNAL_WEIGHT_V1"
QSC_MEAL_DEPOSITION_COVERAGE_LAW_ID = "EXACT_LOCAL_PACKET_SET_ONCE_V1"
QSC_MEAL_DEPOSITION_CLOSURE_LAW_ID = (
    "EXACT_RATIONAL_PACKET_TO_LAYER_TO_TRAY_WITH_COUNTED_BINARY64_ULP_BOUND_V1"
)
QSC_MEAL_DEPOSITION_PER_PARTICLE_LAW_ID = (
    "DECLARED_PACKET_WORK_DIVIDED_BY_T_N_REPRESENTATIVE_PARTICLES_V1"
)
QSC_MEAL_DEPOSITION_EXPECTED_LAW_DEFINITION_DIGEST = (
    "sha256:2c6a9c619330a1ae9188bebc2392637de669a2e7a1e9ede5daf093e87f53b59c"
)


class QSCMealMechanicalWorkDepositionError(ValueError):
    """A deposition authority, bundle, pin, or accepted packet refused closed."""


def _require_identity(name: str, *values: str) -> None:
    if any(type(value) is not str or not value or value != value.strip() for value in values):
        raise QSCMealMechanicalWorkDepositionError(
            f"{name} must contain nonblank exact strings without surrounding whitespace"
        )


def _require_sha256(name: str, value: str) -> None:
    if (
        type(value) is not str
        or len(value) != 71
        or not value.startswith("sha256:")
        or any(character not in "0123456789abcdef" for character in value[7:])
    ):
        raise QSCMealMechanicalWorkDepositionError(
            f"{name} must be a canonical lowercase sha256 identity"
        )


def _require_float(
    name: str,
    *values: float,
    positive: bool = False,
    nonnegative: bool = False,
) -> None:
    if any(type(value) is not float or not math.isfinite(value) for value in values):
        raise QSCMealMechanicalWorkDepositionError(
            f"{name} must contain finite exact binary64 values"
        )
    if positive and any(value <= 0.0 for value in values):
        raise QSCMealMechanicalWorkDepositionError(f"{name} must be strictly positive")
    if nonnegative and any(value < 0.0 for value in values):
        raise QSCMealMechanicalWorkDepositionError(f"{name} must be nonnegative")


def _digest_part(value: str | int | float | bool) -> tuple[bytes, bytes]:
    if type(value) is str:
        return b"s", value.encode("utf-8")
    if type(value) is bool:
        return b"b", b"1" if value else b"0"
    if type(value) is int:
        return b"i", str(value).encode("ascii")
    if type(value) is float:
        return b"f", struct.pack(">d", value)
    raise TypeError(f"unsupported meal-deposition digest part: {type(value).__name__}")


def _digest(domain: str, parts: Iterable[str | int | float | bool]) -> str:
    digest = hashlib.sha256()
    digest.update(b"DTDC-QSC-MEAL-MECHANICAL-WORK-DEPOSITION-V1\0")
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


def _fraction_parts(value: Fraction) -> tuple[int, int]:
    return value.numerator, value.denominator


def _canonical_fraction(numerator: int, denominator: int, *, name: str) -> Fraction:
    if type(numerator) is not int or type(denominator) is not int or denominator <= 0:
        raise QSCMealMechanicalWorkDepositionError(
            f"{name} must be a canonical exact integer ratio with positive denominator"
        )
    value = Fraction(numerator, denominator)
    if (value.numerator, value.denominator) != (numerator, denominator):
        raise QSCMealMechanicalWorkDepositionError(
            f"{name} exact ratio must be in canonical lowest terms"
        )
    return value


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMealMechanicalWorkDepositionLawAuthority:
    """Fixed interpretation of a nominal C2 meal-work allocation."""

    schema_id: str
    schema_revision: int
    source_runtime_definition_digest: str
    source_allocator_law_definition_digest: str
    source_allocation_schema_id: str
    source_allocation_schema_revision: int
    boundary_partition: runtime.ShaftWorkBoundaryPartition
    within_meal_distribution: runtime.WithinMealDistribution
    meal_slot_declaration: str
    packet_basis_law_id: str
    coverage_law_id: str
    closure_law_id: str
    per_particle_law_id: str
    source_temperature_defined: bool

    declared_engineering_law: ClassVar[bool] = True
    production_wired: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    industry_feasible_operating_envelope: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if (
            self.schema_id != QSC_MEAL_DEPOSITION_LAW_SCHEMA_ID
            or type(self.schema_revision) is not int
            or self.schema_revision != QSC_MEAL_DEPOSITION_LAW_SCHEMA_REVISION
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition law has a foreign schema identity or revision"
            )
        _require_sha256("source runtime definition", self.source_runtime_definition_digest)
        _require_sha256(
            "source allocator-law definition", self.source_allocator_law_definition_digest
        )
        if (
            self.source_runtime_definition_digest
            != runtime.QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST
            or self.source_allocator_law_definition_digest
            != allocation.QSC_SHAFT_WORK_EXPECTED_ALLOCATOR_LAW_DEFINITION_DIGEST
            or self.source_allocation_schema_id != allocation.QSC_SHAFT_WORK_SCHEMA_ID
            or type(self.source_allocation_schema_revision) is not int
            or self.source_allocation_schema_revision != allocation.QSC_SHAFT_WORK_SCHEMA_REVISION
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition law binds a foreign runtime, allocator, or allocation schema"
            )
        expected = (
            (
                self.boundary_partition,
                runtime.ShaftWorkBoundaryPartition.ALL_TO_MEAL,
            ),
            (
                self.within_meal_distribution,
                runtime.WithinMealDistribution.WET_HOLDUP_PROPORTIONAL_BY_TRAY,
            ),
            (self.meal_slot_declaration, QSC_MEAL_DEPOSITION_SLOT_DECLARATION),
            (self.packet_basis_law_id, QSC_MEAL_DEPOSITION_PACKET_BASIS_LAW_ID),
            (self.coverage_law_id, QSC_MEAL_DEPOSITION_COVERAGE_LAW_ID),
            (self.closure_law_id, QSC_MEAL_DEPOSITION_CLOSURE_LAW_ID),
            (self.per_particle_law_id, QSC_MEAL_DEPOSITION_PER_PARTICLE_LAW_ID),
        )
        if any(observed != ruled for observed, ruled in expected):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition law differs from the nominal Q_other_s ruling"
            )
        if self.source_temperature_defined is not False:
            raise QSCMealMechanicalWorkDepositionError(
                "mechanical work must not acquire a fictitious source temperature"
            )

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-meal-mechanical-work-deposition-law-authority",
            (
                self.schema_id,
                self.schema_revision,
                self.source_runtime_definition_digest,
                self.source_allocator_law_definition_digest,
                self.source_allocation_schema_id,
                self.source_allocation_schema_revision,
                self.boundary_partition.value,
                self.within_meal_distribution.value,
                self.meal_slot_declaration,
                self.packet_basis_law_id,
                self.coverage_law_id,
                self.closure_law_id,
                self.per_particle_law_id,
                self.source_temperature_defined,
            ),
        )


QSC_MEAL_MECHANICAL_WORK_DEPOSITION_LAW_AUTHORITY = QSCMealMechanicalWorkDepositionLawAuthority(
    schema_id=QSC_MEAL_DEPOSITION_LAW_SCHEMA_ID,
    schema_revision=QSC_MEAL_DEPOSITION_LAW_SCHEMA_REVISION,
    source_runtime_definition_digest=runtime.QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST,
    source_allocator_law_definition_digest=(
        allocation.QSC_SHAFT_WORK_EXPECTED_ALLOCATOR_LAW_DEFINITION_DIGEST
    ),
    source_allocation_schema_id=allocation.QSC_SHAFT_WORK_SCHEMA_ID,
    source_allocation_schema_revision=allocation.QSC_SHAFT_WORK_SCHEMA_REVISION,
    boundary_partition=runtime.ShaftWorkBoundaryPartition.ALL_TO_MEAL,
    within_meal_distribution=(runtime.WithinMealDistribution.WET_HOLDUP_PROPORTIONAL_BY_TRAY),
    meal_slot_declaration=QSC_MEAL_DEPOSITION_SLOT_DECLARATION,
    packet_basis_law_id=QSC_MEAL_DEPOSITION_PACKET_BASIS_LAW_ID,
    coverage_law_id=QSC_MEAL_DEPOSITION_COVERAGE_LAW_ID,
    closure_law_id=QSC_MEAL_DEPOSITION_CLOSURE_LAW_ID,
    per_particle_law_id=QSC_MEAL_DEPOSITION_PER_PARTICLE_LAW_ID,
    source_temperature_defined=False,
)


def validate_qsc_meal_mechanical_work_deposition_law_authority(
    authority: QSCMealMechanicalWorkDepositionLawAuthority,
) -> str:
    if type(authority) is not QSCMealMechanicalWorkDepositionLawAuthority:
        raise TypeError("meal-deposition law must be an exact typed authority")
    QSCMealMechanicalWorkDepositionLawAuthority.__post_init__(authority)
    digest = authority.definition_digest
    if digest != QSC_MEAL_DEPOSITION_EXPECTED_LAW_DEFINITION_DIGEST:
        raise QSCMealMechanicalWorkDepositionError(
            "meal-deposition law differs from its fixed expected definition digest"
        )
    live = (
        runtime.QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST,
        allocation.QSC_SHAFT_WORK_EXPECTED_ALLOCATOR_LAW_DEFINITION_DIGEST,
        allocation.QSC_SHAFT_WORK_SCHEMA_ID,
        allocation.QSC_SHAFT_WORK_SCHEMA_REVISION,
        QSC_MEAL_DEPOSITION_SLOT_DECLARATION,
        QSC_MEAL_DEPOSITION_PACKET_BASIS_LAW_ID,
        QSC_MEAL_DEPOSITION_COVERAGE_LAW_ID,
        QSC_MEAL_DEPOSITION_CLOSURE_LAW_ID,
        QSC_MEAL_DEPOSITION_PER_PARTICLE_LAW_ID,
    )
    pinned = (
        authority.source_runtime_definition_digest,
        authority.source_allocator_law_definition_digest,
        authority.source_allocation_schema_id,
        authority.source_allocation_schema_revision,
        authority.meal_slot_declaration,
        authority.packet_basis_law_id,
        authority.coverage_law_id,
        authority.closure_law_id,
        authority.per_particle_law_id,
    )
    if live != pinned:
        raise QSCMealMechanicalWorkDepositionError(
            "operative meal-deposition or source identity was substituted"
        )
    return digest


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMealMechanicalWorkValidationContext:
    """Caller-held authorities and exact pins required on evaluate and commit."""

    runtime_authority: runtime.QSCMechanicalRuntimeAuthority
    allocator_law_authority: allocation.QSCShaftWorkAllocatorLawAuthority
    deposition_law_authority: QSCMealMechanicalWorkDepositionLawAuthority
    expected_runtime_definition_digest: str
    expected_allocator_law_definition_digest: str
    expected_deposition_law_definition_digest: str
    expected_allocation_request_definition_digest: str
    expected_allocation_definition_digest: str
    #: Owner-approved B9 (2026-08-31): a caller that ran the FULL runtime
    #: validation this interval threads its digest here; every context
    #: re-validation then pins the digest (value-keyed re-hash catches any
    #: mutation) instead of repeating the structural re-walk.  ``None``
    #: keeps the sealed full behavior.  Excluded from definition_digest.
    preverified_runtime_definition_digest: str | None = None

    physically_qualifying: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False

    def __post_init__(self) -> None:
        for name in (
            "expected_runtime_definition_digest",
            "expected_allocator_law_definition_digest",
            "expected_deposition_law_definition_digest",
            "expected_allocation_request_definition_digest",
            "expected_allocation_definition_digest",
        ):
            _require_sha256(name, getattr(self, name))
        if self.preverified_runtime_definition_digest is not None:
            _require_sha256(
                "preverified_runtime_definition_digest",
                self.preverified_runtime_definition_digest,
            )
        observed_runtime = runtime.validate_qsc_mechanical_runtime_authority(
            self.runtime_authority,
            preverified_definition_digest=self.preverified_runtime_definition_digest,
        )
        observed_allocator = allocation.validate_qsc_shaft_work_allocator_law_authority(
            self.allocator_law_authority
        )
        observed_deposition = validate_qsc_meal_mechanical_work_deposition_law_authority(
            self.deposition_law_authority
        )
        if (
            observed_runtime != self.expected_runtime_definition_digest
            or observed_allocator != self.expected_allocator_law_definition_digest
            or observed_deposition != self.expected_deposition_law_definition_digest
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "caller-held runtime, allocator, or deposition-law pin differs"
            )

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-meal-mechanical-work-validation-context",
            (
                self.expected_runtime_definition_digest,
                self.expected_allocator_law_definition_digest,
                self.expected_deposition_law_definition_digest,
                self.expected_allocation_request_definition_digest,
                self.expected_allocation_definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMealPacketMechanicalWorkDeposition:
    physical_tray_id: str
    vertical_layer_id: VerticalLayerId
    packet_id: str
    t_n_payload_identity_digest: str
    t_n_packet_definition_digest: str
    t_n_representative_particles: float
    declared_mechanical_work_j: float
    mechanical_work_per_particle_j: float
    closing_packet_remainder: bool

    source_temperature_defined: ClassVar[bool] = False
    is_heat_source: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("meal-deposition packet identity", self.physical_tray_id, self.packet_id)
        if self.physical_tray_id not in runtime.QSC_MECH01_SOLID_TRAY_ORDER:
            raise QSCMealMechanicalWorkDepositionError("meal-deposition row names a foreign tray")
        if type(self.vertical_layer_id) is not VerticalLayerId:
            raise TypeError("meal-deposition layer must be an exact VerticalLayerId")
        if self.vertical_layer_id.value not in (1, 2):
            raise QSCMealMechanicalWorkDepositionError("meal-deposition layer must be K1 or K2")
        _require_sha256("t_n payload identity", self.t_n_payload_identity_digest)
        _require_sha256("t_n packet definition", self.t_n_packet_definition_digest)
        _require_float(
            "meal-deposition packet work and basis",
            self.t_n_representative_particles,
            self.declared_mechanical_work_j,
            self.mechanical_work_per_particle_j,
            positive=True,
        )
        if self.declared_mechanical_work_j / self.t_n_representative_particles != (
            self.mechanical_work_per_particle_j
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "per-particle mechanical work differs from the ruled t_n weight division"
            )
        if type(self.closing_packet_remainder) is not bool:
            raise TypeError("meal-deposition closing marker must be an exact bool")

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-meal-packet-mechanical-work-deposition",
            (
                self.physical_tray_id,
                self.vertical_layer_id.value,
                self.packet_id,
                self.t_n_payload_identity_digest,
                self.t_n_packet_definition_digest,
                self.t_n_representative_particles,
                self.declared_mechanical_work_j,
                self.mechanical_work_per_particle_j,
                self.closing_packet_remainder,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMealMechanicalWorkDepositionLedger:
    physical_tray_id: str
    energy_datum_id: str
    meal_slot_declaration: str
    packet_basis_law_id: str
    coverage_law_id: str
    closure_law_id: str
    covered_packet_ids: tuple[str, ...]
    layer_mechanical_work_j: tuple[float, float]
    tray_declared_mechanical_work_j: float
    emitted_packet_mechanical_work_j: float
    binary64_closure_residual_j: float
    exact_layer_1_numerator: int
    exact_layer_1_denominator: int
    exact_layer_2_numerator: int
    exact_layer_2_denominator: int
    exact_tray_numerator: int
    exact_tray_denominator: int
    exact_emitted_numerator: int
    exact_emitted_denominator: int
    exact_residual_numerator: int
    exact_residual_denominator: int
    closure_ulp_limit: int
    binary64_roundoff_bound_j: float
    authoritative_external_energy_input: bool

    exact_packet_layer_tray_closure: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity(
            "meal-deposition ledger identity",
            self.physical_tray_id,
            self.energy_datum_id,
            self.meal_slot_declaration,
            self.packet_basis_law_id,
            self.coverage_law_id,
            self.closure_law_id,
        )
        if (
            self.meal_slot_declaration != QSC_MEAL_DEPOSITION_SLOT_DECLARATION
            or self.packet_basis_law_id != QSC_MEAL_DEPOSITION_PACKET_BASIS_LAW_ID
            or self.coverage_law_id != QSC_MEAL_DEPOSITION_COVERAGE_LAW_ID
            or self.closure_law_id != QSC_MEAL_DEPOSITION_CLOSURE_LAW_ID
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition ledger has a foreign slot or arithmetic law"
            )
        if (
            type(self.covered_packet_ids) is not tuple
            or not self.covered_packet_ids
            or any(type(item) is not str or not item for item in self.covered_packet_ids)
            or len(set(self.covered_packet_ids)) != len(self.covered_packet_ids)
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition ledger requires one unique nonempty packet tuple"
            )
        if (
            type(self.layer_mechanical_work_j) is not tuple
            or len(self.layer_mechanical_work_j) != 2
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition ledger requires exactly two K-layer totals"
            )
        _require_float(
            "meal-deposition ledger work",
            *self.layer_mechanical_work_j,
            self.tray_declared_mechanical_work_j,
            self.emitted_packet_mechanical_work_j,
            positive=True,
        )
        _require_float("meal-deposition binary64 closure", self.binary64_closure_residual_j)
        _require_float(
            "meal-deposition binary64 bound", self.binary64_roundoff_bound_j, positive=True
        )
        layer_1 = _canonical_fraction(
            self.exact_layer_1_numerator,
            self.exact_layer_1_denominator,
            name="exact K1 work",
        )
        layer_2 = _canonical_fraction(
            self.exact_layer_2_numerator,
            self.exact_layer_2_denominator,
            name="exact K2 work",
        )
        tray = _canonical_fraction(
            self.exact_tray_numerator,
            self.exact_tray_denominator,
            name="exact tray work",
        )
        emitted = _canonical_fraction(
            self.exact_emitted_numerator,
            self.exact_emitted_denominator,
            name="exact emitted work",
        )
        residual = _canonical_fraction(
            self.exact_residual_numerator,
            self.exact_residual_denominator,
            name="exact closure residual",
        )
        if (
            float(layer_1) != self.layer_mechanical_work_j[0]
            or float(layer_2) != self.layer_mechanical_work_j[1]
            or tray != Fraction.from_float(self.tray_declared_mechanical_work_j)
            or float(emitted) != self.emitted_packet_mechanical_work_j
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition exact ratios differ from their binary64 fields"
            )
        if layer_1 + layer_2 != emitted or residual != emitted - tray:
            raise QSCMealMechanicalWorkDepositionError(
                "packet-to-layer-to-tray exact work closure does not reconstruct"
            )
        if self.binary64_closure_residual_j != (
            self.emitted_packet_mechanical_work_j - self.tray_declared_mechanical_work_j
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "binary64 meal-deposition closure residual does not reconstruct"
            )
        expected_ulp_limit = 2 * len(self.covered_packet_ids) + 2
        if type(self.closure_ulp_limit) is not int or self.closure_ulp_limit != (
            expected_ulp_limit
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition ULP limit differs from the typed counted-sum law"
            )
        expected_bound = self.closure_ulp_limit * math.ulp(self.tray_declared_mechanical_work_j)
        if self.binary64_roundoff_bound_j != expected_bound:
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition roundoff bound differs from the typed ULP law"
            )
        if abs(residual) > Fraction.from_float(expected_bound):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition exact closure exceeds its disclosed ULP bound"
            )
        if self.authoritative_external_energy_input is not True:
            raise QSCMealMechanicalWorkDepositionError(
                "Q_other_s must remain an authoritative external common-energy input"
            )

    @property
    def passed(self) -> bool:
        try:
            self.__post_init__()
        except TypeError, QSCMealMechanicalWorkDepositionError:
            return False
        return True

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-meal-mechanical-work-deposition-ledger",
            (
                self.physical_tray_id,
                self.energy_datum_id,
                self.meal_slot_declaration,
                self.packet_basis_law_id,
                self.coverage_law_id,
                self.closure_law_id,
                *self.covered_packet_ids,
                *self.layer_mechanical_work_j,
                self.tray_declared_mechanical_work_j,
                self.emitted_packet_mechanical_work_j,
                self.binary64_closure_residual_j,
                self.exact_layer_1_numerator,
                self.exact_layer_1_denominator,
                self.exact_layer_2_numerator,
                self.exact_layer_2_denominator,
                self.exact_tray_numerator,
                self.exact_tray_denominator,
                self.exact_emitted_numerator,
                self.exact_emitted_denominator,
                self.exact_residual_numerator,
                self.exact_residual_denominator,
                self.closure_ulp_limit,
                self.binary64_roundoff_bound_j,
                self.authoritative_external_energy_input,
            ),
        )


_QSC_MEAL_DEPOSITION_BUNDLE_SEAL = object()


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMealMechanicalWorkDepositionBundle:
    schema_id: str
    schema_revision: int
    deposition_law_definition_digest: str
    validation_context_definition_digest: str
    runtime_definition_digest: str
    allocator_law_definition_digest: str
    allocation_request_definition_digest: str
    allocation_definition_digest: str
    accepted_step_identity_digest: str
    prior_tray_state_digest: str
    physical_tray_id: str
    t_n_s: float
    end_time_s: float
    energy_datum_id: str
    allocation_request: allocation.QSCMechanicalAllocationRequest
    mechanical_work_allocation: allocation.QSCMechanicalWorkAllocation
    packet_depositions: tuple[QSCMealPacketMechanicalWorkDeposition, ...]
    ledger: QSCMealMechanicalWorkDepositionLedger
    _seal: object = field(repr=False, compare=False)

    bounded_engineering_deposition: ClassVar[bool] = True
    source_temperature_defined: ClassVar[bool] = False
    production_wired: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    f2: ClassVar[bool] = False
    f3: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    industry_feasible_operating_envelope: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self._seal is not _QSC_MEAL_DEPOSITION_BUNDLE_SEAL:
            raise TypeError("meal-deposition bundles are issued only by the validated builder")
        if (
            self.schema_id != QSC_MEAL_DEPOSITION_SCHEMA_ID
            or type(self.schema_revision) is not int
            or self.schema_revision != QSC_MEAL_DEPOSITION_SCHEMA_REVISION
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition bundle has a foreign schema identity or revision"
            )
        for name in (
            "deposition_law_definition_digest",
            "validation_context_definition_digest",
            "runtime_definition_digest",
            "allocator_law_definition_digest",
            "allocation_request_definition_digest",
            "allocation_definition_digest",
            "accepted_step_identity_digest",
            "prior_tray_state_digest",
        ):
            _require_sha256(name, getattr(self, name))
        _require_identity("meal-deposition tray/datum", self.physical_tray_id, self.energy_datum_id)
        _require_float("meal-deposition endpoint", self.t_n_s, self.end_time_s, nonnegative=True)
        if self.end_time_s <= self.t_n_s:
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition endpoint must advance accepted time"
            )
        if type(self.allocation_request) is not allocation.QSCMechanicalAllocationRequest:
            raise TypeError("meal-deposition bundle contains a foreign C2 request")
        if type(self.mechanical_work_allocation) is not allocation.QSCMechanicalWorkAllocation:
            raise TypeError("meal-deposition bundle contains a foreign C2 allocation")
        allocation.QSCMechanicalAllocationRequest.__post_init__(self.allocation_request)
        allocation.QSCMechanicalWorkAllocation.__post_init__(self.mechanical_work_allocation)
        if (
            self.allocation_request.definition_digest != self.allocation_request_definition_digest
            or self.mechanical_work_allocation.definition_digest
            != self.allocation_definition_digest
            or self.accepted_step_identity_digest
            != self.allocation_request.accepted_step_identity_digest
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition bundle differs from its embedded C2 identities"
            )
        if type(self.packet_depositions) is not tuple or not self.packet_depositions:
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition bundle requires a nonempty exact packet tuple"
            )
        for row in self.packet_depositions:
            if type(row) is not QSCMealPacketMechanicalWorkDeposition:
                raise TypeError("meal-deposition bundle contains a foreign packet row")
            QSCMealPacketMechanicalWorkDeposition.__post_init__(row)
        if sum(row.closing_packet_remainder for row in self.packet_depositions) != 1 or not (
            self.packet_depositions[-1].closing_packet_remainder
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "only the final canonical tray packet may close the C2 remainder"
            )
        if type(self.ledger) is not QSCMealMechanicalWorkDepositionLedger:
            raise TypeError("meal-deposition bundle contains a foreign ledger")
        QSCMealMechanicalWorkDepositionLedger.__post_init__(self.ledger)
        if (
            self.ledger.physical_tray_id != self.physical_tray_id
            or self.ledger.energy_datum_id != self.energy_datum_id
            or self.ledger.covered_packet_ids
            != tuple(row.packet_id for row in self.packet_depositions)
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition bundle and local closure ledger disagree"
            )

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-meal-mechanical-work-deposition-bundle",
            (
                self.schema_id,
                self.schema_revision,
                self.deposition_law_definition_digest,
                self.validation_context_definition_digest,
                self.runtime_definition_digest,
                self.allocator_law_definition_digest,
                self.allocation_request_definition_digest,
                self.allocation_definition_digest,
                self.accepted_step_identity_digest,
                self.prior_tray_state_digest,
                self.physical_tray_id,
                self.t_n_s,
                self.end_time_s,
                self.energy_datum_id,
                *(row.definition_digest for row in self.packet_depositions),
                self.ledger.definition_digest,
            ),
        )

    def packet_row(self, packet_id: str) -> QSCMealPacketMechanicalWorkDeposition:
        rows = tuple(row for row in self.packet_depositions if row.packet_id == packet_id)
        if len(rows) != 1:
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition packet lookup did not resolve exactly once"
            )
        return rows[0]


def _validate_context_and_c2(
    *,
    context: QSCMealMechanicalWorkValidationContext,
    request: allocation.QSCMechanicalAllocationRequest,
    work_allocation: allocation.QSCMechanicalWorkAllocation,
) -> tuple[str, str, str]:
    if type(context) is not QSCMealMechanicalWorkValidationContext:
        raise TypeError("meal-deposition validation requires an exact caller-held context")
    QSCMealMechanicalWorkValidationContext.__post_init__(context)
    if request.definition_digest != context.expected_allocation_request_definition_digest:
        raise QSCMealMechanicalWorkDepositionError(
            "C2 allocation request differs from the caller-held request pin"
        )
    allocation_digest = allocation.validate_qsc_shaft_work_allocation(
        work_allocation,
        request=request,
        runtime_authority=context.runtime_authority,
        allocator_law_authority=context.allocator_law_authority,
        expected_request_definition_digest=(context.expected_allocation_request_definition_digest),
        expected_allocator_law_definition_digest=(context.expected_allocator_law_definition_digest),
        preverified_runtime_definition_digest=(context.preverified_runtime_definition_digest),
    )
    if allocation_digest != context.expected_allocation_definition_digest:
        raise QSCMealMechanicalWorkDepositionError(
            "C2 allocation differs from the caller-held allocation pin"
        )
    if (
        request.boundary_partition is not runtime.ShaftWorkBoundaryPartition.ALL_TO_MEAL
        or request.tray_distribution
        is not runtime.WithinMealDistribution.WET_HOLDUP_PROPORTIONAL_BY_TRAY
        or work_allocation.ledger.boundary_partition
        is not runtime.ShaftWorkBoundaryPartition.ALL_TO_MEAL
    ):
        raise QSCMealMechanicalWorkDepositionError(
            "C3a admits only the ruled nominal ALL_TO_MEAL wet-holdup allocation"
        )
    return (
        context.expected_runtime_definition_digest,
        context.expected_allocator_law_definition_digest,
        allocation_digest,
    )


def _build_bundle(
    *,
    accepted_tray: tray_host.AcceptedKCellTrayState,
    end_time_s: float,
    request: allocation.QSCMechanicalAllocationRequest,
    work_allocation: allocation.QSCMechanicalWorkAllocation,
    context: QSCMealMechanicalWorkValidationContext,
) -> QSCMealMechanicalWorkDepositionBundle:
    if type(accepted_tray) is not tray_host.AcceptedKCellTrayState:
        raise TypeError("meal-deposition requires an exact accepted K-cell tray")
    tray_host.AcceptedKCellTrayState.__post_init__(accepted_tray)
    _require_float("meal-deposition end time", end_time_s, nonnegative=True)
    if end_time_s <= accepted_tray.time_s:
        raise QSCMealMechanicalWorkDepositionError(
            "meal-deposition endpoint must advance accepted tray time"
        )
    runtime_digest, allocator_digest, allocation_digest = _validate_context_and_c2(
        context=context,
        request=request,
        work_allocation=work_allocation,
    )
    dt = end_time_s - accepted_tray.time_s
    if request.dt_s != dt:
        raise QSCMealMechanicalWorkDepositionError(
            "C2 allocation duration differs from the exact local tray macro-step"
        )
    if request.energy_datum_id != accepted_tray.energy_datum_id:
        raise QSCMealMechanicalWorkDepositionError(
            "C2 allocation and accepted tray use different common-energy datums"
        )
    tray_id = accepted_tray.physical_tray_id
    host_packets = tuple(
        sorted(
            accepted_tray.packets,
            key=lambda packet: (
                packet.owner_key.vertical_layer_id.value,
                packet.packet_id,
            ),
        )
    )
    request_packets = tuple(
        packet for packet in request.packets if packet.physical_tray_id == tray_id
    )
    allocated_rows = tuple(
        row for row in work_allocation.meal_packet_inputs if row.physical_tray_id == tray_id
    )
    if (
        not host_packets
        or len(host_packets) != len(request_packets)
        or len(host_packets) != len(allocated_rows)
    ):
        raise QSCMealMechanicalWorkDepositionError(
            "C2 allocation does not cover the complete accepted local packet set exactly"
        )
    packet_depositions: list[QSCMealPacketMechanicalWorkDeposition] = []
    for host_packet, request_packet, allocated in zip(
        host_packets, request_packets, allocated_rows, strict=True
    ):
        expected_key = (
            tray_id,
            host_packet.owner_key.vertical_layer_id,
            host_packet.packet_id,
        )
        if (
            request_packet.physical_tray_id,
            request_packet.vertical_layer_id,
            request_packet.packet_id,
        ) != expected_key or (
            allocated.physical_tray_id,
            allocated.vertical_layer_id,
            allocated.packet_id,
        ) != expected_key:
            raise QSCMealMechanicalWorkDepositionError(
                "C2 allocation packet identity/order differs from the accepted t_n packet set"
            )
        if (
            request_packet.particle_state_identity_digest != host_packet.payload_identity
            or request_packet.inventory != host_packet.entry.inventory
            or request_packet.weight != host_packet.entry.weight
            or allocated.packet_definition_digest != request_packet.definition_digest
            or allocated.accepted_weighted_wet_mass_kg
            != request_packet.accepted_weighted_wet_mass_kg
        ):
            raise QSCMealMechanicalWorkDepositionError(
                "C2 packet payload, inventory, weight, or definition differs from t_n"
            )
        per_particle = (
            allocated.declared_mechanical_work_j / host_packet.entry.weight.representative_particles
        )
        _require_float("per-particle Q_other_s input", per_particle, positive=True)
        packet_depositions.append(
            QSCMealPacketMechanicalWorkDeposition(
                physical_tray_id=tray_id,
                vertical_layer_id=host_packet.owner_key.vertical_layer_id,
                packet_id=host_packet.packet_id,
                t_n_payload_identity_digest=host_packet.payload_identity,
                t_n_packet_definition_digest=request_packet.definition_digest,
                t_n_representative_particles=(host_packet.entry.weight.representative_particles),
                declared_mechanical_work_j=allocated.declared_mechanical_work_j,
                mechanical_work_per_particle_j=per_particle,
                closing_packet_remainder=allocated.closing_packet_remainder,
            )
        )
    rows = tuple(packet_depositions)
    tray_row = next(
        row for row in work_allocation.tray_allocations if row.physical_tray_id == tray_id
    )
    exact_layers = tuple(
        sum(
            (
                Fraction.from_float(row.declared_mechanical_work_j)
                for row in rows
                if row.vertical_layer_id.value == layer
            ),
            Fraction(),
        )
        for layer in (1, 2)
    )
    if any(value <= 0 for value in exact_layers):
        raise QSCMealMechanicalWorkDepositionError(
            "nominal C3a requires positive mechanical work in both occupied K layers"
        )
    layer_values = tuple(float(value) for value in exact_layers)
    exact_emitted = sum(
        (Fraction.from_float(row.declared_mechanical_work_j) for row in rows),
        Fraction(),
    )
    emitted = float(exact_emitted)
    exact_tray = Fraction.from_float(tray_row.declared_mechanical_work_j)
    exact_residual = exact_emitted - exact_tray
    closure_ulp_limit = 2 * len(rows) + 2
    roundoff_bound = closure_ulp_limit * math.ulp(tray_row.declared_mechanical_work_j)
    layer_1_parts = _fraction_parts(exact_layers[0])
    layer_2_parts = _fraction_parts(exact_layers[1])
    tray_parts = _fraction_parts(exact_tray)
    emitted_parts = _fraction_parts(exact_emitted)
    residual_parts = _fraction_parts(exact_residual)
    ledger = QSCMealMechanicalWorkDepositionLedger(
        physical_tray_id=tray_id,
        energy_datum_id=accepted_tray.energy_datum_id,
        meal_slot_declaration=QSC_MEAL_DEPOSITION_SLOT_DECLARATION,
        packet_basis_law_id=QSC_MEAL_DEPOSITION_PACKET_BASIS_LAW_ID,
        coverage_law_id=QSC_MEAL_DEPOSITION_COVERAGE_LAW_ID,
        closure_law_id=QSC_MEAL_DEPOSITION_CLOSURE_LAW_ID,
        covered_packet_ids=tuple(row.packet_id for row in rows),
        layer_mechanical_work_j=(layer_values[0], layer_values[1]),
        tray_declared_mechanical_work_j=tray_row.declared_mechanical_work_j,
        emitted_packet_mechanical_work_j=emitted,
        binary64_closure_residual_j=(emitted - tray_row.declared_mechanical_work_j),
        exact_layer_1_numerator=layer_1_parts[0],
        exact_layer_1_denominator=layer_1_parts[1],
        exact_layer_2_numerator=layer_2_parts[0],
        exact_layer_2_denominator=layer_2_parts[1],
        exact_tray_numerator=tray_parts[0],
        exact_tray_denominator=tray_parts[1],
        exact_emitted_numerator=emitted_parts[0],
        exact_emitted_denominator=emitted_parts[1],
        exact_residual_numerator=residual_parts[0],
        exact_residual_denominator=residual_parts[1],
        closure_ulp_limit=closure_ulp_limit,
        binary64_roundoff_bound_j=roundoff_bound,
        authoritative_external_energy_input=True,
    )
    # Retain the exact C2 result as an independent check; the local ledger
    # reports rounded layer sums, while this check operates on emitted rows.
    if abs(exact_residual) > Fraction.from_float(roundoff_bound):
        raise QSCMealMechanicalWorkDepositionError(
            "local C2 tray allocation exceeds the counted exact-rational ULP closure"
        )
    law_digest = validate_qsc_meal_mechanical_work_deposition_law_authority(
        context.deposition_law_authority
    )
    return QSCMealMechanicalWorkDepositionBundle(
        schema_id=QSC_MEAL_DEPOSITION_SCHEMA_ID,
        schema_revision=QSC_MEAL_DEPOSITION_SCHEMA_REVISION,
        deposition_law_definition_digest=law_digest,
        validation_context_definition_digest=context.definition_digest,
        runtime_definition_digest=runtime_digest,
        allocator_law_definition_digest=allocator_digest,
        allocation_request_definition_digest=request.definition_digest,
        allocation_definition_digest=allocation_digest,
        accepted_step_identity_digest=request.accepted_step_identity_digest,
        prior_tray_state_digest=accepted_tray.state_digest,
        physical_tray_id=tray_id,
        t_n_s=accepted_tray.time_s,
        end_time_s=end_time_s,
        energy_datum_id=accepted_tray.energy_datum_id,
        allocation_request=request,
        mechanical_work_allocation=work_allocation,
        packet_depositions=rows,
        ledger=ledger,
        _seal=_QSC_MEAL_DEPOSITION_BUNDLE_SEAL,
    )


def build_qsc_meal_mechanical_work_deposition_bundle(
    *,
    accepted_tray: tray_host.AcceptedKCellTrayState,
    end_time_s: float,
    allocation_request: allocation.QSCMechanicalAllocationRequest,
    mechanical_work_allocation: allocation.QSCMechanicalWorkAllocation,
    validation_context: QSCMealMechanicalWorkValidationContext,
) -> QSCMealMechanicalWorkDepositionBundle:
    """Build one sealed local Q_other_s bundle from exact C2 and t_n authorities."""

    return _build_bundle(
        accepted_tray=accepted_tray,
        end_time_s=end_time_s,
        request=allocation_request,
        work_allocation=mechanical_work_allocation,
        context=validation_context,
    )


def validate_qsc_meal_mechanical_work_deposition_bundle(
    bundle: QSCMealMechanicalWorkDepositionBundle,
    *,
    accepted_tray: tray_host.AcceptedKCellTrayState,
    end_time_s: float,
    validation_context: QSCMealMechanicalWorkValidationContext,
    preverified_definition_digest: str | None = None,
) -> str:
    """Reconstruct a bundle against caller pins and exact accepted t_n packets.

    B7-family memo under the owner's 2026-08-31 rulings: a caller that
    already ran the FULL reconstruction this interval may thread the
    returned digest back as ``preverified_definition_digest``.  The pin
    branch re-derives the bundle's digest from its stored fields (a
    value-keyed walk — a mutated bundle re-hashes and mismatches), so
    tamper detection is preserved while the per-Picard-iteration full
    rebuild collapses.  Default ``None`` keeps the sealed behavior.
    """

    if (
        type(bundle) is not QSCMealMechanicalWorkDepositionBundle
        or bundle._seal is not _QSC_MEAL_DEPOSITION_BUNDLE_SEAL
    ):
        raise TypeError("meal-deposition validation requires a builder-issued exact bundle")
    QSCMealMechanicalWorkDepositionBundle.__post_init__(bundle)
    if preverified_definition_digest is not None:
        _require_sha256("preverified bundle definition", preverified_definition_digest)
        if bundle.definition_digest != preverified_definition_digest:
            raise QSCMealMechanicalWorkDepositionError(
                "meal-deposition bundle differs from its interval-preverified digest"
            )
        return preverified_definition_digest
    expected = _build_bundle(
        accepted_tray=accepted_tray,
        end_time_s=end_time_s,
        request=bundle.allocation_request,
        work_allocation=bundle.mechanical_work_allocation,
        context=validation_context,
    )
    if bundle != expected:
        raise QSCMealMechanicalWorkDepositionError(
            "meal-deposition bundle differs from exact C2/t_n reconstruction"
        )
    return bundle.definition_digest


def packet_q_other_s_roundoff_bound_j(
    row: QSCMealPacketMechanicalWorkDeposition,
    *,
    per_particle_energy_before_j: float,
    per_particle_energy_after_j: float,
) -> float:
    """Typed binary64 bound for divide/add/weight operations on one packet."""

    if type(row) is not QSCMealPacketMechanicalWorkDeposition:
        raise TypeError("Q_other_s roundoff requires an exact packet-deposition row")
    QSCMealPacketMechanicalWorkDeposition.__post_init__(row)
    _require_float(
        "Q_other_s packet energies",
        per_particle_energy_before_j,
        per_particle_energy_after_j,
    )
    weight = row.t_n_representative_particles
    per_particle = row.mechanical_work_per_particle_j
    return math.fsum(
        (
            weight
            * math.fsum(
                (
                    math.ulp(per_particle_energy_before_j),
                    math.ulp(per_particle),
                    math.ulp(per_particle_energy_after_j),
                )
            ),
            abs(per_particle) * math.ulp(weight),
            abs(weight) * math.ulp(per_particle),
            math.ulp(row.declared_mechanical_work_j),
        )
    )


__all__ = (
    "QSCMealMechanicalWorkDepositionBundle",
    "QSCMealMechanicalWorkDepositionError",
    "QSCMealMechanicalWorkDepositionLawAuthority",
    "QSCMealMechanicalWorkDepositionLedger",
    "QSCMealMechanicalWorkValidationContext",
    "QSCMealPacketMechanicalWorkDeposition",
    "QSC_MEAL_DEPOSITION_CLOSURE_LAW_ID",
    "QSC_MEAL_DEPOSITION_COVERAGE_LAW_ID",
    "QSC_MEAL_DEPOSITION_EXPECTED_LAW_DEFINITION_DIGEST",
    "QSC_MEAL_DEPOSITION_LAW_SCHEMA_ID",
    "QSC_MEAL_DEPOSITION_LAW_SCHEMA_REVISION",
    "QSC_MEAL_DEPOSITION_PACKET_BASIS_LAW_ID",
    "QSC_MEAL_DEPOSITION_PER_PARTICLE_LAW_ID",
    "QSC_MEAL_DEPOSITION_SCHEMA_ID",
    "QSC_MEAL_DEPOSITION_SCHEMA_REVISION",
    "QSC_MEAL_DEPOSITION_SLOT_DECLARATION",
    "QSC_MEAL_MECHANICAL_WORK_DEPOSITION_LAW_AUTHORITY",
    "build_qsc_meal_mechanical_work_deposition_bundle",
    "packet_q_other_s_roundoff_bound_j",
    "validate_qsc_meal_mechanical_work_deposition_bundle",
    "validate_qsc_meal_mechanical_work_deposition_law_authority",
)
