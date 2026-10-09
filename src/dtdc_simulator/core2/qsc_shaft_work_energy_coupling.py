"""Detached MECH-01C2 shaft-work allocation evidence for the QSC tower.

The source-pinned MECH-01A runtime authority owns the three D2 post-gearbox
power endpoints and both D3 allocation brackets.  This module does one bounded
piece of arithmetic beneath that authority: it converts one accepted interval's
declared mechanical work into either one input per physical tray wall or one
input per accepted K-layer packet.

Tray allocation uses current accepted weighted wet-mass sums; the bottom-heavy
bracket multiplies those sums by MECH-01A's pinned solid-order rank factors.
The first five tray values are direct binary64 products and the sixth is the
closing remainder.  For the all-to-meal endpoint, a tray's work is allocated by
accepted weighted wet mass.  Wet mass is exactly the sum of dry matter, attached
and internal hexane, external water, and retained water.  The residual-oil label
is already contained within the dry-matter material basis and is deliberately
excluded.

The common-energy ledger records both binary64 closure and the exact rational
sum of all emitted binary64 work values.  A disclosed ULP allowance recognizes
unavoidable multiplication and summation rounding; it is not a tunable physical
tolerance.  Validation reconstructs the complete allocation from an externally
pinned request digest, closing simple tamper and coherent replay-like request
substitution.

This module does not mutate a particle, wall, pressure-drop state, K identity,
or host ledger.  It introduces no second mixing carrier and supplies no finite
source temperature for mechanical work.  It is detached engineering allocation
evidence only: shaft-work wiring, QSC-10, physical qualification, plant
prediction, and industrial feasibility all remain false.
"""

from __future__ import annotations

import functools
import hashlib
import math
import re
import struct
from dataclasses import dataclass, replace as dataclass_replace
from fractions import Fraction
from typing import ClassVar, Iterable

from .packet_weight_scaling_contract import (
    PacketTotals,
    PacketWeight,
    PerParticleInventory,
    ScalingContractError,
    derive_packet_totals,
)
from .qsc_mechanical_runtime_authority import (
    QSC_MECH01_CONFIG_CANONICAL_LF_SHA256,
    QSC_MECH01_CONFIG_RELATIVE_PATH,
    QSC_MECH01_EXPECTED_COMMON_SHAFT_DEFINITION_DIGEST,
    QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST,
    QSC_MECH01_SOLID_TRAY_ORDER,
    QSCMechanicalRuntimeAuthority,
    ShaftWorkBoundaryPartition,
    WithinMealDistribution,
    make_engineering_runtime_authority,
    validate_qsc_mechanical_runtime_authority,
)
from .tray_particle_stateful import VerticalLayerId

QSC_SHAFT_WORK_SCHEMA_ID = "dtdc-core2-qsc-shaft-work-allocation-v1"
QSC_SHAFT_WORK_SCHEMA_REVISION = 1
QSC_SHAFT_WORK_ALLOCATOR_LAW_SCHEMA_ID = "dtdc-core2-qsc-shaft-work-allocator-law-v1"
QSC_SHAFT_WORK_ALLOCATOR_LAW_SCHEMA_REVISION = 1
QSC_SHAFT_WORK_EXPECTED_ALLOCATOR_LAW_DEFINITION_DIGEST = (
    "sha256:a32764e876f71416c239c74319d6765cf78de24bdf334fbba38fadda17ab3969"
)
_SHA256_PATTERN = re.compile(r"sha256:[0-9a-f]{64}\Z")
_RAW_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}\Z")
_TRAY_RANK = {tray_id: rank for rank, tray_id in enumerate(QSC_MECH01_SOLID_TRAY_ORDER)}


class QSCMechanicalAllocationError(ValueError):
    """A detached shaft-work request or allocation refused closed."""


def _negative_zero(value: float) -> bool:
    return struct.pack(">d", value) == b"\x80\x00\x00\x00\x00\x00\x00\x00"


def _require_float(
    name: str,
    *values: float,
    positive: bool = False,
    nonnegative: bool = False,
) -> None:
    if any(type(value) is not float for value in values):
        raise TypeError(f"{name} must contain exact binary64 values")
    if any(not math.isfinite(value) for value in values):
        raise QSCMechanicalAllocationError(f"{name} must contain finite values")
    if any(_negative_zero(value) for value in values):
        raise QSCMechanicalAllocationError(f"{name} rejects negative zero")
    if positive and any(value <= 0.0 for value in values):
        raise QSCMechanicalAllocationError(f"{name} must contain strictly positive values")
    if nonnegative and any(value < 0.0 for value in values):
        raise QSCMechanicalAllocationError(f"{name} must contain non-negative values")


def _require_identity(name: str, *values: str) -> None:
    if any(type(value) is not str for value in values):
        raise TypeError(f"{name} must contain exact strings")
    if any(not value or value != value.strip() for value in values):
        raise QSCMechanicalAllocationError(f"{name} must contain canonical nonblank strings")


def _require_sha256(name: str, value: str) -> None:
    if type(value) is not str:
        raise TypeError(f"{name} must be an exact string")
    if _SHA256_PATTERN.fullmatch(value) is None:
        raise QSCMechanicalAllocationError(f"{name} must be canonical lowercase sha256")


def _derive_packet_totals_typed(
    inventory: PerParticleInventory,
    weight: PacketWeight,
) -> PacketTotals:
    try:
        return derive_packet_totals(inventory, weight)
    except (ScalingContractError, OverflowError) as error:
        raise QSCMechanicalAllocationError(
            "packet-weight scaling produced a nonrepresentable extensive"
        ) from error


def _max_scaled_positive_weights(
    weights: tuple[float, ...],
    *,
    name: str,
) -> tuple[float, ...]:
    if type(weights) is not tuple or not weights:
        raise QSCMechanicalAllocationError(f"{name} requires a nonempty exact tuple")
    _require_float(name, *weights, positive=True)
    scale = max(weights)
    scaled = tuple(weight / scale for weight in weights)
    if any(value == 0.0 for value in scaled):
        raise QSCMechanicalAllocationError(
            f"{name} dynamic range underflows during finite max-rescaling"
        )
    _require_float(f"max-scaled {name}", *scaled, positive=True)
    return scaled


def _digest_part(value: str | int | float | bool) -> tuple[bytes, bytes]:
    if type(value) is str:
        return b"s", value.encode("utf-8")
    if type(value) is bool:
        return b"b", b"1" if value else b"0"
    if type(value) is int:
        return b"i", str(value).encode("ascii")
    if type(value) is float:
        return b"f", struct.pack(">d", value)
    raise TypeError(f"unsupported shaft-work digest part: {type(value).__name__}")


def _digest(domain: str, parts: Iterable[str | int | float | bool]) -> str:
    digest = hashlib.sha256()
    digest.update(b"DTDC-QSC-SHAFT-WORK-ALLOCATION-V1\0")
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


def _canonical_fraction(numerator: int, denominator: int, *, name: str) -> Fraction:
    if type(numerator) is not int or type(denominator) is not int:
        raise TypeError(f"{name} ratio must contain exact integers")
    if denominator <= 0:
        raise QSCMechanicalAllocationError(f"{name} denominator must be positive")
    value = Fraction(numerator, denominator)
    if value.numerator != numerator or value.denominator != denominator:
        raise QSCMechanicalAllocationError(f"{name} ratio must be in canonical lowest terms")
    return value


def _fraction_parts(value: Fraction) -> tuple[int, int]:
    return value.numerator, value.denominator


def _within_disclosed_roundoff(
    observed: Iterable[float],
    target: float,
) -> bool:
    values = tuple(observed)
    exact_observed = sum((Fraction.from_float(value) for value in values), Fraction())
    exact_target = Fraction.from_float(target)
    allowance = (2 * len(values) + 2) * Fraction.from_float(math.ulp(target))
    return abs(exact_observed - exact_target) <= allowance


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCShaftWorkAllocatorLawAuthority:
    """Source-pinned D2/D3 arithmetic law, separate from accepted state."""

    schema_id: str
    schema_revision: int
    source_runtime_definition_digest: str
    source_common_shaft_definition_digest: str
    source_configuration_relative_path: str
    source_configuration_canonical_lf_sha256: str
    solid_tray_order: tuple[str, ...]
    d2_power_endpoints_w: tuple[float, ...]
    bottom_heavy_rank_multipliers: tuple[float, ...]
    weighted_wet_mass_total_indices: tuple[int, ...]
    weighted_wet_mass_component_names: tuple[str, ...]
    excluded_residual_oil_label_name: str
    tray_fraction_law_id: str
    bottom_heavy_law_id: str
    within_tray_meal_law_id: str
    normalization_law_id: str
    tray_closing_law_id: str
    packet_closing_law_id: str
    exact_ledger_law_id: str
    roundoff_ulp_coefficient: int
    roundoff_ulp_intercept: int
    #: Ruled 2026-08-30 (GT_PS2_V3_EXTENSION_AND_RPM_AUTHORITY_RULING): the
    #: declared constant operating speed of the bound runtime source.  A
    #: projection cross-checked against the source digest below — not
    #: independent state, so it does not enter the definition digest.
    source_declared_constant_shaft_rpm: float = 11.0

    declared_engineering_law: ClassVar[bool] = True
    authenticated_machine_data: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self.schema_id != QSC_SHAFT_WORK_ALLOCATOR_LAW_SCHEMA_ID:
            raise QSCMechanicalAllocationError("allocator law has the wrong schema identity")
        if self.schema_revision != QSC_SHAFT_WORK_ALLOCATOR_LAW_SCHEMA_REVISION:
            raise QSCMechanicalAllocationError("allocator law has the wrong schema revision")
        if self.source_declared_constant_shaft_rpm == 11.0:
            if self.source_runtime_definition_digest != (
                QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST
            ):
                raise QSCMechanicalAllocationError("allocator law binds a foreign runtime source")
        elif self.source_runtime_definition_digest != (
            make_engineering_runtime_authority(
                self.source_declared_constant_shaft_rpm
            ).definition_digest
        ):
            raise QSCMechanicalAllocationError(
                "allocator law binds a foreign runtime source at its declared speed"
            )
        if self.source_common_shaft_definition_digest != (
            QSC_MECH01_EXPECTED_COMMON_SHAFT_DEFINITION_DIGEST
        ):
            raise QSCMechanicalAllocationError("allocator law binds a foreign common shaft")
        if self.source_configuration_relative_path != QSC_MECH01_CONFIG_RELATIVE_PATH:
            raise QSCMechanicalAllocationError("allocator law binds a foreign source path")
        if (
            type(self.source_configuration_canonical_lf_sha256) is not str
            or _RAW_SHA256_PATTERN.fullmatch(self.source_configuration_canonical_lf_sha256) is None
            or self.source_configuration_canonical_lf_sha256
            != QSC_MECH01_CONFIG_CANONICAL_LF_SHA256
        ):
            raise QSCMechanicalAllocationError("allocator law binds foreign source bytes")
        if type(self.solid_tray_order) is not tuple or self.solid_tray_order != (
            "PD1",
            "PD2",
            "PD3",
            "MN1",
            "MN2",
            "SP1",
        ):
            raise QSCMechanicalAllocationError("allocator law has a foreign solid tray order")
        if type(self.d2_power_endpoints_w) is not tuple or self.d2_power_endpoints_w != (
            25_000.0,
            50_000.0,
            165_000.0,
        ):
            raise QSCMechanicalAllocationError("allocator law has foreign D2 endpoints")
        if type(
            self.bottom_heavy_rank_multipliers
        ) is not tuple or self.bottom_heavy_rank_multipliers != (1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
            raise QSCMechanicalAllocationError("allocator law has foreign D3 rank multipliers")
        if (
            type(self.weighted_wet_mass_total_indices) is not tuple
            or self.weighted_wet_mass_total_indices != (0, 2, 3, 4, 5)
            or any(type(value) is not int for value in self.weighted_wet_mass_total_indices)
        ):
            raise QSCMechanicalAllocationError("allocator law has a foreign wet-mass mapping")
        if (
            type(self.weighted_wet_mass_component_names) is not tuple
            or self.weighted_wet_mass_component_names
            != (
                "dry_matter_kg",
                "attached_hexane_kg",
                "internal_hexane_kg",
                "external_water_kg",
                "retained_water_kg",
            )
            or self.excluded_residual_oil_label_name != "residual_oil_label_kg"
        ):
            raise QSCMechanicalAllocationError("allocator law has foreign wet-mass semantics")
        expected_law_ids = (
            "CURRENT_ACCEPTED_WEIGHTED_WET_MASS_BY_PHYSICAL_TRAY",
            "CURRENT_ACCEPTED_WEIGHTED_WET_MASS_TIMES_SOLID_ORDER_RANK",
            "CURRENT_ACCEPTED_WEIGHTED_WET_MASS_BY_K_PACKET",
            "FINITE_MAX_RESCALED_POSITIVE_WEIGHTS_V1",
            "FIRST_FIVE_PRODUCTS_SIXTH_BINARY64_REMAINDER",
            "FIRST_N_MINUS_ONE_PRODUCTS_FINAL_BINARY64_REMAINDER",
            "EXACT_RATIONAL_EMITTED_BINARY64_COMMON_ENERGY_LEDGER",
        )
        observed_law_ids = (
            self.tray_fraction_law_id,
            self.bottom_heavy_law_id,
            self.within_tray_meal_law_id,
            self.normalization_law_id,
            self.tray_closing_law_id,
            self.packet_closing_law_id,
            self.exact_ledger_law_id,
        )
        _require_identity("allocator law identity", *observed_law_ids)
        if observed_law_ids != expected_law_ids:
            raise QSCMechanicalAllocationError("allocator law identifiers drifted")
        if (
            type(self.roundoff_ulp_coefficient) is not int
            or type(self.roundoff_ulp_intercept) is not int
            or (self.roundoff_ulp_coefficient, self.roundoff_ulp_intercept) != (2, 2)
        ):
            raise QSCMechanicalAllocationError("allocator law roundoff formula drifted")

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-shaft-work-allocator-law-authority",
            (
                self.schema_id,
                self.schema_revision,
                self.source_runtime_definition_digest,
                self.source_common_shaft_definition_digest,
                self.source_configuration_relative_path,
                self.source_configuration_canonical_lf_sha256,
                *self.solid_tray_order,
                *self.d2_power_endpoints_w,
                *self.bottom_heavy_rank_multipliers,
                *self.weighted_wet_mass_total_indices,
                *self.weighted_wet_mass_component_names,
                self.excluded_residual_oil_label_name,
                self.tray_fraction_law_id,
                self.bottom_heavy_law_id,
                self.within_tray_meal_law_id,
                self.normalization_law_id,
                self.tray_closing_law_id,
                self.packet_closing_law_id,
                self.exact_ledger_law_id,
                self.roundoff_ulp_coefficient,
                self.roundoff_ulp_intercept,
            ),
        )


QSC_SHAFT_WORK_ALLOCATOR_LAW_AUTHORITY = QSCShaftWorkAllocatorLawAuthority(
    schema_id=QSC_SHAFT_WORK_ALLOCATOR_LAW_SCHEMA_ID,
    schema_revision=QSC_SHAFT_WORK_ALLOCATOR_LAW_SCHEMA_REVISION,
    source_runtime_definition_digest=QSC_MECH01_EXPECTED_RUNTIME_DEFINITION_DIGEST,
    source_common_shaft_definition_digest=QSC_MECH01_EXPECTED_COMMON_SHAFT_DEFINITION_DIGEST,
    source_configuration_relative_path=QSC_MECH01_CONFIG_RELATIVE_PATH,
    source_configuration_canonical_lf_sha256=QSC_MECH01_CONFIG_CANONICAL_LF_SHA256,
    solid_tray_order=QSC_MECH01_SOLID_TRAY_ORDER,
    d2_power_endpoints_w=(25_000.0, 50_000.0, 165_000.0),
    bottom_heavy_rank_multipliers=(1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
    weighted_wet_mass_total_indices=(0, 2, 3, 4, 5),
    weighted_wet_mass_component_names=(
        "dry_matter_kg",
        "attached_hexane_kg",
        "internal_hexane_kg",
        "external_water_kg",
        "retained_water_kg",
    ),
    excluded_residual_oil_label_name="residual_oil_label_kg",
    tray_fraction_law_id="CURRENT_ACCEPTED_WEIGHTED_WET_MASS_BY_PHYSICAL_TRAY",
    bottom_heavy_law_id="CURRENT_ACCEPTED_WEIGHTED_WET_MASS_TIMES_SOLID_ORDER_RANK",
    within_tray_meal_law_id="CURRENT_ACCEPTED_WEIGHTED_WET_MASS_BY_K_PACKET",
    normalization_law_id="FINITE_MAX_RESCALED_POSITIVE_WEIGHTS_V1",
    tray_closing_law_id="FIRST_FIVE_PRODUCTS_SIXTH_BINARY64_REMAINDER",
    packet_closing_law_id="FIRST_N_MINUS_ONE_PRODUCTS_FINAL_BINARY64_REMAINDER",
    exact_ledger_law_id="EXACT_RATIONAL_EMITTED_BINARY64_COMMON_ENERGY_LEDGER",
    roundoff_ulp_coefficient=2,
    roundoff_ulp_intercept=2,
)


@functools.lru_cache(maxsize=8)
def make_engineering_allocator_law_authority(
    shaft_rpm: float = 11.0,
) -> QSCShaftWorkAllocatorLawAuthority:
    """The D2/D3 allocator law bound to a ruled constant-speed runtime variant.

    Ruled 2026-08-30 (GT_PS2_V3_EXTENSION_AND_RPM_AUTHORITY_RULING): the
    arithmetic law is speed-invariant; only the bound runtime source digest
    changes, so a variant law validates and digests against the exact
    config-built runtime at its declared speed.
    """

    if shaft_rpm == 11.0:
        return QSC_SHAFT_WORK_ALLOCATOR_LAW_AUTHORITY
    return dataclass_replace(
        QSC_SHAFT_WORK_ALLOCATOR_LAW_AUTHORITY,
        source_runtime_definition_digest=(
            make_engineering_runtime_authority(shaft_rpm).definition_digest
        ),
        source_declared_constant_shaft_rpm=shaft_rpm,
    )


# Public compatibility views are checked against the typed authority before
# every build and validation.  Monkeypatching either view therefore refuses;
# neither tuple is an independent operative authority.
QSC_SHAFT_WORK_POWER_ENDPOINTS_W = QSC_SHAFT_WORK_ALLOCATOR_LAW_AUTHORITY.d2_power_endpoints_w
QSC_SHAFT_WORK_BOTTOM_HEAVY_RANK_MULTIPLIERS = (
    QSC_SHAFT_WORK_ALLOCATOR_LAW_AUTHORITY.bottom_heavy_rank_multipliers
)


def validate_qsc_shaft_work_allocator_law_authority(
    authority: QSCShaftWorkAllocatorLawAuthority,
) -> str:
    """Accept only the exact source-pinned MECH-01C2 allocator law."""

    if type(authority) is not QSCShaftWorkAllocatorLawAuthority:
        raise TypeError("allocator law must be an exact authority record")
    QSCShaftWorkAllocatorLawAuthority.__post_init__(authority)
    digest = authority.definition_digest
    expected = (
        QSC_SHAFT_WORK_EXPECTED_ALLOCATOR_LAW_DEFINITION_DIGEST
        if authority.source_declared_constant_shaft_rpm == 11.0
        else make_engineering_allocator_law_authority(
            authority.source_declared_constant_shaft_rpm
        ).definition_digest
    )
    if digest != expected:
        raise QSCMechanicalAllocationError("allocator law differs from its source-pinned digest")
    if QSC_SHAFT_WORK_POWER_ENDPOINTS_W != authority.d2_power_endpoints_w:
        raise QSCMechanicalAllocationError("runtime D2 endpoint view was substituted")
    if QSC_SHAFT_WORK_BOTTOM_HEAVY_RANK_MULTIPLIERS != (authority.bottom_heavy_rank_multipliers):
        raise QSCMechanicalAllocationError("runtime D3 rank-multiplier view was substituted")
    if QSC_MECH01_SOLID_TRAY_ORDER != authority.solid_tray_order or _TRAY_RANK != {
        tray_id: rank for rank, tray_id in enumerate(authority.solid_tray_order)
    }:
        raise QSCMechanicalAllocationError("runtime solid-tray order view was substituted")
    return digest


#: Delegated-authority memo store for AcceptedQSCMechanicalPacket digests
#: (owner grant 2026-08-31): keyed on the exact identity + content bytes.
_PACKET_DIGEST_CACHE: dict[tuple, str] = {}
_PACKET_DIGEST_CACHE_LIMIT = 500_000


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedQSCMechanicalPacket:
    """Accepted packet identity and exact external weight used by MECH-01C2."""

    physical_tray_id: str
    vertical_layer_id: VerticalLayerId
    packet_id: str
    particle_state_identity_digest: str
    inventory: PerParticleInventory
    weight: PacketWeight

    detached_engineering_evidence: ClassVar[bool] = True
    mutates_particle_state: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("packet identity", self.physical_tray_id, self.packet_id)
        if self.physical_tray_id not in _TRAY_RANK:
            raise QSCMechanicalAllocationError("packet names a foreign physical tray")
        if type(self.vertical_layer_id) is not VerticalLayerId:
            raise TypeError("packet vertical layer must be an exact VerticalLayerId")
        if type(self.vertical_layer_id.value) is not int or self.vertical_layer_id.value < 1:
            raise QSCMechanicalAllocationError(
                "MECH-01C2 packet layer must be an exact positive accepted K layer"
            )
        _require_sha256("particle-state identity", self.particle_state_identity_digest)
        if type(self.inventory) is not PerParticleInventory:
            raise TypeError("packet inventory must be an exact PerParticleInventory")
        if type(self.weight) is not PacketWeight:
            raise TypeError("packet weight must be an exact PacketWeight")
        # Re-run the exact types' validators so post-construction hostile mutation
        # cannot smuggle negative or nonfinite masses into an allocation.
        try:
            PerParticleInventory.__post_init__(self.inventory)
            PacketWeight.__post_init__(self.weight)
        except ScalingContractError as error:
            raise QSCMechanicalAllocationError(
                "packet inventory or weight failed exact scaling-contract revalidation"
            ) from error
        if self.inventory.energy_datum_id != self.weight.energy_datum_id:
            raise QSCMechanicalAllocationError("packet inventory and weight bind foreign datums")
        if self.accepted_weighted_wet_mass_kg <= 0.0:
            raise QSCMechanicalAllocationError("packet accepted weighted wet mass must be positive")

    @property
    def accepted_weighted_wet_mass_kg(self) -> float:
        totals = _derive_packet_totals_typed(self.inventory, self.weight).totals
        try:
            # The typed allocator law binds the same exact mapping.  Keeping the
            # tuple literal here prevents a mutable module alias from becoming a
            # second authority during request construction.
            wet_mass = math.fsum(totals[index] for index in (0, 2, 3, 4, 5))
        except OverflowError as error:
            raise QSCMechanicalAllocationError(
                "packet accepted weighted wet mass is not representable"
            ) from error
        _require_float("packet accepted weighted wet mass", wet_mass, positive=True)
        return wet_mass

    @property
    def definition_digest(self) -> str:
        # Delegated-authority memo (owner grant 2026-08-31, epsilon-rerule
        # record; TEMPORARY-PENDING-OWNER-REVIEW): value-keyed on the exact
        # content bytes - the same round-3/B7 pattern - so a mutated packet
        # re-keys and recomputes.  Evidence: bit-identical march A/B in the
        # landing commit.
        key = (
            self.physical_tray_id,
            self.vertical_layer_id.value,
            self.packet_id,
            self.particle_state_identity_digest,
            self.inventory.energy_datum_id,
            self.weight.energy_datum_id,
            struct.pack(
                ">8d",
                *self.inventory.as_tuple(),
                self.weight.representative_particles,
            ),
        )
        hit = _PACKET_DIGEST_CACHE.get(key)
        if hit is not None:
            return hit
        totals = _derive_packet_totals_typed(self.inventory, self.weight)
        value = _digest(
            "accepted-qsc-mechanical-packet",
            (
                self.physical_tray_id,
                self.vertical_layer_id.value,
                self.packet_id,
                self.particle_state_identity_digest,
                self.inventory.energy_datum_id,
                *self.inventory.as_tuple(),
                self.weight.representative_particles,
                totals.contract_digest,
                self.accepted_weighted_wet_mass_kg,
            ),
        )
        if len(_PACKET_DIGEST_CACHE) < _PACKET_DIGEST_CACHE_LIMIT:
            _PACKET_DIGEST_CACHE[key] = value
        return value


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMechanicalAllocationRequest:
    """Externally pinned accepted-step request for detached allocation."""

    accepted_step_identity_digest: str
    allocator_law_definition_digest: str
    energy_datum_id: str
    dt_s: float
    shaft_power_w: float
    boundary_partition: ShaftWorkBoundaryPartition
    tray_distribution: WithinMealDistribution
    packets: tuple[AcceptedQSCMechanicalPacket, ...]

    detached_engineering_evidence: ClassVar[bool] = True
    persistent_replay_protection_implemented: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_sha256("accepted-step identity", self.accepted_step_identity_digest)
        _require_sha256("allocator-law definition", self.allocator_law_definition_digest)
        _require_identity("common energy datum", self.energy_datum_id)
        _require_float("accepted-step duration", self.dt_s, positive=True)
        _require_float("D2 shaft power", self.shaft_power_w, positive=True)
        # B4 (owner-ruled GO 2026-08-31): the D2 Coulombic law is linear in
        # shaft speed, so the ruled bracket is a CONTINUUM admission; the
        # three endpoint values remain law data pinned in the allocator-law
        # authority.  Every previously admitted request is unchanged.
        if not 25_000.0 <= self.shaft_power_w <= 165_000.0:
            raise QSCMechanicalAllocationError("shaft power left the ruled D2 bracket [25, 165] kW")
        if type(self.boundary_partition) is not ShaftWorkBoundaryPartition:
            raise TypeError("boundary partition must be an exact MECH-01A enum")
        if type(self.tray_distribution) is not WithinMealDistribution:
            raise TypeError("tray distribution must be an exact MECH-01A enum")
        if type(self.packets) is not tuple or not self.packets:
            raise QSCMechanicalAllocationError("allocation requires a nonempty exact packet tuple")
        if any(type(packet) is not AcceptedQSCMechanicalPacket for packet in self.packets):
            raise TypeError("allocation request contains a foreign packet type")
        for packet in self.packets:
            AcceptedQSCMechanicalPacket.__post_init__(packet)
            if packet.inventory.energy_datum_id != self.energy_datum_id:
                raise QSCMechanicalAllocationError("packet binds a foreign common energy datum")

        key = lambda packet: (  # noqa: E731 - local canonical key is clearer inline
            _TRAY_RANK[packet.physical_tray_id],
            packet.vertical_layer_id.value,
            packet.packet_id,
        )
        if self.packets != tuple(sorted(self.packets, key=key)):
            raise QSCMechanicalAllocationError("packets must use canonical tray/K/packet order")
        identities = tuple(
            (packet.physical_tray_id, packet.vertical_layer_id.value, packet.packet_id)
            for packet in self.packets
        )
        if len(set(identities)) != len(identities):
            raise QSCMechanicalAllocationError("allocation request repeats a packet identity")

        packet_masses = tuple(packet.accepted_weighted_wet_mass_kg for packet in self.packets)
        _max_scaled_positive_weights(
            packet_masses,
            name="request accepted weighted wet masses",
        )
        for tray_id in QSC_MECH01_SOLID_TRAY_ORDER:
            masses = tuple(
                packet.accepted_weighted_wet_mass_kg
                for packet in self.packets
                if packet.physical_tray_id == tray_id
            )
            if not masses:
                raise QSCMechanicalAllocationError(
                    "every physical tray requires at least one accepted packet"
                )

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-mechanical-allocation-request",
            (
                self.accepted_step_identity_digest,
                self.allocator_law_definition_digest,
                self.energy_datum_id,
                self.dt_s,
                self.shaft_power_w,
                self.boundary_partition.value,
                self.tray_distribution.value,
                *(packet.definition_digest for packet in self.packets),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCTrayMechanicalWorkAllocation:
    physical_tray_id: str
    distribution_fraction: float
    declared_mechanical_work_j: float
    closing_sixth_remainder: bool

    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("tray-work identity", self.physical_tray_id)
        if self.physical_tray_id not in _TRAY_RANK:
            raise QSCMechanicalAllocationError("tray-work row names a foreign tray")
        _require_float("tray-work fraction", self.distribution_fraction, positive=True)
        _require_float("tray mechanical work", self.declared_mechanical_work_j, positive=True)
        if self.distribution_fraction > 1.0:
            raise QSCMechanicalAllocationError("tray-work fraction must not exceed one")
        if type(self.closing_sixth_remainder) is not bool:
            raise TypeError("closing-remainder marker must be an exact bool")

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-tray-mechanical-work",
            (
                self.physical_tray_id,
                self.distribution_fraction,
                self.declared_mechanical_work_j,
                self.closing_sixth_remainder,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class MealPacketMechanicalWorkInput:
    physical_tray_id: str
    vertical_layer_id: VerticalLayerId
    packet_id: str
    packet_definition_digest: str
    accepted_weighted_wet_mass_kg: float
    declared_mechanical_work_j: float
    closing_packet_remainder: bool

    is_heat_source: ClassVar[bool] = False
    source_temperature_defined: ClassVar[bool] = False
    mutates_particle_state: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("meal-work packet identity", self.physical_tray_id, self.packet_id)
        if self.physical_tray_id not in _TRAY_RANK:
            raise QSCMechanicalAllocationError("meal-work row names a foreign tray")
        if type(self.vertical_layer_id) is not VerticalLayerId:
            raise TypeError("meal-work layer must be an exact VerticalLayerId")
        if type(self.vertical_layer_id.value) is not int or self.vertical_layer_id.value < 1:
            raise QSCMechanicalAllocationError(
                "meal-work layer must be an exact positive accepted K layer"
            )
        _require_sha256("meal-work packet definition", self.packet_definition_digest)
        _require_float("meal-work wet mass", self.accepted_weighted_wet_mass_kg, positive=True)
        _require_float("meal mechanical work", self.declared_mechanical_work_j, positive=True)
        if type(self.closing_packet_remainder) is not bool:
            raise TypeError("packet closing-remainder marker must be an exact bool")

    @property
    def definition_digest(self) -> str:
        return _digest(
            "meal-packet-mechanical-work-input",
            (
                self.physical_tray_id,
                self.vertical_layer_id.value,
                self.packet_id,
                self.packet_definition_digest,
                self.accepted_weighted_wet_mass_kg,
                self.declared_mechanical_work_j,
                self.closing_packet_remainder,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PhysicalTrayWallMechanicalWorkInput:
    physical_tray_id: str
    declared_mechanical_work_j: float

    is_heat_source: ClassVar[bool] = False
    source_temperature_defined: ClassVar[bool] = False
    one_input_per_physical_tray: ClassVar[bool] = True
    mutates_wall_state: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("wall-work tray identity", self.physical_tray_id)
        if self.physical_tray_id not in _TRAY_RANK:
            raise QSCMechanicalAllocationError("wall-work row names a foreign tray")
        _require_float("wall mechanical work", self.declared_mechanical_work_j, positive=True)

    @property
    def definition_digest(self) -> str:
        return _digest(
            "physical-tray-wall-mechanical-work-input",
            (self.physical_tray_id, self.declared_mechanical_work_j),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCCommonEnergyAllocationLedger:
    energy_datum_id: str
    runtime_definition_digest: str
    allocator_law_definition_digest: str
    accepted_step_identity_digest: str
    request_definition_digest: str
    boundary_partition: ShaftWorkBoundaryPartition
    shaft_power_w: float
    dt_s: float
    declared_mechanical_work_j: float
    meal_mechanical_work_j: float
    wall_mechanical_work_j: float
    allocated_mechanical_work_j: float
    binary64_closure_residual_j: float
    exact_declared_numerator: int
    exact_declared_denominator: int
    exact_allocated_numerator: int
    exact_allocated_denominator: int
    exact_residual_numerator: int
    exact_residual_denominator: int
    closure_ulp_limit: int
    authoritative_conservation: bool

    exact_common_energy_ledger: ClassVar[bool] = True
    roundoff_aware: ClassVar[bool] = True
    physically_qualifying: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("ledger energy datum", self.energy_datum_id)
        for name in (
            "runtime_definition_digest",
            "allocator_law_definition_digest",
            "accepted_step_identity_digest",
            "request_definition_digest",
        ):
            _require_sha256(name, getattr(self, name))
        if type(self.boundary_partition) is not ShaftWorkBoundaryPartition:
            raise TypeError("ledger partition must be an exact MECH-01A enum")
        _require_float("ledger shaft power", self.shaft_power_w, positive=True)
        _require_float("ledger duration", self.dt_s, positive=True)
        _require_float("ledger declared work", self.declared_mechanical_work_j, positive=True)
        _require_float(
            "ledger allocated work",
            self.meal_mechanical_work_j,
            self.wall_mechanical_work_j,
            self.allocated_mechanical_work_j,
            nonnegative=True,
        )
        _require_float("ledger closure residual", self.binary64_closure_residual_j)
        if self.shaft_power_w * self.dt_s != self.declared_mechanical_work_j:
            raise QSCMechanicalAllocationError("ledger declared work differs from power times dt")
        if math.fsum((self.meal_mechanical_work_j, self.wall_mechanical_work_j)) != (
            self.allocated_mechanical_work_j
        ):
            raise QSCMechanicalAllocationError("ledger partition totals do not reconstruct")
        if self.binary64_closure_residual_j != (
            self.allocated_mechanical_work_j - self.declared_mechanical_work_j
        ):
            raise QSCMechanicalAllocationError("ledger binary64 residual does not reconstruct")
        if self.boundary_partition is ShaftWorkBoundaryPartition.ALL_TO_MEAL:
            if self.wall_mechanical_work_j != 0.0:
                raise QSCMechanicalAllocationError("all-to-meal ledger contains wall work")
        elif self.meal_mechanical_work_j != 0.0:
            raise QSCMechanicalAllocationError("all-to-wall ledger contains meal work")

        declared = _canonical_fraction(
            self.exact_declared_numerator,
            self.exact_declared_denominator,
            name="exact declared work",
        )
        allocated = _canonical_fraction(
            self.exact_allocated_numerator,
            self.exact_allocated_denominator,
            name="exact allocated work",
        )
        residual = _canonical_fraction(
            self.exact_residual_numerator,
            self.exact_residual_denominator,
            name="exact closure residual",
        )
        if declared != Fraction.from_float(self.declared_mechanical_work_j):
            raise QSCMechanicalAllocationError("exact declared-work ratio differs from binary64")
        if residual != allocated - declared:
            raise QSCMechanicalAllocationError("exact closure-residual ratio does not reconstruct")
        if type(self.closure_ulp_limit) is not int or self.closure_ulp_limit <= 0:
            raise QSCMechanicalAllocationError("closure ULP limit must be a positive exact integer")
        ulp = Fraction.from_float(math.ulp(self.declared_mechanical_work_j))
        if abs(residual) > self.closure_ulp_limit * ulp:
            raise QSCMechanicalAllocationError("mechanical-work closure exceeds disclosed roundoff")
        if self.authoritative_conservation is not True:
            raise QSCMechanicalAllocationError("common-energy conservation must be authoritative")

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-common-energy-allocation-ledger",
            (
                self.energy_datum_id,
                self.runtime_definition_digest,
                self.allocator_law_definition_digest,
                self.accepted_step_identity_digest,
                self.request_definition_digest,
                self.boundary_partition.value,
                self.shaft_power_w,
                self.dt_s,
                self.declared_mechanical_work_j,
                self.meal_mechanical_work_j,
                self.wall_mechanical_work_j,
                self.allocated_mechanical_work_j,
                self.binary64_closure_residual_j,
                self.exact_declared_numerator,
                self.exact_declared_denominator,
                self.exact_allocated_numerator,
                self.exact_allocated_denominator,
                self.exact_residual_numerator,
                self.exact_residual_denominator,
                self.closure_ulp_limit,
                self.authoritative_conservation,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class QSCMechanicalWorkAllocation:
    schema_id: str
    schema_revision: int
    runtime_definition_digest: str
    allocator_law_definition_digest: str
    common_shaft_definition_digest: str
    accepted_history_definition_digest: str
    request_definition_digest: str
    tray_allocations: tuple[QSCTrayMechanicalWorkAllocation, ...]
    meal_packet_inputs: tuple[MealPacketMechanicalWorkInput, ...]
    wall_inputs: tuple[PhysicalTrayWallMechanicalWorkInput, ...]
    ledger: QSCCommonEnergyAllocationLedger

    detached_engineering_allocation_evidence: ClassVar[bool] = True
    shaft_work_wired: ClassVar[bool] = False
    mutates_wall_state: ClassVar[bool] = False
    mutates_pressure_drop_state: ClassVar[bool] = False
    rekeys_k: ClassVar[bool] = False
    independent_mixing_carrier_introduced: ClassVar[bool] = False
    qsc_10: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    industry_feasible_operating_envelope: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self.schema_id != QSC_SHAFT_WORK_SCHEMA_ID:
            raise QSCMechanicalAllocationError("allocation has the wrong schema identity")
        if self.schema_revision != QSC_SHAFT_WORK_SCHEMA_REVISION:
            raise QSCMechanicalAllocationError("allocation has the wrong schema revision")
        for name in (
            "runtime_definition_digest",
            "allocator_law_definition_digest",
            "common_shaft_definition_digest",
            "accepted_history_definition_digest",
            "request_definition_digest",
        ):
            _require_sha256(name, getattr(self, name))
        if type(self.tray_allocations) is not tuple or len(self.tray_allocations) != 6:
            raise QSCMechanicalAllocationError("allocation requires six exact tray rows")
        if any(type(row) is not QSCTrayMechanicalWorkAllocation for row in self.tray_allocations):
            raise TypeError("allocation contains a foreign tray-work row")
        for row in self.tray_allocations:
            QSCTrayMechanicalWorkAllocation.__post_init__(row)
        if tuple(row.physical_tray_id for row in self.tray_allocations) != (
            QSC_MECH01_SOLID_TRAY_ORDER
        ):
            raise QSCMechanicalAllocationError("tray-work rows differ from solid tray order")
        if tuple(row.closing_sixth_remainder for row in self.tray_allocations) != (
            False,
            False,
            False,
            False,
            False,
            True,
        ):
            raise QSCMechanicalAllocationError("only the sixth tray may be the closing remainder")

        if type(self.meal_packet_inputs) is not tuple or type(self.wall_inputs) is not tuple:
            raise TypeError("allocation outputs must be exact tuples")
        if any(type(row) is not MealPacketMechanicalWorkInput for row in self.meal_packet_inputs):
            raise TypeError("allocation contains a foreign meal-work row")
        if any(type(row) is not PhysicalTrayWallMechanicalWorkInput for row in self.wall_inputs):
            raise TypeError("allocation contains a foreign wall-work row")
        for row in self.meal_packet_inputs:
            MealPacketMechanicalWorkInput.__post_init__(row)
        for row in self.wall_inputs:
            PhysicalTrayWallMechanicalWorkInput.__post_init__(row)
        if type(self.ledger) is not QSCCommonEnergyAllocationLedger:
            raise TypeError("allocation ledger has a foreign exact type")
        QSCCommonEnergyAllocationLedger.__post_init__(self.ledger)
        if (
            self.ledger.runtime_definition_digest != self.runtime_definition_digest
            or self.ledger.allocator_law_definition_digest != self.allocator_law_definition_digest
            or self.ledger.request_definition_digest != self.request_definition_digest
        ):
            raise QSCMechanicalAllocationError("ledger binds a foreign runtime, law, or request")

        if not _within_disclosed_roundoff(
            (row.declared_mechanical_work_j for row in self.tray_allocations),
            self.ledger.declared_mechanical_work_j,
        ):
            raise QSCMechanicalAllocationError("tray rows do not close to declared shaft work")
        if self.ledger.boundary_partition is ShaftWorkBoundaryPartition.ALL_TO_MEAL:
            if not self.meal_packet_inputs or self.wall_inputs:
                raise QSCMechanicalAllocationError("all-to-meal output shape is not exclusive")
            for tray_row in self.tray_allocations:
                rows = tuple(
                    row
                    for row in self.meal_packet_inputs
                    if row.physical_tray_id == tray_row.physical_tray_id
                )
                if not rows or sum(row.closing_packet_remainder for row in rows) != 1:
                    raise QSCMechanicalAllocationError(
                        "each meal tray requires exactly one closing packet remainder"
                    )
                if not rows[-1].closing_packet_remainder:
                    raise QSCMechanicalAllocationError("only the final tray packet may close")
                if not _within_disclosed_roundoff(
                    (row.declared_mechanical_work_j for row in rows),
                    tray_row.declared_mechanical_work_j,
                ):
                    raise QSCMechanicalAllocationError("meal packet rows do not close by tray")
        else:
            if self.meal_packet_inputs or len(self.wall_inputs) != 6:
                raise QSCMechanicalAllocationError("all-to-wall output shape is not exclusive")
            if tuple(row.physical_tray_id for row in self.wall_inputs) != (
                QSC_MECH01_SOLID_TRAY_ORDER
            ):
                raise QSCMechanicalAllocationError("wall inputs differ from solid tray order")
            if any(
                wall.declared_mechanical_work_j != tray.declared_mechanical_work_j
                for wall, tray in zip(self.wall_inputs, self.tray_allocations, strict=True)
            ):
                raise QSCMechanicalAllocationError("wall inputs differ from tray allocations")

        emitted = tuple(
            row.declared_mechanical_work_j for row in (*self.meal_packet_inputs, *self.wall_inputs)
        )
        if math.fsum(emitted) != self.ledger.allocated_mechanical_work_j:
            raise QSCMechanicalAllocationError("emitted work differs from common-energy ledger")
        exact_emitted = sum((Fraction.from_float(value) for value in emitted), Fraction())
        if exact_emitted != Fraction(
            self.ledger.exact_allocated_numerator,
            self.ledger.exact_allocated_denominator,
        ):
            raise QSCMechanicalAllocationError("exact emitted work differs from ledger")

    @property
    def definition_digest(self) -> str:
        return _digest(
            "qsc-mechanical-work-allocation",
            (
                self.schema_id,
                self.schema_revision,
                self.runtime_definition_digest,
                self.allocator_law_definition_digest,
                self.common_shaft_definition_digest,
                self.accepted_history_definition_digest,
                self.request_definition_digest,
                *(row.definition_digest for row in self.tray_allocations),
                *(row.definition_digest for row in self.meal_packet_inputs),
                *(row.definition_digest for row in self.wall_inputs),
                self.ledger.definition_digest,
            ),
        )


def _closed_positive_allocations(total: float, fractions: tuple[float, ...]) -> tuple[float, ...]:
    _require_float("allocation total", total, positive=True)
    _require_float("allocation fractions", *fractions, positive=True)
    if not fractions or not _within_disclosed_roundoff(fractions, 1.0):
        raise QSCMechanicalAllocationError("allocation fractions exceed disclosed roundoff")
    values = tuple(total * fraction for fraction in fractions[:-1])
    closing = total - math.fsum(values)
    result = (*values, closing)
    _require_float("closed allocation", *result, positive=True)
    if not _within_disclosed_roundoff(result, total):
        raise QSCMechanicalAllocationError("closed allocation exceeds disclosed roundoff")
    return result


def _normalized_positive_fractions(weights: tuple[float, ...]) -> tuple[float, ...]:
    scaled = _max_scaled_positive_weights(weights, name="normalization weights")
    try:
        total = math.fsum(scaled)
    except OverflowError as error:  # pragma: no cover - bounded by tuple size in practice
        raise QSCMechanicalAllocationError(
            "max-scaled normalization total is not representable"
        ) from error
    _require_float("max-scaled normalization weight total", total, positive=True)
    fractions = tuple(weight / total for weight in scaled[:-1])
    closing = 1.0 - math.fsum(fractions)
    result = (*fractions, closing)
    _require_float("normalized fractions", *result, positive=True)
    if not _within_disclosed_roundoff(result, 1.0):
        raise QSCMechanicalAllocationError("normalized fractions exceed disclosed roundoff")
    return result


def _packet_wet_mass_under_law(
    packet: AcceptedQSCMechanicalPacket,
    allocator_law: QSCShaftWorkAllocatorLawAuthority,
) -> float:
    totals = _derive_packet_totals_typed(packet.inventory, packet.weight).totals
    try:
        mass = math.fsum(totals[index] for index in allocator_law.weighted_wet_mass_total_indices)
    except OverflowError as error:
        raise QSCMechanicalAllocationError(
            "law-bound packet accepted weighted wet mass is not representable"
        ) from error
    _require_float("law-bound packet accepted weighted wet mass", mass, positive=True)
    if mass != packet.accepted_weighted_wet_mass_kg:
        raise QSCMechanicalAllocationError("packet wet mass differs from allocator-law mapping")
    return mass


def _tray_fractions(
    request: QSCMechanicalAllocationRequest,
    distribution: WithinMealDistribution,
    allocator_law: QSCShaftWorkAllocatorLawAuthority,
) -> tuple[float, ...]:
    packet_masses = tuple(
        _packet_wet_mass_under_law(packet, allocator_law) for packet in request.packets
    )
    scaled_packet_masses = _max_scaled_positive_weights(
        packet_masses,
        name="law-bound packet wet masses",
    )
    accepted_wet_masses = tuple(
        math.fsum(
            scaled_mass
            for packet, scaled_mass in zip(
                request.packets,
                scaled_packet_masses,
                strict=True,
            )
            if packet.physical_tray_id == tray_id
        )
        for tray_id in allocator_law.solid_tray_order
    )
    if distribution is WithinMealDistribution.WET_HOLDUP_PROPORTIONAL_BY_TRAY:
        weights = accepted_wet_masses
    else:
        weights = tuple(
            mass * multiplier
            for mass, multiplier in zip(
                accepted_wet_masses,
                allocator_law.bottom_heavy_rank_multipliers,
                strict=True,
            )
        )
    return _normalized_positive_fractions(weights)


def build_qsc_shaft_work_allocation(
    request: QSCMechanicalAllocationRequest,
    *,
    runtime_authority: QSCMechanicalRuntimeAuthority,
    allocator_law_authority: QSCShaftWorkAllocatorLawAuthority,
    preverified_runtime_definition_digest: str | None = None,
) -> QSCMechanicalWorkAllocation:
    """Build one detached allocation from the exact MECH-01A runtime graph."""

    runtime_digest = validate_qsc_mechanical_runtime_authority(
        runtime_authority,
        preverified_definition_digest=preverified_runtime_definition_digest,
    )
    allocator_law_digest = validate_qsc_shaft_work_allocator_law_authority(allocator_law_authority)
    if type(request) is not QSCMechanicalAllocationRequest:
        raise TypeError("request must be an exact QSCMechanicalAllocationRequest")
    QSCMechanicalAllocationRequest.__post_init__(request)
    if request.allocator_law_definition_digest != allocator_law_digest:
        raise QSCMechanicalAllocationError("request binds a foreign allocator law")
    if runtime_authority.definition_digest != (
        allocator_law_authority.source_runtime_definition_digest
    ):
        raise QSCMechanicalAllocationError("allocator law binds a foreign runtime authority")
    if runtime_authority.shaft_power_bracket_w != (allocator_law_authority.d2_power_endpoints_w):
        raise QSCMechanicalAllocationError("runtime authority contains foreign D2 endpoints")
    # B4 (owner-ruled GO 2026-08-31): the MECH-01A bracket admits the D2
    # continuum; the pinned endpoint tuple itself stays law data (digest
    # unchanged) and supplies the interval bounds.
    if not (
        runtime_authority.shaft_power_bracket_w[0]
        <= request.shaft_power_w
        <= runtime_authority.shaft_power_bracket_w[-1]
    ):
        raise QSCMechanicalAllocationError("request power left the MECH-01A D2 bracket")
    declared_work = request.shaft_power_w * request.dt_s
    _require_float("declared mechanical work", declared_work, positive=True)

    fractions = _tray_fractions(
        request,
        request.tray_distribution,
        allocator_law_authority,
    )
    tray_work = _closed_positive_allocations(declared_work, fractions)
    tray_rows = tuple(
        QSCTrayMechanicalWorkAllocation(
            physical_tray_id=tray_id,
            distribution_fraction=fraction,
            declared_mechanical_work_j=work,
            closing_sixth_remainder=index == 5,
        )
        for index, (tray_id, fraction, work) in enumerate(
            zip(allocator_law_authority.solid_tray_order, fractions, tray_work, strict=True)
        )
    )

    meal_rows: list[MealPacketMechanicalWorkInput] = []
    wall_rows: list[PhysicalTrayWallMechanicalWorkInput] = []
    if request.boundary_partition is ShaftWorkBoundaryPartition.ALL_TO_MEAL:
        for tray_id, work in zip(allocator_law_authority.solid_tray_order, tray_work, strict=True):
            packets = tuple(
                packet for packet in request.packets if packet.physical_tray_id == tray_id
            )
            masses = tuple(
                _packet_wet_mass_under_law(packet, allocator_law_authority) for packet in packets
            )
            packet_fractions = _normalized_positive_fractions(masses)
            packet_work = _closed_positive_allocations(work, packet_fractions)
            meal_rows.extend(
                MealPacketMechanicalWorkInput(
                    physical_tray_id=packet.physical_tray_id,
                    vertical_layer_id=packet.vertical_layer_id,
                    packet_id=packet.packet_id,
                    packet_definition_digest=packet.definition_digest,
                    accepted_weighted_wet_mass_kg=mass,
                    declared_mechanical_work_j=packet_energy,
                    closing_packet_remainder=index == len(packets) - 1,
                )
                for index, (packet, mass, packet_energy) in enumerate(
                    zip(packets, masses, packet_work, strict=True)
                )
            )
    else:
        wall_rows.extend(
            PhysicalTrayWallMechanicalWorkInput(
                physical_tray_id=tray_id,
                declared_mechanical_work_j=work,
            )
            for tray_id, work in zip(
                allocator_law_authority.solid_tray_order, tray_work, strict=True
            )
        )

    meal_tuple = tuple(meal_rows)
    wall_tuple = tuple(wall_rows)
    emitted = tuple(row.declared_mechanical_work_j for row in (*meal_tuple, *wall_tuple))
    meal_total = math.fsum(row.declared_mechanical_work_j for row in meal_tuple)
    wall_total = math.fsum(row.declared_mechanical_work_j for row in wall_tuple)
    allocated_total = math.fsum(emitted)
    exact_declared = Fraction.from_float(declared_work)
    exact_allocated = sum((Fraction.from_float(value) for value in emitted), Fraction())
    exact_residual = exact_allocated - exact_declared
    closure_ulp_limit = (
        allocator_law_authority.roundoff_ulp_coefficient * len(emitted)
        + allocator_law_authority.roundoff_ulp_intercept
    )
    if abs(exact_residual) > (closure_ulp_limit * Fraction.from_float(math.ulp(declared_work))):
        raise QSCMechanicalAllocationError("emitted allocation exceeds roundoff-aware closure")
    request_digest = request.definition_digest
    common = runtime_authority.common_shaft
    ledger = QSCCommonEnergyAllocationLedger(
        energy_datum_id=request.energy_datum_id,
        runtime_definition_digest=runtime_digest,
        allocator_law_definition_digest=allocator_law_digest,
        accepted_step_identity_digest=request.accepted_step_identity_digest,
        request_definition_digest=request_digest,
        boundary_partition=request.boundary_partition,
        shaft_power_w=request.shaft_power_w,
        dt_s=request.dt_s,
        declared_mechanical_work_j=declared_work,
        meal_mechanical_work_j=meal_total,
        wall_mechanical_work_j=wall_total,
        allocated_mechanical_work_j=allocated_total,
        binary64_closure_residual_j=allocated_total - declared_work,
        exact_declared_numerator=_fraction_parts(exact_declared)[0],
        exact_declared_denominator=_fraction_parts(exact_declared)[1],
        exact_allocated_numerator=_fraction_parts(exact_allocated)[0],
        exact_allocated_denominator=_fraction_parts(exact_allocated)[1],
        exact_residual_numerator=_fraction_parts(exact_residual)[0],
        exact_residual_denominator=_fraction_parts(exact_residual)[1],
        closure_ulp_limit=closure_ulp_limit,
        authoritative_conservation=True,
    )
    return QSCMechanicalWorkAllocation(
        schema_id=QSC_SHAFT_WORK_SCHEMA_ID,
        schema_revision=QSC_SHAFT_WORK_SCHEMA_REVISION,
        runtime_definition_digest=runtime_digest,
        allocator_law_definition_digest=allocator_law_digest,
        common_shaft_definition_digest=common.definition_digest,
        accepted_history_definition_digest=common.accepted_shaft_history.definition_digest,
        request_definition_digest=request_digest,
        tray_allocations=tray_rows,
        meal_packet_inputs=meal_tuple,
        wall_inputs=wall_tuple,
        ledger=ledger,
    )


def validate_qsc_shaft_work_allocation(
    allocation: QSCMechanicalWorkAllocation,
    *,
    request: QSCMechanicalAllocationRequest,
    runtime_authority: QSCMechanicalRuntimeAuthority,
    allocator_law_authority: QSCShaftWorkAllocatorLawAuthority,
    expected_request_definition_digest: str,
    expected_allocator_law_definition_digest: str,
    preverified_runtime_definition_digest: str | None = None,
) -> str:
    """Reconstruct against an external request pin and return the allocation digest."""

    if type(allocation) is not QSCMechanicalWorkAllocation:
        raise TypeError("allocation must be an exact QSCMechanicalWorkAllocation")
    QSCMechanicalWorkAllocation.__post_init__(allocation)
    if type(request) is not QSCMechanicalAllocationRequest:
        raise TypeError("request must be an exact QSCMechanicalAllocationRequest")
    QSCMechanicalAllocationRequest.__post_init__(request)
    _require_sha256("expected request definition", expected_request_definition_digest)
    _require_sha256("expected allocator-law definition", expected_allocator_law_definition_digest)
    allocator_law_digest = validate_qsc_shaft_work_allocator_law_authority(allocator_law_authority)
    if allocator_law_digest != expected_allocator_law_definition_digest:
        raise QSCMechanicalAllocationError(
            "allocator law differs from the externally pinned law definition"
        )
    if request.allocator_law_definition_digest != expected_allocator_law_definition_digest:
        raise QSCMechanicalAllocationError("request binds a foreign externally pinned law")
    if request.definition_digest != expected_request_definition_digest:
        raise QSCMechanicalAllocationError(
            "request differs from the externally pinned accepted-step allocation request"
        )
    expected = build_qsc_shaft_work_allocation(
        request,
        runtime_authority=runtime_authority,
        allocator_law_authority=allocator_law_authority,
        preverified_runtime_definition_digest=preverified_runtime_definition_digest,
    )
    if allocation != expected:
        raise QSCMechanicalAllocationError(
            "allocation differs from authoritative reconstruction; tamper or substitution refused"
        )
    return allocation.definition_digest


__all__ = (
    "AcceptedQSCMechanicalPacket",
    "MealPacketMechanicalWorkInput",
    "PhysicalTrayWallMechanicalWorkInput",
    "QSCCommonEnergyAllocationLedger",
    "QSCMechanicalAllocationError",
    "QSCMechanicalAllocationRequest",
    "QSCMechanicalWorkAllocation",
    "QSCTrayMechanicalWorkAllocation",
    "QSCShaftWorkAllocatorLawAuthority",
    "QSC_SHAFT_WORK_ALLOCATOR_LAW_AUTHORITY",
    "QSC_SHAFT_WORK_ALLOCATOR_LAW_SCHEMA_ID",
    "QSC_SHAFT_WORK_ALLOCATOR_LAW_SCHEMA_REVISION",
    "QSC_SHAFT_WORK_EXPECTED_ALLOCATOR_LAW_DEFINITION_DIGEST",
    "QSC_SHAFT_WORK_POWER_ENDPOINTS_W",
    "QSC_SHAFT_WORK_BOTTOM_HEAVY_RANK_MULTIPLIERS",
    "QSC_SHAFT_WORK_SCHEMA_ID",
    "QSC_SHAFT_WORK_SCHEMA_REVISION",
    "build_qsc_shaft_work_allocation",
    "validate_qsc_shaft_work_allocation",
    "validate_qsc_shaft_work_allocator_law_authority",
)
