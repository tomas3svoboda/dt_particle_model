"""Detached PR-08 common-shaft mechanical-authority schema foundation.

The owner-approved PR-08 ruling selects one versioned common-shaft authority
schema before any PR-09 arm or PR-11 shaft-work coefficient campaign.  This
module implements only that immutable, fail-closed data boundary.  It does not
contain a reference-machine instance, infer a value from legacy configuration,
evaluate an arm law, calculate shaft work, authorize solid motion, or wire a
production caller.

Every manufactured numerical field carries its exact unit, finite uncertainty
interval, and source provenance.  One top-level accepted-history reference is
shared by all trays.  Tray records deliberately contain geometry and draw-drive
classification only; they cannot carry a per-tray RPM or shaft history.
"""

from __future__ import annotations

import enum
import functools
import hashlib
import math
import re
import unicodedata
from dataclasses import dataclass
from datetime import datetime
from typing import ClassVar, Iterable

MECHANICAL_AUTHORITY_SCHEMA_ID = "dtdc-core2-pr08-common-shaft-authority-v1"

_SHA256_PATTERN = re.compile(r"sha256:[0-9a-f]{64}\Z")
_UTC_PATTERN = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z\Z")


class MechanicalEvidenceClass(enum.Enum):
    """Declared role of a source record; never an authentication result."""

    MANUFACTURED_TEST = "MANUFACTURED_TEST"
    PROVISIONAL_ENGINEERING = "PROVISIONAL_ENGINEERING"
    SOURCE_RECORD = "SOURCE_RECORD"


class MechanicalUnit(enum.Enum):
    COUNT = "1"
    DIMENSIONLESS = "1.0"
    RPM = "rpm"
    RPM_PER_SECOND = "rpm/s"
    METRE = "m"
    RADIAN = "rad"
    CUBIC_METRE_PER_REVOLUTION = "m3/rev"
    WATT = "W"
    NEWTON_METRE = "N*m"


class MechanicalQuantity(enum.Enum):
    SHAFT_SPEED_ZERO = "shaft_speed_zero"
    SHAFT_SPEED_CRAWL = "shaft_speed_crawl"
    SHAFT_SPEED_MINIMUM_STABLE = "shaft_speed_minimum_stable"
    SHAFT_SPEED_NOMINAL = "shaft_speed_nominal"
    SHAFT_SPEED_MAXIMUM_CONTINUOUS = "shaft_speed_maximum_continuous"
    SHAFT_SPEED_TRIP = "shaft_speed_trip"
    SHAFT_SPEED_RESTART = "shaft_speed_restart"
    SHAFT_RAMP_UP_LIMIT = "shaft_ramp_up_limit"
    SHAFT_RAMP_DOWN_LIMIT = "shaft_ramp_down_limit"
    ARM_COUNT = "arm_count"
    ARM_INNER_RADIUS = "arm_inner_radius"
    ARM_OUTER_RADIUS = "arm_outer_radius"
    BLADE_HEIGHT = "blade_height"
    BLADE_WIDTH = "blade_width"
    BLADE_PITCH = "blade_pitch"
    BLADE_RAKE = "blade_rake"
    TIP_CLEARANCE = "tip_clearance"
    SWEPT_VOLUME_PER_REVOLUTION = "swept_volume_per_revolution"
    GEARBOX_RATIO = "gearbox_ratio"
    RATED_POST_GEARBOX_POWER = "rated_post_gearbox_power"
    RATED_POST_GEARBOX_TORQUE = "rated_post_gearbox_torque"
    DRIVETRAIN_EFFICIENCY = "drivetrain_efficiency"


class SolidsConveyanceDomain(enum.Enum):
    """Mechanical domain that conveys solids through a draw device."""

    COMMON_SHAFT_SWEEP_CONVEYED = "COMMON_SHAFT_SWEEP_CONVEYED"
    INDEPENDENT_DRIVE = "INDEPENDENT_DRIVE"
    PASSIVE_GRAVITY = "PASSIVE_GRAVITY"


class DrawActuationDomain(enum.Enum):
    """Authority domain for opening or actuating a draw device."""

    NO_SEPARATE_ACTUATOR = "NO_SEPARATE_ACTUATOR"
    INDEPENDENT_ACTUATOR = "INDEPENDENT_ACTUATOR"


class DrawDeviceKind(enum.Enum):
    INTERTRAY_SWEEP_PORT = "INTERTRAY_SWEEP_PORT"
    CONTROLLED_GATE = "CONTROLLED_GATE"
    ROTARY_AIRLOCK = "ROTARY_AIRLOCK"


class ShaftHistoryStatus(enum.Enum):
    """Causal status of a referenced history, not a cryptographic attestation."""

    ACCEPTED_MEASURED = "ACCEPTED_MEASURED"
    ACCEPTED_ESTIMATED = "ACCEPTED_ESTIMATED"
    OFFERED_OR_SETPOINT_ONLY = "OFFERED_OR_SETPOINT_ONLY"


class ShaftRotationSense(enum.Enum):
    """Rotation sense relative to an explicitly identified positive shaft axis."""

    POSITIVE_ABOUT_DECLARED_AXIS = "POSITIVE_ABOUT_DECLARED_AXIS"
    NEGATIVE_ABOUT_DECLARED_AXIS = "NEGATIVE_ABOUT_DECLARED_AXIS"
    REVERSIBLE_ABOUT_DECLARED_AXIS = "REVERSIBLE_ABOUT_DECLARED_AXIS"


class ShaftHistorySampleSemantics(enum.Enum):
    """Frozen interpretation of accepted RPM history samples."""

    NONNEGATIVE_RPM_MAGNITUDE_IN_DECLARED_SENSE = "NONNEGATIVE_RPM_MAGNITUDE_IN_DECLARED_SENSE"
    SIGNED_RPM_ABOUT_DECLARED_AXIS = "SIGNED_RPM_ABOUT_DECLARED_AXIS"


class ShaftTripSpeedSemantics(enum.Enum):
    TRIP_THRESHOLD_MAGNITUDE = "TRIP_THRESHOLD_MAGNITUDE"
    POST_TRIP_TARGET_MAGNITUDE = "POST_TRIP_TARGET_MAGNITUDE"


class ShaftRestartSpeedSemantics(enum.Enum):
    RESTART_ENABLE_THRESHOLD_MAGNITUDE = "RESTART_ENABLE_THRESHOLD_MAGNITUDE"
    POST_RESTART_TARGET_MAGNITUDE = "POST_RESTART_TARGET_MAGNITUDE"


class GearboxRatioConvention(enum.Enum):
    MOTOR_SPEED_DIVIDED_BY_POST_GEARBOX_SHAFT_SPEED = (
        "MOTOR_SPEED_DIVIDED_BY_POST_GEARBOX_SHAFT_SPEED"
    )
    POST_GEARBOX_SHAFT_SPEED_DIVIDED_BY_MOTOR_SPEED = (
        "POST_GEARBOX_SHAFT_SPEED_DIVIDED_BY_MOTOR_SPEED"
    )


class DriveTrainEfficiencyBoundary(enum.Enum):
    ELECTRICAL_INPUT_TO_POST_GEARBOX_MECHANICAL_OUTPUT = (
        "ELECTRICAL_INPUT_TO_POST_GEARBOX_MECHANICAL_OUTPUT"
    )
    MOTOR_SHAFT_INPUT_TO_POST_GEARBOX_MECHANICAL_OUTPUT = (
        "MOTOR_SHAFT_INPUT_TO_POST_GEARBOX_MECHANICAL_OUTPUT"
    )


_EXPECTED_UNITS: dict[MechanicalQuantity, MechanicalUnit] = {
    MechanicalQuantity.SHAFT_SPEED_ZERO: MechanicalUnit.RPM,
    MechanicalQuantity.SHAFT_SPEED_CRAWL: MechanicalUnit.RPM,
    MechanicalQuantity.SHAFT_SPEED_MINIMUM_STABLE: MechanicalUnit.RPM,
    MechanicalQuantity.SHAFT_SPEED_NOMINAL: MechanicalUnit.RPM,
    MechanicalQuantity.SHAFT_SPEED_MAXIMUM_CONTINUOUS: MechanicalUnit.RPM,
    MechanicalQuantity.SHAFT_SPEED_TRIP: MechanicalUnit.RPM,
    MechanicalQuantity.SHAFT_SPEED_RESTART: MechanicalUnit.RPM,
    MechanicalQuantity.SHAFT_RAMP_UP_LIMIT: MechanicalUnit.RPM_PER_SECOND,
    MechanicalQuantity.SHAFT_RAMP_DOWN_LIMIT: MechanicalUnit.RPM_PER_SECOND,
    MechanicalQuantity.ARM_COUNT: MechanicalUnit.COUNT,
    MechanicalQuantity.ARM_INNER_RADIUS: MechanicalUnit.METRE,
    MechanicalQuantity.ARM_OUTER_RADIUS: MechanicalUnit.METRE,
    MechanicalQuantity.BLADE_HEIGHT: MechanicalUnit.METRE,
    MechanicalQuantity.BLADE_WIDTH: MechanicalUnit.METRE,
    MechanicalQuantity.BLADE_PITCH: MechanicalUnit.RADIAN,
    MechanicalQuantity.BLADE_RAKE: MechanicalUnit.RADIAN,
    MechanicalQuantity.TIP_CLEARANCE: MechanicalUnit.METRE,
    MechanicalQuantity.SWEPT_VOLUME_PER_REVOLUTION: MechanicalUnit.CUBIC_METRE_PER_REVOLUTION,
    MechanicalQuantity.GEARBOX_RATIO: MechanicalUnit.DIMENSIONLESS,
    MechanicalQuantity.RATED_POST_GEARBOX_POWER: MechanicalUnit.WATT,
    MechanicalQuantity.RATED_POST_GEARBOX_TORQUE: MechanicalUnit.NEWTON_METRE,
    MechanicalQuantity.DRIVETRAIN_EFFICIENCY: MechanicalUnit.DIMENSIONLESS,
}

_NONNEGATIVE_FLOAT_QUANTITIES = frozenset(
    {
        MechanicalQuantity.SHAFT_SPEED_ZERO,
        MechanicalQuantity.SHAFT_SPEED_TRIP,
        MechanicalQuantity.SHAFT_SPEED_RESTART,
        MechanicalQuantity.ARM_INNER_RADIUS,
        MechanicalQuantity.TIP_CLEARANCE,
    }
)

_STRICTLY_POSITIVE_FLOAT_QUANTITIES = frozenset(
    {
        MechanicalQuantity.SHAFT_SPEED_CRAWL,
        MechanicalQuantity.SHAFT_SPEED_MINIMUM_STABLE,
        MechanicalQuantity.SHAFT_SPEED_NOMINAL,
        MechanicalQuantity.SHAFT_SPEED_MAXIMUM_CONTINUOUS,
        MechanicalQuantity.SHAFT_RAMP_UP_LIMIT,
        MechanicalQuantity.SHAFT_RAMP_DOWN_LIMIT,
        MechanicalQuantity.ARM_OUTER_RADIUS,
        MechanicalQuantity.BLADE_HEIGHT,
        MechanicalQuantity.BLADE_WIDTH,
        MechanicalQuantity.SWEPT_VOLUME_PER_REVOLUTION,
        MechanicalQuantity.GEARBOX_RATIO,
        MechanicalQuantity.RATED_POST_GEARBOX_POWER,
        MechanicalQuantity.RATED_POST_GEARBOX_TORQUE,
        MechanicalQuantity.DRIVETRAIN_EFFICIENCY,
    }
)


def _require_nonblank(name: str, *values: str) -> None:
    if any(type(value) is not str for value in values):
        raise TypeError(f"{name} must contain exact strings")
    if any(not value.strip() for value in values):
        raise ValueError(f"{name} must not contain blank values")
    if any(value != value.strip() for value in values):
        raise ValueError(f"{name} must not contain surrounding whitespace")
    if any(
        any(ord(character) < 32 or ord(character) == 127 for character in value) for value in values
    ):
        raise ValueError(f"{name} must not contain control characters")
    for value in values:
        try:
            value.encode("utf-8", errors="strict")
        except UnicodeEncodeError as exc:
            raise ValueError(f"{name} must be valid UTF-8 text") from exc
        if unicodedata.normalize("NFC", value) != value:
            raise ValueError(f"{name} must use NFC-normalized Unicode")


def _require_sha256(name: str, value: str) -> None:
    if type(value) is not str:
        raise TypeError(f"{name} must be an exact string")
    if _SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{name} must be canonical lowercase sha256:<64-hex>")


def _parse_utc(name: str, value: str) -> datetime:
    if type(value) is not str:
        raise TypeError(f"{name} must be an exact string")
    if _UTC_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{name} must use canonical whole-second UTC form")
    try:
        return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ")
    except ValueError as exc:
        raise ValueError(f"{name} is not a valid UTC timestamp") from exc


def _require_binary64(name: str, *values: float) -> None:
    if any(type(value) is not float for value in values):
        raise TypeError(f"{name} must contain exact binary64 values")
    if any(not math.isfinite(value) for value in values):
        raise ValueError(f"{name} must contain only finite values")


def _is_negative_zero(value: float) -> bool:
    return value == 0.0 and math.copysign(1.0, value) < 0.0


def _canonical_digest(
    *,
    domain: str,
    parts: Iterable[str | int | float | bool],
) -> str:
    digest = hashlib.sha256()
    digest.update(b"DTDC-PR08-MECHANICAL-AUTHORITY-V1\0")
    domain_bytes = domain.encode("utf-8")
    digest.update(len(domain_bytes).to_bytes(8, "big"))
    digest.update(domain_bytes)
    material = tuple(parts)
    digest.update(len(material).to_bytes(8, "big"))
    for part in material:
        if type(part) is bool:
            tag = b"b"
            payload = b"1" if part else b"0"
        elif type(part) is int:
            tag = b"i"
            payload = str(part).encode("ascii")
        elif type(part) is float:
            if not math.isfinite(part):
                raise ValueError("mechanical identity floats must be finite")
            tag = b"f"
            payload = part.hex().encode("ascii")
        elif type(part) is str:
            tag = b"s"
            try:
                payload = part.encode("utf-8", errors="strict")
            except UnicodeEncodeError as exc:
                raise ValueError("mechanical identity strings must be valid UTF-8") from exc
            if unicodedata.normalize("NFC", part) != part:
                raise ValueError("mechanical identity strings must use NFC-normalized Unicode")
        else:
            raise TypeError("unsupported mechanical identity part")
        digest.update(tag)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class MechanicalFieldProvenance:
    evidence_class: MechanicalEvidenceClass
    source_id: str
    source_content_sha256: str
    source_locator: str
    acquired_at_utc: str
    uncertainty_basis_id: str

    def __post_init__(self) -> None:
        if type(self.evidence_class) is not MechanicalEvidenceClass:
            raise TypeError("evidence_class must be a MechanicalEvidenceClass")
        _require_nonblank(
            "mechanical provenance strings",
            self.source_id,
            self.source_locator,
            self.uncertainty_basis_id,
        )
        _require_sha256("source_content_sha256", self.source_content_sha256)
        _parse_utc("acquired_at_utc", self.acquired_at_utc)

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="field-provenance",
            parts=(
                self.evidence_class.value,
                self.source_id,
                self.source_content_sha256,
                self.source_locator,
                self.acquired_at_utc,
                self.uncertainty_basis_id,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class MechanicalScalarField:
    quantity: MechanicalQuantity
    unit: MechanicalUnit
    value: float
    lower_bound: float
    upper_bound: float
    provenance: MechanicalFieldProvenance

    def __post_init__(self) -> None:
        if type(self.quantity) is not MechanicalQuantity:
            raise TypeError("quantity must be a MechanicalQuantity")
        if type(self.unit) is not MechanicalUnit:
            raise TypeError("unit must be a MechanicalUnit")
        if self.quantity is MechanicalQuantity.ARM_COUNT:
            raise ValueError("ARM_COUNT requires MechanicalIntegerField")
        if self.unit is not _EXPECTED_UNITS[self.quantity]:
            raise ValueError("mechanical scalar unit does not match its quantity")
        _require_binary64(
            "mechanical scalar and bounds",
            self.value,
            self.lower_bound,
            self.upper_bound,
        )
        if any(
            _is_negative_zero(item) for item in (self.value, self.lower_bound, self.upper_bound)
        ):
            raise ValueError("mechanical scalar values and bounds reject negative zero")
        if not self.lower_bound <= self.value <= self.upper_bound:
            raise ValueError("mechanical scalar value must lie inside its uncertainty interval")
        if self.quantity in _NONNEGATIVE_FLOAT_QUANTITIES and self.lower_bound < 0.0:
            raise ValueError("non-negative mechanical quantity has a negative lower bound")
        if self.quantity in _STRICTLY_POSITIVE_FLOAT_QUANTITIES and self.lower_bound <= 0.0:
            raise ValueError("positive mechanical quantity requires a positive lower bound")
        if self.quantity is MechanicalQuantity.SHAFT_SPEED_ZERO and (
            self.value != 0.0 or self.lower_bound != 0.0 or self.upper_bound != 0.0
        ):
            raise ValueError("the zero shaft-speed authority must be exact zero")
        if self.quantity is MechanicalQuantity.DRIVETRAIN_EFFICIENCY and (self.upper_bound > 1.0):
            raise ValueError("drivetrain efficiency must not exceed one")
        if type(self.provenance) is not MechanicalFieldProvenance:
            raise TypeError("mechanical scalar provenance has the wrong type")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="scalar-field",
            parts=(
                self.quantity.value,
                self.unit.value,
                self.value,
                self.lower_bound,
                self.upper_bound,
                self.provenance.definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class MechanicalIntegerField:
    quantity: MechanicalQuantity
    unit: MechanicalUnit
    value: int
    lower_bound: int
    upper_bound: int
    provenance: MechanicalFieldProvenance

    def __post_init__(self) -> None:
        if self.quantity is not MechanicalQuantity.ARM_COUNT:
            raise ValueError("the only PR-08 integer field is ARM_COUNT")
        if self.unit is not MechanicalUnit.COUNT:
            raise ValueError("arm count must use the exact count unit")
        if any(type(item) is not int for item in (self.value, self.lower_bound, self.upper_bound)):
            raise TypeError("mechanical integer values and bounds must be exact integers")
        if self.lower_bound <= 0 or not self.lower_bound <= self.value <= self.upper_bound:
            raise ValueError("arm count and its uncertainty interval must be positive and ordered")
        if type(self.provenance) is not MechanicalFieldProvenance:
            raise TypeError("mechanical integer provenance has the wrong type")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="integer-field",
            parts=(
                self.quantity.value,
                self.unit.value,
                self.value,
                self.lower_bound,
                self.upper_bound,
                self.provenance.definition_digest,
            ),
        )


def _require_scalar_quantity(
    name: str,
    field: MechanicalScalarField,
    quantity: MechanicalQuantity,
) -> None:
    if type(field) is not MechanicalScalarField:
        raise TypeError(f"{name} must be a MechanicalScalarField")
    if field.quantity is not quantity:
        raise ValueError(f"{name} is bound to the wrong mechanical quantity")


@dataclass(frozen=True, slots=True, kw_only=True)
class MachineConfigurationIdentity:
    unit_id: str
    unit_class_id: str
    serial_id: str
    configuration_id: str
    configuration_revision: str
    provenance: MechanicalFieldProvenance

    def __post_init__(self) -> None:
        _require_nonblank(
            "machine configuration identity",
            self.unit_id,
            self.unit_class_id,
            self.serial_id,
            self.configuration_id,
            self.configuration_revision,
        )
        if type(self.provenance) is not MechanicalFieldProvenance:
            raise TypeError("machine identity provenance has the wrong type")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="machine-configuration",
            parts=(
                self.unit_id,
                self.unit_class_id,
                self.serial_id,
                self.configuration_id,
                self.configuration_revision,
                self.provenance.definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ShaftDirectionAuthority:
    """One machine-bound rotation sense relative to a declared shaft axis."""

    unit_id: str
    configuration_id: str
    configuration_revision: str
    direction_reference_frame_id: str
    rotation_sense: ShaftRotationSense
    provenance: MechanicalFieldProvenance

    def __post_init__(self) -> None:
        _require_nonblank(
            "shaft direction authority",
            self.unit_id,
            self.configuration_id,
            self.configuration_revision,
            self.direction_reference_frame_id,
        )
        if type(self.rotation_sense) is not ShaftRotationSense:
            raise TypeError("rotation_sense must be a ShaftRotationSense")
        if type(self.provenance) is not MechanicalFieldProvenance:
            raise TypeError("shaft direction provenance has the wrong type")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="shaft-direction-authority",
            parts=(
                self.unit_id,
                self.configuration_id,
                self.configuration_revision,
                self.direction_reference_frame_id,
                self.rotation_sense.value,
                self.provenance.definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class AcceptedShaftHistoryReference:
    unit_id: str
    configuration_id: str
    configuration_revision: str
    history_id: str
    history_status: ShaftHistoryStatus
    sample_semantics: ShaftHistorySampleSemantics
    speed_unit: MechanicalUnit
    direction_authority_digest: str
    sensor_channel_id: str
    controller_channel_id: str
    history_schema_id: str
    history_content_length_bytes: int
    history_content_sha256: str
    accepted_interval_start_utc: str
    accepted_interval_end_utc: str
    provenance: MechanicalFieldProvenance

    def __post_init__(self) -> None:
        _require_nonblank(
            "accepted shaft-history identity",
            self.unit_id,
            self.configuration_id,
            self.configuration_revision,
            self.history_id,
            self.sensor_channel_id,
            self.controller_channel_id,
            self.history_schema_id,
        )
        if type(self.history_status) is not ShaftHistoryStatus:
            raise TypeError("history_status must be a ShaftHistoryStatus")
        if self.history_status is ShaftHistoryStatus.OFFERED_OR_SETPOINT_ONLY:
            raise ValueError("offered or setpoint-only RPM cannot be an accepted shaft history")
        if self.speed_unit is not MechanicalUnit.RPM:
            raise ValueError("accepted shaft history must declare exact RPM samples")
        if type(self.sample_semantics) is not ShaftHistorySampleSemantics:
            raise TypeError("sample_semantics must be a ShaftHistorySampleSemantics")
        _require_sha256("direction_authority_digest", self.direction_authority_digest)
        if type(self.history_content_length_bytes) is not int:
            raise TypeError("shaft-history content length must be an exact integer")
        if self.history_content_length_bytes <= 0:
            raise ValueError("shaft-history content must be nonempty")
        _require_sha256("history_content_sha256", self.history_content_sha256)
        start = _parse_utc("accepted_interval_start_utc", self.accepted_interval_start_utc)
        end = _parse_utc("accepted_interval_end_utc", self.accepted_interval_end_utc)
        if not start < end:
            raise ValueError("accepted shaft-history interval must have positive duration")
        if type(self.provenance) is not MechanicalFieldProvenance:
            raise TypeError("shaft-history provenance has the wrong type")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="accepted-shaft-history",
            parts=(
                self.unit_id,
                self.configuration_id,
                self.configuration_revision,
                self.history_id,
                self.history_status.value,
                self.sample_semantics.value,
                self.speed_unit.value,
                self.direction_authority_digest,
                self.sensor_channel_id,
                self.controller_channel_id,
                self.history_schema_id,
                self.history_content_length_bytes,
                self.history_content_sha256,
                self.accepted_interval_start_utc,
                self.accepted_interval_end_utc,
                self.provenance.definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ShaftSpeedEnvelope:
    trip_speed_semantics: ShaftTripSpeedSemantics
    restart_speed_semantics: ShaftRestartSpeedSemantics
    semantics_provenance: MechanicalFieldProvenance
    zero_rpm: MechanicalScalarField
    crawl_rpm: MechanicalScalarField
    minimum_stable_rpm: MechanicalScalarField
    nominal_rpm: MechanicalScalarField
    maximum_continuous_rpm: MechanicalScalarField
    trip_rpm: MechanicalScalarField
    restart_rpm: MechanicalScalarField
    ramp_up_limit_rpm_per_s: MechanicalScalarField
    ramp_down_limit_rpm_per_s: MechanicalScalarField

    def __post_init__(self) -> None:
        if type(self.trip_speed_semantics) is not ShaftTripSpeedSemantics:
            raise TypeError("trip_speed_semantics must be a ShaftTripSpeedSemantics")
        if type(self.restart_speed_semantics) is not ShaftRestartSpeedSemantics:
            raise TypeError("restart_speed_semantics must be a ShaftRestartSpeedSemantics")
        if type(self.semantics_provenance) is not MechanicalFieldProvenance:
            raise TypeError("shaft-speed semantics provenance has the wrong type")
        expected = (
            ("zero_rpm", self.zero_rpm, MechanicalQuantity.SHAFT_SPEED_ZERO),
            ("crawl_rpm", self.crawl_rpm, MechanicalQuantity.SHAFT_SPEED_CRAWL),
            (
                "minimum_stable_rpm",
                self.minimum_stable_rpm,
                MechanicalQuantity.SHAFT_SPEED_MINIMUM_STABLE,
            ),
            ("nominal_rpm", self.nominal_rpm, MechanicalQuantity.SHAFT_SPEED_NOMINAL),
            (
                "maximum_continuous_rpm",
                self.maximum_continuous_rpm,
                MechanicalQuantity.SHAFT_SPEED_MAXIMUM_CONTINUOUS,
            ),
            ("trip_rpm", self.trip_rpm, MechanicalQuantity.SHAFT_SPEED_TRIP),
            ("restart_rpm", self.restart_rpm, MechanicalQuantity.SHAFT_SPEED_RESTART),
            (
                "ramp_up_limit_rpm_per_s",
                self.ramp_up_limit_rpm_per_s,
                MechanicalQuantity.SHAFT_RAMP_UP_LIMIT,
            ),
            (
                "ramp_down_limit_rpm_per_s",
                self.ramp_down_limit_rpm_per_s,
                MechanicalQuantity.SHAFT_RAMP_DOWN_LIMIT,
            ),
        )
        for name, field, quantity in expected:
            _require_scalar_quantity(name, field, quantity)
        ordered = (
            self.crawl_rpm.value,
            self.minimum_stable_rpm.value,
            self.nominal_rpm.value,
            self.maximum_continuous_rpm.value,
        )
        if not ordered[0] <= ordered[1] <= ordered[2] <= ordered[3]:
            raise ValueError("crawl/minimum-stable/nominal/maximum-continuous RPM must be ordered")
        ordered_fields = (
            self.crawl_rpm,
            self.minimum_stable_rpm,
            self.nominal_rpm,
            self.maximum_continuous_rpm,
        )
        if any(
            lower.upper_bound >= upper.lower_bound
            for lower, upper in zip(ordered_fields[:-1], ordered_fields[1:], strict=True)
        ):
            raise ValueError("ordinary operating RPM uncertainty intervals must be disjoint")
        if self.trip_speed_semantics is ShaftTripSpeedSemantics.TRIP_THRESHOLD_MAGNITUDE and (
            self.trip_rpm.lower_bound <= self.maximum_continuous_rpm.upper_bound
        ):
            raise ValueError("trip threshold must be above maximum-continuous RPM")
        if (
            self.trip_speed_semantics is ShaftTripSpeedSemantics.POST_TRIP_TARGET_MAGNITUDE
            and self.trip_rpm.upper_bound > self.maximum_continuous_rpm.upper_bound
        ):
            raise ValueError("post-trip target must lie below maximum-continuous RPM")
        if (
            self.restart_speed_semantics
            is ShaftRestartSpeedSemantics.RESTART_ENABLE_THRESHOLD_MAGNITUDE
            and self.restart_rpm.upper_bound >= self.minimum_stable_rpm.lower_bound
        ):
            raise ValueError("restart-enable threshold must be below minimum-stable RPM")
        if (
            self.restart_speed_semantics is ShaftRestartSpeedSemantics.POST_RESTART_TARGET_MAGNITUDE
            and (
                self.restart_rpm.lower_bound < self.crawl_rpm.lower_bound
                or self.restart_rpm.upper_bound > self.maximum_continuous_rpm.upper_bound
            )
        ):
            raise ValueError("post-restart target must lie inside the operating envelope")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="shaft-speed-envelope",
            parts=(
                self.trip_speed_semantics.value,
                self.restart_speed_semantics.value,
                self.semantics_provenance.definition_digest,
                *(
                    field.definition_digest
                    for field in (
                        self.zero_rpm,
                        self.crawl_rpm,
                        self.minimum_stable_rpm,
                        self.nominal_rpm,
                        self.maximum_continuous_rpm,
                        self.trip_rpm,
                        self.restart_rpm,
                        self.ramp_up_limit_rpm_per_s,
                        self.ramp_down_limit_rpm_per_s,
                    )
                ),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ArmGeometryAuthority:
    physical_tray_id: str
    arm_geometry_id: str
    coordinate_frame_id: str
    identity_provenance: MechanicalFieldProvenance
    arm_count: MechanicalIntegerField
    inner_radius_m: MechanicalScalarField
    outer_radius_m: MechanicalScalarField
    blade_height_m: MechanicalScalarField
    blade_width_m: MechanicalScalarField
    blade_pitch_rad: MechanicalScalarField
    blade_rake_rad: MechanicalScalarField
    tip_clearance_m: MechanicalScalarField
    swept_volume_per_revolution_m3: MechanicalScalarField

    def __post_init__(self) -> None:
        _require_nonblank(
            "arm geometry identity",
            self.physical_tray_id,
            self.arm_geometry_id,
            self.coordinate_frame_id,
        )
        if type(self.identity_provenance) is not MechanicalFieldProvenance:
            raise TypeError("arm geometry identity provenance has the wrong type")
        if type(self.arm_count) is not MechanicalIntegerField:
            raise TypeError("arm_count must be a MechanicalIntegerField")
        expected = (
            ("inner_radius_m", self.inner_radius_m, MechanicalQuantity.ARM_INNER_RADIUS),
            ("outer_radius_m", self.outer_radius_m, MechanicalQuantity.ARM_OUTER_RADIUS),
            ("blade_height_m", self.blade_height_m, MechanicalQuantity.BLADE_HEIGHT),
            ("blade_width_m", self.blade_width_m, MechanicalQuantity.BLADE_WIDTH),
            ("blade_pitch_rad", self.blade_pitch_rad, MechanicalQuantity.BLADE_PITCH),
            ("blade_rake_rad", self.blade_rake_rad, MechanicalQuantity.BLADE_RAKE),
            ("tip_clearance_m", self.tip_clearance_m, MechanicalQuantity.TIP_CLEARANCE),
            (
                "swept_volume_per_revolution_m3",
                self.swept_volume_per_revolution_m3,
                MechanicalQuantity.SWEPT_VOLUME_PER_REVOLUTION,
            ),
        )
        for name, field, quantity in expected:
            _require_scalar_quantity(name, field, quantity)
        if self.inner_radius_m.upper_bound >= self.outer_radius_m.lower_bound:
            raise ValueError("arm inner and outer radius uncertainty intervals must be disjoint")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="arm-geometry",
            parts=(
                self.physical_tray_id,
                self.arm_geometry_id,
                self.coordinate_frame_id,
                self.identity_provenance.definition_digest,
                self.arm_count.definition_digest,
                self.inner_radius_m.definition_digest,
                self.outer_radius_m.definition_digest,
                self.blade_height_m.definition_digest,
                self.blade_width_m.definition_digest,
                self.blade_pitch_rad.definition_digest,
                self.blade_rake_rad.definition_digest,
                self.tip_clearance_m.definition_digest,
                self.swept_volume_per_revolution_m3.definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class DrawGeometryReference:
    """Identity-only draw hardware; it selects no footprint or flow law.

    Solids conveyance and device actuation are deliberately distinct.  A gate
    may meter a common-shaft or passive-gravity stream while its actuator is an
    independent device.  A bottom rotary airlock instead conveys solids with
    its own drive.  Neither classification authorizes a discharge rate.
    """

    physical_tray_id: str
    draw_device_id: str
    device_kind: DrawDeviceKind
    conveyance_domain: SolidsConveyanceDomain
    actuation_domain: DrawActuationDomain
    draw_geometry_schema_id: str
    draw_geometry_content_length_bytes: int
    draw_geometry_content_sha256: str
    draw_geometry_uncertainty_content_sha256: str
    independent_conveyance_authority_digest: str | None
    independent_actuation_authority_digest: str | None
    geometry_provenance: MechanicalFieldProvenance
    conveyance_domain_provenance: MechanicalFieldProvenance
    actuation_domain_provenance: MechanicalFieldProvenance

    def __post_init__(self) -> None:
        _require_nonblank(
            "draw geometry identity",
            self.physical_tray_id,
            self.draw_device_id,
            self.draw_geometry_schema_id,
        )
        if type(self.device_kind) is not DrawDeviceKind:
            raise TypeError("device_kind must be a DrawDeviceKind")
        if type(self.conveyance_domain) is not SolidsConveyanceDomain:
            raise TypeError("conveyance_domain must be a SolidsConveyanceDomain")
        if type(self.actuation_domain) is not DrawActuationDomain:
            raise TypeError("actuation_domain must be a DrawActuationDomain")
        if type(self.draw_geometry_content_length_bytes) is not int:
            raise TypeError("draw geometry length must be an exact integer")
        if self.draw_geometry_content_length_bytes <= 0:
            raise ValueError("draw geometry content must be nonempty")
        _require_sha256("draw_geometry_content_sha256", self.draw_geometry_content_sha256)
        _require_sha256(
            "draw_geometry_uncertainty_content_sha256",
            self.draw_geometry_uncertainty_content_sha256,
        )
        if self.device_kind is DrawDeviceKind.INTERTRAY_SWEEP_PORT and (
            self.conveyance_domain is not SolidsConveyanceDomain.COMMON_SHAFT_SWEEP_CONVEYED
            or self.actuation_domain is not DrawActuationDomain.NO_SEPARATE_ACTUATOR
        ):
            raise ValueError(
                "an inter-tray sweep port requires common-shaft conveyance and no actuator"
            )
        if self.device_kind is DrawDeviceKind.CONTROLLED_GATE and (
            self.conveyance_domain is SolidsConveyanceDomain.INDEPENDENT_DRIVE
            or self.actuation_domain is not DrawActuationDomain.INDEPENDENT_ACTUATOR
        ):
            raise ValueError(
                "a controlled gate requires an independent actuator and cannot propel solids"
            )
        if self.device_kind is DrawDeviceKind.ROTARY_AIRLOCK and (
            self.conveyance_domain is not SolidsConveyanceDomain.INDEPENDENT_DRIVE
            or self.actuation_domain is not DrawActuationDomain.INDEPENDENT_ACTUATOR
        ):
            raise ValueError("a rotary airlock requires independent conveyance and actuation")
        if self.conveyance_domain is SolidsConveyanceDomain.INDEPENDENT_DRIVE:
            if type(self.independent_conveyance_authority_digest) is not str:
                raise TypeError("independent solids conveyance requires an authority digest")
            _require_sha256(
                "independent_conveyance_authority_digest",
                self.independent_conveyance_authority_digest,
            )
        elif self.independent_conveyance_authority_digest is not None:
            raise ValueError("only independent solids conveyance may bind a conveyance authority")
        if self.actuation_domain is DrawActuationDomain.INDEPENDENT_ACTUATOR:
            if type(self.independent_actuation_authority_digest) is not str:
                raise TypeError("an independent actuator requires an authority digest")
            _require_sha256(
                "independent_actuation_authority_digest",
                self.independent_actuation_authority_digest,
            )
        elif self.independent_actuation_authority_digest is not None:
            raise ValueError("only an independent actuator may bind an actuator authority")
        provenance_values = (
            self.geometry_provenance,
            self.conveyance_domain_provenance,
            self.actuation_domain_provenance,
        )
        if any(type(item) is not MechanicalFieldProvenance for item in provenance_values):
            raise TypeError("draw geometry, conveyance, and actuation provenance must be explicit")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="draw-geometry-reference",
            parts=(
                self.physical_tray_id,
                self.draw_device_id,
                self.device_kind.value,
                self.conveyance_domain.value,
                self.actuation_domain.value,
                self.draw_geometry_schema_id,
                self.draw_geometry_content_length_bytes,
                self.draw_geometry_content_sha256,
                self.draw_geometry_uncertainty_content_sha256,
                self.independent_conveyance_authority_digest or "NO-INDEPENDENT-CONVEYANCE",
                self.independent_actuation_authority_digest or "NO-INDEPENDENT-ACTUATOR",
                self.geometry_provenance.definition_digest,
                self.conveyance_domain_provenance.definition_digest,
                self.actuation_domain_provenance.definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class TrayMechanicalAuthority:
    physical_tray_id: str
    arm_geometry: ArmGeometryAuthority
    draw_geometry: DrawGeometryReference

    def __post_init__(self) -> None:
        _require_nonblank("tray mechanical identity", self.physical_tray_id)
        if type(self.arm_geometry) is not ArmGeometryAuthority:
            raise TypeError("arm_geometry must be an ArmGeometryAuthority")
        if type(self.draw_geometry) is not DrawGeometryReference:
            raise TypeError("draw_geometry must be a DrawGeometryReference")
        if self.arm_geometry.physical_tray_id != self.physical_tray_id:
            raise ValueError("arm geometry belongs to another physical tray")
        if self.draw_geometry.physical_tray_id != self.physical_tray_id:
            raise ValueError("draw geometry belongs to another physical tray")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="tray-mechanical-authority",
            parts=(
                self.physical_tray_id,
                self.arm_geometry.definition_digest,
                self.draw_geometry.definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class DriveTrainAuthority:
    unit_id: str
    configuration_id: str
    configuration_revision: str
    motor_id: str
    gearbox_id: str
    drive_configuration_id: str
    gearbox_ratio_convention: GearboxRatioConvention
    efficiency_boundary: DriveTrainEfficiencyBoundary
    identity_provenance: MechanicalFieldProvenance
    convention_provenance: MechanicalFieldProvenance
    gearbox_ratio: MechanicalScalarField
    rated_post_gearbox_power_w: MechanicalScalarField
    rated_post_gearbox_torque_n_m: MechanicalScalarField
    efficiency: MechanicalScalarField

    def __post_init__(self) -> None:
        _require_nonblank(
            "drive-train identity",
            self.unit_id,
            self.configuration_id,
            self.configuration_revision,
            self.motor_id,
            self.gearbox_id,
            self.drive_configuration_id,
        )
        if type(self.identity_provenance) is not MechanicalFieldProvenance:
            raise TypeError("drive-train identity provenance has the wrong type")
        if type(self.gearbox_ratio_convention) is not GearboxRatioConvention:
            raise TypeError("gearbox_ratio_convention has the wrong type")
        if type(self.efficiency_boundary) is not DriveTrainEfficiencyBoundary:
            raise TypeError("efficiency_boundary has the wrong type")
        if type(self.convention_provenance) is not MechanicalFieldProvenance:
            raise TypeError("drive-train convention provenance has the wrong type")
        expected = (
            ("gearbox_ratio", self.gearbox_ratio, MechanicalQuantity.GEARBOX_RATIO),
            (
                "rated_post_gearbox_power_w",
                self.rated_post_gearbox_power_w,
                MechanicalQuantity.RATED_POST_GEARBOX_POWER,
            ),
            (
                "rated_post_gearbox_torque_n_m",
                self.rated_post_gearbox_torque_n_m,
                MechanicalQuantity.RATED_POST_GEARBOX_TORQUE,
            ),
            ("efficiency", self.efficiency, MechanicalQuantity.DRIVETRAIN_EFFICIENCY),
        )
        for name, field, quantity in expected:
            _require_scalar_quantity(name, field, quantity)

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="drive-train-authority",
            parts=(
                self.unit_id,
                self.configuration_id,
                self.configuration_revision,
                self.motor_id,
                self.gearbox_id,
                self.drive_configuration_id,
                self.gearbox_ratio_convention.value,
                self.efficiency_boundary.value,
                self.identity_provenance.definition_digest,
                self.convention_provenance.definition_digest,
                self.gearbox_ratio.definition_digest,
                self.rated_post_gearbox_power_w.definition_digest,
                self.rated_post_gearbox_torque_n_m.definition_digest,
                self.efficiency.definition_digest,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class CommonShaftMechanicalAuthority:
    schema_id: str
    authority_revision: int
    machine: MachineConfigurationIdentity
    shaft_direction: ShaftDirectionAuthority
    accepted_shaft_history: AcceptedShaftHistoryReference
    speed_envelope: ShaftSpeedEnvelope
    drive_train: DriveTrainAuthority
    trays: tuple[TrayMechanicalAuthority, ...]

    architecture_only: ClassVar[bool] = True
    owner_approved_pr08_schema_requirement_applied: ClassVar[bool] = True
    reference_machine_data_bound: ClassVar[bool] = False
    accepted_history_authenticated: ClassVar[bool] = False
    source_authenticated: ClassVar[bool] = False
    complete_transitive_authentication: ClassVar[bool] = False
    security_attestation: ClassVar[bool] = False
    machine_coefficients_identified: ClassVar[bool] = False
    production_mechanical_authority_implemented: ClassVar[bool] = False
    production_host_wired: ClassVar[bool] = False
    production_draw_geometry_selected: ClassVar[bool] = False
    production_arm_coefficient_selected: ClassVar[bool] = False
    production_shaft_work_magnitude_selected: ClassVar[bool] = False
    arm_law_implemented: ClassVar[bool] = False
    shaft_work_model_implemented: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if self.schema_id != MECHANICAL_AUTHORITY_SCHEMA_ID:
            raise ValueError("mechanical authority schema_id is not the frozen V1 schema")
        if type(self.authority_revision) is not int:
            raise TypeError("authority_revision must be an exact integer")
        if self.authority_revision <= 0:
            raise ValueError("authority_revision must be positive")
        if type(self.machine) is not MachineConfigurationIdentity:
            raise TypeError("machine must be a MachineConfigurationIdentity")
        if type(self.shaft_direction) is not ShaftDirectionAuthority:
            raise TypeError("shaft_direction must be a ShaftDirectionAuthority")
        direction_subject = (
            self.shaft_direction.unit_id,
            self.shaft_direction.configuration_id,
            self.shaft_direction.configuration_revision,
        )
        machine_subject = (
            self.machine.unit_id,
            self.machine.configuration_id,
            self.machine.configuration_revision,
        )
        if direction_subject != machine_subject:
            raise ValueError("shaft direction belongs to another machine configuration")
        if type(self.accepted_shaft_history) is not AcceptedShaftHistoryReference:
            raise TypeError("accepted_shaft_history has the wrong type")
        if self.accepted_shaft_history.unit_id != self.machine.unit_id:
            raise ValueError("accepted shaft history belongs to another machine unit")
        if self.accepted_shaft_history.configuration_id != self.machine.configuration_id:
            raise ValueError("accepted shaft history belongs to another machine configuration")
        if (
            self.accepted_shaft_history.configuration_revision
            != self.machine.configuration_revision
        ):
            raise ValueError("accepted shaft history uses another configuration revision")
        if (
            self.accepted_shaft_history.direction_authority_digest
            != self.shaft_direction.definition_digest
        ):
            raise ValueError("accepted shaft history binds another shaft direction authority")
        reversible = (
            self.shaft_direction.rotation_sense is ShaftRotationSense.REVERSIBLE_ABOUT_DECLARED_AXIS
        )
        expected_sample_semantics = (
            ShaftHistorySampleSemantics.SIGNED_RPM_ABOUT_DECLARED_AXIS
            if reversible
            else ShaftHistorySampleSemantics.NONNEGATIVE_RPM_MAGNITUDE_IN_DECLARED_SENSE
        )
        if self.accepted_shaft_history.sample_semantics is not expected_sample_semantics:
            raise ValueError("shaft-history sample semantics disagree with rotation authority")
        if type(self.speed_envelope) is not ShaftSpeedEnvelope:
            raise TypeError("speed_envelope has the wrong type")
        if type(self.drive_train) is not DriveTrainAuthority:
            raise TypeError("drive_train has the wrong type")
        drive_subject = (
            self.drive_train.unit_id,
            self.drive_train.configuration_id,
            self.drive_train.configuration_revision,
        )
        if drive_subject != machine_subject:
            raise ValueError("drive train belongs to another machine configuration")
        if type(self.trays) is not tuple:
            raise TypeError("trays must be an exact immutable tuple")
        if not self.trays:
            raise ValueError("at least one physical DT tray is required")
        if any(type(tray) is not TrayMechanicalAuthority for tray in self.trays):
            raise TypeError("trays must contain TrayMechanicalAuthority values")
        tray_ids = tuple(tray.physical_tray_id for tray in self.trays)
        if len(set(tray_ids)) != len(tray_ids):
            raise ValueError("mechanical authority contains duplicate physical trays")
        if tray_ids != tuple(sorted(tray_ids)):
            raise ValueError("mechanical authority trays must use canonical physical-tray order")
        common_shaft_authorities = {
            self.machine.definition_digest,
            self.shaft_direction.definition_digest,
            self.accepted_shaft_history.definition_digest,
            self.speed_envelope.definition_digest,
            self.drive_train.definition_digest,
        }
        for tray in self.trays:
            draw = tray.draw_geometry
            independent_digests = (
                draw.independent_conveyance_authority_digest,
                draw.independent_actuation_authority_digest,
            )
            if any(digest in common_shaft_authorities for digest in independent_digests):
                raise ValueError(
                    "an independent draw authority must not alias the common shaft authority"
                )

    @property
    def definition_digest(self) -> str:
        # Owner-approved B9 (2026-08-31): value-keyed memo through the
        # frozen dataclass hash — mutation re-walks; the constant common
        # shaft's repeated recursive walks collapse to one.
        return _common_shaft_definition_digest(self)


@functools.lru_cache(maxsize=16)
def _common_shaft_definition_digest(authority: CommonShaftMechanicalAuthority) -> str:
    return _canonical_digest(
        domain="common-shaft-mechanical-authority",
        parts=(
            authority.schema_id,
            authority.authority_revision,
            authority.machine.definition_digest,
            authority.shaft_direction.definition_digest,
            authority.accepted_shaft_history.definition_digest,
            authority.speed_envelope.definition_digest,
            authority.drive_train.definition_digest,
            len(authority.trays),
            *(tray.definition_digest for tray in authority.trays),
        ),
    )


def _rebuild_provenance(value: MechanicalFieldProvenance) -> MechanicalFieldProvenance:
    if type(value) is not MechanicalFieldProvenance:
        raise TypeError("nested provenance has the wrong exact type")
    if type(value.evidence_class) is not MechanicalEvidenceClass:
        raise TypeError("nested evidence_class has the wrong exact type")
    return MechanicalFieldProvenance(
        evidence_class=value.evidence_class,
        source_id=value.source_id,
        source_content_sha256=value.source_content_sha256,
        source_locator=value.source_locator,
        acquired_at_utc=value.acquired_at_utc,
        uncertainty_basis_id=value.uncertainty_basis_id,
    )


def _rebuild_scalar(value: MechanicalScalarField) -> MechanicalScalarField:
    if type(value) is not MechanicalScalarField:
        raise TypeError("nested scalar field has the wrong exact type")
    if type(value.quantity) is not MechanicalQuantity or type(value.unit) is not MechanicalUnit:
        raise TypeError("nested scalar quantity or unit has the wrong exact type")
    return MechanicalScalarField(
        quantity=value.quantity,
        unit=value.unit,
        value=value.value,
        lower_bound=value.lower_bound,
        upper_bound=value.upper_bound,
        provenance=_rebuild_provenance(value.provenance),
    )


def _rebuild_integer(value: MechanicalIntegerField) -> MechanicalIntegerField:
    if type(value) is not MechanicalIntegerField:
        raise TypeError("nested integer field has the wrong exact type")
    if type(value.quantity) is not MechanicalQuantity or type(value.unit) is not MechanicalUnit:
        raise TypeError("nested integer quantity or unit has the wrong exact type")
    return MechanicalIntegerField(
        quantity=value.quantity,
        unit=value.unit,
        value=value.value,
        lower_bound=value.lower_bound,
        upper_bound=value.upper_bound,
        provenance=_rebuild_provenance(value.provenance),
    )


def _rebuild_machine(value: MachineConfigurationIdentity) -> MachineConfigurationIdentity:
    if type(value) is not MachineConfigurationIdentity:
        raise TypeError("nested machine identity has the wrong exact type")
    return MachineConfigurationIdentity(
        unit_id=value.unit_id,
        unit_class_id=value.unit_class_id,
        serial_id=value.serial_id,
        configuration_id=value.configuration_id,
        configuration_revision=value.configuration_revision,
        provenance=_rebuild_provenance(value.provenance),
    )


def _rebuild_direction(value: ShaftDirectionAuthority) -> ShaftDirectionAuthority:
    if type(value) is not ShaftDirectionAuthority:
        raise TypeError("nested shaft direction has the wrong exact type")
    if type(value.rotation_sense) is not ShaftRotationSense:
        raise TypeError("nested rotation sense has the wrong exact type")
    return ShaftDirectionAuthority(
        unit_id=value.unit_id,
        configuration_id=value.configuration_id,
        configuration_revision=value.configuration_revision,
        direction_reference_frame_id=value.direction_reference_frame_id,
        rotation_sense=value.rotation_sense,
        provenance=_rebuild_provenance(value.provenance),
    )


def _rebuild_history(value: AcceptedShaftHistoryReference) -> AcceptedShaftHistoryReference:
    if type(value) is not AcceptedShaftHistoryReference:
        raise TypeError("nested shaft-history reference has the wrong exact type")
    if type(value.history_status) is not ShaftHistoryStatus:
        raise TypeError("nested history status has the wrong exact type")
    if type(value.speed_unit) is not MechanicalUnit:
        raise TypeError("nested history speed unit has the wrong exact type")
    if type(value.sample_semantics) is not ShaftHistorySampleSemantics:
        raise TypeError("nested history sample semantics has the wrong exact type")
    return AcceptedShaftHistoryReference(
        unit_id=value.unit_id,
        configuration_id=value.configuration_id,
        configuration_revision=value.configuration_revision,
        history_id=value.history_id,
        history_status=value.history_status,
        sample_semantics=value.sample_semantics,
        speed_unit=value.speed_unit,
        direction_authority_digest=value.direction_authority_digest,
        sensor_channel_id=value.sensor_channel_id,
        controller_channel_id=value.controller_channel_id,
        history_schema_id=value.history_schema_id,
        history_content_length_bytes=value.history_content_length_bytes,
        history_content_sha256=value.history_content_sha256,
        accepted_interval_start_utc=value.accepted_interval_start_utc,
        accepted_interval_end_utc=value.accepted_interval_end_utc,
        provenance=_rebuild_provenance(value.provenance),
    )


def _rebuild_speed_envelope(value: ShaftSpeedEnvelope) -> ShaftSpeedEnvelope:
    if type(value) is not ShaftSpeedEnvelope:
        raise TypeError("nested shaft-speed envelope has the wrong exact type")
    if type(value.trip_speed_semantics) is not ShaftTripSpeedSemantics:
        raise TypeError("nested trip semantics has the wrong exact type")
    if type(value.restart_speed_semantics) is not ShaftRestartSpeedSemantics:
        raise TypeError("nested restart semantics has the wrong exact type")
    return ShaftSpeedEnvelope(
        trip_speed_semantics=value.trip_speed_semantics,
        restart_speed_semantics=value.restart_speed_semantics,
        semantics_provenance=_rebuild_provenance(value.semantics_provenance),
        zero_rpm=_rebuild_scalar(value.zero_rpm),
        crawl_rpm=_rebuild_scalar(value.crawl_rpm),
        minimum_stable_rpm=_rebuild_scalar(value.minimum_stable_rpm),
        nominal_rpm=_rebuild_scalar(value.nominal_rpm),
        maximum_continuous_rpm=_rebuild_scalar(value.maximum_continuous_rpm),
        trip_rpm=_rebuild_scalar(value.trip_rpm),
        restart_rpm=_rebuild_scalar(value.restart_rpm),
        ramp_up_limit_rpm_per_s=_rebuild_scalar(value.ramp_up_limit_rpm_per_s),
        ramp_down_limit_rpm_per_s=_rebuild_scalar(value.ramp_down_limit_rpm_per_s),
    )


def _rebuild_arm(value: ArmGeometryAuthority) -> ArmGeometryAuthority:
    if type(value) is not ArmGeometryAuthority:
        raise TypeError("nested arm geometry has the wrong exact type")
    return ArmGeometryAuthority(
        physical_tray_id=value.physical_tray_id,
        arm_geometry_id=value.arm_geometry_id,
        coordinate_frame_id=value.coordinate_frame_id,
        identity_provenance=_rebuild_provenance(value.identity_provenance),
        arm_count=_rebuild_integer(value.arm_count),
        inner_radius_m=_rebuild_scalar(value.inner_radius_m),
        outer_radius_m=_rebuild_scalar(value.outer_radius_m),
        blade_height_m=_rebuild_scalar(value.blade_height_m),
        blade_width_m=_rebuild_scalar(value.blade_width_m),
        blade_pitch_rad=_rebuild_scalar(value.blade_pitch_rad),
        blade_rake_rad=_rebuild_scalar(value.blade_rake_rad),
        tip_clearance_m=_rebuild_scalar(value.tip_clearance_m),
        swept_volume_per_revolution_m3=_rebuild_scalar(value.swept_volume_per_revolution_m3),
    )


def _rebuild_draw(value: DrawGeometryReference) -> DrawGeometryReference:
    if type(value) is not DrawGeometryReference:
        raise TypeError("nested draw geometry has the wrong exact type")
    if type(value.device_kind) is not DrawDeviceKind:
        raise TypeError("nested draw device kind has the wrong exact type")
    if type(value.conveyance_domain) is not SolidsConveyanceDomain:
        raise TypeError("nested conveyance domain has the wrong exact type")
    if type(value.actuation_domain) is not DrawActuationDomain:
        raise TypeError("nested actuation domain has the wrong exact type")
    return DrawGeometryReference(
        physical_tray_id=value.physical_tray_id,
        draw_device_id=value.draw_device_id,
        device_kind=value.device_kind,
        conveyance_domain=value.conveyance_domain,
        actuation_domain=value.actuation_domain,
        draw_geometry_schema_id=value.draw_geometry_schema_id,
        draw_geometry_content_length_bytes=value.draw_geometry_content_length_bytes,
        draw_geometry_content_sha256=value.draw_geometry_content_sha256,
        draw_geometry_uncertainty_content_sha256=(value.draw_geometry_uncertainty_content_sha256),
        independent_conveyance_authority_digest=(value.independent_conveyance_authority_digest),
        independent_actuation_authority_digest=value.independent_actuation_authority_digest,
        geometry_provenance=_rebuild_provenance(value.geometry_provenance),
        conveyance_domain_provenance=_rebuild_provenance(value.conveyance_domain_provenance),
        actuation_domain_provenance=_rebuild_provenance(value.actuation_domain_provenance),
    )


def _rebuild_tray(value: TrayMechanicalAuthority) -> TrayMechanicalAuthority:
    if type(value) is not TrayMechanicalAuthority:
        raise TypeError("nested tray authority has the wrong exact type")
    return TrayMechanicalAuthority(
        physical_tray_id=value.physical_tray_id,
        arm_geometry=_rebuild_arm(value.arm_geometry),
        draw_geometry=_rebuild_draw(value.draw_geometry),
    )


def _rebuild_drive_train(value: DriveTrainAuthority) -> DriveTrainAuthority:
    if type(value) is not DriveTrainAuthority:
        raise TypeError("nested drive-train authority has the wrong exact type")
    if type(value.gearbox_ratio_convention) is not GearboxRatioConvention:
        raise TypeError("nested gearbox-ratio convention has the wrong exact type")
    if type(value.efficiency_boundary) is not DriveTrainEfficiencyBoundary:
        raise TypeError("nested efficiency boundary has the wrong exact type")
    return DriveTrainAuthority(
        unit_id=value.unit_id,
        configuration_id=value.configuration_id,
        configuration_revision=value.configuration_revision,
        motor_id=value.motor_id,
        gearbox_id=value.gearbox_id,
        drive_configuration_id=value.drive_configuration_id,
        gearbox_ratio_convention=value.gearbox_ratio_convention,
        efficiency_boundary=value.efficiency_boundary,
        identity_provenance=_rebuild_provenance(value.identity_provenance),
        convention_provenance=_rebuild_provenance(value.convention_provenance),
        gearbox_ratio=_rebuild_scalar(value.gearbox_ratio),
        rated_post_gearbox_power_w=_rebuild_scalar(value.rated_post_gearbox_power_w),
        rated_post_gearbox_torque_n_m=_rebuild_scalar(value.rated_post_gearbox_torque_n_m),
        efficiency=_rebuild_scalar(value.efficiency),
    )


def _rebuild_authority(
    authority: CommonShaftMechanicalAuthority,
) -> CommonShaftMechanicalAuthority:
    if type(authority) is not CommonShaftMechanicalAuthority:
        raise TypeError("authority must be a CommonShaftMechanicalAuthority")
    if type(authority.trays) is not tuple:
        raise TypeError("trays must be an exact immutable tuple")
    return CommonShaftMechanicalAuthority(
        schema_id=authority.schema_id,
        authority_revision=authority.authority_revision,
        machine=_rebuild_machine(authority.machine),
        shaft_direction=_rebuild_direction(authority.shaft_direction),
        accepted_shaft_history=_rebuild_history(authority.accepted_shaft_history),
        speed_envelope=_rebuild_speed_envelope(authority.speed_envelope),
        drive_train=_rebuild_drive_train(authority.drive_train),
        trays=tuple(_rebuild_tray(tray) for tray in authority.trays),
    )


def validate_common_shaft_mechanical_authority(
    authority: CommonShaftMechanicalAuthority,
) -> str:
    """Return the canonical identity after recursive exact-type validation."""

    return _rebuild_authority(authority).definition_digest


__all__ = [
    "AcceptedShaftHistoryReference",
    "ArmGeometryAuthority",
    "CommonShaftMechanicalAuthority",
    "DriveTrainEfficiencyBoundary",
    "DrawActuationDomain",
    "DrawDeviceKind",
    "DrawGeometryReference",
    "DriveTrainAuthority",
    "GearboxRatioConvention",
    "MECHANICAL_AUTHORITY_SCHEMA_ID",
    "MachineConfigurationIdentity",
    "MechanicalEvidenceClass",
    "MechanicalFieldProvenance",
    "MechanicalIntegerField",
    "MechanicalQuantity",
    "MechanicalScalarField",
    "MechanicalUnit",
    "ShaftSpeedEnvelope",
    "ShaftDirectionAuthority",
    "ShaftHistoryStatus",
    "ShaftHistorySampleSemantics",
    "ShaftRestartSpeedSemantics",
    "ShaftRotationSense",
    "ShaftTripSpeedSemantics",
    "SolidsConveyanceDomain",
    "TrayMechanicalAuthority",
    "validate_common_shaft_mechanical_authority",
]
