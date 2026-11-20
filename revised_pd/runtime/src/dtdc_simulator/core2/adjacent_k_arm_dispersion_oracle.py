"""Detached manufactured PR-09/PR-10 arm-response oracle.

This module implements only the owner-selected mathematical structure needed
to qualify an adjacent-K finite-volume arm-dispersion law.  All numerical
inputs are explicit manufactured fixtures bound by generic SHA-256 identities.
No reference-machine record, coefficient identification, transfer acceptance,
RTD response, production configuration, or production caller is supplied.

The finite-volume convention is

    D = 0.5 * phi_ex * delta_z_ex**2 * (accepted_rpm / 60) * f_state
    g_f = D * rho_dm_bulk,f * A_f / distance_f
    q_left_to_right = g_f / M_dm,left
    q_right_to_left = g_f / M_dm,right

where ``g_f`` is the shared dry-matter conductance.  Applying the paired jump
rates approximates the same shared-face flux in binary64; the direct face-flux
evaluator exposes the finite-volume expression independently.  No face
averaging or gradient mapping rule is hidden: face area, dry-matter bulk
density, and gradient distance are separately pinned geometry inputs.

The PR-10 helper is a veto classifier only.  It never scales a reported flow
and never authorizes material motion.  This module deliberately accepts no RTD
order, moment, residence, or exit-selection input.
"""

from __future__ import annotations

import enum
import hashlib
import math
import re
import sys
import unicodedata
from dataclasses import dataclass
from fractions import Fraction
from typing import ClassVar, Iterable, Sequence

from .joint_carrier_oracle import ArmLayerTransition, ArmRedistributionAuthority
from .tray_particle_stateful import VerticalLayerId


ADJACENT_K_ARM_DISPERSION_SCHEMA_ID = (
    "dtdc-core2-pr09-pr10-adjacent-k-arm-dispersion-manufactured-v1"
)
MAX_MANUFACTURED_LAYER_COUNT = 4096
MAX_MANUFACTURED_SEGMENT_COUNT = 4096
MAX_DENSE_GENERATOR_LAYER_COUNT = 256
MAX_IDENTITY_UTF8_BYTES = 4096

_SHA256_PATTERN = re.compile(r"sha256:[0-9a-f]{64}\Z")


def _require_identity(name: str, *values: str) -> None:
    for value in values:
        if type(value) is not str:
            raise TypeError(f"{name} must contain exact strings")
        if not value or value != value.strip():
            raise ValueError(f"{name} must not be blank or have surrounding whitespace")
        if unicodedata.normalize("NFC", value) != value:
            raise ValueError(f"{name} must use NFC-normalized Unicode")
        if any(unicodedata.category(character) == "Cc" for character in value):
            raise ValueError(f"{name} must not contain control characters")
        try:
            encoded = value.encode("utf-8", errors="strict")
        except UnicodeEncodeError as error:
            raise ValueError(f"{name} must be valid strict UTF-8") from error
        if len(encoded) > MAX_IDENTITY_UTF8_BYTES:
            raise ValueError(f"{name} exceeds the manufactured identity resource limit")


def _require_sha256(name: str, value: str) -> None:
    if type(value) is not str:
        raise TypeError(f"{name} must be an exact string")
    if _SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{name} must be canonical lowercase sha256:<64-hex>")


def _is_negative_zero(value: float) -> bool:
    return value == 0.0 and math.copysign(1.0, value) < 0.0


def _require_binary64(name: str, *values: float, reject_negative_zero: bool = True) -> None:
    if any(type(value) is not float for value in values):
        raise TypeError(f"{name} must contain exact binary64 floats")
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"{name} must contain only finite values")
    if reject_negative_zero and any(_is_negative_zero(value) for value in values):
        raise ValueError(f"{name} rejects negative zero")


def _require_positive_normal(name: str, value: float) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and strictly positive")
    if value < sys.float_info.min:
        raise ValueError(f"{name} underflowed the normal binary64 domain")
    return value


def _checked_fsum(name: str, values: Iterable[float]) -> float:
    try:
        result = math.fsum(values)
    except OverflowError as error:
        raise ValueError(f"{name} is not representable in binary64") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} is not representable in binary64")
    return result


def _checked_positive_product_ratio(
    name: str,
    *,
    numerators: tuple[float, ...],
    denominators: tuple[float, ...] = (),
) -> float:
    """Correctly round a positive product/ratio and require a normal result."""

    if type(numerators) is not tuple or type(denominators) is not tuple:
        raise TypeError(f"{name} factors must be exact tuples")
    if not numerators:
        raise ValueError(f"{name} requires at least one numerator")
    _require_binary64(name, *numerators, *denominators)
    if any(value <= 0.0 for value in (*numerators, *denominators)):
        raise ValueError(f"{name} factors must be strictly positive")
    exact = Fraction(1)
    for value in numerators:
        exact *= Fraction.from_float(value)
    for value in denominators:
        exact /= Fraction.from_float(value)
    try:
        result = float(exact)
    except OverflowError as error:
        raise ValueError(f"{name} is not representable in binary64") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} is not representable in binary64")
    return _require_positive_normal(name, result)


def _canonical_digest(*, domain: str, parts: Iterable[str | int | float | bool]) -> str:
    _require_identity("digest domain", domain)
    digest = hashlib.sha256()
    digest.update(b"DTDC-ADJACENT-K-ARM-ORACLE-V1\0")
    domain_bytes = domain.encode("utf-8")
    digest.update(len(domain_bytes).to_bytes(8, "big"))
    digest.update(domain_bytes)
    material = tuple(parts)
    digest.update(len(material).to_bytes(8, "big"))
    for part in material:
        if isinstance(part, bool):
            tag, payload = b"b", b"1" if part else b"0"
        elif type(part) is int:
            tag, payload = b"i", str(part).encode("ascii")
        elif type(part) is float:
            _require_binary64("canonical digest floats", part)
            tag, payload = b"f", part.hex().encode("ascii")
        elif type(part) is str:
            _require_identity("canonical digest strings", part)
            tag, payload = b"s", part.encode("utf-8")
        else:  # pragma: no cover - closed union guarded by callers
            raise TypeError(f"unsupported canonical digest part: {type(part).__name__}")
        digest.update(tag)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return "sha256:" + digest.hexdigest()


class ArmDispersionMode(enum.Enum):
    ACTIVE_BASELINE = "ACTIVE_BASELINE"
    INACTIVE_BED = "INACTIVE_BED"
    STRATIFIED_LIMIT = "STRATIFIED_LIMIT"


class ShaftBoundaryKind(enum.Enum):
    """Exact semantics of a shared piecewise-history boundary."""

    CONTINUOUS = "CONTINUOUS"
    INSTANTANEOUS_STEP = "INSTANTANEOUS_STEP"


class PR10ConveyanceDomain(enum.Enum):
    COMMON_SHAFT_SWEEP_CONVEYED = "COMMON_SHAFT_SWEEP_CONVEYED"
    INDEPENDENT_DRIVE = "INDEPENDENT_DRIVE"
    PASSIVE_GRAVITY = "PASSIVE_GRAVITY"


class PR10ConveyanceDisposition(enum.Enum):
    VETO_ZERO_RPM_COMMON_SHAFT_SWEEP = "VETO_ZERO_RPM_COMMON_SHAFT_SWEEP"
    NO_PR10_VETO_POSITIVE_COMMON_SHAFT_SWEEP = "NO_PR10_VETO_POSITIVE_COMMON_SHAFT_SWEEP"
    OUTSIDE_COMMON_SHAFT_SWEEP_DOMAIN = "OUTSIDE_COMMON_SHAFT_SWEEP_DOMAIN"


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedAcceptedShaftSample:
    common_shaft_definition_digest: str
    accepted_history_definition_digest: str
    sample_id: str
    sample_time_s: float
    shaft_speed_rpm: float

    source_authenticated: ClassVar[bool] = False
    accepted_history_authenticated: ClassVar[bool] = False
    production_authority: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_sha256("common_shaft_definition_digest", self.common_shaft_definition_digest)
        _require_sha256(
            "accepted_history_definition_digest", self.accepted_history_definition_digest
        )
        _require_identity("accepted shaft sample identity", self.sample_id)
        _require_binary64("accepted shaft sample", self.sample_time_s, self.shaft_speed_rpm)
        if self.sample_time_s < 0.0:
            raise ValueError("sample_time_s must be non-negative")
        if self.shaft_speed_rpm < 0.0:
            raise ValueError("shaft_speed_rpm must be non-negative; reversal is not selected")
        if self.sample_time_s > 0.0:
            _require_positive_normal("positive accepted sample time", self.sample_time_s)
        if self.shaft_speed_rpm > 0.0:
            _require_positive_normal("positive accepted shaft speed", self.shaft_speed_rpm)

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="manufactured-accepted-shaft-sample",
            parts=(
                self.common_shaft_definition_digest,
                self.accepted_history_definition_digest,
                self.sample_id,
                self.sample_time_s,
                self.shaft_speed_rpm,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedAcceptedShaftSegment:
    common_shaft_definition_digest: str
    accepted_history_definition_digest: str
    segment_id: str
    start_event_id: str
    end_event_id: str
    start_event_kind: ShaftBoundaryKind
    end_event_kind: ShaftBoundaryKind
    time_start_s: float
    time_end_s: float
    rpm_start: float
    rpm_end: float

    source_authenticated: ClassVar[bool] = False
    accepted_history_authenticated: ClassVar[bool] = False
    production_authority: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_sha256("common_shaft_definition_digest", self.common_shaft_definition_digest)
        _require_sha256(
            "accepted_history_definition_digest", self.accepted_history_definition_digest
        )
        _require_identity(
            "accepted shaft segment identity",
            self.segment_id,
            self.start_event_id,
            self.end_event_id,
        )
        if type(self.start_event_kind) is not ShaftBoundaryKind:
            raise TypeError("start_event_kind must be an exact ShaftBoundaryKind")
        if type(self.end_event_kind) is not ShaftBoundaryKind:
            raise TypeError("end_event_kind must be an exact ShaftBoundaryKind")
        _require_binary64(
            "accepted shaft segment",
            self.time_start_s,
            self.time_end_s,
            self.rpm_start,
            self.rpm_end,
        )
        if self.time_start_s < 0.0 or self.time_end_s <= self.time_start_s:
            raise ValueError("shaft segment times must be non-negative and strictly ordered")
        if self.rpm_start < 0.0 or self.rpm_end < 0.0:
            raise ValueError("shaft segment RPM must be non-negative; reversal is not selected")
        _require_positive_normal("positive shaft segment duration", self.duration_s)
        if self.rpm_start > 0.0:
            _require_positive_normal("positive shaft segment start RPM", self.rpm_start)
        if self.rpm_end > 0.0:
            _require_positive_normal("positive shaft segment end RPM", self.rpm_end)
        if self.start_event_id == self.end_event_id:
            raise ValueError("a positive-duration shaft segment needs distinct boundary events")

    @property
    def duration_s(self) -> float:
        return self.time_end_s - self.time_start_s

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="manufactured-accepted-shaft-segment",
            parts=(
                self.common_shaft_definition_digest,
                self.accepted_history_definition_digest,
                self.segment_id,
                self.start_event_id,
                self.end_event_id,
                self.start_event_kind.value,
                self.end_event_kind.value,
                self.time_start_s,
                self.time_end_s,
                self.rpm_start,
                self.rpm_end,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedArmClosure:
    closure_id: str
    closure_source_digest: str
    tray_geometry_definition_digest: str
    phi_ex: float
    delta_z_ex_m: float
    f_state: float

    source_authenticated: ClassVar[bool] = False
    machine_coefficients_identified: ClassVar[bool] = False
    state_modifier_identified: ClassVar[bool] = False
    directed_conveyance_selected: ClassVar[bool] = False
    saturation_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("manufactured arm closure identity", self.closure_id)
        _require_sha256("closure_source_digest", self.closure_source_digest)
        _require_sha256("tray_geometry_definition_digest", self.tray_geometry_definition_digest)
        _require_binary64("manufactured arm closure", self.phi_ex, self.delta_z_ex_m, self.f_state)
        if not 0.0 < self.phi_ex <= 1.0:
            raise ValueError("phi_ex must lie in the manufactured per-pass interval (0, 1]")
        if self.delta_z_ex_m <= 0.0:
            raise ValueError("delta_z_ex_m must be strictly positive")
        if self.f_state <= 0.0:
            raise ValueError("f_state must be explicit and strictly positive")
        _require_positive_normal("phi_ex", self.phi_ex)
        _require_positive_normal("delta_z_ex_m", self.delta_z_ex_m)
        _require_positive_normal("f_state", self.f_state)

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="manufactured-arm-closure",
            parts=(
                self.closure_id,
                self.closure_source_digest,
                self.tray_geometry_definition_digest,
                self.phi_ex,
                self.delta_z_ex_m,
                self.f_state,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PhysicalHeightCell:
    vertical_layer_id: VerticalLayerId
    z_lower_m: float
    z_upper_m: float
    dry_matter_holdup_kg: float

    def __post_init__(self) -> None:
        if type(self.vertical_layer_id) is not VerticalLayerId:
            raise TypeError("vertical_layer_id must be an exact VerticalLayerId")
        if type(self.vertical_layer_id.value) is not int or self.vertical_layer_id.value <= 0:
            raise ValueError("vertical_layer_id must contain a positive exact integer")
        _require_binary64(
            "physical-height cell",
            self.z_lower_m,
            self.z_upper_m,
            self.dry_matter_holdup_kg,
        )
        if self.z_lower_m < 0.0 or self.z_upper_m <= self.z_lower_m:
            raise ValueError("physical-height cell bounds must be non-negative and ordered")
        if self.dry_matter_holdup_kg <= 0.0:
            raise ValueError("cell dry-matter holdup must be strictly positive")
        _require_positive_normal("cell thickness", self.z_upper_m - self.z_lower_m)
        _require_positive_normal("cell dry-matter holdup", self.dry_matter_holdup_kg)
        _ = self.centroid_m

    @property
    def centroid_m(self) -> float:
        centroid = _checked_fsum(
            "physical-height cell centroid",
            (0.5 * self.z_lower_m, 0.5 * self.z_upper_m),
        )
        if not math.isfinite(centroid):
            raise ValueError("physical-height cell centroid is not representable")
        return centroid

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="physical-height-cell",
            parts=(
                self.vertical_layer_id.value,
                self.z_lower_m,
                self.z_upper_m,
                self.dry_matter_holdup_kg,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PhysicalHeightFace:
    lower_vertical_layer_id: VerticalLayerId
    upper_vertical_layer_id: VerticalLayerId
    face_area_m2: float
    dry_matter_bulk_density_kg_m3: float
    gradient_distance_m: float

    def __post_init__(self) -> None:
        endpoints = (self.lower_vertical_layer_id, self.upper_vertical_layer_id)
        if any(type(value) is not VerticalLayerId for value in endpoints):
            raise TypeError("physical-height face endpoints must be exact VerticalLayerId values")
        if any(type(value.value) is not int or value.value <= 0 for value in endpoints):
            raise ValueError("physical-height face endpoints must be positive exact integers")
        if self.upper_vertical_layer_id.value != self.lower_vertical_layer_id.value + 1:
            raise ValueError("arm-dispersion faces must connect exactly adjacent K layers")
        _require_binary64(
            "physical-height face",
            self.face_area_m2,
            self.dry_matter_bulk_density_kg_m3,
            self.gradient_distance_m,
        )
        if (
            self.face_area_m2 <= 0.0
            or self.dry_matter_bulk_density_kg_m3 <= 0.0
            or self.gradient_distance_m <= 0.0
        ):
            raise ValueError(
                "face area, dry-matter bulk density, and gradient distance must be positive"
            )
        _require_positive_normal("face area", self.face_area_m2)
        _require_positive_normal("face dry-matter bulk density", self.dry_matter_bulk_density_kg_m3)
        _require_positive_normal("face gradient distance", self.gradient_distance_m)

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="physical-height-face",
            parts=(
                self.lower_vertical_layer_id.value,
                self.upper_vertical_layer_id.value,
                self.face_area_m2,
                self.dry_matter_bulk_density_kg_m3,
                self.gradient_distance_m,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PhysicalHeightMesh:
    physical_tray_id: str
    mesh_id: str
    tray_geometry_definition_digest: str
    cells: tuple[PhysicalHeightCell, ...]
    faces: tuple[PhysicalHeightFace, ...]

    numerical_resolution_only: ClassVar[bool] = True
    production_geometry_selected: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("physical-height mesh identity", self.physical_tray_id, self.mesh_id)
        _require_sha256("tray_geometry_definition_digest", self.tray_geometry_definition_digest)
        if type(self.cells) is not tuple or type(self.faces) is not tuple:
            raise TypeError("physical-height cells and faces must be exact immutable tuples")
        if not 2 <= len(self.cells) <= MAX_MANUFACTURED_LAYER_COUNT:
            raise ValueError("manufactured physical-height K is outside the resource domain")
        if len(self.faces) != len(self.cells) - 1:
            raise ValueError("a closed adjacent-K mesh requires exactly K-1 internal faces")
        if any(type(cell) is not PhysicalHeightCell for cell in self.cells):
            raise TypeError("mesh cells must be exact PhysicalHeightCell values")
        if any(type(face) is not PhysicalHeightFace for face in self.faces):
            raise TypeError("mesh faces must be exact PhysicalHeightFace values")
        rebuilt_cells = tuple(_rebuild_cell(cell) for cell in self.cells)
        rebuilt_faces = tuple(_rebuild_face(face) for face in self.faces)
        if self.cells != rebuilt_cells or self.faces != rebuilt_faces:
            raise ValueError("mesh children differ from their reconstructed definitions")
        expected_ids = tuple(range(1, len(self.cells) + 1))
        if tuple(cell.vertical_layer_id.value for cell in self.cells) != expected_ids:
            raise ValueError(
                "mesh cells must use canonical contiguous K identities starting at one"
            )
        for lower, upper in zip(self.cells[:-1], self.cells[1:], strict=True):
            if lower.z_upper_m != upper.z_lower_m:
                raise ValueError("physical-height cells must be exactly contiguous")
        expected_faces = tuple((index, index + 1) for index in range(1, len(self.cells)))
        actual_faces = tuple(
            (face.lower_vertical_layer_id.value, face.upper_vertical_layer_id.value)
            for face in self.faces
        )
        if actual_faces != expected_faces:
            raise ValueError("mesh faces must use canonical complete adjacent-K order")

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="physical-height-mesh",
            parts=(
                self.physical_tray_id,
                self.mesh_id,
                self.tray_geometry_definition_digest,
                len(self.cells),
                *(cell.definition_digest for cell in self.cells),
                len(self.faces),
                *(face.definition_digest for face in self.faces),
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedArmFaceRate:
    lower_vertical_layer_id: VerticalLayerId
    upper_vertical_layer_id: VerticalLayerId
    gradient_distance_m: float
    diffusivity_m2_s: float
    dry_matter_conductance_kg_s: float
    lower_to_upper_rate_s_inv: float
    upper_to_lower_rate_s_inv: float

    def __post_init__(self) -> None:
        endpoints = (self.lower_vertical_layer_id, self.upper_vertical_layer_id)
        if any(type(value) is not VerticalLayerId for value in endpoints):
            raise TypeError("face-rate endpoints must be exact VerticalLayerId values")
        if any(type(value.value) is not int or value.value <= 0 for value in endpoints):
            raise ValueError("face-rate endpoints must contain positive exact integers")
        if self.upper_vertical_layer_id.value != self.lower_vertical_layer_id.value + 1:
            raise ValueError("face rates must connect exactly adjacent K layers")
        values = (
            self.gradient_distance_m,
            self.diffusivity_m2_s,
            self.dry_matter_conductance_kg_s,
            self.lower_to_upper_rate_s_inv,
            self.upper_to_lower_rate_s_inv,
        )
        _require_binary64("manufactured arm face rate", *values)
        if self.gradient_distance_m <= 0.0:
            raise ValueError("face-rate gradient distance must be strictly positive")
        _require_positive_normal("face-rate gradient distance", self.gradient_distance_m)
        rate_values = values[1:]
        if any(value < 0.0 for value in rate_values):
            raise ValueError("manufactured arm face rates must be non-negative")
        zero_pattern = tuple(value == 0.0 for value in rate_values)
        if any(zero_pattern) and not all(zero_pattern):
            raise ValueError("an arm face must be either exactly zero or wholly positive")
        if all(value > 0.0 for value in rate_values):
            for value in rate_values:
                _require_positive_normal("positive manufactured arm face rate", value)

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="manufactured-arm-face-rate",
            parts=(
                self.lower_vertical_layer_id.value,
                self.upper_vertical_layer_id.value,
                self.gradient_distance_m,
                self.diffusivity_m2_s,
                self.dry_matter_conductance_kg_s,
                self.lower_to_upper_rate_s_inv,
                self.upper_to_lower_rate_s_inv,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedAdjacentKArmOperator:
    schema_id: str
    mesh: PhysicalHeightMesh
    closure: ManufacturedArmClosure
    accepted_sample: ManufacturedAcceptedShaftSample
    mode: ArmDispersionMode
    face_rates: tuple[ManufacturedArmFaceRate, ...]
    transitions: tuple[ArmLayerTransition, ...]

    architecture_only: ClassVar[bool] = True
    source_authenticated: ClassVar[bool] = False
    accepted_history_authenticated: ClassVar[bool] = False
    machine_data_authenticated: ClassVar[bool] = False
    machine_coefficients_identified: ClassVar[bool] = False
    production_implemented: ClassVar[bool] = False
    production_host_wired: ClassVar[bool] = False
    transfer_authority: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_identity("adjacent-K arm operator schema identity", self.schema_id)
        if self.schema_id != ADJACENT_K_ARM_DISPERSION_SCHEMA_ID:
            raise ValueError("adjacent-K arm operator has the wrong schema identity")
        if type(self.mesh) is not PhysicalHeightMesh:
            raise TypeError("mesh must be an exact PhysicalHeightMesh")
        if type(self.closure) is not ManufacturedArmClosure:
            raise TypeError("closure must be an exact ManufacturedArmClosure")
        if type(self.accepted_sample) is not ManufacturedAcceptedShaftSample:
            raise TypeError("accepted_sample must be an exact ManufacturedAcceptedShaftSample")
        if type(self.mode) is not ArmDispersionMode:
            raise TypeError("mode must be an exact ArmDispersionMode")
        if type(self.face_rates) is not tuple or type(self.transitions) is not tuple:
            raise TypeError("face rates and transitions must be exact immutable tuples")
        if len(self.face_rates) != len(self.mesh.faces):
            raise ValueError("operator requires one face-rate record per mesh face")
        if any(type(value) is not ManufacturedArmFaceRate for value in self.face_rates):
            raise TypeError("operator face_rates contain a wrong exact type")
        if any(type(value) is not ArmLayerTransition for value in self.transitions):
            raise TypeError("operator transitions contain a wrong exact type")
        rebuilt_face_rates = tuple(_rebuild_face_rate(value) for value in self.face_rates)
        rebuilt_transitions = tuple(_rebuild_transition(value) for value in self.transitions)
        if self.face_rates != rebuilt_face_rates or self.transitions != rebuilt_transitions:
            raise ValueError("operator records differ from their reconstructed definitions")
        if self.closure.tray_geometry_definition_digest != (
            self.mesh.tray_geometry_definition_digest
        ):
            raise ValueError("closure and mesh bind different tray geometry identities")
        if _rebuild_mesh(self.mesh) != self.mesh:
            raise ValueError("operator mesh differs from its reconstructed definition")
        if _rebuild_closure(self.closure) != self.closure:
            raise ValueError("operator closure differs from its reconstructed definition")
        if _rebuild_sample(self.accepted_sample) != self.accepted_sample:
            raise ValueError("operator shaft sample differs from its reconstructed definition")
        should_be_zero = (
            self.accepted_sample.shaft_speed_rpm == 0.0
            or self.mode is not ArmDispersionMode.ACTIVE_BASELINE
        )
        if should_be_zero:
            if self.transitions or any(rate.diffusivity_m2_s != 0.0 for rate in self.face_rates):
                raise ValueError("zero/inactive/stratified arm operator must be exactly zero")
        else:
            if len(self.transitions) != 2 * len(self.mesh.faces):
                raise ValueError("active arm operator requires a paired transition per face")
            if any(rate.diffusivity_m2_s <= 0.0 for rate in self.face_rates):
                raise ValueError("active arm operator must be connected-positive")
        canonical_transitions = tuple(
            sorted(
                self.transitions,
                key=lambda value: (
                    value.source_vertical_layer_id.value,
                    value.destination_vertical_layer_id.value,
                ),
            )
        )
        if self.transitions != canonical_transitions:
            raise ValueError("operator transitions must use canonical directed-edge order")
        expected_face_rates, expected_transitions = _derive_face_records(
            mesh=self.mesh,
            closure=self.closure,
            accepted_sample=self.accepted_sample,
            mode=self.mode,
        )
        if self.face_rates != expected_face_rates or self.transitions != expected_transitions:
            raise ValueError("operator rates differ from the frozen finite-volume formula")

    @property
    def diffusivity_m2_s(self) -> float:
        return self.face_rates[0].diffusivity_m2_s

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="manufactured-adjacent-k-arm-operator",
            parts=(
                self.schema_id,
                self.mesh.definition_digest,
                self.closure.definition_digest,
                self.accepted_sample.definition_digest,
                self.mode.value,
                len(self.face_rates),
                *(rate.definition_digest for rate in self.face_rates),
                len(self.transitions),
                *(
                    part
                    for transition in self.transitions
                    for part in (
                        transition.source_vertical_layer_id.value,
                        transition.destination_vertical_layer_id.value,
                        transition.rate_s_inv,
                    )
                ),
            ),
        )

    def generator_matrix(self) -> tuple[tuple[float, ...], ...]:
        """Return the conservative arm generator in column/source convention."""

        count = len(self.mesh.cells)
        if count > MAX_DENSE_GENERATOR_LAYER_COUNT:
            raise ValueError("dense arm generator exceeds the manufactured resource limit")
        contributions: list[list[list[float]]] = [[[] for _ in range(count)] for _ in range(count)]
        outgoing: list[list[float]] = [[] for _ in range(count)]
        for transition in self.transitions:
            source = transition.source_vertical_layer_id.value - 1
            destination = transition.destination_vertical_layer_id.value - 1
            contributions[destination][source].append(transition.rate_s_inv)
            outgoing[source].append(transition.rate_s_inv)
        for source, values in enumerate(outgoing):
            contributions[source][source].append(-_checked_fsum("outgoing arm rate", values))
        return tuple(
            tuple(_checked_fsum("arm generator entry", entry) for entry in row)
            for row in contributions
        )

    def evaluate_extensive_rates(self, inventories: tuple[float, ...]) -> tuple[float, ...]:
        """Apply the same paired incidence/rate operator to one extensive."""

        if type(inventories) is not tuple:
            raise TypeError("extensive inventories must be an exact tuple")
        if len(inventories) != len(self.mesh.cells):
            raise ValueError("one extensive inventory is required per K cell")
        _require_binary64("extensive inventories", *inventories, reject_negative_zero=False)
        contributions: list[list[float]] = [[] for _ in inventories]
        for transition in self.transitions:
            source = transition.source_vertical_layer_id.value - 1
            destination = transition.destination_vertical_layer_id.value - 1
            amount = transition.rate_s_inv * inventories[source]
            if not math.isfinite(amount):
                raise ValueError("extensive arm rate is not representable")
            if inventories[source] != 0.0 and (amount == 0.0 or abs(amount) < sys.float_info.min):
                raise ValueError("nonzero extensive arm rate underflowed binary64")
            contributions[source].append(-amount)
            contributions[destination].append(amount)
        return tuple(_checked_fsum("extensive arm node rate", values) for values in contributions)

    def evaluate_direct_face_fluxes(self, inventories: tuple[float, ...]) -> tuple[float, ...]:
        """Evaluate the declared FV concentration flux without the jump adapter.

        The returned values use the lower-to-upper sign convention.  This path
        exposes the independent face formula and lets qualification quantify
        binary64 differences introduced by ``q=g/M`` followed by ``q*I``.
        """

        if type(inventories) is not tuple:
            raise TypeError("extensive inventories must be an exact tuple")
        if len(inventories) != len(self.mesh.cells):
            raise ValueError("one extensive inventory is required per K cell")
        _require_binary64("extensive inventories", *inventories, reject_negative_zero=False)
        cells = {cell.vertical_layer_id.value: cell for cell in self.mesh.cells}
        fluxes: list[float] = []
        for face_rate in self.face_rates:
            if face_rate.dry_matter_conductance_kg_s == 0.0:
                fluxes.append(0.0)
                continue
            lower_index = face_rate.lower_vertical_layer_id.value - 1
            upper_index = face_rate.upper_vertical_layer_id.value - 1
            lower_mass = cells[face_rate.lower_vertical_layer_id.value].dry_matter_holdup_kg
            upper_mass = cells[face_rate.upper_vertical_layer_id.value].dry_matter_holdup_kg
            lower_concentration = inventories[lower_index] / lower_mass
            upper_concentration = inventories[upper_index] / upper_mass
            if not math.isfinite(lower_concentration) or not math.isfinite(upper_concentration):
                raise ValueError("direct FV concentration is not representable")
            for inventory, concentration in (
                (inventories[lower_index], lower_concentration),
                (inventories[upper_index], upper_concentration),
            ):
                if inventory != 0.0 and (
                    concentration == 0.0 or abs(concentration) < sys.float_info.min
                ):
                    raise ValueError("nonzero direct FV concentration underflowed binary64")
            difference = _checked_fsum(
                "direct FV concentration difference",
                (lower_concentration, -upper_concentration),
            )
            flux = face_rate.dry_matter_conductance_kg_s * difference
            if not math.isfinite(flux):
                raise ValueError("direct FV face flux is not representable")
            if (
                difference != 0.0
                and face_rate.dry_matter_conductance_kg_s != 0.0
                and (flux == 0.0 or abs(flux) < sys.float_info.min)
            ):
                raise ValueError("nonzero direct FV face flux underflowed binary64")
            fluxes.append(flux)
        return tuple(fluxes)

    def transition_stationarity_residual(self) -> tuple[float, ...]:
        """Expose jump-adapter roundoff for the dry-matter equilibrium vector."""

        return self.evaluate_extensive_rates(
            tuple(cell.dry_matter_holdup_kg for cell in self.mesh.cells)
        )

    def transition_conservation_residual(self, inventories: tuple[float, ...]) -> float:
        """Expose the binary64 global residual after node aggregation."""

        return _checked_fsum(
            "transition conservation residual", self.evaluate_extensive_rates(inventories)
        )

    def to_joint_carrier_arm_authority(self) -> ArmRedistributionAuthority:
        """Adapt only the manufactured arm-rate boundary to the existing oracle."""

        validate_manufactured_adjacent_k_arm_operator(self)
        if self.mode is ArmDispersionMode.STRATIFIED_LIMIT:
            raise ValueError(
                "the stratified limiting oracle cannot be relabelled as an inactive host bed"
            )

        return ArmRedistributionAuthority(
            physical_tray_id=self.mesh.physical_tray_id,
            shaft_history_id=self.accepted_sample.accepted_history_definition_digest,
            shaft_speed_rpm=self.accepted_sample.shaft_speed_rpm,
            tray_geometry_id=self.mesh.tray_geometry_definition_digest,
            closure_id=self.closure.definition_digest,
            source_identity=self.definition_digest,
            active_bed=self.mode is ArmDispersionMode.ACTIVE_BASELINE,
            transitions=self.transitions,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class ManufacturedArmExposure:
    segments: tuple[ManufacturedAcceptedShaftSegment, ...]
    closure: ManufacturedArmClosure
    mode: ArmDispersionMode
    integrated_revolutions: float
    integrated_diffusivity_m2: float

    architecture_only: ClassVar[bool] = True
    accepted_history_authenticated: ClassVar[bool] = False
    production_integrator_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.segments) is not tuple:
            raise TypeError("exposure segments must be an exact immutable tuple")
        if type(self.closure) is not ManufacturedArmClosure:
            raise TypeError("exposure closure must be an exact ManufacturedArmClosure")
        if type(self.mode) is not ArmDispersionMode:
            raise TypeError("exposure mode must be an exact ArmDispersionMode")
        _require_binary64(
            "manufactured arm exposure",
            self.integrated_revolutions,
            self.integrated_diffusivity_m2,
        )
        rebuilt_segments, rebuilt_closure, expected_revolutions, expected_diffusivity = (
            _derive_manufactured_arm_exposure(
                segments=self.segments,
                closure=self.closure,
                mode=self.mode,
            )
        )
        if self.segments != rebuilt_segments or self.closure != rebuilt_closure:
            raise ValueError("exposure inputs differ from their reconstructed definitions")
        if self.integrated_revolutions != expected_revolutions:
            raise ValueError("integrated revolutions disagree with the bound shaft segments")
        if self.integrated_diffusivity_m2 != expected_diffusivity:
            raise ValueError("integrated diffusivity disagrees with the bound segments and closure")

    @property
    def common_shaft_definition_digest(self) -> str:
        return self.segments[0].common_shaft_definition_digest

    @property
    def accepted_history_definition_digest(self) -> str:
        return self.segments[0].accepted_history_definition_digest

    @property
    def closure_definition_digest(self) -> str:
        return self.closure.definition_digest

    @property
    def segment_definition_digests(self) -> tuple[str, ...]:
        return tuple(segment.definition_digest for segment in self.segments)

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="manufactured-arm-exposure",
            parts=(
                self.common_shaft_definition_digest,
                self.accepted_history_definition_digest,
                self.closure_definition_digest,
                self.mode.value,
                len(self.segment_definition_digests),
                *self.segment_definition_digests,
                self.integrated_revolutions,
                self.integrated_diffusivity_m2,
            ),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PR10ConveyanceGateResult:
    accepted_sample: ManufacturedAcceptedShaftSample
    physical_tray_id: str
    draw_geometry_definition_digest: str
    conveyance_domain: PR10ConveyanceDomain
    reported_transfer_kg_s: float
    disposition: PR10ConveyanceDisposition
    pr10_veto: bool
    anomalous_nonzero_report: bool

    architecture_only: ClassVar[bool] = True
    authorizes_transfer: ClassVar[bool] = False
    production_implemented: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.accepted_sample) is not ManufacturedAcceptedShaftSample:
            raise TypeError("accepted_sample must be an exact ManufacturedAcceptedShaftSample")
        if _rebuild_sample(self.accepted_sample) != self.accepted_sample:
            raise ValueError("PR-10 sample differs from its reconstructed definition")
        _require_identity("PR-10 tray identity", self.physical_tray_id)
        _require_sha256("draw_geometry_definition_digest", self.draw_geometry_definition_digest)
        if type(self.conveyance_domain) is not PR10ConveyanceDomain:
            raise TypeError("conveyance_domain must be an exact PR10ConveyanceDomain")
        if type(self.disposition) is not PR10ConveyanceDisposition:
            raise TypeError("disposition must be an exact PR10ConveyanceDisposition")
        if type(self.pr10_veto) is not bool or type(self.anomalous_nonzero_report) is not bool:
            raise TypeError("PR-10 gate booleans must be exact bool values")
        _require_binary64("PR-10 reported transfer", self.reported_transfer_kg_s)
        if self.reported_transfer_kg_s < 0.0:
            raise ValueError("reported transfer must be non-negative")
        if self.conveyance_domain is not PR10ConveyanceDomain.COMMON_SHAFT_SWEEP_CONVEYED:
            expected_disposition = PR10ConveyanceDisposition.OUTSIDE_COMMON_SHAFT_SWEEP_DOMAIN
        elif self.accepted_sample.shaft_speed_rpm == 0.0:
            expected_disposition = PR10ConveyanceDisposition.VETO_ZERO_RPM_COMMON_SHAFT_SWEEP
        else:
            expected_disposition = (
                PR10ConveyanceDisposition.NO_PR10_VETO_POSITIVE_COMMON_SHAFT_SWEEP
            )
        if self.disposition is not expected_disposition:
            raise ValueError("PR-10 disposition disagrees with the bound speed and domain")
        expected_veto = expected_disposition is (
            PR10ConveyanceDisposition.VETO_ZERO_RPM_COMMON_SHAFT_SWEEP
        )
        if self.pr10_veto is not expected_veto:
            raise ValueError("PR-10 gate disposition and veto flag disagree")
        if self.anomalous_nonzero_report is not (
            expected_veto and self.reported_transfer_kg_s > 0.0
        ):
            raise ValueError("PR-10 anomalous-report flag is inconsistent")

    @property
    def common_shaft_definition_digest(self) -> str:
        return self.accepted_sample.common_shaft_definition_digest

    @property
    def accepted_sample_definition_digest(self) -> str:
        return self.accepted_sample.definition_digest

    @property
    def accepted_shaft_speed_rpm(self) -> float:
        return self.accepted_sample.shaft_speed_rpm

    @property
    def definition_digest(self) -> str:
        return _canonical_digest(
            domain="pr10-conveyance-gate-result",
            parts=(
                self.common_shaft_definition_digest,
                self.accepted_sample_definition_digest,
                self.accepted_shaft_speed_rpm,
                self.physical_tray_id,
                self.draw_geometry_definition_digest,
                self.conveyance_domain.value,
                self.reported_transfer_kg_s,
                self.disposition.value,
                self.pr10_veto,
                self.anomalous_nonzero_report,
            ),
        )


def _rebuild_sample(value: ManufacturedAcceptedShaftSample) -> ManufacturedAcceptedShaftSample:
    if type(value) is not ManufacturedAcceptedShaftSample:
        raise TypeError("accepted_sample must be an exact ManufacturedAcceptedShaftSample")
    return ManufacturedAcceptedShaftSample(
        common_shaft_definition_digest=value.common_shaft_definition_digest,
        accepted_history_definition_digest=value.accepted_history_definition_digest,
        sample_id=value.sample_id,
        sample_time_s=value.sample_time_s,
        shaft_speed_rpm=value.shaft_speed_rpm,
    )


def _rebuild_segment(value: ManufacturedAcceptedShaftSegment) -> ManufacturedAcceptedShaftSegment:
    if type(value) is not ManufacturedAcceptedShaftSegment:
        raise TypeError("shaft segments must be exact ManufacturedAcceptedShaftSegment values")
    return ManufacturedAcceptedShaftSegment(
        common_shaft_definition_digest=value.common_shaft_definition_digest,
        accepted_history_definition_digest=value.accepted_history_definition_digest,
        segment_id=value.segment_id,
        start_event_id=value.start_event_id,
        end_event_id=value.end_event_id,
        start_event_kind=value.start_event_kind,
        end_event_kind=value.end_event_kind,
        time_start_s=value.time_start_s,
        time_end_s=value.time_end_s,
        rpm_start=value.rpm_start,
        rpm_end=value.rpm_end,
    )


def _rebuild_closure(value: ManufacturedArmClosure) -> ManufacturedArmClosure:
    if type(value) is not ManufacturedArmClosure:
        raise TypeError("closure must be an exact ManufacturedArmClosure")
    return ManufacturedArmClosure(
        closure_id=value.closure_id,
        closure_source_digest=value.closure_source_digest,
        tray_geometry_definition_digest=value.tray_geometry_definition_digest,
        phi_ex=value.phi_ex,
        delta_z_ex_m=value.delta_z_ex_m,
        f_state=value.f_state,
    )


def _rebuild_cell(value: PhysicalHeightCell) -> PhysicalHeightCell:
    if type(value) is not PhysicalHeightCell:
        raise TypeError("mesh cells must be exact PhysicalHeightCell values")
    return PhysicalHeightCell(
        vertical_layer_id=value.vertical_layer_id,
        z_lower_m=value.z_lower_m,
        z_upper_m=value.z_upper_m,
        dry_matter_holdup_kg=value.dry_matter_holdup_kg,
    )


def _rebuild_face(value: PhysicalHeightFace) -> PhysicalHeightFace:
    if type(value) is not PhysicalHeightFace:
        raise TypeError("mesh faces must be exact PhysicalHeightFace values")
    return PhysicalHeightFace(
        lower_vertical_layer_id=value.lower_vertical_layer_id,
        upper_vertical_layer_id=value.upper_vertical_layer_id,
        face_area_m2=value.face_area_m2,
        dry_matter_bulk_density_kg_m3=value.dry_matter_bulk_density_kg_m3,
        gradient_distance_m=value.gradient_distance_m,
    )


def _rebuild_face_rate(value: ManufacturedArmFaceRate) -> ManufacturedArmFaceRate:
    if type(value) is not ManufacturedArmFaceRate:
        raise TypeError("operator face rates must be exact ManufacturedArmFaceRate values")
    return ManufacturedArmFaceRate(
        lower_vertical_layer_id=value.lower_vertical_layer_id,
        upper_vertical_layer_id=value.upper_vertical_layer_id,
        gradient_distance_m=value.gradient_distance_m,
        diffusivity_m2_s=value.diffusivity_m2_s,
        dry_matter_conductance_kg_s=value.dry_matter_conductance_kg_s,
        lower_to_upper_rate_s_inv=value.lower_to_upper_rate_s_inv,
        upper_to_lower_rate_s_inv=value.upper_to_lower_rate_s_inv,
    )


def _rebuild_transition(value: ArmLayerTransition) -> ArmLayerTransition:
    if type(value) is not ArmLayerTransition:
        raise TypeError("operator transitions must be exact ArmLayerTransition values")
    endpoints = (value.source_vertical_layer_id, value.destination_vertical_layer_id)
    if any(type(endpoint) is not VerticalLayerId for endpoint in endpoints):
        raise TypeError("operator transition endpoints must be exact VerticalLayerId values")
    if any(type(endpoint.value) is not int for endpoint in endpoints):
        raise TypeError("operator transition endpoint values must be exact integers")
    _require_binary64("operator transition rate", value.rate_s_inv)
    return ArmLayerTransition(
        source_vertical_layer_id=value.source_vertical_layer_id,
        destination_vertical_layer_id=value.destination_vertical_layer_id,
        rate_s_inv=value.rate_s_inv,
    )


def _rebuild_mesh(value: PhysicalHeightMesh) -> PhysicalHeightMesh:
    if type(value) is not PhysicalHeightMesh:
        raise TypeError("mesh must be an exact PhysicalHeightMesh")
    if type(value.cells) is not tuple or type(value.faces) is not tuple:
        raise TypeError("mesh cells and faces must remain exact tuples")
    return PhysicalHeightMesh(
        physical_tray_id=value.physical_tray_id,
        mesh_id=value.mesh_id,
        tray_geometry_definition_digest=value.tray_geometry_definition_digest,
        cells=tuple(_rebuild_cell(cell) for cell in value.cells),
        faces=tuple(_rebuild_face(face) for face in value.faces),
    )


def _derive_diffusivity_per_revolution(closure: ManufacturedArmClosure) -> float:
    """Derive the shared blade-scale prefactor with one correctly rounded stage."""

    return _checked_positive_product_ratio(
        "per-revolution arm diffusivity",
        numerators=(
            0.5,
            closure.phi_ex,
            closure.delta_z_ex_m,
            closure.delta_z_ex_m,
            closure.f_state,
        ),
    )


def _derive_diffusivity(
    *,
    closure: ManufacturedArmClosure,
    shaft_speed_rpm: float,
    mode: ArmDispersionMode,
) -> float:
    if shaft_speed_rpm == 0.0 or mode is not ArmDispersionMode.ACTIVE_BASELINE:
        return 0.0
    diffusivity_per_revolution = _derive_diffusivity_per_revolution(closure)
    revolutions_per_second = _checked_positive_product_ratio(
        "instantaneous shaft revolutions per second",
        numerators=(shaft_speed_rpm,),
        denominators=(60.0,),
    )
    return _checked_positive_product_ratio(
        "derived arm diffusivity",
        numerators=(diffusivity_per_revolution, revolutions_per_second),
    )


def _derive_face_records(
    *,
    mesh: PhysicalHeightMesh,
    closure: ManufacturedArmClosure,
    accepted_sample: ManufacturedAcceptedShaftSample,
    mode: ArmDispersionMode,
) -> tuple[tuple[ManufacturedArmFaceRate, ...], tuple[ArmLayerTransition, ...]]:
    """Rebuild every formula-owned face record and directed transition."""

    diffusivity = _derive_diffusivity(
        closure=closure,
        shaft_speed_rpm=accepted_sample.shaft_speed_rpm,
        mode=mode,
    )
    cells = {cell.vertical_layer_id.value: cell for cell in mesh.cells}
    face_rates: list[ManufacturedArmFaceRate] = []
    transitions: list[ArmLayerTransition] = []
    for face in mesh.faces:
        lower = cells[face.lower_vertical_layer_id.value]
        upper = cells[face.upper_vertical_layer_id.value]
        distance = face.gradient_distance_m
        if diffusivity == 0.0:
            conductance = lower_to_upper = upper_to_lower = 0.0
        else:
            conductance = _checked_positive_product_ratio(
                "derived dry-matter face conductance",
                numerators=(
                    diffusivity,
                    face.dry_matter_bulk_density_kg_m3,
                    face.face_area_m2,
                ),
                denominators=(distance,),
            )
            lower_to_upper = _checked_positive_product_ratio(
                "derived lower-to-upper arm rate",
                numerators=(conductance,),
                denominators=(lower.dry_matter_holdup_kg,),
            )
            upper_to_lower = _checked_positive_product_ratio(
                "derived upper-to-lower arm rate",
                numerators=(conductance,),
                denominators=(upper.dry_matter_holdup_kg,),
            )
            transitions.extend(
                (
                    ArmLayerTransition(
                        source_vertical_layer_id=face.lower_vertical_layer_id,
                        destination_vertical_layer_id=face.upper_vertical_layer_id,
                        rate_s_inv=lower_to_upper,
                    ),
                    ArmLayerTransition(
                        source_vertical_layer_id=face.upper_vertical_layer_id,
                        destination_vertical_layer_id=face.lower_vertical_layer_id,
                        rate_s_inv=upper_to_lower,
                    ),
                )
            )
        face_rates.append(
            ManufacturedArmFaceRate(
                lower_vertical_layer_id=face.lower_vertical_layer_id,
                upper_vertical_layer_id=face.upper_vertical_layer_id,
                gradient_distance_m=distance,
                diffusivity_m2_s=diffusivity,
                dry_matter_conductance_kg_s=conductance,
                lower_to_upper_rate_s_inv=lower_to_upper,
                upper_to_lower_rate_s_inv=upper_to_lower,
            )
        )
    canonical_transitions = tuple(
        sorted(
            transitions,
            key=lambda value: (
                value.source_vertical_layer_id.value,
                value.destination_vertical_layer_id.value,
            ),
        )
    )
    outgoing_by_layer: list[list[float]] = [[] for _ in mesh.cells]
    for transition in canonical_transitions:
        outgoing_by_layer[transition.source_vertical_layer_id.value - 1].append(
            transition.rate_s_inv
        )
    for values in outgoing_by_layer:
        outgoing_sum = _checked_fsum(
            "per-cell outgoing arm-rate sum",
            values,
        )
        if outgoing_sum > 0.0:
            _require_positive_normal("per-cell outgoing arm-rate sum", outgoing_sum)
    return tuple(face_rates), canonical_transitions


def build_manufactured_adjacent_k_arm_operator(
    *,
    mesh: PhysicalHeightMesh,
    closure: ManufacturedArmClosure,
    accepted_sample: ManufacturedAcceptedShaftSample,
    mode: ArmDispersionMode,
) -> ManufacturedAdjacentKArmOperator:
    """Build the detached physical-height finite-volume arm-rate oracle."""

    rebuilt_mesh = _rebuild_mesh(mesh)
    rebuilt_closure = _rebuild_closure(closure)
    rebuilt_sample = _rebuild_sample(accepted_sample)
    if type(mode) is not ArmDispersionMode:
        raise TypeError("mode must be an exact ArmDispersionMode")
    if rebuilt_mesh.tray_geometry_definition_digest != (
        rebuilt_closure.tray_geometry_definition_digest
    ):
        raise ValueError("closure and mesh bind different tray geometry identities")
    face_rates, transitions = _derive_face_records(
        mesh=rebuilt_mesh,
        closure=rebuilt_closure,
        accepted_sample=rebuilt_sample,
        mode=mode,
    )
    return ManufacturedAdjacentKArmOperator(
        schema_id=ADJACENT_K_ARM_DISPERSION_SCHEMA_ID,
        mesh=rebuilt_mesh,
        closure=rebuilt_closure,
        accepted_sample=rebuilt_sample,
        mode=mode,
        face_rates=face_rates,
        transitions=transitions,
    )


def validate_manufactured_adjacent_k_arm_operator(
    operator: ManufacturedAdjacentKArmOperator,
) -> str:
    """Recursively reconstruct the operator and return its deterministic identity."""

    if type(operator) is not ManufacturedAdjacentKArmOperator:
        raise TypeError("operator must be an exact ManufacturedAdjacentKArmOperator")
    ManufacturedAdjacentKArmOperator.__post_init__(operator)
    expected = build_manufactured_adjacent_k_arm_operator(
        mesh=operator.mesh,
        closure=operator.closure,
        accepted_sample=operator.accepted_sample,
        mode=operator.mode,
    )
    if operator != expected or operator.definition_digest != expected.definition_digest:
        raise ValueError("adjacent-K arm operator differs from its reconstructed definition")
    return expected.definition_digest


def validate_common_shaft_manufactured_operators(
    operators: Sequence[ManufacturedAdjacentKArmOperator],
) -> None:
    """Require one manufactured shaft sample across distinct physical trays."""

    if type(operators) is not tuple:
        raise TypeError("common-shaft operators must be an exact tuple")
    if not operators:
        raise ValueError("at least one manufactured arm operator is required")
    if any(type(operator) is not ManufacturedAdjacentKArmOperator for operator in operators):
        raise TypeError("common-shaft entries must be exact manufactured arm operators")
    for operator in operators:
        validate_manufactured_adjacent_k_arm_operator(operator)
    tray_ids = tuple(operator.mesh.physical_tray_id for operator in operators)
    if len(set(tray_ids)) != len(tray_ids):
        raise ValueError("common-shaft operators must identify distinct physical trays")
    first = operators[0].accepted_sample
    expected = first.definition_digest
    if any(operator.accepted_sample.definition_digest != expected for operator in operators[1:]):
        raise ValueError("all trays must bind one accepted common-shaft sample")


def _derive_manufactured_arm_exposure(
    *,
    segments: tuple[ManufacturedAcceptedShaftSegment, ...],
    closure: ManufacturedArmClosure,
    mode: ArmDispersionMode,
) -> tuple[
    tuple[ManufacturedAcceptedShaftSegment, ...],
    ManufacturedArmClosure,
    float,
    float,
]:
    """Reconstruct a history chain and derive its analytic trapezoidal integral in binary64."""

    if type(segments) is not tuple:
        raise TypeError("shaft segments must be an exact tuple")
    if not 1 <= len(segments) <= MAX_MANUFACTURED_SEGMENT_COUNT:
        raise ValueError("shaft segment count is outside the manufactured resource domain")
    rebuilt_segments = tuple(_rebuild_segment(segment) for segment in segments)
    rebuilt_closure = _rebuild_closure(closure)
    if type(mode) is not ArmDispersionMode:
        raise TypeError("mode must be an exact ArmDispersionMode")
    segment_ids = tuple(segment.segment_id for segment in rebuilt_segments)
    if len(set(segment_ids)) != len(segment_ids):
        raise ValueError("shaft segment identities must be unique within an exposure")
    first = rebuilt_segments[0]
    event_ids = (first.start_event_id, *(segment.end_event_id for segment in rebuilt_segments))
    if len(set(event_ids)) != len(event_ids):
        raise ValueError("logical shaft boundary event identities must not be reused")
    for left, right in zip(rebuilt_segments[:-1], rebuilt_segments[1:], strict=True):
        if left.time_end_s != right.time_start_s:
            raise ValueError("shaft segments must be exactly contiguous without gaps or overlap")
        if left.end_event_id != right.start_event_id:
            raise ValueError("shaft segment boundary event identities must match exactly")
        if left.end_event_kind is not right.start_event_kind:
            raise ValueError("shaft segment boundary event kinds must match exactly")
        if left.end_event_kind is ShaftBoundaryKind.CONTINUOUS and left.rpm_end != right.rpm_start:
            raise ValueError("continuous shaft boundaries require equal one-sided RPM")
        if left.common_shaft_definition_digest != right.common_shaft_definition_digest:
            raise ValueError("shaft segments bind different common-shaft definitions")
        if left.accepted_history_definition_digest != right.accepted_history_definition_digest:
            raise ValueError("shaft segments bind different accepted histories")
    revolutions_by_segment: list[float] = []
    for segment in rebuilt_segments:
        mean_rpm = _checked_fsum(
            "linear-ramp mean RPM",
            (0.5 * segment.rpm_start, 0.5 * segment.rpm_end),
        )
        if mean_rpm == 0.0:
            revolutions = 0.0
        else:
            revolutions = _checked_positive_product_ratio(
                "integrated shaft revolutions",
                numerators=(mean_rpm, segment.duration_s),
                denominators=(60.0,),
            )
        if not math.isfinite(revolutions) or revolutions < 0.0:
            raise ValueError("integrated shaft revolutions are not representable")
        revolutions_by_segment.append(revolutions)
    integrated_revolutions = _checked_fsum("integrated shaft revolutions", revolutions_by_segment)
    if mode is not ArmDispersionMode.ACTIVE_BASELINE or integrated_revolutions == 0.0:
        integrated_diffusivity = 0.0
    else:
        diffusivity_per_revolution = _derive_diffusivity_per_revolution(rebuilt_closure)
        integrated_diffusivity = _checked_positive_product_ratio(
            "integrated arm diffusivity",
            numerators=(diffusivity_per_revolution, integrated_revolutions),
        )
    return (
        rebuilt_segments,
        rebuilt_closure,
        integrated_revolutions,
        integrated_diffusivity,
    )


def integrate_manufactured_arm_exposure(
    *,
    segments: tuple[ManufacturedAcceptedShaftSegment, ...],
    closure: ManufacturedArmClosure,
    mode: ArmDispersionMode,
) -> ManufacturedArmExposure:
    """Integrate the linear-RPM baseline analytically, with disclosed binary64 split roundoff."""

    rebuilt_segments, rebuilt_closure, integrated_revolutions, integrated_diffusivity = (
        _derive_manufactured_arm_exposure(segments=segments, closure=closure, mode=mode)
    )
    return ManufacturedArmExposure(
        segments=rebuilt_segments,
        closure=rebuilt_closure,
        mode=mode,
        integrated_revolutions=integrated_revolutions,
        integrated_diffusivity_m2=integrated_diffusivity,
    )


def validate_manufactured_arm_exposure(exposure: ManufacturedArmExposure) -> str:
    """Recursively reconstruct an exposure and return its deterministic identity."""

    if type(exposure) is not ManufacturedArmExposure:
        raise TypeError("exposure must be an exact ManufacturedArmExposure")
    ManufacturedArmExposure.__post_init__(exposure)
    expected = integrate_manufactured_arm_exposure(
        segments=exposure.segments,
        closure=exposure.closure,
        mode=exposure.mode,
    )
    if exposure != expected or exposure.definition_digest != expected.definition_digest:
        raise ValueError("manufactured arm exposure differs from its reconstructed definition")
    return expected.definition_digest


def evaluate_pr10_conveyance_gate(
    *,
    accepted_sample: ManufacturedAcceptedShaftSample,
    physical_tray_id: str,
    draw_geometry_definition_digest: str,
    conveyance_domain: PR10ConveyanceDomain,
    reported_transfer_kg_s: float,
) -> PR10ConveyanceGateResult:
    """Classify the exact-zero sweep veto without accepting or scaling flow."""

    sample = _rebuild_sample(accepted_sample)
    _require_identity("PR-10 physical tray identity", physical_tray_id)
    _require_sha256("draw_geometry_definition_digest", draw_geometry_definition_digest)
    if type(conveyance_domain) is not PR10ConveyanceDomain:
        raise TypeError("conveyance_domain must be an exact PR10ConveyanceDomain")
    _require_binary64("reported transfer", reported_transfer_kg_s)
    if reported_transfer_kg_s < 0.0:
        raise ValueError("reported transfer must be non-negative")
    if conveyance_domain is not PR10ConveyanceDomain.COMMON_SHAFT_SWEEP_CONVEYED:
        disposition = PR10ConveyanceDisposition.OUTSIDE_COMMON_SHAFT_SWEEP_DOMAIN
    elif sample.shaft_speed_rpm == 0.0:
        disposition = PR10ConveyanceDisposition.VETO_ZERO_RPM_COMMON_SHAFT_SWEEP
    else:
        disposition = PR10ConveyanceDisposition.NO_PR10_VETO_POSITIVE_COMMON_SHAFT_SWEEP
    veto = disposition is PR10ConveyanceDisposition.VETO_ZERO_RPM_COMMON_SHAFT_SWEEP
    return PR10ConveyanceGateResult(
        accepted_sample=sample,
        physical_tray_id=physical_tray_id,
        draw_geometry_definition_digest=draw_geometry_definition_digest,
        conveyance_domain=conveyance_domain,
        reported_transfer_kg_s=reported_transfer_kg_s,
        disposition=disposition,
        pr10_veto=veto,
        anomalous_nonzero_report=veto and reported_transfer_kg_s > 0.0,
    )


def validate_pr10_conveyance_gate_result(result: PR10ConveyanceGateResult) -> str:
    """Recursively reconstruct a PR-10 veto result and return its identity."""

    if type(result) is not PR10ConveyanceGateResult:
        raise TypeError("result must be an exact PR10ConveyanceGateResult")
    PR10ConveyanceGateResult.__post_init__(result)
    expected = evaluate_pr10_conveyance_gate(
        accepted_sample=result.accepted_sample,
        physical_tray_id=result.physical_tray_id,
        draw_geometry_definition_digest=result.draw_geometry_definition_digest,
        conveyance_domain=result.conveyance_domain,
        reported_transfer_kg_s=result.reported_transfer_kg_s,
    )
    if result != expected or result.definition_digest != expected.definition_digest:
        raise ValueError("PR-10 gate result differs from its reconstructed definition")
    return expected.definition_digest


__all__ = [
    "ADJACENT_K_ARM_DISPERSION_SCHEMA_ID",
    "MAX_DENSE_GENERATOR_LAYER_COUNT",
    "MAX_IDENTITY_UTF8_BYTES",
    "MAX_MANUFACTURED_LAYER_COUNT",
    "MAX_MANUFACTURED_SEGMENT_COUNT",
    "ArmDispersionMode",
    "ManufacturedAcceptedShaftSample",
    "ManufacturedAcceptedShaftSegment",
    "ManufacturedAdjacentKArmOperator",
    "ManufacturedArmClosure",
    "ManufacturedArmExposure",
    "ManufacturedArmFaceRate",
    "PR10ConveyanceDisposition",
    "PR10ConveyanceDomain",
    "PR10ConveyanceGateResult",
    "PhysicalHeightCell",
    "PhysicalHeightFace",
    "PhysicalHeightMesh",
    "ShaftBoundaryKind",
    "build_manufactured_adjacent_k_arm_operator",
    "evaluate_pr10_conveyance_gate",
    "integrate_manufactured_arm_exposure",
    "validate_common_shaft_manufactured_operators",
    "validate_manufactured_adjacent_k_arm_operator",
    "validate_manufactured_arm_exposure",
    "validate_pr10_conveyance_gate_result",
]
