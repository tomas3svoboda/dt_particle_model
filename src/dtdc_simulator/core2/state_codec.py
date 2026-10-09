"""Restricted deterministic state codec for future Core2 particle payloads.

This module implements a deliberately small RFC 8949 application profile.  It
is not a general CBOR implementation and it does not claim conformance with the
RFC's core deterministic encoding requirements.  Fields declared ``FLOAT64``
always use binary64, while integer fields remain distinct shortest-form CBOR
integers.  This is an application-specific deterministic profile informed by
RFC 8949 section 4.2; it does not adopt section 4.2.2 Rule 3 globally.

Every encoded state is one definite-length array::

    [profile_magic, profile_version, schema_id, schema_version, variant_tag,
     schema_definition_sha256, field_values]

The registry supplies a fixed, ordered schema for ``field_values``.  Maps,
semantic tags, indefinite-length items, null/undefined values, float16,
float32, mutable containers, subclasses, non-finite values, non-NFC text, and
trailing bytes are rejected.  Finite binary64 subnormals are preserved without
flush-to-zero.  Signed-zero bits are preserved for opaque state; fields
explicitly declared as non-negative reject negative zero.

Canonical bytes are the identity authority.  The state-identity SHA-256 is
applied only after a full decode/re-encode canonicality check.  This foundation does not yet select
or authenticate a complete real-particle variant registry, a production
adapter, persistence, transport, or any physical model.  Those claims remain
false until every accepted particle state is inventoried and qualified.
"""

from __future__ import annotations

import enum
import hashlib
import math
import re
import struct
import unicodedata
from dataclasses import InitVar, dataclass
from functools import lru_cache
from typing import ClassVar, TypeAlias

PROFILE_MAGIC = "GT-PS-2-STATE-CODEC"
PROFILE_VERSION = 1
MAXIMUM_CANONICAL_STATE_BYTES = 16 * 1024 * 1024
MAXIMUM_SCALAR_BYTES = 1024 * 1024
MAXIMUM_ARRAY_ITEMS = 1024 * 1024
MAXIMUM_RECORD_FIELDS = 4096
MAXIMUM_NESTING_DEPTH = 32
MAXIMUM_TOTAL_ITEMS = 262_144
MAXIMUM_UINT64 = (1 << 64) - 1
MINIMUM_CBOR_INTEGER = -(1 << 64)

_SCHEMA_ID = re.compile(r"[A-Za-z][A-Za-z0-9._:/-]{0,127}\Z")
_ISSUE_TOKEN = object()


class StateCodecError(ValueError):
    """Base exception for malformed or unsupported codec values."""


class StateCodecDecodeError(StateCodecError):
    """Raised when encoded bytes are malformed, noncanonical, or unsupported."""


class ScalarKind(enum.Enum):
    """Scalar types admitted by the restricted application profile."""

    UTF8 = "utf8"
    BYTES = "bytes"
    UINT = "uint"
    SINT = "sint"
    BOOL = "bool"
    FLOAT64 = "float64"


def _require_exact_integer(name: str, value: int, *, minimum: int, maximum: int) -> None:
    if type(value) is not int:
        raise TypeError(f"{name} must be an exact integer")
    if not minimum <= value <= maximum:
        raise ValueError(f"{name} must be in [{minimum}, {maximum}]")


def _require_schema_identifier(name: str, value: str) -> None:
    if type(value) is not str or _SCHEMA_ID.fullmatch(value) is None:
        raise ValueError(f"{name} must match the ASCII schema-identifier application profile")


def _require_nfc_text(name: str, value: str) -> bytes:
    if type(value) is not str:
        raise TypeError(f"{name} must be an exact str")
    if unicodedata.normalize("NFC", value) != value:
        raise ValueError(f"{name} must already be NFC-normalized")
    try:
        return value.encode("utf-8", errors="strict")
    except UnicodeEncodeError as error:
        raise ValueError(f"{name} contains text outside strict UTF-8") from error


def _is_negative_zero(value: float) -> bool:
    return struct.pack(">d", value) == b"\x80\x00\x00\x00\x00\x00\x00\x00"


def _require_finite_binary64(name: str, value: float, *, nonnegative: bool) -> None:
    if type(value) is not float:
        raise TypeError(f"{name} must be an exact binary64 float")
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if nonnegative and (value < 0.0 or _is_negative_zero(value)):
        raise ValueError(f"{name} must be non-negative with canonical positive zero")


def _framed_digest(domain: str, parts: tuple[bytes, ...]) -> str:
    digest = hashlib.sha256()
    domain_bytes = domain.encode("ascii")
    for part in (domain_bytes, *parts):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True)
class ScalarField:
    """One exact scalar field in a versioned state schema."""

    name: str
    kind: ScalarKind
    nonnegative: bool = False
    maximum_length: int | None = None

    def __post_init__(self) -> None:
        _require_schema_identifier("scalar field name", self.name)
        if type(self.kind) is not ScalarKind:
            raise TypeError("scalar field kind must be an exact ScalarKind")
        if type(self.nonnegative) is not bool:
            raise TypeError("scalar field nonnegative marker must be an exact bool")
        if self.nonnegative and self.kind is not ScalarKind.FLOAT64:
            raise ValueError("nonnegative is defined only for binary64 fields")
        if self.kind in (ScalarKind.UTF8, ScalarKind.BYTES):
            maximum_length = (
                MAXIMUM_SCALAR_BYTES if self.maximum_length is None else self.maximum_length
            )
            _require_exact_integer(
                "scalar field maximum_length",
                maximum_length,
                minimum=0,
                maximum=MAXIMUM_SCALAR_BYTES,
            )
        elif self.maximum_length is not None:
            raise ValueError("maximum_length is valid only for text and byte fields")

    @property
    def effective_maximum_length(self) -> int | None:
        if self.kind in (ScalarKind.UTF8, ScalarKind.BYTES):
            return MAXIMUM_SCALAR_BYTES if self.maximum_length is None else self.maximum_length
        return None


@dataclass(frozen=True, slots=True)
class ArrayField:
    """One immutable tuple field whose item schema is recursively fixed."""

    name: str
    item: FieldSpec
    length: int | None = None
    maximum_length: int = MAXIMUM_ARRAY_ITEMS

    def __post_init__(self) -> None:
        _require_schema_identifier("array field name", self.name)
        if type(self.item) not in (
            ScalarField,
            ArrayField,
            RecordField,
            EnumField,
            VariantField,
        ):
            raise TypeError("array item must be an exact codec field specification")
        _require_exact_integer(
            "array maximum_length",
            self.maximum_length,
            minimum=0,
            maximum=MAXIMUM_ARRAY_ITEMS,
        )
        if self.length is not None:
            _require_exact_integer(
                "array length",
                self.length,
                minimum=0,
                maximum=self.maximum_length,
            )


@dataclass(frozen=True, slots=True)
class RecordField:
    """One fixed heterogeneous nested record encoded as a CBOR array."""

    name: str
    fields: tuple[FieldSpec, ...]

    def __post_init__(self) -> None:
        _require_schema_identifier("record field name", self.name)
        if type(self.fields) is not tuple:
            raise TypeError("record fields must be an exact immutable tuple")
        if not self.fields:
            raise ValueError("record must declare at least one field")
        if len(self.fields) > MAXIMUM_RECORD_FIELDS:
            raise ValueError("record exceeds the field-count limit")
        if any(
            type(field) not in (ScalarField, ArrayField, RecordField, EnumField, VariantField)
            for field in self.fields
        ):
            raise TypeError("record must contain exact codec field specifications")
        names = tuple(field.name for field in self.fields)
        if len(names) != len(set(names)):
            raise ValueError("record field names must be unique")


@dataclass(frozen=True, slots=True)
class EnumField:
    """One exact schema-bound textual enumeration."""

    name: str
    tokens: tuple[str, ...]

    def __post_init__(self) -> None:
        _require_schema_identifier("enum field name", self.name)
        if type(self.tokens) is not tuple:
            raise TypeError("enum tokens must be an exact immutable tuple")
        if not self.tokens:
            raise ValueError("enum must declare at least one token")
        for token in self.tokens:
            _require_schema_identifier("enum token", token)
        if len(self.tokens) != len(set(self.tokens)):
            raise ValueError("enum tokens must be unique")
        if self.tokens != tuple(sorted(self.tokens)):
            raise ValueError("enum tokens must be in canonical lexical order")


@dataclass(frozen=True, slots=True)
class VariantArm:
    """One explicitly tagged arm of a nested discriminated union."""

    tag: int
    name: str
    fields: tuple[FieldSpec, ...] = ()

    def __post_init__(self) -> None:
        _require_exact_integer("variant arm tag", self.tag, minimum=0, maximum=MAXIMUM_UINT64)
        _require_schema_identifier("variant arm name", self.name)
        if type(self.fields) is not tuple:
            raise TypeError("variant arm fields must be an exact immutable tuple")
        if len(self.fields) > MAXIMUM_RECORD_FIELDS:
            raise ValueError("variant arm exceeds the field-count limit")
        if any(
            type(field) not in (ScalarField, ArrayField, RecordField, EnumField, VariantField)
            for field in self.fields
        ):
            raise TypeError("variant arm must contain exact codec field specifications")
        names = tuple(field.name for field in self.fields)
        if len(names) != len(set(names)):
            raise ValueError("variant arm field names must be unique")


@dataclass(frozen=True, slots=True)
class VariantField:
    """A nested union encoded as ``[tag, fixed-arm-field-array]``."""

    name: str
    arms: tuple[VariantArm, ...]

    def __post_init__(self) -> None:
        _require_schema_identifier("variant field name", self.name)
        if type(self.arms) is not tuple:
            raise TypeError("variant arms must be an exact immutable tuple")
        if not self.arms:
            raise ValueError("variant field must declare at least one arm")
        if any(type(arm) is not VariantArm for arm in self.arms):
            raise TypeError("variant field must contain exact VariantArm values")
        tags = tuple(arm.tag for arm in self.arms)
        names = tuple(arm.name for arm in self.arms)
        if len(tags) != len(set(tags)) or len(names) != len(set(names)):
            raise ValueError("variant arm tags and names must be unique")
        if tags != tuple(sorted(tags)):
            raise ValueError("variant arms must be in canonical tag order")

    def resolve(self, tag: int) -> VariantArm:
        for arm in self.arms:
            if arm.tag == tag:
                return arm
        raise StateCodecError(f"unknown variant tag {tag} for {self.name}")


FieldSpec: TypeAlias = ScalarField | ArrayField | RecordField | EnumField | VariantField
CodecScalar: TypeAlias = str | bytes | int | bool | float
CodecValue: TypeAlias = CodecScalar | tuple["CodecValue", ...]


def _validate_field_tree(
    field: FieldSpec,
    *,
    depth: int,
    active: frozenset[int],
) -> None:
    if depth > MAXIMUM_NESTING_DEPTH:
        raise ValueError("codec schema exceeds the nesting-depth limit")
    if type(field) not in (ScalarField, ArrayField, RecordField, EnumField, VariantField):
        raise TypeError("schema contains a non-codec field specification")
    if id(field) in active:
        raise ValueError("codec schema contains a recursive field cycle")
    if type(field) is ScalarField:
        field.__post_init__()
        return
    if type(field) is EnumField:
        field.__post_init__()
        return
    next_active = active | {id(field)}
    if type(field) is ArrayField:
        field.__post_init__()
        _validate_field_tree(field.item, depth=depth + 1, active=next_active)
        return
    if type(field) is RecordField:
        field.__post_init__()
        for child in field.fields:
            _validate_field_tree(child, depth=depth + 1, active=next_active)
        return
    field.__post_init__()
    for arm in field.arms:
        arm.__post_init__()
        for child in arm.fields:
            _validate_field_tree(child, depth=depth + 1, active=next_active)


def _minimum_field_items(field: FieldSpec) -> int:
    if type(field) in (ScalarField, EnumField):
        return 1
    if type(field) is ArrayField:
        count = 0 if field.length is None else field.length
        return 1 + count * _minimum_field_items(field.item)
    if type(field) is RecordField:
        return 1 + sum(_minimum_field_items(child) for child in field.fields)
    if type(field) is VariantField:
        return 3 + min(
            sum(_minimum_field_items(child) for child in arm.fields) for arm in field.arms
        )
    raise TypeError("field must be an exact codec field specification")


def _field_definition_parts(field: FieldSpec) -> tuple[bytes, ...]:
    if type(field) is ScalarField:
        maximum = field.effective_maximum_length
        return (
            b"scalar",
            field.name.encode("ascii"),
            field.kind.value.encode("ascii"),
            b"1" if field.nonnegative else b"0",
            b"none" if maximum is None else maximum.to_bytes(8, "big"),
        )
    if type(field) is ArrayField:
        length = b"variable" if field.length is None else field.length.to_bytes(8, "big")
        nested = _field_definition_parts(field.item)
        nested_digest = _framed_digest("GT-PS-2/state-codec/field/v1", nested)
        return (
            b"array",
            field.name.encode("ascii"),
            length,
            field.maximum_length.to_bytes(8, "big"),
            nested_digest.encode("ascii"),
        )
    if type(field) is RecordField:
        child_digests = tuple(
            _framed_digest(
                "GT-PS-2/state-codec/field/v1",
                _field_definition_parts(child),
            ).encode("ascii")
            for child in field.fields
        )
        return (b"record", field.name.encode("ascii"), *child_digests)
    if type(field) is EnumField:
        return (
            b"enum",
            field.name.encode("ascii"),
            *(token.encode("ascii") for token in field.tokens),
        )
    if type(field) is VariantField:
        arm_digests = []
        for arm in field.arms:
            child_digests = tuple(
                _framed_digest(
                    "GT-PS-2/state-codec/field/v1",
                    _field_definition_parts(child),
                ).encode("ascii")
                for child in arm.fields
            )
            arm_digests.append(
                _framed_digest(
                    "GT-PS-2/state-codec/variant-arm/v1",
                    (
                        arm.tag.to_bytes(8, "big"),
                        arm.name.encode("ascii"),
                        *child_digests,
                    ),
                ).encode("ascii")
            )
        return (b"variant", field.name.encode("ascii"), *arm_digests)
    raise TypeError("field must be an exact codec field specification")


@dataclass(frozen=True, slots=True)
class StateSchema:
    """One fixed application schema and variant tag."""

    schema_id: str
    schema_version: int
    variant_tag: int
    fields: tuple[FieldSpec, ...]
    maximum_encoded_bytes: int = MAXIMUM_CANONICAL_STATE_BYTES
    # Counts the decoded state-value tree.  Fixed seven-item envelope overhead
    # is bounded separately by maximum_encoded_bytes and is not counted here.
    maximum_items: int = MAXIMUM_TOTAL_ITEMS

    architecture_only: ClassVar[bool] = True
    restricted_application_profile_selected: ClassVar[bool] = True
    core_deterministic_encoding_claimed: ClassVar[bool] = False
    production_variant_registry_complete: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_schema_identifier("schema_id", self.schema_id)
        _require_exact_integer(
            "schema_version",
            self.schema_version,
            minimum=1,
            maximum=MAXIMUM_UINT64,
        )
        _require_exact_integer(
            "variant_tag",
            self.variant_tag,
            minimum=0,
            maximum=MAXIMUM_UINT64,
        )
        if type(self.fields) is not tuple:
            raise TypeError("schema fields must be an exact immutable tuple")
        if not self.fields:
            raise ValueError("schema must declare at least one field")
        if len(self.fields) > MAXIMUM_RECORD_FIELDS:
            raise ValueError("schema exceeds the field-count limit")
        if any(
            type(field) not in (ScalarField, ArrayField, RecordField, EnumField, VariantField)
            for field in self.fields
        ):
            raise TypeError("schema fields must contain exact field specifications")
        names = tuple(field.name for field in self.fields)
        if len(names) != len(set(names)):
            raise ValueError("schema field names must be unique")
        for field in self.fields:
            _validate_field_tree(field, depth=1, active=frozenset())
        _require_exact_integer(
            "maximum_encoded_bytes",
            self.maximum_encoded_bytes,
            minimum=1,
            maximum=MAXIMUM_CANONICAL_STATE_BYTES,
        )
        minimum_items = sum(_minimum_field_items(field) for field in self.fields)
        _require_exact_integer(
            "maximum_items",
            self.maximum_items,
            minimum=minimum_items,
            maximum=MAXIMUM_TOTAL_ITEMS,
        )

    @property
    def identity(self) -> tuple[str, int, int]:
        return self.schema_id, self.schema_version, self.variant_tag

    @property
    def definition_digest(self) -> str:
        # Owner-approved B7 sub-item (2026-08-31): the digest is memoized
        # VALUE-KEYED through the frozen dataclass's own field-wise
        # hash/eq, so a mutated schema hashes differently and recomputes —
        # the anti-mutation property of the walk is preserved while the
        # constant-per-schema SHA walk (measured 6.4k repeats per
        # interval) collapses to one.
        return _schema_definition_digest(self)


@lru_cache(maxsize=4096)
def _schema_definition_digest(schema: StateSchema) -> str:
    field_digests = tuple(
        _framed_digest("GT-PS-2/state-codec/field/v1", _field_definition_parts(field)).encode(
            "ascii"
        )
        for field in schema.fields
    )
    return _framed_digest(
        "GT-PS-2/state-codec/schema/v1",
        (
            PROFILE_MAGIC.encode("ascii"),
            PROFILE_VERSION.to_bytes(8, "big"),
            schema.schema_id.encode("ascii"),
            schema.schema_version.to_bytes(8, "big"),
            schema.variant_tag.to_bytes(8, "big"),
            schema.maximum_encoded_bytes.to_bytes(8, "big"),
            schema.maximum_items.to_bytes(8, "big"),
            *field_digests,
        ),
    )


@dataclass(frozen=True, slots=True)
class StateSchemaRegistry:
    """Immutable registry of the exact variants accepted at one boundary."""

    registry_id: str
    schemas: tuple[StateSchema, ...]

    architecture_only: ClassVar[bool] = True
    source_authenticated: ClassVar[bool] = False
    production_variant_registry_complete: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_schema_identifier("registry_id", self.registry_id)
        if type(self.schemas) is not tuple:
            raise TypeError("registry schemas must be an exact immutable tuple")
        if not self.schemas:
            raise ValueError("registry must contain at least one schema")
        if any(type(schema) is not StateSchema for schema in self.schemas):
            raise TypeError("registry must contain exact StateSchema values")
        identities = tuple(schema.identity for schema in self.schemas)
        if len(identities) != len(set(identities)):
            raise ValueError("registry schema identities must be unique")
        if identities != tuple(sorted(identities)):
            raise ValueError("registry schemas must be in canonical identity order")

    def resolve(self, schema_id: str, schema_version: int, variant_tag: int) -> StateSchema:
        identity = (schema_id, schema_version, variant_tag)
        for schema in self.schemas:
            if schema.identity == identity:
                return schema
        raise StateCodecDecodeError(f"unknown state schema identity {identity!r}")

    @property
    def definition_digest(self) -> str:
        return _framed_digest(
            "GT-PS-2/state-codec/registry/v1",
            (
                self.registry_id.encode("ascii"),
                *(schema.definition_digest.encode("ascii") for schema in self.schemas),
            ),
        )


@dataclass(frozen=True, slots=True)
class DecodedState:
    """Evaluator-issued immutable decoded state and its canonical identity."""

    schema: StateSchema
    values: tuple[CodecValue, ...]
    canonical_bytes: bytes
    identity_digest: str
    _issue_token: InitVar[object | None] = None

    architecture_only: ClassVar[bool] = True
    canonical_bytes_authoritative: ClassVar[bool] = True
    source_authenticated: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False

    def __post_init__(self, _issue_token: object | None) -> None:
        if _issue_token is not _ISSUE_TOKEN:
            raise TypeError("DecodedState values are issued only by the codec")
        if type(self.schema) is not StateSchema:
            raise TypeError("decoded state schema must be an exact StateSchema")
        if type(self.values) is not tuple:
            raise TypeError("decoded state values must be an exact tuple")
        if type(self.canonical_bytes) is not bytes:
            raise TypeError("canonical state bytes must be exact bytes")
        if type(self.identity_digest) is not str:
            raise TypeError("state identity digest must be an exact str")

    def __copy__(self) -> None:
        raise TypeError("DecodedState authority views cannot be copied")

    def __deepcopy__(self, memo: object) -> None:
        del memo
        raise TypeError("DecodedState authority views cannot be deep-copied")

    def __reduce_ex__(self, protocol: int) -> None:
        del protocol
        raise TypeError("DecodedState authority views cannot be serialized")


def _encode_head(major: int, value: int) -> bytes:
    _require_exact_integer("CBOR major type", major, minimum=0, maximum=7)
    _require_exact_integer("CBOR unsigned argument", value, minimum=0, maximum=MAXIMUM_UINT64)
    initial = major << 5
    if value < 24:
        return bytes((initial | value,))
    if value <= 0xFF:
        return bytes((initial | 24, value))
    if value <= 0xFFFF:
        return bytes((initial | 25,)) + value.to_bytes(2, "big")
    if value <= 0xFFFFFFFF:
        return bytes((initial | 26,)) + value.to_bytes(4, "big")
    return bytes((initial | 27,)) + value.to_bytes(8, "big")


class _Encoder:
    __slots__ = ("data", "item_count", "maximum_bytes", "maximum_items")

    def __init__(self, *, maximum_bytes: int, maximum_items: int) -> None:
        self.data = bytearray()
        self.item_count = 0
        self.maximum_bytes = maximum_bytes
        self.maximum_items = maximum_items

    def write(self, payload: bytes) -> None:
        if len(payload) > self.maximum_bytes - len(self.data):
            raise ValueError("encoded state exceeds the schema byte limit")
        self.data.extend(payload)

    def consume_item(self, *, count: int = 1) -> None:
        if count > self.maximum_items - self.item_count:
            raise ValueError("state value exceeds the schema aggregate-item limit")
        self.item_count += count


def _write_text(
    encoder: _Encoder,
    value: str,
    *,
    context: str,
    maximum_length: int,
) -> None:
    if type(value) is not str:
        raise TypeError(f"{context} must be an exact str")
    if len(value) > maximum_length:
        raise ValueError(f"{context} exceeds its UTF-8 byte limit")
    payload = _require_nfc_text(context, value)
    if len(payload) > maximum_length:
        raise ValueError(f"{context} exceeds its UTF-8 byte limit")
    encoder.write(_encode_head(3, len(payload)))
    encoder.write(payload)


def _write_bytes(
    encoder: _Encoder,
    value: bytes,
    *,
    context: str,
    maximum_length: int,
) -> None:
    if type(value) is not bytes:
        raise TypeError(f"{context} must be exact immutable bytes")
    if len(value) > maximum_length:
        raise ValueError(f"{context} exceeds its byte limit")
    encoder.write(_encode_head(2, len(value)))
    encoder.write(value)


def _validate_array_length(field: ArrayField, value: tuple[CodecValue, ...]) -> None:
    if type(value) is not tuple:
        raise TypeError(f"{field.name} must be an exact immutable tuple")
    if len(value) > field.maximum_length:
        raise ValueError(f"{field.name} exceeds its array length limit")
    if field.length is not None and len(value) != field.length:
        raise ValueError(f"{field.name} must contain exactly {field.length} items")


def _encode_field(
    encoder: _Encoder,
    field: FieldSpec,
    value: CodecValue,
    *,
    path: str,
    depth: int,
) -> None:
    if depth > MAXIMUM_NESTING_DEPTH:
        raise ValueError("state value exceeds the nesting-depth limit")
    encoder.consume_item()
    if type(field) is ScalarField:
        if field.kind is ScalarKind.UTF8:
            assert field.effective_maximum_length is not None
            _write_text(
                encoder,
                value,  # type: ignore[arg-type]
                context=path,
                maximum_length=field.effective_maximum_length,
            )
            return
        if field.kind is ScalarKind.BYTES:
            assert field.effective_maximum_length is not None
            _write_bytes(
                encoder,
                value,  # type: ignore[arg-type]
                context=path,
                maximum_length=field.effective_maximum_length,
            )
            return
        if field.kind is ScalarKind.UINT:
            _require_exact_integer(path, value, minimum=0, maximum=MAXIMUM_UINT64)  # type: ignore[arg-type]
            encoder.write(_encode_head(0, value))  # type: ignore[arg-type]
            return
        if field.kind is ScalarKind.SINT:
            _require_exact_integer(
                path,
                value,  # type: ignore[arg-type]
                minimum=MINIMUM_CBOR_INTEGER,
                maximum=MAXIMUM_UINT64,
            )
            if value >= 0:  # type: ignore[operator]
                encoder.write(_encode_head(0, value))  # type: ignore[arg-type]
            else:
                encoder.write(_encode_head(1, -1 - value))  # type: ignore[operator]
            return
        if field.kind is ScalarKind.BOOL:
            if type(value) is not bool:
                raise TypeError(f"{path} must be an exact bool")
            encoder.write(b"\xf5" if value else b"\xf4")
            return
        if field.kind is ScalarKind.FLOAT64:
            _require_finite_binary64(path, value, nonnegative=field.nonnegative)  # type: ignore[arg-type]
            encoder.write(b"\xfb")
            encoder.write(struct.pack(">d", value))
            return
        raise AssertionError("unreachable scalar kind")

    if type(field) is EnumField:
        if type(value) is not str:
            raise TypeError(f"{path} must be an exact enum-token str")
        if value not in field.tokens:
            raise ValueError(f"{path} contains an unknown enum token")
        _write_text(encoder, value, context=path, maximum_length=128)
        return

    if type(field) is ArrayField:
        _validate_array_length(field, value)  # type: ignore[arg-type]
        encoder.write(_encode_head(4, len(value)))  # type: ignore[arg-type]
        for index, item in enumerate(value):  # type: ignore[arg-type]
            _encode_field(
                encoder,
                field.item,
                item,
                path=f"{path}[{index}]",
                depth=depth + 1,
            )
        return

    if type(field) is RecordField:
        if type(value) is not tuple:
            raise TypeError(f"{path} must be an exact immutable tuple")
        if len(value) != len(field.fields):
            raise ValueError(f"{path} differs from its fixed record field count")
        encoder.write(_encode_head(4, len(field.fields)))
        for child, child_value in zip(field.fields, value, strict=True):
            _encode_field(
                encoder,
                child,
                child_value,
                path=f"{path}.{child.name}",
                depth=depth + 1,
            )
        return

    if type(field) is VariantField:
        if type(value) is not tuple or len(value) != 2:
            raise TypeError(f"{path} must be the exact tuple (variant_tag, arm_values)")
        tag, arm_values = value
        _require_exact_integer(
            f"{path} variant tag",
            tag,  # type: ignore[arg-type]
            minimum=0,
            maximum=MAXIMUM_UINT64,
        )
        if type(arm_values) is not tuple:
            raise TypeError(f"{path} variant arm values must be an exact tuple")
        arm = field.resolve(tag)  # type: ignore[arg-type]
        if len(arm_values) != len(arm.fields):
            raise ValueError(f"{path} differs from its selected variant-arm field count")
        encoder.consume_item(count=2)
        encoder.write(_encode_head(4, 2))
        encoder.write(_encode_head(0, tag))  # type: ignore[arg-type]
        encoder.write(_encode_head(4, len(arm.fields)))
        for child, child_value in zip(arm.fields, arm_values, strict=True):
            _encode_field(
                encoder,
                child,
                child_value,
                path=f"{path}.{arm.name}.{child.name}",
                depth=depth + 1,
            )
        return

    raise TypeError("field must be an exact codec field specification")


def encode_state(schema: StateSchema, values: tuple[CodecValue, ...]) -> bytes:
    """Encode one state under the restricted fixed-binary64 profile."""

    if type(schema) is not StateSchema:
        raise TypeError("schema must be an exact StateSchema")
    schema.__post_init__()
    if type(values) is not tuple:
        raise TypeError("state values must be an exact immutable tuple")
    if len(values) != len(schema.fields):
        raise ValueError("state value count differs from the fixed schema")

    encoder = _Encoder(
        maximum_bytes=schema.maximum_encoded_bytes,
        maximum_items=schema.maximum_items,
    )
    schema_definition = bytes.fromhex(schema.definition_digest.removeprefix("sha256:"))
    encoder.write(_encode_head(4, 7))
    _write_text(encoder, PROFILE_MAGIC, context="profile magic", maximum_length=128)
    encoder.write(_encode_head(0, PROFILE_VERSION))
    _write_text(encoder, schema.schema_id, context="schema_id", maximum_length=128)
    encoder.write(_encode_head(0, schema.schema_version))
    encoder.write(_encode_head(0, schema.variant_tag))
    _write_bytes(
        encoder,
        schema_definition,
        context="schema definition digest",
        maximum_length=32,
    )
    encoder.write(_encode_head(4, len(values)))
    for field, value in zip(schema.fields, values, strict=True):
        _encode_field(encoder, field, value, path=field.name, depth=1)
    return bytes(encoder.data)


class _Reader:
    __slots__ = ("data", "position")

    def __init__(self, data: bytes) -> None:
        if type(data) is not bytes:
            raise TypeError("encoded state must be exact immutable bytes")
        if not data:
            raise StateCodecDecodeError("encoded state must not be empty")
        if len(data) > MAXIMUM_CANONICAL_STATE_BYTES:
            raise StateCodecDecodeError("encoded state exceeds the profile byte limit")
        self.data = data
        self.position = 0

    def read(self, count: int, *, context: str) -> bytes:
        end = self.position + count
        if end > len(self.data):
            raise StateCodecDecodeError(f"truncated {context}")
        payload = self.data[self.position : end]
        self.position = end
        return payload

    def byte(self, *, context: str) -> int:
        return self.read(1, context=context)[0]

    def head(self, expected_major: int, *, context: str) -> int:
        initial = self.byte(context=context)
        major = initial >> 5
        additional = initial & 0x1F
        if major != expected_major:
            raise StateCodecDecodeError(
                f"{context} has CBOR major type {major}, expected {expected_major}"
            )
        if additional < 24:
            return additional
        width_by_additional = {24: 1, 25: 2, 26: 4, 27: 8}
        width = width_by_additional.get(additional)
        if width is None:
            raise StateCodecDecodeError(f"{context} uses an unsupported or indefinite head")
        value = int.from_bytes(self.read(width, context=context), "big")
        lower_bounds = {1: 24, 2: 0x100, 4: 0x10000, 8: 0x100000000}
        if value < lower_bounds[width]:
            raise StateCodecDecodeError(f"{context} uses a non-shortest integer or length")
        return value

    def text(self, *, context: str, maximum_length: int) -> str:
        length = self.head(3, context=context)
        if length > maximum_length:
            raise StateCodecDecodeError(f"{context} exceeds its UTF-8 byte limit")
        payload = self.read(length, context=context)
        try:
            value = payload.decode("utf-8", errors="strict")
        except UnicodeDecodeError as error:
            raise StateCodecDecodeError(f"{context} is not strict UTF-8") from error
        if unicodedata.normalize("NFC", value) != value:
            raise StateCodecDecodeError(f"{context} is not NFC-normalized")
        return value


class _DecodeBudget:
    __slots__ = ("item_count", "maximum_items")

    def __init__(self, *, maximum_items: int) -> None:
        self.item_count = 0
        self.maximum_items = maximum_items

    def consume_item(self, *, count: int = 1) -> None:
        if count > self.maximum_items - self.item_count:
            raise StateCodecDecodeError("state value exceeds the schema aggregate-item limit")
        self.item_count += count


def _decode_field(
    reader: _Reader,
    budget: _DecodeBudget,
    field: FieldSpec,
    *,
    path: str,
    depth: int,
) -> CodecValue:
    if depth > MAXIMUM_NESTING_DEPTH:
        raise StateCodecDecodeError("state value exceeds the nesting-depth limit")
    budget.consume_item()
    if type(field) is ScalarField:
        if field.kind is ScalarKind.UTF8:
            assert field.effective_maximum_length is not None
            return reader.text(context=path, maximum_length=field.effective_maximum_length)
        if field.kind is ScalarKind.BYTES:
            assert field.effective_maximum_length is not None
            length = reader.head(2, context=path)
            if length > field.effective_maximum_length:
                raise StateCodecDecodeError(f"{path} exceeds its byte limit")
            return reader.read(length, context=path)
        if field.kind is ScalarKind.UINT:
            return reader.head(0, context=path)
        if field.kind is ScalarKind.SINT:
            initial = reader.data[reader.position] if reader.position < len(reader.data) else None
            if initial is None:
                raise StateCodecDecodeError(f"truncated {path}")
            major = initial >> 5
            if major == 0:
                return reader.head(0, context=path)
            if major == 1:
                return -1 - reader.head(1, context=path)
            raise StateCodecDecodeError(f"{path} is not a signed CBOR integer")
        if field.kind is ScalarKind.BOOL:
            value = reader.byte(context=path)
            if value == 0xF4:
                return False
            if value == 0xF5:
                return True
            raise StateCodecDecodeError(f"{path} is not an exact CBOR boolean")
        if field.kind is ScalarKind.FLOAT64:
            initial = reader.byte(context=path)
            if initial != 0xFB:
                raise StateCodecDecodeError(f"{path} must use the selected binary64 encoding")
            value = struct.unpack(">d", reader.read(8, context=path))[0]
            try:
                _require_finite_binary64(path, value, nonnegative=field.nonnegative)
            except (TypeError, ValueError) as error:
                raise StateCodecDecodeError(str(error)) from error
            return value
        raise AssertionError("unreachable scalar kind")

    if type(field) is EnumField:
        value = reader.text(context=path, maximum_length=128)
        if value not in field.tokens:
            raise StateCodecDecodeError(f"{path} contains an unknown enum token")
        return value

    if type(field) is ArrayField:
        count = reader.head(4, context=path)
        if count > field.maximum_length:
            raise StateCodecDecodeError(f"{path} exceeds its array length limit")
        if field.length is not None and count != field.length:
            raise StateCodecDecodeError(f"{path} must contain exactly {field.length} items")
        return tuple(
            _decode_field(
                reader,
                budget,
                field.item,
                path=f"{path}[{index}]",
                depth=depth + 1,
            )
            for index in range(count)
        )

    if type(field) is RecordField:
        count = reader.head(4, context=path)
        if count != len(field.fields):
            raise StateCodecDecodeError(f"{path} differs from its fixed record field count")
        return tuple(
            _decode_field(
                reader,
                budget,
                child,
                path=f"{path}.{child.name}",
                depth=depth + 1,
            )
            for child in field.fields
        )

    if type(field) is VariantField:
        count = reader.head(4, context=path)
        if count != 2:
            raise StateCodecDecodeError(f"{path} variant must contain exactly two items")
        tag = reader.head(0, context=f"{path} variant tag")
        try:
            arm = field.resolve(tag)
        except StateCodecError as error:
            raise StateCodecDecodeError(str(error)) from error
        arm_count = reader.head(4, context=f"{path}.{arm.name}")
        if arm_count != len(arm.fields):
            raise StateCodecDecodeError(f"{path} differs from its selected variant-arm field count")
        budget.consume_item(count=2)
        arm_values = tuple(
            _decode_field(
                reader,
                budget,
                child,
                path=f"{path}.{arm.name}.{child.name}",
                depth=depth + 1,
            )
            for child in arm.fields
        )
        return tag, arm_values

    raise TypeError("field must be an exact codec field specification")


#: Owner-approved B7 memo (2026-08-31): decode results keyed on the EXACT
#: payload bytes plus the registry instance.  Identical bytes under the
#: same registry were authenticated once by the full gate (parse, schema
#: resolve, canonical re-encode equality, identity digest); any tampered
#: byte is a different key and a full re-decode.  Refusals are never
#: cached.  The value holds a strong reference to the registry, so the
#: id() in the key can never be reused while its entry lives, and the
#: identity check on hit is belt-and-braces.
_DECODE_CACHE: dict[tuple[int, bytes], tuple[StateSchemaRegistry, DecodedState]] = {}
_DECODE_CACHE_ENTRY_LIMIT = 200_000
_decode_cache_hits = 0


def decode_state(data: bytes, registry: StateSchemaRegistry) -> DecodedState:
    """Decode and re-encode one canonical state under an explicit registry.

    Successful decodes are memoized on the exact payload bytes (B7); the
    cached object is the same frozen ``DecodedState`` the uncached gate
    constructed for those exact bytes under that exact registry.
    """

    global _decode_cache_hits
    key: tuple[int, bytes] | None = None
    if type(registry) is StateSchemaRegistry and type(data) is bytes:
        key = (id(registry), data)
        hit = _DECODE_CACHE.get(key)
        if hit is not None and hit[0] is registry:
            _decode_cache_hits += 1
            return hit[1]
    result = _decode_state_uncached(data, registry)
    if key is not None and len(_DECODE_CACHE) < _DECODE_CACHE_ENTRY_LIMIT:
        _DECODE_CACHE[key] = (registry, result)
    return result


def _decode_state_uncached(data: bytes, registry: StateSchemaRegistry) -> DecodedState:
    if type(registry) is not StateSchemaRegistry:
        raise TypeError("registry must be an exact StateSchemaRegistry")
    reader = _Reader(data)
    top_level_count = reader.head(4, context="state envelope")
    if top_level_count != 7:
        raise StateCodecDecodeError("state envelope must contain exactly seven items")
    magic = reader.text(context="profile magic", maximum_length=128)
    if magic != PROFILE_MAGIC:
        raise StateCodecDecodeError("state profile magic does not match")
    profile_version = reader.head(0, context="profile_version")
    if profile_version != PROFILE_VERSION:
        raise StateCodecDecodeError("state profile version does not match")
    schema_id = reader.text(context="schema_id", maximum_length=128)
    try:
        _require_schema_identifier("schema_id", schema_id)
    except (TypeError, ValueError) as error:
        raise StateCodecDecodeError(str(error)) from error
    schema_version = reader.head(0, context="schema_version")
    variant_tag = reader.head(0, context="variant_tag")
    schema = registry.resolve(schema_id, schema_version, variant_tag)
    if len(data) > schema.maximum_encoded_bytes:
        raise StateCodecDecodeError("encoded state exceeds the schema byte limit")
    schema_definition_length = reader.head(2, context="schema definition digest")
    if schema_definition_length != 32:
        raise StateCodecDecodeError("schema definition digest must contain exactly 32 bytes")
    schema_definition = reader.read(32, context="schema definition digest")
    expected_schema_definition = bytes.fromhex(schema.definition_digest.removeprefix("sha256:"))
    if schema_definition != expected_schema_definition:
        raise StateCodecDecodeError("schema definition digest does not match the registry")
    field_count = reader.head(4, context="field values")
    if field_count != len(schema.fields):
        raise StateCodecDecodeError("field value count differs from the fixed schema")
    budget = _DecodeBudget(maximum_items=schema.maximum_items)
    values = tuple(
        _decode_field(reader, budget, field, path=field.name, depth=1) for field in schema.fields
    )
    if reader.position != len(data):
        raise StateCodecDecodeError("trailing bytes follow the state envelope")
    reencoded = encode_state(schema, values)
    if reencoded != data:
        raise StateCodecDecodeError("state bytes are not canonical for the selected profile")
    identity_digest = "sha256:" + hashlib.sha256(data).hexdigest()
    return DecodedState(
        schema=schema,
        values=values,
        canonical_bytes=data,
        identity_digest=identity_digest,
        _issue_token=_ISSUE_TOKEN,
    )


def canonical_state_identity(data: bytes, registry: StateSchemaRegistry) -> str:
    """Return SHA-256 only after validating the complete canonical state."""

    return decode_state(data, registry).identity_digest


__all__ = [
    "ArrayField",
    "CodecValue",
    "DecodedState",
    "EnumField",
    "FieldSpec",
    "MAXIMUM_ARRAY_ITEMS",
    "MAXIMUM_CANONICAL_STATE_BYTES",
    "MAXIMUM_NESTING_DEPTH",
    "MAXIMUM_RECORD_FIELDS",
    "MAXIMUM_SCALAR_BYTES",
    "MAXIMUM_TOTAL_ITEMS",
    "PROFILE_MAGIC",
    "PROFILE_VERSION",
    "RecordField",
    "ScalarField",
    "ScalarKind",
    "StateCodecDecodeError",
    "StateCodecError",
    "StateSchema",
    "StateSchemaRegistry",
    "VariantArm",
    "VariantField",
    "canonical_state_identity",
    "decode_state",
    "encode_state",
]
