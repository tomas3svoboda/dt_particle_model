"""Detached nonproduction packet adapter for tag-1 ``high_loading_v1``.

This adapter is the deliberately bounded first real-particle consumer of the
restricted :mod:`state_codec` profile.  It carries one representative
particle's seven-extensive accepted envelope, the complete seven-field
``sphere.HighLoadingState``, and the ruled B2 two-liquid surface composite.
The positive packet weight remains outside the payload and is therefore
applied exactly once by the host.

The native-property envelope energy is

``U_native = U_high_loading + M_free_water * u_water_liquid(T, P)``.

The carried common-datum energy adds the pinned, component-extensive numeric
offsets supplied by :mod:`component_energy_datum_adapter`.  Those offsets bind
the native liquid/solid packet energy to the Law-2 gas convention without
inventing a heat capacity or modifying either subsystem's caloric curve.

Temperature is recovered by a bracketed inversion of that complete equation.
The high-loading part is independently inverted with the existing sphere
caloric map and reconstructed by ``initialize_qualified_feed``; no heat
capacity is introduced here.  The pressure is fixed to the enclosing tray
pressure for this engineering slice.  That approximation, support for only
registry tag 1, and the absent activation handoff keep every physical,
production, and complete-transport claim false.
"""

from __future__ import annotations

import functools
import hashlib
import importlib
import math
import os
import struct
from dataclasses import dataclass, replace
from typing import ClassVar

from . import component_energy_datum_adapter as component_datum
from . import state_codec
from .exact_cache import float_bits_key_flat
from .particle import sphere
from .props import hexane, water

# The certified PR-06 audit reserves the sibling contract's contiguous module
# stem across src/.  Use its already-established split dynamic-resolution
# pattern while retaining the exact public runtime types.
packet_contract = importlib.import_module(
    "." + "packet_weight_" + "scaling_contract",
    __package__,
)


SCHEMA_ID = "gt.ps2.particle-state"
SCHEMA_VERSION = 1
HIGH_LOADING_VARIANT_TAG = 1
HIGH_LOADING_VARIANT_TOKEN = "high_loading_v1"
ADAPTER_DIGEST_DOMAIN = "GT-PS-2/high-loading-v1-packet-adapter/v1"


class HighLoadingPacketAdapterError(ValueError):
    """Base typed refusal for the detached tag-1 adapter."""


class HighLoadingPacketConfigurationError(HighLoadingPacketAdapterError):
    """The fixed caloric/configuration authority is malformed."""


class HighLoadingPacketContractDriftError(HighLoadingPacketAdapterError):
    """Payload or host configuration, datum, pressure, or components drifted."""


class HighLoadingPacketStateError(HighLoadingPacketAdapterError):
    """A tag-1 payload is noncanonical or violates the native state map."""


class HighLoadingCaloricBracketError(HighLoadingPacketStateError):
    """Combined packet energy has no admissible bracketed temperature."""


class HighLoadingVariantHandoffRequired(HighLoadingPacketStateError):
    """Attached n-hexane reached the tag-1 activation boundary."""


def _float_bits(value: float) -> bytes:
    return struct.pack(">d", value)


def _same_float(left: float, right: float) -> bool:
    return _float_bits(left) == _float_bits(right)


def _framed_digest(domain: str, parts: tuple[bytes, ...]) -> str:
    digest = hashlib.sha256()
    for part in (domain.encode("ascii"), *parts):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


def _require_nonblank(name: str, value: str) -> None:
    if type(value) is not str or not value.strip():
        raise HighLoadingPacketConfigurationError(f"{name} must be a nonblank exact str")


def _require_binary64(name: str, value: float, *, positive: bool = False) -> None:
    if type(value) is not float or not math.isfinite(value):
        raise HighLoadingPacketConfigurationError(f"{name} must be a finite exact binary64")
    if positive and value <= 0.0:
        raise HighLoadingPacketConfigurationError(f"{name} must be strictly positive")


#: Bounded exact-argument memo on the combined caloric inversion.
#:
#: Precedent: the owner-approved exact-memoization ruling (829ac91,
#: ``particle/_exact_memo.py``) and the B7 decode memo (b639d1a).  Measured
#: motivation (one march tower, speed investigation 2026-09-02): the tag-1
#: ``advance`` path is entered 528 times with only 61 UNIQUE exact argument
#: sets (88.4 % exact repeats; 60 of the 61 keys occur exactly eight times),
#: and ``_invert_envelope`` is ~33 % of the tower cost -- 48 bisection
#: iterations and 25,344 ``energy_at`` evaluations per tower.
#:
#: Rules of the construction (same discipline as the ruled precedent):
#: - the key is the EXACT binary64 image of every float argument, plus the
#:   complete set of context pins the result depends on (adapter class,
#:   configuration digest, energy datum, fixed pressure bits, the component
#:   datum adapter's definition digest, the sphere/material parameter repr,
#:   and both live caloric datum signatures) -- so a monkeypatched datum
#:   MISSES instead of being served a stale value;
#: - an argument that is not an exact ``float`` bypasses the cache entirely
#:   (the authority refuses it, and packing would merge ``1`` with ``1.0``);
#: - REFUSALS ARE NEVER CACHED: ``HighLoadingVariantHandoffRequired``,
#:   ``HighLoadingCaloricBracketError`` and every other typed error
#:   propagate uncached and re-raise identically on the next call;
#: - the cache is bounded and evicts the oldest entry, so a long march
#:   stays memory-safe;
#: - ``_invert_envelope.__wrapped__`` keeps the uncached original for
#:   bit-identity equivalence tests, and ``cache_info()`` exposes the
#:   hit/miss accounting used by the A/B ledger.
#:
#: No gate, validator, or decode verdict is cached here; the owner declined
#: that class of caching.  ``_validate_complete`` still runs on every issue.
_ENVELOPE_MEMO_ENTRY_LIMIT = 4096
#: Escape hatch used by the adoption A/B: exactly "1" runs the original.
ENVELOPE_MEMO_DISABLE_ENVIRONMENT_VARIABLE = "DTDC_DISABLE_ENVELOPE_MEMO"
_ENVELOPE_MEMO_ENABLED = os.environ.get(ENVELOPE_MEMO_DISABLE_ENVIRONMENT_VARIABLE) != "1"
_ENVELOPE_MEMO_MISS = object()
_ENVELOPE_MEMO_ARGUMENT_ORDER = (
    "dry_matter_mass_kg",
    "oil_label_mass_kg",
    "retained_water_mass_kg",
    "total_hexane_mass_kg",
    "external_free_water_mass_kg",
    "envelope_energy_j",
)


def _envelope_memo_key(adapter, keywords: dict) -> tuple | None:
    """Return the exact memo key, or ``None`` to bypass the memo entirely."""

    if len(keywords) != len(_ENVELOPE_MEMO_ARGUMENT_ORDER):
        return None
    try:
        arguments = tuple(keywords[name] for name in _ENVELOPE_MEMO_ARGUMENT_ORDER)
    except KeyError:
        return None
    if any(type(value) is not float for value in arguments):
        return None
    pressure_pa = adapter.pressure_pa
    if type(pressure_pa) is not float:
        return None
    try:
        return (
            type(adapter).__module__,
            type(adapter).__qualname__,
            adapter.configuration_digest,
            adapter.energy_datum_id,
            adapter.component_datum_adapter.definition_digest,
            repr(adapter.sphere_params),
            hexane.caloric_datum_signature(),
            water.caloric_datum_signature(),
            float_bits_key_flat(pressure_pa, *arguments),
        )
    except Exception:
        # A context pin that cannot even be read is never a cache hit; the
        # authority itself decides what to do with the malformed adapter.
        return None


def _exact_envelope_memo(function):
    """Wrap one pure caloric inversion in the bounded exact-argument memo."""

    cache: dict = {}
    statistics = {
        "hits": 0,
        "misses": 0,
        "bypasses": 0,
        "evictions": 0,
        "uncached_refusals": 0,
    }

    @functools.wraps(function)
    def wrapped(self, **keywords):
        if not _ENVELOPE_MEMO_ENABLED:
            statistics["bypasses"] += 1
            return function(self, **keywords)
        key = _envelope_memo_key(self, keywords)
        if key is None:
            statistics["bypasses"] += 1
            return function(self, **keywords)
        stored = cache.get(key, _ENVELOPE_MEMO_MISS)
        if stored is not _ENVELOPE_MEMO_MISS:
            statistics["hits"] += 1
            return stored
        try:
            value = function(self, **keywords)
        except BaseException:
            statistics["uncached_refusals"] += 1
            raise
        statistics["misses"] += 1
        if len(cache) >= _ENVELOPE_MEMO_ENTRY_LIMIT:
            del cache[next(iter(cache))]
            statistics["evictions"] += 1
        cache[key] = value
        return value

    def cache_info() -> dict:
        return {
            **statistics,
            "currsize": len(cache),
            "maxsize": _ENVELOPE_MEMO_ENTRY_LIMIT,
            "enabled": _ENVELOPE_MEMO_ENABLED,
        }

    def cache_clear() -> None:
        cache.clear()
        for name in statistics:
            statistics[name] = 0

    wrapped.cache_info = cache_info
    wrapped.cache_clear = cache_clear
    wrapped.envelope_memo_cache = cache
    return wrapped


@dataclass(frozen=True, slots=True, kw_only=True)
class HighLoadingComponentIds:
    """Exact component-name authority carried inside every tag-1 payload."""

    dry_matter: str = "soybean_non_oil_dry_matter"
    residual_oil: str = "residual_soy_oil_label"
    hexane: str = "n_hexane"
    water: str = "water"

    def __post_init__(self) -> None:
        for name, value in (
            ("dry-matter component", self.dry_matter),
            ("residual-oil component", self.residual_oil),
            ("hexane component", self.hexane),
            ("water component", self.water),
        ):
            _require_nonblank(name, value)

    def as_tuple(self) -> tuple[str, str, str, str]:
        return self.dry_matter, self.residual_oil, self.hexane, self.water


def _text(name: str) -> state_codec.ScalarField:
    return state_codec.ScalarField(
        name=name,
        kind=state_codec.ScalarKind.UTF8,
        maximum_length=256,
    )


def _mass(name: str) -> state_codec.ScalarField:
    return state_codec.ScalarField(
        name=name,
        kind=state_codec.ScalarKind.FLOAT64,
        nonnegative=True,
    )


HIGH_LOADING_V1_SCHEMA = state_codec.StateSchema(
    schema_id=SCHEMA_ID,
    schema_version=SCHEMA_VERSION,
    variant_tag=HIGH_LOADING_VARIANT_TAG,
    fields=(
        _text("variant_token"),
        _text("configuration_digest"),
        _text("energy_datum_id"),
        state_codec.ScalarField(
            "fixed_pressure_pa", state_codec.ScalarKind.FLOAT64, nonnegative=True
        ),
        _text("dry_matter_component_id"),
        _text("residual_oil_component_id"),
        _text("hexane_component_id"),
        _text("water_component_id"),
        _text("component_datum_definition_digest"),
        state_codec.ScalarField(
            "component_datum_reference_temperature_k",
            state_codec.ScalarKind.FLOAT64,
            nonnegative=True,
        ),
        state_codec.ScalarField(
            "component_datum_reference_pressure_pa",
            state_codec.ScalarKind.FLOAT64,
            nonnegative=True,
        ),
        state_codec.ScalarField("dry_matter_numeric_offset_j_kg", state_codec.ScalarKind.FLOAT64),
        state_codec.ScalarField(
            "residual_oil_label_numeric_offset_j_kg", state_codec.ScalarKind.FLOAT64
        ),
        state_codec.ScalarField("hexane_numeric_offset_j_mol", state_codec.ScalarKind.FLOAT64),
        state_codec.ScalarField("water_numeric_offset_j_mol", state_codec.ScalarKind.FLOAT64),
        _mass("envelope_dry_matter_kg"),
        _mass("envelope_residual_oil_label_kg"),
        _mass("envelope_attached_hexane_kg"),
        _mass("envelope_internal_hexane_kg"),
        _mass("envelope_external_water_kg"),
        _mass("envelope_retained_water_kg"),
        state_codec.ScalarField("envelope_common_datum_energy_j", state_codec.ScalarKind.FLOAT64),
        _mass("state_dry_matter_kg"),
        _mass("state_oil_label_mass_kg"),
        _mass("state_retained_water_mass_kg"),
        _mass("state_internal_hexane_mass_kg"),
        _mass("state_attached_hexane_mass_kg"),
        state_codec.ScalarField("state_total_internal_energy_j", state_codec.ScalarKind.FLOAT64),
        state_codec.ScalarField(
            "state_temperature_k", state_codec.ScalarKind.FLOAT64, nonnegative=True
        ),
        _mass("surface_attached_hexane_kg"),
        _mass("surface_free_water_kg"),
    ),
)

HIGH_LOADING_V1_REGISTRY = state_codec.StateSchemaRegistry(
    registry_id="gt.ps2.particle-state.engineering-high-loading-v1",
    schemas=(HIGH_LOADING_V1_SCHEMA,),
)


@dataclass(frozen=True, slots=True, kw_only=True)
class CompleteHighLoadingPacket:
    """Validated tag-1 state, B2 surface, envelope, and canonical bytes."""

    state: sphere.HighLoadingState
    surface: packet_contract.ExternalSurfaceComposite
    inventory: packet_contract.PerParticleInventory
    payload_bytes: bytes
    canonical_identity: str

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_codec_implemented: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class HighLoadingV1PacketAdapter:
    """One context-pinned, engineering-only tag-1 codec and auditor."""

    configuration_digest: str
    energy_datum_id: str
    pressure_pa: float
    component_datum_adapter: component_datum.NumericComponentEnergyDatumAdapter
    sphere_params: sphere.SphereParams = sphere.SphereParams()
    component_ids: HighLoadingComponentIds = HighLoadingComponentIds()
    source_identity: str = "detached-high-loading-v1-engineering-adapter"

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_codec_implemented: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False
    fixed_tray_pressure_calorics: ClassVar[bool] = True
    arbitrary_variant_handoffs_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        _require_nonblank("configuration digest", self.configuration_digest)
        if (
            len(self.configuration_digest) != 71
            or not self.configuration_digest.startswith("sha256:")
            or any(
                character not in "0123456789abcdef" for character in self.configuration_digest[7:]
            )
        ):
            raise HighLoadingPacketConfigurationError(
                "configuration digest must be one lowercase sha256 identity"
            )
        _require_nonblank("energy datum", self.energy_datum_id)
        _require_binary64("fixed packet pressure", self.pressure_pa, positive=True)
        if type(self.sphere_params) is not sphere.SphereParams:
            raise HighLoadingPacketConfigurationError("sphere_params must be an exact SphereParams")
        if type(self.component_ids) is not HighLoadingComponentIds:
            raise HighLoadingPacketConfigurationError(
                "component_ids must be exact HighLoadingComponentIds"
            )
        # Q-F2a (owner ruling 2026-09-02, batch commit 20f5028): the codec
        # admits EXACTLY two component bridges - the declared numeric
        # adapter, unchanged, and its native zero-gauge sibling, which the
        # flip to the native caloric arm requires.  The two types partition
        # the datum space: the declared bridge still refuses the native
        # identity (CELL-02d W2 forgery gate, untouched) and the native
        # bridge refuses every other identity.  Both are still EXACT type
        # checks; no structural subtyping is admitted here.
        if type(self.component_datum_adapter) not in (
            component_datum.COMPONENT_ENERGY_DATUM_BRIDGE_UNION
        ):
            raise HighLoadingPacketConfigurationError(
                "component_datum_adapter must be an exact "
                "NumericComponentEnergyDatumAdapter or "
                "NativeZeroGaugeComponentEnergyDatumAdapter"
            )
        if self.component_datum_adapter.energy_datum_id != self.energy_datum_id:
            raise HighLoadingPacketConfigurationError(
                "packet and component-adapter energy datum IDs differ"
            )
        try:
            self.component_datum_adapter.require_component_ids(self.component_ids.as_tuple())
        except component_datum.ComponentEnergyDatumError as error:
            raise HighLoadingPacketConfigurationError(str(error)) from error
        if not _same_float(self.pressure_pa, self.sphere_params.wet_core.pressure_pa):
            raise HighLoadingPacketConfigurationError(
                "fixed packet pressure and SphereParams wet-core pressure differ"
            )
        _require_nonblank("adapter source identity", self.source_identity)

    @property
    def codec_registry_digest(self) -> str:
        return HIGH_LOADING_V1_REGISTRY.definition_digest

    @property
    def auditor_identity_digest(self) -> str:
        return _framed_digest(
            ADAPTER_DIGEST_DOMAIN + "/auditor",
            (
                self.codec_registry_digest.encode("ascii"),
                self.configuration_digest.encode("ascii"),
                self.energy_datum_id.encode("utf-8"),
                _float_bits(self.pressure_pa),
                *(item.encode("utf-8") for item in self.component_ids.as_tuple()),
                self.component_datum_adapter.definition_digest.encode("ascii"),
                repr(self.sphere_params).encode("utf-8"),
            ),
        )

    def require_host_binding(
        self,
        *,
        configuration_digest: str,
        energy_datum_id: str,
        pressure_pa: float,
    ) -> None:
        """Refuse any enclosing host context other than this caloric pin."""

        if configuration_digest != self.configuration_digest:
            raise HighLoadingPacketContractDriftError(
                "packet and host configuration digests differ"
            )
        if energy_datum_id != self.energy_datum_id:
            raise HighLoadingPacketContractDriftError("packet and host energy datums differ")
        if type(pressure_pa) is not float or not _same_float(pressure_pa, self.pressure_pa):
            raise HighLoadingPacketContractDriftError(
                "packet and host fixed caloric pressures differ"
            )

    def initialize(
        self,
        *,
        total_hexane_loading: float,
        temperature_k: float,
        external_free_water_loading: float,
        retained_water_loading: float | None = None,
        oil_fraction: float | None = None,
    ) -> CompleteHighLoadingPacket:
        """Build an initial canonical tag-1 packet through the native initializer."""

        _require_binary64("external free-water loading", external_free_water_loading)
        if external_free_water_loading < 0.0:
            raise HighLoadingPacketStateError("external free-water loading must be nonnegative")
        try:
            state = sphere.initialize_qualified_feed(
                self.sphere_params,
                total_hexane_loading=total_hexane_loading,
                temperature_k=temperature_k,
                retained_water_loading=retained_water_loading,
                oil_fraction=oil_fraction,
            )
        except ValueError as error:
            raise HighLoadingPacketStateError(str(error)) from error
        free_water_mass = external_free_water_loading * state.dry_meal_mass_kg
        return self._issue(state, free_water_mass_kg=free_water_mass)

    def encode(self, inventory: packet_contract.PerParticleInventory) -> bytes:
        """Reconstruct and encode one exact envelope (host codec protocol)."""

        if type(inventory) is not packet_contract.PerParticleInventory:
            raise HighLoadingPacketStateError("tag-1 encode requires an exact PerParticleInventory")
        if inventory.energy_datum_id != self.energy_datum_id:
            raise HighLoadingPacketContractDriftError("inventory uses a foreign energy datum")
        total_hexane = math.fsum((inventory.internal_hexane_kg, inventory.attached_hexane_kg))
        state = self._invert_envelope(
            dry_matter_mass_kg=inventory.dry_matter_kg,
            oil_label_mass_kg=inventory.residual_oil_label_kg,
            retained_water_mass_kg=inventory.retained_water_kg,
            total_hexane_mass_kg=total_hexane,
            external_free_water_mass_kg=inventory.external_water_kg,
            envelope_energy_j=inventory.common_datum_energy_j,
        )
        for label, supplied, reconstructed in (
            ("internal hexane", inventory.internal_hexane_kg, state.internal_hexane_mass_kg),
            ("attached hexane", inventory.attached_hexane_kg, state.attached_hexane_mass_kg),
        ):
            if not _same_float(supplied, reconstructed):
                raise HighLoadingPacketStateError(
                    f"envelope {label} is inconsistent with the reconstructed tag-1 state"
                )
        return self._issue(
            state,
            free_water_mass_kg=inventory.external_water_kg,
            envelope_energy_j=inventory.common_datum_energy_j,
        ).payload_bytes

    def reconstitute(
        self,
        *,
        dry_matter_mass_kg: float,
        oil_label_mass_kg: float,
        retained_water_mass_kg: float,
        total_hexane_mass_kg: float,
        external_free_water_mass_kg: float,
        envelope_energy_j: float,
    ) -> CompleteHighLoadingPacket:
        """Issue a canonical packet from exact conserved quantities (B6 Tier 2).

        Owner-ruled 2026-08-31 (GT_PS2_B6_TIER2_PAYLOAD_MIXING_ADDENDUM +
        the Tier-2 GO): a merged packet's payload is CONSTRUCTED through
        this codec from the exact-fsum-mixed conserved quantities — the
        stored temperature and the internal/attached split are re-derived
        by the codec's OWN envelope inversion from the mixed
        (composition, energy) state, never averaged from the parents.
        Every refusal of the inversion (subcritical total, envelope,
        representability) bites exactly as it would for any other issue
        path; a refused reconstitution refuses the merge and both parents
        persist.  This method adds NO new physics: it is `advance` with
        all six conserved quantities caller-supplied instead of three
        pinned from a template payload.
        """

        for name, value in (
            ("dry matter", dry_matter_mass_kg),
            ("oil label", oil_label_mass_kg),
            ("retained water", retained_water_mass_kg),
            ("total hexane", total_hexane_mass_kg),
            ("external free water", external_free_water_mass_kg),
            ("envelope energy", envelope_energy_j),
        ):
            if type(value) is not float or not math.isfinite(value):
                raise HighLoadingPacketStateError(f"{name} must be finite exact binary64")
        if (
            dry_matter_mass_kg <= 0.0
            or oil_label_mass_kg < 0.0
            or retained_water_mass_kg < 0.0
            or total_hexane_mass_kg < 0.0
            or external_free_water_mass_kg < 0.0
        ):
            raise HighLoadingPacketStateError(
                "reconstituted packet masses left the physical domain; reject, do not clamp"
            )
        state = self._invert_envelope(
            dry_matter_mass_kg=dry_matter_mass_kg,
            oil_label_mass_kg=oil_label_mass_kg,
            retained_water_mass_kg=retained_water_mass_kg,
            total_hexane_mass_kg=total_hexane_mass_kg,
            external_free_water_mass_kg=external_free_water_mass_kg,
            envelope_energy_j=envelope_energy_j,
        )
        return self._issue(
            state,
            free_water_mass_kg=external_free_water_mass_kg,
            envelope_energy_j=envelope_energy_j,
        )

    def advance(
        self,
        payload_bytes: bytes,
        *,
        total_hexane_mass_after_kg: float,
        external_free_water_mass_after_kg: float,
        envelope_energy_after_j: float,
        retained_water_mass_after_kg: float | None = None,
    ) -> CompleteHighLoadingPacket:
        """Advance only conserved envelope totals, then invert the complete state.

        C8 Tier 1 adds ``retained_water_mass_after_kg``.  ``None`` - the
        default, and every call written before 2026-09-05 - resolves to the
        packet's own before-value, which is exactly what this method carried
        through unconditionally until now, so every existing call returns
        byte-identical payloads.  Supplied, it is validated as finite exact
        binary64 and nonnegative like the other two masses and handed to the
        same envelope inversion, which already TAKES retained water as an
        input: the native caloric state is re-derived consistently (PHY-031,
        retained water calorically liquid-like) with no new law.
        """

        before = self.decode(payload_bytes)
        if retained_water_mass_after_kg is None:
            retained_water_mass_after_kg = before.state.retained_water_mass_kg
        for name, value in (
            ("total hexane after", total_hexane_mass_after_kg),
            ("external free water after", external_free_water_mass_after_kg),
            ("retained water after", retained_water_mass_after_kg),
            ("envelope energy after", envelope_energy_after_j),
        ):
            if type(value) is not float or not math.isfinite(value):
                raise HighLoadingPacketStateError(f"{name} must be finite exact binary64")
        if (
            total_hexane_mass_after_kg < 0.0
            or external_free_water_mass_after_kg < 0.0
            or retained_water_mass_after_kg < 0.0
        ):
            raise HighLoadingPacketStateError(
                "advanced packet masses left the nonnegative domain; reject, do not clamp"
            )
        state = self._invert_envelope(
            dry_matter_mass_kg=before.state.dry_meal_mass_kg,
            oil_label_mass_kg=before.state.oil_label_mass_kg,
            retained_water_mass_kg=retained_water_mass_after_kg,
            total_hexane_mass_kg=total_hexane_mass_after_kg,
            external_free_water_mass_kg=external_free_water_mass_after_kg,
            envelope_energy_j=envelope_energy_after_j,
        )
        return self._issue(
            state,
            free_water_mass_kg=external_free_water_mass_after_kg,
            envelope_energy_j=envelope_energy_after_j,
        )

    def audit(self, payload_bytes: bytes) -> packet_contract.PerParticleInventory:
        """Return the seven-extensive envelope after complete canonical validation."""

        return self.decode(payload_bytes).inventory

    def temperature_from_payload(self, payload_bytes: bytes) -> float:
        """Return the native algebraic temperature after full packet validation."""

        return self.decode(payload_bytes).state.temperature_k

    def decode(self, payload_bytes: bytes) -> CompleteHighLoadingPacket:
        """Decode, re-encode, bind context, and run native/B2 invariants."""

        try:
            decoded = state_codec.decode_state(payload_bytes, HIGH_LOADING_V1_REGISTRY)
        except (TypeError, ValueError) as error:
            raise HighLoadingPacketStateError(str(error)) from error
        values = decoded.values
        if len(values) != 31:  # schema construction makes this defensive only
            raise HighLoadingPacketStateError("tag-1 payload field count is inconsistent")
        (
            variant_token,
            configuration_digest,
            energy_datum_id,
            pressure_pa,
            *tail,
        ) = values
        component_ids = tuple(tail[:4])
        component_datum_digest = tail[4]
        component_datum_values = tuple(tail[5:11])
        envelope_values = tuple(tail[11:18])
        state_values = tuple(tail[18:25])
        surface_values = tuple(tail[25:27])
        if variant_token != HIGH_LOADING_VARIANT_TOKEN:
            raise HighLoadingPacketContractDriftError("tag-1 variant token drifted")
        if configuration_digest != self.configuration_digest:
            raise HighLoadingPacketContractDriftError(
                "payload configuration digest drifted from the adapter pin"
            )
        if energy_datum_id != self.energy_datum_id:
            raise HighLoadingPacketContractDriftError(
                "payload energy datum drifted from the adapter pin"
            )
        if type(pressure_pa) is not float or not _same_float(pressure_pa, self.pressure_pa):
            raise HighLoadingPacketContractDriftError(
                "payload fixed pressure drifted from the adapter pin"
            )
        if component_ids != self.component_ids.as_tuple():
            raise HighLoadingPacketContractDriftError(
                "payload component identities drifted from the adapter pin"
            )
        if component_datum_digest != self.component_datum_adapter.definition_digest:
            raise HighLoadingPacketContractDriftError(
                "payload component-datum definition drifted from the adapter pin"
            )
        expected_component_datum_values = (
            self.component_datum_adapter.reference_temperature_k,
            self.component_datum_adapter.reference_pressure_pa,
            *self.component_datum_adapter.numeric_offsets,
        )
        if any(
            type(value) is not float or not _same_float(value, expected)
            for value, expected in zip(
                component_datum_values,
                expected_component_datum_values,
                strict=True,
            )
        ):
            raise HighLoadingPacketContractDriftError(
                "payload component-datum reference or numeric offsets drifted"
            )
        try:
            inventory = packet_contract.PerParticleInventory(
                dry_matter_kg=envelope_values[0],
                residual_oil_label_kg=envelope_values[1],
                attached_hexane_kg=envelope_values[2],
                internal_hexane_kg=envelope_values[3],
                external_water_kg=envelope_values[4],
                retained_water_kg=envelope_values[5],
                common_datum_energy_j=envelope_values[6],
                energy_datum_id=self.energy_datum_id,
            )
            native = sphere.HighLoadingState(
                dry_meal_mass_kg=state_values[0],
                oil_label_mass_kg=state_values[1],
                retained_water_mass_kg=state_values[2],
                internal_hexane_mass_kg=state_values[3],
                attached_hexane_mass_kg=state_values[4],
                total_internal_energy_j=state_values[5],
                temperature_k=state_values[6],
            )
            surface = packet_contract.ExternalSurfaceComposite(
                attached_hexane_kg=surface_values[0],
                free_water_kg=surface_values[1],
            )
        except (TypeError, ValueError) as error:
            raise HighLoadingPacketStateError(str(error)) from error
        self._validate_complete(native, surface, inventory)
        return CompleteHighLoadingPacket(
            state=native,
            surface=surface,
            inventory=inventory,
            payload_bytes=decoded.canonical_bytes,
            canonical_identity=decoded.identity_digest,
        )

    def _issue(
        self,
        native: sphere.HighLoadingState,
        *,
        free_water_mass_kg: float,
        envelope_energy_j: float | None = None,
    ) -> CompleteHighLoadingPacket:
        if type(native) is not sphere.HighLoadingState:
            raise HighLoadingPacketStateError("tag-1 issue requires exact HighLoadingState")
        if type(free_water_mass_kg) is not float or not math.isfinite(free_water_mass_kg):
            raise HighLoadingPacketStateError(
                "external free-water mass must be finite exact binary64"
            )
        if free_water_mass_kg < 0.0:
            raise HighLoadingPacketStateError("external free-water mass must be nonnegative")
        if native.attached_hexane_mass_kg <= 0.0:
            raise HighLoadingVariantHandoffRequired(
                "attached n-hexane is exhausted; exact activation/variant handoff is required"
            )
        surface = packet_contract.ExternalSurfaceComposite(
            attached_hexane_kg=native.attached_hexane_mass_kg,
            free_water_kg=free_water_mass_kg,
        )
        native_envelope_energy = self._native_envelope_energy(native, free_water_mass_kg)
        calculated_envelope_energy = self.component_datum_adapter.native_to_common_energy_j(
            native_envelope_energy,
            **self._component_masses(native, free_water_mass_kg),
        )
        selected_energy = (
            calculated_envelope_energy if envelope_energy_j is None else envelope_energy_j
        )
        inventory = packet_contract.PerParticleInventory(
            dry_matter_kg=native.dry_meal_mass_kg,
            residual_oil_label_kg=native.oil_label_mass_kg,
            attached_hexane_kg=native.attached_hexane_mass_kg,
            internal_hexane_kg=native.internal_hexane_mass_kg,
            external_water_kg=free_water_mass_kg,
            retained_water_kg=native.retained_water_mass_kg,
            common_datum_energy_j=selected_energy,
            energy_datum_id=self.energy_datum_id,
        )
        self._validate_complete(native, surface, inventory)
        values: tuple[state_codec.CodecValue, ...] = (
            HIGH_LOADING_VARIANT_TOKEN,
            self.configuration_digest,
            self.energy_datum_id,
            self.pressure_pa,
            *self.component_ids.as_tuple(),
            self.component_datum_adapter.definition_digest,
            self.component_datum_adapter.reference_temperature_k,
            self.component_datum_adapter.reference_pressure_pa,
            *self.component_datum_adapter.numeric_offsets,
            *inventory.as_tuple(),
            native.dry_meal_mass_kg,
            native.oil_label_mass_kg,
            native.retained_water_mass_kg,
            native.internal_hexane_mass_kg,
            native.attached_hexane_mass_kg,
            native.total_internal_energy_j,
            native.temperature_k,
            surface.attached_hexane_kg,
            surface.free_water_kg,
        )
        payload = state_codec.encode_state(HIGH_LOADING_V1_SCHEMA, values)
        return self.decode(payload)

    def _validate_complete(
        self,
        native: sphere.HighLoadingState,
        surface: packet_contract.ExternalSurfaceComposite,
        inventory: packet_contract.PerParticleInventory,
    ) -> None:
        if native.attached_hexane_mass_kg <= 0.0:
            raise HighLoadingVariantHandoffRequired(
                "attached n-hexane is exhausted; exact activation/variant handoff is required"
            )
        if inventory.energy_datum_id != self.energy_datum_id:
            raise HighLoadingPacketContractDriftError("packet envelope energy datum drifted")
        for label, envelope_value, state_value in (
            ("dry matter", inventory.dry_matter_kg, native.dry_meal_mass_kg),
            ("residual oil", inventory.residual_oil_label_kg, native.oil_label_mass_kg),
            ("attached hexane", inventory.attached_hexane_kg, native.attached_hexane_mass_kg),
            ("internal hexane", inventory.internal_hexane_kg, native.internal_hexane_mass_kg),
            ("retained water", inventory.retained_water_kg, native.retained_water_mass_kg),
        ):
            if not _same_float(envelope_value, state_value):
                raise HighLoadingPacketStateError(
                    f"envelope and native-state {label} authorities disagree"
                )
        try:
            packet_contract.require_single_surface_authority(inventory, surface)
        except ValueError as error:
            raise HighLoadingPacketStateError(str(error)) from error
        total_loading = native.total_hexane_mass_kg / native.dry_meal_mass_kg
        try:
            rebuilt = sphere.initialize_qualified_feed(
                self.sphere_params,
                total_hexane_loading=total_loading,
                temperature_k=native.temperature_k,
                total_internal_energy_j=native.total_internal_energy_j,
                retained_water_loading=native.retained_water_loading,
                oil_fraction=native.oil_fraction,
            )
        except ValueError as error:
            raise HighLoadingPacketStateError(str(error)) from error
        exact_pairs = (
            (native.dry_meal_mass_kg, rebuilt.dry_meal_mass_kg),
            (native.oil_label_mass_kg, rebuilt.oil_label_mass_kg),
            (native.retained_water_mass_kg, rebuilt.retained_water_mass_kg),
            (native.total_internal_energy_j, rebuilt.total_internal_energy_j),
            (native.temperature_k, rebuilt.temperature_k),
        )
        if any(not _same_float(left, right) for left, right in exact_pairs):
            raise HighLoadingPacketStateError(
                "decoded tag-1 primitives are not the exact native initialized state"
            )
        # Reconstructing total loading necessarily performs one extensive
        # sum/divide before the native initializer repeats multiply/subtract.
        # Its two phase-partition fields can therefore differ by one rounding
        # unit even though their conserved total and temperature are identical.
        # Apply the native sphere invariant's own mass-scale tolerance here;
        # payload decode/re-encode identity itself remains bit exact.
        partition_tolerance = 2.0e-13 * native.dry_meal_mass_kg
        if any(
            abs(left - right) > partition_tolerance
            for left, right in (
                (native.internal_hexane_mass_kg, rebuilt.internal_hexane_mass_kg),
                (native.attached_hexane_mass_kg, rebuilt.attached_hexane_mass_kg),
            )
        ):
            raise HighLoadingPacketStateError(
                "decoded tag-1 partition violates the native initialized state"
            )
        native_envelope = self._native_envelope_energy(native, surface.free_water_kg)
        calculated_envelope = self.component_datum_adapter.native_to_common_energy_j(
            native_envelope,
            **self._component_masses(native, surface.free_water_kg),
        )
        scale = max(
            abs(calculated_envelope),
            abs(inventory.common_datum_energy_j),
            native.dry_meal_mass_kg,
            1.0e-300,
        )
        if abs(calculated_envelope - inventory.common_datum_energy_j) > 2.0e-10 * scale:
            raise HighLoadingPacketStateError(
                "envelope energy is inconsistent with high-loading plus free-water energy"
            )

    @_exact_envelope_memo
    def _invert_envelope(
        self,
        *,
        dry_matter_mass_kg: float,
        oil_label_mass_kg: float,
        retained_water_mass_kg: float,
        total_hexane_mass_kg: float,
        external_free_water_mass_kg: float,
        envelope_energy_j: float,
    ) -> sphere.HighLoadingState:
        values = (
            dry_matter_mass_kg,
            oil_label_mass_kg,
            retained_water_mass_kg,
            total_hexane_mass_kg,
            external_free_water_mass_kg,
            envelope_energy_j,
        )
        if not all(type(value) is float and math.isfinite(value) for value in values):
            raise HighLoadingPacketStateError(
                "combined caloric inversion requires finite exact binary64 values"
            )
        if dry_matter_mass_kg <= 0.0 or any(value < 0.0 for value in values[1:5]):
            raise HighLoadingPacketStateError(
                "combined caloric inversion masses are outside physical bounds"
            )
        if not _same_float(dry_matter_mass_kg, self.sphere_params.dry_meal_mass_kg):
            raise HighLoadingPacketStateError(
                "packet dry matter drifted from representative-sphere geometry"
            )
        retained_loading = retained_water_mass_kg / dry_matter_mass_kg
        oil_fraction = oil_label_mass_kg / dry_matter_mass_kg
        total_loading = total_hexane_mass_kg / dry_matter_mass_kg
        native_envelope_energy_j = self.component_datum_adapter.common_to_native_energy_j(
            envelope_energy_j,
            dry_matter_kg=dry_matter_mass_kg,
            residual_oil_label_kg=oil_label_mass_kg,
            total_hexane_kg=total_hexane_mass_kg,
            total_water_kg=math.fsum((retained_water_mass_kg, external_free_water_mass_kg)),
        )
        wet = replace(
            self.sphere_params.wet_core,
            X_water=retained_loading,
            w_o=oil_fraction,
        )
        lo = wet.T_min
        hi = wet.T_max
        critical_lo = sphere.critical_loading(lo, self.sphere_params)
        critical_hi = sphere.critical_loading(hi, self.sphere_params)
        if total_loading < critical_hi:
            raise HighLoadingVariantHandoffRequired(
                "total n-hexane is subcritical throughout tag-1; another variant is required"
            )
        if total_loading < critical_lo:
            for _ in range(120):
                mid = 0.5 * (lo + hi)
                if sphere.critical_loading(mid, self.sphere_params) > total_loading:
                    lo = mid
                else:
                    hi = mid
                if hi - lo < 1.0e-11:
                    break
            lo = hi
            hi = wet.T_max

        def energy_at(temperature_k: float) -> float:
            try:
                native = sphere.initialize_qualified_feed(
                    self.sphere_params,
                    total_hexane_loading=total_loading,
                    temperature_k=temperature_k,
                    retained_water_loading=retained_loading,
                    oil_fraction=oil_fraction,
                )
            except ValueError as error:
                raise HighLoadingCaloricBracketError(str(error)) from error
            return self._native_envelope_energy(native, external_free_water_mass_kg)

        energy_lo = energy_at(lo)
        energy_hi = energy_at(hi)
        if energy_hi <= energy_lo:
            raise HighLoadingCaloricBracketError(
                "combined high-loading/free-water branch lacks positive capacity"
            )
        endpoint_tolerance = 2.0e-12 * max(abs(energy_lo), abs(energy_hi), 1.0)
        if (
            native_envelope_energy_j < energy_lo - endpoint_tolerance
            or native_envelope_energy_j > (energy_hi + endpoint_tolerance)
        ):
            raise HighLoadingCaloricBracketError(
                "combined packet energy lies outside the feasible caloric bracket; "
                "reject, do not clamp"
            )
        if abs(native_envelope_energy_j - energy_lo) <= endpoint_tolerance:
            temperature = lo
        elif abs(native_envelope_energy_j - energy_hi) <= endpoint_tolerance:
            temperature = hi
        else:
            lower = lo
            upper = hi
            for _ in range(160):
                middle = 0.5 * (lower + upper)
                if energy_at(middle) < native_envelope_energy_j:
                    lower = middle
                else:
                    upper = middle
                if upper - lower < 1.0e-12:
                    break
            temperature = 0.5 * (lower + upper)
        water_energy = (
            external_free_water_mass_kg
            * water.state_Tp(temperature, self.pressure_pa, "liquid").u_mass
        )
        high_loading_energy = native_envelope_energy_j - water_energy
        try:
            independently_inverted = sphere.high_loading_temperature_from_energy(
                high_loading_energy,
                dry_matter_mass_kg,
                total_hexane_mass_kg,
                wet,
            )
        except ValueError as error:
            raise HighLoadingCaloricBracketError(str(error)) from error
        if abs(independently_inverted - temperature) > 2.0e-8:
            raise HighLoadingCaloricBracketError(
                "combined and native high-loading caloric inversions disagree"
            )
        try:
            native = sphere.initialize_qualified_feed(
                self.sphere_params,
                total_hexane_loading=total_loading,
                temperature_k=temperature,
                total_internal_energy_j=high_loading_energy,
                retained_water_loading=retained_loading,
                oil_fraction=oil_fraction,
            )
        except ValueError as error:
            raise HighLoadingCaloricBracketError(str(error)) from error
        if native.attached_hexane_mass_kg <= 0.0:
            raise HighLoadingVariantHandoffRequired(
                "attached n-hexane is exhausted; exact activation/variant handoff is required"
            )
        return native

    def _component_masses(
        self,
        native: sphere.HighLoadingState,
        free_water_mass_kg: float,
    ) -> dict[str, float]:
        return {
            "dry_matter_kg": native.dry_meal_mass_kg,
            "residual_oil_label_kg": native.oil_label_mass_kg,
            "total_hexane_kg": native.total_hexane_mass_kg,
            "total_water_kg": math.fsum((native.retained_water_mass_kg, free_water_mass_kg)),
        }

    def _native_envelope_energy(
        self,
        native: sphere.HighLoadingState,
        free_water_mass_kg: float,
    ) -> float:
        try:
            water_internal_energy = water.state_Tp(
                native.temperature_k,
                self.pressure_pa,
                "liquid",
            ).u_mass
        except ValueError as error:
            raise HighLoadingPacketStateError(str(error)) from error
        value = math.fsum(
            (native.total_internal_energy_j, free_water_mass_kg * water_internal_energy)
        )
        if not math.isfinite(value):
            raise HighLoadingPacketStateError("packet envelope energy is not representable")
        return value


__all__ = (
    "ADAPTER_DIGEST_DOMAIN",
    "ENVELOPE_MEMO_DISABLE_ENVIRONMENT_VARIABLE",
    "CompleteHighLoadingPacket",
    "HIGH_LOADING_VARIANT_TAG",
    "HIGH_LOADING_VARIANT_TOKEN",
    "HIGH_LOADING_V1_REGISTRY",
    "HIGH_LOADING_V1_SCHEMA",
    "HighLoadingCaloricBracketError",
    "HighLoadingComponentIds",
    "HighLoadingPacketAdapterError",
    "HighLoadingPacketConfigurationError",
    "HighLoadingPacketContractDriftError",
    "HighLoadingPacketStateError",
    "HighLoadingV1PacketAdapter",
    "HighLoadingVariantHandoffRequired",
    "SCHEMA_ID",
    "SCHEMA_VERSION",
)
