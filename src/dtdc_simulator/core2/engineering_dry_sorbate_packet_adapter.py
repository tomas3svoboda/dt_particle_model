"""TAG-3 packet codec: the dry-sorbate (mobile-liquid exhausted) variant.

Owner-ruled 2026-09-06 morning ("I rule all four as recommended, proceed.");
build specification `GT_PS2_TAG3_DRY_SORBATE_PACKET_BUILD_SPEC_2026-09-06.md`
section 2.  Faithful mirror of the tag-2 `SubcriticalPacketAdapter` contract
on the tag-3 native state (`engineering_dry_sorbate_packet_state`, seam-exact
caloric):

- SAME schema_id and wire-contract generation; **variant_tag 3** per the repo
  convention that the variant tag denotes the PHYSICAL particle variant.
- SAME context pins in every payload (configuration digest, energy datum,
  fixed pressure, component identities, component-datum definition and
  numeric offsets) - drift refuses identically.
- The envelope block keeps all seven extensives with
  ``envelope_attached_hexane_kg`` REQUIRED exactly +0.0 (canonical positive
  zero), so `PerParticleInventory`, the surface composite, and every
  conservation ledger consume tag-3 packets unchanged.  The particle's whole
  hexane inventory is booked in ``envelope_internal_hexane_kg``, which is the
  same conserved extensive tag-2 books there; only the NATIVE block renames it
  to ``state_sorbed_hexane_mass_kg``, because in this band there is no mobile
  pore liquid at all and every kilogram is on a sorption site.
- The native block therefore carries the SORBED loading: a tag-3 particle that
  re-wets (condensation, or heating past its own gate) hands back to tag-2 at
  the DUAL-adapter boundary, never mutates in place.

Deviation from the specification, stated not hidden: section 2 says
``advance`` carries the C8 Tier 1 ``retained_water_mass_after_kg`` keyword
"exactly as tag-2 carries it".  The tag-2 adapter does NOT carry that keyword
(only tag-1 does, `high_loading_packet_adapter.py:591`).  Tag-3 is the band in
which the sorbed-water arm actually binds, so the keyword is carried here in
the tag-1 form: ``None`` resolves to the packet's own before-value, which is
exactly the unconditional behaviour of the tag-2 method.

Nothing in tag-1 or tag-2 is edited by this module.

physically_qualifying: false.  plant_predictive: false.
"""

from __future__ import annotations

import hashlib
import math
import struct
from dataclasses import dataclass
from typing import ClassVar

from . import component_energy_datum_adapter as component_datum
from . import high_loading_packet_adapter as tag1
from . import state_codec
from .engineering_dry_sorbate_packet_state import (
    DrySorbateExhaustedRefusal,
    DrySorbateParticleState,
    DrySorbateReWettingRefusal,
    DrySorbateStateError,
    dry_sorbate_specific_energy,
    dry_sorbate_temperature_bracket,
    initialize_dry_sorbate_state,
)
from .envelope_inversion_memo import exact_envelope_memo
from .particle import sphere
from .props import water
from .packet_weight_scaling_contract import (
    ExternalSurfaceComposite,
    PerParticleInventory,
    require_single_surface_authority,
)

DRY_SORBATE_SCHEMA_VERSION = 1
DRY_SORBATE_VARIANT_TAG = 3
DRY_SORBATE_VARIANT_TOKEN = "engineering_dry_sorbate_v1"
DRY_SORBATE_ADAPTER_DIGEST_DOMAIN = "GT-PS-2/engineering-dry-sorbate-v1-packet-adapter/v1"

#: Keyword order of the tag-3 caloric inversion, and therefore the packing
#: order of the exact-argument memo key (owner ruling 2026-09-06, item 4:
#: "memoize the dry-sorbate caloric envelope inversion exactly the way the
#: tag-1 inversion is memoized").  See ``envelope_inversion_memo`` for the
#: construction and for what the key excludes.
DRY_SORBATE_ENVELOPE_MEMO_ARGUMENT_ORDER = (
    "dry_matter_mass_kg",
    "oil_label_mass_kg",
    "retained_water_mass_kg",
    "sorbed_hexane_mass_kg",
    "external_free_water_mass_kg",
    "envelope_energy_j",
)


class DrySorbatePacketAdapterError(ValueError):
    """Base typed refusal for the detached tag-3 adapter."""


class DrySorbatePacketConfigurationError(DrySorbatePacketAdapterError):
    """The fixed caloric/configuration authority is malformed."""


class DrySorbatePacketContractDriftError(DrySorbatePacketAdapterError):
    """Payload or host configuration, datum, pressure, or components drifted."""


class DrySorbatePacketStateError(DrySorbatePacketAdapterError):
    """A tag-3 payload is noncanonical or violates the native state map."""


class DrySorbateCaloricBracketError(DrySorbatePacketStateError):
    """Combined packet energy has no admissible bracketed temperature."""


class DrySorbateVariantHandoffRequired(DrySorbatePacketStateError):
    """The state left the tag-3 band; the tag-2 (subcritical) variant owns it."""


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


DRY_SORBATE_V1_SCHEMA = state_codec.StateSchema(
    schema_id=tag1.SCHEMA_ID,
    schema_version=DRY_SORBATE_SCHEMA_VERSION,
    variant_tag=DRY_SORBATE_VARIANT_TAG,
    fields=(
        tag1._text("variant_token"),
        tag1._text("configuration_digest"),
        tag1._text("energy_datum_id"),
        state_codec.ScalarField(
            "fixed_pressure_pa", state_codec.ScalarKind.FLOAT64, nonnegative=True
        ),
        tag1._text("dry_matter_component_id"),
        tag1._text("residual_oil_component_id"),
        tag1._text("hexane_component_id"),
        tag1._text("water_component_id"),
        tag1._text("component_datum_definition_digest"),
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
        tag1._mass("envelope_dry_matter_kg"),
        tag1._mass("envelope_residual_oil_label_kg"),
        tag1._mass("envelope_attached_hexane_kg"),
        tag1._mass("envelope_internal_hexane_kg"),
        tag1._mass("envelope_external_water_kg"),
        tag1._mass("envelope_retained_water_kg"),
        state_codec.ScalarField("envelope_common_datum_energy_j", state_codec.ScalarKind.FLOAT64),
        tag1._mass("state_dry_matter_kg"),
        tag1._mass("state_oil_label_mass_kg"),
        tag1._mass("state_retained_water_mass_kg"),
        tag1._mass("state_sorbed_hexane_mass_kg"),
        state_codec.ScalarField("state_total_internal_energy_j", state_codec.ScalarKind.FLOAT64),
        state_codec.ScalarField(
            "state_temperature_k", state_codec.ScalarKind.FLOAT64, nonnegative=True
        ),
        tag1._mass("surface_attached_hexane_kg"),
        tag1._mass("surface_free_water_kg"),
    ),
)

DRY_SORBATE_V1_REGISTRY = state_codec.StateSchemaRegistry(
    registry_id="gt.ps2.particle-state.engineering-dry-sorbate-v1",
    schemas=(DRY_SORBATE_V1_SCHEMA,),
)


@dataclass(frozen=True, slots=True, kw_only=True)
class CompleteDrySorbatePacket:
    """Validated tag-3 state, surface, envelope, and canonical bytes."""

    state: DrySorbateParticleState
    surface: ExternalSurfaceComposite
    inventory: PerParticleInventory
    payload_bytes: bytes
    canonical_identity: str

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class DrySorbatePacketAdapter:
    """One context-pinned, engineering-only tag-3 codec and auditor."""

    configuration_digest: str
    energy_datum_id: str
    pressure_pa: float
    component_datum_adapter: component_datum.NumericComponentEnergyDatumAdapter
    sphere_params: sphere.SphereParams = sphere.SphereParams()
    component_ids: tag1.HighLoadingComponentIds = tag1.HighLoadingComponentIds()
    source_identity: str = "detached-engineering-dry-sorbate-v1-adapter"

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_codec_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        # The context gates are exactly the tag-1 gates; run them through a
        # throwaway tag-1 adapter so the variants can never drift apart.
        try:
            tag1.HighLoadingV1PacketAdapter(
                configuration_digest=self.configuration_digest,
                energy_datum_id=self.energy_datum_id,
                pressure_pa=self.pressure_pa,
                component_datum_adapter=self.component_datum_adapter,
                sphere_params=self.sphere_params,
                component_ids=self.component_ids,
            )
        except tag1.HighLoadingPacketAdapterError as error:
            raise DrySorbatePacketConfigurationError(str(error)) from error
        if type(self.source_identity) is not str or not self.source_identity.strip():
            raise DrySorbatePacketConfigurationError(
                "adapter source identity must be a nonblank exact str"
            )

    @property
    def codec_registry_digest(self) -> str:
        return DRY_SORBATE_V1_REGISTRY.definition_digest

    @property
    def auditor_identity_digest(self) -> str:
        return _framed_digest(
            DRY_SORBATE_ADAPTER_DIGEST_DOMAIN + "/auditor",
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
        if configuration_digest != self.configuration_digest:
            raise DrySorbatePacketContractDriftError("packet and host configuration digests differ")
        if energy_datum_id != self.energy_datum_id:
            raise DrySorbatePacketContractDriftError("packet and host energy datums differ")
        if type(pressure_pa) is not float or not _same_float(pressure_pa, self.pressure_pa):
            raise DrySorbatePacketContractDriftError(
                "packet and host fixed caloric pressures differ"
            )

    def initialize(
        self,
        *,
        sorbed_hexane_loading: float,
        temperature_k: float,
        external_free_water_loading: float,
        retained_water_loading: float | None = None,
        oil_fraction: float | None = None,
    ) -> CompleteDrySorbatePacket:
        if (
            type(external_free_water_loading) is not float
            or not math.isfinite(external_free_water_loading)
            or external_free_water_loading < 0.0
        ):
            raise DrySorbatePacketStateError(
                "external free-water loading must be finite nonnegative binary64"
            )
        try:
            state = initialize_dry_sorbate_state(
                self.sphere_params,
                sorbed_hexane_loading=sorbed_hexane_loading,
                temperature_k=temperature_k,
                retained_water_loading=retained_water_loading,
                oil_fraction=oil_fraction,
            )
        except DrySorbateStateError as error:
            raise DrySorbatePacketStateError(str(error)) from error
        free_water_mass = external_free_water_loading * state.dry_meal_mass_kg
        return self._issue(state, free_water_mass_kg=free_water_mass)

    def reconstitute(
        self,
        *,
        dry_matter_mass_kg: float,
        oil_label_mass_kg: float,
        retained_water_mass_kg: float,
        total_hexane_mass_kg: float,
        external_free_water_mass_kg: float,
        envelope_energy_j: float,
    ) -> CompleteDrySorbatePacket:
        """Issue a canonical tag-3 packet from exact conserved quantities."""

        for name, value in (
            ("dry matter", dry_matter_mass_kg),
            ("oil label", oil_label_mass_kg),
            ("retained water", retained_water_mass_kg),
            ("total hexane", total_hexane_mass_kg),
            ("external free water", external_free_water_mass_kg),
            ("envelope energy", envelope_energy_j),
        ):
            if type(value) is not float or not math.isfinite(value):
                raise DrySorbatePacketStateError(f"{name} must be finite exact binary64")
        if (
            dry_matter_mass_kg <= 0.0
            or oil_label_mass_kg < 0.0
            or retained_water_mass_kg < 0.0
            or total_hexane_mass_kg <= 0.0
            or external_free_water_mass_kg < 0.0
        ):
            raise DrySorbatePacketStateError(
                "reconstituted packet masses left the tag-3 domain; reject, do not clamp"
            )
        state = self._invert_envelope(
            dry_matter_mass_kg=dry_matter_mass_kg,
            oil_label_mass_kg=oil_label_mass_kg,
            retained_water_mass_kg=retained_water_mass_kg,
            sorbed_hexane_mass_kg=total_hexane_mass_kg,
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
    ) -> CompleteDrySorbatePacket:
        """Advance only conserved envelope totals, then invert the state.

        ``retained_water_mass_after_kg`` is the C8 Tier 1 keyword in its tag-1
        form: ``None`` resolves to the packet's own before-value, which is what
        the tag-2 method does unconditionally, so a caller that never passes it
        sees the tag-2 behaviour exactly.
        """

        before = self.decode(payload_bytes)
        if retained_water_mass_after_kg is None:
            retained_water_mass_after_kg = before.state.retained_water_mass_kg
        return self.reconstitute(
            dry_matter_mass_kg=before.state.dry_meal_mass_kg,
            oil_label_mass_kg=before.state.oil_label_mass_kg,
            retained_water_mass_kg=retained_water_mass_after_kg,
            total_hexane_mass_kg=total_hexane_mass_after_kg,
            external_free_water_mass_kg=external_free_water_mass_after_kg,
            envelope_energy_j=envelope_energy_after_j,
        )

    def encode(self, inventory: PerParticleInventory) -> bytes:
        if type(inventory) is not PerParticleInventory:
            raise DrySorbatePacketStateError("encode requires an exact PerParticleInventory")
        if inventory.attached_hexane_kg != 0.0:
            raise DrySorbateVariantHandoffRequired("a positive attached film is tag-1 territory")
        packet = self.reconstitute(
            dry_matter_mass_kg=inventory.dry_matter_kg,
            oil_label_mass_kg=inventory.residual_oil_label_kg,
            retained_water_mass_kg=inventory.retained_water_kg,
            total_hexane_mass_kg=inventory.internal_hexane_kg,
            external_free_water_mass_kg=inventory.external_water_kg,
            envelope_energy_j=inventory.common_datum_energy_j,
        )
        return packet.payload_bytes

    def audit(self, payload_bytes: bytes) -> PerParticleInventory:
        return self.decode(payload_bytes).inventory

    def temperature_from_payload(self, payload_bytes: bytes) -> float:
        return self.decode(payload_bytes).state.temperature_k

    def decode(self, payload_bytes: bytes) -> CompleteDrySorbatePacket:
        try:
            decoded = state_codec.decode_state(payload_bytes, DRY_SORBATE_V1_REGISTRY)
        except (TypeError, ValueError) as error:
            raise DrySorbatePacketStateError(str(error)) from error
        values = decoded.values
        if len(values) != 30:  # schema construction makes this defensive only
            raise DrySorbatePacketStateError("tag-3 payload field count is inconsistent")
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
        state_values = tuple(tail[18:24])
        surface_values = tuple(tail[24:26])
        if variant_token != DRY_SORBATE_VARIANT_TOKEN:
            raise DrySorbatePacketContractDriftError("tag-3 variant token drifted")
        if configuration_digest != self.configuration_digest:
            raise DrySorbatePacketContractDriftError(
                "payload configuration digest drifted from the adapter pin"
            )
        if energy_datum_id != self.energy_datum_id:
            raise DrySorbatePacketContractDriftError(
                "payload energy datum drifted from the adapter pin"
            )
        if type(pressure_pa) is not float or not _same_float(pressure_pa, self.pressure_pa):
            raise DrySorbatePacketContractDriftError(
                "payload fixed pressure drifted from the adapter pin"
            )
        if component_ids != self.component_ids.as_tuple():
            raise DrySorbatePacketContractDriftError(
                "payload component identities drifted from the adapter pin"
            )
        if component_datum_digest != self.component_datum_adapter.definition_digest:
            raise DrySorbatePacketContractDriftError(
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
            raise DrySorbatePacketContractDriftError(
                "payload component-datum reference or numeric offsets drifted"
            )
        try:
            inventory = PerParticleInventory(
                dry_matter_kg=envelope_values[0],
                residual_oil_label_kg=envelope_values[1],
                attached_hexane_kg=envelope_values[2],
                internal_hexane_kg=envelope_values[3],
                external_water_kg=envelope_values[4],
                retained_water_kg=envelope_values[5],
                common_datum_energy_j=envelope_values[6],
                energy_datum_id=self.energy_datum_id,
            )
            native = DrySorbateParticleState(
                dry_meal_mass_kg=state_values[0],
                oil_label_mass_kg=state_values[1],
                retained_water_mass_kg=state_values[2],
                sorbed_hexane_mass_kg=state_values[3],
                total_internal_energy_j=state_values[4],
                temperature_k=state_values[5],
            )
            surface = ExternalSurfaceComposite(
                attached_hexane_kg=surface_values[0],
                free_water_kg=surface_values[1],
            )
        except (TypeError, ValueError) as error:
            raise DrySorbatePacketStateError(str(error)) from error
        self._validate_complete(native, surface, inventory)
        return CompleteDrySorbatePacket(
            state=native,
            surface=surface,
            inventory=inventory,
            payload_bytes=decoded.canonical_bytes,
            canonical_identity=decoded.identity_digest,
        )

    def _issue(
        self,
        native: DrySorbateParticleState,
        *,
        free_water_mass_kg: float,
        envelope_energy_j: float | None = None,
    ) -> CompleteDrySorbatePacket:
        if type(native) is not DrySorbateParticleState:
            raise DrySorbatePacketStateError("tag-3 issue requires exact DrySorbateParticleState")
        if (
            type(free_water_mass_kg) is not float
            or not math.isfinite(free_water_mass_kg)
            or free_water_mass_kg < 0.0
        ):
            raise DrySorbatePacketStateError(
                "external free-water mass must be finite nonnegative binary64"
            )
        surface = ExternalSurfaceComposite(
            attached_hexane_kg=0.0,
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
        inventory = PerParticleInventory(
            dry_matter_kg=native.dry_meal_mass_kg,
            residual_oil_label_kg=native.oil_label_mass_kg,
            attached_hexane_kg=0.0,
            internal_hexane_kg=native.sorbed_hexane_mass_kg,
            external_water_kg=free_water_mass_kg,
            retained_water_kg=native.retained_water_mass_kg,
            common_datum_energy_j=selected_energy,
            energy_datum_id=self.energy_datum_id,
        )
        self._validate_complete(native, surface, inventory)
        values: tuple[state_codec.CodecValue, ...] = (
            DRY_SORBATE_VARIANT_TOKEN,
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
            native.sorbed_hexane_mass_kg,
            native.total_internal_energy_j,
            native.temperature_k,
            surface.attached_hexane_kg,
            surface.free_water_kg,
        )
        payload = state_codec.encode_state(DRY_SORBATE_V1_SCHEMA, values)
        return self.decode(payload)

    def _validate_complete(
        self,
        native: DrySorbateParticleState,
        surface: ExternalSurfaceComposite,
        inventory: PerParticleInventory,
    ) -> None:
        if inventory.attached_hexane_kg != 0.0 or surface.attached_hexane_kg != 0.0:
            raise DrySorbateVariantHandoffRequired("a positive attached film is tag-1 territory")
        if inventory.energy_datum_id != self.energy_datum_id:
            raise DrySorbatePacketContractDriftError("packet envelope energy datum drifted")
        for label, envelope_value, state_value in (
            ("dry matter", inventory.dry_matter_kg, native.dry_meal_mass_kg),
            ("residual oil", inventory.residual_oil_label_kg, native.oil_label_mass_kg),
            ("sorbed hexane", inventory.internal_hexane_kg, native.sorbed_hexane_mass_kg),
            ("retained water", inventory.retained_water_kg, native.retained_water_mass_kg),
        ):
            if not _same_float(envelope_value, state_value):
                raise DrySorbatePacketStateError(
                    f"envelope and native-state {label} authorities disagree"
                )
        try:
            require_single_surface_authority(inventory, surface)
        except ValueError as error:
            raise DrySorbatePacketStateError(str(error)) from error
        try:
            rebuilt = initialize_dry_sorbate_state(
                self.sphere_params,
                sorbed_hexane_loading=native.sorbed_hexane_loading,
                temperature_k=native.temperature_k,
                total_internal_energy_j=native.total_internal_energy_j,
                retained_water_loading=native.retained_water_loading,
                oil_fraction=native.oil_fraction,
            )
        except DrySorbateReWettingRefusal as error:
            raise DrySorbateVariantHandoffRequired(str(error)) from error
        except DrySorbateStateError as error:
            raise DrySorbatePacketStateError(str(error)) from error
        exact_pairs = (
            (native.dry_meal_mass_kg, rebuilt.dry_meal_mass_kg),
            (native.oil_label_mass_kg, rebuilt.oil_label_mass_kg),
            (native.retained_water_mass_kg, rebuilt.retained_water_mass_kg),
            (native.total_internal_energy_j, rebuilt.total_internal_energy_j),
            (native.temperature_k, rebuilt.temperature_k),
        )
        if any(not _same_float(left, right) for left, right in exact_pairs):
            raise DrySorbatePacketStateError(
                "decoded tag-3 primitives are not the exact native initialized state"
            )
        # As in tag-2 there is no phase partition to reconstruct: sorbed hexane
        # IS the total, so the loading survives one divide/multiply.
        partition_tolerance = 2.0e-13 * native.dry_meal_mass_kg
        if abs(native.sorbed_hexane_mass_kg - rebuilt.sorbed_hexane_mass_kg) > (
            partition_tolerance
        ):
            raise DrySorbatePacketStateError(
                "decoded tag-3 hexane inventory violates the native initialized state"
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
            raise DrySorbatePacketStateError(
                "envelope energy is inconsistent with dry-sorbate plus free-water energy"
            )

    @exact_envelope_memo(*DRY_SORBATE_ENVELOPE_MEMO_ARGUMENT_ORDER)
    def _invert_envelope(
        self,
        *,
        dry_matter_mass_kg: float,
        oil_label_mass_kg: float,
        retained_water_mass_kg: float,
        sorbed_hexane_mass_kg: float,
        external_free_water_mass_kg: float,
        envelope_energy_j: float,
    ) -> DrySorbateParticleState:
        if not _same_float(dry_matter_mass_kg, self.sphere_params.dry_meal_mass_kg):
            raise DrySorbatePacketStateError(
                "packet dry matter drifted from representative-sphere geometry"
            )
        retained_loading = retained_water_mass_kg / dry_matter_mass_kg
        oil_fraction = oil_label_mass_kg / dry_matter_mass_kg
        sorbed_loading = sorbed_hexane_mass_kg / dry_matter_mass_kg
        native_envelope_energy_j = self.component_datum_adapter.common_to_native_energy_j(
            envelope_energy_j,
            dry_matter_kg=dry_matter_mass_kg,
            residual_oil_label_kg=oil_label_mass_kg,
            total_hexane_kg=sorbed_hexane_mass_kg,
            total_water_kg=math.fsum((retained_water_mass_kg, external_free_water_mass_kg)),
        )
        try:
            wet = sphere._material_wet_params(self.sphere_params, retained_loading, oil_fraction)
        except ValueError as error:
            raise DrySorbatePacketStateError(str(error)) from error
        try:
            lo, hi = dry_sorbate_temperature_bracket(sorbed_loading, wet)
        except DrySorbateReWettingRefusal as error:
            raise DrySorbateVariantHandoffRequired(str(error)) from error
        except DrySorbateExhaustedRefusal:
            raise
        except DrySorbateStateError as error:
            raise DrySorbateCaloricBracketError(str(error)) from error

        def energy_at(temperature_k: float) -> float:
            particle = dry_matter_mass_kg * dry_sorbate_specific_energy(
                temperature_k, sorbed_loading, wet
            )
            water_energy = (
                external_free_water_mass_kg
                * water.state_Tp(temperature_k, self.pressure_pa, "liquid").u_mass
            )
            return particle + water_energy

        energy_lo = energy_at(lo)
        energy_hi = energy_at(hi)
        if energy_hi <= energy_lo:
            raise DrySorbateCaloricBracketError(
                "combined dry-sorbate/free-water branch lacks positive capacity"
            )
        endpoint_tolerance = 2.0e-12 * max(abs(energy_lo), abs(energy_hi), 1.0)
        if (
            native_envelope_energy_j < energy_lo - endpoint_tolerance
            or native_envelope_energy_j > energy_hi + endpoint_tolerance
        ):
            raise DrySorbateCaloricBracketError(
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
        particle_energy = native_envelope_energy_j - water_energy
        return DrySorbateParticleState(
            dry_meal_mass_kg=dry_matter_mass_kg,
            oil_label_mass_kg=oil_label_mass_kg,
            retained_water_mass_kg=retained_water_mass_kg,
            sorbed_hexane_mass_kg=sorbed_hexane_mass_kg,
            total_internal_energy_j=particle_energy,
            temperature_k=temperature,
        )

    def _component_masses(
        self,
        native: DrySorbateParticleState,
        free_water_mass_kg: float,
    ) -> dict[str, float]:
        return {
            "dry_matter_kg": native.dry_meal_mass_kg,
            "residual_oil_label_kg": native.oil_label_mass_kg,
            "total_hexane_kg": native.sorbed_hexane_mass_kg,
            "total_water_kg": math.fsum((native.retained_water_mass_kg, free_water_mass_kg)),
        }

    def _native_envelope_energy(
        self,
        native: DrySorbateParticleState,
        free_water_mass_kg: float,
    ) -> float:
        try:
            water_internal_energy = water.state_Tp(
                native.temperature_k,
                self.pressure_pa,
                "liquid",
            ).u_mass
        except ValueError as error:
            raise DrySorbatePacketStateError(str(error)) from error
        value = math.fsum(
            (native.total_internal_energy_j, free_water_mass_kg * water_internal_energy)
        )
        if not math.isfinite(value):
            raise DrySorbatePacketStateError("packet envelope energy is not representable")
        return value


def envelope_memo_cache_info() -> dict:
    """The tag-3 memo census, in the SAME field names tag-1 reports.

    The march driver reads the tag-1 census off the class attribute
    (``scripts/run_engineering_mechanism_march.py``, key
    ``envelope_memo_cache_info``); this parallel module-level accessor lets a
    reader collect the tag-3 census the same way without editing the driver.
    """

    return DrySorbatePacketAdapter._invert_envelope.cache_info()  # noqa: SLF001


__all__ = (
    "CompleteDrySorbatePacket",
    "DRY_SORBATE_ADAPTER_DIGEST_DOMAIN",
    "DRY_SORBATE_ENVELOPE_MEMO_ARGUMENT_ORDER",
    "DRY_SORBATE_SCHEMA_VERSION",
    "DRY_SORBATE_V1_REGISTRY",
    "DRY_SORBATE_V1_SCHEMA",
    "DRY_SORBATE_VARIANT_TAG",
    "DRY_SORBATE_VARIANT_TOKEN",
    "DrySorbateCaloricBracketError",
    "DrySorbatePacketAdapter",
    "DrySorbatePacketAdapterError",
    "DrySorbatePacketConfigurationError",
    "DrySorbatePacketContractDriftError",
    "DrySorbatePacketStateError",
    "DrySorbateVariantHandoffRequired",
    "envelope_memo_cache_info",
)
