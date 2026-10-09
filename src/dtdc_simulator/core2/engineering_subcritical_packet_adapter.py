"""TAG-2 packet codec: the subcritical (film-exhausted) particle variant.

Owner-ruled GO 2026-09-01 (morning ruling, disposition 1); design record
`GT_PS2_TAG2_ADAPTER_DESIGN_RECORD_2026-09-01.md`.  Faithful mirror of the
tag-1 `HighLoadingV1PacketAdapter` contract on the tag-2 native state
(`engineering_subcritical_packet_state`, seam-exact caloric):

- SAME schema_id and wire-contract generation; **variant_tag 2** per the
  repo convention that the variant tag denotes the PHYSICAL particle
  variant (the native-datum v2 adapter's own statement).
- SAME context pins in every payload (configuration digest, energy datum,
  fixed pressure, component identities, component-datum definition and
  numeric offsets) - drift refuses identically.
- The envelope block keeps all seven extensives with
  ``envelope_attached_hexane_kg`` REQUIRED exactly +0.0 (canonical
  positive zero), so `PerParticleInventory`, the surface composite, and
  every conservation ledger consume tag-2 packets unchanged.
- The native block drops only the attached-hexane field: a tag-2 particle
  that re-grows a film hands back to tag-1 at the DUAL-adapter boundary,
  never mutates in place.

Nothing here edits tag-1: its schema digest, registry digest, payload
bytes, and canonical identities are untouched (measured in the recon).

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
from .engineering_subcritical_packet_state import (
    SubcriticalParticleState,
    SubcriticalSorptionFloorRefusal,
    SubcriticalStateError,
    SubcriticalSupercriticalRefusal,
    initialize_subcritical_state,
    subcritical_specific_energy,
    subcritical_temperature_bracket,
)
from .envelope_inversion_memo import exact_envelope_memo
from .particle import sphere
from .props import water
from .packet_weight_scaling_contract import (
    ExternalSurfaceComposite,
    PerParticleInventory,
    require_single_surface_authority,
)

SUBCRITICAL_SCHEMA_VERSION = 1
SUBCRITICAL_VARIANT_TAG = 2
SUBCRITICAL_VARIANT_TOKEN = "engineering_subcritical_v1"
SUBCRITICAL_ADAPTER_DIGEST_DOMAIN = "GT-PS-2/engineering-subcritical-v1-packet-adapter/v1"

#: Keyword order of the tag-2 caloric inversion, and therefore the packing
#: order of the exact-argument memo key.  Same construction as tag-1 and
#: tag-3; see ``envelope_inversion_memo``.
SUBCRITICAL_ENVELOPE_MEMO_ARGUMENT_ORDER = (
    "dry_matter_mass_kg",
    "oil_label_mass_kg",
    "retained_water_mass_kg",
    "internal_hexane_mass_kg",
    "external_free_water_mass_kg",
    "envelope_energy_j",
)


class SubcriticalPacketAdapterError(ValueError):
    """Base typed refusal for the detached tag-2 adapter."""


class SubcriticalPacketConfigurationError(SubcriticalPacketAdapterError):
    """The fixed caloric/configuration authority is malformed."""


class SubcriticalPacketContractDriftError(SubcriticalPacketAdapterError):
    """Payload or host configuration, datum, pressure, or components drifted."""


class SubcriticalPacketStateError(SubcriticalPacketAdapterError):
    """A tag-2 payload is noncanonical or violates the native state map."""


class SubcriticalCaloricBracketError(SubcriticalPacketStateError):
    """Combined packet energy has no admissible bracketed temperature."""


class SubcriticalVariantHandoffRequired(SubcriticalPacketStateError):
    """The state left the tag-2 band; the tag-1 (film) variant owns it."""


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


SUBCRITICAL_V1_SCHEMA = state_codec.StateSchema(
    schema_id=tag1.SCHEMA_ID,
    schema_version=SUBCRITICAL_SCHEMA_VERSION,
    variant_tag=SUBCRITICAL_VARIANT_TAG,
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
        tag1._mass("state_internal_hexane_mass_kg"),
        state_codec.ScalarField("state_total_internal_energy_j", state_codec.ScalarKind.FLOAT64),
        state_codec.ScalarField(
            "state_temperature_k", state_codec.ScalarKind.FLOAT64, nonnegative=True
        ),
        tag1._mass("surface_attached_hexane_kg"),
        tag1._mass("surface_free_water_kg"),
    ),
)

SUBCRITICAL_V1_REGISTRY = state_codec.StateSchemaRegistry(
    registry_id="gt.ps2.particle-state.engineering-subcritical-v1",
    schemas=(SUBCRITICAL_V1_SCHEMA,),
)


@dataclass(frozen=True, slots=True, kw_only=True)
class CompleteSubcriticalPacket:
    """Validated tag-2 state, surface, envelope, and canonical bytes."""

    state: SubcriticalParticleState
    surface: ExternalSurfaceComposite
    inventory: PerParticleInventory
    payload_bytes: bytes
    canonical_identity: str

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


@dataclass(frozen=True, slots=True, kw_only=True)
class SubcriticalPacketAdapter:
    """One context-pinned, engineering-only tag-2 codec and auditor."""

    configuration_digest: str
    energy_datum_id: str
    pressure_pa: float
    component_datum_adapter: component_datum.NumericComponentEnergyDatumAdapter
    sphere_params: sphere.SphereParams = sphere.SphereParams()
    component_ids: tag1.HighLoadingComponentIds = tag1.HighLoadingComponentIds()
    source_identity: str = "detached-engineering-subcritical-v1-adapter"

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_codec_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        # The context gates are exactly the tag-1 gates; run them through a
        # throwaway tag-1 adapter so the two variants can never drift apart.
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
            raise SubcriticalPacketConfigurationError(str(error)) from error
        if type(self.source_identity) is not str or not self.source_identity.strip():
            raise SubcriticalPacketConfigurationError(
                "adapter source identity must be a nonblank exact str"
            )

    @property
    def codec_registry_digest(self) -> str:
        return SUBCRITICAL_V1_REGISTRY.definition_digest

    @property
    def auditor_identity_digest(self) -> str:
        return _framed_digest(
            SUBCRITICAL_ADAPTER_DIGEST_DOMAIN + "/auditor",
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
            raise SubcriticalPacketContractDriftError(
                "packet and host configuration digests differ"
            )
        if energy_datum_id != self.energy_datum_id:
            raise SubcriticalPacketContractDriftError("packet and host energy datums differ")
        if type(pressure_pa) is not float or not _same_float(pressure_pa, self.pressure_pa):
            raise SubcriticalPacketContractDriftError(
                "packet and host fixed caloric pressures differ"
            )

    def initialize(
        self,
        *,
        internal_hexane_loading: float,
        temperature_k: float,
        external_free_water_loading: float,
        retained_water_loading: float | None = None,
        oil_fraction: float | None = None,
    ) -> CompleteSubcriticalPacket:
        if (
            type(external_free_water_loading) is not float
            or not math.isfinite(external_free_water_loading)
            or external_free_water_loading < 0.0
        ):
            raise SubcriticalPacketStateError(
                "external free-water loading must be finite nonnegative binary64"
            )
        try:
            state = initialize_subcritical_state(
                self.sphere_params,
                internal_hexane_loading=internal_hexane_loading,
                temperature_k=temperature_k,
                retained_water_loading=retained_water_loading,
                oil_fraction=oil_fraction,
            )
        except SubcriticalStateError as error:
            raise SubcriticalPacketStateError(str(error)) from error
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
    ) -> CompleteSubcriticalPacket:
        """Issue a canonical tag-2 packet from exact conserved quantities."""

        for name, value in (
            ("dry matter", dry_matter_mass_kg),
            ("oil label", oil_label_mass_kg),
            ("retained water", retained_water_mass_kg),
            ("total hexane", total_hexane_mass_kg),
            ("external free water", external_free_water_mass_kg),
            ("envelope energy", envelope_energy_j),
        ):
            if type(value) is not float or not math.isfinite(value):
                raise SubcriticalPacketStateError(f"{name} must be finite exact binary64")
        if (
            dry_matter_mass_kg <= 0.0
            or oil_label_mass_kg < 0.0
            or retained_water_mass_kg < 0.0
            or total_hexane_mass_kg <= 0.0
            or external_free_water_mass_kg < 0.0
        ):
            raise SubcriticalPacketStateError(
                "reconstituted packet masses left the tag-2 domain; reject, do not clamp"
            )
        state = self._invert_envelope(
            dry_matter_mass_kg=dry_matter_mass_kg,
            oil_label_mass_kg=oil_label_mass_kg,
            retained_water_mass_kg=retained_water_mass_kg,
            internal_hexane_mass_kg=total_hexane_mass_kg,
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
    ) -> CompleteSubcriticalPacket:
        """Advance only conserved envelope totals, then invert the state.

        ``retained_water_mass_after_kg`` (owner ruling D8-d, 2026-09-12; the
        C8 Tier 1 keyword tag-1 and tag-3 already carry) is forwarded to
        :meth:`reconstitute` ONLY when supplied; when it is ``None`` the
        template's retained water is pinned exactly as before the ruling, so
        every call written without the keyword is byte-identical by
        construction.  ``reconstitute`` validates the moved value (finite,
        non-negative) and consumes it calorically; nothing is re-derived
        here.
        """

        before = self.decode(payload_bytes)
        return self.reconstitute(
            dry_matter_mass_kg=before.state.dry_meal_mass_kg,
            oil_label_mass_kg=before.state.oil_label_mass_kg,
            retained_water_mass_kg=(
                before.state.retained_water_mass_kg
                if retained_water_mass_after_kg is None
                else retained_water_mass_after_kg
            ),
            total_hexane_mass_kg=total_hexane_mass_after_kg,
            external_free_water_mass_kg=external_free_water_mass_after_kg,
            envelope_energy_j=envelope_energy_after_j,
        )

    def encode(self, inventory: PerParticleInventory) -> bytes:
        if type(inventory) is not PerParticleInventory:
            raise SubcriticalPacketStateError("encode requires an exact PerParticleInventory")
        if inventory.attached_hexane_kg != 0.0:
            raise SubcriticalVariantHandoffRequired("a positive attached film is tag-1 territory")
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

    def decode(self, payload_bytes: bytes) -> CompleteSubcriticalPacket:
        try:
            decoded = state_codec.decode_state(payload_bytes, SUBCRITICAL_V1_REGISTRY)
        except (TypeError, ValueError) as error:
            raise SubcriticalPacketStateError(str(error)) from error
        values = decoded.values
        if len(values) != 30:  # schema construction makes this defensive only
            raise SubcriticalPacketStateError("tag-2 payload field count is inconsistent")
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
        if variant_token != SUBCRITICAL_VARIANT_TOKEN:
            raise SubcriticalPacketContractDriftError("tag-2 variant token drifted")
        if configuration_digest != self.configuration_digest:
            raise SubcriticalPacketContractDriftError(
                "payload configuration digest drifted from the adapter pin"
            )
        if energy_datum_id != self.energy_datum_id:
            raise SubcriticalPacketContractDriftError(
                "payload energy datum drifted from the adapter pin"
            )
        if type(pressure_pa) is not float or not _same_float(pressure_pa, self.pressure_pa):
            raise SubcriticalPacketContractDriftError(
                "payload fixed pressure drifted from the adapter pin"
            )
        if component_ids != self.component_ids.as_tuple():
            raise SubcriticalPacketContractDriftError(
                "payload component identities drifted from the adapter pin"
            )
        if component_datum_digest != self.component_datum_adapter.definition_digest:
            raise SubcriticalPacketContractDriftError(
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
            raise SubcriticalPacketContractDriftError(
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
            native = SubcriticalParticleState(
                dry_meal_mass_kg=state_values[0],
                oil_label_mass_kg=state_values[1],
                retained_water_mass_kg=state_values[2],
                internal_hexane_mass_kg=state_values[3],
                total_internal_energy_j=state_values[4],
                temperature_k=state_values[5],
            )
            surface = ExternalSurfaceComposite(
                attached_hexane_kg=surface_values[0],
                free_water_kg=surface_values[1],
            )
        except (TypeError, ValueError) as error:
            raise SubcriticalPacketStateError(str(error)) from error
        self._validate_complete(native, surface, inventory)
        return CompleteSubcriticalPacket(
            state=native,
            surface=surface,
            inventory=inventory,
            payload_bytes=decoded.canonical_bytes,
            canonical_identity=decoded.identity_digest,
        )

    def _issue(
        self,
        native: SubcriticalParticleState,
        *,
        free_water_mass_kg: float,
        envelope_energy_j: float | None = None,
    ) -> CompleteSubcriticalPacket:
        if type(native) is not SubcriticalParticleState:
            raise SubcriticalPacketStateError("tag-2 issue requires exact SubcriticalParticleState")
        if (
            type(free_water_mass_kg) is not float
            or not math.isfinite(free_water_mass_kg)
            or free_water_mass_kg < 0.0
        ):
            raise SubcriticalPacketStateError(
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
            internal_hexane_kg=native.internal_hexane_mass_kg,
            external_water_kg=free_water_mass_kg,
            retained_water_kg=native.retained_water_mass_kg,
            common_datum_energy_j=selected_energy,
            energy_datum_id=self.energy_datum_id,
        )
        self._validate_complete(native, surface, inventory)
        values: tuple[state_codec.CodecValue, ...] = (
            SUBCRITICAL_VARIANT_TOKEN,
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
            native.total_internal_energy_j,
            native.temperature_k,
            surface.attached_hexane_kg,
            surface.free_water_kg,
        )
        payload = state_codec.encode_state(SUBCRITICAL_V1_SCHEMA, values)
        return self.decode(payload)

    def _validate_complete(
        self,
        native: SubcriticalParticleState,
        surface: ExternalSurfaceComposite,
        inventory: PerParticleInventory,
    ) -> None:
        if inventory.attached_hexane_kg != 0.0 or surface.attached_hexane_kg != 0.0:
            raise SubcriticalVariantHandoffRequired("a positive attached film is tag-1 territory")
        if inventory.energy_datum_id != self.energy_datum_id:
            raise SubcriticalPacketContractDriftError("packet envelope energy datum drifted")
        for label, envelope_value, state_value in (
            ("dry matter", inventory.dry_matter_kg, native.dry_meal_mass_kg),
            ("residual oil", inventory.residual_oil_label_kg, native.oil_label_mass_kg),
            ("internal hexane", inventory.internal_hexane_kg, native.internal_hexane_mass_kg),
            ("retained water", inventory.retained_water_kg, native.retained_water_mass_kg),
        ):
            if not _same_float(envelope_value, state_value):
                raise SubcriticalPacketStateError(
                    f"envelope and native-state {label} authorities disagree"
                )
        try:
            require_single_surface_authority(inventory, surface)
        except ValueError as error:
            raise SubcriticalPacketStateError(str(error)) from error
        try:
            rebuilt = initialize_subcritical_state(
                self.sphere_params,
                internal_hexane_loading=native.internal_hexane_loading,
                temperature_k=native.temperature_k,
                total_internal_energy_j=native.total_internal_energy_j,
                retained_water_loading=native.retained_water_loading,
                oil_fraction=native.oil_fraction,
            )
        except SubcriticalSupercriticalRefusal as error:
            raise SubcriticalVariantHandoffRequired(str(error)) from error
        except SubcriticalStateError as error:
            raise SubcriticalPacketStateError(str(error)) from error
        exact_pairs = (
            (native.dry_meal_mass_kg, rebuilt.dry_meal_mass_kg),
            (native.oil_label_mass_kg, rebuilt.oil_label_mass_kg),
            (native.retained_water_mass_kg, rebuilt.retained_water_mass_kg),
            (native.total_internal_energy_j, rebuilt.total_internal_energy_j),
            (native.temperature_k, rebuilt.temperature_k),
        )
        if any(not _same_float(left, right) for left, right in exact_pairs):
            raise SubcriticalPacketStateError(
                "decoded tag-2 primitives are not the exact native initialized state"
            )
        # Unlike tag-1 there is no phase partition to reconstruct: internal
        # hexane IS the total, so the loading survives one divide/multiply.
        partition_tolerance = 2.0e-13 * native.dry_meal_mass_kg
        if abs(native.internal_hexane_mass_kg - rebuilt.internal_hexane_mass_kg) > (
            partition_tolerance
        ):
            raise SubcriticalPacketStateError(
                "decoded tag-2 hexane inventory violates the native initialized state"
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
            raise SubcriticalPacketStateError(
                "envelope energy is inconsistent with subcritical plus free-water energy"
            )

    @exact_envelope_memo(*SUBCRITICAL_ENVELOPE_MEMO_ARGUMENT_ORDER)
    def _invert_envelope(
        self,
        *,
        dry_matter_mass_kg: float,
        oil_label_mass_kg: float,
        retained_water_mass_kg: float,
        internal_hexane_mass_kg: float,
        external_free_water_mass_kg: float,
        envelope_energy_j: float,
    ) -> SubcriticalParticleState:
        if not _same_float(dry_matter_mass_kg, self.sphere_params.dry_meal_mass_kg):
            raise SubcriticalPacketStateError(
                "packet dry matter drifted from representative-sphere geometry"
            )
        retained_loading = retained_water_mass_kg / dry_matter_mass_kg
        oil_fraction = oil_label_mass_kg / dry_matter_mass_kg
        internal_loading = internal_hexane_mass_kg / dry_matter_mass_kg
        native_envelope_energy_j = self.component_datum_adapter.common_to_native_energy_j(
            envelope_energy_j,
            dry_matter_kg=dry_matter_mass_kg,
            residual_oil_label_kg=oil_label_mass_kg,
            total_hexane_kg=internal_hexane_mass_kg,
            total_water_kg=math.fsum((retained_water_mass_kg, external_free_water_mass_kg)),
        )
        try:
            wet = sphere._material_wet_params(self.sphere_params, retained_loading, oil_fraction)
        except ValueError as error:
            raise SubcriticalPacketStateError(str(error)) from error
        try:
            lo, hi = subcritical_temperature_bracket(internal_loading, wet)
        except SubcriticalSupercriticalRefusal as error:
            raise SubcriticalVariantHandoffRequired(str(error)) from error
        except SubcriticalSorptionFloorRefusal:
            raise
        except SubcriticalStateError as error:
            raise SubcriticalCaloricBracketError(str(error)) from error

        def energy_at(temperature_k: float) -> float:
            particle = dry_matter_mass_kg * subcritical_specific_energy(
                temperature_k, internal_loading, wet
            )
            water_energy = (
                external_free_water_mass_kg
                * water.state_Tp(temperature_k, self.pressure_pa, "liquid").u_mass
            )
            return particle + water_energy

        energy_lo = energy_at(lo)
        energy_hi = energy_at(hi)
        if energy_hi <= energy_lo:
            raise SubcriticalCaloricBracketError(
                "combined subcritical/free-water branch lacks positive capacity"
            )
        endpoint_tolerance = 2.0e-12 * max(abs(energy_lo), abs(energy_hi), 1.0)
        if (
            native_envelope_energy_j < energy_lo - endpoint_tolerance
            or native_envelope_energy_j > energy_hi + endpoint_tolerance
        ):
            raise SubcriticalCaloricBracketError(
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
        return SubcriticalParticleState(
            dry_meal_mass_kg=dry_matter_mass_kg,
            oil_label_mass_kg=oil_label_mass_kg,
            retained_water_mass_kg=retained_water_mass_kg,
            internal_hexane_mass_kg=internal_hexane_mass_kg,
            total_internal_energy_j=particle_energy,
            temperature_k=temperature,
        )

    def _component_masses(
        self,
        native: SubcriticalParticleState,
        free_water_mass_kg: float,
    ) -> dict[str, float]:
        return {
            "dry_matter_kg": native.dry_meal_mass_kg,
            "residual_oil_label_kg": native.oil_label_mass_kg,
            "total_hexane_kg": native.internal_hexane_mass_kg,
            "total_water_kg": math.fsum((native.retained_water_mass_kg, free_water_mass_kg)),
        }

    def _native_envelope_energy(
        self,
        native: SubcriticalParticleState,
        free_water_mass_kg: float,
    ) -> float:
        try:
            water_internal_energy = water.state_Tp(
                native.temperature_k,
                self.pressure_pa,
                "liquid",
            ).u_mass
        except ValueError as error:
            raise SubcriticalPacketStateError(str(error)) from error
        value = math.fsum(
            (native.total_internal_energy_j, free_water_mass_kg * water_internal_energy)
        )
        if not math.isfinite(value):
            raise SubcriticalPacketStateError("packet envelope energy is not representable")
        return value


def envelope_memo_cache_info() -> dict:
    """The tag-2 memo census, in the SAME field names tag-1 reports."""

    return SubcriticalPacketAdapter._invert_envelope.cache_info()  # noqa: SLF001


__all__ = (
    "CompleteSubcriticalPacket",
    "SUBCRITICAL_ADAPTER_DIGEST_DOMAIN",
    "SUBCRITICAL_ENVELOPE_MEMO_ARGUMENT_ORDER",
    "SUBCRITICAL_SCHEMA_VERSION",
    "SUBCRITICAL_V1_REGISTRY",
    "SUBCRITICAL_V1_SCHEMA",
    "SUBCRITICAL_VARIANT_TAG",
    "SUBCRITICAL_VARIANT_TOKEN",
    "SubcriticalCaloricBracketError",
    "SubcriticalPacketAdapter",
    "SubcriticalPacketAdapterError",
    "SubcriticalPacketConfigurationError",
    "SubcriticalPacketContractDriftError",
    "SubcriticalPacketStateError",
    "SubcriticalVariantHandoffRequired",
    "envelope_memo_cache_info",
)
