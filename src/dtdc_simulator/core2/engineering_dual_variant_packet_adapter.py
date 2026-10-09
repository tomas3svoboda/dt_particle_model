"""The variant-dispatching packet codec: one authority for tag-1, tag-2 and tag-3.

Owner-ruled GO 2026-09-01 (tag-2); owner-ruled GO 2026-09-06 (tag-3, "I rule
all four as recommended, proceed."); design records
`GT_PS2_TAG2_ADAPTER_DESIGN_RECORD_2026-09-01.md` piece B and
`GT_PS2_TAG3_DRY_SORBATE_PACKET_BUILD_SPEC_2026-09-06.md` section 3.  A tray
pins ONE codec registry digest and ONE auditor identity across all its packets
(measured recon constraint), so mixed film/subcritical/dry-sorbate populations
require ONE codec object presenting ONE combined registry.  This module is
that object: it holds the three sealed sub-adapters, routes by the payload's
own variant token, and owns the HANDOFFS:

- tag-1 -> tag-2: a conserved-total advance/reconstitute whose state falls
  below pore saturation re-issues under tag-2 from the SAME six conserved
  quantities (the caloric seam is exact by construction - zero latent
  impulse; measured in the tag-2 state bites).
- tag-2 -> tag-1: a state whose pores re-saturate (condensation/heating)
  re-issues under tag-1 (film reborn).
- tag-2 -> tag-3: a state that falls to or below the a_h = 1 sorption gate
  W_eq(1, T) re-issues under tag-3.  Until 2026-09-06 that landing was the
  typed `SubcriticalSorptionFloorRefusal`; it is now the handoff, and the
  seam is again exact by construction (the tag-3 deficit AT the gate is the
  frozen memoized `saturated_binding_deficit` object itself).
- tag-3 -> tag-2: a state that rises strictly above the gate (hexane
  condensed onto a dry packet, or the packet heated past its own gate)
  re-issues under tag-2.  Re-wetting, the mirror crossing.
- tag-3 -> tag-1 DIRECTLY is impossible: a film needs the loading above
  X_c(T), which no single crossing of the sorption gate can reach.  It is
  refused typed (`DualVariantDirectFilmRebirthRefusal`), never routed.

Sealed tag-1-only history is untouched: this codec is ADOPTED by new runs
at fixture construction; every old payload, schema digest, and canonical
identity is byte-identical (measured).

physically_qualifying: false.  plant_predictive: false.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import ClassVar, Union

from . import component_energy_datum_adapter as component_datum
from . import engineering_dry_sorbate_packet_adapter as tag3
from . import engineering_dry_sorbate_packet_state as dry_state
from . import engineering_subcritical_packet_adapter as tag2
from . import engineering_subcritical_packet_state as sub_state
from . import high_loading_packet_adapter as tag1
from . import state_codec
from .packet_weight_scaling_contract import PerParticleInventory
from .particle import sphere

DUAL_VARIANT_ADAPTER_DIGEST_DOMAIN = "GT-PS-2/engineering-dual-variant-packet-adapter/v1"

#: The combined registry.  TAG-3 (owner ruling 2026-09-06) adds a third
#: schema, so both the registry_id and the registry digest move.
#: PROVENANCE, superseded 2026-09-06:
#:   registry_id "gt.ps2.particle-state.engineering-film-and-subcritical-v1"
#:   definition_digest
#:   "sha256:b996b8238e9a33e0cd33786b0d7fc510d0b34a2b28fd2cee2625529454e1a06b"
#: Adding a schema to a registry changes ONLY the registry digest (measured,
#: TAG-2 design record constraint 2): the tag-1 and tag-2 schema digests,
#: payload bytes and canonical payload identities are untouched.
ENGINEERING_DUAL_VARIANT_REGISTRY = state_codec.StateSchemaRegistry(
    registry_id="gt.ps2.particle-state.engineering-film-subcritical-and-dry-sorbate-v1",
    schemas=(
        tag1.HIGH_LOADING_V1_SCHEMA,
        tag2.SUBCRITICAL_V1_SCHEMA,
        tag3.DRY_SORBATE_V1_SCHEMA,
    ),
)

CompletePacket = Union[
    tag1.CompleteHighLoadingPacket,
    tag2.CompleteSubcriticalPacket,
    tag3.CompleteDrySorbatePacket,
]


class DualVariantAdapterError(ValueError):
    """Base typed refusal for the dual-variant codec."""


class DualVariantConfigurationError(DualVariantAdapterError):
    """The two sub-adapters do not share one exact caloric context."""


class DualVariantPayloadError(DualVariantAdapterError):
    """A payload carries no known engineering variant token."""


class DualVariantDirectFilmRebirthRefusal(DualVariantAdapterError):
    """A tag-3 payload was asked to land as tag-1 in one step.

    A film needs the loading above X_c(T); a dry-sorbate packet is below the
    a_h = 1 gate W_eq(1, T), which is two orders of magnitude below X_c.  No
    single lawful crossing spans that, so the request is a defect in the
    caller's conserved totals, not a seam.  Refuse, never route.
    """


def _framed_digest(domain: str, parts: tuple[bytes, ...]) -> str:
    digest = hashlib.sha256()
    for part in (domain.encode("ascii"), *parts):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class EngineeringDualVariantPacketAdapter:
    """One context-pinned codec speaking every engineering particle variant.

    The class NAME is kept (every march-path codec gate admits this exact
    type; TAG-3 adds a variant, not a type).  ``tag3_codec`` is required: a
    tray pins ONE registry digest, and a two-variant object would present a
    different one from the three-variant registry this module publishes.
    """

    tag1_codec: tag1.HighLoadingV1PacketAdapter
    tag2_codec: tag2.SubcriticalPacketAdapter
    tag3_codec: tag3.DrySorbatePacketAdapter

    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False
    production_codec_implemented: ClassVar[bool] = False

    def __post_init__(self) -> None:
        if type(self.tag1_codec) is not tag1.HighLoadingV1PacketAdapter:
            raise DualVariantConfigurationError(
                "tag1_codec must be an exact HighLoadingV1PacketAdapter"
            )
        if type(self.tag2_codec) is not tag2.SubcriticalPacketAdapter:
            raise DualVariantConfigurationError(
                "tag2_codec must be an exact SubcriticalPacketAdapter"
            )
        if type(self.tag3_codec) is not tag3.DrySorbatePacketAdapter:
            raise DualVariantConfigurationError(
                "tag3_codec must be an exact DrySorbatePacketAdapter"
            )
        others = (self.tag2_codec, self.tag3_codec)
        for other in others:
            pairs = (
                (
                    "configuration digest",
                    self.tag1_codec.configuration_digest,
                    other.configuration_digest,
                ),
                ("energy datum", self.tag1_codec.energy_datum_id, other.energy_datum_id),
            )
            for label, left, right in pairs:
                if left != right:
                    raise DualVariantConfigurationError(
                        f"sub-adapters bind different {label} authorities"
                    )
            if self.tag1_codec.pressure_pa != other.pressure_pa:
                raise DualVariantConfigurationError(
                    "sub-adapters bind different fixed caloric pressures"
                )
            if self.tag1_codec.component_ids != other.component_ids:
                raise DualVariantConfigurationError(
                    "sub-adapters bind different component identities"
                )
            if (
                self.tag1_codec.component_datum_adapter.definition_digest
                != other.component_datum_adapter.definition_digest
            ):
                raise DualVariantConfigurationError(
                    "sub-adapters bind different component-datum definitions"
                )
            if repr(self.tag1_codec.sphere_params) != repr(other.sphere_params):
                raise DualVariantConfigurationError(
                    "sub-adapters bind different representative-sphere geometry"
                )

    # -- shared context surface (host codec protocol) ----------------------

    @property
    def configuration_digest(self) -> str:
        return self.tag1_codec.configuration_digest

    @property
    def energy_datum_id(self) -> str:
        return self.tag1_codec.energy_datum_id

    @property
    def pressure_pa(self) -> float:
        return self.tag1_codec.pressure_pa

    @property
    def component_datum_adapter(self) -> component_datum.NumericComponentEnergyDatumAdapter:
        return self.tag1_codec.component_datum_adapter

    @property
    def component_ids(self) -> tag1.HighLoadingComponentIds:
        return self.tag1_codec.component_ids

    @property
    def sphere_params(self) -> sphere.SphereParams:
        return self.tag1_codec.sphere_params

    @property
    def source_identity(self) -> str:
        return "detached-engineering-dual-variant-adapter"

    @property
    def codec_registry_digest(self) -> str:
        return ENGINEERING_DUAL_VARIANT_REGISTRY.definition_digest

    @property
    def auditor_identity_digest(self) -> str:
        # TAG-3 2026-09-06: the third sub-identity joins the frame, so this
        # digest moves with the registry digest.  Both are computed from the
        # module's own code; no known answer anywhere pins either literal.
        return _framed_digest(
            DUAL_VARIANT_ADAPTER_DIGEST_DOMAIN + "/auditor",
            (
                self.codec_registry_digest.encode("ascii"),
                self.tag1_codec.auditor_identity_digest.encode("ascii"),
                self.tag2_codec.auditor_identity_digest.encode("ascii"),
                self.tag3_codec.auditor_identity_digest.encode("ascii"),
            ),
        )

    def require_host_binding(
        self,
        *,
        configuration_digest: str,
        energy_datum_id: str,
        pressure_pa: float,
    ) -> None:
        self.tag1_codec.require_host_binding(
            configuration_digest=configuration_digest,
            energy_datum_id=energy_datum_id,
            pressure_pa=pressure_pa,
        )

    # -- variant routing ---------------------------------------------------

    def _variant_token(self, payload_bytes: bytes) -> str:
        try:
            decoded = state_codec.decode_state(payload_bytes, ENGINEERING_DUAL_VARIANT_REGISTRY)
        except (TypeError, ValueError) as error:
            raise DualVariantPayloadError(str(error)) from error
        token = decoded.values[0]
        if token not in (
            tag1.HIGH_LOADING_VARIANT_TOKEN,
            tag2.SUBCRITICAL_VARIANT_TOKEN,
            tag3.DRY_SORBATE_VARIANT_TOKEN,
        ):
            raise DualVariantPayloadError(
                f"payload carries an unknown engineering variant token {token!r}"
            )
        return token

    def _codec_for(self, token: str):
        if token == tag1.HIGH_LOADING_VARIANT_TOKEN:
            return self.tag1_codec
        if token == tag2.SUBCRITICAL_VARIANT_TOKEN:
            return self.tag2_codec
        return self.tag3_codec

    def decode(self, payload_bytes: bytes) -> CompletePacket:
        return self._codec_for(self._variant_token(payload_bytes)).decode(payload_bytes)

    def audit(self, payload_bytes: bytes) -> PerParticleInventory:
        return self.decode(payload_bytes).inventory

    def temperature_from_payload(self, payload_bytes: bytes) -> float:
        return self.decode(payload_bytes).state.temperature_k

    # -- issue paths with the two handoffs ---------------------------------

    def initialize(self, **kwargs: object) -> CompletePacket:
        """Variant-selecting initializer (fixture/host codec protocol).

        ``total_hexane_loading`` selects the tag-1 (film) initializer,
        ``internal_hexane_loading`` the tag-2 (subcritical) one and
        ``sorbed_hexane_loading`` the tag-3 (dry-sorbate) one - the three
        signatures share no loading keyword, so the selection is exact.
        """

        if "total_hexane_loading" in kwargs:
            return self.tag1_codec.initialize(**kwargs)  # type: ignore[arg-type]
        if "internal_hexane_loading" in kwargs:
            return self.tag2_codec.initialize(**kwargs)  # type: ignore[arg-type]
        if "sorbed_hexane_loading" in kwargs:
            return self.tag3_codec.initialize(**kwargs)  # type: ignore[arg-type]
        raise DualVariantAdapterError("initialize requires a variant-selecting loading argument")

    def encode(self, inventory: PerParticleInventory) -> bytes:
        if type(inventory) is not PerParticleInventory:
            raise DualVariantAdapterError("encode requires an exact PerParticleInventory")
        if inventory.attached_hexane_kg > 0.0:
            return self.tag1_codec.encode(inventory)
        try:
            return self.tag2_codec.encode(inventory)
        except (
            tag2.SubcriticalPacketAdapterError,
            sub_state.SubcriticalSorptionFloorRefusal,
        ) as tag2_refusal:
            # Below the a_h = 1 gate there is no mobile liquid at all: the
            # dry-sorbate continuation owns the state.  A double refusal
            # re-raises the ORIGINAL, so no domain is widened.
            try:
                return self.tag3_codec.encode(inventory)
            except (
                tag3.DrySorbatePacketAdapterError,
                dry_state.DrySorbateStateError,
            ):
                raise tag2_refusal from None

    def reconstitute(
        self,
        *,
        dry_matter_mass_kg: float,
        oil_label_mass_kg: float,
        retained_water_mass_kg: float,
        total_hexane_mass_kg: float,
        external_free_water_mass_kg: float,
        envelope_energy_j: float,
    ) -> CompletePacket:
        """Issue whichever variant the conserved quantities lawfully land on."""

        # The three variants partition the (loading, temperature) plane at two
        # crossings: pore saturation X_c(T) between tag-1 and tag-2, and the
        # a_h = 1 sorption gate W_eq(1, T) between tag-2 and tag-3.  In the
        # OVERLAP bands a variant refuses a genuinely neighbouring state with
        # its CALORIC BRACKET error rather than the typed handoff (its
        # re-bracketed endpoint pins a temperature the conserved energy
        # contradicts) - so BOTH refusal types route onward, and the ENERGY
        # selects the lawful variant.  If no side accepts, the ORIGINAL tag-1
        # refusal re-raises: nothing is widened, the union of the three exact
        # domains is unchanged.
        quantities = dict(
            dry_matter_mass_kg=dry_matter_mass_kg,
            oil_label_mass_kg=oil_label_mass_kg,
            retained_water_mass_kg=retained_water_mass_kg,
            total_hexane_mass_kg=total_hexane_mass_kg,
            external_free_water_mass_kg=external_free_water_mass_kg,
            envelope_energy_j=envelope_energy_j,
        )
        try:
            return self.tag1_codec.reconstitute(**quantities)
        except (
            tag1.HighLoadingVariantHandoffRequired,
            tag1.HighLoadingCaloricBracketError,
        ) as tag1_refusal:
            try:
                return self.tag2_codec.reconstitute(**quantities)
            except (
                tag2.SubcriticalPacketAdapterError,
                sub_state.SubcriticalSorptionFloorRefusal,
            ):
                try:
                    return self.tag3_codec.reconstitute(**quantities)
                except (
                    tag3.DrySorbatePacketAdapterError,
                    dry_state.DrySorbateStateError,
                ):
                    raise tag1_refusal from None

    def advance(
        self,
        payload_bytes: bytes,
        *,
        total_hexane_mass_after_kg: float,
        external_free_water_mass_after_kg: float,
        envelope_energy_after_j: float,
        retained_water_mass_after_kg: float | None = None,
    ) -> CompletePacket:
        """Advance conserved totals; cross the variant seams when they demand it.

        The template's dry matter, oil label, and retained water are pinned
        exactly as each sub-adapter pins them; a handoff re-issues from the
        SAME six conserved quantities on the other side of the seam, so the
        crossing conserves mass and energy bit-for-bit by construction.

        ``retained_water_mass_after_kg`` (C8 Tier 1) is forwarded ONLY when
        supplied, so every call written before 2026-09-05 reaches each
        sub-adapter with the exact argument list it had then.  All three
        sub-adapters carry the keyword: tag-1 and tag-3 since the C8 Tier 1
        build, tag-2 since owner ruling D8-d (2026-09-12), which closed the
        asymmetry under which a tag-1 handoff already landed a moved retained
        water on a tag-2 packet while a tag-2 packet advancing as tag-2
        refused.
        """

        token = self._variant_token(payload_bytes)
        keyword = (
            {}
            if retained_water_mass_after_kg is None
            else {"retained_water_mass_after_kg": retained_water_mass_after_kg}
        )
        before = self._codec_for(token).decode(payload_bytes)
        pinned = dict(
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
        if token == tag1.HIGH_LOADING_VARIANT_TOKEN:
            try:
                return self.tag1_codec.advance(
                    payload_bytes,
                    total_hexane_mass_after_kg=total_hexane_mass_after_kg,
                    external_free_water_mass_after_kg=external_free_water_mass_after_kg,
                    envelope_energy_after_j=envelope_energy_after_j,
                    **keyword,
                )
            except (
                tag1.HighLoadingVariantHandoffRequired,
                tag1.HighLoadingCaloricBracketError,
            ) as tag1_refusal:
                # Same overlap-band routing as reconstitute: the energy
                # selects the side; a triple refusal re-raises the original.
                try:
                    return self.tag2_codec.reconstitute(**pinned)
                except (
                    tag2.SubcriticalPacketAdapterError,
                    sub_state.SubcriticalSorptionFloorRefusal,
                ):
                    try:
                        return self.tag3_codec.reconstitute(**pinned)
                    except (
                        tag3.DrySorbatePacketAdapterError,
                        dry_state.DrySorbateStateError,
                    ):
                        raise tag1_refusal from None
        if token == tag2.SUBCRITICAL_VARIANT_TOKEN:
            try:
                return self.tag2_codec.advance(
                    payload_bytes,
                    total_hexane_mass_after_kg=total_hexane_mass_after_kg,
                    external_free_water_mass_after_kg=external_free_water_mass_after_kg,
                    envelope_energy_after_j=envelope_energy_after_j,
                    **keyword,
                )
            except (
                tag2.SubcriticalVariantHandoffRequired,
                tag2.SubcriticalCaloricBracketError,
                sub_state.SubcriticalSorptionFloorRefusal,
            ) as tag2_refusal:
                # The refusal names the direction: the sorption-floor refusal
                # IS the tag-2 -> tag-3 handoff and is tried DOWNWARD first;
                # everything else is tried upward first.  Both sides are tried
                # either way, and a double refusal re-raises the original.
                downward_first = isinstance(
                    tag2_refusal, sub_state.SubcriticalSorptionFloorRefusal
                )
                order = (
                    (self.tag3_codec, self.tag1_codec)
                    if downward_first
                    else (self.tag1_codec, self.tag3_codec)
                )
                for candidate in order:
                    try:
                        return candidate.reconstitute(**pinned)
                    except (
                        tag1.HighLoadingPacketAdapterError,
                        tag3.DrySorbatePacketAdapterError,
                        dry_state.DrySorbateStateError,
                    ):
                        continue
                raise tag2_refusal from None
        try:
            return self.tag3_codec.advance(
                payload_bytes,
                total_hexane_mass_after_kg=total_hexane_mass_after_kg,
                external_free_water_mass_after_kg=external_free_water_mass_after_kg,
                envelope_energy_after_j=envelope_energy_after_j,
                **keyword,
            )
        except (
            tag3.DrySorbateVariantHandoffRequired,
            tag3.DrySorbateCaloricBracketError,
        ) as tag3_refusal:
            try:
                return self.tag2_codec.reconstitute(**pinned)
            except (
                tag2.SubcriticalPacketAdapterError,
                sub_state.SubcriticalSorptionFloorRefusal,
            ) as tag2_refusal:
                # tag-3 -> tag-1 directly is impossible; if the totals only
                # make sense as a film, say so typed instead of routing.
                if isinstance(tag2_refusal, tag2.SubcriticalVariantHandoffRequired):
                    raise DualVariantDirectFilmRebirthRefusal(
                        "a dry-sorbate packet cannot re-grow an attached film in "
                        "one crossing: the conserved totals land above pore "
                        "saturation X_c(T) from below the a_h = 1 gate"
                    ) from tag3_refusal
                raise tag3_refusal from None


def make_engineering_dual_codec(
    *,
    configuration_digest: str,
    energy_datum_id: str,
    pressure_pa: float,
    component_datum_adapter: component_datum.NumericComponentEnergyDatumAdapter,
    sphere_params: sphere.SphereParams | None = None,
    component_ids: tag1.HighLoadingComponentIds | None = None,
) -> EngineeringDualVariantPacketAdapter:
    """Construct all three sealed sub-adapters from ONE caloric context."""

    kwargs = dict(
        configuration_digest=configuration_digest,
        energy_datum_id=energy_datum_id,
        pressure_pa=pressure_pa,
        component_datum_adapter=component_datum_adapter,
    )
    if sphere_params is not None:
        kwargs["sphere_params"] = sphere_params
    if component_ids is not None:
        kwargs["component_ids"] = component_ids
    return EngineeringDualVariantPacketAdapter(
        tag1_codec=tag1.HighLoadingV1PacketAdapter(**kwargs),
        tag2_codec=tag2.SubcriticalPacketAdapter(**kwargs),
        tag3_codec=tag3.DrySorbatePacketAdapter(**kwargs),
    )


__all__ = (
    "DUAL_VARIANT_ADAPTER_DIGEST_DOMAIN",
    "DualVariantAdapterError",
    "DualVariantConfigurationError",
    "DualVariantDirectFilmRebirthRefusal",
    "DualVariantPayloadError",
    "ENGINEERING_DUAL_VARIANT_REGISTRY",
    "EngineeringDualVariantPacketAdapter",
    "make_engineering_dual_codec",
)
