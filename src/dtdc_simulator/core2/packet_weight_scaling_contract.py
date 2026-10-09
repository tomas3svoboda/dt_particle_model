"""Manufactured packet-weight scaling contract for accepted solid transfers.

This is the B1 gate for PR-06. The T3C3 real-variant inventory states that a
revised, versioned auditor / accepted-state / population-bridge scaling contract
"is required before PR-04 split legality or ``AcceptedSolidTransfer``
implementation"
(``docs/GT_PS2_T3C3_REAL_PARTICLE_STATE_VARIANT_INVENTORY_2026-08-09.md:199-200``).
The owner ruled B1 on 2026-08-10
(``docs/GT_PS2_B1_B2_B3_RULING_RECORD_2026-08-10.md``): the contract is authored
as a NEW versioned detached module in its own packet, and the two frozen seam
oracles stay byte-intact. Nothing here modifies them.

THE DEFECT THIS CONTRACT REPAIRS. At the existing seam the accepted state's
extensive basis is a cohort TOTAL, and the accepted state welds the payload's own
extensives to that total by exact equality. Identical payload bytes therefore
force an identical basis dry matter, with no ratio or multiplicity anywhere to
absorb a difference -- so two PR-04 split children cannot carry different
positive packet weights while sharing bytes, which is exactly what PR-02/PR-04
require. See ``docs/GT_PS2_B1_SCALING_CONTRACT_DESIGN_BASIS_2026-08-10.md`` for
the mechanism and the five downstream seams that re-impose it.

THE DIRECTION IS FIXED BY ARITHMETIC, NOT PREFERENCE. Carrier scaling in the
frozen oracles is elementwise binary64 with no exact-rational arithmetic, so
scale-then-unscale is not a byte identity: a weight-invariant payload cannot be
recovered from a weight-scaled total by division. Therefore this contract stores
PER-PARTICLE values and multiplies OUTWARD to totals. It never stores totals and
divides back, and it never mutates payload bytes.

DERIVED, NOT WRITABLE. Multiplicity is a derived read-only quantity here, never a
second writable authority (inventory ``:194-195``).

DETACHMENT IS MANDATORY HERE, NOT STYLISTIC. Two of the seam modules this
contract reasons about are among the four sibling module names the frozen ``src/``
guard tests scan for, and those scans read raw file text. This module therefore
declares its own seven-extensive type locally, mirroring the frozen field names
and canonical order as value strings, and imports neither of those modules. It
names none of the four guarded module names anywhere, including in prose. The
enumeration lives in ``docs/GT_PS2_B4_B5_RULING_RECORD_2026-08-10.md``.

VALIDATION LADDER. Ruled Tier-B 2026-08-11
(``docs/GT_PS2_VALIDATION_LADDER_RULING_RECORD_2026-08-11.md``): the STRICTER of
the two frozen ladders applies. Negative zero is rejected on the six mass
extensives; ``common_datum_energy_j`` stays sign-free because datum-relative
energy may legitimately be negative; and the phase-partition representability
guard is inherited so overflow is refused at this seam rather than surfacing
downstream as an infinity. This narrows admissibility relative to the permissive
frozen oracle, and the narrowing is declared, not hidden -- see
``NARROWED_VERSUS_PERMISSIVE_ORACLE``.

This is manufactured architecture evidence. It authenticates no source, selects
no production registry, wires no production caller, moves no solid, and closes no
gate. It does not implement the PR-06 transaction, its host lock, its persisted
replay protection, or any physical transport.
"""

from __future__ import annotations

import enum
import hashlib
import math
import struct
from dataclasses import dataclass
from typing import ClassVar

CONTRACT_MAGIC = "GT-PS-2-PACKET-WEIGHT-SCALING-CONTRACT"
CONTRACT_VERSION = 1
CONTRACT_DIGEST_DOMAIN = "GT-PS-2/pr06-packet-weight-scaling/v1"

#: Canonical order of the seven frozen carrier extensives, mirrored as value
#: strings from the frozen carrier authority rather than imported from it.
EXTENSIVE_NAMES: tuple[str, ...] = (
    "dry_matter_kg",
    "residual_oil_label_kg",
    "attached_hexane_kg",
    "internal_hexane_kg",
    "external_water_kg",
    "retained_water_kg",
    "common_datum_energy_j",
)
#: The six mass extensives, i.e. every extensive except the datum-relative
#: energy. Only these are sign-constrained.
MASS_EXTENSIVE_COUNT = 6

#: Declared narrowing versus the permissive frozen carrier authority, per the
#: 2026-08-11 validation-ladder ruling. Recorded so the divergence is auditable.
NARROWED_VERSUS_PERMISSIVE_ORACLE: tuple[str, ...] = (
    "negative_zero_rejected_on_six_mass_extensives",
    "phase_partition_totals_must_be_representable",
)


class ScalingContractError(ValueError):
    """Base exception for a refused scaling-contract operation."""


class WeightAuthorityError(ScalingContractError):
    """Raised when a packet weight is not a strictly positive binary64 mass."""


class ReconciliationError(ScalingContractError):
    """Raised when declared totals do not reconcile with the derived totals."""


class SurfaceAuthorityError(ScalingContractError):
    """Raised when an external-surface duplicate disagrees with the envelope."""


def _is_negative_zero(value: float) -> bool:
    return struct.pack(">d", value) == b"\x80\x00\x00\x00\x00\x00\x00\x00"


def _require_binary64(name: str, *values: float) -> None:
    for value in values:
        if type(value) is not float:
            raise ScalingContractError(f"{name} must contain exact binary64 float values")
        if not math.isfinite(value):
            raise ScalingContractError(f"{name} must contain only finite values")


def _require_nonnegative_mass(name: str, value: float) -> None:
    """Strict ladder: reject both negative values and negative zero."""

    if value < 0.0 or _is_negative_zero(value):
        raise ScalingContractError(f"{name} must be non-negative with canonical positive zero")


def _framed_digest(domain: str, parts: tuple[bytes, ...]) -> str:
    digest = hashlib.sha256()
    for part in (domain.encode("ascii"), *parts):
        digest.update(len(part).to_bytes(8, "big"))
        digest.update(part)
    return "sha256:" + digest.hexdigest()


def _float_part(value: float) -> bytes:
    """Digest a binary64 with -0.0 canonicalized, mirroring the frozen digest."""

    canonical = 0.0 if value == 0.0 else value
    return struct.pack(">d", canonical)


class DeclaredEngineeringAssumption(enum.Enum):
    """Typed declarations; the disposition forbids silent choices."""

    #: Multiplicity is derived from the weight and the per-particle dry matter and
    #: is exposed read-only. No document states a multiplicity representation.
    DERIVED_READ_ONLY_MULTIPLICITY = "derived_read_only_multiplicity_v1"


@dataclass(frozen=True, slots=True, kw_only=True)
class PerParticleInventory:
    """Packet-weight-invariant per-particle extensive inventory.

    These are the values that live inside canonical payload bytes. They are
    per representative particle, never cohort totals, which is what makes two
    split children with identical bytes able to carry different weights.
    """

    dry_matter_kg: float
    residual_oil_label_kg: float
    attached_hexane_kg: float
    internal_hexane_kg: float
    external_water_kg: float
    retained_water_kg: float
    common_datum_energy_j: float
    energy_datum_id: str

    architecture_only: ClassVar[bool] = True
    packet_weight_invariant: ClassVar[bool] = True
    negative_zero_rejected_on_masses: ClassVar[bool] = True

    def __post_init__(self) -> None:
        values = self.as_tuple()
        _require_binary64("per-particle inventory", *values)
        if type(self.energy_datum_id) is not str or not self.energy_datum_id.strip():
            raise ScalingContractError("per-particle inventory energy_datum_id must be nonblank")
        for name, value in zip(EXTENSIVE_NAMES[:MASS_EXTENSIVE_COUNT], values, strict=False):
            _require_nonnegative_mass(name, value)
        if self.residual_oil_label_kg > self.dry_matter_kg:
            raise ScalingContractError("residual-oil label cannot exceed dry matter")
        # Inherited phase-partition representability guard.
        for label, pair in (
            ("hexane", (self.attached_hexane_kg, self.internal_hexane_kg)),
            ("water", (self.external_water_kg, self.retained_water_kg)),
        ):
            try:
                total = math.fsum(pair)
            except OverflowError as error:  # pragma: no cover - defensive
                raise ScalingContractError(
                    f"per-particle {label} partition total is not representable"
                ) from error
            if not math.isfinite(total):
                raise ScalingContractError(
                    f"per-particle {label} partition total is not representable"
                )

    def as_tuple(self) -> tuple[float, ...]:
        """The seven extensives in canonical order; the datum label is excluded."""

        return (
            self.dry_matter_kg,
            self.residual_oil_label_kg,
            self.attached_hexane_kg,
            self.internal_hexane_kg,
            self.external_water_kg,
            self.retained_water_kg,
            self.common_datum_energy_j,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PacketWeight:
    """Host-carried external positive dry-matter packet weight.

    This never lives inside payload bytes. It is the only authority that differs
    between two split children sharing identical bytes.
    """

    representative_particles: float
    energy_datum_id: str

    architecture_only: ClassVar[bool] = True
    host_carried: ClassVar[bool] = True

    def __post_init__(self) -> None:
        _require_binary64("packet weight", self.representative_particles)
        if self.representative_particles <= 0.0:
            raise WeightAuthorityError("packet weight must be strictly positive")
        if type(self.energy_datum_id) is not str or not self.energy_datum_id.strip():
            raise WeightAuthorityError("packet weight energy_datum_id must be nonblank")


@dataclass(frozen=True, slots=True, kw_only=True)
class ExternalSurfaceComposite:
    """Per-particle extensive external-surface composite (B2 binding).

    The certified composite is intensive, on a per-kg-dry basis. B2 binds it into
    the accepted envelope on the per-particle EXTENSIVE basis, so the conversion
    is an exact multiply by the payload's own per-particle dry matter. Surface
    topology stays derived, never stored.
    """

    attached_hexane_kg: float
    free_water_kg: float

    architecture_only: ClassVar[bool] = True
    topology_is_derived: ClassVar[bool] = True

    def __post_init__(self) -> None:
        _require_binary64("external surface composite", self.attached_hexane_kg, self.free_water_kg)
        _require_nonnegative_mass("attached_hexane_kg", self.attached_hexane_kg)
        _require_nonnegative_mass("free_water_kg", self.free_water_kg)

    @property
    def topology(self) -> str:
        """Derived on read, mirroring the certified topology vocabulary."""

        has_hexane = self.attached_hexane_kg > 0.0
        has_water = self.free_water_kg > 0.0
        if has_hexane and has_water:
            return "two_liquid_vlle"
        if has_hexane:
            return "hexane_only"
        if has_water:
            return "water_only"
        return "dry"


@dataclass(frozen=True, slots=True, kw_only=True)
class PacketTotals:
    """Derived cohort totals: per-particle values multiplied outward by weight."""

    totals: tuple[float, ...]
    energy_datum_id: str
    derived_multiplicity: float
    contract_digest: str

    architecture_only: ClassVar[bool] = True
    derived_outward_from_per_particle: ClassVar[bool] = True
    payload_bytes_unmutated: ClassVar[bool] = True
    # Explicitly not supplied by this contract.
    production_codec_implemented: ClassVar[bool] = False
    persisted_replay_protection_selected: ClassVar[bool] = False
    accepted_solid_transfer_implemented: ClassVar[bool] = False
    complete_particle_state_transport_selected: ClassVar[bool] = False
    physically_qualifying: ClassVar[bool] = False
    plant_predictive: ClassVar[bool] = False


def external_surface_from_intensive(
    inventory: PerParticleInventory,
    attached_hexane_kg_per_kg_dry: float,
    free_water_kg_per_kg_dry: float,
) -> ExternalSurfaceComposite:
    """Convert the certified intensive composite to the per-particle basis.

    Refuses at zero dry matter: the intensive basis is undefined there, and the
    conversion must not silently yield zero surface for an unknown state.
    """

    _require_binary64(
        "intensive surface composite",
        attached_hexane_kg_per_kg_dry,
        free_water_kg_per_kg_dry,
    )
    _require_nonnegative_mass("attached_hexane_kg_per_kg_dry", attached_hexane_kg_per_kg_dry)
    _require_nonnegative_mass("free_water_kg_per_kg_dry", free_water_kg_per_kg_dry)
    if inventory.dry_matter_kg <= 0.0:
        raise ScalingContractError(
            "intensive surface conversion requires strictly positive per-particle dry matter"
        )
    return ExternalSurfaceComposite(
        attached_hexane_kg=attached_hexane_kg_per_kg_dry * inventory.dry_matter_kg,
        free_water_kg=free_water_kg_per_kg_dry * inventory.dry_matter_kg,
    )


def require_single_surface_authority(
    inventory: PerParticleInventory, surface: ExternalSurfaceComposite
) -> None:
    """One writable authority: a variant duplicate must equal the envelope exactly.

    B2 admits a variant carrying its own copy of an external quantity only when it
    agrees bit-for-bit with the envelope value; otherwise there would be two
    writable authorities for one physical quantity. BOTH fields of the certified
    two-field composite are covered -- checking only one would leave the other
    silently dual-authority.

    NAME MAPPING, stated explicitly rather than inferred by resemblance (the B1
    design basis flags this as a required mapping): the certified composite's
    ``attached_hexane_kg`` maps to the envelope's ``attached_hexane_kg``, and its
    ``free_water_kg`` -- external surface liquid -- maps to the envelope's
    ``external_water_kg``. ``retained_water_kg`` is internal retained water and is
    NOT a surface quantity, so it is deliberately not part of this binding.
    """

    for label, envelope_value, surface_value in (
        ("attached hexane", inventory.attached_hexane_kg, surface.attached_hexane_kg),
        ("external free water", inventory.external_water_kg, surface.free_water_kg),
    ):
        if struct.pack(">d", envelope_value) != struct.pack(">d", surface_value):
            raise SurfaceAuthorityError(
                f"envelope {label} disagrees with the external-surface composite"
            )


#: Exact-bits memo for ``derive_packet_totals`` (2026-08-30 speed mandate).
#: Keys are the IEEE-754 bytes of every input float plus both datum identity
#: strings, so ``-0.0``/``0.0`` and equal-int/float merges cannot alias and a
#: forged field always misses.  Only successful derivations are stored;
#: refusals re-run and re-raise on the uncached path every time.
_TOTALS_CACHE: dict[tuple[str, str, bytes], PacketTotals] = {}
_TOTALS_CACHE_ENTRY_LIMIT = 500_000
_totals_cache_hits = 0


def derive_packet_totals(inventory: PerParticleInventory, weight: PacketWeight) -> PacketTotals:
    """Multiply per-particle values outward by the packet weight.

    Never divides, never mutates the inventory, and refuses on datum mismatch or
    on any total that is not representable.  Successful derivations are memoized
    on the exact bit patterns of the inputs; the cached object is the same frozen
    ``PacketTotals`` the uncached path constructed for those exact bytes.
    """

    global _totals_cache_hits
    key: tuple[str, str, bytes] | None = None
    if type(inventory) is PerParticleInventory and type(weight) is PacketWeight:
        try:
            key = (
                inventory.energy_datum_id,
                weight.energy_datum_id,
                struct.pack(">8d", *inventory.as_tuple(), weight.representative_particles),
            )
        except struct.error, TypeError:
            key = None
        if key is not None:
            hit = _TOTALS_CACHE.get(key)
            if hit is not None:
                _totals_cache_hits += 1
                return hit
    result = _derive_packet_totals_uncached(inventory, weight)
    if key is not None and len(_TOTALS_CACHE) < _TOTALS_CACHE_ENTRY_LIMIT:
        _TOTALS_CACHE[key] = result
    return result


def _derive_packet_totals_uncached(
    inventory: PerParticleInventory, weight: PacketWeight
) -> PacketTotals:
    if type(inventory) is not PerParticleInventory:
        raise ScalingContractError("inventory must be an exact PerParticleInventory")
    if type(weight) is not PacketWeight:
        raise ScalingContractError("weight must be an exact PacketWeight")
    if inventory.energy_datum_id != weight.energy_datum_id:
        raise ScalingContractError(
            "cannot scale a per-particle inventory against a foreign energy datum"
        )
    factor = weight.representative_particles
    totals = tuple(factor * value for value in inventory.as_tuple())
    _require_binary64("derived packet totals", *totals)
    for name, value in zip(EXTENSIVE_NAMES[:MASS_EXTENSIVE_COUNT], totals, strict=False):
        _require_nonnegative_mass(f"derived {name}", value)
    digest = _framed_digest(
        CONTRACT_DIGEST_DOMAIN,
        (
            CONTRACT_MAGIC.encode("ascii"),
            CONTRACT_VERSION.to_bytes(8, "big"),
            DeclaredEngineeringAssumption.DERIVED_READ_ONLY_MULTIPLICITY.value.encode("ascii"),
            inventory.energy_datum_id.encode("utf-8"),
            *(_float_part(value) for value in inventory.as_tuple()),
            _float_part(factor),
        ),
    )
    return PacketTotals(
        totals=totals,
        energy_datum_id=inventory.energy_datum_id,
        derived_multiplicity=factor,
        contract_digest=digest,
    )


def reconcile_declared_totals(
    inventory: PerParticleInventory,
    weight: PacketWeight,
    declared_totals: tuple[float, ...],
) -> tuple[float, ...]:
    """Reconcile host-declared totals against the derived totals.

    Returns the componentwise residuals. Reconciliation is read-only: it never
    mutates payload bytes and never adjusts either side to agree. A nonzero
    residual is reported, not repaired -- the house rule is that a shortfall is
    "never repaired by a clamp".
    """

    if type(declared_totals) is not tuple or len(declared_totals) != len(EXTENSIVE_NAMES):
        raise ReconciliationError(f"declared totals must be an exact {len(EXTENSIVE_NAMES)}-tuple")
    _require_binary64("declared totals", *declared_totals)
    derived = derive_packet_totals(inventory, weight)
    return tuple(
        math.fsum((declared, -value))
        for declared, value in zip(declared_totals, derived.totals, strict=True)
    )


__all__ = (
    "CONTRACT_DIGEST_DOMAIN",
    "CONTRACT_MAGIC",
    "CONTRACT_VERSION",
    "DeclaredEngineeringAssumption",
    "EXTENSIVE_NAMES",
    "ExternalSurfaceComposite",
    "MASS_EXTENSIVE_COUNT",
    "NARROWED_VERSUS_PERMISSIVE_ORACLE",
    "PacketTotals",
    "PacketWeight",
    "PerParticleInventory",
    "ReconciliationError",
    "ScalingContractError",
    "SurfaceAuthorityError",
    "WeightAuthorityError",
    "derive_packet_totals",
    "external_surface_from_intensive",
    "reconcile_declared_totals",
    "require_single_surface_authority",
)
