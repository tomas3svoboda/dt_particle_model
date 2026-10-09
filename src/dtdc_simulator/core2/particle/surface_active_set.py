r"""Frozen Packet-C external-surface active-set foundation.

This module implements only the authority-complete algebraic layer:

* exact inventory selection of dry, n-hexane-only, water-only, and
  co-located two-liquid VLLE surface topologies;
* the combined external-liquid occupancy and single shared gas-side area;
* pure-liquid versus binary-virial gas fugacity residuals with exact-sign
  one-sided classifications;
* the transparent/impermeable water-film n-hexane-access bracket; and
* one component-flux vector and one common-datum energy flux applied over the
  same area exactly once.

Packet C and PHY-022 do not freeze a complete finite-thickness external-film
coefficient law.  Consequently this foundation does not invent one.  A
caller supplies a coupled film flux vector and explicit interval selections
from its separately qualified film solver.  The selections are retained as
provenance but are not converted here into an unauthorized linear rate law.

Every water-present branch applies the two endpoint access oracles.  The
water-only topology still concerns an *external* film: this module has no
internal-front argument and never applies VLLE between that film and a
spatially internal n-hexane interface.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

from dtdc_simulator.core2.props import binary_gas as bg
from dtdc_simulator.core2.props import hexane as hx
from dtdc_simulator.core2.props import water as wa


class SurfaceTopology(str, Enum):
    DRY = "dry_surface"
    HEXANE_ONLY = "external_hexane_only"
    WATER_ONLY = "external_water_only"
    TWO_LIQUID_VLLE = "external_two_liquid_vlle"


class Component(str, Enum):
    WATER = "water"
    HEXANE = "hexane"


class GasClosure(str, Enum):
    PRESSURE_SECOND_VIRIAL = "pressure_second_virial"
    IDEAL_ABLATION = "ideal_vapor_ablation"


class HexaneAccessOracle(str, Enum):
    TRANSPARENT_INTERFACE = "transparent_interface"
    IMPERMEABLE_WATER_FILM = "impermeable_water_film"


class ResidualSign(str, Enum):
    NEGATIVE = "negative"
    EXACT_ZERO = "exact_zero"
    POSITIVE = "positive"


class PhaseResidualClassification(str, Enum):
    ABSENT_STABLE = "absent_stable"
    ABSENT_ONE_SIDED_CONTACT = "absent_one_sided_contact"
    ABSENT_APPEARANCE_DEMANDED = "absent_appearance_demanded"
    ACTIVE_GAS_UNDERSATURATED = "active_gas_undersaturated"
    ACTIVE_EQUALITY_SATISFIED = "active_equality_satisfied"
    ACTIVE_GAS_SUPERSATURATED = "active_gas_supersaturated"


class InventoryEndpointClassification(str, Enum):
    ABSENT_CONTINUATION = "absent_continuation"
    APPEARANCE = "appearance"
    DISAPPEARANCE = "disappearance"
    PRESENT_CONTINUATION = "present_continuation"


class HexanePrimaryDrainageHistory(str, Enum):
    """Caller-owned irreversible n-hexane drainage-history latch."""

    PRE_CORE_RECESSION = "pre_core_recession"
    CORE_RECESSION_OCCURRED = "core_recession_occurred"
    CORE_DRY_OUT_OCCURRED = "core_dry_out_occurred"


class PrimaryDrainageEnvelopeError(ValueError):
    """A proposed surface state would reverse the qualified drainage history."""


class TransferCoefficientKind(str, Enum):
    BINARY_MAXWELL_STEFAN_FILM = "binary_maxwell_stefan_film"
    HEAT_TRANSFER = "heat_transfer"


@dataclass(frozen=True)
class ExternalSurfaceInventories:
    """External liquid loadings on kg per kg dry-composite-meal basis."""

    attached_hexane_kg_per_kg_dry: float
    free_water_kg_per_kg_dry: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (
            self.attached_hexane_kg_per_kg_dry,
            self.free_water_kg_per_kg_dry,
        )
        if not all(math.isfinite(value) and value >= 0.0 for value in values):
            raise ValueError("external liquid inventories must be finite and non-negative")

    @property
    def topology(self) -> SurfaceTopology:
        has_hexane = self.attached_hexane_kg_per_kg_dry > 0.0
        has_water = self.free_water_kg_per_kg_dry > 0.0
        if has_hexane and has_water:
            return SurfaceTopology.TWO_LIQUID_VLLE
        if has_hexane:
            return SurfaceTopology.HEXANE_ONLY
        if has_water:
            return SurfaceTopology.WATER_ONLY
        return SurfaceTopology.DRY


@dataclass(frozen=True)
class SurfaceGeometry:
    """Frozen bed/particle geometry used by Packet-C occupancy."""

    particle_radius_m: float
    bed_void_fraction: float
    dry_meal_particle_density_kg_m3: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        values = (self.particle_radius_m, self.dry_meal_particle_density_kg_m3)
        if not all(math.isfinite(value) and value > 0.0 for value in values):
            raise ValueError("surface radius and dry-meal density must be positive")
        if (
            not math.isfinite(self.bed_void_fraction)
            or not 0.0 < self.bed_void_fraction < 1.0
        ):
            raise ValueError("bed void fraction must lie strictly inside (0,1)")

    @property
    def dry_meal_concentration_per_bed_volume_kg_m3(self) -> float:
        return (1.0 - self.bed_void_fraction) * self.dry_meal_particle_density_kg_m3

    @property
    def dry_meal_concentration_per_gas_void_kg_m3(self) -> float:
        return (
            self.dry_meal_concentration_per_bed_volume_kg_m3
            / self.bed_void_fraction
        )

    @property
    def geometric_bed_area_m2_m3(self) -> float:
        return (1.0 - self.bed_void_fraction) * 3.0 / self.particle_radius_m


@dataclass(frozen=True)
class ExternalLiquidOccupancy:
    """Exact Packet-C external occupancies and the single shared area."""

    dry_meal_concentration_per_gas_void_kg_m3: float
    liquid_hexane_density_kg_m3: float
    liquid_water_density_kg_m3: float
    hexane_saturation: float
    water_saturation: float
    total_liquid_saturation: float
    geometric_bed_area_m2_m3: float
    shared_gas_side_area_m2_m3: float
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class SurfaceGasThermo:
    """Gas fugacity and partial enthalpy on one selected closure."""

    closure: GasClosure
    temperature_k: float
    pressure_pa: float
    y_water: float
    y_hexane: float
    fugacity_water_pa: float
    fugacity_hexane_pa: float
    water_partial_enthalpy_j_mol: float
    hexane_partial_enthalpy_j_mol: float
    residual_enthalpy_molar_j_mol: float
    below_cross_direct_evidence: bool
    outside_project_pressure: bool
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class ComponentPhaseCondition:
    """Exact-sign pure-liquid/gas fugacity residual for one component."""

    component: Component
    external_liquid_inventory_kg_per_kg_dry: float
    liquid_present: bool
    gas_fugacity_pa: float
    pure_liquid_fugacity_pa: float
    gas_minus_liquid_fugacity_residual_pa: float
    residual_sign: ResidualSign
    classification: PhaseResidualClassification
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class SurfaceActiveSet:
    """One inventory-selected external interface state."""

    inventories: ExternalSurfaceInventories
    geometry: SurfaceGeometry
    topology: SurfaceTopology
    occupancy: ExternalLiquidOccupancy
    gas: SurfaceGasThermo
    water_condition: ComponentPhaseCondition
    hexane_condition: ComponentPhaseCondition
    active_external_liquid_components: tuple[Component, ...]
    hexane_access_oracle: HexaneAccessOracle | None
    hexane_access_multiplier: float
    external_two_liquid_vlle_active: bool
    couples_internal_hexane_front_to_external_water_vlle: bool = field(
        default=False,
        init=False,
    )
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class OneSidedPhaseEvent:
    """Exact zero/nonzero inventory transition with endpoint residual evidence."""

    component: Component
    previous_inventory_kg_per_kg_dry: float
    next_inventory_kg_per_kg_dry: float
    inventory_classification: InventoryEndpointClassification
    endpoint_residual_sign: ResidualSign
    endpoint_phase_classification: PhaseResidualClassification
    physically_qualifying: bool = field(default=False, init=False)


@dataclass(frozen=True)
class SurfaceTransferCoefficientInterval:
    """Caller-owned strictly positive coefficient interval."""

    kind: TransferCoefficientKind
    lower: float
    upper: float
    units: str
    provenance: str
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.kind, TransferCoefficientKind):
            raise TypeError("surface coefficient kind must be explicit")
        if not all(math.isfinite(value) and value > 0.0 for value in (self.lower, self.upper)):
            raise ValueError("surface coefficient interval must be strictly positive")
        if self.lower > self.upper:
            raise ValueError("surface coefficient interval must be ordered")
        if not isinstance(self.units, str) or not self.units.strip():
            raise ValueError("surface coefficient interval needs declared units")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("surface coefficient interval needs provenance")

    def select(self, fraction: float) -> SurfaceTransferCoefficientSelection:
        if not math.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
            raise ValueError("surface coefficient selection fraction must lie in [0,1]")
        value = self.lower + fraction * (self.upper - self.lower)
        return SurfaceTransferCoefficientSelection(self, fraction, value)


@dataclass(frozen=True)
class SurfaceTransferCoefficientSelection:
    """One explicit point within a caller-supplied interval."""

    interval: SurfaceTransferCoefficientInterval
    fraction: float
    value: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.interval, SurfaceTransferCoefficientInterval):
            raise TypeError("surface coefficient selection requires its interval")
        if not math.isfinite(self.fraction) or not 0.0 <= self.fraction <= 1.0:
            raise ValueError("surface coefficient selection fraction must lie in [0,1]")
        expected = self.interval.lower + self.fraction * (
            self.interval.upper - self.interval.lower
        )
        if not math.isfinite(self.value) or self.value != expected:
            raise ValueError(
                "surface coefficient selection value must match its interval fraction"
            )


@dataclass(frozen=True)
class SurfaceTransferSelections:
    """Explicit provenance for the separately supplied mass/heat fluxes."""

    binary_film: SurfaceTransferCoefficientSelection
    heat_transfer: SurfaceTransferCoefficientSelection
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.binary_film, SurfaceTransferCoefficientSelection):
            raise TypeError("binary film coefficient selection is required")
        if not isinstance(self.heat_transfer, SurfaceTransferCoefficientSelection):
            raise TypeError("heat-transfer coefficient selection is required")
        if self.binary_film.interval.kind is not TransferCoefficientKind.BINARY_MAXWELL_STEFAN_FILM:
            raise ValueError("binary-film selection has the wrong coefficient kind")
        if self.heat_transfer.interval.kind is not TransferCoefficientKind.HEAT_TRANSFER:
            raise ValueError("heat-transfer selection has the wrong coefficient kind")


@dataclass(frozen=True)
class SurfaceComponentFluxVector:
    """One outward-positive coupled gas-film component vector, mol/m2/s."""

    water_mol_m2_s: float
    hexane_mol_m2_s: float
    physically_qualifying: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        if not all(math.isfinite(value) for value in (self.water_mol_m2_s, self.hexane_mol_m2_s)):
            raise ValueError("surface component flux vector must remain finite")


@dataclass(frozen=True)
class SurfaceTransferResult:
    """One area-applied component/common-datum energy transfer ledger."""

    active_set: SurfaceActiveSet
    coefficients: SurfaceTransferSelections
    raw_component_flux_density: SurfaceComponentFluxVector
    effective_component_flux_density: SurfaceComponentFluxVector
    conductive_heat_flux_w_m2: float
    shared_area_m2_m3: float
    gas_water_source_mol_m3_s: float
    gas_hexane_source_mol_m3_s: float
    material_water_source_mol_m3_s: float
    material_hexane_source_mol_m3_s: float
    water_enthalpy_flux_w_m2: float
    hexane_enthalpy_flux_w_m2: float
    total_energy_flux_w_m2: float
    gas_energy_source_w_m3: float
    material_energy_source_w_m3: float
    physically_qualifying: bool = field(default=False, init=False)

    @property
    def component_pair_residuals_mol_m3_s(self) -> tuple[float, float]:
        return (
            math.fsum((self.gas_water_source_mol_m3_s, self.material_water_source_mol_m3_s)),
            math.fsum((self.gas_hexane_source_mol_m3_s, self.material_hexane_source_mol_m3_s)),
        )

    @property
    def energy_pair_residual_w_m3(self) -> float:
        return math.fsum((self.gas_energy_source_w_m3, self.material_energy_source_w_m3))


def evaluate_surface_active_set(
    inventories: ExternalSurfaceInventories,
    geometry: SurfaceGeometry,
    *,
    temperature_k: float,
    pressure_pa: float,
    y_water: float,
    y_hexane: float,
    gas_closure: GasClosure = GasClosure.PRESSURE_SECOND_VIRIAL,
    k_wh: float = bg.K_WH_CENTRAL,
    external_water_hexane_access: HexaneAccessOracle | None = None,
) -> SurfaceActiveSet:
    """Evaluate one external active set without an internal-front coupling."""
    if not isinstance(inventories, ExternalSurfaceInventories):
        raise TypeError("external surface inventories are required")
    if not isinstance(geometry, SurfaceGeometry):
        raise TypeError("surface geometry is required")
    if not isinstance(gas_closure, GasClosure):
        raise TypeError("gas closure must be explicit")
    topology = inventories.topology
    external_water_present = inventories.free_water_kg_per_kg_dry > 0.0
    if external_water_present:
        if not isinstance(external_water_hexane_access, HexaneAccessOracle):
            raise ValueError(
                "every external-water topology requires a transparent or impermeable "
                "n-hexane-access oracle"
            )
        access_oracle = external_water_hexane_access
        access_multiplier = (
            1.0
            if access_oracle is HexaneAccessOracle.TRANSPARENT_INTERFACE
            else 0.0
        )
    else:
        if external_water_hexane_access is not None:
            raise ValueError("water-film access oracle requires external water inventory")
        access_oracle = None
        access_multiplier = 1.0

    gas = _surface_gas_thermo(
        temperature_k,
        pressure_pa,
        y_water,
        y_hexane,
        gas_closure,
        k_wh,
    )
    liquid_water = wa.state_Tp(temperature_k, pressure_pa, "liquid")
    liquid_hexane = hx.state_Tp(temperature_k, pressure_pa, "liquid")
    occupancy = external_liquid_occupancy(
        inventories,
        geometry,
        liquid_hexane_density_kg_m3=liquid_hexane.rho_mass,
        liquid_water_density_kg_m3=liquid_water.rho_mass,
    )
    water_condition = classify_phase_condition(
        Component.WATER,
        inventories.free_water_kg_per_kg_dry,
        gas.fugacity_water_pa,
        liquid_water.fugacity,
    )
    hexane_condition = classify_phase_condition(
        Component.HEXANE,
        inventories.attached_hexane_kg_per_kg_dry,
        gas.fugacity_hexane_pa,
        liquid_hexane.fugacity,
    )
    if topology is SurfaceTopology.DRY:
        active_components: tuple[Component, ...] = ()
    elif topology is SurfaceTopology.HEXANE_ONLY:
        active_components = (Component.HEXANE,)
    elif topology is SurfaceTopology.WATER_ONLY:
        active_components = (Component.WATER,)
    else:
        active_components = (Component.WATER, Component.HEXANE)
    return SurfaceActiveSet(
        inventories=inventories,
        geometry=geometry,
        topology=topology,
        occupancy=occupancy,
        gas=gas,
        water_condition=water_condition,
        hexane_condition=hexane_condition,
        active_external_liquid_components=active_components,
        hexane_access_oracle=access_oracle,
        hexane_access_multiplier=access_multiplier,
        external_two_liquid_vlle_active=(
            topology is SurfaceTopology.TWO_LIQUID_VLLE
        ),
    )


def external_liquid_occupancy(
    inventories: ExternalSurfaceInventories,
    geometry: SurfaceGeometry,
    *,
    liquid_hexane_density_kg_m3: float,
    liquid_water_density_kg_m3: float,
) -> ExternalLiquidOccupancy:
    """Return exact ``S_h``, ``S_w``, ``S_L`` and the one shared area."""
    densities = (liquid_hexane_density_kg_m3, liquid_water_density_kg_m3)
    if not all(math.isfinite(value) and value > 0.0 for value in densities):
        raise ValueError("external liquid densities must be positive and finite")
    concentration = geometry.dry_meal_concentration_per_gas_void_kg_m3
    hexane_saturation = (
        concentration
        * inventories.attached_hexane_kg_per_kg_dry
        / liquid_hexane_density_kg_m3
    )
    water_saturation = (
        concentration
        * inventories.free_water_kg_per_kg_dry
        / liquid_water_density_kg_m3
    )
    total_saturation = math.fsum((hexane_saturation, water_saturation))
    if not all(
        math.isfinite(value) and value >= 0.0
        for value in (hexane_saturation, water_saturation, total_saturation)
    ):
        raise ValueError("external liquid occupancy is non-finite or negative")
    if total_saturation >= 1.0:
        raise ValueError(
            f"external liquids fill or exceed bed gas void: S_L={total_saturation!r}"
        )
    geometric_area = geometry.geometric_bed_area_m2_m3
    if total_saturation == 0.0:
        shared_area = geometric_area
    else:
        shared_area = geometric_area * (1.0 - total_saturation) ** (2.0 / 3.0)
    return ExternalLiquidOccupancy(
        dry_meal_concentration_per_gas_void_kg_m3=concentration,
        liquid_hexane_density_kg_m3=liquid_hexane_density_kg_m3,
        liquid_water_density_kg_m3=liquid_water_density_kg_m3,
        hexane_saturation=hexane_saturation,
        water_saturation=water_saturation,
        total_liquid_saturation=total_saturation,
        geometric_bed_area_m2_m3=geometric_area,
        shared_gas_side_area_m2_m3=shared_area,
    )


def external_water_hexane_access_bracket(
    inventories: ExternalSurfaceInventories,
    geometry: SurfaceGeometry,
    *,
    temperature_k: float,
    pressure_pa: float,
    y_water: float,
    y_hexane: float,
    gas_closure: GasClosure = GasClosure.PRESSURE_SECOND_VIRIAL,
    k_wh: float = bg.K_WH_CENTRAL,
) -> tuple[SurfaceActiveSet, SurfaceActiveSet]:
    """Return both access endpoints whenever an external water film exists."""
    if inventories.free_water_kg_per_kg_dry == 0.0:
        raise ValueError("n-hexane access bracket requires external water inventory")
    common = dict(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        y_water=y_water,
        y_hexane=y_hexane,
        gas_closure=gas_closure,
        k_wh=k_wh,
    )
    transparent = evaluate_surface_active_set(
        inventories,
        geometry,
        external_water_hexane_access=HexaneAccessOracle.TRANSPARENT_INTERFACE,
        **common,
    )
    impermeable = evaluate_surface_active_set(
        inventories,
        geometry,
        external_water_hexane_access=HexaneAccessOracle.IMPERMEABLE_WATER_FILM,
        **common,
    )
    return transparent, impermeable


def classify_one_sided_phase_event(
    component: Component,
    previous_inventory_kg_per_kg_dry: float,
    next_inventory_kg_per_kg_dry: float,
    endpoint_condition: ComponentPhaseCondition,
) -> OneSidedPhaseEvent:
    """Classify an exact inventory endpoint without a numerical deadband."""
    if not isinstance(component, Component):
        raise TypeError("phase-event component must be explicit")
    if endpoint_condition.component is not component:
        raise ValueError("phase-event condition belongs to the wrong component")
    values = (previous_inventory_kg_per_kg_dry, next_inventory_kg_per_kg_dry)
    if not all(math.isfinite(value) and value >= 0.0 for value in values):
        raise ValueError("phase-event inventories must be finite and non-negative")
    was_present = previous_inventory_kg_per_kg_dry > 0.0
    is_present = next_inventory_kg_per_kg_dry > 0.0
    if not was_present and not is_present:
        classification = InventoryEndpointClassification.ABSENT_CONTINUATION
    elif not was_present and is_present:
        classification = InventoryEndpointClassification.APPEARANCE
    elif was_present and not is_present:
        classification = InventoryEndpointClassification.DISAPPEARANCE
    else:
        classification = InventoryEndpointClassification.PRESENT_CONTINUATION
    at_zero_inventory_endpoint = classification is not (
        InventoryEndpointClassification.PRESENT_CONTINUATION
    )
    if endpoint_condition.liquid_present is at_zero_inventory_endpoint:
        raise ValueError(
            "phase-event condition must be evaluated on the selected inventory side"
        )
    return OneSidedPhaseEvent(
        component=component,
        previous_inventory_kg_per_kg_dry=previous_inventory_kg_per_kg_dry,
        next_inventory_kg_per_kg_dry=next_inventory_kg_per_kg_dry,
        inventory_classification=classification,
        endpoint_residual_sign=endpoint_condition.residual_sign,
        endpoint_phase_classification=endpoint_condition.classification,
    )


@dataclass(frozen=True)
class CommittedFilmLineageEvidence:
    """A5 (RULED 2026-08-22): the typed lineage evidence for the latch carve-out.

    Authority: ``docs/GT_PS2_A5_LATCH_CARVEOUT_RULING_PACKET_2026-08-22.md``
    section 3, ruled on the frozen PHY-014 text's own carve-out sentence
    ("EXISTING ADMISSIBLE FULL-CORE OR ATTACHED-LIQUID STATES MAY EVOLVE
    NORMALLY", ``release/physics_decisions.yaml#PHY-014``).  The evidence names
    the committed A1/A2-lineage film companion and the exact committed state
    it must be identity-bound to; :func:`require_committed_film_lineage`
    verifies both, plus the bit-for-bit identity of the companion's film with
    the inventory the latch is judging.  Nothing here is an inventory, a
    kinetic equation, or a tolerance.
    """

    film_companion: object
    committed_state: object
    physically_qualifying: bool = field(default=False, init=False)


def _committed_film_intensive(companion: object) -> object:
    """Read the companion's own attached-film intensive value, by lineage shape.

    The A1 sealed record carries it as ``film.after.attached_hexane_kg_per_kg_
    dry`` (the accumulator's committed endpoint); the BC-5 carrier that
    continues it carries the running value as ``attached_hexane_kg_per_kg_dry``
    directly.  Anything else has no committed film to read.
    """

    direct = getattr(companion, "attached_hexane_kg_per_kg_dry", None)
    if direct is not None:
        return direct
    film = getattr(companion, "film", None)
    after = getattr(film, "after", None)
    return getattr(after, "attached_hexane_kg_per_kg_dry", None)


def require_committed_film_lineage(
    evidence: CommittedFilmLineageEvidence,
    attached_hexane_kg_per_kg_dry: float,
) -> None:
    """Verify the A5-ruled lineage or refuse typed; admit nothing loosely.

    The carve-out admits recession under an attached film exactly when the
    film is a lawfully COMMITTED, IDENTITY-BOUND film of the A1/A2 lineage:

    * the companion is the exact A1 sealed record
      (``cut_birth_condensed_film.CommittedBirthAttachedFilm``) or the exact
      BC-5 carrier continuing it
      (``qsc_hexane_post_birth_stepper.PostBirthFilmCarrier`` whose
      ``birth_film_record`` is the exact A1 record);
    * ``is_bound_to``-verified against the declared committed state (BC-4a's
      identity discipline - the record can neither be re-attached to another
      state nor survive a state it does not describe);
    * the companion's film and the judged inventory are ONE number, bit for
      bit, and strictly positive.

    The lineage types are imported lazily: this foundation module stays
    import-free of the birth/stepper layers, and the check runs only when a
    caller claims the carve-out.
    """

    if type(evidence) is not CommittedFilmLineageEvidence:
        raise PrimaryDrainageEnvelopeError(
            "the latch carve-out consumes the exact typed "
            "CommittedFilmLineageEvidence (A5, RULED 2026-08-22); "
            f"got {type(evidence).__name__}"
        )
    companion = evidence.film_companion
    if companion is None or evidence.committed_state is None:
        raise PrimaryDrainageEnvelopeError(
            "the latch carve-out needs both the committed film companion and "
            "the committed state it binds; an absent member is no lineage"
        )
    from dtdc_simulator.core2.particle.cut_birth_condensed_film import (
        CommittedBirthAttachedFilm,
    )

    if type(companion) is CommittedBirthAttachedFilm:
        pass
    else:
        from dtdc_simulator.core2.qsc_hexane_post_birth_stepper import (
            PostBirthFilmCarrier,
        )

        if type(companion) is not PostBirthFilmCarrier:
            raise PrimaryDrainageEnvelopeError(
                "the latch carve-out admits only the A1/A2 lineage - the exact "
                "CommittedBirthAttachedFilm record or the exact "
                "PostBirthFilmCarrier continuing it; got "
                f"{type(companion).__name__}: any attached film WITHOUT the "
                "typed committed lineage keeps the frozen refusal"
            )
        if type(companion.birth_film_record) is not CommittedBirthAttachedFilm:
            raise PrimaryDrainageEnvelopeError(
                "the carrier's birth_film_record is not the exact A1 sealed "
                "CommittedBirthAttachedFilm; the lineage chain is broken and "
                "the frozen refusal stands"
            )
    bound = getattr(companion, "is_bound_to", None)
    if not callable(bound) or bound(evidence.committed_state) is not True:
        raise PrimaryDrainageEnvelopeError(
            "the committed film companion is not identity-bound to the "
            "declared committed state (BC-4a's is_bound_to discipline); an "
            "unbound companion proves no lineage"
        )
    companion_film = _committed_film_intensive(companion)
    if companion_film != attached_hexane_kg_per_kg_dry:
        raise PrimaryDrainageEnvelopeError(
            "the companion's committed film and the judged attached inventory "
            f"must be one number, bit for bit: companion {companion_film!r} vs "
            f"inventory {attached_hexane_kg_per_kg_dry!r}"
        )
    if not attached_hexane_kg_per_kg_dry > 0.0:
        raise PrimaryDrainageEnvelopeError(
            "the latch carve-out is defined only on a strictly positive "
            "committed film; at EXACTLY zero the latch has cleared and the "
            "ordinary machinery owns the state"
        )


def validate_primary_drainage_admission(
    active_set: SurfaceActiveSet,
    history: HexanePrimaryDrainageHistory,
    *,
    committed_film_lineage: CommittedFilmLineageEvidence | None = None,
) -> None:
    """Fail closed on forbidden external n-hexane reverse history.

    The caller, not the surface inventories, owns the irreversible history
    latch.  Once any central-core recession or dry-out has occurred, attached
    n-hexane may not reappear and a supersaturated absent-liquid state may not
    silently proceed into condensation/rewetting.  This check only admits or
    rejects a state; it supplies no nucleation, front, or kinetic equation and
    does not modify the independently evaluated VLLE residuals.

    THE A5 CARVE-OUT (RULED 2026-08-22, ``docs/GT_PS2_A5_LATCH_CARVEOUT_
    RULING_PACKET_2026-08-22.md`` section 3).  The frozen PHY-014 text carries
    its own admissibility sentence - "EXISTING ADMISSIBLE FULL-CORE OR
    ATTACHED-LIQUID STATES MAY EVOLVE NORMALLY" - and the ruling completes
    this latch to it: an attached film that is a lawfully COMMITTED,
    IDENTITY-BOUND film of the A1/A2 lineage (the typed companion, verified
    by :func:`require_committed_film_lineage`) is precisely such a state, and
    recession is ADMISSIBLE while it evolves under the A2 law in both signs.
    The refusal REMAINS byte-identical for (a) any attached film WITHOUT the
    typed lineage - ``committed_film_lineage=None``, every pre-A5 caller -
    and (b) the three named prohibitions (s=R reset, X_f recreation,
    internal-core nucleation), which this carve-out never touches: the
    condensation/rewetting refusal below is enforced exactly as before.
    """
    if not isinstance(active_set, SurfaceActiveSet):
        raise TypeError("surface active set is required")
    if not isinstance(history, HexanePrimaryDrainageHistory):
        raise TypeError("primary-drainage history must be supplied explicitly")
    if history is HexanePrimaryDrainageHistory.PRE_CORE_RECESSION:
        return
    if active_set.inventories.attached_hexane_kg_per_kg_dry > 0.0:
        if committed_film_lineage is None:
            raise PrimaryDrainageEnvelopeError(
                "attached n-hexane after core recession/dry-out is outside the "
                "qualified primary-drainage envelope"
            )
        if history is not HexanePrimaryDrainageHistory.CORE_RECESSION_OCCURRED:
            # The lineage is the BIRTH's: A1 commits the film on the
            # core-recession state, and the A2 axis's own scope says "a
            # dried-out core is outside A1's committed-film scope".  The
            # carve-out therefore never opens on a dry-out history.
            raise PrimaryDrainageEnvelopeError(
                "the A5 carve-out is scoped to the core-recession lineage; a "
                f"{history.value!r} history has no committed A1/A2 film to "
                "admit and the frozen refusal stands"
            )
        require_committed_film_lineage(
            committed_film_lineage,
            active_set.inventories.attached_hexane_kg_per_kg_dry,
        )
    if active_set.hexane_condition.classification is (
        PhaseResidualClassification.ABSENT_APPEARANCE_DEMANDED
    ):
        raise PrimaryDrainageEnvelopeError(
            "n-hexane condensation/rewetting after core recession/dry-out is "
            "outside the qualified primary-drainage envelope"
        )


def assemble_surface_transfer(
    active_set: SurfaceActiveSet,
    raw_component_flux_density: SurfaceComponentFluxVector,
    *,
    conductive_heat_flux_w_m2: float,
    coefficients: SurfaceTransferSelections,
) -> SurfaceTransferResult:
    """Apply one coupled vector and common-datum energy flux over one area.

    Positive flux is from the material/external-liquid side into the gas.
    The impermeable-water-film oracle sets only the n-hexane component to its
    frozen zero-access endpoint; no continuous blocking parameter exists.
    """
    if not isinstance(active_set, SurfaceActiveSet):
        raise TypeError("surface active set is required")
    if not isinstance(raw_component_flux_density, SurfaceComponentFluxVector):
        raise TypeError("one surface component-flux vector is required")
    if not isinstance(coefficients, SurfaceTransferSelections):
        raise TypeError("explicit surface coefficient selections are required")
    if not math.isfinite(conductive_heat_flux_w_m2):
        raise ValueError("surface conductive heat flux must be finite")
    effective_flux = SurfaceComponentFluxVector(
        water_mol_m2_s=raw_component_flux_density.water_mol_m2_s,
        hexane_mol_m2_s=(
            active_set.hexane_access_multiplier
            * raw_component_flux_density.hexane_mol_m2_s
        ),
    )
    water_enthalpy_flux = (
        effective_flux.water_mol_m2_s
        * active_set.gas.water_partial_enthalpy_j_mol
    )
    hexane_enthalpy_flux = (
        effective_flux.hexane_mol_m2_s
        * active_set.gas.hexane_partial_enthalpy_j_mol
    )
    total_energy_flux = math.fsum(
        (conductive_heat_flux_w_m2, water_enthalpy_flux, hexane_enthalpy_flux)
    )
    area = active_set.occupancy.shared_gas_side_area_m2_m3
    gas_water_source = area * effective_flux.water_mol_m2_s
    gas_hexane_source = area * effective_flux.hexane_mol_m2_s
    gas_energy_source = area * total_energy_flux
    return SurfaceTransferResult(
        active_set=active_set,
        coefficients=coefficients,
        raw_component_flux_density=raw_component_flux_density,
        effective_component_flux_density=effective_flux,
        conductive_heat_flux_w_m2=conductive_heat_flux_w_m2,
        shared_area_m2_m3=area,
        gas_water_source_mol_m3_s=gas_water_source,
        gas_hexane_source_mol_m3_s=gas_hexane_source,
        material_water_source_mol_m3_s=-gas_water_source,
        material_hexane_source_mol_m3_s=-gas_hexane_source,
        water_enthalpy_flux_w_m2=water_enthalpy_flux,
        hexane_enthalpy_flux_w_m2=hexane_enthalpy_flux,
        total_energy_flux_w_m2=total_energy_flux,
        gas_energy_source_w_m3=gas_energy_source,
        material_energy_source_w_m3=-gas_energy_source,
    )


def _surface_gas_thermo(
    temperature_k: float,
    pressure_pa: float,
    y_water: float,
    y_hexane: float,
    closure: GasClosure,
    k_wh: float,
) -> SurfaceGasThermo:
    _validate_temperature_pressure_composition(
        temperature_k,
        pressure_pa,
        y_water,
        y_hexane,
    )
    if closure is GasClosure.PRESSURE_SECOND_VIRIAL:
        mixture = bg.state(
            temperature_k,
            pressure_pa,
            y_water,
            y_hexane,
            k_wh=k_wh,
        )
        water_enthalpy = (
            bg.water_ideal_enthalpy_molar(temperature_k)
            + mixture.partial_residual_enthalpy_water_molar
        )
        hexane_enthalpy = (
            hx.h_ideal(temperature_k)
            + mixture.partial_residual_enthalpy_hexane_molar
        )
        return SurfaceGasThermo(
            closure=closure,
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            y_water=mixture.y_water,
            y_hexane=mixture.y_hexane,
            fugacity_water_pa=mixture.fugacity_water_pa,
            fugacity_hexane_pa=mixture.fugacity_hexane_pa,
            water_partial_enthalpy_j_mol=water_enthalpy,
            hexane_partial_enthalpy_j_mol=hexane_enthalpy,
            residual_enthalpy_molar_j_mol=mixture.residual_enthalpy_molar,
            below_cross_direct_evidence=mixture.below_cross_direct_evidence,
            outside_project_pressure=mixture.outside_project_pressure,
        )
    if closure is GasClosure.IDEAL_ABLATION:
        return SurfaceGasThermo(
            closure=closure,
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            y_water=y_water,
            y_hexane=y_hexane,
            fugacity_water_pa=y_water * pressure_pa,
            fugacity_hexane_pa=y_hexane * pressure_pa,
            water_partial_enthalpy_j_mol=bg.water_ideal_enthalpy_molar(temperature_k),
            hexane_partial_enthalpy_j_mol=hx.h_ideal(temperature_k),
            residual_enthalpy_molar_j_mol=0.0,
            below_cross_direct_evidence=False,
            outside_project_pressure=not (
                bg.PROJECT_PRESSURE_MIN_PA
                <= pressure_pa
                <= bg.PROJECT_PRESSURE_MAX_PA
            ),
        )
    raise AssertionError(f"unhandled gas closure: {closure!r}")


def classify_phase_condition(
    component: Component,
    inventory: float,
    gas_fugacity_pa: float,
    liquid_fugacity_pa: float,
) -> ComponentPhaseCondition:
    """Classify one exact fugacity residual with no tolerance or deadband."""
    if not isinstance(component, Component):
        raise TypeError("phase-condition component must be explicit")
    values = (inventory, gas_fugacity_pa, liquid_fugacity_pa)
    if not all(math.isfinite(value) for value in values):
        raise ValueError("phase-condition inputs must be finite")
    if inventory < 0.0:
        raise ValueError("phase-condition inventory cannot be negative")
    if gas_fugacity_pa < 0.0 or liquid_fugacity_pa <= 0.0:
        raise ValueError("phase-condition fugacities are outside their physical domain")
    residual = gas_fugacity_pa - liquid_fugacity_pa
    sign = _exact_sign(residual)
    present = inventory > 0.0
    if present:
        if sign is ResidualSign.NEGATIVE:
            classification = PhaseResidualClassification.ACTIVE_GAS_UNDERSATURATED
        elif sign is ResidualSign.EXACT_ZERO:
            classification = PhaseResidualClassification.ACTIVE_EQUALITY_SATISFIED
        else:
            classification = PhaseResidualClassification.ACTIVE_GAS_SUPERSATURATED
    elif sign is ResidualSign.NEGATIVE:
        classification = PhaseResidualClassification.ABSENT_STABLE
    elif sign is ResidualSign.EXACT_ZERO:
        classification = PhaseResidualClassification.ABSENT_ONE_SIDED_CONTACT
    else:
        classification = PhaseResidualClassification.ABSENT_APPEARANCE_DEMANDED
    return ComponentPhaseCondition(
        component=component,
        external_liquid_inventory_kg_per_kg_dry=inventory,
        liquid_present=present,
        gas_fugacity_pa=gas_fugacity_pa,
        pure_liquid_fugacity_pa=liquid_fugacity_pa,
        gas_minus_liquid_fugacity_residual_pa=residual,
        residual_sign=sign,
        classification=classification,
    )


def _exact_sign(value: float) -> ResidualSign:
    if not math.isfinite(value):
        raise ValueError("phase residual must be finite")
    if value < 0.0:
        return ResidualSign.NEGATIVE
    if value > 0.0:
        return ResidualSign.POSITIVE
    return ResidualSign.EXACT_ZERO


def _validate_temperature_pressure_composition(
    temperature_k: float,
    pressure_pa: float,
    y_water: float,
    y_hexane: float,
) -> None:
    if not math.isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValueError("surface temperature must be positive and finite")
    if not math.isfinite(pressure_pa) or pressure_pa <= 0.0:
        raise ValueError("surface pressure must be positive and finite")
    if not all(math.isfinite(value) and value >= 0.0 for value in (y_water, y_hexane)):
        raise ValueError("surface gas mole fractions must be finite and non-negative")
    if abs(y_water + y_hexane - 1.0) > 1.0e-12:
        raise ValueError("surface gas mole fractions must sum to one")


__all__ = [
    "CommittedFilmLineageEvidence",
    "Component",
    "ComponentPhaseCondition",
    "ExternalLiquidOccupancy",
    "ExternalSurfaceInventories",
    "GasClosure",
    "HexanePrimaryDrainageHistory",
    "HexaneAccessOracle",
    "InventoryEndpointClassification",
    "OneSidedPhaseEvent",
    "PhaseResidualClassification",
    "PrimaryDrainageEnvelopeError",
    "ResidualSign",
    "SurfaceActiveSet",
    "SurfaceComponentFluxVector",
    "SurfaceGeometry",
    "SurfaceTopology",
    "SurfaceTransferCoefficientInterval",
    "SurfaceTransferCoefficientSelection",
    "SurfaceTransferResult",
    "SurfaceTransferSelections",
    "TransferCoefficientKind",
    "assemble_surface_transfer",
    "classify_phase_condition",
    "classify_one_sided_phase_event",
    "evaluate_surface_active_set",
    "external_liquid_occupancy",
    "external_water_hexane_access_bracket",
    "require_committed_film_lineage",
    "validate_primary_drainage_admission",
]
